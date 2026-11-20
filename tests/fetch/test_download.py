import hashlib
import http.client
import os
import re
from pathlib import Path

import pytest
from memory_client import MemoryClient, Route

from seedmatch.errors import FetchError, InputError, SeedMatchError
from seedmatch.fetch import DEFAULT_FETCH_OPTS, FetchOpts, fetch_entry, fetch_resource
from seedmatch.fetch.http import open_url
from seedmatch.resources import Entry, Group, Resource

URL = "https://example.org/data/file.bin"
BODY = bytes(range(256)) * 8


def digest_of(body):
    return hashlib.sha256(body).hexdigest()


def resource_of(body=BODY, url=URL):
    return Resource(url, len(body), digest_of(body))


def opts_of(client, **kwargs):
    return FetchOpts(open_url=client, chunk_size=100, **kwargs)


def names(directory):
    return sorted(path.name for path in directory.iterdir())


def test_fetch_error_is_a_package_error():
    assert issubclass(FetchError, SeedMatchError)
    assert not issubclass(FetchError, ValueError)


def test_the_default_opts_use_the_default_client():
    assert DEFAULT_FETCH_OPTS.open_url is open_url
    assert DEFAULT_FETCH_OPTS.chunk_size == 1 << 20
    assert DEFAULT_FETCH_OPTS.max_redirects == 5


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        ({"chunk_size": 0}, "chunk_size 0 is less than 1; expected 1 or more"),
        ({"chunk_size": True}, "chunk_size True has type bool; expected an integer"),
        ({"chunk_size": "1"}, "chunk_size '1' has type str; expected an integer"),
        ({"max_redirects": -1}, "max_redirects -1 is negative; expected 0 or more"),
        ({"max_redirects": 1.5}, "max_redirects 1.5 has type float; expected an integer"),
    ],
)
def test_opts_reject_a_bad_limit(kwargs, message):
    with pytest.raises(InputError, match=re.escape(message)) as caught:
        FetchOpts(open_url=MemoryClient({}), **kwargs)
    assert str(caught.value) == message


def test_fetch_writes_the_file_and_leaves_no_part_file(tmp_path):
    client = MemoryClient({URL: Route(body=BODY)})
    path = fetch_resource(resource_of(), tmp_path, opts=opts_of(client))
    assert path == tmp_path / "file.bin"
    assert path.read_bytes() == BODY
    assert names(tmp_path) == ["file.bin"]
    assert client.requests == [URL]
    assert client.closed == 1
    assert client.max_read_size == 100


def test_fetch_accepts_the_dest_as_text(tmp_path):
    client = MemoryClient({URL: Route(body=BODY)})
    path = fetch_resource(resource_of(), str(tmp_path), opts=opts_of(client))
    assert path == tmp_path / "file.bin"


def test_fetch_makes_the_dest_directory(tmp_path):
    dest = tmp_path / "a" / "b"
    client = MemoryClient({URL: Route(body=BODY)})
    path = fetch_resource(resource_of(), dest, opts=opts_of(client))
    assert path == dest / "file.bin"
    assert path.read_bytes() == BODY


def bad_dest(kind, root):
    if kind == "file":
        dest = root / "dest"
        dest.write_bytes(b"x")
    elif kind == "under-a-file":
        (root / "file").write_bytes(b"x")
        dest = root / "file" / "dest"
    else:
        dest = root / "dest"
        dest.symlink_to(root / "absent")
    return dest


@pytest.mark.parametrize("kind", ["file", "under-a-file", "dangling-link"])
def test_fetch_rejects_a_dest_that_is_not_a_directory(tmp_path, kind):
    dest = bad_dest(kind, tmp_path)
    before = names(tmp_path)
    client = MemoryClient({URL: Route(body=BODY)})
    message = (
        f"dest {os.fspath(dest)!r} is not a directory; "
        "expected a directory or a path that is absent"
    )
    with pytest.raises(InputError) as caught:
        fetch_resource(resource_of(), dest, opts=opts_of(client))
    assert str(caught.value) == message
    assert client.requests == []
    assert names(tmp_path) == before


def test_fetch_reports_a_dest_that_cannot_be_made(tmp_path, monkeypatch):
    dest = tmp_path / "dest"
    error = PermissionError(13, "Permission denied")

    def failing_mkdir(*_args, **_kwargs):
        raise error

    client = MemoryClient({URL: Route(body=BODY)})
    message = (
        f"dest {os.fspath(dest)!r} cannot be made: {error!r}; "
        "expected a directory or a path that can be made"
    )
    with monkeypatch.context() as patch:
        patch.setattr(Path, "mkdir", failing_mkdir)
        with pytest.raises(FetchError) as caught:
            fetch_resource(resource_of(), dest, opts=opts_of(client))
    assert str(caught.value) == message
    assert caught.value.__cause__ is error
    assert client.requests == []
    assert names(tmp_path) == []


def test_a_present_file_with_the_digest_causes_no_request(tmp_path):
    (tmp_path / "file.bin").write_bytes(BODY)
    client = MemoryClient({URL: Route(body=b"other")})
    path = fetch_resource(resource_of(), tmp_path, opts=opts_of(client))
    assert path == tmp_path / "file.bin"
    assert client.requests == []
    assert path.read_bytes() == BODY


@pytest.mark.parametrize(
    "content", [b"other", BODY[:-1] + b"\x00"], ids=["other-size", "same-size"]
)
def test_a_present_file_with_another_digest_is_an_error_and_stays(tmp_path, content):
    final = tmp_path / "file.bin"
    final.write_bytes(content)
    client = MemoryClient({URL: Route(body=BODY)})
    message = (
        f"{os.fspath(final)!r} has {len(content)} bytes and sha256 {digest_of(content)}; "
        f"expected {len(BODY)} bytes and sha256 {digest_of(BODY)}. "
        "Remove the file to download it again"
    )
    with pytest.raises(FetchError) as caught:
        fetch_resource(resource_of(), tmp_path, opts=opts_of(client))
    assert str(caught.value) == message
    assert final.read_bytes() == content
    assert client.requests == []
    assert names(tmp_path) == ["file.bin"]


def test_a_present_file_that_cannot_be_read_is_an_error(tmp_path):
    final = tmp_path / "file.bin"
    final.mkdir()
    client = MemoryClient({URL: Route(body=BODY)})
    with pytest.raises(FetchError) as caught:
        fetch_resource(resource_of(), tmp_path, opts=opts_of(client))
    assert isinstance(caught.value.__cause__, OSError)
    message = (
        f"{os.fspath(final)!r} cannot be read: {caught.value.__cause__!r}; "
        "expected a readable file or no file"
    )
    assert str(caught.value) == message
    assert client.requests == []


def test_a_digest_mismatch_leaves_no_file(tmp_path):
    wrong = BODY[::-1]
    client = MemoryClient({URL: Route(body=wrong)})
    message = f"{URL!r} has sha256 {digest_of(wrong)}; expected {digest_of(BODY)}"
    with pytest.raises(FetchError) as caught:
        fetch_resource(resource_of(), tmp_path, opts=opts_of(client))
    assert str(caught.value) == message
    assert names(tmp_path) == []
    assert client.closed == 1


def test_a_short_body_leaves_no_file(tmp_path):
    client = MemoryClient({URL: Route(body=BODY[:-1])})
    message = f"{URL!r} sent {len(BODY) - 1} bytes; expected {len(BODY)}"
    with pytest.raises(FetchError) as caught:
        fetch_resource(resource_of(), tmp_path, opts=opts_of(client))
    assert str(caught.value) == message
    assert names(tmp_path) == []
    assert client.closed == 1


@pytest.mark.parametrize("chunk_size", [100, 1 << 20])
def test_a_long_body_stops_at_the_size_plus_one_byte(tmp_path, chunk_size):
    client = MemoryClient({URL: Route(body=BODY * 10)})
    message = f"{URL!r} sent more than {len(BODY)} bytes; expected {len(BODY)}"
    opts = FetchOpts(open_url=client, chunk_size=chunk_size)
    with pytest.raises(FetchError) as caught:
        fetch_resource(resource_of(), tmp_path, opts=opts)
    assert str(caught.value) == message
    assert client.bytes_read == len(BODY) + 1
    assert client.max_read_size == min(chunk_size, len(BODY) + 1)
    assert names(tmp_path) == []
    assert client.closed == 1


@pytest.mark.parametrize(
    "error",
    [ConnectionResetError("connection reset by peer"), http.client.IncompleteRead(b"abc", 7)],
    ids=["os-error", "http-error"],
)
def test_an_interrupted_transfer_leaves_no_file(tmp_path, error):
    body = bytes(range(256)) * 40
    reset_at = len(body) * 6 // 10
    client = MemoryClient({URL: Route(body=body, reset_at=reset_at, error=error)})
    opts = FetchOpts(open_url=client, chunk_size=1000)
    with pytest.raises(FetchError) as caught:
        fetch_resource(resource_of(body), tmp_path, opts=opts)
    assert caught.value.__cause__ is error
    message = f"{URL!r} failed after {reset_at} of {len(body)} bytes: {error!r}"
    assert str(caught.value) == message
    assert names(tmp_path) == []
    assert client.closed == 1
    assert client.max_read_size <= 1000

    client.routes[URL] = Route(body=body)
    path = fetch_resource(resource_of(body), tmp_path, opts=opts)
    assert path.read_bytes() == body
    assert names(tmp_path) == ["file.bin"]


@pytest.mark.parametrize(
    "error",
    [TimeoutError("timed out"), http.client.IncompleteRead(b"ab", 5)],
    ids=["os-error", "http-error"],
)
def test_an_error_when_the_request_opens_is_a_fetch_error(tmp_path, error):
    def failing(_url):
        raise error

    with pytest.raises(FetchError) as caught:
        fetch_resource(resource_of(), tmp_path, opts=FetchOpts(open_url=failing))
    assert caught.value.__cause__ is error
    assert str(caught.value) == f"{URL!r} failed after 0 of {len(BODY)} bytes: {error!r}"
    assert names(tmp_path) == []


def test_a_keyboard_interrupt_leaves_no_file(tmp_path):
    client = MemoryClient({URL: Route(body=BODY)})

    def interrupting(url):
        response = client(url)
        read = response.read

        def read_once(size):
            if client.bytes_read >= 300:
                raise KeyboardInterrupt
            return read(size)

        response.read = read_once
        return response

    opts = FetchOpts(open_url=interrupting, chunk_size=100)
    with pytest.raises(KeyboardInterrupt):
        fetch_resource(resource_of(), tmp_path, opts=opts)
    assert client.bytes_read == 300
    assert client.closed == 1
    assert names(tmp_path) == []


def test_a_stale_part_file_is_replaced(tmp_path):
    (tmp_path / "file.bin.part").write_bytes(b"stale" * 1000)
    client = MemoryClient({URL: Route(body=BODY)})
    path = fetch_resource(resource_of(), tmp_path, opts=opts_of(client))
    assert path.read_bytes() == BODY
    assert names(tmp_path) == ["file.bin"]


@pytest.mark.parametrize("kind", ["directory", "link"])
def test_a_part_path_that_is_not_a_file_is_an_error_before_a_request(tmp_path, kind):
    dest = tmp_path / "dest"
    dest.mkdir()
    outside = tmp_path / "outside"
    outside.write_bytes(b"keep")
    part = dest / "file.bin.part"
    if kind == "directory":
        part.mkdir()
    else:
        part.symlink_to(outside)
    client = MemoryClient({URL: Route(body=BODY)})
    message = (
        f"{os.fspath(part)!r} is not a regular file; expected a regular file or no file. "
        "Remove it to download again"
    )
    with pytest.raises(FetchError) as caught:
        fetch_resource(resource_of(), dest, opts=opts_of(client))
    assert str(caught.value) == message
    assert client.requests == []
    assert names(dest) == ["file.bin.part"]
    assert part.is_dir() == (kind == "directory")
    assert outside.read_bytes() == b"keep"


def test_an_error_in_the_cleanup_does_not_replace_the_first_error(tmp_path, monkeypatch):
    wrong = BODY[::-1]
    client = MemoryClient({URL: Route(body=wrong)})

    def failing_unlink(*_args, **_kwargs):
        raise PermissionError(13, "Permission denied")

    message = f"{URL!r} has sha256 {digest_of(wrong)}; expected {digest_of(BODY)}"
    with monkeypatch.context() as patch:
        patch.setattr(Path, "unlink", failing_unlink)
        with pytest.raises(FetchError) as caught:
            fetch_resource(resource_of(), tmp_path, opts=opts_of(client))
    assert str(caught.value) == message
    assert client.closed == 1


def test_an_http_url_is_an_error_before_a_request(tmp_path):
    url = "http://example.org/data/file.bin"
    client = MemoryClient({url: Route(body=BODY)})
    message = f"URL {url!r} is not https; expected a URL with the scheme https"
    with pytest.raises(FetchError) as caught:
        fetch_resource(resource_of(url=url), tmp_path, opts=opts_of(client))
    assert str(caught.value) == message
    assert client.requests == []
    assert names(tmp_path) == []


@pytest.mark.parametrize("status", [301, 302, 303, 307, 308])
def test_a_redirect_to_https_is_followed(tmp_path, status):
    target = "https://cdn.example.org/store/file.bin"
    client = MemoryClient(
        {
            URL: Route(status=status, headers={"Location": target}, body=b"moved"),
            target: Route(body=BODY),
        }
    )
    path = fetch_resource(resource_of(), tmp_path, opts=opts_of(client))
    assert path.read_bytes() == BODY
    assert client.requests == [URL, target]
    assert client.closed == 2


@pytest.mark.parametrize(
    ("location", "target"),
    [
        ("/other/file.bin", "https://example.org/other/file.bin"),
        ("../up/file.bin", "https://example.org/up/file.bin"),
        ("same.bin", "https://example.org/data/same.bin"),
    ],
)
def test_a_relative_redirect_is_resolved(tmp_path, location, target):
    client = MemoryClient(
        {URL: Route(status=302, headers={"location": location}), target: Route(body=BODY)}
    )
    path = fetch_resource(resource_of(), tmp_path, opts=opts_of(client))
    assert path.read_bytes() == BODY
    assert client.requests == [URL, target]


def test_a_redirect_to_http_is_an_error_before_the_request(tmp_path):
    target = "http://example.org/data/file.bin"
    client = MemoryClient(
        {URL: Route(status=302, headers={"Location": target}), target: Route(body=BODY)}
    )
    message = f"redirect target {target!r} is not https; expected a URL with the scheme https"
    with pytest.raises(FetchError) as caught:
        fetch_resource(resource_of(), tmp_path, opts=opts_of(client))
    assert str(caught.value) == message
    assert client.requests == [URL]
    assert client.closed == 1
    assert names(tmp_path) == []


def test_a_redirect_without_a_location_is_an_error(tmp_path):
    client = MemoryClient({URL: Route(status=302)})
    message = f"{URL!r} answered with status 302 and no Location header; expected a Location header"
    with pytest.raises(FetchError) as caught:
        fetch_resource(resource_of(), tmp_path, opts=opts_of(client))
    assert str(caught.value) == message
    assert client.closed == 1


@pytest.mark.parametrize("location", ["https://[::1/x", "/café/file.bin", "https://exämple.org/x"])
def test_a_location_that_cannot_be_used_is_an_error(tmp_path, location):
    client = MemoryClient({URL: Route(status=302, headers={"Location": location})})
    message = f"{URL!r} answered with Location {location!r}; expected a valid ASCII URL"
    with pytest.raises(FetchError) as caught:
        fetch_resource(resource_of(), tmp_path, opts=opts_of(client))
    assert str(caught.value) == message
    assert client.requests == [URL]
    assert client.closed == 1
    assert names(tmp_path) == []


@pytest.mark.parametrize("max_redirects", [0, 2])
def test_too_many_redirects_is_an_error(tmp_path, max_redirects):
    client = MemoryClient({URL: Route(status=302, headers={"Location": URL})})
    message = f"{URL!r} has more than {max_redirects} redirects; expected at most {max_redirects}"
    with pytest.raises(FetchError) as caught:
        fetch_resource(resource_of(), tmp_path, opts=opts_of(client, max_redirects=max_redirects))
    assert str(caught.value) == message
    assert client.requests == [URL] * (max_redirects + 1)
    assert client.closed == max_redirects + 1


def test_as_many_redirects_as_the_limit_are_followed(tmp_path):
    first = "https://example.org/one/file.bin"
    second = "https://example.org/two/file.bin"
    client = MemoryClient(
        {
            URL: Route(status=302, headers={"Location": first}),
            first: Route(status=302, headers={"Location": second}),
            second: Route(body=BODY),
        }
    )
    path = fetch_resource(resource_of(), tmp_path, opts=opts_of(client, max_redirects=2))
    assert path.read_bytes() == BODY
    assert client.requests == [URL, first, second]


@pytest.mark.parametrize("status", [204, 304, 404, 503])
def test_a_status_other_than_200_is_an_error_with_the_status(tmp_path, status):
    client = MemoryClient({URL: Route(status=status, body=BODY)})
    message = f"{URL!r} answered with status {status}; expected 200"
    with pytest.raises(FetchError) as caught:
        fetch_resource(resource_of(), tmp_path, opts=opts_of(client))
    assert str(caught.value) == message
    assert names(tmp_path) == []
    assert client.closed == 1


def test_fetch_entry_returns_the_paths_in_order(tmp_path):
    second = b"second file"
    url_second = "https://example.org/data/second.bin"
    client = MemoryClient({URL: Route(body=BODY), url_second: Route(body=second)})
    entry = Entry(
        Group.EXAMPLE,
        "pair",
        "two files",
        (resource_of(second, url_second), resource_of()),
    )
    paths = fetch_entry(entry, tmp_path, opts=opts_of(client))
    assert paths == (tmp_path / "second.bin", tmp_path / "file.bin")
    assert [path.read_bytes() for path in paths] == [second, BODY]
    assert client.requests == [url_second, URL]


def test_fetch_entry_stops_at_the_first_error(tmp_path):
    url_second = "https://example.org/data/second.bin"
    client = MemoryClient({URL: Route(status=404), url_second: Route(body=b"x")})
    entry = Entry(
        Group.EXAMPLE, "pair", "two files", (resource_of(), resource_of(b"x", url_second))
    )
    with pytest.raises(FetchError) as caught:
        fetch_entry(entry, tmp_path, opts=opts_of(client))
    assert str(caught.value) == f"{URL!r} answered with status 404; expected 200"
    assert client.requests == [URL]
