import errno
import fcntl
import hashlib
import http.client
import os
import re
import threading
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pytest
from memory_client import MemoryClient, Route

from seedmatch.errors import FetchError, InputError, SeedMatchError
from seedmatch.fetch import DEFAULT_FETCH_OPTS, FetchOpts, fetch_resource
from seedmatch.fetch.http import open_url
from seedmatch.resources import Resource

URL = "https://example.org/data/file.bin"
BODY = bytes(range(256)) * 8
TIMEOUT = 10


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


def test_a_final_path_that_cannot_be_accessed_is_an_error(tmp_path, monkeypatch):
    final = tmp_path / "file.bin"
    error = PermissionError(13, "Permission denied")
    original_exists = Path.exists

    def failing_exists(self, *args, **kwargs):
        if self == final:
            raise error
        return original_exists(self, *args, **kwargs)

    client = MemoryClient({URL: Route(body=BODY)})
    message = (
        f"{os.fspath(final)!r} cannot be accessed: {error!r}; expected a path that can be accessed"
    )
    with monkeypatch.context() as patch:
        patch.setattr(Path, "exists", failing_exists)
        with pytest.raises(FetchError) as caught:
            fetch_resource(resource_of(), tmp_path, opts=opts_of(client))
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


@pytest.mark.parametrize("length", [len(BODY) - 1, len(BODY) + 1, 0])
def test_a_content_length_that_differs_from_the_size_is_an_error_before_the_read(tmp_path, length):
    client = MemoryClient({URL: Route(headers={"content-length": str(length)}, body=BODY)})
    message = f"{URL!r} announced {length} bytes in Content-Length; expected {len(BODY)}"
    with pytest.raises(FetchError) as caught:
        fetch_resource(resource_of(), tmp_path, opts=opts_of(client))
    assert str(caught.value) == message
    assert client.bytes_read == 0
    assert client.closed == 1
    assert names(tmp_path) == []


@pytest.mark.parametrize(
    "value", ["abc", "", "-2048", "2048.0", "2,048", "\u0662\u0660\u0664\u0668"]
)
def test_a_content_length_that_is_not_a_number_is_an_error_before_the_read(tmp_path, value):
    client = MemoryClient({URL: Route(headers={"Content-Length": value}, body=BODY)})
    message = f"{URL!r} sent Content-Length {value!r}; expected a decimal number"
    with pytest.raises(FetchError) as caught:
        fetch_resource(resource_of(), tmp_path, opts=opts_of(client))
    assert str(caught.value) == message
    assert client.bytes_read == 0
    assert client.closed == 1
    assert names(tmp_path) == []


@pytest.mark.parametrize("value", [str(len(BODY)), f" {len(BODY)} ", f"00{len(BODY)}"])
def test_a_content_length_that_agrees_with_the_size_is_downloaded(tmp_path, value):
    client = MemoryClient({URL: Route(headers={"Content-Length": value}, body=BODY)})
    path = fetch_resource(resource_of(), tmp_path, opts=opts_of(client))
    assert path.read_bytes() == BODY


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


@pytest.mark.parametrize("kind", ["directory", "link", "pipe"])
def test_a_part_path_that_is_not_a_file_is_an_error_before_a_request(tmp_path, kind):
    dest = tmp_path / "dest"
    dest.mkdir()
    outside = tmp_path / "outside"
    outside.write_bytes(b"keep")
    part = dest / "file.bin.part"
    if kind == "directory":
        part.mkdir()
    elif kind == "pipe":
        os.mkfifo(part)
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


def paused(client, reached, release):
    """Wrap the client: the first read tells `reached`, then it waits until `release` is set."""

    def open_url(url):
        response = client(url)
        read = response.read

        def read_when_released(size):
            reached.set()
            if not release.wait(TIMEOUT):
                message = "the test did not release the read"
                raise TimeoutError(message)
            return read(size)

        response.read = read_when_released
        return response

    return open_url


def locked_message(part):
    return (
        f"{os.fspath(part)!r} is locked by another run that downloads the file into the same dest; "
        "expected one run for each file in a dest. Run again when the other run ends"
    )


def test_two_runs_into_one_dest_give_one_file_and_one_error(tmp_path):
    reached, release = threading.Event(), threading.Event()
    first_client = MemoryClient({URL: Route(body=BODY)})
    second_client = MemoryClient({URL: Route(body=BODY)})
    first_opts = opts_of(paused(first_client, reached, release))
    with ThreadPoolExecutor(max_workers=1) as pool:
        first = pool.submit(fetch_resource, resource_of(), tmp_path, opts=first_opts)
        try:
            assert reached.wait(TIMEOUT)
            with pytest.raises(FetchError) as caught:
                fetch_resource(resource_of(), tmp_path, opts=opts_of(second_client))
            assert names(tmp_path) == ["file.bin.part"]
        finally:
            release.set()
        path = first.result(timeout=TIMEOUT)
    assert str(caught.value) == locked_message(tmp_path / "file.bin.part")
    assert isinstance(caught.value.__cause__, BlockingIOError)
    assert second_client.requests == []
    assert path == tmp_path / "file.bin"
    assert path.read_bytes() == BODY
    assert names(tmp_path) == ["file.bin"]


def test_a_part_file_locked_by_another_run_is_an_error_and_stays(tmp_path):
    part = tmp_path / "file.bin.part"
    client = MemoryClient({URL: Route(body=BODY)})
    holder = os.open(part, os.O_RDWR | os.O_CREAT, 0o644)
    try:
        fcntl.flock(holder, fcntl.LOCK_EX | fcntl.LOCK_NB)
        os.write(holder, b"bytes of the other run")
        with pytest.raises(FetchError) as caught:
            fetch_resource(resource_of(), tmp_path, opts=opts_of(client))
        assert part.read_bytes() == b"bytes of the other run"
        assert names(tmp_path) == ["file.bin.part"]
    finally:
        os.close(holder)
    assert str(caught.value) == locked_message(part)
    assert client.requests == []

    path = fetch_resource(resource_of(), tmp_path, opts=opts_of(client))
    assert path.read_bytes() == BODY
    assert names(tmp_path) == ["file.bin"]


def test_a_final_file_made_before_the_lock_is_returned_with_no_request(tmp_path, monkeypatch):
    final = tmp_path / "file.bin"
    flock = fcntl.flock

    def flock_after_the_other_run(fd, operation):
        final.write_bytes(BODY)
        flock(fd, operation)

    client = MemoryClient({URL: Route(body=b"other")})
    with monkeypatch.context() as patch:
        patch.setattr(fcntl, "flock", flock_after_the_other_run)
        path = fetch_resource(resource_of(), tmp_path, opts=opts_of(client))
    assert path == final
    assert path.read_bytes() == BODY
    assert client.requests == []
    assert names(tmp_path) == ["file.bin"]


def test_a_final_file_with_another_digest_made_before_the_lock_is_an_error(tmp_path, monkeypatch):
    final = tmp_path / "file.bin"
    other = bytes(reversed(BODY))
    flock = fcntl.flock

    def flock_after_the_other_run(fd, operation):
        final.write_bytes(other)
        flock(fd, operation)

    client = MemoryClient({URL: Route(body=BODY)})
    with monkeypatch.context() as patch:
        patch.setattr(fcntl, "flock", flock_after_the_other_run)
        with pytest.raises(FetchError) as caught:
            fetch_resource(resource_of(), tmp_path, opts=opts_of(client))
    assert str(caught.value) == (
        f"{os.fspath(final)!r} has {len(BODY)} bytes and sha256 {digest_of(other)}; "
        f"expected {len(BODY)} bytes and sha256 {digest_of(BODY)}. "
        "Remove the file to download it again"
    )
    assert client.requests == []
    assert final.read_bytes() == other
    assert names(tmp_path) == ["file.bin"]


def test_a_part_path_that_changes_after_each_lock_is_an_error(tmp_path, monkeypatch):
    part = tmp_path / "file.bin.part"
    lstat = Path.lstat

    def lstat_of_another_file(self):
        result = lstat(self)
        if self != part:
            return result
        return os.stat_result((*result[:1], result.st_ino + 1, *result[2:]))

    client = MemoryClient({URL: Route(body=BODY)})
    with monkeypatch.context() as patch:
        patch.setattr(Path, "lstat", lstat_of_another_file)
        with pytest.raises(FetchError) as caught:
            fetch_resource(resource_of(), tmp_path, opts=opts_of(client))
    assert str(caught.value) == (
        f"{os.fspath(part)!r} names another file after each lock; "
        "expected one run for each file in a dest"
    )
    assert client.requests == []


def test_an_interrupt_after_the_rename_keeps_the_part_file_of_another_run(tmp_path, monkeypatch):
    part = tmp_path / "file.bin.part"
    replace = Path.replace

    def replace_then_another_run_starts(self, target):
        replace(self, target)
        part.write_bytes(b"part file of another run")
        raise KeyboardInterrupt

    client = MemoryClient({URL: Route(body=BODY)})
    with monkeypatch.context() as patch:
        patch.setattr(Path, "replace", replace_then_another_run_starts)
        with pytest.raises(KeyboardInterrupt):
            fetch_resource(resource_of(), tmp_path, opts=opts_of(client))
    assert (tmp_path / "file.bin").read_bytes() == BODY
    assert part.read_bytes() == b"part file of another run"


@pytest.mark.parametrize(
    "replacement", [b"bytes of a third run", None], ids=["replaced", "deleted"]
)
def test_a_part_file_gone_before_the_lock_is_opened_again(tmp_path, monkeypatch, replacement):
    part = tmp_path / "file.bin.part"
    flock = fcntl.flock
    locked = []

    def flock_after_another_run_replaced_the_file(fd, operation):
        if not locked:
            part.unlink()
            if replacement is not None:
                part.write_bytes(replacement)
        locked.append(fd)
        flock(fd, operation)

    client = MemoryClient({URL: Route(body=BODY)})
    with monkeypatch.context() as patch:
        patch.setattr(fcntl, "flock", flock_after_another_run_replaced_the_file)
        path = fetch_resource(resource_of(), tmp_path, opts=opts_of(client))
    assert len(locked) == 2
    assert path.read_bytes() == BODY
    assert client.requests == [URL]
    assert names(tmp_path) == ["file.bin"]


def test_a_part_file_that_cannot_be_locked_is_an_error_before_a_request(tmp_path, monkeypatch):
    part = tmp_path / "file.bin.part"
    error = OSError(errno.ENOLCK, "No locks available")

    def failing_flock(_fd, _operation):
        raise error

    client = MemoryClient({URL: Route(body=BODY)})
    message = (
        f"{os.fspath(part)!r} cannot be locked: {error!r}; "
        "expected a file system that supports locks"
    )
    with monkeypatch.context() as patch:
        patch.setattr(fcntl, "flock", failing_flock)
        with pytest.raises(FetchError) as caught:
            fetch_resource(resource_of(), tmp_path, opts=opts_of(client))
    assert str(caught.value) == message
    assert caught.value.__cause__ is error
    assert client.requests == []


def test_a_part_file_that_cannot_be_opened_is_an_error_before_a_request(tmp_path, monkeypatch):
    part = tmp_path / "file.bin.part"
    error = PermissionError(13, "Permission denied")
    original_open = os.open

    def failing_open(path, *args, **kwargs):
        if Path(path) == part:
            raise error
        return original_open(path, *args, **kwargs)

    client = MemoryClient({URL: Route(body=BODY)})
    message = (
        f"{os.fspath(part)!r} cannot be opened: {error!r}; expected a file that can be written"
    )
    with monkeypatch.context() as patch:
        patch.setattr(os, "open", failing_open)
        with pytest.raises(FetchError) as caught:
            fetch_resource(resource_of(), tmp_path, opts=opts_of(client))
    assert str(caught.value) == message
    assert caught.value.__cause__ is error
    assert client.requests == []
    assert names(tmp_path) == []


def test_an_error_at_the_rename_is_a_fetch_error_and_leaves_no_file(tmp_path, monkeypatch):
    part = tmp_path / "file.bin.part"
    final = tmp_path / "file.bin"
    error = PermissionError(13, "Permission denied")

    def failing_replace(_self, _target):
        raise error

    client = MemoryClient({URL: Route(body=BODY)})
    message = (
        f"{os.fspath(part)!r} cannot be renamed to {os.fspath(final)!r}: {error!r}; "
        "expected a directory that can be written"
    )
    with monkeypatch.context() as patch:
        patch.setattr(Path, "replace", failing_replace)
        with pytest.raises(FetchError) as caught:
            fetch_resource(resource_of(), tmp_path, opts=opts_of(client))
    assert str(caught.value) == message
    assert caught.value.__cause__ is error
    assert client.requests == [URL]
    assert names(tmp_path) == []


class FullDisk:
    """A stream that writes its first chunk and raises `error` at the second."""

    def __init__(self, stream, error):
        self.stream = stream
        self.error = error
        self.writes = 0

    def __enter__(self):
        return self

    def __exit__(self, *_exception):
        self.stream.close()

    def write(self, chunk):
        self.writes += 1
        if self.writes == 2:
            raise self.error
        return self.stream.write(chunk)


def test_a_write_error_during_the_transfer_is_a_fetch_error_and_leaves_no_file(
    tmp_path, monkeypatch
):
    error = OSError(errno.ENOSPC, "No space left on device")
    fdopen = os.fdopen
    client = MemoryClient({URL: Route(body=BODY)})
    message = f"{URL!r} failed after 200 of {len(BODY)} bytes: {error!r}"
    with monkeypatch.context() as patch:
        patch.setattr(
            os, "fdopen", lambda *args, **kwargs: FullDisk(fdopen(*args, **kwargs), error)
        )
        with pytest.raises(FetchError) as caught:
            fetch_resource(resource_of(), tmp_path, opts=opts_of(client))
    assert str(caught.value) == message
    assert caught.value.__cause__ is error
    assert names(tmp_path) == []
    assert client.closed == 1
