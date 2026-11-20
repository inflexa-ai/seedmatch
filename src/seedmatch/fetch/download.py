"""The download of one registry resource into a directory."""

import contextlib
import hashlib
import os
from collections.abc import Callable
from dataclasses import dataclass
from http.client import HTTPException
from pathlib import Path
from typing import Final
from urllib.parse import urljoin

from seedmatch._checks import as_limit, check_at_least, check_not_negative
from seedmatch.errors import FetchError, InputError
from seedmatch.fetch.http import HttpResponse, check_https, open_url
from seedmatch.resources import Entry, Resource

_OK: Final = 200
_REDIRECTS: Final = frozenset({301, 302, 303, 307, 308})


@dataclass(frozen=True)
class FetchOpts:
    """The options of a download.

    `open_url` sends one GET and follows no redirect. `chunk_size` bounds the memory of the read.
    """

    open_url: Callable[[str], HttpResponse]
    chunk_size: int = 1 << 20
    max_redirects: int = 5

    def __post_init__(self) -> None:
        """Store the limits as plain ints and make sure that they are valid."""
        # The dataclass is frozen, so the constructor is the one writer of the normalized values.
        object.__setattr__(self, "chunk_size", as_limit("chunk_size", self.chunk_size))
        object.__setattr__(self, "max_redirects", as_limit("max_redirects", self.max_redirects))
        check_at_least("chunk_size", self.chunk_size, 1)
        check_not_negative("max_redirects", self.max_redirects)


DEFAULT_FETCH_OPTS: Final = FetchOpts(open_url=open_url)


def fetch_resource(
    resource: Resource,
    dest: str | os.PathLike[str],
    *,
    opts: FetchOpts = DEFAULT_FETCH_OPTS,
) -> Path:
    """Download the resource to `dest / resource.name` and return that path.

    `dest` is a directory, made with its parents if it is absent. A file at the final path with the
    size and the digest of the resource is not downloaded again. A file with another digest raises
    `FetchError`, and stays. The body goes to `<name>.part`, and the file gets its final name only
    when its size and its digest agree with the resource. Each failure removes the `.part` file. A
    `.part` path that is not a regular file raises `FetchError` before the request.
    """
    directory = _make_directory(dest)
    final = directory / resource.name
    if _exists(final):
        _check_present(final, resource)
        return final
    part = directory / f"{resource.name}.part"
    _check_part(part)
    try:
        _transfer(resource, part, opts)
        part.replace(final)
    except BaseException:
        # The next run replaces a .part file that stays, and the caller needs the first error.
        with contextlib.suppress(OSError):
            part.unlink(missing_ok=True)
        raise
    return final


def fetch_entry(
    entry: Entry,
    dest: str | os.PathLike[str],
    *,
    opts: FetchOpts = DEFAULT_FETCH_OPTS,
) -> tuple[Path, ...]:
    """Fetch each resource of the entry in order, and return the paths."""
    return tuple(fetch_resource(resource, dest, opts=opts) for resource in entry.resources)


def _make_directory(dest: str | os.PathLike[str]) -> Path:
    directory = Path(dest)
    try:
        directory.mkdir(parents=True, exist_ok=True)
    except (FileExistsError, NotADirectoryError) as error:
        message = (
            f"dest {os.fspath(dest)!r} is not a directory; "
            "expected a directory or a path that is absent"
        )
        raise InputError(message) from error
    except OSError as error:
        message = (
            f"dest {os.fspath(dest)!r} cannot be made: {error!r}; "
            "expected a directory or a path that can be made"
        )
        raise FetchError(message) from error
    return directory


def _exists(final: Path) -> bool:
    try:
        return final.exists()
    except OSError as error:
        message = (
            f"{os.fspath(final)!r} cannot be accessed: {error!r}; "
            "expected a path that can be accessed"
        )
        raise FetchError(message) from error


def _check_part(part: Path) -> None:
    if part.is_symlink() or (part.exists() and not part.is_file()):
        message = (
            f"{os.fspath(part)!r} is not a regular file; expected a regular file or no file. "
            "Remove it to download again"
        )
        raise FetchError(message)


def _check_present(final: Path, resource: Resource) -> None:
    try:
        size = final.stat().st_size
        with final.open("rb") as stream:
            actual = hashlib.file_digest(stream, "sha256").hexdigest()
    except OSError as error:
        message = (
            f"{os.fspath(final)!r} cannot be read: {error!r}; expected a readable file or no file"
        )
        raise FetchError(message) from error
    if size != resource.size or actual != resource.sha256:
        message = (
            f"{os.fspath(final)!r} has {size} bytes and sha256 {actual}; "
            f"expected {resource.size} bytes and sha256 {resource.sha256}. "
            "Remove the file to download it again"
        )
        raise FetchError(message)


def _open_body(start: str, opts: FetchOpts) -> HttpResponse:
    """Return the response with the status 200, after the redirects of the start URL."""
    url, what, redirects = start, "URL", 0
    while True:
        check_https(url, what)
        response = opts.open_url(url)
        if response.status == _OK:
            return response
        try:
            if response.status not in _REDIRECTS:
                message = f"{url!r} answered with status {response.status}; expected {_OK}"
                raise FetchError(message)
            location = response.header("Location")
            if location is None:
                message = (
                    f"{url!r} answered with status {response.status} and no Location header; "
                    "expected a Location header"
                )
                raise FetchError(message)
            if redirects == opts.max_redirects:
                message = (
                    f"{start!r} has more than {opts.max_redirects} redirects; "
                    f"expected at most {opts.max_redirects}"
                )
                raise FetchError(message)
            target = _resolve(url, location)
        finally:
            response.close()
        url, what, redirects = target, "redirect target", redirects + 1


def _resolve(url: str, location: str) -> str:
    target = None
    # http.client encodes the request line as ASCII, and the registry URLs are ASCII, so a
    # non-ASCII Location is an error and not a value to percent-encode.
    if location.isascii():
        with contextlib.suppress(ValueError):  # an invalid IPv6 host raises ValueError
            target = urljoin(url, location)
    if target is None:
        message = f"{url!r} answered with Location {location!r}; expected a valid ASCII URL"
        raise FetchError(message)
    return target


def _transfer(resource: Resource, part: Path, opts: FetchOpts) -> None:
    """Write the body to `part`, and make sure that its size and its digest agree."""
    url, size = resource.url, resource.size
    digest = hashlib.sha256()
    received = 0
    try:
        response = _open_body(url, opts)
        try:
            with part.open("wb") as stream:
                # The read stops at 1 byte above the size: that byte shows a body that is too long.
                while chunk := response.read(min(opts.chunk_size, size + 1 - received)):
                    received += len(chunk)
                    if received > size:
                        message = f"{url!r} sent more than {size} bytes; expected {size}"
                        raise FetchError(message)
                    digest.update(chunk)
                    stream.write(chunk)
        finally:
            response.close()
    except (OSError, HTTPException) as error:
        message = f"{url!r} failed after {received} of {size} bytes: {error!r}"
        raise FetchError(message) from error
    if received != size:
        message = f"{url!r} sent {received} bytes; expected {size}"
        raise FetchError(message)
    if digest.hexdigest() != resource.sha256:
        message = f"{url!r} has sha256 {digest.hexdigest()}; expected {resource.sha256}"
        raise FetchError(message)
