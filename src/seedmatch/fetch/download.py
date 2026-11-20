"""The download of one registry resource into a directory."""

import contextlib
import errno
import fcntl
import hashlib
import os
import stat
from collections.abc import Callable, Iterator
from dataclasses import dataclass
from http.client import HTTPException
from pathlib import Path
from typing import Final
from urllib.parse import urljoin

from seedmatch._checks import as_limit, check_at_least, check_not_negative
from seedmatch.errors import FetchError, InputError
from seedmatch.fetch.http import HttpResponse, check_https, open_url
from seedmatch.resources import Resource

_OK: Final = 200
_REDIRECTS: Final = frozenset({301, 302, 303, 307, 308})
_LOCK_ATTEMPTS: Final = 8


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
    when its size and its digest agree with the resource. Each failure after the lock removes the
    `.part` file. A `.part` path that is not a regular file raises `FetchError` before the request.
    So does a `.part` file that another run holds: that run downloads the same file into the same
    `dest`, and this run leaves its file alone.
    """
    directory = _make_directory(dest)
    final = directory / resource.name
    if _exists(final):
        _check_present(final, resource)
        return final
    part = directory / f"{resource.name}.part"
    # The lock comes before the request. The rename and the delete happen while the run holds the
    # lock, so no run writes, renames, or deletes the .part file of another run.
    with _locked_part(part) as fd:
        try:
            # The run that held the lock before can have made the final file.
            if _exists(final):
                _remove_held(part, fd)
            else:
                _transfer(resource, fd, opts)
                _rename(part, final)
                return final
        except BaseException:
            _remove_held(part, fd)
            raise
    _check_present(final, resource)
    return final


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


def _not_regular(part: Path) -> str:
    return (
        f"{os.fspath(part)!r} is not a regular file; expected a regular file or no file. "
        "Remove it to download again"
    )


def _open_part(part: Path) -> int:
    flags = os.O_RDWR | os.O_CREAT | os.O_NOFOLLOW
    try:
        fd = os.open(part, flags, 0o644)
    except OSError as error:
        if error.errno in {errno.ELOOP, errno.EISDIR}:
            raise FetchError(_not_regular(part)) from error
        message = (
            f"{os.fspath(part)!r} cannot be opened: {error!r}; expected a file that can be written"
        )
        raise FetchError(message) from error
    if not stat.S_ISREG(os.fstat(fd).st_mode):
        os.close(fd)
        raise FetchError(_not_regular(part))
    return fd


def _lock(fd: int, part: Path) -> None:
    try:
        fcntl.flock(fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError as error:
        message = (
            f"{os.fspath(part)!r} is locked by another run that downloads the file into the same "
            "dest; expected one run for each file in a dest. Run again when the other run ends"
        )
        raise FetchError(message) from error
    except OSError as error:
        message = (
            f"{os.fspath(part)!r} cannot be locked: {error!r}; "
            "expected a file system that supports locks"
        )
        raise FetchError(message) from error


def _is_named(part: Path, fd: int) -> bool:
    try:
        named = part.lstat()
    except FileNotFoundError:
        return False
    held = os.fstat(fd)
    return (named.st_dev, named.st_ino) == (held.st_dev, held.st_ino)


def _open_locked(part: Path) -> int:
    for _ in range(_LOCK_ATTEMPTS):
        fd = _open_part(part)
        try:
            _lock(fd, part)
            named = _is_named(part, fd)
        except BaseException:
            os.close(fd)
            raise
        if named:
            return fd
        # The run that held the lock deleted or renamed this file before the lock came. That lock
        # protects nothing, because the path now names another file or none.
        os.close(fd)
    message = (
        f"{os.fspath(part)!r} names another file after each lock; "
        "expected one run for each file in a dest"
    )
    raise FetchError(message)


@contextlib.contextmanager
def _locked_part(part: Path) -> Iterator[int]:
    """Yield the descriptor of the `.part` file with an exclusive lock. The close ends the lock."""
    fd = _open_locked(part)
    try:
        yield fd
    finally:
        os.close(fd)


def _remove_held(part: Path, fd: int) -> None:
    # After the rename, the path can name the .part file of another run, so only the file of the
    # held lock is deleted. The next run replaces a .part file that stays, and the caller needs the
    # first error.
    with contextlib.suppress(OSError):
        if _is_named(part, fd):
            part.unlink()


def _rename(part: Path, final: Path) -> None:
    try:
        part.replace(final)
    except OSError as error:
        message = (
            f"{os.fspath(part)!r} cannot be renamed to {os.fspath(final)!r}: {error!r}; "
            "expected a directory that can be written"
        )
        raise FetchError(message) from error


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


def _check_length(response: HttpResponse, url: str, size: int) -> None:
    """Raise `FetchError` if a Content-Length header is not the size. No header is no error."""
    value = response.header("Content-Length")
    if value is None:
        return
    text = value.strip()
    if not (text.isascii() and text.isdigit()):
        message = f"{url!r} sent Content-Length {value!r}; expected a decimal number"
        raise FetchError(message)
    # The text is not converted to an int, because int() raises ValueError for a very long number.
    if (text.lstrip("0") or "0") != str(size):
        message = f"{url!r} announced {text} bytes in Content-Length; expected {size}"
        raise FetchError(message)


def _transfer(resource: Resource, fd: int, opts: FetchOpts) -> None:
    """Write the body to the file of `fd`, and make sure that its size and its digest agree."""
    url, size = resource.url, resource.size
    digest = hashlib.sha256()
    received = 0
    try:
        os.ftruncate(fd, 0)
        response = _open_body(url, opts)
        try:
            _check_length(response, url, size)
            # closefd=False: the caller owns the descriptor, and its close ends the lock.
            with os.fdopen(fd, "wb", closefd=False) as stream:
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
