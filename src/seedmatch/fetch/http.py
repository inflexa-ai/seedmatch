"""The default HTTP client: one GET with `urllib.request`, and no redirect."""

import urllib.error
import urllib.request
from email.message import Message
from importlib.metadata import version
from typing import IO, Final, Protocol, override
from urllib.parse import urlsplit

from seedmatch.errors import FetchError

_TIMEOUT: Final = 60


class HttpResponse(Protocol):
    """The part of an HTTP response that a download reads."""

    status: int

    def header(self, name: str) -> str | None:
        """Return the value of the header, or `None` if the response has no such header."""
        ...

    def read(self, size: int) -> bytes:
        """Return the next bytes of the body, at most `size`. An empty result is the end."""
        ...

    def close(self) -> None:
        """Release the connection."""
        ...


class _Raw(Protocol):
    @property
    def headers(self) -> Message: ...

    def read(self, size: int, /) -> bytes: ...

    def close(self) -> None: ...


class _Response:
    def __init__(self, raw: _Raw, status: int) -> None:
        self.status = status
        self._raw = raw

    def header(self, name: str) -> str | None:
        value = self._raw.headers.get(name)
        return None if value is None else str(value)

    def read(self, size: int) -> bytes:
        return self._raw.read(size)

    def close(self) -> None:
        self._raw.close()


class _NoRedirect(urllib.request.HTTPRedirectHandler):
    """Let a redirect end as an `HTTPError`, so that the download checks the target first."""

    @override
    def redirect_request(
        self,
        req: urllib.request.Request,
        fp: IO[bytes],
        code: int,
        msg: str,
        headers: Message,
        newurl: str,
    ) -> None:
        return None


def _build_opener(*handlers: urllib.request.BaseHandler) -> urllib.request.OpenerDirector:
    return urllib.request.build_opener(_NoRedirect, *handlers)


def check_https(url: str, what: str = "URL") -> None:
    """Raise `FetchError` if the URL does not have the scheme https."""
    if urlsplit(url).scheme != "https":
        message = f"{what} {url!r} is not https; expected a URL with the scheme https"
        raise FetchError(message)


def _open(url: str, opener: urllib.request.OpenerDirector) -> HttpResponse:
    check_https(url)
    headers = {"User-Agent": f"seedmatch/{version('seedmatch')}"}
    # The check above leaves only https, so urllib cannot open a file: or an ftp: URL here.
    request = urllib.request.Request(url, headers=headers)  # noqa: S310
    try:
        raw = opener.open(request, timeout=_TIMEOUT)
    except urllib.error.HTTPError as error:
        return _Response(error, error.code)
    return _Response(raw, raw.status)


def open_url(url: str) -> HttpResponse:
    """Send one GET to the URL. A redirect and an error status come back as a response."""
    return _open(url, _build_opener())
