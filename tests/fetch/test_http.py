import http.client
import io
import urllib.request
import urllib.response
from importlib.metadata import version

import pytest

from seedmatch.errors import FetchError
from seedmatch.fetch.http import _build_opener, _open, open_url

URL = "https://example.org/data/file.bin"
TARGET = "https://cdn.example.org/store/file.bin"


class CannedHttps(urllib.request.BaseHandler):
    """An `https` handler that answers from memory. The order puts it before the default handler."""

    handler_order = 100

    def __init__(self, routes):
        self.routes = routes
        self.requests = []
        self.bodies = []

    def https_open(self, req):
        self.requests.append(req)
        status, headers, body = self.routes[req.full_url]
        message = http.client.HTTPMessage()
        for name, value in headers.items():
            message[name] = value
        stream = io.BytesIO(body)
        self.bodies.append(stream)
        response = urllib.response.addinfourl(stream, message, req.full_url, status)
        response.msg = http.client.responses.get(status, "")
        return response


def open_with(routes, url=URL):
    canned = CannedHttps(routes)
    return _open(url, _build_opener(canned)), canned


def test_a_200_gives_the_status_the_headers_and_the_body():
    response, canned = open_with({URL: (200, {"Content-Length": "6"}, b"abcdef")})
    try:
        assert response.status == 200
        assert response.header("content-length") == "6"
        assert response.header("ETag") is None
        assert [response.read(4), response.read(4), response.read(4)] == [b"abcd", b"ef", b""]
        assert [body.closed for body in canned.bodies] == [False]
    finally:
        response.close()
    assert [request.full_url for request in canned.requests] == [URL]
    assert [body.closed for body in canned.bodies] == [True]


@pytest.mark.parametrize("status", [301, 302, 303, 307, 308])
def test_a_redirect_is_not_followed(status):
    routes = {URL: (status, {"Location": TARGET}, b"moved"), TARGET: (200, {}, b"target")}
    response, canned = open_with(routes)
    try:
        assert response.status == status
        assert response.header("Location") == TARGET
    finally:
        response.close()
    assert [request.full_url for request in canned.requests] == [URL]
    assert [body.closed for body in canned.bodies] == [True]


@pytest.mark.parametrize("status", [404, 500])
def test_an_error_status_gives_the_status_and_the_body(status):
    response, canned = open_with({URL: (status, {"Content-Type": "text/plain"}, b"not found")})
    try:
        assert response.status == status
        assert response.header("content-type") == "text/plain"
        assert response.read(100) == b"not found"
    finally:
        response.close()
    assert [body.closed for body in canned.bodies] == [True]


def test_the_request_sends_the_user_agent_a_timeout_and_a_get():
    response, canned = open_with({URL: (200, {}, b"")})
    response.close()
    (request,) = canned.requests
    assert request.get_header("User-agent") == f"seedmatch/{version('seedmatch')}"
    assert request.timeout == 60
    assert request.get_method() == "GET"
    assert request.data is None


def test_open_url_uses_the_default_opener(monkeypatch):
    canned = CannedHttps({URL: (200, {}, b"abc")})
    opener = _build_opener(canned)
    monkeypatch.setattr("seedmatch.fetch.http._build_opener", lambda: opener)
    response = open_url(URL)
    try:
        assert (response.status, response.read(10)) == (200, b"abc")
    finally:
        response.close()


@pytest.mark.parametrize(
    "url", ["file:///etc/passwd", "ftp://example.org/x", "http://example.org/x"]
)
def test_open_url_refuses_a_url_that_is_not_https(monkeypatch, url):
    canned = CannedHttps({url: (200, {}, b"abc")})
    opener = _build_opener(canned)
    monkeypatch.setattr("seedmatch.fetch.http._build_opener", lambda: opener)
    message = f"URL {url!r} is not https; expected a URL with the scheme https"
    with pytest.raises(FetchError) as caught:
        open_url(url)
    assert str(caught.value) == message
    assert canned.requests == []
