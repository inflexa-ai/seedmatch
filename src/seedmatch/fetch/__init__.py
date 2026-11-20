"""The download of the registry resources. Only this package uses the network."""

from seedmatch.errors import FetchError
from seedmatch.fetch.download import DEFAULT_FETCH_OPTS, FetchOpts, fetch_resource
from seedmatch.fetch.http import HttpResponse, open_url

__all__ = [
    "DEFAULT_FETCH_OPTS",
    "FetchError",
    "FetchOpts",
    "HttpResponse",
    "fetch_resource",
    "open_url",
]
