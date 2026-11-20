from collections.abc import Mapping
from dataclasses import dataclass, field


def connection_reset():
    return ConnectionResetError("connection reset by peer")


@dataclass(frozen=True)
class Route:
    status: int = 200
    headers: Mapping[str, str] = field(default_factory=dict)
    body: bytes = b""
    reset_at: int | None = None
    error: BaseException = field(default_factory=connection_reset)


class MemoryResponse:
    def __init__(self, client: "MemoryClient", route: Route) -> None:
        self.status = route.status
        self._client = client
        self._route = route
        self._position = 0

    def header(self, name: str) -> str | None:
        wanted = name.lower()
        for key, value in self._route.headers.items():
            if key.lower() == wanted:
                return value
        return None

    def read(self, size: int) -> bytes:
        self._client.max_read_size = max(self._client.max_read_size, size)
        body, reset_at = self._route.body, self._route.reset_at
        if reset_at is not None and self._position >= reset_at:
            raise self._route.error
        limit = len(body) if reset_at is None else min(len(body), reset_at)
        end = min(self._position + size, limit)
        data = body[self._position : end]
        self._position = end
        self._client.bytes_read += len(data)
        return data

    def close(self) -> None:
        self._client.closed += 1


class MemoryClient:
    """An `open_url` function that serves the routes from memory and records what happens."""

    def __init__(self, routes: Mapping[str, Route]) -> None:
        self.routes = dict(routes)
        self.requests: list[str] = []
        self.bytes_read = 0
        self.max_read_size = 0
        self.closed = 0

    def __call__(self, url: str) -> MemoryResponse:
        self.requests.append(url)
        if url not in self.routes:
            message = f"no route for {url}"
            raise AssertionError(message)
        return MemoryResponse(self, self.routes[url])
