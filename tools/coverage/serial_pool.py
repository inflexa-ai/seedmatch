"""A pytest plugin that runs the pool of seedmatch in the thread of the caller.

Use it for the coverage run, with NUMBA_DISABLE_JIT=1. The kernels then run as Python code on the
pool threads, and the coverage tracer of coverage 7.17.3 stops each thread on its data lock, so the
run does not end. A serial executor runs the same lines in the thread of the caller.

The plugin leaves `tests/test_pool.py` on real threads, because that file tests the threads.
"""

from collections.abc import Callable
from concurrent.futures import Future
from typing import Any, Self

import pytest


class _SerialExecutor:
    """The part of `ThreadPoolExecutor` that `seedmatch._pool.map_in_pool` uses."""

    def __init__(self, max_workers: int | None = None) -> None:
        self.max_workers = max_workers

    # The signature of Executor.submit takes the arguments of the submitted function.
    def submit(self, fn: Callable[..., Any], /, *args: Any, **kwargs: Any) -> Future[Any]:  # noqa: ANN401
        future: Future[Any] = Future()
        try:
            result = fn(*args, **kwargs)
        # ThreadPoolExecutor stores each BaseException in the future, so this executor does too.
        except BaseException as error:  # noqa: BLE001
            future.set_exception(error)
        else:
            future.set_result(result)
        return future

    def __enter__(self) -> Self:
        return self

    def __exit__(self, *_: object) -> None:
        return None


@pytest.fixture(autouse=True)
def _serial_pool(request: pytest.FixtureRequest, monkeypatch: pytest.MonkeyPatch) -> None:
    if request.module.__name__.rsplit(".", 1)[-1] != "test_pool":
        # The target is a text, so that seedmatch loads when a test runs, not when pytest loads this
        # plugin. Then coverage measures the import lines of seedmatch.
        monkeypatch.setattr("seedmatch._pool.ThreadPoolExecutor", _SerialExecutor)
