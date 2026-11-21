"""The thread pool of the functions that run many independent searches."""

import os
from collections import deque
from collections.abc import Callable, Iterable, Iterator
from concurrent.futures import Future, ThreadPoolExecutor

from seedmatch._checks import as_limit, check_at_least


def pool_width(workers: int | None) -> int:
    """Return the number of threads: `workers`, or the number of CPUs for `None`."""
    if workers is None:
        # os.cpu_count() is None if the number of CPUs is unknown.
        return os.cpu_count() or 1
    width = as_limit("workers", workers)
    check_at_least("workers", width, 1)
    return width


def map_in_pool[S, T](work: Callable[[S], T], items: Iterable[S], width: int) -> Iterator[T]:
    """Yield the result of `work` for each item, in the order of the items.

    A result comes while the pool works on the next items, thus a caller can write each result into
    an array at once. The function takes the items one at a time from `items`, in the thread of the
    caller. After the caller took k results, it has taken at most k + 2 `width` items. Thus the
    items and the results that wait take O(`width`) memory, also if `items` is long or the caller
    is slow. If `work` raises, the error reaches the caller and the pool starts no new item.
    """
    with ThreadPoolExecutor(max_workers=width) as pool:
        pending: deque[Future[T]] = deque()
        try:
            for item in items:
                pending.append(pool.submit(work, item))
                if len(pending) > 2 * width:
                    yield pending.popleft().result()
            while pending:
                yield pending.popleft().result()
        finally:
            for future in pending:
                future.cancel()
