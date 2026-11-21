import os
import threading
import time
from concurrent.futures import ThreadPoolExecutor

import pytest

from seedmatch._pool import map_in_pool, pool_width
from seedmatch.errors import InputError


def test_map_in_pool_gives_the_results_in_the_order_of_the_items():
    def slow_first(item):
        time.sleep(0.05 if item == 0 else 0)
        return item * item

    assert list(map_in_pool(slow_first, range(20), 4)) == [i * i for i in range(20)]


@pytest.mark.parametrize("width", [1, 3, 64])
def test_map_in_pool_gives_each_result_for_any_width(width):
    assert list(map_in_pool(str, [3, 1, 2], width)) == ["3", "1", "2"]


def test_map_in_pool_of_no_item_gives_no_result():
    assert list(map_in_pool(str, [], 4)) == []


def test_map_in_pool_does_not_run_more_than_two_widths_ahead_of_the_caller():
    width = 2
    started = []
    lock = threading.Lock()

    def work(item):
        with lock:
            started.append(item)
        return item

    for taken, result in enumerate(map_in_pool(work, list(range(100)), width), start=1):
        assert result == taken - 1
        time.sleep(0.001)
        with lock:
            assert len(started) <= taken + 2 * width


def test_map_in_pool_takes_the_items_of_an_iterator_one_at_a_time():
    width = 2
    given = []

    def items():
        for item in range(100):
            given.append(item)
            yield item

    results = map_in_pool(lambda item: item, items(), width)
    assert next(results) == 0
    # The pool holds 2 * width items, and it takes one more item for the first result.
    assert len(given) <= 2 * width + 1
    assert next(results) == 1
    assert len(given) <= 2 * width + 2
    assert list(results) == list(range(2, 100))
    assert given == list(range(100))


def test_map_in_pool_takes_the_items_in_the_thread_of_the_caller():
    threads = set()

    def items():
        for item in range(20):
            threads.add(threading.get_ident())
            yield item

    assert list(map_in_pool(lambda item: item, items(), 3)) == list(range(20))
    assert threads == {threading.get_ident()}


def test_map_in_pool_takes_no_item_after_the_caller_stops():
    width = 2
    given = []

    def items():
        for item in range(100):
            given.append(item)
            yield item

    results = map_in_pool(lambda item: item, items(), width)
    next(results)
    results.close()
    assert len(given) <= 2 * width + 1


@pytest.fixture
def release():
    """Return an event that the pool sets when it shuts down, after it canceled the waiting items.

    A work item that waits for the event keeps its thread busy until the pool has canceled the
    items that wait. Thus the number of items that start does not depend on the time.
    """
    event = threading.Event()
    shutdown = ThreadPoolExecutor.shutdown

    def release_then_shut_down(self, *args, **kwargs):
        event.set()
        return shutdown(self, *args, **kwargs)

    with pytest.MonkeyPatch.context() as patch:
        patch.setattr(ThreadPoolExecutor, "shutdown", release_then_shut_down)
        yield event


def test_map_in_pool_raises_the_error_of_an_item_and_cancels_the_items_that_wait(release):
    width = 2
    started = []

    def work(item):
        started.append(item)
        if item == 0:
            message = "item 0"
            raise ValueError(message)
        assert release.wait(timeout=30)

    with pytest.raises(ValueError, match=r"^item 0$"):
        list(map_in_pool(work, list(range(100)), width))
    # Item 0 frees one thread, and each other started item keeps its thread until the shutdown.
    # Thus the threads start at most width + 1 items, and the error cancels the others.
    assert len(started) <= width + 1


def test_map_in_pool_raises_the_error_of_the_iterable_and_starts_no_new_item(release):
    started = []
    one_started = threading.Event()

    def work(item):
        started.append(item)
        if item == 1:
            one_started.set()
            assert release.wait(timeout=30)

    def items():
        yield from range(3)
        # Item 1 runs, and item 2 waits for the thread.
        assert one_started.wait(timeout=30)
        message = "no more items"
        raise KeyError(message)

    with pytest.raises(KeyError, match=r"no more items"):
        list(map_in_pool(work, items(), 1))
    assert started == [0, 1]


def test_map_in_pool_runs_the_items_on_the_threads_of_the_pool():
    main = threading.get_ident()
    names = set(map_in_pool(lambda _: threading.get_ident(), range(8), 2))
    assert main not in names


def test_pool_width_of_none_is_the_number_of_cpus():
    assert pool_width(None) == (os.cpu_count() or 1)


def test_pool_width_keeps_an_int_of_1_or_more():
    assert pool_width(3) == 3


def test_pool_width_rejects_a_width_below_1_with_the_message_of_check_at_least():
    with pytest.raises(InputError, match=r"^workers 0 is less than 1; expected 1 or more$"):
        pool_width(0)
