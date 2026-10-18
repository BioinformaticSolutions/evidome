"""map_ordered returns input order, for threads, processes and sequential runs."""

import time

import pytest

from evidome.parallel import map_ordered, resolve_threads


def _slow_square(x):
    time.sleep(0.01 * (5 - x % 5))
    return x * x


@pytest.mark.parametrize("kind", ["thread", "process"])
@pytest.mark.parametrize("threads", [1, 4])
def test_results_in_input_order(kind, threads):
    assert map_ordered(_slow_square, range(10), threads=threads, kind=kind) == [x * x for x in range(10)]


def _boom(x):
    if x == 3:
        raise ValueError("task 3 failed")
    return x


def test_task_failure_surfaces():
    with pytest.raises(ValueError, match="task 3"):
        map_ordered(_boom, range(6), threads=3, kind="process")


def test_zero_means_all_cpus():
    assert resolve_threads(0) >= 1 and resolve_threads(5) == 5


def test_empty_input():
    assert map_ordered(_slow_square, [], threads=4) == []
