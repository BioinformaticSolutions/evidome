"""The one way this program runs work in parallel.

It runs on multi-processor machines, so parallelism is the default
throughout, not an optimisation added later.

Rules, inherited from Paralome ("--threads is the total CPU budget"):

- ``--threads`` is the TOTAL budget for a command, never per-task.
- CPU-bound work (decompressing and filtering, statistics, enrichment) runs in
  PROCESSES, because Python threads share one interpreter lock.
- I/O-bound work (downloads) runs in THREADS, and is capped separately so a
  large budget never turns into dozens of simultaneous connections to a public
  server.
- ``threads == 1`` runs sequentially in the calling process: deterministic,
  debuggable, and what the tests exercise alongside the parallel path.
- Results always come back in input order, whatever order tasks finish in, so
  output files never depend on scheduling.

Every parallel call site goes through ``map_ordered``. A second, hand-rolled
pool somewhere else is the "one contract, two implementations" defect. The one
exception is ``database_value._pool``: its workers need an initializer that
loads the shared design matrices once, and the spawn context so each starts
with its maths libraries held to one thread. It keeps the rest of this
contract -- ``threads == 1`` runs in the calling process.
"""

from __future__ import annotations

import os
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from typing import Callable, Iterable, List, Optional, TypeVar

T = TypeVar("T")
R = TypeVar("R")

#: Ceiling on simultaneous downloads from one public server, whatever --threads says.
MAX_DOWNLOADS = 4


def default_threads() -> int:
    """The CPUs this process may run on: the affinity mask that taskset, cgroups and SLURM set, not every CPU in
    the machine."""
    try:
        return max(1, len(os.sched_getaffinity(0)))
    except (AttributeError, OSError):     # not on macOS
        return max(1, os.cpu_count() or 1)


def resolve_threads(threads: Optional[int]) -> int:
    """An explicit --threads wins; else the CPUs SLURM allocated to this job; else the CPUs this process may use."""
    if threads is None or threads <= 0:
        allocated = os.environ.get("SLURM_CPUS_PER_TASK", "")
        if allocated.isdigit() and int(allocated) > 0:
            return int(allocated)
        return default_threads()
    return threads


def map_ordered(func: Callable[[T], R], items: Iterable[T], *, threads: int, kind: str = "process") -> List[R]:
    """``[func(x) for x in items]``, in parallel, results in input order.

    ``kind`` is "process" for CPU-bound work or "thread" for I/O-bound work.
    ``func`` must be a module-level function when ``kind == "process"``.
    The first exception raised by any task is re-raised here.
    """
    work = list(items)
    workers = min(max(1, threads), len(work)) if work else 1
    if workers == 1:
        return [func(item) for item in work]
    if kind == "thread":
        pool = ThreadPoolExecutor(max_workers=workers)
    elif kind == "process":
        pool = ProcessPoolExecutor(max_workers=workers)
    else:
        raise ValueError(f"kind must be 'process' or 'thread', not {kind!r}")
    with pool:
        return list(pool.map(func, work))
