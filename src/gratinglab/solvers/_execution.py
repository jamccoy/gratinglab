"""How a rigorous solver spends the machine: BLAS threads, and workers.

Two things, both measured before they were written (``benchmarks/
solver_speed.py``; the figures are in the changelog).

**BLAS threads.** The rigorous solvers' time is in LAPACK, on matrices of a
few hundred to a few thousand rows. OpenBLAS's default is a thread per core,
and at these sizes the threads spend more time synchronising than computing;
on a machine already busy with other work they also wait for a core each, and
a 122-square inverse that takes microseconds took a tenth of a second. So
:func:`blas_threads` caps the BLAS pools for the length of a solve and puts
them back afterwards. The cap is process-wide, since that is what a BLAS pool
is, so it is reference-counted: two solves overlapping on different threads
share one cap, and the original setting returns when the last one leaves.

**Workers.** The wavelengths of a scan are independent, and LAPACK releases
the GIL, so ``workers`` threads each solving whole wavelengths on one BLAS
thread use the cores far better than one wavelength on many.
:func:`run_rows` does that while keeping the progress contract of
:class:`~gratinglab.solvers.base.Solver` exactly: every callback is made on
the caller's thread, ``(0, total)`` comes first, ``done`` counts finished
rows in order, and an exception from the callback (``SolveCancelled``) is
re-raised unchanged after the rows already running have finished, so no
worker outlives the call and no partial scan is returned.
"""

from __future__ import annotations

import threading
from numbers import Integral
from collections.abc import Callable, Iterator
from concurrent.futures import FIRST_COMPLETED, Future, ThreadPoolExecutor, wait
from contextlib import contextmanager
from typing import TypeVar

from threadpoolctl import threadpool_limits

from .base import Progress

__all__ = ["capped_blas", "check_workers", "check_blas_threads", "run_rows"]

T = TypeVar("T")
S = TypeVar("S")

_lock = threading.Lock()
_depth = 0
_limiter: threadpool_limits | None = None


def check_workers(workers) -> int:
    """``workers`` as a positive int, or ``ValueError``."""
    if isinstance(workers, bool) or not isinstance(workers, Integral) or workers < 1:
        raise ValueError(f"workers must be a positive integer, got {workers!r}")
    return int(workers)


def check_blas_threads(threads) -> int | None:
    """``blas_threads`` as a positive int or ``None``, or ``ValueError``."""
    if threads is None:
        return None
    if isinstance(threads, bool) or not isinstance(threads, Integral) or threads < 1:
        raise ValueError(
            f"blas_threads must be a positive integer or None, got {threads!r}"
        )
    return int(threads)


@contextmanager
def capped_blas(threads: int | None) -> Iterator[None]:
    """Cap every loaded BLAS pool at ``threads`` while the block runs.

    ``None`` leaves the pools alone. Overlapping blocks share the first
    one's cap; the original limits come back when the last block exits,
    however it exits.
    """
    global _depth, _limiter
    if threads is None:
        yield
        return
    with _lock:
        if _depth == 0:
            _limiter = threadpool_limits(limits=threads, user_api="blas")
        _depth += 1
    try:
        yield
    finally:
        with _lock:
            _depth -= 1
            if _depth == 0:
                assert _limiter is not None
                _limiter.restore_original_limits()
                _limiter = None


def run_rows(
    count: int,
    row: Callable[[int, S], T],
    *,
    state: Callable[[], S],
    workers: int = 1,
    progress: Progress | None = None,
) -> list[T]:
    """``[row(i, s) for i in range(count)]``, on ``workers`` threads.

    ``state()`` is called once per thread that runs rows, and its result is
    passed to every row that thread runs: per-thread scratch such as the
    integral solver's kernel memo. The result list is in row order whatever
    order the rows finished in. With ``workers=1`` nothing is threaded, and
    progress is called exactly where the serial loops always called it.
    """
    if workers == 1 or count <= 1:
        local = state()
        results = []
        for index in range(count):
            if progress is not None:
                progress(index, count)
            results.append(row(index, local))
        if progress is not None:
            progress(count, count)
        return results

    if progress is not None:
        progress(0, count)  # a cancellation point before any work
    slot = threading.local()

    def task(index: int) -> T:
        if not hasattr(slot, "state"):
            slot.state = state()
        return row(index, slot.state)

    results: list = [None] * count
    pool = ThreadPoolExecutor(
        max_workers=min(workers, count), thread_name_prefix="gratinglab-row"
    )
    try:
        pending: dict[Future, int] = {
            pool.submit(task, index): index for index in range(count)
        }
        done = 0
        while pending:
            finished, _ = wait(pending, return_when=FIRST_COMPLETED)
            # Several can finish between two wakes; report them one at a
            # time so `done` steps by one, as the serial loop's does.
            for future in finished:
                results[pending.pop(future)] = future.result()
                done += 1
                if progress is not None:
                    progress(done, count)
    finally:
        # On the way out by an exception, drop the rows not yet started and
        # wait for the running ones: nothing computes after the call returns.
        pool.shutdown(wait=True, cancel_futures=True)
    return results
