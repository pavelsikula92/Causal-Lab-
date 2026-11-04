"""
core.stats.parallel
===================

Independent sprinklings on several CPU cores.

parallel_map(fn, arg_list, workers, progress, label) runs fn(*args) for
every tuple of arg_list on `workers` processes (concurrent.futures) and
returns the results IN THE ORDER OF arg_list, never in completion order
(tests that pair sprinklings depend on the order). Every task must be a
pure function of its arguments that seeds its own generators; then the
result is identical to a sequential run (checked by the core test
`parallel_identical`). An exception in any task stops the run: a failed
sprinkling is never silently dropped.

workers: the command line option --workers N, else the environment
variable CAUSALITY_WORKERS, else all CPU cores. workers = 1 runs in the
main process. fn must be defined at module level (pickling; on Windows
the workers import the module afresh).
"""

from __future__ import annotations

import os
from typing import Any, Callable, List, Optional


def default_workers(requested: Optional[int] = None) -> int:
    """--workers, else CAUSALITY_WORKERS, else the number of CPU cores (at least 1)."""
    if requested and requested > 0:
        return int(requested)
    try:
        w = int(os.environ.get("CAUSALITY_WORKERS", "0"))
    except ValueError:
        w = 0
    return max(1, w if w > 0 else (os.cpu_count() or 1))


def parallel_map(fn: Callable, arg_list: List[tuple], workers: int, progress: Callable[[str], None],
                 label: Callable[[int, Any, int, int], str]) -> List[Any]:
    """See the module notes. label(i, result, done, n) -> one progress line."""
    n = len(arg_list)
    w = max(1, min(int(workers), n)) if n else 1
    out: List[Any] = [None] * n
    if w == 1:
        for i, a in enumerate(arg_list):
            out[i] = fn(*a)
            progress(label(i, out[i], i + 1, n))
        return out
    for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ.setdefault(k, "1")          # one BLAS thread per process
    from concurrent.futures import ProcessPoolExecutor, as_completed
    progress(f"  parallel: {n} tasks on {w} processes")
    done = 0
    with ProcessPoolExecutor(max_workers=w) as ex:
        futs = {ex.submit(fn, *a): i for i, a in enumerate(arg_list)}
        for f in as_completed(futs):
            i = futs[f]
            out[i] = f.result()
            done += 1
            progress(label(i, out[i], done, n))
    return out
