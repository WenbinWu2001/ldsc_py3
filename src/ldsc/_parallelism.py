"""Worker budgets and bounded chromosome scheduling for LDSC workflows.

Requests use joblib.cpu_count() directly as the detected available CPU budget.
Native numerical-library thread settings are owned by the workflow workers.
"""

import logging
import argparse
from concurrent.futures import FIRST_COMPLETED, wait
import joblib

LOGGER = logging.getLogger("LDSC.parallelism")


def _run_chromosome_tasks(chromosomes, workers, submit, collect, on_failure):
    """Collect chromosome tasks with at most ``workers`` outstanding futures.

    ``workers`` is the already-resolved positive count. ``submit(chrom)``
    prepares/submits a task and returns its future; ``collect(chrom, result)``
    processes each result. Drain each ready group before submitting more work.
    On an ``Exception``, attempt cancellation, then call
    ``on_failure(error, current_chrom, outstanding_chromosomes)`` and reraise.
    The current chromosome identifies the operation that raised, not the cause
    of a broken pool. The callback may translate the error with explicit cause.

    The caller owns the executor, diagnostics, and scratch. Cancellation does
    not stop dispatched/running work: the caller must finish normal executor
    shutdown before releasing temporary inputs or outputs. BaseException
    interruptions retain the caller's existing handling.
    """
    remaining, pending = iter(chromosomes), {}
    current_chrom = None
    try:
        while True:
            while len(pending) < workers:
                current_chrom = next(remaining, None)
                if current_chrom is None:
                    break
                pending[submit(current_chrom)] = current_chrom
                current_chrom = None
            if not pending:
                return
            ready, _ = wait(pending, return_when=FIRST_COMPLETED)
            for future in ready:
                current_chrom = pending[future]
                collect(current_chrom, future.result())
                del pending[future]
                current_chrom = None
    except Exception as error:
        for future in pending:
            future.cancel()
        on_failure(error, current_chrom, tuple(pending.values()))
        raise


def _validate_workers(workers: int, *, error_type: type[Exception] = ValueError) -> int:
    """Require a nonzero Python integer, excluding booleans, without coercion."""
    if isinstance(workers, bool) or not isinstance(workers, int) or workers == 0:
        raise error_type("workers must be a nonzero integer (1, positive workers, or a negative CPU offset).")
    budget = joblib.cpu_count()
    if workers > budget:
        raise error_type(
            f"Requested workers={workers} exceeds the detected available CPU count={budget} "
            "from joblib.cpu_count(). Choose --workers at or below that count and your job allocation."
        )
    return workers


def _parse_workers(value: str) -> int:
    """Parse CLI integer syntax and apply the shared worker-request validator."""
    try:
        return _validate_workers(int(value))
    except ValueError as error:
        raise argparse.ArgumentTypeError(str(error)) from error


def _available_cpu_count() -> int:
    """Return joblib's detected available CPU count, not idle or allocated CPUs."""
    return joblib.cpu_count()


def _resolve_worker_count(workers: int, task_count: int, *, task_unit: str = "chromosome",
                          batch_size: int | None = None) -> int:
    """Validate and report task concurrency, capped by work after the CPU check.

    Negative ``-k`` requests use ``max(1, joblib.cpu_count() + 1 - k)``.
    An explicit positive request above detected availability is an error even
    when fewer tasks exist. One task uses one managed child, not inline BLAS.
    Scheduler allocations not reflected in supported process restrictions
    remain the user's responsibility; this does not measure idle CPUs.
    """
    _validate_workers(workers)
    budget = _available_cpu_count()
    if workers < 0:
        resolved = budget + 1 + workers
    else:
        resolved = workers
    available_tasks = min(task_count, batch_size) if batch_size is not None else task_count
    resolved = max(1, min(resolved, available_tasks))
    LOGGER.info(
        "Compute workers: requested=%d; detected available CPU count=%d (joblib.cpu_count()); "
        "task unit=%s; task count=%d; batch cap=%s; worker count=%d. "
        "Each worker uses 1 BLAS thread; worker count is not total software threads or guaranteed CPU utilization.",
        workers, budget, task_unit, task_count, batch_size, resolved,
    )
    return resolved


WORKERS_HELP = (
    "Default: 1, one managed child process. Positive N requests at most N workers; "
    "-1 uses the detected available CPU count from joblib.cpu_count(); -2 uses one fewer "
    "(minimum one). Reject zero and positive requests above the detected count, then cap "
    "by available tasks. Each worker uses 1 BLAS thread, in serial and parallel execution. "
    "Detection does not measure idle CPUs or read your scheduler allocation; choose a value "
    "consistent with that allocation. More workers can use more memory."
)
