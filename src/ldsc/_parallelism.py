"""Shared worker-count validation and resolution for LDSC workflows.

Negative requests use process CPU affinity with a machine CPU-count fallback.
Native numerical-library thread settings are owned by the workflow workers.
"""

import os


def _validate_threads(threads: int, *, error_type: type[Exception] = ValueError) -> int:
    """Require a nonzero Python integer, excluding booleans, without coercion."""
    if isinstance(threads, bool) or not isinstance(threads, int) or threads == 0:
        raise error_type("threads must be a nonzero integer (1, positive workers, or a negative CPU offset).")
    return threads


def _parse_threads(value: str) -> int:
    """Parse CLI integer syntax and apply the shared worker-request validator."""
    return _validate_threads(int(value))


def _available_cpu_count() -> int:
    """Return the number of CPUs available to this process.

    Prefer CPU affinity when supported; otherwise use ``os.cpu_count()``.
    Scheduler allocations are respected when reflected in process affinity.
    """
    getaffinity = getattr(os, "sched_getaffinity", None)
    if getaffinity is not None:
        try:
            return len(getaffinity(0)) or 1
        except OSError:
            pass
    return os.cpu_count() or 1


def _resolve_worker_count(threads: int, n_chromosomes: int) -> int:
    """Resolve the effective worker count from the ``threads`` setting.

    Uses the joblib ``n_jobs`` convention: ``1`` is sequential, a positive ``N``
    requests ``N`` workers, ``-1`` requests all available cores, ``-2`` all but
    one, and any negative ``-k`` requests ``n_cpus + 1 - k``. Core counts respect
    CPU affinity (see :func:`_available_cpu_count`). The result is capped at
    ``n_chromosomes`` and floored at ``1`` so a single work item never spawns a
    pool. Query workflows pass their bounded query count as the work limit.
    Positive requests are capped only by work; scheduler allocations influence
    negative requests when reflected in CPU affinity.

    For example, eight available CPUs and three work items resolve requests
    ``4`` and ``-1`` to three workers; ``-7`` resolves to two. A positive request
    is not capped by CPU count, and ``SLURM_CPUS_PER_TASK`` is not read here.
    """
    _validate_threads(threads)
    if n_chromosomes <= 0:
        return 1
    if threads < 0:
        resolved = _available_cpu_count() + 1 + threads
    else:
        resolved = threads
    return max(1, min(resolved, n_chromosomes))
