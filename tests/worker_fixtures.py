"""Local calculation seams for unit tests that intentionally replace kernels.

These fixtures do not test process isolation. Resource, transport, numerical
parity, and worker-failure tests use real spawned children instead.
"""

from contextlib import contextmanager
from unittest.mock import patch

import pandas as pd


@contextmanager
def calculation_body():
    """Exercise the body as a worker would, without changing native settings."""
    with patch("ldsc._runtime.worker_runtime._INITIALIZED", True):
        yield


class LocalQueryWorkers:
    """Substitute only the process boundary while testing mocked model logic."""

    def __init__(self, workspace, inputs, fit, config, workers):
        from ldsc._regression.partitioned_parallel import _stage_frame
        self.snapshot = _stage_frame(inputs.dataset.merged, workspace / "shared")
        self.inputs, self.fit, self.config = inputs, fit, config

    def map(self, tasks, *, progress=None):
        results = []
        for query, column, directory in tasks:
            with column.open() as values:
                result = self.fit(self.inputs, query, pd.DataFrame({query: values}, copy=False),
                                  directory, self.config)
            results.append(result)
            if progress is not None:
                progress.advance(outcome=result[2]["status"])
        return results

    def __enter__(self):
        return self

    def __exit__(self, *args):
        pass


def local_task_pool(*, max_workers, **kwargs):
    """Keep thread/event probes local in tests with artificial numerical kernels."""
    from concurrent.futures import ThreadPoolExecutor
    return ThreadPoolExecutor(max_workers=max_workers)
