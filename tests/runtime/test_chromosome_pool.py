"""A crashed worker becomes one actionable error in every chromosome workflow."""

from concurrent.futures import wait
from concurrent.futures.process import BrokenProcessPool
import os

import pytest

from ldsc._runtime.parallelism import chromosome_pool
from ldsc.errors import LDSCInternalError


@pytest.mark.parametrize("failure_stage", ["result", "submit"])
def test_worker_crash_becomes_internal_error_naming_workers(failure_stage):
    with pytest.raises(LDSCInternalError) as raised:
        with chromosome_pool(1, operation="ldscore chromosome computation") as pool:
            crashed = pool.submit(os._exit, 17)
            if failure_stage == "result":
                crashed.result()
            wait([crashed])
            pool.submit(abs, -1)
    assert isinstance(raised.value.__cause__, BrokenProcessPool)
    message = str(raised.value)
    assert message.startswith("ldscore chromosome computation aborted")
    assert "`--workers`" in message


def test_task_errors_other_than_crashes_propagate_unchanged():
    with pytest.raises(ValueError, match="task failure"):
        with chromosome_pool(1, operation="ldscore chromosome computation") as pool:
            pool.submit(int, "task failure").result()
