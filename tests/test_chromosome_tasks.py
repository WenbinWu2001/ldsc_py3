"""Bounded chromosome scheduling without numerical or executor ownership."""

from concurrent.futures import Future
from concurrent.futures.process import BrokenProcessPool

import pytest

from ldsc._runtime import parallelism as _parallelism


@pytest.mark.parametrize("workers", [1, 2, 4])
def test_bounds_outstanding_work_and_collects_ready_group_before_refill(workers):
    submitted, collected = [], []
    outstanding = set()

    def submit(chrom):
        # Every submitted future is ready, so refill must drain the whole group.
        if len(submitted) >= workers:
            assert not outstanding or len(submitted) % workers != 0
        assert len(outstanding) < workers
        submitted.append(chrom)
        outstanding.add(chrom)
        future = Future()
        future.set_result(f"result-{chrom}")
        return future

    def collect(chrom, result):
        assert result == f"result-{chrom}"
        outstanding.remove(chrom)
        collected.append(chrom)

    def failed(*args):
        pytest.fail(f"Unexpected scheduling failure: {args}")

    _parallelism._run_chromosome_tasks(
        ["1", "2", "3", "4", "5", "6", "7"], workers, submit, collect, failed)
    assert submitted == ["1", "2", "3", "4", "5", "6", "7"]
    assert sorted(collected) == submitted
    assert not outstanding


class CancellableFuture(Future):
    def __init__(self):
        super().__init__()
        self.cancel_attempts = 0

    def cancel(self):
        self.cancel_attempts += 1
        return super().cancel()


@pytest.mark.parametrize("error_type", [ValueError, BrokenProcessPool])
@pytest.mark.parametrize("first_ready", ["1", "2"])
def test_ready_failure_prevents_refill_and_preserves_context(monkeypatch, error_type, first_ready):
    error = error_type("original task failure")
    futures = {chrom: CancellableFuture() for chrom in ["1", "2", "3", "4"]}
    futures["1"].set_result("ok")
    futures["2"].set_exception(error)
    futures["3"].set_running_or_notify_cancel()
    submitted, collected, failures = [], [], []
    # Exercise both legal iteration orders for two futures ready at once.
    order = ["1", "2"] if first_ready == "1" else ["2", "1"]
    monkeypatch.setattr(_parallelism, "wait", lambda *args, **kwargs:
                        ([futures[chrom] for chrom in order], {futures["3"], futures["4"]}))

    def submit(chrom):
        submitted.append(chrom)
        return futures[chrom]

    def failed(original, chrom, outstanding):
        # Cancellation precedes caller diagnostics/shutdown; running work survives.
        assert futures["3"].cancel_attempts == futures["4"].cancel_attempts == 1
        assert futures["3"].running()
        assert futures["4"].cancelled()
        failures.append((original, chrom, outstanding))

    with pytest.raises(error_type) as raised:
        _parallelism._run_chromosome_tasks(
            ["1", "2", "3", "4", "5"], 4, submit,
            lambda chrom, result: collected.append((chrom, result)), failed)
    assert raised.value is error
    assert submitted == ["1", "2", "3", "4"]
    assert collected == ([("1", "ok")] if first_ready == "1" else [])
    assert failures == [(error, "2", ("2", "3", "4") if first_ready == "1" else ("1", "2", "3", "4"))]


@pytest.mark.parametrize("failure_stage", ["initial-submit", "refill-submit", "collect", "wait"])
def test_coordinator_failure_preserves_exception_and_task_context(monkeypatch, failure_stage):
    error = OSError("original coordinator failure")
    first, second = CancellableFuture(), CancellableFuture()
    first.set_result("first")
    submitted, failures = [], []

    def submit(chrom):
        submitted.append(chrom)
        if ((failure_stage == "initial-submit" and chrom == "2")
                or (failure_stage == "refill-submit" and chrom == "3")):
            raise error
        return first if chrom == "1" else second

    def collect(chrom, result):
        if failure_stage == "collect":
            raise error

    if failure_stage == "wait":
        def broken_wait(*args, **kwargs):
            raise error
        monkeypatch.setattr(_parallelism, "wait", broken_wait)

    with pytest.raises(OSError) as raised:
        _parallelism._run_chromosome_tasks(
            ["1", "2", "3", "4"], 2, submit, collect,
            lambda *context: failures.append(context))
    assert raised.value is error
    expected = {
        "initial-submit": ("2", ("1",)),
        "refill-submit": ("3", ("2",)),
        "collect": ("1", ("1", "2")),
        "wait": (None, ("1", "2")),
    }
    assert failures == [(error, *expected[failure_stage])]
    assert submitted == (["1", "2", "3"] if failure_stage == "refill-submit" else ["1", "2"])
    for chrom in expected[failure_stage][1]:
        assert (first if chrom == "1" else second).cancel_attempts == 1


@pytest.mark.parametrize("error_type", [KeyboardInterrupt, SystemExit])
def test_interruptions_do_not_gain_task_failure_handling(error_type):
    error = error_type("interrupted")
    future = CancellableFuture()
    future.set_exception(error)
    failures = []
    with pytest.raises(error_type) as raised:
        _parallelism._run_chromosome_tasks(
            ["1", "2"], 1, lambda chrom: future, lambda *args: None,
            lambda *args: failures.append(args))
    assert raised.value is error
    assert not failures
    assert future.cancel_attempts == 0


def test_empty_chromosome_sequence_does_not_submit():
    def unexpected(*args):
        pytest.fail(f"Unexpected callback for no chromosomes: {args}")
    _parallelism._run_chromosome_tasks([], 1, unexpected, unexpected, unexpected)


def test_ready_chromosome_releases_slot_while_another_remains_running():
    first = Future()
    collected, submitted = [], []

    def submit(chrom):
        submitted.append(chrom)
        if chrom == "1":
            return first
        if chrom == "3":
            assert collected == [("2", "second")]
            first.set_result("first")
        future = Future()
        future.set_result("second" if chrom == "2" else "third")
        return future

    def failed(*args):
        pytest.fail(f"Unexpected scheduling failure: {args}")

    _parallelism._run_chromosome_tasks(
        ["1", "2", "3"], 2, submit, lambda *result: collected.append(result), failed)
    assert submitted == ["1", "2", "3"]
    assert sorted(collected) == [("1", "first"), ("2", "second"), ("3", "third")]
