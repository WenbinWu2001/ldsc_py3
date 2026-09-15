"""Worker policy, immutable transport, and real spawned-worker failure recovery."""

import json
import multiprocessing as mp
import os
import sys

import numpy as np
import pandas as pd
import pytest

from ldsc import RegressionRunner, RegressionConfig, cli
from ldsc import _parallelism as policy
from ldsc import _partitioned_h2_parallel as parallel
from ldsc.errors import LDSCInternalError
from tests.test_regression_streaming import batch_inputs, cli_inputs
from tests.test_partitioned_query_failures import concentrate_queries


@pytest.mark.parametrize("inputs", [
    ["partitioned-h2", "--sumstats-file", "trait", "--ldscore-dir", "ld"],
    ["ldscore"],
    ["build-gene-ldscore-index", "--baseline-annot-sources", "baseline", "--plink-prefix", "panel",
     "--gene-coordinate-file", "genes", "--snp-identifier", "rsid", "--genome-build", "hg19"],
])
def test_cli_thread_contract(inputs, capsys):
    inputs = [*inputs, "--output-dir", "out"]
    parser = cli.build_parser()
    assert parser.parse_args(inputs).threads == 1
    for value in ("1", "4", "-1", "-2", "-20"):
        assert parser.parse_args([*inputs, "--threads", value]).threads == int(value)
    for value in ("0", "True", "False", "1.5", "1.0"):
        with pytest.raises(SystemExit) as error:
            parser.parse_args([*inputs, "--threads", value])
        assert error.value.code == 2
        assert "argument --threads:" in capsys.readouterr().err


@pytest.mark.parametrize("batch_size,expected", [(1000, 3), (2, 2), (1, 1)])
def test_positive_query_workers_follow_existing_work_cap(tmp_path, monkeypatch, batch_size, expected):
    monkeypatch.setattr(policy.os, "sched_getaffinity", lambda pid: {0}, raising=False)
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "1")
    table, source = batch_inputs(tmp_path / "ld")
    RegressionRunner(source.config_snapshot, RegressionConfig(n_blocks=6)).estimate_partitioned_h2_batch(
        table, source, output_dir=tmp_path / "out", threads=4, query_batch_size=batch_size,
    )
    metadata = json.loads((tmp_path / "out/diagnostics/metadata.json").read_text())
    assert metadata["query_workers_effective"] == expected


@pytest.mark.parametrize("reason", ["requested_one", "one_query", "one_batch", "one_cpu", "baseline"])
def test_effective_one_stays_inline_and_preserves_native_settings(tmp_path, monkeypatch, reason):
    from threadpoolctl import threadpool_info

    table, source = batch_inputs(tmp_path / "ld")
    before = threadpool_info()
    monkeypatch.setattr(parallel, "_QueryWorkers", lambda *a, **k: pytest.fail("inline run spawned workers"))
    monkeypatch.setattr(policy.os, "sched_getaffinity", lambda pid: {0}, raising=False)
    threads = 1 if reason == "requested_one" else -1 if reason == "one_cpu" else 4
    result = RegressionRunner(source.config_snapshot, RegressionConfig(n_blocks=6)).estimate_partitioned_h2_batch(
        table, source, output_dir=tmp_path / "out", threads=threads,
        query_columns=[] if reason == "baseline" else ["first"] if reason == "one_query" else None,
        query_batch_size=1 if reason == "one_batch" else 1000,
    )
    assert threadpool_info() == before
    assert len(result.summary) == (1 if reason in {"baseline", "one_query"} else 3)
    assert not list((tmp_path / "out").glob(".ldsc-annotation-*"))


def test_numeric_transport_preserves_dtype_order_and_is_read_only(tmp_path):
    frame = pd.DataFrame({"q": np.array([2., 1., 3.], dtype=np.float32), "N": [10, 20, 30],
                          "SNP": ["rs3", "rs1", "rs2"], "nullable": pd.array([1, None, 3], dtype="Int64")})
    snapshot = parallel._stage_frame(frame, tmp_path / "mapped")
    with snapshot.open() as actual:
        pd.testing.assert_frame_equal(actual, frame)
        assert not actual.q.to_numpy(copy=False).flags.writeable
        with pytest.raises(ValueError, match="read-only"):
            actual.q.to_numpy(copy=False)[0] = 99
    del actual
    with snapshot.numeric["q"].open() as column:
        np.testing.assert_array_equal(column, [2., 1., 3.])
        assert column.dtype == np.float32


def test_nonnumeric_query_keeps_per_query_preparation_failure_policy(tmp_path):
    table, source = batch_inputs(tmp_path / "ld")
    path = tmp_path / "ld/ldscore.query.parquet"
    query = pd.read_parquet(path)
    query["second"] = "invalid query values"
    query.to_parquet(path, index=False)
    runner = RegressionRunner(source.config_snapshot, RegressionConfig(n_blocks=6))
    expected, actual = [runner.estimate_partitioned_h2_batch(
        table, source, output_dir=tmp_path / f"workers-{threads}", threads=threads, continue_on_query_error=True,
    ) for threads in (1, 2)]
    pd.testing.assert_frame_equal(actual.summary, expected.summary)
    pd.testing.assert_frame_equal(actual.query_status, expected.query_status)
    assert actual.query_status.status.tolist() == ["success", "failed", "success"]
    assert actual.query_status.loc[1, "stage"] == "model_preparation"


def _fault_worker(connection, snapshot, inputs_path):
    """Inject failures inside real spawned processes at OS/numerical boundaries."""
    from ldsc._partitioned_h2_parallel import _query_worker
    from threadpoolctl import threadpool_info

    mode = os.environ["LDSC_TEST_WORKER_FAULT"]
    original_solve = np.linalg.solve
    original_replace = os.replace

    def solve(matrix, rhs):
        assert all(item["num_threads"] == 1 for item in threadpool_info())
        if sys.platform == "darwin":
            import ctypes
            blas = ctypes.CDLL("/System/Library/Frameworks/Accelerate.framework/Accelerate")
            if hasattr(blas, "BLASGetThreading"):
                assert blas.BLASGetThreading() == 1
        if matrix.shape == (3, 3) and not matrix[1].any():
            if mode == "crash":
                os._exit(17)
            if mode == "interrupt":
                raise KeyboardInterrupt("worker interruption")
            if mode == "exception":
                raise RuntimeError("injected numerical backend failure")
        return original_solve(matrix, rhs)

    def replace(source, destination):
        if mode == "disk" and "fit-1" in str(destination):
            raise OSError("injected worker output device failure")
        return original_replace(source, destination)

    np.linalg.solve = solve
    os.replace = replace
    _query_worker(connection, snapshot, inputs_path)


@pytest.mark.parametrize("continue_errors", [False, True])
@pytest.mark.parametrize("mode", ["exception", "crash", "interrupt", "disk"])
def test_real_worker_failures_identify_query_and_reap_children(tmp_path, monkeypatch, continue_errors, mode):
    args, _ = cli_inputs(tmp_path)
    if mode != "disk":
        concentrate_queries(tmp_path)
    for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
        monkeypatch.setenv(key, "7")
    monkeypatch.setenv("LDSC_TEST_WORKER_FAULT", mode)
    monkeypatch.setattr(parallel, "_query_worker", _fault_worker)
    before = {child.pid for child in mp.active_children()}
    args.extend(["--threads", "2"])
    if continue_errors:
        args.append("--continue-on-query-error")
    if mode == "exception" and continue_errors:
        result = cli.main(args)
        assert set(result.category) == {"first", "third"}
    else:
        error_type = {"exception": RuntimeError, "crash": LDSCInternalError,
                      "interrupt": KeyboardInterrupt, "disk": OSError}[mode]
        with pytest.raises(error_type) as caught:
            cli.main(args)
        assert "second" in str(caught.value) + " ".join(getattr(caught.value, "__notes__", []))
        assert not (tmp_path / "out/partitioned_h2.tsv").exists()
    assert {child.pid for child in mp.active_children()} == before
    assert not list((tmp_path / "out").glob(".ldsc-annotation-*"))
    if mode == "exception":
        status = pd.read_csv(tmp_path / "out/diagnostics/query_status.tsv", sep="\t")
        assert status.status.tolist() == ["success", "failed", "success"]
        assert status.loc[1, "error_type"] == "RuntimeError"


def test_parent_interrupt_reaps_workers_before_scratch_cleanup(tmp_path, monkeypatch):
    args, _ = cli_inputs(tmp_path)
    before = {child.pid for child in mp.active_children()}

    def stop(*args):
        raise KeyboardInterrupt("parent interruption")

    monkeypatch.setattr(parallel, "wait", stop)
    with pytest.raises(KeyboardInterrupt, match="parent interruption"):
        cli.main([*args, "--threads", "2", "--continue-on-query-error"])
    assert {child.pid for child in mp.active_children()} == before
    assert not list((tmp_path / "out").glob(".ldsc-annotation-*"))
    assert not (tmp_path / "out/partitioned_h2.tsv").exists()
