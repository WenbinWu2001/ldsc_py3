"""Observable worker budgets and isolated native-library settings."""

import os
from concurrent.futures import ProcessPoolExecutor
import multiprocessing as mp

import pytest

from ldsc._runtime import parallelism as policy


@pytest.mark.parametrize("requested,tasks,expected", [(1,22,1),(4,2,2),(-1,22,8),(-2,22,7),(-20,22,1)])
def test_detected_budget_and_task_cap(monkeypatch, requested, tasks, expected, caplog):
    monkeypatch.setattr("joblib.cpu_count", lambda: 8)
    with caplog.at_level("INFO", logger="LDSC"):
        assert policy._resolve_worker_count(requested, tasks) == expected
    assert "joblib.cpu_count()" in caplog.text
    assert "detected available CPU count=8" in caplog.text


def test_explicit_over_budget_request_is_rejected_before_task_cap(monkeypatch):
    monkeypatch.setattr("joblib.cpu_count", lambda: 2)
    with pytest.raises(ValueError, match="detected available CPU count.*2"):
        policy._resolve_worker_count(4, 1)


def test_unverified_native_limits_are_recorded_without_failing(monkeypatch, caplog):
    from ldsc._runtime import worker_runtime as runtime
    monkeypatch.setattr("threadpoolctl.threadpool_info", lambda: [])
    monkeypatch.setattr(runtime, "_ACCELERATE", None)
    monkeypatch.setattr(runtime, "_BUILD_BACKENDS", {"numpy": "unsupported", "scipy": "unsupported"})
    with caplog.at_level("WARNING", logger="LDSC"):
        runtime.log_native_thread_state()
    assert "BLAS target=1" in caplog.text
    assert "unverified or unmet (continuing)" in caplog.text
    assert "unsupported" in caplog.text


def _inspect_child():
    from ldsc._runtime.worker_runtime import initialize_worker, native_thread_state
    initialize_worker()
    return os.getpid(), native_thread_state()


def test_native_limits_are_set_in_child_without_changing_parent(monkeypatch):
    import pyarrow as pa
    from threadpoolctl import threadpool_info
    monkeypatch.setenv("OPENBLAS_NUM_THREADS", "7")
    before = (threadpool_info(), pa.cpu_count(), pa.io_thread_count())
    with ProcessPoolExecutor(max_workers=1, mp_context=mp.get_context("spawn")) as pool:
        pid, state = pool.submit(_inspect_child).result()
        assert pool.submit(_inspect_child).result()[0] == pid
    assert pid != os.getpid()
    assert state["blas_target"] == 1
    assert state["openmp_target"] == 1
    assert state["arrow_cpu"] == state["arrow_io"] == 1
    assert all(item["num_threads"] == 1 for item in state["pools"])
    if state["accelerate"] is not None:
        assert state["accelerate"]["observed"] == "single"
    assert (threadpool_info(), pa.cpu_count(), pa.io_thread_count()) == before
    assert os.environ["OPENBLAS_NUM_THREADS"] == "7"


@pytest.mark.parametrize("workers", [1, 2])
def test_direct_chromosome_logs_append_across_batches_and_reset_on_overwrite(tmp_path, workers, monkeypatch):
    import re
    from dataclasses import replace
    from ldsc import AnnotationBundle, LDScoreCalculator
    from ldsc.outputs import LDScoreOutputConfig
    from ldsc import ldscore_calculator as direct
    from tests.test_ldscore_query_batches import prepared_case
    pools = []
    original_pool = direct.ProcessPoolExecutor
    def pool(**kwargs):
        pools.append(kwargs["max_workers"])
        return original_pool(**kwargs)
    monkeypatch.setattr(direct, "ProcessPoolExecutor", pool)
    metadata, baseline, queries, config, panel, numerical = prepared_case(tmp_path)
    bundle = AnnotationBundle.from_frames(metadata, baseline, queries, config_snapshot=config)
    destination = tmp_path / "out"
    def run(overwrite=False):
        return LDScoreCalculator().run(bundle, panel, replace(numerical, workers=workers, query_batch_size=1),
                                       config, LDScoreOutputConfig(destination, overwrite=overwrite))
    run()
    assert pools == [workers]  # One persistent pool serves all three batches.
    pids = set()
    for chrom in ("1", "2"):
        path = destination / "diagnostics/chromosome-logs" / f"chr{chrom}.log"
        contents = path.read_text()
        assert contents.count(f"Starting chromosome {chrom}; batch") == 3
        assert "BLAS target=1" in contents
        assert "Arrow CPU=1; Arrow I/O=1" in contents
        pids.update(re.findall(r"Worker pid=(\d+)", contents))
        assert str(path) in (destination / "diagnostics/ldscore.log").read_text()
    assert len(pids) <= workers
    assert str(os.getpid()) not in pids
    with pytest.raises(FileExistsError):
        run()
    run(overwrite=True)
    assert pools == [workers, workers]
    assert (destination / "diagnostics/chromosome-logs/chr1.log").read_text().count("Starting chromosome 1; batch") == 3


def _fail_second_chromosome(*args, **kwargs):
    from ldsc.ldscore_calculator import LDScoreCalculator, _compute_one_chromosome
    original = LDScoreCalculator.compute_chromosome
    def fail(self, chrom, **options):
        if chrom == "22":
            raise OSError("simulated chromosome failure")
        return original(self, chrom=chrom, **options)
    LDScoreCalculator.compute_chromosome = fail
    try:
        return _compute_one_chromosome(*args, **kwargs)
    finally:
        LDScoreCalculator.compute_chromosome = original


def test_direct_failure_keeps_logs_but_removes_staged_reference_metadata(tmp_path, monkeypatch):
    from ldsc import ldscore_calculator as direct
    from tests.test_plink_workflow_resolution import write_inputs
    write_inputs(tmp_path)
    monkeypatch.setattr(direct, "_compute_one_chromosome", _fail_second_chromosome)
    args = direct.build_parser().parse_args([
        "--baseline-annot-sources", str(tmp_path / "baseline.annot"),
        "--plink-prefix", str(tmp_path / "1000G.EUR.QC.*"),
        "--snp-identifier", "chr_pos", "--genome-build", "hg19", "--ld-wind-cm", ".05",
        "--regr-snps-exclude-regions", "none", "--regr-snps-file", str(tmp_path / "regression.tsv"),
        "--output-dir", str(tmp_path / "out"), "--workers", "1", "--export-ref-metadata"])
    with pytest.raises(OSError, match="simulated chromosome failure"):
        direct.run_ldscore_from_args(args)
    output = tmp_path / "out"
    assert not (output / "metadata.json").exists()
    assert not (output / "ref_metadata").exists()
    assert not list(output.glob(".ldsc-*"))
    detail = (output / "diagnostics/chromosome-logs/chr22.log").read_text()
    assert "Traceback" in detail and "simulated chromosome failure" in detail
    main = (output / "diagnostics/ldscore.log").read_text()
    assert "chromosome 22" in main and "batch 1/1" in main
    assert "chr22.log" in main


def test_direct_failure_cancels_pending_without_submitting_more(monkeypatch):
    from concurrent.futures import Future
    from types import SimpleNamespace
    from ldsc import ldscore_calculator as direct
    submitted, shutdown = [], []
    failed, pending = Future(), Future()
    failed.set_exception(ValueError("task failure"))
    class Pool:
        def submit(self, function, chrom, *args, **kwargs):
            submitted.append(chrom)
            return failed if chrom == "1" else pending
        def shutdown(self, *, wait):
            shutdown.append(wait)
    monkeypatch.setattr(direct, "_slice_annotation_bundle", lambda bundle, chrom: None)
    with pytest.raises(ValueError, match="task failure"):
        direct.LDScoreCalculator()._run_chromosomes(
            ["1", "2", "3"], None, SimpleNamespace(spec=None), None, None, None, 2, pool=Pool())
    assert submitted == ["1", "2"]
    assert pending.cancelled()
    assert shutdown == [True]


def test_direct_waits_for_running_worker_before_query_batch_cleanup(tmp_path, monkeypatch):
    from concurrent.futures import ThreadPoolExecutor
    from contextlib import contextmanager
    from dataclasses import replace
    from threading import Event
    from ldsc import AnnotationBundle, LDScoreCalculator, _annotation_queries
    from ldsc import ldscore_calculator as direct
    from ldsc.outputs import LDScoreOutputConfig
    from tests.test_ldscore_query_batches import prepared_case

    metadata, baseline, queries, config, panel, numerical = prepared_case(tmp_path)
    bundle = AnnotationBundle.from_frames(metadata, baseline, queries, config_snapshot=config)
    started, release, finished = Event(), Event(), Event()
    events, scratch_paths = [], []
    error = ValueError("first chromosome failed while another was running")
    original_batch = _annotation_queries.execution_query_bundle

    @contextmanager
    def batch(*args, **kwargs):
        with original_batch(*args, **kwargs) as active:
            root = active.workspace.path
            scratch_paths.append(root)
            (root / "borrowed-input.txt").write_text("still needed")
            try:
                yield active
            finally:
                events.append("batch cleanup")
                assert finished.is_set(), "batch cleanup preceded worker shutdown"
                assert (root / "late-output.txt").read_text() == "still needed"

    class Pool(ThreadPoolExecutor):
        def __init__(self, *, max_workers, **kwargs):
            super().__init__(max_workers=max_workers)

        def shutdown(self, wait=True, **kwargs):
            events.append("shutdown")
            assert wait is True
            release.set()
            super().shutdown(wait=wait, **kwargs)
            events.append("workers finished")

    def task(chrom, active, *args, **kwargs):
        if chrom == "1":
            assert started.wait(timeout=10)
            raise error
        started.set()
        assert release.wait(timeout=10)
        root = active.workspace.path
        (root / "late-output.txt").write_text((root / "borrowed-input.txt").read_text())
        finished.set()
        return direct._ChromOutcome(chrom, None, True, "test worker")

    # Real futures/threads make cancellation fail for the running borrower.
    # Scientific calculations and child runtime settings are tested separately.
    monkeypatch.setattr(direct, "ProcessPoolExecutor", Pool)
    monkeypatch.setattr(direct, "_compute_one_chromosome", task)
    monkeypatch.setattr(_annotation_queries, "execution_query_bundle", batch)
    with pytest.raises(ValueError) as raised:
        LDScoreCalculator().run(bundle, panel, replace(numerical, workers=2, query_batch_size=1),
                                config, LDScoreOutputConfig(tmp_path / "out"))
    assert raised.value is error
    assert events.index("workers finished") < events.index("batch cleanup")
    assert len(scratch_paths) == 1
    assert not scratch_paths[0].exists()
    assert not (tmp_path / "out/metadata.json").exists()
    assert (tmp_path / "out/diagnostics/ldscore.log").is_file()


def _crash_direct_task(*args, **kwargs):
    from ldsc._ldscore_diagnostics import chromosome_log
    with chromosome_log(kwargs["log_path"], args[0], kwargs["batch_label"]):
        os._exit(17)


def _crash_indexed_task(*args, **kwargs):
    from ldsc._ldscore_diagnostics import chromosome_log
    with chromosome_log(args[-1], args[1], "indexed query batches"):
        os._exit(17)


@pytest.mark.parametrize("mode", ["direct", "indexed"])
def test_abnormal_ld_termination_reports_unconfirmed_source_and_retains_logs(tmp_path, monkeypatch, mode):
    from ldsc.errors import LDSCInternalError
    from concurrent.futures.process import BrokenProcessPool
    output = tmp_path / "out"
    if mode == "direct":
        from ldsc import ldscore_calculator as direct
        from tests.test_ldscore_parallelism import _write_index_panel, _PANEL_CHROMS, _run_minimal_ldscore
        panel = tmp_path / "panel"
        _write_index_panel(panel, _PANEL_CHROMS)
        monkeypatch.setattr(direct, "_compute_one_chromosome", _crash_direct_task)
        run = lambda: _run_minimal_ldscore(panel, output, 2)
    else:
        from ldsc import _indexed_ldscore_batches as indexed, gene_ldscore_index as core
        from tests.test_gene_index_streaming import two_chromosome_index
        index = two_chromosome_index(tmp_path)
        query = tmp_path / "query.txt"
        query.write_text("G1\nG3\n")
        monkeypatch.setattr(indexed, "_write_chromosome", _crash_indexed_task)
        run = lambda: core.run_indexed_ldscore(index, query_gene_list_sources=[query], output_dir=output, workers=2, _allow_partial_for_tests=True)
    before = {child.pid for child in mp.active_children()}
    with pytest.raises((LDSCInternalError, BrokenProcessPool)):
        run()
    assert {child.pid for child in mp.active_children()} == before
    main = (output / "diagnostics/ldscore.log").read_text()
    assert "Affected in-flight tasks:" in main
    assert "caused termination is unconfirmed" in main
    assert "chromosome-logs" in main
    assert len(list((output / "diagnostics/chromosome-logs").glob("chr*.log"))) == 2
    assert not (output / "metadata.json").exists()
    assert not list(output.glob(".ldsc-*"))


def _fail_indexed_second_batch(*args, **kwargs):
    from ldsc import _indexed_ldscore_batches as indexed
    original = indexed.pq.write_table
    def save(table, path, **options):
        if args[5] == 1 and path.name == "query.parquet":
            raise OSError("indexed batch failure")
        return original(table, path, **options)
    indexed.pq.write_table = save
    try:
        return indexed._write_chromosome(*args, **kwargs)
    finally:
        indexed.pq.write_table = original


def test_indexed_python_failure_reports_batch_and_full_chromosome_traceback(tmp_path, monkeypatch):
    from ldsc import _indexed_ldscore_batches as indexed, gene_ldscore_index as core
    from tests.test_gene_index_streaming import two_chromosome_index
    index = two_chromosome_index(tmp_path)
    queries = [tmp_path / f"query{i}.txt" for i in range(2)]
    for query in queries:
        query.write_text("G1\nG3\n")
    monkeypatch.setattr(indexed, "_write_chromosome", _fail_indexed_second_batch)
    output = tmp_path / "out"
    with pytest.raises(OSError, match="indexed batch failure"):
        core.run_indexed_ldscore(index, query_gene_list_sources=queries, output_dir=output, query_batch_size=1, _allow_partial_for_tests=True)
    main = (output / "diagnostics/ldscore.log").read_text()
    assert "chromosome 21" in main and "batch 2/2" in main and "chr21.log" in main
    detail = (output / "diagnostics/chromosome-logs/chr21.log").read_text()
    assert "Traceback" in detail and "indexed batch failure" in detail
    assert not (output / "metadata.json").exists()
    assert not list(output.glob(".ldsc-*"))


def _fail_index_builder_after_durable_first_chromosome(*args):
    from pathlib import Path
    from ldsc import gene_ldscore_index as builder
    if args[0] == "22":
        shard = args[-1] / "chromosomes/chr21/ldscore_operator.npz"
        assert shard.is_file()
        Path(os.environ["LDSC_TEST_DURABLE_SHARD"]).write_text(str(shard))
        raise RuntimeError("failure after chromosome 21 was staged")
    return builder._build_one_index_chromosome(*args)


def test_real_index_worker_failure_removes_completed_shards(tmp_path, monkeypatch):
    from pathlib import Path
    from ldsc import gene_ldscore_index as builder
    from tests.test_plink_workflow_resolution import write_inputs, index_args
    write_inputs(tmp_path)
    event = tmp_path / "durable-shard.txt"
    monkeypatch.setenv("LDSC_TEST_DURABLE_SHARD", str(event))
    monkeypatch.setattr(builder, "_build_one_index_chromosome", _fail_index_builder_after_durable_first_chromosome)
    args = index_args(tmp_path, str(tmp_path / "1000G.EUR.QC.*"), tmp_path / "index", workers=1)
    with pytest.raises(RuntimeError, match="chromosome 21 was staged"):
        builder.run_build_gene_ldscore_index_from_args(args)
    assert not Path(event.read_text()).exists()
    assert not (tmp_path / "index").exists()
    assert not list(tmp_path.glob(".index.stage-*"))
    log = (tmp_path / ".index.build-state/build-gene-ldscore-index.log").read_text()
    assert "Finished chromosome 21" in log and "chromosome 21 was staged" in log


def test_baseline_only_model_uses_one_child_for_fit_and_summary(tmp_path, monkeypatch):
    from ldsc import RegressionRunner, RegressionConfig
    from ldsc._runtime import serial_calculation as _serial_calculation
    from tests.test_regression_streaming import batch_inputs
    table, source = batch_inputs(tmp_path / "inputs")
    original, calls = _serial_calculation.ProcessPoolExecutor, []
    def pool(**kwargs):
        calls.append(kwargs["max_workers"])
        return original(**kwargs)
    monkeypatch.setattr(_serial_calculation, "ProcessPoolExecutor", pool)
    result = RegressionRunner(source.config_snapshot, RegressionConfig(n_blocks=6)).estimate_partitioned_h2_batch(
        table, source, query_columns=[], output_dir=tmp_path / "out")
    assert not result.summary.empty
    assert calls == [1]
