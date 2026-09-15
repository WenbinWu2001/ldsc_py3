"""Query errors must remain visible without discarding unrelated scan results."""

import json

import numpy as np
import pandas as pd
import pytest

from ldsc import cli, RegressionRunner, RegressionConfig
from ldsc.errors import LDSCInputError
from ldsc._kernel._jackknife import JackknifeIdentifiabilityError, LstsqJackknifeFast
from tests.test_regression_streaming import cli_inputs, batch_inputs


def concentrate_queries(tmp_path, names=("second",)):
    path = tmp_path / "ld/ldscore.query.parquet"
    queries = pd.read_parquet(path)
    for name in names:
        queries[name] = 0.
        queries.loc[22:24, name] = [1., 2., 3.]
    queries.to_parquet(path, index=False, row_group_size=30)


@pytest.mark.parametrize("threads", [1, 2])
def test_continue_publishes_valid_queries_and_logs_singular_block(tmp_path, threads):
    args, _ = cli_inputs(tmp_path)
    args.extend(["--threads", str(threads)])
    before = cli.main([*args, "--summary-sort-by", "category"])
    concentrate_queries(tmp_path)

    result = cli.main([*args, "--continue-on-query-error", "--summary-sort-by", "category", "--overwrite"])

    output = tmp_path / "out"
    assert result.category.tolist() == ["first", "third"]
    pd.testing.assert_frame_equal(result, before.loc[before.category != "second"].reset_index(drop=True))
    assert np.isfinite(result.coefficient).all()
    status = pd.read_csv(output / "diagnostics/query_status.tsv", sep="\t", keep_default_na=False)
    assert status.query_annotation.tolist() == ["first", "second", "third"]
    assert status.status.tolist() == ["success", "unestimable", "success"]
    assert status.loc[1, "stage"] == "estimator"
    assert status.loc[1, "error_type"] == "JackknifeIdentifiabilityError"
    metadata = json.loads((output / "diagnostics/metadata.json").read_text())
    assert metadata["query_error_policy"] == "continue"
    assert metadata["n_queries_failed"] == 1
    assert metadata["n_queries_successful"] == 2
    assert metadata["files"]["query_status"] == "diagnostics/query_status.tsv"
    manifest = pd.read_csv(output / "diagnostics/query_annotations/manifest.tsv", sep="\t")
    assert manifest.query_annotation.tolist() == ["first", "third"]
    log = (output / "diagnostics/partitioned-h2.log").read_text()
    assert "Query 'second' failed at estimator" in log
    assert "block 3/6" in log
    assert "rows [19, 29)" in log
    assert '"CHR": "1", "POS_start": 21, "POS_end": 30' in log
    assert '"second": {"full": 3, "remaining": 0}' in log
    assert "Traceback" in log
    assert "2 successful, 1 failed" in log
    assert not (output / "RUN_FAILED.txt").exists()


@pytest.mark.parametrize("threads", [1, 2])
def test_strict_default_collects_query_failures_without_publishing(tmp_path, threads):
    args, _ = cli_inputs(tmp_path)
    args.extend(["--threads", str(threads)])
    concentrate_queries(tmp_path, ("second", "third"))
    with pytest.raises(JackknifeIdentifiabilityError) as caught:
        cli.main(args)
    assert "2/3 queries failed" in " ".join(caught.value.__notes__)
    output = tmp_path / "out"
    statuses = pd.read_csv(output / "diagnostics/query_status.tsv", sep="\t")
    assert statuses.status.tolist() == ["success", "unestimable", "unestimable"]
    assert not (output / "partitioned_h2.tsv").exists()
    assert not (output / "diagnostics/query_annotations").exists()
    assert not list(output.glob(".ldsc-annotation-*"))


@pytest.mark.parametrize("threads", [1, 2])
def test_continue_with_no_successes_fails_and_keeps_all_statuses(tmp_path, threads):
    args, _ = cli_inputs(tmp_path)
    args.extend(["--threads", str(threads)])
    concentrate_queries(tmp_path, ("first", "second", "third"))
    with pytest.raises(JackknifeIdentifiabilityError) as caught:
        cli.main([*args, "--continue-on-query-error"])
    assert "No query fit succeeded" in " ".join(caught.value.__notes__)
    output = tmp_path / "out"
    statuses = pd.read_csv(output / "diagnostics/query_status.tsv", sep="\t")
    assert statuses.status.tolist() == ["unestimable"] * 3
    assert not (output / "partitioned_h2.tsv").exists()


def test_continue_also_catches_unexpected_numerical_exceptions(tmp_path, monkeypatch):
    args, _ = cli_inputs(tmp_path)
    concentrate_queries(tmp_path)
    original = np.linalg.solve

    def solve(matrix, rhs):
        if matrix.shape == (3, 3) and not matrix[1].any():
            raise RuntimeError("injected numerical backend failure")
        return original(matrix, rhs)

    monkeypatch.setattr(np.linalg, "solve", solve)
    result = cli.main([*args, "--continue-on-query-error"])
    assert set(result.category) == {"first", "third"}
    status = pd.read_csv(tmp_path / "out/diagnostics/query_status.tsv", sep="\t")
    assert status.loc[1, "status"] == "failed"
    assert status.loc[1, "error_type"] == "RuntimeError"
    assert status.loc[1, "error_message"] == "injected numerical backend failure"


def test_nonzero_columns_can_become_dependent_after_deletion():
    t = np.arange(20, dtype=float)
    query = t.copy()
    query[8:12] += [1., 3., 2., 4.]
    design = np.column_stack([np.ones(20), t, query])
    assert np.linalg.matrix_rank(design) == 3
    response = (design @ [1., 2., 3.]).reshape(-1, 1)
    with pytest.raises(JackknifeIdentifiabilityError) as caught:
        LstsqJackknifeFast(design, response, n_blocks=5)
    assert caught.value.failures == [{
        "block_index": 2, "rank": 2, "n_parameters": 3,
        "zero_column_indices": [], "row_start": 8, "row_end": 12,
    }]


def test_all_singular_deletions_are_reported_together():
    design = np.column_stack([np.ones(12), np.r_[np.ones(4), np.zeros(8)], np.r_[np.zeros(4), np.ones(4), np.zeros(4)]])
    assert np.linalg.matrix_rank(design) == 3
    with pytest.raises(JackknifeIdentifiabilityError) as caught:
        LstsqJackknifeFast(design, (design @ [1., 2., 3.]).reshape(-1, 1), n_blocks=3)
    assert [failure["block_index"] for failure in caught.value.failures] == [0, 1, 2]


@pytest.mark.parametrize("error_type", [KeyboardInterrupt, SystemExit])
def test_continue_does_not_swallow_process_control(tmp_path, monkeypatch, error_type):
    args, _ = cli_inputs(tmp_path)

    def stop(*args, **kwargs):
        raise error_type("requested stop")

    monkeypatch.setattr(np.linalg, "solve", stop)
    with pytest.raises(error_type):
        cli.main([*args, "--continue-on-query-error"])
    assert not (tmp_path / "out/partitioned_h2.tsv").exists()


@pytest.mark.parametrize("threads", [1, 2])
def test_continue_does_not_swallow_publication_failure(tmp_path, monkeypatch, threads):
    import os

    args, _ = cli_inputs(tmp_path)
    args.extend(["--threads", str(threads)])
    original = os.replace

    def replace(source, destination):
        if ".query_annotations.tmp." in str(source) and str(destination).endswith("query_annotations"):
            raise OSError("output device unavailable")
        return original(source, destination)

    monkeypatch.setattr(os, "replace", replace)
    with pytest.raises(OSError, match="output device unavailable"):
        cli.main([*args, "--continue-on-query-error"])
    assert not (tmp_path / "out/diagnostics/query_annotations").exists()


@pytest.mark.parametrize("threads", [1, 2])
def test_python_api_records_query_preparation_error_and_baseline_retry_removes_status(tmp_path, threads):
    table, source = batch_inputs(tmp_path / "ld")
    source.count_records[:] = [record for record in source.count_records if record["column"] != "second"]
    runner = RegressionRunner(source.config_snapshot, RegressionConfig(n_blocks=6))
    output = tmp_path / "out"
    result = runner.estimate_partitioned_h2_batch(
        table, source, output_dir=output, continue_on_query_error=True, threads=threads,
    )
    assert set(result.summary.category) == {"first", "third"}
    assert result.query_status.status.tolist() == ["success", "failed", "success"]
    assert result.query_status.loc[1, "stage"] == "model_preparation"
    assert result.query_status.loc[1, "error_type"] == "LDSCInputError"
    assert set(result.per_query_artifacts) == {"first", "third"}
    runner.estimate_partitioned_h2_batch(table, source, output_dir=output, query_columns=[], overwrite=True, threads=threads)
    assert not (output / "diagnostics/query_status.tsv").exists()
    assert not (output / "diagnostics/query_annotations").exists()


def test_continue_does_not_swallow_shared_input_loading_error(tmp_path):
    args, _ = cli_inputs(tmp_path)
    (tmp_path / "ld/ldscore.query.parquet").unlink()
    with pytest.raises(LDSCInputError):
        cli.main([*args, "--continue-on-query-error"])
    assert not (tmp_path / "out/partitioned_h2.tsv").exists()


def test_python_api_requires_explicit_boolean_continuation(tmp_path):
    runner = RegressionRunner()
    with pytest.raises(ValueError, match="continue_on_query_error must be a boolean"):
        runner.estimate_partitioned_h2_batch(None, None, output_dir=tmp_path, continue_on_query_error="False")


def test_parallel_strict_overwrite_preserves_previous_science_and_marks_failure(tmp_path):
    args, _ = cli_inputs(tmp_path)
    cli.main([*args, "--threads", "2"])
    output = tmp_path / "out"
    previous = (output / "partitioned_h2.tsv").read_bytes()
    manifest = (output / "diagnostics/query_annotations/manifest.tsv").read_bytes()
    concentrate_queries(tmp_path)
    with pytest.raises(JackknifeIdentifiabilityError):
        cli.main([*args, "--threads", "2", "--overwrite"])
    assert (output / "partitioned_h2.tsv").read_bytes() == previous
    assert (output / "diagnostics/query_annotations/manifest.tsv").read_bytes() == manifest
    assert (output / "RUN_FAILED.txt").exists()
    assert not list(output.glob(".ldsc-annotation-*"))
