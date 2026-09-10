"""Selective input and fit ownership checks for large pathway batches."""

import weakref
import json

import numpy as np
import pandas as pd
import pytest

from ldsc import regression_runner as workflow
from ldsc import cli
from ldsc.outputs import LDScoreDirectoryWriter, LDScoreOutputConfig
from tests import test_regression_workflow as fixtures
from tests.test_regression_fit_outcomes import h2_dataset
from ldsc.ldscore_calculator import LDScoreResult
from ldsc.overlap_matrix import LDScoreOverlap
from ldsc._kernel.overlap import OverlapContribution


def ld_directory(tmp_path):
    result = fixtures.RegressionWorkflowTest().make_ldscore_result()
    LDScoreDirectoryWriter().write(result, LDScoreOutputConfig(output_dir=tmp_path))
    return result


def test_directory_loading_defers_query_values_and_reads_only_requested_columns(tmp_path, monkeypatch):
    expected = ld_directory(tmp_path)
    original = pd.read_parquet
    query_reads = []

    def read(path, *args, **kwargs):
        if str(path).endswith("ldscore.query.parquet"):
            query_reads.append(kwargs.get("columns"))
        return original(path, *args, **kwargs)

    monkeypatch.setattr(pd, "read_parquet", read)
    source = workflow.load_ldscore_from_dir(str(tmp_path))
    assert query_reads and all(cols is not None for cols in query_reads)
    assert all(not set(cols).intersection(expected.query_columns) for cols in query_reads)
    assert not hasattr(source, "query_table")
    first = source.read_queries(["query2"])
    pd.testing.assert_series_equal(first.query2, expected.query_table.query2.astype(np.float32))
    assert query_reads[-1] == ["query2"]
    ref = weakref.ref(first)
    del first
    assert ref() is None


def test_unused_query_values_are_not_validated_but_missing_schema_still_fails(tmp_path):
    ld_directory(tmp_path)
    path = tmp_path / "ldscore.query.parquet"
    query = pd.read_parquet(path)
    query["query2"] = "unused invalid numeric values"
    query.to_parquet(path, index=False)
    source = workflow.load_ldscore_from_dir(str(tmp_path))
    runner = workflow.RegressionRunner()
    table = fixtures.RegressionWorkflowTest().make_sumstats_table()
    dataset = runner.build_dataset(table, source, query_columns=["query1"])
    assert dataset.retained_ld_columns == ["base", "query1"]
    query.drop(columns="query2").to_parquet(path, index=False)
    with pytest.raises(Exception, match="query2"):
        workflow.load_ldscore_from_dir(str(tmp_path))


def test_shared_preparation_preserves_order_and_model_specific_column_filtering():
    fixture = fixtures.RegressionWorkflowTest()
    result = fixture.make_ldscore_result()
    table = fixture.make_sumstats_table()
    table = workflow.replace(table, data=table.data.iloc[[2, 0, 1]].reset_index(drop=True))
    result.query_table["query2"] = 1.0
    runner = workflow.RegressionRunner()
    prepared = runner.prepare_inputs(table, result)
    first = runner.dataset_from_prepared(prepared, query_columns=["query1"])
    second = runner.dataset_from_prepared(prepared, query_columns=["query2"])
    assert first.merged.SNP.tolist() == ["rs1", "rs2", "rs3"]
    np.testing.assert_array_equal(first.merged.Z, [2., 1., .5])
    np.testing.assert_array_equal(first.merged.query1, [1., 2., 3.])
    assert first.retained_ld_columns == ["base", "query1"]
    assert second.retained_ld_columns == ["base"]
    assert second.dropped_zero_variance_ld_columns == ["query2"]
    np.testing.assert_array_equal(first.reference_snp_count_totals[workflow.COMMON_COUNT_KEY], [8, 18])
    np.testing.assert_array_equal(second.reference_snp_count_totals[workflow.COMMON_COUNT_KEY], [8])
    assert list(prepared.dataset.merged) == list(runner.build_dataset(table, result).merged)


def batch_inputs(tmp_path):
    frame = h2_dataset().merged
    baseline = frame[["SNP", "base"]].assign(
        CHR=["1"] * 30 + ["2"] * 30, POS=np.tile(np.arange(1, 31), 2), regression_ld_scores=1.,
    )
    query = baseline[["CHR", "SNP", "POS"]].assign(
        first=2 + np.sin(frame.base), second=2 + np.cos(frame.base), third=2 + np.sin(.3 * frame.base),
    )
    overlap = LDScoreOverlap.from_contribution(
        OverlapContribution(np.array([[1000., 250., 300., 400.]]), np.array([[1000., 250., 300., 400.]]),
                            np.array([250., 300., 400.]), np.array([250., 300., 400.]), 1000, 1000),
        baseline_columns=["base"], query_columns=["first", "second", "third"],
    )
    source = LDScoreResult(
        baseline, query, [{"group": group, "column": column, "all_reference_snp_count": count,
                          "common_reference_snp_count": count}
                         for group, column, count in [("baseline", "base", 1000.), ("query", "first", 250.),
                                                      ("query", "second", 300.), ("query", "third", 400.)]],
        ["base"], ["first", "second", "third"], frozenset(), frozenset(baseline.SNP), [],
        config_snapshot=workflow.GlobalConfig(snp_identifier="rsid"), overlap=overlap,
    )
    LDScoreDirectoryWriter().write(source, LDScoreOutputConfig(output_dir=tmp_path))
    table = workflow.SumstatsTable(data=frame[["SNP", "Z", "N"]], trait_name="trait", has_alleles=False,
                                   source_path="trait", config_snapshot=source.config_snapshot)
    return table, workflow.load_ldscore_from_dir(str(tmp_path))


@pytest.mark.parametrize("batch_size", [1, 2, 1000])
def test_batch_fits_match_separate_models_and_release_details(tmp_path, monkeypatch, batch_size):
    table, source = batch_inputs(tmp_path / "ld")
    runner = workflow.RegressionRunner(source.config_snapshot, workflow.RegressionConfig(n_blocks=6))
    expected = {}
    for query in source.query_columns:
        dataset = runner.build_dataset(table, source, query_columns=[query])
        outcome = runner._fit_h2_dataset(dataset)
        expected[query] = (
            workflow.summarize_partitioned_h2(outcome.estimator, outcome.dataset, [query]),
            workflow.summarize_partitioned_h2(outcome.estimator, outcome.dataset, outcome.dataset.retained_ld_columns),
            workflow._coefficient_delete_frame(outcome.estimator, outcome.dataset.retained_ld_columns),
        )
    del outcome, dataset
    original_prepare, original_fit = runner.prepare_inputs, runner._fit_h2_dataset
    preparations, live, reads = [], [], []
    original_read = type(source).read_queries

    def read(source, columns):
        reads.append(list(columns))
        return original_read(source, columns)

    def prepare(*args, **kwargs):
        preparations.append(1)
        return original_prepare(*args, **kwargs)

    def fit(*args, **kwargs):
        assert all(ref() is None for ref in live), "previous estimator details remain live"
        outcome = original_fit(*args, **kwargs)
        live.append(weakref.ref(outcome.estimator))
        return outcome

    monkeypatch.setattr(runner, "prepare_inputs", prepare)
    monkeypatch.setattr(runner, "_fit_h2_dataset", fit)
    monkeypatch.setattr(type(source), "read_queries", read)
    result = runner.estimate_partitioned_h2_batch(table, source, output_dir=tmp_path / "out",
                                                query_batch_size=batch_size, summary_sort_by="category")
    assert preparations == [1]
    assert reads == [source.query_columns[i:i + batch_size] for i in range(0, 3, batch_size)]
    assert all(ref() is None for ref in live)
    pd.testing.assert_frame_equal(result.summary, pd.concat([expected[q][0] for q in sorted(expected)], ignore_index=True))
    for query, paths in result.per_query_artifacts.items():
        pd.testing.assert_frame_equal(pd.read_csv(paths.full, sep="\t"), expected[query][1], check_dtype=False)
        pd.testing.assert_frame_equal(pd.read_parquet(paths.delete_values), expected[query][2])
        assert json.loads(paths.metadata.read_text())["n_snps"] == 58
    assert not list((tmp_path / "out").glob(".ldsc-annotation-*"))


def cli_inputs(tmp_path):
    table, source = batch_inputs(tmp_path / "ld")
    path = tmp_path / "trait.parquet"
    fixtures.RegressionWorkflowTest().write_footer_sumstats_parquet(path, frame=table.data, trait_name="MDD")
    return ["partitioned-h2", "--ldscore-dir", str(tmp_path / "ld"), "--sumstats-file", str(path),
            "--output-dir", str(tmp_path / "out"), "--n-blocks", "6", "--query-batch-size", "2"], source


@pytest.mark.parametrize("sort_by", ["auto", *workflow.PARTITIONED_H2_SUMMARY_SORT_COLUMNS])
def test_real_partitioned_command_preserves_sorted_artifacts_metadata_and_overwrite(tmp_path, sort_by):
    args, source = cli_inputs(tmp_path)
    output = tmp_path / "out"
    result = cli.main([*args, "--summary-sort-by", sort_by])
    written = pd.read_csv(output / "partitioned_h2.tsv", sep="\t")
    pd.testing.assert_frame_equal(written, result, check_dtype=False)
    key = "coefficient_p" if sort_by == "auto" else workflow.PARTITIONED_H2_SUMMARY_SORT_COLUMNS[sort_by]
    ascending = sort_by in {"auto", "category", "enrichment-p", "coefficient-p"}
    assert written[key].tolist() == sorted(written[key], reverse=not ascending)
    manifest = pd.read_csv(output / "diagnostics/query_annotations/manifest.tsv", sep="\t")
    assert manifest.query_annotation.tolist() == written.category.tolist()
    assert manifest.folder.tolist() == [f"{i:04d}_{name}" for i, name in enumerate(written.category, 1)]
    metadata = json.loads((output / "diagnostics/metadata.json").read_text())
    assert metadata["trait_name"] == "MDD"
    assert metadata["analysis_type"] == "cell_type_specific"
    assert metadata["headline_metric"] == "coefficient"
    assert metadata["coefficient_p_test"] == "one_sided_greater"
    assert metadata["enrichment_p_test"] == "two_sided_t"
    assert metadata["count_kind"] == "common"
    for row in manifest.itertuples():
        assert (output / row.metadata_path).is_file()
        assert pd.read_csv(output / row.partitioned_h2_full_path, sep="\t").category.tolist() == ["base", row.query_annotation]
        deletes = pd.read_parquet(output / row.coefficient_delete_values_path)
        assert deletes.columns.tolist() == ["delete_block", "base", row.query_annotation]
        assert deletes[row.query_annotation].dtype == np.float64
    before = (output / "partitioned_h2.tsv").read_bytes()
    with pytest.raises(FileExistsError):
        cli.main(args)
    assert (output / "partitioned_h2.tsv").read_bytes() == before
    stale = output / "diagnostics/query_annotations/stale"
    stale.mkdir()
    cli.main([*args, "--overwrite", "--summary-sort-by", sort_by])
    assert not stale.exists()
    assert (output / "diagnostics/partitioned-h2.log").is_file()
    assert not (output / "RUN_FAILED.txt").exists()


def test_functional_command_keeps_complete_model_at_root(tmp_path):
    args, source = cli_inputs(tmp_path)
    baseline = pd.concat([source.baseline_table, source.read_queries(["first"])[["first"]]], axis=1)
    overlap = LDScoreOverlap.from_contribution(
        OverlapContribution(np.array([[1000., 250.], [250., 250.]]), np.array([[1000., 250.], [250., 250.]]),
                            np.array([]), np.array([]), 1000, 1000), baseline_columns=["base", "first"], query_columns=[],
    )
    result = LDScoreResult(baseline, None, [{**r, "group": "baseline"} for r in source.count_records if r["column"] in ("base", "first")],
                          ["base", "first"], [], frozenset(), source.ld_regression_snps, [], config_snapshot=source.config_snapshot,
                          overlap=overlap)
    LDScoreDirectoryWriter().write(result, LDScoreOutputConfig(output_dir=tmp_path / "ld", overwrite=True))
    summary = cli.main(args)
    assert summary.category.tolist() == ["base", "first"]
    metadata = json.loads((tmp_path / "out/diagnostics/metadata.json").read_text())
    assert metadata["analysis_type"] == "functional_category"
    assert metadata["headline_metric"] == "enrichment"
    assert metadata["n_snps"] == 58
    assert metadata["n_blocks_used"] == 6
    assert metadata["retained_ld_columns"] == ["base", "first"]
    assert pd.read_parquet(tmp_path / "out/diagnostics/coefficient_delete_values.parquet").columns.tolist() == ["delete_block", "base", "first"]
    assert not (tmp_path / "out/diagnostics/query_annotations").exists()


@pytest.mark.parametrize("failure", ["fit", "publication"])
def test_late_failures_clean_private_details_and_preserve_existing_failure_contract(tmp_path, monkeypatch, failure):
    from ldsc import outputs

    args, _ = cli_inputs(tmp_path)
    cli.main(args)
    output = tmp_path / "out"
    prior_manifest = (output / "diagnostics/query_annotations/manifest.tsv").read_bytes()
    staged = []
    original_stage = workflow.stage_partitioned_h2_fit

    def stage(*args, **kwargs):
        artifact = original_stage(*args, **kwargs)
        staged.append(artifact)
        return artifact

    monkeypatch.setattr(workflow, "stage_partitioned_h2_fit", stage)
    if failure == "fit":
        original_fit = workflow.RegressionRunner._fit_h2_dataset

        def fit(*args, **kwargs):
            if staged:
                assert staged[0].delete_values.is_file()
                raise RuntimeError("late fitting failure")
            return original_fit(*args, **kwargs)

        monkeypatch.setattr(workflow.RegressionRunner, "_fit_h2_dataset", fit)
    else:
        original_replace = outputs.os.replace

        def replace(source, destination):
            if ".query_annotations.tmp." in str(source) and str(destination).endswith("query_annotations"):
                raise RuntimeError("publication failure")
            return original_replace(source, destination)

        monkeypatch.setattr(outputs.os, "replace", replace)
    with pytest.raises(RuntimeError, match="failure"):
        cli.main([*args, "--overwrite"])
    assert staged
    assert all(not artifact.full.exists() for artifact in staged)
    assert (output / "RUN_FAILED.txt").is_file()
    assert (output / "diagnostics/query_annotations/manifest.tsv").read_bytes() == prior_manifest
    assert not list(output.glob(".ldsc-annotation-*"))


@pytest.mark.parametrize("size", [0, -1])
def test_invalid_query_batch_size_stops_before_input_loading(tmp_path, monkeypatch, size):
    monkeypatch.setattr(workflow, "_load_sumstats_table", lambda *_args: pytest.fail("invalid batch size reached input loading"))
    with pytest.raises(ValueError, match="query_batch_size"):
        cli.main(["partitioned-h2", "--sumstats-file", "missing", "--ldscore-dir", "missing",
                  "--output-dir", str(tmp_path / "out"), "--query-batch-size", str(size)])
