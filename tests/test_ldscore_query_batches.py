"""Scientific and resource contracts for sequential query execution batches."""

import math
import json
import weakref
from dataclasses import replace

import numpy as np
import pandas as pd
import pytest

from ldsc._kernel.ldscore_projection import ArrayAnnotations
from ldsc._kernel.overlap import annotation_statistics
from ldsc.outputs import LDScoreDirectoryWriter, LDScoreOutputConfig
from ldsc.regression_runner import load_ldscore_from_dir
from tests.test_output import make_split_ldscore_result


def prepared_case(tmp_path):
    from ldsc import GlobalConfig, RefPanelConfig, LDScoreConfig
    from ldsc._kernel.ref_panel import ParquetR2RefPanel
    from tests.test_ldscore_parallelism import _PANEL_CHROMS, _write_index_panel
    panel_dir = tmp_path / "panel"
    _write_index_panel(panel_dir, _PANEL_CHROMS)
    metadata = pd.concat([pair[0] for pair in _PANEL_CHROMS.values()], ignore_index=True)
    baseline = pd.DataFrame({"base": np.ones(len(metadata))})
    queries = pd.DataFrame({"first": [1, 0, 1, 0, 1, 0, 0],
                            "second": [0, 1, 0, 1, 0, 1, 0],
                            "third": [.3, -1, .1, 2, .5, .7, -.2]})
    config = GlobalConfig(snp_identifier="rsid")
    panel = ParquetR2RefPanel(config, RefPanelConfig(backend="parquet_r2", r2_dir=str(panel_dir)))
    return metadata, baseline, queries, config, panel, LDScoreConfig(ld_wind_kb=1, whole_chromosome_ok=True)


def test_prepared_single_batch_python_calculation_writes_nothing(tmp_path, monkeypatch):
    import builtins
    import os
    from pathlib import Path
    from ldsc import AnnotationBundle, LDScoreCalculator
    metadata, baseline, queries, config, panel, numerical = prepared_case(tmp_path)
    real_open, real_os_open = builtins.open, os.open
    def open_read_only(file, mode="r", *args, **kwargs):
        assert not any(flag in mode for flag in "wax+"), f"unexpected write: {file}"
        return real_open(file, mode, *args, **kwargs)
    def os_open_read_only(file, flags, *args, **kwargs):
        assert not flags & (os.O_WRONLY | os.O_RDWR | os.O_CREAT | os.O_TRUNC), f"unexpected write: {file}"
        return real_os_open(file, flags, *args, **kwargs)
    def no_mkdir(*args, **kwargs):
        pytest.fail("in-memory calculation created a directory")
    monkeypatch.setattr(builtins, "open", open_read_only)
    monkeypatch.setattr(os, "open", os_open_read_only)
    monkeypatch.setattr(Path, "mkdir", no_mkdir)
    bundle = AnnotationBundle.from_frames(metadata, baseline, queries, config_snapshot=config)
    result = LDScoreCalculator().run(bundle, panel, numerical, config)
    bundle.close()
    assert result.query_columns == list(queries)
    assert result.baseline_table.SNP.tolist() == metadata.SNP.tolist()
    np.testing.assert_allclose([row["all_reference_snp_count"] for row in result.count_records],
                               [7, 3, 3, math.fsum(queries.third.astype(np.float32))], rtol=1e-6, atol=1e-8)
    assert all(hasattr(diagnostics, "frame") for diagnostics in result.identity_drops_by_chrom.values())
    assert all(list(diagnostics.frames()) == [] for diagnostics in result.identity_drops_by_chrom.values())
    assert result.output_paths == {}


@pytest.mark.parametrize("width,workers", [(1, 1), (2, 1), (10, 1), (2, 2)])
def test_direct_batches_match_single_calculation(tmp_path, width, workers):
    from ldsc import AnnotationBundle, LDScoreCalculator
    metadata, baseline, queries, config, panel, numerical = prepared_case(tmp_path)
    bundle = AnnotationBundle.from_frames(metadata, baseline, queries, config_snapshot=config)
    calculator = LDScoreCalculator()
    expected = calculator.run(bundle, panel, numerical, config)
    saved = calculator.run(bundle, panel, replace(numerical, query_batch_size=width, threads=workers), config,
                           output_config=LDScoreOutputConfig(tmp_path / "out"))
    assert not hasattr(saved, "query_table")
    np.testing.assert_allclose(saved.read_queries(list(queries))[queries.columns],
                               expected.query_table[queries.columns], rtol=1e-6, atol=1e-8)
    pd.testing.assert_frame_equal(saved.baseline_table, expected.baseline_table, check_dtype=False)
    assert saved.count_records == expected.count_records
    pd.testing.assert_frame_equal(saved.overlap.baseline_block_all, expected.overlap.baseline_block_all)
    pd.testing.assert_series_equal(saved.overlap.query_diagonal_all, expected.overlap.query_diagonal_all)


def test_multiple_batches_require_output_before_computation(tmp_path):
    from ldsc import AnnotationBundle, LDScoreCalculator
    from ldsc.errors import LDSCUsageError
    metadata, baseline, queries, config, panel, numerical = prepared_case(tmp_path)
    bundle = AnnotationBundle.from_frames(metadata, baseline, queries, config_snapshot=config)
    with pytest.raises(LDSCUsageError, match="Multiple query execution batches require"):
        LDScoreCalculator().run(bundle, panel, replace(numerical, query_batch_size=2), config)


def test_prepared_inputs_preserve_global_duplicates_and_allele_contract():
    from ldsc import AnnotationBundle, GlobalConfig
    from ldsc.errors import LDSCInputError
    metadata = pd.DataFrame({"CHR": [1, 2, 2], "SNP": ["duplicate", "duplicate", "unique"], "POS": [10, 20, 30]})
    baseline = pd.DataFrame({"base": [1., 2., 3.]})
    config = GlobalConfig(snp_identifier="rsid")
    with AnnotationBundle.from_frames(metadata, baseline, config_snapshot=config) as bundle:
        assert bundle.chromosomes == ["2"]
        np.testing.assert_array_equal(bundle.read("2"), [[3.]])
        assert bundle.identity_drops.frame.SNP.tolist() == ["duplicate", "duplicate"]
    with pytest.raises(LDSCInputError, match="both allele columns or neither"):
        AnnotationBundle.from_frames(metadata.assign(A1="A"), baseline, config_snapshot=config)


def test_direct_bed_preparation_releases_each_execution_batch(tmp_path, monkeypatch):
    from ldsc import run_ldscore, get_global_config, set_global_config
    from ldsc._annotation_storage import ColumnStore
    metadata, baseline, _, config, _, _ = prepared_case(tmp_path)
    baseline_path = tmp_path / "baseline.annot"
    pd.concat([metadata, baseline], axis=1).to_csv(baseline_path, sep="\t", index=False)
    query_paths = []
    for ordinal, pos in enumerate([100, 200, 300]):
        path = tmp_path / f"q{ordinal}.bed"
        path.write_text(f"1\t{pos-1}\t{pos}\n2\t{pos-1}\t{pos}\n")
        query_paths.append(str(path))
    restriction = tmp_path / "snps.tsv"
    metadata[["SNP"]].to_csv(restriction, sep="\t", index=False)
    original = ColumnStore.create
    current_parent = None
    previous_paths = []
    def create(path, n_rows, columns, **kwargs):
        nonlocal current_parent
        focal = [column for column in columns if column.startswith("q")]
        if focal:
            assert len(focal) == 1, "all-query preparation exceeds execution width"
            if current_parent is not None and path.parent != current_parent:
                assert all(not old.exists() for old in previous_paths), "completed query scratch remains live"
            current_parent = path.parent
            previous_paths.append(path)
        return original(path, n_rows, columns, **kwargs)
    monkeypatch.setattr(ColumnStore, "create", create)
    previous_config = get_global_config()
    try:
        set_global_config(config)
        saved = run_ldscore(baseline_annot_sources=str(baseline_path), query_annot_bed_sources=query_paths,
            r2_dir=str(tmp_path / "panel"), output_dir=str(tmp_path / "out"),
            regr_snps_file=str(restriction), regr_snps_exclude_regions="none",
            ld_wind_kb=1, yes_really=True, query_batch_size=1)
    finally:
        set_global_config(previous_config)
    assert saved.query_columns == ["q0", "q1", "q2"]
    assert len(saved.query_batches) == 3
    assert all(not path.exists() for path in previous_paths)


@pytest.mark.parametrize("width", [1, 1000])
@pytest.mark.parametrize("layout", ["C", "F"])
def test_quantitative_counts_match_normalized_values_across_batch_layouts(width, layout):
    rng = np.random.default_rng(182)
    signed = -rng.normal(size=6000).astype(np.float32)
    high = (1e6 + np.arange(6000) % 17).astype(np.float32)
    values = np.ones((6000, 1001), dtype=np.float32, order=layout)
    values[:, 1], values[:, 2] = signed, high
    names = ("base", "signed", "high", *[f"query{i}" for i in range(998)])
    common = np.arange(6000) % 3 != 0
    metadata = pd.DataFrame({"MAF": np.where(common, .2, .01)})
    counts, common_counts, overlap, _ = annotation_statistics(
        metadata, ArrayAnnotations(values, names), 1, query_batch_size=width,
    )
    expected = [math.fsum(values[:, j]) for j in range(3)]
    expected_common = [math.fsum(values[common, j]) for j in range(3)]
    np.testing.assert_allclose(counts[:3], expected, rtol=1e-6, atol=1e-8)
    np.testing.assert_allclose(common_counts[:3], expected_common, rtol=1e-6, atol=1e-8)
    np.testing.assert_allclose(overlap.baseline_block_all[0, :3], expected, rtol=1e-6, atol=1e-8)


def test_batch_writer_releases_tables_and_reads_across_files(tmp_path):
    refs = []

    def batches():
        for names, values in [(["first", "second"], [[.5, 1.5], [.7, 1.7]]),
                              (["third"], [[2.5], [2.7]])]:
            assert all(ref() is None for ref in refs)
            result = make_split_ldscore_result()
            table = result.query_table.drop(columns="query")
            table[names] = np.asarray(values)
            result = replace(result, query_table=table, query_columns=names,
                             count_records=[result.count_records[0], *[
                                 {**result.count_records[1], "column": name} for name in names]])
            refs.append(weakref.ref(table))
            yield result
            del result, table

    source = LDScoreDirectoryWriter().write_batches(batches(), LDScoreOutputConfig(tmp_path))
    assert all(ref() is None for ref in refs)
    assert not hasattr(source, "query_table")
    assert source.query_columns == ["first", "second", "third"]
    metadata = json.loads((tmp_path / "metadata.json").read_text())
    assert [entry["file"] for entry in metadata["query_batches"]] == [
        "ldscore.query.batch00001.parquet", "ldscore.query.batch00002.parquet"]
    loaded = load_ldscore_from_dir(str(tmp_path))
    table = loaded.read_queries(["third", "first", "second"])
    assert table.columns.tolist() == ["CHR", "SNP", "POS", "third", "first", "second"]
    np.testing.assert_allclose(table[["third", "first", "second"]], [[2.5, .5, 1.5], [2.7, .7, 1.7]])
    assert loaded.count_records[0]["all_reference_snp_count"] == 10
    assert not list(tmp_path.glob(".ldsc-batches-*"))
def test_batch_write_failure_closes_producer_and_removes_scratch(tmp_path, monkeypatch):
    from ldsc.outputs import LDScoreDirectoryWriter, LDScoreOutputConfig
    from ldsc import _ldscore_batch_output

    closed = []
    def results():
        try:
            yield make_split_ldscore_result()
        finally:
            closed.append(True)
    producer = results()
    def fail(*args, **kwargs):
        raise OSError("simulated disk failure")
    monkeypatch.setattr(_ldscore_batch_output, "_write_chromosome_aligned_parquet", fail)
    with pytest.raises(OSError, match="simulated disk failure"):
        LDScoreDirectoryWriter().write_batches(producer, LDScoreOutputConfig(output_dir=tmp_path / "out"))
    assert closed == [True]
    assert list((tmp_path / "out").iterdir()) == []


def test_batch_publication_is_private_and_overwrite_removes_owned_siblings(tmp_path):
    first = make_split_ldscore_result()
    second = replace(first, query_columns=["second"], query_table=first.query_table.rename(columns={"query": "second"}),
                     count_records=[{**record, "column": "second" if record["group"] == "query" else record["column"]}
                                    for record in first.count_records], overlap=None)
    first = replace(first, overlap=None)
    output = tmp_path / "out"
    def results():
        yield first
        assert not (output / "metadata.json").exists()
        assert not list(output.glob("ldscore.*.parquet"))
        yield second
    writer = LDScoreDirectoryWriter()
    writer.write_batches(results(), LDScoreOutputConfig(output_dir=output))
    unrelated = output / "ldscore.query.batch-not-owned.parquet"
    unrelated.write_text("keep")
    writer.write_batches(iter([first]), LDScoreOutputConfig(output_dir=output, overwrite=True))
    assert (output / "ldscore.query.parquet").exists()
    assert not (output / "ldscore.query.batch00001.parquet").exists()
    writer.write_batches(iter([make_split_ldscore_result(query=False)]), LDScoreOutputConfig(output_dir=output, overwrite=True))
    assert not (output / "ldscore.query.parquet").exists()
    assert unrelated.read_text() == "keep"
    assert load_ldscore_from_dir(str(output)).query_batches == ()


def test_missing_current_query_manifest_is_rejected(tmp_path):
    LDScoreDirectoryWriter().write(make_split_ldscore_result(), LDScoreOutputConfig(output_dir=tmp_path))
    path = tmp_path / "metadata.json"
    payload = json.loads(path.read_text())
    payload.pop("query_batches")
    path.write_text(json.dumps(payload))
    from ldsc.errors import LDSCInputError
    with pytest.raises(LDSCInputError, match="query_batches"):
        load_ldscore_from_dir(str(tmp_path))
