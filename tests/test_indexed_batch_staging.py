"""Indexed query staging is bounded by the active batch, not all queries."""

from pathlib import Path
from dataclasses import replace
from types import SimpleNamespace
import weakref

import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from ldsc import _indexed_ldscore_batches as batches
from tests.test_gene_index_streaming import two_chromosome_index
from tests.test_gene_ldscore_index import run_indexed_ldscore


def _check_worker_serialization_lifetimes(*args, **kwargs):
    from ldsc import gene_ldscore_index as core

    original_load = core.LoadedGeneLDScoreIndex.load_chromosome
    original_assemble = core.assemble_indexed_ld_scores
    original_write = batches.pq.write_table
    arrays, operators = [], []

    def load(*values, **options):
        record = original_load(*values, **options)
        operators.append(weakref.ref(record.operator))
        return record

    def assemble(*values, **options):
        scores = original_assemble(*values, **options)
        arrays.append(weakref.ref(scores))
        return scores

    def write(table, path, **options):
        if Path(path).name == "query.parquet":
            assert all(ref() is None for ref in arrays), "duplicate dense scores survive into serialization"
            assert all(ref() is None for ref in operators), "chromosome operator survives into serialization"
        return original_write(table, path, **options)

    core.LoadedGeneLDScoreIndex.load_chromosome = load
    core.assemble_indexed_ld_scores = assemble
    batches.pq.write_table = write
    try:
        return batches._write_chromosome(*args, **kwargs)
    finally:
        core.LoadedGeneLDScoreIndex.load_chromosome = original_load
        core.assemble_indexed_ld_scores = original_assemble
        batches.pq.write_table = original_write


def test_worker_releases_operator_and_duplicate_scores_before_serialization(tmp_path, monkeypatch):
    index = two_chromosome_index(tmp_path)
    queries = [tmp_path / "first.txt", tmp_path / "second.txt"]
    for query, genes in zip(queries, ["G1\nG3\n", "G2\nG4\n"], strict=True):
        query.write_text(genes)
    monkeypatch.setattr(batches, "_write_chromosome", _check_worker_serialization_lifetimes)
    result = run_indexed_ldscore(index, query_gene_list_sources=queries, output_dir=tmp_path / "out", query_batch_size=2)
    assert list(result.query_columns) == ["first", "second"]


def test_indexed_staging_holds_one_batch_and_reuses_pool(tmp_path, monkeypatch):
    index = two_chromosome_index(tmp_path)
    queries = []
    for ordinal in range(6):
        query = tmp_path / f"query{ordinal}.txt"
        query.write_text("G1\nG3\n" if ordinal % 2 == 0 else "G2\nG4\n")
        queries.append(query)
    original_results = batches.indexed_results
    original_pool = batches.ProcessPoolExecutor
    pools, pending_counts, previous = [], [], []

    def pool(*args, **kwargs):
        pools.append(kwargs["max_workers"])
        return original_pool(*args, **kwargs)

    def results(*args, **kwargs):
        producer = original_results(*args, **kwargs)
        try:
            for result in producer:
                assert all(not path.exists() for path in previous)
                root = Path(args[5]) / "indexed-batches"
                pending = [path for path in root.rglob("query*") if path.is_file()]
                assert len(pending) <= 2, "future query batches accumulated before consumption"
                assert getattr(result, "query_table", None) is None, "genome-wide query array remains live"
                pending_counts.append(len(pending))
                previous[:] = pending
                yield result
        finally:
            producer.close()

    monkeypatch.setattr(batches, "ProcessPoolExecutor", pool)
    monkeypatch.setattr(batches, "indexed_results", results)
    output = tmp_path / "out"
    result = run_indexed_ldscore(index, query_gene_list_sources=queries,
                                output_dir=output, query_batch_size=2, workers=2)
    assert pools == [2]
    assert pending_counts == [2, 2, 2]
    assert all(not path.exists() for path in previous)
    columns = [query.stem for query in queries]
    # G1/G3 select atoms 0,1; G2/G4 select atoms 1,2 in the independent fixture.
    expected = np.array([[.75, -.25], [1.5, 1.1]] * 2, dtype=np.float32)
    np.testing.assert_array_equal(result.read_queries(columns)[columns], np.tile(expected, (1, 3)))
    assert not list(output.glob(".ldsc-*"))


@pytest.mark.parametrize("split", [0, 1, 3, 4])
def test_variation_summaries_preserve_global_float64_and_nan_semantics(split):
    values = np.array([[1., 1., np.nan, np.nan, 0., np.inf],
                       [1., 1., np.nan, np.nan, -0., np.inf],
                       [2., 1. + 2.**-30, np.nan, 0., 0., np.inf],
                       [2., 1. + 2.**-30, np.nan, 0., -0., np.inf]])
    summaries = [SimpleNamespace(first=first, varying=varying)
                 for first, varying in (batches._score_variation(part)
                                        for part in (values[:split], values[split:]))]
    np.testing.assert_array_equal(batches._combined_variation(summaries, values.shape[1]),
                                  [True, True, False, True, False, False])


@pytest.mark.parametrize("width", [1, 2])
def test_query_pruning_precedes_float32_narrowing_and_allows_skipped_batch(tmp_path, width):
    from tests.test_gene_ldscore_index import _artifact_payload, publish_gene_ldscore_index

    record, catalog, identity = _artifact_payload()
    # G1 selects atoms 0,1 and varies only below float32 precision.
    # G2 selects atoms 1,2 and is exactly constant in float64.
    record = replace(record, operator=sparse.csr_matrix([[1., 0., 1.], [1. + 2.**-30, 0., 1.]]))
    index = publish_gene_ldscore_index(tmp_path / "index", index_identity=identity,
                                      gene_catalog=catalog, chromosomes={"22": record}, overwrite=False)
    constant, varying = tmp_path / "constant.txt", tmp_path / "varying.txt"
    constant.write_text("G2\n")
    varying.write_text("G1\n")
    result = run_indexed_ldscore(index, query_gene_list_sources=[constant, varying],
                                output_dir=tmp_path / "out", query_batch_size=width)
    assert list(result.query_columns) == ["varying"]
    assert [(s.query, s.status, s.reason) for s in result.query_statuses] == [
        ("constant", "skipped", "zero_variance_ld_scores"), ("varying", "ok", "")]
    assert len(result.query_batches) == 1
    actual = result.read_queries(["varying"]).varying
    assert actual.dtype == np.float32
    np.testing.assert_array_equal(actual, [1., 1.])
    assert [r["column"] for r in result.count_records] == ["base", "varying"]
    assert list(result.overlap.query_diagonal_all.index) == ["varying"]


def test_fragment_writer_releases_each_chromosome_before_reading_next(tmp_path, monkeypatch):
    from ldsc import _ldscore_fragments as fragments

    index = two_chromosome_index(tmp_path)
    query = tmp_path / "query.txt"
    query.write_text("G1\nG3\n")
    original = fragments.pd.read_parquet
    live = []

    def read(path, *args, **kwargs):
        if Path(path).name == "query.parquet":
            assert all(ref() is None for ref in live), "previous chromosome query frame remains live"
            frame = original(path, *args, **kwargs)
            assert len(frame) == 2
            live.append(weakref.ref(frame))
            return frame
        return original(path, *args, **kwargs)

    monkeypatch.setattr(fragments.pd, "read_parquet", read)
    result = run_indexed_ldscore(index, query_gene_list_sources=[query], output_dir=tmp_path / "out")
    assert len(live) == 2 and all(ref() is None for ref in live)
    assert [entry["chrom"] for entry in result.query_batches[0]["row_groups"]] == ["21", "22"]


def test_fragment_row_mismatch_aborts_publication_and_cleans_scratch(tmp_path, monkeypatch):
    import multiprocessing
    from ldsc.errors import LDSCInputError

    index = two_chromosome_index(tmp_path)
    query = tmp_path / "query.txt"
    query.write_text("G1\nG3\n")
    original = batches.indexed_results

    def corrupt(*args, **kwargs):
        producer = original(*args, **kwargs)
        try:
            for result in producer:
                path = result.query_fragments[-1]
                frame = pd.read_parquet(path)
                frame.loc[0, "SNP"] = "wrong-row"
                frame.to_parquet(path, index=False)
                yield result
        finally:
            producer.close()

    monkeypatch.setattr(batches, "indexed_results", corrupt)
    output = tmp_path / "out"
    before = {p.pid for p in multiprocessing.active_children()}
    with pytest.raises(LDSCInputError, match="SNP row mismatch"):
        run_indexed_ldscore(index, query_gene_list_sources=[query], output_dir=output, workers=2)
    assert {p.pid for p in multiprocessing.active_children()} == before
    assert not (output / "metadata.json").exists()
    assert not list(output.glob("ldscore.*.parquet"))
    assert not list(output.glob(".ldsc-*"))


def test_failed_batch_waits_for_running_chromosome_before_scratch_cleanup(tmp_path, monkeypatch):
    from concurrent.futures import ThreadPoolExecutor
    from threading import Event

    index = two_chromosome_index(tmp_path)
    query = tmp_path / "query.txt"
    query.write_text("G1\nG3\n")
    running, shutting_down, finished = Event(), Event(), Event()

    class Pool(ThreadPoolExecutor):
        def __init__(self, max_workers, **kwargs):
            super().__init__(max_workers=max_workers)

        def shutdown(self, *args, **kwargs):
            shutting_down.set()
            return super().shutdown(*args, **kwargs)

    def task(*args):
        if args[1] == "21":
            assert running.wait(5)
            raise OSError("failed chromosome")
        running.set()
        assert shutting_down.wait(5)
        assert Path(args[7]).is_dir(), "batch scratch disappeared while a worker was running"
        finished.set()

    monkeypatch.setattr(batches, "ProcessPoolExecutor", Pool)
    monkeypatch.setattr(batches, "_write_chromosome", task)
    output = tmp_path / "out"
    with pytest.raises(OSError, match="failed chromosome"):
        run_indexed_ldscore(index, query_gene_list_sources=[query], output_dir=output, workers=2)
    assert finished.is_set()
    assert not list(output.glob(".ldsc-*"))
