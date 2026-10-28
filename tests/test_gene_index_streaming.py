"""Index integrity and assembly retain only one chromosome operator."""

from dataclasses import replace
import weakref

import numpy as np
import pandas as pd
import pytest

from ldsc._ldscore.gene_index import scoring as index_scoring, build as index_module
from ldsc._ldscore.gene_index import artifact as _gene_index_artifact
from ldsc._ldscore.gene_index.artifact import _catalog_identity_records
from tests.test_gene_ldscore_index import (
    _artifact_payload, load_gene_ldscore_index, publish_gene_ldscore_index,
    run_indexed_ldscore,
)


def two_chromosome_index(tmp_path):
    record, catalog, identity = _artifact_payload()
    first_catalog = catalog.assign(
        chrom="21", gene_id=["ENSG3", "ENSG4"], gene_name=["G3", "G4"],
    )
    catalog = pd.concat([first_catalog, catalog], ignore_index=True)
    catalog["gene_index"] = np.arange(4)
    catalog["catalog_line"] = np.arange(2, 6)
    identity = {**identity, "chromosomes": ["21", "22"],
                "catalog": _catalog_identity_records(catalog)}
    first = replace(record, baseline_rows=record.baseline_rows.assign(
        CHR="21", SNP=["rs3", "rs4"]), atom_model=replace(record.atom_model, chromosome="21"))
    return publish_gene_ldscore_index(
        tmp_path / "index", index_identity=identity, gene_catalog=catalog,
        chromosomes={"21": first, "22": record}, overwrite=False,
    )


def test_validated_index_releases_operators_and_keeps_support(tmp_path, monkeypatch):
    path = two_chromosome_index(tmp_path)
    original = _gene_index_artifact._load_index_chromosome
    live = []

    def checked_load(*args, **kwargs):
        assert all(ref() is None for ref in live), "previous chromosome operator remains live"
        record = original(*args, **kwargs)
        live.append(weakref.ref(record.operator))
        return record

    monkeypatch.setattr(index_scoring, "_load_index_chromosome", checked_load)
    index = load_gene_ldscore_index(path)
    assert all(ref() is None for ref in live)
    assert index.chromosomes == ("21", "22")
    np.testing.assert_array_equal(index.gene_support, [2, 1, 2, 1])
    assert not hasattr(index, "index_chromosomes")


def _observed_index_task(*args, **kwargs):
    import json
    import os
    from pathlib import Path
    from ldsc._ldscore.gene_index import batches
    original_load = _gene_index_artifact._load_index_chromosome
    original_assemble = batches.assemble_indexed_ld_scores
    loads, widths, live = [], [], []
    def load(*values, **options):
        assert all(ref() is None for ref in live)
        record = original_load(*values, **options)
        loads.append(values[1])
        live.append(weakref.ref(record.operator))
        return record
    def assemble(operator, selectors):
        if selectors.ndim == 2:
            widths.append(selectors.shape[1])
        return original_assemble(operator, selectors)
    index_scoring._load_index_chromosome = load
    batches.assemble_indexed_ld_scores = assemble
    try:
        result = batches._write_chromosome(*args, **kwargs)
        assert all(ref() is None for ref in live)
        (Path(os.environ["LDSC_TEST_INDEX_EVENTS"]) / f"{args[1]}-{args[5]}.json").write_text(
            json.dumps({"loads": loads, "widths": widths}))
        return result
    finally:
        index_scoring._load_index_chromosome = original_load
        batches.assemble_indexed_ld_scores = original_assemble


@pytest.mark.parametrize("batch_size", [1, 2, 1000])
def test_index_assembly_reloads_and_releases_each_operator_per_batch(tmp_path, monkeypatch, batch_size):
    path = two_chromosome_index(tmp_path)
    queries = []
    for name, genes in [("first", "G1\nG3\n"), ("second", "G2\nG4\n"),
                        ("union", "G1\nG2\nG3\nG4\n")]:
        query = tmp_path / f"{name}.txt"
        query.write_text(genes)
        queries.append(query)
    original_load = _gene_index_artifact._load_index_chromosome
    loads, live = [], []
    def load(*args, **kwargs):
        assert all(ref() is None for ref in live)
        record = original_load(*args, **kwargs)
        live.append(weakref.ref(record.operator))
        loads.append(args[1])
        return record
    from ldsc._ldscore.gene_index import batches
    events = tmp_path / "events"
    events.mkdir()
    monkeypatch.setenv("LDSC_TEST_INDEX_EVENTS", str(events))
    monkeypatch.setattr(index_scoring, "_load_index_chromosome", load)
    monkeypatch.setattr(batches, "_write_chromosome", _observed_index_task)
    result = run_indexed_ldscore(path, query_gene_list_sources=queries,
                                output_dir=tmp_path / "out", query_batch_size=batch_size)
    assert loads == ["21", "22"]
    import json
    total_batches = (len(queries) + batch_size - 1) // batch_size
    observed = [json.loads((events / f"{chrom}-{ordinal}.json").read_text())
                for ordinal in range(total_batches) for chrom in loads]
    assert [item["loads"] for item in observed] == [["21"], ["22"]] * total_batches
    widths = [width for item in observed for width in item["widths"]]
    assert all(ref() is None for ref in live)
    assert max(widths) <= batch_size
    # Stored atoms: G1/G3 select atoms 0,1; G2/G4 select 1,2.
    np.testing.assert_allclose(result.read_queries(["first", "second", "union"])[["first", "second", "union"]],
                               [[.75, -.25, .75], [1.5, 1.1, 1.6]] * 2)
    assert result.baseline_table.SNP.tolist() == ["rs3", "rs4", "rs1", "rs2"]
    assert result.gene_list_batch.audit_path.is_file()


@pytest.mark.parametrize("workers", [1, 2, -1, -2, 10])
def test_parallel_index_batches_preserve_union_scores_and_counts(tmp_path, workers, monkeypatch):
    monkeypatch.setattr("joblib.cpu_count", lambda: 16)
    path = two_chromosome_index(tmp_path)
    queries = []
    for name, genes in [("first", "G1\nG3\n"), ("second", "G2\nG4\n"),
                        ("union", "G1\nG2\nG3\nG4\n")]:
        query = tmp_path / f"{name}.txt"
        query.write_text(genes)
        queries.append(query)
    reference = run_indexed_ldscore(path, query_gene_list_sources=queries,
                                   output_dir=tmp_path / "reference", query_batch_size=1000)
    result = run_indexed_ldscore(path, query_gene_list_sources=queries,
                                output_dir=tmp_path / "out", query_batch_size=1, workers=workers)
    assert not hasattr(result, "query_table")
    assert len(result.query_batches) == 3
    columns = ["union", "second", "first"]
    pd.testing.assert_frame_equal(result.read_queries(columns), reference.read_queries(columns))
    np.testing.assert_allclose(result.read_queries(["first", "second", "union"])[["first", "second", "union"]],
                               [[.75, -.25, .75], [1.5, 1.1, 1.6]] * 2, rtol=1e-6, atol=1e-8)
    assert result.count_records == reference.count_records
    pd.testing.assert_frame_equal(result.overlap.baseline_block_all, reference.overlap.baseline_block_all)
    pd.testing.assert_frame_equal(result.overlap.baseline_block_common, reference.overlap.baseline_block_common)
    assert not list((tmp_path / "out").glob(".ldsc-*"))
    assert not list((tmp_path / "out").glob(".ldsc-annotation-*"))


@pytest.mark.parametrize("batch_size", [0, -1, 1.5, True])
def test_index_batch_size_validation_precedes_io(tmp_path, batch_size):
    with pytest.raises(ValueError, match="query_batch_size"):
        run_indexed_ldscore(tmp_path / "missing", query_gene_list_sources=[],
                            output_dir=tmp_path / "out", query_batch_size=batch_size)
    assert not (tmp_path / "out").exists()
