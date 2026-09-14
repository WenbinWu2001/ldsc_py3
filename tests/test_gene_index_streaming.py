"""Index integrity and assembly retain only one chromosome operator."""

from dataclasses import replace
import weakref

import numpy as np
import pandas as pd
import pytest

from ldsc import gene_ldscore_index as index_module
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
                "catalog": index_module._catalog_identity_records(catalog)}
    first = replace(record, baseline_rows=record.baseline_rows.assign(
        CHR="21", SNP=["rs3", "rs4"]), atom_model=replace(record.atom_model, chromosome="21"))
    return publish_gene_ldscore_index(
        tmp_path / "index", index_identity=identity, gene_catalog=catalog,
        chromosomes={"21": first, "22": record}, overwrite=False,
    )


def test_validated_index_releases_operators_and_keeps_support(tmp_path, monkeypatch):
    path = two_chromosome_index(tmp_path)
    original = index_module._load_index_chromosome
    live = []

    def checked_load(*args, **kwargs):
        assert all(ref() is None for ref in live), "previous chromosome operator remains live"
        record = original(*args, **kwargs)
        live.append(weakref.ref(record.operator))
        return record

    monkeypatch.setattr(index_module, "_load_index_chromosome", checked_load)
    index = load_gene_ldscore_index(path)
    assert all(ref() is None for ref in live)
    assert index.chromosomes == ("21", "22")
    np.testing.assert_array_equal(index.gene_support, [2, 1, 2, 1])
    assert not hasattr(index, "index_chromosomes")


@pytest.mark.parametrize("batch_size", [1, 2, 1000])
def test_index_assembly_loads_each_operator_twice_and_batches_queries(tmp_path, monkeypatch, batch_size):
    path = two_chromosome_index(tmp_path)
    queries = []
    for name, genes in [("first", "G1\nG3\n"), ("second", "G2\nG4\n"),
                        ("union", "G1\nG2\nG3\nG4\n")]:
        query = tmp_path / f"{name}.txt"
        query.write_text(genes)
        queries.append(query)
    original_load = index_module._load_index_chromosome
    original_assemble = index_module.assemble_indexed_ld_scores
    loads, live, widths = [], [], []

    def load(*args, **kwargs):
        assert all(ref() is None for ref in live)
        record = original_load(*args, **kwargs)
        live.append(weakref.ref(record.operator))
        loads.append(args[1])
        return record

    def assemble(operator, selectors):
        if selectors.ndim == 2:
            widths.append(selectors.shape[1])
        return original_assemble(operator, selectors)

    monkeypatch.setattr(index_module, "_load_index_chromosome", load)
    monkeypatch.setattr(index_module, "assemble_indexed_ld_scores", assemble)
    result = run_indexed_ldscore(path, query_gene_list_sources=queries,
                                output_dir=tmp_path / "out", query_batch_size=batch_size)
    assert loads == ["21", "22", "21", "22"]
    assert all(ref() is None for ref in live)
    assert max(widths) <= batch_size
    # Stored atoms: G1/G3 select atoms 0,1; G2/G4 select 1,2.
    np.testing.assert_allclose(result.read_queries(["first", "second", "union"])[["first", "second", "union"]],
                               [[.75, -.25, .75], [1.5, 1.1, 1.6]] * 2)
    assert result.baseline_table.SNP.tolist() == ["rs3", "rs4", "rs1", "rs2"]
    assert result.gene_list_batch.audit_path.is_file()


@pytest.mark.parametrize("threads", [1, 2, -1, -2, 10])
def test_parallel_index_batches_preserve_union_scores_and_counts(tmp_path, threads):
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
                                output_dir=tmp_path / "out", query_batch_size=1, threads=threads)
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
