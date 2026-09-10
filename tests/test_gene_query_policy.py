"""Shared production policy distinguishes missing coverage from measured support."""

import pandas as pd
import pytest

from ldsc.gene_list_resolver import GeneCatalog, resolve_gene_lists
from ldsc.query_annotations import assess_gene_coverage, gene_query_statuses


@pytest.mark.parametrize("policy", ["strict", "resolved-only"])
def test_coverage_uses_unique_selected_genes_and_leaves_support_unknown(tmp_path, policy):
    catalog = tmp_path / "genes.tsv"
    catalog.write_text(
        "gene_id\tgene_name\tchrom\tstart\tend\tgenome_build\n"
        "G1\tONE\t22\t10\t20\thg19\n"
        "G2\tTWO\t1\t10\t20\thg19\n"
    )
    sources = []
    for name, content in [("full", "G1\nONE\n"), ("partial", "G1\nG2\n"),
                          ("none", "G2\n"), ("empty", "")]:
        path = tmp_path / f"{name}.txt"
        path.write_text(content)
        sources.append(path)
    batch = resolve_gene_lists(sources, GeneCatalog.load(catalog), resolution_policy=policy)
    batch, errors = assess_gene_coverage(batch, ["22"])
    assert batch.summary.coverage_status.tolist() == ["full", "partial", "none", "empty"]
    assert batch.summary.selected_genes.tolist() == [1, 2, 1, 0]
    assert batch.summary.covered_genes.tolist() == [1, 1, 0, 0]
    assert len(errors) == 2
    assert batch.summary.zero_support_genes.isna().all()
    assert batch.audit.loc[batch.audit.canonical_gene_id.eq("G2"), "coverage_status"].tolist() == ["uncovered", "uncovered"]
    assert not batch.audit.disposition.eq("unsupported").any()
    statuses = gene_query_statuses(batch)
    assert [s.reason for s in statuses] == ["", "incomplete_chromosome_coverage", "incomplete_chromosome_coverage", "empty_gene_list"]
    assert [s.status for s in statuses] == ["ok", "error", "error", "skipped"]


def test_missing_support_evidence_is_unknown_not_zero(tmp_path):
    catalog = tmp_path / "genes.tsv"
    catalog.write_text("gene_id\tgene_name\tchrom\tstart\tend\tgenome_build\nG1\tONE\t22\t10\t20\thg19\n")
    source = tmp_path / "focal.txt"
    source.write_text("G1\n")
    batch = resolve_gene_lists([source], GeneCatalog.load(catalog)).with_snp_support(pd.Series(dtype="Int64"))
    assert batch.audit.disposition.tolist() == ["retained"]
    assert pd.isna(batch.summary.loc[0, "zero_support_genes"])
    assert pd.isna(batch.summary.loc[0, "genes_with_snp_support"])


def test_exclusions_precede_coverage_and_all_excluded_is_empty(tmp_path):
    catalog = tmp_path / "genes.tsv"
    catalog.write_text("gene_id\tgene_name\tchrom\tstart\tend\tgenome_build\n"
                       "MHC\tMHCGENE\t6\t30000000\t30000100\thg19\n"
                       "G1\tONE\t22\t10\t20\thg19\n")
    excluded, full = tmp_path / "excluded.txt", tmp_path / "full.txt"
    excluded.write_text("MHC\n")
    full.write_text("MHC\nG1\nONE\n")
    batch = resolve_gene_lists([excluded, full], GeneCatalog.load(catalog), gene_exclude_regions="mhc")
    batch, errors = assess_gene_coverage(batch, ["22"])
    assert not errors
    assert batch.summary.coverage_status.tolist() == ["empty", "full"]
    assert batch.summary.selected_genes.tolist() == [0, 1]
    assert batch.summary.excluded_genes.tolist() == [1, 1]
    assert gene_query_statuses(batch)[0].reason == "zero_resolved_genes"
    assert batch.audit.loc[batch.audit.input_gene.eq("MHC"), "coverage_status"].isna().all()
