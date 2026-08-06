from __future__ import annotations

import gzip
from pathlib import Path

import pandas as pd
import pytest

from ldsc.errors import LDSCInputError
from ldsc.gene_list_resolver import GeneCatalog, gene_list_query_name, resolve_gene_list


CATALOG_COLUMNS = [
    "ensgid",
    "gene_name",
    "hg38_chr",
    "hg38_start0",
    "hg38_end",
    "hg38_strand",
    "hg19_chr",
    "hg19_start0",
    "hg19_end",
    "hg19_strand",
]


def _write_catalog(path: Path, rows: list[list[object]]) -> Path:
    pd.DataFrame(rows, columns=CATALOG_COLUMNS).to_csv(path, sep="\t", index=False, compression="gzip")
    return path


def _catalog(tmp_path: Path) -> GeneCatalog:
    path = _write_catalog(
        tmp_path / "catalog.tsv.gz",
        [
            ["ENSG00000000001", "GENEA", "chr1", 10, 20, "+", "chr1", 5, 15, "+"],
            ["ENSG00000000002", "GENEB", "chr2", 30, 45, "-", None, None, None, None],
            ["ENSG00000000003", "DUP", "chr3", 50, 60, "+", "chr3", 40, 50, "+"],
            ["ENSG00000000004", "DUP", "chr4", 70, 80, "+", "chr4", 60, 70, "+"],
            ["ENSG00000000005", "ENSG00000000001", "chr5", 90, 100, "+", "chr5", 80, 90, "+"],
        ],
    )
    return GeneCatalog.load(path)


def test_packaged_catalog_is_valid_and_importable():
    catalog = GeneCatalog.load()

    assert catalog.resource == "protein_coding_genes.tsv.gz"
    assert catalog.release == "GENCODE v49"
    assert len(catalog.frame) == 18401
    assert catalog.frame["ensgid"].is_unique
    assert catalog.content_sha256 == "04e9e1e0603f70af7319d8424e686809e2dd0bbefef44e817ba6345f985f5ff9"


def test_resolve_gene_list_accepts_mixed_ids_versions_names_and_collapses_aliases(tmp_path):
    source = tmp_path / "immune.v2.tsv.gz"
    with gzip.open(source, "wt", encoding="utf-8") as handle:
        handle.write("ENSG00000000001.12\nGENEA\nGENEB\nGENEA\n\n")

    result = resolve_gene_list(source, _catalog(tmp_path), genome_build="hg38", source_ordinal=2)

    assert result.query == "immune.v2"
    assert result.status == "ok"
    assert result.reason == ""
    assert result.canonical_ensembl_ids == ("ENSG00000000001", "ENSG00000000002")
    assert result.intervals == (("1", 10, 20), ("2", 30, 45))
    assert result.counts == {
        "nonblank_input_rows": 4,
        "unique_normalized_input_tokens": 3,
        "repeated_token_rows": 1,
        "matched_input_rows": 4,
        "unique_resolved_canonical_genes": 2,
        "alias_collapsed_rows": 1,
        "blank_rows": 1,
        "unmatched_identifier": 0,
        "invalid_identifier": 0,
        "build_missing": 0,
        "ambiguous_identifier": 0,
        "malformed_input": 0,
        "excluded_gene_region": 0,
    }


def test_resolve_gene_list_reports_every_problem_and_keeps_partial_query(tmp_path):
    source = tmp_path / "partial.txt"
    source.write_text("GENEA\nUNKNOWN\nENSG00000000002.bad\nGENEB\n", encoding="utf-8")

    result = resolve_gene_list(source, _catalog(tmp_path), genome_build="hg19")

    assert result.status == "warning"
    assert result.reason == "partial_resolution"
    assert result.canonical_ensembl_ids == ("ENSG00000000001",)
    assert [(row.line, row.input_gene, row.reason, row.canonical_ensembl_id) for row in result.unresolved] == [
        (2, "UNKNOWN", "unmatched_identifier", None),
        (3, "ENSG00000000002.bad", "invalid_identifier", None),
        (4, "GENEB", "build_missing", "ENSG00000000002"),
    ]


def test_mhc_gene_exclusion_uses_unpadded_half_open_intervals(tmp_path):
    catalog = GeneCatalog.load(
        _write_catalog(
            tmp_path / "mhc-catalog.tsv.gz",
            [
                ["ENSG00000000011", "LEFT", "chr6", 24_999_990, 25_000_000, "+", "chr6", 24_999_990, 25_000_000, "+"],
                ["ENSG00000000012", "OVERLAP_LEFT", "chr6", 24_999_999, 25_000_001, "+", "chr6", 24_999_999, 25_000_001, "+"],
                ["ENSG00000000013", "HLA_TEST", "chr6", 30_000_000, 30_000_100, "+", "chr6", 30_000_000, 30_000_100, "+"],
                ["ENSG00000000014", "RIGHT", "chr6", 35_000_000, 35_000_010, "+", "chr6", 35_000_000, 35_000_010, "+"],
                ["ENSG00000000015", "OTHER_CHROM", "chr7", 30_000_000, 30_000_100, "+", "chr7", 30_000_000, 30_000_100, "+"],
            ],
        )
    )
    source = tmp_path / "genes.txt"
    source.write_text("LEFT\nOVERLAP_LEFT\nHLA_TEST\nRIGHT\nOTHER_CHROM\n", encoding="utf-8")

    result = resolve_gene_list(source, catalog, genome_build="hg19", gene_exclude_regions="mhc")

    assert result.status == "warning"
    assert result.reason == "partial_resolution"
    assert result.canonical_ensembl_ids == (
        "ENSG00000000011",
        "ENSG00000000014",
        "ENSG00000000015",
    )
    assert result.counts["excluded_gene_region"] == 2
    assert [(row.line, row.input_gene, row.reason) for row in result.unresolved] == [
        (2, "OVERLAP_LEFT", "excluded_gene_region"),
        (3, "HLA_TEST", "excluded_gene_region"),
    ]


@pytest.mark.parametrize(
    ("contents", "reason"),
    [
        ("DUP\n", "ambiguous_identifier"),
        ("ENSG00000000001\n", "ambiguous_identifier"),
        ("GENEA\textra\nUNKNOWN\n", "malformed_input"),
    ],
)
def test_structural_gene_list_problems_skip_the_query_after_full_scan(tmp_path, contents, reason):
    source = tmp_path / "bad.list"
    source.write_text(contents, encoding="utf-8")

    result = resolve_gene_list(source, _catalog(tmp_path), genome_build="hg38")

    assert result.status == "skipped"
    assert result.reason == reason
    assert result.unresolved


def test_catalog_rejects_partial_coordinate_triplets(tmp_path):
    path = _write_catalog(
        tmp_path / "bad.tsv.gz",
        [["ENSG00000000001", "GENEA", "chr1", 10, None, "+", "chr1", 5, 15, "+"]],
    )

    with pytest.raises(LDSCInputError, match="coordinate triplet.*partially missing"):
        GeneCatalog.load(path)


@pytest.mark.parametrize(
    ("source", "expected"),
    [
        ("immune_genes.txt.gz", "immune_genes"),
        ("immune.v2.tsv.gz", "immune.v2"),
        ("immune.custom.gz", "immune.custom"),
        ("immune.list", "immune"),
    ],
)
def test_gene_list_query_name_removes_only_recognized_suffixes(source, expected):
    assert gene_list_query_name(source) == expected
