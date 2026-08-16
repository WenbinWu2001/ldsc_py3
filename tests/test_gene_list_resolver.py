from __future__ import annotations

import gzip
from pathlib import Path

import pandas as pd
import pytest

from ldsc import gene_list_resolver
from ldsc.errors import LDSCInputError
from ldsc.gene_list_resolver import (
    AUDIT_COLUMNS,
    CATALOG_ISSUE_COLUMNS,
    SUMMARY_COLUMNS,
    GeneCatalog,
    GeneCatalogValidationError,
    gene_list_query_name,
    resolve_gene_lists,
)


CATALOG_COLUMNS = ["gene_id", "gene_name", "chrom", "start", "end", "genome_build"]


def _write_catalog(path: Path, rows: list[list[object]]) -> Path:
    pd.DataFrame(rows, columns=CATALOG_COLUMNS).to_csv(
        path,
        sep="\t",
        index=False,
        compression="gzip" if path.suffix == ".gz" else None,
    )
    return path


def _catalog(tmp_path: Path) -> GeneCatalog:
    return GeneCatalog.load(
        _write_catalog(
            tmp_path / "catalog.tsv.gz",
            [
                ["ENSG1", "GENEA", "chr1", 11, 20, "GRCh37"],
                ["ENSG2", "DUP", "2", 31, 45, "hg19"],
                ["ENSG3", "DUP", "3", 51, 60, "hg19"],
                ["ENSG4", "ENSG1", "4", 71, 80, "hg19"],
                ["BAD", "BADCOORD", "5", "1.0", 90, "hg19"],
                ["OFF", "OFFCHR", "X", 101, 110, "hg19"],
            ],
        )
    )


def test_catalog_is_required_and_normalizes_one_based_coordinates(tmp_path):
    with pytest.raises(TypeError):
        GeneCatalog.load()

    catalog = _catalog(tmp_path)

    assert catalog.source == "catalog.tsv.gz"
    assert catalog.genome_build == "hg19"
    assert catalog.frame.loc[0, ["gene_id", "chrom", "start", "end", "start0"]].tolist() == [
        "ENSG1",
        "1",
        11,
        20,
        10,
    ]
    assert catalog.issues["reason"].tolist() == [
        "identifier_namespace_conflict",
        "duplicate_gene_name",
        "duplicate_gene_name",
        "identifier_namespace_conflict",
        "invalid_start",
        "invalid_chromosome",
    ]
    assert list(catalog.issues.columns) == list(CATALOG_ISSUE_COLUMNS)


def test_canonical_index_catalog_validation_reports_every_defect(tmp_path):
    path = _write_catalog(
        tmp_path / "bad.tsv",
        [
            ["A", "SHARED", "1", 1, 2, "hg19"],
            ["A", "SHARED", "2", 4, 3, "hg19"],
        ],
    )

    with pytest.raises(GeneCatalogValidationError) as exc_info:
        GeneCatalog.load(path, require_canonical=True)

    issues = exc_info.value.issues
    assert set(issues["reason"]) == {"duplicate_gene_id", "duplicate_gene_name", "end_before_start"}
    assert issues.loc[issues["reason"] == "duplicate_gene_id", "related_catalog_lines"].tolist() == ["2,3", "2,3"]


def test_unique_name_referencing_a_duplicate_catalog_id_is_rejected_with_all_lines(tmp_path):
    catalog = GeneCatalog.load(
        _write_catalog(
            tmp_path / "duplicate-id.tsv",
            [
                ["A", "NAME1", "1", 1, 2, "hg19"],
                ["A", "NAME2", "2", 3, 4, "hg19"],
            ],
        )
    )
    source = tmp_path / "genes.txt"
    source.write_text("NAME1\n", encoding="utf-8")

    result = resolve_gene_lists((source,), catalog, resolution_policy="resolved-only")

    assert result.audit.loc[0, ["match_type", "catalog_lines", "reason"]].tolist() == [
        "gene_name",
        "2,3",
        "catalog_duplicate_id",
    ]


def test_batch_resolution_is_exact_ordered_and_policy_independent(tmp_path):
    focal_a = tmp_path / "pathway-a.txt"
    focal_a.write_text("GENEA\nGENEA\nENSG1.2\nDUP\nENSG1\nBADCOORD\nUNKNOWN\n", encoding="utf-8")
    focal_b = tmp_path / "pathway-b.list.gz"
    with gzip.open(focal_b, "wt", encoding="utf-8") as handle:
        handle.write("ENSG2\n")
    control = tmp_path / "control.tsv"
    control.write_text("ENSG3\nOFFCHR\n", encoding="utf-8")

    strict = resolve_gene_lists(
        (focal_a, focal_b),
        _catalog(tmp_path),
        control_path=control,
        resolution_policy="strict",
    )
    exploratory = resolve_gene_lists(
        (focal_a, focal_b),
        _catalog(tmp_path),
        control_path=control,
        resolution_policy="resolved-only",
    )

    pd.testing.assert_frame_equal(strict.audit, exploratory.audit)
    assert list(strict.audit.columns) == list(AUDIT_COLUMNS)
    assert list(strict.summary.columns) == list(SUMMARY_COLUMNS)
    assert strict.audit[["input_role", "source", "line"]].values.tolist() == [
        ["focal", "pathway-a.txt", 1],
        ["focal", "pathway-a.txt", 2],
        ["focal", "pathway-a.txt", 3],
        ["focal", "pathway-a.txt", 4],
        ["focal", "pathway-a.txt", 5],
        ["focal", "pathway-a.txt", 6],
        ["focal", "pathway-a.txt", 7],
        ["focal", "pathway-b.list.gz", 1],
        ["control", "control.tsv", 1],
        ["control", "control.tsv", 2],
    ]
    assert strict.audit[["input_gene", "disposition", "reason"]].values.tolist() == [
        ["GENEA", "retained", ""],
        ["GENEA", "duplicate", "duplicate_canonical_gene"],
        ["ENSG1.2", "rejected", "unmatched_identifier"],
        ["DUP", "rejected", "ambiguous_gene_name"],
        ["ENSG1", "rejected", "identifier_namespace_conflict"],
        ["BADCOORD", "rejected", "catalog_invalid_coordinates"],
        ["UNKNOWN", "rejected", "unmatched_identifier"],
        ["ENSG2", "retained", ""],
        ["ENSG3", "retained", ""],
        ["OFFCHR", "rejected", "outside_supported_chromosome"],
    ]
    assert strict.has_fatal_gate_a_issues
    assert not exploratory.has_fatal_gate_a_issues
    assert exploratory.selection("focal", 1).intervals == (("1", 10, 20),)
    assert exploratory.selection("control", 0).canonical_gene_ids == ("ENSG3",)
    assert exploratory.summary["unique_resolved_genes"].tolist() == [1, 1, 1]


def test_many_sources_enter_one_bulk_row_resolution_pass(tmp_path, monkeypatch):
    paths = []
    for index in range(120):
        source = tmp_path / f"pathway-{index:03d}.txt"
        source.write_text("GENEA\n" * 25, encoding="utf-8")
        paths.append(source)
    real_resolve_rows = gene_list_resolver._resolve_rows
    calls = 0

    def count_bulk_pass(*args, **kwargs):
        nonlocal calls
        calls += 1
        return real_resolve_rows(*args, **kwargs)

    monkeypatch.setattr(gene_list_resolver, "_resolve_rows", count_bulk_pass)

    result = resolve_gene_lists(paths, _catalog(tmp_path), resolution_policy="resolved-only")

    assert calls == 1
    assert len(result.summary) == 120
    assert len(result.audit) == 3_000
    assert result.summary["unique_resolved_genes"].eq(1).all()


def test_malformed_rows_remain_fatal_under_resolved_only(tmp_path):
    source = tmp_path / "bad.list"
    source.write_text("GENEA\textra\nUNKNOWN\n", encoding="utf-8")

    result = resolve_gene_lists((source,), _catalog(tmp_path), resolution_policy="resolved-only")

    assert result.has_fatal_gate_a_issues
    assert result.audit["reason"].tolist() == ["malformed_input", "unmatched_identifier"]


def test_mhc_exclusion_is_audited_before_duplicate_precedence_without_resolution_failure(tmp_path):
    catalog = GeneCatalog.load(
        _write_catalog(
            tmp_path / "mhc.tsv",
            [["MHC1", "HLA-A", "6", 29_940_000, 29_950_000, "hg19"]],
        )
    )
    source = tmp_path / "mhc-genes.txt"
    source.write_text("HLA-A\nMHC1\n", encoding="utf-8")

    result = resolve_gene_lists(
        (source,),
        catalog,
        gene_exclude_regions="mhc",
    )

    assert result.audit[["disposition", "reason"]].values.tolist() == [
        ["excluded", "excluded_gene_region"],
        ["duplicate", "duplicate_canonical_gene"],
    ]
    assert result.summary.loc[0, ["unique_resolved_genes", "excluded_genes"]].tolist() == [1, 1]
    assert not result.has_fatal_gate_a_issues
    assert result.selection("focal", 1).canonical_gene_ids == ()


def test_duplicate_query_names_are_batched_with_identifier_issues(tmp_path):
    left = tmp_path / "same.txt"
    right_dir = tmp_path / "other"
    right_dir.mkdir()
    right = right_dir / "same.tsv.gz"
    left.write_text("UNKNOWN\n", encoding="utf-8")
    with gzip.open(right, "wt", encoding="utf-8") as handle:
        handle.write("ENSG1\n")

    result = resolve_gene_lists((left, right), _catalog(tmp_path), resolution_policy="resolved-only")

    assert result.has_fatal_gate_a_issues
    assert result.summary["source_reasons"].tolist() == ["duplicate_query_name", "duplicate_query_name"]
    assert result.audit.loc[0, "reason"] == "unmatched_identifier"


def test_repeated_focal_source_is_not_silently_deduplicated(tmp_path):
    source = tmp_path / "same.txt"
    source.write_text("GENEA\n", encoding="utf-8")

    result = resolve_gene_lists(
        (source, source),
        _catalog(tmp_path),
        resolution_policy="resolved-only",
    )

    assert result.summary[["source", "source_ordinal", "source_reasons"]].values.tolist() == [
        ["same.txt", 1, "duplicate_query_name"],
        ["same.txt", 2, "duplicate_query_name"],
    ]
    assert result.audit[["source_ordinal", "input_gene"]].values.tolist() == [
        [1, "GENEA"],
        [2, "GENEA"],
    ]
    assert result.has_fatal_gate_a_issues


def test_focal_globs_expand_lexically_but_unmatched_and_control_globs_are_audited(tmp_path):
    (tmp_path / "set-b.txt").write_text("ENSG2\n", encoding="utf-8")
    (tmp_path / "set-a.txt").write_text("GENEA\n", encoding="utf-8")

    expanded = resolve_gene_lists(
        (tmp_path / "set-*.txt",),
        _catalog(tmp_path),
        resolution_policy="resolved-only",
    )

    assert expanded.summary[["source", "source_ordinal"]].values.tolist() == [
        ["set-a.txt", 1],
        ["set-b.txt", 2],
    ]
    assert not expanded.has_fatal_gate_a_issues

    unmatched = resolve_gene_lists(
        (tmp_path / "missing-*.txt",),
        _catalog(tmp_path),
        resolution_policy="resolved-only",
    )
    assert unmatched.summary.loc[0, ["source_status", "source_reasons"]].tolist() == [
        "error",
        "unreadable_gene_list",
    ]
    assert unmatched.has_fatal_gate_a_issues

    control_glob = resolve_gene_lists(
        (tmp_path / "set-a.txt",),
        _catalog(tmp_path),
        control_path=tmp_path / "set-*.txt",
        resolution_policy="resolved-only",
    )
    assert control_glob.summary.iloc[-1][["input_role", "source_reasons"]].tolist() == [
        "control",
        "unreadable_gene_list",
    ]
    assert control_glob.has_fatal_gate_a_issues


def test_control_path_expands_environment_variables(tmp_path, monkeypatch):
    focal = tmp_path / "focal.txt"
    control = tmp_path / "control.txt"
    focal.write_text("GENEA\n", encoding="utf-8")
    control.write_text("ENSG2\n", encoding="utf-8")
    monkeypatch.setenv("LDSC_TEST_GENE_CONTROL", str(control))

    result = resolve_gene_lists(
        (focal,),
        _catalog(tmp_path),
        control_path="$LDSC_TEST_GENE_CONTROL",
    )

    assert result.summary.iloc[-1][["input_role", "source", "source_status"]].tolist() == [
        "control",
        "control.txt",
        "ok",
    ]


def test_resolved_only_does_not_relax_outside_index_coverage(tmp_path):
    source = tmp_path / "genes.txt"
    source.write_text("GENEA\n", encoding="utf-8")

    result = resolve_gene_lists(
        (source,),
        _catalog(tmp_path),
        resolution_policy="resolved-only",
        index_chromosome_coverage=("22",),
    )

    assert result.audit.loc[0, "reason"] == "outside_index_chromosome_coverage"
    assert result.has_fatal_gate_a_issues


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


def test_missing_catalog_columns_are_structural_and_have_no_resolution_authority(tmp_path):
    path = tmp_path / "catalog.tsv"
    path.write_text("gene_id\tgene_name\nA\tA\n", encoding="utf-8")

    with pytest.raises(GeneCatalogValidationError, match="missing required columns") as exc_info:
        GeneCatalog.load(path)

    assert set(exc_info.value.issues["reason"]) == {
        "missing_required_column",
        "validation_incomplete",
    }


def test_header_only_and_duplicate_header_catalogs_fail_with_repairable_structural_issues(tmp_path):
    header_only = tmp_path / "header-only.tsv"
    header_only.write_text("\t".join(CATALOG_COLUMNS) + "\n", encoding="utf-8")
    duplicate_header = tmp_path / "duplicate-header.tsv"
    duplicate_header.write_text(
        "gene_id\tgene_name\tchrom\tstart\tend\tgenome_build\tgene_id\n"
        "A\tA\t1\t1\t2\thg19\tB\n",
        encoding="utf-8",
    )

    for path in (header_only, duplicate_header):
        with pytest.raises(GeneCatalogValidationError) as exc_info:
            GeneCatalog.load(path)
        assert set(exc_info.value.issues["reason"]) == {
            "catalog_unparseable",
            "validation_incomplete",
        }


@pytest.mark.parametrize("catalog_lines", [[pd.NA, 3], [2, 2], [1, 3]])
def test_embedded_catalog_requires_unique_physical_data_line_numbers(tmp_path, catalog_lines):
    catalog = GeneCatalog.load(
        _write_catalog(
            tmp_path / "catalog.tsv",
            [
                ["A", "GENEA", "1", 1, 2, "hg19"],
                ["B", "GENEB", "2", 3, 4, "hg19"],
            ],
        )
    )
    embedded = catalog.frame.loc[:, [*CATALOG_COLUMNS, "catalog_line"]].copy()
    embedded["catalog_line"] = pd.Series(catalog_lines, dtype="Int64")

    with pytest.raises(GeneCatalogValidationError, match="catalog_line") as exc_info:
        GeneCatalog.from_embedded_frame(embedded)

    assert "catalog_unparseable" in set(exc_info.value.issues["reason"])


def test_empty_embedded_catalog_is_rejected_through_shared_validation():
    empty = pd.DataFrame(columns=[*CATALOG_COLUMNS, "catalog_line"])

    with pytest.raises(GeneCatalogValidationError) as exc_info:
        GeneCatalog.from_embedded_frame(empty)

    assert set(exc_info.value.issues["reason"]) == {
        "catalog_unparseable",
        "validation_incomplete",
    }


def test_live_and_embedded_catalogs_have_identical_resolution_truth(tmp_path):
    live = GeneCatalog.load(
        _write_catalog(
            tmp_path / "catalog.tsv",
            [
                ["G1", "GENE1", "1", 10, 20, "hg19"],
                ["G2", "GENE2", "2", 30, 40, "hg19"],
            ],
        )
    )
    embedded_frame = live.frame.loc[:, CATALOG_COLUMNS].copy()
    embedded_frame["catalog_line"] = live.frame["catalog_line"]
    embedded = GeneCatalog.from_embedded_frame(embedded_frame)
    source = tmp_path / "pathway.txt"
    source.write_text("GENE1\nMISSING\nG2\n", encoding="utf-8")

    live_result = resolve_gene_lists(
        (source,),
        live,
        resolution_policy="resolved-only",
    )
    indexed_result = resolve_gene_lists(
        (source,),
        embedded,
        resolution_policy="resolved-only",
        index_chromosome_coverage=("1", "2"),
    )

    pd.testing.assert_frame_equal(live_result.audit, indexed_result.audit)
    pd.testing.assert_frame_equal(live_result.summary, indexed_result.summary)
