"""Behavioral tests for the LDSC2 munged-sumstats compatibility boundary."""

from __future__ import annotations

import gzip
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from ldsc.config import GlobalConfig, RegressionConfig
from ldsc.errors import LDSCInputError
from ldsc.ldscore_calculator import LDScoreResult
from ldsc.regression_runner import RegressionRunner, _project_legacy_sumstats_to_panel
from ldsc.sumstats_munger import SumstatsTable, load_sumstats


def _panel(frame: pd.DataFrame, mode: str) -> LDScoreResult:
    return LDScoreResult(
        baseline_table=frame,
        query_table=None,
        count_records=[
            {
                "group": "baseline",
                "column": "base",
                "all_reference_snp_count": 100.0,
                "common_reference_snp_count": 80.0,
            }
        ],
        baseline_columns=["base"],
        query_columns=[],
        ld_reference_snps=frozenset(),
        ld_regression_snps=frozenset(frame["SNP"]),
        chromosome_results=[],
        config_snapshot=GlobalConfig(
            snp_identifier=mode,
            genome_build="hg38" if mode.startswith("chr_pos") else None,
        ),
    )


def _legacy(frame: pd.DataFrame, trait: str = "legacy") -> SumstatsTable:
    return SumstatsTable(
        data=frame,
        has_alleles={"A1", "A2"}.issubset(frame.columns),
        source_path=f"{trait}.sumstats.gz",
        trait_name=trait,
        provenance={"source_format": "ldsc2_sumstats", "legacy_ldsc2": True},
        config_snapshot=None,
    )


def test_load_sumstats_marks_only_legacy_text_and_rejects_footerless_parquet(tmp_path: Path) -> None:
    text_path = tmp_path / "trait.sumstats.gz"
    with gzip.open(text_path, "wt", encoding="utf-8") as handle:
        handle.write("SNP A1 A2 Z N\nrs1 A C 2 1000\n")

    table = load_sumstats(text_path)

    assert table.provenance["legacy_ldsc2"] is True
    assert table.config_snapshot is None

    parquet_path = tmp_path / "trait.parquet"
    table.data.to_parquet(parquet_path, index=False)
    with pytest.raises(LDSCInputError, match="missing required LDSC3 identity footer"):
        load_sumstats(parquet_path)


def test_legacy_h2_projects_by_rsid_and_orients_to_one_allele_compatible_panel_row() -> None:
    panel = _panel(
        pd.DataFrame(
            {
                "CHR": ["1", "1", "1"],
                "POS": [10, 11, 20],
                "SNP": ["rs1", "rs1", "rs2"],
                "A1": ["A", "A", "G"],
                "A2": ["C", "G", "T"],
                "regression_ld_scores": [1.0, 1.1, 1.2],
                "base": [1.0, 2.0, 3.0],
            }
        ),
        "chr_pos_allele_aware",
    )
    sumstats = _legacy(
        pd.DataFrame(
            {
                "SNP": ["rs1", "rs2"],
                "A1": ["G", "A"],
                "A2": ["A", "G"],
                "Z": [2.0, 3.0],
                "N": [1000.0, 1000.0],
                "FRQ": [0.25, np.nan],
            }
        )
    )

    dataset = RegressionRunner(regression_config=RegressionConfig()).build_dataset(sumstats, panel)

    # rs1 matches only panel A/G and is swapped, so panel identity is adopted,
    # Z changes sign, and FRQ remains the frequency of the oriented panel A1.
    assert dataset.merged["SNP"].tolist() == ["rs1"]
    assert dataset.merged[["CHR", "POS", "A1", "A2"]].iloc[0].tolist() == ["1", 11, "A", "G"]
    assert dataset.merged.loc[0, "Z"] == -2.0
    assert dataset.merged.loc[0, "FRQ"] == pytest.approx(0.75)
    assert dataset.effective_snp_identifier == "chr_pos_allele_aware"


def test_legacy_h2_allele_unaware_chr_pos_panel_uses_panel_identity_without_reorienting() -> None:
    panel = _panel(
        pd.DataFrame(
            {
                "CHR": ["2", "2"],
                "POS": [100, 200],
                "SNP": ["rs1", "rs2"],
                "regression_ld_scores": [1.0, 1.5],
                "base": [1.0, 2.0],
            }
        ),
        "chr_pos",
    )
    sumstats = _legacy(
        pd.DataFrame(
            {
                "SNP": ["rs2", "rs1"],
                "A1": ["C", "A"],
                "A2": ["A", "C"],
                "Z": [4.0, 2.0],
                "N": [1000.0, 1000.0],
                "FRQ": [0.2, np.nan],
            }
        )
    )

    dataset = RegressionRunner(regression_config=RegressionConfig()).build_dataset(sumstats, panel)

    assert dataset.merged["SNP"].tolist() == ["rs1", "rs2"]
    assert dataset.merged["POS"].tolist() == [100, 200]
    assert dataset.merged["A1"].tolist() == ["A", "C"]
    assert dataset.merged["Z"].tolist() == [2.0, 4.0]
    assert pd.isna(dataset.merged.loc[0, "FRQ"])


def test_legacy_allele_aware_projection_covers_all_four_safe_orientations() -> None:
    panel = _panel(
        pd.DataFrame(
            {
                "CHR": ["1"] * 4,
                "POS": [10, 20, 30, 40],
                "SNP": ["direct", "complement", "swap", "swap_complement"],
                "A1": ["A"] * 4,
                "A2": ["C"] * 4,
                "regression_ld_scores": [1.0, 1.1, 1.2, 1.3],
                "base": [1.0, 2.0, 3.0, 4.0],
            }
        ),
        "rsid_allele_aware",
    )
    sumstats = _legacy(
        pd.DataFrame(
            {
                "SNP": ["direct", "complement", "swap", "swap_complement"],
                "A1": ["A", "T", "C", "G"],
                "A2": ["C", "G", "A", "T"],
                "Z": [1.0, 2.0, 3.0, 4.0],
                "N": [1000.0] * 4,
                "FRQ": [0.1, 0.2, 0.3, 0.4],
            }
        )
    )

    dataset = RegressionRunner(regression_config=RegressionConfig()).build_dataset(sumstats, panel)

    assert dataset.merged["Z"].tolist() == [1.0, 2.0, -3.0, -4.0]
    assert dataset.merged["FRQ"].tolist() == pytest.approx([0.1, 0.2, 0.7, 0.6])
    assert dataset.merged["A1"].tolist() == ["A"] * 4
    assert dataset.merged["A2"].tolist() == ["C"] * 4


def test_legacy_allele_aware_projection_matches_scalar_compatibility_snapshot() -> None:
    """Lock the scalar implementation's finalized output and audit semantics."""
    panel = _panel(
        pd.DataFrame(
            {
                "CHR": ["1"] * 15,
                "POS": [10, 20, 30, 40, 50, 51, 60, 61, 70, 80, 90, 100, 110, 120, 130],
                "SNP": [
                    "direct",
                    "complement",
                    "swap",
                    "swap_complement",
                    "disambiguated",
                    "disambiguated",
                    "ambiguous",
                    "ambiguous",
                    "incompatible",
                    "invalid",
                    "missing_allele",
                    "palindromic",
                    "duplicate",
                    "invalid_frq",
                    "missing_frq",
                ],
                "A1": ["A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A"],
                "A2": ["C", "C", "C", "C", "C", "G", "C", "C", "C", "C", "C", "T", "C", "C", "C"],
                "regression_ld_scores": np.ones(15),
                "base": np.arange(1.0, 16.0),
            }
        ),
        "rsid_allele_aware",
    )
    source = pd.DataFrame(
        {
            "SNP": [
                "direct",
                "complement",
                "swap",
                "swap_complement",
                "disambiguated",
                "ambiguous",
                "incompatible",
                "invalid",
                "missing_allele",
                "palindromic",
                "missing_panel",
                "duplicate",
                "duplicate",
                "invalid_frq",
                "missing_frq",
            ],
            "A1": ["a", "t", "c", "g", "a", "a", "a", "a", None, "a", "a", "a", "a", "a", "a"],
            "A2": ["c", "g", "a", "t", "g", "c", "g", "n", "c", "t", "c", "c", "c", "c", "c"],
            "Z": np.arange(1.0, 16.0),
            "N": np.full(15, 1000.0),
            "FRQ": [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.1, 0.2, 0.3, 0.4, 1.5, None],
            "EXTRA": np.arange(15),
        }
    )
    original = source.copy(deep=True)

    projected, drops = _project_legacy_sumstats_to_panel(_legacy(source), panel)

    expected = pd.DataFrame(
        {
            "SNP": ["direct", "complement", "swap", "swap_complement", "disambiguated", "invalid_frq", "missing_frq"],
            "A1": ["A"] * 7,
            "A2": ["C", "C", "C", "C", "G", "C", "C"],
            "Z": [1.0, 2.0, -3.0, -4.0, 5.0, 14.0, 15.0],
            "N": [1000.0] * 7,
            "FRQ": [0.1, 0.2, 0.7, 0.6, 0.5, np.nan, np.nan],
            "EXTRA": [0, 1, 2, 3, 4, 13, 14],
            "CHR": ["1"] * 7,
            "POS": pd.Series([10, 20, 30, 40, 51, 120, 130], dtype="int64"),
        }
    )
    pd.testing.assert_frame_equal(projected.data, expected)
    audit = pd.DataFrame(drops).reset_index(drop=True)
    assert audit["reason"].tolist() == [
        "ambiguous_panel_mapping",
        "incompatible_alleles",
        "invalid_allele",
        "missing_allele",
        "strand_ambiguous",
        "missing_panel_rsid",
        "duplicate_source_rsid",
        "duplicate_source_rsid",
    ]
    assert audit["panel_candidate_count"].tolist() == [2, 1, 1, 1, 1, 0, 1, 1]
    pd.testing.assert_frame_equal(source, original)


def test_legacy_allele_unaware_projection_requires_one_panel_row_per_rsid() -> None:
    panel = _panel(
        pd.DataFrame(
            {
                "CHR": ["1", "1", "1"],
                "POS": [10, 11, 20],
                "SNP": ["ambiguous", "ambiguous", "good"],
                "regression_ld_scores": [1.0, 1.0, 1.0],
                "base": [1.0, 2.0, 3.0],
            }
        ),
        "rsid",
    )
    projected, drops = _project_legacy_sumstats_to_panel(
        _legacy(
            pd.DataFrame(
                {
                    "SNP": ["ambiguous", "good"],
                    "A1": ["A", "G"],
                    "A2": ["C", "T"],
                    "Z": [1.0, 2.0],
                    "N": [1000.0, 1000.0],
                }
            )
        ),
        panel,
    )

    assert projected.data["SNP"].tolist() == ["good"]
    assert pd.DataFrame(drops)[["SNP", "reason", "panel_candidate_count"]].to_dict("records") == [
        {"SNP": "ambiguous", "reason": "ambiguous_panel_mapping", "panel_candidate_count": 2}
    ]


def test_legacy_projection_drops_duplicate_source_clusters_and_records_stable_reasons() -> None:
    panel = _panel(
        pd.DataFrame(
            {
                "CHR": ["1", "1", "1"],
                "POS": [10, 20, 30],
                "SNP": ["dup", "good", "ambiguous"],
                "regression_ld_scores": [1.0, 1.0, 1.0],
                "base": [1.0, 2.0, 3.0],
            }
        ),
        "rsid",
    )
    sumstats = _legacy(
        pd.DataFrame(
            {
                "SNP": ["dup", "dup", "good", "missing"],
                "A1": ["A", "A", "A", "A"],
                "A2": ["C", "C", "C", "C"],
                "Z": [1.0, 2.0, 3.0, 4.0],
                "N": [1000.0] * 4,
            }
        )
    )

    dataset = RegressionRunner(regression_config=RegressionConfig()).build_dataset(sumstats, panel)

    assert dataset.merged["SNP"].tolist() == ["good"]
    reasons = sumstats.provenance.get("legacy_sumstats_drops", [])
    # Projection returns a new working table; the caller-owned table is never mutated.
    assert reasons == []
    projected = dataset.legacy_sumstats_drops
    assert projected["reason"].value_counts().to_dict() == {
        "duplicate_source_rsid": 2,
        "missing_panel_rsid": 1,
    }


def test_legacy_rg_keeps_missing_frq_and_harmonizes_trait_two_without_panel_alleles() -> None:
    panel = _panel(
        pd.DataFrame(
            {
                "CHR": ["1", "1"],
                "POS": [10, 20],
                "SNP": ["rs1", "rs2"],
                "regression_ld_scores": [1.0, 1.0],
                "base": [1.0, 2.0],
            }
        ),
        "rsid",
    )
    trait_1 = _legacy(
        pd.DataFrame(
            {
                "SNP": ["rs1", "rs2"],
                "A1": ["A", "G"],
                "A2": ["C", "T"],
                "Z": [2.0, 3.0],
                "N": [1000.0, 1000.0],
                "FRQ": [np.nan, 0.4],
            }
        ),
        "one",
    )
    trait_2 = _legacy(
        pd.DataFrame(
            {
                "SNP": ["rs1", "rs2"],
                "A1": ["C", "G"],
                "A2": ["A", "T"],
                "Z": [4.0, 5.0],
                "N": [900.0, 900.0],
                "FRQ": [0.2, np.nan],
            }
        ),
        "two",
    )

    dataset = RegressionRunner(regression_config=RegressionConfig()).build_rg_dataset(trait_1, trait_2, panel)

    assert dataset.merged["SNP"].tolist() == ["rs1", "rs2"]
    assert dataset.merged["Z2"].tolist() == [-4.0, 5.0]
    assert dataset.merged.loc[0, "FRQ2"] == pytest.approx(0.8)
    assert pd.isna(dataset.merged.loc[1, "FRQ2"])


def test_legacy_sumstats_without_alleles_are_rejected_for_regression() -> None:
    panel = _panel(
        pd.DataFrame(
            {
                "CHR": ["1", "1"],
                "POS": [10, 20],
                "SNP": ["rs1", "rs2"],
                "regression_ld_scores": [1.0, 1.0],
                "base": [1.0, 2.0],
            }
        ),
        "rsid",
    )
    sumstats = _legacy(pd.DataFrame({"SNP": ["rs1"], "Z": [1.0], "N": [1000.0]}))

    with pytest.raises(LDSCInputError, match="legacy LDSC2 sumstats.*A1.*A2"):
        RegressionRunner(regression_config=RegressionConfig()).build_dataset(sumstats, panel)
