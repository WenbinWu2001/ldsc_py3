"""Behavioral tests for the legacy munged-sumstats compatibility boundary."""

from __future__ import annotations

import gzip
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from ldsc.config import GlobalConfig, RegressionConfig
from ldsc.errors import LDSCInputError
from ldsc._ldscore.calculator import LDScoreResult
from ldsc._regression.runner import RegressionRunner, _project_legacy_sumstats_to_panel
from ldsc._sumstats.munger import SumstatsTable, load_sumstats


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
        provenance={"source_format": "legacy_sumstats", "legacy_ldsc": True},
        config_snapshot=None,
    )


def test_load_sumstats_marks_only_legacy_text_and_rejects_footerless_parquet(tmp_path: Path) -> None:
    text_path = tmp_path / "trait.sumstats.gz"
    with gzip.open(text_path, "wt", encoding="utf-8") as handle:
        handle.write("SNP A1 A2 Z N\nrs1 A C 2 1000\n")

    table = load_sumstats(text_path)

    assert table.provenance["legacy_ldsc"] is True
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

    # The palindromic (A/T) row is kept with its source alleles, Z and FRQ:
    # its strand against the panel is unknown, so it is never reoriented.
    expected = pd.DataFrame(
        {
            "SNP": ["direct", "complement", "swap", "swap_complement", "disambiguated", "palindromic", "invalid_frq",
                    "missing_frq"],
            "A1": ["A"] * 8,
            "A2": ["C", "C", "C", "C", "G", "T", "C", "C"],
            "Z": [1.0, 2.0, -3.0, -4.0, 5.0, 10.0, 14.0, 15.0],
            "N": [1000.0] * 8,
            "FRQ": [0.1, 0.2, 0.7, 0.6, 0.5, 0.1, np.nan, np.nan],
            "EXTRA": [0, 1, 2, 3, 4, 9, 13, 14],
            "CHR": ["1"] * 8,
            "POS": pd.Series([10, 20, 30, 40, 51, 100, 120, 130], dtype="int64"),
        }
    )
    pd.testing.assert_frame_equal(projected.data, expected)
    audit = pd.DataFrame(drops).reset_index(drop=True)
    assert audit["reason"].tolist() == [
        "ambiguous_panel_mapping",
        "incompatible_alleles",
        "invalid_allele",
        "missing_allele",
        "missing_panel_rsid",
        "duplicate_source_rsid",
        "duplicate_source_rsid",
    ]
    assert audit["panel_candidate_count"].tolist() == [2, 1, 1, 1, 0, 1, 1]
    pd.testing.assert_frame_equal(source, original)


def test_legacy_h2_with_a_base_mode_panel_keeps_a_t_and_c_g_rows() -> None:
    panel = _panel(
        pd.DataFrame(
            {
                "CHR": ["1", "1", "1"],
                "POS": [10, 20, 30],
                "SNP": ["rs1", "rs2", "rs3"],
                "regression_ld_scores": [1.0, 1.1, 1.2],
                "base": [1.0, 2.0, 3.0],
            }
        ),
        "chr_pos",
    )
    sumstats = _legacy(
        pd.DataFrame(
            {
                "SNP": ["rs1", "rs2", "rs3"],
                "A1": ["A", "C", "A"],
                "A2": ["T", "G", "G"],
                "Z": [2.0, 3.0, 4.0],
                "N": [1000.0] * 3,
            }
        )
    )

    projected, drops = _project_legacy_sumstats_to_panel(sumstats, panel)

    # h2 fits Z^2, so an unknown strand is no reason to drop a SNP.
    assert projected.data["SNP"].tolist() == ["rs1", "rs2", "rs3"]
    assert projected.data["Z"].tolist() == [2.0, 3.0, 4.0]
    assert drops.empty


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
    # The drop reasons are checked through the h2 command's audit file below.
    assert reasons == []


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


@pytest.mark.parametrize("panel_mode,downgrade", [("rsid_allele_aware", True), ("rsid", False)])
def test_legacy_sumstats_without_alleles_project_by_rsid(panel_mode: str, downgrade: bool) -> None:
    frame = {
        "CHR": ["1", "1", "1"], "POS": [10, 20, 30], "SNP": ["rs1", "rs2", "rs3"],
        "regression_ld_scores": [1.0, 1.0, 1.0], "base": [1.0, 2.0, 3.0],
    }
    if panel_mode.endswith("allele_aware"):
        frame.update(A1=["A", "A", "C"], A2=["G", "C", "T"])
    sumstats = _legacy(pd.DataFrame({"SNP": ["rs1", "rs2", "rs3"], "Z": [1.0, -2.0, 0.5], "N": [1000.0] * 3}))

    dataset = RegressionRunner(regression_config=RegressionConfig()).build_dataset(sumstats, _panel(pd.DataFrame(frame), panel_mode))

    assert dataset.merged["SNP"].tolist() == ["rs1", "rs2", "rs3"]
    np.testing.assert_allclose(dataset.merged["Z"], [1.0, -2.0, 0.5])
    assert dataset.effective_snp_identifier == "rsid"
    assert dataset.identity_downgrade_applied is downgrade


def test_legacy_sumstats_without_alleles_drop_rsids_shared_by_panel_rows() -> None:
    panel = _panel(
        pd.DataFrame(
            {
                "CHR": ["1", "1", "1"], "POS": [10, 10, 20], "SNP": ["multi", "multi", "good"],
                "A1": ["A", "A", "C"], "A2": ["G", "T", "T"],
                "regression_ld_scores": [1.0, 1.0, 1.0], "base": [1.0, 2.0, 3.0],
            }
        ),
        "rsid_allele_aware",
    )
    projected, drops = _project_legacy_sumstats_to_panel(
        _legacy(pd.DataFrame({"SNP": ["multi", "good"], "Z": [1.0, 2.0], "N": [1000.0, 1000.0]})), panel,
    )

    assert projected.data["SNP"].tolist() == ["good"]
    assert projected.has_alleles is False
    assert projected.config_snapshot.snp_identifier == "rsid"
    assert drops["reason"].tolist() == ["ambiguous_panel_mapping"]
    assert drops[["A1", "A2"]].isna().all().all()


def test_legacy_sumstats_with_one_allele_column_are_rejected() -> None:
    panel = _panel(
        pd.DataFrame({"CHR": ["1"], "POS": [10], "SNP": ["rs1"], "regression_ld_scores": [1.0], "base": [1.0]}),
        "rsid",
    )
    sumstats = _legacy(pd.DataFrame({"SNP": ["rs1"], "A1": ["A"], "Z": [1.0], "N": [1000.0]}))

    with pytest.raises(LDSCInputError, match="has A1 but not A2"):
        RegressionRunner(regression_config=RegressionConfig()).build_dataset(sumstats, panel)


def test_h2_command_runs_allele_less_legacy_sumstats_and_records_the_downgrade(tmp_path: Path) -> None:
    from ldsc._regression.runner import run_h2

    n_snps = 60
    rng = np.random.default_rng(11)
    snps = [f"rs{index}" for index in range(n_snps)]
    ldscore_dir = tmp_path / "ldscores"
    ldscore_dir.mkdir()
    pd.DataFrame(
        {
            "CHR": ["1"] * n_snps, "SNP": snps, "POS": np.arange(1, n_snps + 1) * 100,
            "A1": ["A"] * n_snps, "A2": ["C"] * n_snps,
            "regression_ld_scores": rng.uniform(1.0, 5.0, n_snps), "base": rng.uniform(1.0, 5.0, n_snps),
        }
    ).to_parquet(ldscore_dir / "ldscore.baseline.parquet", index=False)
    (ldscore_dir / "metadata.json").write_text(json.dumps({
        "schema_version": 1, "artifact_type": "ldscore", "files": {"baseline": "ldscore.baseline.parquet"},
        "snp_identifier": "rsid_allele_aware", "genome_build": "hg38", "chromosomes": ["1"],
        "baseline_columns": ["base"], "query_columns": [], "query_batches": [],
        "baseline_row_groups": [{"chrom": "1", "row_group_index": 0, "row_offset": 0, "n_rows": n_snps}],
        "counts": [{"group": "baseline", "column": "base", "all_reference_snp_count": 100.0,
                    "common_reference_snp_count": 80.0}],
    }), encoding="utf-8")
    sumstats = tmp_path / "trait.sumstats.gz"
    with gzip.open(sumstats, "wt", encoding="utf-8") as handle:
        handle.write("SNP\tZ\tN\n")
        handle.writelines(f"{snp}\t{rng.normal(0.0, 1.5):.4f}\t1000\n" for snp in snps)

    run_h2(sumstats_file=str(sumstats), ldscore_dir=str(ldscore_dir), output_dir=str(tmp_path / "out"), num_blocks=10)

    metadata = json.loads((tmp_path / "out" / "diagnostics" / "metadata.json").read_text())
    assert metadata["effective_snp_identifier"] == "rsid"
    assert metadata["identity_downgrade_applied"] is True


def test_h2_command_records_legacy_projection_drops_in_its_audit(tmp_path: Path) -> None:
    from ldsc._regression.runner import run_h2

    n_snps = 60
    rng = np.random.default_rng(7)
    snps = [f"rs{index}" for index in range(n_snps)]
    ldscore_dir = tmp_path / "ldscores"
    ldscore_dir.mkdir()
    pd.DataFrame(
        {
            "CHR": ["1"] * n_snps, "SNP": snps, "POS": np.arange(1, n_snps + 1) * 100,
            "A1": ["A"] * n_snps, "A2": ["C"] * n_snps,
            "regression_ld_scores": rng.uniform(1.0, 5.0, n_snps), "base": rng.uniform(1.0, 5.0, n_snps),
        }
    ).to_parquet(ldscore_dir / "ldscore.baseline.parquet", index=False)
    (ldscore_dir / "metadata.json").write_text(json.dumps({
        "schema_version": 1, "artifact_type": "ldscore", "files": {"baseline": "ldscore.baseline.parquet"},
        "snp_identifier": "rsid", "genome_build": "hg38", "chromosomes": ["1"], "baseline_columns": ["base"],
        "query_columns": [], "query_batches": [],
        "baseline_row_groups": [{"chrom": "1", "row_group_index": 0, "row_offset": 0, "n_rows": n_snps}],
        "counts": [{"group": "baseline", "column": "base", "all_reference_snp_count": 100.0,
                    "common_reference_snp_count": 80.0}],
    }), encoding="utf-8")
    # rs0 appears twice and "missing" is absent from the LD scores.
    rows = [(snp, rng.normal(0.0, 1.5)) for snp in snps] + [("rs0", 0.5), ("missing", 1.0)]
    sumstats = tmp_path / "trait.sumstats.gz"
    with gzip.open(sumstats, "wt", encoding="utf-8") as handle:
        handle.write("SNP\tA1\tA2\tZ\tN\n")
        handle.writelines(f"{snp}\tA\tC\t{z:.4f}\t1000\n" for snp, z in rows)

    run_h2(sumstats_file=str(sumstats), ldscore_dir=str(ldscore_dir), output_dir=str(tmp_path / "out"), num_blocks=10)

    audit = pd.read_csv(tmp_path / "out/diagnostics/dropped_snps/legacy_sumstats.tsv.gz", sep="\t")
    assert audit["reason"].value_counts().to_dict() == {"duplicate_source_rsid": 2, "missing_panel_rsid": 1}
    assert sorted(audit["SNP"]) == ["missing", "rs0", "rs0"]
