"""Behavioral tests for explicit LDSC2 LD-score suite conversion."""

from __future__ import annotations

import gzip
import json
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest

import ldsc
from ldsc import cli
from ldsc.errors import LDSCInputError
from ldsc.config import GlobalConfig, RegressionConfig
from ldsc.regression_runner import RegressionRunner
from ldsc.sumstats_munger import SumstatsTable


def _write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.suffix == ".gz":
        with gzip.open(path, "wt", encoding="utf-8") as handle:
            handle.write(text)
    else:
        path.write_text(text, encoding="utf-8")


def _write_unpartitioned_suite(reference: Path, weights: Path, *, missing_m_chrom: int | None = None) -> None:
    for chrom in range(1, 23):
        _write_text(
            reference / f"{chrom}.l2.ldscore.gz",
            f"CHR SNP BP L2\n{chrom} rs{chrom} {chrom * 100} {1.0 + chrom / 100}\n",
        )
        _write_text(
            weights / f"weights.{chrom}.l2.ldscore.gz",
            f"CHR SNP BP L2\n{chrom} rs{chrom} {chrom * 100 + 1} {2.0 + chrom / 100}\n",
        )
        _write_text(reference / f"{chrom}.l2.M_5_50", "10\n")
        if chrom != missing_m_chrom:
            _write_text(reference / f"{chrom}.l2.M", "20\n")


def _write_partitioned_suite(reference: Path, weights: Path, frequencies: Path) -> None:
    for chrom in range(1, 23):
        _write_text(
            reference / f"baseline.{chrom}.l2.ldscore.gz",
            (
                "CHR SNP BP baseL2 annot\n"
                f"{chrom} rs{chrom}a {chrom * 100} 1.0 2.0\n"
                f"{chrom} rs{chrom}b {chrom * 100 + 10} 1.5 2.5\n"
            ),
        )
        _write_text(
            reference / f"baseline.{chrom}.annot.gz",
            (
                "CHR BP SNP CM base annot\n"
                f"{chrom} {chrom * 100} rs{chrom}a 0.0 1 0\n"
                f"{chrom} {chrom * 100 + 10} rs{chrom}b 0.1 1 1\n"
            ),
        )
        _write_text(reference / f"baseline.{chrom}.l2.M", "2 1\n")
        _write_text(reference / f"baseline.{chrom}.l2.M_5_50", "1 0\n")
        _write_text(
            weights / f"weights.{chrom}.l2.ldscore.gz",
            (
                "CHR SNP BP L2\n"
                f"{chrom} rs{chrom}a {chrom * 100} 3.0\n"
                f"{chrom} rs{chrom}b {chrom * 100 + 10} 4.0\n"
            ),
        )
        _write_text(
            frequencies / f"1000G.EUR.QC.{chrom}.frq.gz",
            f"CHR SNP A1 A2 MAF NCHROBS\n{chrom} rs{chrom}a A C 0.10 1000\n{chrom} rs{chrom}b A G 0.01 1000\n",
        )


def test_converter_is_public_and_cli_has_the_approved_minimal_arguments() -> None:
    assert callable(ldsc.convert_ldsc2_ldscores)
    assert ldsc.LegacyLDScoreConverter is not None

    parser = cli.build_parser()
    args = parser.parse_args(
        [
            "convert-ldsc2-ldscores",
            "--legacy-reference-dir",
            "reference",
            "--legacy-weight-dir",
            "weights",
            "--output-dir",
            "converted",
        ]
    )

    assert args.snp_identifier == "rsid"
    assert args.genome_build == "auto"
    assert not hasattr(args, "common_maf_min")


def test_private_kernel_has_no_legacy_ldscore_emitter_or_regression_readers() -> None:
    from ldsc._kernel import formats, ldscore

    for name in (
        "emit_outputs",
        "write_ldscore_file",
        "write_counts",
        "write_annotation_groups",
        "run_ldscore_from_args",
        "run_ldscore",
        "main",
    ):
        assert not hasattr(ldscore, name)
    for name in ("sumstats", "ldscore", "ldscore_fromlist", "M", "M_fromlist", "annot"):
        assert not hasattr(formats, name)


def test_unpartitioned_conversion_writes_reloadable_canonical_suite_with_nullable_all_count(
    tmp_path: Path,
) -> None:
    reference = tmp_path / "eur_w_ld_chr"
    weights = tmp_path / "weights_hm3_no_hla"
    output = tmp_path / "converted"
    _write_unpartitioned_suite(reference, weights, missing_m_chrom=22)
    with gzip.open(reference / "1.l2.ldscore.gz", "rt", encoding="utf-8") as handle:
        _write_text(reference / "1.l2.ldscore", handle.read())

    result = ldsc.convert_ldsc2_ldscores(
        legacy_reference_dir=reference,
        legacy_weight_dir=weights,
        output_dir=output,
    )

    assert result.profile == "unpartitioned"
    assert result.n_rows == 22
    baseline = pd.read_parquet(output / "ldscore.baseline.parquet")
    assert baseline.columns.tolist() == ["CHR", "POS", "SNP", "regression_ld_scores", "base"]
    assert baseline["SNP"].tolist() == [f"rs{chrom}" for chrom in range(1, 23)]
    metadata = json.loads((output / "metadata.json").read_text(encoding="utf-8"))
    assert metadata["count_config"] == {
        "common_reference_snp_maf_min": 0.05,
        "common_reference_snp_maf_operator": ">",
        "common_reference_snp_semantics": "legacy_ldsc2",
    }
    assert metadata["counts"] == [
        {
            "all_reference_snp_count": None,
            "column": "base",
            "common_reference_snp_count": 220.0,
            "group": "baseline",
        }
    ]
    assert metadata["legacy_ldsc2_import"]["profile"] == "unpartitioned"
    assert metadata["legacy_ldsc2_import"]["selected_prefixes"] == {
        "reference": "",
        "weight": "weights.",
    }
    assert metadata["legacy_ldsc2_import"]["ignored_files"] == [str(reference / "1.l2.ldscore")]
    assert (output / "diagnostics" / "conversion_issues.tsv.gz").exists()
    assert (output / "diagnostics" / "convert-ldsc2-ldscores.log").exists()

    loaded = ldsc.load_ldscore_from_dir(str(output))
    assert loaded.baseline_columns == ["base"]
    assert loaded.count_records[0]["all_reference_snp_count"] is None
    issues = pd.read_csv(output / "diagnostics" / "conversion_issues.tsv.gz", sep="\t")
    assert "discarded_identical_plain_duplicate" in set(issues["reason"])


def test_unpartitioned_conversion_inner_joins_by_rsid_and_audits_coordinate_disagreement(tmp_path: Path) -> None:
    reference = tmp_path / "reference"
    weights = tmp_path / "weights"
    output = tmp_path / "converted"
    _write_unpartitioned_suite(reference, weights)
    # Replace chromosome 1 weight row with one shared and one weight-only rsID.
    _write_text(
        weights / "weights.1.l2.ldscore.gz",
        "CHR SNP BP L2\n1 rs1 999 3.0\n1 weight_only 1000 4.0\n",
    )

    ldsc.convert_ldsc2_ldscores(
        legacy_reference_dir=reference,
        legacy_weight_dir=weights,
        output_dir=output,
    )

    baseline = pd.read_parquet(output / "ldscore.baseline.parquet")
    assert baseline.loc[baseline["SNP"] == "rs1", "POS"].item() == 100
    issues = pd.read_csv(output / "diagnostics" / "conversion_issues.tsv.gz", sep="\t")
    assert {"coordinate_disagreement", "weight_only_rsid"}.issubset(set(issues["reason"]))


def test_converter_rejects_missing_required_common_count_shard(tmp_path: Path) -> None:
    reference = tmp_path / "reference"
    weights = tmp_path / "weights"
    _write_unpartitioned_suite(reference, weights)
    (reference / "22.l2.M_5_50").unlink()

    with pytest.raises(LDSCInputError, match=r"M_5_50.*chromosome 22"):
        ldsc.convert_ldsc2_ldscores(
            legacy_reference_dir=reference,
            legacy_weight_dir=weights,
            output_dir=tmp_path / "converted",
        )


def test_chr_pos_conversion_uses_reference_build_inference_and_normalizes_zero_based_positions(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    reference = tmp_path / "reference"
    weights = tmp_path / "weights"
    output = tmp_path / "converted"
    _write_unpartitioned_suite(reference, weights)
    monkeypatch.setattr(
        "ldsc.legacy_ldscore_converter.infer_chr_pos_build",
        lambda *_args, **_kwargs: SimpleNamespace(genome_build="hg19", coordinate_basis="0-based"),
    )

    ldsc.convert_ldsc2_ldscores(
        legacy_reference_dir=reference,
        legacy_weight_dir=weights,
        output_dir=output,
        snp_identifier="chr_pos",
        genome_build="auto",
    )

    baseline = pd.read_parquet(output / "ldscore.baseline.parquet")
    assert baseline.loc[baseline["SNP"] == "rs1", "POS"].item() == 101
    metadata = json.loads((output / "metadata.json").read_text(encoding="utf-8"))
    assert metadata["snp_identifier"] == "chr_pos"
    assert metadata["genome_build"] == "hg19"
    assert metadata["legacy_ldsc2_import"]["coordinate_provenance"] == {
        "coordinate_basis": "0-based",
        "effective_build": "hg19",
        "inference_error": None,
        "inferred_build": "hg19",
        "requested_build": "auto",
    }


def test_baseline_partitioned_conversion_reconstructs_counts_and_full_overlap(tmp_path: Path) -> None:
    reference = tmp_path / "baseline_v1.2"
    weights = tmp_path / "weights_hm3_no_hla"
    frequencies = tmp_path / "1000G_Phase3_frq"
    output = tmp_path / "converted"
    _write_partitioned_suite(reference, weights, frequencies)

    result = ldsc.convert_ldsc2_ldscores(
        legacy_reference_dir=reference,
        legacy_weight_dir=weights,
        legacy_frequency_dir=frequencies,
        output_dir=output,
    )

    assert result.profile == "baseline_partitioned"
    assert result.baseline_columns == ("baseL2", "annot")
    baseline = pd.read_parquet(output / "ldscore.baseline.parquet")
    assert baseline.columns.tolist() == [
        "CHR",
        "POS",
        "SNP",
        "regression_ld_scores",
        "baseL2",
        "annot",
    ]
    metadata = json.loads((output / "metadata.json").read_text(encoding="utf-8"))
    assert metadata["counts"] == [
        {
            "all_reference_snp_count": 44.0,
            "column": "baseL2",
            "common_reference_snp_count": 22.0,
            "group": "baseline",
        },
        {
            "all_reference_snp_count": 22.0,
            "column": "annot",
            "common_reference_snp_count": 0.0,
            "group": "baseline",
        },
    ]
    overlap = pd.read_parquet(output / "ldscore.overlap.parquet")
    all_values = overlap.pivot(index="row_annotation", columns="col_annotation", values="overlap_all_snps")
    common_values = overlap.pivot(index="row_annotation", columns="col_annotation", values="overlap_common_snps")
    assert all_values.loc["baseL2", "baseL2"] == 44.0
    assert all_values.loc["baseL2", "annot"] == 22.0
    assert all_values.loc["annot", "baseL2"] == 22.0
    assert all_values.loc["annot", "annot"] == 22.0
    assert common_values.loc["baseL2", "baseL2"] == 22.0
    assert common_values.loc["annot", "annot"] == 0.0
    assert metadata["overlap_config"]["common_maf_operator"] == ">"
    assert metadata["overlap_config"]["stored_block"] == "baseline_by_baseline"

    loaded = ldsc.load_ldscore_from_dir(str(output))
    assert loaded.overlap is not None
    assert loaded.query_columns == []


def test_baseline_partitioned_conversion_fails_on_conflicting_legacy_count(tmp_path: Path) -> None:
    reference = tmp_path / "baseline"
    weights = tmp_path / "weights"
    frequencies = tmp_path / "frequencies"
    output = tmp_path / "converted"
    _write_partitioned_suite(reference, weights, frequencies)
    _write_text(reference / "baseline.7.l2.M_5_50", "999 0\n")

    with pytest.raises(LDSCInputError, match=r"count conflict.*chromosome 7.*baseL2"):
        ldsc.convert_ldsc2_ldscores(
            legacy_reference_dir=reference,
            legacy_weight_dir=weights,
            legacy_frequency_dir=frequencies,
            output_dir=output,
        )

    assert not (output / "metadata.json").exists()
    issues = pd.read_csv(output / "diagnostics" / "conversion_issues.tsv.gz", sep="\t")
    mismatch = issues.loc[issues["reason"] == "count_mismatch"]
    assert mismatch[["chromosome", "annotation"]].iloc[0].tolist() == [7, "baseL2"]


def test_imported_unpartitioned_suite_rejects_all_count_regression_when_m_is_missing(tmp_path: Path) -> None:
    reference = tmp_path / "reference"
    weights = tmp_path / "weights"
    output = tmp_path / "converted"
    _write_unpartitioned_suite(reference, weights, missing_m_chrom=22)
    ldsc.convert_ldsc2_ldscores(
        legacy_reference_dir=reference,
        legacy_weight_dir=weights,
        output_dir=output,
    )
    panel = ldsc.load_ldscore_from_dir(str(output))
    sumstats = SumstatsTable(
        data=pd.DataFrame(
            {
                "SNP": [f"rs{chrom}" for chrom in range(1, 23)],
                "Z": [1.0] * 22,
                "N": [1000.0] * 22,
            }
        ),
        has_alleles=False,
        source_path="current.parquet",
        trait_name="trait",
        config_snapshot=GlobalConfig(snp_identifier="rsid"),
    )

    with pytest.raises(LDSCInputError, match=r"all-SNP counts are unavailable.*base.*\.M"):
        RegressionRunner(regression_config=RegressionConfig(use_common_counts=False)).build_dataset(sumstats, panel)


def test_loader_rejects_imported_count_and_overlap_threshold_drift(tmp_path: Path) -> None:
    reference = tmp_path / "baseline"
    weights = tmp_path / "weights"
    frequencies = tmp_path / "frequencies"
    output = tmp_path / "converted"
    _write_partitioned_suite(reference, weights, frequencies)
    ldsc.convert_ldsc2_ldscores(
        legacy_reference_dir=reference,
        legacy_weight_dir=weights,
        legacy_frequency_dir=frequencies,
        output_dir=output,
    )
    metadata_path = output / "metadata.json"
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    metadata["overlap_config"]["common_maf_operator"] = ">="
    metadata_path.write_text(json.dumps(metadata), encoding="utf-8")

    with pytest.raises(LDSCInputError, match="count_config and overlap_config disagree"):
        ldsc.load_ldscore_from_dir(str(output))
