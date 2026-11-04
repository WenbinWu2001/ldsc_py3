"""``--ref-panel-snps-file`` applies before identity cleanup in every PLINK route.

Identity cleanup runs within the restricted panel: a position is multi-allelic
only if the restriction keeps more than one of its allele pairs. So
``ldscore --plink-prefix`` and ``build-r2-matrix`` followed by
``ldscore --r2-matrix-dir`` keep the same reference SNPs from one panel and
one restriction.

Each test copies the 32-SNP chromosome-22 fixture (one row per position, all
allele pairs valid) and edits its ``.bim``. Expected SNP sets follow from the
edits, independently of the package code.
"""
import gzip
from pathlib import Path
import shutil

import pandas as pd
import pytest

from ldsc.cli import run_cli

FIXTURE = Path(__file__).resolve().parents[1] / "fixtures" / "minimal_external_resources" / "plink" / "hm3_chr22_subset"
BIM_COLUMNS = ["CHR", "SNP", "CM", "POS", "A1", "A2"]


def read_fixture_bim() -> pd.DataFrame:
    return pd.read_csv(str(FIXTURE) + ".bim", sep=r"\s+", header=None, names=BIM_COLUMNS)


def write_panel(directory: Path, bim: pd.DataFrame) -> Path:
    """Copy the fixture genotypes with an edited ``.bim``; return the PLINK prefix."""
    directory.mkdir(parents=True, exist_ok=True)
    prefix = directory / "panel"
    for suffix in (".bed", ".fam"):
        shutil.copy(str(FIXTURE) + suffix, str(prefix) + suffix)
    bim.to_csv(str(prefix) + ".bim", sep="\t", header=False, index=False)
    bim.loc[:, ["SNP", "CHR", "POS"]].to_csv(directory / "regression_snps.tsv", sep="\t", index=False)
    return prefix


def write_restriction(path: Path, rows: pd.DataFrame) -> Path:
    rows.loc[:, ["CHR", "SNP", "POS", "A1", "A2"]].to_csv(path, sep="\t", index=False)
    return path


def read_drops(path: Path) -> pd.DataFrame:
    with gzip.open(path, "rt") as handle:
        return pd.read_csv(handle, sep="\t")


def run_three_routes(tmp_path: Path, prefix: Path, restriction: Path, mode: str) -> dict[str, list[str]]:
    """Run PLINK ``ldscore``, ``build-r2-matrix`` and R2 ``ldscore``; return each route's SNPs."""
    common = ["--snp-identifier", mode, "--ld-wind-kb", "1000"]
    build = [] if mode.startswith("rsid") else ["--genome-build", "hg38"]
    regression = ["--regr-snps-file", str(prefix.parent / "regression_snps.tsv"), "--regr-snps-exclude-regions", "none"]
    assert run_cli(["ldscore", "--plink-prefix", str(prefix), *build, *common,
                    "--ref-panel-snps-file", str(restriction), *regression,
                    "--output-dir", str(tmp_path / "ld_plink")]) == 0
    assert run_cli(["build-r2-matrix", "--plink-prefix", str(prefix), "--source-genome-build", "hg38", *common,
                    "--ref-panel-snps-file", str(restriction), "--output-dir", str(tmp_path / "r2")]) == 0
    assert run_cli(["ldscore", "--r2-matrix-dir", str(tmp_path / "r2"), *build, *common,
                    *regression, "--output-dir", str(tmp_path / "ld_r2")]) == 0
    sidecar = pd.read_csv(tmp_path / "r2" / "hg38" / "chr22_meta.tsv.gz", sep="\t", comment="#")
    return {
        "plink_ldscore": sorted(pd.read_parquet(tmp_path / "ld_plink" / "ldscore.baseline.parquet")["SNP"]),
        "r2_sidecar": sorted(sidecar["SNP"]),
        "r2_ldscore": sorted(pd.read_parquet(tmp_path / "ld_r2" / "ldscore.baseline.parquet")["SNP"]),
    }


def test_restriction_listing_one_allele_pair_of_a_multi_allelic_position_keeps_that_snp(tmp_path):
    bim = read_fixture_bim()
    # Row 1 moves to row 0's position with another allele pair: one position, two allele pairs.
    bim.loc[1, "POS"] = bim.loc[0, "POS"]
    bim.loc[1, ["A1", "A2"]] = ("G", "C")
    assert tuple(bim.loc[0, ["A1", "A2"]]) == ("G", "A")
    prefix = write_panel(tmp_path / "panel", bim)
    # The restriction lists every row except row 1, so the restricted panel has one row per position.
    restriction = write_restriction(tmp_path / "restrict.tsv", bim.drop(index=1))
    expected = sorted(bim.drop(index=1)["SNP"])

    snps = run_three_routes(tmp_path, prefix, restriction, "chr_pos_allele_aware")

    assert len(expected) == 31 and bim.loc[0, "SNP"] in expected
    assert snps == {"plink_ldscore": expected, "r2_sidecar": expected, "r2_ldscore": expected}
    for route in ("ld_plink", "r2"):
        assert read_drops(tmp_path / route / "diagnostics" / "dropped_snps" / "chr22_dropped.tsv.gz").empty
    funnel = pd.read_csv(tmp_path / "ld_plink" / "diagnostics" / "reference_snp_funnel.tsv", sep="\t")
    reference = funnel.loc[funnel["funnel"].eq("reference_panel")].sort_values("step_index")
    assert reference["step"].tolist()[:3] == ["panel_rows", "ref_panel_snps_file", "identity_cleanup"]
    assert reference["removed"].tolist()[1:3] == [1, 0]


@pytest.mark.parametrize("mode", ["chr_pos_allele_aware", "rsid_allele_aware"])
def test_restriction_tolerates_unusable_panel_alleles_that_cleanup_then_drops(tmp_path, mode):
    original = read_fixture_bim()
    bim = original.copy()
    bim.loc[10, ["A1", "A2"]] = ("AT", "G")  # an indel the restriction lists with the fixture's alleles
    bim.loc[12, ["A1", "A2"]] = ("G", "G")  # an identical pair the restriction leaves out
    prefix = write_panel(tmp_path / "panel", bim)
    restriction = write_restriction(tmp_path / "restrict.tsv", original.drop(index=12))
    expected = sorted(original.drop(index=[10, 12])["SNP"])

    snps = run_three_routes(tmp_path, prefix, restriction, mode)

    assert snps == {"plink_ldscore": expected, "r2_sidecar": expected, "r2_ldscore": expected}
    # Only the listed row reaches cleanup, so only it is reported, as build-r2-matrix reports it.
    for route in ("ld_plink", "r2"):
        drops = read_drops(tmp_path / route / "diagnostics" / "dropped_snps" / "chr22_dropped.tsv.gz")
        assert drops["SNP"].tolist() == [bim.loc[10, "SNP"]]
        assert drops["reason"].tolist() == ["invalid_allele"]
