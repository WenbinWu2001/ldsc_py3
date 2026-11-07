"""``--ref-panel-snps-file`` applies before identity cleanup in every PLINK route.

Identity cleanup runs within the restricted panel: a position is multi-allelic
only if the restriction keeps more than one of its allele pairs. So, without
liftover or duplicate positions, ``ldscore --plink-prefix`` and
``build-r2-matrix`` followed by ``ldscore --r2-matrix-dir`` keep the same
reference SNPs from one panel and one restriction. A ``.bim`` with a
non-integral position stops every route before computation.

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
from tests.gene_index_helpers import narrow_index_chromosomes
from tests.ref_panel_helpers import copy_panel, HM3_CHR22_PLINK, read_bim

FIXTURE = HM3_CHR22_PLINK


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
    bim = read_bim(FIXTURE)
    # Row 1 moves to row 0's position with another allele pair: one position, two allele pairs.
    bim.loc[1, "POS"] = bim.loc[0, "POS"]
    bim.loc[1, ["A1", "A2"]] = ("G", "C")
    assert tuple(bim.loc[0, ["A1", "A2"]]) == ("G", "A")
    prefix = copy_panel(tmp_path / "panel", bim)
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
    original = read_bim(FIXTURE)
    bim = original.copy()
    bim.loc[10, ["A1", "A2"]] = ("AT", "G")  # an indel the restriction lists with the fixture's alleles
    bim.loc[12, ["A1", "A2"]] = ("G", "G")  # an identical pair the restriction leaves out
    prefix = copy_panel(tmp_path / "panel", bim)
    restriction = write_restriction(tmp_path / "restrict.tsv", original.drop(index=12))
    expected = sorted(original.drop(index=[10, 12])["SNP"])

    snps = run_three_routes(tmp_path, prefix, restriction, mode)

    assert snps == {"plink_ldscore": expected, "r2_sidecar": expected, "r2_ldscore": expected}
    # Only the listed row reaches cleanup, so only it is reported, as build-r2-matrix reports it.
    for route in ("ld_plink", "r2"):
        drops = read_drops(tmp_path / route / "diagnostics" / "dropped_snps" / "chr22_dropped.tsv.gz")
        assert drops["SNP"].tolist() == [bim.loc[10, "SNP"]]
        assert drops["reason"].tolist() == ["invalid_allele"]


@pytest.mark.parametrize("command", ["ldscore", "build-r2-matrix", "build-gene-ldscore-index"])
def test_non_integral_bim_position_stops_every_plink_route_before_computation(tmp_path, caplog, monkeypatch, command):
    bim = read_bim(FIXTURE).astype({"POS": float})
    bim.loc[5, "POS"] += 0.5  # SNP row 6
    prefix = copy_panel(tmp_path / "panel", bim)
    output = tmp_path / "out"
    window = ["--snp-identifier", "chr_pos_allele_aware", "--ld-wind-kb", "1000", "--output-dir", str(output)]
    regression = ["--regr-snps-file", str(prefix.parent / "regression_snps.tsv"), "--regr-snps-exclude-regions", "none"]
    if command == "ldscore":
        argv = ["ldscore", "--plink-prefix", str(prefix), "--genome-build", "hg38", *window, *regression]
    elif command == "build-r2-matrix":
        argv = ["build-r2-matrix", "--plink-prefix", str(prefix), "--source-genome-build", "hg38", *window]
    else:
        # The fixture is chromosome 22 only; an index normally needs autosomes 1 to 22.
        narrow_index_chromosomes(monkeypatch, ("22",))
        baseline = tmp_path / "baseline.annot"
        bim.loc[:, ["CHR", "SNP", "POS", "A1", "A2"]].assign(base=1.0).to_csv(baseline, sep="\t", index=False)
        genes = tmp_path / "genes.tsv"
        genes.write_text("gene_id\tgene_name\tchrom\tstart\tend\tgenome_build\n"
                         "ENSG22\tG22\t22\t16000000\t17000000\thg38\n")
        argv = ["build-gene-ldscore-index", "--plink-prefix", str(prefix), "--baseline-annot-sources", str(baseline),
                "--gene-coordinate-file", str(genes), "--genome-build", "hg38", "--padding-bp", "0",
                "--gene-exclude-regions", "none", *window, *regression]

    status = run_cli(argv)

    # Truncating the position or dropping the row would change the reference SNPs, so the run stops.
    assert status == 1
    assert f"SNP row 6 (SNP {bim.loc[5, 'SNP']}, BP {bim.loc[5, 'POS']})" in caplog.text
    assert not list(output.rglob("*.parquet"))
    assert not (output / "metadata.json").exists()


def write_two_chromosome_panel(tmp_path: Path) -> tuple[Path, pd.DataFrame, pd.DataFrame]:
    """Write a chromosome-21 copy of the fixture and a chromosome 22 whose rows 0 and 1 share a position."""
    original = read_bim(FIXTURE)
    panel = tmp_path / "panel"
    panel.mkdir()
    chr21 = original.assign(CHR=21, SNP="chr21_" + original["SNP"])
    chr22 = original.copy()
    chr22.loc[1, "POS"] = chr22.loc[0, "POS"]  # one position listed twice: a duplicate CHR:POS identity
    for chrom, bim in ((21, chr21), (22, chr22)):
        for suffix in (".bed", ".fam"):
            shutil.copy(str(FIXTURE) + suffix, panel / f"panel.{chrom}{suffix}")
        bim.to_csv(panel / f"panel.{chrom}.bim", sep="\t", header=False, index=False)
        annotation = bim.loc[:, ["CHR", "SNP", "POS"]].drop_duplicates("POS").assign(base=1.0)
        annotation.to_csv(tmp_path / f"baseline.{chrom}.annot.gz", sep="\t", index=False)
    pd.concat([chr21, chr22]).loc[:, ["SNP", "CHR", "POS"]].to_csv(tmp_path / "regr.tsv", sep="\t", index=False)
    return panel, chr21, chr22


def run_restricted_ldscore(tmp_path: Path, reference: list[str], keep: pd.DataFrame, *, annotated: bool) -> int:
    """Run chr_pos ``ldscore`` with ``keep`` as the restriction, with or without baseline annotations."""
    keep.loc[:, ["CHR", "POS"]].to_csv(tmp_path / "keep.tsv", sep="\t", index=False)
    annotations = ["--baseline-annot-sources", str(tmp_path / "baseline.*.annot.gz")] if annotated else []
    return run_cli(["ldscore", *annotations, *reference, "--genome-build", "hg38", "--snp-identifier", "chr_pos",
                    "--ld-wind-kb", "1000", "--ref-panel-snps-file", str(tmp_path / "keep.tsv"),
                    "--regr-snps-file", str(tmp_path / "regr.tsv"), "--regr-snps-exclude-regions", "none",
                    "--output-dir", str(tmp_path / "out"), "--workers", "1"])


def assert_chromosome_22_skipped(output: Path) -> dict[str, list[int]]:
    """Check that only chromosome 21 was computed; return chromosome 22's reference-panel funnel."""
    baseline = pd.read_parquet(output / "ldscore.baseline.parquet")
    assert set(baseline["CHR"].astype(str)) == {"21"}
    assert "Skipping chromosome 22" in (output / "diagnostics" / "ldscore.log").read_text()
    funnel = pd.read_csv(output / "diagnostics" / "reference_snp_funnel.tsv", sep="\t")
    rows = funnel.loc[funnel["funnel"].eq("reference_panel") & funnel["chrom"].astype(str).eq("22")]
    return {row.step: [row.removed, row.remaining] for row in rows.itertuples()}


# Without annotation files, ldscore builds the synthetic `base` annotation from the panel itself.
@pytest.mark.parametrize("annotated", [True, False], ids=["baseline", "synthetic_base"])
def test_chromosome_whose_restricted_rows_cleanup_removes_is_skipped(tmp_path, annotated):
    panel, chr21, chr22 = write_two_chromosome_panel(tmp_path)
    # On chromosome 22 the restriction keeps only the two rows that share a position.
    keep = pd.concat([chr21, chr22.loc[[0, 1]]])

    status = run_restricted_ldscore(tmp_path, ["--plink-prefix", str(panel / "panel.*")], keep, annotated=annotated)

    assert status == 0
    funnel = assert_chromosome_22_skipped(tmp_path / "out")
    assert "duplicate_identity=2" in (tmp_path / "out" / "diagnostics" / "ldscore.log").read_text()
    assert funnel["identity_cleanup"] == [2, 0]


@pytest.mark.parametrize("annotated", [True, False], ids=["baseline", "synthetic_base"])
def test_chromosome_the_restriction_empties_is_skipped(tmp_path, annotated):
    panel, chr21, _chr22 = write_two_chromosome_panel(tmp_path)

    status = run_restricted_ldscore(tmp_path, ["--plink-prefix", str(panel / "panel.*")], chr21, annotated=annotated)

    assert status == 0
    funnel = assert_chromosome_22_skipped(tmp_path / "out")
    assert funnel["ref_panel_snps_file"] == [32, 0]


@pytest.mark.parametrize("annotated", [True, False], ids=["baseline", "synthetic_base"])
def test_r2_matrix_chromosome_the_restriction_empties_is_skipped(tmp_path, annotated):
    panel, chr21, _chr22 = write_two_chromosome_panel(tmp_path)
    assert run_cli(["build-r2-matrix", "--plink-prefix", str(panel / "panel.*"), "--source-genome-build", "hg38",
                    "--snp-identifier", "chr_pos", "--ld-wind-kb", "1000", "--output-dir", str(tmp_path / "r2")]) == 0
    sidecar_rows = len(pd.read_csv(tmp_path / "r2" / "hg38" / "chr22_meta.tsv.gz", sep="\t", comment="#"))

    status = run_restricted_ldscore(tmp_path, ["--r2-matrix-dir", str(tmp_path / "r2")], chr21, annotated=annotated)

    assert status == 0
    funnel = assert_chromosome_22_skipped(tmp_path / "out")
    assert funnel["ref_panel_snps_file"] == [sidecar_rows, 0]
