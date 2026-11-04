"""Reference-panel rows whose alleles cannot form an allele-aware SNP identity.

A row with a missing, non-ACGT (such as an indel) or identical allele pair has
no allele-aware SNP identity. In allele-aware modes every reference-panel
command drops such rows before computing anything, lists each one in
``diagnostics/dropped_snps/`` and prints one WARNING with the counts to stderr.
Base modes ignore alleles and keep them.

The 32-SNP chromosome-22 fixture has only valid pairs, so each test copies it
and relabels two SNPs' ``.bim`` alleles: one indel and one identical pair.
"""
import gzip
from pathlib import Path
import shutil

import pandas as pd
import pytest

from ldsc.cli import run_cli

FIXTURE = Path(__file__).resolve().parents[1] / "fixtures" / "minimal_external_resources" / "plink" / "hm3_chr22_subset"
UNUSABLE = {10: ("AT", "G"), 12: ("G", "G")}


def write_panel(directory: Path) -> Path:
    """Copy the fixture with two unusable allele pairs; return the PLINK prefix."""
    directory.mkdir(parents=True, exist_ok=True)
    prefix = directory / "panel"
    for suffix in (".bed", ".fam"):
        shutil.copy(str(FIXTURE) + suffix, str(prefix) + suffix)
    bim = pd.read_csv(str(FIXTURE) + ".bim", sep=r"\s+", header=None, names=["CHR", "SNP", "CM", "POS", "A1", "A2"])
    for row, alleles in UNUSABLE.items():
        bim.loc[row, ["A1", "A2"]] = alleles
    bim.to_csv(str(prefix) + ".bim", sep="\t", header=False, index=False)
    bim.loc[:, ["SNP", "CHR", "POS"]].to_csv(directory / "regression_snps.tsv", sep="\t", index=False)
    return prefix


def read_drops(path: Path) -> pd.DataFrame:
    with gzip.open(path, "rt") as handle:
        return pd.read_csv(handle, sep="\t")


def ldscore_args(prefix: Path, mode: str, output_dir: Path) -> list[str]:
    return ["ldscore", "--plink-prefix", str(prefix), "--genome-build", "hg38", "--snp-identifier", mode,
            "--ld-wind-kb", "1000", "--regr-snps-file", str(prefix.parent / "regression_snps.tsv"),
            "--regr-snps-exclude-regions", "none", "--output-dir", str(output_dir)]


@pytest.mark.parametrize("mode", ["chr_pos_allele_aware", "rsid_allele_aware"])
def test_plink_ldscore_drops_unusable_alleles_and_reports_them(tmp_path, capsys, mode):
    prefix = write_panel(tmp_path / "panel")

    status = run_cli(ldscore_args(prefix, mode, tmp_path / "out"))

    assert status == 0
    drops = read_drops(tmp_path / "out/diagnostics/dropped_snps/chr22_dropped.tsv.gz")
    assert drops["reason"].tolist() == ["invalid_allele", "invalid_allele"]
    assert set(drops["stage"]) == {"plink_reference_identity_cleanup"}
    funnel = pd.read_csv(tmp_path / "out/diagnostics/reference_snp_funnel.tsv", sep="\t")
    step = funnel.loc[funnel["funnel"].eq("reference_panel") & funnel["step"].eq("identity_cleanup")]
    assert step["removed"].tolist() == [2]
    assert len(pd.read_parquet(tmp_path / "out/ldscore.baseline.parquet")) == 30
    stderr = capsys.readouterr().err
    assert "WARNING: ldscore dropped 2 SNP rows" in stderr
    assert "invalid_allele=2" in stderr
    assert "diagnostics/dropped_snps" in stderr


def test_base_mode_keeps_rows_whose_alleles_are_unusable(tmp_path, capsys):
    prefix = write_panel(tmp_path / "panel")

    status = run_cli(ldscore_args(prefix, "chr_pos", tmp_path / "out"))

    assert status == 0
    assert len(pd.read_parquet(tmp_path / "out/ldscore.baseline.parquet")) == 32
    assert "dropped" not in capsys.readouterr().err


def test_build_r2_matrix_reports_unusable_alleles_on_stderr(tmp_path, capsys):
    prefix = write_panel(tmp_path / "panel")

    status = run_cli(["build-r2-matrix", "--plink-prefix", str(prefix), "--source-genome-build", "hg38",
                      "--snp-identifier", "chr_pos_allele_aware", "--ld-wind-kb", "1000",
                      "--output-dir", str(tmp_path / "r2")])

    assert status == 0
    drops = read_drops(tmp_path / "r2/diagnostics/dropped_snps/chr22_dropped.tsv.gz")
    assert drops["reason"].tolist() == ["invalid_allele", "invalid_allele"]
    stderr = capsys.readouterr().err
    assert "WARNING: build-r2-matrix dropped 2 SNP rows" in stderr
    assert "invalid_allele=2" in stderr


def test_annotate_reports_unusable_alleles_on_stderr(tmp_path, capsys):
    prefix = write_panel(tmp_path / "panel")
    bim = pd.read_csv(str(prefix) + ".bim", sep="\t", header=None)
    (tmp_path / "query.bed").write_text(f"chr22\t{int(bim[3].min()) - 1}\t{int(bim[3].max())}\n")

    status = run_cli(["annotate", "--ref-snp-bim-sources", str(prefix) + ".bim", "--genome-build", "hg38",
                      "--snp-identifier", "chr_pos_allele_aware", "--query-annot-bed-sources",
                      str(tmp_path / "query.bed"), "--padding-bp", "0", "--output-dir", str(tmp_path / "annot")])

    assert status == 0
    drops = read_drops(tmp_path / "annot/diagnostics/dropped_snps/dropped.tsv.gz")
    assert drops["reason"].tolist() == ["invalid_allele", "invalid_allele"]
    stderr = capsys.readouterr().err
    assert "WARNING: annotate dropped 2 SNP rows" in stderr
    assert "invalid_allele=2" in stderr


def test_build_gene_ldscore_index_reports_unusable_alleles_on_stderr(tmp_path, capsys, monkeypatch):
    from ldsc._ldscore.gene_index.build import main as build_index
    from tests.gene_index_helpers import narrow_index_chromosomes

    narrow_index_chromosomes(monkeypatch, ["22"])
    prefix = write_panel(tmp_path / "panel")
    bim = pd.read_csv(str(prefix) + ".bim", sep="\t", header=None, names=["CHR", "SNP", "CM", "POS", "A1", "A2"])
    bim.loc[:, ["CHR", "SNP", "POS", "CM", "A1", "A2"]].assign(base=1.0).to_csv(
        tmp_path / "baseline.annot", sep="\t", index=False)
    (tmp_path / "genes.tsv").write_text(
        "gene_id\tgene_name\tchrom\tstart\tend\tgenome_build\n"
        f"G1\tGene1\t22\t{int(bim.POS.min())}\t{int(bim.POS.median())}\thg38\n")

    build_index(["--baseline-annot-sources", str(tmp_path / "baseline.annot"), "--plink-prefix", str(prefix),
                 "--gene-coordinate-file", str(tmp_path / "genes.tsv"), "--genome-build", "hg38",
                 "--snp-identifier", "chr_pos_allele_aware", "--ld-wind-kb", "1000", "--padding-bp", "0",
                 "--gene-exclude-regions", "none", "--regr-snps-exclude-regions", "none",
                 "--regr-snps-file", str(prefix.parent / "regression_snps.tsv"), "--output-dir", str(tmp_path / "index")])

    stderr = capsys.readouterr().err
    assert "WARNING: build-gene-ldscore-index dropped" in stderr
    assert "invalid_allele=" in stderr
