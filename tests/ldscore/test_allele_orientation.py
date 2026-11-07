"""Allele orientation of files derived from PLINK reference-panel genotypes.

Fixture: the deterministic 32-SNP chromosome-22 1KG-30x subset (3202 samples,
hg38 positions; see ``tests/fixtures/minimal_external_resources/README.md``).
Its ``.bim`` ``A1`` is the major allele for 10 SNPs, and no SNP has allele
frequency exactly 0.5. Expected alleles come from the decoded BED dosages.
"""

import numpy as np
import pandas as pd
import pytest

from ldsc import cli
from tests.refpanel.test_r2_matrix_ld_window_rules import _dosages
from tests.ref_panel_helpers import HM3_CHR22_PLINK, read_bim

pytest.importorskip("bitarray")

PLINK = HM3_CHR22_PLINK
BIM = read_bim(PLINK)


@pytest.fixture(scope="module")
def r2_matrix(tmp_path_factory):
    output = tmp_path_factory.mktemp("r2") / "panel"
    assert cli.run_cli(["build-r2-matrix", "--plink-prefix", str(PLINK), "--output-dir", str(output),
                        "--snp-identifier", "chr_pos", "--source-genome-build", "hg38", "--ld-wind-snps", "3"]) == 0
    return output


def _ldscore(output_dir, *args):
    return cli.run_cli(["ldscore", "--output-dir", str(output_dir), "--snp-identifier", "chr_pos",
                        "--genome-build", "hg38", "--ld-wind-snps", "3", *args])


def _minor_alleles():
    """The minor allele as ``A1``, the other as ``A2``, and the ``A1`` frequency, by ``POS``."""
    a1_frequency = _dosages(PLINK).mean(axis=1) / 2
    swap = a1_frequency > 0.5
    assert swap.sum() == 10 and not (a1_frequency == 0.5).any()
    return pd.DataFrame({"A1": np.where(swap, BIM["A2"], BIM["A1"]), "A2": np.where(swap, BIM["A1"], BIM["A2"]),
                         "MAF": np.where(swap, 1 - a1_frequency, a1_frequency)}, index=BIM["POS"].to_numpy())


def test_both_sidecar_writers_and_plink_ld_scores_give_the_minor_allele_as_a1(tmp_path, r2_matrix):
    """``A1`` is the less frequent allele in the decoded dosages and ``MAF`` is its frequency."""
    expected = _minor_alleles()
    output = tmp_path / "out"
    assert _ldscore(output, "--plink-prefix", str(PLINK), "--export-ref-metadata") == 0

    for path in (r2_matrix / "hg38" / "chr22_meta.tsv.gz", output / "ref_metadata" / "chr22_meta.tsv.gz"):
        sidecar = pd.read_csv(path, sep="\t", comment="#")
        assert len(sidecar) == 32, path
        want = expected.loc[sidecar["POS"]]
        assert sidecar["A1"].tolist() == want["A1"].tolist(), path
        assert sidecar["A2"].tolist() == want["A2"].tolist(), path
        np.testing.assert_allclose(sidecar["MAF"].to_numpy(), want["MAF"].to_numpy(), rtol=5e-9, atol=0, err_msg=str(path))
    ldscores = pd.read_parquet(output / "ldscore.baseline.parquet", columns=["POS", "A1", "A2"])
    assert len(ldscores) > 0
    assert ldscores["A1"].tolist() == expected.loc[ldscores["POS"], "A1"].tolist()
    assert ldscores["A2"].tolist() == expected.loc[ldscores["POS"], "A2"].tolist()


@pytest.mark.parametrize("panel", ["r2_matrix", "plink"])
def test_ld_scores_take_alleles_from_the_panel_whatever_the_annotation_alleles(tmp_path, r2_matrix, panel):
    """Annotation alleles only identify SNPs: ``.bim`` order, swapped or absent give one table."""
    expected = _minor_alleles()
    reference = ["--r2-matrix-dir", str(r2_matrix)] if panel == "r2_matrix" else ["--plink-prefix", str(PLINK)]
    annotation = BIM[["CHR", "POS", "SNP"]].assign(base=1)
    variants = {"bim_order": annotation.assign(A1=BIM["A1"], A2=BIM["A2"]),
                "swapped": annotation.assign(A1=BIM["A2"], A2=BIM["A1"]), "no_alleles": annotation}
    tables = {}
    for name, frame in variants.items():
        path = tmp_path / f"{name}.annot"
        frame.to_csv(path, sep="\t", index=False)
        assert _ldscore(tmp_path / name, *reference, "--baseline-annot-sources", str(path)) == 0
        tables[name] = pd.read_parquet(tmp_path / name / "ldscore.baseline.parquet")

    for name, table in tables.items():
        assert {"A1", "A2"} <= set(table.columns), name
        assert table["A1"].tolist() == expected.loc[table["POS"], "A1"].tolist(), name
        assert table["A2"].tolist() == expected.loc[table["POS"], "A2"].tolist(), name
        pd.testing.assert_frame_equal(table, tables["bim_order"], obj=name)


@pytest.mark.parametrize(("mode", "alleles", "logged"), [
    ("chr_pos_allele_aware", False, True), ("chr_pos_allele_aware", True, False), ("chr_pos", False, False)])
def test_ldscore_log_names_annotation_files_matched_by_base_identity(tmp_path, mode, alleles, logged):
    annotation = BIM[["CHR", "POS", "SNP"]].assign(base=1)
    if alleles:
        annotation = annotation.assign(A1=BIM["A1"], A2=BIM["A2"])
    path = tmp_path / "baseline.annot"
    annotation.to_csv(path, sep="\t", index=False)

    output = tmp_path / "out"
    assert _ldscore(output, "--plink-prefix", str(PLINK), "--baseline-annot-sources", str(path),
                    "--snp-identifier", mode) == 0

    notice = ("No annotation file has A1/A2 columns, so annotation SNPs are matched by chr_pos, the base mode of "
              f"chr_pos_allele_aware, without checking alleles: {path}.")
    assert (notice in (output / "diagnostics" / "ldscore.log").read_text()) == logged


def test_ldscore_log_names_a_regression_snp_file_matched_by_base_identity(tmp_path):
    path = tmp_path / "regression.tsv"
    BIM[["CHR", "POS", "SNP"]].to_csv(path, sep="\t", index=False)

    output = tmp_path / "out"
    assert _ldscore(output, "--plink-prefix", str(PLINK), "--regr-snps-file", str(path),
                    "--snp-identifier", "chr_pos_allele_aware") == 0

    notice = (f"SNP restriction file '{path}' has no A1/A2 columns, so it selects SNPs by chr_pos, the base mode of "
              "chr_pos_allele_aware, without checking alleles.")
    assert notice in (output / "diagnostics" / "ldscore.log").read_text()
