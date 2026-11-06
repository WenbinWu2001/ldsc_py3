"""A/T and C/G SNPs are valid SNP identities, so reference panels keep them in every mode.

The 32-SNP chromosome-22 fixture has no A/T or C/G SNP. Each test copies it
and relabels two SNPs' ``.bim`` alleles, one to A/T and one to C/G. The
genotype bytes are unchanged, so LD among the SNPs, and with it every LD
score and SNP count, must equal the ``chr_pos`` result, which ignores alleles.
"""
import json
from pathlib import Path
import shutil

import numpy as np
import pandas as pd
import pytest

from ldsc import GlobalConfig, reset_global_config, set_global_config
from ldsc._ldscore import workflow
from ldsc._refpanel import builder
from ldsc._refpanel.r2_format import R2_QUANT_SCALE

FIXTURE = Path(__file__).resolve().parents[1] / "fixtures" / "minimal_external_resources" / "plink" / "hm3_chr22_subset"
RELABELLED = {3: ("A", "T"), 7: ("C", "G")}
ALLELE_AWARE = ["chr_pos_allele_aware", "rsid_allele_aware"]
ALL_MODES = ["rsid", "rsid_allele_aware", "chr_pos", "chr_pos_allele_aware"]


def write_relabelled_panel(directory: Path) -> Path:
    """Copy the fixture with two SNPs relabelled A/T and C/G; return the PLINK prefix."""
    directory.mkdir(parents=True, exist_ok=True)
    prefix = directory / "panel"
    for suffix in (".bed", ".fam"):
        shutil.copy(str(FIXTURE) + suffix, str(prefix) + suffix)
    bim = pd.read_csv(str(FIXTURE) + ".bim", sep=r"\s+", header=None, names=["CHR", "SNP", "CM", "POS", "A1", "A2"])
    for row, alleles in RELABELLED.items():
        bim.loc[row, ["A1", "A2"]] = alleles
    bim.to_csv(str(prefix) + ".bim", sep="\t", header=False, index=False)
    bim.loc[:, ["SNP", "CHR", "POS", "A1", "A2"]].to_csv(directory / "regression_snps.tsv", sep="\t", index=False)
    return prefix


def plink_ldscores(tmp_path: Path, prefix: Path, mode: str) -> pd.DataFrame:
    set_global_config(GlobalConfig(snp_identifier=mode, genome_build="hg38"))
    try:
        source = workflow.run_ldscore(
            output_dir=str(tmp_path / f"ldscore_{mode}"), plink_prefix=str(prefix), ld_wind_kb=1000,
            regr_snps_file=str(prefix.parent / "regression_snps.tsv"), regr_snps_exclude_regions="none",
        )
    finally:
        reset_global_config()
    return source.baseline_table, source.count_records


@pytest.mark.parametrize("mode", ALLELE_AWARE)
def test_plink_ldscore_keeps_a_t_and_c_g_snps(tmp_path, mode):
    prefix = write_relabelled_panel(tmp_path / "panel")
    expected, expected_counts = plink_ldscores(tmp_path, prefix, "chr_pos")
    observed, observed_counts = plink_ldscores(tmp_path, prefix, mode)

    assert len(observed) == 32
    assert {frozenset("AT"), frozenset("CG")} <= {frozenset(pair) for pair in zip(observed["A1"], observed["A2"])}
    np.testing.assert_array_equal(observed["POS"].to_numpy(), expected["POS"].to_numpy())
    np.testing.assert_array_equal(observed["base"].to_numpy(), expected["base"].to_numpy())
    assert observed_counts == expected_counts


# All four modes: PLINK ldscore, build-r2-matrix and R2 ldscore share one PLINK selection stage.
@pytest.mark.parametrize("mode", ALL_MODES)
def test_r2_matrix_keeps_a_t_and_c_g_snps_and_matches_plink_ldscores(tmp_path, mode):
    prefix = write_relabelled_panel(tmp_path / "panel")
    expected, _ = plink_ldscores(tmp_path, prefix, mode)
    set_global_config(GlobalConfig(snp_identifier=mode, genome_build="hg38"))
    try:
        built = builder.run_build_r2_matrix(plink_prefix=str(prefix), source_genome_build="hg38",
                                            output_dir=str(tmp_path / "r2"), ld_wind_kb=1000)
        sidecar = pd.read_csv(built.output_paths["meta_hg38"][0], sep="\t", comment="#")
        metadata = json.loads(Path(built.output_paths["metadata"][0]).read_text())
        source = workflow.run_ldscore(
            output_dir=str(tmp_path / "ldscore_r2"), r2_matrix_dir=str(tmp_path / "r2" / "hg38"), ld_wind_kb=1000,
            regr_snps_file=str(prefix.parent / "regression_snps.tsv"), regr_snps_exclude_regions="none",
        )
    finally:
        reset_global_config()

    assert metadata["snp_identifier"] == mode
    assert len(sidecar) == 32
    assert len(source.baseline_table) == 32
    np.testing.assert_array_equal(np.sort(sidecar["POS"].to_numpy()), expected["POS"].to_numpy())
    np.testing.assert_array_equal(source.baseline_table["SNP"].to_numpy(), expected["SNP"].to_numpy())
    # The R2 matrix stores each R2 as int16 / R2_QUANT_SCALE, so each of the
    # at most 31 summed pairs per SNP is off by at most half a quantization step.
    tolerance = 31 * 0.5 / R2_QUANT_SCALE
    np.testing.assert_allclose(source.baseline_table["base"].to_numpy(), expected["base"].to_numpy(), atol=tolerance)
