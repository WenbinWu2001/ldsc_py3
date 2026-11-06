"""Genetic-map errors name the command that read the map and link its troubleshooting section."""
from pathlib import Path

import pandas as pd
import pytest

from ldsc import cli
from ldsc.config import GeneLDScoreIndexBuildConfig
from ldsc.errors import LDSCInputError

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures" / "minimal_external_resources"
TROUBLESHOOTING = Path(__file__).resolve().parents[2] / "docs" / "user-guide" / "troubleshooting.md"
LDSCORE_DOC = "docs/user-guide/troubleshooting.md#ldscore-a-genetic-map-cannot-be-used"
GENE_INDEX_DOC = "docs/user-guide/troubleshooting.md#build-gene-ldscore-index-a-genetic-map-cannot-be-used"


def _duplicated_map(tmp_path):
    """The hg38 fixture map with its second row repeated."""
    rows = pd.read_csv(FIXTURES / "genetic_maps/genetic_map_hg38_chr22_subset.txt", sep=r"\s+")
    path = tmp_path / "dup_map.txt"
    pd.concat([rows.iloc[:2], rows.iloc[1:]], ignore_index=True).to_csv(path, sep=" ", index=False)
    return path


def test_ldscore_names_itself_for_a_bad_genetic_map(tmp_path, caplog):
    bim = pd.read_csv(FIXTURES / "plink/hm3_chr22_subset.bim", sep=r"\s+", header=None,
                      names=["CHR", "SNP", "CM", "POS", "A1", "A2"])
    bim.loc[:, ["SNP", "CHR", "POS"]].to_csv(tmp_path / "regr.tsv", sep="\t", index=False)
    status = cli.run_cli(["ldscore", "--plink-prefix", str(FIXTURES / "plink/hm3_chr22_subset"),
                          "--genome-build", "hg38", "--snp-identifier", "chr_pos_allele_aware", "--ld-wind-cm", "1",
                          "--genetic-map-hg38-sources", str(_duplicated_map(tmp_path)),
                          "--regr-snps-file", str(tmp_path / "regr.tsv"), "--regr-snps-exclude-regions", "none",
                          "--output-dir", str(tmp_path / "out")])

    assert status == 1
    assert f"ldscore could not load genetic map '{tmp_path / 'dup_map.txt'}': duplicate CHR/POS" in caplog.text
    assert f"Other causes & fixes: {LDSCORE_DOC}" in caplog.text
    assert "build-r2-matrix" not in caplog.text


def test_gene_index_build_names_itself_for_a_bad_genetic_map(tmp_path):
    from ldsc._ldscore.gene_index.build import _load_builder_genetic_map

    path = _duplicated_map(tmp_path)
    config = GeneLDScoreIndexBuildConfig(
        baseline_annot_sources=("baseline",), plink_prefix="panel", output_dir="index",
        gene_coordinate_file="genes", genome_build="hg38", snp_identifier="rsid", ld_wind_cm=1.0,
    )
    with pytest.raises(LDSCInputError) as raised:
        _load_builder_genetic_map(config, {"hg19": None, "hg38": str(path)})

    message = str(raised.value)
    assert message.startswith(f"build-gene-ldscore-index could not load genetic map '{path}': duplicate CHR/POS")
    assert message.endswith(f"Other causes & fixes: {GENE_INDEX_DOC}")


@pytest.mark.parametrize("anchor", [LDSCORE_DOC, GENE_INDEX_DOC])
def test_genetic_map_troubleshooting_sections_exist(anchor):
    heading = {
        LDSCORE_DOC: "### `ldscore`: a genetic map cannot be used",
        GENE_INDEX_DOC: "### `build-gene-ldscore-index`: a genetic map cannot be used",
    }[anchor]
    assert heading in TROUBLESHOOTING.read_text(encoding="utf-8").splitlines()
