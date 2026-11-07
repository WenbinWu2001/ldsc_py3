"""Genetic-map errors name the command that read the map and link its troubleshooting section."""
from pathlib import Path

import pandas as pd
import pytest

from ldsc import cli
from ldsc.config import GeneLDScoreIndexBuildConfig
from ldsc.errors import LDSCInputError
from tests.ref_panel_helpers import read_bim

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
    bim = read_bim(FIXTURES / "plink/hm3_chr22_subset")
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


def test_gene_index_build_rejects_a_map_without_an_index_chromosome(tmp_path, monkeypatch):
    from ldsc._ldscore.gene_index import build as index
    from tests.gene_index_helpers import narrow_index_chromosomes
    from tests.ref_panel_helpers import ORIENTATION_PANEL, write_counted_plink

    narrow_index_chromosomes(monkeypatch, ["1"])
    write_counted_plink(tmp_path / "panel", ORIENTATION_PANEL)
    rows = pd.DataFrame({"CHR": "1", "SNP": [row[0] for row in ORIENTATION_PANEL], "POS": [10, 20, 30, 40, 50]})
    rows.assign(base=1).to_csv(tmp_path / "baseline.annot", sep="\t", index=False)
    rows.to_csv(tmp_path / "regression.tsv", sep="\t", index=False)
    pd.DataFrame({"gene_id": ["G1"], "gene_name": ["Gene1"], "chrom": ["1"], "start": [1], "end": [5],
                  "genome_build": ["hg19"]}).to_csv(tmp_path / "genes.tsv", sep="\t", index=False)
    genetic_map = tmp_path / "map_chr2.txt"
    pd.DataFrame({"CHR": ["2", "2"], "POS": [1, 100], "CM": [0.0, 1.0]}).to_csv(genetic_map, sep=" ", index=False)
    config = GeneLDScoreIndexBuildConfig(
        baseline_annot_sources=(str(tmp_path / "baseline.annot"),), plink_prefix=str(tmp_path / "panel"),
        output_dir=str(tmp_path / "index"), gene_coordinate_file=str(tmp_path / "genes.tsv"), genome_build="hg19",
        snp_identifier="rsid", regr_snps_file=str(tmp_path / "regression.tsv"), regr_snps_exclude_regions="none",
        gene_exclude_regions="none", ld_wind_cm=1.0)

    with pytest.raises(LDSCInputError) as raised:
        index._build(config, genetic_map_sources={"hg19": str(genetic_map), "hg38": None},
                     overwrite=False, log_level="INFO")

    message = str(raised.value)
    assert "build-gene-ldscore-index could not use genetic map" in message
    assert "no rows for chromosome 1" in message
    assert GENE_INDEX_DOC in message
    assert "ldscore could not" not in message
    # The input gate stops the build before any chromosome task starts.
    assert not index._chromosome_log_dir(tmp_path / "index" / ".build").exists()
