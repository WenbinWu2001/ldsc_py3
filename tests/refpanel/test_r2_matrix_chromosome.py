"""One R2 matrix chromosome, read the same way by ``ldscore`` and ``query-r2-matrix``.

Fixture: the deterministic 32-SNP chromosome-22 1KG-30x subset (hg38 positions;
see ``tests/fixtures/minimal_external_resources/README.md``), built into an R2
matrix once per module. Both commands open a chromosome's R2 Parquet file and
sidecar through ``R2MatrixChromosome``, so a damaged sidecar fails through
either command with the same reason, named after the command that read it.
"""
import gzip
import shutil
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from ldsc import cli
from ldsc._refpanel.r2_format import r2_matrix_file_paths
from tests.ref_panel_helpers import write_r2_matrix_chromosome

pytest.importorskip("bitarray")

PLINK = Path(__file__).resolve().parents[1] / "fixtures" / "minimal_external_resources" / "plink" / "hm3_chr22_subset"
ANCHORS = {
    "ldscore": "docs/user-guide/troubleshooting.md#ldscore-r2-matrix-is-incompatible",
    "query-r2-matrix": "docs/user-guide/troubleshooting.md#query-r2-matrix-r2-matrix-is-incompatible",
}


@pytest.fixture(scope="module")
def built(tmp_path_factory):
    root = tmp_path_factory.mktemp("r2") / "matrix"
    assert cli.run_cli(["build-r2-matrix", "--plink-prefix", str(PLINK), "--output-dir", str(root),
                        "--snp-identifier", "chr_pos", "--source-genome-build", "hg38", "--ld-wind-snps", "3"]) == 0
    return root / "hg38"


@pytest.fixture
def matrix(tmp_path, built):
    copy = tmp_path / "matrix"
    shutil.copytree(built, copy)
    return copy


def _rewrite_body(path, edit):
    lines = gzip.decompress(path.read_bytes()).decode().splitlines(keepends=True)
    comments = [line for line in lines if line.startswith("#")]
    path.write_bytes(gzip.compress("".join(comments + edit(lines[len(comments):])).encode()))


def rename_header(path):
    _rewrite_body(path, lambda body: [body[0].replace("CHR\tPOS\t", "CHROM\tBP\t", 1)] + body[1:])


def reorder_rows(path):
    _rewrite_body(path, lambda body: [body[0], body[2], body[1]] + body[3:])


def remove(path):
    path.unlink()


def _run(command, matrix, tmp_path, caplog):
    output = tmp_path / "out"
    if command == "ldscore":
        argv = ["ldscore", "--r2-matrix-dir", str(matrix), "--output-dir", str(output), "--snp-identifier", "chr_pos",
                "--genome-build", "hg38", "--ld-wind-snps", "3"]
    else:
        positions = pd.read_csv(PLINK.with_suffix(".bim"), sep=r"\s+", header=None)[3]
        pairs = tmp_path / "pairs.tsv"
        pairs.write_text(f"CHR_1\tPOS_1\tCHR_2\tPOS_2\n22\t{positions[0]}\t22\t{positions[1]}\n")
        argv = ["query-r2-matrix", "--r2-matrix-dir", str(matrix), "--pairs-file", str(pairs),
                "--output-dir", str(output)]
    code = cli.run_cli(argv)
    return code, caplog.text, output


@pytest.mark.parametrize("command", ["ldscore", "query-r2-matrix"])
@pytest.mark.parametrize(("damage", "reason"), [
    (rename_header, "CHROM, BP, SNP"),
    (remove, "is missing"),
    (reorder_rows, "ldsc:sidecar_identity_sha256"),
])
def test_damaged_sidecar_fails_the_same_way_through_both_commands(tmp_path, matrix, caplog, command, damage, reason):
    _, sidecar = r2_matrix_file_paths(matrix, "22")
    damage(sidecar)

    code, message, output = _run(command, matrix, tmp_path, caplog)

    assert code != 0
    assert f"{command} could not use R2 matrix chromosome 22" in message
    assert str(sidecar) in message
    assert reason in message
    assert ANCHORS[command] in message
    assert not (output / "ldscore.baseline.parquet").exists()
    assert not (output / "query_r2_matrix.tsv").exists()


def _tied_panel(build_dir):
    """A 4-SNP rsID panel whose rows 1 and 2 (rs3, then rs2) share position 200."""
    panel = pd.DataFrame({"CHR": ["1"] * 4, "POS": [100, 200, 200, 300], "SNP": ["rs1", "rs3", "rs2", "rs4"],
                          "A1": ["A"] * 4, "A2": ["C", "C", "G", "C"], "MAF": [0.3, 0.25, 0.4, 0.35]})
    pairs = [{"i": i, "j": j, "R2": r2, "sign": "+"} for i, j, r2 in
             [(0, 1, 0.5), (0, 2, 0.1), (1, 2, 0.4), (1, 3, 0.2), (2, 3, 0.6)]]
    write_r2_matrix_chromosome(build_dir, panel, pairs, snp_identifier="rsid", ld_window=("snps", 1.0))
    return panel


def test_preparing_a_chromosome_parses_its_sidecar_once(tmp_path, monkeypatch):
    """Reference filtering, the sidecar binding check and the pair reader share one parse.

    The annotation lists the two same-position SNPs in the other order, so
    preparation also rebuilds the pair reader in R2 matrix row order.
    """
    from ldsc import GlobalConfig, LDScoreConfig, RefPanelConfig
    from ldsc._refpanel import r2_format
    from ldsc._refpanel.adapters import R2MatrixRefPanel
    from ldsc._refpanel.prepared import ChromosomeAnnotations

    panel = _tied_panel(tmp_path / "panel")
    parsed = []
    original = r2_format.read_r2_sidecar
    monkeypatch.setattr(r2_format, "read_r2_sidecar", lambda *args, **kwargs: parsed.append(args) or
                        original(*args, **kwargs))
    ref_panel = R2MatrixRefPanel(GlobalConfig(snp_identifier="rsid"),
                                 RefPanelConfig(backend="r2_matrix", r2_matrix_dir=str(tmp_path / "panel"), maf_min=0.1))
    metadata = panel[["CHR", "POS", "SNP", "A1", "A2"]].iloc[[0, 2, 1, 3]].reset_index(drop=True)
    annotations = ChromosomeAnnotations(metadata, pd.DataFrame({"base": np.ones(4)}), ["base"], [])

    with ref_panel.prepare_chromosome("1", annotations, LDScoreConfig(ld_wind_snps=1)) as prepared:
        assert prepared.reader.panel_rows.tolist() == [0, 1, 2, 3]

    assert len(parsed) == 1


def test_genome_build_checks_find_files_named_by_the_r2_matrix_naming_rule(tmp_path):
    """``_genome`` sits below ``_refpanel`` and spells the R2 matrix file names itself.

    Files named by ``r2_matrix_file_paths`` must be found by its directory scan,
    and the sidecar next to each R2 Parquet file must be read as a declaration.
    """
    from ldsc._genome.build_inspection import inspect_genome_build
    from ldsc._genome.build_policy import GenomeBuildCheck
    from ldsc.errors import LDSCInputError

    panel = pd.DataFrame({"CHR": ["1"] * 2, "POS": [100, 200], "SNP": ["rs1", "rs2"], "A1": ["A", "A"],
                          "A2": ["C", "G"], "MAF": [0.3, 0.2]})
    _, sidecar = write_r2_matrix_chromosome(tmp_path / "hg38", panel, genome_build="hg38")

    assert "hg38" in inspect_genome_build(str(tmp_path))
    check = GenomeBuildCheck("auto", context="test")
    check.r2("--r2-matrix-dir", tmp_path)
    assert check.resolve() == "hg38"

    text = gzip.decompress(sidecar.read_bytes()).decode()
    sidecar.write_bytes(gzip.compress(text.replace("# ldsc:genome_build=hg38", "# ldsc:genome_build=hg19").encode()))
    check = GenomeBuildCheck("auto", context="test")
    check.r2("--r2-matrix-dir", tmp_path / "hg38")
    with pytest.raises(LDSCInputError, match="Conflicting genome-build evidence") as failure:
        check.resolve()
    assert sidecar.name in str(failure.value)
    with pytest.raises(LDSCInputError, match=sidecar.name):
        inspect_genome_build(str(tmp_path))
