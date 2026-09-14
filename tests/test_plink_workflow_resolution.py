"""PLINK token forms preserve chromosome routing and numerical workflow results."""

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from ldsc import gene_ldscore_index, ldscore_calculator, ref_panel_builder
from ldsc.errors import LDSCInputError


def write_inputs(root):
    """Write separated-cM SNPs whose self-only LD scores equal one."""
    metadata = []
    for chrom in ("21", "22"):
        prefix = root / f"1000G.EUR.QC.{chrom}"
        frame = pd.DataFrame({"CHR": [chrom] * 4, "SNP": [f"rs{chrom}_{i}" for i in range(4)],
                              "CM": np.arange(4) * .1, "POS": np.arange(1, 5) * 100,
                              "A1": "A", "A2": "C"})
        frame.to_csv(str(prefix) + ".bim", sep="\t", header=False, index=False)
        Path(str(prefix) + ".fam").write_text("".join(f"F{i} I{i} 0 0 0 -9\n" for i in range(4)))
        Path(str(prefix) + ".bed").write_bytes(b"\x6c\x1b\x01" + b"\xf0" * 4)
        metadata.append(frame)
    frame = pd.concat(metadata, ignore_index=True)
    frame.assign(base=1.0).to_csv(root / "baseline.annot", sep="\t", index=False)
    frame.assign(query=np.tile([1., 0., 0., 0.], 2)).to_csv(root / "query.annot", sep="\t", index=False)
    frame[["SNP", "CHR", "POS"]].to_csv(root / "regression.tsv", sep="\t", index=False)
    (root / "genes.tsv").write_text("gene_id\tgene_name\tchrom\tstart\tend\tgenome_build\n"
                                   "ENSG21\tG21\t21\t100\t400\thg19\n"
                                   "ENSG22\tG22\t22\t100\t400\thg19\n")


def index_args(root, token, output):
    args = gene_ldscore_index.build_parser().parse_args([
        "--baseline-annot-sources", str(root / "baseline.annot"), "--plink-prefix", token,
        "--gene-coordinate-file", str(root / "genes.tsv"), "--genome-build", "hg19",
        "--snp-identifier", "chr_pos", "--ld-wind-cm", ".05", "--snp-batch-size", "1", "--padding-bp", "0",
        "--gene-exclude-regions", "none", "--regr-snps-exclude-regions", "none",
        "--regr-snps-file", str(root / "regression.tsv"), "--threads", "2", "--output-dir", str(output),
    ])
    args._test_chromosomes = ("21", "22")
    return args


@pytest.mark.parametrize("form", ["plain", "glob", "at", "exact"])
def test_gene_index_and_ldscore_share_plink_resolution_and_scores(tmp_path, form):
    write_inputs(tmp_path)
    if form == "exact":
        prefix = tmp_path / "joint"
        prefix.with_suffix(".bim").write_bytes(b"".join((tmp_path / f"1000G.EUR.QC.{c}.bim").read_bytes() for c in (21, 22)))
        prefix.with_suffix(".fam").write_bytes((tmp_path / "1000G.EUR.QC.21.fam").read_bytes())
        prefix.with_suffix(".bed").write_bytes(b"\x6c\x1b\x01" + b"\xf0" * 8)
        token = str(prefix)
    else:
        token = str(tmp_path / ("1000G.EUR.QC." + {"plain": "", "glob": "*", "at": "@"}[form]))
    output = gene_ldscore_index.run_build_gene_ldscore_index_from_args(index_args(tmp_path, token, tmp_path / "index"))
    index = gene_ldscore_index._load_gene_ldscore_index(output, _allow_partial_for_tests=True)
    assert index.chromosomes == ("21", "22")
    np.testing.assert_array_equal(index.gene_support, [4, 4])
    for chrom in index.chromosomes:
        record = index.load_chromosome(chrom)
        np.testing.assert_allclose(record.baseline_rows["base"], 1, atol=1e-7)
        np.testing.assert_allclose(record.operator.toarray(), np.ones((4, 1)), atol=1e-7)

    # Direct @ declares all 22 chromosomes; use the equivalent glob for
    # this deliberately partial numerical fixture.
    direct_token = token if form != "at" else str(tmp_path / "1000G.EUR.QC.*")
    args = ldscore_calculator.build_parser().parse_args([
        "--baseline-annot-sources", str(tmp_path / "baseline.annot"), "--plink-prefix", direct_token,
        "--query-annot-sources", str(tmp_path / "query.annot"),
        "--snp-identifier", "chr_pos", "--genome-build", "hg19", "--ld-wind-cm", ".05", "--snp-batch-size", "1",
        "--regr-snps-exclude-regions", "none", "--regr-snps-file", str(tmp_path / "regression.tsv"),
        "--output-dir", str(tmp_path / "ldscore"), "--threads", "2",
    ])
    result = ldscore_calculator.run_ldscore_from_args(args)
    np.testing.assert_allclose(result.baseline_table["base"], 1, atol=1e-7)
    assert len(result.baseline_table) == 8
    np.testing.assert_allclose(result.read_queries(["query"])["query"], np.tile([1., 0., 0., 0.], 2), atol=1e-7)


def test_gene_index_reports_every_plink_problem_before_projection(tmp_path):
    write_inputs(tmp_path)
    (tmp_path / "1000G.EUR.QC.21.bed").unlink()
    (tmp_path / "1000G.EUR.QC.22.fam").unlink()
    args = index_args(tmp_path, str(tmp_path / "1000G.EUR.QC."), tmp_path / "index")
    with pytest.raises(LDSCInputError, match="complete .bed/.bim/.fam"):
        gene_ldscore_index.run_build_gene_ldscore_index_from_args(args)
    audit = pd.read_csv(tmp_path / ".index.build-state/plink_input_issues.tsv", sep="\t")
    missing_paths = audit.source.dropna().tolist()
    assert any(path.endswith("21.bed") for path in missing_paths)
    assert any(path.endswith("22.fam") for path in missing_paths)
    assert not (tmp_path / "index").exists()
    assert not list(tmp_path.glob(".index.stage-*"))


def test_r2_builder_accepts_plain_dotted_prefix(tmp_path):
    write_inputs(tmp_path)
    args = ref_panel_builder.build_parser().parse_args([
        "--plink-prefix", str(tmp_path / "1000G.EUR.QC."), "--source-genome-build", "hg19",
        "--snp-identifier", "chr_pos", "--ld-wind-kb", ".05", "--output-dir", str(tmp_path / "r2"),
    ])
    result = ref_panel_builder.run_build_ref_panel_from_args(args)
    assert len(result.output_paths["r2_hg19"]) == 2
    for path in result.output_paths["r2_hg19"]:
        assert Path(path).is_file()
