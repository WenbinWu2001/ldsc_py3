"""Content-based direct input scope and consolidated input failures."""

from argparse import Namespace

import pytest

from ldsc.config import GlobalConfig
from ldsc._ldscore_preflight import inspect_direct_inputs


def write_inputs(tmp_path):
    baseline = tmp_path / "baseline.22.annot"
    baseline.write_text("CHR\tSNP\tPOS\tbase\n22\trs1\t10\t1\n")
    prefix = tmp_path / "panel22"
    prefix.with_suffix(".bim").write_text("22 rs1 0 10 A C\n")
    prefix.with_suffix(".fam").write_text("F I 0 0 0 -9\n")
    prefix.with_suffix(".bed").write_bytes(b"\x6c\x1b\x01\x00")
    return Namespace(baseline_annot_sources=str(baseline), plink_prefix=str(prefix), r2_dir=None,
                     query_annot_sources=None, query_annot_bed_sources=None)


def test_exact_and_glob_scope_comes_from_contents(tmp_path):
    args = write_inputs(tmp_path)
    args.baseline_annot_sources = str(tmp_path / "*.22.annot")
    result = inspect_direct_inputs(args, GlobalConfig(snp_identifier="rsid"))
    assert result.chromosomes == ["22"]
    assert result.issues.empty
    (tmp_path / "baseline.22.annot").write_text("CHR\tSNP\tPOS\tbase\n21\trs1\t10\t1\n")
    result = inspect_direct_inputs(args, GlobalConfig(snp_identifier="rsid"))
    assert result.scope["baseline_chromosomes"] == ["21"]
    assert "chromosome_set_mismatch" in result.issues.reason.tolist()


def test_scope_scan_logs_ignored_metadata_once_across_chunks(tmp_path, caplog):
    args = write_inputs(tmp_path)
    (tmp_path / "baseline.22.annot").write_text("CHR SNP POS MAF base\n" + "22 rs1 10 .2 1\n" * 65537)
    with caplog.at_level("INFO", logger="LDSC"):
        result = inspect_direct_inputs(args, GlobalConfig(snp_identifier="rsid"))
    assert result.issues.empty
    assert sum("contains CM/MAF" in record.getMessage() for record in caplog.records) == 1


def test_at_requires_all_autosomes_and_collects_missing_files(tmp_path):
    args = write_inputs(tmp_path)
    args.baseline_annot_sources = str(tmp_path / "baseline.@.annot")
    result = inspect_direct_inputs(args, GlobalConfig(snp_identifier="rsid"))
    missing = result.issues.loc[result.issues.reason.eq("missing_required_input")]
    assert missing.chrom.tolist() == list(map(str, range(1, 22)))
    assert result.scope["requires_all_autosomes"] is True


def test_multiple_plink_trios_for_one_chromosome_fail_preflight(tmp_path):
    args = write_inputs(tmp_path)
    for suffix in (".bed", ".bim", ".fam"):
        (tmp_path / ("extra22" + suffix)).write_bytes((tmp_path / ("panel22" + suffix)).read_bytes())
    args.plink_prefix = str(tmp_path / "*22")
    result = inspect_direct_inputs(args, GlobalConfig(snp_identifier="rsid"))
    assert "ambiguous_chromosome_input" in result.issues.reason.tolist()


def test_invalid_required_members_are_all_reported(tmp_path):
    args = write_inputs(tmp_path)
    (tmp_path / "panel22.bed").write_bytes(b"bad")
    (tmp_path / "baseline.22.annot").write_text("CHR\tSNP\tPOS\tbase\n22\trs1\t10\tbad\n")
    result = inspect_direct_inputs(args, GlobalConfig(snp_identifier="rsid"))
    assert {"baseline", "reference"} <= set(result.issues.input_role)
    assert result.issues.reason.eq("invalid_required_input").sum() >= 2


def test_truncated_gzip_does_not_hide_other_input_failures(tmp_path):
    import gzip

    args = write_inputs(tmp_path)
    compressed = tmp_path / "baseline.annot.gz"
    compressed.write_bytes(gzip.compress(b"CHR\tSNP\tPOS\tbase\n22\trs1\t10\t1\n")[:-5])
    args.baseline_annot_sources = str(compressed)
    (tmp_path / "panel22.bed").write_bytes(b"bad")
    result = inspect_direct_inputs(args, GlobalConfig(snp_identifier="rsid"))
    assert {"baseline", "reference"} <= set(result.issues.input_role)
    assert result.issues.reason.eq("invalid_required_input").sum() >= 2


def write_r2_inputs(tmp_path):
    import numpy as np
    import pandas as pd
    from ldsc._kernel.ref_panel_builder import write_runtime_metadata_sidecar, write_r2_parquet
    from ldsc._kernel.snp_identity import sidecar_identity_sha256

    args = write_inputs(tmp_path)
    directory = tmp_path / "r2"
    directory.mkdir()
    metadata = pd.DataFrame({"CHR": ["22", "22", "22"], "SNP": ["rs1", "rs2", "rs3"],
                             "POS": [10, 20, 30], "A1": ["A"] * 3, "A2": ["C"] * 3,
                             "CM": [0.0, 0.1, 0.2], "MAF": [0.2, 0.3, 0.4]})
    write_runtime_metadata_sidecar(metadata, directory / "chr22_meta.tsv.gz", genome_build="hg19", snp_identifier="rsid")
    write_r2_parquet(pair_chunks=[(np.array([0, 1]), np.array([1, 2]), np.array([0.5, 0.25]), np.array([1, 1]))],
                     path=directory / "chr22_r2.parquet", genome_build="hg19", n_samples=100,
                     snp_identifier="rsid", min_r2=0.0, n_snps=3, sidecar_identity_sha256=sidecar_identity_sha256(metadata))
    metadata.assign(base=1).to_csv(args.baseline_annot_sources, sep="\t", index=False)
    args.plink_prefix, args.r2_dir = None, str(directory)
    return args


def test_parquet_scope_validates_the_pair_and_sidecar_contents(tmp_path):
    args = write_r2_inputs(tmp_path)
    result = inspect_direct_inputs(args, GlobalConfig(snp_identifier="rsid"))
    assert result.issues.empty
    assert result.chromosomes == ["22"]
    (tmp_path / "r2/chr22_r2.parquet").write_bytes(b"broken parquet")
    result = inspect_direct_inputs(args, GlobalConfig(snp_identifier="rsid"))
    assert result.issues.input_role.eq("reference").any()
    assert result.chromosomes == []


def test_parquet_orphan_members_are_batched(tmp_path):
    args = write_r2_inputs(tmp_path)
    (tmp_path / "r2/chr21_r2.parquet").write_bytes(b"bad")
    (tmp_path / "r2/chr20_meta.tsv.gz").write_bytes(b"bad")
    result = inspect_direct_inputs(args, GlobalConfig(snp_identifier="rsid"))
    assert set(result.issues.chrom) == {"20", "21"}


@pytest.mark.parametrize("maf_min, expected", [(None, "success"), (0.5, "zero_annotation_snps")])
def test_parquet_gene_workflow_distinguishes_valid_zero_support(tmp_path, maf_min, expected):
    import json
    import pandas as pd
    from ldsc.ldscore_calculator import build_parser, run_ldscore_from_args
    from ldsc.errors import LDSCInputError

    inputs = write_r2_inputs(tmp_path)
    genes, catalog, regression = tmp_path / "focal.txt", tmp_path / "genes.tsv", tmp_path / "regression.tsv"
    genes.write_text("G1\n")
    catalog.write_text("gene_id\tgene_name\tchrom\tstart\tend\tgenome_build\nG1\tONE\t22\t10\t10\thg19\n")
    regression.write_text("SNP\nrs1\nrs2\nrs3\n")
    argv = ["--output-dir", str(tmp_path / "out"), "--baseline-annot-sources", inputs.baseline_annot_sources,
            "--r2-dir", inputs.r2_dir, "--query-annot-gene-list-sources", str(genes),
            "--gene-coordinate-file", str(catalog), "--genome-build", "hg19", "--snp-identifier", "rsid",
            "--padding-bp", "0", "--ld-wind-snps", "2", "--yes-really", "--regr-snps-file", str(regression),
            "--regr-snps-exclude-regions", "none"]
    if maf_min is not None:
        argv.extend(["--maf-min", str(maf_min)])
    if expected == "success":
        result = run_ldscore_from_args(build_parser().parse_args(argv))
        assert result.query_columns == ["focal"]
        assert result.gene_list_batch.summary.genes_with_snp_support.tolist() == [1]
    else:
        with pytest.raises(LDSCInputError, match="all 1 requested query annotations were skipped"):
            run_ldscore_from_args(build_parser().parse_args(argv))
        statuses = pd.read_csv(tmp_path / "out/diagnostics/query_annotation_status.tsv", sep="\t")
        assert statuses.reason.tolist() == [expected]
    scope = json.loads((tmp_path / "out/diagnostics/chromosome_scope.json").read_text())
    assert scope["chromosomes"] == ["22"]
