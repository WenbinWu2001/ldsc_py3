"""Content-based direct input scope and consolidated input failures."""

from argparse import Namespace
from pathlib import Path

import pandas as pd
import pytest

from ldsc import cli
from ldsc.config import GlobalConfig
from ldsc.errors import LDSCInputError
from ldsc._ldscore.preflight import inspect_direct_inputs, inspect_reference_inputs


def inspect(args, config, tmp_path, *, bed_paths=()):
    """Run the scope gate on inputs prepared as direct ``ldscore`` prepares them."""
    from ldsc._annotation.interval_queries import prepare_bed_queries
    from ldsc._annotation.preflight import resolve_annotation_inputs
    from ldsc._annotation.sources import prepare_annotation_sources
    from ldsc._annotation.storage import AnnotationWorkspace

    reference = inspect_reference_inputs(args, config)
    paths, declared, issues = resolve_annotation_inputs(args.baseline_annot_sources)
    with AnnotationWorkspace(tmp_path) as workspace:
        bed_sources = prepare_bed_queries(list(bed_paths), workspace)
        try:
            prepared = prepare_annotation_sources(workspace, paths, [], mode=config.snp_identifier,
                                                  declared_chromosomes=declared, input_issues=issues,
                                                  autosomes_only=True)
        except LDSCInputError as exc:
            return inspect_direct_inputs(args, config, annotation_issues=exc.input_issues.to_dict("records"),
                                         bed_sources=bed_sources, reference_inputs=reference)
        return inspect_direct_inputs(args, config, annotation_sources=prepared, bed_sources=bed_sources,
                                     reference_inputs=reference)


def write_inputs(tmp_path):
    baseline = tmp_path / "baseline.22.annot"
    baseline.write_text("CHR\tSNP\tPOS\tbase\n22\trs1\t10\t1\n")
    prefix = tmp_path / "panel22"
    prefix.with_suffix(".bim").write_text("22 rs1 0 10 A C\n")
    prefix.with_suffix(".fam").write_text("F I 0 0 0 -9\n")
    prefix.with_suffix(".bed").write_bytes(b"\x6c\x1b\x01\x00")
    return Namespace(baseline_annot_sources=str(baseline), plink_prefix=str(prefix), r2_matrix_dir=None,
                     query_annot_sources=None, query_annot_bed_sources=None)


def write_beds(directory):
    """Two BED files with lines outside autosomes 1-22 and one clean file."""
    directory.mkdir()
    lines = {"a.bed": [("chr22", 5), ("chrX", 1), ("X", 50)],
             "b.bed": [("22", 5), ("chrY", 1), ("chrM", 9), ("MT", 3)],
             "c.bed": [("chr22", 5)]}
    for name, rows in lines.items():
        (directory / name).write_text("".join(f"{chrom}\t{start}\t{start + 10}\n" for chrom, start in rows))
    return str(directory / "*.bed")


def non_autosomal_issues(issues):
    flagged = issues.loc[issues.reason.eq("non_autosomal_chromosomes")]
    return {Path(row.source).name: row for row in flagged.itertuples(index=False)}


def check_non_autosomal_issues(issues):
    flagged = non_autosomal_issues(issues)
    assert sorted(flagged) == ["a.bed", "b.bed"]
    assert flagged["a.bed"].chrom == "X"
    assert "2 lines" in flagged["a.bed"].details and "X (2)" in flagged["a.bed"].details
    assert flagged["b.bed"].chrom == "Y,MT,M"
    assert "3 lines" in flagged["b.bed"].details and "Y (1), MT (1), M (1)" in flagged["b.bed"].details
    for row in flagged.values():
        assert "remove" in row.details.lower() and "remove" in row.repair.lower()
    assert not issues.reason.eq("incomplete_chromosome_coverage").any()


def test_prepared_bed_sources_report_lines_outside_autosomes(tmp_path):
    args = write_inputs(tmp_path)
    write_beds(tmp_path / "beds")
    result = inspect(args, GlobalConfig(snp_identifier="chr_pos", genome_build="hg19"), tmp_path,
                     bed_paths=sorted(Path(tmp_path / "beds").glob("*.bed")))
    check_non_autosomal_issues(result.issues)


def test_ldscore_stops_before_computation_on_bed_lines_outside_autosomes(tmp_path):
    args = write_inputs(tmp_path)
    beds = write_beds(tmp_path / "beds")
    output = tmp_path / "out"
    with pytest.raises(LDSCInputError) as caught:
        cli.main(["ldscore", "--plink-prefix", args.plink_prefix, "--baseline-annot-sources", args.baseline_annot_sources,
                  "--query-annot-bed-sources", beds, "--snp-identifier", "chr_pos", "--genome-build", "hg19",
                  "--ld-wind-snps", "1", "--output-dir", str(output)])
    message = str(caught.value)
    for expected in ("a.bed", "b.bed", "X (2)", "Y (1), MT (1), M (1)", "outside autosomes 1-22", "Remove these lines"):
        assert expected in message
    issues = pd.read_csv(output / "diagnostics" / "input_issues.tsv", sep="\t")
    assert sorted(non_autosomal_issues(issues)) == ["a.bed", "b.bed"]
    assert not (output / "metadata.json").exists()


def test_exact_and_glob_scope_comes_from_contents(tmp_path):
    args = write_inputs(tmp_path)
    args.baseline_annot_sources = str(tmp_path / "*.22.annot")
    result = inspect(args, GlobalConfig(snp_identifier="rsid"), tmp_path)
    assert result.chromosomes == ["22"]
    assert result.issues.empty
    (tmp_path / "baseline.22.annot").write_text("CHR\tSNP\tPOS\tbase\n21\trs1\t10\t1\n")
    result = inspect(args, GlobalConfig(snp_identifier="rsid"), tmp_path)
    assert result.scope["baseline_chromosomes"] == ["21"]
    assert "chromosome_set_mismatch" in result.issues.reason.tolist()


def test_at_requires_all_autosomes_and_collects_missing_files(tmp_path):
    from ldsc._annotation.preflight import resolve_annotation_inputs

    write_inputs(tmp_path)
    _paths, _declared, issues = resolve_annotation_inputs(str(tmp_path / "baseline.@.annot"))
    missing = [row for row in issues if row["reason"] == "missing_required_input"]
    assert [row["chrom"] for row in missing] == list(map(str, range(1, 22)))


def test_multiple_plink_trios_for_one_chromosome_fail_preflight(tmp_path):
    args = write_inputs(tmp_path)
    for suffix in (".bed", ".bim", ".fam"):
        (tmp_path / ("extra22" + suffix)).write_bytes((tmp_path / ("panel22" + suffix)).read_bytes())
    args.plink_prefix = str(tmp_path / "*22")
    result = inspect(args, GlobalConfig(snp_identifier="rsid"), tmp_path)
    assert "ambiguous_chromosome_input" in result.issues.reason.tolist()


def test_invalid_required_members_are_all_reported(tmp_path):
    args = write_inputs(tmp_path)
    (tmp_path / "panel22.bed").write_bytes(b"bad")
    (tmp_path / "baseline.22.annot").write_text("CHR\tSNP\tPOS\tbase\n22\trs1\t10\tbad\n")
    result = inspect(args, GlobalConfig(snp_identifier="rsid"), tmp_path)
    assert {"baseline", "reference"} <= set(result.issues.input_role)
    assert result.issues.reason.eq("invalid_required_input").sum() >= 2


def test_truncated_gzip_does_not_hide_other_input_failures(tmp_path):
    import gzip

    args = write_inputs(tmp_path)
    compressed = tmp_path / "baseline.annot.gz"
    compressed.write_bytes(gzip.compress(b"CHR\tSNP\tPOS\tbase\n22\trs1\t10\t1\n")[:-5])
    args.baseline_annot_sources = str(compressed)
    (tmp_path / "panel22.bed").write_bytes(b"bad")
    result = inspect(args, GlobalConfig(snp_identifier="rsid"), tmp_path)
    assert {"baseline", "reference"} <= set(result.issues.input_role)
    assert result.issues.reason.eq("invalid_required_input").sum() >= 2


def write_r2_inputs(tmp_path):
    import pandas as pd
    from tests.ref_panel_helpers import write_r2_matrix_chromosome

    args = write_inputs(tmp_path)
    directory = tmp_path / "r2"
    metadata = pd.DataFrame({"CHR": ["22", "22", "22"], "SNP": ["rs1", "rs2", "rs3"],
                             "POS": [10, 20, 30], "A1": ["A"] * 3, "A2": ["C"] * 3,
                             "CM": [0.0, 0.1, 0.2], "MAF": [0.2, 0.3, 0.4]})
    write_r2_matrix_chromosome(directory, metadata, [{"i": 0, "j": 1, "R2": 0.5, "sign": "+"},
                                                     {"i": 1, "j": 2, "R2": 0.25, "sign": "+"}],
                               snp_identifier="rsid", ld_window=("snps", 10.0))
    metadata.assign(base=1).to_csv(args.baseline_annot_sources, sep="\t", index=False)
    args.plink_prefix, args.r2_matrix_dir = None, str(directory)
    return args


def test_parquet_scope_validates_the_pair_and_sidecar_contents(tmp_path):
    args = write_r2_inputs(tmp_path)
    result = inspect(args, GlobalConfig(snp_identifier="rsid"), tmp_path)
    assert result.issues.empty
    assert result.chromosomes == ["22"]
    (tmp_path / "r2/chr22_r2.parquet").write_bytes(b"broken parquet")
    result = inspect(args, GlobalConfig(snp_identifier="rsid"), tmp_path)
    assert result.issues.input_role.eq("reference").any()
    assert result.chromosomes == []


def test_r2_matrix_built_in_another_mode_gets_an_r2_matrix_repair(tmp_path):
    args = write_r2_inputs(tmp_path)
    issues = inspect_reference_inputs(args, GlobalConfig(snp_identifier="rsid_allele_aware")).issues

    assert issues.reason.tolist() == ["invalid_required_input"]
    assert "snp_identifier mismatch" in issues.details[0]
    assert "--snp-identifier" in issues.repair[0] and "build-r2-matrix" in issues.repair[0]
    assert "'@'" not in issues.repair[0]


def test_parquet_orphan_members_are_batched(tmp_path):
    args = write_r2_inputs(tmp_path)
    (tmp_path / "r2/chr21_r2.parquet").write_bytes(b"bad")
    (tmp_path / "r2/chr20_meta.tsv.gz").write_bytes(b"bad")
    result = inspect(args, GlobalConfig(snp_identifier="rsid"), tmp_path)
    assert set(result.issues.chrom) == {"20", "21"}


@pytest.mark.parametrize("maf_min, expected", [(None, "success"), (0.5, "zero_annotation_snps")])
def test_parquet_gene_workflow_distinguishes_valid_zero_support(tmp_path, maf_min, expected):
    import json
    import pandas as pd
    from ldsc._ldscore.calculator import build_parser, run_ldscore_from_args
    from ldsc.errors import LDSCInputError

    inputs = write_r2_inputs(tmp_path)
    genes, catalog, regression = tmp_path / "focal.txt", tmp_path / "genes.tsv", tmp_path / "regression.tsv"
    genes.write_text("G1\n")
    catalog.write_text("gene_id\tgene_name\tchrom\tstart\tend\tgenome_build\nG1\tONE\t22\t10\t10\thg19\n")
    regression.write_text("SNP\nrs1\nrs2\nrs3\n")
    argv = ["--output-dir", str(tmp_path / "out"), "--baseline-annot-sources", inputs.baseline_annot_sources,
            "--r2-matrix-dir", inputs.r2_matrix_dir, "--query-annot-gene-list-sources", str(genes),
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
