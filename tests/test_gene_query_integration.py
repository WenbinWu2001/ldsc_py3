"""Whole-workflow gene viability with real PLINK and immutable index artifacts."""

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from ldsc import gene_ldscore_index as index_workflow
from ldsc.ldscore_calculator import build_parser, run_ldscore_from_args
from ldsc.errors import LDSCInputError


@pytest.fixture
def inputs(tmp_path):
    fixture = Path(__file__).parent / "fixtures/plink/plink"
    prefix = tmp_path / "panel"
    for suffix in (".bed", ".fam"):
        prefix.with_suffix(suffix).write_bytes(fixture.with_suffix(suffix).read_bytes())
    bim = pd.read_csv(fixture.with_suffix(".bim"), sep=r"\s+", header=None,
                      names=["CHR", "SNP", "CM", "POS", "A1", "A2"])
    bim["CM"] = np.arange(len(bim)) * 0.2
    bim["A1"], bim["A2"] = "A", "C"
    bim.to_csv(prefix.with_suffix(".bim"), sep="\t", header=False, index=False)
    baseline = tmp_path / "baseline.annot"
    bim.assign(base=1).to_csv(baseline, sep="\t", index=False)
    regression = tmp_path / "regression.tsv"
    bim[["SNP"]].to_csv(regression, sep="\t", index=False)
    catalog = tmp_path / "genes.tsv"
    catalog.write_text("gene_id\tgene_name\tchrom\tstart\tend\tgenome_build\n"
                       "G1\tONE\t1\t5\t6\thg19\n"
                       "G2\tTWO\t1\t6\t8\thg19\n"
                       "G0\tZERO\t1\t100\t110\thg19\n"
                       "GX\tOUTSIDE\t22\t10\t20\thg19\n")
    return prefix, baseline, regression, catalog


def direct_args(tmp_path, inputs, sources, *, policy="strict", control=None, name="direct"):
    prefix, baseline, regression, catalog = inputs
    argv = ["--output-dir", str(tmp_path / name), "--baseline-annot-sources", str(baseline),
            "--query-annot-gene-list-sources", ",".join(map(str, sources)), "--gene-coordinate-file", str(catalog),
            "--padding-bp", "0",
            "--plink-prefix", str(prefix), "--genome-build", "hg19", "--snp-identifier", "rsid",
            "--ld-wind-cm", "0.2", "--regr-snps-file", str(regression), "--regr-snps-exclude-regions", "none"]
    if policy == "resolved-only":
        argv += ["--allow-unresolved-genes"]
    if control:
        argv += ["--control-gene-list-file", str(control)]
    return build_parser().parse_args(argv)


def build_index(tmp_path, inputs):
    prefix, baseline, regression, catalog = inputs
    args = index_workflow.build_parser().parse_args([
        "--output-dir", str(tmp_path / "index"), "--baseline-annot-sources", str(baseline),
        "--plink-prefix", str(prefix), "--gene-coordinate-file", str(catalog),
        "--genome-build", "hg19", "--snp-identifier", "rsid", "--ld-wind-cm", "0.2",
        "--gene-exclude-regions", "none", "--regr-snps-exclude-regions", "none",
        "--regr-snps-file", str(regression),
    ])
    args._test_chromosomes = ("1",)
    return index_workflow.run_build_gene_ldscore_index_from_args(args)


@pytest.mark.parametrize("name", ["base", "SNP"])
def test_direct_gate_a_reports_naming_and_identifier_issues_together(tmp_path, inputs, name):
    named, invalid = tmp_path / f"{name}.txt", tmp_path / "invalid.txt"
    named.write_text("G1\n")
    invalid.write_text("UNKNOWN\n")
    args = direct_args(tmp_path, inputs, [named, invalid])
    with pytest.raises(LDSCInputError):
        run_ldscore_from_args(args)
    summary = pd.read_csv(tmp_path / "direct/diagnostics/gene_list_resolution_summary.tsv", sep="\t").set_index("query")
    assert summary.loc[name, "source_status"] == "error"
    assert "annotation_name_collision" in summary.loc[name, "source_reasons"]
    assert summary.loc["invalid", "rejected_rows"] == 1
    audit = pd.read_csv(tmp_path / "direct/diagnostics/gene_list_audit.tsv.gz", sep="\t")
    assert len(audit) == 2


def test_direct_log_has_preparation_milestones_and_one_exclusion_summary(tmp_path, inputs):
    catalog = inputs[3]
    with catalog.open("a") as stream:
        stream.write("MHC1\tFIRST\t6\t30000000\t30000100\thg19\n"
                     "MHC2\tSECOND\t6\t31000000\t31000100\thg19\n")
    focal = tmp_path / "focal.txt"
    focal.write_text("MHC1\n\nSECOND\nG1\n")
    args = direct_args(tmp_path, inputs, [focal])
    args.gene_exclude_regions = "mhc"
    result = run_ldscore_from_args(args)
    log = (tmp_path / "direct/diagnostics/ldscore.log").read_text()
    assert log.count("Reading annotation inputs:") == 1
    assert log.count("Checking SNP identities and preparing chromosome annotations.") == 1
    assert log.count("Annotation preparation complete:") == 1
    assert log.count("contains CM/MAF") == 1
    assert log.count("intentionally excluded by region policy:") == 1
    assert "source=focal.txt count=2 line:gene=[1:MHC1, 3:SECOND->MHC2]" in log
    assert result.query_columns == ["focal"]
    audit = pd.read_csv(tmp_path / "direct/diagnostics/gene_list_audit.tsv.gz", sep="\t")
    assert audit.disposition.tolist() == ["excluded", "excluded", "retained"]
    assert audit.line.tolist() == [1, 3, 4]


def test_multichromosome_contents_cover_a_single_chromosome_pathway(tmp_path, inputs):
    prefix, baseline, regression, catalog = inputs
    metadata = pd.read_csv(baseline, sep="\t")
    metadata["CHR"] = [21] * 6 + [22] * 2
    metadata[["CHR", "SNP", "CM", "POS", "A1", "A2"]].to_csv(
        prefix.with_suffix(".bim"), sep="\t", header=False, index=False)
    misleading = tmp_path / "baseline.22.annot"
    metadata.to_csv(misleading, sep="\t", index=False)
    catalog.write_text("gene_id\tgene_name\tchrom\tstart\tend\tgenome_build\nG22\tFOCAL\t22\t7\t7\thg19\n")
    focal = tmp_path / "focal.txt"
    focal.write_text("G22\n")
    args = direct_args(tmp_path, (prefix, misleading, regression, catalog), [focal])
    args.baseline_annot_sources = str(tmp_path / "*.22.annot")
    args.yes_really = True
    result = run_ldscore_from_args(args)
    assert result.chromosome_scope["analysis_chromosomes"] == ["21", "22"]
    assert set(result.baseline_table.CHR.astype(str)) == {"21", "22"}
    assert result.gene_list_batch.summary.coverage_status.tolist() == ["full"]
    assert result.query_columns == ["focal"]


def test_plink_glob_routes_by_validated_contents_not_names(tmp_path, inputs):
    from tests.ref_panel_helpers import write_tiny_plink

    _, _, regression, catalog = inputs
    metadata = pd.DataFrame({"CHR": [21, 21, 22, 22], "SNP": ["rs1", "rs2", "rs3", "rs4"],
                             "POS": [1, 2, 1, 2], "base": [1] * 4})
    for chrom, filename in ((21, "swap.22"), (22, "swap.21")):
        write_tiny_plink(tmp_path / (filename + ".panel"), metadata[metadata.CHR.eq(chrom)])
    baseline = tmp_path / "baseline.22.annot"
    metadata.to_csv(baseline, sep="\t", index=False)
    metadata[["SNP"]].to_csv(regression, sep="\t", index=False)
    catalog.write_text("gene_id\tgene_name\tchrom\tstart\tend\tgenome_build\nG22\tFOCAL\t22\t1\t1\thg19\n")
    focal = tmp_path / "focal.txt"
    focal.write_text("G22\n")
    args = direct_args(tmp_path, (tmp_path / "swap.*", baseline, regression, catalog), [focal])
    args.ld_wind_cm, args.ld_wind_snps = None, 1
    args.yes_really = True
    result = run_ldscore_from_args(args)
    assert result.chromosome_scope["analysis_chromosomes"] == ["21", "22"]
    assert result.query_columns == ["focal"]
    assert result.gene_list_batch.summary.genes_with_snp_support.tolist() == [1]


@pytest.mark.parametrize("policy", ["strict", "resolved-only"])
def test_real_direct_and_indexed_outputs_agree(tmp_path, inputs, policy):
    mixed, empty, zero, control = [tmp_path / f"{name}.txt" for name in ("mixed", "empty", "zero", "control")]
    mixed.write_text("G1\nONE\nG0\n" + ("UNKNOWN\n" if policy == "resolved-only" else ""))
    empty.write_text("")
    zero.write_text("G0\n")
    control.write_text("G2\n")
    sources = [mixed, empty, zero]
    direct = run_ldscore_from_args(direct_args(tmp_path, inputs, sources, policy=policy, control=control))
    indexed = index_workflow.run_indexed_ldscore(
        build_index(tmp_path, inputs), query_gene_list_sources=sources, control_gene_list_file=control,
        gene_list_resolution_policy=policy, output_dir=tmp_path / "indexed", _allow_partial_for_tests=True,
    )
    assert direct.query_columns == indexed.query_columns == ["mixed"]
    assert [s.reason for s in direct.query_statuses] == [s.reason for s in indexed.query_statuses] == [
        "partial_gene_resolution" if policy == "resolved-only" else "partial_snp_support", "empty_gene_list", "zero_annotation_snps"]
    np.testing.assert_allclose(direct.read_queries(["mixed"])["mixed"], indexed.read_queries(["mixed"])["mixed"], rtol=1e-6, atol=1e-7)
    np.testing.assert_allclose(direct.baseline_table[["base", "gene_control", "regression_ld_scores"]],
                               indexed.baseline_table[["base", "gene_control", "regression_ld_scores"]], rtol=1e-6, atol=1e-7)
    assert direct.count_records == indexed.count_records
    pd.testing.assert_frame_equal(direct.overlap.baseline_block_all, indexed.overlap.baseline_block_all)
    pd.testing.assert_frame_equal(direct.gene_list_batch.summary, indexed.gene_list_batch.summary)
    for name in ("direct", "indexed"):
        scope = json.loads((tmp_path / name / "diagnostics/chromosome_scope.json").read_text())
        assert scope["chromosomes"] == ["1"]
        metadata = json.loads((tmp_path / name / "metadata.json").read_text())
        assert metadata["query_columns"] == ["mixed"]
        assert all((tmp_path / name / path).exists() for path in metadata["files"].values())
    log = (tmp_path / "direct/diagnostics/ldscore.log").read_text()
    assert "Chromosomes resolved and entering the analysis: 1" in log


@pytest.mark.parametrize("policy", ["strict", "resolved-only"])
def test_incomplete_pathways_and_control_fail_together_with_unknown_support(tmp_path, inputs, policy):
    partial, absent, full, control = [tmp_path / f"{name}.txt" for name in ("partial", "absent", "full", "control")]
    partial.write_text("G1\nGX\n")
    absent.write_text("GX\n")
    full.write_text("G1\n")
    control.write_text("G2\nGX\n")
    with pytest.raises(LDSCInputError, match="coverage preflight failed") as caught:
        run_ldscore_from_args(direct_args(tmp_path, inputs, [partial, absent, full], policy=policy, control=control))
    assert all(name in str(caught.value) for name in ("partial", "absent", "control"))
    log = (tmp_path / "direct/diagnostics/ldscore.log").read_text()
    for message in (str(caught.value), log):
        assert "'@' requires all chromosomes 1-22" in message
        assert "use quoted '*' patterns" in message
        assert "matching baseline and PLINK chromosome sets" in message
        assert "does not filter genes" in message
    summary = pd.read_csv(tmp_path / "direct/diagnostics/gene_list_resolution_summary.tsv", sep="\t")
    assert summary.coverage_status.tolist() == ["partial", "none", "full", "partial"]
    assert summary.zero_support_genes.isna().all()
    assert not (tmp_path / "direct/metadata.json").exists()


@pytest.mark.parametrize("control_text", ["", "G0\n"])
def test_empty_or_zero_support_control_is_fatal(tmp_path, inputs, control_text):
    full, control = tmp_path / "full.txt", tmp_path / "control.txt"
    full.write_text("G1\n")
    control.write_text(control_text)
    with pytest.raises(LDSCInputError, match="control gene list"):
        run_ldscore_from_args(direct_args(tmp_path, inputs, [full], control=control))
    assert (tmp_path / "direct/diagnostics/gene_list_resolution_summary.tsv").exists()
    assert not (tmp_path / "direct/metadata.json").exists()


def test_all_zero_support_queries_write_diagnostics_and_fail(tmp_path, inputs):
    zero = tmp_path / "zero.txt"
    zero.write_text("G0\n")
    with pytest.raises(LDSCInputError, match="all 1 requested query annotations were skipped"):
        run_ldscore_from_args(direct_args(tmp_path, inputs, [zero]))
    summary = pd.read_csv(tmp_path / "direct/diagnostics/gene_list_resolution_summary.tsv", sep="\t")
    assert summary.coverage_status.tolist() == ["full"]
    assert summary.zero_support_genes.tolist() == [1]
    assert not (tmp_path / "direct/metadata.json").exists()


@pytest.mark.parametrize("role", ["focal", "control"])
def test_zero_variance_is_assessed_on_regression_rows(tmp_path, inputs, role):
    prefix, baseline, regression, catalog = inputs
    regression.write_text("SNP\nrs_4\n")
    focal, control = tmp_path / "focal.txt", tmp_path / "control.txt"
    focal.write_text("G1\n")
    control.write_text("G2\n")
    with pytest.raises(LDSCInputError, match="zero.variance|zero_variance"):
        run_ldscore_from_args(direct_args(tmp_path, inputs, [focal], control=control if role == "control" else None))
    status = pd.read_csv(tmp_path / "direct/diagnostics/query_annotation_status.tsv", sep="\t")
    assert status.reason.tolist() == ["zero_variance_ld_scores"]
    assert status.n_annotation_snps.iloc[0] > 0
    assert not (tmp_path / "direct/metadata.json").exists()


def test_support_uses_baseline_intersection_and_keeps_valid_zero_chromosome(tmp_path, inputs):
    prefix, baseline, regression, catalog = inputs
    frame = pd.read_csv(baseline, sep="\t")
    frame.loc[frame.POS.ge(7)].to_csv(baseline, sep="\t", index=False)
    focal = tmp_path / "focal.txt"
    focal.write_text("G1\n")
    with pytest.raises(LDSCInputError, match="all 1 requested query annotations were skipped"):
        args = direct_args(tmp_path, inputs, [focal])
        args.yes_really = True
        run_ldscore_from_args(args)
    summary = pd.read_csv(tmp_path / "direct/diagnostics/gene_list_resolution_summary.tsv", sep="\t")
    assert summary.coverage_status.tolist() == ["full"]
    assert summary.zero_support_genes.tolist() == [1]


@pytest.mark.parametrize("policy", ["strict", "resolved-only"])
def test_indexed_incomplete_coverage_remains_a_separate_fatal_gate(tmp_path, inputs, policy):
    index = build_index(tmp_path, inputs)
    focal, control = tmp_path / "focal.txt", tmp_path / "control.txt"
    focal.write_text("G1\nGX\n")
    control.write_text("GX\n")
    with pytest.raises(LDSCInputError, match="incomplete chromosome coverage"):
        index_workflow.run_indexed_ldscore(index, query_gene_list_sources=[focal], control_gene_list_file=control,
                                         gene_list_resolution_policy=policy, output_dir=tmp_path / "indexed", _allow_partial_for_tests=True)
    summary = pd.read_csv(tmp_path / "indexed/diagnostics/gene_list_resolution_summary.tsv", sep="\t")
    assert summary.coverage_status.tolist() == ["partial", "none"]
    assert summary.rejected_rows.tolist() == [0, 0]
    assert summary.zero_support_genes.isna().all()


def test_indexed_missing_components_are_reported_together(tmp_path, inputs):
    index = build_index(tmp_path, inputs)
    for name in ("gene_to_atom.npz", "atoms.parquet"):
        (index / "chromosomes/chr1" / name).unlink()
    focal = tmp_path / "focal.txt"
    focal.write_text("G1\n")
    with pytest.raises(LDSCInputError, match="preflight"):
        index_workflow.run_indexed_ldscore(index, query_gene_list_sources=[focal], output_dir=tmp_path / "indexed", _allow_partial_for_tests=True)
    issues = pd.read_csv(tmp_path / "indexed/diagnostics/input_issues.tsv", sep="\t")
    assert len(issues) == 2
    assert not (tmp_path / "indexed/diagnostics/gene_list_audit.tsv.gz").exists()
