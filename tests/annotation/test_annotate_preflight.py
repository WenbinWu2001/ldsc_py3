"""Standalone gene-list failures must report independent preflight issues."""

import subprocess
import sys

import pandas as pd
import pytest

from ldsc import GlobalConfig, run_annotate
from ldsc.errors import LDSCInputError
from tests.annotation.test_annotate_streaming import inputs


def grid_inputs(tmp_path, kind):
    baseline, _, _, catalog = inputs(tmp_path, sharded=True)
    if kind == "annot":
        return {"ref_snp_annot_sources": baseline}, catalog
    metadata = pd.concat([pd.read_csv(path, sep="\t") for path in baseline])
    bim = tmp_path / "grid.bim"
    metadata.assign(CM=0, A1="A", A2="C")[
        ["CHR", "SNP", "CM", "POS", "A1", "A2"]
    ].to_csv(bim, sep="\t", header=False, index=False)
    return {"ref_snp_bim_sources": [bim]}, catalog


def assert_no_projection_outputs(output):
    assert not list(output.glob("*.annot.gz"))
    assert not list(output.glob(".ldsc-annotation-*"))
    for name in ("query_annotation_status.tsv", "chromosome_scope.json", "metadata.json"):
        assert not (output / "diagnostics" / name).exists()


def test_annotate_reports_name_and_identifier_issues_together(tmp_path):
    baseline, _, _, catalog = inputs(tmp_path)
    genes = tmp_path / "base.txt"
    genes.write_text("G1\nUNKNOWN\n")
    output = tmp_path / "out"

    with pytest.raises(LDSCInputError) as failure:
        run_annotate(
            ref_snp_annot_sources=baseline,
            query_annot_gene_list_sources=[genes],
            gene_coordinate_file=catalog,
            padding_bp=0,
            genome_build="hg19",
            global_config=GlobalConfig(snp_identifier="rsid"),
            output_dir=output,
        )

    message = str(failure.value)
    assert "UNKNOWN" in message and "unmatched_identifier" in message
    assert "base.txt" in message and "annotation_name_collision" in message
    assert "Rename" in message
    summary = pd.read_csv(output / "diagnostics/gene_list_resolution_summary.tsv", sep="\t")
    assert summary.source_status.tolist() == ["error"]
    assert summary.source_reasons.tolist() == ["annotation_name_collision"]
    assert summary.rejected_rows.tolist() == [1]
    audit = pd.read_csv(output / "diagnostics/gene_list_audit.tsv.gz", sep="\t")
    assert audit.input_gene.tolist() == ["G1", "UNKNOWN"]
    assert audit.disposition.tolist() == ["retained", "rejected"]
    assert audit.reason.iloc[1] == "unmatched_identifier"
    assert_no_projection_outputs(output)


def test_strict_rejection_message_counts_sources_and_guides_repair(tmp_path):
    baseline, _, _, catalog = inputs(tmp_path)
    many, few = tmp_path / "many.txt", tmp_path / "few.txt"
    many.write_text("".join(f"MISSING{index:02d}\n" for index in range(12)))
    few.write_text("G1\nUNKNOWN\n")
    output = tmp_path / "out"

    with pytest.raises(LDSCInputError) as failure:
        run_annotate(
            ref_snp_annot_sources=baseline,
            query_annot_gene_list_sources=[many, few],
            gene_coordinate_file=catalog,
            padding_bp=0,
            genome_build="hg19",
            global_config=GlobalConfig(snp_identifier="rsid"),
            output_dir=output,
        )

    message = str(failure.value)
    assert "rejected 13 of 14 submitted row(s)" in message
    assert "many.txt: 12/12 rejected; few.txt: 1/2 rejected" in message
    assert "3 more rejected row(s) are in the audit" in message
    assert message.index("3 more rejected row(s)") < message.index("has zero usable genes")
    assert "Suggested audit filter: source == 'many.txt' & disposition == 'rejected'" in message
    assert "--overwrite" in message
    assert "docs/user-guide/troubleshooting.md#annotate-gene-list-preflight" in message


def test_strict_rejection_examples_name_their_source(tmp_path):
    baseline, _, _, catalog = inputs(tmp_path)
    first, second = tmp_path / "first.txt", tmp_path / "second.txt"
    first.write_text("G1\nBAD1\n")
    second.write_text("G2\nBAD2\n")

    with pytest.raises(LDSCInputError) as failure:
        run_annotate(
            ref_snp_annot_sources=baseline,
            query_annot_gene_list_sources=[first, second],
            gene_coordinate_file=catalog,
            padding_bp=0,
            genome_build="hg19",
            global_config=GlobalConfig(snp_identifier="rsid"),
            output_dir=tmp_path / "out",
        )

    message = str(failure.value)
    assert "first.txt line 2, 'BAD1', unmatched_identifier" in message
    assert "second.txt line 2, 'BAD2', unmatched_identifier" in message


@pytest.mark.parametrize("kind,name,rows", [
    ("annot", "focal", "UNKNOWN\n"),
    ("annot", "base", "G1\n"),
    ("annot", "base", "UNKNOWN\n"),
    ("annot", "gene_control", "G1\nUNKNOWN\n"),
    ("bim", "focal", "UNKNOWN\n"),
    ("bim", "SNP", "G1\n"),
    ("bim", "SNP", "G1\nUNKNOWN\n"),
])
def test_cli_preflight_message_and_artifacts(tmp_path, kind, name, rows):
    grid, catalog = grid_inputs(tmp_path, kind)
    genes = tmp_path / f"{name}.txt"
    genes.write_text(rows)
    output = tmp_path / "out"
    result = subprocess.run([
        sys.executable, "-m", "ldsc", "annotate",
        f"--ref-snp-{kind}-sources", *map(str, next(iter(grid.values()))),
        "--query-annot-gene-list-sources", str(genes),
        "--gene-coordinate-file", str(catalog), "--padding-bp", "0",
        "--genome-build", "hg19", "--snp-identifier", "rsid",
        "--output-dir", str(output),
    ], capture_output=True, text=True)

    assert result.returncode != 0, result.stdout + result.stderr
    assert genes.name in result.stderr
    rejected = int("UNKNOWN" in rows)
    assert ("unmatched_identifier" in result.stderr) == bool(rejected)
    collision = name != "focal"
    assert ("annotation_name_collision" in result.stderr) == collision
    if collision:
        assert "Rename conflicting query source files" in result.stderr
    for filename in ("gene_list_audit.tsv.gz", "gene_list_resolution_summary.tsv"):
        assert f"diagnostics/{filename}" in result.stderr
    summary = pd.read_csv(output / "diagnostics/gene_list_resolution_summary.tsv", sep="\t").fillna("")
    assert summary.source_status.tolist() == ["error" if collision else "ok"]
    assert summary.source_reasons.tolist() == ["annotation_name_collision" if collision else ""]
    assert summary.rejected_rows.tolist() == [rejected]
    audit = pd.read_csv(output / "diagnostics/gene_list_audit.tsv.gz", sep="\t")
    assert audit.input_gene.tolist() == rows.splitlines()
    assert audit.disposition.eq("rejected").sum() == rejected
    assert "annotation_snp_count" not in audit
    assert_no_projection_outputs(output)


@pytest.mark.parametrize("kind", ["annot", "bim"])
@pytest.mark.parametrize("issue", ["identifier", "naming", "both"])
def test_gene_gate_stops_before_grid_preparation(tmp_path, monkeypatch, kind, issue):
    from ldsc._annotation import annotate as annotate_workflow

    grid, catalog = grid_inputs(tmp_path, kind)
    name = "focal" if issue == "identifier" else ("base" if kind == "annot" else "SNP")
    genes = tmp_path / f"{name}.txt"
    genes.write_text("G1\n" if issue == "naming" else "G1\nUNKNOWN\n")

    def unexpected_preparation(*args, **kwargs):
        pytest.fail("Fatal source/name preflight must stop before build inference or grid preparation")

    for operation in ("_projection_build", "prepare_annotation_sources", "prepare_bim_grid", "build_query_shards"):
        monkeypatch.setattr(annotate_workflow, operation, unexpected_preparation)
    output = tmp_path / "out"
    with pytest.raises(LDSCInputError, match="Gene-list catalog preflight") as failure:
        run_annotate(**grid, query_annot_gene_list_sources=[genes], gene_coordinate_file=catalog,
                     padding_bp=0, genome_build="hg19", global_config=GlobalConfig(snp_identifier="rsid"),
                     output_dir=output)
    assert "()" not in str(failure.value)
    assert_no_projection_outputs(output)


def test_duplicate_names_keep_resolver_reasons_and_identifier_audit(tmp_path):
    grid, catalog = grid_inputs(tmp_path, "annot")
    left, right = tmp_path / "left", tmp_path / "right"
    left.mkdir()
    right.mkdir()
    sources = [left / "same.txt", right / "same.tsv"]
    sources[0].write_text("G1\n")
    sources[1].write_text("UNKNOWN\n")
    output = tmp_path / "out"
    with pytest.raises(LDSCInputError) as failure:
        run_annotate(**grid, query_annot_gene_list_sources=sources, gene_coordinate_file=catalog,
                     padding_bp=0, genome_build="hg19", global_config=GlobalConfig(snp_identifier="rsid"),
                     output_dir=output)
    assert "duplicate_query_name" in str(failure.value)
    assert "UNKNOWN" in str(failure.value)
    summary = pd.read_csv(output / "diagnostics/gene_list_resolution_summary.tsv", sep="\t")
    assert summary.source_reasons.tolist() == ["duplicate_query_name", "duplicate_query_name"]
    audit = pd.read_csv(output / "diagnostics/gene_list_audit.tsv.gz", sep="\t")
    assert audit.input_gene.tolist() == ["G1", "UNKNOWN"]
    assert_no_projection_outputs(output)


@pytest.mark.parametrize("kind", ["annot", "bim"])
@pytest.mark.parametrize("policy", ["strict", "resolved-only"])
@pytest.mark.parametrize("rows", ["", "UNKNOWN\n", "G1\nUNKNOWN\n"])
def test_names_are_checked_even_for_empty_or_unresolved_queries(tmp_path, kind, policy, rows):
    grid, catalog = grid_inputs(tmp_path, kind)
    name = "base" if kind == "annot" else "SNP"
    genes, sibling = tmp_path / f"{name}.txt", tmp_path / "usable.txt"
    genes.write_text(rows)
    sibling.write_text("G1\n")
    output = tmp_path / "out"
    with pytest.raises(LDSCInputError, match="annotation_name_collision"):
        run_annotate(**grid, query_annot_gene_list_sources=[genes, sibling], gene_coordinate_file=catalog,
                     padding_bp=0, genome_build="hg19", gene_list_resolution_policy=policy,
                     global_config=GlobalConfig(snp_identifier="rsid"), output_dir=output)
    summary = pd.read_csv(output / "diagnostics/gene_list_resolution_summary.tsv", sep="\t").set_index("query")
    assert summary.loc[name, "source_status"] == "error"
    assert summary.loc[name, "source_reasons"] == "annotation_name_collision"
    assert summary.loc["usable", "source_status"] == "ok"
    assert_no_projection_outputs(output)


@pytest.mark.parametrize("kind", ["annot", "bim"])
@pytest.mark.parametrize("policy", ["strict", "resolved-only"])
def test_valid_and_resolved_only_outputs_keep_membership_and_precision(tmp_path, kind, policy):
    grid, catalog = grid_inputs(tmp_path, kind)
    genes = tmp_path / "focal.txt"
    genes.write_text("G1\nUNKNOWN\n" if policy == "resolved-only" else "G1\n")
    output = tmp_path / "out"
    with run_annotate(**grid, query_annot_gene_list_sources=[genes], gene_coordinate_file=catalog,
                      padding_bp=0, genome_build="hg19", gene_list_resolution_policy=policy,
                      global_config=GlobalConfig(snp_identifier="rsid"), output_dir=output) as bundle:
        assert bundle.query_columns == ["focal"]
        assert bundle.baseline_columns == (["base"] if kind == "annot" else [])
        # G1 is one-based inclusive [101,110]; positions 100 and 111 are outside.
        values = bundle.read("1", columns=["focal"])
        assert str(values.dtype) == "float32"
        assert values[:, 0].tolist() == [0, 1, 1, 0]
        assert bundle.read("2", columns=["focal"])[:, 0].tolist() == [0]
    for chrom in ("1", "2"):
        query = pd.read_csv(output / f"query.{chrom}.annot.gz", sep="\t")
        if kind == "annot":
            baseline = pd.read_csv(output / f"baseline.{chrom}.annot.gz", sep="\t")
            pd.testing.assert_frame_equal(query.drop(columns="focal"), baseline.drop(columns="base"))
            assert baseline.base.eq(1).all()
        else:
            assert not list(output.glob("baseline.*.annot.gz"))
    summary = pd.read_csv(output / "diagnostics/gene_list_resolution_summary.tsv", sep="\t")
    assert summary.source_status.tolist() == ["ok"]
    assert summary.rejected_rows.tolist() == [int(policy == "resolved-only")]
    audit = pd.read_csv(output / "diagnostics/gene_list_audit.tsv.gz", sep="\t")
    assert audit.input_gene.tolist() == genes.read_text().splitlines()
    assert audit.disposition.tolist() == (["retained", "rejected"] if policy == "resolved-only" else ["retained"])


def test_invalid_headers_stop_before_dependent_gene_checks(tmp_path, monkeypatch):
    from ldsc._annotation import annotate as annotate_workflow

    grid, catalog = grid_inputs(tmp_path, "annot")
    pd.DataFrame({"unusable": [1]}).to_csv(grid["ref_snp_annot_sources"][0], sep="\t", index=False)
    genes = tmp_path / "base.txt"
    genes.write_text("UNKNOWN\n")

    def unexpected_resolution(*args, **kwargs):
        pytest.fail("Unusable annotation headers must stop dependent preflight")

    monkeypatch.setattr(annotate_workflow, "resolve_gene_lists_staged", unexpected_resolution)
    output = tmp_path / "out"
    with pytest.raises(LDSCInputError, match="Input preflight"):
        run_annotate(**grid, query_annot_gene_list_sources=[genes], gene_coordinate_file=catalog,
                     padding_bp=0, genome_build="hg19", global_config=GlobalConfig(snp_identifier="rsid"),
                     output_dir=output)
    assert (output / "diagnostics/input_issues.tsv").is_file()
    assert not (output / "diagnostics/gene_list_resolution_summary.tsv").exists()
    assert_no_projection_outputs(output)
