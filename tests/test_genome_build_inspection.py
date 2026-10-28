"""Read-only build inspection through the public CLI."""

import json
from pathlib import Path
import subprocess
import sys

import pandas as pd
import pytest
import pyarrow as pa
import pyarrow.parquet as pq

from ldsc import cli
from ldsc._genome.build_inference import load_packaged_reference_table


RSID_NOTICE = "This artifact's workflow uses rsIDs for SNP matching."


def coordinates(build="hg19", chrom=1, n=500, offset=0):
    # Fixed packaged reference coordinates establish the expected assembly/basis.
    reference = load_packaged_reference_table()
    rows = reference.loc[reference.CHR.astype(str).eq(str(chrom))].head(n)
    return pd.DataFrame({"CHR": rows.CHR, "POS": rows[f"{build}_POS"] + offset})


def test_parquet_declaration_and_rsid_notice(tmp_path, capsys):
    path = tmp_path / "sumstats.parquet"
    # The declaration is authoritative even when coordinates are too few to infer.
    table = pa.table({"CHR": [1], "POS": [42], "SNP": ["rs1"]})
    pq.write_table(table.replace_schema_metadata({
        b"ldsc:artifact_type": b"sumstats",
        b"ldsc:genome_build": b"GRCh38",
        b"ldsc:snp_identifier": b"rsid",
    }), path)
    before = path.read_bytes()

    assert cli.run_cli(["infer-genome-build", str(path)]) == 0

    captured = capsys.readouterr()
    assert captured.out == (
        "Genome build: hg38 (declared in artifact metadata).\n"
        f"{RSID_NOTICE}\n"
    )
    assert captured.err == ""
    assert path.read_bytes() == before
    assert list(tmp_path.iterdir()) == [path]


@pytest.mark.parametrize("build", ["hg19", "hg38"])
@pytest.mark.parametrize("offset,basis", [(0, "1-based"), (-1, "0-based")])
@pytest.mark.parametrize("suffix", [".bim", ".tsv.gz", ".csv", ".parquet"])
def test_inference_from_coordinates(tmp_path, capsys, build, offset, basis, suffix):
    path = tmp_path / ("input" + suffix)
    frame = coordinates(build, offset=offset)
    if suffix == ".bim":
        frame.assign(SNP="rs1", CM=0, A1="A", A2="C").loc[
            :, ["CHR", "SNP", "CM", "POS", "A1", "A2"]
        ].to_csv(path, sep="\t", index=False, header=False)
    elif suffix == ".parquet":
        pq.write_table(pa.Table.from_pandas(frame).replace_schema_metadata({
            b"ldsc:genome_build": b"", b"ldsc:snp_identifier": b"rsid_allele_aware",
        }), path)
    else:
        frame.rename(columns={"CHR": "chromosome", "POS": "bp"}).assign(SNP="rs1").to_csv(
            path, sep="," if suffix == ".csv" else "\t", index=False,
        )
    assert cli.run_cli(["infer-genome-build", str(path)]) == 0
    output = capsys.readouterr().out
    assert f"Genome build: {build} (inferred from SNP coordinates; {basis})." in output
    assert (RSID_NOTICE in output) == (suffix == ".parquet")


@pytest.mark.parametrize("name,content", [
    ("small.tsv", "CHR\tPOS\n1\t100\n"),
    ("ids.txt", "SNP\nrs1\nrs2\n"),
    ("intervals.bed", "1\t100\t200\n"),
    ("empty.tsv", "CHR\tPOS\n"),
])
def test_readable_insufficient_inputs_are_unknown(tmp_path, capsys, name, content):
    path = tmp_path / name
    path.write_text(content)
    assert cli.run_cli(["infer-genome-build", str(path)]) == 0
    output = capsys.readouterr()
    assert "Genome build: unknown (" in output.out
    assert "--genome-build" not in output.out
    assert output.err == ""


@pytest.mark.parametrize("other_n,other_build,expected_code", [(2, "hg19", 0), (500, "hg38", 1)])
def test_all_chromosomes_participate(tmp_path, capsys, caplog, other_n, other_build, expected_code):
    path = tmp_path / "multi.tsv"
    pd.concat([coordinates(), coordinates(other_build, chrom=2, n=other_n)]).to_csv(path, sep="\t", index=False)
    assert cli.run_cli(["infer-genome-build", str(path)]) == expected_code
    output = capsys.readouterr().out
    if expected_code == 0:
        assert "Genome build: unknown" in output
    else:
        assert "Conflicting genome-build evidence" in caplog.text
        assert "chr1" in caplog.text and "chr2" in caplog.text


def test_separate_named_build_columns_are_alternatives(tmp_path, capsys):
    path = tmp_path / "both.tsv"
    path.write_text("CHR\thg19_POS\thg38_POS\n1\t100\t200\n")
    assert cli.run_cli(["infer-genome-build", str(path)]) == 0
    output = capsys.readouterr().out
    assert output.index("hg19") < output.index("hg38")
    assert "hg19_POS" in output and "hg38_POS" in output
    assert output.count("Genome build") == 2


def write_metadata(root, payload, relative="diagnostics/metadata.json"):
    path = root / relative
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload))
    return path


@pytest.mark.parametrize("target", ["directory", "member", "metadata"])
def test_annotation_metadata_owns_its_member(tmp_path, capsys, target):
    member = tmp_path / "query.1.annot.gz"
    coordinates("hg38", n=2).to_csv(member, sep="\t", index=False)
    metadata = write_metadata(tmp_path, {
        "artifact_type": "annotation_projection", "genome_build": None,
        "projection_genome_build": "hg19", "snp_identifier": "rsid_allele_aware",
        "files": {"query_annotations": [member.name]},
    })
    path = {"directory": tmp_path, "member": member, "metadata": metadata}[target]
    assert cli.run_cli(["infer-genome-build", str(path)]) == 0
    output = capsys.readouterr().out
    assert "Genome build: hg19 (declared in artifact metadata)." in output
    assert RSID_NOTICE in output


def test_unrelated_neighbor_metadata_is_ignored(tmp_path, capsys):
    path = tmp_path / "unrelated.tsv"
    coordinates("hg38").to_csv(path, sep="\t", index=False)
    write_metadata(tmp_path, {
        "artifact_type": "annotation_projection", "genome_build": "hg19", "snp_identifier": "rsid",
        "files": {"query_annotations": ["query.1.annot.gz"]},
    })
    assert cli.run_cli(["infer-genome-build", str(path)]) == 0
    output = capsys.readouterr().out
    assert "Genome build: hg38 (inferred" in output
    assert RSID_NOTICE not in output


@pytest.mark.parametrize("second_build,second_n,expected", [("hg19", 500, "hg19"), ("hg19", 2, "unknown")])
def test_ldscore_directory_infers_all_declared_members(tmp_path, capsys, second_build, second_n, expected):
    coordinates().to_parquet(tmp_path / "ldscore.baseline.parquet", index=False)
    coordinates(second_build, chrom=2, n=second_n).to_parquet(tmp_path / "ldscore.query.parquet", index=False)
    write_metadata(tmp_path, {
        "artifact_type": "ldscore", "genome_build": None, "snp_identifier": "rsid",
        "files": {"baseline": "ldscore.baseline.parquet", "query": "ldscore.query.parquet",
                  "overlap": "unrelated-missing.parquet"},
    }, "metadata.json")
    assert cli.run_cli(["infer-genome-build", str(tmp_path)]) == 0
    output = capsys.readouterr().out
    assert f"Genome build: {expected}" in output
    assert RSID_NOTICE in output


def test_sumstats_directory_detects_trait_named_footer(tmp_path, capsys):
    table = pa.table({"SNP": ["rs1"]}).replace_schema_metadata({
        b"ldsc:artifact_type": b"sumstats", b"ldsc:genome_build": b"", b"ldsc:snp_identifier": b"rsid",
    })
    pq.write_table(table, tmp_path / "trait.parquet")
    (tmp_path / "notes.txt").write_text("Unrelated notes")
    assert cli.run_cli(["infer-genome-build", str(tmp_path)]) == 0
    output = capsys.readouterr().out
    assert "Genome build: unknown" in output and RSID_NOTICE in output


@pytest.mark.parametrize("artifact", [
    "h2_result", "partitioned_h2_result", "rg_result", "quantile_h2_result", "query_r2_result",
    "plot_result", "h2_scale_conversion_result", "partitioned_h2_query_result", "rg_pair_result",
])
def test_saved_results_use_only_recorded_build(tmp_path, capsys, artifact):
    write_metadata(tmp_path, {
        "artifact_type": artifact, "source_genome_build": "hg19", "snp_identifier": "rsid",
        "ldscore_dir": "/must/not/follow/upstream", "files": {"summary": "missing.tsv"},
    })
    assert cli.run_cli(["infer-genome-build", str(tmp_path)]) == 0
    output = capsys.readouterr().out
    assert "Genome build: unknown" in output and RSID_NOTICE in output


def test_gene_index_root_and_chromosome_metadata(tmp_path, capsys):
    payload = {"artifact_type": "gene_ldscore_index", "genome_build": "hg19", "snp_identifier": "rsid"}
    write_metadata(tmp_path, {**payload, "chromosomes": ["1"], "index_identity": {"genome_build": "hg19"}}, "metadata.json")
    write_metadata(tmp_path, payload, "chromosomes/chr1/metadata.json")
    assert cli.run_cli(["infer-genome-build", str(tmp_path)]) == 0
    assert "Genome build: hg19 (declared" in capsys.readouterr().out


def test_r2_directory_reports_both_output_builds(tmp_path, capsys):
    files = {}
    for build in ("hg19", "hg38"):
        directory = tmp_path / build
        directory.mkdir()
        member = directory / "chr1_r2.parquet"
        pq.write_table(pa.table({"IDX_1": [0]}).replace_schema_metadata({
            b"ldsc:artifact_type": b"ref_panel_r2", b"ldsc:genome_build": build.encode(),
            b"ldsc:snp_identifier": b"rsid",
        }), member)
        files[f"r2_{build}"] = [str(member.relative_to(tmp_path))]
    write_metadata(tmp_path, {
        "artifact_type": "ref_panel", "source_genome_build": "hg19", "genome_build": "auto",
        "emitted_genome_builds": ["hg38", "hg19"], "snp_identifier": "rsid", "files": files,
    })
    for path in (tmp_path, tmp_path / "hg38"):
        assert cli.run_cli(["infer-genome-build", str(path)]) == 0
        output = capsys.readouterr().out
        assert "hg38" in output and RSID_NOTICE in output
        if path == tmp_path:
            assert output.index("hg19") < output.index("hg38")
            assert output.count("Genome build") == 2
        else:
            assert "hg19" not in output


def test_saved_result_member_uses_owning_metadata(tmp_path, capsys):
    member = tmp_path / "h2.tsv"
    member.write_text("h2\th2_se\n0.2\t0.03\n")
    write_metadata(tmp_path, {"artifact_type": "h2_result", "genome_build": "hg38",
                              "snp_identifier": "rsid", "files": {"summary": "h2.tsv"}})
    assert cli.run_cli(["infer-genome-build", str(member)]) == 0
    output = capsys.readouterr().out
    assert "Genome build: hg38 (declared" in output and RSID_NOTICE in output


@pytest.mark.parametrize("name,content", [
    ("broken.parquet", b"not parquet"),
    ("broken.bim", b"1 rs1 0 123 A\n"),
    ("broken.tsv.gz", b"not gzip"),
    ("binary.bed", b"\x6c\x1b\x01\x00"),
    ("binary.txt", b"\x00\x01\x02"),
    ("broken.json", b"{"),
    ("unknown.json", b'{"genome_build":"hg19"}'),
    ("malformed.json", b'{"artifact_type":[]}'),
])
def test_malformed_or_unsupported_inputs_fail(tmp_path, capsys, name, content):
    path = tmp_path / name
    path.write_bytes(content)
    assert cli.run_cli(["infer-genome-build", str(path)]) == 1
    assert capsys.readouterr().out == ""
    assert path.read_bytes() == content


def test_conflicting_footer_and_sidecar_declarations(tmp_path, caplog):
    import gzip

    path = tmp_path / "chr1_r2.parquet"
    pq.write_table(pa.table({"IDX_1": [0]}).replace_schema_metadata({
        b"ldsc:genome_build": b"hg19", b"ldsc:sorted_by_build": b"hg19",
    }), path)
    with gzip.open(tmp_path / "chr1_meta.tsv.gz", "wt") as stream:
        stream.write("# ldsc:genome_build=hg38\nCHR\tPOS\n1\t42\n")
    assert cli.run_cli(["infer-genome-build", str(path)]) == 1
    assert "Conflicting genome-build evidence" in caplog.text
    assert "chr1_r2.parquet" in caplog.text and "chr1_meta.tsv.gz" in caplog.text


def test_low_dominance_is_unknown(tmp_path, capsys):
    path = tmp_path / "ambiguous.tsv"
    pd.concat([coordinates("hg19"), coordinates("hg38")]).to_csv(path, sep="\t", index=False)
    assert cli.run_cli(["infer-genome-build", str(path)]) == 0
    assert "Genome build: unknown" in capsys.readouterr().out


def test_sampling_is_bounded_but_later_chromosomes_are_discovered(tmp_path, capsys):
    path = tmp_path / "bounded.tsv"
    # The compact reference has 500 chr1 SNPs; repeat them to fill 5,000 rows.
    first = pd.concat([coordinates("hg19")] * 10)
    later_same_chrom = coordinates("hg38", n=500)
    pd.concat([first, later_same_chrom]).to_csv(path, sep="\t", index=False)
    assert cli.run_cli(["infer-genome-build", str(path)]) == 0
    assert "Genome build: hg19" in capsys.readouterr().out
    coordinates("hg38", chrom=2).to_csv(path, mode="a", sep="\t", index=False, header=False)
    assert cli.run_cli(["infer-genome-build", str(path)]) == 1


@pytest.mark.parametrize("token", ["missing.tsv", "../outside.tsv", "/must/not/read.tsv"])
def test_selected_manifest_members_must_exist_inside_artifact(tmp_path, token):
    write_metadata(tmp_path, {"artifact_type": "annotation_projection", "genome_build": "hg19",
                              "files": {"query_annotations": [token]}})
    assert cli.run_cli(["infer-genome-build", str(tmp_path)]) == 1


def test_metadata_symlink_cannot_escape_artifact(tmp_path):
    outside = tmp_path / "outside.json"
    outside.write_text('{"artifact_type":"h2_result","genome_build":"hg19"}')
    directory = tmp_path / "artifact"
    directory.mkdir()
    (directory / "metadata.json").symlink_to(outside)
    assert cli.run_cli(["infer-genome-build", str(directory)]) == 1


@pytest.mark.parametrize("entrypoint", ["module", "console"])
def test_cli_entrypoints_streams_exits_and_read_only(tmp_path, entrypoint):
    path = tmp_path / "input.tsv"
    path.write_text("SNP\nrs1\n")
    command = ([sys.executable, "-B", "-m", "ldsc"] if entrypoint == "module"
               else [str(Path(sys.executable).with_name("ldsc"))])
    before = {p.relative_to(tmp_path): (p.read_bytes(), p.stat().st_mtime_ns) for p in tmp_path.rglob("*") if p.is_file()}
    for args, expected in [
        ([str(path)], 0), (["--help"], 0), ([str(tmp_path / "missing")], 1),
        ([], 2), ([str(path), str(path)], 2), ([str(path), "--genome-build", "hg19"], 2),
        ([str(tmp_path / "*.tsv")], 1), ([str(tmp_path / "input.@.tsv")], 1),
    ]:
        result = subprocess.run([*command, "infer-genome-build", *args], cwd=tmp_path, capture_output=True, text=True)
        assert result.returncode == expected, result.stderr
        if expected == 0:
            assert result.stderr == ""
            assert "unknown" in result.stdout
        else:
            assert result.stderr and not result.stdout
    after = {p.relative_to(tmp_path): (p.read_bytes(), p.stat().st_mtime_ns) for p in tmp_path.rglob("*") if p.is_file()}
    assert after == before


def test_internal_failure_returns_two(tmp_path, monkeypatch, caplog):
    from ldsc._genome import build_inspection as _genome_build_inspection

    path = tmp_path / "input.tsv"
    path.write_text("SNP\nrs1\n")
    def fail(*args, **kwargs):
        raise RuntimeError("unexpected internal failure")
    monkeypatch.setattr(_genome_build_inspection, "_inspect_table", fail)
    assert cli.run_cli(["infer-genome-build", str(path)]) == 2
    assert "--log-level" not in caplog.text


def test_rejected_overwrite_options_cannot_create_failure_marker(tmp_path):
    path = tmp_path / "input.tsv"
    path.write_text("SNP\nrs1\n")
    assert cli.run_cli(["infer-genome-build", str(path), "--overwrite", "--output-dir", str(tmp_path)]) == 2
    assert list(tmp_path.iterdir()) == [path]


def test_sumstats_directory_ignores_unrelated_invalid_parquet(tmp_path, capsys):
    pq.write_table(pa.table({"SNP": ["rs1"]}).replace_schema_metadata({
        b"ldsc:artifact_type": b"sumstats", b"ldsc:genome_build": b"hg19",
    }), tmp_path / "trait.parquet")
    (tmp_path / "unrelated.parquet").write_bytes(b"not a parquet artifact")
    assert cli.run_cli(["infer-genome-build", str(tmp_path)]) == 0
    assert "Genome build: hg19" in capsys.readouterr().out


def test_panel_manifest_normalizes_emitted_build_aliases(tmp_path, capsys):
    (tmp_path / "hg38").mkdir()
    pq.write_table(pa.table({"IDX_1": [0]}), tmp_path / "hg38/chr1_r2.parquet")
    write_metadata(tmp_path, {"artifact_type": "ref_panel", "emitted_genome_builds": ["GRCh38"],
                              "files": {"r2_hg38": ["hg38/chr1_r2.parquet"]}})
    assert cli.run_cli(["infer-genome-build", str(tmp_path)]) == 0
    assert "Genome build [hg38]: hg38 (declared" in capsys.readouterr().out


@pytest.mark.parametrize("kind,folder", [("partitioned_h2_query_result", "query_annotations"), ("rg_pair_result", "pairs")])
def test_nested_result_member_uses_aggregate_relative_manifest(tmp_path, capsys, kind, folder):
    directory = tmp_path / "diagnostics" / folder / "example"
    directory.mkdir(parents=True)
    member = directory / "summary.tsv"
    member.write_text("h2\tse\n0.2\t0.03\n")
    write_metadata(directory, {"artifact_type": kind, "genome_build": "hg19", "snp_identifier": "rsid",
                               "files": {"summary": str(member.relative_to(tmp_path))}}, "metadata.json")
    assert cli.run_cli(["infer-genome-build", str(member)]) == 0
    output = capsys.readouterr().out
    assert "Genome build: hg19 (declared" in output and RSID_NOTICE in output


def test_malformed_index_chromosomes_are_an_input_error(tmp_path):
    write_metadata(tmp_path, {"artifact_type": "gene_ldscore_index", "genome_build": "hg19",
                              "chromosomes": None}, "metadata.json")
    assert cli.run_cli(["infer-genome-build", str(tmp_path)]) == 1


def test_paired_sidecar_cannot_escape_its_artifact(tmp_path):
    import gzip

    outside = tmp_path / "outside.tsv.gz"
    with gzip.open(outside, "wt") as stream:
        stream.write("# ldsc:genome_build=hg19\nCHR\tPOS\n1\t123\n")
    directory = tmp_path / "artifact"
    directory.mkdir()
    member = directory / "chr1_r2.parquet"
    pq.write_table(pa.table({"IDX_1": [0]}), member)
    (directory / "chr1_meta.tsv.gz").symlink_to(outside)
    assert cli.run_cli(["infer-genome-build", str(member)]) == 1


def test_unsupported_parquet_encoding_is_an_input_error(tmp_path, monkeypatch):
    member = tmp_path / "unsupported.parquet"
    member.write_bytes(b"PAR1")
    def unsupported(*args, **kwargs):
        raise pa.ArrowNotImplementedError("unsupported Parquet encoding")
    monkeypatch.setattr(pq, "read_schema", unsupported)
    assert cli.run_cli(["infer-genome-build", str(member)]) == 1


def test_real_annotation_workflow_artifact(tmp_path, capsys):
    from ldsc import GlobalConfig, run_annotate

    grid, bed, output = tmp_path / "grid.annot", tmp_path / "regions.bed", tmp_path / "result"
    coordinates().assign(SNP=lambda frame: "rs" + frame.POS.astype(str), base=1).to_csv(grid, sep="\t", index=False)
    bed.write_text("1\t0\t1000000000\n")
    with run_annotate(ref_snp_annot_sources=[grid], query_annot_bed_sources=[bed],
                      global_config=GlobalConfig(snp_identifier="rsid"), genome_build="hg19", output_dir=output):
        pass
    assert cli.run_cli(["infer-genome-build", str(output)]) == 0
    text = capsys.readouterr().out
    assert "Genome build: hg19 (declared" in text and RSID_NOTICE in text
