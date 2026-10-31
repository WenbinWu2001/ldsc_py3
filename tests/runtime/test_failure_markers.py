from __future__ import annotations

from pathlib import Path
import shlex
import shutil
from unittest import mock

import pytest

import ldsc
from ldsc import cli
from ldsc._runtime import workflow_logging
from ldsc._runtime.workflow_logging import materializing_overwrite_guard, overwrite_failure_marker
from ldsc._sumstats import munger as sumstats_munger
from ldsc.errors import LDSCInputError, LDSCUsageError, LDSCUserError

PLINK_FIXTURE = Path(__file__).parents[1] / "fixtures/minimal_external_resources/plink/hm3_chr22_subset"

# One run per materializing command that passes argument parsing and then fails
# inside the workflow, with the marker it leaves relative to "{root}".
FAILING_CLI_RUNS = {
    "annotate": (["--ref-snp-annot-sources", "nope.@.annot.gz", "--query-annot-bed-sources", "nope.bed",
                  "--genome-build", "hg38", "--output-dir", "{root}"], "RUN_FAILED.txt"),
    "ldscore": (["--output-dir", "{root}"], "RUN_FAILED.txt"),
    "build-r2-matrix": (["--plink-prefix", str(PLINK_FIXTURE), "--ld-wind-kb", "1000", "--maf-min", "0.5",
                         "--source-genome-build", "hg19", "--output-dir", "{root}"], "RUN_FAILED.chr22.txt"),
    "build-gene-ldscore-index": (["--snp-identifier", "chr_pos", "--genome-build", "hg38",
                                  "--baseline-annot-sources", "nope.@.annot.gz", "--gene-coordinate-file", "nope.tsv",
                                  "--plink-prefix", "nope.@", "--ld-wind-kb", "100", "--output-dir", "{root}"],
                                 "RUN_FAILED.txt"),
    "convert-legacy-ldscores": (["--legacy-reference-dir", "nope", "--legacy-weight-dir", "nope",
                                 "--output-dir", "{root}"], "RUN_FAILED.txt"),
    "munge-sumstats": (["--raw-sumstats-file", "nope.tsv", "--output-dir", "{root}"], "RUN_FAILED.txt"),
    "h2": (["--sumstats-file", "nope.parquet", "--ldscore-dir", "nope", "--output-dir", "{root}"], "RUN_FAILED.txt"),
    "partitioned-h2": (["--sumstats-file", "nope.parquet", "--ldscore-dir", "nope", "--output-dir", "{root}"],
                       "RUN_FAILED.txt"),
    "rg": (["--sumstats-sources", "nope1.parquet,nope2.parquet", "--ldscore-dir", "nope", "--output-dir", "{root}"],
           "RUN_FAILED.txt"),
    "quantile-h2": (["--partitioned-h2-result-dir", "nope", "--target-annot-sources", "nope.@.annot.gz",
                     "--target-annot-column", "x", "--baseline-annot-sources", "nope.@.annot.gz",
                     "--ref-metadata-sources", "nope.@.tsv", "--output-dir", "{root}"], "RUN_FAILED.txt"),
    "query-r2-matrix": (["--r2-matrix-dir", "nope", "--pairs-file", "nope.tsv", "--output-dir", "{root}"],
                        "RUN_FAILED.txt"),
    "convert-h2-scale": (["--h2-result-dir", "{root}", "--samp-prev", "0.5", "--pop-prev", "0.01"],
                         "postprocessing/liability-scale/RUN_FAILED.txt"),
    "plot": (["--result-dir", "{root}"], "plots/RUN_FAILED.txt"),
}


def _record_marker_writes(monkeypatch) -> list[Path]:
    """Record every failure-marker write; the design writes one marker per failed run."""
    writes = []
    write = workflow_logging._atomic_write_failure_marker

    def recording(path, **kwargs):
        writes.append(Path(path))
        return write(path, **kwargs)

    monkeypatch.setattr(workflow_logging, "_atomic_write_failure_marker", recording)
    return writes


def test_authorized_overwrite_failure_writes_durable_marker_without_rollback(tmp_path):
    existing = tmp_path / "result.tsv"
    existing.write_text("partially replaced", encoding="utf-8")

    with pytest.raises(RuntimeError, match="publication failed"):
        with overwrite_failure_marker(
            tmp_path,
            overwrite=True,
            command="ldsc example --output-dir result --overwrite",
            log_path=tmp_path / "diagnostics" / "example.log",
        ):
            raise RuntimeError("publication failed")

    marker = tmp_path / "RUN_FAILED.txt"
    text = marker.read_text(encoding="utf-8")
    assert "publication failed" in text
    assert "incomplete or mixed" in text
    assert "not restored" in text
    assert "diagnostics/example.log" in text
    assert existing.read_text(encoding="utf-8") == "partially replaced"


def test_successful_retry_removes_marker_and_no_overwrite_failure_does_not_create_one(tmp_path):
    marker = tmp_path / "RUN_FAILED.txt"
    marker.write_text("old failure", encoding="utf-8")

    with overwrite_failure_marker(tmp_path, overwrite=True, command="retry"):
        pass
    assert not marker.exists()

    with pytest.raises(RuntimeError):
        with overwrite_failure_marker(tmp_path, overwrite=False, command="fresh run"):
            raise RuntimeError("collision")
    assert not marker.exists()


@pytest.mark.parametrize("command", FAILING_CLI_RUNS)
def test_cli_failure_leaves_one_marker_naming_the_command(tmp_path, monkeypatch, command):
    monkeypatch.chdir(tmp_path)
    writes = _record_marker_writes(monkeypatch)
    options, marker_name = FAILING_CLI_RUNS[command]
    root = tmp_path / "out"
    # An existing empty directory: plot and convert-h2-scale need their result directory to exist.
    root.mkdir()
    argv = [command, *(token.format(root=root) for token in options), "--overwrite"]

    assert cli.run_cli(argv) == 1

    assert writes == [root / marker_name]
    assert f"Command: {shlex.join(['ldsc', *argv])}\n" in (root / marker_name).read_text(encoding="utf-8")


# Python entry points that call a parsed-argument workflow, each failing inside it.
FAILING_PYTHON_CALLS = {
    "run_h2": lambda root: ldsc.run_h2(
        sumstats_file="nope.parquet", ldscore_dir="nope", output_dir=str(root), overwrite=True,
    ),
    "run_ldscore": lambda root: ldsc.run_ldscore(output_dir=str(root), overwrite=True),
    "run_build_r2_matrix": lambda root: ldsc.run_build_r2_matrix(
        plink_prefix=str(PLINK_FIXTURE), output_dir=str(root), source_genome_build="hg19", ld_wind_kb=1000,
        maf_min=0.5, overwrite=True,
    ),
    "build_gene_ldscore_index": lambda root: ldsc.build_gene_ldscore_index(
        ldsc.GeneLDScoreIndexBuildConfig(
            baseline_annot_sources=("nope.@.annot.gz",), plink_prefix="nope.@", output_dir=str(root),
            gene_coordinate_file="nope.tsv", genome_build="hg19", snp_identifier="chr_pos", ld_wind_kb=100,
        ),
        overwrite=True,
    ),
}


@pytest.mark.parametrize("function", FAILING_PYTHON_CALLS)
def test_python_failure_leaves_one_marker_naming_the_public_call(tmp_path, monkeypatch, function):
    monkeypatch.chdir(tmp_path)
    writes = _record_marker_writes(monkeypatch)

    with pytest.raises(LDSCUserError):
        FAILING_PYTHON_CALLS[function](tmp_path / "out")

    assert len(writes) == 1
    assert f"Command: Python API: {function}(...)\n" in writes[0].read_text(encoding="utf-8")


def test_nested_marker_scopes_write_once_with_the_outermost_label(tmp_path, monkeypatch):
    writes = _record_marker_writes(monkeypatch)

    @materializing_overwrite_guard(lambda output_dir: (output_dir, True), command="inner_call(...)")
    def inner_call(output_dir):
        raise RuntimeError("failed")

    with pytest.raises(RuntimeError, match="failed"):
        with overwrite_failure_marker(tmp_path, overwrite=True, command="Python API: outer_call(...)"):
            inner_call(tmp_path)

    assert writes == [tmp_path / "RUN_FAILED.txt"]
    assert "Command: Python API: outer_call(...)\n" in (tmp_path / "RUN_FAILED.txt").read_text(encoding="utf-8")


def test_a_nested_marker_scope_must_name_the_open_marker(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    with overwrite_failure_marker(tmp_path / "out", overwrite=True, command="outer"):
        with overwrite_failure_marker("out", overwrite=True, command="same marker, relative spelling"):
            pass

    with pytest.raises(RuntimeError, match="different failure marker"):
        with overwrite_failure_marker(tmp_path / "out", overwrite=False, command="outer"):
            with overwrite_failure_marker(tmp_path / "elsewhere", overwrite=False, command="inner"):
                pass


def test_munging_option_errors_leave_a_marker_from_python_too(tmp_path):
    raw = tmp_path / "raw.tsv"
    raw.write_text("SNP A1 A2 P BETA N\nrs1 A G 0.05 0.2 1000\n", encoding="utf-8")
    output_dir = tmp_path / "out"

    with pytest.raises(LDSCUsageError):
        sumstats_munger.main([
            "--raw-sumstats-file", str(raw), "--output-dir", str(output_dir), "--overwrite",
            "--snp-identifier", "rsid", "--source-genome-build", "hg38",
        ])

    assert "Error: LDSCUsageError" in (output_dir / "RUN_FAILED.txt").read_text(encoding="utf-8")


def test_a_failure_after_a_successful_run_leaves_no_marker(tmp_path, monkeypatch):
    raw = tmp_path / "raw.tsv"
    raw.write_text("SNP A1 A2 P BETA N\nrs1 A G 0.05 0.2 1000\nrs2 C T 0.5 -0.1 1000\n", encoding="utf-8")
    output_dir = tmp_path / "out"
    output_dir.mkdir()
    (output_dir / "RUN_FAILED.txt").write_text("previous failed run", encoding="utf-8")
    # The console summary of a completed run fails, as when stdout is piped into `head`.
    monkeypatch.setattr(sumstats_munger, "print", mock.Mock(side_effect=BrokenPipeError("stdout")), raising=False)

    status = cli.run_cli([
        "munge-sumstats", "--raw-sumstats-file", str(raw), "--output-dir", str(output_dir), "--overwrite",
        "--snp-identifier", "rsid", "--no-snp-restriction",
    ])

    assert status == 2
    assert (output_dir / "sumstats.parquet").is_file()
    assert not (output_dir / "RUN_FAILED.txt").exists()


@pytest.mark.parametrize("help_flag", ["--help", "-h"])
@pytest.mark.parametrize("existing", [False, True])
def test_cli_help_with_overwrite_leaves_output_scope_untouched(tmp_path, help_flag, existing):
    output_dir = tmp_path / "out"
    marker = output_dir / "RUN_FAILED.txt"
    if existing:
        output_dir.mkdir()
        marker.write_text("previous failed run", encoding="utf-8")

    status = cli.run_cli(["h2", "--output-dir", str(output_dir), "--overwrite", help_flag])

    assert status == 0
    if existing:
        assert marker.read_text(encoding="utf-8") == "previous failed run"
        assert list(output_dir.iterdir()) == [marker]
    else:
        assert not output_dir.exists()


@pytest.mark.parametrize(
    ("argv", "marker_name"),
    [
        (["h2", "--output-dir", "{root}", "--overwrite", "--bogus-flag"], "RUN_FAILED.txt"),
        (["h2", "--output-dir", "{root}", "--overwrite"], "RUN_FAILED.txt"),
        (["ldscore", "--output-dir", "{root}", "--overwrite", "--bogus-flag"], "RUN_FAILED.txt"),
        (["plot", "--result-dir", "{root}", "--overwrite", "--bogus-flag"], "plots/RUN_FAILED.txt"),
        (["build-r2-matrix", "--plink-prefix", "panel.chr22", "--output-dir", "{root}", "--overwrite",
          "--bogus-flag"], "RUN_FAILED.chr22.txt"),
    ],
)
@pytest.mark.parametrize("existing", [False, True])
def test_cli_argument_error_with_overwrite_leaves_output_scope_untouched(tmp_path, argv, marker_name, existing):
    root = tmp_path / "out"
    marker = root / marker_name
    if existing:
        marker.parent.mkdir(parents=True)
        marker.write_text("previous failed run", encoding="utf-8")

    status = cli.run_cli([token.format(root=root) for token in argv])

    assert status == 2
    if existing:
        assert marker.read_text(encoding="utf-8") == "previous failed run"
        assert [path for path in root.rglob("*") if path.is_file()] == [marker]
    else:
        assert not root.exists()


@pytest.mark.parametrize("code", [0, 2])
def test_marker_helper_leaves_scope_untouched_when_argument_parsing_exits(tmp_path, code):
    output_dir = tmp_path / "out"

    with pytest.raises(SystemExit):
        with overwrite_failure_marker(output_dir, overwrite=True, command="ldsc h2 --bogus-flag"):
            raise SystemExit(code)
    assert not output_dir.exists()

    output_dir.mkdir()
    marker = output_dir / "RUN_FAILED.txt"
    marker.write_text("previous failed run", encoding="utf-8")
    with pytest.raises(SystemExit):
        with overwrite_failure_marker(output_dir, overwrite=True, command="ldsc h2 --bogus-flag"):
            raise SystemExit(code)
    assert marker.read_text(encoding="utf-8") == "previous failed run"


@pytest.mark.parametrize("entrypoint", ["cli", "python"])
@pytest.mark.parametrize("input_state", ["valid", "missing", "omitted"])
@pytest.mark.parametrize("existing", [False, True])
def test_infer_only_with_overwrite_never_mutates_output_scope(
    tmp_path, entrypoint, input_state, existing
):
    from ldsc._sumstats import munger as sumstats_munger

    output_dir = tmp_path / "unused"
    marker = output_dir / "RUN_FAILED.txt"
    if existing:
        output_dir.mkdir()
        marker.write_text("previous failed run", encoding="utf-8")
    raw_path = tmp_path / "raw.tsv"
    if input_state == "valid":
        raw_path.write_text("SNP A1 A2 P BETA N\nrs1 A G 0.05 0.2 1000\n", encoding="utf-8")
    argv = [
        "--output-dir", str(output_dir), "--overwrite", "--infer-only",
        "--snp-identifier", "rsid",
    ]
    if input_state != "omitted":
        argv.extend(["--raw-sumstats-file", str(raw_path)])

    if entrypoint == "cli":
        status = cli.run_cli(["munge-sumstats", *argv])
        assert (status == 0) == (input_state == "valid")
    elif input_state == "valid":
        result = sumstats_munger.main(argv)
        assert result.detected_format == "plain"
    else:
        error = SystemExit if input_state == "omitted" else LDSCInputError
        with pytest.raises(error):
            sumstats_munger.main(argv)

    if existing:
        assert marker.read_text(encoding="utf-8") == "previous failed run"
        assert list(output_dir.iterdir()) == [marker]
    else:
        assert not output_dir.exists()


def test_unknown_cli_command_does_not_claim_an_output_scope(tmp_path):
    output_dir = tmp_path / "out"

    assert cli.run_cli(["not-a-command", "--output-dir", str(output_dir), "--overwrite"]) == 2

    assert not output_dir.exists()


def _copy_panel(prefix: Path, *, chrom: str = "22") -> Path:
    """Copy the chromosome 22 PLINK fixture to ``prefix``, writing ``chrom`` into the BIM file."""
    prefix.parent.mkdir(parents=True, exist_ok=True)
    for suffix in (".bed", ".fam"):
        shutil.copy(f"{PLINK_FIXTURE}{suffix}", f"{prefix}{suffix}")
    rows = Path(f"{PLINK_FIXTURE}.bim").read_text(encoding="utf-8").splitlines()
    Path(f"{prefix}.bim").write_text("".join(f"{chrom}\t{row.split(chr(9), 1)[1]}\n" for row in rows), encoding="utf-8")
    return prefix


def _build_r2_matrix(prefix, output_dir, *extra):
    return cli.run_cli([
        "build-r2-matrix", "--plink-prefix", str(prefix), "--output-dir", str(output_dir),
        "--ld-wind-kb", "1000", "--source-genome-build", "hg19", "--overwrite", *extra,
    ])


@pytest.mark.parametrize("name", ["panel", "EUR_5.chr22", "pan*"])
def test_r2_marker_follows_the_bim_chromosome_not_the_file_name(tmp_path, name):
    _copy_panel(tmp_path / "p" / name.replace("*", "el"))
    output_dir = tmp_path / "out"

    # MAF 0.5 removes every SNP after the inputs are checked, so the run fails inside its marker scope.
    assert _build_r2_matrix(tmp_path / "p" / name, output_dir, "--maf-min", "0.5") == 1

    assert sorted(path.name for path in output_dir.glob("RUN_FAILED*")) == ["RUN_FAILED.chr22.txt"]
    assert (output_dir / "diagnostics" / "build-r2-matrix.chr22.log").is_file()


def test_a_chromosome_job_never_removes_another_chromosomes_marker(tmp_path):
    output_dir = tmp_path / "out"
    assert _build_r2_matrix(_copy_panel(tmp_path / "p" / "panelA", chrom="21"), output_dir, "--maf-min", "0.5") == 1
    assert _build_r2_matrix(_copy_panel(tmp_path / "p" / "panelB"), output_dir) == 0

    assert sorted(path.name for path in output_dir.glob("RUN_FAILED*")) == ["RUN_FAILED.chr21.txt"]
    assert (output_dir / "hg19" / "chr22_meta.tsv.gz").is_file()


@pytest.mark.parametrize("failure", ["missing PLINK input", "no LD window"])
def test_r2_errors_before_the_chromosomes_are_known_leave_no_marker(tmp_path, failure):
    output_dir = tmp_path / "out"
    output_dir.mkdir()
    marker = output_dir / "RUN_FAILED.chr22.txt"
    marker.write_text("previous failed run", encoding="utf-8")
    argv = ["build-r2-matrix", "--plink-prefix", str(tmp_path / "missing.chr22"), "--output-dir", str(output_dir),
            "--source-genome-build", "hg19", "--overwrite"]
    if failure == "missing PLINK input":
        argv += ["--ld-wind-kb", "1000"]

    assert cli.run_cli(argv) == 1

    assert sorted(path.name for path in output_dir.glob("RUN_FAILED*")) == ["RUN_FAILED.chr22.txt"]
    assert marker.read_text(encoding="utf-8") == "previous failed run"
