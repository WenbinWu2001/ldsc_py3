"""Every file a run touches stays inside its output directory.

An empty output directory would otherwise normalize to the working directory,
where an authorized overwrite can replace or remove the user's own files.
"""

from __future__ import annotations

import gzip
from pathlib import Path

import pytest

import ldsc
from ldsc import cli
from ldsc._annotation.storage import AnnotationWorkspace
from ldsc._runtime.paths import ensure_output_directory
from ldsc._runtime.workflow_logging import overwrite_failure_marker
from ldsc.errors import LDSCConfigError
from tests.runtime.test_failure_markers import FAILING_CLI_RUNS

PLINK_FIXTURE = Path(__file__).parents[1] / "fixtures/minimal_external_resources/plink/hm3_chr22_subset"


@pytest.fixture
def empty_cwd(tmp_path, monkeypatch):
    """Run from an empty working directory, with inputs kept in a separate folder."""
    cwd = tmp_path / "cwd"
    cwd.mkdir()
    monkeypatch.chdir(cwd)
    return cwd


@pytest.mark.parametrize("command", FAILING_CLI_RUNS)
@pytest.mark.parametrize("blank", ["", "  "])
def test_cli_rejects_an_empty_destination_before_writing(empty_cwd, command, blank, capsys):
    options, _ = FAILING_CLI_RUNS[command]

    assert cli.run_cli([command, *(token.format(root=blank) for token in options), "--overwrite"]) == 2

    assert list(empty_cwd.iterdir()) == []
    assert "must name a directory" in capsys.readouterr().err


def _legacy_suite(root: Path) -> dict[str, str]:
    for role in ("reference", "weights"):
        (root / role).mkdir(parents=True)
    for chrom in range(1, 23):
        with gzip.open(root / "reference" / f"{chrom}.l2.ldscore.gz", "wt") as handle:
            handle.write(f"CHR SNP BP L2\n{chrom} rs{chrom} {chrom * 100} {1.0 + chrom / 100}\n")
        with gzip.open(root / "weights" / f"weights.{chrom}.l2.ldscore.gz", "wt") as handle:
            handle.write(f"CHR SNP BP L2\n{chrom} rs{chrom} {chrom * 100 + 1} {2.0 + chrom / 100}\n")
        (root / "reference" / f"{chrom}.l2.M_5_50").write_text("10\n")
    return {"legacy_reference_dir": str(root / "reference"), "legacy_weight_dir": str(root / "weights")}


def _python_call(name: str, inputs: Path, blank: str):
    """Return a Python call that would succeed with valid inputs if the empty destination were accepted."""
    inputs.mkdir(exist_ok=True)
    bed = inputs / "query.bed"
    bed.write_text("chr22\t16000000\t17000000\n")
    raw = inputs / "raw.tsv"
    raw.write_text("SNP A1 A2 P BETA N\nrs1 A G 0.05 0.2 1000\nrs2 C T 0.5 -0.1 1000\n")
    calls = {
        "run_annotate": lambda: ldsc.run_annotate(
            output_dir=blank, ref_snp_bim_sources=[f"{PLINK_FIXTURE}.bim"], query_annot_bed_sources=[str(bed)],
            genome_build="hg38",
        ),
        "convert_legacy_ldscores": lambda: ldsc.convert_legacy_ldscores(output_dir=blank, **_legacy_suite(inputs)),
        "SumstatsMunger.run": lambda: ldsc.SumstatsMunger().run(
            ldsc.MungeConfig(raw_sumstats_file=str(raw), output_dir=blank, no_snp_restriction=True),
            global_config=ldsc.GlobalConfig(snp_identifier="rsid"),
        ),
        "AnnotationBuilder.run": lambda: ldsc.AnnotationBuilder(
            ldsc.AnnotationBuildConfig(query_annot_bed_sources=(str(bed),)),
        ).run(output_dir=blank),
        "plot_result": lambda: ldsc.plot_result(blank, overwrite=True),
        "plot_result output_dir": lambda: ldsc.plot_result(str(inputs), output_dir=blank, overwrite=True),
        "convert_h2_scale": lambda: ldsc.convert_h2_scale(blank, samp_prev=0.5, pop_prev=0.01, overwrite=True),
        "convert_h2_scale output_dir": lambda: ldsc.convert_h2_scale(
            str(inputs), samp_prev=0.5, pop_prev=0.01, output_dir=blank, overwrite=True,
        ),
    }
    return calls[name]


@pytest.mark.parametrize("name", [
    "run_annotate", "convert_legacy_ldscores", "SumstatsMunger.run", "AnnotationBuilder.run",
    "plot_result", "plot_result output_dir", "convert_h2_scale", "convert_h2_scale output_dir",
])
def test_python_rejects_an_empty_destination_before_writing(empty_cwd, tmp_path, name):
    call = _python_call(name, tmp_path / "inputs", "")

    with pytest.raises(LDSCConfigError, match="empty"):
        call()

    assert list(empty_cwd.iterdir()) == []


@pytest.mark.parametrize("helper", [
    lambda: ensure_output_directory(""),
    lambda: AnnotationWorkspace(" "),
    lambda: overwrite_failure_marker("", overwrite=True, command="unit").__enter__(),
])
def test_output_helpers_reject_an_empty_destination(empty_cwd, helper):
    with pytest.raises(LDSCConfigError, match="empty"):
        helper()

    assert list(empty_cwd.iterdir()) == []
