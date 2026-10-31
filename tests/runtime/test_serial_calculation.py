"""Error reporting across the serial-calculation process boundary."""

import argparse
from concurrent.futures.process import BrokenProcessPool
import sys
from pathlib import Path

import pytest

import ldsc
from ldsc._regression.quantile_h2 import run_quantile_h2_from_args
from ldsc._runtime import serial_calculation as serial_calculation_module
from ldsc._runtime.workflow_logging import overwrite_failure_marker, workflow_logging
from ldsc._runtime.serial_calculation import serial_calculation
from ldsc.cli import run_cli
from ldsc.errors import LDSCInputError

PLINK_FIXTURE = Path(__file__).parents[1] / "fixtures/minimal_external_resources/plink/hm3_chr22_subset"


@serial_calculation
def _raise_expected_input_error():
    raise LDSCInputError("synthetic expected input error")


@serial_calculation
def _fail_after_opening_a_log(log_path):
    with workflow_logging("unit", log_path, log_level="INFO"):
        raise LDSCInputError("synthetic failure after the log opened")


class _CrashedWorkerPool:
    """Stand-in executor whose only worker dies, as when the system kills it for memory."""

    def __init__(self, *args, **kwargs):
        pass

    def __enter__(self):
        return self

    def __exit__(self, *args):
        return False

    def submit(self, *args, **kwargs):
        raise BrokenProcessPool("A child process terminated abruptly")


def test_failure_marker_names_the_log_a_worker_opened(tmp_path):
    log_path = tmp_path / "diagnostics" / "unit.log"
    log_path.parent.mkdir()

    with pytest.raises(LDSCInputError):
        with overwrite_failure_marker(tmp_path, overwrite=True, command="unit"):
            _fail_after_opening_a_log(log_path)

    assert f"Detailed log: {log_path}\n" in (tmp_path / "RUN_FAILED.txt").read_text(encoding="utf-8")


def test_cli_failure_marker_names_the_log_of_a_worker_run(tmp_path):
    output = tmp_path / "panel"

    assert run_cli([
        "build-r2-matrix", "--plink-prefix", str(PLINK_FIXTURE), "--output-dir", str(output),
        "--ld-wind-kb", "1000", "--maf-min", "0.5", "--source-genome-build", "hg19", "--overwrite",
    ]) == 1

    marker = (output / "RUN_FAILED.chr22.txt").read_text(encoding="utf-8")
    assert f"Detailed log: {output / 'diagnostics/build-r2-matrix.chr22.log'}\n" in marker


@pytest.mark.parametrize("workflow", ["convert-legacy-ldscores", "quantile-h2"])
def test_a_crashed_worker_leaves_a_failure_marker(tmp_path, monkeypatch, workflow):
    monkeypatch.setattr(serial_calculation_module, "ProcessPoolExecutor", _CrashedWorkerPool)
    output = tmp_path / "out"

    with pytest.raises(BrokenProcessPool):
        if workflow == "convert-legacy-ldscores":
            ldsc.convert_legacy_ldscores(
                legacy_reference_dir=tmp_path / "reference", legacy_weight_dir=tmp_path / "weights",
                output_dir=output, overwrite=True,
            )
        else:
            run_quantile_h2_from_args(argparse.Namespace(output_dir=str(output), overwrite=True, log_level="INFO"))

    assert "BrokenProcessPool" in (output / "RUN_FAILED.txt").read_text(encoding="utf-8")


def test_child_traceback_reaches_the_workflow_log_footer_without_an_error_record(tmp_path):
    log_path = tmp_path / "unit.log"

    with pytest.raises(LDSCInputError) as caught:
        with workflow_logging("unit", log_path, log_level="INFO"):
            _raise_expected_input_error()

    assert str(caught.value) == "synthetic expected input error"
    log = log_path.read_text(encoding="utf-8")
    assert "Calculation failed" not in log
    assert "in _raise_expected_input_error" in log


def test_cli_reports_an_expected_worker_error_on_one_line(tmp_path, monkeypatch, capsys):
    output = tmp_path / "panel"
    monkeypatch.setattr(sys, "argv", [
        "ldsc", "build-r2-matrix", "--plink-prefix", str(PLINK_FIXTURE), "--output-dir", str(output),
        "--ld-wind-kb", "1000", "--maf-min", "0.5", "--source-genome-build", "hg19",
    ])

    assert run_cli() == 1

    lines = capsys.readouterr().err.strip().splitlines()
    assert len(lines) == 1
    assert lines[0].startswith("[ERROR] Error: PLINK BED loading retained no SNPs")
    assert "_refpanel/plink.py" in (output / "diagnostics/build-r2-matrix.chr22.log").read_text(encoding="utf-8")
