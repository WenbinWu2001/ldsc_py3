"""Error reporting across the serial-calculation process boundary."""

import sys
from pathlib import Path

import pytest

from ldsc._runtime.workflow_logging import workflow_logging
from ldsc._runtime.serial_calculation import serial_calculation
from ldsc.cli import run_cli
from ldsc.errors import LDSCInputError

PLINK_FIXTURE = Path(__file__).parent / "fixtures/minimal_external_resources/plink/hm3_chr22_subset"


@serial_calculation
def _raise_expected_input_error():
    raise LDSCInputError("synthetic expected input error")


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
        "ldsc", "build-r2-panel", "--plink-prefix", str(PLINK_FIXTURE), "--output-dir", str(output),
        "--ld-wind-kb", "1000", "--maf-min", "0.5", "--source-genome-build", "hg19",
    ])

    assert run_cli() == 1

    lines = capsys.readouterr().err.strip().splitlines()
    assert len(lines) == 1
    assert lines[0].startswith("[ERROR] Error: PLINK BED loading retained no SNPs")
    assert "_refpanel/plink.py" in (output / "diagnostics/build-r2-panel.chr22.log").read_text(encoding="utf-8")
