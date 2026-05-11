from __future__ import annotations

import logging
from pathlib import Path
import sys

SRC = Path(__file__).resolve().parents[1] / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from ldsc._kernel.liftover import LiftoverDropReport, log_liftover_drop_report


def test_liftover_drop_report_logs_counts_at_info_and_examples_at_debug(caplog, tmp_path):
    logger = logging.getLogger("LDSC.test_liftover_logging")
    sidecar_path = tmp_path / "dropped_snps" / "dropped.tsv.gz"
    report = LiftoverDropReport(
        reason="target_collision",
        n_dropped=2,
        examples=[
            {"SNP": "rs1", "CHR": "1", "source_POS": 100, "target_POS": 1000},
            {"SNP": "rs2", "CHR": "1", "source_POS": 200, "target_POS": 1000},
        ],
    )

    with caplog.at_level(logging.DEBUG, logger=logger.name):
        log_liftover_drop_report(
            logger,
            report,
            workflow_label="Summary-statistics liftover",
            sidecar_path=sidecar_path,
        )

    info_records = [record for record in caplog.records if record.levelno == logging.INFO]
    debug_records = [record for record in caplog.records if record.levelno == logging.DEBUG]
    assert len(info_records) == 1
    assert len(debug_records) == 1
    info_message = info_records[0].getMessage()
    debug_message = debug_records[0].getMessage()
    assert "Summary-statistics liftover dropped 2 SNPs for target_collision" in info_message
    assert str(sidecar_path) in info_message
    assert "examples=" not in info_message.lower()
    assert "rs1" not in info_message
    assert "target_collision" in debug_message
    assert "rs1" in debug_message
