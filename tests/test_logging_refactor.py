from pathlib import Path
import logging
import sys
import tempfile
import unittest
from unittest import mock

SRC = Path(__file__).resolve().parents[1] / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))


class LoggingRefactorTest(unittest.TestCase):
    def test_workflow_logging_captures_ldsc_child_records_and_restores_state(self):
        from ldsc._logging import log_inputs, log_outputs, workflow_logging

        ldsc_logger = logging.getLogger("LDSC")
        original_level = ldsc_logger.level
        original_handlers = list(ldsc_logger.handlers)
        with tempfile.TemporaryDirectory() as tmpdir:
            log_path = Path(tmpdir) / "workflow.log"

            with workflow_logging("unit", log_path, log_level="INFO"):
                log_inputs(source="input.tsv", output_dir=tmpdir)
                logging.getLogger("LDSC.some_child").info("child logger message")
                log_outputs(result="output.tsv")

            text = log_path.read_text(encoding="utf-8")

        self.assertEqual(ldsc_logger.level, original_level)
        self.assertEqual(ldsc_logger.handlers, original_handlers)
        self.assertIn("LDSC unit Started", text)
        self.assertIn("Inputs:", text)
        self.assertIn("source", text)
        self.assertIn("child logger message", text)
        self.assertIn("Outputs:", text)
        self.assertIn("Finished", text)
        self.assertIn("Elapsed time:", text)

    def test_workflow_logging_formats_call_over_multiple_lines(self):
        from ldsc._logging import workflow_logging

        argv = [
            "./munge_sumstats.py",
            "--genome-build",
            "auto",
            "--sumstats",
            "/data/raw.tsv",
            "--out",
            "/out/sumstats",
        ]
        with tempfile.TemporaryDirectory() as tmpdir:
            log_path = Path(tmpdir) / "workflow.log"
            with mock.patch.object(sys, "argv", argv):
                with workflow_logging("unit", log_path, log_level="INFO"):
                    pass

            text = log_path.read_text(encoding="utf-8")

        self.assertIn(
            "Call:\n"
            "./munge_sumstats.py \\\n"
            "--genome-build auto \\\n"
            "--sumstats /data/raw.tsv \\\n"
            "--out /out/sumstats\n",
            text,
        )
        self.assertNotIn("Invocation:", text)

    def test_workflow_logging_leaves_blank_line_before_outputs(self):
        from ldsc._logging import log_outputs, workflow_logging

        with tempfile.TemporaryDirectory() as tmpdir:
            log_path = Path(tmpdir) / "workflow.log"
            with workflow_logging("unit", log_path, log_level="INFO"):
                logging.getLogger("LDSC.unit").info("Metadata:")
                logging.getLogger("LDSC.unit").info("Mean chi^2 = 2.009")
                log_outputs(sumstats_parquet="/out/sumstats.parquet")

            text = log_path.read_text(encoding="utf-8")

        self.assertIn("Mean chi^2 = 2.009\n\nOutputs:", text)

    def test_workflow_logging_formats_elapsed_time_with_units(self):
        from ldsc._logging import workflow_logging

        with tempfile.TemporaryDirectory() as tmpdir:
            log_path = Path(tmpdir) / "workflow.log"
            with mock.patch("ldsc._logging.time.monotonic", side_effect=[100.0, 232.4]):
                with workflow_logging("unit", log_path, log_level="INFO"):
                    pass

            text = log_path.read_text(encoding="utf-8")

        self.assertIn("Elapsed time: 2.0min:12s", text)
        self.assertNotIn("Elapsed: 2:12", text)

    def test_workflow_logging_keeps_audit_lines_at_error_level(self):
        from ldsc._logging import workflow_logging

        with tempfile.TemporaryDirectory() as tmpdir:
            log_path = Path(tmpdir) / "workflow.log"
            with workflow_logging("unit", log_path, log_level="ERROR"):
                logging.getLogger("LDSC.threshold").info("suppressed info")
                logging.getLogger("LDSC.threshold").error("kept error")

            text = log_path.read_text(encoding="utf-8")

        self.assertIn("LDSC unit Started", text)
        self.assertIn("kept error", text)
        self.assertNotIn("suppressed info", text)
        self.assertIn("Finished", text)
        self.assertIn("Elapsed time:", text)

    def test_workflow_logging_failed_footer_does_not_record_exception_message(self):
        from ldsc._logging import workflow_logging

        with tempfile.TemporaryDirectory() as tmpdir:
            log_path = Path(tmpdir) / "workflow.log"
            with self.assertRaisesRegex(RuntimeError, "boom"):
                with workflow_logging("unit", log_path, log_level="INFO"):
                    raise RuntimeError("boom")

            text = log_path.read_text(encoding="utf-8")

        self.assertIn("Failed", text)
        self.assertIn("Elapsed time:", text)
        self.assertNotIn("boom", text)

    def test_package_exports_domain_error_hierarchy(self):
        import ldsc
        from ldsc.config import ConfigMismatchError

        self.assertTrue(issubclass(ldsc.LDSCUserError, ldsc.LDSCError))
        self.assertTrue(issubclass(ldsc.LDSCConfigError, ldsc.LDSCUserError))
        self.assertTrue(issubclass(ConfigMismatchError, ldsc.LDSCConfigError))
        self.assertTrue(issubclass(ConfigMismatchError, ValueError))

    def test_run_cli_reports_user_errors_without_traceback(self):
        from ldsc import cli

        with mock.patch.object(cli, "main", side_effect=ValueError("bad input")):
            with self.assertLogs("LDSC.cli", level="ERROR") as caught:
                code = cli.run_cli(["munge-sumstats"])

        self.assertEqual(code, 1)
        self.assertEqual(len(caught.records), 1)
        record = caught.records[0]
        self.assertEqual(record.levelname, "ERROR")
        self.assertIsNone(record.exc_info)
        self.assertIn("bad input", record.getMessage())

    def test_run_cli_logs_unexpected_errors_with_traceback(self):
        from ldsc import cli

        with mock.patch.object(cli, "main", side_effect=RuntimeError("boom")):
            with self.assertLogs("LDSC.cli", level="ERROR") as caught:
                code = cli.run_cli(["ldscore"])

        self.assertEqual(code, 2)
        self.assertEqual(len(caught.records), 1)
        record = caught.records[0]
        self.assertEqual(record.levelname, "ERROR")
        self.assertIs(record.exc_info[0], RuntimeError)
        self.assertIn("Internal error while running ldsc", record.getMessage())
        self.assertIn("boom", record.getMessage())


if __name__ == "__main__":
    unittest.main()
