from pathlib import Path
import sys
import unittest
from unittest import mock

SRC = Path(__file__).resolve().parents[1] / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))


class LoggingRefactorTest(unittest.TestCase):
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
