from pathlib import Path
import sys
import unittest

SRC = Path(__file__).resolve().parents[1] / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))


class PackageLayoutTest(unittest.TestCase):
    def test_public_package_imports(self):
        import ldsc
        import ldsc.__main__

        self.assertTrue(hasattr(ldsc, "LDScoreCalculator"))
        self.assertTrue(hasattr(ldsc, "RegressionRunner"))
        self.assertTrue(hasattr(ldsc, "SumstatsMunger"))
        self.assertTrue(hasattr(ldsc, "AnnotationBuilder"))
        self.assertTrue(hasattr(ldsc.__main__, "main"))

    def test_cli_exposes_expected_subcommands(self):
        from ldsc import cli

        parser = cli.build_parser()
        subparsers_action = next(
            action for action in parser._actions if action.__class__.__name__ == "_SubParsersAction"
        )
        self.assertEqual(
            set(subparsers_action.choices),
            {"annotate", "ldscore", "munge-sumstats", "h2", "partitioned-h2", "rg"},
        )

    def test_ldscore_subcommand_accepts_keep(self):
        from ldsc import cli

        parser = cli.build_parser()
        args = parser.parse_args(
            [
                "ldscore",
                "--out",
                "out/example",
                "--baseline-annot",
                "baseline.annot.gz",
                "--bfile",
                "panel",
                "--ld-wind-snps",
                "10",
                "--keep",
                "samples.keep",
            ]
        )

        self.assertEqual(args.command, "ldscore")
        self.assertEqual(args.keep, "samples.keep")

    def test_ldscore_subcommand_accepts_print_snps(self):
        from ldsc import cli

        parser = cli.build_parser()
        args = parser.parse_args(
            [
                "ldscore",
                "--out",
                "out/example",
                "--baseline-annot",
                "baseline.annot.gz",
                "--bfile",
                "panel",
                "--ld-wind-snps",
                "10",
                "--print-snps",
                "filters/print_snps.txt",
            ]
        )

        self.assertEqual(args.command, "ldscore")
        self.assertEqual(args.print_snps, "filters/print_snps.txt")
