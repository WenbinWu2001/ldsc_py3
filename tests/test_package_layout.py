from pathlib import Path
import runpy
import sys
import unittest
import warnings
from unittest import mock

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
        self.assertTrue(hasattr(ldsc, "ReferencePanelBuilder"))
        self.assertTrue(hasattr(ldsc, "LDScoreOutputConfig"))
        self.assertTrue(hasattr(ldsc, "load_ldscore_from_dir"))
        self.assertTrue(hasattr(ldsc.__main__, "main"))

    def test_cli_exposes_expected_subcommands(self):
        from ldsc import cli

        parser = cli.build_parser()
        subparsers_action = next(
            action for action in parser._actions if action.__class__.__name__ == "_SubParsersAction"
        )
        self.assertEqual(
            set(subparsers_action.choices),
            {"annotate", "ldscore", "build-ref-panel", "munge-sumstats", "h2", "partitioned-h2", "rg"},
        )

    def test_ldscore_subcommand_accepts_unified_path_flags(self):
        from ldsc import cli

        parser = cli.build_parser()
        args = parser.parse_args(
            [
                "ldscore",
                "--output-dir",
                "out/ldscores",
                "--baseline-annot-paths",
                "baseline.annot.gz",
                "--plink-path",
                "panel",
                "--ld-wind-snps",
                "10",
                "--keep-indivs-path",
                "samples.keep",
            ]
        )

        self.assertEqual(args.command, "ldscore")
        self.assertEqual(args.keep_indivs_path, "samples.keep")

    def test_ldscore_subcommand_accepts_regression_snps_path(self):
        from ldsc import cli

        parser = cli.build_parser()
        args = parser.parse_args(
            [
                "ldscore",
                "--output-dir",
                "out/ldscores",
                "--baseline-annot-paths",
                "baseline.annot.gz",
                "--plink-path",
                "panel",
                "--ld-wind-snps",
                "10",
                "--regression-snps-path",
                "filters/hm3.txt",
            ]
        )

        self.assertEqual(args.command, "ldscore")
        self.assertEqual(args.regression_snps_path, "filters/hm3.txt")

    def test_ldscore_subcommand_rejects_removed_print_snps_and_regression_snps_flags(self):
        from ldsc import cli

        parser = cli.build_parser()
        with self.assertRaises(SystemExit):
            parser.parse_args(
                [
                    "ldscore",
                    "--output-dir",
                    "out/ldscores",
                    "--baseline-annot-paths",
                    "baseline.annot.gz",
                    "--plink-path",
                    "panel",
                    "--ld-wind-snps",
                    "10",
                    "--print-snps",
                    "filters/print_snps.txt",
                ]
            )
        with self.assertRaises(SystemExit):
            parser.parse_args(
                [
                    "ldscore",
                    "--output-dir",
                    "out/ldscores",
                    "--baseline-annot",
                    "baseline.annot.gz",
                    "--bfile",
                    "panel",
                    "--ld-wind-snps",
                    "10",
                    "--regression-snps",
                    "filters/hm3.txt",
                ]
            )

    def test_annotate_subcommand_rejects_removed_baseline_dir_flag(self):
        from ldsc import cli

        parser = cli.build_parser()
        with self.assertRaises(SystemExit) as exc:
            parser.parse_args(
                [
                    "annotate",
                    "--query-annot-bed-paths",
                    "beds/*.bed",
                    "--baseline-annot-dir",
                    "annotations/baseline_chr",
                    "--output-dir",
                    "out/annot",
                ]
            )

        self.assertNotEqual(exc.exception.code, 0)

    def test_ldscore_subcommand_rejects_removed_chr_flags(self):
        from ldsc import cli

        parser = cli.build_parser()
        with self.assertRaises(SystemExit) as exc:
            parser.parse_args(
                [
                    "ldscore",
                    "--output-dir",
                    "out/ldscores",
                    "--baseline-annot-chr",
                    "annotations/baseline.@.annot.gz",
                    "--plink-path",
                    "panel",
                    "--ld-wind-snps",
                    "10",
                ]
            )

        self.assertNotEqual(exc.exception.code, 0)

    def test_ldscore_subcommand_rejects_removed_out_flag(self):
        from ldsc import cli

        parser = cli.build_parser()
        with self.assertRaises(SystemExit) as exc:
            parser.parse_args(
                [
                    "ldscore",
                    "--out",
                    "out/example",
                    "--baseline-annot-paths",
                    "baseline.annot.gz",
                    "--plink-path",
                    "panel",
                    "--ld-wind-snps",
                    "10",
                ]
            )

        self.assertNotEqual(exc.exception.code, 0)

    def test_build_ref_panel_subcommand_rejects_removed_chr_flag(self):
        from ldsc import cli

        parser = cli.build_parser()
        with self.assertRaises(SystemExit) as exc:
            parser.parse_args(
                [
                    "build-ref-panel",
                    "--bfile-chr",
                    "data/reference/panel_chr@",
                    "--panel-label",
                    "TEST",
                    "--source-genome-build",
                    "hg38",
                    "--genetic-map-hg19",
                    "maps/hg19.txt",
                    "--genetic-map-hg38",
                    "maps/hg38.txt",
                    "--liftover-chain-hg38-to-hg19",
                    "chains/hg38ToHg19.over.chain",
                    "--ld-wind-kb",
                    "10",
                    "--out",
                    "out/panel",
                ]
            )

        self.assertNotEqual(exc.exception.code, 0)

    def test_build_ref_panel_subcommand_accepts_output_dir_and_no_panel_label(self):
        from ldsc import cli

        parser = cli.build_parser()
        args = parser.parse_args(
            [
                "build-ref-panel",
                "--plink-path",
                "data/reference/panel_chr@",
                "--source-genome-build",
                "hg38",
                "--genetic-map-hg19-path",
                "maps/hg19.txt",
                "--genetic-map-hg38-path",
                "maps/hg38.txt",
                "--liftover-chain-hg38-to-hg19-path",
                "chains/hg38ToHg19.over.chain",
                "--ld-wind-kb",
                "10",
                "--output-dir",
                "out/panel",
                "--keep-indivs-path",
                "samples.keep",
            ]
        )

        self.assertEqual(args.command, "build-ref-panel")
        self.assertEqual(args.plink_path, "data/reference/panel_chr@")
        self.assertEqual(args.output_dir, "out/panel")
        self.assertEqual(args.keep_indivs_path, "samples.keep")

    def test_build_ref_panel_help_fast_path_avoids_scipy_backed_imports(self):
        from ldsc import cli

        with mock.patch.object(cli, "_load_regression_runner", side_effect=AssertionError("regression import should not occur")), \
             mock.patch.object(cli, "_load_sumstats_munger", side_effect=AssertionError("munging import should not occur")):
            with self.assertRaises(SystemExit) as exc:
                cli.main(["build-ref-panel", "--help"])

        self.assertEqual(exc.exception.code, 0)

    def test_python_m_ldsc_does_not_treat_successful_results_as_exit_codes(self):
        sentinel = object()

        with mock.patch("ldsc.cli.main", return_value=sentinel) as mocked_main:
            with warnings.catch_warnings():
                warnings.filterwarnings(
                    "ignore",
                    message=r"'ldsc\.__main__' found in sys\.modules",
                    category=RuntimeWarning,
                )
                runpy.run_module("ldsc", run_name="__main__", alter_sys=True)

        mocked_main.assert_called_once_with()
