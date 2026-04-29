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
        self.assertTrue(hasattr(ldsc, "LDScoreDirectoryWriter"))
        self.assertTrue(hasattr(ldsc, "PartitionedH2OutputConfig"))
        self.assertTrue(hasattr(ldsc, "PartitionedH2DirectoryWriter"))
        self.assertTrue(hasattr(ldsc, "load_ldscore_from_dir"))
        self.assertTrue(hasattr(ldsc.__main__, "main"))
        removed_output_names = [
            "ArtifactConfig",
            "ArtifactOutputConfig",
            "ArtifactProducer",
            "OutputManager",
            "PostProcessor",
            "ResultFormatter",
            "ResultWriter",
        ]
        for name in removed_output_names:
            self.assertFalse(hasattr(ldsc, name), f"{name} should not be a public export")

    def test_top_level_exports_genome_build_inference_api(self):
        import ldsc
        from ldsc import ChrPosBuildInference, infer_chr_pos_build, resolve_chr_pos_table
        from ldsc.genome_build_inference import (
            ChrPosBuildInference as ModuleChrPosBuildInference,
            infer_chr_pos_build as module_infer_chr_pos_build,
            resolve_chr_pos_table as module_resolve_chr_pos_table,
        )

        self.assertIs(ChrPosBuildInference, ModuleChrPosBuildInference)
        self.assertIs(infer_chr_pos_build, module_infer_chr_pos_build)
        self.assertIs(resolve_chr_pos_table, module_resolve_chr_pos_table)
        self.assertIn("ChrPosBuildInference", ldsc.__all__)
        self.assertIn("infer_chr_pos_build", ldsc.__all__)
        self.assertIn("resolve_chr_pos_table", ldsc.__all__)

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
        self.assertNotIn("infer-build", subparsers_action.choices)

    def test_ldscore_subcommand_accepts_unified_path_flags(self):
        from ldsc import cli

        parser = cli.build_parser()
        args = parser.parse_args(
            [
                "ldscore",
                "--output-dir",
                "out/ldscores",
                "--baseline-annot-sources",
                "baseline.annot.gz",
                "--plink-prefix",
                "panel",
                "--ld-wind-snps",
                "10",
                "--keep-indivs-file",
                "samples.keep",
                "--maf-min",
                "0.01",
                "--common-maf-min",
                "0.2",
            ]
        )

        self.assertEqual(args.command, "ldscore")
        self.assertEqual(args.keep_indivs_file, "samples.keep")
        self.assertEqual(args.maf_min, 0.01)
        self.assertEqual(args.common_maf_min, 0.2)

    def test_all_output_subcommands_accept_overwrite_flag(self):
        from ldsc import cli

        parser = cli.build_parser()
        cases = [
            [
                "annotate",
                "--query-annot-bed-sources",
                "beds/query.bed",
                "--baseline-annot-sources",
                "baseline.1.annot.gz",
                "--output-dir",
                "out/annot",
                "--overwrite",
            ],
            [
                "ldscore",
                "--output-dir",
                "out/ldscores",
                "--baseline-annot-sources",
                "baseline.annot.gz",
                "--plink-prefix",
                "panel",
                "--ld-wind-snps",
                "10",
                "--overwrite",
            ],
            [
                "build-ref-panel",
                "--plink-prefix",
                "panel.@",
                "--source-genome-build",
                "hg38",
                "--output-dir",
                "out/panel",
                "--ld-wind-kb",
                "10",
                "--overwrite",
            ],
            [
                "munge-sumstats",
                "--sumstats-file",
                "raw.tsv",
                "--output-dir",
                "out/trait",
                "--overwrite",
            ],
            [
                "h2",
                "--sumstats-file",
                "trait.sumstats.gz",
                "--ldscore-dir",
                "ldscores",
                "--output-dir",
                "out/h2",
                "--overwrite",
            ],
            [
                "partitioned-h2",
                "--sumstats-file",
                "trait.sumstats.gz",
                "--ldscore-dir",
                "ldscores",
                "--output-dir",
                "out/partitioned",
                "--write-per-query-results",
                "--overwrite",
            ],
            [
                "rg",
                "--sumstats-1-file",
                "trait1.sumstats.gz",
                "--sumstats-2-file",
                "trait2.sumstats.gz",
                "--ldscore-dir",
                "ldscores",
                "--output-dir",
                "out/rg",
                "--overwrite",
            ],
        ]

        for argv in cases:
            with self.subTest(command=argv[0]):
                args = parser.parse_args(argv)
                self.assertTrue(args.overwrite)

    def test_ldscore_overwrite_flag_maps_to_output_config(self):
        from ldsc import cli
        from ldsc import ldscore_calculator

        parser = cli.build_parser()
        args = parser.parse_args(
            [
                "ldscore",
                "--output-dir",
                "out/ldscores",
                "--baseline-annot-sources",
                "baseline.annot.gz",
                "--plink-prefix",
                "panel",
                "--ld-wind-snps",
                "10",
                "--overwrite",
            ]
        )

        output_config = ldscore_calculator._output_config_from_args(args)

        self.assertEqual(output_config.output_dir, "out/ldscores")
        self.assertTrue(output_config.overwrite)

    def test_ldscore_subcommand_accepts_regression_snps_file(self):
        from ldsc import cli

        parser = cli.build_parser()
        args = parser.parse_args(
            [
                "ldscore",
                "--output-dir",
                "out/ldscores",
                "--baseline-annot-sources",
                "baseline.annot.gz",
                "--plink-prefix",
                "panel",
                "--ld-wind-snps",
                "10",
                "--regression-snps-file",
                "filters/hm3.txt",
            ]
        )

        self.assertEqual(args.command, "ldscore")
        self.assertEqual(args.regression_snps_file, "filters/hm3.txt")

    def test_ldscore_subcommand_rejects_removed_print_snps_and_regression_snps_flags(self):
        from ldsc import cli

        parser = cli.build_parser()
        with self.assertRaises(SystemExit):
            parser.parse_args(
                [
                    "ldscore",
                    "--output-dir",
                    "out/ldscores",
                    "--baseline-annot-sources",
                    "baseline.annot.gz",
                    "--plink-prefix",
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
                    "--query-annot-bed-sources",
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
                    "--plink-prefix",
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
                    "--baseline-annot-sources",
                    "baseline.annot.gz",
                    "--plink-prefix",
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
                "--plink-prefix",
                "data/reference/panel_chr@",
                "--source-genome-build",
                "hg38",
                "--genetic-map-hg19-sources",
                "maps/hg19.txt",
                "--genetic-map-hg38-sources",
                "maps/hg38.txt",
                "--liftover-chain-hg38-to-hg19-file",
                "chains/hg38ToHg19.over.chain",
                "--ld-wind-kb",
                "10",
                "--output-dir",
                "out/panel",
                "--keep-indivs-file",
                "samples.keep",
                "--maf-min",
                "0.01",
            ]
        )

        self.assertEqual(args.command, "build-ref-panel")
        self.assertEqual(args.plink_prefix, "data/reference/panel_chr@")
        self.assertEqual(args.output_dir, "out/panel")
        self.assertEqual(args.keep_indivs_file, "samples.keep")
        self.assertEqual(args.maf_min, 0.01)

    def test_old_maf_flag_is_rejected_for_reference_panel_workflows(self):
        from ldsc import cli

        parser = cli.build_parser()
        with self.assertRaises(SystemExit):
            parser.parse_args(
                [
                    "ldscore",
                    "--output-dir",
                    "out/ldscores",
                    "--baseline-annot-sources",
                    "baseline.annot.gz",
                    "--plink-prefix",
                    "panel",
                    "--ld-wind-snps",
                    "10",
                    "--maf",
                    "0.01",
                ]
            )
        with self.assertRaises(SystemExit):
            parser.parse_args(
                [
                    "build-ref-panel",
                    "--plink-prefix",
                    "data/reference/panel_chr@",
                    "--source-genome-build",
                    "hg38",
                    "--ld-wind-kb",
                    "10",
                    "--output-dir",
                    "out/panel",
                    "--maf",
                    "0.01",
                ]
            )

    def test_build_ref_panel_help_fast_path_avoids_scipy_backed_imports(self):
        from ldsc import cli

        with mock.patch.object(cli, "_load_regression_runner", side_effect=AssertionError("regression import should not occur")), \
             mock.patch.object(cli, "_load_sumstats_munger", side_effect=AssertionError("munging import should not occur")):
            with self.assertRaises(SystemExit) as exc:
                cli.main(["build-ref-panel", "--help"])

        self.assertEqual(exc.exception.code, 0)

    def test_python_m_ldsc_exits_zero_for_successful_results(self):

        with mock.patch("ldsc.cli.run_cli", return_value=0) as mocked_run_cli:
            with warnings.catch_warnings():
                warnings.filterwarnings(
                    "ignore",
                    message=r"'ldsc\.__main__' found in sys\.modules",
                    category=RuntimeWarning,
                )
                with self.assertRaises(SystemExit) as exc:
                    runpy.run_module("ldsc", run_name="__main__", alter_sys=True)

        self.assertEqual(exc.exception.code, 0)
        mocked_run_cli.assert_called_once_with()
