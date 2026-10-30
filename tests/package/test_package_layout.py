from pathlib import Path
import tempfile
import subprocess
import sys
import unittest
from unittest import mock

import pandas as pd

from tests.cli_parsing import parse_command


class PackageLayoutTest(unittest.TestCase):
    def test_installed_import_cli_and_resources_outside_repository(self):
        probe = """
import sys
from importlib.resources import files
import ldsc
assert 'matplotlib.pyplot' not in sys.modules
resources = files('ldsc').joinpath('data')
with resources.joinpath('hm3_curated_map.tsv.gz').open('rb') as stream:
    assert stream.read(2) == b'\\x1f\\x8b'
assert list(resources.joinpath('regions').iterdir())
print(ldsc.__file__)
"""
        with tempfile.TemporaryDirectory() as outside:
            for command in (
                [sys.executable, "-I", "-c", probe],
                [sys.executable, "-I", "-m", "ldsc", "--help"],
                [str(Path(sys.executable).with_name("ldsc")), "--help"],
            ):
                completed = subprocess.run(command, cwd=outside, capture_output=True, text=True)
                self.assertEqual(completed.returncode, 0, completed.stderr)

    def test_public_package_imports(self):
        import ldsc
        import ldsc.__main__

        self.assertTrue(hasattr(ldsc, "LDScoreCalculator"))
        self.assertTrue(hasattr(ldsc, "RegressionRunner"))
        self.assertTrue(hasattr(ldsc, "SumstatsMunger"))
        self.assertTrue(hasattr(ldsc, "AnnotationBuilder"))
        self.assertTrue(hasattr(ldsc, "R2MatrixBuilder"))
        self.assertTrue(hasattr(ldsc, "H2OutputConfig"))
        self.assertTrue(hasattr(ldsc, "H2DirectoryWriter"))
        self.assertTrue(hasattr(ldsc, "LDScoreOutputConfig"))
        self.assertTrue(hasattr(ldsc, "LDScoreDirectoryWriter"))
        self.assertTrue(hasattr(ldsc, "PartitionedH2OutputConfig"))
        self.assertTrue(hasattr(ldsc, "PartitionedH2DirectoryWriter"))
        self.assertTrue(hasattr(ldsc, "load_ldscore_from_dir"))
        self.assertTrue(hasattr(ldsc.__main__, "main"))

    def test_top_level_exports_genome_build_inference_api(self):
        import ldsc
        from ldsc import ChrPosBuildInference, infer_chr_pos_build, resolve_chr_pos_table
        from ldsc._genome.build_inference import (
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

    def test_top_level_exports_sumstats_inference_api(self):
        import ldsc
        from ldsc import RawSumstatsInference, infer_raw_sumstats
        from ldsc._sumstats.munger import (
            RawSumstatsInference as ModuleRawSumstatsInference,
            infer_raw_sumstats as module_infer_raw_sumstats,
        )

        self.assertIs(RawSumstatsInference, ModuleRawSumstatsInference)
        self.assertIs(infer_raw_sumstats, module_infer_raw_sumstats)
        self.assertIn("RawSumstatsInference", ldsc.__all__)
        self.assertIn("infer_raw_sumstats", ldsc.__all__)

    def test_top_level_exports_hm3_curated_map_loader(self):
        import ldsc
        from ldsc import load_hm3_curated_map

        self.assertIs(load_hm3_curated_map, ldsc.load_hm3_curated_map)
        self.assertIn("load_hm3_curated_map", ldsc.__all__)

    def test_load_hm3_curated_map_returns_canonical_full_schema(self):
        from ldsc import load_hm3_curated_map

        frame = load_hm3_curated_map()

        self.assertIn("CHR", frame.columns)
        self.assertIn("hg19_POS", frame.columns)
        self.assertIn("hg38_POS", frame.columns)
        self.assertIn("SNP", frame.columns)
        self.assertIn("A1", frame.columns)
        self.assertIn("A2", frame.columns)
        self.assertGreater(len(frame), 1_000_000)
        self.assertEqual(str(frame["CHR"].dtype), "string")
        self.assertTrue(pd.api.types.is_integer_dtype(frame["hg19_POS"]))
        self.assertTrue(pd.api.types.is_integer_dtype(frame["hg38_POS"]))
        self.assertEqual(str(frame["SNP"].dtype), "string")
        self.assertEqual(str(frame["A1"].dtype), "string")
        self.assertEqual(str(frame["A2"].dtype), "string")
        observed_alleles = set(frame["A1"].dropna().unique()) | set(frame["A2"].dropna().unique())
        self.assertLessEqual(observed_alleles, {"A", "C", "G", "T"})

    def test_cli_exposes_expected_subcommands(self):
        from ldsc import cli

        parser = cli._build_top_level_parser()
        subparsers_action = next(
            action for action in parser._actions if action.__class__.__name__ == "_SubParsersAction"
        )
        self.assertEqual(
            set(subparsers_action.choices),
            {
                "annotate", "ldscore", "build-r2-matrix", "build-gene-ldscore-index",
                "convert-legacy-ldscores",
                "munge-sumstats", "h2", "partitioned-h2", "quantile-h2", "rg", "query-r2-matrix",
                "convert-h2-scale", "plot", "infer-genome-build",
            },
        )

    def test_materializing_cli_subcommands_use_explicit_or_derived_output_directories(self):
        from ldsc import cli

        for command in cli._SUBCOMMAND_HELP:
            command_parser = cli._command_parser(command)
            with self.subTest(command=command):
                output_actions = [
                    action for action in command_parser._actions if action.dest == "output_dir"
                ]
                if command in {"plot", "convert-h2-scale", "infer-genome-build"}:
                    self.assertEqual(output_actions, [])
                else:
                    self.assertEqual(len(output_actions), 1)
                    self.assertTrue(output_actions[0].required)

    def test_top_level_help_does_not_import_heavy_workflows(self):
        import subprocess
        import textwrap

        # `ldsc --help` only lists subcommands, so it must not import the
        # SciPy-heavy regression and munging workflows. Run in a clean
        # subprocess so an unrelated import in this process cannot mask a
        # regression.
        probe = textwrap.dedent(
            """
            import sys
            from ldsc import cli
            try:
                cli.main(["--help"])
            except SystemExit:
                pass
            heavy = [
                name
                for name in ("ldsc._regression.runner", "ldsc._sumstats.munger", "scipy.stats._stats_py")
                if name in sys.modules
            ]
            print(",".join(heavy))
            """
        )
        completed = subprocess.run(
            [sys.executable, "-c", probe],
            capture_output=True,
            text=True,
            cwd=tempfile.gettempdir(),
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)
        imported_heavy = completed.stdout.splitlines()[-1].strip()
        self.assertEqual(
            imported_heavy,
            "",
            f"`ldsc --help` should not import heavy workflows, but imported: {imported_heavy}",
        )

    def test_annotation_builder_owns_public_workflow_entrypoints(self):
        from ldsc._annotation import builder as annotation_builder

        for name in (
            "AnnotationBuilder",
            "AnnotationBundle",
            "run_annotate",
            "run_annotate_from_args",
            "add_annotate_arguments",
            "build_parser",
            "parse_annotate_args",
            "main",
        ):
            self.assertTrue(hasattr(annotation_builder, name), f"{name} should be public")

    def test_ldscore_subcommand_accepts_unified_path_flags(self):
        args = parse_command(
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
        cases = [
            [
                "annotate",
                "--query-annot-bed-sources",
                "beds/query.bed",
                "--ref-snp-annot-sources",
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
                "build-r2-matrix",
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
                "--raw-sumstats-file",
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
                "--overwrite",
            ],
            [
                "rg",
                "--sumstats-sources",
                "trait1.sumstats.gz",
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
                args = parse_command(argv)
                self.assertTrue(args.overwrite)

    def test_ldscore_overwrite_flag_maps_to_output_config(self):
        from ldsc._ldscore import calculator as ldscore_calculator

        args = parse_command(
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

    def test_ldscore_subcommand_accepts_regr_snps_file_and_rejects_old_name(self):
        args = parse_command(
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
                "--regr-snps-file",
                "filters/hm3.txt",
            ]
        )

        self.assertEqual(args.command, "ldscore")
        self.assertEqual(args.regr_snps_file, "filters/hm3.txt")

    def test_munge_sumstats_subcommand_accepts_sumstats_snps_file(self):
        args = parse_command(
            [
                "munge-sumstats",
                "--raw-sumstats-file",
                "raw.tsv",
                "--output-dir",
                "out/trait",
                "--sumstats-snps-file",
                "filters/hm3.tsv.gz",
            ]
        )

        self.assertEqual(args.command, "munge-sumstats")
        self.assertEqual(args.sumstats_snps_file, "filters/hm3.tsv.gz")

    def test_munge_sumstats_subcommand_accepts_no_snp_restriction(self):
        args = parse_command(
            [
                "munge-sumstats",
                "--raw-sumstats-file",
                "raw.tsv",
                "--output-dir",
                "out/trait",
                "--no-snp-restriction",
            ]
        )

        self.assertEqual(args.command, "munge-sumstats")
        self.assertTrue(args.no_snp_restriction)

    def test_ldscore_subcommand_rejects_legacy_print_snps_flag(self):
        with self.assertRaises(SystemExit):
            parse_command(
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

    def test_ldscore_subcommand_rejects_removed_out_flag(self):
        with self.assertRaises(SystemExit) as exc:
            parse_command(
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

    def test_build_ref_panel_subcommand_accepts_output_dir_and_no_panel_label(self):
        args = parse_command(
            [
                "build-r2-matrix",
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

        self.assertEqual(args.command, "build-r2-matrix")
        self.assertEqual(args.plink_prefix, "data/reference/panel_chr@")
        self.assertEqual(args.output_dir, "out/panel")
        self.assertEqual(args.keep_indivs_file, "samples.keep")
        self.assertEqual(args.maf_min, 0.01)

    def test_old_maf_flag_is_rejected_for_reference_panel_workflows(self):
        with self.assertRaises(SystemExit):
            parse_command(
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
            parse_command(
                [
                    "build-r2-matrix",
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
                cli.main(["build-r2-matrix", "--help"])

        self.assertEqual(exc.exception.code, 0)
