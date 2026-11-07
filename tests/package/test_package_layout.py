from pathlib import Path
import tempfile
import subprocess
import sys
import unittest

import pandas as pd


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

    def test_cli_rejects_a_separator_before_the_command(self):
        import contextlib
        import io
        from ldsc import cli

        stderr = io.StringIO()
        with contextlib.redirect_stderr(stderr):
            self.assertEqual(cli.run_cli(["--", "h2"]), 2)
        self.assertIn("usage: ldsc", stderr.getvalue())
        self.assertIn("put the command first", stderr.getvalue())

    def test_every_cli_command_names_functions_its_module_defines(self):
        from importlib import import_module
        from ldsc import cli

        for name, command in cli._COMMANDS.items():
            with self.subTest(command=name):
                self.assertEqual(command.add_arguments is None, command.run is None)
                module = import_module(command.module, "ldsc")
                functions = (
                    ("build_parser", "main") if command.add_arguments is None
                    else (command.add_arguments, command.run)
                )
                for function in functions:
                    self.assertTrue(callable(getattr(module, function, None)), f"{command.module}.{function}")

    def test_materializing_cli_subcommands_use_explicit_or_derived_output_directories(self):
        from ldsc import cli

        for command in cli._COMMANDS:
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

    def test_help_does_not_import_heavy_workflows(self):
        import subprocess
        import textwrap

        # `ldsc --help` only lists subcommands, and `build-r2-matrix` needs
        # neither regression nor munging, so neither help may import the
        # SciPy-heavy regression and munging workflows. Run in a clean
        # subprocess so an unrelated import in this process cannot mask a
        # regression.
        probe = textwrap.dedent(
            """
            import sys
            from ldsc import cli
            try:
                cli.main(sys.argv[1:])
            except SystemExit:
                pass
            heavy = [
                name
                for name in ("ldsc._regression.runner", "ldsc._regression.commands", "ldsc._sumstats.munger",
                             "scipy.stats._stats_py")
                if name in sys.modules
            ]
            print(",".join(heavy))
            """
        )
        for argv in (["--help"], ["build-r2-matrix", "--help"]):
            with self.subTest(argv=argv):
                completed = subprocess.run(
                    [sys.executable, "-c", probe, *argv],
                    capture_output=True,
                    text=True,
                    cwd=tempfile.gettempdir(),
                )
                self.assertEqual(completed.returncode, 0, completed.stderr)
                imported_heavy = completed.stdout.splitlines()[-1].strip()
                self.assertEqual(
                    imported_heavy,
                    "",
                    f"`ldsc {' '.join(argv)}` should not import heavy workflows, but imported: {imported_heavy}",
                )

    def test_regression_reads_munged_artifacts_without_the_munging_workflow(self):
        # Regression reads munged summary statistics through `_sumstats.table`
        # alone; the munging workflow (parser, input resolution, kernel,
        # liftover, writer) stays unloaded. Run in a clean subprocess so that
        # modules loaded by other tests cannot hide an import.
        probe = (
            "import importlib, pkgutil, sys\n"
            "import ldsc._regression as regression\n"
            "for module in pkgutil.walk_packages(regression.__path__, 'ldsc._regression.'):\n"
            "    importlib.import_module(module.name)\n"
            "print(','.join(sorted(name for name in sys.modules if name.startswith('ldsc._sumstats'))))\n"
        )
        completed = subprocess.run(
            [sys.executable, "-c", probe], capture_output=True, text=True, cwd=tempfile.gettempdir(),
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assertEqual(completed.stdout.splitlines()[-1], "ldsc._sumstats,ldsc._sumstats.table")

    def test_annotation_builder_owns_public_workflow_entrypoints(self):
        from ldsc._annotation import builder as annotation_builder

        for name in (
            "AnnotationBuilder",
            "AnnotationBundle",
            "AnnotationResult",
            "run_annotate",
            "run_annotate_from_args",
            "build_parser",
            "main",
        ):
            self.assertTrue(hasattr(annotation_builder, name), f"{name} should be public")

