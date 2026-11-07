from pathlib import Path
import tempfile
import subprocess
import sys
import unittest

import pandas as pd
import pytest


# Names users import from ``ldsc``; each must be listed in ``ldsc.__all__``.
PUBLIC_NAMES = [
    "LDScoreCalculator", "RegressionRunner", "SumstatsMunger", "AnnotationBuilder", "AnnotationBundle",
    "AnnotationResult", "run_annotate", "R2MatrixBuilder", "H2OutputConfig", "H2DirectoryWriter",
    "LDScoreOutputConfig", "LDScoreDirectoryWriter", "PartitionedH2OutputConfig", "PartitionedH2DirectoryWriter",
    "load_ldscore_from_dir", "ChrPosBuildInference", "infer_chr_pos_build", "resolve_chr_pos_table",
    "RawSumstatsInference", "infer_raw_sumstats", "load_hm3_curated_map", "load_sumstats", "R2Matrix",
    "query_r2_matrix", "unbiased_r2_to_pearson_r",
]


@pytest.mark.parametrize("name", PUBLIC_NAMES)
def test_public_name_is_exported(name):
    import ldsc

    assert name in ldsc.__all__
    assert getattr(ldsc, name) is not None


def test_module_entry_point_exposes_main():
    import ldsc.__main__

    assert callable(ldsc.__main__.main)


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

