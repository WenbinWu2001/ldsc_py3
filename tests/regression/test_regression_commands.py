"""The h2, partitioned-h2 and rg commands: parsers, argument handling, output preflight, and run routes."""

import argparse
import contextlib
from dataclasses import replace
import gzip
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest import mock

import numpy as np
import pandas as pd
import pytest

from ldsc.config import GlobalConfig, set_global_config
from ldsc.errors import LDSCConfigError, LDSCUsageError

from ldsc import cli
from ldsc._regression import commands, partitioned, rg
from ldsc._regression.outputs import LEGACY_SUMSTATS_DROP_COLUMNS, RgDirectoryWriter
from ldsc._regression.results import RgResultFamily
from tests.regression_fixtures import RegressionFixtures, parsed_regression_args, stub_h2_kernel


class RegressionCommandsTest(RegressionFixtures, unittest.TestCase):
    def test_run_h2_from_args_uses_ldscore_dir(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            set_global_config(GlobalConfig(snp_identifier="rsid"))
            with gzip.open(tmpdir / "trait.sumstats.gz", "wt", encoding="utf-8") as handle:
                handle.write("SNP\tA1\tA2\tZ\tN\nrs1\tA\tC\t1.0\t1000\n")
            self.write_sumstats_sidecar(tmpdir / "metadata.json", trait_name="trait")
            ldscore_dir = self.write_ldscore_dir(tmpdir / "ldscores", include_query=True)
            args = parsed_regression_args(
                "h2",
                {
                    "sumstats_file": str(tmpdir / "trait.sumstats.gz"),
                    "trait_name": "trait",
                    "ldscore_dir": str(ldscore_dir),
                    "count_kind": "common",
                    "output_dir": "$LDSC_REGRESSION_OUT/h2_out",
                    "overwrite": False,
                    "log_level": "INFO",
                    "num_blocks": 200,
                    "intercept_h2": None,
                    "two_step_cutoff": None,
                    "chisq_max": None,
                    "log_level": "INFO",
                },
            )

            marker = tmpdir / "h2_out/RUN_FAILED.txt"
            marker.parent.mkdir()
            marker.write_text("earlier failed overwrite")
            with mock.patch.dict("os.environ", {"LDSC_REGRESSION_OUT": str(tmpdir)}), contextlib.chdir(tmpdir), stub_h2_kernel(
                mock.Mock(
                    tot=np.array([0.1]),
                    tot_se=np.array([0.01]),
                    intercept=np.array([1.0]),
                    intercept_se=0.01,
                    mean_chisq=np.array([1.1]),
                    lambda_gc=np.array([1.0]),
                    ratio=0.0,
                    ratio_se=0.0,
                    coef=np.array([0.001]),
                ),
            ) as patched:
                result = commands.run_h2_from_args(args)

            summary = result.summary
            self.assertEqual(result.output_paths["summary"], str(tmpdir / "h2_out" / "h2.tsv"))
            self.assertEqual(result.legacy_sumstats_drops.columns.tolist(), LEGACY_SUMSTATS_DROP_COLUMNS)
            patched.assert_called_once()
            self.assertFalse(marker.exists())
            self.assertFalse((tmpdir / "$LDSC_REGRESSION_OUT").exists())
            self.assertEqual(patched.call_args.args[1].shape[1], 1)
            self.assertEqual(summary.loc[0, "trait_name"], "trait")
            self.assertTrue((tmpdir / "h2_out" / "h2.tsv").exists())
            self.assertFalse((tmpdir / "h2_out" / "metadata.json").exists())
            metadata = json.loads((tmpdir / "h2_out" / "diagnostics" / "metadata.json").read_text(encoding="utf-8"))
            self.assertEqual(
                metadata,
                {
                    "artifact_type": "h2_result",
                    "files": {
                        "summary": "h2.tsv",
                        "ld_score_regression_bins": "diagnostics/ld_score_regression_bins.tsv",
                        "legacy_sumstats_drops": "diagnostics/dropped_snps/legacy_sumstats.tsv.gz",
                    },
                    "trait_name": "trait",
                    "sumstats_file": str(tmpdir / "trait.sumstats.gz"),
                    "ldscore_dir": str(ldscore_dir.resolve()),
                    "effective_snp_identifier": "rsid",
                    "genome_build": None,
                    "identity_downgrade_applied": False,
                    "count_key_used_for_regression": "common_reference_snp_counts",
                    "retained_ld_columns": ["base"],
                    "dropped_zero_variance_ld_columns": [],
                    "dropped_nonfinite_statistic_rows": 0,
                    "n_snps": 1,
                    "n_blocks_used": 1,
                    "effective_chisq_max": None,
                    "effective_two_step_cutoff": 30,
                    "intercept_h2_policy": "free",
                    "samp_prev": None,
                    "pop_prev": None,
                    "scale": "observed",
                },
            )
            log_text = (tmpdir / "h2_out" / "diagnostics" / "h2.log").read_text(encoding="utf-8")
            # The command is the entry point and logs the GlobalConfig once; dataset assembly logs none.
            self.assertEqual(log_text.count("using GlobalConfig("), 1)
            self.assertIn("run_h2_from_args using GlobalConfig(", log_text)
            self.assertTrue((tmpdir / "h2_out" / "diagnostics" / "ld_score_regression_bins.tsv").exists())
            legacy_audit = tmpdir / "h2_out" / "diagnostics" / "dropped_snps" / "legacy_sumstats.tsv.gz"
            self.assertTrue(legacy_audit.exists())
            self.assertEqual(
                pd.read_csv(legacy_audit, sep="\t", compression="gzip").columns.tolist(),
                ["trait_name", "source_path", "SNP", "A1", "A2", "reason", "panel_candidate_count"],
            )

    def _h2_args(self, tmpdir, ldscore_dir, *, samp_prev=None, pop_prev=None):
        return parsed_regression_args(
            "h2",
            {
                "sumstats_file": str(tmpdir / "trait.sumstats.gz"),
                "trait_name": "trait",
                "ldscore_dir": str(ldscore_dir),
                "count_kind": "common",
                "output_dir": str(tmpdir / "h2_out"),
                "overwrite": False,
                "log_level": "INFO",
                "num_blocks": 200,
                "intercept_h2": None,
                "two_step_cutoff": None,
                "chisq_max": None,
                "samp_prev": samp_prev,
                "pop_prev": pop_prev,
            },
        )

    def _write_h2_inputs(self, tmpdir):
        set_global_config(GlobalConfig(snp_identifier="rsid"))
        with gzip.open(tmpdir / "trait.sumstats.gz", "wt", encoding="utf-8") as handle:
            handle.write("SNP\tA1\tA2\tZ\tN\nrs1\tA\tC\t1.0\t1000\n")
        self.write_sumstats_sidecar(tmpdir / "metadata.json", trait_name="trait")
        return self.write_ldscore_dir(tmpdir / "ldscores", include_query=True)

    def _patched_h2(self):
        return stub_h2_kernel(
            mock.Mock(
                tot=np.array([0.2]),
                tot_se=np.array([0.03]),
                intercept=np.array([1.0]),
                intercept_se=0.01,
                mean_chisq=np.array([1.1]),
                lambda_gc=np.array([1.0]),
                ratio=0.0,
                ratio_se=0.0,
                coef=np.array([0.001]),
            ),
        )

    def test_run_h2_from_args_liability_scale(self):
        from ldsc._regression.estimators import liability_conversion_factor

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            ldscore_dir = self._write_h2_inputs(tmpdir)
            args = self._h2_args(tmpdir, ldscore_dir, samp_prev=0.5, pop_prev=0.01)
            with self._patched_h2():
                summary = commands.run_h2_from_args(args).summary
        c = liability_conversion_factor(0.5, 0.01)
        row = summary.iloc[0]
        self.assertAlmostEqual(row["total_h2_liab"], 0.2 * c, places=12)
        self.assertEqual(row["total_h2_obs"], 0.2)
        self.assertEqual(row["samp_prev"], 0.5)
        self.assertEqual(row["pop_prev"], 0.01)

    def test_run_h2_from_args_half_specified_prevalence_raises(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            ldscore_dir = self._write_h2_inputs(tmpdir)
            args = self._h2_args(tmpdir, ldscore_dir, samp_prev=0.5, pop_prev=None)
            with self._patched_h2(), self.assertRaisesRegex(LDSCConfigError, "exactly one of sample/population"):
                commands.run_h2_from_args(args)

    def test_run_h2_from_args_requires_output_dir_before_loading_inputs(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            set_global_config(GlobalConfig(snp_identifier="rsid"))
            with gzip.open(tmpdir / "trait.sumstats.gz", "wt", encoding="utf-8") as handle:
                handle.write("SNP\tA1\tA2\tZ\tN\nrs1\tA\tC\t1.0\t1000\n")
            self.write_sumstats_sidecar(tmpdir / "metadata.json", trait_name="trait")
            ldscore_dir = self.write_ldscore_dir(tmpdir / "ldscores", include_query=True)
            args = parsed_regression_args(
                "h2",
                {
                    "sumstats_file": str(tmpdir / "trait.sumstats.gz"),
                    "trait_name": "trait",
                    "ldscore_dir": str(ldscore_dir),
                    "count_kind": "common",
                    "output_dir": None,
                    "overwrite": False,
                    "log_level": "INFO",
                    "num_blocks": 200,
                    "intercept_h2": None,
                    "two_step_cutoff": None,
                    "chisq_max": None,
                },
            )

            with mock.patch.object(
                commands,
                "_load_sumstats_table",
                side_effect=AssertionError("sumstats should not load"),
            ):
                with self.assertRaisesRegex(LDSCUsageError, "requires `--output-dir`"):
                    commands.run_h2_from_args(args)

            self.assertFalse(list(tmpdir.glob("*.log")))
            self.assertFalse((tmpdir / "h2_out" / "metadata.json").exists())
            self.assertFalse((tmpdir / "h2_out" / "diagnostics").exists())

    def test_run_h2_from_args_labels_unnamed_trait_by_filename(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            set_global_config(GlobalConfig(snp_identifier="rsid"))
            self.write_footer_sumstats_parquet(tmpdir / "trait.parquet", trait_name=None, snp_identifier="rsid")
            ldscore_dir = self.write_ldscore_dir(tmpdir / "ldscores", include_query=True)
            args = parsed_regression_args(
                "h2",
                {
                    "sumstats_file": str(tmpdir / "trait.parquet"),
                    "trait_name": None,
                    "ldscore_dir": str(ldscore_dir),
                    "count_kind": "common",
                    "output_dir": str(tmpdir / "h2_out"),
                    "overwrite": False,
                    "log_level": "INFO",
                    "num_blocks": 200,
                    "intercept_h2": None,
                    "two_step_cutoff": None,
                    "chisq_max": None,
                },
            )

            with stub_h2_kernel(
                mock.Mock(
                    tot=np.array([0.1]),
                    tot_se=np.array([0.01]),
                    intercept=np.array([1.0]),
                    intercept_se=0.01,
                    mean_chisq=np.array([1.1]),
                    lambda_gc=np.array([1.0]),
                    ratio=0.0,
                    ratio_se=0.0,
                    coef=np.array([0.001]),
                ),
            ):
                summary = commands.run_h2_from_args(args).summary

        self.assertEqual(summary.loc[0, "trait_name"], "trait.parquet")

    def test_common_regression_arguments_expose_only_ldscore_dir(self):
        parser = argparse.ArgumentParser()
        commands.add_h2_arguments(parser)

        args = parser.parse_args(
            ["--ldscore-dir", "ldscores", "--sumstats-file", "trait.sumstats.gz", "--output-dir", "out"]
        )

        self.assertEqual(args.ldscore_dir, "ldscores")
        self.assertEqual(args.sumstats_file, "trait.sumstats.gz")
        self.assertEqual(args.count_kind, "common")
        self.assertEqual(args.log_level, "INFO")
        parsed = parser.parse_args(
            [
                "--ldscore-dir",
                "ldscores",
                "--sumstats-file",
                "trait.sumstats.gz",
                "--output-dir",
                "out",
                "--count-kind",
                "all",
                "--log-level",
                "DEBUG",
            ]
        )
        self.assertEqual(parsed.count_kind, "all")
        self.assertEqual(parsed.log_level, "DEBUG")
        with self.assertRaises(SystemExit):
            parser.parse_args(["--ldscore", "x", "--sumstats-file", "trait.sumstats.gz", "--output-dir", "out"])
        with self.assertRaises(SystemExit):
            parser.parse_args(["--ldscore-dir", "ldscores", "--sumstats", "trait.sumstats.gz", "--output-dir", "out"])
        with self.assertRaises(SystemExit):
            parser.parse_args(["--ldscore-dir", "ldscores", "--sumstats-file", "trait.sumstats.gz", "--output-dir", "out",
                               "--count-kind", "m_5_50"])

    def test_partitioned_h2_arguments_accept_summary_sort_by(self):
        parser = argparse.ArgumentParser()
        commands.add_partitioned_h2_arguments(parser)

        default_args = parser.parse_args(
            [
                "--ldscore-dir",
                "ldscores",
                "--sumstats-file",
                "trait.sumstats.gz",
                "--output-dir",
                "out",
            ]
        )
        sorted_args = parser.parse_args(
            [
                "--ldscore-dir",
                "ldscores",
                "--sumstats-file",
                "trait.sumstats.gz",
                "--output-dir",
                "out",
                "--summary-sort-by",
                "enrichment-p",
            ]
        )

        self.assertEqual(default_args.summary_sort_by, "auto")
        self.assertEqual(sorted_args.summary_sort_by, "enrichment-p")
        with self.assertRaises(SystemExit):
            parser.parse_args(
                [
                    "--ldscore-dir",
                    "ldscores",
                    "--sumstats-file",
                    "trait.sumstats.gz",
                    "--output-dir",
                    "out",
                    "--summary-sort-by",
                    "prop_h2",
                ]
            )

    def test_rg_arguments_use_sumstats_sources_and_reject_old_pairwise_flags(self):
        parser = argparse.ArgumentParser()
        commands.add_rg_arguments(parser)

        args = parser.parse_args(
            [
                "--ldscore-dir",
                "ldscores",
                "--sumstats-sources",
                "a.sumstats.gz",
                "b.sumstats.gz",
                "--output-dir",
                "out",
                "--intercept-h2",
                "1.02",
                "--intercept-gencov",
                "0.01",
            ]
        )

        self.assertEqual(args.sumstats_sources, ["a.sumstats.gz", "b.sumstats.gz"])
        self.assertEqual(args.intercept_h2, 1.02)
        self.assertEqual(args.intercept_gencov, 0.01)
        anchored = parser.parse_args(
            [
                "--ldscore-dir",
                "ldscores",
                "--sumstats-sources",
                "a.sumstats.gz",
                "b.sumstats.gz",
                "--output-dir",
                "out",
                "--anchor-trait",
                "MDD",
            ]
        )
        self.assertEqual(anchored.anchor_trait, "MDD")
        with self.assertRaises(SystemExit):
            parser.parse_args(
                [
                    "--ldscore-dir",
                    "ldscores",
                    "--sumstats-sources",
                    "a.sumstats.gz",
                    "b.sumstats.gz",
                    "--output-dir",
                    "out",
                    "--sumstats-1-file",
                    "a.sumstats.gz",
                ]
            )
        with self.assertRaises(SystemExit):
            parser.parse_args(
                [
                    "--ldscore-dir",
                    "ldscores",
                    "--sumstats-sources",
                    "a.sumstats.gz",
                    "b.sumstats.gz",
                    "--output-dir",
                    "out",
                    "--anchor-trait-file",
                    "a.sumstats.gz",
                ]
            )

    def test_run_rg_from_args_rejects_missing_output_dir_before_loading_inputs(self):
        args = parsed_regression_args(
            "rg",
            {
                "sumstats_sources": ["a.sumstats.gz", "b.sumstats.gz"],
                "anchor_trait": None,
                "ldscore_dir": "ldscores",
                "count_kind": "common",
                "output_dir": None,
                "overwrite": False,
                "write_per_pair_detail": True,
                "num_blocks": 200,
                "intercept_h2": None,
                "intercept_gencov": None,
                "two_step_cutoff": None,
                "chisq_max": None,
                "log_level": "INFO",
            },
        )

        with mock.patch.object(
            commands,
            "_load_sumstats_table",
            side_effect=AssertionError("sumstats should not load"),
        ):
            with self.assertRaisesRegex(LDSCUsageError, "requires `--output-dir`"):
                commands.run_rg_from_args(args)

    def test_run_rg_from_args_returns_result_family_and_writes_output(self):
        from ldsc._regression.outputs import RG_FULL_COLUMNS
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            set_global_config(GlobalConfig(snp_identifier="rsid"))
            sumstats_sources = []
            for name in ("a", "b"):
                trait_dir = tmpdir / name
                trait_dir.mkdir()
                source = trait_dir / "sumstats.sumstats.gz"
                with gzip.open(source, "wt", encoding="utf-8") as handle:
                    handle.write("SNP\tA1\tA2\tZ\tN\nrs1\tA\tC\t1.0\t1000\n")
                self.write_sumstats_sidecar(trait_dir / "metadata.json", trait_name=name)
                sumstats_sources.append(str(source))
            ldscore_dir = self.write_ldscore_dir(tmpdir / "ldscores", include_query=True)
            expected = RgResultFamily(
                rg=pd.DataFrame([{"trait_1": "a", "trait_2": "b", "n_snps_used": 1, "rg": 0.1, "rg_se": 0.01, "p": 0.02, "note": ""}]),
                rg_full=pd.DataFrame([{"trait_1": "a", "trait_2": "b", "n_snps_used": 1, "rg": 0.1, "rg_se": 0.01, "z": 10.0, "p": 0.02, "status": "ok", "error": ""}]),
                h2_per_trait=pd.DataFrame([{"trait_name": "a"}, {"trait_name": "b"}]),
                per_pair_metadata=[],
                legacy_sumstats_drops=pd.DataFrame(columns=LEGACY_SUMSTATS_DROP_COLUMNS),
            )
            args = parsed_regression_args(
                "rg",
                {
                    "sumstats_sources": sumstats_sources,
                    "anchor_trait": None,
                    "ldscore_dir": str(ldscore_dir),
                    "count_kind": "common",
                    "output_dir": "$LDSC_REGRESSION_OUT/out",
                    "overwrite": False,
                    "write_per_pair_detail": False,
                    "num_blocks": 200,
                    "intercept_h2": None,
                    "intercept_gencov": None,
                    "two_step_cutoff": None,
                    "chisq_max": None,
                    "log_level": "INFO",
                },
            )

            expected = replace(expected, rg_full=expected.rg_full.reindex(columns=RG_FULL_COLUMNS))
            marker = tmpdir / "out/RUN_FAILED.txt"
            marker.parent.mkdir()
            marker.write_text("earlier failed overwrite")
            stdout = io.StringIO()
            with mock.patch.dict("os.environ", {"LDSC_REGRESSION_OUT": str(tmpdir)}), contextlib.chdir(tmpdir), mock.patch.object(rg, "estimate_rg_pairs", return_value=expected), mock.patch.object(
                RgDirectoryWriter,
                "write",
                wraps=RgDirectoryWriter().write,
            ) as writer, contextlib.redirect_stdout(stdout):
                result = commands.run_rg_from_args(args)
            self.assertFalse(marker.exists())
            self.assertFalse((tmpdir / "$LDSC_REGRESSION_OUT").exists())
            self.assertTrue((tmpdir / "out/rg.tsv").is_file())
            self.assertTrue((tmpdir / "out/diagnostics/dropped_snps/legacy_sumstats.tsv.gz").is_file())

        # The returned family is the estimated one plus the paths the writer reports.
        self.assertIs(result.rg_full, expected.rg_full)
        self.assertEqual(result.output_paths["rg"], str(tmpdir / "out/rg.tsv"))
        self.assertEqual(result.output_paths["legacy_sumstats_drops"],
                         str(tmpdir / "out/diagnostics/dropped_snps/legacy_sumstats.tsv.gz"))
        writer.assert_called_once()
        self.assertEqual(stdout.getvalue(), "")

    def test_run_rg_from_args_recovers_metadata_labels_and_resolves_anchor_by_trait_name(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            set_global_config(GlobalConfig(snp_identifier="rsid"))
            sumstats_sources = []
            for stem, trait_name in (("a", "MDD"), ("b", "SCZ")):
                trait_dir = tmpdir / stem
                trait_dir.mkdir()
                source = trait_dir / "sumstats.parquet"
                self.write_footer_sumstats_parquet(source, trait_name=trait_name, snp_identifier="rsid")
                sumstats_sources.append(str(source))
            ldscore_dir = self.write_ldscore_dir(tmpdir / "ldscores", include_query=True)
            expected = RgResultFamily(
                rg=pd.DataFrame(),
                rg_full=pd.DataFrame(),
                h2_per_trait=pd.DataFrame(),
                per_pair_metadata=[],
            )
            args = parsed_regression_args(
                "rg",
                {
                    "sumstats_sources": sumstats_sources,
                    "anchor_trait": "SCZ",
                    "ldscore_dir": str(ldscore_dir),
                    "count_kind": "common",
                    "output_dir": str(tmpdir / "out"),
                    "overwrite": False,
                    "write_per_pair_detail": False,
                    "num_blocks": 200,
                    "intercept_h2": None,
                    "intercept_gencov": None,
                    "two_step_cutoff": None,
                    "chisq_max": None,
                    "log_level": "INFO",
                },
            )

            with mock.patch.object(rg, "estimate_rg_pairs", return_value=expected) as patched:
                result = commands.run_rg_from_args(args)

        self.assertIs(result.rg, expected.rg)
        tables = patched.call_args.args[0]
        self.assertEqual([table.trait_name for table in tables], ["MDD", "SCZ"])
        self.assertEqual(patched.call_args.kwargs["anchor_index"], 1)

    def test_run_rg_from_args_resolves_anchor_by_path_when_no_trait_name_matches(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            set_global_config(GlobalConfig(snp_identifier="rsid"))
            sumstats_sources = []
            for name in ("a", "b"):
                trait_dir = tmpdir / name
                trait_dir.mkdir()
                source = trait_dir / "sumstats.sumstats.gz"
                with gzip.open(source, "wt", encoding="utf-8") as handle:
                    handle.write("SNP\tA1\tA2\tZ\tN\nrs1\tA\tC\t1.0\t1000\n")
                self.write_sumstats_sidecar(trait_dir / "metadata.json", trait_name=name)
                sumstats_sources.append(str(source))
            ldscore_dir = self.write_ldscore_dir(tmpdir / "ldscores", include_query=True)
            expected = RgResultFamily(
                rg=pd.DataFrame(),
                rg_full=pd.DataFrame(),
                h2_per_trait=pd.DataFrame(),
                per_pair_metadata=[],
                legacy_sumstats_drops=pd.DataFrame(columns=LEGACY_SUMSTATS_DROP_COLUMNS),
            )
            args = parsed_regression_args(
                "rg",
                {
                    "sumstats_sources": sumstats_sources,
                    "anchor_trait": sumstats_sources[1],
                    "ldscore_dir": str(ldscore_dir),
                    "count_kind": "common",
                    "output_dir": str(tmpdir / "out"),
                    "overwrite": False,
                    "write_per_pair_detail": False,
                    "num_blocks": 200,
                    "intercept_h2": None,
                    "intercept_gencov": None,
                    "two_step_cutoff": None,
                    "chisq_max": None,
                    "log_level": "INFO",
                },
            )

            with mock.patch.object(rg, "estimate_rg_pairs", return_value=expected) as patched:
                commands.run_rg_from_args(args)

        self.assertEqual(patched.call_args.kwargs["anchor_index"], 1)

    def test_cli_requires_output_dir_for_rg(self):
        expected = type(
            "Result",
            (),
            {
                "rg": pd.DataFrame(
                    [
                        {
                            "trait_1": "a",
                            "trait_2": "b",
                            "n_snps_used": 0,
                            "rg": np.nan,
                            "rg_se": np.nan,
                            "p": np.nan,
                            "note": "Failed",
                        }
                    ]
                )
            },
        )()

        stdout = io.StringIO()
        with mock.patch.object(commands, "run_rg_from_args", return_value=expected) as patched, contextlib.redirect_stdout(stdout):
            with self.assertRaises(SystemExit):
                cli.main(
                [
                    "rg",
                    "--ldscore-dir",
                    "ldscores",
                    "--sumstats-sources",
                    "a.sumstats.gz",
                    "b.sumstats.gz",
                ]
            )

        patched.assert_not_called()
        self.assertEqual(stdout.getvalue(), "")

    def test_cli_does_not_print_rg_table_when_output_dir_is_supplied(self):
        expected = type("Result", (), {"rg": pd.DataFrame([{"trait_1": "a"}]),
                                       "dropped_orientation_unknown_rows": [0],
                                       "dropped_incompatible_allele_rows": [0]})()

        stdout = io.StringIO()
        with mock.patch.object(commands, "run_rg_from_args", return_value=expected), contextlib.redirect_stdout(stdout):
            result = cli.main(
                [
                    "rg",
                    "--ldscore-dir",
                    "ldscores",
                    "--sumstats-sources",
                    "a.sumstats.gz",
                    "b.sumstats.gz",
                    "--output-dir",
                    "out",
                ]
            )

        self.assertIs(result, expected)
        self.assertEqual(stdout.getvalue(), "")

    def test_cli_requires_output_dir_for_h2(self):
        expected = pd.DataFrame([{"trait_name": "trait", "total_h2_obs": 0.1, "total_h2_obs_se": np.nan}])

        stdout = io.StringIO()
        with mock.patch.object(commands, "run_h2_from_args", return_value=expected) as patched, contextlib.redirect_stdout(stdout):
            with self.assertRaises(SystemExit):
                cli.main(
                [
                    "h2",
                    "--ldscore-dir",
                    "ldscores",
                    "--sumstats-file",
                    "trait.sumstats.gz",
                ]
            )

        patched.assert_not_called()
        self.assertEqual(stdout.getvalue(), "")

    def test_cli_requires_output_dir_for_partitioned_h2(self):
        expected = pd.DataFrame([{"category": "query", "prop_snps": 1.0, "coefficient": np.nan}])

        stdout = io.StringIO()
        with mock.patch.object(
            commands,
            "run_partitioned_h2_from_args",
            return_value=expected,
        ) as patched, contextlib.redirect_stdout(stdout):
            with self.assertRaises(SystemExit):
                cli.main(
                [
                    "partitioned-h2",
                    "--ldscore-dir",
                    "ldscores",
                    "--sumstats-file",
                    "trait.sumstats.gz",
                ]
            )

        patched.assert_not_called()
        self.assertEqual(stdout.getvalue(), "")

    def test_run_partitioned_h2_from_args_refuses_stale_per_query_tree_without_overwrite(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            output_dir = tmpdir / "out"
            stale = output_dir / "diagnostics" / "query_annotations" / "old"
            stale.mkdir(parents=True)
            (stale / "metadata.json").write_text("{}\n", encoding="utf-8")
            args = parsed_regression_args(
                "partitioned-h2",
                {
                    "sumstats_file": str(tmpdir / "trait.sumstats.gz"),
                    "trait_name": "trait",
                    "ldscore_dir": str(tmpdir / "ldscores"),
                    "count_kind": "common",
                    "output_dir": str(output_dir),
                    "overwrite": False,
                    "num_blocks": 200,
                    "intercept_h2": None,
                    "two_step_cutoff": None,
                    "chisq_max": None,
                },
            )

            with mock.patch.object(
                partitioned,
                "estimate_partitioned_h2_batch",
                side_effect=AssertionError("estimation should not run after preflight failure"),
            ):
                with self.assertRaisesRegex(FileExistsError, "overwrite"):
                    commands.run_partitioned_h2_from_args(args)

            self.assertTrue(stale.exists())
            self.assertFalse((output_dir / "partitioned_h2.tsv").exists())
            self.assertFalse((output_dir / "partitioned-h2.log").exists())
            self.assertFalse((output_dir / "diagnostics" / "partitioned-h2.log").exists())

    def test_regression_cli_refuses_existing_result_file_by_default(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            set_global_config(GlobalConfig(snp_identifier="rsid"))
            with gzip.open(tmpdir / "trait.sumstats.gz", "wt", encoding="utf-8") as handle:
                handle.write("SNP\tA1\tA2\tZ\tN\nrs1\tA\tC\t1.0\t1000\n")
            self.write_sumstats_sidecar(tmpdir / "metadata.json", trait_name="trait")
            ldscore_dir = self.write_ldscore_dir(tmpdir / "ldscores", include_query=True)
            output_dir = tmpdir / "out"
            output_dir.mkdir()
            existing = output_dir / "h2.tsv"
            existing.write_text("existing\n", encoding="utf-8")
            args = parsed_regression_args(
                "h2",
                {
                    "sumstats_file": str(tmpdir / "trait.sumstats.gz"),
                    "trait_name": "trait",
                    "ldscore_dir": str(ldscore_dir),
                    "count_kind": "common",
                    "output_dir": str(output_dir),
                    "overwrite": False,
                    "num_blocks": 200,
                    "intercept_h2": None,
                    "two_step_cutoff": None,
                    "chisq_max": None,
                },
            )

            with stub_h2_kernel(
                mock.Mock(
                    tot=np.array([0.1]),
                    tot_se=np.array([0.01]),
                    intercept=np.array([1.0]),
                    intercept_se=0.01,
                    mean_chisq=np.array([1.1]),
                    lambda_gc=np.array([1.0]),
                    ratio=0.0,
                    ratio_se=0.0,
                    coef=np.array([0.001]),
                ),
            ):
                with self.assertRaisesRegex(FileExistsError, "overwrite"):
                    commands.run_h2_from_args(args)

            self.assertEqual(existing.read_text(encoding="utf-8"), "existing\n")
            self.assertFalse((output_dir / "h2.log").exists())
            self.assertFalse((output_dir / "diagnostics" / "h2.log").exists())

    def test_regression_cli_refuses_existing_h2_metadata_by_default_before_estimation(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            output_dir = tmpdir / "out"
            existing = output_dir / "diagnostics" / "metadata.json"
            existing.parent.mkdir(parents=True)
            existing.write_text("{}\n", encoding="utf-8")
            args = parsed_regression_args(
                "h2",
                {
                    "sumstats_file": str(tmpdir / "missing.sumstats.gz"),
                    "trait_name": "trait",
                    "ldscore_dir": str(tmpdir / "missing_ldscores"),
                    "count_kind": "common",
                    "output_dir": str(output_dir),
                    "overwrite": False,
                    "num_blocks": 200,
                    "intercept_h2": None,
                    "two_step_cutoff": None,
                    "chisq_max": None,
                },
            )

            with mock.patch.object(
                commands,
                "_fit_h2_dataset",
                side_effect=AssertionError("estimation should not run after preflight failure"),
            ):
                with self.assertRaisesRegex(FileExistsError, "overwrite"):
                    commands.run_h2_from_args(args)

            self.assertEqual(existing.read_text(encoding="utf-8"), "{}\n")
            self.assertFalse((output_dir / "h2.tsv").exists())
            self.assertFalse((output_dir / "h2.log").exists())
            self.assertFalse((output_dir / "diagnostics" / "h2.log").exists())

    def test_regression_cli_allows_existing_result_file_with_overwrite(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            set_global_config(GlobalConfig(snp_identifier="rsid"))
            with gzip.open(tmpdir / "trait.sumstats.gz", "wt", encoding="utf-8") as handle:
                handle.write("SNP\tA1\tA2\tZ\tN\nrs1\tA\tC\t1.0\t1000\n")
            self.write_sumstats_sidecar(tmpdir / "metadata.json", trait_name="trait")
            ldscore_dir = self.write_ldscore_dir(tmpdir / "ldscores", include_query=True)
            output_dir = tmpdir / "out"
            output_dir.mkdir()
            existing = output_dir / "h2.tsv"
            existing.write_text("existing\n", encoding="utf-8")
            existing_metadata = output_dir / "diagnostics" / "metadata.json"
            existing_metadata.parent.mkdir()
            existing_metadata.write_text('{"old": true}\n', encoding="utf-8")
            args = parsed_regression_args(
                "h2",
                {
                    "sumstats_file": str(tmpdir / "trait.sumstats.gz"),
                    "trait_name": "trait",
                    "ldscore_dir": str(ldscore_dir),
                    "count_kind": "common",
                    "output_dir": str(output_dir),
                    "overwrite": True,
                    "num_blocks": 200,
                    "intercept_h2": None,
                    "two_step_cutoff": None,
                    "chisq_max": None,
                },
            )

            with stub_h2_kernel(
                mock.Mock(
                    tot=np.array([0.1]),
                    tot_se=np.array([0.01]),
                    intercept=np.array([1.0]),
                    intercept_se=0.01,
                    mean_chisq=np.array([1.1]),
                    lambda_gc=np.array([1.0]),
                    ratio=0.0,
                    ratio_se=0.0,
                    coef=np.array([0.001]),
                ),
            ):
                commands.run_h2_from_args(args)

            self.assertIn("total_h2", existing.read_text(encoding="utf-8"))
            self.assertFalse((output_dir / "metadata.json").exists())
            metadata = json.loads(existing_metadata.read_text(encoding="utf-8"))
            self.assertNotIn("old", metadata)
            self.assertEqual(metadata["artifact_type"], "h2_result")

    def test_regression_writer_refuses_each_fixed_summary_filename(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir)
            from ldsc._regression.outputs import H2DirectoryWriter, PartitionedH2DirectoryWriter, RgDirectoryWriter

            for filename, writer in (("h2.tsv", H2DirectoryWriter), ("partitioned_h2.tsv", PartitionedH2DirectoryWriter), ("rg.tsv", RgDirectoryWriter)):
                with self.subTest(filename=filename):
                    existing = output_dir / filename
                    existing.write_text("existing\n", encoding="utf-8")

                    with self.assertRaisesRegex(FileExistsError, "overwrite"):
                        writer.artifact_family(output_dir).preflight(overwrite=False)

                    self.assertEqual(existing.read_text(encoding="utf-8"), "existing\n")


@pytest.mark.parametrize("command, flags, message", [
    ("h2", ["--samp-prev", "1.5", "--pop-prev", "0.01"], r"outside the open interval \(0, 1\)"),
    ("partitioned-h2", ["--num-blocks", "1"], "num_blocks=1 must be greater than 1"),
    ("rg", ["--num-blocks", "1"], "num_blocks=1 must be greater than 1"),
])
def test_regression_settings_are_checked_before_anything_is_written(tmp_path, command, flags, message):
    # The regression settings become a RegressionConfig before the output preflight, so an invalid one
    # stops the run before the output directory, its log or the input checks.
    inputs = ["--sumstats-sources", "a.parquet", "b.parquet"] if command == "rg" else ["--sumstats-file", "a.parquet"]
    with pytest.raises(LDSCConfigError, match=message):
        cli.main([command, *inputs, "--ldscore-dir", "missing", "--output-dir", str(tmp_path / "out"), *flags])
    assert not (tmp_path / "out").exists()
