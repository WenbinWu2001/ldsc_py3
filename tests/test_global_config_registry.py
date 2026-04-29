from argparse import Namespace
from pathlib import Path
import sys
import tempfile
import unittest
from unittest import mock

import numpy as np
import pandas as pd

SRC = Path(__file__).resolve().parents[1] / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

import ldsc
from ldsc import GlobalConfig, set_global_config, reset_global_config
from ldsc.config import ConfigMismatchError, validate_config_compatibility
from ldsc import ldscore_calculator, ref_panel_builder, regression_runner
from ldsc._kernel import annotation as kernel_annotation
from ldsc.ldscore_calculator import LDScoreResult
from ldsc.ref_panel_builder import ReferencePanelBuildConfig, ReferencePanelBuilder
from ldsc.regression_runner import RegressionRunner
from ldsc.sumstats_munger import SumstatsTable


def _make_build_config(tmpdir: Path) -> ReferencePanelBuildConfig:
    return ReferencePanelBuildConfig(
        plink_prefix=tmpdir / "panel.@",
        source_genome_build="hg19",
        genetic_map_hg19_sources=tmpdir / "hg19.map",
        genetic_map_hg38_sources=tmpdir / "hg38.map",
        liftover_chain_hg19_to_hg38_file=tmpdir / "hg19ToHg38.over.chain",
        output_dir=tmpdir / "out",
        ld_wind_kb=1.0,
    )


def _make_sumstats_table() -> SumstatsTable:
    return SumstatsTable(
        data=pd.DataFrame(
            {
                "SNP": ["rs1", "rs2"],
                "Z": [2.0, 1.0],
                "N": [1000.0, 1000.0],
                "A1": ["A", "C"],
                "A2": ["G", "T"],
            }
        ),
        has_alleles=True,
        source_path="sumstats.gz",
        trait_name="trait",
        provenance={},
    )


def _make_ldscore_result() -> LDScoreResult:
    baseline_table = pd.DataFrame(
        {
            "CHR": ["1", "1"],
            "SNP": ["rs1", "rs2"],
            "BP": [10, 20],
            "regr_weight": [3.0, 4.0],
            "base": [1.0, 2.0],
        }
    )
    return LDScoreResult(
        baseline_table=baseline_table,
        query_table=None,
        count_records=[
            {
                "group": "baseline",
                "column": "base",
                "all_reference_snp_count": 10.0,
            }
        ],
        baseline_columns=["base"],
        query_columns=[],
        ld_reference_snps=frozenset(),
        ld_regression_snps=frozenset({"rs1", "rs2"}),
        chromosome_results=[],
    )


class GlobalConfigRegistryTest(unittest.TestCase):
    def tearDown(self) -> None:
        reset_global_config()

    def test_run_bed_to_annot_uses_registered_global_config_and_logs_once(self):
        set_global_config(
            GlobalConfig(
                snp_identifier="rsid",
                log_level="DEBUG",
            )
        )

        with mock.patch.object(
            kernel_annotation.AnnotationBuilder,
            "project_bed_annotations",
            autospec=True,
            return_value=mock.sentinel.bundle,
        ) as patched, mock.patch("builtins.print") as patched_print:
            with self.assertLogs("LDSC.config", level="INFO") as caught:
                result = ldsc.run_bed_to_annot(
                    query_annot_bed_sources="beds/*.bed",
                    baseline_annot_sources="baseline.@.annot.gz",
                    output_dir="out",
                )

        self.assertIs(result, mock.sentinel.bundle)
        builder = patched.call_args.args[0]
        self.assertEqual(builder.global_config.snp_identifier, "rsid")
        self.assertIsNone(builder.global_config.genome_build)
        self.assertEqual(builder.global_config.log_level, "DEBUG")
        patched_print.assert_not_called()
        self.assertEqual(len(caught.records), 1)
        self.assertEqual(caught.records[0].levelname, "INFO")
        self.assertIn("GlobalConfig", caught.records[0].getMessage())

    def test_run_bed_to_annot_rejects_removed_shared_kwargs(self):
        with self.assertRaises(TypeError):
            ldsc.run_bed_to_annot(
                query_annot_bed_sources="beds/*.bed",
                baseline_annot_sources="baseline.@.annot.gz",
                output_dir="out",
                snp_identifier="rsid",
            )

    def test_run_ldscore_uses_registered_global_config_for_python_defaults(self):
        set_global_config(GlobalConfig(snp_identifier="rsid", log_level="DEBUG"))

        def _assert_args(args):
            self.assertEqual(args.snp_identifier, "rsid")
            self.assertIsNone(args.genome_build)
            self.assertEqual(args.log_level, "DEBUG")
            raise RuntimeError("stop after inspecting args")

        with mock.patch.object(ldscore_calculator, "_normalize_run_args", side_effect=_assert_args):
            with self.assertRaisesRegex(RuntimeError, "stop after inspecting args"):
                ldscore_calculator.run_ldscore(
                    output_dir="results/example",
                    baseline_annot_sources="baseline.annot.gz",
                    r2_sources="reference.parquet",
                    metadata_sources="reference.tsv.gz",
                    ld_wind_cm=1.0,
                )

    def test_run_ldscore_rejects_removed_shared_kwargs(self):
        with self.assertRaisesRegex(ValueError, "set_global_config"):
            ldscore_calculator.run_ldscore(
                output_dir="results/example",
                baseline_annot_sources="baseline.annot.gz",
                r2_sources="reference.parquet",
                metadata_sources="reference.tsv.gz",
                ld_wind_cm=1.0,
                snp_identifier="rsid",
            )

    def test_run_ldscore_rejects_retired_io_kwargs(self):
        with self.assertRaisesRegex(ValueError, "baseline_annot_paths"):
            ldscore_calculator.run_ldscore(
                output_dir="results/example",
                baseline_annot_paths="baseline.annot.gz",
                plink_path="plink/panel",
                ld_wind_snps=10,
            )

    def test_run_build_ref_panel_rejects_removed_shared_kwargs(self):
        with self.assertRaisesRegex(ValueError, "set_global_config"):
            ref_panel_builder.run_build_ref_panel(
                plink_prefix="plink/panel",
                source_genome_build="hg19",
                genetic_map_hg19_sources="maps/hg19.txt",
                genetic_map_hg38_sources="maps/hg38.txt",
                liftover_chain_hg19_to_hg38_file="chains/hg19ToHg38.over.chain",
                output_dir="out/panel",
                ld_wind_kb=1.0,
                log_level="DEBUG",
            )

    def test_run_build_ref_panel_rejects_retired_io_kwargs(self):
        with self.assertRaisesRegex(ValueError, "plink_path"):
            ref_panel_builder.run_build_ref_panel(
                plink_path="plink/panel",
                source_genome_build="hg19",
                genetic_map_hg19_path="maps/hg19.txt",
                output_dir="out/panel",
                ld_wind_kb=1.0,
            )

    def test_reference_panel_builder_uses_registered_global_config_and_logs_once(self):
        set_global_config(GlobalConfig(snp_identifier="rsid", log_level="DEBUG"))

        with tempfile.TemporaryDirectory() as tmpdir:
            config = _make_build_config(Path(tmpdir))
            builder = ReferencePanelBuilder()

            with mock.patch("ldsc.ref_panel_builder.resolve_plink_prefix_group", return_value=["panel.1"]), mock.patch.object(
                builder, "_configure_logging"
            ), mock.patch.object(
                builder, "_prepare_build_state", return_value=mock.sentinel.state
            ), mock.patch.object(
                builder, "_discover_prefix_chromosomes", return_value=["1"]
            ), mock.patch.object(
                builder,
                "_build_chromosome",
                return_value={"ann": "ann", "ld": "ld", "meta_hg19": "m19", "meta_hg38": "m38"},
            ), mock.patch("builtins.print") as patched_print:
                with self.assertLogs("LDSC.config", level="INFO") as caught:
                    result = builder.run(config)

        self.assertEqual(builder.global_config, ldsc.get_global_config())
        self.assertEqual(result.chromosomes, ["1"])
        patched_print.assert_not_called()
        self.assertEqual(len(caught.records), 1)
        self.assertEqual(caught.records[0].levelname, "INFO")
        self.assertIn("GlobalConfig", caught.records[0].getMessage())

    def test_regression_runner_uses_registered_global_config_and_logs_once(self):
        set_global_config(GlobalConfig(snp_identifier="rsid"))
        runner = RegressionRunner()

        with mock.patch("builtins.print") as patched_print:
            with self.assertLogs("LDSC.config", level="INFO") as caught:
                dataset = runner.build_dataset(_make_sumstats_table(), _make_ldscore_result())

        self.assertEqual(runner.global_config, ldsc.get_global_config())
        self.assertEqual(dataset.merged["SNP"].tolist(), ["rs1", "rs2"])
        patched_print.assert_not_called()
        self.assertEqual(len(caught.records), 1)
        self.assertEqual(caught.records[0].levelname, "INFO")
        self.assertIn("GlobalConfig", caught.records[0].getMessage())

    def test_cli_normalization_ignores_registered_python_global_config(self):
        set_global_config(GlobalConfig(snp_identifier="chr_pos", genome_build="hg38", log_level="WARNING"))

        args = Namespace(
            snp_identifier="rsid",
            genome_build=None,
            log_level="DEBUG",
            output_dir="results/example",
            keep_indivs_file=None,
        )

        normalized_args, global_config = ldscore_calculator._normalize_run_args(args)

        self.assertEqual(normalized_args.snp_identifier, "rsid")
        self.assertEqual(global_config.snp_identifier, "rsid")
        self.assertIsNone(global_config.genome_build)
        self.assertEqual(global_config.log_level, "DEBUG")


class ConfigCompatibilityTest(unittest.TestCase):
    def test_validate_config_compatibility_rejects_genome_build_mismatch(self):
        left = GlobalConfig(genome_build="hg19", snp_identifier="chr_pos")
        right = GlobalConfig(genome_build="hg38", snp_identifier="chr_pos")

        with self.assertRaisesRegex(ConfigMismatchError, "genome_build mismatch"):
            validate_config_compatibility(left, right, context="test")

    def test_validate_config_compatibility_rejects_snp_identifier_mismatch(self):
        with self.assertWarns(UserWarning):
            cfg = GlobalConfig(genome_build="hg38", snp_identifier="rsid")
        self.assertIsNone(cfg.genome_build)
        left = cfg
        right = GlobalConfig(genome_build="hg38", snp_identifier="chr_pos")

        with self.assertRaisesRegex(ConfigMismatchError, "snp_identifier mismatch"):
            validate_config_compatibility(left, right, context="test")

    def test_validate_config_compatibility_accepts_matching_critical_fields(self):
        shared = GlobalConfig(genome_build="hg38", snp_identifier="chr_pos")

        validate_config_compatibility(shared, shared, context="test")
