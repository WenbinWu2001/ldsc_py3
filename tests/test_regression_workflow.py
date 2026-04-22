from dataclasses import replace
import gzip
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

from ldsc.config import ConfigMismatchError, GlobalConfig, RegressionConfig

try:
    from ldsc.ldscore_calculator import LDScoreResult
    from ldsc import regression_runner
    from ldsc.regression_runner import RegressionRunner
    from ldsc.sumstats_munger import SumstatsTable
    from ldsc.annotation_builder import AnnotationBundle
except ImportError:
    LDScoreResult = None
    RegressionRunner = None
    SumstatsTable = None
    AnnotationBundle = None


@unittest.skipIf(RegressionRunner is None, "regression_runner module is not available")
class RegressionWorkflowTest(unittest.TestCase):
    def test_load_ldscore_from_files_is_public(self):
        from ldsc import load_ldscore_from_files

        self.assertTrue(callable(load_ldscore_from_files))

    def test_load_ldscore_from_files_new_format_no_weight_path(self):
        from ldsc import load_ldscore_from_files

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            ldscore_path = tmpdir / "trait.l2.ldscore.gz"
            counts_path = tmpdir / "trait.l2.M_5_50"
            with gzip.open(ldscore_path, "wt", encoding="utf-8") as handle:
                handle.write("CHR\tSNP\tBP\tbase\tquery\tregr_weight\n1\trs1\t10\t1.0\t2.0\t3.0\n")
            counts_path.write_text("5\t6\n", encoding="utf-8")

            result = load_ldscore_from_files(
                ldscore_path=str(ldscore_path),
                counts_path=str(counts_path),
                count_kind="m_5_50",
                snp_identifier="rsid",
            )

        self.assertIn("regr_weight", result.ldscore_table.columns)
        self.assertEqual(result.ld_reference_snps, frozenset())
        self.assertEqual(result.ld_regression_snps, frozenset({"rs1"}))

    def make_ldscore_result(self):
        ldscore_table = pd.DataFrame(
            {
                "CHR": ["1", "1", "1"],
                "SNP": ["rs1", "rs2", "rs3"],
                "BP": [10, 20, 30],
                "base": [1.0, 1.0, 1.0],
                "query1": [1.0, 2.0, 3.0],
                "query2": [0.5, 1.5, 2.5],
                "regr_weight": [2.0, 2.0, 2.0],
            }
        )
        return LDScoreResult(
            ldscore_table=ldscore_table,
            snp_count_totals={
                "all_reference_snp_counts": np.array([10.0, 20.0, 30.0]),
                "common_reference_snp_counts_maf_gt_0_05": np.array([8.0, 18.0, 28.0]),
            },
            baseline_columns=["base"],
            query_columns=["query1", "query2"],
            ld_reference_snps=frozenset(),
            ld_regression_snps=frozenset({"rs1", "rs2", "rs3"}),
            chromosome_results=[],
            config_snapshot=GlobalConfig(genome_build="hg38", snp_identifier="rsid"),
        )

    def make_sumstats_table(self):
        return SumstatsTable(
            data=pd.DataFrame(
                {
                    "SNP": ["rs1", "rs2", "rs3"],
                    "Z": [2.0, 1.0, 0.5],
                    "N": [1000.0, 1000.0, 1000.0],
                    "A1": ["A", "C", "G"],
                    "A2": ["G", "T", "A"],
                }
            ),
            has_alleles=True,
            source_path="sumstats.gz",
            trait_name="trait",
            provenance={},
            config_snapshot=GlobalConfig(genome_build="hg38", snp_identifier="rsid"),
        )

    def make_annotation_bundle(self):
        metadata = pd.DataFrame(
            {
                "CHR": ["1", "1", "1"],
                "BP": [10, 20, 30],
                "SNP": ["rs1", "rs2", "rs3"],
                "CM": [0.1, 0.2, 0.3],
            }
        )
        return AnnotationBundle(
            metadata=metadata,
            baseline_annotations=pd.DataFrame({"base": [1, 1, 1]}),
            query_annotations=pd.DataFrame({"query1": [1, 0, 0], "query2": [0, 1, 0]}),
            baseline_columns=["base"],
            query_columns=["query1", "query2"],
            chromosomes=["1"],
            source_summary={},
        )

    def test_build_dataset_uses_common_counts_and_drops_zero_variance_columns(self):
        runner = RegressionRunner(GlobalConfig(snp_identifier="rsid"), RegressionConfig())
        dataset = runner.build_dataset(self.make_sumstats_table(), self.make_ldscore_result())
        self.assertEqual(dataset.count_key_used_for_regression, "common_reference_snp_counts_maf_gt_0_05")
        self.assertEqual(dataset.retained_ld_columns, ["query1", "query2"])
        self.assertEqual(dataset.dropped_zero_variance_ld_columns, ["base"])
        self.assertEqual(dataset.config_snapshot, GlobalConfig(genome_build="hg38", snp_identifier="rsid"))

    def test_build_dataset_raises_on_mismatched_sumstats_and_ldscore_snapshots(self):
        runner = RegressionRunner(GlobalConfig(snp_identifier="rsid"), RegressionConfig())
        sumstats = self.make_sumstats_table()
        ldscore_result = replace(
            self.make_ldscore_result(),
            config_snapshot=GlobalConfig(genome_build="hg19", snp_identifier="rsid"),
        )

        with self.assertRaisesRegex(ConfigMismatchError, "genome_build mismatch"):
            runner.build_dataset(sumstats, ldscore_result)

    def test_build_dataset_skips_compatibility_check_for_legacy_none_snapshots(self):
        runner = RegressionRunner(GlobalConfig(snp_identifier="rsid"), RegressionConfig())
        sumstats = replace(self.make_sumstats_table(), config_snapshot=None)
        ldscore_result = replace(self.make_ldscore_result(), config_snapshot=None)

        dataset = runner.build_dataset(sumstats, ldscore_result)

        self.assertIsNone(dataset.config_snapshot)

    def test_estimate_partitioned_h2_batch_loops_over_queries(self):
        runner = RegressionRunner(GlobalConfig(snp_identifier="rsid"), RegressionConfig())
        table = self.make_sumstats_table()
        ldscore_result = self.make_ldscore_result()
        annotation_bundle = self.make_annotation_bundle()
        fake_rows = [
            pd.DataFrame([{"query_annotation": "query1", "coefficient": 1.0}]),
            pd.DataFrame([{"query_annotation": "query2", "coefficient": 2.0}]),
        ]
        with mock.patch.object(runner, "estimate_partitioned_h2", side_effect=fake_rows) as patched:
            result = runner.estimate_partitioned_h2_batch(table, ldscore_result, annotation_bundle)
        self.assertEqual(patched.call_count, 2)
        self.assertEqual(result["query_annotation"].tolist(), ["query1", "query2"])

    def test_run_h2_from_args_resolves_scalar_glob_inputs(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            with gzip.open(tmpdir / "trait.sumstats.gz", "wt", encoding="utf-8") as handle:
                handle.write("SNP\tZ\tN\nrs1\t1.0\t1000\n")
            with gzip.open(tmpdir / "trait.l2.ldscore.gz", "wt", encoding="utf-8") as handle:
                handle.write("CHR\tSNP\tPOS\tCM\tbase\n1\trs1\t10\t0.1\t1.0\n")
            with gzip.open(tmpdir / "trait.w.l2.ldscore.gz", "wt", encoding="utf-8") as handle:
                handle.write("CHR\tSNP\tPOS\tCM\tL2\n1\trs1\t10\t0.1\t2.0\n")
            (tmpdir / "trait.l2.M_5_50").write_text("5\n", encoding="utf-8")
            args = type(
                "Args",
                (),
                {
                    "sumstats": str(tmpdir / "trait*.sumstats.gz"),
                    "trait_name": "trait",
                    "ldscore": str(tmpdir / "trait.l2*.ldscore.gz"),
                    "w_ld": str(tmpdir / "trait.w*.ldscore.gz"),
                    "counts": str(tmpdir / "*.l2.M_5_50"),
                    "count_kind": "m_5_50",
                    "out": None,
                    "n_blocks": 200,
                    "no_intercept": False,
                    "intercept_h2": None,
                    "two_step_cutoff": None,
                    "chisq_max": None,
                },
            )()

            with mock.patch.object(regression_runner.RegressionRunner, "estimate_h2", return_value=mock.Mock(
                tot=np.array([0.1]),
                tot_se=np.array([0.01]),
                intercept=np.array([1.0]),
                intercept_se=0.01,
                mean_chisq=np.array([1.1]),
                lambda_gc=np.array([1.0]),
                ratio=0.0,
                ratio_se=0.0,
            )) as patched:
                summary = regression_runner.run_h2_from_args(args)

            patched.assert_called_once()
            self.assertEqual(summary.loc[0, "trait_name"], "trait")

    def test_load_ldscore_result_from_files_requires_canonical_internal_headers(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            ldscore_path = tmpdir / "trait.l2.ldscore.gz"
            weight_path = tmpdir / "trait.w.l2.ldscore.gz"
            counts_path = tmpdir / "trait.l2.M_5_50"
            with gzip.open(ldscore_path, "wt", encoding="utf-8") as handle:
                handle.write("CHR\tSNP\tBP\tCM\tbase\n1\trs1\t10\t0.1\t1.0\n")
            with gzip.open(weight_path, "wt", encoding="utf-8") as handle:
                handle.write("CHR\tSNP\tBP\tCM\tL2\n1\trs1\t10\t0.1\t2.0\n")
            counts_path.write_text("5\n", encoding="utf-8")

            with self.assertRaises(ValueError):
                regression_runner._load_ldscore_result_from_files(
                    ldscore_path=str(ldscore_path),
                    weight_path=str(weight_path),
                    counts_path=str(counts_path),
                    count_kind="m_5_50",
                )

    def test_run_h2_from_args_without_w_ld_uses_embedded_regr_weight(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            with gzip.open(tmpdir / "trait.sumstats.gz", "wt", encoding="utf-8") as handle:
                handle.write("SNP\tZ\tN\nrs1\t1.0\t1000\n")
            with gzip.open(tmpdir / "trait.l2.ldscore.gz", "wt", encoding="utf-8") as handle:
                handle.write("CHR\tSNP\tBP\tbase\tregr_weight\n1\trs1\t10\t1.0\t2.0\n")
            (tmpdir / "trait.l2.M_5_50").write_text("5\n", encoding="utf-8")
            args = type(
                "Args",
                (),
                {
                    "sumstats": str(tmpdir / "trait.sumstats.gz"),
                    "trait_name": "trait",
                    "ldscore": str(tmpdir / "trait.l2.ldscore.gz"),
                    "w_ld": None,
                    "counts": str(tmpdir / "trait.l2.M_5_50"),
                    "count_kind": "m_5_50",
                    "out": None,
                    "n_blocks": 200,
                    "no_intercept": False,
                    "intercept_h2": None,
                    "two_step_cutoff": None,
                    "chisq_max": None,
                },
            )()

            with mock.patch.object(regression_runner.RegressionRunner, "estimate_h2", return_value=mock.Mock(
                tot=np.array([0.1]),
                tot_se=np.array([0.01]),
                intercept=np.array([1.0]),
                intercept_se=0.01,
                mean_chisq=np.array([1.1]),
                lambda_gc=np.array([1.0]),
                ratio=0.0,
                ratio_se=0.0,
            )) as patched:
                summary = regression_runner.run_h2_from_args(args)

        patched.assert_called_once()
        self.assertEqual(summary.loc[0, "trait_name"], "trait")
