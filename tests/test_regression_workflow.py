from dataclasses import replace
import unittest
from unittest import mock

import numpy as np
import pandas as pd

from ldscore.config import CommonConfig, RegressionConfig

try:
    from ldscore.ldscore_workflow import LDScoreResult
    from ldscore.regression_workflow import RegressionRunner
    from ldscore.sumstats_munger import SumstatsTable
    from ldscore.annotation import AnnotationBundle
except ImportError:
    LDScoreResult = None
    RegressionRunner = None
    SumstatsTable = None
    AnnotationBundle = None


@unittest.skipIf(RegressionRunner is None, "regression_workflow module is not available")
class RegressionWorkflowTest(unittest.TestCase):
    def make_ldscore_result(self):
        metadata = pd.DataFrame(
            {
                "CHR": ["1", "1", "1"],
                "SNP": ["rs1", "rs2", "rs3"],
                "BP": [10, 20, 30],
                "CM": [0.1, 0.2, 0.3],
            }
        )
        return LDScoreResult(
            reference_metadata=metadata,
            ld_scores=pd.DataFrame({"base": [1.0, 1.0, 1.0], "query1": [1.0, 2.0, 3.0], "query2": [0.5, 1.5, 2.5]}),
            regression_metadata=metadata.copy(),
            w_ld=pd.DataFrame({"L2": [2.0, 2.0, 2.0]}),
            snp_count_totals={
                "all_reference_snp_counts": np.array([10.0, 20.0, 30.0]),
                "common_reference_snp_counts_maf_gt_0_05": np.array([8.0, 18.0, 28.0]),
            },
            baseline_columns=["base"],
            query_columns=["query1", "query2"],
            reference_snps={"rs1", "rs2", "rs3"},
            regression_snps={"rs1", "rs2", "rs3"},
            chromosome_results=[],
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
        runner = RegressionRunner(CommonConfig(snp_identifier="rsid"), RegressionConfig())
        dataset = runner.build_dataset(self.make_sumstats_table(), self.make_ldscore_result())
        self.assertEqual(dataset.count_key_used_for_regression, "common_reference_snp_counts_maf_gt_0_05")
        self.assertEqual(dataset.retained_ld_columns, ["query1", "query2"])
        self.assertEqual(dataset.dropped_zero_variance_ld_columns, ["base"])

    def test_estimate_partitioned_h2_batch_loops_over_queries(self):
        runner = RegressionRunner(CommonConfig(snp_identifier="rsid"), RegressionConfig())
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
