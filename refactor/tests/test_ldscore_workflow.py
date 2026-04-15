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

from ldsc.config import CommonConfig
from ldsc.outputs import OutputSpec

try:
    from ldsc import ldscore_calculator as ldscore_workflow
except ImportError:
    ldscore_workflow = None


@unittest.skipIf(ldscore_workflow is None, "ldscore_workflow module is not available")
class LDScoreWorkflowTest(unittest.TestCase):
    def make_chrom_result(self, chrom: str, bp: int, score: float, count: float):
        metadata = pd.DataFrame(
            {
                "CHR": [chrom],
                "SNP": [f"rs{chrom}"],
                "BP": [bp],
                "CM": [0.1],
                "MAF": [0.2],
            }
        )
        return ldscore_workflow.ChromLDScoreResult(
            chrom=chrom,
            reference_metadata=metadata,
            ld_scores=pd.DataFrame({"base": [score], "query": [score + 1.0]}),
            regression_metadata=metadata.copy(),
            w_ld=pd.DataFrame({"L2": [score + 2.0]}),
            snp_count_totals={
                "all_reference_snp_counts": np.array([count, count + 1.0]),
                "common_reference_snp_counts_maf_gt_0_05": np.array([count - 1.0, count]),
            },
            baseline_columns=["base"],
            query_columns=["query"],
            reference_snps={f"rs{chrom}"},
            regression_snps={f"rs{chrom}"},
        )

    def test_aggregate_chromosome_results(self):
        calc = ldscore_workflow.LDScoreCalculator()
        result = calc._aggregate_chromosome_results(
            [
                self.make_chrom_result("2", 20, 2.0, 20.0),
                self.make_chrom_result("1", 10, 1.0, 10.0),
            ],
            common_config=CommonConfig(snp_identifier="rsid"),
        )
        self.assertEqual(result.reference_metadata["CHR"].tolist(), ["1", "2"])
        self.assertEqual(result.ld_scores.columns.tolist(), ["base", "query"])
        self.assertEqual(result.w_ld.columns.tolist(), ["L2"])
        np.testing.assert_allclose(result.snp_count_totals["all_reference_snp_counts"], [30.0, 32.0])
        np.testing.assert_allclose(
            result.snp_count_totals["common_reference_snp_counts_maf_gt_0_05"],
            [28.0, 30.0],
        )

    def test_run_ldscore_from_args_writes_outputs(self):
        fake_legacy_result = ldscore_workflow._LegacyChromResult(
            chrom="1",
            metadata=pd.DataFrame(
                {
                    "CHR": ["1"],
                    "SNP": ["rs1"],
                    "BP": [10],
                    "CM": [0.1],
                    "MAF": [0.2],
                }
            ),
            ld_scores=np.array([[1.0, 2.0]], dtype=np.float32),
            w_ld=np.array([[3.0]], dtype=np.float32),
            M=np.array([5.0, 6.0]),
            M_5_50=np.array([4.0, 5.0]),
            ldscore_columns=["base", "query"],
            baseline_columns=["base"],
            query_columns=["query"],
        )
        args = Namespace(
            out="example",
            output_dir=None,
            query_annot=None,
            query_annot_chr=None,
            baseline_annot="baseline.annot",
            baseline_annot_chr=None,
            bfile="panel",
            bfile_chr=None,
            r2_table=None,
            r2_table_chr=None,
            snp_identifier="rsid",
            genome_build=None,
            r2_bias_mode=None,
            r2_sample_size=None,
            regression_snps=None,
            frqfile=None,
            frqfile_chr=None,
            ld_wind_snps=10,
            ld_wind_kb=None,
            ld_wind_cm=None,
            maf=None,
            chunk_size=50,
            per_chr_output=False,
            yes_really=False,
            log_level="INFO",
        )
        with tempfile.TemporaryDirectory() as tmpdir:
            args.out = str(Path(tmpdir) / "example")
            with mock.patch.object(ldscore_workflow.ldscore_new, "validate_args") as validate_args, mock.patch.object(
                ldscore_workflow.ldscore_new, "chromosome_set_from_annotation_inputs", return_value=["1"]
            ), mock.patch.object(
                ldscore_workflow.ldscore_new,
                "resolve_optional_chr_files",
                return_value=[],
            ), mock.patch.object(
                ldscore_workflow.ldscore_new,
                "resolve_annotation_files",
                return_value=["baseline.annot"],
            ), mock.patch.object(
                ldscore_workflow.ldscore_new,
                "combine_annotation_groups",
                return_value=mock.sentinel.bundle,
            ), mock.patch.object(
                ldscore_workflow.ldscore_new,
                "compute_chrom_from_plink",
                return_value=fake_legacy_result,
            ):
                result = ldscore_workflow.run_ldscore_from_args(args)
            called_args = validate_args.call_args[0][0]
            self.assertEqual(called_args.out, args.out)
            self.assertEqual(called_args.snp_identifier, "rsID")
            self.assertEqual(result.reference_metadata["SNP"].tolist(), ["rs1"])
            self.assertIn("ldscore", result.output_paths)
            self.assertTrue(Path(result.output_paths["ldscore"]).exists())
            self.assertTrue(Path(result.output_paths["w_ld"]).exists())
