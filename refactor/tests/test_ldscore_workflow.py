from argparse import Namespace
import gzip
from pathlib import Path
from types import SimpleNamespace
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
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            baseline = tmpdir / "baseline.annot"
            baseline.write_text("CHR\tBP\tSNP\tCM\tbase\n1\t10\trs1\t0.1\t1\n", encoding="utf-8")
            (tmpdir / "panel.bed").write_text("", encoding="utf-8")
            (tmpdir / "panel.bim").write_text("1 rs1 0.1 10 A G\n", encoding="utf-8")
            (tmpdir / "panel.fam").write_text("f i 0 0 0 -9\n", encoding="utf-8")
            args = Namespace(
                out=str(tmpdir / "example"),
                output_dir=None,
                query_annot=None,
                query_annot_chr=None,
                baseline_annot=str(baseline),
                baseline_annot_chr=None,
                bfile=str(tmpdir / "panel"),
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
            with mock.patch.object(ldscore_workflow.ldscore_new, "validate_args") as validate_args, mock.patch.object(
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

    def test_namespace_from_configs_emits_string_paths(self):
        from ldsc._kernel.ref_panel import RefPanelSpec

        common = CommonConfig(snp_identifier="rsid", genome_build="hg19")
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            r2_path = tmpdir / "r2" / "chr1.parquet"
            meta_path = tmpdir / "meta" / "chr1.tsv.gz"
            r2_path.parent.mkdir(parents=True, exist_ok=True)
            meta_path.parent.mkdir(parents=True, exist_ok=True)
            r2_path.write_text("", encoding="utf-8")
            with gzip.open(meta_path, "wt", encoding="utf-8") as handle:
                handle.write("CHR\tBP\tSNP\tCM\tMAF\n1\t10\trs1\t0.1\t0.2\n")
            spec = RefPanelSpec(
                backend="parquet_r2",
                r2_table_paths=(r2_path,),
                maf_metadata_paths=(meta_path,),
            )
            ref_panel = SimpleNamespace(spec=spec)

            args = ldscore_workflow._namespace_from_configs(
                chrom="1",
                ref_panel=ref_panel,
                ldscore_config=ldscore_workflow.LDScoreConfig(ld_wind_cm=1.0),
                common_config=common,
            )

            self.assertIsInstance(args.r2_table, str)
            self.assertIsInstance(args.frqfile, str)
            self.assertEqual(args.r2_table, str(r2_path))
            self.assertEqual(args.frqfile, str(meta_path))
            self.assertEqual(args.genome_build, "hg19")

    def test_run_ldscore_from_args_resolves_glob_and_suite_tokens_in_workflow_layer(self):
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
            ld_scores=np.array([[1.0]], dtype=np.float32),
            w_ld=np.array([[3.0]], dtype=np.float32),
            M=np.array([5.0]),
            M_5_50=np.array([4.0]),
            ldscore_columns=["base"],
            baseline_columns=["base"],
            query_columns=[],
        )
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            with gzip.open(tmpdir / "baseline.1.annot.gz", "wt", encoding="utf-8") as handle:
                handle.write("CHR\tBP\tSNP\tCM\tbase\n1\t10\trs1\t0.1\t1\n")
            (tmpdir / "panel.1.bed").write_text("", encoding="utf-8")
            (tmpdir / "panel.1.bim").write_text("1 rs1 0.1 10 A G\n", encoding="utf-8")
            (tmpdir / "panel.1.fam").write_text("f i 0 0 0 -9\n", encoding="utf-8")
            args = Namespace(
                out=str(tmpdir / "example"),
                output_dir=None,
                query_annot=None,
                query_annot_chr=None,
                baseline_annot=str(tmpdir / "baseline.@"),
                baseline_annot_chr=None,
                bfile=str(tmpdir / "panel.@"),
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
            with mock.patch.object(ldscore_workflow.ldscore_new, "validate_args"), mock.patch.object(
                ldscore_workflow.ldscore_new,
                "combine_annotation_groups",
                return_value=mock.sentinel.bundle,
            ) as combine_groups, mock.patch.object(
                ldscore_workflow.ldscore_new,
                "compute_chrom_from_plink",
                return_value=fake_legacy_result,
            ):
                ldscore_workflow.run_ldscore_from_args(args)

            combine_kwargs = combine_groups.call_args.kwargs
            self.assertEqual(combine_kwargs["baseline_files"], [str(tmpdir / "baseline.1.annot.gz")])
