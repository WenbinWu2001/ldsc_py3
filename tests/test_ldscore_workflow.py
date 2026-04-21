from __future__ import annotations

from argparse import Namespace
import gzip
import importlib.util
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

from ldsc.config import GlobalConfig
from ldsc.outputs import OutputSpec

try:
    from ldsc import AnnotationBundle, LDScoreConfig, PlinkRefPanel, RefPanelSpec
    from ldsc import ldscore_calculator as ldscore_workflow
    from ldsc._kernel import formats as kernel_formats
    from ldsc._kernel import ldscore as kernel_ldscore
except ImportError:
    AnnotationBundle = None
    LDScoreConfig = None
    PlinkRefPanel = None
    RefPanelSpec = None
    kernel_formats = None
    ldscore_workflow = None
    kernel_ldscore = None


FIXTURES = Path(__file__).resolve().parent / "fixtures" / "legacy"
PLINK_FIXTURES = FIXTURES / "plink_test"
_HAS_PYARROW = importlib.util.find_spec("pyarrow") is not None
_HAS_BITARRAY = importlib.util.find_spec("bitarray") is not None


@unittest.skipIf(ldscore_workflow is None, "ldscore_workflow module is not available")
class LDScoreWorkflowTest(unittest.TestCase):
    def make_chrom_result(self, chrom: str, bp: int, score: float, count: float):
        metadata = pd.DataFrame(
            {
                "CHR": [chrom],
                "SNP": [f"rs{chrom}"],
                "POS": [bp],
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

    def _build_annotation_bundle(self, prefix: Path) -> AnnotationBundle:
        bim = pd.read_csv(
            prefix.with_suffix(".bim"),
            sep=r"\s+",
            header=None,
            names=["CHR", "SNP", "CM", "POS", "A1", "A2"],
        )
        metadata = bim.loc[:, ["CHR", "SNP", "CM", "POS"]].copy()
        metadata["CHR"] = metadata["CHR"].astype(str)
        metadata["SNP"] = metadata["SNP"].astype(str)
        metadata["CM"] = pd.to_numeric(metadata["CM"], errors="coerce")
        metadata["POS"] = pd.to_numeric(metadata["POS"], errors="raise").astype(int)
        baseline = pd.DataFrame({"base": np.ones(len(metadata), dtype=np.float32)})
        query = pd.DataFrame(index=metadata.index)
        return AnnotationBundle(
            metadata=metadata,
            baseline_annotations=baseline,
            query_annotations=query,
            baseline_columns=["base"],
            query_columns=[],
            chromosomes=["1"],
            source_summary={},
        )

    def _copy_plink_fixture_with_distinct_fids(self, tmpdir: Path) -> Path:
        prefix = tmpdir / "panel"
        prefix.with_suffix(".bed").write_bytes((PLINK_FIXTURES / "plink.bed").read_bytes())
        prefix.with_suffix(".bim").write_text((PLINK_FIXTURES / "plink.bim").read_text(encoding="utf-8"), encoding="utf-8")
        fam_lines = []
        for idx, line in enumerate((PLINK_FIXTURES / "plink.fam").read_text(encoding="utf-8").splitlines()):
            fields = line.split()
            fields[0] = f"fam{idx}"
            fam_lines.append(" ".join(fields))
        prefix.with_suffix(".fam").write_text("\n".join(fam_lines) + "\n", encoding="utf-8")
        return prefix

    def _write_keep_file(self, path: Path, iids: list[str]) -> None:
        path.write_text("\n".join(iids) + "\n", encoding="utf-8")

    def _expected_plink_result(self, prefix: Path, keep_path: Path, maf_min: float | None) -> tuple[pd.DataFrame, np.ndarray]:
        fam = kernel_formats.PlinkFAMFile(str(prefix.with_suffix(".fam")))
        keep_indivs = fam.loj(kernel_formats.FilterFile(str(keep_path)).IDList)
        bim = kernel_formats.PlinkBIMFile(str(prefix.with_suffix(".bim")))
        bed = kernel_ldscore.PlinkBEDFile(
            str(prefix.with_suffix(".bed")),
            len(fam.IDList),
            bim,
            keep_indivs=keep_indivs,
            mafMin=maf_min,
        )
        metadata = pd.DataFrame(bed.df, columns=bed.colnames).rename(columns={"BP": "POS"})
        metadata["CHR"] = metadata["CHR"].astype(str)
        metadata["SNP"] = metadata["SNP"].astype(str)
        metadata["CM"] = pd.to_numeric(metadata["CM"], errors="coerce")
        metadata["POS"] = pd.to_numeric(metadata["POS"], errors="raise").astype(int)
        metadata["MAF"] = pd.to_numeric(metadata["MAF"], errors="coerce")
        block_left = kernel_ldscore.getBlockLefts(np.arange(bed.m), 10)
        ld_scores = bed.ldScoreVarBlocks(block_left, 50, annot=np.ones((bed.m, 1), dtype=np.float32))
        return metadata.reset_index(drop=True), np.ravel(ld_scores)

    def test_aggregate_chromosome_results(self):
        calc = ldscore_workflow.LDScoreCalculator()
        result = calc._aggregate_chromosome_results(
            [
                self.make_chrom_result("2", 20, 2.0, 20.0),
                self.make_chrom_result("1", 10, 1.0, 10.0),
            ],
            global_config=GlobalConfig(snp_identifier="rsid"),
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
                    "POS": [10],
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
            self.assertEqual(called_args.snp_identifier, "rsid")
            self.assertEqual(result.reference_metadata["SNP"].tolist(), ["rs1"])
            self.assertIn("ldscore", result.output_paths)
            self.assertTrue(Path(result.output_paths["ldscore"]).exists())
            self.assertTrue(Path(result.output_paths["w_ld"]).exists())

    def test_run_ldscore_from_args_print_snps_filters_only_written_reference_table(self):
        fake_legacy_result = ldscore_workflow._LegacyChromResult(
            chrom="1",
            metadata=pd.DataFrame(
                {
                    "CHR": ["1", "1"],
                    "SNP": ["rs1", "rs2"],
                    "POS": [10, 20],
                    "CM": [0.1, 0.2],
                    "MAF": [0.2, 0.3],
                }
            ),
            ld_scores=np.array([[1.0, 2.0], [3.0, 4.0]], dtype=np.float32),
            w_ld=np.array([[5.0], [6.0]], dtype=np.float32),
            M=np.array([7.0, 8.0]),
            M_5_50=np.array([6.0, 7.0]),
            ldscore_columns=["base", "query"],
            baseline_columns=["base"],
            query_columns=["query"],
        )
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            baseline = tmpdir / "baseline.annot"
            baseline.write_text(
                "\n".join(
                    [
                        "CHR\tBP\tSNP\tCM\tbase",
                        "1\t10\trs1\t0.1\t1",
                        "1\t20\trs2\t0.2\t1",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )
            print_snps = tmpdir / "print_snps.txt"
            print_snps.write_text("rs2\n", encoding="utf-8")
            (tmpdir / "panel.bed").write_text("", encoding="utf-8")
            (tmpdir / "panel.bim").write_text("1 rs1 0.1 10 A G\n1 rs2 0.2 20 C T\n", encoding="utf-8")
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
                print_snps=str(print_snps),
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
            ), mock.patch.object(
                ldscore_workflow.ldscore_new,
                "compute_chrom_from_plink",
                return_value=fake_legacy_result,
            ):
                result = ldscore_workflow.run_ldscore_from_args(args)

            self.assertEqual(result.reference_metadata["SNP"].tolist(), ["rs1", "rs2"])
            self.assertEqual(result.reference_snps, {"rs1", "rs2"})
            ldscore_df = pd.read_csv(result.output_paths["ldscore"], sep="\t")
            self.assertEqual(ldscore_df["SNP"].tolist(), ["rs2"])
            weight_df = pd.read_csv(result.output_paths["w_ld"], sep="\t")
            self.assertEqual(weight_df["SNP"].tolist(), ["rs1", "rs2"])

    def test_namespace_from_configs_emits_string_paths(self):
        from ldsc._kernel.ref_panel import RefPanelSpec

        common = GlobalConfig(snp_identifier="rsid", genome_build="hg19")
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            r2_path = tmpdir / "r2" / "chr1.parquet"
            meta_path = tmpdir / "meta" / "chr1.tsv.gz"
            keep_path = tmpdir / "filters" / "samples.keep"
            r2_path.parent.mkdir(parents=True, exist_ok=True)
            meta_path.parent.mkdir(parents=True, exist_ok=True)
            keep_path.parent.mkdir(parents=True, exist_ok=True)
            r2_path.write_text("", encoding="utf-8")
            with gzip.open(meta_path, "wt", encoding="utf-8") as handle:
                handle.write("CHR\tBP\tSNP\tCM\tMAF\n1\t10\trs1\t0.1\t0.2\n")
            keep_path.write_text("iid1\n", encoding="utf-8")
            spec = RefPanelSpec(
                backend="parquet_r2",
                r2_table_paths=(r2_path,),
                maf_metadata_paths=(meta_path,),
            )
            ref_panel = SimpleNamespace(spec=spec)

            args = ldscore_workflow._namespace_from_configs(
                chrom="1",
                ref_panel=ref_panel,
                ldscore_config=ldscore_workflow.LDScoreConfig(
                    ld_wind_cm=1.0,
                    keep_individuals_path=keep_path,
                ),
                global_config=common,
            )

            self.assertIsInstance(args.r2_table, str)
            self.assertIsInstance(args.frqfile, str)
            self.assertEqual(args.r2_table, str(r2_path))
            self.assertEqual(args.frqfile, str(meta_path))
            self.assertEqual(args.keep, str(keep_path))
            self.assertEqual(args.genome_build, "hg19")

    def test_run_ldscore_from_args_resolves_glob_and_suite_tokens_in_workflow_layer(self):
        fake_legacy_result = ldscore_workflow._LegacyChromResult(
            chrom="1",
            metadata=pd.DataFrame(
                {
                    "CHR": ["1"],
                    "SNP": ["rs1"],
                    "POS": [10],
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
                baseline_annot=str(tmpdir / "baseline.@.annot.gz"),
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

    def test_run_ldscore_from_args_rejects_keep_in_parquet_mode(self):
        args = Namespace(
            out="out/example",
            output_dir=None,
            query_annot=None,
            query_annot_chr=None,
            baseline_annot="baseline.annot.gz",
            baseline_annot_chr=None,
            bfile=None,
            bfile_chr=None,
            r2_table="chr1.parquet",
            r2_table_chr=None,
            snp_identifier="rsid",
            genome_build=None,
            r2_bias_mode="unbiased",
            r2_sample_size=None,
            regression_snps=None,
            frqfile=None,
            frqfile_chr=None,
            keep="samples.keep",
            ld_wind_snps=10,
            ld_wind_kb=None,
            ld_wind_cm=None,
            maf=None,
            chunk_size=50,
            per_chr_output=False,
            yes_really=False,
            log_level="INFO",
        )

        with self.assertRaisesRegex(ValueError, "--keep.*PLINK"):
            ldscore_workflow.run_ldscore_from_args(args)

    @unittest.skipUnless(_HAS_BITARRAY, "bitarray is not installed")
    def test_ldscore_calculator_run_applies_keep_filter_by_fam_iid(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            prefix = self._copy_plink_fixture_with_distinct_fids(tmpdir)
            keep_path = tmpdir / "samples.keep"
            self._write_keep_file(keep_path, ["per0", "per1"])

            expected_metadata, expected_ld = self._expected_plink_result(prefix, keep_path, maf_min=None)
            bundle = self._build_annotation_bundle(prefix)
            common = GlobalConfig(snp_identifier="rsid")
            panel = PlinkRefPanel(common, RefPanelSpec(backend="plink", bfile_prefix=prefix))
            result = ldscore_workflow.LDScoreCalculator().run(
                annotation_bundle=bundle,
                ref_panel=panel,
                ldscore_config=LDScoreConfig(
                    ld_wind_snps=10,
                    keep_individuals_path=keep_path,
                    whole_chromosome_ok=True,
                ),
                global_config=common,
            )

            self.assertEqual(result.reference_metadata["SNP"].tolist(), expected_metadata["SNP"].tolist())
            np.testing.assert_allclose(result.reference_metadata["MAF"].to_numpy(), expected_metadata["MAF"].to_numpy())
            np.testing.assert_allclose(result.ld_scores["base"].to_numpy(), expected_ld)
            np.testing.assert_allclose(result.w_ld["L2"].to_numpy(), expected_ld)

    @unittest.skipUnless(_HAS_BITARRAY, "bitarray is not installed")
    def test_ldscore_calculator_run_filters_individuals_before_maf_in_plink_mode(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            prefix = self._copy_plink_fixture_with_distinct_fids(tmpdir)
            keep_path = tmpdir / "samples.keep"
            freq_path = tmpdir / "panel.frq"
            self._write_keep_file(keep_path, ["per0", "per1"])
            freq_path.write_text(
                "\n".join(
                    [
                        "CHR SNP BP CM MAF",
                        "1 rs_0 1 0 0.0",
                        "1 rs_1 2 0 0.0",
                        "1 rs_2 3 0 0.0",
                        "1 rs_3 4 0 0.0",
                        "1 rs_4 5 0 0.0",
                        "1 rs_5 6 0 0.0",
                        "1 rs_6 7 0 0.0",
                        "1 rs_7 8 0 0.0",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )

            expected_metadata, expected_ld = self._expected_plink_result(prefix, keep_path, maf_min=0.01)
            bundle = self._build_annotation_bundle(prefix)
            common = GlobalConfig(snp_identifier="rsid")
            panel = PlinkRefPanel(
                common,
                RefPanelSpec(
                    backend="plink",
                    bfile_prefix=prefix,
                    maf_metadata_paths=(freq_path,),
                ),
            )
            result = ldscore_workflow.LDScoreCalculator().run(
                annotation_bundle=bundle,
                ref_panel=panel,
                ldscore_config=LDScoreConfig(
                    ld_wind_snps=10,
                    maf_min=0.01,
                    keep_individuals_path=keep_path,
                    whole_chromosome_ok=True,
                ),
                global_config=common,
            )

            self.assertEqual(result.reference_metadata["SNP"].tolist(), expected_metadata["SNP"].tolist())
            np.testing.assert_allclose(result.reference_metadata["MAF"].to_numpy(), expected_metadata["MAF"].to_numpy())
            np.testing.assert_allclose(result.ld_scores["base"].to_numpy(), expected_ld)


@unittest.skipIf(kernel_ldscore is None, "ldscore kernel is not available")
class LDScoreParquetNormalizationTest(unittest.TestCase):
    def test_resolve_parquet_files_accepts_chromosome_specific_resolution(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            path1 = tmpdir / "r2.1.parquet"
            path2 = tmpdir / "r2.2.parquet"
            path1.write_text("", encoding="utf-8")
            path2.write_text("", encoding="utf-8")

            args = Namespace(r2_table=str(tmpdir / "r2.@.parquet"))

            self.assertEqual(
                kernel_ldscore.resolve_parquet_files(args, chrom="1"),
                [str(path1)],
            )

    def test_canonicalize_r2_pairs_renames_bp_aliases_to_pos_columns(self):
        df = pd.DataFrame(
            {
                "chr": ["1"],
                "rsID_1": ["rs2"],
                "rsID_2": ["rs1"],
                "hg38_bp1": [120],
                "hg38_bp2": [100],
                "hg19_bp_1": [20],
                "hg19_bp_2": [10],
                "hg38_Uniq_ID_1": ["1:120"],
                "hg38_Uniq_ID_2": ["1:100"],
                "hg19_Uniq_ID_1": ["1:20"],
                "hg19_Uniq_ID_2": ["1:10"],
                "R2": [0.4],
                "Dprime": [0.5],
                "+/-corr": ["+"],
            }
        )

        out = kernel_ldscore.canonicalize_r2_pairs(df, "GRCh37")

        self.assertEqual(out["chr"].tolist(), ["1"])
        self.assertEqual(out["pos_1"].tolist(), [10])
        self.assertEqual(out["pos_2"].tolist(), [20])
        self.assertIn("hg38_pos_1", out.columns)
        self.assertIn("hg38_pos_2", out.columns)
        self.assertIn("hg19_pos_1", out.columns)
        self.assertIn("hg19_pos_2", out.columns)
        self.assertNotIn("pair_chr", out.columns)

    def test_require_runtime_genome_build_accepts_aliases(self):
        self.assertEqual(kernel_ldscore._require_runtime_genome_build("hg37"), "hg19")
        self.assertEqual(kernel_ldscore._require_runtime_genome_build("GRCh37"), "hg19")
        self.assertEqual(kernel_ldscore._require_runtime_genome_build("GRCh38"), "hg38")

    def test_get_r2_build_columns_accepts_reduced_position_only_schema(self):
        columns = ["chr", "hg19_pos_1", "hg19_pos_2"]
        self.assertEqual(
            kernel_ldscore.get_r2_build_columns("hg19", columns),
            ("hg19_pos_1", "hg19_pos_2"),
        )

    @unittest.skipUnless(_HAS_PYARROW, "pyarrow is required for parquet reader coverage")
    def test_sorted_r2_block_reader_projects_actual_raw_schema_columns(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "raw_chr1.parquet"
            pd.DataFrame(
                {
                    "chr": ["1"],
                    "rsID_1": ["rs2"],
                    "rsID_2": ["rs1"],
                    "hg38_bp1": [120],
                    "hg38_bp2": [100],
                    "hg19_bp_1": [20],
                    "hg19_bp_2": [10],
                    "hg38_Uniq_ID_1": ["1:120"],
                    "hg38_Uniq_ID_2": ["1:100"],
                    "hg19_Uniq_ID_1": ["1:20"],
                    "hg19_Uniq_ID_2": ["1:10"],
                    "R2": [0.4],
                    "Dprime": [0.5],
                    "+/-corr": ["+"],
                }
            ).to_parquet(path, index=False)

            metadata = pd.DataFrame(
                {
                    "CHR": ["1", "1"],
                    "SNP": ["rs1", "rs2"],
                    "POS": [10, 20],
                    "CM": [0.1, 0.2],
                }
            )
            reader = kernel_ldscore.SortedR2BlockReader(
                paths=[str(path)],
                chrom="1",
                metadata=metadata,
                identifier_mode="rsID",
                r2_bias_mode="unbiased",
                r2_sample_size=None,
                genome_build="hg19",
            )

            matrix = reader.within_block_matrix(l_B=0, c=2)
            np.testing.assert_allclose(
                matrix,
                np.array([[1.0, 0.4], [0.4, 1.0]], dtype=np.float32),
            )

    @unittest.skipUnless(_HAS_PYARROW, "pyarrow is required for parquet reader coverage")
    def test_sorted_r2_block_reader_accepts_auto_genome_build(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "raw_chr1.parquet"
            pd.DataFrame(
                {
                    "chr": ["1"],
                    "rsID_1": ["rs2"],
                    "rsID_2": ["rs1"],
                    "hg38_bp1": [5120],
                    "hg38_bp2": [5100],
                    "hg19_bp_1": [120],
                    "hg19_bp_2": [100],
                    "hg38_Uniq_ID_1": ["1:5120"],
                    "hg38_Uniq_ID_2": ["1:5100"],
                    "hg19_Uniq_ID_1": ["1:120"],
                    "hg19_Uniq_ID_2": ["1:100"],
                    "R2": [0.4],
                    "Dprime": [0.5],
                    "+/-corr": ["+"],
                }
            ).to_parquet(path, index=False)

            metadata = pd.DataFrame(
                {
                    "CHR": ["1"] * 250,
                    "SNP": [f"rs{i}" for i in range(250)],
                    "POS": [100 + (idx * 10) for idx in range(250)],
                    "CM": [0.1] * 250,
                }
            )
            with mock.patch(
                "ldsc._kernel.ldscore.load_packaged_reference_table",
                return_value=pd.DataFrame(
                    {
                        "CHR": ["1"] * 250,
                        "hg19_POS": [100 + (idx * 10) for idx in range(250)],
                        "hg38_POS": [5100 + (idx * 10) for idx in range(250)],
                    }
                ),
            ):
                reader = kernel_ldscore.SortedR2BlockReader(
                    paths=[str(path)],
                    chrom="1",
                    metadata=metadata,
                    identifier_mode="chr_pos",
                    r2_bias_mode="unbiased",
                    r2_sample_size=None,
                    genome_build="auto",
                )

            self.assertEqual(reader.genome_build, "hg19")
