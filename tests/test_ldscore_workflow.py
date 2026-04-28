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
import warnings

import numpy as np
import pandas as pd

SRC = Path(__file__).resolve().parents[1] / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from ldsc.config import ConfigMismatchError, GlobalConfig

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
        baseline_table = pd.DataFrame(
            {
                "CHR": [chrom],
                "SNP": [f"rs{chrom}"],
                "BP": [bp],
                "regr_weight": [score + 2.0],
                "base": [score],
            }
        )
        query_table = pd.DataFrame(
            {
                "CHR": [chrom],
                "SNP": [f"rs{chrom}"],
                "BP": [bp],
                "query": [score + 1.0],
            }
        )
        return ldscore_workflow.ChromLDScoreResult(
            chrom=chrom,
            baseline_table=baseline_table,
            query_table=query_table,
            count_records=[
                {
                    "group": "baseline",
                    "column": "base",
                    "all_reference_snp_count": count,
                    "common_reference_snp_count_maf_gt_0_05": count - 1.0,
                },
                {
                    "group": "query",
                    "column": "query",
                    "all_reference_snp_count": count + 1.0,
                    "common_reference_snp_count_maf_gt_0_05": count,
                },
            ],
            baseline_columns=["base"],
            query_columns=["query"],
            ld_reference_snps=frozenset(),
            ld_regression_snps=frozenset({f"rs{chrom}"}),
            snp_count_totals={
                "all_reference_snp_counts": np.array([count, count + 1.0]),
                "common_reference_snp_counts_maf_gt_0_05": np.array([count - 1.0, count]),
            },
            config_snapshot=GlobalConfig(genome_build="hg38", snp_identifier="rsid"),
        )

    def make_annotation_bundle(self, chrom_rows: list[tuple[str, str, int]], genome_build: str = "hg38") -> AnnotationBundle:
        metadata = pd.DataFrame(chrom_rows, columns=["CHR", "SNP", "POS"])
        metadata["CM"] = 0.1
        baseline = pd.DataFrame({"base": np.ones(len(metadata), dtype=np.float32)})
        query = pd.DataFrame(index=metadata.index)
        chromosomes = metadata["CHR"].astype(str).drop_duplicates().tolist()
        return AnnotationBundle(
            metadata=metadata,
            baseline_annotations=baseline,
            query_annotations=query,
            baseline_columns=["base"],
            query_columns=[],
            chromosomes=chromosomes,
            source_summary={},
            config_snapshot=GlobalConfig(genome_build=genome_build, snp_identifier="rsid"),
        )

    def make_ref_panel_stub(self, *, backend: str, genome_build: str = "hg38", bfile_prefix: str | None = None, r2_table_paths: tuple[str, ...] = ()) -> SimpleNamespace:
        return SimpleNamespace(
            spec=SimpleNamespace(
                backend=backend,
                genome_build=genome_build,
                bfile_prefix=bfile_prefix,
                r2_table_paths=r2_table_paths,
                maf_metadata_paths=(),
                r2_bias_mode=None,
                sample_size=None,
                ref_panel_snps_path=None,
            )
        )

    def test_global_config_rejects_removed_ref_panel_snps_path(self):
        with self.assertRaises(TypeError):
            GlobalConfig(snp_identifier="rsid", ref_panel_snps_path="filters/reference.tsv.gz")

    def test_global_config_rejects_removed_regression_snps_path(self):
        with self.assertRaises(TypeError):
            GlobalConfig(snp_identifier="rsid", regression_snps_path="filters/hm3.tsv.gz")

    def test_ldscore_config_accepts_regression_snps_path(self):
        config = LDScoreConfig(ld_wind_snps=10, regression_snps_path="/path/to/snps.txt")
        self.assertEqual(config.regression_snps_path, "/path/to/snps.txt")

    def test_chrom_result_uses_split_table_shape(self):
        chrom_result = ldscore_workflow.ChromLDScoreResult(
            chrom="1",
            baseline_table=pd.DataFrame(
                {
                    "CHR": ["1"],
                    "SNP": ["rs1"],
                    "BP": [10],
                    "regr_weight": [3.0],
                    "base": [1.0],
                }
            ),
            query_table=pd.DataFrame(
                {
                    "CHR": ["1"],
                    "SNP": ["rs1"],
                    "BP": [10],
                    "query": [2.0],
                }
            ),
            count_records=[
                {"group": "baseline", "column": "base", "all_reference_snp_count": 10.0},
                {"group": "query", "column": "query", "all_reference_snp_count": 11.0},
            ],
            baseline_columns=["base"],
            query_columns=["query"],
            ld_reference_snps=frozenset(),
            ld_regression_snps=frozenset({"rs1"}),
            snp_count_totals={"all_reference_snp_counts": np.array([10.0, 11.0])},
            config_snapshot=GlobalConfig(genome_build="hg38", snp_identifier="rsid"),
        )
        self.assertFalse(hasattr(chrom_result, "reference_metadata"))
        self.assertFalse(hasattr(chrom_result, "w_ld"))
        self.assertFalse(hasattr(chrom_result, "ldscore_table"))
        self.assertEqual(chrom_result.baseline_table["regr_weight"].tolist(), [3.0])
        self.assertEqual(chrom_result.query_table["query"].tolist(), [2.0])

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
            config_snapshot=GlobalConfig(genome_build="hg38", snp_identifier="rsid"),
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
        self.assertFalse(hasattr(result, "reference_metadata"))
        self.assertEqual(result.baseline_table["CHR"].tolist(), ["1", "2"])
        self.assertEqual(result.baseline_table.columns.tolist(), ["CHR", "SNP", "BP", "regr_weight", "base"])
        self.assertEqual(result.query_table.columns.tolist(), ["CHR", "SNP", "BP", "query"])
        self.assertEqual(result.count_records[0]["all_reference_snp_count"], 30.0)
        self.assertEqual(result.count_records[1]["common_reference_snp_count_maf_gt_0_05"], 30.0)

    def test_build_parser_accepts_query_annot_bed(self):
        parser = ldscore_workflow.build_parser()
        args = parser.parse_args(
            [
                "--output-dir",
                "out/example",
                "--baseline-annot",
                "baseline.annot.gz",
                "--bfile",
                "panel",
                "--ld-wind-snps",
                "10",
                "--query-annot-bed",
                "query.bed",
            ]
        )
        self.assertEqual(args.query_annot_bed, "query.bed")

    def test_build_parser_rejects_query_annot_and_query_annot_bed_together(self):
        parser = ldscore_workflow.build_parser()
        with self.assertRaises(SystemExit):
            parser.parse_args(
                [
                    "--output-dir",
                    "out/example",
                    "--baseline-annot",
                    "baseline.annot.gz",
                    "--bfile",
                    "panel",
                    "--ld-wind-snps",
                    "10",
                    "--query-annot",
                    "query.annot.gz",
                    "--query-annot-bed",
                    "query.bed",
                ]
            )

    def test_run_rejects_annotation_bundle_snapshot_mismatch(self):
        calc = ldscore_workflow.LDScoreCalculator()
        annotation_bundle = AnnotationBundle(
            metadata=pd.DataFrame({"CHR": ["1"], "SNP": ["rs1"], "POS": [10], "CM": [0.1]}),
            baseline_annotations=pd.DataFrame({"base": [1.0]}),
            query_annotations=pd.DataFrame(index=pd.RangeIndex(1)),
            baseline_columns=["base"],
            query_columns=[],
            chromosomes=["1"],
            source_summary={},
            config_snapshot=GlobalConfig(genome_build="hg19", snp_identifier="rsid"),
        )
        ref_panel = SimpleNamespace(spec=SimpleNamespace(genome_build="hg19", backend="plink"))

        with mock.patch.object(calc, "compute_chromosome") as patched_compute:
            with self.assertRaisesRegex(ConfigMismatchError, "AnnotationBundle and LDScoreCalculator runtime config"):
                calc.run(
                    annotation_bundle=annotation_bundle,
                    ref_panel=ref_panel,
                    ldscore_config=LDScoreConfig(ld_wind_snps=10),
                    global_config=GlobalConfig(genome_build="hg38", snp_identifier="rsid"),
                )

        patched_compute.assert_not_called()

    def test_run_ldscore_from_args_writes_outputs(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            args = Namespace(
                output_dir=str(tmpdir / "ldscore_result"),
                query_annot=None,
                query_annot_chr=None,
                baseline_annot="baseline.annot.gz",
                baseline_annot_chr=None,
                bfile="panel",
                bfile_chr=None,
                r2_table=None,
                r2_table_chr=None,
                snp_identifier="rsid",
                genome_build=None,
                ref_panel_snps_path=None,
                regression_snps_path=None,
                r2_bias_mode=None,
                r2_sample_size=None,
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
            annotation_bundle = self.make_annotation_bundle([("1", "rs1", 10)])
            ref_panel = self.make_ref_panel_stub(backend="plink")
            with mock.patch.object(ldscore_workflow.ldscore_new, "validate_args") as validate_args, mock.patch(
                "ldsc._kernel.annotation.AnnotationBuilder.run",
                autospec=True,
                return_value=annotation_bundle,
            ), mock.patch(
                "ldsc._kernel.ref_panel.RefPanelLoader.load",
                autospec=True,
                return_value=ref_panel,
            ), mock.patch.object(
                ldscore_workflow.LDScoreCalculator,
                "compute_chromosome",
                autospec=True,
                return_value=self.make_chrom_result("1", 10, 1.0, 5.0),
            ):
                result = ldscore_workflow.run_ldscore_from_args(args)
            called_args = validate_args.call_args[0][0]
            self.assertEqual(called_args.output_dir, args.output_dir)
            self.assertEqual(called_args.snp_identifier, "rsid")
            self.assertEqual(result.baseline_table["SNP"].tolist(), ["rs1"])
            self.assertIn("manifest", result.output_paths)
            self.assertIn("baseline", result.output_paths)
            self.assertIn("query", result.output_paths)
            self.assertTrue(Path(result.output_paths["manifest"]).exists())
            self.assertFalse(list((tmpdir / "ldscore_result").glob("*.l2.ldscore.gz")))
            self.assertFalse(list((tmpdir / "ldscore_result").glob("*.M*")))

    def test_run_ldscore_from_args_loads_regression_snps_and_writes_filtered_ldscore(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            regression_snps_path = tmpdir / "regression_snps.txt"
            regression_snps_path.write_text("SNP\nrs2\n", encoding="utf-8")
            args = Namespace(
                output_dir=str(tmpdir / "ldscore_result"),
                query_annot=None,
                query_annot_chr=None,
                baseline_annot="baseline.annot.gz",
                baseline_annot_chr=None,
                bfile="panel",
                bfile_chr=None,
                r2_table=None,
                r2_table_chr=None,
                snp_identifier="rsid",
                genome_build=None,
                ref_panel_snps_path=None,
                regression_snps_path=str(regression_snps_path),
                r2_bias_mode=None,
                r2_sample_size=None,
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

            def _compute(_self, chrom, annotation_bundle, ref_panel, ldscore_config, global_config, regression_snps=None):
                self.assertEqual(chrom, "1")
                self.assertEqual(regression_snps, {"rs2"})
                return ldscore_workflow.ChromLDScoreResult(
                    chrom="1",
                    baseline_table=pd.DataFrame(
                        {
                            "CHR": ["1"],
                            "SNP": ["rs2"],
                            "BP": [20],
                            "regr_weight": [6.0],
                            "base": [3.0],
                        }
                    ),
                    query_table=pd.DataFrame(
                        {
                            "CHR": ["1"],
                            "SNP": ["rs2"],
                            "BP": [20],
                            "query": [4.0],
                        }
                    ),
                    count_records=[
                        {"group": "baseline", "column": "base", "all_reference_snp_count": 7.0, "common_reference_snp_count_maf_gt_0_05": 6.0},
                        {"group": "query", "column": "query", "all_reference_snp_count": 8.0, "common_reference_snp_count_maf_gt_0_05": 7.0},
                    ],
                    baseline_columns=["base"],
                    query_columns=["query"],
                    ld_reference_snps=frozenset(),
                    ld_regression_snps=frozenset({"rs2"}),
                    snp_count_totals={
                        "all_reference_snp_counts": np.array([7.0, 8.0]),
                        "common_reference_snp_counts_maf_gt_0_05": np.array([6.0, 7.0]),
                    },
                    config_snapshot=global_config,
                )

            annotation_bundle = self.make_annotation_bundle([("1", "rs1", 10), ("1", "rs2", 20)])
            ref_panel = self.make_ref_panel_stub(backend="plink")
            with mock.patch.object(ldscore_workflow.ldscore_new, "validate_args"), mock.patch(
                "ldsc._kernel.annotation.AnnotationBuilder.run",
                autospec=True,
                return_value=annotation_bundle,
            ), mock.patch(
                "ldsc._kernel.ref_panel.RefPanelLoader.load",
                autospec=True,
                return_value=ref_panel,
            ), mock.patch.object(
                ldscore_workflow.LDScoreCalculator,
                "compute_chromosome",
                autospec=True,
                side_effect=_compute,
            ):
                result = ldscore_workflow.run_ldscore_from_args(args)

            self.assertEqual(result.baseline_table["SNP"].tolist(), ["rs2"])
            self.assertEqual(result.ld_regression_snps, frozenset({"rs2"}))
            baseline_df = pd.read_parquet(result.output_paths["baseline"])
            self.assertEqual(baseline_df["SNP"].tolist(), ["rs2"])

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

    def test_run_ldscore_from_args_passes_path_tokens_to_builder_and_ref_panel_loader(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            args = Namespace(
                output_dir=str(tmpdir / "ldscore_result"),
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
                ref_panel_snps_path=None,
                regression_snps_path=None,
                r2_bias_mode=None,
                r2_sample_size=None,
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
            annotation_bundle = self.make_annotation_bundle([("1", "rs1", 10)])
            ref_panel = self.make_ref_panel_stub(backend="plink", bfile_prefix=str(tmpdir / "panel.@"))
            with mock.patch.object(ldscore_workflow.ldscore_new, "validate_args"), mock.patch(
                "ldsc._kernel.annotation.AnnotationBuilder.run",
                autospec=True,
                return_value=annotation_bundle,
            ) as patched_builder, mock.patch(
                "ldsc._kernel.ref_panel.RefPanelLoader.load",
                autospec=True,
                return_value=ref_panel,
            ) as patched_loader, mock.patch.object(
                ldscore_workflow.LDScoreCalculator,
                "compute_chromosome",
                autospec=True,
                return_value=self.make_chrom_result("1", 10, 1.0, 5.0),
            ):
                ldscore_workflow.run_ldscore_from_args(args)

            source_spec = patched_builder.call_args.args[1]
            self.assertEqual(source_spec.baseline_annot_paths, (str(tmpdir / "baseline.@.annot.gz"),))
            self.assertEqual(source_spec.query_annot_paths, ())
            self.assertEqual(source_spec.bed_paths, ())
            ref_spec = patched_loader.call_args.args[1]
            self.assertEqual(ref_spec.backend, "plink")
            self.assertEqual(str(ref_spec.bfile_prefix), str(tmpdir / "panel.@"))

    def test_run_ldscore_from_args_rejects_keep_in_parquet_mode(self):
        args = Namespace(
            output_dir="out/example",
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
            ref_panel_snps_path=None,
            regression_snps_path=None,
            r2_bias_mode="unbiased",
            r2_sample_size=None,
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

    def test_run_ldscore_from_args_warns_and_skips_empty_intersection_chromosome(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            args = Namespace(
                output_dir=str(tmpdir / "ldscore_result"),
                query_annot=None,
                query_annot_chr=None,
                baseline_annot=str(tmpdir / "baseline.@.annot.gz"),
                baseline_annot_chr=None,
                bfile=None,
                bfile_chr=None,
                r2_table=str(tmpdir / "r2.@.parquet"),
                r2_table_chr=None,
                snp_identifier="rsid",
                genome_build="hg19",
                ref_panel_snps_path=None,
                regression_snps_path=None,
                r2_bias_mode="unbiased",
                r2_sample_size=None,
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

            def _compute_side_effect(_self, chrom, annotation_bundle, ref_panel, ldscore_config, global_config, regression_snps=None):
                if chrom == "1":
                    raise ValueError("No retained annotation SNPs remain on chromosome 1 after parquet intersection.")
                return ldscore_workflow.ChromLDScoreResult(
                    chrom="22",
                    baseline_table=pd.DataFrame(
                        {
                            "CHR": ["22"],
                            "SNP": ["rs22"],
                            "BP": [220],
                            "regr_weight": [3.0],
                            "base": [2.0],
                        }
                    ),
                    query_table=None,
                    count_records=[
                        {"group": "baseline", "column": "base", "all_reference_snp_count": 5.0, "common_reference_snp_count_maf_gt_0_05": 4.0}
                    ],
                    baseline_columns=["base"],
                    query_columns=[],
                    ld_reference_snps=frozenset(),
                    ld_regression_snps=frozenset({"rs22"}),
                    snp_count_totals={
                        "all_reference_snp_counts": np.array([5.0]),
                        "common_reference_snp_counts_maf_gt_0_05": np.array([4.0]),
                    },
                    config_snapshot=global_config,
                )

            annotation_bundle = self.make_annotation_bundle([("1", "rs1", 10), ("22", "rs22", 220)], genome_build="hg19")
            ref_panel = self.make_ref_panel_stub(backend="parquet_r2", genome_build="hg19", r2_table_paths=(str(tmpdir / "r2.@.parquet"),))
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                with mock.patch.object(ldscore_workflow.ldscore_new, "validate_args"), mock.patch(
                    "ldsc._kernel.annotation.AnnotationBuilder.run",
                    autospec=True,
                    return_value=annotation_bundle,
                ), mock.patch(
                    "ldsc._kernel.ref_panel.RefPanelLoader.load",
                    autospec=True,
                    return_value=ref_panel,
                ), mock.patch.object(
                    ldscore_workflow.LDScoreCalculator,
                    "compute_chromosome",
                    autospec=True,
                    side_effect=_compute_side_effect,
                ):
                    result = ldscore_workflow.run_ldscore_from_args(args)

        self.assertEqual(result.baseline_table["CHR"].tolist(), ["22"])
        self.assertTrue(any("Skipping chromosome 1" in str(item.message) for item in caught))

    def test_compute_chromosome_filters_annotation_bundle_to_ref_panel_metadata_before_kernel_call(self):
        annotation_bundle = self.make_annotation_bundle(
            [("1", "rs1", 10), ("1", "rs2", 20), ("1", "rs3", 30)],
        )
        ref_panel = self.make_ref_panel_stub(backend="parquet_r2")
        ref_panel.load_metadata = mock.Mock(
            return_value=pd.DataFrame(
                {
                    "CHR": ["1", "1"],
                    "SNP": ["rs1", "rs3"],
                    "CM": [0.1, 0.3],
                    "POS": [10, 30],
                }
            )
        )

        def _compute_side_effect(chrom, bundle, args, regression_snps):
            self.assertEqual(chrom, "1")
            self.assertEqual(bundle.metadata["SNP"].tolist(), ["rs1", "rs3"])
            self.assertEqual(bundle.annotations["base"].tolist(), [1.0, 1.0])
            return ldscore_workflow._LegacyChromResult(
                chrom="1",
                metadata=pd.DataFrame(
                    {
                        "CHR": ["1", "1"],
                        "SNP": ["rs1", "rs3"],
                        "POS": [10, 30],
                        "CM": [0.1, 0.3],
                        "MAF": [0.2, 0.2],
                    }
                ),
                ld_scores=np.array([[1.0], [2.0]], dtype=np.float32),
                w_ld=np.array([[3.0], [4.0]], dtype=np.float32),
                M=np.array([2.0]),
                M_5_50=np.array([2.0]),
                ldscore_columns=["base"],
                baseline_columns=["base"],
                query_columns=[],
            )

        with mock.patch.object(
            ldscore_workflow.kernel_ldscore,
            "compute_chrom_from_parquet",
            side_effect=_compute_side_effect,
        ):
            result = ldscore_workflow.LDScoreCalculator().compute_chromosome(
                chrom="1",
                annotation_bundle=annotation_bundle,
                ref_panel=ref_panel,
                ldscore_config=LDScoreConfig(ld_wind_snps=10),
                global_config=GlobalConfig(snp_identifier="rsid"),
            )

        ref_panel.load_metadata.assert_called_once_with("1")
        self.assertEqual(result.baseline_table["SNP"].tolist(), ["rs1", "rs3"])
        self.assertEqual(result.count_records[0]["all_reference_snp_count"], 2.0)

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

            self.assertEqual(result.baseline_table["SNP"].tolist(), expected_metadata["SNP"].tolist())
            np.testing.assert_allclose(result.baseline_table["base"].to_numpy(), expected_ld)
            np.testing.assert_allclose(result.baseline_table["regr_weight"].to_numpy(), expected_ld)

    @unittest.skipUnless(_HAS_BITARRAY, "bitarray is not installed")
    def test_ldscore_calculator_run_applies_ref_panel_snp_restriction_before_plink_compute(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            prefix = self._copy_plink_fixture_with_distinct_fids(tmpdir)
            restrict_path = tmpdir / "restrict.snps"
            restrict_path.write_text("SNP\nrs_1\nrs_3\nrs_6\n", encoding="utf-8")

            bundle = self._build_annotation_bundle(prefix)
            common = GlobalConfig(snp_identifier="rsid")
            panel = PlinkRefPanel(
                common,
                RefPanelSpec(
                    backend="plink",
                    bfile_prefix=prefix,
                    ref_panel_snps_path=restrict_path,
                ),
            )
            result = ldscore_workflow.LDScoreCalculator().run(
                annotation_bundle=bundle,
                ref_panel=panel,
                ldscore_config=LDScoreConfig(
                    ld_wind_snps=10,
                    whole_chromosome_ok=True,
                ),
                global_config=common,
            )

            self.assertEqual(result.baseline_table["SNP"].tolist(), ["rs_1", "rs_3", "rs_6"])
            self.assertEqual(result.ld_regression_snps, frozenset({"rs_1", "rs_3", "rs_6"}))
            self.assertEqual(result.count_records[0]["all_reference_snp_count"], 3.0)

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

            self.assertEqual(result.baseline_table["SNP"].tolist(), expected_metadata["SNP"].tolist())
            np.testing.assert_allclose(result.baseline_table["base"].to_numpy(), expected_ld)

    def test_ldscore_calculator_run_warns_and_skips_empty_intersection_chromosome(self):
        annotation_bundle = AnnotationBundle(
            metadata=pd.DataFrame(
                {
                    "CHR": ["1", "22"],
                    "SNP": ["rs1", "rs22"],
                    "POS": [10, 220],
                    "CM": [0.1, 0.2],
                }
            ),
            baseline_annotations=pd.DataFrame({"base": [1.0, 1.0]}, dtype=np.float32),
            query_annotations=pd.DataFrame(index=pd.RangeIndex(2)),
            baseline_columns=["base"],
            query_columns=[],
            chromosomes=["1", "22"],
            source_summary={},
        )
        def _compute_side_effect(chrom, bundle, args, regression_snps):
            if chrom == "1":
                raise ValueError("No retained annotation SNPs remain on chromosome 1 after parquet intersection.")
            return ldscore_workflow._LegacyChromResult(
                chrom="22",
                metadata=pd.DataFrame(
                    {
                        "CHR": ["22"],
                        "SNP": ["rs22"],
                        "POS": [220],
                        "CM": [0.2],
                        "MAF": [0.3],
                    }
                ),
                ld_scores=np.array([[2.0]], dtype=np.float32),
                w_ld=np.array([[3.0]], dtype=np.float32),
                M=np.array([5.0]),
                M_5_50=np.array([4.0]),
                ldscore_columns=["base"],
                baseline_columns=["base"],
                query_columns=[],
            )

        calculator = ldscore_workflow.LDScoreCalculator()
        ref_panel = SimpleNamespace(spec=SimpleNamespace(backend="parquet_r2"))
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            with mock.patch.object(
                ldscore_workflow.kernel_ldscore,
                "compute_chrom_from_parquet",
                side_effect=_compute_side_effect,
            ):
                result = calculator.run(
                    annotation_bundle=annotation_bundle,
                    ref_panel=SimpleNamespace(
                        spec=SimpleNamespace(backend="parquet_r2"),
                        load_metadata=lambda chrom: pd.DataFrame(
                            {
                                "CHR": [chrom],
                                "SNP": [f"rs{chrom}"],
                                "POS": [10 if chrom == "1" else 220],
                                "CM": [0.1 if chrom == "1" else 0.2],
                            }
                        ),
                    ),
                    ldscore_config=LDScoreConfig(ld_wind_cm=1.0),
                    global_config=GlobalConfig(snp_identifier="rsid"),
                )

        self.assertEqual([chrom_result.chrom for chrom_result in result.chromosome_results], ["22"])
        self.assertEqual(result.baseline_table["CHR"].tolist(), ["22"])
        self.assertTrue(any("Skipping chromosome 1" in str(item.message) for item in caught))


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
