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
    from ldsc import AnnotationBundle, LDScoreConfig, PlinkRefPanel, RefPanelConfig
    from ldsc import ldscore_calculator as ldscore_workflow
    from ldsc._kernel import formats as kernel_formats
    from ldsc._kernel import ldscore as kernel_ldscore
    from ldsc._kernel.snp_identity import RestrictionIdentityKeys, empty_identity_drop_frame
except ImportError:
    AnnotationBundle = None
    LDScoreConfig = None
    PlinkRefPanel = None
    RefPanelConfig = None
    kernel_formats = None
    ldscore_workflow = None
    kernel_ldscore = None
    RestrictionIdentityKeys = None
    empty_identity_drop_frame = None


PLINK_FIXTURES = Path(__file__).resolve().parent / "fixtures" / "plink"
_HAS_PYARROW = importlib.util.find_spec("pyarrow") is not None
_HAS_BITARRAY = importlib.util.find_spec("bitarray") is not None


def _write_minimal_r2_parquet(path: Path, extra_meta: dict[bytes, bytes] | None = None) -> None:
    """Write a zero-row canonical R2 parquet table with optional schema metadata."""
    import pyarrow as pa
    import pyarrow.parquet as pq

    schema = pa.schema(
        [
            pa.field("CHR", pa.string()),
            pa.field("POS_1", pa.int64()),
            pa.field("POS_2", pa.int64()),
            pa.field("R2", pa.float32()),
            pa.field("SNP_1", pa.string()),
            pa.field("SNP_2", pa.string()),
        ]
    )
    meta = {
        b"ldsc:sorted_by_build": b"hg19",
        b"ldsc:row_group_size": b"50000",
    }
    if extra_meta:
        meta.update(extra_meta)
    table = pa.table({col: pa.array([], type=schema.field(col).type) for col in schema.names})
    table = table.replace_schema_metadata(meta)
    path.parent.mkdir(parents=True, exist_ok=True)
    pq.write_table(table, str(path))


class R2SchemaMetaReaderTest(unittest.TestCase):
    """Tests for R2 parquet schema metadata reader and resolver helpers."""

    @unittest.skipUnless(_HAS_PYARROW, "pyarrow is required for parquet schema coverage")
    def test_read_r2_schema_meta_returns_both_fields(self):
        from ldsc._kernel.ref_panel import _read_r2_schema_meta

        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "chr1_r2.parquet"
            _write_minimal_r2_parquet(
                path,
                {b"ldsc:n_samples": b"200", b"ldsc:r2_bias": b"unbiased"},
            )
            meta = _read_r2_schema_meta(str(path))

        self.assertEqual(meta.n_samples, 200)
        self.assertEqual(meta.r2_bias, "unbiased")

    @unittest.skipUnless(_HAS_PYARROW, "pyarrow is required for parquet schema coverage")
    def test_read_r2_schema_meta_legacy_file_returns_none_fields(self):
        from ldsc._kernel.ref_panel import _read_r2_schema_meta

        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "chr1_r2.parquet"
            _write_minimal_r2_parquet(path)
            meta = _read_r2_schema_meta(str(path))

        self.assertIsNone(meta.n_samples)
        self.assertIsNone(meta.r2_bias)

    @unittest.skipUnless(_HAS_PYARROW, "pyarrow is required for parquet schema coverage")
    def test_read_r2_schema_meta_warns_for_n_samples_without_bias(self):
        from ldsc._kernel.ref_panel import _read_r2_schema_meta

        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "chr1_r2.parquet"
            _write_minimal_r2_parquet(path, {b"ldsc:n_samples": b"200"})
            with self.assertLogs("LDSC", level="WARNING"):
                meta = _read_r2_schema_meta(str(path))

        self.assertEqual(meta.n_samples, 200)
        self.assertEqual(meta.r2_bias, "raw")

    def test_resolve_unbiased_no_user_n_returns_unbiased(self):
        from ldsc._kernel.ref_panel import _R2SchemaMeta, _resolve_r2_bias_from_meta

        meta = _R2SchemaMeta(n_samples=100, r2_bias="unbiased")
        mode, n = _resolve_r2_bias_from_meta(None, None, meta)

        self.assertEqual(mode, "unbiased")
        self.assertIsNone(n)

    def test_resolve_unbiased_with_user_n_warns_and_returns_unbiased(self):
        from ldsc._kernel.ref_panel import _R2SchemaMeta, _resolve_r2_bias_from_meta

        meta = _R2SchemaMeta(n_samples=100, r2_bias="unbiased")
        with self.assertLogs("LDSC", level="WARNING"):
            mode, n = _resolve_r2_bias_from_meta(None, 999.0, meta)

        self.assertEqual(mode, "unbiased")
        self.assertIsNone(n)

    def test_resolve_raw_no_user_n_autofills_stored_n(self):
        from ldsc._kernel.ref_panel import _R2SchemaMeta, _resolve_r2_bias_from_meta

        meta = _R2SchemaMeta(n_samples=150, r2_bias="raw")
        mode, n = _resolve_r2_bias_from_meta(None, None, meta)

        self.assertEqual(mode, "raw")
        self.assertEqual(n, 150.0)

    def test_resolve_raw_user_n_overrides_stored_n(self):
        from ldsc._kernel.ref_panel import _R2SchemaMeta, _resolve_r2_bias_from_meta

        meta = _R2SchemaMeta(n_samples=150, r2_bias="raw")
        mode, n = _resolve_r2_bias_from_meta(None, 999.0, meta)

        self.assertEqual(mode, "raw")
        self.assertEqual(n, 999.0)

    def test_resolve_absent_metadata_defaults_to_unbiased(self):
        from ldsc._kernel.ref_panel import _R2SchemaMeta, _resolve_r2_bias_from_meta

        meta = _R2SchemaMeta(n_samples=None, r2_bias=None)
        mode, n = _resolve_r2_bias_from_meta(None, None, meta)

        self.assertEqual(mode, "unbiased")
        self.assertIsNone(n)

    def test_resolve_absent_metadata_raw_mode_with_user_n(self):
        from ldsc._kernel.ref_panel import _R2SchemaMeta, _resolve_r2_bias_from_meta

        meta = _R2SchemaMeta(n_samples=None, r2_bias=None)
        mode, n = _resolve_r2_bias_from_meta("raw", 300.0, meta)

        self.assertEqual(mode, "raw")
        self.assertEqual(n, 300.0)


@unittest.skipIf(kernel_ldscore is None, "ldscore kernel is not available")
class R2AutoLoadCLITest(unittest.TestCase):
    @unittest.skipUnless(_HAS_PYARROW, "pyarrow is required for parquet schema coverage")
    def test_cli_autofills_unbiased_from_schema_when_mode_is_none(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "hg19" / "chr1_r2.parquet"
            _write_minimal_r2_parquet(
                path,
                {b"ldsc:n_samples": b"200", b"ldsc:r2_bias": b"unbiased"},
            )
            args = kernel_ldscore.build_parser().parse_args(
                [
                    "--r2-table",
                    str(Path(tmpdir) / "hg19"),
                    "--snp-identifier",
                    "rsid",
                    "--baseline-annot",
                    "fake",
                    "--out",
                    "fake",
                    "--ld-wind-kb",
                    "1",
                ]
            )
            args.r2_bias_mode = None
            args.r2_sample_size = None

            kernel_ldscore.validate_args(args)

        self.assertEqual(args.r2_bias_mode, "unbiased")
        self.assertIsNone(args.r2_sample_size)

    @unittest.skipUnless(_HAS_PYARROW, "pyarrow is required for parquet schema coverage")
    def test_cli_autofills_raw_and_n_from_schema(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "hg19" / "chr1_r2.parquet"
            _write_minimal_r2_parquet(
                path,
                {b"ldsc:n_samples": b"150", b"ldsc:r2_bias": b"raw"},
            )
            args = kernel_ldscore.build_parser().parse_args(
                [
                    "--r2-table",
                    str(Path(tmpdir) / "hg19"),
                    "--snp-identifier",
                    "rsid",
                    "--baseline-annot",
                    "fake",
                    "--out",
                    "fake",
                    "--ld-wind-kb",
                    "1",
                ]
            )
            args.r2_bias_mode = None
            args.r2_sample_size = None

            kernel_ldscore.validate_args(args)

        self.assertEqual(args.r2_bias_mode, "raw")
        self.assertAlmostEqual(args.r2_sample_size, 150.0)

    @unittest.skipUnless(_HAS_PYARROW, "pyarrow is required for parquet schema coverage")
    def test_parquet_panel_autofills_raw_and_n_from_schema(self):
        from ldsc._kernel.ref_panel import ParquetR2RefPanel

        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "hg19" / "chr1_r2.parquet"
            _write_minimal_r2_parquet(
                path,
                {
                    b"ldsc:n_samples": b"200",
                    b"ldsc:r2_bias": b"raw",
                    b"ldsc:schema_version": b"1",
                    b"ldsc:artifact_type": b"ref_panel_r2",
                    b"ldsc:snp_identifier": b"rsid",
                    b"ldsc:genome_build": b"hg19",
                },
            )
            spec = RefPanelConfig(
                backend="parquet_r2",
                r2_dir=str(Path(tmpdir) / "hg19"),
                r2_bias_mode=None,
                r2_sample_size=None,
            )
            panel = ParquetR2RefPanel(
                GlobalConfig(snp_identifier="rsid"),
                spec,
            )
            metadata = pd.DataFrame(
                {
                    "CHR": ["1"],
                    "POS": [100],
                    "SNP": ["rs1"],
                    "CM": [0.0],
                    "MAF": [0.3],
                }
            )

            reader = panel.build_reader("1", metadata=metadata)

        self.assertEqual(reader.r2_bias_mode, "raw")
        self.assertAlmostEqual(reader.r2_sample_size, 200.0)


@unittest.skipIf(ldscore_workflow is None, "ldscore_workflow module is not available")
class LDScoreWorkflowTest(unittest.TestCase):
    def test_build_parser_genome_build_help_documents_chr_pos_requirement(self):
        help_text = ldscore_workflow.build_parser().format_help()

        self.assertIn("Required when", help_text)
        self.assertIn("Not used for rsid-family modes", help_text)

    def test_build_parser_accepts_r2_dir(self):
        args = ldscore_workflow.build_parser().parse_args(
            [
                "--output-dir",
                "out",
                "--r2-dir",
                "panel/hg38",
                "--r2-bias-mode",
                "unbiased",
            ]
        )

        self.assertEqual(args.r2_dir, "panel/hg38")

    def test_build_parser_defaults_r2_bias_mode_to_unbiased(self):
        args = ldscore_workflow.build_parser().parse_args(["--output-dir", "out", "--r2-dir", "panel/hg38"])

        self.assertEqual(args.r2_bias_mode, "unbiased")

    def test_format_ldscore_start_message_describes_synthetic_base(self):
        bundle = SimpleNamespace(
            baseline_columns=["base"],
            query_columns=[],
            source_summary={"baseline": "synthetic all-ones base annotation from retained reference-panel metadata"},
        )

        message = ldscore_workflow._format_ldscore_start_message(bundle, 22)

        self.assertEqual(
            message,
            "Starting LD-score calculation for 22 chromosomes with synthetic base annotation and no query annotations.",
        )

    def test_format_ldscore_start_message_describes_partitioned_columns(self):
        bundle = SimpleNamespace(
            baseline_columns=["base_a", "base_b"],
            query_columns=["query"],
            source_summary={"baseline": "baseline.@.annot.gz"},
        )

        message = ldscore_workflow._format_ldscore_start_message(bundle, 22)

        self.assertEqual(
            message,
            "Starting LD-score calculation for 22 chromosomes with 2 baseline columns and 1 query columns.",
        )

    def test_validate_run_args_accepts_omitted_r2_bias_mode_as_unbiased(self):
        args = Namespace(
            output_dir="out",
            query_annot_sources=None,
            query_annot_bed_sources=None,
            baseline_annot_sources="baseline.annot.gz",
            plink_prefix=None,
            bfile=None,
            bfile_chr=None,
            r2_dir="panel/hg38",
            r2_table_chr=None,
            snp_identifier="rsid",
            genome_build=None,
            keep_indivs_file=None,
            r2_bias_mode=None,
            r2_sample_size=None,
            ld_wind_snps=10,
            ld_wind_kb=None,
            ld_wind_cm=None,
            maf_min=None,
            common_maf_min=0.05,
            snp_batch_size=50,
        )

        ldscore_workflow._validate_run_args(args)
        self.assertEqual(args.r2_bias_mode, "unbiased")

    def test_validate_run_args_requires_sample_size_for_raw_r2(self):
        args = Namespace(
            output_dir="out",
            query_annot_sources=None,
            query_annot_bed_sources=None,
            baseline_annot_sources="baseline.annot.gz",
            plink_prefix=None,
            bfile=None,
            bfile_chr=None,
            r2_dir="panel/hg38",
            r2_table_chr=None,
            snp_identifier="rsid",
            genome_build=None,
            keep_indivs_file=None,
            r2_bias_mode="raw",
            r2_sample_size=None,
            ld_wind_snps=10,
            ld_wind_kb=None,
            ld_wind_cm=None,
            maf_min=None,
            common_maf_min=0.05,
            snp_batch_size=50,
        )

        with self.assertRaisesRegex(ValueError, "--r2-sample-size is required"):
            ldscore_workflow._validate_run_args(args)

    def test_kernel_build_parser_defaults_r2_bias_mode_to_unbiased(self):
        args = kernel_ldscore.build_parser().parse_args(
            [
                "--out",
                "out",
                "--baseline-annot",
                "baseline.annot.gz",
                "--r2-table",
                "panel.@.parquet",
            ]
        )

        self.assertEqual(args.r2_bias_mode, "unbiased")

    def test_kernel_build_parser_accepts_snp_batch_size_and_rejects_chunk_size(self):
        parser = kernel_ldscore.build_parser()
        args = parser.parse_args(
            [
                "--out",
                "out",
                "--baseline-annot",
                "baseline.annot.gz",
                "--r2-table",
                "panel.@.parquet",
                "--snp-batch-size",
                "64",
            ]
        )
        self.assertEqual(args.snp_batch_size, 64)
        with self.assertRaises(SystemExit):
            parser.parse_args(
                [
                    "--out",
                    "out",
                    "--baseline-annot",
                    "baseline.annot.gz",
                    "--r2-table",
                    "panel.@.parquet",
                    "--chunk-size",
                    "64",
                ]
            )

    def test_kernel_validate_args_accepts_omitted_r2_bias_mode_as_unbiased(self):
        args = Namespace(
            out="out",
            query_annot=None,
            baseline_annot="baseline.annot.gz",
            bfile=None,
            r2_table="panel.@.parquet",
            snp_identifier="rsid",
            genome_build=None,
            r2_bias_mode=None,
            r2_sample_size=None,
            regression_snps_file=None,
            frqfile=None,
            query_annot_chr=None,
            baseline_annot_chr=None,
            bfile_chr=None,
            r2_table_chr=None,
            frqfile_chr=None,
            keep=None,
            ld_wind_snps=10,
            ld_wind_kb=None,
            ld_wind_cm=None,
            maf_min=None,
            common_maf_min=0.05,
            snp_batch_size=50,
            per_chr_output=False,
            yes_really=False,
            log_level="INFO",
        )

        kernel_ldscore.validate_args(args)
        self.assertEqual(args.r2_bias_mode, "unbiased")

    def test_build_parser_help_exposes_only_r2_dir_for_parquet_input(self):
        help_text = ldscore_workflow.build_parser().format_help()

        self.assertIn("--r2-dir", help_text)
        self.assertNotIn("--r2-ref-panel-dir", help_text)
        self.assertNotIn("--ref-panel-dir", help_text)
        self.assertNotIn("--r2-sources", help_text)
        self.assertNotIn("--metadata-sources", help_text)

    def test_build_parser_rejects_removed_parquet_input_flags(self):
        parser = ldscore_workflow.build_parser()
        for flag in ("--r2-ref-panel-dir", "--ref-panel-dir", "--r2-sources", "--metadata-sources"):
            with self.subTest(flag=flag):
                with self.assertRaises(SystemExit):
                    parser.parse_args(["--output-dir", "out", flag, "value"])

    def test_run_ldscore_rejects_removed_parquet_input_kwargs(self):
        for name in ("r2_ref_panel_dir", "ref_panel_dir", "r2_sources", "metadata_sources", "chunk_size"):
            with self.subTest(name=name):
                with self.assertRaisesRegex(ValueError, "no longer accepts"):
                    ldscore_workflow.run_ldscore(output_dir="out", **{name: "value"})

    def test_normalize_run_args_chr_pos_still_requires_genome_build(self):
        args = Namespace(
            output_dir="out",
            query_annot_sources=None,
            baseline_annot_sources=None,
            plink_prefix="plink/panel.@",
            snp_identifier="chr_pos",
            genome_build=None,
            keep_indivs_file=None,
            ref_panel_snps_file=None,
            regression_snps_file=None,
            log_level="INFO",
        )

        with self.assertRaisesRegex(ValueError, "genome_build is required"):
            ldscore_workflow._normalize_run_args(args)

    def test_auto_genome_build_does_not_infer_from_r2_directory_name_without_parquet_metadata(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            r2_dir = Path(tmpdir) / "hg38"
            r2_dir.mkdir()
            args = Namespace(
                output_dir="out",
                query_annot_sources=None,
                baseline_annot_sources=None,
                plink_prefix=None,
                r2_dir=str(r2_dir),
                snp_identifier="chr_pos",
                genome_build="auto",
                keep_indivs_file=None,
                ref_panel_snps_file=None,
                regression_snps_file=None,
                log_level="INFO",
            )

            with self.assertRaisesRegex(ValueError, "Cannot infer --genome-build"):
                ldscore_workflow._normalize_run_args(args)

    @unittest.skipUnless(_HAS_PYARROW, "pyarrow is required for parquet schema coverage")
    def test_auto_genome_build_uses_r2_schema_metadata_not_directory_name(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            r2_dir = Path(tmpdir) / "hg38"
            _write_minimal_r2_parquet(
                r2_dir / "chr1_r2.parquet",
                {b"ldsc:sorted_by_build": b"hg19"},
            )
            args = Namespace(
                output_dir="out",
                query_annot_sources=None,
                baseline_annot_sources=None,
                plink_prefix=None,
                r2_dir=str(r2_dir),
                snp_identifier="chr_pos",
                genome_build="auto",
                keep_indivs_file=None,
                ref_panel_snps_file=None,
                regression_snps_file=None,
                log_level="INFO",
            )

            normalized_args, global_config = ldscore_workflow._normalize_run_args(args)

        self.assertEqual(normalized_args.genome_build, "hg19")
        self.assertEqual(global_config.genome_build, "hg19")

    @unittest.skipUnless(_HAS_PYARROW, "pyarrow is required for parquet schema coverage")
    def test_auto_genome_build_uses_single_child_r2_schema_metadata_not_child_name(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            r2_parent = Path(tmpdir) / "panel"
            _write_minimal_r2_parquet(
                r2_parent / "hg38" / "chr1_r2.parquet",
                {b"ldsc:sorted_by_build": b"hg19"},
            )
            args = Namespace(
                output_dir="out",
                query_annot_sources=None,
                baseline_annot_sources=None,
                plink_prefix=None,
                r2_dir=str(r2_parent),
                snp_identifier="chr_pos",
                genome_build="auto",
                keep_indivs_file=None,
                ref_panel_snps_file=None,
                regression_snps_file=None,
                log_level="INFO",
            )

            normalized_args, global_config = ldscore_workflow._normalize_run_args(args)

        self.assertEqual(normalized_args.genome_build, "hg19")
        self.assertEqual(global_config.genome_build, "hg19")

    @unittest.skipUnless(_HAS_PYARROW, "pyarrow is required for parquet schema coverage")
    def test_auto_genome_build_rejects_conflicting_r2_schema_metadata(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            r2_parent = Path(tmpdir) / "panel"
            _write_minimal_r2_parquet(
                r2_parent / "hg19" / "chr1_r2.parquet",
                {b"ldsc:sorted_by_build": b"hg19"},
            )
            _write_minimal_r2_parquet(
                r2_parent / "hg38" / "chr1_r2.parquet",
                {b"ldsc:sorted_by_build": b"hg38"},
            )
            args = Namespace(
                output_dir="out",
                query_annot_sources=None,
                baseline_annot_sources=None,
                plink_prefix=None,
                r2_dir=str(r2_parent),
                snp_identifier="chr_pos",
                genome_build="auto",
                keep_indivs_file=None,
                ref_panel_snps_file=None,
                regression_snps_file=None,
                log_level="INFO",
            )

            with self.assertRaisesRegex(ValueError, "Conflicting R2 parquet genome-build metadata"):
                ldscore_workflow._normalize_run_args(args)

    def make_chrom_result(self, chrom: str, bp: int, score: float, count: float):
        baseline_table = pd.DataFrame(
            {
                "CHR": [chrom],
                "SNP": [f"rs{chrom}"],
                "POS": [bp],
                "regression_ld_scores": [score + 2.0],
                "base": [score],
            }
        )
        query_table = pd.DataFrame(
            {
                "CHR": [chrom],
                "SNP": [f"rs{chrom}"],
                "POS": [bp],
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
                    "common_reference_snp_count": count - 1.0,
                },
                {
                    "group": "query",
                    "column": "query",
                    "all_reference_snp_count": count + 1.0,
                    "common_reference_snp_count": count,
                },
            ],
            baseline_columns=["base"],
            query_columns=["query"],
            ld_reference_snps=frozenset(),
            ld_regression_snps=frozenset({f"rs{chrom}"}),
            snp_count_totals={
                "all_reference_snp_counts": np.array([count, count + 1.0]),
                "common_reference_snp_counts": np.array([count - 1.0, count]),
            },
            config_snapshot=GlobalConfig(snp_identifier="rsid"),
        )

    def make_baseline_only_chrom_result(self, chrom: str, bp: int, score: float, count: float):
        return ldscore_workflow.ChromLDScoreResult(
            chrom=chrom,
            baseline_table=pd.DataFrame(
                {
                    "CHR": [chrom],
                    "SNP": [f"rs{chrom}"],
                    "POS": [bp],
                    "regression_ld_scores": [score + 2.0],
                    "base": [score],
                }
            ),
            query_table=None,
            count_records=[
                {
                    "group": "baseline",
                    "column": "base",
                    "all_reference_snp_count": count,
                    "common_reference_snp_count": count - 1.0,
                }
            ],
            baseline_columns=["base"],
            query_columns=[],
            ld_reference_snps=frozenset(),
            ld_regression_snps=frozenset({f"rs{chrom}"}),
            snp_count_totals={
                "all_reference_snp_counts": np.array([count]),
                "common_reference_snp_counts": np.array([count - 1.0]),
            },
            config_snapshot=GlobalConfig(snp_identifier="rsid"),
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
            config_snapshot=GlobalConfig(snp_identifier="rsid"),
        )

    def make_ref_panel_stub(
        self,
        *,
        backend: str,
        genome_build: str = "hg38",
        plink_prefix: str | None = None,
        metadata: pd.DataFrame | None = None,
    ) -> SimpleNamespace:
        metadata_by_chrom = {}
        if metadata is not None:
            for chrom, chrom_metadata in metadata.groupby(metadata["CHR"].astype(str), sort=False):
                metadata_by_chrom[str(chrom)] = chrom_metadata.reset_index(drop=True)
        panel = SimpleNamespace(
            spec=SimpleNamespace(
                backend=backend,
                genome_build=genome_build,
                plink_prefix=plink_prefix,
                r2_bias_mode=None,
                sample_size=None,
                ref_panel_snps_file=None,
            )
        )
        if metadata is not None:
            panel.available_chromosomes = lambda: list(metadata_by_chrom)
            panel.load_metadata = lambda chrom: metadata_by_chrom[str(chrom)].copy()
        return panel

    def test_global_config_rejects_removed_ref_panel_snps_file(self):
        with self.assertRaises(TypeError):
            GlobalConfig(snp_identifier="rsid", ref_panel_snps_file="filters/reference.tsv.gz")

    def test_global_config_rejects_removed_regression_snps_file(self):
        with self.assertRaises(TypeError):
            GlobalConfig(snp_identifier="rsid", regression_snps_file="filters/hm3.tsv.gz")

    def test_ldscore_config_accepts_regression_snps_file(self):
        config = LDScoreConfig(ld_wind_snps=10, regression_snps_file="/path/to/snps.txt")
        self.assertEqual(config.regression_snps_file, "/path/to/snps.txt")

    def test_build_parser_accepts_hm3_ref_panel_and_regression_flags(self):
        parser = ldscore_workflow.build_parser()

        args = parser.parse_args(
            [
                "--output-dir",
                "out",
                "--baseline-annot-sources",
                "baseline.annot.gz",
                "--plink-prefix",
                "panel",
                "--ld-wind-snps",
                "10",
                "--use-hm3-ref-panel-snps",
                "--use-hm3-regression-snps",
            ]
        )

        self.assertTrue(args.use_hm3_ref_panel_snps)
        self.assertTrue(args.use_hm3_regression_snps)

    def test_normalize_run_args_rejects_hm3_explicit_file_conflicts(self):
        parser = ldscore_workflow.build_parser()
        args = parser.parse_args(
            [
                "--output-dir",
                "out",
                "--baseline-annot-sources",
                "baseline.annot.gz",
                "--plink-prefix",
                "panel",
                "--ld-wind-snps",
                "10",
                "--snp-identifier",
                "rsid",
                "--use-hm3-regression-snps",
                "--regression-snps-file",
                "custom.tsv",
            ]
        )

        with self.assertRaisesRegex(ValueError, "regression_snps_file.*use_hm3_regression_snps"):
            ldscore_workflow._normalize_run_args(args)

    def test_ldscore_config_from_args_preserves_hm3_regression_flag(self):
        parser = ldscore_workflow.build_parser()
        args = parser.parse_args(
            [
                "--output-dir",
                "out",
                "--baseline-annot-sources",
                "baseline.annot.gz",
                "--plink-prefix",
                "panel",
                "--ld-wind-snps",
                "10",
                "--snp-identifier",
                "rsid",
                "--use-hm3-regression-snps",
            ]
        )
        normalized, _global_config = ldscore_workflow._normalize_run_args(args)

        config = ldscore_workflow._ldscore_config_from_args(normalized)

        self.assertTrue(config.use_hm3_regression_snps)

    def test_chrom_result_uses_split_table_shape(self):
        chrom_result = ldscore_workflow.ChromLDScoreResult(
            chrom="1",
            baseline_table=pd.DataFrame(
                {
                    "CHR": ["1"],
                    "SNP": ["rs1"],
                    "POS": [10],
                    "regression_ld_scores": [3.0],
                    "base": [1.0],
                }
            ),
            query_table=pd.DataFrame(
                {
                    "CHR": ["1"],
                    "SNP": ["rs1"],
                    "POS": [10],
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
            config_snapshot=GlobalConfig(snp_identifier="rsid"),
        )
        self.assertFalse(hasattr(chrom_result, "reference_metadata"))
        self.assertFalse(hasattr(chrom_result, "w_ld"))
        self.assertFalse(hasattr(chrom_result, "ldscore_table"))
        self.assertEqual(chrom_result.baseline_table["regression_ld_scores"].tolist(), [3.0])
        self.assertEqual(chrom_result.query_table["query"].tolist(), [2.0])

    def test_wrap_legacy_result_derives_allele_aware_regression_ids_before_split(self):
        legacy_result = ldscore_workflow._LegacyChromResult(
            chrom="1",
            metadata=pd.DataFrame(
                {
                    "CHR": ["1", "1"],
                    "SNP": ["rs1", "rs2"],
                    "POS": [10, 20],
                    "A1": ["A", "A"],
                    "A2": ["C", "G"],
                    "CM": [0.1, 0.2],
                    "MAF": [0.2, 0.3],
                }
            ),
            ld_scores=np.array([[1.0], [2.0]], dtype=np.float32),
            w_ld=np.array([[3.0], [4.0]], dtype=np.float32),
            M=np.array([2.0]),
            M_5_50=None,
            ldscore_columns=["base"],
            baseline_columns=["base"],
            query_columns=[],
        )

        result = ldscore_workflow.LDScoreCalculator()._wrap_legacy_chrom_result(
            legacy_result,
            global_config=GlobalConfig(snp_identifier="chr_pos_allele_aware"),
        )

        self.assertEqual(result.baseline_table.columns.tolist(), ["CHR", "SNP", "POS", "regression_ld_scores", "base"])
        self.assertEqual(result.ld_regression_snps, frozenset({"1:10:A:C", "1:20:A:G"}))

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
            config_snapshot=GlobalConfig(snp_identifier="rsid"),
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
        self.assertEqual(result.baseline_table.columns.tolist(), ["CHR", "SNP", "POS", "regression_ld_scores", "base"])
        self.assertEqual(result.query_table.columns.tolist(), ["CHR", "SNP", "POS", "query"])
        self.assertEqual(result.count_records[0]["all_reference_snp_count"], 30.0)
        self.assertEqual(result.count_records[1]["common_reference_snp_count"], 30.0)

    def test_build_parser_accepts_query_annot_bed(self):
        parser = ldscore_workflow.build_parser()
        args = parser.parse_args(
            [
                "--output-dir",
                "out/example",
                "--baseline-annot-sources",
                "baseline.annot.gz",
                "--plink-prefix",
                "panel",
                "--ld-wind-snps",
                "10",
                "--query-annot-bed-sources",
                "query.bed",
            ]
        )
        self.assertEqual(args.query_annot_bed_sources, "query.bed")

    def test_build_parser_defaults_snp_batch_size_to_128(self):
        parser = ldscore_workflow.build_parser()
        self.assertEqual(parser.get_default("snp_batch_size"), 128)

    def test_build_parser_accepts_snp_batch_size_and_rejects_chunk_size(self):
        parser = ldscore_workflow.build_parser()
        args = parser.parse_args(["--output-dir", "out", "--snp-batch-size", "64"])
        self.assertEqual(args.snp_batch_size, 64)
        with self.assertRaises(SystemExit):
            parser.parse_args(["--output-dir", "out", "--chunk-size", "64"])

    def test_build_parser_rejects_query_annot_and_query_annot_bed_together(self):
        parser = ldscore_workflow.build_parser()
        with self.assertRaises(SystemExit):
            parser.parse_args(
                [
                    "--output-dir",
                    "out/example",
                    "--baseline-annot-sources",
                    "baseline.annot.gz",
                    "--plink-prefix",
                    "panel",
                    "--ld-wind-snps",
                    "10",
                    "--query-annot-sources",
                    "query.annot.gz",
                    "--query-annot-bed-sources",
                    "query.bed",
                ]
            )

    def test_build_parser_accepts_maf_min_and_common_maf_min_but_rejects_old_maf(self):
        parser = ldscore_workflow.build_parser()
        args = parser.parse_args(
            [
                "--output-dir",
                "out/example",
                "--baseline-annot-sources",
                "baseline.annot.gz",
                "--plink-prefix",
                "panel",
                "--ld-wind-snps",
                "10",
                "--maf-min",
                "0.01",
                "--common-maf-min",
                "0.2",
            ]
        )
        self.assertEqual(args.maf_min, 0.01)
        self.assertEqual(args.common_maf_min, 0.2)
        with self.assertRaises(SystemExit):
            parser.parse_args(
                [
                    "--output-dir",
                    "out/example",
                    "--baseline-annot-sources",
                    "baseline.annot.gz",
                    "--plink-prefix",
                    "panel",
                    "--ld-wind-snps",
                    "10",
                    "--maf",
                    "0.01",
                ]
            )

    def test_compute_counts_uses_inclusive_common_maf_min(self):
        metadata = pd.DataFrame({"MAF": [0.1, 0.2, 0.3]})
        annotations = pd.DataFrame({"base": [1.0, 2.0, 3.0], "query": [0.0, 4.0, 5.0]})

        all_counts, common_counts = kernel_ldscore.compute_counts(
            metadata,
            annotations,
            common_maf_min=0.2,
        )

        np.testing.assert_allclose(all_counts, [6.0, 9.0])
        np.testing.assert_allclose(common_counts, [5.0, 9.0])

    def test_regression_mask_from_keys_uses_identity_restriction_match_kind(self):
        metadata = pd.DataFrame(
            {
                "CHR": ["1", "1", "1"],
                "POS": [100, 100, 200],
                "SNP": ["rs1", "rs1", "rs2"],
                "A1": ["A", "A", "A"],
                "A2": ["C", "G", "C"],
            }
        )
        restriction = RestrictionIdentityKeys(
            keys={"rs1:A:C", "rs2:A:C"},
            match_kind="identity",
            dropped=empty_identity_drop_frame(),
            n_input_rows=2,
            n_retained_keys=2,
        )

        mask = kernel_ldscore.regression_mask_from_keys(metadata, restriction, "rsid_allele_aware")

        np.testing.assert_array_equal(mask, np.array([1.0, 0.0, 1.0], dtype=np.float32))

    def test_frequency_metadata_uses_allele_aware_keys_when_available(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            freq_path = Path(tmpdir) / "freq.tsv"
            freq_path.write_text(
                "CHR\tPOS\tSNP\tA1\tA2\tCM\tMAF\n"
                "1\t10\trs1\tA\tC\t0.1\t0.2\n"
                "1\t20\trs2\tA\tG\t0.2\t0.3\n",
                encoding="utf-8",
            )
            metadata = pd.DataFrame(
                {
                    "CHR": ["1", "1"],
                    "POS": [10, 20],
                    "SNP": ["rs1", "rs2"],
                    "A1": ["A", "A"],
                    "A2": ["C", "G"],
                    "CM": [np.nan, np.nan],
                }
            )

            merged = kernel_ldscore.merge_frequency_metadata(
                metadata,
                Namespace(frqfile=str(freq_path)),
                chrom="1",
                identifier_mode="chr_pos_allele_aware",
            )

        self.assertEqual(merged["CM"].tolist(), [0.1, 0.2])
        self.assertEqual(merged["MAF"].tolist(), [0.2, 0.3])

    def test_frequency_metadata_without_alleles_matches_base_keys(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            freq_path = Path(tmpdir) / "freq.tsv"
            freq_path.write_text(
                "CHR\tPOS\tSNP\tCM\tMAF\n"
                "1\t10\trs1\t0.1\t0.2\n"
                "1\t20\trs2\t0.2\t0.3\n",
                encoding="utf-8",
            )
            metadata = pd.DataFrame(
                {
                    "CHR": ["1", "1"],
                    "POS": [10, 20],
                    "SNP": ["rs1", "rs2"],
                    "A1": ["A", "A"],
                    "A2": ["C", "G"],
                    "CM": [np.nan, np.nan],
                }
            )

            merged = kernel_ldscore.merge_frequency_metadata(
                metadata,
                Namespace(frqfile=str(freq_path)),
                chrom="1",
                identifier_mode="chr_pos_allele_aware",
            )

        self.assertEqual(merged["CM"].tolist(), [0.1, 0.2])
        self.assertEqual(merged["MAF"].tolist(), [0.2, 0.3])

    def test_plink_compute_accepts_restriction_identity_keys(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            prefix = self._copy_plink_fixture_with_distinct_fids(Path(tmpdir))
            workflow_bundle = self._build_annotation_bundle(prefix)
            bundle = kernel_ldscore.AnnotationBundle(
                metadata=workflow_bundle.metadata,
                annotations=workflow_bundle.baseline_annotations,
                baseline_columns=workflow_bundle.baseline_columns,
                query_columns=workflow_bundle.query_columns,
            )
            args = Namespace(
                bfile=str(prefix),
                keep=None,
                maf_min=None,
                maf=None,
                ld_wind_snps=10,
                ld_wind_kb=None,
                ld_wind_cm=None,
                yes_really=True,
                snp_batch_size=50,
                common_maf_min=0.05,
                snp_identifier="rsid",
            )
            restriction = RestrictionIdentityKeys(
                keys={"rs_1", "rs_2"},
                match_kind="base",
                dropped=empty_identity_drop_frame(),
                n_input_rows=2,
                n_retained_keys=2,
            )

            result = kernel_ldscore.compute_chrom_from_plink("1", bundle, args, restriction)

            self.assertEqual(result.chrom, "1")
            self.assertGreater(len(result.metadata), 0)
            self.assertEqual(result.w_ld.shape[0], len(result.metadata))

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
            config_snapshot=GlobalConfig(genome_build="hg19", snp_identifier="chr_pos"),
        )
        ref_panel = SimpleNamespace(spec=SimpleNamespace(genome_build="hg19", backend="plink"))

        with mock.patch.object(calc, "compute_chromosome") as patched_compute:
            with self.assertRaisesRegex(ConfigMismatchError, "AnnotationBundle and LDScoreCalculator runtime config"):
                calc.run(
                    annotation_bundle=annotation_bundle,
                    ref_panel=ref_panel,
                    ldscore_config=LDScoreConfig(ld_wind_snps=10),
                    global_config=GlobalConfig(genome_build="hg38", snp_identifier="chr_pos"),
                )

        patched_compute.assert_not_called()

    def test_run_ldscore_from_args_writes_outputs(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            args = Namespace(
                output_dir=str(tmpdir / "ldscore_result"),
                query_annot_sources=None,
                query_annot_chr=None,
                baseline_annot_sources="baseline.annot.gz",
                baseline_annot_chr=None,
                plink_prefix="panel",
                bfile_chr=None,
                r2_table_chr=None,
                snp_identifier="rsid",
                genome_build=None,
                ref_panel_snps_file=None,
                regression_snps_file=None,
                r2_bias_mode=None,
                r2_sample_size=None,
                frqfile_chr=None,
                ld_wind_snps=10,
                ld_wind_kb=None,
                ld_wind_cm=None,
                maf_min=None,
                common_maf_min=0.05,
                snp_batch_size=50,
                per_chr_output=False,
                yes_really=False,
                log_level="INFO",
            )
            annotation_bundle = self.make_annotation_bundle([("1", "rs1", 10)])
            ref_panel = self.make_ref_panel_stub(backend="plink")
            with mock.patch.object(ldscore_workflow.kernel_ldscore, "validate_args"), mock.patch(
                "ldsc.annotation_builder.AnnotationBuilder.run",
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
            self.assertEqual(result.baseline_table["SNP"].tolist(), ["rs1"])
            self.assertIn("manifest", result.output_paths)
            self.assertIn("baseline", result.output_paths)
            self.assertIn("query", result.output_paths)
            self.assertNotIn("log", result.output_paths)
            self.assertTrue(Path(result.output_paths["manifest"]).exists())
            self.assertTrue((tmpdir / "ldscore_result" / "ldscore.log").exists())
            self.assertFalse(list((tmpdir / "ldscore_result").glob("*.l2.ldscore.gz")))
            self.assertFalse(list((tmpdir / "ldscore_result").glob("*.M*")))

    def test_run_ldscore_from_args_synthesizes_base_when_baseline_and_query_are_omitted(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            args = Namespace(
                output_dir=str(tmpdir / "ldscore_result"),
                query_annot_sources=None,
                query_annot_chr=None,
                query_annot_bed_sources=None,
                baseline_annot_sources=None,
                baseline_annot_chr=None,
                plink_prefix="panel",
                bfile_chr=None,
                r2_table_chr=None,
                snp_identifier="rsid",
                genome_build="hg38",
                ref_panel_snps_file=None,
                regression_snps_file=None,
                r2_bias_mode=None,
                r2_sample_size=None,
                frqfile_chr=None,
                ld_wind_snps=10,
                ld_wind_kb=None,
                ld_wind_cm=None,
                maf_min=None,
                common_maf_min=0.05,
                snp_batch_size=50,
                per_chr_output=False,
                yes_really=False,
                log_level="INFO",
            )
            ref_metadata = pd.DataFrame(
                {
                    "CHR": ["1", "1"],
                    "SNP": ["rs1", "rs2"],
                    "CM": [0.1, 0.2],
                    "POS": [10, 20],
                }
            )
            ref_panel = self.make_ref_panel_stub(backend="plink", metadata=ref_metadata)

            def _compute(_self, chrom, annotation_bundle, ref_panel, ldscore_config, global_config, regression_snps=None):
                self.assertEqual(chrom, "1")
                self.assertEqual(annotation_bundle.metadata["SNP"].tolist(), ["rs1", "rs2"])
                self.assertEqual(annotation_bundle.baseline_columns, ["base"])
                self.assertEqual(annotation_bundle.query_columns, [])
                self.assertEqual(annotation_bundle.baseline_annotations["base"].tolist(), [1.0, 1.0])
                self.assertEqual(len(annotation_bundle.query_annotations), 2)
                return ldscore_workflow.ChromLDScoreResult(
                    chrom="1",
                    baseline_table=pd.DataFrame(
                        {
                            "CHR": ["1", "1"],
                            "SNP": ["rs1", "rs2"],
                            "POS": [10, 20],
                            "regression_ld_scores": [3.0, 4.0],
                            "base": [1.0, 2.0],
                        }
                    ),
                    query_table=None,
                    count_records=[
                        {"group": "baseline", "column": "base", "all_reference_snp_count": 2.0, "common_reference_snp_count": 2.0}
                    ],
                    baseline_columns=["base"],
                    query_columns=[],
                    ld_reference_snps=frozenset(),
                    ld_regression_snps=frozenset({"rs1", "rs2"}),
                    snp_count_totals={
                        "all_reference_snp_counts": np.array([2.0]),
                        "common_reference_snp_counts": np.array([2.0]),
                    },
                    config_snapshot=global_config,
                )

            with mock.patch(
                "ldsc.annotation_builder.AnnotationBuilder.run",
                autospec=True,
                side_effect=AssertionError("AnnotationBuilder should not be called for pseudo-base LD scores."),
            ) as patched_builder, mock.patch(
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

            patched_builder.assert_not_called()
            self.assertEqual(result.baseline_columns, ["base"])
            self.assertEqual(result.query_columns, [])
            self.assertIsNone(result.query_table)
            self.assertNotIn("query", result.output_paths)
            baseline_df = pd.read_parquet(result.output_paths["baseline"])
            self.assertEqual(baseline_df.columns.tolist(), ["CHR", "SNP", "POS", "regression_ld_scores", "base"])
            self.assertEqual(result.count_records[0]["column"], "base")

    def test_run_ldscore_from_args_refuses_stale_query_parquet_without_overwrite(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            output_dir = tmpdir / "ldscore_result"
            output_dir.mkdir()
            stale = output_dir / "ldscore.query.parquet"
            stale.write_text("stale\n", encoding="utf-8")
            args = Namespace(
                output_dir=str(output_dir),
                query_annot_sources=None,
                query_annot_bed_sources=None,
                baseline_annot_sources=None,
                plink_prefix="panel",
                snp_identifier="rsid",
                genome_build="hg38",
                ref_panel_snps_file=None,
                regression_snps_file=None,
                r2_bias_mode=None,
                r2_sample_size=None,
                ld_wind_snps=10,
                ld_wind_kb=None,
                ld_wind_cm=None,
                maf_min=None,
                common_maf_min=0.05,
                snp_batch_size=50,
                yes_really=False,
                overwrite=False,
                log_level="INFO",
            )
            ref_panel = self.make_ref_panel_stub(
                backend="plink",
                metadata=pd.DataFrame({"CHR": ["1"], "SNP": ["rs1"], "CM": [0.1], "POS": [10]}),
            )

            with mock.patch(
                "ldsc._kernel.ref_panel.RefPanelLoader.load",
                autospec=True,
                return_value=ref_panel,
            ), mock.patch.object(
                ldscore_workflow.LDScoreCalculator,
                "compute_chromosome",
                autospec=True,
                side_effect=AssertionError("calculation should not run after preflight failure"),
            ):
                with self.assertRaisesRegex(FileExistsError, "overwrite"):
                    ldscore_workflow.run_ldscore_from_args(args)

            self.assertEqual(stale.read_text(encoding="utf-8"), "stale\n")
            self.assertFalse((output_dir / "manifest.json").exists())
            self.assertFalse((output_dir / "ldscore.log").exists())

    def test_run_ldscore_from_args_overwrite_removes_stale_query_parquet(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            output_dir = tmpdir / "ldscore_result"
            output_dir.mkdir()
            stale = output_dir / "ldscore.query.parquet"
            stale.write_text("stale\n", encoding="utf-8")
            args = Namespace(
                output_dir=str(output_dir),
                query_annot_sources=None,
                query_annot_bed_sources=None,
                baseline_annot_sources=None,
                plink_prefix="panel",
                snp_identifier="rsid",
                genome_build="hg38",
                ref_panel_snps_file=None,
                regression_snps_file=None,
                r2_bias_mode=None,
                r2_sample_size=None,
                ld_wind_snps=10,
                ld_wind_kb=None,
                ld_wind_cm=None,
                maf_min=None,
                common_maf_min=0.05,
                snp_batch_size=50,
                yes_really=False,
                overwrite=True,
                log_level="INFO",
            )
            ref_panel = self.make_ref_panel_stub(
                backend="plink",
                metadata=pd.DataFrame({"CHR": ["1"], "SNP": ["rs1"], "CM": [0.1], "POS": [10]}),
            )

            with mock.patch(
                "ldsc._kernel.ref_panel.RefPanelLoader.load",
                autospec=True,
                return_value=ref_panel,
            ), mock.patch.object(
                ldscore_workflow.LDScoreCalculator,
                "compute_chromosome",
                autospec=True,
                return_value=self.make_baseline_only_chrom_result("1", 10, 1.0, 5.0),
            ):
                result = ldscore_workflow.run_ldscore_from_args(args)

            self.assertTrue((output_dir / "manifest.json").exists())
            self.assertTrue((output_dir / "ldscore.log").exists())
            self.assertFalse(stale.exists())
            self.assertNotIn("query", result.output_paths)

    def test_run_ldscore_from_args_rejects_query_annotations_without_baseline(self):
        args = Namespace(
            output_dir="out/example",
            query_annot_sources="query.annot.gz",
            query_annot_chr=None,
            query_annot_bed_sources=None,
            baseline_annot_sources=None,
            baseline_annot_chr=None,
            plink_prefix="panel",
            bfile_chr=None,
            r2_table_chr=None,
            snp_identifier="rsid",
            genome_build=None,
            ref_panel_snps_file=None,
            regression_snps_file=None,
            r2_bias_mode=None,
            r2_sample_size=None,
            frqfile_chr=None,
            ld_wind_snps=10,
            ld_wind_kb=None,
            ld_wind_cm=None,
            maf_min=None,
            common_maf_min=0.05,
            snp_batch_size=50,
            per_chr_output=False,
            yes_really=False,
            log_level="INFO",
        )

        with self.assertRaisesRegex(ValueError, "Query annotations require baseline annotations.*all-ones `base` baseline annotation"):
            ldscore_workflow.run_ldscore_from_args(args)

    def test_run_ldscore_from_args_rejects_query_bed_without_baseline(self):
        args = Namespace(
            output_dir="out/example",
            query_annot_sources=None,
            query_annot_chr=None,
            query_annot_bed_sources="query.bed",
            baseline_annot_sources=None,
            baseline_annot_chr=None,
            plink_prefix="panel",
            bfile_chr=None,
            r2_table_chr=None,
            snp_identifier="rsid",
            genome_build=None,
            ref_panel_snps_file=None,
            regression_snps_file=None,
            r2_bias_mode=None,
            r2_sample_size=None,
            frqfile_chr=None,
            ld_wind_snps=10,
            ld_wind_kb=None,
            ld_wind_cm=None,
            maf_min=None,
            common_maf_min=0.05,
            snp_batch_size=50,
            per_chr_output=False,
            yes_really=False,
            log_level="INFO",
        )

        with self.assertRaisesRegex(ValueError, "Query annotations require baseline annotations.*all-ones `base` baseline annotation"):
            ldscore_workflow.run_ldscore_from_args(args)

    def test_run_ldscore_from_args_loads_regression_snps_and_writes_filtered_ldscore(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            regression_snps_file = tmpdir / "regression_snps.txt"
            regression_snps_file.write_text("SNP\nrs2\n", encoding="utf-8")
            args = Namespace(
                output_dir=str(tmpdir / "ldscore_result"),
                query_annot_sources=None,
                query_annot_chr=None,
                baseline_annot_sources="baseline.annot.gz",
                baseline_annot_chr=None,
                plink_prefix="panel",
                bfile_chr=None,
                r2_table_chr=None,
                snp_identifier="rsid",
                genome_build=None,
                ref_panel_snps_file=None,
                regression_snps_file=str(regression_snps_file),
                r2_bias_mode=None,
                r2_sample_size=None,
                frqfile_chr=None,
                ld_wind_snps=10,
                ld_wind_kb=None,
                ld_wind_cm=None,
                maf_min=None,
                common_maf_min=0.05,
                snp_batch_size=50,
                per_chr_output=False,
                yes_really=False,
                log_level="INFO",
            )

            def _compute(_self, chrom, annotation_bundle, ref_panel, ldscore_config, global_config, regression_snps=None):
                self.assertEqual(chrom, "1")
                self.assertEqual(regression_snps.keys, {"rs2"})
                self.assertEqual(regression_snps.match_kind, "base")
                return ldscore_workflow.ChromLDScoreResult(
                    chrom="1",
                    baseline_table=pd.DataFrame(
                        {
                            "CHR": ["1"],
                            "SNP": ["rs2"],
                            "POS": [20],
                            "regression_ld_scores": [6.0],
                            "base": [3.0],
                        }
                    ),
                    query_table=pd.DataFrame(
                        {
                            "CHR": ["1"],
                            "SNP": ["rs2"],
                            "POS": [20],
                            "query": [4.0],
                        }
                    ),
                    count_records=[
                        {"group": "baseline", "column": "base", "all_reference_snp_count": 7.0, "common_reference_snp_count": 6.0},
                        {"group": "query", "column": "query", "all_reference_snp_count": 8.0, "common_reference_snp_count": 7.0},
                    ],
                    baseline_columns=["base"],
                    query_columns=["query"],
                    ld_reference_snps=frozenset(),
                    ld_regression_snps=frozenset({"rs2"}),
                    snp_count_totals={
                        "all_reference_snp_counts": np.array([7.0, 8.0]),
                        "common_reference_snp_counts": np.array([6.0, 7.0]),
                    },
                    config_snapshot=global_config,
                )

            annotation_bundle = self.make_annotation_bundle([("1", "rs1", 10), ("1", "rs2", 20)])
            ref_panel = self.make_ref_panel_stub(backend="plink")
            with mock.patch.object(ldscore_workflow.kernel_ldscore, "validate_args"), mock.patch(
                "ldsc.annotation_builder.AnnotationBuilder.run",
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
        from ldsc._kernel.ref_panel import RefPanelConfig

        common = GlobalConfig(snp_identifier="chr_pos", genome_build="hg19")
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
            spec = RefPanelConfig(
                backend="parquet_r2",
                r2_dir=str(tmpdir / "panel" / "hg19"),
            )
            ref_panel = SimpleNamespace(
                spec=spec,
                resolve_r2_paths=lambda chrom, required=False: [str(r2_path)],
                resolve_metadata_paths=lambda chrom: [str(meta_path)],
            )

            args = ldscore_workflow._namespace_from_configs(
                chrom="1",
                ref_panel=ref_panel,
                ldscore_config=ldscore_workflow.LDScoreConfig(
                    ld_wind_cm=1.0,
                    snp_batch_size=64,
                    common_maf_min=0.05,
                ),
                global_config=common,
            )

            self.assertIsInstance(args.r2_table, str)
            self.assertIsInstance(args.frqfile, str)
            self.assertEqual(args.r2_table, str(r2_path))
            self.assertEqual(args.frqfile, str(meta_path))
            self.assertIsNone(args.keep)
            self.assertEqual(args.genome_build, "hg19")
            self.assertEqual(args.snp_batch_size, 64)
            self.assertFalse(hasattr(args, "chunk_size"))

    def test_run_ldscore_from_args_passes_path_tokens_to_builder_and_ref_panel_loader(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            args = Namespace(
                output_dir=str(tmpdir / "ldscore_result"),
                query_annot_sources=None,
                query_annot_chr=None,
                baseline_annot_sources=str(tmpdir / "baseline.@.annot.gz"),
                baseline_annot_chr=None,
                plink_prefix=str(tmpdir / "panel.@"),
                bfile_chr=None,
                r2_table_chr=None,
                snp_identifier="rsid",
                genome_build=None,
                ref_panel_snps_file=None,
                regression_snps_file=None,
                r2_bias_mode=None,
                r2_sample_size=None,
                frqfile_chr=None,
                ld_wind_snps=10,
                ld_wind_kb=None,
                ld_wind_cm=None,
                maf_min=None,
                common_maf_min=0.05,
                snp_batch_size=50,
                per_chr_output=False,
                yes_really=False,
                log_level="INFO",
            )
            annotation_bundle = self.make_annotation_bundle([("1", "rs1", 10)])
            ref_panel = self.make_ref_panel_stub(backend="plink", plink_prefix=str(tmpdir / "panel.@"))
            with mock.patch.object(ldscore_workflow.kernel_ldscore, "validate_args"), mock.patch(
                "ldsc.annotation_builder.AnnotationBuilder.run",
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
            self.assertEqual(source_spec.baseline_annot_sources, (str(tmpdir / "baseline.@.annot.gz"),))
            self.assertEqual(source_spec.query_annot_sources, ())
            self.assertEqual(source_spec.query_annot_bed_sources, ())
            ref_spec = patched_loader.call_args.args[1]
            self.assertEqual(ref_spec.backend, "plink")
            self.assertEqual(str(ref_spec.plink_prefix), str(tmpdir / "panel.@"))

    def test_run_ldscore_from_args_passes_r2_dir_to_ref_panel_loader(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            args = Namespace(
                output_dir=str(tmpdir / "ldscore_result"),
                query_annot_sources=None,
                query_annot_chr=None,
                baseline_annot_sources=str(tmpdir / "baseline.@.annot.gz"),
                baseline_annot_chr=None,
                plink_prefix=None,
                bfile_chr=None,
                r2_table_chr=None,
                r2_dir=str(tmpdir / "panel" / "hg38"),
                snp_identifier="rsid",
                genome_build=None,
                ref_panel_snps_file=None,
                regression_snps_file=None,
                r2_bias_mode="unbiased",
                r2_sample_size=None,
                frqfile_chr=None,
                ld_wind_snps=10,
                ld_wind_kb=None,
                ld_wind_cm=None,
                maf_min=None,
                common_maf_min=0.05,
                snp_batch_size=50,
                per_chr_output=False,
                yes_really=False,
                log_level="INFO",
            )
            annotation_bundle = self.make_annotation_bundle([("1", "rs1", 10)])
            ref_panel = self.make_ref_panel_stub(
                backend="parquet_r2",
                genome_build=None,
            )
            with mock.patch.object(ldscore_workflow.kernel_ldscore, "validate_args"), mock.patch(
                "ldsc.annotation_builder.AnnotationBuilder.run",
                autospec=True,
                return_value=annotation_bundle,
            ), mock.patch(
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

            ref_spec = patched_loader.call_args.args[1]
            self.assertEqual(ref_spec.backend, "parquet_r2")
            self.assertEqual(ref_spec.r2_dir, str(tmpdir / "panel" / "hg38"))

    def test_normalize_run_args_chr_pos_auto_infers_matching_annotation_and_r2_dir_builds(self):
        args = Namespace(
            output_dir="out",
            query_annot_sources=None,
            baseline_annot_sources="baseline.@.annot.gz",
            plink_prefix=None,
            r2_dir="panel/hg19",
            snp_identifier="chr_pos",
            genome_build="auto",
            keep_indivs_file=None,
            ref_panel_snps_file=None,
            regression_snps_file=None,
            log_level="INFO",
        )
        sample_calls = []
        resolve_calls = []

        def fake_sample(tokens, *, context, **_kwargs):
            sample_calls.append((tuple(tokens), context))
            return pd.DataFrame({"CHR": ["1"], "POS": [100]}), f"{context}.sample"

        def fake_resolve(genome_build, snp_identifier, sample_frame, *, context, **_kwargs):
            resolve_calls.append((genome_build, snp_identifier, context, sample_frame["POS"].tolist()))
            return "hg19"

        with mock.patch.object(ldscore_workflow, "sample_frame_from_chr_pattern", side_effect=fake_sample), mock.patch.object(
            ldscore_workflow,
            "resolve_genome_build",
            side_effect=fake_resolve,
        ):
            normalized_args, global_config = ldscore_workflow._normalize_run_args(args)

        self.assertEqual(global_config.genome_build, "hg19")
        self.assertEqual(normalized_args.genome_build, "hg19")
        self.assertEqual(
            sample_calls,
            [
                (("baseline.@.annot.gz",), "LD-score annotation inputs"),
            ],
        )
        self.assertEqual([call[2] for call in resolve_calls], ["LD-score annotation inputs"])

    def test_normalize_run_args_chr_pos_auto_does_not_infer_from_ref_panel_child_directory(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            root = tmpdir / "r2_ref_panel_1kg30x_1cM_hm3"
            (root / "hg38").mkdir(parents=True)
            args = Namespace(
                output_dir="out",
                query_annot_sources=None,
                baseline_annot_sources=None,
                plink_prefix=None,
                r2_dir=str(root),
                snp_identifier="chr_pos",
                genome_build="auto",
                keep_indivs_file=None,
                ref_panel_snps_file=None,
                regression_snps_file=None,
                log_level="INFO",
            )

            with self.assertRaisesRegex(ValueError, "Cannot infer --genome-build"):
                ldscore_workflow._normalize_run_args(args)

    def test_normalize_run_args_rsid_concrete_genome_build_warns_and_nulls(self):
        args = Namespace(
            output_dir="out",
            query_annot_sources=None,
            baseline_annot_sources="baseline.annot.gz",
            plink_prefix="panel",
            snp_identifier="rsid",
            genome_build="hg38",
            keep_indivs_file=None,
            ref_panel_snps_file=None,
            regression_snps_file=None,
            log_level="INFO",
        )

        with self.assertWarnsRegex(UserWarning, "ignored"):
            normalized_args, global_config = ldscore_workflow._normalize_run_args(args)

        self.assertIsNone(normalized_args.genome_build)
        self.assertIsNone(global_config.genome_build)

    @unittest.skipUnless(_HAS_PYARROW, "pyarrow is required for parquet schema coverage")
    def test_normalize_run_args_chr_pos_auto_rejects_annotation_r2_build_mismatch(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            r2_dir = Path(tmpdir) / "panel" / "hg38"
            _write_minimal_r2_parquet(
                r2_dir / "chr1_r2.parquet",
                {b"ldsc:sorted_by_build": b"hg38"},
            )
            args = Namespace(
                output_dir="out",
                query_annot_sources=None,
                baseline_annot_sources="baseline.@.annot.gz",
                plink_prefix=None,
                r2_dir=str(r2_dir),
                snp_identifier="chr_pos",
                genome_build="auto",
                keep_indivs_file=None,
                ref_panel_snps_file=None,
                regression_snps_file=None,
                log_level="INFO",
            )

            with mock.patch.object(
                ldscore_workflow,
                "sample_frame_from_chr_pattern",
                return_value=(pd.DataFrame({"CHR": ["1"], "POS": [100]}), "sample"),
            ), mock.patch.object(
                ldscore_workflow,
                "resolve_genome_build",
                return_value="hg19",
            ):
                with self.assertRaisesRegex(ValueError, "genome build.*disagree"):
                    ldscore_workflow._normalize_run_args(args)

    def test_run_ldscore_from_args_rejects_keep_in_parquet_mode(self):
        args = Namespace(
            output_dir="out/example",
            query_annot_sources=None,
            query_annot_chr=None,
            baseline_annot_sources="baseline.annot.gz",
            baseline_annot_chr=None,
            plink_prefix=None,
            bfile_chr=None,
            r2_dir="panel/hg38",
            r2_table_chr=None,
            snp_identifier="rsid",
            genome_build=None,
            ref_panel_snps_file=None,
            regression_snps_file=None,
            r2_bias_mode="unbiased",
            r2_sample_size=None,
            frqfile_chr=None,
            keep_indivs_file="samples.keep",
            ld_wind_snps=10,
            ld_wind_kb=None,
            ld_wind_cm=None,
            maf_min=None,
            common_maf_min=0.05,
            snp_batch_size=50,
            per_chr_output=False,
            yes_really=False,
            log_level="INFO",
        )

        with self.assertRaisesRegex(ValueError, "--keep-indivs-file.*PLINK"):
            ldscore_workflow.run_ldscore_from_args(args)

    def test_run_ldscore_from_args_warns_and_skips_empty_intersection_chromosome(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            args = Namespace(
                output_dir=str(tmpdir / "ldscore_result"),
                query_annot_sources=None,
                query_annot_chr=None,
                baseline_annot_sources=str(tmpdir / "baseline.@.annot.gz"),
                baseline_annot_chr=None,
                plink_prefix=None,
                bfile_chr=None,
                r2_dir=str(tmpdir / "r2_panel" / "hg19"),
                r2_table_chr=None,
                snp_identifier="rsid",
                genome_build="hg19",
                ref_panel_snps_file=None,
                regression_snps_file=None,
                r2_bias_mode="unbiased",
                r2_sample_size=None,
                frqfile_chr=None,
                ld_wind_snps=10,
                ld_wind_kb=None,
                ld_wind_cm=None,
                maf_min=None,
                common_maf_min=0.05,
                snp_batch_size=50,
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
                            "POS": [220],
                            "regression_ld_scores": [3.0],
                            "base": [2.0],
                        }
                    ),
                    query_table=None,
                    count_records=[
                        {"group": "baseline", "column": "base", "all_reference_snp_count": 5.0, "common_reference_snp_count": 4.0}
                    ],
                    baseline_columns=["base"],
                    query_columns=[],
                    ld_reference_snps=frozenset(),
                    ld_regression_snps=frozenset({"rs22"}),
                    snp_count_totals={
                        "all_reference_snp_counts": np.array([5.0]),
                        "common_reference_snp_counts": np.array([4.0]),
                    },
                    config_snapshot=global_config,
                )

            annotation_bundle = self.make_annotation_bundle([("1", "rs1", 10), ("22", "rs22", 220)], genome_build="hg19")
            ref_panel = self.make_ref_panel_stub(backend="parquet_r2", genome_build="hg19")
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                with mock.patch.object(ldscore_workflow.kernel_ldscore, "validate_args"), mock.patch(
                    "ldsc.annotation_builder.AnnotationBuilder.run",
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

    def test_compute_chromosome_aligns_allele_free_annotations_in_base_family_mode(self):
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
                    "A1": ["A", "A"],
                    "A2": ["C", "G"],
                }
            )
        )

        def _compute_side_effect(chrom, bundle, args, regression_snps):
            self.assertEqual(chrom, "1")
            self.assertEqual(bundle.metadata["SNP"].tolist(), ["rs1", "rs3"])
            self.assertEqual(bundle.metadata["A1"].tolist(), ["A", "A"])
            self.assertEqual(bundle.metadata["A2"].tolist(), ["C", "G"])
            return ldscore_workflow._LegacyChromResult(
                chrom="1",
                metadata=pd.DataFrame(
                    {
                        "CHR": ["1", "1"],
                        "SNP": ["rs1", "rs3"],
                        "POS": [10, 30],
                        "A1": ["A", "A"],
                        "A2": ["C", "G"],
                        "CM": [0.1, 0.3],
                        "MAF": [0.2, 0.2],
                    }
                ),
                ld_scores=np.array([[1.0], [2.0]], dtype=np.float32),
                w_ld=np.array([[3.0], [4.0]], dtype=np.float32),
                M=np.array([2.0]),
                M_5_50=None,
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
                global_config=GlobalConfig(snp_identifier="chr_pos_allele_aware", genome_build="hg38"),
            )

        ref_panel.load_metadata.assert_called_once_with("1")
        self.assertEqual(result.baseline_table["SNP"].tolist(), ["rs1", "rs3"])

    def test_pseudo_base_annotation_preserves_ref_panel_alleles(self):
        ref_panel = self.make_ref_panel_stub(
            backend="parquet_r2",
            metadata=pd.DataFrame(
                {
                    "CHR": ["1", "1"],
                    "SNP": ["rs1", "rs2"],
                    "CM": [0.1, 0.2],
                    "POS": [10, 20],
                    "A1": ["A", "A"],
                    "A2": ["C", "G"],
                }
            ),
        )

        bundle = ldscore_workflow._pseudo_base_annotation_bundle_from_ref_panel(
            ref_panel,
            GlobalConfig(snp_identifier="chr_pos_allele_aware", genome_build="hg38"),
        )

        self.assertEqual(bundle.metadata["A1"].tolist(), ["A", "A"])
        self.assertEqual(bundle.metadata["A2"].tolist(), ["C", "G"])

    @unittest.skipUnless(_HAS_BITARRAY, "bitarray is not installed")
    def test_plink_compute_chromosome_enriches_allele_free_annotations_in_default_mode(self):
        prefix = Path(__file__).resolve().parent / "fixtures" / "minimal_external_resources" / "plink" / "hm3_chr22_subset"
        bim = pd.read_csv(
            prefix.with_suffix(".bim"),
            sep=r"\s+",
            header=None,
            names=["CHR", "SNP", "CM", "POS", "A1", "A2"],
        )
        metadata = bim.loc[:, ["CHR", "SNP", "CM", "POS"]].copy()
        metadata["CHR"] = metadata["CHR"].astype(str)
        baseline = pd.DataFrame({"base": np.ones(len(metadata), dtype=np.float32)})
        bundle = AnnotationBundle(
            metadata=metadata,
            baseline_annotations=baseline,
            query_annotations=pd.DataFrame(index=metadata.index),
            baseline_columns=["base"],
            query_columns=[],
            chromosomes=["22"],
            source_summary={},
            config_snapshot=None,
        )
        common = GlobalConfig(snp_identifier="chr_pos_allele_aware", genome_build="hg38")
        panel = PlinkRefPanel(common, RefPanelConfig(backend="plink", plink_prefix=str(prefix)))

        result = ldscore_workflow.LDScoreCalculator().compute_chromosome(
            chrom="22",
            annotation_bundle=bundle,
            ref_panel=panel,
            ldscore_config=LDScoreConfig(ld_wind_snps=10, whole_chromosome_ok=True),
            global_config=common,
        )

        self.assertGreater(len(result.baseline_table), 0)
        self.assertTrue(all(key.startswith("22:") for key in result.ld_regression_snps))

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
            panel = PlinkRefPanel(common, RefPanelConfig(backend="plink", plink_prefix=prefix, keep_indivs_file=keep_path))
            result = ldscore_workflow.LDScoreCalculator().run(
                annotation_bundle=bundle,
                ref_panel=panel,
                ldscore_config=LDScoreConfig(
                    ld_wind_snps=10,
                    whole_chromosome_ok=True,
                ),
                global_config=common,
            )

            self.assertEqual(result.baseline_table["SNP"].tolist(), expected_metadata["SNP"].tolist())
            np.testing.assert_allclose(result.baseline_table["base"].to_numpy(), expected_ld)
            np.testing.assert_allclose(result.baseline_table["regression_ld_scores"].to_numpy(), expected_ld)

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
                RefPanelConfig(
                    backend="plink",
                    plink_prefix=prefix,
                    ref_panel_snps_file=restrict_path,
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

            self.assertEqual(result.baseline_table["SNP"].tolist(), ["rs_6"])
            self.assertEqual(result.ld_regression_snps, frozenset({"rs_6"}))
            self.assertEqual(result.count_records[0]["all_reference_snp_count"], 1.0)

    @unittest.skipUnless(_HAS_BITARRAY, "bitarray is not installed")
    def test_ldscore_calculator_run_filters_individuals_before_maf_in_plink_mode(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            prefix = self._copy_plink_fixture_with_distinct_fids(tmpdir)
            keep_path = tmpdir / "samples.keep"
            self._write_keep_file(keep_path, ["per0", "per1"])

            expected_metadata, expected_ld = self._expected_plink_result(prefix, keep_path, maf_min=0.01)
            bundle = self._build_annotation_bundle(prefix)
            common = GlobalConfig(snp_identifier="rsid")
            panel = PlinkRefPanel(
                common,
                RefPanelConfig(
                    backend="plink",
                    plink_prefix=prefix,
                    maf_min=0.01,
                    keep_indivs_file=keep_path,
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

    def test_canonicalize_r2_pairs_normalizes_mixed_chromosome_labels(self):
        df = pd.DataFrame(
            {
                "chr": ["chr1.0", "1", "24"],
                "rsID_1": ["rs2", "rs4", "rsY2"],
                "rsID_2": ["rs1", "rs3", "rsY1"],
                "hg38_pos_1": [120, 140, 80],
                "hg38_pos_2": [100, 130, 70],
                "hg19_pos_1": [20, 40, 8],
                "hg19_pos_2": [10, 30, 7],
                "hg38_Uniq_ID_1": ["1:120", "1:140", "Y:80"],
                "hg38_Uniq_ID_2": ["1:100", "1:130", "Y:70"],
                "hg19_Uniq_ID_1": ["1:20", "1:40", "Y:8"],
                "hg19_Uniq_ID_2": ["1:10", "1:30", "Y:7"],
                "R2": [0.4, 0.2, 0.1],
                "Dprime": [0.5, 0.3, 0.2],
                "+/-corr": ["+", "+", "-"],
            }
        )

        out = kernel_ldscore.canonicalize_r2_pairs(df, "hg19")

        self.assertEqual(out["chr"].tolist(), ["1", "1", "Y"])
        self.assertEqual(out["pos_1"].tolist(), [10, 30, 7])
        self.assertEqual(out["pos_2"].tolist(), [20, 40, 8])

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
                identifier_mode="rsid",
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
    def test_sorted_r2_block_reader_rejects_auto_genome_build_at_kernel_boundary(self):
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
            with self.assertRaisesRegex(AssertionError, "concrete"):
                reader = kernel_ldscore.SortedR2BlockReader(
                    paths=[str(path)],
                    chrom="1",
                    metadata=metadata,
                    identifier_mode="chr_pos",
                    r2_bias_mode="unbiased",
                    r2_sample_size=None,
                    genome_build="auto",
                )
