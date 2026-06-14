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
import pytest

SRC = Path(__file__).resolve().parents[1] / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from ldsc.config import ConfigMismatchError, GlobalConfig, set_global_config
from ldsc.errors import LDSCConfigError, LDSCInputError, LDSCInternalError, LDSCUsageError

try:
    from ldsc import (
        AnnotationBuildConfig,
        AnnotationBuilder,
        AnnotationBundle,
        LDScoreConfig,
        PlinkRefPanel,
        RefPanelConfig,
    )
    from ldsc import ldscore_calculator as ldscore_workflow
    from ldsc.ldscore_calculator import (
        build_parser,
        _ref_panel_from_args,
        _normalize_run_args,
        _resolve_exclude_regions_build,
        run_ldscore,
    )
    from ldsc._kernel import formats as kernel_formats
    from ldsc._kernel import ldscore as kernel_ldscore
    from ldsc._kernel.snp_identity import RestrictionIdentityKeys, empty_identity_drop_frame
except ImportError:
    AnnotationBuildConfig = None
    AnnotationBuilder = None
    AnnotationBundle = None
    LDScoreConfig = None
    PlinkRefPanel = None
    RefPanelConfig = None
    kernel_formats = None
    ldscore_workflow = None
    build_parser = None
    _ref_panel_from_args = None
    _normalize_run_args = None
    run_ldscore = None
    kernel_ldscore = None
    RestrictionIdentityKeys = None
    empty_identity_drop_frame = None


PLINK_FIXTURES = Path(__file__).resolve().parent / "fixtures" / "plink"


def dict_chunks(rows):
    """Wrap legacy {i,j,R2,sign} dict rows as the chunk iterable write_r2_parquet expects."""
    rows = list(rows)
    if not rows:
        return iter([])
    i = np.array([r["i"] for r in rows], dtype=np.int64)
    j = np.array([r["j"] for r in rows], dtype=np.int64)
    r2 = np.array([r["R2"] for r in rows], dtype=np.float32)
    sign = np.array([1 if r["sign"] == "+" else -1 for r in rows], dtype=np.int8)
    return [(i, j, r2, sign)]


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


def _write_legacy_r2_parquet(path: Path) -> None:
    """Write a zero-row parquet table without LDSC schema metadata."""
    import pyarrow as pa
    import pyarrow.parquet as pq

    path.parent.mkdir(parents=True, exist_ok=True)
    table = pa.table(
        {
            "POS_1": pa.array([], type=pa.int64()),
            "POS_2": pa.array([], type=pa.int64()),
            "R2": pa.array([], type=pa.float32()),
        }
    )
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
        import pyarrow as pa
        import pyarrow.parquet as pq
        from ldsc._kernel.ref_panel import ParquetR2RefPanel
        from ldsc._kernel.snp_identity import sidecar_identity_sha256

        with tempfile.TemporaryDirectory() as tmpdir:
            build_dir = Path(tmpdir) / "hg19"
            build_dir.mkdir(parents=True)
            path = build_dir / "chr1_r2.parquet"
            panel = pd.DataFrame({"CHR": ["1", "1"], "POS": [100, 120], "SNP": ["rs1", "rs2"],
                                  "A1": ["A", "A"], "A2": ["C", "C"], "CM": [0.0, 0.0], "MAF": [0.3, 0.3]})
            with gzip.open(build_dir / "chr1_meta.tsv.gz", "wt") as handle:
                handle.write("# ldsc:schema_version=2\n")
                panel.to_csv(handle, sep="\t", index=False)
            # A hand-built raw-bias index panel (the package writer only emits
            # unbiased; raw is reserved for future external raw-R2 index panels).
            schema = pa.schema(
                [("IDX_1", pa.int32()), ("IDX_2", pa.int32()), ("R2", pa.float32()), ("SIGN", pa.bool_())]
            ).with_metadata({
                b"ldsc:schema_version": b"1",
                b"ldsc:artifact_type": b"ref_panel_r2",
                b"ldsc:snp_identifier": b"rsid",
                b"ldsc:genome_build": b"hg19",
                b"ldsc:sorted_by_build": b"hg19",
                b"ldsc:n_samples": b"200",
                b"ldsc:r2_bias": b"raw",
                b"ldsc:n_snps": b"2",
                b"ldsc:sidecar_identity_sha256": sidecar_identity_sha256(panel).encode("utf-8"),
            })
            tbl = pa.table(
                {"IDX_1": pa.array([0], pa.int32()), "IDX_2": pa.array([1], pa.int32()),
                 "R2": pa.array([0.5], pa.float32()), "SIGN": pa.array([True], pa.bool_())},
                schema=schema,
            )
            pq.write_table(tbl, str(path))

            spec = RefPanelConfig(
                backend="parquet_r2",
                r2_dir=str(build_dir),
            )
            panel_loader = ParquetR2RefPanel(GlobalConfig(snp_identifier="rsid"), spec)
            reader_meta = pd.DataFrame(
                {"CHR": ["1", "1"], "POS": [100, 120], "SNP": ["rs1", "rs2"], "CM": [0.0, 0.0], "MAF": [0.3, 0.3]}
            )
            reader = panel_loader.build_reader("1", metadata=reader_meta)

        self.assertEqual(reader.r2_bias_mode, "raw")
        self.assertAlmostEqual(reader.r2_sample_size, 200.0)

    @unittest.skipUnless(_HAS_PYARROW, "pyarrow is required for parquet schema coverage")
    def test_external_parquet_metadata_drops_duplicate_identity_clusters(self):
        from ldsc._kernel.ref_panel import ParquetR2RefPanel

        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            _write_legacy_r2_parquet(root / "chr1_r2.parquet")
            with gzip.open(root / "chr1_meta.tsv.gz", "wt", encoding="utf-8") as handle:
                handle.write(
                    "CHR\tPOS\tSNP\tCM\tMAF\n"
                    "1\t10\trs1\t0.1\t0.2\n"
                    "1\t10\trs2\t0.2\t0.3\n"
                    "1\t20\trs3\t0.3\t0.4\n"
                )
            panel = ParquetR2RefPanel(
                GlobalConfig(snp_identifier="chr_pos", genome_build="hg19"),
                RefPanelConfig(backend="parquet_r2", r2_dir=str(root)),
            )

            with self.assertLogs("LDSC.ref_panel", level="WARNING"):
                metadata = panel.load_metadata("1")

        self.assertEqual(metadata["SNP"].tolist(), ["rs3"])
        self.assertEqual(metadata["POS"].tolist(), [20])

    @unittest.skipUnless(_HAS_PYARROW, "pyarrow is required for parquet schema coverage")
    def test_external_parquet_metadata_raises_when_identity_cleanup_drops_all_rows(self):
        from ldsc._kernel.ref_panel import ParquetR2RefPanel

        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            _write_legacy_r2_parquet(root / "chr1_r2.parquet")
            with gzip.open(root / "chr1_meta.tsv.gz", "wt", encoding="utf-8") as handle:
                handle.write(
                    "CHR\tPOS\tSNP\tCM\tMAF\n"
                    "1\t10\trs1\t0.1\t0.2\n"
                    "1\t10\trs2\t0.2\t0.3\n"
                )
            panel = ParquetR2RefPanel(
                GlobalConfig(snp_identifier="chr_pos", genome_build="hg19"),
                RefPanelConfig(backend="parquet_r2", r2_dir=str(root)),
            )

            with self.assertRaisesRegex(LDSCInputError, "retained no parquet metadata rows"):
                panel.load_metadata("1")

    @unittest.skipUnless(_HAS_PYARROW, "pyarrow is required for parquet schema coverage")
    def test_package_parquet_metadata_duplicates_remain_invariant_failures(self):
        from ldsc._kernel.ref_panel import ParquetR2RefPanel

        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            _write_minimal_r2_parquet(
                root / "chr1_r2.parquet",
                {
                    b"ldsc:schema_version": b"1",
                    b"ldsc:artifact_type": b"ref_panel_r2",
                    b"ldsc:snp_identifier": b"chr_pos",
                    b"ldsc:genome_build": b"hg19",
                },
            )
            with gzip.open(root / "chr1_meta.tsv.gz", "wt", encoding="utf-8") as handle:
                handle.write(
                    "# ldsc:schema_version=1\n"
                    "# ldsc:artifact_type=ref_panel_metadata\n"
                    "# ldsc:snp_identifier=chr_pos\n"
                    "# ldsc:genome_build=hg19\n"
                    "CHR\tPOS\tSNP\tCM\tMAF\n"
                    "1\t10\trs1\t0.1\t0.2\n"
                    "1\t10\trs2\t0.2\t0.3\n"
                    "1\t20\trs3\t0.3\t0.4\n"
                )
            panel = ParquetR2RefPanel(
                GlobalConfig(snp_identifier="chr_pos", genome_build="hg19"),
                RefPanelConfig(backend="parquet_r2", r2_dir=str(root)),
            )

            with self.assertRaisesRegex(LDSCInputError, "non-unique SNP identifiers"):
                panel.load_metadata("1")


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
            ]
        )

        self.assertEqual(args.r2_dir, "panel/hg38")

    def test_build_parser_rejects_removed_r2_bias_flags(self):
        parser = ldscore_workflow.build_parser()
        for flag in (["--r2-bias-mode", "raw"], ["--r2-sample-size", "100"]):
            with self.assertRaises(SystemExit):
                parser.parse_args(["--output-dir", "out", "--r2-dir", "panel/hg38", *flag])

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
                with self.assertRaisesRegex(LDSCUsageError, "removed IO argument"):
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

        with self.assertRaisesRegex(LDSCUsageError, "chr_pos-family SNP identifiers.*genome build"):
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

            with self.assertRaisesRegex(LDSCInputError, "could not infer the genome build"):
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

            with self.assertRaisesRegex(LDSCInputError, "conflicting R2 parquet genome-build metadata"):
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

        with self.assertRaisesRegex(LDSCUsageError, "two regression SNP restrictions"):
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

    def test_wrap_legacy_chrom_result_carries_overlap(self):
        from ldsc._kernel.overlap import OverlapContribution

        contribution = OverlapContribution(np.array([[2.0]]), None, np.array([]), None, 2, None)
        legacy_result = ldscore_workflow._LegacyChromResult(
            chrom="1",
            metadata=pd.DataFrame(
                {"CHR": ["1", "1"], "SNP": ["rs1", "rs2"], "POS": [10, 20], "CM": [0.1, 0.2], "MAF": [0.2, 0.3]}
            ),
            ld_scores=np.array([[1.0], [2.0]], dtype=np.float32),
            w_ld=np.array([[3.0], [4.0]], dtype=np.float32),
            M=np.array([2.0]),
            M_5_50=None,
            ldscore_columns=["base"],
            baseline_columns=["base"],
            query_columns=[],
            overlap=contribution,
        )
        result = ldscore_workflow.LDScoreCalculator()._wrap_legacy_chrom_result(
            legacy_result, global_config=GlobalConfig(snp_identifier="rsid")
        )
        self.assertIs(result.overlap, contribution)

    def test_aggregate_chromosome_results_sums_overlap_blocks(self):
        from ldsc._kernel.overlap import OverlapContribution

        cfg = GlobalConfig(snp_identifier="rsid")

        def make_chrom(chrom, snp, contribution):
            return ldscore_workflow.ChromLDScoreResult(
                chrom=chrom,
                baseline_table=pd.DataFrame(
                    {"CHR": [chrom], "SNP": [snp], "POS": [10], "regression_ld_scores": [3.0], "base": [1.0]}
                ),
                query_table=pd.DataFrame({"CHR": [chrom], "SNP": [snp], "POS": [10], "query": [2.0]}),
                count_records=[
                    {"group": "baseline", "column": "base", "all_reference_snp_count": 10.0},
                    {"group": "query", "column": "query", "all_reference_snp_count": 11.0},
                ],
                baseline_columns=["base"],
                query_columns=["query"],
                ld_reference_snps=frozenset(),
                ld_regression_snps=frozenset({snp}),
                snp_count_totals={"all_reference_snp_counts": np.array([10.0, 11.0])},
                overlap=contribution,
                config_snapshot=cfg,
            )

        c1 = OverlapContribution(np.array([[10.0, 4.0]]), None, np.array([6.0]), None, 10, None)
        c2 = OverlapContribution(np.array([[5.0, 3.0]]), None, np.array([4.0]), None, 7, None)
        result = ldscore_workflow.LDScoreCalculator()._aggregate_chromosome_results(
            [make_chrom("1", "rs1", c1), make_chrom("2", "rs2", c2)], global_config=cfg
        )
        self.assertIsNotNone(result.overlap)
        np.testing.assert_allclose(result.overlap.baseline_block_all.to_numpy(), [[15.0, 7.0]])
        np.testing.assert_allclose(result.overlap.query_diagonal_all.to_numpy(), [10.0])
        self.assertEqual(result.overlap.total_all_reference_snps, 17)
        self.assertEqual(list(result.overlap.baseline_block_all.columns), ["base", "query"])

    def test_aggregate_chromosome_results_skips_overlap_for_base_only(self):
        from ldsc._kernel.overlap import OverlapContribution

        cfg = GlobalConfig(snp_identifier="rsid")

        def make_chrom(chrom, snp, contribution):
            return ldscore_workflow.ChromLDScoreResult(
                chrom=chrom,
                baseline_table=pd.DataFrame(
                    {"CHR": [chrom], "SNP": [snp], "POS": [10], "regression_ld_scores": [3.0], "base": [1.0]}
                ),
                query_table=None,
                count_records=[
                    {"group": "baseline", "column": "base", "all_reference_snp_count": 10.0},
                ],
                baseline_columns=["base"],
                query_columns=[],
                ld_reference_snps=frozenset(),
                ld_regression_snps=frozenset({snp}),
                snp_count_totals={"all_reference_snp_counts": np.array([10.0])},
                overlap=contribution,
                config_snapshot=cfg,
            )

        # A single all-ones base column makes the overlap a degenerate 1x1 (the SNP
        # count, already in metadata) that no downstream regression consumes.
        c1 = OverlapContribution(np.array([[10.0]]), None, np.array([]), None, 10, None)
        c2 = OverlapContribution(np.array([[7.0]]), None, np.array([]), None, 7, None)
        result = ldscore_workflow.LDScoreCalculator()._aggregate_chromosome_results(
            [make_chrom("1", "rs1", c1), make_chrom("2", "rs2", c2)], global_config=cfg
        )
        self.assertIsNone(result.overlap)

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

        self.assertEqual(result.baseline_table.columns.tolist(), ["CHR", "SNP", "POS", "A1", "A2", "regression_ld_scores", "base"])
        self.assertEqual(result.ld_regression_snps, frozenset({"1:10:A:C", "1:20:A:G"}))
        with tempfile.TemporaryDirectory() as tmpdir:
            output_paths = ldscore_workflow.LDScoreDirectoryWriter().write(
                result,
                ldscore_workflow.LDScoreOutputConfig(output_dir=tmpdir),
            )
            baseline = pd.read_parquet(output_paths["baseline"])
        self.assertEqual(baseline["A1"].tolist(), ["A", "A"])
        self.assertEqual(baseline["A2"].tolist(), ["C", "G"])

    def test_split_ldscore_table_requires_alleles_in_allele_aware_mode(self):
        ldscore_table = pd.DataFrame(
            {
                "CHR": ["1"],
                "SNP": ["rs1"],
                "POS": [10],
                "regression_ld_scores": [1.0],
                "base": [2.0],
            }
        )

        with self.assertRaisesRegex(LDSCInternalError, "A1/A2"):
            ldscore_workflow._split_ldscore_table(
                ldscore_table,
                baseline_columns=["base"],
                query_columns=[],
                snp_identifier="chr_pos_allele_aware",
            )

    def test_kernel_parse_annotation_file_preserves_optional_alleles_as_metadata(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "alleles.annot"
            path.write_text(
                "CHR\tBP\tSNP\tCM\tA1\tA2\tbase_a\n"
                "1\t10\trs1\t0.1\tA\tC\t1\n",
                encoding="utf-8",
            )

            metadata, annotations = kernel_ldscore.parse_annotation_file(str(path))

        self.assertEqual(list(metadata.columns), ["CHR", "POS", "SNP", "CM", "A1", "A2"])
        self.assertEqual(metadata["A1"].tolist(), ["A"])
        self.assertEqual(metadata["A2"].tolist(), ["C"])
        self.assertEqual(list(annotations.columns), ["base_a"])

    def test_kernel_parse_annotation_file_accepts_file_without_cm_and_ignores_maf(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "no_cm.annot"
            path.write_text("CHR\tBP\tSNP\tMAF\tbase_a\n1\t10\trs1\t0.3\t1\n", encoding="utf-8")
            metadata, annotations = kernel_ldscore.parse_annotation_file(str(path))
        self.assertIn("CM", metadata.columns)  # placeholder retained
        self.assertTrue(metadata["CM"].isna().all())
        self.assertNotIn("MAF", metadata.columns)  # MAF never carried
        self.assertEqual(list(annotations.columns), ["base_a"])

    def test_kernel_parse_annotation_file_rejects_single_allele_column(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "one_allele.annot"
            path.write_text(
                "CHR\tBP\tSNP\tCM\tA2\tbase_a\n"
                "1\t10\trs1\t0.1\tC\t1\n",
                encoding="utf-8",
            )

            with self.assertRaisesRegex(LDSCInputError, "only one allele column"):
                kernel_ldscore.parse_annotation_file(str(path))

    def test_kernel_parse_annotation_file_sorts_and_keeps_values_aligned(self):
        # Rows are intentionally out of genomic order; the parser must sort by
        # position AND keep each SNP glued to its own annotation value. The
        # annotation value equals the SNP's numeric suffix, so a misalignment
        # between metadata and values would be visible.
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "shuffled.annot"
            path.write_text(
                "CHR\tBP\tSNP\tbase_a\n"
                "1\t30\trs3\t3\n"
                "1\t10\trs1\t1\n"
                "1\t20\trs2\t2\n",
                encoding="utf-8",
            )
            metadata, annotations = kernel_ldscore.parse_annotation_file(str(path))
        self.assertEqual(metadata["SNP"].tolist(), ["rs1", "rs2", "rs3"])
        self.assertEqual(metadata["POS"].tolist(), [10, 20, 30])
        self.assertEqual(annotations["base_a"].tolist(), [1.0, 2.0, 3.0])

    def test_kernel_combine_annotation_groups_aligns_allele_free_annotations_in_allele_aware_mode(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            base = tmpdir / "base.annot"
            query = tmpdir / "query.annot"
            rows = "CHR\tBP\tSNP\tCM\t{column}\n1\t10\trs1\t0.1\t{first}\n1\t20\trs2\t0.2\t{second}\n"
            base.write_text(rows.format(column="base_a", first=1, second=0), encoding="utf-8")
            query.write_text(rows.format(column="query_a", first=0, second=1), encoding="utf-8")

            bundle = kernel_ldscore.combine_annotation_groups(
                baseline_files=(str(base),),
                query_files=(str(query),),
                chrom="1",
                identifier_mode="chr_pos_allele_aware",
            )

        self.assertIsNotNone(bundle)
        self.assertEqual(bundle.baseline_columns, ["base_a"])
        self.assertEqual(bundle.query_columns, ["query_a"])
        self.assertNotIn("A1", bundle.metadata.columns)

    def test_kernel_combine_annotation_groups_promotes_later_query_alleles_after_alignment(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            base = tmpdir / "base.annot"
            query = tmpdir / "query.annot"
            base.write_text(
                "CHR\tBP\tSNP\tCM\tbase_a\n"
                "1\t10\trs1\t0.1\t1\n",
                encoding="utf-8",
            )
            query.write_text(
                "CHR\tBP\tSNP\tCM\tA1\tA2\tquery_a\n"
                "1\t10\trs1\t0.1\tA\tG\t0\n",
                encoding="utf-8",
            )

            bundle = kernel_ldscore.combine_annotation_groups(
                baseline_files=(str(base),),
                query_files=(str(query),),
                chrom="1",
                identifier_mode="chr_pos_allele_aware",
            )

        self.assertIsNotNone(bundle)
        self.assertEqual(bundle.metadata["A1"].tolist(), ["A"])
        self.assertEqual(bundle.metadata["A2"].tolist(), ["G"])
        self.assertEqual(bundle.baseline_columns, ["base_a"])
        self.assertEqual(bundle.query_columns, ["query_a"])

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
                "--bed-padding-bp",
                "50000",
            ]
        )
        self.assertEqual(args.query_annot_bed_sources, "query.bed")
        self.assertEqual(args.bed_padding_bp, 50000)

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

    def test_compute_counts_common_mask_is_inclusive_greater_equal(self):
        metadata = pd.DataFrame({"MAF": [0.05, 0.0500001, 0.2, 0.04]})
        annotations = pd.DataFrame({"base": [1.0, 1.0, 1.0, 1.0], "cat": [1.0, 0.0, 1.0, 1.0]})

        M, M_5_50 = kernel_ldscore.compute_counts(metadata, annotations, common_maf_min=0.05)

        # MAF == 0.05 is included (inclusive >=), so rows 0, 1 and 2 are common.
        np.testing.assert_allclose(M, [4.0, 3.0])
        np.testing.assert_allclose(M_5_50, [3.0, 2.0])

    def test_count_config_records_inclusive_operator(self):
        from ldsc.ldscore_calculator import _count_config_from_ldscore_config
        from ldsc.config import LDScoreConfig

        cfg = _count_config_from_ldscore_config(LDScoreConfig(ld_wind_cm=1.0, common_maf_min=0.05))
        self.assertEqual(cfg["common_reference_snp_maf_operator"], ">=")

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

    def test_frequency_metadata_duplicate_keys_drop_all_cluster_values(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            freq_path = Path(tmpdir) / "freq.tsv"
            freq_path.write_text(
                "CHR\tPOS\tSNP\tCM\tMAF\n"
                "1\t10\trs1\t0.1\t0.2\n"
                "1\t10\trs1\t9.9\t0.4\n"
                "1\t20\trs2\t0.2\t0.3\n",
                encoding="utf-8",
            )
            metadata = pd.DataFrame(
                {
                    "CHR": ["1", "1"],
                    "POS": [10, 20],
                    "SNP": ["rs1", "rs2"],
                    "CM": [np.nan, np.nan],
                }
            )

            with self.assertLogs("LDSC.ldscore", level="WARNING"):
                merged = kernel_ldscore.merge_frequency_metadata(
                    metadata,
                    Namespace(frqfile=str(freq_path)),
                    chrom="1",
                    identifier_mode="chr_pos",
                )

        self.assertTrue(pd.isna(merged.loc[0, "CM"]))
        self.assertTrue(pd.isna(merged.loc[0, "MAF"]))
        self.assertEqual(merged.loc[1, "CM"], 0.2)
        self.assertEqual(merged.loc[1, "MAF"], 0.3)

    def test_allele_aware_frequency_metadata_without_alleles_drops_duplicate_base_key_cluster_values(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            freq_path = Path(tmpdir) / "freq.tsv"
            freq_path.write_text(
                "CHR\tPOS\tSNP\tCM\tMAF\n"
                "1\t10\trs1\t0.1\t0.2\n"
                "1\t10\trs1\t9.9\t0.4\n"
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

            with self.assertLogs("LDSC.ldscore", level="WARNING"):
                merged = kernel_ldscore.merge_frequency_metadata(
                    metadata,
                    Namespace(frqfile=str(freq_path)),
                    chrom="1",
                    identifier_mode="chr_pos_allele_aware",
                )

        self.assertTrue(pd.isna(merged.loc[0, "CM"]))
        self.assertTrue(pd.isna(merged.loc[0, "MAF"]))
        self.assertEqual(merged.loc[1, "CM"], 0.2)
        self.assertEqual(merged.loc[1, "MAF"], 0.3)

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

    def test_plink_compute_rejects_unsorted_window_positions(self):
        # geno_meta inherits the annotation bundle's row order via keep_snps, so a
        # bundle whose positions are not coordinate-sorted must trip the window-sort
        # guard before get_block_lefts runs (otherwise LD windows are silently wrong).
        with tempfile.TemporaryDirectory() as tmpdir:
            prefix = self._copy_plink_fixture_with_distinct_fids(Path(tmpdir))
            wb = self._build_annotation_bundle(prefix)
            order = np.argsort(wb.metadata["POS"].to_numpy(), kind="mergesort")[::-1]
            bundle = kernel_ldscore.AnnotationBundle(
                metadata=wb.metadata.iloc[order].reset_index(drop=True),
                annotations=wb.baseline_annotations.iloc[order].reset_index(drop=True),
                baseline_columns=wb.baseline_columns,
                query_columns=wb.query_columns,
            )
            args = Namespace(
                bfile=str(prefix), keep=None, maf_min=None, maf=None,
                ld_wind_snps=10, ld_wind_kb=None, ld_wind_cm=None,
                yes_really=True, snp_batch_size=50, common_maf_min=0.05,
                snp_identifier="rsid",
            )
            with self.assertRaisesRegex(LDSCInternalError, "non-decreasing order"):
                kernel_ldscore.compute_chrom_from_plink("1", bundle, args, None)

    def _plink_cm_args(self, prefix: Path, genetic_map) -> Namespace:
        return Namespace(
            bfile=str(prefix), keep=None, maf_min=None, maf=None,
            ld_wind_snps=None, ld_wind_kb=None, ld_wind_cm=1.0,
            yes_really=True, snp_batch_size=50, common_maf_min=0.05,
            snp_identifier="rsid", genetic_map=genetic_map,
        )

    def test_plink_cm_window_rejects_all_zero_bim_cm_without_map(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            prefix = self._copy_plink_fixture_with_distinct_fids(Path(tmpdir))
            wb = self._build_annotation_bundle(prefix)
            bundle = kernel_ldscore.AnnotationBundle(
                metadata=wb.metadata, annotations=wb.baseline_annotations,
                baseline_columns=wb.baseline_columns, query_columns=wb.query_columns,
            )
            # The fixture .bim CM is all zero; without a map this must error even under
            # --yes-really (yes_really=True below).
            with self.assertRaisesRegex(LDSCInputError, "uninformative|all zero"):
                kernel_ldscore.compute_chrom_from_plink("1", bundle, self._plink_cm_args(prefix, None), None)

    def test_plink_cm_window_uses_genetic_map(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            prefix = self._copy_plink_fixture_with_distinct_fids(Path(tmpdir))
            wb = self._build_annotation_bundle(prefix)
            bundle = kernel_ldscore.AnnotationBundle(
                metadata=wb.metadata, annotations=wb.baseline_annotations,
                baseline_columns=wb.baseline_columns, query_columns=wb.query_columns,
            )
            gmap = pd.DataFrame({"CHR": ["1", "1"], "POS": [1, 8], "CM": [0.0, 4.0]})
            result = kernel_ldscore.compute_chrom_from_plink("1", bundle, self._plink_cm_args(prefix, gmap), None)
            self.assertGreater(len(result.metadata), 0)
            # interpolated CM is non-degenerate -> a real cM window, not whole-chromosome
            self.assertGreaterEqual(result.metadata["CM"].nunique(), 2)

    def test_compute_chromosome_writes_ref_metadata_sidecar(self):
        from ldsc._kernel.ref_panel import RefPanelLoader
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            prefix = self._copy_plink_fixture_with_distinct_fids(tmpdir)
            wb = self._build_annotation_bundle(prefix)
            gc = GlobalConfig(snp_identifier="rsid")
            ref_panel = RefPanelLoader(gc).load(RefPanelConfig(backend="plink", plink_prefix=str(prefix)))
            out = tmpdir / "out"
            ldscore_workflow.LDScoreCalculator().compute_chromosome(
                chrom="1", annotation_bundle=wb, ref_panel=ref_panel,
                ldscore_config=LDScoreConfig(ld_wind_snps=10, export_ref_metadata=True, whole_chromosome_ok=True),
                global_config=gc, regression_snps=None, export_dir=str(out),
            )
            sidecar = out / "ref_metadata" / "chr1_meta.tsv.gz"
            self.assertTrue(sidecar.exists())
            df = pd.read_csv(sidecar, sep="\t")
            self.assertEqual(list(df.columns), ["CHR", "POS", "SNP", "A1", "A2", "CM", "MAF"])
            self.assertTrue(df["MAF"].notna().any())  # genotype-derived MAF present

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
            self.assertIn("metadata", result.output_paths)
            self.assertIn("baseline", result.output_paths)
            self.assertIn("query", result.output_paths)
            self.assertNotIn("log", result.output_paths)
            self.assertTrue(Path(result.output_paths["metadata"]).exists())
            self.assertTrue((tmpdir / "ldscore_result" / "diagnostics" / "ldscore.log").exists())
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

            def _compute(_self, chrom, annotation_bundle, ref_panel, ldscore_config, global_config, regression_snps=None, export_dir=None):
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
            self.assertFalse((output_dir / "metadata.json").exists())
            self.assertFalse((output_dir / "diagnostics" / "ldscore.log").exists())

    def test_run_ldscore_from_args_ignores_legacy_root_log_without_overwrite(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            output_dir = tmpdir / "ldscore_result"
            output_dir.mkdir()
            legacy_log = output_dir / "ldscore.log"
            legacy_log.write_text("legacy log\n", encoding="utf-8")
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
                return_value=self.make_baseline_only_chrom_result("1", 10, 1.0, 5.0),
            ):
                ldscore_workflow.run_ldscore_from_args(args)

            self.assertEqual(legacy_log.read_text(encoding="utf-8"), "legacy log\n")
            self.assertTrue((output_dir / "metadata.json").exists())
            self.assertTrue((output_dir / "diagnostics" / "ldscore.log").exists())

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

            self.assertTrue((output_dir / "metadata.json").exists())
            self.assertTrue((output_dir / "diagnostics" / "ldscore.log").exists())
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

        with self.assertRaisesRegex(LDSCUsageError, "query annotations without baseline annotations.*all-ones `base`"):
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

        with self.assertRaisesRegex(LDSCUsageError, "query annotations without baseline annotations.*all-ones `base`"):
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

            def _compute(_self, chrom, annotation_bundle, ref_panel, ldscore_config, global_config, regression_snps=None, export_dir=None):
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

            with self.assertRaisesRegex(LDSCInputError, "could not infer the genome build"):
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
                with self.assertRaisesRegex(LDSCInputError, "conflicting genome-build evidence"):
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

        with self.assertRaisesRegex(LDSCUsageError, "keep-indivs-file.*parquet R2 mode"):
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

            def _compute_side_effect(_self, chrom, annotation_bundle, ref_panel, ldscore_config, global_config, regression_snps=None, export_dir=None):
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

    def test_compute_chromosome_uses_annotation_alleles_for_allele_aware_matching(self):
        annotation_bundle = AnnotationBundle(
            metadata=pd.DataFrame(
                {
                    "CHR": ["1"],
                    "SNP": ["rs1"],
                    "POS": [10],
                    "CM": [0.1],
                    "A1": ["A"],
                    "A2": ["C"],
                }
            ),
            baseline_annotations=pd.DataFrame({"base": [1.0]}),
            query_annotations=pd.DataFrame(index=pd.RangeIndex(1)),
            baseline_columns=["base"],
            query_columns=[],
            chromosomes=["1"],
            source_summary={},
            config_snapshot=None,
        )
        ref_panel = self.make_ref_panel_stub(backend="parquet_r2")
        ref_panel.load_metadata = mock.Mock(
            return_value=pd.DataFrame(
                {
                    "CHR": ["1"],
                    "SNP": ["rs1"],
                    "CM": [0.1],
                    "POS": [10],
                    "A1": ["A"],
                    "A2": ["G"],
                }
            )
        )

        with mock.patch.object(ldscore_workflow.kernel_ldscore, "compute_chrom_from_parquet") as patched_compute:
            with self.assertRaisesRegex(LDSCInputError, "retained no annotation SNPs"):
                ldscore_workflow.LDScoreCalculator().compute_chromosome(
                    chrom="1",
                    annotation_bundle=annotation_bundle,
                    ref_panel=ref_panel,
                    ldscore_config=LDScoreConfig(ld_wind_snps=10),
                    global_config=GlobalConfig(snp_identifier="chr_pos_allele_aware", genome_build="hg38"),
                )

        patched_compute.assert_not_called()

    def test_compute_chromosome_uses_later_query_alleles_for_allele_aware_matching(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            base = tmpdir / "base.annot"
            query = tmpdir / "query.annot"
            base.write_text(
                "CHR\tBP\tSNP\tCM\tbase\n"
                "1\t10\trs1\t0.1\t1\n",
                encoding="utf-8",
            )
            query.write_text(
                "CHR\tBP\tSNP\tCM\tA1\tA2\tquery\n"
                "1\t10\trs1\t0.1\tA\tG\t1\n",
                encoding="utf-8",
            )
            annotation_bundle = AnnotationBuilder(
                GlobalConfig(snp_identifier="chr_pos_allele_aware", genome_build="hg38")
            ).run(
                AnnotationBuildConfig(
                    baseline_annot_sources=(str(base),),
                    query_annot_sources=(str(query),),
                )
            )

        self.assertEqual(annotation_bundle.metadata["A2"].tolist(), ["G"])
        ref_panel = self.make_ref_panel_stub(backend="parquet_r2")
        ref_panel.load_metadata = mock.Mock(
            return_value=pd.DataFrame(
                {
                    "CHR": ["1"],
                    "SNP": ["rs1"],
                    "CM": [0.1],
                    "POS": [10],
                    "A1": ["A"],
                    "A2": ["C"],
                }
            )
        )

        with mock.patch.object(ldscore_workflow.kernel_ldscore, "compute_chrom_from_parquet") as patched_compute:
            with self.assertRaisesRegex(LDSCInputError, "retained no annotation SNPs"):
                ldscore_workflow.LDScoreCalculator().compute_chromosome(
                    chrom="1",
                    annotation_bundle=annotation_bundle,
                    ref_panel=ref_panel,
                    ldscore_config=LDScoreConfig(ld_wind_snps=10),
                    global_config=GlobalConfig(snp_identifier="chr_pos_allele_aware", genome_build="hg38"),
                )

        patched_compute.assert_not_called()

    def test_compute_chromosome_matches_cleaned_reference_panel_before_annotation_alignment(self):
        annotation_bundle = self.make_annotation_bundle(
            [("1", "rs1", 10), ("1", "rs2", 20)],
        )
        ref_panel = self.make_ref_panel_stub(backend="parquet_r2")
        ref_panel.load_metadata = mock.Mock(
            return_value=pd.DataFrame(
                {
                    "CHR": ["1"],
                    "SNP": ["rs2"],
                    "CM": [0.2],
                    "POS": [20],
                    "A1": ["A"],
                    "A2": ["C"],
                }
            )
        )

        def _compute_side_effect(chrom, bundle, args, regression_snps):
            self.assertEqual(bundle.metadata["SNP"].tolist(), ["rs2"])
            return ldscore_workflow._LegacyChromResult(
                chrom=chrom,
                metadata=bundle.metadata.assign(MAF=0.2),
                ld_scores=np.array([[1.0]], dtype=np.float32),
                w_ld=np.array([[2.0]], dtype=np.float32),
                M=np.array([1.0]),
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

        self.assertEqual(result.baseline_table["SNP"].tolist(), ["rs2"])

    def test_compute_chromosome_drops_reference_multiallelic_sites_before_allele_free_annotation_matching(self):
        annotation_bundle = self.make_annotation_bundle(
            [("1", "rs_multi", 10), ("1", "rs_keep", 20)],
        )
        ref_panel = self.make_ref_panel_stub(backend="parquet_r2")
        ref_panel.load_metadata = mock.Mock(
            return_value=pd.DataFrame(
                {
                    "CHR": ["1", "1", "1"],
                    "SNP": ["rs_multi_c", "rs_multi_g", "rs_keep"],
                    "CM": [0.1, 0.1, 0.2],
                    "POS": [10, 10, 20],
                    "A1": ["A", "A", "A"],
                    "A2": ["C", "G", "C"],
                }
            )
        )

        def _compute_side_effect(chrom, bundle, args, regression_snps):
            self.assertEqual(bundle.metadata["SNP"].tolist(), ["rs_keep"])
            self.assertEqual(bundle.metadata["A1"].tolist(), ["A"])
            self.assertEqual(bundle.metadata["A2"].tolist(), ["C"])
            self.assertEqual(bundle.annotations["base"].tolist(), [1.0])
            return ldscore_workflow._LegacyChromResult(
                chrom=chrom,
                metadata=bundle.metadata.assign(MAF=0.2),
                ld_scores=np.array([[1.0]], dtype=np.float32),
                w_ld=np.array([[2.0]], dtype=np.float32),
                M=np.array([1.0]),
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

        self.assertEqual(result.baseline_table["SNP"].tolist(), ["rs_keep"])

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


def _write_index_sidecar(tmp, df, *, gz: bool = True):
    """Write a panel sidecar (CHR,POS,SNP,A1,A2,CM,MAF) as chr1_meta.tsv.gz."""
    meta = Path(tmp) / "hg19" / "chr1_meta.tsv.gz"
    meta.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(meta, "wt") as h:
        h.write("# ldsc:schema_version=2\n")
        df.to_csv(h, sep="\t", index=False)
    return meta


class IndexBindingTest(unittest.TestCase):
    def _sidecar(self, tmp):
        df = pd.DataFrame({"CHR": ["1", "1", "1"], "POS": [10, 20, 30],
                           "SNP": ["rsA", "rsB", "rsC"], "A1": ["A", "C", "G"],
                           "A2": ["G", "T", "A"], "CM": [0.0, 0.0, 0.0], "MAF": [0.2, 0.3, 0.4]})
        return _write_index_sidecar(tmp, df), df

    def test_load_full_panel_sidecar_reads_canonical_columns(self):
        from ldsc._kernel import ldscore as ls
        with tempfile.TemporaryDirectory() as tmp:
            meta, df = self._sidecar(tmp)
            r2 = meta.with_name("chr1_r2.parquet")
            loaded = ls._load_full_panel_sidecar(str(r2))
            self.assertEqual(list(loaded["CHR"].astype(str)), ["1", "1", "1"])
            self.assertEqual(list(loaded["POS"].astype(int)), [10, 20, 30])

    def test_validate_binding_raises_on_hash_mismatch(self):
        from ldsc._kernel import ldscore as ls
        from ldsc._kernel.snp_identity import sidecar_identity_sha256
        with tempfile.TemporaryDirectory() as tmp:
            meta, df = self._sidecar(tmp)
            good = sidecar_identity_sha256(df)
            ls._validate_index_binding(df, n_snps=3, identity_hash=good, context="t")  # no raise
            with self.assertRaisesRegex(LDSCInputError, "identity hash"):
                ls._validate_index_binding(df, n_snps=3, identity_hash="0" * 64, context="t")
            with self.assertRaisesRegex(LDSCInputError, "n_snps"):
                ls._validate_index_binding(df, n_snps=999, identity_hash=good, context="t")


class IndexRemapTest(unittest.TestCase):
    def test_remap_maps_build_rows_to_retained_indices_chr_pos(self):
        from ldsc._kernel.ldscore import build_index_remap

        full = pd.DataFrame({"CHR": ["1"] * 4, "POS": [10, 20, 30, 40],
                             "SNP": ["a", "b", "c", "d"], "A1": list("ACGT"), "A2": list("GTAC")})
        # retained = rows 2 and 0 of the panel, in matrix order [pos30, pos10]
        retained = full.iloc[[2, 0]].reset_index(drop=True)
        remap, retained_build_idx = build_index_remap(full, retained, "chr_pos")
        # build row 2 (pos30) -> retained idx 0 ; build row 0 (pos10) -> retained idx 1 ; others -1
        self.assertEqual(remap.tolist(), [1, -1, 0, -1])
        # retained_build_idx[matrix_idx] = build row
        self.assertEqual(retained_build_idx.tolist(), [2, 0])

    def test_remap_allele_aware_uses_alleles(self):
        from ldsc._kernel.ldscore import build_index_remap

        # Two multi-allelic variants at the same position with distinct,
        # non-strand-complement allele sets: {A,G} vs {A,C}.
        full = pd.DataFrame({"CHR": ["1", "1"], "POS": [10, 10],
                             "SNP": ["a", "b"], "A1": ["A", "A"], "A2": ["G", "C"]})
        retained = full.iloc[[1]].reset_index(drop=True)  # the A/C variant at pos10
        remap, _ = build_index_remap(full, retained, "chr_pos_allele_aware")
        self.assertEqual(remap.tolist(), [-1, 0])


class IndexReaderDecodeTest(unittest.TestCase):
    def _make_panel(self, tmp):
        from ldsc._kernel import ref_panel_builder as kb
        from ldsc._kernel.snp_identity import sidecar_identity_sha256

        full = pd.DataFrame({"CHR": ["1"] * 4, "POS": [10, 20, 30, 40],
                             "SNP": ["a", "b", "c", "d"], "A1": list("ACGT"),
                             "A2": list("GTAC"), "CM": [0.0] * 4, "MAF": [0.3] * 4})
        meta = _write_index_sidecar(tmp, full)
        r2 = meta.with_name("chr1_r2.parquet")
        rows = [
            {"i": 0, "j": 1, "R2": 0.5, "sign": "+"},
            {"i": 0, "j": 2, "R2": 0.2, "sign": "+"},
            {"i": 1, "j": 2, "R2": 0.4, "sign": "-"},
            {"i": 2, "j": 3, "R2": 0.9, "sign": "+"},
        ]
        kb.write_r2_parquet(pair_chunks=dict_chunks(rows), path=r2, genome_build="hg19", n_samples=100,
                            snp_identifier="chr_pos", min_r2=0.0, n_snps=4,
                            sidecar_identity_sha256=sidecar_identity_sha256(full))
        return r2, full

    @unittest.skipUnless(_HAS_PYARROW, "pyarrow required")
    def test_streaming_drops_unretained_endpoint(self):
        from ldsc._kernel.ldscore import SortedR2BlockReader
        with tempfile.TemporaryDirectory() as tmp:
            r2, full = self._make_panel(tmp)
            # retained = panel rows 0,1,2 (drop pos40); matrix order = [10,20,30]
            retained = full.iloc[[0, 1, 2]][["CHR", "POS", "SNP", "A1", "A2", "CM", "MAF"]].reset_index(drop=True)
            reader = SortedR2BlockReader(
                paths=[str(r2)], chrom="1", metadata=retained,
                identifier_mode="chr_pos", r2_bias_mode="unbiased",
                r2_sample_size=None, genome_build="hg19",
            )
            seen = {}
            for i, j, r2v in reader.iter_all_pairs():
                for a, b, v in zip(i.tolist(), j.tolist(), r2v.tolist()):
                    seen[(a, b)] = v
            # pair (2,3) touches the dropped SNP (pos40) -> excluded by remap
            self.assertEqual(sorted(seen), [(0, 1), (0, 2), (1, 2)])
            # off-diagonal R2 is int16-quantized: tolerance is the half-step (1.5e-5)
            self.assertAlmostEqual(seen[(0, 1)], 0.5, delta=2e-5)
            self.assertAlmostEqual(seen[(0, 2)], 0.2, delta=2e-5)
            self.assertAlmostEqual(seen[(1, 2)], 0.4, delta=2e-5)


class R2PanelWindowValidationTest(unittest.TestCase):
    @unittest.skipUnless(_HAS_PYARROW, "pyarrow required")
    def test_ldscore_accepts_equal_recorded_r2_panel_window(self):
        from ldsc._kernel import ref_panel_builder as kb
        from ldsc._kernel.snp_identity import sidecar_identity_sha256

        panel = pd.DataFrame(
            {
                "CHR": ["1", "1"],
                "POS": [100, 200],
                "SNP": ["rs1", "rs2"],
                "A1": ["A", "A"],
                "A2": ["C", "C"],
                "CM": [0.0, 1.0],
                "MAF": [0.3, 0.3],
            }
        )
        with tempfile.TemporaryDirectory() as tmp:
            meta = _write_index_sidecar(tmp, panel)
            r2 = meta.with_name("chr1_r2.parquet")
            kb.write_r2_parquet(
                pair_chunks=dict_chunks([]),
                path=r2,
                genome_build="hg19",
                n_samples=100,
                snp_identifier="chr_pos",
                min_r2=0.0,
                n_snps=2,
                sidecar_identity_sha256=sidecar_identity_sha256(panel),
                ld_window_mode="cm",
                ld_window_value=1.0,
            )

            args = Namespace(
                ld_wind_snps=None,
                ld_wind_kb=None,
                ld_wind_cm=1.0,
                r2_table=str(r2),
                snp_identifier="chr_pos",
            )

            kernel_ldscore.validate_ldscore_window_within_r2_panel_window(
                args,
                parquet_paths=[str(r2)],
                chrom="1",
            )

    @unittest.skipUnless(_HAS_PYARROW, "pyarrow required")
    def test_ldscore_rejects_window_wider_than_recorded_r2_panel_window(self):
        from ldsc._kernel import ref_panel_builder as kb
        from ldsc._kernel.snp_identity import sidecar_identity_sha256

        panel = pd.DataFrame(
            {
                "CHR": ["1"] * 3,
                "POS": [100, 200, 300],
                "SNP": ["rs1", "rs2", "rs3"],
                "A1": ["A", "A", "A"],
                "A2": ["C", "C", "C"],
                "CM": [0.0, 0.1, 0.2],
                "MAF": [0.3, 0.3, 0.3],
            }
        )
        with tempfile.TemporaryDirectory() as tmp:
            meta = _write_index_sidecar(tmp, panel)
            r2 = meta.with_name("chr1_r2.parquet")
            kb.write_r2_parquet(
                pair_chunks=dict_chunks([{"i": 0, "j": 1, "R2": 0.5, "sign": "+"}]),
                path=r2,
                genome_build="hg19",
                n_samples=100,
                snp_identifier="chr_pos",
                min_r2=0.0,
                n_snps=3,
                sidecar_identity_sha256=sidecar_identity_sha256(panel),
                ld_window_mode="kb",
                ld_window_value=0.1,
            )

            args = Namespace(
                ld_wind_snps=None,
                ld_wind_kb=1.0,
                ld_wind_cm=None,
                r2_table=str(r2),
            )
            with self.assertRaisesRegex(
                LDSCUsageError,
                "user-requested LD window `--ld-wind-kb 1.0` is wider than the input R2 parquet panel window `--ld-wind-kb 0.1`",
            ):
                kernel_ldscore.validate_ldscore_window_within_r2_panel_window(
                    args,
                    parquet_paths=[str(r2)],
                    chrom="1",
                )

    @unittest.skipUnless(_HAS_PYARROW, "pyarrow required")
    def test_ldscore_rejects_any_strictly_larger_same_mode_window(self):
        from ldsc._kernel import ref_panel_builder as kb
        from ldsc._kernel.snp_identity import sidecar_identity_sha256

        panel = pd.DataFrame(
            {
                "CHR": ["1", "1"],
                "POS": [100, 200],
                "SNP": ["rs1", "rs2"],
                "A1": ["A", "A"],
                "A2": ["C", "C"],
                "CM": [0.0, 1.0],
                "MAF": [0.3, 0.3],
            }
        )
        with tempfile.TemporaryDirectory() as tmp:
            meta = _write_index_sidecar(tmp, panel)
            r2 = meta.with_name("chr1_r2.parquet")
            kb.write_r2_parquet(
                pair_chunks=dict_chunks([]),
                path=r2,
                genome_build="hg19",
                n_samples=100,
                snp_identifier="chr_pos",
                min_r2=0.0,
                n_snps=2,
                sidecar_identity_sha256=sidecar_identity_sha256(panel),
                ld_window_mode="cm",
                ld_window_value=1.0,
            )

            args = Namespace(
                ld_wind_snps=None,
                ld_wind_kb=None,
                ld_wind_cm=float(np.nextafter(1.0, 2.0)),
                r2_table=str(r2),
                snp_identifier="chr_pos",
            )
            with self.assertRaisesRegex(LDSCUsageError, "is wider than the input R2 parquet panel window"):
                kernel_ldscore.validate_ldscore_window_within_r2_panel_window(
                    args,
                    parquet_paths=[str(r2)],
                    chrom="1",
                )

    @unittest.skipUnless(_HAS_PYARROW, "pyarrow required")
    def test_ldscore_warns_when_requested_window_mode_differs_from_r2_panel_mode(self):
        from ldsc._kernel import ref_panel_builder as kb
        from ldsc._kernel.snp_identity import sidecar_identity_sha256

        panel = pd.DataFrame(
            {
                "CHR": ["1", "1"],
                "POS": [100, 200],
                "SNP": ["rs1", "rs2"],
                "A1": ["A", "A"],
                "A2": ["C", "C"],
                "CM": [0.0, 1.0],
                "MAF": [0.3, 0.3],
            }
        )
        with tempfile.TemporaryDirectory() as tmp:
            meta = _write_index_sidecar(tmp, panel)
            r2 = meta.with_name("chr1_r2.parquet")
            kb.write_r2_parquet(
                pair_chunks=dict_chunks([]),
                path=r2,
                genome_build="hg19",
                n_samples=100,
                snp_identifier="chr_pos",
                min_r2=0.0,
                n_snps=2,
                sidecar_identity_sha256=sidecar_identity_sha256(panel),
                ld_window_mode="cm",
                ld_window_value=1.0,
            )

            args = Namespace(
                ld_wind_snps=None,
                ld_wind_kb=100.0,
                ld_wind_cm=None,
                r2_table=str(r2),
                snp_identifier="chr_pos",
            )
            with self.assertLogs("LDSC.ldscore", level="WARNING") as caught:
                kernel_ldscore.validate_ldscore_window_within_r2_panel_window(
                    args,
                    parquet_paths=[str(r2)],
                    chrom="1",
                )

        message = "\n".join(caught.output)
        self.assertIn("R2 parquet panel was built with `--ld-wind-cm 1.0`", message)
        self.assertIn("ldscore run requested `--ld-wind-kb 100.0`", message)
        self.assertIn("may require SNP pairs not stored in the parquet", message)


class IndexCrossModeParityTest(unittest.TestCase):
    """One index parquet must produce identical LD scores in all four modes."""

    @unittest.skipUnless(_HAS_PYARROW, "pyarrow required")
    def test_ld_scores_identical_across_modes_from_one_index_parquet(self):
        from ldsc._kernel import ref_panel_builder as kb
        from ldsc._kernel.snp_identity import sidecar_identity_sha256

        panel = pd.DataFrame({
            "CHR": ["1"] * 4, "POS": [100, 120, 140, 160], "SNP": ["rs1", "rs2", "rs3", "rs4"],
            "A1": ["A", "A", "A", "A"], "A2": ["C", "G", "C", "G"], "CM": [0.0] * 4, "MAF": [0.3] * 4,
        })
        pairs = [
            {"i": 0, "j": 1, "R2": 0.4, "sign": "+"},
            {"i": 0, "j": 2, "R2": 0.2, "sign": "-"},
            {"i": 1, "j": 2, "R2": 0.6, "sign": "+"},
            {"i": 2, "j": 3, "R2": 0.5, "sign": "+"},
        ]
        # Dense oracle: symmetric R2 matrix with unit diagonal, times all-ones annotation.
        oracle = np.array([
            [1.0, 0.4, 0.2, 0.0],
            [0.4, 1.0, 0.6, 0.0],
            [0.2, 0.6, 1.0, 0.5],
            [0.0, 0.0, 0.5, 1.0],
        ], dtype=np.float64)
        expected = oracle.sum(axis=1)  # [1.6, 2.0, 2.3, 1.5]

        with tempfile.TemporaryDirectory() as tmp:
            meta = _write_index_sidecar(tmp, panel)
            r2 = meta.with_name("chr1_r2.parquet")
            kb.write_r2_parquet(pair_chunks=dict_chunks(pairs), path=r2, genome_build="hg19", n_samples=100,
                                snp_identifier="chr_pos", min_r2=0.0, n_snps=4,
                                sidecar_identity_sha256=sidecar_identity_sha256(panel))
            mode_scores = {}
            for mode in ("rsid", "chr_pos", "rsid_allele_aware", "chr_pos_allele_aware"):
                reader = kernel_ldscore.SortedR2BlockReader(
                    paths=[str(r2)], chrom="1",
                    metadata=panel[["CHR", "POS", "SNP", "A1", "A2", "CM", "MAF"]].copy(),
                    identifier_mode=mode, r2_bias_mode="unbiased",
                    r2_sample_size=None, genome_build="hg19",
                )
                scores = kernel_ldscore.ld_score_streaming_from_r2_reader(
                    block_left=np.zeros(4, dtype=np.int64),
                    annot=np.ones((4, 1), dtype=np.float32), block_reader=reader,
                )
                mode_scores[mode] = scores[:, 0]
            # All four modes read one quantized parquet -> bit-identical to each other.
            first = mode_scores["rsid"]
            for mode, s in mode_scores.items():
                np.testing.assert_array_equal(s, first, err_msg=f"mode {mode} differs from rsid")
            # Dense oracle uses exact R2; int16 quantization adds <= a few half-steps.
            np.testing.assert_allclose(first, expected, rtol=0, atol=5e-5,
                                       err_msg="quantized LD scores drifted from dense oracle")


class RawSchemaRejectedTest(unittest.TestCase):
    @unittest.skipUnless(_HAS_PYARROW, "pyarrow required")
    def test_legacy_raw_schema_is_rejected(self):
        from ldsc._kernel.ldscore import SortedR2BlockReader
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "chr1_r2.parquet"
            _write_legacy_r2_parquet(path)  # writes old 10-col canonical schema
            meta = pd.DataFrame({"CHR": ["1"], "POS": [10], "SNP": ["a"], "A1": ["A"], "A2": ["G"]})
            with self.assertRaisesRegex(LDSCInputError, "index-format|build-ref-panel"):
                SortedR2BlockReader(paths=[str(path)], chrom="1", metadata=meta,
                                    identifier_mode="chr_pos", r2_bias_mode="unbiased",
                                    r2_sample_size=None, genome_build="hg19")


class R2TransformClipTest(unittest.TestCase):
    def test_raw_transform_upper_clips_at_one(self):
        from ldsc._kernel.ldscore import SortedR2BlockReader
        reader = SortedR2BlockReader.__new__(SortedR2BlockReader)
        reader.r2_bias_mode = "raw"
        reader.r2_sample_size = 100.0
        # raw r2 = 1.0 -> corrected 1.0 - 0/98 = 1.0; values >1 must clip to 1.0.
        out = reader._transform_r2(np.array([1.0, 1.05, 0.0], dtype=np.float32))
        self.assertLessEqual(float(out.max()), 1.0)
        self.assertEqual(float(out[0]), 1.0)


class R2DequantizationTest(unittest.TestCase):
    @unittest.skipUnless(_HAS_PYARROW, "pyarrow required")
    def test_reader_dequantizes_int16_panel_to_float32(self):
        from ldsc._kernel import ref_panel_builder as kb
        from ldsc._kernel import ldscore as kernel_ldscore
        from ldsc._kernel.snp_identity import sidecar_identity_sha256
        panel = pd.DataFrame({"CHR": ["1"] * 3, "POS": [100, 120, 140],
                              "SNP": ["rs1", "rs2", "rs3"], "A1": ["A"] * 3, "A2": ["C"] * 3,
                              "CM": [0.0] * 3, "MAF": [0.3] * 3})
        pairs = [{"i": 0, "j": 1, "R2": 1.0, "sign": "+"},
                 {"i": 0, "j": 2, "R2": 0.2, "sign": "-"}]
        with tempfile.TemporaryDirectory() as tmp:
            meta = _write_index_sidecar(tmp, panel)
            r2 = meta.with_name("chr1_r2.parquet")
            kb.write_r2_parquet(pair_chunks=dict_chunks(pairs), path=r2, genome_build="hg19", n_samples=100,
                                snp_identifier="chr_pos", min_r2=0.0, n_snps=3,
                                sidecar_identity_sha256=sidecar_identity_sha256(panel))
            reader = kernel_ldscore.SortedR2BlockReader(
                paths=[str(r2)], chrom="1", metadata=panel.copy(),
                identifier_mode="chr_pos", r2_bias_mode="unbiased",
                r2_sample_size=None, genome_build="hg19")
            self.assertEqual(reader._r2_scale, 32767.0)
            decoded = reader._decode_index_row_group(0)
            self.assertEqual(decoded.r2.dtype, np.float32)
            # endpoint exact; 0.2 within half-step
            got = dict(zip(zip(decoded.i.tolist(), decoded.j.tolist()), decoded.r2.tolist()))
            self.assertEqual(got[(0, 1)], 1.0)
            self.assertAlmostEqual(got[(0, 2)], 0.2, delta=2e-5)

    @unittest.skipUnless(_HAS_PYARROW, "pyarrow required")
    def test_reader_leaves_float32_panel_unscaled(self):
        # A hand-built float32 R2 index parquet must read with _r2_scale is None.
        import pyarrow as pa
        import pyarrow.parquet as pq
        from ldsc._kernel import ldscore as kernel_ldscore
        from ldsc._kernel.snp_identity import sidecar_identity_sha256
        panel = pd.DataFrame({"CHR": ["1", "1"], "POS": [100, 120], "SNP": ["rs1", "rs2"],
                              "A1": ["A", "A"], "A2": ["C", "G"], "CM": [0.0, 0.0], "MAF": [0.3, 0.3]})
        with tempfile.TemporaryDirectory() as tmp:
            meta = _write_index_sidecar(tmp, panel)
            r2 = meta.with_name("chr1_r2.parquet")
            schema = pa.schema(
                [("IDX_1", pa.int32()), ("IDX_2", pa.int32()),
                 ("R2", pa.float32()), ("SIGN", pa.bool_())]
            ).with_metadata({
                b"ldsc:schema_version": b"1", b"ldsc:artifact_type": b"ref_panel_r2",
                b"ldsc:snp_identifier": b"chr_pos", b"ldsc:genome_build": b"hg19",
                b"ldsc:sorted_by_build": b"hg19", b"ldsc:n_samples": b"100",
                b"ldsc:r2_bias": b"unbiased", b"ldsc:min_r2": b"0.0", b"ldsc:n_snps": b"2",
                b"ldsc:sidecar_identity_sha256": sidecar_identity_sha256(panel).encode(),
                b"ldsc:row_group_size": b"50000",
            })
            tbl = pa.table({"IDX_1": pa.array([0], pa.int32()), "IDX_2": pa.array([1], pa.int32()),
                            "R2": pa.array([0.5], pa.float32()), "SIGN": pa.array([True], pa.bool_())},
                           schema=schema)
            pq.write_table(tbl, str(r2))
            reader = kernel_ldscore.SortedR2BlockReader(
                paths=[str(r2)], chrom="1", metadata=panel.copy(),
                identifier_mode="chr_pos", r2_bias_mode="unbiased",
                r2_sample_size=None, genome_build="hg19")
            self.assertIsNone(reader._r2_scale)
            decoded = reader._decode_index_row_group(0)
            self.assertAlmostEqual(float(decoded.r2[0]), 0.5, delta=1e-7)


def test_ldscore_parser_region_flags_and_default():
    parser = build_parser()
    # Default excludes MHC + centromeres.
    default_args = parser.parse_args(["--output-dir", "o", "--plink-prefix", "p", "--ld-wind-cm", "1"])
    assert default_args.exclude_regions == "mhc-and-centromeres"
    # Single-choice enum; comma lists are rejected.
    args = parser.parse_args(
        ["--output-dir", "o", "--plink-prefix", "p", "--ld-wind-cm", "1",
         "--exclude-regions", "centromeres", "--exclude-regions-build", "hg19",
         "--exclude-regions-bed", "/tmp/a.bed,/tmp/b.bed"]
    )
    assert args.exclude_regions == "centromeres"
    assert args.exclude_regions_build == "hg19"
    with pytest.raises(SystemExit):
        parser.parse_args(["--output-dir", "o", "--plink-prefix", "p", "--ld-wind-cm", "1",
                           "--exclude-regions", "mhc,centromeres"])


def test_ldscore_region_build_inferred_for_chr_pos_required_for_rsid():
    # chr_pos: build is taken from the panel build ldscore operates in.
    assert _resolve_exclude_regions_build(
        Namespace(exclude_regions_build=None),
        GlobalConfig(snp_identifier="chr_pos_allele_aware", genome_build="hg19"),
        ("mhc", "centromeres"),
    ) == "hg19"
    # rsid: no genome build is available, so presets require an explicit build.
    with pytest.raises(LDSCUsageError, match="exclude-regions-build"):
        _resolve_exclude_regions_build(
            Namespace(exclude_regions_build=None),
            GlobalConfig(snp_identifier="rsid"),
            ("mhc",),
        )
    # An explicit build always wins, even in rsid mode.
    assert _resolve_exclude_regions_build(
        Namespace(exclude_regions_build="hg38"),
        GlobalConfig(snp_identifier="rsid"),
        ("mhc",),
    ) == "hg38"
    # No presets => no build needed.
    assert _resolve_exclude_regions_build(
        Namespace(exclude_regions_build=None),
        GlobalConfig(snp_identifier="rsid"),
        (),
    ) is None


_CHR22 = Path(__file__).resolve().parent / "fixtures" / "minimal_external_resources" / "plink" / "hm3_chr22_subset"


def _chr22_available() -> bool:
    return all(Path(str(_CHR22) + ext).exists() for ext in (".bed", ".bim", ".fam"))


def test_ldscore_end_to_end_drops_user_bed_snp(tmp_path):
    if not _chr22_available():
        pytest.skip("chr22 PLINK fixture unavailable; run tests/fixtures/generate_minimal_external_resources.py")
    set_global_config(GlobalConfig(snp_identifier="chr_pos_allele_aware", genome_build="hg38"))

    # Discover a real emitted POS from a no-exclusion panel load.
    base_meta = PlinkRefPanel(
        GlobalConfig(snp_identifier="chr_pos_allele_aware", genome_build="hg38"),
        RefPanelConfig(backend="plink", plink_prefix=str(_CHR22)),
    ).load_metadata("22")
    target_pos = int(base_meta["POS"].iloc[0])

    bed = tmp_path / "exclude.bed"
    bed.write_text(f"22\t{target_pos - 1}\t{target_pos}\n", encoding="utf-8")

    result = run_ldscore(
        plink_prefix=str(_CHR22),
        output_dir=str(tmp_path / "ld"),
        ld_wind_snps=10,
        exclude_regions_bed=(str(bed),),
    )
    assert target_pos not in set(result.baseline_table["POS"])
    # Unpartitioned (synthetic base-only) runs have a degenerate 1x1 overlap that
    # no regression consumes, so the sidecar is suppressed end-to-end.
    assert not (tmp_path / "ld" / "ldscore.overlap.parquet").exists()
    assert result.overlap is None


def test_ldscore_rejects_non_minor_maf():
    bad = pd.DataFrame({"CHR": ["22"], "SNP": ["rs1"], "POS": [1],
                        "A1": ["A"], "A2": ["G"], "MAF": [0.73]})
    with pytest.raises(LDSCInternalError, match="MAF"):
        ldscore_workflow._assert_canonical_maf(bad)


def test_assert_canonical_maf_passes_for_minor():
    ok = pd.DataFrame({"MAF": [0.0, 0.2, 0.5]})
    ldscore_workflow._assert_canonical_maf(ok)  # no raise


def test_merge_frequency_metadata_sidecar_is_authoritative(monkeypatch, tmp_path):
    """Reference-panel sidecar CM/MAF overwrite annotation-provided values."""
    metadata = pd.DataFrame(
        {"CHR": ["1"], "POS": [10], "SNP": ["rs1"], "CM": [99.0], "MAF": [0.99]}
    )
    sidecar = tmp_path / "chr1_meta.tsv.gz"
    pd.DataFrame(
        {"CHR": ["1"], "POS": [10], "SNP": ["rs1"], "CM": [0.5], "MAF": [0.2]}
    ).to_csv(sidecar, sep="\t", index=False)
    monkeypatch.setattr(
        kernel_ldscore, "resolve_frequency_files", lambda args, chrom=None: [str(sidecar)]
    )
    args = Namespace(frqfile=str(sidecar), frqfile_chr=None)
    merged = kernel_ldscore.merge_frequency_metadata(
        metadata, args, chrom="1", identifier_mode="rsid"
    )
    assert merged.loc[0, "CM"] == 0.5  # sidecar wins, not annotation 99.0
    assert merged.loc[0, "MAF"] == 0.2  # sidecar wins (folded)


def test_cm_window_consistent_across_backends():
    """Both backends feed reference CM into the same window builder, so identical
    reference CM yields identical block-lefts (the headline cross-backend invariant)."""
    import argparse
    cm = pd.Series([0.0, 0.4, 0.9, 1.6])
    md = pd.DataFrame({"CM": cm, "POS": [1, 2, 3, 4]})
    args = argparse.Namespace(ld_wind_snps=None, ld_wind_kb=None, ld_wind_cm=1.0)
    coords_a, dist_a = kernel_ldscore.build_window_coordinates(md, args)
    coords_b, dist_b = kernel_ldscore.build_window_coordinates(md.copy(), args)
    np.testing.assert_array_equal(
        kernel_ldscore.get_block_lefts(coords_a, dist_a),
        kernel_ldscore.get_block_lefts(coords_b, dist_b),
    )
    # The kb window over the same physical positions is likewise deterministic.
    args_kb = argparse.Namespace(ld_wind_snps=None, ld_wind_kb=1.0, ld_wind_cm=None)
    coords_kb, dist_kb = kernel_ldscore.build_window_coordinates(md, args_kb)
    np.testing.assert_array_equal(coords_kb, md["POS"].to_numpy(dtype=float))


def test_assert_cm_usable_rejects_all_zero():
    with pytest.raises(LDSCInputError, match="all zero|uninformative"):
        kernel_ldscore.assert_cm_usable(pd.Series([0.0, 0.0, 0.0]), chrom="1")


def test_assert_cm_usable_rejects_constant():
    with pytest.raises(LDSCInputError):
        kernel_ldscore.assert_cm_usable(pd.Series([2.5, 2.5, 2.5]), chrom="1")


def test_assert_cm_usable_accepts_two_distinct_values():
    kernel_ldscore.assert_cm_usable(pd.Series([0.0, 0.1, 0.1]), chrom="1")  # no raise


def test_require_reference_maf_raises_when_all_missing():
    with pytest.raises(LDSCInputError, match="MAF"):
        kernel_ldscore.require_reference_maf(pd.DataFrame({"MAF": [np.nan, np.nan]}), chrom="1")


def test_require_reference_maf_passes_with_any_value():
    kernel_ldscore.require_reference_maf(pd.DataFrame({"MAF": [0.1, np.nan]}), chrom="1")  # no raise


def test_namespace_from_configs_propagates_plink_maf_min():
    """--maf-min reaches the kernel for PLINK (previously hardcoded to None)."""
    spec = SimpleNamespace(
        backend="plink", plink_prefix=None, maf_min=0.4, keep_indivs_file=None, sample_size=None
    )
    ref_panel = SimpleNamespace(spec=spec)
    ns = ldscore_workflow._namespace_from_configs(
        "1", ref_panel, LDScoreConfig(ld_wind_snps=10), GlobalConfig(snp_identifier="rsid")
    )
    assert ns.maf_min == 0.4


def test_build_parser_accepts_genetic_map_and_export_flags():
    parser = ldscore_workflow.build_parser()
    args = parser.parse_args([
        "--output-dir", "out", "--plink-prefix", "panel", "--ld-wind-cm", "1.0",
        "--genetic-map-hg38-sources", "map_hg38.txt", "--export-ref-metadata",
    ])
    assert args.genetic_map_hg38_sources == "map_hg38.txt"
    assert args.export_ref_metadata is True


def test_genetic_map_ignored_for_parquet_with_warning(caplog):
    import logging
    spec = SimpleNamespace(
        backend="parquet_r2", plink_prefix=None, maf_min=None, keep_indivs_file=None,
        sample_size=None, genetic_map_hg19_sources="m19", genetic_map_hg38_sources=None,
    )
    ref_panel = SimpleNamespace(spec=spec)
    with caplog.at_level(logging.WARNING):
        ns = ldscore_workflow._namespace_from_configs(
            "1", ref_panel, LDScoreConfig(ld_wind_cm=1.0),
            GlobalConfig(snp_identifier="chr_pos", genome_build="hg38"),
        )
    assert ns.genetic_map is None  # flags ignored; sidecar CM is authoritative
    assert any("Ignoring --genetic-map" in r.message for r in caplog.records)


def test_validate_window_positions_sorted_accepts_non_decreasing():
    metadata = pd.DataFrame({"CHR": ["1", "1", "1"], "POS": [100, 100, 250]})
    # Non-decreasing (including a tie) must pass without raising.
    kernel_ldscore.validate_window_positions_sorted(metadata, chrom="1")


def test_validate_window_positions_sorted_rejects_unsorted_metadata():
    metadata = pd.DataFrame({"CHR": ["1", "1", "1"], "POS": [100, 300, 250]})
    with pytest.raises(LDSCInternalError, match="non-decreasing order"):
        kernel_ldscore.validate_window_positions_sorted(metadata, chrom="1")


def test_validate_window_positions_sorted_ignores_chromosome_boundary_drop():
    # A position drop across a chromosome boundary is expected, not an error:
    # the window is computed per chromosome.
    metadata = pd.DataFrame({"CHR": ["1", "1", "2"], "POS": [100, 999, 1]})
    kernel_ldscore.validate_window_positions_sorted(metadata, chrom="1")
