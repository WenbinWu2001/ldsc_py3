from __future__ import annotations
from ldsc._runtime.serial_calculation import serial_calculation
from tests.annotation_fixtures import build_annotation_fixture, make_annotation_bundle, fixture_metadata, fixture_values, fixture_ids, replace_fixture_queries

from tests.ref_panel_helpers import compute_plink, write_tiny_plink

from argparse import Namespace
from dataclasses import replace as dataclass_replace
import gzip
import importlib.util
from io import StringIO
import json
from pathlib import Path
from types import SimpleNamespace
import tempfile
import unittest
from unittest import mock
import warnings

import numpy as np
import pandas as pd
import pytest


from ldsc.config import ConfigMismatchError, GlobalConfig
from ldsc.errors import LDSCInputError, LDSCInternalError, LDSCUsageError
from ldsc._annotation.query_outcomes import QueryAnnotationStatus, emit_gene_gate_b_notice

try:
    from ldsc import (
        AnnotationBuildConfig,
        AnnotationBuilder,
        AnnotationBundle,
        LDScoreConfig,
        PlinkRefPanel,
        RefPanelConfig,
    )
    from ldsc._ldscore import calculator as ldscore_workflow
    from ldsc._ldscore.calculator import (
        build_parser,
        _normalize_run_args,
        _resolve_regression_region_build,
    )
    from ldsc._core import plink_files as kernel_formats
    from ldsc._ldscore import compute as kernel_ldscore
    from ldsc._annotation.overlap import compute_overlap
    from ldsc._refpanel import plink as plink_bed
    from ldsc._refpanel import adapters as ref_panel_adapters, r2_format, windows as ld_windows
    from ldsc._refpanel.prepared import ChromosomeAnnotations
    from ldsc._refpanel.r2_format import validate_ldscore_window_within_r2_panel_window
    from ldsc._refpanel.windows import (
        assert_cm_usable,
        validate_window_positions_sorted,
        window_coordinates,
    )
    from ldsc._core.snp_identity import RestrictionIdentityKeys, empty_identity_drop_frame
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
    _normalize_run_args = None
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
    """Write a zero-row pre-index R2 placeholder for tests that never read its pairs."""
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


class R2AutoLoadCLITest(unittest.TestCase):
    @unittest.skipUnless(_HAS_PYARROW, "pyarrow is required for parquet schema coverage")
    def test_package_parquet_metadata_duplicates_remain_invariant_failures(self):
        from ldsc._refpanel.adapters import ParquetR2RefPanel

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
        help_text = " ".join(ldscore_workflow.build_parser().format_help().split())

        self.assertIn("Required for coordinate-family modes", help_text)
        self.assertIn("direct gene-list runs default to auto", help_text)
        self.assertIn("Gene-list runs also use this build for --regr-snps-exclude-regions", help_text)

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
                with self.assertRaisesRegex(LDSCUsageError, "unknown keyword argument"):
                    ldscore_workflow.run_ldscore(output_dir="out", **{name: "value"})

    def test_run_ldscore_requires_output_dir_without_creating_a_placeholder(self):
        with self.assertRaisesRegex(LDSCUsageError, r"missing required keyword argument\(s\) output_dir"):
            ldscore_workflow.run_ldscore(plink_prefix="panel", ld_wind_snps=10)

    def test_calculator_rejects_ref_metadata_export_without_output_config(self):
        with self.assertRaisesRegex(LDSCUsageError, "requires output_config"):
            ldscore_workflow.LDScoreCalculator().run(
                None,
                None,
                LDScoreConfig(ld_wind_snps=10, export_ref_metadata=True),
                GlobalConfig(snp_identifier="rsid"),
            )

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
            regr_snps_file=None,
            log_level="INFO",
        )

        with self.assertRaisesRegex(LDSCUsageError, "coordinate-family SNP identifier mode.*genome build"):
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
                regr_snps_file=None,
                log_level="INFO",
            )

            with self.assertRaisesRegex(LDSCInputError, "Genome-build consistency check failed"):
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
                regr_snps_file=None,
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
                regr_snps_file=None,
                log_level="INFO",
            )

            normalized_args, global_config = ldscore_workflow._normalize_run_args(args)

        self.assertEqual(normalized_args.genome_build, "hg19")
        self.assertEqual(global_config.genome_build, "hg19")

    @unittest.skipUnless(_HAS_PYARROW, "pyarrow is required for parquet schema coverage")
    def test_auto_genome_build_requires_selection_among_r2_alternatives(self):
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
                regr_snps_file=None,
                log_level="INFO",
            )

            with self.assertRaisesRegex(LDSCInputError, "could not be inferred or resolved"):
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
        return make_annotation_bundle(
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

    def test_global_config_rejects_removed_regr_snps_file(self):
        with self.assertRaises(TypeError):
            GlobalConfig(snp_identifier="rsid", regr_snps_file="filters/hm3.tsv.gz")

    def test_ldscore_config_accepts_regr_snps_file(self):
        config = LDScoreConfig(ld_wind_snps=10, regr_snps_file="/path/to/snps.txt")
        self.assertEqual(config.regr_snps_file, "/path/to/snps.txt")

    def test_ldscore_config_rejects_old_regression_snps_file_name(self):
        with self.assertRaises(TypeError):
            LDScoreConfig(ld_wind_snps=10, regression_snps_file="/path/to/snps.txt")

    def test_build_parser_does_not_expose_removed_snp_restriction_flags(self):
        parser = ldscore_workflow.build_parser()
        help_text = parser.format_help()
        for tokens in (
            ("--use-hm3-ref-panel-snps",),
            ("--use-hm3-regression-snps",),
            ("--exclude-regions-bed", "regions.bed"),
            ("--exclude-regions-build", "hg19"),
        ):
            with self.subTest(tokens=tokens):
                self.assertNotIn(tokens[0], help_text)
                with self.assertRaises(SystemExit):
                    parser.parse_args(["--output-dir", "out", "--plink-prefix", "panel", "--ld-wind-snps", "10", *tokens])

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
        from ldsc._annotation.overlap import OverlapContribution

        contribution = OverlapContribution(np.array([[2.0]]), None, np.array([]), None, 2, None)
        legacy_result = kernel_ldscore.ChromComputationResult(
            chrom="1",
            metadata=pd.DataFrame(
                {"CHR": ["1", "1"], "SNP": ["rs1", "rs2"], "POS": [10, 20], "CM": [0.1, 0.2], "MAF": [0.2, 0.3]}
            ),
            ld_scores=np.array([[1.0], [2.0]], dtype=np.float32),
            w_ld=np.array([[3.0], [4.0]], dtype=np.float32),
            M=np.array([2.0]),
            M_5_50=np.array([2.0]),
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
        from ldsc._annotation.overlap import OverlapContribution

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
        from ldsc._annotation.overlap import OverlapContribution

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
        legacy_result = kernel_ldscore.ChromComputationResult(
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
            M_5_50=np.array([2.0]),
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
        return make_annotation_bundle(
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
        bed = plink_bed.PlinkBEDFile(
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
        block_left = ld_windows.get_block_lefts(np.arange(bed.m), 10)
        ld_scores = kernel_ldscore.ld_scores_from_genotypes(bed, block_left, 50, annot=np.ones((bed.m, 1), dtype=np.float32))
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

    def test_finalize_query_statuses_prunes_zero_hit_and_zero_variance_queries(self):
        metadata = pd.DataFrame({"CHR": ["1", "1"], "SNP": ["rs1", "rs2"], "POS": [10, 20]})
        result = ldscore_workflow.LDScoreResult(
            baseline_table=metadata.assign(regression_ld_scores=[1.0, 1.5], base=[1.0, 1.0]),
            query_table=metadata.assign(q_zero=[0.0, 0.0], q_constant=[2.0, 2.0], q_good=[1.0, 3.0]),
            count_records=[
                {"group": "baseline", "column": "base", "all_reference_snp_count": 2.0},
                {"group": "query", "column": "q_zero", "all_reference_snp_count": 0.0},
                {"group": "query", "column": "q_constant", "all_reference_snp_count": 1.0},
                {"group": "query", "column": "q_good", "all_reference_snp_count": 1.0},
            ],
            baseline_columns=["base"],
            query_columns=["q_zero", "q_constant", "q_good"],
            ld_reference_snps=frozenset({"rs1", "rs2"}),
            ld_regression_snps=frozenset({"rs1", "rs2"}),
            chromosome_results=[],
            config_snapshot=GlobalConfig(snp_identifier="rsid"),
        )
        statuses = (
            QueryAnnotationStatus("q_zero", "zero.bed", "bed", "ok"),
            QueryAnnotationStatus("q_constant", "constant.bed", "bed", "ok"),
            QueryAnnotationStatus("q_good", "good.bed", "bed", "ok"),
            QueryAnnotationStatus("q_bad", "bad.bed", "bed", "skipped", "malformed_input"),
        )

        finalized = ldscore_workflow.finalize_query_statuses(result, statuses)

        self.assertEqual(finalized.query_columns, ["q_good"])
        self.assertEqual(finalized.query_table.columns.tolist(), ["CHR", "SNP", "POS", "q_good"])
        self.assertEqual([record["column"] for record in finalized.count_records], ["base", "q_good"])
        self.assertEqual(
            [(row.query, row.status, row.reason, row.n_annotation_snps) for row in finalized.query_statuses],
            [
                ("q_zero", "skipped", "zero_annotation_snps", 0.0),
                ("q_constant", "skipped", "zero_variance_ld_scores", 1.0),
                ("q_good", "ok", "", 1.0),
                ("q_bad", "skipped", "malformed_input", None),
            ],
        )

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
                "--padding-bp",
                "50000",
            ]
        )
        self.assertEqual(args.query_annot_bed_sources, "query.bed")
        self.assertEqual(args.padding_bp, 50000)

    def test_ldscore_rejects_explicit_padding_without_live_interval_queries(self):
        parser = ldscore_workflow.build_parser()
        cases = {
            "no query": ["--output-dir", "out", "--padding-bp", "0"],
            "prebuilt annotation query": [
                "--output-dir", "out",
                "--baseline-annot-sources", "baseline.annot.gz",
                "--query-annot-sources", "query.annot.gz",
                "--padding-bp", "0",
            ],
            "indexed gene query": [
                "--output-dir", "out",
                "--gene-ldscore-index-dir", "index",
                "--query-annot-gene-list-sources", "genes.txt",
                "--padding-bp", "0",
            ],
        }

        for label, argv in cases.items():
            for padding_bp in (0, 100000):
                with self.subTest(label=label, padding_bp=padding_bp):
                    argv[-1] = str(padding_bp)
                    with self.assertRaisesRegex(LDSCUsageError, "Remove.*--padding-bp"):
                        ldscore_workflow.run_ldscore_from_args(parser.parse_args(argv))

    def test_ldscore_accepts_explicit_padding_for_live_interval_queries(self):
        parser = ldscore_workflow.build_parser()
        for query_flag, query_source in (
            ("--query-annot-bed-sources", "query.bed"),
            ("--query-annot-gene-list-sources", "genes.txt"),
        ):
            with self.subTest(query_flag=query_flag):
                args = parser.parse_args(
                    [
                        "--output-dir", "out",
                        "--baseline-annot-sources", "baseline.annot.gz",
                        query_flag, query_source,
                        "--padding-bp", "100000",
                    ]
                )
                ldscore_workflow._validate_padding_usage(args)

    def test_live_gene_lists_require_catalog_and_explicit_padding(self):
        parser = ldscore_workflow.build_parser()
        base = [
            "--output-dir", "out",
            "--baseline-annot-sources", "baseline.annot.gz",
            "--query-annot-gene-list-sources", "genes.txt",
        ]

        with self.assertRaisesRegex(LDSCUsageError, "--padding-bp.*must be chosen"):
            ldscore_workflow.run_ldscore_from_args(parser.parse_args([*base, "--gene-coordinate-file", "genes.tsv"]))
        with self.assertRaisesRegex(LDSCUsageError, "--gene-coordinate-file"):
            ldscore_workflow.run_ldscore_from_args(parser.parse_args([*base, "--padding-bp", "0"]))

    def test_gene_resolution_policy_is_scoped_to_gene_list_modes(self):
        parser = ldscore_workflow.build_parser()
        args = parser.parse_args(
            [
                "--output-dir", "out",
                "--query-annot-bed-sources", "query.bed",
                "--allow-unresolved-genes",
            ]
        )

        with self.assertRaisesRegex(LDSCUsageError, "resolution policy.*gene-list"):
            ldscore_workflow.run_ldscore_from_args(args)

    def test_explicit_default_gene_options_are_rejected_outside_gene_list_mode(self):
        parser = ldscore_workflow.build_parser()
        for option, value in (
            ("--gene-exclude-regions", "none"),
        ):
            with self.subTest(option=option):
                argv = [
                    "--output-dir", "out",
                    "--query-annot-bed-sources", "query.bed",
                    option, value,
                ]
                args = parser.parse_args(argv)
                args._explicit_cli_options = frozenset({"--query-annot-bed-sources", option})
                with self.assertRaisesRegex(LDSCUsageError, "valid only with.*gene-list"):
                    ldscore_workflow.run_ldscore_from_args(args)

    def test_python_run_ldscore_rejects_explicit_padding_without_live_interval_queries(self):
        cases = {
            "no query": {},
            "prebuilt annotation query": {
                "baseline_annot_sources": "baseline.annot.gz",
                "query_annot_sources": "query.annot.gz",
            },
            "indexed gene query": {
                "gene_ldscore_index_dir": "index",
                "query_annot_gene_list_sources": "genes.txt",
            },
        }

        for label, kwargs in cases.items():
            for padding_bp in (0, 100000):
                with self.subTest(label=label, padding_bp=padding_bp):
                    with self.assertRaisesRegex(LDSCUsageError, "Remove.*padding_bp"):
                        ldscore_workflow.run_ldscore(
                            output_dir="out", padding_bp=padding_bp, **kwargs
                        )

    def test_build_parser_accepts_gene_lists_and_keeps_query_routes_mutually_exclusive(self):
        parser = ldscore_workflow.build_parser()
        args = parser.parse_args(
            [
                "--output-dir",
                "out/example",
                "--baseline-annot-sources",
                "baseline.annot.gz",
                "--query-annot-gene-list-sources",
                "immune.txt.gz,brain.list",
                "--gene-coordinate-file",
                "genes.tsv.gz",
                "--padding-bp",
                "0",
                "--allow-unresolved-genes",
            ]
        )

        self.assertEqual(args.query_annot_gene_list_sources, "immune.txt.gz,brain.list")
        self.assertIsNone(args.control_gene_list_file)
        self.assertEqual(args.gene_exclude_regions, "none")
        self.assertEqual(args.gene_coordinate_file, "genes.tsv.gz")
        self.assertEqual(args.gene_list_resolution_policy, "resolved-only")
        with self.assertRaises(SystemExit):
            parser.parse_args(
                [
                    "--output-dir",
                    "out/example",
                    "--query-annot-bed-sources",
                    "query.bed",
                    "--query-annot-gene-list-sources",
                    "immune.txt",
                ]
            )

    def test_annotation_build_config_normalizes_gene_list_sources(self):
        spec = AnnotationBuildConfig(
            query_annot_gene_list_sources="immune.txt,brain.list",
            gene_coordinate_file="genes.tsv",
        )

        self.assertEqual(spec.query_annot_gene_list_sources, ("immune.txt,brain.list",))

    def test_gene_list_auto_build_is_separate_from_rsid_identity_metadata(self):
        args = ldscore_workflow.build_parser().parse_args(
            [
                "--output-dir",
                "out/example",
                "--baseline-annot-sources",
                "baseline.annot.gz",
                "--query-annot-gene-list-sources",
                "immune.txt",
                "--snp-identifier",
                "rsid",
                "--genome-build",
                "auto",
            ]
        )

        with mock.patch.object(
            ldscore_workflow,
            "_resolve_ldscore_chr_pos_genome_build",
            return_value="hg38",
        ) as resolve_build:
            normalized, global_config = _normalize_run_args(args)

        self.assertEqual(normalized.gene_catalog_build, "hg38")
        self.assertIsNone(global_config.genome_build)
        resolve_build.assert_called_once_with(
            normalized,
            "auto",
            gene_catalog_build=None,
        )

    def test_gene_list_auto_build_drives_default_regression_regions_in_rsid_mode(self):
        args = ldscore_workflow.build_parser().parse_args(
            [
                "--output-dir",
                "out/example",
                "--baseline-annot-sources",
                "baseline.annot.gz",
                "--query-annot-gene-list-sources",
                "immune.txt",
                "--snp-identifier",
                "rsid",
            ]
        )

        with mock.patch.object(
            ldscore_workflow,
            "_resolve_ldscore_chr_pos_genome_build",
            return_value="hg38",
        ):
            normalized, global_config = _normalize_run_args(args)

        intervals = ldscore_workflow._regression_region_intervals(normalized, global_config)

        self.assertEqual(
            intervals.source_labels,
            ("preset:mhc[hg38]", "preset:centromeres[hg38]"),
        )

    def test_gene_list_omitted_build_defaults_to_auto_for_chr_pos(self):
        args = ldscore_workflow.build_parser().parse_args(
            [
                "--output-dir",
                "out/example",
                "--baseline-annot-sources",
                "baseline.annot.gz",
                "--query-annot-gene-list-sources",
                "immune.txt",
                "--snp-identifier",
                "chr_pos",
            ]
        )

        with mock.patch.object(
            ldscore_workflow,
            "_resolve_ldscore_chr_pos_genome_build",
            return_value="hg19",
        ) as resolve_build:
            normalized, global_config = _normalize_run_args(args)

        self.assertEqual(normalized.gene_catalog_build, "hg19")
        self.assertEqual(global_config.genome_build, "hg19")
        resolve_build.assert_called_once_with(
            normalized,
            "auto",
            gene_catalog_build=None,
        )

    def test_build_parser_defaults_snp_batch_size_to_512(self):
        parser = ldscore_workflow.build_parser()
        self.assertEqual(parser.get_default("snp_batch_size"), 512)

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

    def test_compute_counts_rejects_metadata_without_complete_maf(self):
        # Both reference adapters always supply MAF, so a gap means a custom or
        # broken adapter; dropping the common counts would silently change the
        # default regression counts.
        annotations = pd.DataFrame({"base": [1.0, 1.0]})
        for label, metadata in {
            "no MAF column": pd.DataFrame({"SNP": ["rs1", "rs2"]}),
            "all missing": pd.DataFrame({"MAF": [np.nan, np.nan]}),
            "one missing": pd.DataFrame({"MAF": [0.2, np.nan]}),
        }.items():
            with self.subTest(label), self.assertRaisesRegex(LDSCInternalError, "MAF"):
                kernel_ldscore.compute_counts(metadata, annotations)

    def test_count_config_records_inclusive_operator(self):
        from ldsc._ldscore.calculator import _count_config_from_ldscore_config
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


    def test_plink_compute_accepts_restriction_identity_keys(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            prefix = self._copy_plink_fixture_with_distinct_fids(Path(tmpdir))
            workflow_bundle = self._build_annotation_bundle(prefix)
            bundle = ChromosomeAnnotations(
                metadata=fixture_metadata(workflow_bundle),
                annotations=fixture_values(workflow_bundle, 'baseline'),
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

            result = compute_plink("1", bundle, args, restriction)

            self.assertEqual(result.chrom, "1")
            self.assertGreater(result.reference_snp_count, 0)
            self.assertEqual(result.metadata.SNP.tolist(), [])
            self.assertEqual(result.w_ld.shape[0], len(result.metadata))

    def test_plink_compute_rejects_unsorted_window_positions(self):
        # geno_meta inherits the annotation bundle's row order via keep_snps, so a
        # bundle whose positions are not coordinate-sorted must trip the window-sort
        # guard before get_block_lefts runs (otherwise LD windows are silently wrong).
        with tempfile.TemporaryDirectory() as tmpdir:
            prefix = self._copy_plink_fixture_with_distinct_fids(Path(tmpdir))
            wb = self._build_annotation_bundle(prefix)
            order = np.argsort(fixture_metadata(wb)["POS"].to_numpy(), kind="mergesort")[::-1]
            bundle = ChromosomeAnnotations(
                metadata=fixture_metadata(wb).iloc[order].reset_index(drop=True),
                annotations=fixture_values(wb, 'baseline').iloc[order].reset_index(drop=True),
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
                compute_plink("1", bundle, args, None)

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
            bundle = ChromosomeAnnotations(
                metadata=fixture_metadata(wb), annotations=fixture_values(wb, 'baseline'),
                baseline_columns=wb.baseline_columns, query_columns=wb.query_columns,
            )
            # The fixture .bim CM is all zero; without a map this must error even under
            # --yes-really (yes_really=True below).
            with self.assertRaisesRegex(LDSCInputError, "uninformative|all zero"):
                compute_plink("1", bundle, self._plink_cm_args(prefix, None), None)

    def test_plink_cm_window_uses_genetic_map(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            prefix = self._copy_plink_fixture_with_distinct_fids(Path(tmpdir))
            wb = self._build_annotation_bundle(prefix)
            bundle = ChromosomeAnnotations(
                metadata=fixture_metadata(wb), annotations=fixture_values(wb, 'baseline'),
                baseline_columns=wb.baseline_columns, query_columns=wb.query_columns,
            )
            gmap = pd.DataFrame({"CHR": ["1", "1"], "POS": [1, 8], "CM": [0.0, 4.0]})
            result = compute_plink("1", bundle, self._plink_cm_args(prefix, gmap), None)
            self.assertGreater(len(result.metadata), 0)
            # interpolated CM is non-degenerate -> a real cM window, not whole-chromosome
            self.assertGreaterEqual(result.metadata["CM"].nunique(), 2)

    def test_compute_chromosome_writes_ref_metadata_sidecar(self):
        from ldsc._refpanel.adapters import RefPanelLoader
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

    @unittest.skipUnless(_HAS_BITARRAY, "bitarray is not installed")
    def test_generated_bool_annotations_produce_identical_ld_scores_to_float32_equivalents(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            prefix = self._copy_plink_fixture_with_distinct_fids(Path(tmpdir))
            base_bundle = self._build_annotation_bundle(prefix)
            generated = np.arange(len(fixture_metadata(base_bundle))) % 2 == 0
            continuous = np.linspace(0.125, 0.875, len(fixture_metadata(base_bundle)), dtype=np.float32)
            bool_bundle = replace_fixture_queries(
                base_bundle,
                query_annotations=pd.DataFrame(
                    {"generated": generated, "continuous": continuous}
                ),
                query_columns=["generated", "continuous"],
            )
            float_bundle = replace_fixture_queries(
                bool_bundle,
                query_annotations=fixture_values(bool_bundle, 'query').astype(np.float32),
            )
            global_config = GlobalConfig(snp_identifier="rsid")
            panel = PlinkRefPanel(
                global_config,
                RefPanelConfig(backend="plink", plink_prefix=str(prefix)),
            )
            ldscore_config = LDScoreConfig(ld_wind_snps=10, whole_chromosome_ok=True)

            bool_result = ldscore_workflow.LDScoreCalculator().compute_chromosome(
                "1", bool_bundle, panel, ldscore_config, global_config
            )
            float_result = ldscore_workflow.LDScoreCalculator().compute_chromosome(
                "1", float_bundle, panel, ldscore_config, global_config
            )

        np.testing.assert_allclose(
            bool_result.baseline_table[["regression_ld_scores", "base"]],
            float_result.baseline_table[["regression_ld_scores", "base"]],
            rtol=0,
            atol=0,
        )
        np.testing.assert_allclose(
            bool_result.query_table[["generated", "continuous"]],
            float_result.query_table[["generated", "continuous"]],
            rtol=0,
            atol=0,
        )
        self.assertEqual(bool_result.count_records, float_result.count_records)

    def test_run_rejects_annotation_bundle_snapshot_mismatch(self):
        calc = ldscore_workflow.LDScoreCalculator()
        annotation_bundle = make_annotation_bundle(
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
                    global_config=GlobalConfig(snp_identifier="rsid"),
                )

        patched_compute.assert_not_called()

    def test_run_ldscore_from_args_writes_outputs(self):
        from tests.test_ldscore_chromosome_scope import write_r2_inputs
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            inputs = write_r2_inputs(root)
            regression = root/'regression.txt'
            regression.write_text('SNP\nrs1\nrs2\nrs3\n')
            argv = ['--output-dir',str(root/'out'),'--r2-dir',inputs.r2_dir,
                    '--baseline-annot-sources',inputs.baseline_annot_sources,
                    '--regr-snps-file',str(regression),'--regr-snps-exclude-regions','none',
                    '--snp-identifier','rsid','--ld-wind-snps','2','--yes-really']
            query = root/'query.annot'
            pd.read_csv(inputs.baseline_annot_sources,sep='\t').drop(columns='base').assign(query=[0,1,0]).to_csv(query,sep='\t',index=False)
            result = ldscore_workflow.run_ldscore_from_args(build_parser().parse_args([*argv,'--query-annot-sources',str(query)]))
            # A parquet --ld-wind-snps run also records its panel/kept SNP counts.
            self.assertEqual(set(result.output_paths),{'baseline','query','metadata','overlap','chromosome_scope',
                                                      'dropped_snps_chr22','ld_window_snp_counts'})
            for path in result.output_paths.values():
                self.assertTrue(Path(path).is_file())
            pd.testing.assert_frame_equal(pd.read_parquet(result.output_paths['query']),result.read_queries(result.query_columns))
            self.assertTrue((root/'out/diagnostics/ldscore.log').is_file())
            self.assertFalse(list((root/'out').glob('.ldsc-annotation-*')))


    def test_run_ldscore_from_args_all_skipped_writes_diagnostics_then_raises(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir) / "ldscore_result"
            args = build_parser().parse_args(
                [
                    "--output-dir",
                    str(output_dir),
                    "--baseline-annot-sources",
                    "baseline.annot",
                    "--query-annot-bed-sources",
                    "bad.bed",
                    "--plink-prefix",
                    "panel",
                    "--ld-wind-snps",
                    "10",
                    "--snp-identifier",
                    "rsid",
                    "--regr-snps-exclude-regions",
                    "none",
                ]
            )
            fixture_root = Path(tmpdir)
            baseline_path = fixture_root / "baseline.annot"
            baseline_path.write_text("CHR\tSNP\tPOS\tbase\n1\trs1\t10\t1\n")
            bed_path = fixture_root / "empty.bed"
            bed_path.write_text("")
            args.baseline_annot_sources = str(baseline_path)
            args.query_annot_bed_sources = str(bed_path)
            args.plink_prefix = str(write_tiny_plink(fixture_root / "panel", pd.DataFrame({"CHR": ["1"], "SNP": ["rs1"], "POS": [10]})))
            annotation_bundle = self.make_annotation_bundle([("1", "rs1", 10)])
            status = QueryAnnotationStatus(
                "bad",
                "bad.bed",
                "bed",
                "skipped",
                "empty_input",
                details="BED source contains no intervals",
            )
            annotation_bundle = dataclass_replace(annotation_bundle, query_statuses=(status,))
            result = ldscore_workflow.LDScoreResult(
                baseline_table=pd.DataFrame(
                    {
                        "CHR": ["1"],
                        "SNP": ["rs1"],
                        "POS": [10],
                        "regression_ld_scores": [1.0],
                        "base": [1.0],
                    }
                ),
                query_table=None,
                count_records=[{"group": "baseline", "column": "base", "all_reference_snp_count": 1.0}],
                baseline_columns=["base"],
                query_columns=[],
                ld_reference_snps=frozenset({"rs1"}),
                ld_regression_snps=frozenset({"rs1"}),
                chromosome_results=[],
                config_snapshot=GlobalConfig(snp_identifier="rsid"),
                query_statuses=(status,),
            )

            with mock.patch(
                "ldsc._annotation.builder.AnnotationBuilder.run",
                autospec=True,
                return_value=annotation_bundle,
            ), mock.patch(
                "ldsc._refpanel.adapters.RefPanelLoader.load",
                autospec=True,
                return_value=self.make_ref_panel_stub(backend="plink"),
            ), mock.patch.object(
                ldscore_workflow.LDScoreCalculator,
                "run",
                autospec=True,
                return_value=result,
            ):
                with self.assertRaisesRegex(LDSCInputError, "all 1 requested query annotations were skipped"):
                    ldscore_workflow.run_ldscore_from_args(args)

            self.assertTrue((output_dir / "diagnostics" / "query_annotation_status.tsv").exists())
            self.assertTrue((output_dir / "diagnostics" / "ldscore.log").exists())
            self.assertFalse((output_dir / "metadata.json").exists())
            self.assertIn("empty_input", (output_dir / "diagnostics" / "ldscore.log").read_text(encoding="utf-8"))

    def test_run_ldscore_from_args_gene_list_writes_query_and_resolution_diagnostics(self):
        # Output ownership is tested here with a deliberately mocked calculator.
        from tests.worker_fixtures import calculation_body, local_task_pool
        self.enterContext(calculation_body())
        self.enterContext(mock.patch.object(ldscore_workflow, "ProcessPoolExecutor", local_task_pool))
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            output_dir = tmpdir / "ldscore_result"
            baseline = tmpdir / "baseline.annot"
            genes = tmpdir / "immune.txt"
            catalog = tmpdir / "genes.tsv"
            baseline.write_text(
                "CHR\tPOS\tSNP\tCM\tbase\n"
                "1\t65419\trs1\t0\t1\n"
                "1\t80000\trs2\t0\t1\n",
                encoding="utf-8",
            )
            genes.write_text("OR4F5\nNOT_A_GENE\n", encoding="utf-8")
            catalog.write_text(
                "gene_id\tgene_name\tchrom\tstart\tend\tgenome_build\n"
                "ENSG00000186092\tOR4F5\t1\t65419\t71585\thg38\n",
                encoding="utf-8",
            )
            args = build_parser().parse_args(
                [
                    "--output-dir",
                    str(output_dir),
                    "--baseline-annot-sources",
                    str(baseline),
                    "--query-annot-gene-list-sources",
                    str(genes),
                    "--gene-coordinate-file",
                    str(catalog),
                    "--allow-unresolved-genes",
                    "--padding-bp",
                    "0",
                    "--plink-prefix",
                    "panel",
                    "--ld-wind-snps",
                    "10",
                    "--snp-identifier",
                    "rsid",
                    "--genome-build",
                    "hg38",
                    "--regr-snps-exclude-regions",
                    "none",
                ]
            )
            args.plink_prefix = str(write_tiny_plink(tmpdir / "panel", pd.read_csv(baseline, sep="\t")))
            args.yes_really = True
            chrom_result = ldscore_workflow.ChromLDScoreResult(
                chrom="1",
                baseline_table=pd.DataFrame(
                    {
                        "CHR": ["1", "1"],
                        "SNP": ["rs1", "rs2"],
                        "POS": [65419, 80000],
                        "regression_ld_scores": [2.0, 2.5],
                        "base": [1.0, 1.5],
                    }
                ),
                query_table=pd.DataFrame(
                    {"CHR": ["1", "1"], "SNP": ["rs1", "rs2"], "POS": [65419, 80000], "immune": [1.0, 2.0]}
                ),
                count_records=[
                    {"group": "baseline", "column": "base", "all_reference_snp_count": 2.0},
                    {"group": "query", "column": "immune", "all_reference_snp_count": 1.0},
                ],
                baseline_columns=["base"],
                query_columns=["immune"],
                ld_reference_snps=frozenset({"rs1", "rs2"}),
                ld_regression_snps=frozenset({"rs1", "rs2"}),
                snp_count_totals={"all_reference_snp_counts": np.array([2.0, 1.0])},
                config_snapshot=GlobalConfig(snp_identifier="rsid"),
            )

            with mock.patch.object(
                ldscore_workflow.LDScoreCalculator,
                "compute_chromosome",
                autospec=True,
                return_value=chrom_result,
            ):
                result = ldscore_workflow.run_ldscore_from_args(args)

            self.assertEqual(result.query_columns, ["immune"])
            self.assertEqual(result.query_statuses[0].status, "warning")
            self.assertTrue((output_dir / "ldscore.query.parquet").exists())
            self.assertEqual(
                pd.read_csv(output_dir / "diagnostics" / "gene_list_audit.tsv.gz", sep="\t")["input_gene"].tolist(),
                ["OR4F5", "NOT_A_GENE"],
            )
            metadata = json.loads((output_dir / "metadata.json").read_text(encoding="utf-8"))
            self.assertEqual(metadata["gene_list_resolution_policy"], "resolved-only")
            self.assertEqual(
                metadata["snp_universe_policy"]["regression_rows_and_weights"]["weight_contributors"],
                "same_filtered_regression_set",
            )
            log_text = (output_dir / "diagnostics" / "ldscore.log").read_text(encoding="utf-8")
            self.assertIn("Gene-list catalog projection build: hg38", log_text)
            self.assertIn("status=warning; reason=partial_gene_resolution", log_text)

    def test_strict_gene_gate_a_batches_focal_and_control_before_reference_load(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            output_dir = tmpdir / "result"
            baseline = tmpdir / "baseline.annot"
            catalog = tmpdir / "genes.tsv"
            focal = tmpdir / "focal.txt"
            control = tmpdir / "control.txt"
            baseline.write_text("CHR\tPOS\tSNP\tCM\tbase\n1\t10\trs1\t0\t1\n", encoding="utf-8")
            catalog.write_text(
                "gene_id\tgene_name\tchrom\tstart\tend\tgenome_build\n"
                "G1\tGENE1\t1\t10\t20\thg19\n",
                encoding="utf-8",
            )
            focal.write_text("GENE1\nMISSING_FOCAL\n", encoding="utf-8")
            control.write_text("MISSING_CONTROL\n", encoding="utf-8")
            from tests.test_gene_ldscore_index import _write_plink_trio
            _write_plink_trio(tmpdir / "panel", ["1"])
            args = build_parser().parse_args(
                [
                    "--output-dir", str(output_dir),
                    "--baseline-annot-sources", str(baseline),
                    "--query-annot-gene-list-sources", str(focal),
                    "--gene-coordinate-file", str(catalog),
                    "--control-gene-list-file", str(control),
                    "--padding-bp", "0",
                    "--plink-prefix", str(tmpdir / "panel"),
                    "--ld-wind-snps", "10",
                    "--snp-identifier", "rsid",
                    "--genome-build", "hg19",
                ]
            )

            with mock.patch(
                "ldsc._refpanel.adapters.RefPanelLoader.load",
                autospec=True,
            ) as load_panel, self.assertRaisesRegex(LDSCInputError, "rejected 2 of 3") as error:
                ldscore_workflow.run_ldscore_from_args(args)

            load_panel.assert_not_called()
            self.assertIn(
                "source == 'focal.txt' & disposition == 'rejected'",
                str(error.exception),
            )
            audit = pd.read_csv(output_dir / "diagnostics" / "gene_list_audit.tsv.gz", sep="\t")
            summary = pd.read_csv(output_dir / "diagnostics" / "gene_list_resolution_summary.tsv", sep="\t")
            self.assertEqual(audit["input_gene"].tolist(), ["GENE1", "MISSING_FOCAL", "MISSING_CONTROL"])
            self.assertEqual(summary["input_role"].tolist(), ["focal", "control"])
            self.assertFalse((output_dir / "diagnostics" / "query_annotation_status.tsv").exists())
            self.assertFalse((output_dir / "metadata.json").exists())
            log_text = (output_dir / "diagnostics" / "ldscore.log").read_text(encoding="utf-8")
            self.assertIn("MISSING_FOCAL", log_text)
            self.assertIn("MISSING_CONTROL", log_text)

    def test_gene_gate_b_batches_zero_support_queries_and_preserves_usable_siblings(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            baseline = tmpdir / "baseline.annot"
            catalog = tmpdir / "genes.tsv"
            mixed = tmpdir / "mixed.txt"
            unsupported = tmpdir / "unsupported.txt"
            control = tmpdir / "control.txt"
            baseline.write_text(
                "CHR\tPOS\tSNP\tCM\tbase\n1\t10\trs1\t0\t1\n1\t100\trs2\t0\t1\n",
                encoding="utf-8",
            )
            catalog.write_text(
                "gene_id\tgene_name\tchrom\tstart\tend\tgenome_build\n"
                "G1\tSUPPORTED\t1\t10\t20\thg19\n"
                "G2\tUNSUPPORTED\t1\t100\t110\thg19\n",
                encoding="utf-8",
            )
            mixed.write_text("G1\nG2\n", encoding="utf-8")
            unsupported.write_text("G2\n", encoding="utf-8")
            control.write_text("G1\n", encoding="utf-8")
            bundle = build_annotation_fixture(AnnotationBuilder(
                GlobalConfig(snp_identifier="rsid"), projection_genome_build="hg19"
            ),
                AnnotationBuildConfig(
                    baseline_annot_sources=(baseline,),
                    query_annot_gene_list_sources=(mixed, unsupported),
                    gene_coordinate_file=catalog,
                    control_gene_list_file=control,
                    gene_list_resolution_policy="resolved-only",
                    padding_bp=0,
                )
            )
            metadata = pd.DataFrame({"CHR": ["1"], "POS": [10], "SNP": ["rs1"]})
            prefix = write_tiny_plink(tmpdir / "panel", metadata)
            panel = PlinkRefPanel(GlobalConfig(snp_identifier="rsid"), RefPanelConfig(backend="plink", plink_prefix=prefix))


            gated, control_error = ldscore_workflow._apply_direct_gene_gate_b(bundle, panel, LDScoreConfig(ld_wind_snps=10, whole_chromosome_ok=True), GlobalConfig(snp_identifier="rsid"))

            self.assertIsNone(control_error)
            self.assertEqual(gated.query_columns, ["mixed"])
            self.assertEqual(
                [(status.query, status.status, status.reason) for status in gated.query_statuses],
                [("mixed", "warning", "partial_snp_support"), ("unsupported", "skipped", "zero_annotation_snps")],
            )
            audit = pd.concat(gated.gene_list_batch.audit_frames(), ignore_index=True)
            self.assertEqual(
                audit.loc[audit["input_gene"].eq("G2"), "disposition"].tolist(),
                ["unsupported", "unsupported"],
            )
            summary = gated.gene_list_batch.summary
            self.assertEqual(summary["zero_support_genes"].tolist(), [1, 1, 0])
            self.assertEqual(summary["genes_with_snp_support"].tolist(), [1, 0, 1])

            stream = StringIO()
            with mock.patch("sys.stderr", stream):
                emit_gene_gate_b_notice(
                    gated.query_statuses,
                    gated.gene_list_batch,
                    command="ldscore",
                )
            notice = stream.getvalue()
            self.assertIn("2 gene-list row(s) have no retained reference-SNP support", notice)
            self.assertIn("1 query annotation(s) retained despite some genes having no retained reference-SNP support (partial_snp_support)", notice)
            self.assertIn("1 query annotation(s) skipped because they overlap no retained reference SNPs (zero_annotation_snps)", notice)
            self.assertIn("Skipped queries are omitted from LD-score outputs.", notice)
            self.assertIn("mixed=partial_snp_support", notice)
            self.assertIn("unsupported=zero_annotation_snps", notice)
            self.assertIn("diagnostics/query_annotation_status.tsv", notice)

    def test_gene_gate_b_uses_the_same_padded_intervals_as_annotation_projection(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            baseline = tmpdir / "baseline.annot"
            catalog = tmpdir / "genes.tsv"
            genes = tmpdir / "padded.txt"
            baseline.write_text(
                "CHR\tPOS\tSNP\tCM\tbase\n1\t10\trs1\t0\t1\n",
                encoding="utf-8",
            )
            catalog.write_text(
                "gene_id\tgene_name\tchrom\tstart\tend\tgenome_build\n"
                "G1\tPADDED\t1\t20\t30\thg19\n",
                encoding="utf-8",
            )
            genes.write_text("G1\n", encoding="utf-8")
            bundle = build_annotation_fixture(AnnotationBuilder(
                GlobalConfig(snp_identifier="rsid"), projection_genome_build="hg19"
            ),
                AnnotationBuildConfig(
                    baseline_annot_sources=(baseline,),
                    query_annot_gene_list_sources=(genes,),
                    gene_coordinate_file=catalog,
                    padding_bp=10,
                )
            )
            metadata = pd.DataFrame({"CHR": ["1"], "POS": [10], "SNP": ["rs1"]})
            prefix = write_tiny_plink(tmpdir / "panel", metadata)
            panel = PlinkRefPanel(GlobalConfig(snp_identifier="rsid"), RefPanelConfig(backend="plink", plink_prefix=prefix))


            gated, control_error = ldscore_workflow._apply_direct_gene_gate_b(bundle, panel, LDScoreConfig(ld_wind_snps=10, whole_chromosome_ok=True), GlobalConfig(snp_identifier="rsid"))

            self.assertIsNone(control_error)
            self.assertEqual(gated.query_columns, ["padded"])
            self.assertEqual(gated.query_statuses[0].status, "ok")
            self.assertEqual(gated.gene_list_batch.summary.loc[0, "genes_with_snp_support"], 1)

    def test_run_ldscore_from_args_synthesizes_base_when_baseline_and_query_are_omitted(self):
        from tests.test_ldscore_chromosome_scope import write_r2_inputs
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            inputs = write_r2_inputs(root)
            regression = root/'regression.txt'
            regression.write_text('SNP\nrs1\nrs2\nrs3\n')
            argv = ['--output-dir',str(root/'out'),'--r2-dir',inputs.r2_dir,
                    '--baseline-annot-sources',inputs.baseline_annot_sources,
                    '--regr-snps-file',str(regression),'--regr-snps-exclude-regions','none',
                    '--snp-identifier','rsid','--ld-wind-snps','2','--yes-really']
            position = argv.index('--baseline-annot-sources')
            del argv[position:position+2]
            result = ldscore_workflow.run_ldscore_from_args(build_parser().parse_args(argv))
            self.assertEqual(result.baseline_columns,['base'])
            self.assertEqual(result.query_columns,[])
            self.assertIsNone(result.read_queries([]))
            self.assertNotIn('query',result.output_paths)
            self.assertEqual(result.count_records[0]['all_reference_snp_count'],3)
            self.assertEqual(result.baseline_table.columns.tolist(),['CHR','SNP','POS','A1','A2','regression_ld_scores','base'])
            np.testing.assert_array_equal(result.baseline_table.base,np.array([1+16384/32767,1+24576/32767,1+8192/32767],dtype=np.float32))


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
                regr_snps_file=None,
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
                "ldsc._refpanel.adapters.RefPanelLoader.load",
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
        # Output ownership is tested here with a deliberately mocked calculator.
        from tests.worker_fixtures import calculation_body, local_task_pool
        self.enterContext(calculation_body())
        self.enterContext(mock.patch.object(ldscore_workflow, "ProcessPoolExecutor", local_task_pool))
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            output_dir = tmpdir / "ldscore_result"
            output_dir.mkdir()
            legacy_log = output_dir / "ldscore.log"
            legacy_log.write_text("legacy log\n", encoding="utf-8")
            from tests.test_gene_ldscore_index import _write_plink_trio
            _write_plink_trio(tmpdir / "panel", ["1"])
            args = Namespace(
                output_dir=str(output_dir),
                query_annot_sources=None,
                query_annot_bed_sources=None,
                baseline_annot_sources=None,
                plink_prefix=str(tmpdir / "panel"),
                snp_identifier="rsid",
                genome_build="hg38",
                ref_panel_snps_file=None,
                regr_snps_file=None,
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
                "ldsc._refpanel.adapters.RefPanelLoader.load",
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
        # Output ownership is tested here with a deliberately mocked calculator.
        from tests.worker_fixtures import calculation_body, local_task_pool
        self.enterContext(calculation_body())
        self.enterContext(mock.patch.object(ldscore_workflow, "ProcessPoolExecutor", local_task_pool))
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            output_dir = tmpdir / "ldscore_result"
            output_dir.mkdir()
            stale = output_dir / "ldscore.query.parquet"
            stale.write_text("stale\n", encoding="utf-8")
            from tests.test_gene_ldscore_index import _write_plink_trio
            _write_plink_trio(tmpdir / "panel", ["1"])
            args = Namespace(
                output_dir=str(output_dir),
                query_annot_sources=None,
                query_annot_bed_sources=None,
                baseline_annot_sources=None,
                plink_prefix=str(tmpdir / "panel"),
                snp_identifier="rsid",
                genome_build="hg38",
                ref_panel_snps_file=None,
                regr_snps_file=None,
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
                "ldsc._refpanel.adapters.RefPanelLoader.load",
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

    def test_run_ldscore_from_args_overwrite_preserves_current_overlap_artifact(self):
        # Output ownership is tested here with a deliberately mocked calculator.
        from tests.worker_fixtures import calculation_body, local_task_pool
        self.enterContext(calculation_body())
        self.enterContext(mock.patch.object(ldscore_workflow, "ProcessPoolExecutor", local_task_pool))
        from ldsc._annotation.overlap import OverlapContribution

        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir) / "ldscore_result"
            argv = [
                "--output-dir",
                str(output_dir),
                "--baseline-annot-sources",
                "baseline.annot.gz",
                "--query-annot-sources",
                "query.annot.gz",
                "--plink-prefix",
                "panel",
                "--snp-identifier",
                "rsid",
                "--regr-snps-exclude-regions",
                "none",
                "--ld-wind-snps",
                "10",
            ]
            fixture_root = Path(tmpdir)
            baseline_path, query_path = fixture_root / "baseline.annot", fixture_root / "query.annot"
            baseline_path.write_text("CHR\tSNP\tPOS\tbase\n1\trs1\t10\t1\n")
            query_path.write_text("CHR\tSNP\tPOS\tquery\n1\trs1\t10\t1\n")
            prefix = write_tiny_plink(fixture_root / "panel", pd.DataFrame({"CHR": ["1"], "SNP": ["rs1"], "POS": [10]}))
            argv[argv.index("--baseline-annot-sources") + 1] = str(baseline_path)
            argv[argv.index("--query-annot-sources") + 1] = str(query_path)
            argv[argv.index("--plink-prefix") + 1] = str(prefix)
            annotation_bundle = replace_fixture_queries(
                self.make_annotation_bundle([("1", "rs1", 10)]),
                query_annotations=pd.DataFrame({"query": [1.0]}),
                query_columns=["query"],
            )
            chrom_result = dataclass_replace(
                self.make_chrom_result("1", 10, 1.0, 5.0),
                overlap=OverlapContribution(
                    baseline_block_all=np.array([[5.0, 4.0]]),
                    baseline_block_common=np.array([[4.0, 3.0]]),
                    query_diagonal_all=np.array([6.0]),
                    query_diagonal_common=np.array([5.0]),
                    n_all=5,
                    n_common=4,
                ),
            )
            ref_panel = self.make_ref_panel_stub(backend="plink")

            with mock.patch(
                "ldsc._annotation.builder.AnnotationBuilder.run",
                autospec=True,
                return_value=annotation_bundle,
            ), mock.patch(
                "ldsc._refpanel.adapters.RefPanelLoader.load",
                autospec=True,
                return_value=ref_panel,
            ), mock.patch.object(
                ldscore_workflow.LDScoreCalculator,
                "compute_chromosome",
                autospec=True,
                return_value=chrom_result,
            ):
                ldscore_workflow.run_ldscore_from_args(build_parser().parse_args(argv))
                result = ldscore_workflow.run_ldscore_from_args(
                    build_parser().parse_args([*argv, "--overwrite"])
                )

            metadata = json.loads((output_dir / "metadata.json").read_text(encoding="utf-8"))
            for relative_path in metadata["files"].values():
                self.assertTrue((output_dir / relative_path).exists(), relative_path)
            self.assertIn("overlap", result.output_paths)
            self.assertEqual(
                [(row["group"], row["column"]) for row in metadata["counts"]],
                [("baseline", "base"), ("query", "query")],
            )
            from ldsc import load_ldscore_from_dir

            loaded = load_ldscore_from_dir(output_dir)
            self.assertIsNotNone(loaded.overlap)
            self.assertEqual(loaded.count_records, result.count_records)

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
            regr_snps_file=None,
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
            regr_snps_file=None,
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
        from tests.test_ldscore_chromosome_scope import write_r2_inputs
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            inputs = write_r2_inputs(root)
            regression = root/'regression.txt'
            regression.write_text('SNP\nrs1\nrs2\nrs3\n')
            argv = ['--output-dir',str(root/'out'),'--r2-dir',inputs.r2_dir,
                    '--baseline-annot-sources',inputs.baseline_annot_sources,
                    '--regr-snps-file',str(regression),'--regr-snps-exclude-regions','none',
                    '--snp-identifier','rsid','--ld-wind-snps','2','--yes-really']
            regression.write_text('SNP\nrs2\n')
            result = ldscore_workflow.run_ldscore_from_args(build_parser().parse_args(argv))
            self.assertEqual(result.baseline_table.SNP.tolist(),['rs2'])
            self.assertEqual(result.ld_regression_snps,frozenset({'rs2'}))
            np.testing.assert_array_equal(result.baseline_table.base,np.array([1+24576/32767],dtype=np.float32))
            np.testing.assert_array_equal(result.baseline_table.regression_ld_scores,[1])
            pd.testing.assert_frame_equal(pd.read_parquet(result.output_paths['baseline']),result.baseline_table)



    def test_run_ldscore_from_args_resolves_annotation_and_plink_suites(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for chrom in range(1,23):
                metadata = pd.DataFrame({'CHR':[str(chrom)],'SNP':[f'rs{chrom}'],'POS':[10]})
                metadata.assign(base=1).to_csv(root/f'baseline.{chrom}.annot.gz',sep='\t',index=False)
                write_tiny_plink(root/f'panel{chrom}',metadata)
            regression = root/'regression.txt'
            regression.write_text('SNP\n'+'\n'.join(f'rs{i}' for i in range(1,23))+'\n')
            args = build_parser().parse_args(['--output-dir',str(root/'out'),
                '--baseline-annot-sources',str(root/'baseline.@.annot.gz'),'--plink-prefix',str(root/'panel@'),
                '--regr-snps-file',str(regression),'--regr-snps-exclude-regions','none',
                '--snp-identifier','rsid','--ld-wind-snps','2','--yes-really'])
            result = ldscore_workflow.run_ldscore_from_args(args)
            self.assertEqual(result.baseline_table.CHR.tolist(),list(map(str,range(1,23))))
            np.testing.assert_allclose(result.baseline_table.base,1,rtol=0,atol=1e-6)


    def test_run_ldscore_from_args_resolves_r2_directory(self):
        from tests.test_ldscore_chromosome_scope import write_r2_inputs
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            inputs = write_r2_inputs(root)
            regression = root/'regression.txt'
            regression.write_text('SNP\nrs1\nrs2\nrs3\n')
            argv = ['--output-dir',str(root/'out'),'--r2-dir',inputs.r2_dir,
                    '--baseline-annot-sources',inputs.baseline_annot_sources,
                    '--regr-snps-file',str(regression),'--regr-snps-exclude-regions','none',
                    '--snp-identifier','rsid','--ld-wind-snps','2','--yes-really']
            result = ldscore_workflow.run_ldscore_from_args(build_parser().parse_args(argv))
            self.assertEqual(result.baseline_table.SNP.tolist(),['rs1','rs2','rs3'])
            np.testing.assert_array_equal(result.baseline_table.base,np.array([1+16384/32767,1+24576/32767,1+8192/32767],dtype=np.float32))


    def test_normalize_run_args_chr_pos_auto_infers_matching_annotation_and_r2_dir_builds(self):
        from ldsc._genome.build_inference import load_packaged_reference_table
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            baseline, r2_dir = root / 'baseline.annot', root / 'panel'
            rows = load_packaged_reference_table().head(500)
            rows[['CHR', 'hg19_POS']].rename(columns={'hg19_POS': 'POS'}).assign(base=1).to_csv(baseline, sep='\t', index=False)
            _write_minimal_r2_parquet(r2_dir / 'chr1_r2.parquet', {b'ldsc:sorted_by_build': b'hg19'})
            args = Namespace(output_dir='out', baseline_annot_sources=str(baseline), query_annot_sources=None,
                             plink_prefix=None, r2_dir=str(r2_dir), snp_identifier='chr_pos', genome_build='auto',
                             keep_indivs_file=None, ref_panel_snps_file=None, regr_snps_file=None, log_level='INFO')
            normalized, config = ldscore_workflow._normalize_run_args(args)
            self.assertEqual(config.genome_build, 'hg19')
            self.assertEqual(normalized.genome_build, 'hg19')

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
                regr_snps_file=None,
                log_level="INFO",
            )

            with self.assertRaisesRegex(LDSCInputError, "Genome-build consistency check failed"):
                ldscore_workflow._normalize_run_args(args)

    def test_normalize_run_args_rsid_concrete_genome_build_warns_and_nulls(self):
        args = Namespace(
            output_dir="out",
            query_annot_sources=None,
            baseline_annot_sources=None,
            plink_prefix="panel",
            snp_identifier="rsid",
            genome_build="hg38",
            keep_indivs_file=None,
            ref_panel_snps_file=None,
            regr_snps_file=None,
            log_level="INFO",
        )

        with self.assertWarnsRegex(UserWarning, "ignored"):
            normalized_args, global_config = ldscore_workflow._normalize_run_args(args)

        self.assertIsNone(normalized_args.genome_build)
        self.assertIsNone(global_config.genome_build)

    @unittest.skipUnless(_HAS_PYARROW, "pyarrow is required for parquet schema coverage")
    def test_normalize_run_args_chr_pos_auto_rejects_annotation_r2_build_mismatch(self):
        from ldsc._genome.build_inference import load_packaged_reference_table
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            baseline, r2_dir = root / 'baseline.annot', root / 'panel'
            rows = load_packaged_reference_table().head(500)
            rows[['CHR', 'hg19_POS']].rename(columns={'hg19_POS': 'POS'}).assign(base=1).to_csv(baseline, sep='\t', index=False)
            _write_minimal_r2_parquet(r2_dir / 'chr1_r2.parquet', {b'ldsc:sorted_by_build': b'hg38'})
            args = Namespace(output_dir='out', baseline_annot_sources=str(baseline), query_annot_sources=None,
                             plink_prefix=None, r2_dir=str(r2_dir), snp_identifier='chr_pos', genome_build='auto',
                             keep_indivs_file=None, ref_panel_snps_file=None, regr_snps_file=None, log_level='INFO')
            with self.assertRaisesRegex(LDSCInputError, 'Conflicting genome-build evidence') as caught:
                ldscore_workflow._normalize_run_args(args)
            self.assertIn(str(baseline), str(caught.exception))
            self.assertIn('chr1_r2.parquet', str(caught.exception))

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
            regr_snps_file=None,
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

        with self.assertRaisesRegex(LDSCUsageError, "keep-indivs-file.*R2 panel"):
            ldscore_workflow.run_ldscore_from_args(args)

    def test_run_ldscore_from_args_warns_and_skips_empty_intersection_chromosome(self):
        from tests.test_ldscore_chromosome_scope import write_r2_inputs
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            inputs = write_r2_inputs(root)
            regression = root/'regression.txt'
            regression.write_text('SNP\nrs1\nrs2\nrs3\n')
            argv = ['--output-dir',str(root/'out'),'--r2-dir',inputs.r2_dir,
                    '--baseline-annot-sources',inputs.baseline_annot_sources,
                    '--regr-snps-file',str(regression),'--regr-snps-exclude-regions','none',
                    '--snp-identifier','rsid','--ld-wind-snps','2','--yes-really']
            from tests.test_reference_preparation import write_panel
            write_panel(Path(inputs.r2_dir),chrom='1')
            metadata = pd.read_csv(inputs.baseline_annot_sources,sep='\t')
            missing = metadata.iloc[[0]].assign(CHR='1',SNP='absent')
            pd.concat([missing,metadata],ignore_index=True).to_csv(inputs.baseline_annot_sources,sep='\t',index=False)
            with self.assertWarnsRegex(UserWarning,'Skipping chromosome 1'):
                result = ldscore_workflow.run_ldscore_from_args(build_parser().parse_args(argv))
            self.assertEqual(result.baseline_table.CHR.tolist(),['22']*3)




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
        bundle = make_annotation_bundle(
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


@unittest.skipIf(kernel_ldscore is None, "ldscore kernel is not available")
class LDScoreParquetNormalizationTest(unittest.TestCase):

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
                reader = r2_format.SortedR2BlockReader(
                    paths=[str(path)],
                    chrom="1",
                    metadata=metadata,
                    identifier_mode="chr_pos",
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
        from ldsc._ldscore import compute as ls
        from ldsc._refpanel.r2_format import _load_full_panel_sidecar
        with tempfile.TemporaryDirectory() as tmp:
            meta, df = self._sidecar(tmp)
            r2 = meta.with_name("chr1_r2.parquet")
            loaded = _load_full_panel_sidecar(str(r2))
            self.assertEqual(list(loaded["CHR"].astype(str)), ["1", "1", "1"])
            self.assertEqual(list(loaded["POS"].astype(int)), [10, 20, 30])

    def test_validate_binding_raises_on_hash_mismatch(self):
        from ldsc._ldscore import compute as ls
        from ldsc._refpanel.r2_format import _validate_index_binding
        from ldsc._core.snp_identity import sidecar_identity_sha256
        with tempfile.TemporaryDirectory() as tmp:
            meta, df = self._sidecar(tmp)
            good = sidecar_identity_sha256(df)
            _validate_index_binding(df, n_snps=3, identity_hash=good, context="t")  # no raise
            with self.assertRaisesRegex(LDSCInputError, "identity hash"):
                _validate_index_binding(df, n_snps=3, identity_hash="0" * 64, context="t")
            with self.assertRaisesRegex(LDSCInputError, "n_snps"):
                _validate_index_binding(df, n_snps=999, identity_hash=good, context="t")


class IndexRemapTest(unittest.TestCase):
    def test_remap_maps_build_rows_to_retained_indices_chr_pos(self):
        from ldsc._refpanel.r2_format import build_index_remap

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
        from ldsc._refpanel.r2_format import build_index_remap

        # Two multi-allelic variants at the same position with distinct,
        # non-strand-complement allele sets: {A,G} vs {A,C}.
        full = pd.DataFrame({"CHR": ["1", "1"], "POS": [10, 10],
                             "SNP": ["a", "b"], "A1": ["A", "A"], "A2": ["G", "C"]})
        retained = full.iloc[[1]].reset_index(drop=True)  # the A/C variant at pos10
        remap, _ = build_index_remap(full, retained, "chr_pos_allele_aware")
        self.assertEqual(remap.tolist(), [-1, 0])


class IndexReaderDecodeTest(unittest.TestCase):
    def _make_panel(self, tmp):
        from ldsc._refpanel import build_kernel as kb
        from ldsc._refpanel import r2_format
        from ldsc._core.snp_identity import sidecar_identity_sha256

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
        r2_format.write_r2_parquet(pair_chunks=dict_chunks(rows), path=r2, genome_build="hg19", n_samples=100,
                            snp_identifier="chr_pos", min_r2=0.0, n_snps=4,
                            sidecar_identity_sha256=sidecar_identity_sha256(full), ld_window_mode="kb", ld_window_value=1000.0)
        return r2, full

    @unittest.skipUnless(_HAS_PYARROW, "pyarrow required")
    def test_streaming_drops_unretained_endpoint(self):
        from ldsc._refpanel.r2_format import SortedR2BlockReader
        with tempfile.TemporaryDirectory() as tmp:
            r2, full = self._make_panel(tmp)
            # retained = panel rows 0,1,2 (drop pos40); matrix order = [10,20,30]
            retained = full.iloc[[0, 1, 2]][["CHR", "POS", "SNP", "A1", "A2", "CM", "MAF"]].reset_index(drop=True)
            reader = SortedR2BlockReader(
                paths=[str(r2)], chrom="1", metadata=retained,
                identifier_mode="chr_pos",
                genome_build="hg19",
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
        from ldsc._refpanel import build_kernel as kb
        from ldsc._refpanel import r2_format
        from ldsc._core.snp_identity import sidecar_identity_sha256

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
            r2_format.write_r2_parquet(
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

            validate_ldscore_window_within_r2_panel_window(
                args,
                parquet_paths=[str(r2)],
                chrom="1",
            )

    @unittest.skipUnless(_HAS_PYARROW, "pyarrow required")
    def test_ldscore_rejects_window_wider_than_recorded_r2_panel_window(self):
        from ldsc._refpanel import build_kernel as kb
        from ldsc._refpanel import r2_format
        from ldsc._core.snp_identity import sidecar_identity_sha256

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
            r2_format.write_r2_parquet(
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
                "user-requested LD window `--ld-wind-kb 1.0` is wider than the input R2 panel window `--ld-wind-kb 0.1`",
            ):
                validate_ldscore_window_within_r2_panel_window(
                    args,
                    parquet_paths=[str(r2)],
                    chrom="1",
                )

    @unittest.skipUnless(_HAS_PYARROW, "pyarrow required")
    def test_ldscore_rejects_any_strictly_larger_same_mode_window(self):
        from ldsc._refpanel import build_kernel as kb
        from ldsc._refpanel import r2_format
        from ldsc._core.snp_identity import sidecar_identity_sha256

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
            r2_format.write_r2_parquet(
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
            with self.assertRaisesRegex(LDSCUsageError, "is wider than the input R2 panel window"):
                validate_ldscore_window_within_r2_panel_window(
                    args,
                    parquet_paths=[str(r2)],
                    chrom="1",
                )

    @unittest.skipUnless(_HAS_PYARROW, "pyarrow required")
    def test_ldscore_rejects_requested_window_mode_differing_from_r2_panel_mode(self):
        from ldsc._refpanel import build_kernel as kb
        from ldsc._refpanel import r2_format
        from ldsc._core.snp_identity import sidecar_identity_sha256

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
            r2_format.write_r2_parquet(
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
            with self.assertRaises(LDSCUsageError) as caught:
                validate_ldscore_window_within_r2_panel_window(
                    args,
                    parquet_paths=[str(r2)],
                    chrom="1",
                )

        message = str(caught.exception)
        self.assertIn("with `--ld-wind-kb 100.0`", message)
        self.assertIn("the panel was built with `--ld-wind-cm 1.0`", message)
        self.assertIn("windows in different units cannot be checked for coverage", message)
        self.assertIn("ldsc build-r2-panel --ld-wind-kb 100.0", message)


class IndexCrossModeParityTest(unittest.TestCase):
    """One index parquet must produce identical LD scores in all four modes."""

    @unittest.skipUnless(_HAS_PYARROW, "pyarrow required")
    def test_ld_scores_identical_across_modes_from_one_index_parquet(self):
        from ldsc._refpanel import build_kernel as kb
        from ldsc._refpanel import r2_format
        from ldsc._core.snp_identity import sidecar_identity_sha256

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
            r2_format.write_r2_parquet(pair_chunks=dict_chunks(pairs), path=r2, genome_build="hg19", n_samples=100,
                                snp_identifier="chr_pos", min_r2=0.0, n_snps=4,
                                sidecar_identity_sha256=sidecar_identity_sha256(panel), ld_window_mode="kb", ld_window_value=1000.0)
            mode_scores = {}
            for mode in ("rsid", "chr_pos", "rsid_allele_aware", "chr_pos_allele_aware"):
                reader = r2_format.SortedR2BlockReader(
                    paths=[str(r2)], chrom="1",
                    metadata=panel[["CHR", "POS", "SNP", "A1", "A2", "CM", "MAF"]].copy(),
                    identifier_mode=mode,
                    genome_build="hg19",
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
        from ldsc._refpanel.r2_format import SortedR2BlockReader
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "chr1_r2.parquet"
            _write_legacy_r2_parquet(path)  # writes old 10-col canonical schema
            meta = pd.DataFrame({"CHR": ["1"], "POS": [10], "SNP": ["a"], "A1": ["A"], "A2": ["G"]})
            with self.assertRaisesRegex(LDSCInputError, "index-format|build-r2-panel"):
                SortedR2BlockReader(paths=[str(path)], chrom="1", metadata=meta,
                                    identifier_mode="chr_pos",
                                    genome_build="hg19")


class R2DequantizationTest(unittest.TestCase):
    @unittest.skipUnless(_HAS_PYARROW, "pyarrow required")
    def test_reader_dequantizes_int16_panel_to_float32(self):
        from ldsc._refpanel import build_kernel as kb
        from ldsc._refpanel import r2_format
        from ldsc._ldscore import compute as kernel_ldscore
        from ldsc._refpanel import r2_format
        from ldsc._core.snp_identity import sidecar_identity_sha256
        panel = pd.DataFrame({"CHR": ["1"] * 3, "POS": [100, 120, 140],
                              "SNP": ["rs1", "rs2", "rs3"], "A1": ["A"] * 3, "A2": ["C"] * 3,
                              "CM": [0.0] * 3, "MAF": [0.3] * 3})
        pairs = [{"i": 0, "j": 1, "R2": 1.0, "sign": "+"},
                 {"i": 0, "j": 2, "R2": 0.2, "sign": "-"}]
        with tempfile.TemporaryDirectory() as tmp:
            meta = _write_index_sidecar(tmp, panel)
            r2 = meta.with_name("chr1_r2.parquet")
            r2_format.write_r2_parquet(pair_chunks=dict_chunks(pairs), path=r2, genome_build="hg19", n_samples=100,
                                snp_identifier="chr_pos", min_r2=0.0, n_snps=3,
                                sidecar_identity_sha256=sidecar_identity_sha256(panel), ld_window_mode="kb", ld_window_value=1000.0)
            reader = r2_format.SortedR2BlockReader(
                paths=[str(r2)], chrom="1", metadata=panel.copy(),
                identifier_mode="chr_pos",
                genome_build="hg19")
            self.assertEqual(reader.contract.r2_scale, 32767.0)
            i, j, r2_values = next(reader.iter_all_pairs())
            self.assertEqual(r2_values.dtype, np.float32)
            # endpoint exact; 0.2 within half-step
            got = dict(zip(zip(i.tolist(), j.tolist()), r2_values.tolist()))
            self.assertEqual(got[(0, 1)], 1.0)
            self.assertAlmostEqual(got[(0, 2)], 0.2, delta=2e-5)

def test_ldscore_parser_regr_snp_region_flag_default_and_removed_alias():
    parser = build_parser()
    # Default excludes MHC + centromeres.
    default_args = parser.parse_args(["--output-dir", "o", "--plink-prefix", "p", "--ld-wind-cm", "1"])
    assert default_args.regr_snps_exclude_regions == "mhc-and-centromeres"
    # Single-choice enum; comma lists are rejected.
    args = parser.parse_args(
        ["--output-dir", "o", "--plink-prefix", "p", "--ld-wind-cm", "1",
         "--regr-snps-exclude-regions", "centromeres", "--genome-build", "hg19"]
    )
    assert args.regr_snps_exclude_regions == "centromeres"
    with pytest.raises(SystemExit):
        parser.parse_args(
            ["--output-dir", "o", "--plink-prefix", "p", "--ld-wind-cm", "1",
             "--exclude-regions", "mhc", "--genome-build", "hg19"]
        )
    assert "--exclude-regions" not in parser.format_help()
    with pytest.raises(SystemExit):
        parser.parse_args(["--output-dir", "o", "--plink-prefix", "p", "--ld-wind-cm", "1",
                           "--regr-snps-exclude-regions", "mhc,centromeres"])


def test_ldscore_region_build_uses_genome_build_and_requires_it_for_rsid():
    # chr_pos: build is taken from the panel build ldscore operates in.
    assert _resolve_regression_region_build(
        Namespace(),
        GlobalConfig(snp_identifier="chr_pos_allele_aware", genome_build="hg19"),
        ("mhc", "centromeres"),
    ) == "hg19"
    # rsid: no genome build is available, so presets require an explicit build.
    with pytest.raises(LDSCUsageError, match="genome-build"):
        _resolve_regression_region_build(
            Namespace(),
            GlobalConfig(snp_identifier="rsid"),
            ("mhc",),
        )
    # The sole build declaration is used even in rsid mode.
    assert _resolve_regression_region_build(
        Namespace(genome_build="hg38"),
        GlobalConfig(snp_identifier="rsid"),
        ("mhc",),
    ) == "hg38"
    # Gene-list projection supplies coordinate evidence even though rsID
    # identity metadata deliberately remains build-independent.
    assert _resolve_regression_region_build(
        Namespace(genome_build=None, gene_catalog_build="hg19"),
        GlobalConfig(snp_identifier="rsid"),
        ("mhc",),
    ) == "hg19"
    # No presets => no build needed.
    assert _resolve_regression_region_build(
        Namespace(),
        GlobalConfig(snp_identifier="rsid"),
        (),
    ) is None


_CHR22 = Path(__file__).resolve().parent / "fixtures" / "minimal_external_resources" / "plink" / "hm3_chr22_subset"


def _chr22_available() -> bool:
    return all(Path(str(_CHR22) + ext).exists() for ext in (".bed", ".bim", ".fam"))


def test_regression_region_mask_leaves_reference_counts_and_ld_contributors_intact():
    from ldsc._core import regions

    metadata = pd.DataFrame(
        {
            "CHR": ["6", "6", "6"],
            "POS": [24_900_000, 25_100_000, 36_000_000],
            "SNP": ["near_mhc", "in_mhc", "outside"],
            "MAF": [0.2, 0.2, 0.2],
        }
    )
    annotations = pd.DataFrame({"base": [1.0, 1.0, 0.0], "query": [0.0, 1.0, 1.0]})
    intervals = regions.load_preset_intervals(("mhc",), "hg19")

    mask = kernel_ldscore.regression_mask_from_keys(
        metadata, {"near_mhc", "in_mhc", "outside"}, "rsid", region_intervals=intervals
    )
    counts, common_counts = kernel_ldscore.compute_counts(metadata, annotations)

    np.testing.assert_array_equal(mask, np.array([1.0, 0.0, 1.0], dtype=np.float32))
    np.testing.assert_array_equal(counts, np.array([2.0, 2.0]))
    np.testing.assert_array_equal(common_counts, np.array([2.0, 2.0]))
    class PairReader:
        def iter_all_pairs(self):
            yield np.array([0]), np.array([1]), np.array([0.5], dtype=np.float32)

    combined = np.c_[annotations.to_numpy(dtype=np.float32), mask]
    scores = kernel_ldscore.ld_score_streaming_from_r2_reader(
        np.zeros(3, dtype=np.int64), combined, PairReader()
    )
    # The retained SNP near MHC still receives the MHC query annotation, while
    # w_ld sees only the identically filtered regression set.
    assert scores[0, 1] == pytest.approx(0.5)
    np.testing.assert_allclose(scores[:, 2], np.array([1.0, 0.5, 1.0], dtype=np.float32))

    legacy = SimpleNamespace(
        chrom="6",
        metadata=metadata.loc[mask.astype(bool)].reset_index(drop=True),
        ld_scores=scores[mask.astype(bool), :2],
        w_ld=scores[mask.astype(bool), 2:],
        M=counts,
        M_5_50=common_counts,
        ldscore_columns=["base", "query"],
        baseline_columns=["base"],
        query_columns=["query"],
        overlap=None,
        reference_snp_count=3, regression_selected_snp_count=3,
        regression_region_removed_snp_count=1,
    )
    wrapped = ldscore_workflow.LDScoreCalculator()._wrap_legacy_chrom_result(
        legacy,
        GlobalConfig(snp_identifier="rsid"),
        {"near_mhc", "in_mhc", "outside"},
        intervals,
    )
    assert wrapped.baseline_table["SNP"].tolist() == ["near_mhc", "outside"]
    np.testing.assert_allclose(wrapped.baseline_table["regression_ld_scores"], [1.0, 1.0])
    assert wrapped.count_records[0]["all_reference_snp_count"] == 2.0
    assert wrapped.regression_region_removed_snp_count == 1


def test_gene_list_mhc_query_keeps_reference_counts_but_not_regression_row(tmp_path):
    from ldsc._core import regions

    baseline = tmp_path / "baseline.annot"
    genes = tmp_path / "hla.txt"
    catalog = tmp_path / "gene-coordinates.tsv"
    baseline.write_text(
        "CHR\tPOS\tSNP\tCM\tbase\n"
        "6\t24900000\tnear_mhc\t0\t1\n"
        "6\t29945000\tin_hla_a\t0\t1\n"
        "6\t36000000\toutside\t0\t1\n",
        encoding="utf-8",
    )
    genes.write_text("HLA-A\n", encoding="utf-8")
    catalog.write_text(
        "gene_id\tgene_name\tchrom\tstart\tend\tgenome_build\n"
        "ENSG_HLA_A\tHLA-A\t6\t29940000\t29950000\thg38\n",
        encoding="utf-8",
    )
    bundle = build_annotation_fixture(AnnotationBuilder(
        GlobalConfig(snp_identifier="rsid"), projection_genome_build="hg38"
    ),
        AnnotationBuildConfig(
            baseline_annot_sources=(baseline,),
            query_annot_gene_list_sources=(genes,),
            gene_coordinate_file=catalog,
        )
    )
    metadata = fixture_metadata(bundle).copy()
    metadata["MAF"] = 0.2
    annotations = fixture_values(bundle, include_query=True)
    intervals = regions.load_preset_intervals(("mhc",), "hg38")
    regression_mask = kernel_ldscore.regression_mask_from_keys(
        metadata,
        {"near_mhc", "in_hla_a", "outside"},
        "rsid",
        region_intervals=intervals,
    )
    counts, common_counts = kernel_ldscore.compute_counts(metadata, annotations)
    overlap = compute_overlap(metadata, annotations, n_baseline=1)

    assert fixture_values(bundle, 'query')["hla"].tolist() == [0.0, 1.0, 0.0]
    np.testing.assert_array_equal(regression_mask, np.array([1.0, 0.0, 1.0], dtype=np.float32))
    np.testing.assert_array_equal(counts, np.array([3.0, 1.0]))
    np.testing.assert_array_equal(common_counts, np.array([3.0, 1.0]))
    np.testing.assert_array_equal(overlap.query_diagonal_all, np.array([1.0]))


def test_ldscore_rejects_non_minor_maf():
    bad = pd.DataFrame({"CHR": ["22"], "SNP": ["rs1"], "POS": [1],
                        "A1": ["A"], "A2": ["G"], "MAF": [0.73]})
    with pytest.raises(LDSCInternalError, match="MAF"):
        ldscore_workflow._assert_canonical_maf(bad)


def test_assert_canonical_maf_passes_for_minor():
    ok = pd.DataFrame({"MAF": [0.0, 0.2, 0.5]})
    ldscore_workflow._assert_canonical_maf(ok)  # no raise


def test_window_coordinates_select_each_window_mode():
    """One window builder serves ldscore and build-r2-panel for every mode."""
    positions = np.array([100, 250, 900, 2_000])
    cm = np.array([0.0, 0.4, 0.9, 1.6])

    coords, dist = window_coordinates(positions, cm, ld_wind_snps=None, ld_wind_kb=None, ld_wind_cm=1.0)
    np.testing.assert_array_equal(coords, cm)
    assert dist == 1.0
    np.testing.assert_array_equal(ld_windows.get_block_lefts(coords, dist), [0, 0, 0, 2])

    coords, dist = window_coordinates(positions, None, ld_wind_snps=None, ld_wind_kb=1.0, ld_wind_cm=None)
    np.testing.assert_array_equal(coords, positions.astype(float))
    assert dist == 1000.0

    coords, dist = window_coordinates(positions, None, ld_wind_snps=2, ld_wind_kb=None, ld_wind_cm=None)
    np.testing.assert_array_equal(coords, [0.0, 1.0, 2.0, 3.0])
    assert dist == 2.0


def test_window_coordinates_count_panel_rows_when_given():
    """``--ld-wind-snps`` counts R2 panel sidecar rows, not retained rows."""
    coords, dist = window_coordinates(
        np.array([10, 20, 30]), None, ld_wind_snps=1, ld_wind_kb=None, ld_wind_cm=None,
        snp_positions=np.array([0, 5, 6]),
    )
    np.testing.assert_array_equal(coords, [0.0, 5.0, 6.0])
    np.testing.assert_array_equal(ld_windows.get_block_lefts(coords, dist), [0, 1, 1])


def test_ldscore_window_rejects_missing_cm():
    metadata = pd.DataFrame({"CHR": ["22"] * 3, "POS": [1, 2, 3], "SNP": ["rs1", "rs2", "rs3"],
                             "CM": [0.0, np.nan, 0.2]})
    config = SimpleNamespace(ld_wind_snps=None, ld_wind_kb=None, ld_wind_cm=1.0, whole_chromosome_ok=True)
    with pytest.raises(LDSCInputError, match="ldscore cannot use `--ld-wind-cm`.*missing CM"):
        ref_panel_adapters._prepare_window(metadata, "22", config, "rsid")


def test_assert_cm_usable_rejects_all_zero():
    with pytest.raises(LDSCInputError, match="all zero|uninformative"):
        assert_cm_usable(pd.Series([0.0, 0.0, 0.0]), chrom="1")


def test_assert_cm_usable_rejects_constant():
    with pytest.raises(LDSCInputError):
        assert_cm_usable(pd.Series([2.5, 2.5, 2.5]), chrom="1")


def test_assert_cm_usable_accepts_two_distinct_values():
    assert_cm_usable(pd.Series([0.0, 0.1, 0.1]), chrom="1")  # no raise


@serial_calculation
def test_public_plink_calculation_prepares_genotypes_once():
    pytest.importorskip("bitarray")
    metadata = pd.read_csv(
        PLINK_FIXTURES / "plink.bim", sep=r"\s+", header=None,
        names=["CHR", "SNP", "CM", "POS", "A1", "A2"],
    )
    bundle = make_annotation_bundle(
        metadata=metadata,
        baseline_annotations=pd.DataFrame({"base": np.ones(len(metadata))}),
        query_annotations=pd.DataFrame(index=metadata.index),
        baseline_columns=["base"], query_columns=[], chromosomes=["1"], source_summary={},
    )
    config = GlobalConfig(snp_identifier="rsid")
    panel = PlinkRefPanel(config, RefPanelConfig(backend="plink", plink_prefix=PLINK_FIXTURES / "plink"))
    readers = []
    reader_class = plink_bed.PlinkBEDFile

    def open_reader(*args, **kwargs):
        reader = reader_class(*args, **kwargs)
        readers.append(reader)
        return reader

    with mock.patch("ldsc._refpanel.adapters.PlinkBEDFile", side_effect=open_reader):
        result = ldscore_workflow.LDScoreCalculator().compute_chromosome(
            "1", bundle, panel, LDScoreConfig(ld_wind_snps=10, whole_chromosome_ok=True), config,
        )

    assert len(readers) == 1
    assert readers[0]._fh.closed
    assert len(result.baseline_table) == 4
    assert result.count_records[0]["all_reference_snp_count"] == 4.0
    assert result.count_records[0]["common_reference_snp_count"] == 4.0


def test_build_parser_accepts_genetic_map_and_export_flags():
    parser = ldscore_workflow.build_parser()
    args = parser.parse_args([
        "--output-dir", "out", "--plink-prefix", "panel", "--ld-wind-cm", "1.0",
        "--genetic-map-hg38-sources", "map_hg38.txt", "--export-ref-metadata",
    ])
    assert args.genetic_map_hg38_sources == "map_hg38.txt"
    assert args.export_ref_metadata is True


def test_validate_window_positions_sorted_accepts_non_decreasing():
    metadata = pd.DataFrame({"CHR": ["1", "1", "1"], "POS": [100, 100, 250]})
    # Non-decreasing (including a tie) must pass without raising.
    validate_window_positions_sorted(metadata, chrom="1")


def test_validate_window_positions_sorted_rejects_unsorted_metadata():
    metadata = pd.DataFrame({"CHR": ["1", "1", "1"], "POS": [100, 300, 250]})
    with pytest.raises(LDSCInternalError, match="non-decreasing order"):
        validate_window_positions_sorted(metadata, chrom="1")


def test_validate_window_positions_sorted_ignores_chromosome_boundary_drop():
    # A position drop across a chromosome boundary is expected, not an error:
    # the window is computed per chromosome.
    metadata = pd.DataFrame({"CHR": ["1", "1", "2"], "POS": [100, 999, 1]})
    validate_window_positions_sorted(metadata, chrom="1")
