import gzip
import importlib.util
import json
import logging
from dataclasses import replace as dataclass_replace
from pathlib import Path
import sys
import tempfile
import unittest
from unittest import mock
from typing import Optional

import numpy as np
import pandas as pd

SRC = Path(__file__).resolve().parents[1] / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from ldsc._kernel import ref_panel_builder as kernel_builder
from ldsc._kernel.snp_identity import IDENTITY_DROP_COLUMNS, RestrictionIdentityKeys, empty_identity_drop_frame
from ldsc.config import GlobalConfig, ReferencePanelBuildConfig
from ldsc import ldscore_calculator, ref_panel_builder, reset_global_config, set_global_config


_HAS_BITARRAY = importlib.util.find_spec("bitarray") is not None
_HAS_PYARROW = importlib.util.find_spec("pyarrow") is not None
_HAS_PYLIFTOVER = importlib.util.find_spec("pyliftover") is not None
MINIMAL_EXTERNAL_FIXTURES = Path(__file__).resolve().parent / "fixtures" / "minimal_external_resources"


def _find_resources_root() -> Optional[Path]:
    for parent in Path(__file__).resolve().parents:
        candidate = parent / "resources"
        if candidate.exists():
            return candidate
    return None


class GeneticMapParserTest(unittest.TestCase):
    def test_load_genetic_map_accepts_alkes_headers_and_normalizes_chr23_to_x(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "map.txt"
            path.write_text(
                "chr position COMBINED_rate(cM/Mb) Genetic_Map(cM)\n"
                "23 10 0.0 0.1\n"
                "23 20 0.0 0.2\n",
                encoding="utf-8",
            )

            genetic_map = kernel_builder.load_genetic_map(path)

            self.assertEqual(genetic_map.columns.tolist(), ["CHR", "POS", "CM"])
            self.assertEqual(genetic_map["CHR"].tolist(), ["X", "X"])
            self.assertEqual(genetic_map["POS"].tolist(), [10, 20])
            self.assertEqual(genetic_map["CM"].tolist(), [0.1, 0.2])

    def test_load_genetic_map_accepts_basic_alias_family(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "map.tsv"
            path.write_text(
                "chromosome bp genetic_map_cm\n"
                "1 100 0.0\n"
                "1 200 1.0\n",
                encoding="utf-8",
            )

            genetic_map = kernel_builder.load_genetic_map(path)

            self.assertEqual(genetic_map["CHR"].tolist(), ["1", "1"])
            self.assertEqual(genetic_map["POS"].tolist(), [100, 200])
            self.assertEqual(genetic_map["CM"].tolist(), [0.0, 1.0])

    def test_load_genetic_map_rejects_numeric_mitochondrial_codes(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "map.tsv"
            path.write_text(
                "chromosome bp genetic_map_cm\n"
                "25 100 0.0\n",
                encoding="utf-8",
            )

            with self.assertRaisesRegex(ValueError, "M' or 'MT"):
                kernel_builder.load_genetic_map(path)

    def test_load_genetic_map_rejects_unsorted_rows(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "map.txt"
            path.write_text(
                "chr position Genetic_Map(cM)\n"
                "1 200 1.0\n"
                "1 100 0.0\n",
                encoding="utf-8",
            )

            with self.assertRaises(ValueError):
                kernel_builder.load_genetic_map(path)

    def test_load_genetic_map_rejects_duplicate_positions(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "map.txt"
            path.write_text(
                "chr position Genetic_Map(cM)\n"
                "1 100 0.0\n"
                "1 100 0.1\n",
                encoding="utf-8",
            )

            with self.assertRaises(ValueError):
                kernel_builder.load_genetic_map(path)


class GeneticMapInterpolationTest(unittest.TestCase):
    def test_interpolate_cm_linearly_and_clamps_endpoints(self):
        genetic_map = pd.DataFrame(
            {
                "CHR": ["1", "1", "1"],
                "POS": [100, 200, 400],
                "CM": [0.0, 1.0, 3.0],
            }
        )

        interpolated = kernel_builder.interpolate_genetic_map_cm(
            chrom="1",
            positions=np.array([50, 150, 300, 500], dtype=np.int64),
            genetic_map=genetic_map,
        )

        np.testing.assert_allclose(interpolated, np.array([0.0, 0.5, 2.0, 3.0]))


class LiftOverTranslatorTest(unittest.TestCase):
    def test_map_positions_filters_unmapped_and_cross_chrom_hits(self):
        class _FakeLiftOver:
            def convert_coordinate(self, chrom, pos0):
                table = {
                    ("chr22", 99): [("chr22", 999, "+", 0)],
                    ("chr22", 199): [("chr21", 1999, "+", 0)],
                    ("chr22", 299): [],
                }
                return table[(chrom, pos0)]

        translator = kernel_builder.LiftOverTranslator.__new__(kernel_builder.LiftOverTranslator)
        translator.source_build = "hg38"
        translator.target_build = "hg19"
        translator.chain_path = None
        translator._identity = False
        translator._liftover = _FakeLiftOver()

        result = translator.map_positions("22", np.array([100, 200, 300], dtype=np.int64))

        np.testing.assert_array_equal(result.translated_positions, np.array([1000], dtype=np.int64))
        np.testing.assert_array_equal(result.keep_mask, np.array([True, False, False]))
        self.assertEqual(result.cross_chrom_count, 1)
        self.assertEqual(result.unmapped_count, 1)

    def test_map_positions_ignores_unsupported_auxiliary_contigs(self):
        class _FakeLiftOver:
            def convert_coordinate(self, chrom, pos0):
                table = {
                    ("chr22", 99): [("chrUn_gl000217", 50, "+", 0), ("chr22", 999, "+", 0)],
                    ("chr22", 199): [("chrUn_gl000217", 60, "+", 0)],
                }
                return table[(chrom, pos0)]

        translator = kernel_builder.LiftOverTranslator.__new__(kernel_builder.LiftOverTranslator)
        translator.source_build = "hg38"
        translator.target_build = "hg19"
        translator.chain_path = None
        translator._identity = False
        translator._liftover = _FakeLiftOver()

        result = translator.map_positions("22", np.array([100, 200], dtype=np.int64))

        np.testing.assert_array_equal(result.translated_positions, np.array([1000], dtype=np.int64))
        np.testing.assert_array_equal(result.keep_mask, np.array([True, False]))
        self.assertEqual(result.cross_chrom_count, 1)
        self.assertEqual(result.unmapped_count, 0)

    def test_map_positions_uses_first_same_chromosome_hit(self):
        class _FakeLiftOver:
            def convert_coordinate(self, chrom, pos0):
                return [("chr21", 499, "+", 0), ("chr22", 999, "+", 0), ("chr22", 1499, "+", 0)]

        translator = kernel_builder.LiftOverTranslator.__new__(kernel_builder.LiftOverTranslator)
        translator.source_build = "hg38"
        translator.target_build = "hg19"
        translator.chain_path = None
        translator._identity = False
        translator._liftover = _FakeLiftOver()

        result = translator.map_positions("22", np.array([100], dtype=np.int64))

        np.testing.assert_array_equal(result.translated_positions, np.array([1000], dtype=np.int64))
        np.testing.assert_array_equal(result.keep_mask, np.array([True]))
        self.assertEqual(result.cross_chrom_count, 0)
        self.assertEqual(result.unmapped_count, 0)

    def test_explicit_chain_path_is_required_for_cross_build_translation(self):
        with self.assertRaises(ValueError) as exc:
            kernel_builder.LiftOverTranslator(
                source_build="hg38",
                target_build="hg19",
                chain_path=None,
                chain_flag_hint="--liftover-chain-file",
                workflow_label="summary-statistics liftover",
            )

        self.assertIn("summary-statistics liftover", str(exc.exception))
        self.assertIn("--liftover-chain-file", str(exc.exception))


class RestrictionModeDetectionTest(unittest.TestCase):
    def test_detect_restriction_identifier_mode_accepts_rsid_header(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "restrict.txt"
            path.write_text("rsid\nrs1\nrs2\n", encoding="utf-8")

            mode = kernel_builder.detect_restriction_identifier_mode(path)

            self.assertEqual(mode, "rsid")

    def test_detect_restriction_identifier_mode_accepts_chr_pos_header(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "restrict.txt"
            path.write_text("CHR\tPOS\n1\t100\n1\t200\n", encoding="utf-8")

            mode = kernel_builder.detect_restriction_identifier_mode(path)

            self.assertEqual(mode, "chr_pos")

    def test_detect_restriction_identifier_mode_rejects_headerless_input(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "restrict.txt"
            path.write_text("rs1\nrs2\n", encoding="utf-8")

            with self.assertRaisesRegex(ValueError, "header"):
                kernel_builder.detect_restriction_identifier_mode(path)


class _SequentialSNPGetter:
    def __init__(self, matrix: np.ndarray):
        self.matrix = matrix
        self.cursor = 0

    def __call__(self, width: int) -> np.ndarray:
        block = self.matrix[:, self.cursor:self.cursor + width]
        self.cursor += width
        return block


class PairwiseEmissionTest(unittest.TestCase):
    def test_iter_pairwise_r2_rows_emits_one_unordered_pair_per_window(self):
        standardized = np.array(
            [
                [-1.0, -1.0, 1.0, 1.0],
                [-1.0, 1.0, -1.0, 1.0],
                [1.0, -1.0, -1.0, 1.0],
                [1.0, 1.0, 1.0, -1.0],
            ]
        )
        block_left = np.array([0, 0, 1, 2], dtype=float)

        rows = kernel_builder.iter_pairwise_r2_rows(
            block_left=block_left,
            snp_batch_size=2,
            standardized_snp_getter=_SequentialSNPGetter(standardized),
            m=4,
            n=4,
        )

        expected_pairs = {(0, 1), (1, 2), (2, 3)}
        observed_pairs = {(row["i"], row["j"]) for row in rows}
        self.assertEqual(observed_pairs, expected_pairs)
        self.assertTrue(all(row["i"] < row["j"] for row in rows))

        brute_force = {}
        for j in range(1, standardized.shape[1]):
            for i in range(int(block_left[j]), j):
                corr = float(np.dot(standardized[:, i], standardized[:, j]) / standardized.shape[0])
                r2 = corr**2 - (1.0 - corr**2) / (standardized.shape[0] - 2)
                brute_force[(i, j)] = (r2, "+" if corr >= 0 else "-")

        for row in rows:
            expected_r2, expected_sign = brute_force[(row["i"], row["j"])]
            self.assertAlmostEqual(row["R2"], expected_r2)
            self.assertEqual(row["sign"], expected_sign)

    def test_yield_pairwise_r2_rows_emits_non_decreasing_left_index(self):
        standardized = np.array(
            [
                [-1.0, -1.0, 1.0, 1.0, -1.0],
                [-1.0, 1.0, -1.0, 1.0, 1.0],
                [1.0, -1.0, -1.0, 1.0, -1.0],
                [1.0, 1.0, 1.0, -1.0, 1.0],
            ]
        )
        block_left = np.array([0, 0, 0, 1, 1], dtype=int)

        rows = kernel_builder.iter_pairwise_r2_rows(
            block_left=block_left,
            snp_batch_size=2,
            standardized_snp_getter=_SequentialSNPGetter(standardized),
            m=5,
            n=4,
        )

        left_indices = [row["i"] for row in rows]
        self.assertEqual(left_indices, sorted(left_indices))


class RetainedSnpOrderingTest(unittest.TestCase):
    def test_sort_retained_snps_uses_requested_build_positions(self):
        keep_snps = np.array([10, 20, 30], dtype=int)
        hg19_lookup = {10: 2000, 20: 1000, 30: 3000}
        hg38_lookup = {10: 10685988, 20: 10685981, 30: 10686000}

        sorted_snps = ref_panel_builder._sort_retained_snps_by_build_position(
            keep_snps,
            genome_build="hg38",
            hg19_lookup=hg19_lookup,
            hg38_lookup=hg38_lookup,
        )

        self.assertEqual(sorted_snps.tolist(), [20, 10, 30])

    def test_kb_window_coordinates_use_requested_build_positions(self):
        metadata = pd.DataFrame(
            {
                "CHR": ["1", "1"],
                "POS": [1_000, 2_000],
                "SNP": ["rs1", "rs2"],
            }
        )
        build_positions = np.array([100, 250], dtype=np.int64)

        coords, max_dist = kernel_builder.build_window_coordinates(
            metadata=ref_panel_builder._metadata_with_build_positions(metadata, build_positions),
            cm_values=None,
            ld_wind_snps=None,
            ld_wind_kb=1.0,
            ld_wind_cm=None,
        )

        np.testing.assert_array_equal(coords, np.array([100.0, 250.0]))
        self.assertEqual(max_dist, 1000.0)


class StandardTableFormattingTest(unittest.TestCase):
    def test_build_reference_snp_table_uses_exact_schema(self):
        metadata = pd.DataFrame(
            {
                "CHR": ["1", "X"],
                "SNP": ["1:100:A:G", "rsX"],
                "POS": [100, 200],
                "MAF": [0.2, 0.3],
                "A1": ["A", "C"],
                "A2": ["G", "T"],
            }
        )

        table = kernel_builder.build_reference_snp_table(
            metadata=metadata,
            hg19_positions=np.array([100, 150], dtype=np.int64),
            hg38_positions=np.array([110, 250], dtype=np.int64),
        )

        self.assertEqual(
            table.columns.tolist(),
            ["chr", "hg19_pos", "hg38_pos", "hg19_Uniq_ID", "hg38_Uniq_ID", "rsID", "MAF", "A1", "A2"],
        )
        self.assertEqual(table["hg19_Uniq_ID"].tolist(), ["1:100:A:G", "X:150:C:T"])
        self.assertEqual(table["hg38_Uniq_ID"].tolist(), ["1:110:A:G", "X:250:C:T"])
        self.assertEqual(table["rsID"].tolist(), ["1:100:A:G", "rsX"])

    def test_build_reference_snp_table_allows_missing_opposite_build_columns(self):
        metadata = pd.DataFrame(
            {
                "CHR": ["1"],
                "SNP": ["rs1"],
                "POS": [100],
                "MAF": [0.2],
                "A1": ["A"],
                "A2": ["G"],
            }
        )

        table = kernel_builder.build_reference_snp_table(
            metadata=metadata,
            hg19_positions=np.array([100], dtype=np.int64),
            hg38_positions=None,
        )

        self.assertEqual(
            table.columns.tolist(),
            ["chr", "hg19_pos", "hg38_pos", "hg19_Uniq_ID", "hg38_Uniq_ID", "rsID", "MAF", "A1", "A2"],
        )
        self.assertEqual(table.loc[0, "hg19_pos"], 100)
        self.assertEqual(table.loc[0, "hg19_Uniq_ID"], "1:100:A:G")
        self.assertTrue(pd.isna(table.loc[0, "hg38_pos"]))
        self.assertTrue(pd.isna(table.loc[0, "hg38_Uniq_ID"]))

    def test_build_standard_r2_table_uses_exact_schema(self):
        reference_snp_table = pd.DataFrame(
            {
                "CHR": ["1", "1"],
                "hg19_pos": [100, 120],
                "hg38_pos": [110, 130],
                "hg19_Uniq_ID": ["1:100:A:G", "1:120:C:T"],
                "hg38_Uniq_ID": ["1:110:A:G", "1:130:C:T"],
                "rsID": ["rs1", "rs2"],
                "MAF": [0.2, 0.3],
                "A1": ["A", "C"],
                "A2": ["G", "T"],
            }
        )
        pair_rows = [{"i": 0, "j": 1, "R2": 0.75, "sign": "-"}]

        table = kernel_builder.build_standard_r2_table(
            pair_rows=pair_rows,
            reference_snp_table=reference_snp_table,
            genome_build="hg19",
        )

        self.assertEqual(
            table.columns.tolist(),
            ["CHR", "POS_1", "POS_2", "SNP_1", "SNP_2", "A1_1", "A2_1", "A1_2", "A2_2", "R2"],
        )
        self.assertEqual(table.loc[0, "CHR"], "1")
        self.assertEqual(table.loc[0, "POS_1"], 100)
        self.assertEqual(table.loc[0, "POS_2"], 120)
        self.assertAlmostEqual(table.loc[0, "R2"], 0.75, places=4)
        self.assertEqual(table.loc[0, "SNP_1"], "rs1")
        self.assertEqual(table.loc[0, "SNP_2"], "rs2")
        self.assertEqual(table.loc[0, "A1_1"], "A")
        self.assertEqual(table.loc[0, "A2_1"], "G")
        self.assertEqual(table.loc[0, "A1_2"], "C")
        self.assertEqual(table.loc[0, "A2_2"], "T")
        self.assertEqual(table["POS_1"].dtype, np.dtype("int64"))
        self.assertEqual(table["POS_2"].dtype, np.dtype("int64"))
        self.assertEqual(table["R2"].dtype, np.dtype("float32"))
        self.assertTrue(pd.api.types.is_string_dtype(table["CHR"]))
        self.assertTrue(pd.api.types.is_string_dtype(table["SNP_1"]))
        self.assertTrue(pd.api.types.is_string_dtype(table["SNP_2"]))

    def test_build_standard_r2_table_preserves_empty_schema_dtypes(self):
        reference_snp_table = pd.DataFrame(
            {
                "CHR": ["1"],
                "hg19_pos": [100],
                "hg38_pos": [110],
                "hg19_Uniq_ID": ["1:100:A:G"],
                "hg38_Uniq_ID": ["1:110:A:G"],
                "rsID": ["rs1"],
                "MAF": [0.2],
                "A1": ["A"],
                "A2": ["G"],
            }
        )

        table = kernel_builder.build_standard_r2_table(
            pair_rows=[],
            reference_snp_table=reference_snp_table,
            genome_build="hg19",
        )

        self.assertEqual(
            table.columns.tolist(),
            ["CHR", "POS_1", "POS_2", "SNP_1", "SNP_2", "A1_1", "A2_1", "A1_2", "A2_2", "R2"],
        )
        self.assertEqual(table["POS_1"].dtype, np.dtype("int64"))
        self.assertEqual(table["POS_2"].dtype, np.dtype("int64"))
        self.assertEqual(table["R2"].dtype, np.dtype("float32"))
        self.assertTrue(pd.api.types.is_string_dtype(table["CHR"]))
        self.assertTrue(pd.api.types.is_string_dtype(table["SNP_1"]))
        self.assertTrue(pd.api.types.is_string_dtype(table["SNP_2"]))

    @unittest.skipIf(_HAS_PYARROW, "pyarrow dependency is installed")
    def test_write_r2_parquet_requires_pyarrow_for_canonical_output(self):
        reference_snp_table = pd.DataFrame(
            {
                "CHR": ["1", "1"],
                "hg19_pos": [100, 120],
                "hg38_pos": [110, 130],
                "hg19_Uniq_ID": ["1:100:A:G", "1:120:C:T"],
                "hg38_Uniq_ID": ["1:110:A:G", "1:130:C:T"],
                "rsID": ["rs1", "rs2"],
                "MAF": [0.2, 0.3],
                "A1": ["A", "C"],
                "A2": ["G", "T"],
            }
        )
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "chr1.parquet"
            with self.assertRaises(ImportError) as ctx:
                kernel_builder.write_r2_parquet(
                    pair_rows=iter([{"i": 0, "j": 1, "R2": 0.75, "sign": "+"}]),
                    reference_snp_table=reference_snp_table,
                    path=path,
                    genome_build="hg19",
                    n_samples=10,
                    snp_identifier="chr_pos",
                )
            self.assertIn("requires pyarrow", str(ctx.exception))
            self.assertNotIn("fastparquet", str(ctx.exception))

    @unittest.skipUnless(_HAS_PYARROW, "pyarrow dependency is not installed")
    def test_write_r2_parquet_asserts_sort_invariant(self):
        reference_snp_table = pd.DataFrame(
            {
                "CHR": ["1", "1", "1"],
                "hg19_pos": [100, 80, 120],
                "hg38_pos": [110, 90, 130],
                "hg19_Uniq_ID": ["1:100:A:G", "1:80:C:T", "1:120:G:A"],
                "hg38_Uniq_ID": ["1:110:A:G", "1:90:C:T", "1:130:G:A"],
                "rsID": ["rs1", "rs2", "rs3"],
                "MAF": [0.2, 0.3, 0.1],
                "A1": ["A", "C", "G"],
                "A2": ["G", "T", "A"],
            }
        )
        pair_rows = [
            {"i": 0, "j": 2, "R2": 0.5, "sign": "+"},
            {"i": 1, "j": 2, "R2": 0.3, "sign": "+"},
        ]
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "chr1.parquet"
            with self.assertRaises(ValueError) as ctx:
                kernel_builder.write_r2_parquet(
                    pair_rows=iter(pair_rows),
                    reference_snp_table=reference_snp_table,
                    path=path,
                    genome_build="hg19",
                    n_samples=10,
                    snp_identifier="chr_pos",
                )
            self.assertIn("POS_1=80", str(ctx.exception))
            self.assertIn("POS_1=100", str(ctx.exception))

    @unittest.skipUnless(_HAS_PYARROW, "pyarrow dependency is not installed")
    def test_write_r2_parquet_writes_schema_metadata(self):
        import pyarrow.parquet as pq

        reference_snp_table = pd.DataFrame(
            {
                "CHR": ["1", "1"],
                "hg19_pos": [100, 120],
                "hg38_pos": [110, 130],
                "hg19_Uniq_ID": ["1:100:A:G", "1:120:C:T"],
                "hg38_Uniq_ID": ["1:110:A:G", "1:130:C:T"],
                "rsID": ["rs1", "rs2"],
                "MAF": [0.2, 0.3],
                "A1": ["A", "C"],
                "A2": ["G", "T"],
            }
        )
        pair_rows = [{"i": 0, "j": 1, "R2": 0.75, "sign": "+"}]
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "chr1.parquet"
            kernel_builder.write_r2_parquet(
                pair_rows=iter(pair_rows),
                reference_snp_table=reference_snp_table,
                path=path,
                genome_build="hg19",
                n_samples=10,
                snp_identifier="rsid_allele_aware",
                row_group_size=50_000,
            )
            pf = pq.ParquetFile(str(path))
            meta = pf.schema_arrow.metadata
            self.assertIsNotNone(meta)
            self.assertIn(b"ldsc:sorted_by_build", meta)
            self.assertEqual(meta[b"ldsc:sorted_by_build"].decode("utf-8"), "hg19")
            self.assertIn(b"ldsc:row_group_size", meta)
            self.assertEqual(meta[b"ldsc:row_group_size"].decode("utf-8"), "50000")
            self.assertEqual(meta[b"ldsc:schema_version"].decode("utf-8"), "1")
            self.assertEqual(meta[b"ldsc:artifact_type"].decode("utf-8"), "ref_panel_r2")
            self.assertEqual(meta[b"ldsc:snp_identifier"].decode("utf-8"), "rsid_allele_aware")
            self.assertEqual(meta[b"ldsc:genome_build"].decode("utf-8"), "hg19")

    @unittest.skipUnless(_HAS_PYARROW, "pyarrow dependency is not installed")
    def test_write_r2_parquet_stores_n_samples_and_r2_bias(self):
        import pyarrow.parquet as pq

        ref_table = kernel_builder.build_reference_snp_table(
            metadata=pd.DataFrame(
                {
                    "CHR": ["1", "1"],
                    "SNP": ["rs1", "rs2"],
                    "A1": ["A", "T"],
                    "A2": ["C", "G"],
                    "MAF": [0.3, 0.4],
                }
            ),
            hg19_positions=np.array([100, 200], dtype=int),
            hg38_positions=None,
        )
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "hg19" / "chr1_r2.parquet"
            kernel_builder.write_r2_parquet(
                pair_rows=iter([]),
                reference_snp_table=ref_table,
                path=path,
                genome_build="hg19",
                n_samples=42,
                snp_identifier="chr_pos",
            )
            meta = pq.read_schema(str(path)).metadata

        self.assertEqual(meta[b"ldsc:n_samples"], b"42")
        self.assertEqual(meta[b"ldsc:r2_bias"], b"unbiased")

    def test_build_runtime_metadata_table_is_build_specific(self):
        metadata = pd.DataFrame(
            {
                "CHR": ["1"],
                "SNP": ["rs1"],
                "A1": ["A"],
                "A2": ["C"],
                "MAF": [0.2],
            }
        )

        table = kernel_builder.build_runtime_metadata_table(
            metadata=metadata,
            positions=np.array([120], dtype=np.int64),
            cm_values=np.array([0.5], dtype=float),
        )

        self.assertEqual(table.columns.tolist(), ["CHR", "POS", "SNP", "A1", "A2", "CM", "MAF"])
        self.assertEqual(table.loc[0, "POS"], 120)
        self.assertEqual(table.loc[0, "CM"], 0.5)

    def test_build_runtime_metadata_table_preserves_alleles(self):
        metadata = pd.DataFrame(
            {
                "CHR": ["1"],
                "SNP": ["rs1"],
                "MAF": [0.2],
                "A1": ["A"],
                "A2": ["C"],
            }
        )

        table = kernel_builder.build_runtime_metadata_table(
            metadata=metadata,
            positions=np.array([120], dtype=np.int64),
            cm_values=np.array([0.5], dtype=float),
        )

        self.assertEqual(table.columns.tolist(), ["CHR", "POS", "SNP", "A1", "A2", "CM", "MAF"])
        self.assertEqual(table.loc[0, "A1"], "A")
        self.assertEqual(table.loc[0, "A2"], "C")

    def test_build_runtime_metadata_table_allows_missing_cm_values(self):
        metadata = pd.DataFrame(
            {
                "CHR": ["1"],
                "SNP": ["rs1"],
                "A1": ["A"],
                "A2": ["C"],
                "MAF": [0.2],
            }
        )

        table = kernel_builder.build_runtime_metadata_table(
            metadata=metadata,
            positions=np.array([120], dtype=np.int64),
            cm_values=None,
        )

        self.assertEqual(table.columns.tolist(), ["CHR", "POS", "SNP", "A1", "A2", "CM", "MAF"])
        self.assertEqual(table.loc[0, "POS"], 120)
        self.assertTrue(pd.isna(table.loc[0, "CM"]))


class ReferencePanelBuildConfigOptionalLiftoverTest(unittest.TestCase):
    def test_python_api_defaults_missing_maps_for_snp_window(self):
        config = ReferencePanelBuildConfig(
            plink_prefix="plink/panel.@",
            source_genome_build="hg38",
            output_dir="out",
            ld_wind_snps=10,
        )

        self.assertIsNone(config.genetic_map_hg19_sources)
        self.assertIsNone(config.genetic_map_hg38_sources)
        self.assertEqual(config.output_dir, "out")

    def test_accepts_missing_chain_and_maps_for_kb_window(self):
        config = ReferencePanelBuildConfig(
            plink_prefix="plink/panel.@",
            source_genome_build="hg19",
            genetic_map_hg19_sources=None,
            genetic_map_hg38_sources=None,
            output_dir="out",
            ld_wind_kb=1.0,
        )

        self.assertIsNone(config.liftover_chain_hg19_to_hg38_file)
        self.assertIsNone(config.genetic_map_hg19_sources)
        self.assertIsNone(config.genetic_map_hg38_sources)

    def test_missing_source_map_raises_for_cm_window(self):
        with self.assertRaisesRegex(ValueError, "hg19 genetic map.*required.*ld_wind_cm"):
            ReferencePanelBuildConfig(
                plink_prefix="plink/panel.@",
                source_genome_build="hg19",
                genetic_map_hg19_sources=None,
                genetic_map_hg38_sources="maps/hg38.map",
                output_dir="out",
                ld_wind_cm=1.0,
            )

    def test_matching_chain_does_not_require_target_map_for_kb_window(self):
        config = ReferencePanelBuildConfig(
            plink_prefix="plink/panel.@",
            source_genome_build="hg19",
            genetic_map_hg19_sources="maps/hg19.map",
            genetic_map_hg38_sources=None,
            liftover_chain_hg19_to_hg38_file="chains/hg19ToHg38.over.chain",
            output_dir="out",
            ld_wind_kb=1.0,
        )

        self.assertEqual(config.liftover_chain_hg19_to_hg38_file, "chains/hg19ToHg38.over.chain")
        self.assertIsNone(config.genetic_map_hg38_sources)

    def test_hm3_quick_liftover_emits_source_and_target_builds(self):
        config = ReferencePanelBuildConfig(
            plink_prefix="plink/panel.@",
            source_genome_build="hg19",
            output_dir="out",
            ld_wind_kb=1.0,
            use_hm3_snps=True,
            use_hm3_quick_liftover=True,
        )

        self.assertEqual(ref_panel_builder._emitted_genome_builds(config), ["hg19", "hg38"])

    def test_hm3_quick_liftover_drops_unmapped_rows_with_sidecar_frame(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            hm3_path = tmpdir / "hm3.tsv"
            hm3_path.write_text("CHR\thg19_POS\thg38_POS\tSNP\n1\t100\t1000\trs1\n", encoding="utf-8")
            chrom_df = pd.DataFrame(
                {
                    "CHR": ["1", "1"],
                    "SNP": ["rs1", "rs_unmapped"],
                    "BP": [100, 200],
                }
            )
            state = ref_panel_builder._BuildState(
                genetic_map_hg19=None,
                genetic_map_hg38=None,
                use_hm3_quick_liftover=True,
                hm3_map_file=str(hm3_path),
            )

            keep, hg19_lookup, hg38_lookup, drop_frame = ref_panel_builder.ReferencePanelBuilder(
                GlobalConfig(snp_identifier="chr_pos", genome_build="hg19")
            )._resolve_mappable_snp_positions(
                build_state=state,
                chrom="1",
                source_build="hg19",
                chrom_df=chrom_df,
                keep_snps=np.asarray([0, 1], dtype=int),
                sidecar_path=tmpdir / "dropped.tsv.gz",
            )

        self.assertEqual(keep.tolist(), [0])
        self.assertEqual(hg19_lookup, {0: 100})
        self.assertEqual(hg38_lookup, {0: 1000})
        self.assertEqual(drop_frame["SNP"].tolist(), ["rs_unmapped"])
        self.assertEqual(drop_frame["reason"].tolist(), ["unmapped_liftover"])

    def test_non_matching_chain_does_not_require_target_map(self):
        config = ReferencePanelBuildConfig(
            plink_prefix="plink/panel.@",
            source_genome_build="hg19",
            genetic_map_hg19_sources="maps/hg19.map",
            genetic_map_hg38_sources=None,
            liftover_chain_hg38_to_hg19_file="chains/hg38ToHg19.over.chain",
            output_dir="out",
            ld_wind_kb=1.0,
        )

        self.assertEqual(config.liftover_chain_hg38_to_hg19_file, "chains/hg38ToHg19.over.chain")


class ReferencePanelBuildConfigNoDuplicatePolicyTest(unittest.TestCase):
    def _base_config(self, **kwargs):
        return ReferencePanelBuildConfig(
            plink_prefix="plink/panel.@",
            source_genome_build="hg19",
            output_dir="out",
            ld_wind_kb=1.0,
            **kwargs,
        )

    def test_duplicate_position_policy_field_is_removed(self):
        config = self._base_config()
        self.assertFalse(hasattr(config, "duplicate_position_policy"))

    def test_duplicate_position_policy_constructor_argument_is_rejected(self):
        with self.assertRaises(TypeError):
            self._base_config(duplicate_position_policy="drop-all")


class ReferencePanelBuildConfigFromArgsTest(unittest.TestCase):
    def test_build_parser_defaults_snp_batch_size_to_128(self):
        parser = ref_panel_builder.build_parser()
        self.assertEqual(parser.get_default("snp_batch_size"), 128)

    def test_build_parser_accepts_snp_batch_size(self):
        parser = ref_panel_builder.build_parser()
        args = parser.parse_args([
            "--plink-prefix", "plink/panel.@",
            "--output-dir", "out",
            "--ld-wind-kb", "1",
            "--snp-batch-size", "64",
        ])
        self.assertEqual(args.snp_batch_size, 64)

    def test_build_parser_accepts_hm3_flags(self):
        parser = ref_panel_builder.build_parser()
        args = parser.parse_args([
            "--plink-prefix", "plink/panel.@",
            "--output-dir", "out",
            "--ld-wind-kb", "1",
            "--use-hm3-snps",
            "--use-hm3-quick-liftover",
        ])

        self.assertTrue(args.use_hm3_snps)
        self.assertTrue(args.use_hm3_quick_liftover)

    def test_build_parser_keeps_chunk_size_as_hidden_alias(self):
        parser = ref_panel_builder.build_parser()
        help_text = parser.format_help()
        self.assertIn("--snp-batch-size", help_text)
        self.assertNotIn("--chunk-size", help_text)
        args = parser.parse_args([
            "--plink-prefix", "plink/panel.@",
            "--output-dir", "out",
            "--ld-wind-kb", "1",
            "--chunk-size", "64",
        ])
        self.assertEqual(args.snp_batch_size, 64)

    def test_config_from_args_passes_snp_batch_size_to_config(self):
        parser = ref_panel_builder.build_parser()
        args = parser.parse_args([
            "--plink-prefix", "plink/panel.@",
            "--output-dir", "out",
            "--ld-wind-kb", "1",
            "--snp-batch-size", "64",
        ])
        build_config, _ = ref_panel_builder.config_from_args(args)
        self.assertEqual(build_config.snp_batch_size, 64)

    def test_config_from_args_passes_hm3_flags_to_config(self):
        parser = ref_panel_builder.build_parser()
        args = parser.parse_args([
            "--plink-prefix", "plink/panel.@",
            "--source-genome-build", "hg19",
            "--output-dir", "out",
            "--ld-wind-kb", "1",
            "--use-hm3-snps",
            "--use-hm3-quick-liftover",
        ])

        build_config, _ = ref_panel_builder.config_from_args(args)

        self.assertTrue(build_config.use_hm3_snps)
        self.assertTrue(build_config.use_hm3_quick_liftover)

    def test_build_parser_help_does_not_include_duplicate_position_policy(self):
        parser = ref_panel_builder.build_parser()
        self.assertNotIn("--duplicate-position-policy", parser.format_help())

    def test_build_parser_rejects_duplicate_position_policy(self):
        parser = ref_panel_builder.build_parser()
        with self.assertRaises(SystemExit):
            parser.parse_args([
                "--plink-prefix", "plink/panel.@",
                "--output-dir", "out",
                "--ld-wind-kb", "1",
                "--duplicate-position-policy", "drop-all",
            ])

    def test_build_parser_does_not_accept_genome_build(self):
        parser = ref_panel_builder.build_parser()
        with self.assertRaises(SystemExit):
            parser.parse_args(
                [
                    "--plink-prefix",
                    "plink/panel.@",
                    "--output-dir",
                    "out",
                    "--ld-wind-kb",
                    "1",
                    "--genome-build",
                    "hg38",
                ]
            )

    def test_config_from_args_uses_parser_default_identifier_for_restriction_file(self):
        parser = ref_panel_builder.build_parser()
        args = parser.parse_args(
            [
                "--plink-prefix",
                "plink/panel.@",
                "--genetic-map-hg19-sources",
                "maps/hg19.map",
                "--output-dir",
                "out",
                "--ld-wind-kb",
                "1",
                "--ref-panel-snps-file",
                "hm3.tsv",
            ]
        )

        set_global_config(GlobalConfig(snp_identifier="chr_pos"))
        try:
            build_config, global_config = ref_panel_builder.config_from_args(args)
        finally:
            reset_global_config()

        self.assertEqual(build_config.ref_panel_snps_file, "hm3.tsv")
        self.assertEqual(build_config.source_genome_build, "auto")
        self.assertEqual(global_config.snp_identifier, "chr_pos_allele_aware")
        self.assertEqual(global_config.genome_build, "auto")

    def test_config_from_args_accepts_explicit_auto_source_genome_build(self):
        parser = ref_panel_builder.build_parser()
        args = parser.parse_args(
            [
                "--plink-prefix",
                "plink/panel.@",
                "--source-genome-build",
                "auto",
                "--genetic-map-hg19-sources",
                "maps/hg19.map",
                "--output-dir",
                "out",
                "--ld-wind-kb",
                "1",
            ]
        )

        build_config, _global_config = ref_panel_builder.config_from_args(args)

        self.assertEqual(build_config.source_genome_build, "auto")

    def test_config_from_args_uses_explicit_snp_identifier_for_restriction_file(self):
        parser = ref_panel_builder.build_parser()
        args = parser.parse_args(
            [
                "--plink-prefix",
                "plink/panel.@",
                "--source-genome-build",
                "hg19",
                "--genetic-map-hg19-sources",
                "maps/hg19.map",
                "--output-dir",
                "out",
                "--ld-wind-kb",
                "1",
                "--ref-panel-snps-file",
                "hm3.tsv",
                "--snp-identifier",
                "rsid",
            ]
        )

        build_config, global_config = ref_panel_builder.config_from_args(args)

        self.assertEqual(build_config.source_genome_build, "hg19")
        self.assertIsNone(global_config.genome_build)
        self.assertEqual(global_config.snp_identifier, "rsid")

    def test_config_from_args_does_not_require_genome_build_for_chr_pos_mode(self):
        parser = ref_panel_builder.build_parser()
        args = parser.parse_args(
            [
                "--plink-prefix",
                "plink/panel.@",
                "--source-genome-build",
                "hg19",
                "--genetic-map-hg19-sources",
                "maps/hg19.map",
                "--output-dir",
                "out",
                "--ld-wind-kb",
                "1",
                "--ref-panel-snps-file",
                "hm3.tsv",
                "--snp-identifier",
                "chr_pos",
            ]
        )

        build_config, global_config = ref_panel_builder.config_from_args(args)

        self.assertEqual(build_config.source_genome_build, "hg19")
        self.assertEqual(global_config.snp_identifier, "chr_pos")
        self.assertEqual(global_config.genome_build, "auto")

    def test_run_build_ref_panel_uses_parser_default_identifier_for_restriction_file(self):
        captured = {}

        def fake_run(self, config):
            captured["global_config"] = self.global_config
            captured["config"] = config
            return ref_panel_builder.ReferencePanelBuildResult(panel_name="out", chromosomes=[])

        set_global_config(GlobalConfig(snp_identifier="chr_pos"))
        try:
            with mock.patch.object(ref_panel_builder.ReferencePanelBuilder, "run", fake_run):
                ref_panel_builder.run_build_ref_panel(
                    plink_prefix="plink/panel.@",
                    source_genome_build="hg19",
                    genetic_map_hg19_sources="maps/hg19.map",
                    output_dir="out",
                    ld_wind_kb=1,
                    ref_panel_snps_file="hm3.tsv",
                )
        finally:
            reset_global_config()

        self.assertEqual(captured["config"].ref_panel_snps_file, "hm3.tsv")
        self.assertEqual(captured["global_config"].snp_identifier, "chr_pos_allele_aware")
        self.assertEqual(captured["global_config"].genome_build, "auto")

    def test_run_build_ref_panel_uses_snp_batch_size_keyword(self):
        captured = {}

        def fake_run(self, config):
            captured["config"] = config
            return ref_panel_builder.ReferencePanelBuildResult(panel_name="out", chromosomes=[])

        with mock.patch.object(ref_panel_builder.ReferencePanelBuilder, "run", fake_run):
            ref_panel_builder.run_build_ref_panel(
                plink_prefix="plink/panel.@",
                source_genome_build="hg19",
                genetic_map_hg19_sources="maps/hg19.map",
                output_dir="out",
                ld_wind_kb=1,
                snp_batch_size=64,
            )

        self.assertEqual(captured["config"].snp_batch_size, 64)

    def test_run_build_ref_panel_maps_legacy_chunk_size_keyword(self):
        captured = {}

        def fake_run(self, config):
            captured["config"] = config
            return ref_panel_builder.ReferencePanelBuildResult(panel_name="out", chromosomes=[])

        with mock.patch.object(ref_panel_builder.ReferencePanelBuilder, "run", fake_run):
            ref_panel_builder.run_build_ref_panel(
                plink_prefix="plink/panel.@",
                source_genome_build="hg19",
                genetic_map_hg19_sources="maps/hg19.map",
                output_dir="out",
                ld_wind_kb=1,
                chunk_size=64,
            )

        self.assertEqual(captured["config"].snp_batch_size, 64)

    def test_run_build_ref_panel_rejects_conflicting_batch_size_keywords(self):
        with self.assertRaisesRegex(ValueError, "chunk_size or snp_batch_size"):
            ref_panel_builder.run_build_ref_panel(
                plink_prefix="plink/panel.@",
                source_genome_build="hg19",
                genetic_map_hg19_sources="maps/hg19.map",
                output_dir="out",
                ld_wind_kb=1,
                chunk_size=64,
                snp_batch_size=128,
            )

    def test_run_build_ref_panel_rejects_explicit_genome_build_keyword(self):
        with self.assertRaisesRegex(ValueError, "set_global_config"):
            ref_panel_builder.run_build_ref_panel(
                plink_prefix="plink/panel.@",
                source_genome_build="hg19",
                genetic_map_hg19_sources="maps/hg19.map",
                output_dir="out",
                ld_wind_kb=1,
                ref_panel_snps_file="hm3.tsv",
                snp_identifier="chr_pos",
                genome_build="hg38",
            )

    def test_run_build_ref_panel_rejects_explicit_snp_identifier_keyword(self):
        with self.assertRaisesRegex(ValueError, "set_global_config"):
            ref_panel_builder.run_build_ref_panel(
                plink_prefix="plink/panel.@",
                source_genome_build="hg19",
                genetic_map_hg19_sources="maps/hg19.map",
                output_dir="out",
                ld_wind_kb=1,
                ref_panel_snps_file="hm3.tsv",
                snp_identifier="rsid",
            )


class DuplicatePositionDropAllTest(unittest.TestCase):
    """Unit tests for _resolve_unique_snp_set()."""

    def _make_chrom_df(self, snp_ids, bps):
        """Build a minimal .bim-style DataFrame indexed by PLINK row indices."""
        return pd.DataFrame(
            {"SNP": snp_ids, "BP": bps},
            index=list(range(len(snp_ids))),
        )

    def test_no_duplicates_returns_keep_snps_unchanged(self):
        chrom_df = self._make_chrom_df(["rs1", "rs2", "rs3"], [100, 200, 300])
        keep = np.array([0, 1, 2])
        hg19 = {0: 100, 1: 200, 2: 300}
        hg38 = {0: 1000, 1: 2000, 2: 3000}

        cleaned, dropped = ref_panel_builder._resolve_unique_snp_set("1", chrom_df, keep, hg19, hg38)

        np.testing.assert_array_equal(cleaned, keep)
        self.assertTrue(dropped.empty)

    def test_source_duplicate_drop_all_removes_cluster(self):
        chrom_df = self._make_chrom_df(["rs1", "rs2", "rs3"], [100, 100, 300])
        keep = np.array([0, 1, 2])

        cleaned, dropped = ref_panel_builder._resolve_unique_snp_set("1", chrom_df, keep, {0: 100, 1: 100, 2: 300}, {})

        np.testing.assert_array_equal(cleaned, np.array([2]))
        self.assertEqual(len(dropped), 2)
        self.assertTrue((dropped["reason"] == "source_duplicate").all())
        self.assertTrue(dropped["target_pos"].isna().all())
        self.assertSetEqual(set(dropped["SNP"]), {"rs1", "rs2"})

    def test_target_collision_drop_all_removes_cluster_from_all_builds(self):
        chrom_df = self._make_chrom_df(["rs1", "rs2", "rs3"], [100, 200, 300])
        keep = np.array([0, 1, 2])
        hg19 = {0: 100, 1: 200, 2: 300}
        hg38 = {0: 5000, 1: 5000, 2: 6000}

        cleaned, dropped = ref_panel_builder._resolve_unique_snp_set("1", chrom_df, keep, hg19, hg38)

        np.testing.assert_array_equal(cleaned, np.array([2]))
        self.assertEqual(len(dropped), 2)
        self.assertTrue((dropped["reason"] == "target_collision").all())
        self.assertSetEqual(set(dropped["SNP"]), {"rs1", "rs2"})
        self.assertFalse(dropped["target_pos"].isna().any())

    def test_target_collision_in_hg19_only_drops_both_snps(self):
        chrom_df = self._make_chrom_df(["rs1", "rs2", "rs3"], [100, 200, 300])
        keep = np.array([0, 1, 2])
        hg19 = {0: 5000, 1: 5000, 2: 6000}
        hg38 = {}

        cleaned, dropped = ref_panel_builder._resolve_unique_snp_set("1", chrom_df, keep, hg19, hg38)

        np.testing.assert_array_equal(cleaned, np.array([2]))
        self.assertEqual(len(dropped), 2)
        self.assertTrue((dropped["reason"] == "target_collision").all())

    def test_collision_in_both_builds_deduplicates_provenance_rows(self):
        chrom_df = self._make_chrom_df(["rs1", "rs2", "rs3"], [100, 200, 300])
        keep = np.array([0, 1, 2])
        hg19 = {0: 5000, 1: 5000, 2: 6000}
        hg38 = {0: 7000, 1: 7000, 2: 8000}

        cleaned, dropped = ref_panel_builder._resolve_unique_snp_set("1", chrom_df, keep, hg19, hg38)

        np.testing.assert_array_equal(cleaned, np.array([2]))
        self.assertEqual(len(dropped), 2)
        self.assertSetEqual(set(dropped["SNP"]), {"rs1", "rs2"})

    def test_all_snps_dropped_returns_empty_keep(self):
        chrom_df = self._make_chrom_df(["rs1", "rs2"], [100, 100])
        keep = np.array([0, 1])

        cleaned, dropped = ref_panel_builder._resolve_unique_snp_set("1", chrom_df, keep, {0: 100, 1: 100}, {})

        self.assertEqual(len(cleaned), 0)
        self.assertEqual(len(dropped), 2)

    def test_provenance_columns_are_correct(self):
        chrom_df = self._make_chrom_df(["rs1", "rs2", "rs3"], [100, 100, 300])
        keep = np.array([0, 1, 2])

        _, dropped = ref_panel_builder._resolve_unique_snp_set("1", chrom_df, keep, {0: 100, 1: 100, 2: 300}, {})

        self.assertListEqual(list(dropped.columns), IDENTITY_DROP_COLUMNS)
        self.assertEqual(dropped["CHR"].iloc[0], "1")
        self.assertEqual(dropped["source_pos"].dtype.name, "Int64")
        self.assertEqual(dropped["target_pos"].dtype.name, "Int64")


class ReferencePanelBuilderWorkflowTest(unittest.TestCase):
    DROPPED_SNP_DTYPES = {
        "CHR": "string",
        "SNP": "string",
        "source_pos": "Int64",
        "target_pos": "Int64",
        "reason": "string",
        "base_key": "string",
        "identity_key": "string",
        "allele_set": "string",
        "stage": "string",
    }

    def _read_dropped_sidecar(self, path: Path) -> pd.DataFrame:
        return pd.read_csv(path, sep="\t", compression="gzip", dtype=self.DROPPED_SNP_DTYPES)

    def _write_dummy_plink_prefix(self, root: Path, stem: str, chrom: str):
        prefix = root / stem
        Path(str(prefix) + ".bed").write_bytes(b"")
        Path(str(prefix) + ".fam").write_text("fam iid 0 0 0 -9\n", encoding="utf-8")
        Path(str(prefix) + ".bim").write_text(
            f"{chrom} rs{chrom} 0.0 100 A G\n",
            encoding="utf-8",
        )
        return prefix

    def _write_plink_prefix_rows(self, root: Path, stem: str, rows: list[tuple]):
        prefix = root / stem
        Path(str(prefix) + ".bed").write_bytes(b"")
        Path(str(prefix) + ".fam").write_text("fam iid 0 0 0 -9\n", encoding="utf-8")
        lines = []
        for row in rows:
            chrom, snp, pos, *alleles = row
            a1, a2 = alleles if alleles else ("A", "G")
            lines.append(f"{chrom} {snp} 0.0 {pos} {a1} {a2}\n")
        Path(str(prefix) + ".bim").write_text("".join(lines), encoding="utf-8")
        return prefix

    def _build_config(self, tmpdir: Path) -> ReferencePanelBuildConfig:
        map_hg19 = tmpdir / "hg19.map"
        map_hg38 = tmpdir / "hg38.map"
        map_hg19.write_text("chr position Genetic_Map(cM)\n1 100 0.0\n2 100 0.0\n", encoding="utf-8")
        map_hg38.write_text("chr position Genetic_Map(cM)\n1 100 0.0\n2 100 0.0\n", encoding="utf-8")
        return ReferencePanelBuildConfig(
            plink_prefix=tmpdir / "panel.@",
            source_genome_build="hg19",
            genetic_map_hg19_sources=map_hg19,
            genetic_map_hg38_sources=map_hg38,
            liftover_chain_hg19_to_hg38_file=tmpdir / "hg19ToHg38.over.chain",
            output_dir=tmpdir / "out",
            ld_wind_kb=1.0,
        )

    def test_drop_all_policy_writes_sidecar_under_dropped_snps(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            self._write_plink_prefix_rows(
                tmpdir,
                "panel.1",
                [
                    ("1", "rs1", 100),
                    ("1", "rs2", 100),
                    ("1", "rs3", 200),
                ],
            )
            config = ReferencePanelBuildConfig(
                plink_prefix=tmpdir / "panel.@",
                source_genome_build="hg19",
                output_dir=tmpdir / "out",
                ld_wind_kb=1.0,
            )
            builder = ref_panel_builder.ReferencePanelBuilder(
                global_config=GlobalConfig(snp_identifier="chr_pos", genome_build="hg19")
            )

            def fake_resolve(*, chrom_df, keep_snps, **_kwargs):
                keep = np.asarray(keep_snps, dtype=int)
                hg19 = {int(idx): int(chrom_df.loc[idx, "BP"]) for idx in keep}
                return keep, hg19, {}, ref_panel_builder._empty_unified_drop_frame()

            with mock.patch.object(
                builder,
                "_resolve_mappable_snp_positions",
                side_effect=fake_resolve,
            ), mock.patch(
                "ldsc.ref_panel_builder.kernel_builder.write_r2_parquet"
            ), mock.patch(
                "ldsc.ref_panel_builder.kernel_builder.write_runtime_metadata_sidecar"
            ), mock.patch(
                "ldsc.ref_panel_builder.kernel_ldscore.PlinkBEDFile"
            ) as mock_bed:
                mock_bed.return_value.kept_snps = [2]
                mock_bed.return_value.maf = np.array([0.3])
                mock_bed.return_value.m = 1
                mock_bed.return_value.n = 1
                mock_bed.return_value.nextSNPs = lambda: iter([np.array([0.0])])

                with self.assertLogs("LDSC.ref_panel_builder", level="INFO") as log_ctx:
                    builder.run(config)

            sidecar = tmpdir / "out" / "diagnostics" / "dropped_snps" / "chr1_dropped.tsv.gz"
            self.assertTrue(sidecar.exists(), "sidecar must be written under dropped_snps")
            with gzip.open(sidecar, "rt") as fh:
                dropped_df = pd.read_csv(fh, sep="\t")
            self.assertEqual(len(dropped_df), 2)
            self.assertTrue((dropped_df["reason"] == "duplicate_identity").all())
            self.assertTrue((dropped_df["stage"] == "ref_panel_source_identity_cleanup").all())
            self.assertTrue(any("chr1_dropped.tsv.gz" in line for line in log_ctx.output))

    def test_drop_all_policy_logs_examples_only_at_debug(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            self._write_plink_prefix_rows(
                tmpdir,
                "panel.1",
                [
                    ("1", "rs1", 100),
                    ("1", "rs2", 100),
                    ("1", "rs3", 200),
                ],
            )
            config = ReferencePanelBuildConfig(
                plink_prefix=tmpdir / "panel.@",
                source_genome_build="hg19",
                output_dir=tmpdir / "out",
                ld_wind_kb=1.0,
            )
            builder = ref_panel_builder.ReferencePanelBuilder(
                global_config=GlobalConfig(
                    snp_identifier="chr_pos",
                    genome_build="hg19",
                    log_level="DEBUG",
                )
            )

            def fake_resolve(*, chrom_df, keep_snps, **_kwargs):
                keep = np.asarray(keep_snps, dtype=int)
                hg19 = {int(idx): int(chrom_df.loc[idx, "BP"]) for idx in keep}
                return keep, hg19, {}, ref_panel_builder._empty_unified_drop_frame()

            with mock.patch.object(
                builder,
                "_resolve_mappable_snp_positions",
                side_effect=fake_resolve,
            ), mock.patch(
                "ldsc.ref_panel_builder.kernel_builder.write_r2_parquet"
            ), mock.patch(
                "ldsc.ref_panel_builder.kernel_builder.write_runtime_metadata_sidecar"
            ), mock.patch(
                "ldsc.ref_panel_builder.kernel_ldscore.PlinkBEDFile"
            ) as mock_bed:
                mock_bed.return_value.kept_snps = [2]
                mock_bed.return_value.maf = np.array([0.3])
                mock_bed.return_value.m = 1
                mock_bed.return_value.n = 1
                mock_bed.return_value.nextSNPs = lambda: iter([np.array([0.0])])

                with self.assertLogs("LDSC.ref_panel_builder", level="DEBUG") as log_ctx:
                    builder.run(config)

            info_text = "\n".join(
                record.getMessage() for record in log_ctx.records if record.levelno == logging.INFO
            )
            debug_text = "\n".join(
                record.getMessage() for record in log_ctx.records if record.levelno == logging.DEBUG
            )
            self.assertIn("duplicate_identity", info_text)
            self.assertNotIn("rs1", info_text)
            self.assertIn("duplicate_identity", debug_text)
            self.assertIn("rs1", debug_text)

    def test_resolve_mappable_snp_positions_returns_unmapped_and_cross_chromosome_drop_frame(self):
        builder = ref_panel_builder.ReferencePanelBuilder(global_config=GlobalConfig(snp_identifier="chr_pos"))
        build_state = ref_panel_builder._BuildState(
            genetic_map_hg19=None,
            genetic_map_hg38=None,
            liftover_chain_paths={("hg19", "hg38"): "chain.over"},
        )
        chrom_df = pd.DataFrame({"SNP": ["rs1", "rs2", "rs3"], "BP": [100, 200, 300]}, index=[0, 1, 2])
        mapping = kernel_builder.LiftOverMappingResult(
            translated_positions=np.asarray([1300], dtype=np.int64),
            keep_mask=np.asarray([False, False, True], dtype=bool),
            unmapped_count=1,
            cross_chrom_count=1,
            unmapped_mask=np.asarray([True, False, False], dtype=bool),
            cross_chrom_mask=np.asarray([False, True, False], dtype=bool),
        )

        with mock.patch.object(builder, "_map_positions", return_value=mapping):
            retained, hg19_lookup, hg38_lookup, drop_frame = builder._resolve_mappable_snp_positions(
                build_state=build_state,
                chrom="1",
                source_build="hg19",
                chrom_df=chrom_df,
                keep_snps=np.array([0, 1, 2], dtype=int),
            )

        self.assertEqual(retained.tolist(), [2])
        self.assertEqual(hg19_lookup, {2: 300})
        self.assertEqual(hg38_lookup, {2: 1300})
        self.assertEqual(drop_frame["reason"].tolist(), ["unmapped_liftover", "cross_chromosome_liftover"])
        self.assertEqual({column: str(dtype) for column, dtype in drop_frame.dtypes.items()}, self.DROPPED_SNP_DTYPES)

    def test_dropped_snps_sidecar_written_as_header_only_when_no_drops(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            self._write_plink_prefix_rows(tmpdir, "panel.1", [("1", "rs1", 100)])
            config = ReferencePanelBuildConfig(
                plink_prefix=tmpdir / "panel.@",
                source_genome_build="hg19",
                output_dir=tmpdir / "out",
                ld_wind_kb=1.0,
            )
            builder = ref_panel_builder.ReferencePanelBuilder(
                global_config=GlobalConfig(snp_identifier="chr_pos", genome_build="hg19")
            )

            with mock.patch("ldsc.ref_panel_builder.kernel_builder.write_r2_parquet"), mock.patch(
                "ldsc.ref_panel_builder.kernel_builder.write_runtime_metadata_sidecar"
            ), mock.patch("ldsc.ref_panel_builder.kernel_ldscore.PlinkBEDFile") as mock_bed:
                mock_bed.return_value.kept_snps = [0]
                mock_bed.return_value.maf = np.array([0.3])
                mock_bed.return_value.m = 1
                mock_bed.return_value.n = 1
                mock_bed.return_value.nextSNPs = lambda: iter([np.zeros((1, 1))])

                builder.run(config)

            sidecar = tmpdir / "out" / "diagnostics" / "dropped_snps" / "chr1_dropped.tsv.gz"
            self.assertTrue(sidecar.exists())
            dropped = self._read_dropped_sidecar(sidecar)
            self.assertEqual(len(dropped), 0)
            self.assertEqual(dropped.columns.tolist(), IDENTITY_DROP_COLUMNS)
            self.assertEqual({column: str(dtype) for column, dtype in dropped.dtypes.items()}, self.DROPPED_SNP_DTYPES)

    def test_dropped_snps_sidecar_overwritten_in_place_on_rerun(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            self._write_plink_prefix_rows(tmpdir, "panel.1", [("1", "rs1", 100)])
            sidecar = tmpdir / "out" / "diagnostics" / "dropped_snps" / "chr1_dropped.tsv.gz"
            sidecar.parent.mkdir(parents=True)
            sidecar.write_text("stale\n", encoding="utf-8")
            config = ReferencePanelBuildConfig(
                plink_prefix=tmpdir / "panel.@",
                source_genome_build="hg19",
                output_dir=tmpdir / "out",
                ld_wind_kb=1.0,
                overwrite=True,
            )
            builder = ref_panel_builder.ReferencePanelBuilder(
                global_config=GlobalConfig(snp_identifier="chr_pos", genome_build="hg19")
            )

            with mock.patch("ldsc.ref_panel_builder.kernel_builder.write_r2_parquet"), mock.patch(
                "ldsc.ref_panel_builder.kernel_builder.write_runtime_metadata_sidecar"
            ), mock.patch("ldsc.ref_panel_builder.kernel_ldscore.PlinkBEDFile") as mock_bed:
                mock_bed.return_value.kept_snps = [0]
                mock_bed.return_value.maf = np.array([0.3])
                mock_bed.return_value.m = 1
                mock_bed.return_value.n = 1
                mock_bed.return_value.nextSNPs = lambda: iter([np.zeros((1, 1))])

                builder.run(config)

            dropped = self._read_dropped_sidecar(sidecar)
            self.assertEqual(len(dropped), 0)

    def test_dropped_snps_sidecar_written_header_only_on_restriction_empty_skip(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            prefix = self._write_plink_prefix_rows(tmpdir, "panel.1", [("1", "rs1", 100)])
            restriction = tmpdir / "restrict.tsv"
            restriction.write_text("SNP\nrs_not_present\n", encoding="utf-8")
            config = ReferencePanelBuildConfig(
                plink_prefix=prefix,
                source_genome_build="hg19",
                output_dir=tmpdir / "out",
                ld_wind_kb=1.0,
                ref_panel_snps_file=restriction,
            )
            builder = ref_panel_builder.ReferencePanelBuilder(
                global_config=GlobalConfig(snp_identifier="rsid")
            )

            with self.assertRaisesRegex(ValueError, "No chromosome artifacts"):
                builder.run(config)

            sidecar = tmpdir / "out" / "diagnostics" / "dropped_snps" / "chr1_dropped.tsv.gz"
            self.assertTrue(sidecar.exists())
            self.assertEqual(len(self._read_dropped_sidecar(sidecar)), 0)

    def test_rsid_source_only_build_logs_coordinate_duplicate_filtering_skipped_once(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            self._write_plink_prefix_rows(
                tmpdir,
                "panel.1",
                [
                    ("1", "rs1", 100),
                    ("1", "rs2", 100),
                    ("1", "rs3", 200),
                ],
            )
            config = ReferencePanelBuildConfig(
                plink_prefix=tmpdir / "panel.@",
                source_genome_build="hg19",
                output_dir=tmpdir / "out",
                ld_wind_kb=1.0,
            )
            builder = ref_panel_builder.ReferencePanelBuilder(
                global_config=GlobalConfig(snp_identifier="rsid", genome_build="hg38")
            )

            def fake_resolve(*, chrom_df, keep_snps, **_kwargs):
                keep = np.asarray(keep_snps, dtype=int)
                hg19 = {int(idx): int(chrom_df.loc[idx, "BP"]) for idx in keep}
                return keep, hg19, {}, ref_panel_builder._empty_unified_drop_frame()

            with mock.patch.object(
                builder,
                "_resolve_mappable_snp_positions",
                side_effect=fake_resolve,
            ), mock.patch(
                "ldsc.ref_panel_builder.kernel_builder.write_r2_parquet"
            ), mock.patch(
                "ldsc.ref_panel_builder.kernel_builder.write_runtime_metadata_sidecar"
            ), mock.patch(
                "ldsc.ref_panel_builder.kernel_ldscore.PlinkBEDFile"
            ) as mock_bed:
                mock_bed.return_value.kept_snps = [0, 1, 2]
                mock_bed.return_value.maf = np.array([0.2, 0.2, 0.3])
                mock_bed.return_value.m = 3
                mock_bed.return_value.n = 1
                mock_bed.return_value.nextSNPs = lambda: iter([np.zeros((1, 3))])

                with self.assertLogs("LDSC.ref_panel_builder", level="INFO") as log_ctx:
                    result = builder.run(config)

            self.assertEqual(result.chromosomes, ["1"])
            self.assertTrue((tmpdir / "out" / "diagnostics" / "dropped_snps" / "chr1_dropped.tsv.gz").exists())
            self.assertEqual(len(self._read_dropped_sidecar(tmpdir / "out" / "diagnostics" / "dropped_snps" / "chr1_dropped.tsv.gz")), 0)
            self.assertEqual(mock_bed.call_args.kwargs["keep_snps"], [0, 1, 2])
            messages = "\n".join(log_ctx.output)
            self.assertEqual(messages.count("Coordinate duplicate filtering applies only for chr_pos-family snp_identifier modes"), 1)

    def test_builder_run_rejects_matching_chain_in_rsid_mode(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            self._write_dummy_plink_prefix(tmpdir, "panel.1", "1")
            config = ReferencePanelBuildConfig(
                plink_prefix=tmpdir / "panel.@",
                source_genome_build="hg19",
                liftover_chain_hg19_to_hg38_file=tmpdir / "hg19ToHg38.over.chain",
                output_dir=tmpdir / "out",
                ld_wind_kb=1.0,
            )
            builder = ref_panel_builder.ReferencePanelBuilder(
                global_config=GlobalConfig(snp_identifier="rsid", genome_build="hg38")
            )

            with mock.patch.object(
                ref_panel_builder.ReferencePanelBuilder,
                "_build_chromosome",
                side_effect=AssertionError("chromosome build should not run"),
            ):
                with self.assertRaisesRegex(ValueError, "chain liftover.*chr_pos-family"):
                    builder.run(config)

    def test_builder_run_rejects_matching_chain_in_rsid_allele_aware_mode(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            self._write_dummy_plink_prefix(tmpdir, "panel.1", "1")
            config = ReferencePanelBuildConfig(
                plink_prefix=tmpdir / "panel.@",
                source_genome_build="hg19",
                liftover_chain_hg19_to_hg38_file=tmpdir / "hg19ToHg38.over.chain",
                output_dir=tmpdir / "out",
                ld_wind_kb=1.0,
            )
            builder = ref_panel_builder.ReferencePanelBuilder(
                global_config=GlobalConfig(snp_identifier="rsid_allele_aware")
            )

            with self.assertRaisesRegex(
                ValueError,
                "Reference-panel chain liftover is only valid for chr_pos-family modes; "
                "omit the matching liftover chain in rsID-family modes.",
            ):
                builder.run(config)

    def test_prepare_build_state_rejects_hm3_quick_liftover_in_rsid_family_modes(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            config = ReferencePanelBuildConfig(
                plink_prefix=tmpdir / "panel.@",
                source_genome_build="hg19",
                output_dir=tmpdir / "out",
                ld_wind_kb=1.0,
                use_hm3_snps=True,
                use_hm3_quick_liftover=True,
            )
            for mode in ("rsid", "rsid_allele_aware"):
                builder = ref_panel_builder.ReferencePanelBuilder(global_config=GlobalConfig(snp_identifier=mode))
                with self.assertRaisesRegex(
                    ValueError,
                    "Reference-panel HM3 quick liftover is only valid for chr_pos-family modes.",
                ):
                    builder._prepare_build_state(config)

    def test_allele_aware_mode_requires_bim_allele_columns(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            prefix = tmpdir / "panel.1"
            Path(str(prefix) + ".bed").write_bytes(b"")
            Path(str(prefix) + ".fam").write_text("fam iid 0 0 0 -9\n", encoding="utf-8")
            Path(str(prefix) + ".bim").write_text("1 rs1 0.0 100\n", encoding="utf-8")
            config = ReferencePanelBuildConfig(
                plink_prefix=tmpdir / "panel.@",
                source_genome_build="hg19",
                output_dir=tmpdir / "out",
                ld_wind_kb=1.0,
            )
            builder = ref_panel_builder.ReferencePanelBuilder(
                global_config=GlobalConfig(snp_identifier="chr_pos_allele_aware")
            )

            with self.assertRaisesRegex(ValueError, "A1.*A2|columns do not match"):
                builder.run(config)

    def test_chr_pos_allele_aware_drops_invalid_allele_rows_before_genotype_loading(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            prefix = self._write_plink_prefix_rows(
                tmpdir,
                "panel.1",
                [
                    ("1", "rs1", 100, "A", "C"),
                    ("1", "rs_bad", 200, "A", "N"),
                    ("1", "rs2", 300, "A", "G"),
                ],
            )
            config = ReferencePanelBuildConfig(
                plink_prefix=tmpdir / "panel.@",
                source_genome_build="hg19",
                output_dir=tmpdir / "out",
                ld_wind_kb=1.0,
            )
            build_state = ref_panel_builder._BuildState(genetic_map_hg19=None, genetic_map_hg38=None)
            builder = ref_panel_builder.ReferencePanelBuilder(
                global_config=GlobalConfig(snp_identifier="chr_pos_allele_aware")
            )
            captured = {}

            def stop_after_identity_cleanup(*, keep_snps, **_kwargs):
                captured["keep_snps"] = list(map(int, keep_snps))
                raise RuntimeError("stop after identity cleanup")

            with mock.patch.object(builder, "_resolve_mappable_snp_positions", side_effect=stop_after_identity_cleanup):
                with self.assertRaisesRegex(RuntimeError, "stop after identity cleanup"):
                    builder._build_chromosome(str(prefix), "1", config, build_state)

            self.assertEqual(captured["keep_snps"], [0, 2])

    def test_allele_aware_modes_drop_multiallelic_base_key_clusters(self):
        cases = [
            ("chr_pos_allele_aware", [("1", "rs1", 100, "A", "C"), ("1", "rs2", 100, "A", "G"), ("1", "rs3", 300, "A", "C")]),
            ("rsid_allele_aware", [("1", "rs1", 100, "A", "C"), ("1", "rs1", 200, "A", "G"), ("1", "rs3", 300, "A", "C")]),
        ]
        for mode, rows in cases:
            with self.subTest(mode=mode):
                with tempfile.TemporaryDirectory() as tmpdir:
                    tmpdir = Path(tmpdir)
                    prefix = self._write_plink_prefix_rows(tmpdir, "panel.1", rows)
                    config = ReferencePanelBuildConfig(
                        plink_prefix=tmpdir / "panel.@",
                        source_genome_build="hg19",
                        output_dir=tmpdir / "out",
                        ld_wind_kb=1.0,
                    )
                    build_state = ref_panel_builder._BuildState(genetic_map_hg19=None, genetic_map_hg38=None)
                    builder = ref_panel_builder.ReferencePanelBuilder(global_config=GlobalConfig(snp_identifier=mode))
                    captured = {}

                    def stop_after_identity_cleanup(*, keep_snps, **_kwargs):
                        captured["keep_snps"] = list(map(int, keep_snps))
                        raise RuntimeError("stop after identity cleanup")

                    with mock.patch.object(builder, "_resolve_mappable_snp_positions", side_effect=stop_after_identity_cleanup):
                        with self.assertRaisesRegex(RuntimeError, "stop after identity cleanup"):
                            builder._build_chromosome(str(prefix), "1", config, build_state)

                    self.assertEqual(captured["keep_snps"], [2])

    def test_rsid_base_mode_drops_duplicate_effective_keys_with_identity_reason(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            prefix = self._write_plink_prefix_rows(
                tmpdir,
                "panel.1",
                [("1", "rs1", 100, "A", "C"), ("1", "rs1", 200, "N", "N"), ("1", "rs2", 300, "A", "C")],
            )
            config = ReferencePanelBuildConfig(
                plink_prefix=tmpdir / "panel.@",
                source_genome_build="hg19",
                output_dir=tmpdir / "out",
                ld_wind_kb=1.0,
            )
            build_state = ref_panel_builder._BuildState(genetic_map_hg19=None, genetic_map_hg38=None)
            builder = ref_panel_builder.ReferencePanelBuilder(global_config=GlobalConfig(snp_identifier="rsid"))
            captured = {}

            def stop_after_identity_cleanup(*, keep_snps, **_kwargs):
                captured["keep_snps"] = list(map(int, keep_snps))
                raise RuntimeError("stop after identity cleanup")

            with mock.patch.object(builder, "_resolve_mappable_snp_positions", side_effect=stop_after_identity_cleanup):
                with self.assertRaisesRegex(RuntimeError, "stop after identity cleanup"):
                    builder._build_chromosome(str(prefix), "1", config, build_state)

            self.assertEqual(captured["keep_snps"], [2])

    def test_all_dropped_chromosome_is_skipped_without_poisoning_other_chromosomes(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            self._write_plink_prefix_rows(tmpdir, "panel.1", [("1", "rs1", 100), ("1", "rs2", 100)])
            self._write_plink_prefix_rows(tmpdir, "panel.2", [("2", "rs3", 200)])
            config = ReferencePanelBuildConfig(
                plink_prefix=tmpdir / "panel.@",
                source_genome_build="hg19",
                output_dir=tmpdir / "out",
                ld_wind_kb=1.0,
            )
            builder = ref_panel_builder.ReferencePanelBuilder(
                global_config=GlobalConfig(snp_identifier="chr_pos", genome_build="hg19")
            )

            def fake_resolve(*, chrom_df, keep_snps, **_kwargs):
                keep = np.asarray(keep_snps, dtype=int)
                hg19 = {int(idx): int(chrom_df.loc[idx, "BP"]) for idx in keep}
                return keep, hg19, {}, ref_panel_builder._empty_unified_drop_frame()

            with mock.patch.object(
                builder,
                "_resolve_mappable_snp_positions",
                side_effect=fake_resolve,
            ), mock.patch(
                "ldsc.ref_panel_builder.kernel_builder.write_r2_parquet"
            ), mock.patch(
                "ldsc.ref_panel_builder.kernel_builder.write_runtime_metadata_sidecar"
            ), mock.patch(
                "ldsc.ref_panel_builder.kernel_ldscore.PlinkBEDFile"
            ) as mock_bed:
                mock_bed.return_value.kept_snps = [0]
                mock_bed.return_value.maf = np.array([0.3])
                mock_bed.return_value.m = 1
                mock_bed.return_value.n = 1
                mock_bed.return_value.nextSNPs = lambda: iter([np.zeros((1, 1))])

                result = builder.run(config)

            self.assertEqual(result.chromosomes, ["2"])
            self.assertTrue((tmpdir / "out" / "diagnostics" / "dropped_snps" / "chr1_dropped.tsv.gz").exists())
            self.assertTrue((tmpdir / "out" / "diagnostics" / "dropped_snps" / "chr2_dropped.tsv.gz").exists())
            self.assertEqual(len(self._read_dropped_sidecar(tmpdir / "out" / "diagnostics" / "dropped_snps" / "chr2_dropped.tsv.gz")), 0)

    def test_builder_run_collects_artifact_paths_from_resolved_suite(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            self._write_dummy_plink_prefix(tmpdir, "panel.1", "1")
            self._write_dummy_plink_prefix(tmpdir, "panel.2", "2")
            config = self._build_config(tmpdir)
            builder = ref_panel_builder.ReferencePanelBuilder(global_config=GlobalConfig(snp_identifier="chr_pos", genome_build="hg38"))

            def fake_build(prefix, chrom, config, build_state):
                out_root = Path(config.output_dir)
                return {
                    "r2_hg19": str(out_root / "hg19" / f"chr{chrom}_r2.parquet"),
                    "r2_hg38": str(out_root / "hg38" / f"chr{chrom}_r2.parquet"),
                    "meta_hg19": str(out_root / "hg19" / f"chr{chrom}_meta.tsv.gz"),
                    "meta_hg38": str(out_root / "hg38" / f"chr{chrom}_meta.tsv.gz"),
                }

            with mock.patch.object(
                ref_panel_builder.ReferencePanelBuilder,
                "_build_chromosome",
                side_effect=fake_build,
            ) as patched:
                result = builder.run(config)

            self.assertEqual(result.panel_name, "out")
            self.assertEqual(result.chromosomes, ["1", "2"])
            self.assertEqual(patched.call_count, 2)
            self.assertNotIn("ann", result.output_paths)
            self.assertNotIn("ld", result.output_paths)
            self.assertEqual(result.output_paths["metadata"], [str(tmpdir / "out" / "diagnostics" / "metadata.json")])
            self.assertFalse((tmpdir / "out" / "metadata.json").exists())
            metadata = json.loads((tmpdir / "out" / "diagnostics" / "metadata.json").read_text(encoding="utf-8"))
            self.assertEqual(metadata["artifact_type"], "ref_panel")
            self.assertEqual(metadata["snp_identifier"], "chr_pos")
            self.assertEqual(metadata["genome_build"], "hg38")
            self.assertEqual(metadata["chromosomes"], ["1", "2"])
            self.assertEqual(metadata["files"]["r2_hg19"], ["hg19/chr1_r2.parquet", "hg19/chr2_r2.parquet"])
            self.assertEqual(metadata["files"]["meta_hg38"], ["hg38/chr1_meta.tsv.gz", "hg38/chr2_meta.tsv.gz"])
            self.assertEqual(
                result.output_paths["r2_hg19"],
                [
                    str(tmpdir / "out" / "hg19" / "chr1_r2.parquet"),
                    str(tmpdir / "out" / "hg19" / "chr2_r2.parquet"),
                ],
            )
            self.assertEqual(
                result.output_paths["meta_hg38"],
                [
                    str(tmpdir / "out" / "hg38" / "chr1_meta.tsv.gz"),
                    str(tmpdir / "out" / "hg38" / "chr2_meta.tsv.gz"),
                ],
            )

    def test_run_build_ref_panel_from_args_writes_workflow_log(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            self._write_dummy_plink_prefix(tmpdir, "panel.1", "1")
            parser = ref_panel_builder.build_parser()
            args = parser.parse_args(
                [
                    "--plink-prefix",
                    str(tmpdir / "panel.@"),
                    "--source-genome-build",
                    "hg19",
                    "--output-dir",
                    str(tmpdir / "out"),
                    "--ld-wind-kb",
                    "1",
                ]
            )

            def fake_build(prefix, chrom, config, build_state):
                out_root = Path(config.output_dir)
                return {
                    "r2_hg19": str(out_root / "hg19" / f"chr{chrom}_r2.parquet"),
                    "meta_hg19": str(out_root / "hg19" / f"chr{chrom}_meta.tsv.gz"),
                }

            with mock.patch.object(
                ref_panel_builder.ReferencePanelBuilder,
                "_build_chromosome",
                side_effect=fake_build,
            ):
                result = ref_panel_builder.run_build_ref_panel_from_args(args)

            self.assertTrue((tmpdir / "out" / "diagnostics" / "build-ref-panel.log").exists())
            self.assertNotIn("log", result.output_paths)

    def test_run_build_ref_panel_from_args_single_concrete_prefix_writes_chromosome_scoped_log(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            self._write_dummy_plink_prefix(tmpdir, "panel.1", "1")
            parser = ref_panel_builder.build_parser()
            args = parser.parse_args(
                [
                    "--plink-prefix",
                    str(tmpdir / "panel.1"),
                    "--source-genome-build",
                    "hg19",
                    "--output-dir",
                    str(tmpdir / "out"),
                    "--ld-wind-kb",
                    "1",
                ]
            )

            def fake_build(prefix, chrom, config, build_state):
                out_root = Path(config.output_dir)
                return {
                    "r2_hg19": str(out_root / "hg19" / f"chr{chrom}_r2.parquet"),
                    "meta_hg19": str(out_root / "hg19" / f"chr{chrom}_meta.tsv.gz"),
                }

            with mock.patch.object(
                ref_panel_builder.ReferencePanelBuilder,
                "_build_chromosome",
                side_effect=fake_build,
            ):
                result = ref_panel_builder.run_build_ref_panel_from_args(args)

            self.assertTrue((tmpdir / "out" / "diagnostics" / "build-ref-panel.chr1.log").exists())
            self.assertFalse((tmpdir / "out" / "diagnostics" / "build-ref-panel.log").exists())
            self.assertNotIn("log", result.output_paths)

    def test_builder_run_allows_source_only_output_paths_when_no_chain_is_available(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            self._write_dummy_plink_prefix(tmpdir, "panel.1", "1")
            map_hg19 = tmpdir / "hg19.map"
            map_hg19.write_text("chr position Genetic_Map(cM)\n1 100 0.0\n", encoding="utf-8")
            config = ReferencePanelBuildConfig(
                plink_prefix=tmpdir / "panel.@",
                source_genome_build="hg19",
                genetic_map_hg19_sources=map_hg19,
                genetic_map_hg38_sources=None,
                output_dir=tmpdir / "out",
                ld_wind_kb=1.0,
            )
            builder = ref_panel_builder.ReferencePanelBuilder(global_config=GlobalConfig(snp_identifier="chr_pos", genome_build="hg38"))

            def fake_build(prefix, chrom, config, build_state):
                out_root = Path(config.output_dir)
                return {
                    "r2_hg19": str(out_root / "hg19" / f"chr{chrom}_r2.parquet"),
                    "meta_hg19": str(out_root / "hg19" / f"chr{chrom}_meta.tsv.gz"),
                }

            with self.assertLogs("LDSC.ref_panel_builder", level="INFO") as logs:
                with mock.patch.object(
                    ref_panel_builder.ReferencePanelBuilder,
                    "_build_chromosome",
                    side_effect=fake_build,
                ):
                    result = builder.run(config)

            self.assertIn("meta_hg19", result.output_paths)
            self.assertIn("r2_hg19", result.output_paths)
            self.assertNotIn("ann", result.output_paths)
            self.assertNotIn("ld", result.output_paths)
            self.assertNotIn("meta_hg38", result.output_paths)
            self.assertNotIn("r2_hg38", result.output_paths)
            self.assertTrue(any("source-build-only" in message for message in logs.output))

    def test_builder_run_refuses_existing_candidate_artifact_before_chromosome_build(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            self._write_dummy_plink_prefix(tmpdir, "panel.1", "1")
            config = self._build_config(tmpdir)
            existing = tmpdir / "out" / "hg19" / "chr1_r2.parquet"
            existing.parent.mkdir(parents=True)
            existing.write_text("existing\n", encoding="utf-8")
            builder = ref_panel_builder.ReferencePanelBuilder(global_config=GlobalConfig(snp_identifier="chr_pos", genome_build="hg38"))

            with mock.patch.object(
                ref_panel_builder.ReferencePanelBuilder,
                "_build_chromosome",
                side_effect=AssertionError("chromosome build should not run"),
            ):
                with self.assertRaisesRegex(FileExistsError, "overwrite"):
                    builder.run(config)

            self.assertEqual(existing.read_text(encoding="utf-8"), "existing\n")

    def test_builder_run_allows_existing_candidate_artifact_with_overwrite(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            self._write_dummy_plink_prefix(tmpdir, "panel.1", "1")
            config = dataclass_replace(self._build_config(tmpdir), overwrite=True)
            existing = tmpdir / "out" / "hg19" / "chr1_r2.parquet"
            existing.parent.mkdir(parents=True)
            existing.write_text("existing\n", encoding="utf-8")
            builder = ref_panel_builder.ReferencePanelBuilder(global_config=GlobalConfig(snp_identifier="chr_pos", genome_build="hg38"))

            with mock.patch.object(
                ref_panel_builder.ReferencePanelBuilder,
                "_build_chromosome",
                return_value={"r2_hg19": str(existing), "r2_hg38": "r2", "meta_hg19": "m19", "meta_hg38": "m38"},
            ) as patched:
                result = builder.run(config)

            patched.assert_called_once()
            self.assertEqual(result.chromosomes, ["1"])

    def test_builder_run_ignores_legacy_root_diagnostics_without_overwrite(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            self._write_dummy_plink_prefix(tmpdir, "panel.1", "1")
            config = self._build_config(tmpdir)
            output_dir = tmpdir / "out"
            root_drop_dir = output_dir / "dropped_snps"
            root_drop_dir.mkdir(parents=True)
            (output_dir / "metadata.json").write_text('{"legacy": true}\n', encoding="utf-8")
            (output_dir / "build-ref-panel.log").write_text("legacy log\n", encoding="utf-8")
            (root_drop_dir / "chr1_dropped.tsv.gz").write_text("legacy drops\n", encoding="utf-8")
            builder = ref_panel_builder.ReferencePanelBuilder(
                global_config=GlobalConfig(snp_identifier="chr_pos", genome_build="hg38")
            )
            builder._workflow_log_path = output_dir / "diagnostics" / "build-ref-panel.log"

            with mock.patch.object(
                ref_panel_builder.ReferencePanelBuilder,
                "_build_chromosome",
                return_value={
                    "r2_hg19": str(output_dir / "hg19" / "chr1_r2.parquet"),
                    "r2_hg38": str(output_dir / "hg38" / "chr1_r2.parquet"),
                    "meta_hg19": str(output_dir / "hg19" / "chr1_meta.tsv.gz"),
                    "meta_hg38": str(output_dir / "hg38" / "chr1_meta.tsv.gz"),
                },
            ):
                result = builder.run(config)

            self.assertEqual((output_dir / "metadata.json").read_text(encoding="utf-8"), '{"legacy": true}\n')
            self.assertEqual((output_dir / "build-ref-panel.log").read_text(encoding="utf-8"), "legacy log\n")
            self.assertEqual((root_drop_dir / "chr1_dropped.tsv.gz").read_text(encoding="utf-8"), "legacy drops\n")
            self.assertEqual(result.chromosomes, ["1"])
            self.assertTrue((output_dir / "diagnostics" / "metadata.json").exists())
            self.assertTrue((output_dir / "diagnostics" / "build-ref-panel.log").exists())

    def test_builder_run_refuses_stale_owned_artifact_before_chromosome_build(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            self._write_dummy_plink_prefix(tmpdir, "panel.1", "1")
            config = self._build_config(tmpdir)
            stale = tmpdir / "out" / "hg38" / "chr2_r2.parquet"
            stale.parent.mkdir(parents=True)
            stale.write_text("stale\n", encoding="utf-8")
            builder = ref_panel_builder.ReferencePanelBuilder(global_config=GlobalConfig(snp_identifier="chr_pos", genome_build="hg38"))

            with mock.patch.object(
                ref_panel_builder.ReferencePanelBuilder,
                "_build_chromosome",
                side_effect=AssertionError("chromosome build should not run"),
            ):
                with self.assertRaisesRegex(FileExistsError, "overwrite"):
                    builder.run(config)

            self.assertEqual(stale.read_text(encoding="utf-8"), "stale\n")

    def test_builder_run_overwrite_removes_stale_owned_artifact_after_success(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            self._write_dummy_plink_prefix(tmpdir, "panel.1", "1")
            config = dataclass_replace(self._build_config(tmpdir), overwrite=True)
            stale_r2 = tmpdir / "out" / "hg38" / "chr2_r2.parquet"
            stale_drop = tmpdir / "out" / "diagnostics" / "dropped_snps" / "chr2_dropped.tsv.gz"
            stale_r2.parent.mkdir(parents=True)
            stale_drop.parent.mkdir(parents=True)
            stale_r2.write_text("stale\n", encoding="utf-8")
            stale_drop.write_text("stale\n", encoding="utf-8")
            builder = ref_panel_builder.ReferencePanelBuilder(global_config=GlobalConfig(snp_identifier="chr_pos", genome_build="hg38"))

            with mock.patch.object(
                ref_panel_builder.ReferencePanelBuilder,
                "_build_chromosome",
                return_value={
                    "r2_hg19": str(tmpdir / "out" / "hg19" / "chr1_r2.parquet"),
                    "r2_hg38": str(tmpdir / "out" / "hg38" / "chr1_r2.parquet"),
                    "meta_hg19": str(tmpdir / "out" / "hg19" / "chr1_meta.tsv.gz"),
                    "meta_hg38": str(tmpdir / "out" / "hg38" / "chr1_meta.tsv.gz"),
                },
            ):
                result = builder.run(config)

            self.assertEqual(result.chromosomes, ["1"])
            self.assertFalse(stale_r2.exists())
            self.assertFalse(stale_drop.exists())

    def test_builder_run_allows_hg38_source_only_output_paths(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            self._write_dummy_plink_prefix(tmpdir, "panel.1", "1")
            map_hg38 = tmpdir / "hg38.map"
            map_hg38.write_text("chr position Genetic_Map(cM)\n1 100 0.0\n", encoding="utf-8")
            config = ReferencePanelBuildConfig(
                plink_prefix=tmpdir / "panel.@",
                source_genome_build="hg38",
                genetic_map_hg19_sources=None,
                genetic_map_hg38_sources=map_hg38,
                output_dir=tmpdir / "out",
                ld_wind_kb=1.0,
            )
            builder = ref_panel_builder.ReferencePanelBuilder(global_config=GlobalConfig(snp_identifier="chr_pos", genome_build="hg38"))

            def fake_build(prefix, chrom, config, build_state):
                out_root = Path(config.output_dir)
                return {
                    "r2_hg38": str(out_root / "hg38" / f"chr{chrom}_r2.parquet"),
                    "meta_hg38": str(out_root / "hg38" / f"chr{chrom}_meta.tsv.gz"),
                }

            with self.assertLogs("LDSC.ref_panel_builder", level="INFO") as logs:
                with mock.patch.object(
                    ref_panel_builder.ReferencePanelBuilder,
                    "_build_chromosome",
                    side_effect=fake_build,
                ):
                    result = builder.run(config)

            self.assertIn("meta_hg38", result.output_paths)
            self.assertIn("r2_hg38", result.output_paths)
            self.assertNotIn("meta_hg19", result.output_paths)
            self.assertNotIn("r2_hg19", result.output_paths)
            self.assertTrue(any("source-build-only" in message for message in logs.output))

    def test_builder_run_warns_when_only_non_matching_chain_is_available(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            self._write_dummy_plink_prefix(tmpdir, "panel.1", "1")
            map_hg19 = tmpdir / "hg19.map"
            map_hg19.write_text("chr position Genetic_Map(cM)\n1 100 0.0\n", encoding="utf-8")
            config = ReferencePanelBuildConfig(
                plink_prefix=tmpdir / "panel.@",
                source_genome_build="hg19",
                genetic_map_hg19_sources=map_hg19,
                genetic_map_hg38_sources=None,
                liftover_chain_hg38_to_hg19_file=tmpdir / "hg38ToHg19.over.chain",
                output_dir=tmpdir / "out",
                ld_wind_kb=1.0,
            )
            builder = ref_panel_builder.ReferencePanelBuilder(global_config=GlobalConfig(snp_identifier="chr_pos", genome_build="hg38"))

            with self.assertLogs("LDSC.ref_panel_builder", level="INFO") as logs:
                with mock.patch.object(
                    ref_panel_builder.ReferencePanelBuilder,
                    "_build_chromosome",
                    return_value={"r2_hg19": "r2", "meta_hg19": "m19"},
                ):
                    builder.run(config)

            self.assertTrue(any("ignoring the opposite-direction chain" in message for message in logs.output))

    def test_resolve_mappable_snp_positions_skips_liftover_without_matching_chain(self):
        builder = ref_panel_builder.ReferencePanelBuilder(global_config=GlobalConfig(snp_identifier="chr_pos", genome_build="hg38"))
        build_state = ref_panel_builder._BuildState(
            genetic_map_hg19=pd.DataFrame({"CHR": ["1"], "POS": [100], "CM": [0.0]}),
            genetic_map_hg38=None,
            liftover_chain_paths={("hg19", "hg38"): None, ("hg38", "hg19"): "chains/hg38ToHg19.over.chain"},
        )
        chrom_df = pd.DataFrame({"BP": [100]}, index=[4])

        with mock.patch.object(kernel_builder, "LiftOverTranslator") as patched:
            retained, hg19_lookup, hg38_lookup, drop_frame = builder._resolve_mappable_snp_positions(
                build_state=build_state,
                chrom="1",
                source_build="hg19",
                chrom_df=chrom_df,
                keep_snps=np.array([4], dtype=int),
            )

        patched.assert_not_called()
        self.assertEqual(retained.tolist(), [4])
        self.assertEqual(hg19_lookup, {4: 100})
        self.assertEqual(hg38_lookup, {})
        self.assertEqual(len(drop_frame), 0)

    def test_resolve_mappable_snp_positions_keeps_hg38_source_positions_without_chain(self):
        builder = ref_panel_builder.ReferencePanelBuilder(global_config=GlobalConfig(snp_identifier="chr_pos", genome_build="hg38"))
        build_state = ref_panel_builder._BuildState(
            genetic_map_hg19=None,
            genetic_map_hg38=pd.DataFrame({"CHR": ["1"], "POS": [110], "CM": [0.0]}),
            liftover_chain_paths={("hg19", "hg38"): "chains/hg19ToHg38.over.chain", ("hg38", "hg19"): None},
        )
        chrom_df = pd.DataFrame({"BP": [110]}, index=[7])

        retained, hg19_lookup, hg38_lookup, drop_frame = builder._resolve_mappable_snp_positions(
            build_state=build_state,
            chrom="1",
            source_build="hg38",
            chrom_df=chrom_df,
            keep_snps=np.array([7], dtype=int),
        )

        self.assertEqual(retained.tolist(), [7])
        self.assertEqual(hg19_lookup, {})
        self.assertEqual(hg38_lookup, {7: 110})
        self.assertEqual(len(drop_frame), 0)

    def test_validate_emitted_chr_pos_rejects_liftover_collisions(self):
        metadata = pd.DataFrame(
            {
                "CHR": ["1", "1", "1"],
                "SNP": ["1:6205309:A:C", "1:6205308:A:G", "1:7000000:T:C"],
                "POS": [6205309, 6205308, 7000000],
            }
        )

        with self.assertRaisesRegex(ValueError, "hg19.*1:6265368.*1:6205309:A:C.*1:6205308:A:G"):
            ref_panel_builder._validate_emitted_build_chr_pos_uniqueness(
                metadata=metadata,
                positions=np.array([6265368, 6265368, 7060000], dtype=np.int64),
                genome_build="hg19",
                chrom="1",
                snp_identifier="chr_pos",
            )

    def test_validate_emitted_chr_pos_allows_duplicate_positions_in_rsid_mode(self):
        metadata = pd.DataFrame(
            {
                "CHR": ["1", "1"],
                "SNP": ["rs1", "rs2"],
                "POS": [10, 11],
            }
        )

        ref_panel_builder._validate_emitted_build_chr_pos_uniqueness(
            metadata=metadata,
            positions=np.array([20, 20], dtype=np.int64),
            genome_build="hg19",
            chrom="1",
            snp_identifier="rsid",
        )

    def test_liftover_target_collision_removes_variant_from_all_emitted_builds(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            prefix = self._write_plink_prefix_rows(
                tmpdir,
                "panel.1",
                [("1", "rs1", 100), ("1", "rs2", 200), ("1", "rs3", 300)],
            )
            config = ReferencePanelBuildConfig(
                plink_prefix=tmpdir / "panel.@",
                source_genome_build="hg19",
                liftover_chain_hg19_to_hg38_file=tmpdir / "hg19ToHg38.over.chain",
                output_dir=tmpdir / "out",
                ld_wind_kb=1.0,
            )
            build_state = ref_panel_builder._BuildState(
                genetic_map_hg19=None,
                genetic_map_hg38=None,
                liftover_chain_paths={("hg19", "hg38"): str(tmpdir / "hg19ToHg38.over.chain")},
            )
            builder = ref_panel_builder.ReferencePanelBuilder(
                global_config=GlobalConfig(snp_identifier="chr_pos", genome_build="hg19")
            )

            def fake_mappable(*, keep_snps, **_kwargs):
                keep = np.asarray(keep_snps, dtype=int)
                return keep, {0: 100, 1: 200, 2: 300}, {0: 1000, 1: 1000, 2: 3000}, empty_identity_drop_frame()

            bed_instances = []
            for _ in range(2):
                bed = mock.Mock()
                bed.kept_snps = [2]
                bed.maf = np.array([0.3])
                bed.m = 1
                bed.n = 1
                bed.nextSNPs = lambda _width: np.zeros((1, 1))
                bed_instances.append(bed)

            with mock.patch.object(builder, "_resolve_mappable_snp_positions", side_effect=fake_mappable), mock.patch(
                "ldsc.ref_panel_builder.kernel_builder.write_r2_parquet"
            ), mock.patch(
                "ldsc.ref_panel_builder.kernel_builder.write_runtime_metadata_sidecar"
            ), mock.patch(
                "ldsc.ref_panel_builder.kernel_ldscore.PlinkBEDFile",
                side_effect=bed_instances,
            ) as mock_bed:
                with self.assertLogs("LDSC.ref_panel_builder", level="INFO") as logs:
                    builder._build_chromosome(str(prefix), "1", config, build_state)

            self.assertEqual([call.kwargs["keep_snps"] for call in mock_bed.call_args_list], [[2], [2]])
            dropped = self._read_dropped_sidecar(tmpdir / "out" / "diagnostics" / "dropped_snps" / "chr1_dropped.tsv.gz")
            self.assertEqual(dropped["reason"].tolist(), ["target_collision", "target_collision"])
            self.assertTrue(any("synchronized cross-build drops" in message for message in logs.output))

    def test_builder_run_infers_source_genome_build_before_build_state_preparation(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            self._write_plink_prefix_rows(tmpdir, "panel.1", [("1", "rs1", 100), ("1", "rs2", 200)])
            config = ReferencePanelBuildConfig(
                plink_prefix=tmpdir / "panel.@",
                output_dir=tmpdir / "out",
                ld_wind_kb=1.0,
            )
            builder = ref_panel_builder.ReferencePanelBuilder(global_config=GlobalConfig(snp_identifier="rsid"))
            captured = {}

            def fake_build(prefix, chrom, config, build_state):
                captured["source_genome_build"] = config.source_genome_build
                return {
                    "r2_hg19": str(Path(config.output_dir) / "hg19" / f"chr{chrom}_r2.parquet"),
                    "meta_hg19": str(Path(config.output_dir) / "hg19" / f"chr{chrom}_meta.tsv.gz"),
                }

            def fake_resolve_from_frames(hint, snp_identifier, frames, *, context, logger=None):
                captured["resolve_args"] = (hint, snp_identifier, context, logger)
                captured["evidence_frames"] = list(frames)
                return "hg19"

            with mock.patch.object(
                ref_panel_builder,
                "resolve_genome_build",
                side_effect=AssertionError("build-ref-panel should use chunked CHR/POS evidence"),
            ), mock.patch.object(
                ref_panel_builder,
                "resolve_genome_build_from_chr_pos_frames",
                side_effect=fake_resolve_from_frames,
                create=True,
            ) as patched_resolve:
                with mock.patch.object(ref_panel_builder.ReferencePanelBuilder, "_build_chromosome", side_effect=fake_build):
                    builder.run(config)

        self.assertEqual(captured["source_genome_build"], "hg19")
        self.assertEqual(captured["resolve_args"][:3], ("auto", "chr_pos", "build-ref-panel PLINK .bim"))
        self.assertIs(captured["resolve_args"][3], ref_panel_builder.LOGGER)
        self.assertEqual(len(captured["evidence_frames"]), 1)
        self.assertEqual(captured["evidence_frames"][0]["POS"].tolist(), [100, 200])
        self.assertEqual(patched_resolve.call_count, 1)

    def test_prepare_build_state_reads_source_build_specific_restriction_column(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            restriction = tmpdir / "restrict.tsv"
            restriction.write_text("CHR\thg19_POS\thg38_POS\n1\t100\t110\n", encoding="utf-8")
            config = ReferencePanelBuildConfig(
                plink_prefix=tmpdir / "panel.@",
                source_genome_build="hg19",
                output_dir=tmpdir / "out",
                ld_wind_kb=1.0,
                ref_panel_snps_file=restriction,
            )
            builder = ref_panel_builder.ReferencePanelBuilder(
                global_config=GlobalConfig(snp_identifier="chr_pos", genome_build="hg38")
            )

            with self.assertLogs("LDSC.ref_panel_builder", level="INFO") as logs:
                build_state = builder._prepare_build_state(config)

        self.assertEqual(build_state.restriction_values, {"1:100"})
        self.assertTrue(any("source genome build 'hg19'" in message for message in logs.output))

    def test_prepare_build_state_infers_generic_restriction_pos_build_and_accepts_source_match(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            restriction = tmpdir / "restrict.tsv"
            restriction.write_text(
                "CHR\tPOS\tSNP\n"
                "1\t100\trs1\n"
                "chrUn\t200\tbad_chr\n"
                "1\tbad\tbad_pos\n",
                encoding="utf-8",
            )
            config = ReferencePanelBuildConfig(
                plink_prefix=tmpdir / "panel.@",
                source_genome_build="hg19",
                output_dir=tmpdir / "out",
                ld_wind_kb=1.0,
                ref_panel_snps_file=restriction,
            )
            builder = ref_panel_builder.ReferencePanelBuilder(
                global_config=GlobalConfig(snp_identifier="chr_pos")
            )

            with mock.patch.object(ref_panel_builder, "resolve_genome_build", return_value="hg19") as patched_resolve:
                with self.assertLogs("LDSC.ref_panel_builder", level="WARNING") as logs:
                    build_state = builder._prepare_build_state(config)

        self.assertEqual(build_state.restriction_values, {"1:100"})
        self.assertEqual(patched_resolve.call_args.args[0], "auto")
        self.assertEqual(patched_resolve.call_args.args[2]["POS"].tolist(), [100])
        self.assertTrue(any("Dropped 2 SNPs with invalid or missing CHR/POS" in message for message in logs.output))

    def test_prepare_build_state_rejects_generic_restriction_pos_build_mismatch(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            restriction = tmpdir / "restrict.tsv"
            restriction.write_text("CHR\tPOS\n1\t110\n", encoding="utf-8")
            config = ReferencePanelBuildConfig(
                plink_prefix=tmpdir / "panel.@",
                source_genome_build="hg19",
                output_dir=tmpdir / "out",
                ld_wind_kb=1.0,
                ref_panel_snps_file=restriction,
            )
            builder = ref_panel_builder.ReferencePanelBuilder(
                global_config=GlobalConfig(snp_identifier="chr_pos")
            )

            with mock.patch.object(ref_panel_builder, "resolve_genome_build", return_value="hg38"):
                with self.assertRaisesRegex(ValueError, "restriction.*hg38.*source.*hg19"):
                    builder._prepare_build_state(config)

    def test_chr_pos_restriction_filters_before_liftover(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            prefix = self._write_plink_prefix_rows(tmpdir, "panel.1", [("1", "rs1", 100), ("1", "rs2", 200)])
            config = ReferencePanelBuildConfig(
                plink_prefix=tmpdir / "panel.@",
                source_genome_build="hg19",
                output_dir=tmpdir / "out",
                ld_wind_kb=1.0,
            )
            build_state = ref_panel_builder._BuildState(
                genetic_map_hg19=None,
                genetic_map_hg38=None,
                liftover_chain_paths={("hg19", "hg38"): "chain.over"},
                restriction_mode="chr_pos",
                restriction_values={"1:200"},
            )
            builder = ref_panel_builder.ReferencePanelBuilder(global_config=GlobalConfig(snp_identifier="chr_pos"))
            captured = {}

            def stop_after_restriction(*, keep_snps, **_kwargs):
                captured["keep_snps"] = list(map(int, keep_snps))
                raise RuntimeError("stop after restriction")

            with mock.patch.object(builder, "_resolve_mappable_snp_positions", side_effect=stop_after_restriction):
                with self.assertRaisesRegex(RuntimeError, "stop after restriction"):
                    builder._build_chromosome(str(prefix), "1", config, build_state)

        self.assertEqual(captured["keep_snps"], [1])

    def test_source_duplicates_are_dropped_by_identity_cleanup_before_liftover_mapping(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            prefix = self._write_plink_prefix_rows(
                tmpdir,
                "panel.1",
                [("1", "rs1", 100), ("1", "rs2", 100), ("1", "rs3", 300)],
            )
            config = ReferencePanelBuildConfig(
                plink_prefix=tmpdir / "panel.@",
                source_genome_build="hg19",
                output_dir=tmpdir / "out",
                ld_wind_kb=1.0,
            )
            build_state = ref_panel_builder._BuildState(
                genetic_map_hg19=None,
                genetic_map_hg38=None,
                liftover_chain_paths={("hg19", "hg38"): "chain.over"},
            )
            builder = ref_panel_builder.ReferencePanelBuilder(global_config=GlobalConfig(snp_identifier="chr_pos"))
            captured = {}

            def stop_before_mapping(*, keep_snps, **_kwargs):
                captured["keep_snps"] = list(map(int, keep_snps))
                raise RuntimeError("stop before mapping")

            with mock.patch.object(builder, "_resolve_mappable_snp_positions", side_effect=stop_before_mapping):
                with self.assertRaisesRegex(RuntimeError, "stop before mapping"):
                    builder._build_chromosome(str(prefix), "1", config, build_state)

        self.assertEqual(captured["keep_snps"], [2])

    def test_rsid_restriction_filters_before_liftover_and_ignores_global_build(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            prefix = self._write_plink_prefix_rows(tmpdir, "panel.1", [("1", "rs1", 100), ("1", "rs2", 200)])
            config = ReferencePanelBuildConfig(
                plink_prefix=tmpdir / "panel.@",
                source_genome_build="hg19",
                output_dir=tmpdir / "out",
                ld_wind_kb=1.0,
            )
            build_state = ref_panel_builder._BuildState(
                genetic_map_hg19=None,
                genetic_map_hg38=None,
                liftover_chain_paths={("hg19", "hg38"): "chain.over"},
                restriction_mode="rsid",
                restriction_values={"rs2"},
            )
            builder = ref_panel_builder.ReferencePanelBuilder(
                global_config=GlobalConfig(snp_identifier="rsid", genome_build="hg38")
            )
            captured = {}

            def stop_after_restriction(*, keep_snps, **_kwargs):
                captured["keep_snps"] = list(map(int, keep_snps))
                raise RuntimeError("stop after restriction")

            with mock.patch.object(builder, "_resolve_mappable_snp_positions", side_effect=stop_after_restriction):
                with self.assertRaisesRegex(RuntimeError, "stop after restriction"):
                    builder._build_chromosome(str(prefix), "1", config, build_state)

        self.assertEqual(captured["keep_snps"], [1])

    def test_build_restriction_mask_uses_identity_keys_when_restriction_is_allele_aware(self):
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

        mask = kernel_builder.build_restriction_mask(metadata, restriction, "rsid_allele_aware")

        self.assertEqual(mask.tolist(), [True, False, True])

    def test_allele_aware_restriction_filters_unrelated_invalid_panel_alleles_before_cleanup(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            prefix = self._write_plink_prefix_rows(
                tmpdir,
                "panel.1",
                [
                    ("1", "rs_valid", 100, "A", "C"),
                    ("1", "rs_invalid_unrelated", 200, "A", "N"),
                ],
            )
            restriction = RestrictionIdentityKeys(
                keys={"1:100:A:C"},
                match_kind="identity",
                dropped=empty_identity_drop_frame(),
                n_input_rows=1,
                n_retained_keys=1,
            )
            config = ReferencePanelBuildConfig(
                plink_prefix=tmpdir / "panel.@",
                source_genome_build="hg19",
                output_dir=tmpdir / "out",
                ld_wind_kb=1.0,
            )
            build_state = ref_panel_builder._BuildState(
                genetic_map_hg19=None,
                genetic_map_hg38=None,
                liftover_chain_paths={},
                restriction_mode="chr_pos_allele_aware",
                restriction_values=restriction.keys,
                restriction_keys=restriction,
            )
            builder = ref_panel_builder.ReferencePanelBuilder(
                global_config=GlobalConfig(snp_identifier="chr_pos_allele_aware")
            )
            captured = {}

            def stop_after_cleanup(*, keep_snps, **_kwargs):
                captured["keep_snps"] = list(map(int, keep_snps))
                raise RuntimeError("stop after cleanup")

            with mock.patch.object(builder, "_resolve_mappable_snp_positions", side_effect=stop_after_cleanup):
                with self.assertRaisesRegex(RuntimeError, "stop after cleanup"):
                    builder._build_chromosome(str(prefix), "1", config, build_state)

        self.assertEqual(captured["keep_snps"], [0])

    def test_allele_aware_restriction_selected_invalid_panel_allele_reaches_cleanup_sidecar(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            prefix = self._write_plink_prefix_rows(
                tmpdir,
                "panel.1",
                [
                    ("1", "rs_valid", 100, "A", "C"),
                    ("1", "rs_bad", 200, "A", "N"),
                ],
            )
            restriction = RestrictionIdentityKeys(
                keys={"1:100:A:C", "1:200:A:C"},
                match_kind="identity",
                dropped=empty_identity_drop_frame(),
                n_input_rows=2,
                n_retained_keys=2,
            )
            config = ReferencePanelBuildConfig(
                plink_prefix=tmpdir / "panel.@",
                source_genome_build="hg19",
                output_dir=tmpdir / "out",
                ld_wind_kb=1.0,
            )
            build_state = ref_panel_builder._BuildState(
                genetic_map_hg19=None,
                genetic_map_hg38=None,
                liftover_chain_paths={},
                restriction_mode="chr_pos_allele_aware",
                restriction_values=restriction.keys,
                restriction_keys=restriction,
            )
            builder = ref_panel_builder.ReferencePanelBuilder(
                global_config=GlobalConfig(snp_identifier="chr_pos_allele_aware")
            )

            def fake_resolve(*, chrom_df, keep_snps, **_kwargs):
                keep = np.asarray(keep_snps, dtype=int)
                return (
                    keep,
                    {int(idx): int(chrom_df.loc[idx, "BP"]) for idx in keep},
                    {},
                    ref_panel_builder._empty_unified_drop_frame(),
                )

            with mock.patch.object(
                builder,
                "_resolve_mappable_snp_positions",
                side_effect=fake_resolve,
            ), mock.patch(
                "ldsc.ref_panel_builder.kernel_builder.write_r2_parquet"
            ), mock.patch(
                "ldsc.ref_panel_builder.kernel_builder.write_runtime_metadata_sidecar"
            ), mock.patch(
                "ldsc.ref_panel_builder.kernel_ldscore.PlinkBEDFile"
            ) as mock_bed:
                mock_bed.return_value.kept_snps = [0]
                mock_bed.return_value.maf = np.array([0.3])
                mock_bed.return_value.m = 1
                mock_bed.return_value.n = 1
                mock_bed.return_value.nextSNPs = lambda: iter([np.zeros((1, 1))])

                result = builder._build_chromosome(str(prefix), "1", config, build_state)

            dropped = self._read_dropped_sidecar(tmpdir / "out" / "diagnostics" / "dropped_snps" / "chr1_dropped.tsv.gz")
            self.assertIsNotNone(result)
            self.assertEqual(dropped["SNP"].tolist(), ["rs_bad"])
            self.assertEqual(dropped["reason"].tolist(), ["invalid_allele"])
            self.assertEqual(dropped["stage"].tolist(), ["ref_panel_source_identity_cleanup"])
            self.assertEqual(mock_bed.call_args.kwargs["keep_snps"], [0])

    def test_base_chr_pos_restriction_ignores_single_allele_like_column(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "restrict.tsv"
            path.write_text("CHR\thg19_POS\tREF\n1\t100\tA\n", encoding="utf-8")

            restriction = ref_panel_builder._read_ref_panel_snp_restriction(
                path,
                "chr_pos",
                source_genome_build="hg19",
            )

            self.assertEqual(restriction.match_kind, "base")
            self.assertEqual(restriction.keys, {"1:100"})

    def test_builder_run_rejects_duplicate_chromosome_across_inputs(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            self._write_dummy_plink_prefix(tmpdir, "panel_a", "1")
            self._write_dummy_plink_prefix(tmpdir, "panel_b", "1")
            map_hg19 = tmpdir / "hg19.map"
            map_hg38 = tmpdir / "hg38.map"
            map_hg19.write_text("chr position Genetic_Map(cM)\n1 100 0.0\n", encoding="utf-8")
            map_hg38.write_text("chr position Genetic_Map(cM)\n1 100 0.0\n", encoding="utf-8")
            config = ReferencePanelBuildConfig(
                plink_prefix=str(tmpdir / "panel_*"),
                source_genome_build="hg19",
                genetic_map_hg19_sources=map_hg19,
                genetic_map_hg38_sources=map_hg38,
                liftover_chain_hg19_to_hg38_file=tmpdir / "hg19ToHg38.over.chain",
                output_dir=tmpdir / "out",
                ld_wind_kb=1.0,
            )
            builder = ref_panel_builder.ReferencePanelBuilder(global_config=GlobalConfig(snp_identifier="chr_pos", genome_build="hg38"))

            with mock.patch.object(
                ref_panel_builder.ReferencePanelBuilder,
                "_build_chromosome",
                return_value={"r2_hg19": "r2-19", "r2_hg38": "r2-38", "meta_hg19": "m19", "meta_hg38": "m38"},
            ):
                with self.assertRaises(ValueError):
                    builder.run(config)

    def test_map_positions_passes_explicit_chain_path_to_translator(self):
        builder = ref_panel_builder.ReferencePanelBuilder(global_config=GlobalConfig(snp_identifier="chr_pos", genome_build="hg38"))
        build_state = ref_panel_builder._BuildState(
            genetic_map_hg19=pd.DataFrame({"CHR": ["1"], "POS": [100], "CM": [0.0]}),
            genetic_map_hg38=pd.DataFrame({"CHR": ["1"], "POS": [100], "CM": [0.0]}),
            liftover_chain_paths={("hg38", "hg19"): "chains/hg38ToHg19.over.chain"},
        )
        expected = kernel_builder.LiftOverMappingResult(
            translated_positions=np.array([100], dtype=np.int64),
            keep_mask=np.array([True], dtype=bool),
            unmapped_count=0,
            cross_chrom_count=0,
        )

        with mock.patch.object(kernel_builder, "LiftOverTranslator") as patched:
            patched.return_value.map_positions.return_value = expected
            result = builder._map_positions(
                build_state=build_state,
                chrom="1",
                positions=np.array([100], dtype=np.int64),
                source_build="hg38",
                target_build="hg19",
            )

        self.assertIs(result, expected)
        patched.assert_called_once_with(
            source_build="hg38",
            target_build="hg19",
            chain_path="chains/hg38ToHg19.over.chain",
            chain_flag_hint="--liftover-chain-hg38-to-hg19-file",
            workflow_label="reference-panel liftover",
        )


class RefPanelOutputFamilyTest(unittest.TestCase):
    def test_output_family_includes_existing_owned_artifacts_and_current_run_paths(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            out = Path(tmpdir) / "out"
            existing_r2 = out / "hg19" / "chr1_r2.parquet"
            existing_meta = out / "hg38" / "chr2_meta.tsv.gz"
            existing_drop = out / "diagnostics" / "dropped_snps" / "chr2_dropped.tsv.gz"
            existing_log = out / "diagnostics" / "build-ref-panel.chr2.log"
            produced = out / "hg19" / "chr6_r2.parquet"
            existing_r2.parent.mkdir(parents=True)
            existing_meta.parent.mkdir(parents=True)
            existing_drop.parent.mkdir(parents=True)
            existing_r2.write_text("r2\n", encoding="utf-8")
            existing_meta.write_text("meta\n", encoding="utf-8")
            existing_drop.write_text("drop\n", encoding="utf-8")
            existing_log.write_text("log\n", encoding="utf-8")

            family = ref_panel_builder._ref_panel_output_family(out, [produced])

            self.assertIn(existing_r2, family)
            self.assertIn(existing_meta, family)
            self.assertIn(existing_drop, family)
            self.assertIn(existing_log, family)
            self.assertIn(produced, family)

    def test_output_family_ignores_unknown_subdirectories(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            out = Path(tmpdir) / "out"
            scratch = out / "scratch" / "chr1_r2.parquet"
            scratch.parent.mkdir(parents=True)
            scratch.write_text("unknown\n", encoding="utf-8")

            family = ref_panel_builder._ref_panel_output_family(out, [])

            self.assertNotIn(scratch, family)


class ReferencePanelBuilderSourceOnlySmokeTest(unittest.TestCase):
    def test_hm3_chr22_subset_builds_source_only_without_liftover_chain_or_map(self):
        prefix = MINIMAL_EXTERNAL_FIXTURES / "plink" / "hm3_chr22_subset"
        if not (Path(str(prefix) + ".bed").exists() and Path(str(prefix) + ".bim").exists() and Path(str(prefix) + ".fam").exists()):
            self.skipTest("minimal chr22 PLINK fixture is unavailable; run tests/fixtures/generate_minimal_external_resources.py")

        with tempfile.TemporaryDirectory() as tmpdir:
            build_result = ref_panel_builder.run_build_ref_panel(
                plink_prefix=str(prefix),
                source_genome_build="hg38",
                genetic_map_hg19_sources=None,
                genetic_map_hg38_sources=None,
                output_dir=str(Path(tmpdir) / "panel"),
                ld_wind_snps=10,
                ld_wind_kb=None,
                snp_batch_size=64,
            )

            self.assertEqual(build_result.chromosomes, ["22"])
            self.assertIn("meta_hg38", build_result.output_paths)
            self.assertIn("r2_hg38", build_result.output_paths)
            self.assertNotIn("meta_hg19", build_result.output_paths)
            self.assertNotIn("r2_hg19", build_result.output_paths)
            self.assertNotIn("ann", build_result.output_paths)
            self.assertNotIn("ld", build_result.output_paths)
            self.assertFalse((Path(tmpdir) / "panel" / "parquet").exists())
            self.assertTrue(Path(build_result.output_paths["r2_hg38"][0]).exists())
            self.assertTrue(Path(build_result.output_paths["meta_hg38"][0]).exists())

            meta_hg38 = pd.read_csv(build_result.output_paths["meta_hg38"][0], sep="\t", comment="#")
            self.assertTrue(meta_hg38["CM"].isna().all())


@unittest.skipUnless(
    _HAS_BITARRAY and _HAS_PYARROW and _HAS_PYLIFTOVER,
    "bitarray, pyarrow, and pyliftover are required for the smoke/parity builder test",
)
class ReferencePanelBuilderParityTest(unittest.TestCase):
    def test_hm3_chr22_subset_writes_target_r2_and_metadata_with_missing_cm_when_target_map_is_absent(self):
        resources = _find_resources_root()
        if resources is None:
            self.skipTest("resources directory is not available from this workspace")

        prefix = MINIMAL_EXTERNAL_FIXTURES / "plink" / "hm3_chr22_subset"
        if not (Path(str(prefix) + ".bed").exists() and Path(str(prefix) + ".bim").exists() and Path(str(prefix) + ".fam").exists()):
            self.skipTest("minimal chr22 PLINK fixture is unavailable; run tests/fixtures/generate_minimal_external_resources.py")

        map_hg38 = MINIMAL_EXTERNAL_FIXTURES / "genetic_maps" / "genetic_map_hg38_chr22_subset.txt"
        if not map_hg38.exists():
            self.skipTest("minimal hg38 genetic-map fixture is unavailable; run tests/fixtures/generate_minimal_external_resources.py")

        with tempfile.TemporaryDirectory() as tmpdir:
            build_result = ref_panel_builder.run_build_ref_panel(
                plink_prefix=str(prefix),
                source_genome_build="hg38",
                genetic_map_hg19_sources=None,
                genetic_map_hg38_sources=str(map_hg38),
                liftover_chain_hg38_to_hg19_file=str(resources / "liftover" / "hg38ToHg19.over.chain"),
                output_dir=str(Path(tmpdir) / "panel"),
                ld_wind_snps=10,
                ld_wind_kb=None,
                snp_batch_size=64,
            )

            self.assertIn("meta_hg19", build_result.output_paths)
            self.assertIn("meta_hg38", build_result.output_paths)
            self.assertIn("r2_hg19", build_result.output_paths)
            self.assertIn("r2_hg38", build_result.output_paths)
            meta_hg19 = pd.read_csv(build_result.output_paths["meta_hg19"][0], sep="\t", comment="#")
            meta_hg38 = pd.read_csv(build_result.output_paths["meta_hg38"][0], sep="\t", comment="#")
            self.assertTrue(meta_hg19["CM"].isna().all())
            self.assertFalse(meta_hg38["CM"].isna().any())
            self.assertTrue(Path(build_result.output_paths["r2_hg19"][0]).exists())
            self.assertTrue(Path(build_result.output_paths["r2_hg38"][0]).exists())

    def test_cm_window_requires_target_map_when_liftover_emits_target_build(self):
        resources = _find_resources_root()
        if resources is None:
            self.skipTest("resources directory is not available from this workspace")

        prefix = MINIMAL_EXTERNAL_FIXTURES / "plink" / "hm3_chr22_subset"
        if not (Path(str(prefix) + ".bed").exists() and Path(str(prefix) + ".bim").exists() and Path(str(prefix) + ".fam").exists()):
            self.skipTest("minimal chr22 PLINK fixture is unavailable; run tests/fixtures/generate_minimal_external_resources.py")

        map_hg38 = MINIMAL_EXTERNAL_FIXTURES / "genetic_maps" / "genetic_map_hg38_chr22_subset.txt"
        if not map_hg38.exists():
            self.skipTest("minimal hg38 genetic-map fixture is unavailable; run tests/fixtures/generate_minimal_external_resources.py")

        with tempfile.TemporaryDirectory() as tmpdir:
            with self.assertRaisesRegex(ValueError, "hg19 genetic map.*required"):
                ref_panel_builder.run_build_ref_panel(
                    plink_prefix=str(prefix),
                    source_genome_build="hg38",
                    genetic_map_hg19_sources=None,
                    genetic_map_hg38_sources=str(map_hg38),
                    liftover_chain_hg38_to_hg19_file=str(resources / "liftover" / "hg38ToHg19.over.chain"),
                    output_dir=str(Path(tmpdir) / "panel"),
                    ld_wind_cm=1.0,
                    ld_wind_kb=None,
                    snp_batch_size=64,
                )

    def test_hm3_chr22_subset_runs_direct_and_parquet_ldscore_paths(self):
        resources = _find_resources_root()
        if resources is None:
            self.skipTest("resources directory is not available from this workspace")

        prefix = MINIMAL_EXTERNAL_FIXTURES / "plink" / "hm3_chr22_subset"
        if not (Path(str(prefix) + ".bed").exists() and Path(str(prefix) + ".bim").exists() and Path(str(prefix) + ".fam").exists()):
            self.skipTest("minimal chr22 PLINK fixture is unavailable; run tests/fixtures/generate_minimal_external_resources.py")

        map_hg19 = MINIMAL_EXTERNAL_FIXTURES / "genetic_maps" / "genetic_map_hg19_chr22_subset.txt"
        map_hg38 = MINIMAL_EXTERNAL_FIXTURES / "genetic_maps" / "genetic_map_hg38_chr22_subset.txt"
        if not (map_hg19.exists() and map_hg38.exists()):
            self.skipTest("minimal genetic-map fixtures are unavailable; run tests/fixtures/generate_minimal_external_resources.py")

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            build_result = ref_panel_builder.run_build_ref_panel(
                plink_prefix=str(prefix),
                source_genome_build="hg38",
                genetic_map_hg19_sources=str(map_hg19),
                genetic_map_hg38_sources=str(map_hg38),
                liftover_chain_hg38_to_hg19_file=str(resources / "liftover" / "hg38ToHg19.over.chain"),
                output_dir=str(tmpdir / "panel"),
                ld_wind_snps=10,
                ld_wind_kb=None,
                snp_batch_size=64,
            )

            self.assertEqual(build_result.chromosomes, ["22"])
            r2_path = build_result.output_paths["r2_hg38"][0]
            meta_hg38_path = build_result.output_paths["meta_hg38"][0]
            self.assertTrue(Path(r2_path).exists())
            self.assertTrue(Path(meta_hg38_path).exists())

            with gzip.open(meta_hg38_path, "rt", encoding="utf-8") as handle:
                meta = pd.read_csv(handle, sep="\t", comment="#")
            baseline = tmpdir / "baseline.annot.gz"
            baseline_df = meta.loc[:, ["CHR", "POS", "SNP", "CM"]].rename(columns={"POS": "BP"})
            baseline_df["base"] = 1.0
            with gzip.open(baseline, "wt", encoding="utf-8") as handle:
                baseline_df.to_csv(handle, sep="\t", index=False)

            set_global_config(GlobalConfig(snp_identifier="chr_pos_allele_aware", genome_build="hg38"))
            try:
                direct = ldscore_calculator.run_ldscore(
                    output_dir=str(tmpdir / "direct"),
                    baseline_annot_sources=str(baseline),
                    plink_prefix=str(prefix),
                    ld_wind_snps=10,
                    snp_batch_size=64,
                )
                parquet = ldscore_calculator.run_ldscore(
                    output_dir=str(tmpdir / "parquet"),
                    baseline_annot_sources=str(baseline),
                    r2_dir=str(tmpdir / "panel" / "hg38"),
                    r2_bias_mode="raw",
                    r2_sample_size=3202,
                    ld_wind_snps=10,
                    snp_batch_size=64,
                )
            finally:
                reset_global_config()

            self.assertEqual(
                direct.baseline_table["SNP"].tolist(),
                parquet.baseline_table["SNP"].tolist(),
            )
            self.assertEqual(len(direct.baseline_table), 32)
            self.assertEqual(len(parquet.baseline_table), 32)
            self.assertTrue(np.all(np.isfinite(direct.baseline_table["base"].to_numpy(dtype=float))))
            self.assertTrue(np.all(np.isfinite(parquet.baseline_table["base"].to_numpy(dtype=float))))
            self.assertGreater(float(direct.baseline_table["base"].max()), 1.0)
            self.assertGreater(float(parquet.baseline_table["base"].max()), 1.0)
