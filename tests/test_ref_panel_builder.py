import gzip
import importlib.util
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

    def test_explicit_chain_path_is_required_for_cross_build_translation(self):
        with self.assertRaises(ValueError) as exc:
            kernel_builder.LiftOverTranslator(source_build="hg38", target_build="hg19", chain_path=None)

        self.assertIn("--liftover-chain-hg38-to-hg19", str(exc.exception))


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
            chunk_size=2,
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
            chunk_size=2,
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
            ["chr", "hg19_pos", "hg38_pos", "hg19_Uniq_ID", "hg38_Uniq_ID", "rsID", "MAF", "REF", "ALT"],
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
            ["chr", "hg19_pos", "hg38_pos", "hg19_Uniq_ID", "hg38_Uniq_ID", "rsID", "MAF", "REF", "ALT"],
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
                "REF": ["A", "C"],
                "ALT": ["G", "T"],
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
            ["CHR", "POS_1", "POS_2", "R2", "SNP_1", "SNP_2"],
        )
        self.assertEqual(table.loc[0, "CHR"], "1")
        self.assertEqual(table.loc[0, "POS_1"], 100)
        self.assertEqual(table.loc[0, "POS_2"], 120)
        self.assertAlmostEqual(table.loc[0, "R2"], 0.75, places=4)
        self.assertEqual(table.loc[0, "SNP_1"], "rs1")
        self.assertEqual(table.loc[0, "SNP_2"], "rs2")
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
                "REF": ["A"],
                "ALT": ["G"],
            }
        )

        table = kernel_builder.build_standard_r2_table(
            pair_rows=[],
            reference_snp_table=reference_snp_table,
            genome_build="hg19",
        )

        self.assertEqual(
            table.columns.tolist(),
            ["CHR", "POS_1", "POS_2", "R2", "SNP_1", "SNP_2"],
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
                "REF": ["A", "C"],
                "ALT": ["G", "T"],
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
                "REF": ["A", "C", "G"],
                "ALT": ["G", "T", "A"],
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
                "REF": ["A", "C"],
                "ALT": ["G", "T"],
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
                row_group_size=50_000,
            )
            pf = pq.ParquetFile(str(path))
            meta = pf.schema_arrow.metadata
            self.assertIsNotNone(meta)
            self.assertIn(b"ldsc:sorted_by_build", meta)
            self.assertEqual(meta[b"ldsc:sorted_by_build"].decode("utf-8"), "hg19")
            self.assertIn(b"ldsc:row_group_size", meta)
            self.assertEqual(meta[b"ldsc:row_group_size"].decode("utf-8"), "50000")

    def test_build_runtime_metadata_table_is_build_specific(self):
        metadata = pd.DataFrame(
            {
                "CHR": ["1"],
                "SNP": ["rs1"],
                "MAF": [0.2],
            }
        )

        table = kernel_builder.build_runtime_metadata_table(
            metadata=metadata,
            positions=np.array([120], dtype=np.int64),
            cm_values=np.array([0.5], dtype=float),
        )

        self.assertEqual(table.columns.tolist(), ["CHR", "POS", "SNP", "CM", "MAF"])
        self.assertEqual(table.loc[0, "POS"], 120)
        self.assertEqual(table.loc[0, "CM"], 0.5)

    def test_build_runtime_metadata_table_allows_missing_cm_values(self):
        metadata = pd.DataFrame(
            {
                "CHR": ["1"],
                "SNP": ["rs1"],
                "MAF": [0.2],
            }
        )

        table = kernel_builder.build_runtime_metadata_table(
            metadata=metadata,
            positions=np.array([120], dtype=np.int64),
            cm_values=None,
        )

        self.assertEqual(table.columns.tolist(), ["CHR", "POS", "SNP", "CM", "MAF"])
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


class ReferencePanelBuildConfigFromArgsTest(unittest.TestCase):
    def test_build_parser_defaults_chunk_size_to_128(self):
        parser = ref_panel_builder.build_parser()
        self.assertEqual(parser.get_default("chunk_size"), 128)

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

    def test_config_from_args_uses_registered_identifier_for_restriction_file(self):
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
        self.assertIsNone(build_config.source_genome_build)
        self.assertEqual(global_config.snp_identifier, "chr_pos")
        self.assertIsNone(global_config.genome_build)

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
        self.assertIsNone(global_config.genome_build)

    def test_run_build_ref_panel_uses_registered_identifier_for_restriction_file(self):
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
        self.assertEqual(captured["global_config"].snp_identifier, "chr_pos")
        self.assertIsNone(captured["global_config"].genome_build)

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


class ReferencePanelBuilderWorkflowTest(unittest.TestCase):
    def _write_dummy_plink_prefix(self, root: Path, stem: str, chrom: str):
        prefix = root / stem
        Path(str(prefix) + ".bed").write_bytes(b"")
        Path(str(prefix) + ".fam").write_text("fam iid 0 0 0 -9\n", encoding="utf-8")
        Path(str(prefix) + ".bim").write_text(
            f"{chrom} rs{chrom} 0.0 100 A G\n",
            encoding="utf-8",
        )
        return prefix

    def _write_plink_prefix_rows(self, root: Path, stem: str, rows: list[tuple[str, str, int]]):
        prefix = root / stem
        Path(str(prefix) + ".bed").write_bytes(b"")
        Path(str(prefix) + ".fam").write_text("fam iid 0 0 0 -9\n", encoding="utf-8")
        Path(str(prefix) + ".bim").write_text(
            "".join(f"{chrom} {snp} 0.0 {pos} A G\n" for chrom, snp, pos in rows),
            encoding="utf-8",
        )
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
            retained, hg19_lookup, hg38_lookup = builder._resolve_mappable_snp_positions(
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

    def test_resolve_mappable_snp_positions_keeps_hg38_source_positions_without_chain(self):
        builder = ref_panel_builder.ReferencePanelBuilder(global_config=GlobalConfig(snp_identifier="chr_pos", genome_build="hg38"))
        build_state = ref_panel_builder._BuildState(
            genetic_map_hg19=None,
            genetic_map_hg38=pd.DataFrame({"CHR": ["1"], "POS": [110], "CM": [0.0]}),
            liftover_chain_paths={("hg19", "hg38"): "chains/hg19ToHg38.over.chain", ("hg38", "hg19"): None},
        )
        chrom_df = pd.DataFrame({"BP": [110]}, index=[7])

        retained, hg19_lookup, hg38_lookup = builder._resolve_mappable_snp_positions(
            build_state=build_state,
            chrom="1",
            source_build="hg38",
            chrom_df=chrom_df,
            keep_snps=np.array([7], dtype=int),
        )

        self.assertEqual(retained.tolist(), [7])
        self.assertEqual(hg19_lookup, {})
        self.assertEqual(hg38_lookup, {7: 110})

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

            with mock.patch.object(ref_panel_builder, "resolve_genome_build", return_value="hg19") as patched_resolve:
                with mock.patch.object(ref_panel_builder.ReferencePanelBuilder, "_build_chromosome", side_effect=fake_build):
                    builder.run(config)

        self.assertEqual(captured["source_genome_build"], "hg19")
        self.assertEqual(patched_resolve.call_args.args[0], "auto")
        self.assertEqual(patched_resolve.call_args.args[1], "chr_pos")
        self.assertEqual(patched_resolve.call_args.args[2]["POS"].tolist(), [100, 200])

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
            restriction.write_text("CHR\tPOS\n1\t100\n", encoding="utf-8")
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
                build_state = builder._prepare_build_state(config)

        self.assertEqual(build_state.restriction_values, {"1:100"})
        self.assertEqual(patched_resolve.call_args.args[0], "auto")
        self.assertEqual(patched_resolve.call_args.args[2]["POS"].tolist(), [100])

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
        )


@unittest.skipUnless(
    _HAS_BITARRAY and _HAS_PYARROW,
    "bitarray and pyarrow are required for the source-only builder smoke test",
)
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
                chunk_size=64,
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

            meta_hg38 = pd.read_csv(build_result.output_paths["meta_hg38"][0], sep="\t")
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
                chunk_size=64,
            )

            self.assertIn("meta_hg19", build_result.output_paths)
            self.assertIn("meta_hg38", build_result.output_paths)
            self.assertIn("r2_hg19", build_result.output_paths)
            self.assertIn("r2_hg38", build_result.output_paths)
            meta_hg19 = pd.read_csv(build_result.output_paths["meta_hg19"][0], sep="\t")
            meta_hg38 = pd.read_csv(build_result.output_paths["meta_hg38"][0], sep="\t")
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
                    chunk_size=64,
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
                chunk_size=64,
            )

            self.assertEqual(build_result.chromosomes, ["22"])
            r2_path = build_result.output_paths["r2_hg38"][0]
            meta_hg38_path = build_result.output_paths["meta_hg38"][0]
            self.assertTrue(Path(r2_path).exists())
            self.assertTrue(Path(meta_hg38_path).exists())

            with gzip.open(meta_hg38_path, "rt", encoding="utf-8") as handle:
                meta = pd.read_csv(handle, sep="\t")
            baseline = tmpdir / "baseline.annot.gz"
            baseline_df = meta.loc[:, ["CHR", "POS", "SNP", "CM"]].rename(columns={"POS": "BP"})
            baseline_df["base"] = 1.0
            with gzip.open(baseline, "wt", encoding="utf-8") as handle:
                baseline_df.to_csv(handle, sep="\t", index=False)

            set_global_config(GlobalConfig(snp_identifier="rsid"))
            try:
                direct = ldscore_calculator.run_ldscore(
                    output_dir=str(tmpdir / "direct"),
                    baseline_annot_sources=str(baseline),
                    plink_prefix=str(prefix),
                    ld_wind_snps=10,
                    chunk_size=64,
                )
                parquet = ldscore_calculator.run_ldscore(
                    output_dir=str(tmpdir / "parquet"),
                    baseline_annot_sources=str(baseline),
                    r2_dir=str(tmpdir / "panel" / "hg38"),
                    r2_bias_mode="raw",
                    r2_sample_size=3202,
                    ld_wind_snps=10,
                    chunk_size=64,
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
