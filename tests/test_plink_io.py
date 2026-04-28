from pathlib import Path
import importlib
import sys
import tempfile
import unittest
from unittest import mock

import numpy as np
from numpy.testing import assert_array_almost_equal, assert_array_equal
import pandas as pd

SRC = Path(__file__).resolve().parents[1] / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from ldsc._kernel import formats as ps

try:
    import bitarray as ba
    from ldsc._kernel import ldscore as ld
except ModuleNotFoundError:
    ba = None
    ld = None


PLINK_FIXTURES = Path(__file__).resolve().parent / "fixtures" / "plink"


def _has_module(name: str) -> bool:
    try:
        importlib.import_module(name)
        return True
    except Exception:
        return False


@unittest.skipIf(ld is None or ba is None, "bitarray dependency is not installed")
class LDScoreHelpersTest(unittest.TestCase):
    def test_get_block_lefts(self):
        cases = [
            (np.arange(1, 6), 5, np.zeros(5)),
            (np.arange(1, 6), 0, np.arange(0, 5)),
            ((1, 4, 6, 7, 7, 8), 2, (0, 1, 1, 2, 2, 2)),
        ]
        for coords, max_dist, expected in cases:
            assert_array_equal(ld.getBlockLefts(coords, max_dist), expected)

    def test_block_left_to_right(self):
        cases = [
            ((0, 0, 0, 0, 0), (5, 5, 5, 5, 5)),
            ((0, 1, 2, 3, 4, 5), (1, 2, 3, 4, 5, 6)),
            ((0, 0, 2, 2), (2, 2, 4, 4)),
        ]
        for block_left, expected in cases:
            assert_array_equal(ld.block_left_to_right(block_left), expected)


@unittest.skipIf(ld is None or ba is None, "bitarray dependency is not installed")
class PlinkBedFileTest(unittest.TestCase):
    def setUp(self):
        self.snp_count = 8
        self.sample_count = 5
        self.bim = ps.PlinkBIMFile(str(PLINK_FIXTURES / "plink.bim"))

    def test_bed(self):
        bed = ld.PlinkBEDFile(str(PLINK_FIXTURES / "plink.bed"), self.sample_count, self.bim)
        self.assertEqual(bed.m, 4)
        self.assertEqual(bed.n, self.sample_count)
        self.assertEqual(len(bed.geno), 64)
        expected = np.array([0.6, 0.6, 0.625, 0.625])
        assert_array_almost_equal(bed.freq, expected)

    def test_filter_snps(self):
        bed = ld.PlinkBEDFile(
            str(PLINK_FIXTURES / "plink.bed"),
            self.sample_count,
            self.bim,
            keep_snps=[1, 4],
        )
        self.assertEqual(bed.m, 1)
        self.assertEqual(bed.n, 5)
        self.assertEqual(bed.geno[0:10], ba.bitarray("0001011111"))

    def test_filter_indivs(self):
        bed = ld.PlinkBEDFile(
            str(PLINK_FIXTURES / "plink.bed"),
            self.sample_count,
            self.bim,
            keep_indivs=[0, 1],
        )
        self.assertEqual(bed.m, 2)
        self.assertEqual(bed.n, 2)
        self.assertEqual(bed.geno[0:4], ba.bitarray("0001"))
        self.assertEqual(bed.geno[8:12], ba.bitarray("0001"))

    def test_filter_indivs_and_snps(self):
        bed = ld.PlinkBEDFile(
            str(PLINK_FIXTURES / "plink.bed"),
            self.sample_count,
            self.bim,
            keep_snps=[1, 5],
            keep_indivs=[0, 1],
        )
        self.assertEqual(bed.m, 1)
        self.assertEqual(bed.n, 2)
        self.assertEqual(bed.geno[0:4], ba.bitarray("0001"))

    def test_bad_filename(self):
        with self.assertRaises(ValueError):
            ld.PlinkBEDFile(str(PLINK_FIXTURES / "plink.bim"), 9, self.bim)

    def test_next_snps_errors(self):
        bed = ld.PlinkBEDFile(str(PLINK_FIXTURES / "plink.bed"), self.sample_count, self.bim)
        with self.assertRaises(ValueError):
            bed.nextSNPs(0)
        with self.assertRaises(ValueError):
            bed.nextSNPs(5)

    def test_next_snps(self):
        for width in [1, 2, 3]:
            bed = ld.PlinkBEDFile(
                str(PLINK_FIXTURES / "plink.bed"), self.sample_count, self.bim
            )
            x = bed.nextSNPs(width)
            self.assertEqual(x.shape, (5, width))
            self.assertTrue(np.all(np.abs(np.mean(x, axis=0)) < 0.01))
            self.assertTrue(np.all(np.abs(np.std(x, axis=0) - 1) < 0.01))

    def test_next_snps_maf_ref(self):
        width = 4
        bed = ld.PlinkBEDFile(str(PLINK_FIXTURES / "plink.bed"), self.sample_count, self.bim)
        x = bed.nextSNPs(width)
        bed._currentSNP -= width
        y = bed.nextSNPs(width, minorRef=True)
        assert_array_equal(x, -y)


@unittest.skipIf(ld is None, "ldscore kernel is not available")
@unittest.skipUnless(_has_module("pyarrow"), "pyarrow dependency is not installed")
class RawParquetRuntimeTest(unittest.TestCase):
    def test_raw_parquet_reader_support_hg19(self):
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
                    "BP": [10, 20],
                    "CM": [0.1, 0.2],
                }
            )
            reader = ld.SortedR2BlockReader(
                paths=[str(path)],
                chrom="1",
                metadata=metadata,
                identifier_mode="rsID",
                r2_bias_mode="unbiased",
                r2_sample_size=None,
                genome_build="hg19",
            )
            matrix = reader.within_block_matrix(l_B=0, c=2)
            expected = np.array([[1.0, 0.4], [0.4, 1.0]], dtype=np.float32)
            assert_array_almost_equal(matrix, expected)


@unittest.skipIf(ld is None, "ldscore kernel is not available")
@unittest.skipUnless(_has_module("pyarrow"), "pyarrow dependency is not installed")
class CanonicalParquetRuntimeTest(unittest.TestCase):
    """Tests for the canonical 6-column parquet schema path."""

    def _write_canonical_parquet(self, path: Path, *, genome_build: str = "hg19", n_pairs: int = 5, row_group_size: int = 3) -> None:
        import pyarrow as pa
        import pyarrow.parquet as pq

        pairs = [
            {
                "CHR": "1",
                "POS_1": 100 + i * 10,
                "POS_2": 120 + i * 10,
                "R2": 0.5,
                "SNP_1": f"rs{i + 1}",
                "SNP_2": f"rs{i + 2}",
            }
            for i in range(n_pairs)
        ]
        df = pd.DataFrame(pairs)
        table = pa.Table.from_pandas(df, preserve_index=False)
        meta = {
            b"ldsc:sorted_by_build": genome_build.encode("utf-8"),
            b"ldsc:row_group_size": str(row_group_size).encode("utf-8"),
        }
        enriched = table.schema.with_metadata({**(table.schema.metadata or {}), **meta})
        with pq.ParquetWriter(str(path), enriched) as writer:
            writer.write_table(table.cast(enriched), row_group_size=row_group_size)

    def test_canonical_init_builds_rg_bounds(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "chr1.parquet"
            self._write_canonical_parquet(path, n_pairs=6, row_group_size=2)
            metadata = pd.DataFrame(
                {
                    "CHR": ["1"] * 7,
                    "SNP": [f"rs{i + 1}" for i in range(7)],
                    "BP": [100 + i * 10 for i in range(7)],
                    "CM": [0.0] * 7,
                }
            )
            reader = ld.SortedR2BlockReader(
                paths=[str(path)],
                chrom="1",
                metadata=metadata,
                identifier_mode="rsID",
                r2_bias_mode="unbiased",
                r2_sample_size=None,
                genome_build="hg19",
            )
            self.assertEqual(reader._runtime_layout, "canonical")
            self.assertTrue(hasattr(reader, "_rg_bounds"))
            self.assertGreater(len(reader._rg_bounds), 1)
            for mn, mx, idx in reader._rg_bounds:
                self.assertLessEqual(mn, mx)
                self.assertIsInstance(idx, int)

    def test_canonical_init_excludes_row_groups_without_min_max_stats(self):
        import pyarrow.parquet as pq

        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "chr1.parquet"
            self._write_canonical_parquet(path, n_pairs=1, row_group_size=1)

            real_pf_cls = pq.ParquetFile

            class FakeColumn:
                statistics = None

            class FakeRowGroup:
                def column(self, index):
                    return FakeColumn()

            class FakeMeta:
                num_rows = 1
                num_row_groups = 1

                def row_group(self, index):
                    return FakeRowGroup()

            class FakePF:
                def __init__(self, parquet_path):
                    self._real = real_pf_cls(parquet_path)
                    self.metadata = FakeMeta()
                    self.schema_arrow = self._real.schema_arrow

                def read_row_groups(self, *args, **kwargs):
                    return self._real.read_row_groups(*args, **kwargs)

                def read_row_group(self, *args, **kwargs):
                    return self._real.read_row_group(*args, **kwargs)

            metadata = pd.DataFrame({"CHR": ["1", "1"], "SNP": ["rs1", "rs2"], "BP": [100, 120], "CM": [0.0, 0.0]})
            with mock.patch("pyarrow.parquet.ParquetFile", FakePF):
                reader = ld.SortedR2BlockReader(
                    paths=[str(path)],
                    chrom="1",
                    metadata=metadata,
                    identifier_mode="rsID",
                    r2_bias_mode="unbiased",
                    r2_sample_size=None,
                    genome_build="hg19",
                )
            self.assertEqual(reader._rg_bounds, [])

    def test_canonical_init_raises_on_build_mismatch(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "chr1.parquet"
            self._write_canonical_parquet(path, genome_build="hg19")
            metadata = pd.DataFrame(
                {"CHR": ["1", "1"], "SNP": ["rs1", "rs2"], "BP": [100, 120], "CM": [0.0, 0.0]}
            )
            with self.assertRaises(ValueError) as ctx:
                ld.SortedR2BlockReader(
                    paths=[str(path)],
                    chrom="1",
                    metadata=metadata,
                    identifier_mode="rsID",
                    r2_bias_mode="unbiased",
                    r2_sample_size=None,
                    genome_build="hg38",
                )
            self.assertIn("hg19", str(ctx.exception))
            self.assertIn("hg38", str(ctx.exception))

    def test_canonical_init_warns_on_coarse_row_groups(self):
        import pyarrow as pa
        import pyarrow.parquet as pq

        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "chr1.parquet"
            df = pd.DataFrame(
                {"CHR": ["1"], "POS_1": [100], "POS_2": [120], "R2": [0.5], "SNP_1": ["rs1"], "SNP_2": ["rs2"]}
            )
            table = pa.Table.from_pandas(df, preserve_index=False)
            meta = {b"ldsc:sorted_by_build": b"hg19", b"ldsc:row_group_size": b"50000"}
            enriched = table.schema.with_metadata({**(table.schema.metadata or {}), **meta})
            with pq.ParquetWriter(str(path), enriched) as writer:
                writer.write_table(table.cast(enriched))

            real_pf_cls = pq.ParquetFile

            class FakeMeta:
                num_rows = 600_000
                num_row_groups = 1

                def row_group(self, index):
                    return real_pf_cls(str(path)).metadata.row_group(index)

            class FakePF:
                def __init__(self, parquet_path):
                    self._real = real_pf_cls(parquet_path)
                    self.metadata = FakeMeta()
                    self.schema_arrow = self._real.schema_arrow

                def read_row_groups(self, *args, **kwargs):
                    return self._real.read_row_groups(*args, **kwargs)

                def read_row_group(self, *args, **kwargs):
                    return self._real.read_row_group(*args, **kwargs)

            metadata = pd.DataFrame({"CHR": ["1", "1"], "SNP": ["rs1", "rs2"], "BP": [100, 120], "CM": [0.0, 0.0]})
            with mock.patch("pyarrow.parquet.ParquetFile", FakePF):
                with self.assertLogs("LDSC.new", level="WARNING") as log_ctx:
                    ld.SortedR2BlockReader(
                        paths=[str(path)],
                        chrom="1",
                        metadata=metadata,
                        identifier_mode="rsID",
                        r2_bias_mode="unbiased",
                        r2_sample_size=None,
                        genome_build="hg19",
                    )
            self.assertTrue(any("row group" in msg.lower() for msg in log_ctx.output))

    def test_canonical_init_hard_fails_on_multiple_paths(self):
        import pyarrow as pa
        import pyarrow.parquet as pq

        with tempfile.TemporaryDirectory() as tmpdir:
            path1 = Path(tmpdir) / "chr1_a.parquet"
            path2 = Path(tmpdir) / "chr1_b.parquet"
            for path in (path1, path2):
                df = pd.DataFrame(
                    {"CHR": ["1"], "POS_1": [100], "POS_2": [120], "R2": [0.5], "SNP_1": ["rs1"], "SNP_2": ["rs2"]}
                )
                table = pa.Table.from_pandas(df, preserve_index=False)
                meta = {b"ldsc:sorted_by_build": b"hg19", b"ldsc:row_group_size": b"50000"}
                enriched = table.schema.with_metadata({**(table.schema.metadata or {}), **meta})
                with pq.ParquetWriter(str(path), enriched) as writer:
                    writer.write_table(table.cast(enriched))

            metadata = pd.DataFrame({"CHR": ["1", "1"], "SNP": ["rs1", "rs2"], "BP": [100, 120], "CM": [0.0, 0.0]})
            with self.assertRaises(ValueError) as ctx:
                ld.SortedR2BlockReader(
                    paths=[str(path1), str(path2)],
                    chrom="1",
                    metadata=metadata,
                    identifier_mode="rsID",
                    r2_bias_mode="unbiased",
                    r2_sample_size=None,
                    genome_build="hg19",
                )
            self.assertIn("exactly one file", str(ctx.exception))

    def test_canonical_within_block_matrix_correctness(self):
        import pyarrow as pa
        import pyarrow.parquet as pq

        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "chr1.parquet"
            df = pd.DataFrame(
                {
                    "CHR": ["1", "1", "1"],
                    "POS_1": [100, 100, 120],
                    "POS_2": [120, 140, 140],
                    "R2": [0.4, 0.2, 0.6],
                    "SNP_1": ["rs1", "rs1", "rs2"],
                    "SNP_2": ["rs2", "rs3", "rs3"],
                }
            )
            table = pa.Table.from_pandas(df, preserve_index=False)
            meta = {b"ldsc:sorted_by_build": b"hg19", b"ldsc:row_group_size": b"2"}
            enriched = table.schema.with_metadata({**(table.schema.metadata or {}), **meta})
            with pq.ParquetWriter(str(path), enriched) as writer:
                writer.write_table(table.cast(enriched), row_group_size=2)

            metadata = pd.DataFrame(
                {
                    "CHR": ["1", "1", "1"],
                    "SNP": ["rs1", "rs2", "rs3"],
                    "BP": [100, 120, 140],
                    "CM": [0.0, 0.0, 0.0],
                }
            )
            reader = ld.SortedR2BlockReader(
                paths=[str(path)],
                chrom="1",
                metadata=metadata,
                identifier_mode="rsID",
                r2_bias_mode="unbiased",
                r2_sample_size=None,
                genome_build="hg19",
            )
            self.assertEqual(reader._runtime_layout, "canonical")
            matrix = reader.within_block_matrix(l_B=0, c=3)
            expected = np.array([[1.0, 0.4, 0.2], [0.4, 1.0, 0.6], [0.2, 0.6, 1.0]], dtype=np.float32)
            self.assertEqual(matrix.dtype, np.dtype("float32"))
            assert_array_almost_equal(matrix, expected)

    def test_canonical_query_reads_only_overlapping_row_groups(self):
        import pyarrow as pa
        import pyarrow.parquet as pq

        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "chr1.parquet"
            df = pd.DataFrame(
                {
                    "CHR": ["1"] * 6,
                    "POS_1": [100, 100, 200, 200, 300, 300],
                    "POS_2": [120, 130, 220, 230, 320, 330],
                    "R2": [0.4, 0.3, 0.5, 0.6, 0.7, 0.8],
                    "SNP_1": ["rs1", "rs1", "rs3", "rs3", "rs5", "rs5"],
                    "SNP_2": ["rs2", "rs2b", "rs4", "rs4b", "rs6", "rs6b"],
                }
            )
            table = pa.Table.from_pandas(df, preserve_index=False)
            meta = {b"ldsc:sorted_by_build": b"hg19", b"ldsc:row_group_size": b"2"}
            enriched = table.schema.with_metadata({**(table.schema.metadata or {}), **meta})
            with pq.ParquetWriter(str(path), enriched) as writer:
                writer.write_table(table.cast(enriched), row_group_size=2)

            metadata = pd.DataFrame(
                {
                    "CHR": ["1"] * 7,
                    "SNP": [f"rs{i + 1}" for i in range(7)],
                    "BP": [100, 120, 130, 200, 220, 300, 320],
                    "CM": [0.0] * 7,
                }
            )
            reader = ld.SortedR2BlockReader(
                paths=[str(path)],
                chrom="1",
                metadata=metadata,
                identifier_mode="rsID",
                r2_bias_mode="unbiased",
                r2_sample_size=None,
                genome_build="hg19",
            )
            rg_idxs_used = []
            original_read_row_groups = reader._pf.read_row_groups

            def mock_read_row_groups(indexes, **kwargs):
                rg_idxs_used.extend(indexes)
                return original_read_row_groups(indexes, **kwargs)

            reader._pf.read_row_groups = mock_read_row_groups
            rows = reader._query_union_rows(100, 130)
            self.assertEqual(sorted(set(rg_idxs_used)), [0])
            self.assertEqual(len(rows), 2)
