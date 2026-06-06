from pathlib import Path
import importlib
import sys
import tempfile
import gzip
import unittest
from unittest import mock

import numpy as np
from numpy.testing import assert_array_almost_equal, assert_array_equal
import pandas as pd

SRC = Path(__file__).resolve().parents[1] / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from ldsc._kernel import formats as ps
from ldsc.errors import LDSCConfigError, LDSCInputError, LDSCInternalError

try:
    import bitarray as ba
    from ldsc._kernel import ldscore as ld
except ModuleNotFoundError:
    ba = None
    ld = None


PLINK_FIXTURES = Path(__file__).resolve().parent / "fixtures" / "plink"


def pair_chunk(rows):
    """Build one PairColumns chunk (i, j, r2 float32, sign int8) from dict rows."""
    rows = list(rows)
    i = np.array([r["i"] for r in rows], dtype=np.int64)
    j = np.array([r["j"] for r in rows], dtype=np.int64)
    r2 = np.array([r["R2"] for r in rows], dtype=np.float32)
    sign = np.array([1 if r["sign"] == "+" else -1 for r in rows], dtype=np.int8)
    return (i, j, r2, sign)


def dict_chunks(rows):
    """Wrap legacy dict rows as the chunk iterable write_r2_parquet expects."""
    rows = list(rows)
    return [pair_chunk(rows)] if rows else iter([])


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
        with self.assertRaises(LDSCInputError):
            ld.PlinkBEDFile(str(PLINK_FIXTURES / "plink.bim"), 9, self.bim)

    def test_next_snps_errors(self):
        bed = ld.PlinkBEDFile(str(PLINK_FIXTURES / "plink.bed"), self.sample_count, self.bim)
        with self.assertRaises(LDSCConfigError):
            bed.nextSNPs(0)
        with self.assertRaises(LDSCInternalError):
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

    def test_next_snps_dtype_default_is_float64(self):
        bed = ld.PlinkBEDFile(str(PLINK_FIXTURES / "plink.bed"), self.sample_count, self.bim)
        x = bed.nextSNPs(3)
        self.assertEqual(x.dtype, np.float64)

    def test_next_snps_dtype_float32_matches_float64(self):
        width = 4
        bed = ld.PlinkBEDFile(str(PLINK_FIXTURES / "plink.bed"), self.sample_count, self.bim)
        x64 = bed.nextSNPs(width)
        bed._currentSNP -= width
        x32 = bed.nextSNPs(width, dtype=np.float32)
        self.assertEqual(x32.dtype, np.float32)
        np.testing.assert_allclose(x32, x64, rtol=1e-4, atol=1e-5)


@unittest.skipIf(ld is None, "ldscore kernel is not available")
@unittest.skipUnless(_has_module("pyarrow"), "pyarrow dependency is not installed")


class IndexParquetRuntimeTest(unittest.TestCase):
    """Tests for the 4-column index parquet schema path."""

    def _panel_meta(self, snps, bps, *, a1="A", a2="C"):
        return pd.DataFrame(
            {
                "CHR": ["1"] * len(snps),
                "POS": list(bps),
                "SNP": list(snps),
                "A1": [a1] * len(snps),
                "A2": [a2] * len(snps),
                "CM": [0.0] * len(snps),
                "MAF": [0.3] * len(snps),
            }
        )

    def _reader_meta(self, snps, bps):
        return pd.DataFrame(
            {
                "CHR": ["1"] * len(snps),
                "SNP": list(snps),
                "BP": list(bps),
                "CM": [0.0] * len(snps),
            }
        )

    def _write_index_panel(self, tmpdir, panel_meta, pairs, *, genome_build="hg19", row_group_size=2):
        from ldsc._kernel import ref_panel_builder as kb
        from ldsc._kernel.snp_identity import sidecar_identity_sha256

        r2 = Path(tmpdir) / "chr1_r2.parquet"
        meta = Path(tmpdir) / "chr1_meta.tsv.gz"
        with gzip.open(meta, "wt") as handle:
            handle.write("# ldsc:schema_version=2\n")
            panel_meta.to_csv(handle, sep="\t", index=False)
        kb.write_r2_parquet(
            pair_chunks=dict_chunks(pairs),
            path=r2,
            genome_build=genome_build,
            n_samples=10,
            snp_identifier="chr_pos",
            min_r2=0.0,
            n_snps=len(panel_meta),
            sidecar_identity_sha256=sidecar_identity_sha256(panel_meta),
            row_group_size=row_group_size,
        )
        return r2

    def _reader(self, r2, reader_meta, *, genome_build="hg19"):
        return ld.SortedR2BlockReader(
            paths=[str(r2)],
            chrom="1",
            metadata=reader_meta,
            identifier_mode="rsid",
            r2_bias_mode="unbiased",
            r2_sample_size=None,
            genome_build=genome_build,
        )

    def _four_snp_panel(self):
        return self._panel_meta(["rs1", "rs2", "rs3", "rs4"], [100, 120, 140, 160])

    def _four_snp_reader_meta(self):
        return self._reader_meta(["rs1", "rs2", "rs3", "rs4"], [100, 120, 140, 160])

    def _four_snp_pairs(self):
        return [
            {"i": 0, "j": 1, "R2": 0.4, "sign": "+"},
            {"i": 0, "j": 2, "R2": 0.2, "sign": "+"},
            {"i": 1, "j": 2, "R2": 0.6, "sign": "+"},
            {"i": 2, "j": 3, "R2": 0.5, "sign": "+"},
        ]

    def test_index_init_raises_on_build_mismatch(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            panel = self._panel_meta(["rs1", "rs2"], [100, 120])
            r2 = self._write_index_panel(tmpdir, panel, [{"i": 0, "j": 1, "R2": 0.5, "sign": "+"}], genome_build="hg19")
            with self.assertRaises(LDSCInputError) as ctx:
                self._reader(r2, self._reader_meta(["rs1", "rs2"], [100, 120]), genome_build="hg38")
            self.assertIn("hg19", str(ctx.exception))
            self.assertIn("hg38", str(ctx.exception))

    def test_index_init_warns_on_coarse_row_groups(self):
        import pyarrow.parquet as pq

        with tempfile.TemporaryDirectory() as tmpdir:
            panel = self._panel_meta(["rs1", "rs2"], [100, 120])
            r2 = self._write_index_panel(tmpdir, panel, [{"i": 0, "j": 1, "R2": 0.5, "sign": "+"}], row_group_size=50000)
            real_pf_cls = pq.ParquetFile

            class FakeMeta:
                num_rows = 600_000
                num_row_groups = 1

                def row_group(self, index):
                    return real_pf_cls(str(r2)).metadata.row_group(index)

            class FakePF:
                def __init__(self, parquet_path):
                    self._real = real_pf_cls(parquet_path)
                    self.metadata = FakeMeta()
                    self.schema_arrow = self._real.schema_arrow

                def read_row_group(self, *args, **kwargs):
                    return self._real.read_row_group(*args, **kwargs)

            with mock.patch("pyarrow.parquet.ParquetFile", FakePF):
                with self.assertLogs("LDSC.ldscore", level="WARNING") as log_ctx:
                    self._reader(r2, self._reader_meta(["rs1", "rs2"], [100, 120]))
            self.assertTrue(any("row group" in msg.lower() for msg in log_ctx.output))

    def test_index_init_hard_fails_on_multiple_paths(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            panel = self._panel_meta(["rs1", "rs2"], [100, 120])
            r2a = self._write_index_panel(tmpdir, panel, [{"i": 0, "j": 1, "R2": 0.5, "sign": "+"}])
            r2b = Path(tmpdir) / "chr1_b_r2.parquet"
            import shutil

            shutil.copy(r2a, r2b)
            with self.assertRaises(LDSCInputError) as ctx:
                ld.SortedR2BlockReader(
                    paths=[str(r2a), str(r2b)],
                    chrom="1",
                    metadata=self._reader_meta(["rs1", "rs2"], [100, 120]),
                    identifier_mode="rsid",
                    r2_bias_mode="unbiased",
                    r2_sample_size=None,
                    genome_build="hg19",
                )
            self.assertIn("exactly one file", str(ctx.exception))

    def test_build_index_remap_rejects_collapsing_sidecar(self):
        # Two distinct panel rows share an rsID; under rsid mode they collapse
        # onto the single retained "rsX", which must be a hard error.
        full = pd.DataFrame({"SNP": ["rsX", "rsY", "rsX"]})
        retained = pd.DataFrame({"SNP": ["rsX", "rsY"]})
        with self.assertRaises(LDSCInputError) as ctx:
            ld.build_index_remap(full, retained, "rsid")
        self.assertIn("collapse", str(ctx.exception).lower())

    def test_build_index_remap_injective_builds_retained_build_idx(self):
        # rsB is absent from the retained universe (maps to -1); rsA, rsC are kept.
        full = pd.DataFrame({"SNP": ["rsA", "rsB", "rsC"]})
        retained = pd.DataFrame({"SNP": ["rsA", "rsC"]})
        remap, retained_build_idx = ld.build_index_remap(full, retained, "rsid")
        np.testing.assert_array_equal(remap, [0, -1, 1])
        np.testing.assert_array_equal(retained_build_idx, [0, 2])

    def test_iter_all_pairs_yields_every_stored_pair_once(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            snps = ["rs1", "rs2", "rs3", "rs4"]
            bps = [100, 120, 140, 160]
            panel = self._panel_meta(snps, bps)
            pairs = [
                {"i": 0, "j": 1, "R2": 0.4, "sign": "+"},
                {"i": 0, "j": 2, "R2": 0.2, "sign": "+"},
                {"i": 1, "j": 2, "R2": 0.6, "sign": "+"},
                {"i": 2, "j": 3, "R2": 0.5, "sign": "+"},
            ]
            r2 = self._write_index_panel(tmpdir, panel, pairs, row_group_size=2)
            reader = self._reader(r2, self._reader_meta(snps, bps))
            seen = {}
            for i, j, r2v in reader.iter_all_pairs():
                self.assertEqual(i.dtype, np.dtype("int32"))
                self.assertEqual(r2v.dtype, np.dtype("float32"))
                for a, b, v in zip(i.tolist(), j.tolist(), r2v.tolist()):
                    seen[(a, b)] = v
            self.assertEqual(sorted(seen), [(0, 1), (0, 2), (1, 2), (2, 3)])
            # int16-quantized R2: tolerance is the half-step
            np.testing.assert_allclose(
                [seen[k] for k in [(0, 1), (0, 2), (1, 2), (2, 3)]],
                [0.4, 0.2, 0.6, 0.5], rtol=0, atol=2e-5,
            )

    def test_accumulate_full_window_matches_R_at_annot(self):
        annot = np.ones((4, 1), dtype=np.float32)
        cor_sum = annot.astype(np.float64).copy()  # diagonal R2 = 1
        block_left = np.zeros(4, dtype=np.int64)
        i = np.array([0, 0, 1, 2], dtype=np.int32)
        j = np.array([1, 2, 2, 3], dtype=np.int32)
        r2 = np.array([0.4, 0.2, 0.6, 0.5], dtype=np.float32)
        ld._accumulate_pair_contributions(cor_sum, i, j, r2, annot, block_left)
        # Row sums of symmetric R (unit diagonal):
        np.testing.assert_allclose(cor_sum[:, 0], [1.6, 2.0, 2.3, 1.5], atol=1e-6)

    def test_accumulate_handles_duplicate_left_index(self):
        # SNP 0 pairs with SNP 1 and SNP 2 — both must accumulate into cor_sum[0].
        annot = np.array([[1.0], [2.0], [3.0], [0.0]], dtype=np.float32)
        cor_sum = np.zeros((4, 1), dtype=np.float64)
        block_left = np.zeros(4, dtype=np.int64)
        i = np.array([0, 0], dtype=np.int32)
        j = np.array([1, 2], dtype=np.int32)
        r2 = np.array([0.5, 0.5], dtype=np.float32)
        ld._accumulate_pair_contributions(cor_sum, i, j, r2, annot, block_left)
        np.testing.assert_allclose(cor_sum[0, 0], 0.5 * 2.0 + 0.5 * 3.0)  # 2.5
        np.testing.assert_allclose(cor_sum[1, 0], 0.5 * 1.0)
        np.testing.assert_allclose(cor_sum[2, 0], 0.5 * 1.0)

    def test_accumulate_window_filter_drops_out_of_window_pairs(self):
        annot = np.ones((4, 1), dtype=np.float32)
        cor_sum = np.zeros((4, 1), dtype=np.float64)
        block_left = np.array([0, 0, 1, 2], dtype=np.int64)  # SNP2 window starts at 1
        i = np.array([0, 0, 1, 2], dtype=np.int32)
        j = np.array([1, 2, 2, 3], dtype=np.int32)
        r2 = np.array([0.4, 0.2, 0.6, 0.5], dtype=np.float32)
        # (0,2) dropped: 0 < block_left[2]=1. Others kept.
        ld._accumulate_pair_contributions(cor_sum, i, j, r2, annot, block_left)
        np.testing.assert_allclose(cor_sum[:, 0], [0.4, 1.0, 1.1, 0.5], atol=1e-6)
