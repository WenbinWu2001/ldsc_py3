"""Reading one R2 matrix chromosome through ``ldsc._refpanel.r2_format.R2MatrixChromosome``.

Each test writes a small chromosome with the package writers
(``tests.ref_panel_helpers.write_r2_matrix_chromosome``) and opens it as
``ldscore`` does: the sidecar and its binding to the R2 Parquet file, the map
from R2 matrix rows to retained SNPs, pair decoding and dequantization, and
the check that a requested LD window fits the recorded one.
"""
from pathlib import Path
import tempfile
import unittest

import numpy as np
import pandas as pd

from ldsc.config import LDWindow
from ldsc.errors import LDSCInputError, LDSCUsageError
from ldsc._refpanel import r2_format
from tests.ref_panel_helpers import open_r2_matrix_chromosome


def _set_footer_metadata(path, **values):
    """Overwrite ``ldsc:<key>`` entries in one R2 Parquet file footer."""
    import pyarrow.parquet as pq

    table = pq.read_table(path)
    metadata = dict(table.schema.metadata)
    metadata.update({f"ldsc:{key}".encode(): str(value).encode() for key, value in values.items()})
    pq.write_table(table.replace_schema_metadata(metadata), path)


class LDScoreParquetNormalizationTest(unittest.TestCase):

    def test_sorted_r2_block_reader_rejects_auto_genome_build_at_kernel_boundary(self):
        panel = pd.DataFrame({"CHR": ["1", "1"], "POS": [100, 110], "SNP": ["rs0", "rs1"],
                              "A1": ["A", "A"], "A2": ["C", "C"], "MAF": [0.3, 0.3]})
        with tempfile.TemporaryDirectory() as tmpdir:
            chromosome = open_r2_matrix_chromosome(Path(tmpdir), panel)
            with self.assertRaisesRegex(AssertionError, "concrete"):
                chromosome.pair_reader(panel, "chr_pos", "auto")


class IndexBindingTest(unittest.TestCase):
    PANEL = pd.DataFrame({"CHR": ["1", "1", "1"], "POS": [10, 20, 30],
                          "SNP": ["rsA", "rsB", "rsC"], "A1": ["A", "C", "G"],
                          "A2": ["G", "T", "A"], "CM": [0.0, 0.0, 0.0], "MAF": [0.2, 0.3, 0.4]})

    def test_sidecar_reads_canonical_columns_once(self):
        with tempfile.TemporaryDirectory() as tmp:
            chromosome = open_r2_matrix_chromosome(Path(tmp), self.PANEL)
            sidecar = chromosome.sidecar()
            self.assertEqual(list(sidecar.columns), list(r2_format.SIDECAR_COLUMNS))
            self.assertEqual(sidecar["CHR"].tolist(), ["1", "1", "1"])
            self.assertEqual(sidecar["POS"].dtype, np.int64)
            self.assertEqual(sidecar["POS"].tolist(), [10, 20, 30])
            self.assertEqual(sidecar["MAF"].tolist(), [0.2, 0.3, 0.4])
            self.assertIs(chromosome.sidecar(), sidecar)

    def test_sidecar_binding_rejects_hash_and_row_count_mismatch(self):
        with tempfile.TemporaryDirectory() as tmp:
            chromosome = open_r2_matrix_chromosome(Path(tmp), self.PANEL)
            chromosome.sidecar()  # matched pair: no raise
            for values, reason in (({"sidecar_identity_sha256": "0" * 64}, "ldsc:sidecar_identity_sha256"),
                                   ({"n_snps": 999}, "ldsc:n_snps=999")):
                _set_footer_metadata(chromosome.parquet_path, **values)
                reopened = r2_format.R2MatrixChromosome.open(Path(tmp), "1", command="ldscore")
                with self.assertRaisesRegex(LDSCInputError, reason):
                    reopened.sidecar()


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
    def test_streaming_drops_unretained_endpoint(self):
        full = pd.DataFrame({"CHR": ["1"] * 4, "POS": [10, 20, 30, 40],
                             "SNP": ["a", "b", "c", "d"], "A1": list("ACGT"),
                             "A2": list("GTAC"), "CM": [0.0] * 4, "MAF": [0.3] * 4})
        rows = [
            {"i": 0, "j": 1, "R2": 0.5, "sign": "+"},
            {"i": 0, "j": 2, "R2": 0.2, "sign": "+"},
            {"i": 1, "j": 2, "R2": 0.4, "sign": "-"},
            {"i": 2, "j": 3, "R2": 0.9, "sign": "+"},
        ]
        with tempfile.TemporaryDirectory() as tmp:
            chromosome = open_r2_matrix_chromosome(Path(tmp), full, rows)
            # retained = panel rows 0,1,2 (drop pos40); matrix order = [10,20,30]
            retained = full.iloc[[0, 1, 2]].reset_index(drop=True)
            reader = chromosome.pair_reader(retained, "chr_pos", "hg19")
            seen = {}
            for i, j, r2v in reader.iter_all_pairs():
                for a, b, v in zip(i.tolist(), j.tolist(), r2v.tolist()):
                    seen[(a, b)] = v
            reader.close()
            # pair (2,3) touches the dropped SNP (pos40) -> excluded by remap
            self.assertEqual(sorted(seen), [(0, 1), (0, 2), (1, 2)])
            # off-diagonal R2 is int16-quantized: tolerance is the half-step (1.5e-5)
            self.assertAlmostEqual(seen[(0, 1)], 0.5, delta=2e-5)
            self.assertAlmostEqual(seen[(0, 2)], 0.2, delta=2e-5)
            self.assertAlmostEqual(seen[(1, 2)], 0.4, delta=2e-5)


class R2MatrixWindowValidationTest(unittest.TestCase):
    TWO_SNPS = pd.DataFrame({"CHR": ["1", "1"], "POS": [100, 200], "SNP": ["rs1", "rs2"], "A1": ["A", "A"],
                             "A2": ["C", "C"], "CM": [0.0, 1.0], "MAF": [0.3, 0.3]})

    def _chromosome(self, tmp, ld_window, panel=None, pairs=()):
        return open_r2_matrix_chromosome(Path(tmp), self.TWO_SNPS if panel is None else panel, pairs, ld_window=ld_window)

    def test_ldscore_accepts_equal_recorded_r2_matrix_window(self):
        with tempfile.TemporaryDirectory() as tmp:
            chromosome = self._chromosome(tmp, ("cm", 1.0))
            chromosome.check_ld_window(LDWindow("cm", 1.0))

    def test_ldscore_rejects_window_wider_than_recorded_r2_matrix_window(self):
        panel = pd.DataFrame({"CHR": ["1"] * 3, "POS": [100, 200, 300], "SNP": ["rs1", "rs2", "rs3"],
                              "A1": ["A"] * 3, "A2": ["C"] * 3, "CM": [0.0, 0.1, 0.2], "MAF": [0.3] * 3})
        with tempfile.TemporaryDirectory() as tmp:
            chromosome = self._chromosome(tmp, ("kb", 0.1), panel, [{"i": 0, "j": 1, "R2": 0.5, "sign": "+"}])
            with self.assertRaisesRegex(
                LDSCUsageError,
                "user-requested LD window `--ld-wind-kb 1.0` is wider than the input R2 matrix window `--ld-wind-kb 0.1`",
            ):
                chromosome.check_ld_window(LDWindow("kb", 1.0))

    def test_ldscore_rejects_any_strictly_larger_same_mode_window(self):
        with tempfile.TemporaryDirectory() as tmp:
            chromosome = self._chromosome(tmp, ("cm", 1.0))
            with self.assertRaisesRegex(LDSCUsageError, "is wider than the input R2 matrix window"):
                chromosome.check_ld_window(LDWindow("cm", float(np.nextafter(1.0, 2.0))))

    def test_ldscore_rejects_requested_window_mode_differing_from_r2_matrix_mode(self):
        with tempfile.TemporaryDirectory() as tmp:
            chromosome = self._chromosome(tmp, ("cm", 1.0))
            with self.assertRaises(LDSCUsageError) as caught:
                chromosome.check_ld_window(LDWindow("kb", 100.0))

        message = str(caught.exception)
        self.assertIn("with `--ld-wind-kb 100.0`", message)
        self.assertIn("the R2 matrix was built with `--ld-wind-cm 1.0`", message)
        self.assertIn("windows in different units cannot be checked for coverage", message)
        self.assertIn("ldsc build-r2-matrix --ld-wind-kb 100.0", message)

