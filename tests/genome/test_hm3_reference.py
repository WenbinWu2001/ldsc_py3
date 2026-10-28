from __future__ import annotations

import runpy

import gzip
from pathlib import Path
import tempfile
import unittest

import pandas as pd

TOOL_DIR = Path(__file__).resolve().parents[2] / "tools" / "hm3"



class Hm3ReferenceBuilderTest(unittest.TestCase):
    def test_build_hm3_chr_pos_reference_filters_and_caps_rows(self):
        build_hm3_chr_pos_reference = runpy.run_path(str(TOOL_DIR / "build_hm3_chr_pos_reference.py"))["build_hm3_chr_pos_reference"]

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            source = tmpdir / "hm3.tsv.gz"
            output = tmpdir / "reference.tsv.gz"
            frame = pd.DataFrame(
                {
                    "CHR": ["1", "1", "1", "1", "2"],
                    "hg19_POS": [10, 20, 30, 40, 50],
                    "hg38_POS": [110, 120, 130, 140, 150],
                    "SNP": ["rs1", "rs2", "rs3", "rs4", "rs5"],
                    "A1": ["A", "A", "A", "C", "A"],
                    "A2": ["C", "G", "T", "G", "C"],
                }
            )
            with gzip.open(source, "wt", encoding="utf-8") as handle:
                frame.to_csv(handle, sep="\t", index=False)

            result = build_hm3_chr_pos_reference(source, output, snps_per_chromosome=1)

            self.assertEqual(result["CHR"].tolist(), ["1", "2"])
            self.assertEqual(result["hg19_POS"].tolist(), [10, 50])
            self.assertEqual(result["hg38_POS"].tolist(), [110, 150])
            written = pd.read_csv(output, sep="\t", dtype={"CHR": "string"})
            pd.testing.assert_frame_equal(written, result)

    def test_default_samples_up_to_1000_snps_per_autosome(self):
        build_hm3_chr_pos_reference = runpy.run_path(str(TOOL_DIR / "build_hm3_chr_pos_reference.py"))["build_hm3_chr_pos_reference"]

        n_chr1, n_chr2 = 1200, 5
        positions = list(range(1, n_chr1 + 1)) + list(range(1, n_chr2 + 1))
        frame = pd.DataFrame(
            {
                "CHR": ["1"] * n_chr1 + ["2"] * n_chr2,
                "hg19_POS": [100 * pos for pos in positions],
                "hg38_POS": [100 * pos + 7 for pos in positions],
                "SNP": [f"rs{i}" for i in range(n_chr1 + n_chr2)],
                "A1": ["A"] * (n_chr1 + n_chr2),
                "A2": ["G"] * (n_chr1 + n_chr2),
            }
        )
        with tempfile.TemporaryDirectory() as tmpdir:
            source = Path(tmpdir) / "hm3.tsv.gz"
            frame.to_csv(source, sep="\t", index=False, compression="gzip")

            result = build_hm3_chr_pos_reference(source, Path(tmpdir) / "reference.tsv.gz")

        self.assertEqual(result["CHR"].value_counts().to_dict(), {"1": 1000, "2": 5})
        chr1 = result.loc[result["CHR"] == "1", "hg38_POS"]
        self.assertTrue(chr1.is_monotonic_increasing and chr1.is_unique)
        self.assertEqual((chr1.iloc[0], chr1.iloc[-1]), (107, 100 * n_chr1 + 7))


if __name__ == "__main__":
    unittest.main()
