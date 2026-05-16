from __future__ import annotations

import gzip
from pathlib import Path
import sys
import tempfile
import unittest

import pandas as pd

SRC = Path(__file__).resolve().parents[1] / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))


class Hm3ReferenceBuilderTest(unittest.TestCase):
    def test_build_hm3_chr_pos_reference_filters_and_caps_rows(self):
        from ldsc.hm3_reference import build_hm3_chr_pos_reference

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
                    "MAF": [0.25, 0.10, 0.30, 0.40, 0.50],
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


if __name__ == "__main__":
    unittest.main()
