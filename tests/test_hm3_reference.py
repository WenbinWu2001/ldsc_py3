from __future__ import annotations

import runpy

import gzip
from pathlib import Path
import tempfile
import unittest

import pandas as pd

TOOL_DIR = Path(__file__).resolve().parents[1] / "tools" / "hm3"



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

    def test_filter_reference_candidates_keeps_maf_at_threshold(self):
        _filter_reference_candidates = runpy.run_path(str(TOOL_DIR / "build_hm3_chr_pos_reference.py"))["_filter_reference_candidates"]

        frame = pd.DataFrame(
            {
                "CHR": ["1", "1"],
                "hg19_POS": [10, 20],
                "hg38_POS": [110, 120],
                "SNP": ["rs1", "rs2"],
                "A1": ["A", "A"],
                "A2": ["C", "G"],
                "MAF": [0.2, 0.19],  # 0.2 == min_maf is kept (inclusive >=); 0.19 dropped
            }
        )

        filtered = _filter_reference_candidates(frame, min_maf=0.2)

        self.assertEqual(filtered["SNP"].tolist(), ["rs1"])


if __name__ == "__main__":
    unittest.main()
