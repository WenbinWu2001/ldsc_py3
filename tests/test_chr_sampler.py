import importlib.util
from pathlib import Path
import sys
import tempfile
import unittest

import pandas as pd


SRC = Path(__file__).resolve().parents[1] / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from ldsc._chr_sampler import sample_frame_from_chr_pattern


_HAS_PYARROW = importlib.util.find_spec("pyarrow") is not None


class ChrSamplerTest(unittest.TestCase):
    def test_samples_first_existing_chromosome_token_and_renames_chr_pos(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            path = tmpdir / "baseline.2.annot.gz"
            pd.DataFrame(
                {
                    "chromosome": ["2", "2"],
                    "base_pair": [101, 202],
                    "annotation": [1.0, 0.0],
                }
            ).to_csv(path, sep="\t", index=False)

            frame, sampled_path = sample_frame_from_chr_pattern(
                [str(tmpdir / "baseline.@.annot.gz")],
                chromosomes=("2",),
                nrows=1,
                context="annotation",
            )

            self.assertEqual(sampled_path, str(path))
            self.assertEqual(frame.columns.tolist(), ["CHR", "POS"])
            self.assertEqual(frame.to_dict("records"), [{"CHR": 2, "POS": 101}])

    @unittest.skipUnless(_HAS_PYARROW, "pyarrow is required for parquet sampling coverage")
    def test_samples_parquet_tokens(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            path = tmpdir / "r2.1.parquet"
            pd.DataFrame({"CHR": ["1"], "BP": [12345], "R2": [0.7]}).to_parquet(path, index=False)

            frame, sampled_path = sample_frame_from_chr_pattern(
                [str(tmpdir / "r2.@.parquet")],
                context="reference panel",
            )

            self.assertEqual(sampled_path, str(path))
            self.assertEqual(frame.to_dict("records"), [{"CHR": "1", "POS": 12345}])

    def test_requires_chromosome_placeholder_token(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            path = tmpdir / "baseline.1.annot.gz"
            path.write_text("CHR\tPOS\n1\t100\n", encoding="utf-8")

            with self.assertRaisesRegex(ValueError, "--genome-build"):
                sample_frame_from_chr_pattern([str(path)], context="annotation")

    def test_reports_when_no_candidate_chromosome_file_exists(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)

            with self.assertRaisesRegex(FileNotFoundError, "--genome-build"):
                sample_frame_from_chr_pattern(
                    [str(tmpdir / "baseline.@.annot.gz")],
                    chromosomes=("2",),
                    context="annotation",
                )


if __name__ == "__main__":
    unittest.main()
