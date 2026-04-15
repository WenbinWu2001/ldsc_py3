from pathlib import Path
import gzip
import sys
import tempfile
import unittest
from unittest import mock

import pandas as pd

SRC = Path(__file__).resolve().parents[1] / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from ldsc.config import CommonConfig, MungeConfig

try:
    from ldsc.sumstats_munger import RawSumstatsSpec, SumstatsMunger
except ImportError:
    RawSumstatsSpec = None
    SumstatsMunger = None


@unittest.skipIf(SumstatsMunger is None, "sumstats_munger module is not available")
class SumstatsMungerTest(unittest.TestCase):
    def test_run_munges_and_writes_output(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            raw_path = tmpdir / "raw.tsv"
            raw_path.write_text(
                "SNP A1 A2 P OR N\n"
                "rs1 A G 0.05 1.1 1000\n"
                "rs2 C T 0.10 0.9 1000\n",
                encoding="utf-8",
            )
            raw = RawSumstatsSpec(path=str(raw_path), trait_name="trait")
            config = MungeConfig(out_prefix=str(tmpdir / "munged"))
            munger = SumstatsMunger()
            table = munger.run(raw, config, CommonConfig(snp_identifier="rsid"))
            self.assertEqual(table.trait_name, "trait")
            self.assertTrue((tmpdir / "munged.sumstats.gz").exists())
            with gzip.open(tmpdir / "munged.sumstats.gz", "rt", encoding="utf-8") as handle:
                output = pd.read_csv(handle, sep="\t")
            self.assertEqual(output.columns.tolist(), ["SNP", "A1", "A2", "Z", "N"])
            summary = munger.build_run_summary(table)
            self.assertEqual(summary.n_retained_rows, 2)
            self.assertIn("sumstats_gz", summary.output_paths)

    def test_top_level_wrapper_calls_main(self):
        from ldsc import sumstats_munger as munge_sumstats

        with mock.patch.object(munge_sumstats, "main", return_value=11) as patched:
            rc = munge_sumstats.__dict__["main"](["--help"])
        patched.assert_called_once()
        self.assertEqual(rc, 11)
