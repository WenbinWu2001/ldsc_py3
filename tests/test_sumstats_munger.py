from pathlib import Path
import gzip
import sys
import tempfile
import unittest
import warnings
from unittest import mock

import pandas as pd

SRC = Path(__file__).resolve().parents[1] / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from ldsc.config import GlobalConfig, MungeConfig

try:
    import ldsc
    from ldsc.sumstats_munger import RawSumstatsSpec, SumstatsMunger
except ImportError:
    ldsc = None
    RawSumstatsSpec = None
    SumstatsMunger = None


@unittest.skipIf(SumstatsMunger is None, "sumstats_munger module is not available")
class SumstatsMungerTest(unittest.TestCase):
    def test_load_sumstats_reads_curated_sumstats_gz_with_exact_one_glob(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            sumstats_path = tmpdir / "trait.v1.sumstats.gz"
            with gzip.open(sumstats_path, "wt", encoding="utf-8") as handle:
                handle.write("SNP\tZ\tN\tA1\tA2\tFRQ\nrs1\t1.5\t1000\tA\tG\t0.2\n")

            self.assertTrue(hasattr(ldsc, "load_sumstats"))
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                table = ldsc.load_sumstats(str(tmpdir / "trait*.sumstats.gz"), trait_name="trait")

            self.assertEqual(table.source_path, str(sumstats_path))
            self.assertEqual(table.trait_name, "trait")
            self.assertTrue(table.has_alleles)
            self.assertEqual(table.data.columns.tolist(), ["SNP", "N", "Z", "A1", "A2", "FRQ"])
            self.assertEqual(table.data.loc[0, "SNP"], "rs1")
            self.assertIsNotNone(table.config_snapshot)
            self.assertTrue(any("cannot recover the GlobalConfig" in str(item.message) for item in caught))

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
            table = munger.run(raw, config, GlobalConfig(snp_identifier="rsid"))
            self.assertEqual(table.trait_name, "trait")
            self.assertEqual(table.config_snapshot, GlobalConfig(snp_identifier="rsid"))
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

    def test_run_creates_parent_directory_for_out_prefix(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            raw_path = tmpdir / "raw.tsv"
            raw_path.write_text(
                "SNP A1 A2 P OR N\n"
                "rs1 A G 0.05 1.0 1000\n",
                encoding="utf-8",
            )
            raw = RawSumstatsSpec(path=str(raw_path), trait_name="trait")
            nested_prefix = tmpdir / "nested" / "dir" / "munged"
            config = MungeConfig(out_prefix=str(nested_prefix))

            table = SumstatsMunger().run(raw, config, GlobalConfig(snp_identifier="rsid"))

            self.assertEqual(len(table.data), 1)
            self.assertTrue((tmpdir / "nested" / "dir" / "munged.sumstats.gz").exists())
            self.assertTrue((tmpdir / "nested" / "dir" / "munged.log").exists())

    def test_run_accepts_path_objects_for_input_and_output(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            raw_path = tmpdir / "raw.tsv"
            raw_path.write_text(
                "SNP A1 A2 P OR N\n"
                "rs1 A G 0.05 1.0 1000\n",
                encoding="utf-8",
            )
            raw = RawSumstatsSpec(path=raw_path, trait_name="trait")
            out_prefix = tmpdir / "nested" / "mdd2025"

            table = SumstatsMunger().run(raw, MungeConfig(out_prefix=out_prefix), GlobalConfig(snp_identifier="rsid"))

            self.assertEqual(len(table.data), 1)
            self.assertEqual(table.source_path, str(raw_path))
            self.assertTrue((tmpdir / "nested" / "mdd2025.sumstats.gz").exists())
            self.assertTrue((tmpdir / "nested" / "mdd2025.log").exists())

    def test_run_resolves_glob_pattern_for_single_sumstats_input(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            raw_path = tmpdir / "trait.raw.tsv"
            raw_path.write_text(
                "SNP A1 A2 P OR N\n"
                "rs1 A G 0.05 1.0 1000\n",
                encoding="utf-8",
            )
            raw = RawSumstatsSpec(path=str(tmpdir / "trait.*.tsv"), trait_name="trait")
            out_prefix = tmpdir / "out" / "munged"

            table = SumstatsMunger().run(raw, MungeConfig(out_prefix=out_prefix), GlobalConfig(snp_identifier="rsid"))

            self.assertEqual(table.source_path, str(raw_path))
            self.assertTrue((tmpdir / "out" / "munged.sumstats.gz").exists())

    def test_run_accepts_id_and_ncas_ncon_header_aliases(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            raw_path = tmpdir / "raw.tsv"
            raw_path.write_text(
                "#CHROM\tPOS\tID\tEA\tNEA\tBETA\tSE\tPVAL\tFCAS\tFCON\tIMPINFO\tNEFF\tNCAS\tNCON\tHETI\tHETDF\tHETPVAL\n"
                "1\t123\trs1\tA\tG\t0.1\t0.01\t0.05\t0.2\t0.8\t0.95\t1000\t400\t600\t0\t1\t0.9\n",
                encoding="utf-8",
            )

            table = SumstatsMunger().run(
                RawSumstatsSpec(path=raw_path, trait_name="trait"),
                MungeConfig(out_prefix=tmpdir / "munged"),
                GlobalConfig(snp_identifier="rsid"),
            )

            self.assertEqual(table.data["SNP"].tolist(), ["rs1"])
            self.assertEqual(table.data["N"].tolist(), [1000.0])
