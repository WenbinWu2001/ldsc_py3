from pathlib import Path
import contextlib
import gzip
import importlib.util
import io
import json
import sys
import tempfile
import unittest
import warnings
from unittest import mock

import numpy as np
import pandas as pd
from pandas.testing import assert_series_equal

SRC = Path(__file__).resolve().parents[1] / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

_HAS_PYARROW = importlib.util.find_spec("pyarrow") is not None

from ldsc.config import GlobalConfig, MungeConfig

try:
    import ldsc
    from ldsc import sumstats_munger as sumstats_workflow
    from ldsc._kernel import sumstats_munger as kernel_munge
    from ldsc.sumstats_munger import SumstatsMunger
except ImportError:
    ldsc = None
    kernel_munge = None
    sumstats_workflow = None
    SumstatsMunger = None


@unittest.skipIf(SumstatsMunger is None, "sumstats_munger module is not available")
class SumstatsMungerTest(unittest.TestCase):
    def test_build_parser_defaults_chunksize_to_one_million_rows(self):
        parser = sumstats_workflow.build_parser()
        self.assertEqual(parser.get_default("chunksize"), 1_000_000)

    def test_build_parser_defaults_output_format_to_parquet(self):
        parser = sumstats_workflow.build_parser()
        self.assertEqual(parser.get_default("output_format"), "parquet")
        with self.assertRaises(SystemExit):
            parser.parse_args(["--raw-sumstats-file", "raw.tsv", "--output-dir", "out", "--output-format", "csv"])

    def test_build_parser_defaults_log_level_to_info(self):
        parser = sumstats_workflow.build_parser()
        self.assertEqual(parser.get_default("log_level"), "INFO")
        args = parser.parse_args(["--raw-sumstats-file", "raw.tsv", "--output-dir", "out", "--log-level", "DEBUG"])
        self.assertEqual(args.log_level, "DEBUG")

    def test_build_parser_accepts_trait_name(self):
        parser = sumstats_workflow.build_parser()

        args = parser.parse_args(["--raw-sumstats-file", "raw.tsv", "--output-dir", "out", "--trait-name", "MDD"])

        self.assertEqual(args.trait_name, "MDD")

    def test_build_parser_accepts_liftover_flags(self):
        parser = sumstats_workflow.build_parser()

        args = parser.parse_args(
            [
                "--raw-sumstats-file",
                "raw.tsv",
                "--output-dir",
                "out",
                "--target-genome-build",
                "hg38",
                "--liftover-chain-file",
                "hg19ToHg38.over.chain",
            ]
        )

        self.assertEqual(args.target_genome_build, "hg38")
        self.assertEqual(args.liftover_chain_file, "hg19ToHg38.over.chain")
        self.assertFalse(args.use_hm3_quick_liftover)

        args = parser.parse_args(
            [
                "--raw-sumstats-file",
                "raw.tsv",
                "--output-dir",
                "out",
                "--target-genome-build",
                "hg38",
                "--use-hm3-quick-liftover",
            ]
        )
        self.assertTrue(args.use_hm3_quick_liftover)

    def test_build_parser_accepts_daner_old_and_new_not_legacy_flags(self):
        parser = sumstats_workflow.build_parser()

        old_args = parser.parse_args(["--raw-sumstats-file", "raw.tsv", "--output-dir", "out", "--daner-old"])
        new_args = parser.parse_args(["--raw-sumstats-file", "raw.tsv", "--output-dir", "out", "--daner-new"])

        self.assertTrue(old_args.daner_old)
        self.assertFalse(old_args.daner_new)
        self.assertFalse(new_args.daner_old)
        self.assertTrue(new_args.daner_new)
        with self.assertRaises(SystemExit):
            parser.parse_args(["--raw-sumstats-file", "raw.tsv", "--output-dir", "out", "--daner"])
        with self.assertRaises(SystemExit):
            parser.parse_args(["--raw-sumstats-file", "raw.tsv", "--output-dir", "out", "--daner-n"])

    def test_build_parser_uses_raw_sumstats_file_for_raw_input(self):
        parser = sumstats_workflow.build_parser()

        args = parser.parse_args(
            ["--raw-sumstats-file", "raw.tsv", "--output-dir", "out"]
        )

        self.assertEqual(args.raw_sumstats_file, "raw.tsv")
        with self.assertRaises(SystemExit):
            parser.parse_args(["--sumstats-file", "raw.tsv", "--output-dir", "out"])

    def test_main_delegates_to_run_munge_sumstats_from_args(self):
        with mock.patch.object(sumstats_workflow, "run_munge_sumstats_from_args", return_value=mock.sentinel.result) as patched:
            result = sumstats_workflow.main(["--raw-sumstats-file", "raw.tsv", "--output-dir", "out"])

        self.assertIs(result, mock.sentinel.result)
        self.assertEqual(patched.call_args.args[0].raw_sumstats_file, "raw.tsv")
        self.assertEqual(patched.call_args.args[0].output_dir, "out")

    def test_run_munge_sumstats_from_args_delegates_to_service_with_configs(self):
        args = sumstats_workflow.build_parser().parse_args(
            [
                "--raw-sumstats-file",
                "raw.tsv",
                "--output-dir",
                "out",
                "--sumstats-snps-file",
                "keep.tsv",
                "--overwrite",
                "--snp",
                "variant_id",
                "--chr",
                "chrom",
                "--pos",
                "bp",
                "--N",
                "123",
                "--chunksize",
                "17",
                "--ignore",
                "DROP_ME,ALSO_DROP",
                "--signed-sumstats",
                "BETA,0",
                "--daner-new",
                "--output-format",
                "both",
                "--log-level",
                "DEBUG",
                "--snp-identifier",
                "rsid",
                "--trait-name",
                " MDD ",
            ]
        )

        with mock.patch.object(SumstatsMunger, "run", return_value=mock.sentinel.table) as patched:
            result = sumstats_workflow.run_munge_sumstats_from_args(args)

        self.assertIs(result, mock.sentinel.table)
        raw_config, run_config, global_config = patched.call_args.args
        self.assertEqual(raw_config.raw_sumstats_file, "raw.tsv")
        self.assertEqual(raw_config.trait_name, "MDD")
        self.assertEqual(raw_config.column_hints, {"snp": "variant_id", "chr": "chrom", "pos": "bp"})
        self.assertEqual(run_config.output_dir, "out")
        self.assertEqual(run_config.sumstats_snps_file, "keep.tsv")
        self.assertTrue(run_config.overwrite)
        self.assertEqual(run_config.N, 123)
        self.assertEqual(run_config.chunk_size, 17)
        self.assertEqual(run_config.ignore_columns, ("DROP_ME", "ALSO_DROP"))
        self.assertEqual(run_config.signed_sumstats_spec, "BETA,0")
        self.assertTrue(run_config.daner_new)
        self.assertFalse(run_config.daner_old)
        self.assertEqual(run_config.output_format, "both")
        self.assertEqual(global_config, GlobalConfig(snp_identifier="rsid", log_level="DEBUG"))

    def test_run_munge_sumstats_from_args_passes_liftover_config_in_chr_pos_mode(self):
        args = sumstats_workflow.build_parser().parse_args(
            [
                "--raw-sumstats-file",
                "raw.tsv",
                "--output-dir",
                "out",
                "--snp-identifier",
                "chr_pos",
                "--genome-build",
                "hg19",
                "--target-genome-build",
                "hg38",
                "--use-hm3-quick-liftover",
            ]
        )

        with mock.patch.object(SumstatsMunger, "run", return_value=mock.sentinel.table) as patched:
            result = sumstats_workflow.run_munge_sumstats_from_args(args)

        self.assertIs(result, mock.sentinel.table)
        raw_config, run_config, global_config = patched.call_args.args
        self.assertEqual(raw_config.raw_sumstats_file, "raw.tsv")
        self.assertEqual(run_config.target_genome_build, "hg38")
        self.assertTrue(run_config.use_hm3_quick_liftover)
        self.assertIsNone(run_config.liftover_chain_file)
        self.assertEqual(global_config, GlobalConfig(snp_identifier="chr_pos", genome_build="hg19"))

    def test_munge_config_validates_liftover_flags(self):
        config = MungeConfig(
            output_dir="out",
            target_genome_build="GRCh38",
            liftover_chain_file=Path("liftover") / "hg19ToHg38.over.chain",
        )
        self.assertEqual(config.target_genome_build, "hg38")
        self.assertEqual(config.liftover_chain_file, "liftover/hg19ToHg38.over.chain")

        with self.assertRaisesRegex(ValueError, "target_genome_build"):
            MungeConfig(output_dir="out", use_hm3_quick_liftover=True)
        with self.assertRaisesRegex(ValueError, "mutually exclusive"):
            MungeConfig(
                output_dir="out",
                target_genome_build="hg38",
                liftover_chain_file="chain.over",
                use_hm3_quick_liftover=True,
            )

    def test_kernel_reports_actionable_signed_sumstats_format_error(self):
        args = kernel_munge.parser.parse_args(["--signed-sumstats", "BETA"])

        with self.assertRaisesRegex(ValueError, "Invalid --signed-sumstats value 'BETA'.*BETA,0"):
            kernel_munge.parse_flag_cnames(args)

    def test_kernel_p_to_z_matches_legacy_direction_convention(self):
        z = kernel_munge.p_to_z(pd.Series([0.1, 0.1, 0.1]), pd.Series([1, 2, 3]))
        np.testing.assert_allclose(np.asarray(z), [1.644854, 1.644854, 1.644854], atol=1e-5)

    def test_kernel_filters_reject_invalid_p_info_frq_and_alleles(self):
        args = kernel_munge.parser.parse_args("")
        assert_series_equal(
            kernel_munge.filter_pvals(pd.Series([0, 0.1, 1, 2]), args),
            pd.Series([False, True, True, False]),
        )
        assert_series_equal(
            kernel_munge.filter_info(pd.Series([0.8, 1.0, 1.0]), args),
            pd.Series([False, True, True]),
        )
        assert_series_equal(
            kernel_munge.filter_frq(pd.Series([-1, 0, 0.005, 0.4, 0.6, 0.999, 1, 2]), args),
            pd.Series([False, False, False, True, True, False, False, False]),
        )
        assert_series_equal(
            kernel_munge.filter_alleles(pd.Series(["AC", "AG", "DI", "AAT", "RA"])),
            pd.Series([True, True, False, False, False]),
        )

    def test_load_sumstats_reads_curated_sumstats_gz_with_exact_one_glob(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            sumstats_file = tmpdir / "trait.v1.sumstats.gz"
            with gzip.open(sumstats_file, "wt", encoding="utf-8") as handle:
                handle.write("SNP\tZ\tN\tA1\tA2\tFRQ\nrs1\t1.5\t1000\tA\tG\t0.2\n")

            self.assertTrue(hasattr(ldsc, "load_sumstats"))
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                table = ldsc.load_sumstats(str(tmpdir / "trait*.sumstats.gz"), trait_name="trait")

            self.assertEqual(table.source_path, str(sumstats_file))
            self.assertEqual(table.trait_name, "trait")
            self.assertTrue(table.has_alleles)
            self.assertEqual(table.data.columns.tolist(), ["SNP", "N", "Z", "A1", "A2", "FRQ"])
            self.assertEqual(table.data.loc[0, "SNP"], "rs1")
            self.assertIsNone(table.config_snapshot)
            self.assertTrue(any("cannot recover the GlobalConfig" in str(item.message) for item in caught))

    def test_load_sumstats_reads_uncompressed_curated_sumstats(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            sumstats_file = tmpdir / "trait.sumstats"
            sumstats_file.write_text("SNP\tZ\tN\nrs1\t1.5\t1000\n", encoding="utf-8")

            with warnings.catch_warnings(record=True):
                table = ldsc.load_sumstats(sumstats_file, trait_name="trait")

            self.assertEqual(table.data.columns.tolist(), ["SNP", "N", "Z"])
            self.assertEqual(table.data.loc[0, "SNP"], "rs1")

    def test_load_sumstats_recovers_trait_name_from_sidecar(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            sumstats_file = tmpdir / "trait.sumstats.gz"
            with gzip.open(sumstats_file, "wt", encoding="utf-8") as handle:
                handle.write("SNP\tZ\tN\nrs1\t1.5\t1000\n")
            (tmpdir / "trait.metadata.json").write_text(
                json.dumps(
                    {
                        "format": "ldsc.sumstats.v1",
                        "trait_name": "MDD",
                        "config_snapshot": {"snp_identifier": "rsid", "genome_build": None},
                    }
                ),
                encoding="utf-8",
            )

            table = ldsc.load_sumstats(sumstats_file)

            self.assertEqual(table.trait_name, "MDD")

    def test_load_sumstats_explicit_trait_name_overrides_sidecar(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            sumstats_file = tmpdir / "trait.sumstats.gz"
            with gzip.open(sumstats_file, "wt", encoding="utf-8") as handle:
                handle.write("SNP\tZ\tN\nrs1\t1.5\t1000\n")
            (tmpdir / "trait.metadata.json").write_text(
                json.dumps(
                    {
                        "format": "ldsc.sumstats.v1",
                        "trait_name": "MDD",
                        "config_snapshot": {"snp_identifier": "rsid", "genome_build": None},
                    }
                ),
                encoding="utf-8",
            )

            table = ldsc.load_sumstats(sumstats_file, trait_name=" SCZ ")

            self.assertEqual(table.trait_name, "SCZ")

    def test_load_sumstats_rejects_blank_sidecar_trait_name(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            sumstats_file = tmpdir / "trait.sumstats.gz"
            with gzip.open(sumstats_file, "wt", encoding="utf-8") as handle:
                handle.write("SNP\tZ\tN\nrs1\t1.5\t1000\n")
            (tmpdir / "trait.metadata.json").write_text(
                json.dumps(
                    {
                        "format": "ldsc.sumstats.v1",
                        "trait_name": " ",
                        "config_snapshot": {"snp_identifier": "rsid", "genome_build": None},
                    }
                ),
                encoding="utf-8",
            )

            with self.assertRaisesRegex(ValueError, "trait_name"):
                ldsc.load_sumstats(sumstats_file)

    @unittest.skipUnless(_HAS_PYARROW, "pyarrow is required for sumstats parquet coverage")
    def test_load_sumstats_reads_parquet_with_exact_one_glob_and_sidecar(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            sumstats_file = tmpdir / "trait.parquet"
            pd.DataFrame({"SNP": ["rs1"], "CHR": ["1"], "POS": [100], "Z": [1.5], "N": [1000.0]}).to_parquet(
                sumstats_file,
                index=False,
            )
            (tmpdir / "trait.metadata.json").write_text(
                json.dumps(
                    {
                        "format": "ldsc.sumstats.v1",
                        "config_snapshot": {"snp_identifier": "chr_pos", "genome_build": "hg38"},
                    }
                ),
                encoding="utf-8",
            )

            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                table = ldsc.load_sumstats(str(tmpdir / "trait*.parquet"), trait_name="trait")

            self.assertEqual(table.source_path, str(sumstats_file))
            self.assertEqual(table.data.columns.tolist(), ["SNP", "CHR", "POS", "N", "Z"])
            self.assertEqual(table.config_snapshot, GlobalConfig(snp_identifier="chr_pos", genome_build="hg38"))
            self.assertFalse(any("cannot recover the GlobalConfig" in str(item.message) for item in caught))

    def test_load_sumstats_rejects_sidecar_without_config_snapshot(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            sumstats_file = tmpdir / "trait.sumstats.gz"
            with gzip.open(sumstats_file, "wt", encoding="utf-8") as handle:
                handle.write("SNP\tCHR\tPOS\tZ\tN\nrs1\t1\t100\t1.0\t100.0\n")
            (tmpdir / "trait.metadata.json").write_text(
                json.dumps(
                    {
                        "format": "ldsc.sumstats.v1",
                        "trait_name": "trait",
                    }
                ),
                encoding="utf-8",
            )

            with self.assertRaisesRegex(ValueError, "config_snapshot"):
                ldsc.load_sumstats(sumstats_file)

    def test_load_sumstats_rejects_unknown_suffix(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "trait.csv"
            path.write_text("SNP,Z,N\nrs1,1.5,1000\n", encoding="utf-8")

            with self.assertRaisesRegex(ValueError, "Unsupported munged sumstats format"):
                ldsc.load_sumstats(path)

    @unittest.skipUnless(_HAS_PYARROW, "pyarrow is required for default parquet output")
    def test_run_munges_and_writes_parquet_by_default(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            raw_path = tmpdir / "raw.tsv"
            raw_path.write_text(
                "SNP A1 A2 P OR N\n"
                "rs1 A G 0.05 1.1 1000\n"
                "rs2 C T 0.10 0.9 1000\n",
                encoding="utf-8",
            )
            raw = MungeConfig(raw_sumstats_file=str(raw_path), trait_name="trait")
            config = MungeConfig(output_dir=str(tmpdir / "munged"))
            munger = SumstatsMunger()
            table = munger.run(raw, config, GlobalConfig(snp_identifier="rsid"))
            self.assertEqual(table.trait_name, "trait")
            self.assertEqual(table.config_snapshot, GlobalConfig(snp_identifier="rsid"))
            self.assertTrue((tmpdir / "munged" / "sumstats.parquet").exists())
            self.assertFalse((tmpdir / "munged" / "sumstats.sumstats.gz").exists())
            output = pd.read_parquet(tmpdir / "munged" / "sumstats.parquet")
            self.assertEqual(output.columns.tolist(), ["SNP", "CHR", "POS", "A1", "A2", "Z", "N"])
            self.assertTrue(output["CHR"].isna().all())
            self.assertTrue(output["POS"].isna().all())
            self.assertTrue((tmpdir / "munged" / "sumstats.metadata.json").exists())
            metadata = json.loads((tmpdir / "munged" / "sumstats.metadata.json").read_text(encoding="utf-8"))
            self.assertEqual(metadata["trait_name"], "trait")
            self.assertEqual(
                set(metadata),
                {"format", "trait_name", "config_snapshot"},
            )
            summary = munger.build_run_summary(table)
            self.assertEqual(summary.n_retained_rows, 2)
            self.assertIn("sumstats_parquet", summary.output_paths)
            self.assertNotIn("sumstats_gz", summary.output_paths)
            self.assertNotIn("log", summary.output_paths)
            self.assertIn("metadata_json", summary.output_paths)

    def test_run_writes_tsv_gz_when_requested(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            raw_path = tmpdir / "raw.tsv"
            raw_path.write_text(
                "SNP A1 A2 P OR N\n"
                "rs1 A G 0.05 1.0 1000\n",
                encoding="utf-8",
            )

            table = SumstatsMunger().run(
                MungeConfig(raw_sumstats_file=str(raw_path), trait_name="trait"),
                MungeConfig(output_dir=str(tmpdir / "munged"), output_format="tsv.gz"),
                GlobalConfig(snp_identifier="rsid"),
            )

            self.assertEqual(len(table.data), 1)
            self.assertTrue((tmpdir / "munged" / "sumstats.sumstats.gz").exists())
            self.assertFalse((tmpdir / "munged" / "sumstats.parquet").exists())
            with gzip.open(tmpdir / "munged" / "sumstats.sumstats.gz", "rt", encoding="utf-8") as handle:
                output = pd.read_csv(handle, sep="\t")
            self.assertEqual(output.columns.tolist(), ["SNP", "CHR", "POS", "A1", "A2", "Z", "N"])
            metadata = json.loads((tmpdir / "munged" / "sumstats.metadata.json").read_text(encoding="utf-8"))
            self.assertEqual(metadata["trait_name"], "trait")
            self.assertEqual(set(metadata), {"format", "trait_name", "config_snapshot"})

    @unittest.skipUnless(_HAS_PYARROW, "pyarrow is required for sumstats parquet coverage")
    def test_run_writes_both_formats_and_lists_outputs_in_metadata(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            raw_path = tmpdir / "raw.tsv"
            raw_path.write_text(
                "SNP A1 A2 P OR N\n"
                "rs1 A G 0.05 1.0 1000\n",
                encoding="utf-8",
            )

            SumstatsMunger().run(
                MungeConfig(raw_sumstats_file=str(raw_path), trait_name="trait"),
                MungeConfig(output_dir=str(tmpdir / "munged"), output_format="both"),
                GlobalConfig(snp_identifier="rsid"),
            )

            metadata = json.loads((tmpdir / "munged" / "sumstats.metadata.json").read_text(encoding="utf-8"))
            self.assertTrue((tmpdir / "munged" / "sumstats.parquet").exists())
            self.assertTrue((tmpdir / "munged" / "sumstats.sumstats.gz").exists())
            self.assertEqual(set(metadata), {"format", "trait_name", "config_snapshot"})

    @unittest.skipUnless(_HAS_PYARROW, "pyarrow is required for sumstats parquet coverage")
    def test_run_overwrite_removes_unselected_sumstats_sibling(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            raw_path = tmpdir / "raw.tsv"
            raw_path.write_text("SNP A1 A2 P OR N\nrs1 A G 0.05 1.0 1000\n", encoding="utf-8")
            output_dir = tmpdir / "munged"

            SumstatsMunger().run(
                MungeConfig(raw_sumstats_file=raw_path, trait_name="trait"),
                MungeConfig(output_dir=output_dir, output_format="both"),
                GlobalConfig(snp_identifier="rsid"),
            )
            self.assertTrue((output_dir / "sumstats.sumstats.gz").exists())

            SumstatsMunger().run(
                MungeConfig(raw_sumstats_file=raw_path, trait_name="trait"),
                MungeConfig(output_dir=output_dir, output_format="parquet", overwrite=True),
                GlobalConfig(snp_identifier="rsid"),
            )

            metadata = json.loads((output_dir / "sumstats.metadata.json").read_text(encoding="utf-8"))
            self.assertTrue((output_dir / "sumstats.parquet").exists())
            self.assertFalse((output_dir / "sumstats.sumstats.gz").exists())
            self.assertEqual(set(metadata), {"format", "trait_name", "config_snapshot"})

    def test_run_refuses_unselected_owned_sumstats_sibling_without_overwrite(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            raw_path = tmpdir / "raw.tsv"
            raw_path.write_text("SNP A1 A2 P OR N\nrs1 A G 0.05 1.0 1000\n", encoding="utf-8")
            output_dir = tmpdir / "munged"
            output_dir.mkdir()
            stale = output_dir / "sumstats.sumstats.gz"
            stale.write_text("stale\n", encoding="utf-8")

            with mock.patch.object(kernel_munge, "munge_sumstats", side_effect=AssertionError("kernel should not run")):
                with self.assertRaisesRegex(FileExistsError, "overwrite"):
                    SumstatsMunger().run(
                        MungeConfig(raw_sumstats_file=raw_path, trait_name="trait"),
                        MungeConfig(output_dir=output_dir, output_format="parquet"),
                        GlobalConfig(snp_identifier="rsid"),
                    )

            self.assertEqual(stale.read_text(encoding="utf-8"), "stale\n")
            self.assertFalse((output_dir / "sumstats.log").exists())

    @unittest.skipUnless(_HAS_PYARROW, "pyarrow is required for sumstats parquet coverage")
    def test_write_output_accepts_output_format(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            table = sumstats_workflow.SumstatsTable(
                data=pd.DataFrame({"SNP": ["rs1"], "Z": [1.234567], "N": [1000.123], "A1": ["A"], "A2": ["G"]}),
                has_alleles=True,
                source_path="source.tsv",
                trait_name="trait",
                config_snapshot=GlobalConfig(snp_identifier="rsid"),
            )

            path = SumstatsMunger().write_output(table, tmpdir / "out", output_format="both")

            self.assertEqual(path, str(tmpdir / "out" / "sumstats.parquet"))
            self.assertTrue((tmpdir / "out" / "sumstats.parquet").exists())
            self.assertTrue((tmpdir / "out" / "sumstats.sumstats.gz").exists())

    def test_write_output_refuses_unselected_owned_sumstats_sibling_without_overwrite(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            output_dir = tmpdir / "out"
            output_dir.mkdir()
            stale = output_dir / "sumstats.sumstats.gz"
            stale.write_text("stale\n", encoding="utf-8")
            table = sumstats_workflow.SumstatsTable(
                data=pd.DataFrame({"SNP": ["rs1"], "Z": [1.5], "N": [1000.0]}),
                has_alleles=False,
                source_path="source.tsv",
                trait_name="trait",
            )

            with self.assertRaisesRegex(FileExistsError, "overwrite"):
                SumstatsMunger().write_output(table, output_dir, output_format="parquet")

            self.assertEqual(stale.read_text(encoding="utf-8"), "stale\n")

    @unittest.skipUnless(_HAS_PYARROW, "pyarrow is required for sumstats parquet coverage")
    def test_write_output_overwrite_removes_unselected_sumstats_sibling(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            output_dir = tmpdir / "out"
            output_dir.mkdir()
            stale = output_dir / "sumstats.sumstats.gz"
            stale.write_text("stale\n", encoding="utf-8")
            table = sumstats_workflow.SumstatsTable(
                data=pd.DataFrame({"SNP": ["rs1"], "Z": [1.5], "N": [1000.0]}),
                has_alleles=False,
                source_path="source.tsv",
                trait_name="trait",
                config_snapshot=GlobalConfig(snp_identifier="rsid"),
            )

            SumstatsMunger().write_output(table, output_dir, output_format="parquet", overwrite=True)

            self.assertTrue((output_dir / "sumstats.parquet").exists())
            self.assertTrue((output_dir / "sumstats.metadata.json").exists())
            self.assertFalse(stale.exists())

    @unittest.skipUnless(_HAS_PYARROW, "pyarrow is required for sumstats parquet coverage")
    def test_parquet_output_sorts_by_chr_pos_records_row_groups_and_preserves_precision(self):
        import pyarrow.parquet as pq

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            raw_path = tmpdir / "raw.tsv"
            raw_path.write_text("SNP A1 A2 P OR N\nrs1 A G 0.05 1.1 1000\n", encoding="utf-8")
            returned = pd.DataFrame(
                {
                    "SNP": ["rs2", "rs1", "rs4", "rs3"],
                    "CHR": ["2", "1", pd.NA, "1"],
                    "POS": [200, 100, pd.NA, 50],
                    "A1": ["A", "A", "T", "C"],
                    "A2": ["G", "G", "C", "T"],
                    "Z": [1.123456789, 2.987654321, 0.444444444, -0.333333333],
                    "N": [1000.123456, 2000.987654, 4000.777777, 3000.555555],
                }
            )

            with mock.patch.object(kernel_munge, "munge_sumstats", return_value=returned):
                SumstatsMunger().run(
                    MungeConfig(raw_sumstats_file=str(raw_path), trait_name="trait"),
                    MungeConfig(output_dir=str(tmpdir / "munged")),
                    GlobalConfig(snp_identifier="chr_pos", genome_build="hg38"),
                )

            parquet_path = tmpdir / "munged" / "sumstats.parquet"
            output = pd.read_parquet(parquet_path)
            self.assertEqual(output["SNP"].tolist(), ["rs3", "rs1", "rs2", "rs4"])
            self.assertAlmostEqual(output.loc[0, "Z"], -0.333333333)
            self.assertAlmostEqual(output.loc[1, "N"], 2000.987654)
            self.assertTrue(pd.isna(output.loc[3, "CHR"]))
            self.assertTrue(pd.isna(output.loc[3, "POS"]))
            parquet_file = pq.ParquetFile(parquet_path)
            self.assertEqual(parquet_file.num_row_groups, 3)
            metadata = json.loads((tmpdir / "munged" / "sumstats.metadata.json").read_text(encoding="utf-8"))
            self.assertEqual(set(metadata), {"format", "trait_name", "config_snapshot"})

    @unittest.skipUnless(_HAS_PYARROW, "pyarrow is required for default parquet output")
    def test_run_writes_thin_sidecar_and_logs_coordinate_liftover_provenance(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            raw_path = tmpdir / "raw.tsv"
            raw_path.write_text("CHR POS SNP A1 A2 P BETA N\n1 100 rs1 A G 0.05 0.1 1000\n", encoding="utf-8")

            def fake_munge(args, p=False):
                args._coordinate_metadata = {
                    "format": "ldsc.sumstats.v1",
                    "snp_identifier": "chr_pos",
                    "genome_build": "hg38",
                    "genome_build_inferred": False,
                    "coordinate_basis": "1-based",
                }
                return pd.DataFrame(
                    {"SNP": ["rs1"], "CHR": ["1"], "POS": [100], "A1": ["A"], "A2": ["G"], "Z": [1.0], "N": [1000.0]}
                )

            with mock.patch.object(kernel_munge, "munge_sumstats", side_effect=fake_munge):
                SumstatsMunger().run(
                    MungeConfig(raw_sumstats_file=raw_path, trait_name="trait"),
                    MungeConfig(output_dir=tmpdir / "munged"),
                    GlobalConfig(snp_identifier="chr_pos", genome_build="hg38"),
                )

            metadata = json.loads((tmpdir / "munged" / "sumstats.metadata.json").read_text(encoding="utf-8"))
            self.assertEqual(
                metadata,
                {
                    "format": "ldsc.sumstats.v1",
                    "trait_name": "trait",
                    "config_snapshot": {
                        "snp_identifier": "chr_pos",
                        "genome_build": "hg38",
                        "log_level": "INFO",
                        "fail_on_missing_metadata": False,
                    },
                },
            )
            log_text = (tmpdir / "munged" / "sumstats.log").read_text(encoding="utf-8")
            self.assertIn("Summary-statistics coordinate provenance", log_text)
            self.assertIn('"coordinate_basis": "1-based"', log_text)
            self.assertIn("Summary-statistics liftover report", log_text)
            self.assertEqual(
                set(metadata),
                {"format", "trait_name", "config_snapshot"},
            )

    def test_run_rejects_liftover_request_in_rsid_mode_before_kernel_call(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            raw_path = tmpdir / "raw.tsv"
            raw_path.write_text("SNP A1 A2 P BETA N\nrs1 A G 0.05 0.1 1000\n", encoding="utf-8")

            with mock.patch.object(kernel_munge, "munge_sumstats", side_effect=AssertionError("kernel should not run")):
                with self.assertRaisesRegex(ValueError, "liftover.*chr_pos"):
                    SumstatsMunger().run(
                        MungeConfig(raw_sumstats_file=raw_path, trait_name="trait"),
                        MungeConfig(
                            output_dir=tmpdir / "munged",
                            target_genome_build="hg38",
                            use_hm3_quick_liftover=True,
                        ),
                        GlobalConfig(snp_identifier="rsid"),
                    )

    def test_run_requires_liftover_method_when_target_differs_from_resolved_source(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            raw_path = tmpdir / "raw.tsv"
            raw_path.write_text("CHR POS SNP A1 A2 P BETA N\n1 100 rs1 A G 0.05 0.1 1000\n", encoding="utf-8")

            with self.assertRaisesRegex(ValueError, "liftover method"):
                SumstatsMunger().run(
                    MungeConfig(raw_sumstats_file=raw_path, trait_name="trait"),
                    MungeConfig(output_dir=tmpdir / "munged", target_genome_build="hg38"),
                    GlobalConfig(snp_identifier="chr_pos", genome_build="hg19"),
                )

    def test_sumstats_table_uses_chr_pos_identity_when_configured(self):
        table = sumstats_workflow.SumstatsTable(
            data=pd.DataFrame(
                {
                    "SNP": ["rs_label_1", "1:200"],
                    "CHR": ["1", "1"],
                    "POS": [100, 200],
                    "Z": [1.0, 2.0],
                    "N": [100.0, 100.0],
                }
            ),
            has_alleles=False,
            source_path=None,
            trait_name="trait",
            config_snapshot=GlobalConfig(snp_identifier="chr_pos", genome_build="hg38"),
        )
        metadata = pd.DataFrame({"SNP": ["different_2", "different_1"], "CHR": ["1", "1"], "POS": [200, 100]})

        self.assertEqual(table.snp_identifiers().tolist(), ["1:100", "1:200"])
        self.assertEqual(table.subset_to({"1:200"}).data["SNP"].tolist(), ["1:200"])
        self.assertEqual(table.align_to_metadata(metadata).data["SNP"].tolist(), ["1:200", "rs_label_1"])

    def test_sumstats_table_unknown_config_defaults_to_chr_pos_identity(self):
        table = sumstats_workflow.SumstatsTable(
            data=pd.DataFrame(
                {
                    "SNP": ["rs_label_1", "1:200"],
                    "CHR": ["1", "1"],
                    "POS": [100, 200],
                    "Z": [1.0, 2.0],
                    "N": [100.0, 100.0],
                }
            ),
            has_alleles=False,
            source_path=None,
            trait_name="trait",
            config_snapshot=None,
        )
        metadata = pd.DataFrame({"SNP": ["different_2", "different_1"], "CHR": ["1", "1"], "POS": [200, 100]})

        self.assertEqual(table.snp_identifiers().tolist(), ["1:100", "1:200"])
        self.assertEqual(table.subset_to({"1:200"}).data["SNP"].tolist(), ["1:200"])
        self.assertEqual(table.align_to_metadata(metadata).data["SNP"].tolist(), ["1:200", "rs_label_1"])

    def test_run_accepts_merged_munge_config(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            raw_path = tmpdir / "raw.tsv"
            raw_path.write_text(
                "SNP A1 A2 P OR N\n"
                "rs1 A G 0.05 1.1 1000\n"
                "rs2 C T 0.10 0.9 1000\n",
                encoding="utf-8",
            )
            config = MungeConfig(
                raw_sumstats_file=str(raw_path),
                output_dir=str(tmpdir / "munged"),
                trait_name="trait",
            )

            table = SumstatsMunger().run(config, global_config=GlobalConfig(snp_identifier="rsid"))

            self.assertEqual(table.trait_name, "trait")
            self.assertEqual(table.source_path, str(raw_path))
            self.assertTrue((tmpdir / "munged" / "sumstats.parquet").exists())

    def test_top_level_wrapper_calls_main(self):
        from ldsc import sumstats_munger as munge_sumstats

        with mock.patch.object(munge_sumstats, "main", return_value=11) as patched:
            rc = munge_sumstats.__dict__["main"](["--help"])
        patched.assert_called_once()
        self.assertEqual(rc, 11)

    def test_main_chr_pos_auto_defers_genome_build_inference_to_kernel(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            raw_path = tmpdir / "raw.tsv"
            raw_path.write_text(
                "CHR\tPOS\tSNP\tA1\tA2\tP\tN\n"
                "1\t100\trs1\tA\tG\t0.05\t1000\n",
                encoding="utf-8",
            )

            with mock.patch.object(
                kernel_munge,
                "munge_sumstats",
                return_value=pd.DataFrame(
                    {"SNP": ["rs1"], "CHR": ["1"], "POS": [100], "A1": ["A"], "A2": ["G"], "Z": [1.0], "N": [1000]}
                ),
            ) as patched_munge:
                sumstats_workflow.main(
                    [
                        "--raw-sumstats-file",
                        str(raw_path),
                        "--output-dir",
                        str(tmpdir / "out"),
                        "--snp-identifier",
                        "chr_pos",
                        "--genome-build",
                        "auto",
                    ]
                )

            self.assertEqual(patched_munge.call_args.args[0].genome_build, "auto")

    def test_main_accepts_sumstats_snps_file(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            raw_path = tmpdir / "raw.tsv"
            raw_path.write_text(
                "CHR\tPOS\tSNP\tA1\tA2\tP\tN\n"
                "1\t100\trs1\tA\tG\t0.05\t1000\n",
                encoding="utf-8",
            )
            keep_path = tmpdir / "keep.tsv"
            keep_path.write_text("SNP\nrs1\n", encoding="utf-8")

            with mock.patch.object(
                kernel_munge,
                "munge_sumstats",
                return_value=pd.DataFrame(
                    {"SNP": ["rs1"], "CHR": ["1"], "POS": [100], "A1": ["A"], "A2": ["G"], "Z": [1.0], "N": [1000]}
                ),
            ) as patched_munge:
                sumstats_workflow.main(
                    [
                        "--raw-sumstats-file",
                        str(raw_path),
                        "--output-dir",
                        str(tmpdir / "out"),
                        "--sumstats-snps-file",
                        str(keep_path),
                    ]
                )

            self.assertEqual(patched_munge.call_args.args[0].sumstats_snps, str(keep_path))

    def test_run_creates_output_dir(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            raw_path = tmpdir / "raw.tsv"
            raw_path.write_text(
                "SNP A1 A2 P OR N\n"
                "rs1 A G 0.05 1.0 1000\n",
                encoding="utf-8",
            )
            raw = MungeConfig(raw_sumstats_file=str(raw_path), trait_name="trait")
            output_dir = tmpdir / "nested" / "dir" / "munged"
            config = MungeConfig(output_dir=str(output_dir))

            table = SumstatsMunger().run(raw, config, GlobalConfig(snp_identifier="rsid"))

            self.assertEqual(len(table.data), 1)
            self.assertTrue((output_dir / "sumstats.parquet").exists())
            self.assertTrue((output_dir / "sumstats.log").exists())
            log_text = (output_dir / "sumstats.log").read_text(encoding="utf-8")
            self.assertIn("LDSC munge-sumstats Started", log_text)
            self.assertIn("Munging summary statistics", log_text)

    def test_run_refuses_existing_fixed_outputs_before_kernel_call(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            raw_path = tmpdir / "raw.tsv"
            raw_path.write_text("SNP A1 A2 P OR N\nrs1 A G 0.05 1.0 1000\n", encoding="utf-8")
            output_dir = tmpdir / "munged"
            output_dir.mkdir()
            existing = output_dir / "sumstats.parquet"
            existing.write_text("existing\n", encoding="utf-8")

            with mock.patch.object(kernel_munge, "munge_sumstats", side_effect=AssertionError("kernel should not run")):
                with self.assertRaisesRegex(FileExistsError, "overwrite"):
                    SumstatsMunger().run(
                        MungeConfig(raw_sumstats_file=raw_path, trait_name="trait"),
                        MungeConfig(output_dir=output_dir),
                        GlobalConfig(snp_identifier="rsid"),
                    )

            self.assertEqual(existing.read_text(encoding="utf-8"), "existing\n")
            self.assertFalse((output_dir / "sumstats.log").exists())

    def test_run_allows_existing_fixed_outputs_with_overwrite(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            raw_path = tmpdir / "raw.tsv"
            raw_path.write_text("SNP A1 A2 P OR N\nrs1 A G 0.05 1.0 1000\n", encoding="utf-8")
            output_dir = tmpdir / "munged"
            output_dir.mkdir()
            (output_dir / "sumstats.parquet").write_text("existing\n", encoding="utf-8")
            returned = pd.DataFrame({"SNP": ["rs1"], "N": [1000.0], "Z": [1.5], "A1": ["A"], "A2": ["G"]})

            with mock.patch.object(kernel_munge, "munge_sumstats", return_value=returned) as patched:
                table = SumstatsMunger().run(
                    MungeConfig(raw_sumstats_file=raw_path, trait_name="trait"),
                    MungeConfig(output_dir=output_dir, overwrite=True),
                    GlobalConfig(snp_identifier="rsid"),
                )

            patched.assert_called_once()
            self.assertEqual(table.data["SNP"].tolist(), ["rs1"])

    def test_run_accepts_path_objects_for_input_and_output(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            raw_path = tmpdir / "raw.tsv"
            raw_path.write_text(
                "SNP A1 A2 P OR N\n"
                "rs1 A G 0.05 1.0 1000\n",
                encoding="utf-8",
            )
            raw = MungeConfig(raw_sumstats_file=raw_path, trait_name="trait")
            output_dir = tmpdir / "nested" / "mdd2025"

            table = SumstatsMunger().run(raw, MungeConfig(output_dir=output_dir), GlobalConfig(snp_identifier="rsid"))

            self.assertEqual(len(table.data), 1)
            self.assertEqual(table.source_path, str(raw_path))
            self.assertTrue((output_dir / "sumstats.parquet").exists())
            self.assertTrue((output_dir / "sumstats.log").exists())

    def test_run_resolves_glob_pattern_for_single_sumstats_input(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            raw_path = tmpdir / "trait.raw.tsv"
            raw_path.write_text(
                "SNP A1 A2 P OR N\n"
                "rs1 A G 0.05 1.0 1000\n",
                encoding="utf-8",
            )
            raw = MungeConfig(raw_sumstats_file=str(tmpdir / "trait.*.tsv"), trait_name="trait")
            output_dir = tmpdir / "out" / "munged"

            table = SumstatsMunger().run(raw, MungeConfig(output_dir=output_dir), GlobalConfig(snp_identifier="rsid"))

            self.assertEqual(table.source_path, str(raw_path))
            self.assertTrue((output_dir / "sumstats.parquet").exists())

    def test_kernel_missing_sample_size_error_names_fix_options(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            raw_path = tmpdir / "raw.tsv"
            raw_path.write_text("SNP A1 A2 P OR\nrs1 A G 0.05 1.1\n", encoding="utf-8")
            args = kernel_munge.parser.parse_args(["--sumstats", str(raw_path), "--out", str(tmpdir / "sumstats")])

            with self.assertRaisesRegex(ValueError, "--N.*--N-cas.*--N-con.*N column"):
                kernel_munge.munge_sumstats(args, p=True)

    def test_kernel_error_path_does_not_write_progress_to_stdout(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            raw_path = tmpdir / "raw.tsv"
            raw_path.write_text("SNP A1 A2 P OR\nrs1 A G 0.05 1.1\n", encoding="utf-8")
            args = kernel_munge.parser.parse_args(["--sumstats", str(raw_path), "--out", str(tmpdir / "sumstats")])
            stdout = io.StringIO()

            with contextlib.redirect_stdout(stdout):
                with self.assertRaisesRegex(ValueError, "Could not determine N"):
                    kernel_munge.munge_sumstats(args, p=True)

            self.assertEqual(stdout.getvalue(), "")

    def test_kernel_missing_allele_error_names_no_alleles_fix(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            raw_path = tmpdir / "raw.tsv"
            raw_path.write_text("SNP P OR N\nrs1 0.05 1.1 1000\n", encoding="utf-8")
            args = kernel_munge.parser.parse_args(["--sumstats", str(raw_path), "--out", str(tmpdir / "sumstats")])

            with self.assertRaisesRegex(ValueError, "--a1.*--a2.*--no-alleles"):
                kernel_munge.munge_sumstats(args, p=True)

    def test_kernel_reads_tab_delimited_merge_alleles_with_extra_spaced_columns(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            raw_path = tmpdir / "raw.tsv"
            raw_path.write_text(
                "SNP A1 A2 P BETA N\n"
                "rs1 A G 0.05 0.1 1000\n"
                "rs2 C T 0.10 -0.1 1000\n",
                encoding="utf-8",
            )
            merge_path = tmpdir / "hm3.tsv.gz"
            with gzip.open(merge_path, "wt", encoding="utf-8") as handle:
                handle.write("CHR\thg19_POS\tSNP\tA1\tA2\thm3_alleles\n")
                handle.write("1\t100\trs1\tA\tG\tfrozenset({'A', 'G'})\n")

            args = kernel_munge.parser.parse_args(
                ["--sumstats", str(raw_path), "--out", str(tmpdir / "sumstats"), "--merge-alleles", str(merge_path)]
            )

            table = kernel_munge.munge_sumstats(args, p=True)

            self.assertEqual(table["SNP"].tolist(), ["rs1"])
            self.assertEqual(table["A1"].tolist(), ["A"])
            self.assertEqual(table["A2"].tolist(), ["G"])

    def test_run_restricts_sumstats_snps_file_by_rsid_without_reordering(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            raw_path = tmpdir / "raw.tsv"
            raw_path.write_text(
                "SNP A1 A2 P BETA N\n"
                "rs1 A G 0.05 0.1 1000\n"
                "rs2 C T 0.10 -0.1 1000\n"
                "rs3 G A 0.20 0.2 1000\n",
                encoding="utf-8",
            )
            keep_path = tmpdir / "keep.tsv"
            keep_path.write_text("SNP\nrs3\nrs1\n", encoding="utf-8")

            table = SumstatsMunger().run(
                MungeConfig(raw_sumstats_file=raw_path, trait_name="trait"),
                MungeConfig(output_dir=tmpdir / "munged", sumstats_snps_file=keep_path),
                GlobalConfig(snp_identifier="rsid"),
            )

            self.assertEqual(table.data["SNP"].tolist(), ["rs1", "rs3"])
            log_text = (tmpdir / "munged" / "sumstats.log").read_text(encoding="utf-8")
            self.assertIn("Applying --sumstats-snps-file keep-list", log_text)
            self.assertIn("snp_identifier=rsid", log_text)
            self.assertIn("read 2 keep-list identifiers", log_text)

    def test_run_restricts_sumstats_snps_file_by_chr_pos(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            raw_path = tmpdir / "raw.tsv"
            raw_path.write_text(
                "CHR POS SNP A1 A2 P BETA N\n"
                "1 100 rs1 A G 0.05 0.1 1000\n"
                "1 200 rs2 C T 0.10 -0.1 1000\n",
                encoding="utf-8",
            )
            keep_path = tmpdir / "keep.tsv"
            keep_path.write_text("CHR\tPOS\n1\t200\n", encoding="utf-8")

            table = SumstatsMunger().run(
                MungeConfig(raw_sumstats_file=raw_path, trait_name="trait"),
                MungeConfig(output_dir=tmpdir / "munged", sumstats_snps_file=keep_path),
                GlobalConfig(snp_identifier="chr_pos", genome_build="hg38"),
            )

            self.assertEqual(table.data["SNP"].tolist(), ["rs2"])
            self.assertEqual(table.data["POS"].tolist(), [200])

    def test_run_restricts_sumstats_snps_file_by_chr_pos_logs_missing_coordinate_drops(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            raw_path = tmpdir / "raw.tsv"
            raw_path.write_text(
                "CHR POS SNP A1 A2 P BETA N\n"
                "1 100 rs1 A G 0.05 0.1 1000\n"
                "NA 200 missing_chr C T 0.10 -0.1 1000\n"
                "1 NA missing_pos G A 0.20 0.2 1000\n",
                encoding="utf-8",
            )
            keep_path = tmpdir / "keep.tsv"
            keep_path.write_text("CHR\tPOS\n1\t100\n", encoding="utf-8")

            table = SumstatsMunger().run(
                MungeConfig(raw_sumstats_file=raw_path, trait_name="trait"),
                MungeConfig(output_dir=tmpdir / "munged", sumstats_snps_file=keep_path),
                GlobalConfig(snp_identifier="chr_pos", genome_build="hg38"),
            )

            self.assertEqual(table.data["SNP"].tolist(), ["rs1"])
            log_text = (tmpdir / "munged" / "sumstats.log").read_text(encoding="utf-8")
            self.assertIn("Dropped 2 SNPs with missing CHR/POS", log_text)
            self.assertIn("missing_chr", log_text)

    def test_run_reports_context_when_sumstats_snps_file_removes_everything(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            raw_path = tmpdir / "raw.tsv"
            raw_path.write_text(
                "SNP A1 A2 P BETA N\n"
                "rs1 A G 0.05 0.1 1000\n",
                encoding="utf-8",
            )
            keep_path = tmpdir / "keep.tsv"
            keep_path.write_text("SNP\nrs2\n", encoding="utf-8")

            with self.assertRaisesRegex(
                ValueError,
                "Keep-list file: .*snp_identifier=rsid.*keep-list identifiers=1",
            ):
                SumstatsMunger().run(
                    MungeConfig(raw_sumstats_file=raw_path, trait_name="trait"),
                    MungeConfig(output_dir=tmpdir / "munged", sumstats_snps_file=keep_path),
                    GlobalConfig(snp_identifier="rsid"),
                )

    def test_run_restricts_sumstats_snps_file_by_build_specific_chr_pos_column(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            raw_path = tmpdir / "raw.tsv"
            raw_path.write_text(
                "CHR POS SNP A1 A2 P BETA N\n"
                "1 100 rs1 A G 0.05 0.1 1000\n"
                "1 200 rs2 C T 0.10 -0.1 1000\n",
                encoding="utf-8",
            )
            keep_path = tmpdir / "keep.tsv"
            keep_path.write_text("CHR\thg19_POS\thg38_POS\n1\t200\t999\n", encoding="utf-8")

            table = SumstatsMunger().run(
                MungeConfig(raw_sumstats_file=raw_path, trait_name="trait"),
                MungeConfig(output_dir=tmpdir / "munged", sumstats_snps_file=keep_path),
                GlobalConfig(snp_identifier="chr_pos", genome_build="hg19"),
            )

            self.assertEqual(table.data["SNP"].tolist(), ["rs2"])

    def test_run_restricts_sumstats_snps_file_after_auto_coordinate_normalization(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            raw_path = tmpdir / "raw.tsv"
            raw_path.write_text(
                "CHR POS SNP A1 A2 P BETA N\n"
                "1 99 rs1 A G 0.05 0.1 1000\n"
                "1 199 rs2 C T 0.10 -0.1 1000\n",
                encoding="utf-8",
            )
            keep_path = tmpdir / "keep.tsv"
            keep_path.write_text("CHR\tPOS\n1\t200\n", encoding="utf-8")

            inference = mock.Mock(
                genome_build="hg19",
                coordinate_basis="0-based",
                inspected_snp_count=2,
                match_counts={},
                match_fractions={},
                summary_message="mock inferred hg19 0-based",
            )

            def fake_resolve_chr_pos_table(frame, *, context, logger=None, reference_table=None):
                normalized = frame.copy()
                normalized["CHR"] = normalized["CHR"].astype(str)
                normalized["POS"] = normalized["POS"].astype(int) + 1
                return normalized, inference

            with mock.patch.object(kernel_munge, "resolve_chr_pos_table", side_effect=fake_resolve_chr_pos_table):
                table = SumstatsMunger().run(
                    MungeConfig(raw_sumstats_file=raw_path, trait_name="trait"),
                    MungeConfig(output_dir=tmpdir / "munged", sumstats_snps_file=keep_path),
                    GlobalConfig(snp_identifier="chr_pos", genome_build="auto"),
                )

            self.assertEqual(table.data["SNP"].tolist(), ["rs2"])
            self.assertEqual(table.data["POS"].tolist(), [200])
            self.assertEqual(table.config_snapshot.genome_build, "hg19")

    def test_run_treats_sumstats_snps_file_with_alleles_as_keep_list_only(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            raw_path = tmpdir / "raw.tsv"
            raw_path.write_text(
                "SNP A1 A2 P BETA N\n"
                "rs1 A G 0.05 0.1 1000\n"
                "rs2 C T 0.10 -0.1 1000\n",
                encoding="utf-8",
            )
            keep_path = tmpdir / "keep.tsv"
            keep_path.write_text("SNP\tA1\tA2\nrs1\tT\tC\n", encoding="utf-8")

            table = SumstatsMunger().run(
                MungeConfig(raw_sumstats_file=raw_path, trait_name="trait"),
                MungeConfig(output_dir=tmpdir / "munged", sumstats_snps_file=keep_path),
                GlobalConfig(snp_identifier="rsid"),
            )

            self.assertEqual(table.data["SNP"].tolist(), ["rs1"])
            self.assertEqual(table.data["A1"].tolist(), ["A"])
            self.assertEqual(table.data["A2"].tolist(), ["G"])

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
                MungeConfig(raw_sumstats_file=raw_path, trait_name="trait"),
                MungeConfig(output_dir=tmpdir / "munged"),
                GlobalConfig(snp_identifier="rsid"),
            )

            self.assertEqual(table.data["SNP"].tolist(), ["rs1"])
            self.assertEqual(table.data["N"].tolist(), [1000.0])
            self.assertEqual(table.data["CHR"].tolist(), [1])
            self.assertEqual(table.data["POS"].tolist(), [123])

    def test_run_skips_leading_double_hash_sumstats_metadata_lines(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            raw_path = tmpdir / "raw.tsv"
            raw_path.write_text(
                '##EA="SNP effect allele for ln(OR))"\n'
                '##NEA="SNP non-effect allele"\n'
                '##BETA="ln(Odds Ratio) effect of EA"\n'
                '##PVAL="P-value, uncorrected"\n'
                "## \n"
                "#CHROM\tPOS\tID\tEA\tNEA\tBETA\tPVAL\tNCAS\tNCON\n"
                "1\t753541\trs2073813\tA\tG\t0.004\t0.4614\t310128\t1035355\n",
                encoding="utf-8",
            )

            munger = SumstatsMunger()
            table = munger.run(
                MungeConfig(raw_sumstats_file=raw_path, trait_name="trait"),
                MungeConfig(output_dir=tmpdir / "munged"),
                GlobalConfig(snp_identifier="chr_pos", genome_build="hg38"),
            )

            self.assertEqual(table.data["SNP"].tolist(), ["rs2073813"])
            self.assertEqual(table.data["CHR"].tolist(), ["1"])
            self.assertEqual(table.data["POS"].tolist(), [753541])
            self.assertEqual(munger.build_run_summary(table).n_input_rows, 1)

    def test_run_accepts_chrom_and_bp_aliases_for_coordinates(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            raw_path = tmpdir / "raw.tsv"
            raw_path.write_text(
                "CHROM\tBP\tID\tEA\tNEA\tPVAL\tOR\tN\n"
                "chr1\t123\trs1\tA\tG\t0.05\t1.0\t1000\n",
                encoding="utf-8",
            )

            table = SumstatsMunger().run(
                MungeConfig(raw_sumstats_file=raw_path, trait_name="trait"),
                MungeConfig(output_dir=tmpdir / "munged"),
                GlobalConfig(snp_identifier="chr_pos", genome_build="hg38"),
            )

            self.assertEqual(table.data["CHR"].tolist(), ["1"])
            self.assertEqual(table.data["POS"].tolist(), [123])
            output = pd.read_parquet(tmpdir / "munged" / "sumstats.parquet")
            self.assertEqual(output.columns.tolist(), ["SNP", "CHR", "POS", "A1", "A2", "Z", "N"])

    def test_run_accepts_explicit_chr_and_pos_column_hints(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            raw_path = tmpdir / "raw.tsv"
            raw_path.write_text(
                "chromosome_name\tbase_pair\tvariant\tA1\tA2\tP\tOR\tN\n"
                "2\t200\trs2\tC\tT\t0.05\t1.0\t1000\n",
                encoding="utf-8",
            )

            table = SumstatsMunger().run(
                MungeConfig(
                    raw_sumstats_file=raw_path,
                    trait_name="trait",
                    column_hints={"chr": "chromosome_name", "pos": "base_pair", "snp": "variant"},
                ),
                MungeConfig(output_dir=tmpdir / "munged"),
                GlobalConfig(snp_identifier="chr_pos", genome_build="hg38"),
            )

            self.assertEqual(table.data["SNP"].tolist(), ["rs2"])
            self.assertEqual(table.data["CHR"].tolist(), ["2"])
            self.assertEqual(table.data["POS"].tolist(), [200])

    def test_run_rejects_ambiguous_chromosome_aliases(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            raw_path = tmpdir / "raw.tsv"
            raw_path.write_text(
                "#CHROM\tCHROM\tPOS\tID\tEA\tNEA\tPVAL\tOR\tN\n"
                "1\t1\t123\trs1\tA\tG\t0.05\t1.1\t1000\n",
                encoding="utf-8",
            )

            with self.assertRaisesRegex(ValueError, "CHR"):
                SumstatsMunger().run(
                    MungeConfig(raw_sumstats_file=raw_path, trait_name="trait"),
                    MungeConfig(output_dir=tmpdir / "munged"),
                    GlobalConfig(snp_identifier="chr_pos", genome_build="hg38"),
                )

    def test_load_sumstats_reads_metadata_sidecar_config_snapshot(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            raw_path = tmpdir / "raw.tsv"
            raw_path.write_text(
                "#CHROM\tPOS\tID\tEA\tNEA\tPVAL\tOR\tN\n"
                "1\t123\trs1\tA\tG\t0.05\t1.0\t1000\n",
                encoding="utf-8",
            )
            output_dir = tmpdir / "munged"

            SumstatsMunger().run(
                MungeConfig(raw_sumstats_file=raw_path, trait_name="trait"),
                MungeConfig(output_dir=output_dir),
                GlobalConfig(snp_identifier="chr_pos", genome_build="hg38"),
            )

            metadata = json.loads((output_dir / "sumstats.metadata.json").read_text(encoding="utf-8"))
            self.assertEqual(metadata["format"], "ldsc.sumstats.v1")
            self.assertNotIn("snp_identifier", metadata)
            self.assertNotIn("genome_build", metadata)
            self.assertEqual(metadata["config_snapshot"]["snp_identifier"], "chr_pos")
            self.assertEqual(metadata["config_snapshot"]["genome_build"], "hg38")
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                table = ldsc.load_sumstats(output_dir / "sumstats.parquet", trait_name="trait")
            self.assertEqual(table.config_snapshot, GlobalConfig(snp_identifier="chr_pos", genome_build="hg38"))
            self.assertFalse(any("cannot recover the GlobalConfig" in str(item.message) for item in caught))
