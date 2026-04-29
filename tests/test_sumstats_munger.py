from pathlib import Path
import contextlib
import gzip
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

    def test_kernel_p_to_z_matches_legacy_direction_convention(self):
        z = kernel_munge.p_to_z(pd.Series([0.1, 0.1, 0.1]), pd.Series([1, 2, 3]))
        np.testing.assert_allclose(np.asarray(z), [1.644854, 1.644854, 1.644854], atol=1e-5)

    def test_kernel_filters_reject_invalid_p_info_frq_and_alleles(self):
        log = mock.Mock(log=mock.Mock())
        args = kernel_munge.parser.parse_args("")
        assert_series_equal(
            kernel_munge.filter_pvals(pd.Series([0, 0.1, 1, 2]), log, args),
            pd.Series([False, True, True, False]),
        )
        assert_series_equal(
            kernel_munge.filter_info(pd.Series([0.8, 1.0, 1.0]), log, args),
            pd.Series([False, True, True]),
        )
        assert_series_equal(
            kernel_munge.filter_frq(pd.Series([-1, 0, 0.005, 0.4, 0.6, 0.999, 1, 2]), log, args),
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
            raw = MungeConfig(sumstats_file=str(raw_path), trait_name="trait")
            config = MungeConfig(output_dir=str(tmpdir / "munged"))
            munger = SumstatsMunger()
            table = munger.run(raw, config, GlobalConfig(snp_identifier="rsid"))
            self.assertEqual(table.trait_name, "trait")
            self.assertEqual(table.config_snapshot, GlobalConfig(snp_identifier="rsid"))
            self.assertTrue((tmpdir / "munged" / "sumstats.sumstats.gz").exists())
            with gzip.open(tmpdir / "munged" / "sumstats.sumstats.gz", "rt", encoding="utf-8") as handle:
                output = pd.read_csv(handle, sep="\t")
            self.assertEqual(output.columns.tolist(), ["SNP", "CHR", "POS", "A1", "A2", "Z", "N"])
            self.assertTrue(output["CHR"].isna().all())
            self.assertTrue(output["POS"].isna().all())
            self.assertTrue((tmpdir / "munged" / "sumstats.metadata.json").exists())
            summary = munger.build_run_summary(table)
            self.assertEqual(summary.n_retained_rows, 2)
            self.assertIn("sumstats_gz", summary.output_paths)
            self.assertIn("metadata_json", summary.output_paths)

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
                sumstats_file=str(raw_path),
                output_dir=str(tmpdir / "munged"),
                trait_name="trait",
            )

            table = SumstatsMunger().run(config, global_config=GlobalConfig(snp_identifier="rsid"))

            self.assertEqual(table.trait_name, "trait")
            self.assertEqual(table.source_path, str(raw_path))
            self.assertTrue((tmpdir / "munged" / "sumstats.sumstats.gz").exists())

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
                return_value=pd.DataFrame({"SNP": ["rs1"], "CHR": ["1"], "POS": [100], "Z": [1.0], "N": [1000]}),
            ) as patched_munge:
                sumstats_workflow.main(
                    [
                        "--sumstats-file",
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
                return_value=pd.DataFrame({"SNP": ["rs1"], "CHR": ["1"], "POS": [100], "Z": [1.0], "N": [1000]}),
            ) as patched_munge:
                sumstats_workflow.main(
                    [
                        "--sumstats-file",
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
            raw = MungeConfig(sumstats_file=str(raw_path), trait_name="trait")
            output_dir = tmpdir / "nested" / "dir" / "munged"
            config = MungeConfig(output_dir=str(output_dir))

            table = SumstatsMunger().run(raw, config, GlobalConfig(snp_identifier="rsid"))

            self.assertEqual(len(table.data), 1)
            self.assertTrue((output_dir / "sumstats.sumstats.gz").exists())
            self.assertTrue((output_dir / "sumstats.log").exists())

    def test_run_refuses_existing_fixed_outputs_before_kernel_call(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            raw_path = tmpdir / "raw.tsv"
            raw_path.write_text("SNP A1 A2 P OR N\nrs1 A G 0.05 1.0 1000\n", encoding="utf-8")
            output_dir = tmpdir / "munged"
            output_dir.mkdir()
            existing = output_dir / "sumstats.sumstats.gz"
            existing.write_text("existing\n", encoding="utf-8")

            with mock.patch.object(kernel_munge, "munge_sumstats", side_effect=AssertionError("kernel should not run")):
                with self.assertRaisesRegex(FileExistsError, "overwrite"):
                    SumstatsMunger().run(
                        MungeConfig(sumstats_file=raw_path, trait_name="trait"),
                        MungeConfig(output_dir=output_dir),
                        GlobalConfig(snp_identifier="rsid"),
                    )

            self.assertEqual(existing.read_text(encoding="utf-8"), "existing\n")

    def test_run_allows_existing_fixed_outputs_with_overwrite(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            raw_path = tmpdir / "raw.tsv"
            raw_path.write_text("SNP A1 A2 P OR N\nrs1 A G 0.05 1.0 1000\n", encoding="utf-8")
            output_dir = tmpdir / "munged"
            output_dir.mkdir()
            (output_dir / "sumstats.sumstats.gz").write_text("existing\n", encoding="utf-8")
            returned = pd.DataFrame({"SNP": ["rs1"], "N": [1000.0], "Z": [1.5], "A1": ["A"], "A2": ["G"]})

            with mock.patch.object(kernel_munge, "munge_sumstats", return_value=returned) as patched:
                table = SumstatsMunger().run(
                    MungeConfig(sumstats_file=raw_path, trait_name="trait"),
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
            raw = MungeConfig(sumstats_file=raw_path, trait_name="trait")
            output_dir = tmpdir / "nested" / "mdd2025"

            table = SumstatsMunger().run(raw, MungeConfig(output_dir=output_dir), GlobalConfig(snp_identifier="rsid"))

            self.assertEqual(len(table.data), 1)
            self.assertEqual(table.source_path, str(raw_path))
            self.assertTrue((output_dir / "sumstats.sumstats.gz").exists())
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
            raw = MungeConfig(sumstats_file=str(tmpdir / "trait.*.tsv"), trait_name="trait")
            output_dir = tmpdir / "out" / "munged"

            table = SumstatsMunger().run(raw, MungeConfig(output_dir=output_dir), GlobalConfig(snp_identifier="rsid"))

            self.assertEqual(table.source_path, str(raw_path))
            self.assertTrue((output_dir / "sumstats.sumstats.gz").exists())

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
                MungeConfig(sumstats_file=raw_path, trait_name="trait"),
                MungeConfig(output_dir=tmpdir / "munged", sumstats_snps_file=keep_path),
                GlobalConfig(snp_identifier="rsid"),
            )

            self.assertEqual(table.data["SNP"].tolist(), ["rs1", "rs3"])

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
                MungeConfig(sumstats_file=raw_path, trait_name="trait"),
                MungeConfig(output_dir=tmpdir / "munged", sumstats_snps_file=keep_path),
                GlobalConfig(snp_identifier="chr_pos", genome_build="hg38"),
            )

            self.assertEqual(table.data["SNP"].tolist(), ["rs2"])
            self.assertEqual(table.data["POS"].tolist(), [200])

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
                MungeConfig(sumstats_file=raw_path, trait_name="trait"),
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
                    MungeConfig(sumstats_file=raw_path, trait_name="trait"),
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
                MungeConfig(sumstats_file=raw_path, trait_name="trait"),
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
                MungeConfig(sumstats_file=raw_path, trait_name="trait"),
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
                MungeConfig(sumstats_file=raw_path, trait_name="trait"),
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
                MungeConfig(sumstats_file=raw_path, trait_name="trait"),
                MungeConfig(output_dir=tmpdir / "munged"),
                GlobalConfig(snp_identifier="chr_pos", genome_build="hg38"),
            )

            self.assertEqual(table.data["CHR"].tolist(), ["1"])
            self.assertEqual(table.data["POS"].tolist(), [123])
            with gzip.open(tmpdir / "munged" / "sumstats.sumstats.gz", "rt", encoding="utf-8") as handle:
                output = pd.read_csv(handle, sep="\t")
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
                    sumstats_file=raw_path,
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
                    MungeConfig(sumstats_file=raw_path, trait_name="trait"),
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
                MungeConfig(sumstats_file=raw_path, trait_name="trait"),
                MungeConfig(output_dir=output_dir),
                GlobalConfig(snp_identifier="chr_pos", genome_build="hg38"),
            )

            metadata = json.loads((output_dir / "sumstats.metadata.json").read_text(encoding="utf-8"))
            self.assertEqual(metadata["format"], "ldsc.sumstats.v1")
            self.assertEqual(metadata["snp_identifier"], "chr_pos")
            self.assertEqual(metadata["genome_build"], "hg38")
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                table = ldsc.load_sumstats(output_dir / "sumstats.sumstats.gz", trait_name="trait")
            self.assertEqual(table.config_snapshot, GlobalConfig(snp_identifier="chr_pos", genome_build="hg38"))
            self.assertFalse(any("cannot recover the GlobalConfig" in str(item.message) for item in caught))
