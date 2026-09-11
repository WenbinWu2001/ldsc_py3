from pathlib import Path
import contextlib
import gzip
import importlib.util
import io
import tempfile
import unittest
import warnings
from unittest import mock

import numpy as np
import pandas as pd
from pandas.testing import assert_series_equal


_HAS_PYARROW = importlib.util.find_spec("pyarrow") is not None

from ldsc.config import GlobalConfig, MungeConfig

try:
    import ldsc
    from ldsc import sumstats_munger as sumstats_workflow
    from ldsc._kernel import sumstats_munger as kernel_munge
    from ldsc import _sumstats_input as munge_input
    from ldsc._kernel.snp_identity import IDENTITY_DROP_COLUMNS
    from ldsc.sumstats_munger import SumstatsMunger
except ImportError:
    ldsc = None
    kernel_munge = None
    sumstats_workflow = None
    IDENTITY_DROP_COLUMNS = None
    SumstatsMunger = None


@unittest.skipIf(SumstatsMunger is None, "sumstats_munger module is not available")
class SumstatsMungerTest(unittest.TestCase):

    DROPPED_SNP_DTYPES = {
        "CHR": "string",
        "SNP": "string",
        "source_pos": "Int64",
        "target_pos": "Int64",
        "reason": "string",
        "base_key": "string",
        "identity_key": "string",
        "allele_set": "string",
        "stage": "string",
    }

    def _read_dropped_snps_sidecar(self, path: Path) -> pd.DataFrame:
        return pd.read_csv(
            path,
            sep="\t",
            compression="gzip",
            dtype=self.DROPPED_SNP_DTYPES,
        )

    def _write_raw_sumstats(self, path: Path) -> None:
        path.write_text("SNP A1 A2 P BETA N\nrs1 A G 0.05 0.1 1000\n", encoding="utf-8")

    def _fake_munged_frame(self) -> pd.DataFrame:
        return pd.DataFrame({"SNP": ["rs1"], "A1": ["A"], "A2": ["G"], "Z": [1.0], "N": [1000.0]})

    def _write_footer_parquet(
        self,
        path: Path,
        frame: pd.DataFrame,
        *,
        snp_identifier: str = "rsid",
        genome_build: str | None = None,
        trait_name: str | None = "trait",
    ) -> None:
        """Write a self-describing sumstats parquet with identity footer metadata."""
        footer = sumstats_workflow._sumstats_footer_metadata(
            GlobalConfig(snp_identifier=snp_identifier, genome_build=genome_build), trait_name
        )
        sumstats_workflow._write_sumstats_outputs(
            frame, output_files={"parquet": str(path)}, output_format="parquet", footer_metadata=footer
        )

    def test_build_parser_defaults_chunksize_to_one_million_rows(self):
        parser = sumstats_workflow.build_parser()
        self.assertEqual(parser.get_default("chunksize"), 1_000_000)

    def test_build_parser_defaults_output_format_to_parquet(self):
        parser = sumstats_workflow.build_parser()
        self.assertEqual(parser.get_default("output_format"), "parquet")
        with self.assertRaises(SystemExit):
            parser.parse_args(["--raw-sumstats-file", "raw.tsv", "--output-dir", "out", "--output-format", "csv"])

    def test_build_parser_defaults_source_genome_build_to_auto_and_rejects_old_build_flags(self):
        parser = sumstats_workflow.build_parser()
        self.assertEqual(parser.get_default("source_genome_build"), "auto")

        args = parser.parse_args(["--raw-sumstats-file", "raw.tsv", "--output-dir", "out"])

        self.assertEqual(args.source_genome_build, "auto")
        with self.assertRaises(SystemExit):
            parser.parse_args(["--raw-sumstats-file", "raw.tsv", "--output-dir", "out", "--genome-build", "auto"])
        with self.assertRaises(SystemExit):
            parser.parse_args(["--raw-sumstats-file", "raw.tsv", "--output-dir", "out", "--target-genome-build", "hg38"])

    def test_build_parser_defaults_log_level_to_info(self):
        parser = sumstats_workflow.build_parser()
        self.assertEqual(parser.get_default("log_level"), "INFO")
        args = parser.parse_args(["--raw-sumstats-file", "raw.tsv", "--output-dir", "out", "--log-level", "DEBUG"])
        self.assertEqual(args.log_level, "DEBUG")

    def test_build_parser_accepts_trait_name(self):
        parser = sumstats_workflow.build_parser()

        args = parser.parse_args(["--raw-sumstats-file", "raw.tsv", "--output-dir", "out", "--trait-name", "MDD"])

        self.assertEqual(args.trait_name, "MDD")

    def test_build_parser_accepts_chain_override_and_defaults_to_hm3(self):
        args = sumstats_workflow.build_parser().parse_args([
            "--raw-sumstats-file", "raw.tsv", "--output-dir", "out",
            "--output-genome-build", "hg38", "--liftover-chain-file", "hg19ToHg38.over.chain",
        ])
        self.assertEqual(args.output_genome_build, "hg38")
        self.assertEqual(args.liftover_chain_file, "hg19ToHg38.over.chain")
        self.assertFalse(args.no_snp_restriction)
        self.assertFalse(hasattr(args, "use_hm3_quick_liftover"))

    def test_build_parser_accepts_no_snp_restriction(self):
        args = sumstats_workflow.build_parser().parse_args([
            "--raw-sumstats-file", "raw.tsv", "--output-dir", "out", "--no-snp-restriction",
        ])
        self.assertTrue(args.no_snp_restriction)

    def test_build_parser_defaults_sumstats_format_to_auto(self):
        parser = sumstats_workflow.build_parser()

        args = parser.parse_args(["--raw-sumstats-file", "raw.tsv", "--output-dir", "out"])

        self.assertEqual(args.sumstats_format, "auto")

    def test_build_parser_rejects_vcf_format_name(self):
        parser = sumstats_workflow.build_parser()

        args = parser.parse_args(["--raw-sumstats-file", "raw.tsv", "--output-dir", "out", "--input-format", "plain"])

        self.assertEqual(args.sumstats_format, "plain")
        with self.assertRaises(SystemExit):
            parser.parse_args(["--raw-sumstats-file", "raw.tsv", "--output-dir", "out", "--input-format", "vcf"])

    def test_build_parser_requires_output_dir_for_infer_only(self):
        parser = sumstats_workflow.build_parser()

        with self.assertRaises(SystemExit):
            parser.parse_args(["--raw-sumstats-file", "raw.tsv", "--infer-only"])
        args = parser.parse_args(
            ["--raw-sumstats-file", "raw.tsv", "--output-dir", "unused", "--infer-only"]
        )

        self.assertTrue(args.infer_only)
        self.assertEqual(args.output_dir, "unused")

    def test_build_parser_selects_daner_via_format_and_rejects_legacy_flags(self):
        parser = sumstats_workflow.build_parser()

        format_old_args = parser.parse_args(["--raw-sumstats-file", "raw.tsv", "--output-dir", "out", "--input-format", "daner-old"])
        format_new_args = parser.parse_args(["--raw-sumstats-file", "raw.tsv", "--output-dir", "out", "--input-format", "daner-new"])

        self.assertEqual(format_old_args.sumstats_format, "daner-old")
        self.assertEqual(format_new_args.sumstats_format, "daner-new")
        # The legacy boolean DANER flags are removed; --input-format is the only selector.
        for legacy in ("--daner-old", "--daner-new", "--daner", "--daner-n"):
            with self.assertRaises(SystemExit):
                parser.parse_args(["--raw-sumstats-file", "raw.tsv", "--output-dir", "out", legacy])

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
            result = sumstats_workflow.main(["--no-snp-restriction", "--raw-sumstats-file", "raw.tsv", "--output-dir", "out"])

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
                "--input-format",
                "daner-new",
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

        table = mock.Mock(provenance={"coordinate_provenance": mock.sentinel.metadata})
        with (
            mock.patch.object(SumstatsMunger, "run", return_value=table) as patched,
            mock.patch.object(SumstatsMunger, "build_run_summary", return_value=mock.sentinel.summary),
            mock.patch.object(sumstats_workflow, "_render_munge_summary", return_value="summary"),
            contextlib.redirect_stdout(io.StringIO()),
        ):
            result = sumstats_workflow.run_munge_sumstats_from_args(args)

        self.assertIs(result, table)
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
        self.assertEqual(run_config.sumstats_format, "daner-new")
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
                "--source-genome-build",
                "hg19",
                "--output-genome-build",
                "hg38",
            ]
        )

        table = mock.Mock(provenance={"coordinate_provenance": mock.sentinel.metadata})
        with (
            mock.patch.object(SumstatsMunger, "run", return_value=table) as patched,
            mock.patch.object(SumstatsMunger, "build_run_summary", return_value=mock.sentinel.summary),
            mock.patch.object(sumstats_workflow, "_render_munge_summary", return_value="summary"),
            contextlib.redirect_stdout(io.StringIO()),
        ):
            result = sumstats_workflow.run_munge_sumstats_from_args(args)

        self.assertIs(result, table)
        raw_config, run_config, global_config = patched.call_args.args
        self.assertEqual(raw_config.raw_sumstats_file, "raw.tsv")
        self.assertEqual(run_config.source_genome_build, "hg19")
        self.assertEqual(run_config.output_genome_build, "hg38")
        self.assertFalse(run_config.no_snp_restriction)
        self.assertIsNone(run_config.liftover_chain_file)
        self.assertEqual(global_config, GlobalConfig(snp_identifier="chr_pos", genome_build="hg19"))

    def test_munge_config_validates_liftover_configuration(self):
        config = MungeConfig(
            output_dir="out", source_genome_build="GRCh37", output_genome_build="GRCh38",
            liftover_chain_file=Path("liftover") / "hg19ToHg38.over.chain",
        )
        self.assertEqual(config.source_genome_build, "hg19")
        self.assertEqual(config.output_genome_build, "hg38")
        self.assertEqual(config.liftover_chain_file, "liftover/hg19ToHg38.over.chain")
        with self.assertRaisesRegex(ldsc.LDSCConfigError, "source_genome_build"):
            MungeConfig(source_genome_build=None, output_genome_build="hg38")
        for options in ({"use_hm3_snps": True}, {"use_hm3_quick_liftover": True}):
            with self.assertRaises(TypeError):
                MungeConfig(**options)

    def test_workflow_reports_actionable_signed_sumstats_format_error(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            raw_path = Path(tmpdir) / "raw.tsv"
            self._write_raw_sumstats(raw_path)
            with self.assertRaisesRegex(ldsc.LDSCUsageError, "could not parse --signed-sumstats='BETA'.*BETA,0"):
                SumstatsMunger().run(
                    MungeConfig(no_snp_restriction=True, raw_sumstats_file=raw_path, output_dir=Path(tmpdir) / "out", signed_sumstats_spec="BETA"),
                    global_config=GlobalConfig(snp_identifier="rsid"),
                )

    def test_kernel_p_to_z_matches_legacy_direction_convention(self):
        z = kernel_munge.p_to_z(pd.Series([0.1, 0.1, 0.1]), pd.Series([1, 2, 3]))
        np.testing.assert_allclose(np.asarray(z), [1.644854, 1.644854, 1.644854], atol=1e-5)

    def test_kernel_filters_reject_invalid_p_info_frq_and_alleles(self):
        args = kernel_munge.MungeQC()
        assert_series_equal(
            kernel_munge.filter_pvals(pd.Series([0, 0.1, 1, 2])),
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

    def test_kernel_filter_frq_is_inclusive_at_maf_min(self):
        args = kernel_munge.MungeQC()
        # Folded MAF exactly at the floor is kept (inclusive >=); just below is dropped.
        assert_series_equal(
            kernel_munge.filter_frq(pd.Series([0.0099, 0.01, 0.9901]), args),
            pd.Series([False, True, False]),
        )

    def test_run_filters_sumstats_snps_file_by_rsid_before_process_n(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            raw_path = tmpdir / "raw.tsv"
            raw_path.write_text(
                "SNP A1 A2 P BETA N\n"
                "rs_keep A G 0.05 0.1 1000\n"
                "rs_drop C T 0.10 -0.1 1000\n",
                encoding="utf-8",
            )
            keep_path = tmpdir / "keep.tsv"
            keep_path.write_text("SNP\nrs_keep\n", encoding="utf-8")
            original_process_n = kernel_munge.process_n

            def assert_restricted_before_process_n(dat, args):
                self.assertEqual(dat["SNP"].tolist(), ["rs_keep"])
                return original_process_n(dat, args)

            with mock.patch.object(kernel_munge, "process_n", side_effect=assert_restricted_before_process_n):
                table = SumstatsMunger().run(
                    MungeConfig(no_snp_restriction=True, raw_sumstats_file=raw_path, trait_name="trait"),
                    MungeConfig(output_dir=tmpdir / "munged", sumstats_snps_file=keep_path),
                    GlobalConfig(snp_identifier="rsid"),
                )

            self.assertEqual(table.data["SNP"].tolist(), ["rs_keep"])

    def test_run_filters_sumstats_snps_file_by_chr_pos_before_process_n_without_string_keys(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            raw_path = tmpdir / "raw.tsv"
            raw_path.write_text(
                "CHR POS SNP A1 A2 P BETA N\n"
                "1 10 label_a A G 0.05 0.0 1000\n"
                "1 11 label_b C T 0.10 -0.1 1000\n"
                "2 20 label_c G A 0.20 0.0 1000\n",
                encoding="utf-8",
            )
            restrict_path = Path(tmpdir) / "hm3.tsv"
            restrict_path.write_text(
                "SNP\tCHR\thg19_POS\thg38_POS\n"
                "rs1\t1\t10\t100\n"
                "rs2\t2\t20\t200\n",
                encoding="utf-8",
            )
            original_process_n = kernel_munge.process_n

            def assert_restricted_before_process_n(dat, args):
                self.assertEqual(dat["SNP"].tolist(), ["label_a", "label_c"])
                return original_process_n(dat, args)

            original_packed_keys = kernel_munge.build_packed_chr_pos_series
            packed_key_calls = []

            def assert_packed_keys(chr_values, pos_values, **kwargs):
                keys = original_packed_keys(chr_values, pos_values, **kwargs)
                self.assertEqual(keys.dtype, np.dtype("uint64"))
                packed_key_calls.append(keys)
                return keys

            for mode in ("chr_pos", "chr_pos_allele_aware"):
                with self.subTest(mode=mode):
                    packed_key_calls.clear()
                    with mock.patch.object(kernel_munge, "build_packed_chr_pos_series", side_effect=assert_packed_keys), \
                         mock.patch.object(kernel_munge, "process_n", side_effect=assert_restricted_before_process_n):
                        table = SumstatsMunger().run(
                            MungeConfig(no_snp_restriction=True, raw_sumstats_file=raw_path, trait_name="trait"),
                            MungeConfig(
                                output_dir=tmpdir / f"munged_{mode}",
                                sumstats_snps_file=restrict_path,
                                source_genome_build="hg19",
                                output_genome_build="hg19",
                            ),
                            GlobalConfig(snp_identifier=mode, genome_build="hg19"),
                        )

                    self.assertEqual(table.data["SNP"].tolist(), ["label_a", "label_c"])
                    self.assertEqual(table.data["POS"].tolist(), [10, 20])
                    self.assertTrue(packed_key_calls)

    def test_prepare_sumstats_restriction_uses_rich_reader_for_allele_aware_modes(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            keep_path = tmpdir / "keep.tsv"
            keep_path.write_text("SNP\tA1\tA2\nrs1\tA\tC\n", encoding="utf-8")
            restriction = kernel_munge.prepare_sumstats_restriction(str(keep_path), "rsid_allele_aware", None)

            self.assertEqual(restriction.identity_keys.match_kind, "identity")
            self.assertEqual(restriction.identity_keys.keys, {"rs1:A:C"})

    def test_run_filters_sumstats_snps_file_by_identity_in_allele_aware_mode_before_process_n(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            raw_path = tmpdir / "raw.tsv"
            raw_path.write_text(
                "SNP A1 A2 P BETA N\n"
                "rs1 A C 0.05 0.1 1000\n"
                "rs1 A G 0.10 -0.1 1000\n",
                encoding="utf-8",
            )
            keep_path = tmpdir / "keep.tsv"
            keep_path.write_text("SNP\tA1\tA2\nrs1\tA\tC\n", encoding="utf-8")
            original_process_n = kernel_munge.process_n

            def assert_restricted_before_process_n(dat, args):
                self.assertEqual(dat["SNP"].tolist(), ["rs1"])
                self.assertEqual(dat["A1"].tolist(), ["A"])
                self.assertEqual(dat["A2"].tolist(), ["C"])
                return original_process_n(dat, args)

            with mock.patch.object(kernel_munge, "process_n", side_effect=assert_restricted_before_process_n):
                table = SumstatsMunger().run(
                    MungeConfig(no_snp_restriction=True, raw_sumstats_file=raw_path, trait_name="trait"),
                    MungeConfig(output_dir=tmpdir / "munged", sumstats_snps_file=keep_path),
                    GlobalConfig(snp_identifier="rsid_allele_aware"),
                )

            self.assertEqual(table.data["SNP"].tolist(), ["rs1"])
            self.assertEqual(table.data["A1"].tolist(), ["A"])
            self.assertEqual(table.data["A2"].tolist(), ["C"])

    def test_run_resolves_chr_pos_auto_before_chunk_parsing(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            raw_path = tmpdir / "raw.tsv"
            raw_path.write_text(
                "CHR POS SNP A1 A2 P BETA N\n"
                "1 99 rs1 A G 0.05 0.1 1000\n",
                encoding="utf-8",
            )

            inference = mock.Mock(
                genome_build="hg19",
                coordinate_basis="0-based",
                inspected_snp_count=1,
                match_counts={},
                match_fractions={},
                summary_message="mock inferred hg19 0-based",
            )

            def fake_resolve_chr_pos_table(frame, *, context, logger=None, reference_table=None, coordinate_policy="drop"):
                normalized = frame.copy()
                normalized["CHR"] = normalized["CHR"].astype(str)
                normalized["POS"] = pd.to_numeric(normalized["POS"], errors="raise").astype(int) + 1
                return normalized, inference

            original_parse = kernel_munge.parse_dat

            def assert_auto_resolved_before_parse(dat_gen, request):
                self.assertEqual(request.genome_build, "hg19")
                self.assertTrue(request.coordinate_metadata["genome_build_inferred"])
                self.assertEqual(request.coordinate_basis, "0-based")
                return original_parse(dat_gen, request)

            with mock.patch.object(munge_input, "resolve_chr_pos_table", side_effect=fake_resolve_chr_pos_table), \
                 mock.patch.object(kernel_munge, "parse_dat", side_effect=assert_auto_resolved_before_parse):
                table = SumstatsMunger().run(
                    MungeConfig(no_snp_restriction=True, raw_sumstats_file=raw_path, trait_name="trait"),
                    MungeConfig(no_snp_restriction=True, output_dir=tmpdir / "munged", source_genome_build="auto", output_genome_build="hg19"),
                    GlobalConfig(snp_identifier="chr_pos", genome_build="auto"),
                )

            self.assertEqual(table.config_snapshot.genome_build, "hg19")
            self.assertEqual(table.data["POS"].tolist(), [100])

    @unittest.skipUnless(_HAS_PYARROW, "pyarrow is required for sumstats parquet coverage")
    def test_load_sumstats_reads_footer_parquet_with_exact_one_glob(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            sumstats_file = tmpdir / "trait.v1.parquet"
            self._write_footer_parquet(
                sumstats_file,
                pd.DataFrame({"SNP": ["rs1"], "Z": [1.5], "N": [1000.0], "A1": ["A"], "A2": ["G"], "FRQ": [0.2]}),
                snp_identifier="rsid_allele_aware",
                trait_name="trait",
            )

            self.assertTrue(hasattr(ldsc, "load_sumstats"))
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                table = ldsc.load_sumstats(str(tmpdir / "trait*.parquet"), trait_name="trait")

            self.assertEqual(table.source_path, str(sumstats_file))
            self.assertEqual(table.trait_name, "trait")
            self.assertTrue(table.has_alleles)
            self.assertEqual(table.data.loc[0, "SNP"], "rs1")
            self.assertEqual(table.config_snapshot, GlobalConfig(snp_identifier="rsid_allele_aware"))
            self.assertFalse(any("cannot recover the GlobalConfig" in str(item.message) for item in caught))

    def test_load_sumstats_reads_legacy_sumstats_gz_without_metadata(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            sumstats_file = tmpdir / "trait.sumstats.gz"
            with gzip.open(sumstats_file, "wt", encoding="utf-8") as handle:
                handle.write("SNP\tA1\tA2\tZ\tN\nrs1\tA\tG\t1.5\t1000\n")

            table = ldsc.load_sumstats(sumstats_file, trait_name="trait")

            self.assertIsNone(table.config_snapshot)
            self.assertEqual(table.trait_name, "trait")
            self.assertTrue(table.has_alleles)
            self.assertEqual(table.data.loc[0, "SNP"], "rs1")

    def test_load_sumstats_reads_uncompressed_legacy_sumstats(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            sumstats_file = tmpdir / "trait.sumstats"
            sumstats_file.write_text("SNP\tZ\tN\nrs1\t1.5\t1000\n", encoding="utf-8")

            table = ldsc.load_sumstats(sumstats_file, trait_name="trait")

            self.assertIsNone(table.config_snapshot)
            self.assertEqual(table.data.columns.tolist(), ["SNP", "N", "Z"])
            self.assertEqual(table.data.loc[0, "SNP"], "rs1")

    @unittest.skipUnless(_HAS_PYARROW, "pyarrow is required for sumstats parquet coverage")
    def test_load_sumstats_rejects_footerless_parquet(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            sumstats_file = tmpdir / "plain.parquet"
            pd.DataFrame({"SNP": ["rs1"], "A1": ["A"], "A2": ["G"], "Z": [1.5], "N": [1000.0]}).to_parquet(
                sumstats_file, index=False
            )

            with self.assertRaisesRegex(ldsc.LDSCInputError, "missing required LDSC3 identity footer"):
                ldsc.load_sumstats(sumstats_file)

    @unittest.skipUnless(_HAS_PYARROW, "pyarrow is required for sumstats parquet coverage")
    def test_load_sumstats_recovers_trait_name_from_footer(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            sumstats_file = tmpdir / "trait.parquet"
            self._write_footer_parquet(
                sumstats_file, pd.DataFrame({"SNP": ["rs1"], "Z": [1.5], "N": [1000.0]}), trait_name="MDD"
            )

            table = ldsc.load_sumstats(sumstats_file)

            self.assertEqual(table.trait_name, "MDD")

    @unittest.skipUnless(_HAS_PYARROW, "pyarrow is required for sumstats parquet coverage")
    def test_load_sumstats_explicit_trait_name_overrides_footer(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            sumstats_file = tmpdir / "trait.parquet"
            self._write_footer_parquet(
                sumstats_file, pd.DataFrame({"SNP": ["rs1"], "Z": [1.5], "N": [1000.0]}), trait_name="MDD"
            )

            table = ldsc.load_sumstats(sumstats_file, trait_name=" SCZ ")

            self.assertEqual(table.trait_name, "SCZ")

    @unittest.skipUnless(_HAS_PYARROW, "pyarrow is required for sumstats parquet coverage")
    def test_load_sumstats_reads_chr_pos_footer_parquet_with_exact_one_glob(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            sumstats_file = tmpdir / "trait.parquet"
            self._write_footer_parquet(
                sumstats_file,
                pd.DataFrame({"SNP": ["rs1"], "CHR": ["1"], "POS": [100], "Z": [1.5], "N": [1000.0]}),
                snp_identifier="chr_pos",
                genome_build="hg38",
                trait_name=None,
            )

            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                table = ldsc.load_sumstats(str(tmpdir / "trait*.parquet"), trait_name="trait")

            self.assertEqual(table.source_path, str(sumstats_file))
            self.assertEqual(table.config_snapshot, GlobalConfig(snp_identifier="chr_pos", genome_build="hg38"))
            self.assertFalse(any("cannot recover the GlobalConfig" in str(item.message) for item in caught))

    @unittest.skipUnless(_HAS_PYARROW, "pyarrow is required for sumstats parquet coverage")
    def test_load_sumstats_rejects_allele_aware_artifact_without_alleles(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            sumstats_file = tmpdir / "trait.parquet"
            self._write_footer_parquet(
                sumstats_file,
                pd.DataFrame({"SNP": ["rs1"], "Z": [1.0], "N": [100.0]}),
                snp_identifier="rsid_allele_aware",
                trait_name="trait",
            )

            with self.assertRaisesRegex(
                ldsc.LDSCInputError,
                "Cannot load curated sumstats at .*snp_identifier='rsid_allele_aware' requires A1/A2 columns.*"
                "Re-run munge-sumstats from the raw GWAS file.*"
                "docs/troubleshooting.md#munge-sumstats-curated-artifact-is-malformed-or-outdated",
            ):
                ldsc.load_sumstats(sumstats_file)

    @unittest.skipUnless(_HAS_PYARROW, "pyarrow is required for sumstats parquet coverage")
    def test_load_sumstats_rejects_chr_pos_artifact_with_missing_base_identity(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            sumstats_file = tmpdir / "trait.parquet"
            self._write_footer_parquet(
                sumstats_file,
                pd.DataFrame(
                    {
                        "SNP": ["rs_missing_chr", "rs_missing_pos"],
                        "CHR": [pd.NA, "1"],
                        "POS": [100, pd.NA],
                        "Z": [1.0, 1.0],
                        "N": [100.0, 100.0],
                    }
                ),
                snp_identifier="chr_pos",
                genome_build="hg38",
                trait_name="trait",
            )

            with self.assertRaisesRegex(
                ldsc.LDSCInputError,
                "missing or invalid SNP identity rows.*"
                "docs/troubleshooting.md#munge-sumstats-curated-artifact-is-malformed-or-outdated",
            ):
                ldsc.load_sumstats(sumstats_file)

    def test_load_sumstats_rejects_unknown_suffix(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "trait.csv"
            path.write_text("SNP,Z,N\nrs1,1.5,1000\n", encoding="utf-8")

            with self.assertRaisesRegex(ldsc.LDSCInputError, "unsupported file suffix.*run munge-sumstats first"):
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
            raw = MungeConfig(no_snp_restriction=True, raw_sumstats_file=str(raw_path), trait_name="trait")
            config = MungeConfig(no_snp_restriction=True, output_dir=str(tmpdir / "munged"))
            munger = SumstatsMunger()
            table = munger.run(raw, config, GlobalConfig(snp_identifier="rsid"))
            self.assertEqual(table.trait_name, "trait")
            self.assertEqual(table.config_snapshot, GlobalConfig(snp_identifier="rsid"))
            self.assertTrue((tmpdir / "munged" / "trait.parquet").exists())
            self.assertFalse((tmpdir / "munged" / "trait.sumstats.gz").exists())
            output = pd.read_parquet(tmpdir / "munged" / "trait.parquet")
            self.assertEqual(output.columns.tolist(), ["SNP", "CHR", "POS", "A1", "A2", "Z", "N"])
            self.assertTrue(output["CHR"].isna().all())
            self.assertTrue(output["POS"].isna().all())
            self.assertFalse((tmpdir / "munged" / "metadata.json").exists())
            self.assertFalse((tmpdir / "munged" / "sumstats.metadata.json").exists())
            import pyarrow.parquet as pq

            footer = pq.read_schema(str(tmpdir / "munged" / "trait.parquet")).metadata
            self.assertEqual(footer[b"ldsc:artifact_type"], b"sumstats")
            self.assertEqual(footer[b"ldsc:snp_identifier"], b"rsid")
            self.assertEqual(footer[b"ldsc:genome_build"], b"")
            self.assertEqual(footer[b"ldsc:trait_name"], b"trait")
            self.assertNotIn(b"ldsc:schema_version", footer)
            summary = munger.build_run_summary()
            self.assertEqual(summary.n_retained_rows, 2)
            self.assertIn("sumstats_parquet", summary.output_paths)
            self.assertEqual(summary.output_paths["sumstats_parquet"], str(tmpdir / "munged" / "trait.parquet"))
            self.assertNotIn("sumstats_gz", summary.output_paths)
            self.assertNotIn("log", summary.output_paths)
            self.assertNotIn("metadata_json", summary.output_paths)
            self.assertEqual(summary.inferred_columns["detected_format"], "plain")

    @unittest.skipUnless(_HAS_PYARROW, "pyarrow is required for default parquet output")
    def test_run_base_rsid_without_allele_columns_defaults_to_allele_free_table(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            raw_path = tmpdir / "raw.tsv"
            raw_path.write_text(
                "SNP P BETA N\n"
                "rs1 0.05 0.1 1000\n",
                encoding="utf-8",
            )
            output_dir = tmpdir / "munged"

            table = SumstatsMunger().run(
                MungeConfig(no_snp_restriction=True, raw_sumstats_file=raw_path, trait_name="trait"),
                MungeConfig(no_snp_restriction=True, output_dir=output_dir),
                GlobalConfig(snp_identifier="rsid"),
            )

            self.assertFalse(table.has_alleles)
            self.assertEqual(table.data["SNP"].tolist(), ["rs1"])
            self.assertNotIn("A1", table.data.columns)
            output = pd.read_parquet(output_dir / "trait.parquet")
            self.assertEqual(output["SNP"].tolist(), ["rs1"])

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
                MungeConfig(no_snp_restriction=True, raw_sumstats_file=str(raw_path), trait_name="trait"),
                MungeConfig(no_snp_restriction=True, output_dir=str(tmpdir / "munged"), output_format="tsv.gz"),
                GlobalConfig(snp_identifier="rsid"),
            )

            self.assertEqual(len(table.data), 1)
            self.assertTrue((tmpdir / "munged" / "trait.sumstats.gz").exists())
            self.assertFalse((tmpdir / "munged" / "trait.parquet").exists())
            with gzip.open(tmpdir / "munged" / "trait.sumstats.gz", "rt", encoding="utf-8") as handle:
                output = pd.read_csv(handle, sep="\t")
            self.assertEqual(output.columns.tolist(), ["SNP", "CHR", "POS", "A1", "A2", "Z", "N"])
            # tsv.gz carries no embedded metadata and no metadata.json sidecar is written.
            self.assertFalse((tmpdir / "munged" / "metadata.json").exists())

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
                MungeConfig(no_snp_restriction=True, raw_sumstats_file=str(raw_path), trait_name="trait"),
                MungeConfig(no_snp_restriction=True, output_dir=str(tmpdir / "munged"), output_format="both"),
                GlobalConfig(snp_identifier="rsid"),
            )

            import pyarrow.parquet as pq

            self.assertTrue((tmpdir / "munged" / "trait.parquet").exists())
            self.assertTrue((tmpdir / "munged" / "trait.sumstats.gz").exists())
            self.assertFalse((tmpdir / "munged" / "metadata.json").exists())
            footer = pq.read_schema(str(tmpdir / "munged" / "trait.parquet")).metadata
            self.assertEqual(footer[b"ldsc:artifact_type"], b"sumstats")
            self.assertEqual(footer[b"ldsc:trait_name"], b"trait")

    @unittest.skipUnless(_HAS_PYARROW, "pyarrow is required for sumstats parquet coverage")
    def test_run_overwrite_removes_unselected_sumstats_sibling(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            raw_path = tmpdir / "raw.tsv"
            raw_path.write_text("SNP A1 A2 P OR N\nrs1 A G 0.05 1.0 1000\n", encoding="utf-8")
            output_dir = tmpdir / "munged"

            SumstatsMunger().run(
                MungeConfig(no_snp_restriction=True, raw_sumstats_file=raw_path, trait_name="trait"),
                MungeConfig(no_snp_restriction=True, output_dir=output_dir, output_format="both"),
                GlobalConfig(snp_identifier="rsid"),
            )
            self.assertTrue((output_dir / "trait.sumstats.gz").exists())

            SumstatsMunger().run(
                MungeConfig(no_snp_restriction=True, raw_sumstats_file=raw_path, trait_name="trait"),
                MungeConfig(no_snp_restriction=True, output_dir=output_dir, output_format="parquet", overwrite=True),
                GlobalConfig(snp_identifier="rsid"),
            )

            self.assertTrue((output_dir / "trait.parquet").exists())
            self.assertFalse((output_dir / "trait.sumstats.gz").exists())
            self.assertFalse((output_dir / "metadata.json").exists())

    def test_run_refuses_unselected_owned_sumstats_sibling_without_overwrite(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            raw_path = tmpdir / "raw.tsv"
            raw_path.write_text("SNP A1 A2 P OR N\nrs1 A G 0.05 1.0 1000\n", encoding="utf-8")
            output_dir = tmpdir / "munged"
            output_dir.mkdir()
            stale = output_dir / "trait.sumstats.gz"
            stale.write_text("stale\n", encoding="utf-8")

            with mock.patch.object(kernel_munge, "munge_sumstats", side_effect=AssertionError("kernel should not run")):
                with self.assertRaisesRegex(FileExistsError, "overwrite"):
                    SumstatsMunger().run(
                        MungeConfig(no_snp_restriction=True, raw_sumstats_file=raw_path, trait_name="trait"),
                        MungeConfig(no_snp_restriction=True, output_dir=output_dir, output_format="parquet"),
                        GlobalConfig(snp_identifier="rsid"),
                    )

            self.assertEqual(stale.read_text(encoding="utf-8"), "stale\n")
            self.assertFalse((output_dir / "diagnostics" / "sumstats.log").exists())

    def test_run_ignores_legacy_root_diagnostics_without_overwrite(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            raw_path = tmpdir / "raw.tsv"
            self._write_raw_sumstats(raw_path)
            output_dir = tmpdir / "munged"
            root_drop_dir = output_dir / "dropped_snps"
            root_drop_dir.mkdir(parents=True)
            (output_dir / "sumstats.log").write_text("legacy log\n", encoding="utf-8")
            (root_drop_dir / "dropped.tsv.gz").write_text("legacy drops\n", encoding="utf-8")

            SumstatsMunger().run(
                MungeConfig(no_snp_restriction=True, raw_sumstats_file=raw_path),
                MungeConfig(no_snp_restriction=True, output_dir=output_dir, output_format="tsv.gz"),
                GlobalConfig(snp_identifier="rsid"),
            )

            self.assertEqual((output_dir / "sumstats.log").read_text(encoding="utf-8"), "legacy log\n")
            self.assertEqual((root_drop_dir / "dropped.tsv.gz").read_text(encoding="utf-8"), "legacy drops\n")
            self.assertFalse((output_dir / "metadata.json").exists())
            self.assertTrue((output_dir / "diagnostics" / "sumstats.log").exists())
            self.assertTrue((output_dir / "diagnostics" / "dropped_snps" / "dropped.tsv.gz").exists())

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

            self.assertEqual(path, str(tmpdir / "out" / "trait.parquet"))
            self.assertTrue((tmpdir / "out" / "trait.parquet").exists())
            self.assertTrue((tmpdir / "out" / "trait.sumstats.gz").exists())

    def test_write_output_refuses_unselected_owned_sumstats_sibling_without_overwrite(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            output_dir = tmpdir / "out"
            output_dir.mkdir()
            stale = output_dir / "trait.sumstats.gz"
            stale.write_text("stale\n", encoding="utf-8")
            table = sumstats_workflow.SumstatsTable(
                data=pd.DataFrame({"SNP": ["rs1"], "Z": [1.5], "N": [1000.0]}),
                has_alleles=False,
                source_path="source.tsv",
                trait_name="trait",
                config_snapshot=GlobalConfig(snp_identifier="rsid"),
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
            stale = output_dir / "trait.sumstats.gz"
            stale.write_text("stale\n", encoding="utf-8")
            table = sumstats_workflow.SumstatsTable(
                data=pd.DataFrame({"SNP": ["rs1"], "Z": [1.5], "N": [1000.0]}),
                has_alleles=False,
                source_path="source.tsv",
                trait_name="trait",
                config_snapshot=GlobalConfig(snp_identifier="rsid"),
            )

            SumstatsMunger().write_output(table, output_dir, output_format="parquet", overwrite=True)

            self.assertTrue((output_dir / "trait.parquet").exists())
            self.assertFalse((output_dir / "metadata.json").exists())
            self.assertFalse((output_dir / "sumstats.metadata.json").exists())
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

            with mock.patch.object(kernel_munge, "munge_sumstats", return_value=kernel_munge.MungeResult(
                returned, len(returned), {}, "input_columns",
                {"snp_identifier": "chr_pos", "genome_build": "hg38", "liftover": {"applied": False}}, pd.DataFrame(), pd.DataFrame(),
            )):
                SumstatsMunger().run(
                    MungeConfig(no_snp_restriction=True, raw_sumstats_file=str(raw_path), trait_name="trait"),
                    MungeConfig(no_snp_restriction=True,
                        output_dir=str(tmpdir / "munged"),
                        source_genome_build="hg38",
                        output_genome_build="hg38",
                    ),
                    GlobalConfig(snp_identifier="chr_pos", genome_build="hg38"),
                )

            parquet_path = tmpdir / "munged" / "trait.parquet"
            output = pd.read_parquet(parquet_path)
            self.assertEqual(output["SNP"].tolist(), ["rs3", "rs1", "rs2", "rs4"])
            self.assertAlmostEqual(output.loc[0, "Z"], -0.333333333)
            self.assertAlmostEqual(output.loc[1, "N"], 2000.987654)
            self.assertTrue(pd.isna(output.loc[3, "CHR"]))
            self.assertTrue(pd.isna(output.loc[3, "POS"]))
            parquet_file = pq.ParquetFile(parquet_path)
            self.assertEqual(parquet_file.num_row_groups, 3)
            self.assertFalse((tmpdir / "munged" / "metadata.json").exists())
            footer = parquet_file.schema_arrow.metadata
            self.assertEqual(footer[b"ldsc:artifact_type"], b"sumstats")

    @unittest.skipUnless(_HAS_PYARROW, "pyarrow is required for default parquet output")
    def test_run_writes_thin_sidecar_and_logs_coordinate_liftover_provenance(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            raw_path = tmpdir / "raw.tsv"
            raw_path.write_text("CHR POS SNP A1 A2 P BETA N\n1 100 rs1 A G 0.05 0.1 1000\n", encoding="utf-8")


            SumstatsMunger().run(
                MungeConfig(no_snp_restriction=True, raw_sumstats_file=raw_path, trait_name="trait"),
                MungeConfig(no_snp_restriction=True, output_dir=tmpdir / "munged", source_genome_build="hg38", output_genome_build="hg38"),
                GlobalConfig(snp_identifier="chr_pos", genome_build="hg38"),
            )

            import pyarrow.parquet as pq

            self.assertFalse((tmpdir / "munged" / "metadata.json").exists())
            footer = pq.read_schema(str(tmpdir / "munged" / "trait.parquet")).metadata
            self.assertEqual(footer[b"ldsc:artifact_type"], b"sumstats")
            self.assertEqual(footer[b"ldsc:snp_identifier"], b"chr_pos")
            self.assertEqual(footer[b"ldsc:genome_build"], b"hg38")
            self.assertEqual(footer[b"ldsc:trait_name"], b"trait")
            self.assertNotIn(b"ldsc:schema_version", footer)
            log_text = (tmpdir / "munged" / "diagnostics" / "sumstats.log").read_text(encoding="utf-8")
            self.assertIn("Summary-statistics coordinate provenance", log_text)
            self.assertIn("coordinate_basis=1-based", log_text)
            self.assertIn("Summary-statistics liftover report", log_text)
            self.assertNotIn('"coordinate_basis"', log_text)

    def test_run_rejects_liftover_request_in_rsid_mode_before_kernel_call(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            raw_path = tmpdir / "raw.tsv"
            raw_path.write_text("SNP A1 A2 P BETA N\nrs1 A G 0.05 0.1 1000\n", encoding="utf-8")

            with mock.patch.object(kernel_munge, "munge_sumstats", side_effect=AssertionError("kernel should not run")):
                with self.assertRaisesRegex(ldsc.LDSCUsageError, "liftover.*chr_pos"):
                    SumstatsMunger().run(
                        MungeConfig(no_snp_restriction=True, raw_sumstats_file=raw_path, trait_name="trait"),
                        MungeConfig(no_snp_restriction=True,
                            output_dir=tmpdir / "munged",
                            liftover_chain_file="unused.chain",
                        ),
                        GlobalConfig(snp_identifier="rsid"),
                    )

    def test_run_requires_liftover_method_when_target_differs_from_resolved_source(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            raw_path = tmpdir / "raw.tsv"
            raw_path.write_text("CHR POS SNP A1 A2 P BETA N\n1 100 rs1 A G 0.05 0.1 1000\n", encoding="utf-8")

            with self.assertRaisesRegex(ldsc.LDSCUsageError, "liftover method"):
                SumstatsMunger().run(
                    MungeConfig(no_snp_restriction=True, raw_sumstats_file=raw_path, trait_name="trait"),
                    MungeConfig(no_snp_restriction=True, output_dir=tmpdir / "munged", source_genome_build="hg19", output_genome_build="hg38"),
                    GlobalConfig(snp_identifier="chr_pos", genome_build="hg19"),
                )

    def test_run_writes_header_only_dropped_snps_sidecar_when_no_liftover_drops(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            raw_path = tmpdir / "raw.tsv"
            self._write_raw_sumstats(raw_path)


            SumstatsMunger().run(
                MungeConfig(no_snp_restriction=True, raw_sumstats_file=raw_path),
                MungeConfig(no_snp_restriction=True, output_dir=tmpdir / "munged", output_format="tsv.gz"),
                GlobalConfig(snp_identifier="rsid"),
            )

            sidecar = tmpdir / "munged" / "diagnostics" / "dropped_snps" / "dropped.tsv.gz"
            self.assertTrue(sidecar.exists())
            dropped = self._read_dropped_snps_sidecar(sidecar)
            self.assertEqual(len(dropped), 0)
            self.assertEqual(dropped.columns.tolist(), IDENTITY_DROP_COLUMNS)
            self.assertEqual({column: str(dtype) for column, dtype in dropped.dtypes.items()}, self.DROPPED_SNP_DTYPES)

    def test_run_summary_includes_dropped_snps_sidecar_unconditionally(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            raw_path = tmpdir / "raw.tsv"
            self._write_raw_sumstats(raw_path)
            munger = SumstatsMunger()

            table = munger.run(
                MungeConfig(no_snp_restriction=True, raw_sumstats_file=raw_path),
                MungeConfig(no_snp_restriction=True, output_dir=tmpdir / "munged", output_format="tsv.gz"),
                GlobalConfig(snp_identifier="rsid"),
            )

            summary = munger.build_run_summary()
            self.assertEqual(
                summary.output_paths["dropped_snps_tsv_gz"],
                str(tmpdir / "munged" / "diagnostics" / "dropped_snps" / "dropped.tsv.gz"),
            )

    def test_dropped_snps_sidecar_preflight_blocks_existing_file(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            raw_path = tmpdir / "raw.tsv"
            self._write_raw_sumstats(raw_path)
            sidecar = tmpdir / "munged" / "diagnostics" / "dropped_snps" / "dropped.tsv.gz"
            sidecar.parent.mkdir(parents=True)
            sidecar.write_text("stale\n", encoding="utf-8")

            with mock.patch.object(kernel_munge, "munge_sumstats", side_effect=AssertionError("kernel should not run")):
                with self.assertRaises(FileExistsError):
                    SumstatsMunger().run(
                        MungeConfig(no_snp_restriction=True, raw_sumstats_file=raw_path),
                        MungeConfig(no_snp_restriction=True, output_dir=tmpdir / "munged", output_format="tsv.gz"),
                        GlobalConfig(snp_identifier="rsid"),
                    )

    def test_dropped_snps_sidecar_is_overwritten_when_overwrite_true(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            raw_path = tmpdir / "raw.tsv"
            self._write_raw_sumstats(raw_path)
            sidecar = tmpdir / "munged" / "diagnostics" / "dropped_snps" / "dropped.tsv.gz"
            sidecar.parent.mkdir(parents=True)
            sidecar.write_text("stale\n", encoding="utf-8")

            SumstatsMunger().run(
                MungeConfig(no_snp_restriction=True, raw_sumstats_file=raw_path),
                MungeConfig(no_snp_restriction=True, output_dir=tmpdir / "munged", output_format="tsv.gz", overwrite=True),
                GlobalConfig(snp_identifier="rsid"),
            )

            dropped = self._read_dropped_snps_sidecar(sidecar)
            self.assertEqual(len(dropped), 0)
            with gzip.open(sidecar, "rt", encoding="utf-8") as handle:
                self.assertEqual(handle.readline().strip(), "\t".join(IDENTITY_DROP_COLUMNS))

    def test_sumstats_log_points_to_dropped_snps_sidecar_at_info(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            raw_path = tmpdir / "raw.tsv"
            self._write_raw_sumstats(raw_path)
            drop_frame = pd.DataFrame(
                {
                    "CHR": pd.Series(["1"], dtype="string"),
                    "SNP": pd.Series(["rs_drop"], dtype="string"),
                    "source_pos": pd.Series([100], dtype="Int64"),
                    "target_pos": pd.Series([pd.NA], dtype="Int64"),
                    "reason": pd.Series(["unmapped_liftover"], dtype="string"),
                }
            )

            def fake_munge(request):
                return kernel_munge.MungeResult(
                    self._fake_munged_frame(), 2, {"liftover": 1}, "input_columns",
                    {**request.coordinate_metadata, "liftover": {
                        "applied": True, "method": "hm3_curated", "source_build": "hg19", "target_build": "hg38",
                        "n_input": 2, "n_lifted": 1, "n_dropped": 1, "n_unmapped": 1,
                        "n_missing_chr_pos_dropped": 0, "n_cross_chrom": 0,
                        "n_duplicate_source_dropped": 0, "n_duplicate_target_dropped": 0,
                    }}, drop_frame, pd.DataFrame(),
                )

            with mock.patch.object(kernel_munge, "munge_sumstats", side_effect=fake_munge):
                SumstatsMunger().run(
                    MungeConfig(no_snp_restriction=True, raw_sumstats_file=raw_path),
                    MungeConfig(no_snp_restriction=True, output_dir=tmpdir / "munged", output_format="tsv.gz"),
                    GlobalConfig(snp_identifier="rsid"),
                )

            sidecar = tmpdir / "munged" / "diagnostics" / "dropped_snps" / "dropped.tsv.gz"
            log_text = (tmpdir / "munged" / "diagnostics" / "sumstats.log").read_text(encoding="utf-8")
            dropped = self._read_dropped_snps_sidecar(sidecar)
            self.assertEqual(dropped.columns.tolist(), IDENTITY_DROP_COLUMNS)
            self.assertEqual(dropped["stage"].tolist(), ["liftover"])
            self.assertIn(str(sidecar), log_text)
            self.assertIn("unmapped_liftover=1", log_text)
            self.assertNotIn("rs_drop", log_text)

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
                no_snp_restriction=True,
                raw_sumstats_file=str(raw_path),
                output_dir=str(tmpdir / "munged"),
                trait_name="trait",
            )

            table = SumstatsMunger().run(config, global_config=GlobalConfig(snp_identifier="rsid"))

            self.assertEqual(table.trait_name, "trait")
            self.assertEqual(table.source_path, str(raw_path))
            self.assertTrue((tmpdir / "munged" / "trait.parquet").exists())

    def test_top_level_wrapper_calls_main(self):
        from ldsc import sumstats_munger as munge_sumstats

        with mock.patch.object(munge_sumstats, "main", return_value=11) as patched:
            rc = munge_sumstats.__dict__["main"](["--help"])
        patched.assert_called_once()
        self.assertEqual(rc, 11)

    def test_main_resolves_source_build_before_calling_kernel(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            raw_path = tmpdir / "raw.tsv"
            raw_path.write_text("CHR POS SNP A1 A2 P BETA N\n1 99 rs1 A G .05 .05 1000\n")
            inference = mock.Mock(genome_build="hg38", coordinate_basis="0-based", inspected_snp_count=1,
                                  match_counts={}, match_fractions={}, summary_message="inferred hg38 0-based")
            with mock.patch.object(munge_input, "resolve_chr_pos_table", return_value=(pd.DataFrame(), inference)), \
                 mock.patch.object(kernel_munge, "munge_sumstats", wraps=kernel_munge.munge_sumstats) as kernel:
                sumstats_workflow.main(["--no-snp-restriction",
                    "--raw-sumstats-file", str(raw_path), "--output-dir", str(tmpdir / "out"),
                    "--snp-identifier", "chr_pos", "--source-genome-build", "auto", "--output-genome-build", "hg38",
                ])
            self.assertEqual(kernel.call_args.args[0].genome_build, "hg38")
            table = ldsc.load_sumstats(tmpdir / "out" / "sumstats.parquet")
            self.assertEqual(table.data.POS.tolist(), [100])

    def test_main_accepts_sumstats_snps_file(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            raw_path = tmpdir / "raw.tsv"
            raw_path.write_text("SNP P BETA N\nrs1 .05 .05 1000\nrs2 .05 -.05 1000\n")
            keep_path = tmpdir / "keep.tsv"
            keep_path.write_text("SNP\nrs2\n")
            sumstats_workflow.main([
                "--raw-sumstats-file", str(raw_path), "--output-dir", str(tmpdir / "out"),
                "--snp-identifier", "rsid", "--sumstats-snps-file", str(keep_path),
            ])
            table = ldsc.load_sumstats(tmpdir / "out" / "sumstats.parquet")
            self.assertEqual(table.data.SNP.tolist(), ["rs2"])

    def test_main_requires_output_genome_build_for_default_coordinate_mode(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            raw_path = tmpdir / "raw.tsv"
            raw_path.write_text("CHR POS SNP A1 A2 P N\n1 100 rs1 A G 0.05 1000\n", encoding="utf-8")

            with self.assertRaisesRegex(ldsc.LDSCUsageError, "Pass --output-genome-build"):
                sumstats_workflow.main(
                    ["--no-snp-restriction",
                        "--raw-sumstats-file",
                        str(raw_path),
                        "--output-dir",
                        str(tmpdir / "out"),
                    ]
                )

    def test_run_creates_output_dir(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            raw_path = tmpdir / "raw.tsv"
            raw_path.write_text(
                "SNP A1 A2 P OR N\n"
                "rs1 A G 0.05 1.0 1000\n",
                encoding="utf-8",
            )
            raw = MungeConfig(no_snp_restriction=True, raw_sumstats_file=str(raw_path), trait_name="trait")
            output_dir = tmpdir / "nested" / "dir" / "munged"
            config = MungeConfig(no_snp_restriction=True, output_dir=str(output_dir))

            table = SumstatsMunger().run(raw, config, GlobalConfig(snp_identifier="rsid"))

            self.assertEqual(len(table.data), 1)
            self.assertTrue((output_dir / "trait.parquet").exists())
            self.assertTrue((output_dir / "diagnostics" / "sumstats.log").exists())
            log_text = (output_dir / "diagnostics" / "sumstats.log").read_text(encoding="utf-8")
            self.assertIn("LDSC munge-sumstats Started", log_text)
            self.assertIn("Munging summary statistics", log_text)

    def test_run_refuses_existing_trait_outputs_before_kernel_call(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            raw_path = tmpdir / "raw.tsv"
            raw_path.write_text("SNP A1 A2 P OR N\nrs1 A G 0.05 1.0 1000\n", encoding="utf-8")
            output_dir = tmpdir / "munged"
            output_dir.mkdir()
            existing = output_dir / "trait.parquet"
            existing.write_text("existing\n", encoding="utf-8")

            with mock.patch.object(kernel_munge, "munge_sumstats", side_effect=AssertionError("kernel should not run")):
                with self.assertRaisesRegex(FileExistsError, "overwrite"):
                    SumstatsMunger().run(
                        MungeConfig(no_snp_restriction=True, raw_sumstats_file=raw_path, trait_name="trait"),
                        MungeConfig(no_snp_restriction=True, output_dir=output_dir),
                        GlobalConfig(snp_identifier="rsid"),
                    )

            self.assertEqual(existing.read_text(encoding="utf-8"), "existing\n")
            self.assertFalse((output_dir / "diagnostics" / "sumstats.log").exists())

    def test_run_allows_existing_trait_outputs_with_overwrite(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            raw_path = tmpdir / "raw.tsv"
            raw_path.write_text("SNP A1 A2 P OR N\nrs1 A G 0.05 1.0 1000\n", encoding="utf-8")
            output_dir = tmpdir / "munged"
            output_dir.mkdir()
            (output_dir / "trait.parquet").write_text("existing\n", encoding="utf-8")
            returned = pd.DataFrame({"SNP": ["rs1"], "N": [1000.0], "Z": [1.5], "A1": ["A"], "A2": ["G"]})

            table = SumstatsMunger().run(
                MungeConfig(no_snp_restriction=True, raw_sumstats_file=raw_path, trait_name="trait"),
                MungeConfig(no_snp_restriction=True, output_dir=output_dir, overwrite=True),
                GlobalConfig(snp_identifier="rsid"),
            )

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
            raw = MungeConfig(no_snp_restriction=True, raw_sumstats_file=raw_path, trait_name="trait")
            output_dir = tmpdir / "nested" / "mdd2025"

            table = SumstatsMunger().run(raw, MungeConfig(no_snp_restriction=True, output_dir=output_dir), GlobalConfig(snp_identifier="rsid"))

            self.assertEqual(len(table.data), 1)
            self.assertEqual(table.source_path, str(raw_path))
            self.assertTrue((output_dir / "trait.parquet").exists())
            self.assertTrue((output_dir / "diagnostics" / "sumstats.log").exists())

    def test_run_resolves_glob_pattern_for_single_sumstats_input(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            raw_path = tmpdir / "trait.raw.tsv"
            raw_path.write_text(
                "SNP A1 A2 P OR N\n"
                "rs1 A G 0.05 1.0 1000\n",
                encoding="utf-8",
            )
            raw = MungeConfig(no_snp_restriction=True, raw_sumstats_file=str(tmpdir / "trait.*.tsv"), trait_name="trait")
            output_dir = tmpdir / "out" / "munged"

            table = SumstatsMunger().run(raw, MungeConfig(no_snp_restriction=True, output_dir=output_dir), GlobalConfig(snp_identifier="rsid"))

            self.assertEqual(table.source_path, str(raw_path))
            self.assertTrue((output_dir / "trait.parquet").exists())

    def test_workflow_missing_sample_size_error_names_fix_options(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            raw_path = tmpdir / "raw.tsv"
            raw_path.write_text("SNP A1 A2 P OR\nrs1 A G 0.05 1.1\n", encoding="utf-8")

            with self.assertRaisesRegex(
                ldsc.LDSCInputError,
                "munge-sumstats could not determine sample size \\(N\\).*"
                "--N.*--N-cas.*--N-con.*N column",
            ):
                SumstatsMunger().run(MungeConfig(no_snp_restriction=True, raw_sumstats_file=raw_path, output_dir=tmpdir / "out"), global_config=GlobalConfig(snp_identifier="rsid"))

    def test_sumstats_table_missing_required_columns_uses_input_error(self):
        table = sumstats_workflow.SumstatsTable(
            data=pd.DataFrame({"SNP": ["rs1"], "N": [1000.0]}),
            has_alleles=False,
            source_path="trait.sumstats.gz",
            trait_name="trait",
        )

        with self.assertRaisesRegex(
            ldsc.LDSCInputError,
            "munge-sumstats could not map required column.*Z.*"
            "docs/troubleshooting.md#munge-sumstats-could-not-map-a-required-column",
        ):
            table.validate()

    def test_workflow_error_path_does_not_write_progress_to_stdout(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            raw_path = tmpdir / "raw.tsv"
            raw_path.write_text("SNP A1 A2 P OR\nrs1 A G 0.05 1.1\n", encoding="utf-8")
            stdout = io.StringIO()

            with contextlib.redirect_stdout(stdout):
                with self.assertRaisesRegex(ldsc.LDSCInputError, "munge-sumstats could not determine sample size"):
                    SumstatsMunger().run(MungeConfig(no_snp_restriction=True, raw_sumstats_file=raw_path, output_dir=tmpdir / "out"), global_config=GlobalConfig(snp_identifier="rsid"))

            self.assertEqual(stdout.getvalue(), "")

    def test_internal_sample_size_guard_uses_internal_error(self):
        args = kernel_munge.MungeQC()

        with self.assertRaisesRegex(
            ldsc.LDSCInternalError,
            "munge-sumstats could not derive a sample size \\(N\\).*unreachable.*--debug",
        ):
            kernel_munge.process_n(pd.DataFrame({"P": [0.05]}), args)

    def test_default_allele_aware_mode_rejects_missing_alleles_and_suggests_base_mode(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            raw_path = tmpdir / "raw.tsv"
            raw_path.write_text("CHR POS SNP P BETA N\n1 100 rs1 0.05 0.1 1000\n", encoding="utf-8")

            with self.assertRaisesRegex(
                ldsc.LDSCInputError,
                "snp_identifier='chr_pos_allele_aware'.*--snp-identifier chr_pos",
            ):
                SumstatsMunger().run(MungeConfig(no_snp_restriction=True, raw_sumstats_file=raw_path, output_dir=tmpdir / "out", source_genome_build="hg38", output_genome_build="hg38"), global_config=GlobalConfig(snp_identifier="chr_pos_allele_aware", genome_build="hg38"))

    def test_parser_rejects_removed_no_alleles_flag(self):
        with self.assertRaises(SystemExit):
            sumstats_workflow.build_parser().parse_args(
                [
                    "--raw-sumstats-file",
                    "raw.tsv",
                    "--output-dir",
                    "sumstats",
                    "--no-alleles",
                ]
            )

    def test_workflow_infers_ref_alt_allele_columns(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            raw_path = Path(tmpdir) / "raw.tsv"
            raw_path.write_text("SNP REF ALT P OR N\nrs1 A G 0.05 1.0 1000\n")
            table = SumstatsMunger().run(
                MungeConfig(no_snp_restriction=True, raw_sumstats_file=raw_path, output_dir=Path(tmpdir) / "out"),
                global_config=GlobalConfig(snp_identifier="rsid_allele_aware"),
            )
            self.assertEqual(table.data.A1.tolist(), ["A"])
            self.assertEqual(table.data.A2.tolist(), ["G"])

    def test_rsid_allele_aware_requires_alleles(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            raw_path = tmpdir / "raw.tsv"
            raw_path.write_text("SNP P BETA N\nrs1 0.05 0.1 1000\n", encoding="utf-8")
            with self.assertRaisesRegex(
                ldsc.LDSCInputError,
                "snp_identifier='rsid_allele_aware'.*--snp-identifier rsid",
            ):
                SumstatsMunger().run(MungeConfig(no_snp_restriction=True, raw_sumstats_file=raw_path, output_dir=tmpdir / "out"), global_config=GlobalConfig(snp_identifier="rsid_allele_aware"))

    def test_base_modes_run_without_alleles_without_extra_flag(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            raw_path = tmpdir / "raw.tsv"
            raw_path.write_text("CHR POS SNP P BETA N\n1 100 rs1 0.05 0.1 1000\n")
            for mode in ("rsid", "chr_pos"):
                build = {"source_genome_build": "hg38", "output_genome_build": "hg38"} if mode == "chr_pos" else {}
                table = SumstatsMunger().run(
                    MungeConfig(no_snp_restriction=True, raw_sumstats_file=raw_path, output_dir=tmpdir / mode, **build),
                    global_config=GlobalConfig(snp_identifier=mode, genome_build="hg38" if mode == "chr_pos" else None),
                )
                self.assertEqual(table.data.SNP.tolist(), ["rs1"])

    def test_base_modes_keep_singletons_with_bad_allele_columns(self):
        rows = (
            "CHR POS SNP A1 A2 P BETA N\n"
            "1 100 missing . C 0.05 0.1 1000\n"
            "1 200 ambiguous A T 0.05 0.1 1000\n"
            "1 300 multibase AT C 0.05 0.1 1000\n"
            "1 400 identical G G 0.05 0.1 1000\n"
        )
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            for mode in ("rsid", "chr_pos"):
                raw_path = tmpdir / f"{mode}.tsv"
                raw_path.write_text(rows, encoding="utf-8")
                build = {"source_genome_build": "hg38", "output_genome_build": "hg38"} if mode == "chr_pos" else {}
                munged = SumstatsMunger().run(
                    MungeConfig(no_snp_restriction=True, raw_sumstats_file=raw_path, output_dir=tmpdir / f"{mode}_out", **build),
                    global_config=GlobalConfig(snp_identifier=mode, genome_build="hg38" if mode == "chr_pos" else None),
                ).data

                self.assertEqual(munged["SNP"].tolist(), ["missing", "ambiguous", "multibase", "identical"])

    def test_base_mode_duplicate_drops_use_duplicate_identity(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            raw_path = tmpdir / "raw.tsv"
            raw_path.write_text(
                "SNP A1 A2 P BETA N\n"
                "rs1 A C 0.05 0.1 1000\n"
                "rs1 A G 0.05 0.1 1000\n"
                "rs2 A T 0.05 0.1 1000\n",
                encoding="utf-8",
            )
            output_dir = tmpdir / "munged"

            table = SumstatsMunger().run(
                MungeConfig(no_snp_restriction=True, raw_sumstats_file=raw_path),
                MungeConfig(no_snp_restriction=True, output_dir=output_dir),
                GlobalConfig(snp_identifier="rsid"),
            )

            self.assertEqual(table.data["SNP"].tolist(), ["rs2"])
            dropped = self._read_dropped_snps_sidecar(output_dir / "diagnostics" / "dropped_snps" / "dropped.tsv.gz")
            self.assertEqual(dropped["reason"].tolist(), ["duplicate_identity", "duplicate_identity"])
            self.assertEqual(dropped["stage"].tolist(), ["post_liftover_identity_cleanup", "post_liftover_identity_cleanup"])

    def test_allele_aware_sumstats_snps_file_does_not_raise_on_bad_raw_alleles(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            raw_path = tmpdir / "raw.tsv"
            raw_path.write_text(
                "SNP A1 A2 P BETA N\n"
                "rs1 A C 0.05 0.1 1000\n"
                "rs2 AT C 0.05 0.1 1000\n"
                "rs3 A T 0.05 0.1 1000\n",
                encoding="utf-8",
            )
            keep_path = tmpdir / "keep.tsv"
            keep_path.write_text(
                "SNP\tA1\tA2\n"
                "rs1\tA\tC\n"
                "rs2\tA\tC\n"
                "rs3\tA\tC\n",
                encoding="utf-8",
            )
            output_dir = tmpdir / "munged"

            table = SumstatsMunger().run(
                MungeConfig(no_snp_restriction=True, raw_sumstats_file=raw_path),
                MungeConfig(output_dir=output_dir, sumstats_snps_file=keep_path),
                GlobalConfig(snp_identifier="rsid_allele_aware"),
            )

            self.assertEqual(table.data["SNP"].tolist(), ["rs1"])
            dropped = self._read_dropped_snps_sidecar(output_dir / "diagnostics" / "dropped_snps" / "dropped.tsv.gz")
            self.assertEqual(
                dropped["reason"].tolist(),
                ["invalid_allele", "strand_ambiguous_allele"],
            )
            self.assertEqual(
                dropped["stage"].tolist(),
                ["post_liftover_identity_cleanup", "post_liftover_identity_cleanup"],
            )

    def test_rsid_allele_aware_sumstats_snps_file_filters_bad_raw_alleles_outside_keep_list(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            raw_path = tmpdir / "raw.tsv"
            raw_path.write_text(
                "SNP A1 A2 P BETA N\n"
                "rs_keep A C 0.05 0.1 1000\n"
                "rs_bad_out AT C 0.05 0.1 1000\n",
                encoding="utf-8",
            )
            keep_path = tmpdir / "keep.tsv"
            keep_path.write_text("SNP\tA1\tA2\nrs_keep\tA\tC\n", encoding="utf-8")
            output_dir = tmpdir / "munged"

            table = SumstatsMunger().run(
                MungeConfig(no_snp_restriction=True, raw_sumstats_file=raw_path),
                MungeConfig(output_dir=output_dir, sumstats_snps_file=keep_path),
                GlobalConfig(snp_identifier="rsid_allele_aware"),
            )

            self.assertEqual(table.data["SNP"].tolist(), ["rs_keep"])
            dropped = self._read_dropped_snps_sidecar(output_dir / "diagnostics" / "dropped_snps" / "dropped.tsv.gz")
            self.assertNotIn("rs_bad_out", dropped["SNP"].tolist())

    def test_chr_pos_allele_aware_sumstats_snps_file_filters_bad_raw_alleles_outside_keep_list(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            raw_path = tmpdir / "raw.tsv"
            raw_path.write_text(
                "CHR POS SNP A1 A2 P BETA N\n"
                "1 100 rs_keep A C 0.05 0.1 1000\n"
                "1 200 rs_bad_out AT C 0.05 0.1 1000\n",
                encoding="utf-8",
            )
            keep_path = tmpdir / "keep.tsv"
            keep_path.write_text("CHR\tPOS\tA1\tA2\n1\t100\tA\tC\n", encoding="utf-8")
            output_dir = tmpdir / "munged"

            table = SumstatsMunger().run(
                MungeConfig(no_snp_restriction=True, raw_sumstats_file=raw_path),
                MungeConfig(
                    output_dir=output_dir,
                    sumstats_snps_file=keep_path,
                    source_genome_build="hg38",
                    output_genome_build="hg38",
                ),
                GlobalConfig(snp_identifier="chr_pos_allele_aware", genome_build="hg38"),
            )

            self.assertEqual(table.data["SNP"].tolist(), ["rs_keep"])
            dropped = self._read_dropped_snps_sidecar(output_dir / "diagnostics" / "dropped_snps" / "dropped.tsv.gz")
            self.assertNotIn("rs_bad_out", dropped["SNP"].tolist())

    def test_allele_aware_missing_raw_alleles_are_reported_in_identity_sidecar(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            raw_path = tmpdir / "raw.tsv"
            raw_path.write_text(
                "SNP A1 A2 P BETA N\n"
                "rs1 A C 0.05 0.1 1000\n"
                "rs_missing . C 0.05 0.1 1000\n",
                encoding="utf-8",
            )
            output_dir = tmpdir / "munged"

            table = SumstatsMunger().run(
                MungeConfig(no_snp_restriction=True, raw_sumstats_file=raw_path),
                MungeConfig(no_snp_restriction=True, output_dir=output_dir),
                GlobalConfig(snp_identifier="rsid_allele_aware"),
            )

            self.assertEqual(table.data["SNP"].tolist(), ["rs1"])
            dropped = self._read_dropped_snps_sidecar(output_dir / "diagnostics" / "dropped_snps" / "dropped.tsv.gz")
            self.assertEqual(dropped["SNP"].tolist(), ["rs_missing"])
            self.assertEqual(dropped["reason"].tolist(), ["missing_allele"])
            self.assertEqual(dropped["stage"].tolist(), ["post_liftover_identity_cleanup"])

    def test_allele_aware_all_missing_raw_allele_chunk_reaches_identity_sidecar(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            raw_path = tmpdir / "raw.tsv"
            raw_path.write_text(
                "SNP A1 A2 P BETA N\n"
                "rs_missing_1 . . 0.05 0.1 1000\n"
                "rs_missing_2 . . 0.05 0.1 1000\n",
                encoding="utf-8",
            )
            output_dir = tmpdir / "munged"

            table = SumstatsMunger().run(
                MungeConfig(no_snp_restriction=True, raw_sumstats_file=raw_path),
                MungeConfig(no_snp_restriction=True, output_dir=output_dir),
                GlobalConfig(snp_identifier="rsid_allele_aware"),
            )

            self.assertEqual(table.data["SNP"].tolist(), [])
            dropped = self._read_dropped_snps_sidecar(output_dir / "diagnostics" / "dropped_snps" / "dropped.tsv.gz")
            self.assertEqual(dropped["SNP"].tolist(), ["rs_missing_1", "rs_missing_2"])
            self.assertEqual(dropped["reason"].tolist(), ["missing_allele", "missing_allele"])

    def test_allele_aware_sumstats_snps_file_all_bad_raw_alleles_reaches_identity_sidecar(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            raw_path = tmpdir / "raw.tsv"
            raw_path.write_text(
                "SNP A1 A2 P BETA N\n"
                "rs_bad AT C 0.05 0.1 1000\n"
                "rs_ambiguous A T 0.05 0.1 1000\n",
                encoding="utf-8",
            )
            keep_path = tmpdir / "keep.tsv"
            keep_path.write_text(
                "SNP\tA1\tA2\n"
                "rs_bad\tA\tC\n"
                "rs_ambiguous\tA\tC\n",
                encoding="utf-8",
            )
            output_dir = tmpdir / "munged"

            table = SumstatsMunger().run(
                MungeConfig(no_snp_restriction=True, raw_sumstats_file=raw_path),
                MungeConfig(output_dir=output_dir, sumstats_snps_file=keep_path),
                GlobalConfig(snp_identifier="rsid_allele_aware"),
            )

            self.assertEqual(table.data["SNP"].tolist(), [])
            dropped = self._read_dropped_snps_sidecar(output_dir / "diagnostics" / "dropped_snps" / "dropped.tsv.gz")
            self.assertEqual(dropped["SNP"].tolist(), ["rs_bad", "rs_ambiguous"])
            self.assertEqual(dropped["reason"].tolist(), ["invalid_allele", "strand_ambiguous_allele"])

    @unittest.skipUnless(_HAS_PYARROW, "pyarrow is required for sumstats parquet coverage")
    def test_load_sumstats_rejects_duplicate_effective_keys(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            sumstats_file = tmpdir / "trait.parquet"
            self._write_footer_parquet(
                sumstats_file,
                pd.DataFrame({"SNP": ["rs1", "rs1"], "Z": [1.0, 2.0], "N": [100.0, 100.0]}),
                snp_identifier="rsid",
                trait_name="trait",
            )

            with self.assertRaisesRegex(
                ldsc.LDSCInputError,
                "duplicate or invalid SNP identity rows.*duplicate_identity.*"
                "docs/troubleshooting.md#munge-sumstats-curated-artifact-is-malformed-or-outdated",
            ):
                ldsc.load_sumstats(sumstats_file)

    def test_workflow_missing_sample_size_error_explains_neff_is_not_n(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            raw_path = tmpdir / "raw.tsv"
            raw_path.write_text("SNP A1 A2 P OR NEFF\nrs1 A G 0.05 1.1 1000\n", encoding="utf-8")

            with self.assertRaisesRegex(ldsc.LDSCInputError, "NEFF is not treated as N.*--N-col NEFF"):
                SumstatsMunger().run(MungeConfig(no_snp_restriction=True, raw_sumstats_file=raw_path, output_dir=tmpdir / "out"), global_config=GlobalConfig(snp_identifier="rsid"))

    def test_workflow_missing_signed_stat_error_suggests_likely_effect_column(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            raw_path = tmpdir / "raw.tsv"
            raw_path.write_text("SNP A1 A2 P EFFECT_SIZE N\nrs1 A G 0.05 0.2 1000\n", encoding="utf-8")

            with self.assertRaisesRegex(ldsc.LDSCInputError, "--signed-sumstats EFFECT_SIZE,0"):
                SumstatsMunger().run(MungeConfig(no_snp_restriction=True, raw_sumstats_file=raw_path, output_dir=tmpdir / "out"), global_config=GlobalConfig(snp_identifier="rsid"))

    def test_parser_rejects_removed_merge_alleles(self):
        with self.assertRaises(SystemExit):
            sumstats_workflow.build_parser().parse_args(["--raw-sumstats-file", "raw.tsv", "--output-dir", "sumstats", "--merge-alleles", "hm3.tsv.gz"])

    def test_run_restricts_sumstats_snps_file_by_rsid_without_reordering(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            raw_path = tmpdir / "raw.tsv"
            raw_path.write_text(
                "SNP A1 A2 P BETA N\n"
                "rs1 A G 0.05 0.1 1000\n"
                "rs2 C T 0.10 -0.1 1000\n"
                "rs3 G A 0.20 0.0 1000\n",
                encoding="utf-8",
            )
            keep_path = tmpdir / "keep.tsv"
            keep_path.write_text("SNP\nrs3\nrs1\n", encoding="utf-8")

            table = SumstatsMunger().run(
                MungeConfig(no_snp_restriction=True, raw_sumstats_file=raw_path, trait_name="trait"),
                MungeConfig(output_dir=tmpdir / "munged", sumstats_snps_file=keep_path),
                GlobalConfig(snp_identifier="rsid"),
            )

            self.assertEqual(table.data["SNP"].tolist(), ["rs1", "rs3"])
            log_text = (tmpdir / "munged" / "diagnostics" / "sumstats.log").read_text(encoding="utf-8")
            self.assertIn("Applying SNP keep-list restriction", log_text)
            self.assertIn("snp_identifier=rsid", log_text)
            self.assertIn("read 2 keep-list identifiers", log_text)

    def test_run_restricts_default_hm3_by_rsid(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            raw_path = tmpdir / "raw.tsv"
            raw_path.write_text(
                "SNP A1 A2 P BETA N\n"
                "rs1 A G 0.05 0.1 1000\n"
                "rs2 C T 0.10 -0.1 1000\n",
                encoding="utf-8",
            )
            hm3_path = tmpdir / "hm3.tsv"
            hm3_path.write_text("CHR\thg19_POS\thg38_POS\tSNP\n1\t10\t20\trs2\n", encoding="utf-8")

            with mock.patch("ldsc.sumstats_munger.packaged_hm3_curated_map_path", return_value=str(hm3_path)):
                table = SumstatsMunger().run(
                    MungeConfig(raw_sumstats_file=raw_path, trait_name="trait"),
                    MungeConfig(output_dir=tmpdir / "munged"),
                    GlobalConfig(snp_identifier="rsid"),
                )

        self.assertEqual(table.data["SNP"].tolist(), ["rs2"])

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
                MungeConfig(no_snp_restriction=True, raw_sumstats_file=raw_path, trait_name="trait"),
                MungeConfig(
                    output_dir=tmpdir / "munged",
                    sumstats_snps_file=keep_path,
                    source_genome_build="hg38",
                    output_genome_build="hg38",
                ),
                GlobalConfig(snp_identifier="chr_pos", genome_build="hg38"),
            )

            self.assertEqual(table.data["SNP"].tolist(), ["rs2"])
            self.assertEqual(table.data["POS"].tolist(), [200])

    def test_run_restricts_default_hm3_by_source_build_chr_pos(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            raw_path = tmpdir / "raw.tsv"
            raw_path.write_text(
                "CHR POS SNP A1 A2 P BETA N\n"
                "1 100 rs1 A G 0.05 0.1 1000\n"
                "1 200 rs2 C T 0.10 -0.1 1000\n",
                encoding="utf-8",
            )
            hm3_path = tmpdir / "hm3.tsv"
            hm3_path.write_text("CHR\thg19_POS\thg38_POS\tSNP\n1\t200\t999\trs2\n", encoding="utf-8")

            with mock.patch("ldsc.sumstats_munger.packaged_hm3_curated_map_path", return_value=str(hm3_path)):
                table = SumstatsMunger().run(
                    MungeConfig(raw_sumstats_file=raw_path, trait_name="trait"),
                    MungeConfig(
                        output_dir=tmpdir / "munged",
                        source_genome_build="hg19",
                        output_genome_build="hg19",
                    ),
                    GlobalConfig(snp_identifier="chr_pos", genome_build="hg19"),
                )

        self.assertEqual(table.data["SNP"].tolist(), ["rs2"])
        self.assertEqual(table.data["POS"].tolist(), [200])

    def test_run_restricts_default_hm3_by_rsid_alleles(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            raw_path = tmpdir / "raw.tsv"
            raw_path.write_text(
                "SNP A1 A2 P BETA N\n"
                "rs1 A G 0.05 0.1 1000\n"
                "rs1 A C 0.10 -0.1 1000\n",
                encoding="utf-8",
            )
            hm3_path = tmpdir / "hm3.tsv"
            hm3_path.write_text("CHR\thg19_POS\thg38_POS\tSNP\tA1\tA2\n1\t100\t200\trs1\tA\tC\n", encoding="utf-8")

            with mock.patch("ldsc.sumstats_munger.packaged_hm3_curated_map_path", return_value=str(hm3_path)):
                table = SumstatsMunger().run(
                    MungeConfig(raw_sumstats_file=raw_path, trait_name="trait"),
                    MungeConfig(output_dir=tmpdir / "munged"),
                    GlobalConfig(snp_identifier="rsid_allele_aware"),
                )

        self.assertEqual(table.data["SNP"].tolist(), ["rs1"])
        self.assertEqual(table.data["A1"].tolist(), ["A"])
        self.assertEqual(table.data["A2"].tolist(), ["C"])

    def test_run_restricts_default_hm3_by_source_build_chr_pos_alleles(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            raw_path = tmpdir / "raw.tsv"
            raw_path.write_text(
                "CHR POS SNP A1 A2 P BETA N\n"
                "1 100 rs1 A G 0.05 0.1 1000\n"
                "1 100 rs2 A C 0.10 -0.1 1000\n",
                encoding="utf-8",
            )
            hm3_path = tmpdir / "hm3.tsv"
            hm3_path.write_text("CHR\thg19_POS\thg38_POS\tSNP\tA1\tA2\n1\t100\t200\trs2\tA\tC\n", encoding="utf-8")

            with mock.patch("ldsc.sumstats_munger.packaged_hm3_curated_map_path", return_value=str(hm3_path)):
                table = SumstatsMunger().run(
                    MungeConfig(raw_sumstats_file=raw_path, trait_name="trait"),
                    MungeConfig(
                        output_dir=tmpdir / "munged",
                        source_genome_build="hg19",
                        output_genome_build="hg19",
                    ),
                    GlobalConfig(snp_identifier="chr_pos_allele_aware", genome_build="hg19"),
                )

        self.assertEqual(table.data["SNP"].tolist(), ["rs2"])
        self.assertEqual(table.data["POS"].tolist(), [100])
        self.assertEqual(table.data["A1"].tolist(), ["A"])
        self.assertEqual(table.data["A2"].tolist(), ["C"])

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
                MungeConfig(no_snp_restriction=True, raw_sumstats_file=raw_path, trait_name="trait"),
                MungeConfig(
                    output_dir=tmpdir / "munged",
                    sumstats_snps_file=keep_path,
                    source_genome_build="hg38",
                    output_genome_build="hg38",
                ),
                GlobalConfig(snp_identifier="chr_pos", genome_build="hg38"),
            )

            self.assertEqual(table.data["SNP"].tolist(), ["rs1"])
            log_text = (tmpdir / "munged" / "diagnostics" / "sumstats.log").read_text(encoding="utf-8")
            self.assertIn("Dropped 2 SNPs with invalid or missing CHR/POS", log_text)
            self.assertIn("Removed 0 SNPs with missing values.", log_text)
            self.assertIn("Removed 0 SNPs with out-of-bounds p-values.", log_text)
            self.assertIn("missing_chr", log_text)

    def test_run_drops_invalid_chr_pos_rows_before_writing_outputs(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            raw_path = tmpdir / "raw.tsv"
            raw_path.write_text(
                "CHR POS SNP A1 A2 P BETA N\n"
                "1 100 rs1 A G 0.05 0.1 1000\n"
                "chrUn 200 bad_chr C T 0.10 -0.1 1000\n"
                "1 abc bad_pos G A 0.20 0.2 1000\n",
                encoding="utf-8",
            )

            table = SumstatsMunger().run(
                MungeConfig(no_snp_restriction=True, raw_sumstats_file=raw_path, trait_name="trait"),
                MungeConfig(no_snp_restriction=True, output_dir=tmpdir / "munged", source_genome_build="hg38", output_genome_build="hg38"),
                GlobalConfig(snp_identifier="chr_pos", genome_build="hg38"),
            )

            self.assertEqual(table.data["SNP"].tolist(), ["rs1"])
            self.assertEqual(table.data["POS"].tolist(), [100])
            self.assertEqual(table.provenance["coordinate_provenance"]["n_dropped_invalid_chr_pos"], 2)
            log_text = (tmpdir / "munged" / "diagnostics" / "sumstats.log").read_text(encoding="utf-8")
            self.assertIn("Dropped 2 SNPs with invalid or missing CHR/POS", log_text)
            self.assertIn("Removed 0 SNPs with missing values.", log_text)
            self.assertIn("Removed 0 SNPs with out-of-bounds p-values.", log_text)
            self.assertIn("bad_chr", log_text)

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
                ldsc.LDSCInputError,
                "munge-sumstats no SNPs remain after SNP keep-list restriction.*"
                "snp_identifier=rsid.*keep-list identifiers=1.*"
                "docs/troubleshooting.md#munge-sumstats-no-snps-remain-after-filtering",
            ):
                SumstatsMunger().run(
                    MungeConfig(no_snp_restriction=True, raw_sumstats_file=raw_path, trait_name="trait"),
                    MungeConfig(output_dir=tmpdir / "munged", sumstats_snps_file=keep_path),
                    GlobalConfig(snp_identifier="rsid"),
                )

    def test_run_rejects_empty_sumstats_snps_file_before_chunk_parsing(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            raw_path = tmpdir / "raw.tsv"
            raw_path.write_text("SNP A1 A2 P BETA N\nrs1 A G 0.05 0.1 1000\n", encoding="utf-8")
            keep_path = tmpdir / "keep.tsv"
            keep_path.write_text("SNP\n", encoding="utf-8")

            with mock.patch.object(kernel_munge, "parse_dat", side_effect=AssertionError("raw chunks parsed")):
                with self.assertRaisesRegex(ldsc.LDSCInputError, "keep-list identifiers=0"):
                    SumstatsMunger().run(
                        MungeConfig(no_snp_restriction=True, raw_sumstats_file=raw_path, trait_name="trait"),
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
                MungeConfig(no_snp_restriction=True, raw_sumstats_file=raw_path, trait_name="trait"),
                MungeConfig(
                    output_dir=tmpdir / "munged",
                    sumstats_snps_file=keep_path,
                    source_genome_build="hg19",
                    output_genome_build="hg19",
                ),
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

            with mock.patch.object(munge_input, "resolve_chr_pos_table", side_effect=fake_resolve_chr_pos_table):
                table = SumstatsMunger().run(
                    MungeConfig(no_snp_restriction=True, raw_sumstats_file=raw_path, trait_name="trait"),
                    MungeConfig(
                        output_dir=tmpdir / "munged",
                        sumstats_snps_file=keep_path,
                        source_genome_build="auto",
                        output_genome_build="hg19",
                    ),
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
                MungeConfig(no_snp_restriction=True, raw_sumstats_file=raw_path, trait_name="trait"),
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
                MungeConfig(no_snp_restriction=True, raw_sumstats_file=raw_path, trait_name="trait"),
                MungeConfig(no_snp_restriction=True, output_dir=tmpdir / "munged"),
                GlobalConfig(snp_identifier="rsid"),
            )

            self.assertEqual(table.data["SNP"].tolist(), ["rs1"])
            self.assertEqual(table.data["N"].tolist(), [1000.0])
            self.assertEqual(table.data["CHR"].tolist(), [1])
            self.assertEqual(table.data["POS"].tolist(), [123])

    def test_explicit_n_col_suppresses_inferred_case_control_columns(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            raw_path = tmpdir / "raw.tsv"
            raw_path.write_text(
                "SNP A1 A2 P BETA NEFF NCAS NCON\n"
                "rs1 A G 0.05 0.1 777 400 600\n"
                "rs2 C T 0.10 -0.1 888 400 600\n",
                encoding="utf-8",
            )
            (tmpdir / "munged").mkdir()

            with self.assertWarnsRegex(UserWarning, "--N-col NEFF.*ignored automatically inferred.*NCAS.*NCON"):
                table = SumstatsMunger().run(
                    MungeConfig(no_snp_restriction=True, raw_sumstats_file=raw_path, column_hints={"N_col": "NEFF"}),
                    MungeConfig(no_snp_restriction=True, output_dir=tmpdir / "munged"),
                    GlobalConfig(snp_identifier="rsid"),
                )

            self.assertEqual(table.data["N"].tolist(), [777.0, 888.0])

    def test_explicit_case_control_columns_suppress_inferred_direct_n(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            raw_path = tmpdir / "raw.tsv"
            raw_path.write_text(
                "SNP A1 A2 P BETA N CASES CONTROLS\n"
                "rs1 A G 0.05 0.1 777 400 600\n",
                encoding="utf-8",
            )
            (tmpdir / "munged").mkdir()

            with self.assertWarnsRegex(
                UserWarning,
                "--N-cas-col CASES.*--N-con-col CONTROLS.*ignored automatically inferred.*N",
            ):
                table = SumstatsMunger().run(
                    MungeConfig(no_snp_restriction=True,
                        raw_sumstats_file=raw_path,
                        column_hints={"N_cas_col": "CASES", "N_con_col": "CONTROLS"},
                    ),
                    MungeConfig(no_snp_restriction=True, output_dir=tmpdir / "munged"),
                    GlobalConfig(snp_identifier="rsid"),
                )

            self.assertEqual(table.data["N"].tolist(), [1000.0])

    def test_explicit_direct_and_case_control_column_strategies_conflict(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            raw_path = tmpdir / "raw.tsv"
            raw_path.write_text(
                "SNP A1 A2 P BETA NEFF NCAS NCON\n"
                "rs1 A G 0.05 0.1 777 400 600\n",
                encoding="utf-8",
            )

            with self.assertRaisesRegex(
                ldsc.LDSCUsageError,
                "--N-col.*cannot be combined with --N-cas-col and --N-con-col.*Choose one sample-size strategy",
            ):
                SumstatsMunger().run(
                    MungeConfig(no_snp_restriction=True,
                        raw_sumstats_file=raw_path,
                        column_hints={"N_col": "NEFF", "N_cas_col": "NCAS", "N_con_col": "NCON"},
                    ),
                    MungeConfig(no_snp_restriction=True, output_dir=tmpdir / "munged"),
                    GlobalConfig(snp_identifier="rsid"),
                )

    def test_explicit_case_control_column_hints_must_be_a_pair(self):
        incomplete_hints = (
            ({"N_cas_col": "NCAS"}, "--N-con-col"),
            ({"N_con_col": "NCON"}, "--N-cas-col"),
        )
        for hint, expected_missing in incomplete_hints:
            with self.subTest(hint=hint), tempfile.TemporaryDirectory() as tmpdir:
                tmpdir = Path(tmpdir)
                raw_path = tmpdir / "raw.tsv"
                raw_path.write_text(
                    "SNP A1 A2 P BETA N NCAS NCON\n"
                    "rs1 A G 0.05 0.1 777 400 600\n",
                    encoding="utf-8",
                )

                with self.assertRaisesRegex(
                    ldsc.LDSCUsageError,
                    f"--N-cas-col and --N-con-col must be provided together.*{expected_missing}",
                ):
                    SumstatsMunger().run(
                        MungeConfig(no_snp_restriction=True, raw_sumstats_file=raw_path, column_hints=hint),
                        MungeConfig(no_snp_restriction=True, output_dir=tmpdir / "munged"),
                        GlobalConfig(snp_identifier="rsid"),
                    )

    def test_automatic_direct_and_case_control_inference_is_ambiguous(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            raw_path = tmpdir / "raw.tsv"
            raw_path.write_text(
                "SNP A1 A2 P BETA N NCAS NCON\n"
                "rs1 A G 0.05 0.1 777 400 600\n",
                encoding="utf-8",
            )

            with self.assertRaisesRegex(
                ldsc.LDSCInputError,
                "multiple sample-size strategies.*--N-col N.*--N-cas-col NCAS --N-con-col NCON",
            ):
                SumstatsMunger().run(
                    MungeConfig(no_snp_restriction=True, raw_sumstats_file=raw_path),
                    MungeConfig(no_snp_restriction=True, output_dir=tmpdir / "munged"),
                    GlobalConfig(snp_identifier="rsid"),
                )

    def test_unambiguous_direct_and_case_control_sample_size_workflows(self):
        fixtures = (
            ("SNP A1 A2 P BETA N\nrs1 A G 0.05 0.1 777\n", 777.0),
            ("SNP A1 A2 P BETA NCAS NCON\nrs1 A G 0.05 0.1 400 600\n", 1000.0),
        )
        for raw_text, expected_n in fixtures:
            with self.subTest(expected_n=expected_n), tempfile.TemporaryDirectory() as tmpdir:
                tmpdir = Path(tmpdir)
                raw_path = tmpdir / "raw.tsv"
                raw_path.write_text(raw_text, encoding="utf-8")

                table = SumstatsMunger().run(
                    MungeConfig(no_snp_restriction=True, raw_sumstats_file=raw_path),
                    MungeConfig(no_snp_restriction=True, output_dir=tmpdir / "munged"),
                    GlobalConfig(snp_identifier="rsid"),
                )

                self.assertEqual(table.data["N"].tolist(), [expected_n])

    def test_run_auto_detects_old_daner_sample_sizes(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            raw_path = tmpdir / "daner.tsv"
            raw_path.write_text(
                "CHR SNP BP A1 A2 FRQ_A_40 FRQ_U_60 INFO OR SE P\n"
                "1 rs1 100 A G 0.2 0.3 0.99 1.02 0.1 0.05\n",
                encoding="utf-8",
            )

            table = SumstatsMunger().run(
                MungeConfig(no_snp_restriction=True, raw_sumstats_file=raw_path, trait_name="trait"),
                MungeConfig(no_snp_restriction=True, output_dir=tmpdir / "munged", source_genome_build="hg38", output_genome_build="hg38"),
                GlobalConfig(snp_identifier="chr_pos", genome_build="hg38"),
            )

            self.assertEqual(table.data["N"].tolist(), [100.0])

    def test_run_does_not_treat_neff_as_total_n(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            raw_path = tmpdir / "raw.tsv"
            raw_path.write_text(
                "#CHROM\tPOS\tID\tEA\tNEA\tBETA\tPVAL\tNEFF\n"
                "1\t123\trs1\tA\tG\t0.1\t0.05\t1000\n",
                encoding="utf-8",
            )

            with self.assertRaisesRegex(ldsc.LDSCInputError, "NEFF is not treated as N"):
                SumstatsMunger().run(
                    MungeConfig(no_snp_restriction=True, raw_sumstats_file=raw_path, trait_name="trait"),
                    MungeConfig(no_snp_restriction=True, output_dir=tmpdir / "munged", source_genome_build="hg38", output_genome_build="hg38"),
                    GlobalConfig(snp_identifier="chr_pos", genome_build="hg38"),
                )

    def test_run_auto_detects_new_daner_frequency_and_sample_sizes(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            raw_path = tmpdir / "daner_new.tsv"
            raw_path.write_text(
                "SNP A1 A2 FRQ_U_60 P OR Nca Nco\n"
                "rs1 A G 0.3 0.05 1.02 40 60\n",
                encoding="utf-8",
            )

            table = SumstatsMunger().run(
                MungeConfig(no_snp_restriction=True, raw_sumstats_file=raw_path, trait_name="trait"),
                MungeConfig(no_snp_restriction=True, output_dir=tmpdir / "munged"),
                GlobalConfig(snp_identifier="rsid"),
            )

            self.assertEqual(table.data["N"].tolist(), [100])
            self.assertEqual(table.data["FRQ"].tolist(), [0.3])

    def test_run_auto_treats_comma_separated_impinfo_as_info_list(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            raw_path = tmpdir / "raw.tsv"
            raw_path.write_text(
                "#CHROM\tPOS\tID\tEA\tNEA\tBETA\tPVAL\tIMPINFO\tN\n"
                "1\t123\trs1\tA\tG\t0.1\t0.05\t0.95,0.91,NA\t1000\n",
                encoding="utf-8",
            )

            table = SumstatsMunger().run(
                MungeConfig(no_snp_restriction=True, raw_sumstats_file=raw_path, trait_name="trait"),
                MungeConfig(no_snp_restriction=True, output_dir=tmpdir / "munged", source_genome_build="hg38", output_genome_build="hg38"),
                GlobalConfig(snp_identifier="chr_pos", genome_build="hg38"),
            )

            self.assertEqual(table.data["SNP"].tolist(), ["rs1"])

    def test_run_invalid_comma_separated_impinfo_suggests_ignore_or_info_list(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            raw_path = tmpdir / "raw.tsv"
            raw_path.write_text(
                "#CHROM\tPOS\tID\tEA\tNEA\tBETA\tPVAL\tIMPINFO\tN\n"
                "1\t123\trs1\tA\tG\t0.1\t0.05\t0.95,LOW,0.85\t1000\n",
                encoding="utf-8",
            )

            with self.assertRaisesRegex(ldsc.LDSCInputError, "--ignore IMPINFO.*--info-list IMPINFO"):
                SumstatsMunger().run(
                    MungeConfig(no_snp_restriction=True, raw_sumstats_file=raw_path, trait_name="trait"),
                    MungeConfig(no_snp_restriction=True, output_dir=tmpdir / "munged", source_genome_build="hg38", output_genome_build="hg38"),
                    GlobalConfig(snp_identifier="chr_pos", genome_build="hg38"),
                )

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
                MungeConfig(no_snp_restriction=True, raw_sumstats_file=raw_path, trait_name="trait"),
                MungeConfig(no_snp_restriction=True, output_dir=tmpdir / "munged", source_genome_build="hg38", output_genome_build="hg38"),
                GlobalConfig(snp_identifier="chr_pos", genome_build="hg38"),
            )

            self.assertEqual(table.data["SNP"].tolist(), ["rs2073813"])
            self.assertEqual(table.data["CHR"].tolist(), ["1"])
            self.assertEqual(table.data["POS"].tolist(), [753541])
            self.assertEqual(munger.build_run_summary().n_input_rows, 1)

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
                MungeConfig(no_snp_restriction=True, raw_sumstats_file=raw_path, trait_name="trait"),
                MungeConfig(no_snp_restriction=True, output_dir=tmpdir / "munged", source_genome_build="hg38", output_genome_build="hg38"),
                GlobalConfig(snp_identifier="chr_pos", genome_build="hg38"),
            )

            self.assertEqual(table.data["CHR"].tolist(), ["1"])
            self.assertEqual(table.data["POS"].tolist(), [123])
            output = pd.read_parquet(tmpdir / "munged" / "trait.parquet")
            self.assertEqual(output.columns.tolist(), ["SNP", "CHR", "POS", "A1", "A2", "Z", "N"])

    def test_a1_a2_descriptions_explain_signed_statistic_orientation(self):
        self.assertIn("signed statistic", munge_input.describe_cname["A1"])
        self.assertIn("counterpart", munge_input.describe_cname["A2"])
        self.assertNotIn("ref allele", munge_input.describe_cname["A1"].lower())
        self.assertNotIn("non-ref allele", munge_input.describe_cname["A2"].lower())

    def test_signed_statistic_orients_output_z_relative_to_a1(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            raw_path = tmpdir / "raw.tsv"
            raw_path.write_text(
                "SNP A1 A2 P BETA N\n"
                "rs_positive A G 0.05 0.2 1000\n"
                "rs_negative C T 0.05 -0.2 1000\n",
                encoding="utf-8",
            )

            table = SumstatsMunger().run(
                MungeConfig(no_snp_restriction=True, raw_sumstats_file=raw_path, trait_name="trait"),
                MungeConfig(no_snp_restriction=True, output_dir=tmpdir / "munged"),
                GlobalConfig(snp_identifier="rsid"),
            )

            z_by_snp = dict(zip(table.data["SNP"], table.data["Z"]))
            self.assertGreater(z_by_snp["rs_positive"], 0)
            self.assertLess(z_by_snp["rs_negative"], 0)

    def test_infer_only_reports_detected_format_and_suggested_command(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            raw_path = tmpdir / "raw.tsv"
            raw_path.write_text(
                "CHR SNP BP A1 A2 FRQ_A_40 FRQ_U_60 INFO OR SE P\n"
                "1 rs1 100 A G 0.2 0.3 0.99 1.2 0.1 0.05\n",
                encoding="utf-8",
            )
            stdout = io.StringIO()

            with contextlib.redirect_stdout(stdout):
                result = sumstats_workflow.main(
                    ["--no-snp-restriction",
                        "--raw-sumstats-file",
                        str(raw_path),
                        "--output-dir",
                        str(tmpdir / "unused"),
                        "--infer-only",
                        "--source-genome-build",
                        "hg38",
                        "--output-genome-build",
                        "hg38",
                    ]
                )

            output = stdout.getvalue()
            self.assertEqual(result.detected_format, "daner-old")
            self.assertIn("Detected format: daner-old", output)
            self.assertIn("Runnable: yes", output)
            self.assertIn("--input-format daner-old", output)
            self.assertIn("--output-genome-build hg38", output)

    def test_infer_only_reads_gzip_raw_sumstats(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            raw_path = Path(tmpdir) / "raw.tsv.gz"
            with gzip.open(raw_path, "wt", encoding="utf-8") as handle:
                handle.write(
                    "CHR POS SNP A1 A2 P BETA N\n"
                    "1 100 rs1 A G 0.05 0.1 1000\n"
                )

            result = sumstats_workflow.main(
                ["--no-snp-restriction",
                    "--raw-sumstats-file",
                    str(raw_path),
                    "--output-dir",
                    str(Path(tmpdir) / "unused"),
                    "--infer-only",
                    "--source-genome-build",
                    "hg38",
                    "--output-genome-build",
                    "hg38",
                ]
            )

            self.assertEqual(result.detected_format, "plain")
            self.assertTrue(result.runnable)

    def test_infer_only_requires_chain_for_unrestricted_cross_build_run(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            raw = tmpdir / "raw.tsv"
            self._write_raw_sumstats(raw)
            stdout = io.StringIO()
            with contextlib.redirect_stdout(stdout):
                report = sumstats_workflow.main([
                    "--raw-sumstats-file", str(raw), "--output-dir", str(tmpdir / "unused"),
                    "--infer-only", "--source-genome-build", "hg19", "--output-genome-build", "hg38",
                    "--no-snp-restriction",
                ])
            self.assertFalse(report.runnable)
            self.assertIn("liftover_method", report.missing_fields)
            self.assertEqual(report.liftover_method, "missing; chain file required")
            self.assertIn("--no-snp-restriction", report.suggested_args)
            self.assertIn("--liftover-chain-file <hg19ToHg38.over.chain>", stdout.getvalue())
            self.assertFalse((tmpdir / "unused").exists())

    def test_infer_only_reports_source_build_inference_failure_as_non_runnable(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            raw_path = tmpdir / "raw.tsv"
            raw_path.write_text(
                "CHR POS SNP A1 A2 P BETA N\n"
                "1 100 rs1 A G 0.05 0.1 1000\n",
                encoding="utf-8",
            )
            stdout = io.StringIO()

            with mock.patch("ldsc.sumstats_munger.resolve_genome_build", side_effect=ValueError("insufficient evidence")):
                with contextlib.redirect_stdout(stdout):
                    result = sumstats_workflow.main(
                        ["--no-snp-restriction",
                            "--raw-sumstats-file",
                            str(raw_path),
                            "--output-dir",
                            str(tmpdir / "unused"),
                            "--infer-only",
                            "--output-genome-build",
                            "hg38",
                        ]
                    )

            output = stdout.getvalue()
            self.assertFalse(result.runnable)
            self.assertIn("source_genome_build", result.missing_fields)
            self.assertIn("Unable to infer source genome build", output)
            self.assertIn("--source-genome-build hg19", output)
            self.assertIn("--source-genome-build hg38", output)

    def test_infer_only_coordinate_read_stops_after_sufficient_build_evidence(self):
        from ldsc.genome_build_inference import load_packaged_reference_table

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            raw_path = tmpdir / "raw.tsv"
            reference = load_packaged_reference_table().head(300)
            with raw_path.open("w", encoding="utf-8") as handle:
                handle.write("CHR POS SNP A1 A2 P BETA N\n")
                for idx, row in enumerate(reference.itertuples(index=False)):
                    handle.write(f"{row.CHR} {int(row.hg19_POS)} rs{idx} A G 0.05 0.1 1000\n")
                for idx in range(60_000):
                    handle.write(f"1 {900_000_000 + idx} tail{idx} A G 0.05 0.1 1000\n")

            sample = sumstats_workflow._read_infer_only_coordinate_frame(
                str(raw_path),
                MungeConfig(no_snp_restriction=True, raw_sumstats_file=raw_path),
                MungeConfig(no_snp_restriction=True),
                sumstats_workflow.RawSumstatsInference(detected_format="plain"),
            )

            self.assertLess(len(sample), 10_000)
            self.assertGreaterEqual(len(sample), 300)

    def test_infer_only_checks_liftover_chain_path_and_reports_direction(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            raw_path = tmpdir / "raw.tsv"
            raw_path.write_text(
                "CHR POS SNP A1 A2 P BETA N\n"
                "1 100 rs1 A G 0.05 0.1 1000\n",
                encoding="utf-8",
            )
            chain_path = tmpdir / "custom.over.chain"
            chain_path.write_text("", encoding="utf-8")
            stdout = io.StringIO()

            with mock.patch("ldsc.sumstats_munger.resolve_genome_build", return_value="hg19"):
                with contextlib.redirect_stdout(stdout):
                    result = sumstats_workflow.main(
                        ["--no-snp-restriction",
                            "--raw-sumstats-file",
                            str(raw_path),
                            "--output-dir",
                            str(tmpdir / "unused"),
                            "--infer-only",
                            "--output-genome-build",
                            "hg38",
                            "--liftover-chain-file",
                            str(chain_path),
                        ]
                    )

            output = stdout.getvalue()
            self.assertTrue(result.runnable)
            self.assertIn("Liftover required: yes (method: chain file)", output)
            self.assertIn("Expected chain direction: hg19 -> hg38", output)

    def test_infer_only_reports_missing_liftover_chain_path_with_fix(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            raw_path = tmpdir / "raw.tsv"
            raw_path.write_text(
                "CHR POS SNP A1 A2 P BETA N\n"
                "1 100 rs1 A G 0.05 0.1 1000\n",
                encoding="utf-8",
            )
            missing_chain_path = tmpdir / "missing.over.chain"
            stdout = io.StringIO()

            with mock.patch("ldsc.sumstats_munger.resolve_genome_build", return_value="hg19"):
                with contextlib.redirect_stdout(stdout):
                    result = sumstats_workflow.main(
                        ["--no-snp-restriction",
                            "--raw-sumstats-file",
                            str(raw_path),
                            "--output-dir",
                            str(tmpdir / "unused"),
                            "--infer-only",
                            "--output-genome-build",
                            "hg38",
                            "--liftover-chain-file",
                            str(missing_chain_path),
                        ]
                    )

            output = stdout.getvalue()
            self.assertFalse(result.runnable)
            self.assertIn("liftover_chain_file", result.missing_fields)
            self.assertIn("Could not resolve liftover chain file path from token", output)
            self.assertIn("matched 0 files", output)
            self.assertIn("Most likely the path is misspelled", output)
            self.assertIn("pass one existing chain file for the expected direction (hg19ToHg38.over.chain)", output)
            self.assertIn("Expected chain direction: hg19 -> hg38", output)

    def test_infer_only_reports_hm3_quick_liftover_method(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            raw_path = tmpdir / "raw.tsv"
            raw_path.write_text(
                "CHR POS SNP A1 A2 P BETA N\n"
                "1 100 rs1 A G 0.05 0.1 1000\n",
                encoding="utf-8",
            )
            stdout = io.StringIO()

            with mock.patch("ldsc.sumstats_munger.resolve_genome_build", return_value="hg19"):
                with contextlib.redirect_stdout(stdout):
                    result = sumstats_workflow.main(
                        [
                            "--raw-sumstats-file",
                            str(raw_path),
                            "--output-dir",
                            str(tmpdir / "unused"),
                            "--infer-only",
                            "--output-genome-build",
                            "hg38",
                        ]
                    )

            output = stdout.getvalue()
            self.assertTrue(result.runnable)
            self.assertIn("Liftover required: yes (method: hm3 quick)", output)

    def test_infer_only_ignores_liftover_method_when_source_matches_output(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            raw_path = tmpdir / "raw.tsv"
            raw_path.write_text(
                "CHR POS SNP A1 A2 P BETA N\n"
                "1 100 rs1 A G 0.05 0.1 1000\n",
                encoding="utf-8",
            )
            stdout = io.StringIO()

            with contextlib.redirect_stdout(stdout):
                result = sumstats_workflow.main(
                    ["--no-snp-restriction",
                        "--raw-sumstats-file",
                        str(raw_path),
                        "--output-dir",
                        str(tmpdir / "unused"),
                        "--infer-only",
                        "--source-genome-build",
                        "hg38",
                        "--output-genome-build",
                        "hg38",
                        "--liftover-chain-file",
                        str(tmpdir / "missing.over.chain"),
                    ]
                )

            output = stdout.getvalue()
            self.assertTrue(result.runnable)
            self.assertFalse(result.liftover_required)
            self.assertIn("Liftover required: no (method: none)", output)
            self.assertIn("supplied liftover method will be ignored", output)
            self.assertNotIn("liftover_chain_file", result.missing_fields)

    def test_infer_auto_detects_new_daner_case_control_columns(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            raw_path = tmpdir / "daner_new.tsv"
            raw_path.write_text(
                "SNP A1 A2 P OR Nca Nco\n"
                "rs1 A G 0.05 1.02 40 60\n",
                encoding="utf-8",
            )

            result = sumstats_workflow.infer_raw_sumstats(raw_path)

            self.assertEqual(result.detected_format, "daner-new")
            self.assertFalse(result.column_hints)
            self.assertTrue(result.runnable)

    def test_infer_auto_treats_vcf_style_header_as_plain_with_ref_alt_hints(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            raw_path = tmpdir / "plain_vcf_style.tsv"
            raw_path.write_text(
                "##fileformat=VCFv4.2\n"
                "#CHROM POS ID REF ALT BETA PVAL N\n"
                "1 123 rs1 A G 0.1 0.05 1000\n",
                encoding="utf-8",
            )

            result = sumstats_workflow.infer_raw_sumstats(raw_path)

            self.assertEqual(result.detected_format, "plain")
            self.assertIn("--input-format", result.suggested_args)
            self.assertIn("plain", result.suggested_args)
            self.assertEqual(result.column_hints, {"a1": "REF", "a2": "ALT"})
            self.assertTrue(result.runnable)

    def test_infer_only_suggests_likely_signed_sumstats_flag(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            raw_path = tmpdir / "raw.tsv"
            raw_path.write_text(
                "SNP A1 A2 P EFFECT_SIZE N\n"
                "rs1 A G 0.05 0.2 1000\n",
                encoding="utf-8",
            )
            stdout = io.StringIO()

            with contextlib.redirect_stdout(stdout):
                result = sumstats_workflow.main(
                    ["--no-snp-restriction",
                        "--raw-sumstats-file",
                        str(raw_path),
                        "--output-dir",
                        str(tmpdir / "unused"),
                        "--infer-only",
                        "--source-genome-build",
                        "auto",
                        "--output-genome-build",
                        "hg38",
                    ]
                )

            output = stdout.getvalue()
            self.assertEqual(result.signed_sumstats_spec, "EFFECT_SIZE,0")
            self.assertIn("Signed statistic hint: EFFECT_SIZE,0", output)
            self.assertIn("--signed-sumstats EFFECT_SIZE,0", output)

    def test_infer_only_reports_neff_as_missing_n_not_inferred(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            raw_path = tmpdir / "raw.tsv"
            raw_path.write_text(
                "#CHROM\tPOS\tID\tEA\tNEA\tBETA\tPVAL\tNEFF\n"
                "1\t123\trs1\tA\tG\t0.1\t0.05\t1000\n",
                encoding="utf-8",
            )
            stdout = io.StringIO()

            with contextlib.redirect_stdout(stdout):
                result = sumstats_workflow.main(
                    ["--no-snp-restriction",
                        "--raw-sumstats-file",
                        str(raw_path),
                        "--output-dir",
                        str(tmpdir / "unused"),
                        "--infer-only",
                        "--source-genome-build",
                        "auto",
                        "--output-genome-build",
                        "hg38",
                    ]
                )

            output = stdout.getvalue()
            self.assertFalse(result.runnable)
            self.assertIn("Missing fields: N", output)
            self.assertIn("NEFF is not treated as N", output)

    def test_infer_only_base_mode_does_not_report_missing_alleles(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            raw_path = tmpdir / "raw.tsv"
            raw_path.write_text("CHR POS SNP P BETA N\n1 100 rs1 0.05 0.1 1000\n", encoding="utf-8")
            stdout = io.StringIO()

            with contextlib.redirect_stdout(stdout):
                result = sumstats_workflow.main(
                    ["--no-snp-restriction",
                        "--raw-sumstats-file",
                        str(raw_path),
                        "--output-dir",
                        str(tmpdir / "unused"),
                        "--infer-only",
                        "--snp-identifier",
                        "chr_pos",
                        "--source-genome-build",
                        "hg38",
                        "--output-genome-build",
                        "hg38",
                    ]
                )

            self.assertTrue(result.runnable)
            self.assertNotIn("A1/A2", result.missing_fields)
            self.assertNotIn("A1/A2", stdout.getvalue())

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
                MungeConfig(no_snp_restriction=True,
                    raw_sumstats_file=raw_path,
                    trait_name="trait",
                    column_hints={"chr": "chromosome_name", "pos": "base_pair", "snp": "variant"},
                ),
                MungeConfig(no_snp_restriction=True, output_dir=tmpdir / "munged", source_genome_build="hg38", output_genome_build="hg38"),
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

            with self.assertRaisesRegex(ldsc.LDSCInputError, "CHR"):
                SumstatsMunger().run(
                    MungeConfig(no_snp_restriction=True, raw_sumstats_file=raw_path, trait_name="trait"),
                    MungeConfig(no_snp_restriction=True, output_dir=tmpdir / "munged", source_genome_build="hg38", output_genome_build="hg38"),
                    GlobalConfig(snp_identifier="chr_pos", genome_build="hg38"),
                )

    def test_munge_then_load_roundtrips_identity_through_parquet_footer(self):
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
                MungeConfig(no_snp_restriction=True, raw_sumstats_file=raw_path, trait_name="trait"),
                MungeConfig(no_snp_restriction=True, output_dir=output_dir, source_genome_build="hg38", output_genome_build="hg38"),
                GlobalConfig(snp_identifier="chr_pos", genome_build="hg38"),
            )

            import pyarrow.parquet as pq

            self.assertFalse((output_dir / "metadata.json").exists())
            footer = pq.read_schema(str(output_dir / "trait.parquet")).metadata
            self.assertEqual(footer[b"ldsc:artifact_type"], b"sumstats")
            self.assertEqual(footer[b"ldsc:snp_identifier"], b"chr_pos")
            self.assertEqual(footer[b"ldsc:genome_build"], b"hg38")
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                table = ldsc.load_sumstats(output_dir / "trait.parquet", trait_name="trait")
            self.assertEqual(table.config_snapshot, GlobalConfig(snp_identifier="chr_pos", genome_build="hg38"))
            self.assertFalse(any("cannot recover the GlobalConfig" in str(item.message) for item in caught))


@unittest.skipUnless(_HAS_PYARROW, "pyarrow required")
class SumstatsParquetFooterTest(unittest.TestCase):
    def _frame(self):
        return pd.DataFrame(
            {
                "SNP": ["rs1", "rs2"],
                "CHR": ["1", "1"],
                "POS": [10, 20],
                "A1": ["A", "C"],
                "A2": ["G", "T"],
                "Z": [0.1, 0.2],
                "N": [100, 100],
            }
        )

    def test_footer_encodes_identity_and_trait(self):
        import pyarrow.parquet as pq

        snapshot = GlobalConfig(snp_identifier="chr_pos_allele_aware", genome_build="hg19")
        footer = sumstats_workflow._sumstats_footer_metadata(snapshot, trait_name="height")
        with tempfile.TemporaryDirectory() as tmp:
            out = str(Path(tmp) / "sumstats.parquet")
            sumstats_workflow._write_sumstats_outputs(
                self._frame(), output_files={"parquet": out}, output_format="parquet", footer_metadata=footer
            )
            raw = pq.read_schema(out).metadata
        self.assertEqual(raw[b"ldsc:artifact_type"], b"sumstats")
        self.assertEqual(raw[b"ldsc:snp_identifier"], b"chr_pos_allele_aware")
        self.assertEqual(raw[b"ldsc:genome_build"], b"hg19")
        self.assertEqual(raw[b"ldsc:trait_name"], b"height")
        self.assertNotIn(b"ldsc:schema_version", raw)

    def test_footer_encodes_none_genome_build_and_trait_as_empty(self):
        footer = sumstats_workflow._sumstats_footer_metadata(
            GlobalConfig(snp_identifier="rsid_allele_aware"), trait_name=None
        )
        self.assertEqual(footer[b"ldsc:genome_build"], b"")
        self.assertEqual(footer[b"ldsc:trait_name"], b"")
