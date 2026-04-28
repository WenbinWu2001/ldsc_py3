from pathlib import Path
import importlib
import sys
import tempfile
import unittest
import warnings
from unittest import mock

import pandas as pd

SRC = Path(__file__).resolve().parents[1] / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

import ldsc
from ldsc.config import (
    AnnotationBuildConfig,
    GlobalConfig,
    LDScoreConfig,
    MungeConfig,
    ReferencePanelBuildConfig,
    RefPanelConfig,
    RegressionConfig,
    normalize_genome_build,
)
from ldsc._kernel.identifiers import (
    build_chr_pos_snp_id,
    build_snp_id_series,
    clean_header,
    infer_chr_bp_columns,
    infer_snp_column,
    normalize_chromosome,
    normalize_snp_identifier_mode,
    read_global_snp_restriction,
    validate_unique_snp_ids,
)


class GlobalConfigTest(unittest.TestCase):
    def test_defaults(self):
        config = GlobalConfig(snp_identifier="rsid")
        self.assertEqual(config.snp_identifier, "rsid")
        self.assertIsNone(config.genome_build)
        self.assertEqual(config.log_level, "INFO")
        self.assertFalse(config.fail_on_missing_metadata)

    def test_package_global_registry_round_trip(self):
        original = ldsc.reset_global_config()
        self.assertEqual(original, GlobalConfig(snp_identifier="rsid"))

        configured = GlobalConfig(
            snp_identifier="rsid",
            log_level="DEBUG",
        )
        ldsc.set_global_config(configured)

        self.assertEqual(ldsc.get_global_config(), configured)

    def test_reset_global_config_restores_default(self):
        ldsc.set_global_config(GlobalConfig(snp_identifier="rsid", log_level="DEBUG"))

        reset = ldsc.reset_global_config()

        self.assertEqual(reset, GlobalConfig(snp_identifier="rsid"))
        self.assertEqual(ldsc.get_global_config(), GlobalConfig(snp_identifier="rsid"))

    def test_validates_values(self):
        with self.assertRaises(ValueError):
            GlobalConfig(snp_identifier="bad")
        with self.assertRaises(ValueError):
            GlobalConfig(snp_identifier="chr_pos", genome_build="hg18")
        with self.assertRaises(ValueError):
            GlobalConfig(snp_identifier="rsid", log_level="trace")

    def test_normalizes_genome_build_aliases(self):
        self.assertEqual(GlobalConfig(snp_identifier="chr_pos", genome_build="hg37").genome_build, "hg19")
        self.assertEqual(GlobalConfig(snp_identifier="chr_pos", genome_build="GRCh37").genome_build, "hg19")
        self.assertEqual(GlobalConfig(snp_identifier="chr_pos", genome_build="GRCh38").genome_build, "hg38")
        self.assertEqual(normalize_genome_build("hg37"), "hg19")
        self.assertEqual(normalize_genome_build("GRCh37"), "hg19")
        self.assertEqual(normalize_genome_build("GRCh38"), "hg38")
        self.assertEqual(GlobalConfig(snp_identifier="chr_pos", genome_build="auto").genome_build, "auto")
        self.assertEqual(normalize_genome_build("AUTO"), "auto")


class TestGlobalConfigValidation(unittest.TestCase):
    def test_chr_pos_no_genome_build_raises(self):
        with self.assertRaises(ValueError) as ctx:
            GlobalConfig(snp_identifier="chr_pos")
        self.assertIn("genome_build is required", str(ctx.exception))

    def test_chr_pos_hg38_ok(self):
        cfg = GlobalConfig(snp_identifier="chr_pos", genome_build="hg38")
        self.assertEqual(cfg.genome_build, "hg38")

    def test_chr_pos_auto_ok(self):
        cfg = GlobalConfig(snp_identifier="chr_pos", genome_build="auto")
        self.assertEqual(cfg.genome_build, "auto")

    def test_rsid_auto_raises(self):
        with self.assertRaises(ValueError):
            GlobalConfig(snp_identifier="rsid", genome_build="auto")

    def test_rsid_no_genome_build_ok(self):
        cfg = GlobalConfig(snp_identifier="rsid")
        self.assertIsNone(cfg.genome_build)

    def test_rsid_with_concrete_build_warns_and_nulls(self):
        with self.assertWarns(UserWarning) as ctx:
            cfg = GlobalConfig(snp_identifier="rsid", genome_build="hg38")
        self.assertIn("ignored", str(ctx.warning))
        self.assertIsNone(cfg.genome_build)

    def test_singleton_is_rsid(self):
        from ldsc.config import get_global_config
        cfg = get_global_config()
        self.assertEqual(cfg.snp_identifier, "rsid")
        self.assertIsNone(cfg.genome_build)

    def test_reset_returns_rsid(self):
        from ldsc.config import reset_global_config
        cfg = reset_global_config()
        self.assertEqual(cfg.snp_identifier, "rsid")
        self.assertIsNone(cfg.genome_build)


class WorkflowConfigTest(unittest.TestCase):
    def test_annotation_config_normalizes_paths(self):
        config = AnnotationBuildConfig(
            baseline_annot_paths=[Path("base.annot.gz")],
            query_annot_bed_paths=[Path("query.bed")],
            output_dir=Path("outputs") / "annot",
            overwrite=True,
        )
        self.assertEqual(config.baseline_annot_paths, ("base.annot.gz",))
        self.assertEqual(config.query_annot_bed_paths, ("query.bed",))
        self.assertEqual(config.output_dir, "outputs/annot")
        self.assertTrue(config.batch_mode)
        self.assertTrue(config.overwrite)

    def test_annotation_config_accepts_single_string_tokens_for_plural_fields(self):
        config = AnnotationBuildConfig(
            baseline_annot_paths="baseline.@.annot.gz",
            query_annot_paths="query.*.annot.gz",
            query_annot_bed_paths="beds/*.bed",
        )
        self.assertEqual(config.baseline_annot_paths, ("baseline.@.annot.gz",))
        self.assertEqual(config.query_annot_paths, ("query.*.annot.gz",))
        self.assertEqual(config.query_annot_bed_paths, ("beds/*.bed",))

    def test_removed_compatibility_aliases_are_not_exported(self):
        removed_names_by_module = {
            "ldsc": ("AnnotationSourceSpec", "RawSumstatsSpec", "RefPanelSpec", "OutputSpec"),
            "ldsc.annotation_builder": ("AnnotationSourceSpec",),
            "ldsc._kernel.annotation": ("AnnotationSourceSpec",),
            "ldsc._kernel.ref_panel": ("RefPanelSpec",),
            "ldsc.sumstats_munger": ("RawSumstatsSpec",),
            "ldsc.outputs": ("OutputSpec",),
        }
        for module_name, names in removed_names_by_module.items():
            module = importlib.import_module(module_name)
            exports = getattr(module, "__all__", ())
            lazy_exports = getattr(module, "_LAZY_EXPORTS", {})
            for name in names:
                with self.subTest(module=module_name, name=name):
                    self.assertNotIn(name, vars(module))
                    self.assertNotIn(name, exports)
                    self.assertNotIn(name, lazy_exports)

    def test_ref_panel_config_validates_r2_args(self):
        with self.assertRaises(ValueError):
            RefPanelConfig(backend="bad")
        with self.assertRaises(ValueError):
            RefPanelConfig(r2_bias_mode="bad")
        with self.assertRaises(ValueError):
            RefPanelConfig(r2_sample_size=0)

    def test_ldscore_config_requires_one_window(self):
        config = LDScoreConfig(ld_wind_cm=1.0)
        self.assertEqual(config.ld_wind_cm, 1.0)
        self.assertTrue(config.compute_m5_50)
        with self.assertRaises(ValueError):
            LDScoreConfig()
        with self.assertRaises(ValueError):
            LDScoreConfig(ld_wind_cm=1.0, ld_wind_kb=100.0)

    def test_ldscore_config_normalizes_keep_individuals_path(self):
        config = LDScoreConfig(
            ld_wind_snps=10,
            keep_indivs_path=Path("filters") / "samples.keep",
        )
        self.assertEqual(config.keep_indivs_path, "filters/samples.keep")

    def test_ref_panel_build_config_validates_required_fields(self):
        config = ReferencePanelBuildConfig(
            plink_path="plink/panel.@",
            source_genome_build="GRCh38",
            genetic_map_hg19_path="maps/hg19.txt",
            genetic_map_hg38_path="maps/hg38.txt",
            liftover_chain_hg38_to_hg19_path="liftover/hg38ToHg19.over.chain",
            output_dir="out",
            ld_wind_kb=100.0,
            overwrite=True,
        )
        self.assertEqual(config.source_genome_build, "hg38")
        self.assertEqual(config.plink_path, "plink/panel.@")
        self.assertEqual(config.genetic_map_hg19_path, "maps/hg19.txt")
        self.assertEqual(config.genetic_map_hg38_path, "maps/hg38.txt")
        self.assertEqual(config.liftover_chain_hg38_to_hg19_path, "liftover/hg38ToHg19.over.chain")
        self.assertEqual(config.output_dir, "out")
        self.assertTrue(config.overwrite)

        with self.assertRaisesRegex(ValueError, "ld_wind_cm"):
            ReferencePanelBuildConfig(
                plink_path="plink/panel",
                source_genome_build="hg19",
                genetic_map_hg19_path=None,
                genetic_map_hg38_path="maps/hg38.txt",
                liftover_chain_hg19_to_hg38_path="liftover/hg19ToHg38.over.chain",
                output_dir="out",
                ld_wind_cm=1.0,
            )
        with self.assertRaises(ValueError):
            ReferencePanelBuildConfig(
                plink_path="plink/panel",
                source_genome_build="hg19",
                genetic_map_hg19_path="maps/hg19.txt",
                genetic_map_hg38_path="maps/hg38.txt",
                liftover_chain_hg19_to_hg38_path="liftover/hg19ToHg38.over.chain",
                output_dir="out",
                ld_wind_kb=100.0,
                ld_wind_cm=1.0,
            )
        source_only = ReferencePanelBuildConfig(
            plink_path="plink/panel",
            source_genome_build="hg38",
            genetic_map_hg19_path=None,
            genetic_map_hg38_path="maps/hg38.txt",
            output_dir="out",
            ld_wind_kb=100.0,
        )
        self.assertIsNone(source_only.liftover_chain_hg38_to_hg19_path)
        self.assertIsNone(source_only.genetic_map_hg19_path)

    def test_ref_panel_build_config_normalizes_optional_paths(self):
        config = ReferencePanelBuildConfig(
            plink_path=Path("plink") / "panel",
            source_genome_build="hg19",
            genetic_map_hg19_path=Path("maps") / "hg19.txt",
            genetic_map_hg38_path=Path("maps") / "hg38.txt",
            liftover_chain_hg19_to_hg38_path=Path("liftover") / "hg19ToHg38.over.chain",
            liftover_chain_hg38_to_hg19_path=Path("liftover") / "hg38ToHg19.over.chain",
            output_dir=Path("outputs") / "ref",
            ref_panel_snps_path=Path("restrict") / "snps.txt",
            keep_indivs_path=Path("samples") / "keep.txt",
            ld_wind_snps=500,
        )
        self.assertEqual(config.plink_path, "plink/panel")
        self.assertEqual(config.genetic_map_hg19_path, "maps/hg19.txt")
        self.assertEqual(config.genetic_map_hg38_path, "maps/hg38.txt")
        self.assertEqual(config.liftover_chain_hg19_to_hg38_path, "liftover/hg19ToHg38.over.chain")
        self.assertEqual(config.liftover_chain_hg38_to_hg19_path, "liftover/hg38ToHg19.over.chain")
        self.assertEqual(config.output_dir, "outputs/ref")
        self.assertEqual(config.ref_panel_snps_path, "restrict/snps.txt")
        self.assertEqual(config.keep_indivs_path, "samples/keep.txt")

    def test_munge_config_defaults(self):
        config = MungeConfig(output_dir="out")
        self.assertEqual(config.info_min, 0.9)
        self.assertEqual(config.maf_min, 0.01)
        self.assertEqual(config.chunk_size, int(5e6))
        self.assertFalse(config.overwrite)

    def test_munge_config_normalizes_pathlike_fields(self):
        config = MungeConfig(
            sumstats_path=Path("sumstats") / "trait.tsv.gz",
            output_dir=Path("results") / "trait",
            merge_alleles_path=Path("resources") / "alleles.tsv",
            trait_name="trait",
            overwrite=True,
        )
        self.assertEqual(config.sumstats_path, "sumstats/trait.tsv.gz")
        self.assertEqual(config.output_dir, "results/trait")
        self.assertEqual(config.merge_alleles_path, "resources/alleles.tsv")
        self.assertEqual(config.trait_name, "trait")
        self.assertTrue(config.overwrite)

    def test_munge_config_accepts_source_fields(self):
        raw = MungeConfig(sumstats_path=Path("sumstats") / "trait.tsv.gz", trait_name="trait")
        self.assertIsInstance(raw, MungeConfig)
        self.assertEqual(raw.sumstats_path, "sumstats/trait.tsv.gz")
        self.assertEqual(raw.trait_name, "trait")

    def test_regression_config_defaults(self):
        config = RegressionConfig()
        self.assertEqual(config.n_blocks, 200)
        self.assertTrue(config.use_m_5_50)
        with self.assertRaises(ValueError):
            RegressionConfig(n_blocks=1)

    def test_ref_panel_config_normalizes_path_fields(self):
        config = RefPanelConfig(
            backend="parquet_r2",
            plink_path=Path("plink") / "panel",
            r2_paths=[Path("r2") / "chr1.parquet"],
            metadata_paths=[Path("freq") / "chr@.tsv.gz"],
            genome_build="GRCh37",
            ref_panel_snps_path=Path("filters") / "snps.txt",
        )
        self.assertEqual(config.backend, "parquet_r2")
        self.assertEqual(config.plink_path, "plink/panel")
        self.assertEqual(config.r2_paths, ("r2/chr1.parquet",))
        self.assertEqual(config.metadata_paths, ("freq/chr@.tsv.gz",))
        self.assertEqual(config.genome_build, "hg19")
        self.assertEqual(config.ref_panel_snps_path, "filters/snps.txt")

    def test_ref_panel_config_accepts_runtime_source_fields(self):
        config = RefPanelConfig(backend="plink", plink_path=Path("plink") / "panel")
        self.assertEqual(config.plink_path, "plink/panel")

    def test_ref_panel_config_accepts_single_string_tokens_for_plural_fields(self):
        config = RefPanelConfig(
            r2_paths="r2/reference.@.parquet",
            metadata_paths="meta/reference.*.tsv.gz",
        )
        self.assertEqual(config.r2_paths, ("r2/reference.@.parquet",))
        self.assertEqual(config.metadata_paths, ("meta/reference.*.tsv.gz",))

    def test_public_configs_normalize_pathlike_inputs(self):
        raw = MungeConfig(sumstats_path=Path("sumstats") / "trait.tsv.gz")
        annot = AnnotationBuildConfig(
            baseline_annot_paths=(Path("baseline") / "base.1.annot.gz",),
            query_annot_paths=(Path("query") / "custom.1.annot.gz",),
            query_annot_bed_paths=(Path("beds") / "enhancer.bed",),
        )
        ref = RefPanelConfig(
            backend="parquet_r2",
            plink_path=Path("plink") / "panel",
            r2_paths=(Path("r2") / "chr1.parquet",),
            metadata_paths=(Path("meta") / "chr1.tsv.gz",),
        )
        self.assertEqual(raw.sumstats_path, "sumstats/trait.tsv.gz")
        self.assertEqual(annot.baseline_annot_paths, ("baseline/base.1.annot.gz",))
        self.assertEqual(annot.query_annot_paths, ("query/custom.1.annot.gz",))
        self.assertEqual(annot.query_annot_bed_paths, ("beds/enhancer.bed",))
        self.assertEqual(ref.plink_path, "plink/panel")
        self.assertEqual(ref.r2_paths, ("r2/chr1.parquet",))
        self.assertEqual(ref.metadata_paths, ("meta/chr1.tsv.gz",))

    def test_ref_panel_config_normalizes_genome_build_aliases(self):
        self.assertEqual(RefPanelConfig(backend="parquet_r2", genome_build="hg37").genome_build, "hg19")
        self.assertEqual(RefPanelConfig(backend="parquet_r2", genome_build="GRCh37").genome_build, "hg19")
        self.assertEqual(RefPanelConfig(backend="parquet_r2", genome_build="GRCh38").genome_build, "hg38")
        self.assertEqual(RefPanelConfig(backend="parquet_r2", genome_build="auto").genome_build, "auto")

    def test_public_configs_accept_single_string_tokens_for_plural_fields(self):
        annot = AnnotationBuildConfig(
            baseline_annot_paths="baseline.@.annot.gz",
            query_annot_paths="query.*.annot.gz",
            query_annot_bed_paths="beds/*.bed",
        )
        ref = RefPanelConfig(
            backend="parquet_r2",
            r2_paths="r2/reference.@.parquet",
            metadata_paths="meta/reference.*.tsv.gz",
        )
        self.assertEqual(annot.baseline_annot_paths, ("baseline.@.annot.gz",))
        self.assertEqual(annot.query_annot_paths, ("query.*.annot.gz",))
        self.assertEqual(annot.query_annot_bed_paths, ("beds/*.bed",))
        self.assertEqual(ref.r2_paths, ("r2/reference.@.parquet",))
        self.assertEqual(ref.metadata_paths, ("meta/reference.*.tsv.gz",))

    def test_removed_public_config_fields_are_rejected(self):
        with self.assertRaises(TypeError):
            LDScoreConfig(ld_wind_snps=10, keep_individuals_path="samples.keep")
        with self.assertRaises(TypeError):
            ReferencePanelBuildConfig(
                panel_label="EUR",
                plink_prefix="plink/panel",
                source_genome_build="hg19",
                genetic_map_hg19_path="maps/hg19.txt",
                genetic_map_hg38_path="maps/hg38.txt",
                liftover_chain_hg19_to_hg38_path="liftover/hg19ToHg38.over.chain",
                output_dir="out",
                ld_wind_kb=100.0,
            )
        with self.assertRaises(TypeError):
            MungeConfig(out_prefix="out")
        with self.assertRaises(TypeError):
            RefPanelConfig(plink_prefix="plink/panel")


class IdentifierHelpersTest(unittest.TestCase):
    def test_clean_header(self):
        self.assertEqual(clean_header("foo-bar.foo_BaR"), "FOO_BAR_FOO_BAR")

    def test_normalize_snp_identifier_mode(self):
        for value in ["rsid", "rsID", "SNPID", "snp_id", "snp"]:
            self.assertEqual(normalize_snp_identifier_mode(value), "rsid")
        for value in ["chr_pos", "ChrPos", "chrom_pos"]:
            self.assertEqual(normalize_snp_identifier_mode(value), "chr_pos")

    def test_infer_snp_column_aliases(self):
        for header in [["rsid"], ["rsID"], ["SNPID"], ["snp_id"], ["SNP"], ["id"]]:
            self.assertEqual(infer_snp_column(header), header[0])

    def test_infer_chr_bp_columns(self):
        self.assertEqual(infer_chr_bp_columns(["chromosome", "position"]), ("chromosome", "position"))

    def test_normalize_chromosome(self):
        self.assertEqual(normalize_chromosome("chr01"), "1")
        self.assertEqual(normalize_chromosome("1.0"), "1")
        self.assertEqual(normalize_chromosome("X"), "X")
        self.assertEqual(normalize_chromosome("M"), "M")
        self.assertEqual(normalize_chromosome("MT"), "MT")

    def test_normalize_chromosome_warns_on_sex_chromosome_numeric_inference(self):
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            value = normalize_chromosome("23", context="test-sex-code")

        self.assertEqual(value, "X")
        self.assertEqual(len(caught), 1)
        self.assertIn("sex chromosome", str(caught[0].message).lower())
        self.assertIn("23", str(caught[0].message))
        self.assertIn("X", str(caught[0].message))

    def test_normalize_chromosome_rejects_numeric_mitochondrial_codes(self):
        for value in ["25", "26", "25.0", "26.0"]:
            with self.assertRaisesRegex(ValueError, "M' or 'MT"):
                normalize_chromosome(value, context="test-mito-code")

    def test_build_chr_pos_snp_id(self):
        self.assertEqual(build_chr_pos_snp_id("chr1", 123), "1:123")

    def test_build_snp_id_series(self):
        rsid_df = pd.DataFrame({"SNP": ["rs1", "rs2"]})
        self.assertEqual(build_snp_id_series(rsid_df, "rsid").tolist(), ["rs1", "rs2"])
        chr_pos_df = pd.DataFrame({"CHR": ["chr1", "2"], "POS": [10, 20]})
        self.assertEqual(build_snp_id_series(chr_pos_df, "chr_pos").tolist(), ["1:10", "2:20"])

    def test_build_snp_id_series_accepts_bp_alias(self):
        chr_pos_df = pd.DataFrame({"CHR": ["chr1", "2"], "BP": [10, 20]})
        self.assertEqual(build_snp_id_series(chr_pos_df, "chr_pos").tolist(), ["1:10", "2:20"])

    def test_validate_unique_snp_ids(self):
        df = pd.DataFrame({"SNP": ["rs1", "rs1"]})
        with self.assertRaises(ValueError):
            validate_unique_snp_ids(df, "rsid", context="test")


class RestrictionReadersTest(unittest.TestCase):
    def test_read_global_snp_restriction_rsid_requires_header(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "restrict.txt"
            path.write_text("rs1\nrs2\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "header"):
                read_global_snp_restriction(path, "rsid")

    def test_read_global_snp_restriction_rsid_table(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "restrict.tsv"
            path.write_text("id\tother\nrs1\ta\nrs2\tb\n", encoding="utf-8")
            self.assertEqual(read_global_snp_restriction(path, "rsID"), {"rs1", "rs2"})

    def test_read_global_snp_restriction_chr_pos_table(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "restrict.tsv"
            path.write_text("chrom\tbp\nchr1.0\t10\n2\t20\n", encoding="utf-8")
            self.assertEqual(read_global_snp_restriction(path, "chr_pos"), {"1:10", "2:20"})

    def test_read_global_snp_restriction_chr_pos_fixed_build_skips_auto_inference(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "restrict.tsv"
            path.write_text("CHR\tBP\n1\t20\n2\t30\n", encoding="utf-8")
            with mock.patch("ldsc._kernel.identifiers.resolve_chr_pos_table") as patched:
                self.assertEqual(
                    read_global_snp_restriction(path, "chr_pos", genome_build="hg38"),
                    {"1:20", "2:30"},
                )
            patched.assert_not_called()

    def test_read_global_snp_restriction_chr_pos_uses_build_specific_position_header(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "restrict.tsv"
            path.write_text(
                "rsID\tA1\tchr\thg19_pos\thg38_pos\n"
                "rs1\tA\t1\t10.0\t100.0\n"
                "rs2\tG\tchr2\t20.0\t200.0\n",
                encoding="utf-8",
            )
            self.assertEqual(
                read_global_snp_restriction(path, "chr_pos", genome_build="hg19"),
                {"1:10", "2:20"},
            )

    def test_read_global_snp_restriction_chr_pos_drops_rows_missing_chr_or_pos(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "restrict.tsv"
            path.write_text(
                "rsID\tA1\tchr\thg19_pos\thg38_pos\n"
                "rs1\tA\t1\t10.0\t100.0\n"
                "rs2\tG\t2\t\t200.0\n"
                "rs3\tT\t\t30.0\t300.0\n",
                encoding="utf-8",
            )
            self.assertEqual(
                read_global_snp_restriction(path, "chr_pos", genome_build="hg19"),
                {"1:10"},
            )

    def test_read_global_snp_restriction_chr_pos_requires_header(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "restrict.tsv"
            path.write_text("1\t10\n2\t20\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "header"):
                read_global_snp_restriction(path, "chr_pos")

    def test_read_global_snp_restriction_rsid_missing_alias_raises_clear_error(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "restrict.tsv"
            path.write_text("variant\tother\nrs1\ta\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "SNP"):
                read_global_snp_restriction(path, "rsid")

    def test_read_global_snp_restriction_chr_pos_missing_alias_raises_clear_error(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "restrict.tsv"
            path.write_text("variant\tother\nrs1\ta\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "CHR, POS|POS, CHR"):
                read_global_snp_restriction(path, "chr_pos")
