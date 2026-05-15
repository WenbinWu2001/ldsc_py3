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
from ldsc.genome_build_inference import resolve_genome_build, validate_auto_genome_build_mode
from ldsc._kernel.identifiers import (
    build_chr_pos_snp_id,
    build_snp_id_series,
    clean_header,
    infer_chr_bp_columns,
    infer_snp_column,
    normalize_chromosome,
    normalize_snp_identifier_mode,
    read_global_snp_restriction,
    read_snp_restriction_keys,
    validate_unique_snp_ids,
)
from ldsc.column_inference import A1_COLUMN_SPEC, A2_COLUMN_SPEC, RESTRICTION_ALLELE_SPECS
from ldsc._row_alignment import assert_same_snp_rows


class GlobalConfigTest(unittest.TestCase):
    def test_defaults(self):
        config = GlobalConfig()
        self.assertEqual(config.snp_identifier, "chr_pos_allele_aware")
        self.assertEqual(config.genome_build, "auto")
        self.assertEqual(config.log_level, "INFO")
        self.assertFalse(config.fail_on_missing_metadata)

    def test_package_global_registry_round_trip(self):
        original = ldsc.reset_global_config()
        self.assertEqual(original, GlobalConfig())

        configured = GlobalConfig(
            snp_identifier="rsid",
            log_level="DEBUG",
        )
        ldsc.set_global_config(configured)

        self.assertEqual(ldsc.get_global_config(), configured)

    def test_reset_global_config_restores_default(self):
        ldsc.set_global_config(GlobalConfig(snp_identifier="rsid", log_level="DEBUG"))

        reset = ldsc.reset_global_config()

        self.assertEqual(reset, GlobalConfig())
        self.assertEqual(ldsc.get_global_config(), GlobalConfig())

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
    def test_chr_pos_omitted_genome_build_defaults_to_auto(self):
        for mode in ("chr_pos", "chr_pos_allele_aware"):
            with self.subTest(mode=mode):
                cfg = GlobalConfig(snp_identifier=mode)
                self.assertEqual(cfg.genome_build, "auto")

    def test_chr_pos_explicit_none_raises(self):
        for mode in ("chr_pos", "chr_pos_allele_aware"):
            with self.subTest(mode=mode):
                with self.assertRaisesRegex(ValueError, "Pass genome_build='auto'"):
                    GlobalConfig(snp_identifier=mode, genome_build=None)

    def test_chr_pos_hg38_ok(self):
        cfg = GlobalConfig(snp_identifier="chr_pos", genome_build="hg38")
        self.assertEqual(cfg.genome_build, "hg38")

    def test_chr_pos_auto_ok(self):
        cfg = GlobalConfig(snp_identifier="chr_pos", genome_build="auto")
        self.assertEqual(cfg.genome_build, "auto")

    def test_rsid_auto_raises(self):
        for mode in ("rsid", "rsid_allele_aware"):
            with self.subTest(mode=mode):
                with self.assertRaises(ValueError):
                    GlobalConfig(snp_identifier=mode, genome_build="auto")

    def test_rsid_no_genome_build_ok(self):
        for mode in ("rsid", "rsid_allele_aware"):
            with self.subTest(mode=mode):
                cfg = GlobalConfig(snp_identifier=mode)
                self.assertIsNone(cfg.genome_build)

    def test_rsid_with_concrete_build_warns_and_nulls(self):
        for mode in ("rsid", "rsid_allele_aware"):
            with self.subTest(mode=mode):
                with self.assertWarns(UserWarning) as ctx:
                    cfg = GlobalConfig(snp_identifier=mode, genome_build="hg38")
                self.assertIn("ignored", str(ctx.warning))
                self.assertIsNone(cfg.genome_build)

    def test_singleton_is_chr_pos_allele_aware_auto(self):
        from ldsc.config import get_global_config
        cfg = get_global_config()
        self.assertEqual(cfg.snp_identifier, "chr_pos_allele_aware")
        self.assertEqual(cfg.genome_build, "auto")

    def test_reset_returns_chr_pos_allele_aware_auto(self):
        from ldsc.config import reset_global_config
        cfg = reset_global_config()
        self.assertEqual(cfg.snp_identifier, "chr_pos_allele_aware")
        self.assertEqual(cfg.genome_build, "auto")

    def test_auto_genome_build_validation_uses_identifier_family(self):
        validate_auto_genome_build_mode("chr_pos_allele_aware", "auto")
        with self.assertRaisesRegex(ValueError, "chr_pos"):
            validate_auto_genome_build_mode("rsid_allele_aware", "auto")

    def test_resolve_genome_build_treats_rsid_family_as_buildless(self):
        self.assertIsNone(resolve_genome_build("hg38", "rsid_allele_aware", None, context="test"))

    def test_row_alignment_uses_identifier_family(self):
        left = pd.DataFrame({"CHR": ["1"], "POS": [10], "SNP": ["rs1"]})
        right_chr_pos = pd.DataFrame({"CHR": ["1"], "POS": [10]})
        assert_same_snp_rows(left, right_chr_pos, context="test", snp_identifier="chr_pos_allele_aware")

        right_rsid = pd.DataFrame({"CHR": ["1"], "POS": [10], "SNP": ["rs2"]})
        with self.assertRaisesRegex(ValueError, "SNP"):
            assert_same_snp_rows(left, right_rsid, context="test", snp_identifier="rsid_allele_aware")


class WorkflowConfigTest(unittest.TestCase):
    def test_annotation_config_normalizes_paths(self):
        config = AnnotationBuildConfig(
            baseline_annot_sources=[Path("base.annot.gz")],
            query_annot_bed_sources=[Path("query.bed")],
            output_dir=Path("outputs") / "annot",
            overwrite=True,
        )
        self.assertEqual(config.baseline_annot_sources, ("base.annot.gz",))
        self.assertEqual(config.query_annot_bed_sources, ("query.bed",))
        self.assertEqual(config.output_dir, "outputs/annot")
        self.assertFalse(hasattr(config, "batch_mode"))
        self.assertTrue(config.overwrite)

    def test_annotation_config_accepts_single_string_tokens_for_plural_fields(self):
        config = AnnotationBuildConfig(
            baseline_annot_sources="baseline.@.annot.gz",
            query_annot_sources="query.*.annot.gz",
            query_annot_bed_sources="beds/*.bed",
        )
        self.assertEqual(config.baseline_annot_sources, ("baseline.@.annot.gz",))
        self.assertEqual(config.query_annot_sources, ("query.*.annot.gz",))
        self.assertEqual(config.query_annot_bed_sources, ("beds/*.bed",))

    def test_removed_compatibility_aliases_are_not_exported(self):
        removed_names_by_module = {
            "ldsc": ("AnnotationSourceSpec", "RawSumstatsSpec", "RefPanelSpec", "OutputSpec"),
            "ldsc.annotation_builder": ("AnnotationSourceSpec", "main_bed_to_annot"),
            "ldsc._kernel.annotation": (
                "AnnotationSourceSpec",
                "AnnotationBuilder",
                "AnnotationBundle",
                "main_bed_to_annot",
                "parse_bed_to_annot_args",
                "run_bed_to_annot",
            ),
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
        self.assertEqual(config.snp_batch_size, 128)
        self.assertFalse(hasattr(config, "chunk_size"))
        self.assertEqual(config.common_maf_min, 0.05)
        with self.assertRaises(ValueError):
            LDScoreConfig()
        with self.assertRaises(ValueError):
            LDScoreConfig(ld_wind_cm=1.0, ld_wind_kb=100.0)

    def test_ldscore_config_accepts_snp_batch_size_and_rejects_chunk_size(self):
        config = LDScoreConfig(ld_wind_snps=10, snp_batch_size=64)
        self.assertEqual(config.snp_batch_size, 64)
        with self.assertRaises(TypeError):
            LDScoreConfig(ld_wind_snps=10, chunk_size=64)

    def test_ldscore_config_accepts_hm3_regression_flag_and_rejects_conflict(self):
        config = LDScoreConfig(ld_wind_snps=10, use_hm3_regression_snps=True)

        self.assertTrue(config.use_hm3_regression_snps)
        with self.assertRaisesRegex(ValueError, "regression_snps_file.*use_hm3_regression_snps"):
            LDScoreConfig(
                ld_wind_snps=10,
                regression_snps_file="custom.tsv",
                use_hm3_regression_snps=True,
            )

    def test_ldscore_config_rejects_reference_panel_filter_fields(self):
        with self.assertRaises(TypeError):
            LDScoreConfig(ld_wind_snps=10, maf_min=0.01)
        with self.assertRaises(TypeError):
            LDScoreConfig(ld_wind_snps=10, keep_indivs_file=Path("filters") / "samples.keep")

    def test_ref_panel_config_normalizes_reference_panel_filter_fields(self):
        config = RefPanelConfig(
            backend="plink",
            maf_min=0.01,
            keep_indivs_file=Path("filters") / "samples.keep",
        )
        self.assertEqual(config.maf_min, 0.01)
        self.assertEqual(config.keep_indivs_file, "filters/samples.keep")

    def test_ref_panel_build_config_validates_required_fields(self):
        config = ReferencePanelBuildConfig(
            plink_prefix="plink/panel.@",
            source_genome_build="GRCh38",
            genetic_map_hg19_sources="maps/hg19.txt",
            genetic_map_hg38_sources="maps/hg38.txt",
            liftover_chain_hg38_to_hg19_file="liftover/hg38ToHg19.over.chain",
            output_dir="out",
            ld_wind_kb=100.0,
            overwrite=True,
        )
        self.assertEqual(config.source_genome_build, "hg38")
        self.assertEqual(config.plink_prefix, "plink/panel.@")
        self.assertEqual(config.genetic_map_hg19_sources, "maps/hg19.txt")
        self.assertEqual(config.genetic_map_hg38_sources, "maps/hg38.txt")
        self.assertEqual(config.liftover_chain_hg38_to_hg19_file, "liftover/hg38ToHg19.over.chain")
        self.assertEqual(config.output_dir, "out")
        self.assertEqual(config.snp_batch_size, 128)
        self.assertTrue(config.overwrite)

        with self.assertRaisesRegex(ValueError, "ld_wind_cm"):
            ReferencePanelBuildConfig(
                plink_prefix="plink/panel",
                source_genome_build="hg19",
                genetic_map_hg19_sources=None,
                genetic_map_hg38_sources="maps/hg38.txt",
                liftover_chain_hg19_to_hg38_file="liftover/hg19ToHg38.over.chain",
                output_dir="out",
                ld_wind_cm=1.0,
            )
        with self.assertRaises(ValueError):
            ReferencePanelBuildConfig(
                plink_prefix="plink/panel",
                source_genome_build="hg19",
                genetic_map_hg19_sources="maps/hg19.txt",
                genetic_map_hg38_sources="maps/hg38.txt",
                liftover_chain_hg19_to_hg38_file="liftover/hg19ToHg38.over.chain",
                output_dir="out",
                ld_wind_kb=100.0,
                ld_wind_cm=1.0,
            )
        source_only = ReferencePanelBuildConfig(
            plink_prefix="plink/panel",
            source_genome_build="hg38",
            genetic_map_hg19_sources=None,
            genetic_map_hg38_sources="maps/hg38.txt",
            output_dir="out",
            ld_wind_kb=100.0,
        )
        self.assertIsNone(source_only.liftover_chain_hg38_to_hg19_file)
        self.assertIsNone(source_only.genetic_map_hg19_sources)

    def test_ref_panel_build_config_normalizes_optional_paths(self):
        config = ReferencePanelBuildConfig(
            plink_prefix=Path("plink") / "panel",
            source_genome_build="hg19",
            genetic_map_hg19_sources=Path("maps") / "hg19.txt",
            genetic_map_hg38_sources=Path("maps") / "hg38.txt",
            liftover_chain_hg19_to_hg38_file=Path("liftover") / "hg19ToHg38.over.chain",
            liftover_chain_hg38_to_hg19_file=Path("liftover") / "hg38ToHg19.over.chain",
            output_dir=Path("outputs") / "ref",
            ref_panel_snps_file=Path("restrict") / "snps.txt",
            keep_indivs_file=Path("samples") / "keep.txt",
            ld_wind_snps=500,
        )
        self.assertEqual(config.plink_prefix, "plink/panel")
        self.assertEqual(config.genetic_map_hg19_sources, "maps/hg19.txt")
        self.assertEqual(config.genetic_map_hg38_sources, "maps/hg38.txt")
        self.assertEqual(config.liftover_chain_hg19_to_hg38_file, "liftover/hg19ToHg38.over.chain")
        self.assertEqual(config.liftover_chain_hg38_to_hg19_file, "liftover/hg38ToHg19.over.chain")
        self.assertEqual(config.output_dir, "outputs/ref")
        self.assertEqual(config.ref_panel_snps_file, "restrict/snps.txt")
        self.assertEqual(config.keep_indivs_file, "samples/keep.txt")

    def test_ref_panel_build_config_accepts_hm3_flags_and_rejects_conflicts(self):
        config = ReferencePanelBuildConfig(
            plink_prefix="plink/panel",
            source_genome_build="hg19",
            output_dir="out",
            ld_wind_snps=500,
            use_hm3_snps=True,
            use_hm3_quick_liftover=True,
        )

        self.assertTrue(config.use_hm3_snps)
        self.assertTrue(config.use_hm3_quick_liftover)
        with self.assertRaisesRegex(ValueError, "ref_panel_snps_file.*use_hm3_snps"):
            ReferencePanelBuildConfig(
                plink_prefix="plink/panel",
                source_genome_build="hg19",
                output_dir="out",
                ld_wind_snps=500,
                ref_panel_snps_file="custom.tsv",
                use_hm3_snps=True,
            )
        with self.assertRaisesRegex(ValueError, "use_hm3_snps"):
            ReferencePanelBuildConfig(
                plink_prefix="plink/panel",
                source_genome_build="hg19",
                output_dir="out",
                ld_wind_snps=500,
                use_hm3_quick_liftover=True,
            )
        with self.assertRaisesRegex(ValueError, "mutually exclusive"):
            ReferencePanelBuildConfig(
                plink_prefix="plink/panel",
                source_genome_build="hg19",
                output_dir="out",
                ld_wind_snps=500,
                use_hm3_snps=True,
                use_hm3_quick_liftover=True,
                liftover_chain_hg19_to_hg38_file="chain.over",
            )

    def test_munge_config_defaults(self):
        config = MungeConfig(output_dir="out")
        self.assertEqual(config.info_min, 0.9)
        self.assertEqual(config.maf_min, 0.01)
        self.assertEqual(config.chunk_size, 1_000_000)
        self.assertEqual(config.output_format, "parquet")
        self.assertFalse(config.overwrite)

    def test_munge_config_rejects_unknown_output_format(self):
        with self.assertRaisesRegex(ValueError, "output_format"):
            MungeConfig(output_dir="out", output_format="csv")

    def test_munge_config_normalizes_pathlike_fields(self):
        config = MungeConfig(
            raw_sumstats_file=Path("sumstats") / "trait.tsv.gz",
            output_dir=Path("results") / "trait",
            sumstats_snps_file=Path("resources") / "alleles.tsv",
            trait_name="trait",
            overwrite=True,
        )
        self.assertEqual(config.raw_sumstats_file, "sumstats/trait.tsv.gz")
        self.assertEqual(config.output_dir, "results/trait")
        self.assertEqual(config.sumstats_snps_file, "resources/alleles.tsv")
        self.assertEqual(config.trait_name, "trait")
        self.assertTrue(config.overwrite)

    def test_munge_config_accepts_hm3_flags_and_rejects_conflicts(self):
        config = MungeConfig(
            output_dir="out",
            target_genome_build="hg38",
            use_hm3_snps=True,
            use_hm3_quick_liftover=True,
        )

        self.assertTrue(config.use_hm3_snps)
        self.assertTrue(config.use_hm3_quick_liftover)
        with self.assertRaisesRegex(ValueError, "sumstats_snps_file.*use_hm3_snps"):
            MungeConfig(output_dir="out", sumstats_snps_file="custom.tsv", use_hm3_snps=True)
        with self.assertRaisesRegex(ValueError, "use_hm3_snps"):
            MungeConfig(output_dir="out", target_genome_build="hg38", use_hm3_quick_liftover=True)

    def test_munge_config_normalizes_trait_name(self):
        config = MungeConfig(raw_sumstats_file="sumstats/trait.tsv.gz", trait_name=" MDD ")

        self.assertEqual(config.trait_name, "MDD")

    def test_munge_config_rejects_blank_trait_name(self):
        with self.assertRaisesRegex(ValueError, "trait_name"):
            MungeConfig(raw_sumstats_file="sumstats/trait.tsv.gz", trait_name="  ")

    def test_munge_config_accepts_source_fields(self):
        raw = MungeConfig(raw_sumstats_file=Path("sumstats") / "trait.tsv.gz", trait_name="trait")
        self.assertIsInstance(raw, MungeConfig)
        self.assertEqual(raw.raw_sumstats_file, "sumstats/trait.tsv.gz")
        self.assertEqual(raw.trait_name, "trait")

    def test_regression_config_defaults(self):
        config = RegressionConfig()
        self.assertEqual(config.n_blocks, 200)
        self.assertTrue(config.use_common_counts)
        with self.assertRaises(ValueError):
            RegressionConfig(n_blocks=1)

    def test_ref_panel_config_normalizes_path_fields(self):
        config = RefPanelConfig(
            backend="parquet_r2",
            plink_prefix=Path("plink") / "panel",
            r2_dir=Path("r2_panel") / "hg38",
            ref_panel_snps_file=Path("filters") / "snps.txt",
            keep_indivs_file=Path("filters") / "samples.keep",
            maf_min=0.02,
        )
        self.assertEqual(config.backend, "parquet_r2")
        self.assertEqual(config.plink_prefix, "plink/panel")
        self.assertEqual(config.r2_dir, "r2_panel/hg38")
        self.assertEqual(config.ref_panel_snps_file, "filters/snps.txt")
        self.assertEqual(config.keep_indivs_file, "filters/samples.keep")
        self.assertEqual(config.maf_min, 0.02)

    def test_ref_panel_config_accepts_hm3_flag_and_rejects_conflict(self):
        config = RefPanelConfig(backend="plink", use_hm3_ref_panel_snps=True)

        self.assertTrue(config.use_hm3_ref_panel_snps)
        with self.assertRaisesRegex(ValueError, "ref_panel_snps_file.*use_hm3_ref_panel_snps"):
            RefPanelConfig(
                backend="plink",
                ref_panel_snps_file="custom.tsv",
                use_hm3_ref_panel_snps=True,
            )

    def test_ref_panel_config_accepts_runtime_source_fields(self):
        config = RefPanelConfig(backend="plink", plink_prefix=Path("plink") / "panel")
        self.assertEqual(config.plink_prefix, "plink/panel")

    def test_ref_panel_config_rejects_removed_parquet_source_fields(self):
        for kwargs in (
            {"r2_sources": "r2/reference.@.parquet"},
            {"metadata_sources": "meta/reference.*.tsv.gz"},
        ):
            with self.subTest(kwargs=kwargs):
                with self.assertRaises(TypeError):
                    RefPanelConfig(**kwargs)

    def test_public_configs_normalize_pathlike_inputs(self):
        raw = MungeConfig(raw_sumstats_file=Path("sumstats") / "trait.tsv.gz")
        annot = AnnotationBuildConfig(
            baseline_annot_sources=(Path("baseline") / "base.1.annot.gz",),
            query_annot_sources=(Path("query") / "custom.1.annot.gz",),
            query_annot_bed_sources=(Path("beds") / "enhancer.bed",),
        )
        ref = RefPanelConfig(
            backend="parquet_r2",
            plink_prefix=Path("plink") / "panel",
            r2_dir=Path("r2_panel") / "hg38",
        )
        self.assertEqual(raw.raw_sumstats_file, "sumstats/trait.tsv.gz")
        self.assertEqual(annot.baseline_annot_sources, ("baseline/base.1.annot.gz",))
        self.assertEqual(annot.query_annot_sources, ("query/custom.1.annot.gz",))
        self.assertEqual(annot.query_annot_bed_sources, ("beds/enhancer.bed",))
        self.assertEqual(ref.plink_prefix, "plink/panel")
        self.assertEqual(ref.r2_dir, "r2_panel/hg38")

    def test_ref_panel_config_has_no_genome_build_field(self):
        config = RefPanelConfig(backend="parquet_r2")
        self.assertFalse(hasattr(config, "genome_build"))
        with self.assertRaises(TypeError):
            RefPanelConfig(backend="parquet_r2", genome_build="hg38")

    def test_public_configs_accept_single_string_tokens_for_plural_fields(self):
        annot = AnnotationBuildConfig(
            baseline_annot_sources="baseline.@.annot.gz",
            query_annot_sources="query.*.annot.gz",
            query_annot_bed_sources="beds/*.bed",
        )
        ref = RefPanelConfig(backend="parquet_r2", r2_dir="r2_panel/hg38")
        self.assertEqual(annot.baseline_annot_sources, ("baseline.@.annot.gz",))
        self.assertEqual(annot.query_annot_sources, ("query.*.annot.gz",))
        self.assertEqual(annot.query_annot_bed_sources, ("beds/*.bed",))
        self.assertEqual(ref.r2_dir, "r2_panel/hg38")

    def test_removed_public_config_fields_are_rejected(self):
        with self.assertRaises(TypeError):
            LDScoreConfig(ld_wind_snps=10, keep_individuals_path="samples.keep")
        with self.assertRaises(TypeError):
            ReferencePanelBuildConfig(
                panel_label="EUR",
                plink_prefix="plink/panel",
                source_genome_build="hg19",
                genetic_map_hg19_sources="maps/hg19.txt",
                genetic_map_hg38_sources="maps/hg38.txt",
                liftover_chain_hg19_to_hg38_file="liftover/hg19ToHg38.over.chain",
                output_dir="out",
                ld_wind_kb=100.0,
            )
        with self.assertRaises(TypeError):
            MungeConfig(out_prefix="out")
        with self.assertRaises(TypeError):
            MungeConfig(sumstats_file="sumstats/trait.tsv.gz")
        with self.assertRaises(TypeError):
            RefPanelConfig(plink_path="plink/panel")


class IdentifierHelpersTest(unittest.TestCase):
    def test_clean_header(self):
        self.assertEqual(clean_header("foo-bar.foo_BaR"), "FOO_BAR_FOO_BAR")

    def test_normalize_snp_identifier_mode_accepts_only_public_modes(self):
        for value in ("rsid", "rsid_allele_aware", "chr_pos", "chr_pos_allele_aware"):
            self.assertEqual(normalize_snp_identifier_mode(value), value)
        for value in ("rsID", "SNPID", "snp_id", "snp", "chrpos", "rsid_alleles", "chr_pos_alleles"):
            with self.assertRaises(ValueError):
                normalize_snp_identifier_mode(value)

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

    def test_normalize_chromosome_logs_on_sex_chromosome_numeric_inference(self):
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            with self.assertLogs("LDSC.chromosomes", level="INFO") as captured:
                value = normalize_chromosome("23", context="test-sex-code")

        self.assertEqual(value, "X")
        self.assertEqual(caught, [])
        self.assertEqual(captured.records[0].levelname, "INFO")
        self.assertIn("sex chromosome", captured.output[0].lower())
        self.assertIn("23", captured.output[0])
        self.assertIn("X", captured.output[0])

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

    def test_validate_unique_chr_pos_ids_reports_colliding_snps(self):
        df = pd.DataFrame(
            {
                "CHR": ["1", "1", "1"],
                "POS": [10, 10, 20],
                "SNP": ["rsA", "rsB", "rsC"],
            }
        )
        with self.assertRaisesRegex(ValueError, "1:10.*rsA.*rsB.*prune duplicate-coordinate variants"):
            validate_unique_snp_ids(df, "chr_pos", context="test")


class RestrictionReadersTest(unittest.TestCase):
    def test_restriction_allele_specs_treat_ref_alt_as_external_aliases(self):
        self.assertEqual(A1_COLUMN_SPEC.canonical, "A1")
        self.assertEqual(A2_COLUMN_SPEC.canonical, "A2")
        self.assertEqual(RESTRICTION_ALLELE_SPECS, (A1_COLUMN_SPEC, A2_COLUMN_SPEC))
        self.assertIn("REF", A1_COLUMN_SPEC.aliases)
        self.assertIn("ALT", A2_COLUMN_SPEC.aliases)

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
            self.assertEqual(read_global_snp_restriction(path, "rsid"), {"rs1", "rs2"})

    def test_read_global_snp_restriction_rsid_family_table(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "restrict.tsv"
            path.write_text("id\tA1\tA2\nrs1\tA\tC\nrs2\tG\tT\n", encoding="utf-8")
            self.assertEqual(read_global_snp_restriction(path, "rsid_allele_aware"), {"rs1", "rs2"})

    def test_read_global_snp_restriction_chr_pos_table(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "restrict.tsv"
            path.write_text("chrom\tbp\nchr1.0\t10\n2\t20\n", encoding="utf-8")
            self.assertEqual(read_global_snp_restriction(path, "chr_pos"), {"1:10", "2:20"})

    def test_read_global_snp_restriction_chr_pos_fixed_build_reads_directly(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "restrict.tsv"
            path.write_text("CHR\tBP\n1\t20\n2\t30\n", encoding="utf-8")
            self.assertEqual(
                read_global_snp_restriction(path, "chr_pos", genome_build="hg38"),
                {"1:20", "2:30"},
            )

    def test_read_global_snp_restriction_rejects_auto_genome_build_at_kernel_boundary(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "restrict.tsv"
            path.write_text("CHR\tBP\n1\t20\n2\t30\n", encoding="utf-8")
            with self.assertRaisesRegex(AssertionError, "workflow entry"):
                read_global_snp_restriction(path, "chr_pos", genome_build="auto")

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

    def test_read_global_snp_restriction_chr_pos_treats_dot_as_missing(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "restrict.tsv"
            path.write_text(
                "CHR\tPOS\n"
                "1\t10\n"
                ".\t20\n"
                "2\t.\n",
                encoding="utf-8",
            )
            with self.assertLogs("LDSC.identifiers", level="WARNING") as caught:
                values = read_global_snp_restriction(path, "chr_pos", genome_build="hg19")

            self.assertEqual(values, {"1:10"})
            log_text = "\n".join(caught.output)
            self.assertIn("Dropped 2 SNPs with invalid or missing CHR/POS", log_text)
            self.assertIn("missing CHR=1", log_text)
            self.assertIn("missing POS=1", log_text)

    def test_read_global_snp_restriction_chr_pos_drops_invalid_non_missing_pos(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "restrict.tsv"
            path.write_text("CHR\tPOS\n1\tabc\n", encoding="utf-8")
            with self.assertLogs("LDSC.identifiers", level="WARNING") as caught:
                values = read_global_snp_restriction(path, "chr_pos", genome_build="hg19")

            self.assertEqual(values, set())
            log_text = "\n".join(caught.output)
            self.assertIn("Dropped 1 SNPs with invalid or missing CHR/POS", log_text)
            self.assertIn("invalid POS=1", log_text)

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

    def test_read_snp_restriction_keys_without_alleles_uses_base_match_in_allele_aware_mode(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "restrict.tsv"
            path.write_text("SNP\nrs1\nrs2\nrs1\n", encoding="utf-8")

            restriction = read_snp_restriction_keys(path, "rsid_allele_aware")

            self.assertEqual(restriction.match_kind, "base")
            self.assertEqual(restriction.keys, {"rs1", "rs2"})
            self.assertEqual(restriction.n_input_rows, 3)
            self.assertEqual(restriction.n_retained_keys, 2)
            self.assertTrue(restriction.dropped.empty)

    def test_read_snp_restriction_keys_with_alleles_uses_identity_match_in_allele_aware_mode(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "restrict.tsv"
            path.write_text("SNP\tREF\tALT\nrs1\tA\tC\nrs2\tT\tG\n", encoding="utf-8")

            restriction = read_snp_restriction_keys(path, "rsid_allele_aware")

            self.assertEqual(restriction.match_kind, "identity")
            self.assertEqual(restriction.keys, {"rs1:A:C", "rs2:A:C"})

    def test_read_snp_restriction_keys_requires_both_allele_columns_when_one_resolves(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "restrict.tsv"
            path.write_text("SNP\tREF\nrs1\tA\n", encoding="utf-8")

            with self.assertRaisesRegex(ValueError, "both A1 and A2"):
                read_snp_restriction_keys(path, "rsid_allele_aware")

    def test_read_snp_restriction_keys_drops_invalid_allele_rows_with_warning(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "restrict.tsv"
            path.write_text("SNP\tA1\tA2\nrs1\tA\tC\nrs_bad\t\tG\n", encoding="utf-8")

            with self.assertLogs("LDSC.identifiers", level="WARNING") as caught:
                restriction = read_snp_restriction_keys(path, "rsid_allele_aware")

            self.assertEqual(restriction.match_kind, "identity")
            self.assertEqual(restriction.keys, {"rs1:A:C"})
            self.assertEqual(restriction.dropped["reason"].tolist(), ["missing_allele"])
            self.assertIn("Dropped 1 SNP identity rows", "\n".join(caught.output))

    def test_read_snp_restriction_keys_collapses_duplicate_identity_keys(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "restrict.tsv"
            path.write_text("SNP\tA1\tA2\nrs1\tA\tC\nrs1\tC\tA\n", encoding="utf-8")

            restriction = read_snp_restriction_keys(path, "rsid_allele_aware")

            self.assertEqual(restriction.match_kind, "identity")
            self.assertEqual(restriction.keys, {"rs1:A:C"})
            self.assertEqual(restriction.n_retained_keys, 1)

    def test_read_snp_restriction_keys_drops_multiallelic_restriction_clusters(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "restrict.tsv"
            path.write_text("CHR\tPOS\tA1\tA2\n1\t10\tA\tC\n1\t10\tA\tG\n2\t20\tA\tC\n", encoding="utf-8")

            with self.assertLogs("LDSC.identifiers", level="WARNING"):
                restriction = read_snp_restriction_keys(path, "chr_pos_allele_aware")

            self.assertEqual(restriction.match_kind, "identity")
            self.assertEqual(restriction.keys, {"2:20:A:C"})
            self.assertEqual(restriction.dropped["reason"].tolist(), ["multi_allelic_base_key", "multi_allelic_base_key"])

    def test_read_snp_restriction_keys_base_modes_ignore_bad_allele_columns(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "restrict.tsv"
            path.write_text(
                "CHR\tPOS\tA1\tA2\n"
                "1\t10\tA\tT\n"
                "1\t10\tC\tG\n"
                "2\t20\t\tAT\n",
                encoding="utf-8",
            )

            restriction = read_snp_restriction_keys(path, "chr_pos")

            self.assertEqual(restriction.match_kind, "base")
            self.assertEqual(restriction.keys, {"1:10", "2:20"})
            self.assertTrue(restriction.dropped.empty)
