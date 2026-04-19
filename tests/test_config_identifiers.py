from pathlib import Path
import sys
import tempfile
import unittest

import pandas as pd

SRC = Path(__file__).resolve().parents[1] / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from ldsc.config import (
    AnnotationBuildConfig,
    CommonConfig,
    LDScoreConfig,
    MungeConfig,
    RefPanelConfig,
    RegressionConfig,
    normalize_genome_build,
)
from ldsc.sumstats_munger import RawSumstatsSpec
from ldsc._kernel.annotation import AnnotationSourceSpec
from ldsc._kernel.ref_panel import RefPanelSpec
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


class CommonConfigTest(unittest.TestCase):
    def test_defaults(self):
        config = CommonConfig()
        self.assertEqual(config.snp_identifier, "chr_pos")
        self.assertEqual(config.log_level, "INFO")
        self.assertFalse(config.fail_on_missing_metadata)

    def test_validates_values(self):
        with self.assertRaises(ValueError):
            CommonConfig(snp_identifier="bad")
        with self.assertRaises(ValueError):
            CommonConfig(genome_build="hg18")
        with self.assertRaises(ValueError):
            CommonConfig(log_level="trace")

    def test_normalizes_genome_build_aliases(self):
        self.assertEqual(CommonConfig(genome_build="hg37").genome_build, "hg19")
        self.assertEqual(CommonConfig(genome_build="GRCh37").genome_build, "hg19")
        self.assertEqual(CommonConfig(genome_build="GRCh38").genome_build, "hg38")
        self.assertEqual(normalize_genome_build("hg37"), "hg19")
        self.assertEqual(normalize_genome_build("GRCh37"), "hg19")
        self.assertEqual(normalize_genome_build("GRCh38"), "hg38")


class WorkflowConfigTest(unittest.TestCase):
    def test_annotation_config_normalizes_paths(self):
        config = AnnotationBuildConfig(
            baseline_annotation_paths=[Path("base.annot.gz")],
            query_bed_paths=[Path("query.bed")],
            out_prefix=Path("outputs") / "annot",
        )
        self.assertEqual(config.baseline_annotation_paths, ("base.annot.gz",))
        self.assertEqual(config.query_bed_paths, ("query.bed",))
        self.assertEqual(config.out_prefix, "outputs/annot")
        self.assertTrue(config.batch_mode)

    def test_annotation_config_accepts_single_string_tokens_for_plural_fields(self):
        config = AnnotationBuildConfig(
            baseline_annotation_paths="baseline.@",
            query_annotation_paths="query.*.annot.gz",
            query_bed_paths="beds/*.bed",
        )
        self.assertEqual(config.baseline_annotation_paths, ("baseline.@",))
        self.assertEqual(config.query_annotation_paths, ("query.*.annot.gz",))
        self.assertEqual(config.query_bed_paths, ("beds/*.bed",))

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

    def test_munge_config_defaults(self):
        config = MungeConfig(out_prefix="out")
        self.assertEqual(config.info_min, 0.9)
        self.assertEqual(config.maf_min, 0.01)
        self.assertEqual(config.chunk_size, int(5e6))

    def test_munge_config_normalizes_pathlike_fields(self):
        config = MungeConfig(
            out_prefix=Path("results") / "trait",
            merge_alleles_path=Path("resources") / "alleles.tsv",
        )
        self.assertEqual(config.out_prefix, "results/trait")
        self.assertEqual(config.merge_alleles_path, "resources/alleles.tsv")

    def test_regression_config_defaults(self):
        config = RegressionConfig()
        self.assertEqual(config.n_blocks, 200)
        self.assertTrue(config.use_m_5_50)
        with self.assertRaises(ValueError):
            RegressionConfig(n_blocks=1)

    def test_common_config_normalizes_optional_path(self):
        config = CommonConfig(global_snp_restriction_path=Path("restrict") / "snps.txt")
        self.assertEqual(config.global_snp_restriction_path, "restrict/snps.txt")

    def test_ref_panel_config_normalizes_path_fields(self):
        config = RefPanelConfig(
            plink_prefix=Path("plink") / "panel",
            parquet_r2_paths=[Path("r2") / "chr1.parquet"],
            frequency_paths_chr=[Path("freq") / "chr@"],
        )
        self.assertEqual(config.plink_prefix, "plink/panel")
        self.assertEqual(config.parquet_r2_paths, ("r2/chr1.parquet",))
        self.assertEqual(config.frequency_paths_chr, ("freq/chr@",))

    def test_ref_panel_config_accepts_single_string_tokens_for_plural_fields(self):
        config = RefPanelConfig(
            parquet_r2_paths="r2/reference.@",
            frequency_paths="meta/reference.*.tsv.gz",
        )
        self.assertEqual(config.parquet_r2_paths, ("r2/reference.@",))
        self.assertEqual(config.frequency_paths, ("meta/reference.*.tsv.gz",))

    def test_public_specs_normalize_pathlike_inputs(self):
        raw = RawSumstatsSpec(path=Path("sumstats") / "trait.tsv.gz")
        annot = AnnotationSourceSpec(
            baseline_annot_paths=(Path("baseline") / "base.1.annot.gz",),
            query_annot_paths=(Path("query") / "custom.1.annot.gz",),
            bed_paths=(Path("beds") / "enhancer.bed",),
        )
        ref = RefPanelSpec(
            backend="parquet_r2",
            bfile_prefix=Path("plink") / "panel",
            r2_table_paths=(Path("r2") / "chr1.parquet",),
            maf_metadata_paths=(Path("meta") / "chr1.tsv.gz",),
        )
        self.assertEqual(raw.path, "sumstats/trait.tsv.gz")
        self.assertEqual(annot.baseline_annot_paths, ("baseline/base.1.annot.gz",))
        self.assertEqual(annot.query_annot_paths, ("query/custom.1.annot.gz",))
        self.assertEqual(annot.bed_paths, ("beds/enhancer.bed",))
        self.assertEqual(ref.bfile_prefix, "plink/panel")
        self.assertEqual(ref.r2_table_paths, ("r2/chr1.parquet",))
        self.assertEqual(ref.maf_metadata_paths, ("meta/chr1.tsv.gz",))

    def test_ref_panel_spec_normalizes_genome_build_aliases(self):
        self.assertEqual(RefPanelSpec(backend="parquet_r2", genome_build="hg37").genome_build, "hg19")
        self.assertEqual(RefPanelSpec(backend="parquet_r2", genome_build="GRCh37").genome_build, "hg19")
        self.assertEqual(RefPanelSpec(backend="parquet_r2", genome_build="GRCh38").genome_build, "hg38")

    def test_public_specs_accept_single_string_tokens_for_plural_fields(self):
        annot = AnnotationSourceSpec(
            baseline_annot_paths="baseline.@",
            query_annot_paths="query.*.annot.gz",
            bed_paths="beds/*.bed",
        )
        ref = RefPanelSpec(
            backend="parquet_r2",
            r2_table_paths="r2/reference.@",
            maf_metadata_paths="meta/reference.*.tsv.gz",
        )
        self.assertEqual(annot.baseline_annot_paths, ("baseline.@",))
        self.assertEqual(annot.query_annot_paths, ("query.*.annot.gz",))
        self.assertEqual(annot.bed_paths, ("beds/*.bed",))
        self.assertEqual(ref.r2_table_paths, ("r2/reference.@",))
        self.assertEqual(ref.maf_metadata_paths, ("meta/reference.*.tsv.gz",))


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
        self.assertEqual(normalize_chromosome("X"), "X")

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
    def test_read_global_snp_restriction_rsid_plain(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "restrict.txt"
            path.write_text("rs1\nrs2\n", encoding="utf-8")
            self.assertEqual(read_global_snp_restriction(path, "rsid"), {"rs1", "rs2"})

    def test_read_global_snp_restriction_rsid_table(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "restrict.tsv"
            path.write_text("id\tother\nrs1\ta\nrs2\tb\n", encoding="utf-8")
            self.assertEqual(read_global_snp_restriction(path, "rsID"), {"rs1", "rs2"})

    def test_read_global_snp_restriction_chr_pos_table(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "restrict.tsv"
            path.write_text("chrom\tbp\nchr1\t10\n2\t20\n", encoding="utf-8")
            self.assertEqual(read_global_snp_restriction(path, "chr_pos"), {"1:10", "2:20"})
