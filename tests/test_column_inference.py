from pathlib import Path
import sys
import unittest
import warnings

SRC = Path(__file__).resolve().parents[1] / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from ldsc import column_inference as ci
from ldsc.column_inference import (
    PARQUET_R2_CANONICAL_SPECS,
    POS_COLUMN_SPEC,
    R2_HELPER_COLUMN_SPECS,
    R2_SOURCE_COLUMN_SPECS,
    resolve_restriction_chr_pos_columns,
    resolve_restriction_rsid_column,
    resolve_required_column,
    resolve_required_columns,
)


class ColumnInferenceTest(unittest.TestCase):
    def _spec_by_canonical(self, specs, canonical):
        for spec in specs:
            if spec.canonical == canonical:
                return spec
        self.fail(f"Missing ColumnSpec for canonical field {canonical!r}")

    def test_registry_exposes_shared_normalizers(self):
        self.assertIsNotNone(getattr(ci, "clean_header", None))
        self.assertEqual(ci.clean_header("foo-bar.foo_BaR\n"), "FOO_BAR_FOO_BAR")

        self.assertIsNotNone(getattr(ci, "normalize_snp_identifier_mode", None))
        self.assertEqual(
            tuple(ci.SNP_IDENTIFIER_MODES),
            ("rsid", "rsid_allele_aware", "chr_pos", "chr_pos_allele_aware"),
        )
        for value in ("rsid", "rsid_allele_aware", "chr_pos", "chr_pos_allele_aware"):
            self.assertEqual(ci.normalize_snp_identifier_mode(value), value)
        for value in ("rsID", "SNPID", "snp_id", "snp", "chrpos", "rsid_alleles", "chr_pos_alleles"):
            with self.assertRaises(ValueError):
                ci.normalize_snp_identifier_mode(value)

        self.assertIsNotNone(getattr(ci, "normalize_genome_build", None))
        self.assertEqual(ci.normalize_genome_build("hg37"), "hg19")
        self.assertEqual(ci.normalize_genome_build("GRCh37"), "hg19")
        self.assertEqual(ci.normalize_genome_build("GRCh38"), "hg38")

    def test_raw_sumstats_families_cover_legacy_and_current_aliases(self):
        specs = getattr(ci, "RAW_SUMSTATS_REQUIRED_OR_OPTIONAL_SPECS", None)
        self.assertIsNotNone(specs)

        snp_spec = self._spec_by_canonical(specs, "SNP")
        self.assertEqual(resolve_required_column(["rs_number"], snp_spec), "rs_number")
        self.assertEqual(resolve_required_column(["snp_id"], snp_spec), "snp_id")
        self.assertEqual(resolve_required_column(["id"], snp_spec), "id")

        n_spec = self._spec_by_canonical(specs, "N")
        self.assertEqual(resolve_required_column(["weight"], n_spec), "weight")

        n_cas_spec = self._spec_by_canonical(specs, "N_CAS")
        self.assertEqual(resolve_required_column(["ncas"], n_cas_spec), "ncas")
        self.assertEqual(resolve_required_column(["Nca"], n_cas_spec), "Nca")

        n_con_spec = self._spec_by_canonical(specs, "N_CON")
        self.assertEqual(resolve_required_column(["ncon"], n_con_spec), "ncon")
        self.assertEqual(resolve_required_column(["Nco"], n_con_spec), "Nco")

        info_spec = self._spec_by_canonical(specs, "INFO")
        self.assertEqual(resolve_required_column(["IMPINFO"], info_spec), "IMPINFO")

        frq_spec = self._spec_by_canonical(specs, "FRQ")
        self.assertEqual(resolve_required_column(["frq_u"], frq_spec), "frq_u")

        chr_spec = self._spec_by_canonical(specs, "CHR")
        self.assertEqual(resolve_required_column(["#CHROM"], chr_spec), "#CHROM")
        self.assertEqual(resolve_required_column(["CHROM"], chr_spec), "CHROM")

        pos_spec = self._spec_by_canonical(specs, "POS")
        self.assertEqual(resolve_required_column(["BP"], pos_spec), "BP")

        signed_specs = getattr(ci, "RAW_SUMSTATS_SIGNED_STAT_SPECS", None)
        self.assertIsNotNone(signed_specs)
        beta_spec = self._spec_by_canonical(signed_specs, "BETA")
        self.assertEqual(resolve_required_column(["effects"], beta_spec), "effects")

    def test_external_and_internal_families_have_different_strictness(self):
        annotation_specs = getattr(ci, "ANNOTATION_METADATA_SPECS", None)
        self.assertIsNotNone(annotation_specs)
        annotation_pos = self._spec_by_canonical(annotation_specs, "POS")
        self.assertEqual(resolve_required_column(["BP"], annotation_pos), "BP")

        restriction_specs = getattr(ci, "RESTRICTION_RSID_SPECS", None)
        self.assertIsNotNone(restriction_specs)
        restriction_snp = self._spec_by_canonical(restriction_specs, "SNP")
        self.assertEqual(resolve_required_column(["rs_id"], restriction_snp), "rs_id")
        self.assertEqual(resolve_required_column(["id"], restriction_snp), "id")

        internal_specs = getattr(ci, "INTERNAL_LDSCORE_ARTIFACT_SPECS", None)
        self.assertIsNotNone(internal_specs)
        internal_snp = self._spec_by_canonical(internal_specs, "SNP")
        internal_pos = self._spec_by_canonical(internal_specs, "POS")
        with self.assertRaises(ValueError):
            resolve_required_column(["CHR", "rsid", "POS", "CM"], internal_snp)
        with self.assertRaises(ValueError):
            resolve_required_column(["CHR", "SNP", "BP", "CM"], internal_pos)

    def test_resolve_required_column_logs_when_alias_maps_to_canonical_field(self):
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            with self.assertLogs("LDSC.columns", level="INFO") as captured:
                column = resolve_required_column(
                    ["CHR", "BP", "SNP"],
                    POS_COLUMN_SPEC,
                    context="test-position-alias",
                )

        self.assertEqual(column, "BP")
        self.assertEqual(caught, [])
        self.assertTrue(all(record.levelname == "INFO" for record in captured.records))
        self.assertIn("canonical field 'POS'", captured.output[0])
        self.assertIn("input column 'BP'", captured.output[0])

    def test_resolve_r2_source_columns_accepts_legacy_bp_and_build_aliases(self):
        mapping = resolve_required_columns(
            [
                "chromosome",
                "rs_id_1",
                "rs_id_2",
                "grch38_bp_1",
                "grch38_bp_2",
                "hg37_bp1",
                "hg37_bp2",
                "hg38_uniq_id_1",
                "hg38_uniq_id_2",
                "hg19_uniq_id_1",
                "hg19_uniq_id_2",
                "r2",
                "d_prime",
                "corr",
            ],
            R2_SOURCE_COLUMN_SPECS,
            context="test-r2-source-aliases",
        )

        self.assertEqual(mapping["chr"], "chromosome")
        self.assertEqual(mapping["rsID_1"], "rs_id_1")
        self.assertEqual(mapping["rsID_2"], "rs_id_2")
        self.assertEqual(mapping["hg38_pos_1"], "grch38_bp_1")
        self.assertEqual(mapping["hg38_pos_2"], "grch38_bp_2")
        self.assertEqual(mapping["hg19_pos_1"], "hg37_bp1")
        self.assertEqual(mapping["hg19_pos_2"], "hg37_bp2")
        self.assertEqual(mapping["Dprime"], "d_prime")
        self.assertEqual(mapping["+/-corr"], "corr")

    def test_resolve_r2_helper_columns_accepts_bp_aliases(self):
        mapping = resolve_required_columns(
            ["chr", "bp1", "bp_2"],
            R2_HELPER_COLUMN_SPECS,
            context="test-r2-helper-aliases",
        )

        self.assertEqual(mapping["pos_1"], "bp1")
        self.assertEqual(mapping["pos_2"], "bp_2")

    def test_resolve_canonical_parquet_columns_accepts_aliases(self):
        mapping = resolve_required_columns(
            [
                "chr",
                "bp_1",
                "bp2",
                "rsid_1",
                "rs_2",
                "allele1_1",
                "allele2_1",
                "allele1_2",
                "allele2_2",
                "R2",
            ],
            PARQUET_R2_CANONICAL_SPECS,
            context="test-canonical-parquet-aliases",
        )

        self.assertEqual(mapping["CHR"], "chr")
        self.assertEqual(mapping["POS_1"], "bp_1")
        self.assertEqual(mapping["POS_2"], "bp2")
        self.assertEqual(mapping["SNP_1"], "rsid_1")
        self.assertEqual(mapping["SNP_2"], "rs_2")
        self.assertEqual(mapping["A1_1"], "allele1_1")
        self.assertEqual(mapping["A2_1"], "allele2_1")
        self.assertEqual(mapping["A1_2"], "allele1_2")
        self.assertEqual(mapping["A2_2"], "allele2_2")
        self.assertEqual(mapping["R2"], "R2")

    def test_resolve_restriction_rsid_column_uses_registry_aliases(self):
        self.assertEqual(
            resolve_restriction_rsid_column(["variant", "rs_id", "other"], context="test-restrict-rsid"),
            "rs_id",
        )

    def test_resolve_restriction_chr_pos_columns_prefers_generic_pos_aliases(self):
        self.assertEqual(
            resolve_restriction_chr_pos_columns(
                ["chr", "bp", "hg19_pos"],
                genome_build="hg19",
                context="test-restrict-chr-pos",
            ),
            ("chr", "bp"),
        )

    def test_resolve_restriction_chr_pos_columns_accepts_build_specific_position_alias(self):
        self.assertEqual(
            resolve_restriction_chr_pos_columns(
                ["rsID", "chr", "hg19_pos", "other"],
                genome_build="hg19",
                context="test-restrict-hg19-pos",
            ),
            ("chr", "hg19_pos"),
        )
