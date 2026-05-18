from pathlib import Path
import sys
import unittest

import pandas as pd

SRC = Path(__file__).resolve().parents[1] / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from ldsc._kernel import snp_identity as si


class SNPIdentityModeTest(unittest.TestCase):
    def test_only_exact_public_modes_are_accepted(self):
        self.assertEqual(
            tuple(si.SNP_IDENTIFIER_MODES),
            ("rsid", "rsid_allele_aware", "chr_pos", "chr_pos_allele_aware"),
        )
        for mode in si.SNP_IDENTIFIER_MODES:
            self.assertEqual(si.normalize_snp_identifier_mode(mode), mode)

        rejected = ["rsID", "SNPID", "snp_id", "chrpos", "rsid_alleles", "chr_pos_alleles"]
        for value in rejected:
            with self.subTest(value=value):
                with self.assertRaises(ValueError):
                    si.normalize_snp_identifier_mode(value)

    def test_mode_family_helpers_return_base_family(self):
        self.assertEqual(si.identity_mode_family("rsid"), "rsid")
        self.assertEqual(si.identity_mode_family("rsid_allele_aware"), "rsid")
        self.assertEqual(si.identity_mode_family("chr_pos"), "chr_pos")
        self.assertEqual(si.identity_mode_family("chr_pos_allele_aware"), "chr_pos")

        self.assertEqual(si.identity_base_mode("rsid_allele_aware"), "rsid")
        self.assertEqual(si.identity_base_mode("chr_pos_allele_aware"), "chr_pos")
        self.assertFalse(si.is_allele_aware_mode("rsid"))
        self.assertTrue(si.is_allele_aware_mode("chr_pos_allele_aware"))


class AlleleNormalizationTest(unittest.TestCase):
    def test_unordered_strand_aware_sets_normalize_to_canonical_pair(self):
        for left, right in [("A", "C"), ("C", "A"), ("T", "G"), ("G", "T")]:
            with self.subTest(pair=(left, right)):
                self.assertEqual(si.normalize_allele_set(left, right), "A:C")

    def test_invalid_pairs_return_required_failure_reasons(self):
        examples = [
            (None, "A", "missing_allele"),
            ("", "A", "missing_allele"),
            ("A", "A", "invalid_allele"),
            ("AC", "T", "invalid_allele"),
            ("A", "N", "invalid_allele"),
            ("A", "T", "strand_ambiguous_allele"),
            ("C", "G", "strand_ambiguous_allele"),
        ]
        for a1, a2, reason in examples:
            with self.subTest(pair=(a1, a2)):
                with self.assertRaisesRegex(ValueError, reason):
                    si.normalize_allele_set(a1, a2)
                values, reasons = si.allele_set_series(pd.DataFrame({"A1": [a1], "A2": [a2]}), context="test")
                self.assertTrue(pd.isna(values.iloc[0]))
                self.assertEqual(reasons.iloc[0], reason)


class EffectiveKeyTest(unittest.TestCase):
    def test_effective_keys_for_all_modes(self):
        frame = pd.DataFrame(
            {
                "CHR": [1],
                "POS": [101],
                "SNP": ["rs1"],
                "A1": ["G"],
                "A2": ["T"],
            }
        )

        self.assertEqual(si.effective_merge_key_series(frame, "rsid").tolist(), ["rs1"])
        self.assertEqual(si.effective_merge_key_series(frame, "rsid_allele_aware").tolist(), ["rs1:A:C"])
        self.assertEqual(si.effective_merge_key_series(frame, "chr_pos").tolist(), ["1:101"])
        self.assertEqual(si.effective_merge_key_series(frame, "chr_pos_allele_aware").tolist(), ["1:101:A:C"])

    def test_effective_key_raises_for_invalid_alleles_in_allele_aware_modes(self):
        frame = pd.DataFrame({"CHR": [1], "POS": [101], "SNP": ["rs1"], "A1": ["A"], "A2": ["T"]})

        with self.assertRaisesRegex(ValueError, "strand_ambiguous_allele"):
            si.effective_merge_key_series(frame, "rsid_allele_aware")


class IdentityArtifactCleanupTest(unittest.TestCase):
    def test_base_modes_drop_duplicate_clusters_without_inspecting_alleles(self):
        frame = pd.DataFrame(
            {
                "CHR": [1, 1, 1],
                "POS": [10, 10, 11],
                "SNP": ["rs1", "rs1", "rs2"],
                "A1": ["A", "bad", "A"],
                "A2": ["T", None, None],
            }
        )

        result = si.clean_identity_artifact_table(frame, "rsid", context="test", stage="artifact")

        self.assertEqual(result.cleaned["SNP"].tolist(), ["rs2"])
        self.assertEqual(result.dropped["reason"].tolist(), ["duplicate_identity", "duplicate_identity"])
        self.assertEqual(set(result.dropped["SNP"]), {"rs1"})
        self.assertEqual(result.dropped["source_pos"].tolist(), [10, 10])

    def test_base_modes_keep_singletons_with_invalid_missing_or_ambiguous_alleles(self):
        frame = pd.DataFrame(
            {
                "CHR": [1, 1, 1],
                "POS": [10, 20, 30],
                "SNP": ["rs1", "rs2", "rs3"],
                "A1": [None, "AC", "A"],
                "A2": ["A", "T", "T"],
            }
        )

        result = si.clean_identity_artifact_table(frame, "chr_pos", context="test", stage="artifact")

        self.assertEqual(result.cleaned["SNP"].tolist(), ["rs1", "rs2", "rs3"])
        self.assertTrue(result.dropped.empty)

    def test_allele_aware_mode_drops_multi_allelic_base_key_clusters(self):
        frame = pd.DataFrame(
            {
                "CHR": [1, 1, 1],
                "POS": [10, 10, 20],
                "SNP": ["rs1", "rs1", "rs2"],
                "A1": ["A", "A", "A"],
                "A2": ["C", "G", "C"],
            }
        )

        result = si.clean_identity_artifact_table(frame, "rsid_allele_aware", context="test", stage="artifact")

        self.assertEqual(result.cleaned["SNP"].tolist(), ["rs2"])
        self.assertEqual(result.dropped["reason"].tolist(), ["multi_allelic_base_key", "multi_allelic_base_key"])

    def test_allele_aware_mode_drops_duplicate_identity_clusters(self):
        frame = pd.DataFrame(
            {
                "CHR": [1, 1, 1],
                "POS": [10, 10, 20],
                "SNP": ["rs1", "rs1", "rs2"],
                "A1": ["A", "C", "A"],
                "A2": ["C", "A", "C"],
            }
        )

        result = si.clean_identity_artifact_table(frame, "rsid_allele_aware", context="test", stage="artifact")

        self.assertEqual(result.cleaned["SNP"].tolist(), ["rs2"])
        self.assertEqual(result.dropped["reason"].tolist(), ["duplicate_identity", "duplicate_identity"])


class RestrictionIdentityKeysTest(unittest.TestCase):
    def test_restriction_duplicate_keys_collapse_instead_of_drop_all(self):
        frame = pd.DataFrame(
            {
                "CHR": [1, 1, 1],
                "POS": [10, 10, 10],
                "SNP": ["rs1", "rs1", "rs1"],
                "A1": ["A", "C", "T"],
                "A2": ["C", "A", "G"],
            }
        )

        result = si.collapse_restriction_identity_keys(frame, "rsid_allele_aware", context="restriction")

        self.assertEqual(result.keys, {"rs1:A:C"})
        self.assertEqual(result.match_kind, "identity")
        self.assertEqual(result.n_input_rows, 3)
        self.assertEqual(result.n_retained_keys, 1)
        self.assertTrue(result.dropped.empty)

    def test_restrictions_without_alleles_use_base_keys_without_allele_inspection(self):
        frame = pd.DataFrame({"CHR": [1, 1], "POS": [10, 10], "SNP": ["rs1", "rs1"]})

        result = si.collapse_restriction_identity_keys(frame, "rsid_allele_aware", context="restriction")

        self.assertEqual(result.keys, {"rs1"})
        self.assertEqual(result.match_kind, "base")
        self.assertEqual(result.n_retained_keys, 1)


if __name__ == "__main__":
    unittest.main()
