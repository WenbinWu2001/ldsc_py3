import logging
import unittest

import numpy as np
import pandas as pd


from ldsc._core import snp_identity as si
from ldsc.errors import LDSCConfigError, LDSCInputError, LDSCUsageError


class SNPIdentityModeTest(unittest.TestCase):
    def test_only_exact_public_modes_are_accepted(self):
        self.assertEqual(
            tuple(si.SNP_IDENTIFIER_MODES),
            ("rsid", "rsid_allele_aware", "chr_pos", "chr_pos_allele_aware"),
        )
        for mode in si.SNP_IDENTIFIER_MODES:
            self.assertEqual(si.normalize_snp_identifier_mode(mode), mode)

        rejected = ["rsID", "SNPID", "snp_id", "snp", "chrpos", "rsid_alleles", "chr_pos_alleles"]
        for value in rejected:
            with self.subTest(value=value):
                with self.assertRaises(LDSCConfigError):
                    si.normalize_snp_identifier_mode(value)

    def test_invalid_mode_message_reports_cause_and_fix(self):
        with self.assertRaises(LDSCConfigError) as ctx:
            si.normalize_snp_identifier_mode("rsID")

        message = str(ctx.exception)
        self.assertIn("snp_identifier mode 'rsID' is invalid", message)
        self.assertIn("Most likely", message)
        self.assertIn("Use one of", message)

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

    def test_a_t_and_c_g_pairs_are_valid_identities(self):
        # A base paired with its own complement folds to itself, so these two
        # classes complete the four canonical allele sets A:C, A:G, A:T and C:G.
        for left, right, token in [("A", "T", "A:T"), ("T", "A", "A:T"), ("C", "G", "C:G"), ("G", "C", "C:G")]:
            with self.subTest(pair=(left, right)):
                self.assertEqual(si.normalize_allele_set(left, right), token)
                values, reasons = si.allele_set_series(pd.DataFrame({"A1": [left], "A2": [right]}), context="test")
                self.assertEqual(values.iloc[0], token)
                self.assertTrue(pd.isna(reasons.iloc[0]))

    def test_invalid_pairs_return_required_failure_reasons(self):
        examples = [
            (None, "A", "missing_allele"),
            ("", "A", "missing_allele"),
            ("A", "A", "invalid_allele"),
            ("AC", "T", "invalid_allele"),
            ("A", "N", "invalid_allele"),
        ]
        for a1, a2, reason in examples:
            with self.subTest(pair=(a1, a2)):
                with self.assertRaisesRegex(ValueError, reason):
                    si.normalize_allele_set(a1, a2)
                values, reasons = si.allele_set_series(pd.DataFrame({"A1": [a1], "A2": [a2]}), context="test")
                self.assertTrue(pd.isna(values.iloc[0]))
                self.assertEqual(reasons.iloc[0], reason)


class EffectiveKeyTest(unittest.TestCase):
    def test_effective_keys_for_all_modes_use_the_expected_identity_components(self):
        frame = pd.DataFrame(
            {
                "CHR": [1, 1, 2, 1, 1],
                "POS": [101, 101, 101, 202, 101],
                "SNP": ["rs1", "rs2", "rs1", "rs1", "rs1"],
                "A1": ["G", "G", "G", "G", "A"],
                "A2": ["T", "T", "T", "T", "G"],
            }
        )

        expected = {
            "rsid": ["rs1", "rs2", "rs1", "rs1", "rs1"],
            "rsid_allele_aware": ["rs1:A:C", "rs2:A:C", "rs1:A:C", "rs1:A:C", "rs1:A:G"],
            "chr_pos": ["1:101", "1:101", "2:101", "1:202", "1:101"],
            "chr_pos_allele_aware": [
                "1:101:A:C",
                "1:101:A:C",
                "2:101:A:C",
                "1:202:A:C",
                "1:101:A:G",
            ],
        }

        for mode, keys in expected.items():
            with self.subTest(mode=mode):
                self.assertEqual(si.identity_key_series(frame, mode).tolist(), keys)

    def test_effective_key_raises_for_invalid_alleles_in_allele_aware_modes(self):
        frame = pd.DataFrame({"CHR": [1], "POS": [101], "SNP": ["rs1"], "A1": ["A"], "A2": ["N"]})

        with self.assertRaisesRegex(LDSCInputError, "invalid_allele"):
            si.identity_key_series(frame, "rsid_allele_aware")


# Canonical cell values and the base identity each one gives, written out by
# hand: a normalized chromosome or position, or None for a missing or invalid
# value, whose row gets no key.
_CHR_CELLS = {1: "1", "1": "1", "chr2": "2", "22": "22", "X": "X", "chrX": "X", "chrUn": None, None: None}
_POS_CELLS = {100: 100, "250": 250, 3000.0: 3000, 0: None, -5: None, "abc": None, 1.5: None, None: None}
_SNP_CELLS = {"rs1": "rs1", "rs2": "rs2", "rs3": "rs3", None: None}
_ALLELE_PAIRS = [(a, b) for a in "ACGT" for b in "ACGT" if a != b]
_COMPLEMENT = str.maketrans("ACGT", "TGCA")


def _expected_key(chrom, pos, snp, a1, a2, mode):
    """Row-by-row identity key of the four modes, independent of the package code."""
    if mode.startswith("rsid"):
        base = _SNP_CELLS[snp]
    else:
        chrom, pos = _CHR_CELLS[chrom], _POS_CELLS[pos]
        base = None if chrom is None or pos is None else f"{chrom}:{pos}"
    if base is None or not mode.endswith("allele_aware"):
        return base
    observed = ":".join(sorted((a1, a2)))
    complement = ":".join(sorted((a1.translate(_COMPLEMENT), a2.translate(_COMPLEMENT))))
    return f"{base}:{min(observed, complement)}"


class IdentityKeySeriesTest(unittest.TestCase):
    def test_keys_of_random_canonical_frames_match_the_row_rule_in_all_modes(self):
        # Covers missing and invalid CHR/POS, missing rsIDs and every valid allele pair.
        rng = np.random.default_rng(20261106)
        chr_cells, pos_cells, snp_cells = list(_CHR_CELLS), list(_POS_CELLS), list(_SNP_CELLS)
        for trial in range(20):
            n = 150
            rows = [(chr_cells[rng.integers(len(chr_cells))], pos_cells[rng.integers(len(pos_cells))],
                     snp_cells[rng.integers(len(snp_cells))], *_ALLELE_PAIRS[rng.integers(len(_ALLELE_PAIRS))])
                    for _ in range(n)]
            frame = pd.DataFrame(rows, columns=["CHR", "POS", "SNP", "A1", "A2"], index=rng.permutation(n) + 7)
            for mode in si.SNP_IDENTIFIER_MODES:
                with self.subTest(trial=trial, mode=mode):
                    keys = si.identity_key_series(frame, mode, context="test")
                    expected = [_expected_key(*row, mode) for row in rows]
                    self.assertEqual(keys.index.tolist(), frame.index.tolist())
                    self.assertEqual([None if pd.isna(key) else key for key in keys], expected)

    def test_logger_reports_rows_without_a_coordinate_key(self):
        frame = pd.DataFrame({"CHR": ["1", "chrUn", "2"], "POS": [100, 200, None], "SNP": ["rs1", "rs2", "rs3"]})
        logger = logging.getLogger("LDSC.test_identity_key_series")
        with self.assertLogs(logger, level="WARNING") as logs:
            keys = si.identity_key_series(frame, "chr_pos", context="LD-score SNP matching", logger=logger)
        self.assertEqual(keys.tolist()[0], "1:100")
        self.assertTrue(keys.iloc[1:].isna().all())
        self.assertEqual(len(logs.output), 1)
        self.assertIn("Dropped 2 SNPs with invalid or missing CHR/POS in LD-score SNP matching", logs.output[0])

    def test_aliased_headers_are_rejected_in_every_mode(self):
        # Callers pass canonical tables; an alias such as BP or rsid names no canonical column,
        # so the key builder stops instead of guessing the column.
        frame = pd.DataFrame({"CHR": ["1", "2"], "BP": [10, None], "rsid": ["rs1", "rs2"],
                              "A1": ["A", "C"], "A2": ["G", "T"]})
        for mode in si.SNP_IDENTIFIER_MODES:
            with self.subTest(mode=mode):
                with self.assertRaisesRegex(LDSCInputError, "missing canonical column"):
                    si.identity_key_series(frame, mode, context="test")


class IdentityMetadataValidationTest(unittest.TestCase):
    def test_identity_artifact_metadata_has_no_schema_version(self):
        meta = si.identity_artifact_metadata(
            artifact_type="sumstats",
            snp_identifier="chr_pos_allele_aware",
            genome_build="hg19",
        )
        self.assertNotIn("schema_version", meta)
        self.assertEqual(
            meta,
            {
                "artifact_type": "sumstats",
                "snp_identifier": "chr_pos_allele_aware",
                "genome_build": "hg19",
            },
        )

    def test_validate_accepts_payload_without_schema_version(self):
        meta = {"artifact_type": "sumstats", "snp_identifier": "rsid_allele_aware", "genome_build": None}
        self.assertEqual(
            si.validate_identity_artifact_metadata(meta, expected_artifact_type="sumstats"),
            "rsid_allele_aware",
        )

    def test_validate_rejects_wrong_artifact_type_without_schema_version(self):
        meta = {"artifact_type": "ldscore", "snp_identifier": "rsid", "genome_build": None}
        with self.assertRaises(LDSCInputError):
            si.validate_identity_artifact_metadata(meta, expected_artifact_type="sumstats")

    def test_invalid_artifact_metadata_reports_cause_fix_and_reference(self):
        metadata = {"schema_version": 0, "artifact_type": "ldscore", "snp_identifier": "chr_pos"}

        with self.assertRaises(LDSCInputError) as ctx:
            si.validate_identity_artifact_metadata(metadata, expected_artifact_type="sumstats")

        message = str(ctx.exception)
        self.assertIn("Could not read LDSC3 sumstats artifact metadata", message)
        self.assertIn("Most likely", message)
        self.assertIn("Regenerate the artifact with the current LDSC3 package", message)
        self.assertIn(
            "https://wenbinwu2001.github.io/ldsc_py3/user-guide/troubleshooting.html#common-ldsc3-artifact-schema-or-provenance-is-incompatible",
            message,
        )

    def test_regression_identity_mismatch_reports_cause_and_fix(self):
        with self.assertRaises(LDSCUsageError) as ctx:
            si.resolve_regression_identity_mode(
                "rsid",
                "chr_pos",
                allow_identity_downgrade=False,
            )

        message = str(ctx.exception)
        self.assertIn("Cannot resolve regression SNP identity", message)
        self.assertIn("Most likely", message)
        self.assertIn("Use inputs from the same snp_identifier family", message)


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

    def test_base_modes_keep_singletons_with_invalid_or_missing_alleles(self):
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


    def test_cleanup_keys_equal_the_keys_of_the_cleaned_table(self):
        frame = pd.DataFrame(
            {
                "CHR": [1, 1, 1, 2, None, 3],
                "POS": [10, 10, 12, 5, 7, 9],
                "SNP": ["rs1", "rs1", "rs2", "rs3", "rs4", "rs5"],
                "A1": ["A", "A", "C", "A", "G", "C"],
                "A2": ["G", "G", "T", "T", "C", "A"],
            }
        )
        for mode in si.SNP_IDENTIFIER_MODES:
            with self.subTest(mode=mode):
                cleaned = si.clean_identity_artifact_table(frame, mode, context="test", stage="artifact")
                keys, dropped = si.identity_cleanup_keys(frame, mode, context="test", stage="artifact")
                expected = si.identity_key_series(cleaned.cleaned, mode, context="test")
                pd.testing.assert_series_equal(pd.Series(keys, dtype=object),
                                               expected.astype(object).rename(None), check_names=False)
                pd.testing.assert_frame_equal(dropped, cleaned.dropped)


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


class DroppedSNPSchemaTest(unittest.TestCase):
    DTYPES = {"CHR": "string", "SNP": "string", "source_pos": "Int64", "target_pos": "Int64", "reason": "string",
              "base_key": "string", "identity_key": "string", "allele_set": "string", "stage": "string"}

    def test_rows_without_identity_columns_get_the_full_schema_and_default_stage(self):
        rows = pd.DataFrame({"CHR": [1, 2], "SNP": ["rs1", "rs2"], "source_pos": [10, 20],
                             "target_pos": [pd.NA, 25], "reason": "unmapped_liftover"}, index=[5, 7])

        frame = si.coerce_identity_drop_frame(rows, default_stage="liftover")

        self.assertEqual(frame.columns.tolist(), si.IDENTITY_DROP_COLUMNS)
        self.assertEqual({column: str(dtype) for column, dtype in frame.dtypes.items()}, self.DTYPES)
        self.assertEqual(frame.index.tolist(), [0, 1])
        self.assertEqual(frame["CHR"].tolist(), ["1", "2"])
        self.assertEqual(frame["target_pos"].tolist(), [pd.NA, 25])
        self.assertEqual(frame["stage"].tolist(), ["liftover", "liftover"])
        self.assertTrue(frame["base_key"].isna().all())

    def test_default_stage_keeps_existing_stages(self):
        rows = pd.DataFrame({"reason": ["duplicate_identity", "unmapped_liftover"],
                             "stage": ["post_liftover_identity_cleanup", pd.NA]})

        frame = si.coerce_identity_drop_frame(rows, default_stage="liftover")

        self.assertEqual(frame["stage"].tolist(), ["post_liftover_identity_cleanup", "liftover"])

    def test_a_non_numeric_position_is_kept(self):
        frame = si.coerce_identity_drop_frame(pd.DataFrame({"source_pos": ["10", "bad"]}))

        self.assertEqual(frame["source_pos"].tolist(), ["10", "bad"])

    def test_empty_and_concatenated_frames_share_the_schema(self):
        empty = si.empty_identity_drop_frame()
        rows = pd.DataFrame({"SNP": ["rs1"], "source_pos": [10], "reason": "source_duplicate"})

        combined = si.concat_identity_drop_frames([empty, None, rows, rows], default_stage="liftover")

        self.assertEqual({column: str(dtype) for column, dtype in empty.dtypes.items()}, self.DTYPES)
        self.assertEqual(len(empty), 0)
        self.assertEqual({column: str(dtype) for column, dtype in combined.dtypes.items()}, self.DTYPES)
        self.assertEqual(combined.index.tolist(), [0, 1])
        self.assertEqual(combined["stage"].tolist(), ["liftover", "liftover"])


if __name__ == "__main__":
    unittest.main()
