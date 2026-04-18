from pathlib import Path
import sys
import unittest
import warnings

SRC = Path(__file__).resolve().parents[1] / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from ldsc.column_inference import (
    POS_COLUMN_SPEC,
    R2_HELPER_COLUMN_SPECS,
    R2_SOURCE_COLUMN_SPECS,
    resolve_required_column,
    resolve_required_columns,
)


class ColumnInferenceTest(unittest.TestCase):
    def test_resolve_required_column_warns_when_alias_maps_to_canonical_field(self):
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            column = resolve_required_column(
                ["CHR", "BP", "SNP"],
                POS_COLUMN_SPEC,
                context="test-position-alias",
            )

        self.assertEqual(column, "BP")
        self.assertEqual(len(caught), 1)
        self.assertIn("canonical field 'POS'", str(caught[0].message))
        self.assertIn("input column 'BP'", str(caught[0].message))

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
