import unittest

import numpy as np
import pandas as pd

from ldsc.annotation_semantics import (
    annotation_value_fingerprint,
    classify_annotation_values,
    common_universe_fingerprint,
    require_unique_annotation_names,
)
from ldsc.errors import LDSCInputError


class AnnotationSemanticsTest(unittest.TestCase):
    def test_classification_is_exact_and_advisory(self):
        frame = pd.DataFrame(
            {
                "binary": np.array([0.0, 1.0, 0.0], dtype=np.float32),
                "quantitative": np.array([0.0, 0.5, 1.0], dtype=np.float32),
            }
        )

        self.assertEqual(
            classify_annotation_values(frame),
            {"binary": "binary", "quantitative": "quantitative"},
        )

    def test_duplicate_names_across_groups_are_rejected(self):
        with self.assertRaisesRegex(LDSCInputError, "globally unique.*shared"):
            require_unique_annotation_names(["base", "shared"], ["shared"])

    def test_fingerprints_use_the_versioned_canonical_bytes(self):
        identities = pd.Series(["rs2", "rs1"])
        values = pd.Series([1.0, 0.0])

        self.assertEqual(
            common_universe_fingerprint(identities),
            "2ee29ca633a60847d001f6f7c749ee5c210f02855c5eb6c10329edaa248e2029",
        )
        self.assertEqual(
            annotation_value_fingerprint(identities, values),
            "4eec2756edc870952aad96163a96286df0c6ada0b439a94a85b6048e21b9c28d",
        )


if __name__ == "__main__":
    unittest.main()
