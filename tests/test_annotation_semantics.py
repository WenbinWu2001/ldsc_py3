import unittest

import numpy as np
import pandas as pd

from ldsc.annotation_semantics import (
    classify_annotation_values,
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
        original = frame.copy(deep=True)

        self.assertEqual(
            classify_annotation_values(frame),
            {"binary": "binary", "quantitative": "quantitative"},
        )
        pd.testing.assert_frame_equal(frame, original)

    def test_duplicate_names_across_groups_are_rejected(self):
        with self.assertRaisesRegex(LDSCInputError, "globally unique.*shared"):
            require_unique_annotation_names(["base", "shared"], ["shared"])

if __name__ == "__main__":
    unittest.main()
