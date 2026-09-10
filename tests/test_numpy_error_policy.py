"""Numerical workflows must not change their caller's NumPy error policy."""

import subprocess
import sys
import unittest


class NumpyErrorPolicyTest(unittest.TestCase):
    def test_public_workflow_imports_preserve_caller_error_policy(self):
        for export in ("RegressionRunner", "SumstatsMunger"):
            with self.subTest(export=export):
                completed = subprocess.run(
                    [sys.executable, "-c", (
                        "import numpy as np\n"
                        "np.seterr(divide='warn', invalid='raise', over='raise', under='ignore')\n"
                        "before = np.geterr().copy()\n"
                        f"from ldsc import {export}\n"
                        "assert np.geterr() == before, (before, np.geterr())\n"
                    )],
                    capture_output=True, text=True, check=False,
                )
                self.assertEqual(completed.returncode, 0, completed.stderr)
