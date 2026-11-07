"""Collection markers and shared fixtures for the LDSC test suite."""

import unittest

import pytest


def pytest_collection_modifyitems(items):
    """Attach broad markers without forcing every unittest-style test to import pytest."""
    for item in items:
        path = str(item.path)
        if any(name in path for name in ("test_formats.py", "test_plink_io.py", "test_ldscore_workflow.py", "test_ref_panel.py", "test_r2_matrix_builder.py", "test_r2_format.py", "test_windows.py")):
            item.add_marker(pytest.mark.io)
        if any(name in path for name in ("test_formats.py", "test_plink_io.py")):
            item.add_marker(pytest.mark.file_format_compat)


@pytest.fixture(autouse=True)
def _tmp_path_for_unittest_classes(request):
    """Give each unittest.TestCase test pytest's per-test ``tmp_path`` as ``self.tmp_path``."""
    if isinstance(request.instance, unittest.TestCase):
        request.instance.tmp_path = request.getfixturevalue("tmp_path")
