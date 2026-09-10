"""Collection markers for the LDSC test suite."""

import pytest


def pytest_collection_modifyitems(items):
    """Attach broad markers without forcing every unittest-style test to import pytest."""
    for item in items:
        path = str(item.path)
        if any(name in path for name in ("test_formats.py", "test_plink_io.py", "test_ldscore_workflow.py", "test_ref_panel.py", "test_ref_panel_builder.py")):
            item.add_marker(pytest.mark.io)
        if any(name in path for name in ("test_formats.py", "test_plink_io.py")):
            item.add_marker(pytest.mark.file_format_compat)
