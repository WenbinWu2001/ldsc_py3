"""Pytest fixtures and collection markers for the LDSC test suite."""

from __future__ import annotations

from pathlib import Path

import pytest


TESTS_ROOT = Path(__file__).resolve().parent
FIXTURES_ROOT = TESTS_ROOT / "fixtures"
LEGACY_FIXTURES_ROOT = FIXTURES_ROOT / "legacy"
MINIMAL_EXTERNAL_RESOURCES = FIXTURES_ROOT / "minimal_external_resources"


@pytest.fixture(scope="session")
def fixtures_root() -> Path:
    return FIXTURES_ROOT


@pytest.fixture(scope="session")
def legacy_fixtures_root() -> Path:
    return LEGACY_FIXTURES_ROOT


@pytest.fixture(scope="session")
def minimal_external_resources() -> Path:
    return MINIMAL_EXTERNAL_RESOURCES


def pytest_collection_modifyitems(items):
    """Attach broad markers without forcing every unittest-style test to import pytest."""
    for item in items:
        path = str(item.path)
        if any(name in path for name in ("test_formats.py", "test_plink_io.py", "test_ldscore_workflow.py", "test_ref_panel.py", "test_ref_panel_builder.py")):
            item.add_marker(pytest.mark.io)
        if any(name in path for name in ("test_formats.py", "test_plink_io.py")):
            item.add_marker(pytest.mark.legacy_compat)
