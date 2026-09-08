from pathlib import Path
from unittest import mock

import runpy


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]


def test_matplotlib_is_a_required_distribution_dependency():
    with mock.patch("setuptools.setup") as setup:
        runpy.run_path(str(REPOSITORY_ROOT / "setup.py"), run_name="__main__")

    metadata = setup.call_args.kwargs
    assert "matplotlib>=3.9,<4" in metadata["install_requires"]
    assert "plot" not in metadata["extras_require"]
