"""The Python API reference page documents exactly the exported names.

``docs/reference/python-api.md`` groups ``ldsc.__all__`` into hand-written
``autosummary`` directives; the documentation build generates one page per name
listed there. These checks keep the grouping in step with the package.
"""

from collections import Counter
from pathlib import Path
import re

import ldsc


API_PAGE = Path(__file__).resolve().parents[2] / "docs" / "reference" / "python-api.md"
_BLOCK = re.compile(r"^\.\. autosummary::\n(.*?)^```", re.MULTILINE | re.DOTALL)


def _listed_names() -> list[str]:
    names = []
    for block in _BLOCK.findall(API_PAGE.read_text(encoding="utf-8")):
        for line in block.splitlines():
            line = line.strip()
            if line and not line.startswith(":"):
                names.append(line)
    return names


def test_api_page_lists_every_exported_name():
    assert sorted(set(ldsc.__all__) - set(_listed_names())) == []


def test_api_page_lists_only_exported_names_once():
    listed = _listed_names()
    assert sorted(set(listed) - set(ldsc.__all__)) == []
    assert [name for name, count in Counter(listed).items() if count > 1] == []
