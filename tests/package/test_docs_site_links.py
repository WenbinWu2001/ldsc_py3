"""Messages link to pages and sections of the published documentation site.

Users of an installed package have no ``docs/`` folder, so package text links
to the site instead. Each link must name a built page and one of its headings.
On the site a heading's anchor is the docutils id of its text (ASCII,
lower case, runs of other characters become one hyphen), which differs from
GitHub's anchor when the heading contains ``/`` or ``R2`` with a superscript.
"""

from pathlib import Path
import re
import unicodedata

import pytest


REPO = Path(__file__).resolve().parents[2]
SITE = "https://wenbinwu2001.github.io/ldsc_py3/"
_SITE_LINK = re.compile(re.escape(SITE) + r"([A-Za-z0-9_/-]*?\.html)(?:#([A-Za-z0-9-]+))?")


def _site_anchor(heading: str) -> str:
    """Return the id docutils gives a section with this Markdown heading text."""
    text = re.sub(r"[`*]", "", heading)
    text = unicodedata.normalize("NFKD", text.lower()).encode("ascii", "ignore").decode("ascii")
    text = re.sub(r"[^a-z0-9]+", "-", " ".join(text.split()))
    return re.sub(r"^[-0-9]+|-+$", "", text)


def _source_texts():
    for path in sorted((REPO / "src" / "ldsc").rglob("*.py")):
        yield path, path.read_text(encoding="utf-8")


def _site_links():
    return sorted({(path.relative_to(REPO).as_posix(), *match.groups())
                   for path, text in _source_texts() for match in _SITE_LINK.finditer(text)})


def test_messages_no_longer_point_at_the_troubleshooting_markdown_file():
    stale = [path.relative_to(REPO).as_posix() for path, text in _source_texts()
             if "docs/user-guide/troubleshooting.md" in text]
    assert stale == []


@pytest.mark.parametrize(("source", "page", "anchor"), _site_links())
def test_site_link_names_an_existing_page_and_heading(source, page, anchor):
    markdown = REPO / "docs" / Path(page).with_suffix(".md")
    assert markdown.is_file(), f"{source} links to a missing page {page}"
    if anchor:
        headings = [line.lstrip("#").strip() for line in markdown.read_text(encoding="utf-8").splitlines()
                    if re.match(r"#{1,6} ", line)]
        assert anchor in {_site_anchor(heading) for heading in headings}, f"{source}: no heading for #{anchor}"
