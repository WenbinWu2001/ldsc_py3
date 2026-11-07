"""Sphinx helpers for the LDSC3 documentation site.

The Markdown pages under ``docs/`` are written to read correctly on GitHub.
This extension makes the same files build into the site without edits:

- ``github_slug`` gives headings the anchors GitHub gives them, so links such
  as ``cli.md#ldsc-rg`` resolve in both places.
- ``cli_<command>`` functions expose each subcommand parser to the
  ``argparse`` directive in ``reference/cli.md``.
- A ``source-read`` hook points links to files that are not site pages
  (package source, tests, ``AGENTS.md``, excluded pages) at GitHub instead of
  copying those files into the site as downloads.
"""

from __future__ import annotations

import functools
import importlib
import inspect
import re
import unicodedata
from pathlib import Path

from ldsc.cli import _COMMANDS, _command_parser

DOCS_ROOT = Path(__file__).resolve().parents[1]
REPO_ROOT = DOCS_ROOT.parent
GITHUB_REPO = "https://github.com/WenbinWu2001/ldsc_py3"


# Source links follow the published branch; line anchors may drift from the
# text until the site is rebuilt.
SOURCE_REF = "restructure"


def github_slug(title: str) -> str:
    """Return the anchor GitHub generates for a Markdown heading.

    GitHub lower-cases the text, keeps letters, combining marks, decimal
    digits, spaces, hyphens and underscores, drops everything else (so
    ``R2`` written with a superscript loses the superscript), and turns
    spaces into hyphens.
    """
    kept = (
        char
        for char in title.strip().lower()
        if char in " -_" or unicodedata.category(char)[0] in "LM" or unicodedata.category(char) == "Nd"
    )
    return "".join(kept).replace(" ", "-")


# One zero-argument parser factory per subcommand, for example ``cli_munge_sumstats``.
for _name in _COMMANDS:
    globals()[f"cli_{_name.replace('-', '_')}"] = functools.partial(_command_parser, _name)


def linkcode_resolve(domain: str, info: dict) -> str | None:
    """Return the GitHub line range that defines a documented API object."""
    if domain != "py" or not info.get("module"):
        return None
    obj = importlib.import_module(info["module"])
    for part in info["fullname"].split("."):
        obj = getattr(obj, part, None)
        if obj is None:
            return None
    obj = inspect.unwrap(obj.fget if isinstance(obj, property) else obj)
    try:
        path = Path(inspect.getsourcefile(obj)).resolve()
        lines, start = inspect.getsourcelines(obj)
    except (TypeError, OSError):
        # Data attributes and builtins have no source lines to link.
        return None
    if not path.is_relative_to(REPO_ROOT):
        return None
    end = start + len(lines) - 1
    return f"{GITHUB_REPO}/blob/{SOURCE_REF}/{path.relative_to(REPO_ROOT).as_posix()}#L{start}-L{end}"


def offsite_url(resolved: Path, anchor: str, site_docs: set[str]) -> str | None:
    """Return the GitHub URL for a link target that is not a page of the site.

    Parameters
    ----------
    resolved : pathlib.Path
        Absolute, normalized path the relative link points to, for example
        ``REPO_ROOT / "src/ldsc/_regression/estimators.py"``.
    anchor : str
        Fragment after ``#`` without the ``#``, such as ``"L410"`` or a
        heading anchor; empty when the link has none.
    site_docs : set of str
        Names of the documents built into the site, relative to ``docs/``
        and without suffix, such as ``"reference/cli"``.

    Returns
    -------
    str or None
        The GitHub URL that replaces the link target, or ``None`` to leave the
        link to Sphinx: site pages, and targets outside the repository, which
        Sphinx then reports as broken.
    """
    if not resolved.is_relative_to(REPO_ROOT):
        return None
    if resolved.suffix == ".md" and resolved.is_relative_to(DOCS_ROOT):
        if resolved.relative_to(DOCS_ROOT).with_suffix("").as_posix() in site_docs:
            return None
    # Every other file goes to GitHub, so the site never copies it in as a download.
    kind = "tree" if resolved.is_dir() else "blob"
    url = f"{GITHUB_REPO}/{kind}/{SOURCE_REF}/{resolved.relative_to(REPO_ROOT).as_posix()}"
    return f"{url}#{anchor}" if anchor else url


_FENCE = re.compile(r"^\s*(```|~~~)")
# Inline links, not images: [text](target) or [text](<target>).
_LINK = re.compile(r"(?<!!)(\[[^\]\n]*\]\()(\s*<?)([^)\s>]+)(>?\s*\))")


def _rewrite_links(app, docname: str, source: list[str]) -> None:
    """Replace relative link targets that ``offsite_url`` maps to GitHub."""
    doc_dir = (DOCS_ROOT / docname).parent
    site_docs = set(app.env.found_docs)

    def replace(match: re.Match) -> str:
        target = match.group(3)
        if "://" in target or target.startswith(("#", "mailto:")):
            return match.group(0)
        path, _, anchor = target.partition("#")
        url = offsite_url(Path(doc_dir, path).resolve(), anchor, site_docs)
        if url is None:
            return match.group(0)
        return f"{match.group(1)}{match.group(2)}{url}{match.group(4)}"

    lines = source[0].split("\n")
    in_fence = False
    for index, line in enumerate(lines):
        if _FENCE.match(line):
            in_fence = not in_fence
        elif not in_fence:
            lines[index] = _LINK.sub(replace, line)
    source[0] = "\n".join(lines)


def setup(app):
    from pygments.lexers.special import TextLexer

    # Fences tagged ``tsv`` or ``bed`` name the file format; show them as plain text.
    app.add_lexer("tsv", TextLexer)
    app.add_lexer("bed", TextLexer)
    app.connect("source-read", _rewrite_links)
    return {"parallel_read_safe": True, "parallel_write_safe": True}
