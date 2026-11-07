"""Sphinx configuration for the LDSC3 documentation site.

Build from the repository root with the ``docs`` extra installed::

    sphinx-build -W --keep-going -b html docs docs/_build/html
"""

import importlib.metadata
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent / "_ext"))

from ldsc_docs import linkcode_resolve  # noqa: E402,F401  (read by sphinx.ext.linkcode)

project = "LDSC3"
author = "Wenbin Wu, Anthony Abrantes, Brendan Bulik-Sullivan, and Hilary Finucane"
release = importlib.metadata.version("ldsc")

extensions = [
    "myst_parser",
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.intersphinx",
    "sphinx.ext.linkcode",
    "sphinx.ext.mathjax",
    "numpydoc",
    "sphinxarg.ext",
    "sphinxcontrib.mermaid",
    "sphinx_copybutton",
    "ldsc_docs",
]

source_suffix = {".md": "markdown", ".rst": "restructuredtext"}
root_doc = "index"
# Contributor-only pages stay in the repository; links to them point to GitHub.
exclude_patterns = [
    "_build",
    "**/.DS_Store",
    "wiki",
    "dev/DOCS_VOCABULARY.md",
    "dev/guidelines",
    "dev/validation/legacy-issues-audit",
]

# MyST reads the pages as GitHub renders them: GitHub heading anchors at every
# level, ``$...$``/``$$...$$`` math, and ``mermaid`` fences as diagrams.
myst_heading_anchors = 6
myst_heading_slug_func = "ldsc_docs.github_slug"
myst_enable_extensions = ["dollarmath", "colon_fence"]
myst_fence_as_directive = ["mermaid"]
# Resolve ``page.md#anchor`` links against documents and labels only.
myst_ref_domains = ["std"]

# API pages: one page per exported name, generated from NumPy-style docstrings.
# The autosummary blocks sit in ``eval-rst`` fences because the stub generator
# only finds reStructuredText directives.
autosummary_generate = True
autodoc_typehints = "none"
numpydoc_show_class_members = False
numpydoc_xref_param_type = False

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable", None),
    "pandas": ("https://pandas.pydata.org/docs", None),
    "scipy": ("https://docs.scipy.org/doc/scipy", None),
    "matplotlib": ("https://matplotlib.org/stable", None),
}

html_theme = "pydata_sphinx_theme"
html_title = "LDSC3"
html_static_path = []
html_theme_options = {
    "github_url": "https://github.com/WenbinWu2001/ldsc_py3",
    "navigation_with_keys": False,
    "show_toc_level": 2,
}
