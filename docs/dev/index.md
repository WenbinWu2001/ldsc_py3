# Developer Guide

Last updated on: 2026-11-07

This guide is for contributors to LDSC3.
The working rules for any change, including design confirmation, testing and documentation, are in [`AGENTS.md`](../../AGENTS.md) at the repository root.

## Set up and test

Use the editable development installation from [installation](../user-guide/installation.md#development-version) and run the tests from the repository root:

```bash
pytest
```

Tests do not add `src` to the Python path; they rely on the editable installation.
Golden NPZ files under `tests/fixtures/golden/` are fixed expectations: restore a missing one from version control instead of regenerating it.
The fixtures are described in [`tests/fixtures/README.md`](../../tests/fixtures/README.md).

## Build the documentation site

The Markdown pages under `docs/` are also built into a website with Sphinx, the PyData theme and MyST-Parser.
Install the `docs` extra and build from the repository root:

```bash
python -m pip install -e ".[docs]"
sphinx-build -W --keep-going -b html docs docs/_build/html
```

Open `docs/_build/html/index.html` to read the result; `-W` turns every warning, such as a broken link or anchor, into a failure.
The [`Documentation` workflow](../../.github/workflows/docs.yml) runs the same build on every push and pull request and publishes the site to GitHub Pages from `restructure`.

- **Pages and navigation.** Each section `index.md` ends with hidden `toctree` blocks that set the site navigation; a new page must be added to one of them.
  `exclude_patterns` in [`docs/conf.py`](../conf.py) keeps contributor-only pages, such as the glossary and the writing guidelines, out of the site.
- **Generated reference pages.** The option lists in the [command-line reference](../reference/cli.md) come from the parsers, and the [Python API reference](../reference/python-api.md) generates one page per exported name from its docstring.
  A newly exported name must be added to one of the `autosummary` lists there; `tests/package/test_docs_api_reference.py` fails until it is.
- **Links.** Pages keep relative links that work on GitHub.
  During the build, [`docs/_ext/ldsc_docs.py`](../_ext/ldsc_docs.py) points links to files that are not site pages (package source, tests, excluded pages) at the `restructure` branch on GitHub, and gives headings the anchors GitHub gives them.

## Where to start

| Page | Use it to |
| --- | --- |
| [Code structure](code-structure.md) | Find where a change belongs: module map, change guide, test map. |
| [Architecture](architecture.md) | Understand layers, ownership, cross-cutting concerns and invariants. |
| [Data flow](data-flow.md) | Follow the files each workflow reads and writes. |

## Design guides

| Page | Content |
| --- | --- |
| [Memory and CPU](design/memory-and-cpu.md) | Storage, streaming, reuse and numerical-buffer designs, with their precision and lifetime invariants. |
| [Configuration](design/configuration.md) | Configuration objects and provenance checks. |
| [Plotting](design/plotting.md) | Plot dispatch, liability-scale conversion and derived outputs. |
| [Logging](design/logging.md) | Rules for log, warning and error messages, and the logging API boundary. |

## Writing guidelines

- [Documentation vocabulary](DOCS_VOCABULARY.md): the one glossary; read it before writing any documentation, docstring, help text or message.
- [CLI help](guidelines/cli-help.md): how to write `--help` text and name options.
- [Writing user docs](guidelines/writing-docs.md): conventions for user-guide pages.
- Markdown files carry a `Last updated on:` line and put each prose sentence on its own line ([`AGENTS.md`](../../AGENTS.md)).

## Validation and release

- [Validation](validation/index.md): what has been validated, with which evidence, and what has not.
- [Release](release.md): the maintainer's procedure for releasing `restructure` to `main`.

Decision records go in `docs/decisions/` only when the maintainer asks for one; otherwise a decision's consequences and any rationale that protects a non-obvious constraint go into the page that documents the behavior.

```{toctree}
:hidden:

code-structure
architecture
data-flow
design/memory-and-cpu
design/configuration
design/plotting
design/logging
validation/index
release
```
