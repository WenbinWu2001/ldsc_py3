# Developer Guide

Last updated on: 2026-10-27

This guide is for contributors to LDSC3.
The working rules for any change, including design confirmation, testing and documentation, are in [`AGENTS.md`](../../AGENTS.md) at the repository root.

## Set up and test

Use the editable development installation from [installation](../user-guide/installation.md#development-version) and run the tests from the repository root:

```bash
pytest
python -m unittest discover -s tests -p 'test*.py' -v
```

Run the two runners one after the other, not at the same time.
Tests do not add `src` to the Python path; they rely on the editable installation.
Golden NPZ files under `tests/fixtures/golden/` are fixed expectations: restore a missing one from version control instead of regenerating it.
The fixtures are described in [`tests/fixtures/README.md`](../../tests/fixtures/README.md).

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
