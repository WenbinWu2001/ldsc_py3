> **Top priority:** Keep the user crystal clear about the work, why it is being done, and every procedural step. Never treat AI-generated work as a black box. This transparency overrides other rules and supports rigorous, valid, and tractable code and analysis.

# LDSC Package Guidance

Last updated on: 2026-08-02

`ldsc_py3_Jerry` is a refactored, distributable Python 3 LDSC package, not an analysis repository. Preserve its package layout, public interfaces, CLI contracts, canonical artifact formats, and compatibility boundaries.

## Environment and verification

From the repository root, create the development environment and install the editable package with development extras:

```bash
conda env create -f environment.yml -n ldsc3-dev
conda activate ldsc3-dev
python -m pip install -e ".[dev]"
```

Use the real entry points:

```bash
ldsc --help
python -m ldsc --help
pytest
python -m unittest discover -s tests -p 'test*.py' -v
```

`pytest` is the primary suite. The standard-library unittest command remains a compatibility check during the transition. `pyproject.toml` configures pytest with `tests` as the test path, `src` on `pythonpath`, `-ra`, and the `slow`, `statistical`, `io`, and `file_format_compat` markers. No separate formatter, linter, type checker, documentation builder, or build command is configured; do not invent one.

The package supports Python 3.11 through 3.13. Core dependencies are NumPy, pandas, SciPy, and PyArrow; extras provide PLINK (`bitarray`), BED (`pybedtools` plus external `bedtools` on `PATH`), liftover (`pyliftover`), and tests (`pytest`). Use the constraints in `setup.py`, `requirements.txt`, and `environment.yml` rather than broadening them casually.

## Package structure and public contracts

- `src/ldsc/` is the supported public Python package surface. Import stable user-facing objects from `ldsc`; `src/ldsc/_kernel/` is private low-level numerical and file-format implementation code.
- Keep one CLI surface: `ldsc` with `annotate`, `build-ref-panel`, `ldscore`, `munge-sumstats`, `h2`, `partitioned-h2`, `rg`, and `query-r2`. `ldsc.cli` dispatches; it must not gain numerical logic.
- Public workflow modules own user-facing path resolution, header inference, global configuration, genome-build inference, and output preflight. Kernel modules receive resolved primitive inputs and perform numerical work or low-level parsing. Keep dependencies unidirectional and avoid circular imports.
- Keep alias and identifier normalization centralized in `column_inference.py`; keep hg19/hg38 and 0-based/1-based inference centralized in `genome_build_inference.py`. Prefer existing workflow objects, config dataclasses, path-resolution helpers, and column-inference registries over one-off parsing or normalization.
- Treat public Python exports, CLI flags, package-written schemas, and legacy file formats as compatibility contracts. Public API or file-format changes require a deliberate compatibility decision, not mechanical cleanup.
- The package must remain self-sufficient: do not import from sibling repositories or repository-root wrappers.

## Workflow and artifact invariants

- Public LD-score output is a canonical directory written by `LDScoreDirectoryWriter`: `manifest.json`, `baseline.parquet`, and optional `query.parquet`. Regression consumes this aggregated directory and must not recompute LD scores.
- New artifact-writing workflows use `--output-dir`, `--overwrite`, and `--log-level`, write headline outputs at the directory root, and write `diagnostics/metadata.json` plus `<command>.log` through `workflow_logging`. Reuse the shared `*DirectoryWriter` classes in `outputs.py` with `ensure_output_directory` and `preflight_output_artifact_family`.
- A command may emit a clean, pipeable TSV to stdout when no output directory is supplied, as `query-r2` does. Self-describing Parquet outputs embed provenance in the Parquet footer rather than writing `metadata.json`, as `munge-sumstats` does.
- Legacy LDSC formats are compatibility boundaries, not the public LD-score writer layout: annotation workflows read and write `.annot(.gz)`, munging can write `.sumstats.gz`, and the kernel supports `.l2.ldscore(.gz)`, `.w.l2.ldscore(.gz)`, `.M`, and `.M_5_50`.
- Query annotations require explicit baseline annotations. The synthetic all-ones `base` annotation is only for ordinary unpartitioned LD-score generation.
- Preserve the original LDSC regression default of using `.M_5_50`-style common-SNP counts when available.
- Treat source data and packaged reference resources as immutable inputs unless the user explicitly authorizes changes. Preserve identifier mode, genome build, coordinate basis, allele conventions, provenance, schemas, random seeds, numerical tolerances, and scientific interpretation contracts relevant to the requested change.

## Change discipline

- Keep code concise. Minimize unnecessary validation, guard clauses, and redundant design layers. Document non-obvious logic; do not add comments that merely restate code.
- Preserve unrelated user changes and generated artifacts outside the requested scope. Ask before destructive changes, external publication, replacing an existing rule file, or downloading from an external link.
- Before a structural refactor of a file over 300 lines, remove dead properties, unused imports or exports, and debug logs; keep that cleanup separate from the refactor.
- When the same class of mistake occurs two or more times, or a non-obvious bug requires real investigation, append a one-line summary, root cause, and correction to root `lessons.md`. Skip one-off typos and trivial slips.
- After 10 or more messages, or after resuming from a gap, reread the active plan in `docs/plans/`, `lessons.md` when present, and the source files to be edited.
- After a major change, commit with a meaningful Conventional Commit message, or remind the user to commit: `<type>(<scope>): <description>`. Keep the subject at most 50 characters and body lines at most 72 characters; explain what and why, not how.
- Do not use AI tool names in code comments, commit messages, PR bodies, or authorship. Do not use emojis in documentation, docstrings, Markdown files, reports, tutorials, or manuscripts.
- Do not hard-wrap Markdown or LaTeX prose. Use `\(...\)` for inline math and `$$ ... $$` for displayed math in Markdown and notebook Markdown; use `$...$` for inline math in `.tex` sources.

## Coding style
- Do not preserve backward compatibility. Remove obsolete paths instead of adding compatibility layers, fallbacks, or migrations.
- Choose the simplest implementation that fully meets the current requirements. Avoid speculative abstractions, configuration, and indirection.
- Grow the system in layers. Start from the smallest version that works end to end, and add each new capability on top of a product that already works. Never trade a working product for unfinished complexity.
- Keep components modular and concerns clearly separated.
- Prefer established, well-maintained libraries when they reduce overall complexity or improve reliability. Do not reimplement common functionality without a clear reason.
- Lean on the dependencies already in the project before writing your own implementation or adding packages. Do not assume a library lacks a capability without checking its documentation and types.
- Make architectural decisions for the long term. Do not accept a stopgap that only works for now and is meant to be replaced later.
- Study how established products solve the problem before designing a solution. Adopt their proven patterns and conventions rather than inventing an approach from scratch.


## Workflow skills

- For work beyond a several-line localized edit, use `/grilling` or its appropriate route to resolve non-trivial ambiguities before coding; keep asking until they are resolved.
- Route requests to “keep asking,” “resolve all ambiguities,” or “resolve before implementing” to `grill-me`. Use `grill-with-docs` when the user also asks to record decisions in repository documentation.
- Use `to-spec` when the user asks to write a specification or capture an agreed design as documentation. Use `implementation-plan` when the user asks for an execution or implementation plan. Store their default artifacts at `docs/specs/YYYY-MM-DD_<topic>.md` and `docs/plans/YYYY-MM-DD_<topic>.md`.
- For implementation with observable behavior and a stable test seam, use `tdd`: write the failing test first, then implement and refactor. Numerical tests must use known input/output values, convergence properties, or reference implementations—not only successful execution.
- Use `diagnosing-bugs` for a concrete failure. Diagnose root cause before proposing a fix, check `lessons.md` for related failures, and stop for user direction after three unsuccessful repair attempts.
- Use `two-axis-code-review` only when the user explicitly requests a fixed-diff review.

## Documentation and verification

- `docs/current/` is the active package design and navigation source. Treat its architecture, layer, data-flow, path, configuration, artifact, and workflow contracts as authoritative unless a user-approved change updates them.
- After code or artifact-contract changes, update affected docstrings, relevant `docs/current/` documentation, README workflow/input/output descriptions, and applicable tutorials. Update papers or reports only when requested.
- Every Markdown document created or edited by the agent includes `Last updated on: YYYY-MM-DD` near the top and updates it on later changes.
- Use `scientific-python-docs` whenever a public Python function, class, or module header lacks a docstring or has incomplete or outdated documentation. Use `architecture-doc` when package structure changes significantly, a reusable module is added, or contributors need a refreshed overview.
- For a run-aborting error with three or more distinct causes, keep the message self-contained and add or update the corresponding command section in `docs/troubleshooting.md`; keep any in-code anchor link synchronized.
- Before claiming work is complete or fixed, run the relevant focused test, full suite, CLI command, artifact validation, or notebook check and report the actual result. When no automated check exists, say so and report the concrete checks performed.

## Citations

- When explaining a workflow, pipeline, interface, or mechanism, cite the source file and line number or function/class name.
- When summarizing or quoting a repository document, tutorial, audit, plan, or specification, cite its section, table, figure, or appendix.
- When comparing behavior with legacy LDSC, cite the relevant path and function in `/Users/wenbinwu/Documents_local/Research/SullivanLab/LDSC/repos/ldsc_py2_Bulik_workspace/ldsc_py2_Bulik` as well as the restructured implementation.

## Resources

- `README.md`: installation, package surface, output policy, and user-facing CLI overview.
- `setup.py`, `pyproject.toml`, `requirements.txt`, and `environment.yml`: package metadata, Python/dependency constraints, extras, test configuration, and conda environment.
- `src/ldsc/__init__.py`: supported public Python exports; `src/ldsc/cli.py`: unified CLI; `src/ldsc/_kernel/`: private numerical and file-format implementation.
- `src/ldsc/data/`: packaged HM3 maps and hg19/hg38 region resources used by workflows.
- `tests/`: primary behavioral, numerical, I/O, format-compatibility, and workflow suite; `tests/fixtures/` and `tests/fixtures/minimal_external_resources/`: documented deterministic fixture resources.
- `tutorials/`: package-level Markdown and notebook usage examples.
- `docs/current/`: active architecture, data-flow, configuration, schema, provenance, logging, and workflow-contract documentation.
- `docs/specs/`: agreed design specifications; `docs/plans/`: implementation plans using the dated topic convention; `docs/archive/`: historical context only.
- `docs/audits/legacy-equivalence/`: evidence and progress records for compatibility with the legacy implementation.
- `docs/troubleshooting.md`: command-organized remediation reference; `docs/release.md`: release runbook for maintainers.
- `/Users/wenbinwu/Documents_local/Research/SullivanLab/LDSC/repos/ldsc_py2_Bulik_workspace/ldsc_py2_Bulik`: legacy LDSC codebase used for numerical and compatibility comparisons.
- `/Users/wenbinwu/Documents_local/Research/SullivanLab/LDSC/repos/ldsc_py3_Jerry_workspace/docs/ldsc_papers`: main text and supplements of relevant LD score analysis papers.