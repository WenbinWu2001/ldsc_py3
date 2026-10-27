> **Top priority:** Keep the user crystal clear about the work they are doing, why it is being done, and every procedural step.
> Never treat AI-generated content as a black box.
> Preserve this transparency while following the user's instructions; it supports rigorous, valid, and tractable code and analysis.

# LDSC Package Guidance

Last updated on: 2026-10-27

This repository, `ldsc_py3_restructured`, develops the distributable scientific Python package `ldsc` (distribution version `2.0b0`), with source under `src/ldsc/`.

## User instructions and skill conflicts

The user's instructions take precedence over guidelines provided in a skill. If explicit user instructions conflict with a skill's instructions, prioritize the user's instructions.

If a skill causes you to ask for permission or confirmation, pause, leave requested work unfinished, or diverge from the user's intent, name and link to the exact SKILL.md file you read, quote the relevant instruction, and briefly explain how it applies. Distinguish explicit skill requirements from your interpretation of guidelines.

## Environment and verification

Run from this repository root; environment creation and installation are setup commands, not prerequisites to repeat in an already configured environment:

```bash
conda env create -f environment.yml -n ldsc3-dev
conda activate ldsc3-dev
python -m pip install -e ".[dev]"
ldsc --help
python -m ldsc --help
# Focused package/CLI checks; select scientific tests for the affected workflow.
pytest tests/test_package_layout.py tests/test_cli_help.py
# Broader validation; run the two runners sequentially.
pytest
python -m unittest discover -s tests -p 'test*.py' -v
```

- Python support is `>=3.11,<3.14`; dependency bounds live in `setup.py`, `requirements.txt`, and `environment.yml`.
  Required dependencies are NumPy, pandas, SciPy, PyArrow, Matplotlib, threadpoolctl, and joblib.
  Extras: `plink` (bitarray), `liftover` (pyliftover), `test` (pytest), `all` (workflow extras), and `dev` (all plus tests).
- `pytest` is primary; `pyproject.toml` sets `testpaths = ["tests"]`, `-ra`, and markers `slow`, `statistical`, `io`, `file_format_compat`.
  Tests require the editable installation and do not inject `src`; run pytest and unittest sequentially.
  BED/gene interval projection uses NumPy from the base install, without an optional BED backend.
  Golden NPZs are immutable expectations: restore missing files from version control, never regenerate them as a test prerequisite.
- Packaging uses setuptools in `setup.py` and `MANIFEST.in`; no dedicated build command, CI workflow, linter, formatter, type checker, or documentation builder is configured.
  Do not invent mandatory commands or imply unrun checks passed.
  Before claiming completion, run relevant focused tests, broader tests when warranted, CLI/artifact/notebook checks, or documentation validation and report actual results and omissions.

## Package boundaries and current interfaces

- Import supported user-facing objects from `ldsc`; `src/ldsc/__init__.py` owns exports such as `AnnotationBuilder`, `LDScoreCalculator`, `RegressionRunner`, `SumstatsMunger`, `ReferencePanelBuilder`, `LDScoreSource`, configurations, writers, and convenience functions.
  `src/ldsc/_kernel/` and other underscore-prefixed modules are private implementations, even when selected objects are re-exported publicly.
- `ldsc=ldsc.cli:run_cli` and `python -m ldsc` expose one CLI: `infer-genome-build`, `annotate`, `ldscore`, `build-r2-panel`, `build-gene-ldscore-index`, `convert-ldsc2-ldscores`, `munge-sumstats`, `h2`, `partitioned-h2`, `quantile-h2`, `rg`, `query-r2`, `convert-h2-scale`, and `plot`.
  Keep CLI code limited to parsing/dispatch; numerical and plotting logic belongs to its workflow or kernel.
- Dependency direction is CLI → workflows → shared helpers/private kernels; workflows call output writers, and kernels must not import workflow orchestration or writers.
  Runtime reference adapters currently live in `_kernel/ref_panel.py`; `RefPanel.prepare_chromosome()` owns panel preparation and reader lifetime.
  Numerical projection consumes prepared, resolved inputs; do not add user-token discovery to numerical kernels or create circular imports.
- Reuse `path_resolution.py` for path tokens/output preflight, `column_inference.py` for aliases/identifier normalization, and `genome_build_inference.py` for hg19/hg38 and coordinate-basis inference.
  Exact paths, globs, and `@` suites are supported; PLINK additionally supports plain suite stems through its shared resolver.
  The installed package must be self-sufficient, without sibling-repository imports or root wrappers.
- Regression consumes aggregated LD-score artifacts without recomputing LD scores; plotting and scale conversion consume saved results without refitting.
  Matplotlib is required but lazily imported only for figure-producing paths.

## Design, implementation, and change discipline

- Default against backward compatibility unless explicitly requested: remove obsolete interfaces and update callers, tests, examples, and documentation to the confirmed contract.
  Do not add hypothetical aliases, adapters, migrations, fallback paths, or redundant validation; the established LDSC2 artifact exceptions below are deliberate requirements.
- Choose the simplest readable, maintainable, modular implementation; avoid speculative abstractions, configuration, and indirection unless significant computational gains or technical necessity justify them.
  Retain checks for scientific identity, alignment, convergence, and method assumptions; document non-obvious logic without restating code.
  Build capabilities incrementally on working end-to-end behavior, make durable architectural choices, and avoid disposable stopgaps.
- Prefer established libraries and proven patterns that reduce complexity; examine existing dependency documentation/types before assuming a capability is absent or adding a dependency.
  Ask for guidance if correct library use cannot be established.
- For changes beyond a several-line localized edit, discover repository facts and use `grilling`, `grill-me`, or `grill-with-docs` to resolve all non-trivial ambiguities before coding.
  For significant work, list important decisions and design details and obtain explicit design confirmation before implementation, including public API/file-format changes.
  Use `grill-me` for requests to keep asking; use `grill-with-docs` when authorized decision recording is wanted.
- After design confirmation and implementation authorization, use `tdd` for observable behavior with a stable test seam.
  Write the failing test first; numerical checks need known values, convergence properties, or reference implementations.
  Documentation-only work, exploratory analysis, and tests without meaningful signal follow the skill's exclusions and require appropriate concrete verification instead.
- Use `diagnosing-bugs` for concrete failures, establish the root cause before proposing repairs, and consult `lessons.md`.
  Stop and ask for direction after three unsuccessful repair attempts.
- Preserve unrelated user changes, source data, and generated artifacts.
  Ask before destructive changes, external publication, replacing existing rule files, or downloading external links; read-only browsing is allowed.
- Before structurally refactoring a file over 300 lines, remove dead properties, unused imports/exports, and debug logs and commit that cleanup separately.
  After major changes, commit or suggest a Conventional Commit: `<type>(<scope>): <description>`, subject ≤50 characters and body lines ≤72, explaining what/why with issue or breaking-change footers when applicable.
- After 10 or more messages or a gap, reread the request, confirmed decisions, active plan, `status.md` and `lessons.md` when present, and files to be edited.
  Do not put AI tool names in code comments, commit messages, PR bodies, or authorship, or emojis in any agent-authored content.

## Scientific and artifact contracts

- Make assumptions explicit, support scientific decisions with evidence and empirical validation, preserve reproducibility, and document uncertainty/limitations.
  Source data and packaged references are immutable unless changes are explicitly authorized; preserve relevant provenance, schemas, random seeds, numerical tolerances, and interpretation.
- Identity modes are `rsid`, `rsid_allele_aware`, `chr_pos`, and default `chr_pos_allele_aware`; global build defaults to `auto`.
  Preserve one-based canonical SNP positions, zero-based half-open BED intervals, allele orientation, duplicate-key policies, and hg19/hg38 inference boundaries.
  Annotation rows/columns must align to retained reference identities; `A1` is the signed-statistic allele, not necessarily the genome reference allele.
- Preserve reference-derived CM/MAF, the broad LD-contributor/count/overlap universe, and the separate filtered regression-row/weight universe.
  Native common counts use the inclusive `common_maf_min` threshold (default 0.05); legacy conversion retains strict LDSC2 common-count semantics.
  Regression defaults to common-SNP counts when available; do not silently substitute unavailable requested counts.
- Public annotation reads and LD-score Parquet values are float32; numerical LD accumulation, annotation sums, and overlap products use float64.
  Preserve Boolean interval unions, negative unbiased adjusted-r² values, one diagonal contribution, genomic ordering, and physical BIM/BED row mapping.
  Regression matrices retain `X: (n,p)`, `y/w: (n,1)`, aligned jackknife blocks, and estimator-specific assumptions; current IRWLS uses two reweightings, not a tolerance-based convergence loop.
  Preserve missingness semantics: zero is a valid annotation value, unknown support is not measured zero, and missing fields must not shift whitespace-parsed columns.
- `LDScoreDirectoryWriter` owns `metadata.json`, `ldscore.baseline.parquet`, optional single `ldscore.query.parquet` or numbered `ldscore.query.batchNNNNN.parquet` files, and `ldscore.overlap.parquet` when required.
  Query manifests record ordered files/columns/row groups; output tables span chromosomes, with chromosome row groups rather than public chromosome LD-score shards.
  Writing workflows return `LDScoreSource`; selected reads are explicit and caller-owned.
- Query annotations require explicit baseline annotations; no baseline and no query produces an all-ones `base` for ordinary unpartitioned generation.
  Index construction requires a supplied one-based gene catalog, an explicit hg19 or hg38 build assertion, an explicit SNP identifier (any of the four identity modes), and complete autosomes; indexed scoring inherits the complete immutable index and rejects live overrides or silent fallback.
- Reuse `outputs.py` directory writers, their `artifact_family()` declarations, and `path_resolution.py` preflight helpers before owned writes, including logs.
  Core materializing CLI commands require `--output-dir`; explicit `--overwrite` replaces only owned outputs and reconciles stale owned siblings, preserving unrelated files.
  Derived CLI destinations are `<result-dir>/plots/` and `<h2-result-dir>/postprocessing/liability-scale/`, without CLI overrides.
- Put headline results at the result root and diagnostics/logs under `diagnostics/`; LD-score root metadata is its separate artifact manifest.
  Use Parquet for new internal artifacts and TSV for scientific result tables, with documented legacy/sidecar exceptions.
  Logs and failure markers are not scientific outputs or entries in result `output_paths`.
- Failed authorized overwrites leave `RUN_FAILED.txt` (or `RUN_FAILED.chr<chrom>.txt` for concrete panel shards); successful retries remove the applicable marker.
  Markers add no rollback, quarantine, restoration, or action-order changes and do not enter metadata.
  Help and `munge-sumstats --infer-only` create no directories/logs/artifacts/markers; infer-only still requires `--output-dir`.
  Gene-index publication retains its existing complete atomic transaction, with open logs in hidden sibling build-state and closed successful logs moved into published diagnostics.
- Preserve explicit LDSC2 interoperability: regression accepts `.sumstats(.gz)` with SNP/A1/A2/Z/N and validates allele orientation; munging may export `.sumstats.gz`.
  Annotation workflows read `.annot(.gz)` and write `.annot.gz`; omitting both alleles permits base-identity matching even in allele-aware workflows, while supplying only one allele is invalid.
- Legacy LD-score input requires `convert-ldsc2-ldscores` with complete chromosome 1–22 reference/weight `<prefix><chrom>.l2.ldscore(.gz)` families, required `.l2.M_5_50`, and optional `.l2.M` under the documented count policy.
  Converted identity is allele-unaware; `.w.l2.ldscore(.gz)`, query/thin conversion suites, and direct regression consumption of fragments are unsupported.
  Private kernels do not own legacy sumstats, LD-score, count, or annotation-group output files; full contracts live in `docs/reference/formats/legacy-sumstats.md` and `docs/reference/formats/legacy-ldscores.md`.
- Avoid routine provenance hashes; prefer schemas, source labels, declared metadata, and counts.
  Retain hashes only for necessary semantic identity/integrity guarantees, documenting their purpose; the retained hashes are the R² panel sidecar binding (`docs/reference/formats/r2-panels.md`, section 4) and the gene LD-score index identity: `index_id`, its input fingerprints, and `published_row_metadata_sha256` (`docs/reference/gene-ldscore-index.md`).

## Existing warning, error, and diagnostic policy

- Preserve the repository's current policy; warning-visibility defaults from initialization skills/templates do not override it.
  `src/ldsc/_logging.py` and `docs/current/workflow-logging.md` define file-authoritative workflow logs, CLI stderr errors only during logged workflows, and no package-installed Python API console handler.
  CLI calls without a workflow log retain the existing verbosity fallback; do not introduce blanket warning-display or suppression changes.
- Assume users do not inspect file logs after a successful run, especially under batch schedulers such as SLURM.
  If an unrequested input-resolution, filtering, ambiguity, or compatibility condition can silently change the scientific analysis, fail before substantial computation with a nonzero, actionable error and machine-readable diagnostics.
  If it is unclear whether an operation warrants an error or warning, or whether it could affect the scientific analysis, ask the user to decide.
  Explicitly requested transformations may succeed when their effect is recorded in the result diagnostics and provenance.
- Use exhaustive, stage-specific validation gates that collect every safely discoverable issue across focal/control inputs before aborting.
  An unreadable structural contract may fail immediately, a referenced analysis-changing defect must fail, and an unreferenced latent defect may be diagnosed without failing a correct run.
  Design diagnostics for concrete identification and repair; preserve explicit query-continuation policies and their attempted-query records, without introducing silent recovery.
- Keep run-aborting messages self-contained; errors with three or more distinct causes also require the corresponding command section in `docs/user-guide/troubleshooting.md` and synchronized in-code anchors.
  Guide users toward their objective without requiring knowledge of internal implementation.

## Requested code reviews

These rules apply only when a review is requested; they do not require an extra review after every change.

- Establish and verify the supplied diff, commit range, branch/base or PR, including staged/unstaged/untracked scope for work in progress and the base/merge base for branches.
  State assumptions, report empty diffs, and distinguish reviewed state from later edits.
- Read applicable guidance, surrounding code, callers, and tests; assess concrete defects, regressions, scientific correctness, error handling, tests, and repository consistency alongside fulfillment of the request, missing requirements, unsupported assumptions, contract mismatches, and scope.
  These are distinct questions, not mandated separate passes or a fixed format.
- Ground requirement findings in precise quotations/references to user instructions, confirmed decisions, or governing records; identify uncertain or unavailable requirements rather than inventing them.
  Distinguish defects from preferences and respect confirmed decisions, explaining evidence of contradictions.
- Lead with actionable findings by severity, with file/line, evidence, impact, and remedy; identify open questions, unreviewed areas, and checks not run.
  No findings does not establish correctness of untested paths; review-only authorization permits findings, not edits to code, issues, or PRs.

## Documentation, recording permissions, and evidence

- Maintain affected docstrings, README workflow/input/output descriptions, current contracts, tutorials/wiki guidance, and validation artifacts within authorized change scope; update manuscripts/reports only when requested.
  Core-function documentation should define inputs, outputs, defaults, and caveats; use `scientific-python-docs` for requested or change-affected scientific documentation, not unrelated missing-docstring cleanup.
- Before changing responsibilities, public entry points, dependencies, main data flows, or scientific constraints, read relevant sections of [docs/dev/architecture.md](docs/dev/architecture.md).
  Read the whole architecture for broad relational refactors; internal refactors leaving these facts unchanged require neither a full reread nor an architecture update.
  Update descriptions made inaccurate within authorized scope using `architecture-doc`; preserve this established location and still-valid content.
- Every edited/created Markdown document includes and refreshes `Last updated on: YYYY-MM-DD` near its top.
  Put each Markdown/LaTeX prose sentence on its own source line, preserving required syntax, exact quotations, and protected originals; do not reformat unrelated files.
  Markdown/Quarto/notebook math uses `$...$` inline and `$$...$$` display; LaTeX uses `$...$` inline and `\[...\]` display; chat uses `\(...\)` inline and `$$...$$` display.
- Read [docs/dev/DOCS_VOCABULARY.md](docs/dev/DOCS_VOCABULARY.md) before every writing task, including small edits; newly written documentation, docstrings, CLI help text, and messages use its approved terms.
  Only the user approves, edits, or removes entries; record new candidate terms under its proposed terms, and flag existing violations with location and preferred term instead of rewriting them without authorization.
- Keep working rules here and records in their established destinations; update this documentation map with real paths when adding/moving design documents or tutorials.

| Destination | Responsibility and recording boundary |
| --- | --- |
| `docs/dev/code-structure.md`, `docs/dev/architecture.md`, `docs/dev/data-flow.md`, `docs/current/layer-structure.md` | Contributor navigation, execution/dependency architecture, artifact flows, and ownership; current contracts govern unless a confirmed change supersedes them. |
| `docs/reference/` | Exact command, format, result, and method contracts; current contracts govern unless a confirmed change supersedes them. |
| `docs/dev/design/`, `docs/dev/guidelines/` | Implemented design guides (memory and CPU, configuration, plotting) and contributor writing rules (CLI help, user docs). |
| `docs/dev/DOCS_VOCABULARY.md` | The one glossary: names, approved terms with identifiers and symbols, rules for writers, and proposed terms; the established terminology home for `domain-modeling` and decision records. |
| `docs/current/*-decisions.md` | Preserve existing confirmed decisions; retain dated topic names and do not duplicate records. |
| `docs/specs/`, `docs/plans/` | Temporary folders for in-flight design specifications and implementation plans only. Once the work is implemented, move lasting rationale into `docs/decisions/` (or the relevant `docs/current/*-decisions.md`) and the contract into `docs/current/`, then delete the spec and plan. |
| `docs/decisions/` | New authorized standalone `grill-with-docs` decision records with confirmed rationale; use existing topical records where appropriate, keep terminology in its established home, and use ADR format only under `domain-modeling` criteria. |
| `docs/decisions/2026-10-11-parallel-execution.md` | Confirmed parallel-execution decisions; implemented runtime contract in `docs/reference/parallel-execution.md`. |
| `status.md` | Authorized current work, constraints, blockers, next steps; refresh in place and link decisions/evidence. |
| `history.md` | Authorized meaningful events, appended chronologically with dates; consult on demand, not automatically at startup. |
| `lessons.md` | After a repeated mistake (two or more) or an investigated non-obvious bug, append summary/root cause/correction; create on the first qualifying entry, skipping trivial slips. |
| `docs/user-guide/` (including `troubleshooting.md`), `tutorials/`, `docs/wiki/` | User workflows, runnable Markdown/notebook examples, and command-specific remediation; `tutorials/` and `docs/wiki/` are being merged into `docs/user-guide/`. |
| Validation evidence (no dedicated folder) | Summarize durable measurements inline in the host `docs/reference/`, `docs/dev/`, or `docs/decisions/` document with date, revision, and limits; keep raw logs, fixtures, and probe scripts out of `docs/`. The audits, specs, plans, and archive removed on 2026-10-26 remain in git history at `479fc4d`. |

- Use `project-state` only for useful, authorized status/history updates: these destinations do not authorize completion logs, writes during read-only work/side conversations, or empty-file creation.
  `docs/decisions/` contains authorized decision records; `status.md` exists, while `history.md` remains a convention, currently absent.
  Their use does not bypass lessons or decision-recording rules.
- Anchor substantive workflow/scientific explanations to code paths and lines or functions/classes, and document/paper claims to sections, tables, figures, or appendices.
  Legacy comparisons must cite both implementations, using the legacy repository below.

## Resources

- `README.md`, `setup.py`, `pyproject.toml`, `requirements.txt`, `environment.yml`, and `MANIFEST.in`: user surface, environment, package/test configuration, and distribution contents.
- `src/ldsc/data/hm3_curated_map.tsv.gz`: packaged default HM3 restriction map; `src/ldsc/data/hm3_chr_pos_reference.tsv.gz`: compact inference/test reference; `src/ldsc/data/regions/`: hg19/hg38 MHC and centromere presets.
  `src/ldsc/data/readme.txt` and `src/ldsc/data/ATTRIBUTION.txt` record resource roles and unresolved source/reuse questions; do not invent dataset releases or licenses.
- `tests/fixtures/README.md`, `tests/fixtures/golden/`, and `tests/fixtures/minimal_external_resources/README.md`: immutable expectations and deterministic 32-SNP chromosome-22 1KG-30x subset.
  The minimal fixture guide references optional external `resources/` (including a liftover chain), absent in this checkout; do not claim external-resource checks ran without those inputs.
- `docs/reference/formats/column-schema.md`, `docs/reference/snp-identity-and-genome-build.md`, `docs/reference/ldscore-snp-universes.md`, and `docs/reference/regression-options.md`: schema, identity/build, SNP universes, and estimator contracts.
- `docs/dev/design/memory-and-cpu.md`, `docs/reference/gene-ldscore-index.md`, `docs/reference/methods/gene-ldscore-index-math.md`, `docs/reference/results/partitioned-h2.md`, and `docs/reference/results/quantile-h2.md`: memory ownership, exact indexes, result interpretation, and post-fit quantiles.
- `docs/reference/parallel-execution.md`: managed task workers, joblib CPU budgets, one BLAS thread per worker, and chromosome diagnostics.
- `docs/reference/input-paths.md`, `docs/current/workflow-logging.md`, and `docs/dev/guidelines/cli-help.md`: resolution, output/diagnostic lifecycle, and task-oriented CLI help.
- `tutorials/ld-score-calculation.md`, `tutorials/partitioned-ldsc.md`, `tutorials/cell-specific-ldsc.md`, and `tutorials/convert-legacy-ldscores.md`: main worked workflows, complemented by notebooks and `docs/user-guide/guided-tutorial.md`.
- `docs/reference/ldsc2-flag-map.md`: LDSC2 interface mapping, including options without a current equivalent; no document yet records numerical equivalence with LDSC2. `docs/dev/release.md`: maintainer release runbook, not authorization to release.
  `LICENSE`, `NOTICE`, and `CITATION.cff` preserve licensing, attribution, and citation metadata; retain notices during refactors and packaged resource attribution.
- `/Users/wenbinwu/Documents_local/Research/SullivanLab/LDSC/repos/ldsc_py2_Bulik_workspace/ldsc_py2_Bulik`: legacy numerical/format reference.
- `/Users/wenbinwu/Documents_local/Research/SullivanLab/LDSC/repos/ldsc_py3_Jerry_workspace/docs/ldsc_papers`: LDSC papers and supplements.
