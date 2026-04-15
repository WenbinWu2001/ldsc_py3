# PLANS.md — Purpose and Structure

`PLANS.md` is a single file that tracks all active and future work in the repository. It is the authoritative record of what is being built, why, and how far along each effort is.

**File organization:** Active plans appear first, followed by an `## Archive` section at the bottom for completed or abandoned plans. Do not create separate plan files or folders — everything lives in this one file.

**When to write a plan:** Write an ExecPlan for any new feature, algorithm, or refactor with meaningful unknowns. Skip it for small, self-contained changes (e.g., fixing a typo, bumping a dependency, adding a one-liner utility).

**Goal of an ExecPlan:**

- Guide the model to implement precisely what the plan specifies
- Provide a clear record of decisions and progress the user can follow

---

# Restructure `ldsc_py3_Jerry` Into Layered Workflow Modules

This ExecPlan is a living document. Keep Progress and Decision Log up to date as work proceeds.

## Purpose

Restructure `ldsc_py3_Jerry` into the layered, CLI-first design defined by `architecture.md`, `code_structure.md`, and `class-and-features.md` while preserving the current numerical behavior and output files of the legacy LDSC workflows. After this work, the CLI entry points (`ldsc.py`, `ldsc_new.py`, `munge_sumstats.py`, `make_annot.py`, and `utils/run_bed_to_annot.py`) should become thin wrappers that build config/domain objects and delegate to service modules in `ldscore/`. Verification requires parity tests for LD-score generation, summary-statistics munging, h2/rg/h2-cts workflows, annotation generation, and the new batch partitioned-h2 feature.

## Progress

- [ ] Milestone 0: Freeze legacy behavior with local fixtures and parity tests.
- [ ] Milestone 1: Add shared config, identifier, and domain scaffolding without changing behavior.
- [ ] Milestone 2: Build the flexible output module with registry-based artifact producers.
- [ ] Milestone 3: Build the annotation module and migrate annotation-generation entry points.
- [ ] Milestone 4: Build the reference-panel module and loaders for PLINK and parquet R2.
- [ ] Milestone 5: Build the LD-score calculation module and cut over LD-score entry points.
- [ ] Milestone 6: Build the summary-statistics munging module and cut over `munge_sumstats.py`.
- [ ] Milestone 7: Build the regression module, preserve legacy defaults, and add batch partitioned h2.
- [ ] Milestone 8: Cut over remaining CLIs, clean compatibility wrappers, and run full validation.

## Plan of Work

The restructuring must be incremental and behavior-preserving. Do not perform a big-bang rewrite. First create a local test harness and copied fixtures inside this repository so the refactor does not depend on the sibling `ldsc_py2_Bulik` tree at runtime. Then add new workflow-facing modules alongside the existing scripts and legacy kernels. Only after tests exist for a feature should its CLI or legacy workflow entry point be redirected to the new module.

When the design documents conflict with the current implementation, follow the design documents for public interfaces, class boundaries, config names, and workflow structure. When the current code conflicts only at an internal implementation detail, prefer reusing the current numerical code until parity is proven, then clean it up later if needed. Each milestone below must end with passing tests and one isolated commit before the next milestone begins.

### Milestone 0: Freeze legacy behavior with tests

Create a proper `tests/` tree inside the repository and copy the minimum required fixtures from `../ldsc_py2_Bulik_workspace/ldsc_py2_Bulik/test/` into `tests/fixtures/legacy/`. Port or adapt the old kernel tests for parse, jackknife, IRWLS, regressions, LD-score helpers, and munging so the refactor has a local oracle. Add CLI-level parity tests for:

- `ldsc.py --l2`
- `ldsc_new.py`
- `munge_sumstats.py`
- `ldsc.py --h2`
- `ldsc.py --rg`
- `ldsc.py --h2-cts`
- `make_annot.py`
- `utils/run_bed_to_annot.py`

Use the current `ldsc_py3_Jerry` behavior as the baseline when a direct checked-in oracle is unavailable. Do not start moving code before this test harness exists.

### Milestone 1: Shared config, identifier, and domain scaffolding

Create these new modules:

- `ldscore/config.py`
- `ldscore/identifiers.py`

Add the config dataclasses from `class-and-features.md`:

- `CommonConfig`
- `AnnotationBuildConfig`
- `RefPanelConfig`
- `LDScoreConfig`
- `MungeConfig`
- `RegressionConfig`

Put shared SNP-identifier helpers in `ldscore/identifiers.py`, including:

- canonical `snp_id` builders for `rsid` and `chr_pos`
- flexible SNP-column inference using the current `munge_sumstats.py` alias-cleaning style
- global restriction-list readers
- chromosome normalization helpers

Do not create one giant `models.py`. Keep domain/result classes next to their feature modules later in the plan. The only shared modules at this stage should be config and identifier utilities.

Validation for this milestone:

- unit tests for config defaults and validation
- unit tests for identifier/header inference (`rsid`, `rsID`, `RSID`, `SNPID`, `snp_id`, `SNP`)
- no CLI behavior changes yet

### Milestone 2: Flexible output module

Create:

- `ldscore/output.py`

Implement the output classes and interfaces defined in `class-and-features.md`:

- `OutputSpec`
- `ArtifactConfig`
- `RunSummary`
- `ArtifactProducer`
- `OutputManager`
- `ResultWriter`
- `ResultFormatter`
- `PostProcessor`

Lock these design decisions now:

- `OutputSpec` stays small and CLI-friendly
- `enabled_artifacts=None` means “use the built-in default artifact set”
- future summary, report, and plot outputs are added through registered `ArtifactProducer` implementations
- `OutputManager` must remain orchestration-only rather than becoming a hardcoded dispatcher
- `ArtifactConfig` is reserved for future advanced options so plot/report tuning does not bloat the base output config

Build producers for the current built-in artifacts first:

- LD-score tables
- regression-weight LD-score tables
- count files
- annotation manifest
- summary TSV/JSON
- run metadata dump

Validation for this milestone:

- unit tests with fake result objects
- tests for built-in defaults versus explicit `enabled_artifacts`
- tests proving a newly registered producer can emit an extra artifact without modifying `OutputManager`

### Milestone 3: Annotation module

Create:

- `ldscore/annotation.py`

Define and implement:

- `AnnotationSourceSpec`
- `AnnotationBundle`
- `AnnotationBuilder`

Move or reuse logic from:

- `utils/run_bed_to_annot.py`
- `make_annot.py`
- the annotation-loading pieces of `ldscore/ldscore_new.py`

The annotation module must:

- keep baseline and query annotations separate but aligned
- support chromosome-split and non-split `.annot` inputs
- project BED inputs to SNP-level annotations
- apply the global SNP restriction from `CommonConfig`
- preserve exact row order/alignment invariants

After tests pass, change `make_annot.py` and `utils/run_bed_to_annot.py` to become thin wrappers around the new module.

Validation for this milestone:

- parity tests for `.annot` outputs against legacy fixtures
- tests for mismatched annotation row order failure
- tests for global restriction application
- tests for baseline/query column grouping metadata

### Milestone 4: Reference-panel module

Create:

- `ldscore/ref_panel.py`

Define and implement:

- `RefPanelSpec`
- abstract `RefPanel`
- `PlinkRefPanel`
- `ParquetR2RefPanel`
- `RefPanelLoader`

Reuse existing low-level code from:

- `ldscore.parse`
- `ldscore.ldscore.PlinkBEDFile`
- `ldscore.ldscore_new.SortedR2BlockReader`
- parquet frequency/metadata helpers in `ldscore.ldscore_new`

Keep metadata loading lazy and chromosome-scoped. `RefPanel` should own metadata access and reader creation, but fully materialized retained metadata belongs on final result objects, not on the panel itself.

Validation for this milestone:

- tests for chromosome discovery
- tests for lazy `load_metadata(chrom)`
- tests for `filter_to_snps(chrom, snps)`
- tests for parquet frequency/MAF metadata merge
- tests for ambiguous `chr_pos` matching failure

### Milestone 5: LD-score calculation module

Create:

- `ldscore/ldscore_workflow.py`

Define and implement:

- `ChromLDScoreResult`
- `LDScoreResult`
- `LDScoreCalculator`

Reuse the current LD-score math from:

- `ldscore/ldscore.py`
- `ldscore/ldscore_new.py`
- `ldsc.py` genotype-side workflow

This module must preserve the current behavior while matching the new design:

- chromosome-wise computation first, combined result second
- one global regression-SNP input partitioned internally by chromosome
- `M` and `M_5_50` stored as named count arrays rather than bare fields
- regression defaults still use `M_5_50`-style counts later
- `LDScoreResult` stores aggregated metadata/tables plus `chromosome_results`
- output writing goes through the new output layer

Implementation detail rule:

- it is acceptable in the first cut to reuse the current PLINK and parquet backends as-is internally if the public structure matches the design docs and parity is preserved
- once the calculator is stable, `ldscore/ldscore_new.py`, `ldsc_new.py`, and `ldsc.py --l2` should route through it

Validation for this milestone:

- PLINK-mode parity tests against current `ldsc.py --l2`
- parquet-mode parity tests against current `ldscore/ldscore_new.py`
- tests for aggregated versus per-chromosome outputs
- tests for `w_ld`
- tests for `.M` and `.M_5_50`
- tests for `compute_m5_50` affecting counts only

### Milestone 6: Summary-statistics munging module

Create:

- `ldscore/sumstats_munger.py`

Define and implement:

- `RawSumstatsSpec`
- `SumstatsTable`
- `MungeRunSummary`
- `SumstatsMunger`

Move or reuse logic from `munge_sumstats.py` without changing user-visible behavior:

- header inference
- alias cleaning
- INFO / MAF / p-value / allele filters
- `N`, `N_cas`, `N_con`, `NSTUDY` handling
- merge-alleles behavior
- `.sumstats.gz` schema

After tests pass, rewrite `munge_sumstats.py` as a thin CLI wrapper that builds configs and calls `SumstatsMunger`.

Validation for this milestone:

- parity tests using the copied `munge_test` fixtures
- tests for alias inference and capitalization flexibility
- tests for merge-alleles success/failure
- tests for output schema equality with the legacy script

### Milestone 7: Regression module

Create:

- `ldscore/regression_workflow.py`

Define and implement:

- `RegressionDataset`
- `RegressionRunner`

Reuse the numerical kernels from:

- `ldscore.regressions`
- `ldscore.irwls`
- `ldscore.jackknife`

Keep `ldscore/sumstats.py` as a legacy workflow shim until cutover is complete. The new regression module must:

- build `RegressionDataset` only from already aggregated cross-chromosome LD-score results
- preserve the legacy default of using `M_5_50`-style common-SNP counts in regression unless explicitly overridden by config
- remove zero-variance LD-score columns and expose dropped column names directly
- run h2, partitioned h2, rg, and h2-cts through `RegressionRunner`
- add the new required feature: batch partitioned h2 over multiple query annotations, one model per query annotation, then combine/save the summary table

The new batch partitioned-h2 feature is the only intentional feature addition in this milestone. Everything else should preserve current behavior and outputs.

After parity is proven, route the summary-statistics modes in `ldsc.py` and the old workflow assembly in `ldscore/sumstats.py` through `RegressionRunner`.

Validation for this milestone:

- parity tests for `--h2`
- parity tests for `--rg`
- parity tests for `--h2-cts`
- tests that default regression count selection is `M_5_50`
- tests for `--not-M-5-50`-style override behavior via config
- tests for batch partitioned h2 with multiple query annotations producing one combined summary table with one row per query annotation

### Milestone 8: CLI cutover, compatibility wrappers, and cleanup

Modify these entry points to become thin wrappers:

- `ldsc.py`
- `ldsc_new.py`
- `munge_sumstats.py`
- `make_annot.py`
- `utils/run_bed_to_annot.py`

Compatibility rules:

- keep the repository CLI-first
- preserve existing file-format outputs unless the design docs explicitly require a public change
- when legacy argument names conflict with the design docs, follow the design docs for the new internal config layer and keep wrapper-level translation only where it does not create ambiguity
- keep numerical kernel modules (`ldscore/parse.py`, `ldscore/ldscore.py`, `ldscore/regressions.py`, `ldscore/irwls.py`, `ldscore/jackknife.py`) stable unless a specific failing test requires a change

Finalize the restructure by updating:

- `architecture.md`
- `code_structure.md`
- `README.md` if CLI usage examples move
- `setup.py` only if entry-point wiring actually changes

Do not delete transitional wrappers or compatibility imports until the full validation matrix is green.

## Concrete Steps

Use the repository root for all commands:

    cd /Users/wenbinwu/Documents_local/Research/SullivanLab/LDSC/repos/ldsc_py3_Jerry_workspace/ldsc_py3_Jerry

Baseline and kernel validation:

    pytest tests/test_parse_legacy.py tests/test_irwls_legacy.py tests/test_jackknife_legacy.py tests/test_regressions_legacy.py -q
    # expected: all selected tests pass

Annotation, reference-panel, and LD-score validation:

    pytest tests/test_annotation_module.py tests/test_ref_panel.py tests/test_ldscore_workflow.py tests/test_cli_l2.py -q
    # expected: all selected tests pass

Summary-statistics and regression validation:

    pytest tests/test_munge_workflow.py tests/test_regression_workflow.py tests/test_cli_regression.py -q
    # expected: all selected tests pass

Full repository validation before declaring the restructure stable:

    pytest -q
    # expected: all tests pass

CLI smoke checks after cutover:

    python ldsc.py -h
    python ldsc_new.py -h
    python munge_sumstats.py -h
    python make_annot.py -h
    python utils/run_bed_to_annot.py -h
    # expected: help text is printed and commands exit successfully

## Validation

The restructure is complete only when all of the following are true:

- the new modules exist and the old CLI scripts are thin wrappers around them
- copied legacy fixtures live inside this repository under `tests/fixtures/legacy/`
- the numerical kernels still match the current behavior on the selected fixture suite
- `ldsc.py --l2`, `ldsc_new.py`, `munge_sumstats.py`, `ldsc.py --h2`, `ldsc.py --rg`, `ldsc.py --h2-cts`, `make_annot.py`, and `utils/run_bed_to_annot.py` all pass parity or schema-equivalence tests
- the regression default still uses `M_5_50`-style common-SNP counts
- the new batch partitioned-h2 workflow runs one model per query annotation and writes one combined summary table
- `architecture.md`, `code_structure.md`, and `class-and-features.md` all agree with the implemented structure

## Interfaces and Dependencies

New workflow-facing modules to create:

- `ldscore/config.py`
  Defines `CommonConfig`, `AnnotationBuildConfig`, `RefPanelConfig`, `LDScoreConfig`, `MungeConfig`, and `RegressionConfig`.
- `ldscore/identifiers.py`
  Defines flexible SNP/header inference, chromosome normalization, canonical SNP ID builders, and restriction-list readers.
- `ldscore/output.py`
  Defines `OutputSpec`, `ArtifactConfig`, `RunSummary`, `ArtifactProducer`, `OutputManager`, `ResultWriter`, `ResultFormatter`, and `PostProcessor`.
- `ldscore/annotation.py`
  Defines `AnnotationSourceSpec`, `AnnotationBundle`, and `AnnotationBuilder`.
- `ldscore/ref_panel.py`
  Defines `RefPanelSpec`, abstract `RefPanel`, `PlinkRefPanel`, `ParquetR2RefPanel`, and `RefPanelLoader`.
- `ldscore/ldscore_workflow.py`
  Defines `ChromLDScoreResult`, `LDScoreResult`, and `LDScoreCalculator`.
- `ldscore/sumstats_munger.py`
  Defines `RawSumstatsSpec`, `SumstatsTable`, `MungeRunSummary`, and `SumstatsMunger`.
- `ldscore/regression_workflow.py`
  Defines `RegressionDataset` and `RegressionRunner`.

Existing modules to keep and reuse rather than rewrite from scratch:

- `ldscore/parse.py`
- `ldscore/ldscore.py`
- `ldscore/regressions.py`
- `ldscore/irwls.py`
- `ldscore/jackknife.py`

Existing modules to convert into thin wrappers or compatibility layers:

- `ldscore/ldscore_new.py`
- `ldscore/sumstats.py`
- `ldsc.py`
- `ldsc_new.py`
- `munge_sumstats.py`
- `make_annot.py`
- `utils/run_bed_to_annot.py`

Test layout to create:

- `tests/fixtures/legacy/`
- `tests/test_parse_legacy.py`
- `tests/test_irwls_legacy.py`
- `tests/test_jackknife_legacy.py`
- `tests/test_regressions_legacy.py`
- `tests/test_annotation_module.py`
- `tests/test_ref_panel.py`
- `tests/test_ldscore_workflow.py`
- `tests/test_munge_workflow.py`
- `tests/test_regression_workflow.py`
- `tests/test_cli_l2.py`
- `tests/test_cli_regression.py`

## Decision Log

- Decision: Use one master restructure plan with module-sized milestones instead of many smaller unrelated plans.
  Rationale: The work has strong dependencies across modules, but each milestone can still end in a testable, reviewable state.

- Decision: Add new modules alongside the old code first, then cut over the CLI wrappers after parity tests pass.
  Rationale: This keeps the refactor manageable and avoids breaking every workflow at once.

- Decision: Treat `architecture.md`, `code_structure.md`, and `class-and-features.md` as the higher-priority source of truth for public interfaces and module boundaries.
  Rationale: The user explicitly wants the refactor to follow the design docs even when the current code differs.

- Decision: Preserve the original LDSC regression default of using `M_5_50`-style common-SNP counts unless explicitly overridden.
  Rationale: This is the documented default in the current implementation and was explicitly chosen for the refactor.

- Decision: Keep the output module extensible through `ArtifactProducer` registration rather than hardcoded branches in `OutputManager`.
  Rationale: The user wants future plots, summaries, and reports without making the output layer rigid.

- Decision: Add the new batch partitioned-h2 workflow in the regression milestone rather than earlier.
  Rationale: It depends on the aggregated LD-score result, the output module, and the final regression runner structure.

Plan created on 2026-04-14 to drive the full repository restructuring from the current mixed-script layout to the layered design captured in the Markdown design documents.

## Archive

<!-- Completed or abandoned plans are moved here with a brief note on outcome. -->
