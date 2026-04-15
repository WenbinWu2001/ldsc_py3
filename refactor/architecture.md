# Architecture

`refactor/` is the active Python 3 LDSC package. Its package root is `src/ldsc/`. The codebase is organized around one public workflow module per user-visible function, backed by an internal `_kernel/` layer that contains reused numerical methods and low-level file-format code.

The package reads SNP-level annotations, PLINK or parquet-R2 reference inputs, and munged GWAS summary statistics; computes LD scores chromosome by chromosome; aggregates those results across chromosomes; then runs LDSC regressions on the combined tables. Output writing is routed through a separate artifact layer so the core workflows stay independent of file layout and future post-processing features.

## Package Shape

```text
refactor/
├── setup.py
├── src/
│   └── ldsc/
│       ├── __init__.py
│       ├── cli.py
│       ├── config.py
│       ├── outputs.py
│       ├── annotation_builder.py
│       ├── ldscore_calculator.py
│       ├── sumstats_munger.py
│       ├── regression_runner.py
│       └── _kernel/
│           ├── __init__.py
│           ├── annotation.py
│           ├── formats.py
│           ├── identifiers.py
│           ├── ref_panel.py
│           ├── ldscore.py
│           ├── sumstats_munger.py
│           ├── regression.py
│           ├── _jackknife.py
│           └── _irwls.py
├── tests/
├── tutorials/
├── architecture.md
├── class-and-features.md
├── code_structure.md
├── AGENTS.md
└── PLANS.md
```

## Layer Model

| Layer | Location | Responsibility |
| --- | --- | --- |
| Public CLI | `src/ldsc/cli.py` | single command surface and subcommand dispatch |
| Public workflow modules | `src/ldsc/*.py` | validated user-facing services, dataclasses, and file-driven entry helpers |
| Internal kernel | `src/ldsc/_kernel/*.py` | numerical kernels, low-level readers, and compatibility logic |
| Tests | `tests/` | parity checks, API smoke tests, and workflow coverage |
| Tutorials and docs | `tutorials/`, `*.md` | usage examples and project guidance |

## Public Surface

The public package exports the feature-level services directly from `ldsc`:

- `AnnotationBuilder`
- `LDScoreCalculator`
- `SumstatsMunger`
- `RegressionRunner`
- `OutputManager`
- configuration dataclasses
- reference-panel abstractions

The public CLI is also unified:

- `ldsc annotate`
- `ldsc ldscore`
- `ldsc munge-sumstats`
- `ldsc h2`
- `ldsc partitioned-h2`
- `ldsc rg`

This replaces the previous split between root-level wrapper scripts.

## Feature Mapping

| Feature | Public module | Internal kernel |
| --- | --- | --- |
| SNP-level annotation building | `ldsc.annotation_builder` | `ldsc._kernel.annotation` |
| LD-score calculation | `ldsc.ldscore_calculator` | `ldsc._kernel.ldscore` |
| Summary-statistics munging | `ldsc.sumstats_munger` | `ldsc._kernel.sumstats_munger` |
| Regression assembly and execution | `ldsc.regression_runner` | `ldsc._kernel.regression`, `ldsc._kernel._jackknife`, `ldsc._kernel._irwls` |
| Output writing and post-processing hooks | `ldsc.outputs` | none; this layer is already public and standalone |
| Reference-panel access | public types re-exported from `ldsc` | `ldsc._kernel.ref_panel`, `ldsc._kernel.formats`, `ldsc._kernel.identifiers` |

## Data Flow

### Annotation and LD-score path

1. `ldsc annotate` or `AnnotationBuilder` resolves `.annot`, BED, or gene-set inputs.
2. `ldsc.ldscore_calculator` loads and validates chromosome-specific annotation groups.
3. `ldsc._kernel.ref_panel` opens PLINK or parquet-R2 reference metadata.
4. `ldsc._kernel.ldscore` computes one chromosome at a time.
5. `LDScoreCalculator` aggregates chromosome results into one `LDScoreResult`.
6. `OutputManager` writes `.l2.ldscore`, `.w.l2.ldscore`, `.M`, `.M_5_50`, and manifest artifacts.

### Munging and regression path

1. `ldsc munge-sumstats` or `SumstatsMunger` normalizes raw summary statistics into LDSC-compatible `.sumstats.gz`.
2. `RegressionRunner` loads the munged sumstats, combined LD-score tables, weight tables, and count vector.
3. `RegressionRunner` builds a `RegressionDataset`, drops zero-variance LD-score columns, and selects the count vector used for regression.
4. `ldsc._kernel.regression` runs `Hsq`, `RG`, and related estimator code.
5. Results are formatted as compact summary tables.

## Architectural Invariants

- `src/ldsc/` is the only public Python import surface.
- `src/ldsc/_kernel/` is internal-only and may change without user-facing API guarantees.
- `refactor/` must not import from `legacy/`.
- LD-score calculation remains chromosome-wise; regression consumes only the aggregated cross-chromosome result.
- Regression preserves the original LDSC default of using `.M_5_50` when available.
- Output extensibility goes through artifact producers rather than adding feature-specific branches throughout the workflows.

## Where To Change Code

| Goal | Start here |
| --- | --- |
| add or adjust CLI options | `src/ldsc/cli.py` |
| change annotation loading or BED projection | `src/ldsc/annotation_builder.py`, then `src/ldsc/_kernel/annotation.py` |
| change LD-score orchestration or result objects | `src/ldsc/ldscore_calculator.py` |
| change PLINK/parquet LD-score math | `src/ldsc/_kernel/ldscore.py` |
| change munging API or wrapper behavior | `src/ldsc/sumstats_munger.py` |
| change regression dataset assembly or batch partitioned-h2 behavior | `src/ldsc/regression_runner.py` |
| change estimator math | `src/ldsc/_kernel/regression.py` |
| add new output files or future plot/report hooks | `src/ldsc/outputs.py` |
