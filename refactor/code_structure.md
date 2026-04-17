# Code Structure

This document is a coding-oriented map of the current `refactor/` package after the `src/ldsc/` migration.

## Repository Map

```text
refactor/
├── setup.py
├── src/ldsc/
│   ├── __init__.py
│   ├── cli.py
│   ├── config.py
│   ├── path_resolution.py
│   ├── outputs.py
│   ├── annotation_builder.py
│   ├── ldscore_calculator.py
│   ├── sumstats_munger.py
│   ├── regression_runner.py
│   └── _kernel/
│       ├── annotation.py
│       ├── formats.py
│       ├── identifiers.py
│       ├── ref_panel.py
│       ├── ldscore.py
│       ├── sumstats_munger.py
│       ├── regression.py
│       ├── _jackknife.py
│       └── _irwls.py
├── tests/
├── tutorials/
├── architecture.md
├── class-and-features.md
├── code_structure.md
├── AGENTS.md
└── PLANS.md
```

## Layer Structure

| Layer | Files | Responsibility |
| --- | --- | --- |
| CLI | `src/ldsc/cli.py` | top-level command dispatch |
| Public workflow/config | `src/ldsc/config.py`, `path_resolution.py`, `annotation_builder.py`, `ldscore_calculator.py`, `sumstats_munger.py`, `regression_runner.py`, `outputs.py` | user-visible services, dataclasses, shared path resolution, and output config |
| Internal kernel | `src/ldsc/_kernel/*.py` | low-level readers, numerical kernels, compatibility helpers |
| Verification | `tests/` | parity and API checks |
| User guidance | `tutorials/`, `README.md`, design docs | usage examples and maintenance guidance |

## Feature Inventory

| Feature | Public entry point | Public module | Internal kernel | Main outputs |
| --- | --- | --- | --- | --- |
| BED or gene-set to `.annot` | `ldsc annotate` | `ldsc.annotation_builder` | `ldsc._kernel.annotation` | `.annot`, `.annot.gz` |
| LD-score calculation from PLINK | `ldsc ldscore` | `ldsc.ldscore_calculator` | `ldsc._kernel.ldscore` | `.l2.ldscore(.gz)`, `.w.l2.ldscore(.gz)`, `.M`, optional `.M_5_50` |
| LD-score calculation from parquet R2 | `ldsc ldscore` | `ldsc.ldscore_calculator` | `ldsc._kernel.ldscore` | same as above |
| Summary-statistics munging | `ldsc munge-sumstats` | `ldsc.sumstats_munger` | `ldsc._kernel.sumstats_munger` | `.sumstats.gz`, `.log` |
| Heritability | `ldsc h2` | `ldsc.regression_runner` | `ldsc._kernel.regression` | summary TSV |
| Partitioned heritability | `ldsc partitioned-h2` | `ldsc.regression_runner` | `ldsc._kernel.regression` | summary TSV with one row per query annotation |
| Genetic correlation | `ldsc rg` | `ldsc.regression_runner` | `ldsc._kernel.regression` | summary TSV |
| Output writing / future post-processing | `OutputManager` | `ldsc.outputs` | n/a | tables, manifests, JSON, metadata |

## CLI Surface

| Subcommand | Handler |
| --- | --- |
| `ldsc annotate` | `ldsc.cli._run_annotate()` |
| `ldsc ldscore` | `ldsc.ldscore_calculator.run_ldscore_from_args()` |
| `ldsc munge-sumstats` | `ldsc.sumstats_munger.main()` |
| `ldsc h2` | `ldsc.regression_runner.run_h2_from_args()` |
| `ldsc partitioned-h2` | `ldsc.regression_runner.run_partitioned_h2_from_args()` |
| `ldsc rg` | `ldsc.regression_runner.run_rg_from_args()` |

## Public Modules

### `src/ldsc/config.py`

Purpose:

- hold validated workflow configuration dataclasses

Main classes:

- `CommonConfig`
- `AnnotationBuildConfig`
- `RefPanelConfig`
- `LDScoreConfig`
- `MungeConfig`
- `RegressionConfig`

Path handling contract:

- public dataclasses normalize `PathLike` values to strings
- plural path-bearing fields accept either one token or a sequence of tokens
- tokens stay unresolved here; file discovery happens later in workflow modules

### `src/ldsc/path_resolution.py`

Purpose:

- define the shared public input-token language for filesystem inputs
- expand `~` and environment variables
- resolve exact paths, glob patterns, explicit `@` chromosome suites, and legacy bare prefixes
- keep output destinations literal by separating normalization from discovery

Main helpers:

- `normalize_path_token()`
- `normalize_path_tokens()`
- `split_cli_path_tokens()`
- `resolve_scalar_path()`
- `resolve_file_group()`
- `resolve_chromosome_group()`
- `resolve_plink_prefix()`

### `src/ldsc/annotation_builder.py`

Purpose:

- expose the annotation-building public surface without forcing users into `_kernel`

Public exports:

- `AnnotationSourceSpec`
- `AnnotationBundle`
- `AnnotationBuilder`
- `gene_set_to_bed`
- `make_annot_files`
- `run_bed_to_annot`
- CLI-oriented parser helpers

The actual parsing and projection logic lives in `ldsc._kernel.annotation`.

### `src/ldsc/ldscore_calculator.py`

Purpose:

- define workflow-layer LD-score result objects
- orchestrate chromosome-wise calculation
- aggregate outputs and route them through `OutputManager`

Main classes:

- `_LegacyChromResult`
- `ChromLDScoreResult`
- `LDScoreResult`
- `LDScoreCalculator`

Important functions:

- `build_parser()`
- `run_ldscore_from_args()`
- `run_ldscore(**kwargs)`
- `main(argv=None)`

Input streamline:

1. Parse CLI or keyword arguments.
2. Normalize identifier mode and load optional regression SNP restriction.
3. Resolve public annotation, PLINK, parquet, and frequency tokens into concrete chromosome-specific paths.
4. Build one annotation bundle per chromosome.
5. Dispatch to `_kernel.ldscore.compute_chrom_from_plink()` or `compute_chrom_from_parquet()` with primitive string paths only.
6. Wrap each chromosome as `ChromLDScoreResult`.
7. Aggregate to `LDScoreResult`.
8. Write outputs if requested.

Alignment guarantees:

- `reference_metadata` and `ld_scores` share row order.
- `regression_metadata` and `w_ld` share row order.
- `snp_count_totals` are summed across chromosomes by named key.

### `src/ldsc/sumstats_munger.py`

Purpose:

- expose a workflow-friendly wrapper over the legacy-style munging code

Main classes:

- `RawSumstatsSpec`
- `SumstatsTable`
- `MungeRunSummary`
- `SumstatsMunger`

Important exports:

- `parser`
- legacy header/QC helper functions re-exported from the kernel for compatibility

### `src/ldsc/regression_runner.py`

Purpose:

- assemble regression datasets from combined LD-score outputs and munged sumstats
- expose the user-facing h2, partitioned-h2, and rg entry helpers

Main classes:

- `RegressionDataset`
- `RegressionRunner`

Important functions:

- `add_h2_arguments()`
- `add_partitioned_h2_arguments()`
- `add_rg_arguments()`
- `run_h2_from_args()`
- `run_partitioned_h2_from_args()`
- `run_rg_from_args()`

Current regression data contract:

- inputs are already aggregated across chromosomes
- scalar file inputs use exact-one resolution before loading
- one count vector is selected for regression, defaulting to `.M_5_50`
- zero-variance LD-score columns are dropped before fitting
- batch partitioned-h2 loops over query annotations one at a time and concatenates the summary rows

### `src/ldsc/outputs.py`

Purpose:

- define output config and writing behavior
- keep output extensibility separate from scientific workflows

Main classes:

- `OutputSpec`
- `ArtifactConfig`
- `RunSummary`
- `Artifact`
- `ArtifactProducer`
- `ResultFormatter`
- `ResultWriter`
- `OutputManager`
- `PostProcessor`

Built-in artifact producers:

- `LDScoreTableProducer`
- `WeightLDProducer`
- `CountProducer`
- `AnnotationManifestProducer`
- `SummaryTSVProducer`
- JSON summary and run-metadata producers

## Internal Kernel Modules

### `src/ldsc/_kernel/annotation.py`

Purpose:

- parse `.annot` tables
- enforce row alignment
- apply global SNP restriction
- project BED inputs onto SNP rows

Defines the actual `AnnotationBuilder`, `AnnotationBundle`, and `AnnotationSourceSpec` implementations used by the public wrapper.
It also accepts workflow-resolved input tokens for annotation families and performs chromosome-shard bundling after resolution.

### `src/ldsc/_kernel/formats.py`

Purpose:

- low-level file-family readers and small format utilities reused by kernels

### `src/ldsc/_kernel/identifiers.py`

Purpose:

- normalize SNP identifier mode
- infer flexible SNP headers
- build canonical SNP identifier series
- read global SNP restrictions

### `src/ldsc/_kernel/ref_panel.py`

Purpose:

- reference-panel abstractions and backends

Main classes:

- `RefPanelSpec`
- `RefPanel`
- `PlinkRefPanel`
- `ParquetR2RefPanel`
- `RefPanelLoader`

### `src/ldsc/_kernel/ldscore.py`

Purpose:

- main compute kernel for LD-score calculation

This module consolidates the reused parts of the older PLINK genotype kernel and the newer chromosome/parquet engine.

Key contents:

- annotation combination helpers
- PLINK readers and block LD logic
- parquet R2 block readers and per-chromosome computation
- sorting and aggregation helpers used by the workflow layer

### `src/ldsc/_kernel/sumstats_munger.py`

Purpose:

- raw summary-stat parsing, header inference, QC filtering, allele checks, and file writing reused by the public munger

### `src/ldsc/_kernel/regression.py`

Purpose:

- actual regression estimator code and reusable allele-alignment helpers

Main classes:

- `Jackknife`
- `LstsqJackknifeFast`
- `LstsqJackknifeSlow`
- `RatioJackknife`
- `IRWLS`
- `LD_Score_Regression`
- `Hsq`
- `Gencov`
- `RG`

Helper modules:

- `_jackknife.py`
- `_irwls.py`

## Test Layout

| Test file | Coverage |
| --- | --- |
| `tests/test_package_layout.py` | public package and CLI subcommand shape |
| `tests/test_annotation.py` | annotation building and BED projection wrappers |
| `tests/test_ref_panel.py` | reference-panel abstractions |
| `tests/test_ldscore_workflow.py` | LD-score workflow/result assembly |
| `tests/test_sumstats_munger.py` | public munging wrapper |
| `tests/test_regression_workflow.py` | regression dataset assembly and batch partitioned-h2 |
| `tests/test_path_resolution.py` | shared path token resolution and compatibility rules |
| `tests/*_legacy.py` | parity-oriented coverage for reused kernel logic |

## Dead Code Status

The old root-level wrappers and old `refactor/ldscore/` package were removed. The active codebase no longer depends on sibling legacy modules. The remaining compatibility logic lives only inside `src/ldsc/_kernel/`.
