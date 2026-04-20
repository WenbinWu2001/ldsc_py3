# Classes And Features

This document summarizes the public package surface. For workflow-level file streams, see [data_flow.md](data_flow.md).

## Feature Inventory

| Feature | CLI | Python entry points | Main inputs | Main outputs |
| --- | --- | --- | --- | --- |
| Build query annotations | `ldsc annotate` | `AnnotationBuilder`, `run_bed_to_annot()`, `make_annot_files()` | baseline `.annot(.gz)`, BED or gene-set inputs | `query.<chrom>.annot.gz` or one legacy `.annot(.gz)` |
| Build parquet reference panels | `ldsc build-ref-panel` | `ReferencePanelBuilder`, `run_build_ref_panel()` | PLINK prefix, genetic maps, optional keep/restrict files | per-chromosome `ann.parquet`, `LD.parquet`, `meta_*.tsv.gz` |
| Compute LD scores | `ldsc ldscore` | `LDScoreCalculator`, `run_ldscore()` | annotation shards, PLINK or parquet reference panel, optional frequency metadata | `.l2.ldscore(.gz)`, `.w.l2.ldscore(.gz)`, `.l2.M`, `.l2.M_5_50`, manifests |
| Munge GWAS summary statistics | `ldsc munge-sumstats` | `SumstatsMunger`, `load_sumstats()` | raw sumstats, column hints, QC thresholds | `.sumstats.gz`, `.log` |
| Estimate heritability | `ldsc h2` | `RegressionRunner.estimate_h2()` | munged `.sumstats.gz`, LD-score files, count vector | `.h2.tsv` |
| Estimate partitioned heritability | `ldsc partitioned-h2` | `RegressionRunner.estimate_partitioned_h2_batch()` | munged `.sumstats.gz`, LD-score files, count vector, annotation manifest | `.partitioned_h2.tsv` |
| Estimate genetic correlation | `ldsc rg` | `RegressionRunner.estimate_rg()` | two munged `.sumstats.gz` files, LD-score files, count vector | `.rg.tsv` |

## Public Classes And Result Types

### Configuration

| Type | Role |
| --- | --- |
| `CommonConfig` | shared SNP identifier, genome-build, logging, and restriction settings |
| `AnnotationBuildConfig` | annotation projection and bundle-building settings |
| `RefPanelConfig` | choose and parameterize a runtime reference-panel backend |
| `ReferencePanelBuildConfig` | build a parquet reference panel from PLINK input |
| `LDScoreConfig` | LD-window and retained-SNP settings |
| `MungeConfig` | raw-sumstats munging thresholds and output settings |
| `RegressionConfig` | regression-model settings such as intercept handling and jackknife blocks |
| `OutputSpec` | LD-score artifact layout and emission controls |

### Workflow Services

| Type | Role |
| --- | --- |
| `AnnotationBuilder` | load aligned annotation bundles or project BED inputs to SNP-level annotations |
| `ReferencePanelBuilder` | emit standard parquet reference artifacts from PLINK |
| `RefPanelLoader` | load runtime PLINK or parquet reference-panel adapters |
| `LDScoreCalculator` | run per-chromosome LD-score computation and aggregate outputs |
| `SumstatsMunger` | normalize raw GWAS tables into curated LDSC-ready tables |
| `RegressionRunner` | build regression datasets and run `h2`, partitioned `h2`, and `rg` |
| `OutputManager` | write LD-score artifacts, manifests, and summaries |

### Data And Result Objects

| Type | Role |
| --- | --- |
| `AnnotationSourceSpec` | annotation input token bundle for one run |
| `AnnotationBundle` | aligned SNP metadata plus baseline/query annotation matrices |
| `RefPanelSpec` | runtime description of a PLINK or parquet reference panel |
| `ReferencePanelBuildResult` | summary of parquet panel artifacts written by one build |
| `ChromLDScoreResult` | one chromosome’s LD-score and weight tables |
| `LDScoreResult` | aggregated cross-chromosome LD-score artifacts |
| `RawSumstatsSpec` | one raw summary-statistics input plus optional column hints |
| `SumstatsTable` | validated LDSC-ready summary-statistics table |
| `MungeRunSummary` | compact record of a munging run |
| `RegressionDataset` | merged sumstats plus LD-score matrix used by the estimator |
| `RunSummary` | summary emitted by the output layer |

## Public Path And Header Contracts

- Public inputs accept exact paths, standard globs, explicit chromosome suites using `@`, and some legacy bare prefixes.
- Public dataclasses normalize `PathLike` objects to strings but do not expand inputs immediately.
- Workflow modules resolve input tokens before calling `_kernel`.
- Output destinations such as `out_prefix` are normalized but never glob-expanded.
- Raw user-authored inputs use permissive alias resolution through `column_inference.py`.
- Package-written artifacts use stricter internal headers so downstream workflows reload them deterministically.

## Public Import Boundary

Stable user-facing imports are re-exported from `ldsc.__init__`. That includes the workflow services, config dataclasses, reference-panel abstractions, and convenience helpers such as `run_ldscore()` and `load_sumstats()`. Internal modules under `ldsc._kernel` are implementation details and may change without the same compatibility promise.
