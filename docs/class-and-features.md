# Classes And Features

This document summarizes the public package surface. For workflow-level file streams, see [data_flow.md](data_flow.md).

## Feature Inventory

| Feature | CLI | Python entry points | Main inputs | Main outputs |
| --- | --- | --- | --- | --- |
| Build query annotations | `ldsc annotate` | `AnnotationBuilder`, `run_bed_to_annot()`, `make_annot_files()` | baseline `.annot(.gz)`, BED or gene-set inputs | `query.<chrom>.annot.gz` or one legacy `.annot(.gz)` |
| Build parquet reference panels | `ldsc build-ref-panel` | `ReferencePanelBuilder`, `run_build_ref_panel()` | PLINK prefix, optional liftover chains, conditional genetic maps, optional keep/restrict files; `snp_identifier` only when a SNP restriction file is supplied | per-chromosome `ann.parquet`, `LD.parquet`, emitted `meta_*.tsv.gz` sidecars |
| Compute LD scores | `ldsc ldscore` | `LDScoreCalculator`, `run_ldscore()` | optional baseline annotation shards, optional query annotations only when baseline is explicit, PLINK or parquet reference panel, optional frequency metadata | `manifest.json`, `baseline.parquet`, optional `query.parquet` under `output_dir`; no-annotation runs write synthetic `base` |
| Infer `chr_pos` genome build | workflow flags only: `--genome-build auto`; no standalone CLI command | `infer_chr_pos_build()`, `resolve_genome_build()`, `resolve_chr_pos_table()` | pandas table with `CHR` and `POS`; optional reference table | `ChrPosBuildInference`, resolved `GlobalConfig`, and optionally a normalized 1-based table |
| Munge GWAS summary statistics | `ldsc munge-sumstats` | `SumstatsMunger`, `load_sumstats()` | raw sumstats, column hints, QC thresholds, optional `--chr`/`--pos`, `--genome-build auto` | `sumstats.sumstats.gz`, `sumstats.log`, `sumstats.metadata.json` |
| Estimate heritability | `ldsc h2` | `RegressionRunner.estimate_h2()` | munged `.sumstats.gz` plus sidecar when available, LD-score directory | `h2.tsv` |
| Estimate partitioned heritability | `ldsc partitioned-h2` | `RegressionRunner.estimate_partitioned_h2_batch()` | munged `.sumstats.gz` plus sidecar when available, LD-score directory | `partitioned_h2.tsv` |
| Estimate genetic correlation | `ldsc rg` | `RegressionRunner.estimate_rg()` | two munged `.sumstats.gz` files plus sidecars when available, LD-score directory | `rg.tsv` |

## Public Classes And Result Types

### Configuration

| Type | Role |
| --- | --- |
| `GlobalConfig` | shared SNP identifier, genome-build, logging, and missing-metadata policy |
| `AnnotationBuildConfig` | annotation projection and bundle-building settings |
| `RefPanelConfig` | choose and parameterize a runtime reference-panel backend and its source paths |
| `ReferencePanelBuildConfig` | build a parquet reference panel from PLINK input |
| `LDScoreConfig` | LD-window and retained-SNP settings |
| `MungeConfig` | raw-sumstats input, column hints, munging thresholds, and output settings |
| `RegressionConfig` | regression-model settings such as intercept handling and jackknife blocks |
| `ConfigMismatchError` | explicit failure raised when critical config assumptions disagree |

### Workflow Services

| Type | Role |
| --- | --- |
| `AnnotationBuilder` | load aligned annotation bundles or project BED inputs to SNP-level annotations |
| `ReferencePanelBuilder` | emit standard parquet reference artifacts from PLINK |
| `RefPanelLoader` | load runtime PLINK or parquet reference-panel adapters |
| `LDScoreCalculator` | run per-chromosome LD-score computation and aggregate outputs |
| `SumstatsMunger` | normalize raw GWAS tables into curated LDSC-ready tables |
| `RegressionRunner` | build regression datasets and run `h2`, partitioned `h2`, and `rg` |
| `LDScoreDirectoryWriter` | write the canonical LD-score result directory |

### Data And Result Objects

| Type | Role |
| --- | --- |
| `AnnotationBundle` | aligned SNP metadata plus baseline/query annotation matrices |
| `ReferencePanelBuildResult` | summary of parquet panel artifacts written by one build |
| `ChromLDScoreResult` | one chromosome’s LD-score and weight tables, plus `config_snapshot` provenance |
| `LDScoreResult` | aggregated cross-chromosome LD-score artifacts, plus `config_snapshot` provenance |
| `SumstatsTable` | validated LDSC-ready summary-statistics table with canonical `SNP`, `CHR`, `POS`, `Z`, and `N` when available, plus known or unknown `config_snapshot` provenance |
| `MungeRunSummary` | compact record of a munging run |
| `RegressionDataset` | merged sumstats plus LD-score matrix used by the estimator, plus propagated provenance when available |
| `ChrPosBuildInference` | genome-build and coordinate-basis decision returned by `infer_chr_pos_build()` and `resolve_chr_pos_table()` |

### Global Config Registry

| Helper | Role |
| --- | --- |
| `get_global_config()` | return the package-global Python workflow configuration |
| `set_global_config()` | replace the package-global Python workflow configuration |
| `reset_global_config()` | restore the default package-global configuration: `GlobalConfig(snp_identifier="rsid")` |
| `validate_config_compatibility()` | compare two `GlobalConfig` snapshots and raise or warn on mismatch |

## Public Path And Header Contracts

- Public inputs accept exact paths, standard globs, explicit chromosome suites using `@`, and some legacy bare prefixes.
- Public dataclasses normalize `PathLike` objects to strings but do not expand inputs immediately.
- Workflow modules resolve input tokens before calling `_kernel`.
- `ldsc ldscore` accepts no baseline/query inputs for ordinary unpartitioned LD-score generation; the workflow creates a synthetic all-ones baseline column named exactly `base` from retained reference-panel metadata.
- `query_annot_sources` and `query_annot_bed_sources` require explicit `baseline_annot_sources`; users who want to test query annotations against an all-ones universe must materialize that `base` baseline annotation themselves and run the partitioned workflow.
- Public outputs use fixed workflow filenames under `output_dir`; run identity comes from the directory name.
- Missing output directories are created, existing directories are reused, and
  existing fixed files fail before writing starts unless the caller passes
  `--overwrite` or `overwrite=True`.
- `--overwrite` only permits replacement of the known files for the active
  workflow; it does not remove unrelated files from the output directory.
- Raw user-authored inputs use permissive alias resolution through `column_inference.py`.
- Raw sumstats may begin with `##` metadata/comment lines; these are skipped before header inference, so `#CHROM` remains available as the chromosome header.
- `chr_pos` workflows require an explicit genome build or `--genome-build auto`; `rsid` workflows do not use genome-build metadata.
- Package-written sumstats artifacts include canonical `CHR` and `POS` columns. They are populated from inferred or explicitly flagged raw columns and filled as missing when the raw file lacks coordinates.
- Package-written artifacts use stricter internal headers so downstream workflows reload them deterministically.

## Public Import Boundary

Stable user-facing imports are re-exported from `ldsc.__init__`. That includes the workflow services, config dataclasses, reference-panel abstractions, and convenience helpers such as `run_ldscore()` and `load_sumstats()`. Internal modules under `ldsc._kernel` are implementation details and may change without the same compatibility promise.

The top-level package also re-exports `ConfigMismatchError` and
`validate_config_compatibility()` for notebook and library code that wants to
surface or preflight config compatibility explicitly.

Genome-build inference for `chr_pos` inputs is also part of the top-level Python
API: import `ChrPosBuildInference`, `infer_chr_pos_build()`, and
`resolve_chr_pos_table()` from `ldsc`. The command-line API keeps inference
inside existing workflow flags such as `ldsc annotate --genome-build auto` and
`ldsc ldscore --genome-build auto`; it intentionally does not add a standalone
`infer-build` subcommand.
