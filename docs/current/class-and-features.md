# Classes And Features

This document summarizes the public package surface. For workflow-level file streams, see [data-flow.md](data-flow.md).

## Feature Inventory

| Feature | CLI | Python entry points | Main inputs | Main outputs |
| --- | --- | --- | --- | --- |
| Build query annotations | `ldsc annotate` | `AnnotationBuilder`, `run_bed_to_annot()`, `run_annotate_from_args()`, `annotation_builder.main()` | baseline `.annot(.gz)`, BED inputs | `query.<chrom>.annot.gz`; workflow wrappers also write `annotate.log` |
| Build parquet reference panels | `ldsc build-ref-panel` | `ReferencePanelBuilder`, `run_build_ref_panel()` | PLINK prefix, optional source build inferred from `.bim`, optional liftover chains, conditional genetic maps, optional keep/restrict files; restriction identifier read from `GlobalConfig` and coordinates interpreted in the source build | per-build `chr*_r2.parquet` and `chr*_meta.tsv.gz` artifacts; optional `dropped_snps/chr*_dropped.tsv.gz`; workflow wrappers also write `build-ref-panel.log` |
| Compute LD scores | `ldsc ldscore` | `LDScoreCalculator`, `run_ldscore()` | optional baseline annotation shards, optional query annotations only when baseline is explicit, PLINK or parquet reference panel, optional frequency metadata | `manifest.json`, `ldscore.baseline.parquet`, optional `ldscore.query.parquet` under `output_dir`; parquet row groups are chromosome-aligned; no-annotation runs write synthetic `base`; workflow wrappers also write `ldscore.log` |
| Infer `chr_pos` genome build | workflow flags only: `--genome-build auto`; no standalone CLI command | `infer_chr_pos_build()`, `resolve_genome_build()`, `resolve_chr_pos_table()` | pandas table with `CHR` and `POS`; optional reference table | `ChrPosBuildInference`, resolved `GlobalConfig`, and optionally a normalized 1-based table |
| Munge GWAS summary statistics | `ldsc munge-sumstats` | `SumstatsMunger`, `load_sumstats()` | raw sumstats via `--raw-sumstats-file` or `MungeConfig.raw_sumstats_file`, optional `--trait-name`, column hints, QC thresholds, optional `--chr`/`--pos`, `--genome-build auto`, optional `--daner-old`/`--daner-new` schema handling, optional `--sumstats-snps-file` keep-list, optional `--output-format parquet\|tsv.gz\|both` | `sumstats.parquet` by default, optional `sumstats.sumstats.gz`, `sumstats.log`, `sumstats.metadata.json` |
| Estimate heritability | `ldsc h2` | `RegressionRunner.estimate_h2()` | munged `sumstats.parquet` or `.sumstats.gz` plus sidecar when available, LD-score directory | `h2.tsv`; `h2.log` when `output_dir` is supplied |
| Estimate partitioned heritability | `ldsc partitioned-h2` | `RegressionRunner.estimate_partitioned_h2()`, `RegressionRunner.estimate_partitioned_h2_batch()`, `PartitionedH2DirectoryWriter` | munged `sumstats.parquet` or `.sumstats.gz` plus sidecar when available, LD-score directory with non-empty query LD scores | compact `partitioned_h2.tsv`; optional `query_annotations/manifest.tsv`, per-query `partitioned_h2.tsv`, `partitioned_h2_full.tsv`, and `metadata.json` with `--write-per-query-results`; `partitioned-h2.log` when `output_dir` is supplied |
| Estimate genetic correlation | `ldsc rg` | `RegressionRunner.estimate_rg()`, `RegressionRunner.estimate_rg_pairs()`, `RgDirectoryWriter` | two or more munged `sumstats.parquet` or `.sumstats.gz` files plus sidecars when available, optional `--anchor-trait`, LD-score directory | concise `rg.tsv`; full `rg_full.tsv`; `h2_per_trait.tsv`; optional `pairs/` detail tree; `rg.log` when `output_dir` is supplied |

## Workflow Logging

Artifact-writing workflow wrappers share a per-run logging contract. Logs are
preflighted with the deterministic scientific outputs, opened only after
preflight succeeds, and kept with a `Failed` footer when execution raises after
the log is open. `--log-level` controls ordinary package records; lifecycle
audit lines always appear in the file.

Workflow logs are audit artifacts, not data outputs. Result objects and
`output_paths` mappings exclude log paths, including
`MungeRunSummary.output_paths`. `sumstats.metadata.json["output_files"]`
records curated sumstats artifacts only.

For log filenames and API boundary details, see
[workflow-logging.md](workflow-logging.md).

## Public Classes And Result Types

### Configuration

| Type | Role |
| --- | --- |
| `GlobalConfig` | shared SNP identifier, genome-build, logging, and missing-metadata policy |
| `AnnotationBuildConfig` | annotation projection and bundle-building settings |
| `RefPanelConfig` | choose and parameterize a runtime reference-panel backend and its source paths |
| `ReferencePanelBuildConfig` | build a parquet reference panel from PLINK input |
| `LDScoreConfig` | LD-window settings, common-count threshold, and optional regression row-set restriction |
| `MungeConfig` | `raw_sumstats_file` input, optional `sumstats_snps_file` row restriction, column hints, DANER schema switches, munging thresholds, and output settings |
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
| `LDScoreDirectoryWriter` | write the canonical LD-score result directory, including chromosome-aligned parquet row groups and manifest row-group metadata |

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
| `RgResultFamily` | complete multi-trait genetic-correlation result family: concise rg table, full diagnostic table, per-trait h2 table, and per-pair metadata |
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
- `ldsc partitioned-h2` rejects baseline-only LD-score directories. Use
  `ldsc h2` or `ldsc rg` with synthetic `base` outputs, and use
  `partitioned-h2` only when `ldscore.query.parquet` and non-empty `query_columns`
  are present.
- Public outputs use fixed workflow filenames under `output_dir`; run identity comes from the directory name.
- LD-score `ldscore.baseline.parquet` and `ldscore.query.parquet` stay single flat files, but
  each parquet row group contains one chromosome. `manifest.json` records
  `row_group_layout`, `baseline_row_groups`, and `query_row_groups`.
- Missing output directories are created and existing directories are reused.
  For `munge-sumstats`, `ldscore`, `partitioned-h2`, and `annotate`, existing
  owned siblings from the workflow artifact family fail before writing starts
  unless the caller passes `--overwrite` or `overwrite=True`.
- With overwrite enabled, successful runs remove stale owned siblings that the
  current configuration did not produce. Unrelated files in the output
  directory are preserved.
- `build-ref-panel` keeps an expert-oriented exception: overwrite permits
  replacement of current candidate artifacts, but it does not remove stale
  optional target-build or `dropped_snps` siblings from earlier configurations.
- Raw user-authored inputs use permissive alias resolution through `column_inference.py`.
- Raw sumstats may begin with `##` metadata/comment lines; these are skipped before header inference, so `#CHROM` remains available as the chromosome header.
- `chr_pos` workflows that interpret external coordinates require an explicit genome build or `--genome-build auto`; workflows such as `build-ref-panel` may ignore `GlobalConfig.genome_build` when they own a separate source-build contract. `rsid` workflows do not use genome-build metadata.
- Package-written sumstats artifacts include canonical `CHR` and `POS` columns. They are populated from inferred or explicitly flagged raw columns and filled as missing when the raw file lacks coordinates.
- Package-written artifacts use stricter internal headers so downstream workflows reload them deterministically.

## Public Import Boundary

Stable user-facing imports are re-exported from `ldsc.__init__`. That includes the workflow services, config dataclasses, reference-panel abstractions, and convenience helpers such as `run_bed_to_annot()`, `run_ldscore()`, and `load_sumstats()`. Annotation parser helpers live in `ldsc.annotation_builder`; `annotation_builder.main()` is the supported parser entry point for BED-to-annotation projection. Internal modules under `ldsc._kernel` are implementation details and may change without the same compatibility promise.

The top-level package also re-exports `ConfigMismatchError` and
`validate_config_compatibility()` for notebook and library code that wants to
surface or preflight config compatibility explicitly.

Genome-build inference for `chr_pos` inputs is also part of the top-level Python
API: import `ChrPosBuildInference`, `infer_chr_pos_build()`, and
`resolve_chr_pos_table()` from `ldsc`. The command-line API keeps inference
inside existing workflow flags such as `ldsc annotate --genome-build auto` and
`ldsc ldscore --genome-build auto`; it intentionally does not add a standalone
`infer-build` subcommand.
