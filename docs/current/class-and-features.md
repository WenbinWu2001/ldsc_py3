# Classes And Features

This document summarizes the public package surface. For workflow-level file streams, see [data-flow.md](data-flow.md).

## Feature Inventory

| Feature | CLI | Python entry points | Main inputs | Main outputs |
| --- | --- | --- | --- | --- |
| Build query annotations | `ldsc annotate` | `AnnotationBuilder`, `run_bed_to_annot()`, `run_annotate_from_args()`, `annotation_builder.main()` | baseline `.annot(.gz)`, BED inputs; duplicate effective-key rows are dropped before BED projection | root `query.<chrom>.annot.gz`; diagnostics under `diagnostics/` include `metadata.json`, `dropped_snps/dropped.tsv.gz`, and `annotate.log` |
| Build parquet reference panels | `ldsc build-ref-panel` | `ReferencePanelBuilder`, `run_build_ref_panel()` | PLINK prefix, optional source build inferred from `.bim`, optional liftover chains or HM3 quick liftover in `chr_pos`-family modes, conditional genetic maps, optional keep/restrict files or `--use-hm3-snps`; restriction identifier read from `GlobalConfig` and coordinates interpreted in the source build; restriction files are identity-only filters with duplicate keys collapsed; duplicate coordinate groups drop-all in `chr_pos`-family modes | root per-build `chr*_r2.parquet` and `chr*_meta.tsv.gz` artifacts; diagnostics under `diagnostics/` include `metadata.json`, `dropped_snps/chr*_dropped.tsv.gz`, and build logs |
| Compute LD scores | `ldsc ldscore` | `LDScoreCalculator`, `run_ldscore()` | optional baseline annotation shards, optional query annotations only when baseline is explicit, PLINK or parquet reference panel, optional frequency metadata; duplicate frequency metadata identity clusters are dropped before filling `CM`/`MAF` | root `metadata.json`, `ldscore.baseline.parquet`, optional `ldscore.query.parquet` under `output_dir`; parquet row groups are chromosome-aligned; no-annotation runs write synthetic `base`; workflow wrappers also write `diagnostics/ldscore.log` |
| Infer `chr_pos` genome build | workflow flags only: `--genome-build auto`; no standalone CLI command | `infer_chr_pos_build()`, `resolve_genome_build()`, `resolve_chr_pos_table()` | pandas table with `CHR` and `POS`; optional reference table | `ChrPosBuildInference`, resolved `GlobalConfig`, and optionally a normalized 1-based table |
| Munge GWAS summary statistics | `ldsc munge-sumstats` | `SumstatsMunger`, `infer_raw_sumstats()`, `load_sumstats()` | raw sumstats via `--raw-sumstats-file` or `MungeConfig.raw_sumstats_file`, optional `--trait-name`, default `--format auto`, optional `--infer-only`, column hints only when inference cannot decide safely, QC thresholds, optional `--chr`/`--pos`, `--genome-build auto`, optional explicit DANER/PGC format profile, optional `--sumstats-snps-file` keep-list or `--use-hm3-snps`, optional `chr_pos`-family liftover via `--target-genome-build` plus `--liftover-chain-file` or `--use-hm3-snps --use-hm3-quick-liftover`, optional `--output-format parquet\|tsv.gz\|both` | root `metadata.json` and `sumstats.parquet` by default, optional `sumstats.sumstats.gz`; diagnostics under `diagnostics/` include `sumstats.log` and `dropped_snps/dropped.tsv.gz`; `--infer-only` writes nothing |
| Estimate heritability | `ldsc h2` | `RegressionRunner.estimate_h2()`, `H2DirectoryWriter` | munged `sumstats.parquet` or `.sumstats.gz` plus sidecar, LD-score directory | `h2.tsv` when `output_dir` is supplied; `diagnostics/metadata.json` and `diagnostics/h2.log` are provenance only; without `output_dir`, CLI prints compact TSV to stdout |
| Estimate partitioned heritability | `ldsc partitioned-h2` | `RegressionRunner.estimate_partitioned_h2()`, `RegressionRunner.estimate_partitioned_h2_batch()`, `PartitionedH2DirectoryWriter` | munged `sumstats.parquet` or `.sumstats.gz` plus sidecar, LD-score directory with non-empty query LD scores | compact `partitioned_h2.tsv`; optional `diagnostics/query_annotations/manifest.tsv`, per-query `partitioned_h2.tsv`, `partitioned_h2_full.tsv`, and `metadata.json` with `--write-per-query-results`; diagnostics include `partitioned-h2.log`; without `output_dir`, CLI prints compact TSV to stdout |
| Estimate genetic correlation | `ldsc rg` | `RegressionRunner.estimate_rg()`, `RegressionRunner.estimate_rg_pairs()`, `RgDirectoryWriter` | two or more munged `sumstats.parquet` or `.sumstats.gz` files plus sidecars, optional `--anchor-trait`, LD-score directory | concise `rg.tsv`; full `rg_full.tsv`; `h2_per_trait.tsv`; optional `diagnostics/pairs/` detail tree; diagnostics include `rg.log`; without `output_dir`, CLI prints compact `rg.tsv` to stdout |

## Workflow Logging

Artifact-writing workflow wrappers share a per-run logging contract. Logs are
preflighted with the deterministic scientific outputs, opened only after
preflight succeeds, and kept with a `Failed` footer when execution raises after
the log is open. `--log-level` controls ordinary package records; lifecycle
audit lines always appear in the file.

Workflow logs are audit artifacts, not data outputs. Result objects and
`output_paths` mappings exclude log paths. `MungeRunSummary.output_paths`
includes curated data artifacts and the dropped-SNP audit sidecar, but not
`diagnostics/sumstats.log`. Detailed sumstats provenance and output
bookkeeping are written to `diagnostics/sumstats.log`, row-level liftover drops
are written to `diagnostics/dropped_snps/dropped.tsv.gz`, and
`metadata.json` stays limited to the downstream contract.

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
| `MungeConfig` | `raw_sumstats_file` input, optional `sumstats_format`, optional INFO-list columns, optional `sumstats_snps_file` row restriction, optional `chr_pos` liftover target/method, column hints, legacy DANER schema switches, munging thresholds, and output settings |
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
| `LDScoreDirectoryWriter` | write the canonical LD-score result directory, including chromosome-aligned parquet row groups and root `metadata.json` row-group metadata |
| `H2DirectoryWriter` | write unpartitioned h2 result tables and diagnostic metadata sidecars |
| `PartitionedH2DirectoryWriter` | write compact and optional per-query partitioned-h2 result trees |
| `RgDirectoryWriter` | write the genetic-correlation result family and optional per-pair detail tree |

### Data And Result Objects

| Type | Role |
| --- | --- |
| `AnnotationBundle` | aligned SNP metadata plus baseline/query annotation matrices |
| `ReferencePanelBuildResult` | summary of parquet panel artifacts written by one build |
| `ChromLDScoreResult` | one chromosome’s LD-score and weight tables, plus `config_snapshot` provenance |
| `LDScoreResult` | aggregated cross-chromosome LD-score artifacts, plus `config_snapshot` provenance |
| `SumstatsTable` | validated LDSC-ready summary-statistics table with canonical `SNP`, `CHR`, `POS`, `Z`, and `N` when available, plus known or unknown `config_snapshot` provenance |
| `RawSumstatsInference` | header-level `munge-sumstats` inference report with detected format, safe column hints, missing fields, notes, and suggested command arguments |
| `MungeRunSummary` | compact record of a munging run |
| `RegressionDataset` | merged sumstats plus LD-score matrix used by the estimator, plus propagated provenance when available |
| `RgResultFamily` | complete multi-trait genetic-correlation result family: concise rg table, full diagnostic table, per-trait h2 table, and per-pair metadata |
| `ChrPosBuildInference` | genome-build and coordinate-basis decision returned by `infer_chr_pos_build()` and `resolve_chr_pos_table()` |

### Global Config Registry

| Helper | Role |
| --- | --- |
| `get_global_config()` | return the package-global Python workflow configuration |
| `set_global_config()` | replace the package-global Python workflow configuration |
| `reset_global_config()` | restore the default package-global configuration: `GlobalConfig(snp_identifier="chr_pos_allele_aware", genome_build="auto")` |
| `validate_config_compatibility()` | compare two `GlobalConfig` snapshots and raise or warn on mismatch |

## Public Path And Header Contracts

- Public inputs accept exact paths, standard globs, and explicit chromosome suites using `@`.
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
  each parquet row group contains one chromosome. `metadata.json` records
  `row_group_layout`, `baseline_row_groups`, and `query_row_groups`.
- Missing output directories are created and existing directories are reused.
  For `munge-sumstats`, `build-ref-panel`, `ldscore`, `partitioned-h2`, `rg`,
  and `annotate`, existing owned siblings from the workflow artifact family
  fail before writing starts unless the caller passes `--overwrite` or
  `overwrite=True`.
- With overwrite enabled, successful runs remove stale owned siblings that the
  current configuration did not produce. Unrelated files in the output
  directory are preserved.
- Raw user-authored inputs use permissive alias resolution through `column_inference.py`.
- Raw sumstats may begin with `##` metadata/comment lines; these are skipped before header inference, so `#CHROM` remains available as the chromosome header.
- `chr_pos` workflows that interpret external coordinates require an explicit genome build or `--genome-build auto`; workflows such as `build-ref-panel` may ignore `GlobalConfig.genome_build` when they own a separate source-build contract. `rsid` workflows do not use genome-build metadata.
- Public SNP identifier modes are exactly `rsid`, `rsid_allele_aware`,
  `chr_pos`, and `chr_pos_allele_aware`; the default is
  `chr_pos_allele_aware`. Mode names are exact. Column aliases apply to input
  headers only.
- Base modes are allele-blind: `rsid` uses only `SNP`, and `chr_pos` uses only
  `CHR:POS`. Allele columns may be preserved as passive data but never affect
  base-mode identity, duplicate filtering, retention, or drop reasons.
- Allele-aware modes require usable `A1/A2` on sumstats, reference-panel
  artifacts, R2 parquet endpoints, and LD-score artifacts. They drop missing,
  invalid/non-SNP, identical, strand-ambiguous allele pairs, multi-allelic
  base-key clusters, and duplicate effective-key clusters.
- Restriction files may omit alleles and then match by base key. Allele-bearing
  restrictions, including packaged HM3 restrictions, match by the effective
  allele-aware key in allele-aware modes. Annotation files may omit alleles in
  allele-aware modes because annotations describe genomic membership; if
  annotations include alleles, they participate in allele-aware
  matching.
- External raw R2 parquet inputs are supported only in `rsid` and `chr_pos`.
  Allele-aware modes require package-built canonical R2 parquet with endpoint
  allele columns `A1_1/A2_1/A1_2/A2_2`.
- `--allow-identity-downgrade` is regression-only. It allows same-family
  allele-aware/base mixes to run under the base mode; rsID-family and
  coordinate-family modes never mix.
- Package-written sumstats artifacts include canonical `CHR` and `POS` columns. They are populated from inferred or explicitly flagged raw columns and filled as missing when the raw file lacks coordinates.
- In `chr_pos`-family modes, `SNP` is a label and row identity is based on `CHR/POS`; munger liftover updates only `CHR/POS` and is rejected in `rsid`-family modes.
- `build-ref-panel` chain liftover is also coordinate behavior: matching chains are rejected in `rsid`-family modes, and source-only rsID-family builds skip coordinate duplicate filtering.
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
