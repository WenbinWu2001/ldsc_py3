# Layer Structure

This document maps the refactored `ldsc` package by layer and by main
functionality. Use it to answer two questions quickly:

- What layer owns this object?
- Which objects participate in a given workflow?

The public import boundary is `ldsc` and the public command boundary is
`ldsc.cli`. Modules under `ldsc._kernel` are internal implementation details.

## Layer Matrix

| Layer | Shared CLI, config, path, and schema contracts | Annotation building | Reference-panel building and loading | LD-score calculation | Summary-statistics munging | Regression | Output and persistence |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Public API and command surface | `ldsc.__init__` re-exports stable objects; `ldsc.__main__.main()` runs `python -m ldsc`; `ldsc.cli.build_parser()` and `ldsc.cli.main()` define the `ldsc` subcommands. Top-level imports include `ChrPosBuildInference`, `infer_chr_pos_build()`, and `resolve_chr_pos_table()`; the CLI has no standalone build-inference subcommand. | CLI: `ldsc annotate`; imports: `AnnotationBuilder`, `AnnotationBundle`, `AnnotationBuildConfig`, `run_bed_to_annot()`. | CLI: `ldsc build-ref-panel`; imports: `ReferencePanelBuilder`, `ReferencePanelBuildResult`, `run_build_ref_panel()`, `RefPanel`, `RefPanelLoader`, `RefPanelConfig`, `PlinkRefPanel`, `ParquetR2RefPanel`. | CLI: `ldsc ldscore`; imports: `LDScoreCalculator`, `LDScoreResult`, `ChromLDScoreResult`, `run_ldscore()`. | CLI: `ldsc munge-sumstats`; imports: `MungeConfig`, `SumstatsTable`, `SumstatsMunger`, `MungeRunSummary`, `load_sumstats()`. | CLI: `ldsc h2`, `ldsc partitioned-h2`, `ldsc rg`; imports: `RegressionDataset`, `RegressionRunner`, `load_ldscore_from_dir()`. | Imports: `LDScoreOutputConfig`, `LDScoreDirectoryWriter`. |
| Configuration objects | `GlobalConfig` stores shared SNP identifier, genome build, logging, and missing-metadata policy; `get_global_config()`, `set_global_config()`, `reset_global_config()`, and `validate_config_compatibility()` manage and compare shared assumptions. | `AnnotationBuildConfig` normalizes annotation input tokens, output directory, compression, and missing-query policy. | `ReferencePanelBuildConfig` configures PLINK-to-parquet builds; `RefPanelConfig` configures runtime PLINK or parquet-R2 loading, including backend, source paths, metadata, chromosomes, genome build, R2 bias handling, and SNP restriction. | `LDScoreConfig` configures one LD window, MAF filtering, optional individual keep file, optional regression SNP restriction, chunking, and count emission. | `MungeConfig` configures raw sumstats input, trait name, column hints, QC thresholds, sample-size overrides, output directory, merge-alleles path, and legacy munging switches. | `RegressionConfig` configures jackknife blocks, count-vector choice, intercept handling, two-step cutoff, chi-square filtering, and prevalence settings. | `LDScoreOutputConfig` configures the canonical LD-score result directory. |
| Shared normalization and inference | `path_resolution.py` owns path-token language: exact paths, globs, explicit `@` chromosome suites, output-directory creation, and fixed-artifact collision preflight. `column_inference.py` owns header aliases, SNP identifier modes, genome-build aliases, and strict internal artifact schemas. `chromosome_inference.py` owns chromosome normalization and sort order. `genome_build_inference.py` infers `hg19`/`hg38` and 0-based/1-based basis for `chr_pos` inputs when auto mode is used. | Annotation code uses shared path and column inference to resolve `.annot(.gz)` inputs, BED inputs, required metadata columns, restriction files, and chromosome shards; generated `query.<chrom>.annot.gz` files are preflighted before writing. | Builder and runtime loaders use shared path resolution for PLINK prefixes, optional map and chain inputs, parquet R2 groups, metadata sidecars, keep files, and SNP restrictions. The build-ref-panel CLI/convenience wrapper requires an explicit SNP identifier only for SNP restriction files, and the builder preflights candidate parquet/metadata paths before chromosome processing. | The LD-score workflow uses path resolution for annotation groups, PLINK prefixes, parquet R2 groups, metadata sidecars, keep files, and SNP restrictions; it uses column and genome-build inference before calling kernels, and the canonical result directory is preflighted before any LD-score file is written. | Munging resolves one raw file and one optional merge-alleles file; leading raw `##` metadata lines are skipped, raw headers including `#CHROM`/`CHROM` and `POS`/`BP` are inferred through `column_inference.py`, curated `.sumstats` artifacts are loaded with stricter internal column specs, and fixed sumstats/log/metadata outputs are preflighted before the legacy kernel runs. | Regression treats `ldscore_dir` as a literal result directory and sumstats inputs as scalar paths; it validates config compatibility before merging data when both snapshots are known and chooses `SNP` or normalized `CHR:POS` merge keys from the effective identifier mode. Summary TSV outputs are preflighted before writing. | Output code uses `ensure_output_directory()`, `ensure_output_parent_directory()`, and `ensure_output_paths_available()` from `path_resolution.py`; output paths are literal destinations, not input globs, and overwrite must be explicit. |
| Workflow services | `ldsc.cli` dispatches only; scientific work belongs in workflow services. Public workflow modules convert user objects into resolved primitive inputs before entering `_kernel`. | `AnnotationBuilder` loads aligned baseline/query annotation bundles and projects BED intervals onto a baseline SNP grid; `run_bed_to_annot()` is the convenience wrapper. | `ReferencePanelBuilder` orchestrates PLINK chromosome discovery, optional genetic-map loading, optional liftover, SNP/sample restrictions, and fixed parquet artifact writing. `RefPanelLoader` chooses `PlinkRefPanel` or `ParquetR2RefPanel` at LD-score runtime. | `LDScoreCalculator` slices annotation bundles by chromosome; `run_ldscore()` builds explicit annotation bundles when baseline paths exist, synthesizes all-ones `base` when baseline/query inputs are omitted, asks the reference-panel adapter for metadata/readers, calls the LD-score kernel, aggregates `ChromLDScoreResult` objects into `LDScoreResult`, and optionally writes a canonical result directory. | `SumstatsMunger` builds legacy-compatible munging arguments, delegates QC and normalization to the kernel, reloads the curated output, and records a `MungeRunSummary`. | `RegressionRunner` builds `RegressionDataset` objects from `SumstatsTable` plus `LDScoreResult`, drops zero-variance LD columns, and dispatches h2, partitioned-h2, or rg estimators. | `LDScoreDirectoryWriter` writes canonical LD-score directories. |
| Cross-layer data and result objects | `ConfigMismatchError` marks incompatible shared assumptions. `ColumnSpec` describes a canonical field plus accepted aliases. | `AnnotationBundle` carries aligned metadata, baseline annotation matrix, query annotation matrix, column names, chromosomes, source summary, and config snapshot. | `ReferencePanelBuildResult` summarizes emitted parquet artifacts. `RefPanel` defines `available_chromosomes()`, `load_metadata()`, `build_reader()`, `filter_to_snps()`, and `summary()`. | `ChromLDScoreResult` stores one chromosome's split baseline/query tables, count records, SNP universes, and provenance. `LDScoreResult` stores aggregated cross-chromosome baseline/query tables, count records, SNP universes, chromosome results, output paths, and known or unknown provenance. | `SumstatsTable` carries validated `SNP`, `CHR`, `POS`, `Z`, `N`, optional allele columns, provenance, and known or unknown config snapshot. `MungeRunSummary` records row counts, inferred columns, sample-size rule, and output paths. | `RegressionDataset` carries the merged regression table, retained LD columns, regression weight column, selected count vectors, dropped columns, trait names, chromosomes, and provenance. | `LDScoreOutputConfig` describes where the canonical result directory should be written. |
| Internal compute kernels and adapters | `_kernel.formats` contains legacy file readers. `_kernel.identifiers` builds and validates SNP identifiers and reads SNP restriction files. These modules are internal even when public workflows rely on them. | `_kernel.annotation.AnnotationBuilder`, `_kernel.annotation.AnnotationBundle`, and helpers parse annotation tables, validate row alignment, handle chromosome shards, intersect BED intervals, and write `.annot(.gz)` files. | `_kernel.ref_panel_builder` loads provided genetic maps, performs configured liftover, computes pairwise R2 rows, builds standard annotation/LD/metadata schemas, and writes parquet/sidecar files. `_kernel.ref_panel.PlinkRefPanel` and `_kernel.ref_panel.ParquetR2RefPanel` adapt runtime sources to the `RefPanel` interface. | `_kernel.ldscore` owns PLINK/parquet LD readers, annotation parsing helpers, LD-window logic, chromosome LD-score computation, regression weight computation, count calculation, and legacy-compatible LD-score helpers. | `_kernel.sumstats_munger` preserves the historical LDSC munging behavior: column parsing, filters, allele merge, p-to-z conversion, sample-size handling, and `.sumstats.gz` writing. | `_kernel.regression.Hsq`, `_kernel.regression.RG`, `_kernel.regression.LD_Score_Regression`, `_kernel._jackknife.*`, and `_kernel._irwls.IRWLS` implement the numerical estimators and supporting optimization/jackknife routines. | `_kernel` does not own public output layout decisions, except for legacy-compatible writers embedded in historical kernels and the reference-panel builder's low-level parquet/TSV writers. |
| On-disk artifacts owned by the layer | Public docs, tests, and tutorials define user-visible contracts; the CLI exposes names, not internal files. | Inputs: baseline/query `.annot(.gz)` and BED files. Outputs: generated query `.annot(.gz)` files when projection is requested. | Inputs: PLINK `.bed/.bim/.fam`, conditional genetic maps, optional liftover chains, keep files, SNP restrictions. Outputs: `parquet/ann/chr{chrom}_ann.parquet`, `parquet/ld/chr{chrom}_LD.parquet`, and emitted metadata sidecars such as `parquet/meta/chr{chrom}_meta_hg19.tsv.gz` and/or `parquet/meta/chr{chrom}_meta_hg38.tsv.gz`; omitted maps write `CM=NA` for emitted builds. | Inputs: annotation bundles plus PLINK or parquet-R2 reference sources. Outputs: canonical LD-score result directory with `manifest.json`, `baseline.parquet`, and optional `query.parquet`; parquet row groups are chromosome-aligned. | Inputs: raw GWAS summary-statistics file and optional merge-alleles file. Outputs: `sumstats.sumstats.gz`, `sumstats.log`, and `sumstats.metadata.json` under the chosen output directory. | Inputs: canonical `ldscore_dir` plus one or two munged `.sumstats(.gz)` files, with sidecars used when present. Outputs: `h2.tsv`, `partitioned_h2.tsv`, or `rg.tsv` when an output directory is supplied. | `LDScoreDirectoryWriter` owns `manifest.json`, `baseline.parquet`, and `query.parquet`; its manifest records parquet row-group metadata. Regression summary files are written by `regression_runner`; reference-panel parquet/sidecars are written by the builder kernel. |
| Tests and verification | `tests/test_package_layout.py`, `tests/test_config_identifiers.py`, `tests/test_global_config_registry.py`, `tests/test_path_resolution.py`, `tests/test_column_inference.py`, `tests/test_chromosome_inference.py`, and `tests/test_genome_build_inference.py`. | `tests/test_annotation.py`. | `tests/test_ref_panel.py`, `tests/test_ref_panel_builder.py`. | `tests/test_ldscore_workflow.py`, `tests/test_plink_io.py`, `tests/test_output.py`. | `tests/test_sumstats_munger.py`. | `tests/test_regression_workflow.py`, `tests/test_irwls.py`, `tests/test_jackknife.py`. | `tests/test_output.py` plus workflow tests that assert written artifacts can be consumed downstream. |

## Layer Definitions

### Public API and command surface

This layer is the user-facing contract. `ldsc.__init__` re-exports stable
Python objects, while `ldsc.cli` exposes one command with subcommands grouped
by task. Code here should parse, route, and preserve naming contracts; it
should not contain numerical LDSC logic.

### Configuration and shared normalization

Configuration dataclasses describe user choices with frozen, validated objects.
Shared helper modules normalize path tokens, headers, chromosome labels, SNP
identifier modes, and genome-build names before workflow services call the
private kernels. This keeps user-facing flexibility out of the numerical code.

### Workflow services

Workflow modules are the public implementation boundary for each task. They
resolve inputs, align tables, enforce cross-file contracts, select kernel
backends, aggregate per-chromosome work, and hand output payloads to writer
objects. These modules are where most contributor changes should start.

### Cross-layer data and result objects

These dataclasses are the typed handoff objects between services. They carry
tables, column order, selected SNP universes, output paths, and config
provenance. Objects with a `validate()` method define the shape expected by
downstream layers.

### Internal compute kernels and adapters

`ldsc._kernel` contains private readers, legacy-compatible algorithms, and
numerical routines. Kernel code should receive concrete primitive inputs and
resolved tables. It should not be the first place to add public path-token
rules, CLI flags, or user-facing naming changes.

### Output and persistence

The output layer owns public artifact layout. The current LD-score workflow
writes a canonical result directory:

```text
<ldscore_dir>/
  manifest.json
  baseline.parquet
  query.parquet        # omitted when no query annotations exist
```

Regression commands consume that directory through `--ldscore-dir`. Public
regression workflows no longer require callers to pass separate LD-score,
count, weight, or annotation-manifest files.
The LD-score parquet files are single flat files for compatibility, but
`LDScoreDirectoryWriter` writes one row group per chromosome and records the
layout in `manifest.json` as `row_group_layout`, `baseline_row_groups`, and
`query_row_groups`.

All public writers share the same collision policy. A missing output directory
is created with a warning, an existing directory is reused, and fixed files are
checked before the first write. Existing files raise `FileExistsError` unless
the caller passes `--overwrite` on the CLI or `overwrite=True` through the
Python API. The overwrite flag only permits replacement of known workflow files;
it does not delete unrelated files or clean directories.

## Dependency Rules

- Import from `ldsc` for supported public use; avoid importing from
  `ldsc._kernel` in user code.
- CLI code dispatches to workflow modules and does not own scientific logic.
- Workflow modules may depend on shared config/path/schema helpers and private
  kernels.
- Private kernels should not resolve globs, `@` chromosome suites, or other
  public path-token forms.
- LD-score computation remains chromosome-wise; regression consumes aggregated
  LD-score results from memory or from the canonical result directory.
- Output filenames and directory schemas are public contracts and should be
  changed only with coordinated updates to docs, loaders, and tests.
- Fixed output paths should be precomputed and preflighted before long-running
  kernels start, especially for reference-panel and sumstats workflows.
