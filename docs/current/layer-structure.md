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
| Public API and command surface | `ldsc.__init__` re-exports stable objects; `ldsc.__main__.main()` runs `python -m ldsc`; `ldsc.cli.build_parser()` and `ldsc.cli.main()` define the `ldsc` subcommands. Top-level imports include `ChrPosBuildInference`, `infer_chr_pos_build()`, and `resolve_chr_pos_table()`; the CLI has no standalone build-inference subcommand. | CLI: `ldsc annotate`; imports: `AnnotationBuilder`, `AnnotationBundle`, `AnnotationBuildConfig`, `run_bed_to_annot()`. | CLI: `ldsc build-ref-panel`; imports: `ReferencePanelBuilder`, `ReferencePanelBuildResult`, `run_build_ref_panel()`, `RefPanel`, `RefPanelLoader`, `RefPanelConfig`, `PlinkRefPanel`, `ParquetR2RefPanel`. | CLI: `ldsc ldscore`; imports: `LDScoreCalculator`, `LDScoreResult`, `ChromLDScoreResult`, `run_ldscore()`. | CLI: `ldsc munge-sumstats`; imports: `MungeConfig`, `SumstatsTable`, `SumstatsMunger`, `MungeRunSummary`, `load_sumstats()`. | CLI: `ldsc h2`, `ldsc partitioned-h2`, `ldsc rg`; imports: `RegressionDataset`, `RegressionRunner`, `load_ldscore_from_dir()`. | Imports: `LDScoreOutputConfig`, `LDScoreDirectoryWriter`, `PartitionedH2OutputConfig`, `PartitionedH2DirectoryWriter`. Shared implementation: `ldsc._logging` for workflow audit logs. |
| Configuration objects | `GlobalConfig` stores shared SNP identifier, genome build, logging, and missing-metadata policy; `get_global_config()`, `set_global_config()`, `reset_global_config()`, and `validate_config_compatibility()` manage and compare shared assumptions. | `AnnotationBuildConfig` normalizes annotation input tokens, output directory, compression, and missing-query policy. | `ReferencePanelBuildConfig` configures PLINK-to-parquet builds; `RefPanelConfig` configures runtime PLINK or parquet-R2 loading, including backend, source paths, metadata, chromosomes, R2 bias handling, and SNP restriction. Runtime genome-build assumptions come from `GlobalConfig`. | `LDScoreConfig` configures one LD window, optional regression SNP restriction, LD-score `snp_batch_size`, common-MAF count threshold, and whole-chromosome window override. | `MungeConfig` configures `raw_sumstats_file`, trait name, column hints, QC thresholds, sample-size overrides, output directory, optional sumstats SNP keep-list, and legacy munging switches. | `RegressionConfig` configures jackknife blocks, count-vector choice, intercept handling, two-step cutoff, chi-square filtering, and prevalence settings. | `LDScoreOutputConfig` configures the canonical LD-score result directory. |
| Shared normalization and inference | `path_resolution.py` owns path-token language: exact paths, globs, explicit `@` chromosome suites, output-directory creation, fixed-artifact collision preflight, coherent artifact-family preflight, and stale owned-sibling cleanup. `column_inference.py` owns header aliases, SNP identifier modes, genome-build aliases, and strict internal artifact schemas. `chromosome_inference.py` owns chromosome normalization and sort order. `genome_build_inference.py` infers `hg19`/`hg38` and 0-based/1-based basis for `chr_pos` inputs when auto mode is used. `_logging.py` owns LDSC logger level handling and workflow file-handler lifecycle. | Annotation code uses shared path and column inference to resolve `.annot(.gz)` inputs, BED inputs, required metadata columns, restriction files, and chromosome shards; generated root-level `query.*.annot.gz` files and `annotate.log` are preflighted as a family before writing, and stale query shards are removed after successful overwrites. | Builder and runtime loaders use shared path resolution for PLINK prefixes, optional map and chain inputs, parquet R2 groups, metadata sidecars, keep files, SNP restrictions, and duplicate-drop provenance sidecars. The build-ref-panel CLI constructs a `GlobalConfig` only for the identifier mode, the Python wrapper reads the registered identifier mode, source build is explicit or inferred from `.bim`, and the builder preflights candidate R2/metadata/drop paths plus the selected build-ref-panel workflow log before chromosome processing without cleaning stale optional siblings. | The LD-score workflow uses path resolution for annotation groups, PLINK prefixes, parquet R2 groups, metadata sidecars, keep files, and SNP restrictions; it uses column and genome-build inference before calling kernels, and the complete canonical result family plus `ldscore.log` is preflighted before any LD-score file is written. | Munging resolves one raw file and one optional sumstats SNP keep-list; leading raw `##` metadata lines are skipped, raw headers including `#CHROM`/`CHROM` and `POS`/`BP` are inferred through `column_inference.py`, SNP keep-lists are parsed through `_kernel.identifiers`, curated `sumstats.parquet` and legacy `.sumstats(.gz)` artifacts are loaded with stricter internal column specs, and the full sumstats/log/metadata family is preflighted by the workflow before the kernel runs. | Regression treats `ldscore_dir` as a literal result directory and sumstats inputs as scalar paths; it validates config compatibility before merging data when both snapshots are known and chooses `SNP` or normalized `CHR:POS` merge keys from the effective identifier mode. Summary TSV outputs and command logs are preflighted before writing when `output_dir` is supplied; partitioned-h2 preflights `partitioned_h2.tsv`, optional `query_annotations/`, and `partitioned-h2.log` as one family. | Output code uses `ensure_output_directory()`, `ensure_output_parent_directory()`, `ensure_output_paths_available()`, `preflight_output_artifact_family()`, and `remove_output_artifacts()` from `path_resolution.py`; output paths are literal destinations, not input globs, and overwrite must be explicit. Workflow logs are audit artifacts and are excluded from result `output_paths`. |
| Workflow services | `ldsc.cli` dispatches only; scientific work belongs in workflow services. Public workflow modules convert user objects into resolved primitive inputs before entering `_kernel`, then attach shared workflow logging only around orchestration paths that own output directories. | `AnnotationBuilder` loads aligned baseline/query annotation bundles and projects BED intervals onto a baseline SNP grid; `run_bed_to_annot()` is the convenience wrapper and writes `annotate.log`. | `ReferencePanelBuilder` orchestrates PLINK chromosome discovery, optional genetic-map loading, optional liftover, SNP/sample restrictions, fixed parquet artifact writing, and wrapper-owned build-ref-panel workflow logs. `RefPanelLoader` chooses `PlinkRefPanel` or `ParquetR2RefPanel` at LD-score runtime. | `LDScoreCalculator` slices annotation bundles by chromosome; `run_ldscore()` builds explicit annotation bundles when baseline paths exist, synthesizes all-ones `base` when baseline/query inputs are omitted, asks the reference-panel adapter for metadata/readers, calls the LD-score kernel, aggregates `ChromLDScoreResult` objects into `LDScoreResult`, and writes a canonical result directory plus `ldscore.log`. | `run_munge_sumstats_from_args()` maps CLI args to `MungeConfig`; `SumstatsMunger` builds legacy-compatible munging arguments, owns `sumstats.log`, selected Parquet/TSV outputs, and metadata writing, delegates QC and normalization to the kernel, and records a `MungeRunSummary` whose `output_paths` exclude the log. | `RegressionRunner` builds `RegressionDataset` objects from `SumstatsTable` plus `LDScoreResult`, drops zero-variance LD columns, dispatches h2 and rg from baseline LD scores, and dispatches partitioned-h2 only when query LD scores are present; parsed regression entry points write command logs only when `output_dir` is supplied. | `LDScoreDirectoryWriter` writes canonical LD-score directories; `PartitionedH2DirectoryWriter` writes aggregate and optional per-query partitioned-h2 result directories. |
| Cross-layer data and result objects | `ConfigMismatchError` marks incompatible shared assumptions. `ColumnSpec` describes a canonical field plus accepted aliases. | `AnnotationBundle` carries aligned metadata, baseline annotation matrix, query annotation matrix, column names, chromosomes, source summary, and config snapshot. | `ReferencePanelBuildResult` summarizes emitted parquet artifacts. `RefPanel` defines `available_chromosomes()`, `load_metadata()`, `build_reader()`, `filter_to_snps()`, and `summary()`. | `ChromLDScoreResult` stores one chromosome's split baseline/query tables, count records, SNP universes, and provenance. `LDScoreResult` stores aggregated cross-chromosome baseline/query tables, count records, SNP universes, chromosome results, output paths, and known or unknown provenance. | `SumstatsTable` carries validated `SNP`, `CHR`, `POS`, `Z`, `N`, optional allele columns, provenance, and known or unknown config snapshot. `MungeRunSummary` records row counts, inferred columns, sample-size rule, and output paths. | `RegressionDataset` carries the merged regression table, retained LD columns, regression weight column, selected count vectors, dropped columns, trait names, chromosomes, and provenance. | `LDScoreOutputConfig` describes where the canonical result directory should be written. |
| Internal compute kernels and adapters | `_kernel.formats` contains legacy file readers. `_kernel.identifiers` builds and validates SNP identifiers and reads SNP restriction files. These modules are internal even when public workflows rely on them. | `_kernel.annotation` contains low-level annotation table readers, BED normalization, and BED/SNP-grid intersection helpers. `AnnotationBuilder`, `AnnotationBundle`, CLI parsing, chromosome-shard orchestration, and `.annot(.gz)` output policy live in `ldsc.annotation_builder`. | `_kernel.ref_panel_builder` loads provided genetic maps, performs configured liftover, computes pairwise R2 rows, builds in-memory reference SNP tables, and writes canonical R2 parquet plus metadata sidecars. `_kernel.ref_panel.PlinkRefPanel` and `_kernel.ref_panel.ParquetR2RefPanel` adapt runtime sources to the `RefPanel` interface. | `_kernel.ldscore` owns PLINK/parquet R2 readers, the canonical parquet decoded row-group cache, annotation parsing helpers, LD-window logic, chromosome LD-score computation, regression weight computation, count calculation, and legacy-compatible LD-score helpers. | `_kernel.sumstats_munger` preserves the historical LDSC munging behavior: column parsing, filters, allele merge, p-to-z conversion, sample-size handling, and `.sumstats.gz` writing through package logger messages rather than direct log-file ownership. | `_kernel.regression.Hsq`, `_kernel.regression.RG`, `_kernel.regression.LD_Score_Regression`, `_kernel._jackknife.*`, and `_kernel._irwls.IRWLS` implement the numerical estimators and supporting optimization/jackknife routines. | `_kernel` does not own public output layout decisions, except for legacy-compatible writers embedded in historical kernels and the reference-panel builder's low-level parquet/TSV writers. |
| On-disk artifacts owned by the layer | Public docs, tests, and tutorials define user-visible contracts; the CLI exposes names, not internal files. | Inputs: baseline/query `.annot(.gz)` and BED files. Outputs: generated query `.annot(.gz)` files and `annotate.log` when projection is requested. | Inputs: PLINK `.bed/.bim/.fam`, conditional genetic maps, optional liftover chains, keep files, SNP restrictions. Outputs: build-specific `{build}/chr{chrom}_r2.parquet`, `{build}/chr{chrom}_meta.tsv.gz`, optional `dropped_snps/chr{chrom}_dropped.tsv.gz`, and `build-ref-panel.log` or `build-ref-panel.chr<chrom>.log`; omitted maps write `CM=NA` for emitted SNP/kb-window builds. | Inputs: annotation bundles plus a PLINK prefix or build-specific `r2_dir`. Outputs: canonical LD-score result directory with `manifest.json`, `ldscore.baseline.parquet`, optional `ldscore.query.parquet`, and `ldscore.log`; parquet row groups are chromosome-aligned. | Inputs: raw GWAS summary-statistics file and optional sumstats SNP keep-list. Outputs: `sumstats.parquet` by default, optional `sumstats.sumstats.gz`, `sumstats.log`, and `sumstats.metadata.json` under the chosen output directory. | Inputs: canonical `ldscore_dir` plus one or two munged `sumstats.parquet` or `.sumstats(.gz)` files, with sidecars used when present; partitioned-h2 additionally requires `ldscore.query.parquet` and non-empty `query_columns`. Outputs: `h2.tsv` plus `h2.log`, compact `partitioned_h2.tsv` plus `partitioned-h2.log`, optional `query_annotations/manifest.tsv`, per-query `partitioned_h2.tsv`, `partitioned_h2_full.tsv`, and `metadata.json`, or `rg.tsv` plus `rg.log` when an output directory is supplied. | `LDScoreDirectoryWriter` owns `manifest.json`, `ldscore.baseline.parquet`, and `ldscore.query.parquet`; its manifest records parquet row-group metadata. `PartitionedH2DirectoryWriter` owns compact `partitioned_h2.tsv` and optional per-query result trees; reference-panel parquet/sidecars are written by the builder kernel. Workflow logs are written by `ldsc._logging` contexts in the workflow layer. |
| Tests and verification | `tests/test_package_layout.py`, `tests/test_config_identifiers.py`, `tests/test_global_config_registry.py`, `tests/test_path_resolution.py`, `tests/test_column_inference.py`, `tests/test_chromosome_inference.py`, `tests/test_genome_build_inference.py`, and `tests/test_logging_refactor.py`. | `tests/test_annotation.py`. | `tests/test_ref_panel.py`, `tests/test_ref_panel_builder.py`. | `tests/test_ldscore_workflow.py`, `tests/test_plink_io.py`, `tests/test_output.py`. | `tests/test_sumstats_munger.py`. | `tests/test_regression_workflow.py`, `tests/test_irwls.py`, `tests/test_jackknife.py`. | `tests/test_output.py` plus workflow tests that assert written artifacts can be consumed downstream. |

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
  ldscore.baseline.parquet
  ldscore.query.parquet        # omitted when no query annotations exist
```

Regression commands consume that directory through `--ldscore-dir`. Public
regression workflows no longer require callers to pass separate LD-score,
count, weight, or annotation-manifest files.
The LD-score parquet files are single flat files for compatibility, but
`LDScoreDirectoryWriter` writes one row group per chromosome and records the
layout in `manifest.json` as `row_group_layout`, `baseline_row_groups`, and
`query_row_groups`.

Public result-directory writers use coherent artifact families. A missing
output directory is created with a warning, an existing directory is reused, and
owned workflow artifacts are checked before the first write. Existing owned
siblings raise `FileExistsError` unless the caller passes `--overwrite` on the
CLI or `overwrite=True` through the Python API, even when the current run would
not produce that sibling. With overwrite enabled, successful runs remove stale
owned siblings not produced by the current configuration. Unrelated files are
preserved.

This family policy applies to `munge-sumstats`, `ldscore`, `partitioned-h2`,
and `annotate`. `build-ref-panel` is documented as an exception: its overwrite
mode replaces current candidate artifacts but does not clean stale optional
target-build or `dropped_snps` siblings from earlier expert configurations.

Workflow logs follow that same preflight policy. The log path is checked with
the scientific outputs before the log handler opens the file; if computation
fails after that point, the file is kept with a `Failed` lifecycle footer.
Logs are not included in `output_paths` result mappings.

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
- Fixed output paths, including workflow logs, should be precomputed and
  preflighted before long-running kernels start, especially for reference-panel
  and sumstats workflows.
