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
| Public API and command surface | `ldsc.__init__` re-exports stable objects; `ldsc.__main__.main()` runs `python -m ldsc`; `ldsc.cli.build_parser()` and `ldsc.cli.main()` define the `ldsc` subcommands. Top-level imports include `ChrPosBuildInference`, `infer_chr_pos_build()`, and `resolve_chr_pos_table()`; the CLI has no standalone build-inference subcommand. | CLI: `ldsc annotate`; imports: `AnnotationBuilder`, `AnnotationBundle`, `AnnotationBuildConfig`, `run_bed_to_annot()`. | CLI: `ldsc build-ref-panel`; imports: `ReferencePanelBuilder`, `ReferencePanelBuildResult`, `run_build_ref_panel()`, `RefPanel`, `RefPanelLoader`, `RefPanelConfig`, `PlinkRefPanel`, `ParquetR2RefPanel`. | CLI: `ldsc ldscore`; imports: `LDScoreCalculator`, `LDScoreResult`, `ChromLDScoreResult`, `run_ldscore()`. | CLI: `ldsc munge-sumstats`; imports: `MungeConfig`, `SumstatsTable`, `SumstatsMunger`, `MungeRunSummary`, `load_sumstats()`. | CLI: `ldsc h2`, `ldsc partitioned-h2`, `ldsc rg`; imports: `RegressionDataset`, `RegressionRunner`, `RgResultFamily`, `load_ldscore_from_dir()`. | Imports: `LDScoreOutputConfig`, `LDScoreDirectoryWriter`, `PartitionedH2OutputConfig`, `PartitionedH2DirectoryWriter`, `RgOutputConfig`, `RgDirectoryWriter`. Shared implementation: `ldsc._logging` for workflow audit logs. |
| Configuration objects | `GlobalConfig` stores shared SNP identifier, genome build, logging, and missing-metadata policy; `get_global_config()`, `set_global_config()`, `reset_global_config()`, and `validate_config_compatibility()` manage and compare shared assumptions. | `AnnotationBuildConfig` normalizes annotation input tokens, output directory, compression, and missing-query policy. | `ReferencePanelBuildConfig` configures PLINK-to-parquet builds; `RefPanelConfig` configures runtime PLINK or parquet-R2 loading, including backend, source paths, metadata, chromosomes, R2 bias handling, explicit SNP restriction, and HM3 SNP restriction. Runtime genome-build assumptions come from `GlobalConfig`. | `LDScoreConfig` configures one LD window, optional explicit or HM3 regression SNP restriction, LD-score `snp_batch_size`, common-MAF count threshold, and whole-chromosome window override. | `MungeConfig` configures `raw_sumstats_file`, trait name, raw `sumstats_format`, INFO-list columns, column hints, QC thresholds, sample-size overrides, output directory, optional explicit or HM3 sumstats SNP restriction, optional `chr_pos` liftover target/method, and legacy munging switches. | `RegressionConfig` configures jackknife blocks, count-vector choice, intercept handling, two-step cutoff, chi-square filtering, and prevalence settings. | `LDScoreOutputConfig`, `PartitionedH2OutputConfig`, and `RgOutputConfig` configure fixed result directories. |
| Shared normalization and inference | `path_resolution.py` owns path-token language: exact paths, globs, explicit `@` chromosome suites, output-directory creation, fixed-artifact collision preflight, coherent artifact-family preflight, and stale owned-sibling cleanup. `column_inference.py` owns header aliases, genome-build aliases, and strict internal artifact schemas. `_kernel.snp_identity` owns the exact public SNP identifier modes (`rsid`, `rsid_allele_aware`, `chr_pos`, `chr_pos_allele_aware`), default semantics, effective keys, allele-aware cleanup, restriction matching, and artifact provenance. `chromosome_inference.py` owns chromosome normalization and sort order. `genome_build_inference.py` infers `hg19`/`hg38` and 0-based/1-based basis for `chr_pos`-family inputs when auto mode is used. `hm3.py` owns the public packaged curated HM3 map loader, and `hm3_reference.py` rebuilds the compact HM3 coordinate inference reference; `_kernel.liftover` owns shared chain-file and HM3 coordinate conversion. `_logging.py` owns LDSC logger level handling and workflow file-handler lifecycle. | Annotation code uses shared path and column inference to resolve `.annot(.gz)` inputs, BED inputs, required metadata columns, restriction files, and chromosome shards; generated root-level `query.*.annot.gz` files and diagnostics under `diagnostics/` are preflighted as a family before writing, and stale query shards are removed after successful overwrites. Annotation files may omit alleles in allele-aware modes unless they intend allele-aware matching. | Builder and runtime loaders use shared path resolution for PLINK prefixes, optional map and chain inputs, parquet R2 groups, metadata sidecars, keep files, explicit or HM3 SNP restrictions, and liftover-stage drop audit sidecars. The build-ref-panel CLI constructs a `GlobalConfig` only for the identifier mode, the Python wrapper reads the registered identifier mode, source build is explicit or inferred from `.bim`, and the builder preflights candidate R2/metadata/drop paths plus selected diagnostics before chromosome processing; successful overwrites remove stale current-contract build, chromosome, dropped-SNP, and log siblings. | The LD-score workflow uses path resolution for annotation groups, PLINK prefixes, parquet R2 groups, metadata sidecars, keep files, explicit or HM3 SNP restrictions; it uses column, identity, and genome-build inference before calling kernels, and the complete canonical result family plus `diagnostics/ldscore.log` is preflighted before any LD-score file is written. External raw R2 parquet is base-mode only; allele-aware modes require package-built endpoint alleles. | Munging resolves one raw file and one optional explicit or HM3 sumstats SNP restriction; leading raw `##` metadata lines are skipped, `--format auto` detects common plain text, including VCF-style headers, old DANER, and new DANER, safe aliases are inferred through `column_inference.py`, `--infer-only` reports missing fields and suggested repair flags without writing outputs, SNP restrictions are parsed through `_kernel.identifiers`, optional `chr_pos`-family liftover uses `_kernel.liftover`, curated `sumstats.parquet`, legacy `.sumstats(.gz)`, and dropped-SNP audit artifacts are loaded or written with stricter internal column specs, and the full sumstats root-metadata/data/diagnostics family is preflighted by the workflow before the kernel runs. | Regression treats `ldscore_dir` as a literal result directory and h2/partitioned-h2 sumstats inputs as scalar paths; rg resolves `--sumstats-sources` as a file group and can optionally select one anchor source. It validates config compatibility before merging data when both snapshots are known and chooses the effective SNP identity key for the resolved mode; `--allow-identity-downgrade` is regression-only and same-family only. Summary TSV outputs and command logs are preflighted before writing when `output_dir` is supplied; without `output_dir`, regression CLIs print compact TSV to stdout and write no diagnostics. Partitioned-h2 preflights `partitioned_h2.tsv`, optional `diagnostics/query_annotations/`, and `diagnostics/partitioned-h2.log` as one family; rg preflights `rg.tsv`, `rg_full.tsv`, `h2_per_trait.tsv`, optional `diagnostics/pairs/`, and `diagnostics/rg.log`. | Output code uses `ensure_output_directory()`, `ensure_output_parent_directory()`, `ensure_output_paths_available()`, `preflight_output_artifact_family()`, and `remove_output_artifacts()` from `path_resolution.py`; output paths are literal destinations, not input globs, and overwrite must be explicit. Workflow logs are audit artifacts and are excluded from result `output_paths`. |
| Workflow services | `ldsc.cli` dispatches only; scientific work belongs in workflow services. Public workflow modules convert user objects into resolved primitive inputs before entering `_kernel`, then attach shared workflow logging only around orchestration paths that own output directories. | `AnnotationBuilder` loads aligned baseline/query annotation bundles and projects BED intervals onto a baseline SNP grid; `run_bed_to_annot()` is the convenience wrapper and writes `diagnostics/annotate.log`. | `ReferencePanelBuilder` orchestrates PLINK chromosome discovery, optional genetic-map loading, optional liftover, SNP/sample restrictions, fixed parquet artifact writing, dropped-SNP sidecar writing, diagnostic metadata, and wrapper-owned build-ref-panel workflow logs under `diagnostics/`. `RefPanelLoader` chooses `PlinkRefPanel` or `ParquetR2RefPanel` at LD-score runtime. | `LDScoreCalculator` slices annotation bundles by chromosome; `run_ldscore()` builds explicit annotation bundles when baseline paths exist, synthesizes all-ones `base` when baseline/query inputs are omitted, asks the reference-panel adapter for metadata/readers, calls the LD-score kernel, aggregates `ChromLDScoreResult` objects into `LDScoreResult`, and writes a canonical result directory plus `diagnostics/ldscore.log`. | `run_munge_sumstats_from_args()` maps CLI args to `MungeConfig`; `SumstatsMunger` builds legacy-compatible munging arguments, owns root `metadata.json`, selected Parquet/TSV outputs, diagnostics, and dropped-SNP sidecar writing, delegates QC and normalization to the kernel, and records a `MungeRunSummary` whose `output_paths` exclude the log. | `RegressionRunner` builds `RegressionDataset` and `RGRegressionDataset` objects from `SumstatsTable` plus `LDScoreResult`, drops zero-variance LD columns, dispatches h2 and rg from baseline LD scores, dispatches partitioned-h2 only when query LD scores are present, and exposes `estimate_rg_pairs()` for all-pairs or anchor-vs-rest rg result families; parsed regression entry points write command logs only when `output_dir` is supplied and otherwise print compact TSV to stdout. | `LDScoreDirectoryWriter` writes canonical LD-score directories; `PartitionedH2DirectoryWriter` writes aggregate and optional diagnostic per-query partitioned-h2 result directories; `RgDirectoryWriter` writes the rg result family and optional diagnostic per-pair detail tree. |
| Cross-layer data and result objects | `ConfigMismatchError` marks incompatible shared assumptions. `ColumnSpec` describes a canonical field plus accepted aliases. | `AnnotationBundle` carries aligned metadata, baseline annotation matrix, query annotation matrix, column names, chromosomes, source summary, and config snapshot. | `ReferencePanelBuildResult` summarizes emitted parquet artifacts. `RefPanel` defines `available_chromosomes()`, `load_metadata()`, `build_reader()`, `filter_to_snps()`, and `summary()`. | `ChromLDScoreResult` stores one chromosome's split baseline/query tables, count records, SNP universes, and provenance. `LDScoreResult` stores aggregated cross-chromosome baseline/query tables, count records, SNP universes, chromosome results, output paths, and current identity provenance. | `SumstatsTable` carries validated `SNP`, `CHR`, `POS`, `Z`, `N`, optional allele columns, provenance, and the current config snapshot. `MungeRunSummary` records row counts, inferred columns, sample-size rule, and output paths. | `RegressionDataset` carries the merged regression table, retained LD columns, regression LD-score column, selected count vectors, dropped columns, trait names, chromosomes, and provenance. `RGRegressionDataset` carries the two-trait merge used by one rg kernel call. `RgResultFamily` carries concise rg, full rg diagnostics, per-trait h2, and per-pair metadata. | `LDScoreOutputConfig`, `PartitionedH2OutputConfig`, and `RgOutputConfig` describe fixed result-directory layouts. |
| Internal compute kernels and adapters | `_kernel.formats` contains legacy file readers. `_kernel.identifiers` parses SNP restriction files and delegates exact mode semantics to `_kernel.snp_identity`, which builds base/effective keys, performs allele-aware cleanup, collapses restriction keys, and validates minimal identity provenance. `_kernel.liftover` applies chain-file and HM3 quick liftover, detects duplicate coordinate groups, and returns readable drop-count provenance; the public packaged HM3 loader lives in `ldsc.hm3`. These modules are internal even when public workflows rely on them. | `_kernel.annotation` contains low-level annotation table readers, BED normalization, and BED/SNP-grid intersection helpers. `AnnotationBuilder`, `AnnotationBundle`, CLI parsing, chromosome-shard orchestration, and `.annot(.gz)` output policy live in `ldsc.annotation_builder`. | `_kernel.ref_panel_builder` loads provided genetic maps, performs configured coordinate liftover through the shared liftover helper, computes pairwise R2 rows, builds in-memory reference SNP tables, and writes canonical R2 parquet plus metadata sidecars. `_kernel.ref_panel.PlinkRefPanel` and `_kernel.ref_panel.ParquetR2RefPanel` adapt runtime sources to the `RefPanel` interface. | `_kernel.ldscore` owns PLINK/parquet R2 readers, the canonical parquet decoded row-group cache, annotation parsing helpers, LD-window logic, chromosome LD-score computation, regression-universe LD-score computation, count calculation, and legacy-compatible LD-score helpers. | `_kernel.sumstats_munger` preserves the historical LDSC munging behavior: column parsing, filters, allele merge, p-to-z conversion, sample-size handling, optional post-filter coordinate liftover, identity cleanup, and `.sumstats.gz` writing through package logger messages rather than direct log-file ownership. | `_kernel.regression.Hsq`, `_kernel.regression.RG`, `_kernel.regression.LD_Score_Regression`, `_kernel._jackknife.*`, and `_kernel._irwls.IRWLS` implement the numerical estimators and supporting optimization/jackknife routines. | `_kernel` does not own public output layout decisions, except for legacy-compatible writers embedded in historical kernels and the reference-panel builder's low-level parquet/TSV writers. |
| On-disk artifacts owned by the layer | Public docs, tests, and tutorials define user-visible contracts; the CLI exposes names, not internal files. | Inputs: baseline/query `.annot(.gz)` and BED files. Outputs: generated root query `.annot(.gz)` files plus diagnostics under `diagnostics/` when projection is requested. | Inputs: PLINK `.bed/.bim/.fam`, conditional genetic maps, optional liftover chains, keep files, SNP restrictions. Outputs: build-specific `{build}/chr{chrom}_r2.parquet`, `{build}/chr{chrom}_meta.tsv.gz`, diagnostic metadata, always-written `diagnostics/dropped_snps/chr{chrom}_dropped.tsv.gz` for processed chromosomes, and build logs under `diagnostics/`; omitted maps write `CM=NA` for emitted SNP/kb-window builds. | Inputs: annotation bundles plus a PLINK prefix or build-specific `r2_dir`. Outputs: canonical LD-score result directory with root `metadata.json`, `ldscore.baseline.parquet`, optional `ldscore.query.parquet`, and `diagnostics/ldscore.log`; parquet row groups are chromosome-aligned. | Inputs: raw GWAS summary-statistics file and optional sumstats SNP keep-list. Outputs: root `metadata.json`, `sumstats.parquet` by default, optional `sumstats.sumstats.gz`, `diagnostics/sumstats.log`, and always-written `diagnostics/dropped_snps/dropped.tsv.gz` under the chosen output directory. | Inputs: canonical `ldscore_dir` plus one or more munged `sumstats.parquet` or `.sumstats(.gz)` files, with current root metadata required for package-written artifacts; h2 uses one trait, rg uses two or more traits through `--sumstats-sources`, and partitioned-h2 additionally requires `ldscore.query.parquet` and non-empty `query_columns`. Outputs with `output_dir`: `h2.tsv` plus diagnostics, compact `partitioned_h2.tsv` plus diagnostics, optional `diagnostics/query_annotations/manifest.tsv`, per-query `partitioned_h2.tsv`, `partitioned_h2_full.tsv`, and `metadata.json`, or the rg family (`rg.tsv`, `rg_full.tsv`, `h2_per_trait.tsv`, optional `diagnostics/pairs/`) plus diagnostics. Without `output_dir`, regression CLIs print compact TSV to stdout. | `LDScoreDirectoryWriter` owns `metadata.json`, `ldscore.baseline.parquet`, and `ldscore.query.parquet`; root metadata records parquet row-group metadata. `H2DirectoryWriter` owns `h2.tsv` plus diagnostic metadata. `PartitionedH2DirectoryWriter` owns compact `partitioned_h2.tsv` and optional diagnostic per-query result trees; `RgDirectoryWriter` owns rg data files and optional diagnostic pair trees but not workflow logs; reference-panel parquet/sidecars are written by the builder kernel. Workflow logs are written by `ldsc._logging` contexts in the workflow layer. |
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
  metadata.json
  ldscore.baseline.parquet
  ldscore.query.parquet        # omitted when no query annotations exist
  diagnostics/
    ldscore.log
```

Regression commands consume that directory through `--ldscore-dir`. Public
regression workflows no longer require callers to pass separate LD-score,
count, weight, or annotation-manifest files.
The LD-score parquet files are single flat files for compatibility, but
`LDScoreDirectoryWriter` writes one row group per chromosome and records the
layout in `metadata.json` as `row_group_layout`, `baseline_row_groups`, and
`query_row_groups`.

Public result-directory writers use coherent artifact families. A missing
output directory is created with a warning, an existing directory is reused, and
current-contract owned workflow artifacts are checked before the first write.
Existing owned siblings raise `FileExistsError` unless the caller passes
`--overwrite` on the CLI or `overwrite=True` through the Python API, even when
the current run would not produce that sibling. With overwrite enabled,
successful runs remove stale current-contract owned siblings not produced by the
current configuration. Removed legacy root diagnostic names are ignored, and
unrelated files are preserved.

This family policy applies to `munge-sumstats`, `build-ref-panel`, `ldscore`,
`partitioned-h2`, `rg`, and `annotate`.

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
