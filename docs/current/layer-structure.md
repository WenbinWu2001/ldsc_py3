# Layer Structure

Last updated on: 2026-09-07

This document maps the refactored `ldsc` package by layer and by main
functionality. Use it to answer two questions quickly:

The continuous-annotation extension keeps semantics and post-fit work in the public layers: `annotation_semantics.py` owns advisory classification and compact hashes, `regression_runner.py` persists one delete matrix per fitted model, `quantile_h2.py` owns source reconstruction and projection, and `outputs.py` owns the result family. The regression kernel is unchanged.

Optional result exploration follows the same direction. `plotting` validates
canonical result metadata and owns derived figure output, while its private
`_builders` module alone imports Matplotlib. `h2_scale` loads a canonical h2
result and reuses the regression kernel's conversion factor. Neither path is
called automatically by a core workflow.

- What layer owns this object?
- Which objects participate in a given workflow?

The public import boundary is `ldsc` and the public command boundary is
`ldsc.cli`. Modules under `ldsc._kernel` are internal implementation details.

Gene-list support follows the same boundary. `ldsc.gene_list_resolver` validates
the required or embedded catalog and returns compact batch resolution/audit records;
`ldsc.query_annotations` defines internal BED/gene status records;
`AnnotationBuilder` projects intervals; `LDScoreCalculator` prunes zero-hit and
zero-variance queries; and `LDScoreDirectoryWriter` owns the status and
unresolved-gene diagnostics. No gene identity enters `_kernel`.

Legacy compatibility follows the same boundary. `ldsc.legacy_ldscore_converter`
alone reads selected LDSC2 LD-score/count/annotation/frequency suites and emits
a canonical `LDScoreResult` directory; it is exported as
`LegacyLDScoreConverter` and `convert_ldsc2_ldscores()`. Regression never reads
those fragments. `ldsc.regression_runner` automatically projects legacy LDSC2
sumstats text onto the already-loaded canonical panel and owns the resulting
drop audit. Neither compatibility path adds a legacy reader or writer to
`_kernel`.

## Layer Matrix

| Layer | Shared CLI, config, path, and schema contracts | Annotation building | Reference-panel building and loading | LD-score calculation | Summary-statistics munging | Regression | Output and persistence |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Public API and command surface | `ldsc.__init__` re-exports stable objects; `ldsc.__main__.main()` runs `python -m ldsc`; `ldsc.cli.build_parser()` and `ldsc.cli.main()` define the `ldsc` subcommands. Top-level imports include build inference, R2 query, lazy `plot_result()`/`PlotArtifact`, and lazy `convert_h2_scale()`/`H2ScaleConversionArtifact`; the CLI has no standalone build-inference subcommand. | CLI: `ldsc annotate`; imports: `AnnotationBuilder`, `AnnotationBundle`, `AnnotationBuildConfig`, `run_bed_to_annot()`. | CLI: `ldsc build-ref-panel`, `ldsc query-r2`; imports: `ReferencePanelBuilder`, `ReferencePanelBuildResult`, `run_build_ref_panel()`, `R2Panel`, `query_r2()`, `RefPanel`, `RefPanelLoader`, `RefPanelConfig`, `PlinkRefPanel`, `ParquetR2RefPanel`. | CLI: `ldsc ldscore`; imports: `LDScoreCalculator`, `LDScoreResult`, `ChromLDScoreResult`, `run_ldscore()`. | CLI: `ldsc munge-sumstats`; imports: `MungeConfig`, `SumstatsTable`, `SumstatsMunger`, `MungeRunSummary`, `load_sumstats()`. | CLI: `ldsc h2`, `ldsc partitioned-h2`, `ldsc rg`, `ldsc convert-h2-scale`, `ldsc plot`; imports include regression types plus the lazy conversion and plotting APIs. | Imports: canonical writer configs/classes. `ldsc.plotting` owns fixed nested plot families; `ldsc.h2_scale` owns scale-conversion families; `ldsc._logging` owns audit logs and failure markers. |
| Configuration objects | `GlobalConfig` stores shared SNP identifier, genome build, logging, and missing-metadata policy; `get_global_config()`, `set_global_config()`, `reset_global_config()`, and `validate_config_compatibility()` manage and compare shared assumptions. | `AnnotationBuildConfig` normalizes annotation input tokens, output directory, compression, and missing-query policy. | `ReferencePanelBuildConfig` configures PLINK-to-parquet builds, source-build inference, optional coordinate liftover, explicit SNP/sample restrictions, and optional `min_r2`; `RefPanelConfig` configures runtime PLINK or parquet-R2 loading, including backend, source paths, chromosomes, explicit SNP/sample restrictions, retained-panel MAF filtering, and optional genetic maps. Runtime genome-build assumptions come from `GlobalConfig`. | `LDScoreConfig` configures one LD window, optional explicit or bundled-HM3-default regression SNP restriction, LD-score `snp_batch_size`, common-MAF count threshold, and whole-chromosome window override. | `MungeConfig` configures `raw_sumstats_file`, trait name, raw `sumstats_format`, INFO-list columns, column hints, QC thresholds, sample-size overrides, output directory, optional explicit or HM3 sumstats SNP restriction, optional `chr_pos` liftover target/method, and legacy munging switches. | `RegressionConfig` configures jackknife blocks, count-vector choice, intercept handling, two-step cutoff, chi-square filtering, and prevalence settings. | `LDScoreOutputConfig`, `PartitionedH2OutputConfig`, and `RgOutputConfig` configure fixed result directories. |
| Shared normalization and inference | `path_resolution.py`, schema/build inference, SNP identity, liftover, and `_logging.py` own shared normalization, output preflight, lifecycle logging, and failure-marker I/O. Core CLI workflows require `--output-dir`; `plot` and `convert-h2-scale` derive fixed nested destinations. | Annotation inputs and generated query shards use shared resolution and coherent artifact-family preflight. | Builders/loaders resolve PLINK, R2, map, restriction, and liftover inputs; build outputs are preflighted by chromosome or suite. | LD-score resolution validates annotations, panels, identity/build metadata, and its complete output family before writing. | Munging infers raw formats and columns, applies restrictions/liftover, and preflights its curated artifacts; infer-only reports without writing. | Regression validates compatibility before merging and always preflights the required output directory. Query partitioned-h2 owns the per-query tree by default; baseline-only complete-model artifacts stay at the root. RG pair details remain opt-in. Plot dispatch checks only plotting-required metadata fields. | Output paths are literal or fixed derived destinations, overwrite is explicit, stale owned siblings are cleaned, and workflow logs/markers are excluded from scientific result mappings. |
| Workflow services | `ldsc.cli` dispatches only; public workflow services resolve user inputs, call kernels, and own output/log orchestration. | `AnnotationBuilder.run()` honors `AnnotationBuildConfig.output_dir`; materializing annotation workflows write query shards, diagnostics, and `diagnostics/annotate.log`. `run_bed_to_annot()` uses the same stream. | `ReferencePanelBuilder` writes panel artifacts and diagnostics; `R2Panel` serves required-output pair queries; `RefPanelLoader` chooses runtime backends. | `LDScoreCalculator` computes in memory when no output is requested, but metadata export requires output configuration; `run_ldscore()` requires and writes a canonical output directory. | `run_munge_sumstats_from_args()` maps CLI args to `MungeConfig`; `SumstatsMunger` owns curated outputs and diagnostics. | `RegressionRunner` retains in-memory numerical APIs and the current batch return type. Parsed CLI workflows require output and write command logs. Query batch runs materialize per-query results by default. | Directory writers own canonical result layouts; low-level Python output configs retain their explicit controls. |
| Cross-layer data and result objects | `ConfigMismatchError` marks incompatible shared assumptions. `ColumnSpec` describes a canonical field plus accepted aliases. | `AnnotationBundle` carries aligned metadata, baseline annotation matrix, query annotation matrix, column names, chromosomes, source summary, and config snapshot. | `ReferencePanelBuildResult` summarizes emitted parquet artifacts. `R2Panel` caches per-chromosome sidecar/index state for repeated pair queries. `RefPanel` defines `available_chromosomes()`, `load_metadata()`, `build_reader()`, `filter_to_snps()`, and `summary()`. | `ChromLDScoreResult` stores one chromosome's split baseline/query tables, count records, SNP universes, and provenance. `LDScoreResult` stores aggregated cross-chromosome baseline/query tables, count records, SNP universes, chromosome results, output paths, and current identity provenance. | `SumstatsTable` carries validated `SNP`, `CHR`, `POS`, `Z`, `N`, optional allele columns, provenance, and the current config snapshot. `MungeRunSummary` records row counts, inferred columns, sample-size rule, and output paths. | `RegressionDataset` carries the merged regression table, retained LD columns, regression LD-score column, selected count vectors, dropped columns, trait names, chromosomes, and provenance. `RGRegressionDataset` carries the two-trait merge used by one rg kernel call. `RgResultFamily` carries concise rg, full rg diagnostics, per-trait h2, and per-pair metadata; rg p-values are nominal and uncorrected. | `LDScoreOutputConfig`, `PartitionedH2OutputConfig`, and `RgOutputConfig` describe fixed result-directory layouts. |
| Internal compute kernels and adapters | `_kernel.formats` contains only retained PLINK and identity-list reader primitives; obsolete legacy sumstats, LD-score, count, annotation, and frequency readers were removed. `_kernel.identifiers` parses SNP restriction files and delegates exact mode semantics to `_kernel.snp_identity`, which builds base/effective keys, performs allele-aware cleanup, collapses restriction keys, and validates minimal identity provenance. `_kernel.liftover` applies chain-file and HM3 quick liftover, detects duplicate coordinate groups, and returns readable drop-count provenance; `_kernel.regions` loads packaged/custom exclusion intervals; the public packaged HM3 loader lives in `ldsc.hm3`. These modules are internal even when public workflows rely on them. | `_kernel.annotation` contains low-level annotation table readers, BED normalization, and BED/SNP-grid intersection helpers. `AnnotationBuilder`, `AnnotationBundle`, CLI parsing, chromosome-shard orchestration, and `.annot(.gz)` output policy live in `ldsc.annotation_builder`. | `_kernel.ref_panel_builder` loads provided genetic maps, performs configured coordinate liftover through the shared liftover helper, computes pairwise R2 rows, and writes canonical 4-column index R2 parquet files plus mandatory metadata sidecars. `_kernel.r2_query` performs low-level index-format parquet pair lookup. `_kernel.ref_panel.PlinkRefPanel` and `_kernel.ref_panel.ParquetR2RefPanel` adapt runtime sources to the `RefPanel` interface. | `_kernel.ldscore` owns PLINK/parquet R2 readers, the index-format parquet pair-streaming read path (`iter_all_pairs` + `ld_score_streaming_from_r2_reader`, with binding validation and per-chromosome remap), annotation parsing helpers, LD-window logic, chromosome LD-score computation, regression-universe LD-score computation, count calculation; it has no LDSC2 artifact emitters. | `_kernel.sumstats_munger` preserves the historical LDSC munging behavior: column parsing, filters, allele merge, p-to-z conversion, sample-size handling, optional post-filter coordinate liftover, and identity cleanup; it returns in-memory tables and owns no output files. | `_kernel.regression.Hsq`, `_kernel.regression.RG`, `_kernel.regression.LD_Score_Regression`, `_kernel._jackknife.*`, and `_kernel._irwls.IRWLS` implement the numerical estimators and supporting optimization/jackknife routines. | `_kernel` does not own public output layout decisions; reference-panel low-level parquet/TSV emission remains part of the dedicated builder kernel. |
| On-disk artifacts owned by the layer | Public docs, tests, and tutorials define user-visible contracts; core CLI commands take a required output directory while derived commands take a source result root. | Materializing annotation runs write root `query.<chrom>.annot.gz` shards plus diagnostics and a log. | Panel builds write R2 parquet, metadata, drop audits, and logs. `query-r2` always writes `query_r2.tsv` plus diagnostic metadata and log under its output directory. | LD-score runs write canonical metadata, baseline/query/overlap parquet as applicable, and `diagnostics/ldscore.log`. | Munging writes self-describing parquet, optional legacy TSV, drop audit, and log; infer-only writes nothing. | H2 also writes exact regression-bin diagnostics; partitioned-h2, quantile-h2, and rg write stable result families. `plot` writes one PNG below `plots/`; h2 conversion writes a table and optional sensitivity PNG below `postprocessing/`. | Writers own core layouts; plotting/conversion own derived layouts; workflow logging owns logs and non-scientific `RUN_FAILED` markers. |
| Tests and verification | `tests/test_package_layout.py`, `tests/test_config_identifiers.py`, `tests/test_global_config_registry.py`, `tests/test_path_resolution.py`, `tests/test_column_inference.py`, `tests/test_chromosome_inference.py`, `tests/test_genome_build_inference.py`, `tests/test_regions.py`, and `tests/test_logging_refactor.py`. | `tests/test_annotation.py`. | `tests/test_ref_panel.py`, `tests/test_ref_panel_builder.py`, `tests/test_r2_query.py`. | `tests/test_ldscore_workflow.py`, `tests/test_plink_io.py`, `tests/test_output.py`. | `tests/test_sumstats_munger.py`. | `tests/test_regression_workflow.py`, `tests/test_irwls.py`, `tests/test_jackknife.py`. | `tests/test_output.py` plus workflow tests that assert written artifacts can be consumed downstream. |

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
unrelated files are preserved. Sharded workflows may narrow the owned family to
the shard selected by the current invocation; `build-ref-panel` does this for
concrete chromosome PLINK prefixes, while `@` chromosome-suite prefixes use the
full panel family.

This family policy applies to `munge-sumstats`, `build-ref-panel`, `ldscore`,
`partitioned-h2`, `rg`, and `annotate`.

Workflow logs follow that same preflight policy. The log path is checked with
the scientific outputs before the log handler opens the file; if computation
fails after that point, the file is kept with a `Failed` lifecycle footer.
Lifecycle start/end timestamps and elapsed duration are sampled from paired
entry/exit timepoints, so the footer describes the same interval as the logged
wall-clock timepoints.
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
- Matplotlib remains an optional extra imported only by private plotting
  builders or h2 sensitivity mode. Plotting consumes canonical saved values and
  never refits a regression.
- Fixed output paths, including workflow logs, should be precomputed and
  preflighted before long-running kernels start, especially for reference-panel
  and sumstats workflows.
