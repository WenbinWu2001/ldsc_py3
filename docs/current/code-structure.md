# Code Structure

Last updated on: 2026-09-15

This is the authoritative contributor entry for the `ldsc` package. Start here to locate a change. [Architecture](architecture.md) explains execution boundaries, [data flow](data-flow.md) explains artifact streams, and [layer structure](layer-structure.md) supplies the detailed ownership matrix. These domain references supplement this entry rather than defining competing navigation maps.

## Repository Map

```text
ldsc_py3_restructured/
├── docs/
├── src/ldsc/
│   ├── __init__.py
│   ├── __main__.py
│   ├── cli.py
│   ├── _cli_help.py
│   ├── config.py
│   ├── _result_files.py
│   ├── path_resolution.py
│   ├── _input_preflight.py
│   ├── _progress.py
│   ├── _logging.py
│   ├── column_inference.py
│   ├── chromosome_inference.py
│   ├── genome_build_inference.py
│   ├── gene_ldscore_index.py
│   ├── gene_list_resolver.py
│   ├── query_annotations.py
│   ├── annotation_builder.py
│   ├── annotation_semantics.py
│   ├── ref_panel_builder.py
│   ├── r2_query.py
│   ├── ldscore_calculator.py
│   ├── _ldscore_preflight.py
│   ├── legacy_ldscore_converter.py
│   ├── _sumstats_input.py
│   ├── sumstats_munger.py
│   ├── regression_runner.py
│   ├── h2_scale.py
│   ├── plotting/
│   ├── quantile_h2.py
│   ├── prevalence.py
│   ├── overlap_matrix.py
│   ├── outputs.py
│   └── _kernel/
├── tests/
└── tutorials/
```

## Dependency Direction

- `ldsc.cli` -> public workflow modules
- public workflow modules -> shared helpers in `config.py`, `path_resolution.py`, `column_inference.py`, `chromosome_inference.py`, `genome_build_inference.py`
- public workflow modules -> private `_kernel` modules
- `ldsc.outputs` is called from the workflow layer, not from `_kernel`
- regression reloads written LD-score artifacts; it does not depend on annotation or reference-panel kernels directly
- `ldsc.plotting` and `ldsc.h2_scale` consume canonical result directories; only the private plotting builders and prevalence-range conversion path lazily import required Matplotlib

## Module Map

| Module | Responsibility |
| --- | --- |
| `ldsc.cli` | unified `ldsc` command and subcommand dispatch |
| `ldsc._cli_help` | shared `CLIHelpFormatter` keeps flag names and paths intact when help descriptions wrap; `CHROMOSOME_PATH_HELP` and `SCALAR_PATH_HELP` supply verified path syntax; workflow parsers own their descriptions and argument groups |
| `ldsc.config` | frozen public config dataclasses and basic validation |
| `ldsc.path_resolution` | normalize path tokens, resolve concrete input files, create output directories, preflight fixed output paths, and enforce coherent output artifact families |
| `ldsc._input_preflight` | aggregate declaration, companion, bounded-header, and saved-manifest checks; retain concrete paths and six-column repair diagnostics; leave scientific content validators with their owners |
| `ldsc._progress` | phase boundaries and time-throttled progress through existing logger routes; bounded counters, current objects, and elapsed time |
| `ldsc._logging` | shared workflow logging context, LDSC logger level handling, lifecycle audit lines, CLI console-handler routing (file-authoritative, console error-only) with run-aborting traceback capture, durable authorized-overwrite failure markers without rollback, and log-only formatting helpers |
| `ldsc.column_inference` | resolve header aliases and normalize identifier/build tokens |
| `ldsc.chromosome_inference` | canonical chromosome normalization and ordering |
| `ldsc.genome_build_inference` | public `chr_pos` build and coordinate-basis inference helpers |
| `ldsc.hm3` | public packaged curated HM3 map loader and installed map path helper for workflow internals |
| `ldsc.gene_list_resolver` | required/embedded catalog validation and vectorized focal/control resolution into exact IDs, intervals, audits, and summaries |
| `ldsc.query_annotations` | shared gene coverage, resolution/support precedence, required-control and all-skipped failures, ordered diagnostics, and aligned result pruning; does not import workflow orchestration |
| `tools/hm3/build_hm3_chr_pos_reference.py` | maintenance tool (outside the package) that rebuilds the compact HM3 coordinate reference used by genome-build inference |
| `ldsc._kernel.liftover` | shared hg19/hg38 liftover helpers, chain-file translation, curated HM3 dual-build coordinate conversion, drop-all coordinate collision helpers, and readable drop reports |
| `ldsc._kernel.regions` | packaged and user BED interval loading plus region-exclusion masks |
| `ldsc._kernel.plink_bed` | PLINK genotype reader (`PlinkBEDFile` and its `__GenotypeArrayInMemory__` base, incl. the in-class LD-score block sums): lazy header read, per-SNP selective decode with fused individual filter, and opt-in disk streaming for unrestricted builds; never materializes the whole-chromosome bitarray |
| `ldsc.annotation_builder` | public annotation builder and general standalone parser/API dispatch seam |
| `ldsc.annotation_semantics` | global annotation-name uniqueness and advisory binary/quantitative classification |
| `ldsc._annotation_bundle`, `ldsc._annotation_storage` | explicit dataset ownership, chromosome descriptors, detached selected reads, and output-contained scratch |
| `ldsc._annotation_outputs` | saved query descriptors with original baseline dependencies; deferred preparation on explicit reads and construction-owner detachment |
| `ldsc._annotation_sources`, `ldsc._annotation_identity`, `ldsc._annotation_loading` | bounded input scans, aligned source preparation, global identity cleanup, and shared source-backed construction |
| `ldsc.annotate_workflow`, `ldsc._annotation_queries`, `ldsc._direct_annotation` | standalone BED/gene gates, incremental projection, and direct LD-score integration |
| `ldsc._gene_query_storage` | staged per-source gene selections and complete audit replay with shared catalog state |
| `ldsc.ldscore_source` | validate ordered query manifests and read explicit selections across LD-score batch files |
| `ldsc._annotation_memory` | small prepared annotations and diagnostics without filesystem writes |
| `ldsc._ldscore_batch_output` | consume/release query batches, stage privately, and publish complete LD directories |
| `ldsc._parallelism` | shared nonzero-integer worker validation, CLI parsing, affinity-aware CPU discovery, and work-capped resolution; tested across direct/indexed scoring, gene-index construction, and partitioned query fitting in `tests/test_worker_policy.py` |
| `ldsc._indexed_ldscore_batches` | bounded chromosome workers holding one operator each and deterministic batch-fragment assembly |
| `ldsc._quantile_inputs`, `ldsc._quantile_storage` | bounded common-universe reconstruction, exact global boundaries, sufficient statistics, and streamed alignment diagnostics |
| `ldsc._kernel.ldscore_projection` | float64 output-row accumulators and query batches sharing each LD block |
| `ldsc.ref_panel_builder` | parquet reference-panel build workflow, including source-build inference, optional coordinate liftover, explicit SNP/sample restrictions, and optional `min_r2` pair-emission threshold |
| `ldsc.r2_query` | public `query-r2` CLI/API, `R2Panel`, one-shot `query_r2()`, sidecar-binding validation, endpoint key resolution, sign harmonization, and optional adjusted-R2-to-Pearson-r conversion |
| `ldsc.ldscore_calculator` | LD-score orchestration, catalog-build selection, direct prepared-SNP support measurement, optional synthetic `base`, aggregation, and output routing through shared query finalization |
| `ldsc._ldscore_preflight` | validates direct annotation/reference contents and chromosome equality; enforces complete `@` declarations, authoritative glob matches, and consolidated required-input diagnostics |
| `ldsc.legacy_ldscore_converter` | sole LDSC2 LD-score-suite import boundary: deterministic family discovery, rsID joins, count/overlap validation or reconstruction, source/discovery provenance without source-file hashes, diagnostics, and canonical LDSC3 directory writing |
| `ldsc.sumstats_munger` | raw-sumstats CLI/API orchestration, `--input-format auto` / `--infer-only` header inference, Parquet/TSV curated output writing, self-describing `sumstats.parquet` footer identity metadata, diagnostics under `diagnostics/`, canonical `CHR`/`POS` sumstats output, and curated sumstats loader |
| `ldsc._sumstats_input` | private workflow helper resolving raw schema, DANER/sample-size settings, bounded source-build evidence and keep-lists into `ResolvedMungeInput` |
| `ldsc._kernel.sumstats_munger` | chunk QC and restriction, whole-table N and sign conversion, resolved liftover and global identity cleanup; returns `MungeResult` with counts and provenance |
| `ldsc.regression_runner` | file-driven regression dataset assembly, automatic legacy LDSC2 sumstats rsID-to-panel projection and allele harmonization, active effective identity-key merging (`SNP`, `SNP:<allele_set>`, `CHR:POS`, or `CHR:POS:<allele_set>`), h2/partitioned-h2/rg estimator dispatch (including the two overlap-aware partitioned-h2 regimes), observed/liability-scale summary columns, exact final-fit h2 regression-bin diagnostics, and rg result-family writing |
| `ldsc._partitioned_h2_parallel`, `ldsc._parallelism` | private read-only regression maps, bounded spawned whole-query workers, worker termination/log forwarding, and shared worker-count resolution; statistical calculations remain in `RegressionRunner._fit_partitioned_query()` |
| `ldsc.quantile_h2` | post-fit continuous-target quantile assignment, fitted-source/common-universe reconstruction and verification, vectorized coefficient/delete-value projection, standardized `tau_star`, CLI orchestration, and diagnostics |
| `ldsc.prevalence` | parse and validate binary-trait prevalence inputs (scalar `--samp-prev`/`--pop-prev` for h2/partitioned-h2; comma-separated lists or a `--prevalence-manifest` TSV for rg) into a normalized per-trait `(samp_prev, pop_prev)` structure for observed-to-liability conversion |
| `ldsc.h2_scale` | strict post-fit observed-to-liability conversion from a canonical h2 result, including exact and prevalence-range modes and the fixed nested derived-result family |
| `ldsc.plotting` | sole public metadata-driven plotting dispatcher, lazy plotting-runtime boundary, fixed plot-family output, and live `PlotArtifact` return object |
| `ldsc.plotting._builders` | private Matplotlib-only headless builders for the approved h2, rg, partitioned-h2, and quantile-h2 plots |
| `ldsc.overlap_matrix` | public-layer overlap container (`LDScoreOverlap`), long-form parquet (de)serialization, per-model overlap assembly, the overlap-aware category table (ported `_overlap_output` + augmentation), and the collinearity hard-error check (`model_collinearity_error`) |
| `ldsc.outputs` | per-writer `artifact_family()` declarations shared by workflow preflight, final writes, metadata file entries, and cleanup; artifact naming, LD-score parquet and query-diagnostic layout, partitioned-h2 per-query layout, rg result-family layout, metadata JSON payloads, and serialization |
| `ldsc._kernel.overlap` | low-level annotation overlap-block computation (`OverlapContribution`, `compute_overlap`, `sum_overlap_contributions`) |
| `ldsc._kernel.annotation` | low-level annotation table reading and BED intersection helpers |
| `ldsc._kernel.ref_panel_builder` | optional genetic-map parsing, optional liftover, parquet schemas, pairwise LD emission |
| `ldsc._kernel.ref_panel` | runtime PLINK/parquet adapters; `RefPanel.prepare_chromosome` owns reference filtering, annotation alignment, window/bias policy, and reader lifetime |
| `ldsc._kernel.r2_query` | low-level index-format parquet pair lookup used by `ldsc.r2_query` |
| `ldsc._kernel.ldscore` | `PreparedChromosome` data contract, LD-score projection from prepared state, streaming R2 reader, annotation/window/count primitives; no LDSC2 artifact emitters |
| `ldsc._kernel.sumstats_munger` | legacy-compatible raw summary-statistics QC, normalization, and optional coordinate liftover; returns `MungeResult` with tables, counts and provenance and owns no output files |
| `ldsc._kernel.regression` | LDSC estimators for `Hsq` and `RG` |
| `ldsc._kernel._jackknife`, `ldsc._kernel._irwls` | supporting numerical routines used by regression |
| `ldsc._kernel.formats`, `ldsc._kernel.identifiers` | retained PLINK/list primitives and SNP identifier helpers; obsolete legacy regression-artifact readers are removed |

## Where To Change Code

| Goal | Start here |
| --- | --- |
| change CLI flags or subcommand wiring | workflow `build_parser` / `add_*_arguments` functions, then `src/ldsc/cli.py` for dispatch and copied help groups; follow [CLI help guidelines](cli-help-guidelines.md) |
| change help wrapping or shared log-level wording | `src/ldsc/_cli_help.py:CLIHelpFormatter`, `src/ldsc/_logging.py:LOG_LEVEL_HELP`; verify `tests/test_cli_help.py` |
| change path-token behavior | `src/ldsc/path_resolution.py` |
| change output collision policy | `src/ldsc/path_resolution.py`, then the workflow writer that owns the artifact |
| change workflow log files, lifecycle lines, or log-level handling | `src/ldsc/_logging.py`, then the workflow orchestration function that owns the output directory |
| change header aliases or identifier/build normalization | `src/ldsc/column_inference.py` |
| change automatic `chr_pos` genome-build inference | `src/ldsc/genome_build_inference.py` |
| change annotation loading, `ldsc annotate` behavior, or BED projection | `src/ldsc/annotation_builder.py`, then `src/ldsc/_kernel/annotation.py` |
| change gene catalog validation, identifier resolution, or list naming | `src/ldsc/gene_list_resolver.py`, then `src/ldsc/annotation_builder.py` |
| change BED/gene partial-success statuses or diagnostic schemas | `src/ldsc/query_annotations.py`, `src/ldsc/ldscore_calculator.py`, `src/ldsc/outputs.py` |
| change parquet reference-panel build logic | `src/ldsc/ref_panel_builder.py`, then `src/ldsc/_kernel/ref_panel_builder.py` |
| change named regression-SNP region exclusion for `ldscore` or gene-index construction | `src/ldsc/config.py`, `src/ldsc/_kernel/regions.py`, then `src/ldsc/ldscore_calculator.py` or `src/ldsc/gene_ldscore_index.py` |
| change R2 pair lookup CLI/API behavior | `src/ldsc/r2_query.py`, then `src/ldsc/_kernel/r2_query.py` |
| change runtime PLINK/parquet reference access | `src/ldsc/_kernel/ref_panel.py` |
| change LD-score orchestration, optional-baseline behavior, or output packaging | `src/ldsc/ldscore_calculator.py`, `src/ldsc/outputs.py` |
| change LD-score math | `src/ldsc/_kernel/ldscore.py` |
| change raw sumstats ingestion, format inference, `CHR`/`POS` handling, sumstats SNP keep-list filtering, liftover drop audit sidecars, footer provenance, or curated loading | `src/ldsc/sumstats_munger.py`, then `src/ldsc/_kernel/sumstats_munger.py` |
| change regression dataset assembly or CLI summaries | `src/ldsc/regression_runner.py`, then `src/ldsc/outputs.py`, `docs/current/partitioned-h2-results.md` for partitioned-h2 output layout, `docs/current/partitioned-ldsc-workflow.md` for rg output contracts, and `docs/current/regression-configuration.md` for the tunable estimator parameters and defaults |
| change annotation-name validation or continuous-annotation classification | `src/ldsc/annotation_semantics.py`, then `src/ldsc/ldscore_calculator.py` and `src/ldsc/outputs.py` |
| change target quantiles, standardized coefficients, or post-fit alignment | `src/ldsc/quantile_h2.py`, then `src/ldsc/outputs.py` and `docs/current/continuous-annotation-quantile-h2.md` |
| change LDSC estimators | `src/ldsc/_kernel/regression.py` |
| change binary-trait prevalence parsing or observed-to-liability conversion | `src/ldsc/prevalence.py` (input parsing/validation), then `src/ldsc/regression_runner.py` (summary/metadata wiring) and `src/ldsc/_kernel/regression.py` (`liability_conversion_factor`) |
| change post-fit h2 conversion or sensitivity output | `src/ldsc/h2_scale.py`, then `src/ldsc/_kernel/regression.py` only if the shared numerical primitive itself changes |
| change plot dispatch, supported result contracts, or visual semantics | `src/ldsc/plotting/__init__.py`, `src/ldsc/plotting/_builders.py`, then `docs/current/plotting-module.md` and `tutorials/plotting-results.md` |
| change owned artifact declarations or LD-score result-directory files, parquet row-group layout, partitioned-h2 per-query layout, rg result-family layout, or metadata JSON payloads | `src/ldsc/outputs.py` |

## Architectural Rules That Matter In Practice

- Treat `src/ldsc/` as the only supported Python import surface.
- Do not add user-facing path discovery to `_kernel`; pass concrete files in.
- Keep public file contracts for `.annot(.gz)`, self-describing Parquet munged sumstats (footer identity metadata) plus optional `.sumstats.gz` compatibility output, canonical LD-score result directories, and regression summary directories stable unless the change is intentional and coordinated. LD-score parquet files remain flat files, with chromosome-aligned row groups documented through root `metadata.json`; diagnostic logs, dropped-SNP reports, and diagnostic metadata live under `diagnostics/`. Legacy reference/weight `.l2.ldscore(.gz)`, `.l2.M`, and `.l2.M_5_50` families enter through the explicit [LD-score converter](legacy-ldscore-conversion.md#compatibility-boundary). The `.w.l2.ldscore(.gz)` suffix is unsupported, and private kernel legacy artifact emitters and regression readers have been removed.
- Treat package-built index-format R2 panels as the only public R2 parquet
  format. Both `ldsc ldscore --r2-dir` and `ldsc query-r2` depend on the paired
  `chr*_r2.parquet` / `chr*_meta.tsv.gz` layout and sidecar identity binding.
- Keep optional-baseline behavior in the public LD-score workflow layer: no baseline and no query means a synthetic all-ones `base`; prebuilt, BED, and gene-list queries require explicit baseline annotations. Skipped BED/gene queries never become placeholder scientific columns.
- Keep `ldsc annotate` orchestration in `ldsc.annotation_builder`. The CLI
  registers annotation flags from that module and dispatches parsed namespaces
  to `run_annotate_from_args()`; `_kernel.annotation` must not own parser
  functions, result objects, or public workflow aliases.
- Use the shared output preflight helpers for fixed output files, including
  workflow log files. Public workflows should create missing output
  directories, reuse existing directories, fail on existing owned artifacts by
  default, and require `--overwrite` or `overwrite=True` for replacement.
- Keep coherent output families consistent. For `munge-sumstats`,
  `build-r2-panel`, `ldscore`, `partitioned-h2`, `rg`, and `annotate`,
  no-overwrite mode rejects any current-contract owned sibling from the workflow family.
  Overwrite mode writes the requested current outputs and then removes stale
  current-contract owned siblings that were not produced by the successful run.
  Removed legacy root diagnostic names are not workflow-owned. Unrelated files
  in the output directory must be preserved. Sharded workflows may pass a
  narrowed owned family when shards share an output directory; `build-r2-panel`
  concrete chromosome prefixes own only that chromosome's package, while `@`
  chromosome-suite prefixes own the full panel package.
- Preflight deterministic output paths before expensive or multi-file writes.
  This is especially important for `build-r2-panel`, `munge-sumstats`,
  `annotate`, and summary-table regression commands.
- Treat workflow logs as audit artifacts. Do not include log paths in
  `output_paths` mappings or thin metadata sidecars that downstream code
  interprets as scientific data artifacts. Whole-directory transactional
  publishers must keep open workflow logs outside the replaceable tree; the
  gene-index builder uses hidden sibling `.<index-name>.build-state/` and moves
  only the closed successful log into the published diagnostics.
- Keep regression file-driven: it should be able to rebuild state from written artifacts without recomputing LD scores.
- Keep plotting metadata-driven and post-fit: plot builders consume declared canonical tables and never refit. Matplotlib is required at installation but lazily imported through `ldsc.plotting._builders`.
- Prefer extending shared helpers or the workflow-owned writer over duplicating local parsing or writing logic.

## Test Map

| Area | Main tests |
| --- | --- |
| output layer | `tests/test_output.py` |
| annotation workflow | `tests/test_annotation.py` |
| reference-panel builder | `tests/test_ref_panel_builder.py` |
| R2 pair query | `tests/test_r2_query.py` |
| LD-score workflow | `tests/test_ldscore_workflow.py` |
| gene-list catalog and resolution | `tests/test_gene_list_resolver.py` |
| sumstats munging | `tests/test_sumstats_munger.py` |
| regression workflow | `tests/test_regression_workflow.py` |
| path and config contracts | `tests/test_path_resolution.py`, `tests/test_config_identifiers.py`, `tests/test_column_inference.py`, `tests/test_genome_build_inference.py` |
| region exclusion | `tests/test_regions.py` |

`AnnotationBundle.validate()` validates row alignment and annotation names; identity cleanup occurs during builder construction. It no longer accepts an ignored identifier argument. Column-name inference uses `infer_chr_pos_columns`; numerical window construction uses `get_block_lefts`. Their former forwarding spellings are removed without changing accepted input column aliases or legacy formats.

`ldsc._result_files` owns contained declared-result input paths and atomic JSON publication for `h2_scale` and `plotting`. Each workflow retains its metadata interpretation, table loading and scientific dispatch. `outputs` keeps its separate JSON serialization of NumPy/dataclass payloads.

## Domain contracts and decision history

| Change | Current contract | Main owner |
| --- | --- | --- |
| annotations and interval projection | [BED input](bed-input-format.md), [gene lists](gene-list-input-format.md), [query repair](gene-list-diagnostics-and-repair.md) | `annotation_builder.AnnotationBuilder`; `_kernel.annotation` projects resolved intervals |
| panel preparation, CM/MAF and traversal | [R2 pipeline](parquet-r2-format-and-read-pipeline.md), [pair accumulation](ldscore-parquet-accumulation.md), [window behavior](ld-window-parquet-r2-sidecar-behavior.md) | `_kernel.ref_panel.RefPanel.prepare_chromosome`; `_kernel.ldscore.compute_chromosome` |
| regression rows and region masks | [SNP universe](ldscore-snp-universe-contract.md), [region exclusions](region-exclusion-presets.md) | `ldscore_calculator._regression_region_intervals`; `_kernel.ldscore.regression_mask_from_keys`; `gene_ldscore_index` resolves indexed policy |
| reference construction and pair queries | [pair query](ref-panel-r2-query.md), [path contract](path-specification.md) | `ref_panel_builder.ReferencePanelBuilder`; `r2_query.R2Panel`; `outputs.QueryR2DirectoryWriter` |
| munging and liftover | [munging](munge-sumstats.md), [liftover](liftover-harmonization-decisions.md) | `sumstats_munger.SumstatsMunger`; `_sumstats_input.prepare_munge_input`; `_kernel.sumstats_munger` |
| regression and overlap | [configuration](regression-configuration.md), [partitioned results](partitioned-h2-results.md) | `regression_runner.RegressionRunner`; `overlap_matrix`; `_kernel.regression` |
| indexed gene queries | [index](gene-ldscore-index.md), [mathematics](gene-ldscore-index-mathematics.md) | `gene_ldscore_index`; `_kernel.gene_ldscore_index` |
| derived results | [plotting](plotting-module.md), [quantile projection](continuous-annotation-quantile-h2.md) | `plotting`, `h2_scale`, `quantile_h2`; shared `_result_files` |

`query-r2` requires `--output-dir`; `QueryR2DirectoryWriter` writes its result directory. Region presets affect regression rows and weights, not the LD-reference contributor universe. `build-r2-panel` has no region-pruning fields; use its explicit reference-SNP restriction when pruning is intended. Reference preparation supplies population-specific CM/MAF; annotation CM is an `NA` placeholder. Parquet LD scores use pair streaming, and PLINK filtering happens during adapter preparation.

Historical decisions formerly indexed by root `design_map.md` remain available below. They record decisions at their dates; current contracts above govern current work. The [old map](../archive/design/2026-09-10-design-map.md) and [completed threshold audit](../audits/threshold-comparison/README.md) preserve prior evidence without competing with this entry.

- [2026-05-02-logging-harmonization-implementation-plan](../plans/2026-05-02-logging-harmonization-implementation-plan.md)
- [2026-05-09-batch-rg-implementation-plan](../plans/2026-05-09-batch-rg-implementation-plan.md)
- [2026-05-10-liftover-harmonization](../plans/2026-05-10-liftover-harmonization.md)
- [2026-05-11-ref-panel-stale-class2-warning](../plans/2026-05-11-ref-panel-stale-class2-warning.md)
- [2026-06-04-logging-console-file-routing-plan](../plans/2026-06-04-logging-console-file-routing-plan.md)
- [2026-06-05-ldscore-chromosome-parallelism-plan](../plans/2026-06-05-ldscore-chromosome-parallelism-plan.md)
- [2026-06-06-ref-panel-r2-query-plan](../plans/2026-06-06-ref-panel-r2-query-plan.md)
- [2026-06-06-sumstats-self-describing-artifact-plan](../plans/2026-06-06-sumstats-self-describing-artifact-plan.md)
- [2026-06-11-cm-maf-source-of-truth-plan](../plans/2026-06-11-cm-maf-source-of-truth-plan.md)
- [2026-06-11-overlap-aware-partitioned-h2-plan](../plans/2026-06-11-overlap-aware-partitioned-h2-plan.md)
- [2026-05-02-logging-harmonization-design](../specs/2026-05-02-logging-harmonization-design.md)
- [2026-05-09-batch-rg-design](../specs/2026-05-09-batch-rg-design.md)
- [2026-05-10-liftover-harmonization](../specs/2026-05-10-liftover-harmonization.md)
- [2026-05-11-ref-panel-stale-class2-warning](../specs/2026-05-11-ref-panel-stale-class2-warning.md)
- [2026-06-01-build-ref-panel-memory-optimization-design](../specs/2026-06-01-build-ref-panel-memory-optimization-design.md)
- [2026-06-02-build-ref-panel-reader-streaming-design](../specs/2026-06-02-build-ref-panel-reader-streaming-design.md)
- [2026-06-03-pair-emission-columnar-and-doc-sync-design](../specs/2026-06-03-pair-emission-columnar-and-doc-sync-design.md)
- [2026-06-05-ldscore-chromosome-parallelism-design](../specs/2026-06-05-ldscore-chromosome-parallelism-design.md)
- [2026-06-06-ref-panel-r2-query-design](../specs/2026-06-06-ref-panel-r2-query-design.md)
- [2026-06-06-region-exclusion-design](../specs/2026-06-06-region-exclusion-design.md)
- [2026-06-06-sumstats-self-describing-artifact-design](../specs/2026-06-06-sumstats-self-describing-artifact-design.md)
- [2026-06-11-cm-maf-source-of-truth-design](../specs/2026-06-11-cm-maf-source-of-truth-design.md)
- [2026-06-11-overlap-aware-partitioned-h2-design](../specs/2026-06-11-overlap-aware-partitioned-h2-design.md)
