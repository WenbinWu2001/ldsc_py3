# Code Structure

Last updated on: 2026-10-27

This is the authoritative contributor entry for the `ldsc` package. Start here to locate a change. [Architecture](architecture.md) explains execution boundaries, [data flow](data-flow.md) explains artifact streams, and the [ownership matrix](architecture.md#ownership-by-layer-and-workflow) shows which layer owns each object. These domain references supplement this entry rather than defining competing navigation maps.

For resource-sensitive changes, start with the [memory and CPU efficiency design guide](design/memory-and-cpu.md).
It catalogs implemented storage, streaming, reuse, and numerical-buffer designs, with owners, precision/lifetime invariants, tradeoffs, tests, and remaining large allocations.

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
│   ├── ldscore_source.py
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
└── tools/
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
| `ldsc._kernel.plink_bed` | PLINK genotype reader (`PlinkBEDFile` and its `__GenotypeArrayInMemory__` base, including LD-score block sums): selective SNP reads with fused individual filtering retain only selected packed genotypes; unrestricted panel builds opt into disk streaming instead |
| `ldsc.annotation_builder` | public annotation builder and general standalone parser/API dispatch seam |
| `ldsc.annotation_semantics` | global annotation-name uniqueness and advisory binary/quantitative classification |
| `ldsc._annotation_bundle`, `ldsc._annotation_storage` | explicit dataset ownership, chromosome descriptors, detached selected reads, and output-contained scratch |
| `ldsc._annotation_outputs` | saved baseline/query descriptors with no original-input dependencies; deferred preparation on explicit reads and construction-owner detachment |
| `ldsc._annotation_bim` | metadata-only standalone BIM grid scans using shared identity cleanup; no genotype companions or value synthesis |
| `ldsc._annotation_sources`, `ldsc._annotation_identity`, `ldsc._annotation_loading` | bounded input scans, aligned source preparation, global identity cleanup, and shared source-backed construction |
| `ldsc.annotate_workflow`, `ldsc._annotation_queries`, `ldsc._direct_annotation` | standalone BED/gene gates, incremental projection, and direct LD-score integration |
| `ldsc._gene_query_storage` | staged per-source gene selections and complete audit replay with shared catalog state |
| `ldsc._annotation_memory` | small prepared annotations and diagnostics without filesystem writes |
| `ldsc._ldscore_batch_output` | consume/release materialized or staged query batches, retain explicit shared publication state, and publish complete LD directories |
| `ldsc._ldscore_fragments` | describe producer-owned chromosome fragments and validate/write one query frame at a time against the shared baseline |
| `ldsc._parallelism` | shared `workers` validation, joblib detected CPU budgets, task caps, and execution summaries |
| `ldsc._worker_runtime`, `ldsc._serial_calculation` | fixed 1-BLAS-thread child policy, runtime verification, independent Arrow limits, and controlled standalone numerical calls |
| `ldsc._ldscore_diagnostics` | chromosome logs, batch context, and accurate task/pool failure reporting |
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
| `ldsc.ldscore_source` | canonical directory loading, manifest/schema and table/identity validation shared with writers, resident baseline/count/overlap metadata, and uncached selected-query reads |
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
| `ldsc._kernel.annotation` | concrete annotation table reading and numeric-value validation |
| `ldsc._annotation_projection` | `ChromosomeProjector` maps BED/gene interval unions to chromosome SNP positions; called by `_annotation_queries.build_query_shards()` after shared BED parsing |
| `ldsc._kernel.ref_panel_builder` | optional genetic-map parsing, optional liftover, parquet schemas, pairwise LD emission |
| `ldsc._kernel.ref_panel` | runtime PLINK/parquet adapters; `RefPanel.prepare_chromosome` owns reference filtering, annotation alignment, window/bias policy, and reader lifetime |
| `ldsc._kernel.r2_query` | low-level index-format parquet pair lookup used by `ldsc.r2_query` |
| `ldsc._kernel.ldscore` | `PreparedChromosome` data contract, LD-score projection from prepared state, streaming R2 reader, annotation/window/count primitives; no LDSC2 artifact emitters |
| `ldsc._kernel.regression` | LDSC estimators for `Hsq` and `RG`; `Hsq.summary`, `Gencov.summary`, and `RG.summary` are kept deliberately because `estimate_h2`/`estimate_rg` return these estimators to Python users, even though no internal caller uses the summaries |
| `ldsc._kernel._jackknife`, `ldsc._kernel._irwls` | supporting numerical routines used by regression |
| `ldsc._kernel.formats`, `ldsc._kernel.identifiers` | retained PLINK/list primitives and SNP identifier helpers; obsolete legacy regression-artifact readers are removed |

## Where To Change Code

| Goal | Start here |
| --- | --- |
| change CLI flags or subcommand wiring | workflow `build_parser` / `add_*_arguments` functions, then `src/ldsc/cli.py` for dispatch and copied help groups; follow [CLI help guidelines](guidelines/cli-help.md) |
| change help wrapping or shared log-level wording | `src/ldsc/_cli_help.py:CLIHelpFormatter`, `src/ldsc/_logging.py:LOG_LEVEL_HELP`; verify `tests/test_cli_help.py` |
| change path-token behavior | `src/ldsc/path_resolution.py` |
| change output collision policy | `src/ldsc/path_resolution.py`, then the workflow writer that owns the artifact |
| change workflow log files, lifecycle lines, or log-level handling | `src/ldsc/_logging.py`, then the workflow orchestration function that owns the output directory |
| change header aliases or identifier/build normalization | `src/ldsc/column_inference.py` |
| change read-only artifact build inspection | `src/ldsc/_genome_build_inspection.py`, with CLI dispatch in `src/ldsc/cli.py` and coverage in `tests/test_genome_build_inspection.py` |
| change automatic `chr_pos` genome-build inference | `src/ldsc/genome_build_inference.py` |
| change annotation loading or `ldsc annotate` behavior | `src/ldsc/annotation_builder.py`, then `src/ldsc/_annotation_sources.py` and `src/ldsc/annotate_workflow.py` |
| change BED parsing or interval projection | `src/ldsc/_kernel/regions.py:iter_bed_rows`, then `src/ldsc/_annotation_queries.py` and `src/ldsc/_annotation_projection.py:ChromosomeProjector` |
| change gene catalog validation, identifier resolution, or list naming | `src/ldsc/gene_list_resolver.py`, then `src/ldsc/annotation_builder.py` |
| change BED/gene partial-success statuses or diagnostic schemas | `src/ldsc/query_annotations.py`, `src/ldsc/ldscore_calculator.py`, `src/ldsc/outputs.py` |
| change parquet reference-panel build logic | `src/ldsc/ref_panel_builder.py`, then `src/ldsc/_kernel/ref_panel_builder.py` |
| change named regression-SNP region exclusion for `ldscore` or gene-index construction | `src/ldsc/config.py`, `src/ldsc/_kernel/regions.py`, then `src/ldsc/ldscore_calculator.py` or `src/ldsc/gene_ldscore_index.py` |
| change R2 pair lookup CLI/API behavior | `src/ldsc/r2_query.py`, then `src/ldsc/_kernel/r2_query.py` |
| change runtime PLINK/parquet reference access | `src/ldsc/_kernel/ref_panel.py` |
| change LD-score orchestration, optional-baseline behavior, or output packaging | `src/ldsc/ldscore_calculator.py`, `src/ldsc/outputs.py` |
| change saved LD-score loading or selected-query reads | `src/ldsc/ldscore_source.py`; keep model-specific validation in `src/ldsc/regression_runner.py` |
| change LD-score math | `src/ldsc/_kernel/ldscore.py` |
| change raw sumstats ingestion, format inference, `CHR`/`POS` handling, sumstats SNP keep-list filtering, liftover drop audit sidecars, footer provenance, or curated loading | `src/ldsc/sumstats_munger.py`, then `src/ldsc/_kernel/sumstats_munger.py` |
| change regression dataset assembly or CLI summaries | `src/ldsc/regression_runner.py`, then `src/ldsc/outputs.py`, `docs/reference/results/partitioned-h2.md` for partitioned-h2 output layout, `docs/reference/formats/artifact-metadata.md` for rg output contracts, and `docs/reference/regression-options.md` for the tunable estimator parameters and defaults |
| change annotation-name validation or continuous-annotation classification | `src/ldsc/annotation_semantics.py`, then `src/ldsc/ldscore_calculator.py` and `src/ldsc/outputs.py` |
| change target quantiles, standardized coefficients, or post-fit alignment | `src/ldsc/quantile_h2.py`, then `src/ldsc/outputs.py` and `docs/reference/results/quantile-h2.md` |
| change LDSC estimators | `src/ldsc/_kernel/regression.py` |
| change binary-trait prevalence parsing or observed-to-liability conversion | `src/ldsc/prevalence.py` (input parsing/validation), then `src/ldsc/regression_runner.py` (summary/metadata wiring) and `src/ldsc/_kernel/regression.py` (`liability_conversion_factor`) |
| change post-fit h2 conversion or sensitivity output | `src/ldsc/h2_scale.py`, then `src/ldsc/_kernel/regression.py` only if the shared numerical primitive itself changes |
| change plot dispatch, supported result contracts, or visual semantics | `src/ldsc/plotting/__init__.py`, `src/ldsc/plotting/_builders.py`, then `docs/dev/design/plotting.md` and `docs/user-guide/workflows/plotting.md` |
| change owned artifact declarations or LD-score result-directory files, parquet row-group layout, partitioned-h2 per-query layout, rg result-family layout, or metadata JSON payloads | `src/ldsc/outputs.py` |

## Architectural Rules That Matter In Practice

- Treat `src/ldsc/` as the only supported Python import surface.
- Do not add user-facing path discovery to `_kernel`; pass concrete files in.
- Keep public file contracts for `.annot(.gz)`, self-describing Parquet munged sumstats (footer identity metadata) plus optional `.sumstats.gz` compatibility output, canonical LD-score result directories, and regression summary directories stable unless the change is intentional and coordinated. LD-score parquet files remain flat files, with chromosome-aligned row groups documented through root `metadata.json`; diagnostic logs, dropped-SNP reports, and diagnostic metadata live under `diagnostics/`. Legacy reference/weight `.l2.ldscore(.gz)`, `.l2.M`, and `.l2.M_5_50` families enter through the explicit [LD-score converter](../reference/formats/legacy-ldscores.md#compatibility-boundary). The `.w.l2.ldscore(.gz)` suffix is unsupported, and private kernel legacy artifact emitters and regression readers have been removed.
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
| annotations and interval projection | [BED input](../reference/formats/bed-intervals.md), [gene lists](../reference/formats/gene-lists.md), [query repair](../reference/results/gene-list-diagnostics.md) | `annotation_builder.AnnotationBuilder`; `_annotation_queries.build_query_shards`; `_annotation_projection.ChromosomeProjector` |
| panel preparation, CM/MAF and traversal | [R2 pipeline](../reference/formats/r2-panels.md), [pair accumulation](../reference/methods/ld-score-computation.md), [window behavior](../reference/formats/r2-panels.md#37-ld-windows-and-sidecar-values-in-ldscore) | `_kernel.ref_panel.RefPanel.prepare_chromosome`; `_kernel.ldscore.compute_chromosome` |
| regression rows and region masks | [SNP universe](../reference/ldscore-snp-universes.md), [region exclusions](../reference/formats/region-exclusion-presets.md) | `ldscore_calculator._regression_region_intervals`; `_kernel.ldscore.regression_mask_from_keys`; `gene_ldscore_index` resolves indexed policy |
| reference construction and pair queries | [pair query](../reference/query-r2.md), [path contract](../reference/input-paths.md) | `ref_panel_builder.ReferencePanelBuilder`; `r2_query.R2Panel`; `outputs.QueryR2DirectoryWriter` |
| munging and liftover | [munging](../reference/munge-sumstats.md), [liftover](../reference/snp-identity-and-genome-build.md#liftover) | `sumstats_munger.SumstatsMunger`; `_sumstats_input.prepare_munge_input`; `_kernel.sumstats_munger` |
| regression and overlap | [configuration](../reference/regression-options.md), [partitioned results](../reference/results/partitioned-h2.md) | `regression_runner.RegressionRunner`; `overlap_matrix`; `_kernel.regression` |
| indexed gene queries | [index](../reference/gene-ldscore-index.md), [mathematics](../reference/methods/gene-ldscore-index-math.md) | `gene_ldscore_index`; `_kernel.gene_ldscore_index` |
| derived results | [plotting](design/plotting.md), [quantile projection](../reference/results/quantile-h2.md) | `plotting`, `h2_scale`, `quantile_h2`; shared `_result_files` |

`query-r2` requires `--output-dir`; `QueryR2DirectoryWriter` writes its result directory. Region presets affect regression rows and weights, not the LD-reference contributor universe. `build-r2-panel` has no region-pruning fields; use its explicit reference-SNP restriction when pruning is intended. Reference preparation supplies population-specific CM/MAF; annotation CM is an `NA` placeholder. Parquet LD scores use pair streaming, and PLINK filtering happens during adapter preparation.

Design rationale that protects a non-obvious constraint is stated in the page that documents the constraint; measurements are in [validation](validation/performance.md).
Earlier dated specifications, implementation plans, audits, and the archived design map were removed on 2026-10-26, and the decision records on 2026-10-27; they remain in git history, for example `git show 479fc4d:docs/specs/<file>` or `git show 77a54e2:docs/decisions/2026-10-11-parallel-execution.md`.
