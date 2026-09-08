# Code Structure

Last updated on: 2026-09-08

This is the contributor-facing module map for `ldsc_py3_Jerry`.

## Repository Map

```text
ldsc_py3_Jerry/
├── docs/
├── src/ldsc/
│   ├── __init__.py
│   ├── __main__.py
│   ├── cli.py
│   ├── config.py
│   ├── path_resolution.py
│   ├── _logging.py
│   ├── column_inference.py
│   ├── chromosome_inference.py
│   ├── genome_build_inference.py
│   ├── gene_list_resolver.py
│   ├── query_annotations.py
│   ├── annotation_builder.py
│   ├── annotation_semantics.py
│   ├── ref_panel_builder.py
│   ├── r2_query.py
│   ├── ldscore_calculator.py
│   ├── legacy_ldscore_converter.py
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
| `ldsc.config` | frozen public config dataclasses and basic validation |
| `ldsc.path_resolution` | normalize path tokens, resolve concrete input files, create output directories, preflight fixed output paths, and enforce coherent output artifact families |
| `ldsc._logging` | shared workflow logging context, LDSC logger level handling, lifecycle audit lines, CLI console-handler routing (file-authoritative, console error-only) with run-aborting traceback capture, durable authorized-overwrite failure markers without rollback, and log-only formatting helpers |
| `ldsc.column_inference` | resolve header aliases and normalize identifier/build tokens |
| `ldsc.chromosome_inference` | canonical chromosome normalization and ordering |
| `ldsc.genome_build_inference` | public `chr_pos` build and coordinate-basis inference helpers |
| `ldsc.hm3` | public packaged curated HM3 map loader and installed map path helper for workflow internals |
| `ldsc.gene_list_resolver` | required/embedded catalog validation and vectorized focal/control resolution into exact IDs, intervals, audits, and summaries |
| `ldsc.query_annotations` | internal ordered BED/gene query status record shared by annotation, LD-score, and output layers |
| `tools/hm3/build_hm3_chr_pos_reference.py` | maintenance tool (outside the package) that rebuilds the compact HM3 coordinate reference used by genome-build inference |
| `ldsc._kernel.liftover` | shared hg19/hg38 liftover helpers, chain-file translation, curated HM3 dual-build coordinate conversion, drop-all coordinate collision helpers, and readable drop reports |
| `ldsc._kernel.regions` | packaged and user BED interval loading plus region-exclusion masks |
| `ldsc._kernel.plink_bed` | PLINK genotype reader (`PlinkBEDFile` and its `__GenotypeArrayInMemory__` base, incl. the in-class LD-score block sums): lazy header read, per-SNP selective decode with fused individual filter, and opt-in disk streaming for unrestricted builds; never materializes the whole-chromosome bitarray |
| `ldsc.annotation_builder` | public annotation workflow: path resolution, bundle loading, BED/gene interval projection, query-local source isolation, and query `.annot.gz` writing for `annotate` |
| `ldsc.annotation_semantics` | global annotation-name uniqueness, binary/quantitative classification, and compact common-universe SHA256 semantic fingerprints |
| `ldsc.ref_panel_builder` | parquet reference-panel build workflow, including source-build inference, optional coordinate liftover, explicit SNP/sample restrictions, and optional `min_r2` pair-emission threshold |
| `ldsc.r2_query` | public `query-r2` CLI/API, `R2Panel`, one-shot `query_r2()`, sidecar-binding validation, endpoint key resolution, sign harmonization, and optional adjusted-R2-to-Pearson-r conversion |
| `ldsc.ldscore_calculator` | LD-score orchestration, catalog-build selection, optional synthetic `base`, query-status finalization/pruning, aggregation, and output routing |
| `ldsc.legacy_ldscore_converter` | sole LDSC2 LD-score-suite import boundary: deterministic family discovery, rsID joins, count/overlap validation or reconstruction, provenance hashing, diagnostics, and canonical LDSC3 directory writing |
| `ldsc.sumstats_munger` | raw-sumstats CLI/API orchestration, `--format auto` / `--infer-only` header inference, Parquet/TSV curated output writing, self-describing `sumstats.parquet` footer identity metadata, diagnostics under `diagnostics/`, canonical `CHR`/`POS` sumstats output, and curated sumstats loader |
| `ldsc.regression_runner` | file-driven regression dataset assembly, automatic legacy LDSC2 sumstats rsID-to-panel projection and allele harmonization, active effective identity-key merging (`SNP`, `SNP:<allele_set>`, `CHR:POS`, or `CHR:POS:<allele_set>`), h2/partitioned-h2/rg estimator dispatch (including the two overlap-aware partitioned-h2 regimes), observed/liability-scale summary columns, exact final-fit h2 regression-bin diagnostics, and rg result-family writing |
| `ldsc.quantile_h2` | post-fit continuous-target quantile assignment, fitted-source/common-universe reconstruction and verification, vectorized coefficient/delete-value projection, standardized `tau_star`, CLI orchestration, and diagnostics |
| `ldsc.prevalence` | parse and validate binary-trait prevalence inputs (scalar `--samp-prev`/`--pop-prev` for h2/partitioned-h2; comma-separated lists or a `--prevalence-manifest` TSV for rg) into a normalized per-trait `(samp_prev, pop_prev)` structure for observed-to-liability conversion |
| `ldsc.h2_scale` | strict post-fit observed-to-liability conversion from a canonical h2 result, including exact and prevalence-range modes and the fixed nested derived-result family |
| `ldsc.plotting` | sole public metadata-driven plotting dispatcher, lazy plotting-runtime boundary, fixed plot-family output, and live `PlotArtifact` return object |
| `ldsc.plotting._builders` | private Matplotlib-only headless builders for the approved h2, rg, partitioned-h2, and quantile-h2 plots |
| `ldsc.overlap_matrix` | public-layer overlap container (`LDScoreOverlap`), long-form parquet (de)serialization, per-model overlap assembly, the overlap-aware category table (ported `_overlap_output` + augmentation), and the collinearity hard-error check (`model_collinearity_error`) |
| `ldsc.outputs` | artifact naming, LD-score parquet and query-diagnostic layout, partitioned-h2 per-query layout, rg result-family layout, metadata JSON payloads, and serialization |
| `ldsc._kernel.overlap` | low-level annotation overlap-block computation (`OverlapContribution`, `compute_overlap`, `sum_overlap_contributions`) |
| `ldsc._kernel.annotation` | low-level annotation table reading and BED intersection helpers |
| `ldsc._kernel.ref_panel_builder` | optional genetic-map parsing, optional liftover, parquet schemas, pairwise LD emission |
| `ldsc._kernel.ref_panel` | runtime PLINK/parquet reference-panel adapters |
| `ldsc._kernel.r2_query` | low-level index-format parquet pair lookup used by `ldsc.r2_query` |
| `ldsc._kernel.ldscore` | LD-score math, PLINK/parquet readers, chromosome computation, and count primitives; no LDSC2 artifact emitters |
| `ldsc._kernel.sumstats_munger` | legacy-compatible raw summary-statistics QC, normalization, and optional coordinate liftover; returns in-memory tables and owns no output files |
| `ldsc._kernel.regression` | LDSC estimators for `Hsq` and `RG` |
| `ldsc._kernel._jackknife`, `ldsc._kernel._irwls` | supporting numerical routines used by regression |
| `ldsc._kernel.formats`, `ldsc._kernel.identifiers` | retained PLINK/list primitives and SNP identifier helpers; obsolete legacy regression-artifact readers are removed |

## Where To Change Code

| Goal | Start here |
| --- | --- |
| change CLI flags or subcommand wiring | `src/ldsc/cli.py` |
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
| change raw sumstats ingestion, format inference, `CHR`/`POS` handling, sumstats SNP keep-list filtering, liftover drop audit sidecars, sidecar provenance, or curated loading | `src/ldsc/sumstats_munger.py`, then `src/ldsc/_kernel/sumstats_munger.py` |
| change regression dataset assembly or CLI summaries | `src/ldsc/regression_runner.py`, then `src/ldsc/outputs.py`, `docs/current/partitioned-h2-results.md` for partitioned-h2 output layout, `docs/current/partitioned-ldsc-workflow.md` for rg output contracts, and `docs/current/regression-configuration.md` for the tunable estimator parameters and defaults |
| change continuous-annotation classification or fingerprints | `src/ldsc/annotation_semantics.py`, then `src/ldsc/ldscore_calculator.py` and `src/ldsc/outputs.py` |
| change target quantiles, standardized coefficients, or post-fit alignment | `src/ldsc/quantile_h2.py`, then `src/ldsc/outputs.py` and `docs/current/continuous-annotation-quantile-h2.md` |
| change LDSC estimators | `src/ldsc/_kernel/regression.py` |
| change binary-trait prevalence parsing or observed-to-liability conversion | `src/ldsc/prevalence.py` (input parsing/validation), then `src/ldsc/regression_runner.py` (summary/metadata wiring) and `src/ldsc/_kernel/regression.py` (`liability_conversion_factor`) |
| change post-fit h2 conversion or sensitivity output | `src/ldsc/h2_scale.py`, then `src/ldsc/_kernel/regression.py` only if the shared numerical primitive itself changes |
| change plot dispatch, supported result contracts, or visual semantics | `src/ldsc/plotting/__init__.py`, `src/ldsc/plotting/_builders.py`, then `docs/current/plotting-module.md` and `tutorials/plotting-results.md` |
| change LD-score result-directory files, parquet row-group layout, partitioned-h2 per-query layout, rg result-family layout, or metadata JSON payloads | `src/ldsc/outputs.py` |

## Architectural Rules That Matter In Practice

- Treat `src/ldsc/` as the only supported Python import surface.
- Do not add user-facing path discovery to `_kernel`; pass concrete files in.
- Keep public file contracts for `.annot(.gz)`, self-describing Parquet munged sumstats (footer identity metadata) plus optional `.sumstats.gz` compatibility output, canonical LD-score result directories, and regression summary directories stable unless the change is intentional and coordinated. LD-score parquet files remain flat files, with chromosome-aligned row groups documented through root `metadata.json`; diagnostic logs, dropped-SNP reports, and diagnostic metadata live under `diagnostics/`. Legacy `.l2.ldscore(.gz)`, `.w.l2.ldscore(.gz)`, `.l2.M`, and `.l2.M_5_50` files are compatibility concerns rather than the public LD-score output surface.
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
  `build-ref-panel`, `ldscore`, `partitioned-h2`, `rg`, and `annotate`,
  no-overwrite mode rejects any current-contract owned sibling from the workflow family.
  Overwrite mode writes the requested current outputs and then removes stale
  current-contract owned siblings that were not produced by the successful run.
  Removed legacy root diagnostic names are not workflow-owned. Unrelated files
  in the output directory must be preserved. Sharded workflows may pass a
  narrowed owned family when shards share an output directory; `build-ref-panel`
  concrete chromosome prefixes own only that chromosome's package, while `@`
  chromosome-suite prefixes own the full panel package.
- Preflight deterministic output paths before expensive or multi-file writes.
  This is especially important for `build-ref-panel`, `munge-sumstats`,
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
