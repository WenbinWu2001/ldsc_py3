# Code Structure

Last updated on: 2026-10-30

This is the authoritative contributor entry for LDSC3.
Start here to locate a change.
[Architecture](architecture.md) explains execution boundaries, [data flow](data-flow.md) explains artifact streams, and the [ownership matrix](architecture.md#ownership-by-layer-and-workflow) shows which layer owns each object.
These domain references supplement this entry rather than defining competing navigation maps.

For resource-sensitive changes, start with the [memory and CPU efficiency design guide](design/memory-and-cpu.md).
It catalogs implemented storage, streaming, reuse, and numerical-buffer designs, with owners, precision/lifetime invariants, tradeoffs, tests, and remaining large allocations.

## Repository Map

```text
ldsc_py3_restructured/
├── docs/
├── src/ldsc/
│   ├── __init__.py  __main__.py  cli.py  config.py  errors.py   root shell
│   ├── data/                 packaged HM3 maps and region presets
│   ├── _core/                tier 1: chromosomes, columns, coordinates, SNP identity, regions, PLINK companion files
│   ├── _runtime/             tier 2: paths, input preflight, logging, progress, workers, keyword checks, result-file writers
│   ├── _genome/              tier 4: genome-build inference and policy, build inspection, liftover, HM3 map
│   ├── _annotation/          tier 5: annotation files, gene lists, interval projection, overlap, ldsc annotate
│   ├── _refpanel/            tier 5: PLINK reader, R² matrix format, LD windows, runtime adapters, build-r2-matrix, query-r2-matrix
│   ├── _sumstats/            tier 5: munge-sumstats and the munged-table loader
│   ├── _ldscore/             tier 6: LD-score math and workflow, LD-score directories, legacy conversion
│   │   └── gene_index/       exact gene LD-score indexes
│   ├── _regression/          tier 7: estimators, h2/partitioned-h2/rg, liability scale, quantile-h2
│   └── _plotting/            tier 8: result plots
├── tests/
│   ├── conftest.py  annotation_fixtures.py  ref_panel_helpers.py  worker_fixtures.py   shared helpers
│   ├── fixtures/             small checked-in inputs and golden expectations
│   ├── core/  runtime/  genome/  annotation/  refpanel/  sumstats/  ldscore/  regression/  plotting/
│   │                         one folder per domain package, named without the underscore
│   └── package/              root shell: package layout, import tiers, CLI help, configuration
└── tools/
```

Every package whose name starts with an underscore is private; the public Python surface is the `ldsc` namespace (`src/ldsc/__init__.py`), and the `ldsc` command is `src/ldsc/cli.py`.
A domain package's `__init__.py` holds only a docstring: it imports nothing, so importing one of its modules never loads the whole domain, and the SciPy-backed workflows that `ldsc` exports lazily stay unloaded until used.
Imports between domains name the module, for example `from .._ldscore.source import LDScoreSource`.

## Dependency Direction

Packages are arranged in tiers, and a module may import from its own package and from lower tiers only:

| Tier | Package or module |
| --- | --- |
| 0 | `errors.py` |
| 1 | `_core` |
| 2 | `_runtime` |
| 3 | `config.py` |
| 4 | `_genome` |
| 5 | `_annotation`, `_refpanel`, `_sumstats`; these peers do not import each other |
| 6 | `_ldscore` |
| 7 | `_regression` |
| 8 | `_plotting` |
| 9 | `cli.py`, `__init__.py`, `__main__.py` |

[`tests/package/test_layering.py`](../../tests/package/test_layering.py) enforces the table: it reads every import inside the package, including imports inside functions and under `TYPE_CHECKING`, and fails on an import that points upward or between peers.
It also fails when a module outside the root shell sits at the package root.
Function-level imports remain for optional or heavy dependencies, such as Matplotlib and the SciPy-backed workflows, but never to avoid an import cycle.

- `ldsc.cli` parses arguments and dispatches; workflow logic lives in the domain packages.
- Reference panels hand LD-score computation a `PreparedChromosome` (`_refpanel.prepared`); `_ldscore.compute` consumes it and never reopens sources.
- Each writer lives with the domain that owns its artifact; `_runtime.artifacts` holds the declaration and atomic-write base they share.
- Regression reloads written LD-score directories; it does not depend on annotation or reference-panel code directly.
- `_plotting.dispatch` and `_regression.h2_scale` consume canonical result directories; only the private plotting builders and the prevalence-range conversion path import Matplotlib, lazily.

## Module Map

| Module | Responsibility |
| --- | --- |
| `ldsc.cli` | unified `ldsc` command and subcommand dispatch |
| `ldsc.config` | frozen public config dataclasses and basic validation |
| `ldsc._core.chromosomes` | canonical chromosome normalization and ordering |
| `ldsc._core.columns` | resolve header aliases and normalize identifier/build tokens |
| `ldsc._core.coordinates` | `CHR:POS` keys and coordinate-frame normalization |
| `ldsc._core.snp_identity`, `ldsc._core.identifiers` | the four SNP identifier modes, base/effective identity keys, allele-aware cleanup, restriction files, `identifier_keys` and retained-identity uniqueness checks |
| `ldsc._core.regions` | packaged and user BED interval loading plus region-exclusion masks |
| `ldsc._core.row_alignment` | row-alignment checks between SNP tables and the shared genomic sort |
| `ldsc._core.plink_files` | retained PLINK `.bim`/`.fam` and identifier-list readers; obsolete legacy regression-artifact readers are removed |
| `ldsc._runtime.paths` | normalize path tokens, resolve concrete input files, create output directories, preflight fixed output paths, and enforce coherent output artifact families |
| `ldsc._runtime.preflight` | the declared-input gate (`InputGate`, `inspect_declared_inputs`), issue records, companion and saved-manifest checks; domains plug in their own declaration and header checks, and scientific content validators stay with their owners |
| `ldsc._runtime.artifacts` | `ArtifactFamily` declarations shared by workflow preflight, final writes, metadata file entries and cleanup, plus atomic table and JSON writes |
| `ldsc._runtime.result_files` | contained declared-result input paths and atomic JSON publication for `h2_scale` and plotting |
| `ldsc._runtime.progress` | phase boundaries and time-throttled progress through existing logger routes; bounded counters, current objects, and elapsed time |
| `ldsc._runtime.workflow_logging` | shared workflow logging context, `LDSC` logger level handling, lifecycle audit lines, CLI console-handler routing (file-authoritative, console error-only) with run-aborting traceback capture, durable authorized-overwrite failure markers without rollback, and log-only formatting helpers |
| `ldsc._runtime.cli_help` | shared `CLIHelpFormatter` keeps flag names and paths intact when help descriptions wrap; `CHROMOSOME_PATH_HELP` and `SCALAR_PATH_HELP` supply verified path syntax; workflow parsers own their descriptions and argument groups |
| `ldsc._runtime.workflow_kwargs` | validate Python keyword arguments for `run_ldscore`, `run_build_r2_matrix`, `run_h2`, `run_partitioned_h2`, `run_rg`, and `run_quantile_h2` against the command parser before any write; run the regression and post-fit workflows under a failure marker naming the Python function |
| `ldsc._runtime.parallelism` | shared `workers` validation, joblib detected CPU budgets, task caps, execution summaries, chromosome task scheduling, and the spawn worker pool that turns a crashed worker into `LDSCInternalError` |
| `ldsc._runtime.worker_runtime`, `ldsc._runtime.serial_calculation` | fixed one-BLAS-thread worker policy, runtime verification, independent Arrow limits, and controlled standalone numerical calls |
| `ldsc._genome.build_inference` | public `chr_pos` build and coordinate-basis inference helpers |
| `ldsc._genome.build_policy`, `ldsc._genome.build_inspection` | the build policy workflows apply to their inputs and read-only build inspection of saved artifacts (`infer-genome-build`) |
| `ldsc._genome.liftover` | shared hg19/hg38 liftover helpers, chain-file translation, dual-build HM3 quick liftover, drop-all coordinate collision helpers, and readable drop reports |
| `ldsc._genome.hm3` | public packaged curated HM3 map loader and installed map path helper for workflow internals |
| `tools/hm3/build_hm3_chr_pos_reference.py` | maintenance tool (outside the package) that rebuilds the compact HM3 coordinate reference used by genome-build inference |
| `ldsc._annotation.builder` | public annotation builder and general standalone parser/API dispatch seam |
| `ldsc._annotation.annotate`, `ldsc._annotation.interval_queries` | standalone BED/gene gates and incremental projection |
| `ldsc._annotation.interval_projection` | `ChromosomeProjector` maps BED/gene interval unions to chromosome SNP positions; called by `interval_queries.build_query_shards()` after shared BED parsing |
| `ldsc._annotation.bundle`, `ldsc._annotation.storage` | explicit dataset ownership, chromosome descriptors, detached selected reads, and output-contained scratch |
| `ldsc._annotation.sources`, `ldsc._annotation.identity_cleanup`, `ldsc._annotation.loading` | bounded input scans, aligned source preparation, global identity cleanup, and shared source-backed construction |
| `ldsc._annotation.parsing`, `ldsc._annotation.tables` | annotation chunk normalization and concrete annotation table reading with numeric-value validation |
| `ldsc._annotation.preflight` | annotation declaration resolution and `inspect_annotation_inputs`, the declared-input gate with annotation header checks first |
| `ldsc._annotation.saved_results` | saved baseline/query descriptors with no original-input dependencies; deferred preparation on explicit reads and construction-owner detachment |
| `ldsc._annotation.bim_grid` | metadata-only standalone BIM grid scans using shared identity cleanup; no genotype companions or value synthesis |
| `ldsc._annotation.gene_lists` | required/embedded gene-coordinate-catalog validation and vectorized query/control resolution into exact IDs, intervals, audits, and summaries |
| `ldsc._annotation.gene_query_storage` | staged per-source gene selections and complete audit replay with shared gene-coordinate-catalog state |
| `ldsc._annotation.query_outcomes` | shared gene coverage, resolution/support precedence, required-control and all-skipped failures, ordered diagnostics, and aligned result pruning; does not import workflow orchestration |
| `ldsc._annotation.in_memory` | small prepared annotations and diagnostics without filesystem writes |
| `ldsc._annotation.semantics` | global annotation-name uniqueness and advisory binary/continuous classification |
| `ldsc._annotation.overlap` | annotation overlap-block computation (`OverlapContribution`, `compute_overlap`, `sum_overlap_contributions`); a product of annotations over reference SNPs without LD |
| `ldsc._annotation.outputs` | `AnnotationDirectoryWriter` for the `ldsc annotate` result directory |
| `ldsc._refpanel.plink` | PLINK genotype reader (`PlinkBEDFile` and its `__GenotypeArrayInMemory__` base): selective SNP reads with fused individual filtering retain only selected packed genotypes; unrestricted panel builds opt into disk streaming instead |
| `ldsc._refpanel.r2_format` | the R² matrix format: file names (`r2_matrix_file_paths`, `r2_matrix_chromosomes`, `resolve_r2_build_dir`), schema constants, Parquet and sidecar writers, `R2MatrixChromosome` (footer contract, strict sidecar reading, sidecar binding, LD-window and truncation checks for `ldscore` and `query-r2-matrix`), `SortedR2BlockReader` and the pair lookup used by `query-r2-matrix` |
| `ldsc._refpanel.windows` | LD-window geometry shared by `ldscore` and `build-r2-matrix`: `window_coordinates`, `get_block_lefts`, and window checks |
| `ldsc._refpanel.prepared` | the `PreparedChromosome` contract, the chromosome annotation input `ChromosomeAnnotations`, and the `ArrayAnnotations`/`MappedAnnotations` accessors |
| `ldsc._refpanel.funnel` | the reference-SNP funnel: `ReferenceSNPFunnel` records what each reference-preparation step removed, with readers, writers and the run-level warning text |
| `ldsc._refpanel.adapters` | runtime PLINK/R²-matrix adapters; `RefPanel.prepare_chromosome` owns reference filtering, annotation alignment, window/bias policy, and reader lifetime |
| `ldsc._refpanel.builder`, `ldsc._refpanel.build_kernel` | R²-matrix build workflow (source-build inference, optional liftover, explicit SNP/sample restrictions, optional `r2_min` threshold) and its kernel (genetic maps, pairwise R² emission, restriction masks) |
| `ldsc._refpanel.query` | public `query-r2-matrix` CLI/API, `R2Matrix`, one-shot `query_r2_matrix()`, endpoint key resolution, sign harmonization, and optional unbiased-R²-to-Pearson-r conversion |
| `ldsc._refpanel.preflight` | `inspect_r2_headers` (each R² matrix chromosome's files, footer and sidecar header, for either reading command) and the genetic-map header check for the declared-input gate |
| `ldsc._refpanel.outputs` | `QueryR2MatrixDirectoryWriter` for the `query-r2-matrix` result directory |
| `ldsc._sumstats.munger` | raw-summary-statistics CLI/API orchestration, `--input-format auto` / `--infer-only` header inference, Parquet/TSV munged output writing, self-describing `sumstats.parquet` footer identity metadata, diagnostics under `diagnostics/`, canonical `CHR`/`POS` summary-statistics output, and munged summary-statistics loader |
| `ldsc._sumstats.munge_input` | private workflow helper resolving raw schema, DANER/sample-size settings, bounded source-build evidence and keep-lists into `ResolvedMungeInput` |
| `ldsc._sumstats.munge_kernel` | chunk QC and restriction, whole-table N and sign conversion, resolved liftover and global identity cleanup; returns `MungeResult` with counts and provenance |
| `ldsc._ldscore.compute` | LD-score math over a `PreparedChromosome`: `compute_chromosome`, R² pair streaming, PLINK block sums (`ld_scores_from_genotypes`), the float64 `ProjectionAccumulator`, counts and the regression-SNP mask; no legacy artifact emitters |
| `ldsc._ldscore.calculator` | LD-score orchestration, gene-coordinate-catalog build selection, direct prepared-SNP support measurement, optional synthetic `base`, aggregation, and output routing through shared query finalization |
| `ldsc._ldscore.preflight` | validates direct annotation/reference contents and chromosome equality; enforces complete `@` declarations, authoritative glob matches, and consolidated required-input diagnostics |
| `ldsc._ldscore.direct_annotation`, `ldsc._ldscore.direct_batches` | direct BED/gene query integration into an LD-score run and staged direct query scores |
| `ldsc._ldscore.batch_execution` | the query-batch loop shared by the direct and indexed LD-score producers: one batch in flight through the run's persistent pool, failure shutdown before scratch removal, and query-status finalization |
| `ldsc._ldscore.fragments`, `ldsc._ldscore.batch_publication` | producer-owned chromosome fragments written one query frame at a time against the shared baseline, and publication of complete LD-score directories |
| `ldsc._ldscore.diagnostics`, `ldsc._ldscore.query_variation` | chromosome logs, batch context and task/pool failure reporting, and the write-log-raise sequences both LD-score modes use for rejected gene lists and unusable queries; exact query-column variation across chromosomes |
| `ldsc._ldscore.source` | canonical directory loading, manifest/schema and table/identity validation shared with writers, resident baseline/count/overlap metadata, and uncached selected-query reads |
| `ldsc._ldscore.overlap_table` | `LDScoreOverlap` and its long-frame Parquet read/write for the LD-score directory's overlap table |
| `ldsc._ldscore.outputs` | `LDScoreDirectoryWriter`: LD-score Parquet layout with chromosome row groups, query batch files, query diagnostics and the root `metadata.json` |
| `ldsc._ldscore.legacy_converter` | sole legacy LDSC LD-score-suite import boundary: deterministic family discovery, rsID joins, count/overlap validation or reconstruction, source/discovery provenance without source-file hashes, diagnostics, and canonical LDSC3 directory writing |
| `ldsc._ldscore.gene_index.build` | gene LD-score index building and transactional publication |
| `ldsc._ldscore.gene_index.artifact` | the on-disk index format: chromosome payloads, `index_id` and published-row hashes, staged writes, readers and component inspection |
| `ldsc._ldscore.gene_index.scoring`, `ldsc._ldscore.gene_index.batches` | `load_gene_ldscore_index`, `run_indexed_ldscore`, and bounded chromosome workers holding one operator each |
| `ldsc._ldscore.gene_index.atoms` | sparse disjoint-atom primitives for exact gene LD-score indexes |
| `ldsc._regression.runner` | file-driven regression dataset assembly, automatic legacy summary-statistics rsID-to-panel projection and allele harmonization, active effective identity-key merging (`SNP`, `SNP:<allele_set>`, `CHR:POS`, or `CHR:POS:<allele_set>`), `h2`/`partitioned-h2`/`rg` estimator dispatch (including the two overlap-aware `partitioned-h2` regimes), observed/liability-scale summary columns, and exact final-fit `h2` regression-bin diagnostics |
| `ldsc._regression.partitioned_parallel` | private read-only regression maps, bounded spawned whole-query workers, and worker termination/log forwarding; statistical calculations remain in `RegressionRunner._fit_partitioned_query()` |
| `ldsc._regression.estimators` | LDSC estimators for `Hsq` and `RG` |
| `ldsc._regression.jackknife`, `ldsc._regression.irwls` | supporting numerical routines used by regression |
| `ldsc._regression.overlap_model` | per-model overlap assembly, the overlap-aware annotation table (ported `_overlap_output` with augmentation) and the collinearity check (`model_collinearity_error`) |
| `ldsc._regression.prevalence` | parse and validate binary-trait prevalence inputs (scalar `--samp-prev`/`--pop-prev` for `h2`/`partitioned-h2`; comma-separated lists or a `--prevalence-manifest` TSV for `rg`) into a normalized per-trait `(samp_prev, pop_prev)` structure for observed-to-liability conversion |
| `ldsc._regression.h2_scale` | strict post-fit observed-to-liability conversion from a canonical `h2` result, including exact and prevalence-range modes and the fixed nested derived-result family |
| `ldsc._regression.quantile_h2`, `ldsc._regression.quantile_inputs`, `ldsc._regression.quantile_storage` | post-fit continuous-target quantile assignment, reconstruction and verification of the fitted sources and common reference SNPs, vectorized coefficient/delete-value projection, standardized `tau_star`, and streamed alignment diagnostics |
| `ldsc._regression.outputs` | writers for the `h2`, `partitioned-h2` (per-query layout), `rg` (result family) and `quantile-h2` result directories, with their column constants |
| `ldsc._plotting.dispatch` | sole public metadata-driven plotting dispatcher, lazy plotting-runtime boundary, fixed plot-family output, and live `PlotArtifact` return object |
| `ldsc._plotting.builders` | private Matplotlib-only headless builders for the approved `h2`, `rg`, `partitioned-h2`, and `quantile-h2` plots |

## Where To Change Code

| Goal | Start here |
| --- | --- |
| change CLI flags or subcommand wiring | workflow `build_parser` / `add_*_arguments` functions, then `src/ldsc/cli.py` for dispatch (`main`, `_command_parser`); follow [CLI help guidelines](guidelines/cli-help.md) |
| change help wrapping or shared log-level wording | `src/ldsc/_runtime/cli_help.py:CLIHelpFormatter`, `src/ldsc/_runtime/workflow_logging.py:LOG_LEVEL_HELP`; verify `tests/package/test_cli_help.py` |
| change path-token behavior | `src/ldsc/_runtime/paths.py` |
| change output collision policy | `src/ldsc/_runtime/paths.py`, then the workflow writer that owns the artifact |
| change a keyword-argument Python workflow (`run_ldscore`, `run_build_r2_matrix`, `run_h2`, `run_partitioned_h2`, `run_rg`, `run_quantile_h2`) | the function in its workflow module, keyword checks in `src/ldsc/_runtime/workflow_kwargs.py`; verify `tests/package/test_workflow_keyword_checks.py` and `tests/regression/test_python_regression_workflows.py` |
| change workflow log files, lifecycle lines, or log-level handling | `src/ldsc/_runtime/workflow_logging.py`, then the workflow orchestration function that owns the output directory |
| change header aliases or identifier/build normalization | `src/ldsc/_core/columns.py` |
| change read-only artifact build inspection | `src/ldsc/_genome/build_inspection.py`, with CLI dispatch in `src/ldsc/cli.py` and coverage in `tests/genome/test_genome_build_inspection.py` |
| change automatic `chr_pos` genome-build inference | `src/ldsc/_genome/build_inference.py` |
| change annotation loading or `ldsc annotate` behavior | `src/ldsc/_annotation/builder.py`, then `src/ldsc/_annotation/sources.py` and `src/ldsc/_annotation/annotate.py` |
| change BED parsing or interval projection | `src/ldsc/_core/regions.py:iter_bed_rows`, then `src/ldsc/_annotation/interval_queries.py` and `src/ldsc/_annotation/interval_projection.py:ChromosomeProjector` |
| change gene coordinate catalog validation, identifier resolution, or list naming | `src/ldsc/_annotation/gene_lists.py`, then `src/ldsc/_annotation/builder.py` |
| change BED/gene partial-success statuses or diagnostic schemas | `src/ldsc/_annotation/query_outcomes.py`, `src/ldsc/_ldscore/calculator.py`, `src/ldsc/_ldscore/outputs.py` |
| change R²-matrix build logic | `src/ldsc/_refpanel/builder.py`, then `src/ldsc/_refpanel/build_kernel.py` |
| change named regression-SNP region exclusion for `ldscore` or gene LD-score index construction | `src/ldsc/config.py`, `src/ldsc/_core/regions.py`, then `src/ldsc/_ldscore/calculator.py` or `src/ldsc/_ldscore/gene_index/build.py` |
| change R² pair lookup CLI/API behavior | `src/ldsc/_refpanel/query.py`, then `src/ldsc/_refpanel/r2_format.py:lookup_pairs_in_parquet` |
| change runtime PLINK/R²-matrix reference access | `src/ldsc/_refpanel/adapters.py`, with the R² file format in `src/ldsc/_refpanel/r2_format.py` and LD windows in `src/ldsc/_refpanel/windows.py` |
| change LD-score orchestration, optional-baseline behavior, or output packaging | `src/ldsc/_ldscore/calculator.py`, `src/ldsc/_ldscore/outputs.py` |
| change saved LD-score loading or selected-query reads | `src/ldsc/_ldscore/source.py`; keep model-specific validation in `src/ldsc/_regression/runner.py` |
| change LD-score math | `src/ldsc/_ldscore/compute.py` |
| change raw summary-statistics ingestion, format inference, `CHR`/`POS` handling, summary-statistics SNP keep-list filtering, liftover drop audit sidecars, footer provenance, or munged loading | `src/ldsc/_sumstats/munger.py`, then `src/ldsc/_sumstats/munge_kernel.py` |
| change regression dataset assembly or CLI summaries | `src/ldsc/_regression/runner.py`, then `src/ldsc/_regression/outputs.py`, `docs/reference/results/partitioned-h2.md` for `partitioned-h2` output layout, `docs/reference/formats/artifact-metadata.md` for `rg` output contracts, and `docs/reference/regression-options.md` for the tunable estimator parameters and defaults |
| change annotation-name validation or continuous-annotation classification | `src/ldsc/_annotation/semantics.py`, then `src/ldsc/_ldscore/calculator.py` and `src/ldsc/_ldscore/outputs.py` |
| change target quantiles, standardized coefficients, or post-fit alignment | `src/ldsc/_regression/quantile_h2.py`, then `src/ldsc/_regression/outputs.py` and `docs/reference/results/quantile-h2.md` |
| change LDSC estimators | `src/ldsc/_regression/estimators.py` |
| change binary-trait prevalence parsing or observed-to-liability conversion | `src/ldsc/_regression/prevalence.py` (input parsing/validation), then `src/ldsc/_regression/runner.py` (summary/metadata wiring) and `src/ldsc/_regression/estimators.py` (`liability_conversion_factor`) |
| change post-fit heritability scale conversion or sensitivity output | `src/ldsc/_regression/h2_scale.py`, then `src/ldsc/_regression/estimators.py` only if the shared numerical primitive itself changes |
| change plot dispatch, supported result contracts, or visual semantics | `src/ldsc/_plotting/dispatch.py`, `src/ldsc/_plotting/builders.py`, then `docs/dev/design/plotting.md` and `docs/user-guide/workflows/plotting.md` |
| change owned artifact declarations | `src/ldsc/_runtime/artifacts.py`, then the domain writer: `src/ldsc/_annotation/outputs.py`, `src/ldsc/_refpanel/outputs.py`, `src/ldsc/_ldscore/outputs.py` or `src/ldsc/_regression/outputs.py` |
| change LD-score directory files, Parquet row-group layout, or LD-score metadata | `src/ldsc/_ldscore/outputs.py` |
| change the `partitioned-h2` per-query layout, the `rg` result-family layout, or regression metadata payloads | `src/ldsc/_regression/outputs.py` |
| change which module may import which | the tier table above and `tests/package/test_layering.py` |

## Architectural Rules That Matter In Practice

- Treat `src/ldsc/` as the only supported Python import surface.
- Do not add user-facing path discovery to computation modules such as `_ldscore.compute`, `_refpanel.r2_format` or `_regression.estimators`; pass concrete files in.
- Keep public file contracts for `.annot(.gz)`, self-describing Parquet munged summary statistics (footer identity metadata) plus optional `.sumstats.gz` compatibility output, LD-score directories, and regression summary directories stable unless the change is intentional and coordinated.
  LD-score parquet files remain flat files, with chromosome-aligned row groups documented through root `metadata.json`; diagnostic logs, dropped-SNP reports, and diagnostic metadata live under `diagnostics/`.
  Legacy reference/weight `.l2.ldscore(.gz)`, `.l2.M`, and `.l2.M_5_50` families enter through the explicit [LD-score converter](../reference/formats/legacy-ldscores.md#compatibility-boundary).
  The `.w.l2.ldscore(.gz)` suffix is unsupported, and private kernel legacy artifact emitters and regression readers have been removed.
- Treat package-built R² matrices as the only public R² Parquet
  format. Both `ldsc ldscore --r2-matrix-dir` and `ldsc query-r2-matrix` depend on the paired
  `chr*_r2.parquet` / `chr*_meta.tsv.gz` layout and sidecar identity binding.
- Keep optional-baseline behavior in the public LD-score workflow layer: no baseline and no query means a synthetic all-ones `base`; prebuilt, BED, and gene-list queries require explicit baseline annotations. Skipped BED/gene queries never become placeholder scientific columns.
- Keep `ldsc annotate` orchestration in `ldsc._annotation.builder`. The CLI
  registers annotation flags from that module and dispatches parsed namespaces
  to `run_annotate_from_args()`; `_annotation.tables` must not own parser
  functions, result objects, or public workflow aliases.
- Use the shared output preflight helpers for fixed output files, including
  workflow log files. Public workflows should create missing output
  directories, reuse existing directories, fail on existing owned artifacts by
  default, and require `--overwrite` or `overwrite=True` for replacement.
- Keep coherent output families consistent. For `munge-sumstats`,
  `build-r2-matrix`, `ldscore`, `partitioned-h2`, `rg`, and `annotate`,
  no-overwrite mode rejects any current-contract owned sibling from the workflow family.
  Overwrite mode writes the requested current outputs and then removes stale
  current-contract owned siblings that were not produced by the successful run.
  Removed legacy root diagnostic names are not workflow-owned. Unrelated files
  in the output directory must be preserved. Sharded workflows may pass a
  narrowed owned family when shards share an output directory; `build-r2-matrix`
  concrete chromosome prefixes own only that chromosome's package, while `@`
  chromosome-suite prefixes own the full R² matrix package.
- Preflight deterministic output paths before expensive or multi-file writes.
  This is especially important for `build-r2-matrix`, `munge-sumstats`,
  `annotate`, and summary-table regression commands.
- Treat workflow logs as diagnostic files. Do not include log paths in
  `output_paths` mappings or thin metadata sidecars that downstream code
  interprets as scientific data artifacts. Whole-directory transactional
  publishers must keep open workflow logs outside the replaceable tree; the
  gene LD-score index builder uses hidden sibling `.<index-name>.build-state/` and moves
  only the closed successful log into the published diagnostics.
- Keep regression file-driven: it should be able to rebuild state from written artifacts without recomputing LD scores.
- Keep plotting metadata-driven and post-fit: plot builders consume declared canonical tables and never refit. Matplotlib is required at installation but lazily imported through `ldsc._plotting.builders`.
- Prefer extending shared helpers or the workflow-owned writer over duplicating local parsing or writing logic.

## Test Map

Each test file sits in the folder of the domain package whose code it mainly exercises, such as `tests/ldscore/` for `_ldscore`.
Tests of shared machinery that every workflow uses, such as output writers, workflow logs, failure markers, workers and input preflight, are in `tests/runtime/`.
Tests of the root shell and of the package as a whole are in `tests/package/`.
Shared helpers and `fixtures/` stay at the top of `tests/`; tests import them, and one another, by full module path, for example `from tests.annotation.test_annotate_streaming import inputs`.
Each test folder has an `__init__.py`, which makes those full module paths importable.

| Area | Main tests |
| --- | --- |
| output layer | `tests/runtime/test_output.py` |
| annotation workflow | `tests/annotation/test_annotation.py` |
| reference-panel builder | `tests/refpanel/test_r2_matrix_builder.py` |
| R² pair query | `tests/refpanel/test_r2_matrix_query.py` |
| LD-score workflow | `tests/ldscore/test_ldscore_workflow.py` |
| gene coordinate catalog and gene-list resolution | `tests/annotation/test_gene_list_resolver.py` |
| summary-statistics munging | `tests/sumstats/test_sumstats_munger.py` |
| regression workflow | `tests/regression/test_regression_workflow.py` |
| path and config contracts | `tests/runtime/test_path_resolution.py`, `tests/core/test_config_identifiers.py`, `tests/core/test_column_inference.py`, `tests/genome/test_genome_build_inference.py` |
| package layout and import tiers | `tests/package/test_layering.py`, `tests/package/test_package_layout.py` |
| region exclusion | `tests/core/test_regions.py` |

`AnnotationBundle.validate()` validates row alignment and annotation names; identity cleanup occurs during builder construction. It no longer accepts an ignored identifier argument. Column-name inference uses `infer_chr_pos_columns`; numerical window construction uses `get_block_lefts`. Their former forwarding spellings are removed without changing accepted input column aliases or legacy formats.

`ldsc._runtime.result_files` owns contained declared-result input paths and atomic JSON publication for `h2_scale` and plotting.
Each workflow retains its metadata interpretation, table loading and scientific dispatch.
The result writers keep their separate JSON serialization of NumPy/dataclass payloads in `ldsc._runtime.artifacts`.

## Domain contracts and decision history

| Change | Current contract | Main owner |
| --- | --- | --- |
| annotations and interval projection | [BED input](../reference/formats/bed-intervals.md), [gene lists](../reference/formats/gene-lists.md), [query repair](../reference/results/gene-list-diagnostics.md) | `_annotation.builder.AnnotationBuilder`; `_annotation.interval_queries.build_query_shards`; `_annotation.interval_projection.ChromosomeProjector` |
| panel preparation, CM/MAF and traversal | [R² pipeline](../reference/formats/r2-matrices.md), [pair accumulation](../reference/methods/ld-score-computation.md), [window behavior](../reference/formats/r2-matrices.md#37-ld-windows-and-sidecar-values-in-ldscore) | `_refpanel.adapters.RefPanel.prepare_chromosome`; `_ldscore.compute.compute_chromosome` |
| regression SNPs and region masks | [reference and regression SNPs](../reference/ldscore-snp-universes.md), [region exclusions](../reference/formats/region-exclusion-presets.md) | `_ldscore.calculator._regression_region_intervals`; `_ldscore.compute.regression_mask_from_keys`; `_ldscore.gene_index.build` resolves indexed policy |
| reference construction and pair queries | [pair query](../reference/query-r2-matrix.md), [path contract](../reference/input-paths.md) | `_refpanel.builder.R2MatrixBuilder`; `_refpanel.query.R2Matrix`; `_refpanel.outputs.QueryR2MatrixDirectoryWriter` |
| munging and liftover | [munging](../reference/munge-sumstats.md), [liftover](../reference/snp-identity-and-genome-build.md#liftover) | `_sumstats.munger.SumstatsMunger`; `_sumstats.munge_input.prepare_munge_input`; `_sumstats.munge_kernel` |
| regression and overlap | [configuration](../reference/regression-options.md), [partitioned results](../reference/results/partitioned-h2.md) | `_regression.runner.RegressionRunner`; `_ldscore.overlap_table`; `_regression.overlap_model`; `_regression.estimators` |
| indexed gene queries | [index](../reference/gene-ldscore-index.md), [mathematics](../reference/methods/gene-ldscore-index-math.md) | `_ldscore.gene_index.build`, `.artifact`, `.scoring`; `_ldscore.gene_index.atoms` |
| derived results | [plotting](design/plotting.md), [quantile projection](../reference/results/quantile-h2.md) | `_plotting.dispatch`, `_regression.h2_scale`, `_regression.quantile_h2`; shared `_runtime.result_files` |

`query-r2-matrix` requires `--output-dir`; `QueryR2MatrixDirectoryWriter` writes its result directory.
Region presets affect the regression SNPs and their weights, not the reference SNPs.
`build-r2-matrix` has no region-pruning fields; use its explicit reference-SNP restriction when pruning is intended.
Reference preparation supplies population-specific CM/MAF; annotation CM is an `NA` placeholder.
LD scores from R² matrices use pair streaming, and PLINK filtering happens during adapter preparation.

Design rationale that protects a non-obvious constraint is stated in the page that documents the constraint; measurements are in [validation](validation/performance.md).
Earlier dated specifications, implementation plans, audits, and the archived design map were removed on 2026-10-26, and the decision records on 2026-10-27; they remain in git history, for example `git show 479fc4d:docs/specs/<file>` or `git show 77a54e2:docs/decisions/2026-10-11-parallel-execution.md`.
The specification and plan for the domain-package layout were removed once it was implemented on 2026-10-28; see `git show 869f41e:docs/specs/2026-10-27-domain-subpackages.md`.
