# Code Structure

Last updated on: 2026-11-06

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
│   ├── _runtime/             tier 2: paths, input preflight, run scratch, logging, progress, workers, keyword checks, result-file writers
│   ├── _genome/              tier 4: genome-build inference and policy, build inspection, liftover, HM3 map
│   ├── _annotation/          tier 5: annotation files, gene lists, interval projection, overlap, ldsc annotate
│   ├── _refpanel/            tier 5: PLINK reader, R² matrix format, LD windows, prepared-chromosome contract, reference-SNP funnel, runtime adapters, build-r2-matrix, query-r2-matrix
│   ├── _sumstats/            tier 5: munge-sumstats and the munged-table loader
│   ├── _ldscore/             tier 6: LD-score math and workflow, LD-score directories, legacy conversion
│   │   └── gene_index/       exact gene LD-score indexes
│   ├── _regression/          tier 7: estimators, h2/partitioned-h2/rg, liability scale, quantile-h2
│   └── _plotting/            tier 8: result plots
├── tests/
│   ├── conftest.py  annotation_fixtures.py  gene_index_helpers.py  ldscore_helpers.py  quantile_fixtures.py
│   │   ref_panel_helpers.py  regression_fixtures.py  worker_fixtures.py    shared helpers
│   ├── fixtures/             small checked-in inputs and golden expectations
│   ├── core/  runtime/  genome/  annotation/  refpanel/  sumstats/  ldscore/  regression/  plotting/
│   │                         one folder per domain package, named without the underscore
│   └── package/              root shell: package layout, import tiers, CLI help, configuration
└── tools/                    data-curation and figure scripts outside the package
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
In-package imports sit at module level, so a module's header lists its dependencies.
No modules may import one another, directly or through other modules; the same test fails on any such cycle in the imports Python runs, counting imports inside functions and leaving out imports under `TYPE_CHECKING`.
An in-package import stays inside a function only when moving it to module level would make `import ldsc` or `ldsc --help` load another module, and never to avoid an import cycle.
The Matplotlib plot builders follow their own rule: Matplotlib loads only on paths that draw a figure, so their import stays inside `_plotting.dispatch`, which `import ldsc` and `ldsc --help` do not load but `ldsc plot --help` does.
`DEFERRED_IMPORTS` in `tests/package/test_layering.py` lists each such import with its reason, and the test fails on an unlisted function-level import and on a stale entry.
It also fails when `import ldsc` or `ldsc --help` loads a listed module, when they do not load the importer of a listed import other than the plot builders (so a module-level import would not slow startup), and when `ldsc plot --help` loads the plot builders.
The lazily exported workflows in `ldsc/__init__.py` and the per-command modules of `ldsc.cli` are loaded by name through `importlib`, and optional third-party dependencies such as pyliftover keep their own guarded imports; neither is an in-package import statement.

- `ldsc.cli` parses arguments and dispatches; workflow logic lives in the domain packages.
- Reference panels hand LD-score computation a `PreparedChromosome` (`_refpanel.prepared`); `_ldscore.compute` consumes it and never reopens sources.
- Each writer lives with the domain that owns its artifact; `_runtime.artifacts` holds the declaration and atomic-write base they share.
- Regression reloads written LD-score directories; the fit commands (`h2`, `partitioned-h2`, `rg`) import no annotation or reference-panel module directly.
  They reach annotation code only through the LD-score loader: `_ldscore.source` checks annotation names with `_annotation.semantics`, and `_ldscore.overlap_table` labels `_annotation.overlap.OverlapContribution` values.
  `quantile-h2` is the only regression command that imports annotation code directly: it rebuilds the fitted model's annotations, BED and gene-list query annotations included, with `_annotation.loading.build_annotation_shards()`.
  When that preparation stops at the gene-list Gate A check, it writes the gene-list audit and resolution summary with `_annotation.outputs.AnnotationDirectoryWriter.write_diagnostics()`.
  When it stops on a `_annotation.gene_lists.GeneCatalogValidationError`, it writes the gene coordinate catalog issue table with its own `_regression.outputs.QuantileH2DirectoryWriter.write_gene_catalog_issues()`.
- Run scratch (`ScratchWorkspace`), frame spools and the diagnostic-row handles live in `_runtime.scratch`, so annotation, LD-score, gene-index and regression runs share them without depending on each other.
- `_plotting.dispatch` and `_regression.h2_scale` consume canonical result directories; only the private plotting builders and the prevalence-range conversion path import Matplotlib, lazily.

## Module Map

| Module | Responsibility |
| --- | --- |
| `ldsc.cli` | unified `ldsc` command; the `_COMMANDS` table maps each subcommand to its workflow module and entry functions, including the function that prints a `run_*_from_args` command's successful-run notices (`notices`) |
| `ldsc.config` | frozen public config dataclasses and basic validation, including the one binary-trait prevalence rule (`validate_prevalences()`, vectorized, which `normalize_prevalence_pair()`, `liability_conversion_factor()` and `convert-h2-scale` call) and the one LD-window rule that `LDScoreConfig`, `R2MatrixBuildConfig` and `GeneLDScoreIndexBuildConfig` share (`LDWindow`, returned by their `ld_window` property) |
| `ldsc._core.chromosomes` | canonical chromosome normalization and ordering |
| `ldsc._core.columns` | resolve header aliases and normalize identifier/build tokens |
| `ldsc._core.coordinates` | `CHR:POS` keys and coordinate-frame normalization |
| `ldsc._core.snp_identity`, `ldsc._core.identifiers` | the four SNP identifier modes (`SNP_IDENTIFIER_MODES`) and the default mode (`DEFAULT_SNP_IDENTIFIER`), which the parsers and `GlobalConfig` use, the one allele cleaner (`canonical_allele_series`, which munging and the HM3 map loader also use), `identity_key_series`, which builds the SNP identity of every row of a canonical table, and its base-identity step `base_key_series` (called directly by identity cleanup, restriction matching and R² matrix queries, which handle invalid alleles themselves), allele-aware cleanup, the one base-complement table and the allele-matching tables (`MATCH_ALLELES`, `FLIP_ALLELES`) that regression alignment reads, restriction files, retained-identity uniqueness checks, and `open_identity_table`, which opens package-written text tables past their leading `# ldsc:` identity lines for the R² matrix sidecar reader, the genome-build check of a sidecar's declaration and `quantile-h2` reference metadata |
| `ldsc._core.regions` | packaged and user BED interval loading plus region-exclusion masks |
| `ldsc._core.row_alignment` | row-alignment checks between SNP tables and the shared genomic sort |
| `ldsc._core.plink_files` | retained PLINK `.bim`/`.fam` and identifier-list readers; obsolete legacy regression-artifact readers are removed |
| `ldsc._runtime.paths` | normalize path tokens, expand `@` chromosome suites (`expand_chromosome_suite`), resolve concrete input files, discover PLINK trios (`discover_plink_members`, `resolve_plink_prefix_group`), create output directories, preflight fixed output paths, and enforce coherent output artifact families |
| `ldsc._runtime.preflight` | the declared-input gate (`InputGate`, `inspect_declared_inputs` and its result `DeclaredInputs`: resolved files and check results by role, plus discovered PLINK members), companion and saved-manifest checks; domains plug in their own declaration and header checks, and scientific content validators stay with their owners |
| `ldsc._runtime.input_issues` | the input-issue row schema (`INPUT_ISSUE_COLUMNS`) and its one constructor `input_issue`, whose repair text every caller supplies; used by every input gate, including PLINK discovery in `ldsc._runtime.paths` |
| `ldsc._runtime.scratch` | `ScratchWorkspace`, the owner of one hidden `.ldsc-annotation-*` scratch directory below an output directory for a run or a query batch; `FrameSpool`; and the `DiagnosticRows` protocol (`frames()`, `write_to(path)`) shared by spooled, file-staged (`TsvDiagnostics`) and in-memory (`MemoryDiagnostics`) diagnostic rows |
| `ldsc._runtime.artifacts` | `ArtifactFamily` declarations shared by workflow preflight, final writes, metadata file entries and cleanup; `declared_member()`, the one reader of metadata `files` entries (a relative name without `..`, with symbolic links either allowed anywhere, for LD-score directories and the saved-artifact preflight, or kept inside the directory, for `declared_result_file()`, which `h2_scale`, plotting and `infer-genome-build` use); plus `atomic_replace()`, the one temporary-sibling-then-rename publication, which the atomic table and JSON writers, the failure marker and the legacy summary-statistics audit use |
| `ldsc._runtime.progress` | phase boundaries and time-throttled progress through existing logger routes; bounded counters, current objects, and elapsed time |
| `ldsc._runtime.workflow_logging` | shared workflow logging context, `LDSC` logger level handling (`LOG_LEVEL_CHOICES` and `normalize_log_level()`, the one log-level check that `GlobalConfig` also uses), lifecycle audit lines, CLI console-handler routing (file-authoritative, console error-only) with run-aborting traceback capture, durable authorized-overwrite failure markers without rollback, and log-only formatting helpers |
| `ldsc._runtime.cli_help` | shared `CLIHelpFormatter` keeps flag names and paths intact when help descriptions wrap; `CHROMOSOME_PATH_HELP` and `SCALAR_PATH_HELP` supply verified path syntax; `add_output_dir_argument()` and `add_log_level_argument()` add the shared `--output-dir` and `--log-level` options with each command's help text; workflow parsers own their descriptions and argument groups |
| `ldsc._runtime.workflow_kwargs` | validate Python keyword arguments for `run_ldscore`, `run_build_r2_matrix`, `run_h2`, `run_partitioned_h2`, `run_rg`, and `run_quantile_h2` against the command parser before any write; fill the options a keyword call may not pass (SNP identifier mode, genome build, log level) from the registered `GlobalConfig`, the one place these wrappers read it; record the explicitly supplied options (`_explicit_cli_options`) for the keyword route and, through `parse_command_line()`, for the command lines of these six commands; run the regression and post-fit workflows under a failure marker naming the Python function |
| `ldsc._runtime.parallelism` | shared `workers` and `query_batch_size` validation and their argparse types, joblib detected CPU budgets, task caps, execution summaries, chromosome task scheduling, and the spawn worker pool that turns a crashed worker into `LDSCInternalError` |
| `ldsc._runtime.worker_runtime`, `ldsc._runtime.serial_calculation` | fixed one-BLAS-thread worker policy, runtime verification, independent Arrow limits, and controlled standalone numerical calls |
| `ldsc._genome.build_inference` | public `chr_pos` build and coordinate-basis inference helpers |
| `ldsc._genome.build_policy`, `ldsc._genome.build_inspection` | the build policy workflows apply to their inputs and read-only build inspection of saved artifacts (`infer-genome-build`) |
| `ldsc._genome.build_evidence` | build-evidence readers that the build policy and build inspection share: chunked table reads, delimiter and `# ldsc:` comment sniffing of user-supplied text tables, the bounded per-chromosome coordinate sampler, and R² matrix file and build-child discovery |
| `ldsc._genome.liftover` | liftover mechanics shared by munging and R² matrix building: hg19/hg38 chain-file translation, drop-all coordinate collision helpers, and readable drop reports |
| `ldsc._genome.hm3` | public packaged curated HM3 map loader and installed map path helper for workflow internals |
| `tools/hm3/build_hm3_chr_pos_reference.py` | maintenance tool (outside the package) that rebuilds the compact HM3 coordinate reference used by genome-build inference |
| `ldsc._annotation.builder` | public annotation builder and general standalone parser/API dispatch seam |
| `ldsc._annotation.annotate`, `ldsc._annotation.interval_queries` | standalone BED/gene gates and incremental projection |
| `ldsc._annotation.interval_projection` | `ChromosomeProjector` maps BED/gene interval unions to chromosome SNP positions; called by `interval_queries.build_query_shards()` after shared BED parsing |
| `ldsc._annotation.bundle`, `ldsc._annotation.storage` | explicit dataset ownership, chromosome descriptors, and detached selected reads from packed-binary and dense-float32 column stores |
| `ldsc._annotation.sources`, `ldsc._annotation.identity_cleanup`, `ldsc._annotation.loading` | bounded input scans, aligned source preparation, global identity cleanup (`stage_retained_rows`, shared with BIM grids), the one conversion of prepared sources into a bundle (`PreparedAnnotationSources.bundle`), and shared source-backed construction |
| `ldsc._annotation.parsing` | annotation chunk normalization and numeric-value validation shared by every annotation reader |
| `ldsc._annotation.preflight` | annotation declaration resolution and `inspect_annotation_inputs`, the declared-input gate with annotation header checks first, whose `AnnotationInputs` wraps the gate's `DeclaredInputs` with the annotation paths, `@` declarations, header widths and columns; the query checks that `annotate`, `ldscore` and `build_annotation_shards` share: reserved and colliding BED and gene-list query names, the gene-list Gate A error and the projection genome build |
| `ldsc._annotation.saved_results` | `AnnotationResult`, the record `run_annotate()` returns: saved baseline/query descriptors with no original-input dependencies, construction-owner detachment, and explicit preparation with `prepared_in(directory)` |
| `ldsc._annotation.bim_grid` | metadata-only standalone BIM grid scans using shared identity cleanup; no genotype companions or value synthesis |
| `ldsc._annotation.gene_lists` | required/embedded gene-coordinate-catalog validation, whose `GeneCatalogValidationError` names the first defects (`catalog_defect_message`) and carries the complete issue table for the consuming workflow's writer, and vectorized per-chunk query/control resolution (`resolve_gene_rows`) into exact IDs, intervals and audit rows |
| `ldsc._annotation.gene_query_storage` | gene-list Gate A (`resolve_gene_lists_staged`: bounded resolution, the one cross-chunk duplicate pass, summaries, and the required query-name check), staged per-source gene selections, complete audit replay with shared gene-coordinate-catalog state, and the shared gene coverage check (`with_coverage`) |
| `ldsc._annotation.query_outcomes` | shared resolution/support precedence, required-control and all-skipped failures, ordered diagnostics, and final query statuses from counts and score variation (`finalize_query_measurements`); the LD-score result records prune skipped queries themselves; does not import workflow orchestration |
| `ldsc._annotation.in_memory` | small prepared annotations and diagnostics without filesystem writes |
| `ldsc._annotation.semantics` | global annotation-name uniqueness and advisory binary/continuous classification, combined across chromosomes by `combine_annotation_types` |
| `ldsc._annotation.overlap` | annotation counts, overlap blocks and classification over reference SNPs without LD (`annotation_statistics`, `OverlapContribution`, `sum_overlap_contributions`), shared by direct `ldscore` and the gene-index build |
| `ldsc._annotation.outputs` | `AnnotationDirectoryWriter` for the `ldsc annotate` result directory |
| `ldsc._refpanel.plink` | PLINK genotype reader (`PlinkBEDFile`): selective SNP reads with fused individual filtering retain only selected packed genotypes; unrestricted panel builds opt into disk streaming instead |
| `ldsc._refpanel.r2_format` | the R² matrix format: file names (`r2_matrix_file_paths`, `r2_matrix_chromosomes`, `resolve_r2_build_dir`), schema constants, Parquet and sidecar writers, `R2MatrixChromosome` (footer contract, strict sidecar reading, sidecar binding, LD-window and truncation checks for `ldscore` and `query-r2-matrix`), `SortedR2BlockReader` and the pair lookup used by `query-r2-matrix` |
| `ldsc._refpanel.windows` | LD-window geometry shared by `ldscore`, `build-r2-matrix` and gene-index construction: `window_coordinates`, `get_block_lefts`, the position-order and CM checks, and the whole-chromosome guard (`check_whole_chromosome_window`, `WholeChromosomeWindowError`) |
| `ldsc._refpanel.prepared` | the `PreparedChromosome` contract, the chromosome annotation input `ChromosomeAnnotations`, and the `ArrayAnnotations`/`MappedAnnotations` accessors |
| `ldsc._refpanel.funnel` | the reference-SNP funnel: `ReferenceSNPFunnel` records what each reference-preparation step removed, with readers, writers and the run-level warning text |
| `ldsc._refpanel.plink_source` | `inspect_plink_inputs`, the PLINK content check of every PLINK route (BIM/FAM parsing, whole-number positions, BED header and size, chromosome assignment from `.bim` contents); `select_plink_chromosome`, the one BIM-level selection stage of `ldscore --plink-prefix`, gene-index builds and `build-r2-matrix`: a chromosome's `.bim` rows, the SNP restriction, identity cleanup of the restricted rows, the keep list mapped to `.fam` rows, and the step counts (`PlinkChromosomeSource`); the restriction match `build_restriction_mask`, which the R²-matrix adapter also applies to sidecar rows; and `build_plink_metadata_frame`/`orientation_flip_sign`, the retained SNP metadata with the minor allele as `A1`; kept apart from the GPL-derived reader in `plink` |
| `ldsc._refpanel.genetic_map` | genetic maps for `ldscore`, gene-index builds and `build-r2-matrix`: the position and cM column aliases (`GENETIC_MAP_POS_SPEC`, `GENETIC_MAP_CM_SPEC`), which add the HapMap names `Position(bp)` and `Map(cM)`, the declared-input header check (`GENETIC_MAP_HEADER_CHECKS`), map loading and combining (`load_genetic_map_group`) and cM interpolation (`interpolate_genetic_map_cm`); each command keeps its own map policy |
| `ldsc._refpanel.adapters` | runtime PLINK/R²-matrix adapters; `RefPanel.prepare_chromosome` owns reference filtering, annotation alignment, window/bias policy, and reader lifetime |
| `ldsc._refpanel.builder`, `ldsc._refpanel.build_kernel` | R²-matrix build workflow (source-build inference, optional liftover, explicit SNP/sample restrictions, optional `r2_min` threshold) and its kernel (pairwise R² emission in the exact LD window and genotype re-orientation; `build-r2-matrix` code only) |
| `ldsc._refpanel.query` | public `query-r2-matrix` CLI/API, `R2Matrix`, one-shot `query_r2_matrix()`, endpoint key resolution, sign harmonization, and optional unbiased-R²-to-Pearson-r conversion |
| `ldsc._refpanel.preflight` | `inspect_r2_headers` (each R² matrix chromosome's files, footer and sidecar header, for either reading command) |
| `ldsc._refpanel.outputs` | `R2MatrixBuildDirectoryWriter` for the `build-r2-matrix` output layout (owned files, chromosome-scoped names, metadata, dropped-SNP tables) and `QueryR2MatrixDirectoryWriter` for the `query-r2-matrix` result directory |
| `ldsc._sumstats.munger` | raw-summary-statistics CLI/API orchestration (`SumstatsMunger`), the `--infer-only` and `infer_raw_sumstats()` reports rendered from the munging plan, diagnostics under `diagnostics/`, and the run summary |
| `ldsc._sumstats.table` | the munged artifact: `SumstatsTable`, `load_sumstats()` for munged Parquet and legacy `.sumstats(.gz)` files, and the `ldsc:*` Parquet footer codec; regression imports it without the munging workflow |
| `ldsc._sumstats.outputs` | `curated_sumstats_frame()`, the one output frame of a run (fixed columns, canonical `CHR` labels, genomic row order) that the returned `SumstatsTable` and both data files hold; `SumstatsDirectoryWriter`: the owned munging output family (`artifact_family()`), staged Parquet and `.sumstats.gz` writing with chromosome row groups, and the dropped-SNP table |
| `ldsc._sumstats.munge_input` | private workflow helper for the munging plan (`resolve_munge_plan`, `MungePlan`): build and liftover option checks (`resolve_source_config`), raw schema, DANER/sample-size settings, the keep-list path, `INC_ALLELE`/`DEC_ALLELE` orientation, bounded source-build evidence and the liftover request; `prepare_munge_input` raises the munging plan's problems and loads the keep-list into `ResolvedMungeInput`, reading the packaged HM3 map once for HM3 quick liftover |
| `ldsc._sumstats.munge_kernel` | the raw reader `read_raw_sumstats` that the munging plan shares, chunk QC and restriction (including the rsID-family drop of unsupported `CHR` labels), whole-table N and sign conversion, resolved liftover, global identity cleanup and the stops when a stage retains no row; returns `MungeResult` with counts and provenance |
| `ldsc._sumstats.liftover` | summary-statistics liftover: `SumstatsLiftoverRequest`, HM3 quick liftover, the chain applier, liftover drop rows and the all-rows-dropped error |
| `ldsc._ldscore.compute` | LD-score math over a `PreparedChromosome`: `compute_chromosome`, R² pair streaming, PLINK block sums (`ld_scores_from_genotypes`), the float64 `ProjectionAccumulator`, counts and the regression-SNP mask; no legacy artifact emitters |
| `ldsc._ldscore.workflow` | the `ldscore` command and `run_ldscore()`: the parser, the typed run boundary (`_configs_from_args()` builds every configuration before the output preflight; `_run_direct_ldscore()` and `_run_explicit_indexed_ldscore()` run them), genome-build and gene-coordinate-catalog build selection, regression SNP and region loading, direct prepared-SNP support measurement (Gate B, counted with `ChromosomeProjector.counts`), and optional synthetic `base` |
| `ldsc._ldscore.calculator` | the public LD-score calculation service `LDScoreCalculator` and its `_RunContext`: build resolution and R² matrix checks for prepared inputs, logs and the chromosome pool, the direct LD-score producer's query batches, and output routing through shared query finalization |
| `ldsc._ldscore.preflight` | validates direct annotation/reference contents and chromosome equality from the run's `AnnotationBuildConfig` and `RefPanelConfig`; enforces complete `@` declarations, authoritative glob matches, and consolidated required-input diagnostics |
| `ldsc._ldscore.direct_annotation`, `ldsc._ldscore.direct_batches` | direct BED/gene query integration into an LD-score run; the direct LD-score producer: `_DirectProducer`, the one per-chromosome task `compute_direct_chromosome` of every query batch, staged direct query scores, and chromosome aggregation |
| `ldsc._ldscore.batch_execution` | the query-batch loop shared by the direct and indexed LD-score producers: one batch in flight through the run's persistent pool, failure shutdown before scratch removal, query-status finalization, and the gene-list viability gate after the last batch |
| `ldsc._ldscore.results` | the LD-score result records: the public materialized `LDScoreResult`, its private per-chromosome `ChromLDScoreResult`, and the private file-backed `DirectChromosomeFragment` and `StagedLDScoreBatch`, whose chromosome fragments a direct or indexed LD-score producer owns until the writer has consumed them; imported by the LD-score producers, the writer, regression and the legacy converter |
| `ldsc._ldscore.diagnostics`, `ldsc._ldscore.query_variation` | chromosome logs (and the main `ldscore.log` unless the calling command says with `main_log_open=True` that it already opened it), batch context and task/pool failure reporting, and the write-log-raise sequences both LD-score modes use for rejected gene lists and unusable queries; exact query-column variation across chromosomes |
| `ldsc._ldscore.source` | canonical directory loading (`read_ldscore_manifest()` for the manifest alone, `load_ldscore_from_dir()` for the tables), manifest/schema and table/identity validation shared with writers, the shared schema facts (weight column, count-record fields and `count_record()`, `count_config` fields, index provenance fields), resident baseline/count/overlap metadata, and uncached selected-query reads |
| `ldsc._ldscore.overlap_table` | `LDScoreOverlap`, its query selection and batch joining (`select_queries()`, `concat_query_batches()`), and its long-frame Parquet read/write for the LD-score directory's overlap table |
| `ldsc._ldscore.outputs` | `LDScoreDirectoryWriter`: LD-score Parquet layout with chromosome row groups, query batch files, query diagnostics and the root `metadata.json`, and `write_batches`, the one staged publication of complete LD-score directories, which reads a staged batch's fragments one chromosome query frame at a time against the shared baseline |
| `ldsc._ldscore.legacy_converter` | sole legacy LDSC LD-score-suite import boundary: deterministic family discovery, rsID joins, count/overlap validation or reconstruction, source/discovery provenance without source-file hashes, diagnostics, and LDSC3 LD-score directory writing |
| `ldsc._ldscore.gene_index.build` | gene LD-score index building and transactional publication |
| `ldsc._ldscore.gene_index.artifact` | the on-disk index format, writing and reading: chromosome payloads, `index_id` and published-row hashes, staged writes, component inspection, `load_gene_ldscore_index` and `LoadedGeneLDScoreIndex` |
| `ldsc._ldscore.gene_index.scoring`, `ldsc._ldscore.gene_index.batches` | `run_indexed_ldscore` and bounded chromosome workers holding one operator each |
| `ldsc._ldscore.gene_index.atoms` | sparse disjoint-atom primitives for exact gene LD-score indexes |
| `ldsc._regression.results` | the public regression dataset and result types (`RegressionDataset`, `RegressionFit`, `H2Result`, `PartitionedH2BatchResult`, `RgResultFamily`) and the `PartitionedH2FitArtifacts` the writers publish |
| `ldsc._regression.alignment` | file-driven regression dataset assembly: per-trait preparation once per run (`_PreparedTraitSumstats`: automatic legacy summary-statistics rsID-to-panel projection, the missing-or-infinite `Z`/`N` drop, and keyed rows cached per SNP identifier mode), the genome-build check, per-command allele alignment, active effective identity-key merging (`SNP`, `SNP:<allele_set>`, `CHR:POS`, or `CHR:POS:<allele_set>`), model column and count selection, the preparation spine (`_prepare_h2_dataset()`, `_h2_dataset()`, `_prepare_rg_dataset()`), and the run log's SNP identifier line from the prepared datasets (`_log_regression_identity()`) |
| `ldsc._regression.fitting` | the fit policy, the `h2` and `rg` fits (`_fit_h2_dataset()`, `_fit_rg_dataset()`) and the per-query `partitioned-h2` fit (`fit_partitioned_query()`), the collinearity guard, observed/liability-scale summary columns, exact final-fit `h2` regression-bin diagnostics, result rows and metadata |
| `ldsc._regression.kernel_values` | the one replacement of the estimator kernel's `'NA'` sentinels with NaN, where the fits and the overlap-aware table receive kernel results, and the `RG` flag `h2_out_of_bounds` |
| `ldsc._regression.partitioned`, `ldsc._regression.rg` | the `partitioned-h2` batch over the two overlap-aware regimes, which publishes its result directory, and the `rg` pair run |
| `ldsc._regression.commands` | the `h2`, `partitioned-h2` and `rg` command boundary: argument registration, `_configs_from_args()`, which turns parsed arguments into a `RegressionConfig` before the output preflight, the typed bodies `_run_h2()`, `_run_partitioned_h2()` and `_run_rg()`, and the Python forms `run_h2()`, `run_partitioned_h2()` and `run_rg()` |
| `ldsc._regression.runner` | the public `RegressionRunner` facade, whose methods log the `GlobalConfig` banner when they read inputs and forward its configs to `alignment`, `fitting`, `partitioned` and `rg` |
| `ldsc._regression.partitioned_parallel` | private read-only regression maps, bounded spawned whole-query workers, and worker termination/log forwarding; statistical calculations remain in `_regression.fitting.fit_partitioned_query()`, which the `partitioned-h2` batch passes to the pool |
| `ldsc._regression.estimators` | LDSC estimators for `Hsq` and `RG` |
| `ldsc._regression.jackknife`, `ldsc._regression.irwls` | supporting numerical routines used by regression |
| `ldsc._regression.overlap_model` | per-model overlap assembly, the overlap-aware annotation table (ported `_overlap_output` with augmentation, renamed once to the `partitioned-h2` column names) and the collinearity check (`model_collinearity_error`) |
| `ldsc._regression.prevalence` | parse the `rg` prevalence inputs (comma-separated `--samp-prev`/`--pop-prev` lists or a `--prevalence-manifest` TSV) into a normalized per-trait `(samp_prev, pop_prev)` structure for observed-to-liability conversion, checking each pair with `ldsc.config.normalize_prevalence_pair()`; `RegressionConfig` checks the scalar `h2`/`partitioned-h2` pair |
| `ldsc._regression.h2_scale` | strict post-fit observed-to-liability conversion from a canonical `h2` result, including exact and prevalence-range modes and the fixed nested derived-result family, declared as an `ArtifactFamily` |
| `ldsc._regression.quantile_h2`, `ldsc._regression.quantile_inputs`, `ldsc._regression.quantile_storage` | post-fit continuous-target quantile assignment, reconstruction and verification of the fitted sources and common reference SNPs, vectorized coefficient/delete-value projection, standardized `tau_star`, and streamed alignment diagnostics |
| `ldsc._regression.outputs` | writers for the `h2`, `partitioned-h2` (per-query layout), `rg` (result family) and `quantile-h2` result directories, including the legacy summary-statistics audit, with their column constants, the fitted-model file names and locators (`fitted_model_files()`, `published_query_artifacts()`), and `_publish_staged_tree()`, the staged publication and rollback of the per-query and pair trees |
| `ldsc._plotting.dispatch` | sole public metadata-driven plotting dispatcher, lazy plotting-runtime boundary, fixed plot-family output declared as an `ArtifactFamily`, and live `PlotArtifact` return object |
| `ldsc._plotting.builders` | private Matplotlib-only headless builders for the approved `h2`, `rg`, `partitioned-h2`, and `quantile-h2` plots |

## Where To Change Code

| Goal | Start here |
| --- | --- |
| change CLI flags or subcommand wiring | workflow `build_parser` / `add_*_arguments` functions, then the `_COMMANDS` table in `src/ldsc/cli.py`, which `main` and `_command_parser` read; follow [CLI help guidelines](guidelines/cli-help.md) |
| change help wrapping, the shared `--output-dir`/`--log-level` options or shared log-level wording | `src/ldsc/_runtime/cli_help.py:CLIHelpFormatter`, `add_output_dir_argument()`, `add_log_level_argument()`, `src/ldsc/_runtime/workflow_logging.py:LOG_LEVEL_HELP`; verify `tests/package/test_cli_help.py` |
| change path-token behavior | `src/ldsc/_runtime/paths.py` |
| change output collision policy | `src/ldsc/_runtime/paths.py`, then the workflow writer that owns the artifact |
| change a keyword-argument Python workflow (`run_ldscore`, `run_build_r2_matrix`, `run_h2`, `run_partitioned_h2`, `run_rg`, `run_quantile_h2`) | the function in its workflow module, keyword checks and registered-`GlobalConfig` values in `src/ldsc/_runtime/workflow_kwargs.py`; verify `tests/package/test_workflow_keyword_checks.py`, `tests/package/test_global_config_registry.py` and `tests/regression/test_python_regression_workflows.py` |
| change workflow log files, lifecycle lines, or log-level handling | `src/ldsc/_runtime/workflow_logging.py`, then the workflow orchestration function that owns the output directory |
| change header aliases or identifier/build normalization | `src/ldsc/_core/columns.py` |
| change read-only artifact build inspection | `src/ldsc/_genome/build_inspection.py`, which also holds the command's argument and stdout report, with its read-only `_COMMANDS` row in `src/ldsc/cli.py` and coverage in `tests/genome/test_genome_build_inspection.py` |
| change automatic `chr_pos` genome-build inference | `src/ldsc/_genome/build_inference.py` |
| change annotation loading or `ldsc annotate` behavior | `src/ldsc/_annotation/builder.py`, then `src/ldsc/_annotation/sources.py` and `src/ldsc/_annotation/annotate.py` |
| change BED parsing or interval projection | `src/ldsc/_core/regions.py:iter_bed_rows`, then `src/ldsc/_annotation/interval_queries.py` and `src/ldsc/_annotation/interval_projection.py:ChromosomeProjector` |
| change gene coordinate catalog validation, identifier resolution, or list naming | `src/ldsc/_annotation/gene_lists.py`, then `src/ldsc/_annotation/builder.py` |
| change BED/gene partial-success statuses or diagnostic schemas | `src/ldsc/_annotation/query_outcomes.py`, `src/ldsc/_ldscore/batch_execution.py` (`finalize_batch`), `src/ldsc/_ldscore/results.py` (`select_queries`), `src/ldsc/_ldscore/workflow.py` (Gate B), `src/ldsc/_ldscore/outputs.py` |
| change R²-matrix build logic | `src/ldsc/_refpanel/builder.py`, then `src/ldsc/_refpanel/build_kernel.py` |
| change named regression-SNP region exclusion for `ldscore` or gene LD-score index construction | `src/ldsc/config.py`, `src/ldsc/_core/regions.py`, then `src/ldsc/_ldscore/workflow.py` or `src/ldsc/_ldscore/gene_index/build.py` |
| change R² pair lookup CLI/API behavior | `src/ldsc/_refpanel/query.py`, then `src/ldsc/_refpanel/r2_format.py:lookup_pairs_in_parquet` |
| change runtime PLINK/R²-matrix reference access | `src/ldsc/_refpanel/adapters.py`, with the R² file format in `src/ldsc/_refpanel/r2_format.py` and LD windows in `src/ldsc/_refpanel/windows.py` |
| change which PLINK rows or individuals a chromosome keeps (restriction, identity cleanup, keep list) | `src/ldsc/_refpanel/plink_source.py`, shared by `ldscore --plink-prefix`, gene-index builds and `build-r2-matrix` |
| change how PLINK inputs are discovered or validated before any computation | discovery in `src/ldsc/_runtime/paths.py:discover_plink_members`, content checks in `src/ldsc/_refpanel/plink_source.py:inspect_plink_inputs` |
| change `ldscore` options, input modes, or optional-baseline behavior | `src/ldsc/_ldscore/workflow.py`, then `src/ldsc/_ldscore/direct_annotation.py` |
| change LD-score orchestration or output packaging | `src/ldsc/_ldscore/calculator.py`, `src/ldsc/_ldscore/outputs.py` |
| change the direct per-chromosome LD-score task, chromosome aggregation, or the LD-score result records | `src/ldsc/_ldscore/direct_batches.py`, then `src/ldsc/_ldscore/results.py` and `src/ldsc/_ldscore/batch_execution.py` |
| change saved LD-score loading or selected-query reads | `src/ldsc/_ldscore/source.py`; keep model-specific validation in `src/ldsc/_regression/alignment.py`; `quantile-h2` reads the manifest through `read_ldscore_manifest()` |
| change LD-score math | `src/ldsc/_ldscore/compute.py` |
| change raw summary-statistics ingestion, format inference, `CHR`/`POS` handling, or summary-statistics SNP keep-list filtering | `src/ldsc/_sumstats/munger.py`, then `src/ldsc/_sumstats/munge_input.py` and `src/ldsc/_sumstats/munge_kernel.py` |
| change summary-statistics liftover or its rows in the dropped-SNP table | `src/ldsc/_sumstats/liftover.py`; chain translation shared with `build-r2-matrix` is in `src/ldsc/_genome/liftover.py` |
| change munged output files, footer provenance, or munged loading | `src/ldsc/_sumstats/outputs.py` and `src/ldsc/_sumstats/table.py` |
| change regression dataset assembly, fits or CLI summaries | `src/ldsc/_regression/alignment.py` (datasets), `src/ldsc/_regression/fitting.py` (fits, summaries, metadata), `src/ldsc/_regression/partitioned.py` or `src/ldsc/_regression/rg.py` (batch and pair runs), `src/ldsc/_regression/commands.py` (options and command runs), then `src/ldsc/_regression/outputs.py`, `docs/reference/results/partitioned-h2.md` for `partitioned-h2` output layout, `docs/reference/formats/artifact-metadata.md` for `rg` output contracts, and `docs/reference/regression-options.md` for the tunable estimator parameters and defaults |
| change annotation-name validation or continuous-annotation classification | `src/ldsc/_annotation/semantics.py`, then `src/ldsc/_ldscore/workflow.py`, `src/ldsc/_ldscore/calculator.py` and `src/ldsc/_ldscore/outputs.py` |
| change target quantiles, standardized coefficients, or post-fit alignment | `src/ldsc/_regression/quantile_h2.py`, then `src/ldsc/_regression/outputs.py` and `docs/reference/results/quantile-h2.md` |
| change LDSC estimators | `src/ldsc/_regression/estimators.py` |
| change binary-trait prevalence parsing or observed-to-liability conversion | `src/ldsc/config.py` (`validate_prevalences()`, the rule, and `normalize_prevalence_pair()`, one trait's pair) and `src/ldsc/_regression/prevalence.py` (`rg` list and manifest parsing), then `src/ldsc/_regression/fitting.py` (summary/metadata wiring) and `src/ldsc/_regression/estimators.py` (`liability_conversion_factor`) |
| change post-fit heritability scale conversion or sensitivity output | `src/ldsc/_regression/h2_scale.py`, then `src/ldsc/_regression/estimators.py` only if the shared numerical primitive itself changes |
| change plot dispatch, supported result contracts, or visual semantics | `src/ldsc/_plotting/dispatch.py`, `src/ldsc/_plotting/builders.py`, then `docs/dev/design/plotting.md` and `docs/user-guide/workflows/plotting.md` |
| change owned artifact declarations | `src/ldsc/_runtime/artifacts.py`, then the domain writer: `src/ldsc/_annotation/outputs.py`, `src/ldsc/_refpanel/outputs.py`, `src/ldsc/_sumstats/outputs.py`, `src/ldsc/_ldscore/outputs.py` or `src/ldsc/_regression/outputs.py` |
| change LD-score directory files, Parquet row-group layout, or LD-score metadata | `src/ldsc/_ldscore/outputs.py` |
| change the `partitioned-h2` per-query layout, the `rg` result-family layout, or regression metadata payloads | `src/ldsc/_regression/outputs.py` |
| change which module may import which, or defer an import | the tier table above and `tests/package/test_layering.py` (`DEFERRED_IMPORTS` for imports inside functions) |

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
  to `run_annotate_from_args()`; `_annotation.parsing` must not own parser
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
  interprets as scientific data artifacts. The gene LD-score index
  builder keeps its open workflow log in the index's hidden `.build/` directory and moves
  only the closed successful log into the published diagnostics; nothing is written beside the index.
- Keep regression file-driven: it should be able to rebuild state from written artifacts without recomputing LD scores.
- Keep plotting metadata-driven and post-fit: plot builders consume declared canonical tables and never refit. Matplotlib is required at installation but lazily imported through `ldsc._plotting.builders`.
- Prefer extending shared helpers or the workflow-owned writer over duplicating local parsing or writing logic.

## Test Map

Each test file sits in the folder of the domain package whose code it mainly exercises, such as `tests/ldscore/` for `_ldscore`.
Tests of shared machinery that every workflow uses, such as output writers, workflow logs, failure markers, workers and input preflight, are in `tests/runtime/`.
Tests of the root shell and of the package as a whole are in `tests/package/`.
Shared helpers and `fixtures/` stay at the top of `tests/`; tests import them, and one another, by full module path, for example `from tests.annotation.test_annotate_streaming import inputs`.
Gene LD-score index fixtures are written by `tests/gene_index_helpers.py` through the build's own staging and commit steps; a test that builds or loads an index covering fewer than autosomes 1 to 22 narrows the index chromosome set with its `narrow_index_chromosomes`.
LD-score directory fixtures are written by `tests/ldscore_helpers.py` through `LDScoreDirectoryWriter.write_batches`, so they carry the writer's real manifest; a test of a rejected directory starts from that output and changes it.
The four-SNP quantile-h2 problem, with its writer-built LD scores, hand-written fitted model and source tables, comes from `tests/quantile_fixtures.py`.
Each test folder has an `__init__.py`, which makes those full module paths importable.

| Area | Main tests |
| --- | --- |
| output layer | `tests/runtime/test_output.py` |
| annotation workflow | `tests/annotation/test_annotation.py`, `tests/annotation/test_annotate_grid.py`, `tests/annotation/test_annotate_streaming.py` |
| reference-panel builder | `tests/refpanel/test_r2_matrix_builder.py` |
| runtime reference panels | `tests/refpanel/test_ref_panel.py` |
| R² pair query | `tests/refpanel/test_r2_matrix_query.py` |
| R² matrix chromosome reads and LD windows | `tests/refpanel/test_r2_format.py`, `tests/refpanel/test_windows.py`; also `tests/refpanel/test_plink_io.py` (row remap, pair iteration, `get_block_lefts`), `tests/refpanel/test_r2_parquet_contract.py` (file contract and the int16 scale) and `tests/refpanel/test_r2_matrix_chromosome.py` (sidecar damage through `ldscore` and `query-r2-matrix`) |
| LD-score workflow | `tests/ldscore/test_ldscore_workflow.py` |
| gene coordinate catalog and gene-list resolution | `tests/annotation/test_gene_list_resolver.py` |
| summary-statistics munging | `tests/sumstats/test_sumstats_munger.py`, `tests/sumstats/test_sumstats_outputs.py`, `tests/sumstats/test_sumstats_liftover.py` |
| regression workflow | `tests/regression/test_alignment.py`, `tests/regression/test_fitting.py`, `tests/regression/test_partitioned.py`, `tests/regression/test_rg.py`, `tests/regression/test_regression_commands.py` |
| regression estimators | `tests/regression/test_irwls.py`, `tests/regression/test_jackknife.py` |
| path and config contracts | `tests/runtime/test_path_resolution.py`, `tests/core/test_config_identifiers.py`, `tests/package/test_global_config_registry.py`, `tests/core/test_column_inference.py`, `tests/core/test_chromosome_inference.py`, `tests/genome/test_genome_build_inference.py` |
| workflow logs | `tests/runtime/test_logging_refactor.py` |
| package layout and import tiers | `tests/package/test_layering.py`, `tests/package/test_package_layout.py` |
| region exclusion | `tests/core/test_regions.py` |

`AnnotationBundle.validate()` validates row alignment and annotation names; identity cleanup occurs during builder construction. It no longer accepts an ignored identifier argument. Column-name inference uses `infer_chr_pos_columns`; numerical window construction uses `get_block_lefts`. Their former forwarding spellings are removed without changing accepted input column aliases or legacy formats.

`ldsc._runtime.artifacts` reads declared result files (`declared_result_file()`) beside the `ArtifactFamily` declarations that write them, and publishes JSON for `h2_scale` and plotting as for every writer.
Each workflow retains its metadata interpretation, table loading and scientific dispatch.

## Domain contracts and decision history

| Change | Current contract | Main owner |
| --- | --- | --- |
| annotations and interval projection | [BED input](../reference/formats/bed-intervals.md), [gene lists](../reference/formats/gene-lists.md), [query repair](../reference/results/gene-list-diagnostics.md) | `_annotation.builder.AnnotationBuilder`; `_annotation.interval_queries.build_query_shards`; `_annotation.interval_projection.ChromosomeProjector` |
| panel preparation, CM/MAF and traversal | [R² pipeline](../reference/formats/r2-matrices.md), [pair accumulation](../reference/methods/ld-score-computation.md), [window behavior](../reference/formats/r2-matrices.md#37-ld-windows-and-sidecar-values-in-ldscore) | `_refpanel.adapters.RefPanel.prepare_chromosome`; `_ldscore.compute.compute_chromosome` |
| regression SNPs and region masks | [reference and regression SNPs](../reference/ldscore-snp-universes.md), [region exclusions](../reference/formats/region-exclusion-presets.md) | `_ldscore.workflow._regression_region_intervals`; `_ldscore.compute.regression_mask_from_keys`; `_ldscore.gene_index.build` resolves indexed policy |
| reference construction and pair queries | [pair query](../reference/query-r2-matrix.md), [path contract](../reference/input-paths.md) | `_refpanel.builder.R2MatrixBuilder`; `_refpanel.query.R2Matrix`; `_refpanel.outputs.QueryR2MatrixDirectoryWriter` |
| munging and liftover | [munging](../reference/munge-sumstats.md), [liftover](../reference/snp-identity-and-genome-build.md#liftover) | `_sumstats.munger.SumstatsMunger`; `_sumstats.munge_input.resolve_munge_plan`, `.resolve_raw_schema`, `.prepare_munge_input`; `_sumstats.munge_kernel`; `_sumstats.liftover.apply_sumstats_liftover`; `_sumstats.outputs.curated_sumstats_frame`, `.SumstatsDirectoryWriter` |
| regression and overlap | [configuration](../reference/regression-options.md), [partitioned results](../reference/results/partitioned-h2.md) | `_regression.runner.RegressionRunner`; `_ldscore.overlap_table`; `_regression.overlap_model`; `_regression.estimators` |
| indexed gene queries | [index](../reference/gene-ldscore-index.md), [mathematics](../reference/methods/gene-ldscore-index-math.md) | `_ldscore.gene_index.build`, `.artifact`, `.scoring`, `.batches`; `_ldscore.gene_index.atoms` |
| derived results | [plotting](design/plotting.md), [quantile projection](../reference/results/quantile-h2.md) | `_plotting.dispatch`, `_regression.h2_scale`, `_regression.quantile_h2`; shared `_runtime.artifacts` |

`query-r2-matrix` requires `--output-dir`; `QueryR2MatrixDirectoryWriter` writes its result directory.
Region presets affect the regression SNPs and their weights, not the reference SNPs.
`build-r2-matrix` has no region-pruning fields; use its explicit reference-SNP restriction when pruning is intended.
Reference preparation supplies population-specific CM/MAF; annotation CM is an `NA` placeholder.
LD scores from R² matrices use pair streaming, and PLINK filtering happens during adapter preparation.

Design rationale that protects a non-obvious constraint is stated in the page that documents the constraint; measurements are in [validation](validation/performance.md).
Earlier dated specifications, implementation plans, audits, and the archived design map were removed on 2026-10-26, and the decision records on 2026-10-27; they remain in git history, for example `git show 479fc4d:docs/specs/<file>` or `git show 77a54e2:docs/decisions/2026-10-11-parallel-execution.md`.
The specification and plan for the domain-package layout were removed once it was implemented on 2026-10-28; see `git show 869f41e:docs/specs/2026-10-27-domain-subpackages.md`.
