# Architecture

Last updated on: 2026-11-06

LDSC3 (distribution `ldsc`) is the Python 3 package that replaces legacy LDSC.
It reads optional SNP-level annotations, PLINK reference panels or R² matrices, and GWAS summary statistics; resolves user-facing path and header conventions in the public workflow layer; delegates numerical work to the computation modules of its private domain packages; and writes LDSC3 artifacts that can be chained into later runs.

Start with [the contributor entry](code-structure.md) to locate a change; this document explains the execution architecture. Related domain references:

- [`parallel-execution.md`](../reference/parallel-execution.md): worker budgets, one BLAS thread per worker, process ownership, and chromosome logs
- [`data-flow.md`](data-flow.md): user-visible file streams and flowcharts
- [`python-api.md`](../reference/python-api.md) and [`cli.md`](../reference/cli.md): public Python API and command-line options
- [`code-structure.md`](code-structure.md): module map and change guide
- [`logs-and-outputs.md`](../reference/logs-and-outputs.md): per-run log naming, preflight, and API boundaries
- [`plotting.md`](design/plotting.md): plotting, post-fit heritability scale conversion, dispatch, and derived-output ownership
- [`gene-lists.md`](../reference/formats/gene-lists.md): gene-list parsing, the gene coordinate catalog, projection, and diagnostics
- [`query-r2-matrix.md`](../reference/query-r2-matrix.md): pairwise R² lookup contract for R² matrices
- [`snp-identity-and-genome-build.md`](../reference/snp-identity-and-genome-build.md): genome-build policy, inference, and the shared liftover rules
- [`partitioned-h2.md`](../reference/results/partitioned-h2.md): partitioned-heritability result columns and interpretation
- [`quantile-h2.md`](../reference/results/quantile-h2.md): continuous-annotation post-fit contracts
- [`memory-and-cpu.md`](design/memory-and-cpu.md): consolidated memory and CPU efficiency designs, including bit-packing, genotype/R² streaming, sparse projection, index reuse, workers, regression copies/artifacts, exact quantiles, and evidence limits

## Bird's-Eye View

The [workflow map](../reference/index.md#workflow-map) lists every workflow with its command, Python entry points, main inputs and main outputs, and the [package overview](data-flow.md#package-overview) draws how their artifacts chain.
`munge-sumstats` writes munged summary statistics.
`annotate`, `build-r2-matrix` and `build-gene-ldscore-index` prepare optional inputs of `ldscore`, and `ldscore` or `convert-legacy-ldscores` writes an LD-score directory.
`h2`, `partitioned-h2` and `rg` read munged summary statistics and an LD-score directory and write result directories, which `quantile-h2`, `plot` and `convert-h2-scale` read without refitting.
`query-r2-matrix` looks up SNP pairs in an R² matrix, and `infer-genome-build` reports the genome build of one input without writing files.
Every artifact-writing workflow writes a per-run log inside its output directory ([logs and outputs](../reference/logs-and-outputs.md)).

## Layer Structure

- **CLI Layer**: public command dispatch in `ldsc.cli`
- **Workflow And Preprocessing Layer**: public services in `ldsc._annotation.builder`, `ldsc._refpanel.builder`, `ldsc._refpanel.query`, `ldsc._ldscore.workflow`, `ldsc._ldscore.calculator`, `ldsc._ldscore.legacy_converter`, `ldsc._ldscore.gene_index.build`, `ldsc._sumstats.munger`, the regression workflow modules (`ldsc._regression.commands`, the `RegressionRunner` facade in `ldsc._regression.runner`, and `ldsc._regression.alignment`, `ldsc._regression.fitting`, `ldsc._regression.partitioned` and `ldsc._regression.rg`), and `ldsc._regression.quantile_h2`, plus shared normalization in `ldsc.config`, `ldsc._runtime.paths`, `ldsc._core.columns`, `ldsc._core.chromosomes`, `ldsc._genome.build_inference`, `ldsc._annotation.semantics`, and the internal `ldsc._annotation.gene_lists`
- **Computation**: private file-format and numerical code kept with its domain rather than in one shared kernel layer: annotation value checks in `ldsc._annotation.parsing`, BED parsing in `ldsc._core.regions`, interval projection in `ldsc._annotation.interval_projection`, PLINK and R² matrix access in `ldsc._refpanel.plink` and `ldsc._refpanel.r2_format`, LD-score math in `ldsc._ldscore.compute`, legacy-compatible munging in `ldsc._sumstats.munge_kernel`, and the estimators in `ldsc._regression.estimators`
- **Shared Artifact Access**: LD-score directory opening, schema/identity validation, and selected-query reads in `ldsc._ldscore.source`, consumed by generation and regression
- **Output And Derived-Result Layer**: result writing in each domain's `outputs` module (`ldsc._annotation.outputs`, `ldsc._refpanel.outputs`, `ldsc._sumstats.outputs`, `ldsc._ldscore.outputs`, `ldsc._regression.outputs`) on the shared `ldsc._runtime.artifacts` base; metadata-driven plotting in
  `ldsc._plotting.dispatch`; post-fit scale conversion in `ldsc._regression.h2_scale`

### Ownership by layer and workflow

Each object's owner is recorded once, by module, in the [module map](code-structure.md#module-map) and the module docstring.
The commands, Python entry points, configuration and result objects of each workflow are in the [workflow map](../reference/index.md#workflow-map), the [Python API](../reference/python-api.md) and the [command-line reference](../reference/cli.md), and the files each workflow writes are in [public output layouts](../reference/formats/artifact-metadata.md#public-output-layouts) and [logs and outputs](../reference/logs-and-outputs.md).
The layers above cut across the domain packages, whose imports follow the one-directional [tiers](code-structure.md#dependency-direction) that `tests/package/test_layering.py` enforces.

### Ownership boundaries

The continuous-annotation extension keeps semantics and post-fit work in the public layers: `_annotation/semantics.py` owns advisory classification; routine annotation fingerprints are not computed, `_regression/fitting.py` (each query model) and `_regression/partitioned.py` (the baseline-only model) persist one delete matrix per fitted model, `_regression/quantile_h2.py` owns source reconstruction and projection, and `_regression/outputs.py` owns the result family. The regression kernel is unchanged.

Explicit result exploration follows the same direction.
`_plotting.dispatch` validates canonical result metadata and owns derived figure output, while `_plotting.builders` imports Matplotlib for plots; the only other importer is `_regression.h2_scale`, for its prevalence-range figure.
`_regression.h2_scale` loads a canonical `h2` result and reuses the regression kernel's conversion factor.
Neither path is called automatically by a core workflow.

Gene-list support follows the same boundary.
`ldsc._annotation.gene_lists` validates the required or embedded gene coordinate catalog and returns compact batch resolution/audit records; `ldsc._annotation.query_outcomes` defines internal BED/gene status records; `AnnotationBuilder` projects intervals; LD-score query-batch finalization (`_ldscore.batch_execution.finalize_batch()`) prunes zero-hit and zero-variance queries; and `LDScoreDirectoryWriter` owns the status and unresolved-gene diagnostics.
No gene identity enters the LD-score computation in `_ldscore.compute`.

Legacy compatibility follows the same boundary.
`ldsc._ldscore.legacy_converter` alone reads selected legacy LD-score/count/annotation/frequency suites and writes an LD-score directory, opened with `load_ldscore_from_dir()` as an `LDScoreSource`; it is exported as `convert_legacy_ldscores()`.
Regression never reads those fragments.
`ldsc._regression.alignment` automatically projects legacy summary-statistics text onto the regression SNPs of the already-loaded LD-score directory, once per trait per run, and owns the resulting drop audit.
Neither compatibility path adds a legacy reader or writer to the computation modules.

PLINK discovery lives in `_runtime.paths.discover_plink_members()`, which finds the selected BED/BIM/FAM trios and reports missing members; the declared-input gate calls it, so it stays with path resolution.
PLINK content checks and chromosome assignment live in `_refpanel.plink_source.inspect_plink_inputs()`, beside the PLINK chromosome selection stage.
It uses the low-level BIM/FAM readers for structural inspection (positions, the BED header and size), assigns chromosomes from `.bim` contents, returns combined issue records in a `_runtime.paths.PlinkInputResolution`, and supplies concrete chromosome-prefix mappings to direct LD-score preflight, gene LD-score index construction, and R²-matrix construction.
`PlinkRefPanel` consumes a supplied mapping or resolves it lazily once for direct Python use.
The kernel no longer owns a PLINK path-token wrapper.

## File Tree

The [repository map](code-structure.md#repository-map) shows the source and test layout with each package's import tier.

## Code Map

Each module's own responsibility is in the [module map](code-structure.md#module-map) and its module docstring, and its tests are in the [test map](code-structure.md#test-map).
This map keeps how the modules hand work to one another.

### Command line and shared contracts

`python -m ldsc` (`ldsc.__main__`) and the `ldsc` console script both call `ldsc.cli.run_cli()`, which wraps `main()` in the user-error boundary.
`ldsc.cli` dispatches each subcommand through its `_COMMANDS` row and imports the workflow module only when that command is parsed or run; `_Command` in `src/ldsc/cli.py` describes the row, the read-only flag and who prints successful-run notices.
Workflows resolve path tokens (`_runtime.paths`), header aliases (`_core.columns`), chromosome order and the genome build before they call a kernel, so numerical functions receive resolved inputs.
`_runtime.paths.output_directory_path()` is the one check that rejects an empty or blank destination: `ensure_output_directory()`, the failure-marker scope, `_runtime.scratch.ScratchWorkspace` and the `plots/` and liability-scale roots call it, and `_runtime.cli_help.output_directory_argument()` rejects a blank option while parsing ([output confinement](../reference/logs-and-outputs.md#output-confinement)).

`_runtime.preflight.inspect_declared_inputs()` is the declared-input gate that `annotate`, `ldscore`, `build-gene-ldscore-index`, `build-r2-matrix`, `query-r2-matrix`, `munge-sumstats`, the regression commands and `quantile-h2` run before any large scan.
Each domain plugs in its own checks: `_annotation.preflight.inspect_annotation_inputs()` adds annotation headers, `_runtime.paths.discover_plink_members()` checks PLINK companions, and callers add checks such as `inspect_artifact_paths()` for saved artifacts.
Content gates later in a run use `InputGate` directly; every gate builds its rows with `_runtime.input_issues.input_issue()` and raises one `LDSCInputError` after writing the caller's input-issue table ([staged input validation](../reference/logs-and-outputs.md#staged-input-validation)).

### Genome builds

Workflows and the reference-panel adapters call `_genome.build_policy.GenomeBuildCheck` once, before filtering or numerical work; numerical readers receive the concrete build and do not repeat the decision ([build policy](../reference/snp-identity-and-genome-build.md#genome-build-policy)).
`_genome.build_inspection`, behind the read-only `infer-genome-build`, shares the raw-evidence readers of `_genome.build_evidence` and the HM3 inference of `_genome.build_inference` with the build check but never calls workflow build resolution, opens no failure-marker scope and writes no artifact ([read-only inspection](../reference/logs-and-outputs.md#read-only-inspection-exception)).

### Workflow logs and failure markers

`_runtime.workflow_logging` owns the per-run log, the CLI console handler and the `RUN_FAILED` marker; their contract is in [logs and outputs](../reference/logs-and-outputs.md).
Each run has one marker scope, opened by its workflow (`materializing_overwrite_guard()` or `overwrite_failure_marker()`) from parsed arguments, or, for `build-r2-matrix`, by the `_build()` body that both entry points call, after `R2MatrixBuilder._inspect_inputs()` has assigned chromosomes from BIM contents.
A scope opened inside an open one neither writes nor removes, so the outermost public entry point decides the marker and its label.
The CLI opens no scope; `run_cli()` wraps dispatch in `cli_command()`, and the marker records that command line.
`serial_calculation()` returns the workflow log path opened in its worker, so a marker written in the calling process names that log.

### Plotting and scale conversion

`_plotting.dispatch` and `_regression.h2_scale` read canonical result directories through `_runtime.artifacts.declared_result_file()` and declare their derived families as `ArtifactFamily` values, as the domain writers do; their seams, output families and the `rg` heritability annotations are in [plotting](design/plotting.md).

### Annotations

`run_annotate()` (`_annotation.annotate`), `ldscore` and `_annotation.loading.build_annotation_shards()` (behind `AnnotationBuilder` and `quantile-h2`) share the query checks of `_annotation.preflight`: reserved or colliding BED and gene-list query names stop a run before any grid is read, and the gene-list Gate A error and the projection genome build come from one place, while `ldscore` resolves its own build against the reference panel.
Gate A itself is one call, `_annotation.gene_query_storage.resolve_gene_lists_staged()`, in `annotate`, `ldscore`, indexed scoring and `build_annotation_shards()`; each workflow keeps its own failure handling.
BED and gene intervals reach SNPs through one projector: `_annotation.interval_queries.build_query_shards()` applies `_annotation.interval_projection.ChromosomeProjector` to BED rows from `_core.regions.iter_bed_rows()` and to resolved gene intervals.
Only `annotate` measures gene support in that pass, by passing a `gene_support` object; direct `ldscore` measures it later on the prepared reference SNPs (Gate B).
`_ldscore.direct_annotation` joins the annotation modules to the direct LD-score scope and support gates; storage ownership and lifetimes are in [memory and CPU](design/memory-and-cpu.md#ownership-and-lifetime).

### Reference panels

`_refpanel.plink_source.select_plink_chromosome()` is the one BIM-level selection stage of every PLINK route: `PlinkRefPanel` (`ldscore --plink-prefix` and gene LD-score index construction) records the reference-SNP funnel from the counts it returns, and `R2MatrixBuilder` adds liftover and duplicate-position filtering to the rows it returns.
`RefPanel.prepare_chromosome()` is the boundary between reference panels and LD-score computation: direct `ldscore` chromosome tasks (`compute_direct_chromosome()`) and PLINK gene LD-score index construction call it, and `_ldscore.compute.compute_chromosome()` consumes the returned `PreparedChromosome` without reopening a source.
Preparation keeps only reference SNPs on the annotation SNP grid (`_align_annotations()`) and records each step's removals in the reference-SNP funnel ([reference SNPs on the annotation SNP grid](../reference/ldscore-snp-universes.md#reference-snps-on-the-annotation-snp-grid)).
A panel reads `--ref-panel-snps-file` and `--keep-indivs-file` once, on first use; `_DirectProducer` resolves them on its task panel in the main process (`RefPanel._resolve_run_inputs()`), so the copies sent to chromosome workers carry them.
`R2MatrixRefPanel.check_chromosome()` runs the footer and sidecar-binding checks of one R² matrix chromosome without preparing it; the `ldscore` input preflight (`_ldscore.preflight.inspect_reference_inputs()`) and `R2MatrixRefPanel.check_reference()`, which `LDScoreCalculator` calls, both use it, and preparation opens each chromosome once through `R2MatrixRefPanel.open_chromosome()`.
`_ldscore.compute` picks its kernel from the prepared reader: a `SortedR2BlockReader` streams R² pairs, and a PLINK genotype reader is rewound before its block sums.
`build-r2-matrix` (`_refpanel.builder`) and `query-r2-matrix` (`_refpanel.query`) write and read R² matrices only through `_refpanel.r2_format`; their contracts are in [R² matrices](../reference/formats/r2-matrices.md) and [`query-r2-matrix`](../reference/query-r2-matrix.md).

### LD-score calculation

`ldsc._ldscore.workflow` resolves the inputs of the `ldscore` command and `run_ldscore()` and hands the calculation to `LDScoreCalculator` (`_ldscore.calculator`) with a private `_RunContext`: the resolved `GlobalConfig`, the reference panel built with it, and which steps have already run (output preflight and workflow log, build resolution, R² matrix checks), so `LDScoreCalculator._run_prepared()` repeats none of them.
A direct Python call of `LDScoreCalculator.run()` passes a context that has settled nothing; the calculation then resolves the build against its prepared inputs, checks every R² matrix chromosome, and preflights and opens `diagnostics/ldscore.log` itself.
The parent and every chromosome task build their reference-panel adapter with `_refpanel.adapters.resolved_ref_panel()` from the run's resolved configuration and PLINK chromosome prefixes, so a run resolves its build once; `resolved_ref_panel()` checks no build declaration again, and opening an R² matrix chromosome still compares its declared build with the run's.
`RefPanelLoader.load()`, for direct Python callers, checks the build declarations first and then builds the adapter the same way, as PLINK gene LD-score index construction does.
An option counts as supplied only when `_explicit_cli_options` names it, never by comparison with a parser default, so the registered `GlobalConfig` that a keyword call carries is not a live override in indexed mode.
An indexed run's `run_indexed_ldscore()` relies on the `ldscore` workflow for the output preflight and the open log and says so to `_ldscore.diagnostics.ldscore_logs()` with `main_log_open=True`.

### LD-score producers and query batches

Both LD-score modes run their query batches through `_ldscore.batch_execution.run_query_batches()`, one batch in flight: the direct LD-score producer (`_DirectProducer` in `_ldscore/direct_batches.py`) computes one chromosome from the reference panel per task, and the indexed producer (`_IndexedProducer` in `_ldscore/gene_index/batches.py`) reloads one stored chromosome operator per task.
For direct BED and gene-list queries, `_ldscore.direct_annotation.prepare_direct_annotations()` returns the bundle with the query columns and statuses declared but not projected, plus a `DeferredQueries` record that the calculator hands to the producer; each batch projects its own columns through `_annotation.interval_queries.execution_query_bundle()`.
The producer passes each chromosome's annotation identity-cleanup drops to its first-batch task explicitly, so the bundle carries no LD-score execution state.
A chromosome whose reference preparation leaves no SNP comes back as a skip that carries its reference-SNP funnel (`EmptyReferenceSNPs`); the parent fails the run when the `annotation_snp_grid` step emptied a chromosome (`_refpanel.funnel.annotation_snp_grid_emptied()`) and otherwise logs a WARNING ([empty chromosomes](../reference/ldscore-snp-universes.md#empty-chromosomes)).
A run without annotation files meets such a chromosome earlier, when `prepare_synthetic_base()` reads the panel metadata; the chromosome gets no task, and `_DirectProducer` adds its funnel to the first batch's funnels.

Direct and indexed scoring share every gene-query decision: the coverage check (`with_coverage()` in `_annotation/gene_query_storage.py`), the final query statuses (`_annotation.query_outcomes.finalize_query_measurements()`), the pruning of skipped queries (`select_queries()` in `_ldscore/results.py`), and the stop sequences in `_ldscore/diagnostics.py` (`fail_gene_list_resolution()`, `fail_unusable_queries()`); see [the coverage contract](../reference/results/gene-list-diagnostics.md#chromosome-scope-and-pathway-coverage).
A run that stops before publication describes its diagnostics with `_ldscore.outputs.QueryDiagnostics`, the only input of the writer's diagnostics methods; a gene coordinate catalog that `GeneCatalog.load()` rejects reaches `diagnostics/gene_catalog_issues.tsv` through `QueryDiagnostics.catalog_issues`, and `annotate` and `quantile-h2` write the same table through `AnnotationDirectoryWriter` and `QuantileH2DirectoryWriter`.

### LD-score directory publication

`LDScoreDirectoryWriter.write_batches()` is the only publication path for LD-score directories: direct `ldscore`, indexed scoring and `convert-legacy-ldscores` (as a one-item batch) all publish through it.
It consumes materialized `LDScoreResult` values or private `StagedLDScoreBatch` descriptors behind one interface, and its publication record (`_LDScorePublication.from_first_batch`) keeps the baseline, counts, overlap, provenance and diagnostics but no query values.
The writer performs no cross-batch scientific step: `LDScoreOverlap.concat_query_batches()` in `_ldscore/overlap_table.py` joins the batches' overlap blocks, and `run_query_batches()` applies the gene-list viability gate before the writer publishes, while the writer still writes the query diagnostics of a run the gate stops.
The producer owns its fragments until the writer has consumed them or the workers have shut down; the writer closes the producer before it removes publication staging after a handled failure.
The writer reuses the table and identity validation of `_ldscore.source` and requires each staged fragment's identity columns to equal the shared first-batch baseline rows (`_write_query_fragments()`).
`load_ldscore_from_dir()` reopens the published directory as an `LDScoreSource` for regression, and `quantile-h2` reads only its manifest (`read_ldscore_manifest()`); model-specific validation stays in `_regression.alignment`.
The directory layout is in [artifact metadata](../reference/formats/artifact-metadata.md#ldscore), and the memory design of fragments and batches in [memory and CPU](design/memory-and-cpu.md#direct-ld-score-computation).

### Gene LD-score indexes

`build-gene-ldscore-index` (`_ldscore.gene_index.build`) reuses the direct LD-score path wherever an index must agree with a direct run.
It builds one `PlinkRefPanel` per run, reading the keep list once, and each chromosome worker passes its baseline shard, in canonical (`POS`, `SNP`) order, to `PlinkRefPanel.prepare_chromosome()` as the annotation SNP grid, so the panel's selection, identity cleanup and grid match are those of direct `ldscore`.
Before the workers start, `_validate_baseline_plink_inputs()` compares each chromosome's baseline rows with the rows `select_plink_chromosome()` keeps.
Workers compute baseline LD scores, regression weights, counts, overlaps and annotation types with `_ldscore.compute.compute_chromosome()` before the atom batches, write their chromosome logs through `_ldscore.diagnostics.chromosome_log()`, and take their identity settings from `gene_index.artifact.index_global_config()`.
The coordinator builds root `annotation_types` and `snp_universe_policy` with the helpers direct `ldscore` uses (`_annotation.semantics.combine_annotation_types()`, `_ldscore.outputs.snp_universe_policy_record()`), so indexed metadata matches a direct run.
The build reload-validates its staged and published index with `load_gene_ldscore_index()` (`gene_index.artifact`), the loader indexed scoring uses; `run_indexed_ldscore()` (`gene_index.scoring`) then writes an ordinary LD-score directory through `LDScoreDirectoryWriter` ([gene LD-score index](../reference/gene-ldscore-index.md)).

### Summary statistics

`SumstatsMunger` (`_sumstats.munger`) resolves a run's munging plan with `_sumstats.munge_input.resolve_munge_plan()`, which `--infer-only` and `infer_raw_sumstats()` report from as well, so the previews report the problems the run stops on.
`prepare_munge_input()` then raises the plan's problems and prepares the keep-list; when the plan selects HM3 quick liftover it loads the packaged HM3 map once and hands the same frame, through `ResolvedMungeInput.hm3_map`, to the kernel's liftover (`_sumstats.liftover.apply_sumstats_liftover()`).
`_sumstats.munge_kernel.munge_sumstats()` receives only resolved settings and returns a `MungeResult`; `_sumstats.outputs.curated_sumstats_frame()` turns it into the one frame that the returned `SumstatsTable` and both data files hold, and `SumstatsDirectoryWriter` writes it.
Regression and every other reader import the munged artifact from `_sumstats.table` (`SumstatsTable`, `load_sumstats()`, the footer codec), never from the munging workflow; the command contract is in [`munge-sumstats`](../reference/munge-sumstats.md) and the stage order in [data flow](data-flow.md#5-munge-sumstats-raw-summary-statistics-to-munged-summary-statistics).

### Regression

The regression modules import only the ones listed before them: `results` and the writers in `outputs`, `alignment` (trait preparation and dataset assembly), `fitting` (fit policy, fits, result rows and metadata), `partitioned` and `rg` (the `partitioned-h2` batch and the `rg` pair run), and the two entry layers, `commands` (the three commands and their Python forms) and `runner` (the public `RegressionRunner`), which do not import each other.
Both entry layers and the `rg` pair run go through one preparation spine in `alignment` (`_prepare_h2_dataset()`, `_h2_dataset()`, `_prepare_rg_dataset()`) and the fits in `fitting` (`_fit_h2_dataset()`, `_fit_rg_dataset()`), each called with one explicit `RegressionConfig`; the entry point, not these modules, logs the `GlobalConfig` banner.
Each trait's summary statistics are prepared once per run and shared by every model; the dataset rules (missing statistics, allele alignment, identity keys, collinear columns) are in [regression options](../reference/regression-options.md).
`partitioned.estimate_partitioned_h2_batch()` owns input preparation, batch reads and the publication of its result directory, and passes the module-level `fitting.fit_partitioned_query()` to the query workers in `_regression.partitioned_parallel`, which import neither `partitioned` nor `fitting`; `h2` and `rg` publish from their command bodies.
The Python forms `run_h2()`, `run_partitioned_h2()`, `run_rg()` and `run_quantile_h2()` pass their keyword arguments and the registered `GlobalConfig` to `_runtime.workflow_kwargs.run_from_kwargs()`, which checks them against the command's own parser and calls the same `run_*_from_args()` function as the command line.
`run_partitioned_h2()` checks `workers` and `query_batch_size` before that call, because their argparse types check them only on the command line.
`quantile-h2` (`_regression.quantile_h2`) builds no `GlobalConfig`: it takes the SNP identity and genome build from the fitted model's LD-score directory, which it reads only through `_ldscore.source.read_ldscore_manifest()`, and rebuilds the fitted annotations with `_annotation.loading.build_annotation_shards()`.
Its fitted-model reader and the `partitioned-h2` writer share the file names and locators in `_regression/outputs.py` (`fitted_model_files()`, `published_query_artifacts()`); the `global_metadata` reference and the copying rules are in [partitioned-h2 results](../reference/results/partitioned-h2.md).

### Workers

Direct scoring, indexed scoring and gene LD-score index construction share `_runtime.parallelism.chromosome_pool()`, which turns a crashed worker into `LDSCInternalError`.
`partitioned-h2` query workers live in `_regression.partitioned_parallel` and share the worker-count validation of `_runtime.parallelism`, and standalone numerical entry points run in one controlled child through `_runtime.serial_calculation`.
Worker budgets, the one-BLAS-thread policy and failure handling are in [parallel execution](../reference/parallel-execution.md), and the memory and lifetime designs in [memory and CPU](design/memory-and-cpu.md).

## Cross-Cutting Concerns

- **Input token language**: public inputs accept exact paths, globs and `@` chromosome suites; PLINK inputs also accept plain suite stems through `_runtime.paths.resolve_plink_prefix_group()`, and output paths are literal destinations ([input paths](../reference/input-paths.md)).
- **Column names and SNP identity**: `_core/columns.py` is the single source of truth for input header aliases, internal artifact headers, SNP identifier modes and genome-build aliases, and `_core.snp_identity` for identity keys and allele-aware cleanup.
  The SNP identifier modes, the duplicate and restriction rules and the meaning of `A1` in each file type are in [column schema](../reference/formats/column-schema.md#snp-identifier-modes); why identity cleanup follows the SNP restriction is in [SNP identity and genome build](../reference/snp-identity-and-genome-build.md).
- **Genome builds**: `_genome/build_inference.py` infers hg19/hg38 and the coordinate basis of `chr_pos` tables, and the build policy decides a workflow's build once ([genome build](../reference/snp-identity-and-genome-build.md#genome-build-policy)).
- **Logging**: workflow and kernel records reach one per-run file through the `LDSC` logger; the CLI console shows errors only, and `stdout` carries deliverable command output ([logs and outputs](../reference/logs-and-outputs.md)).
- **Dependencies**: PLINK-backed LD computation needs the `plink` extra and chain-file liftover the `liftover` extra, while BED and gene-list interval projection needs only the base installation ([installation](../user-guide/installation.md#python-versions-and-extras)).
- **Chromosome ordering**: workflows validate chromosome-sharded inputs and reassemble them in stable genomic order.

## Architectural Invariants

These are the rules a change must keep; [code structure](code-structure.md) links here for them.

- Import from `ldsc` for supported public use; every package whose name starts with an underscore is private, even when `ldsc` re-exports a name it defines.
- `ldsc.cli` parses arguments and dispatches; numerical and plotting logic belongs to a workflow or kernel.
- Numerical projection consumes resolved, prepared inputs: user-token resolution belongs to workflows, shared helpers and the runtime adapters in `_refpanel.adapters`, and computation modules such as `_ldscore.compute`, `_refpanel.r2_format` and `_regression.estimators` take concrete files, own no public output layout and carry no workflow diagnostics.
- `_core/columns.py` owns alias families and internal artifact header strictness.
- LD-score computation remains chromosome-wise, global query filtering uses the complete chromosome summaries, and publication follows the success of every required chromosome.
- Regression is file-driven: it rebuilds its state from written artifacts and never recomputes LD scores.
  Plotting and scale conversion read saved results and never refit, and Matplotlib loads only on a figure-producing path.
- Public file contracts stay stable unless a change is intentional and coordinated: `.annot(.gz)`, munged summary statistics (Parquet with its identity footer, optionally also `.sumstats.gz`), R² matrices, LD-score directories (flat Parquet files with one row group per chromosome, described by root `metadata.json`) and regression result directories.
  Emitted Parquet schemas are part of that contract.
- Headline results sit at the result root; logs, dropped-SNP tables and diagnostic metadata live under `diagnostics/`.
- New internal or intermediate artifacts are Parquet and new science-facing result tables are TSV; compatibility text outputs are called out as legacy or sidecar exceptions.
- Legacy reference and weight `.l2.ldscore(.gz)`, `.l2.M` and `.l2.M_5_50` families enter only through `convert-legacy-ldscores` ([compatibility boundary](../reference/formats/legacy-ldscores.md#compatibility-boundary)); the `.w.l2.ldscore(.gz)` suffix is unsupported, and no kernel emits and no regression reader reads legacy LD-score files.
- Package-built R² matrices are the only public R² Parquet format for `ldscore --r2-matrix-dir` and `query-r2-matrix`: both depend on the paired `chr*_r2.parquet` and `chr*_meta.tsv.gz` files and their sidecar identity binding.
- `LDScoreDirectoryWriter` fixes the public LD-score output layout, and `PartitionedH2DirectoryWriter` the compact and full `partitioned-h2` result schemas.
  Public output customization chooses the directory name and an explicit overwrite policy, not per-run filename prefixes.
- Public annotation reads and final LD-score Parquet values are float32; LD accumulation, annotation sums and overlap products use float64.
- Query annotations are valid only with explicit baseline annotations; the LD-score workflow layer builds the synthetic all-ones `base` only for ordinary unpartitioned LD scores and downstream `h2`/`rg`, not for `partitioned-h2`.
- BED and gene-list queries use one partial-success contract: only usable query columns enter scientific artifacts, a skipped query never becomes a placeholder column, and skipped inputs stay visible in `diagnostics/query_annotation_status.tsv`.
- BED projection uses input intervals as provided unless `padding_bp` / `--padding-bp` is set.
  Padding expands both interval ends in base pairs before projection and clips starts at zero; it should not be applied again to BED files already expanded upstream.
  `ldscore` accepts this option only when it is performing live BED or gene-list projection; prebuilt-annotation, no-query and indexed modes reject it.
- The `ldsc annotate` entry points stay in `_annotation.builder` (`build_parser()`, `run_annotate_from_args()`, `main()`), which the CLI dispatches to; `_annotation.parsing` owns no parser functions, result objects or public workflow aliases.
- Every workflow that writes fixed artifacts precomputes its output paths, including its log, and calls the shared output preflight before the first write.
  It creates a missing output directory, reuses an existing one, refuses existing owned artifacts unless `--overwrite` or `overwrite=True` is given, and never deletes unrelated files.
  Coherent output families (`munge-sumstats`, `build-r2-matrix`, `ldscore`, `partitioned-h2`, `rg`, `annotate`) reject any owned sibling without overwrite and remove stale owned siblings after a successful overwrite; a `build-r2-matrix` run with a concrete chromosome prefix owns only that chromosome's files, and an `@` suite run owns the whole R² matrix ([output-family preflight](../reference/logs-and-outputs.md#output-family-preflight)).
- Logs are diagnostic files: result `output_paths` mappings and metadata that downstream code reads as scientific data exclude them.
  The gene LD-score index builder keeps its lock, open log and staged index in the index's hidden `.build/` directory, publishes root `metadata.json` last, and moves only the closed successful log into the published diagnostics.
- Gene LD-score index construction and assembly preserve Boolean interval union, float64 unbiased-R² accumulation (including negative values), one diagonal contribution, and the split between the broad reference SNPs and the filtered regression SNPs.
- Indexed scoring uses only the explicitly named complete gene LD-score index; it never discovers an index or falls back to live computation.
- Extend a shared helper or the workflow-owned writer rather than duplicating local parsing or writing logic.
