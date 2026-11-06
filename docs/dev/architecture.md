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

- **Build query annotations**: project BED or resolved gene intervals onto an annotation SNP grid; standalone `annotate` accepts mutually exclusive BIM or baseline-annotation sources. Entry points: `ldsc annotate`, `ldsc ldscore`, `ldsc.AnnotationBuilder`
- **Build R² matrices**: convert PLINK reference panels into R² matrices.
  Entry points: `ldsc build-r2-matrix`, `ldsc.R2MatrixBuilder`
- **Query R² matrices**: look up unbiased R², `sign_r`, and signed Pearson `r` for SNP pairs in package-built R² matrices.
  Entry points: `ldsc query-r2-matrix`, `ldsc.R2Matrix`, `ldsc.query_r2_matrix()`
- **Build exact gene LD-score indexes**: precompute PLINK-backed disjoint-atom operators for one complete immutable baseline/reference/gene configuration. Entry point: `ldsc build-gene-ldscore-index`
- **Compute LD scores**: align annotations to a live reference panel or explicitly assemble gene-list columns from a validated complete index, then emit the same canonical artifacts. Entry points: `ldsc ldscore`, `ldsc.run_ldscore()`, `ldsc.LDScoreCalculator`
- **Munge raw summary statistics**: normalize raw GWAS tables into Parquet-first munged summary statistics, with optional legacy `.sumstats.gz` output.
  Entry points: `ldsc munge-sumstats`, `ldsc.SumstatsMunger`
- **Run LD score regression**: consume munged summary statistics and LD-score directories to estimate SNP heritability, partitioned heritability, or genetic correlation.
  Entry points: `ldsc h2`, `ldsc partitioned-h2`, `ldsc rg`, `ldsc.run_h2()`, `ldsc.run_partitioned_h2()`, `ldsc.run_rg()`, `ldsc.RegressionRunner`
- **Explore canonical results**: dispatch one approved figure from `h2`, `partitioned-h2`, `quantile-h2`, or `rg` result metadata without refitting.
  Entry points: `ldsc plot`, `ldsc.plot_result()`
- **Convert saved heritability scale**: convert observed-scale SNP heritability at one population prevalence or over a sensitivity grid.
  Entry points: `ldsc convert-h2-scale`, `ldsc.convert_h2_scale()`
- **Project continuous-target quantiles**: verify one fitted partitioned model and resupplied SNP annotations, then calculate joint-model quantile heritability and standardized coefficients without refitting. Entry points: `ldsc quantile-h2`, `ldsc.run_quantile_h2()`
- **Audit workflow runs**: artifact-writing workflow wrappers create deterministic
  per-run logs around their owned work. The gene LD-score index builder writes its log
  to the hidden `.build/` directory inside the index directory, then moves a closed
  successful log into the published `diagnostics/` and removes `.build/`. Logs are diagnostic files and are not
  included in result `output_paths`.

## Layer Structure

- **CLI Layer**: public command dispatch in `ldsc.cli`
- **Workflow And Preprocessing Layer**: public services in `ldsc._annotation.builder`, `ldsc._refpanel.builder`, `ldsc._refpanel.query`, `ldsc._ldscore.workflow`, `ldsc._ldscore.calculator`, `ldsc._sumstats.munger`, the regression workflow modules (`ldsc._regression.commands`, the `RegressionRunner` facade in `ldsc._regression.runner`, and `ldsc._regression.alignment`, `ldsc._regression.fitting`, `ldsc._regression.partitioned` and `ldsc._regression.rg`), and `ldsc._regression.quantile_h2`, plus shared normalization in `ldsc.config`, `ldsc._runtime.paths`, `ldsc._core.columns`, `ldsc._core.chromosomes`, `ldsc._genome.build_inference`, `ldsc._annotation.semantics`, and the internal `ldsc._annotation.gene_lists`
- **Computation**: private file-format and numerical code inside each domain package, such as `ldsc._ldscore.compute`, `ldsc._refpanel.r2_format`, `ldsc._refpanel.plink` and `ldsc._regression.estimators`
- **Shared Artifact Access**: LD-score directory opening, schema/identity validation, and selected-query reads in `ldsc._ldscore.source`, consumed by generation and regression
- **Output And Derived-Result Layer**: result writing in each domain's `outputs` module (`ldsc._annotation.outputs`, `ldsc._refpanel.outputs`, `ldsc._ldscore.outputs`, `ldsc._regression.outputs`) on the shared `ldsc._runtime.artifacts` base; metadata-driven plotting in
  `ldsc._plotting.dispatch`; post-fit scale conversion in `ldsc._regression.h2_scale`

### Ownership by layer and workflow

The matrix answers which layer owns an object and which objects take part in a workflow.
The public import boundary is `ldsc`, the command boundary is `ldsc.cli`, and every package whose name starts with an underscore is internal.
The packages form one-directional tiers, `_core` < `_runtime` < `config` < `_genome` < (`_annotation`, `_refpanel`, `_sumstats`) < `_ldscore` < `_regression` < `_plotting` < the root shell, enforced by `tests/package/test_layering.py`, which also fails on import cycles and on unlisted imports inside functions; the [code structure](code-structure.md#dependency-direction) gives the rules.

| Layer | Shared CLI, config, path, and schema contracts | Annotation building | Reference-panel building and loading | LD-score calculation | Summary-statistics munging | Regression | Output and persistence |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Public API and command surface | `ldsc.__init__` re-exports stable objects; `ldsc.__main__.main()` runs `python -m ldsc`; `ldsc.cli.main()` dispatches the `ldsc` subcommands listed in `ldsc.cli._COMMANDS`, and `ldsc.cli._command_parser()` returns the parser each one uses. Top-level imports include build inference, R² query, lazy `plot_result()`/`PlotArtifact`, and lazy `convert_h2_scale()`/`H2ScaleConversionArtifact`; the CLI also exposes the read-only `infer-genome-build PATH` helper through private `ldsc._genome.build_inspection`, with no new Python export. | CLI: `ldsc annotate`; imports: `AnnotationBuilder`, `AnnotationBundle`, `AnnotationBuildConfig`, `run_annotate()`. | CLI: `ldsc build-r2-matrix`, `ldsc query-r2-matrix`; imports: `R2MatrixBuilder`, `R2MatrixBuildResult`, `run_build_r2_matrix()`, `R2Matrix`, `query_r2_matrix()`, `RefPanel`, `RefPanelLoader`, `RefPanelConfig`, `PlinkRefPanel`, `R2MatrixRefPanel`. | CLI: `ldsc ldscore`; imports: `LDScoreCalculator`, `LDScoreResult`, `LDScoreSource`, `ChromLDScoreResult`, `run_ldscore()`. | CLI: `ldsc munge-sumstats`; imports: `MungeConfig`, `SumstatsTable`, `SumstatsMunger`, `MungeRunSummary`, `load_sumstats()`. | CLI: `ldsc h2`, `ldsc partitioned-h2`, `ldsc quantile-h2`, `ldsc rg`, `ldsc convert-h2-scale`, `ldsc plot`; imports: `run_h2()`, `run_partitioned_h2()`, `run_quantile_h2()`, `run_rg()`, regression types, and the lazy conversion and plotting APIs. | Imports: canonical writer configs/classes. `ldsc._plotting.dispatch` owns fixed nested plot families; `ldsc._regression.h2_scale` owns scale-conversion families; `ldsc._runtime.workflow_logging` owns audit logs and failure markers. |
| Configuration objects | `GlobalConfig` stores the shared SNP identifier mode, genome build, and log level; `get_global_config()`, `set_global_config()`, `reset_global_config()`, and `validate_config_compatibility()` manage and compare shared assumptions. | `AnnotationBuildConfig` holds baseline, query, BED and gene-list sources, the gene coordinate catalog, the control gene list, the gene-list resolution policy, gene exclusion, and padding. | `R2MatrixBuildConfig` configures PLINK-to-R²-matrix builds, source-build inference, optional coordinate liftover, explicit SNP/sample restrictions, and optional `r2_min`; `RefPanelConfig` configures runtime PLINK or R²-matrix loading, including backend, source paths, explicit SNP/sample restrictions, retained-panel MAF filtering, and optional genetic maps; its construction rejects a missing source path and any field of the other backend (`keep_indivs_file`, `plink_prefix` and genetic maps apply only to PLINK, `r2_matrix_dir` only to an R² matrix). Runtime genome-build assumptions come from `GlobalConfig`. | `LDScoreConfig` configures one LD window, optional explicit or bundled-HM3-default regression SNP restriction, `snp_batch_size`, `query_batch_size`, the common-MAF count threshold, the whole-chromosome window override, reference-metadata export, and `workers`. | `MungeConfig` configures `raw_sumstats_file`, trait name, raw `sumstats_format`, INFO-list columns, column hints, QC thresholds, sample-size fallbacks, output directory, optional explicit or HM3 summary-statistics SNP restriction, optional `chr_pos` liftover target/method, and the increasing-allele assumption. | `RegressionConfig` configures the block-jackknife block count, count-vector choice, intercept handling, two-step cutoff, chi-square filtering, and prevalence settings. | `LDScoreOutputConfig`, `H2OutputConfig`, `PartitionedH2OutputConfig`, `QuantileH2OutputConfig`, and `RgOutputConfig` configure fixed result directories. |
| Shared normalization and inference | `_runtime/paths.py`, schema/build inference, SNP identity, liftover, and `_runtime/workflow_logging.py` own shared normalization, output preflight, lifecycle logging, and failure-marker I/O. Core CLI workflows require `--output-dir`; `plot` and `convert-h2-scale` derive fixed nested destinations. | Annotation inputs and generated query shards use shared resolution and coherent artifact-family preflight. | Builders/loaders resolve PLINK, R²-matrix, map, restriction, and liftover inputs; build outputs are preflighted by chromosome or suite. | LD-score resolution validates annotations, panels, identity/build metadata, and its complete output family before writing. | Munging infers raw formats and columns, applies restrictions/liftover, and preflights its munged outputs; infer-only reports without writing. | Regression validates compatibility before merging and always preflights the required output directory. `partitioned-h2` query runs own the per-query tree by default; baseline-only complete-model artifacts stay at the root. `rg` pair details remain opt-in. Plot dispatch checks only plotting-required metadata fields. | Output paths are literal or fixed derived destinations, overwrite is explicit, stale owned siblings are cleaned, and workflow logs/markers are excluded from scientific result mappings. |
| Workflow services | `ldsc.cli` dispatches only; public workflow services resolve user inputs, call kernels, and own output/log orchestration. | `AnnotationBuilder.run(..., output_dir=...)` returns an owned chromosome handle. `run_annotate()` publishes BED/gene query shards incrementally, releases construction scratch, and returns an `AnnotationResult` of persistent output/source references. | `R2MatrixBuilder` writes panel artifacts and diagnostics; `R2Matrix` serves required-output pair queries; `RefPanelLoader` chooses runtime backends. | `LDScoreCalculator` computes in memory when no output is requested, but metadata export requires output configuration; `run_ldscore()` requires and writes a canonical output directory and checks its keyword arguments against the `ldscore` parser in private `_runtime.workflow_kwargs`, which also fills the SNP identifier mode, genome build and log level from the registered `GlobalConfig`. | `run_munge_sumstats_from_args()` maps CLI args to `MungeConfig`; `SumstatsMunger` owns munged outputs and diagnostics. | `RegressionRunner` retains in-memory numerical APIs and the current batch return type. Parsed CLI workflows require output and write command logs. `run_h2()`, `run_partitioned_h2()`, `run_rg()`, and `run_quantile_h2()` check keyword arguments against the command parser in private `_runtime.workflow_kwargs`, then run the same parsed workflows. Query batches stage completed fit details and return persistent file references. | Directory writers own canonical result layouts; low-level Python output configs retain their explicit controls. |
| Cross-layer data and result objects | `ConfigMismatchError` marks incompatible shared assumptions. `ColumnSpec` describes a canonical field plus accepted aliases. | `AnnotationBundle` carries chromosome descriptors, ordered column names, source provenance, configuration, and an explicit private-storage owner. Selected reads return caller-owned arrays. | `R2MatrixBuildResult` summarizes emitted R²-matrix files. `R2MatrixChromosome` opens one chromosome of an R² matrix for `ldscore` and `query-r2-matrix` and keeps its validated sidecar; `R2Matrix` caches per-chromosome identity-key and Parquet state for repeated pair queries. `RefPanel` defines `available_chromosomes()`, `load_metadata()`, `build_reader()`, `check_reference()`, `reference_metadata_sources()`, and `prepare_chromosome()`; `build_reader()` has no caller inside the package and stays for direct Python use. | `PreparedChromosome` binds aligned metadata/annotation rows and LD windows to one owned reader. `ChromLDScoreResult` stores one chromosome's split baseline/query tables, SNP count vectors, and provenance. `LDScoreResult` holds one active genome-wide query batch. Both live in `_ldscore/results.py` with the private file-backed `DirectChromosomeFragment` and `StagedLDScoreBatch`, so regression, the legacy converter and the writer import them without the LD-score orchestrator. Writing workflows and `load_ldscore_from_dir()` in `_ldscore/source.py` return `LDScoreSource` with shared baseline/metadata and explicit selected-query reads across batch files; completed query tables are released. Prepared small inputs can use `AnnotationBundle.from_frames` and return a single in-memory result without filesystem writes. | `SumstatsTable` carries validated `SNP`, `CHR`, `POS`, `Z`, `N`, optional allele columns, provenance, and the current config snapshot. `MungeRunSummary` records row counts, inferred columns, sample-size rule, and output paths. | `RegressionDataset` carries the merged regression table of one trait (`h2`) or one harmonized trait pair (`rg`), retained LD columns, regression LD-score column, selected count vectors, dropped columns, trait names, and provenance. `RegressionFit` carries a fitted kernel estimator with the SNP rows, block count, chi-square cap, two-step cutoff, and diagnostic bins of its final fit. `H2Result` carries the `h2` summary, the paths written, and legacy projection drops. `RgResultFamily` carries the concise genetic-correlation table, full genetic-correlation diagnostics, per-trait SNP heritability, per-pair metadata, legacy projection drops, per-trait non-finite statistic counts, and the paths written; genetic-correlation p-values are nominal and uncorrected. | `ArtifactFamily` carries selected paths and the full owned scope from each directory writer; workflows use that same declaration for collision checks. `LDScoreOutputConfig`, `PartitionedH2OutputConfig`, and `RgOutputConfig` describe result-directory options. |
| Internal compute kernels and adapters | `_core.plink_files` contains only retained PLINK and identity-list reader primitives; obsolete legacy summary-statistics, LD-score, count, annotation, and frequency readers were removed. `_core.identifiers` parses SNP restriction files and opens package-written text tables past their leading `# ldsc:` identity lines (`open_identity_table`, used by the R² matrix sidecar reader, by the genome-build check of an R² matrix sidecar's declaration and by `quantile-h2` reference metadata, which resolves its columns through `_core.columns`; the build check's `_text_options` still reads user-supplied text tables, sniffing their delimiter past any leading comment) and delegates exact mode semantics to `_core.snp_identity`, which builds base/effective keys, performs allele-aware cleanup, collapses restriction keys, validates minimal identity provenance, and coerces every dropped-SNP table to one schema (`coerce_identity_drop_frame`). `_genome.liftover` translates positions through chain files, detects duplicate coordinate groups, and builds readable drop reports for munging and R² matrix building; summary-statistics liftover, including HM3 quick liftover, is in `_sumstats.liftover`; `_core.regions` loads packaged/custom exclusion intervals; the public packaged HM3 loader lives in `ldsc._genome.hm3`. These modules are internal even when public workflows rely on them. | `_annotation.parsing` normalizes annotation chunks and validates their numeric values. Live BED parsing uses `_core.regions.iter_bed_rows`; `_annotation.interval_queries.build_query_shards` projects BED/gene intervals through `_annotation.interval_projection.ChromosomeProjector`. `_annotation.builder` exposes the API/parser, `_annotation.annotate` orchestrates standalone output, and `_annotation.bundle` owns the lazy dataset. | `_refpanel.build_kernel` holds only `build-r2-matrix` code: it computes pairwise R² rows in the exact LD window and re-orients genotype columns; `build-r2-matrix` liftover is in `_refpanel.builder`. `_refpanel.genetic_map` loads, combines and interpolates genetic maps for `ldscore`, gene-index construction and `build-r2-matrix`, with messages that name the running command; each command keeps its own policy for which map it uses. `_refpanel.plink_source.select_plink_chromosome()` is the one BIM-level selection stage of every PLINK route: it reads a chromosome's `.bim` rows, applies the SNP restriction, runs identity cleanup on the rows it keeps, and maps the keep list to `.fam` rows; the same module matches reference-panel SNP restrictions (`build_restriction_mask()`) for the PLINK stage and the R² matrix adapter, and builds PLINK SNP metadata with the minor allele as `A1` (`build_plink_metadata_frame()`). `_refpanel.r2_format` owns the R² matrix format: its file names, the writers of the 4-column R² Parquet files and mandatory sidecars, `R2MatrixChromosome` (one chromosome's footer contract, strict sidecar reading and sidecar binding, parsed once per process), pair streaming through `SortedR2BlockReader`, and low-level pair lookup. `_refpanel.windows` owns LD-window geometry, `_refpanel.prepared` defines `PreparedChromosome`, and `_refpanel.funnel` counts the SNPs each preparation step removes. `_refpanel.adapters.PlinkRefPanel` and `_refpanel.adapters.R2MatrixRefPanel` adapt runtime sources to the `RefPanel` interface and own chromosome preparation, including filtering, alignment, R²-matrix contract and window policy, and reader lifetime. | `_ldscore.compute` owns the LD-score math over a `PreparedChromosome`: R²-matrix pair streaming (`ld_score_streaming_from_r2_reader`), PLINK genotype block sums (`ld_scores_from_genotypes`), the float64 `ProjectionAccumulator`, chromosome LD-score projection, the regression-SNP mask, and count calculation; it has no legacy artifact emitters. | `_sumstats.munge_kernel` preserves the legacy LDSC munging behavior: column parsing, filters, p-to-z conversion, sample-size handling, optional post-filter coordinate liftover, and identity cleanup; it returns `MungeResult` with tables, accounting and provenance and owns no output files. | `_regression.estimators.Hsq`, `_regression.estimators.RG`, `_regression.estimators.LD_Score_Regression`, `_regression.jackknife.*`, and `_regression.irwls.IRWLS` implement the numerical estimators and supporting optimization/jackknife routines. | Computation modules do not own public output layout decisions; R²-matrix Parquet and sidecar emission belongs to `_refpanel.r2_format`. |
| On-disk artifacts owned by the layer | Public docs, tests, and tutorials define user-visible contracts; materializing CLI commands take a required output directory while derived commands take a source result root; `infer-genome-build` accepts one exact input path and writes only to the terminal. | Standalone runs use mutually exclusive BIM/annotation grids and write root `query.<chrom>.annot.gz`, plus aligned `baseline.<chrom>.annot.gz` for annotation-grid input, diagnostics, and a log. | R² matrix builds write R² Parquet files, sidecars, drop audits, and logs. `query-r2-matrix` always writes `query_r2_matrix.tsv` plus diagnostic metadata and log under its output directory. | LD-score runs write canonical metadata, baseline/query/overlap parquet as applicable, and `diagnostics/ldscore.log`. | Munging writes self-describing parquet, optional legacy TSV, drop audit, and log; infer-only writes nothing. | `h2` also writes exact regression-bin diagnostics; `partitioned-h2`, `quantile-h2`, and `rg` write stable result families. `plot` writes one PNG below `plots/`; heritability scale conversion writes a table and optional sensitivity PNG below `postprocessing/`. | Writers own core layouts; plotting/conversion own derived layouts; workflow logging owns logs and non-scientific `RUN_FAILED` markers. |
| Tests and verification | `tests/package/test_package_layout.py`, `tests/core/test_config_identifiers.py`, `tests/package/test_global_config_registry.py`, `tests/runtime/test_path_resolution.py`, `tests/core/test_column_inference.py`, `tests/core/test_chromosome_inference.py`, `tests/genome/test_genome_build_inference.py`, `tests/core/test_regions.py`, and `tests/runtime/test_logging_refactor.py`. | `tests/annotation/test_annotation.py`, `tests/annotation/test_annotate_grid.py`, `tests/annotation/test_annotate_streaming.py`, and `tests/core/test_regions.py`. | `tests/refpanel/test_ref_panel.py`, `tests/refpanel/test_r2_matrix_builder.py`, `tests/refpanel/test_r2_matrix_query.py`. | `tests/ldscore/test_ldscore_workflow.py`, `tests/refpanel/test_plink_io.py`, `tests/runtime/test_output.py`. | `tests/sumstats/test_sumstats_munger.py`. | `tests/regression/test_alignment.py`, `tests/regression/test_fitting.py`, `tests/regression/test_partitioned.py`, `tests/regression/test_rg.py`, `tests/regression/test_regression_commands.py`, `tests/regression/test_irwls.py`, `tests/regression/test_jackknife.py`. | `tests/runtime/test_output.py` plus workflow tests that assert written artifacts can be consumed downstream. |

### Ownership boundaries

The continuous-annotation extension keeps semantics and post-fit work in the public layers: `_annotation/semantics.py` owns advisory classification; routine annotation fingerprints are not computed, `_regression/fitting.py` (each query model) and `_regression/partitioned.py` (the baseline-only model) persist one delete matrix per fitted model, `_regression/quantile_h2.py` owns source reconstruction and projection, and `_regression/outputs.py` owns the result family. The regression kernel is unchanged.

Explicit result exploration follows the same direction.
`_plotting.dispatch` validates canonical result metadata and owns derived figure output, while `_plotting.builders` alone imports Matplotlib.
`_regression.h2_scale` loads a canonical `h2` result and reuses the regression kernel's conversion factor.
Neither path is called automatically by a core workflow.

Gene-list support follows the same boundary.
`ldsc._annotation.gene_lists` validates the required or embedded gene coordinate catalog and returns compact batch resolution/audit records; `ldsc._annotation.query_outcomes` defines internal BED/gene status records; `AnnotationBuilder` projects intervals; LD-score query-batch finalization (`_ldscore.batch_execution.finalize_batch()`) prunes zero-hit and zero-variance queries; and `LDScoreDirectoryWriter` owns the status and unresolved-gene diagnostics.
No gene identity enters the LD-score computation in `_ldscore.compute`.

Legacy compatibility follows the same boundary.
`ldsc._ldscore.legacy_converter` alone reads selected legacy LD-score/count/annotation/frequency suites and emits a canonical `LDScoreResult` directory; it is exported as `LegacyLDScoreConverter` and `convert_legacy_ldscores()`.
Regression never reads those fragments.
`ldsc._regression.alignment` automatically projects legacy summary-statistics text onto the regression SNPs of the already-loaded LD-score directory, once per trait per run, and owns the resulting drop audit.
Neither compatibility path adds a legacy reader or writer to the computation modules.

PLINK discovery lives in `_runtime.paths.discover_plink_members()`, which finds the selected BED/BIM/FAM trios and reports missing members; the declared-input gate calls it, so it stays with path resolution.
PLINK content checks and chromosome assignment live in `_refpanel.plink_source.inspect_plink_inputs()`, beside the PLINK chromosome selection stage.
It uses the low-level BIM/FAM readers for structural inspection (positions, the BED header and size), assigns chromosomes from `.bim` contents, returns combined issue records in a `_runtime.paths.PlinkInputResolution`, and supplies concrete chromosome-prefix mappings to direct LD-score preflight, gene LD-score index construction, and R²-matrix construction.
`PlinkRefPanel` consumes a supplied mapping or resolves it lazily once for direct Python use.
The kernel no longer owns a PLINK path-token wrapper.

## File Tree

```text
ldsc_py3_restructured/
├── docs/
│   ├── user-guide/              # tutorial, concepts, and workflow guides
│   ├── reference/               # command, format, result, and method contracts
│   ├── dev/
│   │   ├── architecture.md      # package-level architecture
│   │   ├── code-structure.md    # contributor entry and module map
│   │   └── data-flow.md         # workflow-level file streams
│   └── assets/
├── src/ldsc/
│   ├── __init__.py  __main__.py  cli.py  config.py  errors.py   # root shell
│   ├── _core/               # chromosomes, columns, coordinates, SNP identity, regions, PLINK companion files
│   ├── _runtime/            # paths, input preflight, workflow logging, progress, workers, keyword checks, artifact writers
│   ├── _genome/             # genome-build inference and policy, build inspection, liftover, HM3 map
│   ├── _annotation/         # annotation files, gene lists, interval projection, overlap, `ldsc annotate`
│   ├── _refpanel/           # PLINK reader, R² matrix format, LD windows, prepared-chromosome contract, reference-SNP funnel, adapters, build and query
│   ├── _sumstats/           # `munge-sumstats` and the munged-table loader
│   ├── _ldscore/            # LD-score math and workflow, directory reader and writer, legacy conversion
│   │   └── gene_index/      # gene LD-score index build, format, scoring, atoms and batches
│   ├── _regression/         # estimators, h2/partitioned-h2/rg, liability scale, quantile-h2, result writers
│   └── _plotting/           # metadata dispatcher and private Matplotlib builders
├── tests/                   # one folder per domain package plus package/ for the root shell; shared helpers and fixtures/ at the top
└── tools/                   # data-curation and figure scripts
```

## Code Map

### `ldsc.cli`

This is the only command-line entry point.
It parses subcommands, preserves the user-facing CLI contract, and hands control to the public workflow modules.
The `_COMMANDS` table lists each subcommand once, with its one-line help, its workflow module, and either the module's own `build_parser()` and `main(argv)` or the `add_*_arguments` and `run_*_from_args` functions for a parser the CLI builds.
A workflow module is imported only when its command is parsed or run, so `ldsc --help` does not import the SciPy-heavy workflows.
Architecture invariant: CLI code should dispatch only; it should not contain numerical logic.

### `ldsc.config`, `ldsc._runtime.paths`, `ldsc._core.columns`, `ldsc._core.chromosomes`, `ldsc._genome.build_inference`

These modules define the package-wide contracts that every workflow shares.
They normalize path tokens, chromosome ordering, genome-build aliases, exact SNP identifier modes, input header aliases, and `chr_pos`-family genome-build inference before the kernel is called.
`ldsc._runtime.paths` also owns output-directory creation, fixed-artifact collision preflight, coherent artifact-family preflight, and stale owned-sibling cleanup helpers.
Its `output_directory_path()` turns a destination into a path and rejects an empty or blank one, which would otherwise mean the working directory; `ensure_output_directory()`, the failure-marker scope, `_runtime.scratch.ScratchWorkspace`, and the `plots/` and liability-scale roots use it, and `cli_help.output_directory_argument` rejects a blank destination option while parsing.
`ldsc._genome.build_inference` is public for Python callers through the top-level `ldsc` exports, while the CLI exposes it through workflow build options and the read-only `infer-genome-build` helper.
Architecture invariant: numerical functions consume resolved inputs.
Shared normalization stays here or in the workflow layer; the runtime reference-panel adapters currently live at `ldsc._refpanel.adapters` and resolve panel-specific inputs before numerical projection.

### `ldsc._genome.build_policy`

`GenomeBuildCheck` collects declared, inferred, unresolved, and unchecked evidence before resolving a workflow coordinate build.
It uses shared path expansion and column aliases, trusts concrete artifact declarations, and bounds coordinate inference per file/chromosome.
All relevant declarations are compared after alias normalization, including labels produced by earlier inference.
Concrete user builds bypass fresh inference while requiring declarations to agree; unlabeled inputs inherit the explicit assertion.
Auto mode reconciles declarations and necessary per-input/per-chromosome inference, failing on conflicts or unresolved required builds.
Build-specific coordinate columns and R² matrix build children represent selectable alternatives, and explicit liftover source/target stages are validated separately.
Workflows and reference adapters call this helper before filtering or numerical work; numerical readers receive a concrete build and do not repeat the decision.
Restriction column selection in `_core/columns.py` prefers the matching named position column, then generic `POS`.
`_genome/build_inference.py` retains the HM3 confidence and coordinate-basis calculations.
See the [build policy](../reference/snp-identity-and-genome-build.md#genome-build-policy) for command scope, BED inheritance, source/output separation, and grouped diagnostics.

### `ldsc._genome.build_inspection`

The private [inspection module](../../src/ldsc/_genome/build_inspection.py) backs `ldsc infer-genome-build PATH` without adding a Python export.
It reads owned artifact metadata and relevant members, reuses shared table/column readers and HM3 inference, and renders declared, inferred, unknown, or separately labeled alternative builds independently of the SNP identifier mode.
It does not call workflow build resolution: unknown evidence is a successful inspection, while contradictory evidence and invalid inputs are errors.
The module also supplies the command's `PATH` argument and writes the report to stdout (`add_inspection_arguments()`, `run_inspection_from_args()`); it opens no failure-marker scope, and the CLI marks the command read-only so its errors give no workflow-log advice.
No artifact writers, workflow logs, or liftover participate.
`tests/genome/test_genome_build_inspection.py` covers artifact families, sampling, rsID reporting, exit codes, and read-only behavior through the CLI.

### `ldsc._runtime.workflow_logging`

This module owns the small shared logging context used by workflow wrappers.
It sets the `LDSC` logger threshold, attaches a per-run file handler to the
`LDSC` ancestor logger, writes lifecycle audit lines, captures a run-aborting
traceback into the log file on failure, restores logger state on exit, and
exposes log-only formatting helpers. It also owns the CLI console-handler
lifecycle: the per-run `.log` file is the authoritative sink, while the console
(stderr) is CLI-only and silent except for errors (it falls back to full
verbosity only for low-level calls without a workflow log). The Python API never
emits console output. Architecture invariant: workflows preflight the
log path with scientific outputs before entering the context, and result
`output_paths` mappings do not include logs.

The same module owns the small `RUN_FAILED` marker helper used at public
materializing boundaries. The helper adds failure visibility for authorized
overwrites but performs no rollback, restoration, quarantine, or action-order
change. See [logs-and-outputs.md](../reference/logs-and-outputs.md).
Each run has one marker scope, opened by its workflow (`materializing_overwrite_guard()` or `overwrite_failure_marker()`) from parsed arguments, or, for `build-r2-matrix`, by the `_build()` body that both entry points call, after `R2MatrixBuilder._inspect_inputs()` has assigned chromosomes from BIM contents.
A scope opened inside an open one neither writes nor removes, so the outermost public entry point decides the marker and its label.
The CLI opens no scope; `run_cli()` wraps dispatch in `cli_command()`, and the marker records that command line.
`serial_calculation()` returns a workflow log path opened in its worker, so a marker written in the calling process names that log.

### `ldsc._plotting.dispatch`, `ldsc._regression.h2_scale`

`ldsc._plotting.dispatch.plot_result()` is the sole public plot dispatcher. It selects
one approved plot from canonical `diagnostics/metadata.json`, follows declared
source paths, and imports private Matplotlib builders only on a plot-producing
path. `ldsc._regression.h2_scale.convert_h2_scale()` reads the observed estimate and SE
from a canonical `h2` result and reuses the kernel liability conversion factor;
only prevalence-range mode imports Matplotlib. Both default to package-owned
nested result roots and allow an advanced Python-only destination override.
Their complete seams are documented in [plotting.md](design/plotting.md).

For `rg`, the dispatcher loads the correlation table and optional metadata-declared single-trait heritability table.
Private builders align observed SNP heritability and SE by trait name for the heatmap diagonal or anchor annotations.
Missing heritability values display `failed`; the plotting layer neither reruns regression nor chooses among pair-specific heritability fits.
See [the annotation contract](design/plotting.md#heritability-annotations-in-rg-plots).

### `ldsc._annotation.builder`

This is the public annotation API and CLI dispatch seam: `AnnotationBuilder`, `AnnotationBundle`, `AnnotationResult`, `run_annotate()`, `parse_annotate_args()`, `run_annotate_from_args()`, and `main()`.
`AnnotationBuilder.run(..., output_dir=...)` prepares a chromosome-backed dataset with explicit ownership.
`_annotation.annotate` owns standalone BED/gene validation and incremental canonical output.
`_annotation.bim_grid` stages metadata-only BIM grids through `_annotation.sources.stage_retained_rows()`, the global identity cleanup that annotation-file preparation also uses.
Annotation-grid input emits aligned baseline/query families; BIM input emits query values only.
`_annotation.saved_results` records these exports in the returned `AnnotationResult`, which owns no files; `prepared_in(directory)` reopens them, without original-input dependencies, as a bundle that owns its private folder in the directory the caller chooses.
`_annotation.preflight` holds the query checks that `run_annotate()`, `ldscore` and `build_annotation_shards()` (behind `AnnotationBuilder` and `quantile-h2`) share: BED and gene-list query names that are reserved or collide stop before any grid is read, and the gene-list Gate A error and the projection genome build come from one place; `ldscore` resolves its own build against the reference panel.
Gene-list Gate A is one call in all four gene-list workflows (`annotate`, `ldscore`, indexed scoring and `build_annotation_shards()`): `_annotation.gene_query_storage.resolve_gene_lists_staged()` resolves the lists and applies the query-name check against its required `existing_columns`, and each workflow keeps its own failure handling.
Each preparation path turns its `PreparedAnnotationSources` into a bundle with `PreparedAnnotationSources.bundle()`.
The other `_annotation` modules own bounded discovery, identity cleanup, separate metadata/value artifacts, projection, and storage lifetime; `_ldscore.direct_annotation` integrates them with direct LD-score scope and support gates.
The complete module and ownership map is in [memory-and-cpu.md](design/memory-and-cpu.md).

For interval projection, `_annotation.interval_queries.prepare_bed_queries()` streams validated BED rows from `_core.regions.iter_bed_rows()` into chromosome spools.
`_annotation.interval_queries.build_query_shards()` uses `_annotation.interval_projection.ChromosomeProjector` for both BED and gene intervals, preserving zero-based half-open membership against one-based SNP positions.
`annotate` measures annotation-row gene support in the same pass by passing a `gene_support` object, which also sets the gene query statuses; the other callers only project.

### `ldsc._annotation.gene_lists`, `ldsc._annotation.query_outcomes`

These internal helpers validate the user-supplied or index-embedded gene coordinate catalog, resolve exact identifiers, and apply ordered query-status rules.
`_annotation.gene_query_storage` shares gene-coordinate-catalog lookups, stages per-source selections and complete audits, and returns compact summaries.
Gene-coordinate-catalog build selection, scientific-query pruning, and final diagnostic publication remain with the consuming workflow.

### `ldsc._refpanel.builder`

This module builds R² matrices from PLINK inputs.
It handles optional genetic-map loading, optional chain-file liftover selection and liftover-stage rows of the dropped-SNP table, then delegates pairwise-R² emission to the kernel.
`_refpanel.outputs.R2MatrixBuildDirectoryWriter` owns the output layout: its `artifact_family()` declares the files a run writes and owns (one chromosome's files, with chromosome-scoped metadata, log and failure-marker names, for a concrete single-chromosome run; the whole directory for an `@` suite run), and it writes the build metadata and the dropped-SNP tables in the shared identity-drop columns.
Each chromosome's SNPs come from `_refpanel.plink_source.select_plink_chromosome()`, the selection stage `ldscore --plink-prefix` also uses (SNP restriction, then identity cleanup of the rows it keeps, and the keep list mapped to `.fam` rows); the builder adds liftover and duplicate-position filtering on the rows it returns.
The SNP restriction and the keep list are read once per run, in `_prepare_build_state()`.
PLINK source build is explicit or inferred locally; a matching chain emits the opposite build only in `chr_pos`-family modes, while `rsid`-family builds are source-build-only.
Coordinate duplicate handling is `chr_pos`-family behavior and always uses drop-all for source/target collision groups; processed chromosomes always write `diagnostics/dropped_snps/chr{chrom}_dropped.tsv.gz`, header-only when no liftover-stage rows dropped.
Named MHC and centromere exclusion is not part of reference-panel construction; `ldscore` and `build-gene-ldscore-index` apply it only to regression/output SNPs.
Architecture invariant: emitted parquet schemas are part of the public file contract for parquet-backed LDSC3 workflows.

### `ldsc._refpanel.query`

This module is the public pair-query surface for package-built R² matrices.
`R2Matrix.open()` accepts a `build-r2-matrix` output directory or its genome-build child and checks every chromosome's files and Parquet footer through `R2MatrixChromosome`; each sidecar and its binding are validated, and per-chromosome index state cached, when a query first touches the chromosome.
`R2Matrix.query_pairs()` and `query_r2_matrix()` accept endpoint-suffixed pair tables and append `r2`, nullable `sign_r`, signed Pearson `r`, and `status`; low-level Parquet lookup lives in `ldsc._refpanel.r2_format`.
Architecture invariant: `query-r2-matrix` reads canonical package-built R² matrices only; it does not define a second external R² format.

### `ldsc._ldscore.workflow`

This module owns the `ldscore` command and the keyword-argument Python workflow `run_ldscore()`: the parser, `main()`, `run_ldscore_from_args()`, and the direct and indexed run bodies.
It resolves annotation and reference-panel inputs, selects the gene-coordinate-catalog projection build, synthesizes the all-ones `base` annotation when an unpartitioned run omits baseline/query inputs, measures gene-list support, and hands the calculation to `LDScoreCalculator`.
`run_ldscore_from_args()` is a thin adapter: `_configs_from_args()` checks the input-mode rules and builds the run's configuration before anything is written, for the command line and `run_ldscore()` alike, and a typed body runs it.
A direct run's `_DirectConfigs` holds the `GlobalConfig` with the requested genome build, `RefPanelConfig`, `LDScoreConfig`, `AnnotationBuildConfig` and `LDScoreOutputConfig`; an indexed run's `_IndexedConfigs` holds its gene lists and execution settings.
An option counts as supplied when the record of explicitly supplied options (`_explicit_cli_options`) names it, never by comparison with a parser default, so the registered `GlobalConfig` that a keyword call carries is not a live override in indexed mode.
The configurations' own range checks run once; `_OPTION_RANGE_ERRORS` words an error that a configuration raises by field name with the command-line option instead.
`_run_direct_ldscore()` preflights the outputs, runs the input gates, resolves the genome build once against every coordinate input (`_resolve_build()`), and passes typed configurations to `_ldscore.preflight` and `_ldscore.direct_annotation`; `_run_explicit_indexed_ldscore()` runs an indexed request.
A direct gene-list run measures each gene's support on the prepared reference SNPs before any LD score is computed (Gate B, `_apply_direct_gene_gate_b()`): it prepares the chromosomes one at a time in the parent and counts with `ChromosomeProjector.counts()`, the padding rule of annotation projection.
The direct run then hands the calculation a private `_RunContext`: the resolved `GlobalConfig`, the reference panel built with it, and which steps have run (output preflight and workflow log, build resolution, R² matrix checks), so `LDScoreCalculator._run_prepared()` repeats none of them.
`LDScoreCalculator.run()` passes a context that has settled nothing; the calculation then resolves the build against its prepared inputs, checks every R² matrix chromosome, and preflights and opens `diagnostics/ldscore.log` itself.
The parent builds the reference-panel adapter with `_refpanel.adapters.resolved_ref_panel(config, spec, chromosome_prefixes=...)` and the resolved configuration, and the direct LD-score producer sends every task one adapter built the same way with the run's configuration and the resolved PLINK chromosome prefixes, so a run resolves its build once; `resolved_ref_panel()` checks no build declaration again, and opening an R² matrix chromosome still compares its declared build with the run's.
`RefPanelLoader.load(spec)`, for direct Python callers, checks the panel and `ref_panel_snps_file` declarations in a coordinate-family mode (resolving an `auto` build) and then builds the adapter through `resolved_ref_panel()`.
PLINK gene LD-score index construction builds its panel through the same `resolved_ref_panel()`.
An indexed run tells `run_indexed_ldscore(main_log_open=True)` the same way that its log is open, so `_ldscore.diagnostics.ldscore_logs()` never has to recognize a caller's log handler.

### `ldsc._ldscore.calculator`

This module holds the public calculation service `LDScoreCalculator` and the `_RunContext` it receives; it orchestrates chromosome-wise LD-score computation over prepared inputs.
`LDScoreCalculator._run_prepared()` resolves the build of prepared inputs (`_resolve_calculator_build()`), checks the reference files (`RefPanel.check_reference()`) and opens the workflow log only where the context says the caller has not.
It then opens the chromosome logs and the chromosome pool, logs the reference-panel CM and MAF sources that `RefPanel.reference_metadata_sources()` reports, and runs the direct LD-score producer.
The direct LD-score producer, `_DirectProducer` in `_ldscore/direct_batches.py`, mirrors `_IndexedProducer` in `_ldscore/gene_index/batches.py`: each task computes one chromosome from the reference panel.
For BED and gene-list queries, `prepare_direct_annotations()` returns the bundle with the query columns and statuses declared but not projected, plus a `DeferredQueries` record (`_annotation/interval_queries.py`) that `LDScoreCalculator._run()` hands to the producer; each batch projects its own columns through `execution_query_bundle()`.
The producer passes each chromosome's annotation identity-cleanup drops to its first-batch task explicitly, so the bundle carries no LD-score execution state.
Every task, in the first and in later query batches, runs `compute_direct_chromosome()`, which prepares the chromosome, projects the batch's annotations with the kernel once, and builds one record from the kernel's `ChromComputationResult`.
The no-output path receives a `ChromLDScoreResult`; a writing run's task stages its post-kernel float32 query scores in a private Parquet fragment and returns a `DirectChromosomeFragment` with baseline data (first batch only), scientific statistics, diagnostic descriptors, and exact variation summaries.
A chromosome whose reference preparation leaves no SNP returns a `_SkippedChromosome` instead.
`LDScoreCalculator.compute_chromosome()` is a thin public wrapper around `compute_direct_chromosome()`.
`_ldscore.batch_execution` runs the query batches; the parent combines statistics and filters queries after every chromosome succeeds, then yields `StagedLDScoreBatch` to the shared writer without constructing a genome-wide query matrix.
Prepared no-output calculations run their one batch through the same LD-score producer tasks and retain the materialized `LDScoreResult` path.
Writing workflows retain `diagnostics/ldscore.log` and per-chromosome logs.
Query-local BED/gene failures are status records; a run continues when at least one query is usable.
Chromosome workers return only the typed `EmptyReferenceSNPs` outcome as a skip, carrying the error's reference-SNP funnel; other failures propagate independently of diagnostic wording.
The parent fails the run, naming every such chromosome, when the `annotation_snp_grid` step emptied a chromosome (`annotation_snp_grid_emptied()` in `_refpanel/funnel.py`), and otherwise logs the skip as a WARNING; the batch-1 funnels become `diagnostics/reference_snp_funnel.tsv`.
Runs without annotation files skip such a chromosome earlier: `_ldscore.direct_annotation.prepare_synthetic_base()` builds the synthetic `base` from `RefPanel.load_metadata()`, which raises `EmptyReferenceSNPs` with the funnel when the panel's filters leave no SNP; the chromosome gets no task, its skip is logged as a WARNING, and `_DirectProducer` adds its funnel to the batch-1 funnels.
Architecture invariant: computation stays chromosome-wise, global filtering uses the complete chromosome summaries, and final publication follows complete success.

### `ldsc._refpanel.adapters`, `ldsc._ldscore.compute`

`RefPanel.prepare_chromosome(chrom, annotations, config)` is the shared preparation boundary for `LDScoreCalculator.compute_chromosome()` and PLINK gene LD-score index construction.
The PLINK/parquet adapters own source resolution, identity and SNP restrictions, sample/MAF filtering, annotation alignment, authoritative reference CM/MAF, LD-window validation, and reader policy.
`PlinkRefPanel._load_source()` selects a chromosome's rows through `_refpanel.plink_source.select_plink_chromosome()`, the stage `build-r2-matrix` uses, and records the reference-SNP funnel from the counts it returns.
Each panel reads `--ref-panel-snps-file` (both backends) and `--keep-indivs-file` (PLINK) once, on first use; `_DirectProducer` resolves them on its task panel in the main process (`RefPanel._resolve_run_inputs()`), so the copies sent to chromosome workers carry them.
They return `PreparedChromosome`: metadata, a selected annotation accessor mapping retained reference SNPs to stored annotation rows, window bounds, and one owned reader.
The annotation rows are the annotation SNP grid: `_align_annotations()` excludes reference-panel SNPs without an annotation row, and `_refpanel.funnel.ReferenceSNPFunnel` records each step's removals as `PreparedChromosome.reference_snp_funnel`, or on `EmptyReferenceSNPs` when no reference SNP remains ([reference SNPs on the annotation SNP grid](../reference/ldscore-snp-universes.md#reference-snps-on-the-annotation-snp-grid)).
Packed binary and dense float32 annotation stores remain behind that accessor.
PLINK preparation preserves physical BIM/BED row indices through sorting and filtering.
R²-matrix preparation opens each chromosome once through `R2MatrixRefPanel.open_chromosome()` and builds its reader from that handle, so the contract check, the single sidecar parse and the sidecar binding follow the same path in production and adapter tests.
`R2MatrixRefPanel.check_chromosome()` runs the footer and sidecar-binding checks without preparing the chromosome, optionally reading every stored pair; the `ldscore` input preflight (`_ldscore.preflight.inspect_reference_inputs()`) and `R2MatrixRefPanel.check_reference()`, which `LDScoreCalculator` calls, both use it.
The calculation picks the kernel from the prepared reader: a `SortedR2BlockReader` streams R² pairs, and a PLINK genotype reader is rewound (`rewind()`) before its genotype block sums.

`_ldscore.compute.compute_chromosome()` projects baseline/query annotations and regression weights together, then derives counts and overlap. It consumes prepared state without reopening sources. Use prepared state in a `with` block: readers close on success or failure, and metadata caches retain tables only. Separate `load_metadata()` inspection may read PLINK data independently; one chromosome calculation itself constructs one BED reader. `tests/refpanel/test_reference_preparation.py` covers known raw/unbiased results, alignment, authoritative metadata, physical row mapping, and failure cleanup; `tests/ldscore/test_reference_snp_funnel.py` covers funnel counts on both backends, the grid-emptied failure, and the published diagnostics; `tests/ldscore/test_ldscore_workflow.py` verifies the public single-reader path.
`tests/ldscore/test_restriction_before_identity_cleanup.py` checks that `ldscore --plink-prefix`, `build-r2-matrix` and `ldscore --r2-matrix-dir` keep the same reference SNPs from one panel and one `--ref-panel-snps-file` when no liftover or duplicate position intervenes, that a `.bim` row whose position is not a whole number stops `ldscore --plink-prefix` and `build-r2-matrix` before computation, and that a chromosome emptied by the restriction or by identity cleanup of the rows it kept is skipped, with or without annotation files.

### `ldsc._ldscore.legacy_converter`

This workflow module is the sole public reader of selected legacy LD-score fragments.
It discovers one coherent autosomal reference family and one weight family, joins them by rsID, accepts either an unpartitioned single-score suite or a complete full baseline suite, validates the fixed legacy common-SNP-count semantics, and writes an ordinary LDSC3 LD-score directory through `LDScoreDirectoryWriter.write_batches()` as a one-item batch, so its files are staged and the published directory is reloaded as in an `ldscore` run.
Baseline conversion reconstructs and validates marginal counts plus the full overlap matrix; it never accepts query or thin annotations.
Inputs remain read-only and provenance includes selected prefixes/files, intersections, coordinate evidence, and count origins.
Routine source hashes are not computed.
Architecture invariant: regression and the computation modules never parse legacy LD-score fragments.

### `ldsc._ldscore.gene_index.build`, `ldsc._ldscore.gene_index.atoms`

The workflow module owns the `build-gene-ldscore-index` command, explicit hg19/hg38 plus any SNP identifier mode, an explicitly selected cM/kb/SNP window, the baseline/PLINK consistency gate, canonical `index_id`, logging, locking and staging in the index's hidden `.build/` directory, incremental chromosome persistence in the staged index, publication with root `metadata.json` written last, and reload validation of the staged and published index.
`run_build_gene_ldscore_index_from_args()` turns the parsed arguments into a `GeneLDScoreIndexBuildConfig`, and `build_gene_ldscore_index()` takes one; both call the private typed body `_build(config, genetic_map_sources=..., overwrite=..., log_level=...)`, because the genetic-map sources, `overwrite` and the log level are not config fields.
`gene_index.artifact` owns the index format, writing and reading: strict mode/build/digest validation of a saved index is its `load_gene_ldscore_index()`, which the build also uses to reload-validate the staged and published index; explicit indexed `ldscore` assembly is `run_indexed_ldscore()` in `gene_index.scoring` with the chromosome workers in `gene_index.batches`.
Each worker atomically installs its payload in a distinct staged chromosome directory and returns compact evidence; the coordinator writes all shared and identity metadata, validates the complete stage, and moves its components into the index directory, root `metadata.json` last.
The open log stays under `<index-dir>/.build/`; after its handler closes, a successful log moves into the published `diagnostics/` and `.build/` is removed.
Workers record their work in `.build/chromosome-logs/chrN.log` through `_ldscore.diagnostics.chromosome_log()`, the per-chromosome log of `ldscore` workers, and take their identity settings from `artifact.index_global_config()`; the coordinator logs each chromosome's completion and evidence-based warnings in the main log, and a successful build moves the chromosome logs to `diagnostics/chromosome-logs/`.
A failure before publication keeps an earlier index loadable, and build-directory or log-placement cleanup after publication is warning-only.
Shared allele-aware cleanup and duplicate-group removal retain audit diagnostics.
Before chromosome workers start, `_validate_baseline_plink_inputs()` compares each chromosome's baseline rows with the PLINK rows `_refpanel.plink_source.select_plink_chromosome()` keeps: contradictory coordinates or usable alleles and empty matches fail with a repair table, cross-chromosome rsID checks use a temporary SQLite join, and the gate counts `chr_pos` label disagreements.
Each worker passes its baseline shard, in canonical (`POS`, `SNP`) order, to `PlinkRefPanel.prepare_chromosome()` as the annotation SNP grid, as direct `ldscore` does, so the panel's selection, identity cleanup and grid match are the only ones; the chromosome's `baseline_plink_intersection` counts come from the reference-SNP funnel and the identity drops (`_baseline_plink_evidence()`), PLINK `CHR/POS/SNP/A1/A2` metadata is authoritative, and physical BIM/BED rows are preserved.
The coordinator builds one `PlinkRefPanel` per run, reading the keep list once, and sends it to the workers.
The private kernel owns disjoint half-open atoms, Boolean gene-to-atom CSR membership, bounded SNP-by-atom blocks, sufficient statistics, and float64 `Y @ z` assembly.
Each index directory is one immutable format-version-6 artifact with no restart, resume, incremental update, component reuse, or legacy-index loader.
The builder computes each chromosome's baseline LD scores, regression weights, counts, overlaps and `annotation_types` with direct `ldscore`'s `_ldscore.compute.compute_chromosome()`, then traverses the genotypes again for each atom batch of the operator; it builds its root `annotation_types` and `snp_universe_policy` with the helpers direct `ldscore` uses (`_annotation.semantics.combine_annotation_types()`, `_ldscore.outputs.snp_universe_policy_record()`), so indexed metadata matches a direct run.
Earlier formats require a complete rebuild.
Indexed output is still a self-contained LD-score directory written by `LDScoreDirectoryWriter`; it inherits identity/build and rejects live overrides.
Architecture invariant: online assembly uses only the explicitly named complete index and never discovers an index or falls back to live computation.

Gene-query eligibility and finalization are shared by direct and indexed computation: the staged gene-list batch's `with_coverage()` in `_annotation/gene_query_storage.py` checks coverage after explicit exclusions, and `_annotation/query_outcomes.py` owns resolution/support reason precedence, required-control/all-skipped viability, and the final query statuses from counts and score variation (`finalize_query_measurements()`); the LD-score result records' `select_queries()` in `_ldscore/results.py` then prune the skipped queries' columns and metadata. `_ldscore/preflight.py` validates direct input chromosome contents, exact baseline/reference set agreement, complete `@` declarations, and authoritative glob matches. Neither helper is a general workflow framework; reference preparation and indexed atom support remain backend-specific. See [the detailed coverage contract](../reference/results/gene-list-diagnostics.md#chromosome-scope-and-pathway-coverage).
Both modes stop through the same sequences in `_ldscore/diagnostics.py`: `fail_gene_list_resolution()` for rejected gene-list identifiers or names, and `fail_unusable_queries()` for unusable queries or controls, which writes the query diagnostics, logs the written paths, gene-list rejections, SNP support and query statuses, and raises; `log_query_outcomes()` logs the same summaries after a successful run.
Runs that stop before publication describe their diagnostics with `_ldscore.outputs.QueryDiagnostics`, the only input of the writer's diagnostics methods; the direct workflow converts its annotation bundle to one before `fail_unusable_queries()`, and `LDScoreDirectoryWriter.artifact_family()` selects outputs by type: `QueryDiagnostics` for diagnostics alone, the private publication record for a complete directory.
After a failed input preflight, `LDScoreDirectoryWriter.write_input_preflight_diagnostics()` writes `input_issues.tsv` and `chromosome_scope.json` and leaves an earlier run's outputs in place, as for every input gate; the other diagnostics writers reconcile stale owned outputs.

### `ldsc._sumstats.munger`

This module wraps the legacy LDSC munging behavior in typed public objects such as `MungeConfig`, `RawSumstatsInference`, and `SumstatsMunger`.
The munged artifact has its own modules: `_sumstats.table` holds `SumstatsTable`, `load_sumstats()` and the Parquet footer codec, which regression imports without the munging workflow, and `_sumstats.outputs` holds `SumstatsDirectoryWriter`.
Raw heterogeneous GWAS input is named explicitly as `raw_sumstats_file` in Python and `--raw-sumstats-file` in the CLI.
`main()` parses only, `run_munge_sumstats_from_args()` maps CLI arguments into `MungeConfig`, and `SumstatsMunger.run()` owns path resolution, raw-format inference, fixed-output preflight, diagnostics under `diagnostics/`, and result construction before delegating primitive munging work to the kernel.
`SumstatsDirectoryWriter.artifact_family()` declares the owned family once for the run and for `SumstatsMunger.write_output()`: both data formats of one filename stem, plus the log, the dropped-SNP table and the input-issue table for a run.
`_sumstats.outputs.curated_sumstats_frame()` puts the kernel's rows into the one frame that the returned `SumstatsTable` and both data files hold (fixed columns, canonical `CHR` labels on complete coordinates, genomic row order); `SumstatsMunger.run()` and `write_output()` build it once before writing.
`SumstatsDirectoryWriter.write()` writes that frame: it stages the Parquet file (only split into chromosome row groups), the `.sumstats.gz` file and the dropped-SNP table inside the output directory, moves them into place after all are complete, and then removes stale owned siblings.
`--input-format auto` is the default and detects common plain text, including VCF-style headers, old DANER, and new DANER; `--infer-only` prints the inference report and suggested minimal command without writing artifacts, while still requiring `--output-dir` for the uniform CLI contract.
The workflow keeps raw-input parsing flexible, skips leading `##` raw metadata lines, writes canonical `CHR` and `POS` columns alongside `SNP`, `Z`, and `N`, loads a headered `sumstats_snps_file` keep-list or packaged HM3 restriction once before parsing, applies that restriction inside each retained chunk through shared identifier readers, and can liftover `CHR`/`POS` after filtering when a `chr_pos`-family mode has a resolved source build, an output build, and exactly one required method.
In `chr_pos`-family modes, chunk filtering uses source-build coordinates and compact packed coordinate keys.
HM3 quick liftover requires the packaged HM3 restriction flag.
Liftover preserves `SNP` labels, drops missing or invalid coordinates before mapping, drops duplicate source-coordinate groups before mapping, drops unmapped/cross-chromosome rows, and removes duplicate target-coordinate groups after mapping.
The `sumstats.parquet` footer persists `ldsc:artifact_type`, `ldsc:snp_identifier`, final output `ldsc:genome_build`, and optional `ldsc:trait_name`; rsID-family summary statistics store an empty `genome_build`.
`_sumstats.table.encode_sumstats_footer()` and `_decode_sumstats_footer()` form one pair, and loading rejects a footer in any form the encoder cannot write.
When the run's QC thresholds are known, the footer also records `ldsc:qc_info_min` and `ldsc:qc_maf_min` as provenance that readers ignore (see `docs/reference/formats/artifact-metadata.md`).
Row-level identity and liftover drops are audited in `diagnostics/dropped_snps/dropped.tsv.gz`; detailed coordinate provenance, liftover reports, HM3 provenance, output files, and row-group layout are readable audit log entries in `diagnostics/sumstats.log`.
DANER-specific schema handling can be inferred or requested explicitly: old DANER reads case/control sample sizes from `FRQ_A_<Ncas>` and `FRQ_U_<Ncon>` headers, while new DANER reads per-SNP case/control counts from `Nca` and `Nco` columns.
`MungeRunSummary.output_paths` lists munged data files and the dropped-SNP audit sidecar, not `diagnostics/sumstats.log`.
Public interface: downstream workflows should consume `SumstatsTable` or the self-describing munged `sumstats.parquet` (or `.sumstats.gz`), not raw heterogeneous GWAS tables.

`_sumstats.munge_input.resolve_munge_plan()` resolves the munging plan, every check and decision a run makes before its munging pass, and returns a `MungePlan`, which `SumstatsMunger.run()`, `--infer-only` and `infer_raw_sumstats()` all use, so the previews report the problems the run stops on before that pass.
Within it, `_sumstats.munge_input.resolve_raw_schema()` owns raw header interpretation: format detection, column aliases and hints, DANER/sample-size selection, the signed statistic, QC thresholds, allele orientation and every header-level issue, returned as a `RawSchema`.
The munging plan then resolves the keep-list path, checks `INC_ALLELE`/`DEC_ALLELE` orientation against the signed statistic, infers or checks the source build and selects the liftover request for that build; it stops at the first of these stages that finds a problem and collects that problem as one more issue, keeping the exception a check raised so that the run re-raises it unchanged.
`resolve_source_config()` checks the build and liftover options against the SNP identifier mode first, and the run also checks the liftover method of an explicit source build before any input IO.
The run's input preflight (`_runtime.preflight.inspect_declared_inputs()`) checks the keep-list path and the chain file of an explicit source build before the munging plan, so the run stops on those paths before any header problem and reaches the munging plan's path issues only for a chain file selected after `auto` resolves; the previews report them all from the munging plan.
`_sumstats.munge_input.prepare_munge_input()` raises the munging plan's issues and prepares the keep-list.
When that request selects HM3 quick liftover, it loads the packaged HM3 map once (`_sumstats.liftover.load_hm3_map_file()`), derives the keep-list keys from it, and passes the same frame to `_sumstats.liftover.apply_sumstats_liftover()` through `ResolvedMungeInput.hm3_map`.
`_sumstats.liftover.apply_sumstats_liftover()` executes the liftover request of the munging plan and makes no option checks of its own: a request without a method is a no-op, and a request whose shape the munging plan never builds is an internal error.
It supplies `ResolvedMungeInput` to `_sumstats.munge_kernel.munge_sumstats()`, which returns `MungeResult` with data, exclusive stage counts, coordinate provenance and separate drop tables. The sole CLI parser is `_sumstats.munger.build_parser()`. Production and tests use the same preparation path; results never travel through added namespace attributes. `tests/sumstats/test_munge_results.py` covers count conservation, chunk independence, bounded build inference before QC, request reuse, and liftover/identity accounting; `tests/sumstats/test_munge_infer_diagnostics.py` and `tests/sumstats/test_munge_allele_policy.py` cover the shared header problems and the allele policy.

### `ldsc._ldscore.source`

`read_ldscore_manifest()` reads and checks `metadata.json` without LD-score values and returns an `LDScoreManifest`: identity, declared files named inside the directory (symbolic links allowed), baseline and query row groups and schemas from Parquet footers, counts, `count_config` and its agreement with `overlap_config`, and a lazy `load_overlap()`; read and decoding failures become `LDSCInputError`.
`load_ldscore_from_dir()` calls it and opens an `LDScoreSource`, retaining baseline values, identities, counts, and overlap while leaving query LD-score values on disk.
This module owns manifest/schema checks, count/overlap consistency, and table/identity validation shared by materialized and staged output.
It also declares the schema facts that the direct and indexed LD-score producers, the legacy converter and the readers share: the weight column `REGRESSION_LD_SCORE_COLUMN`, the count-record fields `ALL_COUNT_FIELD` and `COMMON_COUNT_FIELD` with the one record constructor `count_record()` (an unavailable all-SNP count is recorded as `null`, an unavailable common-SNP count is left out), the `count_config` fields `COMMON_MAF_MIN_FIELD` and `COMMON_MAF_OPERATOR_FIELD`, and the index provenance fields `INDEX_PROVENANCE_FIELDS`.
The diagnostic paths that producers build outside the writer (`INPUT_ISSUES_FILE`, `LD_WINDOW_SNP_COUNTS_FILE`, `REFERENCE_SNP_FUNNEL_FILE`, `chromosome_log_file()`) are declared next to the other file names in `_ldscore/outputs.py`, and the writer takes the `query_diagnostics` paths of `metadata.json` from its own `artifact_family()`.
`validate_ldscore_baseline_table()` validates the complete shared baseline; `validate_ldscore_query_table()` checks a materialized query table against its already-validated baseline rows, and `_ldscore.outputs._write_query_fragments()` requires each staged fragment's identity columns to equal the shared first-batch baseline rows exactly.
`LDScoreSource.read_queries()` returns explicit, uncached selections in requested column order; callers own their memory and must keep the source directory available and unchanged.
The batch output writer and regression workflows both call this loader, and `quantile-h2` reads only the manifest; `_ldscore/source.py` has no dependency on output writers or regression orchestration.
Selected-value validation and alignment to a regression model remain in `ldsc._regression.alignment`.

### `ldsc._regression.runner`, `ldsc._regression.alignment`, `ldsc._regression.fitting`, `ldsc._regression.partitioned`, `ldsc._regression.rg`, `ldsc._regression.commands`

The regression workflow is split into modules that each import only the ones listed before them: `results` (the public dataset and result types), `alignment` (trait preparation, SNP identity resolution and keying, the genome-build check, and h2 and rg dataset assembly), `fitting` (the fit policy, the h2, rg and per-query fits, result rows and metadata), `partitioned` and `rg` (the partitioned-h2 batch and the rg pair run), `commands` (the three commands and their Python forms), and `runner` (the public `RegressionRunner`).
`RegressionRunner` is a facade: it holds a `GlobalConfig`, the SNP identifier fallback for an input without recorded provenance, and a default `RegressionConfig`, and its methods forward both to these modules, which take them as arguments; `commands` calls the same modules directly and does not import `runner`.
The `config or self.regression_config` fallback exists only in the facade.
The `GlobalConfig` banner is an entry-point concern: each command and each facade method that reads inputs logs it once, and the regression modules log none.
The spine that the commands, the facade and the `rg` pair run share, each function with one explicit `RegressionConfig`, is `_prepare_h2_dataset()` (one trait's h2 dataset plus the legacy projection drops and `dropped_nonfinite_statistic_rows` that `h2` results record), `_h2_dataset()` (the same from prepared trait summary statistics), `_prepare_rg_dataset()` (a trait pair) in `alignment`, and `_fit_h2_dataset()` and `_fit_rg_dataset()` in `fitting`.
The standalone fits `_fit_h2_dataset()` and `_fit_rg_dataset()`, the baseline-only partitioned fit and `rg.estimate_rg_pairs()` are module-level `serial_calculation` callables, because the controlled child re-imports a callable by module and qualified name.

Regression accepts saved sources or prepared `LDScoreResult` values and prepares each trait's summary statistics once per run (`_PreparedTraitSumstats`, built by `alignment._prepare_trait_sumstats()`).
That preparation projects legacy summary statistics by rsID onto the regression SNPs in the LD-score directory with allele validation, orientation (A/T and C/G rows keep their source orientation), and a stable drop audit, and drops rows whose `Z` or `N` is missing or infinite with a warning and a recorded count (see [missing or infinite statistics](../reference/regression-options.md#16-missing-or-infinite-statistics)).
A missing or infinite baseline LD score or regression weight stops the run when `_PreparedLDScores` is built, once per run; a damaged query column fails only its model in `_dataset_from_query_values()`.
Each command then applies its own SNP alignment: `h2` and `partitioned-h2` merge by SNP identity only, because they fit $\chi^2 = Z^2$, and `rg` harmonizes trait 2 to trait 1 in `_rg_dataset()`, dropping A/T and C/G SNPs (strand unknown) and SNPs whose traits list different allele pairs, with both counts recorded per pair.
The shared `_select_model_columns()` then drops zero-variance LD-score columns and selects the count vectors on the final rows.
The commands require `output_dir`, write per-command logs under `diagnostics/`, and reach the regression kernel through `fitting` for SNP heritability, partitioned heritability, and genetic correlation.
Regression merges by the effective key for the active mode: `SNP` in `rsid`, `SNP:<allele_set>` in `rsid_allele_aware`, `CHR:POS` in `chr_pos`, and `CHR:POS:<allele_set>` in `chr_pos_allele_aware`.
Regression labels traits from explicit CLI overrides, then the `sumstats.parquet` footer `ldsc:trait_name`, then filenames; `rg` anchor selection uses `--anchor-trait`, matching trait labels before source paths.
`partitioned-h2` produces overlap-aware annotation summaries (legacy LDSC `--overlap-annot` math) and auto-detects two regimes from the LD-score directory: a joint baseline fit over functional annotations when there are no query columns, and a cell-type-specific baseline-plus-query model per query when there are.
Query runs always write the per-query result tree; the retired `--write-per-query-results` flag is rejected.
It requires `ldscore.overlap.parquet` (assembled per model from the stored overlap matrix via `ldsc._ldscore.overlap_table`); both regimes share one public schema.
It also restores legacy LDSC's regression-time collinearity guard as a hard error (`LDSCInputError` when the LD-score design matrix's condition number exceeds 1e5), and logs a WARNING naming any zero-variance LD-score columns it drops before fitting.
`rg` applies the same guard, with the same threshold, to each trait pair's regression rows, as legacy LDSC did; a collinear pair is recorded as a failed pair, like other pair-level errors, while each trait's own `h2` fit stops the run.
Both skip the guard when the LD-score directory stores no overlap matrix or the model has fewer than two LD-score columns.
`--allow-identity-downgrade` is regression-only and permits same-family allele-aware/base mixes to run under the base mode; rsID-family and coordinate-family modes never mix.
Rows keyed for a merge are cached per effective SNP identifier mode and downgrade cleanup on the prepared trait and on `_PreparedLDScores`, at most twice each.
`estimate_rg_pairs()` prepares every trait once and runs one joint genome-build check; each pair keeps its own effective mode and builds its merged dataset inside the pair's error handling, so a pair-level error, including a keying error in a downgraded pair, fails only that pair.
A pair whose jackknife deletion is singular is marked `unestimable`; the rg fit (`_fit_rg_dataset()`, behind `estimate_rg()`) adds the failed component, panel genomic spans and unsupported columns to the error, and the error carries the pair's counts as `fit_summary` for its metadata.
A pair whose heritability estimate is not positive is also `unestimable`, without an exception: the `RG` kernel ported from legacy LDSC flags it with `_negative_hsq` and the string `'NA'` for `rg`, `rg_se`, `z` and `p`, and `_summarize_rg_pair()` reads the flag `h2_out_of_bounds`, keeps the fitted counts and estimates and records the estimate that is not positive in `error`.
`_regression/kernel_values.py` is the one place that replaces the kernel's `'NA'` sentinels with NaN, where the package receives a kernel result: `numeric_ldsc_fit()` (`intercept_se` of a fixed intercept, `ratio` and `ratio_se` when the mean $\chi^2$ is at most 1) and `numeric_rg_fit()` (the `RG` estimates and `h2_out_of_bounds`), which `_fit_h2_dataset()` and `_fit_rg_dataset()` apply, and `numeric_overlap_table()` (`Enrichment_p` of a category containing every SNP), which `overlap_model.overlap_aware_category_table()` applies before it renames the table once to the `partitioned-h2` column names; the estimators themselves are unchanged.
`build_dataset()` and `build_rg_dataset()` are the public dataset builders, and `estimate_h2()` and `estimate_rg()` fit the datasets they return; the alignment that `estimate_partitioned_h2_batch()` keeps across query models (`_prepare_h2_inputs()`, `_dataset_from_prepared()`, `_PreparedRegressionInputs` in `alignment`) is private.
The run log's SNP identifier line comes from the prepared datasets, not from a prediction: `_log_regression_identity()` in `alignment` reports the effective mode and any downgrade after preparation in `h2` and in the `partitioned-h2` batch, and for `rg` per pair from the per-pair metadata, with the mode each input records.
`tests/regression/test_trait_preparation.py` and `tests/regression/test_nonfinite_statistics.py` cover single projection, keying counts, pair equivalence, the per-pair failure contract, and the non-finite drop.
Architecture invariant: regression only consumes aggregated LD-score directories; it does not recompute LD scores.

`run_h2()`, `run_partitioned_h2()`, and `run_rg()` in `commands` are the Python forms of the three commands, and `_regression.quantile_h2.run_quantile_h2()` is the form of `quantile-h2`.
Each passes its keyword arguments and the registered `GlobalConfig` to `_runtime.workflow_kwargs.run_from_kwargs()`, which builds the command's own parser, rejects unknown names, missing required options, values outside declared choices, non-boolean flags, and `log_level`/`snp_identifier`/`genome_build` before any write, takes the log level from the registered `GlobalConfig`, then calls the matching `run_*_from_args()` under a `RUN_FAILED` marker naming the Python function.
The `*_from_args()` functions remain module-level CLI dispatch targets and are not exported from `ldsc`.
Each is a thin adapter: it turns the parsed arguments into typed settings before any output preflight and calls the command's typed body.
`_configs_from_args(args, command)` returns the registered `GlobalConfig`, the SNP identifier fallback for an input without recorded provenance, and the `RegressionConfig`, built with the scalar prevalences for `h2` and `partitioned-h2` and without prevalences for `rg`.
`_run_h2()`, `_run_partitioned_h2()` and `_run_rg()` take those two configs and the other options (inputs, output directory, overwrite, log level, and the command's own options such as workers, query batch size or console notices) as keyword arguments.
`run_partitioned_h2()` checks `workers` and `query_batch_size` before it calls `run_from_kwargs()`, because their argparse types check them only on the command line; like an argparse error, an invalid value then writes nothing, not even the `RUN_FAILED` marker.
`run_h2_from_args()` returns an `H2Result`, `run_partitioned_h2_from_args()` the complete `PartitionedH2BatchResult`, and `run_rg_from_args()` the `RgResultFamily` with the paths written.
The three typed bodies share `_inspect_regression_inputs()`, which checks the declared summary statistics, prevalence manifest and LD-score directory before any input is loaded.
`h2` and `rg` publish from their typed bodies; `estimate_partitioned_h2_batch()` publishes its own result directory, because its query workers stage each fitted model inside it.
`run_quantile_h2_from_args()` follows the same pattern: `_regression.quantile_h2._configs_from_args()` builds the resupplied annotation sources as an `AnnotationBuildConfig` and the target as a `_regression.quantile_inputs.QuantileTarget`, which checks `--num-quantiles` and parses `--target-missing-value` once, and the typed body `_run_quantile_h2()` runs in one worker and passes them to `prepare_quantile_statistics()`; it builds no `GlobalConfig`, because the SNP identifier mode and genome build come from the fitted model's LD-score directory.
That directory is read only through `_ldscore.source.read_ldscore_manifest()`: after loading the fitted model, `_run_quantile_h2()` runs it and, as a separate check of the same input gate, the declared-file check of the regression commands, before any annotation is read, and passes the `LDScoreManifest` on for the SNP identity, common-SNP rule, counts and overlap.
`tests/regression/test_python_regression_workflows.py` checks that the Python and CLI routes write identical result files.

`RegressionFit` values in `_regression/results.py` carry the estimator, fitted dataset, realized block count, effective cap, two-step cutoff, and optional single-annotation diagnostic bins.
The `h2`, `partitioned-h2`, and `rg` reporting paths consume these values directly, and the public `estimate_h2()` and `estimate_rg()` return them.
`_fit_policy()` resolves the fixed intercepts and the requested two-step cutoff from a `RegressionConfig` once, and `_FitPolicy.two_step_cutoff()` applies legacy LDSC's single-annotation default of 30, for the `h2` and `rg` fits and their metadata.
Result metadata is built from the same values: `_dataset_metadata()` and `_fit_metadata()` record a fitted dataset's settings and a fit's counts, cap and cutoff for `h2`, `partitioned-h2` models and `rg` pairs, and `_ldscore_dir_of()` records the LD-score directory that was read.
Prevalences reach the fits checked: `RegressionConfig` and `estimate_rg_pairs()` apply `ldsc.config.normalize_prevalence_pair()`, so a quantitative trait is `(None, None)` throughout.
`estimate_rg_pairs()` reads the per-trait pairs from its `prevalences` argument only and rejects a `RegressionConfig` that sets `samp_prev` or `pop_prev`, which would otherwise be ignored; `ldsc rg` and `run_rg()` resolve their prevalence lists or manifest into that argument and build the config without them.
Output diagnostics remain owned by the workflow rather than being attached dynamically to numerical kernel objects.

Managed partitioned query workers, including serial execution, live in `_regression/partitioned_parallel.py`.
`partitioned.estimate_partitioned_h2_batch()` owns input preparation, batch reads and final publication, and hands the shared module-level `fitting.fit_partitioned_query()` calculation to the private module, which owns read-only numeric transport, spawned processes, bounded dispatch, log forwarding, and worker cleanup and imports neither `partitioned` nor `fitting`.
Workers borrow output-owned maps and stage unique complete-query results.
The parent assembles outcomes in query order before stable summary sorting.
`_runtime.parallelism` shares parsing, validation, and `resolve_worker_count()` with chromosome workflows, using the bounded query count as the work limit.
Tests in `test_regression_streaming.py`, `test_partitioned_query_failures.py`, `test_partitioned_query_workers.py`, and `test_worker_policy.py` cover numerical/artifact parity, CPU and memory ownership, strict/continue policies, and real process failures.
See [query workers](../reference/regression-options.md#43-query-workers-and-memory).

### Managed numerical execution

Regression model preparation borrows aligned read-only columns through `_dataset_from_query_values()` in `_regression/alignment.py`.
`_fit_h2_dataset()` applies the model's SNP filter before materializing its design, and the kernel reuses owned work buffers when dtype promotion permits.
`tests/regression/test_regression_memory.py` covers shared-input immutability, repeated fits, and independent weighted least-squares agreement.

`_runtime.parallelism` validates `workers` against the detected available CPU count from `joblib.cpu_count()`, then caps concurrency by chromosomes or complete query-plus-baseline models.
Its `chromosome_pool()` creates the spawn pool for direct scoring, indexed scoring and index construction, and turns a crashed worker into `LDSCInternalError` in all three.
There is no custom CPU detector, Slurm parsing, CPU binding, or automatic memory budget.
The default of one still uses one worker; direct and indexed LD reuse their pools across query batches, and regression reuses mapped-input workers across model batches.
Gene LD-score index construction now passes staged chromosome descriptors to spawned processes and collects compact persisted-component metadata.
The [focused local transport probe](validation/performance.md#gene-ld-score-index-transport-probe) records the migration evidence and additional RSS cost.

`_runtime.worker_runtime` configures and verifies 1 BLAS thread per worker for serial and parallel calculations, plus separate OpenMP and Arrow CPU/I/O limits of one.
Unsupported numerical controls produce a clear unverified/unmet record and allow computation to continue.
`_runtime.serial_calculation` puts standalone numerical entry points into a worker and reuses the executing worker for nested numerical calls, without creating nested task pools or changing the caller's global settings.
`_ldscore.diagnostics` owns worker-local chromosome file handlers; regression retains its existing pipe-based record forwarding.
All required chromosome work must succeed before aggregation and publication; logs survive failure while temporary scientific results are removed after executor shutdown.
See [parallel execution](../reference/parallel-execution.md) for the complete contract.

### Result writers

Each domain writes its own result directories on the shared `ldsc._runtime.artifacts` base: `ldsc._annotation.outputs` for `annotate`, `ldsc._refpanel.outputs` for `build-r2-matrix` (metadata, dropped-SNP tables and the owned-file declaration; `r2_format` writes the R² Parquet files and sidecars) and `query-r2-matrix`, `ldsc._ldscore.outputs` for LD scores and `ldsc._regression.outputs` for regression results.
`ldsc._ldscore.outputs` owns fixed files inside `output_dir`: `metadata.json`, `ldscore.baseline.parquet`, optional single or numbered query batch files, and an optional `ldscore.overlap.parquet` (the long-form annotation overlap matrix, written only
for runs with two or more annotation columns). Each LD-score file spans the calculated chromosomes with one row group per chromosome; root `metadata.json.query_batches` records the ordered query files, query columns, and row groups. Direct/indexed writing workflows return `LDScoreSource` without retained query tables. Older LD-score directories without the current manifest require regeneration.
`ldsc._regression.outputs` owns the regression results: for `h2`, root `h2.tsv` plus diagnostic
metadata. For `partitioned-h2`, it owns root `partitioned_h2.tsv` plus
the optional diagnostic `diagnostics/query_annotations/` tree containing
`manifest.tsv`, per-query `partitioned_h2.tsv`,
`partitioned_h2_full.tsv`, and `metadata.json`. For `rg`, it owns `rg.tsv`,
`rg_full.tsv`, `h2_per_trait.tsv`, and the optional diagnostic
`diagnostics/pairs/` detail tree; `rg` tables report nominal p-values only.
For a legacy summary-statistics input, each of the three writers also writes `diagnostics/dropped_snps/legacy_sumstats.tsv.gz` before `metadata.json` and lists it in `files`.
It
reuses existing directories but refuses existing canonical files unless
overwrite is supplied. Architecture invariant: public output customization
chooses the directory name and explicit overwrite policy, not per-run filename
prefixes.

Each directory writer exposes `artifact_family()` as the single declaration of its named outputs, conditional siblings, and derived roots. Workflows use it for early collision checks; writers use the final declaration for writing, metadata `files`, and stale cleanup. `ArtifactFamily` (in `_runtime.artifacts`) describes paths and delegates collision checks to `_runtime.paths`; it does not introduce a publication transaction. LD-score preparation no longer predicts a cleanup list that survives computation. See [logs and outputs, Output-Family Preflight](../reference/logs-and-outputs.md#output-family-preflight) and `tests/runtime/test_artifact_declarations.py` for the lifecycle contract and overwrite/reload checks.

### Query execution and result ownership

`stage_partitioned_h2_fit()` marks private regression files as owned so `PartitionedH2DirectoryWriter` can move them into its sorted publication tree.
Persistent `PartitionedH2FitArtifacts` default to copying.
The writer stores run metadata once, and `_regression.quantile_h2.load_fitted_partitioned_model()` resolves the query's relative `global_metadata` reference before applying model-specific fields.
`_regression/outputs.py` spells the fitted-model file names and the `global_metadata` reference once: `fitted_model_files()` locates one fitted model, its expected `artifact_type` and its manifest root for `quantile-h2` (the aggregate root only for a folder under `diagnostics/query_annotations`, so the root never leaves the directory given), and `published_query_artifacts()` names the published per-query files for `PartitionedH2BatchResult`.
`tests/regression/test_partitioned_publication.py` and `tests/regression/test_quantile_streaming.py` cover this writer-to-reader contract, relocation, and publication failure cleanup.

| Module and symbol | Responsibility |
| --- | --- |
| `_annotation.in_memory.bundle_from_frames` | Small prepared annotations and complete diagnostics in RAM, without filesystem staging |
| `_annotation.interval_queries.execution_query_bundle`, `DeferredQueries` | Project and release the active direct BED/gene batch from the producer's `DeferredQueries` record while sharing baseline/control and global preflight state |
| `_ldscore.batch_execution.run_query_batches`, `run_chromosome_batch`, `finalize_batch` | One active query batch at a time for either LD-score producer through the caller's persistent pool: per-batch statuses and scratch, bounded chromosome tasks, failure reporting, pool shutdown before scratch removal, and query-status finalization |
| `LDScoreCalculator.run` | Own the direct run's persistent worker pool and logs; the writing path runs `_DirectProducer` batches through `run_query_batches`, and owned fragment workspace remains live through consumption or worker shutdown |
| `_ldscore.direct_batches.compute_direct_chromosome`, `_DirectProducer` | One chromosome task of any direct query batch: compute once from the reference panel, stage post-kernel float32 query scores, and return baseline plus scientific/diagnostic summaries instead of query matrices; the producer combines them into one batch |
| `_ldscore.gene_index.batches.indexed_results` | Own the indexed run's persistent pool; `_IndexedProducer` reloads and releases one chromosome operator per task and reduces compact float64 statistics in index order |
| `_ldscore.query_variation.summarize_variation`, `combine_variation` | Exact global variation from first values and within-chromosome flags, preserving each LD-score producer's precision |
| `_ldscore.results.StagedLDScoreBatch`, `_ldscore.outputs._write_query_fragments` | Private batch descriptor, and its ordered validation/publication from direct float32 or indexed float64 fragments by the writer; delete each consumed fragment without a genome-wide query array |
| `_ldscore.outputs._LDScorePublication`, `LDScoreDirectoryWriter.write_batches` | Retain explicit shared publication state while releasing each materialized or staged batch; publish the LD-score directory after success and clean owned scratch on handled failure |
| `LDScoreSource.read_queries` | Explicit selected-query reads across manifest files in requested order, without caching or a width cap |

`LDScoreDirectoryWriter.write_batches()` consumes materialized `LDScoreResult` values or private `StagedLDScoreBatch` descriptors behind the same publication interface.
It is the only publication path for LD-score directories: direct `ldscore`, indexed scoring and `convert-legacy-ldscores` all publish through it.
Its publication record retains the baseline, counts, overlap, provenance, and diagnostics, without retaining the LD-score producer's query values or chromosome results.
`_LDScorePublication.from_first_batch` builds that record from the first batch, materialized or staged alike: both batch records carry the same run-level provenance fields, so the writer has no per-type copy.
The producer owns active fragments through consumption or worker shutdown; the writer closes the producer before removing publication staging on handled failure.
Fragment validation consumes tables directly and checks row alignment once per fragment, without constructing partial result objects.

`tests/ldscore/test_ldscore_query_batches.py` verifies numerical precision, zero writes for prepared single-batch calculation, batch lifetimes, publication, and overwrite cleanup.
`tests/ldscore/test_gene_index_streaming.py` covers operator ownership, known union scores, and serial/parallel agreement.
`tests/ldscore/test_batch_execution.py` covers, for any LD-score producer, the shared persistent pool, one batch in flight, per-batch status filtering and finalization, cancellation of pending tasks, and worker shutdown before scratch removal.
`tests/ldscore/test_indexed_batch_staging.py` covers operator and parent-frame lifetimes, exact pre-narrowing variation, skipped batches, and cleanup after row-validation failure.
`tests/ldscore/test_direct_batch_staging.py` covers real-kernel materialized/streamed parity, global filtering, parent-frame release, failure cleanup, and real PLINK gene audits, identity drops, chromosome scope, and reference exports.
`_annotation.query_outcomes.finalize_query_measurements` shares count/variation decisions across materialized direct results and both LD-score producers.
`tests/annotation/test_annotation_format.py` and `tests/annotation/test_annotation_storage.py` cover automatic exact binary classification, packed storage, bounded decoding, and logical column order.
Public annotation reads and final Parquet LD values stay float32; LD accumulation, annotation sums, and overlap products use float64.
Reference SNPs, regression SNPs, and fit-time validation are unchanged.
See [annotation format policy](design/memory-and-cpu.md#annotation-format-policy).

### Computation modules

Each domain keeps its numerical methods and low-level readers with the rest of the domain, rather than in one shared kernel layer.
Annotation table validation is in `_annotation.parsing` and BED parsing in `_core.regions`; PLINK and R²-matrix access in `_refpanel.plink` and `_refpanel.r2_format`; LD-score math in `_ldscore.compute`; legacy-compatible summary-statistics munging in `_sumstats.munge_kernel`; and the regression estimators in `_regression.estimators`.
Live interval projection resides in `_annotation.interval_projection`, called by `_annotation.interval_queries`.
Private boundary: these modules are not the supported import surface; numerical computation receives prepared inputs, while the `_refpanel.adapters` adapters resolve panel-specific sources before preparation.

## Cross-Cutting Concerns

- **Input token language**: public inputs accept exact paths, globs, and `@` chromosome suites.
  PLINK inputs additionally accept plain suite stems through `_runtime.paths.resolve_plink_prefix_group()`; this exception does not apply to other input families.
  Output paths are normalized but remain literal destinations.
- **Column and identifier normalization**: `_core/columns.py` is the single source of truth for raw-input aliases, including `#CHROM`/`CHROM` as `CHR`, internal artifact headers, SNP identifier modes, and genome-build aliases.
  The munging workflow owns higher-level raw-format inference and conservative repair suggestions before the kernel runs.
  In summary statistics `A1` is the effect allele, the allele the signed statistic refers to, not necessarily the genome reference allele; it is never inferred from `REF`/`ALT` or `ALLELE0`/`ALLELE1`/`ALLELE2` names, and `NEFF` is not inferred as `N`.
  `_genome/build_inference.py` owns automatic hg19/hg38 and zero-based/one-based inference for `chr_pos` tables.
- **SNP identity policy**: public `snp_identifier` values are exactly `rsid`, `rsid_allele_aware`, `chr_pos`, and `chr_pos_allele_aware`; the default is `chr_pos_allele_aware`.
  Mode names are exact, while column aliases such as `RSID`, `BP`, and `POSITION` apply only to input headers.
  Base modes are fully allele-blind: `rsid` uses only `SNP`, and `chr_pos` uses only `CHR:POS`; allele columns may be preserved as passive data but never affect base-mode identity, duplicate filtering, retention, or drop reasons.
  Allele-free munging is selected through a base identifier mode, not through a separate skip flag.
  Allele-aware modes require usable `A1/A2` for summary statistics, the R²-matrix sidecar, and LD-score directories; they drop missing, invalid/non-SNP and identical allele pairs, multi-allelic base-key clusters, and duplicate effective-key clusters among the rows that pass any SNP restriction.
  A/T and C/G pairs are valid identities (canonical allele sets `A:T`, `C:G`); only signed comparisons between files (`rg`, `query-r2-matrix` signs) treat their strand as unknown.
  Reference-panel identity cleanup runs after `--ref-panel-snps-file` and before validation in every backend: for PLINK panels in `_refpanel.plink_source.select_plink_chromosome()`, which `PlinkRefPanel._load_source()` and `build-r2-matrix` share, so an R² matrix sidecar and `_align_annotations()` receive a reference with unique identities.
  `select_plink_chromosome()` and the R² matrix adapter match the restriction with `_refpanel.plink_source.build_restriction_mask()`, which keeps a row without an allele-aware SNP identity at a restricted position for cleanup to drop and report, so `ldscore --plink-prefix` and an R² matrix built from the same panel and restriction select the same rows; `build-r2-matrix` liftover and duplicate-position filtering can still remove more SNPs.
  `inspect_plink_inputs()` rejects a `.bim` row whose position is missing, non-numeric, nonpositive or not a whole number in every PLINK route, naming the rows, so no route truncates or drops such a position.
  When cleanup removes every row the restriction keeps on a PLINK chromosome, `PlinkRefPanel._load_source()` raises `EmptyReferenceSNPs`, so the chromosome is skipped like one the restriction empties; without a restriction, a chromosome that cleanup empties raises `LDSCInputError`.
  For artifact-like tables, the duplicate policy is always compute the active effective merge key, then drop all rows in duplicate-key clusters.
- **Restrictions and annotations**: restriction files may omit alleles. Allele-free restrictions match by base key and can retain multiple candidate rows before later artifact cleanup. Allele-bearing restrictions, including packaged HM3 restrictions, match by the effective allele-aware key in allele-aware modes. Restriction files are identity-only filters: duplicate restriction keys collapse to one retained key, and non-identity columns such as `CM`, `MAF`, or other metadata are not carried into downstream artifacts. Annotation files may omit alleles even in allele-aware modes because annotations describe genomic membership, not variant alleles; if annotation files include alleles, those alleles participate in allele-aware matching.
- **Artifact compatibility**: Public downstream chaining uses `.annot.gz`, self-describing Parquet munged summary statistics (identity in the footer) with optional `.sumstats.gz` compatibility output, package-built R² matrices, and LD-score directories with root `metadata.json`.
  The forward format rule is parquet for internal artifacts and TSV for science-facing result tables.
  The Parquet payloads of LD-score directories and munged summary statistics use chromosome-aligned row groups where useful, so full-file readers still work while chromosome-specific readers can skip unrelated row groups.
  Legacy `.sumstats(.gz)` inputs are projected automatically onto the regression SNPs in the LD-score directory.
  Legacy reference/weight `.l2.ldscore(.gz)`, `.l2.M`, and `.l2.M_5_50` families require explicit [LD-score conversion](../reference/formats/legacy-ldscores.md#compatibility-boundary); `.w.l2.ldscore(.gz)` suffixes are unsupported.
  Private kernel legacy artifact emitters and regression readers have been removed.
  Allele-free annotations remain supported even in allele-aware workflows, while converted LD-score suites use allele-unaware identity.
- **Output collision handling**: output directories are literal destinations.
  Most workflows create missing directories and reuse existing ones while
  checking fixed output files, including workflow logs, before writing. The
  gene LD-score index builder keeps its live state in the destination's hidden
  `.build/` directory and publishes with root `metadata.json` last. By
  default an existing artifact raises `FileExistsError`; `--overwrite` or
  `overwrite=True` makes replacement explicit without deleting unrelated files
  or cleaning the directory. Sharded workflows may narrow ownership to the
  current shard; `build-r2-matrix` concrete chromosome runs own only their
  chromosome package, while `@` suite runs own the full R² matrix package.
- **Workflow logging**: `ldsc._runtime.workflow_logging.workflow_logging()` attaches file
  handlers to the `LDSC` logger so workflow and kernel records are captured
  together. The log file is the authoritative sink; the CLI console (stderr,
  installed by `run_cli`) stays silent except for errors, and `stdout` is left
  for deliverable command output. Logs are diagnostic files, not scientific outputs,
  so result `output_paths` mappings exclude them. Lifecycle sections use a
  multi-line `Call:` block, separated `Inputs:`/`Outputs:` audit blocks, and an
  explicit elapsed-time footer such as `Elapsed time: 2.0min:12s`. Start/end
  timestamps and elapsed duration are derived from paired entry/exit timepoints
  so the footer reflects the interval covered by the log. A failed run also
  records the full traceback in the failure footer. Gene LD-score index publication keeps
  its open handler outside the replaceable tree, moves the closed successful
  log into the published diagnostics, and treats cleanup after destination
  reload validation as warning-only garbage collection.
- **Dependency split**: base package dependencies cover core pandas/numpy/SciPy workflows and parquet I/O.
  PLINK-backed LD computation requires the `plink` extra (`bitarray`), and cross-build reference-panel output requires the `liftover` extra (`pyliftover`).
  BED/gene interval projection uses NumPy through `ChromosomeProjector` and requires only the base package.
- **Chromosome ordering**: chromosome-sharded inputs are validated and reassembled in stable genomic order by the workflow layer.
- **Testing approach**: tests under `tests/` cover file contracts, workflow behavior, and legacy compatibility expectations.

## Architectural Invariants

- Import from `ldsc`, not from the private `ldsc._*` packages, for supported public use.
- Numerical projection must consume resolved, prepared inputs; user-token resolution belongs to workflows/shared helpers and the existing runtime adapters in `_refpanel.adapters`.
- `_core/columns.py` owns alias families and internal artifact header strictness.
- LD-score computation remains chromosome-wise; regression consumes only the aggregated artifacts.
- Public LD-score output layout is fixed by `LDScoreDirectoryWriter`;
  `partitioned-h2` result tables use fixed compact/full schemas owned by
  `PartitionedH2DirectoryWriter`.
- Package-built R² matrices are the only public R² Parquet format
  for `ldsc ldscore` and `ldsc query-r2-matrix`; external R² Parquet layouts are not
  supported as public inputs.
- New internal or intermediate artifacts should be parquet by default; new
  science-facing result tables should be TSV. Compatibility text outputs must
  be called out as legacy or sidecar exceptions.
- Query annotations are valid only with explicit baseline annotations; the
  synthetic `base` path is for ordinary unpartitioned LD-score generation and
  downstream `h2`/`rg`, not for `partitioned-h2`.
- BED and gene-list queries use one partial-success contract. Only usable query
  columns enter scientific artifacts; skipped inputs remain visible in
  `diagnostics/query_annotation_status.tsv`.
- BED projection uses input intervals as provided unless `padding_bp` /
  `--padding-bp` is set. Padding expands both interval ends in base pairs
  before projection and clips starts at zero; it should not be applied again to
  BED files already expanded upstream. `ldscore` accepts this option only when
  it is performing live BED or gene-list projection; prebuilt-annotation,
  no-query, and indexed modes reject it.
- Every workflow that writes fixed artifacts must precompute expected output
  paths, including its log path, and call the shared output preflight before
  the first write.
- Gene LD-score index construction and assembly preserve Boolean interval union, float64 unbiased-R² accumulation (including negative values), one diagonal contribution, and the split between the broad reference SNPs and the filtered regression SNPs.
