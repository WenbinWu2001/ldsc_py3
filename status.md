# Status

Last updated on: 2026-10-31

Branch `restructure`; checked on 2026-10-31.
Nothing below is in progress; each open item needs the user's go-ahead.

## Completed: output-directory confinement (2026-10-31)

The user set the rule on 2026-10-30: every file a run creates, writes, renames or deletes stays inside its output directory, and an empty `output_dir` never touches user files.
An audit at `2551648` (record outside the repository) found that an empty destination meant the working directory, where an authorized `annotate` overwrite deleted a user's file, and that the gene LD-score index build and the annotation bundle wrote outside the output directory.
The user approved the designs and authorized the work on 2026-10-30; matplotlib's own home-directory cache stays as it is.

- `450ec21`: `output_directory_path()` in `_runtime/paths.py` rejects a missing or blank destination; `ensure_output_directory`, the failure-marker scope, `AnnotationWorkspace`, the `plots/` and liability-scale roots and the remaining direct derivations use it.
  Every `--output-dir`, `--result-dir` and `--h2-result-dir` option rejects a blank value as an argument-parsing error (`cli_help.output_directory_argument`).
- `a7bcb94`: `plot` and `convert-h2-scale` check their result directory before the marker scope opens, so a misspelled result directory is no longer created.
- `b12e304`: the gene LD-score index build keeps its lock, live log, issue tables, history and staged index in the index's hidden `.build/` directory and removes it after success; publication writes root `metadata.json` last, so the loader never reads an interrupted publication.
  The sibling `.<name>.build-state/` and `.<name>.stage-*/` directories, the backup and the parent-directory recovery are gone.
  Accepted cost: a failure during publication leaves no loadable index; a failure before publication keeps the earlier one.
- `c6f4b2c`: the rsID cross-chromosome SQLite join creates its `SNP` indexes before inserting, so it writes no temporary sort file to the system temporary directory.
  On the 1KG 30x BIM files of all 22 autosomes (28.9 million rows per table) the step took 141.4 s instead of 122.8 s, with a 3.35 GB instead of 3.22 GB database and 1.38 GB instead of 1.45 GB peak memory (one run each).
  `EXPLAIN QUERY PLAN` of all 16 queries of the two other stores (`DiskIdentityIndex`, `QuantileIdentityStore`) showed no temporary B-tree, so they are unchanged.
- `fbe9760`: `AnnotationBundle.prepared_in(workspace)`; `LDScoreCalculator.run` with an output directory prepares a saved `run_annotate` result in its own output directory instead of the annotate output directory.
- This commit: `docs/reference/logs-and-outputs.md` (new "Output confinement" section), `docs/reference/gene-ldscore-index.md`, `docs/user-guide/troubleshooting.md`, `docs/user-guide/workflows/gene-list-ld-scores.md`, `docs/dev/architecture.md`, `docs/dev/code-structure.md` and four other reference pages describe the contract; `AGENTS.md` records the rule; `lessons.md` gains "An empty destination is the working directory".
- Evidence: a before/after comparison of exported copies of `2551648` and the final tree matched exactly in the `--help` output of `ldsc` and all 14 commands, every parser's sample parse, and the modules each help path imports; only the 13 destination options gained the blank-value type.
  Of the 29 marker cases only a `plot` run on a missing result directory changed (it no longer creates the directory and a marker), and an empty `munge-sumstats --output-dir` now stops at parsing (exit 2 instead of 1).
  In sandboxes that list every file anywhere, gene-index builds (fresh, overwrite, failed) left `.index.build-state/` beside the index before and nothing outside it now; a chr_pos annotate-then-LD-score run wrote `.ldsc-annotation-<id>/` into the annotate output before and nothing now.
  Published index, annotation and LD-score files were identical in content (gzip files compared decompressed).
- Full pytest on an exported copy of `fbe9760` with this commit's documentation, the workspace `resources/` linked into its parent folder: 2,424 passed, none skipped.
  Each code commit first ran its focused tests (2,356 for `450ec21`, 426 for `a7bcb94`, 1,091 for `b12e304`, 162 for `c6f4b2c`, 1,776 for `fbe9760`).

### Decisions taken during implementation

- `build-gene-ldscore-index` treats index component names beside `.build/` as an interrupted build (rebuild with `--overwrite`); without `.build/` they stay unrecognized and are refused, so a lone user `metadata.json` is never replaced.
- Publication merges `diagnostics/` entry by entry and keeps unrelated files in the index directory, which the old whole-directory swap deleted.
- A cleanup failure of `.build/` after a successful build is a `UserWarning`, because the workflow log is already closed.
- `estimate_partitioned_h2_batch`, indexed `ldscore` and indexed scoring also route their output directory through `output_directory_path()`, though no public entry point reaches them with an empty value.
- `LDScoreCalculator.compute_chromosome` and a direct `bundle.read()` keep preparing a saved result in the annotate output directory, the bundle's own output, as the brief allows.

## Completed: one failure-marker scope per run (2026-10-30)

Architecture review candidate 10.
Re-checked at `a9131ca`, the review's duplicated chromosome regex was already gone (`cd77462`), but the nested guards remained and two marker defects the review had not reported were measured.
The user accepted two rounds of recommendations on 2026-10-30 and authorized the implementation; the question record and evidence are outside the repository.

- `fea62a8`: a failure marker written in the calling process names a workflow log opened in a `serial_calculation` worker; markers of `build-r2-matrix` (CLI and Python), `convert-legacy-ldscores` (CLI) and `quantile-h2` had said "No workflow log was opened" while the log held the traceback.
  The `quantile-h2` and `convert-legacy-ldscores` guards moved into the calling process around separately named worker functions, so a killed worker still leaves a marker.
- `c8d81da`: one marker scope per run.
  `run_cli` no longer opens a scope from raw arguments (`_cli_failure_marker_scope` and `_raw_option_value` are gone); it records its command line with `cli_command()`, and the workflow's scope writes it.
  A scope opened inside an open one neither writes nor removes and must name the same marker, so a failed run writes one marker instead of two (three for `annotate`, two for the Python `run_h2`).
  `munge-sumstats` opens its scope before its option checks, so an invalid option combination leaves a marker from Python as from the CLI; `--infer-only` and console output after a completed run (munging summary, gene-index funnel notice, plot figure close) are outside every scope.
  `--output-dir ""` with `--overwrite` no longer writes `RUN_FAILED.txt` into the working directory.
- `14a8fb7`: `build-r2-matrix` names its marker from the BIM chromosomes, by the rule its log and metadata already used.
  A prefix without a chromosome in its file name, a one-match glob, or a name such as `EUR_5.chr22` had left `RUN_FAILED.txt`, which a successful run of another chromosome in the same output directory then deleted.
  Input checks run in the calling process before the marker opens, so errors found before the chromosomes are known (options, PLINK files, genetic-map headers, keep and SNP restriction files) leave no marker; their diagnostics files still record them.
- This commit: `docs/reference/logs-and-outputs.md`, `docs/reference/input-paths.md` and `docs/dev/architecture.md` describe the contract; `lessons.md` gains "Per-thread state does not cross the worker boundary" and two recurrences; the `run_build_r2_matrix_from_args` docstring says no scientific output has been replaced, since input checks rewrite diagnostics files.
- Evidence: a before/after comparison of exported copies of `a9131ca` and `14a8fb7` (the candidate 5 harness plus 13 marker scenarios) matched exactly in the `--help` output of `ldsc` and all 14 commands, every parser, and the modules each help path imports.
  Of the 29 earlier marker cases only `run_fail_r2_chr` changed, as designed: its missing PLINK input now leaves diagnostics and no marker.
  The scenarios changed only as designed: one marker write per failure, log pointers for `build-r2-matrix` and `convert-legacy-ldscores`, BIM-named `build-r2-matrix` markers that a later chromosome run leaves in place, no marker for a corrupt PLINK input or an empty `--output-dir`; Python keyword errors, the Python legacy marker and the munging option-error marker were unchanged.
- Full pytest on an exported copy of `14a8fb7` with this commit's documentation, the workspace `resources/` linked into its parent folder: 2,374 passed, none skipped.
  Each code commit first ran its focused tests (1,030 for `fea62a8`, 1,359 for `c8d81da`, 1,189 for `14a8fb7`).

### Decisions taken during implementation

- `build-r2-matrix` option errors (for example no LD window) leave no marker, like its input errors, because they too come before the chromosomes are known.
- Tests that replaced `R2MatrixBuilder.run` now replace `_inspect_inputs` and `_run`, and the PLINK retry test fails inside the run through an SNP restriction that matches nothing.

### Follow-ups closed after candidate 10 (2026-10-30)

The user asked for both items left open by candidates 5 and 10 to be fixed.

- `b2f4bf6`: `ldsc -- h2` stops with a usage error (exit 2) that says to put the command first; before, the top-level parser accepted it and the command exited 0 without running anything.
- `9ac5529`: markers from `run_ldscore`, `run_build_r2_matrix` and `build_gene_ldscore_index` name that function, as the other Python entry points already did, instead of the `run_*_from_args(...)` function they call.
  `run_ldscore` and `build_gene_ldscore_index` open their own scope after their argument checks; `run_build_r2_matrix` passes its name to the parsed workflow, whose marker opens only once the chromosomes are known.
- Evidence: focused tests per commit (197 passed for `b2f4bf6`, 1,193 for `9ac5529`); full pytest on an exported copy of `9ac5529` with the workspace `resources/` linked into its parent folder: 2,378 passed, none skipped.

## Completed: one CLI command table (2026-10-30)

Architecture review candidate 5.
Re-checked at `781d75b`, the review's main reason was gone: `8d305cf` had removed `cli.build_parser()`, so tests already parsed with the parsers users run.
The user accepted two rounds of recommendations on 2026-10-30 for the remaining friction and authorized the implementation; the question record and evidence are outside the repository.

- `45c7965`: six parser descriptions that were overwritten before any help was printed, and three repeated `prog` assignments, are gone.
  `add_h2_arguments`, `add_partitioned_h2_arguments`, `add_rg_arguments` and `add_convert_h2_scale_arguments` now set their own description, so every workflow parser owns its description.
- `82cb4a1`: the `_COMMANDS` table in `src/ldsc/cli.py` lists each subcommand once, with its one-line help, its workflow module, and either the module's `build_parser()` and `main(argv)` or the add-arguments and run functions for a parser the CLI builds.
  `main`, `_command_parser`, the top-level help, the failure-marker scope and the internal-error advice read it.
  A command's module is imported by name only when the command is parsed or run, which replaced six `_load_*` helpers and five eager workflow imports; `cli.py` went from 347 to 276 lines.
  `infer-genome-build`'s argument and stdout report moved to `add_inspection_arguments()` and `run_inspection_from_args()` in `src/ldsc/_genome/build_inspection.py`.
  A new test checks that every row names functions its module defines.
- `ddf0366`: an argument-parsing error with `--overwrite`, such as a misspelled or missing option, no longer creates the output directory or writes a `RUN_FAILED.txt` reading `SystemExit: 2`; `overwrite_failure_marker` re-raises `SystemExit` without touching the output scope.
  The explicit help, `--infer-only` and read-only exclusions stay.
  The contract is in `docs/reference/logs-and-outputs.md`, and `lessons.md` records the recurrence.
- This commit: `docs/dev/architecture.md` and `docs/dev/code-structure.md` describe the table and the inspection module's CLI functions.
- Evidence: a before/after comparison of exported copies of `781d75b` and the final tree matched exactly in the `--help` output of `ldsc` and all 14 commands, every parser's settings and actions with one sample parse per command, and the modules each help path imports; `ldsc --help` and `ldsc build-r2-matrix --help` still import neither the regression nor the munging workflow.
  In a 29-case failure-marker matrix, only the 12 argument-error cases changed, as designed; the help, run-failure, `--infer-only`, `infer-genome-build`, unknown-command and no-overwrite cases matched, including exit codes and messages.
- Full pytest on an exported copy of `ddf0366` with this commit's documentation: 2,361 passed, none skipped; each other commit ran its focused tests (385 for `45c7965`, 393 for `82cb4a1`, 607 for `ddf0366`).

### Left open from candidate 5

- Candidate 10: done; see "Completed: one failure-marker scope per run".
- Errors that workflows raise after parsing for invalid option combinations (`LDSCUsageError`) keep leaving a failure marker, as decided in candidate 10, now also from Python for `munge-sumstats`.
- `ldsc -- h2`: fixed in `b2f4bf6`; see "Follow-ups closed after candidate 10".

## Completed: gene-index fixtures through the build writer (2026-10-30)

Architecture review candidate 4.
The user accepted the updated round-1 recommendations on 2026-10-30, after re-checking the evidence against their own commits through `57de7dd`; the question record and evidence are outside the repository.

- `3d41db6`: `AUTOSOMES` in `src/ldsc/_core/chromosomes.py` names the labels 1 to 22 once, replacing the separate copies in `_ldscore/preflight.py` and `_annotation/preflight.py`; `STANDARD_CHROMOSOMES` derives from it.
  `legacy_converter._CHROMS` stays, because it compares integer shard numbers parsed from file names.
- `6ea4824`: `publish_gene_ldscore_index` and `_write_index_artifact`, which no production code called, are gone, and with them the second copy of the payload-size diagnostics.
  `write_gene_ldscore_index` in `tests/gene_index_helpers.py` publishes fixture chromosome records through the build's own transaction, chromosome staging, finalization and commit.
  The tests of overwrite refusal and replacement, post-commit cleanup warnings, rollback and failure-marker reuse now run `run_build_gene_ldscore_index_from_args`; the rollback test fails the install of the replacement index, so the backup restore inside the commit step is now tested.
  The recovery test still calls `_recover_gene_index_publication` on a hand-built interrupted transaction.
- `53c866a`: `INDEX_CHROMOSOMES` in `src/ldsc/_ldscore/gene_index/artifact.py`, defined as `AUTOSOMES`, replaces `_allow_partial_for_tests` in nine functions and the hidden `args._test_chromosomes`.
  The build's chromosome list and the loader's completeness check read it through the module when they run; tests narrow it with `narrow_index_chromosomes`.
  `load_gene_ldscore_index` and its forwarding twin `_load_gene_ldscore_index` are now one function with the same public signature.
- `b3c6089`: `docs/dev/code-structure.md` describes the helpers, and the gene-list diagnostics reference cites the loader where it lives.
- Evidence: an exact before/after comparison of `57de7dd` and `53c866a` over seven index builds matched in `index_id`, every loaded array and table, the root, component and build-diagnostic metadata (without timings), and every indexed `ldscore` output (gzip files compared decompressed).
  The builds covered the four SNP identifier modes and a non-default `--common-maf-min` on the 8-SNP chromosome 1 fixture, a two-chromosome panel, and the 32-SNP 1KG chromosome 22 subset with SNP and kb LD windows.
- Full pytest on an exported copy of `b3c6089`: 2,349 passed, none skipped; every other commit ran its focused tests (478 for `3d41db6`, 367 for `6ea4824`, 402 for `53c866a`).

## Completed: regression trait preparation (2026-10-30)

Architecture review candidate 3.
The user answered two rounds of design questions on 2026-10-30, confirmed the design, and authorized implementation; the question record and evidence are outside the repository.

- `daede5d`: `_regression/outputs.py` names the fitted-model files and `"../../metadata.json"` once.
  `fitted_model_files()` serves `quantile-h2`, which now also checks `artifact_type`, and `published_query_artifacts()` replaces the batch path's re-spelling.
  The file format is unchanged.
- `1717d74`: `_PreparedTraitSumstats` projects a legacy table once per run and keeps its drops.
  Keyed rows are cached per SNP identifier mode and downgrade cleanup, at most two per trait, and `_PreparedLDScores` caches the baseline LD scores the same way.
  `estimate_rg_pairs()` prepares traits once, checks genome builds once, and keys a downgraded pair's rows inside that pair's error handling.
  The CLI no longer projects; `PartitionedH2BatchResult` and `RgResultFamily` gained `legacy_sumstats_drops`, and the `legacy_projection` provenance entry is gone.
- `a75dec9`: rows whose `Z` or `N` is missing or infinite are dropped per trait with a warning and `dropped_nonfinite_statistic_rows` in result metadata, with reason `nonfinite_statistic` in the legacy audit.
  This is the second recorded exception in `AGENTS.md`; see `docs/reference/regression-options.md` section 1.6.
  It replaced a raw `LinAlgError` in single-annotation `h2`, a silent drop logged as a chi-square outlier under a cap, and an `rg` abort.
- `d4a6119` shares `_select_model_columns()` between `h2` and `rg` after their own alignment, `706610b` merges `RGRegressionDataset` into `RegressionDataset`, `854bd38` returns the public `RegressionFit` from `estimate_h2` and `estimate_rg`, and `e0f4f9a` makes `prepare_inputs`, `dataset_from_prepared` and `PreparedRegressionInputs` private.
  Documentation followed in `374cfce` and `3e33290`.
- Order deviation chosen during implementation: the prepared-trait refactor landed before the non-finite fix, so the fix has one home.
- Evidence: an eight-case harness gave byte-identical outputs for every refactor commit against its predecessor; the fix commit added only the new metadata field.
  The cases covered `h2`, `partitioned-h2` and `rg`, a legacy trait, and a mixed-mode `rg` run whose downgrade cleanup drops multi-allelic panel rows.
  A 330,000-SNP, ten-trait `rg` benchmark ran 165.4 s before and 49.0 s after, with peak RSS 1,447 and 1,608 MiB; see `docs/dev/validation/performance.md`.
  Each revision ran once on a heavily swapping laptop, and HapMap3-sized memory was not measured.
- Full pytest at `a5afc19`: 2,301 passed, none skipped (2,291 at `1c472df`; ten new tests).
  The first full run, at `ea6277c`, failed one genome-build test whose stand-in tables lacked `provenance`; `a5afc19` fixed the stand-in.
  Every other commit ran its focused tests.
- Follow-ups the user confirmed on 2026-10-30: `3c3447e` (`estimate_rg(dataset)` replaces the table-taking form and `_fit_rg_dataset`), `a60a1c6` (unestimable `rg` pairs, and jackknife messages that name the failed fit, block, genomic span and unsupported columns) and `e5bb065` (missing or infinite LD scores and regression weights stop the run; a damaged query fails only that query).
  Full pytest at `e5bb065`: 2,310 passed, none skipped.

### Left open from candidate 3

- The user approved "regression fit" (`RegressionFit`) and "prepared trait summary statistics" on 2026-10-30 (`0c78a59`), and the older "fit outcome" wordings were reworded.
- `rg` still stops the whole run when one trait's own heritability fit fails in its per-trait `h2` phase (`RegressionRunner.estimate_rg_pairs()` in `src/ldsc/_regression/runner.py`), even for an unestimable jackknife; only pair-level errors are recorded per pair.
- `tests/regression/test_regression_streaming.py::test_public_directory_reader_does_not_import_regression` hung twice at interpreter exit inside Arrow's thread-pool shutdown, once blocking a full run for about 30 minutes.
  Its child now has a 120-second timeout, so a hang fails the test; the Arrow shutdown deadlock itself is not investigated.

## Completed: one R² matrix chromosome handle (2026-10-30)

The user accepted the round-1 recommendations for architecture review candidate 2 on 2026-10-30 and authorized the implementation, together with a check for other repeated file-name contracts.

- `11ae004` (`refactor(refpanel)!: add R2MatrixChromosome`): `ldscore` and `query-r2-matrix` open each R² matrix chromosome through `R2MatrixChromosome` in `src/ldsc/_refpanel/r2_format.py`.
  Opening checks that both files exist and validates the Parquet footer; the first `sidecar()` call parses, validates and binds the sidecar and keeps it.
  - The sidecar reader is strict: the writer's identity lines, matching the footer; the header `CHR POS SNP A1 A2 CM MAF` in that order; one chromosome; an integer `POS`, a numeric or `NA` `CM`, and `MAF` in [0, 0.5].
    `read_metadata_table`, `_load_full_panel_sidecar`, the alias matching and the layouts the writer never produced are gone.
  - `r2_matrix_file_paths` and `r2_matrix_chromosomes` own the file names for the builder, the adapter, `R2Matrix` and both input checks; `_genome` keeps its own spelling in one place (`R2_PARQUET_GLOB`, `r2_sidecar_path` in `build_policy.py`), and `tests/refpanel/test_r2_matrix_chromosome.py` pins the agreement.
  - Shared errors name the reading command and link its troubleshooting section; `query-r2-matrix` has a new section, and the misplaced `build-r2-matrix` "reference-panel artifact is incompatible" section was folded into the `ldscore` one (`c8e3525`).
  - Breaking: `R2Matrix.open(validate_binding=...)` is gone; `R2MatrixRefPanel.resolve_r2_paths` and `resolve_metadata_paths` became `open_chromosome`; `load_metadata()` returns the sidecar's column order.
- Evidence, kept outside the repository: an exact before/after comparison over 15 cases (three builds and four windows, `ldscore` in four identifier and window settings, `--ref-panel-snps-file` with `--maf-min`, the same-position reorder case, three `query-r2-matrix` modes and the Python API) matched in every `ldscore`, `query-r2-matrix` and `build-r2-matrix` output.
  Only the column order of `load_metadata()` changed, and a sidecar with renamed columns now fails where it used to be read.
  Per chromosome in one `ldscore` run, sidecar parses fell from 4 to 2 (one per process) and Parquet footer reads from 17 to 11; the other 7 footer reads come from `GenomeBuildCheck`, which this work left out of scope, and a worker builds its reference panel three times.
- Repeated file-name contracts found by the audit and fixed: `abc07fd` (LD-score query-batch, log, dropped-SNP and reference-metadata names, built separately by writer, reader and cleanup), `cd77462` (`FAILURE_MARKER`, the shared `args_output_scope` resolver for eight overwrite guards, and one `build-r2-matrix` marker rule for the workflow and the CLI) and `2afb6f1` (the `annotate` log path, the BED query-name rule, and the annotation scratch shard path).
  The candidate-1 fixture outputs matched exactly between `d0ae477` and `2afb6f1`.
- `25b279d`: removed `REFERENCE_METADATA_SPECS` and `REFERENCE_METADATA_SPEC_MAP`, whose only users were the replaced sidecar readers, and corrected the column schema, which listed sidecars among the alias-accepting inputs.
- Full pytest runs on exported copies of `c8e3525` (2,286 passed; the five resource-dependent builder tests skipped there then passed in a rerun with the workspace `resources/` linked, see `lessons.md`) and of `25b279d` (2,291 passed, none skipped); every other commit ran its focused tests.

### Left open from the naming audit

- `diagnostics/input_issues.tsv` is declared only by the LD-score and annotation directory writers; the regression, `query-r2-matrix`, `build-r2-matrix` and munging writers pass it to `InputGate` without declaring it.
- Lower-tier readers that spell higher-tier names (`_genome` and `_runtime` reading LD-score, regression, gene-index and munged file names or `files` keys) stay as they are, as for the R² names.

## Completed: dead-code and redundant-test cleanup (2026-10-29)

The user confirmed the cleanup inventory from the architecture review and the public-name decisions on 2026-10-29 and authorized the work; the inventory, with evidence per item, stayed outside the repository.

- `ab01daf`: tests import private names and pyarrow directly, so a rename or a broken install fails at collection instead of silently skipping up to 331 tests.
- Unused production code, one commit per domain: `a181036` (`_runtime`, including the never-read `suffixes` parameter of the path resolvers), `0c3a5f9` (`_core`), `1046ade` (`_genome`), `33549f4` (`_annotation`), `aa09540` (`_refpanel`, including the warning for R² row groups that only older LDSC3 versions wrote), `63a9ac5` (`_sumstats`) and `29a895b` (`_regression`).
  Two of these reject files from older LDSC3 versions instead of reading them leniently: a munged Parquet footer without `ldsc:snp_identifier` or `ldsc:genome_build` now raises the curated "re-run munge-sumstats" error, and `quantile-h2` requires `retained_ld_columns` in the fitted-model metadata.
- `66e57e3`: the six `try/except ImportError` guards around pyarrow and Matplotlib, both required dependencies, are gone; Matplotlib stays lazily imported.
- `f0211a2`: the three `_ldscore` leftovers.
  The gene-index destination check accepts only a directory holding just `RUN_FAILED.txt`, not the diagnostics layout of failed builds before `e93f403`; `_wrap_legacy_chrom_result` became `_wrap_chrom_result` without its unused parameters; `_normalize_run_args` no longer sets four attributes nothing reads.
  The `config_snapshot is None` checks in `_ldscore/source.py` and `_ldscore/outputs.py` stay: they raise `LDSCInternalError` when the public LD-score writer receives a result built without a snapshot.
- Two fixes, each reproduced before the change: `7527841` (`convert-legacy-ldscores` accepted abbreviated flags) and `d9817f7` (`RefPanelConfig.backend` is now required, `plink` or `r2_matrix`; with the old `"auto"` default a PLINK panel dropped its genetic map with a log line about R² matrices, and `LDScoreCalculator.run` then failed with a raw `TypeError`).
- Public names removed, all breaking: `RegressionDataset.legacy_sumstats_drops` and `RGRegressionDataset.legacy_sumstats_drops` (`02a7b72`; the `h2` audit file test now covers the drop reasons), the positional `SumstatsMunger.run(config, GlobalConfig)` form (`6a4e1c3`), `RefPanelLoader`'s `ref_panel_config` (`d9817f7`), and `cli.build_parser()` with `_copy_actions` (`8d305cf`).
  `cli._command_parser()` now returns the parser each subcommand uses; `cli.main` and the tests (`tests/cli_parsing.py`) share it.
- Tests: `34485e5` (55 tests that only checked a removed name), `4a4bfb8` (39 duplicates), `18ccede` (15 mocked tests that real runs cover), `f646fba` (tests for files from older LDSC3 versions, two rewritten to keep unrelated-file coverage) and `0b615cf` (skips for committed fixtures become failures, three parser tests that passed for the wrong reason fixed, six unused fixture files and the dead helpers removed).
- Verification: the full pytest suite and `python -m unittest discover -s tests -p 'test*.py'` ran on an exported copy of each tested commit, with its own temporary directory, one run at a time; every run passed.
  Final tree (`8d305cf`): pytest 2,289 passed with no skips, unittest 875 OK; before the work (`b2176a7`), pytest collected 2,428 with one skip and unittest ran 989.
  Thirteen commits had full runs; at the user's request the full runs for `a181036`, `66e57e3`, `7527841`, `02a7b72`, `4a4bfb8`, `f646fba` and `0b615cf` were stopped, and `AGENTS.md` now asks for full runs only on the last commit of each group.

### Kept or left open

- `RefPanel.prepare_chromosome(genetic_map=)` stays, contrary to the inventory: `build-gene-ldscore-index` passes its preloaded map (`_ldscore/gene_index/build.py:1677`).
- Not done: importing the six `munge_input` aliases directly in `_sumstats/munger.py` (cosmetic). The R² sidecar reader rewrite landed with architecture candidate 2 (above).
- Tracked elsewhere: `_allow_partial_for_tests` and `publish_gene_ldscore_index`, one LD-score publication path, and index construction through `compute_chromosome()` (below).
- Kept by decision: one rejection test per legacy LDSC flag; the `plot` rejection of an h2 result without regression bins; the converter's lack of `--common-maf-min`; the BED `empty_input` test; unused parameters that mirror legacy LDSC signatures in `_regression/estimators.py`, `_refpanel/plink.py` and `p_to_z`.
- `AGENTS.md` names pytest as the only test runner (`4f71e9f`); this work also ran unittest discovery because its brief asked for it.

## Completed: one batch execution path for direct and indexed LD scores (2026-10-29)

The user confirmed the design in a grilling session on 2026-10-29 (architecture review candidate 1) and authorized its implementation the same day.

- `9059346` (`fix(ldscore): inherit index settings in output`): indexed output records the index's SNP identifier mode and common-SNP MAF threshold through `LoadedGeneLDScoreIndex.global_config` and `count_config`.
  Before, every allele-aware index was recorded as `chr_pos` and every threshold as 0.05, so `h2` refused allele-aware summary statistics against these LD scores without `--allow-identity-downgrade`, and always failed for an `rsid_allele_aware` index.
- `80ee1f7`: dead-code cleanup, including the public `LDScoreCalculator.run(config_snapshot=...)` parameter (a breaking change) and four duplicated worker-count tests.
- `cc1f3e5`: `_runtime.parallelism.chromosome_pool()` creates the spawn pool for direct scoring, indexed scoring and index construction; a crashed worker raises `LDSCInternalError` in all three.
- `9d77843`: `src/ldsc/_ldscore/batch_execution.py` runs the query batches for both LD-score producers, `_DirectProducer` in `calculator.py` and `_IndexedProducer` in `gene_index/batches.py`.
- `224b0da` and `580f9ad`: the user approved the glossary term "LD-score producer", and docs and docstrings use it.
- `e034c8c`: direct scoring no longer repeats the gene-list name check (Gate A) or the gene-option, coordinate-file and query-requires-baseline argument checks, none of which could fire a second time.
- Deviation from the confirmed design: the direct producer stays in `calculator.py`, because moving it to `direct_batches.py` would create an import cycle.
- Evidence: synthetic fixture outputs recorded after `9059346` (direct prebuilt, BED, gene-list and synthetic-base runs, index construction, indexed scoring and the in-memory path) matched exactly after each later code commit; the comparison script and outputs stayed outside the repository.
  The full pytest and unittest suites passed after each code commit.
- `6452102`: the user confirmed the last two audit items on 2026-10-29, and both modes now share them:
  - the gene-list stop sequences in `_ldscore/diagnostics.py` (`fail_gene_list_resolution()`, `fail_unusable_queries()`, `log_query_outcomes()`), with the superset logging: a failed run logs its written diagnostics, gene-list rejections and SNP support, and query statuses in both modes;
  - `_ldscore.outputs.QueryDiagnostics` replaces the five `SimpleNamespace` stand-ins, and `LDScoreDirectoryWriter.write_input_preflight_diagnostics()` replaces the private writer call; what a failed run deletes is unchanged, so a failed input preflight still keeps an earlier run's outputs.
- Candidate 1 is complete.
- Not done, out of scope: removing `_allow_partial_for_tests` and one LD-score publication path.

## Completed: open defects from the architecture review (2026-10-29)

The user asked for each to be diagnosed and given a failing test first; each was reproduced before its fix.

- `f2d1ccd`: indexed scoring now checks gene-list query names against the index's baseline columns and the reserved names before Gate A, through `validate_gene_query_names()` in `_annotation/query_outcomes.py`.
  Before, a query named like a baseline column ran the whole calculation and then failed, and one named `SNP` crashed with a raw `ArrowInvalid`.
- `e43ca63`: `h2` and `partitioned-h2` record the absolute LD-score directory, including a relative `metadata["ldscore_dir"]` from Python callers, and `quantile-h2` rejects a relative recorded path instead of guessing a base directory.
- `7dbf370`: seven config constructors now catch the `LDSCConfigError` their normalizers raise, so users see the message naming the field and its valid values.
- Fixture outputs still matched the commit-1 baseline exactly; the full pytest and unittest suites passed.

## Completed: review decisions of 2026-10-29

- `1b8830d`: each `rg` trait pair applies the `h2` collinearity check (condition number above 1e5) on its own regression rows, with `h2`'s rule that the check needs a stored overlap matrix and two or more LD-score columns; a collinear pair is recorded as a failed pair.
- `4f71e9f`: pytest is the only test runner; `AGENTS.md` (and `CLAUDE.md`, a link to it) and the developer docs no longer run `python -m unittest discover`.
- Missing reference MAF, verified on 2026-10-29: not a reachable defect.
  The R2 matrix sidecar reader rejects a missing or out-of-range MAF, and PLINK genotype QC removes all-missing and monomorphic SNPs, so neither `compute_counts()` (which raises) nor `annotation_statistics()` (which would omit common counts) ever receives missing MAF.
- `13e95ce`: with the user's confirmation, removed the public names that had no production caller and no use in user-facing docs: `assign_legacy_quantiles` and `QuantileAssignment`, `SumstatsTable.subset_to`, `align_to_metadata` and `snp_identifiers`, `RefPanel.filter_to_snps`, and the `summary()` methods of `RefPanel`, `LDScoreResult` and `ChromLDScoreResult` (a breaking change).
  The legacy quantile-boundary tests now run through `run_quantile_h2_from_args`.
  Public names that user-facing docs use were kept, among them `RegressionRunner.estimate_partitioned_h2`, which `docs/reference/index.md` lists.

## Proposed: index construction through the shared LD-score kernel

Split out of architecture review candidate 1 on 2026-10-29; awaiting the user's go-ahead.

- `build-gene-ldscore-index` assembles each chromosome's LD scores and counts in `build_plink_index_chromosome()` (`src/ldsc/_ldscore/gene_index/build.py`) by calling the kernel primitives `ld_scores_from_genotypes()` and `compute_counts()` itself, instead of going through `compute_chromosome()` in `src/ldsc/_ldscore/compute.py`.
- Unifying them changes numerical code: it needs exact golden comparison of index payloads and indexed LD scores, preferably on a full reference panel, and must keep the shared PLINK projection pass for regression-weight LD scores (see `lessons.md`).

## Completed: domain subpackage reorganization of `src/ldsc`

- The reorganization landed on `restructure` on 2026-10-28, from `9236b80` (dead-code cleanup) to the documentation commit `docs: describe the domain layout`.
  [Code structure](docs/dev/code-structure.md) and [architecture](docs/dev/architecture.md) describe the layout, and `tests/package/test_layering.py` enforces its import tiers.
- Follow-up cleanup landed on 2026-10-28 in the `refactor` commits from `cabe537` to `5b3cede`:
  - `window_coordinates` and `get_block_lefts` in `src/ldsc/_refpanel/windows.py` serve both `ldscore` and `build-r2-matrix`; each command keeps its own missing-CM error.
  - `_runtime.artifacts.atomic_write_json` is the one atomic JSON writer, and every file it writes ends with a newline; `h2_scale` uses the shared `atomic_write_dataframe`.
  - The gene-index build no longer migrates `<index>.build/` logs or takes the standalone lock of earlier releases.
  - `_core.identifiers.open_text` is the one text-file opener; it decodes UTF-8 and treats `.gz` in any letter case as gzip.
  - Four dead or test-only helpers were removed, `_genome.liftover.load_hm3_curated_map(path)` became `load_hm3_map_file`, and three modules stopped re-exporting names they do not define.
  - The 24 underscore-prefixed names that other domains imported lost the underscore (`_runtime.preflight._readable` became `check_readable`), so an underscore again means module-local.
- Look-alike helpers deliberately kept, because each carries a domain-specific error message or return type: the four `_require_columns`, the two `_read_json`, and the two `get_compression`.
- Tests moved into per-domain folders on 2026-10-28: `tests/<domain>/` for each domain package (named without the underscore) and `tests/package/` for the root shell; the [test map](docs/dev/code-structure.md#test-map) describes the layout.
- The glossary approved the terms domain package, tier and root shell on 2026-10-28.

## Annotation SNP grids and LD-score annotation coverage

The earlier handoff for this work was written at `4e2fa3c`, before the implementation; this section replaces it.

### Confirmed decisions

These five decisions are settled; do not reopen them.

1. Annotation definitions (gene lists with their catalog and padding, or BED files) are separate from the annotation SNP grid, which selects the SNPs that receive them.
2. `ldsc annotate` takes exactly one annotation SNP grid: a BIM grid (`--ref-snp-bim-sources`) or a temporarily retained annotation grid (`--ref-snp-annot-sources`).
3. Baseline annotations may supply the annotation grid to `annotate`; downstream they remain the conditioning annotations for LD scores and regression, and choosing a BIM grid does not remove them.
4. The annotation SNP grid defines the candidate SNPs, reference-panel filtering and explicit restrictions select the reference SNPs, and regression SNPs are a separate subset; R² matrices take SNP identities and coordinates from their metadata sidecars.
5. Annotation coverage is validated during `ldscore` reference preparation, before LD-score computation, without zero-filling annotation values or silently shrinking the reference SNPs.

### Implemented

- Decisions 1 to 3 landed in `3ef653a` (`feat(annotate): add BIM grids and aligned exports`, 2026-10-20).
  The CLI uses the mutually exclusive `--ref-snp-bim-sources` and `--ref-snp-annot-sources`; the BIM grid reads `.bim`/`.bim.gz` only and creates no baseline values.
  The contract is [docs/reference/annotate.md](docs/reference/annotate.md); later commits `12cb3fd`, `1b928d2` and `e00fca9` added BED-query and reference-panel failures and that reference page.
- Carryover audit findings now resolved:
  successful partial runs print stderr notices (`emit_resolved_only_notice()` and `emit_gene_gate_b_notice()` in `src/ldsc/_annotation/query_outcomes.py`; contract in `annotate.md`, "Command-line notices");
  gene-list naming collisions are reported together with rejected rows (`naming_error` in `_gate_a_error()`);
  [data-flow.md](docs/dev/data-flow.md) describes both BED and gene-list routes.

### Implemented: decision 5 in `ldscore` (2026-10-28)

The user confirmed the design in a grilling session on 2026-10-28 and authorized its implementation the same day.
The contract is in [reference SNPs on the annotation SNP grid](docs/reference/ldscore-snp-universes.md#reference-snps-on-the-annotation-snp-grid).

- The annotation files define the annotation SNP grid; the reference SNPs are the grid SNPs among the reference-panel SNPs that pass the panel's filters.
  Panel SNPs off the grid are excluded without failing: the user's recorded exception to the AGENTS.md rule that unrequested, analysis-changing filtering must fail, because legacy baseline-LD annotations exist only for 1000 Genomes Phase 3 SNPs.
  With the user's approval, that AGENTS.md rule now carries a one-line pointer to the exception.
  No CLI option was added, and `ldscore` does not read the grid kind recorded by `annotate`.
- The exclusion is never silent: `diagnostics/reference_snp_funnel.tsv` (two funnels per chromosome, fixed rows per backend, `NA` for steps that did not run), one log WARNING and a stderr notice after successful CLI runs, and `LDScoreResult.reference_snp_funnel` for in-memory runs.
  PLINK counts the grid step before genotype QC; counts come from batch 1 only.
- A chromosome emptied at the grid step now fails the run, naming every such chromosome (in Gate B for gene-list runs); one emptied by `--ref-panel-snps-file` or genotype QC is skipped with a logged WARNING instead of a Python `UserWarning`.
- `snp_universe_policy.ld_reference_universe.selection` was removed from `metadata.json`; the reference-panel source is in the log's `Inputs:` section.
- Gene LD-score indexes are format version 4: each chromosome stores `reference_snp_funnel.parquet`, which replaces `pre_qc_rows` and `annotation_intersection_removed` (the latter counted every `.bim` row, not the chromosome's); indexed scoring writes the inherited funnel and repeats the WARNING.
- Baseline and prebuilt query files still need identical rows in identical order; that row-mismatch error now explains the shared annotation SNP grid and links a new troubleshooting section.
- The earlier decision-5 example ("a BIM grid A, B, C with a baseline covering A and B") already failed through row alignment before this work, so it was not the missing check; the tests in `tests/ldscore/test_reference_snp_funnel.py` cover the grid-versus-panel behavior instead.

### Resolved follow-ups (2026-10-28)

- `annotate`'s gene-list preflight error (`_gate_a_error()`) now reuses `gene_list_gate_a_message()`, the `ldscore` wording: per-source rejected/submitted counts, up to ten examples with the omitted count, a suggested audit filter and the `--overwrite` rerun, followed by any naming collision and the troubleshooting link.
  A message with no rejected rows no longer prints an empty `()`.
  Each example names its gene-list file, and the count of examples not shown follows the examples.
- rsID-family runs on an R² matrix no longer print "genome_build is set but will be ignored": `_artifact_global_config()` drops the build that the panel records instead of passing it to `GlobalConfig`.
  The warning remains for a build the caller supplies.

## Reading map

- [annotate contract](docs/reference/annotate.md) and [LD-score SNP universes](docs/reference/ldscore-snp-universes.md).
- `src/ldsc/_refpanel/adapters.py` (`_align_annotations()`, `PlinkRefPanel._load_genotypes()`, `RefPanel.prepare_chromosome()`), `src/ldsc/_refpanel/plink.py` (per-SNP genotype QC), `src/ldsc/_ldscore/preflight.py` (chromosome-coverage gate) and `src/ldsc/_ldscore/calculator.py` (`_snp_universe_policy()`, batch-1 counts, chromosome skips).
- `src/ldsc/_annotation/annotate.py` (`run_annotate()`, `_gate_a_error()`), `src/ldsc/_annotation/sources.py` (`_aligned_metadata()`) and `src/ldsc/_ldscore/gene_index/build.py` (per-chromosome evidence).
- [Lessons](lessons.md), particularly reference preparation, physical BIM/BED row mapping, and output ownership.
