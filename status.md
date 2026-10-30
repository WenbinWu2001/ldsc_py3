# Status

Last updated on: 2026-10-30

Branch `restructure`; checked on 2026-10-30.
Nothing below is in progress; each open item needs the user's go-ahead.

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

### Left open from candidate 3

- The user approved "regression fit" (`RegressionFit`) and "prepared trait summary statistics" on 2026-10-30 (`0c78a59`), and the older "fit outcome" wordings were reworded.
- `rg` pair datasets keep no panel coordinates and `_fit_rg_dataset()` adds no block-jackknife context, unlike `h2`; on 2026-10-30 the user decided not to pursue this.
- Missing or infinite LD scores and regression weights are not covered by the new drop.
  `docs/dev/validation/legacy-issues-audit/bugs.md` (#174) still describes the earlier state.
- `tests/regression/test_regression_streaming.py::test_public_directory_reader_does_not_import_regression` hung once at interpreter exit inside Arrow's thread-pool shutdown; it passed when rerun.

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

- The gene-index build writes `diagnostics/build-gene-ldscore-index.json` with its own copy of the `payload_bytes` logic in `artifact.py`; this belongs to the index-publication candidate.
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
