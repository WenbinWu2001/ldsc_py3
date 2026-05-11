# `munge-sumstats` Liftover Implementation Plan

> **For agentic workers:** Use `superpowers:subagent-driven-development` or
> `superpowers:executing-plans` to implement this plan task by task.

## Goal

Add opt-in genome-build liftover to `ldsc munge-sumstats` for
`snp_identifier = chr_pos` workflows. Liftover updates `CHR/POS`, preserves
`SNP` exactly as a label, writes auditable metadata, and keeps downstream
compatibility on the existing `GlobalConfig` build snapshot.

Also audit package-wide `chr_pos` identity behavior so matching uses `CHR/POS`
everywhere identity is intended. Assume `SNP` contains rsIDs even when it does
not.

## Locked Decisions

- Liftover is only valid for `snp_identifier = chr_pos`.
- Any liftover request in `rsid` mode errors before input IO.
- `SNP` is always treated as a label/rsID field, even if it contains
  coordinate-looking strings.
- In `chr_pos` mode, identity and matching are always based on `CHR` and `POS`,
  never `SNP`.
- Liftover updates `CHR/POS` and never rewrites `SNP`.
- Source build remains `GlobalConfig.genome_build`; do not add
  `MungeConfig.genome_build`.
- Add only `target_genome_build`, `liftover_chain_file`, and
  `use_hm3_quick_liftover` to `MungeConfig`.
- `--sumstats-snps-file` is interpreted in the source build because it runs
  before liftover.
- Method flag without target build is an error.
- Target differs from source without method is an error.
- Target equals source without method is a no-op.
- Target equals source with method is an error.
- Chain-file and HM3 quick methods are mutually exclusive.
- HM3 quick uses packaged `src/ldsc/data/hm3_curated_map.tsv.gz`.
- Do not add `utils/build_hm3_dual_build_map.py`.
- Runtime HM3 map loader hard-errors on duplicate `(CHR, hg19_POS)` or
  `(CHR, hg38_POS)`.
- Missing/unmapped rows are dropped and reported for both HM3 quick and chain
  liftover.
- Chain cross-chromosome hits are dropped and counted separately.
- Chain mapping keeps the first same-chromosome hit; mixed cross/same hit lists
  keep the same-chromosome hit.
- Duplicate source `CHR/POS` groups are dropped before mapping when a sumstats
  liftover request enters the map stage.
- Duplicate target `CHR/POS` groups are removed entirely: if multiple original
  rows map to the same target coordinate, none are retained.
- If liftover drops all rows, error instead of writing an empty artifact.
- Log count summaries at normal verbosity. Example rows, if emitted, appear
  only at `DEBUG`.
- Sidecar metadata stays thin: schema marker, trait label, and full
  `config_snapshot` only.
- Coordinate provenance, liftover reports, HM3 provenance, row counts, and output
  bookkeeping are written to the log instead of the sidecar.
- Metadata method for HM3 quick is `hm3_curated`; CLI flag remains
  `--use-hm3-quick-liftover`.
- Reference-panel PLINK builds keep the current source-build plus optional
  opposite-build emission contract. Do not add public target/method config
  fields for that workflow.
- Reference-panel matching chain liftover is invalid in `rsid` mode.
- Reference-panel `duplicate_position_policy` / `--duplicate-position-policy`
  was removed by the 2026-05-10 harmonization. `drop-all` is the only
  coordinate duplicate behavior.
- Reference-panel coordinate duplicate handling applies only in `chr_pos` mode.
  Source-only `rsid` builds log once that duplicate-position policy is not
  applicable.
- Reference-panel dropped-SNP sidecars remain under `dropped_snps/`, are
  always written for processed chromosomes, and contain the four
  ref-panel-applicable liftover-stage reasons. Runtime metadata TSV/parquet
  schemas stay unchanged.

## Harmonization Addendum

The liftover module is an internal service/helper layer, not a new public API.
Both sumstats munging and reference-panel building use it for chain-file mapping,
workflow-aware chain/dependency errors, drop reports, duplicate-coordinate
detection, and readable examples. Workflow contracts stay separate:

- Sumstats uses explicit source/target/method validation from `GlobalConfig` and
  `MungeConfig`.
- Reference-panel building uses explicit or inferred PLINK source build plus an
  optional matching chain to emit the opposite build.
- HM3 quick liftover remains sumstats-only.
- Provenance details live in workflow `.log` files. Sumstats metadata sidecars
  are compatibility-only (`format`, `trait_name`, `config_snapshot`), and
  existing metadata sidecars that lack `config_snapshot` are invalid rather
  than migrated. Row-level liftover drops are audited in the always-written
  `dropped_snps/dropped.tsv.gz` sidecar.

## Pre-Flight

- [ ] Confirm current worktree state and keep unrelated user changes intact.
- [ ] Confirm the packaged map exists:
  `src/ldsc/data/hm3_curated_map.tsv.gz`.
- [ ] Run the existing suite before implementation:

```bash
pytest -q
```

## Task 1 - Shared Liftover Module

**Files:** `src/ldsc/_kernel/liftover.py`,
`src/ldsc/_kernel/ref_panel_builder.py`,
`tests/_kernel/test_liftover_module.py`

- [ ] Move `LiftOverTranslator` and `LiftOverMappingResult` from
  `_kernel/ref_panel_builder.py` into new `_kernel/liftover.py`.
- [ ] Add `LiftoverError(RuntimeError)` for workflow-level liftover failures.
- [ ] Keep reference-panel behavior unchanged by importing/re-exporting the
  relocated symbols.
- [ ] Generalize missing-chain and missing-optional-dependency messages so they
  are not build-ref-panel-specific.
- [ ] Preserve explicit chain-file semantics; no automatic fallback to HM3.
- [ ] Test module import, identity translation, missing chain path error,
  missing `pyliftover` dependency error, and legacy import compatibility.

## Task 2 - HM3 Curated Map Loader And Lifter

**Files:** `src/ldsc/_kernel/liftover.py`,
`src/ldsc/data/hm3_curated_map.tsv.gz`,
`tests/fixtures/hm3_curated_map.tsv.gz`,
`tests/_kernel/test_liftover_module.py`

- [ ] Implement `load_hm3_curated_map(path)` with cached reads.
- [ ] Loader resolves the logical chromosome, hg19 position, hg38 position, and
  SNP-label columns using the package's existing column-alias inference helpers;
  preferred canonical names are `CHR`, `hg19_POS`, `hg38_POS`, and `SNP`.
- [ ] Ignore extra columns after resolving the four logical columns.
- [ ] Normalize `CHR` and integer position columns.
- [ ] Hard-error if the map has duplicate `(CHR, hg19_POS)` or duplicate
  `(CHR, hg38_POS)` keys.
- [ ] Implement `Hm3DualBuildLifter(source_build, target_build, map_path)`.
- [ ] Lookup is coordinate-only:
  - hg19 to hg38 uses `(CHR, hg19_POS)` and writes `hg38_POS`.
  - hg38 to hg19 uses `(CHR, hg38_POS)` and writes `hg19_POS`.
- [ ] Do not use map `SNP` as a lookup key; retain it only for diagnostics.
- [ ] Return lifted frame plus unmapped row indices/counts.
- [ ] Test success in both directions, identity no-op, unmapped rows, invalid
  builds, column-alias inference, duplicate map keys, and coordinate-only
  behavior with mismatching `SNP` labels.
- [ ] Ensure `setup.py` package data already includes `data/*.tsv.gz`; do not
  add a map-generation script.

## Task 3 - Config And CLI Validation

**Files:** `src/ldsc/config.py`, `src/ldsc/cli.py`,
`src/ldsc/sumstats_munger.py`, `tests/test_munge_config_liftover.py`

- [ ] Add `target_genome_build`, `liftover_chain_file`, and
  `use_hm3_quick_liftover` to `MungeConfig`.
- [ ] Normalize target build aliases to `hg19`/`hg38`; reject target `auto`.
- [ ] Normalize optional chain-file path.
- [ ] Reject simultaneous chain-file and HM3 quick method flags.
- [ ] Do not change `GlobalConfig` rsID-mode normalization.
- [ ] Add CLI flags:
  - `--target-genome-build`
  - `--liftover-chain-file`
  - `--use-hm3-quick-liftover`
- [ ] In workflow/API validation before input IO:
  - reject any liftover request when `global_config.snp_identifier != "chr_pos"`;
  - reject method flag without target build;
  - reject target differs from source without method after source build resolves;
  - reject target equals source with method;
  - allow target equals source without method as no-op.
- [ ] Apply identical validation to CLI and programmatic `SumstatsMunger.run()`.

## Task 4 - Kernel Liftover Stage

**Files:** `src/ldsc/_kernel/sumstats_munger.py`,
`tests/_kernel/test_sumstats_munger_liftover_stage.py`

- [ ] Add request/plan/report structures for source build, target build, method,
  chain file, HM3 map file, and count fields, including
  `n_missing_chr_pos_dropped`.
- [ ] Insert liftover after `_finalize_coordinate_columns()` and
  `filter_sumstats_snps()`, before output writing.
- [ ] Keep `--sumstats-snps-file` filtering in source-build coordinates.
- [ ] If liftover is requested and any retained row lacks complete `CHR/POS`,
  drop that row before mapping and log the count/examples.
- [ ] Keep malformed non-missing `CHR/POS` values as validation errors.
- [ ] Apply HM3 or chain mapping and update only `CHR/POS`; never rewrite `SNP`.
- [ ] Drop every row in duplicated source `CHR/POS` groups before mapping.
- [ ] Drop and count unmapped rows for both methods.
- [ ] Drop and separately count chain cross-chromosome hits.
- [ ] Drop every row in duplicated target `CHR/POS` groups; do not keep the
  first row.
- [ ] If all rows are dropped by liftover, raise an error.
- [ ] Log warnings or info records with up to 5 example rows for
  missing-coordinate, unmapped, cross-chromosome, and duplicate-target-coordinate
  drops.
- [ ] On applied liftover, set logged coordinate provenance `genome_build` to
  the target build and attach the liftover report for logging.

## Task 5 - Thin Sidecar Metadata And Logging

**Files:** `src/ldsc/sumstats_munger.py`,
`tests/test_sumstats_munger_liftover.py`

- [ ] Write only thin sidecar fields: `format`, `trait_name`, and
  `config_snapshot`.
- [ ] Keep `config_snapshot` as the downstream compatibility block with
  `snp_identifier`, `genome_build`, `log_level`, and
  `fail_on_missing_metadata`.
- [ ] Log coordinate provenance, including coordinate source columns, output
  coordinate build, inference status, coordinate basis, missing-coordinate
  counts, and optional build-inference details.
- [ ] Log the liftover report with method, source/target build, path fields when
  present, `n_input`, `n_lifted`, `n_dropped`,
  `n_missing_chr_pos_dropped`, `n_unmapped`, `n_cross_chrom`, and
  `n_duplicate_source_dropped`, `n_duplicate_target_dropped`.
- [ ] Log HM3 provenance when HM3 quick liftover is selected.
- [ ] In no-op liftover reports, use `source_build = null` and
  `target_build = null` because the source/target fields are not applicable.

## Task 6 - Package-Wide chr_pos Identity Audit

**Files:** start with `src/ldsc/sumstats_munger.py`,
`src/ldsc/_row_alignment.py`, `src/ldsc/regression_runner.py`, then audit
annotation, reference-panel, and LD-score modules for identity-sensitive merge
or filter sites.

- [ ] Update `SumstatsTable.snp_identifiers()` to use `CHR/POS` in `chr_pos`
  mode and `SNP` in `rsid` mode.
- [ ] Update `SumstatsTable.subset_to()` to use mode-aware identity keys.
- [ ] Update `SumstatsTable.align_to_metadata()` to align by `CHR/POS` in
  `chr_pos` mode and by `SNP` in `rsid` mode.
- [ ] Make `_row_alignment.assert_same_snp_rows()` identifier-mode aware. In
  `chr_pos` mode, compare `CHR/POS` only and allow `SNP` label mismatches.
- [ ] Add tests proving regression `h2` merges sumstats and LD scores by
  `CHR/POS` in `chr_pos` mode even when `SNP` labels differ.
- [ ] Add tests proving `rg` merges two sumstats and LD scores by `CHR/POS` in
  `chr_pos` mode even when `SNP` labels differ.
- [ ] Audit annotation, reference panel, LD-score, and regression code for
  unconditional `SNP` identity use. Replace identity operations with
  `build_snp_id_series(..., mode)` or temporary `CHR/POS` keys as applicable.
- [ ] Preserve explicitly label/allele-based `SNP` operations where the contract
  is rsID/label matching, such as allele merge inputs.

## Task 7 - End-To-End Liftover Tests

**Files:** `tests/test_sumstats_munger_liftover.py`,
`tests/test_regression_build_mismatch.py`

- [ ] Config and CLI flag validation.
- [ ] HM3 curated map load and duplicate-key rejection.
- [ ] HM3 curated map column-alias inference.
- [ ] HM3 coordinate-only lookup with mismatching `SNP` labels.
- [ ] HM3 hg19 to hg38 and hg38 to hg19.
- [ ] Chain-file path with mapped, unmapped, and cross-chromosome rows.
- [ ] Row dropping and count fields for both methods.
- [ ] Duplicate target-coordinate group removal.
- [ ] All-dropped error.
- [ ] Missing coordinate error.
- [ ] `SNP` unchanged after liftover.
- [ ] No-op metadata schema.
- [ ] Applied metadata schema.
- [ ] Thin sidecar metadata and detailed log provenance.
- [ ] Lifted hg38 `chr_pos` sumstats reject hg19 LD scores through existing
  `ConfigMismatchError`.
- [ ] Lifted hg38 `chr_pos` sumstats proceed with hg38 LD scores.

## Task 8 - User-Facing Docs

**Files:** `docs/current/data_flow.md`, `docs/current/architecture.md`,
`tutorials/munge_sumstats.md` or equivalent

- [ ] Document the new stage after SNP-file filtering.
- [ ] Document that liftover is `chr_pos`-only.
- [ ] Document that `SNP` is a label and `chr_pos` identity uses `CHR/POS`.
- [ ] Document source-build filtering for `--sumstats-snps-file`.
- [ ] Document row-dropping behavior and sidecar counts.
- [ ] Add chain-file example.
- [ ] Add HM3 quick example using `--use-hm3-quick-liftover`.
- [ ] Document the rsID-mode error and rationale.

## Final Verification

Run focused tests:

```bash
pytest tests/_kernel/test_liftover_module.py \
       tests/_kernel/test_sumstats_munger_liftover_stage.py \
       tests/test_munge_config_liftover.py \
       tests/test_sumstats_munger_liftover.py \
       tests/test_regression_build_mismatch.py -v
```

Run package-wide tests:

```bash
pytest -q
```

Run CLI smoke check:

```bash
ldsc munge-sumstats --help | grep -E "target-genome-build|liftover-chain-file|use-hm3-quick-liftover"
```

Expected final behavior:

- Liftover works for `chr_pos` hg19 to hg38 and hg38 to hg19.
- Liftover flags in `rsid` mode fail early.
- `SNP` is unchanged by liftover.
- `chr_pos` matching uses `CHR/POS` package-wide.
- Output sidecars contain only the thin compatibility metadata.
- Existing regression build compatibility rejects mismatched concrete builds.
