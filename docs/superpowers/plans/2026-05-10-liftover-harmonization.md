# Liftover Harmonization Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use
> `superpowers:subagent-driven-development` (recommended) or
> `superpowers:executing-plans` to implement this plan task-by-task. Steps use
> checkbox (`- [ ]`) syntax for tracking. Spec:
> [`docs/superpowers/specs/2026-05-10-liftover-harmonization.md`](../specs/2026-05-10-liftover-harmonization.md).

**Goal.** Harmonize liftover collision handling across `munge-sumstats` and
`build-ref-panel`:

1. Remove `--duplicate-position-policy` (CLI + config field). `drop-all` is
   the only behavior.
2. Both workflows write a `dropped_snps/` sidecar covering the
   workflow-applicable drop reasons. Sumstats covers all five:
   `missing_coordinate`, `source_duplicate`, `unmapped_liftover`,
   `cross_chromosome_liftover`, `target_collision`. Ref-panel covers four
   (everything except `missing_coordinate`, which is unreachable because
   PLINK BIM `BP` is structurally non-null — see spec R3 per-workflow
   coverage table). The vocabulary is shared so consumers can write one
   downstream parser; not every reason fires in every workflow.
3. Per-row example log lines move from `INFO` to `DEBUG`. `INFO` lines stay
   count-only with a sidecar-path pointer.

**Architecture.** Algorithm changes are minimal — the shared
`_kernel/liftover.py` already exposes `duplicate_coordinate_drop_result`,
`liftover_drop_report`, and `mapping_reason_masks`. The work is a workflow-layer
rewiring:

- `_kernel/liftover.py::log_liftover_drop_report` becomes count-only at `INFO`;
  examples emit on a separate `DEBUG` record. A new optional `sidecar_path`
  argument lets callers append a "see ..." pointer; callers that don't know the
  path yet (e.g., ref-panel's `_resolve_mappable_snp_positions`) omit it and
  rely on a separate summary log emitted by the workflow layer once the path
  is known.
- `apply_sumstats_liftover` returns a 3-tuple `(lifted_frame, report_dict,
  drop_frame)`. A dataclass return would force every existing call site
  (1 production + 6 tests) to migrate to attribute access without buying
  anything; the drop_frame does not even ride to the workflow as a return value
  (see transport note below). All 7 unpack sites are migrated to the new
  arity in Task 2.
- The drop_frame travels from kernel to workflow via the existing
  `args._coordinate_metadata` side-channel that already carries
  `liftover_report` between `_apply_liftover_if_requested`
  (`_kernel/sumstats_munger.py:961`) and `SumstatsMunger.run`
  (`sumstats_munger.py:420`). A new `coordinate_metadata['liftover_drop_frame']`
  key holds the DataFrame; the workflow layer pops it and writes
  `dropped_snps/dropped.tsv.gz`.
- `ref_panel_builder._build_chromosome` extends `dropped_frames` to include
  `unmapped_liftover` and `cross_chromosome_liftover` rows produced inside
  `_resolve_mappable_snp_positions`. The reference-panel sidecar will not
  contain `missing_coordinate` rows because PLINK `.bim` `BP` is structurally
  non-null — `chrom_metadata["POS"].astype(int)`
  (`ref_panel_builder.py:455`) and the `to_numpy(dtype=int)` cast at
  line 630 raise on malformed values long before sidecar logic runs. The
  union vocabulary is shared; not every reason fires in every workflow.
  Sidecar column names are unchanged; the dtype contract widens to the
  unified nullable schema described below.
- `ReferencePanelBuildConfig.duplicate_position_policy` is deleted; the
  internal `_resolve_unique_snp_set(..., policy=...)` parameter and the `error`
  branches are removed.

**Sidecar schema (revised — both workflows).** Use pandas nullable types so
sumstats `missing_coordinate` rows (which may have NA `CHR`, NA `POS`, or
both) are representable:

| Column | Pandas dtype | Notes |
| --- | --- | --- |
| `CHR` | `string` | Nullable. NA only for sumstats `missing_coordinate` rows missing `CHR`. |
| `SNP` | `string` | Source label, preserved unchanged. Always present. |
| `source_pos` | `Int64` | Nullable. NA for sumstats `missing_coordinate` rows missing `POS`. |
| `target_pos` | `Int64` | Nullable. NA for `missing_coordinate`, `source_duplicate`, `unmapped_liftover`, `cross_chromosome_liftover`. |
| `reason` | `string` | One of the controlled vocabulary. |

The ref-panel sidecar keeps the same column names, and `source_pos` remains
structurally non-null because BIM `BP` is non-null. The shared dtype contract
still widens both position columns to nullable `Int64`: `target_pos` is empty
for `source_duplicate`, `unmapped_liftover`, and
`cross_chromosome_liftover` rows, and clean runs may produce header-only
sidecars. The canonical TSV is gzip-compressed; the writer passes `na_rep=""`
so nullable-NA cells render as empty strings (matching the spec R2 contract)
which round-trips on read with the explicit `dtype="string"`/`"Int64"` map.

**Always-write contract for the new sidecar (both workflows).** The
sidecar is class-1: always-owned and always-produced. When the current
run drops zero rows, write a header-only gzipped TSV. There is **no
conditional-produced wiring** and **no `Path.unlink()` for this
artifact** in either workflow. A re-run simply overwrites the prior file
in place. This is the simplification confirmed on 2026-05-11.

For `MungeRunSummary.output_paths`: include the sidecar as
`dropped_snps_tsv_gz` unconditionally (always present in the mapping
because the file is always written). Aligns with `metadata_json` being
unconditionally present and matches the file's class-1 status.

**Class-2 (conditionally-produced) artifacts are out of scope.** Sumstats
keeps its existing `preflight_output_artifact_family` +
`remove_output_artifacts` machinery for `output_format` siblings
(sumstats_munger.py:390, 474, 527, 544). Ref-panel keeps
`ensure_output_paths_available` with no auto-cleanup for per-chromosome
and per-build artifacts (ref_panel_builder.py:233). Neither workflow's
class-2 contract changes in this plan. Spec section R7 documents the
asymmetry and the silent-mix failure mode users must manage themselves
for ref-panel.

**Required in-code documentation when implementing.** Each task that
touches preflight wiring must add comments that name and contrast
`owned_paths` (the workflow's stable territory — every file it could
write under any flag combination) and `produced_paths` (only what the
current invocation will actually write). The detailed contract lives in
the function-level docstring at
[path_resolution.py:505](../../src/ldsc/path_resolution.py:505); call-site
comments should be one or two short lines pointing at *why* a particular
path is in one list. For the new sidecar this is straightforward — it
appears in both lists unconditionally because it is always written.
Document the always-written contract at the writer call site so the next
implementer does not "fix" the empty-write as a no-op. This is a
documentation requirement, not a separate task — verify each modified
call site has these comments before marking the relevant task complete.

**Tech Stack.** Python 3.11, pandas, numpy, pytest (the project still accepts
unittest classes; new tests follow the surrounding file's style).

---

## File Map

| File | Change |
|---|---|
| `src/ldsc/_kernel/liftover.py` | `log_liftover_drop_report` → INFO count-only + optional sidecar pointer; examples on a separate `DEBUG` line. `apply_sumstats_liftover` returns per-reason drop frames in addition to the counts dict. |
| `src/ldsc/sumstats_munger.py` | Collect per-reason drop frames from kernel result via `args._coordinate_metadata['liftover_drop_frame']`; **always** write `dropped_snps/dropped.tsv.gz` (header-only when no drops); list it unconditionally in both `produced_paths` and `owned_paths`; surface its path in run summary unconditionally. No `remove_output_artifacts` for this artifact. |
| `src/ldsc/ref_panel_builder.py` | Append `unmapped_liftover`/`cross_chromosome_liftover` rows to `dropped_frames` in `_build_chromosome` (NB: `missing_coordinate` does not apply — PLINK BIM `BP` is structurally non-null). **Always** call `_write_dropped_sidecar` per chromosome processed (header-only when no drops). Remove the conditional write at lines 482-518. Remove `--duplicate-position-policy` CLI flag, the `policy=` plumbing, and the two `error` branches in `_resolve_unique_snp_set`. Update `_write_dropped_sidecar` to handle empty input and summarize all reason categories. (NB: the stale-class-2 `WARNING` helper is split to a follow-up plan; see [`2026-05-11-ref-panel-stale-class2-warning.md`](./2026-05-11-ref-panel-stale-class2-warning.md).) |
| `src/ldsc/config.py` | Delete `ReferencePanelBuildConfig.duplicate_position_policy` field, its docstring entry, and its `__post_init__` validation block. |
| `tests/test_ref_panel_builder.py` | Delete `ReferencePanelBuildConfigDuplicatePolicyTest`, the `error`-policy integration test, the parser-default test, and the rsid "ignored once per run" test (rewrite as "no policy field exists"). Add tests for: sidecar contains all four ref-panel-applicable reasons; header-only sidecar on clean chromosomes; sidecar overwritten in place across re-runs; no-overwrite collision still blocks. (Stale-class-2 warning tests live in the follow-up plan.) |
| `tests/test_sumstats_liftover.py` | Add tests for new `dropped_snps/dropped.tsv.gz` artifact: presence, schema/dtypes, all five reasons, header-only-when-no-drops, no-overwrite collision, in-place overwrite under `overwrite=True`, INFO log points to sidecar path. |
| `tests/test_logging_refactor.py` (or new `tests/test_liftover_logging.py`) | Assert that `log_liftover_drop_report` emits one `INFO` count line and one `DEBUG` examples line; `INFO` includes the sidecar pointer when supplied. |
| `docs/current/liftover-harmonization-decisions.md` | Replace the `duplicate_position_policy` paragraph; document the unified sidecar schema/reason vocabulary; document the sumstats sidecar. |
| `docs/current/data_flow.md`, `docs/current/architecture.md` | Add the new sumstats artifact to the munge stage diagram/text. |
| `tutorials/munge_sumstats.md` (or equivalent) | One paragraph: "auditing dropped SNPs" → sidecar pointer. |
| `docs/superpowers/specs/2026-05-09-munge-sumstats-liftover.md` | Append a `## 2026-05-10 Update` addendum noting the policy collapse and the sumstats sidecar. |
| `docs/superpowers/plans/2026-04-30-duplicate-position-policy.md` | Append note: original plan superseded by `2026-05-10-liftover-harmonization.md`; `duplicate_position_policy` removed. |

---

## Pre-Flight

- [ ] Confirm worktree is clean except for the in-progress harmonization
      changes.
- [ ] Re-read the spec at
      `docs/superpowers/specs/2026-05-10-liftover-harmonization.md`.
- [ ] Run the existing suite as a baseline:

```bash
source /Users/wenbinwu/miniforge3/etc/profile.d/conda.sh && conda activate ldsc3
pytest -q tests/test_ref_panel_builder.py tests/test_sumstats_liftover.py
```

Capture the green baseline so regressions during refactoring are obvious.

---

## Task 1 — Quiet kernel logging; add sidecar pointer

**Files:** `src/ldsc/_kernel/liftover.py`, new
`tests/test_liftover_logging.py`.

The current `log_liftover_drop_report`
(`src/ldsc/_kernel/liftover.py:356-371`) emits the count and the
`examples` payload on the same `INFO` line. Split it.

- [ ] **Step 1: Failing test.** Create
      `tests/test_liftover_logging.py` with one test verifying
      `log_liftover_drop_report` emits:
      - One `INFO` record with the count and reason but **no** `examples=`
        text, and including the `sidecar_path` when supplied.
      - One `DEBUG` record per call containing the `examples` payload.

      Use `caplog` at `DEBUG` level to capture both.

- [ ] **Step 2: Verify it fails.**

```bash
pytest tests/test_liftover_logging.py -v
```

- [ ] **Step 3: Refactor `log_liftover_drop_report`.** Add an optional
      `sidecar_path: str | os.PathLike | None = None` keyword argument.
      Change the body to:
      - `logger.info("%s dropped %d SNPs for %s.%s", workflow_label,
        report.n_dropped, report.reason, sidecar_suffix)` where
        `sidecar_suffix` is `f" See '{sidecar_path}' for row-level
        provenance."` when provided, else empty.
      - `logger.debug("%s example dropped (%s): %s", workflow_label,
        report.reason, report.examples)` only when `report.examples` is
        non-empty.

- [ ] **Step 4: Update existing call sites** to pass `sidecar_path` when
      they know the destination (Tasks 2 and 3 wire this up). For now,
      ensure no caller breaks: the parameter is keyword-only with a default.

- [ ] **Step 5: Re-run tests.** Both new and existing
      sumstats/ref-panel tests must stay green for any caller still passing
      no `sidecar_path`.

```bash
pytest tests/test_liftover_logging.py tests/test_sumstats_liftover.py tests/test_ref_panel_builder.py -v
```

---

## Task 2 — Sumstats `dropped_snps/dropped.tsv.gz` sidecar

**Files:** `src/ldsc/_kernel/liftover.py`, `src/ldsc/sumstats_munger.py`,
`tests/test_sumstats_liftover.py`.

Currently `apply_sumstats_liftover`
(`src/ldsc/_kernel/liftover.py:419-555`) returns `(lifted_df,
report_dict)` where `report_dict` is counts only. We need the per-reason
drop frames so the workflow layer can write a single sidecar at
`{output_dir}/dropped_snps/dropped.tsv.gz`. Per spec R2, the sidecar is
**always written** by every run — header-only TSV when no rows dropped.
There is no conditional-write logic for this artifact.

- [ ] **Step 1: Decide return shape.** Change
      `apply_sumstats_liftover` from `(lifted, report)` to
      `(lifted, report, drop_frame)`. Do **not** wrap in a dataclass —
      `drop_frame` is also stashed on
      `args._coordinate_metadata['liftover_drop_frame']` (Step 5 below)
      because it must travel through `_apply_liftover_if_requested` and
      back to the workflow layer; a dataclass gains nothing here and forces
      every existing call site to migrate to attribute access.

      Empty drop_frame for clean runs is a DataFrame with the unified
      schema columns and zero rows: `CHR: string`, `SNP: string`,
      `source_pos: Int64`, `target_pos: Int64`, `reason: string`. Do not
      return `None` — this lets callers concat unconditionally.

- [ ] **Step 2: Failing tests.** In `tests/test_sumstats_liftover.py`, add:
      - `test_sumstats_liftover_writes_dropped_snps_sidecar_with_all_reasons`:
        seed an input fixture that produces at least one row in each of
        `missing_coordinate`, `source_duplicate`, `unmapped_liftover`,
        `cross_chromosome_liftover`, `target_collision`. Assert the file
        exists at `{output_dir}/dropped_snps/dropped.tsv.gz`, has the
        expected columns, and contains at least one row per reason.
      - `test_sumstats_liftover_writes_header_only_sidecar_when_no_drops`:
        a clean run with no drops produces a `dropped_snps/dropped.tsv.gz`
        that exists, contains exactly the header row, and round-trips via
        `pd.read_csv(path, sep="\t", compression="gzip",
        dtype={"CHR": "string", "SNP": "string", "source_pos": "Int64",
        "target_pos": "Int64", "reason": "string"})` to a 0-row DataFrame
        whose dtypes match the unified sidecar schema. **Pandas cannot
        infer `string`/`Int64` from an empty TSV; the test must pass the
        dtype map explicitly. The writer is responsible only for emitting
        valid headers and proper NA encoding (`na_rep=""`); the consumer
        is responsible for declaring expected dtypes on read.**
      - `test_sumstats_dropped_snps_sidecar_preflight_blocks_existing_file`:
        pre-create the sidecar; run with `overwrite=False`; assert
        `FileExistsError`. (The sidecar is in `owned_paths`, so the
        no-overwrite collision still applies.)
      - `test_sumstats_dropped_snps_sidecar_overwritten_on_overwrite_true`:
        pre-create a sidecar with arbitrary stale content; run with
        `overwrite=True`; assert the file's content matches the current
        run (always overwritten in place).
      - `test_sumstats_log_points_to_dropped_snps_sidecar_at_info`: assert
        the `INFO` log line includes the resolved sidecar path.

- [ ] **Step 3: Verify they fail.**

```bash
pytest tests/test_sumstats_liftover.py -k "dropped_snps_sidecar or points_to_dropped" -v
```

- [ ] **Step 4: Modify `apply_sumstats_liftover`** to:
      - Collect each per-reason drop slice into a list of frames during the
        existing pipeline (missing-coord rows from
        `_normalized_coordinate_frame(..., drop_missing=True)`,
        `source_duplicate_result.report` rows, the unmapped/cross-chrom
        rows produced inside `_apply_chain_liftover` and
        `_apply_hm3_liftover`, and `target_duplicate_result.report` rows).
      - Concatenate them into a single `drop_frame` with columns
        `CHR: string, SNP: string, source_pos: Int64, target_pos: Int64,
        reason: string` (all nullable; `target_pos` is `<NA>` for
        `missing_coordinate`, `source_duplicate`, `unmapped_liftover`, and
        `cross_chromosome_liftover`; both `CHR` and `source_pos` may be
        `<NA>` for `missing_coordinate` rows).
      - Return the `(lifted, report, drop_frame)` 3-tuple. Update the
        signature and docstring.
      - The internal `_apply_chain_liftover` and `_apply_hm3_liftover`
        helpers gain a `collect_drops: bool` flag (or simply always return
        the per-reason mask + frame slice — pick whichever keeps the diff
        smaller).
      - **Migrate the production unpack site** at
        `src/ldsc/_kernel/sumstats_munger.py:965` to
        `dat, liftover_report, liftover_drop_frame =
        apply_sumstats_liftover(...)`. Stash the frame on
        `coordinate_metadata['liftover_drop_frame'] = liftover_drop_frame`
        immediately after `coordinate_metadata['liftover'] =
        liftover_report` so the existing
        `args._coordinate_metadata = coordinate_metadata` assignment
        carries it back to the workflow layer.
      - **Migrate the 6 test unpack sites** in
        `tests/test_sumstats_liftover.py` (around lines 90, 144, 185,
        219, 251, 299) to the 3-tuple. Add a stand-alone test asserting
        the returned `drop_frame` columns and dtypes match the spec.

- [ ] **Step 5: Modify `SumstatsMunger.run`** in
      `src/ldsc/sumstats_munger.py`:
      - Compute `dropped_snps_path = output_dir / "dropped_snps" /
        "dropped.tsv.gz"` once, up front.
      - **Preflight: always-owned, always-produced.** Add
        `dropped_snps_path` to **both** the `produced_paths` and
        `owned_paths` arguments of `preflight_output_artifact_family`
        unconditionally. The sidecar is class-1 (always written, even
        when zero rows dropped), so there is no asymmetry between owned
        and produced for this artifact. With `overwrite=False`, a stale
        sidecar raises `FileExistsError` (consistent with every other
        owned artifact). With `overwrite=True`, the writer overwrites it
        in place — no `remove_output_artifacts` call is needed for this
        path.
      - After the kernel call, read
        `drop_frame = coordinate_metadata.pop("liftover_drop_frame",
        None)` from the args (via the existing
        `getattr(args, "_coordinate_metadata", {})` retrieval at line
        420). If absent (e.g., no liftover requested), use an empty
        DataFrame with the unified schema dtypes — the file still gets
        written as a header-only TSV so the sidecar contract holds for
        every run.
      - **Always** ensure the parent directory exists and write
        `drop_frame.to_csv(dropped_snps_path, sep="\t", index=False,
        compression="gzip", na_rep="")`. The empty-frame case writes
        only the header line. Total disk overhead for an empty sidecar
        is ~30 bytes gzipped.
      - Emit one workflow-level `INFO` line summarizing the run: when
        `len(drop_frame) > 0` it names counts per reason and the
        sidecar path; when empty it states "no SNPs dropped during
        liftover stage; audit sidecar at <path>".
      - Add `dropped_snps_path` to `MungeRunSummary.output_paths` under
        key `dropped_snps_tsv_gz` **unconditionally** (the file is
        always written, so the key is always present in the mapping).
      - Per-reason example log lines from kernel calls inside
        `apply_sumstats_liftover` stay count-only (Task 1's no-pointer
        code path). The single workflow-level summary line above is
        the one that names the sidecar path. This keeps the kernel
        unaware of workflow-layer output paths.

- [ ] **Step 6: Add inline preflight comments at the modified call sites.**
      At `sumstats_munger.py:390` (around the
      `preflight_output_artifact_family` call updated in Step 5), add a
      short comment naming `owned_paths` vs `produced_paths` and noting
      that `dropped_snps/dropped.tsv.gz` appears in **both** lists
      unconditionally because it is always written, even on clean runs.
      Cross-reference the function-level docstring at
      `path_resolution.py:505` rather than restating the contract. At
      the writer call site, add a one-line comment explaining the
      always-write contract so the next reader does not "fix" the
      header-only-write as a no-op.

- [ ] **Step 7: Re-run sumstats tests.**

```bash
pytest tests/test_sumstats_liftover.py tests/test_sumstats_munger.py -v
```

---

## Task 3 — Reference-panel sidecar covers liftover-stage rows

**Files:** `src/ldsc/ref_panel_builder.py`, `tests/test_ref_panel_builder.py`.

Today, ref-panel `_resolve_mappable_snp_positions`
(`src/ldsc/ref_panel_builder.py:619-704`) logs unmapped/cross-chromosome
counts but never adds those rows to `dropped_frames`
(`src/ldsc/ref_panel_builder.py:471`). Extend it so the sidecar contains
every row removed by the liftover stage. Schema is unchanged but the
column dtypes widen to nullable per the unified schema (no behavior change
for existing rows since BIM `BP` is structurally non-null).

**Note on `missing_coordinate`.** Ref-panel sidecars will not contain
`missing_coordinate` rows. PLINK `.bim` `BP` is parsed as `int` at
`ref_panel_builder.py:455` (`chrom_metadata["POS"].astype(int)`) and again
at line 630 (`to_numpy(dtype=int)`); a malformed value raises before any
sidecar logic runs. The reason vocabulary is the union across workflows;
this asymmetry is intentional and documented in the spec/decisions doc.

- [ ] **Step 1: Failing tests.** In `tests/test_ref_panel_builder.py` add:
      - `test_dropped_snps_sidecar_contains_unmapped_and_cross_chromosome_rows`:
        fixture chain produces ≥1 unmapped and ≥1 cross-chromosome row;
        assert the per-chromosome sidecar's `reason` column contains both
        `unmapped_liftover` and `cross_chromosome_liftover` in addition
        to any duplicate rows.
      - `test_dropped_snps_sidecar_written_as_header_only_when_no_drops`:
        clean run with no drops on a chromosome produces
        `dropped_snps/chr{chrom}_dropped.tsv.gz` containing exactly the
        header line; round-trips via `pd.read_csv(path, sep="\t",
        compression="gzip", dtype={"CHR": "string", "SNP": "string",
        "source_pos": "Int64", "target_pos": "Int64", "reason": "string"})`
        to a 0-row DataFrame whose dtypes match. Per the unified sidecar
        contract, consumers must pass an explicit dtype map; pandas
        cannot infer nullable extension types from an empty TSV.
      - `test_dropped_snps_sidecar_overwritten_in_place_on_rerun`:
        Run 1 produces a sidecar with drops; Run 2 over the same dir
        with `--overwrite` and a different input that produces no drops
        on that chromosome; assert the sidecar now contains only the
        header line (overwritten in place; not deleted, not stale).
      - `test_dropped_snps_sidecar_in_preflight_blocks_no_overwrite`:
        pre-create `chr6_dropped.tsv.gz`; run targeting chr6 with
        `overwrite=False`; assert `FileExistsError`.
      - `test_dropped_snps_sidecar_written_header_only_on_restriction_empty_skip`:
        provide a `--ref-panel-snps-file` whose intersection with chr6
        is empty so `_build_chromosome` hits the restriction-empty
        early return (case 1). Assert `chr6_dropped.tsv.gz` exists
        and is header-only (the chromosome was iterated to; the
        always-write contract holds).
      - `test_dropped_snps_sidecar_written_with_target_collision_rows_on_target_dup_empty_skip`:
        seed an input where every retained SNP on chr6 maps to one of
        a small set of target-build positions, so target-collision
        filtering empties the chromosome (case 3). Assert
        `chr6_dropped.tsv.gz` exists, is **non-empty**, and contains
        `target_collision` rows for every dropped SNP (plus any
        upstream `source_duplicate`/`unmapped_liftover`/`cross_chromosome_liftover`
        rows that occurred). This locks in the symmetry with the
        source-duplicate-empty case.
      - `test_dropped_snps_sidecar_written_header_only_on_plink_filter_empty_skip`:
        configure a MAF filter that excludes every SNP on chr6 so
        `_build_chromosome` hits the PLINK-filter-empty early return
        (case 4) inside the per-build loop. Assert
        `chr6_dropped.tsv.gz` exists. Content is whatever the
        upstream always-write produced — header-only when no
        liftover-stage drops occurred, populated otherwise. Use a
        fixture where no liftover-stage drops happen so the assertion
        is "header-only" and deterministic.

- [ ] **Step 2: Verify they fail.**

```bash
pytest tests/test_ref_panel_builder.py -k "dropped_snps_sidecar" -v
```

- [ ] **Step 3: Modify `_resolve_mappable_snp_positions`** to
      additionally return a `liftover_drop_frame: pd.DataFrame` whose
      rows match the unified schema. **Every return branch must
      return a drop frame** — empty with the unified schema for
      branches that did not run liftover, populated with
      `unmapped_liftover` / `cross_chromosome_liftover` rows for
      branches that did.

      The function currently has four return points (`ref_panel_builder.py:619-704`):

      1. **hg19 source, no chain (line 638):** no liftover ran. Return
         an empty drop frame with the unified-schema dtypes
         (`CHR: string, SNP: string, source_pos: Int64, target_pos:
         Int64, reason: string`) and zero rows.
      2. **hg38 source, no chain (line 656):** same as (1) — empty
         unified-schema drop frame.
      3. **hg19 source, chain present (line 704, hg19 branch):** build
         the drop frame from the existing `unmapped_mask` and
         `cross_chrom_mask` via `mapping_reason_masks(mapping)` and the
         `drop_frame` constructed at lines 670-678. Use reason values
         `"unmapped_liftover"` and `"cross_chromosome_liftover"`.
         `target_pos` is `<NA>` for both. If no liftover-stage drops
         occurred (all positions mapped same-chromosome), return an
         empty unified-schema frame.
      4. **hg38 source, chain present (line 704, hg38 branch):** same
         as (3) — populated drop frame on liftover drops, empty
         unified-schema frame when nothing was dropped.

      Add a small private helper `_empty_unified_drop_frame() ->
      pd.DataFrame` that returns the canonical empty frame so the four
      return branches can call it instead of constructing it inline
      four times. The helper is also useful for the empty cases in
      `_build_chromosome` (Step 4) where no `dropped_frames` were
      collected.

      The new return signature becomes
      `tuple[np.ndarray, dict[int, int], dict[int, int], pd.DataFrame]`.
      Update the call site in `_build_chromosome` to unpack the
      fourth value and append it to `dropped_frames` **after the
      `_resolve_mappable_snp_positions` call returns and before the
      target-collision filtering block runs** (i.e., between the two
      duplicate blocks: source duplicates were already filtered before
      mapping; target collisions are about to be filtered using the
      newly-resolved positions). The liftover drops belong with the
      target-collision drops in the dropped-frames pipeline because
      both are post-mapping audit data.

- [ ] **Step 4: Modify `_build_chromosome`** to:
      - Append the liftover-drop frame returned by
        `_resolve_mappable_snp_positions` to `dropped_frames`
        **between the source-duplicate block (which runs before chain
        mapping) and the target-collision block (which runs after
        chain mapping)**. Concretely: insert the append immediately
        after the `_resolve_mappable_snp_positions(...)` call returns
        and before the second `_resolve_unique_snp_set(...)` call that
        does target-collision detection. The duplicate-result frames
        already use the same five-column schema, so concatenation is
        type-stable.
      - **Always** call `_write_dropped_sidecar` for the chromosome,
        regardless of whether `dropped_frames` is empty. When the
        concatenated `dropped_df` is empty, write a header-only TSV
        with the unified-schema dtypes. Eliminates the conditional
        write at lines 482-518 — the file is now always produced for
        every chromosome the current run processes. Total disk overhead
        for an empty per-chromosome sidecar is ~30 bytes gzipped.

      **Define "processed" — early-return cases.** `_build_chromosome`
      has **four** early-return points before the normal end-of-function
      sidecar write. The always-write contract holds for every
      chromosome the run *iterated to* (i.e., reached
      `_build_chromosome` for); each early-return case must therefore
      either write a sidecar before returning, or be downstream of a
      sidecar write that already occurred. The cases form a 2×2 grid:

      | Case | Drained by | Sidecar at return | Why |
      | --- | --- | --- | --- |
      | 1. Restriction-empty (`ref_panel_builder.py:467-469`) | `--ref-panel-snps-file` filtering, before any liftover stage | **Header-only** | Drops here are user restriction, not liftover stage; nothing belongs in `reason` rows. |
      | 2. Source-duplicate-empty (`ref_panel_builder.py:484-489`) | source `CHR/POS` collisions, before chain mapping | **Non-empty** with `source_duplicate` rows | The drops *are* liftover-stage drops; they belong in the audit. |
      | 3. **Target-duplicate-empty (`ref_panel_builder.py:519-521`)** | target `CHR/POS` collisions after chain mapping | **Non-empty** with `source_duplicate` (if any) + `unmapped_liftover` / `cross_chromosome_liftover` (if any) + `target_collision` rows | Symmetric to case 2: target collisions *are* liftover-stage drops; they belong in the audit. |
      | 4. PLINK-filter-empty (`ref_panel_builder.py:551-552`) | PLINK MAF / keep-individuals filtering inside the per-build loop | **Whatever was already written** by the unconditional sidecar write after coordinate duplicate/collision filtering | This case fires *after* that unconditional sidecar write, so the sidecar is already on disk; just return. The sidecar's contents reflect any liftover-stage drops that occurred upstream. |

      **Implementation rules per case:**
      - **Cases 1 (restriction-empty) and 4 (PLINK-filter-empty)** are
        the "drained by non-liftover filtering" cases. Case 1 must
        explicitly call `_write_dropped_sidecar(empty_frame, path,
        chrom)` before `return None`. Case 4 needs no sidecar action
        — by construction it fires after the unconditional sidecar
        write, so the file already exists.
      - **Cases 2 (source-duplicate-empty) and 3 (target-duplicate-empty)**
        are the "drained by liftover-stage filtering" cases. Both must
        call `_write_dropped_sidecar` before `return None` with the
        accumulated `dropped_frames` concatenated. Case 2 has only
        source-duplicate rows; case 3 has the full set of liftover-stage
        rows accumulated up to that point. After Step 4, both go through
        the same `_write_dropped_sidecar` code path — no separate
        branches.

      **Why all four cases write a sidecar (or rely on one already
      written):** the always-write contract says "every chromosome
      `_build_chromosome` reaches produces a sidecar at its expected
      path." If even one early-return path skipped the sidecar, the
      contract would be a lie at exactly the cases users most want to
      audit (chromosomes that emptied unexpectedly). The 2×2 grid
      ensures every drainage scenario leaves a usable artifact behind.

      Document the principle inline at the function header and at each
      early-return site: *"the per-chromosome sidecar is always written
      for any chromosome `_build_chromosome` is invoked on; reason
      rows reflect liftover-stage drops only — restriction-file and
      PLINK MAF drops do not appear because they happened outside the
      liftover stage."*

- [ ] **Step 5: Update `_write_dropped_sidecar`** to:
      - Handle empty input: when `dropped_df` is empty, write the
        header line only (no data rows). Emit a one-line `INFO`
        summary stating "no SNPs dropped on chromosome {chrom}; audit
        sidecar at <path>".
      - For non-empty input, replace the hardcoded counts with a
        `value_counts()` over the `reason` column and format the
        message with the categories present. Emit a one-line `INFO`
        summary naming the sidecar path.
      - **Do not** thread the sidecar path into
        `_resolve_mappable_snp_positions`. That helper runs before
        `_build_chromosome` is positioned to summarize. Per-call kernel
        logging from `_resolve_mappable_snp_positions` stays count-only
        via Task 1's no-pointer code path; the workflow emits exactly
        one summary line per chromosome from `_write_dropped_sidecar`.

- [ ] **Step 6: Update preflight wiring at `ref_panel_builder.py:160`.**
      The `_expected_ref_panel_output_paths` function (line 144)
      currently appends `dropped_snps/chr{chrom}_dropped.tsv.gz`
      unconditionally. After Step 4, that listing is now always
      *accurate* (the file is always written). Confirm no change to
      this function is needed — it already lists the path. Add a one-
      line comment confirming that the sidecar is always written and
      the preflight collision check therefore correctly applies under
      `overwrite=False`.

- [ ] **Step 7: Document the ref-panel class-2 asymmetry inline.** At
      `ref_panel_builder.py:233` (around the `ensure_output_paths_available`
      call), add a short comment that contrasts:
      - **Class 1 (always-written, this work's contract):** the new
        per-chromosome sidecar is always produced, so the preflight
        collision check is the entire contract — no auto-cleanup
        needed.
      - **Class 2 (per-chromosome r2/meta, per-build):** ref-panel
        uses the simpler `ensure_output_paths_available` helper, which
        does not auto-remove stale siblings under `--overwrite`. This
        is the documented expert contract from the module header
        (lines 21-26): users running against the same `output_dir`
        across runs with different chromosome scope or build
        configuration must manually clean stale per-chromosome /
        per-build files. See spec section R7 for the full rationale.

      Cross-reference the `path_resolution.py:505` docstring for
      callers wondering why the sumstats workflow uses a different
      helper for class-2 cleanup.

> **Note:** The originally-planned Step 8 — implementing
> `_warn_stale_class2_artifacts` to surface the silent-mix risk per
> spec R7 — was split out to its own spec/plan to keep this change set
> focused on liftover collision handling. See:
> - Spec: [`docs/superpowers/specs/2026-05-11-ref-panel-stale-class2-warning.md`](../specs/2026-05-11-ref-panel-stale-class2-warning.md)
> - Plan: [`docs/superpowers/plans/2026-05-11-ref-panel-stale-class2-warning.md`](./2026-05-11-ref-panel-stale-class2-warning.md)
>
> Spec R7's three governing principles stay in the parent spec because
> they shape per-workflow contracts; only the warning *implementation*
> moves.

- [ ] **Step 8: Re-run ref-panel tests.**

```bash
pytest tests/test_ref_panel_builder.py -v
```

---

## Task 4 — Remove `--duplicate-position-policy` and the `error` branch

**Files:** `src/ldsc/config.py`, `src/ldsc/ref_panel_builder.py`,
`tests/test_ref_panel_builder.py`.

Breaking change confirmed by user (2026-05-10). Once Tasks 2 and 3 land,
the sidecar plus log counts subsume every diagnostic the `error` policy
provided.

- [ ] **Step 1: Failing tests.** Update
      `tests/test_ref_panel_builder.py`:
      - Delete `ReferencePanelBuildConfigDuplicatePolicyTest` (`tests/test_ref_panel_builder.py:710-723`) and replace with one
        test asserting `ReferencePanelBuildConfig` has no
        `duplicate_position_policy` attribute and constructing one with
        `duplicate_position_policy=...` raises `TypeError`.
      - Delete `test_build_parser_defaults_duplicate_position_policy_to_drop_all`
        and the related parser test (`tests/test_ref_panel_builder.py:765-790`); add one asserting `ldsc build-ref-panel
        --help` output does not contain `--duplicate-position-policy` (use
        `argparse.format_help()` on the parser).
      - Delete the integration test that exercises
        `duplicate_position_policy="error"`
        (`tests/test_ref_panel_builder.py:1168` region) and the
        `drop-all`-explicit one (`...:1200`); a single test using the
        default already covers `drop-all`.
      - Delete `test_rsid_source_only_build_ignores_duplicate_position_policy_once_per_run`
        (`tests/test_ref_panel_builder.py:1239`); replace with
        `test_rsid_source_only_build_logs_duplicate_handling_skipped_once`
        that asserts the existing log message
        ("duplicate_position_policy applies only when
        snp_identifier='chr_pos'") is updated to a new wording (Step 4)
        and emitted exactly once.

- [ ] **Step 2: Verify they fail.**

```bash
pytest tests/test_ref_panel_builder.py -v
```

- [ ] **Step 3: Remove the config field.** In `src/ldsc/config.py`
      (`src/ldsc/config.py:444-516`):
      - Delete the `duplicate_position_policy` parameter docstring entry
        (lines around 444-451).
      - Delete `duplicate_position_policy: str = "drop-all"` (line 469).
      - Delete the `__post_init__` validation block (lines 512-516).

- [ ] **Step 4: Remove the CLI plumbing.** In
      `src/ldsc/ref_panel_builder.py`:
      - Delete the `--duplicate-position-policy` parser entry
        (`src/ldsc/ref_panel_builder.py:1198-1208`).
      - Delete `duplicate_position_policy=args.duplicate_position_policy`
        from `config_from_args` (line 1250).
      - In `_resolve_unique_snp_set` (lines 993-1112), drop the `policy`
        parameter and the two `if policy == "error":` branches (lines
        ~1027-1038 and ~1076-1091). Remove the
        `Use --duplicate-position-policy=drop-all to drop colliding
        clusters` mention.
      - In `_build_chromosome`, drop `policy=config.duplicate_position_policy`
        from the two `_resolve_unique_snp_set` calls
        (`src/ldsc/ref_panel_builder.py:480, 505`).
      - Update the rsid-mode info log
        (`src/ldsc/ref_panel_builder.py:246-249`) to reflect the new
        contract, e.g.:
        ```
        Coordinate duplicate filtering applies only when
        snp_identifier='chr_pos'; keeping duplicate CHR/POS rows in rsid
        source-only reference-panel builds.
        ```
        and align the test assertion in Step 1.
      - Update the **module docstring** at
        `src/ldsc/ref_panel_builder.py:21-26` to drop the
        `duplicate-position policy` mention from the "use a fresh
        output directory" trigger list. New wording:
        > `build-ref-panel` keeps an expert-oriented overwrite
        > contract: `overwrite=True` permits replacing candidate
        > artifacts for the current chromosome/build set but does
        > not remove stale optional target-build or `dropped_snps`
        > siblings from earlier configurations. Use a fresh output
        > directory when changing emitted builds, liftover/coordinate
        > configuration, or chromosome scope.

        This is the wording the spec quotes; keep them in sync.

- [ ] **Step 5: Verify ref-panel tests pass.**

```bash
pytest tests/test_ref_panel_builder.py -v
```

- [ ] **Step 6: CLI smoke check.**

```bash
ldsc build-ref-panel --help | grep -E "duplicate-position-policy" && echo FAIL || echo OK
```

      Must print `OK` (no match).

---

## Task 5 — Documentation

**Files:**
`docs/current/liftover-harmonization-decisions.md`,
`docs/current/data_flow.md`,
`docs/current/architecture.md`,
`tutorials/munge_sumstats.md` (or whichever tutorial covers munging — verify
with `ls tutorials/`),
`docs/superpowers/specs/2026-05-09-munge-sumstats-liftover.md`,
`docs/superpowers/plans/2026-04-30-duplicate-position-policy.md`.

- [ ] **Step 1: Update the harmonization decisions doc.** In
      `docs/current/liftover-harmonization-decisions.md`:
      - Reference-Panel Builder Contract section: replace the
        `duplicate_position_policy` line with "Coordinate duplicate
        groups (source and target) are dropped via `drop-all`. There is
        no public knob."
      - Replace the "duplicate-only sidecar" line with: "Per-chromosome
        sidecars at `dropped_snps/chr{chrom}_dropped.tsv.gz` are
        **always written** for every chromosome processed by the run
        (header-only when no drops). They capture rows removed in the
        liftover stage with `reason ∈ {source_duplicate,
        unmapped_liftover, cross_chromosome_liftover, target_collision}`.
        `missing_coordinate` is not produced by ref-panel because PLINK
        BIM `BP` is structurally non-null. Schema columns: `CHR`, `SNP`,
        `source_pos`, `target_pos`, `reason`; dtypes widened to nullable
        (`string`/`Int64`) to share one schema with sumstats."
      - Sumstats Contract section: add bullet for the new always-written
        `dropped_snps/dropped.tsv.gz` artifact (header-only when no
        drops), same schema, all five reasons may appear.
      - Backward Compatibility section: note (a) removal of
        `duplicate_position_policy` from CLI and Python API, and (b)
        ref-panel sidecar contract change from conditional-write to
        always-write — consumers that used `path.exists()` to detect
        "had drops" must switch to checking the row count.
      - Add a new **Class-2 (conditionally-produced) Artifacts** section
        that:
        - Names the three governing principles from spec R7: don't
          touch what we don't recognize; document the limitation as a
          user contract; inform loudly when divergence is possible.
        - Describes the per-workflow asymmetry kept by this work
          (sumstats auto-cleans `output_format` siblings; ref-panel
          does not auto-clean per-chromosome / per-build siblings;
          ref-panel users must "use a fresh output directory" when
          changing chromosome scope or build configuration).
        - Includes a short worked example of the silent-mix failure
          mode from spec section R7.
        - **Notes that the third principle ("inform loudly")** is
          implemented in a follow-up:
          [`docs/superpowers/specs/2026-05-11-ref-panel-stale-class2-warning.md`](../specs/2026-05-11-ref-panel-stale-class2-warning.md).
          That follow-up's documentation step will extend this Class-2
          Artifacts section with the warning's behavior, log signature,
          and remediation guidance — do not write that prose here.

- [ ] **Step 2: Update data flow / architecture diagrams.** Add the new
      sumstats artifact to the post-liftover output set in
      `docs/current/data_flow.md` and update any artifact list in
      `docs/current/architecture.md`.

- [ ] **Step 3: Tutorial note.** Add one short paragraph to the munging
      tutorial under a "Auditing dropped SNPs" subhead pointing to
      `dropped_snps/dropped.tsv.gz` and listing the five `reason` values.

- [ ] **Step 4: Append addenda to prior superpowers docs.**
      - In `docs/superpowers/specs/2026-05-09-munge-sumstats-liftover.md`,
        add `## 2026-05-10 Update` summarizing the policy collapse and
        sumstats sidecar (link to the new spec).
      - In `docs/superpowers/plans/2026-04-30-duplicate-position-policy.md`,
        prepend or append a note: "Superseded by
        `2026-05-10-liftover-harmonization.md`. The
        `duplicate_position_policy` flag/field has been removed."

---

## Final Verification

Run focused tests:

```bash
pytest tests/test_liftover_logging.py \
       tests/test_sumstats_liftover.py \
       tests/test_sumstats_munger.py \
       tests/test_ref_panel_builder.py -v
```

Run package-wide tests:

```bash
pytest -q
```

CLI smoke checks:

```bash
ldsc build-ref-panel --help | grep -c duplicate-position-policy   # expect 0
ldsc munge-sumstats --help | grep -E "target-genome-build|liftover-chain-file|use-hm3-quick-liftover"
```

End-to-end fixture validation:

1. Run `build-ref-panel` on a small fixture seeded with collisions across
   the four ref-panel-applicable reason categories (`source_duplicate`,
   `unmapped_liftover`, `cross_chromosome_liftover`, `target_collision`).
   Confirm `dropped_snps/chr*_dropped.tsv.gz` contains rows from each
   category; confirm `INFO`-level log shows counts and a sidecar pointer
   but no per-row examples. Confirm a chromosome with no drops produces a
   header-only sidecar (always-write contract).
2. Run `munge-sumstats` with chain liftover on a fixture seeded with all
   five sumstats-applicable categories: missing coordinates, source
   duplicates, unmapped rows, cross-chromosome hits, and target
   collisions. Confirm `dropped_snps/dropped.tsv.gz` contains rows from
   each category. Confirm a clean run produces a header-only sidecar.
3. Re-run (1) and (2) with `--log-level DEBUG`. Confirm example rows now
   appear inline.

Expected final state:

- `--duplicate-position-policy` is gone from CLI and Python API.
- Both workflows write a unified `dropped_snps/` sidecar always (header-
  only when no drops). Sumstats covers all five reasons; ref-panel covers
  four (no `missing_coordinate` — PLINK BIM `BP` is structurally
  non-null).
- Default `INFO` logs are count-only with sidecar pointers; `DEBUG` logs
  expose per-row examples.
- Sumstats `metadata.json` shape, ref-panel runtime metadata schemas, and
  all locked liftover validation rules are unchanged.
