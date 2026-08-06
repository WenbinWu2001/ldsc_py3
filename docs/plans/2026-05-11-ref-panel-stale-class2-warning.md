# Ref-Panel Stale Class-2 Artifact Warning Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use
> `superpowers:subagent-driven-development` (recommended) or
> `superpowers:executing-plans` to implement this plan task-by-task.
> Steps use checkbox (`- [ ]`) syntax for tracking. Spec:
> [`docs/superpowers/specs/2026-05-11-ref-panel-stale-class2-warning.md`](../specs/2026-05-11-ref-panel-stale-class2-warning.md).

**Goal.** Implement the `WARNING` log line that surfaces the silent-mix
risk for stale class-2 artifacts in `build-ref-panel` output
directories. Spawned from spec R7 of the liftover-harmonization work.

**Status.** Implemented on 2026-05-11. The helper is
`src/ldsc/ref_panel_builder.py::_warn_stale_class2_artifacts`, and focused
coverage lives in `tests/test_ref_panel_builder.py::RefPanelStaleClass2WarningTest`.

**Architecture.** A single private helper
`_warn_stale_class2_artifacts(output_dir, expected_paths)` in
`src/ldsc/ref_panel_builder.py`, invoked from `_run` immediately after
`ensure_output_paths_available` (line 233). The helper scans
`{output_dir}/{hg19,hg38}/` for `chr*_r2.parquet` and
`chr*_meta.tsv.gz` files not in the current run's expected set, emits
one `WARNING` record naming up to 10 affected paths, and returns. The
helper is deletion-free, never aborts, and fires under both
`--overwrite` modes. No flag changes, no config changes, no contract
changes for any other artifact.

**Tech Stack.** Python 3.11, pathlib, logging, pytest with `caplog`.

---

## File Map

| File | Change |
|---|---|
| `src/ldsc/ref_panel_builder.py` | Add private helper `_warn_stale_class2_artifacts(output_dir, expected_paths) -> None`. Invoke it from `_run` after `ensure_output_paths_available` at line 233. Helper globs `{output_dir}/{hg19,hg38}/chr*_r2.parquet` and `chr*_meta.tsv.gz`, computes `existing - expected`, emits one `LOGGER.warning(...)` record with the spec-defined wording template and 10-path truncation. |
| `tests/test_ref_panel_builder.py` | Add six tests verifying: warning fires under `--overwrite=True`; warning fires under `--overwrite=False`; warning suppressed when no stale; truncation at 10 with "and N more" tail; class-1 artifacts (`dropped_snps/chr*_dropped.tsv.gz`) excluded from detection; files in non-hg19/hg38 subdirectories excluded. |
| `docs/current/liftover-harmonization-decisions.md` | Extend the Class-2 Artifacts section (added by the liftover harmonization work) with the warning behavior and remediation guidance. |
| `docs/superpowers/specs/2026-05-10-liftover-harmonization.md` | Update R7's "Stale class-2 sibling warning" subsection to mark implementation as landed; add link to this plan. |

---

## Pre-Flight

- [ ] Confirm the liftover-harmonization work (the parent change set
      that added the Class-2 Artifacts section to
      `docs/current/liftover-harmonization-decisions.md`) has landed.
      This spec depends on R7's principles being documented; the
      warning extends those.
- [ ] Confirm worktree is clean.
- [ ] Run the existing suite as a baseline:

```bash
source /Users/wenbinwu/miniforge3/etc/profile.d/conda.sh && conda activate ldsc3
pytest -q tests/test_ref_panel_builder.py
```

Capture the green baseline so regressions during the change are obvious.

---

## Task 1 — Add the warning helper and call site

**Files:** `src/ldsc/ref_panel_builder.py`,
`tests/test_ref_panel_builder.py`.

- [ ] **Step 1: Failing tests.** Add to `tests/test_ref_panel_builder.py`:

      - `test_warn_on_stale_class2_artifacts_logs_warning`: pre-create
        `out/hg19/chr1_r2.parquet` and `out/hg19/chr1_meta.tsv.gz`. Run
        `build-ref-panel` with `--plink-prefix only_chr6.@ --overwrite`
        so the run only produces chr6 paths. Use `caplog.at_level("WARNING")`
        to capture log records. Assert exactly one `WARNING` record
        fires whose message contains both stale paths, the
        "fresh --output-dir" remediation phrase, and the
        "delete the listed files" remediation phrase.

      - `test_warn_on_stale_class2_artifacts_fires_under_no_overwrite_too`:
        same setup but `overwrite=False`. The chr6 paths do not exist
        beforehand, so the existing collision check at line 233 does
        not fire. Assert the new warning still fires for the chr1
        stale paths (warning is independent of `--overwrite`).

      - `test_no_warn_when_run_paths_match_existing`: pre-create chr6
        paths only. Run targeting only chr6 with `--overwrite`. Assert
        no `WARNING` records mention "pre-existing reference-panel
        artifacts" (the warning's signature phrase). Every existing
        class-2 file is in the current run's expected set, so nothing
        is stale.

      - `test_warn_truncates_stale_list_at_ten_paths`: pre-create 15
        stale per-chromosome files (chr1-15) under `out/hg19/`. Run
        targeting chr16. Assert the warning's message lists exactly 10
        paths and contains the literal substring `"... and 5 more"`.

      - `test_warn_ignores_class1_artifacts`: pre-create
        `out/dropped_snps/chr1_dropped.tsv.gz`. Run targeting chr6.
        Assert the warning does not fire (no class-2 stale, only
        class-1 which is always-written).

      - `test_warn_ignores_unknown_subdirectories`: pre-create
        `out/scratch/chr1_r2.parquet` (subdirectory not named hg19 or
        hg38). Run targeting chr6 in hg19. Assert the warning does not
        list `out/scratch/chr1_r2.parquet` (workflow does not own
        arbitrary subdirectories).

- [ ] **Step 2: Verify they fail.**

```bash
source /Users/wenbinwu/miniforge3/etc/profile.d/conda.sh && conda activate ldsc3
pytest tests/test_ref_panel_builder.py -k "warn_on_stale or no_warn_when or warn_truncates or warn_ignores" -v
```

Expected: all six fail with "no WARNING record" assertions because the
helper does not exist yet.

- [ ] **Step 3: Implement the helper.** In
      `src/ldsc/ref_panel_builder.py`, add a private helper above
      `class ReferencePanelBuilder` (e.g., after the
      `_expected_ref_panel_output_paths` definition):

      ```python
      def _warn_stale_class2_artifacts(
          output_dir: Path,
          expected_paths: Iterable[Path],
      ) -> None:
          """Surface the silent-mix risk for class-2 ref-panel artifacts.

          Implements spec R7's "inform loudly when divergence is
          possible" principle without modifying the filesystem.
          See docs/superpowers/specs/2026-05-11-ref-panel-stale-class2-warning.md.
          """
          # Class-2 territory: per-chromosome r2 + meta files under
          # hg19/ or hg38/ subdirs. Class-1 sidecars (dropped_snps,
          # logs) are always-written and cannot be "stale," so they
          # are excluded from this check.
          expected = {str(Path(p)) for p in expected_paths}
          stale: list[Path] = []
          for build_subdir_name in ("hg19", "hg38"):
              build_subdir = output_dir / build_subdir_name
              if not build_subdir.is_dir():
                  continue
              for pattern in ("chr*_r2.parquet", "chr*_meta.tsv.gz"):
                  for path in build_subdir.glob(pattern):
                      if str(path) not in expected:
                          stale.append(path)
          if not stale:
              return
          stale_sorted = sorted(stale, key=lambda p: str(p))
          shown = stale_sorted[:10]
          tail = (
              f"  ... and {len(stale_sorted) - 10} more"
              if len(stale_sorted) > 10
              else None
          )
          path_lines = "\n".join(f"  {p}" for p in shown)
          if tail is not None:
              path_lines = f"{path_lines}\n{tail}"
          LOGGER.warning(
              "Found %d pre-existing reference-panel artifacts in '%s' "
              "that this run will NOT touch. If those came from a previous "
              "run with different configuration (different chromosome scope, "
              "source build, liftover settings, MAF filter, etc.), downstream "
              "tools that load this directory as a complete reference panel "
              "will silently see a mix of two different builds. If this is an "
              "intentional chromosome-by-chromosome batch with consistent "
              "configuration, no action is needed.\n\n"
              "To remove the silent-mix risk, either:\n"
              "  (a) start fresh with a new --output-dir, or\n"
              "  (b) delete the listed files before re-running.\n\n"
              "Pre-existing artifacts not produced by this run:\n%s",
              len(stale_sorted),
              output_dir,
              path_lines,
          )
      ```

      Notes:
      - The `expected` set is built from `str(Path(p))` to match the
        glob output exactly.
      - Sort lexicographically before truncation so the displayed
        first-10 set is deterministic across runs.
      - Use a single `LOGGER.warning(...)` call so the entire message
        appears as one log record (matters for `caplog` assertions and
        for log readers that group by record).

- [ ] **Step 4: Wire the call site.** In
      `ReferencePanelBuilder._run` at `ref_panel_builder.py:233`,
      immediately after the `ensure_output_paths_available(...)` call
      returns, add:

      ```python
      # Spec R7: don't touch what we don't recognize; inform the user
      # loudly about pre-existing class-2 siblings outside the current
      # run's expected set. See
      # docs/superpowers/specs/2026-05-11-ref-panel-stale-class2-warning.md.
      _warn_stale_class2_artifacts(Path(config.output_dir), preflight_paths)
      ```

      Note that `preflight_paths` already contains the class-1 sidecar
      paths — that is fine because the helper inspects only the
      `chr*_r2.parquet` and `chr*_meta.tsv.gz` glob patterns under
      `hg19/` and `hg38/`, so class-1 entries in `expected` are
      harmless extras that simply never match the scan.

- [ ] **Step 5: Re-run tests.**

```bash
pytest tests/test_ref_panel_builder.py -k "warn_on_stale or no_warn_when or warn_truncates or warn_ignores" -v
```

Expected: all six pass.

- [ ] **Step 6: Re-run the full ref-panel suite to catch regressions.**

```bash
pytest tests/test_ref_panel_builder.py -v
```

---

## Task 2 — Documentation

**Files:**
`docs/current/liftover-harmonization-decisions.md`,
`docs/superpowers/specs/2026-05-10-liftover-harmonization.md`.

- [ ] **Step 1: Extend the Class-2 Artifacts section.** In
      `docs/current/liftover-harmonization-decisions.md`, find the
      Class-2 Artifacts section added by the liftover harmonization
      work. Add a new subsection titled "Stale-sibling warning" that:
      - States when the warning fires (after preflight, before
        chromosome processing, under both `--overwrite` modes).
      - Explains why it does not abort (legitimate
        chromosome-by-chromosome batches produce the same on-disk
        pattern as silent-mix).
      - Lists the two remediation options (fresh `--output-dir` or
        manual cleanup).
      - Quotes the warning's signature phrase
        ("Found N pre-existing reference-panel artifacts...") so
        users grepping logs can find it.
      - Cross-references this spec by file path.

- [ ] **Step 2: Update parent spec R7.** In
      `docs/superpowers/specs/2026-05-10-liftover-harmonization.md`,
      find the "Stale class-2 sibling warning (added by this work)"
      subsection inside R7. Replace its forward-looking language
      ("ref-panel emits...") with completed-language plus a link to
      this spec, e.g.:

      > **Stale class-2 sibling warning.** Implementation landed in
      > [`docs/superpowers/specs/2026-05-11-ref-panel-stale-class2-warning.md`](./2026-05-11-ref-panel-stale-class2-warning.md).
      > Behavior summary: ref-panel emits one `WARNING` log record per
      > run when pre-existing class-2 artifacts in `output_dir` are not
      > in the current run's expected set; deletion-free; does not
      > abort; fires under both `--overwrite` modes.

      Adjust adjacent paragraphs as needed for tense agreement.

---

## Final Verification

Run focused tests:

```bash
pytest tests/test_ref_panel_builder.py -v
```

Run package-wide tests:

```bash
pytest -q
```

End-to-end fixture validation:

1. Pre-create fake `out/hg19/chr1_r2.parquet` and
   `out/hg19/chr1_meta.tsv.gz` (just touch them or write small dummy
   parquet/tsv files). Run `ldsc build-ref-panel --plink-prefix
   small_chr6_only.@ --output-dir out --overwrite ...`. Confirm one
   `WARNING` log record names the chr1 files and includes both
   remediation options.
2. Same as (1) without `--overwrite` (chr6 paths fresh). Confirm the
   warning still fires for chr1.
3. Run a clean build into a fresh `--output-dir`. Confirm no warning
   fires.
4. Pre-create 15 stale files under `out/hg19/`; run a build into the
   same `out/`. Confirm the warning lists exactly 10 paths plus an
   "and 5 more" tail.

Expected final state:

- `_warn_stale_class2_artifacts` exists in
  `src/ldsc/ref_panel_builder.py` and is called from `_run` after the
  preflight check.
- Six tests in `tests/test_ref_panel_builder.py` cover the behavior
  matrix (under both `--overwrite` modes, no-warn case, truncation,
  class-1 exclusion, unknown-subdirectory exclusion).
- `docs/current/liftover-harmonization-decisions.md` documents the
  warning under the Class-2 Artifacts section.
- Parent spec R7 in
  `docs/superpowers/specs/2026-05-10-liftover-harmonization.md`
  references this spec as the implementation landing.
- No new flags, no new config fields, no deletion semantics, no
  changes to any other artifact's contract.
