# Design — Ref-Panel Stale Class-2 Artifact Warning

**Status.** Superseded by the current output-family cleanup contract. The
historical warning-only design below is retained for context, but
`build-ref-panel` now owns current-contract ref-panel artifacts, blocks them
without overwrite, removes stale current-contract siblings after successful
overwrites, and ignores removed legacy/root diagnostic names.

Spawned from spec section R7
of [`2026-05-10-liftover-harmonization.md`](./2026-05-10-liftover-harmonization.md)
to keep the liftover-harmonization change set focused on liftover
collision handling. The user-contract framing (R7's three governing
principles) was locked in the parent spec; this spec described an observable
warning-only behavior that is no longer the current implementation.

## Context

`build-ref-panel` produces two artifact classes:

- **Class 1 — always-produced.** Workflow always writes regardless of
  run flags. Examples: `build-ref-panel.log`, the new (per liftover
  harmonization) `dropped_snps/chr{chrom}_dropped.tsv.gz`.
- **Class 2 — conditionally-produced.** What gets written depends on
  per-run flags. Examples: per-chromosome
  `{output_dir}/{build}/chr{N}_r2.parquet` and
  `{output_dir}/{build}/chr{N}_meta.tsv.gz`, only for chromosomes
  processed by the current run; per-build subdirectories `hg19/` vs
  `hg38/`, only when the matching liftover chain is supplied.

Class-1 artifacts have no stale-cleanup problem because they are always
overwritten in place. Class-2 artifacts can leave a directory in an
internally-inconsistent state when a user reuses the same `output_dir`
across runs with different chromosome scope, build configuration, MAF
filter, or PLINK input. The current behavior — documented in the
`ref_panel_builder.py:21-26` module header — is "use a fresh output
directory" as a user contract: the workflow does not auto-clean stale
siblings, but it also does not warn about them. That silence is the
problem this spec addresses.

### The silent-mix failure mode

Worked example. Run 1 builds chr1-22 from `all_chroms.@` PLINK input.
Run 2 reuses the directory with `--plink-prefix only_chr6_7.@
--overwrite`. Run 2's preflight list contains only chr6 and chr7 paths
(loop in `_expected_ref_panel_output_paths` at
`ref_panel_builder.py:152` iterates over chromosomes resolved from the
current PLINK input, not all 22). The collision check passes for those
two — they exist from Run 1 but `--overwrite` allows replacement.
Chr1-5 and chr8-22 from Run 1 are never inspected because they are
not in Run 2's preflight list at all.

Result on disk: `{output_dir}/hg19/` holds chr6,7 from Run 2 and
chr1-5,8-22 from Run 1 side-by-side. A downstream tool that loads the
directory as a "complete reference panel" sees 22 chromosomes — built
from two different PLINK inputs, possibly with different sample sets,
MAF filters, or liftover configurations. The runtime metadata sidecars
(`chr*_meta.tsv.gz`) record the per-chromosome provenance, but most
consumers do not cross-check.

The user has no signal that this happened. The first user gets is a
downstream regression result that quietly diverges from prior runs.

## Recommended Behavior

Add a single `WARNING` log line per run that fires when ref-panel
detects pre-existing class-2 artifacts in `output_dir` that the current
run will not produce.

### Governing principles (carried from spec R7)

1. **Don't touch what we don't recognize.** The workflow remains
   deletion-free. The warning never deletes, renames, or modifies any
   file.
2. **Document the limitation as a user contract.** The user is
   responsible for keeping the output directory internally consistent.
   The workflow's job is to inform; the user's job is to act.
3. **Inform loudly when divergence is possible.** When pre-existing
   class-2 siblings exist outside the current run's expected-output set,
   surface them in a single `WARNING` log record.

### Behavior contract

- **When it fires.** After preflight passes (`ensure_output_paths_available`
  at `ref_panel_builder.py:233`), before chromosome processing begins.
  Firing early matters: a user who sees the warning can Ctrl-C and fix
  the directory before the run wastes hours of CPU.
- **Under both `--overwrite` modes.** The silent-mix risk exists whether
  or not the user opted into overwrite. `--overwrite` controls behavior
  for paths the run *will* touch; the warning concerns paths it *will
  not*.
- **Single record per run.** One `WARNING` log entry, not one per stale
  file. Listing files inside the record is fine; spamming the log is not.
- **No abort.** Chromosome-by-chromosome batch processing on HPC clusters
  is a legitimate workflow that produces the same on-disk pattern as
  the dangerous silent-mix case. The workflow has no way to distinguish
  "intentional batch with consistent config" from "accidental mix with
  divergent config" on its own. The warning is informational; the user
  decides whether to act.
- **Truncated list.** Cap at the first 10 paths sorted lexicographically.
  Append `"  ... and {N-10} more"` when N > 10. Without truncation, a
  22-chromosome × 2-build silent-mix scenario could dump 88 paths into
  a single log record.

### Detection scope

The helper inspects `output_dir/{build}/` for `build ∈ {"hg19", "hg38"}`
and globs `chr*_r2.parquet` and `chr*_meta.tsv.gz`. Paths whose `str`
representation is not in the current run's expected-output set are
collected as "stale class-2 siblings."

**Excluded from detection:**

- Class-1 artifacts: `dropped_snps/chr*_dropped.tsv.gz`,
  `build-ref-panel.log`, log siblings. These are always-written, so by
  construction they cannot be "stale" in the silent-mix sense.
- Files in subdirectories not named `hg19` or `hg38`. The workflow does
  not own arbitrary subdirectories of `output_dir`.
- Files with names that do not match the `chr*_r2.parquet` or
  `chr*_meta.tsv.gz` glob. The workflow does not own files outside its
  documented territory.

This conservative scope upholds principle 1: the workflow only inspects
what it would itself write under some configuration.

### Wording template

A single multi-line `LOGGER.warning(...)` record:

```
Found %d pre-existing reference-panel artifacts in '%s' that this run
will NOT touch. If those came from a previous run with different
configuration (different chromosome scope, source build, liftover
settings, MAF filter, etc.), downstream tools that load this directory
as a complete reference panel will silently see a mix of two different
builds. If this is an intentional chromosome-by-chromosome batch with
consistent configuration, no action is needed.

To remove the silent-mix risk, either:
  (a) start fresh with a new --output-dir, or
  (b) delete the listed files before re-running.

Pre-existing artifacts not produced by this run:
%s
```

The third `%s` is the truncated path list, one path per line indented
two spaces, optionally with the "and N more" tail.

## Why sumstats does not need an equivalent warning

Sumstats class-2 artifacts (`sumstats.parquet` vs `sumstats.sumstats.gz`)
are auto-cleaned by `remove_output_artifacts` after a successful write
(`sumstats_munger.py:474, 544`). When a user switches `--output-format`
between runs in the same directory, the prior format's file is removed
automatically. There is no silent-mix failure mode for sumstats users to
be warned about.

The two workflows' class-2 contracts remain asymmetric (per spec R7,
locked by the liftover-harmonization decision). This spec adds the
observability the asymmetry needs on the ref-panel side without changing
the underlying contract.

## Public Contract Changes

| Change | Magnitude | Approval |
|---|---|---|
| Add `WARNING` log line in `build-ref-panel` when stale class-2 siblings detected | additive — no abort, no deletion, no flag changes | minor: observability only |

No flag changes. No new config fields. No new files written. No
deletion semantics. Existing scripts continue to behave identically;
they just see one new `WARNING` line in the log under the silent-mix
condition.

## Backward Compatibility Impact

**None breaking.** The warning is purely additive. Users who reuse
`output_dir` legitimately (chromosome-by-chromosome batches with
consistent configuration) will see the warning fire on every run after
Run 1; they can either suppress at the log level (`--log-level ERROR`)
or ignore the message. Users who reuse `output_dir` with divergent
configuration get the silent-mix risk surfaced loudly — the intended
behavior.

The warning text explicitly names the legitimate-batch case so users
running per-chromosome workflows do not mistake the warning for an
error.

## Test / Doc Updates Needed

### Tests

- `test_warn_on_stale_class2_artifacts_logs_warning`: pre-create
  `out/hg19/chr1_r2.parquet` and `out/hg19/chr1_meta.tsv.gz`. Run
  `build-ref-panel` with `--plink-prefix only_chr6.@ --overwrite` so
  the run only produces chr6 paths. Assert a single `WARNING` log
  record fires that names both stale paths and includes the
  remediation guidance.
- `test_warn_on_stale_class2_artifacts_fires_under_no_overwrite_too`:
  same setup but `overwrite=False`; the existing collision check does
  not fire (chr6 paths don't exist), but the warning still fires for
  the chr1 stale paths.
- `test_no_warn_when_run_paths_match_existing`: pre-create chr6 paths;
  run targeting only chr6 with `--overwrite`; assert no stale-warning
  fires (every existing class-2 file is in the current run's expected
  set).
- `test_warn_truncates_stale_list_at_ten_paths`: pre-create 15 stale
  per-chromosome files; run targeting one different chromosome; assert
  the warning lists exactly the first 10 paths plus a "… and 5 more"
  tail.
- `test_warn_ignores_class1_artifacts`: pre-create
  `out/dropped_snps/chr1_dropped.tsv.gz` (class-1, always-written,
  never stale); run targeting chr6; assert the warning does not list
  the chr1 dropped-snps file.
- `test_warn_ignores_unknown_subdirectories`: pre-create
  `out/scratch/chr1_r2.parquet` (not under hg19/hg38); run targeting
  chr6; assert the warning does not list the scratch file (workflow
  does not own arbitrary subdirectories).

### Docs

- `docs/current/liftover-harmonization-decisions.md`: extend the
  Class-2 Artifacts section (added by the liftover harmonization
  work) to document the warning behavior and remediation guidance.
- `docs/superpowers/specs/2026-05-10-liftover-harmonization.md`:
  update R7's "Stale class-2 sibling warning" subsection to mark
  the implementation as landed (link this spec).

### Code touchpoints

- `src/ldsc/ref_panel_builder.py`: add private helper
  `_warn_stale_class2_artifacts(output_dir: Path,
  expected_paths: Iterable[Path]) -> None`; invoke it from `_run`
  immediately after `ensure_output_paths_available`.

## Verification Plan (post-implementation)

1. `pytest tests/test_ref_panel_builder.py -v` — all six warning tests
   pass; existing ref-panel tests stay green.
2. End-to-end: pre-create a fake `out/hg19/chr1_r2.parquet` and
   `out/hg19/chr1_meta.tsv.gz`. Run `ldsc build-ref-panel
   --plink-prefix small_chr6_only.@ --output-dir out --overwrite ...`.
   Confirm one `WARNING` log record names the chr1 files and includes
   the remediation guidance.
3. Run the same end-to-end without `--overwrite` (chr6 paths fresh).
   Confirm the warning still fires for chr1.
4. Run a clean fresh-directory build. Confirm no warning fires.
5. Run a build over a directory that contains 15 stale chr*_r2.parquet
   files. Confirm the warning lists 10 paths plus an "and 5 more" tail.

## Confirmed Decisions

User confirmed (2026-05-11):

1. **Split from liftover harmonization.** This warning behavior is
   independent of liftover collision handling and ships as its own
   spec/plan. The liftover work documents the principle (R7) without
   implementing the warning; this spec lands the implementation.
2. **Deletion-free.** The warning never modifies the filesystem.
3. **No abort.** Informational only.
4. **No new flags.** Always fires when the silent-mix condition holds.
