# Implementation Plan — Logging & Error Refinement (ldsc_py3_restructured)

**Date:** 2026-06-04
**Skill driving this work:** `my-skills:logging-refiner`
**Companion files:**
- Pilot template: `docs/superpowers/plans/2026-06-04-logging-refinement-pilot-munge.md`
- Session prompt: `docs/superpowers/plans/2026-06-04-logging-refinement-session-prompt.md`

---

## 1. Goal

Make every error that can **abort an `ldsc` run** point the user at its **cause**
and a **remedy** — **in the terminal message itself**. Back only the few
*multi-cause* errors with a browse-able, command-sectioned `docs/troubleshooting.md`
reference. Deliver the reference **and** the in-code rewrites in the same effort.

Decided scope (confirmed with the user):
- **Errors covered:** run-aborting only — anything that reaches a CLI command or
  public-API boundary and stops the operation. Pure internal invariant guards
  that should never fire are *summarized* under one internal-error pattern, not
  documented individually.
- **Messages are self-contained.** Each carries *what & where* + the most likely
  cause + the top remedy. **No error-ID codes.**
- **Reference doc:** `docs/troubleshooting.md`, organized by command, with
  plain-language `###` headings. An entry exists **only** for errors with ≥3
  distinct causes; the message links to it by the heading's slug. Most errors get
  no entry.
- **Code:** rewrite the in-code messages and exception classes together with the
  doc (not doc-first-then-code).

## 2. Non-goals

- Do **not** change control flow or behavior — only message text, log levels,
  exception *types*, exception *chaining* (`from e`), and the rare doc link.
- Do **not** rewrite the logging *infrastructure* in `_logging.py` (handlers,
  workflow-log-file routing) — it is already sound. Touch it only if a level is
  demonstrably wrong.
- Do **not** document warnings or DEBUG/INFO logs in the reference. (Log *quality*
  fixes — f-string style, level correctness, missing-context — are in scope as a
  secondary pass, but they get no doc entry.)
- Do **not** introduce error-ID codes or an immutable-ID scheme.
- Do **not** work in the `ldsc_py3_Jerry` (`main`) worktree.

## 3. Current state (from 2026-06-04 reconnaissance)

- **Exception architecture is already good.** `errors.py` defines
  `LDSCError → {LDSCUserError → LDSCUsageError, LDSCConfigError, LDSCInputError,
  LDSCDependencyError; LDSCInternalError}`. No redesign needed — the job is to
  *route raises into the right existing subclass*.
- **CLI boundary is already correct** (`cli.py:170 run_cli`): user errors → clean
  `Error: <msg>` + exit 1; `LDSCError` internal + any other `Exception` →
  `_log_internal_error` + exit 2, full traceback in the workflow log file.
- **The real gap is message quality**, at scale: **516 `raise` sites**, 166
  logging calls, 81 `except` blocks across 35 files. Most user-facing raises use
  **bare `ValueError`** with a terse "what" and no "where / why / fix".
- **Two structural smells to fix as we go:**
  1. `_USER_ERROR_TYPES` (`cli.py:34`) blanket-catches **all `ValueError`** as a
     user error. So a genuine *internal* `ValueError` bug is silently shown as a
     clean user message with exit 1. Migrating user raises to the
     `LDSCUserError` family lets us *narrow* this blanket later (removing it is
     out of scope now, but every migrated raise reduces reliance on it).
  2. `REGENERATE_ARTIFACT_MESSAGE` (`_kernel/snp_identity.py:36`) is one message
     reused for 4+ distinct failures — give each call site a message that names
     which check failed, link them all to one reference heading (pilot Part B).
- **Legacy kernel modules** (`_kernel/sumstats_munger.py`, `_kernel/ldscore.py`)
  still use `.format()`/`%`-style strings and leak developer phrasing to users
  (`'... This message indicates a bug.'`). Highest-severity rewrites live here.

## 4. Doc-section map (how the reference is organized)

`docs/troubleshooting.md` is sectioned by the **user-facing command**. Errors
raised in shared/kernel files are documented under the command the user runs when
they hit them; errors that appear across *several* commands go under a `Common`
section. This table maps source files to the doc section their multi-cause errors
land in — it is an *organization* aid, not an ID namespace.

| Doc section (`## ...`) | CLI surface | Primary source files |
|------------------------|-------------|----------------------|
| `munge-sumstats` | `munge-sumstats` | `sumstats_munger.py`, `_kernel/sumstats_munger.py` |
| `build-ref-panel` | `build-ref-panel` | `ref_panel_builder.py`, `_kernel/ref_panel_builder.py`, `_kernel/ref_panel.py`, `_kernel/plink_bed.py`, `_chr_sampler.py` |
| `ldscore` | `ldscore` | `ldscore_calculator.py`, `_kernel/ldscore.py` |
| `annotate` | `annotate` | `annotation_builder.py`, `_kernel/annotation.py` |
| `regression` | `h2`, `partitioned-h2`, `rg` | `regression_runner.py`, `_kernel/regression.py`, `_kernel/_irwls.py`, `_kernel/_jackknife.py` |
| `Common` | shared across commands | `config.py`, `path_resolution.py`, `outputs.py`, `genome_build_inference.py`, `chromosome_inference.py`, `_coordinates.py`, `_kernel/liftover.py`, `_kernel/snp_identity.py`, `_kernel/identifiers.py`, `_row_alignment.py`, `column_inference.py`, `hm3.py`, `hm3_reference.py` |

Put an error in `Common` only if it genuinely surfaces in 2+ commands; otherwise
file it under the single command where users meet it.

## 5. Per-module workflow (apply per command/batch)

For each batch, run the skill's phases on just that batch's raises:

1. **Enumerate run-aborting raises.** `grep -n "raise " <files>`. For each, decide:
   does it reach a CLI/API boundary (in scope) or is it an unreachable internal
   guard (leave as `LDSCInternalError`, summarize once)?
2. **Classify the track** for each: user input/config/usage (`LDSCUserError`
   family) vs. internal state (`LDSCInternalError`). Pick the most specific
   `errors.py` subclass.
3. **Count distinct causes per error.** Several raises that mean the same thing to
   a user are one error. Then:
   - **1–2 causes** → self-contained message only, **no** doc entry.
   - **≥3 distinct causes** → self-contained message for the top cause **plus** a
     `docs/troubleshooting.md` entry (ranked causes/checks/remedies) under the
     command section, linked by the heading slug.
4. **Author the reference entry** (multi-cause only) using the *exact* template in
   the pilot Part B: plain-language `###` heading, `Raised by` / `Exception`,
   `Symptom`, ranked *likely-cause → how-to-check* table, `Remedies`. Every check
   runnable or directly observable; every remedy concrete.
5. **Rewrite the in-code message** to the self-contained shape: *what & where* +
   *one* most-likely cause + *one* top remedy. Add the `#slug` link only for the
   multi-cause errors. Swap the exception class to the chosen subclass. Preserve
   every `from e` / `from None`.
6. **Secondary log pass (no doc entries):** within the same files, fix obvious log
   issues — `.format()`/`%` → f-string, wrong level, missing count/path context,
   `except ...: pass` swallows. Don't expand into new logging coverage unless a
   public I/O function has *zero* logs.
7. **Verify** (see §7) before moving on, and commit per batch.

## 6. Batch ordering

Lowest-level / highest-traffic first so shared messages settle before dependents:

1. **`munge-sumstats`** — already piloted; finish the remaining raises + apply the
   pilot Part C rewrites and fold the pilot Part B section into the doc.
2. **`Common`: config + path resolution** (`config.py` 49 raises, `path_resolution.py`)
   — surface in every command; settle them early.
3. **`Common`: identity / genome-build / liftover** — shared SNP/coordinate layer.
4. **`ldscore`** (`_kernel/ldscore.py` 56 raises — the largest single file).
5. **`build-ref-panel`**.
6. **`regression`**.
7. **`annotate`** + remaining `Common` (outputs, hm3).

Commit after each batch, e.g. `docs(troubleshooting): munge-sumstats section` and
`refactor(munge): self-contained error messages` (or one combined commit/batch).

## 7. Verification (run after every batch — CLAUDE.md "evidence before claims")

```bash
source /Users/wenbinwu/miniforge3/etc/profile.d/conda.sh && conda activate ldsc3-dev
cd /Users/wenbinwu/Documents_local/Research/SullivanLab/LDSC/repos/ldsc_py3_Jerry_workspace/ldsc_py3_restructured

pytest                                   # behavior must be unchanged

# No silent swallows remain:
grep -rn "except Exception: pass\|except Exception:\s*$\|except:\s*$" src --include="*.py" | grep -v __pycache__

# Every doc link in code resolves to a heading slug in the reference:
grep -rhoE "troubleshooting\.md#([a-z0-9-]+)" src --include="*.py" | sed 's/.*#//' | sort -u
grep -oE "^### .+" docs/troubleshooting.md   # eyeball: each code link above appears here (as the heading's slug)
```

Tests must pass with identical *behavior*. If a test asserts on exact error
*text*, update the assertion to the new text in the same commit and note it.

## 8. Reference upkeep rule (install once, in batch 1)

Append to `CLAUDE.md` (the repo convention) under a new `## Troubleshooting
reference` section:

```markdown
## Troubleshooting reference
`docs/troubleshooting.md` is a browse-able, command-sectioned reference for
run-aborting errors that have multiple causes. When you add or change such an error:
- Make the in-code message self-contained: what & where + the most likely cause
  + the top remedy.
- If the error has 3+ distinct causes, add/update its section (plain-language
  heading, ranked causes + checks + remedies) and link the message to that
  heading's slug.
- Keep any in-code `docs/troubleshooting.md#...` link in sync with its heading.
```

## 9. Risks & gotchas

- **Line numbers drift** as edits land. Re-grep before each edit; never trust a
  stale line number from this plan or the pilot.
- **Blanket `ValueError` catch** means a wrong exception class won't *crash* but
  will mis-route. Choosing the right `errors.py` subclass is the point — don't
  leave new raises as bare `ValueError`.
- **Slug drift:** the `#slug` in a message must match the GitHub auto-slug of the
  `###` heading exactly. If you reword a heading a message links to, update the
  link in the same commit (verification in §7 catches mismatches).
- **Over-linking:** resist adding a doc link to simple 1–2-cause errors. The
  message should stand alone; links are the exception, not the rule.
- **Legacy `.format()` strings** in kernels: convert to f-strings but re-check the
  interpolated variable names still exist after refactor.
- **Tests asserting on message text** will break on rewrites — expect to update a
  handful; this is allowed (text, not behavior).
- **`from None`** appears intentionally in a few places — preserve intent.
- Keep messages tight (≤ ~2 sentences inline). Extra causes belong in the doc.

## 10. Definition of done

- [ ] Every run-aborting in-code message is self-contained (what & where + top
      cause + top remedy) and raises the correct `errors.py` subclass.
- [ ] `docs/troubleshooting.md` exists, sectioned by command, with a plain-language
      entry for each ≥3-cause error — and **no** entries for simple errors.
- [ ] Every in-code `docs/troubleshooting.md#...` link resolves to a real heading
      slug (§7 check passes); no orphan links, no over-linking of simple errors.
- [ ] No `except ...: pass` silent swallows remain in `src/`.
- [ ] Secondary log pass done in touched files (f-string style, levels, context).
- [ ] `CLAUDE.md` has the troubleshooting-reference upkeep section.
- [ ] `pytest` passes; any text-assertion updates noted in commits.
