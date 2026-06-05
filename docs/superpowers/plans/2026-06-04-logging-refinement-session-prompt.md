# Session Prompt — paste this to start the implementation session

> Copy everything in the fenced block below into a fresh Claude Code session
> started from the workspace root
> (`/Users/wenbinwu/Documents_local/Research/SullivanLab/LDSC/repos/ldsc_py3_Jerry_workspace`).
> It is self-contained.

```text
Work in the `ldsc_py3_restructured` repo on the `restructure` branch ONLY
(path: .../ldsc_py3_Jerry_workspace/ldsc_py3_restructured). Do NOT touch the
ldsc_py3_Jerry (main) worktree.

Task: refine logging and error messages so every run-aborting error tells the
user — IN THE TERMINAL MESSAGE ITSELF — what & where failed, the most likely
cause, and the top fix. Back only the few errors that have many causes with a
browse-able docs/troubleshooting.md reference. Deliver the reference AND the
in-code rewrites together. NO error-ID codes.

First, invoke the `my-skills:logging-refiner` skill and follow it. Then read,
in this order, and treat them as the source of truth:
  1. docs/superpowers/plans/2026-06-04-logging-refinement-plan.md   (the plan)
  2. docs/superpowers/plans/2026-06-04-logging-refinement-pilot-munge.md (the
     gold-standard template: copy its self-contained in-code message shape
     (Part C) and its command-sectioned reference format (Part B) EXACTLY)
Also read CLAUDE.md (build/test commands, "ask before non-trivial work" rule,
commit conventions) and errors.py (the existing exception hierarchy to route
into — do NOT invent new exception classes).

Scope (already decided — do not re-litigate):
  - Handle ONLY run-aborting errors (reach a CLI/API boundary and stop the run).
    Internal invariant guards that should never fire stay LDSCInternalError and
    are summarized once, not documented individually.
  - Every message is SELF-CONTAINED: what & where + most likely cause + top
    remedy, inline. No error-ID codes.
  - docs/troubleshooting.md is organized by COMMAND with plain-language ###
    headings (section map in plan §4). Add an entry ONLY for an error with 3+
    distinct causes; the message links to it by the heading's slug. Simple 1-2
    cause errors get a self-contained message and NO doc entry — do not over-link.
  - Swap each user-facing raise from bare ValueError to the correct errors.py
    subclass. Behavior must NOT change — only message text, exception type,
    exception chaining, log level, and the rare doc link.

Method: work ONE batch at a time in the order in plan §6 (start with
munge-sumstats, whose pilot is already written — fold the pilot Part B section
into docs/troubleshooting.md, finish the remaining munge raises, then apply the
pilot Part C rewrites). For each batch follow plan §5 (enumerate raises, classify
track, count distinct causes per error, author reference entries for 3+-cause
errors only, rewrite messages self-contained, secondary log-quality pass, verify,
commit).

Before editing each batch: re-grep the raise sites (line numbers in the
plan/pilot will have drifted — never trust stale numbers). Per CLAUDE.md, for
anything beyond a several-line change, ask clarifying questions first and pause
for plan-mode approval before large edits.

Verify after every batch (per CLAUDE.md "evidence before claims"):
  source /Users/wenbinwu/miniforge3/etc/profile.d/conda.sh && conda activate ldsc3-dev
  cd .../ldsc_py3_restructured && pytest
  grep -rn "except Exception: pass\|except Exception:\s*$\|except:\s*$" src --include="*.py" | grep -v __pycache__
  # every doc link in code must resolve to a heading slug in the reference:
  grep -rhoE "troubleshooting\.md#([a-z0-9-]+)" src --include="*.py" | sed 's/.*#//' | sort -u
  grep -oE "^### .+" docs/troubleshooting.md
Tests must pass; if a test asserts exact error TEXT, update that assertion in the
same commit and note it. Commit per batch using Conventional Commits (no AI-tool
names in messages/authorship).

Once (in batch 1): install the troubleshooting-reference upkeep rule into
CLAUDE.md (exact text in plan §8).

Stop and check in with me after the first batch (munge-sumstats) is fully done
and verified, so I can confirm the message tone and reference format before you
scale to the rest.
```
