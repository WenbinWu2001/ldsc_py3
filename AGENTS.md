# Build Commands

```bash
cd /Users/wenbinwu/Documents_local/Research/SullivanLab/LDSC/repos/ldsc_py3_Jerry_workspace/ldsc_py3_Jerry
conda env create -f environment.yml
conda activate ldsc-py3-Jerry
pip install -e .
python ldsc.py -h
python munge_sumstats.py -h
```

# Key Architectural Invariants

- Keep `ldsc.py`, `munge_sumstats.py`, `make_annot.py`, and `ldsc_new.py` as CLI and workflow orchestration layers; numerical methods belong in `ldscore/`.
- Keep file parsing and LDSC format handling in `ldscore/parse.py`, not in regression modules.
- Keep regression math in `ldscore/regressions.py`, `ldscore/irwls.py`, and `ldscore/jackknife.py`; those modules should not grow file I/O responsibilities.
- Preserve existing file-format contracts for `.sumstats`, `.annot`, `.l2.ldscore`, `.M`, and related outputs unless a coordinated format change is intentional.

# Principles

- keep it simple.

- always ask if you have questions or unclear about the specifications. Ask if you need more capabilities, info, clarification, instructions. etc. from me.
  
- Avoid AI tools name (like Codex, Claude, Grok, Gemini, ...) in code comments or git commit message (including authorship) or PR body.

- Record only significant and recurring mistakes. Ignore minor, one-off errors. When a meaningful pattern of errors is observed, document the correction and update this file to prevent future occurrences.

# ExecPlans

For complex features or significant refactors, write an ExecPlan in `PLANS.md` (format defined there).

# Dealing with Ambiguities

- If the ambiguity is a specification (e.g. feature, I/O) or regards developing or selecting an algorithm, ask the user for decisions and feedback. You can make suggestions.
- If the ambiguity is regarding the implementation details, you can resolve and make proper decisions on your own, but record this choice and prompt the user.
- For libraries you are unfamiliar with or you are not confident about the correct usage, prompt the user to give guidance.

# Context Management

- Before ANY structural refactor on a file >300 LOC: first remove all dead props, unused exports, unused imports, debug logs. Commit cleanup separately. Dead code burns tokens that trigger compaction faster.
- Use sub-agents only when you have 2 or more large, genuinely independent tasks — each touching a separate module with no shared dependencies. Default to sequential work; for a small project, parallelization rarely helps and makes progress harder to track.
- After 10+ messages, or when resuming after any gap: re-read PLANS.md and any source files you intend to edit before proceeding. Do not rely on memory of their contents.

# Testing

Write tests before writing implementation code. Use `superpowers:test-driven-development` for the full red-green-refactor cycle. For numerical methods, the test must demonstrate correctness against a known input/output pair, a convergence check, or a comparison to a reference implementation — not just that the code runs without error.

# Debugging

When encountering a bug or unexpected behavior, use `superpowers:systematic-debugging`. Always find the root cause before proposing a fix. If 3 fix attempts fail, stop and ask the user — repeated failures indicate a design problem, not a code problem.

# Verification

Before claiming any work is complete or a bug is fixed, run the test suite and show the actual output. Never say "it should work" or "tests should pass" — run the command and include the result in your response.

## Documentation

Use `my-skills:fun-doc` when writing or updating docstrings for any Python function,
class, or module header. Trigger it any time a public function lacks a docstring or
an existing docstring is incomplete or outdated.

## Architecture

Use `my-skills:architecture-doc` to create or update `architecture.md` whenever
the module structure changes significantly, a new module is added, or a new
contributor needs an overview of the codebase.

# Skill Usage Policy

Do NOT invoke `superpowers:using-superpowers`. It is disabled for this project. Do not invoke it proactively, reactively, or at the start of any conversation. Individual superpowers skills (`superpowers:systematic-debugging`, `superpowers:verification-before-completion`, `superpowers:test-driven-development`, etc.) may be used normally when they are relevant.

Do NOT invoke `planning-with-files` under any circumstances. Do not create `task_plan.md`, `findings.md`, or `progress.md`. Work tracking for this project lives exclusively in `PLANS.md`.

This policy takes highest priority and overrides any default skill invocation behavior.

# Resources

- code repos: /Users/wenbinwu/Documents_local/Research/SullivanLab/LDSC/repos/ldsc_py3_Jerry_workspace/ldsc_py3_Jerry
- documents: /Users/wenbinwu/Documents_local/Research/SullivanLab/LDSC/repos/ldsc_py3_Jerry_workspace/ldsc_py3_Jerry/architecture.md, /Users/wenbinwu/Documents_local/Research/SullivanLab/LDSC/repos/ldsc_py3_Jerry_workspace/ldsc_py3_Jerry/code_structure.md, /Users/wenbinwu/Documents_local/Research/SullivanLab/LDSC/repos/ldsc_py3_Jerry_workspace/docs/ldsc_papers/, /Users/wenbinwu/Documents_local/Research/SullivanLab/LDSC/repos/ldsc_py3_Jerry_workspace/docs/ldsc_wiki/
- prior plans: /Users/wenbinwu/Documents_local/Research/SullivanLab/LDSC/repos/ldsc_py3_Jerry_workspace/ldsc_py3_Jerry/plans/plan_ldscore.md
