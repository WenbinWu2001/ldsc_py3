# Build Commands

```bash
cd /Users/wenbinwu/Documents_local/Research/SullivanLab/LDSC/repos/ldsc_py3_Jerry_workspace/ldsc_py3_restructured
conda env create -f environment.yml      # first-time environment setup
source /Users/wenbinwu/miniforge3/etc/profile.d/conda.sh
conda activate ldsc3
pip install -e .                         # install editable package
ldsc --help                              # installed console entry point
python -m ldsc --help                    # module entry point
pytest                                   # primary test suite
python -m unittest discover -s tests -p 'test*.py' -v
```

For non-interactive agent commands, activate the project environment with:

```bash
source /Users/wenbinwu/miniforge3/etc/profile.d/conda.sh && conda activate ldsc3 && <command>
```

# Key Architectural Invariants

- Always make code, test, and documentation changes in the `restructure`
  branch worktree:
  `/Users/wenbinwu/Documents_local/Research/SullivanLab/LDSC/repos/ldsc_py3_Jerry_workspace/ldsc_py3_restructured`.
- Do not work on the `main` branch worktree
  `/Users/wenbinwu/Documents_local/Research/SullivanLab/LDSC/repos/ldsc_py3_Jerry_workspace/ldsc_py3_Jerry`.
  That worktree is for publish-only operations.
- The only acceptable edits in the `main` worktree are explicit user-requested
  publish/release actions or direct updates to this instruction file.
- `src/ldsc/` is the only supported public package surface. `src/ldsc/_kernel/`
  is private implementation code.
- Keep one CLI surface: `ldsc` with subcommands `annotate`, `build-ref-panel`,
  `ldscore`, `munge-sumstats`, `h2`, `partitioned-h2`, and `rg`.
- Public workflow modules own user-facing path resolution, header inference,
  global config handling, genome-build inference, and output preflight. Kernel
  modules receive resolved primitive inputs and do numerical work or low-level
  parsing.
- Keep alias and identifier normalization centralized in `column_inference.py`.
  Keep hg19/hg38 and 0-based/1-based inference centralized in
  `genome_build_inference.py`.
- Public LD-score output is a canonical result directory written by
  `LDScoreDirectoryWriter`: `manifest.json`, `baseline.parquet`, and optional
  `query.parquet`. Regression consumes that aggregated directory and must not
  recompute LD scores.
- Legacy LDSC formats remain compatibility boundaries, not the public LD-score
  writer layout: annotation workflows read/write `.annot(.gz)`, munging writes
  `.sumstats.gz`, and `_kernel` keeps low-level support for `.l2.ldscore(.gz)`,
  `.w.l2.ldscore(.gz)`, `.M`, and `.M_5_50`.
- Query annotations require explicit baseline annotations. The synthetic
  all-ones `base` annotation is only for ordinary unpartitioned LD-score
  generation.
- Preserve the original LDSC regression default of using `.M_5_50`-style
  common-SNP counts when available.
- Keep the package self-sufficient: no imports from sibling trees or
  repository-root wrappers, and no circular imports.

# Code Writing

- Write concise code. Minimize unnecessary validation and guard clauses. Avoid over-engineering -- no redundant design layers.
- Properly document non-obvious logic; skip comments that merely restate what the code says.
- Prefer existing workflow objects, config dataclasses, path-resolution helpers, and column-inference registries over one-off parsing or normalization.

# Principles

- Keep it simple.

- Be concise in responses.

- After each major change, commit with a meaningful message, or remind the user to commit with a suggested one-line message.

- Avoid AI tools name (like Codex, Claude, Grok, Gemini, ...) in code comments or git commit message (including authorship) or PR body.

- Record only significant and recurring mistakes. Ignore minor, one-off errors. When a meaningful pattern of errors is observed, document the correction and update this file to prevent future occurrences.

## Planning & Ambiguities

For complex features or significant refactors, do not begin coding immediately. Follow this workflow:

**Step 1.** Use `superpowers:brainstorming` to explore design and requirements. Ask clarifying questions -- and keep asking -- until the scope, approach, and every significant design choice is confirmed and fully understood. Resolve all ambiguities before moving on: you may suggest resolutions, but explain your reasoning and get explicit confirmation. Record each confirmed choice. For libraries you are unfamiliar with or not confident about, ask the user for guidance on correct usage.

**Step 2.** Use `superpowers:writing-plans` to produce the implementation plan. Surface any remaining ambiguities to the user before finalizing.

**Step 3.** Execute the plan with `superpowers:executing-plans`.

# Context Management

- Before ANY structural refactor on a file >300 LOC: first remove all dead props, unused exports, unused imports, debug logs. Commit cleanup separately. Dead code burns tokens that trigger compaction faster.
- Use sub-agents only when you have 2 or more large, genuinely independent tasks -- each touching a separate module with no shared dependencies. Default to sequential work; for a small project, parallelization rarely helps and makes progress harder to track.
- After 10+ messages, or when resuming after any gap: re-read the active superpowers plan in `docs/superpowers/` (if any) and any source files you intend to edit before proceeding. Do not rely on memory of their contents.
- Treat `docs/current/architecture.md`, `docs/current/class-and-features.md`, `docs/current/code_structure.md`, `docs/current/data_flow.md`, `docs/current/layer-structure.md`, `docs/current/`, and `docs/superpowers/` as the design source of truth over the current implementation when conflicts arise.

# Testing

Write tests before writing implementation code. Use `superpowers:test-driven-development` for the full red-green-refactor cycle. For numerical methods, the test must demonstrate correctness against a known input/output pair, a convergence check, or a comparison to a reference implementation -- not just that the code runs without error.

- Use the local `tests/` tree. The migrated tests import from `src/ldsc/`.
- Use `pytest` as the primary test command. Keep the standard-library unittest command working while the transition to pytest remains in progress.

# Debugging

When encountering a bug or unexpected behavior, use `superpowers:systematic-debugging`. Always find the root cause before proposing a fix. If 3 fix attempts fail, stop and ask the user -- repeated failures indicate a design problem, not a code problem.

# Verification

Before claiming any work is complete or a bug is fixed, run the test suite and show the actual output. Never say "it should work" or "tests should pass" -- run the command and include the result in your response.

## Documentation

Use `my-skills:fun-doc` when writing or updating docstrings for any Python function,
class, or module header. Trigger it any time a public function lacks a docstring or
an existing docstring is incomplete or outdated.

After any code modification, always update:
- **Docstrings** for any modified functions or classes.
- **Design documents** at `docs/current/` (if applicable).
- **Tutorials** at `tutorials/` (if applicable).

When initializing the project, create `design_map.md` to track the correspondence
between design documents and function/module implementations. After any major change,
update `design_map.md` to keep it aligned with the current implementation.

## Architecture

Use `my-skills:architecture-doc` to create or update `architecture.md` whenever
the module structure changes significantly, a new module is added, or a new
contributor needs an overview of the codebase.

## Explaining Code

When explaining a workflow, processing pipeline, or mechanism, anchor each statement
to the relevant source location (file path and line number, or function name).

## Research and Document Reading

When reading papers or documents, back every major statement or conclusion with a
citation to the source location (e.g., "Section 3.2", "Table 1", "Appendix A").

# Skill Usage Policy

`superpowers:using-superpowers` should only be invoked for non-trivial tasks. Individual superpowers skills may be used normally when relevant, including: `superpowers:brainstorming`, `superpowers:writing-plans`, `superpowers:executing-plans`, `superpowers:systematic-debugging`, `superpowers:verification-before-completion`, `superpowers:test-driven-development`.

No need to invoke skills for trivial tasks (simple typo fixes, obvious one-liners).

This policy takes highest priority and overrides any default skill invocation behavior.

# Resources

- code repos: `/Users/wenbinwu/Documents_local/Research/SullivanLab/LDSC/repos/ldsc_py3_Jerry_workspace/ldsc_py3_restructured`
- public package: `src/ldsc`
- internal kernel: `src/ldsc/_kernel`
- documents: `docs/current/architecture.md`, `docs/current/class-and-features.md`, `docs/current/code_structure.md`, `docs/current/data_flow.md`, `docs/current/layer-structure.md`, `docs/current/`
- tutorials: `tutorials/`
- prior plans: `docs/superpowers/`
- tests: `tests/`
