# Build Commands

```bash
cd /Users/wenbinwu/Documents_local/Research/SullivanLab/LDSC/repos/ldsc_py3_Jerry_workspace/ldsc_py3_restructured
conda env create -f environment.yml -n ldsc3-dev    # first-time development environment setup
source /Users/wenbinwu/miniforge3/etc/profile.d/conda.sh
conda activate ldsc3-dev
python -m pip install -e ".[dev]"                   # install editable package with dev extras
ldsc --help                                         # installed console entry point
python -m ldsc --help                               # module entry point
pytest                                              # primary test suite
python -m unittest discover -s tests -p 'test*.py' -v
```

For non-interactive agent commands, activate the project environment with:

```bash
source /Users/wenbinwu/miniforge3/etc/profile.d/conda.sh && conda activate ldsc3-dev && <command>
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

- Keep it simple. The ultimate goal is to deliver the best user experience.

- Be concise in responses.

- After each major change, commit with a meaningful message, or remind the user to commit with a suggested one-line message. Use Conventional Commits: `<type>(<scope>): <description>`. Subject <= 50 chars; body lines <= 72 chars. Body explains *what* and *why*, not *how*. Footer for issue refs (`Fixes #123`) or breaking changes. Common types: `feat`, `fix`, `docs`, `style`, `refactor`, `test`, `chore`.

- Avoid AI tools name (like Codex, Claude, Grok, Gemini, ...) in code comments or git commit message (including authorship) or PR body.

- When a significant or recurring mistake occurs (same class of error seen 2+ times, or a non-obvious bug whose fix required real investigation), append an entry to `lessons.md` at the repo root. Each entry: one-line summary, root cause, and the correction. Skip one-off typos and trivial slips. Create `lessons.md` on the first such entry if it does not exist yet.

## Planning & Ambiguities

For anything more than a several-line modification, do NOT begin coding immediately:

- **Ask clarifying questions first, and keep asking** until every ambiguity is resolved. Do not start coding while any non-trivial question is open.
- For medium-scope changes (multiple functions, new code paths, or more than one file), prompt the user to turn on plan mode before editing.
- For libraries you are not confident about, ask the user for guidance on correct usage before relying on them.

A "several-line modification" means a localized edit (rename, small bug fix, obvious one-liner). Anything larger falls under the rules above.

### Complex features or significant refactors

For complex features or significant refactors, use the following three steps with `superpowers`. This workflow produces **two separate documents**: a design doc (user-facing specs and concepts) and an implementation plan (codebase nitty-gritties for execution).

**Step 1.** Use `superpowers:brainstorming` to explore design and requirements. Resolve all ambiguities and get explicit user confirmation on every significant design choice before moving on.

**Step 2.** Produce both documents, each named with a date prefix (`YYYY-MM-DD`) plus the topic, sharing the same prefix and topic but differing in suffix and directory:

- **Design doc**: `docs/superpowers/specs/<YYYY-MM-DD>-<topic>-design.md`
  Specs and conceptual, user-facing content. Captures *what* and *why*, not *how*.
- **Implementation plan**: `docs/superpowers/plans/<YYYY-MM-DD>-<topic>-plan.md`
  Codebase nitty-gritties for execution only. Use `superpowers:writing-plans`.

Surface any remaining ambiguities and obtain explicit user approval of both documents before proceeding.

**Step 3.** Execute the implementation plan with `superpowers:executing-plans`.

# Context Management

- Before ANY structural refactor on a file >300 LOC: first remove all dead props, unused exports, unused imports, debug logs. Commit cleanup separately. Dead code burns tokens that trigger compaction faster.
- Use sub-agents only when you have 2 or more large, genuinely independent tasks -- each touching a separate module with no shared dependencies. Default to sequential work; for a small project, parallelization rarely helps and makes progress harder to track.
- After 10+ messages, or when resuming after any gap: re-read the active plan document under `docs/superpowers/` (if any), `lessons.md` (if it exists), and any source files you intend to edit before proceeding. Do not rely on memory of their contents.
- Treat `docs/current/architecture.md`, `docs/current/class-and-features.md`, `docs/current/code-structure.md`, `docs/current/data-flow.md`, `docs/current/layer-structure.md`, `docs/current/`, and `docs/superpowers/` as the design source of truth over the current implementation when conflicts arise.

# Testing

Write tests before writing implementation code. Use `superpowers:test-driven-development` for the full red-green-refactor cycle. For numerical methods, the test must demonstrate correctness against a known input/output pair, a convergence check, or a comparison to a reference implementation -- not just that the code runs without error.

- Use the local `tests/` tree. The migrated tests import from `src/ldsc/`.
- Use `pytest` as the primary test command. Keep the standard-library unittest command working while the transition to pytest remains in progress.

# Debugging

When encountering a bug or unexpected behavior, use `superpowers:systematic-debugging`. Always find the root cause before proposing a fix. Before proposing a fix, check `lessons.md` (if it exists) for prior occurrences of a similar bug. If 3 fix attempts fail, stop and ask the user -- repeated failures indicate a design problem, not a code problem.

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

## Citations

Anchor every claim to its source:

- **Code**: when explaining a workflow, pipeline, or mechanism, cite the file path and line number (or function name).
- **Papers and documents**: when summarizing or quoting, cite the section, table, or appendix (e.g., "Section 3.2", "Table 1", "Appendix A").

# Skill Usage Policy

`superpowers:using-superpowers` should only be invoked for non-trivial tasks. Individual superpowers skills may be used normally when relevant, including: `superpowers:brainstorming`, `superpowers:writing-plans`, `superpowers:executing-plans`, `superpowers:systematic-debugging`, `superpowers:verification-before-completion`, `superpowers:test-driven-development`.

No need to invoke skills for trivial tasks (simple typo fixes, obvious one-liners). 

This policy takes highest priority and overrides any default skill invocation behavior.

# Resources

- code repos: `/Users/wenbinwu/Documents_local/Research/SullivanLab/LDSC/repos/ldsc_py3_Jerry_workspace/ldsc_py3_restructured`
- public package: `src/ldsc`
- internal kernel: `src/ldsc/_kernel`
- documents: `docs/current/architecture.md`, `docs/current/class-and-features.md`, `docs/current/code-structure.md`, `docs/current/data-flow.md`, `docs/current/layer-structure.md`, `docs/current/`
- tutorials: `tutorials/`
- prior plans: `docs/superpowers/`
- tests: `tests/`
