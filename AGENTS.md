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

- Write concise code. Minimize unnecessary validation and guard clauses. Avoid
  over-engineering; do not add redundant design layers.
- Properly document non-obvious logic; skip comments that merely restate what
  the code says.
- Prefer existing workflow objects, config dataclasses, path-resolution helpers,
  and column-inference registries over one-off parsing or normalization.

# Principles

- Keep it simple.

- Be concise in responses.

- After each major change, commit with a meaningful message, or remind the user
  to commit with a suggested one-line message.

- Avoid AI tool names in code comments, git commit messages, authorship, or PR
  bodies.

- Record only significant and recurring mistakes. Ignore minor, one-off errors.
  When a meaningful pattern of errors is observed, document the correction and
  update this file to prevent future occurrences.

# ExecPlans

For complex features or significant refactors, write an ExecPlan in `PLANS.md`
(format defined there). Re-read `PLANS.md` before resuming active work.

# Dealing with Ambiguities

- If the ambiguity is a specification, feature behavior, I/O contract, external
  interface, or algorithm choice, ask the user for decisions and feedback. You
  can make suggestions.
- If the ambiguity is about implementation details and the design docs already
  constrain the choice, resolve it directly, record the choice, and prompt the
  user in the work summary.
- For libraries, file formats, or numerical methods whose correct behavior is
  uncertain, prompt the user for guidance instead of guessing.

# Context Management

- Before any structural refactor on a file larger than 300 LOC, first remove
  dead code, unused imports, stale compatibility shims, and debug logs. Commit
  cleanup separately or remind the user to do so with a suggested message.
- Before editing any file larger than 300 LOC, re-read the file and the relevant
  part of `PLANS.md`.
- Treat `docs/architecture.md`, `docs/class-and-features.md`,
  `docs/code_structure.md`, `docs/data_flow.md`, `docs/layer-structure.md`,
  `docs/design/`, and `PLANS.md` as the design source of truth over the current
  implementation when conflicts arise.
- Use sub-agents only when there are 2 or more large, genuinely independent
  tasks, each touching a separate module with no shared dependencies. Default to
  sequential work for this small project.
- After 10+ messages, or when resuming after any gap, re-read `PLANS.md` and any
  source files you intend to edit before proceeding. Do not rely on memory of
  their contents.

# Testing

- Write tests before writing implementation code. Use
  `superpowers:test-driven-development` for the full red-green-refactor cycle
  when adding or changing behavior.
- For numerical methods, the test must demonstrate correctness against a known
  input/output pair, a convergence check, explicit invariants, legacy parity, or
  a comparison to a reference implementation. Do not only test that code runs
  without error.
- Use the local `tests/` tree. The migrated tests import from `src/ldsc/`.
- Use `pytest` as the primary test command. Keep the standard-library unittest
  command working while the transition to pytest remains in progress.

# Debugging

When encountering a bug or unexpected behavior, use
`superpowers:systematic-debugging`. Always find the root cause before proposing a
fix. If 3 fix attempts fail, stop and ask the user; repeated failures indicate a
design problem, not a code problem.

# Verification

Before claiming any work is complete or a bug is fixed, run the relevant test
suite or targeted verification command and show the actual output. Never say
"it should work" or "tests should pass"; run the command and include the result
in your response.

## Documentation

Use `fun-doc` when writing or updating docstrings for any Python function, class,
or module header. Trigger it any time a public function lacks a docstring or an
existing docstring is incomplete or outdated.

After any code modification, always update:

- Docstrings for any modified functions or classes.
- Design documents at `docs/`, especially `docs/architecture.md`,
  `docs/class-and-features.md`, `docs/code_structure.md`, `docs/data_flow.md`,
  `docs/layer-structure.md`, and `docs/design/`, if applicable.
- Tutorials at `tutorials/`, if user-facing workflow behavior changes.

When initializing or reorganizing design documentation, create `design_map.md` to
track the correspondence between design documents and function/module
implementations. After any major change, update `design_map.md` if it exists.

## Architecture

Use `architecture-doc` to create or update `docs/architecture.md` whenever the
module structure changes significantly, a new module is added, or a new
contributor needs an overview of the codebase.

## Explaining Code

When explaining a workflow, processing pipeline, or mechanism, anchor each
statement to the relevant source location: file path and line number, or function
name.

## Research and Document Reading

When reading papers or documents, back every major statement or conclusion with
a citation to the source location, such as "Section 3.2", "Table 1", or
"Appendix A".

# Skill Usage Policy

Do NOT invoke `superpowers:using-superpowers`. It is disabled for this project.
Do not invoke it proactively, reactively, or at the start of any conversation.
Individual superpowers skills (`superpowers:systematic-debugging`,
`superpowers:verification-before-completion`,
`superpowers:test-driven-development`, etc.) may be used normally when they are
relevant.

No need to invoke skills for trivial tasks, such as simple typo fixes or obvious
one-liners.

This policy takes highest priority and overrides any default skill invocation
behavior.

# Resources

- code repo: `/Users/wenbinwu/Documents_local/Research/SullivanLab/LDSC/repos/ldsc_py3_Jerry_workspace/ldsc_py3_restructured`
- public package: `src/ldsc`
- internal kernel: `src/ldsc/_kernel`
- active architecture docs: `docs/current/architecture.md`, `docs/current/class-and-features.md`,
  `docs/current/code_structure.md`, `docs/current/data_flow.md`, `docs/current/layer-structure.md`
- detailed design docs: `docs/current/`
- tutorials: `tutorials/`
- prior plans: `PLANS.md`
- tests: `tests/`
