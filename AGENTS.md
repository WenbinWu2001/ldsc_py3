# Build Commands

```bash
cd /Users/wenbinwu/Documents_local/Research/SullivanLab/LDSC/repos/ldsc_py3_Jerry_workspace/ldsc_py3_restructured
conda env create -f environment.yml
conda activate ldsc-py3-Jerry
pip install -e .
python -m ldsc.cli -h
python -m unittest discover -s tests -p 'test*.py' -v
```

# Key Architectural Invariants

- Treat `src/ldsc/` as the only public Python package surface.
- Treat `src/ldsc/_kernel/` as internal-only compute and file-format code.
- Keep the public workflow modules narrow and obvious:
  - `ldsc.annotation_builder`
  - `ldsc.ldscore_calculator`
  - `ldsc.sumstats_munger`
  - `ldsc.regression_runner`
  - `ldsc.outputs`
- Keep one CLI surface: `ldsc` with subcommands `annotate`, `ldscore`, `munge-sumstats`, `h2`, `partitioned-h2`, and `rg`.
- Keep the package self-sufficient. Do not add imports that depend on sibling trees or repository-root wrappers.
- Preserve existing file-format contracts for `.annot`, `.sumstats.gz`, `.l2.ldscore(.gz)`, `.w.l2.ldscore(.gz)`, `.M`, and `.M_5_50` unless a coordinated format change is intentional.
- LD-score computation remains chromosome-wise. Regression consumes only the aggregated cross-chromosome result.
- Preserve the original LDSC regression default of using `.M_5_50`-style common-SNP counts when available.

# Principles

- keep it simple
- ask if a specification, algorithm choice, or user-facing behavior is unclear
- do not mention AI tool names in code comments, commits, or PR text
- record only significant, repeated project mistakes in this file

# ExecPlans

For complex features or significant refactors, write or update an ExecPlan in `PLANS.md`.

# Dealing with Ambiguities

- Ask the user about ambiguous feature behavior, algorithm choices, and external interface decisions.
- Resolve local implementation details directly when the design docs already constrain them.
- If a workflow depends on a library or file format whose behavior is uncertain, stop and ask rather than guessing.

# Context Management

- Before editing any file larger than 300 LOC, re-read the file and the relevant part of `PLANS.md`.
- Treat `architecture.md`, `class-and-features.md`, `code_structure.md`, and `PLANS.md` as the design source of truth over the current implementation when conflicts arise.
- Use sub-agents only for large independent tasks with clearly separate write scopes.

# Testing

- Write tests before implementation changes when adding behavior.
- For numerical methods, verify against known outputs, legacy parity, or explicit invariants.
- Use the local `tests/` tree. The migrated tests import from `src/ldsc/`.

# Debugging

- Find the root cause before changing code.
- If three fix attempts fail, stop and ask the user. Repeated failure usually means the design assumption is wrong.

# Verification

- Before claiming work is complete, run the relevant tests and report the actual result.
- Do not say work is done based only on reasoning.

# Documentation

- Use `fun-doc` when writing or updating Python docstrings for public functions, classes, or modules.
- Use `architecture-doc` whenever the package structure changes significantly or `architecture.md` is out of sync with the code.


# Resources

- code repo: /Users/wenbinwu/Documents_local/Research/SullivanLab/LDSC/repos/ldsc_py3_Jerry_workspace/ldsc_py3_restructured
- public package: /Users/wenbinwu/Documents_local/Research/SullivanLab/LDSC/repos/ldsc_py3_Jerry_workspace/ldsc_py3_restructured/src/ldsc
- internal kernel: /Users/wenbinwu/Documents_local/Research/SullivanLab/LDSC/repos/ldsc_py3_Jerry_workspace/ldsc_py3_restructured/src/ldsc/_kernel
- design docs: /Users/wenbinwu/Documents_local/Research/SullivanLab/LDSC/repos/ldsc_py3_Jerry_workspace/ldsc_py3_restructured/architecture.md, /Users/wenbinwu/Documents_local/Research/SullivanLab/LDSC/repos/ldsc_py3_Jerry_workspace/ldsc_py3_restructured/class-and-features.md, /Users/wenbinwu/Documents_local/Research/SullivanLab/LDSC/repos/ldsc_py3_Jerry_workspace/ldsc_py3_restructured/code_structure.md, /Users/wenbinwu/Documents_local/Research/SullivanLab/LDSC/repos/ldsc_py3_Jerry_workspace/ldsc_py3_restructured/PLANS.md
- tests: /Users/wenbinwu/Documents_local/Research/SullivanLab/LDSC/repos/ldsc_py3_Jerry_workspace/ldsc_py3_restructured/tests
