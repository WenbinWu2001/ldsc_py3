# PLANS.md

This file tracks active and future work for the refactored package.

## Completed: Restructure Into `src/ldsc/`

Status: complete

Summary:

- moved the public package surface to `src/ldsc/`
- moved internal compute code to `src/ldsc/_kernel/`
- replaced the old split script surface with one `ldsc` CLI and subcommands
- removed the old root wrappers and the old `refactor/ldscore/` package
- kept the refactor tree self-sufficient with no imports from `legacy/`
- kept the local parity-oriented test suite green

Resulting command surface:

- `ldsc annotate`
- `ldsc ldscore`
- `ldsc munge-sumstats`
- `ldsc h2`
- `ldsc partitioned-h2`
- `ldsc rg`

Resulting package split:

- public modules in `src/ldsc/`
- internal compute kernels in `src/ldsc/_kernel/`

## Current Follow-Up Work

### 1. Stabilize the public CLI

Goals:

- keep subcommand argument names clean and consistent
- reduce leftover legacy argument naming where possible
- add direct smoke tests for installed `ldsc` console entry behavior

Suggested checks:

- `python -m ldsc.cli --help`
- per-subcommand `--help`
- output-path behavior for flat vs per-chrom layouts

### 2. Tighten kernel cleanup

Goals:

- remove any compatibility aliases that only exist for transition support
- trim dead helpers inside `_kernel` once parity coverage proves they are unused
- address minor warnings such as unclosed file handles and `np.linalg.lstsq` future warnings

### 3. Improve tutorials

Goals:

- keep markdown and notebook tutorials aligned with `ldsc` / `ldsc.*`
- add one end-to-end example per main workflow
- keep tutorial outputs separate from tracked fixture data

### 4. Extend the output layer

Goals:

- add richer post-processing producers
- support cleaner query-only partitioned-h2 summaries by default
- leave room for future plotting/report producers without changing workflow modules

## Verification Baseline

Current verification command:

```bash
python -m unittest discover -s tests -p 'test*.py' -v
```

Acceptance standard for future changes:

- keep public imports under `ldsc`
- do not reintroduce imports from `legacy/`
- preserve LDSC-compatible file outputs unless a planned format change is explicit
