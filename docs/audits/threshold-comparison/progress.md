# Threshold Comparison Audit Progress

Last updated on: 2026-09-10

Historical audit, relocated from the repository root. Statements and source line numbers below describe the earlier threshold audit, not the current implementation; see [context and current boundaries](README.md).

- Initialized audit planning files for `ldsc_py3_restructured` on branch `restructure`.
- Confirmed active branch is `restructure`.
- Listed repo files and ran a broad threshold/comparison search across Python sources.
- Ran an AST comparison extraction over `src/ldsc`.
- Inspected threshold implementation sites in sumstats munging, LD-score calculation, reference-panel loading/building, overlap, regression, genome-build inference, and config validation.
- Cross-checked threshold behavior against targeted tests and current docs.
- Added the final threshold matrix and doc/test evidence to `findings.md`.
