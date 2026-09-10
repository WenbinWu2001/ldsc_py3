# Threshold Comparison Audit

Last updated on: 2026-09-10

Historical audit, relocated from the repository root. Statements and source line numbers below describe the earlier threshold audit, not the current implementation; see [context and current boundaries](README.md).

## Goal

Audit threshold-based filtering/comparison semantics across functionality modules and identify whether each threshold uses strict (`>`, `<`) or inclusive (`>=`, `<=`) comparisons.

## Phases

1. **Complete** - Identify functionality modules and relevant CLI/config threshold names.
2. **Complete** - Search implementation code for threshold comparisons and masks.
3. **Complete** - Cross-check tests/docs for intended behavior.
4. **Complete** - Summarize mixed semantics and harmonization candidates.

## Errors Encountered

None yet.
