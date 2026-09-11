# Regression intercept consolidation

Last updated on: 2026-09-11

Implemented the approved [intercept audit recommendation](regression-intercepts.md): `h2`, `partitioned-h2`, and `rg` reject `--no-intercept`. No alias or replacement boolean is retained. The [LDSC2 flag map](../../current/legacy-cli-flag-map.md#regression-intercept-consolidation) records all three replacements and their omission behavior.

| Command | Standard fixed-intercept request |
| --- | --- |
| `h2` | `--intercept-h2 1` |
| `partitioned-h2` | `--intercept-h2 1` |
| `rg` | `--intercept-h2 1 --intercept-gencov 0` |

The regression argument builders and CLI config conversion now use the numeric flags exclusively. The obsolete shortcut-conflict validator and its call sites were removed. Help explains the fixed values, scalar broadcasting in `rg`, omission defaults, and unchanged two-step restrictions. See [`regression_runner._add_common_regression_arguments` / `_runner_from_args`](../../../src/ldsc/regression_runner.py).

The numerical kernel, estimator selection, default cutoff 30, and public `RegressionConfig.use_intercept` behavior are unchanged. In particular, single-annotation covariance-only `rg` still fails when automatic two-step estimation conflicts with a fixed covariance intercept, as explicitly requested. Omitting all intercept flags still estimates the intercepts using the existing defaults.

Updated the current regression reference, IO inventory, relevant wiki explanation, and three Markdown tutorials. Both affected notebook tutorials now use the explicit numeric settings in their CLI and Python examples. Historical audit inventories and dated plans remain historical; their current status is linked from the follow-up audit.

## Verification

Three rejection tests first failed because the old flag was accepted, then passed after removal. Numerical tests compare the real CLI-to-config path against the retained public Python standard-fixed-intercept behavior on deterministic single-annotation h2, multi-annotation h2, and rg fixtures. They check coefficients, covariance matrices, totals, standard errors, jackknife values, and rg statistics, as well as default free fits and the retained two-step rejections. See [`test_cli_help.py`](../../../tests/test_cli_help.py) and [`test_regression_fit_outcomes.py`](../../../tests/test_regression_fit_outcomes.py).

Final verification:

- Full test suite: **1,547 passed, one skipped, 132 subtests passed**, with 189 warnings in 95.03 seconds. The skip applies when the optional `pyarrow` dependency is installed.
- Both affected notebooks completed all five code cells in order, including their real CLI subprocesses. The small rg example emitted a numerical `RuntimeWarning` from the jackknife square root; its smoke run still completed.
- Six real help invocations and six retired-flag rejection invocations passed across `ldsc` and `python -m ldsc`. Help text matched between entry points, omitted the retired flag, and rejected invocations created no output directory.
- Compared 21 parser surfaces with the saved pre-consolidation snapshot: only the three `--no-intercept` entries were removed. All remaining names, destinations, defaults, argument counts, types, choices, required settings, and existing mutual-exclusion groups matched. The CLI now has 212 visible command-option entries.
- Edited Python modules and notebook cells compiled; 123 local documentation links/anchors and dates in 11 documents validated; `git diff --check` passed.

Local verification logs: `/tmp/ldsc-intercept-consolidation-pytest.txt`, `/tmp/ldsc-intercept-notebook-rg.txt`, and `/tmp/ldsc-intercept-notebook-partitioned.txt`. The parser snapshot is `/tmp/ldsc-intercept-consolidation-parser.json`.
