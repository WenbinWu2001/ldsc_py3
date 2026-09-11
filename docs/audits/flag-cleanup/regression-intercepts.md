# Regression intercept flag audit

Last updated on: 2026-09-11

Decision and implementation: the user approved retiring the CLI shortcut and explicitly retaining the covariance-only single-annotation `rg` failure. The consolidation is now implemented; see the [current migration map](../../current/legacy-cli-flag-map.md#regression-intercept-consolidation) and [implementation record](intercept-consolidation.md). The audit below records the pre-consolidation interface and evidence. Its suggested changes to two-step policy and error handling were not included.

Scope: the current `h2`, `partitioned-h2`, and `rg` interfaces on the local `restructure` branch. This is an audit and recommendation; no flags, numerical defaults, or runtime behavior were changed. `quantile-h2` reads saved fitted coefficients and has no intercept-estimation flags; see [`quantile_h2.load_fitted_partitioned_model` / `add_quantile_h2_arguments`](../../../src/ldsc/quantile_h2.py).

## Recommendation

Retire `--no-intercept` in favor of the existing explicit fixed-value flags. Keep `--intercept-h2`, `--intercept-gencov`, and `--two-step-cutoff` separate. This removes one public flag from each regression command without removing an available fitted model. It also removes a misleading name: `--no-intercept` fixes the h2 intercept to 1 rather than making it zero.

| Current command | Current shortcut | Equivalent existing spelling |
| --- | --- | --- |
| `h2` | `--no-intercept` | `--intercept-h2 1` |
| `partitioned-h2` | `--no-intercept` | `--intercept-h2 1` |
| `rg` | `--no-intercept` | `--intercept-h2 1 --intercept-gencov 0` |

Omitting all intercept flags should continue to estimate the intercepts with the existing estimator defaults. No replacement boolean or combined tuple syntax is necessary. This recommendation supersedes the original audit's suggestion to rename the shortcut to `--fix-intercept-defaults`; see [original audit, “Other commands”](README.md#other-commands).

The equivalence follows from [`RegressionRunner._fit_h2_dataset`](../../../src/ldsc/regression_runner.py:824), [`RegressionRunner._fit_rg_dataset`](../../../src/ldsc/regression_runner.py:1063), and [`_select_intercept`](../../../src/ldsc/regression_runner.py:2522). The kernel receives the same fixed intercepts and no automatic two-step cutoff in both spellings. The effective per-pair metadata policies also agree through [`_intercept_policy`](../../../src/ldsc/regression_runner.py:2215). The original command text and raw configuration representation can still differ; this is numerical and effective-policy equivalence, not a claim that all logs are byte-identical.

## Flag inventory and merge decisions

All four options are declared in [`_add_common_regression_arguments`](../../../src/ldsc/regression_runner.py:2961).

| Flag | Commands | Effective behavior | Decision |
| --- | --- | --- | --- |
| `--no-intercept` | All three | Fix h2 intercepts to 1; also fix covariance intercepts to 0 in `rg`. Default: off. Conflicts with either explicit fixed-intercept flag. | Redundant shortcut; absorb into the existing numeric controls. |
| `--intercept-h2 VALUE` | All three | Fix the h2 intercept. In `rg`, the same scalar is passed to both trait-specific h2 fits in every pair. If omitted, estimate it unless the shortcut fixes it. | Keep. This is the complete fixed-h2-intercept control. |
| `--intercept-gencov VALUE` | `rg` | Fix the genetic-covariance intercept to one scalar shared by every pair. If omitted, estimate it unless the shortcut fixes it. | Keep separately from h2: it controls a different fitted response and can have a different value or estimation policy. |
| `--two-step-cutoff VALUE` | All three | Set the inclusive first-stage chi-square threshold for estimating intercepts. Requires a single retained annotation and free intercepts. Omission enables cutoff 30 for single-annotation fits with a free h2 intercept; otherwise no two-step fit. | Keep as an advanced estimator control; it does not specify a fixed intercept value. |

`rg` fits two h2 models using squared Z scores and one covariance model using the product of the two traits' Z scores. Its two intercept controls therefore represent separate model quantities. A single numeric `--intercept` could not express the current default fixed combination `(1, 0)` or keep one quantity free while fixing the other. A tuple or key-value syntax could encode them, but would only relocate the same choices into a more complicated value. Sources: [`_kernel.regression.RG.__init__`](../../../src/ldsc/_kernel/regression.py:830), [`Hsq.__init__`](../../../src/ldsc/_kernel/regression.py:471), and [`Gencov.__init__`](../../../src/ldsc/_kernel/regression.py:680).

Do not merge `--two-step-cutoff` with `--chisq-max`, either. The former selects SNPs for the first-stage intercept estimate, after which the slope uses the retained fit population. The latter removes SNPs from that population before fitting. Sources: [`LD_Score_Regression.__init__`](../../../src/ldsc/_kernel/regression.py:294) and [`RegressionRunner._fit_h2_dataset`](../../../src/ldsc/regression_runner.py:805).

## Interactions that require separate attention

1. **Fixing only the covariance intercept fails in a single-annotation `rg` fit.** With `--intercept-gencov 0` and no h2 override, [`_fit_rg_dataset`](../../../src/ldsc/regression_runner.py:1071) selects cutoff 30 because the h2 intercept is free. [`RG.__init__`](../../../src/ldsc/_kernel/regression.py:836) passes that cutoff to all three fits, and the covariance fit rejects the fixed-intercept/two-step combination. The same covariance-only setting succeeds in a multi-annotation fit, where there is no automatic two-step cutoff. This interaction is not evidence that the two intercept controls should be merged.
2. **There is no explicit way to disable automatic two-step estimation while leaving the h2 intercept free.** Omission means automatic selection, while zero and negative cutoffs are rejected by [`RegressionConfig.__post_init__`](../../../src/ldsc/config.py:1083). Changing that interface or assigning two-step settings separately to the h2 and covariance fits would be an additional estimator-policy decision. Fixing the h2 intercept merely to avoid the error changes the requested model and should not be an automatic workaround.
3. **The fixed-intercept/two-step error can suggest an impossible remedy.** When cutoff 30 was selected automatically, the numerical error says to remove the supplied cutoff even though none was supplied. An eventual fix should identify the effective automatic setting and reject incompatible choices before substantial fitting. See [`LD_Score_Regression.__init__`](../../../src/ldsc/_kernel/regression.py:294).
4. **Current `rg` intercept inputs are scalars, not per-trait lists.** The h2 value is shared by both traits and all pairs; the covariance value is shared by all pairs. This is clear in the float parser actions and kernel call. The wording “value(s)” in [IO argument inventory, `rg`](../../current/io-argument-inventory.md) overstates the CLI's current support. This is a documentation discrepancy, not a reason to collapse the controls.

The first two points remain separate from the proposed redundant-flag removal. Preserve their current behavior unless an additional policy change is approved.

## Compatibility and practical use

The shortcut is used in the CLI examples in [`cross-trait-genetic-correlation.ipynb`](../../../tutorials/cross-trait-genetic-correlation.ipynb) and [`cell-specific-ldsc.ipynb`](../../../tutorials/cell-specific-ldsc.ipynb). It is redundant, not demonstrably unused. Repository searches cannot establish usage in external scripts or by other users.

Legacy LDSC also implements `--no-intercept` by selecting h2=1 and, for `rg`, covariance=0. See [`ldsc.py`, intercept argument declarations](/Users/wenbinwu/Documents_local/Research/SullivanLab/LDSC/repos/ldsc_py2_Bulik_workspace/ldsc_py2_Bulik/ldsc.py:529), [`ldscore/sumstats.py`, `estimate_h2`](/Users/wenbinwu/Documents_local/Research/SullivanLab/LDSC/repos/ldsc_py2_Bulik_workspace/ldsc_py2_Bulik/ldscore/sumstats.py:321), and [`estimate_rg`](/Users/wenbinwu/Documents_local/Research/SullivanLab/LDSC/repos/ldsc_py2_Bulik_workspace/ldsc_py2_Bulik/ldscore/sumstats.py:391), alongside the current resolution functions cited above. Removing the shortcut preserves those numerical models through the explicit replacements but changes a legacy CLI spelling.

If removal is approved, update the shared regression parser and conversion, help dependencies, tutorials, [legacy flag map](../../current/legacy-cli-flag-map.md), [regression configuration reference](../../current/regression-configuration.md), and [IO argument inventory](../../current/io-argument-inventory.md). Keep the original audit as historical evidence. A compatibility alias would preserve old scripts but would not actually retire the accepted flag; that would be a deliberate compatibility choice rather than an automatic extra layer.

The same representational redundancy exists between public `RegressionConfig.use_intercept` and its nullable numeric intercept fields. Removing that Python field would change the public Python interface and requires a separate explicit compatibility decision; it is not required to simplify the CLI. See [`RegressionConfig`](../../../src/ldsc/config.py:1025).

## Verification

- Traced parser declarations, CLI-to-config conversion, effective-intercept resolution, automatic two-step selection, kernel validation, output policy metadata, current documentation, tutorials, and legacy shortcut implementation.
- Ran the existing regression workflow, fit-outcome, and CLI-help suites: **118 tests passed, three subtests passed, 34 warnings**.
- Ran three numerical equivalence probes using the real CLI parser, config conversion, and numerical fit methods on deterministic 60-SNP fixtures. Compared intercepts, coefficients, covariance matrices, total estimates, standard errors, and jackknife delete values exactly; `rg` additionally matched both correlation estimates, correlation standard error, p-value, z-score, and effective intercept-policy labels. All three shortcut/replacement pairs matched.
- Ran 12 additional interaction probes covering defaults, explicit conflicts, fixed intercepts with two-step estimation, multi-annotation rejection, covariance-only `rg` in both annotation regimes, and rejection of cutoff zero. All matched the behavior described above.

Probe inputs reuse [`tests/test_regression_fit_outcomes.py`, `h2_dataset`](../../../tests/test_regression_fit_outcomes.py:16), adding a non-collinear query column for multi-annotation fits and deterministic second-trait data for `rg`. The local reproduction script and results are `/tmp/ldsc-intercept-audit.py` and `/tmp/ldsc-intercept-audit-results.json`; the focused test log is `/tmp/ldsc-intercept-audit-tests.txt`. These temporary files are supporting evidence for this run, not package artifacts.
