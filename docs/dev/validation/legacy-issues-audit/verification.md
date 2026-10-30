# Verification record and limits

Last updated on: 2026-10-28

This record distinguishes source inspection, executed current behavior, and unavailable historical reproductions.
The issue reports remain the source of their original symptoms; a synthetic boundary does not establish every historical dataset's root cause.
The work used local `restructure` at `9d095e66a8fddff02831c644d286ed274876503b`.

## Runtime and focused tests

The interpreter was `/Users/wenbinwu/miniforge3/envs/ldsc3-dev/bin/python`, Python 3.13.13.
Its imported `ldsc` came from this checkout's `src/ldsc/__init__.py`.
Relevant versions were NumPy 2.4.5, pandas 2.3.3, SciPy 1.17.1 and PyArrow 23.0.1.
This was an existing editable development environment, not a fresh installation.
The console command and `python -m ldsc --help` both succeeded and produced identical help text.

From the package repository root, the following focused command passed:

```bash
python -m pytest -q \
  tests/package/test_package_layout.py \
  tests/package/test_package_dependencies.py \
  tests/package/test_cli_help.py \
  tests/annotation/test_annotation_format.py \
  tests/annotation/test_annotation_semantics.py \
  tests/annotation/test_annotate_bed_validation.py \
  tests/annotation/test_annotate_grid.py \
  tests/annotation/test_annotate_preflight.py \
  tests/sumstats/test_sumstats_munger.py \
  tests/sumstats/test_munge_flag_contracts.py \
  tests/sumstats/test_munge_hm3_defaults.py \
  tests/ldscore/test_plink_exact_windows.py \
  tests/ldscore/test_ldscore_projection.py \
  tests/ldscore/test_reference_snp_funnel.py \
  tests/ldscore/test_ldscore_query_batches.py \
  tests/regression/test_kernel_regression.py \
  tests/regression/test_jackknife.py \
  tests/regression/test_overlap_matrix.py \
  tests/regression/test_prevalence.py \
  tests/regression/test_regression_fit_outcomes.py \
  tests/regression/test_h2_diagnostics.py \
  tests/regression/test_partitioned_query_failures.py \
  tests/regression/test_quantile_h2.py
```

Result: **499 passed, 16 warnings, 36 subtests passed in 139.42 seconds**.
There were no reported skips or expected failures in this selected run.
Warnings included the known genome-build setting ignored in rsID-family test fixtures.
This checks the selected package contracts, not all original datasets, every environment, or equivalence to legacy LDSC.
No source tests or immutable golden fixtures were added, edited or regenerated.

## V1: Numeric invariants on direct legacy input

**Executed entry path:** `load_sumstats()` → `RegressionRunner.build_dataset()` → `estimate_h2()`.
The summary statistics were real temporary legacy-format text files; the LD-score input was a synthetic in-memory `LDScoreResult` with 200 SNPs.
All SNPs were on chromosome 22 with A/C alleles, finite positive LD scores and an rsID SNP identifier mode.
NumPy generator seed 128 produced reference weights in [1,3), baseline scores in [1,8), and signed Z values from a normal distribution with SD 2.
Ordinary N was 10,000; one input row was replaced in each independent scenario.
The fit used 10 jackknife blocks.

| Altered field | Loader and dataset result | Fit result |
| --- | --- | --- |
| `Z=inf` | Accepted; all 200 rows retained, including infinity. | Late `LinAlgError`: SVD did not converge. |
| `N=inf` | Accepted; all 200 rows retained, including infinity. | Late `LinAlgError`: SVD did not converge. |
| `N=0` | Accepted; invalid row retained. | Returned a total estimate of approximately 0.014461. |
| `N=-1` | Accepted; invalid row retained. | Returned a total estimate of approximately -0.009870. |

These values characterize the fixture, not valid scientific estimates.
The confirmed gap is missing early finite/positive-N validation on this direct-import path.
The probe did not exercise a complete on-disk LD-score directory or all CLI publication steps.
It does not establish that every raw munging path accepts the same invalid inputs.

The low-level annotation numeric helper also accepted infinity, but the subsequent `classify_annotation_values()` rejects nonfinite values in annotation preparation.
Consequently that helper-only observation is **not** reported as a confirmed annotation-pipeline vulnerability.

**Canonical repair proposal:** validate finite statistics and strictly positive finite sample sizes at the shared accepted-input boundary, with row-level diagnostics and no imputation.
Future regression cases should include positive/negative infinity, zero/negative N, both direct text and native artifacts, and unchanged valid signed statistics.
Sources: [`munger.py`, `load_sumstats()`](../../../../src/ldsc/_sumstats/munger.py), [`runner.py`, dataset construction](../../../../src/ldsc/_regression/runner.py), [`semantics.py`](../../../../src/ldsc/_annotation/semantics.py).

## V2: Genetic-correlation delete-block boundary

**Executed paths:** the real `RG` kernel, then current `_fit_rg_dataset()` and `_summarize_rg_pair()` on an `RGRegressionDataset`.
No estimator was mocked.
This tests current fitting and summary construction; it is not a complete filesystem CLI run.

The deterministic fixture has 80 rows, a single LD-score column $x=\operatorname{linspace}(1,10,80)$, unit regression weights, $N_1=N_2=10000$, $M=1000$, and 10 blocks.
Using `numpy.random.default_rng(0)`, draw an 80-by-1 normal array for the first trait and then another for the second, both with SD 0.8.
Construct $Z_1=\sqrt{\max(0.01,1+0.015x+\epsilon_1)}$ and $Z_2=\sqrt{\max(0.01,1+0.02x+\epsilon_2)}$ elementwise.
The fixture deliberately probes weak-signal ratio arithmetic; it is not a calibration simulation of a GWAS population.

The direct kernel returned positive whole-data heritabilities of approximately 0.00623769 and 0.0000184984.
Five delete-block heritability products were negative.
The direct correlation ratio was approximately 7.01164, with `rg_se=NaN` and an invalid-square-root warning.
The workflow fit applies its current default two-step setting and produced the following summary:

| Field | Observed value |
| --- | --- |
| `rg` | 7.3594863081624275 |
| `rg_se` | `NaN` |
| `p` | `NaN` |
| `status` | `ok` |
| `error` | Empty string |

An independent rerun reproduced these results.
The gap is the missing delete-block validity/uncertainty status, not merely that a finite estimate exceeds one.
Both-negative deleted heritabilities also require individual checks: a positive product does not make either heritability valid.

**Canonical repair proposal:** check each trait's whole-data and delete-block denominators and the resulting uncertainty; label undefined correlation inference explicitly not estimable and report the affected blocks.
Retain valid component heritability/covariance results when feasible.
Do not take absolute values, clip to parameter bounds, silently omit blocks, or substitute zero uncertainty.
Future cases should cover one-negative, both-negative, zero and near-zero deletions, plus valid low-signal fits.
Sources: [`estimators.py:746`](../../../../src/ldsc/_regression/estimators.py#L746), [`runner.py:1252`](../../../../src/ldsc/_regression/runner.py#L1252), [`runner.py:2475`](../../../../src/ldsc/_regression/runner.py#L2475).

## V3: Zero-total derived proportion

Using the same 80-row x, weights, N, M and block count as V2, pass an all-ones chi-square response to `Hsq(..., intercept=1)`.
The fit reaches a derived ratio/proportion with a zero total and raises `FloatingPointError: invalid value encountered in divide`.
The failure occurs before the caller can obtain a regular total-heritability summary.
An independent rerun confirmed it.
The original issue #345 concerns a related small-denominator failure, but its original dataset was unavailable; this boundary reproduction does not identify that dataset's exact cause.

**Canonical repair proposal:** separate estimable totals/coefficients from undefined ratios, and report unavailable derived fields with a reason.
An arbitrary epsilon changes the scientific quantity and is not a canonical repair.
Future cases should exercise exact zero, nonzero totals with a zero deleted total, and sufficiently small finite totals without relaxing numerical tolerances.
Sources: [`estimators.py`, `_prop()` and `_enrichment()`](../../../../src/ldsc/_regression/estimators.py#L410), [`jackknife.py`, ratio deletion](../../../../src/ldsc/_regression/jackknife.py#L599).

## V4: Binary projection and header delimiters

For annotation SNP positions `[1,5,10,15]` and zero-based half-open BED intervals `[(0,10),(4,15)]`, `ChromosomeProjector.project()` returned `[1,1,1,1]`.
Per-SNP overlap multiplicities would instead be `[1,2,2,1]`.
This confirms intentional binary membership and an unsupported count-generation use case, rather than an incorrect binary result.
Sources: [`interval_projection.py:15`](../../../../src/ldsc/_annotation/interval_projection.py#L15), [`legacy-flag-map.md`](../../../reference/legacy-flag-map.md).

A second header case was reproduced through both CLI preview and full munging (issue #227).
The complete one-row fixture is:

```text
"SNP" "A1" "A2" "P" "N" "Z"
"rs1" "A" "C" 0.05 1000 -1.96
```

Using `munge-sumstats --snp-identifier rsid --no-snp-restriction --raw-sumstats-file INPUT --output-dir OUTPUT --infer-only` reported missing N and signed direction.
Removing `--infer-only` failed to recognize the signed-statistic header.
The source splits the raw quoted header without applying the body reader's quote handling.
This is separate from the original thread's incorrect nonzero signed-Z null argument.
The canonical repair is consistent quote-aware header/body parsing, with a focused quoted/unquoted round-trip case.

The header reader returned one combined token for `rsid,bp,ref,alt,beta,pval,N` in a CSV file and seven tokens for the equivalent tab-delimited header.
This reproduces the CSV discoverability/input-format trap discussed in issue #374.
It does not establish a regression in an advertised CSV contract.
The small remedy is an actionable delimiter diagnosis or clearly documented conversion to the accepted format.
Source: [`munge_input.py`, `read_header()`](../../../../src/ldsc/_sumstats/munge_input.py#L97).

## Source, link and peer checks

- Enumerated issue IDs reconcile with independent GitHub search totals of 458/326/132, with `incomplete_results=false` for all three queries.
  Pull requests were excluded, full accessible comments were attached by issue ID, and all 28 FAQ boundaries were checked against raw Markdown.
- Legacy fix claims were checked against thread/event evidence and specific commit/PR state where available.
  An open or closed-but-unmerged PR was not promoted to a merged upstream fix.
- Public HEAD requests checked 57 distinct reference-resource URLs extracted from issue text.
  Results: 23 HTTP 404, 23 HTTP 200, five HTTP 403, one HTTP 400 and five unresolved network/trust failures.
  Six of the HTTP 200 responses are Google console redirects to sign-in.
  These are access barriers, not verified reference-file downloads.
  The [resource report](resources.md) lists every checked URL and the limit of its result.
- The original cross-trait and continuous-annotation papers were checked in local PDF form as well as primary web sources.
  Rendered PDF page 2 of the cross-trait paper and page 8 of the continuous-annotation paper were visually checked for the formulas cited in the methodology report.
- A separate reviewer checked the high-severity numerical claims, reran V2/V3, reviewed the FAQ assessments, and sampled source-based prevention claims across the other issue batches.
  The review corrected annotation-versus-reference alignment wording and missing evidence paths.

## Final documentation checks

Final documentation checks found exactly 458 unique issue assessment anchors and 28 unique FAQ assessment anchors, with no missing or duplicate primary assignments.
An independent parse of the rendered workflow-tag lines reproduced all seven frequency totals.
All seven category summary tables contained one row per individual assessment in nondecreasing severity rank.
All local Markdown targets and anchors resolved, code-line references were within their files, external link targets contained no whitespace, and no template placeholders remained.
The final Git scope contained only this report folder and the validation index; `git diff --check` passed.

## What was not established

No full historical environment matrix, fresh installation trial, Windows/Linux compatibility run, complete legacy numerical-equivalence study, performance benchmark at production scale, or new estimator calibration was performed.
No large resource archives were downloaded or unpacked.
A reachable archive listing does not verify its internal sample IDs, SNP compatibility, biological suitability or reuse terms.
Four API-reported comments and some historical attachments/input files were unavailable; their absence remains attached to the affected conclusions.
Raw source captures, scripts and logs are retained outside the repository at `/private/tmp/ldsc-legacy-issues-audit-20261028/`; durable conclusions and fixture definitions are recorded here.
