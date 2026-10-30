# Legacy LDSC issues and FAQ audit

Last updated on: 2026-10-28

This audit assesses every publicly accessible issue in [bulik/ldsc](https://github.com/bulik/ldsc/issues), including closed issues, and every entry in its [FAQ](https://github.com/bulik/ldsc/wiki/FAQ), against LDSC3.
It is a prevention and user-experience review, not an assertion of numerical equivalence or authorization to implement its recommendations.
Only this report folder and the validation index were changed.

## Audited state and coverage

- LDSC3: local `restructure`, revision `9d095e66a8fddff02831c644d286ed274876503b`; initially clean.
- Legacy LDSC: upstream `master`, revision `2fdeeb3b44379408794154993dbd6101b8946b7e` (2026-01-16).
- Local legacy comparison checkout: `e32c4b484726b8d376bc85b49baa45b33e22733e`; its local documentation commit is not evidence of an upstream fix.
- Sources acquired on 2026-10-28: **458 issues, 326 open and 132 closed; 1,110 accessible comments; 28 FAQ entries**.
- The issue API reports 1,114 comments in total.
  Issues #43, #255, #290 and #357 each expose one fewer comment through both repository-wide and individual comment retrieval.
  All returned comments were reviewed; the four unavailable comments are an explicit evidence gap.
- Issue collection covered six REST pages, including pull-request records that were excluded, followed by an empty page.
  Repository comments covered twelve populated pages and an empty thirteenth page; 668 issue events were retained after excluding pull requests.
- FAQ identifiers `FAQ-01` through `FAQ-28` follow source order, across General, LD Score Estimation, and Genetic Correlation and Heritability.
  The source page reports its last edit as 2017-04-22.
  FAQ subjects in this report are paraphrased; source order identifies the original question.

The category totals below reconcile to 458 issues plus 28 FAQ entries.
Each item has one primary category and its own assessment.
Repeated reports retain individual records; recommendations for their shared mechanism are consolidated.
Full issue bodies, all accessible comment text, later topic changes and relevant resolution evidence informed classification; repeated email quotations were not counted as separate reports.

| Category | Issues | FAQ entries | Report |
| --- | ---: | ---: | --- |
| Correctness and numerical bugs | 83 | 1 | [Correctness and numerical bugs](bugs.md) |
| Installation and dependencies | 60 | 0 | [Installation and dependencies](installation.md) |
| Inputs, workflows, and diagnostics | 83 | 6 | [Inputs, workflows, and diagnostics](workflow.md) |
| Features and compatibility | 45 | 1 | [Features and compatibility](features.md) |
| Resources and data availability | 54 | 0 | [Resources and data availability](resources.md) |
| Methodology and interpretation | 118 | 19 | [Methodology and interpretation](methodology.md) |
| Administrative and insufficient-information records | 15 | 1 | [Administrative and insufficient-information records](administrative.md) |

## Highest-priority current concerns

| Priority | Concern and evidence | Smallest useful next step |
| --- | --- | --- |
| High | Direct legacy summary-statistics ingestion retains infinite `Z`/`N` and nonpositive `N`; a synthetic workflow probe reached fitting, with invalid `N` sometimes producing an estimate ([V1](verification.md#v1-numeric-invariants-on-direct-legacy-input)). | Enforce finite numeric fields and strictly positive sample sizes at the shared ingestion/pre-fit boundary; identify offending rows and abort rather than impute. |
| High | Positive whole-data SNP heritabilities can coexist with invalid genetic-correlation delete-block denominators; the current fit/summary path produced `rg_se=NaN`, `p=NaN`, yet `status='ok'` ([V2](verification.md#v2-genetic-correlation-delete-block-boundary)). | Diagnose invalid trait/block denominators and emit an explicit not-estimable result; preserve valid component estimates without claiming valid correlation uncertainty. |
| High | A1 orientation, conditional coefficient tests, marginal enrichment and standardized coefficients are repeatedly confused, sometimes in suggested solutions. | Add short objective-to-quantity examples and visible result interpretation; emphasize that A1 being the signed-statistic allele does not imply A1 always increases the trait. |
| Moderate | A valid zero-total boundary aborts while computing a derived proportion ([V3](verification.md#v3-zero-total-derived-proportion)). | Preserve estimable total/coefficient results and mark undefined ratios unavailable with a reason; never add an arbitrary denominator epsilon. |
| Moderate | BED/gene projection is binary membership; consuming numeric annotations does not provide interval-overlap-count generation ([V4](verification.md#v4-binary-projection-and-header-delimiters)). | Consider one explicit count-projection option using the existing interval sweep, after documenting the changed annotation meaning. |
| Moderate | Historic resource links and package-name ambiguity can send users to missing data or a different Python distribution. | Provide a small, versioned resource-role table and repository-specific installation recipes, with compatibility requirements beside links. |

No critical finding was established by the bounded checks.
Severity describes potential harm; confidence describes the evidence supporting applicability.
A high-severity conceptual warning is not automatically a confirmed code defect.
The [verification record](verification.md) states exactly which paths ran.

## Frequency by requested workflow

A count is the number of distinct issues whose substantive question or proposed solution concerns that workflow.
Each issue is counted at most once within a group and may belong to several groups.
References to a command in a pasted log, without a substantive concern about that workflow, do not by themselves justify another tag.
The denominator for issue percentages is 458; these percentages are not population incidence, independent bug counts, or estimates of affected users.
Continuous annotations are a cross-cutting group.
Baseline-only partitioned heritability means a joint model of supplied baseline annotations; LDSC-SEG/query testing means an additional query fitted with its controls and baseline.
A partitioning question whose model cannot be determined remains unassigned to those two specific regimes rather than being guessed into one.

| Workflow | All issues | Open | Closed | % of 458 | Membership |
| --- | ---: | ---: | ---: | ---: | --- |
| munge-sumstats | 129 | 92 | 37 | 28.2% | [Issue links](workflow-index.md#munge-sumstats) |
| LD-score calculation | 85 | 60 | 25 | 18.6% | [Issue links](workflow-index.md#ldscore) |
| SNP heritability | 105 | 77 | 28 | 22.9% | [Issue links](workflow-index.md#h2) |
| Genetic correlation | 134 | 95 | 39 | 29.3% | [Issue links](workflow-index.md#rg) |
| Baseline-only partitioned heritability | 65 | 49 | 16 | 14.2% | [Issue links](workflow-index.md#baseline-partitioned-h2) |
| LDSC-SEG / additional query-annotation testing | 68 | 54 | 14 | 14.8% | [Issue links](workflow-index.md#ldsc-seg) |
| Continuous annotations | 11 | 8 | 3 | 2.4% | [Issue links](workflow-index.md#continuous) |

[Workflow membership lists](workflow-index.md) provide every contributing issue link and identify cross-cutting/unassigned issues.
FAQ entries are excluded from the 458-issue percentages.
Their separate counts are:

| Workflow | FAQ entries |
| --- | ---: |
| munge-sumstats | 0 |
| LD-score calculation | 7 |
| SNP heritability | 13 |
| Genetic correlation | 12 |
| Baseline-only partitioned heritability | 9 |
| LDSC-SEG / additional query-annotation testing | 7 |
| Continuous annotations | 4 |

## Reading the assessments

Each category starts with a severity-ordered summary table and then an individual entry for each assigned issue or FAQ item.
Within a severity level, current exposure and recurring mechanisms precede resolved historical concerns.
The source's open/closed state is retained separately from the audit's legacy-resolution assessment.
A closed issue, a referenced commit, or an unmerged PR is not treated as a confirmed fix.

LDSC3 applicability uses these labels: confirmed persistent behavior; apparently prevented by current code; partial coverage; intentionally unsupported; methodological limitation; resource-only concern; or indeterminate.
Source inspection supports qualified conclusions; only a matching executed probe supports a reproduced current failure.
A successful focused test is evidence for its scenario, not proof that every historical input is repaired.
Where an original file, platform or environment is unavailable, the entry identifies the limit rather than inventing a diagnosis.

## Practical prevention priorities

1. Reject invalid numeric inputs before fitting and make undefined estimates visibly not estimable.
2. Explain signed alleles, sample-size conventions, SNP identity/build, and reference versus regression SNPs at the point where users supply them.
3. Select coefficients, enrichment or standardized coefficients by the estimand and null hypothesis, with examples; do not imply that conditional association identifies a causal cell type or pathway.
4. Provide working, repository-specific installation choices and compatible resource combinations; avoid a new resource manager or installation framework.
5. Preserve existing alignment, window, collinearity, query-failure and output-publication safeguards, and add future regression cases at their existing seams.

These are audit recommendations only.
No production code, tests, resources, installation guide, GitHub issue or PR was modified.
The audit used local focused tests and synthetic boundary probes, without HPC jobs, large reference downloads, full historical environment recreation, or a genome-wide numerical-equivalence study.
