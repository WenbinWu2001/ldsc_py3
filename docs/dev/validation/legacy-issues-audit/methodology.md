# Methodology and interpretation

Last updated on: 2026-10-28

[Audit overview](README.md) · [Workflow counts and membership](workflow-index.md) · [Verification](verification.md)

## How to choose the quantity and test

The table concerns the scientific question, not which variable happens to be a table row or column.
The same annotation can require different summaries for different questions.
The implementation references below establish current output semantics; the original papers establish the statistical model and its limits.

| Scientific question | Quantity to read | Test or interpretation | Main interpretation limit |
| --- | --- | --- | --- |
| What SNP heritability does this GWAS support? | `total_h2_obs` and its SE; `total_h2_liab` only with appropriate prevalences | No dedicated headline P-value; report the signed estimate and its SE | The specified population, SNP denominator, summary-statistic calibration and model matter. |
| Does a binary annotation concentrate heritability relative to its size? | `enrichment`, `prop_h2`, `prop_snps` | `enrichment_p`: two-sided inside-versus-complement per-SNP contrast, implemented with a t reference and `n_blocks` degrees of freedom | Overlap-aware shares include contributions through other annotations and need not sum to one across sets. |
| Does a query add a positive contribution conditional on the chosen baseline and controls? | `coefficient`, `coefficient_se`, `coefficient_p` | One-sided normal test of a positive conditional coefficient | Association is conditional on this model, not proof of causal mediation or unique cell-type relevance. |
| How does a continuous annotation relate to per-SNP heritability? | Raw coefficient per annotation unit, or `tau_star` with its stated normalization | The standardized-coefficient table reports two-sided zero tests | Binary-style annotation shares or enrichment are not ordinary set enrichments for arbitrary numeric values. |
| How much fitted heritability lies in low versus high values of a continuous score? | `quantile-h2` totals, shares and quantile enrichment | The quantile contrast is separate from the coefficient test | This projects one saved joint model; it does not fit independent models for each bin or prove a causal dose-response. |
| How strongly do two traits share genome-wide genetic effects? | Genetic correlation, both component heritabilities, covariance and their uncertainty | The reported correlation zero test | Low or undefined component estimates can make the ratio unstable; genetic correlation is not a causal effect. |

Sources: [partitioned result contracts](../../../reference/results/partitioned-h2.md), [quantile result contracts](../../../reference/results/quantile-h2.md), `overlap_aware_category_table()` in [`overlap_model.py`](../../../../src/ldsc/_regression/overlap_model.py), and `RG` in [`estimators.py`](../../../../src/ldsc/_regression/estimators.py).
The published bases are [Finucane et al. 2015](https://doi.org/10.1038/ng.3404), [Finucane et al. 2018](https://doi.org/10.1038/s41588-018-0081-4), [Gazal et al. 2017](https://doi.org/10.1038/ng.3954), and [Bulik-Sullivan et al. 2015, genetic correlations](https://doi.org/10.1038/ng.3406).

### Conditional signal and marginal concentration can disagree

Consider a binary baseline annotation covering 80% of SNPs and a disjoint query covering 10%.
Suppose the fitted per-SNP heritability is proportional to $1+9B+Q$.
The query coefficient is positive, but its SNPs have mean 2 compared with a genome-wide mean of $1+9(0.8)+0.1=8.3$.
Its marginal enrichment is therefore $2/8.3<1$.
A positive conditional coefficient without marginal enrichment is entirely possible; it is not intrinsically contradictory or necessarily rare.
Whether either test is significant also depends on uncertainty.
This algebra illustrates the model rather than reproducing a real GWAS.

The current conceptual draft at [`partitioned-h2-math.md`](../../../user-guide/concepts/partitioned-h2-math.md), especially its four-case interpretation near lines 247–250, needs qualification.
Failure to reject a coefficient or enrichment null does not establish absence of SNP heritability.
A significant positive coefficient does not establish a specific pathological mechanism.
The phrase “de-confounded signal” in the result guide also needs to specify conditioning on measured annotations under the fitted model.
These are proposed documentation corrections; the existing pages were not edited.

### Alleles, extreme probabilities and sample sizes

A1 is the allele to which the signed statistic refers.
The `--a1-inc` assertion is stronger: A1 must be the increasing allele for every retained SNP.
Using that flag simply because A1 is the reported effect allele can destroy direction information and invalidate genetic correlation.
This distinction is already present in command help, but the issue discussions show a discoverability problem.
See [`munge_kernel.py`](../../../../src/ldsc/_sumstats/munge_kernel.py), the `a1_inc` branch, and issue #405 below.

Raw munging reconstructs the magnitude of Z from P and uses the signed statistic for direction.
Providing an additional signed Z column does not bypass an already-underflowed P value in that path.
Direct compatible munged summary statistics use a different entry path.
Sample size must describe the data underlying each statistic; combined discovery/replication N, case-control totals and effective N are not interchangeable conventions.
The audit recommends explicit examples and early invariant checks, not a universal automatic rescaling.

### Statistical limitations and documented future work

- Reference LD, population and marker ascertainment must support the analysis.
  A valid file format or a resource labelled with an ancestry does not establish that ordinary LDSC is calibrated for an admixed, targeted-array, selected-hit or nonhuman analysis.
- Low signal and near-zero denominators limit ratio inference even when the software has no algebraic defect.
  A software improvement should expose non-estimability, rather than conceal it with clipping, absolute values, extra constraints or arbitrary regularization.
- Correlated annotations limit what can be separated, including after deleting a jackknife block.
  A large baseline can reduce misspecification bias but does not guarantee identification or remove all biological confounding.
- The legacy FAQ explicitly identifies targeted-array handling as future work.
  The original single-trait paper identifies admixed populations and more robust all-variant heritability estimation as research directions; these are historical statements, not evidence that every such problem remains unsolved today ([Bulik-Sullivan et al. 2015, Discussion and assumptions](https://pmc.ncbi.nlm.nih.gov/articles/PMC4495769/)).
- Cell-type annotations can be correlated and expression profiles can reflect cell composition.
  The LDSC-SEG paper discusses resolution, tissue coverage and expression-data limitations; a significant conditional annotation is a useful prioritization result, not a mechanistic experiment ([Finucane et al. 2018, Discussion](https://pmc.ncbi.nlm.nih.gov/articles/PMC5896795/)).

For cross-trait fits, the idealized overlap contribution is $\rho N_s/\sqrt{N_1N_2}$, where $\rho$ is the phenotypic correlation among overlapping samples.
The original paper's equation (1) and discussion require more than overlap headcount to justify fixing that intercept; shared confounding must also be addressed.
This was checked against the local paper's PDF page 2 and the [published paper](https://doi.org/10.1038/ng.3406).

Gazal et al.'s equation (5), local PDF page 8, scales the raw coefficient by annotation SD and the mean per-SNP heritability over the SNPs used for that heritability total.
LDSC3's current quantile workflow explicitly uses its full common reference-SNP universe, including target-excluded SNPs for standardization.
Computing SD only on filtered regression SNPs can therefore change the normalization and does not reproduce the current contract.
See `compute_standardized_coefficients()` and `QuantileInputStatistics` in [`quantile_h2.py`](../../../../src/ldsc/_regression/quantile_h2.py) and [`quantile_inputs.py`](../../../../src/ldsc/_regression/quantile_inputs.py).

## Recommendations by kind

- **Code:** numeric input invariants and explicit unestimable states are detailed in the bugs report; a refactor alone cannot repair violated sampling/model assumptions.
- **Documentation:** add the objective-to-quantity table, the signed-allele distinction, the conditional-versus-marginal example, and assumptions beside every recommended interpretation.
- **Features:** favor small diagnostics exposing the retained population and invalid uncertainty over new estimators without independent validation.
- **Resources:** retain ancestry, genome build, marker universe and baseline/control compatibility with every suggested resource.

## Category summary

118 issues and 19 FAQ entries.

| Severity | Source | Concern | LDSC3 disposition | Evidence level |
| --- | --- | --- | --- | --- |
| High | [#98](#issue-98) | One-sided P supplied where two-sided P is expected | Partial coverage | Source inspection |
| High | [#136](#issue-136) | Control gene list means the selection universe | Partial coverage | Source inspection |
| High | [#138](#issue-138) | Single-annotation model assigns it all modeled heritability | Partial coverage | Source inspection |
| High | [#405](#issue-405) | Signed statistic and a1-inc misuse | Partial coverage | Source inspection |
| High | [#425](#issue-425) | Alleles, restriction and overlap | Partial coverage | Source inspection |
| High | [#438](#issue-438) | Continuous enrichment misuse | Partial coverage | Source inspection |
| High | [FAQ-11](#faq-11) | Nonpositive heritability in correlations | Partial coverage | Source inspection |
| High | [#31](#issue-31) | High ratio with weak polygenic signal | Methodological limitation | Source inspection |
| High | [#38](#issue-38) | Very small GWAS yields negative heritability | Methodological limitation | Source inspection |
| High | [#48](#issue-48) | Selected-array meta-analysis and genome-wide extrapolation | Methodological limitation | Source inspection |
| High | [#52](#issue-52) | Cross-ancestry correlation is a different model | Methodological limitation | Source inspection |
| High | [#53](#issue-53) | Effective sample size for related-person GWAS | Methodological limitation | Source inspection |
| High | [#70](#issue-70) | Undoing genomic control | Methodological limitation | Source inspection |
| High | [#71](#issue-71) | Sample overlap does not require a fixed intercept | Methodological limitation | Source inspection |
| High | [#78](#issue-78) | Out-of-range correlation and inconsistent displays | Methodological limitation | Source inspection |
| High | [#82](#issue-82) | GC correction and fixed-intercept P values conflated | Methodological limitation | Source inspection |
| High | [#87](#issue-87) | Renamed covariance-intercept flag mistaken for a recommendation | Methodological limitation | Source inspection |
| High | [#88](#issue-88) | Negative h2 cannot be bypassed with a population definition | Methodological limitation | Source inspection |
| High | [#89](#issue-89) | Large SEs and correlations outside [-1,1] | Methodological limitation | Source inspection |
| High | [#92](#issue-92) | P values from conditional coefficients versus enrichment | Methodological limitation | Source inspection |
| High | [#94](#issue-94) | Unavailable correlation with negative low-signal h2 | Methodological limitation | Source inspection |
| High | [#95](#issue-95) | Case-control total N and effective N conventions | Methodological limitation | Source inspection |
| High | [#102](#issue-102) | Padded annotations and conditional coefficients | Methodological limitation | Source inspection |
| High | [#108](#issue-108) | Wrong test reproduced the wrong cell-type figure | Methodological limitation | Source inspection |
| High | [#123](#issue-123) | Missing chromosomes and invalid genome-wide extrapolation | Methodological limitation | Source inspection |
| High | [#129](#issue-129) | Disease-state gene lists, model differences and resource releases | Methodological limitation | Source inspection |
| High | [#151](#issue-151) | More imputed SNPs need not raise LDSC heritability | Methodological limitation | Source inspection |
| High | [#152](#issue-152) | Shared controls and binary covariance intercept | Methodological limitation | Source inspection |
| High | [#155](#issue-155) | Negative slope and ratio above one | Methodological limitation | Source inspection |
| High | [#156](#issue-156) | Annotation heritability share outside [0,1] | Methodological limitation | Source inspection |
| High | [#189](#issue-189) | Significantly out-of-range genetic correlation | Methodological limitation | Source inspection |
| High | [#191](#issue-191) | Sample size under extreme imbalance | Methodological limitation | Source inspection |
| High | [#195](#issue-195) | Mixed-model GWAS in animal populations | Methodological limitation | Source inspection |
| High | [#200](#issue-200) | ATAC-seq query/control design | Methodological limitation | Source inspection |
| High | [#203](#issue-203) | GC-corrected statistics and low ratio | Methodological limitation | Source inspection |
| High | [#204](#issue-204) | Unknown continuous values versus zero | Methodological limitation | Source inspection |
| High | [#210](#issue-210) | Testing annotation intersections in a genome subset | Methodological limitation | Source inspection |
| High | [#213](#issue-213) | Increasing blocks to avoid singular fits | Methodological limitation | Source inspection |
| High | [#218](#issue-218) | Overlapping samples for genetic correlation | Methodological limitation | Source inspection |
| High | [#222](#issue-222) | Ratio interpretation with weak polygenic signal | Methodological limitation | Source inspection |
| High | [#229](#issue-229) | Out-of-range SNP heritability at N2000 | Methodological limitation | Source inspection |
| High | [#249](#issue-249) | Genetic versus phenotypic correlation | Methodological limitation | Source inspection |
| High | [#251](#issue-251) | METAL Weight as sample size | Methodological limitation | Source inspection |
| High | [#252](#issue-252) | Unbounded estimates with four hundred samples | Methodological limitation | Source inspection |
| High | [#255](#issue-255) | Small imbalanced cohorts and nested meta-analysis | Methodological limitation | Source inspection |
| High | [#257](#issue-257) | Genetic architecture and standardized chi-square derivation | Methodological limitation | Source inspection |
| High | [#259](#issue-259) | Large correlation can be nonsignificant | Methodological limitation | Source inspection |
| High | [#265](#issue-265) | High ratio with weak polygenic component | Methodological limitation | Source inspection |
| High | [#271](#issue-271) | Meaning and danger of --a1-inc | Methodological limitation | Source inspection |
| High | [#283](#issue-283) | Score-test output as signed direction | Methodological limitation | Source inspection |
| High | [#288](#issue-288) | One-sided tests do not constrain estimates | Methodological limitation | Source inspection |
| High | [#299](#issue-299) | QTL summary statistics as GWAS inputs | Methodological limitation | Source inspection |
| High | [#305](#issue-305) | Very small annotations and calibration | Methodological limitation | Source inspection |
| High | [#308](#issue-308) | Constrain intercept to obtain desired heritability | Methodological limitation | Source inspection |
| High | [#310](#issue-310) | Very small first study despite large partner | Methodological limitation | Source inspection |
| High | [#313](#issue-313) | Wide intervals versus UK Biobank | Methodological limitation | Source inspection |
| High | [#318](#issue-318) | Continuous AFR-LD annotation misinterpreted as cross-ancestry transfer | Methodological limitation | Source inspection |
| High | [#328](#issue-328) | Large genetic-correlation estimate with huge SE | Methodological limitation | Source inspection |
| High | [#332](#issue-332) | Sign of binary/continuous genetic correlation | Methodological limitation | Source inspection |
| High | [#334](#issue-334) | BOLT genotype versus imputed output | Methodological limitation | Source inspection |
| High | [#339](#issue-339) | Rare-variant SNP heritability and allele filtering | Methodological limitation | Source inspection |
| High | [#362](#issue-362) | Intercept constraint changes target | Methodological limitation | Source inspection |
| High | [#363](#issue-363) | Overlap intercept identification | Methodological limitation | Source inspection |
| High | [#415](#issue-415) | Repeated-observation sample size | Methodological limitation | Source inspection |
| High | [#419](#issue-419) | Shared-cohort intercept | Methodological limitation | Source inspection |
| High | [#431](#issue-431) | Constraint selected by significance | Methodological limitation | Source inspection |
| High | [#436](#issue-436) | Which p-value across traits | Methodological limitation | Source inspection |
| High | [#444](#issue-444) | Low heritability and covariance intercept | Methodological limitation | Source inspection |
| High | [#452](#issue-452) | Model dependence and subset contrasts | Methodological limitation | Source inspection |
| High | [#458](#issue-458) | Near-boundary rg with weak traits | Methodological limitation | Source inspection |
| High | [#467](#issue-467) | Targeted-array correlation | Methodological limitation | Source inspection |
| High | [#472](#issue-472) | Negative heritability denominator | Methodological limitation | Source inspection |
| High | [#476](#issue-476) | Baseline choice follows analysis aim | Methodological limitation | Source inspection |
| High | [#487](#issue-487) | Contributors, output rows and panel samples | Methodological limitation | Source inspection |
| High | [#488](#issue-488) | Small LD-score row count | Methodological limitation | Source inspection |
| High | [FAQ-10](#faq-10) | Power and sample size | Methodological limitation | Source inspection |
| High | [FAQ-12](#faq-12) | Unbounded correlation estimates | Methodological limitation | Source inspection |
| High | [FAQ-13](#faq-13) | Imprecise genetic correlations | Methodological limitation | Source inspection |
| High | [FAQ-15](#faq-15) | Negative total estimates | Methodological limitation | Source inspection |
| High | [FAQ-17](#faq-17) | Estimates above parameter bounds | Methodological limitation | Source inspection |
| High | [FAQ-18](#faq-18) | Negative annotation contributions | Methodological limitation | Source inspection |
| High | [FAQ-20](#faq-20) | Intercept estimates below unity | Methodological limitation | Source inspection |
| High | [FAQ-21](#faq-21) | Choosing intercept constraints | Methodological limitation | Source inspection |
| High | [FAQ-22](#faq-22) | Small heritability estimates | Methodological limitation | Source inspection |
| High | [FAQ-25](#faq-25) | Targeted genotyping designs | Methodological limitation | Source inspection |
| High | [#388](#issue-388) | Partitioned correlation misconception | Intentionally unsupported | Source inspection |
| High | [#441](#issue-441) | Single-variant PVE versus LDSC expectation | Intentionally unsupported | Source inspection |
| High | [#451](#issue-451) | Per-SNP heritability is not latent output | Intentionally unsupported | Source inspection |
| High | [#498](#issue-498) | Survival GWAS scale | Intentionally unsupported | Source inspection |
| High | [#109](#issue-109) | Base proportion above one and conflicting P values | Indeterminate | Source inspection |
| High | [#430](#issue-430) | Version and uncertainty mismatch | Indeterminate | Source inspection |
| High | [#468](#issue-468) | Correct univariate fits do not validate signs | Indeterminate | Source inspection |
| High | [#100](#issue-100) | Coefficient P value must specify sidedness | Apparently prevented | Source inspection |
| High | [#225](#issue-225) | Single-annotation proportion is always one | Apparently prevented | Source inspection |
| High | [#375](#issue-375) | Effect versus reference allele | Apparently prevented | Source inspection |
| High | [#381](#issue-381) | Observed/liability comparison | Apparently prevented | Source inspection |
| High | [#393](#issue-393) | Rows versus annotated membership | Apparently prevented | Source inspection |
| High | [#439](#issue-439) | Contradictory allele advice | Apparently prevented | Source inspection |
| High | [#447](#issue-447) | Query and all-gene control roles | Apparently prevented | Source inspection |
| High | [#484](#issue-484) | ATAC query conditioning and enrichment | Apparently prevented | Source inspection |
| High | [FAQ-16](#faq-16) | Observed and liability scales | Apparently prevented | Source inspection |
| High | [FAQ-19](#faq-19) | Overlapping annotation interpretation | Apparently prevented | Source inspection |
| Moderate | [#3](#issue-3) | Case-control heritability in weight updates | Methodological limitation | Source inspection |
| Moderate | [#50](#issue-50) | Missing significant SNPs from regression set | Methodological limitation | Source inspection |
| Moderate | [#56](#issue-56) | Correlation, SE, Z and P columns conflated | Methodological limitation | Source inspection |
| Moderate | [#67](#issue-67) | Total-h2 Z score and base proportion | Methodological limitation | Source inspection |
| Moderate | [#74](#issue-74) | Pairwise h2 differs from standalone h2 | Methodological limitation | Source inspection |
| Moderate | [#130](#issue-130) | Self-correlation near one and covariance intercept below one | Methodological limitation | Source inspection |
| Moderate | [#141](#issue-141) | Unequal and small gene-list sizes | Methodological limitation | Source inspection |
| Moderate | [#153](#issue-153) | Why chi-square is approximated by N beta squared | Methodological limitation | Source inspection |
| Moderate | [#154](#issue-154) | Minimum annotation size and power | Methodological limitation | Source inspection |
| Moderate | [#193](#issue-193) | Meaning of genetic-correlation p-value | Methodological limitation | Source inspection |
| Moderate | [#236](#issue-236) | Signed-statistic median sanity check | Methodological limitation | Source inspection |
| Moderate | [#262](#issue-262) | P column absent from munged output | Methodological limitation | Source inspection |
| Moderate | [#284](#issue-284) | Variant count enrichment versus heritability enrichment | Methodological limitation | Source inspection |
| Moderate | [#410](#issue-410) | Unbounded estimate with large uncertainty | Methodological limitation | Source inspection |
| Moderate | [#424](#issue-424) | GWAS frequency versus reference MAF | Methodological limitation | Source inspection |
| Moderate | [#427](#issue-427) | Negative ratio interpretation | Methodological limitation | Source inspection |
| Moderate | [#437](#issue-437) | UKB-specific defaults | Methodological limitation | Source inspection |
| Moderate | [#474](#issue-474) | Negative unbiased LD scores | Methodological limitation | Source inspection |
| Moderate | [#475](#issue-475) | Negative estimate near zero | Methodological limitation | Source inspection |
| Moderate | [#481](#issue-481) | SE and scientific validation | Methodological limitation | Source inspection |
| Moderate | [#494](#issue-494) | Near-zero MAF/LD correlation | Methodological limitation | Source inspection |
| Moderate | [#501](#issue-501) | Which stage needs many SNPs | Methodological limitation | Source inspection |
| Moderate | [FAQ-04](#faq-04) | Negative unbiased LD scores | Methodological limitation | Source inspection |
| Moderate | [FAQ-08](#faq-08) | Numeric annotation sums | Methodological limitation | Source inspection |
| Moderate | [FAQ-23](#faq-23) | Imprecise heritability estimates | Methodological limitation | Source inspection |
| Moderate | [FAQ-26](#faq-26) | Self-correlation numerical identity | Methodological limitation | Source inspection |
| Moderate | [FAQ-27](#faq-27) | Same label across different cohorts | Methodological limitation | Source inspection |
| Moderate | [#390](#issue-390) | Testing correlation against one | Intentionally unsupported | Source inspection |
| Moderate | [#256](#issue-256) | rsIDs versus genome-build coordinates | Apparently prevented | Source inspection |
| Moderate | [#270](#issue-270) | Raw versus standardized coefficient | Apparently prevented | Source inspection |
| Moderate | [#285](#issue-285) | Pair-specific SNP-heritability estimates differ | Apparently prevented | Source inspection |
| Moderate | [#336](#issue-336) | Standardized-coefficient z-score | Apparently prevented | Source inspection |
| Low | [#280](#issue-280) | Small deviations from tutorial numbers | Indeterminate | Source + image inspection |
| Low | [#215](#issue-215) | Parentheses in log summaries | Apparently prevented | Source inspection |
| Low | [FAQ-28](#faq-28) | Efficient block jackknife | Apparently prevented | Source inspection |

## Individual assessments

<a id="issue-98"></a>

### #98: Munge script modifies the Z score values, may be a bug??

[Original issue](https://github.com/bulik/ldsc/issues/98) · Source state: **closed** · Severity: **high** · Accessible comments reviewed: **6**.

Workflow tags: `munge-sumstats`.

**Trigger and harm.** User computes P=normal survival at abs(Z) without multiplying by two.
Reconstructing Z magnitude from P inflates statistics and lambda.

**Legacy disposition.** User identifies and fixes the missing factor of two; maintainer had reproduced correct behavior for a separate dataset.

**LDSC3 assessment — Partial coverage.** LDSC3 intentionally derives Z magnitude from P using a one-degree-of-freedom chi-square and takes only the sign from the supplied statistic.
Therefore inconsistent P and signed Z can still silently change magnitudes.
This is not evidence p_to_z is mathematically wrong, but the interface needs explicit two-sided-test expectations and a useful consistency diagnostic.

**Recommended prevention.** Document P as two-sided one-df association P and consider a sampled P-versus-Z/beta-SE consistency check with declared tolerances; do not preserve conflicting inputs silently.

**Evidence and verification.** source-inspected
[src/ldsc/_sumstats/munge_kernel.py:p_to_z](<../../../../src/ldsc/_sumstats/munge_kernel.py>); [src/ldsc/_sumstats/munge_kernel.py:munge_sumstats](<../../../../src/ldsc/_sumstats/munge_kernel.py>); [docs/reference/munge-sumstats.md](<../../../reference/munge-sumstats.md>).

<a id="issue-136"></a>

### #136: [Cell-type specific analysis] How are the 'control' gene sets determined?

[Original issue](https://github.com/bulik/ldsc/issues/136) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **3**.

Workflow tags: `ldsc-seg`.

**Trigger and harm.** Users confuse control genes with healthy controls or omit them when testing specifically expressed genes.
Enrichment near genes generally can be misattributed to the selected tissue-specific set.

**Legacy disposition.** Maintainer defines controls as all genes considered/passing QC in the expression-selection analysis, including Cahoy source tables.

**LDSC3 assessment — Partial coverage.** LDSC3 supports --control-gene-list-file and adds its annotation to baseline, but the choice is scientifically user-owned and optional input does not verify appropriateness. A generic all-genes list may be wrong if the tested list came from a restricted measurable/QC universe.
Passing one tissue query and obtaining output does not validate the experiment.

**Recommended prevention.** Describe control-list provenance beside gene-list generation and use the exact selection universe; do not hard-code one universal gene background.

**Evidence and verification.** source-inspected
[src/ldsc/_annotation/gene_lists.py](<../../../../src/ldsc/_annotation/gene_lists.py>); [docs/user-guide/workflows/gene-list-ld-scores.md](<../../../user-guide/workflows/gene-list-ld-scores.md>); [Source](<https://github.com/bulik/ldsc/issues/136#issuecomment-452324374>).

<a id="issue-138"></a>

### #138: Prop._h2 always is 1.0 use my own ldscores based on my own annot

[Original issue](https://github.com/bulik/ldsc/issues/138) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **12**.

Workflow tags: `baseline-partitioned-h2`, `ldsc-seg`.

**Trigger and harm.** A custom binary annotation is fit alone and prop_h2 is always one.
Users interpret a model denominator as biological certainty or apply an incorrect overlap formula patch.

**Legacy disposition.** Maintainer explicitly explains the need for coverage outside the annotation and updates documentation; user suggestion replacing M_annot by M_tot is not endorsed.

**LDSC3 assessment — Partial coverage.** LDSC3 requires baseline annotations for query generation and warns about a degenerate single-annotation functional model.
However a supplied baseline is not necessarily an all-ones base, so the coverage assumption still needs validation/documentation.
Adding duplicate all-ones columns to both baseline and query causes a different singularity.

**Recommended prevention.** Show base plus query once per joint model; explain modeled versus genome-wide proportions and explicitly reject the unsupported denominator patch.

**Evidence and verification.** source-inspected
[src/ldsc/_regression/runner.py:_fit_functional_model](<../../../../src/ldsc/_regression/runner.py>); [docs/reference/results/partitioned-h2.md](<../../../reference/results/partitioned-h2.md>); [src/ldsc/_regression/overlap_model.py](<../../../../src/ldsc/_regression/overlap_model.py>).

<a id="issue-405"></a>

### #405: 'Could not find a signed summary statistic column.')

[Original issue](https://github.com/bulik/ldsc/issues/405) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **3**.

Workflow tags: `munge-sumstats`, `rg`.

**Trigger and harm.** A nonstandard beta1 header is not inferred; a user later enables --a1-inc merely because effects refer to A1.
That mistaken option can discard negative signs and invalidate genetic correlation while appearing successful.

**Legacy disposition.** A commenter identifies the column-hint problem and the user reports success, but their explanation of --a1-inc is scientifically incorrect; no code fix is shown.

**LDSC3 assessment — Partial coverage.** LDSC3 supports explicit signed-statistic hints, but a1_inc still skips median checking and sign restoration.
Correct A1 coding alone does not imply every A1 increases the trait.

**Recommended prevention.** Prominently reserve --a1-inc for already trait-increasing A1 coding, show signed beta/Z examples, and diagnose contradictory negative signed inputs if that option is requested.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[src/ldsc/_sumstats/munge_input.py:135](<../../../../src/ldsc/_sumstats/munge_input.py#L135>); [src/ldsc/_sumstats/munge_kernel.py:732](<../../../../src/ldsc/_sumstats/munge_kernel.py#L732>); [src/ldsc/_sumstats/munge_kernel.py:744](<../../../../src/ldsc/_sumstats/munge_kernel.py#L744>).

<a id="issue-425"></a>

### #425: Is it necessary to use the --merge-alleles parameter when doing genetic correlation in ldsc?

[Original issue](https://github.com/bulik/ldsc/issues/425) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **16**.

Workflow tags: `munge-sumstats`, `rg`, `ldscore`.

**Trigger and harm.** Equal row counts or intersection by rsID do not prevent incompatible alleles; later discussion concerns HM3 ancestry and covariance intercepts.
Indels/discordant alleles can fail alignment, and inaccurate community explanations can cause mistaken overlap assumptions.

**Legacy disposition.** Author reports success after removing indels. A later comment incorrectly says LDSC assumes no overlap, then clarifies that the intercept can be estimated.
No code fix is established.

**LDSC3 assessment — Partial coverage.** LDSC3 filters unsupported SNP alleles, matches identity explicitly and aligns effect signs during regression.
HM3 restriction, allele harmonization and ancestry-matched LD are distinct operations; an intercept is not an identified headcount of shared people.

**Recommended prevention.** Teach the three separate steps and give mismatch examples; state that a free covariance intercept accommodates overlap under assumptions but does not by itself identify sample overlap.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[src/ldsc/_sumstats/munge_kernel.py:232](<../../../../src/ldsc/_sumstats/munge_kernel.py#L232>); [src/ldsc/_regression/estimators.py:786](<../../../../src/ldsc/_regression/estimators.py#L786>); [docs/reference/ldscore-snp-universes.md](<../../../reference/ldscore-snp-universes.md>).

<a id="issue-438"></a>

### #438: Prop._h2 is negative

[Original issue](https://github.com/bulik/ldsc/issues/438) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **15**.

Workflow tags: `baseline-partitioned-h2`, `continuous`.

**Trigger and harm.** Weighted continuous annotations yield extreme/negative prop_h2 and enrichment; follow-up confuses raw and standardized coefficients.
Binary-style output columns can be interpreted as impossible biological proportions or used for the wrong test.

**Legacy disposition.** Comments correctly warn against binary interpretation but misdescribe raw coefficient units and recommend a questionable regression-only SD universe; no code fix.

**LDSC3 assessment — Partial coverage.** LDSC3 retains weighted legacy summaries with a log warning and offers quantile_h2 and standardized coefficients.
The documented common reference SNP reconstruction differs from the comment suggesting regression SNPs only. A passing run does not make binary-style summaries interpretable.

**Recommended prevention.** Make continuous-output interpretation visible in headline tables, define raw per-unit versus standardized effects and their reference SNP set, and direct users to coefficient tests or post-fit quantiles.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[docs/reference/results/partitioned-h2.md:28](<../../../reference/results/partitioned-h2.md#two-regimes>); [docs/reference/results/quantile-h2.md](<../../../reference/results/quantile-h2.md>); [src/ldsc/_regression/quantile_h2.py](<../../../../src/ldsc/_regression/quantile_h2.py>); [Gazal et al. 2017, PDF page 8, equation 5](<https://doi.org/10.1038/ng.3954>).

<a id="faq-11"></a>

### FAQ-11: Nonpositive heritability in correlations

[FAQ source](https://github.com/bulik/ldsc/wiki/FAQ) · Source-order identifier: **FAQ-11** · Severity: **high**.

Workflow tags: `h2`, `rg`.

**Legacy disposition.** Historical FAQ guidance; assessed independently rather than treated as a code-fix record.

**LDSC3 assessment — Partial coverage.** The historical sign explanation contains a reversal.
Current RG handles nonpositive whole-data SNP heritability, but positive whole-data values do not guarantee valid delete-block denominators; probe V2 reproduces that separate boundary.

**Recommended prevention.** Report correlation as not estimable with the affected traits/blocks, preserving interpretable covariance and heritability estimates where possible.

**Evidence and verification.** Source-inspected; relevant executed checks are identified in verification.md.
[src/ldsc/_regression/estimators.py:765](<../../../../src/ldsc/_regression/estimators.py#L765>); [src/ldsc/_regression/runner.py:2475](<../../../../src/ldsc/_regression/runner.py#L2475>); [verification.md#v2-genetic-correlation-delete-block-boundary](<verification.md#v2-genetic-correlation-delete-block-boundary>).
[Shared reproduced boundary V2](verification.md#v2-genetic-correlation-delete-block-boundary); the original historical dataset was not reproduced.

<a id="issue-31"></a>

### #31: Large Ratio values in univariate heritability estimations

[Original issue](https://github.com/bulik/ldsc/issues/31) · Source state: **closed** · Severity: **high** · Accessible comments reviewed: **6**.

Workflow tags: `h2`, `rg`, `baseline-partitioned-h2`.

**Trigger and harm.** Mean chi-square is near one and ratio estimates are compared with GCTA or tutorial values. A noisy ratio or differing SNP/model sets can be misread as proven confounding or a faulty estimator.

**Legacy disposition.** Maintainers identify reference-score differences, weak partitioning power and imprecise comparison estimates; no code defect is established.

**LDSC3 assessment — Methodological limitation.** LDSC3 computes the same ratio denominator and offers free/fixed intercepts.
Those controls cannot establish absence of confounding.
The thread shows that phenotype heterogeneity, reference choice and uncertainty matter; selecting a fixed intercept because the answer looks closer to another estimator is not a valid general remedy.

**Recommended prevention.** Explain ratio instability near mean chi-square one, compare confidence intervals on matched targets, and justify fixed intercepts independently of desired results.

**Evidence and verification.** source-inspected
[src/ldsc/_regression/estimators.py:Hsq](<../../../../src/ldsc/_regression/estimators.py>); [docs/user-guide/workflows/heritability.md](<../../../user-guide/workflows/heritability.md>); [docs/user-guide/concepts/partitioned-h2-math.md](<../../../user-guide/concepts/partitioned-h2-math.md>).

<a id="issue-38"></a>

### #38: h2 is negative -1.105

[Original issue](https://github.com/bulik/ldsc/issues/38) · Source state: **closed** · Severity: **high** · Accessible comments reviewed: **1**.

Workflow tags: `h2`, `rg`.

**Trigger and harm.** A study with roughly 300 samples estimates -1.105 with SE 1.57.
Uninformative sampling variation can be mistaken for negative biological heritability or repaired by clipping.

**Legacy disposition.** Maintainer attributes the result to insufficient polygenic signal; no software fix is appropriate.

**LDSC3 assessment — Methodological limitation.** LDSC3 deliberately retains unconstrained heritability estimates.
More explicit input validation cannot create information in low-signal summary statistics, and genetic correlation is especially unstable when the denominator is near zero. A universal sample-size threshold would be misleading because power depends on architecture, phenotype, quality and reference match.

**Recommended prevention.** Describe signal/uncertainty checks before interpretation and direct users with individual-level data toward appropriately validated alternatives; do not truncate estimates to zero.

**Evidence and verification.** source-inspected
[src/ldsc/_regression/estimators.py:Hsq](<../../../../src/ldsc/_regression/estimators.py>); [docs/user-guide/workflows/heritability.md](<../../../user-guide/workflows/heritability.md>); [docs/user-guide/workflows/genetic-correlation.md](<../../../user-guide/workflows/genetic-correlation.md>).

<a id="issue-48"></a>

### #48: Using Meta-analysis Results Containing MetaChip for LDSC

[Original issue](https://github.com/bulik/ldsc/issues/48) · Source state: **closed** · Severity: **high** · Accessible comments reviewed: **3**.

Workflow tags: `h2`, `rg`.

**Trigger and harm.** A GWAS backbone is combined with MetaboChip-enriched variants and larger N.
Genome-wide heritability can be biased by nonrepresentative marker selection and uneven sample-size architecture.

**Legacy disposition.** Maintainer gives explicitly nonrigorous rules of thumb: ignore univariate h2 on selected arrays and prefer the genome-wide stage where possible.

**LDSC3 assessment — Methodological limitation.** LDSC3 can read these summaries but cannot infer the ascertainment design from a header.
The historical claim that selected-array genetic correlation is valid is a heuristic requiring matched designs and assumptions, not a software guarantee. A simple SNP-count threshold cannot diagnose selection bias or rescue genome-wide extrapolation.

**Recommended prevention.** Document genome-wide representative input requirements and distinguish restricted-panel estimands; recommend stage-specific sensitivity comparisons with uncertainty.

**Evidence and verification.** source-inspected
[src/ldsc/_sumstats/munge_kernel.py:process_n](<../../../../src/ldsc/_sumstats/munge_kernel.py>); [docs/user-guide/workflows/heritability.md](<../../../user-guide/workflows/heritability.md>); [docs/user-guide/workflows/genetic-correlation.md](<../../../user-guide/workflows/genetic-correlation.md>); [Source](<https://github.com/bulik/ldsc/issues/48#issuecomment-198450937>).

<a id="issue-52"></a>

### #52: How can I perform genetic correlation of two traits from different population? European versus Asian.

[Original issue](https://github.com/bulik/ldsc/issues/52) · Source state: **closed** · Severity: **high** · Accessible comments reviewed: **3**.

Workflow tags: `rg`.

**Trigger and harm.** One GWAS is European and the other Asian.
Using one population LD model for both can misinterpret cross-population effect similarity.

**Legacy disposition.** Maintainers explicitly say ordinary LDSC is unsuitable and link the subsequently published Popcorn implementation as related work.

**LDSC3 assessment — Methodological limitation.** LDSC3 still fits the ordinary same-reference covariance/correlation model.
Accepting two files does not imply cross-ancestry validity, and matching genome build is insufficient to establish ancestry compatibility.
This is a fundamental model-scope limitation, not an absent parsing flag.
Popcorn is a separate method requiring its own assumptions and resources.

**Recommended prevention.** State same-ancestry/reference assumptions prominently and describe the linked cross-population method as a separate research route, with current primary-source verification before recommendation.

**Evidence and verification.** source-inspected
[src/ldsc/_regression/runner.py:build_rg_dataset](<../../../../src/ldsc/_regression/runner.py>); [src/ldsc/_regression/estimators.py:Gencov](<../../../../src/ldsc/_regression/estimators.py>); [Source](<https://github.com/brielin/Popcorn>).

<a id="issue-53"></a>

### #53: Specifying sample size when GWAS included related individuals

[Original issue](https://github.com/bulik/ldsc/issues/53) · Source state: **closed** · Severity: **high** · Accessible comments reviewed: **2**.

Workflow tags: `munge-sumstats`, `h2`.

**Trigger and harm.** A twin study uses GEE to account for within-family dependence.
Raw participant N may not match the variance scaling required by LDSC.

**Legacy disposition.** Maintainer calls effective N a hunch and explicitly says this design had not been tested.

**LDSC3 assessment — Methodological limitation.** LDSC3 accepts per-SNP or constant N but cannot derive a valid effective size for arbitrary GEE, twin or mixed-model designs.
The case/control effective-N expression is not automatically applicable to familial dependence.
More confident documentation than the original evidence would overstate validity.

**Recommended prevention.** Label family-based effective-N guidance uncertain; require study-specific variance calibration or simulation before treating this design as validated.

**Evidence and verification.** source-inspected
[src/ldsc/_sumstats/munge_kernel.py:process_n](<../../../../src/ldsc/_sumstats/munge_kernel.py>); [docs/reference/munge-sumstats.md](<../../../reference/munge-sumstats.md>); [Source](<https://github.com/bulik/ldsc/issues/53#issuecomment-226207023>).

<a id="issue-70"></a>

### #70: GC corrected data: Undo

[Original issue](https://github.com/bulik/ldsc/issues/70) · Source state: **closed** · Severity: **high** · Accessible comments reviewed: **2**.

Workflow tags: `h2`, `rg`.

**Trigger and harm.** Summary statistics have been scaled by genomic control once or within meta-analysis cohorts.
Naive rescaling can distort SNP heritability, covariance and intercept interpretation.

**Legacy disposition.** Maintainer gives a single-final-GC inverse and points to Proposition 3/formula 4.3 for cohort-specific GC; no built-in legacy correction is claimed.

**LDSC3 assessment — Methodological limitation.** LDSC3 does not infer correction history or accept cohort GC factors as a generic repair.
Scaling an already mixed meta-analysis by one lambda is generally not the stated derivation.
Approximate cancellation in correlation does not restore absolute heritability or justify fixing intercepts after correction.

**Recommended prevention.** Prefer uncorrected released summaries; document known single-scaling inversion separately from unidentifiable cohort histories and cite the original supplement.

**Evidence and verification.** source-inspected
[src/ldsc/_regression/estimators.py](<../../../../src/ldsc/_regression/estimators.py>); [docs/user-guide/workflows/heritability.md](<../../../user-guide/workflows/heritability.md>); [Source](<https://github.com/bulik/ldsc/issues/70#issuecomment-281153170>).

<a id="issue-71"></a>

### #71: Sample overlap

[Original issue](https://github.com/bulik/ldsc/issues/71) · Source state: **closed** · Severity: **high** · Accessible comments reviewed: **7**.

Workflow tags: `rg`.

**Trigger and harm.** Users interpret available intercept flags as recommended adjustments for overlapping participants.
Incorrect fixed values can bias estimates and produce spuriously small SEs.

**Legacy disposition.** Maintainer explains free intercept defaults, h2 null one, covariance null zero, and quantitative/binary phenotype-correlation formulas.

**LDSC3 assessment — Methodological limitation.** LDSC3 estimates intercepts by default and replaces no-intercept with explicit values.
Sample overlap concerns the cross-trait covariance intercept, not an automatic change to each univariate intercept.
Known overlap counts alone do not determine all confounding, and binary ascertainment affects expected covariance-intercept formulas.

**Recommended prevention.** Make free-intercept estimation the main tutorial and show fixed values only as justified sensitivity analyses, clearly separating the two intercepts.

**Evidence and verification.** source-inspected
[src/ldsc/_regression/runner.py:_fit_rg_dataset](<../../../../src/ldsc/_regression/runner.py>); [docs/reference/regression-options.md](<../../../reference/regression-options.md>); [docs/reference/legacy-flag-map.md](<../../../reference/legacy-flag-map.md>).

<a id="issue-78"></a>

### #78: Error: rg out of bounds

[Original issue](https://github.com/bulik/ldsc/issues/78) · Source state: **closed** · Severity: **high** · Accessible comments reviewed: **5**.

Workflow tags: `rg`.

**Trigger and harm.** Legacy log suppresses rg above 1.2 while its summary table prints 1.2452.
Users may interpret display suppression as a calculation failure or clip the result to one.

**Legacy disposition.** Maintainer explains unbounded sampling estimates and overlap-related instability; no methodological defect is demonstrated by an out-of-range estimate alone.

**LDSC3 assessment — Methodological limitation.** LDSC3 reports the unconstrained ratio in saved output with explicit statuses. A value significantly outside the feasible range suggests assumptions, overlap, effect orientation or reference mismatch require investigation; large SEs may simply indicate low information.
Two identical files yielding one is a useful identity check, not general validation of every pair.

**Recommended prevention.** Report estimate and uncertainty consistently and explain why clipping or suppression is not a statistical correction.

**Evidence and verification.** source-inspected
[src/ldsc/_regression/estimators.py:RG](<../../../../src/ldsc/_regression/estimators.py>); [src/ldsc/_regression/runner.py:_summarize_rg_pair](<../../../../src/ldsc/_regression/runner.py>); [docs/user-guide/workflows/genetic-correlation.md](<../../../user-guide/workflows/genetic-correlation.md>).

<a id="issue-82"></a>

### #82: about the p value of genetic correlation for GC correction GWAS

[Original issue](https://github.com/bulik/ldsc/issues/82) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **4**.

Workflow tags: `rg`.

**Trigger and harm.** A corrected and an uncorrected GWAS are compared after also changing intercept constraints.
More significant P values may arise from imposed assumptions rather than improved data.

**Legacy disposition.** Maintainer expects approximate scale cancellation in correlation but explicitly notes limited formal testing and warns against incorrect constraints.

**LDSC3 assessment — Methodological limitation.** LDSC3 has the same correlation-as-covariance-over-h2 structure.
Cancellation under constant scaling is not a guarantee for cohort-specific corrections, changed weights or wrong fixed intercepts.
Choosing constraints only for the corrected run prevents attributing a P-value difference to genomic control alone.

**Recommended prevention.** Document a like-for-like sensitivity comparison holding estimator settings fixed and distinguish absolute h2 bias from approximate correlation robustness.

**Evidence and verification.** source-inspected
[src/ldsc/_regression/estimators.py:RG](<../../../../src/ldsc/_regression/estimators.py>); [docs/reference/regression-options.md](<../../../reference/regression-options.md>); [Source](<https://github.com/bulik/ldsc/issues/82#issuecomment-300568311>).

<a id="issue-87"></a>

### #87: The value after --intercept-gencov flag is the correct intercept term

[Original issue](https://github.com/bulik/ldsc/issues/87) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **2**.

Workflow tags: `rg`.

**Trigger and harm.** User knows phenotype overlap and seeks the old constrain-intercept option. A renamed flag can conceal the strength of the scientific assumption.

**Legacy disposition.** Maintainer confirms --intercept-gencov and recommends checking the unconstrained fit despite known overlap.

**LDSC3 assessment — Methodological limitation.** LDSC3 keeps an explicit covariance-intercept option and documents fixed versus estimated intercepts.
Sample-overlap contribution is one component of the intercept; known Ns and phenotype correlation do not rule out other confounding or model mismatch.
There is no need to restore the obsolete name.

**Recommended prevention.** Use one canonical option name, explain its scale and assumptions, and recommend comparison with the default unconstrained fit.

**Evidence and verification.** source-inspected
[src/ldsc/_regression/runner.py:add_rg_arguments](<../../../../src/ldsc/_regression/runner.py>); [docs/reference/regression-options.md](<../../../reference/regression-options.md>).

<a id="issue-88"></a>

### #88: why the estimated heritablity is negative?

[Original issue](https://github.com/bulik/ldsc/issues/88) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **1**.

Workflow tags: `h2`, `rg`.

**Trigger and harm.** A noisy estimated h2 is negative and user proposes direct computation from marginal betas.
Population formulas and noisy marginal GWAS effects are conflated, producing attenuated or invalid correlation.

**Legacy disposition.** Maintainer explains unbounded estimates and the difference between conditional joint effects and marginal GWAS effects.

**LDSC3 assessment — Methodological limitation.** LDSC3 preserves unconstrained h2 and cannot manufacture a valid correlation denominator from unsupported signal.
Its error/status handling should explain this inferential boundary.
Replacing h2 with its absolute value, zero or a direct beta correlation would change the estimand and ignore LD and sampling error.

**Recommended prevention.** State the population versus estimator distinction and offer covariance reporting when correlation is unavailable, with matched-SNP uncertainty.

**Evidence and verification.** source-inspected
[src/ldsc/_regression/estimators.py:RG](<../../../../src/ldsc/_regression/estimators.py>); [docs/user-guide/workflows/genetic-correlation.md](<../../../user-guide/workflows/genetic-correlation.md>).

<a id="issue-89"></a>

### #89: why are the genetic correlation out of the bound of [-1, 1] and have very large se?

[Original issue](https://github.com/bulik/ldsc/issues/89) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **3**.

Workflow tags: `rg`.

**Trigger and harm.** Small/low-heritability GWAS or dominant loci yield imprecise correlations, including fully overlapping samples.
Users may interpret an out-of-range point estimate as software malfunction or remove overlap unnecessarily.

**Legacy disposition.** Maintainer explains unconstrained sampling estimates and suggests individual-level bivariate GREML when available, with caveats for excluded loci.

**LDSC3 assessment — Methodological limitation.** LDSC3 still has an unconstrained ratio estimator and free covariance intercept.
Full sample overlap is allowed by the model; it does not supply more polygenic information or guarantee a stable denominator.
Extreme-locus exclusions must be reported as target-changing sensitivity analyses.

**Recommended prevention.** Explain uncertainty first, audit retained SNPs and per-trait signal, and describe individual-level alternatives without promising they solve every design.

**Evidence and verification.** source-inspected
[src/ldsc/_regression/estimators.py:RG](<../../../../src/ldsc/_regression/estimators.py>); [docs/user-guide/workflows/genetic-correlation.md](<../../../user-guide/workflows/genetic-correlation.md>); [Source](<https://github.com/bulik/ldsc/issues/89#issuecomment-398883655>).

<a id="issue-92"></a>

### #92: Z-scores New Categories

[Original issue](https://github.com/bulik/ldsc/issues/92) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **3**.

Workflow tags: `ldsc-seg`.

**Trigger and harm.** User plots Z-derived P values for 12 custom epigenomic annotations. A significance plot can be labelled enrichment while testing a conditional contribution.

**Legacy disposition.** Maintainer confirms coefficient tests are appropriate for the stated contribution question, conditional on a sensible baseline.

**LDSC3 assessment — Methodological limitation.** LDSC3 already distinguishes baseline-only enrichment from baseline-plus-one-query coefficient tests.
Each query is fit separately, so these P values do not compare queries directly or establish their mutual independence.
Baseline redundancy and multiple testing across the requested family still matter.

**Recommended prevention.** State the null hypothesis and conditioning model in plot captions, preserve sign and effect magnitude, and report multiplicity across all tested queries.

**Evidence and verification.** source-inspected
[docs/reference/results/partitioned-h2.md:11-30](<../../../reference/results/partitioned-h2.md#two-regimes>); [src/ldsc/_regression/runner.py:_fit_partitioned_query](<../../../../src/ldsc/_regression/runner.py>).

<a id="issue-94"></a>

### #94: Genetic Correlation: nan (nan) (h2  out of bounds)

[Original issue](https://github.com/bulik/ldsc/issues/94) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **3**.

Workflow tags: `h2`, `rg`.

**Trigger and harm.** One trait has h2=-.0057, SE=.0915 and mean chi-square about one.
Users may interpret missing correlation as a file-format failure instead of lack of information.

**Legacy disposition.** Maintainer identifies absent SNP-heritability evidence and recommends checking parsing only if inconsistent with independent GWAS diagnostics.

**LDSC3 assessment — Methodological limitation.** LDSC3 sets correlation unavailable when a full h2 total is nonpositive, while retaining standalone estimates.
This is appropriate inferential scope; the separate delete-block boundary needs improved handling. A valid phenotype can have an uninformative dataset without implying zero biological heritability.

**Recommended prevention.** Explain unavailable correlation with per-trait signal/SE and pair SNP counts; do not recommend absolute-value or constrained-h2 workarounds.

**Evidence and verification.** source-inspected
[src/ldsc/_regression/estimators.py:RG.__init__](<../../../../src/ldsc/_regression/estimators.py>); [src/ldsc/_regression/runner.py:estimate_rg_pairs](<../../../../src/ldsc/_regression/runner.py>).

<a id="issue-95"></a>

### #95: Effective sample size or total N?

[Original issue](https://github.com/bulik/ldsc/issues/95) · Source state: **closed** · Severity: **high** · Accessible comments reviewed: **3**.

Workflow tags: `munge-sumstats`, `h2`, `rg`.

**Trigger and harm.** User proposes effective N to compensate for control imbalance while considering liability conversion.
Mixing N conventions with sample prevalence can double-correct ascertainment and rescale heritability.

**Legacy disposition.** Maintainer says total N plus sample/population prevalence handles the stated imbalance; later question is unanswered.

**LDSC3 assessment — Methodological limitation.** LDSC3 inherits explicit total/case-control N handling and separate prevalence conversion.
This is not a general endorsement for all modern meta-analysis or mixed-model effective-N conventions.
The report should specify the statistic model, N definition and prevalence combination together, because one factor cannot be changed in isolation.

**Recommended prevention.** Add a worked imbalanced case-control example and flag effective-N meta-analysis as requiring source-specific guidance; avoid automatic N reinterpretation.

**Evidence and verification.** source-inspected
[src/ldsc/_sumstats/munge_kernel.py:process_n](<../../../../src/ldsc/_sumstats/munge_kernel.py>); [src/ldsc/_regression/prevalence.py](<../../../../src/ldsc/_regression/prevalence.py>); [docs/reference/munge-sumstats.md](<../../../reference/munge-sumstats.md>).

<a id="issue-102"></a>

### #102: Partitioned heritability by functional category

[Original issue](https://github.com/bulik/ldsc/issues/102) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **1**.

Workflow tags: `baseline-partitioned-h2`, `continuous`.

**Trigger and harm.** An annotation is enriched only with or without a surrounding extension; user misstates 500bp as 500kb.
The biological claim may be assigned to the wrong genomic region or test.

**Legacy disposition.** Maintainer describes +500bp primarily as a control and explains overlapping marginal enrichment versus conditional coefficients, including continuous generalization.

**LDSC3 assessment — Methodological limitation.** LDSC3 exposes explicit padding in base pairs and documents query/baseline tests.
Padding changes annotation membership and the scientific hypothesis; significance of the extended region does not identify the original element as causal. A continuous coefficient is per annotation unit, not an ordinary binary proportion.

**Recommended prevention.** Label exact padding units and tested regions, explain original/extension overlap, and report coefficient tests for conditional questions rather than infer location from a P-value change.

**Evidence and verification.** source-inspected
[docs/reference/annotate.md](<../../../reference/annotate.md>); [docs/reference/results/partitioned-h2.md](<../../../reference/results/partitioned-h2.md>); [docs/user-guide/workflows/continuous-annotations.md](<../../../user-guide/workflows/continuous-annotations.md>).

<a id="issue-108"></a>

### #108: Validation of ldsc results

[Original issue](https://github.com/bulik/ldsc/issues/108) · Source state: **closed** · Severity: **high** · Accessible comments reviewed: **2**.

Workflow tags: `ldsc-seg`.

**Trigger and harm.** User plots enrichment_p when the intended figure tests conditional coefficients.
Many apparently significant cell types can be falsely presented as independent contributions.

**Legacy disposition.** User confirms coefficient Z reproduces the desired pattern.
Their claim enrichment does not account for overlap is imprecise: overlap-aware enrichment and conditional coefficient tests answer different questions.

**LDSC3 assessment — Methodological limitation.** LDSC3 computes overlap-aware proportions and explicitly headlines coefficient_p for query models.
It must not teach that enrichment ignores overlaps; rather, an annotation can inherit marginal signal through overlap while its conditional coefficient adds none. A reproducibility check must match model, resource version, hypothesis and multiple-testing rule.

**Recommended prevention.** Use a concrete conditional-versus-overlap-aware example and correct the causal explanation while preserving the user-confirmed metric distinction.

**Evidence and verification.** source-inspected
[docs/reference/results/partitioned-h2.md:8](<../../../reference/results/partitioned-h2.md#partitioned-heritability-result-columns>); [docs/reference/results/partitioned-h2.md:76-83](<../../../reference/results/partitioned-h2.md#choosing-the-right-column>); [src/ldsc/_regression/overlap_model.py](<../../../../src/ldsc/_regression/overlap_model.py>).

<a id="issue-123"></a>

### #123: SNP-heritability calculation requires consecutive chromosome number?

[Original issue](https://github.com/bulik/ldsc/issues/123) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **1**.

Workflow tags: `h2`, `ldscore`.

**Trigger and harm.** GWAS/genotypes cover nonconsecutive chromosomes or selected genomic regions. A successful partial-input run can extrapolate selected signal to an inappropriate genome-wide SNP count.

**Legacy disposition.** Maintainer explains the h2/M extrapolation and calls region-restricted count substitution untested; GREML is suggested for observed-variant targets.

**LDSC3 assessment — Methodological limitation.** LDSC3 allows explicitly selected chromosome scopes for direct LD scores while full @ suites and gene LD-score indexes require their declared completeness.
These technical contracts do not certify that a restricted set represents the genome.
Count metadata define a target but cannot validate representative sampling or a local-heritability method.

**Recommended prevention.** Label subset analyses by their actual reference/count scope and explicitly warn against interpreting arbitrary chromosome subsets as validated genome-wide or local heritability.

**Evidence and verification.** source-inspected
[docs/reference/ldscore-snp-universes.md](<../../../reference/ldscore-snp-universes.md>); [src/ldsc/_ldscore/preflight.py:validate_direct_scope](<../../../../src/ldsc/_ldscore/preflight.py>); [docs/reference/gene-ldscore-index.md](<../../../reference/gene-ldscore-index.md>).

<a id="issue-129"></a>

### #129: Questions about partitioned heritability, and misc.

[Original issue](https://github.com/bulik/ldsc/issues/129) · Source state: **closed** · Severity: **high** · Accessible comments reviewed: **1**.

Workflow tags: `ldsc-seg`, `baseline-partitioned-h2`, `h2`.

**Trigger and harm.** Users want enrichment of disease differential-expression lists and compare simple versus baseline models.
Association with a derived annotation may be misread as disease causality or a direct comparison between states.

**Legacy disposition.** Maintainer permits arbitrary gene lists but says disease-state lists were untested; explains coefficient-focused cell-type scans and recommends Phase 3 except historical replication.

**LDSC3 assessment — Methodological limitation.** LDSC3 can generate query annotations from gene lists and fit each on baseline controls.
It does not validate the upstream expression contrast, remove disease/treatment confounding, or test one query differs from another.
Different total h2 from richer models may reflect changed assumptions and SNP sets, not a programming error.

**Recommended prevention.** Define gene-selection universe and comparison hypothesis before interpretation, use matched controls, and distinguish capability from validated causal inference.

**Evidence and verification.** source-inspected
[docs/user-guide/workflows/gene-list-ld-scores.md](<../../../user-guide/workflows/gene-list-ld-scores.md>); [src/ldsc/_regression/runner.py:_fit_partitioned_query](<../../../../src/ldsc/_regression/runner.py>); [docs/reference/results/partitioned-h2.md](<../../../reference/results/partitioned-h2.md>).

<a id="issue-151"></a>

### #151: Heritability become less after imputation

[Original issue](https://github.com/bulik/ldsc/issues/151) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **1**.

Workflow tags: `h2`.

**Trigger and harm.** Directly genotyped versus imputed analyses use different regression-SNP sets and h2 estimates .5 versus .3.
Users assume LDSC is an in-sample variance-explained fit monotone in the number of markers.

**Legacy disposition.** Open and unresolved; no evidence of a code bug is supplied.

**LDSC3 assessment — Methodological limitation.** LDSC3 estimates a regression slope scaled by a reference-SNP count, not monotone prediction R² as rows are added.
Imputation changes representativeness, quality, sample size and LD tagging of regression observations. A credible comparison holds reference/count/model policy fixed, examines retained rows, and compares uncertainty rather than raw point estimates.

**Recommended prevention.** Explain the estimand with a before/after SNP-flow table and inspect input calibration; do not treat higher h2 as automatic improvement.

**Evidence and verification.** source-inspected
[src/ldsc/_regression/estimators.py:Hsq](<../../../../src/ldsc/_regression/estimators.py>); [docs/reference/ldscore-snp-universes.md](<../../../reference/ldscore-snp-universes.md>); [docs/user-guide/workflows/heritability.md](<../../../user-guide/workflows/heritability.md>).

<a id="issue-152"></a>

### #152: -interecept-gencov

[Original issue](https://github.com/bulik/ldsc/issues/152) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **1**.

Workflow tags: `rg`.

**Trigger and harm.** Two diseases share controls but no cases, and user wants a constraint to reduce SE.
Using a quantitative overlap formula or choosing the intercept for significance can bias correlation.

**Legacy disposition.** No maintainer answer; formula interpretation remains unresolved in this issue.

**LDSC3 assessment — Methodological limitation.** LDSC3 can estimate or explicitly fix a covariance intercept but cannot infer binary ascertainment from case/control totals alone.
Shared controls induce correlation under specific ascertainment assumptions.
Free intercept estimation remains the practical default; a large SE is evidence of uncertainty, not justification for an arbitrary constraint.

**Recommended prevention.** Explain the binary-trait overlap derivation from the original supplement and require appropriate joint sample counts/assumptions before presenting a fixed-intercept calculation.

**Evidence and verification.** source-inspected
[src/ldsc/_regression/runner.py:add_rg_arguments](<../../../../src/ldsc/_regression/runner.py>); [src/ldsc/_regression/estimators.py:Gencov](<../../../../src/ldsc/_regression/estimators.py>); [docs/reference/regression-options.md](<../../../reference/regression-options.md>).

<a id="issue-155"></a>

### #155: negative h2 and Ratio &gt; 1

[Original issue](https://github.com/bulik/ldsc/issues/155) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **10**.

Workflow tags: `h2`.

**Trigger and harm.** Some examples have weak signal, while another has strongly inflated chi-square/intercept and negative h2.
Treating all cases as sampling noise or automatically excluding MHC can miss major calibration/model problems.

**Legacy disposition.** No maintainer resolution; one user reports MHC exclusion helps and another says it does not.

**LDSC3 assessment — Methodological limitation.** LDSC3 keeps unconstrained h2 and ratio.
Large N or millions of raw SNPs does not ensure a valid positive LD-score relationship.
The strongly negative example needs checks of N, P/statistic consistency, ancestry, imputation and reference matching; small negative estimates with large SE are a different regime.

**Recommended prevention.** Provide a diagnostic sequence distinguishing low signal from clear negative-slope misspecification and report any region-exclusion sensitivity as a changed target.

**Evidence and verification.** source-inspected
[src/ldsc/_regression/estimators.py:Hsq](<../../../../src/ldsc/_regression/estimators.py>); [src/ldsc/_regression/runner.py:_fit_h2_dataset](<../../../../src/ldsc/_regression/runner.py>); [docs/user-guide/workflows/heritability.md](<../../../user-guide/workflows/heritability.md>).

<a id="issue-156"></a>

### #156: Prop._h2 outside [0, 1]

[Original issue](https://github.com/bulik/ldsc/issues/156) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **3**.

Workflow tags: `baseline-partitioned-h2`.

**Trigger and harm.** Overlapping baseline and two extra binary annotations yield negative or above-one estimated shares.
Users may clip noisy estimates or confuse a conditional contribution with an overlap-aware attributable share.

**Legacy disposition.** Users hypothesize narrow annotations/noisy GWAS; no verified diagnosis or code fix.

**LDSC3 assessment — Methodological limitation.** LDSC3 reports unconstrained fitted quantities with jackknife SEs and distinct conditional category_h2 versus overlap-aware prop_h2.
True valid binary-region heritability is nonnegative, but estimates can escape bounds with noise/collinearity; overlap means shares need not sum to one.
Severe stable violations need model/input investigation, not automatic acceptance.

**Recommended prevention.** Explain both quantities and uncertainty, inspect support/collinearity, and never silently bound results or present all out-of-range estimates as equivalent.

**Evidence and verification.** source-inspected
[docs/reference/results/partitioned-h2.md:45-61](<../../../reference/results/partitioned-h2.md#partitioned_h2tsv>); [src/ldsc/_regression/overlap_model.py](<../../../../src/ldsc/_regression/overlap_model.py>).

<a id="issue-189"></a>

### #189: Genetic correlation significantly out of bounds [-1,1], with relatively small SE's

[Original issue](https://github.com/bulik/ldsc/issues/189) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **3**.

Workflow tags: `rg`.

**Trigger and harm.** A proxy-trait comparison returns magnitude1.232 with SE0.064 and changes markedly when intercepts are constrained. A significant nonzero test can coexist with implausible magnitude; arbitrary constraints can conceal misspecification.

**Legacy disposition.** No maintainer resolution or code fix is provided; users report the same phenomenon years later.

**LDSC3 assessment — Methodological limitation.** LDSC3 retains an unbounded ratio estimator.
Its rg p-value tests zero, not whether the parameter equals one or whether assumptions hold.
Allele reversal, ascertainment, phenotype definitions, reference mismatch and intercept assumptions require separate checks; a constrained0.90 estimate is not evidence of correction.

**Recommended prevention.** Explain zero-versus-unit nulls, confidence intervals and justified intercept sensitivity; do not clip estimates or select constraints for plausibility.

**Evidence and verification.** source-inspected
[src/ldsc/_regression/estimators.py:RG](<../../../../src/ldsc/_regression/estimators.py>); [src/ldsc/_regression/runner.py:_fit_rg_dataset](<../../../../src/ldsc/_regression/runner.py>); [docs/reference/regression-options.md](<../../../reference/regression-options.md>).

<a id="issue-191"></a>

### #191: which N to use for large case/control imbalance data

[Original issue](https://github.com/bulik/ldsc/issues/191) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **1**.

Workflow tags: `munge-sumstats`, `h2`, `rg`.

**Trigger and harm.** Disease GWAS with case fraction below0.01 raises total-versus-effective-N choice.
Wrong N rescales SNP heritability and can invalidate liability conversion or weighting.

**Legacy disposition.** A commenter points to issue95; the thread itself supplies no universal validated formula or code fix.

**LDSC3 assessment — Methodological limitation.** LDSC3 now distinguishes direct per-variant N from paired case/control normalization and documents that its legacy rule differs from harmonic effective sample size.
The data producer and association model determine the appropriate N; neither a total count nor NEFF is automatically correct for every imbalanced design.

**Recommended prevention.** Give model-specific examples using producer-documented N and matched prevalence; never silently substitute an effective-N formula.

**Evidence and verification.** source-inspected
[src/ldsc/_sumstats/munge_input.py:_resolve_sample_size_column_strategy](<../../../../src/ldsc/_sumstats/munge_input.py>); [src/ldsc/_sumstats/munge_kernel.py:process_n](<../../../../src/ldsc/_sumstats/munge_kernel.py>); [docs/user-guide/workflows/munge-sumstats.md#sample-size-columns](<../../../user-guide/workflows/munge-sumstats.md#sample-size-columns>).

<a id="issue-195"></a>

### #195: heritability estimation with summary data of mixed model

[Original issue](https://github.com/bulik/ldsc/issues/195) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **0**.

Workflow tags: `h2`.

**Trigger and harm.** Animal GWAS summary statistics come from a single-SNP mixed linear model rather than ordinary regression.
Association-model calibration, relatedness and LD reference differences can violate the regression relationship.

**Legacy disposition.** Open without a resolution or supporting validation study in the thread.

**LDSC3 assessment — Methodological limitation.** LDSC3 can read numeric summary statistics but parsing is not validation of LDSC for an animal mixed-model study.
Human HM3 defaults, hg19/hg38 assumptions and LD-score references make this use case scientifically and mechanically outside the ordinary documented workflow.

**Recommended prevention.** Document human-reference assumptions and distinguish supported file ingestion from validated study designs; seek model-specific primary evidence before use.

**Evidence and verification.** source-inspected
[src/ldsc/_regression/estimators.py:Hsq](<../../../../src/ldsc/_regression/estimators.py>); [docs/user-guide/workflows/heritability.md](<../../../user-guide/workflows/heritability.md>).

<a id="issue-200"></a>

### #200: Celltype/tissue specific analysis on ATAC seq data

[Original issue](https://github.com/bulik/ldsc/issues/200) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **1**.

Workflow tags: `ldsc-seg`.

**Trigger and harm.** Users have MACS2 narrowPeak files and ask about t-statistics, controls and replicate merging.
Treating every peak as cell-type specificity can conflate assay accessibility with the intended conditional biological test.

**Legacy disposition.** Open; no validated answer or code fix is given.

**LDSC3 assessment — Methodological limitation.** LDSC3 can project BED intervals and test each query on explicit baseline annotations; it does not infer a scientifically appropriate control from peak files.
LDSC-SEG expression-specific t-statistics need not be imitated mechanically for ATAC-seq.
Replicate aggregation and an all-accessible-regions control must follow the experimental design.

**Recommended prevention.** Add an ATAC-seq worked example separating peak preparation, query definition, common accessibility control and baseline adjustment.

**Evidence and verification.** source-inspected
[docs/user-guide/workflows/cell-type-enrichment.md](<../../../user-guide/workflows/cell-type-enrichment.md>); [src/ldsc/_annotation/interval_projection.py:ChromosomeProjector.project](<../../../../src/ldsc/_annotation/interval_projection.py>).

<a id="issue-203"></a>

### #203: Low ratio values - are the genetic correlation results reliable?

[Original issue](https://github.com/bulik/ldsc/issues/203) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **2**.

Workflow tags: `h2`, `rg`.

**Trigger and harm.** Both traits are GC corrected, one has low SNP heritability and another an intercept0.8095. A negative ratio is mistaken for proof of reliable genetic correlation despite weak signal and altered calibration.

**Legacy disposition.** Only user speculation is provided; no maintainer guarantee or code fix.

**LDSC3 assessment — Methodological limitation.** LDSC3 retains intercept and ratio summaries and an unbounded genetic-correlation ratio. A low or negative ratio does not validate the analysis when mean chi-square is near1, and low denominator precision can remain consequential even if uniform scaling cancels algebraically in genetic correlation.

**Recommended prevention.** Explain GC correction, ratio denominator instability and SNP-heritability precision; recommend checking uncorrected statistics and justified intercept settings.

**Evidence and verification.** source-inspected
[src/ldsc/_regression/estimators.py:RG](<../../../../src/ldsc/_regression/estimators.py>); [src/ldsc/_regression/runner.py:_fit_rg_dataset](<../../../../src/ldsc/_regression/runner.py>); [docs/reference/regression-options.md](<../../../reference/regression-options.md>).

<a id="issue-204"></a>

### #204: Continuous annotation for a subset of SNP: is it reasonable to set other SNP=0?

[Original issue](https://github.com/bulik/ldsc/issues/204) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **0**.

Workflow tags: `continuous`.

**Trigger and harm.** SDS exists for only four million SNPs; the remainder would be filled with biologically meaningful0.
Replacing unknown scores with zero changes the annotation and can bias interpretation.

**Legacy disposition.** Open without resolution.

**LDSC3 assessment — Methodological limitation.** LDSC3 numeric annotation parsing rejects missing cells rather than silently filling them.
Its quantile target interface permits one explicitly declared missing token and excludes those SNPs from quantile construction.
That post-fit facility does not justify zero-imputing missing values in a fitted annotation.

**Recommended prevention.** Document measured zero versus unknown support; require an explicit scientifically defensible missingness/restriction strategy before generating LD scores.

**Evidence and verification.** source-inspected
[docs/user-guide/workflows/continuous-annotations.md](<../../../user-guide/workflows/continuous-annotations.md>); [docs/reference/results/quantile-h2.md](<../../../reference/results/quantile-h2.md>); [docs/reference/results/partitioned-h2.md](<../../../reference/results/partitioned-h2.md>).

<a id="issue-210"></a>

### #210: doubt with sub-annotations - cell-type enrichment

[Original issue](https://github.com/bulik/ldsc/issues/210) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **0**.

Workflow tags: `ldsc-seg`.

**Trigger and harm.** User restricts baseline toCONC then tests CONC intersect Cahoy and asks which baseline is valid.
Changing the reference-SNP universe and controls changes the scientific estimand and can omit necessary adjustment.

**Legacy disposition.** Open without a supported resolution.

**LDSC3 assessment — Methodological limitation.** LDSC3 separates reference SNPs, regression SNPs and query annotations, but cannot infer whether the target is genome-wide cell-type signal or an interaction withinCONC.
An intersection query with an appropriateCONC control differs from deleting all nonCONC baseline support.

**Recommended prevention.** Specify the contrast first, retain a scientifically justified genome-wide baseline and add explicit region/control terms rather than silently truncating SNP support.

**Evidence and verification.** source-inspected
[docs/reference/results/partitioned-h2.md](<../../../reference/results/partitioned-h2.md>); [src/ldsc/_regression/runner.py:summarize_partitioned_h2](<../../../../src/ldsc/_regression/runner.py>).

<a id="issue-213"></a>

### #213: Effect of --n-blocks on partitioning heritability

[Original issue](https://github.com/bulik/ldsc/issues/213) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **1**.

Workflow tags: `baseline-partitioned-h2`.

**Trigger and harm.** Users increase block count from200to1000or2000 to avoid singular delete-block fits.
An arbitrary numerical workaround changes uncertainty estimates and can create misleading p-values.

**Legacy disposition.** Users observe near-identical point enrichment but changing SE/p-values; no validated universal block-count fix.

**LDSC3 assessment — Methodological limitation.** LDSC3 retains contiguous block jackknife with explicit --num-blocks and detailed failed-block rank/support diagnostics.
More smaller blocks can retain rare annotation support after deletion, but that is not evidence of calibrated uncertainty or adequate effective independent regions.

**Recommended prevention.** Explain block-size/LD dependence and sparse-query support; inspect failed blocks before a prespecified sensitivity analysis.

**Evidence and verification.** source-inspected
[docs/reference/results/partitioned-h2.md](<../../../reference/results/partitioned-h2.md>); [src/ldsc/_regression/runner.py:summarize_partitioned_h2](<../../../../src/ldsc/_regression/runner.py>).

<a id="issue-218"></a>

### #218: Doubt about genetic correlation

[Original issue](https://github.com/bulik/ldsc/issues/218) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **0**.

Workflow tags: `rg`.

**Trigger and harm.** Two phenotype GWAS use the same individuals.
Forcing a zero covariance intercept can bias shared-sample analyses.

**Legacy disposition.** Open with no response or code fix.

**LDSC3 assessment — Methodological limitation.** LDSC3 fits a genetic-covariance intercept and permits explicit constraints. A shared sample does not automatically prohibit genetic-correlation estimation; the overlap contribution also depends on phenotypic covariance and study design.
Knowing the overlap count alone is not a universal justification for an arbitrary fixed intercept.

**Recommended prevention.** Explain how the covariance intercept handles overlap and when external information justifies constraints.

**Evidence and verification.** source-inspected
[src/ldsc/_regression/estimators.py:RG](<../../../../src/ldsc/_regression/estimators.py>); [src/ldsc/_regression/runner.py:_fit_rg_dataset](<../../../../src/ldsc/_regression/runner.py>); [docs/reference/regression-options.md](<../../../reference/regression-options.md>).

<a id="issue-222"></a>

### #222: how to understand attenuation ratio?

[Original issue](https://github.com/bulik/ldsc/issues/222) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **0**.

Workflow tags: `h2`.

**Trigger and harm.** A7000-person case/control GWAS has ratio0.3475 while intercept is1.0291 and mean chi-square1.0839. A modest absolute intercept deviation can yield a large ratio when total inflation is small.

**Legacy disposition.** Open without a validated conclusion.

**LDSC3 assessment — Methodological limitation.** LDSC3 computes ratio=(intercept-1)/(mean chi-square-1).
It quantifies a fitted decomposition under LDSC assumptions, not a direct binary test for population stratification, and becomes unstable near mean chi-square1.
The intercept SE and study calibration remain necessary context.

**Recommended prevention.** Show the ratio formula and denominator sensitivity; avoid declaring acceptable/confounded GWAS from a universal ratio cutoff.

**Evidence and verification.** source-inspected
[src/ldsc/_regression/estimators.py:Hsq](<../../../../src/ldsc/_regression/estimators.py>); [docs/reference/regression-options.md](<../../../reference/regression-options.md>).

<a id="issue-229"></a>

### #229: heritability greater than 1

[Original issue](https://github.com/bulik/ldsc/issues/229) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **1**.

Workflow tags: `h2`.

**Trigger and harm.** One trait estimates0.0325(SE0.2509), another1.1843(SE0.2387).
Very noisy unbounded estimates may be interpreted as biological proportions.

**Legacy disposition.** Open with no established software defect or fix.

**LDSC3 assessment — Methodological limitation.** LDSC3 retains unbounded regression estimates, so values below0or above1 can occur by sampling error.
Here both studies have mean chi-square near1 and low sample size; the uncertainty is central.
Observed-versus-liability scale and N/count errors require separate checks before interpreting an out-of-range estimate.

**Recommended prevention.** Show estimate with SE/confidence interval, flag weak signal and explain scaling; avoid clipping or reporting the point estimate alone.

**Evidence and verification.** source-inspected
[src/ldsc/_regression/estimators.py:Hsq](<../../../../src/ldsc/_regression/estimators.py>); [docs/user-guide/workflows/heritability.md](<../../../user-guide/workflows/heritability.md>).

<a id="issue-249"></a>

### #249: how to interpret genetic correction of two traits that are phenotypically correlated

[Original issue](https://github.com/bulik/ldsc/issues/249) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **0**.

Workflow tags: `rg`.

**Trigger and harm.** Same-cohort BMI/T2D have phenotypic correlation0.4 and genetic correlation0.3; user considers covariate adjustment.
Phenotype conditioning can change the estimand or induce bias rather than isolate genetic sharing.

**Legacy disposition.** Open without a supported answer or fix.

**LDSC3 assessment — Methodological limitation.** LDSC3 estimates shared genetic effects, a quantity distinct from phenotypic correlation.
Same-sample covariance is addressed through the cross-trait intercept under assumptions. A positive phenotypic correlation does not make genetic correlation meaningless, and conditioning one GWAS on the other phenotype is a different scientific analysis.

**Recommended prevention.** Explain variance decomposition, overlap intercept and the consequences of heritable covariate adjustment; avoid recommending residualization as automatic QC.

**Evidence and verification.** source-inspected
[src/ldsc/_regression/estimators.py:RG](<../../../../src/ldsc/_regression/estimators.py>); [src/ldsc/_regression/runner.py:_fit_rg_dataset](<../../../../src/ldsc/_regression/runner.py>); [docs/reference/regression-options.md](<../../../reference/regression-options.md>).

<a id="issue-251"></a>

### #251: How to set N sample size on results from METAL while using munge_sumstats.py

[Original issue](https://github.com/bulik/ldsc/issues/251) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **0**.

Workflow tags: `munge-sumstats`.

**Trigger and harm.** Sample-size-weighted METAL output has a Weight column summing three cohort totals.
Mistaking arbitrary meta-analysis weights for N changes the SNP-heritability scale.

**Legacy disposition.** Open with no established answer or fix.

**LDSC3 assessment — Methodological limitation.** LDSC3 permits --N-col Weight, but cannot verify that METAL Weight represents the N associated with the reported Z statistic.
Under the exact sample-size-weighting specification it may encode the intended participant totals; missing cohorts, case/control imbalance and other METAL schemes require producer-specific interpretation.

**Recommended prevention.** Document the METAL scheme and use per-variant participation/effective-N information when appropriate; do not select a constant from the abstract blindly.

**Evidence and verification.** source-inspected
[src/ldsc/_sumstats/munge_input.py:_resolve_sample_size_column_strategy](<../../../../src/ldsc/_sumstats/munge_input.py>); [src/ldsc/_sumstats/munge_kernel.py:process_n](<../../../../src/ldsc/_sumstats/munge_kernel.py>); [docs/user-guide/workflows/munge-sumstats.md#sample-size-columns](<../../../user-guide/workflows/munge-sumstats.md#sample-size-columns>).

<a id="issue-252"></a>

### #252: how to understand negative heritability or heritability values that are greater than 1 or even 2

[Original issue](https://github.com/bulik/ldsc/issues/252) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **2**.

Workflow tags: `h2`, `rg`.

**Trigger and harm.** SNP heritability estimates range from negative to2.68; later discussion reveals N400.
Users may treat unreliable point estimates or a method switch as a correction.

**Legacy disposition.** Comments cite alternative literature/HDL, but no validated fix of these data is shown.

**LDSC3 assessment — Methodological limitation.** LDSC3 retains the same statistical information limits; a successful fit cannot create power in a very small GWAS.
Negative or above-one regression estimates must be interpreted with uncertainty, scale and model checks.
Alternative methods need their own assumptions and validation, not selection because estimates look plausible.

**Recommended prevention.** Explain low-power limitations and report SE/signal diagnostics; separate alternative-method evidence from an established remedy.

**Evidence and verification.** source-inspected
[src/ldsc/_regression/estimators.py:Hsq](<../../../../src/ldsc/_regression/estimators.py>); [docs/user-guide/workflows/heritability.md](<../../../user-guide/workflows/heritability.md>).

<a id="issue-255"></a>

### #255: rg was out of bounds.

[Original issue](https://github.com/bulik/ldsc/issues/255) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **4**.
The API reports 5 comments; one was unavailable through both retrieval methods.

Workflow tags: `munge-sumstats`, `h2`, `rg`.

**Trigger and harm.** Several low-powered binary GWAS and their meta-analysis return undefined/out-of-bounds correlations and low mean chi-square. N/prevalence mismatch and weak denominators can make comparisons misleading; cohort-versus-meta comparisons also overlap.

**Legacy disposition.** No resolution is visible; GitHub reports one unavailable comment, so full-thread evidence is incomplete.

**LDSC3 assessment — Methodological limitation.** LDSC3 exposes sample-size strategies, prevalence conversion and pairwise diagnostics, but cannot make these low-signal data informative. --samp-prev belongs to regression/scale conversion, not raw munging.
Effective-N0.5 conventions require a consistent derivation; low mean chi-square is not high cross-study correlation.

**Recommended prevention.** Explain model-matched N/prevalence, overlapping meta-analysis samples and weak-signal limits; disclose the unavailable comment.

**Evidence and verification.** source-inspected
[src/ldsc/_regression/estimators.py:RG](<../../../../src/ldsc/_regression/estimators.py>); [src/ldsc/_regression/runner.py:_fit_rg_dataset](<../../../../src/ldsc/_regression/runner.py>); [docs/reference/regression-options.md](<../../../reference/regression-options.md>).

<a id="issue-257"></a>

### #257: 2 architecture questions on LDSC

[Original issue](https://github.com/bulik/ldsc/issues/257) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **1**.

Workflow tags: `h2`, `rg`.

**Trigger and harm.** User conflates phenotypic with genetic correlation and asks why chi-square=N beta squared.
An approximation can be applied to unstandardized effects or interpreted causally.

**Legacy disposition.** Open with follow-up derivation questions; no supported fix.

**LDSC3 assessment — Methodological limitation.** LDSC3 numerical inputs do not imply that arbitrary reported beta satisfies chi-square=N beta squared.
The identity depends on standardized variables and the association model/SE approximation.
Genetic correlation concerns genetic effects and differs from phenotypic correlation even when phenotypes are strongly correlated.

**Recommended prevention.** Provide a compact standardized-model derivation with assumptions and distinguish genetic covariance, phenotypic covariance and causal effects.

**Evidence and verification.** source-inspected
[docs/reference/regression-options.md](<../../../reference/regression-options.md>); [src/ldsc/_regression/estimators.py:LD_Score_Regression](<../../../../src/ldsc/_regression/estimators.py>).

<a id="issue-259"></a>

### #259: Rg &gt;1?   Rg=0.84 but P non-significant

[Original issue](https://github.com/bulik/ldsc/issues/259) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **0**.

Workflow tags: `rg`.

**Trigger and harm.** A user reports rg1.93 and rg0.84 with p0.1768.
Magnitude is confused with precision and validity.

**Legacy disposition.** Open without resolution.

**LDSC3 assessment — Methodological limitation.** LDSC3 retains an unbounded genetic-correlation ratio and tests zero using its SE.
An estimate0.84 can be nonsignificant if uncertainty is large; a value above1 is not a biologically valid population correlation but can arise from noise or misspecification.
No clipping can restore information.

**Recommended prevention.** Display SE/confidence intervals beside magnitude and state the zero-null test; check denominator SNP heritabilities and model assumptions.

**Evidence and verification.** source-inspected
[src/ldsc/_regression/estimators.py:RG](<../../../../src/ldsc/_regression/estimators.py>); [src/ldsc/_regression/runner.py:_fit_rg_dataset](<../../../../src/ldsc/_regression/runner.py>); [docs/reference/regression-options.md](<../../../reference/regression-options.md>).

<a id="issue-265"></a>

### #265: Low polygenicity signal

[Original issue](https://github.com/bulik/ldsc/issues/265) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **2**.

Workflow tags: `h2`.

**Trigger and harm.** More than one million SNPs yield mean chi-square1.0694 and intercept1.0555.
Successful preprocessing and many SNPs are mistaken for sufficient heritable signal.

**Legacy disposition.** Open; reporter says unresolved and no fix is identified.

**LDSC3 assessment — Methodological limitation.** LDSC3 retains these diagnostics but the fitted polygenic component is the small slope-associated part, not total inflation. A high ratio can reflect confounding, model mismatch or weak signal and is not a direct diagnosis by itself; SNP-heritability SE and ancestry/QC evidence remain necessary.

**Recommended prevention.** Explain signal versus inflation and ratio uncertainty; add an example where many valid SNPs still provide insufficient information.

**Evidence and verification.** source-inspected
[src/ldsc/_regression/estimators.py:Hsq](<../../../../src/ldsc/_regression/estimators.py>); [docs/user-guide/workflows/heritability.md](<../../../user-guide/workflows/heritability.md>).

<a id="issue-271"></a>

### #271: interpreting --a1-inc

[Original issue](https://github.com/bulik/ldsc/issues/271) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **0**.

Workflow tags: `munge-sumstats`, `baseline-partitioned-h2`, `ldsc-seg`.

**Trigger and harm.** Two tutorials use the same BMI data with different --a1-inc settings.
All-positive Z is harmless for chi-square-only fitting but invalidates ordinary signed genetic-correlation interpretation.

**Legacy disposition.** Open without clarification.

**LDSC3 assessment — Methodological limitation.** LDSC3 explicitly documents --a1-inc as an assertion that A1 is always increasing and bypasses direction columns.
It produces positive Z from P; it does not infer which allele increases the trait.
Differences in tutorial source allele coding matter, and a file produced this way should not be casually reused for genetic correlation.

**Recommended prevention.** Explain when effect directions are unnecessary and label positive-Z outputs; require genuine signed effects for genetic-correlation workflows.

**Evidence and verification.** source-inspected
[src/ldsc/_sumstats/munge_input.py:prepare_munge_input](<../../../../src/ldsc/_sumstats/munge_input.py>); [src/ldsc/_sumstats/munge_kernel.py:munge_sumstats](<../../../../src/ldsc/_sumstats/munge_kernel.py>).

<a id="issue-283"></a>

### #283: munge_sumstats.py error: questions for input file

[Original issue](https://github.com/bulik/ldsc/issues/283) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **0**.

Workflow tags: `munge-sumstats`, `rg`.

**Trigger and harm.** GMMAT provides SCORE/PVAL/alleles but no beta,ORorZ.
An unsigned test statistic cannot recover effect direction needed for genetic correlation.

**Legacy disposition.** Open without a validated mapping.

**LDSC3 assessment — Methodological limitation.** LDSC3 can accept a user-named signed statistic through --signed-sumstats COLUMN,0, but only if the producer defines its sign relative to A1. A score statistic may be signed or squared; the package cannot infer that semantics. P supplies magnitude, so a validated signed score can suffice without SE.

**Recommended prevention.** Document producer-defined sign and sample-size requirements; do not use --a1-inc merely to bypass an unsigned score.

**Evidence and verification.** source-inspected
[src/ldsc/_sumstats/munge_input.py:prepare_munge_input](<../../../../src/ldsc/_sumstats/munge_input.py>); [src/ldsc/_sumstats/munge_kernel.py:munge_sumstats](<../../../../src/ldsc/_sumstats/munge_kernel.py>).

<a id="issue-288"></a>

### #288: One-sided test with --h2 flag and custom annotations

[Original issue](https://github.com/bulik/ldsc/issues/288) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **1**.

Workflow tags: `baseline-partitioned-h2`, `ldsc-seg`.

**Trigger and harm.** User thinks h2-cts gives only positive estimates because its coefficient test is one-sided.
Changing the p-value tail is mistaken for preventing negative SNP heritability.

**Legacy disposition.** Open without correction; a later user repeats the negative-enrichment concern.

**LDSC3 assessment — Methodological limitation.** LDSC3 reports both two-sided enrichment_p and one-sided coefficient_p with distinct definitions.
Both arise from unbounded fitted coefficients; selecting a one-sided alternative does not force estimates positive.
Marginal overlap-aware enrichment differs from a query conditional contribution beyond the baseline.

**Recommended prevention.** Put explicit hypotheses and allowed estimate ranges beside both columns; do not add positivity constraints merely to match the requested test tail.

**Evidence and verification.** source-inspected
[docs/reference/results/partitioned-h2.md](<../../../reference/results/partitioned-h2.md>); [src/ldsc/_regression/runner.py:summarize_partitioned_h2](<../../../../src/ldsc/_regression/runner.py>).

<a id="issue-299"></a>

### #299: LDSC with QTL summary statistics

[Original issue](https://github.com/bulik/ldsc/issues/299) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **0**.

Workflow tags: Cross-cutting / not assigned to a requested group.

**Trigger and harm.** User asks whether QTL summary statistics can replace GWAS summary statistics.
Cis-restricted or selected variants violate the ordinary genome-wide sampling/model assumptions.

**Legacy disposition.** Open without a validated supported scenario.

**LDSC3 assessment — Methodological limitation.** LDSC3 can parse appropriately formatted numeric association tables, but this does not validate LDSC for every molecular-QTL design. A genome-wide molecular phenotype with adequate sample size/reference may differ fundamentally from a cis-only significant-QTL list; scope and selection must be stated before choosing an estimator.

**Recommended prevention.** Explain genome-wide versus cis/selected designs and distinguish using QTL-derived annotations from treating QTL tests as a full GWAS.

**Evidence and verification.** source-inspected
[src/ldsc/_regression/estimators.py:Hsq](<../../../../src/ldsc/_regression/estimators.py>); [docs/user-guide/workflows/heritability.md](<../../../user-guide/workflows/heritability.md>).

<a id="issue-305"></a>

### #305: How should I run LDSC on a small annotation?

[Original issue](https://github.com/bulik/ldsc/issues/305) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **1**.

Workflow tags: `ldsc-seg`.

**Trigger and harm.** Open-chromatin queries cover below1percent; users consider100kb padding solely to meet a heuristic.
Data-dependent expansion changes the biological question and may not cure invalid small-annotation p-values.

**Legacy disposition.** Open; users cite small-annotation testing work but no supported implementation.

**LDSC3 assessment — Methodological limitation.** LDSC3 detects collinearity and singular delete-block fits, yet successful fitting does not establish calibrated inference for sparse support.
The cited1percent or0.01cM coverage heuristics need their original method definitions, not a universal hard threshold.
Padding chosen to force eligibility is not neutral.

**Recommended prevention.** Explain support across independent genomic blocks and primary-paper calibration; treat specialized small-annotation procedures as separately validated future work.

**Evidence and verification.** source-inspected
[src/ldsc/_regression/jackknife.py](<../../../../src/ldsc/_regression/jackknife.py>); [docs/reference/results/partitioned-h2.md](<../../../reference/results/partitioned-h2.md>); [Source](<https://github.com/bulik/ldsc/issues/305>).

<a id="issue-308"></a>

### #308: Very high Ratio, negative h2, and doubt about constraining the intercept

[Original issue](https://github.com/bulik/ldsc/issues/308) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **2**.

Workflow tags: `h2`.

**Trigger and harm.** A4740-person GWAS has negative SNP heritability, ratio1.59and mean chi-square1.014.
Agreement with GCTA after fixing intercept1 is mistaken for proof of validity.

**Legacy disposition.** Open with repeated users and no validated resolution.

**LDSC3 assessment — Methodological limitation.** LDSC3 preserves free/fixed intercept choices and their different information/bias tradeoffs.
Including PCs and strict QC does not prove the residual intercept equals1.
Weak polygenic signal makes the slope imprecise, and agreement between estimators with different assumptions is not an accuracy reference.

**Recommended prevention.** Explain justified intercept constraints and prespecified sensitivity; report both estimates with uncertainty rather than selecting the plausible one.

**Evidence and verification.** source-inspected
[src/ldsc/_regression/estimators.py:Hsq](<../../../../src/ldsc/_regression/estimators.py>); [docs/user-guide/workflows/heritability.md](<../../../user-guide/workflows/heritability.md>).

<a id="issue-310"></a>

### #310: Can I apply ldsc to small sample size (N&lt;1000) for genetic correlation?

[Original issue](https://github.com/bulik/ldsc/issues/310) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **0**.

Workflow tags: `rg`.

**Trigger and harm.** A600-person GWAS with SNP heritability0.984(SE0.830) is correlated with a very large GWAS. A large partner cannot rescue a near-unidentified denominator from the small study.

**Legacy disposition.** Open without a validated answer.

**LDSC3 assessment — Methodological limitation.** LDSC3 retains uncertainty from both SNP-heritability estimates; low power is not equivalent to systematic inflation, and a numerically returned correlation is not automatically reliable.
The supplied estimate is highly imprecise regardless of its large point value or absence of obvious GWAS inflation.

**Recommended prevention.** Explain denominator precision and distinguish bias from variance; avoid claiming a universal sample-size cutoff or a large partner guarantees validity.

**Evidence and verification.** source-inspected
[src/ldsc/_regression/estimators.py:RG](<../../../../src/ldsc/_regression/estimators.py>); [src/ldsc/_regression/runner.py:_fit_rg_dataset](<../../../../src/ldsc/_regression/runner.py>); [docs/reference/regression-options.md](<../../../reference/regression-options.md>).

<a id="issue-313"></a>

### #313: High ratio, and high h2 confidence interval in comparison with UKBB

[Original issue](https://github.com/bulik/ldsc/issues/313) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **0**.

Workflow tags: `h2`.

**Trigger and harm.** A4069-person LDL GWAS has much wider uncertainty than UKBB; fixing intercept yields a similar point estimate.
Users try parameter tuning to reproduce a larger study rather than acknowledge information limits.

**Legacy disposition.** Open with no validated fix; using own versus European reference does not materially improve the supplied results.

**LDSC3 assessment — Methodological limitation.** LDSC3 cannot make the small cohort as precise as a large biobank.
The original brackets contain SE, not the full confidence interval, and mean chi-square1.008 indicates weak signal.
Intercept constraints trade assumptions for precision and should not be selected by agreement with UKBB.

**Recommended prevention.** Explain sample-size/signal effects, SE-versus-CI and prespecified intercept sensitivity; do not call matching point estimates an improvement.

**Evidence and verification.** source-inspected
[src/ldsc/_regression/estimators.py:Hsq](<../../../../src/ldsc/_regression/estimators.py>); [docs/user-guide/workflows/heritability.md](<../../../user-guide/workflows/heritability.md>).

<a id="issue-318"></a>

### #318: Partitioned heritability: MAF_Adj_LLD_AFR annotation

[Original issue](https://github.com/bulik/ldsc/issues/318) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **0**.

Workflow tags: `baseline-partitioned-h2`, `continuous`.

**Trigger and harm.** A negative MAF_Adj_LLD_AFR enrichment is interpreted as different genetic architecture between Europeans and Africans. A fitted annotation effect is mistaken for trans-ancestry genetic correlation or transferability.

**Legacy disposition.** Open without a supported answer.

**LDSC3 assessment — Methodological limitation.** LDSC3 classifies nonbinary annotations and warns that ordinary binary enrichment interpretation does not apply.
The meaning of this transformed baseline annotation comes from its original definition and conditional model, not the AFR substring. A negative fitted coefficient/enrichment does not directly compare ancestry-specific genetic effects.

**Recommended prevention.** Document the annotation definition and appropriate standardized coefficient; explicitly separate LD-related annotations from cross-ancestry sharing.

**Evidence and verification.** source-inspected
[docs/user-guide/workflows/continuous-annotations.md](<../../../user-guide/workflows/continuous-annotations.md>); [docs/reference/results/quantile-h2.md](<../../../reference/results/quantile-h2.md>); [docs/reference/results/partitioned-h2.md](<../../../reference/results/partitioned-h2.md>).

<a id="issue-328"></a>

### #328: High genetic correlation, but large se and not significant, h2 is weirdly low too

[Original issue](https://github.com/bulik/ldsc/issues/328) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **2**.

Workflow tags: `h2`, `rg`.

**Trigger and harm.** Brain-volume/AD correlation0.966 has SE1.163 while brain-volume SNP heritability is 0.0144(SE0.0261). A large point correlation is mistaken for shared genetic evidence despite an imprecise denominator.

**Legacy disposition.** Comment identifies weak trait heritability and requests clearer output; no code fix.

**LDSC3 assessment — Methodological limitation.** LDSC3 retains own-trait and pairwise estimates with SE and failed-pair notes.
Eighteen-thousand samples do not guarantee adequate trait signal, and a GCTA estimate is not automatically an accuracy reference for a different SNP/reference/model population.
The reported p0.4065 correctly reflects weak evidence against zero.

**Recommended prevention.** Make denominator precision visible near correlation results and explain method/population differences before interpreting disagreement with GCTA.

**Evidence and verification.** source-inspected
[src/ldsc/_regression/estimators.py:RG](<../../../../src/ldsc/_regression/estimators.py>); [src/ldsc/_regression/runner.py:_fit_rg_dataset](<../../../../src/ldsc/_regression/runner.py>); [docs/reference/regression-options.md](<../../../reference/regression-options.md>).

<a id="issue-332"></a>

### #332: Interpretation of genetic correlation between a continuous trait and a binary trait

[Original issue](https://github.com/bulik/ldsc/issues/332) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **1**.

Workflow tags: `rg`.

**Trigger and harm.** User asks whether MMSE-versus-AD positive or negative correlation matches observed phenotype differences.
Effect coding and causal interpretation can be confused with genetic sharing.

**Legacy disposition.** Open; another user repeats the negative-correlation question.

**LDSC3 assessment — Methodological limitation.** LDSC3 aligns signed effects to alleles; a positive correlation means alleles increasing the coded quantitative trait tend to increase the coded binary liability, and negative means opposite directions.
Reversing case coding or the quantitative scale reverses the sign.
It need not equal phenotypic correlation or establish causation.

**Recommended prevention.** Document phenotype/case coding and signed-statistic orientation; interpret the sign only after checking both source GWAS definitions.

**Evidence and verification.** source-inspected
[src/ldsc/_regression/estimators.py:RG](<../../../../src/ldsc/_regression/estimators.py>); [src/ldsc/_regression/runner.py:_fit_rg_dataset](<../../../../src/ldsc/_regression/runner.py>); [docs/reference/regression-options.md](<../../../reference/regression-options.md>).

<a id="issue-334"></a>

### #334: BOLT-LMM v.2.3.5 in combination with LDSC

[Original issue](https://github.com/bulik/ldsc/issues/334) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **0**.

Workflow tags: `munge-sumstats`, `h2`.

**Trigger and harm.** Genotyped .stats retains159418 SNPs, while imputed .bgen.stats retains1.2million HM3 SNPs.
Selecting a sparse genotype-only file can give poor coverage and large uncertainty.

**Legacy disposition.** Open without a validated answer.

**LDSC3 assessment — Methodological limitation.** LDSC3 reads either table when mapped correctly, but HM3 membership does not itself ensure imputation quality.
The appropriate source should contain genome-wide well-calibrated association statistics, model-matched N, INFO/MAF QC and sufficient regression-SNP coverage; mixed-model statistic choice also requires producer guidance.

**Recommended prevention.** Explain the two BOLT output universes and separate HM3 restriction from INFO QC; use suitable imputed genome-wide results when scientifically justified.

**Evidence and verification.** source-inspected
[src/ldsc/_sumstats/munge_input.py:prepare_munge_input](<../../../../src/ldsc/_sumstats/munge_input.py>); [src/ldsc/_sumstats/munge_kernel.py:munge_sumstats](<../../../../src/ldsc/_sumstats/munge_kernel.py>).

<a id="issue-339"></a>

### #339: Variant filtering (not SNPs or strand-ambiguous)

[Original issue](https://github.com/bulik/ldsc/issues/339) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **0**.

Workflow tags: `munge-sumstats`, `h2`.

**Trigger and harm.** A rare non-imputed GWAS loses250117 rows as non-SNP or strand-ambiguous despite no displayed MAF removals.
Users confuse allele validity with rarity and may disable QC for a method outside its reference assumptions.

**Legacy disposition.** Open without a validated rare-variant workflow.

**LDSC3 assessment — Methodological limitation.** LDSC3 separates MAF filtering, identity validity, HM3 restriction and reference counts, but its default common-SNP denominator does not estimate arbitrary rare-variant heritability. A/T/C/G ambiguity is independent of rarity; changing --count-kind toall alone does not validate rare-variant LD/model calibration.

**Recommended prevention.** Explain each filter independently and state the validated common-SNP scope; require rare-variant-specific reference and method evidence before extending claims.

**Evidence and verification.** source-inspected
[docs/reference/ldscore-snp-universes.md](<../../../reference/ldscore-snp-universes.md>); [src/ldsc/_sumstats/munge_kernel.py:parse_dat](<../../../../src/ldsc/_sumstats/munge_kernel.py>); [docs/reference/regression-options.md](<../../../reference/regression-options.md>).

<a id="issue-362"></a>

### #362: h2 estimate and intercept

[Original issue](https://github.com/bulik/ldsc/issues/362) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **0**.

Workflow tags: `h2`.

**Trigger and harm.** An EAS meta-analysis gives much larger heritability when the intercept is forced to one.
Choosing the higher estimate as correct can attribute confounding or model mismatch to polygenic signal.

**Legacy disposition.** Open; no authoritative resolution in the accessible thread.

**LDSC3 assessment — Methodological limitation.** LDSC3 exposes intercept constraints but cannot infer that the scientific assumptions hold.
The ratio depends on mean chi-square minus one, so weak inflation can make it comparatively large without proving a particular cause.
Exact reproduction also depends on SNPs, counts, sample sizes and reference matching.

**Recommended prevention.** Explain free versus constrained intercept assumptions and a matched-input sensitivity comparison; do not select the constraint by the resulting heritability.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[src/ldsc/_regression/estimators.py:517](<../../../../src/ldsc/_regression/estimators.py#L517>); [docs/reference/regression-options.md](<../../../reference/regression-options.md>).

<a id="issue-363"></a>

### #363: Sample overlap

[Original issue](https://github.com/bulik/ldsc/issues/363) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **1**.

Workflow tags: `rg`.

**Trigger and harm.** Asks how to set genetic covariance with two overlapping people across binary and quantitative GWAS.
Treating overlap count alone as the cross-trait intercept can bias genetic covariance and genetic correlation.

**Legacy disposition.** Open; only a follow-up asking the same question.

**LDSC3 assessment — Methodological limitation.** The modeled intercept is on the z1*z2 scale and depends on shared-sample phenotypic covariance and sample-size scaling, not just a count of overlapping subjects.
LDSC3 permits estimating it instead of supplying an unsupported constant.

**Recommended prevention.** Give the intercept formula with its assumptions and recommend estimating the intercept when the shared-sample covariance is unknown.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[src/ldsc/_regression/estimators.py:683](<../../../../src/ldsc/_regression/estimators.py#L683>); [docs/reference/regression-options.md](<../../../reference/regression-options.md>).

<a id="issue-415"></a>

### #415: Sample size with repeated observations

[Original issue](https://github.com/bulik/ldsc/issues/415) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **1**.

Workflow tags: `munge-sumstats`, `h2`, `rg`.

**Trigger and harm.** A meta-analysis repeatedly observes the same individuals at different ages and asks whether N means people or observations.
Either naive count can mis-scale heritability if the association-statistic variance reflects correlated observations.

**Legacy disposition.** A commenter recommends age-specific GWAS; no derivation establishes an effective N for the posted meta-analysis.

**LDSC3 assessment — Methodological limitation.** LDSC3 takes N as supplied and cannot infer the correlation structure or recover an effective sample size from observation counts alone.
The GWAS/meta-analysis model must establish its information scaling before LDSC is applied.

**Recommended prevention.** Document that N must match the variance of the GWAS statistic; request the meta-analysis covariance model rather than choosing people or rows by rule of thumb.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[src/ldsc/_sumstats/munge_kernel.py:611](<../../../../src/ldsc/_sumstats/munge_kernel.py#L611>); [src/ldsc/_regression/estimators.py:238](<../../../../src/ldsc/_regression/estimators.py#L238>).

<a id="issue-419"></a>

### #419: rg with sample overlap

[Original issue](https://github.com/bulik/ldsc/issues/419) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **1**.

Workflow tags: `rg`.

**Trigger and harm.** Both GWAS derive from FinnGen and the user asks for extra restrictions to correct overlap.
Forcing zero covariance intercept can misattribute shared-sample covariance to genetics.

**Legacy disposition.** A commenter explains that an estimated cross-trait intercept absorbs overlap/phenotypic covariance under the model; no code fix is required.

**LDSC3 assessment — Methodological limitation.** LDSC3 can estimate the genetic-covariance intercept.
Cohort identity alone does not quantify the overlap covariance or guarantee every LDSC assumption; no blanket additional restriction follows from this thread.

**Recommended prevention.** Use the free intercept by default and explain the assumptions needed before fixing it; distinguish sample overlap from genetic correlation.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[src/ldsc/_regression/estimators.py:683](<../../../../src/ldsc/_regression/estimators.py#L683>); [docs/reference/regression-options.md](<../../../reference/regression-options.md>).

<a id="issue-431"></a>

### #431: INTERCEPT FLAG

[Original issue](https://github.com/bulik/ldsc/issues/431) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **1**.

Workflow tags: `rg`.

**Trigger and harm.** The user asks whether to constrain the intercept because doing so makes a result significant.
Post hoc constraint selection changes the estimator and invalidates the interpretation of the nominal test.

**Legacy disposition.** A commenter advises against those flags but supplies no case-specific derivation; no code fix is needed.

**LDSC3 assessment — Methodological limitation.** LDSC3 exposes the same scientifically meaningful choice.
Its presence in the CLI does not justify selecting it by the output p-value or shared ancestry alone.

**Recommended prevention.** Explain when an intercept can be treated as known, use a free intercept by default, and prespecify any constrained sensitivity analysis.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[docs/reference/regression-options.md](<../../../reference/regression-options.md>); [src/ldsc/_regression/estimators.py:238](<../../../../src/ldsc/_regression/estimators.py#L238>).

<a id="issue-436"></a>

### #436: Filtration of partition heritability results

[Original issue](https://github.com/bulik/ldsc/issues/436) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **2**.

Workflow tags: `ldsc-seg`.

**Trigger and harm.** Tests one small ATAC annotation across 106 GWAS and asks whether enrichment or coefficient p-values should enter FDR.
Selecting the statistic that yields discoveries or choosing by table orientation can answer the wrong scientific question.

**Legacy disposition.** A commenter recommends enrichment p based on testing traits rather than annotations; that rationale does not identify the estimand.
No code fix is involved.

**LDSC3 assessment — Methodological limitation.** LDSC3 documents separate enrichment and coefficient tests.
Marginal concentration versus additional signal beyond baseline remains the controlling distinction for every trait; the tiny annotation and large SE also limit precision.

**Recommended prevention.** State the null and testing family before FDR, use the same intended test across traits, and report support/SE rather than switch tests to obtain significance.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[docs/reference/results/partitioned-h2.md:49](<../../../reference/results/partitioned-h2.md#partitioned_h2tsv>); [docs/reference/results/partitioned-h2.md:76](<../../../reference/results/partitioned-h2.md#choosing-the-right-column>); [docs/user-guide/concepts/partitioned-h2-math.md](<../../../user-guide/concepts/partitioned-h2-math.md>).

<a id="issue-444"></a>

### #444: Estimating bivariate LDSC intercept when h2 is low

[Original issue](https://github.com/bulik/ldsc/issues/444) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **3**.

Workflow tags: `h2`, `rg`.

**Trigger and harm.** Asks whether a bivariate intercept remains useful when h2 is near/under zero and genetic correlation cannot be estimated.
Conflating intercept, overlap count and ratio failure can either discard useful covariance diagnostics or overclaim overlap identification.

**Legacy disposition.** Commenter notes no obvious derivation-based bias but admits missing simulation evidence; no universal guarantee or code fix is established.

**LDSC3 assessment — Methodological limitation.** LDSC3 can fit covariance components separately from an unestimable rg denominator, but the RG wrapper has delete-value boundary limitations.
An intercept can reflect several non-genetic sources and is not a direct overlap estimator.

**Recommended prevention.** Report covariance/intercept and ratio estimability separately, preserve uncertainty, and require regime-specific validation before using the intercept in downstream overlap models.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[src/ldsc/_regression/estimators.py:626](<../../../../src/ldsc/_regression/estimators.py#L626>); [src/ldsc/_regression/estimators.py:746](<../../../../src/ldsc/_regression/estimators.py#L746>); [docs/user-guide/workflows/genetic-correlation.md](<../../../user-guide/workflows/genetic-correlation.md>).

<a id="issue-452"></a>

### #452: Is enrichment same when using LD scores with different combination of categories?

[Original issue](https://github.com/bulik/ldsc/issues/452) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **3**.

Workflow tags: `ldscore`, `baseline-partitioned-h2`, `ldsc-seg`.

**Trigger and harm.** Asks whether unchanged LD-score columns imply unchanged enrichments, how to compare nested annotations, and whether HM3 rows represent all query SNPs.
Users can treat multivariable results as column-local or shrink the contributor set to the regression rows.

**Legacy disposition.** Commenter explains model dependence; later nested-set and HM3 questions remain unanswered.

**LDSC3 assessment — Methodological limitation.** LDSC3 fits annotation coefficients jointly, so changing covariates changes effects and overlap-derived summaries.
Reference SNPs define annotation coverage while regression SNPs are measurement rows. A subset-versus-parent contrast needs their covariance, not two independent ratios.

**Recommended prevention.** Explain the three distinctions with a nested-set example, retain full baseline adjustment appropriate to the aim, and define any subset contrast before implementing it.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[docs/reference/ldscore-snp-universes.md](<../../../reference/ldscore-snp-universes.md>); [src/ldsc/_regression/overlap_model.py](<../../../../src/ldsc/_regression/overlap_model.py>); [docs/reference/results/partitioned-h2.md](<../../../reference/results/partitioned-h2.md>).

<a id="issue-458"></a>

### #458: Genetic correlation &gt; 1

[Original issue](https://github.com/bulik/ldsc/issues/458) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **0**.

Workflow tags: `rg`, `h2`.

**Trigger and harm.** Two Pan-UKB phenotypes give rg slightly above one, large SE and very weak observed-scale h2.
Users may interpret overlap or liability conversion as a repair for noisy denominator estimates.

**Legacy disposition.** Open; no answer or implementation fix.

**LDSC3 assessment — Methodological limitation.** LDSC3 allows unconstrained estimates and reports scale-invariant rg.
Liability conversion multiplies covariance and heritabilities consistently, so it cannot improve the observed-scale signal-to-noise or force rg into bounds.

**Recommended prevention.** Report uncertainty and weak-signal diagnostics, explain scale invariance and free-intercept overlap adjustment, and avoid clipping or declaring the near-one estimate a code bug.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[src/ldsc/_regression/estimators.py:746](<../../../../src/ldsc/_regression/estimators.py#L746>); [src/ldsc/_regression/runner.py:2235](<../../../../src/ldsc/_regression/runner.py#L2235>); [docs/user-guide/workflows/genetic-correlation.md](<../../../user-guide/workflows/genetic-correlation.md>).

<a id="issue-467"></a>

### #467: Does immunochip summary data with another gwas summary data can be used to LDSC genetic correlation?

[Original issue](https://github.com/bulik/ldsc/issues/467) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **0**.

Workflow tags: `rg`.

**Trigger and harm.** The title asks whether Immunochip summary statistics can be correlated with genome-wide GWAS; no body is provided.
Targeted SNP ascertainment can violate genome-wide regression assumptions even when the files merge.

**Legacy disposition.** Open; no case-specific evidence or resolution.

**LDSC3 assessment — Methodological limitation.** LDSC3 does not validate a targeted array as a representative GWAS solely from its format.
The empty body prevents assessing backbone coverage, imputation or available statistics.

**Recommended prevention.** Explain the targeted-array limitation and require an appropriate genome-wide backbone/validated method; mark this individual application unverified.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[docs/reference/ldscore-snp-universes.md](<../../../reference/ldscore-snp-universes.md>); [Source](<https://github.com/bulik/ldsc/wiki/FAQ>).

<a id="issue-472"></a>

### #472: rg values not in the [-1,1] bound&amp;NAN Values

[Original issue](https://github.com/bulik/ldsc/issues/472) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **2**.

Workflow tags: `h2`, `rg`.

**Trigger and harm.** Very low or strongly negative trait heritability produces NA genetic correlations; the user also asks about out-of-bounds values. A nominally large sample does not rescue a deflated statistic distribution or an invalid ratio denominator.

**Legacy disposition.** A commenter suggests proxy-case bias, but the user does not know whether proxies were used.
The cause is not verified.

**LDSC3 assessment — Methodological limitation.** LDSC3 marks nonpositive whole-data heritability as unestimable for rg.
Proxy phenotype scaling, GC, sample-size errors, model mismatch and signal weakness remain possibilities, not established diagnoses.

**Recommended prevention.** Explain denominator estimability and inspect the actual GWAS design before invoking a proxy correction; preserve negative estimates and uncertainty instead of clipping them.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[src/ldsc/_regression/estimators.py:765](<../../../../src/ldsc/_regression/estimators.py#L765>); [docs/user-guide/workflows/genetic-correlation.md](<../../../user-guide/workflows/genetic-correlation.md>); [Source](<https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1010208>).

<a id="issue-476"></a>

### #476: Baseline LD scores

[Original issue](https://github.com/bulik/ldsc/issues/476) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **2**.

Workflow tags: `baseline-partitioned-h2`, `ldsc-seg`.

**Trigger and harm.** Asks whether baseline and baselineLD models are interchangeable for cell-specific analysis.
The two tests and covariate sets can be conflated, changing conditional meaning and p-values.

**Legacy disposition.** A commenter quotes source recommendations separating baseline v1.2 coefficient tests from baselineLD v2.2 enrichment; a follow-up asks why those goals differ.

**LDSC3 assessment — Methodological limitation.** LDSC3 distinguishes baseline-only enrichment from baseline-plus-query coefficient tests, but does not choose the scientifically appropriate baseline release.
More columns are not automatically better, and continuous baseline annotations need appropriate interpretation.

**Recommended prevention.** Tie baseline recommendations to original method/version and the stated null, with a resource compatibility table; do not present one globally best baseline.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[docs/reference/results/partitioned-h2.md:11](<../../../reference/results/partitioned-h2.md#two-regimes>); [docs/user-guide/workflows/cell-type-enrichment.md](<../../../user-guide/workflows/cell-type-enrichment.md>); [docs/reference/results/quantile-h2.md](<../../../reference/results/quantile-h2.md>).

<a id="issue-487"></a>

### #487: The pre-computed European LD scores included not only hapmap3 SNPs

[Original issue](https://github.com/bulik/ldsc/issues/487) · Source state: **closed** · Severity: **high** · Accessible comments reviewed: **0**.

Workflow tags: `ldscore`.

**Trigger and harm.** Asks why precomputed LD-score rows exceed HM3 and why the EUR panel has489 rather than503 people.
Filtering contributors to output rows changes LD scores; treating differing panel counts as an unexplained typo obscures QC/provenance.

**Legacy disposition.** Closed without an accessible answer or verified historical sample-selection explanation.

**LDSC3 assessment — Methodological limitation.** LDSC3 deliberately separates reference SNPs from regression SNPs and can keep individuals explicitly.
It cannot reconstruct historical ancestry/relatedness exclusions from file counts alone.

**Recommended prevention.** Document panel sample IDs/QC and each SNP set independently; do not invent why14people were excluded or restrict LD contributors to HM3 merely to match output row counts.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[docs/reference/ldscore-snp-universes.md](<../../../reference/ldscore-snp-universes.md>); [src/ldsc/_refpanel/adapters.py](<../../../../src/ldsc/_refpanel/adapters.py>); [docs/reference/formats/legacy-ldscores.md](<../../../reference/formats/legacy-ldscores.md>).

<a id="issue-488"></a>

### #488: SNPs included in the pre-computed LD scores for European are much less than SNPs in 1000 genome phase 1 dataset

[Original issue](https://github.com/bulik/ldsc/issues/488) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **0**.

Workflow tags: `ldscore`.

**Trigger and harm.** Compares chr22 LD-score rows with all Phase1 variants and asks what QC was applied.
Users may infer that only persisted rows contributed LD or reproduce the wrong release/filtering policy.

**Legacy disposition.** Open with no verified historical QC recipe.

**LDSC3 assessment — Methodological limitation.** LDSC3 records the reference-SNP funnel and stores a separate regression-row set. A small output table can coexist with a much larger reference SNP set, and historical Phase1/Phase3 resources must not be conflated.

**Recommended prevention.** Provide a funnel showing source variants, retained reference SNPs and persisted regression SNPs, plus explicit panel release and QC; preserve uncertainty about historical unrecorded steps.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[docs/reference/ldscore-snp-universes.md](<../../../reference/ldscore-snp-universes.md>); [tests/ldscore/test_reference_snp_funnel.py](<../../../../tests/ldscore/test_reference_snp_funnel.py>); [src/ldsc/_ldscore/outputs.py](<../../../../src/ldsc/_ldscore/outputs.py>).

<a id="faq-10"></a>

### FAQ-10: Power and sample size

[FAQ source](https://github.com/bulik/ldsc/wiki/FAQ) · Source-order identifier: **FAQ-10** · Severity: **high**.

Workflow tags: `h2`, `rg`, `baseline-partitioned-h2`, `ldsc-seg`.

**Legacy disposition.** Historical FAQ guidance; assessed independently rather than treated as a code-fix record.

**LDSC3 assessment — Methodological limitation.** Precision depends on signal, sample size, retained regression SNPs and model complexity. A universal sample-size threshold cannot guarantee reliable estimation, and a successful software run does not establish useful power.

**Recommended prevention.** Show SNP retention, SNP-heritability uncertainty and analysis-specific precision; consider individual-level methods where available.

**Evidence and verification.** Source-inspected; relevant executed checks are identified in verification.md.
[src/ldsc/_regression/runner.py:summarize_total_h2](<../../../../src/ldsc/_regression/runner.py>); [docs/user-guide/workflows/heritability.md](<../../../user-guide/workflows/heritability.md>); [Source](<https://pmc.ncbi.nlm.nih.gov/articles/PMC4495769/>).

<a id="faq-12"></a>

### FAQ-12: Unbounded correlation estimates

[FAQ source](https://github.com/bulik/ldsc/wiki/FAQ) · Source-order identifier: **FAQ-12** · Severity: **high**.

Workflow tags: `rg`.

**Legacy disposition.** Historical FAQ guidance; assessed independently rather than treated as a code-fix record.

**LDSC3 assessment — Methodological limitation.** A finite-sample genetic-correlation estimate can exceed its population parameter bounds.
Near-zero heritability and unjustified intercept constraints require scrutiny; clipping to [-1,1] would hide uncertainty and change the estimator.

**Recommended prevention.** Keep estimates and uncertainty distinct from population bounds; flag unstable or undefined cases clearly.

**Evidence and verification.** Source-inspected; relevant executed checks are identified in verification.md.
[src/ldsc/_regression/estimators.py:746](<../../../../src/ldsc/_regression/estimators.py#L746>); [docs/user-guide/workflows/genetic-correlation.md](<../../../user-guide/workflows/genetic-correlation.md>); [Source](<https://www.nature.com/articles/ng.3406>).

<a id="faq-13"></a>

### FAQ-13: Imprecise genetic correlations

[FAQ source](https://github.com/bulik/ldsc/wiki/FAQ) · Source-order identifier: **FAQ-13** · Severity: **high**.

Workflow tags: `rg`.

**Legacy disposition.** Historical FAQ guidance; assessed independently rather than treated as a code-fix record.

**LDSC3 assessment — Methodological limitation.** Weak heritability, limited SNP overlap or unstable jackknife deletions can make correlation uncertainty large or undefined.
These are inferential limitations and may also expose inadequate diagnostic handling.

**Recommended prevention.** Report pair-specific retained SNPs, both heritabilities and valid uncertainty; never replace undefined SE with zero.

**Evidence and verification.** Source-inspected; relevant executed checks are identified in verification.md.
[src/ldsc/_regression/runner.py:_summarize_rg_pair](<../../../../src/ldsc/_regression/runner.py>); [src/ldsc/_regression/estimators.py:772](<../../../../src/ldsc/_regression/estimators.py#L772>).

<a id="faq-15"></a>

### FAQ-15: Negative total estimates

[FAQ source](https://github.com/bulik/ldsc/wiki/FAQ) · Source-order identifier: **FAQ-15** · Severity: **high**.

Workflow tags: `h2`.

**Legacy disposition.** Historical FAQ guidance; assessed independently rather than treated as a code-fix record.

**LDSC3 assessment — Methodological limitation.** The unconstrained weighted regression can yield negative total estimates.
Distinguish a noisy estimate near zero from systematic input/model problems; zero-total derived ratios currently also have a reproducible failure boundary in V3.

**Recommended prevention.** Preserve signed estimates and SEs; investigate reproducibly negative estimates and treat undefined derived quantities explicitly.

**Evidence and verification.** Source-inspected; relevant executed checks are identified in verification.md.
[src/ldsc/_regression/estimators.py:Hsq](<../../../../src/ldsc/_regression/estimators.py>); [verification.md#v3-zero-total-derived-proportion](<verification.md#v3-zero-total-derived-proportion>).

<a id="faq-17"></a>

### FAQ-17: Estimates above parameter bounds

[FAQ source](https://github.com/bulik/ldsc/wiki/FAQ) · Source-order identifier: **FAQ-17** · Severity: **high**.

Workflow tags: `h2`.

**Legacy disposition.** Historical FAQ guidance; assessed independently rather than treated as a code-fix record.

**LDSC3 assessment — Methodological limitation.** Unconstrained estimates can exceed physical parameter bounds, and case-control ascertainment further complicates observed-scale interpretation.
Check sample-size convention, denominator and prevalence before assigning a biological meaning.

**Recommended prevention.** Do not clip estimates as a repair; present their SE and investigate scaling, model fit and input provenance.

**Evidence and verification.** Source-inspected; relevant executed checks are identified in verification.md.
[src/ldsc/_regression/estimators.py:Hsq](<../../../../src/ldsc/_regression/estimators.py>); [src/ldsc/_regression/estimators.py:h2_obs_to_liab](<../../../../src/ldsc/_regression/estimators.py>); [docs/reference/regression-options.md](<../../../reference/regression-options.md>).

<a id="faq-18"></a>

### FAQ-18: Negative annotation contributions

[FAQ source](https://github.com/bulik/ldsc/wiki/FAQ) · Source-order identifier: **FAQ-18** · Severity: **high**.

Workflow tags: `baseline-partitioned-h2`, `ldsc-seg`, `continuous`.

**Legacy disposition.** Historical FAQ guidance; assessed independently rather than treated as a code-fix record.

**LDSC3 assessment — Methodological limitation.** A negative conditional coefficient can coexist with positive overlap-aware annotation heritability.
Sampling noise, correlated predictors and model misspecification must be distinguished from a claim of negative biological variance.

**Recommended prevention.** Identify the reported column and null hypothesis before interpreting sign; avoid equating M_C times coefficient with the full overlap-aware share.

**Evidence and verification.** Source-inspected; relevant executed checks are identified in verification.md.
[src/ldsc/_regression/overlap_model.py:overlap_aware_category_table](<../../../../src/ldsc/_regression/overlap_model.py>); [docs/reference/results/partitioned-h2.md:145](<../../../reference/results/partitioned-h2.md#notes-for-interpretation>); [Source](<https://pmc.ncbi.nlm.nih.gov/articles/PMC4626285/>).

<a id="faq-20"></a>

### FAQ-20: Intercept estimates below unity

[FAQ source](https://github.com/bulik/ldsc/wiki/FAQ) · Source-order identifier: **FAQ-20** · Severity: **high**.

Workflow tags: `h2`.

**Legacy disposition.** Historical FAQ guidance; assessed independently rather than treated as a code-fix record.

**LDSC3 assessment — Methodological limitation.** A below-unity estimate is not by itself proof of a defective implementation.
Summary-statistic correction, calibration and model assumptions matter, while mean chi-square and its relationship to LD provide additional diagnostic context.

**Recommended prevention.** Ask how statistics were corrected and calibrated; do not force an intercept of one solely to make the result look conventional.

**Evidence and verification.** Source-inspected; relevant executed checks are identified in verification.md.
[src/ldsc/_regression/estimators.py:Hsq._ratio](<../../../../src/ldsc/_regression/estimators.py>); [docs/reference/regression-options.md](<../../../reference/regression-options.md>); [Source](<https://pmc.ncbi.nlm.nih.gov/articles/PMC4495769/>).

<a id="faq-21"></a>

### FAQ-21: Choosing intercept constraints

[FAQ source](https://github.com/bulik/ldsc/wiki/FAQ) · Source-order identifier: **FAQ-21** · Severity: **high**.

Workflow tags: `h2`, `rg`.

**Legacy disposition.** Historical FAQ guidance; assessed independently rather than treated as a code-fix record.

**LDSC3 assessment — Methodological limitation.** Free and constrained intercept fits encode different assumptions.
For cross-trait fits, known overlap alone does not universally determine the intercept: overlap-sample phenotypic covariance and relevant confounding enter the model.

**Recommended prevention.** Keep estimation as the default; document the evidence supporting any fixed value and avoid treating a constraint as a generic precision repair.

**Evidence and verification.** Source-inspected; relevant executed checks are identified in verification.md.
[src/ldsc/_regression/runner.py:_select_intercept](<../../../../src/ldsc/_regression/runner.py>); [tests/regression/test_regression_fit_outcomes.py](<../../../../tests/regression/test_regression_fit_outcomes.py>); [Source](<https://www.nature.com/articles/ng.3406>).

<a id="faq-22"></a>

### FAQ-22: Small heritability estimates

[FAQ source](https://github.com/bulik/ldsc/wiki/FAQ) · Source-order identifier: **FAQ-22** · Severity: **high**.

Workflow tags: `h2`.

**Legacy disposition.** Historical FAQ guidance; assessed independently rather than treated as a code-fix record.

**LDSC3 assessment — Methodological limitation.** The estimate refers to a specified SNP/reference universe and population, not automatically family-study heritability.
Sample-size scaling, corrected statistics and reference mismatch can alter it; a low value does not identify which explanation applies.

**Recommended prevention.** Report the SNP denominator and GWAS analysis sample, then inspect data preparation and model fit before comparing studies.

**Evidence and verification.** Source-inspected; relevant executed checks are identified in verification.md.
[docs/reference/ldscore-snp-universes.md](<../../../reference/ldscore-snp-universes.md>); [src/ldsc/_regression/estimators.py:LD_Score_Regression.aggregate](<../../../../src/ldsc/_regression/estimators.py>); [Source](<https://pmc.ncbi.nlm.nih.gov/articles/PMC4495769/>).

<a id="faq-25"></a>

### FAQ-25: Targeted genotyping designs

[FAQ source](https://github.com/bulik/ldsc/wiki/FAQ) · Source-order identifier: **FAQ-25** · Severity: **high**.

Workflow tags: `h2`, `rg`, `baseline-partitioned-h2`.

**Legacy disposition.** Historical FAQ guidance; assessed independently rather than treated as a code-fix record.

**LDSC3 assessment — Methodological limitation.** A targeted marker design can violate the genome-wide representativeness assumed by ordinary LDSC. A refactor cannot repair ascertainment of variants or people merely by accepting the file.
The FAQ explicitly identifies this as future work.

**Recommended prevention.** State the limitation prominently and require a validated, scientifically appropriate design rather than recommend extrapolation from targeted hits.

**Evidence and verification.** Source-inspected; relevant executed checks are identified in verification.md.
[Source](<https://github.com/bulik/ldsc/wiki/FAQ#genetic-correlation-and-heritability>); [Source](<https://pmc.ncbi.nlm.nih.gov/articles/PMC4495769/>); [docs/reference/ldscore-snp-universes.md](<../../../reference/ldscore-snp-universes.md>).

<a id="issue-388"></a>

### #388: Question about partitioned LDSC

[Original issue](https://github.com/bulik/ldsc/issues/388) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **1**.

Workflow tags: `rg`, `baseline-partitioned-h2`.

**Trigger and harm.** Interprets a meta-analysis of annotation enrichments as genetic correlation after selecting SNPs inside each annotation.
Filtering to a functional subset does not generally estimate annotation-specific genetic correlation and changes LD tagging.

**Legacy disposition.** Open; a follow-up repeats the question without a supported implementation.

**LDSC3 assessment — Intentionally unsupported.** LDSC3 separates univariate partitioned heritability from genome-wide rg and does not provide a validated partitioned covariance/correlation estimator.
Listing annotated SNPs cannot supply that missing estimand.

**Recommended prevention.** Explain cross-trait meta-analysis versus cross-trait covariance and prohibit presenting annotation-filtered ordinary rg as partitioned genetic correlation.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[src/ldsc/_regression/runner.py:918](<../../../../src/ldsc/_regression/runner.py#L918>); [src/ldsc/_regression/runner.py:1233](<../../../../src/ldsc/_regression/runner.py#L1233>); [docs/user-guide/concepts/ldsc-vs-partitioned-h2.md](<../../../user-guide/concepts/ldsc-vs-partitioned-h2.md>).

<a id="issue-441"></a>

### #441: Individual variant heritabilities?

[Original issue](https://github.com/bulik/ldsc/issues/441) · Source state: **closed** · Severity: **high** · Accessible comments reviewed: **4**.

Workflow tags: `h2`.

**Trigger and harm.** Requests individual SNP heritability from a genome-wide LDSC result. A fitted annotation-based expected contribution or marginal association PVE can be mistaken for an identified causal per-variant heritability.

**Legacy disposition.** Comments supply b^2 Var(G)/Var(Y) and phenotype-variance requirements; the issue closes without a package feature.

**LDSC3 assessment — Intentionally unsupported.** LDSC3 fits aggregate/annotation models, not independent causal effects for every variant.
The marginal formula needs effect scaling and does not allocate LD-shared genetic variance uniquely.

**Recommended prevention.** Explain marginal PVE, model-implied per-SNP expectation and causal variance allocation as different quantities; avoid adding an apparently precise per-SNP heritability export.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[src/ldsc/_regression/estimators.py:397](<../../../../src/ldsc/_regression/estimators.py#L397>); [docs/user-guide/concepts/partitioned-h2-math.md](<../../../user-guide/concepts/partitioned-h2-math.md>).

<a id="issue-451"></a>

### #451: How to obtain per-SNP heritability by LDSC?

[Original issue](https://github.com/bulik/ldsc/issues/451) · Source state: **closed** · Severity: **high** · Accessible comments reviewed: **4**.

Workflow tags: `h2`.

**Trigger and harm.** Assumes total LDSC heritability was computed by summing already estimated individual SNP contributions.
Users may infer causal SNP effects that the model does not identify.

**Legacy disposition.** Comments supply marginal PVE under genotype/phenotype variance assumptions and explain standardized-effect requirements; no package feature is added.

**LDSC3 assessment — Intentionally unsupported.** LDSC3 estimates aggregate regression coefficients and annotation-based expectations. A squared marginal GWAS coefficient is affected by LD and sampling error and is not uniquely a SNP contribution to the total.

**Recommended prevention.** Explain aggregate estimation versus marginal PVE and model-based annotation expectations; do not manufacture a causal per-SNP table from the total.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[src/ldsc/_regression/estimators.py:397](<../../../../src/ldsc/_regression/estimators.py#L397>); [docs/user-guide/concepts/partitioned-h2-math.md](<../../../user-guide/concepts/partitioned-h2-math.md>).

<a id="issue-498"></a>

### #498: LDSC for survival model

[Original issue](https://github.com/bulik/ldsc/issues/498) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **0**.

Workflow tags: `munge-sumstats`, `h2`.

**Trigger and harm.** Asks whether CoxPH log-hazard-ratio summary statistics from regenie imply a meaningful LDSC heritability. A valid association Z statistic does not establish the ordinary quantitative/liability heritability scale or sample-size normalization under censoring.

**Legacy disposition.** Open with no derivation, validated recipe or code fix.

**LDSC3 assessment — Intentionally unsupported.** LDSC3 can parse numeric beta/P/N, but parsing is not methodological support for survival-trait heritability.
Censoring, event information and association-model scaling need a justified mapping to the LDSC moment equation.

**Recommended prevention.** Document this as unvalidated and seek a survival-specific derivation or simulation before advising a generic liability conversion or effective N formula.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[src/ldsc/_sumstats/munge_input.py](<../../../../src/ldsc/_sumstats/munge_input.py>); [src/ldsc/_regression/estimators.py:238](<../../../../src/ldsc/_regression/estimators.py#L238>); [docs/reference/regression-options.md](<../../../reference/regression-options.md>).

<a id="issue-109"></a>

### #109: how to understand the results of partitioning heritability

[Original issue](https://github.com/bulik/ldsc/issues/109) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **3**.

Workflow tags: `baseline-partitioned-h2`, `ldsc-seg`.

**Trigger and harm.** A purported base has prop_snps=1.64 rather than one, and enrichment/conditional tests disagree.
An incoherent denominator can be mistaken for an interpretation nuance.

**Legacy disposition.** Unanswered; comments suspect a missing all-ones base but do not establish the root cause.

**LDSC3 assessment — Indeterminate.** LDSC3 validates counts/overlaps and distinguishes the two tests, but a label named base does not prove an all-ones column. A genuine binary all-reference-SNP base must have share one; above-one values can indicate mismatched count universes or a continuous column.
The exact input is unavailable, so do not assign a definitive cause.

**Recommended prevention.** Check annotation values and common/all denominators before interpreting P values; report invalid base semantics with actionable count/overlap diagnostics.

**Evidence and verification.** source-inspected
[src/ldsc/_regression/overlap_model.py](<../../../../src/ldsc/_regression/overlap_model.py>); [src/ldsc/_ldscore/source.py](<../../../../src/ldsc/_ldscore/source.py>); [docs/reference/results/partitioned-h2.md](<../../../reference/results/partitioned-h2.md>).

<a id="issue-430"></a>

### #430: p values in version Version 1.0.1 are different from v 2.0.0

[Original issue](https://github.com/bulik/ldsc/issues/430) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **3**.

Workflow tags: `rg`.

**Trigger and harm.** Independent implementations give similar genetic correlation but different significance.
Users may choose a version or intercept setting to recover significance without identifying the changed uncertainty.

**Legacy disposition.** A commenter hypothesizes different jackknife SE and warns that a fixed intercept requires scientific knowledge; no actual cause is verified.

**LDSC3 assessment — Indeterminate.** LDSC3 computes p from estimate/SE, while independent distributions can change row sets, weights, block construction and code.
The report cannot label this a repaired legacy bug without matched inputs and provenance.

**Recommended prevention.** Compare source revisions, resources, retained SNPs, blocks, estimate and SE first; document numerical-equivalence limits rather than promising unchanged p-values.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[src/ldsc/_regression/estimators.py:65](<../../../../src/ldsc/_regression/estimators.py#L65>); [src/ldsc/_regression/estimators.py:779](<../../../../src/ldsc/_regression/estimators.py#L779>); [docs/dev/validation/index.md](<../index.md>).

<a id="issue-468"></a>

### #468: Unexpected genetic correlation (rg) estimate

[Original issue](https://github.com/bulik/ldsc/issues/468) · Source state: **closed** · Severity: **high** · Accessible comments reviewed: **1**.

Workflow tags: `munge-sumstats`, `rg`.

**Trigger and harm.** Two traits match public h2 summaries but their genetic correlation changes magnitude/sign from GWAS Atlas.
Univariate Z-squared statistics can look correct even when effect orientation is wrong for cross-trait Z products.

**Legacy disposition.** Closed without an accessible resolution; a later user repeats the question.

**LDSC3 assessment — Indeterminate.** LDSC3 allele alignment helps but depends on correct signed-source mapping.
Different phenotype definitions, GWAS releases, SNP sets or intercepts are additional possibilities; agreement of h2 is not a sufficient cross-trait validation.

**Recommended prevention.** Check source trait identity and allele/sign coding on shared variants first, then match resources/settings; do not tune until the expected correlation appears.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[src/ldsc/_regression/estimators.py:786](<../../../../src/ldsc/_regression/estimators.py#L786>); [src/ldsc/_sumstats/munge_kernel.py:744](<../../../../src/ldsc/_sumstats/munge_kernel.py#L744>); [docs/user-guide/workflows/genetic-correlation.md](<../../../user-guide/workflows/genetic-correlation.md>).

<a id="issue-100"></a>

### #100: P value for coefficient

[Original issue](https://github.com/bulik/ldsc/issues/100) · Source state: **closed** · Severity: **high** · Accessible comments reviewed: **3**.

Workflow tags: `baseline-partitioned-h2`, `ldsc-seg`.

**Trigger and harm.** User has coefficient and SE but no direct coefficient P output.
An enrichment P value or an unspecified tail can be substituted for the intended conditional test.

**Legacy disposition.** A user and maintainer confirm coefficient/SE can be tested as a normal Z; no specific one- versus two-sided hypothesis is stated in the issue.

**LDSC3 assessment — Apparently prevented.** LDSC3 writes coefficient_z and one-sided coefficient_p for a positive conditional contribution, separately from two-sided enrichment_p. A two-sided departure-from-zero question is different and must not be answered by simply relabelling the positive-tail column.
Negative coefficients require retaining their sign.

**Recommended prevention.** State the null and alternative next to result columns; give formulas for an intentional two-sided coefficient test without conflating it with enrichment.

**Evidence and verification.** source-inspected
[docs/reference/results/partitioned-h2.md:49-53](<../../../reference/results/partitioned-h2.md#partitioned_h2tsv>); [src/ldsc/_regression/overlap_model.py](<../../../../src/ldsc/_regression/overlap_model.py>).

<a id="issue-225"></a>

### #225: Prop._h2 always is 1.0 use my own ldscores based on my own annot

[Original issue](https://github.com/bulik/ldsc/issues/225) · Source state: **closed** · Severity: **high** · Accessible comments reviewed: **0**.

Workflow tags: `ldsc-seg`.

**Trigger and harm.** A sparse custom annotation yields Prop._h2=1 despite containing many zeros. A self-normalized single-component model is mistaken for genome-wide enrichment.

**Legacy disposition.** Closed without answer; refers back to issue138, with no attributable fix.

**LDSC3 assessment — Apparently prevented.** LDSC3 requires explicit baseline annotations when query annotations are supplied.
That prevents the most direct custom-query-only workflow, while baseline-only h2 can still summarize its only component as1 by definition.
Sparsity alone does not create the missing outside-annotation comparison.

**Recommended prevention.** Explain normalization and require an all-SNP baseline plus appropriate covariates for a query-enrichment question.

**Evidence and verification.** source-inspected
[docs/reference/results/partitioned-h2.md](<../../../reference/results/partitioned-h2.md>); [src/ldsc/_regression/runner.py:summarize_partitioned_h2](<../../../../src/ldsc/_regression/runner.py>).

<a id="issue-375"></a>

### #375: About A1 / A2 allele definition in wiki and `munge_sumstats.py`

[Original issue](https://github.com/bulik/ldsc/issues/375) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **6**.

Workflow tags: `munge-sumstats`, `rg`.

**Trigger and harm.** Legacy messages call A1 a reference allele while users interpret that as genome reference rather than signed-statistic allele.
Misoriented effects can reverse cross-trait signals; a blanket both-traits-flipped argument does not repair inconsistent per-SNP coding.

**Legacy disposition.** A commenter explains signed A1 orientation but also conflates genomic reference with coding conventions; no code fix is documented.

**LDSC3 assessment — Apparently prevented.** The LDSC3 glossary explicitly defines summary-statistic A1 as the signed-statistic allele and says it is not necessarily genome reference.
Allele-aware regression still depends on correct source mapping.

**Recommended prevention.** Use effect_allele -> A1 and other_allele -> A2 examples, with a sign-flip example and a warning against inferring coding from REF/ALT alone.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[docs/dev/DOCS_VOCABULARY.md](<../../DOCS_VOCABULARY.md>); [src/ldsc/_regression/estimators.py:786](<../../../../src/ldsc/_regression/estimators.py#L786>); [docs/reference/formats/legacy-sumstats.md](<../../../reference/formats/legacy-sumstats.md>).

<a id="issue-381"></a>

### #381: Help in replicating heritability estimates

[Original issue](https://github.com/bulik/ldsc/issues/381) · Source state: **closed** · Severity: **high** · Accessible comments reviewed: **1**.

Workflow tags: `h2`.

**Trigger and harm.** Replication of published insomnia heritability is roughly half the published number.
Comparing different scales can be mistaken for a numerical software discrepancy.

**Legacy disposition.** The author explicitly identifies missing prevalence conversion and closes the issue; this is an analysis correction, not a code fix.

**LDSC3 assessment — Apparently prevented.** LDSC3 provides labeled observed/liability fields and a saved-result scale-conversion command, but cannot choose the appropriate population prevalence.
The supplied paper inputs were not rerun.

**Recommended prevention.** Require matching scale, ascertainment and prevalence when reproducing a paper; show which output fields can be compared directly.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[src/ldsc/_regression/estimators.py:104](<../../../../src/ldsc/_regression/estimators.py#L104>); [docs/user-guide/workflows/heritability.md](<../../../user-guide/workflows/heritability.md>); [Source](<https://github.com/bulik/ldsc/issues/381#issuecomment-1537044841>).

<a id="issue-393"></a>

### #393: Why is the number of SNPs in a cell type the same to the control's one?

[Original issue](https://github.com/bulik/ldsc/issues/393) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **2**.

Workflow tags: `ldscore`, `ldsc-seg`.

**Trigger and harm.** Confuses equal LD-score row sets with equal numbers of SNPs inside query and control gene annotations.
Users may filter away all outside-query SNPs and alter the regression rather than computing annotation-weighted LD scores.

**Legacy disposition.** A commenter explains that every regression SNP receives a score even when it is outside the query; an ATAC control follow-up remains unanswered.

**LDSC3 assessment — Apparently prevented.** LDSC3 distinguishes reference SNPs, regression SNPs, annotation values and control gene lists.
Its joint computation can preserve aligned rows while query membership differs, but tutorial explanations must make that distinction visible.

**Recommended prevention.** Use a small matrix showing shared SNP rows and different query/control columns; explain the biological purpose of the control for expression versus ATAC annotations.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[docs/reference/ldscore-snp-universes.md](<../../../reference/ldscore-snp-universes.md>); [docs/user-guide/workflows/cell-type-enrichment.md](<../../../user-guide/workflows/cell-type-enrichment.md>); [src/ldsc/_ldscore/calculator.py](<../../../../src/ldsc/_ldscore/calculator.py>).

<a id="issue-439"></a>

### #439: Filtering of munge_sumstats.py

[Original issue](https://github.com/bulik/ldsc/issues/439) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **7**.

Workflow tags: `munge-sumstats`, `rg`.

**Trigger and harm.** A thread about strand-ambiguous SNP filtering turns into conflicting instructions that A1 must be genomic REF.
Users can flip correct effects incorrectly, especially when REF/ALT conventions differ between studies.

**Legacy disposition.** Several community replies contradict the signed-A1 documentation; no verified code fix resolves the conceptual disagreement.

**LDSC3 assessment — Apparently prevented.** LDSC3 explicitly defines A1 as the signed-statistic allele and distinguishes it from genome reference.
Squared-LD scores are invariant to allele reversal; signed cross-trait alignment still matters.
Ambiguous SNP exclusion is a harmonization safeguard, not proof all palindromes are biologically invalid.

**Recommended prevention.** Provide one consistent effect_allele mapping with examples of beta sign and OR inversion; warn readers that issue comments are not authoritative package contracts.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[docs/dev/DOCS_VOCABULARY.md](<../../DOCS_VOCABULARY.md>); [src/ldsc/_sumstats/munge_kernel.py:232](<../../../../src/ldsc/_sumstats/munge_kernel.py#L232>); [src/ldsc/_regression/estimators.py:786](<../../../../src/ldsc/_regression/estimators.py#L786>).

<a id="issue-447"></a>

### #447: Questions Regarding S-LDSC Analysis with --ref-ld-chr-cts

[Original issue](https://github.com/bulik/ldsc/issues/447) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **9**.

Workflow tags: `ldsc-seg`.

**Trigger and harm.** Asks why cell-type output emphasizes coefficients and whether all-gene controls filter SNPs, act as covariates or modify effects.
Users may confuse conditional query effects with raw enrichment, or mistake shared controls for tissue-independence guarantees.

**Legacy disposition.** Comments correctly identify joint covariates but loosely describe raw coefficients as direct annotation changes in chi-square; no code fix is provided.

**LDSC3 assessment — Apparently prevented.** LDSC3 supports control gene lists as baseline annotations and reports full per-query models.
Regression is on annotation-specific LD scores, with coefficients describing conditional per-SNP contributions; an all-gene control adjusts generic gene proximity.

**Recommended prevention.** Include a model equation with baseline, all-gene control and one query, show both hypotheses, and avoid causal or independence claims from conditional association.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[docs/user-guide/workflows/cell-type-enrichment.md](<../../../user-guide/workflows/cell-type-enrichment.md>); [docs/reference/results/partitioned-h2.md](<../../../reference/results/partitioned-h2.md>); [src/ldsc/_annotation/gene_lists.py](<../../../../src/ldsc/_annotation/gene_lists.py>).

<a id="issue-484"></a>

### #484: Clarification on Baseline Model and Obtaining Enrichment Scores for Cell-Type-Specific Heritability (e.g., with scATAC-seq annotations)

[Original issue](https://github.com/bulik/ldsc/issues/484) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **1**.

Workflow tags: `ldsc-seg`.

**Trigger and harm.** Asks whether an ATAC query needs a full baseline and how to obtain enrichment from coefficient-only cts output. A coefficient is not the annotation heritability itself, and base-only fits can confound query specificity with general regulatory annotations.

**Legacy disposition.** A commenter cites the original separate-query/full-baseline analysis and a legacy output function; no code fix is needed.

**LDSC3 assessment — Apparently prevented.** LDSC3 reports both conditional coefficient and overlap-aware enrichment in full per-query models, with distinct tests.
Baseline release and the biological goal still need prespecification, especially for continuous baselineLD annotations.

**Recommended prevention.** Lead the tutorial with the scientific question and matching p-value; retain full baseline adjustment appropriate to that analysis and explain the coefficient units explicitly.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[docs/reference/results/partitioned-h2.md](<../../../reference/results/partitioned-h2.md>); [docs/user-guide/workflows/cell-type-enrichment.md](<../../../user-guide/workflows/cell-type-enrichment.md>); [src/ldsc/_regression/runner.py:918](<../../../../src/ldsc/_regression/runner.py#L918>).

<a id="faq-16"></a>

### FAQ-16: Observed and liability scales

[FAQ source](https://github.com/bulik/ldsc/wiki/FAQ) · Source-order identifier: **FAQ-16** · Severity: **high**.

Workflow tags: `h2`, `rg`.

**Legacy disposition.** Historical FAQ guidance; assessed independently rather than treated as a code-fix record.

**LDSC3 assessment — Apparently prevented.** LDSC3 stores observed and liability values separately and uses trait-specific sample and population prevalences.
Genetic correlation is invariant to these positive rescalings; covariance and heritability are not.

**Recommended prevention.** Require appropriate prevalence inputs and compare the same estimand and scale across studies.

**Evidence and verification.** Source-inspected; relevant executed checks are identified in verification.md.
[src/ldsc/_regression/estimators.py:h2_obs_to_liab](<../../../../src/ldsc/_regression/estimators.py>); [src/ldsc/_regression/runner.py:_summarize_rg_pair](<../../../../src/ldsc/_regression/runner.py>); [tests/regression/test_prevalence.py](<../../../../tests/regression/test_prevalence.py>).

<a id="faq-19"></a>

### FAQ-19: Overlapping annotation interpretation

[FAQ source](https://github.com/bulik/ldsc/wiki/FAQ) · Source-order identifier: **FAQ-19** · Severity: **high**.

Workflow tags: `baseline-partitioned-h2`, `ldsc-seg`.

**Legacy disposition.** Historical FAQ guidance; assessed independently rather than treated as a code-fix record.

**LDSC3 assessment — Apparently prevented.** LDSC3 uses the overlap matrix automatically for partitioned summaries.
Annotation shares need not sum to one across overlapping sets, and the conditional coefficient answers a different question from marginal enrichment.

**Recommended prevention.** Explain both quantities with a small overlapping-set example and select the test by scientific question.

**Evidence and verification.** Source-inspected; relevant executed checks are identified in verification.md.
[src/ldsc/_regression/overlap_model.py](<../../../../src/ldsc/_regression/overlap_model.py>); [docs/reference/results/partitioned-h2.md:6](<../../../reference/results/partitioned-h2.md#partitioned-heritability-result-columns>); [tests/regression/test_overlap_matrix.py](<../../../../tests/regression/test_overlap_matrix.py>).

<a id="issue-3"></a>

### #3: CCA weights

[Original issue](https://github.com/bulik/ldsc/issues/3) · Source state: **closed** · Severity: **moderate** · Accessible comments reviewed: **1**.

Workflow tags: `h2`.

**Trigger and harm.** An observed-scale estimate falls outside [0,1] in an ascertained case-control analysis.
Constraining the reported estimate and stabilizing an iterative weight formula answer different questions.

**Legacy disposition.** The discussion raises a Taylor-approximation/ascertainment question; closure supplies no derivation or verified fix.

**LDSC3 assessment — Methodological limitation.** Hsq.weights clips the provisional heritability used for weights; the reported fitted total is not forced into [0,1].
This preserves a numerical stabilization inherited from legacy LDSC, not a proof that every ascertained phenotype satisfies the working variance model.
Liability conversion separately requires sample and population prevalences.

**Recommended prevention.** Explain the distinction between weight stabilization, an unconstrained estimate, and observed-to-liability conversion; do not present clipping as a statistical correction.

**Evidence and verification.** source-inspected
[src/ldsc/_regression/estimators.py:Hsq.weights](<../../../../src/ldsc/_regression/estimators.py>); [src/ldsc/_regression/prevalence.py](<../../../../src/ldsc/_regression/prevalence.py>); [docs/user-guide/workflows/liability-scale.md](<../../../user-guide/workflows/liability-scale.md>).

<a id="issue-50"></a>

### #50: Question regarding LD Score Estimation

[Original issue](https://github.com/bulik/ldsc/issues/50) · Source state: **closed** · Severity: **moderate** · Accessible comments reviewed: **4**.

Workflow tags: `munge-sumstats`, `ldscore`, `rg`.

**Trigger and harm.** Genome-wide hits are absent from HM3 or precomputed regression LD-score rows.
Users may append selected significant variants and thereby bias marker selection.

**Legacy disposition.** Maintainer explains representative polygenic tagging, separates regression/reference SNP filters and cites Price et al. long-range LD regions.

**LDSC3 assessment — Methodological limitation.** LDSC3 separately exposes summary-statistics restrictions, regression-SNP selection and reference-panel preparation.
The packaged default need not contain each lead SNP to tag polygenic signal; adding only hits is not an unbiased extension.
Changing restriction files also requires corresponding LD scores and weights, not merely extra raw GWAS rows.

**Recommended prevention.** Give one end-to-end custom regression-SNP example and explain that long-range LD exclusions must match genome build and reference assumptions.

**Evidence and verification.** source-inspected
[docs/reference/ldscore-snp-universes.md](<../../../reference/ldscore-snp-universes.md>); [src/ldsc/_ldscore/calculator.py:_snp_universe_policy](<../../../../src/ldsc/_ldscore/calculator.py>); [docs/reference/formats/region-exclusion-presets.md](<../../../reference/formats/region-exclusion-presets.md>).

<a id="issue-56"></a>

### #56: Understanding results from LD regression

[Original issue](https://github.com/bulik/ldsc/issues/56) · Source state: **closed** · Severity: **moderate** · Accessible comments reviewed: **1**.

Workflow tags: `rg`.

**Trigger and harm.** A correlation estimate .91 with SE 2.64 is mistaken for the Z statistic .3454 in a wrapped table.
Users may report a large point estimate as strong evidence despite extreme uncertainty.

**Legacy disposition.** Maintainer links an answer on the user group; the issue itself does not contain its text.

**LDSC3 assessment — Methodological limitation.** LDSC3 writes named TSV columns, preventing terminal wrapping from defining column positions.
The quantities remain distinct: genetic correlation estimates effect similarity, SE measures uncertainty, and a zero-test P value is based on the corresponding statistic. A high point estimate with this SE is not evidence of near-identical genetic architecture.

**Recommended prevention.** Add one annotated example with weak evidence and require readers to use column names rather than visual position in a wrapped log.

**Evidence and verification.** source-inspected
[src/ldsc/_regression/outputs.py:RG_CONCISE_COLUMNS](<../../../../src/ldsc/_regression/outputs.py>); [docs/user-guide/workflows/genetic-correlation.md](<../../../user-guide/workflows/genetic-correlation.md>).

<a id="issue-67"></a>

### #67: Stratified h2

[Original issue](https://github.com/bulik/ldsc/issues/67) · Source state: **closed** · Severity: **moderate** · Accessible comments reviewed: **2**.

Workflow tags: `baseline-partitioned-h2`, `ldsc-seg`.

**Trigger and harm.** User expects every cell annotation to have proportion of heritability one.
Confusing normalization with an annotation-specific result obscures the tested scientific question.

**Legacy disposition.** Maintainer explains h2/SE and that only the all-SNP base has proportion one.

**LDSC3 assessment — Methodological limitation.** LDSC3 exposes total SNP heritability with SE and overlap-aware annotation shares.
Query rows are separate baseline-plus-query models, whereas baseline-only rows share one model.
Proportions need not sum to one for overlapping annotations, and coefficient significance is not derived from total-h2 Z.

**Recommended prevention.** Provide a small base-plus-overlapping-annotations example and distinguish total-h2 signal screening from query significance.

**Evidence and verification.** source-inspected
[docs/reference/results/partitioned-h2.md:11-30](<../../../reference/results/partitioned-h2.md#two-regimes>); [docs/reference/results/partitioned-h2.md:41-61](<../../../reference/results/partitioned-h2.md#partitioned_h2tsv>).

<a id="issue-74"></a>

### #74: Different h2 from the same sumstats.gz file using --h2 flag and --rg flag

[Original issue](https://github.com/bulik/ldsc/issues/74) · Source state: **closed** · Severity: **moderate** · Accessible comments reviewed: **3**.

Workflow tags: `h2`, `rg`.

**Trigger and harm.** Correlation uses SNPs shared by both traits and valid aligned alleles.
Users may treat a pair-specific estimate as interchangeable with a standalone genome-wide estimate.

**Legacy disposition.** Maintainer confirms differing SNP intersections explain the small sign change near zero.

**LDSC3 assessment — Methodological limitation.** LDSC3 retains separate h2_per_trait and pair-specific full outputs.
It estimates covariance and its denominator h2 on the same final pair intersection, with explicit pair filtering.
This is required for the ratio and is not solved by substituting the larger standalone h2 value when a pair becomes unstable.

**Recommended prevention.** Label pairwise SNP counts and h2 prominently and explain which estimate answers the standalone versus correlation question.

**Evidence and verification.** source-inspected
[src/ldsc/_regression/runner.py:build_rg_dataset](<../../../../src/ldsc/_regression/runner.py>); [src/ldsc/_regression/runner.py:estimate_rg_pairs](<../../../../src/ldsc/_regression/runner.py>); [tests/regression/test_regression_workflow.py:test_build_rg_dataset_uses_final_three_way_snp_intersection](<../../../../tests/regression/test_regression_workflow.py>).

<a id="issue-130"></a>

### #130: Genetic correlation on same data is not 1

[Original issue](https://github.com/bulik/ldsc/issues/130) · Source state: **closed** · Severity: **moderate** · Accessible comments reviewed: **2**.

Workflow tags: `rg`.

**Trigger and harm.** Identical summaries give rg .999 and covariance intercept .8919. A useful identity check is mistaken for a requirement that the covariance intercept always equals one.

**Legacy disposition.** Maintainer explains the expected self covariance intercept equals the fitted univariate intercept and acknowledges small estimator-implementation differences.

**LDSC3 assessment — Methodological limitation.** LDSC3 fits h2 and covariance separately with potentially differing weighting/start details.
Roundoff-level self-correlation departures are not automatically defects; a materially larger discrepancy needs matched-filter and weighting inspection. A covariance intercept below one can be consistent with GC or input calibration.

**Recommended prevention.** Document a tolerance-based self-correlation sanity check and compare intercepts on the correct scale, without forcing exact unity or clipping results.

**Evidence and verification.** source-inspected
[src/ldsc/_regression/estimators.py:RG](<../../../../src/ldsc/_regression/estimators.py>); [src/ldsc/_regression/runner.py:_fit_rg_dataset](<../../../../src/ldsc/_regression/runner.py>); [tests/regression/test_kernel_regression.py](<../../../../tests/regression/test_kernel_regression.py>).

<a id="issue-141"></a>

### #141: Questions about genesets' sample size when doing the regression

[Original issue](https://github.com/bulik/ldsc/issues/141) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **1**.

Workflow tags: `ldsc-seg`.

**Trigger and harm.** User wants query lists of 100-200 genes rather than 2,000 genes each. A nonsignificant query may be blamed on a fixed size rule, while sparse genomic support can make uncertainty unestimable.

**Legacy disposition.** Unanswered; no evidence establishes a mandatory 2,000-gene threshold.

**LDSC3 assessment — Methodological limitation.** LDSC3 accepts varying gene-list sizes, records resolution/support, and diagnoses collinearity and singular delete-block fits.
Power depends on SNP coverage, LD, effect architecture, baseline overlap, GWAS signal and block distribution, not gene count alone.
Query coefficient P values do not directly compare lists of different biological selection rules.

**Recommended prevention.** Document support/precision checks and a prespecified sensitivity or simulation design; avoid a hard minimum gene count without evidence.

**Evidence and verification.** source-inspected
[src/ldsc/_annotation/query_outcomes.py](<../../../../src/ldsc/_annotation/query_outcomes.py>); [src/ldsc/_regression/jackknife.py](<../../../../src/ldsc/_regression/jackknife.py>); [docs/user-guide/workflows/gene-list-ld-scores.md](<../../../user-guide/workflows/gene-list-ld-scores.md>).

<a id="issue-153"></a>

### #153: Issue related to method derivation

[Original issue](https://github.com/bulik/ldsc/issues/153) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **0**.

Workflow tags: `h2`.

**Trigger and harm.** Reader treats a standardized-effect derivation as an identity for every raw GWAS beta.
Using N*beta² without matching scaling can produce incorrect association statistics.

**Legacy disposition.** Open with no answer; the question cites the 2015 supplement.

**LDSC3 assessment — Methodological limitation.** LDSC3 munging converts valid one-df P values to chi-square/Z, rather than blindly applying N times any input beta squared.
The derivation assumes standardized genotype/phenotype and an asymptotic sampling variance; raw beta/SE or other association models need their own statistic definition.
This is a teaching gap rather than a demonstrated implementation flaw.

**Recommended prevention.** Add a short assumptions-first derivation connecting Z=beta/SE to sqrt(N)*standardized beta, and state the approximation limits.

**Evidence and verification.** source-inspected
[src/ldsc/_sumstats/munge_kernel.py:p_to_z](<../../../../src/ldsc/_sumstats/munge_kernel.py>); [docs/user-guide/concepts/ldsc-vs-partitioned-h2.md](<../../../user-guide/concepts/ldsc-vs-partitioned-h2.md>).

<a id="issue-154"></a>

### #154: minimum sequence or number of SNPs needed for estimating partitioned heritability in my own annotations?

[Original issue](https://github.com/bulik/ldsc/issues/154) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **0**.

Workflow tags: `ldsc-seg`.

**Trigger and harm.** Custom intervals cover under 1Mb and show no enrichment.
Nonsignificance may be interpreted as no biology or a universal minimum-span violation.

**Legacy disposition.** Unanswered; no data support a hard size threshold.

**LDSC3 assessment — Methodological limitation.** LDSC3 accepts sparse annotations but can only estimate coefficients when the LD-score design and delete-block fits retain support.
Power depends on common reference SNPs, LD spillover, baseline overlap, GWAS signal and spatial distribution, not base-pair length alone.
Query coefficient and marginal enrichment answer different power questions.

**Recommended prevention.** Report support and coefficient uncertainty, avoid arbitrary minimum Mb/gene thresholds, and use prespecified matched simulations if annotation-specific power is scientifically needed.

**Evidence and verification.** source-inspected
[src/ldsc/_annotation/query_outcomes.py](<../../../../src/ldsc/_annotation/query_outcomes.py>); [src/ldsc/_regression/jackknife.py](<../../../../src/ldsc/_regression/jackknife.py>); [docs/reference/results/partitioned-h2.md](<../../../reference/results/partitioned-h2.md>).

<a id="issue-193"></a>

### #193: Meaning of p-value in genetic correlation output (rg)

[Original issue](https://github.com/bulik/ldsc/issues/193) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **0**.

Workflow tags: `rg`.

**Trigger and harm.** The output says only p-value for rg without spelling out the null.
Users may confuse evidence against zero with evidence for equality or a well-calibrated estimate.

**Legacy disposition.** Open without clarification in the thread.

**LDSC3 assessment — Methodological limitation.** LDSC3 computes a two-sided normal test using genetic correlation divided by its block-jackknife standard error.
This tests a zero genetic correlation, conditional on valid model assumptions, and does not test equality to1, equivalence of two studies or absence of confounding.

**Recommended prevention.** State H0:r_g=0 and two-sided alternative beside the p column; provide separate guidance for comparison with1.

**Evidence and verification.** source-inspected
[src/ldsc/_regression/estimators.py:p_z_norm](<../../../../src/ldsc/_regression/estimators.py>); [src/ldsc/_regression/estimators.py:RG](<../../../../src/ldsc/_regression/estimators.py>); [docs/user-guide/workflows/genetic-correlation.md](<../../../user-guide/workflows/genetic-correlation.md>).

<a id="issue-236"></a>

### #236: WARNING: median value of SIGNED_SUMSTATS is 0.5 (should be close to 0.0). This column may be mislabeled

[Original issue](https://github.com/bulik/ldsc/issues/236) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **2**.

Workflow tags: `munge-sumstats`.

**Trigger and harm.** A direction column has median0.5 although the expected beta/Z null is 0.
Wrongly mapped or uniformly positive effects can corrupt direction, especially genetic correlation.

**Legacy disposition.** Open with repeated users but no validated explanation or fix.

**LDSC3 assessment — Methodological limitation.** LDSC3 retains a median-null sanity check and actionable --signed-sumstats guidance.
The check is a heuristic for coding/orientation mistakes, not a theorem that every selected or ascertained dataset has median0.
OR has a null of1 and needs different treatment from beta/Z.

**Recommended prevention.** Explain the null-value argument and common wrong columns; do not encourage centering effects solely to pass the check.

**Evidence and verification.** source-inspected
[src/ldsc/_sumstats/munge_input.py:prepare_munge_input](<../../../../src/ldsc/_sumstats/munge_input.py>); [src/ldsc/_sumstats/munge_kernel.py:munge_sumstats](<../../../../src/ldsc/_sumstats/munge_kernel.py>).

<a id="issue-262"></a>

### #262: RG&gt;1 &amp;&amp; sumsates.gz without column P

[Original issue](https://github.com/bulik/ldsc/issues/262) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **0**.

Workflow tags: `munge-sumstats`, `rg`.

**Trigger and harm.** User suspects missing output P caused genetic correlation above1.
Normal transformed output is mistaken for data loss or a cause of an unbounded estimate.

**Legacy disposition.** Open without an answer or fix.

**LDSC3 assessment — Methodological limitation.** LDSC3 raw munging uses P to derive signed Z and writes an analysis-ready Z,N table; direct legacy regression does not require P.
Its absence is therefore intentional, whereas an out-of-bounds correlation requires checking precision, allele orientation and model assumptions independently.

**Recommended prevention.** Explain the raw-P to signed-Z transformation and distinguish expected output schema from invalid correlation diagnostics.

**Evidence and verification.** source-inspected
[src/ldsc/_sumstats/munge_input.py:prepare_munge_input](<../../../../src/ldsc/_sumstats/munge_input.py>); [src/ldsc/_sumstats/munge_kernel.py:munge_sumstats](<../../../../src/ldsc/_sumstats/munge_kernel.py>).

<a id="issue-284"></a>

### #284: genetic enrichment of 3'UTR

[Original issue](https://github.com/bulik/ldsc/issues/284) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **0**.

Workflow tags: `ldsc-seg`.

**Trigger and harm.** A user asks whether variants are enriched in selected3primeUTRs.
LDSC heritability enrichment can be confused with overrepresentation of associated SNP counts.

**Legacy disposition.** Open without clarification.

**LDSC3 assessment — Methodological limitation.** LDSC3 can test whether SNPs annotated by these UTRs contribute conditional SNP heritability beyond the baseline, with LD adjustment.
That differs from asking whether more GWAS hits fall in a region or whether an individual locus is causal. A matched binary query, broad baseline and sufficient support are needed.

**Recommended prevention.** Start the workflow with the scientific question and distinguish coefficient/enrichment tests from SNP-count enrichment.

**Evidence and verification.** source-inspected
[docs/reference/results/partitioned-h2.md](<../../../reference/results/partitioned-h2.md>); [src/ldsc/_regression/runner.py:summarize_partitioned_h2](<../../../../src/ldsc/_regression/runner.py>).

<a id="issue-410"></a>

### #410: correlation value higher than 1

[Original issue](https://github.com/bulik/ldsc/issues/410) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **1**.

Workflow tags: `rg`.

**Trigger and harm.** An estimated genetic correlation of 1.0316 has SE 1.1478.
Users may treat an unconstrained noisy estimate as proof of impossible biology or a code bug.

**Legacy disposition.** A commenter explains unconstrained estimation and notes consistency with zero; no code change is proposed.

**LDSC3 assessment — Methodological limitation.** LDSC3 retains an unconstrained ratio estimate and its uncertainty.
The observed value alone is not evidence of an implementation defect, and clipping to one would hide instability.

**Recommended prevention.** Show interval/precision guidance and denominator checks alongside the point estimate; explain that a nonzero-null conclusion is unsupported here.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[src/ldsc/_regression/estimators.py:746](<../../../../src/ldsc/_regression/estimators.py#L746>); [docs/user-guide/workflows/genetic-correlation.md](<../../../user-guide/workflows/genetic-correlation.md>).

<a id="issue-424"></a>

### #424: Missing frequency from GWAS

[Original issue](https://github.com/bulik/ldsc/issues/424) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **1**.

Workflow tags: `munge-sumstats`, `baseline-partitioned-h2`.

**Trigger and harm.** Asks whether missing allele frequency in some GWAS changes partitioned results.
Users can conflate optional raw-data QC with the reference-panel MAF that defines counts and LD scores.

**Legacy disposition.** A commenter states regression uses reference MAF, but that does not imply raw GWAS frequency can never affect munging filters.

**LDSC3 assessment — Methodological limitation.** LDSC3 uses reference-derived MAF for LD scores and common-SNP counts.
Its munger can also filter using an optional raw FRQ column, so inclusion can affect which rows survive if QC thresholds are active.

**Recommended prevention.** Explain both roles separately and record whether optional GWAS frequency QC was applied when comparing the ten traits.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[src/ldsc/_sumstats/munge_kernel.py:216](<../../../../src/ldsc/_sumstats/munge_kernel.py#L216>); [docs/reference/ldscore-snp-universes.md](<../../../reference/ldscore-snp-universes.md>); [src/ldsc/_annotation/parsing.py](<../../../../src/ldsc/_annotation/parsing.py>).

<a id="issue-427"></a>

### #427: Causes for Ratio &lt; 0 other than GC?

[Original issue](https://github.com/bulik/ldsc/issues/427) · Source state: **closed** · Severity: **moderate** · Accessible comments reviewed: **2**.

Workflow tags: `h2`, `rg`.

**Trigger and harm.** An intercept slightly below one yields a negative ratio despite no known genomic-control correction.
Users may confuse genomic control with genetic engineering or misdiagnose small random deviations as a package defect.

**Legacy disposition.** A commenter explains the intercept/SE and flags low heritability as the more serious issue; author acknowledges.

**LDSC3 assessment — Methodological limitation.** LDSC3 computes the ratio from intercept minus one over mean chi-square minus one. A negative value alone does not identify its cause and low signal can destabilize genetic correlation.

**Recommended prevention.** Describe the formula and uncertainty, use the full term genomic control, and avoid presenting ratio<0 as a specific diagnosis.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[src/ldsc/_regression/estimators.py:517](<../../../../src/ldsc/_regression/estimators.py#L517>); [docs/user-guide/workflows/heritability.md](<../../../user-guide/workflows/heritability.md>).

<a id="issue-437"></a>

### #437: Partition heritability for UKB data

[Original issue](https://github.com/bulik/ldsc/issues/437) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **0**.

Workflow tags: `ldsc-seg`.

**Trigger and harm.** Asks whether all UKB summary statistics need special partitioned-analysis parameters. A cohort brand can be mistaken for a complete specification of ancestry, GWAS model and sample-size scaling.

**Legacy disposition.** Open; no answer or verified code change.

**LDSC3 assessment — Methodological limitation.** LDSC3 can analyze suitable UKB summary statistics using the same framework, but binary/quantitative scales, effective sample sizes, mixed-model corrections and ancestry must be assessed per GWAS release.

**Recommended prevention.** Provide a concise UKB input checklist tied to the actual association model and fields, not a blanket UKB flag preset.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[src/ldsc/_sumstats/munge_input.py:184](<../../../../src/ldsc/_sumstats/munge_input.py#L184>); [docs/reference/regression-options.md](<../../../reference/regression-options.md>); [docs/user-guide/workflows/cell-type-enrichment.md](<../../../user-guide/workflows/cell-type-enrichment.md>).

<a id="issue-474"></a>

### #474: Negative value of LDscore

[Original issue](https://github.com/bulik/ldsc/issues/474) · Source state: **closed** · Severity: **moderate** · Accessible comments reviewed: **0**.

Workflow tags: `ldscore`.

**Trigger and harm.** Downloaded baselineLD scores include negative values that seem impossible under a raw r-squared definition.
Clamping legitimate unbiased estimates changes the regression predictor and its bias properties.

**Legacy disposition.** Closed without an accessible response or code fix.

**LDSC3 assessment — Methodological limitation.** LDSC3 intentionally uses sample-size-corrected unbiased R2, which can be negative, and sums annotation-weighted values.
The naive nonnegative population definition is not the finite-sample estimator used by the software.

**Recommended prevention.** Explain the estimator formula, retain negatives and distinguish a few expected values from a broad pattern requiring panel/window QC.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[src/ldsc/_refpanel/query.py](<../../../../src/ldsc/_refpanel/query.py>); [docs/dev/DOCS_VOCABULARY.md](<../../DOCS_VOCABULARY.md>); [docs/reference/ldscore-snp-universes.md](<../../../reference/ldscore-snp-universes.md>).

<a id="issue-475"></a>

### #475: LDSC negative heritability

[Original issue](https://github.com/bulik/ldsc/issues/475) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **1**.

Workflow tags: `h2`.

**Trigger and harm.** An EAS stroke GWAS has mean chi-square near one and h2=-0.0024 with SE0.0059.
Users can treat sampling noise as a negative biological variance or automatically alter the model.

**Legacy disposition.** A commenter explains the estimate is not distinguishable from zero; no code fix is needed.

**LDSC3 assessment — Methodological limitation.** LDSC3 preserves unconstrained estimates, which are useful for assessing precision.
The result does not establish a software defect or identify which study-processing factor caused weak signal.

**Recommended prevention.** Present estimate with SE and signal diagnostics, avoid clipping or significance-seeking constraints, and distinguish this case from a precisely large negative estimate.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[src/ldsc/_regression/estimators.py:463](<../../../../src/ldsc/_regression/estimators.py#L463>); [docs/user-guide/workflows/heritability.md](<../../../user-guide/workflows/heritability.md>).

<a id="issue-481"></a>

### #481: what is the value in parenthesis in Total Observed scale h2?

[Original issue](https://github.com/bulik/ldsc/issues/481) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **3**.

Workflow tags: `h2`.

**Trigger and harm.** Asks what the parenthesized value means and how to know a heritability estimate is sound.
Statistical significance can be confused with correct scale, estimand and data processing.

**Legacy disposition.** Comments identify SE and recommend independent estimates; no code fix is involved.

**LDSC3 assessment — Methodological limitation.** LDSC3 labels SE columns, but comparison to family-based or individual-level estimates requires matching variant set, population, trait definition and scale.
Such comparators are not automatically ground truth for common-SNP heritability.

**Recommended prevention.** Explain SE, interval and diagnostics, then give a validation checklist with matched estimands rather than a significance threshold or an unqualified twin-study comparison.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[docs/user-guide/workflows/heritability.md](<../../../user-guide/workflows/heritability.md>); [docs/user-guide/concepts/ldsc-vs-partitioned-h2.md](<../../../user-guide/concepts/ldsc-vs-partitioned-h2.md>); [src/ldsc/_regression/estimators.py:404](<../../../../src/ldsc/_regression/estimators.py#L404>).

<a id="issue-494"></a>

### #494: negative MAF/LD Score Correlation

[Original issue](https://github.com/bulik/ldsc/issues/494) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **0**.

Workflow tags: `ldscore`.

**Trigger and harm.** Some chromosomes have small negative correlations between MAF and estimated LD scores.
Users may treat a weak diagnostic correlation as an automatic pass/fail test and alter scores to force it positive.

**Legacy disposition.** Open with no panel/QC details or resolution.

**LDSC3 assessment — Methodological limitation.** LDSC3 preserves unbiased R2 and panel-specific MAF.
The expected diagnostic depends on population, sample size, LD window and filtering, so a sign alone cannot establish score validity or a software error.

**Recommended prevention.** Inspect panel sample count, polymorphism, CM/map correctness, window span and score distribution jointly; do not clamp scores or require positive correlation by fiat.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[src/ldsc/_refpanel/query.py](<../../../../src/ldsc/_refpanel/query.py>); [docs/reference/ldscore-snp-universes.md](<../../../reference/ldscore-snp-universes.md>); [docs/dev/design/memory-and-cpu.md](<../../design/memory-and-cpu.md>).

<a id="issue-501"></a>

### #501: &lt;200,000 variants

[Original issue](https://github.com/bulik/ldsc/issues/501) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **0**.

Workflow tags: `ldscore`, `rg`.

**Trigger and harm.** Asks whether the historical 200k recommendation applies to LD computation or genetic correlation.
Users may conflate the number of reference variants, stored rows and effective regression information.

**Legacy disposition.** Open with no clarification.

**LDSC3 assessment — Methodological limitation.** LDSC3 separates reference and regression SNPs.
The historical threshold is a regression-coverage diagnostic, not a mathematical minimum for computing an LD score or a universal guarantee of precision above that number.

**Recommended prevention.** Explain the stage, signal/SE and block/model constraints, and show counts by funnel; do not enforce a universal 200k gate on tiny demonstration calculations.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[docs/reference/ldscore-snp-universes.md](<../../../reference/ldscore-snp-universes.md>); [docs/user-guide/workflows/genetic-correlation.md](<../../../user-guide/workflows/genetic-correlation.md>); [Source](<https://github.com/bulik/ldsc/wiki/FAQ>).

<a id="faq-04"></a>

### FAQ-04: Negative unbiased LD scores

[FAQ source](https://github.com/bulik/ldsc/wiki/FAQ) · Source-order identifier: **FAQ-04** · Severity: **moderate**.

Workflow tags: `ldscore`.

**Legacy disposition.** Historical FAQ guidance; assessed independently rather than treated as a code-fix record.

**LDSC3 assessment — Methodological limitation.** The reference implementation deliberately preserves negative unbiased R² contributions. A negative estimated LD score is not automatically an invalid input; widespread negativity calls for checking reference sample size, window geometry and input quality.

**Recommended prevention.** Explain population squared correlations versus finite-reference unbiased estimates; do not clip the latter to zero.

**Evidence and verification.** Source-inspected; relevant executed checks are identified in verification.md.
[src/ldsc/_refpanel/query.py](<../../../../src/ldsc/_refpanel/query.py>); [src/ldsc/_ldscore/compute.py](<../../../../src/ldsc/_ldscore/compute.py>); [docs/dev/DOCS_VOCABULARY.md](<../../DOCS_VOCABULARY.md>).

<a id="faq-08"></a>

### FAQ-08: Numeric annotation sums

[FAQ source](https://github.com/bulik/ldsc/wiki/FAQ) · Source-order identifier: **FAQ-08** · Severity: **moderate**.

Workflow tags: `ldscore`, `continuous`.

**Legacy disposition.** Historical FAQ guidance; assessed independently rather than treated as a code-fix record.

**LDSC3 assessment — Methodological limitation.** LDSC3 accumulates supplied numeric annotation values in float64.
For non-binary annotations these are weighted sums, not literal numbers of SNPs belonging to a set, and binary enrichment language can therefore mislead.

**Recommended prevention.** Name the weighted quantity explicitly; use coefficients, standardized coefficients or post-fit quantiles for continuous annotations.

**Evidence and verification.** Source-inspected; relevant executed checks are identified in verification.md.
[src/ldsc/_annotation/overlap.py](<../../../../src/ldsc/_annotation/overlap.py>); [docs/reference/results/quantile-h2.md](<../../../reference/results/quantile-h2.md>).

<a id="faq-23"></a>

### FAQ-23: Imprecise heritability estimates

[FAQ source](https://github.com/bulik/ldsc/wiki/FAQ) · Source-order identifier: **FAQ-23** · Severity: **moderate**.

Workflow tags: `h2`, `baseline-partitioned-h2`.

**Legacy disposition.** Historical FAQ guidance; assessed independently rather than treated as a code-fix record.

**LDSC3 assessment — Methodological limitation.** More annotations can increase uncertainty because the model must separate correlated effects.
Fixing the intercept may reduce variance only under a defensible constraint and can introduce bias otherwise.

**Recommended prevention.** Show uncertainty and model dimensionality; simplify unsupported models or obtain more information rather than automatically constraining the intercept.

**Evidence and verification.** Source-inspected; relevant executed checks are identified in verification.md.
[src/ldsc/_regression/jackknife.py](<../../../../src/ldsc/_regression/jackknife.py>); [src/ldsc/_regression/overlap_model.py:79](<../../../../src/ldsc/_regression/overlap_model.py#L79>); [Source](<https://pmc.ncbi.nlm.nih.gov/articles/PMC4626285/>).

<a id="faq-26"></a>

### FAQ-26: Self-correlation numerical identity

[FAQ source](https://github.com/bulik/ldsc/wiki/FAQ) · Source-order identifier: **FAQ-26** · Severity: **moderate**.

Workflow tags: `rg`.

**Legacy disposition.** Historical FAQ guidance; assessed independently rather than treated as a code-fix record.

**LDSC3 assessment — Methodological limitation.** Current Hsq and Gencov use distinct weight functions and initialization with a fixed number of reweightings.
Identical inputs do not by themselves justify an exact floating-point equality assertion for the final ratio.

**Recommended prevention.** Use a scientifically matched, well-conditioned self-correlation check with a stated tolerance; do not interpret small algorithmic differences as different biology.

**Evidence and verification.** Source-inspected; relevant executed checks are identified in verification.md.
[src/ldsc/_regression/estimators.py:Hsq.weights](<../../../../src/ldsc/_regression/estimators.py>); [src/ldsc/_regression/estimators.py:Gencov.weights](<../../../../src/ldsc/_regression/estimators.py>); [src/ldsc/_regression/irwls.py](<../../../../src/ldsc/_regression/irwls.py>).

<a id="faq-27"></a>

### FAQ-27: Same label across different cohorts

[FAQ source](https://github.com/bulik/ldsc/wiki/FAQ) · Source-order identifier: **FAQ-27** · Severity: **moderate**.

Workflow tags: `rg`.

**Legacy disposition.** Historical FAQ guidance; assessed independently rather than treated as a code-fix record.

**LDSC3 assessment — Methodological limitation.** A shared trait label does not ensure identical phenotype definitions, ascertainment or genetic architecture across cohorts.
Sampling uncertainty and reference compatibility also affect the comparison.

**Recommended prevention.** Document population and phenotype definitions and interpret the estimated correlation with uncertainty rather than assume a target of exactly one.

**Evidence and verification.** Source-inspected; relevant executed checks are identified in verification.md.
[src/ldsc/_regression/estimators.py:RG](<../../../../src/ldsc/_regression/estimators.py>); [Source](<https://www.nature.com/articles/ng.3406>).

<a id="issue-390"></a>

### #390: How to test the null hypothesis of rg=1

[Original issue](https://github.com/bulik/ldsc/issues/390) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **0**.

Workflow tags: `rg`.

**Trigger and harm.** Asks for a test of complete genetic correlation rather than the default zero-correlation test. A significant default p-value does not reject rg=1 and cannot prove incomplete sharing.

**Legacy disposition.** Open with no answer or code fix.

**LDSC3 assessment — Intentionally unsupported.** LDSC3 reports a zero-null normal test. A separate contrast (rg-1)/SE is an approximate inference that needs appropriate uncertainty and boundary interpretation; it is not obtained by reading the existing p column.

**Recommended prevention.** Document the default null explicitly and, if demanded, add a clearly labeled contrast against a specified null with methodological qualifications rather than repurposing the existing p-value.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[src/ldsc/_regression/estimators.py:65](<../../../../src/ldsc/_regression/estimators.py#L65>); [src/ldsc/_regression/estimators.py:779](<../../../../src/ldsc/_regression/estimators.py#L779>); [docs/user-guide/workflows/genetic-correlation.md](<../../../user-guide/workflows/genetic-correlation.md>).

<a id="issue-256"></a>

### #256: Do I need convert a list of SNP coordinates from hg19 to hg18

[Original issue](https://github.com/bulik/ldsc/issues/256) · Source state: **closed** · Severity: **moderate** · Accessible comments reviewed: **0**.

Workflow tags: `munge-sumstats`.

**Trigger and harm.** User thinks an HM3 rsID list requires converting hg19 coordinates to hg18.
Unnecessary or wrong-direction liftover can destroy valid matching.

**Legacy disposition.** Closed without explanation or fix.

**LDSC3 assessment — Apparently prevented.** LDSC3 explicitly separates rsID-family identity from coordinate-family identity.
An rsID membership list is not itself a set of hg18 positions; coordinate liftover is needed when coordinate-based sources differ in genome build, not merely because the list is old.
It supports hg19/hg38, not a general hg18 workflow.

**Recommended prevention.** Explain which fields define identity and support boundaries; require source/target build only for coordinate conversion.

**Evidence and verification.** source-inspected
[src/ldsc/_genome/build_policy.py](<../../../../src/ldsc/_genome/build_policy.py>); [docs/reference/snp-identity-and-genome-build.md](<../../../reference/snp-identity-and-genome-build.md>).

<a id="issue-270"></a>

### #270: calculation of tau* from LDSC output

[Original issue](https://github.com/bulik/ldsc/issues/270) · Source state: **closed** · Severity: **moderate** · Accessible comments reviewed: **1**.

Workflow tags: `baseline-partitioned-h2`, `continuous`.

**Trigger and harm.** Users ask how Coefficient and SE map to the standardized coefficient formula.
Raw conditional per-unit effects are compared across annotations with different scales.

**Legacy disposition.** Closed without a definitive calculation; later users remain confused.

**LDSC3 assessment — Apparently prevented.** LDSC3 retains raw coefficients and provides standardized_coefficients.tsv through quantile-h2 using the annotation standard deviation and total fitted-model SNP heritability. A coefficient SE alone is insufficient to reconstruct the scale factor, and nonpositive total SNP heritability makes standardization undefined.

**Recommended prevention.** Document the standardization formula, eligible common reference SNPs and unavailable-value policy; distinguish raw coefficient tests from standardized magnitudes.

**Evidence and verification.** source-inspected
[docs/user-guide/workflows/continuous-annotations.md](<../../../user-guide/workflows/continuous-annotations.md>); [docs/reference/results/quantile-h2.md](<../../../reference/results/quantile-h2.md>); [docs/reference/results/partitioned-h2.md](<../../../reference/results/partitioned-h2.md>).

<a id="issue-285"></a>

### #285: Different h2 results for the same phenotype

[Original issue](https://github.com/bulik/ldsc/issues/285) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **0**.

Workflow tags: `h2`, `rg`.

**Trigger and harm.** The same trait has estimates0.081,0.082and0.083 when paired with different partners.
Users interpret expected pairwise-intersection differences as nondeterminism.

**Legacy disposition.** Open without an established code failure.

**LDSC3 assessment — Apparently prevented.** LDSC3 reports h2_per_trait.tsv from each trait own SNPs and pair-specific heritabilities in rg_full.tsv.
Each pair can retain different SNPs and alleles, so the fitted denominator can differ.
The separate outputs make the intended comparison explicit rather than demanding identical estimates across unequal inputs.

**Recommended prevention.** Explain own-trait versus pairwise populations and report retained SNP counts; compare on the same population when investigating numerical differences.

**Evidence and verification.** source-inspected
[docs/user-guide/workflows/genetic-correlation.md](<../../../user-guide/workflows/genetic-correlation.md>); [src/ldsc/_regression/runner.py:build_rg_dataset](<../../../../src/ldsc/_regression/runner.py>).

<a id="issue-336"></a>

### #336: tau* z-scores

[Original issue](https://github.com/bulik/ldsc/issues/336) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **0**.

Workflow tags: `continuous`.

**Trigger and harm.** User asks how raw coefficient SE is used for the standardized-coefficient z-score.
Rescaled effects may be tested with an inconsistent SE or mistaken for different evidence.

**Legacy disposition.** Open without an answer or code fix.

**LDSC3 assessment — Apparently prevented.** LDSC3 fixed-scale standardization multiplies raw coefficient and SE by the same positive annotation-SD/per-SNP-heritability factor, so nondegenerate standardized and raw z-scores agree.
This convention does not propagate a new denominator ratio uncertainty and is distinct from the one-sided query coefficient p-value in partitioned_h2.tsv.

**Recommended prevention.** Show the shared scale algebra and fixed-scale convention; report zero-SD/nonpositive-total cases as undefined and distinguish test tails.

**Evidence and verification.** source-inspected
[src/ldsc/_regression/quantile_h2.py:_standardized_coefficients](<../../../../src/ldsc/_regression/quantile_h2.py>); [docs/reference/results/quantile-h2.md](<../../../reference/results/quantile-h2.md>).

<a id="issue-280"></a>

### #280: Little variations in intercept value

[Original issue](https://github.com/bulik/ldsc/issues/280) · Source state: **open** · Severity: **low** · Accessible comments reviewed: **0**.

Workflow tags: `h2`, `rg`.

**Trigger and harm.** A user obtains slightly different intercepts from historical tutorial screenshots.
Harmless numerical/data-version differences can be mistaken for an invalid installation.

**Legacy disposition.** Open without numerical tables or a reproducible cause; screenshots alone are insufficient to establish a bug.

**LDSC3 assessment — Indeterminate.** Both screenshots were inspected: the displayed mean chi-square differs (1.2336 versus 1.2234), as do heritability, intercept and correlation.
This is evidence that the compared inputs or filtering differ, rather than mere final rounding on an identical fitted population.
LDSC3 has explicit input/fit-population diagnostics, but the original command and retained rows are unavailable, so the particular cause cannot be diagnosed.

**Recommended prevention.** Publish versioned example inputs and numerical tolerances; investigate material differences using retained counts and policies before changing code.

**Evidence and verification.** Original screenshots visually inspected; source inspected; original data not reproduced.
`tests/fixtures/golden/`; [docs/reference/regression-options.md](<../../../reference/regression-options.md>); [docs/user-guide/guided-tutorial.md](<../../../user-guide/guided-tutorial.md>).

<a id="issue-215"></a>

### #215: Value meaning

[Original issue](https://github.com/bulik/ldsc/issues/215) · Source state: **closed** · Severity: **low** · Accessible comments reviewed: **2**.

Workflow tags: `h2`, `rg`.

**Trigger and harm.** Users cannot distinguish a point estimate from its parenthesized standard error.
Basic result uncertainty is misread or ignored.

**Legacy disposition.** Comment explains parentheses are standard errors; reporter acknowledges and closes.

**LDSC3 assessment — Apparently prevented.** LDSC3 scientific TSV outputs separate estimate and *_se columns, reducing dependence on an unexplained log convention.
SNP heritability, intercept and ratio SEs should still be defined near the first example so a novice understands confidence intervals and noisy estimates.

**Recommended prevention.** Include a short labeled estimate/SE example and distinguish SE from a confidence interval.

**Evidence and verification.** source-inspected
[src/ldsc/_regression/runner.py:summarize_total_h2](<../../../../src/ldsc/_regression/runner.py>); [docs/user-guide/workflows/heritability.md](<../../../user-guide/workflows/heritability.md>).

<a id="faq-28"></a>

### FAQ-28: Efficient block jackknife

[FAQ source](https://github.com/bulik/ldsc/wiki/FAQ) · Source-order identifier: **FAQ-28** · Severity: **low**.

Workflow tags: `h2`, `rg`, `baseline-partitioned-h2`, `ldsc-seg`.

**Legacy disposition.** Historical FAQ guidance; assessed independently rather than treated as a code-fix record.

**LDSC3 assessment — Apparently prevented.** Current least-squares jackknife combines block sufficient statistics and solves delete-block systems.
It is computationally efficient without assuming independent SNPs within a block, but block selection and singular deletions remain substantive constraints.

**Recommended prevention.** Keep an explanation of block dependence and compare the fast calculation with an explicit small deletion reference.

**Evidence and verification.** Source-inspected; relevant executed checks are identified in verification.md.
[src/ldsc/_regression/jackknife.py:LstsqJackknifeFast](<../../../../src/ldsc/_regression/jackknife.py>); [tests/regression/test_jackknife.py](<../../../../tests/regression/test_jackknife.py>).

