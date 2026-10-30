# Correctness and numerical bugs

Last updated on: 2026-10-28

[Audit overview](README.md) · [Workflow counts and membership](workflow-index.md) · [Verification](verification.md)

## Recurrent correctness and numerical risks

The most consequential reproduced gaps are early numeric validation and explicit handling of undefined ratio uncertainty.
They are detailed, with fixture definitions and limits, in [V1–V3](verification.md).
The executed focused suite passed, which shows that these boundary cases are not covered by those passing checks.
It does not erase existing protections or establish that all historical bugs persist.

| Demonstrated mechanism | Current evidence | Proposed canonical solution | Future regression scenario |
| --- | --- | --- | --- |
| Invalid direct-import statistics/sample sizes reach fitting | V1 uses the real text loader and dataset construction; invalid N can still produce a returned fit. | Require finite statistics and strictly positive finite N at the accepted-input boundary, with offending-row diagnostics. | Text/native inputs with infinity, zero and negative N, plus valid finite signed statistics. |
| Genetic-correlation SE is undefined but summary says `ok` | V2 uses real fits and current summary construction; invalid delete-block products yield `NaN` uncertainty. | Explicit not-estimable status with trait/block diagnostics; preserve valid components rather than fabricate a ratio. | One/both nonpositive delete-block heritabilities, zero and near-zero denominators, and valid weak-signal controls. |
| Zero total prevents a usable summary | V3 reaches a derived ratio after a zero-total fit. | Return valid primary quantities and explicit unavailable derived fields; avoid epsilon denominators. | Exact-zero total, zero deleted total, and small nonzero totals. |

Other entries below include source-level residuals, apparently prevented historical defects and cases without sufficient original evidence to diagnose.
In particular, distinct weighting policies or intended LD-window changes must not be labelled a refactoring regression merely because the answers differ.
A matched comparison must control SNPs, counts, intercepts, weights, LD window and other scientific settings.

## Recommendations by kind

- **Code:** use the existing reader, fit and summary seams for the smallest repair; avoid imputation, pseudoinverses, silent block deletion and new recovery frameworks.
- **Documentation:** identify what failed, which scientific quantity is unavailable, and which valid results remain interpretable.
- **Features:** expose existing invalid-fit information before adding new estimators or compatibility modes.
- **Resources:** keep alignment/count/window mismatches distinct from numerical-kernel failures.

## Category summary

83 issues and 1 FAQ entry.

| Severity | Source | Concern | LDSC3 disposition | Evidence level |
| --- | --- | --- | --- | --- |
| High | [#40](#issue-40) | Invalid delete-block heritability denominator | Confirmed persistent behavior | Source + shared boundary probe |
| High | [#54](#issue-54) | Weak-signal correlation square-root failure | Confirmed persistent behavior | Source + shared boundary probe |
| High | [#117](#issue-117) | Fixed intercept masks invalid correlation jackknife | Confirmed persistent behavior | Source + shared boundary probe |
| High | [#220](#issue-220) | Negative delete-block genetic-correlation denominator | Confirmed persistent behavior | Source + shared boundary probe |
| High | [#241](#issue-241) | Negative jackknife heritability product | Confirmed persistent behavior | Source + shared boundary probe |
| High | [#290](#issue-290) | Undefined ratio jackknife despite positive standalone heritability | Confirmed persistent behavior | Source + shared boundary probe |
| High | [#58](#issue-58) | Whitespace parsing shifts empty tab fields | Partial coverage | Source inspection |
| High | [#66](#issue-66) | Empty munging result has several causes | Partial coverage | Source inspection |
| High | [#144](#issue-144) | Extremely small P values underflow or parse as text | Partial coverage | Source inspection |
| High | [#165](#issue-165) | Negative computed total variance after a modified baseline | Partial coverage | Source inspection |
| High | [#174](#issue-174) | Nonfinite values reach regression weights | Partial coverage | Source + input-boundary probe |
| High | [#181](#issue-181) | Munging stalls at large chunk size | Partial coverage | Source inspection |
| High | [#207](#issue-207) | Multi-annotation heritability differs inside genetic correlation | Partial coverage | Source inspection |
| High | [#228](#issue-228) | Invalid arithmetic in heritability weight update | Partial coverage | Source + input-boundary probe |
| High | [#242](#issue-242) | Extreme values become nonnumeric or underflow | Partial coverage | Source inspection |
| High | [#258](#issue-258) | Very large PLINK panel memory | Partial coverage | Source inspection |
| High | [#263](#issue-263) | Huge LD block from chunk-size misunderstanding | Partial coverage | Source inspection |
| High | [#268](#issue-268) | Munging takes eight to seventy-two hours | Partial coverage | Source inspection |
| High | [#275](#issue-275) | Extremely small p-values below floating range | Partial coverage | Source inspection |
| High | [#293](#issue-293) | P-values around five times10^-700 | Partial coverage | Source inspection |
| High | [#312](#issue-312) | Twenty-four-hour run reduced by chunk size | Partial coverage | Source inspection |
| High | [#323](#issue-323) | Multi-annotation RG uses different heritability weights | Partial coverage | Source inspection |
| High | [#418](#issue-418) | Extreme p-value underflow | Partial coverage | Source inspection |
| High | [#496](#issue-496) | NA genetic correlation and estimability status | Partial coverage | Source + shared boundary probe |
| High | [FAQ-09](#faq-09) | Collinear annotation designs | Partial coverage | Source inspection |
| High | [#22](#issue-22) | Historical proportion-of-heritability calculation | Indeterminate | Source inspection |
| High | [#292](#issue-292) | Abrupt munging runtime cliff near one million rows | Indeterminate | Source inspection |
| High | [#10](#issue-10) | Genetic-covariance regression variance factor | Apparently prevented | Source inspection |
| High | [#11](#issue-11) | Annotation pruning and overlap/count alignment | Apparently prevented | Source inspection |
| High | [#17](#issue-17) | Missing high-statistic/two-step regression handling | Apparently prevented | Source inspection |
| High | [#24](#issue-24) | Incompatible allele four-mer lookup | Apparently prevented | Source inspection |
| High | [#35](#issue-35) | Fixed intercept produced grossly inflated heritability | Apparently prevented | Source inspection |
| High | [#72](#issue-72) | Fast merge compares differently indexed Series | Apparently prevented | Source inspection |
| High | [#77](#issue-77) | Overlap matrix/count dimensions mismatch | Apparently prevented | Source inspection |
| High | [#84](#issue-84) | Noninteger block slice indices | Apparently prevented | Source inspection |
| High | [#85](#issue-85) | Frequency mask applied to a different SNP table | Apparently prevented | Source inspection |
| High | [#97](#issue-97) | Redundant chromosome/MAF partitions destabilize covariance | Apparently prevented | Source inspection |
| High | [#99](#issue-99) | Wrong prevalence used in multitraits liability table | Apparently prevented | Source inspection |
| High | [#103](#issue-103) | Phase 3 annotations combined with Phase 1 frequencies | Apparently prevented | Source inspection |
| High | [#114](#issue-114) | Uninitialized PLINK padding after individual filtering | Apparently prevented | Source inspection |
| High | [#126](#issue-126) | Duplicate BIM positions expand annotation rows | Apparently prevented | Source inspection |
| High | [#137](#issue-137) | Filtered allele keys not aligned to retained traits | Apparently prevented | Source inspection |
| High | [#139](#issue-139) | Repeated incompatible-allele lookup report | Apparently prevented | Source inspection |
| High | [#140](#issue-140) | Allele filter applied to data but not allele vector | Apparently prevented | Source inspection |
| High | [#158](#issue-158) | Singular models from duplicate bases, empty queries or collinearity | Apparently prevented | Source inspection |
| High | [#186](#issue-186) | cM window produces identical scores | Apparently prevented | Source inspection |
| High | [#214](#issue-214) | Globbing chromosome-prefix regression | Apparently prevented | Source inspection |
| High | [#223](#issue-223) | Empty concatenation after glob update | Apparently prevented | Source inspection |
| High | [#234](#issue-234) | Zero-variance columns misalign overlap | Apparently prevented | Source inspection |
| High | [#306](#issue-306) | Overlapping intervals duplicate annotation rows | Apparently prevented | Source inspection |
| High | [#327](#issue-327) | No surviving chunks after QC/restriction | Apparently prevented | Source inspection |
| High | [#342](#issue-342) | LD-column/count order | Apparently prevented | Source inspection |
| High | [#356](#issue-356) | Integer mask inversion | Apparently prevented | Source inspection |
| High | [#493](#issue-493) | Coordinate-only join duplicates rows | Apparently prevented | Source inspection |
| High | [#502](#issue-502) | Annotation join multiplies rows | Apparently prevented | Source inspection |
| Moderate | [#227](#issue-227) | Quoted header cannot map required SNP | Confirmed persistent behavior | Quoted-header CLI reproduction |
| Moderate | [#2](#issue-2) | Duplicate rsIDs in merged inputs | Partial coverage | Source inspection |
| Moderate | [#187](#issue-187) | Munging never finishes | Partial coverage | Source inspection |
| Moderate | [#216](#issue-216) | Sixteen-million-row munging stalls | Partial coverage | Source inspection |
| Moderate | [#243](#issue-243) | Compressed-file parsing failures | Partial coverage | Source inspection |
| Moderate | [#269](#issue-269) | Gzip suffix/content mismatch | Partial coverage | Source inspection |
| Moderate | [#329](#issue-329) | Nonnumeric token inside numeric data | Partial coverage | Source inspection |
| Moderate | [#345](#issue-345) | Undefined enrichment denominator | Partial coverage | Source + shared boundary probe |
| Moderate | [#211](#issue-211) | Jackknife variance sqrt failure | Indeterminate | Source inspection |
| Moderate | [#19](#issue-19) | Malformed annotation yields an opaque shape error | Apparently prevented | Source inspection |
| Moderate | [#21](#issue-21) | Identical annotation names across LD-score sources | Apparently prevented | Source inspection |
| Moderate | [#28](#issue-28) | No recognized signed statistic causes IndexError | Apparently prevented | Source inspection |
| Moderate | [#34](#issue-34) | Annotation absent on one chromosome | Apparently prevented | Source inspection |
| Moderate | [#47](#issue-47) | Duplicate identifier aliases yield ambiguous Series | Apparently prevented | Source inspection |
| Moderate | [#128](#issue-128) | Annotation columns differ between chromosomes | Apparently prevented | Source inspection |
| Moderate | [#160](#issue-160) | Explicit SNP choice fails to override another alias | Apparently prevented | Source inspection |
| Moderate | [#175](#issue-175) | Boolean inversion produces integer indexes | Apparently prevented | Source inspection |
| Moderate | [#198](#issue-198) | Negative BED start after gene padding | Apparently prevented | Source inspection |
| Moderate | [#232](#issue-232) | Undescribed new glob regression | Apparently prevented | Source inspection |
| Moderate | [#239](#issue-239) | Empty logs on repeated runs | Apparently prevented | Source inspection |
| Moderate | [#245](#issue-245) | Allele-merge Boolean mask KeyError | Apparently prevented | Source inspection |
| Moderate | [#272](#issue-272) | Explicit column mapping should override aliases | Apparently prevented | Source inspection |
| Moderate | [#367](#issue-367) | Per-source chromosome discovery | Apparently prevented | Source inspection |
| Moderate | [#422](#issue-422) | Boolean Series with iloc | Apparently prevented | Source inspection |
| Moderate | [#435](#issue-435) | Comprehension-scope parser bug | Apparently prevented | Source inspection |
| Moderate | [#483](#issue-483) | SNP-list delimiter/header mismatch | Apparently prevented | Source inspection |
| Moderate | [#500](#issue-500) | Gzip header bytes/text | Apparently prevented | Source inspection |
| Low | [#6](#issue-6) | Historical simulation generator returned NAs | Intentionally unsupported | Source inspection |
| Low | [#244](#issue-244) | Unspecified least-squares rcond | Apparently prevented | Source inspection |

## Individual assessments

<a id="issue-40"></a>

### #40: error in rg

[Original issue](https://github.com/bulik/ldsc/issues/40) · Source state: **closed** · Severity: **high** · Accessible comments reviewed: **5**.

Workflow tags: `rg`.

**Trigger and harm.** Full h2 estimates are positive but one or more delete-block h2 products are negative.
Correlation SE becomes NaN or a floating-point error discards an otherwise usable genetic-covariance fit.

**Legacy disposition.** Maintainers diagnosed the unguarded jackknife square root and suggested unavailable correlation plus retained covariance; no committed legacy fix is identified.

**LDSC3 assessment — Confirmed persistent behavior.** LDSC3 RG.__init__ still checks only full totals and then takes sqrt of delete-block products.
Root probe V2 reproduces positive full totals, five negative products and NaN correlation SE.
Both negative h2 deletions also need individual validity checks even though their product is positive.
This is a numerical reporting defect atop a genuine low-signal limitation.
The default-runner probe additionally returns status=ok with empty error despite NaN SE/P, so downstream status is also misleading.

**Recommended prevention.** Validate each delete-block denominator and report unavailable correlation/SE with explicit block diagnostics while retaining genetic covariance; never clip negative deletions to fabricate precision.

**Evidence and verification.** Root probe V2 reproduced invalid square root and NaN rg_se; source-inspected.
[src/ldsc/_regression/estimators.py:RG.__init__](<../../../../src/ldsc/_regression/estimators.py>); [docs/dev/validation/legacy-issues-audit/verification.md#v2-genetic-correlation-delete-block-boundary](<verification.md#v2-genetic-correlation-delete-block-boundary>).

Legacy implementation context: [`ldscore/regressions.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/regressions.py), [`ldscore/sumstats.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/sumstats.py).
[Shared reproduced boundary V2](verification.md#v2-genetic-correlation-delete-block-boundary); the original historical dataset was not reproduced.

<a id="issue-54"></a>

### #54: ERROR computing rg for phenotype 2/2, FloatingPointError: invalid value encountered in sqrt

[Original issue](https://github.com/bulik/ldsc/issues/54) · Source state: **closed** · Severity: **high** · Accessible comments reviewed: **9**.

Workflow tags: `rg`.

**Trigger and harm.** A low-signal sex-stratified trait has positive full h2 but invalid delete-block denominator; some loci are very strong.
The fit loses correlation uncertainty and users misdiagnose positive odds ratios as the cause.

**Legacy disposition.** Maintainer explains low h2 Z-score and extreme-locus sensitivity; removing top loci helped one user, without a verified general code repair.

**LDSC3 assessment — Confirmed persistent behavior.** LDSC3 retains the unguarded delete-block sqrt documented for #40; root probe V2 reproduces the mechanism.
Positive OR values are expected, while signed Z is obtained relative to OR=1.
Excluding loci after seeing results changes the scientific target and must be a reported sensitivity analysis, not a universal preprocessing fix.

**Recommended prevention.** Guard each delete-block denominator, retain covariance output, and explain that valid file format cannot compensate for insufficient polygenic signal.

**Evidence and verification.** Root probe V2; source-inspected; original attachments not independently reproduced.
[src/ldsc/_regression/estimators.py:RG.__init__](<../../../../src/ldsc/_regression/estimators.py>); [src/ldsc/_sumstats/munge_kernel.py:p_to_z](<../../../../src/ldsc/_sumstats/munge_kernel.py>); [docs/dev/validation/legacy-issues-audit/verification.md#v2-genetic-correlation-delete-block-boundary](<verification.md#v2-genetic-correlation-delete-block-boundary>).

Legacy implementation context: [`ldscore/regressions.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/regressions.py), [`ldscore/sumstats.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/sumstats.py).
[Shared reproduced boundary V2](verification.md#v2-genetic-correlation-delete-block-boundary); the original historical dataset was not reproduced.

<a id="issue-117"></a>

### #117: LD score regression trait correlation ERROR:invalid value encountered in sqrt

[Original issue](https://github.com/bulik/ldsc/issues/117) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **10**.

Workflow tags: `rg`, `h2`.

**Trigger and harm.** A free-intercept fit hits negative delete-block h2 while a constrained fit returns a plausible number.
Users may choose constraints to force completion and compare SNP heritability with twin estimates on another scale.

**Legacy disposition.** Maintainer diagnoses the denominator, warns that no overlap does not imply no confounding, and user reports strong-locus exclusion helps.
No general legacy code fix is established.

**LDSC3 assessment — Confirmed persistent behavior.** LDSC3 retains the invalid-square-root path.
Root V2 also confirms default runner output can have NaN SE/P while status is ok and error empty.
Strong full h2 does not rule out a problematic deletion or locus architecture.
Liability and twin-study comparisons concern different targets and scales.

**Recommended prevention.** Guard individual delete-block denominators, propagate unavailable status, preserve covariance, and make locus exclusion a transparent sensitivity analysis rather than choose intercepts for successful completion.

**Evidence and verification.** Root V2 reproduces NaN uncertainty with status ok; source-inspected.
[src/ldsc/_regression/estimators.py:RG.__init__](<../../../../src/ldsc/_regression/estimators.py>); [src/ldsc/_regression/runner.py:_summarize_rg_pair](<../../../../src/ldsc/_regression/runner.py>); [docs/dev/validation/legacy-issues-audit/verification.md#v2-genetic-correlation-delete-block-boundary](<verification.md#v2-genetic-correlation-delete-block-boundary>).

Legacy implementation context: [`ldscore/regressions.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/regressions.py), [`ldscore/sumstats.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/sumstats.py).
[Shared reproduced boundary V2](verification.md#v2-genetic-correlation-delete-block-boundary); the original historical dataset was not reproduced.

<a id="issue-220"></a>

### #220: LD score regression trait correlation ERROR: FloatingPointError: invalid value encountered in sqrt

[Original issue](https://github.com/bulik/ldsc/issues/220) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **1**.

Workflow tags: `munge-sumstats`, `h2`, `rg`.

**Trigger and harm.** A small study has weak SNP heritability; sqrt of the product of delete-block heritabilities fails. A positive full-sample denominator does not guarantee all ratio-jackknife replicates are defined.

**Legacy disposition.** Open; another user reports related multiply/sqrt failures, without a fix.

**LDSC3 assessment — Confirmed persistent behavior.** LDSC3 RG checks full-sample heritabilities but directly takes sqrt of delete-block products.
Thus the numerical boundary remains even when individual full fits are positive.
The separate long-munging complaint shares the chunk-size issue; successful correlation with a third trait does not validate this pair.

**Recommended prevention.** Check each trait's delete-block heritability for finite positivity, then validate the denominator and uncertainty and report affected blocks; never force completion by clipping, dropping blocks or taking absolute values.

**Evidence and verification.** Source-inspected; the shared current boundary is reproduced by V2; original dataset not reproduced.
[src/ldsc/_regression/estimators.py:RG](<../../../../src/ldsc/_regression/estimators.py>); [src/ldsc/_sumstats/munge_kernel.py:parse_dat](<../../../../src/ldsc/_sumstats/munge_kernel.py>).

Legacy implementation context: [`munge_sumstats.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/munge_sumstats.py), [`ldscore/regressions.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/regressions.py), [`ldscore/sumstats.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/sumstats.py).
[Shared reproduced boundary V2](verification.md#v2-genetic-correlation-delete-block-boundary); the original historical dataset was not reproduced.

<a id="issue-241"></a>

### #241: estimate_rg (FloatingPointError: invalid value encountered in sqrt)

[Original issue](https://github.com/bulik/ldsc/issues/241) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **5**.

Workflow tags: `rg`.

**Trigger and harm.** A binary-trait pair fails at sqrt of delete-block SNP-heritability products.
Low-denominator signal produces undefined genetic-correlation uncertainty and a cryptic arithmetic error.

**Legacy disposition.** Comment identifies negative jackknife estimates as the mechanism; no code repair is provided.

**LDSC3 assessment — Confirmed persistent behavior.** LDSC3 RG still computes delete-block square roots after checking only full-sample estimates.
Logistic beta/SE ranges do not guarantee valid denominator replicates.
The scientific limitation is weak SNP-heritability information; the preventable software issue is failing without a direct denominator diagnosis.

**Recommended prevention.** Check each trait's delete-block heritability for finite positivity, then validate the denominator and uncertainty and report affected blocks; never force completion by clipping, dropping blocks or taking absolute values.

**Evidence and verification.** Source-inspected; the shared current boundary is reproduced by V2; original dataset not reproduced.
[src/ldsc/_regression/estimators.py:RG](<../../../../src/ldsc/_regression/estimators.py>); [Source](<https://github.com/bulik/ldsc/issues/241#issuecomment-2225810093>).

Legacy implementation context: [`ldscore/regressions.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/regressions.py), [`ldscore/sumstats.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/sumstats.py).
[Shared reproduced boundary V2](verification.md#v2-genetic-correlation-delete-block-boundary); the original historical dataset was not reproduced.

<a id="issue-290"></a>

### #290: Some results have NA values

[Original issue](https://github.com/bulik/ldsc/issues/290) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **1**.
The API reports 2 comments; one was unavailable through both retrieval methods.

Workflow tags: `rg`.

**Trigger and harm.** Several pairs yield NA and sqrt failures although standalone estimates are positive.
The full estimate does not guarantee all delete-block heritability products are positive.

**Legacy disposition.** Open without fix; one reported comment is unavailable and must be disclosed.

**LDSC3 assessment — Confirmed persistent behavior.** LDSC3 preserves the same unguarded delete-product square root, so the boundary persists.
Pairwise SNP intersection and noise can make a delete-block denominator cross zero even with positive own-trait fits.
An absent uncertainty estimate should not be presented as evidence for no genetic correlation.

**Recommended prevention.** Provide an explicit weak-denominator/block diagnostic and preserve failed-pair status; disclose the missing comment and avoid silent clipping.

**Evidence and verification.** source-inspected
[src/ldsc/_regression/estimators.py:RG](<../../../../src/ldsc/_regression/estimators.py>); [docs/user-guide/workflows/genetic-correlation.md](<../../../user-guide/workflows/genetic-correlation.md>).

Legacy implementation context: [`ldscore/regressions.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/regressions.py), [`ldscore/sumstats.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/sumstats.py).
[Shared reproduced boundary V2](verification.md#v2-genetic-correlation-delete-block-boundary); the original historical dataset was not reproduced.

<a id="issue-58"></a>

### #58: Validate expected dtypes in munge_sumstats

[Original issue](https://github.com/bulik/ldsc/issues/58) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **1**.

Workflow tags: `munge-sumstats`.

**Trigger and harm.** Consecutive tabs encode a missing value, but whitespace splitting collapses them and changes column alignment.
Values may be assigned to the wrong statistic or allele before a late SciPy dtype error.

**Legacy disposition.** Author proposes numeric-dtype checks and maintainer invites a PR; no merged fix is identified in the issue.

**LDSC3 assessment — Partial coverage.** LDSC3 checks numeric input columns before distribution functions, preventing the specific opaque chdtri error.
However numeric validation cannot detect every shifted row if displaced fields remain numeric.
The reader delimiter/missing-field policy therefore needs an explicit fixture and must not equate coercible numeric text with correct alignment.

**Recommended prevention.** Preserve empty delimited fields or reject ragged rows with line numbers; prioritize a public-input probe with an internal empty tab rather than only a dtype unit test.

**Evidence and verification.** source-inspected
[src/ldsc/_sumstats/munge_kernel.py:parse_dat](<../../../../src/ldsc/_sumstats/munge_kernel.py>); [src/ldsc/_sumstats/munge_kernel.py:numeric_cols](<../../../../src/ldsc/_sumstats/munge_kernel.py>); [docs/reference/formats/column-schema.md](<../../../reference/formats/column-schema.md>).

Legacy implementation context: [`munge_sumstats.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/munge_sumstats.py).

<a id="issue-66"></a>

### #66: Error converting summary statistics

[Original issue](https://github.com/bulik/ldsc/issues/66) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **20**.

Workflow tags: `munge-sumstats`, `rg`.

**Trigger and harm.** Header mismatch, an unused all-NaN case-count column, decorated rsIDs, or malformed delimiter removes every row; later correlation hits invalid sqrt. A generic concat error hides the cause and can trigger destructive manual reformatting.

**Legacy disposition.** Multiple users report different successful workarounds; no universal legacy fix is verified.
Later discussion explicitly diagnoses invalid jackknife h2.

**LDSC3 assessment — Partial coverage.** LDSC3 validates sample-size precedence, offers explicit SNP identifier modes, numeric checks and keep-list empty errors.
Early optional columns must not become required solely because their names match aliases.
The correlation subcase remains persistent via V2, so the overall issue cannot be labelled resolved because preprocessing improved.

**Recommended prevention.** Report the exact emptying stage and active columns, use source-verified identifier mapping, and handle unavailable correlation without clipping or rerunning a different target silently.

**Evidence and verification.** source-inspected
[src/ldsc/_sumstats/munge_input.py:resolve_sample_size](<../../../../src/ldsc/_sumstats/munge_input.py>); [src/ldsc/_sumstats/munge_kernel.py:parse_dat](<../../../../src/ldsc/_sumstats/munge_kernel.py>); [src/ldsc/_regression/estimators.py:RG.__init__](<../../../../src/ldsc/_regression/estimators.py>); [docs/dev/validation/legacy-issues-audit/verification.md#v2-genetic-correlation-delete-block-boundary](<verification.md#v2-genetic-correlation-delete-block-boundary>).

Legacy implementation context: [`munge_sumstats.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/munge_sumstats.py), [`ldscore/regressions.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/regressions.py), [`ldscore/sumstats.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/sumstats.py).

<a id="issue-144"></a>

### #144: Error converting summary statistics UKBiobank p-values

[Original issue](https://github.com/bulik/ldsc/issues/144) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **8**.

Workflow tags: `munge-sumstats`.

**Trigger and harm.** A valid P such as 1e-3739 is outside ordinary float64 numeric representation.
Munging drops strong variants or aborts with a misleading dtype error; manual capping changes scientific statistics.

**Legacy disposition.** PR #319 proposes coercion/exclusion but users state it was not merged; maintainer agrees extreme loci are often removed later, not that fabricated P values are safe.

**LDSC3 assessment — Partial coverage.** LDSC3 still parses ordinary numeric P and derives Z through chi2.isf.
Underflowed zero fails the P filter and object-valued extreme strings can fail numeric validation.
Its estimator-specific filters differ, so the historical claim all such rows would always be excluded later is too broad, especially ordinary single-annotation h2.

**Recommended prevention.** Give an explicit underflow diagnostic; prefer verified beta/SE or signed Z/log-P parsing when implemented and validated, and never silently cap P values.

**Evidence and verification.** source-inspected
[src/ldsc/_sumstats/munge_kernel.py:filter_pvals](<../../../../src/ldsc/_sumstats/munge_kernel.py>); [src/ldsc/_sumstats/munge_kernel.py:p_to_z](<../../../../src/ldsc/_sumstats/munge_kernel.py>); [src/ldsc/_regression/runner.py:_resolve_default_chisq_max](<../../../../src/ldsc/_regression/runner.py>); [Source](<https://github.com/bulik/ldsc/pull/319>).

Legacy implementation context: [`munge_sumstats.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/munge_sumstats.py).

<a id="issue-165"></a>

### #165: Cell type specific analysis crashes with specific traits using a modified baseline model

[Original issue](https://github.com/bulik/ldsc/issues/165) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **1**.

Workflow tags: `ldsc-seg`.

**Trigger and harm.** Only some traits fail at sqrt(tot_cov) in a custom baseline/query fit. A cancellation or unstable covariance calculation discards inference and hides whether the model is identifiable.

**Legacy disposition.** Open with no diagnosis; differing baseline behavior is evidence of model sensitivity, not proof of a particular defect.

**LDSC3 assessment — Partial coverage.** LDSC3 _tot still computes tot_cov as the sum of category covariance entries and immediately takes its square root.
It adds full-model collinearity and singular-delete diagnostics, but those do not mathematically prevent small negative roundoff in an otherwise valid covariance contrast or every nonfinite fit.
Original arrays are absent, so exact reproduction is unproved.

**Recommended prevention.** Check finite/rank context and compute uncertainty from the total delete values or a stable covariance contrast; only tolerate negative roundoff under an explicit validated criterion.

**Evidence and verification.** source-inspected
[src/ldsc/_regression/estimators.py:LD_Score_Regression._tot](<../../../../src/ldsc/_regression/estimators.py>); [src/ldsc/_regression/overlap_model.py](<../../../../src/ldsc/_regression/overlap_model.py>); [src/ldsc/_regression/jackknife.py](<../../../../src/ldsc/_regression/jackknife.py>).

Legacy implementation context: [`ldscore/regressions.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/regressions.py).

<a id="issue-174"></a>

### #174: FloatingPointError: invalid value encountered in less_equal

[Original issue](https://github.com/bulik/ldsc/issues/174) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **9**.

Workflow tags: `rg`.

**Trigger and harm.** Genetic-correlation fitting fails at w<=0 after merging to only39528 regression SNPs.
Missing annotation values or nonfinite summary statistics can produce a late opaque numerical failure.

**Legacy disposition.** Commenters found missing annotation cells or NaN/inf Z values; other reporters remained unresolved and no code fix is identified.

**LDSC3 assessment — Partial coverage.** LDSC3 rejects missing/non-numeric annotation cells, classifies fitted annotations as finite, and drops missing pairwise numeric fields.
Imported legacy summary statistics remain a gap: a real-loader/synthetic200-row LD-score fixture retained Z=inf and N=inf and fitting failed late with SVD did not converge.
This is a related input-boundary reproduction, not the unavailable original pair.

**Recommended prevention.** Require finite Z,N,LD scores and weights before fitting, with offending rows/columns; retain distinct low-overlap guidance.

**Evidence and verification.** actual boundary/workflow probe: load_sumstats text file -> build_dataset -> estimate_h2; synthetic LDScoreResult fixture, probe-nonfinite-b.py
[src/ldsc/_annotation/tables.py:_validate_annotation_values](<../../../../src/ldsc/_annotation/tables.py>); [src/ldsc/_regression/runner.py:RegressionRunner.build_rg_dataset](<../../../../src/ldsc/_regression/runner.py>).

Legacy implementation context: [`ldscore/regressions.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/regressions.py), [`ldscore/sumstats.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/sumstats.py).

<a id="issue-181"></a>

### #181: munge_sumstats.py never finished, then I had to use AWK

[Original issue](https://github.com/bulik/ldsc/issues/181) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **3**.

Workflow tags: `munge-sumstats`.

**Trigger and harm.** Five-million-row default chunks run for10-48hours; users also ask whether BETA/SE can replace Z.
Memory pressure blocks preprocessing and motivates unsafe ad hoc input reconstruction.

**Legacy disposition.** Two commenters confirm --chunksize500000 helps; no committed algorithmic fix is established.

**LDSC3 assessment — Partial coverage.** LDSC3 lowers the default chunk size to one million and exposes --chunksize, but parse_dat retains each surviving chunk then concatenates the entire table.
Memory therefore still grows with retained rows. P and signed direction remain required in raw munging; legacy regression accepts Z,N directly.

**Recommended prevention.** Document a smaller-chunk remedy and memory growth; add progress/memory evidence before considering out-of-core redesign.

**Evidence and verification.** source-inspected
[src/ldsc/_sumstats/munge_kernel.py:parse_dat](<../../../../src/ldsc/_sumstats/munge_kernel.py>); [src/ldsc/_sumstats/munger.py:1417](<../../../../src/ldsc/_sumstats/munger.py#L1417>).

Legacy implementation context: [`munge_sumstats.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/munge_sumstats.py).

<a id="issue-207"></a>

### #207: heritability obtain using --h2 and --rg are different.

[Original issue](https://github.com/bulik/ldsc/issues/207) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **2**.

Workflow tags: `h2`, `rg`.

**Trigger and harm.** The same trait has0.12 standalone SNP heritability but-2.33 inside a multi-annotation genetic-correlation fit.
Different SNP populations, filtering and weight policies can produce unstable denominators and incomparable reported estimates.

**Legacy disposition.** A commenter explains pairwise intersection, but does not explain the extreme discrepancy or establish a fix.

**LDSC3 assessment — Partial coverage.** LDSC3 records model defaults and pairwise intersections, so moderate differences can be expected.
Its standalone multi-annotation Hsq path uses legacy weights while RG internally constructs Hsq without old_weights, preserving an additional estimator-policy difference beyond SNP overlap.
Constraining intercepts to obtain plausible values is not a valid general solution.

**Recommended prevention.** Compare on identical SNPs and explicit policies; document or align multi-annotation weight behavior only after numerical validation.

**Evidence and verification.** source-inspected
[src/ldsc/_regression/runner.py:_fit_h2_dataset](<../../../../src/ldsc/_regression/runner.py>); [src/ldsc/_regression/runner.py:_fit_rg_dataset](<../../../../src/ldsc/_regression/runner.py>); [src/ldsc/_regression/estimators.py:RG](<../../../../src/ldsc/_regression/estimators.py>); [docs/reference/regression-options.md](<../../../reference/regression-options.md>).

Legacy implementation context: [`ldscore/regressions.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/regressions.py), [`ldscore/sumstats.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/sumstats.py).

<a id="issue-228"></a>

### #228: FloatingPointError: invalid value encountered in multiply

[Original issue](https://github.com/bulik/ldsc/issues/228) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **0**.

Workflow tags: `rg`.

**Trigger and harm.** One lifespan-healthspan pair fails at hsq*N/M while each trait correlates successfully withT2D.
Nonfinite or degenerate pair-specific fits can cause late low-level failures.

**Legacy disposition.** Open with no original input or verified fix; third-trait success does not establish valid pairwise inputs.

**LDSC3 assessment — Partial coverage.** LDSC3 preserves the weighting arithmetic, while direct legacy summary-statistics loading does not enforce finite positive N or finite Z. A synthetic loader/workflow probe retained infinity through dataset preparation and failed at linear algebra; this demonstrates a related boundary gap rather than reproducing the unavailable original pair.

**Recommended prevention.** Validate finite values and positive N before fitting; report pair-specific counts and denominator conditions with an actionable failure.

**Evidence and verification.** actual boundary/workflow probe: legacy text loader with Z/N infinity accepted; synthetic LDScoreResult fixture; original pair unavailable
[src/ldsc/_regression/estimators.py:Hsq.weights](<../../../../src/ldsc/_regression/estimators.py>); [src/ldsc/_sumstats/munger.py:load_sumstats](<../../../../src/ldsc/_sumstats/munger.py>); [src/ldsc/_regression/runner.py:RegressionDataset.validate](<../../../../src/ldsc/_regression/runner.py>).

Legacy implementation context: [`ldscore/regressions.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/regressions.py), [`ldscore/sumstats.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/sumstats.py).

<a id="issue-242"></a>

### #242: ValueError: Columns ['P'] are expected to be numeric

[Original issue](https://github.com/bulik/ldsc/issues/242) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **1**.

Workflow tags: `munge-sumstats`.

**Trigger and harm.** The munger says P must be numeric; a commenter links extreme p-values/odds ratios.
Strong association signals can be lost or preprocessing blocked by floating-point representation.

**Legacy disposition.** User advice identifies a possible cause, but the original data and a merged fix are absent.

**LDSC3 assessment — Partial coverage.** LDSC3 improves numeric-column diagnostics but raw P is still parsed as ordinary floating point and converted with chi-square inverse survival.
Values smaller than double precision can become0 and be removed; an input Z does not rescue magnitude because P is authoritative.

**Recommended prevention.** Explain supported P precision and underflow removals; consider a deliberate log-P or direct-Z route with stable numerics rather than arbitrary p-value flooring.

**Evidence and verification.** source-inspected
[src/ldsc/_sumstats/munge_input.py:prepare_munge_input](<../../../../src/ldsc/_sumstats/munge_input.py>); [src/ldsc/_sumstats/munge_kernel.py:munge_sumstats](<../../../../src/ldsc/_sumstats/munge_kernel.py>).

Legacy implementation context: [`munge_sumstats.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/munge_sumstats.py).

<a id="issue-258"></a>

### #258: MemoryError

[Original issue](https://github.com/bulik/ldsc/issues/258) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **1**.

Workflow tags: `ldscore`.

**Trigger and harm.** A64GB BED reference panel exhausts230GB RAM.
Dense genotype/block expansion can make reference LD-score generation impractical.

**Legacy disposition.** Open; another user reports the same issue, with no verified fix.

**LDSC3 assessment — Partial coverage.** LDSC3 implements streamed PLINK access, explicit LD windows, query batches and memory ownership controls.
These reduce avoidable allocations but do not prove a64GB panel fits in230GB: genotype sample count, SNP count, LD-window density and concurrent workers determine peak memory.

**Recommended prevention.** Document a measured memory estimate and small-chromosome pilot; use representative reference samples and conservative workers before complex distributed designs.

**Evidence and verification.** source-inspected
[src/ldsc/_refpanel/plink.py](<../../../../src/ldsc/_refpanel/plink.py>); [src/ldsc/_refpanel/build_kernel.py](<../../../../src/ldsc/_refpanel/build_kernel.py>); [docs/dev/design/memory-and-cpu.md](<../../design/memory-and-cpu.md>).

Legacy implementation context: [`ldscore/ldscore.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/ldscore.py).

<a id="issue-263"></a>

### #263: Memory Error for Linkage disequilibrium with annotation

[Original issue](https://github.com/bulik/ldsc/issues/263) · Source state: **closed** · Severity: **high** · Accessible comments reviewed: **1**.

Workflow tags: `ldscore`.

**Trigger and harm.** Legacy LD calculation allocates a b-by500000 correlation block and exhausts300GB.
An oversized compute-block parameter causes excessive allocation; restricting reference SNPs can change the scientific denominator.

**Legacy disposition.** Reporter reduces panel SNPs to summary-statistics SNPs and closes; that is not a generally valid LD-reference fix.

**LDSC3 assessment — Partial coverage.** LDSC3 separates LD compute blocks/query batches from munging chunks and preserves broad reference SNPs versus filtered regression SNPs.
Its explicit windows and streamed backend address the allocation mechanism, but users must not shrink the LD-contributing panel solely to force memory usage down.

**Recommended prevention.** Explain block size units and memory scaling; reduce compute blocks/workers rather than silently redefining reference SNPs.

**Evidence and verification.** source-inspected
[src/ldsc/_refpanel/build_kernel.py](<../../../../src/ldsc/_refpanel/build_kernel.py>); [docs/reference/ldscore-snp-universes.md](<../../../reference/ldscore-snp-universes.md>); [docs/dev/design/memory-and-cpu.md](<../../design/memory-and-cpu.md>).

Legacy implementation context: [`ldscore/ldscore.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/ldscore.py).

<a id="issue-268"></a>

### #268: It took long time to run munge_sumstats.py --

[Original issue](https://github.com/bulik/ldsc/issues/268) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **3**.

Workflow tags: `munge-sumstats`.

**Trigger and harm.** A one-million-SNP example takes8hours and six million SNPs run72hours without output.
The first preprocessing step becomes unusable and invites unnecessary parallelization.

**Legacy disposition.** Two comments report smaller --chunksize500000 helps; no algorithmic fix established.

**LDSC3 assessment — Partial coverage.** LDSC3 has a lower one-million-row default, progress reporting and configurable chunks.
It still accumulates all retained chunks for global operations, so runtime/memory scaling remains a real validation concern.
Parallelizing input parsing would not automatically solve allocation pressure and may worsen peak memory.

**Recommended prevention.** Offer the established smaller-chunk remedy and log row progress; measure before introducing parallel or out-of-core machinery.

**Evidence and verification.** source-inspected
[src/ldsc/_sumstats/munge_kernel.py:parse_dat](<../../../../src/ldsc/_sumstats/munge_kernel.py>); [src/ldsc/_sumstats/munger.py:1417](<../../../../src/ldsc/_sumstats/munger.py#L1417>).

Legacy implementation context: [`munge_sumstats.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/munge_sumstats.py).

<a id="issue-275"></a>

### #275: munge_sumstats error on encountering low p-values. Implications for genetic correlation?

[Original issue](https://github.com/bulik/ldsc/issues/275) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **0**.

Workflow tags: `munge-sumstats`, `rg`.

**Trigger and harm.** Raw P reaches1e-2134 and preprocessing works only after flooring to1e-300.
Clipping changes extreme chi-square magnitudes and potentially regression weights or genetic-correlation estimates.

**Legacy disposition.** No fix supplied; the user reports a workaround and asks about its consequences.

**LDSC3 assessment — Partial coverage.** LDSC3 still treats raw P as floating point, requires P even when Z is supplied and derives magnitude through chi-square isf.
Values outside double precision cannot be represented faithfully and may become0 or object dtype.
Its improved numeric error does not preserve the original signal.

**Recommended prevention.** Provide a stable deliberate direct-Z or log-P input route before promising such support; explain that arbitrary flooring changes analysis and test consequences.

**Evidence and verification.** source-inspected
[src/ldsc/_sumstats/munge_input.py:prepare_munge_input](<../../../../src/ldsc/_sumstats/munge_input.py>); [src/ldsc/_sumstats/munge_kernel.py:munge_sumstats](<../../../../src/ldsc/_sumstats/munge_kernel.py>).

Legacy implementation context: [`munge_sumstats.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/munge_sumstats.py), [`ldscore/regressions.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/regressions.py), [`ldscore/sumstats.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/sumstats.py).

<a id="issue-293"></a>

### #293: LD score regression of data including SNPs with extremely low P-value

[Original issue](https://github.com/bulik/ldsc/issues/293) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **0**.

Workflow tags: `munge-sumstats`, `h2`, `rg`.

**Trigger and harm.** User considers dropping highly significant SNPs or rounding their p-values.
Both workarounds alter signal and can bias SNP-heritability or covariance fitting.

**Legacy disposition.** Open without a validated remedy.

**LDSC3 assessment — Partial coverage.** LDSC3 raw munging still requires floating-point P and derives Z magnitude from it, so values beyond double precision cannot be represented faithfully.
The existing direct legacy Z,N route can preserve already valid finite Z statistics, but manual reconstruction requires valid sign and N and does not justify arbitrary clipping.

**Recommended prevention.** Document precision limits and a scientifically valid direct-Z workflow; consider stable log-P support only with appropriate numerical tests.

**Evidence and verification.** source-inspected
[src/ldsc/_sumstats/munge_input.py:prepare_munge_input](<../../../../src/ldsc/_sumstats/munge_input.py>); [src/ldsc/_sumstats/munge_kernel.py:munge_sumstats](<../../../../src/ldsc/_sumstats/munge_kernel.py>).

Legacy implementation context: [`munge_sumstats.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/munge_sumstats.py), [`ldscore/regressions.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/regressions.py), [`ldscore/sumstats.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/sumstats.py).

<a id="issue-312"></a>

### #312: Duration of execution for munge_sumstats.py

[Original issue](https://github.com/bulik/ldsc/issues/312) · Source state: **closed** · Severity: **high** · Accessible comments reviewed: **1**.

Workflow tags: `munge-sumstats`.

**Trigger and harm.** A729MB file runs24hours with no log progress on32GB RAM. A common onboarding step appears hung despite modest file size.

**Legacy disposition.** Reporter confirms --chunksize500000 finishes in about one minute; issue closes without a code fix.

**LDSC3 assessment — Partial coverage.** LDSC3 lowers the default parse chunk size and adds phase/progress logs, but still retains surviving chunks and concatenates globally.
The dramatic legacy behavior is not a universal computational cost; supported-stack measurements are needed before promising a bound or introducing parallel parsing.

**Recommended prevention.** Document the observed practical chunk workaround and a measurable progress expectation; keep parsing size separate from LD compute blocks.

**Evidence and verification.** source-inspected
[src/ldsc/_sumstats/munge_kernel.py:parse_dat](<../../../../src/ldsc/_sumstats/munge_kernel.py>); [src/ldsc/_sumstats/munger.py:1417](<../../../../src/ldsc/_sumstats/munger.py#L1417>).

Legacy implementation context: [`munge_sumstats.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/munge_sumstats.py).

<a id="issue-323"></a>

### #323: Using stratified model for Rg (old_weight is always FALSE for --rg)

[Original issue](https://github.com/bulik/ldsc/issues/323) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **0**.

Workflow tags: `h2`, `rg`.

**Trigger and harm.** Same variants and same x,y,w,N,M yield negative RG heritability because old_weights is False versus True standalone.
The same Hsq inputs can be fitted by fixed-initial-weight WLS or iterative reweighting, potentially changing the component estimates.

**Legacy disposition.** Open without a fixing commit; reporter localizes a concrete legacy policy discrepancy.

**LDSC3 assessment — Partial coverage.** LDSC3 standalone multi-annotation heritability uses old_weights=True (fixed initial-weight WLS), while RG constructs Hsq with its default False (iterative reweighting).
This source-confirmed policy difference is separate from pairwise SNP intersection.
The original dataset and reported sign reversal were not reproduced, and the difference alone does not establish which estimator policy is scientifically preferable.

**Recommended prevention.** Define and validate the supported multi-annotation genetic-correlation weighting policy; use a matched-input regression test and document intended differences before changing it.

**Evidence and verification.** source-inspected
[src/ldsc/_regression/runner.py:879](<../../../../src/ldsc/_regression/runner.py#L879>); [src/ldsc/_regression/estimators.py:LD_Score_Regression.__init__](<../../../../src/ldsc/_regression/estimators.py>); [src/ldsc/_regression/estimators.py:RG](<../../../../src/ldsc/_regression/estimators.py>).

Legacy implementation context: [`ldscore/regressions.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/regressions.py), [`ldscore/sumstats.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/sumstats.py).

<a id="issue-418"></a>

### #418: ERROR: P&lt;1e-300

[Original issue](https://github.com/bulik/ldsc/issues/418) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **3**.

Workflow tags: `munge-sumstats`.

**Trigger and harm.** Very small p-values such as 1e-400 become floating-point zero during raw munging.
Strong association rows can be dropped or lose usable magnitudes despite finite original Z statistics.

**Legacy disposition.** A commenter recommends finite original Z in a ready summary-statistics file; no upstream stable-log-P fix is established.

**LDSC3 assessment — Partial coverage.** LDSC3 still filters P to (0,1] and derives magnitude through chi2.isf.
Adding signed Z to the same raw input does not avoid P reconstruction.
Direct validated munged Z input is a distinct route, not a repair of the munger.

**Recommended prevention.** Prefer supplied finite signed Z or beta/SE through an explicit validated magnitude path; alternatively support log-P if truly needed.
Do not silently clamp tiny P or describe 1e-400 as mathematically zero.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[src/ldsc/_sumstats/munge_kernel.py:185](<../../../../src/ldsc/_sumstats/munge_kernel.py#L185>); [src/ldsc/_sumstats/munge_kernel.py:652](<../../../../src/ldsc/_sumstats/munge_kernel.py#L652>); [src/ldsc/_sumstats/munge_kernel.py:741](<../../../../src/ldsc/_sumstats/munge_kernel.py#L741>); [docs/reference/formats/legacy-sumstats.md](<../../../reference/formats/legacy-sumstats.md>).

Legacy implementation context: [`munge_sumstats.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/munge_sumstats.py).

<a id="issue-496"></a>

### #496: LDSC Genetic Correlation Analysis Returns All NA Values - Seeking Help

[Original issue](https://github.com/bulik/ldsc/issues/496) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **0**.

Workflow tags: `h2`, `rg`.

**Trigger and harm.** A rare binary trait has almost null inflation and rg returns all NA despite successful univariate fits.
Successful h2 commands do not establish positive stable delete-one denominators; opaque NA output obscures why inference failed.

**Legacy disposition.** Open; no full failed-pair traceback or verified historical fix is supplied.

**LDSC3 assessment — Partial coverage.** LDSC3 still has a confirmed boundary where positive whole-data h2 coexists with negative delete products, yielding rg_se/p NaN while a pair is marked ok.
The root probe demonstrates this failure family, not the cause of this exact unavailable input.

**Recommended prevention.** Validate every delete-one heritability denominator and report an explicit unestimable pair with reason while preserving finite covariance/intercept estimates; do not tune intercept/cutoffs merely to force a number.

**Evidence and verification.** Root V2 deterministic local probe reproduces nonfinite rg uncertainty with status ok; historical GWAS inputs not rerun.
[src/ldsc/_regression/estimators.py:772](<../../../../src/ldsc/_regression/estimators.py#L772>); [src/ldsc/_regression/runner.py:_fit_rg_dataset](<../../../../src/ldsc/_regression/runner.py>); [src/ldsc/_regression/runner.py:_summarize_rg_pair](<../../../../src/ldsc/_regression/runner.py>).

Legacy implementation context: [`ldscore/regressions.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/regressions.py), [`ldscore/sumstats.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/sumstats.py).
[Shared reproduced boundary V2](verification.md#v2-genetic-correlation-delete-block-boundary); the original historical dataset was not reproduced.

<a id="faq-09"></a>

### FAQ-09: Collinear annotation designs

[FAQ source](https://github.com/bulik/ldsc/wiki/FAQ) · Source-order identifier: **FAQ-09** · Severity: **high**.

Workflow tags: `baseline-partitioned-h2`, `ldsc-seg`, `continuous`.

**Legacy disposition.** Historical FAQ guidance; assessed independently rather than treated as a code-fix record.

**LDSC3 assessment — Partial coverage.** Current model assembly checks near-collinearity before fitting, and singular delete-block fits receive detailed diagnostics. A full-data design that passes still need not remain identifiable after every block deletion.

**Recommended prevention.** Retain explicit model failure with offending annotations/blocks; simplify scientifically redundant annotations rather than introduce a silent pseudoinverse.

**Evidence and verification.** Source-inspected; relevant executed checks are identified in verification.md.
[src/ldsc/_regression/overlap_model.py:79](<../../../../src/ldsc/_regression/overlap_model.py#L79>); [src/ldsc/_regression/runner.py:2733](<../../../../src/ldsc/_regression/runner.py#L2733>); [tests/regression/test_partitioned_query_failures.py](<../../../../tests/regression/test_partitioned_query_failures.py>).

<a id="issue-22"></a>

### #22: fix prop in dev

[Original issue](https://github.com/bulik/ldsc/issues/22) · Source state: **closed** · Severity: **high** · Accessible comments reviewed: **0**.

Workflow tags: `baseline-partitioned-h2`.

**Trigger and harm.** Developer flags disabled proportion tests in the dev branch.
An incorrect proportion calculation changes partitioned interpretation and enrichment.

**Legacy disposition.** Closed and referenced 60e633dc7377c1cbee86844a64120ec549a8be19; the short report does not specify the numerical discrepancy, so that event is not a verified fix.

**LDSC3 assessment — Indeterminate.** LDSC3 implements overlap-aware proportion estimates in Hsq and overlap_model with explicit count/column alignment.
No failing input, expected value or identified legacy hunk is supplied here, so the exact historical bug cannot be declared reproduced or cured.
Existing overlap tests support current behavior but are not a retrospective proof of this unspecified failure.

**Recommended prevention.** Keep the issue indeterminate and link hand-worked overlap calculations; avoid inferring resolution from closure or passing unrelated tests.

**Evidence and verification.** source-inspected
[src/ldsc/_regression/estimators.py:Hsq._overlap_output](<../../../../src/ldsc/_regression/estimators.py>); [src/ldsc/_regression/overlap_model.py](<../../../../src/ldsc/_regression/overlap_model.py>); [tests/regression/test_overlap_matrix.py](<../../../../tests/regression/test_overlap_matrix.py>).

Legacy implementation context: [`ldscore/regressions.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/regressions.py).

<a id="issue-292"></a>

### #292: munge_sumstats.py takes forever with file contains more than 1M lines.

[Original issue](https://github.com/bulik/ldsc/issues/292) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **3**.

Workflow tags: `munge-sumstats`.

**Trigger and harm.** Files differing by a few rows change runtime from20seconds to many hours.
Parsing/environment pathology blocks normal preprocessing and is not explained by row count alone.

**Legacy disposition.** One suggestion is --chunksize500000, but another user says it does not help; no proven root cause/fix.

**LDSC3 assessment — Indeterminate.** LDSC3 lowers default chunks and uses a modern stack, but no original file/environment is supplied to reproduce this discontinuity.
It remains wrong to claim every slow-munging issue was fixed solely by a new chunk default; retained table accumulation still consumes memory.

**Recommended prevention.** Offer the practical chunk workaround and progress evidence, then request a minimal boundary reproduction before any algorithm redesign.

**Evidence and verification.** source-inspected
[src/ldsc/_sumstats/munge_kernel.py:parse_dat](<../../../../src/ldsc/_sumstats/munge_kernel.py>); [src/ldsc/_sumstats/munger.py:1417](<../../../../src/ldsc/_sumstats/munger.py#L1417>).

Legacy implementation context: [`munge_sumstats.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/munge_sumstats.py).

<a id="issue-10"></a>

### #10: rg regression weights

[Original issue](https://github.com/bulik/ldsc/issues/10) · Source state: **closed** · Severity: **high** · Accessible comments reviewed: **1**.

Workflow tags: `rg`.

**Trigger and harm.** The product-statistic variance incorrectly doubles the covariance contribution.
Incorrect weights can reduce efficiency and alter genetic-correlation uncertainty.

**Legacy disposition.** Maintainer explicitly reported the factor error fixed in dev for v1.0.0; no exact fixing commit was verified from the discussion.
Current and historical weight formulas were inspected separately.

**LDSC3 assessment — Apparently prevented.** LDSC3 Gencov.weights uses the product of marginal expected chi-square terms plus one squared covariance term, matching the corrected product-of-normal-variables variance.
The method is inherited rather than rederived by the CLI.
This is source evidence of the corrected formula, not an independent numerical calibration of every weight update.

**Recommended prevention.** Retain a small known-value variance-formula test and cite the corrected weight expression when documenting inherited numerical behavior.

**Evidence and verification.** source-inspected
[src/ldsc/_regression/estimators.py:Gencov.weights](<../../../../src/ldsc/_regression/estimators.py>); [tests/regression/test_kernel_regression.py](<../../../../tests/regression/test_kernel_regression.py>).

Legacy implementation context: [`ldscore/regressions.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/regressions.py), [`ldscore/sumstats.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/sumstats.py).

<a id="issue-11"></a>

### #11: --overlap-annot and _keep_ld/_check_variance

[Original issue](https://github.com/bulik/ldsc/issues/11) · Source state: **closed** · Severity: **high** · Accessible comments reviewed: **2**.

Workflow tags: `baseline-partitioned-h2`.

**Trigger and harm.** LD-score columns are removed while annotation overlap matrices or counts retain another order.
Heritability can be assigned to the wrong annotation without a shape error.

**Legacy disposition.** Maintainer cites e2761618996529fcc084dbae865178f5888cb559; local history verifies changes in regressions.py and sumstats.py, but the large patch is not reduced to a single proven hunk here.

**LDSC3 assessment — Apparently prevented.** LDSC3 assembles per-model counts and overlaps by the selected annotation names and checks names/order.
Duplicate names are rejected before generation.
Zero-variance LD columns are removed after the applicable SNP merge, so the audit should preserve explicit checks that every corresponding count/overlap selection follows the retained list.

**Recommended prevention.** Retain name/order invariants and include a pruned-middle-column case in future parity validation.

**Evidence and verification.** source-inspected
[src/ldsc/_regression/overlap_model.py](<../../../../src/ldsc/_regression/overlap_model.py>); [src/ldsc/_regression/runner.py:build_dataset](<../../../../src/ldsc/_regression/runner.py>); [tests/regression/test_overlap_matrix.py](<../../../../tests/regression/test_overlap_matrix.py>); [Source](<https://github.com/bulik/ldsc/commit/e2761618996529fcc084dbae865178f5888cb559>).

Legacy implementation context: [`ldscore/regressions.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/regressions.py).

<a id="issue-17"></a>

### #17: dev currently has no chi^2 filter

[Original issue](https://github.com/bulik/ldsc/issues/17) · Source state: **closed** · Severity: **high** · Accessible comments reviewed: **0**.

Workflow tags: `h2`, `rg`.

**Trigger and harm.** Early dev code lacks the chi-square filtering policy required by its estimator.
Large associations can dominate intercept fitting or make analyses inconsistent across commands.

**Legacy disposition.** Verified local commit 4f2ae1704fcc1bc7a226b4bfca2e7a3562b7da97 explicitly closes #17 and changes estimator/filtering code.

**LDSC3 assessment — Apparently prevented.** LDSC3 defaults to two-step estimation for single-annotation free-intercept h2/rg and an explicit chi-square cap for multivariate partitioned models.
The defaults intentionally differ by estimator and fixed-intercept setting.
Consequently absence of a universal filter is not the original bug; the audit must verify each policy on the final fitted SNP set.

**Recommended prevention.** Explain estimator-specific defaults and preserve focused tests for fixed intercepts, inclusive cutoffs and postfilter counts.

**Evidence and verification.** source-inspected
[src/ldsc/_regression/runner.py:_fit_h2_dataset](<../../../../src/ldsc/_regression/runner.py>); [src/ldsc/_regression/runner.py:_fit_rg_dataset](<../../../../src/ldsc/_regression/runner.py>); [tests/regression/test_regression_fit_outcomes.py](<../../../../tests/regression/test_regression_fit_outcomes.py>); [Source](<https://github.com/bulik/ldsc/commit/4f2ae1704fcc1bc7a226b4bfca2e7a3562b7da97>).

Legacy implementation context: [`ldscore/regressions.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/regressions.py), [`ldscore/sumstats.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/sumstats.py).

<a id="issue-24"></a>

### #24: deal with mismatched alleles in sumstats files not generated with --merge-alleles

[Original issue](https://github.com/bulik/ldsc/issues/24) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **2**.

Workflow tags: `rg`.

**Trigger and harm.** Two valid SNP allele pairs are mutually incompatible after the trait merge. A raw dictionary KeyError hides which SNPs were excluded or why the pair failed.

**Legacy disposition.** Open; maintainer recommends legacy --merge-alleles as prevention, without a verified numerical fix in the issue.

**LDSC3 assessment — Apparently prevented.** LDSC3 validates/harmonizes the merged allele pairs and removes incompatible SNPs before sign alignment; legacy projection reports strand-ambiguous and incompatible mappings.
It therefore does not rely on a user remembering an old munging flag.
An emptied intersection still must yield an actionable failure with the selected identity settings.

**Recommended prevention.** Preserve explicit allele-drop reasons and include both incompatible and strand-complement cases in the issue-to-test map.

**Evidence and verification.** source-inspected
[src/ldsc/_regression/runner.py:build_rg_dataset](<../../../../src/ldsc/_regression/runner.py>); [src/ldsc/_regression/estimators.py:_filter_alleles](<../../../../src/ldsc/_regression/estimators.py>); [tests/regression/test_regression_workflow.py:test_build_rg_dataset_harmonizes_alleles_and_flips_second_trait_z](<../../../../tests/regression/test_regression_workflow.py>).

Legacy implementation context: [`ldscore/regressions.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/regressions.py), [`ldscore/sumstats.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/sumstats.py).

<a id="issue-35"></a>

### #35: no-intercept produces h2=4.84...

[Original issue](https://github.com/bulik/ldsc/issues/35) · Source state: **closed** · Severity: **high** · Accessible comments reviewed: **8**.

Workflow tags: `h2`, `baseline-partitioned-h2`.

**Trigger and harm.** A historical partitioned fit jumps from 0.17 to 4.84 under no-intercept.
An intercept implementation error can silently yield plausible-looking but invalid scientific output.

**Legacy disposition.** The user confirms updating legacy LDSC changed the constrained estimate to 0.2545; referenced 860c95a8b9726a52d7c04844ca7548345e0ee0c8 is unavailable locally and was not verified as the fix.

**LDSC3 assessment — Apparently prevented.** LDSC3 uses explicit --intercept-h2 1 rather than the removed shortcut.
Focused tests cover fixed/free fits and incompatible two-step combinations.
Fixed estimates can still differ substantively from free-intercept estimates under misspecification; matching the old number alone is insufficient proof of correctness.

**Recommended prevention.** Keep known-value fixed-intercept checks and explain the statistical assumption; do not recommend fixing the intercept simply to force heritability into range.

**Evidence and verification.** source-inspected
[src/ldsc/_regression/runner.py:_select_intercept](<../../../../src/ldsc/_regression/runner.py>); [tests/regression/test_regression_fit_outcomes.py:test_explicit_cli_intercepts_preserve_standard_fixed_fit](<../../../../tests/regression/test_regression_fit_outcomes.py>); [docs/reference/legacy-flag-map.md](<../../../reference/legacy-flag-map.md>).

Legacy implementation context: [`ldscore/regressions.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/regressions.py).

<a id="issue-72"></a>

### #72: smart_merge crashes when x.index is different from y.index

[Original issue](https://github.com/bulik/ldsc/issues/72) · Source state: **closed** · Severity: **high** · Accessible comments reviewed: **0**.

Workflow tags: `h2`, `rg`.

**Trigger and harm.** Two SNP tables have equal lengths and matching values but different pandas row labels.
The optimization crashes or risks mistaken row correspondence.

**Legacy disposition.** Verified patch 94b3957e77b4d723b146024eb1abb77dfec914b5, titled Fix #72 and merged through PR #73, adds an x.index==y.index guard before the fast SNP-Series comparison.

**LDSC3 assessment — Apparently prevented.** LDSC3 assembles datasets by explicit effective SNP identity and genomic order instead of the historical smart_merge shortcut.
Tests cover reordering, coordinate-family joins and allele harmonization.
It must still distinguish positional row equality from identity equality; merely resetting indices without matching keys would be unsafe.

**Recommended prevention.** Retain a differently indexed/permuted-input test and validate one-to-one keys before any positional fast path.

**Evidence and verification.** source-inspected
[src/ldsc/_regression/runner.py:build_dataset](<../../../../src/ldsc/_regression/runner.py>); [src/ldsc/_core/row_alignment.py](<../../../../src/ldsc/_core/row_alignment.py>); [tests/regression/test_regression_workflow.py:test_build_dataset_orders_merged_rows_by_ld_score_genomic_order](<../../../../tests/regression/test_regression_workflow.py>); [Source](<https://github.com/bulik/ldsc/commit/94b3957e77b4d723b146024eb1abb77dfec914b5>).

Legacy implementation context: [`ldscore/regressions.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/regressions.py), [`ldscore/sumstats.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/sumstats.py).

<a id="issue-77"></a>

### #77: Error while partitioning heritability with the baseline model

[Original issue](https://github.com/bulik/ldsc/issues/77) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **1**.

Workflow tags: `baseline-partitioned-h2`.

**Trigger and harm.** A filtered baseline supplies 53 overlap values for a one-column fitted model.
Broadcast errors can hide inconsistent annotation selection and counts.

**Legacy disposition.** Maintainer asks for files and points to empty/extra columns or mismatched M; no verified legacy fix.

**LDSC3 assessment — Apparently prevented.** LDSC3 saves named annotations, overlap rows and count metadata together and selects complete per-model overlap blocks.
Source validation rejects missing requested counts and duplicate names rather than broadcasting unnamed arrays.
The actual legacy input is unavailable, so only the failure class is assessed, not its specific root cause.

**Recommended prevention.** Keep structural dimension/name validation before fitting and report every mismatched annotation/count source in the diagnostic.

**Evidence and verification.** source-inspected
[src/ldsc/_regression/overlap_model.py](<../../../../src/ldsc/_regression/overlap_model.py>); [src/ldsc/_ldscore/source.py](<../../../../src/ldsc/_ldscore/source.py>); [tests/regression/test_overlap_matrix.py](<../../../../tests/regression/test_overlap_matrix.py>).

Legacy implementation context: [`ldscore/regressions.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/regressions.py).

<a id="issue-84"></a>

### #84: Question about partition heritability stratified by allele freqency

[Original issue](https://github.com/bulik/ldsc/issues/84) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **1**.

Workflow tags: `ldscore`, `baseline-partitioned-h2`.

**Trigger and harm.** The LD-score kernel slices annotation matrices with a float-valued block boundary.
Ordinary and annotation-weighted LD-score calculations fail.

**Legacy disposition.** Maintainer suspects c42e055c68e1be9a9c4bb27c7e5b544826b44f43 but explicitly could not reproduce; that is a candidate fix, not verified user resolution.

**LDSC3 assessment — Apparently prevented.** LDSC3 has rewritten exact-window projection and explicit integer index geometry, with PLINK and R² backend tests.
MAF-bin annotations are supplied numeric columns; their scientific interpretation still depends on common/all SNP counts and appropriate baseline controls.
No audit evidence establishes the original tmp panel root cause.

**Recommended prevention.** Retain integer-boundary tests and compare exact windows on a small fixture; document MAF-bin denominators instead of adding a special slicing workaround.

**Evidence and verification.** source-inspected
[src/ldsc/_ldscore/compute.py](<../../../../src/ldsc/_ldscore/compute.py>); [src/ldsc/_refpanel/windows.py](<../../../../src/ldsc/_refpanel/windows.py>); [tests/ldscore/test_plink_exact_windows.py](<../../../../tests/ldscore/test_plink_exact_windows.py>).

Legacy implementation context: [`ldscore/ldscore.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/ldscore.py), [`ldscore/regressions.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/regressions.py).

<a id="issue-85"></a>

### #85: Unalignable boolean Series provided as indexer

[Original issue](https://github.com/bulik/ldsc/issues/85) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **2**.

Workflow tags: `baseline-partitioned-h2`.

**Trigger and harm.** Annotation and frequency files contain different variants or pandas indices.
Row-aligned masking can crash or select incorrect reference SNPs.

**Legacy disposition.** The user confirms mismatched frequency files caused the error; no algorithmic patch was required for that input.

**LDSC3 assessment — Apparently prevented.** Native LDSC3 derives MAF from the selected reference panel and aligns annotations by explicit identities.
Legacy conversion validates/reconstructs overlaps and counts according to its contract rather than blindly applying another table index.
Matched lengths alone remain insufficient evidence of matched SNPs or release.

**Recommended prevention.** Require identity-aware joins and release provenance; make mismatch diagnostics identify differing SNPs and the affected annotation family.

**Evidence and verification.** source-inspected
[src/ldsc/_refpanel/adapters.py](<../../../../src/ldsc/_refpanel/adapters.py>); [src/ldsc/_ldscore/legacy_converter.py](<../../../../src/ldsc/_ldscore/legacy_converter.py>); [docs/reference/formats/legacy-ldscores.md](<../../../reference/formats/legacy-ldscores.md>).

Legacy implementation context: [`ldscore/regressions.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/regressions.py).

<a id="issue-97"></a>

### #97: Partitioned heritability inflated SE and FloatingPointError

[Original issue](https://github.com/bulik/ldsc/issues/97) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **0**.

Workflow tags: `baseline-partitioned-h2`.

**Trigger and harm.** A model includes base plus all 22 chromosome indicators, possibly together with exhaustive MAF bins.
Exact linear dependence makes conditional coefficients nonidentifiable and can produce enormous SEs or negative numerical variances.

**Legacy disposition.** Open with no maintainer resolution; the displayed base-plus-all-chromosomes design is structurally redundant if correctly generated.

**LDSC3 assessment — Apparently prevented.** LDSC3 checks model collinearity and reports singular delete-block models instead of using pseudoinverse.
Removing a redundant column or changing the scientific partition can produce an estimable model; merely splitting output files or downgrading libraries does not solve the mathematical dependence.
The original arrays were not supplied, so this remains a design-based diagnosis.

**Recommended prevention.** Explain exhaustive-partition redundancy and show one identifiable parameterization; never force inversion or suppress invalid-square-root errors.

**Evidence and verification.** source-inspected
[src/ldsc/_regression/overlap_model.py:model_collinearity_error](<../../../../src/ldsc/_regression/overlap_model.py>); [src/ldsc/_regression/runner.py:_raise_on_model_collinearity](<../../../../src/ldsc/_regression/runner.py>); [src/ldsc/_regression/jackknife.py:LstsqJackknifeFast](<../../../../src/ldsc/_regression/jackknife.py>).

Legacy implementation context: [`ldscore/regressions.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/regressions.py).

<a id="issue-99"></a>

### #99: Discrepancies in output (results table) for Liability Scale Heritability

[Original issue](https://github.com/bulik/ldsc/issues/99) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **2**.

Workflow tags: `rg`.

**Trigger and harm.** The first pair is correct but later summary rows use another trait prevalence.
Reported liability-scale heritability and SE can be wrong while log values are correct.

**Legacy disposition.** Maintainer recognizes a recurring prevalence-index bug and says earlier log values are correct; no verified fixing commit appears.

**LDSC3 assessment — Apparently prevented.** LDSC3 resolves per-trait prevalence before fitting and passes prev[i] and prev[j] explicitly to pair summaries.
Separate standalone per-trait results use the corresponding table/prevalence tuple.
Existing prevalence tests cover list/manifest validation, but the exact three-distinct-prevalence row-order scenario should remain a release-level check.

**Recommended prevention.** Validate at least three binary traits with distinct prevalences and reversed order against hand conversions; label observed and liability columns consistently.

**Evidence and verification.** source-inspected
[src/ldsc/_regression/runner.py:estimate_rg_pairs](<../../../../src/ldsc/_regression/runner.py>); [src/ldsc/_regression/prevalence.py](<../../../../src/ldsc/_regression/prevalence.py>); [tests/regression/test_prevalence.py](<../../../../tests/regression/test_prevalence.py>).

Legacy implementation context: [`ldscore/regressions.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/regressions.py), [`ldscore/sumstats.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/sumstats.py).

<a id="issue-103"></a>

### #103: Run Partitioned heritability have errors

[Original issue](https://github.com/bulik/ldsc/issues/103) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **2**.

Workflow tags: `baseline-partitioned-h2`.

**Trigger and harm.** Frequency-mask indices do not match baseline annotation SNPs.
The run fails after an already suspicious negative heritability fit, inviting wrong fixes.

**Legacy disposition.** Maintainer identifies the release mismatch and user confirms matching Phase 3 frequency files solve the parsing error.

**LDSC3 assessment — Apparently prevented.** Native LDSC3 derives MAF from its panel and records overlaps/counts together, removing the independent frqfile mask path.
Conversion validates legacy families.
Resolving the file mismatch does not itself validate the prior negative h2 or unusually high intercept; those require separate scientific QC on matched inputs.

**Recommended prevention.** Make release matching explicit and rerun scientific diagnostics after correcting resources; do not treat successful parsing as validity of the inference.

**Evidence and verification.** source-inspected
[src/ldsc/_refpanel/adapters.py](<../../../../src/ldsc/_refpanel/adapters.py>); [src/ldsc/_ldscore/legacy_converter.py](<../../../../src/ldsc/_ldscore/legacy_converter.py>); [docs/reference/formats/legacy-ldscores.md](<../../../reference/formats/legacy-ldscores.md>).

Legacy implementation context: [`ldscore/regressions.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/regressions.py).

<a id="issue-114"></a>

### #114: Incorrect SNP filtering and MAF when using --keep in LD score calculations

[Original issue](https://github.com/bulik/ldsc/issues/114) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **0**.

Workflow tags: `ldscore`.

**Trigger and harm.** Keeping a sample count not divisible by four leaves uninitialized padded genotype bits in the historical buffer.
MAF and SNP filtering can depend on memory contents, including impossible frequencies.

**Legacy disposition.** Issue gives a root cause and candidate one-line fix 60efb37, unavailable in local history; do not claim a merged repair from the proposal alone.

**LDSC3 assessment — Apparently prevented.** LDSC3 selected-SNP reader explicitly allocates the filtered bitarray and calls z.setall(0) before copying individuals.
Thus padding remains zero and cannot contaminate bit-count MAF.
Existing golden/selective-read checks support the rewrite; a non-multiple-of-four retained sample probe is more specific than four-sample parity.

**Recommended prevention.** Keep a 1/2/3-mod-four sample-count fixture checking both MAF and retained SNPs against decoded genotypes; advise recomputation only for genuinely affected historical scores.

**Evidence and verification.** source-inspected
[src/ldsc/_refpanel/plink.py:189-207](<../../../../src/ldsc/_refpanel/plink.py#L189>); [tests/refpanel/test_plink_bed_reader.py](<../../../../tests/refpanel/test_plink_bed_reader.py>); [tests/refpanel/test_plink_io.py](<../../../../tests/refpanel/test_plink_io.py>).

Legacy implementation context: [`ldscore/ldscore.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/ldscore.py).

<a id="issue-126"></a>

### #126: make_annot.py returns different number of results on the same bim file with different annotations

[Original issue](https://github.com/bulik/ldsc/issues/126) · Source state: **closed** · Severity: **high** · Accessible comments reviewed: **1**.

Workflow tags: `ldscore`.

**Trigger and harm.** BED intersection with duplicated genomic positions yields more rows than the BIM grid.
Without SNP identifiers, expanded rows cannot be aligned reliably to genotypes.

**Legacy disposition.** User identifies duplicate BIM BP as the cause; no source fix is established.

**LDSC3 assessment — Apparently prevented.** LDSC3 interval projection returns one value per retained grid row and identity cleanup applies the selected allele-aware/base policy.
It writes identity metadata rather than anonymous thin columns.
Duplicate positions with different alleles are not automatically the same SNP, so simply deduplicating BP would be scientifically unsafe.

**Recommended prevention.** Keep one-row-per-retained-identity tests with duplicate coordinates and distinct alleles, and report dropped duplicate identities explicitly.

**Evidence and verification.** source-inspected
[src/ldsc/_annotation/interval_projection.py:ChromosomeProjector.project](<../../../../src/ldsc/_annotation/interval_projection.py>); [src/ldsc/_annotation/bim_grid.py](<../../../../src/ldsc/_annotation/bim_grid.py>); [docs/reference/annotate.md](<../../../reference/annotate.md>).

Legacy implementation context: [`make_annot.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/make_annot.py).

<a id="issue-137"></a>

### #137: Running error --rg Incompatible alleles

[Original issue](https://github.com/bulik/ldsc/issues/137) · Source state: **closed** · Severity: **high** · Accessible comments reviewed: **1**.

Workflow tags: `rg`.

**Trigger and harm.** A correlation pair passes SNP filtering then lookup of GCGC/GAGT fails. A filtering implementation error is mistaken for bad input that must be discarded wholesale.

**Legacy disposition.** Closed without a fix in this thread; identical reports #139/#140 later identify the unfiltered allele-vector bug.

**LDSC3 assessment — Apparently prevented.** LDSC3 builds the final merged pair, filters compatible pairs and flips signed statistics on that retained frame.
It does not align a filtered trait vector against stale unfiltered allele keys.
Standard standalone h2 success is not a test of this cross-trait operation, so pair-level allele fixtures are the relevant evidence.

**Recommended prevention.** Keep mixed valid/incompatible pair tests with nonconsecutive indices and verify count, retained identities and sign after filtering.

**Evidence and verification.** source-inspected
[src/ldsc/_regression/runner.py:build_rg_dataset](<../../../../src/ldsc/_regression/runner.py>); [tests/regression/test_regression_workflow.py:test_build_rg_dataset_harmonizes_alleles_and_flips_second_trait_z](<../../../../tests/regression/test_regression_workflow.py>).

Legacy implementation context: [`ldscore/regressions.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/regressions.py), [`ldscore/sumstats.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/sumstats.py).

<a id="issue-139"></a>

### #139: Running Error

[Original issue](https://github.com/bulik/ldsc/issues/139) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **1**.

Workflow tags: `rg`.

**Trigger and harm.** Same GCGC failure as #137 occurs after 340,046 compatible SNPs remain.
Users may unnecessarily discard most variants or rewrite munging to bypass a pair-merge bug.

**Legacy disposition.** Open duplicate-like report; no new fix. #140 contains a later precise stale-vector workaround.

**LDSC3 assessment — Apparently prevented.** LDSC3 filters and aligns on the merged pair rather than using a stale prefilter allele Series.
Allele-aware identities can also remove incompatible candidates earlier.
High INFO does not make two allele representations compatible, and omitting HM3 restriction does not remove the need for orientation checks.

**Recommended prevention.** Cross-link #137/#140 as one mechanism while counting each issue once; preserve explicit incompatibility diagnostics rather than advise hand-deleting allele classes.

**Evidence and verification.** source-inspected
[src/ldsc/_regression/runner.py:build_rg_dataset](<../../../../src/ldsc/_regression/runner.py>); [src/ldsc/_regression/estimators.py:_filter_alleles](<../../../../src/ldsc/_regression/estimators.py>); [tests/regression/test_regression_workflow.py](<../../../../tests/regression/test_regression_workflow.py>).

Legacy implementation context: [`ldscore/regressions.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/regressions.py), [`ldscore/sumstats.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/sumstats.py).

<a id="issue-140"></a>

### #140: Heritability have error

[Original issue](https://github.com/bulik/ldsc/issues/140) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **6**.

Workflow tags: `rg`.

**Trigger and harm.** The loop frame is filtered while the allele-key Series still includes incompatible pairs.
The subsequent sign dictionary lookup crashes despite earlier reporting valid alleles.

**Legacy disposition.** A user supplies the exact repair alleles=alleles[li] and surrounding lines; no merged upstream commit is verified.

**LDSC3 assessment — Apparently prevented.** LDSC3 constructs and harmonizes the retained merged data directly, with tests for orientation and effective-key joins.
This avoids stale-Series alignment rather than depending on a munging keep-list to hide the bug.
Treat the suggested removal of A/T or G/C classes carefully: strand ambiguity is distinct from this indexing error.

**Recommended prevention.** Retain one mixed-pair regression fixture demonstrating dropped incompatible SNPs and correct remaining Z signs; do not add legacy patch code to the new workflow.

**Evidence and verification.** source-inspected
[src/ldsc/_regression/runner.py:build_rg_dataset](<../../../../src/ldsc/_regression/runner.py>); [tests/regression/test_regression_workflow.py:test_build_rg_dataset_allele_aware_mode_harmonizes_after_identity_merge](<../../../../tests/regression/test_regression_workflow.py>); [Source](<https://github.com/bulik/ldsc/issues/140#issuecomment-1031683241>).

Legacy implementation context: [`ldscore/regressions.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/regressions.py), [`ldscore/sumstats.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/sumstats.py).

<a id="issue-158"></a>

### #158: LinAlgError: Singular matrix

[Original issue](https://github.com/bulik/ldsc/issues/158) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **10**.

Workflow tags: `baseline-partitioned-h2`, `ldsc-seg`, `continuous`.

**Trigger and harm.** Continuous or cell-type models have redundant/all-zero LD-score columns or lose support in a deletion. A raw singular-matrix error tempts users to tune block count until the run succeeds.

**Legacy disposition.** Several users identify duplicate all-ones bases, repeated control inputs or empty BEDs; no single universal legacy fix exists.

**LDSC3 assessment — Apparently prevented.** LDSC3 rejects duplicate names and collinear models, checks empty/support conditions, and reports failed delete blocks with parameter/support context. --continue-on-query-error can retain explicitly recorded other successes.
More blocks may change a sparse deletion problem but cannot identify coefficients in an exactly redundant full model.

**Recommended prevention.** Diagnose full-rank versus delete-block support first, include base exactly once, and avoid pseudoinverse, forced inversion or arbitrary block-count repairs.

**Evidence and verification.** source-inspected
[src/ldsc/_regression/overlap_model.py:model_collinearity_error](<../../../../src/ldsc/_regression/overlap_model.py>); [src/ldsc/_regression/jackknife.py](<../../../../src/ldsc/_regression/jackknife.py>); [src/ldsc/_annotation/query_outcomes.py](<../../../../src/ldsc/_annotation/query_outcomes.py>); [tests/regression/test_partitioned_query_failures.py](<../../../../tests/regression/test_partitioned_query_failures.py>).

Legacy implementation context: [`ldscore/regressions.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/regressions.py).

<a id="issue-186"></a>

### #186: L2 calculation error

[Original issue](https://github.com/bulik/ldsc/issues/186) · Source state: **closed** · Severity: **high** · Accessible comments reviewed: **0**.

Workflow tags: `ldscore`.

**Trigger and harm.** Changing --ld-wind-cm from1 to1e-6 leaves scores unchanged, with --yes-really enabled.
An unusable genetic map can silently change LD neighborhoods and invalidate reference scores.

**Legacy disposition.** Closed without comments or an attributable fix; the supplied log alone does not establish the root cause.

**LDSC3 assessment — Apparently prevented.** LDSC3 checks usable cM values or derives them from an explicit genetic map, and exact-window tests exercise boundary behavior.
Identical scores can also be legitimate if no pair crosses the changed window; reproduce with informative coordinates before alleging arithmetic failure.

**Recommended prevention.** Reject missing/constant genetic maps for cM windows and report the resolved map/window; preserve an exact boundary test.

**Evidence and verification.** source-inspected
[src/ldsc/_refpanel/adapters.py:prepare_chromosome](<../../../../src/ldsc/_refpanel/adapters.py>); [src/ldsc/_refpanel/adapters.py:assert_cm_usable](<../../../../src/ldsc/_refpanel/adapters.py>); [tests/ldscore/test_plink_exact_windows.py](<../../../../tests/ldscore/test_plink_exact_windows.py>).

Legacy implementation context: [`ldscore/ldscore.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/ldscore.py).

<a id="issue-214"></a>

### #214: Error parsing refernce panel LD Score

[Original issue](https://github.com/bulik/ldsc/issues/214) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **4**.

Workflow tags: `baseline-partitioned-h2`, `continuous`.

**Trigger and harm.** A path prefix valid in1.0.0 yields no objects to concatenate after wildcard expansion changed.
Documented partitioned-heritability workflows fail before regression.

**Legacy disposition.** Verified master fix: aa33296abac9569a6422ee6ba7eb4b902422cc74, PR221, moves glob handling out of _splitp.

**LDSC3 assessment — Apparently prevented.** LDSC3 uses shared exact/glob/chromosome-suite resolvers and explicit LD-score-directory conversion instead of inheriting the buggy legacy _splitp.
It also produces actionable missing-input errors rather than an empty concat.
The fork commits referenced in events are not the authoritative master fix.

**Recommended prevention.** Preserve separate exact-prefix and glob tests and document @ chromosome suites and legacy conversion.

**Evidence and verification.** source-inspected
[src/ldsc/_runtime/paths.py:resolve_file_group](<../../../../src/ldsc/_runtime/paths.py>); [src/ldsc/_ldscore/legacy_converter.py](<../../../../src/ldsc/_ldscore/legacy_converter.py>); [Source](<https://github.com/bulik/ldsc/commit/aa33296abac9569a6422ee6ba7eb4b902422cc74>).

Legacy implementation context: [`ldscore/regressions.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/regressions.py).

<a id="issue-223"></a>

### #223: Partitioned heritability: Error parsing reference panel LD Score.

[Original issue](https://github.com/bulik/ldsc/issues/223) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **4**.

Workflow tags: `baseline-partitioned-h2`.

**Trigger and harm.** Partitioned-heritability tutorial prefixes produce ValueError:No objects to concatenate. A path-expansion regression blocks the advertised workflow.

**Legacy disposition.** Users confirm rollback b02f2a6; authoritative later fix is merged PR221/aa33296abac9569a6422ee6ba7eb4b902422cc74.

**LDSC3 assessment — Apparently prevented.** LDSC3 no longer routes chromosome prefixes through legacy _splitp.
Shared input resolution validates selected files and legacy conversion requires complete matching chromosome families; build-aware resources also address the supplied hg38 frequency-versus-baseline concern separately from this parsing bug.

**Recommended prevention.** Keep actionable empty-selection diagnostics and tests for exact suites, globs and genome-build compatibility.

**Evidence and verification.** source-inspected
[src/ldsc/_runtime/paths.py:resolve_chromosome_group](<../../../../src/ldsc/_runtime/paths.py>); [src/ldsc/_ldscore/legacy_converter.py](<../../../../src/ldsc/_ldscore/legacy_converter.py>); [Source](<https://github.com/bulik/ldsc/commit/aa33296abac9569a6422ee6ba7eb4b902422cc74>).

Legacy implementation context: [`ldscore/regressions.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/regressions.py).

<a id="issue-234"></a>

### #234: annotations with zero variance LD scores

[Original issue](https://github.com/bulik/ldsc/issues/234) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **0**.

Workflow tags: `baseline-partitioned-h2`, `ldsc-seg`.

**Trigger and harm.** Removing constant LD-score columns updates counts but not the overlap matrix.
Model summaries fail or pair the wrong overlap entries with coefficients.

**Legacy disposition.** Reporter proposes two-axis subsetting; no merged fixing commit is identified.

**LDSC3 assessment — Apparently prevented.** LDSC3 records retained/dropped column names, subsets count vectors consistently and builds the overlap matrix by retained names rather than reusing its original positional shape.
This directly addresses the reported alignment mechanism; a wholly constant model still fails explicitly rather than yielding fictitious coefficients.

**Recommended prevention.** Keep a regression case with a constant annotation between two retained columns; disclose omitted coefficients and counts.

**Evidence and verification.** source-inspected
[src/ldsc/_regression/runner.py:_dataset_from_query_values](<../../../../src/ldsc/_regression/runner.py>); [src/ldsc/_regression/overlap_model.py:assemble_model_overlap](<../../../../src/ldsc/_regression/overlap_model.py>).

Legacy implementation context: [`ldscore/regressions.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/regressions.py).

<a id="issue-306"></a>

### #306: Problem with duplicates when making annotation file [make_annot.py]

[Original issue](https://github.com/bulik/ldsc/issues/306) · Source state: **closed** · Severity: **high** · Accessible comments reviewed: **0**.

Workflow tags: `ldscore`.

**Trigger and harm.** A BED/BIM intersection duplicates BP rows and left join creates more rows than the BIM.
Row multiplication can corrupt SNP alignment and LD-score generation.

**Legacy disposition.** PR307 proposed drop_duplicates but closed unmerged; closure of issue306 is not evidence master adopted the patch.

**LDSC3 assessment — Apparently prevented.** LDSC3 ChromosomeProjector projects interval union membership directly onto fixed annotation-row positions, so overlaps never multiply rows.
Duplicate BIM identities are handled separately.
The deliberate binary union also differs from requested overlap-count annotations, which remain a separate unsupported feature.

**Recommended prevention.** Retain overlapping/nested interval tests with exactly one output row per grid SNP; document binary union versus counts explicitly.

**Evidence and verification.** source-inspected
[src/ldsc/_annotation/interval_projection.py:ChromosomeProjector.project](<../../../../src/ldsc/_annotation/interval_projection.py>); [Source](<https://github.com/bulik/ldsc/pull/307>); [docs/reference/annotate.md](<../../../reference/annotate.md>).

Legacy implementation context: [`make_annot.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/make_annot.py).

<a id="issue-327"></a>

### #327: munge_sumstats.py_error "ValueError: No object to concatenate"

[Original issue](https://github.com/bulik/ldsc/issues/327) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **10**.

Workflow tags: `munge-sumstats`.

**Trigger and harm.** Munging ends at pd.concat(dat_list) with no objects despite reading many chunks. A completely removed dataset yields an opaque error instead of showing why no SNPs survived.

**Legacy disposition.** Some users fix extra all-NA columns or change coordinate IDs to rsIDs; chunk-size advice fails for others.

**LDSC3 assessment — Apparently prevented.** LDSC3 reads selected columns, permits coordinate-family identity and handles an empty dat_list before raising retained0-of-N with stage counts.
This directly improves the common failure mechanism.
Tabs are valid whitespace; a comment claiming tab delimiters are categorically invalid is incorrect.

**Recommended prevention.** Use infer-only and stage-wise removal counts to distinguish ID mismatch, all-missing selected fields and QC; never relax scientific thresholds blindly.

**Evidence and verification.** source-inspected
[src/ldsc/_sumstats/munge_input.py:prepare_munge_input](<../../../../src/ldsc/_sumstats/munge_input.py>); [src/ldsc/_sumstats/munge_kernel.py:munge_sumstats](<../../../../src/ldsc/_sumstats/munge_kernel.py>).

Legacy implementation context: [`munge_sumstats.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/munge_sumstats.py).

<a id="issue-342"></a>

### #342: Bug: l2_parser can fail to load ldscores

[Original issue](https://github.com/bulik/ldsc/issues/342) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **0**.

Workflow tags: `baseline-partitioned-h2`.

**Trigger and harm.** Legacy pandas concatenation reorders annotation columns, producing 11 count columns versus 10 LD-score columns.
Count/coefficient misalignment can crash or corrupt annotation attribution.

**Legacy disposition.** Author proposes PR #341; API verification shows the PR remains open and unmerged.
No merged upstream fix is established.

**LDSC3 assessment — Apparently prevented.** LDSC3 uses named columns in LD-score directories and validates M shape against the model matrix.
Legacy conversion remains the relevant import boundary and must retain one explicit ordered annotation schema.

**Recommended prevention.** Verify converter rejection of shuffled/missing annotation columns and retain an order-permutation regression test; inspect PR #341 before attributing an upstream fix.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[src/ldsc/_ldscore/legacy_converter.py](<../../../../src/ldsc/_ldscore/legacy_converter.py>); [src/ldsc/_regression/estimators.py:263](<../../../../src/ldsc/_regression/estimators.py#L263>); [Source](<https://github.com/bulik/ldsc/pull/341>).

Legacy implementation context: [`ldscore/regressions.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/regressions.py).

<a id="issue-356"></a>

### #356: Error converting summary statistics

[Original issue](https://github.com/bulik/ldsc/issues/356) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **1**.

Workflow tags: `munge-sumstats`.

**Trigger and harm.** Legacy allele_merge applies bitwise negation to an integer-valued allele-match Series, yielding -1/-2 row labels.
Valid summary statistics fail after expensive filtering; an unchecked mask could also select incorrect rows.

**Legacy disposition.** Open; the later comment asks whether it was resolved.

**LDSC3 assessment — Apparently prevented.** LDSC3 no longer uses the legacy NaN-padding allele_merge assignment; allele-aware restriction forms explicit identity masks and filters chunks.
This is source evidence against the exact old indexing mechanism, not reproduction on the original GWAS.

**Recommended prevention.** Retain a regression test with matched and mismatched alleles after row filtering and assert Boolean dtype/alignment before selection.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[src/ldsc/_sumstats/munge_kernel.py:458](<../../../../src/ldsc/_sumstats/munge_kernel.py#L458>); [tests/sumstats/test_munge_hm3_defaults.py](<../../../../tests/sumstats/test_munge_hm3_defaults.py>).

Legacy implementation context: [`munge_sumstats.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/munge_sumstats.py).

<a id="issue-493"></a>

### #493: make_annot.py generates duplicate output lines

[Original issue](https://github.com/bulik/ldsc/issues/493) · Source state: **closed** · Severity: **high** · Accessible comments reviewed: **2**.

Workflow tags: `ldscore`.

**Trigger and harm.** Legacy make_annot joins on BP alone, multiplying rows at shared positions or chromosomes. A thin annotation can silently lose one-to-one alignment with genotype rows.

**Legacy disposition.** Commenter requires per-chromosome unique positions; the reporter accepts this restriction.
Closure does not establish an implementation fix.

**LDSC3 assessment — Apparently prevented.** LDSC3 projects directly onto normalized SNP rows and uses explicit identifier modes, rather than a BP-only left merge.
Same-position biallelic variants can remain distinct under allele-aware identity, while duplicate identities follow cleanup policy.

**Recommended prevention.** Retain tests for shared positions across chromosomes and distinct allele pairs; never solve row multiplication by globally discarding scientifically distinct same-position variants.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[src/ldsc/_annotation/interval_projection.py](<../../../../src/ldsc/_annotation/interval_projection.py>); [src/ldsc/_annotation/identity_cleanup.py](<../../../../src/ldsc/_annotation/identity_cleanup.py>); [src/ldsc/_annotation/bim_grid.py](<../../../../src/ldsc/_annotation/bim_grid.py>).

Legacy implementation context: [`ldscore/ldscore.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/ldscore.py).

<a id="issue-502"></a>

### #502: Partitioned LDSC failing for 2 gene sets.

[Original issue](https://github.com/bulik/ldsc/issues/502) · Source state: **closed** · Severity: **high** · Accessible comments reviewed: **1**.

Workflow tags: `ldscore`, `ldsc-seg`.

**Trigger and harm.** Two query gene lists produce more annotation rows than BIM/FRQ rows and fail Boolean frequency filtering.
Row multiplicity can silently destroy alignment before an eventual parser error.

**Legacy disposition.** Author reports applying PR #307 fixed their case.
API verifies PR #307 is closed but unmerged; this is a local workaround, not an upstream merge fix.

**LDSC3 assessment — Apparently prevented.** LDSC3 interval projection operates on SNP rows and identity cleanup avoids a BP-only multiplying join.
Its full annotation format and source/grid alignment checks address the mechanism without requiring thin annotations or blindly dropping duplicate output lines.

**Recommended prevention.** Keep a regression test with repeated positions across SNPs/chromosomes and require one output per retained SNP identity; preserve distinct alleles instead of a broad drop_duplicates workaround.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[src/ldsc/_annotation/interval_projection.py](<../../../../src/ldsc/_annotation/interval_projection.py>); [src/ldsc/_annotation/identity_cleanup.py](<../../../../src/ldsc/_annotation/identity_cleanup.py>); [Source](<https://github.com/bulik/ldsc/pull/307>).

Legacy implementation context: [`make_annot.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/make_annot.py), [`ldscore/regressions.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/regressions.py).

<a id="issue-227"></a>

### #227: Error in using munge_sumstats.py could not find SNP column

[Original issue](https://github.com/bulik/ldsc/issues/227) · Source state: **closed** · Severity: **moderate** · Accessible comments reviewed: **2**.

Workflow tags: `munge-sumstats`.

**Trigger and harm.** The raw file has quoted column names and uses --signed-sumstats Z,6.
Header parsing blocks input, and a nonzero Z null would separately misassign signs.

**Legacy disposition.** A commenter suggests removing quotes; later user reports similar failure without quotes, so closure is not a verified general fix.

**LDSC3 assessment — Confirmed persistent behavior.** An actual LDSC3 CLI probe with a quoted whitespace header reports missing N/direction in --infer-only and fails full munging with an unrecognized signed-statistic header, despite columns SNP/A1/A2/P/N/Z being present. read_header splits raw tokens without quote handling while the pandas body reader handles quotes.
The original Z,6 argument is a separate sign-null error; a valid Z null is 0.

**Recommended prevention.** Use one header parser consistent with the body reader; until then document unquoted whitespace headers and the correct Z,0 null.

**Evidence and verification.** actual CLI probe: python -m ldsc munge-sumstats, --infer-only and full run; quoted-header.sumstats.txt; probe-quoted-header-b.log and probe-quoted-header-full-b.log
[src/ldsc/_sumstats/munge_input.py:read_header](<../../../../src/ldsc/_sumstats/munge_input.py>); [src/ldsc/_sumstats/munge_input.py:prepare_munge_input](<../../../../src/ldsc/_sumstats/munge_input.py>).

Legacy implementation context: [`munge_sumstats.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/munge_sumstats.py).
[Executed header checks V4](verification.md#v4-binary-projection-and-header-delimiters).

<a id="issue-2"></a>

### #2: repeat rsids

[Original issue](https://github.com/bulik/ldsc/issues/2) · Source state: **closed** · Severity: **moderate** · Accessible comments reviewed: **0**.

Workflow tags: `h2`, `rg`.

**Trigger and harm.** Repeated rsIDs occur outside or inside the final regression intersection.
An arbitrary retained duplicate can change the modeled SNP or abort an otherwise usable analysis.

**Legacy disposition.** Closed without an identified fixing commit; the requested behavior was to discard offending intersecting SNPs and warn.

**LDSC3 assessment — Partial coverage.** LDSC3 has effective-identity uniqueness checks and drops duplicate effective keys during explicitly allowed identity downgrade.
Its legacy-summary-statistics projection drops duplicated source rsIDs before projection; this is safer than picking one but can remove duplicates outside the final intersection.
Exact equivalence to the requested intersection-only policy is therefore unproved.

**Recommended prevention.** Retain deterministic drop-all behavior, document the stage and affected counts, and add an intersection-specific audit case before changing policy.

**Evidence and verification.** source-inspected
[src/ldsc/_regression/runner.py:_project_legacy_sumstats_to_panel](<../../../../src/ldsc/_regression/runner.py>); [src/ldsc/_regression/runner.py:1794](<../../../../src/ldsc/_regression/runner.py#L1794>); [tests/regression/test_regression_workflow.py:test_identity_downgrade_drops_duplicate_effective_keys_before_merge_and_logs](<../../../../tests/regression/test_regression_workflow.py>).

Legacy implementation context: [`ldscore/regressions.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/regressions.py), [`ldscore/sumstats.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/sumstats.py).

<a id="issue-187"></a>

### #187: Munge error

[Original issue](https://github.com/bulik/ldsc/issues/187) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **1**.

Workflow tags: `munge-sumstats`.

**Trigger and harm.** The same file finishes on another machine but appears to run indefinitely locally.
Resource-dependent preprocessing failures are difficult to distinguish from a hang.

**Legacy disposition.** Reporter redirects to issue181, where smaller chunks are the confirmed user workaround.

**LDSC3 assessment — Partial coverage.** LDSC3 exposes chunk sizing and logs preparation progress, but retains the full filtered table for global N filtering and identity cleanup.
Smaller chunks reduce peak parsing pressure rather than guaranteeing bounded total memory; no supplied input proves an algorithmic loop.

**Recommended prevention.** Link chunk-size and retained-table memory guidance; report progress with a reproducible row count and runtime.

**Evidence and verification.** source-inspected
[src/ldsc/_sumstats/munge_kernel.py:parse_dat](<../../../../src/ldsc/_sumstats/munge_kernel.py>); [src/ldsc/_sumstats/munger.py:1417](<../../../../src/ldsc/_sumstats/munger.py#L1417>).

Legacy implementation context: [`munge_sumstats.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/munge_sumstats.py).

<a id="issue-216"></a>

### #216: munge_sumstats.py stuck and taking too much time

[Original issue](https://github.com/bulik/ldsc/issues/216) · Source state: **closed** · Severity: **moderate** · Accessible comments reviewed: **2**.

Workflow tags: `munge-sumstats`.

**Trigger and harm.** Munging16.4million rows appears stuck during five-million-row input chunks.
High memory and no visible progress obstruct preprocessing.

**Legacy disposition.** Confirmed user workaround --chunksize500000; reporter says it works and closes.

**LDSC3 assessment — Partial coverage.** LDSC3 defaults to one-million-row chunks with phase progress and HM3 restriction, addressing this common pressure point.
It still concatenates retained rows and computes global N filters, so a smaller chunk is a practical workaround rather than a guarantee of bounded-memory execution.

**Recommended prevention.** Document chunk sizing and progress expectations; use row-retention diagnostics to distinguish a hang from memory pressure.

**Evidence and verification.** source-inspected
[src/ldsc/_sumstats/munge_kernel.py:parse_dat](<../../../../src/ldsc/_sumstats/munge_kernel.py>); [src/ldsc/_sumstats/munger.py:1417](<../../../../src/ldsc/_sumstats/munger.py#L1417>).

Legacy implementation context: [`munge_sumstats.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/munge_sumstats.py).

<a id="issue-243"></a>

### #243: Error for ldsc.py

[Original issue](https://github.com/bulik/ldsc/issues/243) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **4**.

Workflow tags: `rg`.

**Trigger and harm.** Regression cannot parse a .sumstats.gz despite reported successful munging.
Wrong compression or a broken dependency environment blocks analysis with confusing empty-file diagnostics.

**Legacy disposition.** Users report dependency updates or removing a false .gz suffix; neither establishes a universal code fix.

**LDSC3 assessment — Partial coverage.** LDSC3 uses current pandas and a clear legacy-summary-statistics loader, but compression selection still follows filename suffix. A plain text file mislabeled .gz will fail, as it should; error context should distinguish malformed compression from an empty or wrong-schema file.

**Recommended prevention.** Validate actual gzip integrity and identify the source path; do not rename a genuinely corrupted archive as a repair.

**Evidence and verification.** source-inspected
[src/ldsc/_sumstats/munger.py:_read_curated_sumstats_artifact](<../../../../src/ldsc/_sumstats/munger.py>); [docs/reference/formats/legacy-sumstats.md](<../../../reference/formats/legacy-sumstats.md>).

Legacy implementation context: [`ldscore/regressions.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/regressions.py), [`ldscore/sumstats.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/sumstats.py).

<a id="issue-269"></a>

### #269: IOError: Not a gzipped file

[Original issue](https://github.com/bulik/ldsc/issues/269) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **0**.

Workflow tags: `rg`.

**Trigger and harm.** Genetic correlation raises Not a gzipped file for a .sumstats.gz input.
Users cannot tell false suffix, corruption and wrong file selection apart.

**Legacy disposition.** Open; author asserts files are gzipped but supplies no reproducible bytes or fix.

**LDSC3 assessment — Partial coverage.** LDSC3 legacy loading and raw input compression still rely on suffix-based selection. A mislabeled plain file or damaged archive therefore remains an input failure, rather than evidence of a regression estimator bug.
Native Parquet output avoids this particular compression ambiguity.

**Recommended prevention.** Report filename and gzip-integrity context clearly; use a content check for diagnosis without silently changing how the file is interpreted.

**Evidence and verification.** source-inspected
[src/ldsc/_sumstats/munger.py:_read_curated_sumstats_artifact](<../../../../src/ldsc/_sumstats/munger.py>); [src/ldsc/_sumstats/munge_input.py:get_compression](<../../../../src/ldsc/_sumstats/munge_input.py>).

Legacy implementation context: [`ldscore/regressions.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/regressions.py), [`ldscore/sumstats.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/sumstats.py).

<a id="issue-329"></a>

### #329: Munge Error - ValueError: could not convert string to float: OR

[Original issue](https://github.com/bulik/ldsc/issues/329) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **4**.

Workflow tags: `munge-sumstats`.

**Trigger and harm.** Parser cannot convert the tokenOR to float in an effect column. A repeated header or malformed row produces a low-level error before useful row diagnostics.

**Legacy disposition.** Comments recommend numeric-column inspection, with mixed success; no universal legacy fix.

**LDSC3 assessment — Partial coverage.** LDSC3 detects nonnumeric mapped columns, but signed columns are forced to float by pandas while reading chunks, so some malformed tokens can fail before the custom numeric-column message.
Automatically dropping them would risk hiding structural errors; exact source-row context is the simple missing usability improvement.

**Recommended prevention.** Wrap parsing exceptions with source path, column and likely repeated-header/malformed-row guidance; preserve the original exception and fail rather than coercing silently.

**Evidence and verification.** source-inspected
[src/ldsc/_sumstats/munge_kernel.py:munge_sumstats](<../../../../src/ldsc/_sumstats/munge_kernel.py>); [src/ldsc/_sumstats/munge_kernel.py:parse_dat](<../../../../src/ldsc/_sumstats/munge_kernel.py>).

Legacy implementation context: [`munge_sumstats.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/munge_sumstats.py).

<a id="issue-345"></a>

### #345: FloatingPointError: invalid value encountered in divide

[Original issue](https://github.com/bulik/ldsc/issues/345) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **3**.

Workflow tags: `h2`.

**Trigger and harm.** Non-European heritability fit raises invalid division while deriving enrichment. A completed fit can lose its result to an undefined ratio, obscuring whether counts or near-zero total heritability caused it.

**Legacy disposition.** Open; later users report the same exception without a verified fix.

**LDSC3 assessment — Partial coverage.** The shared estimator still computes cat/M divided by tot/M_tot without an explicit zero-total guard.
Current workflow count validation removes some invalid inputs, but this does not establish handling of a legitimate zero total.
Ancestry itself is not a diagnosed cause.

**Recommended prevention.** Preserve finite fitted quantities and explicitly mark undefined derived ratios with a reason; reproduce zero-count and zero-total cases separately before changing numerical policy.

**Evidence and verification.** Root V3 local zero-total Hsq probe raises FloatingPointError in derived-ratio handling.
The historical non-European input was not available to establish its exact cause.
[src/ldsc/_regression/estimators.py:424](<../../../../src/ldsc/_regression/estimators.py#L424>); [src/ldsc/_regression/runner.py](<../../../../src/ldsc/_regression/runner.py>).

Legacy implementation context: [`ldscore/regressions.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/regressions.py).
[Shared zero-total boundary V3](verification.md#v3-zero-total-derived-proportion).

<a id="issue-211"></a>

### #211: Error when following the tutorial

[Original issue](https://github.com/bulik/ldsc/issues/211) · Source state: **closed** · Severity: **moderate** · Accessible comments reviewed: **1**.

Workflow tags: `baseline-partitioned-h2`.

**Trigger and harm.** A tutorial partitioned-heritability run fails when taking the square root of jackknife variance.
Standard errors cannot be produced despite substantial input and memory.

**Legacy disposition.** Reporter succeeds on Google Cloud and closes as a cluster-environment problem; no fixing commit or precise numerical cause is supplied.

**LDSC3 assessment — Indeterminate.** LDSC3 computes jackknife covariance and square roots through modern NumPy, but source inspection alone cannot prove the same environment-specific error is impossible.
The supplied report does not distinguish roundoff, a bad BLAS environment or invalid model inputs.

**Recommended prevention.** Record supported dependency/BLAS versions and reproduce with minimal inputs; diagnose negative variance instead of globally suppressing floating-point errors.

**Evidence and verification.** source-inspected
[src/ldsc/_regression/jackknife.py:Jackknife.jknife](<../../../../src/ldsc/_regression/jackknife.py>); [tests/regression/test_jackknife.py](<../../../../tests/regression/test_jackknife.py>); [docs/user-guide/installation.md](<../../../user-guide/installation.md>).

Legacy implementation context: [`ldscore/regressions.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/regressions.py).

<a id="issue-19"></a>

### #19: misformatted .annot doesn't raise sensible error

[Original issue](https://github.com/bulik/ldsc/issues/19) · Source state: **closed** · Severity: **moderate** · Accessible comments reviewed: **1**.

Workflow tags: `ldscore`.

**Trigger and harm.** An annotation supplies SNP/values but omits required identity metadata.
Late pandas errors hide the source and delay repair.

**Legacy disposition.** Verified commit 4b6f38d3d0e8a9eaff500752443d98f8c2a47965 closes #19; maintainer describes wrapping annotation reads in actionable error handling.

**LDSC3 assessment — Apparently prevented.** LDSC3 declares annotation inputs and checks headers before numerical LD-score work.
Parsing normalizes identity fields and rejects missing/non-numeric annotation values; diagnostics identify input roles and repair actions. A file valid for rsID matching may legitimately need fewer fields than a coordinate-family file, so validation should explain the selected SNP identifier mode rather than impose obsolete fixed columns.

**Recommended prevention.** Retain mode-specific preflight and test a malformed minimal table through the public command, not only private parsing.

**Evidence and verification.** source-inspected
[src/ldsc/_annotation/preflight.py](<../../../../src/ldsc/_annotation/preflight.py>); [src/ldsc/_annotation/parsing.py](<../../../../src/ldsc/_annotation/parsing.py>); [tests/annotation/test_annotation_format.py](<../../../../tests/annotation/test_annotation_format.py>); [Source](<https://github.com/bulik/ldsc/commit/4b6f38d3d0e8a9eaff500752443d98f8c2a47965>).

Legacy implementation context: [`make_annot.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/make_annot.py).

<a id="issue-21"></a>

### #21: --ref-ld-list w/ identical column names

[Original issue](https://github.com/bulik/ldsc/issues/21) · Source state: **closed** · Severity: **moderate** · Accessible comments reviewed: **0**.

Workflow tags: `baseline-partitioned-h2`.

**Trigger and harm.** Multiple LD-score files contain the same annotation column name.
Duplicate names can merge or mislabel coefficients and overlaps.

**Legacy disposition.** Verified commit 2d9d7c564130cdab8f8c1331cbf091a0e845d938 changes two test headers to exercise #21; this demonstrates a fixture update rather than a numerical fix.

**LDSC3 assessment — Apparently prevented.** LDSC3 rejects duplicate annotation names across baseline and query groups before construction. A saved LD-score directory also carries explicit ordered columns and selected reads.
Rejection is a deliberate safer interface than guessing suffixes; users must rename colliding annotation columns with scientifically meaningful labels.

**Recommended prevention.** Keep the error explicit about every colliding source/name and show a rename example in troubleshooting.

**Evidence and verification.** source-inspected
[src/ldsc/_annotation/semantics.py:require_unique_annotation_names](<../../../../src/ldsc/_annotation/semantics.py>); [tests/annotation/test_annotation_semantics.py:test_duplicate_names_across_groups_are_rejected](<../../../../tests/annotation/test_annotation_semantics.py>); [Source](<https://github.com/bulik/ldsc/commit/2d9d7c564130cdab8f8c1331cbf091a0e845d938>).

Legacy implementation context: [`ldscore/regressions.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/regressions.py).

<a id="issue-28"></a>

### #28: Appropriate Error Message for missing singed-sumstats

[Original issue](https://github.com/bulik/ldsc/issues/28) · Source state: **closed** · Severity: **moderate** · Accessible comments reviewed: **1**.

Workflow tags: `munge-sumstats`.

**Trigger and harm.** A raw GWAS header lacks a recognized beta, odds ratio or Z field. A non-actionable list-index traceback blocks preprocessing.

**Legacy disposition.** Verified local patch dadf6e2405d9bb4d7f0ab24746b2b6ee16e57f30, titled close issue #28, adds an explicit len(sign_cnames)==0 ValueError before indexing.
The unrelated fd06bc9 event was unavailable locally.

**LDSC3 assessment — Apparently prevented.** LDSC3 resolves the input schema before QC, rejects ambiguous mappings, and requires a valid signed statistic/null when needed.
The public munger offers infer-only inspection without writing outputs. A P-only file can support unsigned heritability under appropriate settings but cannot supply genetic-correlation signs.

**Recommended prevention.** Make the missing-field error list accepted signed-statistic options and explain the P-only heritability versus correlation boundary.

**Evidence and verification.** source-inspected
[src/ldsc/_sumstats/munge_input.py](<../../../../src/ldsc/_sumstats/munge_input.py>); [docs/reference/munge-sumstats.md](<../../../reference/munge-sumstats.md>); [tests/sumstats/test_munge_infer_diagnostics.py](<../../../../tests/sumstats/test_munge_infer_diagnostics.py>); [Source](<https://github.com/bulik/ldsc/commit/dadf6e2405d9bb4d7f0ab24746b2b6ee16e57f30>).

Legacy implementation context: [`munge_sumstats.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/munge_sumstats.py).

<a id="issue-34"></a>

### #34: error for generate ldscore file when no annotation in that chromosome

[Original issue](https://github.com/bulik/ldsc/issues/34) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **4**.

Workflow tags: `ldscore`, `baseline-partitioned-h2`.

**Trigger and harm.** A valid genome-wide annotation is all zero on chromosome 22 but nonzero elsewhere.
Per-chromosome pruning can produce wrong shapes or discard legitimate sparse annotations.

**Legacy disposition.** Open; maintainer suggested concatenating chromosome scores and summing counts as a workaround, not a verified repair.

**LDSC3 assessment — Apparently prevented.** LDSC3 aggregates chromosome results into one LD-score directory and measures query variation across the complete run.
Zero is a valid annotation value.
Global zero-support queries and empty input sources are handled separately from chromosome-local zeros, while regression can still be unestimable if all supporting SNPs occur in one jackknife block.

**Recommended prevention.** Preserve columns across chromosome-local zeros and add a sparse-one-chromosome case to the audit probes; never manufacture support with nonzero filler.

**Evidence and verification.** source-inspected
[src/ldsc/_ldscore/query_variation.py](<../../../../src/ldsc/_ldscore/query_variation.py>); [src/ldsc/_ldscore/calculator.py:_aggregate_chromosome_results](<../../../../src/ldsc/_ldscore/calculator.py>); [tests/ldscore/test_ldscore_query_batches.py](<../../../../tests/ldscore/test_ldscore_query_batches.py>).

Legacy implementation context: [`make_annot.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/make_annot.py), [`ldscore/regressions.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/regressions.py).

<a id="issue-47"></a>

### #47: munge_sumstats.py error

[Original issue](https://github.com/bulik/ldsc/issues/47) · Source state: **closed** · Severity: **moderate** · Accessible comments reviewed: **3**.

Workflow tags: `munge-sumstats`.

**Trigger and harm.** Both MarkerName and rsID are automatically recognized after an explicit --snp selection.
The munger receives duplicate SNP columns and crashes with a pandas truth-value error.

**Legacy disposition.** Ignoring MarkerName fixed the user input.
Maintainer cites PR #75; local merged history 4da1199468ef6866a61f0f3252e212e2972fe5ed includes duplicate-name/canonical-field diagnostic checks, verified by inspecting the contained patches.

**LDSC3 assessment — Apparently prevented.** LDSC3 resolves explicit column choices and rejects duplicate mappings before chunk QC, with instructions to rename/remove conflicting labels.
This is the correct seam for the failure; coercing the resulting DataFrame to a Series would hide the ambiguity.
Existing sample-size precedence tests do not substitute for an identifier-alias case.

**Recommended prevention.** Preserve explicit-selection precedence and a diagnostic listing every competing source column; map this exact dual-ID case to a focused test.

**Evidence and verification.** source-inspected
[src/ldsc/_sumstats/munge_input.py:prepare_munge_input](<../../../../src/ldsc/_sumstats/munge_input.py>); [tests/sumstats/test_munge_infer_diagnostics.py](<../../../../tests/sumstats/test_munge_infer_diagnostics.py>); [Source](<https://github.com/bulik/ldsc/pull/75>); [Source](<https://github.com/bulik/ldsc/commit/4da1199468ef6866a61f0f3252e212e2972fe5ed>).

Legacy implementation context: [`munge_sumstats.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/munge_sumstats.py).

<a id="issue-128"></a>

### #128: Different number of annotations per chromosome

[Original issue](https://github.com/bulik/ldsc/issues/128) · Source state: **closed** · Severity: **moderate** · Accessible comments reviewed: **4**.

Workflow tags: `ldscore`, `baseline-partitioned-h2`.

**Trigger and harm.** One chromosome lacks a column present in every other shard. M counts and predictors have incompatible dimensions.

**Legacy disposition.** Maintainer suggests concatenating/regenerating aligned chromosome data; no verified general repair.

**LDSC3 assessment — Apparently prevented.** LDSC3 uses a common named annotation schema and publishes aggregated LD scores. A valid globally defined annotation must retain a zero column on unsupported chromosomes, not delete the column.
The global-variation check can remove a wholly absent query separately.
Merging malformed shards without fixing their column meaning is not sufficient.

**Recommended prevention.** Explain missing-column versus all-zero-column semantics and keep a chromosome-local zero example; never infer zero for unknown or unreadable annotations.

**Evidence and verification.** source-inspected
[src/ldsc/_annotation/sources.py](<../../../../src/ldsc/_annotation/sources.py>); [src/ldsc/_ldscore/query_variation.py](<../../../../src/ldsc/_ldscore/query_variation.py>); [src/ldsc/_ldscore/source.py](<../../../../src/ldsc/_ldscore/source.py>).

Legacy implementation context: [`make_annot.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/make_annot.py), [`ldscore/regressions.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/regressions.py).

<a id="issue-160"></a>

### #160: Error in munge_sumstats.py if multiple columns have a "default name"

[Original issue](https://github.com/bulik/ldsc/issues/160) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **1**.

Workflow tags: `munge-sumstats`.

**Trigger and harm.** A file has RSID and SNPID and the user selects one explicitly.
Useful raw schemas require unnecessary destructive renaming.

**Legacy disposition.** Open request for explicit options to override defaults; no verified legacy fix.

**LDSC3 assessment — Apparently prevented.** LDSC3 get_cname_map collects used canonical targets from explicit hints and excludes automatic aliases for those targets.
Thus --snp SNPID prevents RSID from also becoming SNP, while genuine duplicate raw header labels remain rejected.
This is precise source-level prevention, not merely a more informative rejection.

**Recommended prevention.** Retain an explicit dual-ID precedence fixture and document the selected source in infer-only output; do not add another identifier-renaming option.

**Evidence and verification.** source-inspected
[src/ldsc/_sumstats/munge_input.py:get_cname_map](<../../../../src/ldsc/_sumstats/munge_input.py>); [src/ldsc/_core/columns.py](<../../../../src/ldsc/_core/columns.py>).

Legacy implementation context: [`munge_sumstats.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/munge_sumstats.py).

<a id="issue-175"></a>

### #175: KeyError: u"None of [Int64Index([-1, -2, -1, -2, -2, -1, -1, -1, -1, -1,\n

[Original issue](https://github.com/bulik/ldsc/issues/175) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **2**.

Workflow tags: `munge-sumstats`.

**Trigger and harm.** Munging raises KeyError with -1/-2 integer labels under newer pandas.
Normal input cannot be munged when a Boolean mask becomes an integer indexer.

**Legacy disposition.** Merged PR185 fixes allele_merge; merge be1985677b8301e3f9e2524b9708e380186bad29. [Source](<https://github.com/bulik/ldsc/pull/185>)

**LDSC3 assessment — Apparently prevented.** LDSC3 uses current pandas and explicit Boolean identity/allele masks rather than the failing legacy integer-index expression.
The installation distinction is between the interpreter in the environment and the Anaconda distribution, not whether Anaconda3 was installed.

**Recommended prevention.** Retain Boolean-mask regression coverage across supported pandas versions and explain supported Python versions.

**Evidence and verification.** source-inspected
[src/ldsc/_sumstats/munge_kernel.py:_filter_sumstats_chunk_to_restriction](<../../../../src/ldsc/_sumstats/munge_kernel.py>); [tests/sumstats/test_munge_flag_contracts.py](<../../../../tests/sumstats/test_munge_flag_contracts.py>); [Source](<https://github.com/bulik/ldsc/pull/185>).

Legacy implementation context: [`munge_sumstats.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/munge_sumstats.py).

<a id="issue-198"></a>

### #198: make_annot: OverflowError: can't convert negative value to CHRPOS on recent version

[Original issue](https://github.com/bulik/ldsc/issues/198) · Source state: **closed** · Severity: **moderate** · Accessible comments reviewed: **1**.

Workflow tags: `ldscore`, `ldsc-seg`.

**Trigger and harm.** A gene near chromosome start is padded to one-based0 then converted to BED start-1. pybedtools raises OverflowError and annotation generation stops.

**Legacy disposition.** Confirmed merged PR199 clips the one-based start at1 before conversion; merge388e5aedfe325b7f8d932bb066c83367e664aa0d.

**LDSC3 assessment — Apparently prevented.** LDSC3 validates one-based gene catalog coordinates and its NumPy interval projector clips zero-based padded starts at0.
The implementation no longer depends on pybedtools CHRPOS conversion, removing this exact underflow while preserving correct boundary semantics.

**Recommended prevention.** Retain a chromosome-start padding example and boundary regression test; preserve one-based-versus-BED conventions.

**Evidence and verification.** source-inspected
[src/ldsc/_annotation/interval_projection.py:ChromosomeProjector.project](<../../../../src/ldsc/_annotation/interval_projection.py>); [src/ldsc/_annotation/gene_lists.py:_validate_catalog_rows](<../../../../src/ldsc/_annotation/gene_lists.py>); [Source](<https://github.com/bulik/ldsc/pull/199>).

Legacy implementation context: [`make_annot.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/make_annot.py), [`ldscore/regressions.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/regressions.py).

<a id="issue-232"></a>

### #232: Issue with new glob update

[Original issue](https://github.com/bulik/ldsc/issues/232) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **0**.

Workflow tags: `baseline-partitioned-h2`.

**Trigger and harm.** A short report says the new glob update throws an error.
Wildcard/prefix changes can block valid model inputs, but no traceback is supplied.

**Legacy disposition.** The known matching regression was fixed by PR221/aa33296; this thread alone does not demonstrate exact identity with that bug.

**LDSC3 assessment — Apparently prevented.** LDSC3 uses a dedicated input resolver with exact, glob and chromosome-suite semantics rather than the legacy splitter.
The report is consistent with issues214and223, but should retain its lower evidence specificity instead of being claimed as a separately reproduced defect.

**Recommended prevention.** Keep the known glob regression covered and request a concrete failing token/traceback for any distinct failure.

**Evidence and verification.** source-inspected
[src/ldsc/_runtime/paths.py:resolve_file_group](<../../../../src/ldsc/_runtime/paths.py>); [Source](<https://github.com/bulik/ldsc/commit/aa33296abac9569a6422ee6ba7eb4b902422cc74>).

Legacy implementation context: [`ldscore/regressions.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/regressions.py).

<a id="issue-239"></a>

### #239: Empty log files

[Original issue](https://github.com/bulik/ldsc/issues/239) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **1**.

Workflow tags: `h2`, `rg`.

**Trigger and harm.** A reporter reproduces an empty log when repeating the same legacy output prefix.
The only available diagnostic record can be lost, obscuring resource failures.

**Legacy disposition.** Deleting the old log and rerunning is reported as a workaround; no verified logging fix.

**LDSC3 assessment — Apparently prevented.** LDSC3 owns logs through workflow contexts and preflights output collisions; reruns require explicit --overwrite.
This removes silent reuse of a legacy write-binary logger. A scheduler kill before logging or inaccessible output remains a distinct possibility that cannot be diagnosed from the original large input.

**Recommended prevention.** Retain repeat-run/overwrite logging checks and provide a clear startup message naming the authoritative log.

**Evidence and verification.** source-inspected
[src/ldsc/_runtime/workflow_logging.py](<../../../../src/ldsc/_runtime/workflow_logging.py>); [docs/reference/logs-and-outputs.md](<../../../reference/logs-and-outputs.md>).

Legacy implementation context: [`ldscore/regressions.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/regressions.py), [`ldscore/sumstats.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/sumstats.py).

<a id="issue-245"></a>

### #245: Issue with munge_sumstats.py : KeyError: u"None of [Int64Index([-1, -1, -1, -1, -1, -1, -1, -1, -1, -1,\n

[Original issue](https://github.com/bulik/ldsc/issues/245) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **2**.

Workflow tags: `munge-sumstats`.

**Trigger and harm.** --merge-alleles produces a -1/-2 index KeyError while omitting it succeeds.
Users may bypass allele QC or manually replace alleles without adjusting signed effects.

**Legacy disposition.** Same known legacy allele_merge mechanism fixed by merged PR185/be1985677b8301e3f9e2524b9708e380186bad29.

**LDSC3 assessment — Apparently prevented.** LDSC3 performs explicit SNP restriction and allele-compatible identity matching with Boolean masks.
The old flag combined row restriction with allele checks; copying HM3 alleles onto raw rows without orienting Z is not equivalent to a correct merge.

**Recommended prevention.** Explain restriction versus orientation and preserve signed statistics during any allele swap; retain the Boolean-mask regression case.

**Evidence and verification.** source-inspected
[src/ldsc/_sumstats/munge_kernel.py:_filter_sumstats_chunk_to_restriction](<../../../../src/ldsc/_sumstats/munge_kernel.py>); [Source](<https://github.com/bulik/ldsc/pull/185>); [docs/reference/munge-sumstats.md](<../../../reference/munge-sumstats.md>).

Legacy implementation context: [`munge_sumstats.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/munge_sumstats.py).

<a id="issue-272"></a>

### #272: munge_sumstats.py column names specified on the command line should override any defaults

[Original issue](https://github.com/bulik/ldsc/issues/272) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **1**.

Workflow tags: `munge-sumstats`.

**Trigger and harm.** --snp ldsc_id still conflicts with another automatically recognized rs column.
Correctly specified files are rejected or require fragile header renaming.

**Legacy disposition.** A commenter proposes --ignore as a workaround; no merged legacy fix is identified.

**LDSC3 assessment — Apparently prevented.** LDSC3 get_cname_map first maps explicit hints, records used canonical targets and suppresses competing default aliases for those targets.
Thus explicit SNP/frequency selection now has the precedence requested by this issue, while genuinely ambiguous automatic inference still fails.

**Recommended prevention.** Retain tests with both explicit and default aliases; show the resolved mapping in --infer-only.

**Evidence and verification.** source-inspected
[src/ldsc/_sumstats/munge_input.py:get_cname_map](<../../../../src/ldsc/_sumstats/munge_input.py>); [tests/sumstats/test_munge_flag_contracts.py](<../../../../tests/sumstats/test_munge_flag_contracts.py>).

Legacy implementation context: [`munge_sumstats.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/munge_sumstats.py).

<a id="issue-367"></a>

### #367: parse.py error?

[Original issue](https://github.com/bulik/ldsc/issues/367) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **1**.

Workflow tags: `baseline-partitioned-h2`.

**Trigger and harm.** Legacy annot parser computes the chromosome list outside its per-input loop.
Multi-prefix annotations can load the wrong chromosomes or fail before overlap-aware results.

**Legacy disposition.** Reporter proposes moving discovery inside the loop; a commenter mentions an unspecified PR, without verified merge evidence.

**LDSC3 assessment — Apparently prevented.** LDSC3 resolves annotation sources per chromosome and carries explicit source/column metadata, rather than using the old loop.
Compatibility still requires matching chromosome support across supplied files.

**Recommended prevention.** Keep a multi-source fixture with different missing chromosomes and assert a source-specific failure; do not globally reuse the first annotation source discovery result.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[src/ldsc/_annotation/sources.py](<../../../../src/ldsc/_annotation/sources.py>); [src/ldsc/_runtime/paths.py:269](<../../../../src/ldsc/_runtime/paths.py#L269>); [src/ldsc/_ldscore/legacy_converter.py](<../../../../src/ldsc/_ldscore/legacy_converter.py>).

Legacy implementation context: [`ldscore/regressions.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/regressions.py).

<a id="issue-422"></a>

### #422: NotImplementedError: iLocation based boolean indexing on an integer type is not available

[Original issue](https://github.com/bulik/ldsc/issues/422) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **0**.

Workflow tags: `ldscore`.

**Trigger and harm.** A third-party Python-3 distribution passes a pandas Boolean Series to iloc during --print-snps filtering.
LD-score output fails after computation.

**Legacy disposition.** Reporter proposes loc instead of iloc; no verified bulik/ldsc master fix is provided.

**LDSC3 assessment — Apparently prevented.** LDSC3 has separate explicit regression-SNP filtering and does not use the shown legacy df.iloc[df.SNP.isin(...)] path.
Exact row alignment remains important when identities are filtered.

**Recommended prevention.** Retain a nontrivial-index regression-SNP selection test and clearly identify the installed distribution in bug reports.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[src/ldsc/_ldscore/calculator.py](<../../../../src/ldsc/_ldscore/calculator.py>); [src/ldsc/_ldscore/compute.py](<../../../../src/ldsc/_ldscore/compute.py>); [tests/ldscore/test_reference_snp_funnel.py](<../../../../tests/ldscore/test_reference_snp_funnel.py>).

Legacy implementation context: [`ldscore/ldscore.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/ldscore.py).

<a id="issue-435"></a>

### #435: Error parsing .annot file

[Original issue](https://github.com/bulik/ldsc/issues/435) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **7**.

Workflow tags: `ldsc-seg`.

**Trigger and harm.** Python-3 annotation parsing uses fh after a comprehension, where it is no longer bound. A completed fit fails during overlap output; path-only advice fails to address the actual exception.

**Legacy disposition.** A commenter identifies changed comprehension scope and the reporter confirms changing fh to fh_list[0] works locally.
No merged fix is established.

**LDSC3 assessment — Apparently prevented.** LDSC3 builds overlap metadata during scoring and uses explicit per-source resolution.
It does not execute the shown legacy annot parser.
Source-specific chromosome validation is safer than assuming the first source represents every file.

**Recommended prevention.** Retain multi-source/unequal-scope validation and avoid reducing this historical case to a missing trailing dot or a missing file.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[src/ldsc/_annotation/sources.py](<../../../../src/ldsc/_annotation/sources.py>); [src/ldsc/_ldscore/outputs.py](<../../../../src/ldsc/_ldscore/outputs.py>); [src/ldsc/_regression/overlap_model.py](<../../../../src/ldsc/_regression/overlap_model.py>).

Legacy implementation context: [`ldscore/regressions.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/regressions.py).

<a id="issue-483"></a>

### #483: Error with calculating LD score using w_hm3.snplist: After merging with --print-snps, no SNPs remain.

[Original issue](https://github.com/bulik/ldsc/issues/483) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **1**.

Workflow tags: `ldscore`.

**Trigger and harm.** Legacy --print-snps reads a tabular SNP/A1/A2 file with comma defaults as one composite string column.
No identifiers match despite real overlap, and users patch parser code unnecessarily.

**Legacy disposition.** Author and commenter report selecting/exporting the SNP column as a workaround; no accepted upstream parser fix is shown.

**LDSC3 assessment — Apparently prevented.** LDSC3 uses a documented headered identity-file contract for regression SNP selection and its packaged HM3 default.
That contract differs from the legacy headerless one-column list and needs a clear migration example.

**Recommended prevention.** Validate required identity columns and separator before computation; show exact accepted keep-list examples rather than accepting arbitrary filename-based interpretations.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[docs/reference/legacy-flag-map.md:122](<../../../reference/legacy-flag-map.md#established-command-and-input-renames>); [src/ldsc/_core/columns.py](<../../../../src/ldsc/_core/columns.py>); [docs/reference/ldscore-snp-universes.md](<../../../reference/ldscore-snp-universes.md>).

Legacy implementation context: [`ldscore/ldscore.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/ldscore.py).

<a id="issue-500"></a>

### #500: ldsc39 ./munge_sumstats.py  has issue with *.gz file

[Original issue](https://github.com/bulik/ldsc/issues/500) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **0**.

Workflow tags: `munge-sumstats`.

**Trigger and harm.** A third-party Python-3.9 version reads a compressed header as bytes and calls rstrip with a text string.
Compressed input fails while the decompressed equivalent succeeds.

**Legacy disposition.** Open; the posted successful decompressed run is a workaround, not a bulik master fix.

**LDSC3 assessment — Apparently prevented.** LDSC3 explicitly decodes byte headers to UTF-8 before splitting, and supports gzip/bzip2/plain routes.
This targets the exact header mechanism; corrupted compressed data remain a different condition.

**Recommended prevention.** Retain equivalent plain/gzip/bzip2 header and munging tests with identical retained rows and signs; identify the installed distribution in reports.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[src/ldsc/_sumstats/munge_input.py:77](<../../../../src/ldsc/_sumstats/munge_input.py#L77>); [src/ldsc/_sumstats/munge_input.py:97](<../../../../src/ldsc/_sumstats/munge_input.py#L97>); [src/ldsc/_sumstats/munge_input.py:118](<../../../../src/ldsc/_sumstats/munge_input.py#L118>).

Legacy implementation context: [`munge_sumstats.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/munge_sumstats.py).

<a id="issue-6"></a>

### #6: something is wrong with bivariatePointNormal

[Original issue](https://github.com/bulik/ldsc/issues/6) · Source state: **closed** · Severity: **low** · Accessible comments reviewed: **4**.

Workflow tags: Cross-cutting / not assigned to a requested group.

**Trigger and harm.** Calling the developer bivariatePointNormal simulation helper returned invalid draws. A defective simulation generator can invalidate method-validation evidence.

**Legacy disposition.** Maintainer supplied a rewritten mixture sampler; local history contains 370a1bf7adab2d19b11bf2ba3ec358f2a27f746b, titled fixed bivariate point normal.

**LDSC3 assessment — Intentionally unsupported.** LDSC3 does not ship the historical simulate.py or bivariatePointNormal helper.
No production workflow invokes it, so the reported failure has no current execution path.
This does not certify new simulation code: any future correlated sparse-effect generator still needs covariance feasibility and marginal-moment checks.

**Recommended prevention.** Record as an out-of-scope developer-tool issue; validate moments and sparsity if a simulation helper is later added.

**Evidence and verification.** source-inspected
[src/ldsc/__init__.py](<../../../../src/ldsc/__init__.py>); [tests/package/test_package_layout.py](<../../../../tests/package/test_package_layout.py>); [Source](<https://github.com/bulik/ldsc/commit/370a1bf7adab2d19b11bf2ba3ec358f2a27f746b>).

<a id="issue-244"></a>

### #244: Warning message meaning : FutureWarning: `rcond` parameter will change to the default of machine precision times

[Original issue](https://github.com/bulik/ldsc/issues/244) · Source state: **open** · Severity: **low** · Accessible comments reviewed: **1**.

Workflow tags: `h2`, `rg`.

**Trigger and harm.** NumPy warns that the default least-squares rank cutoff will change.
Different dependency versions can change rank handling and clutter every run.

**Legacy disposition.** No fix attributed by this thread.

**LDSC3 assessment — Apparently prevented.** LDSC3 IRWLS explicitly calls np.linalg.lstsq(...,rcond=None), eliminating this version-default warning and making the chosen policy visible.
That choice does not prove equivalence with the historical rcond=-1 threshold for nearly singular designs, which remain a separate conditioning question.

**Recommended prevention.** Retain explicit rcond and conditioning tests; document the numerical policy rather than suppressing the warning globally.

**Evidence and verification.** source-inspected
[src/ldsc/_regression/irwls.py:IRWLS.wls](<../../../../src/ldsc/_regression/irwls.py>).

Legacy implementation context: [`ldscore/regressions.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/regressions.py), [`ldscore/sumstats.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/sumstats.py).

