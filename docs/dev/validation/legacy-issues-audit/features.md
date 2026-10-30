# Features and compatibility

Last updated on: 2026-10-28

[Audit overview](README.md) · [Workflow counts and membership](workflow-index.md) · [Verification](verification.md)

## Practical feature priorities

Existing capability and missing convenience are different findings.
LDSC3 already supports supplied numeric annotations, baseline-only joint partitioning, one-query-at-a-time fits, post-fit quantiles, saved delete values, directory artifacts and managed workers.
Those should be documented and verified before adding another interface.
An intentional unsupported interface is not automatically a numerical bug.

| Candidate | Current support and evidence | Practical recommendation |
| --- | --- | --- |
| BED interval overlap multiplicity | `ChromosomeProjector.project()` thresholds a cumulative overlap count to Boolean membership; V4 confirms overlaps remain one. | A single explicit count-projection option could retain the existing cumulative counts instead of thresholding. Cover duplicates, interval endpoints, SNP order and zero support; preserve binary membership as the default. |
| Continuous annotation interpretation | Numeric annotation input and `quantile-h2` already exist. | Make annotation type and interpretation visible beside saved results, so binary-style enrichment columns are not read without their caveat. |
| Extreme association-statistic input | Raw munging uses P to derive Z magnitude; direct munged input has different requirements. | First explain supported routes and diagnose underflow; consider a direct finite signed-statistic route only if needed, without inventing precision already lost in P. |
| Annotation-conditioned genetic covariance/correlation | Current `rg` is not a complete implementation of every requested stratified correlation method. | Treat this as a separate scientific extension with an appropriate estimator and validation, rather than expose partitioned columns as if already validated for that question. |
| Better retained-SNP and unestimable-result reporting | Existing funnels, query ledgers and pair statuses provide much of the needed structure. | Fill demonstrated gaps at those existing seams before adding a new tracking subsystem. |

The overlap-count proposal covers the use case described by legacy `make_annot.py --nomerge`.
It does not require a new storage format or external interval library.
The current `counts()` helper counts SNPs per interval; it is not the per-SNP number of overlapping intervals.
Supplied finite numeric annotations can already carry externally constructed counts, but their coefficients describe a per-unit count contribution and their interpretation differs from binary membership.
See [`interval_projection.py`](../../../../src/ldsc/_annotation/interval_projection.py), [`legacy-flag-map.md`](../../../reference/legacy-flag-map.md), and [verification V4](verification.md#v4-binary-projection-and-header-delimiters).

## Recommendations by kind

- **Code:** implement none of these proposals during this audit; the bugs report separates actual input/status defects from missing capabilities.
- **Documentation:** surface current equivalents and intentional limitations next to the affected workflow.
- **Features:** prioritize interval counts and demonstrated diagnostic gaps; defer broad new estimators until a concrete scientific use case and validation reference are specified.
- **Resources:** retain enough original annotation/reference metadata to reproduce any optional extension.

## Category summary

45 issues and 1 FAQ entry.

| Severity | Source | Concern | LDSC3 disposition | Evidence level |
| --- | --- | --- | --- | --- |
| High | [#20](#issue-20) | Single-locus domination of genetic correlation | Partial coverage | Source inspection |
| High | [#489](#issue-489) | Unsigned heritability input | Partial coverage | Source inspection |
| High | [#302](#issue-302) | Structural-variant heritability | Intentionally unsupported | Source inspection |
| High | [#314](#issue-314) | Subset-versus-parent enrichment contrast | Intentionally unsupported | Source inspection |
| High | [#340](#issue-340) | CNV heritability | Intentionally unsupported | Source inspection |
| High | [#353](#issue-353) | Nonhuman genome scope | Intentionally unsupported | Source inspection |
| High | [#364](#issue-364) | Custom baseline outside human scope | Intentionally unsupported | Source inspection |
| High | [#383](#issue-383) | Nonhuman chromosome count | Intentionally unsupported | Source inspection |
| High | [#411](#issue-411) | Partitioned genetic correlation | Intentionally unsupported | Source inspection |
| High | [#454](#issue-454) | Crop GWAS and single-cell extension | Intentionally unsupported | Source inspection |
| High | [#457](#issue-457) | Expression-trait LDSC extension | Intentionally unsupported | Source inspection |
| High | [#463](#issue-463) | Functional-group cross-trait analysis | Intentionally unsupported | Source inspection |
| High | [#469](#issue-469) | Local-region heritability | Intentionally unsupported | Source inspection |
| High | [#471](#issue-471) | Nonhuman regression weights | Intentionally unsupported | Source inspection |
| High | [#497](#issue-497) | Functional partition correlation | Intentionally unsupported | Source inspection |
| High | [#166](#issue-166) | Coordinate identifiers do not match rsID resources | Apparently prevented | Source inspection |
| Moderate | [#15](#issue-15) | Per-trait prevalence and intercept lists | Partial coverage | Source inspection |
| Moderate | [#68](#issue-68) | Partition heritability by chromosome | Partial coverage | Source inspection |
| Moderate | [#79](#issue-79) | Universal coordinate-to-rsID lookup | Partial coverage | Source inspection |
| Moderate | [#303](#issue-303) | Select a subset of annotation columns | Partial coverage | Source inspection |
| Moderate | [#394](#issue-394) | Continuous-score generation | Partial coverage | Source inspection |
| Moderate | [#395](#issue-395) | Single-cell tutorial | Partial coverage | Source inspection |
| Moderate | [#434](#issue-434) | P-value underflow reporting | Partial coverage | Source inspection |
| Moderate | [#448](#issue-448) | Reuse unchanged baseline LD scores | Partial coverage | Source inspection |
| Moderate | [#460](#issue-460) | Custom GTEx specificity scores | Partial coverage | Source inspection |
| Moderate | [#286](#issue-286) | Smallest sufficient regression-SNP set | Methodological limitation | Source inspection |
| Moderate | [#32](#issue-32) | Infer per-SNP sample size from effect precision | Intentionally unsupported | Source inspection |
| Moderate | [#107](#issue-107) | Indel-only LDSC analysis | Intentionally unsupported | Source inspection |
| Moderate | [#167](#issue-167) | Which regions drive an annotation enrichment | Intentionally unsupported | Source inspection |
| Moderate | [#172](#issue-172) | Chromosome X inference | Intentionally unsupported | Source inspection |
| Moderate | [#317](#issue-317) | Constrain one trait intercept only | Intentionally unsupported | Source inspection |
| Moderate | [#9](#issue-9) | Observed-to-liability conversion | Apparently prevented | Source inspection |
| Moderate | [#120](#issue-120) | Full per-query quantities beyond coefficient P | Apparently prevented | Source inspection |
| Moderate | [#122](#issue-122) | Uncontrolled native threading | Apparently prevented | Source inspection |
| Moderate | [#298](#issue-298) | All-versus-all trait correlations | Apparently prevented | Source inspection |
| Moderate | [#359](#issue-359) | Python-3 support | Apparently prevented | Source inspection |
| Moderate | [#396](#issue-396) | Separate annotation paths | Apparently prevented | Source inspection |
| Low | [#406](#issue-406) | CM/MAF output metadata | Partial coverage | Source inspection |
| Low | [#413](#issue-413) | Regression diagnostic plot | Partial coverage | Source inspection |
| Low | [#416](#issue-416) | Retained SNP export | Partial coverage | Source inspection |
| Low | [#491](#issue-491) | Per-row LD-window contributor count | Partial coverage | Source inspection |
| Low | [#188](#issue-188) | R implementation | Intentionally unsupported | Source inspection |
| Low | [FAQ-06](#faq-06) | Retired automatic binning failure | Intentionally unsupported | Source inspection |
| Low | [#5](#issue-5) | Compressed LD-score storage | Apparently prevented | Source inspection |
| Low | [#8](#issue-8) | Explicit summary-statistics filenames | Apparently prevented | Source inspection |
| Low | [#18](#issue-18) | LD-score output size and metadata columns | Apparently prevented | Source inspection |

## Individual assessments

<a id="issue-20"></a>

### #20: warn if rg driven by single locus

[Original issue](https://github.com/bulik/ldsc/issues/20) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **0**.

Workflow tags: `rg`.

**Trigger and harm.** Extreme Z1*Z2 observations drive a genome-wide correlation. A precise-looking genome-wide estimate may describe one pleiotropic locus or model misspecification.

**Legacy disposition.** Open request suggests checking extreme Z products; no implemented legacy fix is evidenced.

**LDSC3 assessment — Partial coverage.** LDSC3 accepts an explicit product-statistic chi-square cap and uses block jackknife uncertainty, but no automatic influence warning or leave-one-region summary was found in the correlation workflow. A cap changes the fitted target and should not be applied silently.
Block jackknife alone is not an interpretation guide for locus-driven estimates.

**Recommended prevention.** Add a simple documented sensitivity analysis excluding a prespecified dominant region and report its SNP count; consider an influence diagnostic only after validating it.

**Evidence and verification.** source-inspected
[src/ldsc/_regression/runner.py:_fit_rg_dataset](<../../../../src/ldsc/_regression/runner.py>); [src/ldsc/_regression/estimators.py:RG](<../../../../src/ldsc/_regression/estimators.py>); [docs/user-guide/workflows/genetic-correlation.md](<../../../user-guide/workflows/genetic-correlation.md>).

Legacy implementation context: [`ldscore/regressions.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/regressions.py), [`ldscore/sumstats.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/sumstats.py).

<a id="issue-489"></a>

### #489: Is it possible to run heritability and/or partitioned heritability with just p values and no signed-sumstats? And if yes - how?

[Original issue](https://github.com/bulik/ldsc/issues/489) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **1**.

Workflow tags: `munge-sumstats`, `h2`, `baseline-partitioned-h2`.

**Trigger and harm.** Has SNP identifiers and p-values but no direction and asks whether chi-square regression can use them. A workaround that invents positive directions can later contaminate genetic-correlation analyses.

**Legacy disposition.** A commenter says new code is needed and links heritability fitting; no dedicated unsigned route is implemented upstream.

**LDSC3 assessment — Partial coverage.** Univariate LDSC needs chi-square magnitudes rather than effect directions, but still needs valid N and suitable SNP coverage.
LDSC3 munger has an A1-increasing assertion, which is not a scientifically honest generic unsigned-data flag; rg requires real harmonized signs.

**Recommended prevention.** Consider an explicit unsigned provenance/magnitude input restricted to h2/partitioned-h2, with rg rejection; do not tell users to assert all A1 effects increase the trait just to bypass validation.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[src/ldsc/_sumstats/munge_kernel.py:741](<../../../../src/ldsc/_sumstats/munge_kernel.py#L741>); [src/ldsc/_sumstats/munger.py:1365](<../../../../src/ldsc/_sumstats/munger.py#L1365>); [src/ldsc/_regression/estimators.py:468](<../../../../src/ldsc/_regression/estimators.py#L468>).

Legacy implementation context: [`munge_sumstats.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/munge_sumstats.py), [`ldscore/regressions.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/regressions.py).

<a id="issue-302"></a>

### #302: LDSC for structural variants

[Original issue](https://github.com/bulik/ldsc/issues/302) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **1**.

Workflow tags: Cross-cutting / not assigned to a requested group.

**Trigger and harm.** Users ask about deletions, inversions and structural-variant datasets.
Applying SNP-specific LD/reference assumptions to structural variants may yield invalid inference.

**Legacy disposition.** Open with another user request, no validated method or implementation.

**LDSC3 assessment — Intentionally unsupported.** LDSC3 standard allele-aware handling and packaged references are designed for biallelic SNP workflows; parsing arbitrary labels in a base mode would not supply structural-variant LD or an appropriate statistical model.
This is a substantial method/validation extension, not an omitted input alias.

**Recommended prevention.** State the validated SNP scope and avoid advertising structural-variant support until appropriate references and scientific validation exist.

**Evidence and verification.** source-inspected
[src/ldsc/_core/snp_identity.py](<../../../../src/ldsc/_core/snp_identity.py>); [docs/reference/ldscore-snp-universes.md](<../../../reference/ldscore-snp-universes.md>); [src/ldsc/_regression/estimators.py:_filter_alleles](<../../../../src/ldsc/_regression/estimators.py>).

<a id="issue-314"></a>

### #314: Partitioned heritability - comparing enrichment between two categories

[Original issue](https://github.com/bulik/ldsc/issues/314) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **0**.

Workflow tags: `baseline-partitioned-h2`, `ldsc-seg`.

**Trigger and harm.** User wants a significance test comparing a subset with its parent annotation.
Comparing two enrichment p-values ignores covariance and does not test their difference.

**Legacy disposition.** Open; cites Hujoel2019and LDSC-SEG as possible approaches, without an implemented resolution.

**LDSC3 assessment — Intentionally unsupported.** LDSC3 supplies coefficient estimates/delete values and baseline/query conditional tests, but no direct arbitrary subset-versus-parent enrichment contrast was found. A coefficient beyond a parent control and an average-per-SNP heritability contrast answer related but distinct hypotheses; standardized coefficients are not automatically enrichment differences.

**Recommended prevention.** If needed, add one explicit covariance-aware block-jackknife contrast after defining subset versus complement within the parent; do not infer differences from separate significance.

**Evidence and verification.** source-inspected
[src/ldsc/_regression/runner.py:_coefficient_delete_frame](<../../../../src/ldsc/_regression/runner.py>); [src/ldsc/_regression/overlap_model.py](<../../../../src/ldsc/_regression/overlap_model.py>); [docs/reference/results/partitioned-h2.md](<../../../reference/results/partitioned-h2.md>).

Legacy implementation context: [`ldscore/regressions.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/regressions.py).

<a id="issue-340"></a>

### #340: Using LDSC for CNV-GWAS

[Original issue](https://github.com/bulik/ldsc/issues/340) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **0**.

Workflow tags: `h2`.

**Trigger and harm.** Requests heritability of copy-number variants from CNV-GWAS.
SNP LD and biallelic coding assumptions cannot be transferred automatically to structural variants.

**Legacy disposition.** Open; no answer or implementation fix in the accessible thread.

**LDSC3 assessment — Intentionally unsupported.** LDSC3 prepares biallelic SNP reference panels and a human HM3 regression set; accepting a summary-statistic table would not validate a CNV estimand.
This is an unsupported scientific extension, not an installation failure.

**Recommended prevention.** State the SNP scope prominently and require a separately validated LD/variant model before considering CNV support.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[src/ldsc/_sumstats/munge_kernel.py:232](<../../../../src/ldsc/_sumstats/munge_kernel.py#L232>); [docs/reference/snp-identity-and-genome-build.md](<../../../reference/snp-identity-and-genome-build.md>).

Legacy implementation context: [`ldscore/regressions.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/regressions.py).

<a id="issue-353"></a>

### #353: After merging with reference panel LD, 0 SNPs remain

[Original issue](https://github.com/bulik/ldsc/issues/353) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **1**.

Workflow tags: `h2`, `ldscore`.

**Trigger and harm.** Arabidopsis summary statistics use positional identifiers and only chromosome 2.
Human HM3 defaults and autosomal assumptions can make valid nonhuman data appear to have no overlap.

**Legacy disposition.** Open; a later comment asks for the source summary statistics but does not resolve the error.

**LDSC3 assessment — Intentionally unsupported.** LDSC3 is designed around human hg19/hg38 and the packaged human HM3 map; numeric SNP labels are not sufficient to establish nonhuman compatibility. A single chromosome also supplies limited information for genome-wide inference.

**Recommended prevention.** State the human-genome boundary and require explicit nonhuman chromosome/build/reference abstractions plus methodological validation before supporting this use.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[src/ldsc/_core/chromosomes.py](<../../../../src/ldsc/_core/chromosomes.py>); [src/ldsc/_genome/hm3.py](<../../../../src/ldsc/_genome/hm3.py>); [docs/reference/snp-identity-and-genome-build.md](<../../../reference/snp-identity-and-genome-build.md>).

Legacy implementation context: [`ldscore/ldscore.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/ldscore.py), [`ldscore/regressions.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/regressions.py).

<a id="issue-364"></a>

### #364: how to make my own baseline model?

[Original issue](https://github.com/bulik/ldsc/issues/364) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **1**.

Workflow tags: `baseline-partitioned-h2`.

**Trigger and harm.** Requests a custom baseline; the only comment specifies another species.
An arbitrary baseline or a human baseline transplanted to another species can produce uninterpretable partitioning.

**Legacy disposition.** Open with no construction recipe or verified fix.

**LDSC3 assessment — Intentionally unsupported.** LDSC3 accepts baseline annotation files but its human-genome and HM3 conventions remain.
File acceptance is separate from validating a baseline that controls relevant LD-dependent architecture in a new organism.

**Recommended prevention.** Document custom human baseline construction and scientific adequacy; describe nonhuman support as a separate unvalidated extension.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[src/ldsc/_annotation/builder.py](<../../../../src/ldsc/_annotation/builder.py>); [src/ldsc/_core/chromosomes.py](<../../../../src/ldsc/_core/chromosomes.py>); [docs/user-guide/workflows/partitioned-heritability.md](<../../../user-guide/workflows/partitioned-heritability.md>).

Legacy implementation context: [`ldscore/regressions.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/regressions.py).

<a id="issue-383"></a>

### #383: Could I set the number of chromosomes in LDSC?

[Original issue](https://github.com/bulik/ldsc/issues/383) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **2**.

Workflow tags: `ldscore`, `h2`.

**Trigger and harm.** Cattle analysis reads only chromosomes 1-22 although cattle has 29 autosomes.
Silently omitting or relabeling chromosomes changes the analysis and biological interpretation.

**Legacy disposition.** A commenter suggests concatenating files and --ref-ld; no methodological validation or current master fix is shown.

**LDSC3 assessment — Intentionally unsupported.** LDSC3 human normalization maps numeric 23/24 to X/Y and rejects other unsupported numeric labels.
It cannot safely acquire cattle support merely by changing a loop count or concatenating files.

**Recommended prevention.** State the human-only scope and reject recommendations to renumber nonhuman chromosomes; a species-aware extension requires explicit build and chromosome semantics.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[src/ldsc/_core/chromosomes.py:51](<../../../../src/ldsc/_core/chromosomes.py#L51>); [src/ldsc/_genome/hm3.py](<../../../../src/ldsc/_genome/hm3.py>).

Legacy implementation context: [`ldscore/ldscore.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/ldscore.py), [`ldscore/regressions.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/regressions.py).

<a id="issue-411"></a>

### #411: Partitioned genetic correlation

[Original issue](https://github.com/bulik/ldsc/issues/411) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **7**.

Workflow tags: `rg`, `baseline-partitioned-h2`.

**Trigger and harm.** Requests annotation-specific genetic correlation; later comments propose subsetting SNPs and running ordinary LDSC.
An easy-looking workaround changes the estimand and can produce misleading functional sharing claims.

**Legacy disposition.** Open; no validated implementation is provided, and the subset suggestion is not a demonstrated remedy.

**LDSC3 assessment — Intentionally unsupported.** LDSC3 has univariate partitioned heritability and ordinary genetic correlation, not a validated annotation-specific covariance/correlation workflow.
The mere existence of both commands does not imply their outputs can be combined.

**Recommended prevention.** State this limitation explicitly and evaluate a published stratified covariance method separately if demand justifies it; do not advertise SNP-filtered ordinary rg as equivalent.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[src/ldsc/_regression/runner.py:918](<../../../../src/ldsc/_regression/runner.py#L918>); [src/ldsc/_regression/runner.py:1233](<../../../../src/ldsc/_regression/runner.py#L1233>); [docs/reference/legacy-flag-map.md](<../../../reference/legacy-flag-map.md>).

Legacy implementation context: [`ldscore/regressions.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/regressions.py), [`ldscore/sumstats.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/sumstats.py).

<a id="issue-454"></a>

### #454: Other species

[Original issue](https://github.com/bulik/ldsc/issues/454) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **0**.

Workflow tags: `ldsc-seg`.

**Trigger and harm.** Requests wheat/other-species integration of GWAS and single-cell annotations.
Human genome-build inference, chromosome coding and baseline assumptions can silently misrepresent a new species.

**Legacy disposition.** Open with no implementation or validated recipe.

**LDSC3 assessment — Intentionally unsupported.** LDSC3 explicitly supports human hg19/hg38 and chromosome conventions, so data-format similarity does not establish crop applicability.
Polyploidy and species-specific LD further require methodological consideration.

**Recommended prevention.** State current human scope clearly; evaluate nonhuman support only with species-aware coordinates, reference SNP selection, genetic maps and a validated baseline model.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[src/ldsc/_core/chromosomes.py](<../../../../src/ldsc/_core/chromosomes.py>); [src/ldsc/_genome/build_inference.py](<../../../../src/ldsc/_genome/build_inference.py>); [docs/reference/snp-identity-and-genome-build.md](<../../../reference/snp-identity-and-genome-build.md>).

Legacy implementation context: [`ldscore/regressions.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/regressions.py).

<a id="issue-457"></a>

### #457: Extension of S-LDSC

[Original issue](https://github.com/bulik/ldsc/issues/457) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **0**.

Workflow tags: `baseline-partitioned-h2`, `rg`.

**Trigger and harm.** Requests the Liu et al. expression-summary-statistic extension and aggregate cross-tissue expression correlation.
Ordinary one-trait GWAS inputs cannot automatically reproduce a probe-aggregated or cis/trans-specific estimand.

**Legacy disposition.** Open; the requested published extension is named but no implementation or tutorial is supplied.

**LDSC3 assessment — Intentionally unsupported.** LDSC3 supports ordinary summary-statistic models, not the stated multi-probe extension.
Its presence in a methodological paper is not evidence that this software implements it or that a file-format adaptation suffices.

**Recommended prevention.** Document this boundary and inspect the original extension/code before prioritizing an implementation; require a clear expression estimand and a practical validated use case.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[src/ldsc/_regression/runner.py](<../../../../src/ldsc/_regression/runner.py>); [docs/reference/regression-options.md](<../../../reference/regression-options.md>); [Source](<https://github.com/bulik/ldsc/issues/457>).

Legacy implementation context: [`ldscore/regressions.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/regressions.py), [`ldscore/sumstats.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/sumstats.py).

<a id="issue-463"></a>

### #463: How to calculate stratified-LDSC between two traits in special function groups?

[Original issue](https://github.com/bulik/ldsc/issues/463) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **0**.

Workflow tags: `rg`, `baseline-partitioned-h2`.

**Trigger and harm.** Asks for stratified analysis of SCZ and BMI in 17 functional groups, referring to a paper screenshot.
An ordinary genome-wide rg command cannot establish annotation-specific genetic sharing.

**Legacy disposition.** Open; no executable method description or code fix is supplied.
The screenshot alone was not used as a verified methodological specification.

**LDSC3 assessment — Intentionally unsupported.** LDSC3 does not expose a validated stratified covariance/correlation workflow.
It can generate annotations and ordinary rg, but composing them by SNP filtering changes the statistical problem.

**Recommended prevention.** Request the exact cited method/estimand before implementing; explain the current scope and avoid presenting a shortcut as the paper reproduction.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[src/ldsc/_regression/runner.py:1233](<../../../../src/ldsc/_regression/runner.py#L1233>); [docs/reference/regression-options.md](<../../../reference/regression-options.md>).

Legacy implementation context: [`ldscore/regressions.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/regressions.py), [`ldscore/sumstats.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/sumstats.py).

<a id="issue-469"></a>

### #469: LDSC for each of the 1703 regions defined by ldetech (Joseph K Pickrell 2016)

[Original issue](https://github.com/bulik/ldsc/issues/469) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **0**.

Workflow tags: `ldscore`, `baseline-partitioned-h2`.

**Trigger and harm.** Requests separate heritability for 1703 LD blocks, while make_annot produces one union column.
Treating each small-region ordinary LDSC fit as local heritability can be statistically underidentified and poorly calibrated.

**Legacy disposition.** Open; no validated implementation or recipe.

**LDSC3 assessment — Intentionally unsupported.** LDSC3 BED projection makes one binary query per source, not one column per interval by default, and joint region indicators plus an all-ones base can create dependence.
Local heritability is a separate estimand from ordinary genome-wide enrichment.

**Recommended prevention.** Explain the difference between a union annotation, per-region query models and a validated local-heritability method; do not add 1703 regressions as an automatic workaround.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[src/ldsc/_annotation/interval_projection.py](<../../../../src/ldsc/_annotation/interval_projection.py>); [src/ldsc/_annotation/interval_queries.py](<../../../../src/ldsc/_annotation/interval_queries.py>); [docs/reference/results/partitioned-h2.md](<../../../reference/results/partitioned-h2.md>).

Legacy implementation context: [`ldscore/ldscore.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/ldscore.py), [`ldscore/regressions.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/regressions.py).

<a id="issue-471"></a>

### #471: how to build my own input weight files for --w-ld-chr

[Original issue](https://github.com/bulik/ldsc/issues/471) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **3**.

Workflow tags: `ldscore`.

**Trigger and harm.** Requests weight files for another species and follows up that the whole reference remains unavailable.
Creating a weight file alone does not make the human defaults valid for a new organism.

**Legacy disposition.** Comments point to the weight tutorial and common variants/unrelated reference genotypes; no nonhuman validation is established.

**LDSC3 assessment — Intentionally unsupported.** LDSC3 human chromosome/build/HM3 policies prevent a simple substitution.
It can compute human reference and regression SNP weights but does not claim species-independent scientific applicability.

**Recommended prevention.** Document required scientific inputs and current human scope; defer a species extension until chromosome, build, LD and variant-selection assumptions are explicit.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[src/ldsc/_core/chromosomes.py](<../../../../src/ldsc/_core/chromosomes.py>); [src/ldsc/_genome/hm3.py](<../../../../src/ldsc/_genome/hm3.py>); [docs/reference/ldscore-snp-universes.md](<../../../reference/ldscore-snp-universes.md>).

Legacy implementation context: [`ldscore/ldscore.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/ldscore.py).

<a id="issue-497"></a>

### #497: how to do Partitioned genetic correlation analysis

[Original issue](https://github.com/bulik/ldsc/issues/497) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **0**.

Workflow tags: `rg`, `baseline-partitioned-h2`.

**Trigger and harm.** Requests genetic correlations across named functional annotations.
Annotation-specific heritability and ordinary cross-trait correlation cannot be combined into a validated stratified correlation by command composition.

**Legacy disposition.** Open with no supported implementation or source-method details.

**LDSC3 assessment — Intentionally unsupported.** LDSC3 does not expose stratified genetic covariance/correlation and has no claim of calibrated inference for this estimand. A separate published method and its uncertainty construction are required.

**Recommended prevention.** Document the boundary and collect the precise intended estimand and published reference before prioritizing; avoid an ordinary rg run on selected annotation SNPs.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[src/ldsc/_regression/runner.py:1233](<../../../../src/ldsc/_regression/runner.py#L1233>); [docs/reference/regression-options.md](<../../../reference/regression-options.md>).

Legacy implementation context: [`ldscore/regressions.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/regressions.py), [`ldscore/sumstats.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/sumstats.py).

<a id="issue-166"></a>

### #166: LD-sc using chr:bp instead of rs numbers

[Original issue](https://github.com/bulik/ldsc/issues/166) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **0**.

Workflow tags: `munge-sumstats`, `h2`.

**Trigger and harm.** HRC-imputed variants use chr:bp while a keep list and LD scores use rsIDs.
All rows disappear or a user bypasses matching without a coherent identity contract.

**Legacy disposition.** Open without a legacy implementation or fix.

**LDSC3 assessment — Apparently prevented.** LDSC3 supports coordinate-family SNP identifier modes, hg19/hg38 declarations/inference, HM3 mapping and explicit legacy-summary-statistics projection onto panel coordinates.
Merely naming a chr:bp column rs does not convert it to rsID; coordinate build and alleles must agree.
The pasted command also contains a stray comma in a path, a separate usage error.

**Recommended prevention.** Show a complete coordinate-family workflow with correct genome build and matching LD-score directory; retain explicit no-overlap diagnostics rather than guessing identities.

**Evidence and verification.** source-inspected
[docs/reference/snp-identity-and-genome-build.md](<../../../reference/snp-identity-and-genome-build.md>); [src/ldsc/_genome/hm3.py](<../../../../src/ldsc/_genome/hm3.py>); [src/ldsc/_regression/runner.py:_project_legacy_sumstats_to_panel](<../../../../src/ldsc/_regression/runner.py>).

Legacy implementation context: [`munge_sumstats.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/munge_sumstats.py), [`ldscore/regressions.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/regressions.py).

<a id="issue-15"></a>

### #15: setting --constrain-intercept / --samp-prev / --pop-prev with --rg and a list

[Original issue](https://github.com/bulik/ldsc/issues/15) · Source state: **closed** · Severity: **moderate** · Accessible comments reviewed: **1**.

Workflow tags: `rg`.

**Trigger and harm.** A multi-trait correlation run needs different settings for each trait or pair.
Using one prevalence or intercept indiscriminately can change scales or covariance adjustment.

**Legacy disposition.** Maintainer confirmed multiargument support in dev/v1.0.0. A later referenced commit is not independently treated as the fix.

**LDSC3 assessment — Partial coverage.** LDSC3 supports per-trait prevalence lists or a manifest and all-pairs or anchor correlation runs.
Current intercept-h2 and intercept-gencov CLI controls are scalar, so arbitrary legacy lists of trait-specific/pair-specific intercepts are not reproduced.
Users can run separate fits, but that workaround should be explicit rather than silently reusing a scalar.

**Recommended prevention.** Document the scalar-intercept boundary and recommend separate explicit fits when intercept assumptions differ; add a richer interface only for demonstrated demand.

**Evidence and verification.** source-inspected
[src/ldsc/_regression/prevalence.py](<../../../../src/ldsc/_regression/prevalence.py>); [src/ldsc/_regression/runner.py:add_rg_arguments](<../../../../src/ldsc/_regression/runner.py>); [docs/reference/regression-options.md](<../../../reference/regression-options.md>).

Legacy implementation context: [`ldscore/regressions.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/regressions.py), [`ldscore/sumstats.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/sumstats.py).

<a id="issue-68"></a>

### #68: How can I get the partitioned heritability by chromosome 1-22 other than functional category?

[Original issue](https://github.com/bulik/ldsc/issues/68) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **0**.

Workflow tags: `baseline-partitioned-h2`.

**Trigger and harm.** User requests chromosome-level heritability and enrichment rather than functional annotations.
Running 22 isolated regressions can be confused with jointly partitioning a genome-wide target.

**Legacy disposition.** Open with no maintainer answer or verified implementation.

**LDSC3 assessment — Partial coverage.** LDSC3 consumes arbitrary prebuilt annotation columns, so one binary column per chromosome could define a joint baseline partition subject to rank/support checks.
It does not provide a dedicated chromosome-enrichment command or certify that narrow chromosome-specific fits meet genome-wide LDSC assumptions.
Query mode would instead fit separate models and answer another question.

**Recommended prevention.** Document the distinction and assess estimability before adding convenience generation; avoid a new command unless a validated chromosome-partition use case is needed.

**Evidence and verification.** source-inspected
[src/ldsc/_annotation/parsing.py](<../../../../src/ldsc/_annotation/parsing.py>); [src/ldsc/_regression/runner.py:_fit_functional_model](<../../../../src/ldsc/_regression/runner.py>); [docs/user-guide/workflows/partitioned-heritability.md](<../../../user-guide/workflows/partitioned-heritability.md>).

Legacy implementation context: [`ldscore/regressions.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/regressions.py).

<a id="issue-79"></a>

### #79: Feature Request

[Original issue](https://github.com/bulik/ldsc/issues/79) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **2**.

Workflow tags: `munge-sumstats`.

**Trigger and harm.** Decorated coordinate names are unusable in rsID-only downstream software.
Naive position-only mapping can select another build, allele or multiallelic record.

**Legacy disposition.** Maintainer declines a universal matcher because of identity nuance and limited LDSC need.

**LDSC3 assessment — Partial coverage.** LDSC3 supports coordinate-family modes and an HM3 dual-build map, so LDSC analyses need not first obtain rsIDs.
It does not offer a universal dbSNP/1000G annotation service for every external application.
The narrow bundled mapping and explicit allele-aware identity are preferable to silently guessing rsIDs from arbitrary strings.

**Recommended prevention.** Document coordinate-family workflows and refer broad identifier annotation to dedicated validated tools; keep ambiguous mappings visible.

**Evidence and verification.** source-inspected
[src/ldsc/_core/snp_identity.py](<../../../../src/ldsc/_core/snp_identity.py>); [src/ldsc/_genome/hm3.py](<../../../../src/ldsc/_genome/hm3.py>); [docs/reference/snp-identity-and-genome-build.md](<../../../reference/snp-identity-and-genome-build.md>).

Legacy implementation context: [`munge_sumstats.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/munge_sumstats.py).

<a id="issue-303"></a>

### #303: specifying columns in annotation file

[Original issue](https://github.com/bulik/ldsc/issues/303) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **0**.

Workflow tags: `ldscore`, `baseline-partitioned-h2`.

**Trigger and harm.** User wants selected columns from a many-column annotation file without making separate files.
Extra annotations change conditioning and increase computation/collinearity.

**Legacy disposition.** Open without implementation or workaround.

**LDSC3 assessment — Partial coverage.** LDSC3 Python regression APIs can select query_columns, but baseline columns in a directory define the common fitted baseline and no general public baseline-column selector was found.
Users can prepare a deliberate subset annotation table while preserving metadata and recompute counts/overlap consistently; editing only LD-score columns is unsafe.

**Recommended prevention.** Consider one simple explicit column-selection option if repeated demand warrants it; first document safe input-column preparation and resulting model changes.

**Evidence and verification.** source-inspected
[src/ldsc/_regression/runner.py:estimate_partitioned_h2_batch](<../../../../src/ldsc/_regression/runner.py>); [src/ldsc/_annotation/sources.py:prepare_annotation_sources](<../../../../src/ldsc/_annotation/sources.py>); [docs/user-guide/workflows/annotations.md](<../../../user-guide/workflows/annotations.md>).

Legacy implementation context: [`make_annot.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/make_annot.py), [`ldscore/regressions.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/regressions.py).

<a id="issue-394"></a>

### #394: s-LDSC continuous annotations

[Original issue](https://github.com/bulik/ldsc/issues/394) · Source state: **closed** · Severity: **moderate** · Accessible comments reviewed: **2**.

Workflow tags: `ldscore`, `ldsc-seg`, `continuous`.

**Trigger and harm.** Requests BED intervals with values to create continuous annotations.
Binary interval membership loses the supplied quantitative score or overlap multiplicity.

**Legacy disposition.** Author closes after obtaining scores from an external database; no generation feature is added.

**LDSC3 assessment — Partial coverage.** LDSC3 consumes continuous annotation columns and provides coefficient/quantile summaries, but BED projection is binary interval union. A BED fourth column must not be assumed to become a per-SNP quantitative value.

**Recommended prevention.** Document how externally derived values map to SNPs and distinguish binary union, overlap counts and weighted interval aggregation; add only the explicitly needed projection mode.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[src/ldsc/_annotation/interval_projection.py:15](<../../../../src/ldsc/_annotation/interval_projection.py#L15>); [docs/reference/legacy-flag-map.md:150](<../../../reference/legacy-flag-map.md#legacy-ldsc-options-with-no-current-equivalent>); [docs/user-guide/workflows/continuous-annotations.md](<../../../user-guide/workflows/continuous-annotations.md>).

Legacy implementation context: [`ldscore/ldscore.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/ldscore.py), [`ldscore/regressions.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/regressions.py).

<a id="issue-395"></a>

### #395: tutorial request

[Original issue](https://github.com/bulik/ldsc/issues/395) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **2**.

Workflow tags: `ldsc-seg`.

**Trigger and harm.** Requests an end-to-end analysis from custom single-cell data and downloaded GWAS.
The gap between biological data and a valid query/control annotation makes users copy incomplete pipelines.

**Legacy disposition.** Open; later comments repeat the request and link scLDSC, without a canonical tutorial.

**LDSC3 assessment — Partial coverage.** LDSC3 has gene-list, BED and cell-type enrichment workflows, but does not perform expression normalization or define a scientifically justified cell-specific score automatically.

**Recommended prevention.** Add a compact worked example separating external cell-type score construction, query selection, control genes, LD scores and the coefficient test; avoid bundling a whole single-cell framework.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[docs/user-guide/workflows/cell-type-enrichment.md](<../../../user-guide/workflows/cell-type-enrichment.md>); [src/ldsc/_annotation/gene_lists.py](<../../../../src/ldsc/_annotation/gene_lists.py>); [docs/user-guide/guided-tutorial.md](<../../../user-guide/guided-tutorial.md>).

Legacy implementation context: [`ldscore/regressions.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/regressions.py).

<a id="issue-434"></a>

### #434: Lowest limit of rg p-value reported?

[Original issue](https://github.com/bulik/ldsc/issues/434) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **6**.

Workflow tags: `rg`.

**Trigger and harm.** A very large genetic-correlation test statistic is reported with p=0.
Users may report an impossible exact zero or an arbitrary lower threshold.

**Legacy disposition.** Commenter recommends reporting interval or Z and clarifies the parenthesized quantity is jackknife SE; no log-p feature is supplied.

**LDSC3 assessment — Partial coverage.** LDSC3 still computes p through chi2.sf(Z**2,1), which can underflow while Z and SE remain finite.
This is a numerical representation limit, not proof that a high same-trait correlation is erroneous.

**Recommended prevention.** Keep estimate/SE/Z available and label underflow; an optional log-p from a stable log-survival function is a small useful addition if exact tail reporting is needed.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[src/ldsc/_regression/estimators.py:65](<../../../../src/ldsc/_regression/estimators.py#L65>); [docs/user-guide/workflows/genetic-correlation.md](<../../../user-guide/workflows/genetic-correlation.md>).

Legacy implementation context: [`ldscore/regressions.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/regressions.py), [`ldscore/sumstats.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/sumstats.py).

<a id="issue-448"></a>

### #448: How to update the partitioned LDscore?

[Original issue](https://github.com/bulik/ldsc/issues/448) · Source state: **closed** · Severity: **moderate** · Accessible comments reviewed: **2**.

Workflow tags: `ldscore`, `ldsc-seg`, `continuous`.

**Trigger and harm.** Asks to reuse baseline columns when replacing a query annotation to avoid days of recomputation.
Repeated baseline work wastes compute, but blind column pasting can desynchronize identities, counts and overlaps.

**Legacy disposition.** Commenter approves appending legacy columns; author acknowledges, without a compatibility demonstration.

**LDSC3 assessment — Partial coverage.** LDSC3 supports query batching, reusable R2 matrices and gene indexes, but its directory metadata binds additional count/overlap artifacts.
Manually editing one Parquet column is not a complete supported append operation.

**Recommended prevention.** Recommend supported batching/index reuse first; consider a simple validated append only if demand remains, requiring matching panel/window/identity and recomputed query overlap metadata.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[docs/dev/design/memory-and-cpu.md](<../../design/memory-and-cpu.md>); [docs/reference/gene-ldscore-index.md](<../../../reference/gene-ldscore-index.md>); [src/ldsc/_ldscore/outputs.py](<../../../../src/ldsc/_ldscore/outputs.py>).

Legacy implementation context: [`ldscore/ldscore.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/ldscore.py), [`ldscore/regressions.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/regressions.py).

<a id="issue-460"></a>

### #460: customized expression data for cell specific analysis

[Original issue](https://github.com/bulik/ldsc/issues/460) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **0**.

Workflow tags: `ldsc-seg`.

**Trigger and harm.** Requests code for t-statistic ranking of tissue-specific genes in GTEx v8.
An annotation can be syntactically valid yet represent a different specificity contrast than the published method.

**Legacy disposition.** Open; no score-construction code or validated updated annotation is supplied.

**LDSC3 assessment — Partial coverage.** LDSC3 can turn a supplied gene list into a query and include an all-gene control, but it does not fit the expression-specificity model.
That preprocessing remains a scientifically consequential external step.

**Recommended prevention.** Document the required output and method-specific ranking/control choices, with a small external-preprocessing example rather than silently embedding an expression-analysis pipeline.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[src/ldsc/_annotation/gene_lists.py](<../../../../src/ldsc/_annotation/gene_lists.py>); [docs/user-guide/workflows/cell-type-enrichment.md](<../../../user-guide/workflows/cell-type-enrichment.md>).

Legacy implementation context: [`ldscore/regressions.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/regressions.py).

<a id="issue-286"></a>

### #286: Minimal set of rsIDs for heritability estimation (not an issue)

[Original issue](https://github.com/bulik/ldsc/issues/286) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **0**.

Workflow tags: `h2`.

**Trigger and harm.** A user wants the minimal rsID list permitting accurate SNP-heritability estimation.
Aggressive thinning can reduce precision or create selection bias.

**Legacy disposition.** Open without a validated universal list or feature.

**LDSC3 assessment — Methodological limitation.** LDSC3 allows explicit regression-SNP selection but retains a broad reference-SNP denominator.
The minimum useful set depends on trait signal, LD coverage, ancestry and QC; passing a fixed row threshold does not ensure calibrated inference.
HM3 is a practical default rather than a proof of minimality.

**Recommended prevention.** Document the information/precision tradeoff and avoid presenting an unsupported minimal SNP count as a guarantee.

**Evidence and verification.** source-inspected
[docs/reference/ldscore-snp-universes.md](<../../../reference/ldscore-snp-universes.md>); [docs/reference/regression-options.md](<../../../reference/regression-options.md>).

Legacy implementation context: [`ldscore/regressions.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/regressions.py).

<a id="issue-32"></a>

### #32: infer N from beta [SE or P] and reference panel MAF

[Original issue](https://github.com/bulik/ldsc/issues/32) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **2**.

Workflow tags: `munge-sumstats`.

**Trigger and harm.** A publication gives maximum N but no per-SNP N, beta precision and reference MAF might permit approximation.
Using maximum N for all SNPs can bias scaling when coverage differs.

**Legacy disposition.** Open proposal; discussion calls reference-MAF inference an approximation, not a validated universal estimator.

**LDSC3 assessment — Intentionally unsupported.** LDSC3 supports supplied per-SNP and constant sample sizes with explicit precedence, but does not infer N from beta/SE and external MAF.
Such inference depends on phenotype scale, imputation, model and allele-frequency agreement.
Automatically adding it would replace missing evidence with assumptions and can be particularly wrong for mixed models or case-control meta-analysis.

**Recommended prevention.** Prefer requesting released per-SNP N or documenting a sensitivity analysis; keep any derived-N helper separate and assumption-explicit.

**Evidence and verification.** source-inspected
[src/ldsc/_sumstats/munge_input.py:resolve_sample_size](<../../../../src/ldsc/_sumstats/munge_input.py>); [src/ldsc/_sumstats/munge_kernel.py:process_n](<../../../../src/ldsc/_sumstats/munge_kernel.py>); [docs/reference/munge-sumstats.md](<../../../reference/munge-sumstats.md>).

Legacy implementation context: [`munge_sumstats.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/munge_sumstats.py).

<a id="issue-107"></a>

### #107: LD scores for indels from 1000 genome

[Original issue](https://github.com/bulik/ldsc/issues/107) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **0**.

Workflow tags: `ldscore`, `munge-sumstats`.

**Trigger and harm.** User wants LD scores for indels and asks why munging removes them. A technically readable genotype panel is not evidence that indel-only genome-wide regression is calibrated.

**Legacy disposition.** Open without a maintainer answer or verified extension.

**LDSC3 assessment — Intentionally unsupported.** LDSC3 summary-statistics allele filtering targets strand-unambiguous single-nucleotide variants, while PLINK low-level scoring and annotation inputs have separate representational constraints.
An indel-only analysis is not an advertised validated workflow.
Fewer, nonrepresentative variants and imputation quality affect method assumptions beyond parsing.

**Recommended prevention.** State the supported variant scope clearly; require representation harmonization and method validation before adding indel support rather than removing one filter.

**Evidence and verification.** source-inspected
[src/ldsc/_sumstats/munge_kernel.py:filter_alleles](<../../../../src/ldsc/_sumstats/munge_kernel.py>); [src/ldsc/_core/snp_identity.py](<../../../../src/ldsc/_core/snp_identity.py>); [docs/reference/munge-sumstats.md](<../../../reference/munge-sumstats.md>).

Legacy implementation context: [`munge_sumstats.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/munge_sumstats.py), [`ldscore/ldscore.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/ldscore.py).

<a id="issue-167"></a>

### #167: Estimate regions that contribute to enrichment

[Original issue](https://github.com/bulik/ldsc/issues/167) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **0**.

Workflow tags: `ldsc-seg`.

**Trigger and harm.** User wants genomic localization behind a cell-type enrichment result. A genome-wide annotation association can be misrepresented as identifying causal loci.

**Legacy disposition.** Open with no validated legacy localization feature.

**LDSC3 assessment — Intentionally unsupported.** LDSC3 can inspect annotation membership, LD-score rows and aggregate model outputs, but it does not implement a calibrated local-heritability or causal-region decomposition.
Jackknife deletion effects measure influence, not per-region biological heritability. A simple region-exclusion sensitivity is useful if reported as changing the target.

**Recommended prevention.** Offer transparent influence/sensitivity summaries only with their limits; refer local genetic architecture questions to separately validated methods rather than invent a causal ranking.

**Evidence and verification.** source-inspected
[src/ldsc/_regression/jackknife.py](<../../../../src/ldsc/_regression/jackknife.py>); [docs/user-guide/workflows/cell-type-enrichment.md](<../../../user-guide/workflows/cell-type-enrichment.md>); [docs/reference/results/partitioned-h2.md](<../../../reference/results/partitioned-h2.md>).

Legacy implementation context: [`ldscore/regressions.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/regressions.py).

<a id="issue-172"></a>

### #172: X chromosome LD score regression

[Original issue](https://github.com/bulik/ldsc/issues/172) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **1**.

Workflow tags: `h2`, `rg`.

**Trigger and harm.** User asks whether chromosome X permits polygenicity, SNP heritability and genetic-correlation analyses.
Sex-specific dosage, male hemizygosity and appropriate denominators require more than accepting an X label.

**Legacy disposition.** Open; comments speculate about female-only analysis but supply neither validated support nor a fix.

**LDSC3 assessment — Intentionally unsupported.** LDSC3 normalizes X labels, but LD-score preflight and standard chromosome suites target autosomes; this is not evidence for a validated X-chromosome estimator.
The proposed female-only workaround remains a scientific assumption, not an implemented capability.

**Recommended prevention.** State autosomal scope clearly; consider X only after a dedicated sex-aware methodological validation.

**Evidence and verification.** source-inspected
[src/ldsc/_core/chromosomes.py:normalize_chromosome](<../../../../src/ldsc/_core/chromosomes.py>); [src/ldsc/_ldscore/preflight.py:130](<../../../../src/ldsc/_ldscore/preflight.py#L130>); [docs/reference/input-paths.md](<../../../reference/input-paths.md>).

Legacy implementation context: [`ldscore/regressions.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/regressions.py), [`ldscore/sumstats.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/sumstats.py).

<a id="issue-317"></a>

### #317: rg –intercept-h2: constrain the intercept of a single sumstat

[Original issue](https://github.com/bulik/ldsc/issues/317) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **0**.

Workflow tags: `rg`.

**Trigger and harm.** A user wants one SNP-heritability intercept fixed and the other estimated freely.
The same constraint applied to both traits can encode an unintended assumption.

**Legacy disposition.** Open with no supported legacy resolution.

**LDSC3 assessment — Intentionally unsupported.** The LDSC3 RG kernel has separate intercept_hsq1/intercept_hsq2 parameters, but the public runner selects one scalar intercept_h2 and passes it to both.
Thus the scientific kernel can express the requested fit while the public CLI/API configuration cannot.
Mixed constraints also require an explicit compatible two-step policy.

**Recommended prevention.** Consider a simple per-trait intercept configuration only if justified; retain free intercept defaults and document estimator compatibility.

**Evidence and verification.** source-inspected
[src/ldsc/_regression/runner.py:_fit_rg_dataset](<../../../../src/ldsc/_regression/runner.py>); [src/ldsc/_regression/runner.py:_select_intercept](<../../../../src/ldsc/_regression/runner.py>); [src/ldsc/_regression/estimators.py:RG](<../../../../src/ldsc/_regression/estimators.py>).

Legacy implementation context: [`ldscore/regressions.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/regressions.py), [`ldscore/sumstats.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/sumstats.py).

<a id="issue-9"></a>

### #9: Feature request

[Original issue](https://github.com/bulik/ldsc/issues/9) · Source state: **closed** · Severity: **moderate** · Accessible comments reviewed: **3**.

Workflow tags: `h2`, `rg`.

**Trigger and harm.** A binary-trait estimate may be reported on the wrong scale.
Scientifically incompatible heritability values can be compared or published as the same quantity.

**Legacy disposition.** Maintainer confirmed --samp-prev/--pop-prev support in v1.0.0 and documented nan entries for quantitative traits.

**LDSC3 assessment — Apparently prevented.** LDSC3 supports prevalence-aware regression outputs and post-fit convert-h2-scale.
The dedicated post-fit command consumes a canonical h2 result rather than the arbitrary multirow text file first requested.
It saves observed and liability quantities explicitly; per-trait prevalence manifests support correlation analyses with mixed trait types.

**Recommended prevention.** Make scale labels mandatory in worked interpretations and show mixed binary/quantitative prevalence inputs; no new general-purpose table converter is needed.

**Evidence and verification.** source-inspected
[src/ldsc/_regression/prevalence.py](<../../../../src/ldsc/_regression/prevalence.py>); [src/ldsc/_regression/h2_scale.py](<../../../../src/ldsc/_regression/h2_scale.py>); [docs/user-guide/workflows/liability-scale.md](<../../../user-guide/workflows/liability-scale.md>).

Legacy implementation context: [`ldscore/regressions.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/regressions.py), [`ldscore/sumstats.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/sumstats.py).

<a id="issue-120"></a>

### #120: Calculating additional information from CTS analysis

[Original issue](https://github.com/bulik/ldsc/issues/120) · Source state: **closed** · Severity: **moderate** · Accessible comments reviewed: **2**.

Workflow tags: `ldsc-seg`.

**Trigger and harm.** User wants normalized coefficient, proportion of heritability and SNP share from a cell-type scan.
An aggregate query table may combine quantities from different fitted models.

**Legacy disposition.** Maintainer recommends separate regressions because legacy h2-cts does not load annotations for overlap summaries.

**LDSC3 assessment — Apparently prevented.** LDSC3 stores full per-query baseline-plus-one-query fits with overlap-aware summaries and float64 coefficient delete values. quantile-h2 can compute standardized coefficients from one complete fitted model.
Its aggregate root is not one joint model across all queries, so downstream summaries must select a complete query directory.

**Recommended prevention.** Document the per-query route and distinguish normalized/standardized coefficients from overlap-aware proportions; reuse saved fits rather than add redundant refitting.

**Evidence and verification.** source-inspected
[docs/reference/results/partitioned-h2.md:65-74](<../../../reference/results/partitioned-h2.md#partitioned_h2tsv>); [docs/user-guide/workflows/continuous-annotations.md:20-41](<../../../user-guide/workflows/continuous-annotations.md#1-preserve-one-complete-fitted-model>); [src/ldsc/_regression/quantile_h2.py](<../../../../src/ldsc/_regression/quantile_h2.py>).

Legacy implementation context: [`ldscore/regressions.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/regressions.py).

<a id="issue-122"></a>

### #122: Threads used when estimating LD scores

[Original issue](https://github.com/bulik/ldsc/issues/122) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **1**.

Workflow tags: `ldscore`.

**Trigger and harm.** NumPy/BLAS uses all available CPU cores during one score calculation. A nominally single-task run can oversubscribe shared machines or scheduler allocations.

**Legacy disposition.** Maintainer points to BLAS environment controls, not Python-level parallelism; no legacy CLI fix is identified.

**LDSC3 assessment — Apparently prevented.** LDSC3 separates managed worker count from native numerical threads and enforces one BLAS thread per worker through its runtime helpers.
CPU detection and task caps are documented. A worker budget is still not a guarantee about every external library or total memory, and selected source runtime checks matter on the execution platform.

**Recommended prevention.** Keep one clear workers option and actual detected/used CPU diagnostics; do not add several overlapping threading flags.

**Evidence and verification.** source-inspected
[src/ldsc/_runtime/worker_runtime.py](<../../../../src/ldsc/_runtime/worker_runtime.py>); [src/ldsc/_runtime/serial_calculation.py](<../../../../src/ldsc/_runtime/serial_calculation.py>); [src/ldsc/_runtime/parallelism.py](<../../../../src/ldsc/_runtime/parallelism.py>); [docs/reference/parallel-execution.md](<../../../reference/parallel-execution.md>).

Legacy implementation context: [`ldscore/ldscore.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/ldscore.py).

<a id="issue-298"></a>

### #298: rg all vs all

[Original issue](https://github.com/bulik/ldsc/issues/298) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **6**.

Workflow tags: `rg`.

**Trigger and harm.** Eleven traits require every pair, while legacy --rg fits only against the first.
Shell loops and manual bookkeeping create duplicate/missing comparisons and poor onboarding.

**Legacy disposition.** Reporter supplies nested shell loops; no legacy feature merge identified.

**LDSC3 assessment — Apparently prevented.** LDSC3 implements all-pairs fitting by default plus --anchor-trait for the legacy anchored use case, with per-pair status and scientific TSV outputs.
This directly supplies a practical requested capability without a new orchestration framework.

**Recommended prevention.** Lead with the all-pairs example and report the number of unique pairs; keep failure notes and multiplicity guidance visible.

**Evidence and verification.** source-inspected
[src/ldsc/_regression/runner.py:_iter_rg_pairs](<../../../../src/ldsc/_regression/runner.py>); [docs/user-guide/workflows/genetic-correlation.md](<../../../user-guide/workflows/genetic-correlation.md>).

Legacy implementation context: [`ldscore/regressions.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/regressions.py), [`ldscore/sumstats.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/sumstats.py).

<a id="issue-359"></a>

### #359: LDSC is not available for Python 3

[Original issue](https://github.com/bulik/ldsc/issues/359) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **1**.

Workflow tags: Cross-cutting / not assigned to a requested group.

**Trigger and harm.** Requests a supported Python-3 version; a later comment links another implementation.
Users otherwise rely on obsolete Python or untracked forks.

**Legacy disposition.** Open; the linked cbiit project is an alternative, not evidence of a bulik/ldsc master fix.

**LDSC3 assessment — Apparently prevented.** LDSC3 is explicitly a Python-3 package with modern dependency bounds and one CLI.
That satisfies the requested platform capability but does not establish numerical equivalence with all third-party ports.

**Recommended prevention.** Make the development branch and distribution identity unmistakable, and document supported interpreter versions.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[setup.py](<../../../../setup.py>); [src/ldsc/__main__.py](<../../../../src/ldsc/__main__.py>); [docs/user-guide/installation.md](<../../../user-guide/installation.md>).

<a id="issue-396"></a>

### #396: --overlap-annot does not recognize path in --annot

[Original issue](https://github.com/bulik/ldsc/issues/396) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **0**.

Workflow tags: `baseline-partitioned-h2`.

**Trigger and harm.** Legacy overlap analysis ignores --annot as a separate location and searches beside LD-score files. A fit can complete then fail when annotation files live elsewhere.

**Legacy disposition.** Open; the quoted _read_annot implementation uses reference prefixes, and no accepted fix is shown.

**LDSC3 assessment — Apparently prevented.** LDSC3 computes and stores the overlap matrix in the LD-score directory before regression.
Native regression no longer rediscovers raw annotation files by a shared prefix; conversion exposes explicit source inputs.

**Recommended prevention.** Explain artifact contents and make missing overlap data a pre-fit validation error, preserving separate annotation source paths in conversion diagnostics.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[src/ldsc/_ldscore/outputs.py](<../../../../src/ldsc/_ldscore/outputs.py>); [src/ldsc/_regression/overlap_model.py](<../../../../src/ldsc/_regression/overlap_model.py>); [src/ldsc/_ldscore/legacy_converter.py](<../../../../src/ldsc/_ldscore/legacy_converter.py>).

Legacy implementation context: [`ldscore/regressions.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/regressions.py).

<a id="issue-406"></a>

### #406: Add CM and MAF file to LD Score output .gz

[Original issue](https://github.com/bulik/ldsc/issues/406) · Source state: **open** · Severity: **low** · Accessible comments reviewed: **0**.

Workflow tags: `ldscore`.

**Trigger and harm.** Requests genetic distance and allele frequency alongside LD-score output.
Missing metadata complicates alignment checks and follow-up analyses.

**Legacy disposition.** Open; no implementation answer or fix.

**LDSC3 assessment — Partial coverage.** LDSC3 uses Parquet outputs and reference-derived metadata rather than an extensible legacy gzip table.
The reference data own CM/MAF; not every persisted regression row exposes all reference columns as headline output.

**Recommended prevention.** Document where reference-derived coordinates/frequencies are retained and offer a simple explicit export only if a demonstrated workflow needs them; do not copy annotation-supplied CM/MAF as authoritative.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[src/ldsc/_ldscore/outputs.py](<../../../../src/ldsc/_ldscore/outputs.py>); [docs/reference/formats/artifact-metadata.md](<../../../reference/formats/artifact-metadata.md>); [docs/reference/ldscore-snp-universes.md](<../../../reference/ldscore-snp-universes.md>).

Legacy implementation context: [`ldscore/ldscore.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/ldscore.py).

<a id="issue-413"></a>

### #413: How to draw the visialization of the LDSC regression plot like its NG paper?

[Original issue](https://github.com/bulik/ldsc/issues/413) · Source state: **open** · Severity: **low** · Accessible comments reviewed: **3**.

Workflow tags: `h2`.

**Trigger and harm.** Requests the LD-score versus chi-square visualization used in the method paper.
Users lack an accessible diagnostic to inspect signal and fit behavior.

**Legacy disposition.** A commenter suggests plotting base LD score against squared Z; no integrated upstream implementation is shown.

**LDSC3 assessment — Partial coverage.** LDSC3 provides a plot command, but plots must use aligned regression SNPs and saved model context. A raw scatter from unmatched rows is not a fit diagnostic, and a binned mean line differs from the weighted regression.

**Recommended prevention.** Offer a clearly labeled aligned LD-score/chi-square diagnostic with fitted line and filtering context, reusing saved inputs rather than refitting invisibly.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[src/ldsc/cli.py](<../../../../src/ldsc/cli.py>); [docs/user-guide/workflows/plotting.md](<../../../user-guide/workflows/plotting.md>); [docs/reference/ldscore-snp-universes.md](<../../../reference/ldscore-snp-universes.md>).

Legacy implementation context: [`ldscore/regressions.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/regressions.py).

<a id="issue-416"></a>

### #416: how to get the filtered snps

[Original issue](https://github.com/bulik/ldsc/issues/416) · Source state: **open** · Severity: **low** · Accessible comments reviewed: **0**.

Workflow tags: `ldscore`, `h2`.

**Trigger and harm.** Asks how to identify the SNPs retained for a trait after computation.
Without stage-specific retained sets, users cannot explain differences between runs or reproduce the fitted population.

**Legacy disposition.** Open; the short request has no answer.

**LDSC3 assessment — Partial coverage.** LDSC3 records SNP identity in materialized tables and filtering diagnostics, but LD-score rows and trait-specific regression SNPs are different sets. A single generic filtered-SNP file would be ambiguous.

**Recommended prevention.** Document how to inspect each retained set and expose a simple trait-specific regression-row export if saved outputs do not already provide it.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[docs/reference/ldscore-snp-universes.md](<../../../reference/ldscore-snp-universes.md>); [docs/reference/logs-and-outputs.md](<../../../reference/logs-and-outputs.md>); [src/ldsc/_regression/runner.py:495](<../../../../src/ldsc/_regression/runner.py#L495>).

Legacy implementation context: [`ldscore/ldscore.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/ldscore.py), [`ldscore/regressions.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/regressions.py).

<a id="issue-491"></a>

### #491: Number of variants included for computing LD score at each SNP

[Original issue](https://github.com/bulik/ldsc/issues/491) · Source state: **open** · Severity: **low** · Accessible comments reviewed: **0**.

Workflow tags: `ldscore`.

**Trigger and harm.** Requests how many reference variants contribute to each SNP LD score.
Users cannot distinguish window density/filtering changes from genuine LD differences using a score alone.

**Legacy disposition.** Open with no implemented output or answer.

**LDSC3 assessment — Partial coverage.** LDSC3 knows exact LD windows and records aggregate SNP funnels, but a per-row contributor count is a distinct diagnostic from annotation sums or regression-SNP counts.
The score value cannot identify that count.

**Recommended prevention.** If demanded, expose the count while traversing already resolved windows, with separate all-reference and annotation-weighted meanings; avoid an extra all-pairs recomputation.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[src/ldsc/_refpanel/windows.py](<../../../../src/ldsc/_refpanel/windows.py>); [src/ldsc/_ldscore/compute.py](<../../../../src/ldsc/_ldscore/compute.py>); [docs/reference/ldscore-snp-universes.md](<../../../reference/ldscore-snp-universes.md>).

Legacy implementation context: [`ldscore/ldscore.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/ldscore.py).

<a id="issue-188"></a>

### #188: R implementation

[Original issue](https://github.com/bulik/ldsc/issues/188) · Source state: **open** · Severity: **low** · Accessible comments reviewed: **0**.

Workflow tags: Cross-cutting / not assigned to a requested group.

**Trigger and harm.** A user requests an R implementation of LDSC. R users face language integration friction, but no failure of the estimator is demonstrated.

**Legacy disposition.** Open with no response or linked implementation in the thread.

**LDSC3 assessment — Intentionally unsupported.** LDSC3 is a Python package with a CLI and reusable Python API; the audited repository does not provide a native R implementation.
Scientific TSV outputs and CLI invocation already allow R workflows without maintaining a second numerical implementation.

**Recommended prevention.** Offer a concise R system2 invocation/read.delim example before considering an independent R port.

**Evidence and verification.** source-inspected
[src/ldsc/__init__.py](<../../../../src/ldsc/__init__.py>); [src/ldsc/cli.py](<../../../../src/ldsc/cli.py>); [docs/user-guide/python-usage.md](<../../../user-guide/python-usage.md>).

<a id="faq-06"></a>

### FAQ-06: Retired automatic binning failure

[FAQ source](https://github.com/bulik/ldsc/wiki/FAQ) · Source-order identifier: **FAQ-06** · Severity: **low**.

Workflow tags: `ldscore`, `continuous`.

**Legacy disposition.** Historical FAQ guidance; assessed independently rather than treated as a code-fix record.

**LDSC3 assessment — Intentionally unsupported.** The historical --cts-bin workflow is absent.
Supplied numeric annotations are supported, while quantile-h2 projects a saved fitted model after regression; it does not refit binary bins.

**Recommended prevention.** Offer the explicit numeric-annotation or separately constructed binary-bin workflow, with its different scientific question.

**Evidence and verification.** Source-inspected; relevant executed checks are identified in verification.md.
[docs/reference/legacy-flag-map.md:143](<../../../reference/legacy-flag-map.md#legacy-ldsc-options-with-no-current-equivalent>); [src/ldsc/_regression/quantile_h2.py](<../../../../src/ldsc/_regression/quantile_h2.py>).

<a id="issue-5"></a>

### #5: space

[Original issue](https://github.com/bulik/ldsc/issues/5) · Source state: **closed** · Severity: **low** · Accessible comments reviewed: **3**.

Workflow tags: `ldscore`.

**Trigger and harm.** Many cell-type LD-score files consume hundreds of gigabytes.
Storage and repeated decompression become practical barriers to large annotation collections.

**Legacy disposition.** Maintainer agreed to gzip by default; current legacy LDSC emits .gz LD scores, although the issue does not identify the introducing commit.

**LDSC3 assessment — Apparently prevented.** LDSC3 reads compressed annotation and legacy input tables and writes Parquet LD-score directories with separate baseline/query storage.
Query batches avoid requiring all query columns in memory.
This changes the storage interface and addresses the motivating use case, but it does not establish production-scale disk requirements or eliminate the need to budget retained R² matrices.

**Recommended prevention.** Document estimated output sizes from measured runs and which reusable files can be archived; prefer existing compression and batching over another storage framework.

**Evidence and verification.** source-inspected
[src/ldsc/_ldscore/outputs.py](<../../../../src/ldsc/_ldscore/outputs.py>); [src/ldsc/_ldscore/source.py](<../../../../src/ldsc/_ldscore/source.py>); [docs/dev/validation/performance.md](<../performance.md>).

Legacy implementation context: [`ldscore/ldscore.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/ldscore.py).

<a id="issue-8"></a>

### #8: --h2 w/ suffix

[Original issue](https://github.com/bulik/ldsc/issues/8) · Source state: **closed** · Severity: **low** · Accessible comments reviewed: **2**.

Workflow tags: `h2`, `rg`.

**Trigger and harm.** A command guesses between prefix-based .chisq, .gz and pickle artifacts.
Ambiguous filename precedence can select the wrong input.

**Legacy disposition.** Maintainer explicitly confirmed full filenames in dev for v1.0.0.

**LDSC3 assessment — Apparently prevented.** LDSC3 requires --sumstats-file or --sumstats-sources pointing to concrete files or documented path tokens.
It does not infer a .chisq/pickle family from an old prefix.
Legacy .sumstats(.gz) interoperability is explicit, and quoted glob expansion is documented for multiple traits.
This fulfills the practical request without carrying forward obsolete suffix precedence.

**Recommended prevention.** Use exact filenames in examples and diagnose empty or multiply matched path selections before fitting.

**Evidence and verification.** source-inspected
[src/ldsc/_regression/runner.py:_add_regression_sections](<../../../../src/ldsc/_regression/runner.py>); [docs/reference/formats/legacy-sumstats.md](<../../../reference/formats/legacy-sumstats.md>); [docs/reference/input-paths.md](<../../../reference/input-paths.md>).

Legacy implementation context: [`ldscore/regressions.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/regressions.py), [`ldscore/sumstats.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/sumstats.py).

<a id="issue-18"></a>

### #18: change .l2.ldscore file format

[Original issue](https://github.com/bulik/ldsc/issues/18) · Source state: **closed** · Severity: **low** · Accessible comments reviewed: **1**.

Workflow tags: `ldscore`.

**Trigger and harm.** Legacy text files redundantly store CM/MAF and excessively precise floats.
Redundant fields inflate resources and encourage assumptions about implicit column positions.

**Legacy disposition.** Maintainer confirmed dropping CM/MAF and printing three decimals; a8de524740211cb5bdbdac6a6fb3eceba97bf433 explicitly closes #18.

**LDSC3 assessment — Apparently prevented.** LDSC3 writes typed Parquet LD scores and metadata rather than reproducing the old text layout.
It preserves needed reference information separately and converts supported legacy file families explicitly.
Numerical precision is governed by its float32 public storage/float64 accumulation contract, so copying the legacy three-decimal rule into native output would be a regression.

**Recommended prevention.** Keep schema/version and precision documentation near conversion instructions; avoid another user-selectable output-format matrix.

**Evidence and verification.** source-inspected
[src/ldsc/_ldscore/outputs.py](<../../../../src/ldsc/_ldscore/outputs.py>); [docs/reference/formats/artifact-metadata.md](<../../../reference/formats/artifact-metadata.md>); [Source](<https://github.com/bulik/ldsc/commit/a8de524740211cb5bdbdac6a6fb3eceba97bf433>).

Legacy implementation context: [`ldscore/ldscore.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/ldscore.py).

