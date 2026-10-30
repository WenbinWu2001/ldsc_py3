# Inputs, workflows, and diagnostics

Last updated on: 2026-10-28

[Audit overview](README.md) · [Workflow counts and membership](workflow-index.md) · [Verification](verification.md)

## Prevention at the user boundary

Many historical failures originate in input meaning or composition rather than an incorrect numerical formula.
The recurring boundaries are allele/sign conventions, delimiters and aliases, sample-size conventions, reference/annotation identity, cM units, and resource-role confusion.
LDSC3 already has many explicit preflight and output-ownership checks; an audit should retain those protections while improving gaps exposed by real reports.

Annotation files supplied together still require identical SNP rows/order.
Alignment between prepared annotations and a reference panel is a separate operation based on SNP identity.
The annotation SNP grid defines candidate SNPs, whereas the reference-SNP funnel and regression-SNP filters describe different stages.
A low regression count cannot be repaired merely by padding missing annotation rows with zero or relaxing allele/build checks.

## Recommendations by kind

- **Code:** add small actionable diagnoses at existing boundaries, such as a comma-separated header mistaken for one whitespace field, and explicit partial/unestimable statuses where demonstrated.
- **Documentation:** show a short complete input/output example for each workflow, including quantities and assumptions, and explain what to inspect after unexpected attrition.
- **Features:** prefer discoverability of current flags, output ledgers and supported formats over aliases for every historical command.
- **Resources:** name the role of each input and distinguish native LD-score generation from one-time conversion of a supported legacy suite.

## Category summary

83 issues and 6 FAQ entries.

| Severity | Source | Concern | LDSC3 disposition | Evidence level |
| --- | --- | --- | --- | --- |
| High | [#29](#issue-29) | Morgans versus centimorgans | Partial coverage | Source inspection |
| High | [#43](#issue-43) | Signed effects, duplicate aliases and misleading filtering counts | Partial coverage | Source inspection |
| High | [#96](#issue-96) | Baseline/query regression-SNP lists differ | Partial coverage | Source inspection |
| High | [#116](#issue-116) | Summary statistics omit an allele column | Partial coverage | Source inspection |
| High | [#125](#issue-125) | Median-beta failure with only 27k variants | Partial coverage | Source inspection |
| High | [#320](#issue-320) | Signed columns, regression weights and negative denominator | Partial coverage | Source inspection |
| High | [#421](#issue-421) | Whole-file input versus chromosome prefix | Partial coverage | Source inspection |
| High | [#433](#issue-433) | LD-window unit mistake | Partial coverage | Source inspection |
| High | [#456](#issue-456) | Allele QC versus low signal | Partial coverage | Source inspection |
| High | [FAQ-24](#faq-24) | Too few regression SNPs | Partial coverage | Source inspection |
| High | [#119](#issue-119) | Discordant alleles empty the retained GWAS | Apparently prevented | Source inspection |
| High | [#133](#issue-133) | A1 is effect allele, not genome reference | Apparently prevented | Source inspection |
| High | [#171](#issue-171) | Missing cM collapses an LD window to the whole chromosome | Apparently prevented | Source inspection |
| High | [#226](#issue-226) | Ambiguous and incompatible allele pairs | Apparently prevented | Source inspection |
| High | [#233](#issue-233) | Baseline/query SNP grids differ | Apparently prevented | Source inspection |
| High | [#254](#issue-254) | Successful run writes only a log | Apparently prevented | Source inspection |
| High | [#274](#issue-274) | Cross-ancestry query/baseline LD mismatch | Apparently prevented | Source inspection |
| High | [#304](#issue-304) | Constant N is a fallback, not override | Apparently prevented | Source inspection |
| High | [#315](#issue-315) | BETA direction and explicit SNP alias precedence | Apparently prevented | Source inspection |
| High | [#321](#issue-321) | Reference SNPs confused with printed regression SNPs | Apparently prevented | Source inspection |
| High | [#326](#issue-326) | Missing per-SNP N column | Apparently prevented | Source inspection |
| High | [#385](#issue-385) | EAS query/baseline rows | Apparently prevented | Source inspection |
| High | [#389](#issue-389) | Baseline is not regression weights | Apparently prevented | Source inspection |
| High | [#480](#issue-480) | Null value mistaken for column index | Apparently prevented | Source inspection |
| High | [FAQ-05](#faq-05) | Annotation and reference alignment | Apparently prevented | Source inspection |
| High | [FAQ-07](#faq-07) | Accidental chromosome-wide windows | Apparently prevented | Source inspection |
| Moderate | [#12](#issue-12) | Undocumented frequency-table format | Partial coverage | Source inspection |
| Moderate | [#25](#issue-25) | Pair failures hidden in long correlation logs | Partial coverage | Source inspection |
| Moderate | [#27](#issue-27) | Broken documentation links and reference/weight ambiguity | Partial coverage | Source inspection |
| Moderate | [#37](#issue-37) | Absolute median guard depends on effect units | Partial coverage | Source inspection |
| Moderate | [#65](#issue-65) | Mislabelled OR and headerless external summaries | Partial coverage | Source inspection |
| Moderate | [#131](#issue-131) | Missing SNP header and values parsed in the wrong column | Partial coverage | Source inspection |
| Moderate | [#134](#issue-134) | Strong associations removed by ordinary QC | Partial coverage | Source inspection |
| Moderate | [#145](#issue-145) | Large munging chunks cause extreme slowdowns | Partial coverage | Source inspection |
| Moderate | [#150](#issue-150) | Synthetic data all fail allele/QC checks | Partial coverage | Source inspection |
| Moderate | [#190](#issue-190) | Gene list yields all-zero annotation | Partial coverage | Source inspection |
| Moderate | [#253](#issue-253) | Empty munging concatenation | Partial coverage | Source + image inspection |
| Moderate | [#281](#issue-281) | Munging chunks versus LD compute blocks | Partial coverage | Source inspection |
| Moderate | [#297](#issue-297) | Single-cell-derived gene-list workflow | Partial coverage | Source inspection |
| Moderate | [#316](#issue-316) | Gzip CRC corruption error | Partial coverage | Source inspection |
| Moderate | [#348](#issue-348) | All-zero annotation | Partial coverage | Source inspection |
| Moderate | [#374](#issue-374) | Delimiter versus header name | Partial coverage | CSV/TSV header probe |
| Moderate | [#404](#issue-404) | Inconsistent text-field counts | Partial coverage | Source inspection |
| Moderate | [#420](#issue-420) | Chunk-dependent membership slowdown | Partial coverage | Source inspection |
| Moderate | [FAQ-03](#faq-03) | Historical time and memory estimates | Partial coverage | Source inspection |
| Moderate | [#57](#issue-57) | LD Hub service UX and data-release confusion | Intentionally unsupported | Source inspection |
| Moderate | [#237](#issue-237) | Thin annotation lacks overlap metadata | Intentionally unsupported | Source inspection |
| Moderate | [#279](#issue-279) | Thin annotation read as full metadata file | Intentionally unsupported | Source inspection |
| Moderate | [#478](#issue-478) | Malformed path and interpreter claim | Indeterminate | Source inspection |
| Moderate | [#13](#issue-13) | Munged summary-statistics contract | Apparently prevented | Source inspection |
| Moderate | [#23](#issue-23) | Tutorial and release mismatch | Apparently prevented | Source inspection |
| Moderate | [#76](#issue-76) | P values plus direction can supply signed Z | Apparently prevented | Source inspection |
| Moderate | [#83](#issue-83) | Malformed cell-type manifest has empty lines | Apparently prevented | Source inspection |
| Moderate | [#110](#issue-110) | Annotation rows differ from output LD-score rows | Apparently prevented | Source inspection |
| Moderate | [#111](#issue-111) | How a query joins the baseline | Apparently prevented | Source inspection |
| Moderate | [#143](#issue-143) | Quantitative traits do not require odds ratios | Apparently prevented | Source inspection |
| Moderate | [#146](#issue-146) | Leading padded rows mistaken for failed munging | Apparently prevented | Source inspection |
| Moderate | [#147](#issue-147) | Relative cell-type manifest paths and missing downloads | Apparently prevented | Source inspection |
| Moderate | [#149](#issue-149) | Missing LD-score prefix mistaken for invalid resources | Apparently prevented | Source inspection |
| Moderate | [#170](#issue-170) | Padded HM3 output interpreted as missing processing | Apparently prevented | Source inspection |
| Moderate | [#183](#issue-183) | Absolute and manifest-relative paths | Apparently prevented | Source inspection |
| Moderate | [#192](#issue-192) | Combining annotation files and model roles | Apparently prevented | Source inspection |
| Moderate | [#238](#issue-238) | Multiple signed-statistic columns | Apparently prevented | Source inspection |
| Moderate | [#240](#issue-240) | Missing CHR in LD-score input | Apparently prevented | Source inspection |
| Moderate | [#250](#issue-250) | No chromosome files resolved | Apparently prevented | Source inspection |
| Moderate | [#260](#issue-260) | LD Hub downstream file missing | Apparently prevented | Source inspection |
| Moderate | [#278](#issue-278) | LD-score command expects missing annotation | Apparently prevented | Source inspection |
| Moderate | [#291](#issue-291) | Missing null in signed-statistic option | Apparently prevented | Source inspection |
| Moderate | [#295](#issue-295) | Computing a custom reference panel | Apparently prevented | Source inspection |
| Moderate | [#309](#issue-309) | Incomplete LD-score prefix | Apparently prevented | Source inspection |
| Moderate | [#333](#issue-333) | Matching SNP list still drops ambiguous alleles | Apparently prevented | Source inspection |
| Moderate | [#335](#issue-335) | N-min documentation contradicts code | Apparently prevented | Source inspection |
| Moderate | [#349](#issue-349) | Missing weight prefix | Apparently prevented | Source inspection |
| Moderate | [#366](#issue-366) | Opaque empty prefix discovery | Apparently prevented | Source inspection |
| Moderate | [#392](#issue-392) | Wrong annotation path | Apparently prevented | Source inspection |
| Moderate | [#401](#issue-401) | Adding a query to baseline | Apparently prevented | Source inspection |
| Moderate | [#408](#issue-408) | Identity dtype mismatch | Apparently prevented | Source inspection |
| Moderate | [#426](#issue-426) | Oversized genotype batch | Apparently prevented | Source inspection |
| Moderate | [#461](#issue-461) | Custom cell-type annotation recipe | Apparently prevented | Source inspection |
| Moderate | [#466](#issue-466) | Chromosome number embedded twice | Apparently prevented | Source inspection |
| Moderate | [#477](#issue-477) | Constant N versus N column | Apparently prevented | Source inspection |
| Moderate | [#482](#issue-482) | Missing path separator | Apparently prevented | Source inspection |
| Low | [FAQ-14](#faq-14) | Regression execution time | Partial coverage | Source inspection |
| Low | [#4](#issue-4) | Module/script naming collision | Apparently prevented | Source inspection |
| Low | [#7](#issue-7) | Truncated trait names in correlation summaries | Apparently prevented | Source inspection |
| Low | [#41](#issue-41) | Tutorial output-name and line-continuation mismatch | Apparently prevented | Source inspection |
| Low | [#135](#issue-135) | P-value alias confused with service input format | Apparently prevented | Source inspection |
| Low | [#387](#issue-387) | Shell flag as command | Apparently prevented | Source inspection |
| Low | [FAQ-02](#faq-02) | Chromosome-level parallel execution | Apparently prevented | Source inspection |

## Individual assessments

<a id="issue-29"></a>

### #29: Option for map column to be in Morgans or cM

[Original issue](https://github.com/bulik/ldsc/issues/29) · Source state: **closed** · Severity: **high** · Accessible comments reviewed: **1**.

Workflow tags: `ldscore`.

**Trigger and harm.** The BIM genetic-distance column is in Morgans while the requested LD window is interpreted in cM. A hundredfold LD-window mismatch silently changes scores and downstream estimates.

**Legacy disposition.** Maintainer retained cM expectations and suggested scaling the window for a Morgan-valued BIM; no automatic-unit fix was promised.

**LDSC3 assessment — Partial coverage.** LDSC3 calls the field CM and uses it directly for PLINK cM windows.
It has kb/SNP alternatives and explicit genetic-map inputs for R² construction, but no reliable automatic Morgan detector was found.
Plausible numeric magnitudes alone cannot prove units across small genomic regions.

**Recommended prevention.** Require documentation to say BIM column 3 must be cM for cM windows, show explicit conversion to a new file, and avoid heuristic unit guessing.

**Evidence and verification.** source-inspected
[src/ldsc/_refpanel/windows.py:window_coordinates](<../../../../src/ldsc/_refpanel/windows.py>); [src/ldsc/_refpanel/adapters.py](<../../../../src/ldsc/_refpanel/adapters.py>); [docs/reference/methods/ld-score-computation.md](<../../../reference/methods/ld-score-computation.md>).

Legacy implementation context: [`ldscore/ldscore.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/ldscore.py).

<a id="issue-43"></a>

### #43: fatal error while running munge_sumstats.py

[Original issue](https://github.com/bulik/ldsc/issues/43) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **23**.
The API reports 24 comments; one was unavailable through both retrieval methods.

Workflow tags: `munge-sumstats`, `h2`, `rg`.

**Trigger and harm.** A user makes negative betas positive to bypass a median check, then encounters duplicate allele headers and confusing HM3 output counts.
Destroying effect signs invalidates genetic correlation; ambiguous mappings and count logs conceal the analysis population.

**Legacy disposition.** Thread identifies explicit/automatic alias collisions and a double-counted MAF log; it also corrects percentages passed as prevalences.
No single verified legacy fix covers all subtopics.
One of 24 expected comments is unavailable.

**LDSC3 assessment — Partial coverage.** LDSC3 resolves duplicate aliases, validates prevalence in (0,1), and outputs retained SNPs rather than padding every HM3 row.
The absolute median guard remains scale-sensitive and dropped-SNP sidecars do not enumerate every early QC removal.
Its HM3 default is a regression-SNP restriction, not evidence that only those SNPs contribute LD.

**Recommended prevention.** Explicitly prohibit taking absolute betas as repair, document signed-Z/units, and distinguish mutually exclusive removal counts from cumulative counts and padded rows.

**Evidence and verification.** source-inspected
[src/ldsc/_sumstats/munge_input.py](<../../../../src/ldsc/_sumstats/munge_input.py>); [src/ldsc/_sumstats/munge_kernel.py:parse_dat](<../../../../src/ldsc/_sumstats/munge_kernel.py>); [src/ldsc/_regression/prevalence.py](<../../../../src/ldsc/_regression/prevalence.py>); [docs/reference/ldscore-snp-universes.md](<../../../reference/ldscore-snp-universes.md>).

Legacy implementation context: [`munge_sumstats.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/munge_sumstats.py), [`ldscore/regressions.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/regressions.py), [`ldscore/sumstats.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/sumstats.py).

<a id="issue-96"></a>

### #96: Baseline Model for Partitioning Analysis

[Original issue](https://github.com/bulik/ldsc/issues/96) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **7**.

Workflow tags: `ldscore`, `ldsc-seg`.

**Trigger and harm.** Two public HapMap3 lists and baseline resources have slightly different output rows.
Users may recompute inconsistent scores or restrict reference SNPs when intending only output rows.

**Legacy disposition.** Comment identifies w_hm3.snplist SNP column as the matching historical list; no universal legacy patch is required.

**LDSC3 assessment — Partial coverage.** LDSC3 computes baseline and query scores on one selected annotation SNP grid and stores aligned rows together, with separate regression-SNP restriction.
Its bundled HM3 map is a newly curated resource and should not be assumed byte-equivalent to either historical list.
Legacy mixed precomputed suites are constrained by conversion support.

**Recommended prevention.** Document grid versus regression-SNP roles and require a coherent baseline/query generation recipe; compare matched identities rather than only row counts.

**Evidence and verification.** source-inspected
[src/ldsc/_ldscore/calculator.py](<../../../../src/ldsc/_ldscore/calculator.py>); [docs/reference/ldscore-snp-universes.md](<../../../reference/ldscore-snp-universes.md>); [docs/user-guide/workflows/migrating-from-legacy-ldsc.md](<../../../user-guide/workflows/migrating-from-legacy-ldsc.md>).

Legacy implementation context: [`ldscore/ldscore.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/ldscore.py), [`ldscore/regressions.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/regressions.py).

<a id="issue-116"></a>

### #116: Incomplete data - A1/A2 column is not available

[Original issue](https://github.com/bulik/ldsc/issues/116) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **1**.

Workflow tags: `munge-sumstats`, `h2`, `rg`.

**Trigger and harm.** A2 is absent, but user wants heritability or genetic correlation.
Inventing the missing allele can flip signs or confuse multiallelic variants.

**Legacy disposition.** Maintainer distinguishes allele-free heritability from correlation requiring both alleles and careful reference lookup.

**LDSC3 assessment — Partial coverage.** LDSC3 offers explicit base SNP identifier modes for appropriate unsigned heritability workflows and requires alleles for correlation.
Default allele-aware modes deliberately reject missing allele metadata.
Reference lookup needs genome build and multiallelic/strand checks, not a guessed complement; no universal imputation helper is warranted.

**Recommended prevention.** Provide separate P-only/base-mode heritability and signed-allele correlation recipes; fail with the scientific reason when correlation lacks alleles.

**Evidence and verification.** source-inspected
[src/ldsc/_regression/runner.py:_require_rg_alleles](<../../../../src/ldsc/_regression/runner.py>); [docs/reference/snp-identity-and-genome-build.md](<../../../reference/snp-identity-and-genome-build.md>); [docs/reference/munge-sumstats.md](<../../../reference/munge-sumstats.md>).

Legacy implementation context: [`munge_sumstats.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/munge_sumstats.py), [`ldscore/regressions.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/regressions.py), [`ldscore/sumstats.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/sumstats.py).

<a id="issue-125"></a>

### #125: Error converting summary statistic

[Original issue](https://github.com/bulik/ldsc/issues/125) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **3**.

Workflow tags: `munge-sumstats`.

**Trigger and harm.** Raw beta median .41 fails the guard; N and SNP count both equal 27,130 in the pasted command. A selected/small marker set or mistaken sample size can undermine inference beyond the visible median error.

**Legacy disposition.** No maintainer resolution; a later user suggests signed Z from beta/SE, which is valid only with correct source semantics.

**LDSC3 assessment — Partial coverage.** LDSC3 retains the absolute median guard, so this trigger can remain.
Correct signed Z may avoid unit dependence, but it does not fix wrong N, nonrepresentative SNP selection or insufficient signal.
The identical N/row counts warrant checking, not a definitive claim that the user substituted SNP count for people.

**Recommended prevention.** Explain N is participants and request source-verified units/direction; do not recenter effects or relax QC solely to obtain output.

**Evidence and verification.** source-inspected
[src/ldsc/_sumstats/munge_kernel.py:check_median](<../../../../src/ldsc/_sumstats/munge_kernel.py>); [src/ldsc/_sumstats/munge_kernel.py:process_n](<../../../../src/ldsc/_sumstats/munge_kernel.py>); [docs/reference/munge-sumstats.md](<../../../reference/munge-sumstats.md>).

Legacy implementation context: [`munge_sumstats.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/munge_sumstats.py).

<a id="issue-320"></a>

### #320: Which --signed-sumstats?  And NA genetic correlation

[Original issue](https://github.com/bulik/ldsc/issues/320) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **2**.

Workflow tags: `munge-sumstats`, `rg`.

**Trigger and harm.** A user asks which beta/Z/OR to select, why reference paths repeat and why genetic correlation isNA.
Three distinct input/method concepts are conflated, encouraging arbitrary settings.

**Legacy disposition.** Reporter later supplies a working command, but its changed paths/N do not identify a single general fix.

**LDSC3 assessment — Partial coverage.** LDSC3 distinguishes P-derived magnitude from direction, centralizes regression weights in LD-score directories and reports pair-specific SNP heritability/failed status.
Equivalent correctly coded beta/Z/OR directions should agree; negative denominator information cannot be repaired by choosing a different column solely for plausible output.

**Recommended prevention.** Provide one task-oriented diagnostic sequence: inspect column mapping, validate reference/weights, then assess denominator signal and uncertainty.

**Evidence and verification.** source-inspected
[docs/user-guide/workflows/munge-sumstats.md](<../../../user-guide/workflows/munge-sumstats.md>); [docs/reference/ldscore-snp-universes.md](<../../../reference/ldscore-snp-universes.md>); [docs/user-guide/workflows/genetic-correlation.md](<../../../user-guide/workflows/genetic-correlation.md>).

Legacy implementation context: [`munge_sumstats.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/munge_sumstats.py), [`ldscore/regressions.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/regressions.py), [`ldscore/sumstats.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/sumstats.py).

<a id="issue-421"></a>

### #421: lsdc :IndexError: list index out of range

[Original issue](https://github.com/bulik/ldsc/issues/421) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **7**.

Workflow tags: `rg`.

**Trigger and harm.** An AFR LD-score table is passed to legacy chromosome-split flags; the two traits may also have different ancestries. A path error blocks analysis, while a suggested EUR replacement can make a run execute with the wrong scientific reference.

**Legacy disposition.** A commenter identifies --ref-ld/--w-ld as the relevant legacy flags; no confirmed final fix.

**LDSC3 assessment — Partial coverage.** LDSC3 removes these parallel prefix modes in favor of an explicit LD-score directory.
Its supported legacy converter requires complete standard families, so an arbitrary AFR table is not a drop-in input, and cross-ancestry correlation needs separate methodological justification.

**Recommended prevention.** Diagnose the input format first and preserve ancestry matching; never treat a successful EUR-reference workaround as validation of an AFR analysis.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[docs/reference/formats/legacy-ldscores.md](<../../../reference/formats/legacy-ldscores.md>); [src/ldsc/_ldscore/legacy_converter.py](<../../../../src/ldsc/_ldscore/legacy_converter.py>); [docs/reference/input-paths.md](<../../../reference/input-paths.md>).

Legacy implementation context: [`ldscore/regressions.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/regressions.py), [`ldscore/sumstats.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/sumstats.py).

<a id="issue-433"></a>

### #433: African data analysis

[Original issue](https://github.com/bulik/ldsc/issues/433) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **1**.

Workflow tags: `ldscore`, `h2`.

**Trigger and harm.** AFR LD scores are computed with a 1-kb window while the user expects a conventional broader LD window. A syntactically valid short window can omit relevant LD and change regression behavior.

**Legacy disposition.** Comment recommends a broader kb window and chromosome parallelism, without matching the target reference scientifically.

**LDSC3 assessment — Partial coverage.** LDSC3 requires an explicit window and preserves unit/value metadata, but cannot infer that 1kb was unintended. A universal replacement of 1cM by 1000kb is also not scientifically exact.

**Recommended prevention.** Show units in examples and preflight summaries, distinguish cM from physical distance, and match the intended panel/model before optimizing runtime.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[docs/reference/ldscore-snp-universes.md](<../../../reference/ldscore-snp-universes.md>); [src/ldsc/_refpanel/windows.py](<../../../../src/ldsc/_refpanel/windows.py>); [docs/reference/parallel-execution.md](<../../../reference/parallel-execution.md>).

Legacy implementation context: [`ldscore/ldscore.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/ldscore.py), [`ldscore/regressions.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/regressions.py).

<a id="issue-456"></a>

### #456: valid alleles QC

[Original issue](https://github.com/bulik/ldsc/issues/456) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **1**.

Workflow tags: `munge-sumstats`, `rg`.

**Trigger and harm.** An EAS rg run drops about 81k rows at allele QC and then has negative h2 under fixed intercepts.
Explaining only input row-count differences misses both the allele-loss stage and an unestimable genetic-correlation denominator.

**Legacy disposition.** A commenter discusses the intersection count, not the later allele-QC loss; no cause or fix is verified.

**LDSC3 assessment — Partial coverage.** LDSC3 aligns alleles explicitly and records filtering, but cannot infer the cause of these lost rows without inspecting pairs.
It should not rescue a negative denominator by clipping or an unjustified intercept constraint.

**Recommended prevention.** Separate overlap loss, invalid allele, discordance and ambiguous-strand counts; explain low-signal estimability before changing intercepts or allowing incompatible alleles.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[src/ldsc/_regression/estimators.py:781](<../../../../src/ldsc/_regression/estimators.py#L781>); [src/ldsc/_regression/estimators.py:765](<../../../../src/ldsc/_regression/estimators.py#L765>); [docs/user-guide/workflows/genetic-correlation.md](<../../../user-guide/workflows/genetic-correlation.md>).

Legacy implementation context: [`munge_sumstats.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/munge_sumstats.py), [`ldscore/regressions.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/regressions.py), [`ldscore/sumstats.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/sumstats.py).

<a id="faq-24"></a>

### FAQ-24: Too few regression SNPs

[FAQ source](https://github.com/bulik/ldsc/wiki/FAQ) · Source-order identifier: **FAQ-24** · Severity: **high**.

Workflow tags: `h2`, `rg`, `baseline-partitioned-h2`, `ldsc-seg`.

**Legacy disposition.** Historical FAQ guidance; assessed independently rather than treated as a code-fix record.

**LDSC3 assessment — Partial coverage.** Input matching and filtering must be reviewed stage by stage. A retained-SNP count alone is not a universal validity criterion; the effective genome-wide information and intended population matter.

**Recommended prevention.** Make unexpected attrition easy to find and repair build/identity/resource mismatches instead of bypassing warnings by lowering thresholds.

**Evidence and verification.** Source-inspected; relevant executed checks are identified in verification.md.
[docs/reference/ldscore-snp-universes.md](<../../../reference/ldscore-snp-universes.md>); [src/ldsc/_regression/runner.py:RegressionRunner.prepare_inputs](<../../../../src/ldsc/_regression/runner.py>); [src/ldsc/_regression/runner.py:RegressionRunner.build_rg_dataset](<../../../../src/ldsc/_regression/runner.py>).

<a id="issue-119"></a>

### #119: Error converting summary statistics munge_sumstat

[Original issue](https://github.com/bulik/ldsc/issues/119) · Source state: **closed** · Severity: **high** · Accessible comments reviewed: **3**.

Workflow tags: `munge-sumstats`.

**Trigger and harm.** Reported rsIDs have allele pairs incompatible with the HM3 list, including multiallelic ambiguity.
No retained data can be misdiagnosed as a software/environment failure.

**Legacy disposition.** Maintainer identifies allele discordance from example rows and user confirms that cause.

**LDSC3 assessment — Apparently prevented.** LDSC3 validates allele-aware identities and records projection/restriction outcomes rather than trusting a matching rsID alone.
The correct remedy is a source-verified representation, not changing allele labels until rows match.
An empty restriction should identify the active mode, candidate overlap and allele rejection reasons.

**Recommended prevention.** Keep explicit incompatible/multiallelic diagnostics and show how to inspect representative rejected records without auto-rewriting scientific inputs.

**Evidence and verification.** source-inspected
[src/ldsc/_core/snp_identity.py](<../../../../src/ldsc/_core/snp_identity.py>); [src/ldsc/_sumstats/munge_kernel.py:_raise_sumstats_restriction_empty](<../../../../src/ldsc/_sumstats/munge_kernel.py>); [src/ldsc/_regression/runner.py:_project_legacy_sumstats_to_panel](<../../../../src/ldsc/_regression/runner.py>).

Legacy implementation context: [`munge_sumstats.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/munge_sumstats.py).

<a id="issue-133"></a>

### #133: Question about A1/A2 alleles for munge_sumstat.py

[Original issue](https://github.com/bulik/ldsc/issues/133) · Source state: **closed** · Severity: **high** · Accessible comments reviewed: **3**.

Workflow tags: `munge-sumstats`, `rg`.

**Trigger and harm.** User maps A1 to genome ref although Z is relative to alt. A sign reversal can reverse genetic-correlation direction without a parsing error.

**Legacy disposition.** Maintainer explains overloaded reference terminology and confirms the frequency-option label typo had no functional effect; no source commit is identified.

**LDSC3 assessment — Apparently prevented.** LDSC3 input descriptions explicitly define A1 as the allele the signed statistic is relative to and A2 as counterpart.
Downstream harmonization cannot repair a wrongly declared signed-statistic allele because it trusts the input semantics.
The approved vocabulary also distinguishes this from genomic reference and reference-panel minor allele.

**Recommended prevention.** Show one swapped-allele/signed-Z example and clearly separate effect, minor and genome-reference alleles; never flip beta independently of its allele.

**Evidence and verification.** source-inspected
[src/ldsc/_sumstats/munge_input.py:describe_cname](<../../../../src/ldsc/_sumstats/munge_input.py>); [docs/dev/DOCS_VOCABULARY.md](<../../DOCS_VOCABULARY.md>); [src/ldsc/_regression/runner.py:build_rg_dataset](<../../../../src/ldsc/_regression/runner.py>).

Legacy implementation context: [`munge_sumstats.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/munge_sumstats.py), [`ldscore/regressions.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/regressions.py), [`ldscore/sumstats.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/sumstats.py).

<a id="issue-171"></a>

### #171: Error when calculating LD scores

[Original issue](https://github.com/bulik/ldsc/issues/171) · Source state: **closed** · Severity: **high** · Accessible comments reviewed: **3**.

Workflow tags: `ldscore`.

**Trigger and harm.** BIM CM values are all zero while --ld-wind-cm 1 is requested.
Bypassing a generic whole-chromosome warning can launch costly scientifically meaningless scoring.

**Legacy disposition.** User identifies missing cM; later comment suggests physical kb windows.
No algorithm defect is needed to explain the report.

**LDSC3 assessment — Apparently prevented.** LDSC3 assert_cm_usable explicitly rejects fewer than two distinct CM values and suggests a real genetic map or kb/SNP windows. --yes-really does not bypass missing/uninformative CM. A physical window is a valid deliberate alternative, not a numerical conversion preserving a 1cM scientific target.

**Recommended prevention.** Keep the actionable early failure and explain window units; validate real map coverage rather than fill CM with zero or blindly authorize chromosome-wide work.

**Evidence and verification.** source-inspected
[src/ldsc/_refpanel/windows.py:assert_cm_usable](<../../../../src/ldsc/_refpanel/windows.py>); [src/ldsc/_refpanel/windows.py:check_whole_chromosome_window](<../../../../src/ldsc/_refpanel/windows.py>); [docs/reference/methods/ld-score-computation.md](<../../../reference/methods/ld-score-computation.md>).

Legacy implementation context: [`ldscore/ldscore.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/ldscore.py).

<a id="issue-226"></a>

### #226: What's Incompatible alleles mean?

[Original issue](https://github.com/bulik/ldsc/issues/226) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **3**.

Workflow tags: `rg`.

**Trigger and harm.** Genetic correlation reports an incompatible ATAT allele combination; later comments discuss multiallelic rsIDs.
Wrong allele alignment can reverse genetic effects or match distinct variants.

**Legacy disposition.** Comments suggest removing ambiguous or incompatible alleles; no specific merged fix is established.

**LDSC3 assessment — Apparently prevented.** LDSC3 validates allele-compatible rsID projection and uses explicit Boolean match/flip masks; ambiguous strand pairs and multiallelic matches are diagnosed rather than blindly indexed in FLIP_ALLELES.
Coordinate-family allele-aware identity also distinguishes variants sharing an rsID, but does not replace orientation checks.

**Recommended prevention.** Explain ambiguity versus a true mismatch, preserve row-level drop reasons, and discourage manually flipping all alleles without checking signed effects.

**Evidence and verification.** source-inspected
[src/ldsc/_regression/runner.py:_orient_sumstats_z_to_reference_alleles](<../../../../src/ldsc/_regression/runner.py>); [src/ldsc/_regression/estimators.py:_filter_alleles](<../../../../src/ldsc/_regression/estimators.py>); [docs/reference/formats/legacy-sumstats.md](<../../../reference/formats/legacy-sumstats.md>).

Legacy implementation context: [`ldscore/regressions.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/regressions.py), [`ldscore/sumstats.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/sumstats.py).

<a id="issue-233"></a>

### #233: ValueError: LD Scores for concatenation must have identical SNP columns.

[Original issue](https://github.com/bulik/ldsc/issues/233) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **9**.

Workflow tags: `ldscore`, `ldsc-seg`.

**Trigger and harm.** Custom peak LD-score files have98691 SNPs while distributed baseline files have98643.
Combining unmatched rows can misalign annotations; rebuilding baseline from incomplete BED sources can change the model.

**Legacy disposition.** Commenters suggest the original print-SNP list or regeneration and later Zenodo8292725; no code fix is shown.

**LDSC3 assessment — Apparently prevented.** LDSC3 calculates baseline and query LD scores on one validated annotation SNP grid and separates reference SNPs from regression SNPs.
Native directories avoid manually concatenating independently computed chromosome files.
Missing continuous baseline sources cannot be recreated faithfully from binary BED files alone.

**Recommended prevention.** Recompute queries and baseline from a compatible panel/grid, record regression-SNP selection, and avoid blind row intersections or incomplete baseline reconstruction.

**Evidence and verification.** source-inspected
[docs/reference/ldscore-snp-universes.md](<../../../reference/ldscore-snp-universes.md>); [src/ldsc/_ldscore/preflight.py](<../../../../src/ldsc/_ldscore/preflight.py>); [src/ldsc/_ldscore/direct_annotation.py](<../../../../src/ldsc/_ldscore/direct_annotation.py>).

Legacy implementation context: [`ldscore/ldscore.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/ldscore.py), [`ldscore/regressions.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/regressions.py).

<a id="issue-254"></a>

### #254: No results file from partitioned heritability

[Original issue](https://github.com/bulik/ldsc/issues/254) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **2**.

Workflow tags: `ldsc-seg`.

**Trigger and harm.** A single custom annotation without baseline or --overlap-annot yields no .results file. A conditional enrichment question is analyzed with an insufficient model and users scrape incomplete logs.

**Legacy disposition.** Comment identifies legacy results-file emission depends on --overlap-annot; no merged change is supplied.

**LDSC3 assessment — Apparently prevented.** LDSC3 separates h2 from partitioned-h2, requires a baseline for queries and always has fixed-name scientific result outputs for the selected workflow.
Stored overlap replaces the hidden output flag. A native query-plus-baseline model also prevents the single-component proportion1 interpretation.

**Recommended prevention.** Guide users by scientific question and required result table; do not recommend parsing log text to reconstruct scientific outputs.

**Evidence and verification.** source-inspected
[docs/reference/results/partitioned-h2.md](<../../../reference/results/partitioned-h2.md>); [src/ldsc/_regression/runner.py:summarize_partitioned_h2](<../../../../src/ldsc/_regression/runner.py>).

Legacy implementation context: [`ldscore/regressions.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/regressions.py).

<a id="issue-274"></a>

### #274: ValueError: Missing some LD scores from cts files. Are you sure all SNPs in ref-ld-chr are also in ref-ld-chr-cts

[Original issue](https://github.com/bulik/ldsc/issues/274) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **0**.

Workflow tags: `ldscore`, `ldsc-seg`.

**Trigger and harm.** African/Asian custom scores are combined with distributed cell-type scores missing some reference SNPs.
Changing a join to inner can misalign counts/matrices while retaining population-mismatched LD.

**Legacy disposition.** Open; hand-editing sumstats.py produces another dimension error and is not a verified fix.

**LDSC3 assessment — Apparently prevented.** LDSC3 computes query and baseline scores together from one ancestry-appropriate reference and validated grid.
Biological gene-list definitions may be transferred, but their LD scores cannot generally be borrowed across ancestries.
Alignment of rows alone does not establish compatible reference counts or LD.

**Recommended prevention.** Regenerate all LD-dependent components consistently; document why a manual inner join is not a scientific compatibility repair.

**Evidence and verification.** source-inspected
[src/ldsc/_ldscore/preflight.py](<../../../../src/ldsc/_ldscore/preflight.py>); [docs/reference/ldscore-snp-universes.md](<../../../reference/ldscore-snp-universes.md>).

Legacy implementation context: [`ldscore/ldscore.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/ldscore.py), [`ldscore/regressions.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/regressions.py).

<a id="issue-304"></a>

### #304: --N flag ignored by munge_sumstats.py

[Original issue](https://github.com/bulik/ldsc/issues/304) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **0**.

Workflow tags: `munge-sumstats`.

**Trigger and harm.** Legacy help claims --N overrides an existing N column, but code uses the column.
Users believe a sample-size transformation occurred when it did not.

**Legacy disposition.** Reporter recommends correcting documentation; no fixing commit linked.

**LDSC3 assessment — Apparently prevented.** LDSC3 explicitly documents input sample-size columns taking precedence and describes --N as a constant fallback.
The code still follows that policy, so silent reinterpretation is avoided by an accurate contract rather than changing an established normalization rule.
Explicit direct-N and case/control choices are also separated.

**Recommended prevention.** Show the selected N source in preview/provenance; explain how to deliberately exclude a column when a trusted constant is intended.

**Evidence and verification.** source-inspected
[src/ldsc/_sumstats/munge_input.py:_resolve_sample_size_column_strategy](<../../../../src/ldsc/_sumstats/munge_input.py>); [src/ldsc/_sumstats/munge_kernel.py:process_n](<../../../../src/ldsc/_sumstats/munge_kernel.py>); [docs/user-guide/workflows/munge-sumstats.md#sample-size-columns](<../../../user-guide/workflows/munge-sumstats.md#sample-size-columns>).

Legacy implementation context: [`munge_sumstats.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/munge_sumstats.py).

<a id="issue-315"></a>

### #315: IndexError: list index out of range | ValueError: Found 2 different SNP column

[Original issue](https://github.com/bulik/ldsc/issues/315) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **0**.

Workflow tags: `munge-sumstats`.

**Trigger and harm.** User blames BETA,0 for an unrelated IndexError and reports duplicate SNP aliases despite --snp rsid.
Wrong diagnosis obscures path/input errors, while arbitrary allele swaps can reverse effects.

**Legacy disposition.** Open without a resolution; the reported IndexError lacks a traceback to localize it.

**LDSC3 assessment — Apparently prevented.** LDSC3 accepts BETA,0 as direction with P-derived magnitude, and explicit column hints suppress competing default aliases through used_targets. A1/A2 counts remaining unchanged after swapping columns is not evidence that signed effect orientation remains valid.
IndexError cannot be attributed to the signed-statistic choice from this report.

**Recommended prevention.** Show accepted sign/null examples, preview column mapping and explain that swapping A1/A2 requires the corresponding sign change.

**Evidence and verification.** source-inspected
[src/ldsc/_sumstats/munge_input.py:prepare_munge_input](<../../../../src/ldsc/_sumstats/munge_input.py>); [src/ldsc/_sumstats/munge_kernel.py:munge_sumstats](<../../../../src/ldsc/_sumstats/munge_kernel.py>).

Legacy implementation context: [`munge_sumstats.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/munge_sumstats.py).

<a id="issue-321"></a>

### #321: How to replicate pre-computed LD score

[Original issue](https://github.com/bulik/ldsc/issues/321) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **4**.

Workflow tags: `ldscore`.

**Trigger and harm.** User applies --extract HM3 while trying to reproduce EUR/EAS precomputed scores.
Filtering LD-contributing SNPs differs from selecting output/regression rows and changes scores/counts.

**Legacy disposition.** Comments suggest different filtering histories; a recalled private email is not reproducible evidence of equivalence.

**LDSC3 assessment — Apparently prevented.** LDSC3 explicitly separates --ref-panel-snps-file from --regr-snps-file and records a reference-SNP funnel.
Reproducing a historical score still needs exact ancestry/sample QC, LD window, MAF and SNP-selection provenance; similar filenames or small row differences do not establish numerical identity.

**Recommended prevention.** Show contributor-versus-regression selection with a small example and list the information required to reproduce a distributed reference.

**Evidence and verification.** source-inspected
[docs/reference/ldscore-snp-universes.md](<../../../reference/ldscore-snp-universes.md>); [src/ldsc/_ldscore/diagnostics.py](<../../../../src/ldsc/_ldscore/diagnostics.py>); [docs/reference/cli.md](<../../../reference/cli.md>).

Legacy implementation context: [`ldscore/ldscore.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/ldscore.py).

<a id="issue-326"></a>

### #326: N column missing from summary statistics

[Original issue](https://github.com/bulik/ldsc/issues/326) · Source state: **closed** · Severity: **high** · Accessible comments reviewed: **1**.

Workflow tags: `munge-sumstats`, `rg`.

**Trigger and harm.** Some raw GWAS files lack sample size for genetic-correlation preparation.
An invented N may permit execution while miscalibrating SNP heritability and weights.

**Legacy disposition.** Reporter resolves using a constant --N in munging; that is valid only when a trustworthy constant applies.

**LDSC3 assessment — Apparently prevented.** LDSC3 supports a documented constant-N fallback when input N/case-control columns are absent, and requires sample size by the time regression runs.
It does not estimate unknown N from the statistics.
Meta-analysis participation or mixed sample sizes can make a single abstract-level N inappropriate.

**Recommended prevention.** Show the fallback syntax and require documented study-specific N; prefer per-variant N where it varies.

**Evidence and verification.** source-inspected
[src/ldsc/_sumstats/munge_input.py:_resolve_sample_size_column_strategy](<../../../../src/ldsc/_sumstats/munge_input.py>); [src/ldsc/_sumstats/munge_kernel.py:process_n](<../../../../src/ldsc/_sumstats/munge_kernel.py>); [docs/user-guide/workflows/munge-sumstats.md#sample-size-columns](<../../../user-guide/workflows/munge-sumstats.md#sample-size-columns>).

Legacy implementation context: [`munge_sumstats.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/munge_sumstats.py), [`ldscore/regressions.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/regressions.py), [`ldscore/sumstats.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/sumstats.py).

<a id="issue-385"></a>

### #385: Error: Missing some LD scores from cts files. Are you sure all SNPs in ref-ld-chr are also in ref-ld-chr-cts

[Original issue](https://github.com/bulik/ldsc/issues/385) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **4**.

Workflow tags: `ldscore`, `ldsc-seg`.

**Trigger and harm.** Custom EAS query LD scores omit SNPs present in the supplied EAS baseline, possibly due to --print-snps restriction.
Manual intersection may hide inconsistent reference/annotation support or unnecessarily reduce the regression set.

**Legacy disposition.** A user reports manually intersecting files; another suggests an HM3-list mismatch.
Neither establishes a canonical fix.

**LDSC3 assessment — Apparently prevented.** LDSC3 computes baseline and query scores together on declared reference and regression SNP sets, preventing independently assembled cts tables from silently drifting.
Missing annotation support must still not be filled with zero.

**Recommended prevention.** Generate baseline and query from one matched grid and panel; diagnose lost SNPs by stage rather than recommending blind intersections.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[docs/reference/ldscore-snp-universes.md](<../../../reference/ldscore-snp-universes.md>); [src/ldsc/_ldscore/direct_annotation.py](<../../../../src/ldsc/_ldscore/direct_annotation.py>); [src/ldsc/_annotation/builder.py](<../../../../src/ldsc/_annotation/builder.py>).

Legacy implementation context: [`ldscore/ldscore.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/ldscore.py), [`ldscore/regressions.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/regressions.py).

<a id="issue-389"></a>

### #389: Error calculating heritability

[Original issue](https://github.com/bulik/ldsc/issues/389) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **1**.

Workflow tags: `h2`, `baseline-partitioned-h2`.

**Trigger and harm.** A multicolumn baselineLD suite is supplied as --w-ld-chr.
The wrong resource role causes a startup error and invites incorrect single-column extraction.

**Legacy disposition.** A commenter directs the user to the separate HM3-no-MHC regression-weight archive; no code fix is needed.

**LDSC3 assessment — Apparently prevented.** LDSC3 stores weight and predictor roles explicitly in an LD-score directory and the converter takes separate reference/weight suites.
Filename similarity still makes the conceptual distinction worth documenting.

**Recommended prevention.** Show one table contrasting baseline annotations, predictor LD scores and regression weights; do not suggest selecting an arbitrary baseline column as a weight.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[src/ldsc/_ldscore/legacy_converter.py](<../../../../src/ldsc/_ldscore/legacy_converter.py>); [docs/reference/ldscore-snp-universes.md](<../../../reference/ldscore-snp-universes.md>).

Legacy implementation context: [`ldscore/regressions.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/regressions.py).

<a id="issue-480"></a>

### #480: munge_sumstats.py - median value of SIGNED_SUMSTATS is 0.02 (should be close to 6.0). This column may be mislabeled

[Original issue](https://github.com/bulik/ldsc/issues/480) · Source state: **closed** · Severity: **high** · Accessible comments reviewed: **0**.

Workflow tags: `munge-sumstats`.

**Trigger and harm.** Passes --signed-sumstats Effect_Beta,6 because Effect_Beta is the sixth column. A median guard prevents incorrect sign thresholding; bypassing it could invert many valid effects.

**Legacy disposition.** Closed without an accessible resolution or code fix.

**LDSC3 assessment — Apparently prevented.** LDSC3 explicitly defines the comma argument as the null value, with beta/Z null zero and OR null one.
The same median validation is retained and should reject six for an ordinary beta.

**Recommended prevention.** Put the null-value meaning in examples and error remediation, explicitly saying it is not a column number; retain the guard rather than weakening it.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[src/ldsc/_sumstats/munge_input.py:135](<../../../../src/ldsc/_sumstats/munge_input.py#L135>); [src/ldsc/_sumstats/munge_kernel.py:657](<../../../../src/ldsc/_sumstats/munge_kernel.py#L657>); [docs/reference/cli.md:212](<../../../reference/cli.md#ldsc-munge-sumstats>).

Legacy implementation context: [`munge_sumstats.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/munge_sumstats.py).

<a id="faq-05"></a>

### FAQ-05: Annotation and reference alignment

[FAQ source](https://github.com/bulik/ldsc/wiki/FAQ) · Source-order identifier: **FAQ-05** · Severity: **high**.

Workflow tags: `ldscore`.

**Legacy disposition.** Historical FAQ guidance; assessed independently rather than treated as a code-fix record.

**LDSC3 assessment — Apparently prevented.** LDSC3 aligns prepared annotation rows to reference-panel SNP identity, so legacy annotation-versus-BIM physical row order is not the same contract.
However, annotation files read together still require the same SNP rows in the same order; reference alignment does not remove that inter-file requirement.
Genome build, alleles and excluded reference SNPs remain important.

**Recommended prevention.** Document identity-based alignment and inspect excluded reference SNPs; do not repair mismatches by blindly sorting or padding zero rows.

**Evidence and verification.** Source-inspected; relevant executed checks are identified in verification.md.
[src/ldsc/_annotation/parsing.py](<../../../../src/ldsc/_annotation/parsing.py>); [src/ldsc/_refpanel/adapters.py](<../../../../src/ldsc/_refpanel/adapters.py>); [docs/reference/ldscore-snp-universes.md](<../../../reference/ldscore-snp-universes.md>); [src/ldsc/_annotation/sources.py:130](<../../../../src/ldsc/_annotation/sources.py#L130>).

<a id="faq-07"></a>

### FAQ-07: Accidental chromosome-wide windows

[FAQ source](https://github.com/bulik/ldsc/wiki/FAQ) · Source-order identifier: **FAQ-07** · Severity: **high**.

Workflow tags: `ldscore`.

**Legacy disposition.** Historical FAQ guidance; assessed independently rather than treated as a code-fix record.

**LDSC3 assessment — Apparently prevented.** Whole-chromosome LD windows have a dedicated guard. A cM input with no useful genetic distances also has an unconditional check that is not bypassed by --yes-really.

**Recommended prevention.** Keep both checks and explain the appropriate SNP, kb or cM window rather than recommending the override first.

**Evidence and verification.** Source-inspected; relevant executed checks are identified in verification.md.
[src/ldsc/_refpanel/windows.py:70](<../../../../src/ldsc/_refpanel/windows.py#L70>); [src/ldsc/_refpanel/windows.py:112](<../../../../src/ldsc/_refpanel/windows.py#L112>); [tests/ldscore/test_plink_exact_windows.py](<../../../../tests/ldscore/test_plink_exact_windows.py>).

<a id="issue-12"></a>

### #12: document .frq file format

[Original issue](https://github.com/bulik/ldsc/issues/12) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **1**.

Workflow tags: `baseline-partitioned-h2`.

**Trigger and harm.** Users need a frequency file for counts but cannot inspect a hosted example.
Ambiguous schema and inaccessible examples block reproducible resource preparation.

**Legacy disposition.** Still open; a 2023 comment reports GCP access barriers.
No legacy documentation fix is established.

**LDSC3 assessment — Partial coverage.** Native LDSC3 derives reference MAF from its prepared reference panel and uses it for common-SNP counts.
Legacy conversion has a separate input contract rather than treating any frequency table as interchangeable.
Existing format and reference-data pages help, but provenance, sample membership, genome build and frequency versus MAF need explicit worked examples.

**Recommended prevention.** Explain when an external frequency file is needed, show a tiny schema example, and link a maintained source with access instructions.

**Evidence and verification.** source-inspected
[docs/reference/formats/legacy-ldscores.md](<../../../reference/formats/legacy-ldscores.md>); [docs/reference/ldscore-snp-universes.md](<../../../reference/ldscore-snp-universes.md>); [docs/user-guide/reference-data.md](<../../../user-guide/reference-data.md>).

Legacy implementation context: [`ldscore/regressions.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/regressions.py).

<a id="issue-25"></a>

### #25: better error messages for --rg

[Original issue](https://github.com/bulik/ldsc/issues/25) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **0**.

Workflow tags: `rg`.

**Trigger and harm.** A multi-pair run continues after one pair fails and prints missing values in its summary.
Users may mistake partial completion for a fully successful requested analysis.

**Legacy disposition.** Open request for a completed-with-errors footer and repeated error summary; no legacy implementation is verified.

**LDSC3 assessment — Partial coverage.** LDSC3 records pair status and compact errors in saved correlation tables and logs each exception. estimate_rg_pairs intentionally continues after pair-level failures.
The CLI completion log reports the number of trait pairs, but this does not ensure a novice notices failed pairs; successful process completion is not proof every pair was estimated.

**Recommended prevention.** Add a concise final success/failure count and point to failed rows; preserve explicit partial-result statuses instead of imputing estimates.

**Evidence and verification.** source-inspected
[src/ldsc/_regression/runner.py:estimate_rg_pairs](<../../../../src/ldsc/_regression/runner.py>); [src/ldsc/_regression/runner.py:run_rg_from_args](<../../../../src/ldsc/_regression/runner.py>); [src/ldsc/_regression/outputs.py:RG_CONCISE_COLUMNS](<../../../../src/ldsc/_regression/outputs.py>).

Legacy implementation context: [`ldscore/regressions.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/regressions.py), [`ldscore/sumstats.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/sumstats.py).

<a id="issue-27"></a>

### #27: Missing tutorial/FAQ links

[Original issue](https://github.com/bulik/ldsc/issues/27) · Source state: **closed** · Severity: **moderate** · Accessible comments reviewed: **2**.

Workflow tags: `h2`, `rg`.

**Trigger and harm.** README links target removed local paths and tutorial descriptions swap reference and weight directories.
Readers can choose a different LD-score model while believing they reproduced a tutorial.

**Legacy disposition.** Maintainer explains the intended use of HM3-restricted scores for both historical inputs; the issue does not identify a documentation fixing commit.

**LDSC3 assessment — Partial coverage.** LDSC3 has local linked user guides and one LD-score directory carrying predictors, weights and counts.
This removes legacy flag ambiguity but does not make their SNP sets identical: reference SNPs and regression SNPs still have different roles.
Conversion preserves distinct reference and weight families.

**Recommended prevention.** Check documentation links and explain both SNP sets with one concrete example; do not repeat the historical heuristic as a universal resource recommendation.

**Evidence and verification.** source-inspected
[docs/reference/ldscore-snp-universes.md](<../../../reference/ldscore-snp-universes.md>); [docs/reference/formats/legacy-ldscores.md](<../../../reference/formats/legacy-ldscores.md>); [docs/user-guide/workflows/migrating-from-legacy-ldsc.md](<../../../user-guide/workflows/migrating-from-legacy-ldsc.md>).

Legacy implementation context: [`ldscore/regressions.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/regressions.py), [`ldscore/sumstats.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/sumstats.py).

<a id="issue-37"></a>

### #37: median value issue for munge_sumstats.py

[Original issue](https://github.com/bulik/ldsc/issues/37) · Source state: **closed** · Severity: **moderate** · Accessible comments reviewed: **14**.

Workflow tags: `munge-sumstats`.

**Trigger and harm.** Raw quantitative beta coefficients have a median outside the fixed 0.1 tolerance; recent comments also involve N=629 and skewed effects. A valid unit rescaling can pass or fail preprocessing without changing signs or P values.

**Legacy disposition.** An early user confirms division by phenotype SD passed; later discussion warns low sample size/outliers and does not endorse median centering.

**LDSC3 assessment — Partial coverage.** LDSC3 retains check_median with absolute tolerance and therefore the same scale sensitivity.
Its actionable message recommends verifying column/null labels, but that does not address valid beta units.
Scaling must preserve the GWAS statistic and allele sign; arbitrary median centering can change signs and invalidate genetic correlation.

**Recommended prevention.** Explain why the guard exists and how to verify units or provide signed Z; consider a scale-aware diagnostic after validating it, rather than encouraging recentering.

**Evidence and verification.** source-inspected
[src/ldsc/_sumstats/munge_kernel.py:check_median](<../../../../src/ldsc/_sumstats/munge_kernel.py>); [src/ldsc/_sumstats/munge_input.py](<../../../../src/ldsc/_sumstats/munge_input.py>); [Source](<https://github.com/bulik/ldsc/issues/37#issuecomment-3968866251>).

Legacy implementation context: [`munge_sumstats.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/munge_sumstats.py).

<a id="issue-65"></a>

### #65: Problems with external data source

[Original issue](https://github.com/bulik/ldsc/issues/65) · Source state: **closed** · Severity: **moderate** · Accessible comments reviewed: **3**.

Workflow tags: `munge-sumstats`.

**Trigger and harm.** A personality GWAS column called OR is centered near zero and its downloaded file lacks headers.
Wrong null value or an added header with too many fields can remove every row or reverse interpretation.

**Legacy disposition.** Maintainer suspects header/data mismatch and empty columns; no input-level resolution is documented.

**LDSC3 assessment — Partial coverage.** LDSC3 accepts explicit signed-statistic column/null settings and validates required numeric fields, but cannot determine biological meaning from a misleading OR label.
Header inference must remain reviewable. A syntactically valid file may still encode beta under an odds-ratio name, so source documentation is needed.

**Recommended prevention.** Show how to inspect inferred columns and override a verified null; do not transform values merely to satisfy a median guard.
External image was not read as an independent data source.

**Evidence and verification.** source-inspected
[src/ldsc/_sumstats/munge_input.py](<../../../../src/ldsc/_sumstats/munge_input.py>); [src/ldsc/_sumstats/munge_kernel.py:check_median](<../../../../src/ldsc/_sumstats/munge_kernel.py>); [docs/reference/munge-sumstats.md](<../../../reference/munge-sumstats.md>).

Legacy implementation context: [`munge_sumstats.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/munge_sumstats.py).

<a id="issue-131"></a>

### #131: Error converting summary statistics

[Original issue](https://github.com/bulik/ldsc/issues/131) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **3**.

Workflow tags: `munge-sumstats`.

**Trigger and harm.** A minimal report expands into unrecognized headers, unsupported compression and an rsID converted as a signed statistic.
Users may coerce identifier columns to numbers instead of fixing schema alignment.

**Legacy disposition.** Maintainer explains allowed aliases, whitespace assumptions and decompression; user reports a float conversion worked without enough detail to prove the cause.

**LDSC3 assessment — Partial coverage.** LDSC3 header inference, explicit mapping and numeric validation improve diagnosis, but unknown delimiter/header mismatch must still be resolved from original source rows.
Numeric coercion cannot make an rsID a valid beta.
The issue lacks a complete original input, so resolution is not evidence that arbitrary conversion is safe.

**Recommended prevention.** Surface the inferred map and first offending row/column; document supported compression/delimiters and avoid recommending wholesale float conversion.

**Evidence and verification.** source-inspected
[src/ldsc/_sumstats/munge_input.py](<../../../../src/ldsc/_sumstats/munge_input.py>); [src/ldsc/_sumstats/munge_kernel.py:parse_dat](<../../../../src/ldsc/_sumstats/munge_kernel.py>); [docs/reference/formats/column-schema.md](<../../../reference/formats/column-schema.md>).

Legacy implementation context: [`munge_sumstats.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/munge_sumstats.py).

<a id="issue-134"></a>

### #134: Why large Z score are removed

[Original issue](https://github.com/bulik/ldsc/issues/134) · Source state: **closed** · Severity: **moderate** · Accessible comments reviewed: **1**.

Workflow tags: `munge-sumstats`.

**Trigger and harm.** User interprets loss of high-Z SNPs as a statistic-magnitude filter.
They may disable quality filters or append only significant variants.

**Legacy disposition.** Maintainer lists INFO/MAF/N/allele/restriction reasons and explains representative genome-wide coverage.

**LDSC3 assessment — Partial coverage.** LDSC3 logs stage-specific removals and has explicit HM3/custom restriction options.
Its dropped-SNP sidecar does not include every early QC-removal row, so users cannot always attribute one particular hit using the final sidecar alone.
Preserving every strong SNP is not the objective of a calibrated regression-SNP set.

**Recommended prevention.** Explain the actual reason hierarchy and how to inspect a specific rejected SNP; consider targeted rejection inspection rather than exporting all raw rows by default.

**Evidence and verification.** source-inspected
[src/ldsc/_sumstats/munge_kernel.py:parse_dat](<../../../../src/ldsc/_sumstats/munge_kernel.py>); [src/ldsc/_sumstats/munge_kernel.py:munge_sumstats](<../../../../src/ldsc/_sumstats/munge_kernel.py>); [docs/reference/munge-sumstats.md](<../../../reference/munge-sumstats.md>).

Legacy implementation context: [`munge_sumstats.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/munge_sumstats.py).

<a id="issue-145"></a>

### #145: FIXED - Fail to converst summary statistic in .sumstats format: munge_sumstats is taking hours

[Original issue](https://github.com/bulik/ldsc/issues/145) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **13**.

Workflow tags: `munge-sumstats`.

**Trigger and harm.** A five-million-row chunk stalls for hours while 500,000 rows finishes quickly.
Users mistake a performance problem for corrupted data or installation failure.

**Legacy disposition.** Many independent comments confirm smaller --chunksize helps; no profile or universal causal fix is provided.

**LDSC3 assessment — Partial coverage.** LDSC3 retains a chunk-size control with a lower one-million-row default and stage progress.
It still combines retained chunks for global operations, so smaller input chunks do not guarantee bounded total memory.
The pasted example also uses a schizophrenia N/output label for a bipolar file; performance success does not validate these scientific settings.

**Recommended prevention.** Document a smaller-chunk troubleshooting example and total-memory limits, then profile a repeatable workload before adding adaptive chunking or claiming a speedup.

**Evidence and verification.** source-inspected
[src/ldsc/config.py:MungeConfig](<../../../../src/ldsc/config.py>); [src/ldsc/_sumstats/munge_kernel.py:parse_dat](<../../../../src/ldsc/_sumstats/munge_kernel.py>); [docs/reference/munge-sumstats.md](<../../../reference/munge-sumstats.md>).

Legacy implementation context: [`munge_sumstats.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/munge_sumstats.py).

<a id="issue-150"></a>

### #150: munge_sumstats.py throwing error on simulated dataset

[Original issue](https://github.com/bulik/ldsc/issues/150) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **2**.

Workflow tags: `munge-sumstats`.

**Trigger and harm.** Every displayed simulated record uses the strand-ambiguous pair T/A and invented rsIDs/coordinates. A concat error is blamed on simulation identifiers while another filter may remove every row.

**Legacy disposition.** No proven resolution; comment speculates about an absent merge list and another user bypasses the program.

**LDSC3 assessment — Partial coverage.** LDSC3 permits explicit rsID base mode and no-SNP-restriction for simulations, while coordinate-family inputs must satisfy supported genome-build policy.
The displayed all-T/A input would fail ordinary allele-aware ambiguity requirements independently of made-up rsIDs.
Stage-specific empty diagnostics are needed to distinguish those causes.

**Recommended prevention.** Use a synthetic fixture with nonambiguous alleles and deliberate restriction/build settings; report the actual emptying stage rather than infer a hidden default from unrelated logs.

**Evidence and verification.** source-inspected
[src/ldsc/_sumstats/munge_kernel.py:parse_dat](<../../../../src/ldsc/_sumstats/munge_kernel.py>); [src/ldsc/_sumstats/munge_kernel.py:_raise_sumstats_restriction_empty](<../../../../src/ldsc/_sumstats/munge_kernel.py>); [docs/reference/munge-sumstats.md](<../../../reference/munge-sumstats.md>).

Legacy implementation context: [`munge_sumstats.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/munge_sumstats.py).

<a id="issue-190"></a>

### #190: Make_annot.py issue and windows size to choose

[Original issue](https://github.com/bulik/ldsc/issues/190) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **2**.

Workflow tags: `ldscore`, `ldsc-seg`.

**Trigger and harm.** MSigDB identifiers, an Ensembl coordinate table and100kb padding yield only zero annotation values. A valid-looking empty query may hide identifier/build mismatch and cannot test the intended biology.

**Legacy disposition.** Thread remains open; an external rewritten script is offered without a demonstrated root-cause resolution.

**LDSC3 assessment — Partial coverage.** LDSC3 resolves gene identifiers through a validated gene coordinate catalog and records unmatched/ambiguous entries.
This improves detection of naming failures, but a resolved interval can still legitimately contain no annotation-row SNPs; padding is a scientific choice rather than a repair for every empty query.

**Recommended prevention.** Report matched genes and covered reference SNPs separately; explain identifier namespaces, genome build and biologically justified padding.

**Evidence and verification.** source-inspected
[src/ldsc/_annotation/gene_lists.py:GeneCatalog](<../../../../src/ldsc/_annotation/gene_lists.py>); [src/ldsc/_annotation/interval_projection.py:ChromosomeProjector.project](<../../../../src/ldsc/_annotation/interval_projection.py>); [docs/reference/annotate.md](<../../../reference/annotate.md>).

Legacy implementation context: [`make_annot.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/make_annot.py), [`ldscore/regressions.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/regressions.py).

<a id="issue-253"></a>

### #253: ValueError: No objects to concatenate

[Original issue](https://github.com/bulik/ldsc/issues/253) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **0**.

Workflow tags: `munge-sumstats`.

**Trigger and harm.** A screenshot and brief text report failure while converting GWAS summary statistics.
The user cannot identify whether all rows were removed or input was malformed.

**Legacy disposition.** Open without a textual traceback or reproducible file; no fix attributable.

**LDSC3 assessment — Partial coverage.** The screenshot was inspected: legacy parse_dat calls pandas.concat on an empty chunk list.
LDSC3 returns an empty frame from parse_dat and raises an explicit retained-0-of-N message with removal counts, apparently preventing that cryptic empty-concatenation error.
The original input is unavailable, so the cause of retaining no chunks remains unproved.

**Recommended prevention.** Record insufficient original evidence and point users to per-stage QC counts and a small textual reproduction.

**Evidence and verification.** Original screenshot visually inspected; current source inspected; original file not reproduced.
[src/ldsc/_sumstats/munge_input.py:prepare_munge_input](<../../../../src/ldsc/_sumstats/munge_input.py>); [src/ldsc/_sumstats/munge_kernel.py:munge_sumstats](<../../../../src/ldsc/_sumstats/munge_kernel.py>).

Legacy implementation context: [`munge_sumstats.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/munge_sumstats.py).

<a id="issue-281"></a>

### #281: mungestats takes a long time

[Original issue](https://github.com/bulik/ldsc/issues/281) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **4**.

Workflow tags: `munge-sumstats`, `ldscore`.

**Trigger and harm.** Users confuse --chunksize500000 with the LD --chunk-size default50.
Applying a parsing parameter to an LD matrix block can cause excessive allocation.

**Legacy disposition.** Users report identical munged outputs at1e5,5e5and1e6 rows; no universal performance fix.

**LDSC3 assessment — Partial coverage.** LDSC3 retains a raw-input --chunksize and separates LD computation settings.
Global N filtering happens after chunk assembly, helping chunk invariance, but input dtype inference and whole-table memory still warrant checks.
The existing user comparison supports these data, not every pathological input.

**Recommended prevention.** Document each parameter unit and purpose and show the established smaller-munging-chunk remedy without changing scientific filters.

**Evidence and verification.** source-inspected
[src/ldsc/_sumstats/munge_kernel.py:parse_dat](<../../../../src/ldsc/_sumstats/munge_kernel.py>); [src/ldsc/_sumstats/munger.py:1417](<../../../../src/ldsc/_sumstats/munger.py#L1417>).

Legacy implementation context: [`munge_sumstats.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/munge_sumstats.py), [`ldscore/ldscore.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/ldscore.py).

<a id="issue-297"></a>

### #297: How to use own gene set obtained from single cell data to do cell type specific analysis?

[Original issue](https://github.com/bulik/ldsc/issues/297) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **0**.

Workflow tags: `ldsc-seg`.

**Trigger and harm.** The title asks how to test a gene list obtained from single-cell data; body is empty.
Users need a path from biological selection to a conditional query, but exact design is unknown.

**Legacy disposition.** Open with no details or response.

**LDSC3 assessment — Partial coverage.** LDSC3 supports gene-list query sources, a validated gene coordinate catalog and explicit control gene lists.
This supplies the mechanics requested by the title; it does not establish how cell specificity, top-gene selection or an assay-matched control should be defined from unspecified data.

**Recommended prevention.** Link a short end-to-end gene-list example and state that list/control construction must follow the scientific question.

**Evidence and verification.** source-inspected
[docs/user-guide/workflows/gene-list-ld-scores.md](<../../../user-guide/workflows/gene-list-ld-scores.md>); [docs/user-guide/workflows/cell-type-enrichment.md](<../../../user-guide/workflows/cell-type-enrichment.md>).

Legacy implementation context: [`ldscore/regressions.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/regressions.py).

<a id="issue-316"></a>

### #316: IOError: CRC check failed 0xf157293c != 0xb6a8ec8fL

[Original issue](https://github.com/bulik/ldsc/issues/316) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **0**.

Workflow tags: `munge-sumstats`.

**Trigger and harm.** A user sees an unexplained CRC check failed message while reading compressed data.
File corruption is confused with an analysis/model failure.

**Legacy disposition.** Open without reproducible bytes, named input or a verified fix.

**LDSC3 assessment — Partial coverage.** LDSC3 continues to rely on standard gzip readers for legacy text, so corrupted compressed bytes remain unreadable; accepting them would be unsafe.
Native Parquet avoids that specific wrapper but not general integrity failures.
The improvement needed is file-specific explanatory context, not suppressing the check.

**Recommended prevention.** Explain CRC as archive corruption/truncation, identify the file and recommend integrity verification or reacquisition from the source.

**Evidence and verification.** source-inspected
[src/ldsc/_sumstats/munge_input.py:get_compression](<../../../../src/ldsc/_sumstats/munge_input.py>); [src/ldsc/_sumstats/munger.py:_read_curated_sumstats_artifact](<../../../../src/ldsc/_sumstats/munger.py>).

Legacy implementation context: [`munge_sumstats.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/munge_sumstats.py).

<a id="issue-348"></a>

### #348: Getting Nan in the MAF/LD Score Correlation Matrix

[Original issue](https://github.com/bulik/ldsc/issues/348) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **2**.

Workflow tags: `ldscore`.

**Trigger and harm.** Annotation sums to zero on most chromosomes, so LD scores and correlation denominators are zero.
Users can mistake an undefined correlation for an LD algorithm failure or continue with an annotation that has no support.

**Legacy disposition.** A commenter correctly identifies zero annotation variance; no code fix is reported.

**LDSC3 assessment — Partial coverage.** Binary interval projection and query-status diagnostics help expose empty support, while a per-chromosome zero column can be legitimate for a localized annotation.
Global support and model estimability, not per-chromosome positivity, determine usability.

**Recommended prevention.** Explain NaN correlation for constant annotations and show genome-wide support totals; reject only unusable fitted columns rather than requiring every chromosome to contain query SNPs.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[src/ldsc/_annotation/query_outcomes.py](<../../../../src/ldsc/_annotation/query_outcomes.py>); [src/ldsc/_annotation/interval_projection.py](<../../../../src/ldsc/_annotation/interval_projection.py>); [tests/annotation/test_gene_query_policy.py](<../../../../tests/annotation/test_gene_query_policy.py>).

Legacy implementation context: [`make_annot.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/make_annot.py).

<a id="issue-374"></a>

### #374: Why do I keep reporting an error: Could not find SNP column？

[Original issue](https://github.com/bulik/ldsc/issues/374) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **25**.

Workflow tags: `munge-sumstats`.

**Trigger and harm.** Explicit SNP hints still fail because a comma-delimited file is parsed as a whitespace header; other commenters have distinct unresolved inputs.
The missing-column message sends users toward repeated renaming instead of correcting the separator.

**Legacy disposition.** A user confirms replacing commas with whitespace solved their CSV case; other later reports remain unresolved.
Some suggested commands wrongly use --N N.

**LDSC3 assessment — Partial coverage.** LDSC3 read_header still splits on whitespace.
Explicit column hints cannot rescue an entire comma-separated header token, so the same CSV usability trap remains unless converted first.
Header mismatch and wrong reference-path errors are separate possibilities.

**Recommended prevention.** Detect a single header token containing commas and give a TSV-conversion instruction; distinguish --N numeric from --N-col column without adding speculative format inference.

**Evidence and verification.** Local CSV/TSV header probe: CSV returns one combined token; equivalent TSV returns seven distinct fields.
Original GWAS not rerun.
[src/ldsc/_sumstats/munge_input.py:97](<../../../../src/ldsc/_sumstats/munge_input.py#L97>); [src/ldsc/_sumstats/munge_input.py:287](<../../../../src/ldsc/_sumstats/munge_input.py#L287>).

Legacy implementation context: [`munge_sumstats.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/munge_sumstats.py).
[Executed header checks V4](verification.md#v4-binary-projection-and-header-delimiters).

<a id="issue-404"></a>

### #404: ，

[Original issue](https://github.com/bulik/ldsc/issues/404) · Source state: **closed** · Severity: **moderate** · Accessible comments reviewed: **1**.

Workflow tags: `munge-sumstats`.

**Trigger and harm.** Munging encounters a file with one parsed field in the header and two in a later row.
The parser fails before producing usable summary statistics.

**Legacy disposition.** Closed without a posted resolution; a later commenter repeats the error.

**LDSC3 assessment — Partial coverage.** LDSC3 retains whitespace-delimited raw-input parsing, with special handling for leading ## metadata.
That does not make arbitrary comma headers, malformed lines or extra preambles valid.
The missing source prevents identifying which case occurred.

**Recommended prevention.** Show detected header/separator and first invalid row in diagnostics; require format repair rather than silently skipping malformed data.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[src/ldsc/_sumstats/munge_input.py:84](<../../../../src/ldsc/_sumstats/munge_input.py#L84>); [src/ldsc/_sumstats/munge_input.py:97](<../../../../src/ldsc/_sumstats/munge_input.py#L97>); [src/ldsc/_sumstats/munge_kernel.py:709](<../../../../src/ldsc/_sumstats/munge_kernel.py#L709>).

Legacy implementation context: [`munge_sumstats.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/munge_sumstats.py).

<a id="issue-420"></a>

### #420: munge_sumstats.py taking forever to complete analysis

[Original issue](https://github.com/bulik/ldsc/issues/420) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **5**.

Workflow tags: `munge-sumstats`.

**Trigger and harm.** Legacy SNP membership filtering runs for days with the default huge chunk, but 50,000-row chunks complete quickly.
An unexpectedly pathological preprocessing step wastes compute and appears hung.

**Legacy disposition.** Two users confirm the smaller-chunk workaround; no diagnosed library-level root cause or upstream commit fix is supplied.

**LDSC3 assessment — Partial coverage.** LDSC3 still performs chunkwise restriction but builds explicit identity sets and uses modern pandas.
The original large-file timing was not reproduced, and final retained-table concatenation can still set a memory floor.

**Recommended prevention.** Document munging chunk size separately from LD-score SNP batch size, report progress by stage, and benchmark the membership step before choosing a more elaborate engine.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[src/ldsc/_sumstats/munge_kernel.py:369](<../../../../src/ldsc/_sumstats/munge_kernel.py#L369>); [src/ldsc/_sumstats/munge_kernel.py:424](<../../../../src/ldsc/_sumstats/munge_kernel.py#L424>); [src/ldsc/_sumstats/munge_kernel.py:697](<../../../../src/ldsc/_sumstats/munge_kernel.py#L697>).

Legacy implementation context: [`munge_sumstats.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/munge_sumstats.py).

<a id="faq-03"></a>

### FAQ-03: Historical time and memory estimates

[FAQ source](https://github.com/bulik/ldsc/wiki/FAQ) · Source-order identifier: **FAQ-03** · Severity: **moderate**.

Workflow tags: `ldscore`.

**Legacy disposition.** Historical FAQ guidance; assessed independently rather than treated as a code-fix record.

**LDSC3 assessment — Partial coverage.** The historical hardware estimate is not a performance guarantee for LDSC3.
Input size, LD window, backend, annotation count and output batching govern the current workload; the repository already records bounded performance evidence.

**Recommended prevention.** Give representative measured configurations and limits, with a small pilot recommendation, rather than copying an old runtime number.

**Evidence and verification.** Source-inspected; relevant executed checks are identified in verification.md.
[docs/dev/validation/performance.md](<../performance.md>); [docs/dev/design/memory-and-cpu.md](<../../design/memory-and-cpu.md>).

<a id="issue-57"></a>

### #57: LD Hub website documentation issues

[Original issue](https://github.com/bulik/ldsc/issues/57) · Source state: **closed** · Severity: **moderate** · Accessible comments reviewed: **5**.

Workflow tags: `munge-sumstats`, `rg`.

**Trigger and harm.** Web uploads expire, schemas are unclear, resource matrices differ in age, and support ownership is unclear.
Users cannot tell input rejection, service failure and outdated scientific data apart.

**Legacy disposition.** Legacy maintainers only fixed support/contact documentation via PR #60; they explicitly did not own the LD Hub service.

**LDSC3 assessment — Intentionally unsupported.** LDSC3 is a local package and does not inherit web sessions or upload code.
Transferable lessons apply to explicit schema inference, progress, durable failure/status records, trait labels, and release provenance. A new upload service or data portal would be a large unrelated feature, not a remedy for the local package.

**Recommended prevention.** Keep local diagnostics self-contained and distinguish support for the package from support for external services; do not add a web front end to address this issue.

**Evidence and verification.** source-inspected
[src/ldsc/_runtime/progress.py](<../../../../src/ldsc/_runtime/progress.py>); [src/ldsc/_sumstats/munge_input.py](<../../../../src/ldsc/_sumstats/munge_input.py>); [src/ldsc/_regression/runner.py:_disambiguate_trait_names](<../../../../src/ldsc/_regression/runner.py>); [Source](<https://github.com/bulik/ldsc/pull/60>).

Legacy implementation context: [`munge_sumstats.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/munge_sumstats.py), [`ldscore/regressions.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/regressions.py), [`ldscore/sumstats.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/sumstats.py).

<a id="issue-237"></a>

### #237: --overlap-annot does not work when reading in "thin" annot files?

[Original issue](https://github.com/bulik/ldsc/issues/237) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **0**.

Workflow tags: `baseline-partitioned-h2`.

**Trigger and harm.** Legacy --overlap-annot reads thin files without CHR/BP/CM and fails.
Users receive apparently contradictory thin-file recommendations between LD-score construction and regression.

**Legacy disposition.** Reporter succeeds after restoring metadata from the matching BIM; no implementation fix verified.

**LDSC3 assessment — Intentionally unsupported.** LDSC3 deliberately requires metadata-bearing annotations and stores computed overlap in the LD-score directory.
Regression no longer rereads thin files for overlap, but legacy thin-file generation/conversion is unsupported rather than transparently compatible.

**Recommended prevention.** State the format boundary and provide a validated metadata-restoration recipe using the original BIM order, then generate native LD-score directories.

**Evidence and verification.** source-inspected
[docs/reference/results/partitioned-h2.md](<../../../reference/results/partitioned-h2.md>); [src/ldsc/_regression/runner.py:summarize_partitioned_h2](<../../../../src/ldsc/_regression/runner.py>).

Legacy implementation context: [`ldscore/regressions.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/regressions.py).

<a id="issue-279"></a>

### #279: While calculating LD scores over annotation file got IndexError: positional indexers out of bounds

[Original issue](https://github.com/bulik/ldsc/issues/279) · Source state: **closed** · Severity: **moderate** · Accessible comments reviewed: **1**.

Workflow tags: `ldscore`.

**Trigger and harm.** make_annot emits a one-column ANNOT file but ldscore is invoked without --thin-annot. A positional-index error gives no guidance about the actual format mismatch.

**Legacy disposition.** Reporter adds --thin-annot and confirms success; no implementation fix.

**LDSC3 assessment — Intentionally unsupported.** LDSC3 deliberately accepts metadata-bearing annotation files and rejects thin files rather than inheriting this flag-dependent ambiguity.
Its annotate writes the annotation SNP grid metadata, and header-based resolution avoids depending on a fixed positional CHR/SNP/BP column count.

**Recommended prevention.** Explain the unsupported thin format and generate full native annotation files from the original grid.

**Evidence and verification.** source-inspected
[src/ldsc/_annotation/parsing.py:normalize_annotation_chunk](<../../../../src/ldsc/_annotation/parsing.py>); [docs/reference/annotate.md](<../../../reference/annotate.md>); [docs/reference/legacy-flag-map.md](<../../../reference/legacy-flag-map.md>).

Legacy implementation context: [`make_annot.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/make_annot.py).

<a id="issue-478"></a>

### #478: ldsc python 3  H2 estimation question

[Original issue](https://github.com/bulik/ldsc/issues/478) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **1**.

Workflow tags: `h2`.

**Trigger and harm.** A posted command contains ../ ADHD.susmstats.gz and missing continuation while attributing an AttributeError to Python3.
Shell tokenization can hide the actual input and confound comparisons across machines.

**Legacy disposition.** A commenter points out the embedded space; no traceback or successful verification is provided.

**LDSC3 assessment — Indeterminate.** LDSC3 Python3 support does not establish the cause of this unspecified AttributeError.
Its explicit input path checks help when parsing reaches the workflow, but cannot repair an incorrectly split shell command.

**Recommended prevention.** Use a quoted one-line command with a full traceback and interpreter/source provenance; do not diagnose a Python-version regression from this record.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[src/ldsc/_runtime/paths.py:182](<../../../../src/ldsc/_runtime/paths.py#L182>); [src/ldsc/cli.py](<../../../../src/ldsc/cli.py>); [docs/user-guide/workflows/heritability.md](<../../../user-guide/workflows/heritability.md>).

Legacy implementation context: [`ldscore/regressions.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/regressions.py).

<a id="issue-13"></a>

### #13: document .sumstats file format (replacement for .chisq)

[Original issue](https://github.com/bulik/ldsc/issues/13) · Source state: **closed** · Severity: **moderate** · Accessible comments reviewed: **1**.

Workflow tags: `munge-sumstats`, `h2`, `rg`.

**Trigger and harm.** Users cannot tell which fields replace the obsolete .chisq format.
Incorrect headers, signed alleles or sample sizes propagate into regression.

**Legacy disposition.** Maintainer states the .sumstats format was documented in the wiki; issue closed.

**LDSC3 assessment — Apparently prevented.** LDSC3 has a dedicated legacy-summary-statistics format page and a self-describing Parquet output.
Legacy interoperability requires SNP/A1/A2/Z/N, while Parquet retains SNP identifier mode and genome-build metadata.
The user guide separates raw inputs from munged outputs and explains that A1 is the allele for the signed statistic.

**Recommended prevention.** Keep one short input/output example near installation and link it from validation errors; avoid duplicating schemas across many pages.

**Evidence and verification.** source-inspected
[docs/reference/formats/legacy-sumstats.md](<../../../reference/formats/legacy-sumstats.md>); [docs/reference/munge-sumstats.md](<../../../reference/munge-sumstats.md>); [src/ldsc/_sumstats/munger.py](<../../../../src/ldsc/_sumstats/munger.py>).

Legacy implementation context: [`munge_sumstats.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/munge_sumstats.py), [`ldscore/regressions.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/regressions.py), [`ldscore/sumstats.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/sumstats.py).

<a id="issue-23"></a>

### #23: Can't find munge_sumstats.py script

[Original issue](https://github.com/bulik/ldsc/issues/23) · Source state: **closed** · Severity: **moderate** · Accessible comments reviewed: **4**.

Workflow tags: `munge-sumstats`, `rg`.

**Trigger and harm.** Tutorial names a script not present on the downloaded branch, followed by an incompatible reference-resource download. A beginner can execute an internal module without producing results or combine the wrong resource format.

**Legacy disposition.** Maintainer identified dev-versus-master documentation mismatch, supplied an appropriate EUR LD-score URL, and the user confirmed success.

**LDSC3 assessment — Apparently prevented.** LDSC3 exposes one installed CLI and documents it with its own checkout.
Installation smoke checks exercise ldsc and python -m ldsc.
Resource conversion is explicit, so the historical RDS file is not a regression input.
This remains a release-documentation maintenance obligation rather than a reason to add script aliases.

**Recommended prevention.** Tie documentation examples to a release/revision and distinguish reference LD scores, weights and source genotypes at the download point.

**Evidence and verification.** source-inspected
[src/ldsc/cli.py](<../../../../src/ldsc/cli.py>); [docs/user-guide/installation.md](<../../../user-guide/installation.md>); [docs/user-guide/reference-data.md](<../../../user-guide/reference-data.md>); [tests/package/test_cli_help.py](<../../../../tests/package/test_cli_help.py>).

Legacy implementation context: [`munge_sumstats.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/munge_sumstats.py), [`ldscore/regressions.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/regressions.py), [`ldscore/sumstats.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/sumstats.py).

<a id="issue-76"></a>

### #76: only p values no beta

[Original issue](https://github.com/bulik/ldsc/issues/76) · Source state: **closed** · Severity: **moderate** · Accessible comments reviewed: **2**.

Workflow tags: `munge-sumstats`, `rg`.

**Trigger and harm.** Only two-sided P values and numeric direction indicators are available, without beta magnitude.
Users may unnecessarily predict betas from another study or discard valid correlation inputs.

**Legacy disposition.** Maintainer recommends numeric +1/-1 as signed statistic with null zero; user later closes the issue.

**LDSC3 assessment — Apparently prevented.** LDSC3 converts P to Z magnitude and uses the selected statistic relative to its null for sign.
Numeric direction therefore suffices when it is correctly aligned to A1 and P represents the intended one-degree-of-freedom test.
It does not make unsigned P alone usable for genetic correlation or validate arbitrary multi-degree tests.

**Recommended prevention.** Add a concise numeric-direction example and verify allele orientation rather than introducing a beta-imputation feature.

**Evidence and verification.** source-inspected
[src/ldsc/_sumstats/munge_kernel.py:p_to_z](<../../../../src/ldsc/_sumstats/munge_kernel.py>); [src/ldsc/_sumstats/munge_kernel.py:munge_sumstats](<../../../../src/ldsc/_sumstats/munge_kernel.py>); [docs/reference/munge-sumstats.md](<../../../reference/munge-sumstats.md>).

Legacy implementation context: [`munge_sumstats.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/munge_sumstats.py), [`ldscore/regressions.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/regressions.py), [`ldscore/sumstats.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/sumstats.py).

<a id="issue-83"></a>

### #83: Issues with ldscores from 2 LDSC_Seg annotations

[Original issue](https://github.com/bulik/ldsc/issues/83) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **2**.

Workflow tags: `ldsc-seg`.

**Trigger and harm.** An edited .ldcts file no longer has two fields per line.
Cell-type analysis fails with an opaque unpacking error.

**Legacy disposition.** User confirms reformatting the path manifest fixed Franke/ImmGen inputs.

**LDSC3 assessment — Apparently prevented.** LDSC3 does not consume the old .ldcts manifest; it builds query annotations and one LD-score directory, with path/header preflight before fitting.
This prevents the specific tuple-unpacking path.
Migrating users still need an explicit explanation that legacy cell-type manifest files are not accepted directly.

**Recommended prevention.** Show how to convert or regenerate a query LD-score workflow; reject a .ldcts file with an objective-oriented migration message.

**Evidence and verification.** source-inspected
[docs/reference/legacy-flag-map.md:124](<../../../reference/legacy-flag-map.md#established-command-and-input-renames>); [src/ldsc/_ldscore/preflight.py](<../../../../src/ldsc/_ldscore/preflight.py>); [docs/user-guide/workflows/cell-type-enrichment.md](<../../../user-guide/workflows/cell-type-enrichment.md>).

Legacy implementation context: [`ldscore/regressions.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/regressions.py).

<a id="issue-110"></a>

### #110: partitioning heritability

[Original issue](https://github.com/bulik/ldsc/issues/110) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **0**.

Workflow tags: `ldscore`, `baseline-partitioned-h2`.

**Trigger and harm.** User expects annotation, predictor and weight files to have identical row counts.
They may shrink reference SNPs to output SNPs and change the LD model.

**Legacy disposition.** Open without an answer.

**LDSC3 assessment — Apparently prevented.** LDSC3 explicitly separates reference SNPs contributing to scores/counts/overlap from filtered regression SNPs appearing as LD-score rows.
The native directory keeps predictor and weight alignment together.
Annotation and score tables need not have the same number of rows; reference and regression sets have different purposes.

**Recommended prevention.** Put a small two-set example in the beginner workflow and explain identity joins rather than requiring equal raw row counts.

**Evidence and verification.** source-inspected
[docs/reference/ldscore-snp-universes.md](<../../../reference/ldscore-snp-universes.md>); [src/ldsc/_ldscore/compute.py:compute_chromosome](<../../../../src/ldsc/_ldscore/compute.py>).

Legacy implementation context: [`make_annot.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/make_annot.py), [`ldscore/regressions.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/regressions.py).

<a id="issue-111"></a>

### #111: To analyze a new annotation, you can add it to the full baseline model.

[Original issue](https://github.com/bulik/ldsc/issues/111) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **0**.

Workflow tags: `ldscore`, `ldsc-seg`.

**Trigger and harm.** User asks whether to append query columns physically to the full baseline annotation file.
They may overwrite distributed baselines or fit the query without its controls.

**Legacy disposition.** Open without a maintainer answer.

**LDSC3 assessment — Apparently prevented.** LDSC3 accepts separate --baseline-annot-sources and --query-annot-sources/BED/gene-list inputs, requires explicit baseline annotations for query scoring, and fits one query plus baseline per model.
Names and SNP identities are validated across inputs; physical concatenation and alteration of canonical baseline files are unnecessary.

**Recommended prevention.** Show separate baseline/query files on the same annotation SNP grid and explain how baseline-only joint modeling differs from a query scan.

**Evidence and verification.** source-inspected
[src/ldsc/_ldscore/calculator.py:run_ldscore_from_args](<../../../../src/ldsc/_ldscore/calculator.py>); [docs/user-guide/workflows/partitioned-heritability.md](<../../../user-guide/workflows/partitioned-heritability.md>).

Legacy implementation context: [`ldscore/ldscore.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/ldscore.py), [`ldscore/regressions.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/regressions.py).

<a id="issue-143"></a>

### #143: Can i use ldsc to calculate quantitative traits.

[Original issue](https://github.com/bulik/ldsc/issues/143) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **0**.

Workflow tags: `munge-sumstats`, `h2`, `rg`.

**Trigger and harm.** User cannot compute an OR for a continuous phenotype.
Artificially converting a linear effect to OR can misstate null values and signed statistics.

**Legacy disposition.** Open with no maintainer answer.

**LDSC3 assessment — Apparently prevented.** LDSC3 accepts beta, log odds, Z or a verified signed statistic with explicit null.
Quantitative traits use their linear association statistics and require no liability conversion.
The pasted hg18 coordinate label is outside the supported hg19/hg38 inference scope and needs source verification/liftover; matching rsIDs alone does not validate the coordinates.

**Recommended prevention.** Show one quantitative beta/SE or signed-Z example, define N and effect allele, and explicitly address unsupported genome-build labels.

**Evidence and verification.** source-inspected
[docs/reference/munge-sumstats.md](<../../../reference/munge-sumstats.md>); [src/ldsc/_sumstats/munge_input.py:null_values](<../../../../src/ldsc/_sumstats/munge_input.py>); [docs/reference/snp-identity-and-genome-build.md](<../../../reference/snp-identity-and-genome-build.md>).

Legacy implementation context: [`munge_sumstats.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/munge_sumstats.py), [`ldscore/regressions.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/regressions.py), [`ldscore/sumstats.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/sumstats.py).

<a id="issue-146"></a>

### #146: Partitioning heritability tutorial --merge-alleles issue

[Original issue](https://github.com/bulik/ldsc/issues/146) · Source state: **closed** · Severity: **moderate** · Accessible comments reviewed: **2**.

Workflow tags: `munge-sumstats`.

**Trigger and harm.** head shows identifiers but empty other fields after a legacy HM3 merge.
Users disable alignment/restriction even though most retained data are later in the file.

**Legacy disposition.** A later user confirms rows are padded to the merge list and scrolling reveals data; no defect is established.

**LDSC3 assessment — Apparently prevented.** LDSC3 writes retained munged rows rather than padding every HM3 record, with explicit n_input/n_output and exclusive removals.
Legacy sumstats loading still handles missing rows according to its contract.
This is an intentional format/UX change that migration guidance should mention, not a requirement to preserve empty placeholders.

**Recommended prevention.** Show an output preview with actual retained counts and explain why LDSC3 row counts differ from legacy HM3-padded files.

**Evidence and verification.** source-inspected
[src/ldsc/_sumstats/munge_kernel.py:munge_sumstats](<../../../../src/ldsc/_sumstats/munge_kernel.py>); [docs/reference/formats/legacy-sumstats.md](<../../../reference/formats/legacy-sumstats.md>); [docs/reference/munge-sumstats.md](<../../../reference/munge-sumstats.md>).

Legacy implementation context: [`munge_sumstats.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/munge_sumstats.py).

<a id="issue-147"></a>

### #147: Cell type analysis: IOError: Could not open Cahoy.1.l2.ldscore[./gz/bz2]

[Original issue](https://github.com/bulik/ldsc/issues/147) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **10**.

Workflow tags: `ldsc-seg`.

**Trigger and harm.** Cahoy prefixes inside .ldcts resolve from the calling directory; later users also report a directory typo and cloud access issues.
Users recompute scores or change software when simple path/access repair is needed.

**Legacy disposition.** Maintainer explains current-directory resolution; one user confirms ldscore versus ld_score typo. A Polyfun workaround is not a verified general fix.

**LDSC3 assessment — Apparently prevented.** LDSC3 uses explicit validated path tokens and saved LD-score directories, and does not consume .ldcts manifests.
Missing input errors should show the resolved concrete path and required companions.
External Cahoy availability remains a resource issue even though path resolution improves.

**Recommended prevention.** Provide migration guidance and distinguish missing file, unreadable file and access-restricted download; do not recommend recomputation before inspecting paths.

**Evidence and verification.** source-inspected
[src/ldsc/_runtime/paths.py](<../../../../src/ldsc/_runtime/paths.py>); [docs/reference/legacy-flag-map.md](<../../../reference/legacy-flag-map.md>); [docs/user-guide/workflows/cell-type-enrichment.md](<../../../user-guide/workflows/cell-type-enrichment.md>).

Legacy implementation context: [`ldscore/regressions.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/regressions.py).

<a id="issue-149"></a>

### #149: Estimate heritability

[Original issue](https://github.com/bulik/ldsc/issues/149) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **7**.

Workflow tags: `h2`.

**Trigger and harm.** The log opens eur_w_ld/1.l2.ldscore although commands mention another directory.
Users undertake unnecessary genotype/LD-score recomputation.

**Legacy disposition.** User says self-computed scores work but provides no root-cause proof; subsequent users remain confused.

**LDSC3 assessment — Apparently prevented.** LDSC3 takes a validated LD-score directory with declared members rather than guessing chromosome suffixes from legacy prefixes.
Its converter has explicit family discovery and missing-file diagnostics.
Successful installation/tests do not establish correct analysis paths, and 128k retained SNPs also require scientific coverage assessment after the path is repaired.

**Recommended prevention.** Diagnose the actual resolved path and manifest first; do not prescribe recomputation or assume the downloaded LD scores are scientifically wrong.

**Evidence and verification.** source-inspected
[src/ldsc/_ldscore/source.py](<../../../../src/ldsc/_ldscore/source.py>); [src/ldsc/_ldscore/legacy_converter.py](<../../../../src/ldsc/_ldscore/legacy_converter.py>); [docs/reference/input-paths.md](<../../../reference/input-paths.md>).

Legacy implementation context: [`ldscore/regressions.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/regressions.py).

<a id="issue-170"></a>

### #170: Empty rows in output after running munge_sumstats.py

[Original issue](https://github.com/bulik/ldsc/issues/170) · Source state: **closed** · Severity: **moderate** · Accessible comments reviewed: **2**.

Workflow tags: `munge-sumstats`.

**Trigger and harm.** Legacy output has all HM3 identifiers but only 1,065,283 populated records.
Users perform avoidable postprocessing or treat normal placeholders as allele errors.

**Legacy disposition.** User correctly identifies output padding, although calls it a recurring problem; no algorithmic fix is established.

**LDSC3 assessment — Apparently prevented.** LDSC3 outputs retained rows and records exclusive drop counts; it does not emit all-HM3 placeholder rows.
Legacy interoperability readers still remove missing required values appropriately.
More output lines do not imply more regression information, and retaining identifiers of excluded SNPs is a presentation convention rather than added scientific observations.

**Recommended prevention.** Document retained-row semantics and show how logs reconcile to data rows; preserve originals rather than edit legacy padded files in place.

**Evidence and verification.** source-inspected
[src/ldsc/_sumstats/munge_kernel.py:munge_sumstats](<../../../../src/ldsc/_sumstats/munge_kernel.py>); [docs/reference/formats/legacy-sumstats.md](<../../../reference/formats/legacy-sumstats.md>).

Legacy implementation context: [`munge_sumstats.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/munge_sumstats.py).

<a id="issue-183"></a>

### #183: IOError could not open

[Original issue](https://github.com/bulik/ldsc/issues/183) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **1**.

Workflow tags: `ldsc-seg`.

**Trigger and harm.** A missing leading slash or a relative path embedded in a cell-type manifest causes IOError.
Users confuse inaccessible input with an LD-score calculation defect.

**Legacy disposition.** Reporter fixes the absolute paths; another changes working directory for manifest entries.
No numerical fix is involved.

**LDSC3 assessment — Apparently prevented.** LDSC3 resolves explicit path tokens and validates missing inputs before fitting; native LD-score directories and query manifests replace the legacy manually paired chromosome prefixes.
Relative paths still need an explicit documented base, especially in copied scripts.

**Recommended prevention.** Show absolute or clearly based paths and name the unresolved token; explain manifest-relative path rules where applicable.

**Evidence and verification.** source-inspected
[src/ldsc/_runtime/paths.py:resolve_file_group](<../../../../src/ldsc/_runtime/paths.py>); [docs/reference/input-paths.md](<../../../reference/input-paths.md>).

Legacy implementation context: [`ldscore/regressions.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/regressions.py).

<a id="issue-192"></a>

### #192: How to merge annotation files

[Original issue](https://github.com/bulik/ldsc/issues/192) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **3**.

Workflow tags: `ldscore`, `baseline-partitioned-h2`, `ldsc-seg`.

**Trigger and harm.** Users want several custom annotations in columns and ask whether to merge files or fit each separately.
Joint annotation adjustment and one-query-per-baseline tests answer different questions.

**Legacy disposition.** Commenters confirm a comma-separated legacy --ref-ld-chr list works; fitting each with baseline is an alternative model, not identical analysis.

**LDSC3 assessment — Apparently prevented.** LDSC3 accepts multiple baseline annotation sources and multiple query sources with explicit roles.
Baseline-only partitioned heritability fits jointly; query annotations are tested individually on the common baseline.
This removes a need for manual pasted LD-score columns while making the model distinction essential.

**Recommended prevention.** Show side-by-side joint-baseline and individual-query examples; align annotations by SNP identity rather than linewise paste.

**Evidence and verification.** source-inspected
[docs/user-guide/workflows/annotations.md](<../../../user-guide/workflows/annotations.md>); [docs/user-guide/workflows/partitioned-heritability.md](<../../../user-guide/workflows/partitioned-heritability.md>); [docs/user-guide/workflows/cell-type-enrichment.md](<../../../user-guide/workflows/cell-type-enrichment.md>).

Legacy implementation context: [`make_annot.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/make_annot.py), [`ldscore/regressions.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/regressions.py).

<a id="issue-238"></a>

### #238: alueError: Too many signed sumstat columns. Specify which to ignore with the --ignore flag.

[Original issue](https://github.com/bulik/ldsc/issues/238) · Source state: **closed** · Severity: **moderate** · Accessible comments reviewed: **2**.

Workflow tags: `munge-sumstats`.

**Trigger and harm.** Raw data contain both OR and Z and munging refuses to choose.
Ambiguous effect columns block input but prevent silent sign misinterpretation.

**Legacy disposition.** Commenter confirms --ignore OR works; reporter acknowledges and closes.

**LDSC3 assessment — Apparently prevented.** LDSC3 requires one direction statistic or an explicit --signed-sumstats choice and suppresses competing default signed aliases when selected.
This is an intentional ambiguity error rather than a numerical defect. P still determines Z magnitude even when the direction column is named Z.

**Recommended prevention.** Prefer --signed-sumstats Z,0 or OR,1 and show the selected mapping in --infer-only output.

**Evidence and verification.** source-inspected
[src/ldsc/_sumstats/munge_input.py:prepare_munge_input](<../../../../src/ldsc/_sumstats/munge_input.py>); [src/ldsc/_sumstats/munge_kernel.py:munge_sumstats](<../../../../src/ldsc/_sumstats/munge_kernel.py>).

Legacy implementation context: [`munge_sumstats.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/munge_sumstats.py).

<a id="issue-240"></a>

### #240: LDSC KeyError: "CHR"

[Original issue](https://github.com/bulik/ldsc/issues/240) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **0**.

Workflow tags: `h2`.

**Trigger and harm.** A legacy --ref-ld/--w-ld input causes a bare KeyError:CHR.
An incomplete or wrongly selected reference file cannot be diagnosed from the message.

**Legacy disposition.** Open without input files or a verified cause/fix.

**LDSC3 assessment — Apparently prevented.** LDSC3 requires a self-describing LD-score directory or explicit legacy conversion with full chromosome families and canonical metadata.
Header resolution and input validation identify missing fields before fitting.
The exact malformed input is unavailable, so no numerical defect should be inferred.

**Recommended prevention.** Name the offending file and required CHR/POS/SNP fields; distinguish LD-score files from SNP lists or raw annotations.

**Evidence and verification.** source-inspected
[src/ldsc/_ldscore/legacy_converter.py](<../../../../src/ldsc/_ldscore/legacy_converter.py>); [src/ldsc/_core/columns.py:resolve_required_column](<../../../../src/ldsc/_core/columns.py>); [docs/reference/formats/legacy-ldscores.md](<../../../reference/formats/legacy-ldscores.md>).

Legacy implementation context: [`ldscore/regressions.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/regressions.py).

<a id="issue-250"></a>

### #250: IndexError: list index out of range

[Original issue](https://github.com/bulik/ldsc/issues/250) · Source state: **closed** · Severity: **moderate** · Accessible comments reviewed: **9**.

Workflow tags: `baseline-partitioned-h2`, `rg`.

**Trigger and harm.** Legacy prefix paths, missing slashes or unextracted archives yield chrs[0] IndexError. A basic input-location problem appears to be an internal regression bug.

**Legacy disposition.** Reporter confirms extracting folders and using correct paths works; comments also identify passing a SNP list as weights.

**LDSC3 assessment — Apparently prevented.** LDSC3 validates file groups and complete legacy suites before conversion/regression, so an empty discovery becomes an actionable input error.
Compressed individual text files can be valid; extracting an archive directory is distinct from uncompressing every .gz member.

**Recommended prevention.** Show archive-versus-member compression, exact paths and weight-file roles; fail with the missing selection instead of indexing an empty list.

**Evidence and verification.** source-inspected
[src/ldsc/_runtime/paths.py:resolve_file_group](<../../../../src/ldsc/_runtime/paths.py>); [docs/reference/input-paths.md](<../../../reference/input-paths.md>).

Legacy implementation context: [`ldscore/regressions.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/regressions.py), [`ldscore/sumstats.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/sumstats.py).

<a id="issue-260"></a>

### #260: IOError: [Errno 2] No such file or directory: '/var/www/ldhub/ldhub_web/static/sumstats/ldhub_100_sumstat_full_HC_comp1.txt.5fed644d-592e-48c9-8c27-addf36e7e3d0.sumstats.gz'

[Original issue](https://github.com/bulik/ldsc/issues/260) · Source state: **closed** · Severity: **moderate** · Accessible comments reviewed: **3**.

Workflow tags: `munge-sumstats`, `h2`.

**Trigger and harm.** LD Hub tries regression after munging failed, so the expected uploaded .sumstats.gz does not exist.
Users see a downstream IOError rather than the actual preprocessing TypeError.

**Legacy disposition.** A suggestion to limit uploads does not resolve a later user; closure alone is not a fix.

**LDSC3 assessment — Apparently prevented.** LDSC3 is a local package rather than the LD Hub service and validates numeric N during raw munging before a returned table is available.
Its separate commands still need scripts that stop on preprocessing failure; creating expected filenames manually would conceal the root cause.

**Recommended prevention.** Document checking the first failed stage and shell exit codes; distinguish LD Hub service support from local package diagnostics.

**Evidence and verification.** source-inspected
[src/ldsc/_sumstats/munge_kernel.py:parse_dat](<../../../../src/ldsc/_sumstats/munge_kernel.py>); [src/ldsc/_runtime/paths.py:resolve_scalar_path](<../../../../src/ldsc/_runtime/paths.py>).

Legacy implementation context: [`munge_sumstats.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/munge_sumstats.py), [`ldscore/regressions.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/regressions.py).

<a id="issue-278"></a>

### #278: annot.gz file not generated when running LD Score Estimation

[Original issue](https://github.com/bulik/ldsc/issues/278) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **0**.

Workflow tags: `ldscore`.

**Trigger and harm.** A custom make_annot script did not create the file subsequently passed as --annot. A downstream missing-file traceback obscures a failed or incomplete preceding command.

**Legacy disposition.** Open without the custom script, its stderr or a demonstrated fix.

**LDSC3 assessment — Apparently prevented.** LDSC3 supplies standalone annotate plus direct BED/gene query forms for ldscore, with output and input preflight.
It cannot guarantee an external custom script succeeds; the annotation file must exist before it is read. A log alone is not a scientific annotation output.

**Recommended prevention.** Use the native annotation workflow and stop a shell pipeline on the first nonzero exit; name missing expected outputs explicitly.

**Evidence and verification.** source-inspected
[docs/reference/annotate.md](<../../../reference/annotate.md>); [src/ldsc/_runtime/paths.py:resolve_scalar_path](<../../../../src/ldsc/_runtime/paths.py>); [docs/user-guide/workflows/annotations.md](<../../../user-guide/workflows/annotations.md>).

Legacy implementation context: [`make_annot.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/make_annot.py).

<a id="issue-291"></a>

### #291: munge_sumstats: --signed-sumstats error

[Original issue](https://github.com/bulik/ldsc/issues/291) · Source state: **closed** · Severity: **moderate** · Accessible comments reviewed: **4**.

Workflow tags: `munge-sumstats`.

**Trigger and harm.** --signed-sumstats b is supplied without the required comma/null value. A raw unpacking traceback makes a simple argument mistake difficult to fix.

**Legacy disposition.** Comments give beta,0 and OR,1; no implementation fix is needed to clarify valid syntax.

**LDSC3 assessment — Apparently prevented.** LDSC3 _column_hints catches the malformed argument and raises a usage error with concrete examples, and CLI help names COLUMN,NULL.
The number denotes the no-effect null, not SE or the desired direction; this distinction matters for OR versus beta.

**Recommended prevention.** Retain actionable examples and parse validation, and explain the statistical null in the raw-input tutorial.

**Evidence and verification.** source-inspected
[src/ldsc/_sumstats/munge_input.py:_column_hints](<../../../../src/ldsc/_sumstats/munge_input.py>); [src/ldsc/_sumstats/munger.py:1357](<../../../../src/ldsc/_sumstats/munger.py#L1357>).

Legacy implementation context: [`munge_sumstats.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/munge_sumstats.py).

<a id="issue-295"></a>

### #295: Can I calculate ld scores from my own reference samples?

[Original issue](https://github.com/bulik/ldsc/issues/295) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **0**.

Workflow tags: `ldscore`.

**Trigger and harm.** User asks which script calculates LD scores from their own samples.
Private implementation scripts are mistaken for the intended public workflow.

**Legacy disposition.** Open without an answer or code fix.

**LDSC3 assessment — Apparently prevented.** LDSC3 provides the public ldscore command and LDScoreCalculator with PLINK input or an R² matrix.
Users must select representative ancestry, QC, LD window and regression SNPs; merely running a command on any samples does not validate the resulting reference.

**Recommended prevention.** Show one minimal custom-panel command with required PLINK extra, map/window assumptions and a reference-SNP funnel.

**Evidence and verification.** source-inspected
[src/ldsc/_refpanel/adapters.py:PlinkRefPanel.prepare_chromosome](<../../../../src/ldsc/_refpanel/adapters.py>); [src/ldsc/_ldscore/calculator.py:LDScoreCalculator](<../../../../src/ldsc/_ldscore/calculator.py>); [docs/user-guide/workflows/ld-scores.md](<../../../user-guide/workflows/ld-scores.md>).

Legacy implementation context: [`ldscore/ldscore.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/ldscore.py).

<a id="issue-309"></a>

### #309: error

[Original issue](https://github.com/bulik/ldsc/issues/309) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **10**.

Workflow tags: `rg`, `h2`, `ldsc-seg`.

**Trigger and harm.** Missing trailing slash/full weight prefix leaves no chromosome files and chrs[0] raises IndexError.
Repeated vague failures push users toward unrelated methods instead of correcting input selection.

**Legacy disposition.** Comments report full prefixes or different resource downloads work; no universal code fix is demonstrated.

**LDSC3 assessment — Apparently prevented.** LDSC3 validates exact/glob/@ selections and native LD-score directories before fitting.
Absolute paths remove unnecessary requirements to put all files in one directory. A correct basename is still necessary, and changing the estimator or interpreter is not a repair for an empty file selection.

**Recommended prevention.** Name the missing prefix/suite and show one matching filename; distinguish directories, chromosome prefixes and embedded manifest paths.

**Evidence and verification.** source-inspected
[src/ldsc/_runtime/paths.py:resolve_file_group](<../../../../src/ldsc/_runtime/paths.py>); [docs/reference/input-paths.md](<../../../reference/input-paths.md>).

Legacy implementation context: [`ldscore/regressions.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/regressions.py), [`ldscore/sumstats.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/sumstats.py).

<a id="issue-333"></a>

### #333: problems on --merge-alleles

[Original issue](https://github.com/bulik/ldsc/issues/333) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **1**.

Workflow tags: `munge-sumstats`.

**Trigger and harm.** A custom SNP/allele list agrees with input but rows are removed as non-SNP or strand-ambiguous.
Users assume membership matching should disable allele validity checks.

**Legacy disposition.** Open; author clarifies no indels, with no demonstrated fix.

**LDSC3 assessment — Apparently prevented.** LDSC3 distinguishes SNP restriction, identity validity and pairwise orientation. A/T or C/G pairs can be strand-ambiguous even when identical strings occur in both lists; invalid/repeated alleles are another reason.
Matching counts are not a guarantee that a signed cross-study comparison is identifiable.

**Recommended prevention.** Explain separate filter reasons with counts/examples and do not bypass ambiguity solely because a supplied list agrees.

**Evidence and verification.** source-inspected
[src/ldsc/_regression/runner.py:_orient_sumstats_z_to_reference_alleles](<../../../../src/ldsc/_regression/runner.py>); [src/ldsc/_regression/estimators.py:_filter_alleles](<../../../../src/ldsc/_regression/estimators.py>); [docs/reference/formats/legacy-sumstats.md](<../../../reference/formats/legacy-sumstats.md>).

Legacy implementation context: [`munge_sumstats.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/munge_sumstats.py).

<a id="issue-335"></a>

### #335: `--n-min` default in args description does not reflect underlying code

[Original issue](https://github.com/bulik/ldsc/issues/335) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **0**.

Workflow tags: `munge-sumstats`.

**Trigger and harm.** Legacy help says90th-percentile N divided by2 while implementation divides by1.5.
Users cannot reproduce filtering and may unintentionally change retained SNPs.

**Legacy disposition.** Open with precise documentation correction but no linked legacy fix.

**LDSC3 assessment — Apparently prevented.** LDSC3 help explicitly states division by1.5 and the kernel uses that value.
It also documents zero as selecting the default rule and assigns constant N after this filter.
This resolves the contract mismatch without altering the scientific default.

**Recommended prevention.** Keep help/reference/default-contract checks synchronized and report the realized N threshold in provenance.

**Evidence and verification.** source-inspected
[src/ldsc/_sumstats/munger.py:1283](<../../../../src/ldsc/_sumstats/munger.py#L1283>); [src/ldsc/_sumstats/munge_kernel.py:process_n](<../../../../src/ldsc/_sumstats/munge_kernel.py>).

Legacy implementation context: [`munge_sumstats.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/munge_sumstats.py).

<a id="issue-349"></a>

### #349: Cell type specific analyses index error

[Original issue](https://github.com/bulik/ldsc/issues/349) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **9**.

Workflow tags: `ldsc-seg`.

**Trigger and harm.** Legacy cell-type analysis reaches chrs[0] when a weight or annotation prefix finds no files.
An opaque IndexError delays routine path repair and obscures which resource is absent.

**Legacy disposition.** Comments suggest absolute paths; the traceback specifically occurs while reading weights, so the ldcts explanation is not proven for the original command.

**LDSC3 assessment — Apparently prevented.** LDSC3 replaces prefix-based cts loading with an LD-score directory and gives zero-match path diagnostics.
Relative paths still depend on the working directory and need clear examples.

**Recommended prevention.** Keep input-specific zero-match errors and show complete resource paths; avoid claiming that editing ldcts alone solves this traceback.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[src/ldsc/_runtime/paths.py:116](<../../../../src/ldsc/_runtime/paths.py#L116>); [src/ldsc/_runtime/paths.py:202](<../../../../src/ldsc/_runtime/paths.py#L202>); [docs/reference/input-paths.md](<../../../reference/input-paths.md>).

Legacy implementation context: [`ldscore/regressions.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/regressions.py).

<a id="issue-366"></a>

### #366: IndexError: list index out of range

[Original issue](https://github.com/bulik/ldsc/issues/366) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **13**.

Workflow tags: `h2`, `ldsc-seg`.

**Trigger and harm.** A missing trailing separator, typo, unavailable resource, or wrong compression makes legacy chromosome discovery empty.
Many users spend hours on a bare chrs[0] IndexError and unrelated suggested repairs.

**Legacy disposition.** Multiple comments identify path/resource workarounds; no code fix is established.
Later cell-type examples have the same discovery symptom at a different input.

**LDSC3 assessment — Apparently prevented.** LDSC3 reports unmatched path tokens and converts legacy suites before regression.
This removes the precise chrs[0] mechanism but does not make arbitrary renamed or miscompressed inputs valid.

**Recommended prevention.** Retain separate diagnostics for missing paths, unreadable compression, incomplete chromosome suites and schema errors; provide one exact working conversion example.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[src/ldsc/_runtime/paths.py:99](<../../../../src/ldsc/_runtime/paths.py#L99>); [src/ldsc/_ldscore/legacy_converter.py](<../../../../src/ldsc/_ldscore/legacy_converter.py>); [docs/reference/input-paths.md](<../../../reference/input-paths.md>).

Legacy implementation context: [`ldscore/regressions.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/regressions.py).

<a id="issue-392"></a>

### #392: IndexError: list index out of range

[Original issue](https://github.com/bulik/ldsc/issues/392) · Source state: **closed** · Severity: **moderate** · Accessible comments reviewed: **2**.

Workflow tags: `ldsc-seg`.

**Trigger and harm.** Baseline-plus-query analysis encounters the legacy empty-chromosome IndexError.
An invalid annotation prefix is hidden behind an unrelated indexing failure.

**Legacy disposition.** The author confirms the cause was the wrong annotation path and closes the issue; no code fix is identified.

**LDSC3 assessment — Apparently prevented.** LDSC3 resolves annotation path tokens before scoring and requires explicit baseline/query roles, avoiding comma-concatenated legacy-prefix ambiguity.
Current path failures identify zero matches.

**Recommended prevention.** Retain source-specific missing-path diagnostics and show quoted globs/exact files in a minimal baseline-plus-query recipe.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[src/ldsc/_runtime/paths.py:116](<../../../../src/ldsc/_runtime/paths.py#L116>); [src/ldsc/_annotation/preflight.py](<../../../../src/ldsc/_annotation/preflight.py>); [docs/reference/input-paths.md](<../../../reference/input-paths.md>).

Legacy implementation context: [`ldscore/regressions.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/regressions.py).

<a id="issue-401"></a>

### #401: how to add annotation to the baseline model???

[Original issue](https://github.com/bulik/ldsc/issues/401) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **0**.

Workflow tags: `ldsc-seg`.

**Trigger and harm.** The user cannot find how to add an annotation to the baseline.
They may replace the baseline or fit every query jointly when their aim is a separate conditional test.

**Legacy disposition.** Open; the short request has no answer or code fix.

**LDSC3 assessment — Apparently prevented.** LDSC3 explicitly separates baseline and query annotation inputs and fits baseline-plus-one-query models.
The relevant distinction is already implemented; this record primarily motivates discoverable examples.

**Recommended prevention.** Link a minimal baseline-plus-query command directly from the main tutorial and contrast it with a joint baseline-only fit.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[docs/user-guide/workflows/cell-type-enrichment.md](<../../../user-guide/workflows/cell-type-enrichment.md>); [docs/reference/results/partitioned-h2.md:11](<../../../reference/results/partitioned-h2.md#two-regimes>); [src/ldsc/_regression/runner.py:918](<../../../../src/ldsc/_regression/runner.py#L918>).

Legacy implementation context: [`ldscore/regressions.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/regressions.py).

<a id="issue-408"></a>

### #408: Error running LD score regression

[Original issue](https://github.com/bulik/ldsc/issues/408) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **0**.

Workflow tags: `h2`.

**Trigger and harm.** Legacy regression attempts to merge float64 and object identity columns.
The analysis aborts, and the generic pandas suggestion to concatenate would be scientifically unsafe.

**Legacy disposition.** Open; no source table, answer or verified fix is provided.

**LDSC3 assessment — Apparently prevented.** LDSC3 normalizes SNP identity centrally and validates supported identifier modes before joins.
Actual missing/numeric identifier content still needs examination; changing join type is not a valid repair.

**Recommended prevention.** Report source columns, normalized key examples and missing values on join failure; retain dtype-mixed identifier tests and never recommend pd.concat as an alignment substitute.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[src/ldsc/_core/columns.py](<../../../../src/ldsc/_core/columns.py>); [docs/reference/snp-identity-and-genome-build.md](<../../../reference/snp-identity-and-genome-build.md>); [tests/regression/test_legacy_sumstats_compatibility.py](<../../../../tests/regression/test_legacy_sumstats_compatibility.py>).

Legacy implementation context: [`ldscore/regressions.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/regressions.py).

<a id="issue-426"></a>

### #426: MemoryError: When running LDSCORE analysis

[Original issue](https://github.com/bulik/ldsc/issues/426) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **4**.

Workflow tags: `ldscore`.

**Trigger and harm.** Fourteen million variants are decoded with --chunk-size 500000 and cause MemoryError. A user-set batch can require enormous dense genotype buffers even on a large-memory cluster.

**Legacy disposition.** Comments recommend checking actual memory allocation and leaving the default chunk size; no code fix is identified.

**LDSC3 assessment — Apparently prevented.** LDSC3 uses a much smaller SNP batch default and exact windows independent of batch size.
Large user-requested batches and high concurrent worker counts can still exceed memory; a cluster label is not an allocation guarantee.

**Recommended prevention.** Explain batch-by-sample memory scaling and the worker multiplier; recommend a small batch and measured allocation before a larger redesign.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[docs/dev/design/memory-and-cpu.md](<../../design/memory-and-cpu.md>); [docs/reference/parallel-execution.md](<../../../reference/parallel-execution.md>); [tests/ldscore/test_plink_exact_windows.py](<../../../../tests/ldscore/test_plink_exact_windows.py>).

Legacy implementation context: [`ldscore/ldscore.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/ldscore.py).

<a id="issue-461"></a>

### #461: How to generate one own cts annotation file?

[Original issue](https://github.com/bulik/ldsc/issues/461) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **0**.

Workflow tags: `ldscore`, `ldsc-seg`.

**Trigger and harm.** Requests how to create query annotations for a tissue absent from the precomputed set.
Without a recipe users assume historical annotation archives are mandatory.

**Legacy disposition.** Open; no answer or code fix.

**LDSC3 assessment — Apparently prevented.** LDSC3 has BED and gene-list projection, explicit gene coordinate catalog/padding, baseline and optional control gene lists, followed by one scoring and partitioned-fit workflow.
Biological query construction still belongs to the researcher.

**Recommended prevention.** Make the custom-data workflow easy to find and show a minimal gene-list and ATAC-BED example with the appropriate control rationale.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[src/ldsc/_annotation/annotate.py](<../../../../src/ldsc/_annotation/annotate.py>); [docs/user-guide/workflows/cell-type-enrichment.md](<../../../user-guide/workflows/cell-type-enrichment.md>); [docs/reference/formats/gene-lists.md](<../../../reference/formats/gene-lists.md>).

Legacy implementation context: [`make_annot.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/make_annot.py), [`ldscore/regressions.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/regressions.py).

<a id="issue-466"></a>

### #466: Index error/ list index out of range

[Original issue](https://github.com/bulik/ldsc/issues/466) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **4**.

Workflow tags: `baseline-partitioned-h2`.

**Trigger and harm.** A chr3-specific prefix is passed to a legacy flag that appends all chromosome numbers. A malformed prefix produces empty discovery rather than the intended genome-wide fit.

**Legacy disposition.** A commenter removes the embedded chromosome and the reporter confirms success; later advice renames files.
No code fix is identified.

**LDSC3 assessment — Apparently prevented.** LDSC3 uses explicit path tokens for generation and one LD-score directory for regression.
It distinguishes per-chromosome computation from a complete fitted model; renaming files cannot fix their scientific contents.

**Recommended prevention.** Show chromosome-suite expansion and a genome-wide regression command side by side; preserve source-specific missing-input diagnostics.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[src/ldsc/_runtime/paths.py](<../../../../src/ldsc/_runtime/paths.py>); [docs/reference/input-paths.md](<../../../reference/input-paths.md>); [docs/user-guide/workflows/partitioned-heritability.md](<../../../user-guide/workflows/partitioned-heritability.md>).

Legacy implementation context: [`ldscore/regressions.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/regressions.py).

<a id="issue-477"></a>

### #477: munge_sumstats.py: error: argument --N: invalid float value: 'N'

[Original issue](https://github.com/bulik/ldsc/issues/477) · Source state: **closed** · Severity: **moderate** · Accessible comments reviewed: **1**.

Workflow tags: `munge-sumstats`.

**Trigger and harm.** Passes --N N and expects the literal N to name a numeric column.
Changing the file dtype cannot repair a command-line scalar parsing error.
The example also reverses effect and non-effect allele assignments.

**Legacy disposition.** A commenter correctly prescribes --N-col; no code fix is required.

**LDSC3 assessment — Apparently prevented.** LDSC3 keeps separate explicit scalar and column options with validation.
The interface still needs a novice-facing example because a column named N looks plausible after --N.

**Recommended prevention.** Show --N 100000 versus --N-col N side by side and map Effect_Allele to A1 in the same example, rather than fixing only the parse error.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[src/ldsc/_sumstats/munge_input.py:184](<../../../../src/ldsc/_sumstats/munge_input.py#L184>); [src/ldsc/_sumstats/munger.py](<../../../../src/ldsc/_sumstats/munger.py>); [docs/reference/cli.md](<../../../reference/cli.md>).

Legacy implementation context: [`munge_sumstats.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/munge_sumstats.py).

<a id="issue-482"></a>

### #482: ldsc: IndexError: list index out of range

[Original issue](https://github.com/bulik/ldsc/issues/482) · Source state: **closed** · Severity: **moderate** · Accessible comments reviewed: **2**.

Workflow tags: `rg`.

**Trigger and harm.** Ten-trait rg analysis fails at empty LD-score chromosome discovery because a directory path lacks its separator.
Users copy incompatible -chr/non-chr flag advice across different input layouts.

**Legacy disposition.** Author confirms adding the missing slash solved it; closure is an input correction, not a code fix.

**LDSC3 assessment — Apparently prevented.** LDSC3 regression consumes one explicit directory, removing legacy prefix concatenation, while its path validator names missing inputs.
All-pairs versus anchor default behavior also needs attention in this multi-trait use case.

**Recommended prevention.** Show one valid LD-score directory and explain all-pairs/anchor selection; preserve explicit zero-match errors rather than a generic IndexError.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[src/ldsc/_runtime/paths.py:116](<../../../../src/ldsc/_runtime/paths.py#L116>); [src/ldsc/_regression/runner.py:1233](<../../../../src/ldsc/_regression/runner.py#L1233>); [docs/reference/input-paths.md](<../../../reference/input-paths.md>).

Legacy implementation context: [`ldscore/regressions.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/regressions.py), [`ldscore/sumstats.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/sumstats.py).

<a id="faq-14"></a>

### FAQ-14: Regression execution time

[FAQ source](https://github.com/bulik/ldsc/wiki/FAQ) · Source-order identifier: **FAQ-14** · Severity: **low**.

Workflow tags: `h2`, `rg`, `baseline-partitioned-h2`, `ldsc-seg`.

**Legacy disposition.** Historical FAQ guidance; assessed independently rather than treated as a code-fix record.

**LDSC3 assessment — Partial coverage.** LDSC3 separates file preparation, fitting and result publication and can reuse prepared LD-score inputs.
Runtime depends on I/O and model count as well as regression size; historical timing cannot certify present throughput.

**Recommended prevention.** Explain input reuse and query-worker budgets, using measured performance evidence with its limits.

**Evidence and verification.** Source-inspected; relevant executed checks are identified in verification.md.
[docs/dev/validation/performance.md](<../performance.md>); [docs/reference/parallel-execution.md](<../../../reference/parallel-execution.md>).

<a id="issue-4"></a>

### #4: rename ldsc or ldsc.py

[Original issue](https://github.com/bulik/ldsc/issues/4) · Source state: **closed** · Severity: **low** · Accessible comments reviewed: **1**.

Workflow tags: Cross-cutting / not assigned to a requested group.

**Trigger and harm.** Importing a helper from ldsc.py resolves the sibling ldsc package instead.
Ambiguous import names complicate reuse and shell entry points.

**Legacy disposition.** Maintainers proposed renaming the internal package; final legacy LDSC uses ldscore/ alongside ldsc.py. A later referenced commit is not itself proof of the original fix.

**LDSC3 assessment — Apparently prevented.** LDSC3 uses a src-layout package, a setuptools console entry point, and python -m ldsc.
Supported objects are exported from the ldsc namespace; private domain modules are distinct from the command dispatcher.
There is no root ldsc.py wrapper to collide with the package.

**Recommended prevention.** Keep installation and Python examples on the installed namespace and console command; do not reintroduce root wrappers.

**Evidence and verification.** source-inspected
[src/ldsc/__init__.py](<../../../../src/ldsc/__init__.py>); [src/ldsc/__main__.py](<../../../../src/ldsc/__main__.py>); [setup.py](<../../../../setup.py>); [tests/package/test_package_layout.py](<../../../../tests/package/test_package_layout.py>).

<a id="issue-7"></a>

### #7: table at the end of --rg truncates long file names

[Original issue](https://github.com/bulik/ldsc/issues/7) · Source state: **closed** · Severity: **low** · Accessible comments reviewed: **1**.

Workflow tags: `rg`.

**Trigger and harm.** Long input paths are printed through a pandas display table.
Users can confuse traits whose visible prefixes match.

**Legacy disposition.** Maintainer recommended pd.set_option max_colwidth=1000; no exact production fixing commit is established.

**LDSC3 assessment — Apparently prevented.** LDSC3 writes complete trait names to rg.tsv and related tables instead of using console rendering as the result.
Duplicate trait names are disambiguated and their source paths are retained in diagnostics.
File contents therefore avoid pandas display truncation, although users still need meaningful trait labels when filenames are opaque.

**Recommended prevention.** Keep full names in saved tables and demonstrate naming two similarly named GWAS inputs.

**Evidence and verification.** source-inspected
[src/ldsc/_regression/runner.py:_disambiguate_trait_names](<../../../../src/ldsc/_regression/runner.py>); [src/ldsc/_regression/outputs.py:RgDirectoryWriter](<../../../../src/ldsc/_regression/outputs.py>); [tests/regression/test_regression_workflow.py:test_cli_does_not_print_rg_table_when_output_dir_is_supplied](<../../../../tests/regression/test_regression_workflow.py>).

Legacy implementation context: [`ldscore/regressions.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/regressions.py), [`ldscore/sumstats.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/sumstats.py).

<a id="issue-41"></a>

### #41: Incorrect command in tutorial? --out to wrong file?

[Original issue](https://github.com/bulik/ldsc/issues/41) · Source state: **closed** · Severity: **low** · Accessible comments reviewed: **1**.

Workflow tags: `h2`.

**Trigger and harm.** A tutorial writes one output prefix but tells users to inspect another and omits whitespace before a shell continuation.
Copy-and-paste examples fail or display stale results.

**Legacy disposition.** Maintainer confirms the tutorial was edited accordingly.

**LDSC3 assessment — Apparently prevented.** LDSC3 uses output directories with fixed result filenames and maintains worked CLI workflow pages.
This reduces output-prefix ambiguity, but example syntax still needs review as part of documentation validation.
CLI help smoke checks cannot prove every multiline tutorial is executable against real resources.

**Recommended prevention.** Check commands and subsequent inspect paths together; use one declared output root per worked example.

**Evidence and verification.** source-inspected
[docs/user-guide/guided-tutorial.md](<../../../user-guide/guided-tutorial.md>); [docs/user-guide/workflows/heritability.md](<../../../user-guide/workflows/heritability.md>); [src/ldsc/_regression/outputs.py:H2DirectoryWriter](<../../../../src/ldsc/_regression/outputs.py>).

Legacy implementation context: [`ldscore/regressions.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/ldscore/regressions.py).

<a id="issue-135"></a>

### #135: Include 'P-value'

[Original issue](https://github.com/bulik/ldsc/issues/135) · Source state: **open** · Severity: **low** · Accessible comments reviewed: **1**.

Workflow tags: `munge-sumstats`.

**Trigger and harm.** User thinks P-value is an unsupported header.
Unnecessary manual renaming can mask a delimiter or web-service parsing problem.

**Legacy disposition.** Maintainer confirms legacy support and asks for the actual error; no demonstrated alias bug.

**LDSC3 assessment — Apparently prevented.** LDSC3 clean_header replaces hyphens/dots with underscores, and the P aliases include P_VALUE.
Thus P-value is recognized by the native munger.
Header acceptance does not certify quoted CSV or an LD Hub upload, which use separate format/service rules.

**Recommended prevention.** Show accepted aliases via infer-only output and diagnose delimiter/header structure before suggesting a new synonym or contacting external-service support.

**Evidence and verification.** source-inspected
[src/ldsc/_core/columns.py:clean_header](<../../../../src/ldsc/_core/columns.py>); [src/ldsc/_core/columns.py:RAW_SUMSTATS_REQUIRED_OR_OPTIONAL_SPECS](<../../../../src/ldsc/_core/columns.py>).

Legacy implementation context: [`munge_sumstats.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/munge_sumstats.py).

<a id="issue-387"></a>

### #387: merge-alleles: command not found

[Original issue](https://github.com/bulik/ldsc/issues/387) · Source state: **open** · Severity: **low** · Accessible comments reviewed: **0**.

Workflow tags: `munge-sumstats`.

**Trigger and harm.** The script path is missing and --merge-alleles is then entered as a standalone shell command.
The user cannot start munging despite having an environment prompt.

**Legacy disposition.** Open; no answer or software fix in the thread.

**LDSC3 assessment — Apparently prevented.** LDSC3 installs ldsc as a console command and accepts options as part of one invocation.
Shell continuation mistakes can still happen regardless of package internals.

**Recommended prevention.** Provide a one-line munging command before multiline examples and explain that flags are arguments, not executable commands.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[setup.py](<../../../../setup.py>); [src/ldsc/cli.py](<../../../../src/ldsc/cli.py>); [README.md](<../../../../README.md>).

Legacy implementation context: [`munge_sumstats.py`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/munge_sumstats.py).

<a id="faq-02"></a>

### FAQ-02: Chromosome-level parallel execution

[FAQ source](https://github.com/bulik/ldsc/wiki/FAQ) · Source-order identifier: **FAQ-02** · Severity: **low**.

Workflow tags: `ldscore`.

**Legacy disposition.** Historical FAQ guidance; assessed independently rather than treated as a code-fix record.

**LDSC3 assessment — Apparently prevented.** Managed LD-score work is already divided by chromosome, with one native numerical thread per worker.
This is compatible with within-chromosome LD estimation; worker count is a resource budget rather than a change to the estimand.

**Recommended prevention.** Point users to --workers and distinguish chromosome tasks from BLAS threads.

**Evidence and verification.** Source-inspected; relevant executed checks are identified in verification.md.
[docs/reference/parallel-execution.md](<../../../reference/parallel-execution.md>); [src/ldsc/_runtime/parallelism.py](<../../../../src/ldsc/_runtime/parallelism.py>).

