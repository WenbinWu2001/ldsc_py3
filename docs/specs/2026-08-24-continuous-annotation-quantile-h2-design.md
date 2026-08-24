# Continuous Annotation Support and Quantile Heritability Specification

Last updated on: 2026-08-24

## Status and references

Status: implemented and validated.

This specification captures the confirmed design for quantitative annotations in partitioned LDSC and the post-fit `ldsc quantile-h2` workflow. Scientific behavior is grounded in Gazal et al. (`docs/ldsc_papers/paper_continuous_annot.pdf` and `docs/ldsc_papers/supplementary_continous_annot.pdf`). Legacy operational behavior is grounded in `docs/wiki/partitioned-ldsc-continuous-annotations-ldsc2.md` and the original `ldsc_py2_Bulik/ContinuousAnnotations/quantile_M.pl` and `quantile_h2g.r` scripts.

## Problem and goal

Partitioned LDSC already fits quantitative annotations correctly through annotation-specific LD scores, but its category-size, proportion, and enrichment summaries use annotation sums as though they were SNP counts. Those summaries retain numerical values for compatibility but do not have the ordinary binary-category interpretation for quantitative rows. The current package also does not persist coefficient delete values or provide the legacy post-fit quantile projection and standardized coefficient summaries needed to interpret quantitative annotations.

The feature succeeds when:

- native LD-score artifacts classify and minimally fingerprint annotations without changing LD-score or regression numerics;
- partitioned-h2 retains its current numerical tables, logs quantitative-row interpretation warnings, and persists coefficient delete-one-block values for each fitted model;
- `ldsc quantile-h2` accepts one canonical LDSC3 fitted model plus resupplied annotation and reference metadata, verifies them, and reports joint-model heritability by target quantile and paper-compatible standardized coefficients;
- binary-only workflows remain numerically unchanged; and
- failures that could change the scientific result abort with actionable messages and machine-readable SNP diagnostics.

## Domain language

- **Binary annotation:** every retained annotation value is exactly `0` or `1`.
- **Quantitative annotation:** at least one retained value is not exactly `0` or `1`, including discrete nonbinary annotations.
- **Unknown annotation type:** a legacy-converted or older artifact whose original values were unavailable for classification.
- **Common reference-SNP universe:** reference SNPs that pass the inherited common-MAF threshold and operator.
- **Fitted model:** one complete joint coefficient vector. A baseline-only result is one fitted model; each cell-type per-query directory is a separate fitted baseline-plus-one-query model.
- **Target annotation:** the annotation used only to assign SNPs to quantiles. It may or may not be part of the fitted model.

## Observable behavior and acceptance criteria

### Native LD-score construction

1. Each retained annotation is classified once over retained reference-SNP values and recorded as `binary` or `quantitative`.
2. Classification is advisory and does not alter LD scores, overlap values, counts, coefficients, uncertainty, P values, or result rows.
3. Annotation names are globally unique across baseline and query groups. Duplicate names fail before computation and when a persisted artifact is loaded.
4. Missing, nonnumeric, NaN, positive-infinite, and negative-infinite fitted annotation values are rejected.
5. Native artifacts record the approved common-reference-SNP-universe and per-annotation fingerprints. No SNP-by-annotation matrix is persisted.

### Partitioned-h2

1. Binary-only results remain numerically identical to current behavior.
2. Quantitative rows retain numerical `prop_snps`, `prop_h2`, `enrichment`, uncertainty, and overlap-derived summaries for legacy compatibility.
3. The log lists quantitative annotations and warns that those weighted summaries do not have the ordinary binary-category interpretation. It identifies `tau`, its SE, and its zero-test as interpretable and directs users to `ldsc quantile-h2`.
4. The warning is emitted once per run. It separately explains that binary quantile annotations followed by refitting are required to estimate distinct conditional bin coefficients.
5. Every written fitted model has a coefficient delete-value artifact in exact fitted coefficient order.

### Fitted-model selection

1. A baseline-only partitioned-h2 result root is accepted as one fitted model.
2. One `diagnostics/query_annotations/<query>/` directory is accepted as one fitted baseline-plus-query model.
3. An aggregate multi-query root is rejected because its rows come from separate regressions. Coefficients are never concatenated across query runs.
4. A result without `coefficient_delete_values.parquet` is rejected with instructions to rerun `partitioned-h2`; LD scores do not need to be recomputed.

### Target, missingness, and quantiles

1. The target is supplied explicitly even when it is fitted. An external target defines quantile membership but contributes no coefficient, annotation sum, or standardized coefficient to the fitted-model projection.
2. The target must contain every SNP in the common reference-SNP universe. Omitted rows are errors rather than implicit missingness.
3. `--target-missing-value` is the sole exclusion mechanism and defaults to no exclusion. Finite numeric tokens use numeric equality; nonnumeric tokens use trimmed, exact, case-sensitive raw-field matching before numeric conversion.
4. After explicit exclusions, every target value must be numeric and finite. Zero is valid unless explicitly selected as the missing token.
5. `--num-quantiles` defaults to 5 and must be at least 2.
6. Quantiles follow the legacy algorithm: choose boundaries at `round(i * (n - 1) / Q)`, assign the first interval as `[minimum, first upper boundary]`, assign later intervals as `(lower, upper]`, and keep boundary ties in the lower-valued quantile.
7. Empty quantiles are fatal and never silently merged. Rows are ordered low to high and report actual boundaries and SNP counts.

### Failures and degenerate estimates

1. Duplicate effective identities, missing required target rows, missing or nonfinite intersected MAF, failed allele inference, fitted-source coverage defects, and provenance mismatches are fatal.
2. Reference-metadata-only and annotation-only SNPs are nonfatal exclusions when reconstructed aggregate contracts still match the fitted artifact.
3. Negative quantile heritability estimates remain visible. A negative but nonzero target-eligible total retains proportion and enrichment calculations with a prominent warning.
4. A zero target-eligible total retains quantile heritability but writes `NaN` for ratio-derived fields.
5. A nonpositive fitted-model total over the complete common reference-SNP universe retains raw `tau` fields but writes `NaN` for all `tau_star` fields.
6. If any delete replicate required by a metric has an invalid denominator, that metric's SE and P value are `NaN`; delete replicates are never selectively dropped.

## Scientific contracts

### Regression and overlap

Annotation-specific LD scores retain the quantitative definition

$$
l(j,c)=\sum_k a_c(k)r_{jk}^2.
$$

The overlap matrix is not used to estimate `tau` or its covariance. It remains useful after the joint fit. For a binary set (B), including when other fitted annotations are quantitative,

$$
\widehat h_g^2(B)=\sum_c \widehat\tau_c\sum_{j\in B}a_c(j).
$$

For a quantitative row the same weighted cross-products are defined but do not describe a literal SNP category. The implementation warns without suppressing or changing those numerical outputs.

### Quantile projection and jackknife

Let (Q_q) be target quantile (q), let (C) be the retained fitted annotations, and define

$$
S_{cq}=\sum_{j\in Q_q}a_c(j), \qquad \widehat h_g^2(Q_q)=\sum_{c\in C}\widehat\tau_cS_{cq}.
$$

If (S) has shape (C\times Q), `tau` has length (C), and `T_delete` has shape (B\times C), compute

$$
\mathbf h=S^\mathsf{T}\widehat{\boldsymbol\tau}, \qquad H_{delete}=T_{delete}S.
$$

The sufficient-statistic matrix (S) is accumulated once. Projection and uncertainty calculations use vectorized matrix operations.

For fixed SNP counts (N_q),

$$
\mathrm{prop\_h2}_q=\frac{h_q}{\sum_rh_r}, \qquad \mathrm{enrichment}_q=\frac{\mathrm{prop\_h2}_q}{N_q/\sum_rN_r}.
$$

`enrichment_p` is the two-sided normal P value for

$$
\frac{h_q}{N_q}-\frac{h_{all}-h_q}{N_{all}-N_q},
$$

recomputed from every delete coefficient vector. It is not `enrichment / enrichment_se` and does not test the target coefficient. Block-jackknife SEs use the package's existing pseudovalue convention and the actual number of persisted blocks. `enrichment_se` equals `prop_h2_se / prop_snps` because SNP counts are fixed.

### Standardized coefficients

For fitted annotation (c), let (M) be the number of SNPs in the complete common reference-SNP universe and use population SD (`ddof=0`) over that universe:

$$
\tau_c^*=\frac{\widehat\tau_c\,\mathrm{SD}(a_c)}{\widehat h_g^2/M}, \qquad \mathrm{SE}(\tau_c^*)=\frac{\mathrm{SE}(\widehat\tau_c)\,\mathrm{SD}(a_c)}{\widehat h_g^2/M}.
$$

The denominator is not ratio-jackknifed. `tau_p` and `tau_star_p` are the same two-sided normal zero-test. Standardized coefficients use the complete common reference-SNP universe and do not change when target missingness excludes SNPs from quantile summaries.

Existing observed-to-liability conversion applies to absolute quantile heritability and its SE. Proportions, enrichment, P values, `tau`, and `tau_star` are scale invariant.

## Public interface

The command is `ldsc quantile-h2`. `--output-dir` is required because the command writes multiple headline artifacts and diagnostics.

```text
--partitioned-h2-result-dir
--baseline-annot-sources
--query-annot-sources
--query-annot-bed-sources
--query-annot-gene-list-sources
--target-annot-sources
--target-annotation
--ref-metadata-sources
--target-missing-value
--num-quantiles
--output-dir
--overwrite
--log-level
```

The three query-source forms are mutually exclusive and retain existing BED/gene projection options. Whole-genome inputs and explicit `@` chromosome suites use the existing path-token language. Fitted sources may contain extra columns; the fitted result selects and reorders exactly its retained annotation names. If the target name matches a fitted annotation, the vectors must match or the external target must be renamed.

`--ref-metadata-sources` accepts whole-genome or `@` chromosome-suite tables containing at least `CHR`, `POS`, `SNP`, and `MAF`, plus alleles when required. Parquet-R2 users may supply existing `chr*_meta.tsv.gz` sidecars. PLINK users may supply metadata produced by `ldscore --export-ref-metadata`. The command inherits the fitted threshold and operator, exposes no common-MAF override, and never recomputes MAF from genotypes.

## Identity, alignment, and provenance

The first baseline source defines the annotation SNP grid. Other prebuilt baseline and query sources cover exactly that grid by effective identity, although row order may differ. BED and gene-list queries are projected onto it. The fitted-model reference universe is the intersection of this grid with `--ref-metadata-sources`, matching LD-score construction. Reference metadata order supplies deterministic processing order. In allele-aware modes, omitted annotation alleles may be inferred only from unambiguous reference metadata.

Minimal aggregate verification requires exact identity-mode and genome-build compatibility, globally unique annotation names, one occurrence of every retained fitted annotation, reconstructed all-reference and common reference-SNP universe sizes matching `overlap_config`, annotation sums matching stored count records, and available fitted-model cross-products matching the overlap artifact. Aggregates use deterministic float64 accumulation and `numpy.isclose(rtol=1e-6, atol=1e-8)`; mismatches are fatal and name the annotation and statistic.

Native LD-score metadata records:

```json
{
  "annotation_fingerprints": {
    "algorithm": "sha256",
    "canonicalization": "ldsc_common_annotation_v1",
    "common_reference_snp_universe": "<sha256>",
    "annotation_values": {"base": "<sha256>"}
  }
}
```

The universe fingerprint hashes sorted effective identities in the common reference-SNP universe. Each annotation fingerprint hashes a canonical two-column frame of `effective_snp_id` and normalized float32 annotation value over the same universe. The frame is sorted by identity, formats values with 17 significant digits, and serializes as headered UTF-8 TSV with `\n` line endings. Paths, compression, source-file bytes, MAF changes that do not change membership, rare-SNP values, and unrelated columns are excluded.

Native matches record `verification_level=exact_common_values`. Older LDSC3 and converted LDSC2 LD-score artifacts without fingerprints remain usable with aggregate checks, `verification_level=aggregate_only`, and a warning. Raw LDSC2 regression result prefixes are not accepted.

## Artifact contracts

### Partitioned-h2 additions

Baseline-only results write `diagnostics/coefficient_delete_values.parquet`; per-query results write `diagnostics/query_annotations/<query>/coefficient_delete_values.parquet`. The wide parquet contains zero-based integer `delete_block` followed by float64 retained annotation columns in coefficient order. Query artifacts are written only with the per-query tree. Containing metadata exposes `files.coefficient_delete_values`, block count, and retained order.

### Quantile-h2 result family

```text
<output-dir>/
  quantile_h2.tsv
  standardized_coefficients.tsv
  diagnostics/
    metadata.json
    quantile-h2.log
    snp_alignment_issues.tsv.gz
```

`quantile_h2.tsv` columns are:

```text
quantile target_value_lower target_value_upper n_snps prop_snps
h2_obs h2_obs_se h2_liab h2_liab_se
prop_h2 prop_h2_se enrichment enrichment_se enrichment_p
```

`quantile` is an integer from 1 through (Q). Liability columns remain present and contain `NaN` when prevalence is unavailable.

`standardized_coefficients.tsv` columns are:

```text
annotation annotation_type annotation_sd
tau tau_se tau_z tau_p
tau_star tau_star_se tau_star_z tau_star_p
```

`annotation_type` is descriptive and never changes numerical fields.

Successful `diagnostics/metadata.json` uses `artifact_type=quantile_h2_result` and records relative files, selected model path/type, linked LD-score directory, identity/build, retained order, source and target provenance, common-MAF rule, SNP counts, quantile and boundary rules, missing-token rule/count, prevalence, P-value conventions, fingerprint canonicalization, and verification level.

`diagnostics/snp_alignment_issues.tsv.gz` is always written, including header-only on a clean run, with columns:

```text
source_role source CHR POS SNP A1 A2 effective_snp_id issue action details
```

Core issues are `duplicate_identity`, `missing_reference_metadata`, `missing_baseline_annotation`, `missing_query_annotation`, `missing_target_annotation`, `missing_maf`, and `allele_inference_failed`; `action` is `excluded` or `fatal`. Fatal validation retains this file and the workflow log but writes no headline TSVs or successful metadata.

## Logging contract

`ldscore` logs one compact classification block. `partitioned-h2` repeats persisted classification and warns once per run. `quantile-h2` logs the fitted model, target, inherited common-MAF rule, common reference-SNP universe size, missing exclusion count, quantile boundaries/counts, provenance checks, and numerical warnings.

## Validation strategy

Unit seams cover annotation validation and classification, global-name uniqueness, fingerprints, legacy quantile boundaries and ties, missing-token behavior, matrix projection, jackknife statistics, `tau_star`, alignment issues, aggregate tolerances, and verification levels.

Integration seams cover LD-score metadata round trips, baseline-only and per-query delete artifacts, fitted and external targets, whole-genome and `@`-suite inputs, numeric and string missing tokens, observed and liability outputs, fatal diagnostics, and regression equivalence. A deterministic legacy fixture compares the seven legacy quantile statistics with the corresponding new columns.

Milestones run focused pytest modules, then full `pytest`, standard-library unittest compatibility, `ldsc --help`, `ldsc quantile-h2 --help`, and real baseline-only and per-query smoke runs.

## Scientific and operational constraints

- This is a post-fit projection, not a nonlinear test or refit of quantile bins.
- Inputs must already have appropriate scaling, transformation, missingness, and confounder adjustment. The command does not normalize annotations or add MAF bins.
- Common reference-SNP universe, build, identity, alleles, prevalence, and coefficient order are immutable fitted-model provenance.
- Inputs are read-only; outputs follow standard preflight and overwrite policy.
- Full annotation matrices are not persisted. SNP scans build common reference-SNP universe columns and sufficient statistics once; uncertainty uses matrix operations.
- Existing dependencies are sufficient.
- Implementation must update workflow, CLI, artifact, troubleshooting, interpretation, and example documentation, explicitly covering Parquet-R2 sidecars and PLINK `--export-ref-metadata`.

## Out of scope

- Changing the regression estimator, coefficient scaling, overlap formulas, or existing partitioned-h2 numerical rows.
- Constructing binary quantile annotations or refitting quantile-bin models.
- Accepting raw LDSC2 `.results` and `.part_delete` prefixes.
- Persisting full annotation matrices or recomputing MAF from PLINK genotypes.
- Adding a common-MAF override, automatic annotation normalization, plotting, or meta-analysis.

## Risks and open questions

No product, scientific, or architectural questions remain open. Implementation risks are coefficient/annotation order drift, fitted-universe reconstruction errors, aggregate tolerance behavior, failure-diagnostic atomicity, and output-schema drift. Any discovery requiring a contract change returns to design review.
