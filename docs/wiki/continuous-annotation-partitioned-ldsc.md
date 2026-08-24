# Continuous annotations in partitioned LDSC

Last updated on: 2026-08-24

This page starts after the partitioned-LDSC steps in the [guided tutorial](guided-tutorial.md). The regression itself does not use quantile bins: fit the continuous annotation directly, then use `ldsc quantile-h2` to summarize how the complete fitted joint model distributes heritability from low to high target values.

## What to interpret

For a quantitative annotation, interpret its conditional coefficient `tau`, SE, and coefficient p-value. The numerical `prop_snps`, `prop_h2`, and `enrichment` fields remain in `partitioned_h2.tsv` for compatibility, but they are weighted summaries and do not have the ordinary binary-category meaning.

The post-fit command adds two complementary views:

- `quantile_h2.tsv` reports joint-model heritability in target-defined quantiles.
- `standardized_coefficients.tsv` reports `tau_star`, the coefficient per one annotation SD divided by average per-SNP heritability.

The target score may be fitted or external. It only assigns SNPs to quantiles; heritability is always calculated from the annotations and coefficients in the selected fitted model.

## 1. Preserve one complete fitted model

For a baseline-only partitioned model, use the `partitioned-h2` output root. For a cell-type-specific run, include `--write-per-query-results` when fitting and later select one query directory:

```text
partitioned-h2/<trait>/diagnostics/query_annotations/0001_<query>/
```

Do not pass the aggregate multi-query root: its rows belong to different baseline-plus-one-query regressions.

## 2. Prepare the resupplied sources

You need every annotation source used to fit the selected model, one target source, and reference metadata. LDSC3 stores compact fingerprints and regression delete values, not copies of the large annotation matrices.

- Parquet-R2 users can supply the existing `chr*_meta.tsv.gz` sidecars.
- PLINK users should first run `ldsc ldscore --export-ref-metadata`; use the resulting `ref_metadata/chr@_meta.tsv.gz` files.

Annotation names must be unique across all baseline and query sources in one LD-score artifact.

## 3. Run the post-fit projection

Baseline-only example:

```bash
ldsc quantile-h2 \
  --partitioned-h2-result-dir "${PARTITIONED_H2_DIR}" \
  --baseline-annot-sources "${BASELINE_ANNOT_DIR}/baseline.@.annot.gz" \
  --target-annot-sources "${TARGET_ANNOT_DIR}/target.@.annot.gz" \
  --target-annotation "continuous_score" \
  --ref-metadata-sources "${R2_DIR}/chr@_meta.tsv.gz" \
  --num-quantiles 5 \
  --output-dir "${OUTPUT_ROOT}/continuous-score-quantiles" \
  --overwrite
```

For one cell-type-specific model, point `--partitioned-h2-result-dir` to its per-query directory and resupply that fitted query source as well:

```bash
ldsc quantile-h2 \
  --partitioned-h2-result-dir "${PARTITIONED_H2_DIR}/diagnostics/query_annotations/0001_my_query" \
  --baseline-annot-sources "${BASELINE_ANNOT_DIR}/baseline.@.annot.gz" \
  --query-annot-sources "${QUERY_ANNOT_DIR}/my_query.@.annot.gz" \
  --target-annot-sources "${TARGET_ANNOT_DIR}/external_score.@.annot.gz" \
  --target-annotation "external_score" \
  --ref-metadata-sources "${R2_DIR}/chr@_meta.tsv.gz" \
  --output-dir "${OUTPUT_ROOT}/external-score-quantiles"
```

Use the same query input form used for LD scores: `--query-annot-sources`, `--query-annot-bed-sources`, or `--query-annot-gene-list-sources` with its gene options.

By default, every target value must be finite and numeric. If exactly one token means missingness, name it explicitly, for example `--target-missing-value NaN`. Those SNPs are excluded from quantile construction. Zero is retained unless you explicitly choose zero as the missing token.

## 4. Read the outputs

Rows in `quantile_h2.tsv` run from the lowest to the highest target values.

| Column group | Interpretation |
| --- | --- |
| `quantile`, `target_value_lower`, `target_value_upper`, `n_snps`, `prop_snps` | Realized target interval and eligible-SNP share. |
| `h2_obs`, `h2_obs_se` | Joint fitted-model heritability total and block-jackknife SE in the quantile. |
| `h2_liab`, `h2_liab_se` | Liability-scale counterparts when the fitted model has both prevalences. |
| `prop_h2`, `prop_h2_se` | Quantile heritability divided by the total across quantiles, with denominator uncertainty. |
| `enrichment`, `enrichment_se` | Heritability share divided by eligible-SNP share. |
| `enrichment_p` | Two-sided inside-versus-complement test of average per-SNP heritability. |

`standardized_coefficients.tsv` has one row per fitted annotation. `tau_star` uses the complete common reference-SNP universe, including SNPs excluded only because the target value was missing. If total fitted-model heritability on that universe is nonpositive, raw `tau` remains available and `tau_star` is reported as `NaN`.

The diagnostics directory contains the run log, full provenance, and a compressed table of excluded or problematic SNP identities. If you need coefficients for separate binary quantile categories rather than a post-fit projection, construct binary bins, recompute their LD scores, and fit those bins as annotations.
