# Partitioned LDSC Result Columns

This page describes the tabular outputs from `ldsc partitioned-h2`. The root
summary table is designed for scanning query annotation results. The optional
per-query folders add the full category-level table for each fitted model.

## `partitioned_h2.tsv`

The root `partitioned_h2.tsv` has one row per query annotation. Columns are
written in this order. The file is generated only for LD-score directories that
contain query annotations; baseline-only LD-score directories should be analyzed
with `ldsc h2` or `ldsc rg`, not `ldsc partitioned-h2`.

| Column | Meaning |
| --- | --- |
| `Category` | Query annotation tested in the baseline-plus-query model. |
| `Prop._SNPs` | Proportion of reference SNPs assigned to the category. |
| `Prop._h2` | Proportion of total SNP heritability attributed to the category. |
| `Enrichment` | Heritability enrichment, computed as `Prop._h2 / Prop._SNPs`. |
| `Enrichment_p` | Two-sided p-value testing whether per-SNP heritability in the category differs from the remaining model categories. |
| `Coefficient` | LDSC regression coefficient for the category. |
| `Coefficient_p` | Two-sided p-value testing whether the coefficient differs from zero. |

The same seven-column schema is used for each one-row
`query_annotations/<folder>/partitioned_h2.tsv` file.

## `partitioned_h2_full.tsv`

When `--write-per-query-results` is supplied, each query folder also contains
`partitioned_h2_full.tsv`. This table has one row for every retained category in
the fitted baseline-plus-query model, including baseline categories and the
query category. Columns are written in this order:

| Column | Meaning |
| --- | --- |
| `Category` | Model category represented by the row. |
| `Prop._SNPs` | Proportion of reference SNPs assigned to the category. |
| `Category_h2` | Estimated SNP heritability contribution for the category. |
| `Category_h2_std_error` | Jackknife standard error for `Category_h2`. |
| `Prop._h2` | Proportion of total SNP heritability attributed to the category. |
| `Prop._h2_std_error` | Jackknife standard error for `Prop._h2`. |
| `Enrichment` | Heritability enrichment, computed as `Prop._h2 / Prop._SNPs`. |
| `Enrichment_std_error` | Jackknife standard error for `Enrichment`. |
| `Enrichment_p` | Two-sided p-value testing whether per-SNP heritability in the category differs from the remaining model categories. |
| `Coefficient` | LDSC regression coefficient for the category. |
| `Coefficient_std_error` | Jackknife standard error for `Coefficient`. |
| `Coefficient_p` | Two-sided p-value testing whether the coefficient differs from zero. |

## Notes For Interpretation

`Prop._h2`, `Category_h2`, and `Enrichment` are regression estimates and can be
negative or greater than one. This can happen when estimates are noisy or when
annotations are correlated. P-values and standard errors may be missing when a
valid comparison cannot be formed, for example when there is no non-category
complement.

The `query_annotations/manifest.tsv` file is an index for the per-query result
tree. It records each original query annotation name, the sanitized folder name,
and the relative path columns `summary_path`, `partitioned_h2_full_path`, and
`metadata_path`. It is not a scientific result table.
