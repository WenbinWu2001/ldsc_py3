# Partitioned LDSC Result Columns

This page describes the tabular outputs from `ldsc partitioned-h2`. The summaries
are **overlap-aware**: they reproduce the legacy `--overlap-annot` math, in which
each category's heritability is propagated through the annotation overlap matrix
`O = AᵀA` rather than assuming disjoint categories. The required overlap matrix
is read from `ldscore.overlap.parquet` in the LD-score directory; a directory
produced by an older `ldsc ldscore` is rejected with a regenerate message
(`h2` and `rg` do not need it).

## Two regimes

`partitioned-h2` auto-detects one of two regimes from the LD-score directory,
with no user flag:

- **Functional-category regime** — the directory has **no** query columns. A
  single joint fit of all baseline annotations yields one row per baseline
  category. The headline is `Enrichment` (with two-sided `Enrichment_p`). This
  reproduces the joint functional-enrichment analysis (Finucane et al. 2015 /
  legacy `--overlap-annot`).
- **Cell-type-specific regime** — the directory **has** query columns. For each
  query annotation a `baseline + one query` model is fit, yielding one row per
  query. The headline is `Coefficient` (with one-sided `Coefficient_p`). This is
  the cell-type-specific analysis (Finucane et al. 2018): the query's conditional
  contribution beyond the baseline.

Both regimes write the **same** `partitioned_h2.tsv` column schema; they differ
only in which rows appear (baseline categories vs. queries) and the default sort.

## `partitioned_h2.tsv`

One uniform schema, written in this column order:

| Column | Meaning |
| --- | --- |
| `Category` | Baseline category (functional regime) or query annotation (cell-type regime). |
| `Prop._SNPs` | `M_c / M_tot`: the category's share of the reference SNP universe. |
| `Category_h2_obs` | **Observed-scale conditional** category heritability `M_c · tau_c` (this category's own coefficient times its SNPs). Can be negative under overlap. |
| `Category_h2_obs_std_error` | Jackknife standard error for `Category_h2_obs`. |
| `Category_h2_liab` | **Liability-scale** conditional category heritability: `Category_h2_obs` times the observed-to-liability factor `c(samp_prev, pop_prev)`. `NaN` unless both `--samp-prev` and `--pop-prev` are supplied. |
| `Category_h2_liab_std_error` | Jackknife standard error for `Category_h2_liab` (the observed SE scaled by the same `c`). |
| `samp_prev` | Sample (case) prevalence applied for liability conversion, or `NaN` for an observed-scale run. |
| `pop_prev` | Population prevalence applied for liability conversion, or `NaN`. |
| `Prop._h2` | **Marginal**, overlap-aware proportion of total heritability explained by the SNPs in this category: `Σⱼ (O[c,j]/Mⱼ)·propⱼ`. Folds in every category those SNPs also belong to. |
| `Prop._h2_std_error` | Jackknife standard error for `Prop._h2`. |
| `Enrichment` | `Prop._h2 / Prop._SNPs`. The headline of the functional regime. |
| `Enrichment_std_error` | Standard error of `Enrichment`. |
| `Enrichment_p` | **Two-sided** t-test (df = `n_blocks`) of whether the category's mean per-SNP heritability differs from the rest of the model (the legacy overlap contrast). `NaN` for a category that contains every SNP (e.g. `base`). |
| `Coefficient` | Conditional per-SNP contribution `tau_c` (controls for the other model annotations). The headline of the cell-type regime. |
| `Coefficient_std_error` | Jackknife standard error for `Coefficient`. |
| `Coefficient_z` | `Coefficient / Coefficient_std_error`. |
| `Coefficient_p` | **One-sided** normal p-value testing `Coefficient > 0` (does this annotation add heritability beyond the baseline?). |
| `overlap_aware` | Whether the fitted model's annotations overlap (any off-diagonal of `O` > 0). With a `base` (all-ones) category this is effectively always true. |

The same schema is used for each one-row
`diagnostics/query_annotations/<folder>/partitioned_h2.tsv` file.

### Choosing the right column

- **Functional enrichment** ("are this category's SNPs enriched for
  heritability?"): read `Enrichment` and `Enrichment_p`. `Enrichment > 1` means
  the category's SNPs carry a larger heritability share than their SNP share; a
  small `Enrichment_p` means the enrichment is significantly different from 1.
- **Cell-type-specific contribution** ("does this annotation add heritability
  beyond the baseline?"): read `Coefficient` together with the one-sided
  `Coefficient_p`. A positive `Coefficient` with a small `Coefficient_p` is the
  de-confounded signal. The `Enrichment` of a query is confounded by its overlap
  with generally-enriched baseline categories, so prefer `Coefficient_p` here.

### Sorting and signaling

`--summary-sort-by` defaults to `auto`, which resolves to `coefficient-p` in the
cell-type regime (most significant query first) and `category` in the functional
regime (baseline-annotation order). Explicit keys are `category`, `prop-snps`,
`prop-h2`, `enrichment`, `enrichment-p`, `coefficient`, and `coefficient-p`;
p-value keys sort ascending, other numeric keys descending, ties preserve input
order, and missing values are placed last. When per-query diagnostics are
enabled, the manifest order and folder ordinals follow the sorted summary order.

The run logs a one-block regime banner naming the column to focus on, and records
`analysis_type` (`functional_category` | `cell_type_specific`), `headline_metric`
(`enrichment` | `coefficient`), `enrichment_p_test` (`two_sided_t`), and
`coefficient_p_test` (`one_sided_greater`) in `diagnostics/metadata.json`.

Output directories follow the coherent artifact-family policy. The root
`partitioned_h2.tsv`, diagnostic metadata, optional
`diagnostics/query_annotations/` tree, and `diagnostics/partitioned-h2.log` are
checked together. Without overwrite, any existing owned sibling rejects the run.
With overwrite, a successful aggregate-only run removes a stale
`diagnostics/query_annotations/` tree from an earlier per-query run. If
`--output-dir` is omitted, the CLI prints the table to stdout and writes no
diagnostics.

## `partitioned_h2_full.tsv`

When `--write-per-query-results` is supplied (cell-type regime), each query
folder also contains `partitioned_h2_full.tsv`. This table uses the **same**
column schema as `partitioned_h2.tsv`, with one row for every retained category
in that query's fitted `baseline + query` model (all baseline categories plus the
query). Its baseline rows are the joint functional enrichments conditional on the
added query.

## Notes For Interpretation

`Prop._h2`, `Category_h2_obs` (and its liability counterpart `Category_h2_liab`),
and `Enrichment` are regression estimates and can be negative or greater than one
when estimates are noisy or annotations are correlated; `Category_h2_obs` is
especially prone to this under overlap because it is the conditional `M_c · tau_c`.
Only the `Category_h2_*` columns differ between scales; `Prop._SNPs`, `Prop._h2`,
`Enrichment`, and the coefficients are scale-invariant. P-values and standard errors may be missing when a
valid comparison cannot be formed — for example `Enrichment_p` is `NaN` for a
category that contains every SNP, and `Coefficient_p` is `NaN` when the
coefficient standard error is zero.

The `diagnostics/query_annotations/manifest.tsv` file is an index for the
per-query result tree. It records each original query annotation name, the
sanitized folder name, and the relative path columns `summary_path`,
`partitioned_h2_full_path`, and `metadata_path`. It is not a scientific result
table.
