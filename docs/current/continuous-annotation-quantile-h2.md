# Continuous-Annotation Quantile Heritability: Technical Contract

Last updated on: 2026-09-15

`ldsc quantile-h2` is a post-fit projection workflow for interpreting a continuous target annotation. It consumes exactly one fitted joint `partitioned-h2` model, reconstructs that model's common reference-SNP universe, assigns eligible SNPs to target-value quantiles, and projects the saved whole-data and delete-one-block coefficient vectors onto those quantiles. It does not refit LDSC.

## Model and interpretation

For fitted annotations `c`, the model is `Var(beta_j) = sum_c a_c(j) tau_c`. Annotation-specific LD scores, rather than the overlap matrix, identify `tau` and its jackknife uncertainty. The overlap-derived `prop_snps`, `prop_h2`, and `enrichment` values remain numerical output for legacy compatibility, but they do not have the usual category interpretation for a quantitative annotation. The raw coefficient, its SE, and its coefficient zero-test remain meaningful.

The quantile projection uses only fitted annotations:

`h2(Q_q) = sum_c tau_c sum_{j in Q_q} a_c(j)`.

The target annotation only defines `Q_q`. It may be one of the fitted annotations or an external score. An external target receives no coefficient and contributes no annotation sum to the projection.

The standardized coefficient is

`tau_star_c = tau_c SD(a_c) / (h2_g / M)`,

where `SD(a_c)` uses `ddof=0`, and both `h2_g` and `M` use the complete common reference-SNP universe before target-missing exclusions. The scale is fixed; it is not ratio-jackknifed. If the fitted-model total is nonpositive, only `tau_star` fields are `NaN`.

## Accepted fitted models

`--partitioned-h2-result-dir` accepts one complete fitted joint model:

- a baseline-only `partitioned-h2` output root; or
- one `diagnostics/query_annotations/<query>/` directory containing a baseline-plus-one-query fit.

An aggregate cell-type root is rejected because its rows came from separate regressions. Raw LDSC2 prefixes are not accepted. Current `partitioned-h2` writes float64 delete-one-block coefficients to `diagnostics/coefficient_delete_values.parquet` for a baseline-only model and to `coefficient_delete_values.parquet` inside each per-query directory.

## Required resupplied data

The command deliberately does not copy original annotation matrices into regression results. Users resupply:

- every original fitted annotation source with `--baseline-annot-sources` and, for a per-query model, the matching query source flag;
- the target source with `--target-annot-sources` and its unique column name with `--target-annotation`; and
- reference metadata with `--ref-metadata-sources`.

Path arguments accept exact files, globs, `@` chromosome placeholders, and existing path-token conventions. Parquet-R2 users can supply the existing `chr*_meta.tsv.gz` sidecars. PLINK users should supply metadata previously produced with `ldscore --export-ref-metadata`. Reference metadata requires `CHR`, `POS` (or `BP`), `SNP`, and `MAF`; allele-aware identity also uses `A1/A2`, with inference allowed only when the base identity is unique.

Annotation names must be globally unique within one LD-score artifact. The annotation classifier labels a column `binary` only when all finite values are exactly `0` or `1`; every other finite numeric column is `quantitative`. Classification affects logs and metadata only, never coefficients or other numerical results.

### Gene-list resolution

When reconstructing a gene-list query, supply `--query-annot-gene-list-sources` and reproduce the original catalog, padding, exclusions, optional control, and resolution policy. Resolution is strict by default. Add the value-free `--allow-unresolved-genes` flag if the original LD-score run deliberately used the `resolved-only` subset; do not append a policy value. The old `--gene-list-resolution-policy` option is rejected without an alias. Other source, coverage, and model-consistency checks still apply. See [resolution choices and migration](gene-list-input-format.md#strict-and-exploratory-resolution) and `add_quantile_h2_arguments()` in [quantile_h2.py](../../src/ldsc/quantile_h2.py).

## Missing target values

No target value is excluded by default. A nonnumeric, `NaN`, or infinite target value is an error that recommends `--target-missing-value`. When supplied, that one token is excluded before quantile construction. Finite numeric placeholders use numeric equality; nonnumeric and nonfinite spellings such as `NaN` use trimmed, case-sensitive string equality. Zero is an ordinary target value unless the user explicitly selects `--target-missing-value 0`.

## Reconstruction and provenance checks

The workflow inherits the linked LD-score artifact's inclusive common-MAF threshold; it has no common-MAF override. It validates effective SNP identity, duplicate rows, reference/annotation intersections, recorded all-reference and common reference-SNP universe sizes, recorded common-SNP fitted annotation sums, and available overlap cross-products. Aggregate comparisons retain `rtol=1e-6` and `atol=1e-8`. A target whose name matches a fitted annotation must match that resupplied fitted column's float32 values on target-eligible SNPs.

These checks establish alignment among resupplied sources and agreement with the available stored aggregate statistics. They do not establish the original annotation value at every SNP: value reassignments that preserve the checked sums and cross-products can pass and change the quantile projections. Resupply the original sources even when alternative sources pass validation. See `quantile_h2._prepare_quantile_inputs()` for the validation sequence.

Rows outside an otherwise valid intersection are recorded as exclusions. Duplicate identities, missing target coverage on the common universe, invalid MAF, and provenance mismatches are fatal. Row-addressable issues are written to `diagnostics/snp_alignment_issues.tsv.gz` even when the run fails.

## Quantiles and jackknife calculations

The default is five quantiles. Boundaries use the LDSC2 rule `floor(i * (n - 1) / Q + 0.5)` on sorted eligible target values. Intervals are low-to-high; ties at an internal boundary remain in the lower quantile. An empty realized quantile is a fatal error.

Within-quantile fitted-annotation sums are computed once. Matrix multiplication produces whole-data and all delete-one-block quantile totals. Ratio metrics are recomputed from every delete vector so denominator uncertainty is included. `enrichment_p` is the two-sided normal test of the inside-versus-complement per-SNP contrast; it is not `enrichment / enrichment_se` and does not test the target's coefficient.

When that contrast has zero or missing jackknife SE, `enrichment_p` is `NaN`. `compute_quantile_h2` skips the undefined division, including under a caller's strict NumPy error policy, while retaining the other quantile summaries.

## Output family

```text
<output_dir>/
  quantile_h2.tsv
  standardized_coefficients.tsv
  diagnostics/
    metadata.json
    quantile-h2.log
    snp_alignment_issues.tsv.gz
```

`quantile_h2.tsv` contains `quantile`, target bounds, `n_snps`, `prop_snps`, observed/liability `h2` and SE, `prop_h2` and SE, `enrichment` and SE, and `enrichment_p`. Liability fields are `NaN` unless the fitted regression recorded both prevalences. `standardized_coefficients.tsv` contains every fitted annotation's classification, SD, `tau`, SE, z, two-sided p-value, and corresponding `tau_star` fields.

`diagnostics/metadata.json` records the selected model, linked LD-score directory, target and reference sources, inherited common-MAF rule, common and eligible SNP counts, missing-token policy, quantile rule, and statistic definitions.

## Bounded exact reconstruction

The implementation validates the global common SNP universe using bounded annotation/reference/target reads and output-owned identity locators. A first pass accumulates full-common sums, cross-products, stable centered variance statistics, and eligible float64 targets. It computes the existing exact global quantile boundaries once, releases the target vector, and reuses those boundaries in a second chromosome pass. Quantile sums/counts require no whole-genome fitted-annotation matrix or dense SNP-by-quantile indicators. Ties, raw missing tokens, aggregate tolerances, and full-common standardization remain unchanged. Complete alignment diagnostics stream to a persistent file referenced by `QuantileH2Result.alignment_issues_path`. See [the memory design](annotation-memory-design.md#exact-global-quantiles).
