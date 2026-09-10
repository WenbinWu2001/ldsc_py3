# LD-score SNP-universe contract

Last updated on: 2026-09-10

This table distinguishes the SNP universes used by canonical LD Score
Regression (LDSC) and stratified LDSC. It is especially important for MHC and
pericentromeric exclusions.

| Quantity | Appropriate SNP universe |
|---|---|
| Baseline/query LD-score contributors \(k\) | Broad retained 1000 Genomes/PLINK reference universe, including MHC and pericentromeric SNPs |
| Written/regression rows \(j\) | HapMap3 intersected with the complement of MHC, pericentromeric, and other declared regression-row exclusions |
| Regression-weight LD score `w_ld` | Computed within the filtered HapMap3 regression-SNP set; excluded-region SNPs are neither rows nor `w_ld` contributors |
| `M`, `M_5_50`, and baseline-overlap counts | Broad retained reference/causal universe (with the ordinary MAF threshold for `M_5_50`), not only regression rows |
| Final `h2`, `rg`, and `partitioned-h2` observations | Intersection of the filtered regression rows, written LD-score rows, regression weights, and munged summary statistics |

## Chromosome coverage versus SNP support

Direct query chromosome coverage is established from validated baseline/reference inputs before SNP, sample, or MAF filtering. Their chromosome sets must match exactly. `@` declares autosomes 1–22; exact paths and globs select actual artifacts whose contents establish scope. Every selected gene must lie within scope after explicit gene exclusions. Coverage failures abort the batch, whereas a valid chromosome whose retained computational SNP universe becomes empty supplies a measured zero-support outcome. Direct gene support uses the prepared baseline/reference intersection; indexed support uses immutable gene-to-atom counts. Unknown or unevaluated support stays blank. Query LD-score variance is assessed separately on regression rows. See [coverage diagnostics](gene-list-diagnostics-and-repair.md#chromosome-scope-and-pathway-coverage).

## Projection and traversal contract

For each chromosome, baseline/control annotations, query annotations, and the binary regression-SNP mask share one numerical LD traversal. PLINK computes each genotype-correlation block once; parquet-R2 decodes each stored pair chunk once. Baseline/control and regression-weight projections are computed once per block/chunk, and every query batch is projected against that same block/chunk before it is released. Query batching must not restart the chromosome kernel or repeat genotype/R2 traversal. The synthetic all-ones `base` annotation follows the same rule in unpartitioned runs.

The approved memory refactor permits separate bounded annotation reads for query batches instead of requiring one combined all-query matrix. It resolves output rows before score allocation and accumulates scores only for those rows, while retaining every eligible reference-SNP contributor and the existing all/common counts and overlap universes. Non-output endpoints still contribute to output endpoints. `w_ld` keeps its distinct filtered-regression contributor set. A caller requesting all output rows remains supported.

Implementation status: at commit `a505c45`, direct computation still materializes the combined annotation matrix and full-reference-row score buffers before filtering output rows. The preceding batching/output-row contract is approved but not implemented; see [memory decisions](annotation-memory-decisions.md#confirmed-direct-ld-score-batching-and-output-row-accumulation) and the [refactor plan](../plans/2026-09-10-annotation-workflow-memory.md). The public LD-score output remains one aggregated baseline Parquet file and one optional aggregated query Parquet file, not separate chromosome files.

The gene-index builder applies the same rule to its fixed common payload:
supplied baseline columns and the regression mask share one PLINK traversal.
Its separate PLINK calls for bounded gene-atom column batches are intentional:
they construct the stored atom operator without materializing the full
SNP-by-atom matrix and are not a second computation of the regression weights.

## Exact gene LD-score index

The v1 gene index defaults written/regression rows to bundled HapMap3 SNPs
minus `mhc-and-centromeres`, intersected with retained reference rows. A custom
`--regr-snps-file` replaces the HapMap3 candidate set, and the selected
`--regr-snps-exclude-regions` policy is then subtracted. This is a **regression-row
policy**, not an HM3-only reference panel:

- baseline and focal LD scores still receive contributions from the broad
  retained PLINK universe;
- non-HM3, MHC, and pericentromeric SNPs remain eligible contributors;
- counts and overlap sufficient statistics use the broad retained universe;
- only `regression_ld_scores`/`w_ld` uses the filtered regression set as both
  rows and contributors.

The policy is fixed when the index is built. Indexed `ldscore` does not accept
an online regression-SNP or region override; build a complete new index for a
different policy.

## Consequence for region exclusions

For canonical LDSC-compatible artifacts, MHC and pericentromeric exclusions
are applied to the regression-SNP/output-row universe. They must not remove
SNPs from the baseline/query LD-score contributor universe, because doing so
changes the reference/causal universe and the numerical LD scores.

The `w_ld` exception is deliberate: it is an overcounting weight calculated
within regression SNPs, so it follows the filtered regression-SNP universe.

The original LDSC paper removes long-range-LD and pericentromeric variants
from LD-score regression while retaining the former in LD-score estimation.
The stratified-LDSC supplement separately defines HapMap3 regression SNPs and
a broad 1000 Genomes reference/causal universe. See:

- [Bulik-Sullivan et al. 2015](https://pmc.ncbi.nlm.nih.gov/articles/PMC4495769/)
- [Finucane et al. 2015](https://www.nature.com/articles/ng.3404)
- [official LDSC partitioned-heritability documentation](https://github.com/bulik/ldsc/wiki/Partitioned-Heritability)

The legacy implementation confirms the mechanism: it computes LD scores over
the full genotype/annotation matrix, applies `--print-snps` only when writing
rows, and then calculates `M` and `M_5_50` from the full annotation matrix.
