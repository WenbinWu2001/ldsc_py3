# LD-score SNP-universe contract

Last updated on: 2026-08-06

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

## Projection and traversal contract

For each chromosome, ordinary direct `ldscore` constructs one combined
annotation matrix whose columns are ordered as baseline, query/control, then
the binary regression-SNP mask. The PLINK backend computes genotype-correlation
blocks once and projects all of those columns in the same `ldScoreVarBlocks`
call; the parquet-R2 backend likewise streams stored pairs once over the same
combined matrix. The result is split afterward into partitioned LD scores and
`w_ld`. This also applies to an unpartitioned run, whose synthetic all-ones
`base` column is projected together with the regression mask.

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
