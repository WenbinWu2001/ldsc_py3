# LD-score SNP-universe contract

Last updated on: 2026-08-03

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

## Exact gene LD-score index

The v1 gene index fixes written/regression rows to bundled HapMap3 SNPs minus
`mhc-and-centromeres`, intersected with retained PLINK rows. This is a fixed
**regression-row policy**, not an HM3-only reference panel:

- baseline and focal LD scores still receive contributions from the broad
  retained PLINK universe;
- non-HM3, MHC, and pericentromeric SNPs remain eligible contributors;
- counts and overlap sufficient statistics use the broad retained universe;
- only `regression_ld_scores`/`w_ld` uses the filtered regression set as both
  rows and contributors.

Indexed mode does not accept a custom regression-SNP file. A different
regression-row policy requires direct `ldsc ldscore` mode.

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
