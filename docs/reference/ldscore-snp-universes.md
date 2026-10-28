# Reference and regression SNPs

Last updated on: 2026-10-27

This table distinguishes the SNP sets used by canonical LD score regression (LDSC) and by stratified LD score regression (S-LDSC; Finucane et al. 2015).
It is especially important for MHC and pericentromeric exclusions.

| Quantity | Appropriate SNP set |
|---|---|
| Reference SNPs $k$ that contribute to baseline and query LD scores | Broad retained 1000 Genomes/PLINK reference SNPs, including MHC and pericentromeric SNPs |
| Written LD-score rows and regression SNPs $j$ | HapMap3 intersected with the complement of MHC, pericentromeric, and other declared regression-SNP exclusions |
| Regression-weight LD score `w_ld` | Computed within the filtered HapMap3 regression SNPs; excluded-region SNPs neither enter the regression nor contribute to `w_ld` |
| All-SNP counts, common-SNP counts, and baseline-overlap counts | Broad retained reference SNPs, which are also the assumed causal SNPs (with the ordinary MAF threshold for common-SNP counts), not only regression SNPs |
| Final `h2`, `rg`, and `partitioned-h2` observations | Intersection of the filtered regression SNPs, written LD-score rows, regression weights, and munged summary statistics |

## Chromosome coverage versus SNP support

Direct query chromosome coverage is established from validated baseline/reference inputs before SNP, sample, or MAF filtering.
Their chromosome sets must match exactly.
`@` declares autosomes 1–22; exact paths and globs select actual artifacts whose contents establish scope.
Every selected gene must lie within scope after explicit gene exclusions.
Coverage failures abort the batch, whereas a valid chromosome whose retained SNPs for computation become empty supplies a measured zero-support outcome.
Direct gene support uses the prepared baseline/reference intersection; indexed support uses immutable gene-to-atom counts.
Unknown or unevaluated support stays blank.
Query LD-score variance is assessed separately on the regression SNPs.
See [coverage diagnostics](results/gene-list-diagnostics.md#chromosome-scope-and-pathway-coverage).

## Projection and traversal contract

Within each direct query execution batch, baseline/control annotations, the active query columns, and the binary regression-SNP mask share one numerical LD traversal per chromosome.
PLINK computes each genotype-correlation block once for that batch; an R² panel decodes bounded pair chunks.
The next execution batch repeats reference preparation and genotype/R² traversal.
This repeated work is deliberate: completed query batches are written and released before the next batch is prepared.
Shared baseline counts and overlap blocks are retained once, not added repeatedly across execution batches.
The synthetic all-ones `base` annotation uses one no-query calculation in unpartitioned runs.

The implementation reads bounded annotation tiles and resolves output rows before score allocation.
It accumulates scores only for output rows while retaining every eligible reference SNP and the established SNP sets for all-SNP counts, common-SNP counts and overlap.
Non-output endpoints still contribute to output endpoints.
`w_ld` keeps its distinct set of contributing SNPs, the filtered regression SNPs.
A caller requesting all output rows remains supported.
Counts sum binary 0/1 and normalized float32 continuous annotation values in float64, preserving existing validation tolerances across query widths and array layouts.
Packed storage and selected decoding follow the [annotation format policy](../dev/design/memory-and-cpu.md#annotation-format-policy).

The [sequential query-batch design](../dev/design/memory-and-cpu.md#direct-ld-score-computation) supersedes the earlier requirement for one traversal across all queries. Public output contains one shared genome-wide baseline file and one or more genome-wide query files, with chromosome row groups and the required ordered `query_batches` manifest. Generation batch width and chromosome workers affect execution memory, not SNP membership or scientific definitions. Sources: `LDScoreCalculator.run` in [_ldscore/calculator.py](../../src/ldsc/_ldscore/calculator.py), `ProjectionAccumulator` in [_kernel/ldscore_projection.py](../../src/ldsc/_kernel/ldscore_projection.py), and `annotation_statistics` in [_annotation/overlap.py](../../src/ldsc/_annotation/overlap.py). See the [batch verification](../dev/design/memory-and-cpu.md#diagnostics-and-verification).

The gene LD-score index builder applies the same rule to its fixed common payload:
supplied baseline columns and the regression mask share one PLINK traversal.
Its separate PLINK calls for bounded gene-atom column batches are intentional:
they construct the stored atom operator without materializing the full
SNP-by-atom matrix and are not a second computation of the regression weights.

## Exact gene LD-score index

The gene LD-score index defaults the written rows and regression SNPs to the bundled HapMap3 SNPs minus `mhc-and-centromeres`, intersected with the retained reference SNPs.
A custom `--regr-snps-file` replaces the HapMap3 candidate set, and the selected `--regr-snps-exclude-regions` policy is then subtracted.
This is a **regression-SNP policy**, not an HM3-only reference panel:

- baseline and query LD scores still receive contributions from the broad retained reference SNPs in the PLINK panel;
- non-HM3, MHC, and pericentromeric SNPs remain eligible reference SNPs;
- counts and overlap sufficient statistics use the broad retained reference SNPs;
- only `regression_ld_scores`/`w_ld` uses the filtered regression SNPs both as rows and as contributing SNPs.

The policy is fixed when the index is built. Indexed `ldscore` does not accept
an online regression-SNP or region override; build a complete new index for a
different policy.

## Consequence for region exclusions

For canonical LDSC-compatible artifacts, MHC and pericentromeric exclusions are applied to the regression SNPs and output rows.
They must not remove SNPs from the reference SNPs that contribute to baseline and query LD scores, because doing so changes the reference (assumed causal) SNPs and the numerical LD scores.

The `w_ld` exception is deliberate: it is an overcounting weight calculated within regression SNPs, so it follows the filtered regression SNPs.

The first LD score regression paper removes long-range-LD and pericentromeric variants from LD score regression while retaining the former in LD-score estimation.
The S-LDSC supplement separately defines HapMap3 regression SNPs and broad 1000 Genomes reference SNPs, which are also the assumed causal SNPs.
See:

- [Bulik-Sullivan et al. 2015](https://pmc.ncbi.nlm.nih.gov/articles/PMC4495769/)
- [Finucane et al. 2015](https://www.nature.com/articles/ng.3404)
- [legacy LDSC partitioned-heritability wiki](https://github.com/bulik/ldsc/wiki/Partitioned-Heritability)

The legacy implementation confirms the mechanism: it computes LD scores over
the full genotype/annotation matrix, applies `--print-snps` only when writing
rows, and then calculates the all-SNP and common-SNP counts (`M` and `M_5_50`)
from the full annotation matrix.
