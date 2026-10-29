# Reference and regression SNPs

Last updated on: 2026-10-29

This table distinguishes the SNP sets used by canonical LD score regression (LDSC) and by stratified LD score regression (S-LDSC; Finucane et al. 2015).
It is especially important for MHC and pericentromeric exclusions.

| Quantity | Appropriate SNP set |
|---|---|
| Reference SNPs $k$ that contribute to baseline and query LD scores | Broad retained 1000 Genomes/PLINK reference SNPs, including MHC and pericentromeric SNPs |
| Written LD-score rows and regression SNPs $j$ | HapMap3 intersected with the complement of MHC, pericentromeric, and other declared regression-SNP exclusions |
| Regression-weight LD score `w_ld` | Computed within the filtered HapMap3 regression SNPs; excluded-region SNPs neither enter the regression nor contribute to `w_ld` |
| All-SNP counts, common-SNP counts, and baseline-overlap counts | Broad retained reference SNPs, which are also the assumed causal SNPs (with the ordinary MAF threshold for common-SNP counts), not only regression SNPs |
| Final `h2`, `rg`, and `partitioned-h2` observations | Intersection of the filtered regression SNPs, written LD-score rows, regression weights, and munged summary statistics |

## Reference SNPs on the annotation SNP grid

The annotation files given to `ldscore` define the annotation SNP grid: their SNP rows are the candidate SNPs.
Baseline and prebuilt query files must list the same SNP rows in the same order, and BED and gene-list queries are projected onto the baseline rows, so all annotation files of one run share one grid.
The reference SNPs are the grid SNPs among the reference-panel SNPs that pass the panel's filters: identity cleanup, `--ref-panel-snps-file`, `--keep-indivs-file`, genotype QC and `--maf-min`.
Reference-panel SNPs without an annotation row are excluded at the `annotation_snp_grid` step, and grid rows without a reference-panel SNP are dropped at the same step.
No annotation value is filled in or imputed; a missing annotation value is rejected when the file is read.

The exclusion does not stop the run.
It is a deliberate exception, confirmed by the maintainer on 2026-10-28, to the repository rule that unrequested filtering that can change the analysis must fail ([AGENTS.md](../../AGENTS.md), "Existing warning, error, and diagnostic policy").
Supplying annotation files is the explicit choice of candidate SNPs.
Legacy baseline-LD continuous annotations exist only for 1000 Genomes Phase 3 SNPs, so a newer panel's extra SNPs cannot be annotated, and requiring an option on every such run was rejected.
Legacy LDSC instead required the `.annot` rows to equal the `.bim` rows (legacy `ldsc.py`, lines 156-158).

The exclusion is never silent:

- `diagnostics/reference_snp_funnel.tsv` records every chromosome's reference-SNP funnel, including skipped chromosomes;
- when the `annotation_snp_grid` step removes at least one reference-panel SNP, `diagnostics/ldscore.log` receives one WARNING with the total, the number of chromosomes and the funnel path, and a successful command-line run prints the same line to stderr;
- Python callers of `LDScoreCalculator.run` without an output directory find the same table in `LDScoreResult.reference_snp_funnel`, and the WARNING in the package log.

The warning has no threshold.
`metadata.json` records the number of reference SNPs on the grid as `snp_universe_policy.ld_reference_universe.retained_snp_count`, and the funnel describes how they were selected.

### Reference-SNP funnel

Each chromosome has two funnels in `diagnostics/reference_snp_funnel.tsv`, one row per step, with the columns `funnel`, `chrom`, `backend`, `step_index`, `step`, `removed` and `remaining`.
Steps follow the backend's own execution order.
Each backend always writes the same rows.
A step that did not run, because its option was not set, the identifier mode skips it, or preparation stopped earlier, has `removed` equal to `NA` and leaves `remaining` unchanged.

| Funnel | Backend | Steps in order |
| --- | --- | --- |
| `reference_panel` | PLINK | `panel_rows` (BIM rows on the chromosome), `identity_cleanup` (`rsid` and `chr_pos` modes), `ref_panel_snps_file`, `allele_aware_identity_cleanup` (allele-aware modes), `annotation_snp_grid`, `genotype_qc` (monomorphic or all-missing SNPs after `--keep-indivs-file`), `maf_floor` |
| `reference_panel` | R² matrix | `panel_rows` (sidecar rows), `ref_panel_snps_file`, `maf_floor`, `allele_aware_identity_cleanup`, `annotation_snp_grid` |
| `annotation_snp_grid` | PLINK | `grid_rows`, `annotation_snp_grid` (grid rows with no remaining reference-panel SNP), `genotype_qc`, `maf_floor` |
| `annotation_snp_grid` | R² matrix | `grid_rows`, `annotation_snp_grid` |

For PLINK panels the grid step runs before genotype QC, so a reference-panel SNP off the grid is counted even when genotype QC would have removed it, and the genotype-QC and MAF counts cover grid SNPs only.
For R² matrices the MAF floor comes from the sidecar and runs before the grid step.
Both funnels of a chromosome end at the same count: the reference SNPs that enter LD scores, SNP counts and overlap.
The counts come from one reference preparation per chromosome per run; later query batches and the gene-list support pass add nothing.
Sources: `ReferenceSNPFunnel` in [_refpanel/funnel.py](../../src/ldsc/_refpanel/funnel.py), `_align_annotations()` and each panel's `prepare_chromosome()` in [_refpanel/adapters.py](../../src/ldsc/_refpanel/adapters.py), and `_DirectProducer.reduce()` in [_ldscore/calculator.py](../../src/ldsc/_ldscore/calculator.py).

### Empty chromosomes

If the `annotation_snp_grid` step removes every reference-panel SNP left on a chromosome, `ldscore` fails and names every such chromosome with its grid-row and reference-panel SNP counts.
A whole chromosome without one shared SNP usually indicates a genome-build, allele or SNP-identifier mismatch; to exclude a chromosome on purpose, leave it out of the annotation files.
Gene-list runs raise this error in the gene-support pass, before LD-score computation; other runs raise it after the first query batch's reference preparation.
A chromosome emptied by an earlier step, such as `--ref-panel-snps-file`, or by genotype QC is skipped with a WARNING in `diagnostics/ldscore.log`, and its funnel shows the step that removed its SNPs.
A run in which no chromosome keeps a reference SNP fails.

## Chromosome coverage versus SNP support

Direct query chromosome coverage is established from validated baseline/reference inputs before SNP, sample, or MAF filtering.
Their chromosome sets must match exactly.
`@` declares autosomes 1–22; exact paths and globs select actual artifacts whose contents establish scope.
Every selected gene must lie within scope after explicit gene exclusions.
Coverage failures abort the batch.
A valid chromosome emptied by `--ref-panel-snps-file` or QC supplies a measured zero-support outcome; one emptied at the annotation SNP grid step fails, as described under [empty chromosomes](#empty-chromosomes).
Direct gene support counts the reference SNPs on the annotation SNP grid; indexed support uses immutable gene-to-atom counts.
Unknown or unevaluated support stays blank.
Query LD-score variance is assessed separately on the regression SNPs.
See [coverage diagnostics](results/gene-list-diagnostics.md#chromosome-scope-and-pathway-coverage).

## Projection and traversal contract

Within each direct query execution batch, baseline/control annotations, the active query columns, and the binary regression-SNP mask share one numerical LD traversal per chromosome.
PLINK computes each genotype-correlation block once for that batch; an R² matrix decodes bounded pair chunks.
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

The [sequential query-batch design](../dev/design/memory-and-cpu.md#direct-ld-score-computation) supersedes the earlier requirement for one traversal across all queries. Public output contains one shared genome-wide baseline file and one or more genome-wide query files, with chromosome row groups and the required ordered `query_batches` manifest. Generation batch width and chromosome workers affect execution memory, not SNP membership or scientific definitions. Sources: `LDScoreCalculator.run` in [_ldscore/calculator.py](../../src/ldsc/_ldscore/calculator.py), `ProjectionAccumulator` in [_ldscore/compute.py](../../src/ldsc/_ldscore/compute.py), and `annotation_statistics` in [_annotation/overlap.py](../../src/ldsc/_annotation/overlap.py). See the [batch verification](../dev/design/memory-and-cpu.md#diagnostics-and-verification).

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

The index build applies the annotation SNP grid rule above to its baseline annotations and stores each chromosome's reference-SNP funnel in the index.
Its `annotation_snp_grid` funnel starts from every baseline row, before the build's baseline/PLINK identifier intersection.
The build logs and prints the same WARNING, pointing to `diagnostics/build-gene-ldscore-index.json`.
Indexed `ldscore` copies the stored funnel to its own `diagnostics/reference_snp_funnel.tsv` and repeats the WARNING; it does not recompute the funnel.

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
