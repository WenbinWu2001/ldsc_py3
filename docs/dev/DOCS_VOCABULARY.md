# Documentation Vocabulary

Last updated on: 2026-10-28

This is the glossary of the LDSC3 documentation.
Each term below has one meaning and one wording across the user guide, the reference, the developer guide, the README, docstrings, command help and messages.
Commands, flags, arguments, columns and file names appear in code font exactly as the code spells them, even where the prose term differs.
Rules for writers and terms still under proposal follow the term tables.

## Names

| Term | Definition | Avoid | API/CLI identifier |
| --- | --- | --- | --- |
| LD score regression (LDSC) | The method: regression of GWAS $\chi^2$ statistics on LD scores, used to estimate SNP heritability, genetic correlation and partitioned heritability. In prose, "LDSC" means the method only. | "LD-score regression", "LD Score Regression" | |
| LDSC3 | This Python 3 package (distribution `ldsc`, version 2.0b0); a provisional name. | "the refactored package", "ldsc3", "the new LDSC" | import `ldsc`; command `ldsc` |
| legacy LDSC | The original Python 2 software (bulik/ldsc) that LDSC3 replaces. Before a file or format, "legacy" means one that legacy LDSC reads or writes, as in "legacy LD-score files". | "LDSC2", "original LDSC", "upstream LDSC", "Python 2 LDSC", "the legacy codebase" | `convert-ldsc2-ldscores` |
| partitioned heritability | The analysis that decomposes SNP heritability across annotations fitted jointly with the baseline, and its results. | "partitioned LDSC", "pLDSC", "stratified LDSC", "partitioned h2" in prose | `partitioned-h2`; `RegressionRunner.estimate_partitioned_h2` |
| stratified LD score regression (S-LDSC) | The published method that estimates partitioned heritability (Finucane et al. 2015); used only when naming that method. | "S-LDSC" before the full name | |
| cell-type enrichment | The analysis that fits one query annotation per cell type and compares their heritability contributions. | "cell-specific LDSC" | `partitioned-h2` |
| LDSC-SEG | The published method for cell-type enrichment from specifically expressed genes (Finucane et al. 2018); used only when naming that method. | | |

The digit in LDSC3 is the Python version, not the software version.
Write "legacy LDSC", not a version number, for the original software; "legacy" alone refers only to legacy LDSC and its files, never to an earlier LDSC3 version.
`ldsc` in code font is the command and the import name, never a prose name.
The software title in `README.md` and `CITATION.cff` is the distribution name, `ldsc`.
On each page that names S-LDSC, the first mention gives the full name and citation, for example "partitioned heritability, estimated with stratified LD score regression (S-LDSC; Finucane et al. 2015)", so that searches for "S-LDSC" find the page.
`partitioned-h2` in code font names the command; prose uses "partitioned heritability".

## Approved terms

### SNP identity and genome build

| Term | Definition | Avoid | API/CLI identifier |
| --- | --- | --- | --- |
| SNP identifier mode | The setting that decides which fields identify a SNP when inputs are matched: `rsid`, `rsid_allele_aware`, `chr_pos`, or `chr_pos_allele_aware` (default). | "identity mode", "SNP identity mode", "matching mode" | `--snp-identifier`; `GlobalConfig.snp_identifier` |
| SNP identity | The key a SNP receives under the SNP identifier mode. | "SNP key"; "SNP ID" when the key is meant | |
| base mode | `rsid` or `chr_pos`: identity without alleles. | "allele-free mode", "allele-blind mode" | |
| allele-aware mode | `rsid_allele_aware` or `chr_pos_allele_aware`: identity that also uses the two alleles. | "allele mode" | |
| coordinate-family mode | `chr_pos` or `chr_pos_allele_aware`. | "coordinate mode", "position mode" | |
| rsID-family mode | `rsid` or `rsid_allele_aware`. | "rsID mode" when both are meant | |
| genome build | The reference assembly that positions refer to: hg19 (GRCh37) or hg38 (GRCh38); `auto` infers it. | "assembly", "build version", "reference build" | `--genome-build`; `GlobalConfig.genome_build`; `ldsc infer-genome-build` |
| liftover | Converting positions between hg19 and hg38 with a chain file, or with the HM3 map during munging. | "lift-over", "lift over", "coordinate conversion" | `--liftover-chain-file`; `--output-genome-build` |
| A1, A2 | In summary statistics, A1 is the allele the signed statistic refers to; in reference-panel, LD-score and annotation files, A1 is the minor allele. A2 is the other allele. A1 is not necessarily the genome reference allele. | "reference allele" and "ref/alt" as names for A1 and A2 | columns `A1`, `A2`; `--a1`, `--a2` |
| one-based position | A SNP position that counts the first base as 1. | "1-indexed" | column `POS` |
| zero-based half-open interval | A BED interval whose start counts from 0 and whose end is excluded. | "0-indexed interval" | BED files |

### Reference data and LD

| Term | Definition | Avoid | API/CLI identifier | Symbol |
| --- | --- | --- | --- | --- |
| reference panel | Genotype data (PLINK) or an R² panel from which LD scores are computed. | "ref panel", "reference genotypes", "LD panel" | `--plink-prefix`; `--r2-dir`; `RefPanelConfig` | |
| R² panel | Package-built Parquet files of pairwise unbiased R² with one sidecar per chromosome, written by `build-r2-panel`. | "parquet R² reference", "R² table", "index-format panel", "R² parquet panel" | `build-r2-panel`; `--r2-dir`; `--panel-dir`; `R2Panel` | |
| sidecar | The per-chromosome metadata file of an R² panel that records SNP identity, CM and MAF. | "metadata file", "meta file" | `chr<N>_meta.tsv.gz` | |
| R² | Squared correlation between two SNPs' genotypes; written "R²" (or $R^2$ in math) in Markdown prose, and "R2" in code text (docstrings, CLI help, messages, comments) and identifiers. | "R2" in Markdown prose, "r2" | column `R2` | |
| unbiased R² | The sample-size-corrected squared correlation $r^2 - (1 - r^2)/(n - 2)$ (`src/ldsc/_refpanel/query.py`); it can be negative. | "adjusted R²", "bias-adjusted R²" | column `R2`; `ldsc:r2_bias="unbiased"` | $\tilde r^2$ |
| LD window | The maximum distance between two SNPs whose R² enters an LD score, in SNPs, kb or cM. | "window size" | `--ld-wind-snps`; `--ld-wind-kb`; `--ld-wind-cm` | |
| LD score | The sum of unbiased R² between a SNP and the reference SNPs in its LD window, weighted by annotation values. Write "LD score" as a noun and "LD-score" before another noun. | "LDscore", "L2" | `ldscore`; `ldscore.baseline.parquet` | $\ell(j, C)$ |
| LD-score directory | The output of `ldscore` or `convert-ldsc2-ldscores`: `metadata.json` with baseline, query and overlap Parquet files; the input to regression. | "LD-score artifact", "LD-score result directory", "canonical LD-score directory" | `--ldscore-dir`; `LDScoreSource` | |
| gene LD-score index | Precomputed data that give exact gene-list LD scores for one baseline, reference panel and gene coordinate catalog. | "gene index", "exact index" | `build-gene-ldscore-index`; `--gene-ldscore-index-dir` | |
| HM3 map | The packaged HapMap3 SNP table with hg19 and hg38 positions; by default it defines the regression SNPs and restricts munged summary statistics, and munging uses it for quick liftover. Genome-build inference uses the separate `hm3_chr_pos_reference.tsv.gz`. | "HM3 list", "HapMap3 file" | `src/ldsc/data/hm3_curated_map.tsv.gz`; `load_hm3_curated_map()` | |
| genetic map | A table of per-position genetic distances in cM, used for cM LD windows. | "recombination map" | `--genetic-map-hg19-sources`; `--genetic-map-hg38-sources` | |
| region-exclusion preset | A named set of genomic regions (MHC, centromeres) removed from the regression SNPs. | "blacklist", "masked regions" | `--regr-snps-exclude-regions` | |

### SNP sets and counts

| Term | Definition | Avoid | API/CLI identifier | Symbol |
| --- | --- | --- | --- | --- |
| reference SNPs | The reference-panel SNPs retained after filtering; they contribute to LD scores, SNP counts and annotation overlap. | "LD contributors", "contributor universe", "count universe" | | |
| reference-SNP funnel | The per-chromosome count of SNPs that each reference-preparation step removed, from the reference panel's rows to the reference SNPs. | | `diagnostics/reference_snp_funnel.tsv`; `LDScoreResult.reference_snp_funnel`; gene-index payload `reference_snp_funnel.parquet` | |
| regression SNPs | The SNPs whose rows enter the regression and its weights (the HM3 map by default, after region exclusion). | "regression rows", "regression-row universe", "weight SNPs" | `--regr-snps-file`; `--regr-snps-exclude-regions` | |
| common-SNP count | The number of reference SNPs with MAF at least `common_maf_min` (default 0.05, inclusive), per annotation; the default heritability denominator. | "M_5_50", "common-variant count", "common counts" | `--count-kind common`; `--common-maf-min` | $M_C$ |
| all-SNP count | The number of all reference SNPs per annotation. | "full count", "M count" | `--count-kind all` | $M_C$ |

### Annotations

| Term | Definition | Avoid | API/CLI identifier | Symbol |
| --- | --- | --- | --- | --- |
| annotation | A per-SNP value column used to partition heritability; binary (0/1) or continuous. | "feature", "track"; "category" except as the result column `category` | | $a_C(j)$ |
| baseline annotation | An annotation that is in every fitted model, such as the baseline-LD set. | "baseline category", "background annotation" | `--baseline-annot-sources` | |
| query annotation | An annotation tested one at a time on top of the baseline. | "focal annotation", "test annotation" | `--query-annot-sources`; `--query-annot-bed-sources`; `--query-annot-gene-list-sources` | |
| continuous annotation | An annotation with non-binary values. | "quantitative annotation" | `quantile-h2` | |
| annotation SNP grid | The SNP rows on which annotations are defined; they are the candidate SNPs. `annotate` takes them from a BIM grid (PLINK `.bim` files) or an annotation grid (existing annotation files) and writes them as its annotation rows; `ldscore` and `build-gene-ldscore-index` read them from the annotation files. It should list the reference panel's SNPs: the reference SNPs are the grid SNPs kept after reference-panel filtering, and reference-panel SNPs not on it are excluded and counted in the reference-SNP funnel. Before another noun, write "annotation rows" or "annotation row SNPs". | "SNP grid" alone; "annotation grid" except for the kind selected by `--ref-snp-annot-sources`; "SNP universe" for this meaning | `--ref-snp-bim-sources`; `--ref-snp-annot-sources`; `annotation_snp_grid` (funnel name and step) | |
| gene list | A file of gene identifiers that becomes one query annotation. | "gene set" for the file, "pathway file" | `--query-annot-gene-list-sources` | |
| gene coordinate catalog | The table of one-based gene coordinates that resolves gene identifiers to intervals. | "gene catalog", "coordinate catalog", "gene annotation file" | `--gene-coordinate-file` | |
| control gene list | A gene list added to the baseline as column `gene_control`. | "background genes" | `--control-gene-list-file` | |
| padding | Base pairs added on both sides of each gene interval. | "flank", "window" | `--padding-bp` | |
| query batch | The query annotations computed together in one LD-score pass. | "chunk" | `--query-batch-size` | |

### Summary statistics

| Term | Definition | Avoid | API/CLI identifier | Symbol |
| --- | --- | --- | --- | --- |
| summary statistics | Per-SNP GWAS association results; written in full in prose. | "sumstats" in prose | `munge-sumstats` | |
| raw summary statistics | The GWAS table given to `munge-sumstats`. | "raw sumstats", "input GWAS" | `--raw-sumstats-file` | |
| munged summary statistics | The output of `munge-sumstats`: Parquet with an identity footer, optionally also `.sumstats.gz`. | "curated summary statistics", "cleaned summary statistics" | `--sumstats-file`; `--sumstats-sources` | |
| sample size | The number of individuals behind each SNP's statistic. | | column `N`; `--N`; `--N-col` | $N$ |

### Regression quantities

| Term | Definition | Avoid | API/CLI identifier | Symbol |
| --- | --- | --- | --- | --- |
| SNP heritability | The proportion of phenotypic variance explained by the SNPs counted in the model's denominator. | "h2" in prose, "SNP-heritability", "chip heritability" | `h2`; `total_h2_obs`; `h2_obs` | $h^2_g$ |
| observed scale, liability scale | Scales of heritability for a binary trait; conversion uses sample and population prevalence. | "observed-scale h2" in prose | `*_obs`, `*_liab`; `convert-h2-scale` | |
| sample prevalence, population prevalence | The case fraction in the GWAS sample, and the trait prevalence in the population. | | `--samp-prev`; `--pop-prev` | $P$, $K$ |
| intercept | The regression intercept; it absorbs confounding and, for genetic covariance, sample overlap. | "LDSC constant" | `intercept`; `--intercept-h2`; `--intercept-gencov` | |
| ratio | $(\text{intercept} - 1)/(\text{mean } \chi^2 - 1)$ (`src/ldsc/_regression/estimators.py`). | "attenuation ratio" | `ratio` | |
| genetic correlation | The correlation of genetic effects between two traits. | "rg" in prose, "genetic correlation coefficient" | `rg` | $r_g$ |
| genetic covariance | The covariance of genetic effects between two traits. | "gencov" in prose | `gencov` | $\rho_g$ |
| coefficient | The per-annotation regression coefficient: the per-SNP heritability contributed by one unit of the annotation, conditional on the other annotations. | "tau" in prose, "effect size" | `coefficient` in `partitioned_h2.tsv`; `tau` in the standardized-coefficient table | $\tau_C$ |
| standardized coefficient | The coefficient scaled by the annotation standard deviation and the per-SNP heritability (Gazal et al. 2017). | "tau star" | `tau_star` | $\tau^*_C$ |
| proportion of heritability | The share of SNP heritability attributed to an annotation. | "heritability fraction" | `prop_h2` | |
| proportion of SNPs | The annotation's share of the SNP count. | "SNP fraction" | `prop_snps` | |
| enrichment | Proportion of heritability divided by proportion of SNPs (`src/ldsc/_regression/estimators.py`). | "fold enrichment", "enrichment ratio" | `enrichment` | |
| block jackknife | Standard-error estimation by deleting contiguous blocks of SNPs, 200 by default. | "jackknife blocks" for the method | `--n-blocks` | |

### Inputs, outputs and execution

| Term | Definition | Avoid | API/CLI identifier |
| --- | --- | --- | --- |
| path token | A string that names inputs: an exact path, a glob, or a chromosome suite. | "path pattern" for all three | |
| chromosome suite | A path containing `@`, expanded to chromosomes 1 to 22. | "chromosome pattern", "@ token" | `@` |
| output directory | The directory into which a workflow writes its fixed-name outputs. | "output prefix", "out dir" | `--output-dir` |
| diagnostics | Files under `diagnostics/` that explain filtering, issues and failures; not headline results. | "audit files" | `diagnostics/` |
| failure marker | The file left after a failed overwrite, removed by a successful retry. | "error file" | `RUN_FAILED.txt` |
| worker | A managed child process that runs one task (a chromosome or a model); each uses one BLAS thread. | "thread", "job", "process" for this unit | `--workers`; `workers=` |

## Identifier differences

These identifiers differ from the prose term or from each other; they are recorded here, not renamed.

- The coefficient is the column `coefficient` in `partitioned_h2.tsv` and `tau` in the standardized-coefficient table (`src/ldsc/_regression/outputs.py`).
- Prose "R²" corresponds to the identifiers `R2`, `r2_dir` and `build-r2-panel`.
- `--count-kind common|all` on the command line corresponds to `RegressionConfig.use_common_counts` in Python.
- Prose "gene coordinate catalog" corresponds to the flag `--gene-coordinate-file`.
- Prose "legacy LDSC" corresponds to the identifiers that contain `ldsc2`, such as `convert-ldsc2-ldscores`.

## Rules for writers

1. **One concept, one wording:** use the approved term every time; do not vary wording for style.
2. **One glossary:** this file is the only glossary for the documentation, README, docstrings, CLI help text and messages; add terms here and do not start another glossary.
3. **Read before writing:** read this file before every writing task, including small edits.
4. **Entry content:** each approved entry gives the term, a short definition, wordings to avoid, the API/CLI identifier when readers need to connect the term to code, and a symbol for quantities.
5. **Precedence:** approved terms take precedence over general style preferences; keep technically different concepts and verbs distinct.
6. **Approved versus proposed:** approved and proposed terms are kept apart; when a concept appears under several names, record a proposal with the preferred term and the variants with their locations.
7. **Author control:** only the maintainer approves proposals or edits and removes approved entries.
8. **Existing text:** violations in existing text are flagged with location and preferred term; rewriting them needs separate authorization; all newly written prose follows the approved entries.
9. **Scope:** consulting the glossary does not authorize broader rewriting.
10. **Identifiers stay as spelled:** commands, flags, arguments, columns and file names appear in code font exactly as the code spells them, even where the prose term differs.
11. **ASCII-only code:** every script (package source, tests, tools), including docstrings, comments, CLI help text and messages, uses ASCII characters only.
    Write `R2`, `h2`, `chi^2`, `+/-` and `->` there; Unicode symbols such as R², ± and → appear only in Markdown documentation.
    When a figure or printed output must show a symbol, produce it from ASCII source: matplotlib mathtext such as `R$^2$`, or a `\uXXXX` escape.
12. **No emoji:** emoji are not used anywhere, including documentation, code, messages, outputs and commit messages.
13. **Figure labels:** text drawn inside a figure may use condensed phrasings, such as "sumstats", where space is limited; captions and the surrounding prose use the approved terms.
14. **Verbatim copies:** documentation that reproduces code text verbatim, such as the CLI help tables in `docs/reference/cli.md`, the docstring summaries in `docs/reference/python-api.md`, and quoted messages, keeps the code's ASCII spelling (for example `R2`).

## Proposed terms

Add a row when a concept appears under several names: the concept, the proposed term, and each variant with its file and line.

| Concept | Proposed term | Variants and locations |
| --- | --- | --- |
| A private package under `src/ldsc/` that holds one area of the package, such as `_annotation` or `_ldscore` | domain package | "domain subpackage" (the removed reorganization specification); "domain packages" (`docs/dev/code-structure.md:36`, `docs/dev/data-flow.md:15`, `docs/dev/architecture.md:52`) |
| A level in the one-directional import order of the domain packages, enforced by `tests/test_layering.py` | tier | "tiers" (`docs/dev/code-structure.md:21`, `docs/dev/architecture.md:63`); "layer" already names the CLI, workflow and output layers in `docs/dev/architecture.md` and should not be used for this concept |
| The modules at the package root that belong to no domain package: `__init__`, `__main__`, `cli`, `config` and `errors` | root shell | "public shell" (`docs/dev/code-structure.md:19`, `docs/dev/architecture.md:113`); "shell" (`docs/dev/code-structure.md:57`); "root shell" (`tests/test_layering.py:8`) |
