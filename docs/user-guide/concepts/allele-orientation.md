# Allele orientation

Last updated on: 2026-11-02

Every file with `A1` and `A2` columns names the two alleles of a SNP, but which allele is `A1` depends on the file.
In summary statistics `A1` is the effect allele; in files the package derives from reference-panel genotypes it is the minor allele; in inputs and in `annotate` output it keeps the order of the source.
This page states the four rules behind these meanings and how each command applies them.
The meaning in each file type is listed in [What A1 means in each file type](../../reference/formats/column-schema.md#what-a1-means-in-each-file-type).

## 1. Matching never depends on allele order

- In allele-aware modes (`rsid_allele_aware`, `chr_pos_allele_aware`), SNP identity uses the two alleles as an unordered pair, folded over the DNA strand: `A`/`G`, `G`/`A`, `T`/`C` and `C`/`T` give one identity (`normalize_allele_set` in `src/ldsc/_core/snp_identity.py`).
  Swapping `A1` and `A2` in a file therefore never changes which SNPs match.
- Allele order matters only for signed values, and each signed use aligns the alleles first:
  - in allele-aware modes, `h2` and `partitioned-h2` reverse the sign of `Z` where the summary statistics' `A1` is the LD-score `A2`;
  - `rg` aligns the second trait's `Z` to the first trait's `A1`;
  - `query-r2-matrix` reports `sign_r` and `r` for the dosages of each endpoint's `A1` in the pair file, which need not be the minor allele.
- LD scores use only R² and the folded MAF, so the rules below change no LD score, SNP count or regression result.

## 2. What A1 means depends on whether the package knows allele frequencies

- **Summary statistics:** `A1` is the effect allele, the allele the signed statistic refers to.
  `munge-sumstats` never guesses or reorients it (see [allele columns](../../reference/munge-sumstats.md#allele-columns)).
- **Files derived from PLINK reference-panel genotypes:** `A1` is the minor allele in the retained reference individuals (after `--keep-indivs-file`, counting non-missing calls), so `MAF` is the frequency of `A1` and at most 0.5.
  When both alleles have frequency exactly 0.5, the `.bim` order is kept.
  - `build-r2-matrix` writes its sidecars and correlation signs in this orientation.
  - `ldscore` writes the reference panel's `A1` and `A2` in every LD-score directory, whether the LD scores come from a PLINK reference panel, an R² matrix or a gene LD-score index, and in every SNP identifier mode; the alleles of annotation files take part only in matching.
  - `ldscore --plink-prefix ... --export-ref-metadata` writes `ref_metadata/chr<N>_meta.tsv.gz` in this orientation, so these files match the sidecars that `build-r2-matrix` writes for the same panel.
  - `build-gene-ldscore-index` stores its SNP rows in this orientation; indexes built before format version 5 must be rebuilt.
- **Inputs and `annotate` output:** `A1` and `A2` keep the order of the source and only identify SNPs.
  This covers PLINK `.bim` files, annotation files, SNP restriction files, the packaged HM3 map, `query-r2-matrix` pair files and `quantile-h2` inputs.
  `annotate` reads no genotypes, so it writes the alleles of its annotation SNP grid in their source order.
- LD-score directories written before these rules may carry `.bim` or annotation-file allele order.
  They still load, and because matching ignores order, regression results from them are unchanged.
- The name `A1` keeps both meanings; the [glossary](../../dev/DOCS_VOCABULARY.md) defines it for each kind of file.

## 3. A frequency column is always the frequency of A1

- `FRQ` in summary statistics is the frequency of the effect allele and may exceed 0.5; it is never folded.
- `MAF` in sidecars and `--export-ref-metadata` files is the frequency of `A1`, the minor allele, so it is at most 0.5.

## 4. Missing alleles are handled by the input's role

- Inputs that attach values or filter SNPs fall back to base identity in allele-aware modes:
  - When no annotation file has `A1` and `A2` (in `ldscore`, `annotate --ref-snp-annot-sources`, `build-gene-ldscore-index` baselines and `quantile-h2` resupplied annotations), annotation SNPs are matched by the base mode, and the workflow log names the files and the base mode.
    When at least one of the annotation files read together has both columns, the shared annotation SNP grid takes its alleles and matching stays allele-aware.
  - A `quantile-h2` target annotation file without alleles is matched the same way, with a log line naming it.
  - SNP restriction files (`--ref-panel-snps-file`, `--regr-snps-file`, munging keep-lists and gene-index restrictions) without allele columns filter by base identity, and the log names each such file and the base mode; a single allele column is an error.
- Inputs that define the query or the SNP universe must carry alleles in allele-aware modes:
  - A `query-r2-matrix` pair file without `A1` and `A2` for both endpoints stops the run; `--snp-identifier chr_pos` or `rsid` returns `r2` with empty `sign_r` and `r` instead.
  - `quantile-h2` reference metadata without `A1` and `A2` stops the run for an allele-aware fit; the R² matrix sidecars and the `--export-ref-metadata` files carry both.
    In base modes these alleles are optional and ignored, and duplicate base identities stop the run.
- Regression combines artifacts that declare different SNP identifier modes only with `--allow-identity-downgrade`; allele-less legacy `.sumstats(.gz)` files in `h2` and `partitioned-h2` are the exception (see [allele-less legacy files](../../reference/formats/legacy-sumstats.md#allele-less-legacy-files)).
