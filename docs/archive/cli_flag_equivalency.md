# Old vs Refactored CLI Flag Equivalency

This document compares the legacy CLI surface in `munge_sumstats.py` and `ldsc.py` against the current refactored codebase.

Important scope note:

- Only LD-score estimation currently has a shipped refactored CLI: `ldsc_new.py` -> `ldscore/ldscore_new.py`.
- For `munge sumstats`, `heritability`, `partitioned LDSC`, `cross-trait regression`, and `cell type specific LDSC`, the refactor currently defines design targets in `class-and-features.md`, but it does not yet ship a concrete replacement parser with stable new flags.
- In those sections, `Not yet implemented` means there is no concrete refactored CLI flag in code today.
- `-` in the `old flag` column means the option is new in the refactored CLI and has no legacy counterpart.

## Munge Sumstats

Refactor status: planned as `ldsc munge-sumstats ...`, but no refactored CLI flags are implemented yet.

| old flag | new flag | meaning |
| --- | --- | --- |
| `--sumstats` | `Not yet implemented` | Raw GWAS summary-statistics input file. |
| `--N` | `Not yet implemented` | Override or supply total sample size. |
| `--N-cas` | `Not yet implemented` | Override or supply number of cases. |
| `--N-con` | `Not yet implemented` | Override or supply number of controls. |
| `--out` | `Not yet implemented` | Output prefix for munged `.sumstats.gz` and logs. |
| `--info-min` | `Not yet implemented` | Minimum INFO threshold used for QC filtering. |
| `--maf-min` | `Not yet implemented` | Minimum MAF threshold used for QC filtering. |
| `--daner-old` | `Not yet implemented` | Parse older `daner*` summary-statistics format with case/control N encoded in `FRQ_A_<Ncas>` and `FRQ_U_<Ncon>` headers. |
| `--daner-new` | `Not yet implemented` | Parse newer `daner*` format with `Nca` and `Nco` columns. |
| `--no-alleles` | `Not yet implemented` | Allow munging without allele columns; mainly useful for `h2` / partitioned `h2`, not `rg`. |
| `--merge-alleles` | `Not yet implemented` | Harmonize alleles to a reference `SNP A1 A2` file and drop mismatches. |
| `--n-min` | `Not yet implemented` | Minimum per-variant sample size threshold. |
| `--chunksize` | `Not yet implemented` | Chunk size for streaming large summary-statistics files. |
| `--snp` | `Not yet implemented` | Explicit raw-column name for SNP identifier. |
| `--N-col` | `Not yet implemented` | Explicit raw-column name for total sample size. |
| `--N-cas-col` | `Not yet implemented` | Explicit raw-column name for case count. |
| `--N-con-col` | `Not yet implemented` | Explicit raw-column name for control count. |
| `--a1` | `Not yet implemented` | Explicit raw-column name for effect allele. |
| `--a2` | `Not yet implemented` | Explicit raw-column name for non-effect allele. |
| `--p` | `Not yet implemented` | Explicit raw-column name for p-value. |
| `--frq` | `Not yet implemented` | Explicit raw-column name for allele frequency / MAF. |
| `--signed-sumstats` | `Not yet implemented` | Signed effect-stat column plus null value, for example `Z,0` or `OR,1`. |
| `--info` | `Not yet implemented` | Explicit raw-column name for INFO score. |
| `--info-list` | `Not yet implemented` | Comma-separated INFO columns to average before filtering. |
| `--nstudy` | `Not yet implemented` | Explicit raw-column name for number of studies. |
| `--nstudy-min` | `Not yet implemented` | Minimum number-of-studies threshold. |
| `--ignore` | `Not yet implemented` | Comma-separated raw columns to ignore during inference. |
| `--a1-inc` | `Not yet implemented` | Declare `A1` to be the increasing allele for signed statistics. |
| `--keep-maf` | `Not yet implemented` | Preserve MAF in the munged LDSC output if available. |

## LD Score Calculation

Refactor status: implemented in `ldsc_new.py` / `ldscore/ldscore_new.py`.

| old flag | new flag | meaning |
| --- | --- | --- |
| `--l2` | `implicit in ldsc_new.py` | Enter LD-score-estimation mode. The refactored CLI is a dedicated LD-score command, so no separate mode flag is needed. |
| `--out` | `--out` | Output prefix for `.l2.ldscore.gz`, `.w.l2.ldscore.gz`, `.l2.M`, `.l2.M_5_50`, and annotation-group manifest files. |
| `--bfile` | `--bfile` | Single PLINK prefix for the LD reference panel. |
| `-` | `--bfile-chr` | Chromosome-pattern PLINK prefix for per-chromosome reference-panel inputs. |
| `--extract` | `Not yet implemented` | Restrict the reference SNP set before LD-score computation. |
| `--keep` | `Not yet implemented` | Restrict the PLINK individual set before LD-score computation. |
| `--ld-wind-snps` | `--ld-wind-snps` | LD window size in number of SNPs. |
| `--ld-wind-kb` | `--ld-wind-kb` | LD window size in kilobases. |
| `--ld-wind-cm` | `--ld-wind-cm` | LD window size in centiMorgans. |
| `--print-snps` | `Not yet implemented` | Restrict which SNP rows are written to the legacy `.l2.ldscore` output. This is not the same as new `--regression-snps`. |
| `--annot` | `--baseline-annot`, `--query-annot`, `--baseline-annot-chr`, `--query-annot-chr` | Supply SNP-level annotation inputs. The refactor splits legacy annotation input into baseline and query groups and supports explicit chromosome-pattern inputs. |
| `--thin-annot` | `Not yet implemented` | Legacy shortcut for annotation-only matrices without `CHR/BP/SNP/CM`. The refactored parser expects full SNP-level metadata. |
| `--cts-bin` | `Not yet implemented` | Build annotations on the fly by binning continuous variables during LD-score computation. |
| `--cts-breaks` | `Not yet implemented` | Breakpoint definition for `--cts-bin`. |
| `--cts-names` | `Not yet implemented` | Display names for variables supplied through `--cts-bin`. |
| `--per-allele` | `Not yet implemented` | Compute per-allele LD scores weighted by `p(1-p)`. |
| `--pq-exp` | `Not yet implemented` | Compute LD scores with `(p(1-p))^a` weighting. |
| `--no-print-annot` | `Not yet implemented` | Suppress printing the generated annotation matrix. |
| `--maf` | `--maf` | MAF threshold on retained SNPs when MAF is available. |
| `--chunk-size` | `--chunk-size` | Chunk size for the legacy PLINK-backed block computations reused by the refactor. |
| `--yes-really` | `--yes-really` | Allow whole-chromosome LD windows. |
| `--frqfile` | `--frqfile` | Optional frequency / metadata file used to populate MAF and CM and to enable `.M_5_50` when possible. |
| `--frqfile-chr` | `--frqfile-chr` | Chromosome-pattern frequency / metadata inputs. |
| `--pickle` | `Not yet implemented` | Legacy option to write pickled LD-score files instead of gzipped text. |
| `-` | `--r2-table` | Use sorted parquet pairwise `R2` input instead of PLINK genotypes. |
| `-` | `--r2-table-chr` | Chromosome-pattern sorted parquet `R2` inputs. |
| `-` | `--snp-identifier` | Select SNP matching mode, either `chr_pos` or `rsID`. |
| `-` | `--genome-build` | Genome build used for parquet / `chr_pos` matching. |
| `-` | `--r2-bias-mode` | Declare whether parquet `R2` values are raw sample `r^2` or already unbiased. |
| `-` | `--r2-sample-size` | LD reference sample size used to convert raw parquet `R2` values to unbiased `r^2`. |
| `-` | `--regression-snps` | Optional SNP list defining the regression SNP set used to compute `w_ld`. This is a new concept and does not replace legacy `--print-snps`. |
| `-` | `--per-chr-output` | Write per-chromosome outputs instead of one aggregated result set. |
| `-` | `--log-level` | Set refactored logging verbosity. |

## Heritability Estimates

Refactor status: planned as `ldsc regress h2 ...`, but no refactored CLI flags are implemented yet.

| old flag | new flag | meaning |
| --- | --- | --- |
| `--h2` | `Not yet implemented` | Run single-trait LD Score regression from one munged `.sumstats` input. |
| `--out` | `Not yet implemented` | Output prefix for heritability logs and optional result files. |
| `--ref-ld` | `Not yet implemented` | Reference LD-score predictor files for regression. |
| `--ref-ld-chr` | `Not yet implemented` | Chromosome-split reference LD-score predictor files. |
| `--w-ld` | `Not yet implemented` | Regression-weight LD-score file set. |
| `--w-ld-chr` | `Not yet implemented` | Chromosome-split regression-weight LD-score file set. |
| `--no-intercept` | `Not yet implemented` | Constrain the heritability intercept to 1. |
| `--intercept-h2` | `Not yet implemented` | Provide an explicit constrained intercept for `h2`. |
| `--M` | `Not yet implemented` | Override annotation SNP counts instead of reading `.l2.M` / `.l2.M_5_50`. |
| `--two-step` | `Not yet implemented` | Use the two-step estimator with the supplied chi-square cutoff. |
| `--chisq-max` | `Not yet implemented` | Maximum chi-square threshold used to filter large statistics. |
| `--print-cov` | `Not yet implemented` | Write covariance matrix of regression estimates. |
| `--print-delete-vals` | `Not yet implemented` | Write block-jackknife delete values. |
| `--invert-anyway` | `Not yet implemented` | Force inversion of ill-conditioned LD-score matrices. |
| `--n-blocks` | `Not yet implemented` | Number of jackknife blocks. |
| `--not-M-5-50` | `Not yet implemented` | Use `.l2.M` rather than `.l2.M_5_50` for regression counts. |
| `--samp-prev` | `Not yet implemented` | Sample prevalence for liability-scale conversion of binary traits. |
| `--pop-prev` | `Not yet implemented` | Population prevalence for liability-scale conversion of binary traits. |

## Partitioned LDSC

Refactor status: planned as `ldsc regress partitioned-h2 ...`, but no refactored CLI flags are implemented yet.

| old flag | new flag | meaning |
| --- | --- | --- |
| `--h2` | `Not yet implemented` | Run partitioned heritability through the same legacy `--h2` entrypoint when multiple LD-score annotations are supplied. |
| `--out` | `Not yet implemented` | Output prefix for partitioned LDSC logs and result tables. |
| `--ref-ld` | `Not yet implemented` | Reference LD-score predictor files containing multiple annotation columns. |
| `--ref-ld-chr` | `Not yet implemented` | Chromosome-split partitioned LD-score predictor files. |
| `--w-ld` | `Not yet implemented` | Regression-weight LD-score file set. |
| `--w-ld-chr` | `Not yet implemented` | Chromosome-split regression-weight LD-score file set. |
| `--overlap-annot` | `Not yet implemented` | Declare that annotation categories overlap, which changes output interpretation. |
| `--print-coefficients` | `Not yet implemented` | Print coefficients in addition to heritability summaries when annotations overlap. |
| `--frqfile` | `Not yet implemented` | Frequency file used to prune to common SNPs for overlap-aware summaries unless `--not-M-5-50` is set. |
| `--frqfile-chr` | `Not yet implemented` | Chromosome-split version of `--frqfile`. |
| `--M` | `Not yet implemented` | Override annotation SNP counts instead of reading `.l2.M` / `.l2.M_5_50`. |
| `--two-step` | `Not yet implemented` | Two-step estimator cutoff. |
| `--chisq-max` | `Not yet implemented` | Maximum chi-square threshold used to filter large statistics. |
| `--no-intercept` | `Not yet implemented` | Constrain the partitioned `h2` intercept to 1. |
| `--intercept-h2` | `Not yet implemented` | Provide an explicit constrained intercept. |
| `--print-cov` | `Not yet implemented` | Write covariance matrix of regression estimates. |
| `--print-delete-vals` | `Not yet implemented` | Write block-jackknife delete values, including partitioned delete values. |
| `--invert-anyway` | `Not yet implemented` | Force inversion of ill-conditioned LD-score matrices. |
| `--n-blocks` | `Not yet implemented` | Number of jackknife blocks. |
| `--not-M-5-50` | `Not yet implemented` | Use `.l2.M` rather than `.l2.M_5_50` for regression counts. |
| `--samp-prev` | `Not yet implemented` | Sample prevalence for liability-scale conversion of binary traits. |
| `--pop-prev` | `Not yet implemented` | Population prevalence for liability-scale conversion of binary traits. |

## Cross-Trait Regression

Refactor status: planned as `ldsc regress rg ...`, but no refactored CLI flags are implemented yet.

| old flag | new flag | meaning |
| --- | --- | --- |
| `--rg` | `Not yet implemented` | Run genetic correlation from two or more munged `.sumstats` files. |
| `--out` | `Not yet implemented` | Output prefix for genetic-correlation logs and optional artifact files. |
| `--ref-ld` | `Not yet implemented` | Reference LD-score predictor files shared by all traits in the run. |
| `--ref-ld-chr` | `Not yet implemented` | Chromosome-split reference LD-score predictor files. |
| `--w-ld` | `Not yet implemented` | Regression-weight LD-score file set. |
| `--w-ld-chr` | `Not yet implemented` | Chromosome-split regression-weight LD-score file set. |
| `--no-intercept` | `Not yet implemented` | Constrain trait heritability intercepts to 1 and genetic-covariance intercepts to 0. |
| `--intercept-h2` | `Not yet implemented` | Trait-specific constrained heritability intercepts. |
| `--intercept-gencov` | `Not yet implemented` | Trait-pair constrained genetic-covariance intercepts. |
| `--M` | `Not yet implemented` | Override annotation SNP counts instead of reading `.l2.M` / `.l2.M_5_50`. |
| `--two-step` | `Not yet implemented` | Two-step estimator cutoff. |
| `--chisq-max` | `Not yet implemented` | Maximum chi-square threshold used to filter large statistics. |
| `--print-cov` | `Not yet implemented` | Write covariance matrices for `rg` estimates. |
| `--print-delete-vals` | `Not yet implemented` | Write block-jackknife delete values for `rg` estimates. |
| `--invert-anyway` | `Not yet implemented` | Force inversion of ill-conditioned LD-score matrices. |
| `--n-blocks` | `Not yet implemented` | Number of jackknife blocks. |
| `--not-M-5-50` | `Not yet implemented` | Use `.l2.M` rather than `.l2.M_5_50` for regression counts. |
| `--return-silly-things` | `Not yet implemented` | Force return of otherwise rejected extreme / nonsensical genetic-correlation estimates. |
| `--no-check-alleles` | `Not yet implemented` | Skip cross-trait allele-matching checks. |
| `--samp-prev` | `Not yet implemented` | Trait-specific sample prevalences for liability-scale reporting. |
| `--pop-prev` | `Not yet implemented` | Trait-specific population prevalences for liability-scale reporting. |

## Cell Type Specific LDSC

Refactor status: not yet implemented as a concrete refactored CLI. The design docs describe it as a later specialization of the regression stack, but they do not currently define replacement flags.

| old flag | new flag | meaning |
| --- | --- | --- |
| `--h2-cts` | `Not yet implemented` | Run cell-type-specific LDSC from one munged `.sumstats` file. |
| `--out` | `Not yet implemented` | Output prefix for `.cell_type_results.txt` and logs. |
| `--ref-ld-chr` | `Not yet implemented` | Baseline chromosome-split LD-score predictors used in every cell-type regression. |
| `--ref-ld-chr-cts` | `Not yet implemented` | File listing chromosome-split LD-score prefixes for each cell-type-specific predictor set. |
| `--w-ld` | `Not yet implemented` | Regression-weight LD-score file set. |
| `--w-ld-chr` | `Not yet implemented` | Chromosome-split regression-weight LD-score file set. |
| `--no-intercept` | `Not yet implemented` | Constrain the intercept to 1. |
| `--intercept-h2` | `Not yet implemented` | Provide an explicit constrained intercept. |
| `--M` | `Not yet implemented` | Override annotation SNP counts instead of reading `.l2.M` / `.l2.M_5_50`. |
| `--chisq-max` | `Not yet implemented` | Maximum chi-square threshold used to filter large statistics. |
| `--print-all-cts` | `Not yet implemented` | Print all coefficients for multi-column cell-type LD-score inputs rather than only the first coefficient. |
| `--invert-anyway` | `Not yet implemented` | Force inversion of ill-conditioned LD-score matrices. |
| `--n-blocks` | `Not yet implemented` | Number of jackknife blocks. |
| `--not-M-5-50` | `Not yet implemented` | Use `.l2.M` rather than `.l2.M_5_50` for regression counts. |
