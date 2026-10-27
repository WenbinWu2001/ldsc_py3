# LDSC2 to Current CLI Flag Map

Last updated on: 2026-10-27

This document maps changed LDSC2 flags to the supported `ldsc` commands and records the flag-cleanup decisions.
LDSC2 source links point to `bulik/ldsc` at commit `aa33296`. “Removed” means the parser rejects the old option; there are no compatibility aliases. A retired option whose behavior is now automatic has that behavior described explicitly below. The [IO argument inventory](cli.md) lists the complete current surface, and [munge-sumstats](munge-sumstats.md) describes its input and artifact contracts.

## Changes in the flag cleanup

Some replaced spellings came from earlier LDSC3 revisions rather than LDSC2. They are distinguished here so that the mapping does not invent a legacy equivalent.

| LDSC2 flag | Earlier LDSC3 spelling | Current command/flag | Current behavior |
| --- | --- | --- | --- |
| No equivalent LDSC2 flag | `--gene-list-resolution-policy {strict,resolved-only}` | `annotate`, `ldscore`, and `quantile-h2`: `--allow-unresolved-genes` | Value-free opt-in flag. Omission retains strict resolution; supplying it selects the existing `resolved-only` policy and omission diagnostics. Remove the entire `--gene-list-resolution-policy strict` option/value pair from existing commands; replace `--gene-list-resolution-policy resolved-only` with the new flag. The old option is rejected without an alias. Python policy arguments and metadata values remain unchanged. |
| `munge_sumstats.py --keep-maf` | `munge-sumstats --keep-maf` | Removed; no replacement toggle | Preserve a recognized/selected input frequency column as `FRQ` whenever present; omit `FRQ` when absent. Preserve its values without converting them to MAF, including frequencies above 0.5. `--maf-min` still filters on folded MAF. There is no `--keep-frequency` option. |
| `munge_sumstats.py --daner` | `munge-sumstats --format daner-old` | `ldsc munge-sumstats --input-format daner-old` | Read constant case/control counts from `FRQ_A_<Ncas>` and `FRQ_U_<Ncon>` headers. Default `--input-format auto` can detect this profile. |
| `munge_sumstats.py --daner-n` | `munge-sumstats --format daner-new` | `ldsc munge-sumstats --input-format daner-new` | Retain new-DANER selection, using the same shared preparation, aliases, optional fields, validation, and reader as automatically detected new-DANER input. `Nca/Nco`, `NCAS/NCON`, and case variants use the central alias registry. `FRQ_U_*` is optional. |
| No generic LDSC2 format selector | `munge-sumstats --format` | `ldsc munge-sumstats --input-format` | Values remain `auto`, `plain`, `daner-old`, and `daner-new`; default remains `auto`. This selects the input profile, independently of `--output-format`. Earlier LDSC3 booleans `--daner-old` and `--daner-new` are also unsupported. |
| No equivalent LDSC2 flag | `ldscore --exclude-regions` | `ldsc ldscore --regr-snps-exclude-regions` | The hidden alias is removed. Default remains `mhc-and-centromeres`; choices remain `none`, `mhc`, `centromeres`, and `mhc-and-centromeres`. This changes regression/output rows, not the LD-score contributor or count universe. |
| No LDSC2 gene-index command | `build-gene-ldscore-index --exclude-regions` | `ldsc build-gene-ldscore-index --regr-snps-exclude-regions` | Same alias removal and unchanged default region policy. |
| `ldsc.py --no-intercept` with `--h2` | `h2 --no-intercept` | `ldsc h2 --intercept-h2 1` | Removed shortcut; explicitly fix the h2 intercept to 1. Omit the replacement to estimate it using the existing defaults. |
| `ldsc.py --no-intercept` in partitioned `--h2` or `--h2-cts` runs | `partitioned-h2 --no-intercept` | `ldsc partitioned-h2 --intercept-h2 1` | Removed shortcut; explicitly fix the h2 intercept to 1 in each fitted model. Omit the replacement to estimate it using the existing defaults. |
| `ldsc.py --no-intercept` with `--rg` | `rg --no-intercept` | `ldsc rg --intercept-h2 1 --intercept-gencov 0` | Removed shortcut; fix both traits' h2 intercepts to 1 in every pair and each genetic-covariance intercept to 0. Omit both replacements to estimate the intercepts using the existing defaults. |

Direct/indexed `ldscore`, `build-gene-ldscore-index`, and batch regression now use `--workers`.
The former `--threads` spelling is retired without an alias.
The default is one managed child; every worker uses 1 BLAS thread in serial and parallel execution.
See [parallel execution](parallel-execution.md) for joblib CPU detection, validation, and task caps.

Python counterparts also retain `LDScoreConfig.workers`, `run_ldscore(..., workers=...)`, and `GeneLDScoreIndexBuildConfig.workers`. `MungeConfig.keep_maf` is removed because preserving frequency is unconditional. The Python input-profile field remains `MungeConfig.sumstats_format`; its CLI spelling is now `--input-format`.

Frequency preservation applies to the in-memory result and both output formats. Parquet retains numeric precision. Gzip TSV preserves the parsed frequency without the three-decimal rounding used for other floating columns; missing fields are written as `NA` so whitespace readers do not shift column positions. An explicitly ignored frequency column is not selected and therefore is absent from output. Multiple competing frequency columns still follow existing column-selection/ambiguity rules; this change does not choose among them silently.

Sources: current [`sumstats_munger.build_parser`, `_apply_raw_sumstats_inference`, `_write_sumstats_tsv_gz`](../../src/ldsc/sumstats_munger.py), [`_sumstats_input.prepare_munge_input`](../../src/ldsc/_sumstats_input.py), [`_kernel.sumstats_munger.parse_dat` / `filter_frq`](../../src/ldsc/_kernel/sumstats_munger.py), [`ldscore_calculator.build_parser` / `_resolve_worker_count`](../../src/ldsc/ldscore_calculator.py), and [`gene_ldscore_index.build_parser`](../../src/ldsc/gene_ldscore_index.py). LDSC2 declarations and behavior: [`munge_sumstats.py`, parser and `parse_dat`](https://github.com/bulik/ldsc/blob/aa33296abac9569a6422ee6ba7eb4b902422cc74/munge_sumstats.py#L452); [`ldsc.py`, parser](https://github.com/bulik/ldsc/blob/aa33296abac9569a6422ee6ba7eb4b902422cc74/ldsc.py#L423).

## Regression intercept consolidation

The LDSC2 `--no-intercept` flag disabled intercept estimation by choosing standard fixed values. Its name did not mean an h2 intercept of zero. The current CLI rejects that flag, with no hidden alias. Use the following explicit replacements:

| Current command | Flags replacing LDSC2 `--no-intercept` |
| --- | --- |
| `h2` | `--intercept-h2 1` |
| `partitioned-h2` | `--intercept-h2 1` |
| `rg` | `--intercept-h2 1 --intercept-gencov 0` |

These disable intercept estimation by fixing h2 intercepts to 1 and, for `rg`, covariance intercepts to 0. In `rg`, one scalar h2 value applies to both traits in every pair, and one scalar covariance value applies to every pair. Omitting these flags continues to estimate the intercepts with the existing defaults; fixed intercepts have not become automatic. Other numeric fixed values remain accepted through the same flags. `quantile-h2` reads a saved fitted model and has no intercept-estimation flags.

Single-annotation `rg --intercept-gencov 0` with a free h2 intercept continues to fail: automatic two-step estimation selects cutoff 30, which conflicts with the fixed covariance intercept. Explicit `--two-step-cutoff` remains incompatible with fixed intercepts and multi-annotation fits. This consolidation adds no way to disable automatic two-step estimation.

The public Python `RegressionConfig.use_intercept=False` option remains available and unchanged. The CLI now expresses its standard fixed-intercept behavior through `--intercept-h2` and `--intercept-gencov` only.

Sources: LDSC2 [`ldsc.py`, intercept flags](https://github.com/bulik/ldsc/blob/aa33296abac9569a6422ee6ba7eb4b902422cc74/ldsc.py#L528) and [`ldscore/sumstats.py`, `estimate_h2` / `estimate_rg`](https://github.com/bulik/ldsc/blob/aa33296abac9569a6422ee6ba7eb4b902422cc74/ldscore/sumstats.py#L321); current [`regression_runner._add_common_regression_arguments`, `_runner_from_args`, `_fit_h2_dataset`, and `_fit_rg_dataset`](../../src/ldsc/regression_runner.py), [`RegressionConfig`](../../src/ldsc/config.py), and [regression configuration, “Requesting fixed intercepts”](regression-options.md#15-requesting-fixed-intercepts).

## Legacy behavior explicitly retained

| LDSC2 flag | Current flag on `ldsc munge-sumstats` | Behavior |
| --- | --- | --- |
| `--N` | `--N` | Constant fallback. Input per-variant N or paired case/control columns take precedence. With no per-variant strategy, supplied `--N` takes precedence over the case/control constants. |
| `--N-cas`, `--N-con` | Same names | Both remain accepted. When used as the fallback, N is their sum. These are distinct from the paired per-variant column strategy. |
| `--N-col`, `--N-cas-col`, `--N-con-col` | Same names | Keep explicit per-variant strategy selection. Direct N and paired case/control selections remain mutually exclusive under the current schema contract. |
| `--n-min` | `--n-min` | Default is the 90th percentile of per-variant N divided by **1.5**. Omitted and zero values both invoke that default. Keep rows with `N >= n_min`. Constant N is assigned after the filter and bypasses it, as in LDSC2. |
| `--nstudy`, `--nstudy-min` | Same names | Study-count filtering applies when per-variant N is absent. Omitted/zero thresholds use the maximum study count. A per-variant N strategy takes precedence over this study-count filter. |
| `--chunksize` | `--chunksize` | Spelling retained. Current default is 1,000,000 input rows per chunk; LDSC2's parser default is 5,000,000. This cleanup does not change the current default. It is separate from the old LD-score `--chunk-size` flag below. |
| `--maf-min` | `--maf-min` | Compare folded `min(FRQ, 1-FRQ)` with the inclusive threshold. The LDSC2 default 0.01 applies only when a frequency column is selected; without one the default filter is skipped, and an explicit value is an error (LDSC2 skipped it silently). This does not overwrite the stored frequency. |

The existing per-variant case/control normalization is also unchanged: define total T and case fraction p per row, then `N_i = T_i * p_i / mean(p among rows with maximum T)`. It must not be replaced with simply summing the paired per-variant columns or with an effective-N formula as part of flag cleanup. See [“Sample-Size Column Selection”](munge-sumstats.md#sample-size-column-selection).

Legacy numerical source: [`munge_sumstats.py:321`, `process_n`](https://github.com/bulik/ldsc/blob/aa33296abac9569a6422ee6ba7eb4b902422cc74/munge_sumstats.py#L321); current numerical source: [`_kernel.sumstats_munger.process_n`](../../src/ldsc/_kernel/sumstats_munger.py). The legacy help text itself incorrectly describes the divisor as 2; the implementation is the authority for the retained 1.5 rule. No new default or override policy is introduced here.

## Legacy behavior intentionally changed

### Exact PLINK LD windows

LDSC2 computed PLINK LD scores in chunks of `--chunk-size` SNPs (default 50) and widened each SNP's window to a multiple of the chunk, so the chunk size changed the LD scores.
LDSC2's own code comment notes that these windows differ from the definition in the paper.
The current PLINK backend keeps the chunked computation but drops every pair outside the requested window, so `--snp-batch-size` in `ldsc ldscore` and `ldsc build-gene-ldscore-index` changes only speed and memory.
Its LD scores equal LDSC2 run with `--chunk-size 1` and use the same window rule as the parquet-R2 backend, which always used exact windows.
The two backends still differ slightly, because parquet R² is stored as int16 and a parquet `--ld-wind-snps` window counts panel rows rather than retained SNPs.
A cM window also differs when the panel's `CM`, interpolated from a genetic map, differs from the `.bim` `CM` the PLINK run uses.
Gene LD-score indexes built with the earlier rounded windows (`format_version` 2) are rejected and must be rebuilt.
Published LDSC2 LD scores remain usable unchanged through `ldsc convert-ldsc2-ldscores`.

Measurement (2026-10-27; former LDSC2-ported kernel at `c67b40e` versus the exact kernel introduced with this section), on 1000G.EUR.QC chromosome 22 (489 individuals, 141,123 SNPs), `--ld-wind-cm 1`, 17,489 HapMap3 output SNPs:

- The former kernel at batch 50 reproduced the published Phase 3 baseline `baseL2` within its three printed decimals.
- Relative to the exact window, batch 50 raised the mean base LD score by 0.05% (largest per-SNP difference 6.7; 2.2% of SNPs differed by more than 0.5), and the former default batch 128 raised it by 0.13% (largest difference 7.6; 8.5% of SNPs).
- Chromosome-22 SNP heritability for the PGC3 schizophrenia summary statistics (wave 3, European) was 0.008406 at batch 50, 0.008413 at batch 128, and 0.008411 with exact windows (standard error 0.00114).
- The exact kernel took 18 s at batch 128 and gave the same written values as the former kernel at batch 1, which took 264 s.

Limits: one chromosome, panel, trait, and window were measured.
Smaller `--ld-wind-snps` windows and sparse panels change more; on the 32-SNP test fixture, `--ld-wind-snps 3` at batch 128 formerly covered the whole chromosome.

Sources: LDSC2 [`ldscore/ldscore.py:167-168`, `__corSumVarBlocks__`](https://github.com/bulik/ldsc/blob/aa33296abac9569a6422ee6ba7eb4b902422cc74/ldscore/ldscore.py#L167-L168); current [`_kernel/plink_bed.py`, `__corSumVarBlocks__` and `_zero_pairs_outside_window`](../../src/ldsc/_kernel/plink_bed.py), tested in [`test_ldscore_projection.py`](../../tests/test_ldscore_projection.py) and [`test_plink_exact_windows.py`](../../tests/test_plink_exact_windows.py).

## Established command and input renames

These mappings describe current workflows and are not additional changes made by this cleanup. Where workflow semantics differ, a replacement is described as a workflow rather than a one-to-one alias.

| LDSC2 command/flag | Current entry point | Meaning or required adaptation |
| --- | --- | --- |
| `munge_sumstats.py --sumstats FILE` | `ldsc munge-sumstats --raw-sumstats-file FILE` | Raw input; keep it distinct from regression's curated `--sumstats-file`. |
| `--out PREFIX` | `--output-dir DIRECTORY` on materializing commands | Directory of fixed artifacts, not an output prefix. `plot` and `convert-h2-scale` derive fixed destinations from the input result directory. |
| `ldsc.py --l2` | `ldsc ldscore` | LD-score computation is a subcommand. |
| `ldsc.py --bfile PREFIX` | `ldsc ldscore --plink-prefix PREFIX` | PLINK BED/BIM/FAM input; current prefix resolution also handles chromosome suites. |
| `ldsc.py --extract FILE` | `ldsc ldscore --ref-panel-snps-file FILE` | Select the reference contributor universe using the current headered identity-file contract. |
| `ldsc.py --keep FILE` | `ldsc ldscore --keep-indivs-file FILE` | Restrict PLINK individuals; invalid with precomputed R2. |
| `ldsc.py --annot SOURCES` | `--baseline-annot-sources`, optionally `--query-annot-sources` | Explicitly separate baseline annotations from focal queries; current annotation inputs include identity metadata. |
| `ldsc.py --maf VALUE` | `ldsc ldscore --maf-min VALUE` | Reference-panel MAF filtering; `--common-maf-min` separately controls common-SNP count vectors. |
| `ldsc.py --chunk-size N` | `ldsc ldscore --snp-batch-size N` | PLINK genotype computation batch size, current default 512. Unlike `--chunk-size`, it does not widen LD windows; see [Exact PLINK LD windows](#exact-plink-ld-windows). This is not munging's `--chunksize`. |
| `ldsc.py --print-snps FILE` | `--regr-snps-file FILE` with a deliberate workflow adaptation | The current set determines both persisted LD-score rows and regression-weight contributors. It is not merely the old write-only filter; default HM3 selection and region exclusions also apply. |
| `ldsc.py --h2 FILE` | `ldsc h2 --sumstats-file FILE` or `ldsc partitioned-h2 --sumstats-file FILE` | Select the intended workflow; both consume a canonical LD-score directory. |
| `ldsc.py --h2-cts FILE`, `--ref-ld-chr-cts` | `ldsc partitioned-h2 --sumstats-file FILE --ldscore-dir DIRECTORY` | Build query LD scores first and use the canonical artifact; old cell-type manifests are not accepted as regression inputs. |
| `ldsc.py --rg FILES` | `ldsc rg --sumstats-sources FILES` | Current default is all pairs; use `--anchor-trait` for an anchor-versus-rest run. |
| `--ref-ld`, `--ref-ld-chr`, `--w-ld`, `--w-ld-chr` | `--ldscore-dir DIRECTORY` on regression commands | Predictors, weights, and count metadata come from one canonical artifact. Convert legacy fragments with `convert-ldsc2-ldscores` first. |
| `--two-step CUTOFF` | `--two-step-cutoff CUTOFF` | Estimator cutoff; current inclusive threshold and model/intercept restrictions still apply. |
| `--not-M-5-50` | `--count-kind all` | Opt into all-SNP counts. Default remains `--count-kind common`. |
| `munge_sumstats.py --no-alleles` | `--snp-identifier rsid` or `--snp-identifier chr_pos` | Choose allele-blind identity explicitly; default current identity is `chr_pos_allele_aware`. |
| `munge_sumstats.py --merge-alleles FILE` | No exact alias; packaged HM3 restriction is the default, `--sumstats-snps-file` replaces it, and `--no-snp-restriction` disables it | Current identity matching and downstream allele alignment have separate roles. A keep-list does not imply allele rewriting or liftover. |

Sources: LDSC2 [`ldsc.py`, parser](https://github.com/bulik/ldsc/blob/aa33296abac9569a6422ee6ba7eb4b902422cc74/ldsc.py#L423) and [`munge_sumstats.py`, parser](https://github.com/bulik/ldsc/blob/aa33296abac9569a6422ee6ba7eb4b902422cc74/munge_sumstats.py#L449); current [`cli.build_parser`](../../src/ldsc/cli.py), [`ldscore_calculator.build_parser`](../../src/ldsc/ldscore_calculator.py), [`regression_runner` argument builders](../../src/ldsc/regression_runner.py), and the [IO inventory, command sections](cli.md). Unlisted unsupported legacy options should not be assumed to have an automatic replacement.

## LDSC2 options with no current equivalent

These LDSC2 options have no counterpart in any current subcommand's `--help`; the parser rejects them.
Line references are to the LDSC2 repository.

| LDSC2 option | Current situation |
| --- | --- |
| `ldsc.py --per-allele`, `--pq-exp` (`ldsc.py:477`, `:481`) | No per-allele or $p(1-p)$-exponent LD-score scaling exists, and there is no workaround. |
| `ldsc.py --cts-bin`, `--cts-breaks`, `--cts-names`, `--no-print-annot` (`ldsc.py:460-485`) | No on-the-fly binning of a SNP-to-value file. Continuous-valued annotation files are accepted as query annotations, and `quantile-h2` summarizes a continuous target after fitting. Binning by hand remains possible but is lossy. |
| `ldsc.py --M` (`ldsc.py:538`) | No manual SNP-count override. Counts come from the LD-score directory; `--count-kind` selects common or all-SNP counts. |
| `ldsc.py --print-cov` (`ldsc.py:549`) | No covariance-matrix output. |
| `ldsc.py --print-delete-vals` (`ldsc.py:552`) | No flag; partitioned fits always write `diagnostics/coefficient_delete_values.parquet` ([partitioned-h2 results](results/partitioned-h2.md)). |
| `ldsc.py --invert-anyway` (`ldsc.py:564`) | No forced inversion of ill-conditioned matrices. |
| `ldsc.py --return-silly-things` (`ldsc.py:570`) | No flag; the default matches LDSC2's default behavior. |
| `ldsc.py --pickle`, `--print-all-cts` (`ldsc.py:560`, `:546`) | No equivalent output options. |
| `make_annot.py --nomerge` (`make_annot.py:43`) | Annotation projection is membership-only; there is no mode proportional to the number of overlapping intervals. |

## Naming audit boundary

Current README, documentation under `docs/`, troubleshooting guidance, tutorials, public configuration docstrings, code comments, and tests use the current names. The shared `--workers` spelling remains current for both LD-score computation and gene-index building. References to the abandoned `--threads` spelling are limited to this decision record, rejection tests, and git history.
