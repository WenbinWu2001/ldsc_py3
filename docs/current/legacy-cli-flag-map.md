# LDSC2 to Current CLI Flag Map

Last updated on: 2026-09-11

This document maps changed LDSC2 flags to the supported `ldsc` commands and records the flag-cleanup decisions. “Removed” means the parser rejects the old option; there are no compatibility aliases. A retired option whose behavior is now automatic has that behavior described explicitly below. The [IO argument inventory](io-argument-inventory.md) lists the complete current surface, and [munge-sumstats](munge-sumstats.md) describes its input and artifact contracts.

## Changes in the flag cleanup

Some replaced spellings came from earlier LDSC3 revisions rather than LDSC2. They are distinguished here so that the mapping does not invent a legacy equivalent.

| LDSC2 flag | Earlier LDSC3 spelling | Current command/flag | Current behavior |
| --- | --- | --- | --- |
| `munge_sumstats.py --keep-maf` | `munge-sumstats --keep-maf` | Removed; no replacement toggle | Preserve a recognized/selected input frequency column as `FRQ` whenever present; omit `FRQ` when absent. Preserve its values without converting them to MAF, including frequencies above 0.5. `--maf-min` still filters on folded MAF. There is no `--keep-frequency` option. |
| `munge_sumstats.py --daner` | `munge-sumstats --format daner-old` | `ldsc munge-sumstats --input-format daner-old` | Read constant case/control counts from `FRQ_A_<Ncas>` and `FRQ_U_<Ncon>` headers. Default `--input-format auto` can detect this profile. |
| `munge_sumstats.py --daner-n` | `munge-sumstats --format daner-new` | `ldsc munge-sumstats --input-format daner-new` | Retain new-DANER selection, using the same shared preparation, aliases, optional fields, validation, and reader as automatically detected new-DANER input. `Nca/Nco`, `NCAS/NCON`, and case variants use the central alias registry. `FRQ_U_*` is optional. |
| No generic LDSC2 format selector | `munge-sumstats --format` | `ldsc munge-sumstats --input-format` | Values remain `auto`, `plain`, `daner-old`, and `daner-new`; default remains `auto`. This selects the input profile, independently of `--output-format`. Earlier LDSC3 booleans `--daner-old` and `--daner-new` are also unsupported. |
| No equivalent LDSC2 flag | `ldscore --exclude-regions` | `ldsc ldscore --regr-snps-exclude-regions` | The hidden alias is removed. Default remains `mhc-and-centromeres`; choices remain `none`, `mhc`, `centromeres`, and `mhc-and-centromeres`. This changes regression/output rows, not the LD-score contributor or count universe. |
| No LDSC2 gene-index command | `build-gene-ldscore-index --exclude-regions` | `ldsc build-gene-ldscore-index --regr-snps-exclude-regions` | Same alias removal and unchanged default region policy. |
| `ldsc.py --no-intercept` with `--h2` | `h2 --no-intercept` | `ldsc h2 --intercept-h2 1` | Removed shortcut; explicitly fix the h2 intercept to 1. Omit the replacement to estimate it using the existing defaults. |
| `ldsc.py --no-intercept` in partitioned `--h2` or `--h2-cts` runs | `partitioned-h2 --no-intercept` | `ldsc partitioned-h2 --intercept-h2 1` | Removed shortcut; explicitly fix the h2 intercept to 1 in each fitted model. Omit the replacement to estimate it using the existing defaults. |
| `ldsc.py --no-intercept` with `--rg` | `rg --no-intercept` | `ldsc rg --intercept-h2 1 --intercept-gencov 0` | Removed shortcut; fix both traits' h2 intercepts to 1 in every pair and each genetic-covariance intercept to 0. Omit both replacements to estimate the intercepts using the existing defaults. |

Both `ldscore` and `build-gene-ldscore-index` retain `--threads`; the proposed `--workers` rename was reverted. In direct `ldscore`, it controls chromosome worker processes. In `build-gene-ldscore-index`, it controls a chromosome thread pool. Both default to `1` (sequential); positive N requests N concurrent chromosomes, `-1` uses available cores, and `-2` leaves one core free, with effective concurrency capped at the chromosome count. Computation and concurrency policies are unchanged. BLAS thread limits are separate internal controls.

Python counterparts also retain `LDScoreConfig.threads`, `run_ldscore(..., threads=...)`, and `GeneLDScoreIndexBuildConfig.threads`. `MungeConfig.keep_maf` is removed because preserving frequency is unconditional. The Python input-profile field remains `MungeConfig.sumstats_format`; its CLI spelling is now `--input-format`.

Frequency preservation applies to the in-memory result and both output formats. Parquet retains numeric precision. Gzip TSV preserves the parsed frequency without the three-decimal rounding used for other floating columns; missing fields are written as `NA` so whitespace readers do not shift column positions. An explicitly ignored frequency column is not selected and therefore is absent from output. Multiple competing frequency columns still follow existing column-selection/ambiguity rules; this change does not choose among them silently.

Sources: current [`sumstats_munger.build_parser`, `_apply_raw_sumstats_inference`, `_write_sumstats_tsv_gz`](../../src/ldsc/sumstats_munger.py), [`_sumstats_input.prepare_munge_input`](../../src/ldsc/_sumstats_input.py), [`_kernel.sumstats_munger.parse_dat` / `filter_frq`](../../src/ldsc/_kernel/sumstats_munger.py), [`ldscore_calculator.build_parser` / `_resolve_worker_count`](../../src/ldsc/ldscore_calculator.py), and [`gene_ldscore_index.build_parser`](../../src/ldsc/gene_ldscore_index.py). LDSC2 declarations and behavior: [`munge_sumstats.py`, parser and `parse_dat`](</Users/wenbinwu/Documents_local/Research/SullivanLab/LDSC/repos/ldsc_py2_Bulik_workspace/ldsc_py2_Bulik/munge_sumstats.py:459>); [`ldsc.py`, parser](</Users/wenbinwu/Documents_local/Research/SullivanLab/LDSC/repos/ldsc_py2_Bulik_workspace/ldsc_py2_Bulik/ldsc.py:424>).

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

Sources: LDSC2 [`ldsc.py`, intercept flags](</Users/wenbinwu/Documents_local/Research/SullivanLab/LDSC/repos/ldsc_py2_Bulik_workspace/ldsc_py2_Bulik/ldsc.py:529>) and [`ldscore/sumstats.py`, `estimate_h2` / `estimate_rg`](</Users/wenbinwu/Documents_local/Research/SullivanLab/LDSC/repos/ldsc_py2_Bulik_workspace/ldsc_py2_Bulik/ldscore/sumstats.py:321>); current [`regression_runner._add_common_regression_arguments`, `_runner_from_args`, `_fit_h2_dataset`, and `_fit_rg_dataset`](../../src/ldsc/regression_runner.py), [`RegressionConfig`](../../src/ldsc/config.py), and [regression configuration, “Requesting fixed intercepts”](regression-configuration.md#15-requesting-fixed-intercepts).

## Legacy behavior explicitly retained

| LDSC2 flag | Current flag on `ldsc munge-sumstats` | Behavior |
| --- | --- | --- |
| `--N` | `--N` | Constant fallback. Input per-variant N or paired case/control columns take precedence. With no per-variant strategy, supplied `--N` takes precedence over the case/control constants. |
| `--N-cas`, `--N-con` | Same names | Both remain accepted. When used as the fallback, N is their sum. These are distinct from the paired per-variant column strategy. |
| `--N-col`, `--N-cas-col`, `--N-con-col` | Same names | Keep explicit per-variant strategy selection. Direct N and paired case/control selections remain mutually exclusive under the current schema contract. |
| `--n-min` | `--n-min` | Default is the 90th percentile of per-variant N divided by **1.5**. Omitted and zero values both invoke that default. Keep rows with `N >= n_min`. Constant N is assigned after the filter and bypasses it, as in LDSC2. |
| `--nstudy`, `--nstudy-min` | Same names | Study-count filtering applies when per-variant N is absent. Omitted/zero thresholds use the maximum study count. A per-variant N strategy takes precedence over this study-count filter. |
| `--chunksize` | `--chunksize` | Spelling retained. Current default is 1,000,000 input rows per chunk; LDSC2's parser default is 5,000,000. This cleanup does not change the current default. It is separate from the old LD-score `--chunk-size` flag below. |
| `--maf-min` | `--maf-min` | Current QC behavior is unchanged: compare folded `min(FRQ, 1-FRQ)` with the inclusive threshold, default 0.01. This does not overwrite the stored frequency. |

The existing per-variant case/control normalization is also unchanged: define total T and case fraction p per row, then `N_i = T_i * p_i / mean(p among rows with maximum T)`. It must not be replaced with simply summing the paired per-variant columns or with an effective-N formula as part of flag cleanup. See [“Sample-Size Column Selection”](munge-sumstats.md#sample-size-column-selection).

Legacy numerical source: [`munge_sumstats.py:327`, `process_n`](</Users/wenbinwu/Documents_local/Research/SullivanLab/LDSC/repos/ldsc_py2_Bulik_workspace/ldsc_py2_Bulik/munge_sumstats.py:327>); current numerical source: [`_kernel.sumstats_munger.process_n`](../../src/ldsc/_kernel/sumstats_munger.py). The legacy help text itself incorrectly describes the divisor as 2; the implementation is the authority for the retained 1.5 rule. No new default or override policy is introduced here.

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
| `ldsc.py --chunk-size N` | `ldsc ldscore --snp-batch-size N` | PLINK genotype computation batch size, current default 128. This is not munging's `--chunksize`. |
| `ldsc.py --print-snps FILE` | `--regr-snps-file FILE` with a deliberate workflow adaptation | The current set determines both persisted LD-score rows and regression-weight contributors. It is not merely the old write-only filter; default HM3 selection and region exclusions also apply. |
| `ldsc.py --h2 FILE` | `ldsc h2 --sumstats-file FILE` or `ldsc partitioned-h2 --sumstats-file FILE` | Select the intended workflow; both consume a canonical LD-score directory. |
| `ldsc.py --h2-cts FILE`, `--ref-ld-chr-cts` | `ldsc partitioned-h2 --sumstats-file FILE --ldscore-dir DIRECTORY` | Build query LD scores first and use the canonical artifact; old cell-type manifests are not accepted as regression inputs. |
| `ldsc.py --rg FILES` | `ldsc rg --sumstats-sources FILES` | Current default is all pairs; use `--anchor-trait` for an anchor-versus-rest run. |
| `--ref-ld`, `--ref-ld-chr`, `--w-ld`, `--w-ld-chr` | `--ldscore-dir DIRECTORY` on regression commands | Predictors, weights, and count metadata come from one canonical artifact. Convert legacy fragments with `convert-ldsc2-ldscores` first. |
| `--two-step CUTOFF` | `--two-step-cutoff CUTOFF` | Estimator cutoff; current inclusive threshold and model/intercept restrictions still apply. |
| `--not-M-5-50` | `--count-kind all` | Opt into all-SNP counts. Default remains `--count-kind common`. |
| `munge_sumstats.py --no-alleles` | `--snp-identifier rsid` or `--snp-identifier chr_pos` | Choose allele-blind identity explicitly; default current identity is `chr_pos_allele_aware`. |
| `munge_sumstats.py --merge-alleles FILE` | No exact alias; use `--sumstats-snps-file` or `--use-hm3-snps` for restriction | Current identity matching and downstream allele alignment have separate roles. A keep-list does not imply allele rewriting or liftover. |

Sources: LDSC2 [`ldsc.py`, parser](</Users/wenbinwu/Documents_local/Research/SullivanLab/LDSC/repos/ldsc_py2_Bulik_workspace/ldsc_py2_Bulik/ldsc.py:424>) and [`munge_sumstats.py`, parser](</Users/wenbinwu/Documents_local/Research/SullivanLab/LDSC/repos/ldsc_py2_Bulik_workspace/ldsc_py2_Bulik/munge_sumstats.py:454>); current [`cli.build_parser`](../../src/ldsc/cli.py), [`ldscore_calculator.build_parser`](../../src/ldsc/ldscore_calculator.py), [`regression_runner` argument builders](../../src/ldsc/regression_runner.py), and the [IO inventory, command sections](io-argument-inventory.md). Unlisted unsupported legacy options should not be assumed to have an automatic replacement.

## Naming audit boundary

Current README, `docs/current`, `docs/wiki`, troubleshooting guidance, tutorials, public configuration docstrings, code comments, tests, and the annotation-memory benchmark use the current names. The shared `--threads` spelling remains current for both LD-score computation and gene-index building. References to the abandoned `--workers` spelling are limited to this decision record, historical records, rejection tests, and the annotation-memory benchmark's independent orchestration option. Dated `docs/audits`, `docs/archive`, `docs/plans`, and `docs/specs` may retain historical names and are not the current CLI reference. The original [flag audit](../audits/flag-cleanup/README.md) is a snapshot of behavior before this cleanup.
