# Threshold Comparison Findings

## Initial Scope

- Active target repo: `ldsc_py3_restructured` on branch `restructure`.
- Main implementation package: `src/ldsc`.
- Broad threshold search surfaced explicit MAF filtering in `src/ldsc/_kernel/ldscore.py`, metadata contract entries in `src/ldsc/outputs.py`, and validation/inference thresholds in `src/ldsc/genome_build_inference.py`.
- Need distinguish user-facing thresholds from ordinary loop/index/bounds comparisons before reporting.

## Confirmed Implementation Semantics

- Sumstats P-value filter keeps `0 < P <= 1`.
- Sumstats INFO keeps `INFO >= info_min` (single column) and average-like multi-column threshold `sum(INFO) >= info_min * n_columns`.
- Sumstats MAF keeps folded `FRQ > maf_min`.
- Sumstats sample-size filters keep `N >= n_min` and `NSTUDY >= nstudy_min`.
- LD-score retained SNP MAF filters keep `MAF > maf_min`; common-count/overlap masks keep `MAF > common_maf_min`.
- Reference-panel MAF filtering keeps `MAF > maf_min`; PLINK backend also keeps folded genotype MAF `> mafMin`.
- HM3 coordinate-reference builder keeps `MAF >= min_maf`, unlike the other MAF filters.
- Build-ref-panel pair filter emits pairs with `r2 >= min_r2` when `min_r2 > 0`; `min_r2 <= 0` disables the filter.
- Regression `chisq_max` keeps `chisq < chisq_max`; cross-trait two-step keeps both `z^2 < twostep`.
- Genome-build inference requires `informative_count >= MIN_INSPECTED_REFERENCE_SNPS` and `best_fraction >= DOMINANCE_THRESHOLD`.

## Notes

- Config validation commonly accepts threshold parameters inclusively within allowed numeric ranges, for example `0 <= maf_min <= 0.5`; these validations are not row-retention masks.
- `ldscore_calculator.py` docstrings call the common-MAF threshold "inclusive", but the active implementation and output metadata use strict `>`.

## Final Matrix

| Area | User/internal threshold | Retained / accepted condition | Equality kept? | Primary implementation |
|---|---|---:|---:|---|
| `munge-sumstats` | P-value range | `P > 0` and `P <= 1` | lower no, upper yes | `src/ldsc/_kernel/sumstats_munger.py:210` |
| `munge-sumstats` | `--info-min` | `INFO >= info_min`; multi-INFO uses `sum(INFO) >= info_min * n_cols` | yes | `src/ldsc/_kernel/sumstats_munger.py:222`, `:226` |
| `munge-sumstats` | `--maf-min` | folded `FRQ > maf_min` | no | `src/ldsc/_kernel/sumstats_munger.py:256` |
| `munge-sumstats` | `--n-min` | `N >= n_min` | yes | `src/ldsc/_kernel/sumstats_munger.py:670` |
| `munge-sumstats` | `--nstudy-min` | `NSTUDY >= nstudy_min` | yes | `src/ldsc/_kernel/sumstats_munger.py:677` |
| `ldscore` | `--maf-min` for annotation/reference metadata | `MAF > maf_min` | no | `src/ldsc/_kernel/ldscore.py:972` |
| `ldscore` | `--maf-min` with PLINK backend | folded genotype MAF `> mafMin` and heterozygote/missing count `< n_eff` | no for MAF | `src/ldsc/_kernel/plink_bed.py:317` |
| `ldscore` | `--common-maf-min` counts | `MAF > common_maf_min` | no | `src/ldsc/_kernel/ldscore.py:1572` |
| `ldscore` overlap | common universe | `MAF > common_maf_min` | no | `src/ldsc/_kernel/overlap.py:61` |
| `runtime ref panel` | `RefPanelConfig.maf_min` | `MAF > maf_min` | no | `src/ldsc/_kernel/ref_panel.py:260` |
| `build-ref-panel` | `ReferencePanelBuildConfig.maf_min` / `--maf-min` | PLINK folded genotype MAF `> mafMin` | no | `src/ldsc/_kernel/plink_bed.py:317` via `src/ldsc/ref_panel_builder.py:926` |
| `build-ref-panel` | `--min-r2` | if `min_r2 > 0`, keep `r2 >= min_r2`; `min_r2 <= 0` disables filtering | yes when enabled | `src/ldsc/_kernel/ref_panel_builder.py:426-428` |
| `build-ref-panel` / `ldscore` | LD window distance | drops while distance `> max_dist`, so retained pairs are `distance <= max_dist` | yes | `src/ldsc/_kernel/ref_panel_builder.py:328`, `src/ldsc/_kernel/ldscore.py:355` |
| region exclusion | BED interval | 1-based `p` excluded iff `start < p <= end` | excludes upper endpoint, not lower | `src/ldsc/_kernel/regions.py:273-275`, `:307-311` |
| annotation BED projection | BED overlap | SNP `[POS-1, POS)` overlaps BED with count `> 0` | BED half-open semantics | `src/ldsc/_kernel/annotation.py:81-88`, `:132` |
| `hm3_reference` builder | `min_maf` | `MAF >= min_maf` | yes | `src/ldsc/hm3_reference.py:54` |
| regression h2 | `--chisq-max` | `chisq < chisq_max` | no | `src/ldsc/regression_runner.py:607` |
| regression rg | `--two-step-cutoff` | `z1**2 < twostep` and `z2**2 < twostep` | no | `src/ldsc/_kernel/regression.py:632` |
| genome-build inference | minimum inspected reference SNPs | fails when inspected `< MIN_INSPECTED_REFERENCE_SNPS`, so accepts `>=` | yes | `src/ldsc/genome_build_inference.py:248`, `:460` |
| genome-build inference | dominance threshold | `best_fraction >= DOMINANCE_THRESHOLD` | yes | `src/ldsc/genome_build_inference.py:251`, `:472` |
| coordinate validation | minimum base-pair position | invalid when position `< min_position`, so accepts `>=` | yes | `src/ldsc/_coordinates.py:95`, `:309` |
| overlap diagnostics | condition-number threshold | warning only when `condition_number > threshold`; returns clean at `<=` | equality is clean | `src/ldsc/overlap_matrix.py:207` |
| allele canonicalization | minor-allele flip | flips when frequency `> 0.5`; ties stay ordered | equality not flipped | `src/ldsc/_kernel/snp_identity.py:642` |

## Direct Test/Doc Evidence

- `tests/test_sumstats_munger.py:342-354` pins P, INFO, and FRQ boundary behavior.
- `tests/test_ldscore_workflow.py:1276-1291` pins strict `common_maf_min` and metadata operator `">"`.
- `tests/test_overlap_matrix.py:18-31` pins `MAF == 0.05` as not common.
- `tests/test_ref_panel_builder.py:334-379` pins `min_r2=0.0` as unfiltered and positive `min_r2` as an emission floor.
- `docs/current/io-argument-inventory.md:183` says `--common-maf-min` is strict `MAF > common_maf_min`.
- `docs/current/ld-window-parquet-r2-sidecar-behavior.md:78` conflicts with code/tests by saying common counts use `MAF >= common_maf_min`.
