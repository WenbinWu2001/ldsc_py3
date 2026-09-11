# CLI Flag Cleanup Audit

Last updated on: 2026-09-11

Historical snapshot: this report and its inventory/probes describe the pre-cleanup interface. Approved changes are now documented in the [current legacy CLI flag map](../../current/legacy-cli-flag-map.md); the original probe results below are retained as evidence.

The [regression intercept follow-up audit](regression-intercepts.md) revisits the shortcut's redundancy and records the two-step interactions. The approved replacement of `--no-intercept` with existing fixed-value flags is now implemented; see the [implementation record](intercept-consolidation.md). The covariance-only single-annotation `rg` failure remains intentional.

Audited the local `ldsc_py3_restructured` working tree on `restructure`, based on commit `69081b2`, including pre-existing uncommitted changes. No package code, scientific defaults, public interfaces, or existing documentation was changed by this audit.

The parser exposes **218 command-option entries across 13 commands**, excluding built-in help and including two hidden aliases. `munge-sumstats` has **40**. There is little evidence for wholesale deletion: the clearest obsolete surface is the hidden `--exclude-regions` alias in two commands. Most munging options perform real work, but several have misleading names, incorrect help, conditional non-effects, or broken combinations.

## Scope and evidence

- [Complete inventory](inventory.md): all 218 entries, proposed disposition, declaration location, and literal-flag documentation/test coverage.
- [Machine-readable inventory](inventory.json): parser defaults, choices, requirements, help text, declaration locations, and exact matching file lists.
- [Reproduction script](probe_flags.py) and [observed results](probe-results.json): 22 small probes using temporary synthetic inputs, predominantly through the public munging CLI dispatcher.
- Source tracing covers parsing, configuration conversion, input preparation, kernel consumers, output writers, and relevant examples/tests. The command surface comes from [`cli.build_parser`](../../../src/ldsc/cli.py), and the active IO contract is [`docs/current/io-argument-inventory.md`, “Format Policy” and command sections](../../current/io-argument-inventory.md).

“Unused in practice” cannot be established from this repository alone: there is no user-execution telemetry, and no HPC jobs, external scripts, shell history, or unrelated directories were inspected. Literal mentions are evidence of discoverability, not execution or scientific necessity. Python API tests often exercise a feature without mentioning its CLI spelling. The audit therefore distinguishes **confirmed redundancy**, **conditional non-effect**, **behavior/help defects**, and **optional consolidation**.

## Highest-priority findings

| Priority | Flags | Verified behavior and consequence | Recommendation |
| --- | --- | --- | --- |
| High | munge `--N`, `--N-cas`, `--N-con` | Help says supplied constants take priority over input columns. The implementation uses constants only after failing to obtain per-variant N. Input N=1000 plus `--N 9999` produces N=1000. | State the fallback rule accurately and reject contradictory explicit strategies, or deliberately choose new precedence. Do not silently change numerical behavior as a rename. |
| High | munge `--n-min` | Help says default is the 90th percentile divided by 2; code divides by **1.5**. Zero is treated as omitted. `--n-min 2000 --N 1000` does not filter constant-N rows. | Correct help; distinguish `None` from zero; define applicability to constant N explicitly. |
| High | munge `--nstudy-min` | NSTUDY is removed from the read schema when per-variant N or a complete case/control pair is present. An explicit study-count threshold then has no effect. Zero also requests the maximum-study default. | Retain as an advanced constant-N/meta-analysis filter; reject an explicit ineffective request or deliberately make it independent of N. |
| High | munge `--format daner-new` | `auto` succeeds on a table with `Nca/Nco` and no `FRQ_U_*`; explicit `daner-new` raises `IndexError`. Explicit mode also rejects `NCAS/NCON` aliases that automatic inference recognizes. | Unify explicit and inferred preparation. Strong candidate to retire the dedicated new-DANER kernel branch while retaining normal column/count handling. |
| High | munge `--infer-only` | Reports `Runnable: yes` for missing P and ambiguous N strategies. Suggested commands omit supplied options such as `--n-min`, `--keep-maf`, `--trait-name`, and rsID identity mode. It may suggest explicit `daner-new`, triggering the failure above. | Retain this useful mode, but share structural validation and preserve effective configuration when rendering commands. A diagnostic suggestion must not change the requested analysis. |
| Medium | munge `--info-list` | A single column containing values such as `0.99,NA` works. Two column names map to duplicate canonical INFO fields and are rejected before the multi-column kernel path can run. | Either make it explicitly one list-valued column (`--info-list-col`) or implement a deliberate multi-column aggregation contract. |
| Medium | munge `--info-min`, `--maf-min` | Explicit thresholds succeed without filtering anything when the corresponding input field is absent. Probes with `--info-min 10` and absent INFO, or `--maf-min .5` and absent frequency, retain all valid rows. | Keep both controls. Distinguish an omitted default from an explicitly requested filter and fail clearly when the requested field is missing. |
| Medium | munge `--keep-maf` | Preserves original FRQ, including 0.8 and 0.7. Folding to minor-allele frequency is local to the QC mask. | Rename to `--keep-frequency`, or remove the flag and always preserve supplied frequency after an explicit output-schema decision. |
| Medium | munge `--a1-inc` | Suppresses inferred effect columns and emits positive Z from P. A negative input BETA becomes positive Z when the option is supplied. | Rename to `--a1-is-increasing`; keep only as an advanced explicit input assumption. It does not infer the increasing allele or swap A1/A2. |
| Low | ldscore and gene-index `--exclude-regions` | Hidden alias writes exactly the same destination as `--regr-snps-exclude-regions`; both accept the same choices. | Remove both aliases. The explicit spelling communicates which SNP universe is affected. |

Evidence: [`sumstats_munger.build_parser`, lines 1046–1107](../../../src/ldsc/sumstats_munger.py); [`_kernel.sumstats_munger.process_n`, line 567](../../../src/ldsc/_kernel/sumstats_munger.py); [`_sumstats_input.prepare_munge_input`, DANER handling and NSTUDY removal](../../../src/ldsc/_sumstats_input.py); [`infer_raw_sumstats`, `_apply_raw_sumstats_inference`, `_render_inference_report`](../../../src/ldsc/sumstats_munger.py); [`parse_dat`, `filter_frq`, `munge_sumstats`](../../../src/ldsc/_kernel/sumstats_munger.py). Every behavioral claim in the table except allele-swapping absence and conditional case/control details has a corresponding direct probe in [probe-results.json](probe-results.json); those remaining details were traced in source.

The 90th-percentile example uses N values 550 and 1000. The percentile is 955, so division by 2 would retain both rows, while the implemented division by 1.5 removes the 550 row. The probe observes the latter. A numerical-default decision should precede any attempt to “fix” this by changing the denominator.

## All 40 munge-sumstats flags

The following groups cover the complete munging surface. Names below are proposals, not newly available flags.

| Existing flags | Actual role | Disposition |
| --- | --- | --- |
| `--raw-sumstats-file`, `--output-dir`, `--overwrite`, `--log-level` | Resolve one input, own a result directory, authorize replacement, route logs. | Keep. Uniform directory/logging controls are useful. |
| `--trait-name` | Biological trait label carried into artifacts and regression summaries. | Keep; downstream labels are meaningful. |
| `--format` | Raw-input profile, separate from output serialization. | Rename `--input-format`; repair DANER consistency. Keep old-DANER header-count extraction unless supported inputs are deliberately reduced. |
| `--infer-only` | Diagnose input and render a proposed run without writing. | Keep and repair. Requiring `--output-dir` is currently deliberate; making it optional is a workflow-contract choice, not dead-code cleanup. |
| `--output-format` | `parquet`, legacy `tsv.gz`, or both. | Keep interoperability. Consider always writing canonical Parquet and replacing this with an optional legacy-export switch; this would remove the metadata-poor TSV-only mode, not the export capability. |
| `--sumstats-snps-file`, `--use-hm3-snps` | Custom keep-list versus packaged HM3 restriction. | Keep both capabilities. The packaged flag avoids requiring a user to locate internal resources. `--keep-snps-file` is a clearer candidate spelling for the former. |
| `--source-genome-build`, `--output-genome-build` | Identify input coordinate build and select final build. | Keep distinct. Combining them would obscure liftover and source interpretation. |
| `--liftover-chain-file`, `--use-hm3-quick-liftover` | Chain liftover versus packaged HM3 coordinate mapping. | Keep distinct methods. Consider `--liftover-hm3`; document its required HM3 restriction in help. Do not make filtering or liftover implicit without a deliberate design decision. |
| `--snp-identifier` | Controls identity, duplicate groups, keep-list matching, and provenance. | Keep. A clearer cross-command name could be `--snp-identity-mode`, but the change would be broad. |
| `--N` | Constant sample-size fallback. | Keep as `--sample-size`; clarify precedence and positivity/finite-value expectations. |
| `--N-cas`, `--N-con` | Constant fallback counts are added to produce one N. | Strong consolidation candidates: numerically equivalent to a supplied total N in this branch. Preserve count provenance explicitly if removing them. |
| `--N-col` | Explicit direct per-variant N selection. | Keep; normalize capitalization to `--n-col`, or use `--sample-size-col`. |
| `--N-cas-col`, `--N-con-col` | Per-variant case/control strategy; both required. | Keep as `--case-count-col` and `--control-count-col`. These are not interchangeable with constant-N flags or direct effective-N columns. |
| `--snp`, `--chr`, `--pos`, `--a1`, `--a2`, `--p` | Explicit column names, not values or filters. | Keep as `--snp-col`, `--chr-col`, `--pos-col`, `--a1-col`, `--a2-col`, `--p-col`. Automatic inference does not eliminate legitimate overrides. |
| `--frq`, `--info`, `--nstudy` | Frequency, INFO, study-count column selectors. | Keep as `--frequency-col`, `--info-col`, `--study-count-col`. Study-count mapping is advanced/conditional. |
| `--signed-sumstats` | One column plus its null value, used for sign and median checking; P supplies Z magnitude. | Keep; rename to `--signed-stat COLUMN,NULL` and explain the separation of sign and magnitude. Splitting this into two flags would increase the count without necessarily improving use. |
| `--info-list` | Parse numeric/NA comma-separated values within a selected INFO cell. | Keep the supported input capability; repair or narrow the advertised plural-column behavior. |
| `--ignore` | Exclude headers from inference, including explicit hints; not a SNP filter. | Keep as `--ignore-cols`. It remains useful for duplicate aliases and unwanted QC fields; explicit N strategy selection has removed one former workaround use. |
| `--info-min`, `--maf-min` | QC thresholds when the corresponding data exist. | Keep; fix explicit missing-field requests. `--maf-min` correctly refers to folded MAF even though stored FRQ is not necessarily MAF. |
| `--n-min`, `--nstudy-min` | Whole-table N filtering or fallback study-count filtering. | Keep as advanced controls with the behavior corrections above; consider `--study-count-min` for consistency. |
| `--chunksize` | Number of raw rows read per chunk. | Keep as `--chunk-size`, with units and memory trade-off in help. Retained chunks are concatenated; this is not a bound on total retained-table memory. |
| `--a1-inc` | Assert that A1 is already the increasing allele for every row. | Advanced niche candidate; rename first. Absence from normal examples does not prove it is safe to remove. |
| `--keep-maf` | Carry FRQ through output rather than drop it after QC. | Rename or remove by always preserving frequency. It does not control the MAF filter. |

The per-variant case/control path computes total count T and case fraction p, then uses `N_i = T_i * p_i / mean(p_j among rows with maximum T)`. Consequently, the paired **column** flags must not be removed on the grounds that the paired **constant** flags merely add two numbers. See [`process_n`](../../../src/ldsc/_kernel/sumstats_munger.py) and [“Sample-Size Column Selection”](../../current/munge-sumstats.md).

`--keep-maf` is not needed by the regression estimator: FRQ is carried and may be reoriented during allele alignment, while fitted numerical inputs use Z, N, LD scores, weights, and counts. Frequency can still be useful to downstream consumers outside this package. See [`RegressionRunner._fit_h2_dataset`, `_fit_rg_dataset`, and legacy-frequency handling](../../../src/ldsc/regression_runner.py).

The main help burden can be reduced without deleting supported input formats: group flags into input/output, column mapping, sample-size strategy, QC, identity/liftover, and advanced execution. Twelve simple column selectors alone belong together. Add clear metavar values and defaults instead of a flat list of 40 equally prominent options.

## Other commands

Common `--output-dir`, `--overwrite`, and `--log-level` flags remain active. The lack of `--output-dir` on `plot` and `convert-h2-scale` is an approved fixed-destination contract, not missing functionality. The [complete inventory](inventory.md) enumerates every flag individually; the table below records the command-specific conclusions.

| Command | Entries | Audit conclusion |
| --- | ---: | --- |
| `annotate` | 12 | No confirmed dead flag. Keep BED/gene alternatives, baseline sources, gene catalog/policy/exclusion, padding, and identity/build controls. Clarify `--snp-identifier` beyond “bundle validation”; defaults and permitted gene-only options deserve grouped help. |
| `ldscore` | 34 | Remove hidden `--exclude-regions`. Rename `--yes-really` to `--allow-whole-chromosome-window` and `--threads` to `--workers` (it uses processes). Separate direct PLINK, stored R2, and indexed-gene options in help. Keep `--ref-panel-snps-file` distinct from `--regr-snps-file`: contributors versus written/regression rows. Keep `--maf-min` distinct from `--common-maf-min`: filtering versus count-vector threshold. |
| `build-ref-panel` | 18 | No confirmed dead flag. SNP/kb/cM windows specify different units. Two directional chain flags allow different emitted builds and are not simple duplicates. `--min-r2` is scientifically consequential sparsification of unbiased R2, not a generic read/QC filter; clearer name: `--min-stored-r2`. Keep PLINK batching, individual restrictions, reference SNP universe, maps, build and identity controls. |
| `build-gene-ldscore-index` | 21 | Remove hidden `--exclude-regions`. `--genome-build` has only the choice `hg19`, but it is a deliberate required assertion, not an accidentally dead parameter. Keep SNP and atom batching: they bound different work. `--threads` uses an actual thread pool here; `--workers` would give both commands a consistent user-level concept. Gene exclusion defaults to MHC here but none in live annotation workflows; make that visible. |
| `convert-ldsc2-ldscores` | 8 | No confirmed dead flag. Reference and weight directories are distinct inputs; optional frequencies supply metadata. Add missing help and explain accepted legacy suite layouts and build assertions. |
| `h2` | 15 | Rename `--no-intercept`: it fixes the intercept to **1**, not zero or an absent model term. It duplicates `--intercept-h2 1` numerically, but is a useful shared policy spelling. Improve `--n-blocks` to `--jackknife-blocks`; `--count-kind` to `--reference-snp-count-kind`. Keep explicit intercept, cutoff, prevalence, identity-downgrade, and count choices. Update sumstats help to lead with canonical Parquet. |
| `partitioned-h2` | 17 | Same regression naming issues. `--two-step-cutoff` is invalid for actual multi-annotation fits; label this restriction and reject before expensive work. Single-column models can still reach the estimator, so this is not universally unread code. Keep `--query-batch-size` and `--summary-sort-by`: loading bound versus presentation. Query batching fits each query separately, not a joint query model. |
| `quantile-h2` | 18 | No confirmed dead flag. Baseline/query sources reconstruct the fitted model; target sources define quantile assignment; reference metadata define the reference universe. These superficially repetitive paths are not interchangeable. `--target-annotation` is a column name; consider `--target-annotation-col`. `--target-missing-value` is an exclusion sentinel, not an imputation value: use `--target-missing-token`. Many flags currently have no help. Align explanation of gene padding with the other workflows. |
| `rg` | 18 | `--no-intercept` fixes h2 intercepts to **1** and covariance intercept to **0**. `--chisq-max c` filters `Z1^2 * Z2^2 <= c^2`, equivalently `abs(Z1*Z2) <= c`, whereas h2 applies `Z^2 <= c`; the same name hides different operations. Retain but document the per-command meaning or use an explicit product name for rg. The ordered prevalence lists overlap with `--prevalence-manifest`; prefer the named manifest for batches, retain lists as a small-run convenience unless deliberate consolidation is desired. Keep anchor selection and per-pair detail. |
| `query-r2` | 7 | No confirmed dead flag. `--pairs` should be `--pairs-file` (keep `-` for stdin). `--snp-identifier` can deliberately rekey queries, not merely assert the panel's stored identity; explain or rename to `--query-identity-mode`. `--genome-build` selects a build subdirectory, not liftover. Align `--panel-dir`/`--r2-dir` terminology only after making accepted root/subdirectory behavior explicit. |
| `convert-h2-scale` | 7 | `--num-points` only affects prevalence ranges and is ignored for exact `--pop-prev`; the probe even accepts a negative count in exact mode. Keep and gate it by mode. Rename `--samp-prev`, `--pop-prev`, and `--pop-prev-range` to full sample/population-prevalence names consistently across regression and conversion commands. Existing inline prevalence conversion and this postprocessing workflow serve different stages. |
| `plot` | 3 | All three flags are useful: result selection, overwrite, and logging. No cleanup beyond consistent help is indicated. |

Parser and consumer references: [`annotation_builder.add_annotate_arguments` / `run_annotate_from_args`](../../../src/ldsc/annotation_builder.py); [`ldscore_calculator.build_parser`, `_run_explicit_indexed_ldscore`, `_ldscore_config_from_args`, `_resolve_worker_count`](../../../src/ldsc/ldscore_calculator.py); [`ref_panel_builder.build_parser` / `config_from_args`](../../../src/ldsc/ref_panel_builder.py); [`gene_ldscore_index.build_parser` / `run_build_gene_ldscore_index_from_args`](../../../src/ldsc/gene_ldscore_index.py); [`legacy_ldscore_converter.main`](../../../src/ldsc/legacy_ldscore_converter.py); [`regression_runner` argument builders, `_validate_intercept_conflicts`, `_fit_h2_dataset`, `_fit_rg_dataset`](../../../src/ldsc/regression_runner.py); [`_kernel.regression.LD_Score_Regression`, two-step validation](../../../src/ldsc/_kernel/regression.py); [`_quantile_inputs._model_annotations` / `prepare_quantile_statistics`](../../../src/ldsc/_quantile_inputs.py); [`_quantile_storage.target_chunks`](../../../src/ldsc/_quantile_storage.py); [`r2_query.R2Panel.open` / `_chrom_state`](../../../src/ldsc/r2_query.py); [`h2_scale._population_prevalence_grid`](../../../src/ldsc/h2_scale.py); [`plotting.add_plot_arguments` / `run_plot_from_args`](../../../src/ldsc/plotting/__init__.py).

## Conditional non-effects versus invalid combinations

These are important for deciding what to remove. A backend-specific option can be useful without being meaningful in every mode.

- **Stored-R2 ldscore:** `--snp-batch-size` does not affect pair streaming; genetic-map options cannot replace authoritative sidecar CM; `--export-ref-metadata` does not create an export when metadata already come from a Parquet sidecar. An explicit unsupported request should fail or be visibly explained. `--keep-indivs-file` already raises an error in this backend, so it is not a silent no-op. See [`ParquetR2RefPanel.load_metadata`, `prepare_chromosome`, `_resolve_genetic_map`](../../../src/ldsc/_kernel/ref_panel.py) and [`LDScoreCalculator.run`, export-directory selection](../../../src/ldsc/ldscore_calculator.py).
- **Indexed-gene ldscore:** the workflow already rejects direct reference, build, window, map, and many resource knobs, including explicitly supplied values equal to defaults. Do not list those as dead; make the accepted mode matrix discoverable. See [`_run_explicit_indexed_ldscore`](../../../src/ldsc/ldscore_calculator.py).
- **Munging liftover:** same-build runs warn and ignore an otherwise supplied method; rsID runs reject concrete build and liftover settings. `--source-genome-build auto` is accepted in rsID mode despite broader documentation wording saying source-build flags are rejected. See [`_resolve_main_global_config`, `_validate_liftover_request_before_io`](../../../src/ldsc/sumstats_munger.py).
- **Regression:** fixed intercept plus two-step estimation and multi-annotation plus two-step estimation are already rejected in the numerical layer. `--no-intercept` plus explicit intercept overrides is already rejected by the workflow. Improve where and how users discover these restrictions; do not label rejected combinations silently ignored.
- **Prevalence conversion:** exact and range modes share a parser, but grid size matters only to ranges. This is a small scope/help correction, not reason to remove range functionality.

## Stale guidance and naming conventions

Several recovery messages currently recommend flags that do not exist: `--snp-col` and `--a1-col` in `SumstatsTable.validate`, `--daner-new` in input preparation, and `--debug` in kernel errors. A kernel log also calls the input `--sumstats`. The actual current spellings are `--snp`, `--a1`, `--format daner-new`, `--log-level DEBUG`, and `--raw-sumstats-file`. Rename implementation and all recovery advice together; do not add compatibility aliases just to make old messages work. Sources: [`SumstatsTable.validate`](../../../src/ldsc/sumstats_munger.py), [`prepare_munge_input`](../../../src/ldsc/_sumstats_input.py), [`filter_info`, `parse_dat`, `process_n`](../../../src/ldsc/_kernel/sumstats_munger.py).

Useful conventions for a deliberate cleanup:

1. Use lowercase hyphenated names; distinguish constant values from `*-col` selectors.
2. Spell out ambiguous abbreviations when doing so improves interpretation: sample/population prevalence, individual keep file, jackknife blocks, workers.
3. Preserve distinctions among reference SNP contributors, regression/output SNPs, and raw-sumstats keep-lists. Similar names here reflect different scientific universes.
4. Keep one required output-directory convention for materializing commands and the documented fixed destinations for derived commands.
5. Make mode restrictions, defaults, units, and accepted list/path syntax visible. For example, `*-sources` accepts `nargs='+'` in some commands and one comma-separated argument in ldscore/gene-index commands; harmonize parsing deliberately rather than merely changing labels.

## Recommended cleanup decisions

**Remove with high confidence:** the two hidden `--exclude-regions` aliases. They add no capability and retain a less precise spelling.

**Consolidate after a deliberate interface decision:** fixed `--N-cas`/`--N-con` into total N; dedicated new-DANER preparation into shared aliases/hints; `--keep-maf` into automatic frequency preservation; optional TSV-only munging into canonical Parquet plus export. Each has a numerical, provenance, input-support, or artifact consequence even when it reduces flag count.

**Retain but make advanced:** explicit column overrides, INFO-list cells, study-count filtering, increasing-allele input, chunk/batch sizes, custom SNP universes, regression estimator controls. The repository does not establish that these capabilities are unused outside its own examples.

**Repair before changing defaults:** sample-size precedence/help, N threshold denominator and zero handling, explicit ineffective QC requests, DANER dispatch, and inference-only validation/command generation. These can affect results or make a suggested run fail, so cosmetic cleanup should not obscure them.

Implementation remains a separate step. The repository treats public flags and artifacts as deliberate contracts while also directing removal of obsolete compatibility paths. This audit supplies the concrete evidence needed to choose the removals and behavioral rules; it does not assume approval for a new public interface.

## Verification

- The installed development Python imports this workspace's `src/ldsc`, not the sibling main-branch repository.
- All 13 real `python -m ldsc <command> --help` invocations exited 0.
- Runtime parser inventory: 218 entries; all 218 matched to AST-located argument declarations.
- `python -m pytest tests/test_sumstats_munger.py tests/test_munge_results.py -q`: **150 passed, 3 warnings, 6 subtests passed**. All three reported warnings concern pandas dtype assignment deprecations.
- `python docs/audits/flag-cleanup/probe_flags.py`: **22 observations** saved in `probe-results.json`, including expected successes, two controlled public input errors, and the explicit-DANER `IndexError`. These are observation probes, not a regression suite claiming current behavior is correct.
- No full numerical suite or production GWAS/HPC run was performed for this documentation-only audit. Conditional findings outside the targeted probes are based on source tracing.
