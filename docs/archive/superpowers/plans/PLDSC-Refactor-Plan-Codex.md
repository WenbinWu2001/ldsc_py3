# Dedicated Partitioned-LDSC Refactor

## Summary
Implement the refactor described in the April 21 plan, with the live-code adjustments discovered during inspection:

- Extend scope to the coupled modules that the current plan undercounts: `outputs.py` and `ref_panel_builder.py`.
- Treat the `LDScoreResult` / `ChromLDScoreResult` change as a hard public API break.
- Make per-chromosome LD-score output the only supported output mode and remove aggregate/per-chrom toggles rather than keeping compatibility knobs.
- Keep `build-ref-panel` consistent with the new config boundary in the same pass.

## Key Changes
### Config and reference-panel boundary
- Remove `ref_panel_snps_file` and `regression_snps_file` from `GlobalConfig`; keep only identifier/build/logging/shared runtime settings there.
- Add `regression_snps_file` to `LDScoreConfig`.
- Add `ref_panel_snps_file` and `r2_bias_mode` to `RefPanelConfig`.
- Move reference SNP restriction ownership into `_kernel/ref_panel.py` so `load_metadata()` applies the restriction from `RefPanelConfig`, not from `GlobalConfig`.
- Update `ref_panel_builder.py` so `build-ref-panel` also stops depending on `GlobalConfig.ref_panel_snps_file` and instead carries restriction state through its build config/spec path consistently.

### Annotation and LD-score workflow
- Remove gene-set TSV support entirely; BED is the only raw query format.
- Add in-memory BED projection to `AnnotationBuilder.run()` and make `run_bed_to_annot()` / `project_bed_annotations()` return `AnnotationBundle`, with optional `.annot.gz` writing only for the explicit annotate/Python escape hatch.
- Rewrite `run_ldscore_from_args()` to always go through `AnnotationBuilder` + `LDScoreCalculator`; delete the old direct workflow helper path once parity is verified.
- Add `--query-annot-bed` to `ldsc ldscore` and make it mutually exclusive with `--query-annot`.
- Remove aggregate/per-chrom output mode selection from the workflow/output layer; `ldsc ldscore` always writes `<out_prefix>.<chrom>.l2.ldscore.gz` plus count artifacts.

### Public result and output shape
- Replace the split public result shape:
  `reference_metadata` + `ld_scores` + `regression_metadata` + `w_ld`
  with a single merged `ldscore_table`.
- Normalize both `ChromLDScoreResult` and `LDScoreResult` to:
  `ldscore_table`, `snp_count_totals`, `baseline_columns`, `query_columns`, `ld_reference_snps`, `ld_regression_snps`, `output_paths`, `config_snapshot`.
- Make `ldscore_table` the file-equivalent public table:
  rows = `ld_regression_snps`
  columns = `[CHR, SNP, BP, <annotation LD columns>, regr_weight]`
- Set `ld_reference_snps = frozenset()` on normalized/public results because the full reference universe is not recoverable from written rows.
- Keep `chromosome_results` as a list of per-chromosome normalized results so summaries and writers still have explicit per-chrom objects.
- Update `outputs.py` to consume `ldscore_table` directly, emit merged `.l2.ldscore.gz`, stop writing `.w.l2.ldscore.gz`, and compute summaries from the normalized shape.

### Regression loading and consumption
- Rename `_load_ldscore_result_from_files` to public `load_ldscore_from_files`.
- Support both:
  new merged format with embedded `regr_weight`
  legacy separate weight-file format via optional `weight_path`
- Reconstruct public results from disk as normalized results:
  `ldscore_table` from file rows, `ld_reference_snps = frozenset()`, `ld_regression_snps` from `build_snp_id_series(...)`.
- Update `RegressionRunner` and helpers to build datasets from `ldscore_table` instead of separate reference/weight tables.
- Make `--w-ld` optional on regression subcommands and resolve weights from embedded `regr_weight` first when present.

## Test Plan
- Annotation tests:
  BED projection through `AnnotationBuilder.run()` produces binary query columns with no implicit `.annot.gz` writes.
  `run_bed_to_annot()` returns `AnnotationBundle` and only writes annotations when `output_dir` is explicitly provided.
  gene-set APIs/CLI paths are gone.
- Config/ref-panel tests:
  `GlobalConfig` no longer exposes `ref_panel_snps_file` or `regression_snps_file`.
  `LDScoreConfig` accepts `regression_snps_file`.
  `RefPanelConfig` accepts `ref_panel_snps_file` and `r2_bias_mode`.
  `build-ref-panel` still honors SNP restriction through the new boundary.
- LD-score workflow tests:
  direct `AnnotationBuilder` + `LDScoreCalculator` matches `run_ldscore_from_args()`.
  parser accepts `--query-annot-bed` and rejects it with `--query-annot`.
  outputs are per-chromosome only, include `regr_weight`, and never emit `.w.l2.ldscore.gz`.
  public results expose `ldscore_table` and no longer expose split-table fields.
- Regression tests:
  `load_ldscore_from_files` is public.
  merged-format loads work without `weight_path`.
  legacy-format loads still work with `weight_path`.
  `partitioned-h2` works without `--w-ld` when `regr_weight` is embedded.
  existing regression summaries still compute after the dataset build path is switched to `ldscore_table`.

## Assumptions and Defaults
- Hard break: no compatibility properties or deprecation shims for old public result fields.
- Per-chromosome-only output is enforced by removing output-mode switches, not by silently ignoring them.
- `chromosome_results` remains public, but each chromosome result uses the same normalized merged-table shape as the aggregate result.
- The kernel compute layer remains behaviorally unchanged; all normalization happens in the workflow/output layers above it.
