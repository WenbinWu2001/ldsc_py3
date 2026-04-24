# Partitioned-LDSC Refactor: Code Review Checklist

Review against:
- `docs/design/partitioned-ldsc-workflow.md` — SNP universe semantics, output formats, CLI reference
- `docs/superpowers/plans/2026-04-21-partitioned-ldsc-refactor.md` — task-by-task implementation plan
- `docs/superpowers/plans/PLDSC-Refactor-Plan-Codex.md` — extended scope and hard constraints

Each item targets a specific constraint that is easy to get wrong or to leave half-done.

---

## 1. Config boundary (`config.py`, `GlobalConfig`, `LDScoreConfig`, `RefPanelSpec`)

- [ ] `GlobalConfig` does **not** expose `ref_panel_snps_path` — field removed entirely.
- [ ] `GlobalConfig` does **not** expose `regression_snps_path` — field removed entirely.
- [ ] `LDScoreConfig` has a `regression_snps_path` field (moved from `GlobalConfig`).
- [ ] `RefPanelSpec` has a `ref_panel_snps_path` field (moved from `GlobalConfig`).
- [ ] `RefPanelSpec` has an `r2_bias_mode: str | None = None` field.
- [ ] No callers of `GlobalConfig` still read `.ref_panel_snps_path` or `.regression_snps_path` — grep confirms zero uses.
- [ ] `_normalize_run_args()` in `ldscore_calculator.py` routes `--ref-panel-snps` into `RefPanelSpec`, not `GlobalConfig`.
- [ ] `_normalize_run_args()` routes `--regression-snps` into `LDScoreConfig`, not `GlobalConfig`.
- [ ] `RefPanel.load_metadata()` applies restriction from `RefPanelSpec.ref_panel_snps_path`, not from `GlobalConfig`.
- [ ] `ref_panel_builder.py` passes SNP restriction through its build config / spec path — no dependency on `GlobalConfig.ref_panel_snps_path`.

---

## 2. Gene-set TSV removal (`_kernel/annotation.py`, `annotation_builder.py`, `cli.py`)

- [ ] `AnnotationSourceSpec` has no `gene_set_paths` field.
- [ ] `gene_set_to_bed`, `make_annot_files`, `parse_make_annot_args`, `main_make_annot` are absent from `_kernel/annotation.py`.
- [ ] `annotation_builder.py` `__all__` exports no gene-set symbols.
- [ ] `cli.py` `_run_annotate()` has no gene-set dispatch branch.
- [ ] `cli.py` `_add_annotate_arguments()` registers none of: `--gene-set-file`, `--gene-coord-file`, `--windowsize`, `--bed-file`, `--nomerge`, `--bimfile`, `--annot-file`.

---

## 3. In-memory BED projection (`_kernel/annotation.py`)

- [ ] `AnnotationBuilder.run()` handles non-empty `AnnotationSourceSpec.bed_paths` — no silent ignore.
- [ ] `_compute_bed_query_columns()` returns a `pd.DataFrame`. No `.annot.gz` files are written as a side-effect.
- [ ] BED column names are taken from file stems. Duplicate stems raise `ValueError` before any computation.
- [ ] Stems that clash with loaded baseline/query column names raise `ValueError`.
- [ ] `run_bed_to_annot()` returns `AnnotationBundle`, not `Path`.
- [ ] `project_bed_annotations()` returns `AnnotationBundle`. File writing only occurs when `output_dir` is explicitly provided.
- [ ] `_write_bundle_query_as_annot_files()` is called only from `project_bed_annotations()` / `run_bed_to_annot()` — not from `run_ldscore_from_args()` or any ldscore CLI path.

---

## 4. SNP universe filtering — `ref_panel_snps_path`

- [ ] `AnnotationBuilder._run_single_universe()` keeps the full annotation-file universe `B`; it does **not** apply `ref_panel_snps_path` on the annotation side.
- [ ] `RefPanel.load_metadata()` applies `RefPanelSpec.ref_panel_snps_path` at load time — not deferred to the caller.
- [ ] `LDScoreCalculator.compute_chromosome()` aligns each chromosome-local annotation bundle to `ref_panel.load_metadata(chrom)` before the legacy kernel call, so the kernel sees `B ∩ A'`.
- [ ] `_filter_annotation_bundle_to_reference_universe()` in `ldscore_calculator.py` is **deleted** (dead code after convergence).
- [ ] `_chromosome_set_from_annotation_inputs()` in `ldscore_calculator.py` is **deleted**.

---

## 5. SNP universe fields on results (`ldscore_calculator.py`, `regression_runner.py`)

- [ ] `ChromLDScoreResult` has `ld_reference_snps: frozenset[str]` — old name `reference_snps` is gone.
- [ ] `ChromLDScoreResult` has `ld_regression_snps: frozenset[str]` — old name `regression_snps` is gone.
- [ ] Same rename on `LDScoreResult`.
- [ ] `_wrap_legacy_chrom_result()` derives the regression row mask from the compute-time reference metadata, so `ld_regression_snps` is a subset of the internal reference universe by construction before normalization.
- [ ] On normalized / public results produced in memory or reconstructed from disk via `load_ldscore_from_files`, `ld_reference_snps = frozenset()` — the full reference universe is intentionally not recoverable from normalized row tables.
- [ ] `ld_regression_snps` on disk-loaded results is reconstructed via `build_snp_id_series(ldscore_table, snp_identifier)`.
- [ ] No callers still reference `.reference_snps` or `.regression_snps` — grep confirms zero uses.

---

## 6. Public result shape — `ldscore_table` (`ldscore_calculator.py`, `regression_runner.py`)

- [ ] `ChromLDScoreResult` has a single `ldscore_table: pd.DataFrame` field. Columns: `[CHR, SNP, BP, <annotation>_L2 ..., regr_weight]`. Rows = `ld_regression_snps`.
- [ ] `LDScoreResult` has the same `ldscore_table` shape (aggregate across chromosomes).
- [ ] Old split fields (`reference_metadata`, `ld_scores`, `regression_metadata`, `w_ld`) are **absent** from both dataclasses.
- [ ] `ldscore_table` always contains a `regr_weight` column — never absent.
- [ ] `snp_count_totals` is computed over `ld_reference_snps`, not `ld_regression_snps`.
- [ ] `chromosome_results: list[ChromLDScoreResult]` is preserved; each element uses the same normalized merged-table shape.
- [ ] `_subset_ldscore_result()` in `regression_runner.py` operates on `ldscore_table` columns, not split tables.

---

## 7. Output files (`outputs.py`, `ldscore_calculator.py`)

- [ ] `ldsc ldscore` writes `<out_prefix>.<chrom>.l2.ldscore.gz` — per-chromosome only; no aggregate output mode.
- [ ] The per-chromosome file contains a `regr_weight` column. No `.w.l2.ldscore.gz` is written.
- [ ] `.M` and `.M_5_50` counts are accumulated over `ld_reference_snps` — the count loop does not use `ld_regression_snps` rows.
- [ ] Output-mode selection switches (`per_chr_output`, `aggregate_output`, or equivalent flags) are removed, not silently defaulted.
- [ ] `outputs.py` reads from `ldscore_table` directly, not from split tables.

---

## 8. `run_ldscore_from_args()` convergence (`ldscore_calculator.py`)

- [ ] `run_ldscore_from_args()` builds an `AnnotationSourceSpec` and calls `AnnotationBuilder(global_config).run(source_spec)` — no inline annotation loading.
- [ ] `run_ldscore_from_args()` builds a `RefPanel` via `_ref_panel_from_args()` using `RefPanelSpec`.
- [ ] `run_ldscore_from_args()` builds an `LDScoreConfig` via `_ldscore_config_from_args()`.
- [ ] `run_ldscore_from_args()` calls `LDScoreCalculator().run(...)` — no direct `kernel_ldscore.combine_annotation_groups()` call.
- [ ] `_ref_panel_from_args()` constructs `RefPanelSpec` with `ref_panel_snps_path` sourced from the args.
- [ ] `_ldscore_config_from_args()` constructs `LDScoreConfig` with `regression_snps_path`.

---

## 9. CLI (`ldscore_calculator.py`, `cli.py`)

- [ ] `build_parser()` uses a `mutually_exclusive_group` for `--query-annot` and `--query-annot-bed`.
- [ ] Passing both `--query-annot` and `--query-annot-bed` at the same time causes `SystemExit` at parse time.
- [ ] `--query-annot-bed` help text states mutual exclusivity with `--query-annot` and that no `.annot.gz` files are written.
- [ ] `_normalize_run_args()` defaults `query_annot_bed` to `None` when the attribute is absent.
- [ ] `cli.py` `_run_annotate()` requires `--bed-files` and has no gene-set dispatch branch.
- [ ] `ldsc ldscore --help` output includes `--query-annot-bed`.

---

## 10. Regression loading (`regression_runner.py`, `__init__.py`)

- [ ] `load_ldscore_from_files` is a public function (no leading underscore).
- [ ] `load_ldscore_from_files` is exported from `src/ldsc/__init__.py`.
- [ ] New-format load (no `weight_path`): detects `regr_weight` column and uses it as regression weights.
- [ ] Legacy-format load: `weight_path` provided → loads weights from separate `.w.l2.ldscore.gz`.
- [ ] When `weight_path` is absent and no `regr_weight` column is present, `ValueError` is raised.
- [ ] `--w-ld` is optional on `partitioned-h2`, `h2`, and `rg` subcommands; regression works without it when `regr_weight` is embedded.
- [ ] `RegressionRunner` builds its dataset from `ldscore_table` — not from `reference_metadata` + `ld_scores` + `regression_metadata` + `w_ld`.

---

## 11. Public API surface (`__init__.py`)

- [ ] `from ldsc import AnnotationBuilder, AnnotationBundle, AnnotationSourceSpec` works.
- [ ] `from ldsc import LDScoreCalculator, LDScoreResult` works.
- [ ] `from ldsc import load_ldscore_from_files` works.
- [ ] `from ldsc import RegressionRunner` works.
- [ ] `from ldsc import run_bed_to_annot` works.
- [ ] No gene-set symbols are reachable from the top-level `ldsc` package.

---

## 12. Test coverage

- [ ] `AnnotationBuilder.run()` with `bed_paths` returns binary query columns and writes no `.annot.gz` files.
- [ ] `run_bed_to_annot()` returns `AnnotationBundle` and writes files only when `output_dir` is given.
- [ ] `LDScoreResult` / `ChromLDScoreResult` expose `ld_reference_snps` and `ld_regression_snps` as `frozenset`.
- [ ] A workflow test verifies that `compute_chromosome()` filters the chromosome-local annotation bundle to the restricted reference-panel metadata before calling the kernel.
- [ ] `build_parser()` rejects `--query-annot` + `--query-annot-bed` simultaneously.
- [ ] `run_ldscore_from_args()` result matches a direct `AnnotationBuilder` + `LDScoreCalculator` run on the same inputs (parity / regression test).
- [ ] `load_ldscore_from_files` is importable from `ldsc` top-level.
- [ ] New-format regression (embedded `regr_weight`, no `--w-ld`) produces consistent heritability estimates.
- [ ] Legacy-format regression (separate weight file via `--w-ld`) still works.
- [ ] `GlobalConfig` constructor raises `TypeError` or `AttributeError` when passed `ref_panel_snps_path` or `regression_snps_path`.
- [ ] `RefPanelSpec` accepts `ref_panel_snps_path` and `r2_bias_mode`; both default to `None`.

---

## Quick grep commands

Run these to catch stale references before marking the review complete:

```bash
# Old SNP universe field names on results
grep -rn "\.reference_snps\b\|\.regression_snps\b" src/ tests/

# Old split-table fields
grep -rn "reference_metadata\b\|regression_metadata\b\|\.ld_scores\b\|\.w_ld\b" src/ tests/

# GlobalConfig must not carry restriction paths
grep -rn "ref_panel_snps_path\|regression_snps_path" src/ldsc/config.py

# w.l2.ldscore.gz must not be written
grep -rn "w\.l2\.ldscore" src/

# Gene-set symbols must be gone
grep -rn "gene_set_to_bed\|make_annot_files\|main_make_annot\|parse_make_annot\|gene_set_paths" src/

# combine_annotation_groups must not be called from the workflow layer
grep -rn "combine_annotation_groups" src/ldsc/ldscore_calculator.py

# Dead helpers must be gone
grep -rn "_filter_annotation_bundle_to_reference_universe\|_chromosome_set_from_annotation_inputs" src/
```
