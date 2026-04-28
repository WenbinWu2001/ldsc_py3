# Partitioned-LDSC Pipeline Refactoring Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Refactor annotation processing and the partitioned-h2 pipeline to support in-memory BED → `AnnotationBundle` construction, a unified `run_ldscore_from_args` code path through `AnnotationBuilder` + `LDScoreCalculator`, a `--query-annot-bed` flag on `ldsc ldscore`, a single merged public `ldscore_table` result shape with embedded `regr_weight`, per-chromosome-only `.l2.ldscore.gz` output, and a cleaner config boundary where `ref_panel_snps_path` lives on `RefPanelConfig` and `regression_snps_path` lives on `LDScoreConfig` — both removed from `GlobalConfig`.

**Architecture:** Nine sequential tasks. Each builds on the previous: (1) remove gene-set TSV support; (2) normalize public `LDScoreResult` / `ChromLDScoreResult` to the merged `ldscore_table` shape; (3) restructure SNP-restriction config — move `ref_panel_snps_path` to `RefPanelConfig` and `regression_snps_path` to `LDScoreConfig`, remove both from `GlobalConfig`, add `r2_bias_mode` to `RefPanelConfig`; (4) add in-memory BED projection to `AnnotationBuilder.run()`; (5) refactor `project_bed_annotations` to return an `AnnotationBundle`; (6) converge `run_ldscore_from_args` through the workflow layer and make per-chromosome output the only output mode; (7) add `--query-annot-bed` to the ldscore CLI; (8) make `load_ldscore_from_files` public; (9) make `--w-ld` optional when `regr_weight` is embedded.

**Tech Stack:** Python 3.10+, pandas, numpy, pybedtools, argparse. Tests: `python -m unittest discover -s tests -p 'test*.py' -v`

---

## File Map

| File | Change type | Purpose |
|---|---|---|
| `src/ldsc/_kernel/annotation.py` | Modify (large) | Remove gene-set code; remove `_filter_metadata_by_global_restriction()`; add `_compute_bed_query_columns()`; modify `_run_single_universe()`; refactor `project_bed_annotations()`; add `_write_bundle_query_as_annot_files()` |
| `src/ldsc/annotation_builder.py` | Modify | Update `__all__` to remove gene-set exports |
| `src/ldsc/config.py` | Modify | Remove `ref_panel_snps_path` and `regression_snps_path` from `GlobalConfig`; add `regression_snps_path` to `LDScoreConfig` |
| `src/ldsc/_kernel/ref_panel.py` | Modify | Add `ref_panel_snps_path` and `r2_bias_mode` to `RefPanelConfig`; rename `_filter_metadata_by_global_restriction()` → `_apply_snp_restriction()` reading from spec |
| `src/ldsc/ldscore_calculator.py` | Modify (large) | Normalize public results to merged `ldscore_table` shape; add `_ref_panel_from_args()`, `_ldscore_config_from_args()`; rewrite `run_ldscore_from_args()`; enforce per-chromosome output; add `--query-annot-bed` (mutually exclusive with `--query-annot`) to `build_parser()`; remove `--ref-panel-snps` and `--regression-snps` from passing through `GlobalConfig` |
| `src/ldsc/regression_runner.py` | Modify | Consume `ldscore_table`; rename `_load_ldscore_result_from_files` → `load_ldscore_from_files` |
| `src/ldsc/__init__.py` | Modify | Export `load_ldscore_from_files` |
| `src/ldsc/cli.py` | Modify | Remove gene-set CLI arguments; simplify `_run_annotate()`; remove `--ref-panel-snps` and `--regression-snps` from non-ldscore subcommands |
| `tests/test_annotation.py` | Modify | Update for new API and BED in-memory behavior |
| `tests/test_ldscore_workflow.py` | Modify | Update for merged `ldscore_table` result shape; add GlobalConfig/LDScoreConfig tests |
| `tests/test_regression_workflow.py` | Modify | Update for merged `ldscore_table` result shape |
| `tests/test_ref_panel.py` | Modify | Test `r2_bias_mode` and `ref_panel_snps_path` fields |

---

## Task 1: Remove gene-set TSV support

**Files:**
- Modify: `src/ldsc/_kernel/annotation.py`
- Modify: `src/ldsc/annotation_builder.py`
- Modify: `src/ldsc/cli.py`
- Test: `tests/test_annotation.py`
- Test: `tests/test_package_layout.py`

- [ ] **Step 1: Write the failing test**

In `tests/test_annotation.py`, add:
```python
def test_annotation_source_spec_has_no_gene_set_paths(self):
    """AnnotationBuildConfig must not expose gene_set_paths after removal."""
    spec = AnnotationBuildConfig(baseline_annot_paths=("baseline.1.annot.gz",))
    self.assertFalse(hasattr(spec, "gene_set_paths"))

def test_gene_set_functions_not_exported(self):
    """Removed gene-set functions must not appear on the annotation_builder module."""
    import ldsc.annotation_builder as ab
    for name in ("gene_set_to_bed", "make_annot_files", "main_make_annot", "parse_make_annot_args"):
        self.assertFalse(hasattr(ab, name), f"{name} should not be exported")
```

- [ ] **Step 2: Run test to verify it fails**

```
python -m unittest tests.test_annotation.TestAnnotationBuildConfig.test_annotation_source_spec_has_no_gene_set_paths -v
```
Expected: FAIL (field still exists)

- [ ] **Step 3: Remove gene-set code from `_kernel/annotation.py`**

In `AnnotationBuildConfig`, remove the `gene_set_paths` field and its `__post_init__` line:
```python
# DELETE these two lines from AnnotationBuildConfig:
gene_set_paths: ... = field(default_factory=tuple)
# in __post_init__:
object.__setattr__(self, "gene_set_paths", _normalize_path_tuple(self.gene_set_paths))
```

In `_build_bundle()`, remove from `source_summary` dict:
```python
# DELETE:
"gene_set_paths": list(source_spec.gene_set_paths),
```

In `_run_sharded_inputs()`, remove from `chrom_source_spec` construction:
```python
# DELETE:
gene_set_paths=source_spec.gene_set_paths,
```

Remove the four module-level functions entirely:
- `gene_set_to_bed(args)` (lines ~747-759)
- `make_annot_files(args, bed_for_annot)` (lines ~762-765)
- `parse_make_annot_args(argv)` (lines ~813-833)
- `main_make_annot(argv)` (lines ~836-847)

- [ ] **Step 4: Update `annotation_builder.py` exports**

```python
# In src/ldsc/annotation_builder.py, change the import and __all__ to:
from ._kernel.annotation import (
    AnnotationBuilder,
    AnnotationBundle,
    AnnotationBuildConfig,
    main_bed_to_annot,
    parse_bed_to_annot_args,
    run_bed_to_annot,
)

__all__ = [
    "AnnotationBuilder",
    "AnnotationBundle",
    "AnnotationBuildConfig",
    "main_bed_to_annot",
    "parse_bed_to_annot_args",
    "run_bed_to_annot",
]
```

- [ ] **Step 5: Update `cli.py` to remove gene-set dispatch and arguments**

In `_run_annotate()`, remove the `main_make_annot` branch entirely. The function becomes:
```python
def _run_annotate(args: argparse.Namespace):
    """Dispatch the ``annotate`` subcommand."""
    if not args.bed_files:
        raise SystemExit("ldsc annotate requires --bed-files and --baseline-annot.")
    return annotation_builder.main_bed_to_annot(
        _namespace_to_argv(args, exclude={"command"})
    )
```

In `_add_annotate_arguments()`, remove these argument registrations:
- `--gene-set-file`
- `--gene-coord-file`
- `--windowsize`
- `--bed-file`
- `--nomerge`
- `--bimfile`
- `--annot-file`

- [ ] **Step 6: Run tests to verify they pass**

```
python -m unittest tests.test_annotation tests.test_package_layout -v
```
Expected: all PASS

- [ ] **Step 7: Commit**

```bash
git add src/ldsc/_kernel/annotation.py src/ldsc/annotation_builder.py src/ldsc/cli.py tests/test_annotation.py tests/test_package_layout.py
git commit -m "$(cat <<'EOF'
Remove gene-set TSV annotation support

BED files are now the only supported raw query annotation format.
Remove gene_set_to_bed, make_annot_files, main_make_annot,
parse_make_annot_args, and gene_set_paths from AnnotationBuildConfig.
Simplify ldsc annotate to require --bed-files.

Co-Authored-By: Claude Sonnet 4.6 <noreply@anthropic.com>
EOF
)"
```

---

## Task 2: Normalize public LDScoreResult / ChromLDScoreResult to merged ldscore_table shape

**Files:**
- Modify: `src/ldsc/ldscore_calculator.py`
- Modify: `src/ldsc/regression_runner.py`
- Modify: `tests/test_ldscore_workflow.py`
- Modify: `tests/test_regression_workflow.py`

- [ ] **Step 1: Write the failing test**

In `tests/test_ldscore_workflow.py`, add to an existing `LDScoreResult`-producing test:
```python
def test_ldscore_result_uses_single_table_shape(self):
    result = _build_minimal_ldscore_result()  # helper that returns an LDScoreResult
    self.assertTrue(hasattr(result, "ldscore_table"))
    self.assertTrue(hasattr(result, "ld_reference_snps"))
    self.assertTrue(hasattr(result, "ld_regression_snps"))
    # Split-table public fields must be gone
    self.assertFalse(hasattr(result, "reference_metadata"))
    self.assertFalse(hasattr(result, "regression_metadata"))
    self.assertFalse(hasattr(result, "w_ld"))
    # Public result mirrors written .l2.ldscore.gz rows
    self.assertIn("regr_weight", result.ldscore_table.columns)
    self.assertGreater(len(result.ldscore_table), 0)
    # ld_reference_snps is not recoverable from normalized row tables
    self.assertIsInstance(result.ld_reference_snps, frozenset)
    self.assertEqual(result.ld_reference_snps, frozenset())
    self.assertIsInstance(result.ld_regression_snps, frozenset)
```

- [ ] **Step 2: Run test to verify it fails**

```
python -m unittest tests.test_ldscore_workflow -v 2>&1 | head -30
```
Expected: FAIL (old split-table fields still present)

- [ ] **Step 3: Replace the split-table public dataclass shape**

In `src/ldsc/ldscore_calculator.py`, change both dataclass field definitions:
```python
# In ChromLDScoreResult, replace the public row tables:
ldscore_table: pd.DataFrame
ld_reference_snps: frozenset[str]
ld_regression_snps: frozenset[str]

# Remove public fields:
# - reference_metadata
# - regression_metadata
# - w_ld
#
# Same shape change in LDScoreResult. The public result object becomes a
# normalized/file-equivalent view:
# - ldscore_table: one DataFrame with [CHR, SNP, BP, <annot columns>, regr_weight]
# - snp_count_totals: counts already computed over the full ld_reference_snps
# - ld_reference_snps: frozenset() on normalized/public results
# - ld_regression_snps: reconstructed from ldscore_table rows
```

- [ ] **Step 4: Update `_wrap_legacy_chrom_result()` to build the merged normalized row table**

```python
def _wrap_legacy_chrom_result(self, legacy_result, global_config, regression_snps=None):
    reference_metadata = legacy_result.metadata.reset_index(drop=True).copy()
    ld_scores = pd.DataFrame(legacy_result.ld_scores, columns=list(legacy_result.ldscore_columns))
    internal_ld_reference_snps = frozenset(
        build_snp_id_series(reference_metadata, global_config.snp_identifier)
    )
    internal_ld_regression_snps = (
        internal_ld_reference_snps if regression_snps is None
        else frozenset(internal_ld_reference_snps.intersection(regression_snps))
    )
    regression_keep = build_snp_id_series(
        reference_metadata, global_config.snp_identifier
    ).isin(internal_ld_regression_snps)
    ldscore_table = pd.concat(
        [
            reference_metadata.loc[regression_keep, ["CHR", "SNP", "BP"]].reset_index(drop=True),
            ld_scores.loc[regression_keep].reset_index(drop=True),
            pd.DataFrame(
                {"regr_weight": np.asarray(legacy_result.w_ld, dtype=np.float32)[regression_keep.to_numpy()]}
            ).reset_index(drop=True),
        ],
        axis=1,
    )
    count_map = {"all_reference_snp_counts": np.asarray(legacy_result.M, dtype=np.float64)}
    if legacy_result.M_5_50 is not None:
        count_map["common_reference_snp_counts_maf_gt_0_05"] = np.asarray(legacy_result.M_5_50, dtype=np.float64)
    result = ChromLDScoreResult(
        chrom=str(legacy_result.chrom),
        ldscore_table=ldscore_table,
        snp_count_totals=count_map,
        baseline_columns=list(legacy_result.baseline_columns),
        query_columns=list(legacy_result.query_columns),
        ld_reference_snps=frozenset(),
        ld_regression_snps=frozenset(
            build_snp_id_series(ldscore_table, global_config.snp_identifier)
        ),
        output_paths={},
        config_snapshot=global_config,
    )
    result.validate()
    return result
```

- [ ] **Step 5: Update `_aggregate_chromosome_results()` to concatenate `ldscore_table`**

```python
ldscore_table=pd.concat(
    [result.ldscore_table for result in chromosome_results],
    ignore_index=True,
)
ld_reference_snps=frozenset(),
ld_regression_snps=frozenset().union(
    *(result.ld_regression_snps for result in chromosome_results)
),
```

- [ ] **Step 6: Update `_replace_result_output_paths()` to use new field names**

```python
# In _replace_result_output_paths(), change:
ldscore_table=result.ldscore_table,
ld_reference_snps=result.ld_reference_snps,
ld_regression_snps=result.ld_regression_snps,
```

- [ ] **Step 7: Update `regression_runner.py` consumers to use `ldscore_table`**

In `_subset_ldscore_result()`:
```python
# Subset ldscore_result.ldscore_table on the retained SNP rows and carry:
ldscore_table=subset_table.reset_index(drop=True),
ld_reference_snps=frozenset(),
ld_regression_snps=frozenset(ldscore_result.ld_regression_snps),
```

In `_load_ldscore_result_from_files()`:
```python
# Return:
ldscore_table=file_rows_df.reset_index(drop=True),
ld_reference_snps=frozenset(),
ld_regression_snps=frozenset(build_snp_id_series(file_rows_df, snp_identifier)),
```

- [ ] **Step 8: Update `summary()` / `validate()` methods for the new public shape**

Search and replace in `ldscore_calculator.py`:
```bash
grep -n "reference_metadata\|regression_metadata\|w_ld\|reference_snps\|regression_snps" src/ldsc/ldscore_calculator.py
```
Fix any remaining references (method bodies in `summary()`, `validate()`, etc.). Validation
should check the normalized/file-equivalent contract:
- `ldscore_table` contains `CHR`, `SNP`, `BP`, and `regr_weight`
- `ld_reference_snps == frozenset()`
- `ld_regression_snps` is derived from `ldscore_table` rows

- [ ] **Step 9: Update test files**

In `tests/test_ldscore_workflow.py` and `tests/test_regression_workflow.py`, replace all:
- `.reference_metadata` / `.regression_metadata` / `.w_ld` assertions with `ldscore_table`
- `.reference_snps` → `.ld_reference_snps`
- `.regression_snps` → `.ld_regression_snps`

```bash
grep -rn "\.reference_metadata\|\.regression_metadata\|\.w_ld\|\.reference_snps\|\.regression_snps" tests/
```

- [ ] **Step 10: Run tests**

```
python -m unittest tests.test_ldscore_workflow tests.test_regression_workflow -v
```
Expected: all PASS

- [ ] **Step 11: Commit**

```bash
git add src/ldsc/ldscore_calculator.py src/ldsc/regression_runner.py tests/test_ldscore_workflow.py tests/test_regression_workflow.py
git commit -m "$(cat <<'EOF'
Normalize public LDScoreResult to merged ldscore_table shape

Replace the public split-table shape (reference_metadata /
regression_metadata / w_ld) with a single normalized ldscore_table that
matches the written .l2.ldscore.gz rows and embeds regr_weight. Public
ld_reference_snps is now frozenset() because the full reference universe
is not recoverable from normalized row tables.

Co-Authored-By: Claude Sonnet 4.6 <noreply@anthropic.com>
EOF
)"
```

---

## Task 3: Restructure SNP-restriction config: move ref_panel_snps_path and regression_snps_path out of GlobalConfig

`ref_panel_snps_path` and `regression_snps_path` are LD-score-specific parameters with no meaning during annotation building, sumstats munging, or regression. This task moves them to the correct homes: `ref_panel_snps_path` → `RefPanelConfig` (the reference panel owns its own SNP restriction and applies it inside `load_metadata()`); `regression_snps_path` → `LDScoreConfig`; both removed from `GlobalConfig`. Also adds `r2_bias_mode` to `RefPanelConfig` (needed by Task 6).

**Key consequence for the SNP universe:** `AnnotationBuilder` no longer filters its annotation bundle by `ref_panel_snps_path`. The bundle keeps all annotation-file SNPs. The kernel intersection `annotation_snps ∩ reference_panel_snps` naturally yields the correct `ld_reference_snps` because the reference panel is already restricted by `RefPanelConfig.ref_panel_snps_path` inside `load_metadata()`.

**`--ref-panel-snps` and `--regression-snps` CLI flags** are only valid on `ldsc ldscore`. They must be removed from all other subcommand parsers (`ldsc annotate`, `ldsc h2`, `ldsc partitioned-h2`, `ldsc rg`, `ldsc munge-sumstats`).

**Files:**

- Modify: `src/ldsc/config.py`
- Modify: `src/ldsc/_kernel/ref_panel.py`
- Modify: `src/ldsc/_kernel/annotation.py`
- Modify: `tests/test_ref_panel.py`
- Modify: `tests/test_ldscore_workflow.py`

- [ ] **Step 1: Write the failing tests**

In `tests/test_ref_panel.py`:
```python
def test_ref_panel_spec_accepts_r2_bias_mode(self):
    from ldsc._kernel.ref_panel import RefPanelConfig
    spec = RefPanelConfig(backend="parquet_r2", r2_bias_mode="unbiased")
    self.assertEqual(spec.r2_bias_mode, "unbiased")

def test_ref_panel_spec_r2_bias_mode_defaults_to_none(self):
    from ldsc._kernel.ref_panel import RefPanelConfig
    spec = RefPanelConfig(backend="parquet_r2")
    self.assertIsNone(spec.r2_bias_mode)

def test_ref_panel_spec_accepts_ref_panel_snps_path(self):
    from ldsc._kernel.ref_panel import RefPanelConfig
    spec = RefPanelConfig(backend="plink", ref_panel_snps_path="/path/to/snps.txt")
    self.assertEqual(spec.ref_panel_snps_path, "/path/to/snps.txt")

def test_ref_panel_spec_ref_panel_snps_path_defaults_to_none(self):
    from ldsc._kernel.ref_panel import RefPanelConfig
    spec = RefPanelConfig(backend="plink")
    self.assertIsNone(spec.ref_panel_snps_path)
```

In `tests/test_ldscore_workflow.py`:
```python
def test_global_config_has_no_ref_panel_snps_path(self):
    from ldsc.config import GlobalConfig
    config = GlobalConfig(snp_identifier="rsid")
    self.assertFalse(hasattr(config, "ref_panel_snps_path"))

def test_global_config_has_no_regression_snps_path(self):
    from ldsc.config import GlobalConfig
    config = GlobalConfig(snp_identifier="rsid")
    self.assertFalse(hasattr(config, "regression_snps_path"))

def test_ldscore_config_accepts_regression_snps_path(self):
    from ldsc.config import LDScoreConfig
    config = LDScoreConfig(regression_snps_path="/path/to/snps.txt")
    self.assertEqual(config.regression_snps_path, "/path/to/snps.txt")
```

- [ ] **Step 2: Run tests to verify they fail**

```
python -m unittest tests.test_ref_panel tests.test_ldscore_workflow -v 2>&1 | head -30
```
Expected: FAIL (AttributeError on the new fields / field still present on GlobalConfig)

- [ ] **Step 3: Remove ref_panel_snps_path and regression_snps_path from GlobalConfig in config.py**

In `src/ldsc/config.py`, in the `GlobalConfig` dataclass, delete:
```python
# DELETE from GlobalConfig field definitions:
ref_panel_snps_path: str | PathLike[str] | None = None
regression_snps_path: str | PathLike[str] | None = None

# DELETE from GlobalConfig.__post_init__:
object.__setattr__(self, "ref_panel_snps_path", _normalize_optional_path(self.ref_panel_snps_path))
object.__setattr__(self, "regression_snps_path", _normalize_optional_path(self.regression_snps_path))
```

Also remove the comparison blocks for both fields from `_assert_global_configs_compatible()`:
```python
# DELETE the if-blocks that compare ref_panel_snps_path and regression_snps_path
```

- [ ] **Step 4: Add regression_snps_path to LDScoreConfig in config.py**

In the `LDScoreConfig` dataclass, add:
```python
regression_snps_path: str | PathLike[str] | None = None
```

Add to `LDScoreConfig.__post_init__` (create if not present):
```python
object.__setattr__(self, "regression_snps_path", _normalize_optional_path(self.regression_snps_path))
```

- [ ] **Step 5: Add ref_panel_snps_path and r2_bias_mode to RefPanelConfig; rename _filter_metadata_by_global_restriction()**

In `src/ldsc/_kernel/ref_panel.py`, add to `RefPanelConfig` after `genome_build`:
```python
r2_bias_mode: str | None = None
ref_panel_snps_path: str | PathLike[str] | None = None
```

Add to `RefPanelConfig.__post_init__`:
```python
object.__setattr__(self, "ref_panel_snps_path", _normalize_optional_path(self.ref_panel_snps_path))
```

Rename `_filter_metadata_by_global_restriction()` → `_apply_snp_restriction()` and change it to read from `self.spec` instead of `self.global_config`:
```python
def _apply_snp_restriction(self, metadata: pd.DataFrame) -> pd.DataFrame:
    """Filter metadata rows to RefPanelConfig.ref_panel_snps_path if set."""
    restrict_path = self.spec.ref_panel_snps_path
    if restrict_path is None or len(metadata) == 0:
        return metadata
    restrict_ids = read_global_snp_restriction(
        restrict_path,
        self.global_config.snp_identifier,
        genome_build=self.global_config.genome_build,
        logger=LOGGER,
    )
    keep = build_snp_id_series(metadata, self.global_config.snp_identifier).isin(restrict_ids)
    return metadata.loc[keep].reset_index(drop=True)
```

Update both call sites in `PlinkRefPanel.load_metadata()` and `ParquetR2RefPanel.load_metadata()`:

```python
# Change:
metadata = self._filter_metadata_by_global_restriction(metadata)
# To:
metadata = self._apply_snp_restriction(metadata)
```

- [ ] **Step 6: Remove _filter_metadata_by_global_restriction() from AnnotationBuilder**

In `src/ldsc/_kernel/annotation.py`, find and delete the `_filter_metadata_by_global_restriction()` method on `AnnotationBuilder` and every call to it:

```bash
grep -n "_filter_metadata_by_global_restriction" src/ldsc/_kernel/annotation.py
```
Delete the method body and all call sites. The annotation bundle now keeps all annotation-file SNPs; the reference panel's `_apply_snp_restriction()` ensures `ld_reference_snps` is correctly restricted.

- [ ] **Step 7: Run tests**

```
python -m unittest tests.test_ref_panel tests.test_ldscore_workflow tests.test_annotation -v
```
Expected: all PASS

- [ ] **Step 8: Commit**

```bash
git add src/ldsc/config.py src/ldsc/_kernel/ref_panel.py src/ldsc/_kernel/annotation.py tests/test_ref_panel.py tests/test_ldscore_workflow.py
git commit -m "$(cat <<'EOF'
Move ref_panel_snps_path and regression_snps_path out of GlobalConfig

ref_panel_snps_path moves to RefPanelConfig — the reference panel applies it
in _apply_snp_restriction() inside load_metadata(). regression_snps_path
moves to LDScoreConfig. Both removed from GlobalConfig. AnnotationBuilder
no longer filters its annotation bundle; the kernel intersection with the
already-restricted reference panel produces the correct ld_reference_snps.
Also adds r2_bias_mode and ref_panel_snps_path to RefPanelConfig.

Co-Authored-By: Claude Sonnet 4.6 <noreply@anthropic.com>
EOF
)"
```

---

## Task 4: Add in-memory BED projection to AnnotationBuilder.run()

This task adds `_compute_bed_query_columns()` and wires it into `_run_single_universe()`.
After this task, `AnnotationBuilder.run(AnnotationBuildConfig(baseline_annot_paths=..., bed_paths=[...]))` returns a bundle with BED-projected binary query columns. No intermediate `.annot.gz` files are written.

**Files:**
- Modify: `src/ldsc/_kernel/annotation.py`
- Modify: `tests/test_annotation.py`

- [ ] **Step 1: Write the failing test**

```python
def test_run_with_bed_paths_returns_bundle_with_binary_query_columns(self):
    """AnnotationBuilder.run() with bed_paths returns binary query annotation columns."""
    # Uses the test fixtures: baseline.1.annot.gz and a simple BED file
    baseline_path = str(FIXTURE_DIR / "baseline.1.annot.gz")
    bed_path = str(FIXTURE_DIR / "test_query.bed")
    config = GlobalConfig(snp_identifier="chr_pos")
    builder = AnnotationBuilder(config)
    spec = AnnotationBuildConfig(
        baseline_annot_paths=(baseline_path,),
        bed_paths=(bed_path,),
    )
    bundle = builder.run(spec)
    # Query columns are named by BED file stem
    self.assertIn("test_query", bundle.query_columns)
    # All query values are 0 or 1
    q = bundle.query_annotations["test_query"]
    self.assertTrue(q.isin([0.0, 1.0]).all())
    # Row counts must match
    self.assertEqual(len(bundle.metadata), len(bundle.query_annotations))

def test_run_with_bed_paths_no_disk_annot_written(self):
    """AnnotationBuilder.run() with bed_paths does not create .annot.gz files."""
    import os, tempfile
    with tempfile.TemporaryDirectory() as tmpdir:
        baseline_path = str(FIXTURE_DIR / "baseline.1.annot.gz")
        bed_path = str(FIXTURE_DIR / "test_query.bed")
        config = GlobalConfig(snp_identifier="chr_pos")
        builder = AnnotationBuilder(config)
        spec = AnnotationBuildConfig(
            baseline_annot_paths=(baseline_path,),
            bed_paths=(bed_path,),
        )
        builder.run(spec)
        # No .annot.gz files should have been written anywhere
        written = [f for f in os.listdir(tmpdir) if f.endswith(".annot.gz")]
        self.assertEqual(written, [])
```

Note: you may need to add a minimal BED fixture `tests/fixtures/test_query.bed` covering a few SNP positions from the baseline fixture. A file with content like:
```
chr1	10000	20000
chr1	50000	60000
```

- [ ] **Step 2: Run tests to verify they fail**

```
python -m unittest tests.test_annotation.TestAnnotationBuilderBED -v
```
Expected: FAIL (AttributeError or ValueError — bed_paths not yet handled by run())

- [ ] **Step 3: Add `_compute_bed_query_columns()` module-level function in `annotation.py`**

Add after `_compute_bed_overlap_mask()` (around line 1125):

```python
def _compute_bed_query_columns(
    metadata: pd.DataFrame,
    bed_paths: Sequence[Path],
    tempdir: Path,
) -> pd.DataFrame:
    """Compute binary BED-overlap annotation columns from SNP metadata in memory.

    One column per BED file is returned, named by BED file stem, aligned
    row-for-row to ``metadata``. Uses pybedtools temp files internally but
    does not write any LDSC-format ``.annot.gz`` intermediates.
    """
    pybedtools = _get_pybedtools()
    rows = [
        _BaselineRow(
            chrom=str(row.CHR),
            pos=int(row.POS),
            snp=str(row.SNP),
            cm=str(row.CM),
        )
        for row in metadata.itertuples(index=False)
    ]
    baseline_bed_path = _write_baseline_bed(rows, tempdir / "bed_query_baseline.bed")
    baseline_bed = pybedtools.BedTool(str(baseline_bed_path))
    annotation_names = [path.stem for path in bed_paths]
    masks = [_compute_bed_overlap_mask(rows, baseline_bed, bed_path) for bed_path in bed_paths]
    pybedtools.cleanup(remove_all=True)
    return pd.DataFrame(
        {name: np.array(mask, dtype=np.float32) for name, mask in zip(annotation_names, masks)},
        index=metadata.index,
    )
```

- [ ] **Step 4: Add BED projection block to `_run_single_universe()`**

In `AnnotationBuilder._run_single_universe()`, add the following block after the `if metadata is None: raise` check and **before** the `baseline_df = self._concat_or_empty(...)` line:

```python
        if metadata is None:
            raise ValueError("No annotation rows were loaded from the supplied sources.")

        # Project BED-file paths into binary query annotation columns, in memory.
        if source_spec.bed_paths:
            resolved_bed_paths = [
                Path(p) for p in resolve_file_group(source_spec.bed_paths, label="BED file")
            ]
            stems = [p.stem for p in resolved_bed_paths]
            duplicate_stems = sorted({s for s in stems if stems.count(s) > 1})
            if duplicate_stems:
                raise ValueError(
                    "BED basenames must be unique (they become annotation column names). "
                    f"Duplicates: {', '.join(duplicate_stems)}"
                )
            overlap_with_existing = sorted(set(stems) & seen_columns)
            if overlap_with_existing:
                raise ValueError(
                    "BED file stems clash with loaded annotation column names: "
                    f"{overlap_with_existing}"
                )
            with tempfile.TemporaryDirectory(prefix="bed2annot_") as tmpdir:
                tempdir_path = Path(tmpdir)
                normalized_beds: list[Path] = []
                for bed_path in resolved_bed_paths:
                    normalized_path = tempdir_path / bed_path.name
                    _write_normalized_bed(bed_path, normalized_path)
                    normalized_beds.append(normalized_path)
                bed_df = _compute_bed_query_columns(metadata, normalized_beds, tempdir_path)
            query_blocks.append(bed_df)
            query_columns.extend(bed_df.columns.tolist())
            seen_columns.update(bed_df.columns)

        baseline_df = self._concat_or_empty(baseline_blocks, index=metadata.index)
        query_df = self._concat_or_empty(query_blocks, index=metadata.index)
        # ... rest unchanged ...
```

- [ ] **Step 5: Run tests**

```
python -m unittest tests.test_annotation -v
```
Expected: all PASS

- [ ] **Step 6: Commit**

```bash
git add src/ldsc/_kernel/annotation.py tests/test_annotation.py tests/fixtures/
git commit -m "$(cat <<'EOF'
Add in-memory BED projection to AnnotationBuilder.run()

Add _compute_bed_query_columns() which converts SNP metadata rows into
binary annotation columns via pybedtools overlap, without writing any
LDSC-format .annot.gz intermediate files. Wire it into _run_single_universe()
so AnnotationBuildConfig.bed_paths is now fully handled by run().

Co-Authored-By: Claude Sonnet 4.6 <noreply@anthropic.com>
EOF
)"
```

---

## Task 5: Refactor project_bed_annotations to return AnnotationBundle

After this task: `run_bed_to_annot()` returns an `AnnotationBundle` instead of a `Path`. Writing `.annot.gz` files is optional via an optional `output_dir` argument.

**Scoping note:** `_write_bundle_query_as_annot_files()` is used exclusively by `project_bed_annotations()` (and transitively `run_bed_to_annot()`). It is **not** called from `run_ldscore_from_args()` or anywhere in the `ldsc ldscore` CLI path. Saving annotation files to disk is an explicit, opt-in choice for Python API callers who need the `.annot.gz` artifacts — it is not part of the normal LD score computation workflow. The `ldsc ldscore` CLI has no `--save-annot` flag; the LD score files written by `--out` are the only disk output of that command.

**Files:**
- Modify: `src/ldsc/_kernel/annotation.py`
- Modify: `tests/test_annotation.py`

- [ ] **Step 1: Write the failing test**

```python
def test_run_bed_to_annot_returns_annotation_bundle(self):
    """run_bed_to_annot() must return an AnnotationBundle, not a Path."""
    from ldsc.annotation_builder import run_bed_to_annot
    baseline_path = str(FIXTURE_DIR / "baseline.1.annot.gz")
    bed_path = str(FIXTURE_DIR / "test_query.bed")
    result = run_bed_to_annot(
        bed_files=[bed_path],
        baseline_annot_paths=[baseline_path],
    )
    self.assertIsInstance(result, AnnotationBundle)
    self.assertIn("test_query", result.query_columns)

def test_run_bed_to_annot_writes_files_when_output_dir_given(self):
    """run_bed_to_annot() with output_dir writes per-chromosome .annot.gz files."""
    import tempfile, os
    from ldsc.annotation_builder import run_bed_to_annot
    baseline_path = str(FIXTURE_DIR / "baseline.1.annot.gz")
    bed_path = str(FIXTURE_DIR / "test_query.bed")
    with tempfile.TemporaryDirectory() as tmpdir:
        result = run_bed_to_annot(
            bed_files=[bed_path],
            baseline_annot_paths=[baseline_path],
            output_dir=tmpdir,
        )
        self.assertIsInstance(result, AnnotationBundle)
        written = [f for f in os.listdir(tmpdir) if f.endswith(".annot.gz")]
        self.assertTrue(len(written) >= 1)
        # Each written file must contain the BED column name
        import gzip, pandas as pd
        for fname in written:
            df = pd.read_csv(os.path.join(tmpdir, fname), sep="\t", compression="gzip")
            self.assertIn("test_query", df.columns)
```

- [ ] **Step 2: Run tests to verify they fail**

```
python -m unittest tests.test_annotation.TestRunBedToAnnot -v
```
Expected: FAIL (returns Path, not AnnotationBundle)

- [ ] **Step 3: Add `_write_bundle_query_as_annot_files()` helper function in `annotation.py`**

Add after `_write_annot_file()` (around line 1136):

```python
def _write_bundle_query_as_annot_files(bundle: "AnnotationBundle", output_dir: Path) -> list[Path]:
    """Write query annotations from a bundle as per-chromosome .annot.gz files.

    Each file is named ``query.<chrom>.annot.gz`` and contains the metadata
    columns (CHR/POS/SNP/CM) followed by the query annotation columns.
    """
    output_paths = []
    for chrom in bundle.chromosomes:
        chrom_mask = bundle.metadata["CHR"].astype(str) == str(chrom)
        chrom_meta = bundle.metadata.loc[chrom_mask].reset_index(drop=True)
        chrom_query = bundle.query_annotations.loc[chrom_mask].reset_index(drop=True)
        out_path = output_dir / f"query.{chrom}.annot.gz"
        out_path.parent.mkdir(parents=True, exist_ok=True)
        combined = pd.concat([chrom_meta, chrom_query], axis=1)
        with gzip.open(out_path, "wt") as handle:
            combined.to_csv(handle, sep="\t", index=False)
        output_paths.append(out_path)
    return output_paths
```

- [ ] **Step 4: Refactor `project_bed_annotations()` to return `AnnotationBundle`**

Replace the current body of `project_bed_annotations()` with:

```python
def project_bed_annotations(
    self,
    bed_files: str | PathLike[str] | Sequence[str | PathLike[str]],
    baseline_annot_paths: str | PathLike[str] | Sequence[str | PathLike[str]],
    output_dir: str | Path | None = None,
    batch: bool | None = None,
    log_level: str | None = None,
) -> "AnnotationBundle":
    """Convert BED annotations into an in-memory AnnotationBundle.

    If ``output_dir`` is provided, also writes per-chromosome
    ``query.<chrom>.annot.gz`` files containing the projected annotation
    columns. The ``batch`` parameter is retained for backward compatibility
    but is ignored (all BED columns are combined in a single file per chromosome).
    """
    _configure_logging(log_level or self.global_config.log_level)
    source_spec = AnnotationBuildConfig(
        baseline_annot_paths=baseline_annot_paths,
        bed_paths=bed_files,
    )
    bundle = self.run(source_spec)
    if output_dir is not None:
        output_path = ensure_output_directory(output_dir, label="output directory")
        _write_bundle_query_as_annot_files(bundle, output_path)
    return bundle
```

- [ ] **Step 5: Update `run_bed_to_annot()` and `_run_bed_to_annot_with_global_config()`**

```python
def run_bed_to_annot(
    bed_files: str | PathLike[str] | Sequence[str | PathLike[str]],
    baseline_annot_paths: str | PathLike[str] | Sequence[str | PathLike[str]],
    output_dir: str | Path | None = None,
    batch: bool = True,
) -> "AnnotationBundle":
    """Build an in-memory AnnotationBundle from BED files and baseline annotations.

    If ``output_dir`` is provided, per-chromosome .annot.gz files are also
    written. Python workflows read identifier, genome-build, and SNP-restriction
    settings from the registered GlobalConfig.
    """
    return _run_bed_to_annot_with_global_config(
        bed_files=bed_files,
        baseline_annot_paths=baseline_annot_paths,
        output_dir=output_dir,
        batch=batch,
        global_config=get_global_config(),
        entrypoint="run_bed_to_annot",
    )


def _run_bed_to_annot_with_global_config(
    bed_files,
    baseline_annot_paths,
    output_dir,
    *,
    batch: bool,
    global_config: GlobalConfig,
    entrypoint: str,
) -> "AnnotationBundle":
    print_global_config_banner(entrypoint, global_config)
    build_config = AnnotationBuildConfig(query_bed_paths=bed_files, batch_mode=batch)
    return AnnotationBuilder(global_config, build_config).project_bed_annotations(
        bed_files=bed_files,
        baseline_annot_paths=baseline_annot_paths,
        output_dir=output_dir,
        batch=batch,
        log_level=global_config.log_level,
    )
```

- [ ] **Step 6: Run tests**

```
python -m unittest tests.test_annotation -v
```
Expected: all PASS

- [ ] **Step 7: Commit**

```bash
git add src/ldsc/_kernel/annotation.py tests/test_annotation.py
git commit -m "$(cat <<'EOF'
Refactor project_bed_annotations to return AnnotationBundle

run_bed_to_annot() and project_bed_annotations() now return an
AnnotationBundle. Disk writes are optional via an optional output_dir
argument. Callers that required a Path return value must update.
Add _write_bundle_query_as_annot_files() for the file-writing path.

Co-Authored-By: Claude Sonnet 4.6 <noreply@anthropic.com>
EOF
)"
```

---

## Task 6: Converge run_ldscore_from_args through AnnotationBuilder + LDScoreCalculator

After this task, `run_ldscore_from_args()` uses `AnnotationBuilder.run()` to build the annotation bundle and `LDScoreCalculator.run()` for computation, removing the inline `kernel_ldscore.combine_annotation_groups()` loop and the `_filter_annotation_bundle_to_reference_universe()` dead helper. The `ref_panel_snps_path` restriction is **not** applied inside `AnnotationBuilder` — it is applied by `RefPanel._apply_snp_restriction()` inside `load_metadata()` (Task 3). The annotation bundle keeps all annotation-file SNPs; the kernel intersection with the already-restricted reference panel produces the correct `ld_reference_snps`. This task also makes per-chromosome `.l2.ldscore.gz` output the only output mode; any aggregate-vs-per-chrom output switch is removed or ignored.

**Files:**
- Modify: `src/ldsc/ldscore_calculator.py`
- Modify: `tests/test_ldscore_workflow.py`

- [ ] **Step 1: Write regression test**

In `tests/test_ldscore_workflow.py`, add a parity test that runs both the old and new code path on the same PLINK fixtures and checks that key numerical results are identical:

```python
def test_run_ldscore_from_args_produces_same_result_as_calculator_run(self):
    """run_ldscore_from_args result must match a direct AnnotationBuilder + LDScoreCalculator run."""
    from ldsc.ldscore_calculator import run_ldscore_from_args, LDScoreCalculator
    from ldsc._kernel.annotation import AnnotationBuilder, AnnotationBuildConfig
    from ldsc._kernel.ref_panel import PlinkRefPanel, RefPanelConfig
    from ldsc.config import GlobalConfig, LDScoreConfig
    import argparse, numpy as np

    baseline = str(FIXTURE_DIR / "baseline.@.annot.gz")
    bfile = str(FIXTURE_DIR / "plink.@")  # adjust to actual fixture prefix
    out_prefix = str(FIXTURE_DIR / "test_parity_out")
    config = GlobalConfig(snp_identifier="rsid")

    args = argparse.Namespace(
        baseline_annot=baseline, query_annot=None, query_annot_bed=None,
        bfile=bfile, r2_table=None, out=out_prefix,
        snp_identifier="rsid", genome_build=None,
        ld_wind_snps=None, ld_wind_kb=1.0, ld_wind_cm=None,
        maf=None, chunk_size=50,
        yes_really=False, r2_bias_mode=None, r2_sample_size=None,
        ref_panel_snps_path=None, regression_snps_path=None,
        frqfile=None, keep=None, log_level="WARNING",
    )
    result = run_ldscore_from_args(args)
    self.assertIsNotNone(result)
    self.assertGreater(len(result.ldscore_table), 0)
    self.assertIn("regr_weight", result.ldscore_table.columns)
```

- [ ] **Step 2: Run test to verify it passes with old code (baseline)**

```
python -m unittest tests.test_ldscore_workflow.TestRunLdscoreFromArgs.test_run_ldscore_from_args_produces_same_result_as_calculator_run -v
```
Expected: PASS (verify it passes NOW before the rewrite)

- [ ] **Step 3: Add `_ref_panel_from_args()` and `_ldscore_config_from_args()` helpers**

Add both below `_load_regression_snps()` in `ldscore_calculator.py`:

```python
def _ref_panel_from_args(args: argparse.Namespace, global_config: GlobalConfig):
    """Build a RefPanel adapter from normalized CLI arguments."""
    from ._kernel.ref_panel import PlinkRefPanel, ParquetR2RefPanel, RefPanelConfig
    freq_tokens = split_cli_path_tokens(getattr(args, "frqfile", None))
    ref_panel_snps_path = getattr(args, "ref_panel_snps_path", None)
    if getattr(args, "r2_table", None) is not None:
        r2_tokens = split_cli_path_tokens(args.r2_table)
        spec = RefPanelConfig(
            backend="parquet_r2",
            r2_paths=tuple(r2_tokens) if r2_tokens else (),
            metadata_paths=tuple(freq_tokens) if freq_tokens else (),
            r2_bias_mode=getattr(args, "r2_bias_mode", None),
            sample_size=getattr(args, "r2_sample_size", None),
            genome_build=global_config.genome_build,
            ref_panel_snps_path=ref_panel_snps_path,
        )
        return ParquetR2RefPanel(global_config, spec)
    spec = RefPanelConfig(
        backend="plink",
        bfile_prefix=getattr(args, "bfile", None),
        metadata_paths=tuple(freq_tokens) if freq_tokens else (),
        genome_build=global_config.genome_build,
        ref_panel_snps_path=ref_panel_snps_path,
    )
    return PlinkRefPanel(global_config, spec)


def _ldscore_config_from_args(args: argparse.Namespace) -> LDScoreConfig:
    """Build an LDScoreConfig from normalized CLI arguments."""
    return LDScoreConfig(
        ld_wind_snps=getattr(args, "ld_wind_snps", None),
        ld_wind_kb=getattr(args, "ld_wind_kb", None),
        ld_wind_cm=getattr(args, "ld_wind_cm", None),
        chunk_size=getattr(args, "chunk_size", 50),
        maf_min=getattr(args, "maf", None),
        whole_chromosome_ok=getattr(args, "yes_really", False),
        keep_individuals_path=getattr(args, "keep", None),
        regression_snps_path=getattr(args, "regression_snps_path", None),
    )
```

- [ ] **Step 4: Rewrite `run_ldscore_from_args()`**

Replace the body of `run_ldscore_from_args()` with:

```python
def run_ldscore_from_args(args: argparse.Namespace) -> LDScoreResult:
    """Run LD-score calculation from a parsed CLI namespace.

    Resolves annotation inputs through AnnotationBuilder so that BED files,
    glob patterns, and @ chromosome-suite tokens are all handled uniformly
    before computation.
    """
    from ._kernel.annotation import AnnotationBuilder, AnnotationBuildConfig
    normalized_args, global_config = _normalize_run_args(args)
    print_global_config_banner("run_ldscore_from_args", global_config)
    kernel_ldscore.validate_args(normalized_args)

    # regression_snps_path lives on LDScoreConfig, not GlobalConfig.
    ldscore_config = _ldscore_config_from_args(normalized_args)
    regression_snps = _load_regression_snps(ldscore_config.regression_snps_path, global_config)

    baseline_tokens = split_cli_path_tokens(normalized_args.baseline_annot)
    query_tokens = split_cli_path_tokens(normalized_args.query_annot)
    bed_tokens = split_cli_path_tokens(getattr(normalized_args, "query_annot_bed", None))

    source_spec = AnnotationBuildConfig(
        baseline_annot_paths=tuple(baseline_tokens) if baseline_tokens else (),
        query_annot_paths=tuple(query_tokens) if query_tokens else (),
        bed_paths=tuple(bed_tokens) if bed_tokens else (),
    )
    annotation_bundle = AnnotationBuilder(global_config).run(source_spec)
    # ref_panel_snps_path is passed via RefPanelConfig; no separate filtering step needed.
    ref_panel = _ref_panel_from_args(normalized_args, global_config)
    output_spec = _output_spec_from_args(normalized_args)

    calculator = LDScoreCalculator()
    return calculator.run(
        annotation_bundle=annotation_bundle,
        ref_panel=ref_panel,
        ldscore_config=ldscore_config,
        global_config=global_config,
        output_spec=output_spec,
        regression_snps=regression_snps,
    )
```

- [ ] **Step 5: Remove `_filter_annotation_bundle_to_reference_universe()` and `_chromosome_set_from_annotation_inputs()` since they are no longer called**

These two private helpers in `ldscore_calculator.py` are dead code after the rewrite. Delete them (lines ~792-832). Verify with:
```bash
grep -n "_filter_annotation_bundle_to_reference_universe\|_chromosome_set_from_annotation_inputs" src/ldsc/ldscore_calculator.py
```

- [ ] **Step 6: Update the output writer in `LDScoreCalculator.run()` to produce the merged format**

**Output format contract:**
- Rows = `ld_regression_snps` (B ∩ A' ∩ C). When `regression_snps_path` is absent, C = A', so rows = `ld_reference_snps`.
- Columns = `[CHR, SNP, BP, <baseline_L2 cols>, <query_L2 cols>, regr_weight]`
- `regr_weight` is the per-SNP total LD score used as regression weight (the value that was previously written to `.w.l2.ldscore.gz`)
- No separate `.w.l2.ldscore.gz` file is emitted
- `.M` and `.M_5_50` count files are unchanged: they track `ld_reference_snps` (B ∩ A') per annotation column
- Public `LDScoreResult` / `ChromLDScoreResult` objects are normalized to this same merged row-table shape
- Per-chromosome output is the only output mode; write `<out_prefix>.<chrom>.l2.ldscore.gz` for each chromosome

In `LDScoreCalculator.run()` (or its per-chromosome output assembly), replace the two-file write with a single-file write that appends `regr_weight` as the last column before writing:

```python
# After computing per-chromosome ld_scores (rows = ld_reference_snps) and
# regression_weights (rows = ld_regression_snps):
# Filter ld_scores to regression rows and append weight column.
regression_mask = build_snp_id_series(reference_metadata, global_config.snp_identifier).isin(
    ld_regression_snps
)
output_df = ld_scores.loc[regression_mask].copy()
output_df["regr_weight"] = regression_weights.values
# Write output_df to <out_prefix>.<chrom>.l2.ldscore.gz
```

Add a test:
```python
def test_ldscore_output_has_regr_weight_column_and_no_w_ld_file(self):
    import os, pandas as pd
    # Run ldsc ldscore on fixture data and verify output
    result = run_ldscore_from_args(args)
    # Output file has regr_weight column
    out_path = str(FIXTURE_DIR / "test_parity_out.1.l2.ldscore.gz")
    df = pd.read_csv(out_path, sep="\t", compression="gzip")
    self.assertIn("regr_weight", df.columns)
    # No separate w_ld file written
    w_ld_path = str(FIXTURE_DIR / "test_parity_out.1.w.l2.ldscore.gz")
    self.assertFalse(os.path.exists(w_ld_path))
    # Row count matches ld_regression_snps
    self.assertEqual(len(df), len(result.ld_regression_snps))
```

- [ ] **Step 7: Run tests**

```
python -m unittest tests.test_ldscore_workflow -v
```
Expected: all PASS, including the parity test from Step 1 and the output-format test from Step 6.

- [ ] **Step 8: Commit**

```bash
git add src/ldsc/ldscore_calculator.py tests/test_ldscore_workflow.py
git commit -m "$(cat <<'EOF'
Converge run_ldscore_from_args; merge regr_weight into output ldscore file

Replace the per-chromosome combine_annotation_groups + direct kernel call
loop with AnnotationBuilder.run() + LDScoreCalculator.run(). Add
_ref_panel_from_args() and _ldscore_config_from_args() helpers. Remove
_filter_annotation_bundle_to_reference_universe() and
_chromosome_set_from_annotation_inputs() (dead code). Output .l2.ldscore.gz
now has rows = ld_regression_snps and a regr_weight column; no separate
.w.l2.ldscore.gz is written.

Co-Authored-By: Claude Sonnet 4.6 <noreply@anthropic.com>
EOF
)"
```

---

## Task 7: Add --query-annot-bed to ldsc ldscore CLI; enforce mutual exclusivity with --query-annot

After this task, `ldsc ldscore --query-annot-bed my_annotation.bed ...` is a valid command that runs in-memory BED projection and LD score computation in a single step. `--query-annot-bed` and `--query-annot` are mutually exclusive: the former accepts raw BED files and projects them in memory; the latter accepts pre-built `.annot.gz` files. Supplying both at once is an error enforced at parse time by argparse.

**Files:**
- Modify: `src/ldsc/ldscore_calculator.py` (`build_parser()`)
- Modify: `tests/test_ldscore_workflow.py`

- [ ] **Step 1: Write the failing tests**

```python
def test_build_parser_has_query_annot_bed_argument(self):
    from ldsc.ldscore_calculator import build_parser
    parser = build_parser()
    args = parser.parse_args(["--out", "test", "--query-annot-bed", "my.bed"])
    self.assertEqual(args.query_annot_bed, "my.bed")

def test_build_parser_query_annot_and_query_annot_bed_are_mutually_exclusive(self):
    from ldsc.ldscore_calculator import build_parser
    parser = build_parser()
    with self.assertRaises(SystemExit):
        parser.parse_args([
            "--out", "test",
            "--query-annot", "q.1.annot.gz",
            "--query-annot-bed", "my.bed",
        ])
```

- [ ] **Step 2: Run tests to verify they fail**

```
python -m unittest tests.test_ldscore_workflow.TestBuildParser.test_build_parser_has_query_annot_bed_argument tests.test_ldscore_workflow.TestBuildParser.test_build_parser_query_annot_and_query_annot_bed_are_mutually_exclusive -v
```
Expected: FAIL (unrecognized argument)

- [ ] **Step 3: Replace the standalone `--query-annot` `add_argument` call with a `mutually_exclusive_group` in `build_parser()`**

Find the existing `parser.add_argument("--query-annot", ...)` line and replace it (and only it) with:
```python
_query_group = parser.add_mutually_exclusive_group()
_query_group.add_argument(
    "--query-annot",
    default=None,
    help="Pre-built query annotation files (.annot[.gz]). Comma-separated path tokens.",
)
_query_group.add_argument(
    "--query-annot-bed",
    default=None,
    help=(
        "Comma-separated BED file path tokens (exact paths or globs) to project "
        "in memory as binary query annotations. No intermediate .annot.gz files "
        "are written. Requires --baseline-annot. Mutually exclusive with --query-annot."
    ),
)
```

Also add `query_annot_bed` to `_normalize_run_args()` default handling so callers that build an `argparse.Namespace` by hand are not required to include the field:
```python
# In _normalize_run_args(), add to the attribute-defaulting block:
if not hasattr(normalized_args, "query_annot_bed"):
    setattr(normalized_args, "query_annot_bed", None)
normalized_args.query_annot_bed = getattr(args, "query_annot_bed", None)
```

The attribute is already consumed by `run_ldscore_from_args()` via `getattr(normalized_args, "query_annot_bed", None)` (wired in Task 6).

- [ ] **Step 4: Run tests**

```
python -m unittest tests.test_ldscore_workflow -v
```
Expected: all PASS

- [ ] **Step 5: Commit**

```bash
git add src/ldsc/ldscore_calculator.py tests/test_ldscore_workflow.py
git commit -m "$(cat <<'EOF'
Add --query-annot-bed flag to ldsc ldscore; enforce mutual exclusivity with --query-annot

Users can now pass raw BED files directly to ldsc ldscore without a
separate ldsc annotate step. BED projection runs in memory via
AnnotationBuilder.run() before LD score computation. --query-annot-bed
and --query-annot are in a mutually_exclusive_group so argparse rejects
both being set at the same time.

Co-Authored-By: Claude Sonnet 4.6 <noreply@anthropic.com>
EOF
)"
```

---

## Task 8: Make load_ldscore_from_files public

**Files:**
- Modify: `src/ldsc/regression_runner.py`
- Modify: `src/ldsc/__init__.py`
- Modify: `tests/test_regression_workflow.py`

- [ ] **Step 1: Write the failing tests**

```python
def test_load_ldscore_from_files_is_public(self):
    """load_ldscore_from_files must be importable from the top-level ldsc package."""
    from ldsc import load_ldscore_from_files
    self.assertTrue(callable(load_ldscore_from_files))

def test_load_ldscore_from_files_new_format_no_weight_path(self):
    """New merged format: regr_weight column present; weight_path omitted."""
    from ldsc import load_ldscore_from_files
    # Fixture produced by the new writer: has regr_weight column, rows = ld_regression_snps
    ldscore_path = str(FIXTURE_DIR / "test_merged.l2.ldscore.gz")
    m_path = str(FIXTURE_DIR / "test.l2.M_5_50")
    result = load_ldscore_from_files(
        ldscore_path=ldscore_path,
        counts_path=m_path,
        count_kind="m_5_50",
        snp_identifier="rsid",
    )
    self.assertIsNotNone(result)
    self.assertGreater(len(result.ldscore_table), 0)
    self.assertIn("regr_weight", result.ldscore_table.columns)
    # ld_reference_snps is not recoverable from disk; must be empty frozenset
    self.assertEqual(result.ld_reference_snps, frozenset())
    # ld_regression_snps reconstructed from file rows via build_snp_id_series
    self.assertIsInstance(result.ld_regression_snps, frozenset)
    self.assertGreater(len(result.ld_regression_snps), 0)

def test_load_ldscore_from_files_legacy_format_separate_weight_file(self):
    """Legacy format: no regr_weight column; weights loaded from separate weight_path."""
    from ldsc import load_ldscore_from_files
    ldscore_path = str(FIXTURE_DIR / "test_legacy.l2.ldscore.gz")
    w_ld_path = str(FIXTURE_DIR / "test.w.l2.ldscore.gz")
    m_path = str(FIXTURE_DIR / "test.l2.M_5_50")
    result = load_ldscore_from_files(
        ldscore_path=ldscore_path,
        counts_path=m_path,
        count_kind="m_5_50",
        weight_path=w_ld_path,
        snp_identifier="rsid",
    )
    self.assertIsNotNone(result)
    self.assertGreater(len(result.ldscore_table), 0)
    self.assertIn("regr_weight", result.ldscore_table.columns)
```

- [ ] **Step 2: Run test to verify it fails**

```
python -m unittest tests.test_regression_workflow.TestLoadLdscoreFromFiles -v
```
Expected: FAIL (ImportError — function not yet exported)

- [ ] **Step 3: Rename the function in `regression_runner.py` and update its signature**

Rename `_load_ldscore_result_from_files` → `load_ldscore_from_files` (remove the leading `_`) and make `weight_path` optional:

```python
def load_ldscore_from_files(
    ldscore_path: str,
    counts_path: str,
    count_kind: str,
    weight_path: str | None = None,
    annotation_manifest_path: str | None = None,
    explicit_query_columns: str | None = None,
    snp_identifier: str = "rsid",
) -> LDScoreResult:
    """Rebuild an LDScoreResult from previously written LD-score artifacts.

    Supports two file formats:
    - New merged format: `.l2.ldscore.gz` has a `regr_weight` column; `weight_path` omitted.
    - Legacy format: weights in a separate `.w.l2.ldscore.gz`; pass via `weight_path`.
    """
```

**SNP-set reconstruction contract:**

```python
df = pd.read_csv(ldscore_path, sep="\t", compression="gzip")
if "regr_weight" in df.columns:
    regression_weights = df["regr_weight"].values
    file_rows_df = df
elif weight_path is not None:
    w_df = pd.read_csv(weight_path, sep="\t", compression="gzip")
    regression_weights = w_df["L2"].values
    file_rows_df = df
else:
    regression_weights = None
    file_rows_df = df

# file_rows_df is the normalized/public row table:
# - new format: df already includes regr_weight
# - legacy format: append regr_weight from the separate weight file
if "regr_weight" not in file_rows_df.columns and regression_weights is not None:
    file_rows_df = file_rows_df.assign(regr_weight=regression_weights)

# ld_reference_snps is not recoverable from disk: column sums (M/M_5_50) were
# computed over ld_reference_snps at write time, but the file rows are
# ld_regression_snps (a strict subset). Setting to empty frozenset is correct;
# regression never reads ld_reference_snps after M/M_5_50 are loaded.
ld_reference_snps = frozenset()

# ld_regression_snps: reconstruct from normalized row table using the
# snp_identifier mode so that chr_pos vs rsid is applied consistently.
ld_regression_snps = frozenset(build_snp_id_series(file_rows_df, snp_identifier))
```

Update the three callers inside `regression_runner.py` that reference the old name:

```python
# In run_h2_from_args, run_partitioned_h2_from_args, run_rg_from_args:
# Change _load_ldscore_result_from_files(...) → load_ldscore_from_files(...)
# Pass snp_identifier=global_config.snp_identifier
```

- [ ] **Step 4: Export from `src/ldsc/__init__.py`**

Add to the imports and `__all__` in `__init__.py`:
```python
from .regression_runner import load_ldscore_from_files
# and add "load_ldscore_from_files" to __all__
```

- [ ] **Step 5: Run all tests**

```
python -m unittest discover -s tests -p 'test*.py' -v
```
Expected: all PASS

- [ ] **Step 6: Commit**

```bash
git add src/ldsc/regression_runner.py src/ldsc/__init__.py tests/test_regression_workflow.py
git commit -m "$(cat <<'EOF'
Make load_ldscore_from_files a public API function

Rename _load_ldscore_result_from_files to load_ldscore_from_files and
export it from the top-level ldsc package. Enables Python callers to
reconstruct an LDScoreResult from disk artifacts and pass it directly
to RegressionRunner without re-running LD score computation.

Co-Authored-By: Claude Sonnet 4.6 <noreply@anthropic.com>
EOF
)"
```

---

## Task 9: Make --w-ld optional on regression subcommands; fall back to regr_weight column

After this task, `ldsc partitioned-h2` (and `h2`, `rg`) no longer require `--w-ld` when the
input `.l2.ldscore.gz` was produced by the new writer with an embedded `regr_weight` column.
`--w-ld` remains a valid optional argument for backward compatibility with pre-computed weight
files (e.g., EUR `w_ld_chr/` downloads). When `--w-ld` is absent the regression runner detects
the `regr_weight` column and uses it as the per-SNP weight.

**Files:**

- Modify: `src/ldsc/regression_runner.py`
- Modify: `tests/test_regression_workflow.py`

- [ ] **Step 1: Write the failing test**

```python
def test_run_partitioned_h2_from_args_without_w_ld_uses_regr_weight(self):
    """regression subcommand works when --w-ld is absent and ldscore file has regr_weight."""
    import argparse
    from ldsc.regression_runner import run_partitioned_h2_from_args
    args = argparse.Namespace(
        ldscore_path=str(FIXTURE_DIR / "test_merged.l2.ldscore.gz"),
        w_ld=None,          # absent — must fall back to regr_weight column
        sumstats=str(FIXTURE_DIR / "test.sumstats.gz"),
        ref_ld_chr=str(FIXTURE_DIR / "test_merged.@"),
        # ...other required args...
    )
    result = run_partitioned_h2_from_args(args)
    self.assertIsNotNone(result)
```

- [ ] **Step 2: Run test to verify it fails**

```
python -m unittest tests.test_regression_workflow.TestRunPartitionedH2FromArgs.test_run_partitioned_h2_from_args_without_w_ld_uses_regr_weight -v
```

Expected: FAIL (error about missing weights)

- [ ] **Step 3: Add `_resolve_regression_weights()` helper in `regression_runner.py`**

Add this helper above `run_h2_from_args`:

```python
def _resolve_regression_weights(
    ldscore_result: LDScoreResult,
    w_ld_path: str | None,
    global_config: GlobalConfig,
) -> pd.DataFrame:
    """Return a regression-weight DataFrame from either a separate file or the embedded column.

    Priority:
    1. If w_ld_path is provided, load from file (backward compatibility).
    2. If ldscore_result has a 'regr_weight' column in ldscore_table, extract it.
    3. Raise ValueError — weights are required for regression.
    """
    if w_ld_path is not None:
        return pd.read_csv(w_ld_path, sep="\t", compression="infer")
    table = ldscore_result.ldscore_table
    if "regr_weight" in table.columns:
        return table[["CHR", "SNP", "BP", "regr_weight"]].rename(
            columns={"regr_weight": "L2"}
        )
    raise ValueError(
        "Regression weights unavailable: provide --w-ld or use an .l2.ldscore.gz "
        "produced by this tool (which embeds a regr_weight column)."
    )
```

- [ ] **Step 4: Wire `_resolve_regression_weights()` into each regression subcommand**

In `run_h2_from_args`, `run_partitioned_h2_from_args`, and `run_rg_from_args`, replace
the existing unconditional weight-file load with:

```python
# Before (remove this):
w_ld_df = pd.read_csv(args.w_ld, sep="\t", compression="infer")

# After (in each function, after loading ldscore_result):
w_ld_df = _resolve_regression_weights(ldscore_result, getattr(args, "w_ld", None), global_config)
```

`--w-ld` is already optional (default `None`) in `cli.py` for all three subcommands; no CLI
change is needed if the argument already exists. Confirm with:

```bash
grep -n "w.ld" src/ldsc/cli.py
```

If `--w-ld` has `required=True`, change it to `default=None`.

- [ ] **Step 5: Run tests**

```
python -m unittest tests.test_regression_workflow -v
```

Expected: all PASS, including the new no-`--w-ld` test and any existing `--w-ld` tests.

- [ ] **Step 6: Commit**

```bash
git add src/ldsc/regression_runner.py tests/test_regression_workflow.py
git commit -m "$(cat <<'EOF'
Make --w-ld optional on regression subcommands; fall back to regr_weight column

When the input .l2.ldscore.gz has a regr_weight column (new merged format),
regression subcommands no longer require --w-ld. _resolve_regression_weights()
checks for the embedded column before raising an error. --w-ld remains valid
for backward compatibility with pre-computed EUR w_ld files.

Co-Authored-By: Claude Sonnet 4.6 <noreply@anthropic.com>
EOF
)"
```

---

## Final Verification

- [ ] **Run the full test suite**

```
python -m unittest discover -s tests -p 'test*.py' -v 2>&1 | tail -20
```
Expected: all tests PASS, 0 errors.

- [ ] **Smoke-test the CLI**

```bash
python -m ldsc --help
python -m ldsc annotate --help
python -m ldsc ldscore --help | grep -A2 "query-annot-bed"
python -m ldsc partitioned-h2 --help
```

- [ ] **Verify public API surface**

```python
python -c "
from ldsc import (
    AnnotationBuilder, AnnotationBundle, AnnotationBuildConfig,
    LDScoreCalculator, LDScoreResult,
    RegressionRunner, load_ldscore_from_files,
    run_bed_to_annot,
)
print('Public API OK')
"
```

- [ ] **Update PLANS.md**

Add a completed section for this refactoring under "Completed" in `PLANS.md`.

---

## Self-Review Checklist

**Spec coverage:**
- [x] Q1 — `project_bed_annotations` / `run_bed_to_annot` return `AnnotationBundle`, `output_dir` optional (Tasks 4–5)
- [x] Q2 — Public `LDScoreResult` / `ChromLDScoreResult` normalized to one merged `ldscore_table`; `ld_reference_snps = frozenset()` on normalized/public results (Task 2)
- [x] Q3 — Resumability already present via `_load_ldscore_result_from_files`; now public with weight_path optional and SNP reconstruction via `build_snp_id_series` (Task 8)
- [x] Q4 — Python API and CLI changes: in-memory BED path, `--query-annot-bed` flag mutually exclusive with `--query-annot` (Tasks 4–5, 7)
- [x] Q5 — Gene-set TSV removed; BED only (Task 1)
- [x] Q6 — `run_ldscore_from_args` converged through `AnnotationBuilder` + `LDScoreCalculator` (Task 6)
- [x] Merged output format — `regr_weight` embedded in `.l2.ldscore.gz`; rows = `ld_regression_snps`; no `.w.l2.ldscore.gz` emitted (Task 6)
- [x] `--w-ld` optional on regression subcommands — falls back to `regr_weight` column; backward compatible with pre-computed w_ld files (Task 9)
- [x] Config restructure — `ref_panel_snps_path` on `RefPanelConfig`; `regression_snps_path` on `LDScoreConfig`; both removed from `GlobalConfig`; `--ref-panel-snps` and `--regression-snps` only on `ldsc ldscore` (Task 3)

**Config restructure scope (Task 3) — all callers updated:**

- `GlobalConfig`: fields `ref_panel_snps_path` and `regression_snps_path` deleted
- `LDScoreConfig`: field `regression_snps_path` added
- `RefPanelConfig`: fields `ref_panel_snps_path` and `r2_bias_mode` added
- `RefPanel._filter_metadata_by_global_restriction()` → `_apply_snp_restriction()` reading from `self.spec`
- `AnnotationBuilder._filter_metadata_by_global_restriction()` removed entirely
- `_ref_panel_from_args()`: passes `ref_panel_snps_path` to `RefPanelConfig`
- `_ldscore_config_from_args()`: passes `regression_snps_path` to `LDScoreConfig`
- `run_ldscore_from_args()`: reads regression SNPs from `ldscore_config`, not `global_config`

**Public result-shape scope (Task 2) — all callers updated:**
- `ChromLDScoreResult.ldscore_table`, `ld_reference_snps`, `ld_regression_snps`
- `LDScoreResult.ldscore_table`, `ld_reference_snps`, `ld_regression_snps`
- `_wrap_legacy_chrom_result()` normalizes to file-equivalent rows with embedded `regr_weight`
- `_aggregate_chromosome_results()`
- `_replace_result_output_paths()`
- `_subset_ldscore_result()` in regression_runner.py
- `_load_ldscore_result_from_files()` → `load_ldscore_from_files()`
- All test files

**BED in-memory — no intermediate .annot.gz files:**
- `_compute_bed_query_columns()` uses pybedtools temp files (internal) but returns a DataFrame
- `project_bed_annotations()` calls `self.run(source_spec)` which calls `_run_single_universe()` which calls `_compute_bed_query_columns()`
- File writing only when `output_dir` is explicitly provided

**Dependency order confirmed:**
Task 1 (gene-set removal) → Task 4 (BED in run()) since `_run_sharded_inputs` no longer passes `gene_set_paths`
Task 3 (RefPanelConfig.r2_bias_mode) → Task 6 (convergence uses RefPanelConfig)
Task 4+5 (BED bundle) → Task 6 (run_ldscore_from_args uses AnnotationBuilder with bed_paths)
Task 6 (convergence) → Task 7 (--query-annot-bed consumed by run_ldscore_from_args; mutual exclusivity with --query-annot enforced at parse time)
