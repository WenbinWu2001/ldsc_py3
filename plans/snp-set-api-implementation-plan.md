# Implementation Plan: Two-Field SNP-Set API

## Background and Motivation

`GlobalConfig.restrict_snps_path` currently overloads two scientifically distinct
concepts:

1. **Reference-panel / annotation universe** — which SNPs get physical rows in
   `.annot.gz` files and in the LD computation (corresponds to LDSC `--bfile` universe).
2. **Regression SNP subset** — which SNPs appear in the output `.l2.ldscore.gz` table
   and therefore enter the LDSC regression (corresponds to LDSC `--print-snps`).

As a result the BED-to-annot step applies the restriction as a **zero-mask** on
annotation values (rows are still written; restricted SNPs get annotation value 0),
while `AnnotationBuilder._filter_aligned_tables_by_global_restriction()` applies it
as a **row filter** at load time. The intermediate `.annot.gz` artifacts and the
in-memory `AnnotationBundle` objects therefore express the restriction in two
incompatible representations.

This plan replaces the single field with two explicit fields:

- `GlobalConfig.ref_panel_snps_path` — row-level filter on baseline template;
  controls which SNPs exist in annotation artifacts and the LD computation.
- `GlobalConfig.regression_snps_path` — single source of truth for the regression SNP
  set; controls **both** which SNPs contribute to the weight LD kernel **and** which
  rows appear in the written `.l2.ldscore.gz` and `.w.l2.ldscore.gz` output files.

It also fixes the BED-to-annot artifact contract so that when `ref_panel_snps_path`
is set, the output `.annot.gz` files physically contain only reference-universe rows
rather than full-template rows with zeroed query values.

---

## Scope and File Change Summary

| File | Change type | Description |
| --- | --- | --- |
| `src/ldsc/config.py` | Rename + add | Replace `restrict_snps_path` with `ref_panel_snps_path` and `regression_snps_path` in `GlobalConfig`; update `validate_config_compatibility()`; rename field in `ReferencePanelBuildConfig` |
| `src/ldsc/_kernel/annotation.py` | Behavior + rename | `_process_baseline_file`: row-filter before BED projection; `AnnotationBuilder._filter_aligned_tables_by_global_restriction`: read `ref_panel_snps_path`; `run_bed_to_annot` argparse: rename CLI flag to `--ref-panel-snps-path` |
| `src/ldsc/_kernel/ref_panel.py` | Rename | `RefPanel._filter_metadata_by_global_restriction`: read `ref_panel_snps_path` |
| `src/ldsc/ref_panel_builder.py` | Rename | `ReferencePanelBuilder.run()`, `run_build_ref_panel()`, `build_parser()`: rename to `ref_panel_snps_path` |
| `src/ldsc/outputs.py` | Remove field + fix filter | Remove `OutputSpec.print_snps_path`; delete `_load_print_snps()`; `LDScoreTableProducer.build()` filters `.l2.ldscore.gz` using `result.regression_snps` (already a `set[str]` from the kernel) |
| `src/ldsc/ldscore_calculator.py` | Remove flags + add flags + wire | Remove `--regression-snps` and `--print-snps`; add `--ref-panel-snps-path` and `--regression-snps-path`; update `_normalize_run_args()` and `run_ldscore()` to propagate both `GlobalConfig` fields; wire `global_config.regression_snps_path` to the kernel |
| `src/ldsc/cli.py` | Rename | Rename `--restrict-snps-path` to `--ref-panel-snps-path` in annotate subcommand |
| `src/ldsc/__init__.py` | No change needed | `GlobalConfig` is already exported; field rename is transparent to callers |
| `docs/config-design.md` | Update | Document two-field model, updated critical/advisory classification, CLI flag table |
| `tests/test_config_identifiers.py` | Update | Rename field references; add `regression_snps_path` normalization test |
| `tests/test_global_config_registry.py` | Update | Rename field references; `ref_panel_snps_path` mismatch → `ConfigMismatchError`; `regression_snps_path` mismatch → warning |
| `tests/test_annotation.py` | Update | Rename field references; add row-filter behavior test |
| `tests/test_ref_panel.py` | Update | Rename field references |
| `tests/test_output.py` | Update | Rename `print_snps_path` → `regression_snps_path`; rename test methods |
| `tests/test_ldscore_workflow.py` | Update | Replace `print_snps=` kwarg with GlobalConfig-based approach |
| `tests/test_package_layout.py` | Update | Remove print_snps test; add ref_panel_snps_path and regression_snps_path tests |

**Files that do NOT need changes:**

- `src/ldsc/sumstats_munger.py` — no restriction logic.
- `src/ldsc/regression_runner.py` — no restriction logic.

---

## Step 0 — Read Before Writing Any Code

Read these files in full before editing anything:

- `src/ldsc/config.py` — complete file; understand `GlobalConfig.__post_init__`,
  `ReferencePanelBuildConfig`, and `validate_config_compatibility()`
- `src/ldsc/_kernel/annotation.py` — focus on:
  - `_process_baseline_file()` (~line 1160)
  - `_build_restrict_mask()` (~line 1110)
  - `_combine_masks()` (~line 1128)
  - `_RestrictResource` dataclass (~line 202)
  - `_build_restrict_resource()` (~line 1064)
  - `_filter_rows_to_restriction()` — does not exist yet; you will write it
  - `AnnotationBuilder._filter_aligned_tables_by_global_restriction()` (~line 678)
  - `AnnotationBuilder.run_bed_to_annot()` / the internal dispatch path (~line 530)
  - The `main_bed_to_annot()` argparse section (~line 770)
- `src/ldsc/_kernel/ref_panel.py` — focus on `RefPanel._filter_metadata_by_global_restriction()` (~line 123)
- `src/ldsc/ref_panel_builder.py` — focus on `ReferencePanelBuilder.run()` (~line 120), `config_from_args()` (~line 440), `run_build_ref_panel()` (~line 464)
- `src/ldsc/ldscore_calculator.py` — focus on the `--print-snps` argument addition (~line 416) and `run_ldscore_from_args()` (~line 640) where `print_snps_path` is resolved into `OutputSpec`
- `src/ldsc/outputs.py` — focus on `OutputSpec` (~line 30), `_load_print_snps()` (~line 489), and `_filter_ldscore_table()` (~line 500)

---

## Step 1 — Update `GlobalConfig` in `src/ldsc/config.py`

### 1a. Replace `restrict_snps_path` with two fields

**Location:** `GlobalConfig` dataclass, approximately lines 78–115.

Remove:

```python
restrict_snps_path: str | PathLike[str] | None = None
```

Add two fields in its place (both optional, both default `None`):

```python
ref_panel_snps_path: str | PathLike[str] | None = None
regression_snps_path: str | PathLike[str] | None = None
```

Update the docstring to describe both fields:

```text
ref_panel_snps_path : str or os.PathLike[str] or None, optional
    Path to a SNP list that defines the annotation and LD-computation universe.
    When set, annotation files are physically limited to these SNPs (rows outside
    the universe are dropped, not zeroed). Applied as a row filter at
    AnnotationBundle load time and at BED-to-annot write time. Corresponds
    conceptually to the SNP universe of LDSC's ``--bfile`` reference panel.
    Default is ``None`` (use the full baseline template universe).
regression_snps_path : str or os.PathLike[str] or None, optional
    Path to a SNP list that restricts which SNPs appear in the written
    ``.l2.ldscore.gz`` output table and therefore enter the LDSC regression.
    Does not affect which rows are in annotation files or which SNPs contribute
    to LD sums. Corresponds to LDSC's ``--print-snps`` flag. Must be a subset of
    the annotation universe (validated at LD-score write time). Default is ``None``
    (write all annotation-universe SNPs to the LD-score table).
```

Update `__post_init__` to normalize both fields:

```python
object.__setattr__(
    self, "ref_panel_snps_path", _normalize_optional_path(self.ref_panel_snps_path)
)
object.__setattr__(
    self, "regression_snps_path", _normalize_optional_path(self.regression_snps_path)
)
```

Remove the old `restrict_snps_path` normalization line.

### 1b. Update `validate_config_compatibility()`

**Location:** `validate_config_compatibility()` function, approximately lines 144–169.

Replace the current `restrict_snps_path` advisory comparison block with:

```python
# ref_panel_snps_path is critical: it defines which rows exist in annotation
# artifacts. Two results computed on different annotation universes cannot be
# safely combined.
if a.ref_panel_snps_path != b.ref_panel_snps_path:
    raise ConfigMismatchError(
        f"ref_panel_snps_path mismatch{prefix}: {a.ref_panel_snps_path!r} vs "
        f"{b.ref_panel_snps_path!r}. These objects were computed on different "
        "annotation universes and cannot be safely merged."
    )
# regression_snps_path is advisory: it controls the output filter, not the
# computation space. Results with different regression SNP sets are combinable.
if a.regression_snps_path != b.regression_snps_path:
    warnings.warn(
        f"regression_snps_path differs{prefix}: {a.regression_snps_path!r} vs "
        f"{b.regression_snps_path!r}. Results are combinable but regression SNP "
        "sets differ.",
        UserWarning,
        stacklevel=2,
    )
```

**Caveat:** `ref_panel_snps_path` is now **critical** (raises `ConfigMismatchError`), not advisory. The existing `config-design.md` table and tests must be updated to reflect this.

### 1c. Rename `restrict_snps_path` in `ReferencePanelBuildConfig`

**Location:** `ReferencePanelBuildConfig` dataclass, approximately lines 328–395.

Rename the field:

```python
# Before:
restrict_snps_path: str | PathLike[str] | None = None

# After:
ref_panel_snps_path: str | PathLike[str] | None = None
```

Update `__post_init__` normalization accordingly:

```python
object.__setattr__(self, "ref_panel_snps_path", _normalize_optional_path(self.ref_panel_snps_path))
```

**Note:** `ReferencePanelBuildConfig.ref_panel_snps_path` is a build-time filter that
restricts which SNPs are retained in the parquet reference panel during construction.
It is distinct from `GlobalConfig.ref_panel_snps_path` but the two should point to
the same file in a consistent workflow. They are separate because `ReferencePanelBuildConfig`
is a build-once config while `GlobalConfig` is a runtime analysis config.

---

## Step 2 — Fix BED-to-Annot Artifact Contract in `src/ldsc/_kernel/annotation.py`

This is the most significant behavioral change in the plan.

### Current behavior (to be replaced)

`_process_baseline_file()` calls `_build_restrict_mask()` which produces a boolean
mask over all baseline rows, then `_combine_masks(overlap_mask, restrict_mask)` zeroes
out annotation values for restricted SNPs while still writing all rows. The output
`.annot.gz` therefore has the full baseline-template row count with 0 annotation values
for out-of-restriction SNPs.

### New behavior

Apply the restriction as a **row filter** on `rows` before the BED projection step.
The output `.annot.gz` contains only reference-universe rows. Annotation values are
derived purely from BED overlap; the restriction is never visible as a spurious zero.

### 2a. Add `_filter_rows_to_restriction()`

Insert a new helper function near the other row-level utilities (after `_build_restrict_mask`,
approximately line 1120):

```python
def _filter_rows_to_restriction(
    rows: list[_BaselineRow],
    restrict_resource: _RestrictResource,
) -> list[_BaselineRow]:
    """Return only the rows whose SNP identity satisfies the restriction set.

    Parameters
    ----------
    rows : list of _BaselineRow
        All rows from one baseline annotation shard.
    restrict_resource : _RestrictResource
        Normalized restriction input (rsid set or BED-intersection resource).

    Returns
    -------
    list of _BaselineRow
        Subset of rows passing the restriction. Order is preserved.
    """
    if restrict_resource.mode == "rsid":
        assert restrict_resource.snp_ids is not None
        return [row for row in rows if row.snp in restrict_resource.snp_ids]
    # chr_pos mode: use the same BED-intersection approach but keep only
    # rows whose (chrom, pos) falls inside the restriction BED.
    assert restrict_resource.bed_path is not None
    pybedtools = _get_pybedtools()
    all_bed = pybedtools.BedTool(
        "\n".join(
            f"{_format_chrom_for_bed(row.chrom)}\t{row.pos - 1}\t{row.pos}"
            for row in rows
        ),
        from_string=True,
    )
    intersected = all_bed.intersect(str(restrict_resource.bed_path), u=True)
    keep_positions: set[tuple[str, int]] = {
        (normalize_chromosome(interval.chrom), int(interval.end))
        for interval in intersected
    }
    return [
        row for row in rows
        if (normalize_chromosome(row.chrom), row.pos) in keep_positions
    ]
```

**Caveat:** The BED-intersection path here constructs a temporary BedTool from a
string. For very large baseline shards (>1M rows) this may be slower than the
previous approach. If performance is a concern, a sorted-merge approach using
`baseline_bed.intersect(restrict_bed, u=True)` on the already-written baseline BED
file is preferable. In that case, re-use the `baseline_bed_path` already constructed
in `_process_baseline_file()` and pass it into this helper.

**Alternative cleaner implementation** (preferred): Since `_process_baseline_file()`
already writes `baseline_bed_path` and constructs `baseline_bed`, pass the already-
materialized BedTool in rather than re-constructing it:

```python
def _filter_rows_to_restriction(
    rows: list[_BaselineRow],
    baseline_bed: Any,             # pybedtools.BedTool already constructed
    restrict_resource: _RestrictResource,
) -> list[_BaselineRow]:
    if restrict_resource.mode == "rsid":
        assert restrict_resource.snp_ids is not None
        return [row for row in rows if row.snp in restrict_resource.snp_ids]
    assert restrict_resource.bed_path is not None
    results = baseline_bed.intersect(str(restrict_resource.bed_path), c=True, wa=True)
    keep_mask = _validate_and_convert_intersection(results, rows, sanity_mode="chr_pos")
    return [row for row, keep in zip(rows, keep_mask) if keep]
```

This reuses the existing `_validate_and_convert_intersection()` helper and avoids
materializing a second BedTool.

### 2b. Update `_process_baseline_file()`

**Location:** `_process_baseline_file()`, approximately lines 1160–1191.

```python
def _process_baseline_file(
    baseline_path: Path,
    bed_paths: Sequence[Path],
    output_dir: Path,
    batch: bool,
    restrict_resource: _RestrictResource | None,
    tempdir: Path,
) -> None:
    """Project all requested BED annotations onto one baseline template file."""
    rows = _read_baseline_annot(baseline_path)
    baseline_bed_path = _write_baseline_bed(rows, tempdir / f"{baseline_path.name}.snps.bed")
    pybedtools = _get_pybedtools()
    baseline_bed = pybedtools.BedTool(str(baseline_bed_path))

    # NEW: apply restriction as a row filter BEFORE BED projection.
    # This replaces the old zero-masking approach. The output file will contain
    # only reference-universe rows; annotation values reflect BED geometry only.
    if restrict_resource is not None:
        rows = _filter_rows_to_restriction(rows, baseline_bed, restrict_resource)
        # Re-write the baseline BED to match the filtered row set.
        baseline_bed_path = _write_baseline_bed(
            rows, tempdir / f"{baseline_path.name}.filtered.snps.bed"
        )
        baseline_bed = pybedtools.BedTool(str(baseline_bed_path))

    if batch:
        annotation_names = [path.stem for path in bed_paths]
        masks = []
        for bed_path in bed_paths:
            overlap_mask = _compute_bed_overlap_mask(rows, baseline_bed, bed_path)
            masks.append(overlap_mask)          # no _combine_masks needed
        output_name = _query_output_name(baseline_path)
        _write_annot_file(output_dir / output_name, rows, annotation_names, masks)
    else:
        for bed_path in bed_paths:
            overlap_mask = _compute_bed_overlap_mask(rows, baseline_bed, bed_path)
            bed_output_dir = output_dir / bed_path.stem
            output_name = _query_output_name(baseline_path)
            _write_annot_file(bed_output_dir / output_name, rows, [bed_path.stem], [overlap_mask])

    pybedtools.cleanup(remove_all=True)
```

Remove `_build_restrict_mask` and `_combine_masks` calls entirely from this function.
`_build_restrict_mask()` and `_combine_masks()` can remain in the file as dead code
initially (for safety during review), then be deleted in a follow-up cleanup.

### 2c. Update `AnnotationBuilder.run_bed_to_annot()` / internal dispatch

**Location:** The block that constructs `restrict_resource` inside
`AnnotationBuilder.run_bed_to_annot()` (or the equivalent internal helper), approximately
lines 574–607.

Change the field read from:

```python
if self.global_config.restrict_snps_path is not None:
    restrict_path = Path(resolve_scalar_path(
        self.global_config.restrict_snps_path, label="global SNP restriction"
    ))
```

to:

```python
if self.global_config.ref_panel_snps_path is not None:
    restrict_path = Path(resolve_scalar_path(
        self.global_config.ref_panel_snps_path, label="reference-panel SNP universe"
    ))
```

### 2d. Update `_filter_aligned_tables_by_global_restriction()`

**Location:** `AnnotationBuilder._filter_aligned_tables_by_global_restriction()`,
approximately lines 678–702.

Change the field read:

```python
# Before:
restrict_path = self.global_config.restrict_snps_path

# After:
restrict_path = self.global_config.ref_panel_snps_path
```

No other changes to this method are needed. Row-filtering behavior is already correct.

### 2e. Update the argparse section in `main_bed_to_annot()`

**Location:** argparse argument definitions near line 779.

```python
# Before:
parser.add_argument(
    "--restrict-snps-path",
    default=None,
    help="Optional SNP restriction file matched by --snp-identifier.",
)

# After:
parser.add_argument(
    "--ref-panel-snps-path",
    default=None,
    help=(
        "Optional path to a SNP list defining the reference-panel universe. "
        "When set, output .annot.gz files contain only rows for these SNPs. "
        "The restriction is a row filter, not a zero-mask."
    ),
)
```

Update the `GlobalConfig` construction that reads from parsed args (~line 796):

```python
# Before:
restrict_snps_path=None if args.restrict_snps_path is None else str(args.restrict_snps_path),

# After:
ref_panel_snps_path=None if args.ref_panel_snps_path is None else str(args.ref_panel_snps_path),
```

**Caveat:** `args.restrict_snps_path` → `args.ref_panel_snps_path` because argparse
converts `--ref-panel-snps-path` to `ref_panel_snps_path` via the hyphen-to-underscore
rule. Verify with `vars(args)` in a test before finalizing.

---

## Step 3 — Update `src/ldsc/_kernel/ref_panel.py`

**Location:** `RefPanel._filter_metadata_by_global_restriction()`, approximately lines 123–135.

```python
def _filter_metadata_by_global_restriction(self, metadata: pd.DataFrame) -> pd.DataFrame:
    """Filter metadata rows to the configured reference-panel SNP universe."""
    restrict_path = self.global_config.ref_panel_snps_path   # renamed field
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

Only the field name changes; the rest of the method body is unchanged.

---

## Step 4 — Update `src/ldsc/ref_panel_builder.py`

### 4a. `ReferencePanelBuilder.run()` field read

**Location:** approximately lines 125–130.

```python
# Before:
if self.global_config.restrict_snps_path:
    restriction_path = resolve_scalar_path(
        self.global_config.restrict_snps_path, ...
    )

# After:
if self.global_config.ref_panel_snps_path:
    restriction_path = resolve_scalar_path(
        self.global_config.ref_panel_snps_path, ...
    )
```

### 4b. `config_from_args()` — reading from parsed args

**Location:** approximately lines 440–453.

```python
# Before:
restrict_snps_path=args.restrict_snps_path,
...
restrict_snps_path=build_config.restrict_snps_path,

# After:
ref_panel_snps_path=args.ref_panel_snps_path,
...
ref_panel_snps_path=build_config.ref_panel_snps_path,
```

### 4c. `run_build_ref_panel()` convenience wrapper

**Location:** approximately lines 464–496.

```python
# Before:
forbidden = sorted({"log_level", "restrict_snps_path"} & set(kwargs))
...
defaults["restrict_snps_path"] = global_config.restrict_snps_path

# After:
forbidden = sorted({"log_level", "ref_panel_snps_path"} & set(kwargs))
...
defaults["ref_panel_snps_path"] = global_config.ref_panel_snps_path
```

Also update the argparse flag in `ref_panel_builder.build_parser()` (wherever
`--restrict-snps-path` is defined for the `build-ref-panel` subcommand):

```python
# Before:
parser.add_argument("--restrict-snps-path", ...)

# After:
parser.add_argument("--ref-panel-snps-path", ...)
```

---

## Step 5 — Remove `--regression-snps` and `--print-snps`; Wire `regression_snps_path` as Single Source of Truth

**Background:** The existing codebase has two separate SNP-list flags in
`ldscore_calculator.py` that both sounded like "regression SNP filtering" but
operated at different stages:

- `--regression-snps` (line 408) → `_load_regression_snps()` → feeds the `regression_snps`
  set into the **weight LD computation kernel**, determining which SNPs populate
  `regression_metadata` / `w_ld`.
- `--print-snps` (line 416) → `OutputSpec.print_snps_path` → **output filter only** —
  restricted which rows appeared in the written `.l2.ldscore.gz` from a separate file.

The design problem: a user could set one without the other, producing paired output
files (`.l2.ldscore.gz` and `.w.l2.ldscore.gz`) with different row sets. The fix is
to collapse both into `GlobalConfig.regression_snps_path` as the single authority:

1. `GlobalConfig.regression_snps_path` is loaded once via `_load_regression_snps()`
   and passed to the weight LD kernel — replacing `--regression-snps`.
2. The kernel stores the retained intersection as `LDScoreResult.regression_snps: set[str]`.
3. The output layer reads `result.regression_snps` directly (already a resolved
   `set[str]`, no second file-read) and filters `.l2.ldscore.gz` rows to that set.
4. `.w.l2.ldscore.gz` is written from `result.regression_metadata / w_ld` which the
   kernel already built with only regression-SNP rows.

Result: `.l2.ldscore.gz` and `.w.l2.ldscore.gz` are guaranteed to have exactly the
same row set, driven by one config field. `OutputSpec` no longer carries any
regression SNP path.

### 5a. Remove `--regression-snps` and `--print-snps`; add `--ref-panel-snps-path` and `--regression-snps-path`

**Location:** `build_parser()` in `src/ldsc/ldscore_calculator.py`, approximately
lines 405–431.

Remove both existing flags:

```python
# Remove both:
parser.add_argument("--regression-snps", ...)
parser.add_argument("--print-snps", ...)
```

Add the two new GlobalConfig-aligned flags in their place:

```python
parser.add_argument(
    "--ref-panel-snps-path",
    default=None,
    help=(
        "Path to a SNP list defining the reference-panel universe. "
        "Annotation rows outside this set are dropped before LD-score computation. "
        "Accepts rsid list or chr_pos table (see --snp-identifier)."
    ),
)
parser.add_argument(
    "--regression-snps-path",
    default=None,
    help=(
        "Path to a SNP list defining the regression SNP set. "
        "Controls both which SNPs contribute to weight LD-score computation "
        "and which rows appear in the written .l2.ldscore.gz and .w.l2.ldscore.gz "
        "files. Accepts rsid list or chr_pos table (see --snp-identifier)."
    ),
)
```

Argparse converts `--ref-panel-snps-path` → `ref_panel_snps_path` and
`--regression-snps-path` → `regression_snps_path` as namespace attribute names.

### 5b. Update `_normalize_run_args()` to include both new fields in `GlobalConfig`

**Location:** `_normalize_run_args()` in `src/ldsc/ldscore_calculator.py`,
approximately lines 618–635.

```python
default_config = GlobalConfig()   # reads module-level _GLOBAL_CONFIG
_ref_panel = getattr(args, "ref_panel_snps_path", None) or default_config.ref_panel_snps_path
_regression = getattr(args, "regression_snps_path", None) or default_config.regression_snps_path

global_config = GlobalConfig(
    snp_identifier=normalized_mode,
    genome_build=default_config.genome_build if getattr(args, "genome_build", None) is None else getattr(args, "genome_build"),
    ref_panel_snps_path=normalize_optional_path_token(_ref_panel),
    regression_snps_path=normalize_optional_path_token(_regression),
    log_level=getattr(args, "log_level", "INFO"),
)
```

This mirrors how `snp_identifier` and `genome_build` are propagated: explicit CLI arg
wins; falls back to the module-level global config. The `or` short-circuit is safe
because `normalize_optional_path_token(None)` returns `None` and paths are never
empty strings after normalization.

### 5c. Wire `global_config.regression_snps_path` to the kernel in `run_ldscore_from_args()`

**Location:** `run_ldscore_from_args()`, approximately line 447.

```python
# Before:
regression_snps = _load_regression_snps(normalized_args.regression_snps, global_config)

# After:
regression_snps = _load_regression_snps(global_config.regression_snps_path, global_config)
```

`_load_regression_snps()` already uses `read_global_snp_restriction()` internally
(same format support as `ref_panel_snps_path`). Only the source of the path changes —
from a now-removed CLI arg to the `GlobalConfig` field. No changes to the function
body itself.

### 5d. Update `run_ldscore()` Python convenience wrapper

**Location:** `run_ldscore()`, approximately lines 592–608.

```python
# Before:
defaults["snp_identifier"] = global_config.snp_identifier
defaults["genome_build"] = global_config.genome_build
defaults["log_level"] = global_config.log_level
forbidden = sorted({"snp_identifier", "genome_build", "log_level"} & set(kwargs))

# After:
defaults["snp_identifier"] = global_config.snp_identifier
defaults["genome_build"] = global_config.genome_build
defaults["ref_panel_snps_path"] = global_config.ref_panel_snps_path
defaults["regression_snps_path"] = global_config.regression_snps_path
defaults["log_level"] = global_config.log_level
forbidden = sorted(
    {"snp_identifier", "genome_build", "ref_panel_snps_path", "regression_snps_path", "log_level"}
    & set(kwargs)
)
```

### 5e. Update `_output_spec_from_args()` — remove `print_snps_path`

**Location:** `_output_spec_from_args()` in `src/ldsc/ldscore_calculator.py`,
approximately lines 649–661.

```python
# Before:
return OutputSpec(
    ...
    print_snps_path=normalize_optional_path_token(getattr(args, "print_snps", None)),
    ...
)

# After:
return OutputSpec(
    ...
    # print_snps_path removed — row filtering is driven by result.regression_snps
    ...
)
```

The function signature stays `_output_spec_from_args(args: argparse.Namespace)` —
no `global_config` parameter needed since `OutputSpec` no longer carries any
regression SNP path.

### 5f. Update `src/ldsc/outputs.py` — remove `print_snps_path`; filter via `result.regression_snps`

**Location:** `OutputSpec` dataclass, approximately line 39. Remove the field and its
`__post_init__` normalization line:

```python
# Remove:
print_snps_path: str | PathLike[str] | None = None
```

**Location:** `LDScoreTableProducer.build()`, approximately line 199. Replace the
`_load_print_snps()` call with a direct read of `result.regression_snps`:

```python
# Before:
print_snps = _load_print_snps(output_spec.print_snps_path)
table = _filter_ldscore_table(table, print_snps)

# After:
_cfg = getattr(result, "config_snapshot", None) or get_global_config()
regression_filter = set(result.regression_snps) if _cfg.regression_snps_path else None
table = _filter_ldscore_table(table, regression_filter)
```

Apply the same `regression_filter` variable to the per-chromosome branch (line 210)
and update the guard and error message at the end:

```python
# Before:
if print_snps is not None and not artifacts:
    raise ValueError("After merging with --print-snps, no SNPs remain.")

# After:
if regression_filter is not None and not artifacts:
    raise ValueError("After applying regression_snps_path filter, no SNPs remain.")
```

**Delete `_load_print_snps()`** (approximately lines 489–497) entirely.

**`WeightLDProducer.build()` — no changes needed.** The `regression_metadata / w_ld`
tables are produced by the kernel using the regression SNP set. They already have
only regression-SNP rows. The output layer writes them as-is.

**Note on `_filter_ldscore_table()`:** Keep name and signature. Rename its internal
`print_snps` parameter to `regression_snps` for clarity; behavior is unchanged.

**Import:** Add `get_global_config` from `src/ldsc/config.py` to `outputs.py` if
not already imported (needed as fallback when `result.config_snapshot` is `None`).
Remove any import of `read_global_snp_restriction` added for the old approach (it is
no longer needed in `outputs.py`).

### 5g. Update tests

**`tests/test_package_layout.py`**

- Remove `test_ldscore_subcommand_accepts_print_snps` — `--print-snps` no longer exists.
- Remove `test_ldscore_subcommand_accepts_regression_snps` — `--regression-snps` no
  longer exists.
- Add two new tests:

  ```python
  def test_ldscore_subcommand_accepts_ref_panel_snps_path(self): ...
  def test_ldscore_subcommand_accepts_regression_snps_path(self): ...
  ```

**`tests/test_output.py`**

Remove `print_snps_path=` from all `OutputSpec` constructor calls (the field no
longer exists). The tests must be rewritten to drive filtering via `result.regression_snps`
and a `config_snapshot` with `regression_snps_path` set, rather than via `OutputSpec`:

- `test_print_snps_filters_only_reference_ldscore_artifact` →
  `test_regression_snps_path_filters_only_reference_ldscore_artifact`
- `test_print_snps_filters_per_chrom_reference_outputs_only` →
  `test_regression_snps_path_filters_per_chrom_reference_outputs_only`
- `test_print_snps_raises_when_no_reference_rows_remain` →
  `test_regression_snps_path_raises_when_no_reference_rows_remain`; update expected
  error message from `"After merging with --print-snps, no SNPs remain."` to
  `"After applying regression_snps_path filter, no SNPs remain."`.

For each test, construct the mock `LDScoreResult` with:

- `regression_snps=<the filtered SNP set>`
- `config_snapshot=GlobalConfig(regression_snps_path="path/to/snps.txt")`

**`tests/test_ldscore_workflow.py`**

Update `test_run_ldscore_from_args_print_snps_filters_only_written_reference_table`:

- Rename the test to `test_regression_snps_path_produces_aligned_paired_outputs`.
- Replace `print_snps=str(...)` kwarg with `regression_snps_path=str(...)` in the
  `run_ldscore()` call (or set via `set_global_config()`).
- Assert that both `.l2.ldscore.gz` and `.w.l2.ldscore.gz` have the same row set
  (not just that `.l2.ldscore.gz` is filtered).

---

## Step 6 — Update `src/ldsc/cli.py`

### 6a. Rename the annotate subcommand flag

**Location:** `_add_annotate_arguments()`, approximately line 151.

```python
# Before:
parser.add_argument(
    "--restrict-snps-path",
    default=None,
    help="Optional global SNP restriction file.",
)

# After:
parser.add_argument(
    "--ref-panel-snps-path",
    default=None,
    help=(
        "Path to a SNP list defining the reference-panel universe for annotation. "
        "Rows outside this set are dropped from output .annot.gz files."
    ),
)
```

### 6b. Update `_run_annotate()` exclusion set

**Location:** `_run_annotate()`, approximately line 137.

```python
# Before:
_namespace_to_argv(args, exclude={..., "restrict_snps_path", ...})

# After:
_namespace_to_argv(args, exclude={..., "ref_panel_snps_path", ...})
```

---

## Step 7 — Update Tests

### 7a. `tests/test_config_identifiers.py`

- Replace all `restrict_snps_path=...` with `ref_panel_snps_path=...` in `GlobalConfig`
  construction calls (approximately lines 56 and 210–246).
- Add a test verifying `regression_snps_path` is normalized the same way:

  ```python
  config = GlobalConfig(regression_snps_path=Path("output") / "hm3.txt")
  self.assertEqual(config.regression_snps_path, "output/hm3.txt")
  ```

### 7b. `tests/test_global_config_registry.py`

- Replace `restrict_snps_path="restrict.txt"` with `ref_panel_snps_path="restrict.txt"`
  in construction calls (approximately lines 89, 110, 245–246).
- Update test `test_validate_config_compatibility_warns_on_restrict_snps_difference`
  (approximately line 244):
  - Rename to `test_validate_config_compatibility_raises_on_ref_panel_snps_difference`
  - Assert `ConfigMismatchError` is raised (not just a warning):

    ```python
    left = GlobalConfig(genome_build="hg38", snp_identifier="chr_pos", ref_panel_snps_path="left.txt")
    right = GlobalConfig(genome_build="hg38", snp_identifier="chr_pos", ref_panel_snps_path="right.txt")
    with self.assertRaises(ConfigMismatchError, msg="ref_panel_snps_path mismatch"):
        validate_config_compatibility(left, right)
    ```

- Add a new test for `regression_snps_path` difference emitting a warning (not error):

  ```python
  left = GlobalConfig(genome_build="hg38", regression_snps_path="a.txt")
  right = GlobalConfig(genome_build="hg38", regression_snps_path="b.txt")
  with self.assertWarns(UserWarning):
      validate_config_compatibility(left, right)
  ```

### 7c. `tests/test_annotation.py`

- Replace `restrict_snps_path=...` with `ref_panel_snps_path=...` (approximately
  lines 85 and 108).
- Add a test asserting that BED-to-annot output with `ref_panel_snps_path` set
  produces only reference-universe rows:

  ```python
  # Confirm that the output .annot.gz has N_restricted rows, not N_baseline rows.
  # Confirm that no row has annotation value 0 due solely to restriction
  # (all 0s should come from BED non-overlap only).
  ```

### 7d. `tests/test_ref_panel.py`

- Replace `restrict_snps_path=...` with `ref_panel_snps_path=...` (approximately
  line 74).

### 7e. New test for two-field behavior

Add a test (in `test_config_identifiers.py` or a new `test_snp_set_api.py`) verifying
that setting both fields on the same `GlobalConfig` is valid and that they are treated
independently:

```python
cfg = GlobalConfig(
    ref_panel_snps_path="universe/hm3_extended.txt",
    regression_snps_path="filters/hm3.txt",
)
self.assertEqual(cfg.ref_panel_snps_path, "universe/hm3_extended.txt")
self.assertEqual(cfg.regression_snps_path, "filters/hm3.txt")
```

---

## Step 8 — Update `docs/config-design.md`

Update the **Critical vs. Advisory Fields** table to reflect the new classification:

| Field | Class | Mismatch behavior |
| --- | --- | --- |
| `genome_build` | Critical | Hard error (`ConfigMismatchError`) |
| `snp_identifier` | Critical | Hard error (`ConfigMismatchError`) |
| `ref_panel_snps_path` | Critical | Hard error (`ConfigMismatchError`) |
| `regression_snps_path` | Advisory | Warning only |
| `log_level` | Advisory | Ignored during compatibility checks |
| `fail_on_missing_metadata` | Advisory | Ignored during compatibility checks |

Add a new section **Two-Level SNP-Set Model** (see the companion design doc
`docs/config-design.md` for the full text).

---

## Known Caveats and Edge Cases

### C1. Breaking API change for users of `restrict_snps_path`

Renaming `GlobalConfig.restrict_snps_path` is a **breaking public API change**. Any
user code that passes `restrict_snps_path=...` to `GlobalConfig()` will raise
`TypeError` at construction time. Options:

1. **Hard break (recommended):** Remove the old field entirely. Users are on a
   pre-release internal package; the API contract is not yet frozen.
2. **Deprecation shim:** Add a `restrict_snps_path` property that emits
   `DeprecationWarning` and redirects to `ref_panel_snps_path`. Remove after one
   release cycle.

If the hard-break option is chosen, grep for all call sites outside `src/` and
`tests/` (notebooks, tutorial files) and update them as part of this PR.

### C2. Intermediate artifact row-count change

When `ref_panel_snps_path` is set, `query.<chrom>.annot.gz` files will have fewer
rows than the baseline template. Any downstream workflow that assumes the output
`.annot.gz` and the baseline `.annot.gz` have identical row counts will break.

This is **intentional and correct**, but it means:

- The AnnotationBundle loading path must accept baseline and query files with
  different row counts when one or both are already restricted. Verify that
  `AnnotationBuilder.run()` does not hard-fail when a query file has fewer rows
  than the raw baseline template.
- If users supply these restricted query files to the ldscore step alongside an
  unrestricted baseline, alignment will fail. Document that the same `ref_panel_snps_path`
  must be used for both the annotate and ldscore steps.

### C3. `ReferencePanelBuildConfig.ref_panel_snps_path` vs `GlobalConfig.ref_panel_snps_path`

These two fields serve the same concept but at different stages:

- `ReferencePanelBuildConfig.ref_panel_snps_path` restricts which SNPs are retained
  in the parquet files at panel-build time.
- `GlobalConfig.ref_panel_snps_path` restricts which SNPs are loaded at annotation
  and LD-computation time.

In a consistent workflow, both should point to the same file. There is currently no
validation that enforces this. A future enhancement could warn when the two paths
differ, but it requires reading parquet metadata at runtime and is out of scope here.

### C4. `regression_snps_path` is honored automatically through the kernel result

`GlobalConfig.regression_snps_path` is loaded in `run_ldscore_from_args()` and
passed to the weight LD kernel. The kernel stores the retained SNP intersection in
`LDScoreResult.regression_snps`. The output layer then reads `result.regression_snps`
directly — no `OutputSpec` field is involved. This means:

- Users who call `run_ldscore()` or `run_ldscore_from_args()` get correct filtering
  automatically.
- Users who construct an `LDScoreResult` directly (e.g. in tests) must populate
  `regression_snps` and `config_snapshot.regression_snps_path` explicitly to trigger
  output filtering.
- `OutputSpec` no longer has any regression SNP path field; do not add one.

### C5. Both `--print-snps` and `--regression-snps` are fully removed — no backward-compat shim

Both flags are removed from the `ldscore` subcommand. The replacement is the single
`--regression-snps-path` flag. Users who previously passed either flag on the command
line must migrate:

- `--print-snps path/to/snps.txt` → `--regression-snps-path path/to/snps.txt`
- `--regression-snps path/to/snps.txt` → `--regression-snps-path path/to/snps.txt`

From Python, use `set_global_config(GlobalConfig(regression_snps_path=...))` before
calling `run_ldscore()`.

Grep for `print-snps`, `print_snps`, `regression-snps`, and `regression_snps` (the
old standalone form) in notebooks and tutorial files outside `src/` and `tests/` and
update them as part of this PR.

### C6. Zero-mask removal may break existing user-inspected artifacts

Users who currently depend on the behavior where out-of-restriction SNPs appear
with annotation value `0` in `query.*.annot.gz` will see a different artifact after
this change. Those SNPs will simply not appear as rows. If users were relying on the
zero values to identify "which SNPs were excluded," they must now inspect the row
set difference between the baseline template and the output query file. Log a message
at `INFO` level during BED-to-annot projection when `ref_panel_snps_path` is set,
reporting the number of rows retained and dropped.

### C7. `_build_restrict_mask()` and `_combine_masks()` become dead code

After Step 2, these functions are no longer called in `_process_baseline_file()`.
They can be deleted, but do a final grep before removing to confirm no other caller
uses them.

---

## Execution Order

The steps are mostly independent but should be executed in this order to keep tests
green at each checkpoint:

1. Step 1 (`config.py`) — makes the new field names available; all downstream tests will break until Step 7.
2. Step 2 (`_kernel/annotation.py`) — most invasive; write `_filter_rows_to_restriction()` first, then update `_process_baseline_file()`.
3. Step 3 (`_kernel/ref_panel.py`) — simple field rename.
4. Step 4 (`ref_panel_builder.py`) — rename cascade including `build_parser()` CLI flag.
5. Step 5 (`outputs.py` + `ldscore_calculator.py`) — remove `--regression-snps` and `--print-snps`; add `--ref-panel-snps-path` and `--regression-snps-path` (5a); update `_normalize_run_args()` (5b); wire kernel feed (5c); update `run_ldscore()` (5d); clean `_output_spec_from_args()` (5e); remove `OutputSpec.print_snps_path` and filter via `result.regression_snps` (5f). Do `outputs.py` first since `ldscore_calculator.py` references `OutputSpec`.
6. Step 6 (`cli.py`) — rename annotate subcommand flag.
7. Step 7 (tests) — update all renamed field references; delete removed tests; add new tests.
8. Step 8 (docs) — `config-design.md` updated; grep for `restrict_snps_path`, `print_snps_path`, `print-snps`, `regression-snps` (old standalone flag) in all docs and tutorials outside `src/` and `tests/`.

Run `python -m unittest discover -s tests -p 'test*.py' -v` after Steps 1–2 and
again after Step 7 to confirm no regressions.
