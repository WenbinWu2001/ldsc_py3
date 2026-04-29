# Genome-Build Configuration Refactor: Design Document

## Status

Implemented.

This document is a historical design record. Sections that say "current"
describe the codebase at the time the refactor was proposed; the implemented
state is summarized in §4 and §6.

---

## 1. Problem Statement

`genome_build` is threaded through the codebase as a single field that conflates
four distinct concepts:

| Concept | Current home | Problem |
|---|---|---|
| **Source build** — what build the raw input files are in | `ReferencePanelBuildConfig.source_genome_build` | No problem; explicit and required. |
| **Runtime build** — what build the user declares for this run | `GlobalConfig.genome_build` | Defaults to `"hg38"`, so omitting `--genome-build` silently applies hg38 assumptions to all chr_pos reads. |
| **Inferred build** — what build was detected by scoring CHR/POS values | Resolved locally inside `ref_panel.py`, `ldscore.py`, `identifiers.py`, `annotation.py` | Inference happens lazily and independently at four sites; the resolved value is never written back, so each site may re-infer or conflict. |
| **File-declared build** — what build is stored in parquet row-group metadata | Validated inside `_init_canonical_path()` | No problem; used as a cross-check, not a config source. |

The three concrete defects this causes:

1. **Silent hg19-as-hg38 substitution.** Any user who passes hg19 chr_pos data
   without `--genome-build hg19` gets `"hg38"` baked into metadata paths and
   parquet column selection — no warning, no error.

2. **rsid workflows inherit a meaningless default.** `GlobalConfig.genome_build`
   defaults to `"hg38"` even when `snp_identifier="rsid"`, where it is never
   consulted. The value creates a false impression that the rsid run made a
   genome-build commitment.

3. **`RefPanelConfig.genome_build` is a silent redundant copy.** It is always
   set from `GlobalConfig.genome_build` at `ldscore_calculator.py:848,856`, and
   metadata loading reads `GlobalConfig` directly anyway (`ref_panel.py:341`).
   The field is present but not authoritative.

---

## 2. Goals

1. A chr_pos workflow that omits `--genome-build` must fail with an actionable
   error — not silently use hg38.
2. `"auto"` inference runs once at the workflow entry point, logs the result,
   and produces a concrete build that flows through the rest of the run.
3. Internal components (kernel modules, metadata loaders, parquet readers) never
   see `"auto"` — they assert or type-check for a concrete build.
4. rsid workflows have `genome_build=None` throughout and never touch
   genome-build logic.
5. `RefPanelConfig.genome_build` is removed; reference-panel code reads
   `GlobalConfig.genome_build` directly.

### Non-goals

- Changing the scoring algorithm in `genome_build_inference.py`.
- Supporting per-panel genome-build mixing within a single run.
- Adding a `--genome-build` override that ignores parquet file-declared metadata.

---

## 3. Proposed Design

### 3.1 `GlobalConfig.genome_build` default: `None`

**File:** `src/ldsc/config.py`

Change the `genome_build` field default from `"hg38"` to `None`:

```python
# Before
genome_build: GenomeBuildInput | None = "hg38"

# After
genome_build: GenomeBuildInput | None = None
```

Replace the existing `normalize_genome_build` call block in `__post_init__` with
the following (keep the `normalize_genome_build` call, then add the new checks
immediately after it):

```python
def __post_init__(self) -> None:
    if self.snp_identifier not in {"rsid", "chr_pos"}:
        raise ValueError("snp_identifier must be 'rsid' or 'chr_pos'.")
    object.__setattr__(self, "genome_build", normalize_genome_build(self.genome_build))
    object.__setattr__(self, "log_level", _normalize_log_level(self.log_level))
    # New validation block:
    if self.snp_identifier == "chr_pos" and self.genome_build is None:
        raise ValueError(
            "genome_build is required when snp_identifier='chr_pos'. "
            "Pass 'auto' to infer from data, or 'hg19'/'hg38' explicitly."
        )
    if self.snp_identifier == "rsid" and self.genome_build == "auto":
        raise ValueError(
            "genome_build='auto' is not valid for snp_identifier='rsid'."
        )
    if self.snp_identifier == "rsid" and self.genome_build is not None:
        import warnings
        warnings.warn(
            "genome_build is set but will be ignored in rsid mode.",
            UserWarning,
            stacklevel=2,
        )
        object.__setattr__(self, "genome_build", None)
```

**Module singleton and reset** — update both lines that construct `GlobalConfig()`
with no arguments; the new default raises for chr_pos:

```python
# config.py line 115 — was: GlobalConfig()
_GLOBAL_CONFIG: GlobalConfig = GlobalConfig(snp_identifier="rsid")

# config.py line 138 inside reset_global_config() — was: GlobalConfig()
_GLOBAL_CONFIG = GlobalConfig(snp_identifier="rsid")
```

**Impact on `validate_config_compatibility`** — no change required. The existing
check `a.genome_build != b.genome_build` handles all cases correctly:

- Both `None` (rsid + rsid) → passes.
- Both concrete and matching → passes.
- Both concrete and mismatching → raises `ConfigMismatchError`.
- `None` vs concrete — only possible if `snp_identifier` also differs, which
  the `snp_identifier` check already catches first.

---

### 3.2 Eager `"auto"` resolution: `resolve_genome_build()`

**File:** `src/ldsc/genome_build_inference.py`

Add one new public function. This is the **only** place in the codebase that
calls `infer_chr_pos_build()` for build inference after this refactor. All
workflow entry points call it once before constructing `GlobalConfig`.

#### Function specification

```python
def resolve_genome_build(
    hint: str | None,
    snp_identifier: str,
    sample_frame: pd.DataFrame | None,
    *,
    context: str,
    logger=None,
) -> str | None:
    """Resolve a genome-build hint to a concrete build or None.

    Parameters
    ----------
    hint : {"auto", "hg19", "hg38", None}
        Raw value from CLI or config. May be "auto", a concrete build string,
        or None. Alias inputs (e.g. "hg37", "GRCh38") are accepted and
        normalized via normalize_genome_build().
    snp_identifier : {"rsid", "chr_pos"}
        Already-normalized SNP identifier mode.
    sample_frame : pd.DataFrame or None
        DataFrame with columns "CHR" and "POS" representing a small sample of
        the user's input data. Required when hint="auto" and
        snp_identifier="chr_pos"; ignored (and may be None) otherwise.
    context : str
        Human-readable label for the data source, used in log and error
        messages. Example: "annotation chr1 of baseline.@.annot.gz".
    logger : logging.Logger-like, optional
        Receives inference summary messages.

    Returns
    -------
    str or None
        "hg19", "hg38", or None. Never returns "auto".
        Returns None when snp_identifier="rsid" regardless of hint.

    Raises
    ------
    ValueError
        - hint="auto" and sample_frame is None.
        - hint="auto" and sample_frame has insufficient HM3 overlap
          (fewer than MIN_INSPECTED_REFERENCE_SNPS=200 uniquely informative
          SNPs). Error message directs user to pass --genome-build explicitly.
        - hint="auto" and snp_identifier="rsid" (validated at GlobalConfig
          construction; this function re-checks defensively).
    """
```

#### Behaviour table

| `snp_identifier` | `hint` | `sample_frame` | Return value | Side effect |
|---|---|---|---|---|
| `"rsid"` | any | any | `None` | None |
| `"chr_pos"` | `"hg19"` or `"hg38"` | any | hint (normalized) | None |
| `"chr_pos"` | `"auto"` | valid DataFrame | inferred build | logs summary |
| `"chr_pos"` | `"auto"` | `None` | — | raises `ValueError` |
| `"chr_pos"` | `"auto"` | insufficient overlap | — | raises `ValueError` |

#### Log output

When inference runs (`hint="auto"`, chr_pos), the function logs at INFO level:

```
Inferred genome build for <context>: <ChrPosBuildInference.summary_message>
```

If `coordinate_basis="0-based"` was detected, log at WARNING level instead
(matching the existing `resolve_chr_pos_table` convention).

#### Exports

Add `resolve_genome_build` to:
- `__all__` in `genome_build_inference.py`
- the public re-export block in `src/ldsc/__init__.py`

---

### 3.3 Per-workflow sampling helper: `sample_frame_from_chr_pattern()`

**File:** `src/ldsc/_chr_sampler.py` (new private module, not kernel)

This helper is used by `ldscore_calculator.py` and `annotation_builder.py`. It
lives outside the kernel because it performs I/O and workflow-level path
resolution. It does not belong in `genome_build_inference.py`, which must remain
pure (Q4 decision).

#### Function specification

```python
def sample_frame_from_chr_pattern(
    tokens: Sequence[str],
    chromosomes: Sequence[str] | None,
    *,
    logger=None,
    nrows: int = 5000,
) -> tuple[pd.DataFrame, str]:
    """Load a small CHR/POS sample from the first resolvable chromosome file.

    Parameters
    ----------
    tokens : sequence of str
        Path tokens (already split by split_cli_path_tokens). Each token may
        contain "@" as a chromosome placeholder, or be a literal path.
    chromosomes : sequence of str or None
        Chromosome list from args.chromosomes. Used as fallback when chr1 is
        not found. If None or empty, only "1" is tried.
    logger : logging.Logger-like, optional
        Receives the pre-read INFO message.
    nrows : int
        Maximum number of rows to read. Default 5000.

    Returns
    -------
    frame : pd.DataFrame
        Columns ["CHR", "POS"] (normalized). Exactly nrows rows or fewer.
    resolved_path : str
        The file path that was actually read (used in context strings for
        resolve_genome_build()).

    Raises
    ------
    ValueError
        If no token with "@" exists in tokens, or no chromosome substitution
        produces an existing file.
    """
```

#### Algorithm

1. Find the first token that contains `"@"`. Tokens without `"@"` are literal
   paths and are not suitable for chromosome substitution — skip them.
   If no `"@"` token exists, raise `ValueError`:
   ```
   No per-chromosome pattern (@) found in tokens; cannot sample for
   genome-build inference. Pass --genome-build hg19 or --genome-build hg38.
   ```
2. Build the candidate chromosome list: `["1"] + list(chromosomes or [])`.
   Deduplicate while preserving order.
3. For each candidate chromosome `c`, call
   `substitute_chromosome(token, c)` (from `path_resolution.py`).
   If the resulting path exists (`os.path.exists`), use it. Stop at the first
   match.
4. If no candidate produces an existing file, raise `ValueError`:
   ```
   Could not find a readable chromosome file for pattern <token>.
   Tried chromosomes: 1, <chromosomes>. Pass --genome-build explicitly.
   ```
5. Log at INFO level **before** reading:
   ```
   Auto-inferring genome build from first <nrows> rows of <resolved_path>
   (chr<c> of <token>). Pass --genome-build hg19 or --genome-build hg38 to
   skip inference.
   ```
6. Read the file:
   - `.annot.gz` / `.annot` files: use `pd.read_csv(..., sep="\t",
     compression="infer", nrows=nrows)`.
   - `.parquet` files: use `pyarrow.parquet.read_table(...,
     columns=None).to_pandas().head(nrows)` — then select CHR and position
     columns from the result.
7. Identify the position column from the file's header using
   `infer_chr_pos_columns(df.columns, context=resolved_path)` from
   `column_inference.py`. This returns `(chr_col, pos_col)`.
8. Return `df.rename(columns={chr_col: "CHR", pos_col: "POS"})[["CHR", "POS"]]`
   and `resolved_path`.

#### Imports needed in `_chr_sampler.py`

```python
import os
from collections.abc import Sequence
import pandas as pd
from .column_inference import infer_chr_pos_columns
from .path_resolution import substitute_chromosome
```

---

### 3.4 Workflow entry-point integration

Each workflow's entry point runs the following sequence **before** constructing
`GlobalConfig`, when `snp_identifier == "chr_pos"`:

#### ldscore (`ldscore_calculator.py:_normalize_run_args`, ~line 788)

`ldscore` is the only workflow that has both annotation and reference panel
inputs. When `hint == "auto"`, **both** sources are sampled and their inferred
builds are cross-checked.

```python
# After deriving normalized_mode from args.snp_identifier:

if normalized_mode == "chr_pos":
    hint = getattr(args, "genome_build", None)   # raw CLI value or None

    if hint == "auto":
        # --- Annotation source ---
        annot_tokens = split_cli_path_tokens(getattr(args, "baseline_annot_sources", None))
        annot_sample, annot_path = sample_frame_from_chr_pattern(
            annot_tokens, getattr(args, "chromosomes", None), logger=logger
        )
        annot_build = resolve_genome_build(
            hint, normalized_mode, annot_sample,
            context=f"annotation chr1 of {annot_tokens[0]}",
            logger=logger,
        )

        # --- Reference panel source ---
        ref_tokens = split_cli_path_tokens(getattr(args, "r2_sources", None))
        ref_sample, ref_path = sample_frame_from_chr_pattern(
            ref_tokens, getattr(args, "chromosomes", None), logger=logger
        )
        ref_build = resolve_genome_build(
            hint, normalized_mode, ref_sample,
            context=f"reference panel chr1 of {ref_tokens[0]}",
            logger=logger,
        )

        # --- Cross-check ---
        if annot_build != ref_build:
            raise ValueError(
                f"Genome build mismatch: annotation files inferred as "
                f"{annot_build!r} but reference panel files inferred as "
                f"{ref_build!r}. Ensure both inputs use the same build, "
                "or pass --genome-build hg19 or --genome-build hg38 explicitly."
            )
        resolved_build = annot_build

    elif hint is not None:
        # Concrete hint: trust the user; parquet cross-check runs at read time.
        resolved_build = normalize_genome_build(hint)
    else:
        raise ValueError(
            "genome_build is required when snp_identifier='chr_pos'. "
            "Pass --genome-build auto, --genome-build hg19, or --genome-build hg38."
        )

    global_config = GlobalConfig(
        snp_identifier=normalized_mode,
        genome_build=resolved_build,
        log_level=getattr(args, "log_level", "INFO"),
    )
else:
    # rsid mode: genome_build not needed; None enforced by GlobalConfig.__post_init__
    global_config = GlobalConfig(
        snp_identifier=normalized_mode,
        log_level=getattr(args, "log_level", "INFO"),
    )
```

Remove the old `default_config = GlobalConfig()` fallback block at lines 814–817.

**Edge case:** when only one source is present (annotation-only or ref-panel-only
run), skip the cross-check and call `resolve_genome_build()` once against the
available source. Detect which source is missing by checking whether
`split_cli_path_tokens(...)` returns an empty list.

#### annotate (`annotation_builder.py`)

Same dual-sampling + cross-check pattern as ldscore. Annotation source:
`--baseline-annot` tokens. Reference panel source: `--r2-table` tokens (if any;
otherwise annotation-only, no cross-check).

#### munge-sumstats (`sumstats_munger.py`)

No reference panel input; no cross-check.

```python
if normalized_mode == "chr_pos" and hint == "auto":
    sample_frame = pd.read_csv(
        sumstats_file, sep="\t", compression="infer", nrows=5000,
        usecols=lambda col: col.upper() in {"CHR", "BP", "POS", "POSITION"},
    )
    chr_col, pos_col = infer_chr_pos_columns(sample_frame.columns, context=sumstats_file)
    sample_frame = sample_frame.rename(columns={chr_col: "CHR", pos_col: "POS"})[["CHR", "POS"]]
    logger.info(
        "Auto-inferring genome build from first 5000 rows of %s. "
        "Pass --genome-build hg19 or --genome-build hg38 to skip.",
        sumstats_file,
    )
    resolved_build = resolve_genome_build(
        hint, normalized_mode, sample_frame, context=sumstats_file, logger=logger
    )
```

#### build-ref-panel (`ref_panel_builder.py`)

`source_genome_build` is always a concrete build (required field, no `"auto"`
accepted). Call `resolve_genome_build()` with `sample_frame=None` — the function
returns the concrete hint unchanged without touching the sample frame.

```python
resolved_build = resolve_genome_build(
    source_genome_build, "chr_pos", None, context="build-ref-panel source", logger=logger
)
```

---

### 3.5 Remove `RefPanelConfig.genome_build`

**Files:** `src/ldsc/config.py`, `src/ldsc/ldscore_calculator.py`,
`src/ldsc/_kernel/ref_panel.py`

- Remove the field `genome_build: GenomeBuildInput | None = None` from
  `RefPanelConfig`.
- Remove `object.__setattr__(self, "genome_build", normalize_genome_build(self.genome_build))`
  from `RefPanelConfig.__post_init__`.
- Remove `genome_build=global_config.genome_build` from both `RefPanelConfig(...)`
  constructions in `_ref_panel_from_args()` (`ldscore_calculator.py:848, 856`).
- Replace `self.spec.genome_build or self.global_config.genome_build` at
  `ref_panel.py:269` with `self.global_config.genome_build`.

---

### 3.6 Remove scattered lazy `"auto"` blocks

After §3.4 is in place, the following blocks become dead code and are deleted:

| File | Lines | What to delete |
|---|---|---|
| `_kernel/ref_panel.py` | 341–347 | `if snp_identifier == "chr_pos" and global_config.genome_build == "auto":` block + `resolve_chr_pos_table` call. Replace with `assert global_config.genome_build in {"hg19", "hg38", None}`. |
| `_kernel/ldscore.py` | 1257–1265 | `SortedR2BlockReader.__init__()` `"auto"` branch. Replace with `assert genome_build in {"hg19", "hg38", None}`. |
| `_kernel/identifiers.py` | 230 | `infer_build = (genome_build == "auto")` and the `if infer_build:` branch in `_finalize_chr_pos_restriction_frame`. |
| `_kernel/annotation.py` | 497 | auto-inference trigger. |

Remove `validate_auto_genome_build_mode()` calls from all kernel-internal sites.
Keep the call in `RefPanel.__init__` as a temporary defensive assert; remove it
in a follow-up once the Step 4 assertions have been validated by a full test run.

---

### 3.7 CLI changes

**Files:** `src/ldsc/ldscore_calculator.py` (parser), `src/ldsc/cli.py`

`--genome-build` in the ldscore and annotate parsers currently defaults to
`None`. After the change, `None` passed at runtime is the sentinel that triggers
the `ValueError` in §3.4. The parser `default=None` is correct; no change to the
parser default is needed.

Update the help string:

```
--genome-build {auto,hg19,hg37,hg38,GRCh37,GRCh38}
    Genome build for chr_pos inputs. Required when --snp-identifier chr_pos
    (the default). Use 'auto' to infer hg19/hg38 and 0-based/1-based
    coordinates from data. Not used when --snp-identifier rsid.
```

`build-ref-panel` already uses `--source-genome-build` (required, no `"auto"`).
No change needed there.

---

## 4. Invariants After the Refactor

| State | `snp_identifier` | `genome_build` |
|---|---|---|
| After `GlobalConfig` construction, chr_pos | `"chr_pos"` | `"hg19"`, `"hg38"`, or `"auto"` |
| After `resolve_genome_build()` at workflow entry | `"chr_pos"` | `"hg19"` or `"hg38"` |
| After `GlobalConfig` construction, rsid | `"rsid"` | `None` |
| After `resolve_genome_build()`, rsid | `"rsid"` | `None` |
| Inside any kernel module | any | never `"auto"` |

---

## 5. What Does Not Change

- `ReferencePanelBuildConfig.source_genome_build` — required and explicit.
- The inference algorithm in `infer_chr_pos_build()` and `resolve_chr_pos_table()`
  — unchanged.
- `validate_auto_genome_build_mode()` — kept as a utility function; its
  scattered kernel-internal call sites are removed (the `GlobalConfig.__post_init__`
  checks now enforce the same constraint at construction time), but the function
  itself is not deleted.
- Parquet file-declared build cross-check in `_init_canonical_path()` — unchanged;
  it remains a hard-error validation step comparing the file's embedded build
  against `GlobalConfig.genome_build`.
- Legacy format compatibility boundaries.
- The packaged HM3 reference file at `src/ldsc/data/hm3_chr_pos_reference.tsv.gz`
  (11 000 SNPs, 500 per chromosome) — already sufficient and unchanged.

---

## 6. Design Decisions

**Q1 — `GlobalConfig()` no-args → Option C (strictest).**
`GlobalConfig()` raises because `snp_identifier` defaults to `"chr_pos"` and
`genome_build` defaults to `None`, which fails the new validation. The module
singleton and `reset_global_config()` become `GlobalConfig(snp_identifier="rsid")`.
All CLI workflows construct `GlobalConfig` with an explicit resolved `genome_build`.

Complete callsite inventory requiring update:

*Source code:*
- `config.py:115` — singleton init → `GlobalConfig(snp_identifier="rsid")`
- `config.py:138` — `reset_global_config()` → same
- `__init__.py:26` — docstring example `GlobalConfig().snp_identifier` → update
- `ldscore_calculator.py:814–817` — `default_config = GlobalConfig()` fallback → remove entirely (§3.4 adds the explicit construction)
- `_kernel/annotation.py:719` — same removal

*Tests:*
- `test_config_identifiers.py:41,49,65,66` — bare `GlobalConfig()` → `GlobalConfig(snp_identifier="rsid")`
- `test_config_identifiers.py:72` — `GlobalConfig(genome_build="hg18")` → `GlobalConfig(snp_identifier="chr_pos", genome_build="hg18")`
- `test_config_identifiers.py:74` — `GlobalConfig(log_level="trace")` → `GlobalConfig(snp_identifier="rsid", log_level="trace")`
- `test_config_identifiers.py:77–79,83` — `GlobalConfig(genome_build="hg37")` etc. → add `snp_identifier="chr_pos"`
- `test_config_identifiers.py:263,267` — stale; `GlobalConfig(ref_panel_snps_file=...)` / `GlobalConfig(regression_snps_file=...)` — fields no longer on `GlobalConfig`; **delete these tests**
- `test_ref_panel_builder.py` lines 739, 789, 819, 839, 865, 901, 914, 937, 975, 986 — bare `GlobalConfig()` → `GlobalConfig(snp_identifier="chr_pos", genome_build="hg38")`
- `test_annotation.py:127` — `GlobalConfig(snp_identifier="chr_pos")` missing `genome_build` → add `genome_build="hg38"`

*Tutorials/docs:*
- `tutorials/build-parquet-reference-panel-from-plink.md:247` — `GlobalConfig(log_level="INFO")` → `GlobalConfig(snp_identifier="chr_pos", genome_build="hg38", log_level="INFO")`

**Q2 — `RefPanelConfig.genome_build` removed immediately** (no deprecation).
All usage is internal to `_ref_panel_from_args()`.

**Q3 — rsid + `genome_build` set → `UserWarning`.**
`"genome_build is set but will be ignored in rsid mode."` Emitted in
`GlobalConfig.__post_init__`. The field is then forced to `None`.

**Q4 — `resolve_genome_build()` is pure; accepts `pd.DataFrame | None`.**
I/O lives in `sample_frame_from_chr_pattern()` in `src/ldsc/_chr_sampler.py`.

**Q5 — Fail immediately** with an actionable error when the sample has fewer than
`MIN_INSPECTED_REFERENCE_SNPS` (200) HM3-overlapping SNPs. No auto-expansion.
Error message: `"Insufficient overlap with HM3 reference to infer genome build
for <context> (matched N SNPs; need ≥200). Pass --genome-build hg19 or
--genome-build hg38 explicitly."`
