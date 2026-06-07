# Chromosome-Region Exclusion for `build-ref-panel` and `ldscore`

**Date:** 2026-06-06
**Status:** Approved design (pre-implementation)
**Topic:** `region-exclusion`
**Companion plan:** `docs/superpowers/plans/2026-06-06-region-exclusion-plan.md` (to be written)

---

## 1. Motivation

Certain genomic regions carry LD structure that distorts LD-score regression and
downstream heritability/genetic-correlation estimates: the MHC/HLA block on
chr6 (extreme, long-range LD), and centromeric assembly gaps (sparse, unreliable
genotyping). Standard practice is to **remove all SNPs in these regions before
any LD computation**, so their LD never enters the emitted artifacts.

This feature adds a way to name such regions and have every SNP inside them
dropped:

- in **`build-ref-panel`**, before pairwise-R2 is computed and written to the
  reference-panel parquet artifacts; and
- in **`ldscore`**, before LD scores are computed from either the PLINK or the
  parquet-R2 reference backend.

The two modules are independent runs; the feature is configured per run (not via
`GlobalConfig`), mirroring the existing per-run `ref_panel_snps_file`
restriction.

## 2. Scope

### In scope
- Two curated, named region **presets**: `mhc` and `centromeres`, packaged as
  standard per-build BED files.
- A **user-supplied BED** escape hatch for arbitrary region exclusion.
- Wiring both into the existing single SNP-restriction chokepoint in each module.
- Build-resolution rules (see §6), validation, logging, and tests.

### Out of scope (deferrable without rework — see §11)
- Additional presets (long-range-LD/Price 2008, telomeres): future rows of
  curated data, no new code paths.
- Inline CLI interval syntax (e.g. `--exclude-region 6:25000000-34000000`).
- Recording region-excluded SNPs in the `build-ref-panel` dropped-SNP sidecar
  (consistent with how `ref_panel_snps_file` drops are already kept out of that
  audit vocabulary; see `ref_panel_builder.py` `_build_chromosome` docstring).

## 3. Design principles & rationale

The design rests on three observations about the existing code:

1. **Both modules already have exactly one SNP-restriction chokepoint keyed on
   CHR/POS.** Region exclusion is structurally an *exclude-by-interval* filter —
   the dual of the existing *keep-by-identity* `ref_panel_snps_file` filter.
   Reusing those seams means **no numerical kernel changes**: excluded SNPs never
   reach R2 streaming or LD accumulation, so zero compute is spent on them.

2. **In `build-ref-panel`, coordinates are in the *source* build at the filter
   point**, before liftover. Excluding on source-build positions therefore
   removes those SNPs from *every* emitted build with a single region set. This
   matches the module's existing "synchronized cross-build drop" behavior, where
   a SNP dropped in any emitted build is dropped from all.

3. **In `ldscore`, every path funnels through `RefPanel.load_metadata` →
   `_apply_snp_restriction`** (PLINK backend, parquet backend, the synthetic
   all-ones `base` annotation, and annotation alignment). One filter there covers
   the whole module.

A region preset is **sugar over an interval list**: presets and user BEDs both
reduce to `(chrom, start, end)` intervals fed to one mask engine. This is why
adding future presets requires only curated data, not code.

## 4. Region data files

### 4.1 Layout

Standard 3-column BED, one file **per preset per build**, under
`src/ldsc/data/regions/`:

```
src/ldsc/data/regions/mhc.hg19.bed
src/ldsc/data/regions/mhc.hg38.bed
src/ldsc/data/regions/centromeres.hg19.bed
src/ldsc/data/regions/centromeres.hg38.bed
```

Rationale for four standard BEDs over two files with a `build` column: a `build`
column produces a bespoke schema only our loader understands, whereas plain BED
stays inspectable/validatable in IGV, bedtools, etc. — these are coordinate
tables a human will re-curate against UCSC. The loader becomes a one-liner with
no row-filtering branch, and a missing/typo'd build fails loudly (file not found)
rather than silently dropping a build. The cost is two extra ~1 KB files.

### 4.2 Format

- **Standard BED**: tab-delimited, **0-based, half-open** `[start, end)`.
- Columns: `chrom  start  end` (a 4th `name` column is permitted and ignored).
- No header line (BED convention). Lines beginning with `#`, `track`, or
  `browser` are skipped.
- `chrom` uses bare tokens (`6`, `X`), but the loader normalizes any `chr`
  prefix and `X`/`Y` casing via `normalize_chromosome` so user-edited files with
  `chr6` still match.

### 4.3 Contents

- **`mhc.*.bed`** — single MHC interval on chr6. Coordinates differ by build:
  - hg19: `6  25000000  35000000` (the conventional broad MHC exclusion window).
  - hg38: `6  28477797  33448354` (lifted MHC core; final exact bounds fixed
    during implementation against the curation source).
- **`centromeres.*.bed`** — one interval per chromosome (autosomes + X),
  sourced from the UCSC `gap` table (`type = centromere`) for the matching build
  (~24 rows). Exact rows are materialized during implementation from UCSC and
  recorded with their provenance in a comment block at the top of each file.

> **Implementation note:** the precise centromere/MHC coordinates and their
> UCSC source URLs/dates are pinned in the implementation plan and embedded as
> `#`-comment provenance headers in each BED file so future curators can
> reproduce them.

### 4.4 Packaging

`setup.py` currently declares `package_data={"ldsc": ["data/*.tsv.gz"]}`. Extend
to include the BEDs:

```python
package_data={"ldsc": ["data/*.tsv.gz", "data/regions/*.bed"]},
```

Files are resolved at runtime exactly like the HM3 map (`hm3.py`
`packaged_hm3_curated_map_path`):

```python
from importlib import resources
resources.files("ldsc").joinpath("data", "regions", f"{name}.{build}.bed")
```

## 5. Kernel helper — `src/ldsc/_kernel/regions.py` (new, private)

A small, dependency-light module owning interval loading and masking. It only
ever sees primitive inputs (names, a build, a metadata frame) per the
public-workflow/private-kernel split.

### 5.1 Public-to-kernel surface

```python
# Preset menu — single source of truth for valid names.
REGION_PRESETS: frozenset[str] = frozenset({"mhc", "centromeres"})

@dataclass(frozen=True)
class RegionIntervals:
    """Per-chromosome exclusion intervals in 0-based half-open BED coordinates.

    intervals : dict[str, np.ndarray]
        Maps normalized chromosome token -> (k, 2) int64 array of [start, end)
        rows, sorted by start. Empty mapping means "exclude nothing".
    source_labels : tuple[str, ...]
        Human-readable provenance for logging, e.g. ("preset:mhc[hg19]",
        "bed:/path/blacklist.bed").
    """
    intervals: dict[str, np.ndarray]
    source_labels: tuple[str, ...]

def load_preset_intervals(names: Sequence[str], build: str) -> RegionIntervals:
    """Read packaged BEDs for the requested presets in the given build.

    Validates names against REGION_PRESETS and build against {"hg19","hg38"};
    raises LDSCUsageError / LDSCConfigError with self-contained messages
    otherwise. Unions intervals across presets.
    """

def load_bed_intervals(paths: Sequence[str | PathLike[str]]) -> RegionIntervals:
    """Read one or more user BED files, applied as-is on panel CHR/POS.

    Minimal validation: >= 3 whitespace/tab columns, integer start/end,
    start < end. Raises LDSCInputError with the offending file + line on
    malformed rows. Build-agnostic by construction (no translation).
    """

def merge_intervals(*groups: RegionIntervals) -> RegionIntervals:
    """Union several RegionIntervals (presets + user BEDs) into one."""

def region_exclusion_keep_mask(
    metadata: pd.DataFrame,
    intervals: RegionIntervals,
    *,
    chr_col: str = "CHR",
    pos_col: str = "POS",
) -> np.ndarray:
    """Return a boolean KEEP mask over metadata rows.

    A 1-based SNP POS is EXCLUDED iff some interval [start, end) (0-based
    half-open) contains it, i.e. start < POS <= end. Vectorized per chromosome
    via np.searchsorted on sorted interval bounds. Rows whose chromosome has no
    intervals are always kept. Returns all-True when intervals is empty.
    """
```

### 5.2 Coordinate convention (the one subtle invariant)

Panel metadata `POS`/`BP` is **1-based** (PLINK convention). BED intervals are
**0-based half-open**. A 1-based position `p` lies in a BED interval
`[start, end)` iff:

```
start <= (p - 1) < end      ⟺      start < p <= end
```

`region_exclusion_keep_mask` implements the right-hand inequality directly. This
must be covered by explicit boundary tests (a SNP exactly at `start`, at
`start+1`, at `end`, at `end+1`).

## 6. Build-resolution rules (the core asymmetry)

Presets are build-specific; user BEDs are build-agnostic. Each module resolves
the **preset** build differently, and this asymmetry is intentional and safe.

| Module | Preset build source | User-BED build |
|---|---|---|
| `build-ref-panel` | The already-resolved `source_genome_build` (concrete `hg19`/`hg38` by the filter point). **No build flag.** | Applied as-is on source-build CHR/POS. |
| `ldscore` | An explicit `--exclude-regions-build {hg19,hg38}`, **required whenever any preset is named.** | Applied as-is on panel CHR/POS. |

### Why `build-ref-panel` needs no `--exclude-regions-build`

`ReferencePanelBuilder._run` calls `_resolve_source_genome_build`
(`ref_panel_builder.py:425`) *before* `_prepare_build_state`
(`ref_panel_builder.py:426`), so `source_genome_build` is always a concrete
`hg19`/`hg38` when intervals load. That resolved build **is** the coordinate
system every panel SNP is addressed in at the filter point, so it is the only
correct preset build — there is no second build to disambiguate. Adding a flag
would *introduce* a contradiction state (`--source-genome-build hg19
--exclude-regions-build hg38`) we would then have to detect and reject. Omitting
it makes that state unrepresentable.

**Where auto-build resolution happens.** The preset build *can* be auto-resolved,
but only as a side effect of the panel's own source-build inference — the preset
itself performs no inference. When `--source-genome-build auto` (the default),
`_resolve_source_genome_build` infers `hg19`/`hg38` from the PLINK `.bim`
coordinates (HM3-overlap coordinate inference). Because that runs before
`_prepare_build_state`, `load_preset_intervals(..., source_build)`
(`ref_panel_builder.py:700`) always receives a concrete build. So
`build-ref-panel` is the *only* place an auto-inferred build selects a preset
BED, and it inherits whatever the panel resolved to — there is no separate
region-build inference step.

### Why `ldscore` does need it

`ldscore` may address the panel by rsID (`rsid` / `rsid_allele_aware`
identifier modes), in which case `GlobalConfig.genome_build` is `None` and no
build is in scope — yet the panel's `POS` column is still in *some* build. A
coordinate preset cannot know which packaged BED to apply without being told.
We therefore **always require** `--exclude-regions-build` when a preset is named
(even in `chr_pos` modes where a build is otherwise known), so the safe,
explicit path is the only path and we never silently rely on an auto-inferred
build that might disagree with the panel. The flag accepts `hg19`/`hg38` only —
there is deliberately **no `auto` choice** — so on the `ldscore` side the preset
build is never inferred; the caller always states it. `load_preset_intervals`
enforces this, raising `LDSCConfigError` for any build outside `{hg19, hg38}`.

User BEDs are exempt: the user owns those coordinates and supplies them in the
panel's build, the same trust model as `--ref-panel-snps-file`. The `-build`
flag governs **preset selection only**.

## 7. Configuration changes (`src/ldsc/config.py`)

### 7.1 `ReferencePanelBuildConfig` (`config.py:513`)

Add:
```python
exclude_regions: tuple[str, ...] = ()        # subset of REGION_PRESETS
exclude_regions_bed: tuple[str, ...] = ()     # user BED path tokens
```
`__post_init__` additions:
- Normalize each `exclude_regions_bed` token via `_normalize_optional_path`
  (drop empties).
- Validate `set(exclude_regions) <= REGION_PRESETS`; raise `LDSCConfigError`
  naming the unknown preset(s) and listing the valid menu.
- No build field — `source_genome_build` supplies it.

### 7.2 `RefPanelConfig` (`config.py:360`)

Add:
```python
exclude_regions: tuple[str, ...] = ()
exclude_regions_bed: tuple[str, ...] = ()
exclude_regions_build: Literal["hg19", "hg38"] | None = None
```
`__post_init__` additions:
- Same preset-name and path normalization as above.
- If `exclude_regions` is non-empty and `exclude_regions_build is None`, raise
  `LDSCConfigError` instructing the user to pass `--exclude-regions-build`.
- Validate `exclude_regions_build in {None, "hg19", "hg38"}`.

> `LDScoreConfig` is **not** touched: region exclusion narrows the *reference
> panel* universe, which `RefPanelConfig` owns (it already owns
> `ref_panel_snps_file`). Keeping it there means the existing
> `_apply_snp_restriction` seam is the natural home.

## 8. CLI changes

Both parsers gain `--exclude-regions` and `--exclude-regions-bed`; only
`ldscore` gains `--exclude-regions-build`. Multi-valued flags accept a
comma-separated list, reusing the existing `split_cli_path_tokens` idiom for
BED paths.

### 8.1 `build-ref-panel` (`ref_panel_builder.py` `build_parser`)

```
--exclude-regions {mhc,centromeres}[,...]
    Named curated regions to exclude (drop all SNPs inside them) before R2 is
    computed. Resolved in the panel's source genome build; the same SNPs are
    removed from every emitted build.
--exclude-regions-bed PATH[,PATH...]
    User BED file(s) of regions to exclude, applied as-is on source-build
    CHR/POS (0-based half-open). No build translation is performed.
```
Mapped in `config_from_args` (`ref_panel_builder.py:1798`) into the two new
fields.

### 8.2 `ldscore` (`ldscore_calculator.py` `build_parser`)

```
--exclude-regions {mhc,centromeres}[,...]
    Named curated regions to exclude before LD scores are computed. Requires
    --exclude-regions-build.
--exclude-regions-build {hg19,hg38}
    Genome build of the panel coordinates, used to select the preset BEDs.
    Required whenever --exclude-regions is given.
--exclude-regions-bed PATH[,PATH...]
    User BED file(s) of regions to exclude, applied as-is on panel CHR/POS
    (0-based half-open).
```
Mapped through `_normalize_run_args` / `_ref_panel_from_args`
(`ldscore_calculator.py:1432`) into `RefPanelConfig`.

`choices` validation in argparse enumerates the preset menu in `--help`,
covering the discoverability gap of a single menu flag versus separate booleans.

## 9. Wiring (the two chokepoints)

### 9.1 `build-ref-panel`

- **Load intervals once** in `_prepare_build_state` (`ref_panel_builder.py:555`),
  where `config.source_genome_build` is already concrete. Merge presets (in
  source build) with any user BEDs into a single `RegionIntervals`, store on
  `_BuildState` (new field `region_intervals: RegionIntervals`). Log the
  resolved `source_labels`.
- **Apply per chromosome** in `_build_chromosome` (`ref_panel_builder.py:734`),
  folded into the existing keep logic next to `build_restriction_mask`
  (`ref_panel_builder.py:789`): compute `region_exclusion_keep_mask` on
  `chrom_metadata` (which has `CHR`/`POS`) and intersect with `keep_snps`. This
  runs **before** identity cleanup and liftover, on source coordinates.
- Emit an INFO log with the per-source dropped count for the chromosome. If a
  chromosome empties out, reuse the existing "no SNPs remain after …" skip path.

### 9.2 `ldscore`

- **Apply** inside `RefPanel._apply_snp_restriction`
  (`_kernel/ref_panel.py:267`), immediately after the existing keep-restriction,
  so the single seam covers both backends, the synthetic `base` path, and
  annotation alignment.
- `RefPanel` (or its `spec`/`RefPanelConfig`) already carries the three new
  fields; intervals are loaded lazily on first `load_metadata` call (cache on
  the instance) using `exclude_regions_build` for presets and the raw paths for
  user BEDs.
- Emit an INFO log with per-source dropped counts (aggregated across
  chromosomes, or per chromosome — match the surrounding logging granularity).

## 10. Diagnostics & error handling

- **Logging:** both modules log resolved region sources and dropped SNP counts
  at INFO. Example: `Excluded 1,832 SNPs via preset:mhc[hg19] on chromosome 6.`
- **Empty result:** if exclusion removes every SNP on a chromosome (or all
  chromosomes), the existing empty-intersection / empty-build error paths fire
  unchanged, now also citing region exclusion as a possible cause.
- **Errors (self-contained, per repo convention — what + where + likely cause +
  remedy):**
  - Unknown preset name → `LDSCConfigError` listing the valid menu.
  - `ldscore` preset without `--exclude-regions-build` → `LDSCConfigError`
    telling the user to add the flag.
  - Malformed user BED row → `LDSCInputError` naming the file and line.
  - Missing packaged BED for a `(name, build)` pair → `LDSCInternalError`
    (packaging bug, not user error).
- No new `docs/troubleshooting.md` section is required unless an error grows to
  3+ distinct causes; the single-cause messages above are self-contained.

## 11. Extensibility

Because the engine only consumes `(chrom, start, end)` intervals:
- **New presets** (long-range-LD, telomeres) = add `name.hg19.bed` /
  `name.hg38.bed` and one entry to `REGION_PRESETS`. No new flags, config
  fields, or wiring.
- **Inline ranges** or **per-build user BEDs** could later feed the same
  `RegionIntervals` without touching the mask engine or chokepoints.

## 12. Testing plan (TDD — tests written first)

Numerical/behavioral tests, not just "runs without error":

1. **Mask correctness (`regions.py`):** known intervals, with explicit boundary
   cases on the half-open edge — SNP at `start` (kept), `start+1` (excluded),
   `end` (excluded), `end+1` (kept). Multi-interval and multi-chromosome cases.
2. **Preset loading:** `load_preset_intervals(["mhc"], "hg19")` and `"hg38"`
   return the expected chr6 interval; `centromeres` returns one interval per
   chromosome; unknown name and unknown build raise.
3. **User BED parsing:** `chr6`→`6` normalization; malformed row (2 columns,
   non-integer, `start >= end`) raises with file+line.
4. **`ldscore` build-flag guard:** preset without `--exclude-regions-build`
   raises; with it, the correct preset build is selected.
5. **`build-ref-panel` two-build emission:** a SNP inside the source-build MHC is
   absent from **both** the hg19 and hg38 emitted parquet/metadata artifacts
   (small synthetic PLINK + liftover fixture).
6. **Integration (both modules):** excluded SNPs do not appear in the emitted
   artifacts, and LD scores for retained SNPs are unaffected by the exclusion of
   distant excluded SNPs outside their LD window.

## 13. Affected files (implementation checklist)

| File | Change |
|---|---|
| `src/ldsc/data/regions/{mhc,centromeres}.{hg19,hg38}.bed` | New curated BED data (4 files) |
| `setup.py` | Extend `package_data` glob to `data/regions/*.bed` |
| `src/ldsc/_kernel/regions.py` | New: interval loading + keep-mask engine |
| `src/ldsc/config.py` | New fields + validation on `ReferencePanelBuildConfig`, `RefPanelConfig` |
| `src/ldsc/ref_panel_builder.py` | Parser flags, `config_from_args`, `_BuildState`, `_prepare_build_state`, `_build_chromosome` wiring |
| `src/ldsc/ldscore_calculator.py` | Parser flags, `_normalize_run_args`/`_ref_panel_from_args` wiring |
| `src/ldsc/_kernel/ref_panel.py` | Apply mask in `_apply_snp_restriction` |
| `tests/…` | New tests per §12 |
| `docs/current/…`, `design_map.md` | Update if module structure/data-flow docs reference the restriction chokepoints |

---

*End of design.*
