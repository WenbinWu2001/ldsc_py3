# Chromosome-Region Exclusion for `build-ref-panel` and `ldscore`

**Date:** 2026-06-06 (revised post-implementation)
**Status:** Implemented; §6 and §8 revised to match the shipped behavior.
**Topic:** `region-exclusion`
**Companion plan:** `docs/superpowers/plans/2026-06-06-region-exclusion-plan.md`

> **Evolution since the original design.** Three things changed after the first
> implementation; the sections below reflect the current code:
> - `--exclude-regions` is a **fixed single-choice enum**
>   (`none`, `mhc`, `centromeres`, `mhc-and-centromeres`) via
>   `EXCLUDE_REGIONS_CHOICES`, **defaulting to `mhc-and-centromeres` (exclusion ON
>   by default)** — not a free comma list (§8).
> - `--exclude-regions-build` is **inferred** from the panel build in chr_pos
>   modes and required only in rsID modes — not "always required" (§6).
> - The active `mhc` preset is the **broad chr6:25-35Mb** window and the
>   `centromeres` preset is the **pericentromeric ±3 cM** region (LDSC parity);
>   the narrower `mhc_core` / `centromeres_core` definitions are preserved as
>   Python-API-only reference presets (§4.3). See also
>   `docs/current/region-exclusion-presets.md` for the full preset/coordinate
>   table.

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
src/ldsc/data/regions/mhc.hg19.bed                  # active: broad chr6:25-35Mb
src/ldsc/data/regions/mhc.hg38.bed                  # active: broad chr6:25-35Mb
src/ldsc/data/regions/mhc_core.hg19.bed             # classical-HLA core, reference
src/ldsc/data/regions/mhc_core.hg38.bed             # classical-HLA core, reference
src/ldsc/data/regions/centromeres.hg19.bed          # active: pericentromeric +/-3 cM
src/ldsc/data/regions/centromeres.hg38.bed          # active: pericentromeric +/-3 cM
src/ldsc/data/regions/centromeres_core.hg19.bed     # raw gap, reference
src/ldsc/data/regions/centromeres_core.hg38.bed     # raw gap, reference
```

> The full preset/build/coordinate/provenance table is maintained at
> `docs/current/region-exclusion-presets.md`.

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

- **`mhc.*.bed`** (the **active** `mhc` preset) — single broad GWAS exclusion
  window on chr6, **build-consistent**: `6  25000000  35000000` in **both**
  builds (hg19 25-35Mb lifts to hg38 ~25.00-35.03Mb, so the same numeric window
  covers the same biology). Pinned constant.
- **`mhc_core.*.bed`** — the narrow **classical-HLA core**, reference only (not
  wired to any CLI choice): the GWAS-Catalog / GRC MHC region in each build's
  native coordinates — hg19 `6  28477797  33448354`, GRCh38
  `6  28510120  33480577`.
- **`centromeres.*.bed`** (the **active** `centromeres` preset) — one interval
  per chromosome (autosomes + X), the **pericentromeric region: each centromere
  span padded by ±3 cM**, matching the LD Score regression exclusion of
  Bulik-Sullivan et al. 2015 Nat Genet (Online Methods). Padding is computed at
  curation time: cM at the centromere endpoints is read from the Alkes-group
  recombination map (`genetic_map_{build}_withX`), extended ±3 cM, and inverted
  back to base pairs. The padded interval is **unioned with the centromere span**
  so it always contains it — this also gives the correct conservative behavior on
  the acrocentric chromosomes (13/14/15/21/22), whose p-arm is absent from the
  genetic map. The genetic maps are a curation-time input only (workspace
  `resources/`, overridable via `LDSC_GENETIC_MAP_DIR`), not bundled.
- **`centromeres_core.*.bed`** — the **raw** UCSC centromere span (hg19 `gap`
  track `type = centromere`; hg38 `centromeres` track merged per chromosome).
  Kept for reference and loadable via the `centromeres_core` preset (Python API),
  but **not wired to any CLI choice**. Retained for later expert review of the
  centromere-vs-pericentromere definition.

> **Provenance.** Each BED's first `#` line records its source (UCSC track or the
> Alkes-group map + the 3 cM parameter). All four are regenerable by
> `tools/regions/build_region_beds.py`.

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

## 6. Build-resolution rules

> **Updated to current behavior.** The original design required
> `--exclude-regions-build` whenever a preset was named. The shipped behavior
> *infers* it where possible; this section documents what the code does.

Presets are build-specific; user BEDs are build-agnostic. Each module resolves
the **preset** build differently.

| Module | Preset build source | User-BED build |
|---|---|---|
| `build-ref-panel` | The already-resolved `source_genome_build` (concrete `hg19`/`hg38` by the filter point). **No build flag.** | Applied as-is on source-build CHR/POS. |
| `ldscore` | `--exclude-regions-build` if given; otherwise **inferred** from the panel build in chr_pos modes; **required** in rsID modes. | Applied as-is on panel CHR/POS. |

### `ldscore` — `--exclude-regions-build` is optional, resolved by `_resolve_exclude_regions_build`

`ldscore_calculator._resolve_exclude_regions_build(args, global_config, presets)`
fills the build *before* `RefPanelConfig` is constructed, in this order:

1. **Explicit flag wins** — if `--exclude-regions-build` is passed, use it.
2. **No presets active** (`--exclude-regions none`, or only a user BED) — build
   stays `None`; nothing to resolve.
3. **chr_pos-family identifier modes** — reuse `global_config.genome_build` (the
   build `ldscore` already operates the panel in). If that is not a concrete
   `hg19`/`hg38`, raise `LDSCUsageError` asking for an explicit
   `--exclude-regions-build` or `--exclude-regions none`.
4. **rsID-family modes** — `global_config.genome_build` is `None` (rsID
   identifiers carry no build), so an explicit `--exclude-regions-build` is
   required; otherwise raise `LDSCUsageError`.

`RefPanelConfig` still records a concrete build for presets (its `__post_init__`
rejects presets with a `None` build), but the *workflow* supplies that build, so
the user only has to type the flag when it cannot be inferred (rsID modes, or an
unresolved panel build). `load_preset_intervals` enforces the final
`{hg19, hg38}` domain, raising `LDSCConfigError` otherwise.

### `build-ref-panel` — no `--exclude-regions-build`

`ReferencePanelBuilder._run` calls `_resolve_source_genome_build` *before*
`_prepare_build_state`, so `source_genome_build` is always a concrete
`hg19`/`hg38` when intervals load — and that resolved build is the coordinate
system every panel SNP is addressed in at the filter point. The preset performs
no inference of its own: when `--source-genome-build auto` (the default), the
panel's own `.bim`-based inference produces the concrete build that
`load_preset_intervals(..., source_build)` then receives. Adding a separate
region-build flag could only *introduce* a contradiction
(`--source-genome-build hg19 --exclude-regions-build hg38`), so it is omitted.

User BEDs are exempt in both modules: the user owns those coordinates and
supplies them in the panel's build, the same trust model as
`--ref-panel-snps-file`. The build resolution governs **preset selection only**.

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

> **Updated to current behavior.** `--exclude-regions` is a single-choice enum
> (not a free comma list) and **defaults to `mhc-and-centromeres` — exclusion is
> ON by default**. The choice → preset-name expansion lives in
> `regions.EXCLUDE_REGIONS_CHOICES` / `exclude_regions_choice_to_presets`.

The vocabulary is shared by both modules:

```
--exclude-regions {none,mhc,centromeres,mhc-and-centromeres}
    Curated regions to exclude before LD computation / R2 emission.
    Default: mhc-and-centromeres. Use 'none' to keep all regions.
    `mhc` is the broad chr6:25-35Mb window; `centromeres` is the
    pericentromeric +/-3 cM region (LDSC). The `*_core` reference presets
    (mhc_core, centromeres_core) are Python-API only and not offered here.
--exclude-regions-bed PATH[,PATH...]
    User BED file(s) to exclude, applied as-is on panel/source CHR/POS
    (0-based half-open). No build translation.
```

### 8.1 `build-ref-panel`

`--exclude-regions` + `--exclude-regions-bed`; **no `--exclude-regions-build`**
(reuses `source_genome_build`). `config_from_args` maps the choice via
`exclude_regions_choice_to_presets` and `split_cli_path_tokens`.

### 8.2 `ldscore`

`--exclude-regions` + `--exclude-regions-bed`, plus an **optional**
`--exclude-regions-build {hg19,hg38}` (see §6: inferred in chr_pos modes,
required in rsID modes when presets are active). `_ref_panel_from_args` maps the
choice and calls `_resolve_exclude_regions_build` before constructing
`RefPanelConfig`.

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
