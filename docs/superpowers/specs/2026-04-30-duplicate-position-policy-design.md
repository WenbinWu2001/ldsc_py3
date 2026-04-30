# Duplicate-Position Policy for PLINK-Based Reference Panel Builder

**Date:** 2026-04-30
**Scope:** `src/ldsc/ref_panel_builder.py` (PLINK-based `build-ref-panel` workflow only)

---

## Problem

In `chr_pos` mode, LDSC uses only `CHR:POS` as the SNP key. Two retained
variants sharing the same emitted `CHR:POS` make downstream lookup ambiguous.
There are two distinct origins:

- **Source-build duplicates:** the source `.bim` already contains multiple
  variants at the same `CHR:POS` (e.g. multiallelic or adjacent-normalized
  sites). Detectable before liftover.
- **Target-build collisions:** distinct source positions collapse to the same
  emitted position after liftover. Cannot be anticipated by the user; only
  detectable after the liftover chain is applied.

The current PLINK builder has no explicit detection or policy for either case.

---

## Design Decisions

| Question | Decision |
|---|---|
| Scope | PLINK builder only (`ref_panel_builder.py`); TOP-LD pipeline already has its own fix |
| Default policy | `error` — fail fast, report clusters |
| Opt-in policy | `drop-all` — drop every SNP in a colliding cluster |
| Keep-one rule | Not supported; any deterministic rule is arbitrary in `chr_pos` mode |
| Cross-build consistency | Drops apply to all emitted builds when a collision is found in any build |
| Provenance | Sidecar TSV at panel root + WARNING log pointing to it |

---

## CLI Interface

New flag on `build-ref-panel`:

```
--duplicate-position-policy {error,drop-all}
    How to handle SNPs that share a CHR:POS key in any emitted build.

    error     Abort and report all duplicate clusters. (default)
    drop-all  Drop every SNP in each colliding cluster; write a provenance
              sidecar to {output_dir}/chr{chrom}_dropped.tsv.gz.
```

### Config dataclass

`ReferencePanelBuildConfig` gains one field:

```python
duplicate_position_policy: str = "error"  # "error" | "drop-all"
```

`config_from_args()` passes `args.duplicate_position_policy` into the dataclass
directly; `choices=` in `build_parser()` enforces the two legal values.

---

## Detection & Drop Logic

### New helper: `_resolve_unique_snp_set()`

```python
def _resolve_unique_snp_set(
    chrom: str,
    chrom_df: pd.DataFrame,       # source .bim rows, indexed by PLINK row index
    keep_snps: np.ndarray,         # restriction-filtered PLINK indices
    hg19_lookup: dict[int, int],   # idx -> hg19 position (may be empty)
    hg38_lookup: dict[int, int],   # idx -> hg38 position (may be empty)
    policy: str,                   # "error" | "drop-all"
) -> tuple[np.ndarray, pd.DataFrame]:
```

**Returns:** cleaned `keep_snps` array + provenance `DataFrame` (empty if nothing
dropped). The DataFrame has columns `CHR, SNP, source_pos, target_pos, reason`.
One row per dropped SNP (not per cluster): a two-SNP collision produces two
provenance rows.

**Two detection passes, executed in order:**

1. **Source-build duplicates** — check uniqueness of `chrom_df.loc[keep_snps, "BP"]`.
   Variants sharing a source `CHR:POS` form a cluster. `target_pos` is `pd.NA`
   in provenance rows.

2. **Target-build collisions** — for each non-empty lookup (hg19, hg38) — empty
   lookups occur in source-only builds and are skipped — check uniqueness of
   position values among the surviving `keep_snps`. A SNP is a collider if its
   translated position is shared by any other retained SNP in either lookup.
   Colliders are removed from both lookups to enforce cross-build consistency.

**Policy application (same for both passes):**

- `"error"`: raise `ValueError` with a formatted table of all duplicate clusters
  (CHR, positions, SNP IDs). No sidecar is written.
- `"drop-all"`: remove all SNPs in any cluster from `keep_snps`; append
  provenance rows.

---

## Sidecar File

**Path:** `{output_dir}/chr{chrom}_dropped.tsv.gz` — panel root, not under any
`{build}/` subdirectory, because drops apply to all emitted builds.

**Columns:**

| column | type | notes |
|---|---|---|
| `CHR` | str | chromosome label, e.g. `"1"` |
| `SNP` | str | rsID from .bim |
| `source_pos` | int | position in source build |
| `target_pos` | Int64 (nullable) | lifted-over position; `NA` for source-build duplicates |
| `reason` | str | `"source_duplicate"` or `"target_collision"` |

Written with `df.to_csv(path, sep="\t", index=False, compression="gzip")`,
matching the convention of `chr{chrom}_meta.tsv.gz`.

**Written only when** provenance DataFrame is non-empty and policy is `"drop-all"`.
If nothing is dropped, no sidecar is written.

**Log message** (WARNING level):

```
Dropped {n} SNPs on chromosome {chrom} due to duplicate positions
({k} source-build duplicates, {m} target-build collisions).
Provenance written to '{output_dir}/chr{chrom}_dropped.tsv.gz'.
```

---

## Integration into `_build_chromosome()`

One insertion point, after `_resolve_mappable_snp_positions()` returns and before
the per-build emit loop:

```python
keep_snps, hg19_lookup, hg38_lookup = self._resolve_mappable_snp_positions(...)

keep_snps, dropped_df = _resolve_unique_snp_set(
    chrom=chrom,
    chrom_df=chrom_df,
    keep_snps=keep_snps,
    hg19_lookup=hg19_lookup,
    hg38_lookup=hg38_lookup,
    policy=config.duplicate_position_policy,
)
if not dropped_df.empty:
    _write_dropped_sidecar(dropped_df, Path(config.output_dir) / f"chr{chrom}_dropped.tsv.gz")

if len(keep_snps) == 0:
    LOGGER.info(f"Skipping chromosome {chrom}: no SNPs remain after duplicate-position filtering.")
    return None
```

`_write_dropped_sidecar()` is a small private helper: writes the TSV and emits
the WARNING log message.

`_expected_ref_panel_output_paths()` does **not** include the sidecar — it is
conditional and optional, so it is not pre-checked in the output-path preflight.

The per-build emit loop receives the already-cleaned `keep_snps` and both
lookups with no further changes.

---

## Testing Strategy

### Unit tests for `_resolve_unique_snp_set()`

| scenario | expected behaviour |
|---|---|
| All positions unique, both builds | `keep_snps` unchanged, empty DataFrame |
| Source-build duplicate, `"error"` | raises `ValueError` listing the cluster |
| Source-build duplicate, `"drop-all"` | removes both SNPs; `reason="source_duplicate"`, `target_pos=NA` |
| Target-build collision (hg38 only), `"error"` | raises `ValueError` |
| Target-build collision (hg38 only), `"drop-all"` | removes both SNPs from `keep_snps` (cross-build); `reason="target_collision"` |
| Collision in both builds simultaneously | all colliders removed; both reasons appear in provenance |
| All SNPs dropped | returns empty `keep_snps`, non-empty DataFrame |

### Integration tests via `ReferencePanelBuilder.run()`

- `--duplicate-position-policy error` with a `.bim` containing source-build
  duplicates → `ValueError` raised, no output files written.
- `--duplicate-position-policy drop-all` with target-build collisions → sidecar
  written at panel root, both build parquets contain only unique positions, log
  includes the WARNING with sidecar path.
