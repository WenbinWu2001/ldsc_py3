# LD Window Behavior for Parquet R2 Panels

Last updated on: 2026-10-27

This document records how `ldsc ldscore` interprets LD-window flags when the
reference panel backend is the canonical index-format parquet R2 pair.

## Paired artifacts

A complete package-built parquet R2 chromosome has exactly two artifacts, both
mandatory:

```text
chr{chrom}_r2.parquet   — 4-column index format: IDX_1, IDX_2, R2, SIGN_R
chr{chrom}_meta.tsv.gz  — per-SNP metadata sidecar (CHR, POS, SNP, A1, A2, CM, MAF)
```

The R2 parquet stores pairwise LD as sidecar-row indices — it carries no SNP
identity of its own and is meaningless without the exact matching sidecar. The
metadata sidecar defines the full panel SNP universe and the index space.

The sidecar format:

```text
# ldsc:artifact_type=ref_panel_metadata
# ldsc:snp_identifier=chr_pos_allele_aware
# ldsc:genome_build=hg38
CHR  POS  SNP  A1  A2  CM  MAF
```

At LD-score runtime, pairwise R2 values come from `chr{chrom}_r2.parquet`.
Window coordinates, retained SNP order, identifiers, `CM`, and `MAF` come from
the runtime metadata table loaded from `chr{chrom}_meta.tsv.gz`.

## Sidecar is mandatory

The metadata sidecar is a hard requirement. If it is absent, `ParquetR2RefPanel`
raises immediately with an actionable message — there is no synthesize-from-endpoints
fallback. Regenerate the panel with `ldsc build-r2-panel` to produce a paired
parquet + sidecar.

## LD-window flag behavior

| Flag | Coordinate source | Behavior |
|---|---|---|
| `--ld-wind-snps` | sidecar row positions | Counts panel rows, not the SNPs this run keeps. |
| `--ld-wind-kb` | sidecar `POS` | Normal operation. |
| `--ld-wind-cm` | sidecar `CM` | Normal when `CM` is non-missing. Fails when any retained SNP has missing `CM`. |

The shared LD-window resolver requires exactly one of `--ld-wind-snps`,
`--ld-wind-kb`, or `--ld-wind-cm`.

### Window rules against the panel's build window

A panel stores only the SNP pairs inside the window it was built with, recorded as `ldsc:ld_window_mode` and `ldsc:ld_window_value`.
`ldsc ldscore` therefore accepts only windows the panel is known to cover from those settings:

| Requested window | Result |
|---|---|
| Same unit, no wider than the build window | Accepted. |
| Same unit, wider than the build window | `LDSCUsageError` ("is wider than the input R2 parquet panel window"). |
| Different unit (for example a kb panel and `--ld-wind-snps`) | `LDSCUsageError` ("windows in different units cannot be checked for coverage"). |

The CLI checks every selected chromosome file in the header preflight, before pair scans, annotation preparation, or any LD-score computation.
`ParquetR2RefPanel.prepare_chromosome` repeats the check per chromosome because the `LDScoreCalculator` Python API skips that preflight.
A different-unit window is rejected even when the panel happens to cover it, because coverage is not inferred from per-SNP coordinates.
See `validate_ldscore_window_within_r2_panel_window` in `src/ldsc/_kernel/ldscore.py` and `inspect_r2_headers` in `src/ldsc/_input_preflight.py`.

### `--ld-wind-snps` counts panel rows

For a parquet panel, `--ld-wind-snps K` spans `K` sidecar rows on each side, whatever SNPs the run keeps.
A window no wider than the panel's then never needs a pair the panel did not store.
When the run keeps fewer SNPs than the panel (ldscore `--maf-min`, `--ref-panel-snps-file`, or annotations covering part of the panel), its SNPs have fewer than `K` contributors per side.
The PLINK backend counts the SNPs the run keeps, as LDSC2 did, so the two backends then differ for the same flags.
Each such run logs one warning with the number of affected chromosomes, one example with its kept and panel SNP counts, and the location of the per-chromosome counts.
Every parquet `--ld-wind-snps` run writes `diagnostics/ld_window_snp_counts.tsv` with `CHR`, `panel_snps`, `kept_snps`, `ld_wind_snps`, and `counting_basis` (`panel_rows`), listed in `metadata.json` `files`.
Kept SNPs must be in sidecar row order; several SNPs at one position listed in a different order from the sidecar stop the run with an `LDSCInputError`.
See `build_window_coordinates` in `src/ldsc/_kernel/ldscore.py` and `_ld_window_snp_counts` in `src/ldsc/ldscore_calculator.py`.

**`CM`/`MAF` are always reference-panel-sourced.** Annotation `CM`/`MAF` are
ignored; for the parquet backend the sidecar is authoritative, and the PLINK
backend uses `.bim` `CM` (or an interpolated genetic map) and genotype `MAF`.
For PLINK, `--ld-wind-cm` additionally rejects an *unusable* `.bim` `CM` (fewer
than two distinct finite values per chromosome) with a dedicated error that
`--yes-really` does not bypass; supply `--genetic-map-hg19-sources` /
`--genetic-map-hg38-sources` to interpolate `CM` at the `.bim` positions. See
`docs/troubleshooting.md#ldscore-unusable-cm-for-ld-wind-cm`.

For `--ld-wind-cm`, the check is strict: every retained SNP must have a
non-missing `CM`. If any retained row has missing `CM`, chromosome computation
raises:

```text
--ld-wind-cm requires non-missing CM values for all retained SNPs.
```

Package-built sidecars from source-only runs (no genetic map supplied) write
`CM=NA` for every SNP.
Such a panel was built with a kb or SNP-count window, so the window-unit rule above rejects `--ld-wind-cm` before this check; use the panel's own unit, or rebuild the panel with a genetic map and `--ld-wind-cm`.

## Annotation interaction

When no baseline annotations are supplied, `ldscore` synthesizes an all-ones
`base` annotation from `ref_panel.load_metadata()`.

When baseline annotations are supplied, the reference-panel metadata is used to
intersect annotation rows with the retained reference-panel universe. For
LD-score calculation, the **sidecar is authoritative** for `CM` (and `MAF`):
annotation-provided `CM`/`MAF` are ignored, and `ParquetR2RefPanel.prepare_chromosome` supplies the sidecar values for every retained SNP. (Annotation files
no longer carry meaningful `CM`/`MAF` — `CM` is a NaN placeholder and `MAF` is
not carried — so in practice the sidecar is the only source.)

## MAF and counts

The paired sidecar must supply usable `MAF`; the loader rejects missing MAF rather than substituting annotation frequencies.

- `--maf-min` filtering works normally.
- Common-SNP count vectors (and the common-universe overlap matrix) use
  inclusive `MAF >= common_maf_min` from the sidecar.
- All-SNP count vectors are computed from retained annotation rows.

## Memory

Build and read workflows have different working sets. Streaming avoids materializing all pairs, but does not make the builder's genotype window independent of the requested LD window.

### Build side (`build-r2-panel`)

Unrestricted builds stream genotypes from disk. With `--ref-panel-snps-file`, the selected genotype payload is read into RAM. Individual filtering is applied during the read. Pair computation retains the current batch and window-spanning carry-over columns, plus pending pair batches; memory therefore depends on retained individuals, selected SNPs, window density, and `--snp-batch-size`. Wider windows also increase computation and output size. A positive `--min-r2` reduces emitted pairs, and `ldscore` rejects such a panel (`query-r2` still reads it). See `ReferencePanelBuilder._build_chromosome` in `src/ldsc/ref_panel_builder.py` and `yield_pairwise_r2_rows` in `src/ldsc/_kernel/ref_panel_builder.py`.

### Read side (`ldscore`)

The LD-score reader streams every stored pair once per execution batch and projects bounded pair chunks through float64 sparse products.
The main score array covers selected output rows and active annotation columns, including shared baseline and weight columns.
Other allocations include a bounded decoded Arrow batch, pair buffers and CSR construction, annotation/product tiles, sidecar mappings, and output conversion.
The read-time window filters contributions before sparse construction without allocating a dense window matrix; it can affect both work and the retained sparse operator.
Each chromosome worker needs its own working set.
See the [read-side memory model](parquet-r2-format-and-read-pipeline.md#36-read-side-memory-streaming) and the [package resource design](annotation-memory-design.md).

## Practical contract

Keep `chr{chrom}_meta.tsv.gz` alongside `chr{chrom}_r2.parquet` for every
chromosome and emitted build. The parquet is bound to its sidecar by a SHA-256
identity hash recorded in the parquet metadata; the reader hard-fails if the
sidecar is wrong, reordered, or edited. See
`docs/current/parquet-r2-format-and-read-pipeline.md` §4 for the binding
specification.
