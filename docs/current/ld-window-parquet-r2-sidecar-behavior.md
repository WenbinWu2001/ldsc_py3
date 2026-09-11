# LD Window Behavior for Parquet R2 Panels

Last updated on: 2026-09-11

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
| `--ld-wind-snps` | sidecar row order | Normal operation. |
| `--ld-wind-kb` | sidecar `POS` | Normal operation. |
| `--ld-wind-cm` | sidecar `CM` | Normal when `CM` is non-missing. Fails when any retained SNP has missing `CM`. |

The shared LD-window resolver requires exactly one of `--ld-wind-snps`,
`--ld-wind-kb`, or `--ld-wind-cm`.

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
`CM=NA` for every SNP. Use `--ld-wind-snps` or `--ld-wind-kb` in that case, or
regenerate with a genetic-map input.

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

Unrestricted builds stream genotypes from disk. With `--ref-panel-snps-file`, the selected genotype payload is read into RAM. Individual filtering is applied during the read. Pair computation retains the current batch and window-spanning carry-over columns, plus pending pair batches; memory therefore depends on retained individuals, selected SNPs, window density, and `--snp-batch-size`. Wider windows also increase computation and output size. A positive `--min-r2` reduces emitted pairs at the cost of scientific completeness. See `ReferencePanelBuilder._build_chromosome` in `src/ldsc/ref_panel_builder.py` and `yield_pairwise_r2_rows` in `src/ldsc/_kernel/ref_panel_builder.py`.

### Read side (`ldscore`)

The LD-score reader streams every stored pair once and accumulates LD scores using bounded chunks of sparse matrix multiplication. Memory includes the float64 accumulator (`m · n_a · 8` bytes), a bounded pair/CSR chunk, one decoded row group, and workflow overhead. The read-time LD window filters contributions without allocating a dense window matrix. Cross-chromosome workers each need their own working set. See the [read-side memory model](parquet-r2-format-and-read-pipeline.md#36-read-side-memory-streaming).

## Practical contract

Keep `chr{chrom}_meta.tsv.gz` alongside `chr{chrom}_r2.parquet` for every
chromosome and emitted build. The parquet is bound to its sidecar by a SHA-256
identity hash recorded in the parquet metadata; the reader hard-fails if the
sidecar is wrong, reordered, or edited. See
`docs/current/parquet-r2-format-and-read-pipeline.md` §4 for the binding
specification.
