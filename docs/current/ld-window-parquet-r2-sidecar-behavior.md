# LD Window Behavior for Parquet R2 Panels

This document records how `ldsc ldscore` interprets LD-window flags when the
reference panel backend is the canonical index-format parquet R2 pair.

## Paired artifacts

A complete package-built parquet R2 chromosome has exactly two artifacts, both
mandatory:

```text
chr{chrom}_r2.parquet   — 4-column index format: IDX_1, IDX_2, R2, SIGN
chr{chrom}_meta.tsv.gz  — per-SNP metadata sidecar (CHR, POS, SNP, A1, A2, CM, MAF)
```

The R2 parquet stores pairwise LD as sidecar-row indices — it carries no SNP
identity of its own and is meaningless without the exact matching sidecar. The
metadata sidecar defines the full panel SNP universe and the index space.

The sidecar format:

```text
# ldsc:schema_version=1
# ldsc:artifact_type=ref_panel_metadata
# ldsc:snp_identifier=chr_pos_allele_aware
# ldsc:genome_build=hg38
CHR  POS  SNP  CM  MAF  A1  A2
```

At LD-score runtime, pairwise R2 values come from `chr{chrom}_r2.parquet`.
Window coordinates, retained SNP order, identifiers, `CM`, and `MAF` come from
the runtime metadata table loaded from `chr{chrom}_meta.tsv.gz`.

## Sidecar is mandatory

The metadata sidecar is a hard requirement. If it is absent, `ParquetR2RefPanel`
raises immediately with an actionable message — there is no synthesize-from-endpoints
fallback. Regenerate the panel with `ldsc build-ref-panel` to produce a paired
parquet + sidecar.

## LD-window flag behavior

| Flag | Coordinate source | Behavior |
|---|---|---|
| `--ld-wind-snps` | sidecar row order | Normal operation. |
| `--ld-wind-kb` | sidecar `POS` | Normal operation. |
| `--ld-wind-cm` | sidecar `CM` | Normal when `CM` is non-missing. Fails when any retained SNP has missing `CM`. |

The shared LD-window resolver requires exactly one of `--ld-wind-snps`,
`--ld-wind-kb`, or `--ld-wind-cm`.

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
LD-score calculation, annotation-provided `CM` is the first source; sidecar `CM`
only fills missing values from the annotation side.

## MAF and counts

`MAF` is available when present in the sidecar (package builds always write it).

- `--maf-min` filtering works normally.
- Common-SNP count vectors use `MAF >= common_maf_min` from the sidecar.
- All-SNP count vectors are computed from retained annotation rows.

## Memory

Build and read have **opposite** memory profiles with respect to the LD window,
so they are described separately.

### Build side (`build-ref-panel`)

The builder never loads the whole chromosome's genotypes into RAM. When a SNP
restriction is supplied (`--ref-panel-snps-file` / `--use-hm3-snps`) it reads only
the kept SNP blocks from the `.bed`; the default unrestricted build streams a
sliding window directly from disk. Individual filtering (`--keep-indivs-file`) is
fused into the per-SNP read, so the raw and filtered bitarrays never coexist.

As a result, **build** peak RSS is governed by that bounded genotype read plus the
workflow/import floor — **not** by `--snp-batch-size`, `--ld-wind-*`, `--min-r2`,
or `--maf-min`. Those control speed and output size; the window/min-r2 options
also change the pair count and pending-pair working set, but they are not build
peak-RSS levers. R2 pairs are emitted as columnar batches; the on-disk parquet
format is unchanged.

### Read side (`ldscore`)

Reading a parquet R2 panel for LD-score computation is the reverse: peak RSS is
dominated by the reader's **decoded row-group cache**, not a bounded genotype
read. That cache is auto-sized per chromosome, and its size scales with **LD
window width × local SNP density**. So on the read side `--ld-wind-*` *is* a
peak-RSS lever, and dense, wide-window chromosomes (notably the chr6 MHC) set the
high-water mark; under cross-chromosome parallelism each worker holds its own
chromosome's cache, so aggregate RSS scales with the worker count. See
`docs/current/parquet-r2-format-and-read-pipeline.md` §3.6 for the cache-sizing
formula, worked chr6/chr22 numbers, and per-worker budgeting.

## Practical contract

Keep `chr{chrom}_meta.tsv.gz` alongside `chr{chrom}_r2.parquet` for every
chromosome and emitted build. The parquet is bound to its sidecar by a SHA-256
identity hash recorded in the parquet metadata; the reader hard-fails if the
sidecar is wrong, reordered, or edited. See
`docs/current/parquet-r2-format-and-read-pipeline.md` §4 for the binding
specification.
