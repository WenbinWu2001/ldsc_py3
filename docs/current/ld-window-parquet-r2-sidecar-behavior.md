# LD Window Behavior for Parquet R2 Panels

This document records how `ldsc ldscore` interprets LD-window flags when the
reference panel backend is package-style parquet R2, especially when the
per-chromosome metadata sidecar is absent.

## Normal Parquet R2 Inputs

A complete package-built parquet R2 chromosome has two artifacts:

```text
chr{chrom}_r2.parquet
chr{chrom}_meta.tsv.gz
```

The R2 parquet stores pairwise LD rows:

```text
CHR  POS_1  POS_2  R2  SNP_1  SNP_2
```

The metadata sidecar stores the retained SNP universe:

```text
CHR  POS  SNP  CM  MAF
```

At LD-score runtime, pairwise R2 values come from `chr{chrom}_r2.parquet`.
Window coordinates, retained SNP order, identifiers, `CM`, and `MAF` come from
the runtime metadata table, normally loaded from `chr{chrom}_meta.tsv.gz`.

## Missing Sidecar Fallback

If the metadata sidecar is missing and
`GlobalConfig.fail_on_missing_metadata` is false, `ParquetR2RefPanel` falls
back to scanning R2 pair endpoints. It reads the pair parquet columns and
synthesizes minimal SNP metadata:

```text
CHR  POS  SNP  CM
```

In this fallback:

- `CHR`, `POS`, and `SNP` come from the union of `SNP_1`/`SNP_2` and
  `POS_1`/`POS_2` endpoints.
- duplicate endpoint rows are removed.
- rows are sorted by chromosome, position, and SNP.
- `CM` is present but missing for every row.
- `MAF` is unavailable.
- SNPs with no emitted R2 pair cannot be recovered and are absent from the
  runtime reference universe.

If `GlobalConfig.fail_on_missing_metadata` is true, missing sidecars raise
immediately instead of using this fallback.

## Window Flag Behavior Without Sidecar

| Flag | Coordinate source without sidecar | Behavior |
|---|---|---|
| `--ld-wind-snps` | synthesized row order | Can run, but only over SNPs recovered from R2 pair endpoints. |
| `--ld-wind-kb` | synthesized `POS` | Can run, but only over SNPs recovered from R2 pair endpoints. |
| `--ld-wind-cm` | synthesized `CM`, which is all missing | Fails unless complete `CM` is supplied by annotation metadata. |

The shared LD-window resolver requires exactly one of `--ld-wind-snps`,
`--ld-wind-kb`, or `--ld-wind-cm`.

For `--ld-wind-cm`, the check is strict: every retained SNP must have a
non-missing `CM`. If any retained row has missing `CM`, chromosome computation
raises:

```text
--ld-wind-cm requires non-missing CM values for all retained SNPs.
```

## Annotation Interaction

When no baseline annotations are supplied, `ldscore` synthesizes an all-ones
`base` annotation from `ref_panel.load_metadata()`. If the sidecar is missing,
that synthetic annotation is built from the fallback endpoint metadata, so the
runtime SNP universe is limited to SNPs present in at least one R2 pair.

When baseline annotations are supplied, the reference-panel metadata is used to
intersect the annotation rows with the retained reference-panel universe. The
annotation metadata then enters chromosome computation. For LD-score
calculation, annotation-provided `CM` is the first source; sidecar metadata only
fills missing `CM` values later.

This means `--ld-wind-cm` can still work without a sidecar only when the
annotation inputs provide complete non-missing `CM` values for all retained
rows. Without annotation-side `CM`, the fallback metadata has missing `CM`, and
the cM window check raises the same missing-`CM` error.

## MAF and Counts

The fallback endpoint metadata has no `MAF`.

Consequences:

- `--maf-min` cannot be applied; the code logs a warning and keeps rows.
- common-SNP count vectors are unavailable because `MAF >= common_maf_min`
  cannot be evaluated.
- all-SNP count vectors are still computed from the retained annotation rows.

## Practical Contract

Package-built parquet R2 panels should keep `chr{chrom}_meta.tsv.gz` alongside
`chr{chrom}_r2.parquet` for every chromosome and emitted build.

The missing-sidecar fallback is a compatibility path for limited external or
legacy inputs. It is acceptable for some SNP- and kb-window runs, but it is not
equivalent to a complete reference panel because it cannot recover singleton
SNPs, SNPs with no emitted pairs, `CM`, or `MAF`.
