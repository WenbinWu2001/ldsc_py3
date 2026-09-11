# Convert legacy LD-score files

Last updated on: 2026-09-11

`ldsc convert-ldsc2-ldscores` converts a complete legacy reference/weight suite into a directory accepted by LDSC3 regression through `--ldscore-dir`. Use it for ordinary unpartitioned scores or a baseline-only partitioned suite.

## Standard filenames

Place files directly inside the directory passed to the corresponding flag. Each reference, weight, and frequency family needs all chromosomes **1–22** with one constant prefix. `<chrom>` is an unpadded number; `<prefix>` may be empty or include a separator, such as `baseline.`. `(.gz)` means either plain text or gzip. Counts are plain text.

| Role | Standard filename | Example | Directory flag |
| --- | --- | --- | --- |
| Reference LD scores | `<prefix><chrom>.l2.ldscore(.gz)` | `baseline.1.l2.ldscore.gz` | `--legacy-reference-dir` |
| Regression weights | `<prefix><chrom>.l2.ldscore(.gz)` | `weights.1.l2.ldscore.gz` | `--legacy-weight-dir` |
| Required common-SNP counts | `<reference-prefix><chrom>.l2.M_5_50` | `baseline.1.l2.M_5_50` | `--legacy-reference-dir` |
| Optional all-SNP counts | `<reference-prefix><chrom>.l2.M` | `baseline.1.l2.M` | `--legacy-reference-dir` |
| Full baseline annotations | `<reference-prefix><chrom>.annot(.gz)` | `baseline.1.annot.gz` | `--legacy-reference-dir` |
| Baseline frequencies | `<prefix><chrom>.frq(.gz)` | `1000G.EUR.QC.1.frq.gz` | `--legacy-frequency-dir` |

Annotations and counts use exactly the reference LD-score prefix. Weight and frequency prefixes can differ. Keep different releases in separate directories; pass directories rather than individual files, prefixes, or globs.

## Commands

Replace these illustrative `/data/...` paths with your directories. For an ordinary unpartitioned suite, run:

```bash
ldsc convert-ldsc2-ldscores \
  --legacy-reference-dir /data/legacy/reference \
  --legacy-weight-dir /data/legacy/weights \
  --output-dir /data/converted/unpartitioned
```

If one suitable unpartitioned suite supplies both reference and weight scores, pass the same directory to both input flags.

For a baseline-partitioned suite, include the frequency directory:

```bash
ldsc convert-ldsc2-ldscores \
  --legacy-reference-dir /data/legacy/baseline \
  --legacy-weight-dir /data/legacy/weights \
  --legacy-frequency-dir /data/legacy/frequencies \
  --output-dir /data/converted/baseline
```

| Flag | Accepted values and usage |
| --- | --- |
| `--legacy-reference-dir` | Required directory of reference scores and counts; baseline conversion also needs full annotations there |
| `--legacy-weight-dir` | Required directory of one-column weight LD scores |
| `--legacy-frequency-dir` | Required for baseline conversion; omit for unpartitioned conversion |
| `--output-dir` | Required converted-output directory |
| `--snp-identifier` | `rsid` (default) or `chr_pos` |
| `--genome-build` | `auto` (default), `hg19`, or `hg38`; `chr_pos` requires a known or successfully inferred build |
| `--log-level` | `DEBUG`, `INFO` (default), `WARNING`, or `ERROR` for workflow-log detail |
| `--overwrite` | Switch with no value; use to replace existing converter-owned output files |

For a known hg19 coordinate-based output, add `--snp-identifier chr_pos --genome-build hg19`. The default `rsid` output does not require successful automatic build inference. The converter has no profile, prefix, chromosome-subset, or MAF-threshold flag.

See the [detailed user guide](../../../tutorials/convert-legacy-ldscores.md) for complete layouts, file contents, count policies, and troubleshooting examples. The names and flags above follow [`legacy_ldscore_converter.build_parser` and its discovery helpers](../../../src/ldsc/legacy_ldscore_converter.py).
