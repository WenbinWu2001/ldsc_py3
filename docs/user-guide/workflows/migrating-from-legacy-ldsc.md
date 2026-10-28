# Migrating from legacy LDSC

Last updated on: 2026-10-28

LDSC3 implements the same method as legacy LDSC, LD score regression (LDSC), behind one `ldsc` command with subcommands and new file formats.
Build new commands from `ldsc COMMAND --help` and the workflow pages rather than translating old ones line by line.
The [legacy LDSC flag map](../../reference/legacy-flag-map.md) lists the replacement for each changed legacy LDSC flag, the removed flags, and the behaviors that changed on purpose, such as exact PLINK LD windows.
This page summarizes what changed, which legacy files you can reuse, and how to convert a legacy LD-score suite.

## What changed

| Topic | LDSC3 behavior |
| --- | --- |
| SNP identity | The default SNP identifier mode is `chr_pos_allele_aware`; `chr_pos`, `rsid` and `rsid_allele_aware` are explicit alternatives where supported. Use the same SNP identifier mode across inputs. |
| Exceptions to the default | `build-gene-ldscore-index` requires an explicit SNP identifier mode (any of the four), hg19 or hg38, and exactly one cM, kb or SNP-count LD window; indexed LD scoring inherits the index settings. `convert-legacy-ldscores` defaults to `rsid` and accepts only the base modes. |
| Genome build | Build labels in the inputs must agree with an explicitly given build after alias normalization; unlabeled inputs use that build without inference or conversion. With `auto`, LDSC3 reconciles declarations and inference across inputs and chromosomes. See [SNP identity and genome build](../concepts/snp-identity-and-genome-build.md). |
| Baseline annotations | LD scores for query annotations require explicit baseline annotations. Only ordinary unpartitioned LD scores can use the synthetic all-ones `base` annotation. Indexed scoring uses the baseline stored in the gene LD-score index. |
| Query annotations | Supply prebuilt annotation files, BED files or gene lists as mutually exclusive routes. `ldscore` projects BED files and gene lists directly; `annotate` optionally writes reusable query annotation files. |
| LD scores | Regression reads one LD-score directory with `metadata.json`, `ldscore.baseline.parquet` and any query and overlap Parquet files, through `--ldscore-dir`. It never recomputes LD scores. |
| Regression SNPs | New LD scores default to the packaged HM3 map as regression SNPs, followed by MHC and centromere exclusion. This does not restrict the reference SNPs or the annotation counts. Munging separately restricts summary statistics to the HM3 map by default. |
| Per-query results | `partitioned-h2` always writes one complete baseline-plus-query result per query; the removed `--write-per-query-results` flag is rejected. |

## Reuse legacy files

| Legacy file | Use in LDSC3 |
| --- | --- |
| Munged `.sumstats` or `.sumstats.gz` with `SNP`, `A1`, `A2`, `Z` and `N` | Pass directly to `h2`, `partitioned-h2` or `rg` with `--sumstats-file` or `--sumstats-sources`; munging again is optional. LDSC3 matches these files by rsID to the LD-score directory and checks allele orientation; see [legacy summary statistics](../../reference/formats/legacy-sumstats.md). |
| Raw GWAS tables | Munge with [`munge-sumstats`](munge-sumstats.md). |
| Complete one-column LD-score suite with regression weights | Convert once with `convert-legacy-ldscores` (below); use the result with [`h2`](heritability.md) and [`rg`](genetic-correlation.md). |
| Complete baseline-only partitioned suite | Convert once with `convert-legacy-ldscores`; use the result with [`partitioned-h2`](partitioned-heritability.md) for the joint baseline model. |
| Cell-type or other query LD scores, thin annotations, `--ref-ld-chr-cts` manifests | Not convertible. Compute query LD scores with LDSC3 from baseline annotations; see [partitioned and cell-type analyses](#partitioned-and-cell-type-analyses). |

Compatibility runs one way: LDSC3 reads these legacy files, but its own outputs are not promised to work in legacy LDSC.
The `BP` header in the `.annot.gz` files that LDSC3 writes is a narrow convenience, and LDSC3 accepts either `BP` or `POS` when reading annotation files.

## Convert a legacy LD-score suite

`ldsc convert-legacy-ldscores` imports a complete legacy reference and regression-weight suite into an LD-score directory, reusing the supplied LD scores without recomputing them.
Regression commands then read the converted directory through `--ldscore-dir`.
The source files are never changed.

The converter detects the suite type from the files, with no type flag:

- **Unpartitioned**: reference LD scores with one value column and no annotation files; the converted directory has one baseline column named `base`, for `h2` and `rg`.
- **Baseline partitioned**: reference LD scores with several columns plus complete full annotation files; the converted directory keeps the baseline columns, for the baseline-only form of `partitioned-h2`.

Query and cell-type LD scores, thin annotations and legacy summary statistics are not inputs to this command.

### Arrange the input files

Place the files directly inside the directory passed to each flag; subdirectories are not searched.
Each reference, weight and frequency family must contain all chromosomes 1 to 22 under one constant prefix, with unpadded chromosome numbers placed immediately before `.l2.ldscore`, `.annot`, `.l2.M`, `.l2.M_5_50` or `.frq`.
The prefix may be empty; include any separator in it, so `baseline.` gives `baseline.1.l2.ldscore.gz`.
In the table and diagrams, `(.gz)` means plain text or gzip, `{1..22}` means 22 separate files, and `<prefix>` and `<chrom>` are not literal text.

| Input | Directory flag | Filename pattern | Examples |
| --- | --- | --- | --- |
| Reference LD scores | `--legacy-reference-dir` | `<prefix><chrom>.l2.ldscore(.gz)` | `1.l2.ldscore.gz`, `baseline.1.l2.ldscore` |
| Regression-weight LD scores | `--legacy-weight-dir` | `<prefix><chrom>.l2.ldscore(.gz)` | `1.l2.ldscore.gz`, `weights.1.l2.ldscore.gz` |
| Common-SNP counts (required) | `--legacy-reference-dir` | `<reference-prefix><chrom>.l2.M_5_50` | `1.l2.M_5_50`, `baseline.1.l2.M_5_50` |
| All-SNP counts (optional) | `--legacy-reference-dir` | `<reference-prefix><chrom>.l2.M` | `1.l2.M`, `baseline.1.l2.M` |
| Full baseline annotations | `--legacy-reference-dir` | `<reference-prefix><chrom>.annot(.gz)` | `baseline.1.annot.gz` |
| Frequencies (baseline suites) | `--legacy-frequency-dir` | `<prefix><chrom>.frq(.gz)` | `1000G.EUR.QC.1.frq.gz`, `1.frq` |

Annotation and count files must use exactly the reference LD-score prefix; weight and frequency prefixes can differ.
Count files are plain text only.
Supply one complete family per role: if a directory holds two complete releases, move each release to its own directory, because there is no prefix-selection flag.
A plain and a gzip copy of the same file are accepted only when their decompressed contents agree; the gzip copy is then used.
Matching names do not make two resources compatible; take every family from the same release.

### Convert an unpartitioned suite

```text
/path/to/legacy/
├── reference/
│   ├── {1..22}.l2.ldscore.gz
│   ├── {1..22}.l2.M_5_50
│   └── {1..22}.l2.M                 # optional
└── weights/
    └── weights.{1..22}.l2.ldscore.gz
```

```bash
ldsc convert-legacy-ldscores \
  --legacy-reference-dir /path/to/legacy/reference \
  --legacy-weight-dir /path/to/legacy/weights \
  --output-dir /path/to/converted/unpartitioned
```

Each reference and weight table needs `CHR SNP BP` followed by exactly one numeric LD-score column, usually `L2`; allele columns are not needed, and the weight directory needs no count files:

```text
CHR SNP BP L2
1 rs1 100 1.01
```

The `.l2.M` and `.l2.M_5_50` files hold headerless counts, one value for the single reference column.
They count the reference SNPs, which can be more numerous than the rows of the LD-score table, so do not create them by counting table rows.
If any chromosome lacks `.l2.M`, conversion warns and records the all-SNP count as unavailable: regression with the default common-SNP counts still works, but a later `--count-kind all` fails.
If one suite supplies both the reference and the weight LD scores, pass the same directory to `--legacy-reference-dir` and `--legacy-weight-dir` instead of copying its files.

### Convert a baseline-partitioned suite

```text
/path/to/legacy/
├── baseline/
│   ├── baseline.{1..22}.l2.ldscore.gz
│   ├── baseline.{1..22}.annot.gz
│   ├── baseline.{1..22}.l2.M_5_50
│   └── baseline.{1..22}.l2.M        # optional; reconstructed if absent
├── weights/
│   └── weights.{1..22}.l2.ldscore.gz
└── frequencies/
    └── 1000G.EUR.QC.{1..22}.frq.gz
```

```bash
ldsc convert-legacy-ldscores \
  --legacy-reference-dir /path/to/legacy/baseline \
  --legacy-weight-dir /path/to/legacy/weights \
  --legacy-frequency-dir /path/to/legacy/frequencies \
  --output-dir /path/to/converted/baseline
```

The reference LD-score tables can hold several columns.
The full annotations need `CHR BP SNP CM` followed by numeric annotation columns, and every LD-score column must map to exactly one annotation column: the names match exactly, or the LD-score name adds one `L2` suffix.
Output columns keep the LD-score spelling.
For example, these headers correspond:

```text
LD-score table:  CHR SNP BP baseL2 codingL2
Annotation:      CHR BP SNP CM base coding
Frequency:       SNP MAF
```

Frequency files need `SNP` and exactly one usable `MAF` or `FRQ` column, and every annotation SNP needs one usable frequency match by rsID.
The converter rebuilds the annotation counts and overlaps from the annotations and frequencies, using legacy LDSC's strict common-SNP rule (`MAF > 0.05`, or `0.05 < FRQ < 0.95`), not the inclusive threshold of new LDSC3 LD scores.
Supplied count files are checked against this reconstruction: a missing `.l2.M` is reconstructed, a missing `.l2.M_5_50` is an error, and a conflicting count is an error.
The count policies are detailed in [legacy LD-score conversion](../../reference/formats/legacy-ldscores.md#counts).

The reference and weight tables need not contain the same SNPs: conversion keeps their rsID intersection, and the coordinates come from the reference LD-score table.
Duplicate rsIDs, non-finite values, incompatible annotation or count layouts, and an empty intersection are rejected.
The result is a reusable baseline model; it is not a way to import a complete historical baseline-plus-cell-type analysis.

### Flags

| Flag | Value and purpose | If omitted |
| --- | --- | --- |
| `--legacy-reference-dir` | Directory with the reference LD scores and counts, plus full annotations for a baseline suite | Required |
| `--legacy-weight-dir` | Directory with one-column regression-weight LD scores; may equal the reference directory for an unpartitioned suite | Required |
| `--legacy-frequency-dir` | Directory with a chromosome 1 to 22 `.frq(.gz)` family | Required for a baseline suite; not read for an unpartitioned suite |
| `--output-dir` | Destination for the converted LD-score directory and its diagnostics | Required |
| `--snp-identifier` | `rsid` or `chr_pos` | `rsid` |
| `--genome-build` | `auto`, `hg19` or `hg38`: build of the reference coordinates | `auto` |
| `--log-level` | `DEBUG`, `INFO`, `WARNING` or `ERROR` | `INFO` |
| `--overwrite` | Replace this command's earlier outputs | Existing outputs cause an error |

Pass directories, not files, prefixes, `@` placeholders or globs.
There is no output-format, suite-type, prefix, chromosome-subset or common-MAF flag.

### Choose the SNP identifier mode and genome build

The defaults, `--snp-identifier rsid --genome-build auto`, suit rsID-based conversion.
If the build cannot be inferred, rsID conversion warns and continues, because rsID identity does not depend on the build; the output's build is then recorded as empty.

For coordinate identity, use `--snp-identifier chr_pos` with `--genome-build auto` when inference should succeed, or with the known `hg19` or `hg38`.
Coordinate identity requires a resolved build and unique coordinates; an explicit build that conflicts with clear evidence, or with build labels in the reference, weight or annotation files, is rejected.
Conversion never lifts over coordinates or infers alleles, so the converted LD scores are allele-unaware in both modes.

Munge new summary statistics with a SNP identifier mode of the same family: `rsid` against an `rsid` conversion, `chr_pos` against a `chr_pos` conversion.
Allele-aware munged summary statistics need `--allow-identity-downgrade` at regression, and rsID-family and coordinate-family modes never mix; see [regression provenance and downgrade](../../reference/snp-identity-and-genome-build.md#regression-provenance-hierarchy-and-downgrade).
Legacy `.sumstats(.gz)` files are matched by rsID automatically.

### Outputs

```text
converted/
├── metadata.json
├── ldscore.baseline.parquet
├── ldscore.overlap.parquet          # baseline-partitioned conversion only
└── diagnostics/
    ├── convert-legacy-ldscores.log
    └── conversion_issues.tsv.gz
```

`ldscore.baseline.parquet` holds the retained regression SNPs, the reference LD-score columns and the regression weights.
`metadata.json` records the SNP identity, the counts, the selected source files and any ignored alternatives.
There is no query file.
Pass the whole directory to the regression command's `--ldscore-dir`.

### When conversion fails

Read the error message and `diagnostics/convert-legacy-ldscores.log`; naming errors state the accepted pattern and a repair, and the log keeps this text even at `--log-level ERROR`.
`diagnostics/conversion_issues.tsv.gz` records the error and any per-file or per-row issues.
A failed directory holds only diagnostics and cannot be used for regression; after fixing the inputs, rerun with a new output directory or add `--overwrite`.
Errors that occur before logging starts, such as a missing input directory or an output collision, appear only in the terminal.

An unused `.w` file next to a complete recognized family does not block conversion: it is ignored with a warning, recorded as `discarded_unsupported_weight_filename`, and listed in `legacy_ldscore_conversion.ignored_files`.

| Problem | Example | Fix |
| --- | --- | --- |
| `.w` between the chromosome and `.l2.ldscore` | Only `1.w.l2.ldscore.gz` to `22.w.l2.ldscore.gz` | Copy them as `weights.1.l2.ldscore.gz` to `weights.22.l2.ldscore.gz` into a new directory and pass it to `--legacy-weight-dir` |
| Padded chromosome numbers | `weights.01.l2.ldscore.gz` to `weights.22.l2.ldscore.gz` | Rename to unpadded `weights.1.l2.ldscore.gz` to `weights.22.l2.ldscore.gz` |
| Wrong directory level | Files are in `/path/to/legacy/weights/`, but the flag points to `/path/to/legacy/` | Pass the directory that directly contains the files |
| Missing chromosome | A family has chromosomes 1 to 21 only | Restore chromosome 22 from the same release; single-chromosome conversion is not supported |
| Several complete releases | `releaseA.{1..22}.l2.ldscore.gz` and `releaseB.{1..22}.l2.ldscore.gz` in one directory | Put each release in its own directory and select one |
| Wrong frequency suffix | `1000G.EUR.QC.1.freq.gz` | If the contents are a valid frequency table, rename copies to `1000G.EUR.QC.1.frq.gz` consistently across chromosomes |
| Annotation prefix differs from the reference | Reference `baseline.1.l2.ldscore.gz`, annotation `other.1.annot.gz` | For annotations from the same suite, use `baseline.1.annot.gz`; apply the same prefix to the counts |
| Wrong common-count suffix | `baseline.1.M_5_50` or `baseline.1.l2.M_5_50.gz` | Supply the plain-text count as `baseline.1.l2.M_5_50`; decompress a gzip file before removing `.gz` |
| Missing baseline frequencies | Baseline suite without `--legacy-frequency-dir` | Supply the matching complete `.frq(.gz)` directory |
| Unsupported table structure | Thin annotations; several weight columns; several reference columns without full annotations | Obtain a supported, coherent suite; renaming cannot repair the contents |

A `.w.l2.ldscore.gz` file may hold valid regression weights; only its name is not recognized.
To rename without touching the source suite, copy into a new directory; for example, this Bash loop copies 22 gzip files without changing their bytes:

```bash
mkdir -p /path/to/organized/weights
for chrom in {1..22}; do
  cp -n "/path/to/legacy/weights/${chrom}.w.l2.ldscore.gz" \
    "/path/to/organized/weights/weights.${chrom}.l2.ldscore.gz"
done
```

Then pass `--legacy-weight-dir /path/to/organized/weights`.
The same applies to plain files without `.gz`; changing an extension does not compress or decompress a file.
When reorganizing a reference suite, rename its annotations and counts consistently and keep the original release.
See [troubleshooting](../troubleshooting.md#convert-legacy-ldscores) for more repairs.

## Partitioned and cell-type analyses

- Legacy LDSC computes new LD scores only for a query annotation and reuses precomputed baseline LD scores, annotations, frequency files and regression weights.
  LDSC3 direct mode computes LD scores for the baseline and query annotations in the same run and derives the counts and weights from the reference panel; indexed mode takes the fixed baseline from a gene LD-score index and assembles the query LD scores from precomputed parts.
- Precomputed baseline LD scores do not replace the baseline annotations: new query LD scores need the baseline annotation files, and a converted baseline directory cannot be extended with queries.
- Legacy LDSC `--h2-cts` with `--ref-ld-chr-cts` becomes `ldscore` with query annotations followed by `partitioned-h2`; old cell-type manifests are not accepted.
  See [partitioned heritability](partitioned-heritability.md) and [cell-type enrichment](cell-type-enrichment.md).
- LDSC3 starts from the packaged HM3 map as regression SNPs and removes MHC and centromere SNPs from them by default.
  Legacy LDSC uses the regression SNPs you supply, such as `w_hm3.print_snps`, and has no built-in MHC or centromere exclusion step.

## Reference

- [legacy LDSC flag map](../../reference/legacy-flag-map.md)
- [Legacy LD-score conversion](../../reference/formats/legacy-ldscores.md)
- [Legacy summary statistics](../../reference/formats/legacy-sumstats.md)
- [SNP identity and genome build](../../reference/snp-identity-and-genome-build.md)
- [Artifact metadata](../../reference/formats/artifact-metadata.md)
- [Troubleshooting](../troubleshooting.md#convert-legacy-ldscores)

For the Python API, see [Python usage](../python-usage.md#migrating-from-legacy-ldsc).
