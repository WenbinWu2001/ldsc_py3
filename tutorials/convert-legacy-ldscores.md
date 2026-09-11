# Convert a legacy LD-score suite

Last updated on: 2026-09-11

Use `ldsc convert-ldsc2-ldscores` to import an existing LDSC2 reference/weight LD-score suite into a canonical LDSC3 directory. This reuses the supplied LD scores without recomputing them. Subsequent regression commands consume the converted directory through `--ldscore-dir`.

The converter supports an ordinary one-column unpartitioned suite for `h2` or `rg`, or a complete baseline-only partitioned suite for `partitioned-h2`. It selects the profile from the files; there is no profile flag. Query/cell-type LD-score models and thin annotations are outside this conversion contract. Legacy sumstats use a separate automatic regression input path and are not inputs to this command.

All `/data/...` paths below are illustrative. Replace them with your own directories. In directory diagrams, `{1..22}` means 22 separate chromosome files; braces, `<prefix>`, and `<chrom>` are not literal filename text.

## Accepted filename patterns

Files must be directly inside the directory passed to the corresponding flag. Discovery is not recursive. Each independently discovered reference, weight, or frequency family must contain all chromosomes 1 through 22 under one constant prefix. Use unpadded chromosome numbers (`1`, `2`, ..., `22`). The chromosome number appears immediately before `.l2.ldscore`, `.annot`, `.l2.M`, `.l2.M_5_50`, or `.frq`.

`<prefix>` may be empty. If you use a separator, include it in the prefix: `baseline.` gives `baseline.1.l2.ldscore.gz`. Reference, weight, and frequency prefixes may differ. Associated annotations and count files must use the selected **reference prefix** exactly. `(.gz)` means either plain text or gzip is recognized; count files have no gzip alternative.

| Input | Directory flag | Recognized pattern | Examples |
| --- | --- | --- | --- |
| Reference LD scores | `--legacy-reference-dir` | `<prefix><chrom>.l2.ldscore(.gz)` | `1.l2.ldscore.gz`, `baseline.1.l2.ldscore` |
| Regression-weight LD scores | `--legacy-weight-dir` | `<prefix><chrom>.l2.ldscore(.gz)` | `1.l2.ldscore.gz`, `weights.1.l2.ldscore.gz` |
| Common-SNP counts; required | `--legacy-reference-dir` | `<reference-prefix><chrom>.l2.M_5_50` | `1.l2.M_5_50`, `baseline.1.l2.M_5_50` |
| All-SNP counts; optional under the count policy below | `--legacy-reference-dir` | `<reference-prefix><chrom>.l2.M` | `1.l2.M`, `baseline.1.l2.M` |
| Full baseline annotations | `--legacy-reference-dir` | `<reference-prefix><chrom>.annot(.gz)` | `baseline.1.annot.gz` |
| Frequencies for baseline conversion | `--legacy-frequency-dir` | `<prefix><chrom>.frq(.gz)` | `1000G.EUR.QC.1.frq.gz`, `1.frq` |

Supply one complete family per role. If two releases each have a complete chromosome family in the same directory, place the releases in separate directories and select the intended directory. There is no prefix-selection flag. A plain/gzip duplicate is accepted only when its decompressed bytes agree; the gzip copy is then selected. Matching naming alone does not establish that the biological resources or table contents are compatible.

## Canonical flags and values

| Flag | Value and purpose | Omission behavior |
| --- | --- | --- |
| `--legacy-reference-dir` | Existing directory containing reference LD scores, counts, and baseline annotations when applicable | Required |
| `--legacy-weight-dir` | Existing directory containing one-column regression-weight LD scores | Required |
| `--legacy-frequency-dir` | Existing directory containing a chromosome 1–22 `.frq(.gz)` family | Required for baseline conversion; omit for unpartitioned conversion |
| `--output-dir` | Destination for the canonical LDSC3 artifact and diagnostics | Required |
| `--snp-identifier` | `rsid` for rsID identity; `chr_pos` for coordinate identity | `rsid` |
| `--genome-build` | `auto`, `hg19`, or `hg38` for reference build inference/declaration | `auto` |
| `--log-level` | `DEBUG`, `INFO`, `WARNING`, or `ERROR` for workflow-log detail; each includes more severe messages | `INFO`; completion/failure status is always recorded |
| `--overwrite` | Switch with no value; permits replacing converter-owned destination files | Off; existing owned files cause an error |

Pass directories, not individual files, filename prefixes, `@` placeholders, or globs. The command has no output-format, profile, prefix, chromosome-subset, or common-MAF-threshold flag. Allele-aware identity choices are not accepted.

## Working example: unpartitioned conversion

```text
/data/legacy/
├── reference/
│   ├── {1..22}.l2.ldscore.gz
│   ├── {1..22}.l2.M_5_50
│   └── {1..22}.l2.M                 # optional
└── weights/
    └── weights.{1..22}.l2.ldscore.gz
```

```bash
ldsc convert-ldsc2-ldscores \
  --legacy-reference-dir /data/legacy/reference \
  --legacy-weight-dir /data/legacy/weights \
  --output-dir /data/converted/unpartitioned \
  --snp-identifier rsid \
  --genome-build auto \
  --log-level INFO
```

Each reference and weight table requires `CHR SNP BP` followed by exactly one numeric LD-score column, commonly `L2`. No allele columns are required. The weight directory does not need count files. For example, a reference chromosome table may begin:

```text
CHR SNP BP L2
1 rs1 100 1.01
```

The `.l2.M` and `.l2.M_5_50` files contain headerless counts, one value for this single reference column. They describe the LD-reference SNP universe, which can be larger than the SNP rows in the LD-score table; do not manufacture counts by counting those output rows. If any chromosome lacks `.l2.M`, conversion records the all-SNP count as unavailable and warns. Common-count regression remains usable; requesting unavailable all-SNP counts later fails.

One suitable unpartitioned suite can supply both reference and weight scores. To use `/data/legacy/reference` for both roles, pass that same directory to both `--legacy-reference-dir` and `--legacy-weight-dir`; do not duplicate or rename its files solely to assign the second role.

## Working example: baseline-partitioned conversion

```text
/data/legacy/
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
ldsc convert-ldsc2-ldscores \
  --legacy-reference-dir /data/legacy/baseline \
  --legacy-weight-dir /data/legacy/weights \
  --legacy-frequency-dir /data/legacy/frequencies \
  --output-dir /data/converted/baseline \
  --snp-identifier rsid \
  --genome-build auto \
  --log-level INFO
```

These reference LD-score tables can contain multiple scientific columns. Full annotations require `CHR BP SNP CM` plus numeric annotation columns. Every LD-score column must map unambiguously to one annotation column: either the names match exactly or the LD-score name appends one `L2` suffix. Scientific output column names retain the input LD-score spelling. For example, these headers correspond:

```text
LD-score table:  CHR SNP BP baseL2 codingL2
Annotation:     CHR BP SNP CM base coding
Frequency:      SNP MAF
```

Frequency files must supply `SNP` and exactly one usable `MAF` or `FRQ` column. Every annotation SNP needs one usable frequency match by rsID. Common-SNP counts use the fixed legacy rule `MAF > 0.05`, or `0.05 < FRQ < 0.95`. Baseline count vectors have one value per reference LD-score column in its column order. Supplied counts are validated against the full annotation/frequency reconstruction. Missing `.l2.M` is reconstructed; missing `.l2.M_5_50` is an error. Conflicting counts are errors, not an invitation to replace them with arbitrary values.

The reference and weight tables need not have identical SNP sets: conversion retains their rsID intersection. Duplicate rsIDs, non-finite scientific values, incompatible annotation/count schemas, and an empty final intersection are rejected. Coordinates in the reference LD-score table supply the output coordinates. These are reusable baseline-model inputs, not a conversion route for a historical baseline-plus-cell-type analysis.

## Choosing identity and build

The default `--snp-identifier rsid --genome-build auto` is sufficient for rsID-based conversion. Failure to infer the build can warn and continue because canonical rsID identity does not depend on a genome build; its output build metadata remains null.

For coordinate identity, use `--snp-identifier chr_pos` with `--genome-build auto` when inference is expected to succeed, or an explicit known `hg19`/`hg38` build. Coordinate identity requires a resolved build and unique canonical coordinates. An explicit build that conflicts with decisive inference is rejected. Conversion does not lift over coordinates or infer alleles. Converted LD scores remain allele-unaware in both modes.

## Non-working suites and practical repairs

| Problem | Example | Practice that fixes it |
| --- | --- | --- |
| Weight chromosome separated from `.l2.ldscore` by `.w` | Only `1.w.l2.ldscore.gz` through `22.w.l2.ldscore.gz` | Rename copies to `weights.1.l2.ldscore.gz` through `weights.22.l2.ldscore.gz`; point `--legacy-weight-dir` at that directory |
| Padded chromosome names | `weights.01.l2.ldscore.gz` through `weights.22.l2.ldscore.gz` | Normalize to unpadded `weights.1.l2.ldscore.gz` through `weights.22.l2.ldscore.gz` |
| Wrong directory level | Files are under `/data/legacy/weights/`, but the flag points to `/data/legacy/` | Pass the directory directly containing the files; discovery does not recurse |
| Missing chromosome | Reference or weight family has chromosomes 1–21 only | Restore chromosome 22 from the same resource release; single-chromosome conversion is not supported |
| Multiple complete releases | Both `releaseA.{1..22}.l2.ldscore.gz` and `releaseB.{1..22}.l2.ldscore.gz` in one directory | Organize each release in a separate directory and select one |
| Wrong frequency suffix | `1000G.EUR.QC.1.freq.gz` | If the contents are a valid frequency table, rename the copy to `1000G.EUR.QC.1.frq.gz`, consistently across chromosomes |
| Associated prefix mismatch | Reference `baseline.1.l2.ldscore.gz`, annotation `other.1.annot.gz` | For annotations from that same suite, use `baseline.1.annot.gz`; apply the same reference prefix to counts |
| Wrong common-count suffix | `baseline.1.M_5_50` or `baseline.1.l2.M_5_50.gz` | Supply the matching plain-text count as `baseline.1.l2.M_5_50`; decompress gzip contents before removing `.gz` |
| Missing baseline frequencies | Baseline suite without `--legacy-frequency-dir` | Supply the matching complete `.frq(.gz)` directory |
| Unsupported table structure | Thin annotations; multiple weight value columns; multi-column reference scores without full annotations | Obtain a supported coherent source suite; renaming alone cannot repair these contents |

A `.w.l2.ldscore.gz` file may contain valid regression-weight data: the limitation is filename discovery. To organize copies while preserving an existing source suite, use a new destination directory. For example, this Bash sequence copies 22 gzip files without changing their bytes:

```bash
mkdir -p /data/organized/weights
for chrom in {1..22}; do
  cp -n "/data/legacy/weights/${chrom}.w.l2.ldscore.gz" \
    "/data/organized/weights/weights.${chrom}.l2.ldscore.gz"
done
```

Then use `--legacy-weight-dir /data/organized/weights`. The same naming correction applies to plain files without `.gz`; changing an extension does not compress or decompress the contents. When organizing a reference suite, rename its associated annotations and counts consistently, and retain the correct source release. Standard names describe roles and chromosome grouping; they do not make unrelated source datasets interchangeable.

## Outputs and failed-run diagnostics

On success the destination contains:

```text
converted/
├── metadata.json
├── ldscore.baseline.parquet
├── ldscore.overlap.parquet          # baseline-partitioned conversion only
└── diagnostics/
    ├── convert-ldsc2-ldscores.log
    └── conversion_issues.tsv.gz
```

The baseline table contains retained regression SNPs, reference LD-score columns, and regression weights. Metadata records identity, counts, selected source files, and ignored alternatives. Baseline conversion also supplies the overlap matrix. These profiles produce no query Parquet file. Pass the whole converted directory to the subsequent regression's `--ldscore-dir`.

For validation failures after output preflight, read the exception and `diagnostics/convert-ldsc2-ldscores.log`. Naming errors include the accepted pattern and a repair practice; the shared failure footer records this text even at `--log-level ERROR`. `diagnostics/conversion_issues.tsv.gz` records the conversion error and any available per-file or per-row issues. A diagnostics-only failed directory cannot be used for regression. After correcting the inputs, rerun with a new output directory or add `--overwrite` to replace the previous converter-owned files. Errors before logging starts, such as an invalid directory argument or output collision, are reported to the caller/terminal and may not create a diagnostic log.

An unused `.w` file alongside a complete recognized family does not invalidate that family: the file is ignored with a warning, audited as `discarded_unsupported_weight_filename`, and listed in `legacy_ldsc2_import.ignored_files`.

## References

- [Conversion contract](../docs/current/legacy-ldscore-conversion.md): counts, identity, schemas, discovery, and provenance.
- [Troubleshooting](../docs/troubleshooting.md#convert-ldsc2-ldscores): repair guidance for rejected suites.
- [Concise utility wiki](../docs/wiki/utility-functionalities/convert-ldsc2-ldscores.md): standard names and flags.
- [Implementation](../src/ldsc/legacy_ldscore_converter.py): `build_parser`, `_discover_ldscore_family`, `_discover_associated_family`, `_discover_frequency_family`, and `LegacyLDScoreConverter.convert` define the behavior described here.
