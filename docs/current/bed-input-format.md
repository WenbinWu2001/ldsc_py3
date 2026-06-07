# BED Input Format

This document defines the package contract for text BED interval inputs. It
applies to user-supplied query annotation BEDs and region-exclusion BEDs:

- `--query-annot-bed-sources`
- `AnnotationBuildConfig(query_annot_bed_sources=...)`
- `AnnotationBuilder.project_bed_annotations(...)`
- `run_bed_to_annot(...)`
- `--exclude-regions-bed`
- `RefPanelConfig(exclude_regions_bed=...)`
- `ReferencePanelBuildConfig(exclude_regions_bed=...)`

This contract does not apply to PLINK binary `.bed` genotype files. PLINK
inputs are accepted through `plink_prefix` / `--plink-prefix` and must be part
of a complete `.bed/.bim/.fam` trio.

## Standard BED

Text BED inputs must use BED3 or wider format:

```text
chrom  start  end  [optional extra columns...]
```

The first three fields are required:

| Field | Meaning | Requirement |
| --- | --- | --- |
| `chrom` | Chromosome label | Accepted with or without `chr` prefix, for example `1`, `chr1`, `X`, or `chrX`. |
| `start` | Interval start | Integer, 0-based, inclusive. |
| `end` | Interval end | Integer, 0-based, exclusive. Must be greater than `start`. |

All intervals are interpreted as 0-based half-open intervals: `[start, end)`.
For a 1-based SNP position `p`, the corresponding BED coordinate is `p - 1`;
the SNP overlaps an interval when `start <= p - 1 < end`.

Extra columns are allowed. LDSC uses only interval membership and ignores extra
BED columns during annotation projection and region exclusion.

## Skipped Lines

BED readers should ignore non-interval lines before validating data rows:

- blank lines
- comment lines whose first non-whitespace character is `#`
- UCSC `track ...` lines
- UCSC `browser ...` lines
- one optional leading header row whose first three fields are, ignoring case:
  - `chrom start end`
  - `chrom chromStart chromEnd`

The optional `chrom start end` or `chrom chromStart chromEnd` header is skipped
only when it appears before the first data row. A header-like line after data
rows is treated as malformed input rather than silently ignored.

Examples:

```bed
# enhancer intervals
chrom	start	end	name
chr1	1000	1250	enhancer_A
chr1	5000	5400	enhancer_B
```

```bed
browser position chr6:25000000-35000000
track name="mhc_exclusion" description="MHC exclusion interval"
chrom	chromStart	chromEnd
chr6	25000000	35000000
```

```bed
# centromere intervals
browser hide all
track name="centromeres"
1	121535434	124535434
2	92326171	95326171
```

## Delimiters

BED is conventionally tab-delimited. The package should prefer tab-delimited
examples and documentation, but accept any whitespace between fields in user
inputs. CSV input is not part of the BED contract.

## Validation

After skipped lines are removed, every data row must satisfy the standard BED
requirements:

- at least three fields
- integer `start` and `end`
- non-negative `start`
- `start < end`

Files that violate these requirements should raise a user-facing input error
that identifies the file and line where parsing failed.

## Workflow-Specific Handling

The shared syntax contract is the same for all text BED interval inputs, but
workflow-specific transforms still happen after parsing:

- query annotation BEDs may be expanded by `bed_padding_bp` /
  `--bed-padding-bp`; starts are clipped at zero after padding
- region-exclusion BEDs are coalesced by chromosome before masking reference
  panel SNPs
- preset region exclusions (`--exclude-regions mhc,centromeres`) load packaged
  BED files under `src/ldsc/data/regions/` and should obey the same syntax
  rules

## Compression

Both query annotation BEDs and region-exclusion BEDs should accept plain text
`.bed` files and gzip-compressed `.bed.gz` files. Compression support should be
consistent across the two text BED interval input paths.
