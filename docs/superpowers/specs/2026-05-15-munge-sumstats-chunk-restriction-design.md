# Munge-Sumstats Chunk Restriction Design

## Purpose

`munge-sumstats` reads raw GWAS summary statistics, applies LDSC-compatible
quality control, and writes curated sumstats artifacts. When a HapMap3
restriction or `--sumstats-snps-file` keep-list is supplied, the munger should
avoid retaining non-kept rows through the full in-memory concatenation step.

This design moves SNP restriction into the chunk-processing stage while keeping
the user-visible workflow and log messages close to the existing behavior.

## Workflow

1. The public workflow resolves input paths, output paths, output format, and
   the active `snp_identifier` / `genome_build` configuration.
2. The kernel reads the raw sumstats header, applies user column hints, and maps
   raw columns to canonical LDSC columns such as `SNP`, `CHR`, `POS`, `P`,
   allele columns, and sample-size columns.
3. For `chr_pos`-family modes with `genome_build=auto`, the source build is
   resolved before chunk parsing. The munger reads only lightweight raw `CHR`
   and `POS` evidence, calls the shared genome-build inference code, and records
   the inferred source build and coordinate basis.
4. If HM3 filtering or `--sumstats-snps-file` is active, the keep-list is loaded
   once before chunk parsing. Base `rsid` mode keeps canonical SNP labels. Base
   `chr_pos` mode uses the resolved source build to choose the relevant position
   column and stores compact coordinate keys. Allele-aware modes use the shared
   restriction identity reader; allele-free restrictions match by base key and
   allele-bearing restrictions match by effective key.
5. The raw sumstats file is streamed in chunks.
6. For each chunk, the munger drops rows with missing required non-coordinate
   fields, coerces configured INFO-list columns, renames columns to their
   canonical names, and applies the existing INFO, frequency, p-value, and
   allele filters.
7. For `chr_pos`-family modes, each chunk also normalizes valid `CHR` and `POS`
   values and drops rows whose coordinates are missing or invalid. These
   coordinate drops are logged separately from the standard missing-value and QC
   counters.
8. If a keep-list is active, the chunk is restricted before it is appended to
   the retained chunk list. Non-kept rows do not survive to the full concat.
9. Only retained chunk rows are concatenated. Post-concat work remains limited
   to steps that need the retained table: shared identity cleanup, sample-size
   filtering, p-value to Z conversion, signed-statistic direction, optional
   source-build liftover, metadata, and output writing.

## Identifier Behavior

`rsid` mode treats `SNP` labels as row identity. Keep-list filtering happens in
chunk order, so retained rows preserve the input order rather than the keep-list
order. Duplicate effective-key clusters are dropped by the shared identity
cleanup after chunk filtering.

`chr_pos` mode treats normalized source-build `CHR` and `POS` as row identity.
The base-mode munger path filters on compact packed coordinate keys rather than
materializing `CHR:POS` strings. If the keep-list contains build-specific
position columns such as `hg19_POS` and `hg38_POS`, the resolved source build
selects the column. Allele-aware modes require usable `A1/A2` in the raw
sumstats artifact path and use effective identity keys when allele-bearing
restriction files are supplied.

## Genome Build and Liftover

Automatic genome-build inference uses the shared package decision rule: compare
input coordinates against the packaged HapMap3 hypotheses, require at least 200
informative matches, and require the best hypothesis to explain at least 99% of
informative matches. If a 0-based source is inferred, the munger converts
positions to canonical 1-based coordinates before keep-list matching and output.

Keep-list filtering is source-build filtering. Optional summary-stat liftover
runs after source-build SNP filtering and updates coordinates only after the
retained source rows have been selected.

## Empty Results

A chunk may produce zero retained rows. That is normal and does not stop the
run. If all chunks together produce zero retained rows after the keep-list is
applied, the run raises the same user-facing keep-list error class of message:
no SNPs remain, with the keep-list path, identifier mode, genome build, usable
row identifiers, and keep-list size included.

If a keep-list file is valid but contains zero usable identifiers, the munger
fails before scanning the raw sumstats file.
