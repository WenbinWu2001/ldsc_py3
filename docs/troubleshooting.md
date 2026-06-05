# Troubleshooting

This reference explains `ldsc` errors that can **abort a run** and have more than
one likely cause. It is organized by command. Each entry lists the likely causes
(ranked most-probable first), how to confirm each, and how to fix it.

Most errors are self-explanatory from their terminal message alone and are not
repeated here. When a message says `... see docs/troubleshooting.md#<section>`,
that slug is a heading below — jump to it.

## Common

### Common: input path did not resolve to one file

**Raised by:** `path_resolution.resolve_scalar_path()` and group/path-prefix resolvers
· **Exception:** `LDSCInputError`
**Symptom:** `Could not resolve <label> path from token '<token>': matched 0 files` / `matched N files`

**Likely causes & how to check** (most probable first):

| # | Likely cause | How to check |
|---|--------------|--------------|
| 1 | The path is misspelled or relative to a different working directory | `pwd`; then `ls <path>` from the same shell |
| 2 | A glob is too broad for an input that must be one file | `python -c "import glob; print(glob.glob('<pattern>'))"` |
| 3 | A chromosome-suite token is missing the explicit `@` placeholder | Check whether the token looks like `chr@` / `.@.` rather than a bare prefix |
| 4 | The file exists only with a suffix LDSC does not infer | `ls <prefix>*`; pass the full filename including suffix |
| 5 | A PLINK prefix is incomplete | Confirm all three files exist: `<prefix>.bed`, `<prefix>.bim`, `<prefix>.fam` |

**Remedies:**

1. Pass the exact existing file path when the command expects a single file.
2. Narrow broad globs so they match exactly the intended file.
3. For chromosome suites, use an explicit `@` token such as `baseline.@.annot.gz`.

### Common: output artifact already exists

**Raised by:** `path_resolution.ensure_output_paths_available()` and
`path_resolution.preflight_output_artifact_family()` · **Exception:** `FileExistsError`
**Symptom:** `Cannot write <artifact>: existing output artifact already exists at ...`

**Likely causes & how to check** (most probable first):

| # | Likely cause | How to check |
|---|--------------|--------------|
| 1 | The output directory contains results from an earlier run | `ls <output-dir>` |
| 2 | A previous run used a different output-format option and left stale sibling files | Compare existing files against the current command's `--output-format` / workflow mode |
| 3 | The target directory is shared by two concurrent or interrupted runs | Check job logs and file modification times with `ls -l <output-dir>` |
| 4 | The output path points at a file created manually or by another workflow | Inspect the listed path before overwriting it |

**Remedies:**

1. Use a fresh output directory for independent runs.
2. If replacing prior results is intended, pass `--overwrite` on the CLI or
   `overwrite=True` in Python.
3. Do not share one output directory across concurrent runs unless the workflow
   explicitly supports chromosome-sharded output ownership.

### Common: genome build could not be inferred

**Raised by:** `genome_build_inference.resolve_genome_build()` and
`genome_build_inference.resolve_chr_pos_table()` · **Exception:** `LDSCInputError`
**Symptom:** `Could not infer genome build for <context>...`

**Likely causes & how to check** (most probable first):

| # | Likely cause | How to check |
|---|--------------|--------------|
| 1 | The CHR/POS input has too little overlap with HapMap3 reference SNPs | Count matched/reference-like variants; small custom SNP lists often cannot infer a build |
| 2 | The input coordinates are from another build, assembly, or custom reference | Inspect several known SNP coordinates against hg19/hg38 in a genome browser or source metadata |
| 3 | CHR/POS columns were parsed from the wrong columns or wrong delimiter | Print the header and first rows; confirm `CHR` and `POS` values look like chromosome labels and base-pair positions |
| 4 | Coordinates are mixed across builds or coordinate bases | Check whether some rows match hg19 and others match hg38/0-based positions |

**Remedies:**

1. Pass an explicit build: `--genome-build hg19` or `--genome-build hg38`.
2. Fix column mapping or delimiter issues before using `--genome-build auto`.
3. Rebuild the input so all CHR/POS rows use one genome build and one coordinate basis.

### Common: LDSC artifact schema or provenance is incompatible

**Raised by:** `_kernel.snp_identity.validate_identity_artifact_metadata()` and
artifact reload guards · **Exception:** `LDSCInputError`
**Symptom:** `Could not read LDSC <artifact> artifact metadata...`

**Likely causes & how to check** (most probable first):

| # | Likely cause | How to check |
|---|--------------|--------------|
| 1 | The artifact was written by an older LDSC package version | Inspect the artifact metadata sidecar or parquet metadata for `schema_version` |
| 2 | The artifact type does not match the loader | Confirm the file was produced by the command you are now trying to load from |
| 3 | The metadata sidecar was copied without its matching data file, or vice versa | Compare file modification times and paths for the artifact plus sidecar |
| 4 | SNP identity provenance was hand-edited or corrupted | Inspect `snp_identifier`, `genome_build`, and identity metadata fields |

**Remedies:**

1. Regenerate the artifact with the current LDSC package and the same command family.
2. Keep LDSC-generated data files and their sidecars together; do not hand-edit them.
3. If you need a different `snp_identifier` or `genome_build`, regenerate from the upstream input.

## munge-sumstats

### munge-sumstats: could not map a required column

**Raised by:** `_kernel/sumstats_munger` column resolution · `SumstatsTable.validate()`
· **Exception:** `LDSCInputError`
**Symptom:** `munge-sumstats could not map the required column '<field>' from the header...`

**Likely causes & how to check** (most probable first):

| # | Likely cause | How to check |
|---|--------------|--------------|
| 1 | The header uses a name the auto-mapper doesn't recognize | `zcat <file> \| head -1` (or `head -1`); compare against the recognized aliases in `column_inference.py` |
| 2 | The column exists under a synonym you must declare | Re-run with an explicit hint, e.g. `--snp-col MarkerName --a1-col Allele1` |
| 3 | Wrong delimiter, so the whole header parsed as one column | `head -1 <file> \| cat -A` — look for one field with embedded tabs/commas |
| 4 | Wrong `--format`, so expected columns differ | Confirm the format flag matches the file (e.g. drop `--daner-new` for a non-PGC file) |

**Remedies:**

1. Pass explicit column hints for the unmapped field(s): `--<field>-col <name>`.
2. Rename the columns to recognized names, or add the alias to `column_inference.py`
   if broadly useful.
3. Verify the delimiter is whitespace/tab as expected for `.sumstats`/`.txt` inputs.

### munge-sumstats: no SNPs remain after filtering

**Raised by:** `_kernel/sumstats_munger` (post-filter and keep-list paths)
· **Exception:** `LDSCInputError`
**Symptom:** `munge-sumstats removed every SNP...` / `...no SNPs remain after applying --sumstats-snps-file`

**Likely causes & how to check** (most probable first):

| # | Likely cause | How to check |
|---|--------------|--------------|
| 1 | `--sumstats-snps-file` keep-list uses a different SNP-id space (rsID vs chr:pos) | Compare the first IDs of each file; confirm both use `snp_identifier=<mode>` |
| 2 | Keep-list and sumstats are on different genome builds | Compare the `genome_build` in the error line against the keep-list build |
| 3 | INFO / MAF / N thresholds removed every row | Re-run with relaxed `--info-min` / `--maf-min`; inspect the dropped-SNP sidecar |
| 4 | Input is effectively empty after a delimiter/format mis-parse | `zcat <file> \| wc -l` |

**Remedies:**

1. Regenerate the keep-list in the same `snp_identifier` mode and build as the sumstats.
2. Relax quality filters and re-run; review the dropped-SNP audit sidecar to see
   which filter removed the rows.

### munge-sumstats: curated artifact is malformed or outdated

**Raised by:** `sumstats_munger.load_sumstats()` and its metadata helpers
· **Exception:** `LDSCInputError`
**Symptom:** `Cannot load curated sumstats at '<path>': ...` (sidecar missing / old
schema / bad provenance / missing A1-A2 / duplicate identity rows)

**Likely causes & how to check** (most probable first):

| # | Likely cause | How to check |
|---|--------------|--------------|
| 1 | Artifact predates the current schema (no/old `metadata.json` sidecar) | Check for `metadata.json` beside the artifact; inspect its `schema_version` |
| 2 | The `metadata.json` sidecar is missing or unreadable | `python -m json.tool metadata.json` |
| 3 | Identity provenance in the sidecar is invalid/corrupt | Inspect the sidecar's identity fields against the current contract |
| 4 | An allele-aware `snp_identifier` artifact lacks A1/A2 columns | `python -c "import pandas; print(pandas.read_parquet('<f>').columns)"` |
| 5 | Duplicate/invalid SNP-identity rows survived in the artifact | Re-munge from raw input; the loader reports the dropped-row reasons |

**Remedies:**

1. Re-run `ldsc munge-sumstats` from the **raw** GWAS file to regenerate the
   artifact with the current schema.
2. Do not hand-edit curated `.sumstats`/`.parquet` artifacts; treat them as outputs.

### munge-sumstats: liftover dropped all rows

**Raised by:** `_kernel.liftover.apply_sumstats_liftover()`
· **Exception:** `LDSCInputError`
**Symptom:** `Summary-statistics liftover from <source> to <target> using <method> dropped all rows...`

**Likely causes & how to check** (most probable first):

| # | Likely cause | How to check |
|---|--------------|--------------|
| 1 | The declared source build does not match the input CHR/POS coordinates | Compare several input positions against the declared `--source-genome-build` |
| 2 | CHR/POS columns are missing, malformed, or parsed from the wrong fields | Inspect the dropped-SNP sidecar and the first rows of the raw file |
| 3 | The HM3 quick-liftover map cannot cover this non-HM3 SNP set | Check whether the input SNPs are HapMap3-like; use chain-file liftover for broader coverage |
| 4 | Source or target coordinates collide after liftover | Inspect the drop sidecar for `source_duplicate` or `target_collision` reasons |
| 5 | The chain file is for the wrong direction or assembly pair | Confirm the chain filename/source-target pair matches the command flags |

**Remedies:**

1. Fix the source CHR/POS coordinates or pass the correct `--source-genome-build`.
2. Use an explicit chain file for non-HM3 SNP sets instead of HM3 quick liftover.
3. Review the dropped-SNP sidecar to identify whether missing coordinates,
   unmapped variants, or duplicate coordinates removed the rows.
