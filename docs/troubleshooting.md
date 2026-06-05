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
