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

## ldscore

### ldscore: no annotation SNPs remain after reference-panel intersection

**Raised by:** `ldscore_calculator._align_annotation_bundle_to_ref_panel()`,
`_kernel.ldscore.compute_chrom_from_parquet()`, and
`_kernel.ldscore.compute_chrom_from_plink()` · **Exception:** `LDSCInputError`
**Symptom:** `ldscore retained no annotation SNPs on chromosome <chrom> after ... intersection`

**Likely causes & how to check** (most probable first):

| # | Likely cause | How to check |
|---|--------------|--------------|
| 1 | Annotation SNP identifiers use a different identity space than the reference panel | Compare a few annotation SNP IDs against reference-panel metadata; confirm both use the same `--snp-identifier` |
| 2 | Annotation coordinates and reference-panel coordinates are on different genome builds | Check `--genome-build` and the build recorded in R2 parquet or reference-panel metadata |
| 3 | Allele-aware mode was requested but annotation/reference alleles do not match | Inspect A1/A2 columns in annotation and reference metadata for a few expected overlapping SNPs |
| 4 | A reference-panel SNP restriction removed all overlapping annotation SNPs | Re-run without `--ref-panel-snps-file` / HM3 restriction to confirm overlap returns |
| 5 | Chromosome labels or chromosome-sharded path tokens resolved to mismatched chromosomes | Print CHR values from the annotation and reference metadata for the failing chromosome |

**Remedies:**

1. Regenerate annotation and reference-panel artifacts with the same
   `--snp-identifier` and `--genome-build`.
2. Use an allele-aware SNP identifier only when both annotation and reference
   metadata carry matching A1/A2 columns.
3. Relax or rebuild the reference-panel SNP restriction so it overlaps the
   annotation SNP universe.

### ldscore: parquet R2 input is incompatible

**Raised by:** `_kernel.ldscore.SortedR2BlockReader`, sidecar-binding guards,
and parquet row decoders · **Exception:** `LDSCInputError`
**Symptom:** `ldscore could not use R2 parquet ...` / `... sidecar ... missing`

**Likely causes & how to check** (most probable first):

| # | Likely cause | How to check |
|---|--------------|--------------|
| 1 | The R2 parquet was written by an older LDSC version or another tool | Inspect columns; current files contain `IDX_1`, `IDX_2`, `R2`, and `SIGN` |
| 2 | The matching `chrN_meta.tsv.gz` sidecar is missing or was copied from another panel | Confirm each `chrN_r2.parquet` has a same-directory `chrN_meta.tsv.gz` with matching timestamps/provenance |
| 3 | Parquet schema metadata was stripped or edited | Inspect parquet metadata for `ldsc:sorted_by_build`, `ldsc:n_snps`, and `ldsc:sidecar_identity_sha256` |
| 4 | The R2 directory mixes chromosomes or builds from different reference-panel runs | List files in the R2 directory and compare recorded `genome_build` metadata |
| 5 | Duplicate or conflicting pair rows survived in the R2 artifact | Regenerate from a deduplicated reference-panel source and current `ldsc build-ref-panel` |

**Remedies:**

1. Regenerate the reference panel with the current `ldsc build-ref-panel`.
2. Keep each `chrN_r2.parquet` with its matching `chrN_meta.tsv.gz` sidecar; do
   not hand-edit or independently copy sidecars.
3. Use one consistent R2 directory per genome build and reference-panel run.

### ldscore: genome build could not be resolved consistently

**Raised by:** `ldscore_calculator._resolve_ldscore_chr_pos_genome_build()` and
`ldscore_calculator._infer_r2_dir_genome_build()` · **Exception:** `LDSCInputError`
**Symptom:** `ldscore could not infer the genome build...` / `ldscore found conflicting genome-build evidence...`

**Likely causes & how to check** (most probable first):

| # | Likely cause | How to check |
|---|--------------|--------------|
| 1 | `--genome-build auto` had no annotation sample or R2 parquet build metadata to inspect | Confirm annotation path tokens resolve and R2 parquet metadata contains `ldsc:sorted_by_build` |
| 2 | Annotation inputs and R2 parquet metadata were generated on different genome builds | Compare annotation coordinate build against the R2 parquet `ldsc:sorted_by_build` value |
| 3 | The R2 directory contains parquet files from multiple builds | Inspect build metadata across `chr*_r2.parquet` files in the directory |
| 4 | Annotation CHR/POS columns were parsed from the wrong fields | Print the annotation header and first rows; confirm CHR/POS look like genomic coordinates |

**Remedies:**

1. Pass an explicit build: `--genome-build hg19` or `--genome-build hg38`.
2. Regenerate annotation and R2 reference-panel artifacts on the same genome build.
3. Keep only one build's R2 parquet files in a given `--r2-dir`.

## build-ref-panel

### build-ref-panel: no reference-panel artifacts were produced

**Raised by:** `ref_panel_builder.ReferencePanelBuilder.run()` and reference-panel
metadata cleanup paths · **Exception:** `LDSCInputError`
**Symptom:** `build-ref-panel produced no chromosome artifacts...` / `...retained no parquet metadata rows...`

**Likely causes & how to check** (most probable first):

| # | Likely cause | How to check |
|---|--------------|--------------|
| 1 | The SNP restriction does not overlap the PLINK source panel | Compare the first IDs/coordinates in `--ref-panel-snps-file` against the `.bim` file |
| 2 | The restriction or PLINK panel is on the wrong genome build | Check `--source-genome-build` and any build-specific restriction position columns |
| 3 | Identity cleanup dropped all rows because identifiers are missing or duplicated | Inspect `diagnostics/dropped_snps/chr*_dropped.tsv.gz` for `reason` values |
| 4 | Liftover or duplicate-position filtering removed every retained SNP | Inspect the dropped-SNP sidecar for `unmapped_liftover`, `source_duplicate`, or `target_collision` |
| 5 | MAF or individual filtering removed all SNPs before artifact writing | Relax `--maf-min` or check the keep-individual file against the `.fam` file |

**Remedies:**

1. Build the SNP restriction from the same PLINK source build and SNP identifier mode.
2. Review `diagnostics/dropped_snps/` to identify the first filter that removed rows.
3. Relax filters or regenerate the PLINK/reference inputs so at least one SNP remains
   per chromosome.

### build-ref-panel: liftover or genetic-map configuration is incomplete

**Raised by:** `ref_panel_builder.ReferencePanelBuilder._prepare_build_state()`,
chromosome liftover setup, and genetic-map interpolation helpers · **Exception:** `LDSCUsageError` / `LDSCInputError`
**Symptom:** `build-ref-panel cannot emit <build>...` / `...genetic map...`

**Likely causes & how to check** (most probable first):

| # | Likely cause | How to check |
|---|--------------|--------------|
| 1 | `--ld-wind-cm` was requested without the source-build genetic map | Confirm the matching `--genetic-map-hg19-sources` or `--genetic-map-hg38-sources` option is present |
| 2 | Liftover enabled an opposite-build output without that build's genetic map | Check whether a chain file or HM3 quick liftover emits the other build |
| 3 | The genetic map is for the wrong build or lacks the failing chromosome | Inspect CHR values in the map and compare them to `--source-genome-build` |
| 4 | Genetic map shards overlap or are not sorted by chromosome/position | Sort and deduplicate map rows by CHR/POS |
| 5 | A chain file is missing, reversed, or incompatible with the source/target pair | Check the chain filename and source-target CLI flag direction |

**Remedies:**

1. For cM windows, provide every genetic map needed by the emitted build(s).
2. Use `--ld-wind-kb` or `--ld-wind-snps` when genetic maps are unavailable.
3. Use liftover chains only in chr_pos-family modes and make the chain direction
   match the source and target builds.

### build-ref-panel: SNP restriction does not match the source panel

**Raised by:** `ref_panel_builder._read_ref_panel_snp_restriction()` and
restriction build/column readers · **Exception:** `LDSCInputError`
**Symptom:** `build-ref-panel SNP restriction does not match the source panel build...`

**Likely causes & how to check** (most probable first):

| # | Likely cause | How to check |
|---|--------------|--------------|
| 1 | The restriction's generic POS column is on a different build than the PLINK source | Infer or inspect several CHR/POS rows against hg19/hg38 |
| 2 | The restriction file has ambiguous build-specific position columns | Print the header and keep only one `hg19_POS` or `hg38_POS` source-build column |
| 3 | CHR/POS columns were parsed from the wrong delimiter or header names | Inspect the header with `head -1 <file> | cat -A` |
| 4 | Allele-aware restriction rows are missing A1/A2 values | Check whether every row has both allele columns when using allele-aware modes |
| 5 | The file contains too little HM3 overlap for automatic build inference | Use explicit source-build-specific position columns instead of generic `POS` |

**Remedies:**

1. Prepare the restriction file on the same build as the PLINK source panel.
2. Prefer source-build-specific position columns such as `hg19_POS` or `hg38_POS`.
3. Fix ragged rows or delimiter/header issues before rerunning.

### build-ref-panel: reference-panel artifact is incompatible

**Raised by:** `_kernel.ref_panel.ParquetR2RefPanel` and metadata sidecar readers
· **Exception:** `LDSCInputError`
**Symptom:** `Reference-panel metadata sidecar is missing...` / `...required identity/provenance keys are missing...`

**Likely causes & how to check** (most probable first):

| # | Likely cause | How to check |
|---|--------------|--------------|
| 1 | The R2 parquet was copied without its matching `chrN_meta.tsv.gz` sidecar | Confirm each `chrN_r2.parquet` has a same-directory `chrN_meta.tsv.gz` |
| 2 | The artifact was written by an older LDSC version | Inspect parquet or sidecar metadata for `schema_version` and `artifact_type` |
| 3 | A parent directory with multiple build children was passed without a concrete build | Check for `hg19/` and `hg38/` children under `--r2-dir` |
| 4 | The metadata sidecar lacks required SNP identity columns for the active mode | Inspect sidecar columns for SNP or CHR/POS plus A1/A2 in allele-aware modes |
| 5 | The selected chromosome/build was not generated | List `chr*_r2.parquet` files in the selected build directory |

**Remedies:**

1. Keep R2 parquet files and metadata sidecars together as one artifact family.
2. Regenerate the reference panel with the current `ldsc build-ref-panel`.
3. Pass a concrete build-specific R2 directory or set the matching genome build.
