# Troubleshooting

Last updated on: 2026-10-29

This reference expands run-aborting `ldsc` errors with at least three distinct plausible causes, following the [current logging and error-message guide](../dev/design/logging.md#diagnostics-and-troubleshooting).
It also retains guidance for failed-run artifacts, documented query omissions, and cleanup warnings that users need to interpret.
Simple validation errors should explain their cause and repair in the primary message rather than require a separate reference entry.
Entries are organized by the command you run.

Start with the symptom, the affected input or query, and any counts in the actual error.
The listed causes are alternatives to check, not a diagnosis of your run.
Choose a repair supported by that evidence; do not assume that the first listed cause applies.
Module names and exception types appear in supplementary developer references after the user guidance.

If the error names a diagnostic table, inspect its reported inputs, reasons, and repairs.
A file may not exist when failure occurs before that diagnostic or workflow log is written; an older file in the same directory may describe a different run.
Use the available error message first rather than assume every failure produced all the files mentioned here.
Known input problems should be understandable without selecting `DEBUG`; saved tracebacks and additional `DEBUG` detail support investigation of unexpected failures.

Use a fresh output directory for an independent retry.
Use `--overwrite` only when replacing existing owned outputs is intended.
Changing SNP restrictions, filters, model inputs, or query-continuation settings can change the analysis or the set of published results; review that consequence before treating such a change as a repair.

When a message includes `docs/user-guide/troubleshooting.md#<section>`, its fragment identifies a heading below.

## `infer-genome-build`: input cannot be inspected

The helper reports readable but insufficient build evidence as `unknown` with exit code 0.
An inspection error instead means the input cannot be read as a supported artifact, or available evidence is contradictory.

| Symptom | Check and repair |
| --- | --- |
| Missing or unrecognized path | Supply one existing exact file or recognized LDSC3 artifact directory; patterns and PLINK prefixes are not expanded. |
| Unreadable, corrupt, binary, or malformed input | Check permissions, compression, and format; use the `.bim` rather than a binary PLINK `.bed` for coordinate inspection. |
| Unrecognized metadata or directory | Supply a supported table directly or restore the canonical artifact metadata and layout. |
| Missing member or member outside the artifact | Restore the referenced member inside the artifact directory; the helper does not follow upstream paths or escaping links. |
| Conflicting builds | Inspect the paths, chromosomes, and metadata fields named in the error and supply consistent coordinates/declarations for that coordinate set. |

Build declarations are reported without coordinate verification; a successful report is not a full artifact-integrity check.
Deliberately separate hg19/hg38 coordinate alternatives are reported separately.
There is no override flag, output directory, or log file for this command.
See the [helper contract](../reference/snp-identity-and-genome-build.md#standalone-build-inspection).

## Staged input validation

**Symptom:** a combined error lists multiple input objects or chromosomes before staging or numerical work proceeds.

| Cause | Check and repair |
| --- | --- |
| Missing path or PLINK companion | A path-resolution failure means no file matched on the executing node; verify the filename, location, and visibility there. Read every `source` and `chrom` in `input_issues.tsv`; restore each selected BED/BIM/FAM trio. Quote globs and preserve dotted prefixes. |
| Unreadable input | Use the exception type and operating-system error number in `details` to check the named path and access permissions on the executing node. |
| Wrong header, delimiter, or compressed/Parquet file | Inspect the named file's header and format; restore canonical fields and readable compression. Annotation checks use the existing alias and identity rules. |
| Missing or incompatible index component | Restore the complete immutable index with matching chromosome/root metadata; never combine components from separate indexes. |
| Declared chromosome coverage differs from content | Use matching annotation/reference suites and a covered catalog selection. An explicit `@` requires the workflow's documented chromosome suite; a glob selects actual matches. |
| Later alignment, identity, or support failure | These checks need scanned rows or filtered intersections. Use the existing chromosome-scope, gene-list, and SNP-drop audits; passing the earlier path gate cannot establish them. |

Writing workflows retain the six fields `input_role`, `source`, `chrom`, `reason`, `details`, and `repair` in `diagnostics/input_issues.tsv`; gene-index construction uses `.<index>.build-state/input_issues.tsv`. Python exceptions also expose the repair frame as `input_issues`. Existing specialized diagnostics remain authoritative for gene rows, PLINK contents, and alignment. Correct all reported defects at the current gate before rerunning; later content gates may expose defects that could not be established from headers.

For `munge-sumstats`, a missing `--raw-sumstats-file` match is recorded as `missing_required_input` before summary-statistic reading begins.
The message gives a path-specific repair rather than suggesting column or delimiter changes for this failure.
Captured input exceptions remain chained in the workflow log, including all captured causes when multiple inputs fail.
Failed declaration-validation phases report the number of distinct inputs with issues, while unknown phase totals are omitted from progress denominators.

Output collisions are a separate error: passing `--overwrite` authorizes owned output replacement but never bypasses input validation. Early input failures preserve existing scientific files and retain the established authorized-overwrite `RUN_FAILED.txt` contract. INFO phase records distinguish validation, staging, computation, and publication; periodic progress occurs when bounded chunks return, not during a blocking third-party call.

**Developer reference:** `_runtime.preflight.InputGate` and workflow content gates.

## Common

### Common: `RUN_FAILED` is present in an output directory

The latest authorized overwrite failed. Read the marker first, then inspect the
detailed log it names when a log was opened. The directory may contain
incomplete or mixed artifacts because LDSC3 preserves each workflow's existing
write order and performs no rollback or restoration. Correct the underlying
error and rerun the same materializing command with `--overwrite`; a successful
retry removes the applicable marker. Do not treat the marker as a scientific
result file.

Markers, logs, audits, and scientific outputs use the same expanded destination. A marker alone does not prevent retry. If an older version placed logs or markers under a literal `$VARIABLE` or `~` directory, update the package and rerun using the intended destination; existing misplaced directories are not migrated or deleted automatically. Regression coverage is in [`test_output_directory_retries.py`](../../tests/runtime/test_output_directory_retries.py).

### Common: input path did not resolve to one file

**Symptom:** `Could not resolve <label> path from token '<token>': matched 0 files` / `matched N files`

**Likely causes & how to check** (most probable first):

| # | Likely cause | How to check |
|---|--------------|--------------|
| 1 | The path is misspelled or relative to a different working directory | `pwd`; then `ls <path>` from the same shell |
| 2 | A glob is too broad for an input that must be one file | `python3 -c "import glob; print(glob.glob('<pattern>'))"` |
| 3 | An ordinary file input uses an unsupported bare prefix (PLINK prefixes are supported) | Check whether the token looks like `chr@` / `.@.` rather than a bare prefix |
| 4 | The file exists only with a suffix LDSC3 does not infer | `ls <prefix>*`; pass the full filename including suffix |
| 5 | A PLINK prefix is incomplete | Confirm all three files exist: `<prefix>.bed`, `<prefix>.bim`, `<prefix>.fam` |

**Remedies:**

1. Pass the exact existing file path when the command expects a single file.
2. Narrow broad globs so they match exactly the intended file.
3. For annotation chromosome suites, use an explicit `@` token such as `baseline.@.annot.gz` or a quoted glob with the full file suffix.

**PLINK resolution (`ldscore`, `build-gene-ldscore-index`, `build-r2-matrix`):** A plain stem such as `1000G.EUR.QC.` is supported and preserves numeric suffixes such as `.22`. These commands share `_runtime.paths.inspect_plink_inputs()`. If `Could not resolve PLINK inputs` appears, read every listed issue: restore each missing BED/BIM/FAM member, replace malformed or truncated trios, verify BIM chromosome contents against any explicit `@` declaration, and narrow inputs when multiple trios contain the same chromosome. Renaming a file does not change its chromosome assignment. Gene-index construction requires chromosomes 1-22 even without `@`; direct query LD-score full-suite declarations also retain their coverage checks. Direct LD-score issues appear in its input diagnostics; gene-index issues are retained in `.<output-name>.build-state/plink_input_issues.tsv`; R²-builder issues appear in `diagnostics/plink_input_issues.tsv`. See the [shared resolution contract](../reference/input-paths.md#plink-prefix-resolution).

**Developer reference:** `_runtime.paths.resolve_scalar_path()` and group/path-prefix resolvers · **Exception:** `LDSCInputError`

### Common: output artifact already exists

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

**Developer reference:** `_runtime.paths.preflight_output_artifact_family()` · **Exception:** `FileExistsError`

### Common: genome build could not be inferred

**Symptom:** `Could not infer genome build for <context>...`

**Likely causes & how to check** (most probable first):

| # | Likely cause | How to check |
|---|--------------|--------------|
| 1 | The CHR/POS input has too little overlap with HapMap3 reference SNPs | Count matched/reference-like variants; small custom SNP lists often cannot infer a build |
| 2 | The input coordinates are from another genome build or a custom reference | Inspect several known SNP coordinates against hg19/hg38 in a genome browser or source metadata |
| 3 | CHR/POS columns were parsed from the wrong columns or wrong delimiter | Print the header and first rows; confirm `CHR` and `POS` values look like chromosome labels and base-pair positions |
| 4 | Coordinates are mixed across builds or coordinate bases | Check whether some rows match hg19 and others match hg38/0-based positions |
| 5 | A dense input (for example a full imputed GWAS or 1000 Genomes panel) has chance coordinate coincidences with the wrong hypotheses | Read the `Hypothesis support` list in the error, as described below |

**Reading the hypothesis support:** when the best hypothesis falls below the 99% threshold, the error lists matched counts and fractions for all four genome build and coordinate-basis hypotheses (`hg19/1-based`, `hg19/0-based`, `hg38/1-based`, `hg38/0-based`).
If one hypothesis is just below 99% and the remainder is spread thinly across the other three (for example 98.8%, 0.4%, 0.5%, 0.3%), the shortfall is most likely chance coordinate coincidences in a dense input, and the best hypothesis is likely the true genome build.
If a second hypothesis carries substantial support (for example 60% and 40%), the file most likely mixes genome builds or coordinate bases, or was assembled from differently built sources.
In either case, confirm the genome build from the source metadata before declaring it.

**Remedies:**

1. If the input build is known, declare it with the option for that command.
   Use `--source-genome-build hg19` or `hg38` for coordinate-based `munge-sumstats` inputs and `build-r2-matrix` source panels; use `--genome-build hg19` or `hg38` for `annotate` and direct `ldscore`.
   Do not guess a build solely to pass validation.
2. Fix column mapping or delimiter issues before retrying automatic inference.
3. Rebuild the input so all CHR/POS rows use one genome build and one coordinate basis.

**Developer reference:** `_genome.build_inference.resolve_genome_build()` and `_genome.build_inference.resolve_chr_pos_table()` · **Exception:** `LDSCInputError`

### Common: genome-build consistency check failed

This illustrative report follows the [genome-build policy](../reference/snp-identity-and-genome-build.md#genome-build-policy).
The workflow examines all relevant inputs before raising a consistency error and includes each selected file or shard in the grouped findings.

```text
Genome-build consistency check failed for ldscore.

Build evidence across inputs:
  hg19:
    --plink-prefix: panel/chr1.bim (coordinate inference)
    --baseline-annot-sources: annotations/base.1.annot.gz (coordinate inference)
  hg38:
    --baseline-annot-sources: annotations/base.2.annot.gz (coordinate inference)
  unresolved:
    --regr-snps-file: keep.tsv (too few informative coordinates)
  unchecked:
    --query-annot-bed-sources: regions.bed (ordinary BED; inherits workflow build)

The reference input and baseline annotation shards contain conflicting
hg19/hg38 build evidence. The build of keep.tsv could not be resolved.
regions.bed was not independently checked and does not establish a build.

Check the listed inputs and ensure their coordinates use a consistent build.
For verified unlabeled inputs, specify --genome-build hg19 or --genome-build hg38.
Explicit builds cannot override conflicting declarations. Supply matching coordinates,
correct erroneous metadata only after verification, or convert coordinates explicitly.
No coordinates or build labels were changed.
```

The same input flag can appear in multiple groups when its files disagree.
Actual reports identify the files and evidence found in that run, explain each unresolved input, and include matching inputs as well as conflicting ones.
Metadata declarations must be labeled as declarations rather than coordinate inference.
An unchecked ordinary BED file is allowed to inherit a successfully resolved workflow build; the unchecked label is not itself an error.
A generic coordinate keep-list that lacks sufficient inference evidence is an error under `auto`, even if all resolved inputs agree; it does not inherit their build.
rsID-only keep-lists used for rsID filtering need no build validation and can be listed as `unchecked` with the reason "rsID filtering; no coordinate build required".
Valid build metadata and matching build-specific position columns supply build evidence without coordinate re-inference.
When a matching build-specific column is present, it takes precedence over generic `POS`; an unused generic column does not create a conflict.

Override advice must use the option for the affected command.
For `build-r2-matrix` and coordinate-based `munge-sumstats`, the source-build example is `--source-genome-build hg19` or `--source-genome-build hg38`.
`--output-genome-build` specifies munging's target build and does not override source evidence.
For commands without a build override, direct the user to consistent upstream artifacts instead of suggesting an unsupported flag.
Explicit build declarations remain authoritative within their existing scope; they do not repair unreadable files or other structural input errors.

### Common: keep-list build could not be resolved

This illustrative example shows insufficient evidence in a generic coordinate keep-list.
A small list is still an error under `auto`, even when the reference panel is resolved.
For the full distinction between rsID filtering, declared coordinates, and generic coordinates, see the [SNP keep-list rules](../reference/snp-identity-and-genome-build.md#snp-keep-list-rules).

```text
Genome-build consistency check failed for ldscore.

Build evidence across inputs:
  hg19:
    --plink-prefix: panel/chr1.bim (coordinate inference)
  hg38: none
  unresolved:
    --regr-snps-file: keep.tsv (40 informative HM3 matches; at least 200 required)
  unchecked: none

The generic CHR/POS keep-list keep.tsv could not be assigned a genome build.
The reference panel's build does not establish the keep-list's build.
Check that the keep-list coordinates use the intended workflow build, then
explicitly specify --genome-build hg19 or --genome-build hg38 for this run.
The explicit declaration is authoritative and does not convert coordinates.
```

The error and workflow log must identify the actual keep-list flag, file, and reason inference failed, with all other inputs included in the grouped report.
The explicit build option allows the user to declare a known build for a small coordinate list without supplying additional inference evidence.
For panel building or munging, replace the example override with the command's existing `--source-genome-build` option.
Inference uses bounded samples from every relevant file/chromosome and retains the existing confidence criteria; an exhaustive coordinate scan is not required.

### Common: LDSC3 artifact schema or provenance is incompatible

**Symptom:** `Could not read LDSC <artifact> artifact metadata...`

**Likely causes & how to check** (most probable first):

| # | Likely cause | How to check |
|---|--------------|--------------|
| 1 | The artifact was written by an older LDSC3 version | Inspect the artifact metadata sidecar or parquet footer for `artifact_type`, `snp_identifier`, and `genome_build` |
| 2 | The artifact type does not match the loader | Confirm the file was produced by the command you are now trying to load from |
| 3 | The metadata sidecar was copied without its matching data file, or vice versa | Compare file modification times and paths for the artifact plus sidecar |
| 4 | SNP identity provenance was hand-edited or corrupted | Inspect `snp_identifier`, `genome_build`, and identity metadata fields |

**Remedies:**

1. Regenerate the artifact with the current LDSC3 version and the same command family.
2. Keep LDSC-generated data files and their sidecars together; do not hand-edit them.
3. If you need a different `snp_identifier` or `genome_build`, regenerate from the upstream input.

**Developer reference:** `_core.snp_identity.validate_identity_artifact_metadata()` and artifact reload guards · **Exception:** `LDSCInputError`

## `munge-sumstats`

### `munge-sumstats`: could not map a required column

**Symptom:** `munge-sumstats could not map the required column '<field>' from the header...`

**Likely causes & how to check** (most probable first):

| # | Likely cause | How to check |
|---|--------------|--------------|
| 1 | The header uses a name LDSC3 does not recognize | Read the decompressed header and compare it with the field named in the error; `ldsc munge-sumstats --help` lists column-hint options |
| 2 | The column exists under a synonym you must declare | Re-run with an explicit hint, e.g. `--snp MarkerName --a1 Allele1` |
| 3 | Wrong delimiter, so the whole header parsed as one column | Inspect the decompressed header; check whether column names are separated as the selected input format expects |
| 4 | Wrong `--input-format`, so expected columns differ | Confirm `--input-format` matches the file; `auto` and an explicit profile share aliases, optional-field rules, and validation |

**Remedies:**

1. Map the missing field to its actual header name, for example `--snp MarkerName`, `--chr CHROM`, `--pos BP`, `--a1 Allele1`, `--a2 Allele2`, or `--p PVAL`.
   Sample-size columns use `--N-col` or the paired `--N-cas-col` / `--N-con-col` options; signed statistics use `--signed-sumstats COLUMN,NULL`.
2. Check the same command's mapping with `--infer-only`, which requires `--raw-sumstats-file` and `--output-dir` but creates no outputs.
   Coordinate-based modes also require `--output-genome-build`.
   See [input inspection](../reference/munge-sumstats.md#--infer-only).
3. If the header or delimiter is malformed, prepare a corrected input copy and select the appropriate `--input-format`.
   Explicit hints cannot repair a file whose fields were not parsed into separate columns.

**Developer reference:** `_sumstats.munge_kernel` column resolution · `SumstatsTable.validate()` · **Exception:** `LDSCInputError`

### `munge-sumstats`: multiple or conflicting sample-size column strategies

**Symptom:** `munge-sumstats found multiple sample-size strategies...`, an incomplete `--N-cas-col`/`--N-con-col` pair, or a conflict with `--N-col`.

**Likely causes & how to check:**

| # | Likely cause | How to check |
|---|--------------|--------------|
| 1 | The header auto-maps both direct N and case/control counts | Inspect the header for `N` (or another direct-N alias) together with aliases such as `NCAS` and `NCON` |
| 2 | Both explicit strategies were supplied | Check whether the command contains `--N-col` plus `--N-cas-col`/`--N-con-col` |
| 3 | Only one case/control column flag was supplied | Confirm both `--N-cas-col` and `--N-con-col` are present |

**Remedies:**

1. Choose direct N with `--N-col <column>`; inferred case/control columns are
   suppressed with a warning.
2. Or choose case/control N with
   `--N-cas-col <cases> --N-con-col <controls>`; inferred direct N is
   suppressed with a warning.
3. Do not combine the two strategies. For `NEFF + NCAS + NCON`, use
   `--N-col NEFF` when those exact effective-N values are intended; no
   `--ignore` flag is required.

**Developer reference:** `_sumstats.munge_kernel` sample-size column resolution · **Exception:** `LDSCInputError` or `LDSCUsageError`

### `munge-sumstats`: no SNPs remain after filtering

**Symptom:** `munge-sumstats removed every SNP...` / `...no SNPs remain after SNP keep-list restriction`

Munging restricts to packaged HM3 by default.
First use the error and log counts to identify whether QC, the SNP restriction, or a later stage removed the rows.
To intentionally process all input SNPs subject to ordinary QC, use `--no-snp-restriction`; to intentionally replace HM3 with a custom list, use `--sumstats-snps-file FILE`.
Both choices change the requested SNP population.
For cross-build output, either override requires `--liftover-chain-file`; default HM3 uses automatic quick liftover unless a chain is supplied.
Coordinate-based identity requires an explicit `--output-genome-build`; if source-build inference fails, declare the known input build with `--source-genome-build`.
Concrete source/output build options do not apply to rsID-family munging.

**Likely causes & how to check** (most probable first):

| # | Likely cause | How to check |
|---|--------------|--------------|
| 1 | The input has no overlap with the default HM3 set, or a custom keep-list uses a different SNP-id space (rsID vs chr:pos) | Compare the first IDs of each file; confirm both use `snp_identifier=<mode>` |
| 2 | Keep-list and summary statistics are on different genome builds | Compare the `genome_build` in the error line against the keep-list build |
| 3 | Imputation-quality (INFO), minor-allele-frequency (MAF), or sample-size (N) filters removed every row | Inspect filter counts in the error/log and verify the column meanings and intended thresholds; earlier QC removals are not in the dropped-SNP sidecar |
| 4 | Input is effectively empty after a delimiter/format mis-parse | `zcat <file> \| wc -l` |

**Remedies:**

1. Regenerate the keep-list in the same `snp_identifier` mode and genome build as the summary statistics.
2. For QC removals, inspect the error/log counts and confirm the mapped quality, frequency, and sample-size columns.
   Correct a mistaken mapping or threshold; change an intended QC threshold only after deciding that the new criterion is appropriate for the analysis.
3. Earlier QC and SNP-restriction removals have counts rather than individual rows in the dropped-SNP sidecar.
   Use `diagnostics/dropped_snps/dropped.tsv.gz` for liftover or final identity-cleanup removals only when that file was written.
   See [whole-run row accounting](../reference/munge-sumstats.md#preparation-and-run-accounting).

**Developer reference:** `_sumstats.munge_kernel` (post-filter and keep-list paths) · **Exception:** `LDSCInputError`

### `munge-sumstats`: munged summary statistics are malformed or outdated

**Symptom:** `Cannot load munged summary statistics at '<path>': ...` (old schema / bad provenance / missing A1-A2 / duplicate identity rows)

**Likely causes & how to check** (most probable first):

| # | Likely cause | How to check |
|---|--------------|--------------|
| 1 | Artifact predates the current self-describing parquet schema | Inspect the parquet footer for `ldsc:artifact_type`, `ldsc:snp_identifier`, and `ldsc:genome_build` |
| 2 | The artifact is footerless Parquet | Re-munge from the raw GWAS input to produce self-describing `sumstats.parquet`; footerless Parquet is not a legacy compatibility format |
| 3 | Identity provenance in the parquet footer is invalid/corrupt | Inspect the footer identity fields against the current contract |
| 4 | An allele-aware `snp_identifier` artifact lacks A1/A2 columns | `python3 -c "import pandas; print(pandas.read_parquet('<f>').columns)"` |
| 5 | Duplicate/invalid SNP-identity rows survived in the artifact | Re-munge from raw input; the loader reports the dropped-row reasons |

**Remedies:**

1. Re-run `ldsc munge-sumstats` from the **raw** GWAS file to regenerate the
   artifact with the current schema.
2. A genuine legacy `.sumstats` or `.sumstats.gz` text artifact is accepted
   directly by regression when it contains `SNP`, `A1`, `A2`, `Z`, and `N`;
   see `docs/reference/formats/legacy-sumstats.md`.
3. Do not hand-edit curated `.sumstats`/`.parquet` artifacts; treat them as outputs.

**Developer reference:** `_sumstats.munger.load_sumstats()` and its metadata helpers · **Exception:** `LDSCInputError`

### `munge-sumstats`: liftover dropped all rows

**Symptom:** `Summary-statistics liftover from <source> to <target> using <method> dropped all rows...`

**Likely causes & how to check** (most probable first):

| # | Likely cause | How to check |
|---|--------------|--------------|
| 1 | The declared source build does not match the input CHR/POS coordinates | Compare several input positions against the declared `--source-genome-build` |
| 2 | CHR/POS columns are missing, malformed, or parsed from the wrong fields | Inspect the count-bearing error and the first rows of the raw file |
| 3 | The selected mapping cannot resolve the retained source coordinates | Check the selected method and source/output build pair; verify the chain direction if overriding packaged HM3 mapping |
| 4 | Source or target coordinates collide after liftover | Inspect duplicate-source and duplicate-target counts in the error/log |
| 5 | The chain file is for the wrong direction or genome-build pair | Confirm the chain filename/source-target pair matches the command flags |

**Remedies:**

1. Fix the source CHR/POS coordinates or pass the correct `--source-genome-build`.
2. For broader SNP coverage, set `--no-snp-restriction` or a custom `--sumstats-snps-file`, and supply a source-to-output chain for cross-build conversion. Supplying a chain alone changes the mapping method while retaining the default HM3 restriction.
3. Review the error/log counts to identify whether missing coordinates, unmapped variants, or duplicate coordinates removed the rows. After a successful run with partial removals, use the dropped-SNP sidecar to inspect the affected rows.

Successful runs report the method and counts in stdout and `diagnostics/sumstats.log`, including at `--log-level ERROR`; failed runs retain the error and available counts in the log and do not print a success summary. An all-dropped liftover fails before writing the dropped-SNP sidecar, so inspect its count-bearing error first. Mapping input is measured after earlier QC and keep-list filtering. See the [liftover behavior table](../reference/munge-sumstats.md#liftover-rules).

**Developer reference:** `_genome.liftover.apply_sumstats_liftover()` · **Exception:** `LDSCInputError`

## `annotate`

### `annotate`: grid selection and downstream alignment

Use exactly one of `--ref-snp-annot-sources` or `--ref-snp-bim-sources`.
The old standalone `--baseline-annot-sources` flag is removed; downstream `ldscore` retains it.
BIM input is a six-column headerless `.bim`/`.bim.gz` file with positive integer one-based BP, not a PLINK prefix, and requires no BED/FAM companions.
For an `@` BIM suite, provide all 22 autosomal members with matching chromosome contents; use exact paths or globs for a deliberately smaller scope.
BIM format and member failures are recorded in `diagnostics/input_issues.tsv`.

If `ldscore` reports row mismatch after annotation generation, pair the exported `baseline.<chrom>.annot.gz` files with the exported `query.<chrom>.annot.gz` files through `--baseline-annot-sources` and `--query-annot-sources`; see [`ldscore`: annotation SNP rows do not match across files](#ldscore-annotation-snp-rows-do-not-match-across-files).
Original baselines may have additional rows removed by annotation identity cleanup.
Inspect the input/output/dropped counts in `diagnostics/annotate.log` and `diagnostics/metadata.json` (`snp_rows`).
BIM mode exports query values only; using them as downstream queries still requires actual aligned baseline annotations.

### `annotate`: input preflight

Missing files, invalid annotation values or headers, and incorrect chromosome members are recorded in `diagnostics/input_issues.tsv`. Exact paths and globs select their actual inputs. For gene-list annotation, `@` explicitly requires autosomes 1–22; each member must contain its declared chromosome. Repair all listed inputs before rerunning. Input defects are not interpreted as zero pathway support.

### `annotate`: gene-list preflight

Use exactly one BED or gene-list query route. Gene lists require a readable gene coordinate catalog and explicit nonnegative `--padding-bp`, including `0` for gene bodies. Check `diagnostics/gene_catalog_issues.tsv`, `diagnostics/gene_list_audit.tsv.gz`, and `diagnostics/gene_list_resolution_summary.tsv` for catalog defects, malformed or unreadable lists, naming collisions, and unresolved identifiers. `resolved-only` permits only the established identifier-resolution omission allowlist; it does not bypass malformed rows, missing sources, or incomplete chromosome coverage.

The error message, shared with gene-list `ldscore` runs, gives each source's rejected and submitted row counts; up to ten examples with their file and line, then the number not shown; a suggested filter for `diagnostics/gene_list_audit.tsv.gz`; and the rerun with `--overwrite`.

Standalone `annotate` checks query names from source filenames against validated baseline headers and reserved metadata/control names before preparing the annotation SNP grid.
For example, if baseline column `base` conflicts with `base.txt` and that list contains `UNKNOWN`, the same terminal error reports both problems.
The source summary records `annotation_name_collision`, while the row audit records the unresolved identifier; duplicate source-derived query names retain `duplicate_query_name`.
Rename the conflicting focal source files to change their query names, and repair rejected identifiers when using strict resolution.
`--allow-unresolved-genes` does not permit naming collisions, even for empty lists or lists with no resolved genes.
Invalid headers or an unusable catalog may stop an earlier gate because dependent checks cannot be performed reliably.

After resolution, every selected gene must be covered by validated grid contents.
`diagnostics/chromosome_scope.json` and the gene audit explain missing coverage.
Measured support is then reported as `annotation_snp_count` on cleaned grid SNP rows; reference-panel support remains unevaluated.
Empty and globally unsupported focal queries are skipped, usable siblings continue, and an all-skipped batch fails without publishing a new query family.
All-one annotation columns remain valid.

An authorized overwrite can leave `RUN_FAILED.txt` after failure. Inspect diagnostics and rerun with corrected inputs and `--overwrite`; no rollback to earlier results is promised. Private staging is removed on handled failure or bundle closure. Python callers should use `with run_annotate(...):` or call `close()` after consuming the returned handle; canonical `query.<chrom>.annot.gz` files remain available.

### `annotate`: no annotation SNP rows remain

**Symptom:** `Annotation file '<path>' has no SNP rows...` / `annotate retained no annotation rows after SNP identity cleanup...`

**Likely causes & how to check** (most probable first):

| # | Likely cause | How to check |
|---|--------------|--------------|
| 1 | The selected chromosome is absent from the annotation shards | List the resolved annotation files and confirm their chromosome token or filename contains the requested chromosome |
| 2 | Annotation files are empty or contain only headers | `zcat <annot.gz> \| head` (or `head <annot>`); confirm data rows exist after the header |
| 3 | SNP identity cleanup dropped every row because identities are missing, duplicated, or incompatible | Inspect `diagnostics/dropped_snps/` or the run log for identity-cleanup drop reasons |
| 4 | Baseline and query annotation inputs use mismatched shard sets or a mixed sharded/unsharded layout | Compare resolved baseline and query files; each chromosome should have the same shard structure |
| 5 | Allele-aware mode was requested but annotation rows lack usable A1/A2 values | Inspect the annotation header and first rows for both allele columns and non-missing allele values |

**Remedies:**

1. Use baseline and query annotation inputs that cover the same chromosome shards.
2. Regenerate nonempty annotation files with the identity columns required by the selected mode and at least one numeric annotation column.
   See [annotation column requirements](../reference/formats/column-schema.md); reference-panel CM/MAF, rather than supplied annotation CM/MAF, governs LD-score calculation.
3. Match `--snp-identifier` and `--genome-build` to the annotation identity columns, then rerun.

**Developer reference:** `_annotation.preflight.inspect_annotation_inputs()`, `_annotation.sources._scan()`, and `_annotation.sources.prepare_annotation_sources()` · **Exception:** `LDSCInputError`

## `build-gene-ldscore-index`

### `build-gene-ldscore-index`: output directory is nonempty but invalid

A regular root `RUN_FAILED.txt` file is an owned diagnostic, not an index artifact.
If a failed build leaves only that marker, rerun the corrected command in the same output directory.
Preflight preserves the marker, and a successful build removes it.
Versions before this retry fix incorrectly rejected marker-only directories even with `--overwrite`; update the package rather than deleting the directory.

If the current version still reports this error, the directory contains an incomplete index or other unrecognized contents. The marker does not grant permission to replace those contents. Preserve them for inspection and choose a fresh output directory, or explicitly resolve the conflicting files. A valid existing index still requires `--overwrite`.

### `build-gene-ldscore-index`: baseline/PLINK identifier intersection fails

**Symptom:** the builder reports an empty `rsid` or `chr_pos` baseline/PLINK
intersection before genotype QC, or warns that duplicate groups were dropped.

The builder uses the same identifier-key inner intersection as direct
PLINK-backed `ldscore`. Baseline-only and PLINK-only rows are allowed, dropped,
and counted. For a duplicate effective identity, every row in that source group
is dropped; no representative is selected. Inspect
`<index-dir>/diagnostics/dropped_snps/chrN_dropped.tsv.gz` after success (or the
live log during a failing build). An empty result leaves no reference SNPs for LD scores.
Under rsID matching, coordinate disagreements stop the build with a repair table.
Under coordinate matching, differing SNP labels are reported and the PLINK label is published.
During
a failed build, check
`<parent>/.<index-name>.build-state/build-gene-ldscore-index.log`; after a
successful build, check `<index-dir>/diagnostics/build-gene-ldscore-index.log`.
Confirm that all coordinate sources match the explicit hg19/hg38 build, that `--snp-identifier` matches the intended key, and that the inputs overlap after shared identity cleanup.
There is no build inference or liftover; declarations must agree and cannot be replaced by the command-line assertion.
A failed first build leaves the index destination absent or
empty and can be retried directly; the previous log is archived under
`<parent>/.<index-name>.build-state/history/`.

### `build-gene-ldscore-index`: input conflicts

Build-label conflicts identify the affected catalog, annotation, restriction, or map and its declared/requested builds.
Names such as GRCh37/hg19 and GRCh38/hg38 are equivalent.
Supply matching coordinates, correct an erroneous label only after verifying its coordinates, or perform an explicit conversion before index construction.
The builder never converts coordinates or overwrites a conflicting label.

SNP contradictions produce `.<index-name>.build-state/baseline_plink_issues.tsv` before chromosome workers start.
The table identifies shared keys, baseline/PLINK sources, coordinates and usable normalized allele pairs, and a reason such as `coordinate_mismatch`, `allele_mismatch`, or `cross_chromosome_coordinate_mismatch`.
Verify the reference releases and repair the indicated records; ordinary unmatched SNPs and shared-policy identity cleanup remain separate diagnostics.

**Developer reference:** `_ldscore.gene_index.build._validate_baseline_plink_inputs()` and `_genome.build_policy.GenomeBuildCheck` · **Exception:** `LDSCInputError`

### `build-gene-ldscore-index`: a partial hidden stage remains after failure

**Symptom:** a marked `.<index-name>.stage-<run-id>/` sibling contains one or
more `chromosomes/chrN` directories, but the public index is absent, empty, or
still contains its prior complete version.

`Finished chromosome N` means that shard was durably written inside the private
run transaction; it never means partial chromosome coverage was published. A
graceful failure removes the transaction best-effort. Abrupt termination or a
filesystem cleanup error can retain it, but it is not a restart checkpoint and
must not be copied into the public destination. Retry the full command. The
next invocation removes a recognized computation-only stage before starting a
new transaction. If cleanup warns again, confirm that no builder for the same
destination is active, then remove only the exact marked path named in the
warning. Unrecognized neighboring directories are never cleanup candidates.

### `build-gene-ldscore-index`: transaction cleanup remains after success

**Symptom:** the command finishes successfully but warns that a builder-owned
`.stage-*` transaction path could not be removed, commonly with `Directory not
empty` or `Device or resource busy` on a shared filesystem.

The destination was already replaced and reload-validated. It is a successful,
loadable index; cleanup is garbage collection and the warning reports the exact
retained path. Do not rerun the expensive build solely for this warning. After
confirming no build for the same destination is active, remove the reported
builder-owned path manually or let the next invocation retry recognized
cleanup. Never delete an unmarked neighboring directory based only on a similar
name.

Earlier releases kept logs in a visible `<index-dir>.build/` directory and used a standalone hidden lock, `.<index-name>.build-gene-ldscore-index.lock`.
The current builder ignores both, so it does not wait for a build started by an earlier release.
Let any such build finish before building the same destination, then remove the old directory and lock manually if you no longer need their logs.

### `build-gene-ldscore-index`: required identity or build option is missing

**Symptom:** parsing reports a missing `--snp-identifier`, `--genome-build`, or LD-window option.

Choose one of the four SNP identifier modes, explicitly supply hg19 or hg38, and set exactly one of `--ld-wind-cm`, `--ld-wind-kb`, or `--ld-wind-snps`.
There are no defaults for those decisions.
Whole-chromosome windows are rejected; use a smaller window because this builder has no override.
These settings belong to index construction.
Remove them from `ldsc ldscore --gene-ldscore-index-dir ...`, which inherits its immutable settings.

## `ldscore`

### `ldscore`: an explicit gene LD-score index is missing, corrupt, or incompatible

**Symptom:** indexed LD-score calculation stops because the requested index is missing, unreadable, or incompatible with the command or current reader.

Start with the named file or failed requirement in the error.
Restore missing components from the same complete index rather than combining files from different builds.

Pass the one complete index directory containing root `metadata.json` and
`chromosomes/`. Do not add
live baseline, reference, build, padding, window, map, or region settings: those
belong to the immutable index. A failed validation never falls back to direct
mode. Remove `--padding-bp` entirely rather than passing a zero value. Validation
is completed before canonical scientific output publication. Restore
or rebuild the index, or remove `--gene-ldscore-index-dir` and supply the
full direct-mode inputs explicitly. Older gene-index metadata contracts are not
loaded by the current strict reader and must be rebuilt; ordinary canonical
LD-score directories already produced from them are unaffected.
An `Unsupported gene LD-score index format version 2` or `3` error means the index predates the current format: version 3 introduced exact PLINK LD windows and version 4 stores each chromosome's reference-SNP funnel.
Rebuild it with the current `ldsc build-gene-ldscore-index`.

**Developer reference:** integrity checks include metadata IDs, mode/build agreement, the ordered published-row metadata digest, chromosome coverage, duplicate effective identities, Parquet rows, NPZ members, CSR structure/dtypes, and dimensions.

### `ldscore`: gene-list preflight stopped or a focal query is missing

**Symptom:** a requested gene-list query is absent from LD-score outputs or downstream `partitioned-h2` results, or the run reports that every query was skipped.

Strict gene-list resolution stops before LD-score work if any focal/control row
cannot resolve uniquely. Start with
`diagnostics/gene_list_resolution_summary.tsv`, then filter
`diagnostics/gene_list_audit.tsv.gz` by the affected `source` and
`disposition == 'rejected'`.
Direct and indexed runs also stop here when a query name, taken from its source filename, matches a baseline annotation column, a reserved metadata column, or `gene_control`; the summary records `annotation_name_collision`, and renaming the source file fixes it.
After gene-list validation succeeds, `diagnostics/query_annotation_status.tsv` distinguishes empty lists, lists with no resolved genes, queries with no SNP support, and queries whose LD scores are constant on the SNP rows used for regression.
Support on retained reference SNPs and variation on regression SNPs are different checks: a fully supported query can still have constant LD scores on the regression SNPs.
Skipped queries are omitted from LD-score outputs; usable sibling queries continue, and a requested control never disappears silently.
Use the diagnostic reason and affected source to correct the list, catalog, or incompatible input.
Then choose a fresh output directory or `--overwrite` when replacing the failed run's owned diagnostics is intended.
When all focal queries are skipped, no root scientific result is written.

For the complete reason vocabulary, column definitions, filters, and prioritized
repair steps, see [Gene-list diagnostics and repair](../reference/results/gene-list-diagnostics.md).

**Developer reference:** query annotation status handling

### `ldscore`: chromosome coverage preflight

**Symptom:** the batch fails before publishing scientific outputs.

Read `diagnostics/input_issues.tsv` and `diagnostics/chromosome_scope.json`. For gene lists, inspect `coverage_status`, `selected_genes`, `covered_genes`, `missing_chromosomes`, and `uncovered_gene_ids` in `gene_list_resolution_summary.tsv`; audit rows with `coverage_status == 'uncovered'` identify affected input lines. Fix every reported issue before rerunning with `--overwrite`.

- Missing, unreadable, malformed, or mismatched PLINK trio / R²-sidecar / baseline inputs: restore valid matched artifacts. `@` requires every autosome 1–22. All safely discoverable independent issues are reported at the current gate.
- If a chromosome subset was intended, use quoted `*` patterns or exact paths selecting matching baseline and PLINK chromosome sets, for example `"baseline.*.annot.gz"` and `"panel.*"`. Use `*` to select available files without declaring a complete 1-22 suite. Check the matched files and their contents; a broad wildcard can select unintended chromosomes. For a full-suite run, restore missing files instead of narrowing the declaration. `--r2-matrix-dir` takes a literal directory containing the matching chromosome set, not a glob.
- Baseline/reference chromosome-set disagreement: supply exactly matching validated sets. Glob matches are authoritative; filenames do not establish content. A missing file may be undetectable when the remaining groups consistently cover the same subset.
- Incomplete focal/control coverage: supply matching inputs covering every selected gene, or deliberately revise the lists. Neither `resolved-only` nor chromosome-subset inputs authorize dropping cross-chromosome genes or pathways.
- Invalid immutable index: restore or rebuild the complete public index. Resolution policy never repairs an index or changes its scientific configuration.

Coverage follows unique identifier resolution and explicit gene exclusions, before SNP filtering. Unevaluated support remains blank. Fully covered genes with no retained computational reference SNPs have valid zero-support measurements; focal queries can be skipped for zero support/variance, whereas unusable controls and all-focal-skipped batches fail. See [the complete diagnostic vocabulary](../reference/results/gene-list-diagnostics.md).

Replacing `@` with `*` changes file selection only; it does not remove genes outside the selected chromosome set. The raised error, `diagnostics/ldscore.log`, and input-issue repair text describe these alternatives. Source: [`_ldscore.preflight.validate_direct_scope` / `inspect_direct_inputs`](../../src/ldsc/_ldscore/preflight.py).

For BED queries, `--padding-bp` extends the supplied intervals at both ends. Set it to 0 when your BED files already contain the intended padded regions, to avoid double padding.

**Developer reference:** direct input/coverage validation or immutable index validation.

### `ldscore`: annotation values are malformed

**Symptom:** `Could not parse annotation file ... annotation value columns must be numeric and non-missing`

**Likely causes & how to check** (most probable first):

| # | Likely cause | How to check |
|---|--------------|--------------|
| 1 | A data row has fewer fields than the header | Compare the reported data row with the header and count fields using the file's delimiter |
| 2 | An empty metadata placeholder was collapsed by whitespace parsing | Inspect adjacent delimiters around `CM`; package-written missing metadata uses the literal `NA` token |
| 3 | Delimiters are inconsistent or a row/file is truncated | Inspect the reported row and neighboring rows for mixed tabs/spaces or an incomplete final record |
| 4 | An annotation value is missing or non-numeric | Inspect the reported annotation columns and rows; missing tokens are allowed only in metadata columns |

**Remedies:**

1. Ensure every data row has exactly one field per header column.
2. Use `NA` for missing metadata such as `CM`, never an empty field or a numeric sentinel.
3. Supply numeric, non-missing values in every annotation column, then rerun.

**Developer reference:** `_annotation.tables._validate_annotation_values()` through the annotation and LD-score loaders · **Exception:** `LDSCInputError`

### `ldscore`: annotation SNP rows do not match across files

**Symptom:** `Annotation SNP rows do not match across files: ... row count mismatch` or `... SNP row mismatch at row ...`

Baseline and prebuilt query annotation files read together share one annotation SNP grid: they must list the same SNP rows in the same order.

**Likely causes & how to check** (most probable first):

| # | Likely cause | How to check |
|---|--------------|--------------|
| 1 | The query was built by `annotate --ref-snp-bim-sources`, but the baseline comes from another grid, such as a legacy baseline-LD file | Read `grid_kind` in the query's `annotate` output `diagnostics/metadata.json`; compare row counts of the baseline and query files for one chromosome |
| 2 | The original baseline files were paired with `annotate` queries built on them, but SNP identity cleanup removed rows from the exported grid | Compare `snp_rows` in the `annotate` `diagnostics/metadata.json` with the original baseline row counts |
| 3 | The files are sharded differently, sorted differently, or were written with different SNP identifier modes | Compare the first rows and `CHR` values of the two files named in the error |

**Remedies:**

1. Build the query on the baseline's grid with `ldsc annotate --ref-snp-annot-sources` and pass the exported `baseline.<chrom>.annot.gz` and `query.<chrom>.annot.gz` files together.
2. For a query built on a BIM grid, supply baseline annotations built on that same BIM grid.
3. Regenerate the pair from one `annotate` run with the same `--snp-identifier`.

**Developer reference:** `_aligned_metadata()` in `src/ldsc/_annotation/sources.py` · **Exception:** `LDSCInputError`

### `ldscore`: no reference SNPs remain on the annotation SNP grid

**Symptom:** `ldscore found no reference SNPs on the annotation SNP grid for chromosome <chrom> (<n> annotation rows, <m> reference-panel SNPs)`, or `ldscore retained no reference SNPs on chromosome <chrom>: none of its <n> annotation rows matches ...`

The annotation files define the candidate SNPs; reference-panel SNPs without an annotation row are excluded ([reference SNPs on the annotation SNP grid](../reference/ldscore-snp-universes.md#reference-snps-on-the-annotation-snp-grid)).
A chromosome whose remaining panel SNPs are all off the grid stops the run, because a whole chromosome without one shared SNP usually signals an input mismatch.
Inspect `diagnostics/reference_snp_funnel.tsv` when it was written, to see which step removed the SNPs.

**Likely causes & how to check** (most probable first):

| # | Likely cause | How to check |
|---|--------------|--------------|
| 1 | Annotation SNP identifiers use a different identity space than the reference panel | Compare a few annotation SNP identifiers against reference-panel metadata; confirm both use the same `--snp-identifier` |
| 2 | Annotation coordinates and reference-panel coordinates are on different genome builds | Check `--genome-build` and the build recorded in R² matrix or reference-panel metadata |
| 3 | Allele-aware mode was requested but annotation/reference alleles do not match | Inspect A1/A2 columns in annotation and reference metadata for a few expected overlapping SNPs |
| 4 | Chromosome labels or chromosome-sharded path tokens resolved to mismatched chromosomes | Print CHR values from the annotation and reference metadata for the failing chromosome |

**Remedies:**

1. Regenerate annotation and reference-panel artifacts with the same
   `--snp-identifier` and `--genome-build`.
2. Use an allele-aware SNP identifier only when both annotation and reference
   metadata carry matching A1/A2 columns.
3. To exclude a chromosome on purpose, leave it out of the annotation files.

A chromosome emptied by `--ref-panel-snps-file` or genotype QC is skipped with a WARNING instead; check its funnel rows if that skip was not intended.

**Developer reference:** `_refpanel.adapters._align_annotations()` during `RefPanel.prepare_chromosome()`, the gene-support pass `_apply_direct_gene_gate_b()`, and the batch check in `_ldscore.calculator._DirectProducer.reduce()` · **Exception:** `LDSCInputError`

<a id="ldscore-r2-matrix-is-incompatible"></a>

### `ldscore`: an R² matrix is incompatible

**Symptom:** `R2 parquet '...' does not match the format written by the current ldsc build-r2-matrix: ...` / `ldscore could not use R2 parquet ...` / `ldscore cannot use R2 parquet ... omits SNP pairs with adjusted R2 below ...` / `... sidecar ... missing`

**Likely causes & how to check** (most probable first):

| # | Likely cause | How to check |
|---|--------------|--------------|
| 1 | The R² matrix was written by an older LDSC3 version or another tool | Inspect columns; current files contain `IDX_1`, `IDX_2`, `R2`, and `SIGN_R`. The former `SIGN` column is no longer accepted; rebuild according to the [current format requirements](../reference/formats/r2-matrices.md#5-caveats-and-constraints). |
| 2 | The matching `chrN_meta.tsv.gz` sidecar is missing or was copied from another R² matrix | Confirm each `chrN_r2.parquet` has a same-directory `chrN_meta.tsv.gz` with matching timestamps/provenance |
| 3 | Parquet schema metadata was stripped or edited | The error lists each missing or invalid key; current R² matrices record `ldsc:r2_encoding`, `ldsc:r2_scale`, `ldsc:r2_bias=unbiased`, `ldsc:n_samples`, `ldsc:ld_window_mode`/`ldsc:ld_window_value`, `ldsc:r2_min`, `ldsc:sorted_by_build`, `ldsc:n_snps`, `ldsc:sidecar_identity_sha256`, and identity keys, with an int16 `R2` column |
| 4 | The R² matrix directory mixes chromosomes or builds from different reference-panel runs | List files in the R² matrix directory and compare recorded `genome_build` metadata |
| 5 | Duplicate or conflicting pair rows survived in the R² matrix | Regenerate from a deduplicated reference-panel source and current `ldsc build-r2-matrix` |
| 6 | The R² matrix was built with a positive `--r2-min`, so it omits low-R² pairs | The error names `ldsc:r2_min`; rebuild the R² matrix without `--r2-min` for `ldscore`, and keep the truncated R² matrix only for `ldsc query-r2-matrix` |
| 7 | The LD window is in a different unit from the R² matrix's build window, or wider than it (`windows in different units cannot be checked for coverage` / `is wider than the input R2 matrix window`) | The error names both windows; rerun with the R² matrix's unit and a value no wider, or rebuild the R² matrix with the requested window. See [LD windows for R² matrices](../reference/formats/r2-matrices.md#37-ld-windows-and-sidecar-values-in-ldscore) |

**Remedies:**

1. Regenerate the reference panel with the current `ldsc build-r2-matrix`.
2. Keep each `chrN_r2.parquet` with its matching `chrN_meta.tsv.gz` sidecar; do
   not hand-edit or independently copy sidecars.
3. Use one consistent R² matrix directory per genome build and reference-panel run.

**Developer reference:** `_refpanel.r2_format.read_r2_parquet_contracts`, `_refpanel.r2_format.SortedR2BlockReader`, sidecar-binding guards, and `validate_ldscore_window_within_r2_matrix_window` · **Exception:** `LDSCInputError`; `LDSCUsageError` for LD-window mismatches (cause 7)

### `ldscore`: genome build could not be resolved consistently

**Symptom:** `Genome-build consistency check failed for ldscore.`

**Likely causes & how to check** (most probable first):

| # | Likely cause | How to check |
|---|--------------|--------------|
| 1 | A coordinate input, including a small keep list, lacks sufficient build evidence | Read the `unresolved` group for the input flag, concrete path, and reason |
| 2 | Annotation, PLINK BIM, gene coordinate catalog, or R² matrix declarations identify different builds | Compare the complete `hg19` and `hg38` groups, including every selected shard |
| 3 | R² Parquet and sidecar declarations disagree | Keep matching `chrN_r2.parquet` and `chrN_meta.tsv.gz` artifacts together |
| 4 | CHR/POS columns were parsed from the wrong fields | Inspect the listed input's header; matching build-specific position columns take precedence over generic POS |

**Remedies:**

1. Check the listed inputs and supply coordinates on one build.
2.
For verified unlabeled inputs, assert the build with `--genome-build hg19` or `--genome-build hg38`; this bypasses fresh inference but cannot override existing declarations or convert coordinates.
3. Exact paths, globs, and `@` suites all support auto inference, including PLINK-only input and suites without chromosome 1.

See the [complete grouped report](#common-genome-build-consistency-check-failed) and [keep-list guidance](#common-keep-list-build-could-not-be-resolved).

**Developer reference:** `_ldscore.calculator._resolve_ldscore_chr_pos_genome_build()` and `_genome.build_policy.GenomeBuildCheck` · **Exception:** `LDSCInputError`

### `ldscore`: unusable CM for ld-wind-cm

`--ld-wind-cm` needs a genetic-map coordinate (`CM`) that orders SNPs along a
chromosome. The run aborts when the reference panel's `CM` is all zero, constant,
or missing — it cannot define a genetic-distance window, and `--yes-really` does
**not** override this (it only authorizes whole-chromosome windows from *valid*
`CM`).

| Likely cause | Check |
|---|---|
| PLINK `.bim` has an all-zero `CM` column (the common PLINK default) | Inspect the third `.bim` column; it is `0` for every SNP |
| Genome build for the genetic map could not be determined (rsID-family modes) | Confirm whether `--genome-build` was passed |

**Remedies:**

1. Add a genetic map for the panel's build:
   `--genetic-map-hg38-sources <file>` or `--genetic-map-hg19-sources <file>`.
   ldscore interpolates `CM` at the `.bim` positions (PLINK backend).
2. Provide a `.bim` whose third column carries real genetic-map positions, or a
   R² matrix built with `ldsc build-r2-matrix --genetic-map-<build>-sources --ld-wind-cm`.
3. Use a physical-distance window instead: `--ld-wind-kb 1000` or `--ld-wind-snps`.
   With an R² matrix, use the R² matrix's own window unit, no wider than its build window.
4. In rsID identifier modes, pass `--genome-build hg19` / `--genome-build hg38`
   so the matching genetic map can be selected.

## `build-r2-matrix`

### `build-r2-matrix`: no reference-panel artifacts were produced

**Symptom:** `Reference-panel SNP selection failed: N input issue(s)...` listing chromosomes where `all K SNPs were removed at <stage>` (SNP restriction, identity cleanup, or liftover and duplicate-position filtering), with `diagnostics/input_issues.tsv` / `build-r2-matrix stopped at chromosome C: PLINK genotype filtering removed all ...` / `build-r2-matrix produced no chromosome artifacts...` / `...retained no parquet metadata rows...`

SNP selection is checked for every chromosome before any genotypes are read, and every autosome that loses all its SNPs is reported in one error; an autosome emptied by genotype filtering stops the run when it is reached. The causes below apply per chromosome.

**Likely causes & how to check** (most probable first):

| # | Likely cause | How to check |
|---|--------------|--------------|
| 1 | The SNP restriction does not overlap the PLINK source panel | Compare the first IDs/coordinates in `--ref-panel-snps-file` against the `.bim` file |
| 2 | The restriction or PLINK panel is on the wrong genome build | Check `--source-genome-build` and any build-specific restriction position columns |
| 3 | Identity cleanup dropped all rows because identifiers are missing or duplicated | Inspect `diagnostics/dropped_snps/chr*_dropped.tsv.gz` for `reason` values |
| 4 | Liftover or duplicate-position filtering removed every retained SNP | Inspect the dropped-SNP sidecar for `unmapped_liftover`, `source_duplicate`, or `target_collision` |
| 5 | MAF or individual filtering removed all SNPs of a chromosome when its genotypes were read | Check the logged counts, `--maf-min`, and the keep-individual file against the `.fam` file |

**Remedies:**

1. Build the SNP restriction from the same PLINK source build and SNP identifier mode.
2. Use the error/log counts to identify the failing stage; inspect a named dropped-SNP sidecar when it was written and covers that stage.
3. Correct an unintended threshold, restriction, or individual selection, or regenerate compatible PLINK/reference inputs.
   If the intended filters legitimately retain no usable SNPs, do not weaken them solely to make the command finish; any revised filter changes the R² matrix being constructed.

**Developer reference:** `_refpanel.builder.R2MatrixBuilder.run()` and reference-panel metadata cleanup paths · **Exception:** `LDSCInputError`

### `build-r2-matrix`: liftover or genetic-map configuration is incomplete

**Symptom:** `build-r2-matrix cannot emit <build>...` / `...genetic map...`

**Likely causes & how to check** (most probable first):

| # | Likely cause | How to check |
|---|--------------|--------------|
| 1 | `--ld-wind-cm` was requested without the source-build genetic map | Confirm the matching `--genetic-map-hg19-sources` or `--genetic-map-hg38-sources` option is present |
| 2 | Liftover enabled an opposite-build output without that build's genetic map | Check whether a liftover chain file emits the other build |
| 3 | The genetic map is for the wrong build or lacks the failing chromosome | Inspect CHR values in the map and compare them to `--source-genome-build` |
| 4 | Genetic map shards overlap or are not sorted by chromosome/position | Sort and deduplicate map rows by CHR/POS |
| 5 | A chain file is missing, reversed, or incompatible with the source/target pair | Check the chain filename and source-target CLI flag direction |

**Remedies:**

1. For cM windows, provide every genetic map needed by the emitted build(s).
2. Use `--ld-wind-kb` or `--ld-wind-snps` when genetic maps are unavailable.
3. Use liftover chains only in chr_pos-family modes and make the chain direction
   match the source and target builds.

**Developer reference:** `_refpanel.builder.R2MatrixBuilder._prepare_build_state()`, chromosome liftover setup, and genetic-map interpolation helpers · **Exception:** `LDSCUsageError` / `LDSCInputError`

### `build-r2-matrix`: SNP restriction does not match the source panel

**Symptom:** `Genome-build consistency check failed for build-r2-matrix...`, with `--ref-panel-snps-file` listed under a conflicting build or `unresolved`.

**Likely causes & how to check** (most probable first):

| # | Likely cause | How to check |
|---|--------------|--------------|
| 1 | The restriction's generic POS column is on a different build than the PLINK source | Infer or inspect several CHR/POS rows against hg19/hg38 |
| 2 | More than one alias describes the same build-specific position | Keep one column for each build; hg19 and hg38 columns may coexist as alternatives |
| 3 | CHR/POS columns were parsed from the wrong delimiter or header names | Inspect the header with `head -1 <file> | cat -A` |
| 4 | Allele-aware restriction rows are missing A1/A2 values | Check whether every row has both allele columns when using allele-aware modes |
| 5 | The file contains too little HM3 overlap for automatic build inference | Verify its build, then specify `--source-genome-build hg19` or `hg38` |

**Remedies:**

1. Prepare the restriction file on the same build as the PLINK source panel.
2. A matching `hg19_POS` or `hg38_POS` is preferred over generic `POS` without re-inference.
   Otherwise generic positions need enough evidence in auto mode.
   A concrete `--source-genome-build` bypasses fresh inference and must agree with declarations; it does not convert positions.
3. Fix ragged rows or delimiter/header issues before rerunning.

**Developer reference:** `_refpanel.builder.R2MatrixBuilder._resolve_source_genome_build()`, `_genome.build_policy.GenomeBuildCheck`, and `_core.columns.resolve_restriction_chr_pos_columns()` · **Exception:** `LDSCInputError`

### `build-r2-matrix`: reference-panel artifact is incompatible

**Symptom:** `Reference-panel metadata sidecar is missing...` / `...required identity/provenance keys are missing...` / `...has missing, non-numeric, or out-of-range MAF for ... SNP rows...`

**Likely causes & how to check** (most probable first):

| # | Likely cause | How to check |
|---|--------------|--------------|
| 1 | The R² Parquet file was copied without its matching `chrN_meta.tsv.gz` sidecar | Confirm each `chrN_r2.parquet` has a same-directory `chrN_meta.tsv.gz` |
| 2 | The artifact was written by an older LDSC3 version | Inspect parquet or sidecar metadata for `artifact_type`, `snp_identifier`, and `genome_build` |
| 3 | A parent directory with multiple build children was passed without a concrete build | Check for `hg19/` and `hg38/` children under `--r2-matrix-dir` |
| 4 | The metadata sidecar lacks required SNP identity columns for the active mode | Inspect sidecar columns for SNP or CHR/POS plus A1/A2 in allele-aware modes |
| 5 | The selected chromosome/build was not generated | List `chr*_r2.parquet` files in the selected build directory |
| 6 | The sidecar was edited or produced by another tool, leaving `MAF` missing, non-numeric, or outside [0, 0.5] | The error names the count and first SNPs; `build-r2-matrix` writes the minor-allele (`A1`) frequency for every row |

**Remedies:**

1. Keep R² Parquet files and their sidecars together as one artifact family.
2. Regenerate the reference panel with the current `ldsc build-r2-matrix`.
3. Pass a concrete build-specific R² matrix directory or set the matching genome build.

**Developer reference:** `_refpanel.adapters.R2MatrixRefPanel` and metadata sidecar readers · **Exception:** `LDSCInputError`

## `convert-legacy-ldscores`

The [conversion user guide](workflows/migrating-from-legacy-ldsc.md) gives full layouts, table requirements, flags, and copy-and-rename commands. After logging starts, conversion exceptions and their remedies are recorded in `diagnostics/convert-legacy-ldscores.log` even at `--log-level ERROR`, with available issues in `diagnostics/conversion_issues.tsv.gz`. Invalid directory arguments or output collisions may fail before a log opens.

### standard filenames and directory organization

Use these patterns directly inside the selected input directories. `<chrom>` means unpadded `1` through `22`; `<prefix>` is constant within each family and can be empty. Include separators in the prefix, such as `baseline.`. `(.gz)` denotes a plain or gzip alternative; counts must be plain text.

| Input | Standard pattern | Example | Repair practice |
| --- | --- | --- | --- |
| Reference or weight scores | `<prefix><chrom>.l2.ldscore(.gz)` | `weights.1.l2.ldscore.gz` | Restore missing chromosomes or rename valid copies consistently; select the directory containing the files directly |
| Baseline annotations | `<reference-prefix><chrom>.annot(.gz)` | `baseline.1.annot.gz` | Use the same prefix as the reference scores and restore all 22 matching shards |
| Required common-SNP counts | `<reference-prefix><chrom>.l2.M_5_50` | `baseline.1.l2.M_5_50` | Restore the matching source count or rename a misnamed copy; `.l2.M` cannot substitute |
| Optional all-SNP counts | `<reference-prefix><chrom>.l2.M` | `baseline.1.l2.M` | Supply matching original counts if available; see the profile-specific missing-count policy |
| Baseline frequencies | `<prefix><chrom>.frq(.gz)` | `1000G.EUR.QC.1.frq.gz` | Rename valid `.freq` copies to `.frq` if needed, restore missing chromosomes, and select one complete family |

Place different releases in separate directories when discovery finds multiple complete families. Do not point at a parent directory, filename prefix, individual file, or glob. Counts and annotations must match the selected reference prefix; weight and frequency prefixes can differ. Renaming repairs discovery only: retain the source release and valid table contents, and decompress gzip contents before changing a compressed count file to its plain-text standard name.

### unsupported weight filenames

**Symptom:** family discovery reports unsupported `.w.l2.ldscore(.gz)` filenames, such as `1.w.l2.ldscore.gz`.

Both reference and weight suites require `<prefix><chrom>.l2.ldscore(.gz)` names across chromosomes 1–22. Rename copies in a separate input directory, for example `1.w.l2.ldscore.gz -> weights.1.l2.ldscore.gz`, and pass that directory to `--legacy-weight-dir`. Use the same correction without `.gz` for plain files. Each affected file in the inspected directory is recorded in `diagnostics/conversion_issues.tsv.gz` as `discarded_unsupported_weight_filename`. If a complete recognized family exists alongside these files, conversion uses that family and records the unused alternatives with a warning and in `legacy_ldscore_conversion.ignored_files`. After a failed conversion, rerun with a new output directory or `--overwrite` to replace its diagnostics.

### conversion rejected the legacy suite

Inspect `diagnostics/conversion_issues.tsv.gz`; conversion failures intentionally
leave diagnostics but no canonical `metadata.json` or Parquet result.

Common causes are a missing or ambiguous chromosome 1-22 family, duplicate
rsIDs within or across shards, non-finite LD-score/annotation/frequency values,
a partial or thin annotation family, missing `.l2.M_5_50`, non-bijective
annotation-to-score names, or baseline `.M`/`.M_5_50` values that disagree with
the full annotation/frequency reconstruction. The converter reports causal SNPs
or annotation names in the exception and full available issue rows in the
sidecar. Do not mix releases or edit the source suite in place; supply coherent
reference, weight, and (for baseline) frequency directories and rerun with
`--overwrite` because the first failure already owns its diagnostic files.

For `chr_pos`, `--genome-build auto` must find decisive reference evidence. Use
an explicit build only when it is known; an explicit declaration that conflicts
with decisive evidence is rejected. See
`docs/reference/formats/legacy-ldscores.md` for the complete contract.

**Developer reference:** `_ldscore.legacy_converter.LegacyLDScoreConverter` · **Exception:** `LDSCInputError`

## regression

### regression: no SNPs remain after merging regression inputs

**Symptom:** `h2 regression retained no overlapping...` / `rg regression retained no overlapping...`

**Likely causes & how to check** (most probable first):

| # | Likely cause | How to check |
|---|--------------|--------------|
| 1 | Summary statistics and LD-score directories use different SNP identifier modes | Compare `snp_identifier` in each artifact's metadata |
| 2 | Coordinate-family inputs are on different genome builds | Compare `genome_build` metadata for the summary statistics and the LD-score directory |
| 3 | One input was generated allele-aware and another was generated allele-blind | Check whether each artifact has A1/A2 columns and allele-aware provenance |
| 4 | The regression SNPs or filters removed all shared SNPs | Inspect the retained SNP counts and any dropped-SNP sidecars from upstream commands |
| 5 | For `rg`, the two traits use incompatible allele conventions | Confirm both munged summary statistics were generated with harmonized A1/A2 columns |

**Remedies:**

1. Regenerate the munged summary statistics and LD scores with compatible `--snp-identifier` settings.
   For coordinate-based identity, the munged `--output-genome-build` must match the LD-score build selected with `--genome-build`; the raw-input `--source-genome-build` can differ when liftover was requested.
2. Use `--allow-identity-downgrade` only for same-family allele-aware/base mixes
   when duplicate base identities have been cleaned.
3. Check whether intended SNP selection leaves a usable intersection before changing it.
   Rebuild mismatched inputs; broaden an intended SNP selection only if that change is justified for the analysis.

**Developer reference:** `_regression.runner.RegressionRunner.build_dataset()` and `_regression.runner.RegressionRunner.build_rg_dataset()` · **Exception:** `LDSCInputError`

### `partitioned-h2`: a query regression fails

**Symptom:** one or more query models could not be fitted or summarized, or their jackknife uncertainty could not be estimated.

Inspect `diagnostics/query_status.tsv` for attempted-query outcomes when it was written, and the saved `diagnostics/partitioned-h2.log` for the affected query and failed stage.
For a failed jackknife calculation, removing one SNP block can leave too little information to estimate all coefficients even when the full-data model can be fitted.
The log retains the failed block, SNP support counts, numerical ranks, and traceback for deeper investigation.

| Cause | Interpretation and action |
| --- | --- |
| Query LD scores occur in one jackknife block | Removing that block erases the query coefficient's information. Mark the query unestimable under this design; a full-data fit cannot provide the missing jackknife uncertainty. |
| Remaining columns become dependent | Deletion can make baseline/query columns dependent even when every column remains nonzero. Inspect the conditional model and the failed block. |
| Other query preparation, estimator, or summary exception | Use the logged stage, exception type, and full traceback to investigate the specific error. |

Default strict mode collects per-query failures and publishes no new scientific results if any fail.
If retaining only the successful query models is acceptable for the intended analysis, explicitly select `--continue-on-query-error`; failed queries remain in the diagnostic ledger and are not valid fitted results.
Shared input/output errors, process interruption, and all-query-failed scans still fail.
Do not substitute a pseudoinverse, fabricate uncertainty, or automatically move block boundaries.
See [the result contract](../reference/results/partitioned-h2.md#query-failures-and-continuation).

**Developer reference:** `RegressionRunner.estimate_partitioned_h2_batch()` and `LstsqJackknifeFast.block_values_to_delete_values()`.
Singular deletions raise `JackknifeIdentifiabilityError`, a NumPy `LinAlgError` subclass; other per-query errors retain their original exception type.

### `partitioned-h2`: a query worker exits or cannot write

**Symptom:** query fitting stopped because a worker exited, communication failed, or a temporary result could not be written.

This policy applies at every worker count, including `--workers 1`.
Ordinary model exceptions follow the query error policy above; worker, communication, and output failures abort the run even with `--continue-on-query-error`.
Inspect the main log for the affected query or in-flight queries and any available worker traceback.
An in-flight query identifies interrupted work, not necessarily the query that caused process termination.
The final query-status table may not have been written.

For a filesystem exception, inspect the recorded path and resolve disk-space or access problems before retrying.
For an unexplained worker exit, check scheduler/job termination messages and memory-limit evidence rather than assume that missing traceback means out-of-memory.
Reducing `--workers` reduces simultaneous model workspaces; reducing `--query-batch-size` reduces loaded query values.
`--workers 1` still uses one worker, not in-process execution.
Python scripts must guard workflow calls with `if __name__ == "__main__":`, including one-worker calls.
The parent reaps workers before removing owned scratch during handled failure; forced parent termination may leave temporary files.
The scan has no automatic resume.
See [CPU and memory policy](../reference/regression-options.md#43-query-workers-and-memory).

### `partitioned-h2`: missing overlap matrix

**Symptom:** `partitioned-h2 needs the annotation overlap matrix, but the LD-score directory has no ldscore.overlap.parquet...`

Overlap-aware partitioned heritability (both the functional and cell-type
regimes) needs `ldscore.overlap.parquet`, written by `ldsc ldscore` alongside
the baseline/query parquet files whenever the run has **two or more** annotation
columns. A baseline-only directory (multiple baseline annotations, no query) is
**not** an error — it runs the functional regime. An unpartitioned,
single-annotation run (e.g. the synthetic `base`) deliberately omits the sidecar
because its overlap collapses to a SNP count already in `metadata.json`, and
such a directory cannot be partitioned.

**Likely causes & how to check** (most probable first):

| # | Likely cause | How to check |
|---|--------------|--------------|
| 1 | The directory is an unpartitioned single-annotation run (e.g. base-only) | Inspect `metadata.json`: `baseline_columns` + `query_columns` total fewer than 2 |
| 2 | The LD-score directory predates the overlap sidecar (older `ldsc ldscore`) | Inspect `metadata.json` for a missing `files.overlap` / `overlap_config` despite >=2 annotation columns |
| 3 | The directory was copied without `ldscore.overlap.parquet` | Confirm the file exists next to `ldscore.baseline.parquet` |

**Remedies:**

1. Partitioned heritability requires >=2 annotation columns; regenerate the LD-score
   directory with explicit baseline (and optional query) annotations, not the
   synthetic `base`.
2. If the directory predates the sidecar, regenerate it with the current
   `ldsc ldscore`.
3. Keep `metadata.json`, the baseline/query parquet, and `ldscore.overlap.parquet`
   from the same run together.
4. `h2` and `rg` do not require the overlap matrix (their collinearity
   check, for `h2` and for each `rg` trait pair, uses it only when present),
   so they run on unpartitioned or older directories.

**Developer reference:** `_regression.runner.RegressionRunner.estimate_partitioned_h2_batch()` / `summarize_partitioned_h2()` · **Exception:** `LDSCInputError`

## `quantile-h2`

### `quantile-h2`: shared fitted-model metadata is missing

**Symptom:** loading a fitted query reports missing shared metadata or an invalid `global_metadata` reference.

Current query results store run-wide provenance once in the aggregate result's `diagnostics/metadata.json`.
Each query's `metadata.json` points to it with `"global_metadata": "../../metadata.json"`.
Copying only `diagnostics/query_annotations/<query>/` breaks that reference.

Restore the complete result directory, or restore the enclosing `diagnostics/metadata.json` and query metadata from the same fitted run while preserving their relative locations.
Then select the successful query folder within that tree as `--partitioned-h2-result-dir`.
Use the current LDSC3 reader, which resolves shared provenance and gives fit-specific fields precedence.
Restoring missing metadata from the original run does not require recomputing LD scores or refitting the model.
Do not substitute metadata from another run or reconstruct provenance from unrelated results.
See [keeping fitted results together](workflows/partitioned-heritability.md#reading-the-outputs) and the [output contract](../reference/results/partitioned-h2.md).

**Developer reference:** `load_fitted_partitioned_model()` in [_regression/quantile_h2.py](../../src/ldsc/_regression/quantile_h2.py) · **Exception:** `LDSCInputError`

### `quantile-h2`: fitted model or resupplied sources do not verify

**Symptom:** `quantile-h2` cannot verify that the saved fit, annotation inputs, and reference metadata describe the same fitted model.

**Likely causes & how to check** (most probable first):

| # | Likely cause | How to check |
|---|---|---|
| 1 | An aggregate cell-type root was supplied | Select one `diagnostics/query_annotations/<query>/` directory; aggregate rows are separate fits. |
| 2 | The fitted result predates coefficient-delete persistence | Check for `coefficient_delete_values.parquet`; rerun `partitioned-h2` with current LDSC3. |
| 3 | Not every original fitted annotation was resupplied | Compare `retained_ld_columns` in model metadata with source headers. |
| 4 | Reference metadata or annotations come from a different panel/release | Inspect the mismatch text for SNP counts, annotation sums, or overlaps. |
| 5 | Effective SNP identities are duplicated or target coverage is incomplete | Inspect `diagnostics/snp_alignment_issues.tsv.gz`. |

Resupply the exact original annotation sources and matching reference metadata.
R² matrix users can use the existing `chr*_meta.tsv.gz` sidecars; PLINK users should regenerate LD scores with `--export-ref-metadata`.
Do not concatenate coefficient tables from different query directories.

`quantile-h2` requires an inclusive common-SNP MAF threshold (`>=`).
If the saved fit records a different comparison, use a compatible fitted result.
Changing this SNP-count policy requires regenerating LD scores and refitting `partitioned-h2`; it changes the fitted SNP population and must be intentional.
Keep the saved comparison metadata faithful to the fit that produced it.

**Developer reference:** `_regression.quantile_h2.load_fitted_partitioned_model()` / `_regression.quantile_inputs.prepare_quantile_statistics()` · **Exception:** `LDSCInputError`

### `quantile-h2`: target values are invalid or quantiles are empty

**Symptom:** the target annotation contains an invalid value, or the requested quantile split leaves an empty group.

If the target table cannot be parsed, check its text encoding, whitespace-delimited header, chromosome labels, and numeric base-pair positions.
Confirm that `--target-annot-sources` selects the intended annotation table and `--target-annot-column` names its intended value column.
Restore or correct the source table before retrying; changing its filename does not change its coordinates or values.

By default, every target value must be numeric and finite.
If exactly one token denotes missingness, pass it with `--target-missing-value`; zero is retained unless explicitly selected.
If boundary ties produce an empty quantile, reduce `--num-quantiles` or use a less discrete target, recognizing that this changes the requested grouping.
Ties are intentionally kept in the lower-valued quantile for legacy compatibility.

**Developer reference:** `_regression.quantile_storage.target_chunks()` / `_regression.quantile_inputs.prepare_quantile_statistics()` · **Exception:** `LDSCInputError`

## `plot`

### `plot`: result directory is unsupported or incomplete

**Symptom:** `plot` cannot recognize the result directory or find the saved artifacts required for its plots.

| # | Likely cause | How to check |
|---|---|---|
| 1 | An internal per-query `partitioned-h2` directory was supplied | Pass the aggregate `partitioned-h2` result directory containing `diagnostics/metadata.json`. |
| 2 | The result predates a required current artifact | For `h2`, check for `diagnostics/ld_score_regression_bins.tsv`; rerun `ldsc h2` with the current package. |
| 3 | Metadata and result tables came from different runs | Inspect `diagnostics/metadata.json` and verify every relevant `files` entry exists relative to the same root. |
| 4 | The artifact type or scientific regime has no approved plot | Compare the input with the supported table in [the plotting manual](workflows/plotting.md#which-result-directory-produces-which-plot). |

The plot command does not infer from loose TSV files or reconstruct missing
diagnostics. Matplotlib is a required package dependency. If an existing
environment reports that it is unavailable, repair that environment with
`python -m pip install "matplotlib>=3.9,<4"` and rerun the command.

**Developer reference:** `_plotting.dispatch.plot_result()` · **Exception:** `LDSCInputError`

### `plot`: `rg` heritability source is invalid

**Symptom:** `plot` cannot read or validate the declared per-trait heritability table for an `rg` result.

| Likely cause | How to check and fix |
|---|---|
| The declared heritability table is unreadable or malformed | Check `files.h2_per_trait` in `diagnostics/metadata.json`. Restore the original readable TSV from the same `rg` run, with consistent field counts and valid text encoding. |
| Required columns or trait identifiers are missing | The table must contain `trait_name`, `total_h2_obs`, and `total_h2_obs_se`, and each row must have a nonempty trait name. Restore the canonical table. |
| A trait has multiple rows | Check repeated `trait_name` values and restore one single-trait fit per trait; do not substitute pair-specific fits. |
| The metadata path is absolute, escapes the result root, or points to a directory | Use a relative path to a regular file inside the same result directory. Symlinks must also resolve inside that root. |

Missing declarations/files or absent trait rows are allowed and display `failed`. Missing, nonfinite, or nonnumeric heritability or SE values and negative SEs also display `failed`; these do not abort plotting. Restore valid saved inputs and rerun `ldsc plot --result-dir RESULT_DIR --overwrite` to refresh existing figures. See [heritability annotation behavior](../dev/design/plotting.md#heritability-annotations-in-rg-plots).

**Developer reference:** `_plotting.dispatch.plot_result()` and `_plotting.builders._trait_h2_labels()` · **Exception:** `LDSCInputError`
