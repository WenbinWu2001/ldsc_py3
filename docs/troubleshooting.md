# Troubleshooting

Last updated on: 2026-09-14

This reference explains `ldsc` errors that can **abort a run** and have more than
one likely cause. It is organized by command. Each entry lists the likely causes
(ranked most-probable first), how to confirm each, and how to fix it.

Most errors are self-explanatory from their terminal message alone and are not
repeated here. When a message says `... see docs/troubleshooting.md#<section>`,
that slug is a heading below — jump to it.

## Staged input validation

**Raised by:** `_input_preflight.InputGate` and workflow content gates. **Symptom:** a combined error lists multiple input objects or chromosomes before staging or numerical work proceeds.

| Cause | Check and repair |
| --- | --- |
| Missing path or PLINK companion | Read every `source` and `chrom` in `input_issues.tsv`; restore each selected BED/BIM/FAM trio. Quote globs and preserve dotted prefixes. |
| Wrong header, delimiter, or compressed/Parquet file | Inspect the named file's header and format; restore canonical fields and readable compression. Annotation checks use the existing alias and identity rules. |
| Missing or incompatible index component | Restore the complete immutable index with matching chromosome/root metadata; never combine components from separate indexes. |
| Declared chromosome coverage differs from content | Use matching annotation/reference suites and a covered catalog selection. An explicit `@` requires the workflow's documented chromosome suite; a glob selects actual matches. |
| Later alignment, identity, or support failure | These checks need scanned rows or filtered intersections. Use the existing chromosome-scope, gene-list, and SNP-drop audits; passing the earlier path gate cannot establish them. |

Writing workflows retain the six fields `input_role`, `source`, `chrom`, `reason`, `details`, and `repair` in `diagnostics/input_issues.tsv`; gene-index construction uses `.<index>.build-state/input_issues.tsv`. Python exceptions also expose the repair frame as `input_issues`. Existing specialized diagnostics remain authoritative for gene rows, PLINK contents, and alignment. Correct all reported defects at the current gate before rerunning; later content gates may expose defects that could not be established from headers.

Output collisions are a separate error: passing `--overwrite` authorizes owned output replacement but never bypasses input validation. Early input failures preserve existing scientific files and retain the established authorized-overwrite `RUN_FAILED.txt` contract. INFO phase records distinguish validation, staging, computation, and publication; periodic progress occurs when bounded chunks return, not during a blocking third-party call.

## Common

### Common: `RUN_FAILED` is present in an output directory

The latest authorized overwrite failed. Read the marker first, then inspect the
detailed log it names when a log was opened. The directory may contain
incomplete or mixed artifacts because LDSC preserves each workflow's existing
write order and performs no rollback or restoration. Correct the underlying
error and rerun the same materializing command with `--overwrite`; a successful
retry removes the applicable marker. Do not treat the marker as a scientific
result file.

Markers, logs, audits, and scientific outputs use the same expanded destination. A marker alone does not prevent retry. If an older version placed logs or markers under a literal `$VARIABLE` or `~` directory, update the package and rerun using the intended destination; existing misplaced directories are not migrated or deleted automatically. The [output-directory audit](audits/2026-09-14-output-directory-retries.md) records the affected paths and regression coverage.

### Common: input path did not resolve to one file

**Raised by:** `path_resolution.resolve_scalar_path()` and group/path-prefix resolvers
· **Exception:** `LDSCInputError`
**Symptom:** `Could not resolve <label> path from token '<token>': matched 0 files` / `matched N files`

**Likely causes & how to check** (most probable first):

| # | Likely cause | How to check |
|---|--------------|--------------|
| 1 | The path is misspelled or relative to a different working directory | `pwd`; then `ls <path>` from the same shell |
| 2 | A glob is too broad for an input that must be one file | `python -c "import glob; print(glob.glob('<pattern>'))"` |
| 3 | An ordinary file input uses an unsupported bare prefix (PLINK prefixes are supported) | Check whether the token looks like `chr@` / `.@.` rather than a bare prefix |
| 4 | The file exists only with a suffix LDSC does not infer | `ls <prefix>*`; pass the full filename including suffix |
| 5 | A PLINK prefix is incomplete | Confirm all three files exist: `<prefix>.bed`, `<prefix>.bim`, `<prefix>.fam` |

**Remedies:**

1. Pass the exact existing file path when the command expects a single file.
2. Narrow broad globs so they match exactly the intended file.
3. For annotation chromosome suites, use an explicit `@` token such as `baseline.@.annot.gz` or a quoted glob with the full file suffix.

**PLINK resolution (`ldscore`, `build-gene-ldscore-index`, `build-r2-panel`):** A plain stem such as `1000G.EUR.QC.` is supported and preserves numeric suffixes such as `.22`. These commands share `path_resolution.inspect_plink_inputs()`. If `Could not resolve PLINK inputs` appears, read every listed issue: restore each missing BED/BIM/FAM member, replace malformed or truncated trios, verify BIM chromosome contents against any explicit `@` declaration, and narrow inputs when multiple trios contain the same chromosome. Renaming a file does not change its chromosome assignment. Gene-index construction requires chromosomes 1-22 even without `@`; direct query LD-score full-suite declarations also retain their coverage checks. Direct LD-score issues appear in its input diagnostics; gene-index issues are retained in `.<output-name>.build-state/plink_input_issues.tsv`; R²-builder issues appear in `diagnostics/plink_input_issues.tsv`. See the [shared resolution contract](current/path-specification.md#plink-prefix-resolution).

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
| 1 | The artifact was written by an older LDSC package version | Inspect the artifact metadata sidecar or parquet footer for `artifact_type`, `snp_identifier`, and `genome_build` |
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
| 2 | The column exists under a synonym you must declare | Re-run with an explicit hint, e.g. `--snp MarkerName --a1 Allele1` |
| 3 | Wrong delimiter, so the whole header parsed as one column | `head -1 <file> \| cat -A` — look for one field with embedded tabs/commas |
| 4 | Wrong `--input-format`, so expected columns differ | Confirm `--input-format` matches the file; `auto` and an explicit profile share aliases, optional-field rules, and validation |

**Remedies:**

1. Pass explicit column hints for the unmapped field(s): `--<field>-col <name>`.
2. Rename the columns to recognized names, or add the alias to `column_inference.py`
   if broadly useful.
3. Verify the delimiter is whitespace/tab as expected for `.sumstats`/`.txt` inputs.

### munge-sumstats: multiple or conflicting sample-size column strategies

**Raised by:** `_kernel/sumstats_munger` sample-size column resolution
· **Exception:** `LDSCInputError` or `LDSCUsageError`
**Symptom:** `munge-sumstats found multiple sample-size strategies...`, an
incomplete `--N-cas-col`/`--N-con-col` pair, or a conflict with `--N-col`.

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

### munge-sumstats: no SNPs remain after filtering

**Raised by:** `_kernel/sumstats_munger` (post-filter and keep-list paths)
· **Exception:** `LDSCInputError`
**Symptom:** `munge-sumstats removed every SNP...` / `...no SNPs remain after SNP keep-list restriction`

Munging now restricts to packaged HM3 by default. To process all input SNPs subject to ordinary QC, use `--no-snp-restriction`; to replace HM3 with a custom list, use `--sumstats-snps-file FILE`. For cross-build output, both override modes require `--liftover-chain-file`; default HM3 uses automatic quick liftover unless a chain is supplied. Source-build inference failure requires an explicit `--source-genome-build`; output builds always require an explicit choice.

**Likely causes & how to check** (most probable first):

| # | Likely cause | How to check |
|---|--------------|--------------|
| 1 | The input has no overlap with the default HM3 set, or a custom keep-list uses a different SNP-id space (rsID vs chr:pos) | Compare the first IDs of each file; confirm both use `snp_identifier=<mode>` |
| 2 | Keep-list and sumstats are on different genome builds | Compare the `genome_build` in the error line against the keep-list build |
| 3 | INFO / MAF / N thresholds removed every row | Inspect filter counts in the error/log and verify the column meanings and intended thresholds; earlier QC removals are not in the dropped-SNP sidecar |
| 4 | Input is effectively empty after a delimiter/format mis-parse | `zcat <file> \| wc -l` |

**Remedies:**

1. Regenerate the keep-list in the same `snp_identifier` mode and build as the sumstats.
2. Relax quality filters and re-run; review the dropped-SNP audit sidecar to see
   which filter removed the rows.

### munge-sumstats: curated artifact is malformed or outdated

**Raised by:** `sumstats_munger.load_sumstats()` and its metadata helpers
· **Exception:** `LDSCInputError`
**Symptom:** `Cannot load curated sumstats at '<path>': ...` (old schema /
bad provenance / missing A1-A2 / duplicate identity rows)

**Likely causes & how to check** (most probable first):

| # | Likely cause | How to check |
|---|--------------|--------------|
| 1 | Artifact predates the current self-describing parquet schema | Inspect the parquet footer for `ldsc:artifact_type`, `ldsc:snp_identifier`, and `ldsc:genome_build` |
| 2 | The artifact is footerless Parquet | Re-munge from the raw GWAS input to produce self-describing `sumstats.parquet`; footerless Parquet is not an LDSC2 compatibility format |
| 3 | Identity provenance in the parquet footer is invalid/corrupt | Inspect the footer identity fields against the current contract |
| 4 | An allele-aware `snp_identifier` artifact lacks A1/A2 columns | `python -c "import pandas; print(pandas.read_parquet('<f>').columns)"` |
| 5 | Duplicate/invalid SNP-identity rows survived in the artifact | Re-munge from raw input; the loader reports the dropped-row reasons |

**Remedies:**

1. Re-run `ldsc munge-sumstats` from the **raw** GWAS file to regenerate the
   artifact with the current schema.
2. A genuine LDSC2 `.sumstats` or `.sumstats.gz` text artifact is accepted
   directly by regression when it contains `SNP`, `A1`, `A2`, `Z`, and `N`;
   see `docs/current/legacy-sumstats-compatibility.md`.
3. Do not hand-edit curated `.sumstats`/`.parquet` artifacts; treat them as outputs.

### munge-sumstats: liftover dropped all rows

**Raised by:** `_kernel.liftover.apply_sumstats_liftover()`
· **Exception:** `LDSCInputError`
**Symptom:** `Summary-statistics liftover from <source> to <target> using <method> dropped all rows...`

**Likely causes & how to check** (most probable first):

| # | Likely cause | How to check |
|---|--------------|--------------|
| 1 | The declared source build does not match the input CHR/POS coordinates | Compare several input positions against the declared `--source-genome-build` |
| 2 | CHR/POS columns are missing, malformed, or parsed from the wrong fields | Inspect the count-bearing error and the first rows of the raw file |
| 3 | The selected mapping cannot resolve the retained source coordinates | Check the selected method and source/output build pair; verify the chain direction if overriding packaged HM3 mapping |
| 4 | Source or target coordinates collide after liftover | Inspect duplicate-source and duplicate-target counts in the error/log |
| 5 | The chain file is for the wrong direction or assembly pair | Confirm the chain filename/source-target pair matches the command flags |

**Remedies:**

1. Fix the source CHR/POS coordinates or pass the correct `--source-genome-build`.
2. For broader SNP coverage, set `--no-snp-restriction` or a custom `--sumstats-snps-file`, and supply a source-to-output chain for cross-build conversion. Supplying a chain alone changes the mapping method while retaining the default HM3 restriction.
3. Review the error/log counts to identify whether missing coordinates, unmapped variants, or duplicate coordinates removed the rows. After a successful run with partial removals, use the dropped-SNP sidecar to inspect the affected rows.

Successful runs report the method and counts in stdout and `diagnostics/sumstats.log`, including at `--log-level ERROR`; failed runs retain the error and available counts in the log and do not print a success summary. An all-dropped liftover fails before writing the dropped-SNP sidecar, so inspect its count-bearing error first. Mapping input is measured after earlier QC and keep-list filtering. See the [liftover behavior table](current/munge-sumstats.md#liftover-rules).

## annotate

### annotate: input preflight

Missing files, invalid annotation values or headers, and incorrect chromosome members are recorded in `diagnostics/input_issues.tsv`. Exact paths and globs select their actual inputs. For gene-list annotation, `@` explicitly requires autosomes 1–22; each member must contain its declared chromosome. Repair all listed inputs before rerunning. Input defects are not interpreted as zero pathway support.

### annotate: gene-list preflight

Use exactly one BED or gene-list query route. Gene lists require a readable coordinate catalog and explicit nonnegative `--padding-bp`, including `0` for gene bodies. Check `diagnostics/gene_catalog_issues.tsv`, `diagnostics/gene_list_audit.tsv.gz`, and `diagnostics/gene_list_resolution_summary.tsv` for catalog defects, malformed or unreadable lists, naming collisions, and unresolved identifiers. `resolved-only` permits only the established identifier-resolution omission allowlist; it does not bypass malformed rows, missing sources, or incomplete chromosome coverage.

After resolution, every selected gene must be covered by validated baseline contents. `diagnostics/chromosome_scope.json` and the gene audit explain missing coverage. Measured support is then reported as `annotation_snp_count` on cleaned baseline SNP rows; reference-panel support remains unevaluated. Empty and globally unsupported focal queries are skipped, usable siblings continue, and an all-skipped batch fails without publishing a new query family. All-one annotation columns remain valid.

An authorized overwrite can leave `RUN_FAILED.txt` after failure. Inspect diagnostics and rerun with corrected inputs and `--overwrite`; no rollback to earlier results is promised. Private staging is removed on handled failure or bundle closure. Python callers should use `with run_annotate(...):` or call `close()` after consuming the returned handle; canonical `query.<chrom>.annot.gz` files remain available.

### annotate: no annotation SNP rows remain

**Raised by:** `annotation_builder.AnnotationBuilder._run_single_universe()`,
`annotation_builder.AnnotationBuilder._run_sharded_inputs()`, and
`annotation_builder.AnnotationBuilder._apply_identity_cleanup()` · **Exception:** `LDSCInputError`
**Symptom:** `annotate loaded no SNP rows...` / `...no annotation rows remain after SNP identity cleanup...`

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
2. Regenerate annotation files with `CHR`, `POS`, `SNP`, `CM`, and at least one annotation column.
3. Match `--snp-identifier` and `--genome-build` to the annotation identity columns, then rerun.

## ldscore

### build-gene-ldscore-index: output directory is nonempty but invalid

A regular root `RUN_FAILED.txt` file is an owned diagnostic, not an index artifact. If a failed build leaves only that marker (optionally with recognized legacy build diagnostics), rerun the corrected command in the same output directory. Preflight preserves the marker, and a successful build removes it. Versions before this retry fix incorrectly rejected marker-only directories even with `--overwrite`; update the package rather than deleting the directory.

If the current version still reports this error, the directory contains an incomplete index or other unrecognized contents. The marker does not grant permission to replace those contents. Preserve them for inspection and choose a fresh output directory, or explicitly resolve the conflicting files. A valid existing index still requires `--overwrite`.

### build-gene-ldscore-index: baseline/PLINK identifier intersection fails

**Symptom:** the builder reports an empty `rsid` or `chr_pos` baseline/PLINK
intersection before genotype QC, or warns that duplicate groups were dropped.

The builder uses the same identifier-key inner intersection as direct
PLINK-backed `ldscore`. Baseline-only and PLINK-only rows are allowed, dropped,
and counted. For a duplicate effective identity, every row in that source group
is dropped; no representative is selected. Inspect
`<index-dir>/diagnostics/dropped_snps/chrN_dropped.tsv.gz` after success (or the
live log during a failing build). An empty result cannot define an LD universe.
Under rsID matching, coordinate disagreements warn and use PLINK coordinates.
Under coordinate matching, differing SNP labels warn and the PLINK label is
published. During
a failed build, check
`<parent>/.<index-name>.build-state/build-gene-ldscore-index.log`; after a
successful build, check `<index-dir>/diagnostics/build-gene-ldscore-index.log`.
Confirm that all sources are hg19, that `--snp-identifier` matches the intended
key, and that the inputs overlap after complete duplicate groups are removed.
There is no build inference or liftover; `--genome-build hg19` is an advanced-user
provenance assertion. A failed first build leaves the index destination absent or
empty and can be retried directly; the previous log is archived under
`<parent>/.<index-name>.build-state/history/`.

### build-gene-ldscore-index: a partial hidden stage remains after failure

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

### build-gene-ldscore-index: transaction cleanup remains after success

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

The preceding release used a visible `<index-dir>.build/` log directory and a
standalone hidden lock. On the next invocation, recognized logs/history migrate
to `.<index-name>.build-state/`, and the obsolete standalone lock participates
in that transition run before it is removed. Unrecognized files are preserved
and reported rather than deleted.

### ldscore: an explicit gene index is missing, corrupt, or incompatible

**Symptom:** indexed mode rejects metadata IDs, mode/build agreement, the ordered published-row metadata digest, chromosome coverage, duplicate effective identities, Parquet rows, NPZ members, CSR structure/dtypes, or dimensions.

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

The retired `effective_identity_sha256` field alone does not require rebuilding: otherwise-valid existing indexes still load. New indexes omit this redundant field and require an updated reader; upgrade the reader if an older installation reports a missing effective-identity digest. The retained `published_row_metadata_sha256`, semantic `index_id`, and input fingerprints remain required under their existing contracts.

### build-gene-ldscore-index: required identity or build option is missing

**Symptom:** parsing stops with `--snp-identifier is required; choose rsid or
chr_pos.` or `--genome-build is required; choose hg19.`

New construction has no default or inference for either decision. Pass exactly
one of `--snp-identifier rsid` or `--snp-identifier chr_pos`, plus
`--genome-build hg19`. Do not use `auto`, hg38, or an allele-aware mode. These
options belong only to Stage 1 construction. Remove both options from Stage 2
`ldsc ldscore --gene-ldscore-index-dir ...`, which inherits them from the index.

### ldscore: gene-list preflight stopped or a focal query is missing

**Raised by:** query annotation status handling · **Symptom:** one requested
query is absent from `ldscore.query.parquet` or downstream partitioned-h2 rows,
or the run reports that every query was skipped.

Strict gene-list resolution stops before LD-score work if any focal/control row
cannot resolve uniquely. Start with
`diagnostics/gene_list_resolution_summary.tsv`, then filter
`diagnostics/gene_list_audit.tsv.gz` by the affected `source` and
`disposition == 'rejected'`. After Gate A succeeds,
`diagnostics/query_annotation_status.tsv` distinguishes empty/zero-resolved,
zero-annotation-SNP, and zero-variance focal queries. Usable sibling queries
continue; a requested control never disappears silently. Correct the list or
catalog evidence and rerun with `--overwrite` because diagnostics are owned
artifacts. When all focal queries are skipped, no root scientific result is
written.

For the complete reason vocabulary, column definitions, filters, and prioritized
repair steps, see [Gene-list diagnostics and repair](current/gene-list-diagnostics-and-repair.md).

### ldscore: chromosome coverage preflight

**Raised by:** direct input/coverage validation or immutable index validation. **Symptom:** the batch fails before publishing scientific outputs.

Read `diagnostics/input_issues.tsv` and `diagnostics/chromosome_scope.json`. For gene lists, inspect `coverage_status`, `selected_genes`, `covered_genes`, `missing_chromosomes`, and `uncovered_gene_ids` in `gene_list_resolution_summary.tsv`; audit rows with `coverage_status == 'uncovered'` identify affected input lines. Fix every reported issue before rerunning with `--overwrite`.

- Missing, unreadable, malformed, or mismatched PLINK trio / R²-sidecar / baseline inputs: restore valid matched artifacts. `@` requires every autosome 1–22. All safely discoverable independent issues are reported at the current gate.
- If a chromosome subset was intended, use quoted `*` patterns or exact paths selecting matching baseline and PLINK chromosome sets, for example `"baseline.*.annot.gz"` and `"panel.*"`. Use `*` to select available files without declaring a complete 1-22 suite. Check the matched files and their contents; a broad wildcard can select unintended chromosomes. For a full-suite run, restore missing files instead of narrowing the declaration. `--r2-dir` takes a literal directory containing the matching chromosome set, not a glob.
- Baseline/reference chromosome-set disagreement: supply exactly matching validated sets. Glob matches are authoritative; filenames do not establish content. A missing file may be undetectable when the remaining groups consistently cover the same subset.
- Incomplete focal/control coverage: supply matching inputs covering every selected gene, or deliberately revise the lists. Neither `resolved-only` nor chromosome-subset inputs authorize dropping cross-chromosome genes or pathways.
- Invalid immutable index: restore or rebuild the complete public index. Resolution policy never repairs an index or changes its scientific configuration.

Coverage follows unique identifier resolution and explicit gene exclusions, before SNP filtering. Unevaluated support remains blank. Fully covered genes with no retained computational reference SNPs have valid zero-support measurements; focal queries can be skipped for zero support/variance, whereas unusable controls and all-focal-skipped batches fail. See [the complete diagnostic vocabulary](current/gene-list-diagnostics-and-repair.md).

Replacing `@` with `*` changes file selection only; it does not remove genes outside the selected chromosome set. The raised error, `diagnostics/ldscore.log`, and input-issue repair text describe these alternatives. Source: [`_ldscore_preflight.validate_direct_scope` / `inspect_direct_inputs`](../src/ldsc/_ldscore_preflight.py).

For BED queries, `--padding-bp` extends the supplied intervals at both ends. Set it to 0 when your BED files already contain the intended padded regions, to avoid double padding.

### ldscore: annotation values are malformed

**Raised by:** `_kernel.annotation._validate_annotation_values()` through the
annotation and LD-score loaders · **Exception:** `LDSCInputError`
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
| 1 | The R2 parquet was written by an older LDSC version or another tool | Inspect columns; current files contain `IDX_1`, `IDX_2`, `R2`, and `SIGN_R`. The former `SIGN` column is no longer accepted; rebuild according to the [current format requirements](current/parquet-r2-format-and-read-pipeline.md#5-caveats-and-constraints). |
| 2 | The matching `chrN_meta.tsv.gz` sidecar is missing or was copied from another panel | Confirm each `chrN_r2.parquet` has a same-directory `chrN_meta.tsv.gz` with matching timestamps/provenance |
| 3 | Parquet schema metadata was stripped or edited | Inspect parquet metadata for `ldsc:sorted_by_build`, `ldsc:n_snps`, and `ldsc:sidecar_identity_sha256` |
| 4 | The R2 directory mixes chromosomes or builds from different reference-panel runs | List files in the R2 directory and compare recorded `genome_build` metadata |
| 5 | Duplicate or conflicting pair rows survived in the R2 artifact | Regenerate from a deduplicated reference-panel source and current `ldsc build-r2-panel` |

**Remedies:**

1. Regenerate the reference panel with the current `ldsc build-r2-panel`.
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

### ldscore: unusable CM for `--ld-wind-cm`

`--ld-wind-cm` needs a genetic-map coordinate (`CM`) that orders SNPs along a
chromosome. The run aborts when the reference panel's `CM` is all zero, constant,
or missing — it cannot define a genetic-distance window, and `--yes-really` does
**not** override this (it only authorizes whole-chromosome windows from *valid*
`CM`).

| Likely cause | Check |
|---|---|
| PLINK `.bim` has an all-zero `CM` column (the common PLINK default) | Inspect the third `.bim` column; it is `0` for every SNP |
| Parquet panel built without a genetic map (sidecar `CM=NA`) | Inspect `chr*_meta.tsv.gz` `CM` column |
| Genome build for the genetic map could not be determined (rsID modes) | Confirm whether `--genome-build` was passed |

**Remedies:**

1. Add a genetic map for the panel's build:
   `--genetic-map-hg38-sources <file>` or `--genetic-map-hg19-sources <file>`.
   ldscore interpolates `CM` at the `.bim` positions (PLINK backend).
2. Provide a `.bim` whose third column carries real genetic-map positions, or a
   parquet panel built with `ldsc build-r2-panel --genetic-map-<build>-sources`.
3. Use a physical-distance window instead: `--ld-wind-kb 1000` or `--ld-wind-snps`.
4. In rsID identifier modes, pass `--genome-build hg19` / `--genome-build hg38`
   so the matching genetic map can be selected.

## build-r2-panel

### build-r2-panel: no reference-panel artifacts were produced

**Raised by:** `ref_panel_builder.ReferencePanelBuilder.run()` and reference-panel
metadata cleanup paths · **Exception:** `LDSCInputError`
**Symptom:** `build-r2-panel produced no chromosome artifacts...` / `...retained no parquet metadata rows...`

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

### build-r2-panel: liftover or genetic-map configuration is incomplete

**Raised by:** `ref_panel_builder.ReferencePanelBuilder._prepare_build_state()`,
chromosome liftover setup, and genetic-map interpolation helpers · **Exception:** `LDSCUsageError` / `LDSCInputError`
**Symptom:** `build-r2-panel cannot emit <build>...` / `...genetic map...`

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

### build-r2-panel: SNP restriction does not match the source panel

**Raised by:** `ref_panel_builder._read_ref_panel_snp_restriction()` and
restriction build/column readers · **Exception:** `LDSCInputError`
**Symptom:** `build-r2-panel SNP restriction does not match the source panel build...`

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

### build-r2-panel: reference-panel artifact is incompatible

**Raised by:** `_kernel.ref_panel.ParquetR2RefPanel` and metadata sidecar readers
· **Exception:** `LDSCInputError`
**Symptom:** `Reference-panel metadata sidecar is missing...` / `...required identity/provenance keys are missing...`

**Likely causes & how to check** (most probable first):

| # | Likely cause | How to check |
|---|--------------|--------------|
| 1 | The R2 parquet was copied without its matching `chrN_meta.tsv.gz` sidecar | Confirm each `chrN_r2.parquet` has a same-directory `chrN_meta.tsv.gz` |
| 2 | The artifact was written by an older LDSC version | Inspect parquet or sidecar metadata for `artifact_type`, `snp_identifier`, and `genome_build` |
| 3 | A parent directory with multiple build children was passed without a concrete build | Check for `hg19/` and `hg38/` children under `--r2-dir` |
| 4 | The metadata sidecar lacks required SNP identity columns for the active mode | Inspect sidecar columns for SNP or CHR/POS plus A1/A2 in allele-aware modes |
| 5 | The selected chromosome/build was not generated | List `chr*_r2.parquet` files in the selected build directory |

**Remedies:**

1. Keep R2 parquet files and metadata sidecars together as one artifact family.
2. Regenerate the reference panel with the current `ldsc build-r2-panel`.
3. Pass a concrete build-specific R2 directory or set the matching genome build.

## convert-ldsc2-ldscores

The [conversion user guide](../tutorials/convert-legacy-ldscores.md) gives full layouts, table requirements, flags, and copy-and-rename commands. After logging starts, conversion exceptions and their remedies are recorded in `diagnostics/convert-ldsc2-ldscores.log` even at `--log-level ERROR`, with available issues in `diagnostics/conversion_issues.tsv.gz`. Invalid directory arguments or output collisions may fail before a log opens.

### standard filenames and directory organization

Use these patterns directly inside the selected input directories. `<chrom>` means unpadded `1` through `22`; `<prefix>` is constant within each family and can be empty. Include separators in the prefix, such as `baseline.`. `(.gz)` denotes a plain or gzip alternative; counts must be plain text.

| Input | Standard pattern | Example | Repair practice |
| --- | --- | --- | --- |
| Reference or weight scores | `<prefix><chrom>.l2.ldscore(.gz)` | `weights.1.l2.ldscore.gz` | Restore missing chromosomes or rename valid copies consistently; select the directory containing the files directly |
| Baseline annotations | `<reference-prefix><chrom>.annot(.gz)` | `baseline.1.annot.gz` | Use the same prefix as the reference scores and restore all 22 matching shards |
| Required common counts | `<reference-prefix><chrom>.l2.M_5_50` | `baseline.1.l2.M_5_50` | Restore the matching source count or rename a misnamed copy; `.l2.M` cannot substitute |
| Optional all-SNP counts | `<reference-prefix><chrom>.l2.M` | `baseline.1.l2.M` | Supply matching original counts if available; see the profile-specific missing-count policy |
| Baseline frequencies | `<prefix><chrom>.frq(.gz)` | `1000G.EUR.QC.1.frq.gz` | Rename valid `.freq` copies to `.frq` if needed, restore missing chromosomes, and select one complete family |

Place different releases in separate directories when discovery finds multiple complete families. Do not point at a parent directory, filename prefix, individual file, or glob. Counts and annotations must match the selected reference prefix; weight and frequency prefixes can differ. Renaming repairs discovery only: retain the source release and valid table contents, and decompress gzip contents before changing a compressed count file to its plain-text standard name.

### unsupported weight filenames

**Symptom:** family discovery reports unsupported `.w.l2.ldscore(.gz)` filenames, such as `1.w.l2.ldscore.gz`.

Both reference and weight suites require `<prefix><chrom>.l2.ldscore(.gz)` names across chromosomes 1–22. Rename copies in a separate input directory, for example `1.w.l2.ldscore.gz -> weights.1.l2.ldscore.gz`, and pass that directory to `--legacy-weight-dir`. Use the same correction without `.gz` for plain files. Each affected file in the inspected directory is recorded in `diagnostics/conversion_issues.tsv.gz` as `discarded_unsupported_weight_filename`. If a complete recognized family exists alongside these files, conversion uses that family and records the unused alternatives with a warning and in `legacy_ldsc2_import.ignored_files`. After a failed conversion, rerun with a new output directory or `--overwrite` to replace its diagnostics.

### conversion rejected the legacy suite

**Raised by:** `legacy_ldscore_converter.LegacyLDScoreConverter`
· **Exception:** `LDSCInputError`

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
`docs/current/legacy-ldscore-conversion.md` for the complete contract.

## regression

### regression: no SNPs remain after merging regression inputs

**Raised by:** `regression_runner.RegressionRunner.build_dataset()` and
`regression_runner.RegressionRunner.build_rg_dataset()` · **Exception:** `LDSCInputError`
**Symptom:** `h2 regression retained no overlapping...` / `rg regression retained no overlapping...`

**Likely causes & how to check** (most probable first):

| # | Likely cause | How to check |
|---|--------------|--------------|
| 1 | Sumstats and LD-score artifacts use different SNP identifier modes | Compare `snp_identifier` in each artifact's metadata |
| 2 | Coordinate-mode inputs are on different genome builds | Compare `genome_build` metadata for sumstats and the LD-score directory |
| 3 | One input was generated allele-aware and another was generated allele-blind | Check whether each artifact has A1/A2 columns and allele-aware provenance |
| 4 | The regression SNP universe or filters removed all shared SNPs | Inspect the retained SNP counts and any dropped-SNP sidecars from upstream commands |
| 5 | For `rg`, the two traits use incompatible allele conventions | Confirm both munged sumstats were generated with harmonized A1/A2 columns |

**Remedies:**

1. Regenerate sumstats and LD scores with the same `--snp-identifier` and
   `--genome-build`.
2. Use `--allow-identity-downgrade` only for same-family allele-aware/base mixes
   when duplicate base identities have been cleaned.
3. Broaden or rebuild the regression SNP universe so the traits and LD-score
   directory share retained SNPs.

### partitioned-h2: missing overlap matrix

**Raised by:** `regression_runner.estimate_partitioned_h2_batch()` /
`summarize_partitioned_h2()` · **Exception:** `LDSCInputError`
**Symptom:** `partitioned-h2 needs the annotation overlap matrix, but the LD-score directory has no ldscore.overlap.parquet...`

Overlap-aware partitioned heritability (both the functional and cell-type
regimes) needs `ldscore.overlap.parquet`, written by `ldsc ldscore` alongside
the baseline/query parquet files whenever the run has **two or more** annotation
columns. A baseline-only directory (multiple baseline annotations, no query) is
**not** an error — it runs the functional-category regime. An unpartitioned,
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

1. Partitioned-h2 requires >=2 annotation columns; regenerate the LD-score
   directory with explicit baseline (and optional query) annotations, not the
   synthetic `base`.
2. If the directory predates the sidecar, regenerate it with the current
   `ldsc ldscore`.
3. Keep `metadata.json`, the baseline/query parquet, and `ldscore.overlap.parquet`
   from the same run together.
4. `h2` and `rg` do not require the overlap matrix (the shared `h2` collinearity
   guard uses it only when present), so they run on unpartitioned or older
   directories.

## quantile-h2

### quantile-h2: fitted model or resupplied sources do not verify

**Raised by:** `quantile_h2.load_fitted_partitioned_model()` / `quantile_h2._prepare_quantile_inputs()` · **Exception:** `LDSCInputError`

**Likely causes & how to check** (most probable first):

| # | Likely cause | How to check |
|---|---|---|
| 1 | An aggregate cell-type root was supplied | Select one `diagnostics/query_annotations/<query>/` directory; aggregate rows are separate fits. |
| 2 | The fitted result predates coefficient-delete persistence | Check for `coefficient_delete_values.parquet`; rerun `partitioned-h2` with current LDSC3. |
| 3 | Not every original fitted annotation was resupplied | Compare `retained_ld_columns` in model metadata with source headers. |
| 4 | Reference metadata or annotations come from a different panel/release | Inspect SNP-universe size, annotation-sum, or overlap mismatch text. |
| 5 | Effective SNP identities are duplicated or target coverage is incomplete | Inspect `diagnostics/snp_alignment_issues.tsv.gz`. |

Resupply the exact original annotation sources and matching reference metadata. Parquet-R2 users can use the existing `chr*_meta.tsv.gz` sidecars; PLINK users should regenerate LD scores with `--export-ref-metadata`. Do not concatenate coefficient tables from different query directories.

### quantile-h2: target values are invalid or quantiles are empty

**Raised by:** `quantile_h2._read_target_annotation()` / `quantile_h2.assign_legacy_quantiles()` · **Exception:** `LDSCInputError`

By default, every target value must be numeric and finite. If exactly one token denotes missingness, pass it with `--target-missing-value`; zero is retained unless explicitly selected. If boundary ties produce an empty quantile, reduce `--num-quantiles` or use a less discrete target. Ties are intentionally kept in the lower-valued quantile for LDSC2 compatibility.

## plot

### plot: result directory is unsupported or incomplete

**Raised by:** `plotting.plot_result()` · **Exception:** `LDSCInputError`

| # | Likely cause | How to check |
|---|---|---|
| 1 | An internal per-query partitioned-h2 directory was supplied | Pass the aggregate partitioned-h2 root containing `diagnostics/metadata.json`. |
| 2 | The result predates a required current artifact | For h2, check for `diagnostics/ld_score_regression_bins.tsv`; rerun `ldsc h2` with the current package. |
| 3 | Metadata and result tables came from different runs | Inspect `diagnostics/metadata.json` and verify every relevant `files` entry exists relative to the same root. |
| 4 | The artifact type or scientific regime has no approved plot | Compare the input with the supported table in [the plotting manual](../tutorials/plotting-results.md#which-result-directory-produces-which-plot). |

The plot command does not infer from loose TSV files or reconstruct missing
diagnostics. Matplotlib is a required package dependency. If an existing
environment reports that it is unavailable, repair that environment with
`python -m pip install "matplotlib>=3.9,<4"` and rerun the command.

### plot: rg heritability source is invalid

**Raised by:** `plotting.plot_result()` and `plotting._builders._trait_h2_labels()` · **Exception:** `LDSCInputError`

| Likely cause | How to check and fix |
|---|---|
| The declared heritability table is unreadable or malformed | Check `files.h2_per_trait` in `diagnostics/metadata.json`. Restore the original readable TSV from the same rg run, with consistent field counts and valid text encoding. |
| Required columns or trait identifiers are missing | The table must contain `trait_name`, `total_h2_obs`, and `total_h2_obs_se`, and each row must have a nonempty trait name. Restore the canonical table. |
| A trait has multiple rows | Check repeated `trait_name` values and restore one single-trait fit per trait; do not substitute pair-specific fits. |
| The metadata path is absolute, escapes the result root, or points to a directory | Use a relative path to a regular file inside the same result directory. Symlinks must also resolve inside that root. |

Missing declarations/files or absent trait rows are allowed and display `failed`. Missing, nonfinite, or nonnumeric h2/SE values and negative SEs also display `failed`; these do not abort plotting. Restore valid saved inputs and rerun `ldsc plot --result-dir RESULT_DIR --overwrite` to refresh existing figures. See [heritability annotation behavior](current/plotting-module.md#heritability-annotations-in-rg-plots).
