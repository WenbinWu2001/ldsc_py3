# Inputs and Outputs

Last updated on: 2026-10-27

This page explains how LDSC3 commands find their input files, check them, and write their results.
The exact rules are in [input paths](../../reference/input-paths.md) and [logs and outputs](../../reference/logs-and-outputs.md).

## Naming input files

Options that take files accept path tokens:

- an exact path;
- a glob such as `annotations/*.annot.gz`;
- a chromosome suite such as `baseline.@.annot.gz`, where `@` stands for chromosomes 1 to 22;
- for PLINK input, an exact prefix or member file, a plain suite stem such as `1000G.EUR.QC.`, a glob, or an `@` pattern.

For PLINK input, the chromosomes come from the `.bim` contents, every selected `.bed`/`.bim`/`.fam` trio must be complete, and two trios cannot claim the same chromosome; see [PLINK resolution](../../reference/input-paths.md#plink-prefix-resolution).
All path tokens are resolved before any computation starts.

## How inputs are checked

Commands check their inputs in stages before substantial computation, and report every independent problem found at a stage in one error with a machine-readable list of issues.
Content checks such as row alignment, SNP identity cleanup and support after filtering happen later, so passing the path checks does not by itself make a dataset valid.
See [staged validation](../../reference/logs-and-outputs.md#staged-input-validation) and its [troubleshooting entry](../troubleshooting.md#staged-input-validation).

## Where results go

Commands that write results take `--output-dir`, which is always a directory:

- a missing directory is created, and an existing one is reused;
- each command writes fixed file names inside it, so you name a run by naming its directory;
- if any file the command owns already exists, the run stops before writing anything, including its log;
- `--overwrite` (Python: `overwrite=True`) replaces only the files the command owns, never unrelated files.

`plot` and `convert-h2-scale` write to fixed folders inside the result directory they read: `<result-dir>/plots/` and `<h2-result-dir>/postprocessing/liability-scale/`.
Headline results are at the top of the output directory, and logs and other diagnostics are under `diagnostics/`.
A failed run that was allowed to overwrite leaves `RUN_FAILED.txt`, which the next successful run removes; see [failure markers](../../reference/logs-and-outputs.md#failed-overwrite-markers).

## Following a run

Each command writes a workflow log under `diagnostics/`, marking validation, preparation, computation and publication with elapsed times.
The terminal shows errors and a few deliberate notices, for example when queries were skipped; everything else is in the log.
See [logs and outputs](../../reference/logs-and-outputs.md) for log names and what each level records.
