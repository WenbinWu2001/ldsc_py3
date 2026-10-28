# Using LDSC3 from Python

Last updated on: 2026-10-27

This page collects the Python API for each workflow in the user guide.
Each section below assumes the setup described here and links back to the workflow page that explains the analysis.

## Setup

**Import from `ldsc`.**
Import every name from the top-level package, for example `from ldsc import GlobalConfig, run_ldscore`.
Modules and names that start with an underscore, such as `ldsc._kernel`, are private and can change without notice.

**The CLI and the Python API run the same workflows.**
Each `ldsc` command parses its flags and calls the same workflow code as its Python counterpart, so both produce the same scientific outputs.
Python workflow functions take the command's option names without the leading `--` (for example `--ld-wind-cm` becomes `ld_wind_cm`) and return result objects with the output paths, instead of printing a report.

**Set shared settings once with `GlobalConfig`.**
The SNP identifier mode, the genome build and the log level are shared by all workflows and are set once per Python session:

```python
from ldsc import GlobalConfig, get_global_config, set_global_config

set_global_config(GlobalConfig(snp_identifier="chr_pos_allele_aware", genome_build="hg19"))
print(get_global_config())
```

`GlobalConfig()` defaults to `snp_identifier="chr_pos_allele_aware"`, `genome_build="auto"` and `log_level="INFO"`, and `reset_global_config()` restores these defaults.
Workflow functions such as `run_ldscore` read the registered configuration and reject `snp_identifier`, `genome_build` and `log_level` as keyword arguments; `run_annotate`, `RegressionRunner` and `SumstatsMunger.run` also accept an explicit `global_config=`.
Saved LD-score directories and munged summary statistics keep the SNP identity and build recorded when they were written, and regression uses those records rather than the session setting.
The regression functions `run_h2`, `run_partitioned_h2`, `run_rg` and `run_quantile_h2` therefore take only the log level from the registered configuration.
These four functions, `run_ldscore` and `run_build_ref_panel` also reject, before writing anything, a keyword name the command does not have, a missing required option, a value outside the command's choices, and a non-boolean value for an on/off option.
`plot_result`, `convert_h2_scale`, `convert_ldsc2_ldscores` and `build_gene_ldscore_index` take an explicit `log_level=` for one call; when it is omitted they use the registered log level.
See [SNP identity and genome build](concepts/snp-identity-and-genome-build.md).

**Run workflows under `if __name__ == "__main__":`.**
Most numerical workflows, including LD-score computation, R² panel and index construction, regression and LD-score conversion, run in managed workers started with `spawn`, even with one worker.
A Python script must therefore run its workflow calls inside a main guard:

```python
from ldsc import GlobalConfig, run_ldscore, set_global_config

def main():
    set_global_config(GlobalConfig(snp_identifier="rsid"))
    run_ldscore(
        plink_prefix="/path/to/plink/1000G.EUR.QC.",
        ld_wind_cm=1.0,
        output_dir="results/ldscores",
    )

if __name__ == "__main__":
    main()
```

From a notebook that cannot start such workers, run the CLI in a subprocess, as shown under [Cell-type enrichment](#cell-type-enrichment).
Worker counts and process start-up are described in [parallel execution](../reference/parallel-execution.md).

## Partitioned heritability

The Python workflow mirrors the two CLI steps in [partitioned heritability](workflows/partitioned-heritability.md): `run_ldscore` writes the LD-score directory, and `run_partitioned_h2` fits one baseline-plus-query model per query annotation and writes the same result directory as `ldsc partitioned-h2`.
Set the SNP identifier mode and genome build once with `set_global_config`; `run_ldscore` rejects them as keyword arguments.

```python
from ldsc import GlobalConfig, run_ldscore, run_partitioned_h2, set_global_config

def main():
    set_global_config(GlobalConfig(snp_identifier="chr_pos_allele_aware", genome_build="hg19"))
    ldscore_dir = "results/partitioned_ldscores"

    source = run_ldscore(
        baseline_annot_sources="/path/to/baseline/baseline.@.annot.gz",
        query_annot_bed_sources="/path/to/beds/*.bed",
        r2_dir="/path/to/r2-panel/hg19",
        ld_wind_cm=1.0,
        query_batch_size=100,
        workers=1,
        output_dir=ldscore_dir,
    )

    result = run_partitioned_h2(
        sumstats_file="results/trait/trait.parquet",
        ldscore_dir=ldscore_dir,
        workers=1,
        query_batch_size=100,
        output_dir="results/partitioned_h2",
    )
    print(result.summary)
    print(result.query_status)

if __name__ == "__main__":
    main()
```

`run_ldscore` accepts the `ldscore` option names without the leading `--`, for example `ref_panel_snps_file`, `regr_snps_file`, `plink_prefix` or `query_annot_sources`.
The returned `LDScoreSource` holds the baseline LD scores and the saved query file paths, not the query LD scores; regression reads only the query columns it needs.
To reopen a saved directory in a later session, use `load_ldscore_from_dir(ldscore_dir)`.

`run_partitioned_h2` accepts the `partitioned-h2` option names without the leading `--`, for example `trait_name`, `count_kind="all"` for all-SNP counts, `intercept_h2=1.0` for a fixed intercept, `summary_sort_by`, `continue_on_query_error` and `overwrite`; see [regression options](../reference/regression-options.md).
It writes `partitioned_h2.tsv`, the per-query results and `diagnostics/partitioned-h2.log`, which `plot_result()` and `run_quantile_h2` read later.
It returns a `PartitionedH2BatchResult` with `summary`, `output_paths`, `per_query_artifacts` and `query_status`.
The in-memory `summary` has the scientific columns only; the saved `partitioned_h2.tsv` adds `query_result_dir`.

### Fit loaded inputs

`RegressionRunner.estimate_partitioned_h2_batch` fits the same models from a loaded `SumstatsTable` and `LDScoreSource`:

```python
from ldsc import RegressionConfig, RegressionRunner, load_sumstats

# Inside main(), after run_ldscore(...) returns source.
sumstats = load_sumstats("results/trait/trait.parquet", trait_name="trait")
runner = RegressionRunner(regression_config=RegressionConfig())
result = runner.estimate_partitioned_h2_batch(
    sumstats, source,
    workers=1,
    query_batch_size=100,
    output_dir="results/partitioned_h2_loaded",
    metadata={"ldscore_dir": "results/partitioned_ldscores", "trait_name": "trait"},
)
print(result.per_query_artifacts)  # Persistent per-query table, delete-value and metadata paths.
```

`load_sumstats` reads munged summary statistics and recovers their SNP identity and genome build from the Parquet footer; it also accepts legacy `.sumstats(.gz)` files.
`RegressionConfig` holds the regression settings, for example `RegressionConfig(intercept_h2=1.0)` or `use_common_counts=False` (CLI `--count-kind all`).
`estimate_partitioned_h2_batch` also takes `query_columns` to fit a subset of queries, or an empty list to fit the baseline annotations as one joint model, as well as `summary_sort_by` and `overwrite`.
It writes no workflow log and records only the `metadata` it is given.

### Inspect query values

To look at query LD scores, read them explicitly and release them when done:

```python
# Inside main(), after run_ldscore(...) returns.
print(source.output_paths)
print(source.query_columns)
selected = source.read_queries(source.query_columns[:2])
if selected is not None:
    print(selected.head())
del selected
```

`read_queries` can span several query files and returns columns in the requested order; it has no cache or width limit, so the caller owns its memory.

### Keep successful fits when queries fail

Batch regression is strict by default, matching the CLI.
`continue_on_query_error=True` is the equivalent of `--continue-on-query-error`; inspect `result.query_status` for every attempted query:

```python
# Inside main(), replacing the strict call above.
result = run_partitioned_h2(
    sumstats_file="results/trait/trait.parquet",
    ldscore_dir=ldscore_dir,
    output_dir="results/partitioned_h2",
    continue_on_query_error=True,
)
failed = result.query_status.loc[result.query_status["status"] != "success"]
print(failed[["query_annotation", "stage", "error_type", "error_message"]])
```

`run_partitioned_h2` writes the query tracebacks to `diagnostics/partitioned-h2.log`, as the command does.
`RegressionRunner.estimate_partitioned_h2_batch(..., continue_on_query_error=True)` installs no log file and sends them to the `LDSC` logger; both write `diagnostics/query_status.tsv`.

### Reuse query annotations and check SNP retention

`run_annotate` writes reusable annotation files and returns a bundle that must be used inside a `with` block.
On an annotation grid it exports aligned baseline and query files; pass both lists to `run_ldscore`:

```python
import pandas as pd
from ldsc import run_annotate, run_ldscore

def main():
    # Assumes set_global_config(...) was called as above.
    with run_annotate(
        query_annot_bed_sources="/path/to/beds/*.bed",
        ref_snp_annot_sources="/path/to/baseline/baseline.@.annot.gz",
        output_dir="results/query_annotations",
    ) as bundle:
        print(bundle.summary())
        annotation_rows = {chrom: shard.n_rows for chrom, shard in bundle.shards.items()}
        baseline_files = bundle.output_paths["baseline_annotations"]
        query_files = bundle.output_paths["query_annotations"]

    source = run_ldscore(
        baseline_annot_sources=baseline_files,
        query_annot_sources=query_files,
        r2_dir="/path/to/r2-panel/hg19",
        ref_panel_snps_file="/path/to/filters/reference_snps.tsv.gz",
        regr_snps_file="/path/to/filters/regression_snps.tsv.gz",
        ld_wind_cm=1.0,
        output_dir="results/partitioned_ldscores",
    )

    # Compare annotation rows with the regression SNPs written per chromosome.
    before = pd.Series(annotation_rows).sort_index()
    after = source.baseline_table["CHR"].astype(str).value_counts().sort_index()
    comparison = (
        before.rename("annotation_rows").to_frame()
        .join(after.rename("regression_snps"), how="outer")
        .fillna(0)
        .astype(int)
    )
    comparison["not_written"] = comparison["annotation_rows"] - comparison["regression_snps"]
    print(comparison)

if __name__ == "__main__":
    main()
```

For a BIM grid, replace `ref_snp_annot_sources` with `ref_snp_bim_sources`; that route writes query files only and supplies no baseline values.
`ref_panel_snps_file` restricts the reference SNPs, and `regr_snps_file` independently chooses the regression SNPs; omit either to use all reference SNPs or the packaged HM3 map.
Both files are identity-only filters: duplicate keys collapse, and columns such as `CM` and `MAF` do not override the reference panel.

The comparison above counts annotation rows against the regression SNPs written to `ldscore.baseline.parquet`, which every query file shares.
Rows are not written when they are absent from the reference panel after its restrictions and MAF filter, absent from the regression SNP candidates, or inside the excluded regions.
It does not count how many reference SNPs contributed to the LD scores; those stay in the broader reference set used for LD scores, counts and overlaps.
See [reference and regression SNPs](../reference/ldscore-snp-universes.md).

### Parallel fits from a notebook

Workers are started with `spawn`, so Python scripts must run workflow calls under `if __name__ == "__main__":` at every worker count (see [parallel execution](../reference/parallel-execution.md)).
`run_partitioned_h2(..., workers=4, query_batch_size=100)` fits four models at once from an importable script.
From a notebook that cannot start importable workers, run the CLI in a subprocess; see [Cell-type enrichment](#cell-type-enrichment).

## Cell-type enrichment

[Cell-type enrichment](workflows/cell-type-enrichment.md) uses the same Python calls as [partitioned heritability](#partitioned-heritability): pass the cell-type BED files, prebuilt annotation files or gene lists to `run_ldscore`, then call `run_partitioned_h2`, which fits each cell-type query in its own baseline-plus-query model.

### Run parallel fits from a notebook

Workers are started with `spawn`, so `run_partitioned_h2(..., workers=4)` must run from an importable script under `if __name__ == "__main__":` (see [parallel execution](../reference/parallel-execution.md)).
A notebook that cannot start such workers can run the CLI in a subprocess instead:

```python
import subprocess
import sys

output_dir = "results/cell_type_h2"
command = [
    sys.executable, "-m", "ldsc", "partitioned-h2",
    "--sumstats-file", "results/trait/trait.parquet",
    "--ldscore-dir", "results/cell_type_ldscores",
    "--workers", "4",
    "--query-batch-size", "100",
    "--output-dir", output_dir,
]
completed = subprocess.run(command, capture_output=True, text=True)
if completed.returncode != 0:
    print(completed.stdout)
    print(completed.stderr)
    raise RuntimeError(f"ldsc partitioned-h2 failed with exit code {completed.returncode}")
```

`sys.executable -m ldsc` runs the `ldsc` installed in the notebook's own environment.

### Read the saved results

The result directory written by either route can be read with pandas:

```python
import json
from pathlib import Path

import pandas as pd

result_dir = Path("results/cell_type_h2")
summary = pd.read_csv(result_dir / "partitioned_h2.tsv", sep="\t")
status = pd.read_csv(result_dir / "diagnostics" / "query_status.tsv", sep="\t", keep_default_na=False)
metadata = json.loads((result_dir / "diagnostics" / "metadata.json").read_text())
print(metadata["query_workers_requested"], metadata["query_workers_effective"])

# query_result_dir is relative to the result directory and stays with its row after re-sorting.
for row in summary.itertuples(index=False):
    full_model = pd.read_csv(result_dir / row.query_result_dir / "partitioned_h2_full.tsv", sep="\t")
    print(row.category, row.coefficient, row.coefficient_p, len(full_model))
```

`keep_default_na=False` keeps empty `error_message` cells of successful queries as empty strings.
`diagnostics/metadata.json` records the requested and effective worker counts of the run.

## Continuous annotations

This section is the Python counterpart of [continuous annotations](workflows/continuous-annotations.md).
`run_quantile_h2` is the Python form of `ldsc quantile-h2`: it reads one saved partitioned-heritability model, groups SNPs into quantiles of a continuous target annotation and writes `quantile_h2.tsv` and `standardized_coefficients.tsv` without refitting.

```python
from ldsc import plot_result, run_quantile_h2

def main():
    result = run_quantile_h2(
        partitioned_h2_result_dir="results/partitioned_h2",
        baseline_annot_sources="/path/to/baseline/baseline.@.annot.gz",
        target_annot_sources="/path/to/target/target.@.annot.gz",
        target_annotation="continuous_score",
        ref_metadata_sources="/path/to/r2-panel/hg19/chr@_meta.tsv.gz",
        num_quantiles=5,
        output_dir="results/continuous_score_quantiles",
    )
    print(result.quantile_h2[["quantile", "n_snps", "h2_obs", "enrichment"]])
    print(result.standardized_coefficients)
    plot_result("results/continuous_score_quantiles")

if __name__ == "__main__":
    main()
```

Pass the annotation sources used for the fitted model; for one query model, point `partitioned_h2_result_dir` to its folder under `diagnostics/query_annotations/` and add its query source.
For a gene-list query, also pass the original gene settings; `gene_list_resolution_policy="resolved-only"` corresponds to `--allow-unresolved-genes`.
`run_quantile_h2` returns a `QuantileH2Result` with `quantile_h2`, `standardized_coefficients`, `metadata` and `output_paths`.

## Migrating from legacy LDSC

`convert_ldsc2_ldscores` is the Python form of `ldsc convert-ldsc2-ldscores` described in [migrating from legacy LDSC](workflows/migrating-from-legacy-ldsc.md).
Its keyword arguments match the flags without the leading `--`, with the same defaults (`snp_identifier="rsid"`, `genome_build="auto"`, `overwrite=False`), except that an omitted `log_level` uses the registered `GlobalConfig.log_level`; like the command, it has no common-MAF threshold.
Conversion runs in a managed worker, so call it under `if __name__ == "__main__":` (see [parallel execution](../reference/parallel-execution.md)).

```python
from ldsc import convert_ldsc2_ldscores, load_ldscore_from_dir

def main():
    result = convert_ldsc2_ldscores(
        legacy_reference_dir="/path/to/legacy/baseline",
        legacy_weight_dir="/path/to/legacy/weights",
        legacy_frequency_dir="/path/to/legacy/frequencies",  # Omit for an unpartitioned suite.
        output_dir="/path/to/converted/baseline",
    )
    print(result.profile, result.n_rows, result.baseline_columns)
    print(result.output_paths)

    source = load_ldscore_from_dir("/path/to/converted/baseline")

if __name__ == "__main__":
    main()
```

The returned `LegacyLDScoreConversionResult` reports the detected suite type (`profile`), the number of retained rows, the baseline columns and the output paths, including the issue table and log.
`load_ldscore_from_dir` opens the converted directory as an `LDScoreSource` for `RegressionRunner`; see [Partitioned heritability](#partitioned-heritability).
Legacy `.sumstats(.gz)` files need no conversion: `load_sumstats` recognizes them, and `RegressionRunner` matches them to the LD-score directory automatically.

## Munge summary statistics

`SumstatsMunger.run` is the Python form of `ldsc munge-sumstats` described in [munge summary statistics](workflows/munge-sumstats.md); the CLI runs the same method.
One `MungeConfig` holds the raw file, the output settings and the QC thresholds, using the flag names without the leading `--` where they exist.
The SNP identifier mode comes from the global configuration (default `chr_pos_allele_aware`); set it with `set_global_config(GlobalConfig(snp_identifier="chr_pos"))`, for example, or pass `global_config=` to `run`.

```python
from ldsc import MungeConfig, SumstatsMunger

munger = SumstatsMunger()
sumstats = munger.run(
    MungeConfig(
        raw_sumstats_file="/path/to/raw/trait.tsv.gz",
        trait_name="trait",
        source_genome_build="auto",
        output_genome_build="hg38",
        output_dir="results/trait",
        # sumstats_snps_file="/path/to/keep_list.tsv.gz",  # Replaces the packaged HM3 map.
        # no_snp_restriction=True,  # Disables keep-list filtering.
        # liftover_chain_file="/path/to/hg19ToHg38.over.chain",  # Replaces quick liftover.
        # sumstats_format="daner-old",  # CLI --input-format; "auto" detects it.
        # column_hints={"N_col": "NEFF"},  # Only if NEFF is the intended sample size.
        # output_format="both",  # Also write trait.sumstats.gz.
        # overwrite=True,  # Replace earlier outputs for this trait.
    ),
)
```

`output_genome_build` must be chosen explicitly in coordinate-family modes, and rsID-family modes reject the build and liftover fields.
`output_format` takes `"parquet"` (default), `"tsv.gz"` or `"both"`.
The liftover options follow the same rules as the CLI; see the [liftover rules](../reference/munge-sumstats.md#liftover-rules).

### Inspect the run summary

The CLI prints a `Munge-sumstats summary:` block after a successful run; a Python run writes the same block only to `diagnostics/sumstats.log`.
Inspect the returned table and the accounting from the most recent run instead:

```python
summary = munger.build_run_summary()
print(sumstats.summary())  # Rows, alleles, trait name and source path.
print(sumstats.provenance["coordinate_provenance"]["liftover"])  # Liftover method and mapping counts.
print(summary.n_input_rows, summary.n_retained_rows)
print(summary.drop_counts)  # Rows removed at each stage; their sum is input minus retained rows.
print(summary.output_paths)  # Data files and the dropped-SNP table; the log is not listed.
```

`drop_counts` gives exclusive removals in stage order: `NA`, `coordinates`, `INFO`, `FRQ`, `P`, `sumstats_snps` (keep-list), `N`, `NSTUDY`, `liftover` and `identity` (final SNP identity cleanup).
Liftover counts start after keep-list filtering and earlier QC.
The summary also records `used_n_rule` (the sample-size strategy actually used) and `qc_thresholds`, where `None` means the INFO or MAF filter did not run because no such column was selected.
`build_run_summary()` raises `LDSCUsageError` if no run has completed.

### Preview without writing

`infer_raw_sumstats` is the Python counterpart of `--infer-only`: it reads the header and the first data row and writes nothing.

```python
from ldsc import MungeConfig, infer_raw_sumstats

preview = infer_raw_sumstats(
    "/path/to/raw/trait.tsv.gz",
    MungeConfig(raw_sumstats_file="/path/to/raw/trait.tsv.gz", output_genome_build="hg38"),
)
print(preview.detected_format, preview.column_hints, preview.missing_fields)
print(preview.runnable, preview.suggested_args)
```

## LD scores

`run_ldscore()` is the Python form of `ldsc ldscore`.
It takes the command-line option names without the leading `--` and with underscores, for example `baseline_annot_sources`, `query_annot_bed_sources`, `plink_prefix`, `r2_dir`, `ld_wind_cm`, `workers` and `output_dir`.
The SNP identifier mode, genome build and log level are not keyword arguments: set them once with `set_global_config(GlobalConfig(...))`, and `run_ldscore()` rejects `snp_identifier`, `genome_build` or `log_level` passed directly.
It also rejects, before writing anything, a name that is not an `ldscore` option (for example a misspelled `ld_wind_cn`), a value outside the command's choices and a non-boolean value for an on/off option such as `overwrite`.
Per-run SNP controls such as `ref_panel_snps_file` and `regr_snps_file` stay ordinary keyword arguments.
Put workflow calls in a script under `if __name__ == "__main__":`, as explained in [parallel execution](../reference/parallel-execution.md).

### Compute an LD-score directory

```python
from ldsc import GlobalConfig, run_ldscore, set_global_config

if __name__ == "__main__":
    set_global_config(GlobalConfig(snp_identifier="chr_pos_allele_aware", genome_build="hg38"))

    source = run_ldscore(
        output_dir="results/ldscores/baseline",
        baseline_annot_sources="/path/to/baseline/baseline.@.annot.gz",
        r2_dir="/path/to/r2-panel/hg38",
        ld_wind_cm=1.0,
        # overwrite=True,  # replace existing LD-score files and remove stale query files
    )

    print(source.baseline_columns)
    print(source.baseline_table.loc[:, ["CHR", "SNP", "POS", "regression_ld_scores"]].head())
    print(source.output_paths["baseline"])
    print(source.output_paths["metadata"])
    print(source.config_snapshot)
```

Omit `baseline_annot_sources` to compute the unpartitioned all-ones `base` LD score.
Add `query_annot_bed_sources="beds/*.bed"` to compute BED query annotations; their names are then in `source.query_columns`.

`run_ldscore()` always writes an LD-score directory and returns an `LDScoreSource`.
The returned object holds the shared baseline table, counts and metadata, but not the query LD scores.
`output_paths` lists the scientific files and diagnostics, but not `diagnostics/ldscore.log`.

### Read selected query columns

```python
from ldsc import load_ldscore_from_dir

source = load_ldscore_from_dir("results/ldscores/bed_queries")
queries = source.read_queries(["enhancers", "promoters"])
```

`read_queries()` loads the named query columns, even when they are stored in different query batch files, and returns them in the requested order after the SNP columns.
It does not cache values or limit how many columns you request, so the memory for the returned table is yours to budget.
`load_ldscore_from_dir()` opens any existing LD-score directory, including one written by the command line.

### Read one chromosome

Ordinary readers such as `pandas.read_parquet("ldscore.baseline.parquet")` read whole files.
For large outputs, the row-group records in `metadata.json` let you read one chromosome:

```python
import json
from pathlib import Path

import pyarrow.parquet as pq

ldscore_dir = Path("results/ldscores/baseline")
metadata = json.loads((ldscore_dir / "metadata.json").read_text())

if metadata["row_group_layout"] != "one_per_chromosome":
    raise ValueError(f"Unsupported row-group layout: {metadata['row_group_layout']!r}")

baseline_row_group = {
    entry["chrom"]: entry["row_group_index"] for entry in metadata["baseline_row_groups"]
}
parquet_file = pq.ParquetFile(ldscore_dir / "ldscore.baseline.parquet")
baseline_chr22 = parquet_file.read_row_group(baseline_row_group["22"]).to_pandas()
print(baseline_chr22.head())
```

For query files, iterate over `metadata["query_batches"]`: each entry has `file`, `query_columns` and `row_groups`, with row-group records shaped like the baseline ones.
`query_batches` is empty when the directory has no query annotations.
To read selected query columns genome-wide, `read_queries()` is simpler.

### Calculate small LD scores without writing files

For a small calculation that fits in one query batch, `LDScoreCalculator.run(..., output_config=None)` computes LD scores from prepared in-memory inputs and writes no annotation, diagnostic, log or output files; it may still read existing reference-panel files.

```python
from ldsc import AnnotationBundle, GlobalConfig, LDScoreCalculator, LDScoreConfig

if __name__ == "__main__":
    config = GlobalConfig(snp_identifier="rsid", genome_build="hg19")
    # metadata, baseline_values and query_values are aligned in-memory DataFrames.
    # ref_panel is an already-configured PlinkRefPanel or ParquetR2RefPanel.
    with AnnotationBundle.from_frames(
        metadata, baseline_values, query_values, config_snapshot=config,
    ) as annotations:
        result = LDScoreCalculator().run(
            annotations, ref_panel, LDScoreConfig(ld_wind_cm=1, query_batch_size=100),
            config, output_config=None,
        )
    # result.query_table and the in-memory diagnostics remain available after the block.
```

- `metadata` must contain `CHR`, a one-based integer `POS`, and `SNP` for an rsID-family mode; supply both `A1` and `A2` or neither.
- The frames align by position; annotation values must be finite numbers, and annotation names must be unique and differ from the metadata columns.
- `from_frames()` stores columns of exact 0/1 values as binary and copies other columns as float32; reads always return float32.
- It drops every member of a group of SNPs with the same SNP identity, does not infer the genome build and performs no liftover, so configure the reference panel for the same SNP identifier mode and genome build.
- The calculator takes already-resolved `regression_snps` and `regression_regions` arguments; omitting them uses every retained reference SNP as a regression SNP and excludes no regions, unlike `run_ldscore()`, which applies the HM3 map and the MHC and centromere exclusion by default.
  Pass the same restrictions when you compare the two.
- More than one query batch requires an output directory, `LDScoreOutputConfig(output_dir=...)`.

The result is an in-memory `LDScoreResult` rather than an `LDScoreSource`.

### Check a CHR/POS table's genome build

`infer_chr_pos_build()` reports the inferred genome build of a `CHR`/`POS` table; `resolve_chr_pos_table()` also returns the table with canonical chromosome labels and one-based positions.

```python
import pandas as pd
from ldsc import infer_chr_pos_build, resolve_chr_pos_table

# The table needs enough SNPs that overlap HapMap3 to support inference.
table = pd.read_csv("filters/snps_chr_pos.tsv.gz", sep="\t")

normalized, inference = resolve_chr_pos_table(table, context="my SNP table")
decision = infer_chr_pos_build(normalized.loc[:, ["CHR", "POS"]], context="normalized SNP table")

print(inference.genome_build)
print(inference.coordinate_basis)
print(decision.summary_message)
```

Workflows run the same check themselves when the genome build is `auto`.

### Configuration and provenance

`GlobalConfig` carries only the shared settings `snp_identifier`, `genome_build` and `log_level`.
Every result keeps a frozen `config_snapshot` of the settings it was computed with; changing the registered `GlobalConfig` later does not alter existing results.
Combining results computed under incompatible SNP identifier modes or genome builds raises `ConfigMismatchError`.

## Gene-list LD scores

The gene-list workflow uses the same `run_ldscore()` call as other LD scores; its keyword arguments are the command-line option names with underscores.
Put workflow calls in a script under `if __name__ == "__main__":`, as explained in [parallel execution](../reference/parallel-execution.md).

### Compute gene-list LD scores directly

```python
from ldsc import GlobalConfig, run_ldscore, set_global_config

if __name__ == "__main__":
    set_global_config(GlobalConfig(snp_identifier="rsid", genome_build="hg19"))

    source = run_ldscore(
        output_dir="results/ldscores/gene-lists",
        baseline_annot_sources="/path/to/baseline/baseline.@.annot.gz",
        plink_prefix="/path/to/plink/1000G.EUR.QC.",
        query_annot_gene_list_sources="/path/to/gene-lists/*.txt",
        gene_coordinate_file="/path/to/genes/gene-coordinates.hg19.tsv.gz",
        padding_bp=100000,
        gene_exclude_regions="mhc",
        ld_wind_cm=1.0,
        # gene_list_resolution_policy="resolved-only",  # exploratory screens only
    )

    print(source.query_columns)
    print(source.query_statuses)
```

Python selects the resolution policy with `gene_list_resolution_policy="strict"` (the default) or `"resolved-only"`, the equivalent of `--allow-unresolved-genes`; result metadata and diagnostic tables use the same two names.
`padding_bp` must be given for gene lists; passing it with prebuilt query annotations, without query annotations, or with an index raises an error even when it is 0.
Python calls return the query statuses and diagnostic paths and print no console notice.

### Build a gene LD-score index

`build_gene_ldscore_index()` is the Python form of `ldsc build-gene-ldscore-index`; its configuration fields carry the command-line option names, and the genetic-map sources are separate keyword arguments.

```python
from ldsc import GeneLDScoreIndexBuildConfig, build_gene_ldscore_index

if __name__ == "__main__":
    index_dir = build_gene_ldscore_index(
        GeneLDScoreIndexBuildConfig(
            baseline_annot_sources=("/path/to/baseline/baseline.@.annot.gz",),
            plink_prefix="/path/to/plink/1000G.EUR.QC.",
            output_dir="/path/to/indexes/baseline_100kb",
            gene_coordinate_file="/path/to/genes/gene-coordinates.hg19.tsv.gz",
            genome_build="hg19",
            snp_identifier="rsid",
            ld_wind_cm=1.0,
            padding_bp=100000,
        ),
        # genetic_map_hg19_sources="/path/to/maps/genetic_map_chr1.txt,...",
        # overwrite=True,  # replace an existing complete index
    )
```

### Open a gene LD-score index

```python
from ldsc import load_gene_ldscore_index

index = load_gene_ldscore_index("/path/to/indexes/baseline_100kb")
print(index.index_id)
print(index.chromosomes)
print(index.snp_identifier, index.genome_build)
```

Opening an index validates each chromosome and releases its data before loading the next.
The returned object holds the shared metadata, the embedded gene coordinate catalog (`gene_catalog`), per-gene support counts and the paths of the chromosome files, not the stored operators.

### Assemble gene-list LD scores from an index

`run_ldscore(gene_ldscore_index_dir=...)` is the Python form of `ldsc ldscore --gene-ldscore-index-dir`, described in [Calculate LD scores from gene lists](workflows/gene-list-ld-scores.md#step-2-assemble-gene-list-ld-scores).

```python
from ldsc import run_ldscore

if __name__ == "__main__":
    source = run_ldscore(
        output_dir="results/ldscores/gene-lists",
        gene_ldscore_index_dir="/path/to/indexes/baseline_100kb",
        query_annot_gene_list_sources="/path/to/gene-lists/*.txt",
        # control_gene_list_file="/path/to/gene-lists/control-genes.txt",
        # workers=8,
    )

    print(source.query_columns)
    print(source.config_snapshot)
```

The gene LD-score index supplies the SNP identifier mode and genome build, so the registered `GlobalConfig.snp_identifier` and `GlobalConfig.genome_build` are ignored for this call, and `source.config_snapshot` records the index values.
The keyword arguments that indexed mode accepts are those of the command: `control_gene_list_file`, `gene_list_resolution_policy`, `workers`, `query_batch_size`, `overwrite` and `output_dir`; the scientific options stored in the index are rejected.
The resulting directory opens with `load_ldscore_from_dir()` like any other LD-score directory.

## Annotations

`run_annotate()` is the Python form of `ldsc annotate`.
Its keyword arguments carry the command-line option names with underscores: exactly one of `ref_snp_annot_sources` or `ref_snp_bim_sources`, exactly one of `query_annot_bed_sources` or `query_annot_gene_list_sources`, and `output_dir`.
The SNP identifier mode and genome build come from `global_config=GlobalConfig(...)` or the registered global configuration, and `genome_build=` overrides the genome build for one call.
The former keyword `baseline_annot_sources` is not accepted; use `ref_snp_annot_sources` for an annotation grid.

### Project BED files onto a SNP grid

```python
from ldsc import GlobalConfig, run_annotate

with run_annotate(
    ref_snp_bim_sources="panel/reference.@.bim",
    query_annot_bed_sources="baseline_beds/*.bed",
    global_config=GlobalConfig(snp_identifier="chr_pos_allele_aware", genome_build="hg19"),
    output_dir="annotations/from_bim",
    # overwrite=True,  # replace existing annotation files and remove stale ones
) as bundle:
    print(bundle.query_columns)
    print(bundle.chromosomes)
    print(bundle.output_paths["query_annotations"])
```

For an annotation grid, replace `ref_snp_bim_sources` with `ref_snp_annot_sources`.
That result also has `bundle.output_paths["baseline_annotations"]`; pass both path lists to `run_ldscore(baseline_annot_sources=..., query_annot_sources=..., ...)`.

### Project gene lists

```python
from ldsc import GlobalConfig, run_annotate

with run_annotate(
    ref_snp_annot_sources="/path/to/baseline/baseline.@.annot.gz",
    query_annot_gene_list_sources="pathways/*.txt",
    gene_coordinate_file="/path/to/genes/gene-coordinates.hg19.tsv.gz",
    padding_bp=100000,
    gene_exclude_regions="none",
    global_config=GlobalConfig(snp_identifier="rsid", genome_build="hg19"),
    output_dir="annotations/pathways",
    # gene_list_resolution_policy="resolved-only",  # exploratory screens only
) as bundle:
    print(bundle.query_columns)
    print(bundle.query_statuses)
```

`gene_list_resolution_policy=None` (the default) means strict resolution; `"resolved-only"` is the equivalent of `--allow-unresolved-genes`.
`run_annotate()` prints no console notice; it returns the query statuses and output paths instead.

### Bundle lifetime

The returned bundle refers to the exported annotation files (only query files for a BIM grid), not to your original inputs.
Construction scratch is released before the call returns; reading a chromosome from the bundle can prepare private data again under the output directory.
Use `with run_annotate(...) as bundle:` or call `bundle.close()` when you are done; closing keeps every exported file and your original inputs.

<a id="r2-panels"></a>

## R² panels

This section is the Python counterpart of [Build and query R² panels](workflows/r2-panels.md).

### Build a panel

The Python builder keeps the names `ReferencePanelBuilder`, `ReferencePanelBuildConfig` and `run_build_ref_panel`; the command is `build-r2-panel`.
`run_build_ref_panel()` takes the command's options as keyword arguments.
It reads the SNP identifier mode and log level from the registered `GlobalConfig`, so call `set_global_config()` first; passing `snp_identifier`, `genome_build` or `log_level` to it is an error, and `GlobalConfig.genome_build` is ignored: `source_genome_build` gives the build of the PLINK input.
Like `run_ldscore()`, it rejects names that are not `build-r2-panel` options and values outside the command's choices before writing anything.

```python
from ldsc import GlobalConfig, run_build_ref_panel, set_global_config

if __name__ == "__main__":
    set_global_config(GlobalConfig(snp_identifier="chr_pos_allele_aware", log_level="INFO"))
    result = run_build_ref_panel(
        plink_prefix="tests/fixtures/minimal_external_resources/plink/hm3_chr22_subset",
        source_genome_build="hg38",
        ld_wind_snps=10,
        output_dir="tutorial_outputs/r2_panel_python",
    )
    print(result.chromosomes)
    print(result.output_paths["r2_hg38"], result.output_paths["meta_hg38"])
```

`result.output_paths` maps `r2_<build>` and `meta_<build>` to the per-chromosome files, plus `metadata`; it does not list the log.
`run_build_ref_panel()` writes the workflow log under `diagnostics/`.
`ReferencePanelBuilder(global_config=...).run(ReferencePanelBuildConfig(...))` builds the same files without a workflow log.

### Inspect a panel

Read one row group instead of loading a large pair table in full.
A direct Pandas or PyArrow read returns the stored 16-bit integers; divide `R2` by the footer's `ldsc:r2_scale`.

```python
import pandas as pd
import pyarrow.parquet as pq

root = "tutorial_outputs/r2_panel_chr22/hg38"
meta = pd.read_csv(f"{root}/chr22_meta.tsv.gz", sep="\t", comment="#")
parquet = pq.ParquetFile(f"{root}/chr22_r2.parquet")
pairs = parquet.read_row_group(0).to_pandas()
scale = float(parquet.schema_arrow.metadata[b"ldsc:r2_scale"])
pairs["r2_decoded"] = pairs["R2"].astype("float32") / scale
print(meta.head())
print(pairs.head())
```

Never sort or edit the sidecar apart from its pair table: its row order defines `IDX_1` and `IDX_2`.

### Query pairs

`query_r2()` opens a panel, queries a DataFrame of pairs and returns the input rows with `r2`, `sign_r`, `r` and `status`.
The pair columns are the same as for `ldsc query-r2`.
This example queries the first two SNPs of the sidecar read above, in the sidecar's allele orientation, and saves the pairs for the command-line query:

```python
from ldsc import query_r2

query = pd.DataFrame({
    f"{column}_{endpoint}": [meta.iloc[row][column]]
    for endpoint, row in ((1, 0), (2, 1))
    for column in ("CHR", "POS", "A1", "A2")
})
result = query_r2(query, panel_dir=root, genome_build="hg38")
print(result[["r2", "sign_r", "r", "status"]])
query.to_csv("tutorial_outputs/r2_panel_chr22/pairs.tsv", sep="\t", index=False)
```

In Python, `sign_r` is a nullable `Int8` column; the saved TSV writes missing values as `NaN`.
A missing `panel_dir` is an error.

For repeated queries, open the panel once.
`R2Panel.open()` checks every chromosome's Parquet footer immediately and loads each chromosome's sidecar and lookup state on first use, then reuses them:

```python
from ldsc import R2Panel

panel = R2Panel.open(root, genome_build="hg38")
first = panel.query_pairs(query)
print(panel.chromosomes, panel.snp_identifier, panel.n_samples)
```

`unbiased_r2_to_pearson_r(r2_adj, n, sign=None)` converts unbiased R² values to Pearson r for reference sample size `n`; without `sign` it returns $|r|$.
`query_r2()` and `R2Panel.query_pairs()` apply it with the panel's `ldsc:n_samples` and the aligned `sign_r`.
See the [query reference](../reference/query-r2.md#python-api).

## SNP heritability

This section is the Python counterpart of [Estimate SNP heritability](workflows/heritability.md).

Compute or reuse unpartitioned LD scores, then fit with `run_h2`, the Python form of `ldsc h2`:

```python
from ldsc import GlobalConfig, convert_h2_scale, plot_result, run_h2, run_ldscore, set_global_config

def main():
    set_global_config(GlobalConfig(snp_identifier="chr_pos_allele_aware", genome_build="auto"))

    run_ldscore(
        output_dir="results/unpartitioned_ldscores",
        r2_dir="/path/to/resources/r2_panel/hg38",
        ld_wind_cm=1.0,
    )

    summary = run_h2(
        sumstats_file="results/trait_sumstats/trait.parquet",
        ldscore_dir="results/unpartitioned_ldscores",
        output_dir="results/trait_h2",
    )
    print(summary[["trait_name", "total_h2_obs", "total_h2_obs_se", "intercept"]])

    plot_result("results/trait_h2")
    convert_h2_scale("results/trait_h2", samp_prev=0.5, pop_prev=0.01)

if __name__ == "__main__":
    main()
```

`run_h2` accepts the `h2` option names without the leading `--`, for example `trait_name`, `count_kind`, `intercept_h2`, `two_step_cutoff`, `chisq_max`, `n_blocks`, `samp_prev`, `pop_prev` and `overwrite`; see [regression options](../reference/regression-options.md).
It writes `h2.tsv`, the regression-bin diagnostic, `diagnostics/metadata.json` and `diagnostics/h2.log`, and returns the one-row summary saved in `h2.tsv`.
`plot_result()` and `convert_h2_scale()` read that directory; the conversion above writes the liability-scale table under `results/trait_h2/postprocessing/liability-scale/`.
`run_ldscore()` without annotations computes the all-ones `base` annotation.
LD-score directories from older versions that lack the `query_batches` manifest in `metadata.json` must be regenerated.
To munge in the same script, run `SumstatsMunger().run(MungeConfig(...))` and pass the munged Parquet file it writes as `sumstats_file`; see [Munge summary statistics](#munge-summary-statistics).

### Fit loaded inputs

`RegressionRunner` fits the same model from a loaded `SumstatsTable` and LD scores, and writes nothing:

```python
from ldsc import RegressionConfig, RegressionRunner, load_ldscore_from_dir, load_sumstats

# Inside main().
sumstats = load_sumstats("results/trait_sumstats/trait.parquet", trait_name="trait")
ldscore_result = load_ldscore_from_dir("results/unpartitioned_ldscores")
runner = RegressionRunner(regression_config=RegressionConfig())
dataset = runner.build_dataset(sumstats, ldscore_result)
h2 = runner.estimate_h2(dataset)
print(float(h2.tot), float(h2.tot_se), float(h2.intercept))
```

`RegressionConfig` holds the regression options: `n_blocks`, `use_common_counts` (the CLI's `--count-kind`), `intercept_h2`, `two_step_cutoff`, `chisq_max`, `samp_prev`, `pop_prev` and `allow_identity_downgrade`.
`RegressionConfig(use_intercept=False)` is equivalent to fixing the intercept at 1.
`SumstatsMunger().run(MungeConfig(...))` returns the same `SumstatsTable` that `load_sumstats()` reads back.

**Provenance checks.**
Munged Parquet files store the `GlobalConfig` used to make them in their footer, together with the trait label; `load_sumstats()` restores both as `SumstatsTable.config_snapshot` and `trait_name`.
`build_dataset()` checks that the summary statistics and the LD scores agree on settings such as `snp_identifier` and `genome_build`, and raises an error instead of merging incompatible inputs.
Inspect `sumstats.config_snapshot` and `dataset.config_snapshot` to see the settings used.
A legacy `.sumstats.gz` file has no stored settings and is matched to the LD scores by rsID; a Parquet file without the footer is rejected.

**Weights.**
The LD-score table includes a `regression_ld_scores` column (legacy LDSC's `w_ld`); the final regression weights are computed from it during fitting.

**Saved results.**
`estimate_h2()` returns the kernel estimator and writes nothing.
`plot_result()` and `convert_h2_scale()` need a result directory written by `run_h2` or `ldsc h2`.

## Genetic correlation

This section is the Python counterpart of [Estimate genetic correlation](workflows/genetic-correlation.md).

`run_rg` is the Python form of `ldsc rg`: it fits two or more traits against one LD-score directory and writes the same result directory, which `plot_result()` reads:

```python
from ldsc import plot_result, run_rg

def main():
    result = run_rg(
        sumstats_sources=[
            "results/traits/trait_1/trait_1.parquet",
            "results/traits/trait_2/trait_2.parquet",
        ],
        ldscore_dir="results/unpartitioned_ldscores",
        output_dir="results/trait_1_trait_2_rg",
    )
    print(result.rg)
    print(result.h2_per_trait[["trait_name", "total_h2_obs", "total_h2_obs_se"]])

    plot = plot_result("results/trait_1_trait_2_rg")
    print(plot.path)

if __name__ == "__main__":
    main()
```

`run_rg` accepts the `rg` option names without the leading `--`; `sumstats_sources` takes a list of paths or quoted `*` patterns.
The result has four parts: `result.rg`, `result.rg_full` and `result.h2_per_trait` match `rg.tsv`, `rg_full.tsv` and `h2_per_trait.tsv`, and `result.per_pair_metadata` holds one record per pair.
A pair whose fit fails is recorded with `NaN` estimates and `status="failed"`, and the other pairs are still fitted.
Every trait needs `A1` and `A2`; inputs without them stop the call before any fit.
Inputs that share a trait label are relabeled `label@<suffix>` from their paths, as in [inputs that share a trait label](workflows/genetic-correlation.md#inputs-that-share-a-trait-label).
The `h2_per_trait` columns `total_h2_obs` and `total_h2_obs_se` are the values the plot shows on the heatmap diagonal or beside the anchor plot.

**Anchor, intercept and prevalence options.**
`anchor_trait="trait_1"` (a trait name or input path) fits only the pairs with that trait, and `plot_result()` then draws the anchor forest plot.
`intercept_h2=1.0` and `intercept_gencov=0.0` fix the intercepts, and `write_per_pair_detail=True` also writes `diagnostics/pairs/`.
For liability-scale estimates, `samp_prev=[0.5, float("nan")]` and `pop_prev=[0.01, float("nan")]` give one value per trait in resolved `sumstats_sources` order, with `nan` for a quantitative trait; `prevalence_manifest` matches prevalences by trait name instead.

### Fit loaded tables

`RegressionRunner.estimate_rg_pairs()` fits loaded `SumstatsTable` objects and returns the same tables; save them with `RgDirectoryWriter` to get a directory that `plot_result()` can read:

```python
from ldsc import (
    RegressionConfig, RegressionRunner, RgDirectoryWriter, RgOutputConfig,
    load_ldscore_from_dir, load_sumstats,
)

# Inside main().
trait_1 = load_sumstats("results/traits/trait_1/trait_1.parquet", trait_name="trait_1")
trait_2 = load_sumstats("results/traits/trait_2/trait_2.parquet", trait_name="trait_2")
ldscore_result = load_ldscore_from_dir("results/unpartitioned_ldscores")

runner = RegressionRunner(regression_config=RegressionConfig(intercept_h2=1.0, intercept_gencov=0.0))
result = runner.estimate_rg_pairs([trait_1, trait_2], ldscore_result, anchor_index=0)
RgDirectoryWriter().write(result, RgOutputConfig(output_dir="results/trait_1_trait_2_rg_loaded"))
```

`anchor_index` is the input position of the anchor trait.
Tables that share a trait label are relabeled from their `source_path`; tables loaded with `load_sumstats()` carry their path, and tables that share a label without a path must be given distinct `trait_name` values.
Exporting the tables with `to_csv()` does not create a plotting input; `plot_result()` needs the directory and metadata written by `RgDirectoryWriter` or `run_rg`.
`SumstatsMunger().run(MungeConfig(...))` returns a `SumstatsTable` that can be passed to `estimate_rg_pairs()` directly; see [Munge summary statistics](#munge-summary-statistics).
Tables made by the munger, and Parquet files read by `load_sumstats()`, carry the `GlobalConfig` used to make them; their settings are checked against the LD scores before fitting.
A legacy `.sumstats.gz` file is matched to the LD scores by rsID; a Parquet file without the footer is rejected.

`RegressionRunner.estimate_rg()` fits one pair and returns the low-level estimator object; prefer `estimate_rg_pairs()` or `run_rg`, which return the same tables as the command.
