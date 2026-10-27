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

The Python workflow mirrors the two CLI steps in [partitioned heritability](workflows/partitioned-heritability.md): `run_ldscore` writes the LD-score directory and returns an `LDScoreSource`, and `RegressionRunner.estimate_partitioned_h2_batch` fits one baseline-plus-query model per query annotation.
Set the SNP identifier mode and genome build once with `set_global_config`; `run_ldscore` rejects them as keyword arguments.

```python
from ldsc import (
    GlobalConfig, RegressionConfig, RegressionRunner,
    load_sumstats, run_ldscore, set_global_config,
)

def main():
    global_config = GlobalConfig(snp_identifier="chr_pos_allele_aware", genome_build="hg19")
    set_global_config(global_config)
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

    sumstats = load_sumstats("results/trait/trait.parquet", trait_name="trait")
    runner = RegressionRunner(global_config=global_config, regression_config=RegressionConfig())
    result = runner.estimate_partitioned_h2_batch(
        sumstats, source,
        workers=1,
        query_batch_size=100,
        output_dir="results/partitioned_h2",
        metadata={"ldscore_dir": ldscore_dir, "trait_name": "trait"},
    )
    print(result.summary)
    print(result.per_query_artifacts)  # Persistent per-query table, delete-value and metadata paths.

if __name__ == "__main__":
    main()
```

`run_ldscore` accepts the `ldscore` option names without the leading `--`, for example `ref_panel_snps_file`, `regr_snps_file`, `plink_prefix` or `query_annot_sources`.
The returned `LDScoreSource` holds the baseline LD scores and the saved query file paths, not the query LD scores; regression reads only the query columns it needs.
To reopen a saved directory in a later session, use `load_ldscore_from_dir(ldscore_dir)`.
`load_sumstats` reads munged summary statistics and recovers their SNP identity and genome build from the Parquet footer; it also accepts LDSC2 `.sumstats(.gz)` files.

`RegressionConfig` holds the regression settings, for example `RegressionConfig(intercept_h2=1.0)` for a fixed intercept or `use_common_counts=False` for all-SNP counts (CLI `--count-kind all`).
`estimate_partitioned_h2_batch` also takes `query_columns` to fit a subset of queries, or an empty list to fit the baseline annotations as one joint model, as well as `summary_sort_by` and `overwrite`.
It returns a `PartitionedH2BatchResult` with `summary`, `output_paths`, `per_query_artifacts` and `query_status`.
The in-memory `summary` has the scientific columns only; the saved `partitioned_h2.tsv` adds `query_result_dir`.

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
result = runner.estimate_partitioned_h2_batch(
    sumstats, source,
    output_dir="results/partitioned_h2",
    continue_on_query_error=True,
)
failed = result.query_status.loc[result.query_status["status"] != "success"]
print(failed[["query_annotation", "stage", "error_type", "error_message"]])
```

Direct Python calls send errors to the `LDSC` logger and do not install the `diagnostics/partitioned-h2.log` file handler; they still write `diagnostics/query_status.tsv`.

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
`estimate_partitioned_h2_batch(..., workers=4, query_batch_size=100)` fits four models at once from an importable script.
From a notebook that cannot start importable workers, run the CLI in a subprocess; see [Cell-type enrichment](#cell-type-enrichment).

## Cell-type enrichment

[Cell-type enrichment](workflows/cell-type-enrichment.md) uses the same Python calls as [partitioned heritability](#partitioned-heritability): pass the cell-type BED files, prebuilt annotation files or gene lists to `run_ldscore`, then call `RegressionRunner.estimate_partitioned_h2_batch`, which fits each cell-type query in its own baseline-plus-query model.

### Run parallel fits from a notebook

Workers are started with `spawn`, so `estimate_partitioned_h2_batch(..., workers=4)` must run from an importable script under `if __name__ == "__main__":` (see [parallel execution](../reference/parallel-execution.md)).
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

## Migrating from LDSC2

`convert_ldsc2_ldscores` is the Python form of `ldsc convert-ldsc2-ldscores` described in [migrating from LDSC2](workflows/migrating-from-ldsc2.md).
Its keyword arguments match the flags without the leading `--`, with the same defaults (`snp_identifier="rsid"`, `genome_build="auto"`, `overwrite=False`, `log_level="INFO"`); like the command, it has no common-MAF threshold.
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
LDSC2 `.sumstats(.gz)` files need no conversion: `load_sumstats` recognizes them, and `RegressionRunner` matches them to the LD-score directory automatically.

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
