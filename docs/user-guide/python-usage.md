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
