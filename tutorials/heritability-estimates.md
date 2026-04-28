# Heritability Estimates

Goal: estimate SNP heritability for one trait with the refactored package, starting from raw summary statistics and LD scores built from an R2-table reference panel.

The examples below assume chromosome-pattern inputs such as `annotations/baseline.1.annot.gz`, `r2/reference.1.parquet`, and `r2/reference_metadata.1.tsv.gz`.
For parquet-backed LD scores, `r2/reference.@.parquet` is the canonical pair table
and `r2/reference_metadata.@.tsv.gz` is the required per-SNP sidecar. The sidecar
defines the reference SNP universe and supplies `MAF`/`CM`; the parquet is used
only for row-group-pruned LD window queries.

Path-token rules:

- scalar inputs such as raw sumstats must resolve to exactly one file
- group inputs may use globs or chromosome-suite tokens such as `baseline.@.annot.gz`
- output directories are literal destinations

Resolution behavior:

- there is no separate per-chromosome argument for LD-score inputs anymore; `baseline_annot_paths`, `r2_paths`, and `metadata_paths` all accept the unified token format
- a group token may resolve to multiple files
- when multiple matched files are chromosome-sharded and the chromosome can be inferred from the filenames, only the files for the active chromosome are used in that chromosome pass
- otherwise the matched files are read and filtered by `CHR` internally

## Python API

The Python workflow is the most direct end-to-end path because `run_ldscore(...)` returns one merged in-memory `LDScoreResult` that `RegressionRunner` can consume immediately.

```python
from ldsc import (
    GlobalConfig,
    MungeConfig,
    RegressionConfig,
    RegressionRunner,
    SumstatsMunger,
    load_sumstats,
    run_ldscore,
    set_global_config,
)

GLOBAL_CONFIG = GlobalConfig(
    snp_identifier="rsid",
    genome_build="hg38",
)
set_global_config(GLOBAL_CONFIG)

sumstats = SumstatsMunger().run(
    MungeConfig(
        sumstats_path="data/trait*.tsv.gz",
        trait_name="trait",
        column_hints={"snp": "SNP", "a1": "A1", "a2": "A2", "p": "P", "N_col": "N"},
        output_dir="tutorial_outputs/trait",
        signed_sumstats_spec="BETA,0",
    ),
    global_config=GLOBAL_CONFIG,
)

# If you already have a curated .sumstats.gz artifact on disk, set the intended
# GlobalConfig first, then load it directly. The loaded table will warn because
# the original munge-time config is not recoverable from disk.
# sumstats = load_sumstats("tutorial_outputs/trait/sumstats.sumstats.gz", trait_name="trait")

ldscore_result = run_ldscore(
    output_dir="tutorial_outputs/trait_ldscores",
    baseline_annot_paths="annotations/baseline.@.annot.gz",
    r2_paths="r2/reference.@.parquet",
    metadata_paths="r2/reference_metadata.@.tsv.gz",
    regression_snps_path="filters/hapmap3.txt",
    ld_wind_cm=1.0,
)

runner = RegressionRunner(global_config=GLOBAL_CONFIG, regression_config=RegressionConfig())
dataset = runner.build_dataset(sumstats, ldscore_result)
h2 = runner.estimate_h2(dataset)

print("h2 =", float(h2.tot))
print("h2_se =", float(h2.tot_se))
print("intercept =", float(h2.intercept))
print("ldscore_rows =", len(ldscore_result.baseline_table))
print("embedded_weight_column =", "regr_weight" in ldscore_result.baseline_table.columns)
print("sumstats_config =", sumstats.config_snapshot)
print("dataset_config =", dataset.config_snapshot)
```

`SumstatsMunger` captures `GlobalConfig` provenance into `SumstatsTable.config_snapshot`, and `RegressionRunner.build_dataset()` checks that the sumstats and LD-score results agree on critical settings such as `snp_identifier` and `genome_build`. If they do not, the workflow raises `ConfigMismatchError` instead of silently merging inconsistent artifacts.

When you set `regression_snps_path` during LD-score calculation, the same regression SNP subset defines the rows of the normalized `baseline_table`. Regression weights come from the embedded `regr_weight` column by default.

## CLI

The regression CLI reads the canonical LD-score result directory written by
`ldsc ldscore`.

```bash
ldsc munge-sumstats \
  --sumstats-path "data/trait*.tsv.gz" \
  --snp SNP \
  --a1 A1 \
  --a2 A2 \
  --p P \
  --N-col N \
  --signed-sumstats BETA,0 \
  --output-dir tutorial_outputs/trait

ldsc h2 \
  --sumstats-path tutorial_outputs/trait/sumstats.sumstats.gz \
  --ldscore-dir tutorial_outputs/trait_ldscores \
  --count-kind m_5_50 \
  --output-dir tutorial_outputs/trait_h2
```

The command writes `tutorial_outputs/trait_h2/h2.tsv`.
