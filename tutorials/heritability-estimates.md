# Heritability Estimates

Goal: estimate SNP heritability for one trait with the refactored package, starting from raw summary statistics and an R2-table reference panel.

The examples below assume chromosome-pattern inputs such as `annotations/baseline.1.annot.gz`, `r2/reference.1.parquet`, and `r2/reference_metadata.1.tsv.gz`.

Path-token rules:

- scalar inputs such as raw sumstats must resolve to exactly one file
- group inputs may use globs or chromosome-suite tokens such as `baseline.@.annot.gz`
- output prefixes are literal destinations

Resolution behavior:

- there is no separate per-chromosome argument for LD-score inputs anymore; `baseline_annot`, `r2_table`, and `frqfile` all accept the unified token format
- a group token may resolve to multiple files
- when multiple matched files are chromosome-sharded and the chromosome can be inferred from the filenames, only the files for the active chromosome are used in that chromosome pass
- otherwise the matched files are read and filtered by `CHR` internally

## Python API

```python
from ldsc import (
    GlobalConfig,
    MungeConfig,
    RawSumstatsSpec,
    RegressionConfig,
    RegressionRunner,
    SumstatsMunger,
    load_sumstats,
    run_ldscore,
    set_global_config,
)

GLOBAL_CONFIG = GlobalConfig(snp_identifier="rsid")
set_global_config(GLOBAL_CONFIG)

sumstats = SumstatsMunger().run(
    RawSumstatsSpec(
        path="data/trait*.tsv.gz",
        trait_name="trait",
        column_hints={"snp": "SNP", "a1": "A1", "a2": "A2", "p": "P", "N_col": "N"},
    ),
    MungeConfig(
        out_prefix="tutorial_outputs/trait",
        signed_sumstats_spec="BETA,0",
    ),
    global_config=GLOBAL_CONFIG,
)

# If you already have a curated .sumstats.gz artifact on disk, set the intended
# GlobalConfig first, then load it directly. The loaded table will warn because
# the original munge-time config is not recoverable from disk.
# sumstats = load_sumstats("tutorial_outputs/trait.sumstats.gz", trait_name="trait")

ldscore_result = run_ldscore(
    out="tutorial_outputs/trait_ldscores",
    baseline_annot="annotations/baseline.@.annot.gz",
    r2_table="r2/reference.@.parquet",
    frqfile="r2/reference_metadata.@.tsv.gz",
    ld_wind_cm=1.0,
)

runner = RegressionRunner(global_config=GLOBAL_CONFIG, regression_config=RegressionConfig())
dataset = runner.build_dataset(sumstats, ldscore_result)
h2 = runner.estimate_h2(dataset)

print("h2 =", float(h2.tot))
print("h2_se =", float(h2.tot_se))
print("intercept =", float(h2.intercept))
print("ldscore_file =", ldscore_result.output_paths["ldscore"])
print("weight_file =", ldscore_result.output_paths["w_ld"])
print("sumstats_config =", sumstats.config_snapshot)
print("dataset_config =", dataset.config_snapshot)
```

`SumstatsMunger` now captures `GlobalConfig` provenance into
`SumstatsTable.config_snapshot`, and `RegressionRunner.build_dataset()` checks
that the sumstats and LD-score results agree on critical settings such as
`snp_identifier` and `genome_build`. If they do not, the workflow raises
`ConfigMismatchError` instead of silently merging inconsistent artifacts.

## CLI

```bash
ldsc munge-sumstats \
  --sumstats "data/trait*.tsv.gz" \
  --snp SNP \
  --a1 A1 \
  --a2 A2 \
  --p P \
  --N-col N \
  --signed-sumstats BETA,0 \
  --out tutorial_outputs/trait

ldsc ldscore \
  --out tutorial_outputs/trait_ldscores \
  --baseline-annot "annotations/baseline.@.annot.gz" \
  --r2-table "r2/reference.@.parquet" \
  --frqfile "r2/reference_metadata.@.tsv.gz" \
  --snp-identifier rsid \
  --ld-wind-cm 1.0

ldsc h2 \
  --sumstats tutorial_outputs/trait.sumstats.gz \
  --ldscore tutorial_outputs/trait_ldscores.l2.ldscore.gz \
  --w-ld tutorial_outputs/trait_ldscores.w.l2.ldscore.gz \
  --counts tutorial_outputs/trait_ldscores.M_5_50 \
  --count-kind m_5_50 \
  --out tutorial_outputs/trait_h2
```
