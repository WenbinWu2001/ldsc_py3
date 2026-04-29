# Heritability Estimates

Goal: estimate SNP heritability for one trait with the refactored package, starting from raw summary statistics and ordinary unpartitioned LD scores built from an R2-table reference panel.

The examples below assume chromosome-pattern inputs such as `annotations/baseline.1.annot.gz`, `r2/reference.1.parquet`, and `r2/reference_metadata.1.tsv.gz`.
For parquet-backed LD scores, `r2/reference.@.parquet` is the canonical pair table
and `r2/reference_metadata.@.tsv.gz` is the required per-SNP sidecar. The sidecar
defines the reference SNP universe and supplies `MAF`/`CM`; the parquet is used
only for row-group-pruned LD window queries.

Path-token rules:

- scalar inputs such as raw sumstats must resolve to exactly one file
- group inputs may use globs or chromosome-suite tokens such as `baseline.@.annot.gz`
- output directories are literal destinations
- missing output directories are created and existing directories are reused
- existing fixed files are refused before writing starts unless you pass
  `--overwrite` or `overwrite=True`

Resolution behavior:

- there is no separate per-chromosome argument for LD-score inputs anymore; `baseline_annot_sources`, `r2_sources`, and `metadata_sources` all accept the unified token format
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
    snp_identifier="chr_pos",
    genome_build="auto",
)
set_global_config(GLOBAL_CONFIG)

sumstats = SumstatsMunger().run(
    MungeConfig(
        sumstats_file="data/trait*.tsv.gz",
        trait_name="trait",
        column_hints={
            "snp": "ID",
            "chr": "#CHROM",
            "pos": "POS",
            "a1": "EA",
            "a2": "NEA",
            "p": "PVAL",
            "N_col": "NEFF",
            "info": "IMPINFO",
        },
        # sumstats_snps_file="filters/hapmap3.tsv.gz",  # optional row keep-list
        output_dir="tutorial_outputs/trait",
        signed_sumstats_spec="BETA,0",
        # overwrite=True,  # enable only when intentionally replacing sumstats/log files
    ),
    global_config=GLOBAL_CONFIG,
)

# If you already have a curated .sumstats.gz artifact on disk, load it directly.
# Current artifacts recover config_snapshot from sumstats.metadata.json. Older
# files without that sidecar warn and use config_snapshot=None.
# sumstats = load_sumstats("tutorial_outputs/trait/sumstats.sumstats.gz", trait_name="trait")

ldscore_result = run_ldscore(
    output_dir="tutorial_outputs/trait_ldscores",
    r2_sources="r2/reference.@.parquet",
    metadata_sources="r2/reference_metadata.@.tsv.gz",
    r2_bias_mode="unbiased",
    regression_snps_file="filters/hapmap3.txt",
    common_maf_min=0.05,
    ld_wind_cm=1.0,
    # overwrite=True,  # enable only when intentionally replacing LD-score outputs
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

`SumstatsMunger` captures `GlobalConfig` provenance into `SumstatsTable.config_snapshot` and writes the same settings to `sumstats.metadata.json`. `RegressionRunner.build_dataset()` checks that known sumstats and LD-score snapshots agree on critical settings such as `snp_identifier` and `genome_build`. If both snapshots are known and incompatible, the workflow raises `ConfigMismatchError` instead of silently merging inconsistent artifacts. Sumstats loaded from current `.sumstats.gz` artifacts recover that snapshot from the sidecar; older files without the sidecar have `config_snapshot=None`, so regression treats their config provenance as unknown and validates the LD-score side only.

In `chr_pos` mode, regression merges sumstats and LD scores by normalized `CHR:POS`
coordinates. The raw munger accepts common coordinate headers such as `#CHROM`,
`CHROM`, `CHR`, `POS`, and `BP`; use `--chr` and `--pos` or `column_hints` when
the header is ambiguous. Leading raw `##` metadata lines are skipped before the
real header is parsed.

When you set `regression_snps_file` during LD-score calculation, the same regression SNP subset defines the rows of the normalized `baseline_table`. Regression weights come from the embedded `regr_weight` column by default.

Because this is ordinary unpartitioned heritability, `run_ldscore(...)` does
not need baseline annotations. With no baseline and no query inputs, it writes
a synthetic all-ones baseline column named `base`; `h2` reads that directory
through the same manifest path as any other LD-score result. The written
`baseline.parquet` file is flat, but its row groups are chromosome-aligned and
listed in `manifest.json` for chromosome-scoped inspection.

## CLI

The regression CLI reads the canonical LD-score result directory written by
`ldsc ldscore`.

```bash
ldsc munge-sumstats \
  --sumstats-file "data/trait*.tsv.gz" \
  --snp ID \
  --chr '#CHROM' \
  --pos POS \
  --a1 EA \
  --a2 NEA \
  --p PVAL \
  --N-col NEFF \
  --info IMPINFO \
  --sumstats-snps-file filters/hapmap3.tsv.gz \
  --signed-sumstats BETA,0 \
  --snp-identifier chr_pos \
  --genome-build auto \
  --output-dir tutorial_outputs/trait

ldsc ldscore \
  --output-dir tutorial_outputs/trait_ldscores \
  --r2-sources "r2/reference.@.parquet" \
  --metadata-sources "r2/reference_metadata.@.tsv.gz" \
  --r2-bias-mode unbiased \
  --regression-snps-file filters/hapmap3.txt \
  --snp-identifier chr_pos \
  --genome-build auto \
  --common-maf-min 0.05 \
  --ld-wind-cm 1.0

ldsc h2 \
  --sumstats-file tutorial_outputs/trait/sumstats.sumstats.gz \
  --ldscore-dir tutorial_outputs/trait_ldscores \
  --count-kind common \
  --output-dir tutorial_outputs/trait_h2
```

The command writes `tutorial_outputs/trait_h2/h2.tsv`.
If `sumstats.sumstats.gz`, `sumstats.log`, `sumstats.metadata.json`, or `h2.tsv`
already exists from a previous run, the relevant command fails before writing.
Add `--overwrite` only for an intentional rerun.
