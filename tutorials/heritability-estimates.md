# Heritability Estimates

Goal: estimate SNP heritability for one trait with the refactored package, starting from raw summary statistics and ordinary unpartitioned LD scores built from an R2-table reference panel.

The examples below assume chromosome-pattern annotation inputs such as
`annotations/baseline.1.annot.gz` and a package-built R2 directory such as
`r2_ref_panel_1kg30x_1cM_hm3/hg38`. For parquet-backed LD scores, each
`chr*_r2.parquet` file is the canonical pair table and the matching
`chr*_meta.tsv.gz` sidecar is optional but strongly recommended because it
defines the complete reference SNP universe and supplies `MAF`/`CM`.
Package-built R2 parquet files carry `ldsc:r2_bias` and `ldsc:n_samples` in
schema metadata, so the examples omit R2 bias and sample-size arguments.

Path-token rules:

- scalar inputs such as raw sumstats must resolve to exactly one file
- group inputs may use globs or chromosome-suite tokens such as `baseline.@.annot.gz`
- output directories are literal destinations
- missing output directories are created and existing directories are reused
- existing owned workflow artifacts are refused before writing starts unless
  you pass `--overwrite` or `overwrite=True`; successful overwrites remove
  stale owned siblings that the current run did not produce and preserve
  unrelated files

Resolution behavior:

- there is no separate per-chromosome argument for LD-score annotation inputs anymore; `baseline_annot_sources` accepts the unified token format, and the parquet R2 reference panel is supplied as one build directory with `r2_dir`
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
    snp_identifier="chr_pos_allele_aware",
    genome_build="auto",
)
set_global_config(GLOBAL_CONFIG)

sumstats = SumstatsMunger().run(
    MungeConfig(
        raw_sumstats_file="data/trait*.tsv.gz",
        trait_name="trait",
        # Most common columns are inferred automatically. If this file has
        # NEFF and you want to use it as the munger's N, pass this explicitly:
        # column_hints={"N_col": "NEFF"},
        # use_hm3_snps=True,  # packaged HM3 row restriction
        # target_genome_build="hg38",
        # use_hm3_quick_liftover=True,  # requires use_hm3_snps=True
        output_dir="tutorial_outputs/trait",
        # overwrite=True,  # also removes stale unselected sumstats sibling formats
    ),
    global_config=GLOBAL_CONFIG,
)

# If you already have a curated sumstats artifact on disk, load it directly.
# Current parquet and .sumstats.gz artifacts recover config_snapshot from
# sumstats.metadata.json. Old package-written files without current provenance
# must be regenerated with the current LDSC package.
# sumstats = load_sumstats("tutorial_outputs/trait/sumstats.parquet", trait_name="trait")

ldscore_result = run_ldscore(
    output_dir="tutorial_outputs/trait_ldscores",
    r2_dir="r2_ref_panel_1kg30x_1cM_hm3/hg38",
    use_hm3_regression_snps=True,
    common_maf_min=0.05,
    ld_wind_cm=1.0,
    # overwrite=True,  # also removes stale ldscore.query.parquet if this run is baseline-only
)

runner = RegressionRunner(global_config=GLOBAL_CONFIG, regression_config=RegressionConfig())
dataset = runner.build_dataset(sumstats, ldscore_result)
h2 = runner.estimate_h2(dataset)

print("h2 =", float(h2.tot))
print("h2_se =", float(h2.tot_se))
print("intercept =", float(h2.intercept))
print("ldscore_rows =", len(ldscore_result.baseline_table))
print("embedded_regression_ld_scores_column =", "regression_ld_scores" in ldscore_result.baseline_table.columns)
print("sumstats_config =", sumstats.config_snapshot)
print("dataset_config =", dataset.config_snapshot)
```

`SumstatsMunger` captures `GlobalConfig` provenance into `SumstatsTable.config_snapshot` and writes the same settings to `sumstats.metadata.json`. If you pass `trait_name=` in Python or `--trait-name` in the CLI, the sidecar also stores that biological label for downstream result tables. `RegressionRunner.build_dataset()` checks that current sumstats and LD-score snapshots agree on critical settings such as `snp_identifier` and `genome_build`. If current snapshots are incompatible, the workflow raises `ConfigMismatchError` instead of silently merging inconsistent artifacts. Sumstats loaded from current `sumstats.parquet` or `.sumstats.gz` artifacts recover that snapshot and trait label from the sidecar; old package-written files without current provenance must be regenerated with the current LDSC package.

`SumstatsMunger.run()` is also the shared implementation path behind
`ldsc munge-sumstats`: the CLI parser maps arguments into `MungeConfig`, then
the workflow writes `sumstats.parquet` by default, optional
`sumstats.sumstats.gz` compatibility output when requested, `sumstats.log`,
`sumstats.metadata.json`, and `dropped_snps/dropped.tsv.gz` under the selected
output directory. The log is an audit file; `MungeRunSummary.output_paths` lists
curated sumstats artifacts and the dropped-SNP audit sidecar, not
`sumstats.log`.

In the default `chr_pos_allele_aware` mode, regression merges sumstats and LD
scores by normalized `CHR:POS:<allele_set>` identity. `A1/A2` are required and
are used only to make the merge key safer; the `SNP` column is a label and may
contain rsIDs even when coordinate mode is active. To run coordinate identity
without allele-aware matching, set `snp_identifier="chr_pos"` or pass
`--snp-identifier chr_pos`; for rsID-only identity, use `snp_identifier="rsid"`
or `--snp-identifier rsid`. The removed `--no-alleles` flag is not accepted.
The raw munger accepts common coordinate headers
such as `#CHROM`, `CHROM`, `CHR`, `POS`, and `BP`; use `--chr` and `--pos` or
`column_hints` when the header is ambiguous. Leading raw `##` metadata lines are
skipped before the real header is parsed. `--format auto` is the default and
detects plain text, old DANER, new DANER, and PGC VCF-style headers. Use
`ldsc munge-sumstats --raw-sumstats-file raw.txt --infer-only` to inspect
format, column, INFO-list, and missing-field decisions without writing output.
`A1` is the allele that the signed statistic is relative to, not necessarily
the genome reference allele; `NEFF` is not inferred as `N` unless you opt in
with `--N-col NEFF` or `column_hints={"N_col": "NEFF"}`. Optional munger liftover runs after
the source-build SNP filter, drops duplicate source/target coordinate
groups, changes only `CHR`/`POS`, and requires
`--target-genome-build` plus either `--liftover-chain-file` or
`--use-hm3-snps --use-hm3-quick-liftover`. Drop counts are written to `sumstats.log`,
examples appear only at `DEBUG`, and row-level drops are audited in
`dropped_snps/dropped.tsv.gz`; the metadata sidecar stays limited to current
artifact provenance.

When you set `regression_snps_file` or `use_hm3_regression_snps` during LD-score calculation, the same regression SNP subset defines the rows of the normalized `baseline_table`. Regression uses the embedded `regression_ld_scores` column as the historical `w_ld` component; final model-dependent weights are computed later in the regression kernel.

Because this is ordinary unpartitioned heritability, `run_ldscore(...)` does
not need baseline annotations. With no baseline and no query inputs, it writes
a synthetic all-ones baseline column named `base`; `h2` reads that directory
through the same manifest path as any other LD-score result. The written
`ldscore.baseline.parquet` file is flat, but its row groups are
chromosome-aligned and listed in `manifest.json` for chromosome-scoped
inspection.

## CLI

The regression CLI reads the canonical LD-score result directory written by
`ldsc ldscore`.

```bash
ldsc munge-sumstats \
  --raw-sumstats-file "data/trait*.tsv.gz" \
  --trait-name trait \
  --use-hm3-snps \
  --snp-identifier chr_pos \
  --genome-build auto \
  --output-dir tutorial_outputs/trait

# Add explicit repair flags only when --infer-only reports that they are needed,
# for example --N-col NEFF if that is appropriate for the analysis.

# Optional chr_pos-family liftover when the resolved source build differs:
#   --target-genome-build hg38 \
#   --use-hm3-snps \
#   --use-hm3-quick-liftover

ldsc ldscore \
  --output-dir tutorial_outputs/trait_ldscores \
  --r2-dir "r2_ref_panel_1kg30x_1cM_hm3/hg38" \
  --use-hm3-regression-snps \
  --snp-identifier chr_pos \
  --genome-build auto \
  --common-maf-min 0.05 \
  --ld-wind-cm 1.0

ldsc h2 \
  --sumstats-file tutorial_outputs/trait/sumstats.parquet \
  --ldscore-dir tutorial_outputs/trait_ldscores \
  --count-kind common \
  --output-dir tutorial_outputs/trait_h2
```

The command writes `tutorial_outputs/trait_h2/h2.tsv` and
`tutorial_outputs/trait_h2/h2.log`.
If any owned output-family file already exists from a previous run, the
relevant command fails before writing. For `munge-sumstats`, that family is
`sumstats.parquet`, `sumstats.sumstats.gz`, `sumstats.metadata.json`,
`dropped_snps/dropped.tsv.gz`, and `sumstats.log`, even if the current
`--output-format` would write only one data format. Add `--overwrite` only for an intentional rerun; when used, stale
sumstats sibling formats not produced by the successful run are removed.
