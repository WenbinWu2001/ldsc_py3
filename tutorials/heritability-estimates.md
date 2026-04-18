# Heritability Estimates

Goal: estimate SNP heritability for one trait with the refactored package, starting from raw summary statistics and an R2-table reference panel.

The examples below assume chromosome-pattern inputs such as `annotations/baseline.1.annot.gz`, `r2/reference.1.parquet`, and `r2/reference_metadata.1.tsv.gz`.

Path-token rules:

- scalar inputs such as raw sumstats must resolve to exactly one file
- group inputs may use globs or chromosome-suite tokens such as `baseline.@`
- output prefixes are literal destinations

## Python API

```python
from ldsc import CommonConfig, MungeConfig, RawSumstatsSpec, RegressionConfig, RegressionRunner, SumstatsMunger, run_ldscore

common = CommonConfig(snp_identifier="rsid")

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
    common,
)

ldscore_result = run_ldscore(
    out="tutorial_outputs/trait_ldscores",
    baseline_annot="annotations/baseline.@",
    r2_table="r2/reference.@",
    frqfile="r2/reference_metadata.@",
    snp_identifier="rsid",
    ld_wind_cm=1.0,
)

runner = RegressionRunner(common, RegressionConfig())
dataset = runner.build_dataset(sumstats, ldscore_result)
h2 = runner.estimate_h2(dataset)

print("h2 =", float(h2.tot))
print("h2_se =", float(h2.tot_se))
print("intercept =", float(h2.intercept))
print("ldscore_file =", ldscore_result.output_paths["ldscore"])
print("weight_file =", ldscore_result.output_paths["w_ld"])
```

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
  --baseline-annot "annotations/baseline.@" \
  --r2-table "r2/reference.@" \
  --frqfile "r2/reference_metadata.@" \
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
