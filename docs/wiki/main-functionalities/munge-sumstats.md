# Munge summary statistics

Last updated on: 2026-08-10

`ldsc munge-sumstats` converts a raw GWAS table into LDSC3-ready summary statistics. See the current [munge-sumstats guide](../../current/munge-sumstats.md) for the full workflow and output contract.

## Recommended workflow

First run `--infer-only` to preview how LDSC3 interprets the input file. After reviewing the inferred configuration and any suggestions, run the command again without `--infer-only` to create the munged summary-statistics file.

### 1. Preview the inferred configuration

```bash
ldsc munge-sumstats \
  --snp-identifier chr_pos \
  --source-genome-build auto \
  --output-genome-build hg19 \
  --raw-sumstats-file "${RAW_SUMSTATS_FILE}" \
  --use-hm3-snps \
  --trait-name mdd2025 \
  --infer-only
```

`--infer-only` performs a dry run and writes no output artifacts. It reports the detected file format, inferred column mappings, inferred source genome build, missing required fields, whether liftover is needed, and a suggested command for the full run.

LDSC3 recognizes common aliases for required columns. If a column cannot be inferred correctly, specify it explicitly as described in [Override column-name inference](#override-column-name-inference).

### 2. Munge the summary statistics

```bash
ldsc munge-sumstats \
  --snp-identifier chr_pos \
  --source-genome-build auto \
  --output-genome-build hg19 \
  --raw-sumstats-file "${RAW_SUMSTATS_FILE}" \
  --use-hm3-snps \
  --use-hm3-quick-liftover \
  --trait-name mdd2025 \
  --output-dir "${OUTPUT_DIR}" \
  --overwrite
```

### Flags used in this command

- `--snp-identifier chr_pos` identifies SNPs by chromosome and base-pair position.
- `--source-genome-build auto` asks LDSC3 to infer whether the input coordinates use hg19 or hg38.
- `--output-genome-build hg19` writes SNP coordinates using hg19.
- `--raw-sumstats-file` specifies the raw GWAS summary-statistics file.
- `--use-hm3-snps` restricts the output to HapMap3 SNPs.
- `--use-hm3-quick-liftover` uses the packaged dual-build HapMap3 mapping to convert coordinates when the inferred source build differs from the output build. This option requires `--use-hm3-snps`; it is unnecessary when the input already uses hg19.
- `--trait-name mdd2025` assigns the trait name recorded in the output metadata. Replace `mdd2025` with a descriptive name for your trait.
- `--output-dir` specifies the directory for the munged summary statistics and diagnostic artifacts.
- `--overwrite` permits existing artifacts in `--output-dir` to be replaced. Use it with caution; without this flag, the command returns an error rather than overwriting existing files.

### Defaults used by this command

The following options use their default values and are therefore omitted. Specify them only when you need to change the defaults:

- `--format auto` automatically detects the input file format.
- `--output-format parquet` writes the munged summary statistics in Parquet format.

## Override column-name inference

If `--infer-only` does not identify the required columns correctly, use the legacy-compatible column flags to specify them explicitly. For example:

```bash
ldsc munge-sumstats \
  --snp MarkerID \
  --a1 Allele2 \
  --a2 Allele1 \
  --N-cas-col N_case \
  --N-con-col N_ctrl \
  ...
```

The case and control sample-size flags must be supplied together. If the file instead contains only an effective sample-size column such as `NEFF`, use `--N-col NEFF` when that quantity is scientifically appropriate for the intended analysis. Do not combine `--N-col` with `--N-cas-col` and `--N-con-col`; see the detailed sample-size guidance below.

## Sample-size columns

Use one per-variant sample-size strategy:

- Direct N: `--N-col <column>`.
- Case/control counts: `--N-cas-col <cases> --N-con-col <controls>`.

The case and control flags are a required pair, and they cannot be combined with `--N-col`. An explicit direct-N choice suppresses automatically inferred case/control columns with a warning; an explicit case/control choice suppresses an automatically inferred direct-N column with a warning.

### How the case/control pair becomes N

LDSC3 preserves the LDSC2 legacy case/control normalization. For variant `i`, let `T_i = NCAS_i + NCON_i` and `p_i = NCAS_i / T_i`. Let `p_ref` be the mean `p_i` among variants having the maximum `T_i`. The munger calculates:

```text
N_i = T_i * p_i / p_ref
```

This is **not** the commonly used harmonic effective sample size `4 / (1/NCAS_i + 1/NCON_i)`. If the case fraction is constant, the legacy calculation reduces to `NCAS_i + NCON_i`; if it varies, N is scaled relative to the case fraction among the maximum-total-count variants.

Choose the strategy according to the scientific meaning of the source fields:

- Use `--N-col NEFF` when the data producer documents `NEFF` as the effective sample size corresponding to the reported association statistic—for example, after accounting for case/control imbalance, per-variant missingness, meta-analysis participation, or analysis weights—and those exact values are intended for LDSC.
- Use `--N-cas-col NCAS --N-con-col NCON` when the per-variant counts are the trusted inputs and LDSC2-compatible legacy case/control normalization is the intended rule.

`NEFF` definitions vary between producers. Check the study documentation and compare `NEFF` with the count-derived values rather than assuming the two strategies are interchangeable.

If the header automatically maps both strategies, such as `N + NCAS + NCON`, LDSC3 stops rather than allowing one value to overwrite the other. Rerun with either:

```bash
--N-col N
```

or:

```bash
--N-cas-col NCAS --N-con-col NCON
```

For `NEFF + NCAS + NCON`, use `--N-col NEFF` when effective N is the intended quantity. The inferred `NCAS`/`NCON` columns are then suppressed automatically; `--ignore` is unnecessary.
