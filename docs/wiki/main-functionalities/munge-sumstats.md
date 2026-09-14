# Munge summary statistics

Last updated on: 2026-09-14

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
  --trait-name mdd2025 \
  --output-dir "${OUTPUT_DIR}" \
  --infer-only
```

`--infer-only` inspects headers and sample rows and writes no output artifacts, logs, or directories; `--output-dir` is still required. It reports the detected file format, inferred column mappings, inferred source genome build, missing required fields, SNP restriction, selected liftover method, and a suggested command for the full run. It does not perform full filtering or mapping and therefore cannot report whole-run drop counts.

LDSC3 recognizes common aliases for required columns. If a column cannot be inferred correctly, specify it explicitly as described in [Override column-name inference](#override-column-name-inference).

### 2. Munge the summary statistics

```bash
ldsc munge-sumstats \
  --snp-identifier chr_pos \
  --source-genome-build auto \
  --output-genome-build hg19 \
  --raw-sumstats-file "${RAW_SUMSTATS_FILE}" \
  --trait-name mdd2025 \
  --output-dir "${OUTPUT_DIR}" \
  --overwrite
```

### Flags used in this command

- `--snp-identifier chr_pos` identifies SNPs by chromosome and base-pair position.
- `--source-genome-build auto` asks LDSC3 to infer whether the input coordinates use hg19 or hg38.
- `--output-genome-build hg19` writes SNP coordinates using hg19.
- `--raw-sumstats-file` specifies the raw GWAS summary-statistics file.
- `--trait-name mdd2025` stores the trait label in Parquet metadata and names the data file `mdd2025.parquet`. Unsafe filename characters are sanitized without changing the metadata label. Replace `mdd2025` with a descriptive name for your trait.
- `--output-dir` specifies the directory for the munged summary statistics and diagnostic artifacts.
- `--overwrite` permits existing artifacts in `--output-dir` to be replaced. Use it with caution; without this flag, the command returns an error rather than overwriting existing files.

### Defaults used by this command

The following options use their default values and are therefore omitted. Specify them only when you need to change the defaults:

- `--input-format auto` automatically detects the input file format.
- `--output-format parquet` writes the munged summary statistics in Parquet format.

Packaged HapMap3 restriction is also the default, and ordinary QC still applies. Use `--sumstats-snps-file FILE` to replace HM3 with a custom headered keep-list or `--no-snp-restriction` to disable keep-list filtering. These two overrides are mutually exclusive. No HM3 enable flag is needed.

The packaged HM3 reference includes `A1/A2`; allele-aware modes match its allele-aware keys as well as the base identifier. The example above explicitly chooses base `chr_pos`, which matches coordinates without alleles. See [HM3 filtering details](../../current/munge-sumstats.md#hm3-filter-and-quick-liftover).

## Genome-build conversion

Always choose `--output-genome-build` explicitly for coordinate-based identity. The package may infer the source build, but it never chooses the output build for you.

| Situation | Behavior |
| --- | --- |
| Source and requested output builds match | No liftover; a supplied chain is ignored |
| Builds differ, explicit `--liftover-chain-file FILE` supplied | Use that chain; disable quick liftover while retaining the selected SNP restriction |
| Builds differ, packaged HM3 restriction (default), no chain supplied | Automatically use quick liftover from package-bundled reference HM3 metadata |
| Builds differ, custom list or unrestricted SNPs, no chain supplied | Stop and require a chain file, even if the custom list contains only HM3 SNPs |
| Source build cannot be resolved | Stop and request `--source-genome-build hg19` or `hg38` |

To opt out of automatic quick liftover, add `--liftover-chain-file /path/to/sourceToOutput.over.chain` in the source-to-output direction. This changes the mapping method; use `--no-snp-restriction` separately if you also want to disable the HM3 keep-list. rsID-based identity rejects build and liftover options. See the [current liftover contract](../../current/munge-sumstats.md#liftover-rules).

## Outputs and run summary

The example writes `mdd2025.parquet`, `diagnostics/sumstats.log`, and `diagnostics/dropped_snps/dropped.tsv.gz`. Use `--output-format both` to additionally write `mdd2025.sumstats.gz`. Without `--trait-name`, the data filenames are `sumstats.parquet` and `sumstats.gz`; Parquet remains the default.

After success, stdout and the log show a `Munge-sumstats summary:` block at every log level. It names the restriction and selected mapping method and reports mapping input, mapped/retained, dropped rows, and drop reasons. Whole-run row totals and exclusive per-stage drops are separate: `sumstats_snps` counts keep-list removals, `liftover` counts mapping removals, and `identity` counts final identity cleanup. Mapping input excludes rows already removed by earlier QC or the keep-list. The dropped-SNP sidecar covers liftover and identity removals, not every QC stage. See [count interpretation](../../current/munge-sumstats.md#liftover-rules).

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

## Minimal raw input fields

For an ordinary headered, whitespace-delimited raw GWAS table, supply the fields below. Column names may use the aliases in the next table or explicit column flags. Each row represents one variant; matching the schema does not guarantee that the row survives QC or the default HM3 restriction.

| Field | When needed | Meaning and alternatives |
| --- | --- | --- |
| `SNP` | Required by the current raw-input reader in every identity mode. | A variant label. In rsID modes it supplies the matching identifier; coordinate modes still require this column even though matching uses coordinates. |
| `CHR`, `POS` | Required for `chr_pos` and `chr_pos_allele_aware`. | Chromosome and positive, one-based base-pair position in the source build. Optional for rsID modes. Choose the output build separately with `--output-genome-build`. |
| `A1`, `A2` | Required for allele-aware modes, including the default `chr_pos_allele_aware`. | `A1` is the allele relative to which the signed statistic is defined; `A2` is the other allele. Base `chr_pos` and `rsid` modes can munge without allele columns. |
| `P` | Required, including when input `Z` is available. | Association p-value in `(0, 1]`. LDSC derives the output Z magnitude from P; it does not compute P from a supplied effect or Z column. |
| One signed statistic: `Z`, `BETA`, `LOG_ODDS`, or `OR` | Required unless `--a1-inc` is explicitly used. | Supplies the direction relative to A1. The null is `0` for Z/beta/log odds and `1` for OR. Choose one explicitly with `--signed-sumstats COLUMN,NULL` if several are present. `--a1-inc` asserts that A1 is always the increasing allele and produces positive Z; use it only when that assertion is justified. |
| `N`, or both `N_CAS` and `N_CON` | Required unless sample size is supplied as a constant. | Use per-variant N or a case/control pair as described in [Sample-size columns](#sample-size-columns). Alternatively, use `--N VALUE` or both `--N-cas VALUE --N-con VALUE`. Input sample-size columns take precedence over these fallback constants. |

For the default allele-aware coordinate mode with per-variant N, a minimal header is:

```text
SNP CHR POS A1 A2 P BETA N
```

For explicit `--snp-identifier rsid`, a minimal header is `SNP P BETA N`; allele-aware rsID mode additionally needs A1/A2. `INFO`, `FRQ`, and `NSTUDY` are optional. Standard errors are not required because this workflow derives Z from P and the selected direction statistic. Already-munged legacy `SNP A1 A2 Z N` files can be passed directly to regression; the raw munging route still requires P. See [legacy summary-statistics compatibility](../../current/legacy-sumstats-compatibility.md).

Sources: `prepare_munge_input` in [_sumstats_input.py](../../../src/ldsc/_sumstats_input.py) validates required fields and sample-size choices; `munge_sumstats` in [_kernel/sumstats_munger.py](../../../src/ldsc/_kernel/sumstats_munger.py) converts P to signed Z.

## Automatically inferred column aliases

*<u>[Clean up this section.]</u>*

The table lists the raw-sumstats registry, including canonical spellings. Matching is case-insensitive; leading/trailing whitespace is stripped, and periods and hyphens are converted to underscores. For example, `p.value` matches `P_VALUE`. This reader uses cleaned exact aliases, not arbitrary substring or suffix matching. If multiple columns map to one field, choose the intended column explicitly or exclude the extra columns with `--ignore`.

| Canonical field | Recognized headers | Explicit override / role |
| --- | --- | --- |
| `SNP` | `SNP`, `MARKERNAME`, `MARKERID`, `SNPID`, `SNP_ID`, `RS`, `RSID`, `RS_ID`, `ID`, `RS_NUMBER`, `RS_NUMBERS`, `MARKER` | `--snp COLUMN` |
| `CHR` | `CHR`, `#CHROM`, `CHROM`, `CHROMOSOME` | `--chr COLUMN` |
| `POS` | `POS`, `BP`, `POSITION`, `BASE_PAIR`, `BASEPAIR` | `--pos COLUMN` |
| `NSTUDY` | `NSTUDY`, `N_STUDY`, `NSTUDIES`, `N_STUDIES` | `--nstudy COLUMN`; optional study-count filtering |
| `P` | `P`, `PVALUE`, `P_VALUE`, `PVAL`, `P_VAL`, `GC_PVALUE` | `--p COLUMN` |
| `A1` | `A1`, `ALLELE1`, `ALLELE_1`, `EFFECT_ALLELE`, `REFERENCE_ALLELE`, `INC_ALLELE`, `EA` | `--a1 COLUMN`; signed-statistic allele |
| `A2` | `A2`, `ALLELE2`, `ALLELE_2`, `OTHER_ALLELE`, `NON_EFFECT_ALLELE`, `DEC_ALLELE`, `NEA` | `--a2 COLUMN`; other allele |
| `N` | `N`, `WEIGHT` | `--N-col COLUMN` |
| `N_CAS` | `N_CAS`, `NCASE`, `CASES_N`, `N_CASE`, `N_CASES`, `NCAS`, `Nca` | `--N-cas-col COLUMN`; pair with control count |
| `N_CON` | `N_CON`, `N_CONTROLS`, `NCONTROL`, `CONTROLS_N`, `N_CONTROL`, `N_CTRL`, `NCON`, `Nco` | `--N-con-col COLUMN`; pair with case count |
| `INFO` | `INFO`, `IMPINFO` | `--info COLUMN`; optional imputation-quality filter |
| `FRQ` | `FRQ`, `EAF`, `MAF`, `FRQ_U`, `F_U` | `--frq COLUMN`; optional frequency/MAF filtering; selected values are preserved |
| `Z` | `Z`, `ZSCORE`, `Z-SCORE`, `GC_ZSCORE` | `--signed-sumstats COLUMN,0` |
| `OR` | `OR` | `--signed-sumstats COLUMN,1` |
| `BETA` | `BETA`, `B`, `EFFECTS`, `EFFECT` | `--signed-sumstats COLUMN,0` |
| `LOG_ODDS` | `LOG_ODDS` | `--signed-sumstats COLUMN,0` |
| `SIGNED_SUMSTAT` | `SIGNED_SUMSTAT` | Recognized field name, but its null is not inferred; specify `--signed-sumstats SIGNED_SUMSTAT,NULL` |

`NEFF` is deliberately not an automatic alias for N. Use `--N-col NEFF` only when its scientific meaning matches the intended sample-size input. A separate plain-format heuristic maps `REF` to A1 and `ALT` to A2 when both are present and none of `A1`, `A2`, `EA`, or `NEA` is present. Those hints are applied automatically unless overridden, although REF/ALT are not entries in the raw alias registry. Verify the effect orientation and use explicit `--a1`/`--a2` flags if the statistic refers to a different allele. Names such as `EFFECT_SIZE`, `LOGOR`, and `BETA_HAT` can produce suggestions rather than automatic mappings; supply `--signed-sumstats COLUMN,NULL` after checking their meaning. Frequency aliases do not establish allele orientation: a source field named `MAF` remains folded minor-allele frequency, not necessarily A1 frequency.

For DANER inputs, format-specific header handling additionally recognizes encoded count/frequency fields. The table describes ordinary raw-column aliases, not every DANER header convention. `--infer-only` previews mappings and suggestions; it does not replace the full run's required-field and QC checks.

Source: `RAW_SUMSTATS_REQUIRED_OR_OPTIONAL_SPECS`, `RAW_SUMSTATS_SIGNED_STAT_SPECS`, and `clean_header` in [column_inference.py](../../../src/ldsc/column_inference.py); `get_cname_map` in [_sumstats_input.py](../../../src/ldsc/_sumstats_input.py); `infer_raw_sumstats` in [sumstats_munger.py](../../../src/ldsc/sumstats_munger.py).

## Flags for specifying input fields and sample sizes

### Required fields and sample-size alternatives

| Field | Flag | Example |
| --- | --- | --- |
| SNP identifier | `--snp COLUMN` | `--snp variant_id` |
| Chromosome | `--chr COLUMN` | `--chr chromosome` |
| Position | `--pos COLUMN` | `--pos base_pair_location` |
| Effect allele, A1 | `--a1 COLUMN` | `--a1 effect_allele` |
| Other allele, A2 | `--a2 COLUMN` | `--a2 other_allele` |
| P-value | `--p COLUMN` | `--p p_value` |
| Signed effect statistic | `--signed-sumstats COLUMN,NULL` | `--signed-sumstats effect_size,0` |
| Per-variant sample size | `--N-col COLUMN` | `--N-col sample_size` |
| Case count, alternative to N | `--N-cas-col COLUMN` | `--N-cas-col cases` |
| Control count, paired with cases | `--N-con-col COLUMN` | `--N-con-col controls` |

Field requirements depend on identity mode and sample-size strategy; see [Minimal raw input fields](#minimal-raw-input-fields). For signed statistics, use null `0` for BETA, Z, or log odds; `1` for OR. For per-variant sample sizes, supply either `--N-col` or both case/control column flags.

**Constant sample sizes:** these flags take numbers instead of column names and supply the same sample size for every retained variant when per-variant sample-size columns are absent.

| Field | Flag | Example |
| --- | --- | --- |
| Constant total sample size | `--N VALUE` | `--N 100000` |
| Constant case count | `--N-cas VALUE` | `--N-cas 20000` |
| Constant control count | `--N-con VALUE` | `--N-con 80000` |

Supply either `--N` or both `--N-cas` and `--N-con`; the paired constants are added, so the example gives `N = 100000`. These are fallbacks: input N or paired case/control columns take precedence, followed by `--N`, then the paired constants. To use a constant instead of existing sample-size columns, exclude those columns with `--ignore` and omit their column-mapping flags. Constants are assigned after sample-size filtering and are not filtered by `--n-min`. Source: sample-size options in `build_parser` in [sumstats_munger.py](../../../src/ldsc/sumstats_munger.py) and `process_n` in [_kernel/sumstats_munger.py](../../../src/ldsc/_kernel/sumstats_munger.py).

### Optional fields and ignored columns

| Field | Flag | Example |
| --- | --- | --- |
| Allele frequency or MAF | `--frq COLUMN` | `--frq effect_allele_frequency` |
| Scalar imputation INFO | `--info COLUMN` | `--info imputation_quality` |
| INFO values stored as comma-separated lists within cells | `--info-list COLUMNS` | `--info-list per_study_info` |
| Number of contributing studies | `--nstudy COLUMN` | `--nstudy study_count` |
| Columns to exclude from inference and reading | `--ignore COL1,COL2` | `--ignore unused_beta,unused_frequency` |

Source: column options in `build_parser` in [sumstats_munger.py](../../../src/ldsc/sumstats_munger.py).
