# Supported Column Alias Inference

Column aliases are centralized in `src/ldsc/column_inference.py`. The registry
is used by the main workflows (`munge-sumstats`, annotation loading,
`build-ref-panel`, `ldscore`, regression SNP restrictions, and package-written
artifacts) before data are renamed to canonical LDSC fields.

Alias matching is case-insensitive and punctuation-insensitive. For example,
`p-value`, `P_VALUE`, and `pvalue` all match the same alias family. Ambiguous
matches fail fast instead of guessing.

---

## Central registry

### Shared metadata columns

Used for annotation metadata, reference metadata, restriction files, and
coordinate inference.

| Canonical field | Accepted aliases |
|---|---|
| `SNP` | `SNP`, `SNPID`, `SNP_ID`, `RSID`, `RS_ID`, `RS`, `ID`, `MARKERNAME`, `MARKER` |
| `CHR` | `CHR`, `#CHROM`, `CHROM`, `CHROMOSOME` |
| `POS` | `POS`, `BP`, `POSITION`, `BASE_PAIR`, `BASEPAIR` |
| `CM` | `CM`, `CMBP`, `CENTIMORGAN` |
| `MAF` | `MAF`, `FRQ`, `FREQ`, `FREQUENCY` |

Restriction files in `chr_pos` mode also accept build-specific position
columns. At `build-ref-panel` time the source PLINK build chooses the column;
runtime LD-score and munging restrictions use the workflow's resolved
`GlobalConfig.genome_build`:

| Build | Accepted position aliases |
|---|---|
| `hg19` | `hg19_*`, `hg37_*`, or `grch37_*` with `pos`, `bp`, `position`, `base_pair`, or `basepair` |
| `hg38` | `hg38_*` or `grch38_*` with `pos`, `bp`, `position`, `base_pair`, or `basepair` |

### Raw summary statistics

Used by `src/ldsc/_kernel/sumstats_munger.py` through the same central
registry. The public `ldsc.sumstats_munger` workflow now runs a lightweight
header inference pass before delegating to the kernel. The default raw format
profile is `--format auto`, which detects plain whitespace text, old DANER, new
DANER, and PGC VCF-style headers. Users can inspect the decision without
writing outputs:

```bash
ldsc munge-sumstats --raw-sumstats-file raw.txt --infer-only
```

`--infer-only` reports the detected format, safe column hints, missing required
fields, INFO-list columns, notes, and a suggested minimal command. Normal runs
still use the same minimal form when headers are inferable:

```bash
ldsc munge-sumstats --raw-sumstats-file raw.txt --output-dir out
```

| Canonical field | Accepted aliases |
|---|---|
| `SNP` | `SNP`, `MARKERNAME`, `SNPID`, `SNP_ID`, `RS`, `RSID`, `RS_ID`, `ID`, `RS_NUMBER`, `RS_NUMBERS`, `MARKER` |
| `CHR` | `CHR`, `#CHROM`, `CHROM`, `CHROMOSOME` |
| `POS` | `POS`, `BP`, `POSITION`, `BASE_PAIR`, `BASEPAIR` |
| `P` | `P`, `PVALUE`, `P_VALUE`, `PVAL`, `P_VAL`, `GC_PVALUE` |
| `A1` | `A1`, `ALLELE1`, `ALLELE_1`, `EFFECT_ALLELE`, `REFERENCE_ALLELE`, `INC_ALLELE`, `EA` |
| `A2` | `A2`, `ALLELE2`, `ALLELE_2`, `OTHER_ALLELE`, `NON_EFFECT_ALLELE`, `DEC_ALLELE`, `NEA` |
| `N` | `N`, `WEIGHT` |
| `N_CAS` | `NCASE`, `CASES_N`, `N_CASE`, `N_CASES`, `N_CAS`, `NCAS`, `Nca` |
| `N_CON` | `N_CONTROLS`, `N_CON`, `NCONTROL`, `CONTROLS_N`, `N_CONTROL`, `NCON`, `Nco` |
| `INFO` | `INFO`, `IMPINFO` |
| `FRQ` | `EAF`, `FRQ`, `MAF`, `FRQ_U`, `F_U` |
| `NSTUDY` | `NSTUDY`, `N_STUDY`, `NSTUDIES`, `N_STUDIES` |

`A1` and `A2` are legacy LDSC output column names and remain the canonical
written schema for compatibility. Their semantics are:

- `A1` is the allele that the signed statistic is relative to; in practice this
  is usually the effect or increasing allele.
- `A2` is the counterpart allele.
- Positive `Z`, positive `BETA`, positive `LOG_ODDS`, and `OR > 1` mean the
  trait or risk is increasing with respect to `A1`.
- Alias names containing "reference" mean reference for the signed summary
  statistic, not the genome reference allele.

Do not add genome-oriented allele aliases to the global A1/A2 registry unless
their signed-statistic semantics are known. `REF`/`ALT` are handled only through
format-aware inference, currently for PGC VCF-style inputs. `ALT_ALLELE` is not
an accepted alias. `BEAA` is not an accepted alias because it is a privacy-masked
form of `BETA`, not a public header convention.

`NEFF` is intentionally not an alias for `N`. Effective sample size and total
sample size are different quantities. If a user decides a particular analysis
should use `NEFF` as the munger's `N`, they must pass `--N-col NEFF`
explicitly.

Signed statistic aliases:

| Canonical field | Accepted aliases |
|---|---|
| `Z` | `ZSCORE`, `Z-SCORE`, `GC_ZSCORE`, `Z` |
| `OR` | `OR` |
| `BETA` | `B`, `BETA`, `EFFECTS`, `EFFECT` |
| `LOG_ODDS` | `LOG_ODDS` |
| `SIGNED_SUMSTAT` | `SIGNED_SUMSTAT` |

Explicit column hints still take priority over inferred raw-sumstats aliases.
The useful hints are `--snp`, `--chr`, `--pos`, `--a1`, `--a2`, `--p`, and
signed-statistic options such as `--signed-sumstats`.

INFO handling accepts ordinary scalar values such as `IMPINFO=0.97`. If an
INFO-like column contains comma-separated per-study values with numeric and
missing tokens, for example `IMPINFO=0.852,0.113,0.842,0.88,NA`, the munger
can treat it as `--info-list IMPINFO` and filter on the mean of the numeric
non-missing values. A mixed list such as `IMPINFO=0.95,LOW,0.88` is rejected
instead of coerced; the error suggests either `--ignore IMPINFO` or an explicit
`--info-list` only when the list is genuinely numeric/NA.

Error messages are expected to suggest exact repair flags when likely columns
exist: missing alleles with `REF`/`ALT` suggests `--a1 REF --a2 ALT` with the
signed-statistic caveat; missing signed statistics with likely effect columns
suggests `--signed-sumstats <col>,0`; missing `N` with `NEFF` explains that
`NEFF` is not inferred automatically.

### Parquet R2 input

Canonical parquet R2 files are read as the logical fields
`CHR`, `POS_1`, `POS_2`, `SNP_1`, `SNP_2`, and `R2`. Accepted load-time aliases:

| Canonical field | Accepted aliases |
|---|---|
| `CHR` | `CHR`, `#CHROM`, `CHROM`, `CHROMOSOME` |
| `POS_1` | `POS_1`, `POS1`, `BP_1`, `BP1`, `POSITION_1`, `POSITION1` |
| `POS_2` | `POS_2`, `POS2`, `BP_2`, `BP2`, `POSITION_2`, `POSITION2` |
| `SNP_1` | numbered SNP aliases such as `SNP_1`, `SNP1`, `RSID_1`, `RSID1`, `ID_1`, `ID1` |
| `SNP_2` | numbered SNP aliases such as `SNP_2`, `SNP2`, `RSID_2`, `RSID2`, `ID_2`, `ID2` |
| `R2` | `R2` |

Legacy raw R2 parquet files also accept build-specific pair columns such as
`hg19_pos_1`, `hg19_pos_2`, `hg38_pos_1`, `hg38_pos_2`, their `bp`/`position`
variants, numbered unique-ID columns, `Dprime` / `D_PRIME`, and
`+/-corr` / `CORR` / `SIGNCORR` / `SIGN_CORR`.

### Internal artifacts

Package-written artifacts are intentionally strict. Curated `.sumstats`,
annotation, and LD-score artifacts are expected to use canonical field names
only (`SNP`, `CHR`, `POS`, `N`, `Z`, `A1`, `A2`, `FRQ`, `CM`, `MAF`, as
applicable). This catches schema drift in files produced by the package.

### Non-column token aliases

The same registry also normalizes these user-facing tokens:

| Setting | Accepted aliases |
|---|---|
| `snp_identifier=rsid` | `rsid`, `rs_id`, `rs`, `snp`, `snpid`, `snp_id` |
| `snp_identifier=chr_pos` | `chr_pos`, `chrpos`, `chrom_pos`, `chromosome_position`, `position` |
| `genome_build=hg19` | `hg19`, `hg37`, `grch37` |
| `genome_build=hg38` | `hg38`, `grch38` |
| `genome_build=auto` | `auto` |

---

## Functionality-specific inference

Most aliases should be added to `src/ldsc/column_inference.py`. The current
local exception is genetic-map loading in
`src/ldsc/_kernel/ref_panel_builder.py`, where the centiMorgan column accepts:

`CM`, `GENETIC_MAP_CM`, `GENETICMAPCM`, `GENETIC_MAP(CM)`,
`GENETICMAP(CM)`.

To find any future local exceptions, search for `ColumnSpec(` outside
`column_inference.py`.
