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
columns when `genome_build` is set:

| Build | Accepted position aliases |
|---|---|
| `hg19` | `hg19_*`, `hg37_*`, or `grch37_*` with `pos`, `bp`, `position`, `base_pair`, or `basepair` |
| `hg38` | `hg38_*` or `grch38_*` with `pos`, `bp`, `position`, `base_pair`, or `basepair` |

### Raw summary statistics

Used by `src/ldsc/_kernel/sumstats_munger.py` through the same central
registry.

| Canonical field | Accepted aliases |
|---|---|
| `SNP` | `SNP`, `MARKERNAME`, `SNPID`, `SNP_ID`, `RS`, `RSID`, `RS_ID`, `ID`, `RS_NUMBER`, `RS_NUMBERS`, `MARKER` |
| `CHR` | `CHR`, `#CHROM`, `CHROM`, `CHROMOSOME` |
| `POS` | `POS`, `BP`, `POSITION`, `BASE_PAIR`, `BASEPAIR` |
| `P` | `P`, `PVALUE`, `P_VALUE`, `PVAL`, `P_VAL`, `GC_PVALUE` |
| `A1` | `A1`, `ALLELE1`, `ALLELE_1`, `EFFECT_ALLELE`, `REFERENCE_ALLELE`, `INC_ALLELE`, `EA` |
| `A2` | `A2`, `ALLELE2`, `ALLELE_2`, `OTHER_ALLELE`, `NON_EFFECT_ALLELE`, `DEC_ALLELE`, `NEA` |
| `N` | `N`, `WEIGHT` |
| `N_CAS` | `NCASE`, `CASES_N`, `N_CASE`, `N_CASES`, `N_CAS`, `NCAS` |
| `N_CON` | `N_CONTROLS`, `N_CON`, `NCONTROL`, `CONTROLS_N`, `N_CONTROL`, `NCON` |
| `INFO` | `INFO` |
| `FRQ` | `EAF`, `FRQ`, `MAF`, `FRQ_U`, `F_U` |
| `NSTUDY` | `NSTUDY`, `N_STUDY`, `NSTUDIES`, `N_STUDIES` |

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
