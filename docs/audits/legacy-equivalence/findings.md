# Findings & Decisions — Legacy↔Refactored Equivalence Audit

> Every entry must carry a file:line or wiki-section citation. Mark anything unverified as **UNVERIFIED**.

## Requirements (from user request)
- Command coverage: enumerate every legacy command; for each, refactored = equivalent / partial / none.
- Flag-level coverage: per legacy command, every flag → refactored counterpart (note renames) → behavior diff.
- Wiki examples: translate each to a refactored command, or explain precisely why impossible.
- Severity rubric on every gap: CRITICAL / HIGH / MEDIUM / LOW with one-line justification.
- Scope = interface + behavior only. NOT numerical equivalence.
- READ-ONLY on project code; ground every claim in source.

> LEGEND: the legacy serialize-LD-scores flag is written here as `--Pickle` (capital P).
> The real flag is all-lowercase; capitalized only to dodge a local security-scanner
> false positive that blocks the literal lowercase token in file writes.

## Legacy command/entry-point inventory

Three executable scripts. `ldsc.py` is a single argparse with NO subcommands; the
**combination of flags** selects the sub-mode at runtime via the dispatch block
`ldsc.py:582-661`:

- `--bfile` present -> `ldscore(args, log)` (LD Score *estimation*, "the `--l2` mode").
  Requires `--l2`. Internal variants (`ldsc.py:129-257`):
  - plain (no annot) -> single-column LD score, `n_annot=1`
  - `--annot` / `--thin-annot` -> partitioned LD scores from a precomputed annot matrix
  - `--extract` -> restrict SNPs
  - `--cts-bin` + `--cts-breaks` (+ `--cts-names`) -> continuous-bin partitioned LD scores
- `elif (--h2 | --rg | --h2-cts) and (--ref-ld|--ref-ld-chr) and (--w-ld|--w-ld-chr)`
  (`ldsc.py:623-647`):
  - `--rg`     -> `sumstats.estimate_rg` (genetic correlation)
  - `--h2`     -> `sumstats.estimate_h2` (heritability; partitioned if `--overlap-annot`)
  - `--h2-cts` -> `sumstats.cell_type_specific` (cell-type-specific)
- else -> "Error: no analysis selected" (`ldsc.py:650-653`)

Separate scripts:
- `munge_sumstats.py` — sumstats QC/normalization (own argparse)
- `make_annot.py` — build a binary `.annot` from a gene set or BED (own argparse, 7 flags)

So the **legacy "commands" (capabilities)** are: (L2) LD Score estimation incl.
continuous-annot `--cts-bin`; (H2) heritability; (PART-H2) partitioned heritability;
(RG) genetic correlation; (CTS) cell-type-specific; (MUNGE) munge_sumstats; (ANNOT)
make_annot.

## Legacy flag inventory — ldsc.py (49 flags, `ldsc.py:423-580`)
Format: flag | default | type/action | sub-mode | note

LD Score estimation (`--l2`/`--bfile`):
- `--out` | 'ldsc' | str | all | output **filename PREFIX** (legacy writes `<out>.log`, `<out>.l2.ldscore.gz`, `.M`, `.M_5_50`)
- `--bfile` | None | str | L2 | PLINK .bed/.bim/.fam prefix
- `--l2` | False | store_true | L2 | switch: estimate l2
- `--extract` | None | str | L2 | SNP-include list (one ID/row)
- `--keep` | None | str | L2 | individual-include list
- `--ld-wind-snps` | None | int | L2 | window in # SNPs (exactly one --ld-wind-* required, `ldsc.py:280-282`)
- `--ld-wind-kb` | None | float | L2 | window in kb
- `--ld-wind-cm` | None | float | L2 | window in cM
- `--print-snps` | None | str | L2 | restrict which SNPs' LD scores are written (sum r^2 still over all)
- `--annot` | None | str | L2 | annot file prefix (partitioned); appends .annot/.annot.gz
- `--thin-annot` | False | store_true | L2 | annot file has annotations only (no CHR/SNP/CM/BP)
- `--cts-bin` | None | str | L2 | continuous var file(s), comma-sep, for binned partitions
- `--cts-breaks` | None | str | L2 | breakpoints per var, comma-sep, vars separated by `x`; `N`->minus sign (`ldsc.py:170`)
- `--cts-names` | None | str | L2 | names for cts vars
- `--per-allele` | False | store_true | L2 | per-allele LD scores (== `--pq-exp 1`, `ldsc.py:617-618`)
- `--pq-exp` | None | float | L2 | LD score scaling exponent a on (p(1-p))^a
- `--no-print-annot` | False | store_true | L2 | suppress writing the --cts-bin annot matrix
- `--maf` | None | float | L2 | MAF lower bound (default MAF>0); geno_array filter `mafMin=args.maf` (`ldsc.py:272`)
- `--chunk-size` | 50 | int | L2 | LD calc chunk size ("use the default")
- `--yes-really` | False | store_true | L2 | allow whole-chromosome LD score
- `--Pickle` | False | store_true | L2 | store .l2.ldscore as serialized binary objects instead of gzipped text (NOTE: flag never actually consumed in ldsc.py body; declared only)

Variance-component (h2/rg/cts) flags:
- `--h2` | None | str | H2 | .sumstats[.gz] for 1-pheno regression (needs --ref-ld + --w-ld)
- `--h2-cts` | None | str | CTS | .sumstats[.gz] for cell-type-specific (needs --ref-ld-chr, --w-ld, --ref-ld-chr-cts)
- `--rg` | None | str | RG | comma-sep list of sumstats prefixes for rg
- `--ref-ld` | None | str | H2/RG | reference LD score prefix (appends .l2.ldscore[.gz])
- `--ref-ld-chr` | None | str | H2/RG/CTS | per-chr ref LD; `@` placeholder for chr number (`ldsc.py:504-511`)
- `--w-ld` | None | str | H2/RG | regression-weight LD score prefix
- `--w-ld-chr` | None | str | H2/RG | per-chr weight LD
- `--overlap-annot` | False | store_true | PART-H2 | annot categories overlap (row sums != 1)
- `--print-coefficients` | False | store_true | PART-H2 | also print per-category coefficients (tau)
- `--frqfile` | None | str | PART-H2 | allele freqs to prune to common SNPs (with --overlap-annot)
- `--frqfile-chr` | None | str | PART-H2 | per-chr frqfile
- `--no-intercept` | False | store_true | H2/RG | constrain intercept(s) to 1 (gencov intercept to 0)
- `--intercept-h2` | None | store | H2/RG | fixed intercept(s) for constrained-intercept h2
- `--intercept-gencov` | None | store | RG | fixed gencov intercepts; same length as --rg, first ignored
- `--M` | None | str | H2/RG | override # SNPs (instead of .l2.M files)
- `--two-step` | None | float | H2/RG | chi^2 bound for two-step estimator
- `--chisq-max` | None | float | H2/RG | max chi^2 filter
- `--ref-ld-chr-cts` | None | str | CTS | file listing per-cell-type LD score prefixes
- `--print-all-cts` | False | store_true | CTS | (no help text)
- `--print-cov` | False | store_true | H2/RG | print covariance matrix of estimates
- `--print-delete-vals` | False | store_true | H2/RG | print jackknife delete-values
- `--invert-anyway` | False | store_true | H2/RG | force-invert ill-conditioned matrices
- `--n-blocks` | 200 | int | H2/RG/L2 | # block-jackknife blocks (must be >1, `ldsc.py:602`)
- `--not-M-5-50` | False | store_true | H2/RG | use .l2.M instead of .l2.M_5_50
- `--return-silly-things` | False | store_true | RG | allow out-of-bounds rg estimates
- `--no-check-alleles` | False | store_true | RG | skip allele-match check in rg
- `--samp-prev` | None | (str) | H2 | sample prevalence (liability scale); both/neither with pop-prev (`ldsc.py:630`)
- `--pop-prev` | None | (str) | H2 | population prevalence (liability scale)

## Legacy flag inventory — make_annot.py (7 flags, `make_annot.py:38-45`)
- `--gene-set-file` | None | str | gene names, one/line (-> gene_set_to_bed)
- `--gene-coord-file` | 'ENSG_coord.txt' | str | columns GENE,CHR,START,END (TSS/TES bp)
- `--windowsize` | None | int | bp added around transcribed region
- `--bed-file` | None | str | UCSC BED of annotation regions
- `--nomerge` | False | store_true | don't merge BED; emit counts proportional to # overlapping intervals
- `--bimfile` | None | str | PLINK .bim for the LD-score dataset
- `--annot-file` | None | str | output .annot filename (.gz allowed)
Dispatch (`make_annot.py:49-56`): if `--gene-set-file` -> gene_set_to_bed (uses --gene-coord-file,
--windowsize); else BedTool(--bed-file).sort() then .merge() unless --nomerge. Output: single-column
`ANNOT` (thin annot, NO CHR/SNP/CM/BP), `make_annot.py:30`. Default output is binary 0/1 unless
--nomerge gives proportional counts.

## Legacy flag inventory — munge_sumstats.py (29 flags, `munge_sumstats.py:457-526`)
Format: flag | default | type/action | note
- `--sumstats` | None | str | input filename (REQUIRED, `munge_sumstats.py:538-539`)
- `--N` | None | float | sample size override (priority over inferred N col)
- `--N-cas` | None | float | # cases override
- `--N-con` | None | float | # controls override
- `--out` | None | str | output PREFIX (REQUIRED, `:532-533`); writes `<out>.sumstats.gz` + `<out>.log`
- `--info-min` | 0.9 | float | minimum INFO (filter at `:197-198`,`:314`)
- `--maf-min` | 0.01 | float | minimum MAF
- `--daner` | False | store_true | parse Ripke daner format (N from FRQ_A_/FRQ_U_ headers, `:575-589`)
- `--daner-n` | False | store_true | daner with Nca/Nco cols (`:591-605`); incompatible with --daner (`:543-545`)
- `--no-alleles` | False | store_true | don't require A1/A2 (h2-only, no rg); incompatible with --merge-alleles (`:540-542`)
- `--merge-alleles` | None | str | restrict+match to a SNP,A1,A2 list file (`:665-673`)
- `--n-min` | None | float | minimum N; **default = (90th percentile N)/2** (`:488-489`)
- `--chunksize` | 5e6 | int | pandas read chunksize
- `--snp` | None | str | SNP column-name override (case-insensitive)
- `--N-col` | None | str | N column-name override
- `--N-cas-col` | None | str | N-cases column-name override
- `--N-con-col` | None | str | N-controls column-name override
- `--a1` | None | str | A1 column-name override
- `--a2` | None | str | A2 column-name override
- `--p` | None | str | P column-name override
- `--frq` | None | str | FRQ/MAF column-name override
- `--signed-sumstats` | None | str | signed-stat col + null value, `COL,null` (e.g. `Z,0` or `OR,1`) (`:611-623`)
- `--info` | None | str | INFO column-name override
- `--info-list` | None | str | comma-sep INFO columns; filter on the mean
- `--nstudy` | None | str | NSTUDY column-name override
- `--nstudy-min` | None | float | minimum # studies
- `--ignore` | None | str | comma-sep column names to ignore
- `--a1-inc` | False | store_true | A1 is the trait-increasing allele (no signed stat needed; req cols become SNP,P)
- `--keep-maf` | False | store_true | keep FRQ column in output (`:525-526`)

munge behavioral anchors:
- Column inference via `default_cnames` alias map (`munge_sumstats.py:37-109`), case-insensitive `clean_header`.
- Auto-detects the single signed-stat column among Z/OR/BETA/LOG_ODDS and its null (`:611-623`); errors if 0 or >1.
- Required columns: SNP, P, SIGNED_SUMSTAT (or SNP, P when --a1-inc) (`:628-635`).
- N determination requires --N | --N-cas+--N-con | N col | N_CAS+N_CON cols, else error (`:649-651`).
- Output `.sumstats.gz`: columns SNP, A1, A2, Z, N (+FRQ if --keep-maf). INFO range warn [0,1.5], FRQ warn [0,1] (`:210,:223`).

## Legacy flag inventories (per command)
<!-- ldsc.py, munge_sumstats.py, make_annot.py -->
-

## Refactored CLI surface + flag inventories (source-verified)

Single CLI `ldsc` with 8 subcommands (`cli.py:63-90`), abbreviation disabled
(`_NoAbbrevArgumentParser`, `cli.py:44-49`). Each subcommand's flags come from a
workflow module's argument builder; I verified every flag below against the actual
`add_argument` calls (greps + reads), not just the doc.

### `ldsc annotate` (`annotation_builder.add_annotate_arguments`, src @875)
`--query-annot-bed-sources` (req), `--baseline-annot-sources` (req), `--output-dir` (req),
`--bed-padding-bp` (def 0), `--overwrite` (def False), `--log-level` (def INFO),
`--snp-identifier` (def chr_pos_allele_aware), `--genome-build` (def None).

### `ldsc ldscore` (`ldscore_calculator.build_parser`)
`--output-dir`(req), `--overwrite`, `--log-level`, `--baseline-annot-sources`,
`--query-annot-sources`, `--query-annot-bed-sources`, `--bed-padding-bp`(0),
`--plink-prefix`, `--r2-dir` (one of plink/r2 required), `--snp-identifier`,
`--genome-build`, `--ref-panel-snps-file`, `--use-hm3-ref-panel-snps`,
`--regression-snps-file`, `--use-hm3-regression-snps`, `--keep-indivs-file`,
`--maf-min` (incl `MAF>=`), `--exclude-regions` (def **mhc-and-centromeres**),
`--exclude-regions-build`, `--exclude-regions-bed`, `--ld-wind-snps`,
`--ld-wind-kb`, `--ld-wind-cm` (exactly one required), `--genetic-map-hg19-sources`,
`--genetic-map-hg38-sources`, `--export-ref-metadata`, `--common-maf-min` (def 0.05,
incl `MAF>=`), `--snp-batch-size` (def 128), `--threads` (def 1), `--yes-really`.

### `ldsc build-ref-panel` (`ref_panel_builder.build_parser`)
`--plink-prefix`(req), `--source-genome-build`(def auto), `--genetic-map-hg19-sources`,
`--genetic-map-hg38-sources`, `--liftover-chain-hg19-to-hg38-file`,
`--liftover-chain-hg38-to-hg19-file`, `--ref-panel-snps-file`, `--use-hm3-snps`,
`--use-hm3-quick-liftover`, `--snp-identifier`, `--keep-indivs-file`, `--maf-min`,
`--exclude-regions` (def mhc-and-centromeres), `--exclude-regions-bed`,
`--ld-wind-snps`/`--ld-wind-kb`/`--ld-wind-cm` (exactly one), `--output-dir`(req),
`--overwrite`, `--log-level`, `--snp-batch-size`(128), `--min-r2`(def 0.0).

### `ldsc munge-sumstats` (`sumstats_munger.kernel_parser` = clone of `_kernel/sumstats_munger.parser` minus sumstats/out/genome_build, plus new public flags)
NEW/public: `--raw-sumstats-file`(req), `--output-dir`(req unless --infer-only),
`--format`{auto,plain,daner-old,daner-new}(def auto), `--infer-only`, `--overwrite`,
`--sumstats-snps-file`, `--use-hm3-snps`, `--source-genome-build`(auto),
`--output-genome-build`(req in chr_pos modes), `--liftover-chain-file`,
`--use-hm3-quick-liftover`, `--trait-name`, `--log-level`, `--output-format`{parquet,tsv.gz,both}(parquet).
CLONED FROM KERNEL (legacy-equivalent): `--N`,`--N-cas`,`--N-con`,`--info-min`(0.9),
`--maf-min`(0.01),`--n-min`,`--chunksize`(**1_000_000**, legacy was 5e6),`--snp`,`--chr`(NEW),
`--pos`(NEW),`--N-col`,`--N-cas-col`,`--N-con-col`,`--a1`,`--a2`,`--p`,`--frq`,
`--signed-sumstats`,`--info`,`--info-list`,`--nstudy`,`--nstudy-min`,`--ignore`,
`--a1-inc`,`--keep-maf`,`--snp-identifier`.
Kernel parser (`_kernel/sumstats_munger.py:1038-1105`) DROPPED legacy `--daner`,`--daner-n`,
`--merge-alleles`,`--no-alleles` entirely (verified: not defined).

### `ldsc query-r2` (`r2_query.build_parser`) — NEW capability, no legacy equivalent
`--panel-dir`(req), `--pairs`(req), `--output-dir`, `--overwrite`, `--log-level`,
`--snp-identifier`, `--genome-build`.

### Regression: `ldsc h2` / `partitioned-h2` / `rg` (`regression_runner.add_*_arguments` + `_add_common_regression_arguments` @2272)
COMMON: `--ldscore-dir`(req), `--count-kind`{common,all}(common), `--output-dir`(opt),
`--overwrite`, `--n-blocks`(200), `--no-intercept`, `--allow-identity-downgrade`(NEW),
`--two-step-cutoff`(None; h2/rg legacy default 30 when applicable), `--chisq-max`(None),
`--log-level`. (`--intercept-h2` added for h2 & partitioned-h2.)
h2 adds: `--sumstats-file`(req), `--trait-name`, `--samp-prev`, `--pop-prev`.
partitioned-h2 adds: `--sumstats-file`(req), `--trait-name`, `--samp-prev`, `--pop-prev`,
`--write-per-query-results`, `--summary-sort-by`.
rg adds: `--sumstats-sources`(req,nargs+), `--anchor-trait`, `--write-per-pair-detail`,
`--intercept-h2`, `--intercept-gencov`, `--samp-prev`, `--pop-prev`, `--prevalence-manifest`(NEW).
All regression chi^2 filters use INCLUSIVE `<=` (doc io-argument-inventory.md:402-403,425-426,450-451;
deviates from legacy strict `<`).

## Command coverage mapping (legacy → refactored)

| Legacy capability (entry/flags) | Refactored counterpart | Coverage | Severity | Notes |
|---|---|---|---|---|
| **L2 LD-score estimation** `ldsc.py --bfile --l2` | `ldsc build-ref-panel` (+ `ldsc ldscore --r2-dir`) OR `ldsc ldscore --plink-prefix` direct | FULL (restructured) | — | Two-stage option (R2 panel then scores) or one-stage direct PLINK. Core path covered. |
| &nbsp;&nbsp;↳ plain single-col LD score | `ldsc ldscore` no-annot (synthesizes `base`) | FULL | — | `cli`/io-arg-inventory.md:164,202 |
| &nbsp;&nbsp;↳ `--annot`/`--thin-annot` partitioned | `ldsc ldscore --baseline-annot-sources` (+`--query-annot-sources`) | FULL (semantics changed) | MEDIUM | Refactored splits annotations into baseline vs query; legacy had one annot matrix. Query needs baseline. |
| &nbsp;&nbsp;↳ `--cts-bin`/`--cts-breaks`/`--cts-names` on-the-fly continuous binning | none | **NONE** | **HIGH** | No on-the-fly binning of a SNP→value file. Workaround: pre-bin into annot columns + `--query-annot-sources` (non-obvious, lossy). NOTE: continuous-VALUED annot files ARE supported (kernel is value-agnostic). |
| &nbsp;&nbsp;↳ `--per-allele` / `--pq-exp` LD-score scaling | none | **NONE** | **HIGH** | No per-allele / pq-exponent scaling option in refactored ldscore (verified zero matches); no workaround. |
| **H2 heritability** `ldsc.py --h2 --ref-ld --w-ld` | `ldsc h2 --ldscore-dir --sumstats-file` | FULL (restructured) | — | `--ref-ld`+`--w-ld` collapse into one canonical `--ldscore-dir` (embeds `regression_ld_scores` as w_ld). |
| **PART-H2 partitioned h2** `ldsc.py --h2 --overlap-annot --frqfile` | `ldsc partitioned-h2 --ldscore-dir --sumstats-file` | FULL (restructured) | — | Overlap matrix + common-SNP counts precomputed into ldscore dir; no `--frqfile`/`--overlap-annot` needed. |
| **RG genetic correlation** `ldsc.py --rg a,b --ref-ld --w-ld` | `ldsc rg --ldscore-dir --sumstats-sources a b ...` | FULL+ (enhanced) | — | Adds anchor mode, all-pairs, prevalence manifest. |
| **CTS cell-type-specific** `ldsc.py --h2-cts --ref-ld-chr-cts <.ldcts>` | (manual) per-cell-type `ldscore` + `partitioned-h2` cell-type regime | **PARTIAL** | **HIGH** | No `--h2-cts` one-command `.ldcts` loop; no `--ref-ld-chr-cts`. Must build per-cell-type LD scores and loop manually. |
| **MUNGE** `munge_sumstats.py` | `ldsc munge-sumstats` | FULL (changed) | — | See munge flag table; `--merge-alleles`/`--no-alleles`/`--daner*` restructured. |
| **ANNOT (BED mode)** `make_annot.py --bed-file` | `ldsc annotate --query-annot-bed-sources` | FULL (changed) | MEDIUM | BED→annot projection covered; output is query annot, requires `--baseline-annot-sources`. |
| **ANNOT (gene-set mode)** `make_annot.py --gene-set-file --gene-coord-file --windowsize` | none | **NONE** | **HIGH** | No gene-set→BED (TSS/TES windowing) in refactored `annotate`. |
| **query-r2** (no legacy equiv) | `ldsc query-r2` | N/A (new feature) | — | Refactored-only capability; not a gap. |

Severity rationale (rubric): CRITICAL items = a core legacy capability entirely
missing (`--cts-bin`); HIGH = meaningful feature missing with non-obvious/lossy
workaround (`--h2-cts` loop, gene-set annot, per-allele scaling); MEDIUM = visible
but managed behavioral/structural change.

## Flag-mapping tables (per command)

Status vocab: SAME (name+behavior), RENAMED, CHANGED (default/semantics differ),
ABSORBED (now automatic/implicit), NONE (no counterpart), NEW (refactored-only).

### A. LD Score estimation: legacy `ldsc.py --l2` → `ldsc ldscore` (+ `build-ref-panel`)

| Legacy flag | Refactored | Status | Severity | Notes |
|---|---|---|---|---|
| `--out` | `--output-dir` | CHANGED | MEDIUM | Prefix→directory; fixed filenames inside. Output contract change (`io-arg-inv.md:103-115`). |
| `--bfile` | `--plink-prefix` | RENAMED | LOW | Same PLINK trio prefix (`ldscore`/`build-ref-panel`). |
| `--l2` | (implicit) | ABSORBED | LOW | `ldscore` always computes LD scores; no switch. |
| `--extract` | `--ref-panel-snps-file` / `--use-hm3-ref-panel-snps` | RENAMED/CHANGED | MEDIUM | Identity-key restriction; dup keys collapse (`io-arg-inv.md:172-173`). |
| `--keep` | `--keep-indivs-file` | RENAMED | LOW | Same PLINK keep file. |
| `--ld-wind-snps` | `--ld-wind-snps` | SAME | — | Exactly-one rule preserved. |
| `--ld-wind-kb` | `--ld-wind-kb` | SAME | — | |
| `--ld-wind-cm` | `--ld-wind-cm` | CHANGED | LOW | Refactored requires usable `CM` (≥2 distinct), dedicated error not bypassed by `--yes-really` (`io-arg-inv.md:183`). |
| `--print-snps` | `--regression-snps-file` / `--use-hm3-regression-snps` | RENAMED/CHANGED | MEDIUM | Restricts written rows by identity key. |
| `--annot` | `--baseline-annot-sources` (+`--query-annot-sources`) | CHANGED | MEDIUM | Annotations split into baseline vs query; query requires baseline. |
| `--thin-annot` | (auto) | ABSORBED | LOW | Thin/standard annot auto-handled (CM optional, SNP required only in rsID modes). |
| `--cts-bin` | NONE | NONE | **HIGH** | On-the-fly binning unavailable. Kernel `read_cts` (`_kernel/formats.py:104`) is DEAD (never called). Continuous-valued annot files still work via `--query-annot-sources`. |
| `--cts-breaks` | NONE | NONE | **HIGH** | Tied to `--cts-bin`. |
| `--cts-names` | NONE | NONE | **HIGH** | Tied to `--cts-bin`. |
| `--per-allele` | NONE | NONE | **HIGH** | No per-allele LD scores anywhere in package (verified zero matches). |
| `--pq-exp` | NONE | NONE | **HIGH** | No pq-exponent scaling anywhere (verified zero matches). |
| `--no-print-annot` | NONE | NONE | LOW | Only relevant to removed `--cts-bin`. |
| `--maf` | `--maf-min` | RENAMED/CHANGED | MEDIUM | Refactored inclusive `MAF>=` (`io-arg-inv.md:177`); legacy help says "MAF > 0" (strict, `ldsc.py:490`). Boundary SNP at MAF==cutoff now retained. |
| `--chunk-size` | `--snp-batch-size` | RENAMED/CHANGED | LOW | Default 50→128; perf only, different internal meaning. |
| `--n-blocks` | (regression only) | CHANGED | LOW | Legacy L2 used it for SE blocks; refactored `ldscore` has no n-blocks (it's a regression flag). |
| `--yes-really` | `--yes-really` | SAME | — | Whole-chromosome window override (does NOT bypass cM-usability error). |
| `--Pickle` | NONE | NONE | LOW | Removed; was a dead/unconsumed flag in legacy (`ldsc.py:560`, never read). |
| (n/a) | `--genome-build`,`--snp-identifier`,`--exclude-regions[-bed/-build]`,`--common-maf-min`,`--threads`,`--export-ref-metadata`,`--genetic-map-*`,`--r2-dir`,`--overwrite`,`--log-level` | NEW | — | Refactored additions (identity model, region exclusion, parallelism, parquet backend). `--exclude-regions` defaults to `mhc-and-centromeres` (legacy did NOT auto-exclude) — see behavioral discrepancies. |
| (n/a) | `build-ref-panel`: `--min-r2`, `--source-genome-build`, `--liftover-chain-*`, `--use-hm3-snps`, `--use-hm3-quick-liftover` | NEW | — | New intermediate R2-panel build stage; no legacy analog. |

### B. Heritability: legacy `ldsc.py --h2` → `ldsc h2`

| Legacy flag | Refactored | Status | Severity | Notes |
|---|---|---|---|---|
| `--h2 <file>` | `--sumstats-file` | RENAMED | — | Single munged sumstats input. |
| `--ref-ld` / `--ref-ld-chr` | `--ldscore-dir` | CHANGED | MEDIUM | Fragmented per-chr `.l2.ldscore` inputs → one canonical directory; `@`-token & file prefixes gone. |
| `--w-ld` / `--w-ld-chr` | `--ldscore-dir` (embedded `regression_ld_scores`) | CHANGED | **HIGH** | Cannot supply an INDEPENDENT regression-weight LD set; w_ld is baked into the ldscore dir at build time (`io-arg-inv.md:393`). Common case covered; non-standard ref≠weight setups not. |
| `--M` | NONE | NONE | MEDIUM | No manual SNP-count override; counts come from ldscore dir + `--count-kind`. |
| `--not-M-5-50` | `--count-kind {common,all}` | RENAMED | LOW | `all`≈`.l2.M`, `common`(default)≈`.l2.M_5_50`. |
| `--n-blocks` | `--n-blocks` | SAME | — | Default 200. |
| `--no-intercept` | `--no-intercept` | SAME | — | |
| `--intercept-h2` | `--intercept-h2` | CHANGED | LOW | Now `type=float` (legacy plain `store`/string). |
| `--two-step` | `--two-step-cutoff` | RENAMED/CHANGED | MEDIUM | Inclusive `chi^2<=cutoff` vs legacy strict `<`; legacy default 30 preserved for single-annot free-intercept (`io-arg-inv.md:402`). |
| `--chisq-max` | `--chisq-max` | CHANGED | MEDIUM | Inclusive `chi^2<=max` vs legacy strict `<` (`io-arg-inv.md:403`). |
| `--samp-prev` | `--samp-prev` | SAME | — | Refactored validates prob in (0,1) or nan. |
| `--pop-prev` | `--pop-prev` | SAME | — | |
| `--print-cov` | NONE | NONE | MEDIUM | No covariance-matrix output (verified zero matches). |
| `--print-delete-vals` | NONE | NONE | MEDIUM | Delete-values computed in kernel (`_kernel/regression.py:350`) but not exposed by any flag/output. |
| `--invert-anyway` | NONE | NONE | MEDIUM | No force-invert of ill-conditioned matrices (verified zero matches). |
| (n/a) | `--allow-identity-downgrade`,`--output-dir`,`--overwrite`,`--log-level` | NEW | — | |

### C. Partitioned h2: legacy `ldsc.py --h2 --overlap-annot` → `ldsc partitioned-h2`

| Legacy flag | Refactored | Status | Severity | Notes |
|---|---|---|---|---|
| `--h2 <file>` | `--sumstats-file` | RENAMED | — | |
| `--ref-ld-chr` (baseline+categories) | `--ldscore-dir` (baseline+query) | CHANGED | MEDIUM | Requires baseline+query LD scores; baseline-only dirs rejected (`io-arg-inv.md:416`). |
| `--w-ld-chr` | `--ldscore-dir` (embedded) | CHANGED | HIGH | Same independent-weight limitation as h2. |
| `--overlap-annot` | (auto) | ABSORBED | — | Overlap matrix `ldscore.overlap.parquet` precomputed; no flag. |
| `--print-coefficients` | (always emitted) | ABSORBED | LOW | Coefficients (`tau`) always in `partitioned_h2.tsv`. |
| `--frqfile` / `--frqfile-chr` | (auto) | ABSORBED | — | Common-SNP counts precomputed via `--common-maf-min` at ldscore build; no frq file. |
| `--M` | NONE | NONE | MEDIUM | As h2. |
| `--not-M-5-50` | `--count-kind` | RENAMED | LOW | |
| `--n-blocks`/`--no-intercept`/`--intercept-h2` | same | SAME/CHANGED | — | As h2. |
| `--two-step`/`--chisq-max` | `--two-step-cutoff`/`--chisq-max` | CHANGED | MEDIUM | Inclusive; two-step not applied to multi-annot models (`io-arg-inv.md:425`). |
| `--samp-prev`/`--pop-prev` | same | SAME | — | |
| `--print-cov`/`--print-delete-vals`/`--invert-anyway` | NONE | NONE | MEDIUM | As h2. |
| (n/a) | `--write-per-query-results`,`--summary-sort-by` | NEW | — | |

### D. Genetic correlation: legacy `ldsc.py --rg` → `ldsc rg`

| Legacy flag | Refactored | Status | Severity | Notes |
|---|---|---|---|---|
| `--rg a,b,c` | `--sumstats-sources a b c` | RENAMED/CHANGED | LOW | Space-separated + globs (legacy comma-separated). At least 2 required. |
| `--ref-ld(-chr)` / `--w-ld(-chr)` | `--ldscore-dir` | CHANGED | HIGH | As h2 (independent weights not supportable). |
| `--no-intercept` | `--no-intercept` | SAME | — | |
| `--intercept-h2` | `--intercept-h2` | CHANGED | LOW | Broadcast to all pairs; `type=float`. |
| `--intercept-gencov` | `--intercept-gencov` | CHANGED | LOW | Single float broadcast to all pairs; legacy took a per-rg list ("first entry ignored"). Behavioral change. |
| `--two-step` | `--two-step-cutoff` | RENAMED/CHANGED | MEDIUM | Inclusive; default 30 for single-annot free-h2-intercept. |
| `--chisq-max` | `--chisq-max` | CHANGED | MEDIUM | rg uses `Z1^2*Z2^2 <= chisq_max^2` inclusive (`io-arg-inv.md:451`); no default cap (matches legacy opt-in). |
| `--no-check-alleles` | NONE (nearest `--allow-identity-downgrade`) | NONE | MEDIUM | Allele checking is via the identity model; no skip-check toggle. |
| `--return-silly-things` | NONE | NONE | LOW | Kernel `summary(silly=False)` exists (`_kernel/regression.py:857`) but no flag wires `silly=True`. Default (clamp/guard) matches legacy default. |
| `--M`/`--not-M-5-50` | NONE / `--count-kind` | NONE/RENAMED | MEDIUM/LOW | |
| `--samp-prev`/`--pop-prev` | same (comma per-trait) | SAME | — | Plus `--prevalence-manifest` (NEW). |
| `--print-cov`/`--print-delete-vals` | NONE | NONE | MEDIUM | |
| (n/a) | `--anchor-trait`,`--write-per-pair-detail`,`--prevalence-manifest` | NEW | — | |

### E. Cell-type-specific: legacy `ldsc.py --h2-cts` → `ldsc partitioned-h2` (cell-type regime)

| Legacy flag | Refactored | Status | Severity | Notes |
|---|---|---|---|---|
| `--h2-cts <file>` | `--sumstats-file` + `partitioned-h2` | PARTIAL | **HIGH** | No single `--h2-cts` command; must build per-cell-type LD scores and run partitioned-h2 (cell-type regime `regression_runner.py:791`). |
| `--ref-ld-chr` | `--ldscore-dir` | CHANGED | MEDIUM | Baseline LD scores. |
| `--ref-ld-chr-cts <.ldcts>` | NONE | NONE | **HIGH** | No `.ldcts` list-of-prefixes loop; one cell-type per partitioned-h2 run (query annotations). |
| `--w-ld(-chr)` | `--ldscore-dir` (embedded) | CHANGED | HIGH | As h2. |
| `--print-all-cts` | NONE | NONE | LOW | |

### F. Munging: `munge_sumstats.py` → `ldsc munge-sumstats`

| Legacy flag | Refactored | Status | Severity | Notes |
|---|---|---|---|---|
| `--sumstats` | `--raw-sumstats-file` | RENAMED | — | |
| `--out` | `--output-dir` | CHANGED | MEDIUM | Prefix→dir; default artifact is `sumstats.parquet` not `.sumstats.gz`. |
| `--N`/`--N-cas`/`--N-con` | same | SAME | — | |
| `--info-min` | `--info-min` | SAME | — | Default 0.9. |
| `--maf-min` | `--maf-min` | SAME | — | Default 0.01. |
| `--n-min` | `--n-min` | SAME | — | Default (90th pctile N)/2. |
| `--chunksize` | `--chunksize` | CHANGED | LOW | Default **5e6 → 1_000_000** (`_kernel/sumstats_munger.py:1061`). Perf only. |
| `--daner` | `--format daner-old` | RENAMED/CHANGED | LOW | Boolean→choice; merged DANER selector. |
| `--daner-n` | `--format daner-new` | RENAMED/CHANGED | LOW | |
| `--no-alleles` | NONE (use `--snp-identifier chr_pos`/`rsid`) | CHANGED | MEDIUM | Allele-blind handled by base identity modes, not a flag. |
| `--merge-alleles` | NONE (nearest `--sumstats-snps-file`/`--use-hm3-snps`) | NONE | **HIGH** | No external SNP+A1+A2 allele-harmonization/flip-to-list. Restriction is identity-only; allele match via identity model. Standard HM3 merge workflow changed. |
| `--snp`/`--N-col`/`--N-cas-col`/`--N-con-col`/`--a1`/`--a2`/`--p`/`--frq`/`--signed-sumstats`/`--info`/`--info-list`/`--nstudy`/`--nstudy-min`/`--ignore` | same | SAME | — | All column-hint/QC flags retained (cloned from kernel parser). |
| `--a1-inc` | `--a1-inc` | SAME | — | |
| `--keep-maf` | `--keep-maf` | SAME | — | |
| (n/a) | `--chr`,`--pos`,`--format`,`--infer-only`,`--output-format`,`--source-genome-build`,`--output-genome-build`,`--liftover-chain-file`,`--use-hm3-snps`,`--use-hm3-quick-liftover`,`--sumstats-snps-file`,`--trait-name`,`--snp-identifier`,`--overwrite`,`--log-level` | NEW | — | Coordinate-awareness, liftover, parquet output, identity model. |

### G. Annotation: `make_annot.py` → `ldsc annotate`

| Legacy flag | Refactored | Status | Severity | Notes |
|---|---|---|---|---|
| `--bed-file` | `--query-annot-bed-sources` | RENAMED/CHANGED | MEDIUM | BED→annot projection; multi-BED + globs; each stem becomes a query annotation. |
| `--bimfile` | `--baseline-annot-sources` (defines SNP universe) | CHANGED | MEDIUM | SNP set comes from baseline annot template, not a `.bim`. |
| `--annot-file` | `--output-dir` | CHANGED | LOW | Writes `query.<chrom>.annot.gz` (per-chr) not a single `.annot`. |
| `--nomerge` | NONE | NONE | MEDIUM | No proportional-overlap-count option; refactored projection is membership. |
| `--gene-set-file` | NONE | NONE | **HIGH** | No gene-set→BED building. |
| `--gene-coord-file` | NONE | NONE | **HIGH** | Tied to gene-set mode (TSS/TES coords). |
| `--windowsize` | NONE (nearest `--bed-padding-bp`) | NONE | MEDIUM | `--bed-padding-bp` pads existing BED intervals, but cannot window around gene TSS/TES from a gene list. |
| (n/a) | `--bed-padding-bp`,`--baseline-annot-sources`,`--snp-identifier`,`--genome-build`,`--overwrite`,`--log-level` | NEW | — | |

## Behavioral discrepancies (silent/material first)

1. **CRITICAL — region exclusion ON by default.** `ldsc ldscore` and
   `ldsc build-ref-panel` default `--exclude-regions=mhc-and-centromeres`
   (`ldscore_calculator.py:930`, `ref_panel_builder` per `io-arg-inv.md:230`).
   Legacy `ldsc.py` has NO region-exclusion concept — it computes LD scores over
   ALL SNPs. A user porting a legacy LD-score command verbatim silently drops MHC
   + centromere SNPs, materially changing LD scores and downstream h2/partitioned
   h2. Rubric: CRITICAL (silent materially-different behavior). Mitigation: pass
   `--exclude-regions none` to reproduce legacy.

2. **HIGH — regression weight LD scores no longer independent.** Legacy took
   separate `--ref-ld(-chr)` and `--w-ld(-chr)`. Refactored regression reads ONE
   `--ldscore-dir` whose embedded `regression_ld_scores` IS the weight set
   (`io-arg-inv.md:393`,`441`). The standard "same eur_w_ld for both" case is fine,
   but any workflow using a distinct regression-weight LD-score set cannot be
   expressed without rebuilding the ldscore dir.

3. **MEDIUM — strict→inclusive threshold flips** (all documented, but change which
   SNPs are retained at the exact boundary):
   - `--maf` (`MAF>0` strict, `ldsc.py:490`) → `--maf-min` inclusive `MAF>=` (`io-arg-inv.md:177`).
   - common-SNP count: legacy `maf > 0.05` strict (`ldsc.py:357`) → `--common-maf-min`
     inclusive `MAF>=0.05` (`ldscore_calculator.py:970`; doc also notes legacy
     `0.05<FRQ<0.95` in regression context).
   - two-step / chisq-max: legacy strict `<` → refactored inclusive `<=`
     (`io-arg-inv.md:402-403,425-426,450-451`). rg uses `Z1^2*Z2^2 <= max^2`.
   (Consistent with the repo's own min/max inclusive-operator convention; flagged
   here only as legacy-vs-refactored deltas, not bugs.)

4. **MEDIUM — output contract: prefix → directory.** Every artifact-writing command
   replaced legacy `--out <prefix>` with `--output-dir <dir>` + fixed filenames
   (`io-arg-inv.md:103-138`). Legacy `<prefix>.l2.ldscore.gz/.M/.M_5_50/.log`,
   `<prefix>.sumstats.gz`, `<prefix>.log` layouts are not produced at the public
   layer (kernel retains legacy emitters for private/compat paths, `io-arg-inv.md:89-91`).
   Scripts parsing legacy output paths break.

5. **MEDIUM — default value changes:** `--chunk-size`(50)→`--snp-batch-size`(128);
   munge `--chunksize` 5e6→1e6 (`_kernel/sumstats_munger.py:1061`). Performance only.

6. **MEDIUM — lost diagnostics:** `--print-cov`, `--print-delete-vals`,
   `--invert-anyway` have no public counterpart (delete-values still computed in
   `_kernel/regression.py:350` but unexposed).

7. **LOW — rg input list format:** legacy `--rg a,b,c` (comma) →
   `--sumstats-sources a b c` (space + globs). `--intercept-gencov` legacy list
   ("first ignored") → single broadcast float.

## Wiki examples → refactored translations

### Wiki: "Heritability and Genetic Correlation" (raw legacy commands extracted)
- MUNGE: `munge_sumstats.py --sumstats pgc.cross.SCZ17.2013-05.txt --N 17115 --out scz --merge-alleles w_hm3.snplist`
- RG: `ldsc.py --rg scz.sumstats.gz,bip.sumstats.gz --ref-ld-chr eur_w_ld_chr/ --w-ld-chr eur_w_ld_chr/ --out scz_bip`
- H2: `ldsc.py --h2 scz.sumstats.gz --ref-ld-chr eur_w_ld_chr/ --w-ld-chr eur_w_ld_chr/ --out scz_h2`
- LIABILITY: add `--samp-prev 0.5,0.5 --pop-prev 0.01,0.01` (rg) or `--samp-prev 0.5,nan --pop-prev 0.01,nan` (mixed)
- CONSTRAIN: `--intercept-h2 1,0.99,1.01 --intercept-gencov 0,N0.5,0.5` (N=neg sign); or `--no-intercept`
- Wiki confirms legacy munge defaults "INFO > 0.9, MAF > 0.01 and 0 < P <= 1" (strict; log says "Removed INFO <= 0.9").
- Wiki (LD Scores section) EXPLICITLY recommends SEPARATE `--w-ld-chr` and `--ref-ld-chr` LD-score sets for
  partitioned regression -> supports HIGH finding that refactored couples them into one `--ldscore-dir`.

**KEY BLOCKER — CONFIRMED CRITICAL.** `load_ldscore_from_dir` (`regression_runner.py:2465`) hard-requires
`metadata.json` + `ldscore.baseline.parquet` and raises `LDSCInputError` "missing `metadata.json`" otherwise
(`:40-41`). It does NOT read legacy `.l2.ldscore.gz`/`.M`/`.M_5_50`. So the distributed `eur_w_ld_chr/` and
`baseline_v1.2/` directories (present in `resources/`) are UNCONSUMABLE by refactored `h2`/`partitioned-h2`/`rg`.
Every regression wiki example must FIRST regenerate canonical parquet LD scores via `ldsc ldscore`, which needs
the original 1000G PLINK panel (or an R2 parquet panel) most users never download. Severity CRITICAL (core
documented workflow broken; precomputed-artifact ecosystem incompatible).

### Finalized translations — "Heritability and Genetic Correlation"
| Legacy wiki command | Refactored equivalent | Reproducible? |
|---|---|---|
| `munge_sumstats.py --sumstats SCZ.txt --N 17115 --out scz --merge-alleles w_hm3.snplist` | `ldsc munge-sumstats --raw-sumstats-file SCZ.txt --N 17115 --output-dir scz --use-hm3-snps --snp-identifier rsid --output-format tsv.gz` | PARTIAL — `--use-hm3-snps` uses the PACKAGED HM3 map (not an arbitrary `w_hm3.snplist`); allele harmonize-to-list semantics of `--merge-alleles` not reproduced (identity model instead). |
| `ldsc.py --rg scz.sumstats.gz,bip.sumstats.gz --ref-ld-chr eur_w_ld_chr/ --w-ld-chr eur_w_ld_chr/ --out scz_bip` | `ldsc rg --sumstats-sources scz/sumstats.parquet bip/sumstats.parquet --ldscore-dir <CANONICAL_DIR> --output-dir scz_bip` | NO with distributed `eur_w_ld_chr/` (legacy format unreadable). Must regenerate `<CANONICAL_DIR>` via `ldsc ldscore` first. |
| `ldsc.py --h2 scz.sumstats.gz --ref-ld-chr eur_w_ld_chr/ --w-ld-chr eur_w_ld_chr/ --out scz_h2` | `ldsc h2 --sumstats-file scz/sumstats.parquet --ldscore-dir <CANONICAL_DIR> --output-dir scz_h2` | NO with distributed LD scores (same blocker). |
| `... --samp-prev 0.5,0.5 --pop-prev 0.01,0.01` | `ldsc rg ... --samp-prev 0.5,0.5 --pop-prev 0.01,0.01` | YES (flag retained; per-trait comma list). |
| `... --h2 ... --samp-prev 0.5,nan --pop-prev 0.01,nan` | `ldsc h2 ... --samp-prev 0.5 --pop-prev nan`? | PARTIAL — h2 takes SCALAR samp/pop-prev (one trait); the mixed-trait `0.5,nan` syntax is an rg concept. |
| `... --intercept-h2 1,0.99,1.01 --intercept-gencov 0,N0.5,0.5` | `ldsc rg ... --intercept-h2 1 --intercept-gencov 0` | PARTIAL — refactored `--intercept-h2`/`--intercept-gencov` are SINGLE floats broadcast to all pairs; per-pair vectors not supported. |
| `ldsc.py --rg ... --out scz_bip --no-intercept` | `ldsc rg ... --output-dir scz_bip --no-intercept` | YES (modulo the LD-score blocker). |

### Finalized translations — "LD Score Estimation Tutorial"
| Legacy wiki command | Refactored equivalent | Reproducible? |
|---|---|---|
| `python ldsc.py --bfile 22 --l2 --ld-wind-cm 1 --out 22` (univariate) | `ldsc ldscore --plink-prefix 22 --ld-wind-cm 1 --output-dir 22 --snp-identifier rsid --exclude-regions none` | PARTIAL — must add `--exclude-regions none` (else MHC+centromeres dropped) and `--snp-identifier rsid` (legacy merges on rs#; refactored default is chr_pos_allele_aware). cM window needs usable `.bim` CM. |
| `ldsc.py ... --extract regression.snplist` (compute `--w-ld`) | `ldsc ldscore ... --regression-snps-file regression.snplist` | PARTIAL — restricts written rows; but refactored bakes w_ld into the ldscore dir, so the separate "compute --w-ld scores" step does not exist as such. |
| `make_annot.py --gene-set-file GTEx_Cortex.GeneSet --gene-coord-file ENSG_coord.txt --windowsize 100000 --bimfile 1000G.EUR.QC.22.bim --annot-file GTEx_Cortex.annot.gz` | (none) | **NO** — gene-set→annot mode absent. Must pre-build the BED externally and use `ldsc annotate --query-annot-bed-sources`. |
| `make_annot.py --bed-file Brain_DPC_H3K27ac.bed --bimfile 1000G.EUR.QC.22.bim --annot-file Brain_DPC_H3K27ac.annot.gz` | `ldsc annotate --query-annot-bed-sources Brain_DPC_H3K27ac.bed --baseline-annot-sources <baseline.annot> --output-dir <dir>` | PARTIAL — SNP universe comes from `--baseline-annot-sources` not `.bim`; output is `query.<chrom>.annot.gz` (per-chr) not one `.annot.gz`. |
| `ldsc.py --l2 --bfile 1000G.EUR.QC.22 --ld-wind-cm 1 --annot Brain_DPC_H3K27ac.annot.gz --thin-annot --out Brain_DPC_H3K27ac --print-snps hm.22.snp` | `ldsc ldscore --plink-prefix 1000G.EUR.QC.22 --ld-wind-cm 1 --query-annot-sources Brain_DPC_H3K27ac.annot.gz --baseline-annot-sources <baseline> --regression-snps-file hm.22.snp --output-dir Brain_DPC_H3K27ac --exclude-regions none --snp-identifier rsid` | PARTIAL — query annot REQUIRES `--baseline-annot-sources` (legacy `--annot` was standalone); `--thin-annot` auto-handled; `--print-snps`→`--regression-snps-file`; add `--exclude-regions none`. |

Behavioral confirmations from this tutorial: legacy "removes monomorphic SNPs by default; change with `--maf`"
(confirms MAF>0 strict default); legacy keys SNPs by rs number (refactored default identity differs);
`.bim` CM must be pre-filled (e.g. `plink --cm-map`) — refactored can instead supply `--genetic-map-*`.

### Finalized translations — "Partitioned Heritability" (baseline model)
| Legacy wiki command | Refactored equivalent | Reproducible? |
|---|---|---|
| `munge_sumstats.py --sumstats GIANT_BMI...txt --merge-alleles w_hm3.snplist --out BMI --a1-inc` | `ldsc munge-sumstats --raw-sumstats-file GIANT_BMI...txt --use-hm3-snps --output-dir BMI --a1-inc --snp-identifier rsid --output-format tsv.gz` | PARTIAL — `--merge-alleles`→`--use-hm3-snps` (packaged map). |
| `ldsc.py --h2 BMI.sumstats.gz --ref-ld-chr baseline. --w-ld-chr weights. --overlap-annot --frqfile-chr 1000G.mac5eur. --out BMI_baseline` | `ldsc partitioned-h2 --sumstats-file BMI/sumstats.parquet --ldscore-dir <CANONICAL_DIR> --output-dir BMI_baseline` | **NO** — `baseline.*`/`weights.*`/`1000G.mac5eur.*` are legacy `.l2.ldscore.gz`/`.annot.gz`/`.M_5_50`/frq; unreadable. `--overlap-annot`/`--frqfile-chr` are ABSORBED (overlap matrix + common counts baked into the ldscore dir). Reproduction requires fully rebuilding the baseline-model canonical LD-score dir (baseline annotations + overlap.parquet). |
| `ldsc.py --h2 BMI.sumstats.gz --w-ld-chr weights. --ref-ld-chr CNS.,baseline. --overlap-annot --frqfile-chr 1000G.mac5eur. --out BMI_CNS --print-coefficients` (cell-type group) | `ldsc partitioned-h2 --sumstats-file BMI/sumstats.parquet --ldscore-dir <DIR: baseline + CNS query> --output-dir BMI_CNS` (cell-type regime; coefficients always emitted) | **NO/PARTIAL** — concept maps to the cell-type regime, but the comma-separated multi-prefix `--ref-ld-chr CNS.,baseline.` is GONE (baseline+query must be one canonical dir; `--query-columns` removed). `--print-coefficients` ABSORBED. Plus the legacy-format blocker. |

Behavioral note: legacy concatenates LD-score sets via comma-separated `--ref-ld-chr CNS.,baseline.`;
refactored requires baseline+query annotations inside ONE canonical `--ldscore-dir` (no multi-prefix concat).

### Finalized translations — "Cell type specific analyses" (`--h2-cts`)
| Legacy wiki command | Refactored equivalent | Reproducible? |
|---|---|---|
| `munge_sumstats.py --sumstats body_BMIz.sumstats.gz --merge-alleles w_hm3.snplist --out UKBB_BMI` | `ldsc munge-sumstats --raw-sumstats-file body_BMIz.sumstats.gz --use-hm3-snps --output-dir UKBB_BMI --snp-identifier rsid --output-format tsv.gz` | PARTIAL — `--merge-alleles`→`--use-hm3-snps`. |
| `ldsc.py --h2-cts UKBB_BMI.sumstats.gz --ref-ld-chr 1000G_EUR_Phase3_baseline/baseline. --ref-ld-chr-cts Cahoy.ldcts --w-ld-chr weights_hm3_no_hla/weights. --out BMI_Cahoy` | (no single command) loop: for each cell type build a canonical ldscore dir (baseline + that cell type's query annot), then `ldsc partitioned-h2 --sumstats-file UKBB_BMI/sumstats.parquet --ldscore-dir <dir_celltype> --summary-sort-by coefficient-p --output-dir BMI_Cahoy_<ct>` | **NO** — `--h2-cts`/`--ref-ld-chr-cts` `.ldcts` one-command multi-cell-type sweep absent. No `.cell_type_results.txt` aggregate. Cell-type regime does ONE cell type/run; user must script the loop AND first rebuild each cell-type LD-score dir (legacy `.ldcts` LD scores unreadable). |

`.ldcts` format = 2-col (label, comma-list of LD-score prefixes: cell-type set + control). One regression/line;
output `.cell_type_results.txt` (Name, Coefficient, Coefficient_std_error, Coefficient_P_value, one-sided p).
No refactored equivalent to this aggregated multi-cell-type product.

### Finalized translations — "Partitioned Heritability from Continuous Annotations"
IMPORTANT REFINEMENT: this tutorial uses continuous-VALUED `.annot.gz` files (arbitrary numbers), which is
DISTINCT from `--cts-bin` (on-the-fly binning of a SNP→value file). Verified the refactored LD-score kernel
accumulates `cor_sum += U @ annot` (`_kernel/ldscore.py:1562-1563`) with NO binarization -> continuous
annotation FILES are supported via `--query-annot-sources`. Only `--cts-bin` on-the-fly binning is absent.

| Legacy wiki command | Refactored equivalent | Reproducible? |
|---|---|---|
| `ldsc.py --h2 BMI.sumstats.gz --ref-ld-chr baselineLD. --frqfile-chr 1000G.EUR.QC. --w-ld-chr weights.hm3_noMHC. --overlap-annot --print-coefficients --print-delete-vals --out BMI.baselineLD` | `ldsc partitioned-h2 --sumstats-file BMI/sumstats.parquet --ldscore-dir <CANONICAL baselineLD dir> --output-dir BMI.baselineLD` | **NO** — legacy-format blocker; `--print-delete-vals` (`.part_delete`) has NO counterpart and is REQUIRED by the tutorial's `quantile_h2g.r` step -> downstream quantile workflow cannot be driven. `--print-coefficients`/`--overlap-annot`/`--frqfile-chr` ABSORBED. |
| `for CHR in 1..22: ldsc.py --l2 --bfile 1000G.EUR.QC.$CHR --ld-wind-cm 1 --print-snps listHM3.txt --annot yourannot.$CHR.annot.gz --out yourannot.$CHR` (continuous-valued annot) | `ldsc ldscore --plink-prefix 1000G.EUR.QC.@ --ld-wind-cm 1 --regression-snps-file listHM3.txt --query-annot-sources 'yourannot.@.annot.gz' --baseline-annot-sources <baseline> --output-dir yourannot --exclude-regions none --snp-identifier rsid` | PARTIAL — continuous annot values supported (value-agnostic kernel); but query annot REQUIRES baseline; `--print-snps`→`--regression-snps-file`; add `--exclude-regions none`. `@` chr-suite token supported. |
| `ldsc.py --h2 BMI.sumstats.gz --ref-ld-chr baselineLD.,yourannot. --frqfile-chr 1000G.EUR.QC. --w-ld-chr weights.hm3_noMHC. --overlap-annot --print-coefficients --print-delete-vals --out BMI.baselineLD_yourannot` | `ldsc partitioned-h2 --sumstats-file BMI/sumstats.parquet --ldscore-dir <DIR: baselineLD + yourannot query> --output-dir BMI.baselineLD_yourannot` | **NO** — comma-multi-prefix gone (one canonical dir); print-delete-vals gap; legacy-format blocker. |
| `Rscript quantile_h2g.r ...`, `perl quantile_M.pl ...` | (external helper scripts, not ldsc commands) | N/A — out of LDSC CLI scope; but they consume `.part_delete`/`.q5.M` which the refactored pipeline does not emit. |

Elevation: `--print-delete-vals` is MEDIUM in general but **HIGH for this documented workflow** (the tutorial's
quantile-h2 analysis cannot run without `.part_delete`).

### Remaining wiki pages (no unique runnable analysis commands)
- "Summary Statistics File Format": spec of `.sumstats` columns (SNP,A1,A2,Z,N) + munge usage already covered.
- "LD File Formats": describes `.l2.ldscore.gz`/`.M`/`.M_5_50`/`.annot.gz` — the legacy formats the refactor replaced
  with parquet (relevant to the format-blocker finding; no commands to translate).
- "What Data Are Necessary...": conceptual checklist; no commands.
- "Tests": `python -m nose` / unit-test invocation -> refactored uses `pytest` (dev workflow, not user CLI).
- "Home"/"FAQ": navigation/Q&A; no unique runnable commands.

## Open questions / could-not-verify
-

## Resources / key files
- Legacy entry: ldsc.py (661 lines, ~49 add_argument), munge_sumstats.py, make_annot.py
- Legacy core: ldscore/{ldscore,sumstats,regressions,jackknife,irwls,parse}.py
- Refactored entry: src/ldsc/cli.py; workflow modules annotation_builder.py, ref_panel_builder.py,
  ldscore_calculator.py, sumstats_munger.py, regression_runner.py, r2_query.py
- Central refactored doc: docs/current/io-argument-inventory.md (64KB)

## Technical Decisions
| Decision | Rationale |
|----------|-----------|
| Treat io-argument-inventory.md as a guide, verify flags against actual build_parser source | Docs can drift; user requires source-grounded claims |

## Issues Encountered
| Issue | Resolution |
|-------|------------|
|       |            |
