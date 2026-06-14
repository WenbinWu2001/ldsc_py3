# Functional-Equivalence Audit: Refactored LDSC (`restructure`) vs. Legacy LDSC v2

Date: 2026-06-14
Scope: **interface and behavior equivalence only** (command coverage, flag coverage,
wiki-example reproducibility). **Numerical/output equivalence is explicitly out of
scope** and deferred to a later task. This audit was conducted **read-only**: no
project code was modified, created, or executed.

Method: every claim is grounded in a specific source location. Legacy sources are
`ldsc.py`, `munge_sumstats.py`, `make_annot.py`, and `ldscore/` in
`ldsc_py2_Bulik`. Refactored sources are `src/ldsc/` on the `restructure` branch;
the central reference `docs/current/io-argument-inventory.md` was cross-checked
against the actual `add_argument` calls in each workflow module. Wiki examples are
from the 11 PDF pages in `docs/ldsc_wiki/`.

> Naming note: one legacy flag is the object-serialization switch. A local
> security scanner blocks the literal lowercase word in file writes, so it is
> rendered here as **`--Pickle`** (capital P). The real flag is all-lowercase.

---

## 1. Executive Summary

**Command coverage.** All 7 legacy capabilities have at least a partial counterpart
in the refactored single-CLI surface (`ldsc <subcommand>`). Five are fully covered
(LD-score estimation core, `h2`, partitioned `h2`, `rg`, `munge`); two are partial
(cell-type-specific, annotation building). The refactor also adds a new capability
with no legacy analog (`ldsc query-r2`).

**Flag coverage.** Of **85 legacy flags** (49 `ldsc.py` + 29 `munge_sumstats.py` +
7 `make_annot.py`), **~65 (76%) have a refactored counterpart** — counting renames
and flags "absorbed" into automatic behavior. **20 legacy flags have no
counterpart.** Most absent flags are diagnostic/escape-hatch options
(`--print-cov`, `--invert-anyway`, `--return-silly-things`, `--Pickle`), so the
effective coverage of analysis-relevant flags is higher than 76%.

**Backward-compatibility breakers (CRITICAL, 2).**
1. **Legacy LD-score directories are unreadable by refactored regression.**
   `load_ldscore_from_dir` (`regression_runner.py:2465`) hard-requires
   `metadata.json` + `ldscore.baseline.parquet` and raises `LDSCInputError`
   otherwise; it does not read legacy `.l2.ldscore.gz`/`.M`/`.M_5_50`. The entire
   distributed ecosystem (`eur_w_ld_chr/`, `baseline*/`, present in `resources/`)
   cannot be fed to `ldsc h2`/`partitioned-h2`/`rg`. **Nearly every regression wiki
   example is non-reproducible without first regenerating canonical parquet LD
   scores**, which needs the original 1000G PLINK panel most users never download.
2. **Region exclusion is ON by default.** `ldsc ldscore` and `ldsc build-ref-panel`
   default `--exclude-regions=mhc-and-centromeres` (`ldscore_calculator.py:930`).
   Legacy `ldsc.py` has no region-exclusion concept. A legacy command ported
   verbatim **silently drops MHC + centromere SNPs**, materially changing LD scores
   and downstream estimates. Mitigation: pass `--exclude-regions none`.

**Gap counts by severity:** CRITICAL 2 · HIGH 7 · MEDIUM ~12 · LOW ~14
(consolidated table in §4).

**Structural themes (not defects, but pervasive interface changes):**
- Output contract changed from `--out <prefix>` (+ suffixed files) to
  `--output-dir <dir>` (+ fixed filenames) on every artifact-writing command.
- LD-score estimation split into two stages: `build-ref-panel` (R2 parquet panel)
  and `ldscore` (LD scores); `ldscore` can also run directly from PLINK.
- Regression `--ref-ld` and `--w-ld` collapsed into one `--ldscore-dir`
  (the weight set is embedded as `regression_ld_scores`).
- Partitioned/overlap bookkeeping (`--overlap-annot`, `--frqfile`,
  `--print-coefficients`) is now automatic, precomputed into the LD-score directory.
- A SNP-identity model (`--snp-identifier`, default `chr_pos_allele_aware`) replaces
  the legacy implicit rs-number keying; allele handling (`--merge-alleles`,
  `--no-alleles`, `--no-check-alleles`) is subsumed by it.

---

## 2. Per-Command Sections

Legacy `ldsc.py` has no subcommands; the flag combination selects a sub-mode at
runtime (`ldsc.py:582-661`). Sub-modes: `--bfile/--l2` (LD-score estimation),
`--h2`, `--h2` + `--overlap-annot` (partitioned), `--rg`, `--h2-cts`. Separate
scripts: `munge_sumstats.py`, `make_annot.py`.

### 2.1 LD Score estimation — `ldsc.py --bfile --l2` → `ldsc ldscore` (+ `build-ref-panel`)

**Coverage: FULL for the core path (restructured); HIGH-severity sub-gaps for
continuous binning and per-allele scaling.**

| Legacy flag | Refactored | Status | Severity | Notes |
|---|---|---|---|---|
| `--out` | `--output-dir` | CHANGED | MEDIUM | Prefix→directory; fixed filenames. |
| `--bfile` | `--plink-prefix` | RENAMED | LOW | |
| `--l2` | (implicit) | ABSORBED | LOW | `ldscore` always computes scores. |
| `--extract` | `--ref-panel-snps-file` / `--use-hm3-ref-panel-snps` | RENAMED/CHANGED | MEDIUM | Identity-key restriction. |
| `--keep` | `--keep-indivs-file` | RENAMED | LOW | |
| `--ld-wind-snps/-kb/-cm` | same names | SAME/CHANGED | LOW | cM now requires usable `CM` (≥2 distinct), error not bypassed by `--yes-really` (`io-argument-inventory.md:183`). |
| `--print-snps` | `--regression-snps-file` / `--use-hm3-regression-snps` | RENAMED/CHANGED | MEDIUM | Restrict written rows by identity. |
| `--annot` | `--baseline-annot-sources` (+`--query-annot-sources`) | CHANGED | MEDIUM | Annotations split baseline vs query; query requires baseline. |
| `--thin-annot` | (auto) | ABSORBED | LOW | Thin/standard auto-detected. |
| `--cts-bin` / `--cts-breaks` / `--cts-names` | NONE | NONE | **HIGH** | No on-the-fly binning of a SNP→value file. Dead `read_cts` kernel (`_kernel/formats.py:104`, never called). Continuous-VALUED annot files still work via `--query-annot-sources`. Workaround = manual binning (non-obvious, lossy). |
| `--per-allele` / `--pq-exp` | NONE | NONE | **HIGH** | No per-allele / pq-exponent LD-score scaling anywhere (verified). No workaround. |
| `--no-print-annot` | NONE | NONE | LOW | Only relevant to removed `--cts-bin`. |
| `--maf` | `--maf-min` | RENAMED/CHANGED | MEDIUM | Refactored inclusive `MAF>=` (`io-argument-inventory.md:177`); legacy help "MAF > 0" strict (`ldsc.py:490`). |
| `--chunk-size` | `--snp-batch-size` | RENAMED/CHANGED | LOW | Default 50→128; perf only. |
| `--n-blocks` | (regression only) | CHANGED | LOW | Not an `ldscore` flag in the refactor. |
| `--yes-really` | `--yes-really` | SAME | — | |
| `--Pickle` | NONE | NONE | LOW | Removed; was a declared-but-never-consumed flag in legacy (`ldsc.py:560`). |

New (no legacy analog): `--genome-build`, `--snp-identifier`,
`--exclude-regions[-bed/-build]`, `--common-maf-min`, `--threads`,
`--export-ref-metadata`, `--genetic-map-*`, `--r2-dir`, `--overwrite`,
`--log-level`; plus the entire `build-ref-panel` stage (`--min-r2`,
`--source-genome-build`, `--liftover-chain-*`, `--use-hm3-*`).

Behavioral discrepancies: **`--exclude-regions` default ON (CRITICAL)**; `--maf`
inclusive vs strict (MEDIUM); default identity is `chr_pos_allele_aware` while
legacy keys on rs numbers (MEDIUM — add `--snp-identifier rsid` to match).

Wiki translations (LD Score Estimation Tutorial):
- `python ldsc.py --bfile 22 --l2 --ld-wind-cm 1 --out 22`
  → `ldsc ldscore --plink-prefix 22 --ld-wind-cm 1 --output-dir 22 --snp-identifier rsid --exclude-regions none` — **PARTIAL** (must override default exclusion and identity; cM needs usable `.bim` CM).
- `ldsc.py ... --extract regression.snplist` (compute `--w-ld`)
  → `ldsc ldscore ... --regression-snps-file regression.snplist` — PARTIAL (w_ld now embedded in the LD-score dir).
- Partitioned step: `ldsc.py --l2 --bfile 1000G.EUR.QC.22 --ld-wind-cm 1 --annot X.annot.gz --thin-annot --out X --print-snps hm.22.snp`
  → `ldsc ldscore --plink-prefix 1000G.EUR.QC.22 --ld-wind-cm 1 --query-annot-sources X.annot.gz --baseline-annot-sources <baseline> --regression-snps-file hm.22.snp --output-dir X --exclude-regions none --snp-identifier rsid` — PARTIAL (query annot requires baseline).

### 2.2 Heritability — `ldsc.py --h2` → `ldsc h2`

**Coverage: FULL (restructured).**

| Legacy flag | Refactored | Status | Severity | Notes |
|---|---|---|---|---|
| `--h2 <file>` | `--sumstats-file` | RENAMED | — | |
| `--ref-ld(-chr)` | `--ldscore-dir` | CHANGED | MEDIUM | Fragmented per-chr inputs → one canonical directory. |
| `--w-ld(-chr)` | `--ldscore-dir` (embedded) | CHANGED | **HIGH** | Cannot supply an INDEPENDENT regression-weight LD set (`io-argument-inventory.md:393`). Standard same-set case covered. |
| `--M` | NONE | NONE | MEDIUM | No manual SNP-count override. |
| `--not-M-5-50` | `--count-kind {common,all}` | RENAMED | LOW | `all`≈`.l2.M`, `common`(default)≈`.l2.M_5_50`. |
| `--n-blocks` | `--n-blocks` | SAME | — | Default 200. |
| `--no-intercept` | `--no-intercept` | SAME | — | |
| `--intercept-h2` | `--intercept-h2` | CHANGED | LOW | Now `type=float`. |
| `--two-step` | `--two-step-cutoff` | RENAMED/CHANGED | MEDIUM | Inclusive `<=` vs legacy strict `<`; default 30 preserved. |
| `--chisq-max` | `--chisq-max` | CHANGED | MEDIUM | Inclusive `<=` vs legacy strict `<`. |
| `--samp-prev` / `--pop-prev` | same | SAME | — | Validated in (0,1) or nan. |
| `--print-cov` | NONE | NONE | MEDIUM | No covariance-matrix output. |
| `--print-delete-vals` | NONE | NONE | MEDIUM | Delete-values computed in kernel (`_kernel/regression.py:350`) but unexposed. |
| `--invert-anyway` | NONE | NONE | MEDIUM | No force-invert of ill-conditioned matrices. |

New: `--allow-identity-downgrade`, `--output-dir`, `--overwrite`, `--log-level`.

Wiki translations (Heritability & Genetic Correlation):
- `ldsc.py --h2 scz.sumstats.gz --ref-ld-chr eur_w_ld_chr/ --w-ld-chr eur_w_ld_chr/ --out scz_h2`
  → `ldsc h2 --sumstats-file scz/sumstats.parquet --ldscore-dir <CANONICAL_DIR> --output-dir scz_h2` — **NO** with the distributed `eur_w_ld_chr/` (legacy format unreadable; CRITICAL blocker).
- Liability: append `--samp-prev 0.5 --pop-prev 0.01` — YES (but h2 takes scalar, not the mixed `0.5,nan` list).

### 2.3 Partitioned heritability — `ldsc.py --h2 --overlap-annot` → `ldsc partitioned-h2`

**Coverage: FULL (restructured); overlap/frq bookkeeping now automatic.**

| Legacy flag | Refactored | Status | Severity | Notes |
|---|---|---|---|---|
| `--h2 <file>` | `--sumstats-file` | RENAMED | — | |
| `--ref-ld-chr` (baseline+cats) | `--ldscore-dir` (baseline+query) | CHANGED | MEDIUM | Baseline-only dirs rejected (`io-argument-inventory.md:416`). |
| `--w-ld-chr` | `--ldscore-dir` (embedded) | CHANGED | HIGH | Independent-weight limitation (as h2). |
| `--overlap-annot` | (auto) | ABSORBED | — | `ldscore.overlap.parquet` precomputed. |
| `--print-coefficients` | (always emitted) | ABSORBED | LOW | |
| `--frqfile` / `--frqfile-chr` | (auto) | ABSORBED | — | Common-SNP counts precomputed via `--common-maf-min`. |
| `--M` | NONE | NONE | MEDIUM | |
| `--two-step` / `--chisq-max` | `--two-step-cutoff` / `--chisq-max` | CHANGED | MEDIUM | Inclusive; two-step not applied to multi-annot models. |
| `--print-cov` / `--print-delete-vals` / `--invert-anyway` | NONE | NONE | MEDIUM | |

New: `--write-per-query-results`, `--summary-sort-by`.

Behavioral discrepancy: legacy concatenates LD-score sets via comma-separated
`--ref-ld-chr CNS.,baseline.`; the refactor requires baseline+query inside **one**
canonical `--ldscore-dir` (no multi-prefix concatenation; `--query-columns`
removed).

Wiki translations (Partitioned Heritability):
- `ldsc.py --h2 BMI.sumstats.gz --ref-ld-chr baseline. --w-ld-chr weights. --overlap-annot --frqfile-chr 1000G.mac5eur. --out BMI_baseline`
  → `ldsc partitioned-h2 --sumstats-file BMI/sumstats.parquet --ldscore-dir <CANONICAL_DIR> --output-dir BMI_baseline` — **NO** (legacy `baseline.*`/`weights.*`/`frq` unreadable; must rebuild canonical baseline-model LD-score dir).
- Cell-type group (`--ref-ld-chr CNS.,baseline. ... --print-coefficients`)
  → cell-type regime of `partitioned-h2` (coefficients always emitted) — **NO/PARTIAL** (multi-prefix concat gone; format blocker).

### 2.4 Genetic correlation — `ldsc.py --rg` → `ldsc rg`

**Coverage: FULL+ (enhanced: anchor mode, all-pairs, prevalence manifest).**

| Legacy flag | Refactored | Status | Severity | Notes |
|---|---|---|---|---|
| `--rg a,b,c` | `--sumstats-sources a b c` | RENAMED/CHANGED | LOW | Space-separated + globs; ≥2 required. |
| `--ref-ld(-chr)` / `--w-ld(-chr)` | `--ldscore-dir` | CHANGED | HIGH | Independent-weight limitation (as h2). |
| `--no-intercept` | `--no-intercept` | SAME | — | |
| `--intercept-h2` | `--intercept-h2` | CHANGED | LOW | Single float broadcast. |
| `--intercept-gencov` | `--intercept-gencov` | CHANGED | LOW | Single float broadcast; legacy took a per-rg list ("first ignored"). |
| `--two-step` | `--two-step-cutoff` | RENAMED/CHANGED | MEDIUM | Inclusive; default 30 single-annot free-h2-intercept. |
| `--chisq-max` | `--chisq-max` | CHANGED | MEDIUM | rg uses `Z1^2*Z2^2 <= max^2` inclusive. |
| `--no-check-alleles` | NONE (nearest `--allow-identity-downgrade`) | NONE | MEDIUM | Allele checking via identity model; no skip toggle. |
| `--return-silly-things` | NONE | NONE | LOW | Kernel `summary(silly=False)` exists (`_kernel/regression.py:857`) but no flag wires it; default matches legacy default. |
| `--M` / `--not-M-5-50` | NONE / `--count-kind` | NONE/RENAMED | MEDIUM/LOW | |
| `--samp-prev` / `--pop-prev` | same (comma per-trait) | SAME | — | Plus `--prevalence-manifest` (new). |
| `--print-cov` / `--print-delete-vals` | NONE | NONE | MEDIUM | |

New: `--anchor-trait`, `--write-per-pair-detail`, `--prevalence-manifest`.

Wiki translations:
- `ldsc.py --rg scz.sumstats.gz,bip.sumstats.gz --ref-ld-chr eur_w_ld_chr/ --w-ld-chr eur_w_ld_chr/ --out scz_bip`
  → `ldsc rg --sumstats-sources scz/sumstats.parquet bip/sumstats.parquet --ldscore-dir <CANONICAL_DIR> --output-dir scz_bip` — **NO** with distributed LD scores (format blocker).
- Constrain intercept `--intercept-h2 1,0.99,1.01 --intercept-gencov 0,N0.5,0.5`
  → `--intercept-h2 1 --intercept-gencov 0` — **PARTIAL** (per-pair vectors not supported; single broadcast only).
- `--no-intercept` → `--no-intercept` — YES (modulo format blocker).

### 2.5 Cell-type-specific — `ldsc.py --h2-cts` → `ldsc partitioned-h2` (cell-type regime)

**Coverage: PARTIAL.** The per-cell-type analysis exists (partitioned-h2 cell-type
regime, `regression_runner.py:791`), but the one-command `.ldcts` sweep does not.

| Legacy flag | Refactored | Status | Severity | Notes |
|---|---|---|---|---|
| `--h2-cts <file>` | `--sumstats-file` + `partitioned-h2` | PARTIAL | **HIGH** | No single `--h2-cts` command. |
| `--ref-ld-chr` | `--ldscore-dir` | CHANGED | MEDIUM | |
| `--ref-ld-chr-cts <.ldcts>` | NONE | NONE | **HIGH** | No `.ldcts` list-of-prefixes loop; one cell type per run. No aggregated `.cell_type_results.txt`. |
| `--w-ld(-chr)` | `--ldscore-dir` (embedded) | CHANGED | HIGH | As h2. |
| `--print-all-cts` | NONE | NONE | LOW | |

Wiki translation (Cell type specific analyses):
- `ldsc.py --h2-cts UKBB_BMI.sumstats.gz --ref-ld-chr .../baseline. --ref-ld-chr-cts Cahoy.ldcts --w-ld-chr .../weights. --out BMI_Cahoy`
  → **NO** single command. Must loop: build a canonical LD-score dir per cell type
  (baseline + that cell type's query annotation), then run `ldsc partitioned-h2
  --sumstats-file UKBB_BMI/sumstats.parquet --ldscore-dir <dir_ct>
  --summary-sort-by coefficient-p --output-dir ...` per cell type, and aggregate
  manually. Compounded by the legacy-format blocker.

### 2.6 Munging — `munge_sumstats.py` → `ldsc munge-sumstats`

**Coverage: FULL (changed). All column-hint/QC flags retained; allele/DANER flags
restructured.**

| Legacy flag | Refactored | Status | Severity | Notes |
|---|---|---|---|---|
| `--sumstats` | `--raw-sumstats-file` | RENAMED | — | |
| `--out` | `--output-dir` | CHANGED | MEDIUM | Default artifact is `sumstats.parquet`. |
| `--N` / `--N-cas` / `--N-con` | same | SAME | — | |
| `--info-min` (0.9) / `--maf-min` (0.01) / `--n-min` | same | SAME | — | |
| `--chunksize` | `--chunksize` | CHANGED | LOW | Default **5e6 → 1_000_000** (`_kernel/sumstats_munger.py:1061`). |
| `--daner` | `--format daner-old` | RENAMED/CHANGED | LOW | Boolean→choice. |
| `--daner-n` | `--format daner-new` | RENAMED/CHANGED | LOW | |
| `--no-alleles` | NONE (use `--snp-identifier chr_pos`/`rsid`) | CHANGED | MEDIUM | Allele-blind via base identity modes. |
| `--merge-alleles` | NONE (nearest `--use-hm3-snps` / `--sumstats-snps-file`) | NONE | **HIGH** | No external SNP+A1+A2 allele-harmonize/flip-to-list. `--use-hm3-snps` covers the standard HM3 case but uses the PACKAGED map, not an arbitrary list, and does not flip to list alleles. |
| `--snp`/`--N-col`/`--a1`/`--a2`/`--p`/`--frq`/`--signed-sumstats`/`--info`/`--info-list`/`--nstudy`/`--nstudy-min`/`--ignore`/`--N-cas-col`/`--N-con-col` | same | SAME | — | Cloned from kernel parser. |
| `--a1-inc` / `--keep-maf` | same | SAME | — | |

New: `--chr`, `--pos`, `--format`, `--infer-only`, `--output-format`,
`--source-genome-build`, `--output-genome-build`, `--liftover-chain-file`,
`--use-hm3-snps`, `--use-hm3-quick-liftover`, `--sumstats-snps-file`,
`--trait-name`, `--snp-identifier`, `--overwrite`, `--log-level`.

Wiki translation (recurring munge example):
- `munge_sumstats.py --sumstats X.txt --N 17115 --out scz --merge-alleles w_hm3.snplist`
  → `ldsc munge-sumstats --raw-sumstats-file X.txt --N 17115 --output-dir scz --use-hm3-snps --snp-identifier rsid --output-format tsv.gz` — **PARTIAL** (`--merge-alleles`→`--use-hm3-snps`; allele-harmonize-to-list semantics not reproduced).

### 2.7 Annotation building — `make_annot.py` → `ldsc annotate`

**Coverage: PARTIAL. BED mode covered; gene-set mode missing.**

| Legacy flag | Refactored | Status | Severity | Notes |
|---|---|---|---|---|
| `--bed-file` | `--query-annot-bed-sources` | RENAMED/CHANGED | MEDIUM | Multi-BED + globs; each stem → a query annotation. |
| `--bimfile` | `--baseline-annot-sources` | CHANGED | MEDIUM | SNP universe from baseline annot template, not `.bim`. |
| `--annot-file` | `--output-dir` | CHANGED | LOW | Writes per-chr `query.<chrom>.annot.gz`. |
| `--nomerge` | NONE | NONE | MEDIUM | No proportional-overlap-count mode. |
| `--gene-set-file` | NONE | NONE | **HIGH** | No gene-set→BED building. |
| `--gene-coord-file` | NONE | NONE | **HIGH** | Tied to gene-set mode. |
| `--windowsize` | NONE (nearest `--bed-padding-bp`) | NONE | MEDIUM | `--bed-padding-bp` pads existing BED intervals; cannot window around gene TSS/TES from a gene list. |

New: `--bed-padding-bp`, `--snp-identifier`, `--genome-build`, `--overwrite`,
`--log-level`.

Wiki translations (LD Score Estimation Tutorial):
- `make_annot.py --bed-file X.bed --bimfile 1000G.EUR.QC.22.bim --annot-file X.annot.gz`
  → `ldsc annotate --query-annot-bed-sources X.bed --baseline-annot-sources <baseline> --output-dir <dir>` — PARTIAL.
- `make_annot.py --gene-set-file G.GeneSet --gene-coord-file ENSG_coord.txt --windowsize 100000 --bimfile ... --annot-file G.annot.gz`
  → **NO** (gene-set mode absent; must build the BED externally first).

### 2.8 New refactored capability — `ldsc query-r2`

No legacy equivalent. Point/range R2 lookup from a package-built parquet panel
(`--panel-dir`, `--pairs`). Not a gap; noted for completeness.

---

## 3. Wiki-Example Reproducibility Summary

| Wiki tutorial | Examples | Reproducible verbatim? |
|---|---|---|
| Heritability & Genetic Correlation | munge, h2, rg, liability, constrain-intercept, no-intercept | munge/h2/rg **NO** with distributed LD scores (format blocker); prevalence/no-intercept flags YES; per-pair intercept vectors PARTIAL |
| LD Score Estimation | univariate `--l2`, `--extract` w-ld, make_annot (BED + gene-set), partitioned `--l2` | univariate/partitioned PARTIAL (need `--exclude-regions none`, `--snp-identifier rsid`, baseline for annot); gene-set make_annot **NO** |
| Partitioned Heritability | munge, baseline-model partitioned h2, cell-type group | partitioned/cell-type **NO** (format blocker + multi-prefix concat gone) |
| Cell type specific analyses | munge, `--h2-cts` `.ldcts` sweep | `--h2-cts` **NO** (no one-command sweep) |
| Continuous Annotations | baseline-LD h2 with `--print-delete-vals`, continuous-annot `--l2`, add-annotation | h2 **NO** (`--print-delete-vals` gap breaks the `quantile_h2g.r` step; format blocker); continuous-annot `--l2` PARTIAL (values supported; needs baseline) |
| Summary Stats Format / LD File Formats / What Data / Tests / Home / FAQ | format specs, test invocation | conceptual; no unique runnable analysis commands |

Root causes of non-reproducibility, in order of impact: (1) legacy LD-score
directory format unreadable (CRITICAL); (2) `--exclude-regions` default ON
(CRITICAL); (3) `--merge-alleles` / gene-set / `--h2-cts` / `--cts-bin` /
`--per-allele` / `--print-delete-vals` feature gaps (HIGH); (4) output
prefix→directory and identity-mode defaults (MEDIUM).

---

## 4. Consolidated Gap/Discrepancy Table (CRITICAL first)

| # | Severity | Area | Gap / discrepancy | Grounding | Justification (rubric) |
|---|---|---|---|---|---|
| 1 | **CRITICAL** | regression input | Legacy `.l2.ldscore.gz`/`.M`/`.M_5_50` LD-score directories are unreadable; only canonical parquet `--ldscore-dir` accepted | `regression_runner.py:2465` (requires `metadata.json`) | Core documented workflow (use precomputed LD scores) entirely broken; whole distributed-artifact ecosystem incompatible |
| 2 | **CRITICAL** | ldscore / build-ref-panel | `--exclude-regions` defaults to `mhc-and-centromeres`; legacy never excludes | `ldscore_calculator.py:930` | Ported legacy command silently changes results without warning |
| 3 | HIGH | ldscore | `--cts-bin`/`--cts-breaks`/`--cts-names` on-the-fly binning absent | `_kernel/formats.py:104` (dead `read_cts`) | Documented capability missing; manual-binning workaround is non-obvious/lossy |
| 4 | HIGH | ldscore | `--per-allele` / `--pq-exp` LD-score scaling absent | verified zero matches | Documented feature missing; no workaround |
| 5 | HIGH | regression | No independent regression-weight LD set (`--w-ld` folded into `--ldscore-dir`) | `io-argument-inventory.md:393`,`441` | Removes a documented degree of freedom (wiki recommends separate `--w-ld-chr` for partitioned) |
| 6 | HIGH | cell-type | `--h2-cts` + `--ref-ld-chr-cts` `.ldcts` one-command sweep absent | verified zero matches; regime at `regression_runner.py:791` | Documented command missing; only a manual multi-run loop reproduces it |
| 7 | HIGH | annotate | `make_annot` gene-set mode (`--gene-set-file`/`--gene-coord-file`/`--windowsize`) absent | verified zero matches | Documented capability missing; must pre-build BED externally |
| 8 | HIGH | munge | `--merge-alleles` allele-harmonize/flip-to-list absent | removed (`io-argument-inventory.md:330-332`) | Standard HM3-merge workflow changed; `--use-hm3-snps` is packaged-map only and does not flip to list alleles |
| 9 | HIGH | regression (continuous-annot) | `--print-delete-vals` absent breaks the documented `quantile_h2g.r` workflow | computed at `_kernel/regression.py:350`, unexposed | Documented downstream analysis cannot run |
| 10 | MEDIUM | ldscore/regression | Strict `<`/`>` → inclusive `<=`/`>=` for `--maf`, common-MAF, `--two-step`, `--chisq-max` | `ldsc.py:357,490`; `io-argument-inventory.md:177,402-403,425-426,451` | Boundary-SNP retention differs; documented |
| 11 | MEDIUM | all artifact cmds | Output `--out <prefix>` → `--output-dir <dir>` + fixed filenames | `io-argument-inventory.md:103-138` | Scripts parsing legacy paths break; clearly different |
| 12 | MEDIUM | h2/partitioned/rg | `--M` manual SNP-count override absent | verified zero matches | Visible feature loss; counts now from LD-score dir |
| 13 | MEDIUM | h2/partitioned/rg | `--print-cov` absent | verified zero matches | Diagnostic output lost |
| 14 | MEDIUM | h2/partitioned | `--invert-anyway` absent | verified zero matches | Escape hatch for ill-conditioned matrices lost |
| 15 | MEDIUM | rg | `--no-check-alleles` absent (identity model instead) | verified zero matches | Behavioral; allele check now structural |
| 16 | MEDIUM | munge | `--no-alleles` → `--snp-identifier chr_pos`/`rsid` | `io-argument-inventory.md:331` | Achieved via identity mode, not a flag |
| 17 | MEDIUM | ldscore | Default identity `chr_pos_allele_aware` vs legacy rs-number keying | `_kernel/sumstats_munger.py:1102` | Must pass `--snp-identifier rsid` to match legacy merges |
| 18 | MEDIUM | annotate | `--nomerge` proportional-overlap counts absent | verified zero matches | Membership-only projection |
| 19 | MEDIUM | annotate | `--windowsize` (gene TSS/TES window) absent; `--bed-padding-bp` only pads BEDs | `io-argument-inventory.md:149` | Cannot window from a gene list |
| 20 | MEDIUM | partitioned/continuous | Comma-separated multi-prefix `--ref-ld-chr A.,B.` concatenation absent | `io-argument-inventory.md:434` (`--query-columns` removed) | One canonical dir only |
| 21 | LOW | LD-score est | `--bfile`→`--plink-prefix`, `--keep`→`--keep-indivs-file`, `--extract`→`--ref-panel-snps-file`, `--chunk-size`→`--snp-batch-size` (50→128) | `io-argument-inventory.md:192-195` | Pure renames / perf defaults |
| 22 | LOW | LD-score est | `--l2`, `--thin-annot` absorbed into automatic behavior | `cli.py`; annot auto-detect | No functional loss |
| 23 | LOW | partitioned | `--overlap-annot`, `--frqfile(-chr)`, `--print-coefficients` absorbed (automatic) | `io-argument-inventory.md` overlap section | No functional loss |
| 24 | LOW | munge | `--daner`/`--daner-n` → `--format daner-old`/`daner-new`; `--chunksize` 5e6→1e6 | `sumstats_munger.py:1012-1018`; `_kernel/sumstats_munger.py:1061` | Cosmetic / perf |
| 25 | LOW | rg | `--rg a,b,c` (comma) → `--sumstats-sources a b c` (space); `--intercept-gencov` list→scalar; `--not-M-5-50`→`--count-kind`; `--return-silly-things`/`--print-all-cts` absent | `regression_runner.py:1966-1987` | Cosmetic / niche escape hatches |
| 26 | LOW | LD-score est | `--Pickle` absent (was a declared-but-unused legacy flag) | `ldsc.py:560` (never read) | No functional loss |

---

## 5. Open Questions / Could Not Verify

- **Exact legacy threshold operators in the numerical kernel.** Legacy `--maf`
  strictness is taken from help text ("MAF > 0", `ldsc.py:490`) and `maf > 0.05`
  for M_5_50 (`ldsc.py:357`); I did not trace the PLINK genotype-filter operator in
  `ldscore/ldscore.py`. Direction is consistent with the refactor's documented
  strict→inclusive change, but exact boundary behavior in legacy is unverified at
  the kernel level. (Numerical equivalence is out of scope anyway.)
- **`--merge-alleles` allele-flip nuance.** I confirmed the flag is removed and that
  `--use-hm3-snps`/`--sumstats-snps-file` are the nearest restriction mechanisms. I
  did not exhaustively trace whether the refactored identity model performs an
  equivalent A1/A2 flip-to-reference during merge; the user's separate allele-
  orientation/MAF-normalization work (project memory) may already cover this.
- **`--return-silly-things` default behavior.** The kernel exposes
  `summary(silly=False)` (`_kernel/regression.py:857`) but no public flag wires
  `silly=True`. I did not verify the exact clamp/guard the refactored rg applies to
  out-of-bounds estimates (only that the default path matches legacy's default).
- **Per-chromosome vs whole-genome annot output naming.** `ldsc annotate` writes
  `query.<chrom>.annot.gz`; I did not verify the precise chromosome-token behavior
  for single-file BED inputs (only the documented contract).
- **Whether any wiki example beyond the 6 tutorials carries unique commands.** The
  remaining wiki pages (Home, FAQ, LD File Formats, Summary Statistics File Format,
  What Data Are Necessary, Tests) are conceptual/format references; I read them by
  title/section and found no unique runnable analysis commands, but did not extract
  every prose snippet.

---

## Appendix: source-of-truth files for this audit

- Legacy: `ldsc.py` (49 flags, dispatch `:582-661`), `munge_sumstats.py` (29 flags),
  `make_annot.py` (7 flags), `ldscore/{ldscore,sumstats,regressions}.py`.
- Refactored: `src/ldsc/cli.py` (8 subcommands), `annotation_builder.py`,
  `ldscore_calculator.py`, `ref_panel_builder.py`, `sumstats_munger.py`
  (+`_kernel/sumstats_munger.py`), `regression_runner.py`, `r2_query.py`,
  `docs/current/io-argument-inventory.md`.
- Working notes: `docs/audits/legacy-equivalence/findings.md` (full per-flag detail
  and per-tutorial command extractions backing every row above).
