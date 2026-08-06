# Exact gene LD-score index chromosome-22 gate

Date: 2026-08-03
Last updated on: 2026-08-03

Status: passed locally for both approved baseline suites. This gate permits a
separate production whole-genome construction; no production whole-genome suite
was built by this implementation task.

## Inputs and settings

- PLINK: local `1000G.EUR.QC.22` bed/bim/fam suite.
- Baselines: `1000G_EUR_Phase3_baseline/baseline.22.annot.gz` and
  `1000G_Phase3_baselineLD_v2.2_ldscores/baselineLD.22.annot.gz`.
- Profile: hg19, rsID, 1 cM, bundled HM3 minus `mhc-and-centromeres`, no
  explicit retained-reference MAF filter, inclusive common MAF 0.05, 100 kb
  padding, MHC gene exclusion, one worker, SNP batch 128, atom batch 64.
- Queries: overlapping/nested chromosome-22 catalog genes, repeated identifiers,
  a versioned Ensembl alias, and an adjacent overlapping set. The default
  all-protein-coding control was included.

Both baseline sources and the BIM contained exactly 141,123 pre-QC
`CHR/POS/SNP` rows. Genotype preparation retained all 141,123 rows; bundled HM3
minus region exclusions produced 17,380 persisted regression rows.

## Construction measurements

| Baseline suite | Wall time | Peak RSS | Atoms | `nnz(Y)` | CSR operator bytes | Payload bytes |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| baseline | 103.40 s | 1,533,739,008 bytes | 838 | 614,523 | 7,443,800 | 10,533,234 |
| baselineLD v2.2 | 119.45 s | 1,758,478,336 bytes | 838 | 614,523 | 7,443,800 | 14,992,166 |

The final measurements come from the builder's own diagnostics after publishing
and reloading each artifact with the current strict loader. The baseline process
was also timed with macOS `/usr/bin/time -l`; its scientific build completed,
but sandbox denial of the tool's final `sysctl` made that command return a
nonzero status after publication. Both final artifacts independently reload.
Each selected all 489 FAM individuals in genotype order and recorded the same
selected-content SHA-256. Both measured peaks are below the 4–8 GB target.

## Direct versus indexed results

For both suites, regression rows and ordering, baseline/query columns,
canonical Parquet schemas, count records, statuses, and diagnostics matched.
Every final baseline/control/query LD-score value had absolute difference 0
after the normal output cast. Baseline-suite overlaps were exactly identical;
baselineLD v2.2 continuous-annotation overlap accumulation had maximum absolute
difference `7.654307410120964e-08`, with identical labels and shape.

The empirical acceptance tolerance is therefore absolute `1e-7` for overlap
sufficient-statistic accumulation and exact equality after canonical output
casting for LD-score columns under this configuration. A larger discrepancy is
not accepted without investigation.

A deterministic synthetic chromosome-22 sumstats table with all 17,380 rows was
fed to unchanged `partitioned-h2` for matching direct and indexed baseline
outputs. Category ordering, coefficients, standard errors, enrichment,
p-values, and all other numeric fields were exactly identical (maximum absolute
difference 0).

## Remaining external evidence

Only chromosome 22 of the approved PLINK suite is available locally.
Deterministic unit tests cover chromosome-6 separation of gene MHC filtering
from the broad LD-score contributor universe. A full chromosome-6 resource
measurement and production chromosomes 1–22 construction remain HPC/distribution
operations, not prerequisites hidden by this local gate.

## Human-readable logging gate

The fresh baseline chromosome-22 run at
`/private/tmp/codex_gene_index_logging_fresh_chr22_20260803_baseline` confirms
that the operational log is written through shared workflow logging and is
reloaded beside the structured JSON summary:

- [fresh log](/private/tmp/codex_gene_index_logging_fresh_chr22_20260803_baseline/profiles/padding-100000bp-mhc/diagnostics/build-gene-ldscore-index.log)
- [fresh JSON](/private/tmp/codex_gene_index_logging_fresh_chr22_20260803_baseline/profiles/padding-100000bp-mhc/diagnostics/build-gene-ldscore-index.json)

The log contains the shared lifecycle header/footer, resolved hg19/rsID and
1 cM configuration, baseline/PLINK/map/sample resolution counts, the broad
retained-PLINK versus filtered-HM3 SNP-universe explanation, and explicit
`Starting chromosome 22`/`Finished chromosome 22` progress. The chromosome
summary records 426 protein-coding genes after MHC exclusion and 100 kb padding,
141,123 pre-QC and retained-reference rows, 17,380 regression rows, 838 atoms,
and `nnz(Y)=614,523`. The JSON carries the same values; its log is 6,190 bytes
versus 454 bytes for the previous summary-only log, while the scientific build
completed successfully and the final profile reloaded.

Representative excerpt:

```text
Starting chromosome 22.
Finished chromosome 22: protein-coding genes=426, retained-reference=141123,
regression-rows=17380, atoms=838, nnz(Y)=614523, elapsed=104.858s.
SNP universe: broad retained PLINK SNPs are LD-score contributors and count/overlap members;
filtered HM3 SNPs are persisted regression rows; w_ld uses the filtered regression set as rows and contributors.
```

The focused failure test also confirms that strict baseline/BIM failures retain
the `Failed` footer, phase, exception, and traceback without publishing suite or
profile metadata; diagnostics-only output is recoverable by a subsequent retry.
