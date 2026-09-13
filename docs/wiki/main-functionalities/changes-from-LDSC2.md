# Moving from LDSC2

Last updated on: 2026-09-13

Build commands from the current `ldsc COMMAND --help` and workflow examples. Use the [legacy CLI flag map](../../current/legacy-cli-flag-map.md) to check replacements and removed flags.

| Topic | Current contract |
| --- | --- |
| SNP identity | The general default is `chr_pos_allele_aware`. Base `chr_pos`, `rsid`, and `rsid_allele_aware` are explicit alternatives where supported. Use matching identities across inputs. |
| Identity exceptions | `build-gene-ldscore-index` requires explicit base `rsid` or `chr_pos` with hg19. Indexed LD scoring inherits the index settings. Legacy LD-score conversion defaults to `rsid` and supports only base identity modes. |
| Baseline annotations | Direct query LD scoring requires explicit baseline annotations. Only ordinary unpartitioned scoring can synthesize an all-ones `base`. Indexed scoring uses the stored baseline. |
| Query annotations | Supply prebuilt annotation columns, BED files, or gene lists as mutually exclusive routes. Direct BED/gene projection is integrated into `ldscore`; standalone `annotate` optionally creates reusable query files. |
| LD-score artifacts | Regression reads a canonical directory with `metadata.json`, `ldscore.baseline.parquet`, and applicable query/overlap Parquet files. It does not recompute LD scores. |
| Legacy sumstats | Regression accepts `.sumstats`/`.sumstats.gz` with `SNP`, `A1`, `A2`, `Z`, and `N` directly. See [compatibility requirements](../../current/legacy-sumstats-compatibility.md). |
| Legacy LD scores | Use [explicit conversion](../utility-functionalities/convert-ldsc2-ldscores.md) for complete chromosome 1–22 unpartitioned or baseline-only suites. Precomputed scores do not replace the baseline annotations needed for new query LD scoring. |
| Regression SNPs | New LD scoring defaults to packaged HapMap3 candidates, followed by MHC-and-centromere exclusion from output/regression rows. This does not restrict reference contributors or annotation counts. Munging independently defaults to HM3 restriction. |
| Per-query results | `partitioned-h2` automatically writes one complete baseline-plus-query result per query. The removed `--write-per-query-results` flag is not accepted. |

See [global configuration](global-config.md), [LD scoring](ldscore.md), and [partitioned heritability](partitioned-h2.md) for examples. Sources: `GlobalConfig` in [config.py](../../../src/ldsc/config.py), `build_parser` in [ldscore_calculator.py](../../../src/ldsc/ldscore_calculator.py), and `add_partitioned_h2_arguments` in [regression_runner.py](../../../src/ldsc/regression_runner.py).
