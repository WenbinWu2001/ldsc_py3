# Complete CLI Flag Inventory

Last updated on: 2026-09-11

Snapshot of the working tree on `restructure` (base commit `69081b2`). Excludes built-in help; includes hidden aliases. All 218 command-option entries are listed. Counts in D/T are the number of files containing the literal flag in current documentation/tutorials and tests, across all commands. They are navigation evidence, not execution telemetry or proof of missing Python-API coverage. Full file lists, defaults, choices, help text, and declaration locations are in [inventory.json](inventory.json).

Recommendations are proposals; package behavior has not changed. See [the audit](README.md) for severity, scientific distinctions, and reproduction evidence.

## annotate (12)

| Flag | Default / required | Disposition | D/T | Declaration |
| --- | --- | --- | --- | --- |
| `--query-annot-bed-sources` | `None` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 12/3 | [annotation_builder.py:86](../../../src/ldsc/annotation_builder.py) |
| `--query-annot-gene-list-sources` | `None` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 12/5 | [annotation_builder.py:91](../../../src/ldsc/annotation_builder.py) |
| `--gene-coordinate-file` | `None` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 11/5 | [annotation_builder.py:92](../../../src/ldsc/annotation_builder.py) |
| `--gene-list-resolution-policy` | `None` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 8/2 | [annotation_builder.py:93](../../../src/ldsc/annotation_builder.py) |
| `--gene-exclude-regions` | `None` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 9/3 | [annotation_builder.py:94](../../../src/ldsc/annotation_builder.py) |
| `--baseline-annot-sources` | `required` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 19/8 | [annotation_builder.py:95](../../../src/ldsc/annotation_builder.py) |
| `--output-dir` | `required` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 37/18 | [annotation_builder.py:101](../../../src/ldsc/annotation_builder.py) |
| `--padding-bp` | `None` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 17/5 | [annotation_builder.py:106](../../../src/ldsc/annotation_builder.py) |
| `--overwrite` | `False` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 36/10 | [annotation_builder.py:112](../../../src/ldsc/annotation_builder.py) |
| `--snp-identifier` | `chr_pos_allele_aware` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 26/10 | [annotation_builder.py:118](../../../src/ldsc/annotation_builder.py) |
| `--genome-build` | `None` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 24/12 | [annotation_builder.py:124](../../../src/ldsc/annotation_builder.py) |
| `--log-level` | `INFO` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 7/2 | [annotation_builder.py:134](../../../src/ldsc/annotation_builder.py) |

## ldscore (34)

| Flag | Default / required | Disposition | D/T | Declaration |
| --- | --- | --- | --- | --- |
| `--output-dir` | `required` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 37/18 | [ldscore_calculator.py:937](../../../src/ldsc/ldscore_calculator.py) |
| `--gene-ldscore-index-dir` | `None` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 9/2 | [ldscore_calculator.py:938](../../../src/ldsc/ldscore_calculator.py) |
| `--overwrite` | `False` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 36/10 | [ldscore_calculator.py:943](../../../src/ldsc/ldscore_calculator.py) |
| `--query-annot-sources` | `None` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 7/1 | [ldscore_calculator.py:950](../../../src/ldsc/ldscore_calculator.py) |
| `--query-annot-bed-sources` | `None` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 12/3 | [ldscore_calculator.py:955](../../../src/ldsc/ldscore_calculator.py) |
| `--query-annot-gene-list-sources` | `None` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 12/5 | [ldscore_calculator.py:960](../../../src/ldsc/ldscore_calculator.py) |
| `--padding-bp` | `None` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 17/5 | [ldscore_calculator.py:969](../../../src/ldsc/ldscore_calculator.py) |
| `--gene-coordinate-file` | `None` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 11/5 | [ldscore_calculator.py:980](../../../src/ldsc/ldscore_calculator.py) |
| `--gene-list-resolution-policy` | `strict` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 8/2 | [ldscore_calculator.py:988](../../../src/ldsc/ldscore_calculator.py) |
| `--gene-exclude-regions` | `none` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 9/3 | [ldscore_calculator.py:997](../../../src/ldsc/ldscore_calculator.py) |
| `--control-gene-list-file` | `None` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 11/3 | [ldscore_calculator.py:1003](../../../src/ldsc/ldscore_calculator.py) |
| `--baseline-annot-sources` | `None` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 19/8 | [ldscore_calculator.py:1008](../../../src/ldsc/ldscore_calculator.py) |
| `--plink-prefix` | `None` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 12/7 | [ldscore_calculator.py:1013](../../../src/ldsc/ldscore_calculator.py) |
| `--r2-dir` | `None` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 13/3 | [ldscore_calculator.py:1021](../../../src/ldsc/ldscore_calculator.py) |
| `--snp-identifier` | `chr_pos_allele_aware` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 26/10 | [ldscore_calculator.py:1027](../../../src/ldsc/ldscore_calculator.py) |
| `--genome-build` | `None` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 24/12 | [ldscore_calculator.py:1033](../../../src/ldsc/ldscore_calculator.py) |
| `--ref-panel-snps-file` | `None` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 12/1 | [ldscore_calculator.py:1045](../../../src/ldsc/ldscore_calculator.py) |
| `--regr-snps-exclude-regions` | `mhc-and-centromeres` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 10/7 | [ldscore_calculator.py:1054](../../../src/ldsc/ldscore_calculator.py) |
| `--exclude-regions` | `suppressed alias` | Remove: hidden alias of --regr-snps-exclude-regions. | 4/3 | [ldscore_calculator.py:1060](../../../src/ldsc/ldscore_calculator.py) |
| `--regr-snps-file` | `None` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 9/6 | [ldscore_calculator.py:1067](../../../src/ldsc/ldscore_calculator.py) |
| `--keep-indivs-file` | `None` | Rename: --keep-individuals-file; preserve IID contract. | 5/1 | [ldscore_calculator.py:1075](../../../src/ldsc/ldscore_calculator.py) |
| `--ld-wind-snps` | `None` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 4/4 | [ldscore_calculator.py:1080](../../../src/ldsc/ldscore_calculator.py) |
| `--ld-wind-kb` | `None` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 4/3 | [ldscore_calculator.py:1081](../../../src/ldsc/ldscore_calculator.py) |
| `--ld-wind-cm` | `None` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 17/3 | [ldscore_calculator.py:1082](../../../src/ldsc/ldscore_calculator.py) |
| `--maf-min` | `None` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 6/3 | [ldscore_calculator.py:1083](../../../src/ldsc/ldscore_calculator.py) |
| `--common-maf-min` | `0.05` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 11/3 | [ldscore_calculator.py:1084](../../../src/ldsc/ldscore_calculator.py) |
| `--genetic-map-hg19-sources` | `None` | Keep PLINK map; ignored by parquet backend (sidecar CM authoritative). | 7/3 | [ldscore_calculator.py:1085](../../../src/ldsc/ldscore_calculator.py) |
| `--genetic-map-hg38-sources` | `None` | Keep PLINK map; ignored by parquet backend (sidecar CM authoritative). | 6/3 | [ldscore_calculator.py:1086](../../../src/ldsc/ldscore_calculator.py) |
| `--export-ref-metadata` | `False` | Keep PLINK export; no export in parquet-sidecar mode. | 5/2 | [ldscore_calculator.py:1087](../../../src/ldsc/ldscore_calculator.py) |
| `--snp-batch-size` | `128` | Keep PLINK batching; ineffective in ldscore parquet mode. | 4/2 | [ldscore_calculator.py:1088](../../../src/ldsc/ldscore_calculator.py) |
| `--query-batch-size` | `1000` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 9/1 | [ldscore_calculator.py:1089](../../../src/ldsc/ldscore_calculator.py) |
| `--threads` | `1` | Rename --workers: processes for ldscore, threads for gene-index builder. | 8/2 | [ldscore_calculator.py:1090](../../../src/ldsc/ldscore_calculator.py) |
| `--yes-really` | `False` | Rename: --allow-whole-chromosome-window. | 3/3 | [ldscore_calculator.py:1091](../../../src/ldsc/ldscore_calculator.py) |
| `--log-level` | `INFO` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 7/2 | [ldscore_calculator.py:1092](../../../src/ldsc/ldscore_calculator.py) |

## build-ref-panel (18)

| Flag | Default / required | Disposition | D/T | Declaration |
| --- | --- | --- | --- | --- |
| `--plink-prefix` | `required` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 12/7 | [ref_panel_builder.py:1632](../../../src/ldsc/ref_panel_builder.py) |
| `--source-genome-build` | `auto` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 14/3 | [ref_panel_builder.py:1640](../../../src/ldsc/ref_panel_builder.py) |
| `--genetic-map-hg19-sources` | `None` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 7/3 | [ref_panel_builder.py:1646](../../../src/ldsc/ref_panel_builder.py) |
| `--genetic-map-hg38-sources` | `None` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 6/3 | [ref_panel_builder.py:1651](../../../src/ldsc/ref_panel_builder.py) |
| `--liftover-chain-hg19-to-hg38-file` | `None` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 2/0 | [ref_panel_builder.py:1656](../../../src/ldsc/ref_panel_builder.py) |
| `--liftover-chain-hg38-to-hg19-file` | `None` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 3/2 | [ref_panel_builder.py:1664](../../../src/ldsc/ref_panel_builder.py) |
| `--output-dir` | `required` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 37/18 | [ref_panel_builder.py:1672](../../../src/ldsc/ref_panel_builder.py) |
| `--overwrite` | `False` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 36/10 | [ref_panel_builder.py:1673](../../../src/ldsc/ref_panel_builder.py) |
| `--ld-wind-snps` | `None` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 4/4 | [ref_panel_builder.py:1682](../../../src/ldsc/ref_panel_builder.py) |
| `--ld-wind-kb` | `None` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 4/3 | [ref_panel_builder.py:1683](../../../src/ldsc/ref_panel_builder.py) |
| `--ld-wind-cm` | `None` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 17/3 | [ref_panel_builder.py:1684](../../../src/ldsc/ref_panel_builder.py) |
| `--maf-min` | `None` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 6/3 | [ref_panel_builder.py:1685](../../../src/ldsc/ref_panel_builder.py) |
| `--ref-panel-snps-file` | `None` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 12/1 | [ref_panel_builder.py:1686](../../../src/ldsc/ref_panel_builder.py) |
| `--snp-identifier` | `chr_pos_allele_aware` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 26/10 | [ref_panel_builder.py:1695](../../../src/ldsc/ref_panel_builder.py) |
| `--keep-indivs-file` | `None` | Rename: --keep-individuals-file; preserve IID contract. | 5/1 | [ref_panel_builder.py:1705](../../../src/ldsc/ref_panel_builder.py) |
| `--snp-batch-size` | `128` | Keep PLINK batching; ineffective in ldscore parquet mode. | 4/2 | [ref_panel_builder.py:1706](../../../src/ldsc/ref_panel_builder.py) |
| `--min-r2` | `0.0` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 3/0 | [ref_panel_builder.py:1713](../../../src/ldsc/ref_panel_builder.py) |
| `--log-level` | `INFO` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 7/2 | [ref_panel_builder.py:1726](../../../src/ldsc/ref_panel_builder.py) |

## build-gene-ldscore-index (21)

| Flag | Default / required | Disposition | D/T | Declaration |
| --- | --- | --- | --- | --- |
| `--baseline-annot-sources` | `required` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 19/8 | [gene_ldscore_index.py:124](../../../src/ldsc/gene_ldscore_index.py) |
| `--plink-prefix` | `required` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 12/7 | [gene_ldscore_index.py:125](../../../src/ldsc/gene_ldscore_index.py) |
| `--output-dir` | `required` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 37/18 | [gene_ldscore_index.py:126](../../../src/ldsc/gene_ldscore_index.py) |
| `--gene-coordinate-file` | `required` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 11/5 | [gene_ldscore_index.py:127](../../../src/ldsc/gene_ldscore_index.py) |
| `--genome-build` | `required` | Keep explicit hg19 assertion under current contract; single allowed value is deliberate. | 24/12 | [gene_ldscore_index.py:132](../../../src/ldsc/gene_ldscore_index.py) |
| `--snp-identifier` | `required` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 26/10 | [gene_ldscore_index.py:138](../../../src/ldsc/gene_ldscore_index.py) |
| `--padding-bp` | `0` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 17/5 | [gene_ldscore_index.py:144](../../../src/ldsc/gene_ldscore_index.py) |
| `--gene-exclude-regions` | `mhc` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 9/3 | [gene_ldscore_index.py:145](../../../src/ldsc/gene_ldscore_index.py) |
| `--ld-wind-cm` | `1.0` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 17/3 | [gene_ldscore_index.py:146](../../../src/ldsc/gene_ldscore_index.py) |
| `--maf-min` | `None` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 6/3 | [gene_ldscore_index.py:147](../../../src/ldsc/gene_ldscore_index.py) |
| `--common-maf-min` | `0.05` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 11/3 | [gene_ldscore_index.py:148](../../../src/ldsc/gene_ldscore_index.py) |
| `--keep-indivs-file` | `None` | Rename: --keep-individuals-file; preserve IID contract. | 5/1 | [gene_ldscore_index.py:149](../../../src/ldsc/gene_ldscore_index.py) |
| `--regr-snps-file` | `None` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 9/6 | [gene_ldscore_index.py:150](../../../src/ldsc/gene_ldscore_index.py) |
| `--regr-snps-exclude-regions` | `mhc-and-centromeres` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 10/7 | [gene_ldscore_index.py:158](../../../src/ldsc/gene_ldscore_index.py) |
| `--exclude-regions` | `suppressed alias` | Remove: hidden alias of --regr-snps-exclude-regions. | 4/3 | [gene_ldscore_index.py:167](../../../src/ldsc/gene_ldscore_index.py) |
| `--genetic-map-hg19-sources` | `None` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 7/3 | [gene_ldscore_index.py:174](../../../src/ldsc/gene_ldscore_index.py) |
| `--snp-batch-size` | `128` | Keep PLINK batching; ineffective in ldscore parquet mode. | 4/2 | [gene_ldscore_index.py:175](../../../src/ldsc/gene_ldscore_index.py) |
| `--atom-batch-size` | `64` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 1/0 | [gene_ldscore_index.py:176](../../../src/ldsc/gene_ldscore_index.py) |
| `--threads` | `1` | Rename --workers: processes for ldscore, threads for gene-index builder. | 8/2 | [gene_ldscore_index.py:177](../../../src/ldsc/gene_ldscore_index.py) |
| `--overwrite` | `False` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 36/10 | [gene_ldscore_index.py:178](../../../src/ldsc/gene_ldscore_index.py) |
| `--log-level` | `INFO` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 7/2 | [gene_ldscore_index.py:179](../../../src/ldsc/gene_ldscore_index.py) |

## convert-ldsc2-ldscores (8)

| Flag | Default / required | Disposition | D/T | Declaration |
| --- | --- | --- | --- | --- |
| `--legacy-reference-dir` | `required` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 2/1 | [legacy_ldscore_converter.py:689](../../../src/ldsc/legacy_ldscore_converter.py) |
| `--legacy-weight-dir` | `required` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 2/1 | [legacy_ldscore_converter.py:690](../../../src/ldsc/legacy_ldscore_converter.py) |
| `--legacy-frequency-dir` | `None` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 2/0 | [legacy_ldscore_converter.py:691](../../../src/ldsc/legacy_ldscore_converter.py) |
| `--output-dir` | `required` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 37/18 | [legacy_ldscore_converter.py:692](../../../src/ldsc/legacy_ldscore_converter.py) |
| `--snp-identifier` | `rsid` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 26/10 | [legacy_ldscore_converter.py:693](../../../src/ldsc/legacy_ldscore_converter.py) |
| `--genome-build` | `auto` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 24/12 | [legacy_ldscore_converter.py:694](../../../src/ldsc/legacy_ldscore_converter.py) |
| `--overwrite` | `False` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 36/10 | [legacy_ldscore_converter.py:695](../../../src/ldsc/legacy_ldscore_converter.py) |
| `--log-level` | `INFO` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 7/2 | [legacy_ldscore_converter.py:696](../../../src/ldsc/legacy_ldscore_converter.py) |

## munge-sumstats (40)

| Flag | Default / required | Disposition | D/T | Declaration |
| --- | --- | --- | --- | --- |
| `--raw-sumstats-file` | `required` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 13/3 | [sumstats_munger.py:980](../../../src/ldsc/sumstats_munger.py) |
| `--output-dir` | `required` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 37/18 | [sumstats_munger.py:981](../../../src/ldsc/sumstats_munger.py) |
| `--format` | `auto` | Repair auto/explicit DANER divergence; rename --input-format; explicit daner-new is a consolidation candidate. | 11/1 | [sumstats_munger.py:982](../../../src/ldsc/sumstats_munger.py) |
| `--infer-only` | `False` | Keep; repair false runnable reports and incomplete suggested commands. | 16/2 | [sumstats_munger.py:989](../../../src/ldsc/sumstats_munger.py) |
| `--overwrite` | `False` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 36/10 | [sumstats_munger.py:995](../../../src/ldsc/sumstats_munger.py) |
| `--sumstats-snps-file` | `None` | Keep; name could become --keep-snps-file, consistent with other universes. | 5/2 | [sumstats_munger.py:1001](../../../src/ldsc/sumstats_munger.py) |
| `--use-hm3-snps` | `False` | Keep; packaged convenience filter, distinct from liftover. | 13/3 | [sumstats_munger.py:1009](../../../src/ldsc/sumstats_munger.py) |
| `--source-genome-build` | `auto` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 14/3 | [sumstats_munger.py:1015](../../../src/ldsc/sumstats_munger.py) |
| `--output-genome-build` | `None` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 15/1 | [sumstats_munger.py:1021](../../../src/ldsc/sumstats_munger.py) |
| `--liftover-chain-file` | `None` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 9/2 | [sumstats_munger.py:1027](../../../src/ldsc/sumstats_munger.py) |
| `--use-hm3-quick-liftover` | `False` | Keep method; naming candidate --liftover-hm3; depends on HM3 restriction. | 9/1 | [sumstats_munger.py:1032](../../../src/ldsc/sumstats_munger.py) |
| `--trait-name` | `None` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 11/1 | [sumstats_munger.py:1038](../../../src/ldsc/sumstats_munger.py) |
| `--log-level` | `INFO` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 7/2 | [sumstats_munger.py:1039](../../../src/ldsc/sumstats_munger.py) |
| `--output-format` | `parquet` | Keep export ability; candidate: always Parquet plus optional legacy export. | 8/1 | [sumstats_munger.py:1040](../../../src/ldsc/sumstats_munger.py) |
| `--N` | `None` | Keep fallback N; fix incorrect override help; rename --sample-size. | 1/1 | [sumstats_munger.py:1046](../../../src/ldsc/sumstats_munger.py) |
| `--N-cas` | `None` | Consolidation candidate: fixed counts are summed; keep provenance deliberately. | 1/1 | [sumstats_munger.py:1050](../../../src/ldsc/sumstats_munger.py) |
| `--N-con` | `None` | Consolidation candidate: fixed counts are summed; keep per-variant count columns. | 1/1 | [sumstats_munger.py:1054](../../../src/ldsc/sumstats_munger.py) |
| `--info-min` | `0.9` | Keep; only applies when INFO is present; explicit missing-field request should fail. | 1/0 | [sumstats_munger.py:1058](../../../src/ldsc/sumstats_munger.py) |
| `--maf-min` | `0.01` | Keep; only applies when FRQ is present; explicit missing-field request should fail. | 6/3 | [sumstats_munger.py:1060](../../../src/ldsc/sumstats_munger.py) |
| `--n-min` | `None` | Keep; fix documented default and zero/constant-N behavior. | 1/0 | [sumstats_munger.py:1062](../../../src/ldsc/sumstats_munger.py) |
| `--chunksize` | `1000000` | Keep advanced memory control; rename --chunk-size. | 1/1 | [sumstats_munger.py:1064](../../../src/ldsc/sumstats_munger.py) |
| `--snp` | `None` | Keep explicit schema override; rename --snp-col. | 5/1 | [sumstats_munger.py:1066](../../../src/ldsc/sumstats_munger.py) |
| `--chr` | `None` | Keep explicit schema override; rename --chr-col. | 6/1 | [sumstats_munger.py:1068](../../../src/ldsc/sumstats_munger.py) |
| `--pos` | `None` | Keep explicit schema override; rename --pos-col. | 6/1 | [sumstats_munger.py:1070](../../../src/ldsc/sumstats_munger.py) |
| `--N-col` | `None` | Keep; rename --n-col; per-variant strategy is scientifically distinct. | 7/1 | [sumstats_munger.py:1072](../../../src/ldsc/sumstats_munger.py) |
| `--N-cas-col` | `None` | Keep; rename --case-count-col; per-variant strategy is scientifically distinct. | 4/1 | [sumstats_munger.py:1075](../../../src/ldsc/sumstats_munger.py) |
| `--N-con-col` | `None` | Keep; rename --control-count-col; per-variant strategy is scientifically distinct. | 4/1 | [sumstats_munger.py:1078](../../../src/ldsc/sumstats_munger.py) |
| `--a1` | `None` | Keep explicit schema override; rename --a1-col. | 4/0 | [sumstats_munger.py:1081](../../../src/ldsc/sumstats_munger.py) |
| `--a2` | `None` | Keep explicit schema override; rename --a2-col. | 4/0 | [sumstats_munger.py:1083](../../../src/ldsc/sumstats_munger.py) |
| `--p` | `None` | Keep explicit schema override; rename --p-col. | 2/0 | [sumstats_munger.py:1085](../../../src/ldsc/sumstats_munger.py) |
| `--frq` | `None` | Keep explicit schema override; rename --frequency-col. | 1/0 | [sumstats_munger.py:1087](../../../src/ldsc/sumstats_munger.py) |
| `--signed-sumstats` | `None` | Rename --signed-stat COLUMN,NULL; supplies sign, P supplies magnitude. | 2/1 | [sumstats_munger.py:1089](../../../src/ldsc/sumstats_munger.py) |
| `--info` | `None` | Keep explicit schema override; rename --info-col. | 1/0 | [sumstats_munger.py:1091](../../../src/ldsc/sumstats_munger.py) |
| `--info-list` | `None` | Repair multi-column rejection; distinguish list-valued cells from multiple INFO columns. | 3/1 | [sumstats_munger.py:1093](../../../src/ldsc/sumstats_munger.py) |
| `--nstudy` | `None` | Keep explicit schema override; rename --study-count-col. | 1/0 | [sumstats_munger.py:1095](../../../src/ldsc/sumstats_munger.py) |
| `--nstudy-min` | `None` | Advanced conditional filter; ignored with N columns; zero selects default. | 1/0 | [sumstats_munger.py:1097](../../../src/ldsc/sumstats_munger.py) |
| `--ignore` | `None` | Keep; rename --ignore-cols; overrides inference/hints. | 5/1 | [sumstats_munger.py:1100](../../../src/ldsc/sumstats_munger.py) |
| `--a1-inc` | `False` | Advanced assumption; rename --a1-is-increasing; ignores signed-statistic direction. | 1/0 | [sumstats_munger.py:1102](../../../src/ldsc/sumstats_munger.py) |
| `--keep-maf` | `False` | Rename --keep-frequency, or always preserve supplied FRQ; output is not necessarily MAF. | 1/0 | [sumstats_munger.py:1104](../../../src/ldsc/sumstats_munger.py) |
| `--snp-identifier` | `chr_pos_allele_aware` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 26/10 | [sumstats_munger.py:1106](../../../src/ldsc/sumstats_munger.py) |

## h2 (15)

| Flag | Default / required | Disposition | D/T | Declaration |
| --- | --- | --- | --- | --- |
| `--ldscore-dir` | `required` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 16/3 | [regression_runner.py:2912](../../../src/ldsc/regression_runner.py) |
| `--count-kind` | `common` | Clarify/rename: --reference-snp-count-kind; common/all are distinct. | 9/1 | [regression_runner.py:2913](../../../src/ldsc/regression_runner.py) |
| `--output-dir` | `required` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 37/18 | [regression_runner.py:2919](../../../src/ldsc/regression_runner.py) |
| `--overwrite` | `False` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 36/10 | [regression_runner.py:2920](../../../src/ldsc/regression_runner.py) |
| `--n-blocks` | `200` | Clarify/rename: --jackknife-blocks. | 5/1 | [regression_runner.py:2921](../../../src/ldsc/regression_runner.py) |
| `--no-intercept` | `False` | Rename: --fix-intercept-defaults (h2=1, covariance=0). | 5/0 | [regression_runner.py:2922](../../../src/ldsc/regression_runner.py) |
| `--allow-identity-downgrade` | `False` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 8/0 | [regression_runner.py:2923](../../../src/ldsc/regression_runner.py) |
| `--intercept-h2` | `None` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 3/1 | [regression_runner.py:2615](../../../src/ldsc/regression_runner.py) |
| `--two-step-cutoff` | `None` | Keep advanced estimator control; fixed-intercept/multi-annotation combinations are invalid. | 3/0 | [regression_runner.py:2935](../../../src/ldsc/regression_runner.py) |
| `--chisq-max` | `None` | Keep; document h2 Z^2 bound versus rg |Z1*Z2| bound. | 4/1 | [regression_runner.py:2936](../../../src/ldsc/regression_runner.py) |
| `--log-level` | `INFO` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 7/2 | [regression_runner.py:2937](../../../src/ldsc/regression_runner.py) |
| `--sumstats-file` | `required` | Keep; help is stale: canonical Parquet is accepted and should lead. | 12/4 | [regression_runner.py:2567](../../../src/ldsc/regression_runner.py) |
| `--trait-name` | `None` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 11/1 | [regression_runner.py:2568](../../../src/ldsc/regression_runner.py) |
| `--samp-prev` | `None` | Rename: --sample-prevalence. | 10/1 | [regression_runner.py:2548](../../../src/ldsc/regression_runner.py) |
| `--pop-prev` | `None` | Rename: --population-prevalence. | 10/1 | [regression_runner.py:2555](../../../src/ldsc/regression_runner.py) |

## partitioned-h2 (17)

| Flag | Default / required | Disposition | D/T | Declaration |
| --- | --- | --- | --- | --- |
| `--query-batch-size` | `1000` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 9/1 | [regression_runner.py:2574](../../../src/ldsc/regression_runner.py) |
| `--ldscore-dir` | `required` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 16/3 | [regression_runner.py:2912](../../../src/ldsc/regression_runner.py) |
| `--count-kind` | `common` | Clarify/rename: --reference-snp-count-kind; common/all are distinct. | 9/1 | [regression_runner.py:2913](../../../src/ldsc/regression_runner.py) |
| `--output-dir` | `required` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 37/18 | [regression_runner.py:2919](../../../src/ldsc/regression_runner.py) |
| `--overwrite` | `False` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 36/10 | [regression_runner.py:2920](../../../src/ldsc/regression_runner.py) |
| `--n-blocks` | `200` | Clarify/rename: --jackknife-blocks. | 5/1 | [regression_runner.py:2921](../../../src/ldsc/regression_runner.py) |
| `--no-intercept` | `False` | Rename: --fix-intercept-defaults (h2=1, covariance=0). | 5/0 | [regression_runner.py:2922](../../../src/ldsc/regression_runner.py) |
| `--allow-identity-downgrade` | `False` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 8/0 | [regression_runner.py:2923](../../../src/ldsc/regression_runner.py) |
| `--intercept-h2` | `None` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 3/1 | [regression_runner.py:2615](../../../src/ldsc/regression_runner.py) |
| `--two-step-cutoff` | `None` | Keep advanced estimator control; fixed-intercept/multi-annotation combinations are invalid. | 3/0 | [regression_runner.py:2935](../../../src/ldsc/regression_runner.py) |
| `--chisq-max` | `None` | Keep; document h2 Z^2 bound versus rg |Z1*Z2| bound. | 4/1 | [regression_runner.py:2936](../../../src/ldsc/regression_runner.py) |
| `--log-level` | `INFO` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 7/2 | [regression_runner.py:2937](../../../src/ldsc/regression_runner.py) |
| `--sumstats-file` | `required` | Keep; help is stale: canonical Parquet is accepted and should lead. | 12/4 | [regression_runner.py:2567](../../../src/ldsc/regression_runner.py) |
| `--trait-name` | `None` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 11/1 | [regression_runner.py:2568](../../../src/ldsc/regression_runner.py) |
| `--samp-prev` | `None` | Rename: --sample-prevalence. | 10/1 | [regression_runner.py:2548](../../../src/ldsc/regression_runner.py) |
| `--pop-prev` | `None` | Rename: --population-prevalence. | 10/1 | [regression_runner.py:2555](../../../src/ldsc/regression_runner.py) |
| `--summary-sort-by` | `auto` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 5/2 | [regression_runner.py:2580](../../../src/ldsc/regression_runner.py) |

## quantile-h2 (18)

| Flag | Default / required | Disposition | D/T | Declaration |
| --- | --- | --- | --- | --- |
| `--partitioned-h2-result-dir` | `required` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 3/0 | [quantile_h2.py:446](../../../src/ldsc/quantile_h2.py) |
| `--baseline-annot-sources` | `required` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 19/8 | [quantile_h2.py:447](../../../src/ldsc/quantile_h2.py) |
| `--query-annot-sources` | `None` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 7/1 | [quantile_h2.py:449](../../../src/ldsc/quantile_h2.py) |
| `--query-annot-bed-sources` | `None` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 12/3 | [quantile_h2.py:450](../../../src/ldsc/quantile_h2.py) |
| `--query-annot-gene-list-sources` | `None` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 12/5 | [quantile_h2.py:451](../../../src/ldsc/quantile_h2.py) |
| `--gene-coordinate-file` | `None` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 11/5 | [quantile_h2.py:452](../../../src/ldsc/quantile_h2.py) |
| `--control-gene-list-file` | `None` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 11/3 | [quantile_h2.py:453](../../../src/ldsc/quantile_h2.py) |
| `--gene-list-resolution-policy` | `strict` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 8/2 | [quantile_h2.py:454](../../../src/ldsc/quantile_h2.py) |
| `--gene-exclude-regions` | `none` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 9/3 | [quantile_h2.py:455](../../../src/ldsc/quantile_h2.py) |
| `--padding-bp` | `0` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 17/5 | [quantile_h2.py:456](../../../src/ldsc/quantile_h2.py) |
| `--target-annot-sources` | `required` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 3/0 | [quantile_h2.py:457](../../../src/ldsc/quantile_h2.py) |
| `--target-annotation` | `required` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 3/0 | [quantile_h2.py:458](../../../src/ldsc/quantile_h2.py) |
| `--ref-metadata-sources` | `required` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 4/0 | [quantile_h2.py:459](../../../src/ldsc/quantile_h2.py) |
| `--target-missing-value` | `None` | Clarify/rename: --target-missing-token; exclusion sentinel, not imputation. | 3/0 | [quantile_h2.py:460](../../../src/ldsc/quantile_h2.py) |
| `--num-quantiles` | `5` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 2/0 | [quantile_h2.py:461](../../../src/ldsc/quantile_h2.py) |
| `--output-dir` | `required` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 37/18 | [quantile_h2.py:462](../../../src/ldsc/quantile_h2.py) |
| `--overwrite` | `False` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 36/10 | [quantile_h2.py:463](../../../src/ldsc/quantile_h2.py) |
| `--log-level` | `INFO` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 7/2 | [quantile_h2.py:464](../../../src/ldsc/quantile_h2.py) |

## rg (18)

| Flag | Default / required | Disposition | D/T | Declaration |
| --- | --- | --- | --- | --- |
| `--ldscore-dir` | `required` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 16/3 | [regression_runner.py:2912](../../../src/ldsc/regression_runner.py) |
| `--count-kind` | `common` | Clarify/rename: --reference-snp-count-kind; common/all are distinct. | 9/1 | [regression_runner.py:2913](../../../src/ldsc/regression_runner.py) |
| `--output-dir` | `required` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 37/18 | [regression_runner.py:2919](../../../src/ldsc/regression_runner.py) |
| `--overwrite` | `False` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 36/10 | [regression_runner.py:2920](../../../src/ldsc/regression_runner.py) |
| `--n-blocks` | `200` | Clarify/rename: --jackknife-blocks. | 5/1 | [regression_runner.py:2921](../../../src/ldsc/regression_runner.py) |
| `--no-intercept` | `False` | Rename: --fix-intercept-defaults (h2=1, covariance=0). | 5/0 | [regression_runner.py:2922](../../../src/ldsc/regression_runner.py) |
| `--allow-identity-downgrade` | `False` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 8/0 | [regression_runner.py:2923](../../../src/ldsc/regression_runner.py) |
| `--two-step-cutoff` | `None` | Keep advanced estimator control; fixed-intercept/multi-annotation combinations are invalid. | 3/0 | [regression_runner.py:2935](../../../src/ldsc/regression_runner.py) |
| `--chisq-max` | `None` | Keep; document h2 Z^2 bound versus rg |Z1*Z2| bound. | 4/1 | [regression_runner.py:2936](../../../src/ldsc/regression_runner.py) |
| `--log-level` | `INFO` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 7/2 | [regression_runner.py:2937](../../../src/ldsc/regression_runner.py) |
| `--sumstats-sources` | `required` | Keep; help is stale: canonical Parquet is accepted and should lead. | 8/2 | [regression_runner.py:2595](../../../src/ldsc/regression_runner.py) |
| `--anchor-trait` | `None` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 8/1 | [regression_runner.py:2604](../../../src/ldsc/regression_runner.py) |
| `--write-per-pair-detail` | `False` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 4/0 | [regression_runner.py:2609](../../../src/ldsc/regression_runner.py) |
| `--intercept-h2` | `None` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 3/1 | [regression_runner.py:2615](../../../src/ldsc/regression_runner.py) |
| `--intercept-gencov` | `None` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 2/1 | [regression_runner.py:2616](../../../src/ldsc/regression_runner.py) |
| `--samp-prev` | `None` | Rename: --sample-prevalence. | 10/1 | [regression_runner.py:2548](../../../src/ldsc/regression_runner.py) |
| `--pop-prev` | `None` | Rename: --population-prevalence. | 10/1 | [regression_runner.py:2555](../../../src/ldsc/regression_runner.py) |
| `--prevalence-manifest` | `None` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 4/0 | [regression_runner.py:2631](../../../src/ldsc/regression_runner.py) |

## query-r2 (7)

| Flag | Default / required | Disposition | D/T | Declaration |
| --- | --- | --- | --- | --- |
| `--panel-dir` | `required` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 2/1 | [r2_query.py:510](../../../src/ldsc/r2_query.py) |
| `--pairs` | `required` | Rename: --pairs-file; retain stdin support. | 2/1 | [r2_query.py:511](../../../src/ldsc/r2_query.py) |
| `--output-dir` | `required` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 37/18 | [r2_query.py:512](../../../src/ldsc/r2_query.py) |
| `--overwrite` | `False` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 36/10 | [r2_query.py:513](../../../src/ldsc/r2_query.py) |
| `--snp-identifier` | `None` | Keep query rekeying capability; clarify that override selects query matching mode. | 26/10 | [r2_query.py:514](../../../src/ldsc/r2_query.py) |
| `--genome-build` | `None` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 24/12 | [r2_query.py:515](../../../src/ldsc/r2_query.py) |
| `--log-level` | `INFO` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 7/2 | [r2_query.py:516](../../../src/ldsc/r2_query.py) |

## convert-h2-scale (7)

| Flag | Default / required | Disposition | D/T | Declaration |
| --- | --- | --- | --- | --- |
| `--h2-result-dir` | `required` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 3/2 | [h2_scale.py:257](../../../src/ldsc/h2_scale.py) |
| `--samp-prev` | `required` | Rename: --sample-prevalence. | 10/1 | [h2_scale.py:258](../../../src/ldsc/h2_scale.py) |
| `--pop-prev` | `None` | Rename: --population-prevalence. | 10/1 | [h2_scale.py:260](../../../src/ldsc/h2_scale.py) |
| `--pop-prev-range` | `None` | Rename: --population-prevalence-range. | 3/0 | [h2_scale.py:261](../../../src/ldsc/h2_scale.py) |
| `--num-points` | `201` | Keep range mode only; ignored for exact --pop-prev. | 2/0 | [h2_scale.py:268](../../../src/ldsc/h2_scale.py) |
| `--overwrite` | `False` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 36/10 | [h2_scale.py:269](../../../src/ldsc/h2_scale.py) |
| `--log-level` | `INFO` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 7/2 | [h2_scale.py:270](../../../src/ldsc/h2_scale.py) |

## plot (3)

| Flag | Default / required | Disposition | D/T | Declaration |
| --- | --- | --- | --- | --- |
| `--result-dir` | `required` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 6/2 | [__init__.py:162](../../../src/ldsc/plotting/__init__.py) |
| `--overwrite` | `False` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 36/10 | [__init__.py:163](../../../src/ldsc/plotting/__init__.py) |
| `--log-level` | `INFO` | Keep: active workflow/configuration/artifact control; see report for mode restrictions. | 7/2 | [__init__.py:164](../../../src/ldsc/plotting/__init__.py) |
