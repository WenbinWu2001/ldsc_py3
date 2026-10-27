# Resource navigation

Last updated on: 2026-10-27

Select resources that match the reference population, genome build, and intended annotation universe. The example `/path/to/...` and project-specific paths in tutorials must be replaced with your own inputs.

| Resource | Preparation and contract |
| --- | --- |
| Packaged HapMap3 maps and exclusion regions | See [packaged resource notes](../../src/ldsc/data/readme.txt), [attribution](../../src/ldsc/data/ATTRIBUTION.txt), and [region presets](../reference/formats/region-exclusion-presets.md). |
| PLINK or R² reference panel | Supply a complete BED/BIM/FAM genotype trio or [build reusable R² tables](../wiki/main-functionalities/build-r2-panel.md). Keep each R² Parquet with its matching metadata sidecar. |
| Baseline annotations | Supply the actual annotation files for direct query LD scoring or index construction. See [LD scoring](workflows/ld-scores.md). |
| Legacy LD-score suites | Use [explicit conversion](workflows/migrating-from-ldsc2.md); precomputed scores are a different input from baseline annotation files. |
| Query BEDs and gene lists | See [BED preparation](workflows/annotations.md) and [gene-list inputs](workflows/gene-list-ld-scores.md#gene-lists-and-coordinate-authority). Direct gene lists require a caller-supplied coordinate catalog; no catalog fallback is packaged. |
| Small local examples | [Deterministic test resources](../../tests/fixtures/minimal_external_resources/README.md) support smoke tests and demonstrations, not production inference. |

For exact paths, glob selection, and chromosome-suite requirements, use the [path specification](../reference/input-paths.md#pattern-support-in-command-help).
