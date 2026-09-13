# Resource navigation

Last updated on: 2026-09-13

Select resources that match the reference population, genome build, and intended annotation universe. The example `/path/to/...` and project-specific paths in tutorials must be replaced with your own inputs.

| Resource | Preparation and contract |
| --- | --- |
| Packaged HapMap3 maps and exclusion regions | See [packaged resource notes](../../../src/ldsc/data/readme.txt), [attribution](../../../src/ldsc/data/ATTRIBUTION.txt), and [region presets](../../current/region-exclusion-presets.md). |
| PLINK or R² reference panel | Supply a complete BED/BIM/FAM genotype trio or [build reusable R² tables](../main-functionalities/build-r2-panel.md). Keep each R² Parquet with its matching metadata sidecar. |
| Baseline annotations | Supply the actual annotation files for direct query LD scoring or index construction. See [LD scoring](../main-functionalities/ldscore.md). |
| Legacy LD-score suites | Use [explicit conversion](convert-ldsc2-ldscores.md); precomputed scores are a different input from baseline annotation files. |
| Query BEDs and gene lists | See [BED preparation](how-to-customize-your-bed-files.md) and [gene-list inputs](../main-functionalities/ldscore-from-gene-list.md#gene-lists-and-coordinate-authority). Direct gene lists require a caller-supplied coordinate catalog; no catalog fallback is packaged. |
| Small local examples | [Deterministic test resources](../../../tests/fixtures/minimal_external_resources/README.md) support smoke tests and demonstrations, not production inference. |

For exact paths, glob selection, and chromosome-suite requirements, use the [path specification](../../current/path-specification.md#pattern-support-in-command-help).
