# Validation

Last updated on: 2026-10-27

This page records what has been validated in LDSC3, with which evidence, and what has not.
Each evidence page states its date, revision and limits; raw logs and probe scripts are kept outside the repository.

## Recorded evidence

| Evidence | Page | Status |
| --- | --- | --- |
| Coordinates of the bundled HM3 map against UCSC, Ensembl and NCBI | [HM3 map coordinates](hm3-map-coordinates.md) | Describes the map as of 2026-08-05; the map was rebuilt from dbSNP build 157 on 2026-10-27 and has not been through the full protocol again. |
| Parallel execution, regression and LD-score resource measurements | [Performance](performance.md) | Local and Longleaf pilots with stated limits; no measurement is a guarantee for other data or systems. |
| Scientific equivalence across implementation changes | The test suite (`tests/`), including fixed golden NPZ files under `tests/fixtures/golden/` | Covers the package's own consistency; see the test map in [code structure](../code-structure.md). |

## Not yet validated

- **Numerical equivalence with LDSC2.** No document records a comparison of LDSC3 LD scores, SNP heritability, genetic correlation or partitioned heritability with LDSC2 on the same inputs; [migrating from LDSC2](../../user-guide/workflows/migrating-from-ldsc2.md) describes the intended differences.
- **Region-exclusion presets.** The bundled coordinates have not yet been independently checked against their sources ([region-exclusion presets](../../reference/formats/region-exclusion-presets.md)).
- **Production-scale resources.** Genome-wide memory needs, LD-score scaling across chromosomes and input and output performance on HPC storage have not been measured ([performance](performance.md#limits)).

Add new evidence here as a page with its date, revision and limits, and summarize its result in the table above.
