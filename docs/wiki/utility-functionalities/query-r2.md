# Query Pairwise R² and Signed r

Last updated on: 2026-09-11

`ldsc query-r2` looks up specific SNP pairs in a [build-r2-panel](../main-functionalities/build-r2-panel.md) output directory. It returns bias-adjusted R² and, in allele-aware modes, Pearson r in your query allele orientation. This lets you inspect local LD without loading the full pair table.

## Input and minimal command

Provide a TSV whose endpoint columns have `_1` and `_2` suffixes. For `chr_pos_allele_aware`, a row looks like this (illustrative coordinates):

```text
CHR_1	POS_1	A1_1	A2_1	CHR_2	POS_2	A1_2	A2_2
22	100000	A	G	22	101000	C	T
```

Use coordinates from the selected panel build and the alleles whose dosage correlation you want. rsID modes use `SNP_1` and `SNP_2`; allele-aware rsID mode also requires both allele pairs. Registered column aliases are accepted.

```bash
ldsc query-r2 \
  --panel-dir results/r2_panel/hg38 \
  --genome-build hg38 \
  --snp-identifier chr_pos_allele_aware \
  --pairs pairs.tsv \
  --output-dir results/pair_query
```

The panel must contain current `SIGN_R` Parquet files and their matching `chrN_meta.tsv.gz` sidecars. Individual Parquet files are not accepted as `--panel-dir`.

## Results

`query_r2.tsv` preserves input rows and adds:

| Column | Meaning |
| --- | --- |
| `r2` | Dequantized, bias-adjusted R²; can be negative |
| `sign_r` | +1/-1 in query A1 dosage orientation, or missing |
| `r` | Signed Pearson correlation, recovered using reference sample size |
| `status` | Empty for available R²; otherwise `not_in_panel`, `cross_chromosome`, or `absent` |

The Python result uses nullable Int8 for `sign_r`; the TSV writes missing values as `NaN`. There is no `sign` output alias. Stored `SIGN_R` is a boolean in panel allele orientation, whereas returned `sign_r` is harmonized to query alleles. Swapping exactly one endpoint's allele order reverses `sign_r` and `r` without changing `r2`.

`absent` means the resolved pair was not stored, for example because of the build window or a positive `--min-r2` threshold. It does not establish zero correlation. A resolved diagonal has R²=1. Base `rsid`/`chr_pos` modes ignore alleles and return missing `sign_r` and `r`; signed r also needs the panel's sample-size metadata.

The command writes `diagnostics/metadata.json` with provenance and status counts, plus `diagnostics/query-r2.log`. `--overwrite` replaces owned outputs; `--log-level` changes log-file verbosity.

## Options and repeated queries

| Option | Meaning |
| --- | --- |
| `--panel-dir DIR` | Required panel output root or build child |
| `--pairs FILE` | Required TSV or `.csv`; `-` reads TSV from stdin |
| `--output-dir DIR` | Required result directory |
| `--snp-identifier MODE` | Defaults to panel provenance; override to select a public identity mode |
| `--genome-build hg19\|hg38` | Select the coordinate build when needed |
| `--overwrite` | Replace existing command-owned outputs; default off |
| `--log-level LEVEL` | DEBUG, INFO (default), WARNING, or ERROR |

Small queries prune Parquet row groups; larger batches use a streaming scan. Both return the same values. A reusable Python `R2Panel.open(...)` handle caches chromosome state; use `query_r2(...)` for one-off calls. See the [API and format reference](../../current/ref-panel-r2-query.md) and the [runnable tutorial](../../../tutorials/build-parquet-reference-panel-from-plink.md#query-pairs).

Source: [`R2Panel.query_pairs`, CLI arguments, and metadata](../../../src/ldsc/r2_query.py), and [`QueryR2DirectoryWriter`](../../../src/ldsc/outputs.py).
