# R² Matrix Pair Query

Last updated on: 2026-10-28

Query the unbiased R² stored in an R² matrix built by `ldsc build-r2-matrix` for a list of SNP pairs, and convert it to a signed Pearson correlation `r`.
This is the random-access counterpart to the LD-score read path: the LD-score workflow streams every stored pair once, while this tool looks up specific pairs by SNP identity.

See the [R² matrix format](formats/r2-matrices.md) and the [query wiki guide](../user-guide/workflows/r2-matrices.md) for a worked command example.

## Python API

Three public objects are exported from `ldsc`:

- `R2Matrix.open(...)` — a reusable handle. Open once, query many times.
  `open()` validates every chromosome's parquet footer up front; each chromosome's sidecar, sidecar-binding check, identity-key map, and parquet handle are loaded lazily on first touch and cached.
- `query_r2_matrix(pairs, ...)` — a one-shot wrapper that opens an R² matrix, queries, and
  returns the result.
- `unbiased_r2_to_pearson_r(r2_adj, n, sign=None)` — a pure converter.

```python
import pandas as pd
from ldsc import query_r2_matrix, unbiased_r2_to_pearson_r

pairs = pd.DataFrame({
    "CHR_1": [1, 1], "POS_1": [752721, 776546], "A1_1": ["A", "G"], "A2_1": ["G", "A"],
    "CHR_2": [1, 1], "POS_2": [777122, 798959], "A1_2": ["A", "T"], "A2_2": ["C", "G"],
})
result = query_r2_matrix(pairs, r2_matrix_dir="r2_matrix", genome_build="hg38")
print(result[["r2", "sign_r", "r", "status"]])
```

The signed Pearson `r` column is always included (it is all-NaN in base modes).

For repeated queries against the same R² matrix, keep the handle:

```python
from ldsc import R2Matrix

matrix = R2Matrix.open("r2_matrix", genome_build="hg38")
out1 = matrix.query_pairs(pairs_a)
out2 = matrix.query_pairs(pairs_b)
print(matrix.chromosomes, matrix.snp_identifier, matrix.n_samples)
```

A `build-r2-matrix` output directory is the only R² matrix input mode; omitting
`r2_matrix_dir` is a usage error.

## Pair input format

`pairs` is a `DataFrame` (the CLI reads a TSV/CSV into one). Each row is one
queried pair; the two endpoints are columns suffixed `_1` and `_2`:

| Endpoint 1 | Endpoint 2 | Required when |
| --- | --- | --- |
| `CHR_1, POS_1` | `CHR_2, POS_2` | `chr_pos*` modes |
| `A1_1, A2_1` | `A1_2, A2_2` | allele-aware modes (matching + sign) |
| `SNP_1` | `SNP_2` | `rsid*` modes |

The SNP identifier mode defaults to the R² matrix's recorded `ldsc:snp_identifier`.
`snp_identifier=` / `--snp-identifier` accepts only that mode or its base form (`rsid` or `chr_pos`); any other mode fails at `R2Matrix.open()`.

**Base modes ignore alleles.** In `rsid` / `chr_pos` mode only the matching
columns are read; any `A1_*`/`A2_*` columns are dropped before resolution and
never used. A base-mode query behaves identically with or without alleles
supplied. Allele columns are read only in allele-aware modes.

## Output schema

The result preserves input rows and columns and writes `r2`, `sign_r`, `r`, and `status` in that order for fresh input. Existing columns with those names are replaced. No `sign` compatibility column is emitted. The persisted `query_r2_matrix.tsv` uses the same names and writes missing values as `NaN`.

| Column | Meaning |
| --- | --- |
| `r2` | unbiased R² as stored, dequantized to float; `NaN` when not available |
| `sign_r` | harmonized sign `+1`/`-1` (allele-aware modes); `NA` in base modes or where `r2` is `NaN` |
| `status` | empty (`""`) for a valid `r2`; otherwise the cause of the `NaN` |
| `r` | signed Pearson `r` (always emitted); `NaN` where `r2` is `NaN` or `sign_r` is `NA` |

`status` is always present and carries one of a closed vocabulary, only for
`NaN`-`r2` rows:

- `not_in_matrix` — one or both endpoints did not resolve to an R² matrix SNP under the
  active mode (including invalid/ambiguous query alleles in allele-aware mode).
- `cross_chromosome` — both endpoints resolved, but to different chromosomes.
- `absent` — both endpoints are on the same chromosome but the pair is not stored
  (outside the R² matrix's LD window or below its `r2_min` floor).

A pair whose two endpoints resolve to the same R² matrix SNP (the diagonal) returns
`r2 = 1.0` with empty `status`.

## Sign convention

The stored `SIGN_R` is boolean (`True` means Pearson r ≥ 0) in the R² matrix's A1/A2 orientation. The output `sign_r` is nullable Int8 (+1/-1), not that stored boolean. In **allele-aware** modes
the returned `sign_r` is harmonized to *your* allele coding: for each endpoint the
query alleles are classified as aligned or swapped against the R² matrix (allowing
strand complement), and the pair sign is flipped when an odd number of endpoints
are swapped. In **base** modes no allele information is consulted, so `sign_r` is
always `NA`. Because base modes carry no sign, the always-emitted `r` column is
all-`NaN` there (with a warning); use an allele-aware R² matrix/mode for signed `r`.

## R²→r conversion

`unbiased_r2_to_pearson_r` inverts the unbiased correction
`r2_adj = r2_raw - (1 - r2_raw)/(n - 2)` to recover the biased squared
correlation, then takes the root and applies the sign. It is vectorized and pure;
`sign=None` returns the magnitude `|r|`. The query path supplies the R² matrix's
`ldsc:n_samples` and the harmonized `sign_r` automatically.

## Lookup strategy

Pairs are matched by an int64 key against the row groups of the R² Parquet file.
The lookup auto-selects between two internal paths by query size (not user-configurable):

- **random-access** — prune row groups using `IDX_1` footer statistics; used for
  small/interactive queries (≤ 50 000 pairs).
- **streaming** — scan every row group once; used for large batches.

Both produce identical results.

## CLI

```
ldsc query-r2-matrix --r2-matrix-dir DIR
              --pairs-file FILE --output-dir DIR [--overwrite] [--log-level LEVEL]
              [--snp-identifier MODE] [--genome-build {hg19,hg38}]
```

`--pairs-file` is a TSV/CSV with the `_1`/`_2` endpoint columns (`-` reads stdin). The
result echoes the input columns plus `r2`, `sign_r`, `r`, and `status`.

`query-r2-matrix` requires `--output-dir` and follows the canonical result-directory
layout shared by the other commands: `query_r2_matrix.tsv` at the root, plus
`diagnostics/metadata.json` (provenance and per-run status counts) and
`diagnostics/query-r2-matrix.log` (workflow log). `--overwrite` replaces the command-owned result files and log; `--log-level` sets the diagnostics-log verbosity.

```bash
# Canonical result directory
ldsc query-r2-matrix --r2-matrix-dir r2_matrix --genome-build hg38 \
    --pairs-file pairs.tsv --output-dir pairs_r2
```
