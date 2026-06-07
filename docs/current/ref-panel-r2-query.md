# Reference-Panel R² Pair Query

Query the adjusted (unbiased) R² stored in a `ldsc build-ref-panel` index-format
panel for a list of SNP pairs, and convert it to a signed Pearson correlation
`r`. This is the random-access counterpart to the LD-score read path: the
LD-score workflow streams every stored pair once, while this tool looks up
specific pairs by SNP identity.

See the design spec
`docs/superpowers/specs/2026-06-06-ref-panel-r2-query-design.md` for the full
semantics, and `docs/current/parquet-r2-format-and-read-pipeline.md` for the
panel format.

## Python API

Three public objects are exported from `ldsc`:

- `R2Panel.open(...)` — a reusable handle. Open once, query many times; each
  chromosome is loaded lazily on first touch and cached (sidecar, schema
  metadata, sidecar-binding check, identity-key map, parquet handle).
- `query_r2(pairs, ...)` — a one-shot wrapper that opens a panel, queries, and
  returns the result.
- `unbiased_r2_to_pearson_r(r2_adj, n, sign=None)` — a pure converter.

```python
import pandas as pd
from ldsc import query_r2, unbiased_r2_to_pearson_r

pairs = pd.DataFrame({
    "CHR_1": [1, 1], "POS_1": [752721, 776546], "A1_1": ["A", "G"], "A2_1": ["G", "A"],
    "CHR_2": [1, 1], "POS_2": [777122, 798959], "A1_2": ["A", "T"], "A2_2": ["C", "G"],
})
result = query_r2(pairs, panel_dir="ref_panel", genome_build="hg38", with_r=True)
print(result[["r2", "sign", "status", "r"]])
```

For repeated queries against the same panel, keep the handle:

```python
from ldsc import R2Panel

panel = R2Panel.open("ref_panel", genome_build="hg38")
out1 = panel.query_pairs(pairs_a)
out2 = panel.query_pairs(pairs_b, with_r=True)
print(panel.chromosomes, panel.snp_identifier, panel.n_samples)
```

A single chromosome can be opened by explicit paths instead of a directory:

```python
panel = R2Panel.open(meta_path="chr1_meta.tsv.gz", parquet_path="chr1_r2.parquet")
```

Supplying neither `panel_dir` nor an explicit `meta_path`+`parquet_path` pair —
or supplying both — is a usage error.

## Pair input format

`pairs` is a `DataFrame` (the CLI reads a TSV/CSV into one). Each row is one
queried pair; the two endpoints are columns suffixed `_1` and `_2`:

| Endpoint 1 | Endpoint 2 | Required when |
| --- | --- | --- |
| `CHR_1, POS_1` | `CHR_2, POS_2` | `chr_pos*` modes |
| `A1_1, A2_1` | `A1_2, A2_2` | allele-aware modes (matching + sign) |
| `SNP_1` | `SNP_2` | `rsid*` modes |

The SNP identifier mode defaults to the panel's recorded `ldsc:snp_identifier`
and can be overridden with `snp_identifier=` / `--snp-identifier`.

**Base modes ignore alleles.** In `rsid` / `chr_pos` mode only the matching
columns are read; any `A1_*`/`A2_*` columns are dropped before resolution and
never used. A base-mode query behaves identically with or without alleles
supplied. Allele columns are read only in allele-aware modes.

## Output schema

The result echoes the input columns (in input row order) plus:

| Column | Meaning |
| --- | --- |
| `r2` | adjusted (unbiased) R² as stored, dequantized to float; `NaN` when not available |
| `sign` | harmonized sign `+1`/`-1` (allele-aware modes); `NA` in base modes or where `r2` is `NaN` |
| `status` | empty (`""`) for a valid `r2`; otherwise the cause of the `NaN` |
| `r` | signed Pearson `r` (only when `with_r=True` / `--with-r`); `NaN` where `r2` is `NaN` or `sign` is `NA` |

`status` is always present and carries one of a closed vocabulary, only for
`NaN`-`r2` rows:

- `not_in_panel` — one or both endpoints did not resolve to a panel SNP under the
  active mode (including invalid/ambiguous query alleles in allele-aware mode).
- `cross_chromosome` — both endpoints resolved, but to different chromosomes.
- `absent` — both endpoints are on the same chromosome but the pair is not stored
  (outside the panel's LD window or below its `min_r2` floor).

A pair whose two endpoints resolve to the same panel SNP (the diagonal) returns
`r2 = 1.0` with empty `status`.

## Sign convention

The stored sign is in the panel's A1/A2 orientation. In **allele-aware** modes
the returned `sign` is harmonized to *your* allele coding: for each endpoint the
query alleles are classified as aligned or swapped against the panel (allowing
strand complement), and the pair sign is flipped when an odd number of endpoints
are swapped. In **base** modes no allele information is consulted, so `sign` is
always `NA`. Because base modes carry no sign, `--with-r` there yields an
all-`NaN` `r` column (with a warning); use an allele-aware panel/mode for signed
`r`.

## R²→r conversion

`unbiased_r2_to_pearson_r` inverts the unbiased correction
`r2_adj = r2_raw - (1 - r2_raw)/(n - 2)` to recover the biased squared
correlation, then takes the root and applies the sign. It is vectorized and pure;
`sign=None` returns the magnitude `|r|`. The handle/`with_r` path supplies the
panel's `ldsc:n_samples` and the harmonized `sign` automatically.

## Lookup strategy

Pairs are matched by an int64 key against the parquet's row groups. Two
strategies share that primitive:

- **random-access** — prune row groups using `IDX_1` footer statistics; best for
  small/interactive queries.
- **streaming** — scan every row group once; best for large batches.

`strategy="auto"` (the default) chooses random-access when the query has at most
`strategy_threshold` pairs (default 50 000), otherwise streaming. Both are
overridable via `strategy=` / `--strategy` and `strategy_threshold=` /
`--strategy-threshold`; they produce identical results.

## CLI

```
ldsc query-r2 (--panel-dir DIR | --meta M --parquet P)
              --pairs FILE [--out FILE]
              [--snp-identifier MODE] [--genome-build {hg19,hg38}]
              [--with-r]
              [--strategy {auto,random,stream}] [--strategy-threshold N]
```

`--pairs` is a TSV/CSV with the `_1`/`_2` endpoint columns (`-` reads stdin). The
output TSV echoes the input columns plus `r2`, `sign`, `status` (and `r` with
`--with-r`), written to `--out` or stdout.

```bash
ldsc query-r2 --panel-dir ref_panel --genome-build hg38 \
    --pairs pairs.tsv --out pairs_r2.tsv --with-r
```
