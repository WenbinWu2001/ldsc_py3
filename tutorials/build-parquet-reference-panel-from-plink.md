# Build an R² Panel from PLINK

Last updated on: 2026-09-14

`ldsc build-r2-panel` computes pairwise, bias-adjusted R² from PLINK genotypes and writes reusable chromosome tables plus SNP metadata. Use these outputs with `ldsc ldscore` to aggregate LD scores, or with `ldsc query-r2` to inspect particular SNP pairs. The current command has no `build-ref-panel` alias.

## Run a small example

Run this from the installed package repository root. The small hg38 chromosome-22 fixture is included under `tests/fixtures/minimal_external_resources/`; it is for demonstration, not a population reference for scientific analysis.

```bash
ldsc build-r2-panel \
  --plink-prefix tests/fixtures/minimal_external_resources/plink/hm3_chr22_subset \
  --source-genome-build hg38 \
  --ld-wind-snps 10 \
  --output-dir tutorial_outputs/r2_panel_chr22
```

The same arguments work with `python -m ldsc build-r2-panel`. Use a new output directory, or explicitly add `--overwrite` to replace the command's existing artifacts. This example supplies no genetic map or liftover chain, so it writes only hg38 and records missing `CM` values in the SNP sidecar.

For your own data, provide the prefix of a complete `.bed`, `.bim`, `.fam` trio. A quoted chromosome stem, `@` pattern, or ordinary glob can resolve multiple trios; see the [path rules](../docs/current/path-specification.md). The builder processes the resolved chromosomes sequentially. Separate chromosome invocations can share an output root because their diagnostics are scoped by chromosome.

## Choose the SNP universe and LD window

The `--plink-prefix` argument also accepts a plain stem such as `reference/1000G.EUR.QC.`. Its shared resolver discovers selected BED/BIM/FAM trios and reads BIM chromosomes; numeric filename suffixes remain intact and `@` is optional. Missing trio members, malformed contents, and conflicting chromosome sources fail together, with an audit in `diagnostics/plink_input_issues.tsv`. See the [shared PLINK resolution contract](../docs/current/path-specification.md#plink-prefix-resolution).

By default, all otherwise eligible reference SNPs contribute. PLINK filtering removes unusable variants; explicit `--maf-min`, `--ref-panel-snps-file`, and `--keep-indivs-file` further restrict variants or individuals. In coordinate modes, duplicate-coordinate groups are dropped using the fixed `drop-all` policy. There is no duplicate-position-policy CLI switch.

Choose exactly one window: `--ld-wind-snps N`, `--ld-wind-kb KB`, or `--ld-wind-cm CM`. The SNP-count window counts positions in the retained SNP sequence. A cM window requires a genetic map for every emitted build; SNP-count and kb windows do not. Provided maps populate sidecar `CM`; the builder does not use the BIM cM placeholder as the authoritative genetic map.

`--min-r2` defaults to zero, which keeps every pair within the window, including negative adjusted values. A positive threshold omits smaller values. Missing pairs contribute zero in LD-score accumulation, while `query-r2` reports `r2=NaN` and `status=absent`. Do not interpret an absent query pair as measured zero LD.

The reference individuals determine the population whose LD is estimated. An IID restriction file contains one individual ID per row. SNP restrictions are headered identity tables in the PLINK source build; use `SNP` for rsID modes or `CHR`/`POS` for coordinate modes, with allele columns when available. Allele-free restrictions match by base identity even in allele-aware modes.

## Output files

The concrete chromosome-22 example writes:

```text
tutorial_outputs/r2_panel_chr22/
├── hg38/
│   ├── chr22_r2.parquet
│   └── chr22_meta.tsv.gz
└── diagnostics/
    ├── metadata.chr22.json
    ├── build-r2-panel.chr22.log
    └── dropped_snps/
        └── chr22_dropped.tsv.gz
```

A suite invocation uses `diagnostics/metadata.json` and `diagnostics/build-r2-panel.log`. Diagnostics record provenance and processing decisions. The dropped-SNP audit is header-only when it has no rows to report. `--log-level` controls the workflow log's verbosity; it does not enable ordinary terminal progress. Failed authorized overwrites leave a `RUN_FAILED` marker; see the [logging contract](../docs/current/workflow-logging.md).

### Pair table

Each row of `chr22_r2.parquet` is one unordered, off-diagonal SNP pair within the requested window that survives the optional R² threshold. `IDX_1 < IDX_2`; rows are ordered by nondecreasing `IDX_1`. Diagonal entries are implicit and have R²=1.

| Column | Arrow type | Meaning |
| --- | --- | --- |
| `IDX_1` | int32 | First SNP's zero-based row index in the matching sidecar |
| `IDX_2` | int32 | Second SNP's zero-based row index |
| `R2` | int16 | Bias-adjusted R², quantized with scale 32767 |
| `SIGN_R` | bool | True when Pearson r ≥ 0 in panel A1 dosage orientation |

For reference sample size \(N>2\), the adjusted estimate is \(r^2-(1-r^2)/(N-2)\), upper-clipped to 1. Negative adjusted values are valid and retained by default. Decode the stored integer by dividing by the footer's `ldsc:r2_scale` (32767 for package output). `SIGN_R` describes the sign of Pearson r, not the sign of adjusted R². It refers to the sidecar's alleles, where A1 is the minor allele.

### SNP sidecar

`chr22_meta.tsv.gz` is a gzip TSV with leading `# ldsc:*` provenance comments and columns `CHR, POS, SNP, A1, A2, CM, MAF`. Rows are ordered by the emitted build's position. Coordinates are one-based; MAF refers to retained reference individuals. `SNP` preserves the BIM label, which need not be an rsID.

The sidecar is mandatory: its row order defines the pair indices, and its identity hash binds it to the Parquet. Keep the files together. There is no reconstruction from pair endpoints when it is missing. See the [format specification](../docs/current/parquet-r2-format-and-read-pipeline.md) for storage and binding details.

## Inspect the output safely

Read one row group rather than loading a large pair table in full:

```python
import pandas as pd
import pyarrow.parquet as pq

root = "tutorial_outputs/r2_panel_chr22/hg38"
meta = pd.read_csv(f"{root}/chr22_meta.tsv.gz", sep="\t", comment="#")
parquet = pq.ParquetFile(f"{root}/chr22_r2.parquet")
pairs = parquet.read_row_group(0).to_pandas()
scale = float(parquet.schema_arrow.metadata[b"ldsc:r2_scale"])
pairs["r2_decoded"] = pairs["R2"].astype("float32") / scale
print(meta.head())
print(pairs.head())
```

A direct Pandas/PyArrow read returns quantized integers; Parquet itself does not apply the R² scale. Never reorder or edit the sidecar independently of the pair table.

## Query pairs

Create a query from the first two sidecar rows, using the sidecar allele orientation:

```python
from ldsc import query_r2

query = pd.DataFrame({
    f"{column}_{endpoint}": [meta.iloc[row][column]]
    for endpoint, row in ((1, 0), (2, 1))
    for column in ("CHR", "POS", "A1", "A2")
})
result = query_r2(query, panel_dir=root, genome_build="hg38")
print(result[["r2", "sign_r", "r", "status"]])
query.to_csv("tutorial_outputs/r2_panel_chr22/pairs.tsv", sep="\t", index=False)
```

Or query the saved TSV through the CLI:

```bash
ldsc query-r2 \
  --panel-dir tutorial_outputs/r2_panel_chr22/hg38 \
  --genome-build hg38 \
  --pairs tutorial_outputs/r2_panel_chr22/pairs.tsv \
  --output-dir tutorial_outputs/r2_query_chr22
```

The result is `query_r2.tsv`, with `r2`, nullable `sign_r`, `r`, and `status`. `sign_r` is +1/-1 in the query's allele orientation; it is not the boolean `SIGN_R` storage column. Base identity modes return missing `sign_r` and `r`. See the [query reference](../docs/current/ref-panel-r2-query.md) for all statuses and allele handling.

## Compute LD scores

For the small SNP-window example, use the same or a narrower SNP-count window:

```bash
ldsc ldscore \
  --r2-dir tutorial_outputs/r2_panel_chr22/hg38 \
  --genome-build hg38 \
  --snp-identifier chr_pos_allele_aware \
  --ld-wind-snps 10 \
  --output-dir tutorial_outputs/ldscores_chr22
```

For scientific analyses, build the intended genome-wide reference population and window. Ordinary unpartitioned `ldscore` can synthesize an all-ones baseline; partitioned/query LD scores require explicit baseline annotations. A downstream cM window requires usable sidecar CM values. A wider window cannot recover pairs omitted during panel construction.

## Build from Python

The Python builder names remain `ReferencePanelBuilder`, `ReferencePanelBuildConfig`, and `run_build_ref_panel`; the CLI rename does not rename those APIs.

```python
from ldsc import GlobalConfig, run_build_ref_panel, set_global_config

set_global_config(GlobalConfig(snp_identifier="chr_pos_allele_aware", log_level="INFO"))
result = run_build_ref_panel(
    plink_prefix="tests/fixtures/minimal_external_resources/plink/hm3_chr22_subset",
    source_genome_build="hg38",
    ld_wind_snps=10,
    output_dir="tutorial_outputs/r2_panel_python",
)
print(result.output_paths["r2_hg38"])
```

The convenience wrapper creates the workflow log. Direct `ReferencePanelBuilder.run(config)` returns data paths without creating a workflow log by default.

## Advanced controls

| Option | Purpose and default |
| --- | --- |
| `--plink-prefix` | Required complete PLINK prefix or chromosome pattern |
| `--output-dir` | Required destination directory |
| `--source-genome-build` | `auto` by default; accepts hg19/hg37/GRCh37 or hg38/GRCh38 |
| `--snp-identifier` | `chr_pos_allele_aware` by default; also `chr_pos`, `rsid`, `rsid_allele_aware` |
| `--ld-wind-snps`, `--ld-wind-kb`, `--ld-wind-cm` | Exactly one window is required |
| `--genetic-map-hg19-sources`, `--genetic-map-hg38-sources` | Matching maps for sidecar CM; required for cM windows in each emitted build |
| `--liftover-chain-hg19-to-hg38-file`, `--liftover-chain-hg38-to-hg19-file` | Matching chain emits the opposite build in coordinate modes; invalid in rsID modes |
| `--ref-panel-snps-file` | Optional headered SNP restriction in the source build |
| `--keep-indivs-file` | Optional one-IID-per-row restriction |
| `--maf-min` | Optional MAF cutoff in [0, 0.5] |
| `--min-r2` | Optional adjusted-R² sparsification threshold; default 0 keeps all within-window pairs |
| `--snp-batch-size` | Genotype computation batch size; default 128 |
| `--overwrite` | Replace owned outputs and remove stale owned artifacts after success |
| `--log-level` | Workflow log verbosity: DEBUG, INFO (default), WARNING, ERROR |

Without SNP restriction, genotypes are streamed from disk. With an explicit restriction, the selected genotype payload is held in memory. Pair computation retains the current batch and window-spanning columns, so memory depends on retained individuals, SNPs, window density, and batch size. Wider windows also increase pair counts and runtime. There is no universal runtime estimate; benchmark one representative chromosome before scheduling a genome-wide build. See `ldsc build-r2-panel --help` and the [argument inventory](../docs/current/io-argument-inventory.md) for accepted path formats.

Sources: [`ReferencePanelBuilder.run` and CLI arguments](../src/ldsc/ref_panel_builder.py), [`write_r2_parquet` and `yield_pairwise_r2_rows`](../src/ldsc/_kernel/ref_panel_builder.py), and [`R2Panel.query_pairs`](../src/ldsc/r2_query.py).
