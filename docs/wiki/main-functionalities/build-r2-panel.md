# Build an R² Panel

Last updated on: 2026-09-11

`ldsc build-r2-panel` turns PLINK reference genotypes into reusable pairwise R² tables. These describe LD in the reference individuals; they are inputs to [LD-score calculation](ldscore.md) and [pair queries](../utility-functionalities/query-r2.md).

## Minimal command

Supply a complete `.bed/.bim/.fam` prefix and exactly one LD window:

```bash
ldsc build-r2-panel \
  --plink-prefix /data/reference/panel_chr22 \
  --source-genome-build hg38 \
  --ld-wind-kb 1000 \
  --output-dir results/r2_panel
```

This writes hg38 only, with missing genetic-map CM values unless a matching map is supplied. For a cM window, use `--ld-wind-cm` with `--genetic-map-hg38-sources`. A matching liftover chain in coordinate identity modes adds the opposite genome build; each emitted build needs its own map for a cM window.

## Outputs and interpretation

| Artifact | Contents |
| --- | --- |
| `hg38/chr22_r2.parquet` | One row per stored off-diagonal pair: `IDX_1`, `IDX_2`, `R2`, `SIGN_R` |
| `hg38/chr22_meta.tsv.gz` | Required index-to-SNP table: `CHR, POS, SNP, A1, A2, CM, MAF` |
| `diagnostics/metadata.chr22.json` | Run provenance for this concrete chromosome invocation |
| `diagnostics/build-r2-panel.chr22.log` | Workflow log |
| `diagnostics/dropped_snps/chr22_dropped.tsv.gz` | Dropped-SNP audit |

`IDX_1` and `IDX_2` are zero-based sidecar row numbers, with `IDX_1 < IDX_2`. Decode the int16 `R2` column by dividing by the footer's scale (32767). The value is bias-adjusted and can be negative. `SIGN_R` is boolean: true means nonnegative Pearson correlation in the panel's A1 dosage orientation. Keep each Parquet with its exact sidecar; it cannot be interpreted without that table.

A suite invocation writes one pair of scientific files per chromosome and uses unscoped `diagnostics/metadata.json` and `diagnostics/build-r2-panel.log`. `--log-level` controls log-file detail. Existing owned outputs require `--overwrite`.

## Controls and caveats

- `--ref-panel-snps-file`, `--keep-indivs-file`, and `--maf-min` explicitly restrict SNPs, individuals, or MAF. Match the reference individuals to the population whose LD you need.
- `--min-r2 0` keeps negative and zero estimates within the window. A positive threshold omits smaller values; this changes subsequent LD-score sums.
- Unstored pairs are unavailable in `query-r2`, not measured zero LD. Diagonal R²=1 is implicit.
- SNP restrictions use the PLINK source build. Coordinate duplicate groups use fixed `drop-all` handling; no duplicate-policy flag exists.
- Unrestricted builds stream genotypes; explicitly restricted builds keep the selected genotype payload in memory. Window density, sample count, and batch size affect resources. Benchmark a chromosome before a large run.
- The current names are `build-r2-panel` and `SIGN_R`; no old-name aliases are supported.

For every flag, a runnable fixture example, and Python usage, see the [build tutorial](../../../tutorials/build-parquet-reference-panel-from-plink.md) and [argument inventory](../../current/io-argument-inventory.md). The [storage specification](../../current/parquet-r2-format-and-read-pipeline.md) defines types, quantization, and binding.

Source: [`ReferencePanelBuilder.run` and command arguments](../../../src/ldsc/ref_panel_builder.py), and [`write_r2_parquet`](../../../src/ldsc/_kernel/ref_panel_builder.py).
