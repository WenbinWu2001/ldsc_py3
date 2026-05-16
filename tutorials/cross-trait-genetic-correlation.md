# Cross-Trait Genetic Correlation

Goal: estimate genetic correlation for two or more traits from munged summary statistics and one matched LD-score reference.

The regression step expects:

- two or more curated `sumstats.parquet` files, or explicit legacy `.sumstats.gz`
  compatibility files, with at least `SNP`, `Z`, and `N`; current
  package-written artifacts also include `CHR` and `POS`, and `A1`/`A2` are
  recommended so each tested pair can be allele-aligned
- one canonical LD-score result directory containing `metadata.json` and `ldscore.baseline.parquet`

Use the same `GlobalConfig` assumptions for both traits and the LD-score reference. For cross-trait rg, the LD scores should usually be the baseline, non-cell-specific LD scores used for single-trait h2.
Current LD-score directories keep `ldscore.baseline.parquet` as one flat file
with chromosome-aligned row groups listed in root `metadata.json`; regression reads
the full table, while inspection code can use the row-group metadata.

Output directories are literal destinations. Missing directories are created,
existing directories are reused, and existing owned workflow artifacts are
refused before writing unless you pass `--overwrite` or `overwrite=True`.
For each trait's sumstats directory, that family includes `sumstats.parquet`,
`sumstats.sumstats.gz`, root `metadata.json`, `diagnostics/sumstats.log`, and
`diagnostics/dropped_snps/dropped.tsv.gz`, even when the current run would not write every
format. Successful overwrites remove
stale owned siblings not produced by the current configuration and preserve
unrelated files.

## Python API

```python
import pandas as pd

from ldsc import (
    GlobalConfig,
    MungeConfig,
    RegressionConfig,
    RegressionRunner,
    SumstatsMunger,
    load_ldscore_from_dir,
    load_sumstats,
    set_global_config,
)

GLOBAL_CONFIG = GlobalConfig(
    snp_identifier="chr_pos_allele_aware",
    genome_build="hg38",
)
set_global_config(GLOBAL_CONFIG)

# Option A: munge raw sumstats in the same workflow.
trait_1 = SumstatsMunger().run(
    MungeConfig(
        raw_sumstats_file="data/trait_1.tsv.gz",
        trait_name="trait_1",
        # Common SNP, allele, p-value, N, and signed-statistic headers are
        # inferred automatically. Use column_hints only for ambiguous files.
        # use_hm3_snps=True,  # packaged HM3 row restriction
        # target_genome_build="hg38",
        # use_hm3_quick_liftover=True,  # requires use_hm3_snps=True
        output_dir="tutorial_outputs/trait_1",
        # overwrite=True,  # also removes stale unselected sumstats sibling formats
    ),
    global_config=GLOBAL_CONFIG,
)

trait_2 = SumstatsMunger().run(
    MungeConfig(
        raw_sumstats_file="data/trait_2.tsv.gz",
        trait_name="trait_2",
        # Common SNP, allele, p-value, N, and signed-statistic headers are
        # inferred automatically. Use column_hints only for ambiguous files.
        # use_hm3_snps=True,  # packaged HM3 row restriction
        # target_genome_build="hg38",
        # use_hm3_quick_liftover=True,  # requires use_hm3_snps=True
        output_dir="tutorial_outputs/trait_2",
        # overwrite=True,  # also removes stale unselected sumstats sibling formats
    ),
    global_config=GLOBAL_CONFIG,
)

# Option B: load existing curated artifacts instead. Current artifacts recover
# config_snapshot from root metadata.json. Old package-written artifacts
# without current provenance must be regenerated with the current LDSC package.
# trait_1 = load_sumstats("tutorial_outputs/trait_1/sumstats.parquet", trait_name="trait_1")
# trait_2 = load_sumstats("tutorial_outputs/trait_2/sumstats.parquet", trait_name="trait_2")

ldscore_result = load_ldscore_from_dir(
    "tutorial_outputs/baseline_ldscores",
)

runner = RegressionRunner(
    global_config=GLOBAL_CONFIG,
    regression_config=RegressionConfig(),
)
result = runner.estimate_rg_pairs([trait_1, trait_2], ldscore_result)

result.rg.to_csv("tutorial_outputs/trait_1_trait_2_rg.tsv", sep="\t", index=False)
result.rg_full.to_csv("tutorial_outputs/trait_1_trait_2_rg_full.tsv", sep="\t", index=False)
result.h2_per_trait.to_csv("tutorial_outputs/trait_1_trait_2_h2_per_trait.tsv", sep="\t", index=False)
print(result.rg)
```

`RegressionRunner.estimate_rg()` remains available for low-level pairwise code
that needs the raw kernel object. For user-facing analyses, prefer
`estimate_rg_pairs()` because it returns the same output family as the CLI.

When both traits are produced by `SumstatsMunger.run()` in the same workflow,
their known `GlobalConfig` snapshots are checked against the LD-score snapshot
before regression. Current disk artifacts recover that snapshot from
root `metadata.json`; old package-written `.sumstats.gz` files without
current provenance must be regenerated rather than loaded with inferred
provenance. In `chr_pos_allele_aware` mode, both trait tables and the
LD-score table merge by normalized `CHR:POS:<allele_set>` identity, and
`SNP` is treated as a label. To run coordinate identity without allele-aware
matching, set `snp_identifier="chr_pos"` or pass `--snp-identifier chr_pos`;
use the base `rsid` mode for rsID-only identity. The removed `--no-alleles`
flag is not accepted.
The munger defaults to `--format auto`, and `--infer-only` can report
missing fields or exact repair flags without writing outputs. `A1` is the
allele that the signed statistic is relative to; `A2` is its counterpart.
`NEFF` is not inferred as total `N` automatically. Optional munger liftover is valid for chr_pos-family modes; use
`target_genome_build` with either a chain file or `use_hm3_snps=True` plus HM3 quick liftover when both
traits need to be converted to the LD-score build. Liftover drops duplicate
source/target coordinate groups, writes count summaries to `diagnostics/sumstats.log`, and
audits row-level drops in `diagnostics/dropped_snps/dropped.tsv.gz`; examples appear only at
`DEBUG`, not in the compatibility sidecar.

## CLI

```bash
ldsc munge-sumstats \
  --raw-sumstats-file data/trait_1.tsv.gz \
  --trait-name trait_1 \
  --use-hm3-snps \
  --output-dir tutorial_outputs/trait_1

ldsc munge-sumstats \
  --raw-sumstats-file data/trait_2.tsv.gz \
  --trait-name trait_2 \
  --use-hm3-snps \
  --output-dir tutorial_outputs/trait_2

ldsc rg \
  --sumstats-sources tutorial_outputs/trait_1/sumstats.parquet tutorial_outputs/trait_2/sumstats.parquet \
  --ldscore-dir tutorial_outputs/baseline_ldscores \
  --count-kind common \
  --output-dir tutorial_outputs/trait_1_trait_2
```

The command writes `tutorial_outputs/trait_1_trait_2/rg.tsv`,
`rg_full.tsv`, `h2_per_trait.tsv`, and `diagnostics/rg.log`.
If either fixed output already exists, the command fails before writing; add
`--overwrite` only when replacing the previous summary is intentional.

For a trait panel, pass a glob. Without an anchor this computes all unordered
pairs in input order:

```bash
ldsc rg \
  --sumstats-sources tutorial_outputs/traits/*.parquet \
  --ldscore-dir tutorial_outputs/baseline_ldscores \
  --output-dir tutorial_outputs/panel_rg \
  --write-per-pair-detail
```

To compare every trait against one anchor, add `--anchor-trait`. The anchor
matches a recovered trait label first, then a resolved file path:

```bash
ldsc rg \
  --sumstats-sources tutorial_outputs/traits/*.parquet \
  --anchor-trait trait_1 \
  --ldscore-dir tutorial_outputs/baseline_ldscores \
  --output-dir tutorial_outputs/trait_1_anchor_rg
```
