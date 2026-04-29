# Cross-Trait Genetic Correlation

Goal: estimate genetic correlation between two traits from munged summary statistics and one matched LD-score reference.

The regression step expects:

- two curated `.sumstats.gz` files with at least `SNP`, `Z`, and `N`; current
  package-written artifacts also include `CHR` and `POS`, and `A1`/`A2` are
  recommended so the second trait can be allele-aligned
- one canonical LD-score result directory containing `manifest.json` and `baseline.parquet`

Use the same `GlobalConfig` assumptions for both traits and the LD-score reference. For cross-trait rg, the LD scores should usually be the baseline, non-cell-specific LD scores used for single-trait h2.
Current LD-score directories keep `baseline.parquet` as one flat file with
chromosome-aligned row groups listed in `manifest.json`; regression reads the
full table, while inspection code can use the row-group metadata.

Output directories are literal destinations. Missing directories are created,
existing directories are reused, and existing fixed files such as
`sumstats.sumstats.gz`, `sumstats.log`, `sumstats.metadata.json`, or `rg.tsv`
are refused before writing unless you pass `--overwrite` or `overwrite=True`.

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
    snp_identifier="rsid",
    genome_build="hg38",
)
set_global_config(GLOBAL_CONFIG)

# Option A: munge raw sumstats in the same workflow.
trait_1 = SumstatsMunger().run(
    MungeConfig(
        sumstats_file="data/trait_1.tsv.gz",
        trait_name="trait_1",
        column_hints={"snp": "SNP", "a1": "A1", "a2": "A2", "p": "P", "N_col": "N"},
        # sumstats_snps_file="filters/hapmap3.tsv.gz",  # optional row keep-list
        output_dir="tutorial_outputs/trait_1",
        signed_sumstats_spec="BETA,0",
        # overwrite=True,  # enable only when intentionally replacing trait_1 outputs
    ),
    global_config=GLOBAL_CONFIG,
)

trait_2 = SumstatsMunger().run(
    MungeConfig(
        sumstats_file="data/trait_2.tsv.gz",
        trait_name="trait_2",
        column_hints={"snp": "SNP", "a1": "A1", "a2": "A2", "p": "P", "N_col": "N"},
        # sumstats_snps_file="filters/hapmap3.tsv.gz",  # optional row keep-list
        output_dir="tutorial_outputs/trait_2",
        signed_sumstats_spec="BETA,0",
        # overwrite=True,  # enable only when intentionally replacing trait_2 outputs
    ),
    global_config=GLOBAL_CONFIG,
)

# Option B: load existing curated artifacts instead. Current artifacts recover
# config_snapshot from sumstats.metadata.json. Older artifacts without that
# sidecar warn and use config_snapshot=None.
# trait_1 = load_sumstats("tutorial_outputs/trait_1/sumstats.sumstats.gz", trait_name="trait_1")
# trait_2 = load_sumstats("tutorial_outputs/trait_2/sumstats.sumstats.gz", trait_name="trait_2")

ldscore_result = load_ldscore_from_dir(
    "tutorial_outputs/baseline_ldscores",
)

runner = RegressionRunner(
    global_config=GLOBAL_CONFIG,
    regression_config=RegressionConfig(),
)
rg = runner.estimate_rg(trait_1, trait_2, ldscore_result)

summary = pd.DataFrame(
    [
        {
            "trait_1": trait_1.trait_name,
            "trait_2": trait_2.trait_name,
            "rg": getattr(rg, "rg_ratio", None),
            "rg_se": getattr(rg, "rg_se", None),
            "z": getattr(rg, "z", None),
            "p": getattr(rg, "p", None),
        }
    ]
)
summary.to_csv("tutorial_outputs/trait_1_trait_2_rg.tsv", sep="\t", index=False)
print(summary)
```

When both traits are produced by `SumstatsMunger.run()` in the same workflow,
their known `GlobalConfig` snapshots are checked against the LD-score snapshot
before regression. Current disk artifacts recover that snapshot from
`sumstats.metadata.json`; older `.sumstats.gz` files without the sidecar have
unknown provenance and skip sumstats-side compatibility validation rather than
fabricating a snapshot. In `chr_pos` mode, both trait tables and the LD-score
table merge by normalized `CHR:POS` coordinates.

## CLI

```bash
ldsc munge-sumstats \
  --sumstats-file data/trait_1.tsv.gz \
  --snp SNP \
  --a1 A1 \
  --a2 A2 \
  --p P \
  --N-col N \
  --sumstats-snps-file filters/hapmap3.tsv.gz \
  --signed-sumstats BETA,0 \
  --output-dir tutorial_outputs/trait_1

ldsc munge-sumstats \
  --sumstats-file data/trait_2.tsv.gz \
  --snp SNP \
  --a1 A1 \
  --a2 A2 \
  --p P \
  --N-col N \
  --sumstats-snps-file filters/hapmap3.tsv.gz \
  --signed-sumstats BETA,0 \
  --output-dir tutorial_outputs/trait_2

ldsc rg \
  --sumstats-1-file tutorial_outputs/trait_1/sumstats.sumstats.gz \
  --sumstats-2-file tutorial_outputs/trait_2/sumstats.sumstats.gz \
  --trait-name-1 trait_1 \
  --trait-name-2 trait_2 \
  --ldscore-dir tutorial_outputs/baseline_ldscores \
  --count-kind common \
  --output-dir tutorial_outputs/trait_1_trait_2
```

The command writes `tutorial_outputs/trait_1_trait_2/rg.tsv`.
If that TSV already exists, the command fails before writing; add `--overwrite`
only when replacing the previous summary is intentional.
