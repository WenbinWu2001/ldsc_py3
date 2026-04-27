# Cross-Trait Genetic Correlation

Goal: estimate genetic correlation between two traits from munged summary statistics and one matched LD-score reference.

The regression step expects:

- two curated `.sumstats.gz` files with at least `SNP`, `Z`, and `N`; `A1`/`A2` are recommended so the second trait can be allele-aligned
- one merged `.l2.ldscore.gz` table with an embedded `regr_weight` column, or a legacy separate `.w.l2.ldscore.gz` file
- one count vector such as `.l2.M_5_50`

Use the same `GlobalConfig` assumptions for both traits and the LD-score reference. For cross-trait rg, the LD scores should usually be the baseline, non-cell-specific LD scores used for single-trait h2.

## Python API

```python
import pandas as pd

from ldsc import (
    GlobalConfig,
    MungeConfig,
    RawSumstatsSpec,
    RegressionConfig,
    RegressionRunner,
    SumstatsMunger,
    load_ldscore_from_files,
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
    RawSumstatsSpec(
        path="data/trait_1.tsv.gz",
        trait_name="trait_1",
        column_hints={"snp": "SNP", "a1": "A1", "a2": "A2", "p": "P", "N_col": "N"},
    ),
    MungeConfig(
        out_prefix="tutorial_outputs/trait_1",
        signed_sumstats_spec="BETA,0",
    ),
    global_config=GLOBAL_CONFIG,
)

trait_2 = SumstatsMunger().run(
    RawSumstatsSpec(
        path="data/trait_2.tsv.gz",
        trait_name="trait_2",
        column_hints={"snp": "SNP", "a1": "A1", "a2": "A2", "p": "P", "N_col": "N"},
    ),
    MungeConfig(
        out_prefix="tutorial_outputs/trait_2",
        signed_sumstats_spec="BETA,0",
    ),
    global_config=GLOBAL_CONFIG,
)

# Option B: load existing curated artifacts instead.
# trait_1 = load_sumstats("tutorial_outputs/trait_1.sumstats.gz", trait_name="trait_1")
# trait_2 = load_sumstats("tutorial_outputs/trait_2.sumstats.gz", trait_name="trait_2")

ldscore_result = load_ldscore_from_files(
    ldscore_path="tutorial_outputs/baseline_ldscores.merged.l2.ldscore.gz",
    counts_path="tutorial_outputs/baseline_ldscores.l2.M_5_50",
    count_kind="m_5_50",
    snp_identifier=GLOBAL_CONFIG.snp_identifier,
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

If your LD-score table does not contain `regr_weight`, add `weight_path="tutorial_outputs/baseline_ldscores.w.l2.ldscore.gz"` to `load_ldscore_from_files(...)`.

## CLI

```bash
ldsc munge-sumstats \
  --sumstats data/trait_1.tsv.gz \
  --snp SNP \
  --a1 A1 \
  --a2 A2 \
  --p P \
  --N-col N \
  --signed-sumstats BETA,0 \
  --out tutorial_outputs/trait_1

ldsc munge-sumstats \
  --sumstats data/trait_2.tsv.gz \
  --snp SNP \
  --a1 A1 \
  --a2 A2 \
  --p P \
  --N-col N \
  --signed-sumstats BETA,0 \
  --out tutorial_outputs/trait_2

ldsc rg \
  --sumstats-1 tutorial_outputs/trait_1.sumstats.gz \
  --sumstats-2 tutorial_outputs/trait_2.sumstats.gz \
  --trait-name-1 trait_1 \
  --trait-name-2 trait_2 \
  --ldscore tutorial_outputs/baseline_ldscores.merged.l2.ldscore.gz \
  --counts tutorial_outputs/baseline_ldscores.l2.M_5_50 \
  --count-kind m_5_50 \
  --out tutorial_outputs/trait_1_trait_2
```

The command writes `tutorial_outputs/trait_1_trait_2.rg.tsv`. Add `--w-ld <path>` only for legacy LD-score outputs that keep regression weights in a separate file.
