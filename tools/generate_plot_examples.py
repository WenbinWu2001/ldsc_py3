#!/usr/bin/env python3
"""Generate deterministic example LDSC result directories and all approved plots."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from ldsc import convert_h2_scale, plot_result


def _write_result(root: Path, table_name: str, table: pd.DataFrame, metadata: dict) -> None:
    root.mkdir(parents=True, exist_ok=True)
    (root / "diagnostics").mkdir(exist_ok=True)
    table.to_csv(root / table_name, sep="\t", index=False, na_rep="NaN")
    (root / "diagnostics" / "metadata.json").write_text(
        json.dumps(metadata, indent=2) + "\n",
        encoding="utf-8",
    )


def generate_examples(output_dir: Path) -> list[Path]:
    """Generate small canonical inputs and return the seven saved PNG paths."""
    output_dir.mkdir(parents=True, exist_ok=True)
    saved: list[Path] = []

    h2_root = output_dir / "h2"
    bins = pd.DataFrame(
        {
            "bin": np.arange(1, 13),
            "n_snps": [100] * 12,
            "ld_score_min": np.linspace(1.0, 19.0, 12),
            "ld_score_max": np.linspace(2.5, 22.0, 12),
            "mean_ld_score": np.linspace(1.8, 20.5, 12),
            "mean_chi_square": [1.03, 1.05, 1.07, 1.09, 1.12, 1.13, 1.17, 1.18, 1.21, 1.25, 1.27, 1.31],
            "sd_chi_square": [0.08] * 12,
            "mean_sample_size": [100_000] * 12,
            "mean_fitted_chi_square": np.linspace(1.025, 1.305, 12),
            "mean_regression_weight": np.linspace(0.2, 1.0, 12) ** 2,
        }
    )
    summary = pd.DataFrame(
        {
            "trait_name": ["Example case-control trait"],
            "total_h2_obs": [0.105],
            "total_h2_obs_se": [0.018],
        }
    )
    _write_result(
        h2_root,
        "h2.tsv",
        summary,
        {
            "artifact_type": "h2_result",
            "trait_name": "Example case-control trait",
            "files": {
                "summary": "h2.tsv",
                "ld_score_regression_bins": "diagnostics/ld_score_regression_bins.tsv",
            },
        },
    )
    bins.to_csv(h2_root / "diagnostics" / "ld_score_regression_bins.tsv", sep="\t", index=False)
    h2_plot = plot_result(h2_root)
    saved.append(h2_plot.path)
    conversion = convert_h2_scale(
        h2_root,
        samp_prev=0.50,
        pop_prev_range=(0.01, 0.20),
        num_points=201,
    )
    assert conversion.plot_path is not None
    saved.append(conversion.plot_path)

    functional_root = output_dir / "functional_partitioned_h2"
    _write_result(
        functional_root,
        "partitioned_h2.tsv",
        pd.DataFrame(
            {
                "category": ["Coding", "Conserved", "Promoter", "Enhancer", "H3K27ac", "Intron"],
                "enrichment": [6.8, 4.7, 2.2, 1.7, 1.4, 0.8],
                "enrichment_se": [1.1, 0.7, 0.4, 0.3, 0.25, 0.15],
            }
        ),
        {
            "artifact_type": "partitioned_h2_result",
            "analysis_type": "functional_category",
            "headline_metric": "enrichment",
            "files": {"summary": "partitioned_h2.tsv"},
        },
    )
    functional_plot = plot_result(functional_root)
    saved.append(functional_plot.path)

    query_root = output_dir / "cell_type_partitioned_h2"
    _write_result(
        query_root,
        "partitioned_h2.tsv",
        pd.DataFrame(
            {
                "category": ["Microglia", "Excitatory neuron", "Oligodendrocyte", "Astrocyte", "Endothelial"],
                "coefficient_p": [2e-8, 6e-5, 0.012, 0.09, 0.42],
                "coefficient_z": [5.49, 3.85, 2.26, 1.34, 0.20],
            }
        ),
        {
            "artifact_type": "partitioned_h2_result",
            "analysis_type": "cell_type_specific",
            "headline_metric": "coefficient",
            "files": {"summary": "partitioned_h2.tsv"},
        },
    )
    query_plot = plot_result(query_root)
    saved.append(query_plot.path)

    quantile_root = output_dir / "quantile_h2"
    _write_result(
        quantile_root,
        "quantile_h2.tsv",
        pd.DataFrame(
            {
                "quantile": [1, 2, 3, 4, 5],
                "target_value_lower": [0.00, 0.20, 0.40, 0.60, 0.80],
                "target_value_upper": [0.20, 0.40, 0.60, 0.80, 1.00],
                "enrichment": [0.62, 0.83, 1.06, 1.48, 2.15],
                "enrichment_se": [0.11, 0.12, 0.14, 0.20, 0.29],
            }
        ),
        {
            "artifact_type": "quantile_h2_result",
            "target_annotation": "Conservation score",
            "files": {"quantile_h2": "quantile_h2.tsv"},
        },
    )
    quantile_plot = plot_result(quantile_root)
    saved.append(quantile_plot.path)

    traits = ["Trait A", "Trait B", "Trait C", "Trait D"]
    rg_root = output_dir / "rg_all_pairs"
    _write_result(
        rg_root,
        "rg.tsv",
        pd.DataFrame(
            {
                "trait_1": ["Trait A", "Trait A", "Trait A", "Trait B", "Trait B", "Trait C"],
                "trait_2": ["Trait B", "Trait C", "Trait D", "Trait C", "Trait D", "Trait D"],
                "rg": [0.62, -0.31, 0.18, -0.12, np.nan, 0.47],
                "rg_se": [0.09, 0.08, 0.07, 0.10, np.nan, 0.11],
            }
        ),
        {
            "artifact_type": "rg_result",
            "pair_kind": "all_pairs",
            "trait_names": traits,
            "files": {"rg": "rg.tsv"},
        },
    )
    heatmap = plot_result(rg_root)
    saved.append(heatmap.path)

    anchor_root = output_dir / "rg_anchor"
    _write_result(
        anchor_root,
        "rg.tsv",
        pd.DataFrame(
            {
                "trait_1": ["Trait A", "Trait A", "Trait A"],
                "trait_2": ["Trait B", "Trait C", "Trait D"],
                "rg": [0.62, -0.31, 0.18],
                "rg_se": [0.09, 0.08, 0.07],
            }
        ),
        {
            "artifact_type": "rg_result",
            "pair_kind": "anchor",
            "trait_names": traits,
            "files": {"rg": "rg.tsv"},
        },
    )
    anchor_plot = plot_result(anchor_root)
    saved.append(anchor_plot.path)
    return saved


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output_dir", type=Path, help="Empty destination for deterministic example result trees.")
    args = parser.parse_args()
    for path in generate_examples(args.output_dir):
        print(path)


if __name__ == "__main__":
    main()
