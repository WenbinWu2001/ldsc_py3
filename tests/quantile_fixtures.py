"""A four-SNP quantile-h2 problem with hand-worked answers.

The fitted model has the baseline annotations ``base`` (1 on every SNP) and
``cont`` (0, 1, 2, 3 on rs1 to rs4), with coefficients 0.05 and 0.10. The
target annotation repeats ``cont``, so two quantiles hold rs1-rs3 and rs4,
with observed-scale heritability 0.45 and 0.35. The LD-score directory is
written by the production writer; the fitted model and source tables are
written by hand.
"""

from __future__ import annotations

from argparse import Namespace
import json
from pathlib import Path

import numpy as np
import pandas as pd

from tests.ldscore_helpers import baseline_overlap, write_ldscore_dir


def write_quantile_ldscores(root, *, snp_identifier="rsid", genome_build=None, overlap=True, cont_common_count=6.0,
                            overwrite=False) -> Path:
    """Write ``root/ldscores`` for the baseline annotations ``base`` and ``cont`` of the four-SNP problem.

    The common reference-SNP counts are 4 for ``base`` and
    ``cont_common_count`` for ``cont``. With ``overlap``, the directory holds
    the cross-products ``[[4, 6], [6, 14]]`` and both SNP totals equal 4.
    Pass ``overwrite`` to replace an earlier fixture directory.
    """
    baseline = pd.DataFrame({
        "CHR": ["1"] * 4, "SNP": [f"rs{k}" for k in range(1, 5)], "POS": [10, 20, 30, 40],
        "A1": ["A"] * 4, "A2": ["C"] * 4, "regression_ld_scores": [1.0] * 4,
        "base": [1.0] * 4, "cont": [1.0] * 4,
    })
    products = pd.DataFrame([[4.0, 6.0], [6.0, 14.0]], index=["base", "cont"], columns=["base", "cont"])
    return write_ldscore_dir(
        Path(root) / "ldscores", baseline, baseline_columns=["base", "cont"],
        counts={"base": (4.0, 4.0), "cont": (cont_common_count, cont_common_count)},
        snp_identifier=snp_identifier, genome_build=genome_build,
        overlap=baseline_overlap(products, total_all_reference_snps=4, total_common_reference_snps=4)
        if overlap else None,
        annotation_types={"base": "binary", "cont": "quantitative"}, overwrite=overwrite,
    )


def write_quantile_inputs(root) -> Namespace:
    """Write the LD scores, fitted model and source tables under ``root`` and return the quantile-h2 arguments.

    The arguments write their result to ``root/quantile``.
    """
    root = Path(root)
    ldscore_dir = write_quantile_ldscores(root)
    fitted_dir = root / "fitted"
    (fitted_dir / "diagnostics").mkdir(parents=True)
    (fitted_dir / "diagnostics" / "metadata.json").write_text(
        json.dumps(
            {
                "artifact_type": "partitioned_h2_result",
                "analysis_type": "functional_category",
                "ldscore_dir": str(ldscore_dir),
                "retained_ld_columns": ["base", "cont"],
            }
        ),
        encoding="utf-8",
    )
    pd.DataFrame(
        {
            "category": ["base", "cont"],
            "coefficient": [0.05, 0.10],
            "coefficient_se": [0.01, 0.02],
            "samp_prev": [np.nan, np.nan],
            "pop_prev": [np.nan, np.nan],
        }
    ).to_csv(fitted_dir / "partitioned_h2.tsv", sep="\t", index=False, na_rep="NaN")
    pd.DataFrame(
        {
            "delete_block": [0, 1],
            "base": np.array([0.04, 0.06], dtype=np.float64),
            "cont": np.array([0.09, 0.11], dtype=np.float64),
        }
    ).to_parquet(fitted_dir / "diagnostics" / "coefficient_delete_values.parquet", index=False)
    annot = root / "annotations.tsv"
    annot.write_text(
        "CHR\tPOS\tSNP\tCM\tbase\tcont\n"
        "1\t10\trs1\t0\t1\t0\n"
        "1\t20\trs2\t0\t1\t1\n"
        "1\t30\trs3\t0\t1\t2\n"
        "1\t40\trs4\t0\t1\t3\n",
        encoding="utf-8",
    )
    target = root / "target.tsv"
    target.write_text(
        "CHR\tPOS\tSNP\ttarget\n"
        "1\t10\trs1\t0\n"
        "1\t20\trs2\t1\n"
        "1\t30\trs3\t2\n"
        "1\t40\trs4\t3\n",
        encoding="utf-8",
    )
    reference = root / "ref.tsv"
    reference.write_text(
        "CHR\tPOS\tSNP\tMAF\n"
        "1\t10\trs1\t0.2\n"
        "1\t20\trs2\t0.2\n"
        "1\t30\trs3\t0.2\n"
        "1\t40\trs4\t0.2\n",
        encoding="utf-8",
    )
    return Namespace(
        partitioned_h2_result_dir=str(fitted_dir),
        baseline_annot_sources=[str(annot)],
        query_annot_sources=None,
        query_annot_bed_sources=None,
        query_annot_gene_list_sources=None,
        gene_coordinate_file=None,
        control_gene_list_file=None,
        gene_list_resolution_policy="strict",
        gene_exclude_regions="none",
        padding_bp=0,
        target_annot_sources=[str(target)],
        target_annot_column="target",
        ref_metadata_sources=[str(reference)],
        target_missing_value=None,
        num_quantiles=2,
        output_dir=str(root / "quantile"),
        overwrite=False,
        log_level="INFO",
    )
