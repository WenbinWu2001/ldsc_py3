"""LD-score directory fixtures written through the production directory writer."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from ldsc.config import GlobalConfig
from ldsc._ldscore.outputs import LDScoreDirectoryWriter, LDScoreOutputConfig, common_count_config
from ldsc._ldscore.overlap_table import LDScoreOverlap
from ldsc._ldscore.results import LDScoreResult
from ldsc._ldscore.source import count_record


def edit_json(path, edit) -> None:
    """Apply ``edit`` to the parsed JSON document at ``path`` (an artifact's metadata) and save it."""
    path = Path(path)
    payload = json.loads(path.read_text(encoding="utf-8"))
    edit(payload)
    path.write_text(json.dumps(payload), encoding="utf-8")


def write_ldscore_dir(
    output_dir,
    baseline: pd.DataFrame,
    *,
    baseline_columns,
    counts,
    snp_identifier: str,
    genome_build: str | None = None,
    query: pd.DataFrame | None = None,
    query_columns=(),
    overlap: LDScoreOverlap | None = None,
    annotation_types=None,
    common_maf_min: float = 0.05,
    overwrite: bool = False,
) -> Path:
    """Publish one LD-score batch with ``LDScoreDirectoryWriter.write_batches`` and return the directory.

    ``baseline`` holds ``CHR``, ``SNP``, ``POS``, the alleles when the mode
    needs them, ``regression_ld_scores`` and every baseline column, in genomic
    order; ``query`` holds the same SNP rows with the query columns.
    ``counts`` maps each annotation to its ``(all, common)`` reference-SNP
    counts. The writer narrows float64 values to float32, records
    ``common_count_config(common_maf_min)`` and the overlap totals, and
    reloads the directory, so the fixture carries the real manifest.
    """
    baseline_columns, query_columns = list(baseline_columns), list(query_columns)
    records = [count_record("baseline" if column in baseline_columns else "query", column, *counts[column])
               for column in [*baseline_columns, *query_columns]]
    result = LDScoreResult(
        baseline_table=baseline,
        query_table=query,
        count_records=records,
        baseline_columns=baseline_columns,
        query_columns=query_columns,
        count_config=common_count_config(common_maf_min),
        config_snapshot=GlobalConfig(snp_identifier=snp_identifier, genome_build=genome_build),
        overlap=overlap,
        annotation_types=dict(annotation_types or {}),
    )
    LDScoreDirectoryWriter().write_batches([result], LDScoreOutputConfig(output_dir, overwrite=overwrite))
    return Path(output_dir)


def baseline_overlap(products, *, total_all_reference_snps, total_common_reference_snps) -> LDScoreOverlap:
    """Return a baseline-only overlap whose all- and common-SNP cross-products are both ``products``.

    ``products`` is a square frame indexed and labeled by the baseline annotations.
    """
    return LDScoreOverlap(
        baseline_block_all=products, baseline_block_common=products.copy(),
        query_diagonal_all=pd.Series(dtype="float64"), query_diagonal_common=pd.Series(dtype="float64"),
        total_all_reference_snps=float(total_all_reference_snps),
        total_common_reference_snps=float(total_common_reference_snps),
    )


def write_one_snp_ldscore_dir(root: Path, *, include_query: bool = True) -> Path:
    """Write a one-SNP rsID LD-score directory with alleles, with or without the query column ``query``.

    The baseline column ``base`` counts 5 all and 5 common reference SNPs;
    ``query`` counts 6 and 6. The directory records genome build ``None``,
    as every rsID-family LD-score directory does.
    """
    identity = {"CHR": ["1"], "SNP": ["rs1"], "POS": [10], "A1": ["A"], "A2": ["C"]}
    baseline = pd.DataFrame({**identity, "regression_ld_scores": [2.0], "base": [1.0]})
    query = pd.DataFrame({**identity, "query": [2.0]}) if include_query else None
    counts = {"base": (5.0, 5.0), "query": (6.0, 6.0)}
    return write_ldscore_dir(
        root, baseline, baseline_columns=["base"], query=query, query_columns=["query"] if include_query else [],
        counts=counts, snp_identifier="rsid",
    )
