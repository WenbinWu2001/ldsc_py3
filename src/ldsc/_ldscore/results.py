"""LD-score result records shared by the producers, the writer, regression and the converter.

``ChromLDScoreResult`` and ``LDScoreResult`` are the public, materialized
per-chromosome and cross-chromosome results. ``DirectChromosomeFragment`` and
``StagedLDScoreBatch`` are their private file-backed counterparts: a writing
run's chromosome tasks stage query scores in Parquet fragments, and the staged
batch carries the ordered fragment paths to the writer. These records hold
data and check their own shape; computing, aggregating and writing them
belong to the LD-score producers (``_ldscore.direct_batches``,
``_ldscore.gene_index.batches``) and ``_ldscore.outputs``.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING, Any

import numpy as np
import pandas as pd

from .._core.row_alignment import assert_same_snp_rows
from ..errors import LDSCInputError
from .source import validate_ldscore_schemas

if TYPE_CHECKING:
    from .._annotation.overlap import OverlapContribution
    from .._annotation.query_outcomes import QueryAnnotationStatus
    from .._annotation.storage import TsvDiagnostics
    from ..config import GlobalConfig
    from .overlap_table import LDScoreOverlap
    from .query_variation import ScoreVariation


@dataclass(frozen=True)
class ChromLDScoreResult:
    """Normalized per-chromosome LD-score output in split baseline/query form.

    Parameters
    ----------
    chrom : str
        Chromosome label.
    baseline_table : pandas.DataFrame
        Table with ``CHR``, ``SNP``, ``POS``, ``regression_ld_scores``, and baseline
        LD-score columns.
    query_table : pandas.DataFrame or None
        Optional table with ``CHR``, ``SNP``, ``POS``, and query LD-score
        columns.
    baseline_columns, query_columns : list of str
        Ordered annotation LD-score columns in the baseline and query tables.
    snp_count_totals : dict of str to numpy.ndarray, optional
        All-SNP (``all_reference_snp_counts``) and common-SNP
        (``common_reference_snp_counts``) count vectors over the reference
        SNPs, one entry per baseline column followed by each query column.
    annotation_types : dict of str to str, optional
        Per-column ``binary`` or ``quantitative`` classification used only for
        interpretation messages and persisted provenance.
    reference_snp_funnel : pandas.DataFrame or None, optional
        SNPs removed by each reference-preparation step on this chromosome;
        see ``ldsc._refpanel.funnel``.
    """
    chrom: str
    baseline_table: pd.DataFrame
    query_table: pd.DataFrame | None
    baseline_columns: list[str]
    query_columns: list[str]
    snp_count_totals: dict[str, np.ndarray] = field(default_factory=dict, repr=False)
    config_snapshot: GlobalConfig | None = None
    overlap: OverlapContribution | None = field(default=None, repr=False)
    reference_snp_count: int = 0
    regression_selected_snp_count: int = 0
    regression_region_removed_snp_count: int = 0
    identity_drops: TsvDiagnostics | None = field(default=None, repr=False)
    annotation_types: dict[str, str] = field(default_factory=dict)
    reference_snp_funnel: pd.DataFrame | None = field(default=None, repr=False)


@dataclass(frozen=True)
class LDScoreResult:
    """Materialized LD scores for one cross-chromosome calculation batch.

    ``LDScoreCalculator.run`` returns this object for a single prepared-input
    batch without an output configuration. Writing workflows return
    ``LDScoreSource`` so completed query tables can be released.

    Parameters
    ----------
    baseline_table : pandas.DataFrame
        Cross-chromosome table persisted as ``ldscore.baseline.parquet`` when
        the result is written. Required columns are ``CHR``, ``SNP``, ``POS``,
        ``regression_ld_scores``, and every entry in ``baseline_columns``.
    query_table : pandas.DataFrame or None
        Optional materialized cross-chromosome query table. A standalone write
        saves it as ``ldscore.query.parquet``; sequential output uses the
        applicable numbered batch filename. Required columns are ``CHR``,
        ``SNP``, ``POS``, and every entry in ``query_columns``.
    count_records : list of dict
        Manifest-ready count records. Each record names an annotation column and
        its all-SNP and optional common-SNP counts.
    count_config : dict, optional
        Manifest-ready metadata describing how optional common-SNP counts were
        computed, including the inclusive common-MAF threshold.
    baseline_columns, query_columns : list of str
        Ordered LD-score columns available to regression workflows.
    output_paths : dict, optional
        Paths written by ``LDScoreDirectoryWriter`` when output was requested.
    config_snapshot : GlobalConfig or None, optional
        Shared configuration active when the result was computed.
    query_statuses : tuple of QueryAnnotationStatus, optional
        Ordered per-source outcomes for BED and gene-list query annotations.
    gene_list_batch : object or None, optional
        Compact gene-list summaries and a streamed audit. Completed command
        results reference persistent diagnostics.
    identity_drops_by_chrom : dict, optional
        Chromosome diagnostic artifacts with explicit ``path`` and ``frames()``
        access. Published results reference persistent files. A calculator run
        without output configuration owns complete in-memory diagnostics,
        which remain available after the input bundle is closed.
    chromosome_scope : dict, optional
        Validated input chromosome coverage and effective analysis scope.
        Empty for workflows that do not perform query coverage preflight.
    legacy_ldscore_conversion : dict or None, optional
        Explicit converter provenance for a legacy LD-score suite. ``None``
        for natively computed LDSC3 results.
    ld_window : dict or None, optional
        The run's LD window and whole-chromosome override, recorded in
        ``metadata.json``: ``unit``, ``value`` and ``whole_chromosome_ok``.
        Indexed runs take the index's window. ``None`` for converted legacy
        LD scores, whose window is unknown.
    annotation_types : dict of str to str, optional
        Per-column semantic classification used for interpretation only.
    ld_window_snp_counts : pandas.DataFrame or None, optional
        For an R2 matrix with ``--ld-wind-snps``, one row per chromosome
        with ``CHR``, ``panel_snps`` (sidecar rows), ``kept_snps`` (SNPs this
        run kept), ``ld_wind_snps``, and ``counting_basis="panel_rows"``.
        ``None`` for kb/cM windows and PLINK reference panels.
    reference_snp_funnel : pandas.DataFrame or None, optional
        Reference-SNP funnel of every chromosome, including skipped ones: the
        SNPs each reference-preparation step removed, in the backend's order.
        Columns are ``funnel`` (``reference_panel`` or ``annotation_snp_grid``),
        ``chrom``, ``backend``, ``step_index``, ``step``, ``removed`` and
        ``remaining``; ``removed`` is missing for a step that did not run.
        Writing runs save it as ``diagnostics/reference_snp_funnel.tsv``.
        ``None`` when the reference panel records no funnel.
    """
    baseline_table: pd.DataFrame
    query_table: pd.DataFrame | None
    count_records: list[dict[str, Any]]
    baseline_columns: list[str]
    query_columns: list[str]
    output_paths: dict[str, str] = field(default_factory=dict)
    count_config: dict[str, Any] = field(default_factory=dict)
    config_snapshot: GlobalConfig | None = None
    overlap: "LDScoreOverlap | None" = field(default=None, repr=False)
    query_statuses: tuple[QueryAnnotationStatus, ...] = ()
    gene_list_batch: Any | None = field(default=None, repr=False)
    snp_universe_policy: dict[str, Any] | None = None
    index_provenance: dict[str, str] | None = None
    legacy_ldscore_conversion: dict[str, Any] | None = None
    ld_window: dict[str, Any] | None = None
    identity_drops_by_chrom: dict[str, TsvDiagnostics] = field(default_factory=dict, repr=False)
    annotation_types: dict[str, str] = field(default_factory=dict)
    chromosome_scope: dict[str, Any] = field(default_factory=dict)
    ld_window_snp_counts: pd.DataFrame | None = field(default=None, repr=False)
    reference_snp_funnel: pd.DataFrame | None = field(default=None, repr=False)

    def validate(self) -> None:
        """Check the normalized public contract for aggregated results."""
        validate_ldscore_schemas(
            self.baseline_columns, self.query_columns, self.baseline_table.columns,
            None if self.query_table is None else self.query_table.columns,
        )
        if self.query_table is not None:
            assert_same_snp_rows(
                self.baseline_table,
                self.query_table,
                context="query rows must match baseline rows on CHR/SNP/POS",
                snp_identifier=getattr(self.config_snapshot, "snp_identifier", "chr_pos_allele_aware"),
            )

    def read_queries(self, columns):
        """Select named queries from this materialized result without file I/O.

        Parameters
        ----------
        columns : sequence of str
            Query names in the desired output column order.

        Returns
        -------
        pandas.DataFrame or None
            SNP metadata followed by selected query values, in existing row
            order with a reset index. An empty selection returns metadata
            only, or ``None`` when this result has no query table.

        Raises
        ------
        LDSCInputError
            Nonempty query selection is requested from a baseline-only result.
        KeyError
            A selected name is absent from the materialized table.
        """
        if self.query_table is None:
            if columns:
                raise LDSCInputError("LD-score result has no query table.")
            return None
        metadata = [name for name in ("CHR", "SNP", "POS", "A1", "A2") if name in self.query_table]
        selected = self.query_table.loc[:, metadata + list(columns)]
        if not isinstance(selected.index, pd.RangeIndex) or selected.index.start != 0:
            selected = selected.reset_index(drop=True)
        return selected


@dataclass(frozen=True)
class DirectChromosomeFragment:
    """Direct chromosome statistics with a private query path instead of scores."""

    chrom: str
    baseline_table: pd.DataFrame | None
    query_fragment: Path | None
    baseline_columns: list[str]
    query_columns: list[str]
    snp_count_totals: dict
    config_snapshot: Any
    overlap: Any
    reference_snp_count: int
    regression_selected_snp_count: int
    regression_region_removed_snp_count: int
    identity_drops: Any
    annotation_types: dict
    variation: ScoreVariation
    reference_snp_funnel: pd.DataFrame | None = None


@dataclass(frozen=True)
class StagedLDScoreBatch:
    """Metadata and ordered query paths for one private LD-score batch.

    ``query_fragments`` retain direct post-kernel float32 scores or indexed
    float64 scores before public narrowing, in chromosome/SNP order.
    ``query_columns`` are finalized using variation summaries in the respective
    LD-score producer's precision; reads select only these retained columns.
    This internal descriptor is neither a materialized ``LDScoreResult`` nor
    a persistent ``LDScoreSource`` and must not escape its producer lifetime.
    """

    baseline_table: pd.DataFrame
    query_fragments: tuple[Path, ...]
    baseline_columns: list[str]
    query_columns: list[str]
    count_records: list[dict]
    count_config: dict
    config_snapshot: Any
    overlap: Any
    query_statuses: tuple
    gene_list_batch: Any
    chromosome_scope: dict
    index_provenance: dict
    annotation_types: dict = field(default_factory=dict)
    identity_drops_by_chrom: dict = field(default_factory=dict)
    snp_universe_policy: dict | None = None
    ld_window: dict | None = None
    ld_window_snp_counts: pd.DataFrame | None = None
    reference_snp_funnel: pd.DataFrame | None = None
