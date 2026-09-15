"""regression_runner.py

Core functionality:
    Load normalized LD-score artifacts, assemble regression-ready datasets, and
    dispatch heritability or genetic-correlation estimators.

Overview
--------
This module is the workflow-layer boundary for the LDSC regression commands.
It consumes the canonical LD-score result directory introduced by the IO
refactor: ``metadata.json`` plus split baseline/query parquet tables. Public
callers pass one ``--ldscore-dir`` input instead of individual LD-score,
weight, count, and annotation-manifest files.

Legacy LDSC2 ``.sumstats`` and ``.sumstats.gz`` inputs are accepted at this
boundary. Their rsIDs are lookup keys: each trait is projected onto the
canonical LD-score panel, allele orientation is validated and harmonized, and
rejected rows are recorded in a stable diagnostic audit. Legacy LD-score
fragments are never read here and require explicit conversion first.

Regression chooses either the metadata ``common_reference_snp_counts`` vector or
``all_reference_snp_counts`` through ``--count-kind common|all``. The default
``common`` mode falls back to all-SNP counts when common counts are unavailable.
Regression merges by the effective key for the resolved SNP identifier mode.
Base modes are allele-blind; allele-aware/base mixes can run under the base mode
only when ``--allow-identity-downgrade`` or
``RegressionConfig.allow_identity_downgrade`` is set, and only within the same
rsID or coordinate family.
Partitioned-h2 fits the complete baseline model for baseline-only inputs and
one baseline-plus-query model per query annotation otherwise. Query runs share
trait/baseline preparation, load query columns in bounded batches, and fit
inline by default or in bounded whole-query processes via ``threads``. Completed
category/delete-value/metadata files are staged immediately. The output-layer
``PartitionedH2DirectoryWriter`` publishes them after final summary sorting.
Strict batches require every query to succeed. Explicit continuation records
failed queries and publishes successful models with an attempted-query ledger.
Genetic correlation accepts two or more munged summary-statistic sources and
returns the full rg output family:
the concise headline table, a diagnostic full table, per-trait h2 summaries,
and optional per-pair metadata for filesystem detail outputs.

Regression commands require an ``output_dir`` and create per-run logs under
``diagnostics/``. Individual estimator methods remain in-memory APIs; batch
partitioned fitting requires an output directory and returns persistent detail
paths. Final HM3 LD tables remain aggregate files.
"""

from __future__ import annotations

from contextlib import contextmanager, nullcontext
import hashlib
import json
import logging
import math
import time
from dataclasses import dataclass, field, replace
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Any, Sequence

import numpy as np
import pandas as pd

from ._cli_help import CLIHelpFormatter, SCALAR_PATH_HELP
from ._logging import LOG_LEVEL_HELP
from ._coordinates import CHR_POS_KEY_COLUMN, build_chr_pos_key_frame
from .config import (
    ConfigMismatchError,
    GlobalConfig,
    RegressionConfig,
    get_global_config,
    print_global_config_banner,
    suppress_global_config_banner,
)
from .path_resolution import ensure_output_directory, normalize_path_token
from ._logging import log_inputs, log_outputs, materializing_overwrite_guard, workflow_logging
from ._kernel import regression as reg
from ._kernel._jackknife import JackknifeIdentifiabilityError
from ._kernel.identifiers import build_snp_id_series
from ._kernel.snp_identity import (
    clean_identity_artifact_table,
    effective_merge_key_series,
    identity_base_mode,
    identity_mode_family,
    is_allele_aware_mode,
    resolve_regression_identity_mode,
    validate_identity_artifact_metadata,
)
from ._row_alignment import assert_same_snp_rows
from ._parallelism import _parse_threads, _resolve_worker_count, _validate_threads
from .column_inference import infer_chr_pos_columns, normalize_snp_identifier_mode
from .genome_build_inference import infer_chr_pos_build
from .ldscore_calculator import LDScoreResult
from .ldscore_source import LDScoreSource
from .outputs import (
    H2DirectoryWriter,
    H2OutputConfig,
    H2_REGRESSION_BIN_COLUMNS,
    PARTITIONED_H2_COLUMNS,
    PartitionedH2DirectoryWriter,
    PartitionedH2OutputConfig,
    PartitionedH2FitArtifacts,
    stage_partitioned_h2_fit,
    REGRESSION_LD_SCORE_COLUMN,
    RG_CONCISE_COLUMNS,
    RG_FULL_COLUMNS,
    RgDirectoryWriter,
    RgOutputConfig,
    _validate_ldscore_allele_columns,
)
from ._annotation_storage import AnnotationWorkspace
from .sumstats_munger import SumstatsTable, load_sumstats
from .errors import LDSCInputError, LDSCInternalError, LDSCUsageError, LDSCUserError


COMMON_COUNT_KEY = "common_reference_snp_counts"
ALL_COUNT_KEY = "all_reference_snp_counts"
LOGGER = logging.getLogger("LDSC.regression_runner")
PARTITIONED_H2_SUMMARY_SORT_COLUMNS = {
    "category": "category",
    "prop-snps": "prop_snps",
    "prop-h2": "prop_h2",
    "enrichment": "enrichment",
    "enrichment-p": "enrichment_p",
    "coefficient": "coefficient",
    "coefficient-p": "coefficient_p",
}
PARTITIONED_H2_SUMMARY_ASCENDING_SORTS = {
    "enrichment-p",
    "coefficient-p",
}
PARTITIONED_H2_REQUIRES_OVERLAP_MESSAGE = (
    "partitioned-h2 needs the annotation overlap matrix, but the LD-score directory has no "
    "ldscore.overlap.parquet. Most likely it is an unpartitioned LD-score directory (a single "
    "annotation such as the synthetic `base`), which has no overlap matrix and cannot be "
    "partitioned -- regenerate it with two or more annotation columns. It may also predate the "
    "overlap sidecar; if so, regenerate with the current `ldsc ldscore`. "
    "Other causes & fixes: docs/troubleshooting.md#partitioned-h2-missing-overlap-matrix"
)
FAILED_RG_NOTE = "Failed; see rg_full.tsv error column; use --output-dir for diagnostics/rg.log."
REGRESSION_IDENTITY_KEY_COLUMN = "_ldsc_regression_identity_key"
LEGACY_SUMSTATS_DROP_COLUMNS = [
    "trait_name",
    "source_path",
    "SNP",
    "A1",
    "A2",
    "reason",
    "panel_candidate_count",
]
_REGRESSION_NO_OVERLAP_DOC = (
    "docs/troubleshooting.md#regression-no-snps-remain-after-merging-regression-inputs"
)
_REGRESSION_SCHEMA_DOC = "docs/troubleshooting.md#common-ldsc-artifact-schema-or-provenance-is-incompatible"


@contextmanager
def _log_phase_timing(phase: str):
    """Log elapsed wall time for one observable workflow phase."""
    started = time.perf_counter()
    LOGGER.info("Phase start: %s.", phase)
    try:
        yield
    except Exception:
        LOGGER.info("Phase timing: %s failed after %.3fs.", phase, time.perf_counter() - started)
        raise
    LOGGER.info("Phase timing: %s completed in %.3fs.", phase, time.perf_counter() - started)


@dataclass(frozen=True)
class RegressionDataset:
    """Merged regression-ready dataset built from sumstats and LD-score tables.

    ``legacy_sumstats_drops`` carries stable row-level projection diagnostics
    when the source was LDSC2 text; it is empty for native LDSC3 inputs.
    """
    merged: pd.DataFrame
    ref_ld_columns: list[str]
    weight_column: str
    reference_snp_count_totals: dict[str, np.ndarray]
    count_key_used_for_regression: str
    retained_ld_columns: list[str]
    dropped_zero_variance_ld_columns: list[str]
    trait_names: list[str]
    chromosomes_aggregated: list[str]
    config_snapshot: GlobalConfig | None = None
    effective_snp_identifier: str = "chr_pos_allele_aware"
    identity_downgrade_applied: bool = False
    ldscore_overlap: "LDScoreOverlap | None" = None
    legacy_sumstats_drops: pd.DataFrame = field(
        default_factory=lambda: pd.DataFrame(columns=LEGACY_SUMSTATS_DROP_COLUMNS)
    )

    def validate(self) -> None:
        """Validate that the merged table contains the required LDSC columns."""
        required = {"SNP", self.weight_column, "Z", "N"}
        missing = required - set(self.merged.columns)
        if missing:
            raise LDSCInternalError(
                f"h2 regression dataset validation failed: merged table is missing columns {sorted(missing)}. "
                "Most likely regression preprocessing dropped required sumstats or LD-score columns. "
                "Re-run with `--log-level DEBUG` and report the traceback."
            )


@dataclass(frozen=True)
class PreparedRegressionInputs:
    """Shared trait alignment before model-dependent column/count filtering.

    ``dataset`` carries baseline LD values and aligned traits in genomic order.
    ``source_rows`` maps these rows back to the original LD-score table, so
    query batches reuse exactly the same join and allele orientation.
    """

    dataset: RegressionDataset
    source_rows: np.ndarray
    ldscore: LDScoreResult | LDScoreSource
    ldscore_identity_mode: str
    counts_by_column: dict[str, dict]


@dataclass(frozen=True)
class _PartitionedModelInputs:
    """Aligned model state without input readers or the full query source."""

    dataset: RegressionDataset
    counts_by_column: dict[str, dict]


@dataclass(frozen=True)
class RGRegressionDataset:
    """Merged genetic-correlation dataset built from two traits and LD scores.

    ``legacy_sumstats_drops`` combines projection diagnostics for both traits
    before pairwise allele harmonization.
    """
    merged: pd.DataFrame
    ref_ld_columns: list[str]
    weight_column: str
    reference_snp_count_totals: dict[str, np.ndarray]
    count_key_used_for_regression: str
    retained_ld_columns: list[str]
    dropped_zero_variance_ld_columns: list[str]
    trait_names: list[str]
    chromosomes_aggregated: list[str]
    config_snapshot: GlobalConfig | None = None
    effective_snp_identifier: str = "chr_pos_allele_aware"
    identity_downgrade_applied: bool = False
    legacy_sumstats_drops: pd.DataFrame = field(
        default_factory=lambda: pd.DataFrame(columns=LEGACY_SUMSTATS_DROP_COLUMNS)
    )

    def validate(self) -> None:
        """Validate that the merged table contains the required RG columns."""
        required = {"SNP", self.weight_column, "Z1", "N1", "Z2", "N2"}
        missing = required - set(self.merged.columns)
        if missing:
            raise LDSCInternalError(
                f"rg regression dataset validation failed: merged table is missing columns {sorted(missing)}. "
                "Most likely regression preprocessing dropped required trait or LD-score columns. "
                "Re-run with `--log-level DEBUG` and report the traceback."
            )


@dataclass(frozen=True)
class _FitOutcome:
    """Estimator and the realized population supplied to its final fit.

    The fitted dataset preserves SNP order, annotation order, reference counts
    and identity provenance. It reflects the workflow's chi-square filter;
    the two-step estimator's first-stage subset is not the final population.
    Reporting consumes these facts instead of reconstructing selection rules.
    """

    estimator: reg.Hsq | reg.RG
    dataset: RegressionDataset | RGRegressionDataset
    n_blocks: int
    effective_chisq_max: float | None
    diagnostic_bins: pd.DataFrame | None = None

    @property
    def n_snps(self) -> int:
        """Number of rows actually supplied to the final estimator."""
        return len(self.dataset.merged)


@dataclass(frozen=True)
class PartitionedH2BatchResult:
    """Compact batch summary and valid persistent artifacts for every fitted model.

    Category tables and coefficient delete values are written after each fit.
    Only the summary and path descriptors remain when the batch returns.
    ``query_status`` contains one row per attempted query in input order, with
    status, failed stage, exception type, and message. Failed queries have no
    scientific summary or model artifacts. Baseline-only runs leave it empty.
    Its columns are ``query_annotation``, ``status``, ``stage``, ``error_type``,
    and ``error_message``. Status is ``success``, ``unestimable`` for singular
    jackknife deletions, or ``failed`` for other per-query exceptions. Successful
    query scans also persist this table as ``diagnostics/query_status.tsv``.
    Worker completion order does not affect status order, stable summary
    sorting, or the final per-query folder order. Worker-count provenance is
    recorded in root and successful per-query metadata, not scientific tables.
    """

    summary: pd.DataFrame
    output_paths: dict[str, str] = field(default_factory=dict)
    per_query_artifacts: dict[str, PartitionedH2FitArtifacts] = field(default_factory=dict)
    query_status: pd.DataFrame = field(default_factory=pd.DataFrame)


@dataclass(frozen=True)
class RgResultFamily:
    """Result family returned by multi-trait genetic-correlation workflows.

    Parameters
    ----------
    rg : pandas.DataFrame
        Concise one-row-per-pair table with the publication-oriented rg schema:
        trait names, SNP count, :math:`r_g`, standard error, nominal p-value,
        and a failure note.
    rg_full : pandas.DataFrame
        Comprehensive diagnostic table with one row per attempted pair. Failed
        pairs keep their row, set numeric fields to NaN, and store
        ``status='failed'`` plus an ``error`` message.
    h2_per_trait : pandas.DataFrame
        One total-heritability summary per input trait, computed once per trait
        on the trait's single-trait LDSC regression dataset.
    per_pair_metadata : list of dict
        Ordered metadata records aligned to ``rg_full`` rows. Writers use these
        records for optional ``diagnostics/pairs/`` detail output.
    """

    rg: pd.DataFrame
    rg_full: pd.DataFrame
    h2_per_trait: pd.DataFrame
    per_pair_metadata: list[dict[str, object]]


class RegressionRunner:
    """Assemble LDSC regression datasets and dispatch estimator kernels."""

    def __init__(
        self,
        global_config: GlobalConfig | None = None,
        regression_config: RegressionConfig | None = None,
    ) -> None:
        """Initialize the runner with shared defaults for regression workflows."""
        self.global_config = global_config or get_global_config()
        self.regression_config = regression_config or RegressionConfig()

    def build_dataset(self, sumstats_table, ldscore_result, config=None, query_columns=None) -> RegressionDataset:
        """Align inputs and apply the selected model's count and variance rules.

        Use ``prepare_inputs`` and ``dataset_from_prepared`` to share alignment
        across separate pathway models. This convenience method forms one fit.
        """
        prepared = self.prepare_inputs(sumstats_table, ldscore_result, config=config)
        return self.dataset_from_prepared(prepared, query_columns=query_columns, config=config)

    def prepare_inputs(
        self,
        sumstats_table: SumstatsTable,
        ldscore_result: LDScoreSource | LDScoreResult,
        config: RegressionConfig | None = None,
    ) -> PreparedRegressionInputs:
        """Align traits and shared baseline LD scores once for a model batch.

        If both inputs carry known ``GlobalConfig`` snapshots, their critical
        settings are checked before merging. Unknown-provenance inputs, such as
        disk-loaded sumstats, are allowed through this compatibility boundary.
        The merge uses the effective key for the resolved mode: ``SNP`` in
        ``rsid``, ``SNP:<allele_set>`` in ``rsid_allele_aware``, ``CHR:POS`` in
        ``chr_pos``, and ``CHR:POS:<allele_set>`` in
        ``chr_pos_allele_aware``. Base modes never inspect alleles for
        identity. ``RegressionConfig.allow_identity_downgrade`` permits
        same-family allele-aware/base mixes to run under the base mode.

        This preparation preserves every baseline column. Model-dependent
        variance and count selection occur in ``dataset_from_prepared``;
        chi-square filtering, weights, and jackknife work occur inside each fit.
        """
        print_global_config_banner(type(self).__name__, self.global_config)
        config = config or self.regression_config
        legacy_drops = pd.DataFrame(columns=LEGACY_SUMSTATS_DROP_COLUMNS)
        if _is_legacy_sumstats(sumstats_table):
            sumstats_table, legacy_drops = _project_legacy_sumstats_to_panel(sumstats_table, ldscore_result)
        weight_column = REGRESSION_LD_SCORE_COLUMN
        ref_ld_columns = list(ldscore_result.baseline_columns)
        identity = _resolve_h2_identity(sumstats_table, ldscore_result, self.global_config, config)
        identifier_mode = identity.effective_mode
        _regression_genome_build(
            sumstats_table,
            ldscore_result,
            identifier_mode,
            context=sumstats_table.source_path or sumstats_table.trait_name or "sumstats",
        )
        if sumstats_table.config_snapshot is not None and ldscore_result.config_snapshot is not None:
            _validate_regression_config_compatibility(
                sumstats_table.config_snapshot,
                ldscore_result.config_snapshot,
                context="SumstatsTable and LDScoreResult",
            )
        ldscore_mode = _ldscore_identity_mode(ldscore_result, self.global_config, fallback_mode=identifier_mode)
        ldscore_frame = _assemble_regression_ldscore_table(
            ldscore_result,
            [],
            snp_identifier=ldscore_mode,
        )
        ldscore_frame["_ldsc_source_row"] = np.arange(len(ldscore_frame), dtype=np.int64)
        sumstats_frame = sumstats_table.data
        dropped_identity_rows = 0
        if identity.downgrade_applied:
            sumstats_frame, sumstats_dropped = _prepare_regression_identity_table(
                sumstats_frame,
                identifier_mode,
                context=sumstats_table.source_path or "sumstats",
                logger=LOGGER,
            )
            ldscore_frame, ldscore_dropped = _prepare_regression_identity_table(
                ldscore_frame,
                identifier_mode,
                context="LD-score table",
                logger=LOGGER,
            )
            dropped_identity_rows = int(len(sumstats_dropped) + len(ldscore_dropped))
            sumstats_mode = _sumstats_identity_mode(sumstats_table, self.global_config, fallback_mode=ldscore_mode)
            LOGGER.warning(
                f"Identity downgrade enabled: LD-score mode {ldscore_mode}, sumstats mode {sumstats_mode}; "
                f"running regression with effective snp_identifier={identifier_mode!r}. "
                f"Dropped {dropped_identity_rows} duplicate effective-key rows before merge."
            )
        if is_allele_aware_mode(identifier_mode):
            sumstats_keyed = _with_effective_identity_key(
                sumstats_frame,
                identifier_mode,
                context=sumstats_table.source_path or "sumstats",
            )
            ldscore_keyed = _with_effective_identity_key(ldscore_frame, identifier_mode, context="LD-score table")
            ldscore_columns = [REGRESSION_IDENTITY_KEY_COLUMN, *ref_ld_columns, weight_column, "_ldsc_source_row"]
            ldscore_payload = ldscore_keyed.loc[:, ldscore_columns].copy()
            if {"A1", "A2"}.issubset(ldscore_keyed.columns):
                ldscore_payload["A1_ld"] = ldscore_keyed["A1"]
                ldscore_payload["A2_ld"] = ldscore_keyed["A2"]
            # LD-score is the left frame so the inner merge inherits its genomic
            # order; the block jackknife needs genomically contiguous rows for valid
            # SEs, and the sumstats order alone (legacy .sumstats.gz is unsorted) does
            # not guarantee that. Values still align by key regardless of side.
            merged = pd.merge(
                ldscore_payload.reset_index(drop=True),
                sumstats_keyed,
                how="inner",
                on=REGRESSION_IDENTITY_KEY_COLUMN,
                sort=False,
            )
            if {"A1", "A2", "A1_ld", "A2_ld"}.issubset(merged.columns):
                merged = _orient_sumstats_z_to_reference_alleles(merged)
        elif identity_mode_family(identifier_mode) == "rsid":
            merged = pd.merge(
                ldscore_frame.loc[:, ["SNP", *ref_ld_columns, weight_column, "_ldsc_source_row"]].reset_index(drop=True),
                sumstats_frame,
                how="inner",
                on="SNP",
                sort=False,
            )
        else:
            sumstats_keyed = _with_chr_pos_key(sumstats_frame, context=sumstats_table.source_path or "sumstats")
            ldscore_keyed = _with_chr_pos_key(ldscore_frame, context="LD-score table")
            merged = pd.merge(
                ldscore_keyed.loc[:, [CHR_POS_KEY_COLUMN, *ref_ld_columns, weight_column, "_ldsc_source_row"]].reset_index(drop=True),
                sumstats_keyed,
                how="inner",
                on=CHR_POS_KEY_COLUMN,
                sort=False,
            )
        if merged.empty:
            source = sumstats_table.source_path or sumstats_table.trait_name or "sumstats"
            raise LDSCInputError(
                f"h2 regression retained no overlapping {identifier_mode} SNPs after merging sumstats "
                f"'{source}' with {len(ldscore_frame)} LD-score rows. Most likely the sumstats and "
                "LD-score directory use different SNP identifiers or genome builds. Regenerate both inputs "
                "with the same `--snp-identifier` and `--genome-build`, or use `--allow-identity-downgrade` "
                f"for same-family allele-aware/base mixes. Other causes & fixes: {_REGRESSION_NO_OVERLAP_DOC}"
            )

        source_rows = merged.pop("_ldsc_source_row").to_numpy(dtype=np.int64)
        # Keep panel coordinates for delete-block diagnostics even when rsID
        # summary statistics omit coordinates or contain missing placeholders.
        try:
            coordinate_columns = infer_chr_pos_columns(ldscore_result.baseline_table.columns, context="LD-score diagnostics")
        except LDSCInputError:
            coordinate_columns = ()
        for coordinate, source_column in zip(("CHR", "POS"), coordinate_columns):
            merged[f"_ldsc_panel_{coordinate}"] = ldscore_result.baseline_table[source_column].iloc[source_rows].to_numpy()
        dataset = RegressionDataset(
            merged=merged.reset_index(drop=True),
            ref_ld_columns=ref_ld_columns,
            weight_column=weight_column,
            reference_snp_count_totals={},
            count_key_used_for_regression="",
            retained_ld_columns=ref_ld_columns,
            dropped_zero_variance_ld_columns=[],
            trait_names=[name for name in [sumstats_table.trait_name] if name],
            chromosomes_aggregated=[result.chrom for result in ldscore_result.chromosome_results],
            config_snapshot=ldscore_result.config_snapshot,
            effective_snp_identifier=identifier_mode,
            identity_downgrade_applied=identity.downgrade_applied,
            ldscore_overlap=ldscore_result.overlap,
            legacy_sumstats_drops=legacy_drops,
        )
        dataset.validate()
        return PreparedRegressionInputs(dataset, source_rows, ldscore_result, ldscore_mode,
                                        {str(record["column"]): record for record in ldscore_result.count_records})

    def dataset_from_prepared(self, prepared, *, query_columns=None, config=None) -> RegressionDataset:
        """Select one complete genome-wide model and apply its count/variance rules.

        Query identities are validated when this fit selects query columns.
        Shared preparation never decides model-dependent filters or weights.
        """
        selected = list(query_columns or [])
        values = self._read_query_batch(prepared, selected) if selected else None
        return self._dataset_from_query_values(prepared, selected, values, config=config)

    def _read_query_batch(self, prepared, columns):
        """Validate the batch's stored rows once and select the shared trait rows."""
        frame = prepared.ldscore.read_queries(columns)
        assert_same_snp_rows(
            prepared.ldscore.baseline_table, frame,
            context="query rows must match baseline rows on CHR/SNP/POS",
            snp_identifier=prepared.ldscore_identity_mode,
        )
        return frame.loc[:, columns].iloc[prepared.source_rows].reset_index(drop=True)

    def _dataset_from_query_values(self, prepared, selected, query_values, *, config=None):
        """Form one model from already aligned, bounded query values."""
        config = config or self.regression_config
        ref_ld_columns = list(prepared.dataset.ref_ld_columns) + list(selected)
        merged = prepared.dataset.merged
        if selected:
            weight_position = merged.columns.get_loc(prepared.dataset.weight_column)
            merged = pd.concat([merged.iloc[:, :weight_position], query_values.loc[:, selected], merged.iloc[:, weight_position:]], axis=1)
        retained_ld_columns = list(ref_ld_columns)
        dropped_ld_columns: list[str] = []
        if retained_ld_columns:
            variances = merged.loc[:, retained_ld_columns].var()
            dropped_ld_columns = variances.index[variances == 0].tolist()
            retained_ld_columns = [column for column in retained_ld_columns if column not in dropped_ld_columns]
            if retained_ld_columns:
                merged = merged.loc[:, [column for column in merged.columns if column not in dropped_ld_columns]]
                if dropped_ld_columns:
                    LOGGER.warning(
                        "Dropped %d zero-variance LD-score column(s) after merging with sumstats: %s. "
                        "These annotations are constant on the regression SNP set and receive no coefficient.",
                        len(dropped_ld_columns),
                        ", ".join(dropped_ld_columns),
                    )
            else:
                raise LDSCInputError(
                    "h2 regression cannot run because all retained LD-score columns have zero variance "
                    "after merging with sumstats. Most likely the LD-score directory contains a constant "
                    "annotation on the retained SNP set. Use LD scores with at least one varying annotation "
                    "or broaden the regression SNP set."
                )

        count_totals = _count_totals_for_columns(
            [prepared.counts_by_column[column] for column in ref_ld_columns if column in prepared.counts_by_column],
            ref_ld_columns,
        )
        count_key = _select_count_key(count_totals, config.use_common_counts, ref_ld_columns)
        if dropped_ld_columns:
            dropped_index = [ref_ld_columns.index(column) for column in dropped_ld_columns]
            keep_index = [idx for idx in range(len(ref_ld_columns)) if idx not in dropped_index]
            for key, values in list(count_totals.items()):
                count_totals[key] = np.asarray(values)[keep_index]

        dataset = replace(
            prepared.dataset, merged=merged.reset_index(drop=True), ref_ld_columns=ref_ld_columns,
            reference_snp_count_totals=count_totals, count_key_used_for_regression=count_key,
            retained_ld_columns=retained_ld_columns, dropped_zero_variance_ld_columns=dropped_ld_columns,
        )
        dataset.validate()
        return dataset

    def build_rg_dataset(
        self,
        sumstats_table_1: SumstatsTable,
        sumstats_table_2: SumstatsTable,
        ldscore_result: LDScoreSource | LDScoreResult,
        config: RegressionConfig | None = None,
    ) -> RGRegressionDataset:
        """Build the complete genetic-correlation preprocessing dataset.

        This helper merges trait 1, trait 2, and baseline LD scores on the
        resolved SNP identifier, drops missing rows, harmonizes alleles when
        possible, and only then removes zero-variance LD-score columns on the
        final RG SNP set. Allele harmonization exists, but it is skipped if
        allele columns are absent.
        """
        print_global_config_banner(type(self).__name__, self.global_config)
        config = config or self.regression_config
        legacy_drop_frames: list[pd.DataFrame] = []
        if _is_legacy_sumstats(sumstats_table_1):
            sumstats_table_1, drops = _project_legacy_sumstats_to_panel(sumstats_table_1, ldscore_result)
            legacy_drop_frames.append(drops)
        if _is_legacy_sumstats(sumstats_table_2):
            sumstats_table_2, drops = _project_legacy_sumstats_to_panel(sumstats_table_2, ldscore_result)
            legacy_drop_frames.append(drops)
        legacy_drops = (
            pd.concat(legacy_drop_frames, ignore_index=True)
            if legacy_drop_frames
            else pd.DataFrame(columns=LEGACY_SUMSTATS_DROP_COLUMNS)
        )
        identity = _resolve_rg_identity(sumstats_table_1, sumstats_table_2, ldscore_result, self.global_config, config)
        identifier_mode = identity.effective_mode
        for table, label in ((sumstats_table_1, "trait 1"), (sumstats_table_2, "trait 2")):
            _regression_genome_build(
                table,
                ldscore_result,
                identifier_mode,
                context=table.source_path or table.trait_name or label,
            )
            if table.config_snapshot is not None and ldscore_result.config_snapshot is not None:
                _validate_regression_config_compatibility(
                    table.config_snapshot,
                    ldscore_result.config_snapshot,
                    context=f"{label} SumstatsTable and LDScoreResult",
                )

        weight_column = REGRESSION_LD_SCORE_COLUMN
        ref_ld_columns = list(ldscore_result.baseline_columns)
        ldscore_mode = _ldscore_identity_mode(ldscore_result, self.global_config, fallback_mode=identifier_mode)
        trait_1_mode = _sumstats_identity_mode(sumstats_table_1, self.global_config, fallback_mode=ldscore_mode)
        trait_2_mode = _sumstats_identity_mode(sumstats_table_2, self.global_config, fallback_mode=ldscore_mode)
        ldscore_frame = _assemble_regression_ldscore_table(ldscore_result, [], snp_identifier=ldscore_mode)
        left_frame = sumstats_table_1.data
        right_frame = sumstats_table_2.data
        if identity.downgrade_applied:
            left_frame, left_dropped = _prepare_regression_identity_table(
                left_frame,
                identifier_mode,
                context=sumstats_table_1.source_path or "trait 1 sumstats",
                logger=LOGGER,
            )
            right_frame, right_dropped = _prepare_regression_identity_table(
                right_frame,
                identifier_mode,
                context=sumstats_table_2.source_path or "trait 2 sumstats",
                logger=LOGGER,
            )
            ldscore_frame, ldscore_dropped = _prepare_regression_identity_table(
                ldscore_frame,
                identifier_mode,
                context="LD-score table",
                logger=LOGGER,
            )
            dropped_identity_rows = int(len(left_dropped) + len(right_dropped) + len(ldscore_dropped))
            LOGGER.warning(
                f"Identity downgrade enabled: LD-score mode {ldscore_mode}, trait 1 mode {trait_1_mode}, "
                f"trait 2 mode {trait_2_mode}; running rg with effective snp_identifier={identifier_mode!r}. "
                f"Dropped {dropped_identity_rows} duplicate effective-key rows before merge."
            )

        left = left_frame.rename(columns={"N": "N1", "Z": "Z1"})
        left = left.rename(columns={"FRQ": "FRQ1"})
        right = right_frame.rename(
            columns={"A1": "A1x", "A2": "A2x", "N": "N2", "Z": "Z2", "FRQ": "FRQ2"}
        )
        right_payload = [column for column in ["A1x", "A2x", "N2", "Z2", "FRQ2"] if column in right.columns]
        if is_allele_aware_mode(identifier_mode):
            left_keyed = _with_effective_identity_key(left_frame, identifier_mode, context=sumstats_table_1.source_path or "sumstats")
            left = left_keyed.rename(columns={"N": "N1", "Z": "Z1", "FRQ": "FRQ1"})
            right_keyed = _with_effective_identity_key(right_frame, identifier_mode, context=sumstats_table_2.source_path or "sumstats")
            right = right_keyed.rename(
                columns={"A1": "A1x", "A2": "A2x", "N": "N2", "Z": "Z2", "FRQ": "FRQ2"}
            )
            right_payload = [column for column in ["A1x", "A2x", "N2", "Z2", "FRQ2"] if column in right.columns]
            ldscore_keyed = _with_effective_identity_key(ldscore_frame, identifier_mode, context="LD-score table")
            # LD-score leads the first merge so the result (and the final merge below)
            # inherits its genomic order for valid jackknife block SEs. See build_dataset.
            left_with_ld = pd.merge(
                ldscore_keyed.loc[:, [REGRESSION_IDENTITY_KEY_COLUMN, *ref_ld_columns, weight_column]].reset_index(drop=True),
                left,
                how="inner",
                on=REGRESSION_IDENTITY_KEY_COLUMN,
                sort=False,
            )
            merged = pd.merge(
                left_with_ld,
                right.loc[:, [REGRESSION_IDENTITY_KEY_COLUMN, *right_payload]],
                how="inner",
                on=REGRESSION_IDENTITY_KEY_COLUMN,
                sort=False,
            )
        elif identity_mode_family(identifier_mode) == "rsid":
            left_with_ld = pd.merge(
                ldscore_frame.loc[:, ["SNP", *ref_ld_columns, weight_column]].reset_index(drop=True),
                left,
                how="inner",
                on="SNP",
                sort=False,
            )
            merged = pd.merge(
                left_with_ld,
                right.loc[:, ["SNP", *right_payload]],
                how="inner",
                on="SNP",
                sort=False,
            )
        else:
            left_keyed = _with_chr_pos_key(left, context=sumstats_table_1.source_path or "sumstats")
            right_keyed = _with_chr_pos_key(right, context=sumstats_table_2.source_path or "sumstats")
            ldscore_keyed = _with_chr_pos_key(ldscore_frame, context="LD-score table")
            left_with_ld = pd.merge(
                ldscore_keyed.loc[:, [CHR_POS_KEY_COLUMN, *ref_ld_columns, weight_column]].reset_index(drop=True),
                left_keyed,
                how="inner",
                on=CHR_POS_KEY_COLUMN,
                sort=False,
            )
            merged = pd.merge(
                left_with_ld,
                right_keyed.loc[:, [CHR_POS_KEY_COLUMN, *right_payload]],
                how="inner",
                on=CHR_POS_KEY_COLUMN,
                sort=False,
            )
        required_numeric = ["Z1", "N1", "Z2", "N2", weight_column, *ref_ld_columns]
        merged = merged.dropna(subset=[column for column in required_numeric if column in merged.columns]).reset_index(
            drop=True
        )
        if merged.empty:
            raise LDSCInputError(
                f"rg regression retained no overlapping {identifier_mode} SNPs after merging both sumstats "
                f"tables with {len(ldscore_frame)} LD-score rows. Most likely one trait was munged with a "
                "different SNP identifier mode, genome build, or allele convention than the LD-score directory. "
                "Regenerate all inputs with matching identity settings, or use `--allow-identity-downgrade` "
                f"for same-family allele-aware/base mixes. Other causes & fixes: {_REGRESSION_NO_OVERLAP_DOC}"
            )

        if {"A1", "A2", "A1x", "A2x"}.issubset(merged.columns):
            alleles = merged["A1"] + merged["A2"] + merged["A1x"] + merged["A2x"]
            keep = reg._filter_alleles(alleles)
            kept_alleles = alleles.loc[keep].reset_index(drop=True)
            merged = merged.loc[keep].reset_index(drop=True)
            if merged.empty:
                raise LDSCInputError(
                    "rg regression retained no allele-compatible SNPs after harmonizing the two sumstats tables. "
                    "Most likely the traits use different effect-allele conventions or one input was not allele "
                    "harmonized during munging. Regenerate both munged sumstats with A1/A2 columns and matching "
                    "SNP identifier mode."
                )
            merged["Z2"] = reg._align_alleles(merged["Z2"].copy(), kept_alleles)
            if "FRQ2" in merged.columns:
                flip_index = kept_alleles.map(reg.FLIP_ALLELES).astype(bool)
                valid_frq = pd.to_numeric(merged["FRQ2"], errors="coerce").between(0.0, 1.0, inclusive="both")
                transform = flip_index & valid_frq
                merged.loc[transform, "FRQ2"] = 1.0 - pd.to_numeric(
                    merged.loc[transform, "FRQ2"], errors="coerce"
                )

        retained_ld_columns = list(ref_ld_columns)
        dropped_ld_columns: list[str] = []
        if retained_ld_columns:
            variances = merged.loc[:, retained_ld_columns].var()
            dropped_ld_columns = variances.index[variances == 0].tolist()
            retained_ld_columns = [column for column in retained_ld_columns if column not in dropped_ld_columns]
            if retained_ld_columns:
                merged = merged.loc[:, [column for column in merged.columns if column not in dropped_ld_columns]]
                if dropped_ld_columns:
                    LOGGER.warning(
                        "Dropped %d zero-variance LD-score column(s) after merging with sumstats: %s. "
                        "These annotations are constant on the regression SNP set and receive no coefficient.",
                        len(dropped_ld_columns),
                        ", ".join(dropped_ld_columns),
                    )
            else:
                raise LDSCInputError(
                    "rg regression cannot run because all retained LD-score columns have zero variance "
                    "after merging the two traits with LD scores. Most likely the LD-score directory contains "
                    "constant annotations on the retained SNP set. Use LD scores with at least one varying "
                    "annotation or broaden the regression SNP set."
                )

        count_totals = _count_totals_for_columns(ldscore_result.count_records, ref_ld_columns)
        if dropped_ld_columns:
            dropped_index = [ref_ld_columns.index(column) for column in dropped_ld_columns]
            keep_index = [idx for idx in range(len(ref_ld_columns)) if idx not in dropped_index]
            for key, values in list(count_totals.items()):
                count_totals[key] = np.asarray(values)[keep_index]
        count_key = _select_count_key(count_totals, config.use_common_counts, ref_ld_columns)

        dataset = RGRegressionDataset(
            merged=merged.reset_index(drop=True),
            ref_ld_columns=ref_ld_columns,
            weight_column=weight_column,
            reference_snp_count_totals=count_totals,
            count_key_used_for_regression=count_key,
            retained_ld_columns=retained_ld_columns,
            dropped_zero_variance_ld_columns=dropped_ld_columns,
            trait_names=[name for name in [sumstats_table_1.trait_name, sumstats_table_2.trait_name] if name],
            chromosomes_aggregated=[result.chrom for result in ldscore_result.chromosome_results],
            config_snapshot=ldscore_result.config_snapshot,
            effective_snp_identifier=identifier_mode,
            identity_downgrade_applied=identity.downgrade_applied,
            legacy_sumstats_drops=legacy_drops,
        )
        dataset.validate()
        return dataset

    def estimate_h2(
        self,
        dataset: RegressionDataset,
        config: RegressionConfig | None = None,
    ):
        """Estimate single-trait heritability from a prepared dataset.

        Partitioned (multi-annotation) datasets apply the legacy default
        chi-square cap ``max(0.001 * N.max(), 80)`` when ``chisq_max`` is unset,
        down-weighting outlier SNPs. Single-annotation datasets stay uncapped and
        default the two-step cutoff to ``30`` instead. See
        :func:`_resolve_default_chisq_max`.
        """
        return self._fit_h2_dataset(dataset, config=config).estimator

    def _fit_h2_dataset(
        self,
        dataset: RegressionDataset,
        config: RegressionConfig | None = None,
    ) -> _FitOutcome:
        """Fit h2 once and retain its exact SNP population and diagnostic bins."""
        config = config or self.regression_config
        merged = dataset.merged
        n_snp = len(merged)
        n_blocks = min(n_snp, config.n_blocks)
        x = np.asarray(merged[dataset.retained_ld_columns])
        chisq = np.asarray(merged["Z"] ** 2).reshape((n_snp, 1))
        # N.max() is read pre-filter so the partitioned default cap matches legacy.
        chisq_max = _resolve_default_chisq_max(
            config.chisq_max, len(dataset.retained_ld_columns), merged["N"].max()
        )
        if chisq_max is not None:
            keep = np.ravel(chisq <= chisq_max)
            n_removed = n_snp - int(keep.sum())
            merged = merged.loc[keep].reset_index(drop=True)
            n_snp = len(merged)
            n_blocks = min(n_snp, config.n_blocks)
            x = np.asarray(merged[dataset.retained_ld_columns])
            chisq = np.asarray(merged["Z"] ** 2).reshape((n_snp, 1))
            if n_removed:
                LOGGER.info(
                    "Removed %d SNPs with chi^2 > %g (%d SNPs remain).",
                    n_removed,
                    chisq_max,
                    n_snp,
                )
        intercept = None
        if not config.use_intercept:
            intercept = 1
        elif config.intercept_h2 is not None and not isinstance(config.intercept_h2, list):
            intercept = float(config.intercept_h2)
        two_step = config.two_step_cutoff
        if two_step is None and intercept is None and len(dataset.retained_ld_columns) == 1:
            two_step = 30
        old_weights = len(dataset.retained_ld_columns) > 1
        _raise_on_model_collinearity(dataset, x)
        reference_snp_counts = np.asarray(
            dataset.reference_snp_count_totals[dataset.count_key_used_for_regression]
        ).reshape((1, -1))
        try:
            hsq = reg.Hsq(
                chisq,
                x,
                np.asarray(merged[[dataset.weight_column]]),
                np.asarray(merged[["N"]]),
                reference_snp_counts,
                n_blocks=n_blocks,
                intercept=intercept,
                twostep=two_step,
                old_weights=old_weights,
            )
        except JackknifeIdentifiabilityError as error:
            # Two-step jackknives may use a different row population. The
            # partitioned model uses these exact post-filter rows throughout.
            if two_step is None:
                _add_jackknife_model_context(error, merged, dataset.retained_ld_columns, intercept)
            raise
        diagnostic_bins = None
        if len(dataset.retained_ld_columns) == 1:
            diagnostic_bins = summarize_ld_score_regression_bins(
                hsq,
                ld_score=x[:, 0],
                chi_square=chisq,
                sample_size=np.asarray(merged[["N"]]),
                regression_ld_score=np.asarray(merged[[dataset.weight_column]]),
                reference_snp_count=float(reference_snp_counts[0, 0]),
            )
        return _FitOutcome(
            estimator=hsq, dataset=replace(dataset, merged=merged), n_blocks=n_blocks,
            effective_chisq_max=chisq_max, diagnostic_bins=diagnostic_bins,
        )

    def estimate_partitioned_h2(
        self,
        sumstats_table: SumstatsTable,
        ldscore_result: LDScoreSource | LDScoreResult,
        *,
        query_column: str,
        config: RegressionConfig | None = None,
    ) -> pd.DataFrame:
        """Estimate partitioned heritability for one query annotation.

        Parameters
        ----------
        sumstats_table : SumstatsTable
            Munged single-trait summary statistics.
        ldscore_result : LDScoreSource or LDScoreResult
            Canonical LD-score result containing baseline LD-score columns and
            at least one query LD-score column.
        query_column : str
            Query annotation to test in the baseline-plus-query model. The
            column must be present in ``ldscore_result.query_columns``.
        config : RegressionConfig, optional
            Regression settings. Defaults to the runner's config.

        Returns
        -------
        pandas.DataFrame
            One-row compact partitioned-h2 summary for ``query_column``.

        Raises
        ------
        ValueError
            If the LD-score result has no query annotations or ``query_column``
            is not one of the available query annotations.
        """
        selected_queries = _validate_partitioned_query_columns(ldscore_result, [query_column])
        dataset = self.build_dataset(sumstats_table, ldscore_result, config=config, query_columns=selected_queries)
        hsq = self.estimate_h2(dataset, config=config)
        samp_prev, pop_prev = _config_prevalence(config or self.regression_config)
        return summarize_partitioned_h2(hsq, dataset, selected_queries, samp_prev=samp_prev, pop_prev=pop_prev)

    def _fit_partitioned_query(self, prepared, query, values, directory, *, config):
        """Fit, summarize and privately stage one complete genome-wide query.

        Both inline and worker execution use this boundary. Only model errors
        become query outcomes; staging errors and interrupts still abort runs.
        Returned objects contain no estimator, SNP matrix or traceback frames.
        """
        stage = "model_preparation"
        try:
            with _log_phase_timing("query model preparation"):
                dataset = self._dataset_from_query_values(prepared, [query], values, config=config)
            stage = "estimator"
            with _log_phase_timing("estimator execution"):
                outcome = self._fit_h2_dataset(dataset, config=config)
            hsq, dataset = outcome.estimator, outcome.dataset
            stage = "summary"
            with _log_phase_timing("query summary"):
                samp_prev, pop_prev = _config_prevalence(config)
                focal = summarize_partitioned_h2(hsq, dataset, [query], samp_prev=samp_prev, pop_prev=pop_prev)
                full = summarize_partitioned_h2(hsq, dataset, dataset.retained_ld_columns,
                                               samp_prev=samp_prev, pop_prev=pop_prev)
                delete_values = _coefficient_delete_frame(hsq, dataset.retained_ld_columns)
        except Exception as error:
            status = "unestimable" if isinstance(error, JackknifeIdentifiabilityError) else "failed"
            record = {"query_annotation": query, "status": status, "stage": stage,
                      "error_type": type(error).__name__, "error_message": str(error)}
            LOGGER.exception("Query %r failed at %s; marked %s: %s", query, stage, status, error)
            if isinstance(error, JackknifeIdentifiabilityError):
                LOGGER.error("Query %r jackknife diagnostics: %s", query, json.dumps(error.failures))
            error.add_note(f"partitioned-h2 query {query!r} failed at {stage}.")
            error.__traceback__ = error.__cause__ = error.__context__ = None
            return None, None, record, error
        with _log_phase_timing("query output staging"):
            artifacts = stage_partitioned_h2_fit(
                directory, full, delete_values,
                {"dropped_zero_variance_ld_columns": list(dataset.dropped_zero_variance_ld_columns),
                 "n_snps": outcome.n_snps, "n_blocks_used": outcome.n_blocks,
                 "effective_chisq_max": outcome.effective_chisq_max,
                 "effective_snp_identifier": dataset.effective_snp_identifier,
                 "identity_downgrade_applied": dataset.identity_downgrade_applied},
            )
        return focal, artifacts, {"query_annotation": query, "status": "success", "stage": "complete",
                                  "error_type": "", "error_message": ""}, None

    def estimate_partitioned_h2_batch(
        self, sumstats_table, ldscore_result, *, output_dir, query_columns=None,
        config=None, query_batch_size=1000, overwrite=False, summary_sort_by="auto", metadata=None,
        continue_on_query_error=False, threads=1,
    ) -> PartitionedH2BatchResult:
        """Fit and write separate pathway models using shared SNP preparation.

        Parameters
        ----------
        sumstats_table : SumstatsTable
            Munged trait statistics, aligned once to shared baseline LD scores.
        ldscore_result : LDScoreSource or LDScoreResult
            Canonical shared baseline, metadata, and explicit query reads.
        output_dir : path-like
            Required destination for canonical results and owned private staging.
            User and environment tokens are expanded before either is created.
        query_columns : sequence of str or None, optional
            Focal models to fit; None selects all declared queries. An empty
            sequence fits the complete baseline as one functional model.
        config : RegressionConfig, optional
            Statistical settings. Model-dependent filtering and weights are
            recomputed for every fit over the complete aligned genome.
        query_batch_size : int, optional
            Maximum query columns loaded at once, default 1000. Pathways in a
            batch are still fitted separately against baseline categories.
            Also caps concurrent query workers; a width of one runs inline.
            Independent of the batch width used to generate the LD scores.
        threads : int, optional
            Query worker processes, default 1 (inline). Positive values request
            that many workers; -1 uses available CPUs, -2 leaves one CPU free.
            Uses the same validation and resolution as ldscore and gene-index
            construction: negative values use CPU affinity (machine count
            fallback), and all requests are capped by query count and batch
            size. Zero, booleans and non-integers are invalid. Positive requests
            must fit the caller's CPU allocation; SLURM_CPUS_PER_TASK is not
            read separately. Parallel workers use one native numerical
            thread each; inline execution preserves the caller's thread policy.
            On macOS before 15, launch Python with VECLIB_MAXIMUM_THREADS=1
            for parallel execution; newer macOS uses the Accelerate runtime API.
        overwrite : bool, optional
            Replace workflow-owned outputs, default False.
        summary_sort_by : str, optional
            Final summary ordering; auto uses coefficient-p for focal models.
        metadata : dict, optional
            Additional diagnostic provenance copied into each model's metadata.
        continue_on_query_error : bool, optional
            Default False: any query regression error prevents publication.
            If True, catch ordinary exceptions during per-query preparation,
            fitting, or summarization, record their status and traceback, and
            publish successful fits. Shared input and output errors, interrupts,
            and runs without any successful fit still fail.

        Returns
        -------
        PartitionedH2BatchResult
            Summary and persistent detail paths after successful publication.
            ``query_status`` records every attempted query in input order;
            scientific summaries and model artifacts include successful fits
            only. Strict runs publish diagnostics alone when any query fails.

        Notes
        -----
        Parallel execution uses spawned processes and read-only numeric maps
        inside private output scratch. Each worker owns one complete model;
        batch maps are removed between batches and workers exit before shared
        scratch cleanup. Python scripts using multiple workers must launch
        this call under ``if __name__ == "__main__":``. Worker death and transport
        failures abort the run and identify the affected in-flight queries.
        Query statuses follow requested query order; stable summary sorting
        determines the published row and folder order regardless of completion
        order. Numeric maps reduce transport copies, but each concurrent model
        still needs its own filtered matrices and estimator workspace.

        ``continue_on_query_error`` is the Python equivalent of the CLI flag
        ``--continue-on-query-error``. Omit the flag or leave this argument
        False for strict publication; enable it to skip failed query fits.
        It does not change gene-list resolution, SNP filtering, weights, block
        boundaries, or the numerical solver. Failed queries have no inferred
        coefficients, standard errors, p-values, or result folders.

        Errors are sent to the configured LDSC logger. The CLI installs the
        ``diagnostics/partitioned-h2.log`` file handler; this direct Python
        method does not install a log file handler. The status TSV is written
        in either interface after the query scan, including strict failures.
        Shared loading/publication errors and interrupts can stop the scan
        before that final diagnostic table is written.
        """
        _validate_threads(threads)
        if isinstance(query_batch_size, bool) or not isinstance(query_batch_size, int) or query_batch_size < 1:
            raise ValueError("query_batch_size must be a positive integer.")
        if not isinstance(continue_on_query_error, bool):
            raise ValueError("continue_on_query_error must be a boolean.")
        config = config or self.regression_config
        if ldscore_result.overlap is None:
            raise LDSCInputError(PARTITIONED_H2_REQUIRES_OVERLAP_MESSAGE)
        queries = list(ldscore_result.query_columns if query_columns is None else query_columns)
        queries = _validate_partitioned_query_columns(ldscore_result, queries)
        worker_count = _resolve_worker_count(threads, min(len(queries), query_batch_size))
        sort_by = _resolve_summary_sort(summary_sort_by, has_queries=bool(queries))
        output_dir = normalize_path_token(output_dir)
        writer = PartitionedH2DirectoryWriter()
        writer.artifact_family(output_dir, write_per_query_results=bool(queries),
                               coefficient_delete_values=not queries, query_status=bool(queries)).preflight(overwrite=overwrite)
        samp_prev, pop_prev = _config_prevalence(config)
        metadata = dict(metadata or {})
        defaults = {
            "trait_name": sumstats_table.trait_name,
            "count_kind": "common" if config.use_common_counts else "all",
            "analysis_type": "cell_type_specific" if queries else "functional_category",
            "headline_metric": "coefficient" if queries else "enrichment",
            "enrichment_p_test": "two_sided_t", "coefficient_p_test": "one_sided_greater",
            "annotation_types": dict(ldscore_result.annotation_types),
        }
        source_metadata = ldscore_result.output_paths.get("metadata")
        if source_metadata:
            defaults["ldscore_dir"] = str(Path(source_metadata).resolve().parent)
        for name, value in defaults.items():
            metadata.setdefault(name, value)
        with AnnotationWorkspace(output_dir) as workspace:
            with _log_phase_timing("shared regression preparation"):
                prepared = self.prepare_inputs(sumstats_table, ldscore_result, config=config)
            if not queries:
                dataset = self.dataset_from_prepared(prepared, config=config)
                if len(dataset.retained_ld_columns) == 1:
                    LOGGER.warning(
                        "partitioned-h2 functional regime retained a single annotation (%s); the fit is "
                        "degenerate and collapses to single-annotation h2 behavior (two-step estimator, no "
                        "default chi-square cap). Most likely the LD-score directory has only one baseline "
                        "annotation. Provide multiple baseline annotations for a meaningful partitioned analysis.",
                        dataset.retained_ld_columns[0],
                    )
                outcome = self._fit_h2_dataset(dataset, config=config)
                summary = summarize_partitioned_h2(outcome.estimator, outcome.dataset, outcome.dataset.retained_ld_columns,
                                                   samp_prev=samp_prev, pop_prev=pop_prev)
                summary = _sort_partitioned_h2_summary(summary, sort_by)
                with _log_phase_timing("output writing"):
                    paths = writer.write(
                        summary, PartitionedH2OutputConfig(output_dir=output_dir, overwrite=overwrite),
                        metadata={**metadata, "n_snps": outcome.n_snps, "effective_chisq_max": outcome.effective_chisq_max,
                                  "n_blocks_used": outcome.n_blocks},
                        coefficient_delete_values=_coefficient_delete_frame(outcome.estimator, outcome.dataset.retained_ld_columns),
                    )
                return PartitionedH2BatchResult(summary, paths)

            rows, private_artifacts, statuses = [], {}, []
            first_error = None
            LOGGER.info("Query workers: requested=%d, effective=%d; native threads=%s; query batch size=%d.",
                        threads, worker_count, "1 per worker" if worker_count > 1 else "caller settings", query_batch_size)
            metadata.update(query_workers_requested=threads, query_workers_effective=worker_count,
                            query_worker_native_threads=1 if worker_count > 1 else None)
            models = _PartitionedModelInputs(prepared.dataset, prepared.counts_by_column)
            if worker_count > 1:
                from ._partitioned_h2_parallel import _QueryWorkers, _stage_query_batch
                with _log_phase_timing("query worker preparation"):
                    workers = _QueryWorkers(workspace.path, models, self, config, worker_count)
            else:
                workers = nullcontext(None)
            with workers as pool:
                for start in range(0, len(queries), query_batch_size):
                    batch_columns = queries[start:start + query_batch_size]
                    with _log_phase_timing("query batch loading"):
                        batch = self._read_query_batch(prepared, batch_columns)
                    if pool is None:
                        outcomes = (self._fit_partitioned_query(
                            models, query, batch, workspace.path / f"fit-{start + index}", config=config,
                        ) for index, query in enumerate(batch_columns))
                    else:
                        with TemporaryDirectory(prefix="query-batch-", dir=workspace.path) as batch_root:
                            with _log_phase_timing("query batch mapping"):
                                mapped = _stage_query_batch(batch, Path(batch_root) / "values")
                            batch = None
                            tasks = [(query, mapped[query], workspace.path / f"fit-{start + index}")
                                     for index, query in enumerate(batch_columns)]
                            with _log_phase_timing("query batch execution"):
                                outcomes = pool.map(tasks)
                    for focal, artifacts, status, error in outcomes:
                        statuses.append(status)
                        if error is not None:
                            first_error = first_error or error
                        else:
                            rows.append(focal)
                            private_artifacts[status["query_annotation"]] = artifacts
                    del outcomes, batch
            del models
            del prepared
            query_status = pd.DataFrame(statuses)
            n_failed = len(queries) - len(rows)
            metadata.update(query_error_policy="continue" if continue_on_query_error else "raise",
                            n_queries_requested=len(queries), n_queries_successful=len(rows), n_queries_failed=n_failed)
            if first_error is not None:
                LOGGER.warning("Partitioned-h2 query scan: %d successful, %d failed. See diagnostics/query_status.tsv "
                               "and diagnostics/partitioned-h2.log.", len(rows), n_failed)
                if not continue_on_query_error or not rows:
                    writer.write_query_status(output_dir, query_status)
                    first_error.add_note(
                        f"{n_failed}/{len(queries)} queries failed; no results published. "
                        + ("No query fit succeeded." if not rows else
                           "Use --continue-on-query-error to publish successful query fits.")
                    )
                    raise first_error
            summary = _sort_partitioned_h2_summary(pd.concat(rows, ignore_index=True), sort_by)
            with _log_phase_timing("output writing"):
                paths = writer.write(summary, PartitionedH2OutputConfig(output_dir=output_dir, overwrite=overwrite,
                                                                       write_per_query_results=True),
                                     per_query_artifacts=private_artifacts, metadata=metadata, query_status=query_status)
            root = Path(paths["per_query_root"])
            artifacts = {
                str(record["query_annotation"]): PartitionedH2FitArtifacts(
                    root / str(record["folder"]) / "partitioned_h2_full.tsv",
                    root / str(record["folder"]) / "coefficient_delete_values.parquet",
                    root / str(record["folder"]) / "metadata.json",
                ) for record in writer._query_records(summary)
            }
            return PartitionedH2BatchResult(summary, paths, artifacts, query_status)

    def estimate_rg(
        self,
        sumstats_table_1: SumstatsTable,
        sumstats_table_2: SumstatsTable,
        ldscore_result: LDScoreSource | LDScoreResult,
        config: RegressionConfig | None = None,
    ):
        """Estimate genetic correlation between two munged summary-stat tables.

        RG preprocessing is delegated to :meth:`build_rg_dataset`, which keeps
        the final trait1/trait2/LD-score SNP set explicit before arrays are
        passed to the kernel. Allele harmonization exists, but it is skipped if
        allele columns are absent.
        """
        config = config or self.regression_config
        dataset = self.build_rg_dataset(sumstats_table_1, sumstats_table_2, ldscore_result, config=config)
        return self._fit_rg_dataset(dataset, config=config).estimator

    def _fit_rg_dataset(
        self,
        dataset: RGRegressionDataset,
        config: RegressionConfig | None = None,
    ) -> _FitOutcome:
        """Fit rg and retain the population selected by its explicit product cap."""
        config = config or self.regression_config
        merged = dataset.merged
        n_snp = len(merged)
        if config.chisq_max is not None:
            # Legacy rg filter (sumstats.py `_rg`): drop SNPs with
            # Z1^2 * Z2^2 > chisq_max^2, kept inclusive (<=) per the -max convention.
            keep = np.ravel((merged["Z1"] ** 2) * (merged["Z2"] ** 2) <= config.chisq_max ** 2)
            merged = merged.loc[keep].reset_index(drop=True)
            n_snp = len(merged)
        n_blocks = min(n_snp, config.n_blocks)
        intercept_hsq = _select_intercept(config.intercept_h2, use_intercept=config.use_intercept, default_when_disabled=1)
        intercept_gencov = _select_intercept(
            config.intercept_gencov,
            use_intercept=config.use_intercept,
            default_when_disabled=0,
        )
        # Legacy estimate_rg defaults two-step to 30 for a single-annotation fit
        # with a free h2 intercept and no explicit cutoff (sumstats.py:400).
        two_step = config.two_step_cutoff
        if two_step is None and intercept_hsq is None and len(dataset.retained_ld_columns) == 1:
            two_step = 30
        fitted = reg.RG(
            np.asarray(merged[["Z1"]]),
            np.asarray(merged[["Z2"]]),
            np.asarray(merged[dataset.retained_ld_columns]),
            np.asarray(merged[[dataset.weight_column]]),
            np.asarray(merged[["N1"]]),
            np.asarray(merged[["N2"]]),
            np.asarray(dataset.reference_snp_count_totals[dataset.count_key_used_for_regression]).reshape((1, -1)),
            intercept_hsq1=intercept_hsq,
            intercept_hsq2=intercept_hsq,
            intercept_gencov=intercept_gencov,
            n_blocks=n_blocks,
            twostep=two_step,
        )

        return _FitOutcome(
            estimator=fitted, dataset=replace(dataset, merged=merged), n_blocks=n_blocks,
            effective_chisq_max=config.chisq_max,
        )

    def estimate_rg_pairs(
        self,
        sumstats_tables: Sequence[SumstatsTable],
        ldscore_result: LDScoreSource | LDScoreResult,
        *,
        anchor_index: int | None = None,
        config: RegressionConfig | None = None,
        prevalences: Sequence[tuple[float | None, float | None]] | None = None,
    ) -> RgResultFamily:
        """Estimate genetic correlations for all requested trait pairs.

        With no anchor, all unordered pairs are fit in input order. With an
        anchor, only anchor-vs-rest pairs are fit, preserving the original input
        order of non-anchor traits. Pair-level failures are recorded as NaN rows
        with ``status='failed'`` and do not stop later pairs.

        Parameters
        ----------
        sumstats_tables : sequence of SumstatsTable
            Munged summary-statistic tables for two or more traits. Trait names
            should already be unique; the CLI wrapper performs filename-based
            disambiguation before calling this method.
        ldscore_result : LDScoreSource or LDScoreResult
            Canonical LD-score result. Only baseline LD scores are used for rg.
        anchor_index : int or None, optional
            If provided, index of the anchor trait. The method computes
            ``anchor`` against every other trait and skips non-anchor pairs.
            Default is ``None``, which computes all unordered pairs.
        config : RegressionConfig or None, optional
            Regression settings. Defaults to the runner's config.

        Returns
        -------
        RgResultFamily
            Concise rg table, full diagnostic table, per-trait h2 table, and
            ordered per-pair metadata.

        Raises
        ------
        ValueError
            If fewer than two traits are supplied or ``anchor_index`` is out of
            bounds. Whole-run input failures, such as incompatible h2 datasets,
            also propagate as ordinary exceptions.
        """
        config = config or self.regression_config
        tables = list(sumstats_tables)
        if len(tables) < 2:
            raise LDSCUserError(
                "rg regression requires at least two sumstats inputs. Most likely `--sumstats-sources` "
                "resolved to fewer than two files. Pass two or more munged `.sumstats` artifacts."
            )
        if anchor_index is not None and not 0 <= anchor_index < len(tables):
            raise LDSCUsageError(
                f"rg regression received anchor_index={anchor_index}, but valid indices are 0 through {len(tables) - 1}. "
                "Most likely the anchor trait was resolved against a different input list. "
                "Pass an anchor trait name or path that belongs to the current `--sumstats-sources` list."
            )

        prev = list(prevalences) if prevalences is not None else [(None, None)] * len(tables)

        h2_rows = []
        for table, (samp_prev, pop_prev) in zip(tables, prev):
            h2_dataset = self.build_dataset(table, ldscore_result, config=config)
            outcome = self._fit_h2_dataset(h2_dataset, config=config)
            hsq, h2_dataset = outcome.estimator, outcome.dataset
            n_snps_used = outcome.n_snps
            h2_rows.append(
                summarize_total_h2(
                    hsq, h2_dataset, trait_name=_trait_label(table), n_snps_used=n_snps_used,
                    samp_prev=samp_prev, pop_prev=pop_prev,
                )
            )
        h2_per_trait = pd.concat(h2_rows, axis=0, ignore_index=True) if h2_rows else pd.DataFrame()

        rg_full_rows: list[dict[str, object]] = []
        rg_rows: list[dict[str, object]] = []
        per_pair_metadata: list[dict[str, object]] = []
        pair_kind = "anchor" if anchor_index is not None else "all_pairs"
        for i, j in _iter_rg_pairs(len(tables), anchor_index):
            trait_1 = _trait_label(tables[i])
            trait_2 = _trait_label(tables[j])
            try:
                dataset = self.build_rg_dataset(tables[i], tables[j], ldscore_result, config=config)
                outcome = self._fit_rg_dataset(dataset, config=config)
                fitted, dataset = outcome.estimator, outcome.dataset
                n_snps_used = outcome.n_snps
                full_row = _summarize_rg_pair(
                    fitted, dataset, trait_1=trait_1, trait_2=trait_2, pair_kind=pair_kind,
                    n_snps_used=n_snps_used,
                    samp_prev_1=prev[i][0], pop_prev_1=prev[i][1],
                    samp_prev_2=prev[j][0], pop_prev_2=prev[j][1],
                )
                metadata = _rg_pair_metadata(tables[i], tables[j], outcome, full_row, config, pair_kind)
            except Exception as exc:
                error = _format_exception(exc)
                LOGGER.warning(f"Failed rg for pair '{trait_1}' vs '{trait_2}': {error}", exc_info=True)
                full_row = _failed_rg_full_row(trait_1=trait_1, trait_2=trait_2, pair_kind=pair_kind, error=error)
                metadata = _failed_rg_pair_metadata(tables[i], tables[j], full_row, config, pair_kind)
            rg_full_rows.append(full_row)
            rg_rows.append(_concise_rg_row(full_row))
            per_pair_metadata.append(metadata)

        rg_full = pd.DataFrame(rg_full_rows, columns=RG_FULL_COLUMNS)
        rg = pd.DataFrame(rg_rows, columns=RG_CONCISE_COLUMNS)
        return RgResultFamily(rg=rg, rg_full=rg_full, h2_per_trait=h2_per_trait, per_pair_metadata=per_pair_metadata)


def _validate_regression_config_compatibility(
    a: GlobalConfig,
    b: GlobalConfig,
    context: str = "",
) -> None:
    """Validate regression provenance fields outside SNP identity mode."""
    prefix = f" when combining {context}" if context else ""
    if a.genome_build != b.genome_build:
        raise ConfigMismatchError(
            f"genome_build mismatch{prefix}: {a.genome_build!r} vs "
            f"{b.genome_build!r}. These objects were computed under different "
            "genome-build assumptions and cannot be safely merged."
        )


def _regression_genome_build(
    sumstats_table: SumstatsTable,
    ldscore_result: LDScoreSource | LDScoreResult,
    identifier_mode: str,
    *,
    context: str,
) -> None:
    """Verify a chr_pos sumstats shares the LD-score panel's genome build.

    rsID-family runs skip the check entirely because their merge identity is
    independent of coordinate build. For chr_pos-family runs the sumstats build
    is taken from its footer metadata when present, otherwise inferred from its
    coordinates. A build that disagrees with the panel, or coordinates that
    cannot be confidently dated, is a hard error directing the user to liftover
    or re-munge.
    """
    if identity_mode_family(identifier_mode) != "chr_pos":
        return
    panel_build = getattr(ldscore_result.config_snapshot, "genome_build", None)
    snapshot = sumstats_table.config_snapshot
    sumstats_build = getattr(snapshot, "genome_build", None) if snapshot is not None else None
    if sumstats_build is None:
        if not {"CHR", "POS"}.issubset(sumstats_table.data.columns):
            raise LDSCInputError(
                f"Regression could not verify the genome build of sumstats '{context}': it carries no build "
                "metadata and has no CHR/POS columns to infer one from. Most likely it is a legacy rsID-keyed "
                "file being run against a chr_pos LD-score panel. Re-munge it with the current `ldsc munge-sumstats` "
                "so it records its build, or liftover it to match the panel, then re-run."
            )
        try:
            sumstats_build = infer_chr_pos_build(sumstats_table.data, context=context).genome_build
        except LDSCInputError as exc:
            raise LDSCInputError(
                f"Regression could not verify the genome build of sumstats '{context}': its coordinates did not "
                "match the reference closely enough to infer a build, and it carries no build metadata. "
                "Re-munge it with the current `ldsc munge-sumstats` so it records its build, or liftover it to "
                "match the LD-score panel, then re-run."
            ) from exc
        LOGGER.info(f"Inferred genome build '{sumstats_build}' for sumstats '{context}'.")
    if panel_build is not None and sumstats_build != panel_build:
        raise LDSCInputError(
            f"Regression genome-build mismatch: sumstats '{context}' is {sumstats_build!r} but the LD-score "
            f"panel is {panel_build!r}. Coordinate merges require the same build. Liftover the sumstats to "
            f"{panel_build!r} and re-run."
        )
    LOGGER.info(f"Regression genome build for sumstats '{context}': {sumstats_build!r} (panel {panel_build!r}).")


def _snapshot_identity_mode(snapshot: GlobalConfig | None) -> str | None:
    """Return a normalized SNP identity mode from ``snapshot`` if present."""
    value = getattr(snapshot, "snp_identifier", None)
    if value is None:
        return None
    return normalize_snp_identifier_mode(value)


def _is_legacy_sumstats(table: SumstatsTable) -> bool:
    """Return whether ``table`` came from an LDSC2 munged text artifact."""
    return bool(table.provenance.get("legacy_ldsc2", False)) and "legacy_projection" not in table.provenance


def _has_legacy_sumstats_source(table: SumstatsTable) -> bool:
    """Return whether ``table`` originated from legacy LDSC2 text."""
    return bool(table.provenance.get("legacy_ldsc2", False))


def _project_legacy_sumstats_to_panel(
    table: SumstatsTable,
    ldscore_result: LDScoreSource | LDScoreResult,
) -> tuple[SumstatsTable, pd.DataFrame]:
    """Project one LDSC2 munged-sumstats table onto canonical panel identity.

    Legacy ``SNP`` values are treated only as rsID lookup keys. The supplied
    LD-score panel owns output coordinates and identifier mode. Allele-aware
    panels also own allele order; allele-unaware panels retain the source
    alleles so two traits can still be harmonized for genetic correlation.

    Parameters
    ----------
    table : SumstatsTable
        Legacy text table marked with ``provenance['legacy_ldsc2']``.
    ldscore_result : LDScoreSource or LDScoreResult
        Canonical LDSC3 LD-score panel used by the pending regression.

    Returns
    -------
    projected : SumstatsTable
        New working table under the panel's identity configuration. The input
        object and source file are not mutated.
    dropped : pandas.DataFrame
        Stable row-level audit records for rejected source rows.

    Raises
    ------
    LDSCInputError
        If required legacy columns or canonical panel provenance are missing,
        or no source row can be projected safely.
    """
    required = {"SNP", "A1", "A2", "Z", "N"}
    missing = sorted(required - set(table.data.columns))
    if missing:
        raise LDSCInputError(
            f"Cannot use legacy LDSC2 sumstats '{table.source_path or table.trait_name or 'sumstats'}': "
            f"required columns are missing: {missing}. Legacy regression compatibility requires SNP, A1, A2, Z, "
            "and N; allele-less legacy inputs are unsupported for h2, partitioned-h2, and rg."
        )
    panel_snapshot = ldscore_result.config_snapshot
    panel_mode = _snapshot_identity_mode(panel_snapshot)
    if panel_snapshot is None or panel_mode is None:
        raise LDSCInputError(
            "Legacy LDSC2 sumstats projection requires a canonical LD-score panel with recorded "
            "snp_identifier metadata. Regenerate or explicitly convert the LD-score directory first."
        )
    panel_source = ldscore_result.baseline_table
    missing_panel = sorted({"CHR", "POS", "SNP"} - set(panel_source.columns))
    if missing_panel:
        raise LDSCInputError(
            f"Legacy LDSC2 sumstats projection cannot use the LD-score panel because it is missing {missing_panel}. "
            "Regenerate the canonical LD-score directory."
        )
    allele_aware = is_allele_aware_mode(panel_mode)
    if allele_aware and not {"A1", "A2"}.issubset(panel_source.columns):
        raise LDSCInputError(
            f"Legacy LDSC2 sumstats projection requires panel A1/A2 for snp_identifier='{panel_mode}', "
            "but the canonical LD-score table does not contain both allele columns."
        )

    source = table.data.copy().reset_index(drop=True)
    source["SNP"] = source["SNP"].astype("string").str.strip()
    source["A1"] = source["A1"].astype("string").str.upper().str.strip()
    source["A2"] = source["A2"].astype("string").str.upper().str.strip()
    source_rsids = set(source["SNP"].dropna().astype(str)) - {""}
    panel_snp = panel_source["SNP"].astype("string").str.strip()
    panel_columns = ["CHR", "POS", "SNP", *(["A1", "A2"] if allele_aware else [])]
    panel = panel_source.loc[panel_snp.isin(source_rsids), panel_columns].reset_index(drop=True).copy()
    panel["SNP"] = panel["SNP"].astype("string").str.strip()
    if allele_aware:
        panel["A1"] = panel["A1"].astype("string").str.upper().str.strip()
        panel["A2"] = panel["A2"].astype("string").str.upper().str.strip()

    invalid_frq_count = 0
    if "FRQ" in source.columns:
        raw_frq = source["FRQ"]
        numeric_frq = pd.to_numeric(raw_frq, errors="coerce")
        invalid_frq = raw_frq.notna() & (~numeric_frq.between(0.0, 1.0, inclusive="both"))
        invalid_frq_count = int(invalid_frq.sum())
        source["FRQ"] = numeric_frq.where(numeric_frq.between(0.0, 1.0, inclusive="both"))
        if invalid_frq_count:
            LOGGER.warning(
                "Legacy sumstats '%s' has %d nonnumeric, infinite, or out-of-range FRQ value(s); "
                "they were retained as missing because regression does not consume FRQ.",
                table.source_path or table.trait_name or "sumstats",
                invalid_frq_count,
            )

    candidate_counts = panel["SNP"].value_counts(sort=False, dropna=False)
    source_candidate_count = source["SNP"].map(candidate_counts).fillna(0).astype("int64")
    duplicate_source = source["SNP"].notna() & source["SNP"].duplicated(keep=False)
    valid_bases = {"A", "C", "G", "T"}
    reasons = pd.Series(pd.NA, index=source.index, dtype="string")

    def assign_reason(mask: pd.Series, reason: str) -> None:
        reasons.loc[reasons.isna() & mask.fillna(False)] = reason

    assign_reason(duplicate_source, "duplicate_source_rsid")
    assign_reason(
        source["A1"].isna() | source["A2"].isna() | source["A1"].eq("") | source["A2"].eq(""),
        "missing_allele",
    )
    assign_reason(
        ~source["A1"].isin(valid_bases) | ~source["A2"].isin(valid_bases) | source["A1"].eq(source["A2"]),
        "invalid_allele",
    )
    assign_reason(
        ((source["A1"] == "A") & (source["A2"] == "T"))
        | ((source["A1"] == "T") & (source["A2"] == "A"))
        | ((source["A1"] == "C") & (source["A2"] == "G"))
        | ((source["A1"] == "G") & (source["A2"] == "C")),
        "strand_ambiguous",
    )
    assign_reason(source_candidate_count.eq(0), "missing_panel_rsid")

    eligible = reasons.isna()
    selected_panel: pd.DataFrame
    flip_by_source = pd.Series(False, index=source.index, dtype=bool)
    if allele_aware:
        source_candidates = source.loc[eligible, ["SNP", "A1", "A2"]].copy()
        source_candidates["_legacy_source_row"] = source_candidates.index
        panel_candidates = panel.rename(
            columns={column: f"_legacy_panel_{column}" for column in ["CHR", "POS", "SNP", "A1", "A2"]}
        )
        pairs = source_candidates.merge(
            panel_candidates,
            how="inner",
            left_on="SNP",
            right_on="_legacy_panel_SNP",
            sort=False,
            validate="many_to_many",
        )
        complement_a1 = pairs["A1"].map({"A": "T", "T": "A", "C": "G", "G": "C"})
        complement_a2 = pairs["A2"].map({"A": "T", "T": "A", "C": "G", "G": "C"})
        direct = pairs["A1"].eq(pairs["_legacy_panel_A1"]) & pairs["A2"].eq(pairs["_legacy_panel_A2"])
        complemented = complement_a1.eq(pairs["_legacy_panel_A1"]) & complement_a2.eq(
            pairs["_legacy_panel_A2"]
        )
        swapped = pairs["A1"].eq(pairs["_legacy_panel_A2"]) & pairs["A2"].eq(
            pairs["_legacy_panel_A1"]
        )
        swapped_complement = complement_a1.eq(pairs["_legacy_panel_A2"]) & complement_a2.eq(
            pairs["_legacy_panel_A1"]
        )
        compatible = direct | complemented | swapped | swapped_complement
        compatible_pairs = pairs.loc[compatible].copy()
        compatible_pairs["_legacy_flip"] = (swapped | swapped_complement).loc[compatible_pairs.index].to_numpy()
        compatible_counts = compatible_pairs.groupby("_legacy_source_row", sort=False).size()
        compatible_count_by_source = source.index.to_series().map(compatible_counts).fillna(0).astype("int64")
        assign_reason(eligible & compatible_count_by_source.eq(0), "incompatible_alleles")
        assign_reason(eligible & compatible_count_by_source.gt(1), "ambiguous_panel_mapping")
        accepted_rows = reasons[reasons.isna()].index
        selected_panel = compatible_pairs.loc[
            compatible_pairs["_legacy_source_row"].isin(accepted_rows)
        ].set_index("_legacy_source_row")
        selected_panel = selected_panel.loc[accepted_rows]
        flip_by_source.loc[accepted_rows] = selected_panel["_legacy_flip"].to_numpy()
    else:
        assign_reason(eligible & source_candidate_count.ne(1), "ambiguous_panel_mapping")
        accepted_rows = reasons[reasons.isna()].index
        source_matches = source.loc[accepted_rows, ["SNP"]].copy()
        source_matches["_legacy_source_row"] = source_matches.index
        unique_panel = panel.loc[panel["SNP"].isin(source_matches["SNP"])]
        selected_panel = source_matches.merge(unique_panel, how="left", on="SNP", sort=False, validate="one_to_one")
        selected_panel = selected_panel.set_index("_legacy_source_row").loc[accepted_rows]

    accepted_rows = reasons[reasons.isna()].index
    if accepted_rows.empty:
        counts = reasons.dropna().value_counts().to_dict()
        raise LDSCInputError(
            f"Legacy LDSC2 sumstats '{table.source_path or table.trait_name or 'sumstats'}' retained no SNPs "
            f"after projection to the canonical LD-score panel (drop counts: {counts})."
        )
    projected_frame = source.loc[accepted_rows].copy()
    panel_column = (lambda column: f"_legacy_panel_{column}") if allele_aware else (lambda column: column)
    projected_frame["CHR"] = selected_panel.loc[accepted_rows, panel_column("CHR")].to_numpy()
    projected_frame["POS"] = selected_panel.loc[accepted_rows, panel_column("POS")].to_numpy()
    projected_frame["SNP"] = selected_panel.loc[accepted_rows, panel_column("SNP")].to_numpy()
    if allele_aware:
        projected_frame["A1"] = selected_panel.loc[accepted_rows, panel_column("A1")].to_numpy()
        projected_frame["A2"] = selected_panel.loc[accepted_rows, panel_column("A2")].to_numpy()
        flipped_rows = flip_by_source.loc[accepted_rows]
        projected_frame.loc[flipped_rows, "Z"] = -pd.to_numeric(
            projected_frame.loc[flipped_rows, "Z"], errors="raise"
        )
        if "FRQ" in projected_frame.columns:
            valid_flipped_frq = flipped_rows & projected_frame["FRQ"].notna()
            projected_frame.loc[valid_flipped_frq, "FRQ"] = 1.0 - projected_frame.loc[valid_flipped_frq, "FRQ"]
    projected_frame = projected_frame.reset_index(drop=True)
    projected_frame["POS"] = pd.to_numeric(projected_frame["POS"], errors="raise").astype("int64")
    dropped_mask = reasons.notna()
    drops = source.loc[dropped_mask, ["SNP", "A1", "A2"]].copy()
    drops.insert(0, "source_path", table.source_path)
    drops.insert(0, "trait_name", table.trait_name or Path(table.source_path or "sumstats").name)
    drops["reason"] = reasons.loc[dropped_mask].to_numpy()
    drops["panel_candidate_count"] = source_candidate_count.loc[dropped_mask].to_numpy()
    drops = drops.loc[:, LEGACY_SUMSTATS_DROP_COLUMNS].reset_index(drop=True)
    drop_counts = drops["reason"].value_counts(sort=False)
    if not drop_counts.empty:
        LOGGER.warning(
            "Legacy sumstats '%s' projection dropped %d row(s): %s.",
            table.source_path or table.trait_name or "sumstats",
            len(drops),
            ", ".join(f"{reason}={int(count)}" for reason, count in drop_counts.items()),
        )
    LOGGER.info(
        "Projected %d legacy sumstats row(s) from '%s' onto panel identity mode %s; dropped %d.",
        len(projected_frame),
        table.source_path or table.trait_name or "sumstats",
        panel_mode,
        len(drops),
    )
    projected_table = replace(
        table,
        data=projected_frame,
        has_alleles=True,
        provenance={
            **table.provenance,
            "legacy_projection": {
                "panel_snp_identifier": panel_mode,
                "retained_rows": len(projected_frame),
                "dropped_rows": len(drops),
                "invalid_frq_values": invalid_frq_count,
            },
        },
        config_snapshot=panel_snapshot,
    )
    return projected_table, drops


def _sumstats_identity_mode(
    sumstats_table: SumstatsTable,
    runner_config: GlobalConfig,
    *,
    fallback_mode: str,
) -> str:
    """Return the mode represented by one sumstats input for regression."""
    return _snapshot_identity_mode(sumstats_table.config_snapshot) or normalize_snp_identifier_mode(fallback_mode)


def _ldscore_identity_mode(
    ldscore_result: LDScoreSource | LDScoreResult,
    runner_config: GlobalConfig,
    *,
    fallback_mode: str | None = None,
) -> str:
    """Return the mode represented by one LD-score input for regression."""
    return (
        _snapshot_identity_mode(ldscore_result.config_snapshot)
        or normalize_snp_identifier_mode(fallback_mode or runner_config.snp_identifier)
    )


def _resolve_h2_identity(
    sumstats_table: SumstatsTable,
    ldscore_result: LDScoreSource | LDScoreResult,
    runner_config: GlobalConfig,
    regression_config: RegressionConfig,
):
    """Resolve one sumstats table and one LD-score result for h2 merging."""
    ldscore_mode = _ldscore_identity_mode(ldscore_result, runner_config)
    sumstats_mode = _sumstats_identity_mode(sumstats_table, runner_config, fallback_mode=ldscore_mode)
    return resolve_regression_identity_mode(
        sumstats_mode,
        ldscore_mode,
        allow_identity_downgrade=regression_config.allow_identity_downgrade,
    )


def _resolve_rg_identity(
    sumstats_table_1: SumstatsTable,
    sumstats_table_2: SumstatsTable,
    ldscore_result: LDScoreSource | LDScoreResult,
    runner_config: GlobalConfig,
    regression_config: RegressionConfig,
):
    """Resolve trait 1, trait 2, and LD-score identity modes for rg."""
    ldscore_mode = _ldscore_identity_mode(ldscore_result, runner_config)
    trait_1_mode = _sumstats_identity_mode(sumstats_table_1, runner_config, fallback_mode=ldscore_mode)
    trait_2_mode = _sumstats_identity_mode(sumstats_table_2, runner_config, fallback_mode=ldscore_mode)
    pairwise = [
        resolve_regression_identity_mode(
            trait_1_mode,
            ldscore_mode,
            allow_identity_downgrade=regression_config.allow_identity_downgrade,
        ),
        resolve_regression_identity_mode(
            trait_2_mode,
            ldscore_mode,
            allow_identity_downgrade=regression_config.allow_identity_downgrade,
        ),
        resolve_regression_identity_mode(
            trait_1_mode,
            trait_2_mode,
            allow_identity_downgrade=regression_config.allow_identity_downgrade,
        ),
    ]
    if len({trait_1_mode, trait_2_mode, ldscore_mode}) == 1:
        return pairwise[0]
    effective_mode = identity_base_mode(trait_1_mode)
    return replace(pairwise[0], effective_mode=effective_mode, downgrade_applied=True)


def _prepare_regression_identity_table(
    frame: pd.DataFrame,
    mode: str,
    *,
    context: str,
    logger: logging.Logger,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Clean duplicate effective identities before downgraded regression merge."""
    cleanup = clean_identity_artifact_table(
        frame,
        mode,
        context=context,
        stage="regression_identity_downgrade",
        logger=logger,
    )
    return cleanup.cleaned, cleanup.dropped


def _with_effective_identity_key(frame: pd.DataFrame, mode: str, *, context: str) -> pd.DataFrame:
    """Return a copy with the active effective regression identity key."""
    keyed = frame.copy()
    keyed[REGRESSION_IDENTITY_KEY_COLUMN] = effective_merge_key_series(keyed, mode, context=context)
    return keyed.loc[keyed[REGRESSION_IDENTITY_KEY_COLUMN].notna()].reset_index(drop=True)


def _orient_sumstats_z_to_reference_alleles(merged: pd.DataFrame) -> pd.DataFrame:
    """Orient h2 sumstats Z values to LD-score/reference-panel allele order."""
    sumstats_a1 = merged["A1"].astype("string").str.upper()
    sumstats_a2 = merged["A2"].astype("string").str.upper()
    ld_a1 = merged["A1_ld"].astype("string").str.upper()
    ld_a2 = merged["A2_ld"].astype("string").str.upper()
    allele_pairs = sumstats_a1 + sumstats_a2 + ld_a1 + ld_a2
    compatible = allele_pairs.isin(reg.MATCH_ALLELES)
    incompatible = ~compatible
    if bool(incompatible.any()):
        LOGGER.warning(
            "Dropping %d SNPs with incompatible sumstats and LD-score allele order.",
            int(incompatible.sum()),
        )
    oriented = merged.loc[~incompatible].copy()
    flip_index = allele_pairs.loc[~incompatible].map(reg.FLIP_ALLELES).astype(bool)
    oriented.loc[flip_index.index[flip_index], "Z"] *= -1
    return oriented.reset_index(drop=True)


def _with_chr_pos_key(frame: pd.DataFrame, *, context: str) -> pd.DataFrame:
    """Return a copy with a private canonical ``CHR:POS`` merge key."""
    chr_col, pos_col = infer_chr_pos_columns(frame.columns, context=context)
    keyed, _report = build_chr_pos_key_frame(
        frame,
        context=context,
        chr_col=chr_col,
        pos_col=pos_col,
        drop_missing=True,
        logger=LOGGER,
    )
    return keyed


def _select_count_key(
    count_totals: dict[str, np.ndarray],
    use_common_counts: bool,
    columns: Sequence[str] = (),
) -> str:
    """Pick the regression count vector key, preferring common-SNP counts by default."""
    if use_common_counts and COMMON_COUNT_KEY in count_totals:
        return COMMON_COUNT_KEY
    if use_common_counts and ALL_COUNT_KEY in count_totals:
        return ALL_COUNT_KEY
    if not use_common_counts and ALL_COUNT_KEY in count_totals:
        return ALL_COUNT_KEY
    if not use_common_counts:
        raise LDSCInputError(
            "Regression requested all-SNP counts, but all-SNP counts are unavailable for LD-score column(s) "
            f"{list(columns)}. This is expected for an imported unpartitioned LDSC2 suite when one or more "
            "chromosome .M files were missing. Use the default common counts or convert a suite with complete .M files."
        )
    raise LDSCInputError(
        "Regression could not find a complete reference-SNP count vector for the requested LD-score columns. "
        "Regenerate or reconvert the LD-score directory."
    )


def _validate_partitioned_query_columns(ldscore_result: LDScoreSource | LDScoreResult, query_columns: Sequence[str]) -> list[str]:
    """Validate that requested cell-type query columns exist in the LD-score result.

    The functional regime (no query columns) is decided by the caller; this
    helper is only used on the cell-type path, so it just checks that every
    requested query column is available.
    """
    requested = list(query_columns)
    available = list(ldscore_result.query_columns)
    missing = [column for column in requested if column not in available]
    if missing:
        raise LDSCInputError(
            f"partitioned-h2 could not find requested query annotation columns {missing} in the LD-score directory. "
            f"Most likely the query annotation name was misspelled or the LD-score directory was generated from "
            f"different query inputs. Choose one of the available query annotations: {available}."
        )
    return requested


def _assemble_regression_ldscore_table(
    ldscore_result: LDScoreSource | LDScoreResult,
    query_columns: Sequence[str],
    *,
    snp_identifier: str | None = None,
) -> pd.DataFrame:
    """Build the table used for one regression dataset from split LD-score tables."""
    baseline_table = ldscore_result.baseline_table.reset_index(drop=True)
    if not query_columns:
        return baseline_table.copy()
    query_table = ldscore_result.read_queries(query_columns)
    assert_same_snp_rows(
        baseline_table,
        query_table,
        context="query rows must match baseline rows on CHR/SNP/POS",
        snp_identifier=snp_identifier or getattr(ldscore_result.config_snapshot, "snp_identifier", "chr_pos_allele_aware"),
    )
    return pd.concat([baseline_table, query_table.loc[:, list(query_columns)]], axis=1)


def _count_totals_for_columns(count_records: Sequence[dict[str, Any]], columns: Sequence[str]) -> dict[str, np.ndarray]:
    """Return LDSC count vectors aligned to ``columns`` using metadata records."""
    records_by_column = {str(record["column"]): record for record in count_records}
    missing = [column for column in columns if column not in records_by_column]
    if missing:
        raise LDSCInputError(
            f"Regression cannot align LD-score count metadata: count records are missing for columns {missing}. "
            "Most likely the LD-score metadata was written by an older version or edited separately from "
            "the parquet tables. Regenerate the LD-score directory with the current `ldsc ldscore`."
        )
    count_totals: dict[str, np.ndarray] = {}
    all_values = [records_by_column[column].get("all_reference_snp_count") for column in columns]
    if all(value is not None and pd.notna(value) for value in all_values):
        all_counts = np.asarray([float(value) for value in all_values], dtype=np.float64)
        if not np.isfinite(all_counts).all():
            raise LDSCInputError("Regression found non-finite all-SNP count metadata in the LD-score directory.")
        count_totals[ALL_COUNT_KEY] = all_counts
    if all("common_reference_snp_count" in records_by_column[column] for column in columns):
        common_values = [records_by_column[column].get("common_reference_snp_count") for column in columns]
        if any(value is None or pd.isna(value) for value in common_values):
            raise LDSCInputError(
                "Regression found a malformed common-SNP count record with a missing value. "
                "Reconvert or regenerate the LD-score directory; common counts never fall back silently."
            )
        common_counts = np.asarray([float(value) for value in common_values], dtype=np.float64)
        if not np.isfinite(common_counts).all():
            raise LDSCInputError("Regression found non-finite common-SNP count metadata in the LD-score directory.")
        count_totals[COMMON_COUNT_KEY] = common_counts
    return count_totals


def _is_missing_prevalence(value) -> bool:
    """True when a prevalence is unset (None) or a NaN quantitative-trait marker."""
    return value is None or (isinstance(value, float) and math.isnan(value))


def _liability_pair(obs, obs_se, samp_prev, pop_prev):
    """Return (liab, liab_se) for an absolute h2 estimate, or (nan, nan).

    Liability is defined only when both prevalences are finite probabilities;
    a quantitative trait (None/NaN) or unset prevalence yields NaN. ``obs`` and
    ``obs_se`` are the observed-scale point estimate and standard error.
    """
    if _is_missing_prevalence(samp_prev) or _is_missing_prevalence(pop_prev):
        return float("nan"), float("nan")
    c = float(reg.liability_conversion_factor(samp_prev, pop_prev))
    return obs * c, obs_se * c


def _prev_cell(value):
    """Normalize a prevalence for a table cell: None -> NaN, else float."""
    return float("nan") if value is None else float(value)


def _config_prevalence(config: "RegressionConfig | None") -> tuple[float | None, float | None]:
    """Scalar (samp_prev, pop_prev) from a regression config, or (None, None)."""
    if config is None:
        return None, None
    return getattr(config, "samp_prev", None), getattr(config, "pop_prev", None)


def _nan_to_none(value):
    """Map a NaN prevalence cell back to None for clean JSON metadata."""
    return None if (isinstance(value, float) and math.isnan(value)) else value


def _scale_label(*prevalences) -> str:
    """Return 'liability' if any supplied prevalence is finite, else 'observed'."""
    return "liability" if any(p is not None for p in prevalences) else "observed"


def summarize_total_h2(
    hsq,
    dataset: RegressionDataset,
    trait_name: str | None = None,
    n_snps_used: int | None = None,
    samp_prev: float | None = None,
    pop_prev: float | None = None,
) -> pd.DataFrame:
    """Build the one-row total-heritability summary from a fitted ``Hsq`` result.

    Reports total heritability on both the observed (``total_h2_obs``) and
    liability (``total_h2_liab``) scales, plus the prevalences applied. Liability
    columns are NaN unless both ``samp_prev`` and ``pop_prev`` are finite
    probabilities in (0, 1). ``n_snps_used`` is the post-chi-square-filter SNP
    count from the internal fit outcome; when omitted it falls back
    to the supplied dataset's row count. Workflows supply the fitted dataset.
    """
    obs = _scalar(hsq.tot)
    obs_se = _scalar(hsq.tot_se)
    liab, liab_se = _liability_pair(obs, obs_se, samp_prev, pop_prev)
    return pd.DataFrame(
        [
            {
                "trait_name": trait_name,
                "n_snps": int(n_snps_used) if n_snps_used is not None else len(dataset.merged),
                "total_h2_obs": obs,
                "total_h2_obs_se": obs_se,
                "total_h2_liab": liab,
                "total_h2_liab_se": liab_se,
                "samp_prev": _prev_cell(samp_prev),
                "pop_prev": _prev_cell(pop_prev),
                "intercept": _scalar_or_value(hsq.intercept),
                "intercept_se": getattr(hsq, "intercept_se", None),
                "mean_chisq": _scalar_or_value(hsq.mean_chisq),
                "lambda_gc": _scalar_or_value(hsq.lambda_gc),
                "ratio": getattr(hsq, "ratio", None),
                "ratio_se": getattr(hsq, "ratio_se", None),
            }
        ]
    )


def summarize_ld_score_regression_bins(
    hsq,
    *,
    ld_score,
    chi_square,
    sample_size,
    regression_ld_score,
    reference_snp_count: float,
    max_bins: int = 50,
) -> pd.DataFrame:
    """Summarize the exact fitted SNP population into LD Score rank bins.

    Parameters
    ----------
    hsq : ldsc._kernel.regression.Hsq-like object
        Final unpartitioned fit exposing ``coef``, ``tot``, and ``intercept``.
    ld_score, chi_square, sample_size, regression_ld_score : array_like
        Equal-length per-SNP values after the h2 chi-square filter. The order is
        preserved among tied LD Scores.
    reference_snp_count : float
        Reference-SNP count used by the fitted h2 model.
    max_bins : int, optional
        Maximum number of nonempty equal-count rank bins. Default is 50.

    Returns
    -------
    pandas.DataFrame
        One row per ascending LD Score bin in
        :data:`ldsc.outputs.H2_REGRESSION_BIN_COLUMNS` order. Chi-square SD uses
        ``ddof=1`` and is NaN for a singleton bin.

    Raises
    ------
    LDSCInternalError
        If the fitted arrays are empty or differ in length, or if ``max_bins``
        is less than one.

    Notes
    -----
    Fitted expectations use the final LDSC slope and intercept with each SNP's
    own sample size. Weights are the existing heritability regression weights
    evaluated at the final model; this function does not fit another model.
    """
    arrays = [
        np.ravel(np.asarray(values, dtype=np.float64))
        for values in (ld_score, chi_square, sample_size, regression_ld_score)
    ]
    lengths = {len(values) for values in arrays}
    if len(lengths) != 1 or not lengths or next(iter(lengths)) == 0:
        raise LDSCInternalError(
            "h2 regression diagnostics require nonempty, equal-length LD Score, chi-square, "
            "sample-size, and regression-weight LD Score arrays."
        )
    if max_bins < 1:
        raise LDSCInternalError("h2 regression diagnostics require max_bins >= 1.")
    ld_values, chisq_values, n_values, w_ld_values = arrays
    n_snps = len(ld_values)
    n_bins = min(int(max_bins), n_snps)
    slope = _scalar(hsq.coef)
    total_h2 = _scalar(hsq.tot)
    intercept = _scalar(hsq.intercept)
    fitted = intercept + np.multiply(np.multiply(n_values, slope), ld_values)
    weights = np.ravel(
        reg.Hsq.weights(
            ld_values.reshape((-1, 1)),
            w_ld_values.reshape((-1, 1)),
            n_values.reshape((-1, 1)),
            float(reference_snp_count),
            total_h2,
            intercept,
        )
    )
    order = np.argsort(ld_values, kind="stable")
    rows: list[dict[str, float | int]] = []
    for bin_index, indices in enumerate(np.array_split(order, n_bins), start=1):
        bin_chisq = chisq_values[indices]
        rows.append(
            {
                "bin": bin_index,
                "n_snps": len(indices),
                "ld_score_min": float(np.min(ld_values[indices])),
                "ld_score_max": float(np.max(ld_values[indices])),
                "mean_ld_score": float(np.mean(ld_values[indices])),
                "mean_chi_square": float(np.mean(bin_chisq)),
                "sd_chi_square": float(np.std(bin_chisq, ddof=1)) if len(indices) > 1 else float("nan"),
                "mean_sample_size": float(np.mean(n_values[indices])),
                "mean_fitted_chi_square": float(np.mean(fitted[indices])),
                "mean_regression_weight": float(np.mean(weights[indices])),
            }
        )
    return pd.DataFrame(rows, columns=H2_REGRESSION_BIN_COLUMNS)


def summarize_partitioned_h2(
    hsq,
    dataset: RegressionDataset,
    annotation_columns: Sequence[str],
    samp_prev: float | None = None,
    pop_prev: float | None = None,
) -> pd.DataFrame:
    """Build overlap-aware partitioned-h2 rows from a fitted ``Hsq`` result.

    Assembles the model overlap matrix for the retained LD-score columns, runs
    the ported legacy ``_overlap_output`` (augmented with the one-sided
    coefficient p-value, conditional category heritability, and the overlap
    flag), and returns the requested ``annotation_columns`` rows in the single
    lowercase ``PARTITIONED_H2_COLUMNS`` schema. The whole-model total heritability
    is carried on every row (``total_h2_obs``, constant across one fitted model)
    so each row's ``category_h2_obs`` can be read against it. Per-category and
    total heritability are reported on both the observed (``*_obs``) and liability
    (``*_liab``) scales; liability columns are NaN unless both prevalences are
    finite probabilities in (0, 1). Proportions, enrichment, and coefficients are
    scale-invariant. ``annotation_columns`` must be a subset of
    ``dataset.retained_ld_columns``.
    """
    from .overlap_matrix import assemble_model_overlap, overlap_aware_category_table

    overlap = dataset.ldscore_overlap
    if overlap is None:
        raise LDSCInputError(PARTITIONED_H2_REQUIRES_OVERLAP_MESSAGE)
    use_common = dataset.count_key_used_for_regression == COMMON_COUNT_KEY
    model_overlap = assemble_model_overlap(overlap, list(dataset.retained_ld_columns), use_common)
    m_annot = np.asarray(
        dataset.reference_snp_count_totals[dataset.count_key_used_for_regression], dtype=np.float64
    )
    m_tot = overlap.total_common_reference_snps if use_common else overlap.total_all_reference_snps
    table = overlap_aware_category_table(
        hsq, model_overlap, m_annot, float(m_tot), dataset.retained_ld_columns
    )
    table = table.rename(
        columns={
            "Category": "category",
            "Prop._SNPs": "prop_snps",
            "Prop._h2": "prop_h2",
            "Prop._h2_std_error": "prop_h2_se",
            "Enrichment": "enrichment",
            "Enrichment_std_error": "enrichment_se",
            "Enrichment_p": "enrichment_p",
            "Coefficient": "coefficient",
            "Coefficient_std_error": "coefficient_se",
            "Coefficient_z": "coefficient_z",
            "Coefficient_p": "coefficient_p",
            "overlap_aware": "overlap_annot",
            "Category_h2": "category_h2_obs",
            "Category_h2_std_error": "category_h2_obs_se",
        }
    )
    c = (
        float(reg.liability_conversion_factor(samp_prev, pop_prev))
        if not _is_missing_prevalence(samp_prev) and not _is_missing_prevalence(pop_prev)
        else float("nan")
    )
    table["category_h2_liab"] = table["category_h2_obs"] * c
    table["category_h2_liab_se"] = table["category_h2_obs_se"] * c
    # Whole-model total heritability: a single scalar broadcast to every row of one
    # fitted model (constant across the functional regime; each query's own model
    # total in the cell-type regime).
    total_obs = _scalar(hsq.tot)
    total_obs_se = _scalar(hsq.tot_se)
    table["total_h2_obs"] = total_obs
    table["total_h2_obs_se"] = total_obs_se
    table["total_h2_liab"] = total_obs * c
    table["total_h2_liab_se"] = total_obs_se * c
    table["samp_prev"] = _prev_cell(samp_prev)
    table["pop_prev"] = _prev_cell(pop_prev)
    rows = table.set_index("category").loc[list(annotation_columns)].reset_index()
    return rows.loc[:, PARTITIONED_H2_COLUMNS].reset_index(drop=True)


def _iter_rg_pairs(n_traits: int, anchor_index: int | None) -> list[tuple[int, int]]:
    """Return rg pair indices in the public output order."""
    if anchor_index is not None:
        return [(anchor_index, idx) for idx in range(n_traits) if idx != anchor_index]
    return [(i, j) for i in range(n_traits) for j in range(i + 1, n_traits)]


def _summarize_rg_pair(
    rg_result,
    dataset: RGRegressionDataset,
    *,
    trait_1: str,
    trait_2: str,
    pair_kind: str,
    n_snps_used: int | None = None,
    samp_prev_1: float | None = None,
    pop_prev_1: float | None = None,
    samp_prev_2: float | None = None,
    pop_prev_2: float | None = None,
) -> dict[str, object]:
    """Build the full public rg row from one fitted kernel result.

    Per-trait heritability and the pair's genetic covariance are reported on both
    the observed (``*_obs``) and liability (``*_liab``) scales. The rg ratio
    (``rg``/``rg_se``/``z``/``p``) is scale-invariant and unconverted. A trait's
    ``h2_*_liab`` is NaN unless its prevalence is finite; ``gencov_liab`` is NaN
    only when no prevalence was supplied for the run. ``n_snps_used`` is the
    post-product-filter SNP count from the internal fit outcome; it falls
    back to the supplied dataset's row count. Workflows supply the fitted dataset.
    """
    h2_1_obs = _numeric_attr(getattr(rg_result, "hsq1", None), "tot", "h2_1_obs")
    h2_1_obs_se = _numeric_attr(getattr(rg_result, "hsq1", None), "tot_se", "h2_1_obs_se")
    h2_2_obs = _numeric_attr(getattr(rg_result, "hsq2", None), "tot", "h2_2_obs")
    h2_2_obs_se = _numeric_attr(getattr(rg_result, "hsq2", None), "tot_se", "h2_2_obs_se")
    gencov_obs = _numeric_attr(getattr(rg_result, "gencov", None), "tot", "gencov_obs")
    gencov_obs_se = _numeric_attr(getattr(rg_result, "gencov", None), "tot_se", "gencov_obs_se")
    h2_1_liab, h2_1_liab_se = _liability_pair(h2_1_obs, h2_1_obs_se, samp_prev_1, pop_prev_1)
    h2_2_liab, h2_2_liab_se = _liability_pair(h2_2_obs, h2_2_obs_se, samp_prev_2, pop_prev_2)
    if any(v is not None for v in (samp_prev_1, pop_prev_1, samp_prev_2, pop_prev_2)):
        gfac = reg.gencov_obs_to_liab(1.0, samp_prev_1, samp_prev_2, pop_prev_1, pop_prev_2)
        gencov_liab, gencov_liab_se = gencov_obs * gfac, gencov_obs_se * gfac
    else:
        gencov_liab, gencov_liab_se = float("nan"), float("nan")
    return {
        "trait_1": trait_1,
        "trait_2": trait_2,
        "n_snps_used": int(n_snps_used) if n_snps_used is not None else int(len(dataset.merged)),
        "rg": _required_numeric_scalar(getattr(rg_result, "rg_ratio", None), "rg"),
        "rg_se": _required_numeric_scalar(getattr(rg_result, "rg_se", None), "rg_se"),
        "z": _required_numeric_scalar(getattr(rg_result, "z", None), "z"),
        "p": _required_numeric_scalar(getattr(rg_result, "p", None), "p"),
        "h2_1_obs": h2_1_obs,
        "h2_1_obs_se": h2_1_obs_se,
        "h2_1_liab": h2_1_liab,
        "h2_1_liab_se": h2_1_liab_se,
        "h2_2_obs": h2_2_obs,
        "h2_2_obs_se": h2_2_obs_se,
        "h2_2_liab": h2_2_liab,
        "h2_2_liab_se": h2_2_liab_se,
        "gencov_obs": gencov_obs,
        "gencov_obs_se": gencov_obs_se,
        "gencov_liab": gencov_liab,
        "gencov_liab_se": gencov_liab_se,
        "samp_prev_1": _prev_cell(samp_prev_1),
        "pop_prev_1": _prev_cell(pop_prev_1),
        "samp_prev_2": _prev_cell(samp_prev_2),
        "pop_prev_2": _prev_cell(pop_prev_2),
        "intercept_h2_1": _numeric_attr(getattr(rg_result, "hsq1", None), "intercept", "intercept_h2_1"),
        "intercept_h2_1_se": _numeric_attr(getattr(rg_result, "hsq1", None), "intercept_se", "intercept_h2_1_se"),
        "intercept_h2_2": _numeric_attr(getattr(rg_result, "hsq2", None), "intercept", "intercept_h2_2"),
        "intercept_h2_2_se": _numeric_attr(getattr(rg_result, "hsq2", None), "intercept_se", "intercept_h2_2_se"),
        "intercept_gencov": _numeric_attr(getattr(rg_result, "gencov", None), "intercept", "intercept_gencov"),
        "intercept_gencov_se": _numeric_attr(
            getattr(rg_result, "gencov", None),
            "intercept_se",
            "intercept_gencov_se",
        ),
        "ratio_1": _numeric_attr(getattr(rg_result, "hsq1", None), "ratio", "ratio_1"),
        "ratio_1_se": _numeric_attr(getattr(rg_result, "hsq1", None), "ratio_se", "ratio_1_se"),
        "ratio_2": _numeric_attr(getattr(rg_result, "hsq2", None), "ratio", "ratio_2"),
        "ratio_2_se": _numeric_attr(getattr(rg_result, "hsq2", None), "ratio_se", "ratio_2_se"),
        "lambda_gc_1": _numeric_attr(getattr(rg_result, "hsq1", None), "lambda_gc", "lambda_gc_1"),
        "lambda_gc_2": _numeric_attr(getattr(rg_result, "hsq2", None), "lambda_gc", "lambda_gc_2"),
        "mean_chisq_1": _numeric_attr(getattr(rg_result, "hsq1", None), "mean_chisq", "mean_chisq_1"),
        "mean_chisq_2": _numeric_attr(getattr(rg_result, "hsq2", None), "mean_chisq", "mean_chisq_2"),
        "pair_kind": pair_kind,
        "status": "ok",
        "error": "",
    }


def _failed_rg_full_row(*, trait_1: str, trait_2: str, pair_kind: str, error: str) -> dict[str, object]:
    """Build a full rg row for a failed pair."""
    row = {column: math.nan for column in RG_FULL_COLUMNS}
    row.update(
        {
            "trait_1": trait_1,
            "trait_2": trait_2,
            "n_snps_used": 0,
            "pair_kind": pair_kind,
            "status": "failed",
            "error": error,
        }
    )
    return row


def _concise_rg_row(full_row: dict[str, object]) -> dict[str, object]:
    """Build the concise public rg row from a full row."""
    status = str(full_row.get("status", ""))
    return {
        "trait_1": full_row["trait_1"],
        "trait_2": full_row["trait_2"],
        "n_snps_used": full_row["n_snps_used"],
        "rg": full_row["rg"],
        "rg_se": full_row["rg_se"],
        "p": full_row["p"],
        "note": "" if status == "ok" else FAILED_RG_NOTE,
    }


def _h2_metadata(
    args,
    sumstats_table: SumstatsTable,
    outcome: _FitOutcome,
    samp_prev: float | None = None,
    pop_prev: float | None = None,
) -> dict[str, object]:
    """Build h2 metadata from realized fit selection and validated prevalences."""
    dataset = outcome.dataset
    config_snapshot = dataset.config_snapshot
    return {
        "artifact_type": "h2_result",
        "trait_name": sumstats_table.trait_name,
        "sumstats_file": getattr(args, "sumstats_file", None),
        "ldscore_dir": getattr(args, "ldscore_dir", None),
        "effective_snp_identifier": dataset.effective_snp_identifier,
        "genome_build": None if config_snapshot is None else getattr(config_snapshot, "genome_build", None),
        "identity_downgrade_applied": dataset.identity_downgrade_applied,
        "count_key_used_for_regression": dataset.count_key_used_for_regression,
        "retained_ld_columns": list(dataset.retained_ld_columns),
        "dropped_zero_variance_ld_columns": list(dataset.dropped_zero_variance_ld_columns),
        "n_snps": outcome.n_snps,
        "effective_chisq_max": outcome.effective_chisq_max,
        "n_blocks_used": outcome.n_blocks,
        "samp_prev": samp_prev,
        "pop_prev": pop_prev,
        "scale": _scale_label(samp_prev, pop_prev),
    }


def _rg_prevalence_metadata(full_row: dict[str, object]) -> dict[str, object]:
    """Per-pair prevalence + scale provenance pulled from a full rg row."""
    prev = {
        key: _nan_to_none(full_row.get(key))
        for key in ("samp_prev_1", "pop_prev_1", "samp_prev_2", "pop_prev_2")
    }
    return {**prev, "scale": _scale_label(*prev.values())}


def _rg_pair_metadata(
    table_1: SumstatsTable,
    table_2: SumstatsTable,
    outcome: _FitOutcome,
    full_row: dict[str, object],
    config: RegressionConfig,
    pair_kind: str,
) -> dict[str, object]:
    """Build per-pair metadata for an estimated rg pair."""
    dataset = outcome.dataset
    return {
        "trait_1": full_row["trait_1"],
        "trait_2": full_row["trait_2"],
        "source_1": table_1.source_path,
        "source_2": table_2.source_path,
        "pair_kind": pair_kind,
        "status": full_row["status"],
        "error": full_row["error"],
        "n_snps_used": outcome.n_snps,
        "n_blocks_used": outcome.n_blocks,
        "effective_chisq_max": outcome.effective_chisq_max,
        **_rg_prevalence_metadata(full_row),
        "count_key_used_for_regression": dataset.count_key_used_for_regression,
        "retained_ld_columns": list(dataset.retained_ld_columns),
        "dropped_zero_variance_ld_columns": list(dataset.dropped_zero_variance_ld_columns),
        "effective_snp_identifier": dataset.effective_snp_identifier,
        "identity_downgrade_applied": dataset.identity_downgrade_applied,
        "intercept_h2_policy": _intercept_policy(config.intercept_h2, config.use_intercept, default_when_disabled=1),
        "intercept_gencov_policy": _intercept_policy(
            config.intercept_gencov,
            config.use_intercept,
            default_when_disabled=0,
        ),
    }


def _failed_rg_pair_metadata(
    table_1: SumstatsTable,
    table_2: SumstatsTable,
    full_row: dict[str, object],
    config: RegressionConfig,
    pair_kind: str,
) -> dict[str, object]:
    """Build per-pair metadata when a pair failed before a dataset was available."""
    return {
        "trait_1": full_row["trait_1"],
        "trait_2": full_row["trait_2"],
        "source_1": table_1.source_path,
        "source_2": table_2.source_path,
        "pair_kind": pair_kind,
        "status": "failed",
        "error": full_row["error"],
        "n_snps_used": 0,
        "n_blocks_used": 0,
        "effective_chisq_max": config.chisq_max,
        **_rg_prevalence_metadata(full_row),
        "intercept_h2_policy": _intercept_policy(config.intercept_h2, config.use_intercept, default_when_disabled=1),
        "intercept_gencov_policy": _intercept_policy(
            config.intercept_gencov,
            config.use_intercept,
            default_when_disabled=0,
        ),
    }


def _intercept_policy(value: float | None, use_intercept: bool, default_when_disabled: float) -> str:
    """Describe the effective intercept policy in per-pair metadata."""
    if not use_intercept:
        return f"fixed:{default_when_disabled:g}"
    if value is None:
        return "free"
    return f"fixed:{float(value):g}"


def _trait_label(table: SumstatsTable) -> str:
    """Return the public trait label for a sumstats table."""
    if table.trait_name:
        return str(table.trait_name)
    if table.source_path:
        return Path(table.source_path).name
    return "trait"


def _source_key(table: SumstatsTable) -> str:
    """Return a stable source string for duplicate trait-name disambiguation."""
    return str(Path(table.source_path).resolve()) if table.source_path else _trait_label(table)


def _numeric_attr(obj: object | None, attr: str, field: str) -> float:
    """Return a numeric scalar attribute or NaN when unavailable."""
    if obj is None:
        return math.nan
    return _numeric_scalar(getattr(obj, attr, None), field)


def _required_numeric_scalar(value: object, field: str) -> float:
    """Return a numeric scalar or raise for kernel string/invalid headline values."""
    return _numeric_scalar(value, field, fail_on_non_numeric=True)


def _numeric_scalar(value: object, field: str, *, fail_on_non_numeric: bool = False) -> float:
    """Return the first numeric scalar from common kernel return shapes."""
    if value is None:
        return math.nan
    if isinstance(value, str):
        stripped = value.strip()
        if stripped.upper() == "NA":
            if fail_on_non_numeric:
                raise LDSCInternalError(
                    f"Regression summary field {field} is the non-numeric kernel value 'NA'. "
                    "Most likely the estimator returned an unavailable headline statistic where a number was required. "
                    "Re-run with `--log-level DEBUG` and report the traceback."
                )
            return math.nan
        try:
            return float(stripped)
        except ValueError as exc:
            if fail_on_non_numeric:
                raise LDSCInternalError(
                    f"Regression summary field {field} is non-numeric kernel value {value!r}. "
                    "Most likely the estimator returned an unavailable headline statistic where a number was required. "
                    "Re-run with `--log-level DEBUG` and report the traceback."
                ) from exc
            return math.nan
    if hasattr(value, "__array__") or isinstance(value, (list, tuple)):
        array = np.ravel(value)
        if array.size == 0:
            return math.nan
        return _numeric_scalar(array[0].item() if hasattr(array[0], "item") else array[0], field, fail_on_non_numeric=fail_on_non_numeric)
    try:
        if pd.isna(value):
            return math.nan
    except TypeError:
        pass
    try:
        return float(value)
    except (TypeError, ValueError) as exc:
        if fail_on_non_numeric:
            raise LDSCInternalError(
                f"Regression summary field {field} is non-numeric kernel value {value!r}. "
                "Most likely the estimator returned an unavailable headline statistic where a number was required. "
                "Re-run with `--log-level DEBUG` and report the traceback."
            ) from exc
        return math.nan


def _format_exception(exc: Exception) -> str:
    """Return compact user-facing error text for a pair failure."""
    message = str(exc).strip()
    return f"{exc.__class__.__name__}: {message}" if message else exc.__class__.__name__


def _scalar(value) -> float:
    """Return the first scalar numeric value from a possibly array-like object."""
    return float(np.ravel(value)[0])


def _coefficient_delete_frame(hsq, annotation_columns: Sequence[str]) -> pd.DataFrame:
    """Return float64 delete-one-block coefficient values in fitted order.

    ``Hsq.part_delete_values`` already contains annotation coefficients on the
    per-SNP heritability scale. The intercept delete column is intentionally
    excluded because post-fit annotation projection consumes only ``tau``.
    """
    raw_values = getattr(hsq, "part_delete_values", None)
    if not isinstance(raw_values, (np.ndarray, list, tuple, pd.DataFrame)):
        raw_values = np.empty((0, len(annotation_columns)))
    values = np.asarray(raw_values, dtype=np.float64)
    if values.ndim != 2 or values.shape[1] != len(annotation_columns):
        raise LDSCInternalError(
            "partitioned-h2 coefficient delete values do not match the retained annotation order. "
            f"Expected {len(annotation_columns)} columns, got shape {values.shape}. "
            "Re-run with `--log-level DEBUG` and report the traceback."
        )
    frame = pd.DataFrame(values, columns=list(annotation_columns), dtype=np.float64)
    frame.insert(0, "delete_block", np.arange(len(frame), dtype=np.int64))
    return frame


def _scalar_or_value(value):
    """Return a scalar float for array-like values, otherwise preserve the original value."""
    if hasattr(value, "__array__") or isinstance(value, (list, tuple)):
        return _scalar(value)
    return value


def _resolve_summary_sort(sort_by: str, *, has_queries: bool) -> str:
    """Resolve the ``auto`` summary sort to a regime-appropriate column key.

    Cell-type runs (query annotations present) surface the most significant query
    first by ``coefficient-p``; functional runs preserve baseline order with
    ``category``. Any explicit choice is returned unchanged.
    """
    if sort_by != "auto":
        return sort_by
    return "coefficient-p" if has_queries else "category"


def _log_partitioned_h2_regime(ldscore_result: LDScoreSource | LDScoreResult, has_queries: bool) -> None:
    """Log a one-block banner naming the regime and the column to focus on."""
    if has_queries:
        LOGGER.info(
            "Cell-type-specific regime: baseline + one query per model (%d queries). "
            "Focus on `coefficient` together with the one-sided `coefficient_p` "
            "(a one-sided test of whether `coefficient` > 0): a positive `coefficient` with a small "
            "`coefficient_p` means the query annotation contributes heritability beyond the baseline "
            "annotations. The `enrichment` column is confounded by the query annotation's overlap with "
            "the baseline annotations, so use `coefficient_p` to judge whether the additional "
            "contribution is significant.",
            len(ldscore_result.query_columns),
        )
    else:
        LOGGER.info(
            "Functional-category regime: joint fit of %d baseline categories. "
            "Focus on `enrichment` and `enrichment_p`: `enrichment > 1` means the category's SNPs "
            "explain a larger share of heritability than their share of SNPs (`< 1` means a smaller "
            "share); a small `enrichment_p` indicates the enrichment is significantly different from 1, "
            "i.e. significantly larger or smaller.",
            len(ldscore_result.baseline_columns),
        )


def _log_quantitative_annotation_interpretation(ldscore_result: LDScoreResult) -> None:
    """Warn once when a fitted model contains quantitative annotations."""
    quantitative = [
        name for name, annotation_type in ldscore_result.annotation_types.items()
        if annotation_type == "quantitative"
    ]
    if not quantitative:
        return
    LOGGER.warning(
        "Quantitative annotation(s) detected: %s. Numerical prop_snps, prop_h2, and enrichment summaries "
        "remain visible for compatibility, but they are weighted summaries and do not have the ordinary "
        "binary-category interpretation. Interpret tau (`coefficient`), its SE, and its coefficient zero-test; "
        "use `ldsc quantile-h2` for post-fit quantile enrichment. To estimate distinct conditional coefficients "
        "for bins, construct binary quantile annotations and refit partitioned-h2.",
        ", ".join(quantitative),
    )


def _log_effective_regression_identity(
    sumstats_tables: Sequence[SumstatsTable],
    ldscore_result: LDScoreSource | LDScoreResult,
    runner_config: GlobalConfig,
    regression_config: RegressionConfig,
) -> None:
    """Log the SNP-identity mode regression will actually use, and its source.

    The GlobalConfig banner reports ``snp_identifier`` as a CLI-level default.
    Regression instead derives the effective mode from each input artifact's
    recorded provenance and only falls back to that default when an artifact
    omits it. This line names the resolved mode, the per-input sources, and any
    allele-awareness downgrade so the banner default is not mistaken for the
    mode in force. It is informational only: incompatible inputs are left for
    the estimation step to reject with its authoritative error.
    """
    default_mode = normalize_snp_identifier_mode(runner_config.snp_identifier)
    used_default = False

    def describe(label: str, snapshot: GlobalConfig | None) -> tuple[str, str]:
        nonlocal used_default
        recorded = _snapshot_identity_mode(snapshot)
        if recorded is None:
            used_default = True
            return default_mode, f"{label}={default_mode} (no provenance recorded; using GlobalConfig default)"
        return recorded, f"{label}={recorded} (from provenance)"

    modes: list[str] = []
    sources: list[str] = []
    ld_mode, ld_source = describe("LD-score", ldscore_result.config_snapshot)
    modes.append(ld_mode)
    sources.append(ld_source)
    for table in sumstats_tables:
        mode, source = describe(_trait_label(table), table.config_snapshot)
        modes.append(mode)
        sources.append(source)

    effective = modes[0]
    downgrade = False
    try:
        for mode in modes[1:]:
            compat = resolve_regression_identity_mode(
                effective, mode, allow_identity_downgrade=regression_config.allow_identity_downgrade
            )
            effective = compat.effective_mode
            downgrade = downgrade or compat.downgrade_applied
    except LDSCUsageError:
        # Inputs are incompatible; the estimation step raises the authoritative
        # error. Still surface per-input modes so the conflict is visible early.
        LOGGER.info(
            "Regression SNP identifier could not be resolved: the inputs disagree. "
            "Per-input modes: %s. The GlobalConfig snp_identifier %r shown above is "
            "applied only to inputs that have no recorded provenance.",
            "; ".join(sources),
            default_mode,
        )
        return

    # The GlobalConfig snp_identifier in the banner above is only a fallback. State
    # plainly whether it was used at all, so it is never mistaken for the mode in force.
    if used_default:
        default_clause = (
            f"The GlobalConfig snp_identifier {default_mode!r} shown above was used only as a "
            f"fallback for inputs without recorded provenance"
        )
    else:
        default_clause = (
            f"The GlobalConfig snp_identifier {default_mode!r} shown above was NOT used; "
            f"this mode comes entirely from input provenance"
        )
    downgrade_clause = (
        " An allele-awareness downgrade was applied to reach this shared base mode."
        if downgrade
        else ""
    )
    LOGGER.info(
        "SNP identifier used for this regression: %s.%s %s. Per-input modes: %s.",
        effective,
        downgrade_clause,
        default_clause,
        "; ".join(sources),
    )


def _add_jackknife_model_context(error, merged, columns, intercept):
    """Attach filtered-row support and genomic spans to failed deletions."""
    names = [*columns, *(["intercept"] if intercept is None else [])]
    support = np.count_nonzero(merged[columns].to_numpy(), axis=0)
    for failure in error.failures:
        start, end = failure["row_start"], failure["row_end"]
        deleted = merged.iloc[start:end]
        remaining = support - np.count_nonzero(deleted[columns].to_numpy(), axis=0)
        failure["columns"] = names
        failure["zero_columns"] = [names[i] for i in failure["zero_column_indices"]]
        failure["support_counts"] = {
            name: {"full": int(total), "remaining": int(left)}
            for name, total, left in zip(columns, support, remaining)
        }
        chr_column = "_ldsc_panel_CHR" if "_ldsc_panel_CHR" in deleted else "CHR"
        pos_column = "_ldsc_panel_POS" if "_ldsc_panel_POS" in deleted else "POS"
        if {chr_column, pos_column}.issubset(deleted.columns):
            failure["genomic_spans"] = [
                {"CHR": str(chrom), "POS_start": int(rows[pos_column].min()), "POS_end": int(rows[pos_column].max())}
                for chrom, rows in deleted.dropna(subset=[pos_column]).groupby(chr_column, sort=False)
            ]


def _raise_on_model_collinearity(dataset: RegressionDataset, x: np.ndarray) -> None:
    """Abort the fit when the LD-score design matrix is near-collinear.

    Matches legacy LDSC's hard stop on an ill-conditioned design matrix: rather
    than warn and emit a partition the user must second-guess, collinear
    annotations raise ``LDSCInputError`` so the only results produced are
    trustworthy ones.
    """
    overlap = getattr(dataset, "ldscore_overlap", None)
    if overlap is None or len(dataset.retained_ld_columns) < 2:
        return
    from .overlap_matrix import assemble_model_overlap, model_collinearity_error

    use_common = dataset.count_key_used_for_regression == COMMON_COUNT_KEY
    try:
        model_overlap = assemble_model_overlap(overlap, list(dataset.retained_ld_columns), use_common)
    except (ValueError, KeyError):
        return
    message = model_collinearity_error(x, list(dataset.retained_ld_columns), model_overlap)
    if message:
        raise LDSCInputError(message)


def _sort_partitioned_h2_summary(summary: pd.DataFrame, sort_by: str = "category") -> pd.DataFrame:
    """Return ``summary`` ordered by a public partitioned-h2 summary column."""
    if sort_by == "category":
        return summary.reset_index(drop=True)
    column = PARTITIONED_H2_SUMMARY_SORT_COLUMNS.get(sort_by)
    if column is None:
        valid = ", ".join(PARTITIONED_H2_SUMMARY_SORT_COLUMNS)
        raise LDSCUsageError(
            f"partitioned-h2 summary sort key {sort_by!r} is not supported. "
            f"Choose one of: {valid}."
        )
    if column not in summary.columns:
        raise LDSCInternalError(
            f"partitioned-h2 summary cannot be sorted by {sort_by!r} because required column {column!r} "
            "is missing. Most likely the regression workflow returned a result table with the wrong schema. "
            "Re-run with DEBUG logging and report the traceback."
        )
    return summary.sort_values(
        by=column,
        ascending=sort_by in PARTITIONED_H2_SUMMARY_ASCENDING_SORTS,
        na_position="last",
        kind="mergesort",
    ).reset_index(drop=True)


def _select_intercept(value: float | None, use_intercept: bool, default_when_disabled: float):
    """Resolve fixed-intercept config into the scalar expected by the kernel."""
    if not use_intercept:
        return default_when_disabled
    if value is None:
        return None
    return float(value)


def _resolve_default_chisq_max(chisq_max: float | None, n_annot: int, n_max: float) -> float | None:
    """Effective chi-square ceiling for an h2 fit, matching legacy ``estimate_h2``.

    Legacy LDSC (``ldscore/sumstats.py`` ``estimate_h2``) caps multi-annotation
    (partitioned) models at ``max(0.001 * N.max(), 80)`` when ``--chisq-max`` is
    unset, to down-weight outlier SNPs whose extreme chi-square would otherwise
    dominate the regression. Single-annotation models stay uncapped and rely on
    the two-step estimator instead. An explicit ``chisq_max`` always takes
    precedence over the default.
    """
    if chisq_max is not None:
        return float(chisq_max)
    if n_annot > 1:
        return max(0.001 * float(n_max), 80.0)
    return None


def _add_scalar_prevalence_arguments(parser) -> None:
    """Add paired case fractions for optional liability-scale estimates."""
    parser.add_argument(
        '--samp-prev', type=float, default=None, metavar='VALUE',
        help=(
            'Sample case fraction, strictly between 0 and 1, for liability-scale estimates. Requires '
            '--pop-prev; omit both to report observed-scale estimates only.'
        ),
    )
    parser.add_argument(
        '--pop-prev', type=float, default=None, metavar='VALUE',
        help=(
            'Population case fraction, strictly between 0 and 1, for liability-scale estimates. Requires '
            '--samp-prev; omit both to report observed-scale estimates only.'
        ),
    )


def add_h2_arguments(parser) -> None:
    """Register grouped h2 command-line options."""
    parser.prog = "ldsc h2"
    parser.formatter_class = CLIHelpFormatter
    inputs, model, scale, advanced = _add_common_regression_arguments(parser, include_h2_intercept=True)
    inputs.add_argument(
        '--trait-name', default=None,
        help=(
            'Trait label used in summaries. If omitted, use the stored trait name or derive one from the '
            'input filename.'
        ),
    )
    _add_scalar_prevalence_arguments(scale)
    runtime = parser.add_argument_group("Output, performance, and logging")
    runtime.add_argument(
        '--overwrite', action='store_true', default=False,
        help=(
            "Replace this command's existing output files and remove obsolete outputs from an earlier run. "
            'Default: off; stop if output files already exist.'
        ),
    )
    runtime.add_argument(
        '--log-level', default='INFO', choices=('DEBUG', 'INFO', 'WARNING', 'ERROR'),
        help=LOG_LEVEL_HELP,
    )


def add_partitioned_h2_arguments(parser) -> None:
    """Register grouped partitioned-h2 command-line options."""
    parser.prog = "ldsc partitioned-h2"
    parser.formatter_class = CLIHelpFormatter
    inputs, model, scale, advanced = _add_common_regression_arguments(parser, include_h2_intercept=True)
    inputs.add_argument(
        '--trait-name', default=None,
        help=(
            'Trait label used in summaries. If omitted, use the stored trait name or derive one from the '
            'input filename.'
        ),
    )
    _add_scalar_prevalence_arguments(scale)
    model.add_argument(
        '--summary-sort-by', default='auto', choices=('auto', *PARTITIONED_H2_SUMMARY_SORT_COLUMNS),
        help=(
            'Order partitioned_h2.tsv by the selected column. Default: auto, coefficient-p for per-query '
            'fits and category for baseline-only fits. P-values sort ascending; other numeric values sort '
            'descending.'
        ),
    )
    runtime = parser.add_argument_group("Output, performance, and logging")
    runtime.add_argument(
        '--threads', type=_parse_threads, default=1, metavar='N',
        help=(
            'Whole-query worker processes: 1 runs inline (default); positive N requests N workers; '
            '-1 uses available cores and -2 leaves one free; zero is invalid. '
            'Same worker-count rule as ldscore and build-gene-ldscore-index: negative values use '
            'CPU affinity with a machine CPU-count fallback; all values are capped by query count '
            'and query batch size. Choose positive N within your CPU allocation. '
            'Parallel workers use one native numerical thread each; inline execution keeps caller settings.'
        ),
    )
    runtime.add_argument(
        '--query-batch-size', type=int, default=1000, metavar='N',
        help=(
            'Maximum number of query annotations loaded at once; each query is fitted separately with the '
            'shared baseline. Default: 1000; a smaller positive value reduces memory and caps active query workers.'
        ),
    )
    runtime.add_argument(
        '--continue-on-query-error', action='store_true', default=False,
        help=(
            'Skip and record queries whose model preparation, regression, or result calculation fails; '
            'publish successful fits and write diagnostics/query_status.tsv plus error tracebacks in the log. '
            'Default: off; any query error prevents publication. Shared input/output failures and runs '
            'with no successful queries still fail.'
        ),
    )
    runtime.add_argument(
        '--overwrite', action='store_true', default=False,
        help=(
            "Replace this command's existing output files and remove obsolete outputs from an earlier run. "
            'Default: off; stop if output files already exist.'
        ),
    )
    runtime.add_argument(
        '--log-level', default='INFO', choices=('DEBUG', 'INFO', 'WARNING', 'ERROR'),
        help=LOG_LEVEL_HELP,
    )


def add_rg_arguments(parser) -> None:
    """Register grouped rg command-line options."""
    parser.prog = "ldsc rg"
    parser.formatter_class = CLIHelpFormatter
    inputs, model, scale, advanced = _add_common_regression_arguments(parser, include_h2_intercept=False)
    model.add_argument(
        '--anchor-trait', default=None,
        help=(
            'Trait name or input file path to compare against every other trait. Must identify one of '
            '--sumstats-sources. If omitted, estimate every unique trait pair.'
        ),
    )
    scale.add_argument(
        '--samp-prev', default=None,
        help=(
            'Comma-separated sample case fractions in resolved --sumstats-sources order, one per trait; '
            'each in (0, 1) or nan for a quantitative trait. Requires --pop-prev; cannot be combined with '
            '--prevalence-manifest. If neither form is supplied, report observed-scale estimates only.'
        ),
    )
    scale.add_argument(
        '--pop-prev', default=None,
        help=(
            'Comma-separated population case fractions in resolved --sumstats-sources order, one per trait; '
            'each in (0, 1) or nan. Requires --samp-prev; cannot be combined with --prevalence-manifest. If '
            'neither form is supplied, report observed-scale estimates only.'
        ),
    )
    scale.add_argument(
        '--prevalence-manifest', default=None,
        help=(
            'TSV with trait_name, samp_prev, and pop_prev columns, matched to exact input trait names. '
            'Cannot be combined with --samp-prev or --pop-prev. If omitted, use paired prevalence lists '
            'when supplied; otherwise report observed-scale estimates only.'
        ),
    )
    runtime = parser.add_argument_group("Output, performance, and logging")
    runtime.add_argument(
        '--write-per-pair-detail', action='store_true', default=False,
        help=(
            'Also write detailed results for each tested pair under diagnostics/pairs/. Default: off; write '
            'aggregate pair and trait tables only.'
        ),
    )
    runtime.add_argument(
        '--overwrite', action='store_true', default=False,
        help=(
            "Replace this command's existing output files and remove obsolete outputs from an earlier run. "
            'Default: off; stop if output files already exist.'
        ),
    )
    runtime.add_argument(
        '--log-level', default='INFO', choices=('DEBUG', 'INFO', 'WARNING', 'ERROR'),
        help=LOG_LEVEL_HELP,
    )


@materializing_overwrite_guard(
    lambda args: (
        (getattr(args, "output_dir"), getattr(args, "overwrite", False), "RUN_FAILED.txt")
        if getattr(args, "output_dir", None)
        else None
    ),
    command="run_h2_from_args(...)",
)
def run_h2_from_args(args):
    """Run single-trait heritability estimation from parsed CLI arguments.

    The workflow requires ``args.output_dir`` and preflights ``h2.tsv``, the
    LD Score regression-bin diagnostic, metadata, and ``diagnostics/h2.log``
    before loading inputs. A successful overwrite removes default plots and
    liability-scale post-processing derived from the superseded h2 result.
    """
    output_dir, log_path = _preflight_regression_outputs(args, "h2", H2DirectoryWriter)
    with workflow_logging("h2", log_path, log_level=getattr(args, "log_level", "INFO")):
        from ._input_preflight import inspect_declared_inputs, inspect_artifact_paths
        inspect_declared_inputs(
            files=[('sumstats', getattr(args, 'sumstats_file', None) or getattr(args, 'sumstats_sources', None)),
                   ('prevalence manifest', getattr(args, 'prevalence_manifest', None))],
            checks=[('LD scores', args.ldscore_dir, lambda: inspect_artifact_paths(args.ldscore_dir))],
            issues_path=Path(output_dir)/'diagnostics/input_issues.tsv')
        runner, config = _runner_from_args(args)
        print_global_config_banner("run_h2_from_args", runner.global_config)
        log_inputs(sumstats_file=args.sumstats_file, ldscore_dir=args.ldscore_dir, output_dir=output_dir)
        LOGGER.info(f"Starting h2 regression for '{args.sumstats_file}' using LD-score directory '{args.ldscore_dir}'.")
        sumstats_table = _load_sumstats_table(args.sumstats_file, getattr(args, "trait_name", None))
        ldscore_result = load_ldscore_from_dir(args.ldscore_dir)
        legacy_drops = pd.DataFrame(columns=LEGACY_SUMSTATS_DROP_COLUMNS)
        if _is_legacy_sumstats(sumstats_table):
            sumstats_table, legacy_drops = _project_legacy_sumstats_to_panel(sumstats_table, ldscore_result)
        _log_effective_regression_identity([sumstats_table], ldscore_result, runner.global_config, config)
        with suppress_global_config_banner():
            dataset = runner.build_dataset(sumstats_table, ldscore_result, config=config)
        if not legacy_drops.empty:
            dataset = replace(dataset, legacy_sumstats_drops=legacy_drops)
        outcome = runner._fit_h2_dataset(dataset, config=config)
        hsq, dataset = outcome.estimator, outcome.dataset
        n_snps_used = outcome.n_snps
        summary = summarize_total_h2(
            hsq, dataset, trait_name=sumstats_table.trait_name, n_snps_used=n_snps_used,
            samp_prev=config.samp_prev, pop_prev=config.pop_prev,
        )
        output_dir_arg = getattr(args, "output_dir", None)
        if output_dir_arg:
            written = H2DirectoryWriter().write(
                summary,
                H2OutputConfig(
                    output_dir=output_dir,
                    overwrite=getattr(args, "overwrite", False),
                ),
                metadata=_h2_metadata(
                    args,
                    sumstats_table,
                    outcome,
                    samp_prev=config.samp_prev,
                    pop_prev=config.pop_prev,
                ),
                diagnostic_bins=outcome.diagnostic_bins,
            )
            audit_path = _write_or_remove_legacy_sumstats_audit(
                Path(output_dir),
                legacy_used=_has_legacy_sumstats_source(sumstats_table),
                drops=legacy_drops,
                overwrite=getattr(args, "overwrite", False),
            )
            if audit_path is not None:
                written["legacy_sumstats_drops"] = str(audit_path)
            log_outputs(**written)
        LOGGER.info(f"Finished h2 regression with {n_snps_used} regression SNPs.")
    return summary


@materializing_overwrite_guard(
    lambda args: (
        (getattr(args, "output_dir"), getattr(args, "overwrite", False), "RUN_FAILED.txt")
        if getattr(args, "output_dir", None)
        else None
    ),
    command="run_partitioned_h2_from_args(...)",
)
def run_partitioned_h2_from_args(args):
    """Run batch partitioned heritability from parsed CLI arguments.

    The workflow requires ``args.output_dir`` and preflights
    ``partitioned_h2.tsv``, ``diagnostics/query_annotations/``, and
    ``diagnostics/partitioned-h2.log`` before loading inputs. Query-annotation
    runs always write per-query results; baseline-only runs keep the complete
    fitted-model artifacts at the result root. A successful overwrite removes
    the default plot root derived from the superseded result.
    ``args.continue_on_query_error`` explicitly permits query regression
    failures: successful models are published and every attempted query is
    recorded in ``diagnostics/query_status.tsv``. Strict mode, the default,
    collects query errors and publishes only the diagnostic ledger on failure.
    ``args.threads`` uses the shared LD-score/index worker-count policy, with
    query count and ``args.query_batch_size`` as the work limits. The default
    and every effective count of one run inline. Parallel processes each fit
    one complete query model with one native numerical thread; baseline-only
    input remains one inline joint fit. Worker death and transport/output
    failures abort under either query-error policy. Requested/effective counts
    and phase timings are recorded in the workflow log.
    """
    _validate_threads(getattr(args, "threads", 1))
    query_batch_size = getattr(args, "query_batch_size", 1000)
    if isinstance(query_batch_size, bool) or not isinstance(query_batch_size, int) or query_batch_size < 1:
        raise ValueError("query_batch_size must be a positive integer.")
    output_dir, log_path = _preflight_regression_outputs(args, "partitioned-h2", PartitionedH2DirectoryWriter)
    with workflow_logging("partitioned-h2", log_path, log_level=getattr(args, "log_level", "INFO")):
        from ._input_preflight import inspect_declared_inputs, inspect_artifact_paths
        inspect_declared_inputs(
            files=[('sumstats', getattr(args, 'sumstats_file', None) or getattr(args, 'sumstats_sources', None)),
                   ('prevalence manifest', getattr(args, 'prevalence_manifest', None))],
            checks=[('LD scores', args.ldscore_dir, lambda: inspect_artifact_paths(args.ldscore_dir))],
            issues_path=Path(output_dir)/'diagnostics/input_issues.tsv')
        runner, config = _runner_from_args(args)
        print_global_config_banner("run_partitioned_h2_from_args", runner.global_config)
        log_inputs(sumstats_file=args.sumstats_file, ldscore_dir=args.ldscore_dir, output_dir=output_dir)
        LOGGER.info(
            f"Starting partitioned-h2 regression for '{args.sumstats_file}' using LD-score directory '{args.ldscore_dir}'."
        )
        with _log_phase_timing("sumstats loading"):
            sumstats_table = _load_sumstats_table(args.sumstats_file, getattr(args, "trait_name", None))
        with _log_phase_timing("LD-score loading"):
            ldscore_result = load_ldscore_from_dir(args.ldscore_dir)
        legacy_drops = pd.DataFrame(columns=LEGACY_SUMSTATS_DROP_COLUMNS)
        with _log_phase_timing("legacy projection"):
            if _is_legacy_sumstats(sumstats_table):
                sumstats_table, legacy_drops = _project_legacy_sumstats_to_panel(sumstats_table, ldscore_result)
        has_queries = bool(ldscore_result.query_columns)
        _log_effective_regression_identity([sumstats_table], ldscore_result, runner.global_config, config)
        _log_partitioned_h2_regime(ldscore_result, has_queries)
        _log_quantitative_annotation_interpretation(ldscore_result)
        with suppress_global_config_banner():
            result = runner.estimate_partitioned_h2_batch(
                sumstats_table, ldscore_result, config=config, output_dir=output_dir,
                query_batch_size=getattr(args, "query_batch_size", 1000),
                threads=getattr(args, "threads", 1),
                continue_on_query_error=getattr(args, "continue_on_query_error", False),
                overwrite=getattr(args, "overwrite", False),
                summary_sort_by=getattr(args, "summary_sort_by", "auto"),
                metadata={
                    "trait_name": sumstats_table.trait_name,
                    "count_kind": getattr(args, "count_kind", "common"),
                    "ldscore_dir": getattr(args, "ldscore_dir", None),
                    "analysis_type": "cell_type_specific" if has_queries else "functional_category",
                    "headline_metric": "coefficient" if has_queries else "enrichment",
                    "enrichment_p_test": "two_sided_t", "coefficient_p_test": "one_sided_greater",
                    "annotation_types": dict(ldscore_result.annotation_types),
                },
            )
        summary, written = result.summary, dict(result.output_paths)
        audit_path = _write_or_remove_legacy_sumstats_audit(
            Path(output_dir), legacy_used=_has_legacy_sumstats_source(sumstats_table),
            drops=legacy_drops, overwrite=getattr(args, "overwrite", False),
        )
        if audit_path is not None:
            written["legacy_sumstats_drops"] = str(audit_path)
        log_outputs(**written)
        LOGGER.info(
            f"Finished partitioned-h2 regression for {len(ldscore_result.query_columns)} query annotations "
            f"and {len(summary)} summary rows."
        )
    return summary


@materializing_overwrite_guard(
    lambda args: (
        (getattr(args, "output_dir"), getattr(args, "overwrite", False), "RUN_FAILED.txt")
        if getattr(args, "output_dir", None)
        else None
    ),
    command="run_rg_from_args(...)",
)
def run_rg_from_args(args):
    """Run multi-trait genetic-correlation estimation from parsed CLI args.

    The workflow requires ``--output-dir`` and writes the rg result family:
    ``rg.tsv``, ``rg_full.tsv``, ``h2_per_trait.tsv``, optional
    ``diagnostics/pairs/``, and workflow-owned ``diagnostics/rg.log``. A
    successful overwrite removes the default plot root derived from the
    superseded result.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed rg options. Required fields are ``sumstats_sources`` and
        ``ldscore_dir`` and ``output_dir``. Optional fields include ``anchor_trait``,
        ``write_per_pair_detail``, intercept settings, and
        common regression options.

    Returns
    -------
    RgResultFamily
        Complete in-memory result family. The function itself never prints.

    Raises
    ------
    ValueError
        If fewer than two sumstats paths resolve or ``output_dir`` is absent.
    """
    from .prevalence import resolve_rg_prevalences

    output_dir, log_path = _preflight_regression_outputs(args, "rg", RgDirectoryWriter)
    with workflow_logging("rg", log_path, log_level=getattr(args, "log_level", "INFO")):
        from ._input_preflight import inspect_declared_inputs, inspect_artifact_paths
        declared_inputs = inspect_declared_inputs(
            files=[('sumstats', getattr(args, 'sumstats_file', None) or getattr(args, 'sumstats_sources', None)),
                   ('prevalence manifest', getattr(args, 'prevalence_manifest', None))],
            checks=[('LD scores', args.ldscore_dir, lambda: inspect_artifact_paths(args.ldscore_dir))],
            issues_path=Path(output_dir)/'diagnostics/input_issues.tsv')
        sumstats_paths = declared_inputs.files['sumstats']
        if len(sumstats_paths) < 2:
            raise LDSCUserError(
                "rg requires at least two sumstats inputs. Most likely `--sumstats-sources` "
                "resolved to fewer than two files. Pass two or more munged `.sumstats` artifacts."
            )
        runner, config = _runner_from_args(args)
        print_global_config_banner("run_rg_from_args", runner.global_config)
        log_inputs(
            sumstats_sources=[str(path) for path in sumstats_paths],
            anchor_trait=getattr(args, "anchor_trait", None) or "none",
            ldscore_dir=args.ldscore_dir,
            output_dir=output_dir,
        )
        LOGGER.info(
            f"Starting rg regression for {len(sumstats_paths)} sumstats files "
            f"using LD-score directory '{args.ldscore_dir}'."
        )
        from ._input_preflight import InputGate
        from ._progress import PhaseProgress
        gate = InputGate('Summary-statistic content validation', Path(output_dir)/'diagnostics/input_issues.tsv')
        sumstats_tables = []
        with PhaseProgress(LOGGER, 'validation', 'summary-statistic content', len(sumstats_paths)) as progress:
            for path in sumstats_paths:
                table = gate.check('sumstats', path, lambda path=path: _load_sumstats_table(str(path), None))
                if table is not None:
                    sumstats_tables.append(table)
                progress.advance(object=str(path))
            gate.finish()
        # Resolve prevalence against the original munged names so the manifest
        # duplicate-name guard fires before disambiguation renames collisions.
        prevalences = resolve_rg_prevalences(
            samp_prev=getattr(args, "samp_prev", None),
            pop_prev=getattr(args, "pop_prev", None),
            manifest_path=getattr(args, "prevalence_manifest", None),
            trait_names=[_trait_label(table) for table in sumstats_tables],
            trait_paths=[str(path) for path in sumstats_paths],
        )
        sumstats_tables = _disambiguate_trait_names(sumstats_tables)
        if prevalences is not None:
            LOGGER.info(
                "Liability-scale prevalences applied: "
                + ", ".join(
                    f"{_trait_label(table)}=(P={sp}, K={pp})"
                    for table, (sp, pp) in zip(sumstats_tables, prevalences)
                )
            )
        anchor_index = _resolve_anchor_index(getattr(args, "anchor_trait", None), sumstats_paths, sumstats_tables)
        ldscore_result = load_ldscore_from_dir(args.ldscore_dir)
        legacy_drop_frames: list[pd.DataFrame] = []
        projected_tables: list[SumstatsTable] = []
        for table in sumstats_tables:
            if _is_legacy_sumstats(table):
                table, drops = _project_legacy_sumstats_to_panel(table, ldscore_result)
                legacy_drop_frames.append(drops)
            projected_tables.append(table)
        sumstats_tables = projected_tables
        legacy_drops = (
            pd.concat(legacy_drop_frames, ignore_index=True)
            if legacy_drop_frames
            else pd.DataFrame(columns=LEGACY_SUMSTATS_DROP_COLUMNS)
        )
        _log_effective_regression_identity(sumstats_tables, ldscore_result, runner.global_config, config)
        with suppress_global_config_banner():
            result = runner.estimate_rg_pairs(
                sumstats_tables,
                ldscore_result,
                anchor_index=anchor_index,
                config=config,
                prevalences=prevalences,
            )
        output_dir_arg = getattr(args, "output_dir", None)
        if output_dir_arg:
            written = RgDirectoryWriter().write(
                result,
                RgOutputConfig(
                    output_dir=output_dir,
                    overwrite=getattr(args, "overwrite", False),
                    write_per_pair_detail=getattr(args, "write_per_pair_detail", False),
                ),
            )
            audit_path = _write_or_remove_legacy_sumstats_audit(
                Path(output_dir),
                legacy_used=any(_has_legacy_sumstats_source(table) for table in sumstats_tables),
                drops=legacy_drops,
                overwrite=getattr(args, "overwrite", False),
            )
            if audit_path is not None:
                written["legacy_sumstats_drops"] = str(audit_path)
            log_outputs(**written)
        LOGGER.info(
            f"Finished rg regression for {len(result.rg)} trait pairs "
            f"and {len(result.h2_per_trait)} traits."
        )
    return result


def _add_common_regression_arguments(parser, include_h2_intercept: bool):
    """Add shared regression groups and return them for command-specific options."""
    parser.allow_abbrev = False
    inputs = parser.add_argument_group("Inputs and output")
    model = parser.add_argument_group("Model settings")
    scale = parser.add_argument_group("Liability scale (optional)")
    advanced = parser.add_argument_group("Advanced estimation and SNP matching")
    if include_h2_intercept:
        inputs.add_argument(
            '--sumstats-file', required=True, metavar='FILE',
            help=(
            'Required cleaned summary-statistics Parquet or LDSC2 .sumstats(.gz) file. SNP identity must be '
            'compatible with --ldscore-dir. '
            + SCALAR_PATH_HELP
        ),
        )
    else:
        inputs.add_argument(
            '--sumstats-sources', nargs='+', required=True, metavar='SOURCES',
            help=(
            "Required cleaned Parquet or LDSC2 .sumstats(.gz) files; accepts exact paths or quoted '*' patterns "
            'matching filename text and requires at least two files. SNP identity must be compatible with '
            "--ldscore-dir. '@' is not expanded."
        ),
        )
    inputs.add_argument(
        '--ldscore-dir', required=True, metavar='DIR',
        help=(
            'Required LD-score result directory from ldscore or convert-ldsc2-ldscores, containing '
            'predictors, regression weights, and SNP counts.'
        ),
    )
    inputs.add_argument(
        '--output-dir', required=True, metavar='DIR',
        help=(
            'Required destination for regression estimates and diagnostics.'
        ),
    )
    model.add_argument(
        '--count-kind', choices=('common', 'all'), default='common',
        help=(
            'Reference SNP counts used to scale regression estimates: common uses common-variant counts, '
            'all uses all-reference counts. Default: common, falling back to all counts when common counts '
            'are unavailable. Requested counts must be present.'
        ),
    )
    if include_h2_intercept:
        advanced.add_argument(
            '--intercept-h2', type=float, default=None, metavar='VALUE',
            help=(
                'Fix the heritability intercept to this value; use 1 for the standard fixed intercept. '
                'Cannot be combined with an explicit --two-step-cutoff. If omitted, estimate the intercept.'
            ),
        )
    else:
        advanced.add_argument(
            '--intercept-h2', type=float, default=None, metavar='VALUE',
            help=(
                'Fix the heritability intercept to the same value for both traits in every pair; use 1 for '
                'the standard fixed intercept. Cannot be combined with an explicit --two-step-cutoff. '
                'If omitted, estimate the intercepts.'
            ),
        )
        advanced.add_argument(
            '--intercept-gencov', type=float, default=None, metavar='VALUE',
            help=(
                'Fix the genetic-covariance intercept to the same value for every trait pair; use 0 for '
                'the standard fixed intercept. Cannot be combined with --two-step-cutoff. Single-annotation '
                'fits also require --intercept-h2 because automatic two-step estimation otherwise conflicts '
                'with this fixed value. If omitted, estimate the covariance intercepts.'
            ),
        )
    if include_h2_intercept:
        advanced.add_argument(
            '--two-step-cutoff', type=float, default=None, metavar='VALUE',
            help=(
                'First-step chi-square cutoff, including the boundary. Requires a single-annotation model with '
                'a free intercept; cannot be combined with --intercept-h2. If omitted, use 30 '
                'for such models; otherwise disable two-step estimation.'
            ),
        )
        advanced.add_argument(
            '--chisq-max', type=float, default=None, metavar='VALUE',
            help=(
                'Keep SNPs with chi-square at or below this value. If omitted, use max(0.001 * maximum N, 80) '
                'for multi-annotation models; single-annotation models have no chi-square cap.'
            ),
        )
    else:
        advanced.add_argument(
            '--two-step-cutoff', type=float, default=None, metavar='VALUE',
            help=(
                'First-step chi-square cutoff, including the boundary. Requires a single-annotation model with '
                'free intercepts; cannot be combined with --intercept-h2 or '
                '--intercept-gencov. If omitted, use 30 with a free h2 intercept in single-annotation fits; '
                'otherwise disable it.'
            ),
        )
        advanced.add_argument(
            '--chisq-max', type=float, default=None, metavar='VALUE',
            help=(
                'Keep SNPs satisfying Z1^2 * Z2^2 <= VALUE^2 for each trait pair. If omitted, apply no '
                'chi-square product cap. This differs from the single-trait chi-square filter.'
            ),
        )
    advanced.add_argument(
        '--n-blocks', type=int, default=200, metavar='N',
        help=(
            'Number of SNP blocks for jackknife standard errors; must be at least 2. Default: 200, capped '
            'at the number of retained SNPs.'
        ),
    )
    advanced.add_argument(
        '--allow-identity-downgrade', action='store_true', default=False,
        help=(
            'Allow allele-aware and allele-unaware inputs in the same identity family by matching without '
            'alleles. rsid and chr_pos families still cannot be mixed. Default: off; require compatible '
            'identity modes.'
        ),
    )
    return inputs, model, scale, advanced


def _preflight_regression_outputs(
    args,
    workflow_name: str,
    writer,
) -> tuple[str, Path]:
    """Preflight regression outputs and return normalized output dir plus log path."""
    output_dir_arg = getattr(args, "output_dir", None)
    if not output_dir_arg:
        raise LDSCUsageError(
            f"{workflow_name} requires `--output-dir`. Most likely a programmatic namespace omitted output_dir. "
            "Pass an explicit result directory."
        )
    output_dir = ensure_output_directory(output_dir_arg, label="output directory")
    log_path = output_dir / "diagnostics" / f"{workflow_name}.log"
    writer.artifact_family(output_dir).preflight(
        overwrite=getattr(args, "overwrite", False),
        additional_paths=[log_path, _legacy_sumstats_audit_path(output_dir)],
    )
    log_path.parent.mkdir(parents=True, exist_ok=True)
    return str(output_dir), log_path


def _legacy_sumstats_audit_path(output_dir: Path) -> Path:
    """Declare the workflow-owned legacy projection audit shared by h2 and rg."""
    return output_dir / "diagnostics" / "dropped_snps" / "legacy_sumstats.tsv.gz"


def _write_or_remove_legacy_sumstats_audit(
    output_dir: Path,
    *,
    legacy_used: bool,
    drops: pd.DataFrame,
    overwrite: bool,
) -> Path | None:
    """Write the stable legacy projection audit, or remove a stale owned copy."""
    path = _legacy_sumstats_audit_path(output_dir)
    if not legacy_used:
        if overwrite and path.exists():
            path.unlink()
        return None
    path.parent.mkdir(parents=True, exist_ok=True)
    drops.reindex(columns=LEGACY_SUMSTATS_DROP_COLUMNS).to_csv(
        path,
        sep="\t",
        index=False,
        na_rep="",
        compression={"method": "gzip", "mtime": 0},
    )
    return path


def _disambiguate_trait_names(tables: Sequence[SumstatsTable]) -> list[SumstatsTable]:
    """Return tables with deterministic unique trait names for rg outputs."""
    labels = [_trait_label(table) for table in tables]
    duplicates = {label for label in labels if labels.count(label) > 1}
    if not duplicates:
        return [replace(table, trait_name=label) for table, label in zip(tables, labels)]

    proposed: list[str] = []
    for table, label in zip(tables, labels):
        if label not in duplicates:
            proposed.append(label)
            continue
        source = Path(table.source_path) if table.source_path else Path(label)
        parent = source.parent.name or "source"
        proposed.append(f"{label}@{parent}")

    seen: dict[str, int] = {}
    unique: list[str] = []
    for table, label in zip(tables, proposed):
        if label not in seen:
            seen[label] = 1
            unique.append(label)
            continue
        source_hash = hashlib.sha1(_source_key(table).encode("utf-8")).hexdigest()[:8]
        final = f"{label}@{source_hash}"
        while final in seen:
            source_hash = hashlib.sha1(f"{_source_key(table)}:{seen[label]}".encode("utf-8")).hexdigest()[:8]
            final = f"{label}@{source_hash}"
        seen[label] += 1
        seen[final] = 1
        unique.append(final)

    for original, resolved in zip(labels, unique):
        if original != resolved:
            LOGGER.info(f"Disambiguated duplicate rg trait name '{original}' as '{resolved}'.")
    return [replace(table, trait_name=label) for table, label in zip(tables, unique)]


def _resolve_anchor_index(
    anchor_trait: str | None,
    sumstats_paths: Sequence[str],
    sumstats_tables: Sequence[SumstatsTable],
) -> int | None:
    """Resolve ``--anchor-trait`` by trait label first, then source path."""
    if not anchor_trait:
        return None
    token = normalize_path_token(anchor_trait)
    trait_matches = {idx for idx, table in enumerate(sumstats_tables) if table.trait_name == token}
    if len(trait_matches) == 1:
        return next(iter(trait_matches))
    if len(trait_matches) > 1:
        available = [table.trait_name or Path(path).name for table, path in zip(sumstats_tables, sumstats_paths)]
        raise LDSCUsageError(
            f"rg could not resolve `--anchor-trait={anchor_trait}` because it matched "
            f"{len(trait_matches)} trait names. Most likely multiple inputs have the same trait label. "
            f"Use a unique input path or one of these resolved trait names: {available}."
        )

    path_matches: set[int] = set()
    try:
        anchor_path = Path(token).resolve(strict=False)
        for idx, path in enumerate(sumstats_paths):
            if Path(path).resolve(strict=False) == anchor_path:
                path_matches.add(idx)
    except OSError:
        pass  # The actionable unknown-anchor error below includes valid inputs.
    if len(path_matches) != 1:
        available = [table.trait_name or Path(path).name for table, path in zip(sumstats_tables, sumstats_paths)]
        raise LDSCUsageError(
            f"rg could not resolve `--anchor-trait={anchor_trait}` because it matched "
            f"{len(path_matches)} input paths. Most likely the anchor token is misspelled or not part "
            f"of `--sumstats-sources`. Use one of these trait names or exact input paths: {available}."
        )
    return next(iter(path_matches))


def _runner_from_args(args) -> tuple[RegressionRunner, RegressionConfig]:
    """Build the regression workflow objects from parsed CLI arguments.

    For ``h2`` / ``partitioned-h2`` the scalar ``--samp-prev`` / ``--pop-prev``
    floats are validated and stored on the config. For ``rg`` these arguments are
    comma-separated strings resolved later in :func:`run_rg_from_args` against the
    resolved trait list, so they are left off the (per-run) config here.
    """
    from .prevalence import parse_scalar_prevalence

    count_kind = getattr(args, "count_kind", "common")
    samp_prev = getattr(args, "samp_prev", None)
    pop_prev = getattr(args, "pop_prev", None)
    if isinstance(samp_prev, str) or isinstance(pop_prev, str):
        scalar_samp, scalar_pop = None, None  # rg comma-lists: resolved in run_rg_from_args
    else:
        scalar_samp, scalar_pop = parse_scalar_prevalence(samp_prev, pop_prev)
    config = RegressionConfig(
        n_blocks=args.n_blocks,
        use_common_counts=(count_kind == "common"),
        intercept_h2=args.intercept_h2,
        intercept_gencov=getattr(args, "intercept_gencov", None),
        two_step_cutoff=args.two_step_cutoff,
        chisq_max=args.chisq_max,
        samp_prev=scalar_samp,
        pop_prev=scalar_pop,
        allow_identity_downgrade=getattr(args, "allow_identity_downgrade", False),
    )
    runner = RegressionRunner(get_global_config(), config)
    return runner, config


def _load_sumstats_table(path: str, trait_name: str | None) -> SumstatsTable:
    """Load one curated sumstats artifact through the public workflow helper."""
    return load_sumstats(path, trait_name=trait_name)


def load_ldscore_from_dir(
    ldscore_dir: str,
    snp_identifier: str | None = None,
) -> LDScoreSource:
    """Open a canonical LD-score directory with query values read on demand.

    Parameters
    ----------
    ldscore_dir : str
        Directory containing ``metadata.json``, ``ldscore.baseline.parquet``,
        and the query files declared by the current ``query_batches`` manifest.
        Older directories without this manifest must be regenerated.
    snp_identifier : {"rsid", "rsid_allele_aware", "chr_pos", "chr_pos_allele_aware"} or None, optional
        Identifier mode used to reconstruct the public regression SNP set from
        the baseline table. When omitted, the metadata value is used. An
        explicit value must agree with the saved identity configuration.

    Returns
    -------
    LDScoreSource
        Shared baseline values, SNP metadata, counts, and overlap statistics.
        Call ``read_queries(columns)`` for explicit query selections, including
        selections spanning saved files. Reads preserve requested column order,
        are uncached, and have no width limit. The caller owns their memory
        cost. The original directory must remain available and unchanged.

    Raises
    ------
    LDSCInputError
        The directory or root metadata is missing, or the manifest, identity,
        count, or overlap contract is invalid.
    LDSCInternalError
        A declared baseline/query schema lacks required columns.
    ConfigMismatchError
        An explicit identity mode disagrees with the saved configuration.

    Notes
    -----
    Opening checks all declared schemas and allele metadata without loading
    query LD-score values. Selected-query numerical checks and alignment to
    baseline rows are performed by the consuming regression workflow.
    """
    root = Path(normalize_path_token(ldscore_dir))
    if not root.is_dir():
        raise LDSCInputError(
            f"Regression could not load LD-score directory '{root}': path is not an existing directory. "
            "Most likely `--ldscore-dir` points at the wrong location. Pass the directory written by `ldsc ldscore`."
        )
    metadata_path = root / "metadata.json"
    if not metadata_path.exists():
        raise LDSCInputError(
            f"Regression could not load LD-score directory '{root}': missing `metadata.json`. "
            "Most likely this is not a canonical `ldsc ldscore` output directory or the diagnostics file "
            "was copied without the data files. Pass the complete LD-score output directory."
        )
    LOGGER.info(f"Loading LD-score result directory from '{root}'.")
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    config_snapshot = _global_config_from_metadata(metadata)
    files = metadata.get("files", {})
    baseline_rel = files.get("baseline")
    if not baseline_rel:
        raise LDSCInputError(
            f"Regression could not load LD-score directory '{root}': metadata is missing `files.baseline`. "
            "Most likely the LD-score artifact was written by an older version or edited by hand. "
            f"Regenerate it with the current `ldsc ldscore`. Other causes & fixes: {_REGRESSION_SCHEMA_DOC}"
        )
    baseline_table = pd.read_parquet(root / baseline_rel)
    baseline_columns = [str(column) for column in metadata.get("baseline_columns", [])]
    query_columns = [str(column) for column in metadata.get("query_columns", [])]
    count_records = [dict(record) for record in metadata.get("counts", [])]
    from .ldscore_source import load_query_manifest
    query_batches = load_query_manifest(root, metadata, baseline_columns, baseline_table.columns)
    _validate_count_overlap_config(metadata, root)
    overlap = None
    overlap_rel = files.get("overlap")
    if overlap_rel:
        from .overlap_matrix import overlap_from_long_frame
        overlap_config = metadata.get("overlap_config") or {}
        overlap = overlap_from_long_frame(
            pd.read_parquet(root / overlap_rel),
            baseline_columns=baseline_columns,
            query_columns=query_columns,
            total_all_reference_snps=overlap_config.get("total_all_reference_snps"),
            total_common_reference_snps=overlap_config.get("total_common_reference_snps"),
        )
    metadata_identifier = config_snapshot.snp_identifier if config_snapshot is not None else metadata.get("snp_identifier")
    if snp_identifier is not None:
        requested_identifier = normalize_snp_identifier_mode(snp_identifier)
        if metadata_identifier is not None and requested_identifier != metadata_identifier:
            raise ConfigMismatchError(
                "snp_identifier mismatch when loading LD-score directory: "
                f"override {requested_identifier!r} vs metadata {metadata_identifier!r}. "
                "Use an LD-score artifact written for the requested SNP-identifier mode."
            )
        effective_identifier = requested_identifier
    else:
        effective_identifier = metadata_identifier or "chr_pos_allele_aware"
    _validate_ldscore_allele_columns(
        baseline_table,
        table_name="baseline_table",
        snp_identifier=effective_identifier,
    )
    query_table = None
    for entry in query_batches:
        identities = pd.read_parquet(entry["path"], columns=[
            name for name in ("CHR", "SNP", "POS", "A1", "A2") if name in entry["schema"]])
        _validate_ldscore_allele_columns(
            identities,
            table_name="query_table",
            snp_identifier=effective_identifier,
        )
        if query_table is None:
            query_table = identities
        del identities
    result = LDScoreSource(
        baseline_table=baseline_table,
        query_metadata=query_table,
        query_batches=query_batches,
        count_records=count_records,
        baseline_columns=baseline_columns,
        query_columns=query_columns,
        ld_reference_snps=frozenset(),
        ld_regression_snps=frozenset(build_snp_id_series(baseline_table, effective_identifier)),
        chromosome_results=[],
        output_paths={
            "metadata": str(metadata_path),
            "baseline": str(root / baseline_rel),
            **{key: str(root / value) for key, value in files.items()},
        },
        count_config=dict(metadata.get("count_config", {})),
        config_snapshot=config_snapshot,
        overlap=overlap,
        annotation_types={str(key): str(value) for key, value in (metadata.get("annotation_types") or {}).items()},
        snp_identifier=effective_identifier,
        chromosome_scope=metadata.get("chromosome_scope") or {},
        snp_universe_policy=metadata.get("snp_universe_policy"),
        index_provenance=({key: metadata[key] for key in ("index_id", "index_snp_identifier", "index_genome_build")}
                          if metadata.get("index_id") is not None else None),
        legacy_ldsc2_import=metadata.get("legacy_ldsc2_import"),
    )
    query_rows = 0 if query_table is None else len(query_table)
    LOGGER.info(
        f"Loaded LD-score directory '{root}' with {len(baseline_table)} baseline rows, "
        f"{query_rows} query rows, and config provenance {'present' if config_snapshot is not None else 'unknown'}."
    )
    return result


def _validate_count_overlap_config(metadata: dict[str, Any], root: Path) -> None:
    """Require count and overlap metadata to describe one common reference-SNP universe."""
    overlap_config = metadata.get("overlap_config")
    if not overlap_config:
        return
    count_config = metadata.get("count_config") or {}
    count_threshold = count_config.get("common_reference_snp_maf_min")
    overlap_threshold = overlap_config.get("common_maf_min")
    count_operator = count_config.get("common_reference_snp_maf_operator")
    overlap_operator = overlap_config.get("common_maf_operator")
    thresholds_agree = False
    try:
        thresholds_agree = float(count_threshold) == float(overlap_threshold)
    except (TypeError, ValueError):
        thresholds_agree = False
    if not thresholds_agree or count_operator != overlap_operator:
        raise LDSCInputError(
            f"Regression cannot load LD-score directory '{root}': count_config and overlap_config disagree "
            f"about common-SNP semantics (count threshold/operator {count_threshold!r}/{count_operator!r}, "
            f"overlap {overlap_threshold!r}/{overlap_operator!r}). Reconvert or regenerate the directory so "
            "counts and overlap use the same threshold and operator."
        )


def _global_config_from_metadata(metadata: dict[str, Any]) -> GlobalConfig | None:
    """Recreate a GlobalConfig snapshot when LD-score metadata contains one."""
    if "artifact_type" not in metadata:
        raise LDSCInputError(
            "Regression could not read LD-score artifact provenance: metadata lacks `artifact_type`. "
            "Most likely the LD-score directory was written by an older LDSC version. "
            f"Regenerate it with the current `ldsc ldscore`. Other causes & fixes: {_REGRESSION_SCHEMA_DOC}"
        )
    mode = validate_identity_artifact_metadata(metadata, expected_artifact_type="ldscore")
    return GlobalConfig(
        snp_identifier=mode,
        genome_build=metadata.get("genome_build"),
        log_level="INFO",
    )
