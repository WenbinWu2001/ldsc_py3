"""runner.py

Core functionality:
    Consume saved LD-score sources, assemble regression-ready datasets, and
    dispatch heritability or genetic-correlation estimators.

Overview
--------
This module is the workflow-layer boundary for the LDSC3 regression commands.
It consumes the LD-score directory introduced by the IO
refactor: ``metadata.json`` plus split baseline/query parquet tables. Public
callers pass one ``--ldscore-dir`` input instead of individual LD-score,
weight, count, and annotation-manifest files. Directory opening and artifact
validation belong to ``_ldscore.source.load_ldscore_from_dir``; this module
owns selected-value validation, trait alignment, and model construction.

Legacy ``.sumstats`` and ``.sumstats.gz`` inputs are accepted at this
boundary. Their rsIDs are lookup keys: each trait is projected onto the
regression SNPs in the LD-score directory, allele orientation is validated and harmonized, and
rejected rows are recorded in a stable diagnostic audit. Legacy LD-score
fragments are never read here and require explicit conversion first.

Regression chooses either the metadata ``common_reference_snp_counts`` vector or
``all_reference_snp_counts`` through ``--count-kind common|all`` (default
``common``). Requested counts must be present; the other kind is never substituted.
Regression merges by the effective key for the resolved SNP identifier mode.
Base modes ignore alleles; allele-aware/base mixes can run under the base mode
only when ``--allow-identity-downgrade`` or
``RegressionConfig.allow_identity_downgrade`` is set, and only within
rsID-family or within coordinate-family modes.
The ``partitioned-h2`` workflow fits the complete baseline model for baseline-only inputs and
one baseline-plus-query model per query annotation otherwise. Query runs share
trait/baseline preparation, load query columns in bounded batches, and fit
in whole-query workers via ``workers`` (default one). Completed
annotation-table/delete-value/metadata files are staged immediately. The output-layer
``PartitionedH2DirectoryWriter`` publishes them after final summary sorting.
Strict batches require every query to succeed. Explicit continuation records
failed queries and publishes successful models with an attempted-query ledger.
Genetic correlation accepts two or more munged summary-statistic sources and
returns the full ``rg`` output family:
the concise headline table, a diagnostic full table, per-trait SNP heritability summaries,
and optional per-pair metadata for filesystem detail outputs.

Regression commands require an ``output_dir`` and create per-run logs under
``diagnostics/``. Their Python forms ``run_h2``, ``run_partitioned_h2`` and
``run_rg`` take the command option names as keyword arguments and write the
same result directories. Individual estimator methods remain in-memory APIs; batch
partitioned fitting requires an output directory and returns persistent detail
paths. Final HM3 LD tables remain aggregate files.
"""

from __future__ import annotations

from contextlib import contextmanager
import json
import logging
import math
import numbers
import sys
import time
from dataclasses import dataclass, field, replace
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
from typing import Any, Iterable, Sequence

import numpy as np
import pandas as pd

from .._runtime.serial_calculation import serial_calculation
from .._runtime.cli_help import CLIHelpFormatter, SCALAR_PATH_HELP, output_directory_argument
from .._runtime.workflow_logging import args_output_scope, LOG_LEVEL_HELP
from .._core.coordinates import CHR_POS_KEY_COLUMN, build_chr_pos_key_frame
from ..config import (
    GlobalConfig,
    RegressionConfig,
    get_global_config,
    normalize_prevalence_pair,
    print_global_config_banner,
    suppress_global_config_banner,
)
from .._runtime.paths import ensure_output_directory, normalize_path_token, output_directory_path
from .._runtime.workflow_logging import log_inputs, log_outputs, materializing_overwrite_guard, workflow_logging
from . import estimators as reg
from . import partitioned_parallel
from .jackknife import JackknifeIdentifiabilityError
from .._core.snp_identity import (
    clean_identity_artifact_table,
    effective_merge_key_series,
    identity_base_mode,
    identity_mode_family,
    is_allele_aware_mode,
    resolve_regression_identity_mode,
)
from .._core.row_alignment import assert_same_snp_rows, first_identity_mismatch
from .._runtime.parallelism import (
    WORKERS_HELP, parse_query_batch_size, parse_workers, resolve_worker_count, validate_query_batch_size,
    validate_workers,
)
from .._core.columns import infer_chr_pos_columns
from .._core.snp_identity import normalize_snp_identifier_mode
from .._ldscore.results import LDScoreResult
from .._ldscore.source import LDScoreSource, load_ldscore_from_dir
from .outputs import (
    H2DirectoryWriter,
    H2OutputConfig,
    H2_REGRESSION_BIN_COLUMNS,
    LEGACY_SUMSTATS_DROP_COLUMNS,
    PARTITIONED_H2_COLUMNS,
    PartitionedH2DirectoryWriter,
    PartitionedH2OutputConfig,
    PartitionedH2FitArtifacts,
    published_query_artifacts,
    stage_partitioned_h2_fit,
    RG_CONCISE_COLUMNS,
    RG_FULL_COLUMNS,
    RgDirectoryWriter,
    RgOutputConfig,
)
from .._ldscore.outputs import REGRESSION_LD_SCORE_COLUMN
from .._annotation.storage import AnnotationWorkspace
from .._sumstats.table import SumstatsTable, load_sumstats
from ..errors import LDSCInputError, LDSCInternalError, LDSCUsageError, LDSCUserError
from .._runtime.progress import PhaseProgress
from .._genome.build_policy import GenomeBuildCheck
from .overlap_model import assemble_model_overlap, overlap_aware_category_table, model_collinearity_error
from .prevalence import resolve_rg_prevalences
from .._runtime.preflight import InputGate, inspect_artifact_paths, inspect_declared_inputs
from .._runtime.workflow_kwargs import command_parser, run_from_kwargs


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
    "Other causes & fixes: docs/user-guide/troubleshooting.md#partitioned-h2-missing-overlap-matrix"
)
FAILED_RG_NOTE = "Failed; see the rg_full.tsv error column and diagnostics/rg.log."
_NONFINITE_LD_DOC = "docs/user-guide/troubleshooting.md#regression-ld-scores-or-regression-weights-are-missing-or-infinite"
UNESTIMABLE_RG_NOTE = (
    "Unestimable: a jackknife block deletion left the fit unidentifiable; "
    "see the rg_full.tsv error column and diagnostics/rg.log."
)
H2_OUT_OF_BOUNDS_RG_NOTE = (
    "Unestimable: a SNP heritability estimate is not positive (h2 out of bounds); "
    "see the rg_full.tsv error column and diagnostics/rg.log."
)
_RG_H2_OUT_OF_BOUNDS_DOC = "docs/user-guide/troubleshooting.md#rg-a-snp-heritability-estimate-is-not-positive"
REGRESSION_IDENTITY_KEY_COLUMN = "_ldsc_regression_identity_key"
_REGRESSION_NO_OVERLAP_DOC = (
    "docs/user-guide/troubleshooting.md#regression-no-snps-remain-after-merging-regression-inputs"
)
_RG_ALLELES_DOC = "docs/user-guide/troubleshooting.md#rg-a-trait-has-no-a1a2-columns"


@contextmanager
def _log_phase_timing(phase: str, *, level=logging.INFO):
    """Log elapsed wall time for one observable workflow phase."""
    started = time.perf_counter()
    LOGGER.log(level, "Phase start: %s.", phase)
    try:
        yield
    except Exception:
        LOGGER.log(level, "Phase timing: %s failed after %.3fs.", phase, time.perf_counter() - started)
        raise
    LOGGER.log(level, "Phase timing: %s completed in %.3fs.", phase, time.perf_counter() - started)


@dataclass(frozen=True)
class RegressionDataset:
    """Merged regression-ready dataset built from summary statistics and LD-score tables.

    One trait name marks an h2 dataset with ``Z`` and ``N``; two mark an rg
    pair with ``Z1``, ``N1``, ``Z2`` and ``N2``, trait 2 already harmonized to
    trait 1. ``ldscore_overlap`` lets either fit apply the collinearity check.
    An rg pair also counts the SNPs that harmonization dropped:
    ``dropped_orientation_unknown_rows`` are A/T or C/G SNPs, one SNP in both
    traits whose strand, and so the sign of its statistic, cannot be aligned;
    ``dropped_incompatible_allele_rows`` are SNPs whose two traits list
    different allele pairs, usually an input error. Both are zero for h2.
    """
    merged: pd.DataFrame
    ref_ld_columns: list[str]
    weight_column: str
    reference_snp_count_totals: dict[str, np.ndarray]
    count_key_used_for_regression: str
    retained_ld_columns: list[str]
    dropped_zero_variance_ld_columns: list[str]
    trait_names: list[str]
    config_snapshot: GlobalConfig | None = None
    effective_snp_identifier: str = "chr_pos_allele_aware"
    identity_downgrade_applied: bool = False
    ldscore_overlap: "LDScoreOverlap | None" = None
    dropped_orientation_unknown_rows: int = 0
    dropped_incompatible_allele_rows: int = 0

    def validate(self) -> None:
        """Validate that the merged table contains the columns its estimator reads."""
        command, statistics = ("h2", ["Z", "N"]) if len(self.trait_names) == 1 else ("rg", ["Z1", "N1", "Z2", "N2"])
        missing = {"SNP", self.weight_column, *statistics} - set(self.merged.columns)
        if missing:
            raise LDSCInternalError(
                f"{command} regression dataset validation failed: merged table is missing columns {sorted(missing)}. "
                "Most likely regression preprocessing dropped required summary-statistics or LD-score columns. "
                "Re-run with `--log-level DEBUG` and report the traceback."
            )


@dataclass(frozen=True)
class _PreparedRegressionInputs:
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
    ld_variances: dict = field(default_factory=dict)


@dataclass(frozen=True)
class _PartitionedModelInputs:
    """Aligned model state without input readers or the full query source."""

    dataset: RegressionDataset
    counts_by_column: dict[str, dict]
    ld_variances: dict = field(default_factory=dict)


@dataclass(frozen=True)
class RegressionFit:
    """One fitted regression: the estimator and the SNP rows of its final fit.

    Returned by :meth:`RegressionRunner.estimate_h2` and
    :meth:`RegressionRunner.estimate_rg`; the command workflows report from
    the same object, so Python callers see what the result files record.

    Attributes
    ----------
    estimator : Hsq or RG
        The fitted kernel estimator (``Hsq`` for h2, ``RG`` for rg), with
        estimates such as ``tot``, ``tot_se`` and ``intercept``, or
        ``rg_ratio`` and ``rg_se``.
    dataset : RegressionDataset
        The exact rows supplied to the final fit, in SNP order, after the
        chi-square filter, with annotation order, count vectors and SNP
        identifier provenance.
    n_blocks : int
        Jackknife block count actually used.
    effective_chisq_max : float or None
        Chi-square cap actually applied; ``None`` when the fit was uncapped.
    diagnostic_bins : pandas.DataFrame or None
        LD score regression bins of a single-annotation h2 fit, over the same
        rows; ``None`` for rg and partitioned models.
    effective_two_step_cutoff : float or None
        Step-1 chi-square cutoff of the two-step estimator, requested or the
        single-annotation default of 30; ``None`` when the fit used no
        two-step estimator.

    Notes
    -----
    ``n_snps`` counts the rows of ``dataset``. The two-step estimator fits
    its first stage on a subset of these rows; that subset is not the final
    population.
    """

    estimator: reg.Hsq | reg.RG
    dataset: RegressionDataset
    n_blocks: int
    effective_chisq_max: float | None
    diagnostic_bins: pd.DataFrame | None = None
    effective_two_step_cutoff: float | None = None

    @property
    def n_snps(self) -> int:
        """Number of rows actually supplied to the final estimator."""
        return len(self.dataset.merged)


@dataclass(frozen=True)
class H2Result:
    """Result of one ``h2`` run: its summary and the files it wrote.

    Returned by :func:`run_h2`. ``summary`` is the one-row table saved as
    ``h2.tsv``. ``output_paths`` names the written files as
    ``H2DirectoryWriter.write`` returns them: ``summary``, ``metadata``,
    ``ld_score_regression_bins`` for a single-annotation fit, and
    ``legacy_sumstats_drops`` for a legacy input. ``legacy_sumstats_drops``
    holds the rows the legacy ``.sumstats(.gz)`` projection dropped, with the
    columns of ``diagnostics/dropped_snps/legacy_sumstats.tsv.gz``; it is
    ``None`` when the input was not a legacy table. The fitted estimator is
    available from :meth:`RegressionRunner.estimate_h2`.
    """

    summary: pd.DataFrame
    output_paths: dict[str, str] = field(default_factory=dict)
    legacy_sumstats_drops: pd.DataFrame | None = None


@dataclass(frozen=True)
class PartitionedH2BatchResult:
    """Compact batch summary and valid persistent artifacts for every fitted model.

    Annotation tables and coefficient delete values are written after each fit.
    Only the summary and path descriptors remain when the batch returns.
    ``summary`` contains scientific columns only. The saved root TSV appends
    ``query_result_dir`` relative to the output root, or ``NaN`` for a
    baseline-only run; it does not change the original ``category`` labels.
    ``query_status`` contains one row per attempted query in input order, with
    status, failed stage, exception type, and message. Failed queries have no
    scientific summary or model artifacts. Baseline-only runs leave it empty.
    Its columns are ``query_annotation``, ``status``, ``stage``, ``error_type``,
    and ``error_message``. Status is ``success``, ``unestimable`` for singular
    jackknife deletions, or ``failed`` for other per-query exceptions. Successful
    query scans also persist this table as ``diagnostics/query_status.tsv``.
    Worker completion order does not affect status order, stable summary
    sorting, or the final per-query folder order. Worker-count provenance is
    recorded once in root metadata and referenced by successful query metadata.
    ``legacy_sumstats_drops`` holds the rows the legacy ``.sumstats(.gz)``
    projection dropped, with the columns of
    ``diagnostics/dropped_snps/legacy_sumstats.tsv.gz``; it is ``None`` when
    the input was not a legacy table.
    """

    summary: pd.DataFrame
    output_paths: dict[str, str] = field(default_factory=dict)
    per_query_artifacts: dict[str, PartitionedH2FitArtifacts] = field(default_factory=dict)
    query_status: pd.DataFrame = field(default_factory=pd.DataFrame)
    legacy_sumstats_drops: pd.DataFrame | None = None


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
        pairs keep their row and the traits' input prevalences, set the
        estimates to NaN, and store
        ``status='failed'``, or ``status='unestimable'`` when a jackknife block
        deletion left the fit unidentifiable, plus an ``error`` message. A pair
        whose heritability estimate for either trait is not positive is also
        ``unestimable``: it keeps its SNP count and its heritability, genetic
        covariance and intercept estimates, and only ``rg``, ``rg_se``, ``z``
        and ``p`` are NaN.
    h2_per_trait : pandas.DataFrame
        One total-heritability summary per input trait, computed once per trait
        on the trait's single-trait LDSC regression dataset, with ``status``
        (``'ok'``, ``'failed'``, or ``'unestimable'`` when a jackknife block
        deletion left the fit unidentifiable) and an ``error`` message. A
        failed trait keeps its row with NaN estimates, and its pairs are still
        fitted.
    per_pair_metadata : list of dict
        Ordered metadata records aligned to ``rg_full`` rows. Writers use these
        records for optional ``diagnostics/pairs/`` detail output.
    sumstats_files : list of str or None, optional
        Each trait's source path, aligned with ``h2_per_trait`` rows; ``None``
        for an in-memory table. Recorded in the root ``metadata.json``.
    ldscore_dir : str or None, optional
        Source LD-score directory when the LD scores were read from disk.
    legacy_sumstats_drops : pandas.DataFrame or None, optional
        Rows the legacy ``.sumstats(.gz)`` projection dropped, for all legacy
        traits in input order, with the columns of
        ``diagnostics/dropped_snps/legacy_sumstats.tsv.gz``; ``None`` when no
        input was a legacy table.
    dropped_nonfinite_statistic_rows : list of int or None, optional
        Each trait's rows dropped because ``Z`` or ``N`` was missing or
        infinite, aligned with ``h2_per_trait`` rows. Recorded in the root
        ``metadata.json``.
    dropped_orientation_unknown_rows : list of int or None, optional
        Each pair's A/T or C/G SNPs that harmonization dropped because their
        strand, and so the sign of their statistic, cannot be aligned; an
        expected loss of the method. Aligned with ``rg_full`` rows and
        recorded in the root ``metadata.json``.
    dropped_incompatible_allele_rows : list of int or None, optional
        Each pair's SNPs that harmonization dropped because the two traits
        list different allele pairs, usually an input error that only a base
        SNP identifier mode lets through. Aligned with ``rg_full`` rows and
        recorded in the root ``metadata.json``.
    output_paths : dict of str to str, optional
        Files the run wrote, keyed as ``RgDirectoryWriter.write`` returns
        them, including ``legacy_sumstats_drops`` for a legacy input. Empty
        for a family that ``estimate_rg_pairs`` returns without writing.
    """

    rg: pd.DataFrame
    rg_full: pd.DataFrame
    h2_per_trait: pd.DataFrame
    per_pair_metadata: list[dict[str, object]]
    sumstats_files: list[str | None] | None = None
    ldscore_dir: str | None = None
    legacy_sumstats_drops: pd.DataFrame | None = None
    dropped_nonfinite_statistic_rows: list[int] | None = None
    dropped_orientation_unknown_rows: list[int] | None = None
    dropped_incompatible_allele_rows: list[int] | None = None
    output_paths: dict[str, str] = field(default_factory=dict)


@dataclass(frozen=True)
class _PreparedTraitSumstats:
    """Prepared trait summary statistics: one trait, prepared once for every model that uses it.

    For a legacy ``.sumstats(.gz)`` input, ``table`` is its projection onto the
    LD-score panel and ``legacy_drops`` records the rows that projection dropped;
    otherwise ``table`` is the input and ``legacy_drops`` is ``None``. Either way,
    rows whose ``Z`` or ``N`` is missing or infinite are removed and counted in
    ``dropped_nonfinite_statistic_rows``; for a legacy input they are also
    ``legacy_drops`` rows with reason ``nonfinite_statistic``. Rows keyed
    for a merge are built on first use and cached per effective SNP identifier
    mode and downgrade cleanup, so a trait shared by many rg pairs is keyed at
    most twice: once for merges in the LD-score directory's mode and once for
    downgraded merges, which run in its base mode after cleanup.
    """

    table: SumstatsTable
    legacy_drops: pd.DataFrame | None = None
    dropped_nonfinite_statistic_rows: int = 0
    _keyed: dict[tuple[str, bool], tuple[pd.DataFrame, int]] = field(default_factory=dict, repr=False, compare=False)

    @property
    def context(self) -> str:
        """Name the trait in messages."""
        return self.table.source_path or self.table.trait_name or "sumstats"

    def keyed_rows(self, mode: str, *, cleaned: bool) -> tuple[pd.DataFrame, int]:
        """Return rows keyed for a ``mode`` merge and the rows a downgrade cleanup dropped."""
        if (mode, cleaned) not in self._keyed:
            self._keyed[mode, cleaned] = _key_regression_rows(self.table.data, mode, cleaned=cleaned, context=self.context)
        return self._keyed[mode, cleaned]


@dataclass(frozen=True)
class _PreparedLDScores:
    """Baseline LD scores keyed for merges, cached like ``_PreparedTraitSumstats``.

    Keyed frames carry ``_ldsc_source_row``, each row's position in the
    baseline table, so query batches reuse exactly the same join.
    """

    source: LDScoreSource | LDScoreResult
    _keyed: dict[tuple[str, bool], tuple[pd.DataFrame, int]] = field(default_factory=dict, repr=False, compare=False)

    def __post_init__(self):
        """Stop before any fit if a baseline LD score or regression weight is missing or infinite.

        The whole table is checked once per run: every row is a regression
        SNP, so a bad value means the directory itself is damaged, even where
        no trait matches that SNP.
        """
        table = self.source.baseline_table
        problems = _nonfinite_columns(table, [*self.source.baseline_columns, REGRESSION_LD_SCORE_COLUMN], table["SNP"])
        if problems:
            metadata = self.source.output_paths.get("metadata")
            where = f"LD-score directory '{Path(metadata).parent}'" if metadata else "the supplied LD scores"
            raise LDSCInputError(
                f"Regression cannot use {where}: missing or infinite values in {'; '.join(problems)}. "
                "Every LD score and regression weight must be finite; dropping these SNPs would change the "
                "regression without notice. Most likely the directory was edited or only partly written after it "
                "was created. Regenerate it with `ldsc ldscore` or `ldsc convert-legacy-ldscores`, or restore the "
                f"original files; do not replace missing values with zero. Other causes & fixes: {_NONFINITE_LD_DOC}"
            )

    def keyed_rows(self, mode: str, *, cleaned: bool) -> tuple[pd.DataFrame, int]:
        """Return baseline rows keyed for a ``mode`` merge and the rows a downgrade cleanup dropped."""
        if (mode, cleaned) not in self._keyed:
            frame = _assemble_regression_ldscore_table(self.source)
            frame["_ldsc_source_row"] = np.arange(len(frame), dtype=np.int64)
            self._keyed[mode, cleaned] = _key_regression_rows(frame, mode, cleaned=cleaned, context="LD-score table")
        return self._keyed[mode, cleaned]


class RegressionRunner:
    """Assemble regression datasets and fit in workers with 1 BLAS thread each.

    Standalone fits use one worker. Batch partitioned fits reuse their
    requested workers across complete query-plus-baseline models.
    """

    def __init__(
        self,
        global_config: GlobalConfig | None = None,
        regression_config: RegressionConfig | None = None,
    ) -> None:
        """Initialize the runner with shared defaults for regression workflows."""
        self.global_config = global_config or get_global_config()
        self.regression_config = regression_config or RegressionConfig()

    def build_dataset(self, sumstats_table, ldscore_result, config=None, query_columns=None) -> RegressionDataset:
        """Align one trait with the baseline LD scores and select one model's columns and counts.

        If both inputs carry known ``GlobalConfig`` snapshots, their critical
        settings are checked before merging. Unknown-provenance inputs, such as
        disk-loaded summary statistics, are allowed through this compatibility boundary.
        The merge uses the effective key for the resolved mode: ``SNP`` in
        ``rsid``, ``SNP:<allele_set>`` in ``rsid_allele_aware``, ``CHR:POS`` in
        ``chr_pos``, and ``CHR:POS:<allele_set>`` in
        ``chr_pos_allele_aware``. Base modes never inspect alleles for
        identity. ``RegressionConfig.allow_identity_downgrade`` permits
        same-family allele-aware/base mixes to run under the base mode.

        ``query_columns`` adds query LD scores to the baseline model.
        Zero-variance columns and count vectors are selected for that model;
        chi-square filtering, weights, and jackknife work occur inside the fit.
        ``estimate_partitioned_h2_batch`` shares one alignment across query models.
        """
        config = config or self.regression_config
        prepared = self._prepare_h2_inputs(_prepare_trait_sumstats(sumstats_table, ldscore_result),
                                           _PreparedLDScores(ldscore_result), config)
        return self._dataset_from_prepared(prepared, query_columns=query_columns, config=config)

    def _prepare_h2_inputs(
        self,
        trait: _PreparedTraitSumstats,
        ldscores: _PreparedLDScores,
        config: RegressionConfig,
        *,
        build_checked: bool = False,
    ) -> _PreparedRegressionInputs:
        """Merge prepared trait summary statistics with the baseline LD scores.

        This alignment keeps every baseline column, so several models can share
        it; ``_dataset_from_prepared`` selects each model's columns and counts.
        ``build_checked`` skips the genome-build check when the caller already
        checked this trait and LD-score directory together.
        """
        print_global_config_banner(type(self).__name__, self.global_config)
        sumstats_table, ldscore_result = trait.table, ldscores.source
        weight_column = REGRESSION_LD_SCORE_COLUMN
        ref_ld_columns = list(ldscore_result.baseline_columns)
        identity = _resolve_h2_identity(sumstats_table, ldscore_result, self.global_config, config)
        identifier_mode = identity.effective_mode
        if not build_checked:
            _regression_genome_build(sumstats_table, ldscore_result, identifier_mode, context=trait.context)
        ldscore_mode = _ldscore_identity_mode(ldscore_result, self.global_config, fallback_mode=identifier_mode)
        sumstats_keyed, sumstats_dropped = trait.keyed_rows(identifier_mode, cleaned=identity.downgrade_applied)
        ldscore_keyed, ldscore_dropped = ldscores.keyed_rows(identifier_mode, cleaned=identity.downgrade_applied)
        if identity.downgrade_applied:
            sumstats_mode = _sumstats_identity_mode(sumstats_table, fallback_mode=ldscore_mode)
            LOGGER.warning(
                f"Identity downgrade enabled: LD-score mode {ldscore_mode}, summary-statistics mode {sumstats_mode}; "
                f"running regression with effective snp_identifier={identifier_mode!r}. "
                f"Dropped {sumstats_dropped + ldscore_dropped} duplicate effective-key rows before merge."
            )
        key = _merge_key_column(identifier_mode)
        payload = {name: ldscore_keyed[name] for name in [key, *ref_ld_columns, weight_column, "_ldsc_source_row"]}
        # LD-score is the left frame so the inner merge inherits its genomic
        # order; the block jackknife needs genomically contiguous rows for valid
        # SEs, and the sumstats order alone (legacy .sumstats.gz is unsorted) does
        # not guarantee that. Values still align by key regardless of side.
        # The payload borrows the keyed columns; the merge makes the one copy.
        # The fit uses chi^2 = Z^2, so Z keeps the summary statistics' A1 orientation.
        merged = pd.merge(pd.DataFrame(payload, copy=False), sumstats_keyed, how="inner", on=key, sort=False)
        if merged.empty:
            raise LDSCInputError(
                f"h2 regression retained no overlapping {identifier_mode} SNPs after merging summary statistics "
                f"'{trait.context}' with {len(ldscore_keyed)} LD-score rows. Most likely the summary statistics and "
                "LD-score directory use different SNP identifiers or genome builds. Regenerate both inputs "
                "with the same `--snp-identifier` and `--genome-build`, or use `--allow-identity-downgrade` "
                f"for same-family allele-aware/base mixes. Other causes & fixes: {_REGRESSION_NO_OVERLAP_DOC}"
            )

        source_rows = _add_panel_coordinates(merged, ldscore_result)
        dataset = RegressionDataset(
            merged=merged,
            ref_ld_columns=ref_ld_columns,
            weight_column=weight_column,
            reference_snp_count_totals={},
            count_key_used_for_regression="",
            retained_ld_columns=ref_ld_columns,
            dropped_zero_variance_ld_columns=[],
            trait_names=[_trait_label(sumstats_table)],
            config_snapshot=ldscore_result.config_snapshot,
            effective_snp_identifier=identifier_mode,
            identity_downgrade_applied=identity.downgrade_applied,
            ldscore_overlap=ldscore_result.overlap,
        )
        dataset.validate()
        # Baseline variances do not depend on the query, so every model reuses them.
        return _PreparedRegressionInputs(dataset, source_rows, ldscore_result, ldscore_mode,
                                        {str(record["column"]): record for record in ldscore_result.count_records},
                                        dict(_ld_column_variances(merged, ref_ld_columns)))

    def _dataset_from_prepared(self, prepared, *, query_columns=None, config=None) -> RegressionDataset:
        """Select one complete genome-wide model and apply its count/variance rules.

        Query identities are validated when this fit selects query columns.
        Shared preparation never decides model-dependent filters or weights.
        The returned frame borrows aligned columns from ``prepared`` and the
        selected query batch. Treat these inputs as read-only during fitting.
        """
        selected = list(query_columns or [])
        values = self._read_query_batch(prepared, selected) if selected else None
        return _dataset_from_query_values(prepared, selected, values, config or self.regression_config)

    def _read_query_batch(self, prepared, columns):
        """Validate the batch's stored rows once and select the shared trait rows.

        Query rows whose identity columns equal the baseline's exactly need no
        identifier-key comparison; others are compared by key as before.
        """
        frame = prepared.ldscore.read_queries(columns)
        baseline = prepared.ldscore.baseline_table
        identity = [name for name in ("CHR", "SNP", "POS", "A1", "A2") if name in frame.columns]
        if not set(identity) <= set(baseline.columns) or first_identity_mismatch(baseline, frame, identity) is not None:
            assert_same_snp_rows(
                baseline, frame,
                context="query rows must match baseline rows on CHR/SNP/POS",
                snp_identifier=prepared.ldscore_identity_mode,
            )
        rows = prepared.source_rows
        return pd.DataFrame({column: frame[column].to_numpy()[rows] for column in columns}, copy=False)

    def build_rg_dataset(
        self,
        sumstats_table_1: SumstatsTable,
        sumstats_table_2: SumstatsTable,
        ldscore_result: LDScoreSource | LDScoreResult,
        config: RegressionConfig | None = None,
    ) -> RegressionDataset:
        """Build the complete genetic-correlation preprocessing dataset.

        This helper merges trait 1, trait 2, and baseline LD scores on the
        resolved SNP identifier, drops missing rows, aligns trait 2's signed
        statistic to trait 1's ``A1``, and only then removes zero-variance
        LD-score columns on the final RG SNP set. Both traits must carry
        ``A1`` (the effect allele, the allele the signed statistic refers to)
        and ``A2``; otherwise ``LDSCInputError`` is raised.
        """
        _require_rg_alleles([sumstats_table_1, sumstats_table_2])
        return self._rg_dataset(
            _prepare_trait_sumstats(sumstats_table_1, ldscore_result),
            _prepare_trait_sumstats(sumstats_table_2, ldscore_result),
            _PreparedLDScores(ldscore_result),
            config or self.regression_config,
        )

    def _rg_dataset(
        self,
        trait_1: _PreparedTraitSumstats,
        trait_2: _PreparedTraitSumstats,
        ldscores: _PreparedLDScores,
        config: RegressionConfig,
        *,
        build_checked: bool = False,
    ) -> RegressionDataset:
        """Build one pair's dataset from prepared trait summary statistics; see ``build_rg_dataset``.

        ``build_checked`` skips the genome-build check when the caller already
        checked both traits and the LD-score directory together.
        """
        print_global_config_banner(type(self).__name__, self.global_config)
        sumstats_table_1, sumstats_table_2, ldscore_result = trait_1.table, trait_2.table, ldscores.source
        identity = _resolve_rg_identity(sumstats_table_1, sumstats_table_2, ldscore_result, self.global_config, config)
        identifier_mode = identity.effective_mode
        if not build_checked:
            _regression_genome_build([sumstats_table_1, sumstats_table_2], ldscore_result,
                                     identifier_mode, context="rg")

        weight_column = REGRESSION_LD_SCORE_COLUMN
        ref_ld_columns = list(ldscore_result.baseline_columns)
        left_keyed, left_dropped = trait_1.keyed_rows(identifier_mode, cleaned=identity.downgrade_applied)
        right_keyed, right_dropped = trait_2.keyed_rows(identifier_mode, cleaned=identity.downgrade_applied)
        ldscore_keyed, ldscore_dropped = ldscores.keyed_rows(identifier_mode, cleaned=identity.downgrade_applied)
        if identity.downgrade_applied:
            ldscore_mode = _ldscore_identity_mode(ldscore_result, self.global_config, fallback_mode=identifier_mode)
            trait_1_mode = _sumstats_identity_mode(sumstats_table_1, fallback_mode=ldscore_mode)
            trait_2_mode = _sumstats_identity_mode(sumstats_table_2, fallback_mode=ldscore_mode)
            LOGGER.warning(
                f"Identity downgrade enabled: LD-score mode {ldscore_mode}, trait 1 mode {trait_1_mode}, "
                f"trait 2 mode {trait_2_mode}; running rg with effective snp_identifier={identifier_mode!r}. "
                f"Dropped {left_dropped + right_dropped + ldscore_dropped} duplicate effective-key rows before merge."
            )

        key = _merge_key_column(identifier_mode)
        # The merge inputs borrow the keyed columns; each merge makes one copy.
        left = left_keyed.rename(columns={"N": "N1", "Z": "Z1", "FRQ": "FRQ1"}, copy=False)
        right = right_keyed.rename(columns={"A1": "A1x", "A2": "A2x", "N": "N2", "Z": "Z2", "FRQ": "FRQ2"}, copy=False)
        right_payload = [column for column in ["A1x", "A2x", "N2", "Z2", "FRQ2"] if column in right.columns]
        ld_payload = {name: ldscore_keyed[name] for name in [key, *ref_ld_columns, weight_column, "_ldsc_source_row"]}
        # LD-score leads the first merge so the result (and the final merge below)
        # inherits its genomic order for valid jackknife block SEs. See build_dataset.
        left_with_ld = pd.merge(pd.DataFrame(ld_payload, copy=False), left, how="inner", on=key, sort=False)
        merged = pd.merge(left_with_ld, pd.DataFrame({name: right[name] for name in [key, *right_payload]}, copy=False),
                          how="inner", on=key, sort=False)
        del left_with_ld
        _add_panel_coordinates(merged, ldscore_result)
        # Z and N are finite after trait preparation, and LD scores and weights
        # after _PreparedLDScores, so no row needs dropping here.
        if merged.empty:
            raise LDSCInputError(
                f"rg regression retained no overlapping {identifier_mode} SNPs after merging both summary-statistics "
                f"tables with {len(ldscore_keyed)} LD-score rows. Most likely one trait was munged with a "
                "different SNP identifier mode, genome build, or allele convention than the LD-score directory. "
                "Regenerate all inputs with matching identity settings, or use `--allow-identity-downgrade` "
                f"for same-family allele-aware/base mixes. Other causes & fixes: {_REGRESSION_NO_OVERLAP_DOC}"
            )

        compatible, flip = _allele_pair_flags([merged["A1"], merged["A2"], merged["A1x"], merged["A2x"]])
        orientation_unknown = ~compatible & _same_strand_ambiguous_pair(merged)
        drops = (int(orientation_unknown.sum()), int((~compatible & ~orientation_unknown).sum()))
        if drops[0] or drops[1]:
            LOGGER.warning(
                f"rg dropped {drops[0]} A/T or C/G SNPs whose strand cannot be aligned between "
                f"'{_trait_label(sumstats_table_1)}' and '{_trait_label(sumstats_table_2)}', and {drops[1]} SNPs whose "
                "two traits list different allele pairs."
            )
            merged = merged.take(np.flatnonzero(compatible))
            merged.index = pd.RangeIndex(len(merged))
        flip = flip[compatible]
        if merged.empty:
            error = LDSCInputError(
                "rg regression retained no allele-compatible SNPs after harmonizing the two summary-statistics tables "
                f"({drops[0]} A/T or C/G SNPs whose strand cannot be aligned, {drops[1]} SNPs whose two traits list "
                "different allele pairs). Most likely the two traits' A1/A2 columns (A1 is the effect allele, the "
                "allele the signed statistic refers to) hold different allele pairs for the shared SNPs, for example "
                "because one trait was munged from columns that are not its alleles. Re-run `ldsc munge-sumstats` for "
                "both traits with --a1 <effect-allele column> --a2 <other column> and the same SNP identifier mode."
            )
            error.allele_alignment_drops = drops
            raise error
        # The integer signs keep the legacy result dtype (float32 Z widens to float64).
        merged["Z2"] = merged["Z2"].to_numpy() * np.where(flip, -1, 1)
        if "FRQ2" in merged.columns:
            flip_index = pd.Series(flip, index=merged.index)
            valid_frq = pd.to_numeric(merged["FRQ2"], errors="coerce").between(0.0, 1.0, inclusive="both")
            transform = flip_index & valid_frq
            merged.loc[transform, "FRQ2"] = 1.0 - pd.to_numeric(
                merged.loc[transform, "FRQ2"], errors="coerce"
            )

        # Model columns are chosen on the final pair SNP set, after allele harmonization.
        dataset = RegressionDataset(
            ref_ld_columns=ref_ld_columns,
            weight_column=weight_column,
            **_select_model_columns(merged, ref_ld_columns, ldscore_result.count_records,
                                    config.use_common_counts, command="rg"),
            trait_names=[_trait_label(sumstats_table_1), _trait_label(sumstats_table_2)],
            config_snapshot=ldscore_result.config_snapshot,
            effective_snp_identifier=identifier_mode,
            identity_downgrade_applied=identity.downgrade_applied,
            ldscore_overlap=ldscore_result.overlap,
            dropped_orientation_unknown_rows=drops[0],
            dropped_incompatible_allele_rows=drops[1],
        )
        dataset.validate()
        return dataset

    def estimate_h2(
        self,
        dataset: RegressionDataset,
        config: RegressionConfig | None = None,
    ) -> RegressionFit:
        """Estimate single-trait heritability from a prepared dataset.

        Partitioned (multi-annotation) datasets apply the legacy default
        chi-square cap ``max(0.001 * N.max(), 80)`` when ``chisq_max`` is unset,
        down-weighting outlier SNPs. Single-annotation datasets stay uncapped and
        default the two-step cutoff to ``30`` instead. See
        :func:`_resolve_default_chisq_max`.

        Returns
        -------
        RegressionFit
            The ``Hsq`` estimator as ``estimator``, with the SNP rows, block
            count, chi-square cap and diagnostic bins of the fit.
        """
        return _fit_h2_dataset(dataset, config or self.regression_config)

    @serial_calculation
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
            LD-score directory containing baseline LD-score columns and
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
        fit = self.estimate_h2(dataset, config=config)
        samp_prev, pop_prev = _config_prevalence(config or self.regression_config)
        return summarize_partitioned_h2(fit.estimator, fit.dataset, selected_queries, samp_prev=samp_prev,
                                        pop_prev=pop_prev)

    @serial_calculation
    def _fit_functional_model(self, dataset, *, config):
        """Fit and summarize the complete baseline-only model in one child."""
        outcome = _fit_h2_dataset(dataset, config)
        samp_prev, pop_prev = _config_prevalence(config)
        summary = summarize_partitioned_h2(
            outcome.estimator, outcome.dataset, outcome.dataset.retained_ld_columns,
            samp_prev=samp_prev, pop_prev=pop_prev,
        )
        return outcome, summary

    def estimate_partitioned_h2_batch(
        self, sumstats_table, ldscore_result, *, output_dir, query_columns=None,
        config=None, query_batch_size=100, overwrite=False, summary_sort_by="auto", metadata=None,
        continue_on_query_error=False, workers=1,
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
            Maximum query columns loaded at once, default 100; a positive
            integer, otherwise ``LDSCUsageError`` is raised before anything is
            written. Pathways in a
            batch are still fitted separately against the baseline annotations.
            Also caps concurrent query workers; a width of one allows one worker.
            Independent of the batch width used to generate the LD scores.
        workers : int, optional
            Maximum concurrent workers, default 1, including on HPC.
            Use joblib.cpu_count() as the detected available CPU budget.
            Positive requests above that budget fail before capping by task
            count; negative -k uses max(1, budget + 1 - k). Zero, booleans,
            and non-integers are invalid. One worker still runs in a separate process.
            Every worker uses 1 BLAS thread in serial and parallel execution.
            Numerical OpenMP and Arrow CPU/I/O limits are also one; unsupported
            numerical limits are logged as unverified and computation continues.
            This does not detect idle CPUs, guarantee scheduler allocation
            compliance, or ensure sufficient memory.
        overwrite : bool, optional
            Replace workflow-owned outputs, default False.
        summary_sort_by : str, optional
            Final summary ordering; auto uses coefficient-p for focal models.
        metadata : dict, optional
            Additional diagnostic provenance stored once in the root metadata
            and referenced by each model's metadata.
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
            The in-memory summary retains scientific columns only; the saved
            root TSV adds ``query_result_dir`` relative to ``output_dir``, or
            ``NaN`` when no per-query folders are written.
            ``query_status`` records every attempted query in input order;
            scientific summaries and model artifacts include successful fits
            only. Strict runs publish diagnostics alone when any query fails.

        Notes
        -----
        Serial and parallel execution use spawned processes and read-only numeric maps
        inside private output scratch. Each worker owns one complete model;
        batch maps are removed between batches and workers exit before shared
        scratch cleanup. Python scripts at every worker count must launch
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
        validate_workers(workers)
        validate_query_batch_size(query_batch_size, error_type=LDSCUsageError)
        if not isinstance(continue_on_query_error, bool):
            raise ValueError("continue_on_query_error must be a boolean.")
        config = config or self.regression_config
        if ldscore_result.overlap is None:
            raise LDSCInputError(PARTITIONED_H2_REQUIRES_OVERLAP_MESSAGE)
        queries = list(ldscore_result.query_columns if query_columns is None else query_columns)
        queries = _validate_partitioned_query_columns(ldscore_result, queries)
        # Check requested counts before fitting: baseline counts are shared by
        # every model; query counts too unless continuation was requested, which
        # records a query lacking its counts as a failed query instead.
        model_columns = [*ldscore_result.baseline_columns, *([] if continue_on_query_error else queries)]
        count_key = _select_count_key(_count_totals_for_columns(ldscore_result.count_records, model_columns),
                                      config.use_common_counts, model_columns)
        worker_count = resolve_worker_count(workers, max(1, len(queries)),
                                             task_unit="query-plus-baseline model", batch_size=query_batch_size)
        sort_by = _resolve_summary_sort(summary_sort_by, has_queries=bool(queries))
        output_dir = str(output_directory_path(output_dir))
        writer = PartitionedH2DirectoryWriter()
        writer.artifact_family(output_dir, write_per_query_results=bool(queries),
                               coefficient_delete_values=not queries, query_status=bool(queries)).preflight(overwrite=overwrite)
        metadata = dict(metadata or {})
        defaults = {
            "trait_name": _trait_label(sumstats_table),
            "count_kind": "common" if config.use_common_counts else "all",
            "analysis_type": "cell_type_specific" if queries else "functional_category",
            "headline_metric": "coefficient" if queries else "enrichment",
            "enrichment_p_test": "two_sided_t", "coefficient_p_test": "one_sided_greater",
            "annotation_types": dict(ldscore_result.annotation_types),
        }
        ldscore_dir = _ldscore_dir_of(ldscore_result)
        if ldscore_dir:
            defaults["ldscore_dir"] = ldscore_dir
        for name, value in defaults.items():
            metadata.setdefault(name, value)
        if metadata.get("ldscore_dir"):
            metadata["ldscore_dir"] = _absolute_ldscore_dir(metadata["ldscore_dir"])
        with AnnotationWorkspace(output_dir) as workspace:
            with _log_phase_timing("shared regression preparation"):
                trait = _prepare_trait_sumstats(sumstats_table, ldscore_result)
                prepared = self._prepare_h2_inputs(trait, _PreparedLDScores(ldscore_result), config)
                # The trait's provenance, which h2 metadata records as well.
                metadata.update(
                    sumstats_file=sumstats_table.source_path, genome_build=_genome_build_of(prepared.dataset),
                    effective_snp_identifier=prepared.dataset.effective_snp_identifier,
                    identity_downgrade_applied=prepared.dataset.identity_downgrade_applied,
                    count_key_used_for_regression=count_key,
                    dropped_nonfinite_statistic_rows=trait.dropped_nonfinite_statistic_rows,
                    intercept_h2_policy=_intercept_policies(config)["intercept_h2_policy"],
                )
                # Fits need only the aligned rows; release the trait's keyed copy.
                legacy_drops, trait = trait.legacy_drops, None
            if not queries:
                dataset = self._dataset_from_prepared(prepared, config=config)
                if len(dataset.retained_ld_columns) == 1:
                    LOGGER.warning(
                        "partitioned-h2 functional regime retained a single annotation (%s); the fit is "
                        "degenerate and collapses to single-annotation h2 behavior (two-step estimator, no "
                        "default chi-square cap). Most likely the LD-score directory has only one baseline "
                        "annotation. Provide multiple baseline annotations for a meaningful partitioned analysis.",
                        dataset.retained_ld_columns[0],
                    )
                outcome, summary = self._fit_functional_model(dataset, config=config)
                summary = _sort_partitioned_h2_summary(summary, sort_by)
                with _log_phase_timing("output writing"):
                    paths = writer.write(
                        summary, PartitionedH2OutputConfig(output_dir=output_dir, overwrite=overwrite),
                        metadata={**metadata, **_fit_metadata(outcome)},
                        coefficient_delete_values=_coefficient_delete_frame(outcome.estimator, outcome.dataset.retained_ld_columns),
                        legacy_sumstats_drops=legacy_drops,
                    )
                return PartitionedH2BatchResult(summary, paths, legacy_sumstats_drops=legacy_drops)

            rows, private_artifacts, statuses = [], {}, []
            first_error = None
            metadata.update(query_workers_requested=workers, query_workers_effective=worker_count,
                            query_worker_native_threads=1)
            # Workers read the model's numeric columns, SNP and the panel coordinates only.
            fit_columns = {*prepared.dataset.ref_ld_columns, prepared.dataset.weight_column, "Z", "N", "SNP",
                           "_ldsc_panel_CHR", "_ldsc_panel_POS"}
            shared = prepared.dataset.merged
            models = _PartitionedModelInputs(
                replace(prepared.dataset, merged=pd.DataFrame(
                    {name: shared[name] for name in shared.columns if name in fit_columns}, copy=False)),
                prepared.counts_by_column, prepared.ld_variances)
            del shared
            with _log_phase_timing("query worker preparation"):
                worker_pool = partitioned_parallel._QueryWorkers(
                    workspace.path, models, fit_partitioned_query, config, worker_count)
            total_batches = (len(queries) + query_batch_size - 1) // query_batch_size
            del models
            # Loading a batch needs only the LD-score source and row map, not the aligned frame.
            prepared = replace(prepared, dataset=replace(prepared.dataset, merged=pd.DataFrame()))
            with worker_pool as pool, PhaseProgress(LOGGER, "computation", "query fitting", len(queries)) as progress:
                for start in range(0, len(queries), query_batch_size):
                    ordinal = start // query_batch_size + 1
                    with PhaseProgress(LOGGER, "batch", f"query fitting {ordinal}/{total_batches}",
                                       min(query_batch_size, len(queries) - start)) as batch_progress:
                        batch_columns = queries[start:start + query_batch_size]
                        with _log_phase_timing("query batch loading", level=logging.DEBUG):
                            batch = self._read_query_batch(prepared, batch_columns)
                        with TemporaryDirectory(prefix="query-batch-", dir=workspace.path) as batch_root:
                            with _log_phase_timing("query batch mapping", level=logging.DEBUG):
                                mapped = partitioned_parallel._stage_query_batch(batch, Path(batch_root) / "values")
                            batch = None
                            tasks = [(query, mapped[query], workspace.path / f"fit-{start + index}")
                                     for index, query in enumerate(batch_columns)]
                            with _log_phase_timing("query batch execution", level=logging.DEBUG):
                                outcomes = pool.map(tasks, progress=progress)
                        for focal, artifacts, status, error in outcomes:
                            statuses.append(status)
                            if error is not None:
                                first_error = first_error or error
                            else:
                                rows.append(focal)
                                private_artifacts[status["query_annotation"]] = artifacts
                        del outcomes, batch
                        batch_progress.completed = len(batch_columns)
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
                                     per_query_artifacts=private_artifacts, metadata=metadata, query_status=query_status,
                                     legacy_sumstats_drops=legacy_drops)
            return PartitionedH2BatchResult(summary, paths, published_query_artifacts(output_dir, summary), query_status,
                                            legacy_drops)

    @serial_calculation
    def estimate_rg(
        self,
        dataset: RegressionDataset,
        config: RegressionConfig | None = None,
    ) -> RegressionFit:
        """Estimate genetic correlation from a prepared pair dataset.

        Build the dataset with :meth:`build_rg_dataset`, which keeps the final
        trait 1/trait 2/LD-score SNP set explicit and harmonizes trait 2's
        alleles to trait 1. Pass ``chisq_max`` to drop SNPs whose
        ``Z1^2 * Z2^2`` exceeds ``chisq_max^2``; there is no default cap. A
        single-annotation fit with a free heritability intercept defaults the
        two-step cutoff to ``30``, as legacy LDSC does.

        Returns
        -------
        RegressionFit
            The ``RG`` estimator as ``estimator``, with the pair's SNP rows
            after the chi-square filter, block count and the cap applied.

        Raises
        ------
        JackknifeIdentifiabilityError
            If deleting a jackknife block leaves a component fit (trait 1 or
            trait 2 heritability, or the genetic covariance) unidentifiable.
            Each failure names its ``component``; without the two-step
            estimator it also names the unsupported columns and the deleted
            block's genomic span from the LD-score panel.

        Any error raised by the fit carries ``fit_summary``: the pair's SNP and
        block counts, the cap and the dataset settings.
        """
        config = config or self.regression_config
        merged = dataset.merged
        n_snp = len(merged)
        if config.chisq_max is not None:
            # Legacy rg filter (sumstats.py `_rg`): drop SNPs with
            # Z1^2 * Z2^2 > chisq_max^2, kept inclusive (<=) per the -max convention.
            keep = np.ravel((merged["Z1"] ** 2) * (merged["Z2"] ** 2) <= config.chisq_max ** 2)
            if not keep.all():
                merged = merged.take(np.flatnonzero(keep))
                merged.index = pd.RangeIndex(len(merged))
            n_snp = len(merged)
        n_blocks = min(n_snp, config.num_blocks)
        policy = _fit_policy(config)
        intercept_hsq, intercept_gencov = policy.intercept_h2, policy.intercept_gencov
        two_step = policy.two_step_cutoff(len(dataset.retained_ld_columns))
        x = np.asarray(merged[dataset.retained_ld_columns])
        try:
            _raise_on_model_collinearity(dataset, x)
            fitted = reg.RG(
                np.asarray(merged[["Z1"]]),
                np.asarray(merged[["Z2"]]),
                x,
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
        except Exception as error:
            # As for h2, a two-step first stage fits a subset of these rows, so
            # spans and support counts are added only without it.
            if isinstance(error, JackknifeIdentifiabilityError) and two_step is None:
                _add_jackknife_model_context(error, merged, dataset.retained_ld_columns)
            # The pair's rows and settings, for failed-pair reporting.
            error.fit_summary = {"n_snps": n_snp, "n_blocks": n_blocks, "effective_chisq_max": config.chisq_max,
                                 "effective_two_step_cutoff": two_step, **_dataset_metadata(dataset)}
            raise

        return RegressionFit(
            estimator=fitted, dataset=replace(dataset, merged=merged), n_blocks=n_blocks,
            effective_chisq_max=config.chisq_max, effective_two_step_cutoff=two_step,
        )

    @serial_calculation
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
        with ``status='failed'``, or ``status='unestimable'`` for a singular
        jackknife deletion, and do not stop later pairs. A pair whose
        heritability estimate for either trait is not positive is
        ``unestimable`` with NaN ``rg``, ``rg_se``, ``z`` and ``p``, and keeps
        its fitted counts and other estimates. A failure of a trait's
        own heritability fit is recorded the same way in ``h2_per_trait`` and
        does not stop the run.

        Parameters
        ----------
        sumstats_tables : sequence of SumstatsTable
            Munged summary-statistic tables for two or more traits, each with
            ``A1`` (the effect allele) and ``A2``. Inputs sharing a trait label are relabeled
            ``label@<path suffix>`` with the shortest source-path suffix that
            separates them, as in ``ldsc rg``; a shared label without source
            paths raises ``LDSCInputError``.
        ldscore_result : LDScoreSource or LDScoreResult
            LD-score directory. Only baseline LD scores are used for genetic correlation.
        anchor_index : int or None, optional
            If provided, index of the anchor trait. The method computes
            ``anchor`` against every other trait and skips non-anchor pairs.
            Default is ``None``, which computes all unordered pairs.
        config : RegressionConfig or None, optional
            Regression settings. Defaults to the runner's config. Its
            ``samp_prev`` and ``pop_prev`` must be unset: they hold one trait's
            pair, and rg takes the pairs from ``prevalences``.
        prevalences : sequence of (float or None, float or None), optional
            One ``(samp_prev, pop_prev)`` pair per input trait, in input order,
            for liability-scale columns; ``(None, None)`` or NaN marks a
            quantitative trait. Default ``None`` reports every trait on the
            observed scale.

        Returns
        -------
        RgResultFamily
            Concise rg table, full diagnostic table, per-trait SNP heritability table, and
            ordered per-pair metadata.

        Raises
        ------
        LDSCUserError
            Before any fit, if fewer than two traits are supplied,
            ``anchor_index`` is out of bounds, ``prevalences`` does not hold
            one valid pair per trait (``LDSCConfigError`` for an invalid pair),
            or the config sets ``samp_prev`` or ``pop_prev``
            (``LDSCUsageError``). Whole-run input failures, such as
            incompatible h2 datasets, also propagate as ordinary exceptions.
        """
        config = config or self.regression_config
        if config.samp_prev is not None or config.pop_prev is not None:
            raise LDSCUsageError(
                f"rg regression does not read the prevalences of RegressionConfig (samp_prev={config.samp_prev!r}, "
                f"pop_prev={config.pop_prev!r}), which hold one trait's pair; ignoring them would report observed-scale "
                "estimates. Pass one (samp_prev, pop_prev) pair per trait, in input order, with prevalences=..., "
                "and leave samp_prev and pop_prev unset in the config."
            )
        tables = list(sumstats_tables)
        if len(tables) < 2:
            raise LDSCUserError(
                "rg regression requires at least two summary-statistics inputs. Most likely `--sumstats-sources` "
                "resolved to fewer than two files. Pass two or more munged `.sumstats` artifacts."
            )
        if anchor_index is not None and not 0 <= anchor_index < len(tables):
            raise LDSCUsageError(
                f"rg regression received anchor_index={anchor_index}, but valid indices are 0 through {len(tables) - 1}. "
                "Most likely the anchor trait was resolved against a different input list. "
                "Pass an anchor trait name or path that belongs to the current `--sumstats-sources` list."
            )
        if prevalences is not None and len(prevalences) != len(tables):
            raise LDSCUsageError(
                f"rg regression received {len(prevalences)} prevalence pairs for {len(tables)} traits. "
                "Most likely the prevalence list was built for a different input list. "
                "Pass one (samp_prev, pop_prev) pair per trait, in input order."
            )
        prev = ([normalize_prevalence_pair(*pair) for pair in prevalences] if prevalences is not None
                else [(None, None)] * len(tables))
        # Reject allele-free traits before any fit; a pair-level error would be
        # recorded as a failed pair instead of stopping the run.
        _require_rg_alleles(tables)
        tables = _disambiguate_trait_names(tables)
        # Run-level preparation: each trait is projected once, and one joint
        # genome-build check covers every pair, whose effective SNP identifier
        # mode stays in the LD-score directory's family.
        traits = [_prepare_trait_sumstats(table, ldscore_result) for table in tables]
        tables = [trait.table for trait in traits]
        _regression_genome_build(tables, ldscore_result,
                                 _ldscore_identity_mode(ldscore_result, self.global_config), context="rg")
        ldscores = _PreparedLDScores(ldscore_result)

        h2_rows = []
        with suppress_global_config_banner(), PhaseProgress(LOGGER, "computation", "per-trait h2", len(tables)) as progress:
            for trait, (samp_prev, pop_prev) in zip(traits, prev):
                trait_name = _trait_label(trait.table)
                h2_dataset = None
                try:
                    # As for pairs, an error fails only this trait's own fit; its pairs are still fitted.
                    h2_dataset = self._dataset_from_prepared(
                        self._prepare_h2_inputs(trait, ldscores, config, build_checked=True), config=config,
                    )
                    outcome = _fit_h2_dataset(h2_dataset, config, log_level=logging.DEBUG)
                    row = summarize_total_h2(
                        outcome.estimator, outcome.dataset, trait_name=trait_name, n_snps_used=outcome.n_snps,
                        samp_prev=samp_prev, pop_prev=pop_prev,
                    ).assign(status="ok", error="")
                except Exception as exc:
                    status = "unestimable" if isinstance(exc, JackknifeIdentifiabilityError) else "failed"
                    LOGGER.warning(f"SNP heritability of trait '{trait_name}' marked {status}: {_format_exception(exc)}",
                                   exc_info=True)
                    if status == "unestimable":
                        LOGGER.error("Trait %r jackknife diagnostics: %s", trait_name, json.dumps(exc.failures))
                    row = _failed_h2_row(trait_name=trait_name, error=_format_exception(exc), status=status,
                                         n_snps=0 if h2_dataset is None else len(h2_dataset.merged),
                                         samp_prev=samp_prev, pop_prev=pop_prev)
                h2_rows.append(row)
                progress.advance(outcome=row["status"].iloc[0])
                h2_dataset = outcome = None  # release this trait before the next is built
        h2_per_trait = pd.concat(h2_rows, axis=0, ignore_index=True) if h2_rows else pd.DataFrame()

        rg_full_rows: list[dict[str, object]] = []
        rg_rows: list[dict[str, object]] = []
        per_pair_metadata: list[dict[str, object]] = []
        alignment_drops: list[tuple[int, int]] = []
        pair_kind = "anchor" if anchor_index is not None else "all_pairs"
        pair_count = len(tables) - 1 if anchor_index is not None else len(tables) * (len(tables) - 1) // 2
        with suppress_global_config_banner(), PhaseProgress(LOGGER, "computation", "rg pairs", pair_count) as progress:
            for i, j in _iter_rg_pairs(len(tables), anchor_index):
                trait_1 = _trait_label(tables[i])
                trait_2 = _trait_label(tables[j])
                drops = (0, 0)
                try:
                    # Pair-specific keyed rows are built here, so an error fails only this pair.
                    dataset = self._rg_dataset(traits[i], traits[j], ldscores, config, build_checked=True)
                    drops = (dataset.dropped_orientation_unknown_rows, dataset.dropped_incompatible_allele_rows)
                    outcome = self.estimate_rg(dataset, config=config)
                    fitted, dataset = outcome.estimator, outcome.dataset
                    n_snps_used = outcome.n_snps
                    full_row = _summarize_rg_pair(
                        fitted, dataset, trait_1=trait_1, trait_2=trait_2, pair_kind=pair_kind,
                        n_snps_used=n_snps_used,
                        samp_prev_1=prev[i][0], pop_prev_1=prev[i][1],
                        samp_prev_2=prev[j][0], pop_prev_2=prev[j][1],
                    )
                    note = ""
                    if full_row["status"] == "unestimable":
                        note = H2_OUT_OF_BOUNDS_RG_NOTE
                        LOGGER.warning(f"rg for pair '{trait_1}' vs '{trait_2}' marked unestimable: {full_row['error']}")
                    metadata = _rg_pair_metadata(tables[i], tables[j], full_row, config, pair_kind, (prev[i], prev[j]),
                                                 {**_fit_metadata(outcome, "n_snps_used"), **_dataset_metadata(dataset)})
                except Exception as exc:
                    drops = getattr(exc, "allele_alignment_drops", drops)
                    error = _format_exception(exc)
                    status = "unestimable" if isinstance(exc, JackknifeIdentifiabilityError) else "failed"
                    note = UNESTIMABLE_RG_NOTE if status == "unestimable" else FAILED_RG_NOTE
                    LOGGER.warning(f"rg for pair '{trait_1}' vs '{trait_2}' marked {status}: {error}", exc_info=True)
                    if status == "unestimable":
                        LOGGER.error("Pair %r vs %r jackknife diagnostics: %s", trait_1, trait_2, json.dumps(exc.failures))
                    fit_summary = getattr(exc, "fit_summary", None)
                    full_row = _failed_rg_full_row(trait_1=trait_1, trait_2=trait_2, pair_kind=pair_kind, error=error,
                                                   status=status, n_snps_used=fit_summary["n_snps"] if fit_summary else 0,
                                                   samp_prev_1=prev[i][0], pop_prev_1=prev[i][1],
                                                   samp_prev_2=prev[j][0], pop_prev_2=prev[j][1])
                    fit = None if fit_summary is None else {
                        "n_snps_used": fit_summary["n_snps"], "n_blocks_used": fit_summary["n_blocks"],
                        **{key: value for key, value in fit_summary.items() if key not in {"n_snps", "n_blocks"}},
                    }
                    metadata = _rg_pair_metadata(tables[i], tables[j], full_row, config, pair_kind, (prev[i], prev[j]),
                                                 fit)
                rg_full_rows.append(full_row)
                rg_rows.append(_concise_rg_row(full_row, note))
                per_pair_metadata.append(metadata)
                alignment_drops.append(drops)
                progress.advance(outcome=full_row["status"])
                dataset = outcome = fitted = None  # release this pair before the next is built

        rg_full = pd.DataFrame(rg_full_rows, columns=RG_FULL_COLUMNS)
        rg = pd.DataFrame(rg_rows, columns=RG_CONCISE_COLUMNS)
        legacy_drops = [trait.legacy_drops for trait in traits if trait.legacy_drops is not None]
        return RgResultFamily(
            rg=rg, rg_full=rg_full, h2_per_trait=h2_per_trait, per_pair_metadata=per_pair_metadata,
            sumstats_files=[table.source_path for table in tables],
            ldscore_dir=_ldscore_dir_of(ldscore_result),
            legacy_sumstats_drops=pd.concat(legacy_drops, ignore_index=True) if legacy_drops else None,
            dropped_nonfinite_statistic_rows=[trait.dropped_nonfinite_statistic_rows for trait in traits],
            dropped_orientation_unknown_rows=[unknown for unknown, _ in alignment_drops],
            dropped_incompatible_allele_rows=[incompatible for _, incompatible in alignment_drops],
        )


def _dataset_from_query_values(prepared, selected, query_values, config):
    """Borrow aligned columns for one model without copying the baseline.

    Preparation and query loading already provide matching RangeIndexes.
    Borrowed columns must remain unmodified; numerical fitting allocates its
    own design array after the model's chi-square filter has selected rows.
    """
    ref_ld_columns = list(prepared.dataset.ref_ld_columns) + list(selected)
    merged = prepared.dataset.merged
    # A query's damage fails only its own model, inside the query's error handling.
    problems = _nonfinite_columns(query_values, selected, merged["SNP"]) if selected else []
    if problems:
        raise LDSCInputError(
            f"Regression cannot use query annotation {'; '.join(problems)}: its LD scores are missing or "
            "infinite on regression SNPs, and every LD score must be finite. Most likely the LD-score directory "
            "was edited or only partly written after it was created. Regenerate it with `ldsc ldscore`, or "
            f"restore the original query files. Other causes & fixes: {_NONFINITE_LD_DOC}"
        )
    if selected:
        weight_position = merged.columns.get_loc(prepared.dataset.weight_column)
        order = [*merged.columns[:weight_position], *selected, *merged.columns[weight_position:]]
        merged = pd.DataFrame({name: (query_values if name in selected else merged)[name].to_numpy(copy=False)
                               for name in order}, copy=False)
    dataset = replace(prepared.dataset, ref_ld_columns=ref_ld_columns, **_select_model_columns(
        merged, ref_ld_columns, prepared.counts_by_column, config.use_common_counts, command="h2",
        known_variances=prepared.ld_variances,
    ))
    dataset.validate()
    return dataset


@serial_calculation
def _fit_h2_dataset(
    dataset: RegressionDataset,
    config: RegressionConfig,
    *, log_level=logging.INFO, retain_rows=True,
) -> RegressionFit:
    """Fit SNP heritability once and retain its exact SNP population and diagnostic bins.

    The aligned frame can borrow read-only worker maps. A multi-annotation
    design is filled column by column from the rows the chi-square filter
    keeps, so the frame itself is not copied before the fit; the filtered
    frame is built afterwards, for the returned dataset. With
    ``retain_rows=False`` (per-query partitioned fits, which read only the
    row count) the returned dataset's ``merged`` is an empty frame with one
    row per fitted SNP. The estimator never modifies the borrowed input columns.
    """
    merged = dataset.merged
    n_snp = len(merged)
    n_blocks = min(n_snp, config.num_blocks)
    chisq = np.asarray(merged["Z"] ** 2).reshape((n_snp, 1))
    # N.max() is read pre-filter so the partitioned default cap matches legacy.
    chisq_max = _resolve_default_chisq_max(
        config.chisq_max, len(dataset.retained_ld_columns), merged["N"].max()
    )
    keep = None
    if chisq_max is not None:
        mask = np.ravel(chisq <= chisq_max)
        n_removed = n_snp - int(mask.sum())
        if n_removed:
            keep = mask
            n_snp = int(mask.sum())
            n_blocks = min(n_snp, config.num_blocks)
            chisq = chisq[keep]
            LOGGER.log(log_level,
                "Removed %d SNPs with chi^2 > %g (%d SNPs remain).",
                n_removed,
                chisq_max,
                n_snp,
            )
    columns = dataset.retained_ld_columns
    if keep is not None and len(columns) == 1:
        # Single-annotation fits read the filtered frame for their diagnostic bins.
        merged, keep = _kept_rows(merged, keep), None
    x = _ld_design(merged, columns, keep)
    policy = _fit_policy(config)
    intercept, two_step = policy.intercept_h2, policy.two_step_cutoff(len(columns))
    old_weights = len(columns) > 1
    _raise_on_model_collinearity(dataset, x)
    reference_snp_counts = np.asarray(
        dataset.reference_snp_count_totals[dataset.count_key_used_for_regression]
    ).reshape((1, -1))
    try:
        hsq = reg.Hsq(
            chisq,
            x,
            np.asarray(_kept_rows(merged[[dataset.weight_column]], keep)),
            np.asarray(_kept_rows(merged[["N"]], keep)),
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
            _add_jackknife_model_context(error, _kept_rows(merged, keep), columns)
        raise
    diagnostic_bins = None
    if len(columns) == 1:
        diagnostic_bins = summarize_ld_score_regression_bins(
            hsq,
            ld_score=x[:, 0],
            chi_square=chisq,
            sample_size=np.asarray(merged[["N"]]),
            regression_ld_score=np.asarray(merged[[dataset.weight_column]]),
            reference_snp_count=float(reference_snp_counts[0, 0]),
        )
    del x
    merged = _kept_rows(merged, keep) if retain_rows else pd.DataFrame(index=pd.RangeIndex(n_snp))
    return RegressionFit(
        estimator=hsq, dataset=replace(dataset, merged=merged), n_blocks=n_blocks,
        effective_chisq_max=chisq_max, diagnostic_bins=diagnostic_bins, effective_two_step_cutoff=two_step,
    )


def fit_partitioned_query(prepared, query, values, directory, config):
    """Fit, summarize and privately stage one complete genome-wide query.

    Every query worker calls this module-level function, which spawn pickles
    by name; only model errors become query outcomes, while staging errors
    and interrupts still abort runs. Returned objects contain no estimator,
    SNP matrix or traceback frames.
    """
    stage = "model_preparation"
    try:
        with _log_phase_timing("query model preparation", level=logging.DEBUG):
            dataset = _dataset_from_query_values(prepared, [query], values, config)
        stage = "estimator"
        with _log_phase_timing("estimator execution", level=logging.DEBUG):
            outcome = _fit_h2_dataset(dataset, config, log_level=logging.DEBUG, retain_rows=False)
        hsq, dataset = outcome.estimator, outcome.dataset
        stage = "summary"
        with _log_phase_timing("query summary", level=logging.DEBUG):
            samp_prev, pop_prev = _config_prevalence(config)
            full = summarize_partitioned_h2(hsq, dataset, dataset.retained_ld_columns,
                                           samp_prev=samp_prev, pop_prev=pop_prev)
            # The summary computes every retained row; the focal row is one of them.
            focal = (full.set_index("category").loc[[query]].reset_index()
                     .loc[:, PARTITIONED_H2_COLUMNS].reset_index(drop=True))
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
    with _log_phase_timing("query output staging", level=logging.DEBUG):
        artifacts = stage_partitioned_h2_fit(
            directory, full, delete_values,
            {"dropped_zero_variance_ld_columns": list(dataset.dropped_zero_variance_ld_columns),
             **_fit_metadata(outcome),
             "effective_snp_identifier": dataset.effective_snp_identifier,
             "identity_downgrade_applied": dataset.identity_downgrade_applied},
        )
    return focal, artifacts, {"query_annotation": query, "status": "success", "stage": "complete",
                              "error_type": "", "error_message": ""}, None


def _regression_genome_build(
    sumstats_tables: SumstatsTable | Sequence[SumstatsTable],
    ldscore_result: LDScoreSource | LDScoreResult,
    identifier_mode: str,
    *,
    context: str,
) -> None:
    """Check all coordinate-family artifacts before fitting any model.

    Valid snapshots establish builds without coordinate re-inference.
    In-memory objects lacking a concrete snapshot use bounded coordinate
    samples. rsID joins need no coordinate build. Regression has no build
    override; failures direct users to matching upstream artifacts.
    """
    if identity_mode_family(identifier_mode) != "chr_pos":
        return

    check = GenomeBuildCheck('auto', context='regression', option=None, logger=LOGGER)
    tables = sumstats_tables if isinstance(sumstats_tables, (list, tuple)) else [sumstats_tables]
    role = '--sumstats-sources' if isinstance(sumstats_tables, (list, tuple)) else '--sumstats-file'
    for table in tables:
        source = table.source_path or table.trait_name or context
        build = getattr(table.config_snapshot, 'genome_build', None)
        if _has_legacy_sumstats_source(table):
            check.unchecked(role, [source], 'legacy rsID projection inherits panel coordinates')
        elif build in {'hg19', 'hg38'}:
            check.declared(role, source, build)
        else:
            check.frames(role, source, [table.data])
    panel_source = str(ldscore_result.output_paths.get('metadata', 'LD-score directory'))
    build = getattr(ldscore_result.config_snapshot, 'genome_build', None)
    if build in {'hg19', 'hg38'}:
        check.declared('--ldscore-dir', panel_source, build)
    else:
        check.frames('--ldscore-dir', panel_source, [ldscore_result.baseline_table])
    check.resolve()


def _snapshot_identity_mode(snapshot: GlobalConfig | None) -> str | None:
    """Return a normalized SNP identifier mode from ``snapshot`` if present."""
    value = getattr(snapshot, "snp_identifier", None)
    if value is None:
        return None
    return normalize_snp_identifier_mode(value)


def _has_legacy_sumstats_source(table: SumstatsTable) -> bool:
    """Return whether ``table`` originated from legacy text."""
    return bool(table.provenance.get("legacy_ldsc", False))


def _is_allele_less_legacy(table: SumstatsTable) -> bool:
    """Return whether ``table`` is legacy text without A1/A2, which projects by rsID under the panel's base mode."""
    return _has_legacy_sumstats_source(table) and not {"A1", "A2"}.intersection(table.data.columns)


def _prepare_trait_sumstats(
    table: SumstatsTable,
    ldscore_result: LDScoreSource | LDScoreResult,
) -> _PreparedTraitSumstats:
    """Prepare one trait's summary statistics once for every model that uses it.

    A legacy ``.sumstats(.gz)`` table is projected onto the LD-score panel under
    the panel's recorded SNP identifier mode; other tables are used as given.
    Rows whose ``Z`` or ``N`` is missing or infinite cannot enter a fit and are
    dropped with a warning instead of failing the run (a recorded exception to
    the fail-on-unrequested-filtering rule, as in legacy LDSC).
    """
    context = table.source_path or table.trait_name or "sumstats"
    legacy_drops = None
    if _has_legacy_sumstats_source(table):
        table, legacy_drops = _project_legacy_sumstats_to_panel(table, ldscore_result)
        dropped = int(legacy_drops["reason"].eq("nonfinite_statistic").sum())
    else:
        finite = _finite_statistics(table.data)
        dropped = int((~finite).sum())
        if dropped:
            kept = table.data.take(np.flatnonzero(finite.to_numpy()))
            kept.index = pd.RangeIndex(len(kept))
            table = replace(table, data=kept)
    if dropped:
        LOGGER.warning(
            "Dropped %d row(s) of summary statistics '%s' whose Z or N is missing or infinite; such rows cannot "
            "enter a regression. The count is recorded as dropped_nonfinite_statistic_rows in the result metadata.",
            dropped, context,
        )
    return _PreparedTraitSumstats(table, legacy_drops, dropped)


def _finite_statistics(frame: pd.DataFrame) -> pd.Series:
    """Mark rows whose ``Z`` and ``N`` are both finite numbers."""
    return (np.isfinite(pd.to_numeric(frame["Z"], errors="coerce"))
            & np.isfinite(pd.to_numeric(frame["N"], errors="coerce")))


def _project_legacy_sumstats_to_panel(
    table: SumstatsTable,
    ldscore_result: LDScoreSource | LDScoreResult,
) -> tuple[SumstatsTable, pd.DataFrame]:
    """Project one legacy munged summary-statistics table onto the regression SNPs in the LD-score directory.

    Legacy ``SNP`` values are treated only as rsID lookup keys. The supplied
    LD-score panel owns output coordinates and identifier mode. Allele-aware
    panels also own allele order, except for A/T and C/G SNPs, whose strand
    against the panel is unknown: they keep their source alleles and ``Z``.
    Allele-unaware panels retain the source alleles so two traits can still
    be harmonized for genetic correlation.

    Parameters
    ----------
    table : SumstatsTable
        Legacy text table marked with ``provenance['legacy_ldsc']``.
    ldscore_result : LDScoreSource or LDScoreResult
        LDSC3 LD-score directory used by the pending regression.

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
        If required legacy columns or the LD-score directory's SNP identifier metadata are missing,
        or no source row can be projected safely.
    """
    label = table.source_path or table.trait_name or 'sumstats'
    missing = sorted({"SNP", "Z", "N"} - set(table.data.columns))
    if missing:
        raise LDSCInputError(
            f"Cannot use legacy summary-statistics file '{label}': required columns are missing: {missing}. "
            "Legacy regression compatibility requires SNP, Z and N, and A1 and A2 as a pair when present."
        )
    allele_columns = [column for column in ("A1", "A2") if column in table.data.columns]
    if len(allele_columns) == 1:
        present = allele_columns[0]
        absent = "A2" if present == "A1" else "A1"
        raise LDSCInputError(
            f"Cannot use legacy summary-statistics file '{label}': it has {present} but not {absent}. "
            "A1 (the effect allele, the allele the signed statistic refers to) and A2 (the other allele) must come "
            "as a pair. Re-create the file with both allele columns, or without either for h2 and partitioned-h2."
        )
    has_alleles = len(allele_columns) == 2
    panel_snapshot = ldscore_result.config_snapshot
    panel_mode = _snapshot_identity_mode(panel_snapshot)
    if panel_snapshot is None or panel_mode is None:
        raise LDSCInputError(
            "Legacy summary-statistics projection requires an LD-score directory with recorded "
            "snp_identifier metadata. Regenerate or explicitly convert the LD-score directory first."
        )
    panel_source = ldscore_result.baseline_table
    missing_panel = sorted({"CHR", "POS", "SNP"} - set(panel_source.columns))
    if missing_panel:
        raise LDSCInputError(
            f"Legacy summary-statistics projection cannot use the LD-score panel because it is missing {missing_panel}. "
            "Regenerate the LD-score directory."
        )
    # Without alleles, rows project by rsID alone, under the panel's base mode.
    allele_aware = is_allele_aware_mode(panel_mode) and has_alleles
    projected_snapshot = panel_snapshot if has_alleles else replace(
        panel_snapshot, snp_identifier=identity_base_mode(panel_mode))
    if allele_aware and not {"A1", "A2"}.issubset(panel_source.columns):
        raise LDSCInputError(
            f"Legacy summary-statistics projection requires panel A1/A2 for snp_identifier='{panel_mode}', "
            "but the canonical LD-score table does not contain both allele columns."
        )

    source = table.data.copy(deep=False)  # normalized columns are replaced, not modified in place
    source.index = pd.RangeIndex(len(source))
    source["SNP"] = source["SNP"].astype("string").str.strip()
    if has_alleles:
        source["A1"] = source["A1"].astype("string").str.upper().str.strip()
        source["A2"] = source["A2"].astype("string").str.upper().str.strip()
    source_rsids = set(source["SNP"].dropna().astype(str)) - {""}
    panel_snp = panel_source["SNP"].astype("string").str.strip()
    panel_columns = ["CHR", "POS", "SNP", *(["A1", "A2"] if allele_aware else [])]
    in_source = np.flatnonzero(panel_snp.isin(source_rsids).to_numpy())
    panel = panel_source.loc[:, panel_columns].take(in_source)
    panel.index = pd.RangeIndex(len(panel))
    panel["SNP"] = panel_snp.take(in_source).to_numpy()
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
                "Legacy summary-statistics file '%s' has %d nonnumeric, infinite, or out-of-range FRQ value(s); "
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
    if has_alleles:
        assign_reason(
            source["A1"].isna() | source["A2"].isna() | source["A1"].eq("") | source["A2"].eq(""),
            "missing_allele",
        )
        assign_reason(
            ~source["A1"].isin(valid_bases) | ~source["A2"].isin(valid_bases) | source["A1"].eq(source["A2"]),
            "invalid_allele",
        )
    assign_reason(source_candidate_count.eq(0), "missing_panel_rsid")

    eligible = reasons.isna()
    selected_panel: pd.DataFrame
    flip_by_source = pd.Series(False, index=source.index, dtype=bool)
    if allele_aware:
        source_candidates = source.loc[eligible, ["SNP", "A1", "A2"]]
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
        complement_a1 = pairs["A1"].map(reg.COMPLEMENT)
        complement_a2 = pairs["A2"].map(reg.COMPLEMENT)
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
        compatible_pairs = pairs.take(np.flatnonzero(compatible.to_numpy()))
        # An A/T or C/G pair is its own strand complement, so it matches both
        # ways round and its orientation is unknown: it is never flipped.
        flip = (swapped | swapped_complement) & complement_a1.ne(pairs["A2"])
        compatible_pairs["_legacy_flip"] = flip.loc[compatible_pairs.index].to_numpy()
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
        source_matches = source.loc[accepted_rows, ["SNP"]]
        source_matches["_legacy_source_row"] = source_matches.index
        unique_panel = panel.loc[panel["SNP"].isin(source_matches["SNP"])]
        selected_panel = source_matches.merge(unique_panel, how="left", on="SNP", sort=False, validate="one_to_one")
        selected_panel = selected_panel.set_index("_legacy_source_row").loc[accepted_rows]
    # Checked last, so a row with an earlier reason keeps it.
    assign_reason(~_finite_statistics(source), "nonfinite_statistic")

    accepted_rows = reasons[reasons.isna()].index
    if accepted_rows.empty:
        counts = reasons.dropna().value_counts().to_dict()
        raise LDSCInputError(
            f"Legacy summary-statistics file '{table.source_path or table.trait_name or 'sumstats'}' retained no SNPs "
            f"after projection onto the regression SNPs in the LD-score directory (drop counts: {counts})."
        )
    projected_frame = source.take(accepted_rows.to_numpy())
    panel_column = (lambda column: f"_legacy_panel_{column}") if allele_aware else (lambda column: column)
    projected_frame["CHR"] = selected_panel.loc[accepted_rows, panel_column("CHR")].to_numpy()
    projected_frame["POS"] = selected_panel.loc[accepted_rows, panel_column("POS")].to_numpy()
    projected_frame["SNP"] = selected_panel.loc[accepted_rows, panel_column("SNP")].to_numpy()
    if allele_aware:
        # A/T and C/G rows keep their source alleles, the orientation their Z refers to.
        oriented = projected_frame["A1"].map(reg.COMPLEMENT).ne(projected_frame["A2"]).to_numpy()
        for column in ("A1", "A2"):
            panel_alleles = selected_panel.loc[accepted_rows, panel_column(column)].to_numpy()
            projected_frame[column] = np.where(oriented, panel_alleles, projected_frame[column].to_numpy())
        flipped_rows = flip_by_source.loc[accepted_rows]
        projected_frame.loc[flipped_rows, "Z"] = -pd.to_numeric(
            projected_frame.loc[flipped_rows, "Z"], errors="raise"
        )
        if "FRQ" in projected_frame.columns:
            valid_flipped_frq = flipped_rows & projected_frame["FRQ"].notna()
            projected_frame.loc[valid_flipped_frq, "FRQ"] = 1.0 - projected_frame.loc[valid_flipped_frq, "FRQ"]
    projected_frame.index = pd.RangeIndex(len(projected_frame))
    projected_frame["POS"] = pd.to_numeric(projected_frame["POS"], errors="raise").astype("int64")
    dropped_mask = reasons.notna()
    drops = source.reindex(columns=["SNP", "A1", "A2"]).loc[dropped_mask]
    drops.insert(0, "source_path", table.source_path)
    drops.insert(0, "trait_name", table.trait_name or Path(table.source_path or "sumstats").name)
    drops["reason"] = reasons.loc[dropped_mask].to_numpy()
    drops["panel_candidate_count"] = source_candidate_count.loc[dropped_mask].to_numpy()
    drops = drops.loc[:, LEGACY_SUMSTATS_DROP_COLUMNS]
    drops.index = pd.RangeIndex(len(drops))
    drop_counts = drops["reason"].value_counts(sort=False)
    if not drop_counts.empty:
        LOGGER.warning(
            "Legacy summary-statistics file '%s' projection dropped %d row(s): %s.",
            table.source_path or table.trait_name or "sumstats",
            len(drops),
            ", ".join(f"{reason}={int(count)}" for reason, count in drop_counts.items()),
        )
    LOGGER.info(
        "Projected %d legacy summary-statistics row(s) from '%s' onto panel SNP identifier mode %s; dropped %d.",
        len(projected_frame),
        table.source_path or table.trait_name or "sumstats",
        projected_snapshot.snp_identifier,
        len(drops),
    )
    return replace(table, data=projected_frame, has_alleles=has_alleles, config_snapshot=projected_snapshot), drops


def _sumstats_identity_mode(
    sumstats_table: SumstatsTable,
    *,
    fallback_mode: str,
) -> str:
    """Return the mode represented by one summary-statistics input for regression."""
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
    """Resolve one summary-statistics table and one LD-score result for h2 merging."""
    ldscore_mode = _ldscore_identity_mode(ldscore_result, runner_config)
    sumstats_mode = _sumstats_identity_mode(sumstats_table, fallback_mode=ldscore_mode)
    # Allele-less legacy text has no allele-aware alternative, so its projection is downgraded without the flag.
    return resolve_regression_identity_mode(
        sumstats_mode,
        ldscore_mode,
        allow_identity_downgrade=regression_config.allow_identity_downgrade or _is_allele_less_legacy(sumstats_table),
    )


def _resolve_rg_identity(
    sumstats_table_1: SumstatsTable,
    sumstats_table_2: SumstatsTable,
    ldscore_result: LDScoreSource | LDScoreResult,
    runner_config: GlobalConfig,
    regression_config: RegressionConfig,
):
    """Resolve trait 1, trait 2, and LD-score SNP identifier modes for rg."""
    ldscore_mode = _ldscore_identity_mode(ldscore_result, runner_config)
    trait_1_mode = _sumstats_identity_mode(sumstats_table_1, fallback_mode=ldscore_mode)
    trait_2_mode = _sumstats_identity_mode(sumstats_table_2, fallback_mode=ldscore_mode)
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


def _key_regression_rows(frame: pd.DataFrame, mode: str, *, cleaned: bool, context: str) -> tuple[pd.DataFrame, int]:
    """Key one side of a regression merge for ``mode``.

    ``cleaned`` first applies the identity-downgrade cleanup. Allele-aware modes
    add the effective identity key and coordinate-family base modes the
    ``CHR:POS`` key; ``rsid`` merges on ``SNP`` as given. Returns the keyed rows
    and the number of duplicate effective-key rows the cleanup dropped.
    """
    dropped = 0
    if cleaned:
        frame, dropped_rows = _prepare_regression_identity_table(frame, mode, context=context, logger=LOGGER)
        dropped = len(dropped_rows)
    if is_allele_aware_mode(mode):
        frame = _with_effective_identity_key(frame, mode, context=context)
    elif identity_mode_family(mode) != "rsid":
        frame = _with_chr_pos_key(frame, context=context)
    return frame, dropped


def _merge_key_column(mode: str) -> str:
    """Return the column ``_key_regression_rows`` merges on under ``mode``."""
    if is_allele_aware_mode(mode):
        return REGRESSION_IDENTITY_KEY_COLUMN
    return "SNP" if identity_mode_family(mode) == "rsid" else CHR_POS_KEY_COLUMN


def _with_effective_identity_key(frame: pd.DataFrame, mode: str, *, context: str) -> pd.DataFrame:
    """Return the rows with an effective regression identity key, with that key added.

    The result shares the input's columns unless rows are dropped; callers
    only read it (merges copy).
    """
    keyed = frame.copy(deep=False)
    keyed[REGRESSION_IDENTITY_KEY_COLUMN] = effective_merge_key_series(keyed, mode, context=context)
    present = keyed[REGRESSION_IDENTITY_KEY_COLUMN].notna().to_numpy()
    if not present.all():
        keyed = keyed.take(np.flatnonzero(present))
    keyed.index = pd.RangeIndex(len(keyed))
    return keyed


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


def _kept_rows(frame: pd.DataFrame, keep: np.ndarray | None) -> pd.DataFrame:
    """Return the rows ``keep`` marks with a fresh RangeIndex, or ``frame`` itself when ``keep`` is None."""
    if keep is None:
        return frame
    kept = frame.loc[keep]
    kept.index = pd.RangeIndex(len(kept))
    return kept


def _ld_design(merged: pd.DataFrame, columns: Sequence[str], keep: np.ndarray | None) -> np.ndarray:
    """Return the LD-score design matrix of ``columns`` over the kept rows.

    A multi-column float design is filled column by column into the
    column-major array ``DataFrame.to_numpy()`` returns, without first copying
    the frame; other designs use ``to_numpy()`` on the kept rows.
    """
    columns = list(columns)
    dtypes = [merged[column].dtype for column in columns]
    if len(columns) < 2 or not all(isinstance(dtype, np.dtype) and dtype.kind == "f" for dtype in dtypes):
        return np.asarray(_kept_rows(merged, keep)[columns])
    n_rows = len(merged) if keep is None else int(keep.sum())
    design = np.empty((n_rows, len(columns)), dtype=np.result_type(*dtypes), order="F")
    for index, column in enumerate(columns):
        values = merged[column].to_numpy()
        # Assignment widens float32 columns exactly, as to_numpy() would.
        design[:, index] = values if keep is None else values[keep]
    return design


def _ld_column_variances(merged: pd.DataFrame, columns: Sequence[str], known: dict | None = None) -> pd.Series:
    """Return each column's variance, reusing ``known`` values for unchanged columns.

    Each variance is the pandas column variance, which equals the value the
    column receives within ``DataFrame.var()``.
    """
    known = known or {}
    return pd.Series({column: known[column] if column in known else merged[column].var() for column in columns})


def _allele_pair_flags(columns: Sequence[pd.Series]) -> tuple[np.ndarray, np.ndarray]:
    """Return per-row (compatible, flip) flags for concatenated allele pairs.

    The pair string is the four alleles concatenated; it is compatible when
    it is in ``MATCH_ALLELES`` and flips when ``FLIP_ALLELES`` says so. The string is built once per distinct allele
    combination with the same pandas operations a per-row build would use.
    """
    codes, sizes = [], []
    for column in columns:
        column_codes, uniques = pd.factorize(column, use_na_sentinel=True)
        codes.append(column_codes + 1)
        sizes.append(len(uniques) + 1)
    if np.prod(sizes, dtype=object) < 2 ** 62:
        combined = np.zeros(len(codes[0]), dtype=np.int64)
        for column_codes, size in zip(codes, sizes):
            combined = combined * size + column_codes
        _, first, inverse = np.unique(combined, return_index=True, return_inverse=True)
    else:
        _, first, inverse = np.unique(np.stack(codes, axis=1), axis=0, return_index=True, return_inverse=True)
    inverse = inverse.reshape(-1)
    parts = [column.iloc[first].reset_index(drop=True) for column in columns]
    pairs = parts[0] + parts[1] + parts[2] + parts[3]
    compatible = pairs.isin(reg.MATCH_ALLELES).to_numpy(dtype=bool)
    flip = np.zeros(len(pairs), dtype=bool)
    flip[compatible] = pairs[compatible].map(reg.FLIP_ALLELES).astype(bool).to_numpy()
    return compatible[inverse], flip[inverse]


def _same_strand_ambiguous_pair(merged: pd.DataFrame) -> np.ndarray:
    """Mark rg rows where both traits list the same A/T or C/G pair, in either order."""
    a1, a2, a1x, a2x = (merged[column] for column in ("A1", "A2", "A1x", "A2x"))
    same_pair = (a1.eq(a1x) & a2.eq(a2x)) | (a1.eq(a2x) & a2.eq(a1x))
    return (a1.map(reg.COMPLEMENT).eq(a2) & same_pair).to_numpy(dtype=bool, na_value=False)


def _nonfinite_columns(frame: pd.DataFrame, columns: Sequence[str], snps: pd.Series) -> list[str]:
    """Describe each column with missing or infinite values: its name, row count and up to five example SNPs."""
    problems = []
    for column in columns:
        values = frame[column]
        # Float columns are tested as stored; widening them would not change finiteness.
        values = (values.to_numpy() if isinstance(values.dtype, np.dtype) and values.dtype.kind == "f"
                  else values.to_numpy(dtype=np.float64))
        bad = ~np.isfinite(values)
        if bad.any():
            count = int(bad.sum())
            examples = ", ".join(snps.astype(str).to_numpy()[bad][:5])
            problems.append(f"'{column}' ({count} row{'' if count == 1 else 's'}; for example {examples})")
    return problems


def _select_model_columns(
    merged: pd.DataFrame,
    ref_ld_columns: Sequence[str],
    count_records: Iterable[dict[str, Any]] | dict[str, dict[str, Any]],
    use_common_counts: bool,
    *,
    command: str,
    known_variances: dict | None = None,
) -> dict[str, object]:
    """Drop zero-variance LD-score columns and select the requested count vectors.

    Each command calls this after its own SNP alignment, so variances are
    taken on the rows that enter the fit. A dropped column's values are
    constant there and identify no coefficient; its name is kept for the
    fitted-model metadata. Returns the ``RegressionDataset`` fields
    ``merged``, ``retained_ld_columns``, ``dropped_zero_variance_ld_columns``,
    ``reference_snp_count_totals`` and ``count_key_used_for_regression``.
    ``merged`` borrows the input's columns; callers must not modify them.
    ``known_variances`` supplies variances already computed on these rows.
    """
    variances = _ld_column_variances(merged, ref_ld_columns, known_variances)
    dropped = variances.index[variances == 0].tolist()
    retained = [column for column in ref_ld_columns if column not in dropped]
    if not retained:
        raise LDSCInputError(
            f"{command} regression cannot run because all retained LD-score columns have zero variance on the "
            "merged regression SNPs. Most likely the LD-score directory contains constant annotations on the "
            "retained SNP set. Use LD scores with at least one varying annotation or broaden the regression SNP set."
        )
    if dropped:
        merged = pd.DataFrame({name: merged[name].to_numpy(copy=False) for name in merged if name not in dropped},
                              copy=False)
        LOGGER.warning(
            "Dropped %d zero-variance LD-score column(s) after merging with summary statistics. "
            "These annotations are constant on the regression SNP set and receive no coefficient; "
            "names are retained in fitted-model metadata (dropped_zero_variance_ld_columns).",
            len(dropped),
        )
    count_totals = _count_totals_for_columns(count_records, ref_ld_columns)
    count_key = _select_count_key(count_totals, use_common_counts, ref_ld_columns)
    keep = [index for index, column in enumerate(ref_ld_columns) if column in retained]
    return {
        "merged": merged,
        "retained_ld_columns": retained,
        "dropped_zero_variance_ld_columns": dropped,
        "reference_snp_count_totals": {key: np.asarray(values)[keep] for key, values in count_totals.items()},
        "count_key_used_for_regression": count_key,
    }


def _select_count_key(
    count_totals: dict[str, np.ndarray],
    use_common_counts: bool,
    columns: Sequence[str] = (),
) -> str:
    """Return the requested count vector key; the other count kind is never substituted."""
    if use_common_counts:
        if COMMON_COUNT_KEY in count_totals:
            return COMMON_COUNT_KEY
        raise LDSCInputError(
            "Regression requested common-SNP counts (`--count-kind common`, the default), but the LD-score "
            "directory does not record a common-SNP count for every requested LD-score column "
            f"({_preview_names(columns)}). All-SNP counts are never substituted. Most likely the metadata.json "
            "was written by an older LDSC3 version or edited by hand. Regenerate or reconvert the LD-score "
            "directory, or pass `--count-kind all` (`RegressionConfig(use_common_counts=False)`) to scale "
            "estimates by all-SNP counts."
        )
    if ALL_COUNT_KEY in count_totals:
        return ALL_COUNT_KEY
    raise LDSCInputError(
        "Regression requested all-SNP counts, but all-SNP counts are unavailable for LD-score column(s) "
        f"{list(columns)}. This is expected for an imported unpartitioned legacy suite when one or more "
        "chromosome .M files were missing. Use the default common-SNP counts or convert a suite with complete .M files."
    )


def _preview_names(names: Sequence[str]) -> str:
    """Join the first five names for a message and count the rest."""
    names = list(names)
    preview = ", ".join(names[:5])
    return f"{preview}, ... ({len(names) - 5} more)" if len(names) > 5 else preview


def _validate_partitioned_query_columns(ldscore_result: LDScoreSource | LDScoreResult, query_columns: Sequence[str]) -> list[str]:
    """Validate that requested cell-type query columns exist in the LD-score result.

    The functional regime (no query columns) is decided by the caller; this
    helper is only used on the cell-type path, so it just checks that every
    requested query column is available.
    """
    requested = list(query_columns)
    available = list(ldscore_result.query_columns)
    available_set = set(available)
    missing = [column for column in requested if column not in available_set]
    if missing:
        raise LDSCInputError(
            f"partitioned-h2 could not find requested query annotation columns {missing} in the LD-score directory. "
            f"Most likely the query annotation name was misspelled or the LD-score directory was generated from "
            f"different query inputs. Choose one of the available query annotations: {available}."
        )
    return requested


def _assemble_regression_ldscore_table(ldscore_result: LDScoreSource | LDScoreResult) -> pd.DataFrame:
    """Return the baseline LD-score table as a shallow copy with a RangeIndex.

    Callers add columns to it but never modify the shared values in place.
    """
    table = ldscore_result.baseline_table.copy(deep=False)
    table.index = pd.RangeIndex(len(table))
    return table


def _count_totals_for_columns(count_records, columns: Sequence[str]) -> dict[str, np.ndarray]:
    """Return LDSC count vectors aligned to ``columns`` using metadata records.

    ``count_records`` is a sequence of records or a mapping from column name to record.
    """
    records_by_column = (count_records if isinstance(count_records, dict)
                         else {str(record["column"]): record for record in count_records})
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
            raise LDSCInputError(
                "Regression found non-finite all-SNP count metadata in the LD-score directory."
                " Restore consistent finite count metadata from the same LD-score run, or regenerate the "
                "LD-score directory. Do not replace missing or invalid counts with zero."
            )
        count_totals[ALL_COUNT_KEY] = all_counts
    if all("common_reference_snp_count" in records_by_column[column] for column in columns):
        common_values = [records_by_column[column].get("common_reference_snp_count") for column in columns]
        if any(value is None or pd.isna(value) for value in common_values):
            raise LDSCInputError(
                "Regression found a malformed common-SNP count record with a missing value. "
                "Reconvert or regenerate the LD-score directory; common-SNP counts never fall back silently."
            )
        common_counts = np.asarray([float(value) for value in common_values], dtype=np.float64)
        if not np.isfinite(common_counts).all():
            raise LDSCInputError(
                "Regression found non-finite common-SNP count metadata in the LD-score directory."
                " Restore consistent finite count metadata from the same LD-score run, or regenerate the "
                "LD-score directory. Do not replace missing or invalid counts with zero."
            )
        count_totals[COMMON_COUNT_KEY] = common_counts
    return count_totals


def _liability_pair(obs, obs_se, samp_prev, pop_prev):
    """Return (liab, liab_se) for an absolute SNP heritability estimate, or (nan, nan).

    The prevalences come checked from ``RegressionConfig`` or
    ``estimate_rg_pairs``, where a quantitative trait is ``(None, None)``;
    it yields NaN. ``obs`` and ``obs_se`` are the observed-scale point
    estimate and standard error.
    """
    if samp_prev is None or pop_prev is None:
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
    probabilities in (0, 1). ``intercept_se`` is NaN for a fixed intercept, and
    ``ratio``/``ratio_se`` are NaN when the mean chi-square is at most 1.
    ``n_snps_used`` is the post-chi-square-filter SNP
    count from the regression fit; when omitted it falls back
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
                "intercept_se": _numeric_scalar(getattr(hsq, "intercept_se", None), "intercept_se"),
                "mean_chisq": _scalar_or_value(hsq.mean_chisq),
                "lambda_gc": _scalar_or_value(hsq.lambda_gc),
                "ratio": _numeric_scalar(getattr(hsq, "ratio", None), "ratio"),
                "ratio_se": _numeric_scalar(getattr(hsq, "ratio_se", None), "ratio_se"),
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
    """Summarize the exact fitted SNP population into LD-score rank bins.

    Parameters
    ----------
    hsq : ldsc._regression.estimators.Hsq-like object
        Final unpartitioned fit exposing ``coef``, ``tot``, and ``intercept``.
    ld_score, chi_square, sample_size, regression_ld_score : array_like
        Equal-length per-SNP values after the h2 chi-square filter. The order is
        preserved among tied LD scores.
    reference_snp_count : float
        Reference-SNP count used by the fitted h2 model.
    max_bins : int, optional
        Maximum number of nonempty equal-count rank bins. Default is 50.

    Returns
    -------
    pandas.DataFrame
        One row per ascending LD-score bin in
        :data:`ldsc._regression.outputs.H2_REGRESSION_BIN_COLUMNS` order. Chi-square SD uses
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
            "h2 regression diagnostics require nonempty, equal-length LD-score, chi-square, "
            "sample-size, and regression-weight LD-score arrays."
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


@serial_calculation
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
    coefficient p-value, conditional annotation heritability, and the overlap
    flag), and returns the requested ``annotation_columns`` rows in the single
    lowercase ``PARTITIONED_H2_COLUMNS`` schema. The whole-model total heritability
    is carried on every row (``total_h2_obs``, constant across one fitted model)
    so each row's ``category_h2_obs`` can be read against it. Per-annotation and
    total heritability are reported on both the observed (``*_obs``) and liability
    (``*_liab``) scales; liability columns are NaN unless both prevalences are
    finite probabilities in (0, 1). Proportions, enrichment, and coefficients are
    scale-invariant. ``annotation_columns`` must be a subset of
    ``dataset.retained_ld_columns``.
    """
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
        if samp_prev is not None and pop_prev is not None
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
    dataset: RegressionDataset,
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
    the observed (``*_obs``) and liability (``*_liab``) scales. The genetic correlation
    (``rg``/``rg_se``/``z``/``p``) is scale-invariant and unconverted. A trait's
    ``h2_*_liab`` is NaN unless its prevalence is finite; ``gencov_liab`` is NaN
    only when no prevalence was supplied for the run. ``n_snps_used`` is the
    post-product-filter SNP count from the regression fit; it falls
    back to the supplied dataset's row count. Workflows supply the fitted dataset.

    ``status`` is ``'ok'``, or ``'unestimable'`` with an ``error`` naming the
    estimate that is not positive when either heritability estimate is;
    that row keeps its count and its heritability, genetic covariance and
    intercept estimates, and its ``rg``/``rg_se``/``z``/``p`` are NaN.
    """
    n_snps = int(n_snps_used) if n_snps_used is not None else int(len(dataset.merged))
    estimate, h2_out_of_bounds = _kernel_rg_estimate(rg_result)
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
        "n_snps_used": n_snps,
        **estimate,
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
        "status": "unestimable" if h2_out_of_bounds else "ok",
        "error": _h2_out_of_bounds_error(trait_1, trait_2, h2_1_obs, h2_2_obs, n_snps) if h2_out_of_bounds else "",
    }


def _kernel_rg_estimate(rg_result) -> tuple[dict[str, float], bool]:
    """Return a fitted pair's ``rg``, ``rg_se``, ``z`` and ``p``, and whether a heritability estimate is not positive.

    When either trait's heritability estimate is at most 0, the ``RG`` kernel
    ported from legacy LDSC sets ``_negative_hsq`` and stores the string
    ``'NA'`` in these four fields (legacy ``RG.summary()`` reports the pair as
    "h2 out of bounds"); they are returned as NaN with the flag set.
    """
    fields = {"rg": "rg_ratio", "rg_se": "rg_se", "z": "z", "p": "p"}
    if rg_result._negative_hsq:
        return dict.fromkeys(fields, math.nan), True
    return {column: _required_numeric_scalar(getattr(rg_result, attr), column) for column, attr in fields.items()}, False


def _h2_out_of_bounds_error(trait_1: str, trait_2: str, h2_1_obs: float, h2_2_obs: float, n_snps: int) -> str:
    """Explain an rg pair whose genetic correlation is undefined because a heritability estimate is not positive."""
    values = " and ".join(
        f"h2_{index}_obs (trait '{trait}') is {value:.3g}"
        for index, trait, value in ((1, trait_1, h2_1_obs), (2, trait_2, h2_2_obs)) if value <= 0
    )
    return (
        f"SNP heritability estimate not positive: on the pair's {n_snps:,} SNPs, {values}. "
        "The genetic correlation divides by the square root of the product of the two heritability estimates, "
        "so rg, rg_se, z and p are NaN; the heritability and genetic covariance estimates are kept. "
        "Most likely the trait has too little heritability signal on these SNPs for its sample size; "
        "a munging error, such as a wrong N or signed-statistic column, can also cause it. "
        f"Other causes & fixes: {_RG_H2_OUT_OF_BOUNDS_DOC}"
    )


def _failed_h2_row(
    *, trait_name: str, error: str, status: str, n_snps: int, samp_prev=None, pop_prev=None,
) -> pd.DataFrame:
    """Build the ``h2_per_trait`` row of a trait whose own ``rg`` heritability fit failed.

    The columns are those of :func:`summarize_total_h2` plus ``status`` and
    ``error``. ``n_snps`` is the trait's SNP count when its dataset was built,
    else 0; the estimates are NaN and the prevalences are the trait's inputs.
    """
    estimates = ("total_h2_obs", "total_h2_obs_se", "total_h2_liab", "total_h2_liab_se")
    diagnostics = ("intercept", "intercept_se", "mean_chisq", "lambda_gc", "ratio", "ratio_se")
    return pd.DataFrame([{
        "trait_name": trait_name, "n_snps": n_snps, **dict.fromkeys(estimates, math.nan),
        "samp_prev": _prev_cell(samp_prev), "pop_prev": _prev_cell(pop_prev),
        **dict.fromkeys(diagnostics, math.nan), "status": status, "error": error,
    }])


def _failed_rg_full_row(
    *, trait_1: str, trait_2: str, pair_kind: str, error: str, status: str = "failed", n_snps_used: int = 0,
    samp_prev_1=None, pop_prev_1=None, samp_prev_2=None, pop_prev_2=None,
) -> dict[str, object]:
    """Build a full rg row for a failed or unestimable pair.

    ``n_snps_used`` is the pair's SNP count when its dataset was built, else 0.
    The prevalences are the traits' inputs, which a failed pair keeps although
    it has no estimates.
    """
    row = {column: math.nan for column in RG_FULL_COLUMNS}
    row.update(
        {
            "trait_1": trait_1,
            "trait_2": trait_2,
            "n_snps_used": n_snps_used,
            "samp_prev_1": _prev_cell(samp_prev_1),
            "pop_prev_1": _prev_cell(pop_prev_1),
            "samp_prev_2": _prev_cell(samp_prev_2),
            "pop_prev_2": _prev_cell(pop_prev_2),
            "pair_kind": pair_kind,
            "status": status,
            "error": error,
        }
    )
    return row


def _concise_rg_row(full_row: dict[str, object], note: str) -> dict[str, object]:
    """Build the concise public rg row from a full row and its note: empty for an ``ok`` pair."""
    return {
        "trait_1": full_row["trait_1"],
        "trait_2": full_row["trait_2"],
        "n_snps_used": full_row["n_snps_used"],
        "rg": full_row["rg"],
        "rg_se": full_row["rg_se"],
        "p": full_row["p"],
        "note": note,
    }


def _absolute_ldscore_dir(value) -> str | None:
    """Return the absolute LD-score directory to record, resolving a relative value against the working directory.

    Saved results link to their LD-score directory, which ``quantile-h2``
    reopens from any working directory, so a relative path would be ambiguous.
    """
    return None if value is None else str(Path(normalize_path_token(str(value))).resolve())


def _ldscore_dir_of(ldscore_result: LDScoreSource | LDScoreResult) -> str | None:
    """Return the absolute LD-score directory a result was read from, or ``None`` for one built in memory."""
    source_metadata = ldscore_result.output_paths.get("metadata")
    return str(Path(source_metadata).parent.resolve()) if source_metadata else None


def _genome_build_of(dataset: RegressionDataset) -> str | None:
    """Return the genome build of a dataset's LD scores, or ``None`` when they record none."""
    return None if dataset.config_snapshot is None else getattr(dataset.config_snapshot, "genome_build", None)


def _dataset_metadata(dataset: RegressionDataset) -> dict[str, object]:
    """Record a fitted dataset's SNP identity, count vector and LD-score columns in result metadata."""
    return {
        "effective_snp_identifier": dataset.effective_snp_identifier,
        "identity_downgrade_applied": dataset.identity_downgrade_applied,
        "count_key_used_for_regression": dataset.count_key_used_for_regression,
        "retained_ld_columns": list(dataset.retained_ld_columns),
        "dropped_zero_variance_ld_columns": list(dataset.dropped_zero_variance_ld_columns),
    }


def _fit_metadata(fit: RegressionFit, snps_key: str = "n_snps") -> dict[str, object]:
    """Record a fit's SNP and block counts, chi-square cap and two-step cutoff.

    ``rg`` pairs name the SNP count ``n_snps_used``.
    """
    return {snps_key: fit.n_snps, "n_blocks_used": fit.n_blocks, "effective_chisq_max": fit.effective_chisq_max,
            "effective_two_step_cutoff": fit.effective_two_step_cutoff}


def _h2_metadata(
    sumstats_table: SumstatsTable,
    ldscore_dir: str | None,
    outcome: RegressionFit,
    config: RegressionConfig,
    dropped_nonfinite_statistic_rows: int = 0,
) -> dict[str, object]:
    """Build h2 metadata from the realized fit and the config's checked prevalences.

    ``sumstats_file`` is the resolved path the table was read from, as in
    ``partitioned-h2`` and ``rg`` metadata.
    """
    return {
        "artifact_type": "h2_result",
        "trait_name": _trait_label(sumstats_table),
        "sumstats_file": sumstats_table.source_path,
        "ldscore_dir": ldscore_dir,
        "genome_build": _genome_build_of(outcome.dataset),
        **_dataset_metadata(outcome.dataset),
        "dropped_nonfinite_statistic_rows": dropped_nonfinite_statistic_rows,
        **_fit_metadata(outcome),
        "intercept_h2_policy": _intercept_policies(config)["intercept_h2_policy"],
        "samp_prev": config.samp_prev,
        "pop_prev": config.pop_prev,
        "scale": _scale_label(config.samp_prev, config.pop_prev),
    }


def _rg_prevalence_metadata(prevalences) -> dict[str, object]:
    """Per-pair prevalences and scale from the pair's two checked ``(samp_prev, pop_prev)`` inputs."""
    (samp_prev_1, pop_prev_1), (samp_prev_2, pop_prev_2) = prevalences
    return {"samp_prev_1": samp_prev_1, "pop_prev_1": pop_prev_1, "samp_prev_2": samp_prev_2,
            "pop_prev_2": pop_prev_2, "scale": _scale_label(samp_prev_1, pop_prev_1, samp_prev_2, pop_prev_2)}


def _rg_pair_metadata(
    table_1: SumstatsTable,
    table_2: SumstatsTable,
    full_row: dict[str, object],
    config: RegressionConfig,
    pair_kind: str,
    prevalences,
    fit: dict[str, object] | None = None,
) -> dict[str, object]:
    """Build the metadata of one rg pair; ``prevalences`` holds the two traits' checked pairs.

    ``fit`` holds the pair's SNP and block counts, cap and dataset settings:
    from its ``RegressionFit``, or for a failed pair from the ``fit_summary``
    its error carries. Without it the pair failed before its dataset was
    built, and the counts are 0 with the requested cap.
    """
    metadata = {
        "trait_1": full_row["trait_1"],
        "trait_2": full_row["trait_2"],
        "source_1": table_1.source_path,
        "source_2": table_2.source_path,
        "pair_kind": pair_kind,
        "status": full_row["status"],
        "error": full_row["error"],
        "n_snps_used": 0,
        "n_blocks_used": 0,
        "effective_chisq_max": config.chisq_max,
        "effective_two_step_cutoff": None,
        **_rg_prevalence_metadata(prevalences),
        **_intercept_policies(config),
    }
    if fit is not None:
        metadata.update(fit)
    return metadata


def _intercept_policies(config: RegressionConfig) -> dict[str, str]:
    """Describe the heritability and genetic-covariance intercepts: ``free`` or ``fixed:<value>``."""
    policy = _fit_policy(config)
    return {
        f"intercept_{name}_policy": "free" if value is None else f"fixed:{float(value):g}"
        for name, value in (("h2", policy.intercept_h2), ("gencov", policy.intercept_gencov))
    }


def _require_rg_alleles(tables: Sequence[SumstatsTable]) -> None:
    """Reject traits without ``A1``/``A2``; rg needs them to align signed statistics."""
    missing = [table for table in tables if not {"A1", "A2"}.issubset(table.data.columns)]
    if missing:
        listed = ", ".join(f"{_trait_label(table)} ('{table.source_path or 'in-memory table'}')" for table in missing)
        raise LDSCInputError(
            "rg needs A1 (the effect allele, the allele the signed statistic refers to) and A2 in every summary-"
            "statistics input to align the signed statistics of the two traits, but these inputs have no A1/A2 "
            f"columns: {listed}. Most likely `ldsc munge-sumstats` was not given --a1/--a2 and could not infer "
            "them, for example because the raw allele columns were REF/ALT. Re-run `ldsc munge-sumstats` for them "
            "with --a1 <effect-allele column> --a2 <other column>. Other causes & fixes: "
            f"{_RG_ALLELES_DOC}"
        )


def _trait_label(table: SumstatsTable) -> str:
    """Return the public trait label for a summary-statistics table."""
    if table.trait_name:
        return str(table.trait_name)
    if table.source_path:
        return Path(table.source_path).name
    return "trait"


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
            "Functional-category regime: joint fit of %d baseline annotations. "
            "Focus on `enrichment` and `enrichment_p`: `enrichment > 1` means the annotation's SNPs "
            "explain a larger share of heritability than their share of SNPs (`< 1` means a smaller "
            "share); a small `enrichment_p` indicates the enrichment is significantly different from 1, "
            "i.e. significantly larger or smaller.",
            len(ldscore_result.baseline_columns),
        )


def _log_quantitative_annotation_interpretation(ldscore_result: LDScoreResult) -> None:
    """Warn once when a fitted model contains continuous annotations."""
    quantitative = [
        name for name, annotation_type in ldscore_result.annotation_types.items()
        if annotation_type == "quantitative"
    ]
    if not quantitative:
        return
    LOGGER.warning(
        "Continuous annotation(s) detected: count=%d; names are recorded in LD-score annotation_types metadata. Numerical prop_snps, prop_h2, and enrichment summaries "
        "remain visible for compatibility, but they are weighted summaries and do not have the ordinary "
        "binary-annotation interpretation. Interpret the coefficient (`coefficient`), its SE, and its coefficient zero-test; "
        "use `ldsc quantile-h2` for post-fit quantile enrichment. To estimate distinct conditional coefficients "
        "for bins, construct binary quantile annotations and refit partitioned-h2.",
        len(quantitative),
    )


def _log_effective_regression_identity(
    sumstats_tables: Sequence[SumstatsTable],
    ldscore_result: LDScoreSource | LDScoreResult,
    runner_config: GlobalConfig,
    regression_config: RegressionConfig,
) -> None:
    """Log the SNP identifier mode regression will actually use, and its source.

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
        if _is_allele_less_legacy(table):
            # Without alleles, projection gives a legacy table the base form of the LD-score directory's mode.
            mode = identity_base_mode(ld_mode)
            source = f"{_trait_label(table)}={mode} (legacy projection by rsID, without alleles)"
        elif _has_legacy_sumstats_source(table):
            # Projection gives a legacy table the LD-score directory's mode.
            mode, source = ld_mode, f"{_trait_label(table)}={ld_mode} (legacy projection onto the LD-score panel)"
        else:
            mode, source = describe(_trait_label(table), table.config_snapshot)
        modes.append(mode)
        sources.append(source)

    source_details = "; ".join(sources) if len(sources) <= 10 else f"{len(sources)} input artifacts; see input provenance and per-pair metadata"
    if len(sources) > 10 and LOGGER.isEnabledFor(logging.DEBUG):
        LOGGER.debug("Per-input SNP modes: %s", "; ".join(sources))
    effective = modes[0]
    downgrade = False
    try:
        for mode in modes[1:]:
            compat = resolve_regression_identity_mode(
                effective, mode, allow_identity_downgrade=regression_config.allow_identity_downgrade
                or any(_is_allele_less_legacy(table) for table in sumstats_tables),
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
            source_details,
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
        source_details,
    )


def _add_jackknife_model_context(error, merged, columns):
    """Attach filtered-row support and genomic spans to failed deletions.

    A failure with one more parameter than ``columns`` had a free intercept,
    whose column follows the LD-score columns. An rg pair's three fits can
    differ here: the genetic covariance intercept may be fixed while the
    heritability intercepts are free.
    """
    support = np.count_nonzero(merged[columns].to_numpy(), axis=0)
    for failure in error.failures:
        names = [*columns, *(["intercept"] if failure["n_parameters"] == len(columns) + 1 else [])]
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


def _add_panel_coordinates(merged: pd.DataFrame, ldscore_result: LDScoreSource | LDScoreResult) -> np.ndarray:
    """Replace ``_ldsc_source_row`` with the LD-score panel's coordinates and return the source rows.

    Delete-block diagnostics report genomic spans from ``_ldsc_panel_CHR`` and
    ``_ldsc_panel_POS`` even when summary statistics omit coordinates (rsID
    inputs) or carry missing placeholders. These columns never reach outputs.
    """
    source_rows = merged.pop("_ldsc_source_row").to_numpy(dtype=np.int64)
    try:
        coordinate_columns = infer_chr_pos_columns(ldscore_result.baseline_table.columns, context="LD-score diagnostics")
    except LDSCInputError:
        coordinate_columns = ()
    for coordinate, source_column in zip(("CHR", "POS"), coordinate_columns):
        merged[f"_ldsc_panel_{coordinate}"] = ldscore_result.baseline_table[source_column].iloc[source_rows].to_numpy()
    return source_rows


def _raise_on_model_collinearity(dataset: RegressionDataset, x: np.ndarray) -> None:
    """Abort the fit when the LD-score design matrix is near-collinear.

    Matches legacy LDSC's hard stop on an ill-conditioned design matrix: rather
    than warn and emit a partition the user must second-guess, collinear
    annotations raise ``LDSCInputError`` so the only results produced are
    trustworthy ones. h2 fits and each rg trait pair apply it on their own
    regression rows; it is skipped when the LD-score directory stores no
    overlap or the model has fewer than two LD-score columns.
    """
    overlap = getattr(dataset, "ldscore_overlap", None)
    if overlap is None or len(dataset.retained_ld_columns) < 2:
        return

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


@dataclass(frozen=True)
class _FitPolicy:
    """Intercepts and two-step setting of a regression run, resolved once from its config.

    ``None`` marks an intercept that the fit estimates.
    """

    intercept_h2: float | None
    intercept_gencov: float | None
    requested_two_step_cutoff: float | None

    def two_step_cutoff(self, n_ld_columns: int) -> float | None:
        """Two-step cutoff of a model with ``n_ld_columns`` LD-score columns, or ``None`` for no two-step fit.

        A requested cutoff is used as given. Otherwise a single-annotation
        model with an estimated heritability intercept uses 30, as legacy
        LDSC's ``estimate_h2`` and ``estimate_rg`` do (``sumstats.py``).
        """
        if self.requested_two_step_cutoff is not None:
            return self.requested_two_step_cutoff
        if self.intercept_h2 is None and n_ld_columns == 1:
            return 30
        return None


def _fit_policy(config: RegressionConfig) -> _FitPolicy:
    """Resolve a config's intercepts and two-step request for ``h2``, ``partitioned-h2`` and ``rg`` fits.

    ``use_intercept=False`` fixes the heritability intercept at 1 and the
    genetic-covariance intercept at 0; otherwise ``intercept_h2`` and
    ``intercept_gencov`` fix them when set.
    """
    if not config.use_intercept:
        return _FitPolicy(1, 0, config.two_step_cutoff)
    return _FitPolicy(
        None if config.intercept_h2 is None else float(config.intercept_h2),
        None if config.intercept_gencov is None else float(config.intercept_gencov),
        config.two_step_cutoff,
    )


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
            '--pop-prev; omit both, or pass nan for both, to report observed-scale estimates only.'
        ),
    )
    parser.add_argument(
        '--pop-prev', type=float, default=None, metavar='VALUE',
        help=(
            'Population case fraction, strictly between 0 and 1, for liability-scale estimates. Requires '
            '--samp-prev; omit both, or pass nan for both, to report observed-scale estimates only.'
        ),
    )


def add_h2_arguments(parser) -> None:
    """Register grouped h2 command-line options."""
    parser.prog = "ldsc h2"
    parser.formatter_class = CLIHelpFormatter
    parser.description = "Estimate heritability from munged summary statistics and LD scores."
    sections = _add_regression_sections(parser, include_h2_intercept=True)
    sections.inputs.add_argument(
        '--trait-name', default=None,
        help=(
            'Trait label used in summaries. If omitted, use the stored trait name or derive one from the '
            'input filename.'
        ),
    )
    _add_regression_output_arguments(sections.output)
    _add_scalar_prevalence_arguments(sections.scale)
    _add_regression_model_arguments(sections, include_h2_intercept=True)
    sections.logs.add_argument(
        '--log-level', default='INFO', choices=('DEBUG', 'INFO', 'WARNING', 'ERROR'),
        help=LOG_LEVEL_HELP,
    )


def add_partitioned_h2_arguments(parser) -> None:
    """Register grouped partitioned-h2 command-line options."""
    parser.prog = "ldsc partitioned-h2"
    parser.formatter_class = CLIHelpFormatter
    parser.description = "Estimate partitioned heritability with independent query models."
    sections = _add_regression_sections(parser, include_h2_intercept=True)
    sections.inputs.add_argument(
        '--trait-name', default=None,
        help=(
            'Trait label used in summaries. If omitted, use the stored trait name or derive one from the '
            'input filename.'
        ),
    )
    _add_regression_output_arguments(sections.output)
    sections.output.add_argument(
        '--continue-on-query-error', action='store_true', default=False,
        help=(
            'Skip and record queries whose model preparation, regression, or result calculation fails; '
            'publish successful fits and write diagnostics/query_status.tsv plus error tracebacks in the log. '
            'Default: off; any query error prevents publication. Shared input/output failures and runs '
            'with no successful queries still fail.'
        ),
    )
    sections.output.add_argument(
        '--summary-sort-by', default='auto', choices=('auto', *PARTITIONED_H2_SUMMARY_SORT_COLUMNS),
        help=(
            'Order partitioned_h2.tsv by the selected column. Default: auto, coefficient-p for per-query '
            'fits and category for baseline-only fits. P-values sort ascending; other numeric values sort '
            'descending.'
        ),
    )
    _add_scalar_prevalence_arguments(sections.scale)
    _add_regression_model_arguments(sections, include_h2_intercept=True)
    sections.computing.add_argument(
        '--workers', type=parse_workers, default=1, metavar='N',
        help='Maximum complete query-plus-baseline models running simultaneously; also capped by query batch size. ' + WORKERS_HELP,
    )
    sections.computing.add_argument(
        '--query-batch-size', type=parse_query_batch_size, default=100, metavar='N',
        help=(
            'Maximum number of query annotations loaded at once; each query is fitted separately with the '
            'shared baseline. Default: 100; a smaller positive value reduces memory and caps active query workers.'
        ),
    )
    sections.logs.add_argument(
        '--log-level', default='INFO', choices=('DEBUG', 'INFO', 'WARNING', 'ERROR'),
        help=LOG_LEVEL_HELP,
    )


def add_rg_arguments(parser) -> None:
    """Register grouped rg command-line options."""
    parser.prog = "ldsc rg"
    parser.formatter_class = CLIHelpFormatter
    parser.description = "Estimate genetic correlation."
    sections = _add_regression_sections(parser, include_h2_intercept=False)
    sections.inputs.add_argument(
        '--anchor-trait', default=None,
        help=(
            'Trait name or input file path to compare against every other trait. Must identify one of '
            '--sumstats-sources. If omitted, estimate every unique trait pair.'
        ),
    )
    _add_regression_output_arguments(sections.output)
    sections.output.add_argument(
        '--write-per-pair-detail', action='store_true', default=False,
        help=(
            'Also write detailed results for each tested pair under diagnostics/pairs/. Default: off; write '
            'aggregate pair and trait tables only.'
        ),
    )
    sections.scale.add_argument(
        '--samp-prev', default=None,
        help=(
            'Comma-separated sample case fractions in resolved --sumstats-sources order, one per trait; '
            'each in (0, 1) or nan for a quantitative trait. Requires --pop-prev; cannot be combined with '
            '--prevalence-manifest. If neither form is supplied, report observed-scale estimates only.'
        ),
    )
    sections.scale.add_argument(
        '--pop-prev', default=None,
        help=(
            'Comma-separated population case fractions in resolved --sumstats-sources order, one per trait; '
            'each in (0, 1) or nan. Requires --samp-prev; cannot be combined with --prevalence-manifest. If '
            'neither form is supplied, report observed-scale estimates only.'
        ),
    )
    sections.scale.add_argument(
        '--prevalence-manifest', default=None,
        help=(
            'TSV with trait_name, samp_prev, and pop_prev columns, matched to exact input trait names. '
            'Cannot be combined with --samp-prev or --pop-prev. If omitted, use paired prevalence lists '
            'when supplied; otherwise report observed-scale estimates only.'
        ),
    )
    _add_regression_model_arguments(sections, include_h2_intercept=False)
    sections.logs.add_argument(
        '--log-level', default='INFO', choices=('DEBUG', 'INFO', 'WARNING', 'ERROR'),
        help=LOG_LEVEL_HELP,
    )


@materializing_overwrite_guard(
    args_output_scope,
    command="run_h2_from_args(...)",
)
def run_h2_from_args(args):
    """Run single-trait heritability estimation from parsed CLI arguments.

    The regression settings become a ``RegressionConfig`` before anything is
    written. The workflow requires ``args.output_dir`` and preflights
    ``h2.tsv``, the LD score regression-bin diagnostic, metadata, and
    ``diagnostics/h2.log`` before loading inputs. A successful overwrite
    removes default plots and liability-scale post-processing derived from the
    superseded h2 result. Returns the :class:`H2Result` of the run.
    """
    global_config, config = _configs_from_args(args, "h2")
    return _run_h2(
        global_config, config, sumstats_file=args.sumstats_file, ldscore_dir=args.ldscore_dir,
        trait_name=args.trait_name, output_dir=args.output_dir, overwrite=args.overwrite, log_level=args.log_level,
    )


def _run_h2(
    global_config: GlobalConfig,
    config: RegressionConfig,
    *,
    sumstats_file: str,
    ldscore_dir: str,
    trait_name: str | None,
    output_dir: str,
    overwrite: bool,
    log_level: str,
) -> H2Result:
    """Estimate SNP heritability from typed settings and write the ``h2`` result directory.

    ``global_config`` is the registered configuration, the SNP identifier
    fallback for an input without recorded provenance; the caller opens the
    failure-marker scope.
    """
    output_dir, log_path = _preflight_regression_outputs(output_dir, overwrite, "h2", H2DirectoryWriter)
    with workflow_logging("h2", log_path, log_level=log_level):
        _inspect_regression_inputs(output_dir, sumstats=sumstats_file, ldscore_dir=ldscore_dir)
        runner = RegressionRunner(global_config, config)
        print_global_config_banner("run_h2_from_args", runner.global_config)
        log_inputs(sumstats_file=sumstats_file, ldscore_dir=ldscore_dir, output_dir=output_dir)
        LOGGER.info(f"Starting h2 regression for '{sumstats_file}' using LD-score directory '{ldscore_dir}'.")
        sumstats_table = _load_sumstats_table(sumstats_file, trait_name)
        ldscore_result = load_ldscore_from_dir(ldscore_dir)
        _log_effective_regression_identity([sumstats_table], ldscore_result, runner.global_config, config)
        trait = _prepare_trait_sumstats(sumstats_table, ldscore_result)
        with suppress_global_config_banner():
            dataset = runner._dataset_from_prepared(
                runner._prepare_h2_inputs(trait, _PreparedLDScores(ldscore_result), config), config=config,
            )
        # The fit needs only the aligned rows; release the trait's keyed copy.
        legacy_drops, nonfinite_rows, trait = trait.legacy_drops, trait.dropped_nonfinite_statistic_rows, None
        outcome = _fit_h2_dataset(dataset, config)
        hsq, dataset = outcome.estimator, outcome.dataset
        n_snps_used = outcome.n_snps
        summary = summarize_total_h2(
            hsq, dataset, trait_name=_trait_label(sumstats_table), n_snps_used=n_snps_used,
            samp_prev=config.samp_prev, pop_prev=config.pop_prev,
        )
        written = H2DirectoryWriter().write(
            summary,
            H2OutputConfig(
                output_dir=output_dir,
                overwrite=overwrite,
            ),
            metadata=_h2_metadata(
                sumstats_table, _ldscore_dir_of(ldscore_result), outcome, config,
                dropped_nonfinite_statistic_rows=nonfinite_rows,
            ),
            diagnostic_bins=outcome.diagnostic_bins,
            legacy_sumstats_drops=legacy_drops,
        )
        log_outputs(**written)
        LOGGER.info(f"Finished h2 regression with {n_snps_used} regression SNPs.")
    return H2Result(summary, written, legacy_drops)


@materializing_overwrite_guard(
    args_output_scope,
    command="run_partitioned_h2_from_args(...)",
)
def run_partitioned_h2_from_args(args):
    """Run batch partitioned heritability from parsed CLI arguments.

    The regression settings become a ``RegressionConfig`` before the output
    preflight; the argparse types of ``--workers`` and ``--query-batch-size``,
    or :func:`run_partitioned_h2`, have already checked those two options.
    The workflow requires ``args.output_dir`` and preflights
    ``partitioned_h2.tsv``, ``diagnostics/query_annotations/``, and
    ``diagnostics/partitioned-h2.log`` before loading inputs. Query-annotation
    runs always write per-query results; baseline-only runs keep the complete
    fitted-model artifacts at the result root. It returns the
    ``PartitionedH2BatchResult``, whose ``output_paths`` also name the legacy
    summary-statistics audit when one is written. The saved root TSV appends
    ``query_result_dir`` relative to ``args.output_dir`` (``NaN`` for baseline-only
    runs); the returned summary and per-query tables retain scientific columns
    only. A successful overwrite removes
    the default plot root derived from the superseded result.
    ``args.continue_on_query_error`` explicitly permits query regression
    failures: successful models are published and every attempted query is
    recorded in ``diagnostics/query_status.tsv``. Strict mode, the default,
    collects query errors and publishes only the diagnostic ledger on failure.
    ``args.workers`` uses the shared LD-score/index worker-count policy, with
    query count and ``args.query_batch_size`` as the work limits. The default
    and every effective count of one use one worker. Each worker fits
    one complete query model with one native numerical thread; baseline-only
    input remains one joint fit in one worker. Worker death and transport/output
    failures abort under either query-error policy. Requested/effective counts
    and phase timings are recorded in the workflow log.
    """
    global_config, config = _configs_from_args(args, "partitioned-h2")
    return _run_partitioned_h2(
        global_config, config, sumstats_file=args.sumstats_file, ldscore_dir=args.ldscore_dir,
        trait_name=args.trait_name, output_dir=args.output_dir, overwrite=args.overwrite, log_level=args.log_level,
        workers=args.workers, query_batch_size=args.query_batch_size,
        continue_on_query_error=args.continue_on_query_error, summary_sort_by=args.summary_sort_by,
    )


def _run_partitioned_h2(
    global_config: GlobalConfig,
    config: RegressionConfig,
    *,
    sumstats_file: str,
    ldscore_dir: str,
    trait_name: str | None,
    output_dir: str,
    overwrite: bool,
    log_level: str,
    workers: int,
    query_batch_size: int,
    continue_on_query_error: bool,
    summary_sort_by: str,
) -> PartitionedH2BatchResult:
    """Fit partitioned heritability from typed settings and write the ``partitioned-h2`` result directory.

    ``global_config`` is the SNP identifier fallback, as in :func:`_run_h2`;
    the caller opens the failure-marker scope.
    """
    output_dir, log_path = _preflight_regression_outputs(
        output_dir, overwrite, "partitioned-h2", PartitionedH2DirectoryWriter,
    )
    with workflow_logging("partitioned-h2", log_path, log_level=log_level):
        _inspect_regression_inputs(output_dir, sumstats=sumstats_file, ldscore_dir=ldscore_dir)
        runner = RegressionRunner(global_config, config)
        print_global_config_banner("run_partitioned_h2_from_args", runner.global_config)
        log_inputs(sumstats_file=sumstats_file, ldscore_dir=ldscore_dir, output_dir=output_dir)
        LOGGER.info(
            f"Starting partitioned-h2 regression for '{sumstats_file}' using LD-score directory '{ldscore_dir}'."
        )
        with _log_phase_timing("sumstats loading"):
            sumstats_table = _load_sumstats_table(sumstats_file, trait_name)
        with _log_phase_timing("LD-score loading"):
            ldscore_result = load_ldscore_from_dir(ldscore_dir)
        has_queries = bool(ldscore_result.query_columns)
        _log_effective_regression_identity([sumstats_table], ldscore_result, runner.global_config, config)
        _log_partitioned_h2_regime(ldscore_result, has_queries)
        _log_quantitative_annotation_interpretation(ldscore_result)
        with suppress_global_config_banner():
            result = runner.estimate_partitioned_h2_batch(
                sumstats_table, ldscore_result, config=config, output_dir=output_dir,
                query_batch_size=query_batch_size,
                workers=workers,
                continue_on_query_error=continue_on_query_error,
                overwrite=overwrite,
                summary_sort_by=summary_sort_by,
            )
        log_outputs(**result.output_paths)
        LOGGER.info(
            f"Finished partitioned-h2 regression for {len(ldscore_result.query_columns)} query annotations "
            f"and {len(result.summary)} summary rows."
        )
    return result


@materializing_overwrite_guard(
    args_output_scope,
    command="run_rg_from_args(...)",
)
def run_rg_from_args(args):
    """Run multi-trait genetic-correlation estimation from parsed CLI args.

    The regression settings become a ``RegressionConfig`` before anything is
    written. The workflow requires ``--output-dir`` and writes the rg result
    family: ``rg.tsv``, ``rg_full.tsv``, ``h2_per_trait.tsv``, optional
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
        Complete in-memory result family, with ``output_paths`` naming the
        files written. From the command line, a run whose harmonization
        dropped SNPs also prints one WARNING with the counts to stderr; the
        Python wrapper ``run_rg`` never prints.

    Raises
    ------
    ValueError
        If fewer than two summary-statistics paths resolve or ``output_dir`` is absent.
    """
    global_config, config = _configs_from_args(args, "rg")
    return _run_rg(
        global_config, config, sumstats_sources=args.sumstats_sources, ldscore_dir=args.ldscore_dir,
        anchor_trait=args.anchor_trait, samp_prev=args.samp_prev, pop_prev=args.pop_prev,
        prevalence_manifest=args.prevalence_manifest, output_dir=args.output_dir, overwrite=args.overwrite,
        write_per_pair_detail=args.write_per_pair_detail, log_level=args.log_level,
        # run_from_kwargs turns the notice off for Python callers; the command line prints it.
        emit_console_notices=getattr(args, "_emit_console_notices", True),
    )


def _run_rg(
    global_config: GlobalConfig,
    config: RegressionConfig,
    *,
    sumstats_sources: Sequence[str],
    ldscore_dir: str,
    anchor_trait: str | None,
    samp_prev: str | None,
    pop_prev: str | None,
    prevalence_manifest: str | None,
    output_dir: str,
    overwrite: bool,
    write_per_pair_detail: bool,
    log_level: str,
    emit_console_notices: bool,
) -> RgResultFamily:
    """Estimate genetic correlations from typed settings and write the ``rg`` result family.

    ``samp_prev`` and ``pop_prev`` are the comma-separated per-trait lists,
    resolved here against the loaded traits; ``config`` carries no
    prevalences. ``global_config`` is the SNP identifier fallback, as in
    :func:`_run_h2`; the caller opens the failure-marker scope.
    """
    output_dir, log_path = _preflight_regression_outputs(output_dir, overwrite, "rg", RgDirectoryWriter)
    with workflow_logging("rg", log_path, log_level=log_level):
        declared_inputs = _inspect_regression_inputs(
            output_dir, sumstats=sumstats_sources, ldscore_dir=ldscore_dir, prevalence_manifest=prevalence_manifest,
        )
        sumstats_paths = declared_inputs.files['sumstats']
        if len(sumstats_paths) < 2:
            raise LDSCUserError(
                "rg requires at least two summary-statistics inputs. Most likely `--sumstats-sources` "
                "resolved to fewer than two files. Pass two or more munged `.sumstats` artifacts."
            )
        runner = RegressionRunner(global_config, config)
        print_global_config_banner("run_rg_from_args", runner.global_config)
        log_inputs(
            sumstats_sources=f"{len(sumstats_paths)} files; source paths are recorded in diagnostics/metadata.json",
            anchor_trait=anchor_trait or "none",
            ldscore_dir=ldscore_dir,
            output_dir=output_dir,
        )
        LOGGER.info(
            f"Starting rg regression for {len(sumstats_paths)} summary-statistics files "
            f"using LD-score directory '{ldscore_dir}'."
        )
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
            samp_prev=samp_prev,
            pop_prev=pop_prev,
            manifest_path=prevalence_manifest,
            trait_names=[_trait_label(table) for table in sumstats_tables],
            trait_paths=[str(path) for path in sumstats_paths],
        )
        sumstats_tables = _disambiguate_trait_names(sumstats_tables)
        if prevalences is not None:
            LOGGER.info("Liability-scale prevalences applied to %d traits; values are recorded in diagnostics/metadata.json.", len(sumstats_tables))
            if LOGGER.isEnabledFor(logging.DEBUG):
                LOGGER.debug(
                    "Liability-scale prevalences applied: "
                    + ", ".join(
                        f"{_trait_label(table)}=(P={sp}, K={pp})"
                        for table, (sp, pp) in zip(sumstats_tables, prevalences)
                    )
                )
        anchor_index = _resolve_anchor_index(anchor_trait, sumstats_paths, sumstats_tables)
        ldscore_result = load_ldscore_from_dir(ldscore_dir)
        _log_effective_regression_identity(sumstats_tables, ldscore_result, runner.global_config, config)
        with suppress_global_config_banner():
            result = runner.estimate_rg_pairs(
                sumstats_tables,
                ldscore_result,
                anchor_index=anchor_index,
                config=config,
                prevalences=prevalences,
            )
        written = RgDirectoryWriter().write(
            result,
            RgOutputConfig(
                output_dir=output_dir,
                overwrite=overwrite,
                write_per_pair_detail=write_per_pair_detail,
            ),
        )
        log_outputs(**written)
        LOGGER.info(
            f"Finished rg regression for {len(result.rg)} trait pairs "
            f"and {len(result.h2_per_trait)} traits."
        )
    if emit_console_notices:
        _emit_rg_alignment_notice(result)
    return replace(result, output_paths=written)


def _emit_rg_alignment_notice(result: RgResultFamily) -> None:
    """Print one stderr WARNING when harmonization dropped SNPs from any rg pair.

    The workflow log already holds each pair's line; the notice exists because
    batch-scheduler users rarely read logs after a successful run.
    """
    unknown = sum(result.dropped_orientation_unknown_rows or [])
    incompatible = sum(result.dropped_incompatible_allele_rows or [])
    if not unknown and not incompatible:
        return
    print(
        f"WARNING: rg could not align {unknown + incompatible:,} SNPs between traits and dropped them: "
        f"{unknown:,} A/T or C/G SNPs whose strand, and so the sign of their statistic, is unknown (expected), "
        f"and {incompatible:,} SNPs whose two traits list different allele pairs (check the inputs if nonzero). "
        "Per-pair counts are dropped_orientation_unknown_rows and dropped_incompatible_allele_rows in "
        "diagnostics/metadata.json.",
        file=sys.stderr,
    )


def run_h2(**kwargs) -> H2Result:
    """Estimate SNP heritability and save the ``h2`` result directory.

    This is the Python form of ``ldsc h2`` and writes the same files. Keyword
    arguments are the command's option names without the leading ``--``:
    ``sumstats_file``, ``ldscore_dir`` and ``output_dir`` are required;
    ``trait_name``, ``count_kind``, ``intercept_h2``, ``two_step_cutoff``,
    ``chisq_max``, ``num_blocks``, ``allow_identity_downgrade``, ``samp_prev``,
    ``pop_prev`` and ``overwrite`` are optional and default to the command's
    defaults. Run the call under ``if __name__ == "__main__":`` because the
    fit runs in a spawned worker.

    The SNP identifier mode and genome build come from the munged summary
    statistics and the LD-score directory; the registered ``GlobalConfig``
    supplies only the level of ``diagnostics/h2.log``.

    Returns
    -------
    H2Result
        ``summary``, the one-row table saved as ``h2.tsv``; ``output_paths``,
        the files written; and ``legacy_sumstats_drops`` for a legacy input.

    Raises
    ------
    LDSCUsageError
        Before anything is written, when a keyword is unknown, a required one
        is missing, a value is outside the command's choices, a flag is not a
        bool, or ``log_level``, ``snp_identifier`` or ``genome_build`` is passed.

    Notes
    -----
    The result directory holds ``h2.tsv``, the regression-bin diagnostic and
    the metadata that ``plot_result()`` and ``convert_h2_scale()`` read.
    Existing outputs stop the run unless ``overwrite=True``; a failed
    authorized overwrite leaves ``RUN_FAILED.txt`` naming this call.
    """
    return run_from_kwargs(
        "run_h2", command_parser(add_h2_arguments), run_h2_from_args, kwargs,
        registered=get_global_config(),
    )


def run_partitioned_h2(**kwargs) -> PartitionedH2BatchResult:
    """Estimate partitioned heritability and save the ``partitioned-h2`` result directory.

    This is the Python form of ``ldsc partitioned-h2`` and writes the same
    files. Keyword arguments are the command's option names without the
    leading ``--``: ``sumstats_file``, ``ldscore_dir`` and ``output_dir`` are
    required; ``trait_name``, ``count_kind``, ``intercept_h2``,
    ``two_step_cutoff``, ``chisq_max``, ``num_blocks``,
    ``allow_identity_downgrade``, ``samp_prev``, ``pop_prev``,
    ``summary_sort_by``, ``workers``, ``query_batch_size``,
    ``continue_on_query_error`` and ``overwrite`` are optional and default to
    the command's defaults. Run the call under ``if __name__ == "__main__":``
    because fits run in spawned workers.

    An LD-score directory with query annotations gives one baseline-plus-query
    model per query annotation; one without them gives a single joint model of
    the baseline annotations. The SNP identifier mode and genome build come
    from the munged summary statistics and the LD-score directory; the
    registered ``GlobalConfig`` supplies only the level of
    ``diagnostics/partitioned-h2.log``.

    Returns
    -------
    PartitionedH2BatchResult
        ``summary`` with the scientific columns of the saved
        ``partitioned_h2.tsv``, ``query_status`` with one row per attempted
        query (empty for a baseline-only run), ``per_query_artifacts`` and
        ``output_paths``.

    Raises
    ------
    LDSCUsageError
        Before anything is written, when a keyword is unknown, a required one
        is missing, a value is outside the command's choices, a flag is not a
        bool, ``query_batch_size`` is not a positive integer, or
        ``log_level``, ``snp_identifier`` or ``genome_build`` is passed.
    ValueError
        Before anything is written, when ``workers`` is not a nonzero integer
        or exceeds the detected available CPU count.

    Notes
    -----
    Query errors stop publication unless ``continue_on_query_error=True``.
    Existing outputs stop the run unless ``overwrite=True``; a failed
    authorized overwrite leaves ``RUN_FAILED.txt`` naming this call.
    """
    # The argparse types check these two options on the command line, before any run starts. A keyword
    # call is checked here, before run_from_kwargs opens the failure-marker scope, so it writes nothing.
    if "workers" in kwargs:
        validate_workers(kwargs["workers"])
    if "query_batch_size" in kwargs:
        validate_query_batch_size(kwargs["query_batch_size"], error_type=LDSCUsageError)
    return run_from_kwargs(
        "run_partitioned_h2", command_parser(add_partitioned_h2_arguments), run_partitioned_h2_from_args, kwargs,
        registered=get_global_config(),
    )


def run_rg(**kwargs) -> RgResultFamily:
    """Estimate genetic correlation and save the ``rg`` result directory.

    This is the Python form of ``ldsc rg`` and writes the same files. Keyword
    arguments are the command's option names without the leading ``--``:
    ``sumstats_sources`` (two or more munged summary-statistics paths or
    ``*`` patterns), ``ldscore_dir`` and ``output_dir`` are required;
    ``anchor_trait``, ``count_kind``, ``intercept_h2``, ``intercept_gencov``,
    ``two_step_cutoff``, ``chisq_max``, ``num_blocks``,
    ``allow_identity_downgrade``, ``samp_prev``, ``pop_prev``,
    ``prevalence_manifest``, ``write_per_pair_detail`` and ``overwrite`` are
    optional and default to the command's defaults. Run the call under
    ``if __name__ == "__main__":`` because fits run in a spawned worker.

    ``samp_prev`` and ``pop_prev`` take one value per trait in resolved
    ``sumstats_sources`` order, as a list such as ``[0.5, float("nan")]``
    (``nan`` for a quantitative trait) or as the command's comma-separated
    string; ``prevalence_manifest`` matches prevalences by trait name instead.
    Every trait needs ``A1`` (the effect allele, the allele the signed
    statistic refers to) and ``A2``. Traits that share a label are
    relabeled ``label@<path suffix>``. The SNP identifier mode and genome
    build come from the inputs; the registered ``GlobalConfig`` supplies only
    the level of ``diagnostics/rg.log``.

    Returns
    -------
    RgResultFamily
        ``rg``, ``rg_full`` and ``h2_per_trait`` with the contents of the
        saved ``rg.tsv``, ``rg_full.tsv`` and ``h2_per_trait.tsv``, one
        ``per_pair_metadata`` record per pair, and ``output_paths`` naming
        the files written.

    Raises
    ------
    LDSCUsageError
        Before anything is written, when a keyword is unknown, a required one
        is missing, a value is outside the command's choices, a flag is not a
        bool, ``samp_prev`` or ``pop_prev`` is a single number, or
        ``log_level``, ``snp_identifier`` or ``genome_build`` is passed.

    Notes
    -----
    A pair whose fit fails keeps its row with ``status="failed"``, or
    ``status="unestimable"`` for a singular jackknife deletion or a
    heritability estimate that is not positive; the other pairs are still
    fitted. Existing outputs stop the run unless
    ``overwrite=True``; a failed authorized overwrite leaves
    ``RUN_FAILED.txt`` naming this call.
    """
    for name in ("samp_prev", "pop_prev"):
        value = kwargs.get(name)
        if isinstance(value, numbers.Real):
            raise LDSCUsageError(
                f"run_rg() cannot run: {name}={value!r} must give one value per trait, "
                f"for example {name}=[0.5, float('nan')]."
            )
        if value is not None and not isinstance(value, str):
            kwargs[name] = ",".join(str(item) for item in value)
    return run_from_kwargs(
        "run_rg", command_parser(add_rg_arguments), run_rg_from_args, kwargs,
        registered=get_global_config(),
    )


def _add_regression_sections(parser, include_h2_intercept: bool) -> SimpleNamespace:
    """Create the regression help sections in display order and add the shared inputs.

    Callers add the remaining options section by section, in display order,
    so the usage line follows the help sections. Sections that receive no
    options are omitted from the help.
    """
    parser.allow_abbrev = False
    sections = SimpleNamespace(
        inputs=parser.add_argument_group("Step 1: Munged summary statistics and LD scores"),
        output=parser.add_argument_group("Step 2: Output"),
        scale=parser.add_argument_group("Optional: liability scale"),
        model=parser.add_argument_group("Optional: model settings"),
        advanced=parser.add_argument_group("Optional: advanced estimation and SNP matching"),
        computing=parser.add_argument_group("Computing resources"),
        logs=parser.add_argument_group("Logging"),
    )
    if include_h2_intercept:
        sections.inputs.add_argument(
            '--sumstats-file', required=True, metavar='FILE',
            help=(
            'Required munged summary-statistics Parquet or legacy .sumstats(.gz) file. SNP identity must be '
            'compatible with --ldscore-dir. '
            + SCALAR_PATH_HELP
        ),
        )
    else:
        sections.inputs.add_argument(
            '--sumstats-sources', nargs='+', required=True, metavar='SOURCES',
            help=(
            "Required munged summary-statistics Parquet or legacy .sumstats(.gz) files; accepts exact paths or quoted '*' patterns "
            'matching filename text and requires at least two files. SNP identity must be compatible with '
            "--ldscore-dir. '@' is not expanded."
        ),
        )
    sections.inputs.add_argument(
        '--ldscore-dir', required=True, metavar='DIR',
        help=(
            'Required LD-score directory from ldscore or convert-legacy-ldscores, containing '
            'predictors, regression weights, and SNP counts.'
        ),
    )
    return sections


def _add_regression_output_arguments(output) -> None:
    """Add the regression output directory and overwrite policy."""
    output.add_argument(
        '--output-dir', required=True, metavar='DIR', type=output_directory_argument,
        help=(
            'Required destination for regression estimates and diagnostics.'
        ),
    )
    output.add_argument(
        '--overwrite', action='store_true', default=False,
        help=(
            "Replace this command's existing output files and remove obsolete outputs from an earlier run. "
            'Default: off; stop if output files already exist.'
        ),
    )


def _add_regression_model_arguments(sections: SimpleNamespace, include_h2_intercept: bool) -> None:
    """Add SNP-count, intercept, chi-square, jackknife, and SNP-matching settings."""
    sections.model.add_argument(
        '--count-kind', choices=('common', 'all'), default='common',
        help=(
            'Reference SNP counts used to scale regression estimates: common uses common-SNP counts, '
            'all uses all-SNP counts. Default: common. The requested counts must be in the LD-score '
            'directory; the other kind is never substituted.'
        ),
    )
    if include_h2_intercept:
        sections.advanced.add_argument(
            '--intercept-h2', type=float, default=None, metavar='VALUE',
            help=(
                'Fix the heritability intercept to this value; use 1 for the standard fixed intercept. '
                'Cannot be combined with an explicit --two-step-cutoff. If omitted, estimate the intercept.'
            ),
        )
    else:
        sections.advanced.add_argument(
            '--intercept-h2', type=float, default=None, metavar='VALUE',
            help=(
                'Fix the heritability intercept to the same value for both traits in every pair; use 1 for '
                'the standard fixed intercept. Cannot be combined with an explicit --two-step-cutoff. '
                'If omitted, estimate the intercepts.'
            ),
        )
        sections.advanced.add_argument(
            '--intercept-gencov', type=float, default=None, metavar='VALUE',
            help=(
                'Fix the genetic-covariance intercept to the same value for every trait pair; use 0 for '
                'the standard fixed intercept. Cannot be combined with --two-step-cutoff. Single-annotation '
                'fits also require --intercept-h2 because automatic two-step estimation otherwise conflicts '
                'with this fixed value. If omitted, estimate the covariance intercepts.'
            ),
        )
    if include_h2_intercept:
        sections.advanced.add_argument(
            '--two-step-cutoff', type=float, default=None, metavar='VALUE',
            help=(
                'First-step chi-square cutoff, including the boundary. Requires a single-annotation model with '
                'a free intercept; cannot be combined with --intercept-h2. If omitted, use 30 '
                'for such models; otherwise disable two-step estimation.'
            ),
        )
        sections.advanced.add_argument(
            '--chisq-max', type=float, default=None, metavar='VALUE',
            help=(
                'Keep SNPs with chi-square at or below this value. If omitted, use max(0.001 * maximum N, 80) '
                'for multi-annotation models; single-annotation models have no chi-square cap.'
            ),
        )
    else:
        sections.advanced.add_argument(
            '--two-step-cutoff', type=float, default=None, metavar='VALUE',
            help=(
                'First-step chi-square cutoff, including the boundary. Requires a single-annotation model with '
                'free intercepts; cannot be combined with --intercept-h2 or '
                '--intercept-gencov. If omitted, use 30 with a free heritability intercept in single-annotation fits; '
                'otherwise disable it.'
            ),
        )
        sections.advanced.add_argument(
            '--chisq-max', type=float, default=None, metavar='VALUE',
            help=(
                'Keep SNPs satisfying Z1^2 * Z2^2 <= VALUE^2 for each trait pair. If omitted, apply no '
                'chi-square product cap. This differs from the single-trait chi-square filter.'
            ),
        )
    sections.advanced.add_argument(
        '--num-blocks', type=int, default=200, metavar='N',
        help=(
            'Number of SNP blocks for jackknife standard errors; must be at least 2. Default: 200, capped '
            'at the number of retained SNPs.'
        ),
    )
    sections.advanced.add_argument(
        '--allow-identity-downgrade', action='store_true', default=False,
        help=(
            'Allow allele-aware and allele-unaware inputs of the same mode family by matching without '
            'alleles. rsID-family and coordinate-family modes still cannot be mixed. Default: off; require '
            'compatible SNP identifier modes.'
        ),
    )


def _preflight_regression_outputs(
    output_dir: str | None,
    overwrite: bool,
    workflow_name: str,
    writer,
) -> tuple[str, Path]:
    """Preflight regression outputs and return normalized output dir plus log path."""
    if not output_dir:
        raise LDSCUsageError(
            f"{workflow_name} requires `--output-dir`. Most likely a programmatic namespace omitted output_dir. "
            "Pass an explicit result directory."
        )
    output_path = ensure_output_directory(output_dir, label="output directory")
    log_path = output_path / "diagnostics" / f"{workflow_name}.log"
    writer.artifact_family(output_path).preflight(
        overwrite=overwrite,
        additional_paths=[log_path],
    )
    log_path.parent.mkdir(parents=True, exist_ok=True)
    return str(output_path), log_path


def _inspect_regression_inputs(output_dir, *, sumstats, ldscore_dir, prevalence_manifest=None):
    """Check the declared summary statistics, prevalence manifest and LD-score directory before any fit.

    ``sumstats`` is the ``--sumstats-file`` path or the ``--sumstats-sources``
    tokens. Every issue found is recorded in ``diagnostics/input_issues.tsv``
    before the run stops; otherwise returns the resolved inputs.
    """
    return inspect_declared_inputs(
        files=[('sumstats', sumstats), ('prevalence manifest', prevalence_manifest)],
        checks=[('LD scores', ldscore_dir, lambda: inspect_artifact_paths(ldscore_dir))],
        issues_path=Path(output_dir)/'diagnostics/input_issues.tsv')


def _distinguishing_path_suffixes(paths: Sequence[Path]) -> list[str]:
    """Return the shortest trailing path parts that tell distinct paths apart.

    Parent directories are tried first, so files sharing a name show only their
    directories. Files in the same directory fall back to suffixes ending in the
    file name. One suffix length is used for the whole group, which keeps the
    labels independent of input order.
    """
    for parts in ([path.parent.parts for path in paths], [path.parts for path in paths]):
        for depth in range(1, max(map(len, parts)) + 1):
            suffixes = [Path(*item[-depth:]).as_posix() for item in parts]
            if len(set(suffixes)) == len(suffixes):
                return suffixes
    return [path.as_posix() for path in paths]


def _disambiguate_trait_names(tables: Sequence[SumstatsTable]) -> list[SumstatsTable]:
    """Return tables with deterministic unique trait names for rg outputs.

    A label shared by several inputs becomes ``label@<path suffix>``, using the
    shortest source-path suffix that separates that group.
    """
    labels = [_trait_label(table) for table in tables]
    unique = list(labels)
    for label in dict.fromkeys(labels):
        members = [index for index, item in enumerate(labels) if item == label]
        if len(members) == 1:
            continue
        if any(not tables[index].source_path for index in members):
            raise LDSCInputError(
                f"rg cannot give each input a unique trait label: '{label}' labels {len(members)} inputs, and "
                "at least one has no source path to tell them apart. Give each SumstatsTable a distinct "
                "trait_name."
            )
        suffixes = _distinguishing_path_suffixes([Path(tables[index].source_path) for index in members])
        for index, suffix in zip(members, suffixes):
            unique[index] = f"{label}@{suffix}"
    collisions = [label for label in dict.fromkeys(unique) if unique.count(label) > 1]
    if collisions:
        sources = ", ".join(f"'{table.source_path}'" for table, label in zip(tables, unique) if label == collisions[0])
        raise LDSCInputError(
            f"rg cannot give each input a unique trait label: '{collisions[0]}' would label {sources}. "
            "Most likely one trait name contains '@' and equals the label generated to separate "
            "inputs that share a trait name. Re-run `ldsc munge-sumstats` for one of these inputs "
            "with a distinct --trait-name."
        )

    renamed = sum(original != resolved for original, resolved in zip(labels, unique))
    if renamed:
        LOGGER.info("Renamed %d duplicate rg trait labels; resolved names and source paths are recorded in diagnostics/metadata.json.", renamed)
    for original, resolved in zip(labels, unique):
        if original != resolved:
            LOGGER.debug("Disambiguated duplicate rg trait name '%s' as '%s'.", original, resolved)
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


def _configs_from_args(args, command: str) -> tuple[GlobalConfig, RegressionConfig]:
    """Return the ``GlobalConfig`` and ``RegressionConfig`` of parsed ``command`` arguments.

    ``command`` is ``"h2"``, ``"partitioned-h2"`` or ``"rg"``. Regression
    reads the SNP identifier mode and genome build from its saved inputs; the
    registered ``GlobalConfig`` is the fallback for an input without recorded
    provenance. The scalar ``--samp-prev``/``--pop-prev`` of ``h2`` and
    ``partitioned-h2`` go on the config, which validates them. The per-trait
    lists of ``rg`` are resolved against the loaded traits in :func:`_run_rg`,
    so its config has no prevalences.
    """
    rg = command == "rg"
    config = RegressionConfig(
        num_blocks=args.num_blocks,
        use_common_counts=(args.count_kind == "common"),
        intercept_h2=args.intercept_h2,
        intercept_gencov=args.intercept_gencov if rg else None,
        two_step_cutoff=args.two_step_cutoff,
        chisq_max=args.chisq_max,
        samp_prev=None if rg else args.samp_prev,
        pop_prev=None if rg else args.pop_prev,
        allow_identity_downgrade=args.allow_identity_downgrade,
    )
    return get_global_config(), config


def _load_sumstats_table(path: str, trait_name: str | None) -> SumstatsTable:
    """Load one munged summary-statistics file through the public workflow helper."""
    return load_sumstats(path, trait_name=trait_name)
