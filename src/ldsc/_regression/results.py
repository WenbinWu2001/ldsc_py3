"""Result and dataset types of the regression workflows.

``RegressionDataset`` is the merged, regression-ready table one fit reads;
``RegressionFit`` is one fitted regression. ``H2Result``,
``PartitionedH2BatchResult`` and ``RgResultFamily`` are what ``h2``,
``partitioned-h2`` and ``rg`` return, and ``PartitionedH2FitArtifacts`` names
the files of one fitted ``partitioned-h2`` model. The workflow modules build
these values and ``outputs`` writes them; this module imports neither.
``import ldsc`` loads it through the writers, so the estimator classes are
named for type checkers only, which keeps SciPy out of package start-up.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

from ..config import GlobalConfig
from ..errors import LDSCInternalError

if TYPE_CHECKING:
    from . import estimators as reg


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
class PartitionedH2FitArtifacts:
    """Paths to one complete fitted model and its publication ownership.

    ``owned=False`` preserves persistent files by copying them at publication.
    ``owned=True`` transfers the full table and delete values from disposable
    staging on the output filesystem; those source paths are consumed even if
    later publication fails. The caller remains responsible for scratch cleanup.
    """

    full: Path
    delete_values: Path
    metadata: Path
    owned: bool = False


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
