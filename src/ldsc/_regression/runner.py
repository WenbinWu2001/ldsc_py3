"""The public ``RegressionRunner``: in-memory h2, partitioned-h2 and rg calls.

``RegressionRunner`` holds a ``GlobalConfig`` (the SNP identifier fallback
for an input without recorded provenance) and a default
``RegressionConfig``; each method takes an optional ``config`` that
overrides the default for that call. Its methods are entry points: each one
that reads inputs logs the ``GlobalConfig`` once, and all of them forward to
the regression modules, which take both configs explicitly and log no
banner:

- ``alignment``: trait preparation, including the legacy ``.sumstats(.gz)``
  projection onto the LD-score panel, SNP identity resolution and keying,
  the genome-build check, and h2 and rg dataset assembly.
- ``fitting``: the fit policy, the h2, rg and per-query partitioned fits,
  result rows and metadata.
- ``partitioned`` and ``rg``: the partitioned-h2 batch, which publishes its
  result directory, and the rg pair run.
- ``commands``: the ``ldsc h2``, ``ldsc partitioned-h2`` and ``ldsc rg``
  entry points and their Python forms ``run_h2``, ``run_partitioned_h2``
  and ``run_rg``, which call the same modules.

Standalone fits run in one controlled child with one BLAS thread; batch
partitioned fits run each query model in a query worker.
"""

from __future__ import annotations

from typing import Sequence

import pandas as pd

from .._ldscore.results import LDScoreResult
from .._ldscore.source import LDScoreSource
from .._runtime.serial_calculation import serial_calculation
from .._sumstats.table import SumstatsTable
from ..config import GlobalConfig, RegressionConfig, get_global_config, print_global_config_banner
from . import partitioned, rg
from .results import PartitionedH2BatchResult, RegressionDataset, RegressionFit, RgResultFamily
from .alignment import _prepare_h2_dataset, _prepare_rg_dataset, _validate_partitioned_query_columns
from .fitting import _config_prevalence, _fit_h2_dataset, _fit_rg_dataset, summarize_partitioned_h2


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
        print_global_config_banner("RegressionRunner", self.global_config)
        return _prepare_h2_dataset(sumstats_table, ldscore_result, config or self.regression_config,
                                   self.global_config, query_columns=query_columns).dataset

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
        print_global_config_banner("RegressionRunner", self.global_config)
        return _prepare_rg_dataset(sumstats_table_1, sumstats_table_2, ldscore_result,
                                   config or self.regression_config, self.global_config)

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
        print_global_config_banner("RegressionRunner", self.global_config)
        return partitioned.estimate_partitioned_h2_batch(
            sumstats_table, ldscore_result, config or self.regression_config, self.global_config,
            output_dir=output_dir, query_columns=query_columns, query_batch_size=query_batch_size,
            overwrite=overwrite, summary_sort_by=summary_sort_by, metadata=metadata,
            continue_on_query_error=continue_on_query_error, workers=workers,
        )

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
        return _fit_rg_dataset(dataset, config or self.regression_config)

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
        print_global_config_banner("RegressionRunner", self.global_config)
        return rg.estimate_rg_pairs(
            sumstats_tables, ldscore_result, config or self.regression_config, self.global_config,
            anchor_index=anchor_index, prevalences=prevalences,
        )
