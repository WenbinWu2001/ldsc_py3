"""regression_runner.py

Core functionality:
    Load normalized LD-score artifacts, assemble regression-ready datasets, and
    dispatch heritability or genetic-correlation estimators.

Overview
--------
This module is the workflow-layer boundary for the LDSC regression commands.
It consumes the canonical LD-score result directory introduced by the IO
refactor: ``manifest.json`` plus split baseline/query parquet tables. Public
callers pass one ``--ldscore-dir`` input instead of individual LD-score,
weight, count, and annotation-manifest files.

Regression chooses either the manifest ``common_reference_snp_counts`` vector or
``all_reference_snp_counts`` through ``--count-kind common|all``. The default
``common`` mode falls back to all-SNP counts when common counts are unavailable.
Partitioned-h2 requires non-empty query LD-score columns and fits one
baseline-plus-query model per query annotation. It can optionally retain and
write the full category table for each fitted model through the output-layer
``PartitionedH2DirectoryWriter``. Genetic correlation accepts a list of two or
more munged summary-statistic sources and returns the full rg output family:
the concise headline table, a diagnostic full table, per-trait h2 summaries,
and optional per-pair metadata for filesystem detail outputs.

Regression commands create per-run logs only when an ``output_dir`` is supplied:
``h2.log``, ``partitioned-h2.log``, or ``rg.log``. In-memory regression calls
and CLI invocations without an output directory remain console-only.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
import math
import warnings
from dataclasses import dataclass, replace
from pathlib import Path
from types import SimpleNamespace
from typing import Any, Sequence

import numpy as np
import pandas as pd
from scipy import stats

from .config import (
    GlobalConfig,
    RegressionConfig,
    get_global_config,
    print_global_config_banner,
    suppress_global_config_banner,
    validate_config_compatibility,
)
from .path_resolution import (
    ensure_output_directory,
    ensure_output_paths_available,
    normalize_path_token,
    preflight_output_artifact_family,
    resolve_file_group,
)
from ._logging import log_inputs, log_outputs, workflow_logging
from ._kernel import regression as reg
from ._kernel.identifiers import build_snp_id_series
from ._row_alignment import assert_same_snp_rows
from .column_inference import infer_chr_pos_columns, normalize_snp_identifier_mode
from .ldscore_calculator import LDScoreResult
from .outputs import (
    LDSCORE_RESULT_FORMAT,
    PartitionedH2DirectoryWriter,
    PartitionedH2OutputConfig,
    REGRESSION_LD_SCORE_COLUMN,
    RG_CONCISE_COLUMNS,
    RG_FULL_COLUMNS,
    RgDirectoryWriter,
    RgOutputConfig,
)
from .sumstats_munger import SumstatsTable, load_sumstats


COMMON_COUNT_KEY = "common_reference_snp_counts"
ALL_COUNT_KEY = "all_reference_snp_counts"
CHR_POS_KEY_COLUMN = "_ldsc_chr_pos_key"
LOGGER = logging.getLogger("LDSC.regression_runner")
PARTITIONED_H2_AGGREGATE_COLUMNS = [
    "Category",
    "Prop._SNPs",
    "Prop._h2",
    "Enrichment",
    "Enrichment_p",
    "Coefficient",
    "Coefficient_p",
]
PARTITIONED_H2_FULL_COLUMNS = [
    "Category",
    "Prop._SNPs",
    "Category_h2",
    "Category_h2_std_error",
    "Prop._h2",
    "Prop._h2_std_error",
    "Enrichment",
    "Enrichment_std_error",
    "Enrichment_p",
    "Coefficient",
    "Coefficient_std_error",
    "Coefficient_p",
]
PARTITIONED_H2_REQUIRES_QUERY_ANNOTATIONS_MESSAGE = (
    "partitioned-h2 requires query annotations in --ldscore-dir. "
    "Rerun `ldsc ldscore` with --query-annot-sources or --query-annot-bed-sources plus explicit baseline annotations."
)
FAILED_RG_NOTE = "Failed; see rg_full.tsv error column; use --output-dir for rg.log."


@dataclass(frozen=True)
class RegressionDataset:
    """Merged regression-ready dataset built from sumstats and LD-score tables."""
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

    def validate(self) -> None:
        """Validate that the merged table contains the required LDSC columns."""
        required = {"SNP", self.weight_column, "Z", "N"}
        missing = required - set(self.merged.columns)
        if missing:
            raise ValueError(f"RegressionDataset is missing required columns: {sorted(missing)}")


@dataclass(frozen=True)
class RGRegressionDataset:
    """Merged genetic-correlation dataset built from two traits and LD scores."""
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

    def validate(self) -> None:
        """Validate that the merged table contains the required RG columns."""
        required = {"SNP", self.weight_column, "Z1", "N1", "Z2", "N2"}
        missing = required - set(self.merged.columns)
        if missing:
            raise ValueError(f"RGRegressionDataset is missing required columns: {sorted(missing)}")


@dataclass(frozen=True)
class PartitionedH2BatchResult:
    """Batch partitioned-h2 summaries plus optional per-query detail tables.

    ``summary`` is the compact one-row-per-query table written as
    ``partitioned_h2.tsv``. The optional detail dictionaries are keyed by
    original query annotation name and are populated only when callers request
    the full baseline-plus-query tables that are written as
    ``partitioned_h2_full.tsv``.
    """
    summary: pd.DataFrame
    per_query_category_tables: dict[str, pd.DataFrame]
    per_query_metadata: dict[str, dict[str, object]]


@dataclass(frozen=True)
class RgResultFamily:
    """Result family returned by multi-trait genetic-correlation workflows.

    Parameters
    ----------
    rg : pandas.DataFrame
        Concise one-row-per-pair table with the publication-oriented rg schema:
        trait names, SNP count, :math:`r_g`, standard error, p-value,
        BH-adjusted p-value, and a failure note.
    rg_full : pandas.DataFrame
        Comprehensive diagnostic table with one row per attempted pair. Failed
        pairs keep their row, set numeric fields to NaN, and store
        ``status='failed'`` plus an ``error`` message.
    h2_per_trait : pandas.DataFrame
        One total-heritability summary per input trait, computed once per trait
        on the trait's single-trait LDSC regression dataset.
    per_pair_metadata : list of dict
        Ordered metadata records aligned to ``rg_full`` rows. Writers use these
        records for optional ``pairs/`` detail output.
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

    def build_dataset(
        self,
        sumstats_table: SumstatsTable,
        ldscore_result: LDScoreResult,
        config: RegressionConfig | None = None,
        query_columns: Sequence[str] | None = None,
    ) -> RegressionDataset:
        """Merge sumstats, LD scores, and weights into a regression dataset.

        If both inputs carry known ``GlobalConfig`` snapshots, their critical
        settings are checked before merging. Unknown-provenance inputs, such as
        disk-loaded sumstats, are allowed through this compatibility boundary.
        In ``rsid`` mode the merge uses the literal ``SNP`` column. In
        ``chr_pos`` mode it uses normalized private ``CHR:POS`` keys and keeps
        the original sumstats ``SNP`` value as metadata.

        Zero-variance LD-score columns are dropped here so the estimator kernel
        receives only informative regressors. The selected count vector is
        carried alongside the merged table for later use by ``Hsq``.
        """
        print_global_config_banner(type(self).__name__, self.global_config)
        config = config or self.regression_config
        if sumstats_table.config_snapshot is not None and ldscore_result.config_snapshot is not None:
            validate_config_compatibility(
                sumstats_table.config_snapshot,
                ldscore_result.config_snapshot,
                context="SumstatsTable and LDScoreResult",
            )
        weight_column = REGRESSION_LD_SCORE_COLUMN
        selected_query_columns = list(query_columns or [])
        ref_ld_columns = list(ldscore_result.baseline_columns) + selected_query_columns
        ldscore_frame = _assemble_regression_ldscore_table(ldscore_result, selected_query_columns)
        identifier_mode = _effective_snp_identifier_mode(sumstats_table, ldscore_result, self.global_config)
        if identifier_mode == "rsid":
            merged = pd.merge(
                sumstats_table.data,
                ldscore_frame.loc[:, ["SNP", *ref_ld_columns, weight_column]].reset_index(drop=True),
                how="inner",
                on="SNP",
                sort=False,
            )
        else:
            sumstats_keyed = _with_chr_pos_key(sumstats_table.data, context=sumstats_table.source_path or "sumstats")
            ldscore_keyed = _with_chr_pos_key(ldscore_frame, context="LD-score table")
            merged = pd.merge(
                sumstats_keyed,
                ldscore_keyed.loc[:, [CHR_POS_KEY_COLUMN, *ref_ld_columns, weight_column]].reset_index(drop=True),
                how="inner",
                on=CHR_POS_KEY_COLUMN,
                sort=False,
            )
        if merged.empty:
            source = sumstats_table.source_path or sumstats_table.trait_name or "sumstats"
            raise ValueError(
                f"No overlapping {identifier_mode} SNPs remain after merging sumstats '{source}' "
                f"with {len(ldscore_frame)} LD-score rows. Check that snp_identifier and genome_build match. "
                f"Active config: {self.global_config!r}."
            )

        retained_ld_columns = list(ref_ld_columns)
        dropped_ld_columns: list[str] = []
        if retained_ld_columns:
            variances = merged.loc[:, retained_ld_columns].var()
            dropped_ld_columns = variances.index[variances == 0].tolist()
            retained_ld_columns = [column for column in retained_ld_columns if column not in dropped_ld_columns]
            if retained_ld_columns:
                merged = merged.loc[:, [column for column in merged.columns if column not in dropped_ld_columns]]
            else:
                raise ValueError("All LD-score columns have zero variance.")

        count_totals = _count_totals_for_columns(ldscore_result.count_records, ref_ld_columns)
        count_key = _select_count_key(count_totals, config.use_common_counts)
        if dropped_ld_columns:
            dropped_index = [ref_ld_columns.index(column) for column in dropped_ld_columns]
            keep_index = [idx for idx in range(len(ref_ld_columns)) if idx not in dropped_index]
            for key, values in list(count_totals.items()):
                count_totals[key] = np.asarray(values)[keep_index]

        dataset = RegressionDataset(
            merged=merged.reset_index(drop=True),
            ref_ld_columns=ref_ld_columns,
            weight_column=weight_column,
            reference_snp_count_totals=count_totals,
            count_key_used_for_regression=count_key,
            retained_ld_columns=retained_ld_columns,
            dropped_zero_variance_ld_columns=dropped_ld_columns,
            trait_names=[name for name in [sumstats_table.trait_name] if name],
            chromosomes_aggregated=[result.chrom for result in ldscore_result.chromosome_results],
            config_snapshot=ldscore_result.config_snapshot,
        )
        dataset.validate()
        return dataset

    def build_rg_dataset(
        self,
        sumstats_table_1: SumstatsTable,
        sumstats_table_2: SumstatsTable,
        ldscore_result: LDScoreResult,
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
        for table, label in ((sumstats_table_1, "trait 1 SumstatsTable"), (sumstats_table_2, "trait 2 SumstatsTable")):
            if table.config_snapshot is not None and ldscore_result.config_snapshot is not None:
                validate_config_compatibility(
                    table.config_snapshot,
                    ldscore_result.config_snapshot,
                    context=f"{label} and LDScoreResult",
                )

        weight_column = REGRESSION_LD_SCORE_COLUMN
        ref_ld_columns = list(ldscore_result.baseline_columns)
        ldscore_frame = _assemble_regression_ldscore_table(ldscore_result, [])
        identifier_mode = _effective_snp_identifier_mode(sumstats_table_1, ldscore_result, self.global_config)

        left = sumstats_table_1.data.rename(columns={"N": "N1", "Z": "Z1"})
        right = sumstats_table_2.data.rename(columns={"A1": "A1x", "A2": "A2x", "N": "N2", "Z": "Z2"})
        right_payload = [column for column in ["A1x", "A2x", "N2", "Z2"] if column in right.columns]
        if identifier_mode == "rsid":
            left_with_ld = pd.merge(
                left,
                ldscore_frame.loc[:, ["SNP", *ref_ld_columns, weight_column]].reset_index(drop=True),
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
                left_keyed,
                ldscore_keyed.loc[:, [CHR_POS_KEY_COLUMN, *ref_ld_columns, weight_column]].reset_index(drop=True),
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
        merged = merged.dropna(how="any").reset_index(drop=True)
        if merged.empty:
            raise ValueError(
                f"No overlapping {identifier_mode} SNPs remain after merging both sumstats tables "
                f"with {len(ldscore_frame)} LD-score rows."
            )

        if {"A1", "A2", "A1x", "A2x"}.issubset(merged.columns):
            alleles = merged["A1"] + merged["A2"] + merged["A1x"] + merged["A2x"]
            keep = reg._filter_alleles(alleles)
            kept_alleles = alleles.loc[keep].reset_index(drop=True)
            merged = merged.loc[keep].reset_index(drop=True)
            if merged.empty:
                raise ValueError("No allele-compatible SNPs remain after harmonizing the two sumstats tables.")
            merged["Z2"] = reg._align_alleles(merged["Z2"].copy(), kept_alleles)

        retained_ld_columns = list(ref_ld_columns)
        dropped_ld_columns: list[str] = []
        if retained_ld_columns:
            variances = merged.loc[:, retained_ld_columns].var()
            dropped_ld_columns = variances.index[variances == 0].tolist()
            retained_ld_columns = [column for column in retained_ld_columns if column not in dropped_ld_columns]
            if retained_ld_columns:
                merged = merged.loc[:, [column for column in merged.columns if column not in dropped_ld_columns]]
            else:
                raise ValueError("All LD-score columns have zero variance.")

        count_totals = _count_totals_for_columns(ldscore_result.count_records, ref_ld_columns)
        if dropped_ld_columns:
            dropped_index = [ref_ld_columns.index(column) for column in dropped_ld_columns]
            keep_index = [idx for idx in range(len(ref_ld_columns)) if idx not in dropped_index]
            for key, values in list(count_totals.items()):
                count_totals[key] = np.asarray(values)[keep_index]
        count_key = _select_count_key(count_totals, config.use_common_counts)

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
        )
        dataset.validate()
        return dataset

    def estimate_h2(
        self,
        dataset: RegressionDataset,
        config: RegressionConfig | None = None,
    ):
        """Estimate single-trait heritability from a prepared dataset."""
        config = config or self.regression_config
        merged = dataset.merged
        n_snp = len(merged)
        n_blocks = min(n_snp, config.n_blocks)
        x = np.asarray(merged[dataset.retained_ld_columns])
        chisq = np.asarray(merged["Z"] ** 2).reshape((n_snp, 1))
        if config.chisq_max is not None:
            keep = np.ravel(chisq < config.chisq_max)
            merged = merged.loc[keep].reset_index(drop=True)
            n_snp = len(merged)
            n_blocks = min(n_snp, config.n_blocks)
            x = np.asarray(merged[dataset.retained_ld_columns])
            chisq = np.asarray(merged["Z"] ** 2).reshape((n_snp, 1))
        intercept = None
        if not config.use_intercept:
            intercept = 1
        elif config.intercept_h2 is not None and not isinstance(config.intercept_h2, list):
            intercept = float(config.intercept_h2)
        two_step = config.two_step_cutoff
        if two_step is None and intercept is None and len(dataset.retained_ld_columns) == 1:
            two_step = 30
        old_weights = len(dataset.retained_ld_columns) > 1
        return reg.Hsq(
            chisq,
            x,
            np.asarray(merged[[dataset.weight_column]]),
            np.asarray(merged[["N"]]),
            np.asarray(dataset.reference_snp_count_totals[dataset.count_key_used_for_regression]).reshape((1, -1)),
            n_blocks=n_blocks,
            intercept=intercept,
            twostep=two_step,
            old_weights=old_weights,
        )

    def estimate_partitioned_h2(
        self,
        sumstats_table: SumstatsTable,
        ldscore_result: LDScoreResult,
        *,
        query_column: str,
        config: RegressionConfig | None = None,
    ) -> pd.DataFrame:
        """Estimate partitioned heritability for one query annotation.

        Parameters
        ----------
        sumstats_table : SumstatsTable
            Munged single-trait summary statistics.
        ldscore_result : LDScoreResult
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
        return summarize_partitioned_h2(hsq, dataset, selected_queries)

    def estimate_partitioned_h2_batch(
        self,
        sumstats_table: SumstatsTable,
        ldscore_result: LDScoreResult,
        annotation_bundle,
        config: RegressionConfig | None = None,
        include_full_partitioned_h2: bool = False,
        include_model_categories: bool | None = None,
    ) -> pd.DataFrame | PartitionedH2BatchResult:
        """Estimate one baseline-plus-query model per query annotation.

        Partitioned-h2 is only defined here for explicit query annotations.
        Baseline-only LD-score directories should be analyzed with
        :meth:`estimate_h2` or :meth:`estimate_rg`.

        Parameters
        ----------
        sumstats_table : SumstatsTable
            Munged single-trait summary statistics.
        ldscore_result : LDScoreResult
            Canonical LD-score result containing baseline and query columns.
        annotation_bundle : object
            Object exposing ordered ``query_columns``. Disk-loaded CLI runs use
            a lightweight namespace derived from the LD-score manifest.
        config : RegressionConfig, optional
            Regression settings. Defaults to the runner's config.
        include_full_partitioned_h2 : bool, optional
            If ``True``, return ``PartitionedH2BatchResult`` with per-query
            full category tables for output writers. If ``False``, return the
            aggregate dataframe.
        include_model_categories : bool, optional
            Backward-compatible alias for ``include_full_partitioned_h2``.

        Returns
        -------
        pandas.DataFrame or PartitionedH2BatchResult
            Aggregate summary by default, or aggregate plus per-query detail
            tables when ``include_full_partitioned_h2`` is enabled.

        Raises
        ------
        ValueError
            If no query annotations are available or any requested query column
            is absent from ``ldscore_result.query_columns``.
        """
        if include_model_categories is not None:
            include_full_partitioned_h2 = include_model_categories

        query_columns = _validate_partitioned_query_columns(ldscore_result, annotation_bundle.query_columns)
        rows = []
        per_query_category_tables: dict[str, pd.DataFrame] = {}
        per_query_metadata: dict[str, dict[str, object]] = {}
        for query_column in query_columns:
            dataset = self.build_dataset(sumstats_table, ldscore_result, config=config, query_columns=[query_column])
            hsq = self.estimate_h2(dataset, config=config)
            query_row = summarize_partitioned_h2(hsq, dataset, [query_column])
            rows.append(query_row)
            if include_full_partitioned_h2:
                per_query_category_tables[query_column] = summarize_partitioned_h2(
                    hsq,
                    dataset,
                    dataset.retained_ld_columns,
                    include_full_columns=True,
                )
                per_query_metadata[query_column] = {
                    "dropped_zero_variance_ld_columns": list(dataset.dropped_zero_variance_ld_columns),
                    "retained_ld_columns": list(dataset.retained_ld_columns),
                    "n_snps": len(dataset.merged),
                }
        if not rows:
            summary = pd.DataFrame(columns=PARTITIONED_H2_AGGREGATE_COLUMNS)
        else:
            summary = pd.concat(rows, axis=0, ignore_index=True)
        if include_full_partitioned_h2:
            return PartitionedH2BatchResult(
                summary=summary,
                per_query_category_tables=per_query_category_tables,
                per_query_metadata=per_query_metadata,
            )
        return summary

    def estimate_rg(
        self,
        sumstats_table_1: SumstatsTable,
        sumstats_table_2: SumstatsTable,
        ldscore_result: LDScoreResult,
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
        return self._fit_rg_dataset(dataset, config=config)

    def _fit_rg_dataset(
        self,
        dataset: RGRegressionDataset,
        config: RegressionConfig | None = None,
    ):
        """Fit one prepared genetic-correlation dataset with the legacy kernel."""
        config = config or self.regression_config
        merged = dataset.merged
        n_snp = len(merged)
        n_blocks = min(n_snp, config.n_blocks)
        intercept_hsq = _select_intercept(config.intercept_h2, use_intercept=config.use_intercept, default_when_disabled=1)
        intercept_gencov = _select_intercept(
            config.intercept_gencov,
            use_intercept=config.use_intercept,
            default_when_disabled=0,
        )
        return reg.RG(
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
            twostep=config.two_step_cutoff,
        )

    def estimate_rg_pairs(
        self,
        sumstats_tables: Sequence[SumstatsTable],
        ldscore_result: LDScoreResult,
        *,
        anchor_index: int | None = None,
        config: RegressionConfig | None = None,
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
        ldscore_result : LDScoreResult
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
            raise ValueError("--sumstats-sources must resolve to at least two sumstats files.")
        if anchor_index is not None and not 0 <= anchor_index < len(tables):
            raise ValueError(f"anchor_index must be in [0, {len(tables) - 1}], got {anchor_index}.")

        h2_rows = []
        for table in tables:
            h2_dataset = self.build_dataset(table, ldscore_result, config=config)
            hsq = self.estimate_h2(h2_dataset, config=config)
            h2_rows.append(summarize_total_h2(hsq, h2_dataset, trait_name=_trait_label(table)))
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
                fitted = self._fit_rg_dataset(dataset, config=config)
                full_row = _summarize_rg_pair(fitted, dataset, trait_1=trait_1, trait_2=trait_2, pair_kind=pair_kind)
                metadata = _rg_pair_metadata(tables[i], tables[j], dataset, full_row, config, pair_kind)
            except Exception as exc:
                error = _format_exception(exc)
                LOGGER.warning("Failed rg for pair '%s' vs '%s': %s", trait_1, trait_2, error, exc_info=True)
                full_row = _failed_rg_full_row(trait_1=trait_1, trait_2=trait_2, pair_kind=pair_kind, error=error)
                metadata = _failed_rg_pair_metadata(tables[i], tables[j], full_row, config, pair_kind)
            rg_full_rows.append(full_row)
            rg_rows.append(_concise_rg_row(full_row))
            per_pair_metadata.append(metadata)

        rg_full = pd.DataFrame(rg_full_rows, columns=RG_FULL_COLUMNS)
        rg = pd.DataFrame(rg_rows, columns=RG_CONCISE_COLUMNS)
        _apply_rg_multiple_testing(rg, rg_full)
        return RgResultFamily(rg=rg, rg_full=rg_full, h2_per_trait=h2_per_trait, per_pair_metadata=per_pair_metadata)


def _effective_snp_identifier_mode(
    sumstats_table: SumstatsTable,
    ldscore_result: LDScoreResult,
    runner_config: GlobalConfig,
) -> str:
    """Resolve the SNP identifier mode for a regression merge."""
    for snapshot in (ldscore_result.config_snapshot, sumstats_table.config_snapshot, runner_config):
        value = getattr(snapshot, "snp_identifier", None)
        if value is not None:
            return normalize_snp_identifier_mode(value)
    return "rsid"


def _coordinate_missing_mask(series: pd.Series) -> pd.Series:
    """Return rows with missing or NA-like coordinate tokens."""
    tokens = series.astype("string")
    return tokens.isna() | tokens.str.strip().str.lower().isin({"", "na", "nan", "none"})


def _with_chr_pos_key(frame: pd.DataFrame, *, context: str) -> pd.DataFrame:
    """Return a copy with a private canonical ``CHR:POS`` merge key."""
    chr_col, pos_col = infer_chr_pos_columns(frame.columns, context=context)
    chr_missing = _coordinate_missing_mask(frame[chr_col])
    pos_missing = _coordinate_missing_mask(frame[pos_col])
    pos_numeric = pd.to_numeric(frame[pos_col], errors="coerce")
    invalid_pos = (~pos_missing) & pos_numeric.isna()
    if invalid_pos.any():
        bad_value = frame.loc[invalid_pos, pos_col].iloc[0]
        raise ValueError(f"POS values in {context} must be numeric; got {bad_value!r}.")
    complete = ~(chr_missing | pos_missing)
    keyed = frame.loc[complete].copy()
    if keyed.empty:
        keyed[CHR_POS_KEY_COLUMN] = pd.Series(dtype="string")
        return keyed
    keyed[pos_col] = pos_numeric.loc[complete]
    non_integral = (keyed[pos_col] % 1) != 0
    if non_integral.any():
        bad_value = keyed.loc[non_integral, pos_col].iloc[0]
        raise ValueError(f"POS values in {context} must be integer base-pair positions; got {bad_value!r}.")
    keyed[pos_col] = keyed[pos_col].astype("int64")
    keyed[CHR_POS_KEY_COLUMN] = build_snp_id_series(keyed, "chr_pos")
    return keyed


def _select_count_key(count_totals: dict[str, np.ndarray], use_common_counts: bool) -> str:
    """Pick the regression count vector key, preferring common-SNP counts by default."""
    if use_common_counts and COMMON_COUNT_KEY in count_totals:
        return COMMON_COUNT_KEY
    if ALL_COUNT_KEY in count_totals:
        return ALL_COUNT_KEY
    return sorted(count_totals.keys())[0]


def _validate_partitioned_query_columns(ldscore_result: LDScoreResult, query_columns: Sequence[str]) -> list[str]:
    """Return requested query columns after enforcing the partitioned-h2 contract."""
    requested = list(query_columns)
    if not ldscore_result.query_columns or not requested:
        raise ValueError(PARTITIONED_H2_REQUIRES_QUERY_ANNOTATIONS_MESSAGE)
    available = list(ldscore_result.query_columns)
    missing = [column for column in requested if column not in available]
    if missing:
        raise ValueError(
            f"Unknown query annotation requested for partitioned-h2: {missing}. "
            f"Available query annotations: {available}"
        )
    return requested


def _assemble_regression_ldscore_table(ldscore_result: LDScoreResult, query_columns: Sequence[str]) -> pd.DataFrame:
    """Build the table used for one regression dataset from split LD-score tables."""
    baseline_table = ldscore_result.baseline_table.reset_index(drop=True)
    if not query_columns:
        return baseline_table.copy()
    if ldscore_result.query_table is None:
        raise ValueError("LD-score result does not contain query annotations.")
    missing = [column for column in query_columns if column not in ldscore_result.query_columns]
    if missing:
        raise ValueError(f"Unknown query LD-score columns requested: {missing}")
    query_table = ldscore_result.query_table.reset_index(drop=True)
    assert_same_snp_rows(
        baseline_table,
        query_table,
        context="query rows must match baseline rows on CHR/SNP/POS",
    )
    return pd.concat([baseline_table, query_table.loc[:, list(query_columns)]], axis=1)


def _count_totals_for_columns(count_records: Sequence[dict[str, Any]], columns: Sequence[str]) -> dict[str, np.ndarray]:
    """Return LDSC count vectors aligned to ``columns`` using manifest records."""
    records_by_column = {str(record["column"]): record for record in count_records}
    missing = [column for column in columns if column not in records_by_column]
    if missing:
        raise ValueError(f"LD-score manifest is missing count records for columns: {missing}")
    all_counts = [float(records_by_column[column]["all_reference_snp_count"]) for column in columns]
    count_totals = {ALL_COUNT_KEY: np.asarray(all_counts, dtype=np.float64)}
    if all("common_reference_snp_count" in records_by_column[column] for column in columns):
        count_totals[COMMON_COUNT_KEY] = np.asarray(
            [float(records_by_column[column]["common_reference_snp_count"]) for column in columns],
            dtype=np.float64,
        )
    return count_totals


def summarize_total_h2(hsq, dataset: RegressionDataset, trait_name: str | None = None) -> pd.DataFrame:
    """Build the one-row total-heritability summary from a fitted ``Hsq`` result."""
    return pd.DataFrame(
        [
            {
                "trait_name": trait_name,
                "n_snps": len(dataset.merged),
                "total_h2": _scalar(hsq.tot),
                "total_h2_se": _scalar(hsq.tot_se),
                "intercept": _scalar_or_value(hsq.intercept),
                "intercept_se": getattr(hsq, "intercept_se", None),
                "mean_chisq": _scalar_or_value(hsq.mean_chisq),
                "lambda_gc": _scalar_or_value(hsq.lambda_gc),
                "ratio": getattr(hsq, "ratio", None),
                "ratio_se": getattr(hsq, "ratio_se", None),
            }
        ]
    )


def summarize_partitioned_h2(
    hsq,
    dataset: RegressionDataset,
    annotation_columns: Sequence[str],
    *,
    include_full_columns: bool = False,
) -> pd.DataFrame:
    """Build public partitioned-h2 rows from a fitted ``Hsq`` result.

    ``annotation_columns`` must be a subset of ``dataset.retained_ld_columns``.
    Returned rows use legacy-style public column names. By default the result is
    the compact query-summary schema; ``include_full_columns=True`` also emits
    category h2 estimates and standard errors for full per-query model tables.
    """
    rows = []
    coefficients = np.ravel(hsq.coef)
    coefficient_covariance = _coefficient_covariance_matrix(hsq, coefficients)
    coefficient_ses = np.ravel(hsq.coef_se)
    category_h2 = np.ravel(hsq.cat)
    category_h2_ses = np.ravel(hsq.cat_se)
    proportions = np.ravel(hsq.prop)
    proportion_ses = np.ravel(hsq.prop_se)
    enrichments = np.ravel(hsq.enrichment) if getattr(hsq, "enrichment", None) is not None else None
    snp_counts = _retained_snp_counts(dataset)
    snp_proportions = _safe_vector_divide(snp_counts, float(np.sum(snp_counts)))
    n_blocks = getattr(hsq, "n_blocks", None)
    for annotation_column in annotation_columns:
        annotation_index = dataset.retained_ld_columns.index(annotation_column)
        coefficient = float(coefficients[annotation_index])
        coefficient_se = float(coefficient_ses[annotation_index])
        proportion_snps = float(snp_proportions[annotation_index])
        proportion_h2_se = float(proportion_ses[annotation_index])
        rows.append(
            {
                "Category": annotation_column,
                "Prop._SNPs": proportion_snps,
                "Category_h2": float(category_h2[annotation_index]),
                "Category_h2_std_error": float(category_h2_ses[annotation_index]),
                "Prop._h2": float(proportions[annotation_index]),
                "Prop._h2_std_error": proportion_h2_se,
                "Enrichment": math.nan if enrichments is None else float(enrichments[annotation_index]),
                "Enrichment_std_error": _safe_divide(proportion_h2_se, proportion_snps),
                "Enrichment_p": _enrichment_p_value(
                    annotation_index,
                    coefficients,
                    coefficient_covariance,
                    snp_counts,
                    n_blocks,
                ),
                "Coefficient": coefficient,
                "Coefficient_std_error": coefficient_se,
                "Coefficient_p": _coefficient_p_value(coefficient, coefficient_se),
            }
        )
    columns = PARTITIONED_H2_FULL_COLUMNS if include_full_columns else PARTITIONED_H2_AGGREGATE_COLUMNS
    return pd.DataFrame(rows, columns=columns)


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
) -> dict[str, object]:
    """Build the full public rg row from one fitted kernel result."""
    return {
        "trait_1": trait_1,
        "trait_2": trait_2,
        "n_snps_used": int(len(dataset.merged)),
        "rg": _required_numeric_scalar(getattr(rg_result, "rg_ratio", None), "rg"),
        "rg_se": _required_numeric_scalar(getattr(rg_result, "rg_se", None), "rg_se"),
        "z": _required_numeric_scalar(getattr(rg_result, "z", None), "z"),
        "p": _required_numeric_scalar(getattr(rg_result, "p", None), "p"),
        "p_fdr_bh": math.nan,
        "p_bonferroni": math.nan,
        "h2_1": _numeric_attr(getattr(rg_result, "hsq1", None), "tot", "h2_1"),
        "h2_1_se": _numeric_attr(getattr(rg_result, "hsq1", None), "tot_se", "h2_1_se"),
        "h2_2": _numeric_attr(getattr(rg_result, "hsq2", None), "tot", "h2_2"),
        "h2_2_se": _numeric_attr(getattr(rg_result, "hsq2", None), "tot_se", "h2_2_se"),
        "gencov": _numeric_attr(getattr(rg_result, "gencov", None), "tot", "gencov"),
        "gencov_se": _numeric_attr(getattr(rg_result, "gencov", None), "tot_se", "gencov_se"),
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
        "p_fdr_bh": full_row["p_fdr_bh"],
        "note": "" if status == "ok" else FAILED_RG_NOTE,
    }


def _apply_rg_multiple_testing(rg: pd.DataFrame, rg_full: pd.DataFrame) -> None:
    """Apply local BH-FDR and Bonferroni adjustments, ignoring NaN p-values."""
    p_values = pd.to_numeric(rg_full.get("p", pd.Series(dtype=float)), errors="coerce").to_numpy(dtype=float)
    fdr = _benjamini_hochberg(p_values)
    bonferroni = _bonferroni(p_values)
    for frame in (rg, rg_full):
        if "p_fdr_bh" not in frame.columns:
            frame["p_fdr_bh"] = math.nan
        frame.loc[:, "p_fdr_bh"] = fdr
    if "p_bonferroni" not in rg_full.columns:
        rg_full["p_bonferroni"] = math.nan
    rg_full.loc[:, "p_bonferroni"] = bonferroni


def _benjamini_hochberg(p_values: np.ndarray) -> np.ndarray:
    """Return BH-adjusted p-values with NaNs preserved."""
    p = np.asarray(p_values, dtype=np.float64)
    adjusted = np.full(p.shape, math.nan, dtype=np.float64)
    finite = np.isfinite(p)
    if not finite.any():
        return adjusted
    finite_indices = np.flatnonzero(finite)
    finite_p = p[finite]
    order = np.argsort(finite_p)
    ranked = finite_p[order] * len(finite_p) / np.arange(1, len(finite_p) + 1, dtype=np.float64)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    adjusted[finite_indices[order]] = np.clip(ranked, 0.0, 1.0)
    return adjusted


def _bonferroni(p_values: np.ndarray) -> np.ndarray:
    """Return Bonferroni-adjusted p-values with NaNs preserved."""
    p = np.asarray(p_values, dtype=np.float64)
    adjusted = np.full(p.shape, math.nan, dtype=np.float64)
    finite = np.isfinite(p)
    if not finite.any():
        return adjusted
    adjusted[finite] = np.clip(p[finite] * int(finite.sum()), 0.0, 1.0)
    return adjusted


def _rg_pair_metadata(
    table_1: SumstatsTable,
    table_2: SumstatsTable,
    dataset: RGRegressionDataset,
    full_row: dict[str, object],
    config: RegressionConfig,
    pair_kind: str,
) -> dict[str, object]:
    """Build per-pair metadata for an estimated rg pair."""
    n_snps = int(len(dataset.merged))
    return {
        "trait_1": full_row["trait_1"],
        "trait_2": full_row["trait_2"],
        "source_1": table_1.source_path,
        "source_2": table_2.source_path,
        "pair_kind": pair_kind,
        "status": full_row["status"],
        "error": full_row["error"],
        "n_snps_used": n_snps,
        "n_blocks_used": min(n_snps, int(config.n_blocks)),
        "count_key_used_for_regression": dataset.count_key_used_for_regression,
        "retained_ld_columns": list(dataset.retained_ld_columns),
        "dropped_zero_variance_ld_columns": list(dataset.dropped_zero_variance_ld_columns),
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
                raise ValueError(f"{field} is non-numeric kernel value 'NA'.")
            return math.nan
        try:
            return float(stripped)
        except ValueError as exc:
            if fail_on_non_numeric:
                raise ValueError(f"{field} is non-numeric kernel value {value!r}.") from exc
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
            raise ValueError(f"{field} is non-numeric kernel value {value!r}.") from exc
        return math.nan


def _format_exception(exc: Exception) -> str:
    """Return compact user-facing error text for a pair failure."""
    message = str(exc).strip()
    return f"{exc.__class__.__name__}: {message}" if message else exc.__class__.__name__


def _retained_snp_counts(dataset: RegressionDataset) -> np.ndarray:
    """Return reference SNP counts aligned to retained LD-score columns."""
    counts = np.asarray(
        dataset.reference_snp_count_totals[dataset.count_key_used_for_regression],
        dtype=np.float64,
    ).reshape(-1)
    if len(counts) != len(dataset.retained_ld_columns):
        raise ValueError(
            "Reference SNP count vector is not aligned to retained LD-score columns."
        )
    return counts


def _coefficient_covariance_matrix(hsq, coefficients: np.ndarray) -> np.ndarray:
    """Return a coefficient covariance matrix, falling back to diagonal SEs."""
    covariance = getattr(hsq, "coef_cov", None)
    if covariance is not None:
        covariance = np.asarray(covariance, dtype=np.float64)
        if covariance.shape == (len(coefficients), len(coefficients)):
            return covariance
    coefficient_ses = np.ravel(getattr(hsq, "coef_se", np.full(len(coefficients), math.nan)))
    return np.diag(np.square(coefficient_ses))


def _coefficient_p_value(coefficient: float, coefficient_se: float) -> float:
    """Return a two-sided normal p-value for a coefficient estimate."""
    if not np.isfinite(coefficient) or not np.isfinite(coefficient_se) or coefficient_se == 0:
        return math.nan
    return float(2 * stats.norm.sf(abs(coefficient / coefficient_se)))


def _enrichment_p_value(
    annotation_index: int,
    coefficients: np.ndarray,
    coefficient_covariance: np.ndarray,
    snp_counts: np.ndarray,
    n_blocks,
) -> float:
    """Return the legacy-style category-vs-complement enrichment p-value.

    The contrast compares the selected coefficient with the SNP-count-weighted
    mean coefficient of all remaining retained categories. A missing complement,
    unavailable jackknife block count, or non-positive variance yields NaN.
    """
    total_count = float(np.sum(snp_counts))
    category_count = float(snp_counts[annotation_index])
    complement_count = total_count - category_count
    if (
        total_count <= 0
        or category_count <= 0
        or complement_count <= 0
        or n_blocks is None
        or int(n_blocks) <= 0
    ):
        return math.nan
    contrast = -snp_counts / complement_count
    contrast[annotation_index] = 1.0
    contrast_estimate = float(np.dot(contrast, coefficients))
    contrast_variance = float(np.dot(np.dot(contrast, coefficient_covariance), contrast.T))
    if not np.isfinite(contrast_variance) or contrast_variance <= 0:
        return math.nan
    contrast_se = math.sqrt(contrast_variance)
    return float(2 * stats.t.sf(abs(contrast_estimate / contrast_se), int(n_blocks)))


def _safe_divide(numerator: float, denominator: float) -> float:
    """Return ``numerator / denominator`` or NaN for invalid ratios."""
    if not np.isfinite(numerator) or not np.isfinite(denominator) or denominator == 0:
        return math.nan
    return float(numerator / denominator)


def _safe_vector_divide(numerator: np.ndarray, denominator: float) -> np.ndarray:
    """Vectorized safe division for scalar denominators."""
    if not np.isfinite(denominator) or denominator == 0:
        return np.full_like(numerator, math.nan, dtype=np.float64)
    return numerator / denominator


def _scalar(value) -> float:
    """Return the first scalar numeric value from a possibly array-like object."""
    return float(np.ravel(value)[0])


def _scalar_or_value(value):
    """Return a scalar float for array-like values, otherwise preserve the original value."""
    if hasattr(value, "__array__") or isinstance(value, (list, tuple)):
        return _scalar(value)
    return value


def _select_intercept(value: float | None, use_intercept: bool, default_when_disabled: float):
    """Resolve fixed-intercept config into the scalar expected by the kernel."""
    if not use_intercept:
        return default_when_disabled
    if value is None:
        return None
    return float(value)


def add_h2_arguments(parser) -> None:
    """Register heritability CLI arguments on ``parser``."""
    _add_common_regression_arguments(parser, include_h2_intercept=True)
    parser.add_argument("--sumstats-file", required=True, help="Munged .sumstats(.gz) file.")
    parser.add_argument("--trait-name", default=None, help="Optional trait label for summaries.")


def add_partitioned_h2_arguments(parser) -> None:
    """Register partitioned-heritability CLI arguments on ``parser``."""
    _add_common_regression_arguments(parser, include_h2_intercept=True)
    parser.add_argument("--sumstats-file", required=True, help="Munged .sumstats(.gz) file.")
    parser.add_argument("--trait-name", default=None, help="Optional trait label for summaries.")
    parser.add_argument(
        "--write-per-query-results",
        action="store_true",
        default=False,
        help="Also write one sanitized result folder per query annotation under output_dir/query_annotations.",
    )


def add_rg_arguments(parser) -> None:
    """Register genetic-correlation CLI arguments on ``parser``."""
    _add_common_regression_arguments(parser, include_h2_intercept=False)
    parser.add_argument(
        "--sumstats-sources",
        nargs="+",
        required=True,
        help=(
            "Munged .sumstats(.gz) files or glob patterns. At least two resolved files are required; "
            "use --output-dir for the complete rg output family."
        ),
    )
    parser.add_argument(
        "--anchor-trait-file",
        default=None,
        help="Optional anchor sumstats file or trait name; compute only anchor-vs-rest correlations.",
    )
    parser.add_argument(
        "--write-per-pair-detail",
        action="store_true",
        default=False,
        help="Also write one sanitized result folder per tested trait pair under output_dir/pairs.",
    )
    parser.add_argument("--intercept-h2", type=float, default=None, help="Fixed h2 intercept broadcast to every rg pair.")
    parser.add_argument("--intercept-gencov", type=float, default=None, help="Fixed genetic-covariance intercept for every rg pair.")


def run_h2_from_args(args):
    """Run single-trait heritability estimation from parsed CLI arguments.

    When ``args.output_dir`` is provided, the workflow preflights ``h2.tsv`` and
    ``h2.log`` before loading inputs. Without an output directory, it returns
    the summary table without creating a log file.
    """
    output_dir, log_path = _preflight_regression_outputs(args, "h2", ["h2.tsv"])
    with workflow_logging("h2", log_path, log_level=getattr(args, "log_level", "INFO")):
        runner, config = _runner_from_args(args)
        print_global_config_banner("run_h2_from_args", runner.global_config)
        log_inputs(sumstats_file=args.sumstats_file, ldscore_dir=args.ldscore_dir, output_dir=output_dir or "none")
        LOGGER.info(f"Starting h2 regression for '{args.sumstats_file}' using LD-score directory '{args.ldscore_dir}'.")
        sumstats_table = _load_sumstats_table(args.sumstats_file, getattr(args, "trait_name", None))
        ldscore_result = load_ldscore_from_dir(args.ldscore_dir)
        with suppress_global_config_banner():
            dataset = runner.build_dataset(sumstats_table, ldscore_result, config=config)
        hsq = runner.estimate_h2(dataset, config=config)
        summary = summarize_total_h2(hsq, dataset, trait_name=sumstats_table.trait_name)
        _maybe_write_dataframe(
            summary,
            getattr(args, "output_dir", None),
            "h2.tsv",
            overwrite=getattr(args, "overwrite", False),
        )
        if output_dir is not None:
            log_outputs(summary=str(Path(output_dir) / "h2.tsv"))
        LOGGER.info(f"Finished h2 regression with {len(dataset.merged)} regression SNPs.")
    return summary


def run_partitioned_h2_from_args(args):
    """Run batch partitioned heritability from parsed CLI arguments.

    When ``args.output_dir`` is provided, the workflow preflights
    ``partitioned_h2.tsv``, optional ``query_annotations/``, and
    ``partitioned-h2.log`` before loading inputs. With overwrite enabled,
    successful aggregate-only runs remove stale ``query_annotations/`` trees.
    Without an output directory, it returns the summary table without creating
    a log file.
    """
    preflight_names = ["partitioned_h2.tsv"]
    if getattr(args, "write_per_query_results", False):
        preflight_names.append("query_annotations")
    output_dir, log_path = _preflight_regression_outputs(
        args,
        "partitioned-h2",
        preflight_names,
        owned_output_names=["partitioned_h2.tsv", "query_annotations"],
    )
    with workflow_logging("partitioned-h2", log_path, log_level=getattr(args, "log_level", "INFO")):
        runner, config = _runner_from_args(args)
        print_global_config_banner("run_partitioned_h2_from_args", runner.global_config)
        log_inputs(sumstats_file=args.sumstats_file, ldscore_dir=args.ldscore_dir, output_dir=output_dir or "none")
        LOGGER.info(
            f"Starting partitioned-h2 regression for '{args.sumstats_file}' using LD-score directory '{args.ldscore_dir}'."
        )
        sumstats_table = _load_sumstats_table(args.sumstats_file, getattr(args, "trait_name", None))
        ldscore_result = load_ldscore_from_dir(args.ldscore_dir)
        query_bundle = SimpleNamespace(
            query_columns=_validate_partitioned_query_columns(ldscore_result, ldscore_result.query_columns)
        )
        with suppress_global_config_banner():
            result = runner.estimate_partitioned_h2_batch(
                sumstats_table,
                ldscore_result,
                query_bundle,
                config=config,
                include_full_partitioned_h2=getattr(args, "write_per_query_results", False),
            )
        if isinstance(result, PartitionedH2BatchResult):
            summary = result.summary
            per_query_category_tables = result.per_query_category_tables
            per_query_metadata = result.per_query_metadata
        else:
            summary = result
            per_query_category_tables = None
            per_query_metadata = None
        output_dir_arg = getattr(args, "output_dir", None)
        if output_dir_arg:
            written = PartitionedH2DirectoryWriter().write(
                summary,
                PartitionedH2OutputConfig(
                    output_dir=output_dir_arg,
                    overwrite=getattr(args, "overwrite", False),
                    write_per_query_results=getattr(args, "write_per_query_results", False),
                ),
                per_query_category_tables=per_query_category_tables,
                metadata={
                    "trait_name": sumstats_table.trait_name,
                    "count_kind": getattr(args, "count_kind", "common"),
                    "ldscore_dir": getattr(args, "ldscore_dir", None),
                },
                per_query_metadata=per_query_metadata,
            )
            log_outputs(**written)
        LOGGER.info(
            f"Finished partitioned-h2 regression for {len(ldscore_result.query_columns)} query annotations "
            f"and {len(summary)} summary rows."
        )
    return summary


def run_rg_from_args(args):
    """Run multi-trait genetic-correlation estimation from parsed CLI args.

    With ``--output-dir``, the workflow writes the rg result family:
    ``rg.tsv``, ``rg_full.tsv``, ``h2_per_trait.tsv``, optional ``pairs/``,
    and workflow-owned ``rg.log``. Without ``--output-dir``, it returns the same
    in-memory result family and the CLI prints only ``rg.tsv`` to stdout.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed rg options. Required fields are ``sumstats_sources`` and
        ``ldscore_dir``. Optional fields include ``anchor_trait_file``,
        ``output_dir``, ``write_per_pair_detail``, intercept settings, and
        common regression options.

    Returns
    -------
    RgResultFamily
        Complete in-memory result family. The function itself never prints; CLI
        dispatch owns stdout behavior for no-output-dir runs.

    Raises
    ------
    ValueError
        If fewer than two sumstats paths resolve, per-pair detail is requested
        without an output directory, or fixed-intercept flags conflict with
        ``--no-intercept``.
    """
    _validate_intercept_conflicts(args)
    if getattr(args, "write_per_pair_detail", False) and not getattr(args, "output_dir", None):
        raise ValueError("--write-per-pair-detail requires --output-dir.")
    output_names = ["rg.tsv", "rg_full.tsv", "h2_per_trait.tsv"]
    if getattr(args, "write_per_pair_detail", False):
        output_names.append("pairs")
    output_dir, log_path = _preflight_regression_outputs(
        args,
        "rg",
        output_names,
        owned_output_names=["rg.tsv", "rg_full.tsv", "h2_per_trait.tsv", "pairs"],
    )
    sumstats_paths = resolve_file_group(getattr(args, "sumstats_sources", ()), label="sumstats sources")
    if len(sumstats_paths) < 2:
        raise ValueError("--sumstats-sources must resolve to at least two sumstats files.")
    with workflow_logging("rg", log_path, log_level=getattr(args, "log_level", "INFO")):
        runner, config = _runner_from_args(args)
        print_global_config_banner("run_rg_from_args", runner.global_config)
        log_inputs(
            sumstats_sources=[str(path) for path in sumstats_paths],
            anchor_trait_file=getattr(args, "anchor_trait_file", None) or "none",
            ldscore_dir=args.ldscore_dir,
            output_dir=output_dir or "none",
        )
        LOGGER.info(
            f"Starting rg regression for {len(sumstats_paths)} sumstats files "
            f"using LD-score directory '{args.ldscore_dir}'."
        )
        sumstats_tables = [_load_sumstats_table(str(path), None) for path in sumstats_paths]
        sumstats_tables = _disambiguate_trait_names(sumstats_tables)
        anchor_index = _resolve_anchor_index(getattr(args, "anchor_trait_file", None), sumstats_paths, sumstats_tables)
        ldscore_result = load_ldscore_from_dir(args.ldscore_dir)
        with suppress_global_config_banner():
            result = runner.estimate_rg_pairs(
                sumstats_tables,
                ldscore_result,
                anchor_index=anchor_index,
                config=config,
            )
        output_dir_arg = getattr(args, "output_dir", None)
        if output_dir_arg:
            written = RgDirectoryWriter().write(
                result,
                RgOutputConfig(
                    output_dir=output_dir_arg,
                    overwrite=getattr(args, "overwrite", False),
                    write_per_pair_detail=getattr(args, "write_per_pair_detail", False),
                ),
            )
            log_outputs(**written)
        LOGGER.info(
            f"Finished rg regression for {len(result.rg)} trait pairs "
            f"and {len(result.h2_per_trait)} traits."
        )
    return result


def _add_common_regression_arguments(parser, include_h2_intercept: bool) -> None:
    """Add the file and model options shared by all regression subcommands."""
    parser.allow_abbrev = False
    parser.add_argument("--ldscore-dir", required=True, help="Canonical LD-score result directory written by `ldsc ldscore`.")
    parser.add_argument(
        "--count-kind",
        choices=("common", "all"),
        default="common",
        help="Reference SNP count vector used by regression.",
    )
    parser.add_argument("--output-dir", default=None, help="Optional output directory for workflow result files; strongly recommended.")
    parser.add_argument("--overwrite", action="store_true", default=False, help="Replace existing workflow output artifacts.")
    parser.add_argument("--n-blocks", type=int, default=200)
    parser.add_argument("--no-intercept", action="store_true", default=False, help="Fix the intercept to the LDSC default.")
    if include_h2_intercept:
        parser.add_argument("--intercept-h2", type=float, default=None, help="Fixed h2 intercept for single-trait runs.")
    parser.add_argument("--two-step-cutoff", type=float, default=None)
    parser.add_argument("--chisq-max", type=float, default=None)
    parser.add_argument("--log-level", default="INFO", choices=("DEBUG", "INFO", "WARNING", "ERROR"), help="Logging verbosity.")


def _preflight_regression_outputs(
    args,
    workflow_name: str,
    output_names: list[str],
    owned_output_names: list[str] | None = None,
) -> tuple[str | None, Path | None]:
    """Preflight regression outputs and return normalized output dir plus log path."""
    output_dir_arg = getattr(args, "output_dir", None)
    if not output_dir_arg:
        return None, None
    output_dir = ensure_output_directory(output_dir_arg, label="output directory")
    paths = [output_dir / name for name in output_names]
    log_path = output_dir / f"{workflow_name}.log"
    if owned_output_names is None:
        ensure_output_paths_available(
            [*paths, log_path],
            overwrite=getattr(args, "overwrite", False),
            label="regression output artifact",
        )
    else:
        preflight_output_artifact_family(
            [*paths, log_path],
            [*(output_dir / name for name in owned_output_names), log_path],
            overwrite=getattr(args, "overwrite", False),
            label="regression output artifact",
        )
    return str(output_dir), log_path


def _validate_intercept_conflicts(args) -> None:
    """Reject contradictory intercept options before inputs are loaded."""
    if not getattr(args, "no_intercept", False):
        return
    conflicting = []
    if getattr(args, "intercept_h2", None) is not None:
        conflicting.append("--intercept-h2")
    if getattr(args, "intercept_gencov", None) is not None:
        conflicting.append("--intercept-gencov")
    if conflicting:
        raise ValueError(f"--no-intercept cannot be combined with {', '.join(conflicting)}.")


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
            LOGGER.info("Disambiguated duplicate rg trait name '%s' as '%s'.", original, resolved)
    return [replace(table, trait_name=label) for table, label in zip(tables, unique)]


def _resolve_anchor_index(
    anchor_trait_file: str | None,
    sumstats_paths: Sequence[str],
    sumstats_tables: Sequence[SumstatsTable],
) -> int | None:
    """Resolve ``--anchor-trait-file`` against source paths or trait names."""
    if not anchor_trait_file:
        return None
    token = normalize_path_token(anchor_trait_file)
    matches: set[int] = set()
    try:
        anchor_path = Path(token).resolve(strict=False)
        for idx, path in enumerate(sumstats_paths):
            if Path(path).resolve(strict=False) == anchor_path:
                matches.add(idx)
    except OSError:
        pass
    for idx, table in enumerate(sumstats_tables):
        if table.trait_name == token:
            matches.add(idx)
    if len(matches) != 1:
        available = [table.trait_name or Path(path).name for table, path in zip(sumstats_tables, sumstats_paths)]
        raise ValueError(
            f"--anchor-trait-file must match exactly one input path or trait name; "
            f"got {len(matches)} matches for {anchor_trait_file!r}. Available traits: {available}"
        )
    return next(iter(matches))


def _runner_from_args(args) -> tuple[RegressionRunner, RegressionConfig]:
    """Build the regression workflow objects from parsed CLI arguments."""
    _validate_intercept_conflicts(args)
    count_kind = getattr(args, "count_kind", "common")
    config = RegressionConfig(
        n_blocks=args.n_blocks,
        use_common_counts=(count_kind == "common"),
        use_intercept=not args.no_intercept,
        intercept_h2=args.intercept_h2,
        intercept_gencov=getattr(args, "intercept_gencov", None),
        two_step_cutoff=args.two_step_cutoff,
        chisq_max=args.chisq_max,
    )
    runner = RegressionRunner(get_global_config(), config)
    return runner, config


def _load_sumstats_table(path: str, trait_name: str | None) -> SumstatsTable:
    """Load one curated sumstats artifact through the public workflow helper."""
    return load_sumstats(path, trait_name=trait_name)


def load_ldscore_from_dir(
    ldscore_dir: str,
    snp_identifier: str | None = None,
) -> LDScoreResult:
    """Load a canonical LD-score result directory.

    Parameters
    ----------
    ldscore_dir : str
        Directory containing ``manifest.json``, ``ldscore.baseline.parquet``,
        and optional ``ldscore.query.parquet`` files written by the public
        LD-score writer.
    snp_identifier : {"rsid", "chr_pos"} or None, optional
        Identifier mode used to reconstruct the public regression SNP set from
        the baseline table. When omitted, the manifest value is used, falling
        back to ``"rsid"`` for legacy directories.

    Returns
    -------
    LDScoreResult
        Disk-loaded LD-score result. Missing or invalid manifest
        ``config_snapshot`` provenance emits a warning and returns
        ``config_snapshot=None`` so legacy directories can still be used.

    Raises
    ------
    FileNotFoundError
        If ``manifest.json`` is absent.
    NotADirectoryError
        If ``ldscore_dir`` is not an existing directory.
    ValueError
        If the manifest format or required file records are unsupported.
    """
    root = Path(normalize_path_token(ldscore_dir))
    if not root.is_dir():
        raise NotADirectoryError(f"LD-score directory does not exist or is not a directory: {root}")
    manifest_path = root / "manifest.json"
    if not manifest_path.exists():
        raise FileNotFoundError(f"LD-score directory is missing manifest.json: {root}")
    LOGGER.info(f"Loading LD-score result directory from '{root}'.")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("format") != LDSCORE_RESULT_FORMAT:
        raise ValueError(f"Unsupported LD-score directory format: {manifest.get('format')!r}")
    files = manifest.get("files", {})
    baseline_rel = files.get("baseline")
    if not baseline_rel:
        raise ValueError("LD-score manifest is missing files.baseline.")
    baseline_table = pd.read_parquet(root / baseline_rel)
    query_rel = files.get("query")
    query_table = pd.read_parquet(root / query_rel) if query_rel else None
    baseline_columns = [str(column) for column in manifest.get("baseline_columns", [])]
    query_columns = [str(column) for column in manifest.get("query_columns", [])]
    count_records = [dict(record) for record in manifest.get("counts", [])]
    effective_identifier = snp_identifier or manifest.get("snp_identifier") or "rsid"
    config_snapshot = _global_config_from_manifest(manifest)
    result = LDScoreResult(
        baseline_table=baseline_table,
        query_table=query_table,
        count_records=count_records,
        baseline_columns=baseline_columns,
        query_columns=query_columns,
        ld_reference_snps=frozenset(),
        ld_regression_snps=frozenset(build_snp_id_series(baseline_table, effective_identifier)),
        chromosome_results=[],
        output_paths={
            "manifest": str(manifest_path),
            "baseline": str(root / baseline_rel),
            **({"query": str(root / query_rel)} if query_rel else {}),
        },
        count_config=dict(manifest.get("count_config", {})),
        config_snapshot=config_snapshot,
    )
    result.validate(require_query_alignment=False)
    query_rows = 0 if query_table is None else len(query_table)
    LOGGER.info(
        f"Loaded LD-score directory '{root}' with {len(baseline_table)} baseline rows, "
        f"{query_rows} query rows, and config provenance {'present' if config_snapshot is not None else 'unknown'}."
    )
    return result


def _global_config_from_manifest(manifest: dict[str, Any]) -> GlobalConfig | None:
    """Recreate a GlobalConfig snapshot when the manifest contains one."""
    snapshot = manifest.get("config_snapshot")
    if not isinstance(snapshot, dict):
        warnings.warn(
            "LD-score manifest GlobalConfig provenance is missing; treating config compatibility as unknown.",
            UserWarning,
            stacklevel=3,
        )
        return None
    try:
        _recovered_snp = snapshot.get("snp_identifier") or manifest.get("snp_identifier") or "rsid"
        return GlobalConfig(
            snp_identifier=_recovered_snp,
            genome_build=(
                snapshot.get("genome_build")
                or manifest.get("genome_build")
                or ("auto" if _recovered_snp == "chr_pos" else None)
            ),
            log_level=snapshot.get("log_level", "INFO"),
            fail_on_missing_metadata=bool(snapshot.get("fail_on_missing_metadata", False)),
        )
    except Exception as exc:
        LOGGER.debug(f"Invalid LD-score manifest GlobalConfig provenance: {exc}", exc_info=True)
        warnings.warn(
            f"LD-score manifest GlobalConfig provenance is invalid ({exc}); "
            "treating config compatibility as unknown.",
            UserWarning,
            stacklevel=3,
        )
        return None


def _maybe_write_dataframe(
    df: pd.DataFrame,
    output_dir: str | None,
    filename: str,
    overwrite: bool = False,
) -> None:
    """Write a fixed-name summary table only when an output directory is requested.

    Existing ``h2.tsv``, ``partitioned_h2.tsv``, or ``rg.tsv`` files are refused
    before writing unless ``overwrite`` is true.
    """
    if not output_dir:
        return
    path = ensure_output_directory(output_dir, label="output directory") / filename
    ensure_output_paths_available([path], overwrite=overwrite, label="regression output artifact")
    df.to_csv(path, sep="\t", index=False)
