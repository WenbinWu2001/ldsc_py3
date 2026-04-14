"""Workflow-layer regression adapters and dataset assembly."""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any

import numpy as np
import pandas as pd
from scipy import stats

from .config import CommonConfig, RegressionConfig
from . import regressions as reg
from . import sumstats as legacy_sumstats
from .ldscore_workflow import LDScoreResult
from .sumstats_munger import SumstatsTable


COMMON_COUNT_KEY = "common_reference_snp_counts_maf_gt_0_05"
ALL_COUNT_KEY = "all_reference_snp_counts"


@dataclass(frozen=True)
class RegressionDataset:
    merged: pd.DataFrame
    ref_ld_columns: list[str]
    weight_column: str
    reference_snp_count_totals: dict[str, np.ndarray]
    count_key_used_for_regression: str
    retained_ld_columns: list[str]
    dropped_zero_variance_ld_columns: list[str]
    trait_names: list[str]
    chromosomes_aggregated: list[str]

    def validate(self) -> None:
        required = {"SNP", self.weight_column, "Z", "N"}
        missing = required - set(self.merged.columns)
        if missing:
            raise ValueError(f"RegressionDataset is missing required columns: {sorted(missing)}")


class RegressionRunner:
    """Build regression datasets and dispatch legacy estimator kernels."""

    def __init__(
        self,
        common_config: CommonConfig | None = None,
        regression_config: RegressionConfig | None = None,
    ) -> None:
        self.common_config = common_config or CommonConfig()
        self.regression_config = regression_config or RegressionConfig()

    def build_dataset(
        self,
        sumstats_table: SumstatsTable,
        ldscore_result: LDScoreResult,
        config: RegressionConfig | None = None,
    ) -> RegressionDataset:
        config = config or self.regression_config
        ref_ld_frame = pd.concat(
            [
                ldscore_result.reference_metadata.loc[:, ["SNP"]].reset_index(drop=True),
                ldscore_result.ld_scores.reset_index(drop=True),
            ],
            axis=1,
        )
        weight_frame = pd.concat(
            [
                ldscore_result.regression_metadata.loc[:, ["SNP"]].reset_index(drop=True),
                ldscore_result.w_ld.reset_index(drop=True),
            ],
            axis=1,
        )
        weight_column = ldscore_result.w_ld.columns[0]
        merged = pd.merge(sumstats_table.data, ref_ld_frame, how="inner", on="SNP", sort=False)
        merged = pd.merge(merged, weight_frame, how="inner", on="SNP", sort=False)
        ref_ld_columns = list(ldscore_result.ld_scores.columns)

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

        count_key = _select_count_key(ldscore_result.snp_count_totals, config.use_m_5_50)
        count_totals = {
            key: np.asarray(values, dtype=np.float64)
            for key, values in ldscore_result.snp_count_totals.items()
        }
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
        )
        dataset.validate()
        return dataset

    def estimate_h2(
        self,
        dataset: RegressionDataset,
        config: RegressionConfig | None = None,
    ):
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
        annotation_bundle,
        config: RegressionConfig | None = None,
    ) -> pd.DataFrame:
        del annotation_bundle  # grouping metadata is carried by the selected columns
        dataset = self.build_dataset(sumstats_table, ldscore_result, config=config)
        hsq = self.estimate_h2(dataset, config=config)
        query_column = ldscore_result.query_columns[-1] if ldscore_result.query_columns else ldscore_result.baseline_columns[-1]
        query_index = dataset.retained_ld_columns.index(query_column)
        coefficient = float(hsq.coef[query_index])
        coefficient_se = float(hsq.coef_se[query_index])
        z_value = math.nan if coefficient_se == 0 else coefficient / coefficient_se
        p_value = math.nan if coefficient_se == 0 else 2 * stats.norm.sf(abs(z_value))
        enrichment = float(np.ravel(hsq.enrichment)[query_index]) if getattr(hsq, "enrichment", None) is not None else math.nan
        return pd.DataFrame(
            [
                {
                    "query_annotation": query_column,
                    "coefficient": coefficient,
                    "coefficient_se": coefficient_se,
                    "coefficient_z": z_value,
                    "coefficient_p": p_value,
                    "category_h2": float(np.ravel(hsq.cat)[query_index]),
                    "category_h2_se": float(np.ravel(hsq.cat_se)[query_index]),
                    "proportion_h2": float(np.ravel(hsq.prop)[query_index]),
                    "proportion_h2_se": float(np.ravel(hsq.prop_se)[query_index]),
                    "enrichment": enrichment,
                }
            ]
        )

    def estimate_partitioned_h2_batch(
        self,
        sumstats_table: SumstatsTable,
        ldscore_result: LDScoreResult,
        annotation_bundle,
        config: RegressionConfig | None = None,
    ) -> pd.DataFrame:
        rows = []
        for query_column in annotation_bundle.query_columns:
            subset_result = _subset_ldscore_result(ldscore_result, baseline_columns=ldscore_result.baseline_columns, query_columns=[query_column])
            query_row = self.estimate_partitioned_h2(sumstats_table, subset_result, annotation_bundle, config=config)
            rows.append(query_row)
        if not rows:
            return pd.DataFrame(columns=["query_annotation"])
        return pd.concat(rows, axis=0, ignore_index=True)

    def estimate_rg(
        self,
        sumstats_table_1: SumstatsTable,
        sumstats_table_2: SumstatsTable,
        ldscore_result: LDScoreResult,
        config: RegressionConfig | None = None,
    ):
        config = config or self.regression_config
        dataset_1 = self.build_dataset(sumstats_table_1, ldscore_result, config=config)
        merged = pd.merge(
            dataset_1.merged.rename(columns={"N": "N1", "Z": "Z1", "A1": "A1", "A2": "A2"}),
            sumstats_table_2.data.rename(columns={"A1": "A1x", "A2": "A2x", "N": "N2", "Z": "Z2"}),
            how="inner",
            on="SNP",
            sort=False,
        ).dropna(how="any")
        if {"A1", "A2", "A1x", "A2x"}.issubset(merged.columns):
            alleles = merged["A1"] + merged["A2"] + merged["A1x"] + merged["A2x"]
            keep = legacy_sumstats._filter_alleles(alleles)
            merged = merged.loc[keep].reset_index(drop=True)
            merged["Z2"] = legacy_sumstats._align_alleles(merged["Z2"], alleles.loc[keep])
        n_snp = len(merged)
        n_blocks = min(n_snp, config.n_blocks)
        intercept_hsq1 = _select_intercept(config.intercept_h2, index=0, use_intercept=config.use_intercept, default_when_disabled=1)
        intercept_hsq2 = _select_intercept(config.intercept_h2, index=1, use_intercept=config.use_intercept, default_when_disabled=1)
        intercept_gencov = _select_intercept(config.intercept_gencov, index=1, use_intercept=config.use_intercept, default_when_disabled=0)
        return reg.RG(
            np.asarray(merged[["Z1"]]),
            np.asarray(merged[["Z2"]]),
            np.asarray(merged[dataset_1.retained_ld_columns]),
            np.asarray(merged[[dataset_1.weight_column]]),
            np.asarray(merged[["N1"]]),
            np.asarray(merged[["N2"]]),
            np.asarray(dataset_1.reference_snp_count_totals[dataset_1.count_key_used_for_regression]).reshape((1, -1)),
            intercept_hsq1=intercept_hsq1,
            intercept_hsq2=intercept_hsq2,
            intercept_gencov=intercept_gencov,
            n_blocks=n_blocks,
            twostep=config.two_step_cutoff,
        )

    def run_legacy_from_args(self, args):
        if getattr(args, "rg", None):
            return legacy_sumstats.estimate_rg(args)
        if getattr(args, "h2", None):
            return legacy_sumstats.estimate_h2(args)
        if getattr(args, "h2_cts", None):
            return legacy_sumstats.cell_type_specific(args)
        raise ValueError("No regression analysis selected.")


def _select_count_key(count_totals: dict[str, np.ndarray], use_m_5_50: bool) -> str:
    if use_m_5_50 and COMMON_COUNT_KEY in count_totals:
        return COMMON_COUNT_KEY
    if ALL_COUNT_KEY in count_totals:
        return ALL_COUNT_KEY
    return sorted(count_totals.keys())[0]


def _select_intercept(value, index: int, use_intercept: bool, default_when_disabled: float):
    if not use_intercept:
        return default_when_disabled
    if value is None:
        return None
    if isinstance(value, list):
        return float(value[index])
    return float(value)


def _subset_ldscore_result(ldscore_result: LDScoreResult, baseline_columns: list[str], query_columns: list[str]) -> LDScoreResult:
    selected_columns = baseline_columns + query_columns
    selected_index = [ldscore_result.ld_scores.columns.get_loc(column) for column in selected_columns]
    snp_count_totals = {
        key: np.asarray(values)[selected_index]
        for key, values in ldscore_result.snp_count_totals.items()
    }
    return LDScoreResult(
        reference_metadata=ldscore_result.reference_metadata.copy(),
        ld_scores=ldscore_result.ld_scores.loc[:, selected_columns].copy(),
        regression_metadata=ldscore_result.regression_metadata.copy(),
        w_ld=ldscore_result.w_ld.copy(),
        snp_count_totals=snp_count_totals,
        baseline_columns=list(baseline_columns),
        query_columns=list(query_columns),
        reference_snps=set(ldscore_result.reference_snps),
        regression_snps=set(ldscore_result.regression_snps),
        chromosome_results=list(ldscore_result.chromosome_results),
        output_paths=dict(ldscore_result.output_paths),
    )
