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
"""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass
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
from .path_resolution import ensure_output_parent_directory, normalize_path_token
from ._kernel import regression as reg
from ._kernel.identifiers import build_snp_id_series
from .ldscore_calculator import LDScoreResult
from .outputs import LDSCORE_RESULT_FORMAT
from .sumstats_munger import SumstatsTable, load_sumstats


COMMON_COUNT_KEY = "common_reference_snp_counts_maf_gt_0_05"
ALL_COUNT_KEY = "all_reference_snp_counts"


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

        Zero-variance LD-score columns are dropped here so the estimator kernel
        receives only informative regressors. The selected count vector is
        carried alongside the merged table for later use by ``Hsq`` and ``RG``.
        """
        print_global_config_banner(type(self).__name__, self.global_config)
        config = config or self.regression_config
        if sumstats_table.config_snapshot is not None and ldscore_result.config_snapshot is not None:
            validate_config_compatibility(
                sumstats_table.config_snapshot,
                ldscore_result.config_snapshot,
                context="SumstatsTable and LDScoreResult",
            )
        weight_column = "regr_weight"
        selected_query_columns = list(query_columns or [])
        ref_ld_columns = list(ldscore_result.baseline_columns) + selected_query_columns
        ldscore_frame = _assemble_regression_ldscore_table(ldscore_result, selected_query_columns)
        merged = pd.merge(
            sumstats_table.data,
            ldscore_frame.loc[:, ["SNP", *ref_ld_columns, weight_column]].reset_index(drop=True),
            how="inner",
            on="SNP",
            sort=False,
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
        count_key = _select_count_key(count_totals, config.use_m_5_50)
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
        annotation_bundle,
        config: RegressionConfig | None = None,
    ) -> pd.DataFrame:
        """Estimate partitioned heritability for one selected query annotation."""
        del annotation_bundle  # grouping metadata is carried by the selected columns
        query_column = ldscore_result.query_columns[-1] if ldscore_result.query_columns else ldscore_result.baseline_columns[-1]
        selected_queries = [query_column] if query_column in ldscore_result.query_columns else []
        dataset = self.build_dataset(sumstats_table, ldscore_result, config=config, query_columns=selected_queries)
        hsq = self.estimate_h2(dataset, config=config)
        return summarize_partitioned_h2(hsq, dataset, [query_column])

    def estimate_partitioned_h2_batch(
        self,
        sumstats_table: SumstatsTable,
        ldscore_result: LDScoreResult,
        annotation_bundle,
        config: RegressionConfig | None = None,
    ) -> pd.DataFrame:
        """Loop over query annotations and concatenate partitioned-h2 summaries."""
        rows = []
        for query_column in annotation_bundle.query_columns:
            dataset = self.build_dataset(sumstats_table, ldscore_result, config=config, query_columns=[query_column])
            hsq = self.estimate_h2(dataset, config=config)
            query_row = summarize_partitioned_h2(hsq, dataset, [query_column])
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
        """Estimate genetic correlation between two munged summary-stat tables."""
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
            keep = reg._filter_alleles(alleles)
            merged = merged.loc[keep].reset_index(drop=True)
            merged["Z2"] = reg._align_alleles(merged["Z2"], alleles.loc[keep])
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

def _select_count_key(count_totals: dict[str, np.ndarray], use_m_5_50: bool) -> str:
    """Pick the regression count vector key, preferring LDSC's common-SNP default."""
    if use_m_5_50 and COMMON_COUNT_KEY in count_totals:
        return COMMON_COUNT_KEY
    if ALL_COUNT_KEY in count_totals:
        return ALL_COUNT_KEY
    return sorted(count_totals.keys())[0]


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
    baseline_keys = baseline_table.loc[:, ["CHR", "SNP", "BP"]]
    query_keys = query_table.loc[:, ["CHR", "SNP", "BP"]]
    if not baseline_keys.equals(query_keys):
        raise ValueError("query rows must match baseline rows on CHR/SNP/BP.")
    return pd.concat([baseline_table, query_table.loc[:, list(query_columns)]], axis=1)


def _count_totals_for_columns(count_records: Sequence[dict[str, Any]], columns: Sequence[str]) -> dict[str, np.ndarray]:
    """Return LDSC count vectors aligned to ``columns`` using manifest records."""
    records_by_column = {str(record["column"]): record for record in count_records}
    missing = [column for column in columns if column not in records_by_column]
    if missing:
        raise ValueError(f"LD-score manifest is missing count records for columns: {missing}")
    all_counts = [float(records_by_column[column]["all_reference_snp_count"]) for column in columns]
    count_totals = {ALL_COUNT_KEY: np.asarray(all_counts, dtype=np.float64)}
    if all("common_reference_snp_count_maf_gt_0_05" in records_by_column[column] for column in columns):
        count_totals[COMMON_COUNT_KEY] = np.asarray(
            [float(records_by_column[column]["common_reference_snp_count_maf_gt_0_05"]) for column in columns],
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


def summarize_partitioned_h2(hsq, dataset: RegressionDataset, annotation_columns: Sequence[str]) -> pd.DataFrame:
    """Build category-level heritability rows from a fitted ``Hsq`` result."""
    rows = []
    coefficients = np.ravel(hsq.coef)
    coefficient_ses = np.ravel(hsq.coef_se)
    category_h2 = np.ravel(hsq.cat)
    category_h2_ses = np.ravel(hsq.cat_se)
    proportions = np.ravel(hsq.prop)
    proportion_ses = np.ravel(hsq.prop_se)
    enrichments = np.ravel(hsq.enrichment) if getattr(hsq, "enrichment", None) is not None else None
    for annotation_column in annotation_columns:
        annotation_index = dataset.retained_ld_columns.index(annotation_column)
        coefficient = float(coefficients[annotation_index])
        coefficient_se = float(coefficient_ses[annotation_index])
        z_value = math.nan if coefficient_se == 0 else coefficient / coefficient_se
        rows.append(
            {
                "query_annotation": annotation_column,
                "coefficient": coefficient,
                "coefficient_se": coefficient_se,
                "coefficient_z": z_value,
                "coefficient_p": (
                    math.nan if coefficient_se == 0 else 2 * stats.norm.sf(abs(z_value))
                ),
                "category_h2": float(category_h2[annotation_index]),
                "category_h2_se": float(category_h2_ses[annotation_index]),
                "proportion_h2": float(proportions[annotation_index]),
                "proportion_h2_se": float(proportion_ses[annotation_index]),
                "enrichment": math.nan if enrichments is None else float(enrichments[annotation_index]),
            }
        )
    return pd.DataFrame(rows)


def _scalar(value) -> float:
    """Return the first scalar numeric value from a possibly array-like object."""
    return float(np.ravel(value)[0])


def _scalar_or_value(value):
    """Return a scalar float for array-like values, otherwise preserve the original value."""
    if hasattr(value, "__array__") or isinstance(value, (list, tuple)):
        return _scalar(value)
    return value


def _select_intercept(value, index: int, use_intercept: bool, default_when_disabled: float):
    """Resolve fixed-intercept config into the scalar expected by the kernel."""
    if not use_intercept:
        return default_when_disabled
    if value is None:
        return None
    if isinstance(value, list):
        return float(value[index])
    return float(value)


def add_h2_arguments(parser) -> None:
    """Register heritability CLI arguments on ``parser``."""
    _add_common_regression_arguments(parser, include_h2_intercept=True)
    parser.add_argument("--sumstats", required=True, help="Munged .sumstats(.gz) file.")
    parser.add_argument("--trait-name", default=None, help="Optional trait label for summaries.")


def add_partitioned_h2_arguments(parser) -> None:
    """Register partitioned-heritability CLI arguments on ``parser``."""
    _add_common_regression_arguments(parser, include_h2_intercept=True)
    parser.add_argument("--sumstats", required=True, help="Munged .sumstats(.gz) file.")
    parser.add_argument("--trait-name", default=None, help="Optional trait label for summaries.")


def add_rg_arguments(parser) -> None:
    """Register genetic-correlation CLI arguments on ``parser``."""
    _add_common_regression_arguments(parser, include_h2_intercept=False)
    parser.add_argument("--sumstats-1", required=True, help="First munged .sumstats(.gz) file.")
    parser.add_argument("--sumstats-2", required=True, help="Second munged .sumstats(.gz) file.")
    parser.add_argument("--trait-name-1", default=None, help="Optional label for the first trait.")
    parser.add_argument("--trait-name-2", default=None, help="Optional label for the second trait.")
    parser.add_argument("--intercept-h2", nargs=2, type=float, default=None, metavar=("H2_1", "H2_2"))
    parser.add_argument("--intercept-gencov", type=float, default=None)


def run_h2_from_args(args):
    """Run single-trait heritability estimation from parsed CLI arguments."""
    runner, config = _runner_from_args(args)
    print_global_config_banner("run_h2_from_args", runner.global_config)
    sumstats_table = _load_sumstats_table(args.sumstats, getattr(args, "trait_name", None))
    ldscore_result = load_ldscore_from_dir(args.ldscore_dir)
    with suppress_global_config_banner():
        dataset = runner.build_dataset(sumstats_table, ldscore_result, config=config)
    hsq = runner.estimate_h2(dataset, config=config)
    summary = summarize_total_h2(hsq, dataset, trait_name=sumstats_table.trait_name)
    _maybe_write_dataframe(summary, args.out, ".h2.tsv")
    return summary


def run_partitioned_h2_from_args(args):
    """Run batch partitioned heritability from parsed CLI arguments."""
    runner, config = _runner_from_args(args)
    print_global_config_banner("run_partitioned_h2_from_args", runner.global_config)
    sumstats_table = _load_sumstats_table(args.sumstats, getattr(args, "trait_name", None))
    ldscore_result = load_ldscore_from_dir(args.ldscore_dir)
    if not ldscore_result.query_columns:
        raise ValueError("partitioned-h2 requires query annotations in --ldscore-dir.")
    query_bundle = SimpleNamespace(query_columns=list(ldscore_result.query_columns))
    with suppress_global_config_banner():
        summary = runner.estimate_partitioned_h2_batch(sumstats_table, ldscore_result, query_bundle, config=config)
    _maybe_write_dataframe(summary, args.out, ".partitioned_h2.tsv")
    return summary


def run_rg_from_args(args):
    """Run genetic-correlation estimation from parsed CLI arguments."""
    runner, config = _runner_from_args(args)
    print_global_config_banner("run_rg_from_args", runner.global_config)
    sumstats_table_1 = _load_sumstats_table(args.sumstats_1, getattr(args, "trait_name_1", None))
    sumstats_table_2 = _load_sumstats_table(args.sumstats_2, getattr(args, "trait_name_2", None))
    ldscore_result = load_ldscore_from_dir(args.ldscore_dir)
    with suppress_global_config_banner():
        rg_result = runner.estimate_rg(sumstats_table_1, sumstats_table_2, ldscore_result, config=config)
    summary = pd.DataFrame(
        [
            {
                "trait_1": sumstats_table_1.trait_name,
                "trait_2": sumstats_table_2.trait_name,
                "rg": getattr(rg_result, "rg_ratio", None),
                "rg_se": getattr(rg_result, "rg_se", None),
                "z": getattr(rg_result, "z", None),
                "p": getattr(rg_result, "p", None),
            }
        ]
    )
    _maybe_write_dataframe(summary, args.out, ".rg.tsv")
    return summary


def _add_common_regression_arguments(parser, include_h2_intercept: bool) -> None:
    """Add the file and model options shared by all regression subcommands."""
    parser.add_argument("--ldscore-dir", required=True, help="Canonical LD-score result directory written by `ldsc ldscore`.")
    parser.add_argument(
        "--count-kind",
        choices=("m_5_50", "all"),
        default="m_5_50",
        help="Interpretation of the supplied count vector. Default matches LDSC's M_5_50 behavior.",
    )
    parser.add_argument("--out", default=None, help="Optional output prefix for summary tables.")
    parser.add_argument("--n-blocks", type=int, default=200)
    parser.add_argument("--no-intercept", action="store_true", default=False, help="Fix the intercept to the LDSC default.")
    if include_h2_intercept:
        parser.add_argument("--intercept-h2", type=float, default=None, help="Fixed h2 intercept for single-trait runs.")
    parser.add_argument("--two-step-cutoff", type=float, default=None)
    parser.add_argument("--chisq-max", type=float, default=None)


def _runner_from_args(args) -> tuple[RegressionRunner, RegressionConfig]:
    """Build the regression workflow objects from parsed CLI arguments."""
    count_kind = getattr(args, "count_kind", "m_5_50")
    config = RegressionConfig(
        n_blocks=args.n_blocks,
        use_m_5_50=(count_kind == "m_5_50"),
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
    """Load a canonical LD-score result directory."""
    root = Path(normalize_path_token(ldscore_dir))
    if not root.is_dir():
        raise NotADirectoryError(f"LD-score directory does not exist or is not a directory: {root}")
    manifest_path = root / "manifest.json"
    if not manifest_path.exists():
        raise FileNotFoundError(f"LD-score directory is missing manifest.json: {root}")
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
        config_snapshot=config_snapshot,
    )
    result.validate()
    return result


def _global_config_from_manifest(manifest: dict[str, Any]) -> GlobalConfig | None:
    """Recreate a GlobalConfig snapshot when the manifest contains one."""
    snapshot = manifest.get("config_snapshot")
    if not isinstance(snapshot, dict):
        return None
    try:
        return GlobalConfig(
            snp_identifier=snapshot.get("snp_identifier") or manifest.get("snp_identifier") or "rsid",
            genome_build=snapshot.get("genome_build") or manifest.get("genome_build"),
            log_level=snapshot.get("log_level", "INFO"),
            fail_on_missing_metadata=bool(snapshot.get("fail_on_missing_metadata", False)),
        )
    except Exception:
        return None


def _maybe_write_dataframe(df: pd.DataFrame, out_prefix: str | None, suffix: str) -> None:
    """Write a summary table only when the caller requested an output prefix."""
    if not out_prefix:
        return
    path = Path(normalize_path_token(out_prefix) + suffix)
    path = ensure_output_parent_directory(path, label="out_prefix")
    df.to_csv(path, sep="\t", index=False)
