"""Workflow-layer regression adapters and dataset assembly.

Regression file inputs stay at the workflow boundary. Public callers may pass
literal paths or exact-one input tokens such as globs; those tokens are
resolved before tabular artifacts are loaded, while output prefixes remain
literal destinations.
"""

from __future__ import annotations

import argparse
import math
from dataclasses import dataclass, field
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import numpy as np
import pandas as pd
from scipy import stats

from .config import CommonConfig, RegressionConfig
from .path_resolution import normalize_path_token, resolve_scalar_path
from ._kernel import regression as reg
from .ldscore_calculator import LDScoreResult
from .sumstats_munger import SumstatsTable


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
        """Merge sumstats, LD scores, and weights into a regression dataset.

        Zero-variance LD-score columns are dropped here so the estimator kernel
        receives only informative regressors. The selected count vector is
        carried alongside the merged table for later use by ``Hsq`` and ``RG``.
        """
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
        """Loop over query annotations and concatenate partitioned-h2 summaries."""
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


def _select_intercept(value, index: int, use_intercept: bool, default_when_disabled: float):
    """Resolve fixed-intercept config into the scalar expected by the kernel."""
    if not use_intercept:
        return default_when_disabled
    if value is None:
        return None
    if isinstance(value, list):
        return float(value[index])
    return float(value)


def _subset_ldscore_result(ldscore_result: LDScoreResult, baseline_columns: list[str], query_columns: list[str]) -> LDScoreResult:
    """Return a copy of ``ldscore_result`` restricted to the requested columns."""
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
    parser.add_argument(
        "--annotation-manifest",
        default=None,
        help="Annotation-group manifest with columns 'column' and 'group'.",
    )
    parser.add_argument(
        "--query-columns",
        default=None,
        help="Comma-separated query annotation columns when no manifest is available.",
    )


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
    sumstats_table = _load_sumstats_table(args.sumstats, getattr(args, "trait_name", None))
    ldscore_result = _load_ldscore_result_from_files(
        ldscore_path=args.ldscore,
        weight_path=args.w_ld,
        counts_path=args.counts,
        count_kind=args.count_kind,
        annotation_manifest_path=getattr(args, "annotation_manifest", None),
    )
    dataset = runner.build_dataset(sumstats_table, ldscore_result, config=config)
    hsq = runner.estimate_h2(dataset, config=config)
    summary = pd.DataFrame(
        [
            {
                "trait_name": sumstats_table.trait_name,
                "n_snps": len(dataset.merged),
                "total_h2": float(np.ravel(hsq.tot)[0]),
                "total_h2_se": float(np.ravel(hsq.tot_se)[0]),
                "intercept": float(np.ravel(hsq.intercept)[0]) if hasattr(hsq.intercept, "__array__") else hsq.intercept,
                "intercept_se": getattr(hsq, "intercept_se", None),
                "mean_chisq": float(np.ravel(hsq.mean_chisq)[0]) if hasattr(hsq.mean_chisq, "__array__") else hsq.mean_chisq,
                "lambda_gc": float(np.ravel(hsq.lambda_gc)[0]) if hasattr(hsq.lambda_gc, "__array__") else hsq.lambda_gc,
                "ratio": getattr(hsq, "ratio", None),
                "ratio_se": getattr(hsq, "ratio_se", None),
            }
        ]
    )
    _maybe_write_dataframe(summary, args.out, ".h2.tsv")
    return summary


def run_partitioned_h2_from_args(args):
    """Run batch partitioned heritability from parsed CLI arguments."""
    runner, config = _runner_from_args(args)
    sumstats_table = _load_sumstats_table(args.sumstats, getattr(args, "trait_name", None))
    ldscore_result = _load_ldscore_result_from_files(
        ldscore_path=args.ldscore,
        weight_path=args.w_ld,
        counts_path=args.counts,
        count_kind=args.count_kind,
        annotation_manifest_path=getattr(args, "annotation_manifest", None),
        explicit_query_columns=getattr(args, "query_columns", None),
    )
    if not ldscore_result.query_columns:
        raise ValueError("partitioned-h2 requires query annotation columns via --annotation-manifest or --query-columns.")
    query_bundle = SimpleNamespace(query_columns=list(ldscore_result.query_columns))
    summary = runner.estimate_partitioned_h2_batch(sumstats_table, ldscore_result, query_bundle, config=config)
    _maybe_write_dataframe(summary, args.out, ".partitioned_h2.tsv")
    return summary


def run_rg_from_args(args):
    """Run genetic-correlation estimation from parsed CLI arguments."""
    runner, config = _runner_from_args(args)
    sumstats_table_1 = _load_sumstats_table(args.sumstats_1, getattr(args, "trait_name_1", None))
    sumstats_table_2 = _load_sumstats_table(args.sumstats_2, getattr(args, "trait_name_2", None))
    ldscore_result = _load_ldscore_result_from_files(
        ldscore_path=args.ldscore,
        weight_path=args.w_ld,
        counts_path=args.counts,
        count_kind=args.count_kind,
    )
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
    parser.add_argument("--ldscore", required=True, help="Reference LD-score table (.l2.ldscore[.gz]).")
    parser.add_argument("--w-ld", required=True, help="Regression-weight LD-score table (.w.l2.ldscore[.gz]).")
    parser.add_argument("--counts", required=True, help="Annotation count vector file (.M or .M_5_50).")
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
    config = RegressionConfig(
        n_blocks=args.n_blocks,
        use_m_5_50=(args.count_kind == "m_5_50"),
        use_intercept=not args.no_intercept,
        intercept_h2=args.intercept_h2,
        intercept_gencov=getattr(args, "intercept_gencov", None),
        two_step_cutoff=args.two_step_cutoff,
        chisq_max=args.chisq_max,
    )
    runner = RegressionRunner(CommonConfig(), config)
    return runner, config


def _read_table(path: str) -> pd.DataFrame:
    """Read a whitespace-delimited LDSC artifact table with gzip support."""
    resolved = resolve_scalar_path(path, label="tabular artifact")
    compression = "gzip" if str(resolved).endswith(".gz") else "infer"
    return pd.read_csv(resolved, sep=r"\s+", compression=compression)


def _read_count_vector(path: str) -> np.ndarray:
    """Read one count-vector artifact such as ``.l2.M`` or ``.l2.M_5_50``."""
    resolved = resolve_scalar_path(path, label="count vector")
    text = Path(resolved).read_text(encoding="utf-8").strip()
    if not text:
        raise ValueError(f"Count file is empty: {resolved}")
    return np.asarray([float(value) for value in text.split()], dtype=np.float64)


def _load_sumstats_table(path: str, trait_name: str | None) -> SumstatsTable:
    """Load one munged sumstats file and wrap it in ``SumstatsTable``."""
    resolved = resolve_scalar_path(path, label="munged sumstats")
    df = _read_table(resolved)
    table = SumstatsTable(
        data=df.reset_index(drop=True),
        has_alleles={"A1", "A2"}.issubset(df.columns),
        source_path=resolved,
        trait_name=trait_name or Path(resolved).name,
    )
    table.validate()
    return table


def _load_ldscore_result_from_files(
    ldscore_path: str,
    weight_path: str,
    counts_path: str,
    count_kind: str,
    annotation_manifest_path: str | None = None,
    explicit_query_columns: str | None = None,
) -> LDScoreResult:
    """
    Rebuild an ``LDScoreResult`` from previously written LD-score artifacts.

    This keeps the regression CLI file-driven while reusing the in-memory
    workflow objects defined by the refactored package.
    """
    ld_table = _read_table(ldscore_path)
    weight_table = _read_table(weight_path)
    metadata_columns = [column for column in ["CHR", "SNP", "POS", "BP", "CM", "MAF"] if column in ld_table.columns]
    weight_metadata_columns = [column for column in ["CHR", "SNP", "POS", "BP", "CM", "MAF"] if column in weight_table.columns]
    if "SNP" not in metadata_columns or "SNP" not in weight_metadata_columns:
        raise ValueError("Both ldscore and w-ld tables must include an SNP column.")

    reference_metadata = ld_table.loc[:, metadata_columns].reset_index(drop=True)
    regression_metadata = weight_table.loc[:, weight_metadata_columns].reset_index(drop=True)
    if "BP" in reference_metadata.columns and "POS" not in reference_metadata.columns:
        reference_metadata = reference_metadata.rename(columns={"BP": "POS"})
    if "BP" in regression_metadata.columns and "POS" not in regression_metadata.columns:
        regression_metadata = regression_metadata.rename(columns={"BP": "POS"})
    ld_scores = ld_table.drop(columns=metadata_columns).reset_index(drop=True)
    w_ld = weight_table.drop(columns=weight_metadata_columns).reset_index(drop=True)
    if w_ld.shape[1] != 1:
        raise ValueError("The regression-weight table must contain exactly one non-metadata column.")

    baseline_columns, query_columns = _resolve_annotation_groups(
        ld_columns=list(ld_scores.columns),
        annotation_manifest_path=annotation_manifest_path,
        explicit_query_columns=explicit_query_columns,
    )
    count_key = COMMON_COUNT_KEY if count_kind == "m_5_50" else ALL_COUNT_KEY
    count_totals = {count_key: _read_count_vector(counts_path)}
    return LDScoreResult(
        reference_metadata=reference_metadata,
        ld_scores=ld_scores,
        regression_metadata=regression_metadata,
        w_ld=w_ld,
        snp_count_totals=count_totals,
        baseline_columns=baseline_columns,
        query_columns=query_columns,
        reference_snps=set(reference_metadata["SNP"].astype(str)),
        regression_snps=set(regression_metadata["SNP"].astype(str)),
        chromosome_results=[],
        output_paths={},
    )


def _resolve_annotation_groups(
    ld_columns: list[str],
    annotation_manifest_path: str | None,
    explicit_query_columns: str | None,
) -> tuple[list[str], list[str]]:
    """Resolve baseline/query columns from a manifest or explicit query list."""
    if annotation_manifest_path:
        manifest = _read_table(resolve_scalar_path(annotation_manifest_path, label="annotation manifest"))
        if not {"column", "group"}.issubset(manifest.columns):
            raise ValueError("Annotation manifest must contain 'column' and 'group' columns.")
        baseline = manifest.loc[manifest["group"] == "baseline", "column"].astype(str).tolist()
        query = manifest.loc[manifest["group"] == "query", "column"].astype(str).tolist()
        return [column for column in baseline if column in ld_columns], [column for column in query if column in ld_columns]
    if explicit_query_columns:
        query = [column.strip() for column in explicit_query_columns.split(",") if column.strip()]
        baseline = [column for column in ld_columns if column not in query]
        return baseline, query
    return list(ld_columns), []


def _maybe_write_dataframe(df: pd.DataFrame, out_prefix: str | None, suffix: str) -> None:
    """Write a summary table only when the caller requested an output prefix."""
    if not out_prefix:
        return
    path = Path(normalize_path_token(out_prefix) + suffix)
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False)
