"""Workflow-layer LD-score orchestration and normalized result objects.

This module is the public boundary for LD-score path handling and result
normalization. Callers may provide exact paths, glob patterns, or explicit
chromosome-suite tokens using ``@``. The workflow resolves those tokens into
concrete per-chromosome files, aligns each chromosome-local annotation bundle
to the prepared reference-panel metadata returned by ``ref_panel.load_metadata``,
and only then dispatches to the primitive-only kernel.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass, field
from typing import Any, Sequence
import warnings

import numpy as np
import pandas as pd

from .column_inference import normalize_snp_identifier_mode
from .config import (
    ConfigMismatchError,
    GlobalConfig,
    AnnotationBuildConfig,
    LDScoreConfig,
    RefPanelConfig,
    get_global_config,
    print_global_config_banner,
    validate_config_compatibility,
)
from .genome_build_inference import validate_auto_genome_build_mode
from .outputs import LDScoreDirectoryWriter, LDScoreOutputConfig
from .path_resolution import (
    FREQUENCY_SUFFIXES,
    PARQUET_SUFFIXES,
    normalize_optional_path_token,
    normalize_path_token,
    resolve_chromosome_group,
    resolve_plink_prefix,
    resolve_scalar_path,
    split_cli_path_tokens,
)
from ._kernel import ldscore as kernel_ldscore
from ._kernel.identifiers import build_snp_id_series, read_global_snp_restriction


# Temporary compatibility alias while tests and call sites finish moving from the
# old wrapper naming.
ldscore_new = kernel_ldscore


_LDSCORE_SUFFIX_COLUMNS = ("CHR", "SNP", "POS", "BP", "CM", "MAF")


@dataclass(frozen=True)
class _LegacyChromResult:
    """Compatibility record for chromosome outputs produced by the kernel."""
    chrom: str
    metadata: pd.DataFrame
    ld_scores: np.ndarray
    w_ld: np.ndarray
    M: np.ndarray
    M_5_50: np.ndarray | None
    ldscore_columns: list[str]
    baseline_columns: list[str]
    query_columns: list[str]


@dataclass(frozen=True)
class ChromLDScoreResult:
    """Normalized per-chromosome LD-score output in split baseline/query form.

    Parameters
    ----------
    chrom : str
        Chromosome label.
    baseline_table : pandas.DataFrame
        Table with ``CHR``, ``SNP``, ``BP``, ``regr_weight``, and baseline
        LD-score columns.
    query_table : pandas.DataFrame or None
        Optional table with ``CHR``, ``SNP``, ``BP``, and query LD-score
        columns.
    count_records : list of dict
        Manifest-ready count records keyed by annotation column.
    baseline_columns, query_columns : list of str
        Ordered annotation LD-score columns in the baseline and query tables.
    ld_reference_snps, ld_regression_snps : frozenset of str
        SNP universes used for count records and regression rows.
    """
    chrom: str
    baseline_table: pd.DataFrame
    query_table: pd.DataFrame | None
    count_records: list[dict[str, Any]]
    baseline_columns: list[str]
    query_columns: list[str]
    ld_reference_snps: frozenset[str]
    ld_regression_snps: frozenset[str]
    snp_count_totals: dict[str, np.ndarray] = field(default_factory=dict, repr=False)
    output_paths: dict[str, str] = field(default_factory=dict)
    config_snapshot: GlobalConfig | None = None

    def validate(self) -> None:
        """Check the normalized public contract for chromosome-level results."""
        required = {"CHR", "SNP", "BP", "regr_weight", *self.baseline_columns}
        missing = required - set(self.baseline_table.columns)
        if missing:
            raise ValueError(f"baseline_table is missing required columns: {sorted(missing)}")
        if self.query_columns and self.query_table is None:
            raise ValueError("query_table is required when query_columns are present.")
        if self.query_table is not None:
            missing_query = set(self.query_columns) - set(self.query_table.columns)
            if missing_query:
                raise ValueError(f"query_table is missing required columns: {sorted(missing_query)}")
            baseline_keys = self.baseline_table.loc[:, ["CHR", "SNP", "BP"]].reset_index(drop=True)
            query_keys = self.query_table.loc[:, ["CHR", "SNP", "BP"]].reset_index(drop=True)
            if not baseline_keys.equals(query_keys):
                raise ValueError("query rows must match baseline rows on CHR/SNP/BP.")

    def summary(self) -> dict[str, Any]:
        """Return a compact summary of chromosome-level retained rows."""
        return {
            "chrom": self.chrom,
            "n_rows": len(self.baseline_table),
            "count_columns": [record["column"] for record in self.count_records],
        }


@dataclass(frozen=True)
class LDScoreResult:
    """Aggregated cross-chromosome LD-score result in split persisted form.

    Parameters
    ----------
    baseline_table : pandas.DataFrame
        Cross-chromosome table persisted as ``baseline.parquet`` when the result
        is written. Required columns are ``CHR``, ``SNP``, ``BP``,
        ``regr_weight``, and every entry in ``baseline_columns``.
    query_table : pandas.DataFrame or None
        Optional cross-chromosome table persisted as ``query.parquet``. Required
        columns are ``CHR``, ``SNP``, ``BP``, and every entry in
        ``query_columns``.
    count_records : list of dict
        Manifest-ready count records. Each record names an annotation column and
        its all-SNP and optional common-SNP counts.
    baseline_columns, query_columns : list of str
        Ordered LD-score columns available to regression workflows.
    ld_reference_snps, ld_regression_snps : frozenset of str
        SNP universes used for count records and persisted regression rows.
    chromosome_results : list of ChromLDScoreResult
        Per-chromosome components used to build the aggregate result.
    output_paths : dict, optional
        Paths written by ``LDScoreDirectoryWriter`` when output was requested.
    config_snapshot : GlobalConfig or None, optional
        Shared configuration active when the result was computed.
    """
    baseline_table: pd.DataFrame
    query_table: pd.DataFrame | None
    count_records: list[dict[str, Any]]
    baseline_columns: list[str]
    query_columns: list[str]
    ld_reference_snps: frozenset[str]
    ld_regression_snps: frozenset[str]
    chromosome_results: list[ChromLDScoreResult]
    output_paths: dict[str, str] = field(default_factory=dict)
    config_snapshot: GlobalConfig | None = None

    def validate(self) -> None:
        """Check the normalized public contract for aggregated results."""
        required = {"CHR", "SNP", "BP", "regr_weight", *self.baseline_columns}
        missing = required - set(self.baseline_table.columns)
        if missing:
            raise ValueError(f"baseline_table is missing required columns: {sorted(missing)}")
        if self.query_columns and self.query_table is None:
            raise ValueError("query_table is required when query_columns are present.")
        if not self.query_columns and self.query_table is not None:
            raise ValueError("query_table was provided but query_columns is empty.")
        if self.query_table is not None:
            missing_query = set(self.query_columns) - set(self.query_table.columns)
            if missing_query:
                raise ValueError(f"query_table is missing required columns: {sorted(missing_query)}")
            baseline_keys = self.baseline_table.loc[:, ["CHR", "SNP", "BP"]].reset_index(drop=True)
            query_keys = self.query_table.loc[:, ["CHR", "SNP", "BP"]].reset_index(drop=True)
            if not baseline_keys.equals(query_keys):
                raise ValueError("query rows must match baseline rows on CHR/SNP/BP.")

    def summary(self) -> dict[str, Any]:
        """Return a compact cross-chromosome summary."""
        return {
            "n_rows": len(self.baseline_table),
            "chromosomes": [result.chrom for result in self.chromosome_results],
            "count_columns": [record["column"] for record in self.count_records],
        }


class LDScoreCalculator:
    """Orchestrate chromosome-wise LD-score calculation.

    This service assembles annotation and reference-panel inputs, delegates the
    heavy computation to the internal LD-score kernel, aggregates chromosome
    outputs, and optionally hands the result to the output layer. The calculator
    treats ``ref_panel`` as the owner of the reference-panel SNP universe:
    chromosome bundles are intersected with ``ref_panel.load_metadata(chrom)``
    before the kernel runs so ``RefPanelConfig.ref_panel_snps_path`` is honored
    without leaking that setting into the calculator interface. For LD-score
    calculation, the annotation file's ``CM`` is the first source. The sidecar
    metadata only fills missing ``CM`` values.
    """

    def __init__(self, output_writer: LDScoreDirectoryWriter | None = None) -> None:
        """Initialize the calculator with the directory writer used for LD-score outputs."""
        self.output_writer = output_writer or LDScoreDirectoryWriter()

    def run(
        self,
        annotation_bundle,
        ref_panel,
        ldscore_config: LDScoreConfig,
        global_config: GlobalConfig,
        output_config: LDScoreOutputConfig | None = None,
        regression_snps: set[str] | None = None,
        config_snapshot: dict[str, Any] | None = None,
    ) -> LDScoreResult:
        """Compute and aggregate LD scores across all chromosomes.

        Parameters
        ----------
        annotation_bundle : AnnotationBundle
            Aligned SNP-level baseline and query annotations.
        ref_panel : RefPanel
            Reference-panel adapter that supplies chromosome readers and
            metadata. Any ``RefPanelConfig.ref_panel_snps_path`` restriction is
            already applied when the workflow calls ``ref_panel.load_metadata``.
            Reference-panel sidecar ``CM`` values are used only to fill missing
            annotation ``CM`` values during LD-score calculation.
        ldscore_config : LDScoreConfig
            LD-window and retained-SNP settings.
        global_config : GlobalConfig
            Shared identifier and restriction settings.
        output_config : LDScoreOutputConfig or None, optional
            If provided, write the canonical LD-score result directory after
            the aggregate result is built.
            Default is ``None``, which keeps the result in memory only.
        regression_snps : set of str or None, optional
            Optional regression SNP universe used to define the weight table.
            Default is ``None``, which uses the retained reference SNP universe.
        config_snapshot : dict or None, optional
            Optional run metadata forwarded to the output layer. Default is
            ``None``.

        Returns
        -------
        LDScoreResult
            Aggregated cross-chromosome result with aligned metadata and output
            paths if writing was requested.
        """
        print_global_config_banner(type(self).__name__, global_config)
        if annotation_bundle.config_snapshot is not None:
            validate_config_compatibility(
                annotation_bundle.config_snapshot,
                global_config,
                context="AnnotationBundle and LDScoreCalculator runtime config",
            )
        ref_panel_build = getattr(getattr(ref_panel, "spec", None), "genome_build", None)
        if ref_panel_build is not None and global_config.genome_build not in (None, "auto"):
            if ref_panel_build != global_config.genome_build:
                raise ConfigMismatchError(
                    f"genome_build mismatch between RefPanelConfig ({ref_panel_build!r}) "
                    f"and active GlobalConfig ({global_config.genome_build!r})."
                )
        chromosome_results: list[ChromLDScoreResult] = []
        for chrom in _chromosomes_from_bundle(annotation_bundle):
            chrom_bundle = _slice_annotation_bundle(annotation_bundle, chrom)
            try:
                chromosome_result = self.compute_chromosome(
                    chrom=chrom,
                    annotation_bundle=chrom_bundle,
                    ref_panel=ref_panel,
                    ldscore_config=ldscore_config,
                    global_config=global_config,
                    regression_snps=regression_snps,
                )
            except ValueError as exc:
                if _warn_and_skip_empty_intersection(exc, chrom):
                    continue
                raise
            chromosome_results.append(chromosome_result)
        if not chromosome_results:
            raise ValueError("No chromosome results were produced after intersecting annotations with the reference panel.")
        result = self._aggregate_chromosome_results(chromosome_results, global_config=global_config)
        if output_config is not None:
            output_paths = self.output_writer.write(result, output_config)
            result = _replace_result_output_paths(result, output_paths)
        return result

    def compute_chromosome(
        self,
        chrom: str,
        annotation_bundle,
        ref_panel,
        ldscore_config: LDScoreConfig,
        global_config: GlobalConfig,
        regression_snps: set[str] | None = None,
    ) -> ChromLDScoreResult:
        """Compute normalized LD-score outputs for one chromosome.

        The workflow first asks ``ref_panel`` for the prepared chromosome
        metadata, which already reflects any ``ref_panel_snps_path`` restriction.
        It then restricts the chromosome-local ``AnnotationBundle`` to the same
        identifier set so the legacy kernel sees ``B_chrom ∩ A'_chrom`` rather
        than the raw annotation universe ``B_chrom``. ``regression_snps`` is
        applied later when the normalized row table is formed.
        """
        backend = getattr(getattr(ref_panel, "spec", None), "backend", None)
        if backend == "parquet_r2" and ldscore_config.keep_indivs_path is not None:
            raise ValueError("keep_indivs_path/--keep-indivs-path is only supported for PLINK reference panels.")
        annotation_bundle = _align_annotation_bundle_to_ref_panel(
            annotation_bundle=annotation_bundle,
            ref_panel=ref_panel,
            chrom=chrom,
            global_config=global_config,
        )
        args = _namespace_from_configs(
            chrom=chrom,
            ref_panel=ref_panel,
            ldscore_config=ldscore_config,
            global_config=global_config,
        )
        legacy_bundle = kernel_ldscore.AnnotationBundle(
            metadata=annotation_bundle.metadata.copy(),
            annotations=annotation_bundle.annotation_matrix(include_query=True).copy(),
            baseline_columns=list(annotation_bundle.baseline_columns),
            query_columns=list(annotation_bundle.query_columns),
        )
        if getattr(ref_panel.spec, "backend", None) == "parquet_r2":
            legacy_result = kernel_ldscore.compute_chrom_from_parquet(chrom, legacy_bundle, args, regression_snps)
        else:
            legacy_result = kernel_ldscore.compute_chrom_from_plink(chrom, legacy_bundle, args, regression_snps)
        return self._wrap_legacy_chrom_result(legacy_result, global_config=global_config, regression_snps=regression_snps)

    def _wrap_legacy_chrom_result(
        self,
        legacy_result: _LegacyChromResult | Any,
        global_config: GlobalConfig,
        regression_snps: set[str] | None = None,
    ) -> ChromLDScoreResult:
        """Convert one kernel chromosome result into the typed public result."""
        reference_metadata = legacy_result.metadata.reset_index(drop=True).copy()
        ld_scores = pd.DataFrame(legacy_result.ld_scores, columns=list(legacy_result.ldscore_columns))
        reference_ids = frozenset(build_snp_id_series(reference_metadata, global_config.snp_identifier))
        retained_regression_snps = (
            reference_ids if regression_snps is None else frozenset(reference_ids.intersection(regression_snps))
        )
        regression_keep = build_snp_id_series(reference_metadata, global_config.snp_identifier).isin(retained_regression_snps)
        bp_column = "BP" if "BP" in reference_metadata.columns else "POS"
        regression_weights = np.asarray(legacy_result.w_ld, dtype=np.float32).reshape(-1)
        ldscore_table = pd.concat(
            [
                reference_metadata.loc[regression_keep, ["CHR", "SNP", bp_column]].rename(columns={bp_column: "BP"}).reset_index(drop=True),
                ld_scores.loc[regression_keep].reset_index(drop=True),
                pd.DataFrame({"regr_weight": regression_weights[regression_keep.to_numpy()]}).reset_index(drop=True),
            ],
            axis=1,
        )
        ldscore_table = kernel_ldscore.sort_frame_by_genomic_position(ldscore_table)
        count_map = {"all_reference_snp_counts": np.asarray(legacy_result.M, dtype=np.float64)}
        if legacy_result.M_5_50 is not None:
            count_map["common_reference_snp_counts_maf_gt_0_05"] = np.asarray(legacy_result.M_5_50, dtype=np.float64)
        baseline_table, query_table = _split_ldscore_table(
            ldscore_table,
            baseline_columns=list(legacy_result.baseline_columns),
            query_columns=list(legacy_result.query_columns),
        )
        result = ChromLDScoreResult(
            chrom=str(legacy_result.chrom),
            baseline_table=baseline_table,
            query_table=query_table,
            count_records=_count_records_from_totals(
                baseline_columns=list(legacy_result.baseline_columns),
                query_columns=list(legacy_result.query_columns),
                count_totals=count_map,
            ),
            baseline_columns=list(legacy_result.baseline_columns),
            query_columns=list(legacy_result.query_columns),
            ld_reference_snps=frozenset(),
            ld_regression_snps=frozenset(build_snp_id_series(baseline_table, global_config.snp_identifier)),
            snp_count_totals=count_map,
            config_snapshot=global_config,
        )
        result.validate()
        return result

    def _aggregate_chromosome_results(
        self,
        chromosome_results: Sequence[ChromLDScoreResult],
        global_config: GlobalConfig,
    ) -> LDScoreResult:
        """Concatenate and sum per-chromosome results into one aggregate object."""
        if not chromosome_results:
            raise ValueError("At least one chromosome result is required.")
        snapshots = [result.config_snapshot for result in chromosome_results if result.config_snapshot is not None]
        if snapshots:
            for snapshot in snapshots[1:]:
                validate_config_compatibility(
                    snapshots[0],
                    snapshot,
                    context="ChromLDScoreResult aggregation",
                )

        count_keys = sorted({key for result in chromosome_results for key in result.snp_count_totals})
        count_totals = {
            key: np.sum(np.vstack([result.snp_count_totals[key] for result in chromosome_results if key in result.snp_count_totals]), axis=0)
            for key in count_keys
        }
        merged_table = pd.concat(
            [_join_split_tables(result.baseline_table, result.query_table, result.query_columns) for result in chromosome_results],
            axis=0,
            ignore_index=True,
        )
        merged_table = kernel_ldscore.sort_frame_by_genomic_position(merged_table)
        baseline_table, query_table = _split_ldscore_table(
            merged_table,
            baseline_columns=list(chromosome_results[0].baseline_columns),
            query_columns=list(chromosome_results[0].query_columns),
        )
        result = LDScoreResult(
            baseline_table=baseline_table,
            query_table=query_table,
            count_records=_count_records_from_totals(
                baseline_columns=list(chromosome_results[0].baseline_columns),
                query_columns=list(chromosome_results[0].query_columns),
                count_totals=count_totals,
            ),
            baseline_columns=list(chromosome_results[0].baseline_columns),
            query_columns=list(chromosome_results[0].query_columns),
            ld_reference_snps=frozenset(),
            ld_regression_snps=frozenset().union(*(result.ld_regression_snps for result in chromosome_results)),
            chromosome_results=list(chromosome_results),
            config_snapshot=snapshots[0] if snapshots else None,
        )
        result.validate()
        return result

    def write_outputs(
        self,
        result: LDScoreResult,
        output_config: LDScoreOutputConfig,
        config_snapshot: dict[str, Any] | None = None,
    ):
        """Write a previously computed result as a canonical LD-score directory.

        Parameters
        ----------
        result : LDScoreResult
            Aggregate LD-score result to serialize.
        output_config : LDScoreOutputConfig
            Directory path and overwrite/compression controls.
        config_snapshot : dict or None, optional
            Optional metadata recorded in the emitted run summary. Default is
            ``None``.

        Returns
        -------
        dict
            Resolved output paths keyed by artifact name.
        """
        del config_snapshot
        return self.output_writer.write(result, output_config)


def _split_ldscore_table(
    ldscore_table: pd.DataFrame,
    *,
    baseline_columns: list[str],
    query_columns: list[str],
) -> tuple[pd.DataFrame, pd.DataFrame | None]:
    """Split a merged LD-score table into baseline and optional query tables."""
    metadata_columns = ["CHR", "SNP", "BP"]
    baseline_order = [*metadata_columns, "regr_weight", *baseline_columns]
    query_order = [*metadata_columns, *query_columns]
    missing_baseline = [column for column in baseline_order if column not in ldscore_table.columns]
    if missing_baseline:
        raise ValueError(f"LD-score table is missing baseline columns: {missing_baseline}")
    baseline_table = ldscore_table.loc[:, baseline_order].reset_index(drop=True).copy()
    query_table = None
    if query_columns:
        missing_query = [column for column in query_order if column not in ldscore_table.columns]
        if missing_query:
            raise ValueError(f"LD-score table is missing query columns: {missing_query}")
        query_table = ldscore_table.loc[:, query_order].reset_index(drop=True).copy()
    return baseline_table, query_table


def _join_split_tables(
    baseline_table: pd.DataFrame,
    query_table: pd.DataFrame | None,
    query_columns: Sequence[str],
) -> pd.DataFrame:
    """Join split LD-score tables for sorting or regression assembly."""
    if query_table is None:
        return baseline_table.copy()
    baseline_keys = baseline_table.loc[:, ["CHR", "SNP", "BP"]].reset_index(drop=True)
    query_keys = query_table.loc[:, ["CHR", "SNP", "BP"]].reset_index(drop=True)
    if not baseline_keys.equals(query_keys):
        raise ValueError("query rows must match baseline rows on CHR/SNP/BP.")
    query_values = query_table.loc[:, list(query_columns)].reset_index(drop=True)
    return pd.concat([baseline_table.reset_index(drop=True), query_values], axis=1)


def _count_records_from_totals(
    *,
    baseline_columns: list[str],
    query_columns: list[str],
    count_totals: dict[str, np.ndarray],
) -> list[dict[str, Any]]:
    """Convert positional count vectors into manifest-friendly column records."""
    columns = [*baseline_columns, *query_columns]
    groups = ["baseline"] * len(baseline_columns) + ["query"] * len(query_columns)
    all_counts = np.asarray(count_totals.get("all_reference_snp_counts"), dtype=np.float64)
    if all_counts.size != len(columns):
        raise ValueError("all_reference_snp_counts length does not match annotation columns.")
    common_raw = count_totals.get("common_reference_snp_counts_maf_gt_0_05")
    common_counts = None if common_raw is None else np.asarray(common_raw, dtype=np.float64)
    if common_counts is not None and common_counts.size != len(columns):
        raise ValueError("common_reference_snp_counts_maf_gt_0_05 length does not match annotation columns.")
    records: list[dict[str, Any]] = []
    for idx, (group, column) in enumerate(zip(groups, columns)):
        record: dict[str, Any] = {
            "group": group,
            "column": column,
            "all_reference_snp_count": float(all_counts[idx]),
        }
        if common_counts is not None:
            record["common_reference_snp_count_maf_gt_0_05"] = float(common_counts[idx])
        records.append(record)
    return records


def build_parser() -> argparse.ArgumentParser:
    """Build the feature parser for LD-score calculation."""
    parser = argparse.ArgumentParser(
        description="Estimate LDSC-compatible LD scores from SNP-level annotation files using PLINK or sorted parquet R2 input.",
        allow_abbrev=False,
    )
    parser.add_argument("--output-dir", required=True, help="Output directory for the canonical LD-score result.")
    query_group = parser.add_mutually_exclusive_group()
    query_group.add_argument(
        "--query-annot-paths",
        default=None,
        help="Comma-separated query annotation path tokens: exact paths, globs, or explicit @ suite tokens.",
    )
    query_group.add_argument(
        "--query-annot-bed-paths",
        default=None,
        help="Comma-separated BED file path tokens projected in memory as query annotations.",
    )
    parser.add_argument("--baseline-annot-paths", default=None, help="Comma-separated baseline annotation path tokens: exact paths, globs, or explicit @ suite tokens.")
    parser.add_argument("--plink-path", default=None, help="PLINK prefix token for the reference panel.")
    parser.add_argument("--r2-paths", default=None, help="Comma-separated parquet R2 path tokens: exact paths, globs, or explicit @ suite tokens.")
    parser.add_argument("--snp-identifier", default="chr_pos", help="Identifier mode used to match annotations to the reference panel.")
    parser.add_argument(
        "--genome-build",
        choices=("auto", "hg19", "hg37", "GRCh37", "hg38", "GRCh38"),
        default=None,
        help="Genome build for chr_pos matching. Use 'auto' to infer hg19/hg38 and 0-based/1-based coordinates.",
    )
    parser.add_argument("--r2-bias-mode", choices=("raw", "unbiased"), default=None, help="Whether parquet R2 values are raw or already unbiased.")
    parser.add_argument("--r2-sample-size", default=None, type=float, help="LD reference sample size used to correct raw parquet R2 values.")
    parser.add_argument(
        "--ref-panel-snps-path",
        default=None,
        help="Optional SNP list defining the retained reference-panel universe A'; the workflow intersects each chromosome annotation bundle with this prepared panel before LD computation.",
    )
    parser.add_argument("--regression-snps-path", default=None, help="Optional SNP list defining the regression SNP set and the written LD-score row set.")
    parser.add_argument("--metadata-paths", default=None, help="Optional frequency or metadata path tokens for MAF and CM.")
    parser.add_argument(
        "--keep-indivs-path",
        default=None,
        help="File with individuals to include in LD Score estimation. The file should contain one IID per row.",
    )
    parser.add_argument("--ld-wind-snps", default=None, type=int, help="LD window size in SNPs.")
    parser.add_argument("--ld-wind-kb", default=None, type=float, help="LD window size in kilobases.")
    parser.add_argument("--ld-wind-cm", default=None, type=float, help="LD window size in centiMorgans.")
    parser.add_argument("--maf", default=None, type=float, help="Optional MAF filter for retained SNPs when MAF is available.")
    parser.add_argument("--chunk-size", default=50, type=int, help="Chunk size for legacy PLINK block computations.")
    parser.add_argument("--yes-really", default=False, action="store_true", help="Allow whole-chromosome LD windows.")
    parser.add_argument("--log-level", default="INFO", choices=("DEBUG", "INFO", "WARNING", "ERROR"), help="Logging verbosity.")
    return parser


def run_ldscore_from_args(args: argparse.Namespace) -> LDScoreResult:
    """Run LD-score calculation from a parsed CLI namespace.

    The workflow resolves unified path tokens for baseline annotations, optional
    query annotations or BED files, and the reference panel. For each
    chromosome it intersects the annotation rows with
    ``ref_panel.load_metadata(chrom)`` before calling the kernel, then returns
    the normalized public ``LDScoreResult`` with split baseline/query tables.
    """
    from ._kernel.annotation import AnnotationBuilder

    normalized_args, global_config = _normalize_run_args(args)
    print_global_config_banner("run_ldscore_from_args", global_config)
    kernel_ldscore.validate_args(normalized_args)
    ldscore_config = _ldscore_config_from_args(normalized_args)
    regression_snps = _load_regression_snps(ldscore_config.regression_snps_path, global_config)
    source_spec = AnnotationBuildConfig(
        baseline_annot_paths=tuple(split_cli_path_tokens(normalized_args.baseline_annot_paths)),
        query_annot_paths=tuple(split_cli_path_tokens(normalized_args.query_annot_paths)),
        query_annot_bed_paths=tuple(split_cli_path_tokens(getattr(normalized_args, "query_annot_bed_paths", None))),
    )
    annotation_bundle = AnnotationBuilder(global_config).run(source_spec)
    ref_panel = _ref_panel_from_args(normalized_args, global_config)
    calculator = LDScoreCalculator()
    return calculator.run(
        annotation_bundle=annotation_bundle,
        ref_panel=ref_panel,
        ldscore_config=ldscore_config,
        global_config=global_config,
        output_config=_output_config_from_args(normalized_args),
        regression_snps=regression_snps,
    )


def run_ldscore(**kwargs) -> LDScoreResult:
    """Run LD-score calculation from Python using public CLI-style names.

    Keyword arguments are interpreted as CLI-equivalent option names without
    leading ``--``; for example ``baseline_annot_paths``, ``query_annot_paths``,
    ``query_annot_bed_paths``, ``plink_path``, ``r2_paths``,
    ``metadata_paths``, ``keep_indivs_path``, and ``output_dir``. Shared
    runtime assumptions such as ``snp_identifier`` and ``genome_build`` must be
    supplied through ``set_global_config(...)`` first, while per-run controls
    such as ``ref_panel_snps_path`` and ``regression_snps_path`` remain ordinary
    keyword arguments here.

    Returns
    -------
    LDScoreResult
        Aggregated result with ``baseline_table``, optional ``query_table``,
        manifest count records, and canonical output paths.
    """
    forbidden = sorted({"snp_identifier", "genome_build", "log_level"} & set(kwargs))
    if forbidden:
        joined = ", ".join(forbidden)
        raise ValueError(f"Python run_ldscore() no longer accepts {joined}; call set_global_config(...) first.")
    removed = sorted(
        {
            "out",
            "baseline_annot",
            "query_annot",
            "query_annot_bed",
            "bfile",
            "r2_table",
            "frqfile",
            "keep",
        }
        & set(kwargs)
    )
    if removed:
        joined = ", ".join(removed)
        raise ValueError(f"Python run_ldscore() no longer accepts removed IO argument(s): {joined}.")
    parser = build_parser()
    defaults = vars(parser.parse_args(["--output-dir", "placeholder"]))
    global_config = get_global_config()
    defaults["snp_identifier"] = global_config.snp_identifier
    defaults["genome_build"] = global_config.genome_build
    defaults["log_level"] = global_config.log_level
    defaults.update(kwargs)
    args = argparse.Namespace(**defaults)
    return run_ldscore_from_args(args)


def main(argv: Sequence[str] | None = None) -> LDScoreResult:
    """Command-line entry point for the LD-score workflow."""
    parser = build_parser()
    args = parser.parse_args(argv)
    return run_ldscore_from_args(args)


def _normalize_run_args(args: argparse.Namespace) -> tuple[argparse.Namespace, GlobalConfig]:
    """Normalize CLI-style args and derive the shared ``GlobalConfig`` object."""
    normalized_mode = normalize_snp_identifier_mode(args.snp_identifier)
    validate_auto_genome_build_mode(normalized_mode, getattr(args, "genome_build", None))
    normalized_args = argparse.Namespace(**vars(args))
    for attr in ("query_annot_chr", "baseline_annot_chr", "bfile_chr", "r2_table_chr", "frqfile_chr"):
        if not hasattr(normalized_args, attr):
            setattr(normalized_args, attr, None)
    for attr in ("query_annot_paths", "baseline_annot_paths", "plink_path", "r2_paths", "metadata_paths", "query_annot_bed_paths", "keep_indivs_path"):
        if not hasattr(normalized_args, attr):
            setattr(normalized_args, attr, None)
    for attr in ("ref_panel_snps_path", "regression_snps_path"):
        if not hasattr(normalized_args, attr):
            setattr(normalized_args, attr, None)
    normalized_args.snp_identifier = normalized_mode
    normalized_args.output_dir = normalize_path_token(args.output_dir)
    normalized_args.keep_indivs_path = normalize_optional_path_token(getattr(args, "keep_indivs_path", None))
    normalized_args.ref_panel_snps_path = normalize_optional_path_token(getattr(args, "ref_panel_snps_path", None))
    normalized_args.regression_snps_path = normalize_optional_path_token(getattr(args, "regression_snps_path", None))
    # The numerical kernel still consumes the historical namespace shape.
    normalized_args.query_annot = normalized_args.query_annot_paths
    normalized_args.baseline_annot = normalized_args.baseline_annot_paths
    normalized_args.bfile = normalized_args.plink_path
    normalized_args.r2_table = normalized_args.r2_paths
    normalized_args.frqfile = normalized_args.metadata_paths
    normalized_args.keep = normalized_args.keep_indivs_path
    default_config = GlobalConfig()
    global_config = GlobalConfig(
        snp_identifier=normalized_mode,
        genome_build=default_config.genome_build if getattr(args, "genome_build", None) is None else getattr(args, "genome_build"),
        log_level=getattr(args, "log_level", "INFO"),
    )
    return normalized_args, global_config


def _load_regression_snps(path: str | None, global_config: GlobalConfig) -> set[str] | None:
    """Load ``LDScoreConfig.regression_snps_path`` using the active identifier mode."""
    if not path:
        return None
    return read_global_snp_restriction(
        resolve_scalar_path(path, label="regression SNP list"),
        global_config.snp_identifier,
        genome_build=global_config.genome_build,
    )


def _ref_panel_from_args(args: argparse.Namespace, global_config: GlobalConfig):
    """Build the reference-panel adapter that owns the ``A -> A'`` restriction."""
    from ._kernel.ref_panel import RefPanelLoader

    metadata_tokens = split_cli_path_tokens(getattr(args, "metadata_paths", None))
    ref_panel_snps_path = normalize_optional_path_token(getattr(args, "ref_panel_snps_path", None))
    if getattr(args, "r2_paths", None) is not None:
        r2_tokens = split_cli_path_tokens(args.r2_paths)
        spec = RefPanelConfig(
            backend="parquet_r2",
            r2_paths=tuple(r2_tokens),
            metadata_paths=tuple(metadata_tokens),
            r2_bias_mode=getattr(args, "r2_bias_mode", None),
            sample_size=getattr(args, "r2_sample_size", None),
            genome_build=global_config.genome_build,
            ref_panel_snps_path=ref_panel_snps_path,
        )
    else:
        spec = RefPanelConfig(
            backend="plink",
            plink_path=getattr(args, "plink_path", None),
            metadata_paths=tuple(metadata_tokens),
            genome_build=global_config.genome_build,
            ref_panel_snps_path=ref_panel_snps_path,
        )
    return RefPanelLoader(global_config).load(spec)


def _ldscore_config_from_args(args: argparse.Namespace) -> LDScoreConfig:
    """Build an ``LDScoreConfig`` from normalized LD-score arguments."""
    return LDScoreConfig(
        ld_wind_snps=getattr(args, "ld_wind_snps", None),
        ld_wind_kb=getattr(args, "ld_wind_kb", None),
        ld_wind_cm=getattr(args, "ld_wind_cm", None),
        maf_min=getattr(args, "maf", None),
        keep_indivs_path=getattr(args, "keep_indivs_path", None),
        regression_snps_path=getattr(args, "regression_snps_path", None),
        chunk_size=getattr(args, "chunk_size", 50),
        whole_chromosome_ok=getattr(args, "yes_really", False),
    )


def _output_config_from_args(args: argparse.Namespace) -> LDScoreOutputConfig:
    """Translate LD-score CLI arguments into the canonical directory output config."""
    return LDScoreOutputConfig(output_dir=normalize_path_token(args.output_dir))


def _replace_result_output_paths(result: LDScoreResult, output_paths: dict[str, str]) -> LDScoreResult:
    """Return ``result`` with updated artifact-path metadata after writing outputs."""
    return LDScoreResult(
        baseline_table=result.baseline_table,
        query_table=result.query_table,
        count_records=result.count_records,
        baseline_columns=result.baseline_columns,
        query_columns=result.query_columns,
        ld_reference_snps=result.ld_reference_snps,
        ld_regression_snps=result.ld_regression_snps,
        chromosome_results=result.chromosome_results,
        output_paths=dict(output_paths),
        config_snapshot=result.config_snapshot,
    )


def _chromosomes_from_bundle(annotation_bundle) -> list[str]:
    """Infer the chromosome processing order from an annotation bundle."""
    chromosomes = getattr(annotation_bundle, "chromosomes", None)
    if chromosomes:
        return list(chromosomes)
    return sorted(annotation_bundle.metadata["CHR"].astype(str).unique().tolist())


def _slice_annotation_bundle(annotation_bundle, chrom: str):
    """Return the per-chromosome view of an annotation bundle."""
    keep = annotation_bundle.metadata["CHR"].astype(str) == str(chrom)
    return type(annotation_bundle)(
        metadata=annotation_bundle.metadata.loc[keep].reset_index(drop=True),
        baseline_annotations=annotation_bundle.baseline_annotations.loc[keep].reset_index(drop=True),
        query_annotations=annotation_bundle.query_annotations.loc[keep].reset_index(drop=True),
        baseline_columns=list(annotation_bundle.baseline_columns),
        query_columns=list(annotation_bundle.query_columns),
        chromosomes=[str(chrom)],
        source_summary=dict(getattr(annotation_bundle, "source_summary", {})),
        config_snapshot=getattr(annotation_bundle, "config_snapshot", None),
    )


def _align_annotation_bundle_to_ref_panel(annotation_bundle, ref_panel, chrom: str, global_config: GlobalConfig):
    """Restrict one chromosome bundle to the prepared reference-panel universe.

    ``AnnotationBuilder`` defines the annotation universe ``B``. The reference
    panel owns the optional ``A -> A'`` restriction through
    ``RefPanelConfig.ref_panel_snps_path``. This helper materializes the intended
    compute-time universe ``B_chrom ∩ A'_chrom`` immediately before the legacy
    kernel call. It restricts annotation rows but does not replace annotation
    metadata; for LD-score calculation, the annotation file's ``CM`` is the
    first source, and sidecar metadata only fills missing ``CM`` values later in
    the kernel.
    """
    reference_metadata = ref_panel.load_metadata(chrom)
    reference_ids = set(build_snp_id_series(reference_metadata, global_config.snp_identifier))
    keep = build_snp_id_series(annotation_bundle.metadata, global_config.snp_identifier).isin(reference_ids)
    if not bool(keep.any()):
        backend = getattr(getattr(ref_panel, "spec", None), "backend", None)
        intersection = "parquet" if backend == "parquet_r2" else "PLINK"
        raise ValueError(f"No retained annotation SNPs remain on chromosome {chrom} after {intersection} intersection.")
    if bool(keep.all()):
        return annotation_bundle
    return type(annotation_bundle)(
        metadata=annotation_bundle.metadata.loc[keep].reset_index(drop=True),
        baseline_annotations=annotation_bundle.baseline_annotations.loc[keep].reset_index(drop=True),
        query_annotations=annotation_bundle.query_annotations.loc[keep].reset_index(drop=True),
        baseline_columns=list(annotation_bundle.baseline_columns),
        query_columns=list(annotation_bundle.query_columns),
        chromosomes=list(getattr(annotation_bundle, "chromosomes", [str(chrom)])),
        source_summary=dict(getattr(annotation_bundle, "source_summary", {})),
        config_snapshot=getattr(annotation_bundle, "config_snapshot", None),
    )


def _warn_and_skip_empty_intersection(error: ValueError, chrom: str) -> bool:
    """Warn and signal skip when a chromosome loses all SNPs after reference intersection."""
    message = str(error)
    if not message.startswith("No retained annotation SNPs remain on chromosome "):
        return False
    if not any(suffix in message for suffix in (" after parquet intersection.", " after PLINK intersection.")):
        return False
    warnings.warn(f"Skipping chromosome {chrom}: {message}", UserWarning, stacklevel=3)
    return True


def _namespace_from_configs(chrom: str, ref_panel, ldscore_config: LDScoreConfig, global_config: GlobalConfig) -> argparse.Namespace:
    """Build the legacy-kernel namespace expected by the chromosome backends."""
    spec = getattr(ref_panel, "spec", None)
    backend = getattr(spec, "backend", None)
    bfile = None
    r2_table = None
    frqfile = None
    if backend == "plink" and getattr(spec, "plink_path", None) is not None:
        bfile = resolve_plink_prefix(spec.plink_path, chrom=chrom)
    if backend == "parquet_r2" and getattr(spec, "r2_paths", ()):
        r2_table = ",".join(
            resolve_chromosome_group(
                getattr(spec, "r2_paths", ()),
                chrom=chrom,
                suffixes=PARQUET_SUFFIXES,
                label="parquet R2",
                required=False,
            )
        ) or None
    if getattr(spec, "metadata_paths", ()):
        frqfile = ",".join(
            resolve_chromosome_group(
                getattr(spec, "metadata_paths", ()),
                chrom=chrom,
                suffixes=FREQUENCY_SUFFIXES,
                label="frequency or metadata",
                required=False,
            )
        ) or None
    return argparse.Namespace(
        out=None,
        query_annot=None,
        query_annot_chr=None,
        baseline_annot=None,
        baseline_annot_chr=None,
        bfile=bfile,
        bfile_chr=None,
        r2_table=r2_table,
        r2_table_chr=None,
        snp_identifier=global_config.snp_identifier,
        genome_build=(getattr(spec, "genome_build", None) or global_config.genome_build),
        r2_bias_mode=getattr(spec, "r2_bias_mode", None),
        r2_sample_size=getattr(spec, "sample_size", None),
        frqfile=frqfile,
        frqfile_chr=None,
        keep=ldscore_config.keep_indivs_path,
        ld_wind_snps=ldscore_config.ld_wind_snps,
        ld_wind_kb=ldscore_config.ld_wind_kb,
        ld_wind_cm=ldscore_config.ld_wind_cm,
        maf=ldscore_config.maf_min,
        chunk_size=ldscore_config.chunk_size,
        per_chr_output=False,
        yes_really=ldscore_config.whole_chromosome_ok,
        log_level=global_config.log_level,
    )
