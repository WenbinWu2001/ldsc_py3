"""Workflow-layer LD-score orchestration and result objects.

This module is the public boundary for LD-score path handling. Callers may
provide exact paths, glob patterns, explicit chromosome-suite tokens using
``@``, or legacy bare prefixes. The workflow resolves those tokens into
concrete per-chromosome files before calling the primitive-only kernel.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Sequence

import numpy as np
import pandas as pd

from .column_inference import normalize_snp_identifier_mode
from .config import CommonConfig, LDScoreConfig
from .outputs import OutputManager, OutputSpec
from .path_resolution import (
    ANNOTATION_SUFFIXES,
    FREQUENCY_SUFFIXES,
    PARQUET_SUFFIXES,
    normalize_optional_path_token,
    normalize_path_token,
    resolve_chromosome_group,
    resolve_file_group,
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
    """Chromosome-specific LD-score outputs with aligned metadata tables.

    The LD-score pipeline is executed chromosome by chromosome for efficiency.
    This dataclass keeps the reference LD scores, regression-weight LD scores,
    and count vectors for one chromosome while preserving the row-order
    guarantees needed by the aggregated result.
    """
    chrom: str
    reference_metadata: pd.DataFrame
    ld_scores: pd.DataFrame
    regression_metadata: pd.DataFrame
    w_ld: pd.DataFrame
    snp_count_totals: dict[str, np.ndarray]
    baseline_columns: list[str]
    query_columns: list[str]
    reference_snps: set[str]
    regression_snps: set[str]
    output_paths: dict[str, str] = field(default_factory=dict)

    def validate(self) -> None:
        """Check metadata/value alignment for the chromosome result."""
        if len(self.reference_metadata) != len(self.ld_scores):
            raise ValueError("reference_metadata and ld_scores must have the same number of rows.")
        if len(self.regression_metadata) != len(self.w_ld):
            raise ValueError("regression_metadata and w_ld must have the same number of rows.")

    def to_ldscore_table(self) -> pd.DataFrame:
        """Return ``reference_metadata`` and ``ld_scores`` as one table."""
        return pd.concat([self.reference_metadata.reset_index(drop=True), self.ld_scores.reset_index(drop=True)], axis=1)

    def to_weight_table(self) -> pd.DataFrame:
        """Return ``regression_metadata`` and ``w_ld`` as one table."""
        return pd.concat([self.regression_metadata.reset_index(drop=True), self.w_ld.reset_index(drop=True)], axis=1)

    def summary(self) -> dict[str, Any]:
        """Return a compact summary of chromosome-level retained rows."""
        return {
            "chrom": self.chrom,
            "n_reference_rows": len(self.reference_metadata),
            "n_regression_rows": len(self.regression_metadata),
            "count_keys": sorted(self.snp_count_totals.keys()),
        }


@dataclass(frozen=True)
class LDScoreResult:
    """Aggregated cross-chromosome LD-score result used by regression.

    Parameters
    ----------
    reference_metadata : pandas.DataFrame
        Metadata aligned row-for-row with ``ld_scores``.
    ld_scores : pandas.DataFrame
        Aggregated LD-score columns for retained reference SNPs.
    regression_metadata : pandas.DataFrame
        Metadata aligned row-for-row with ``w_ld``.
    w_ld : pandas.DataFrame
        Regression-weight LD-score table.
    snp_count_totals : dict of str to ndarray
        Named count vectors such as all-reference counts and common-SNP counts.
    baseline_columns, query_columns : list of str
        Annotation-column names carried through the run.
    reference_snps, regression_snps : set of str
        Retained SNP identifiers in the reference and regression universes.
    chromosome_results : list of ChromLDScoreResult
        Per-chromosome intermediate results used to build this aggregate.
    output_paths : dict of str to str, optional
        Paths written by the output layer, when available. Default is an empty
        dict.
    """
    reference_metadata: pd.DataFrame
    ld_scores: pd.DataFrame
    regression_metadata: pd.DataFrame
    w_ld: pd.DataFrame
    snp_count_totals: dict[str, np.ndarray]
    baseline_columns: list[str]
    query_columns: list[str]
    reference_snps: set[str]
    regression_snps: set[str]
    chromosome_results: list[ChromLDScoreResult]
    output_paths: dict[str, str] = field(default_factory=dict)

    def validate(self) -> None:
        """Check metadata/value alignment for the aggregated result."""
        if len(self.reference_metadata) != len(self.ld_scores):
            raise ValueError("reference_metadata and ld_scores must have the same number of rows.")
        if len(self.regression_metadata) != len(self.w_ld):
            raise ValueError("regression_metadata and w_ld must have the same number of rows.")

    def to_ldscore_table(self) -> pd.DataFrame:
        """Return the aggregated reference LD-score table."""
        return pd.concat([self.reference_metadata.reset_index(drop=True), self.ld_scores.reset_index(drop=True)], axis=1)

    def to_weight_table(self) -> pd.DataFrame:
        """Return the aggregated regression-weight table."""
        return pd.concat([self.regression_metadata.reset_index(drop=True), self.w_ld.reset_index(drop=True)], axis=1)

    def summary(self) -> dict[str, Any]:
        """Return a compact cross-chromosome summary."""
        return {
            "n_reference_rows": len(self.reference_metadata),
            "n_regression_rows": len(self.regression_metadata),
            "chromosomes": [result.chrom for result in self.chromosome_results],
            "count_keys": sorted(self.snp_count_totals.keys()),
        }


class LDScoreCalculator:
    """Orchestrate chromosome-wise LD-score calculation.

    This service assembles annotation and reference-panel inputs, delegates the
    heavy computation to the internal LD-score kernel, aggregates chromosome
    outputs, and optionally hands the result to the output layer.
    """

    def __init__(self, output_manager: OutputManager | None = None) -> None:
        self.output_manager = output_manager or OutputManager()

    def run(
        self,
        annotation_bundle,
        ref_panel,
        ldscore_config: LDScoreConfig,
        common_config: CommonConfig,
        output_spec: OutputSpec | None = None,
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
            metadata.
        ldscore_config : LDScoreConfig
            LD-window and retained-SNP settings.
        common_config : CommonConfig
            Shared identifier and restriction settings.
        output_spec : OutputSpec or None, optional
            If provided, write outputs after the aggregate result is built.
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
        chromosome_results: list[ChromLDScoreResult] = []
        for chrom in _chromosomes_from_bundle(annotation_bundle):
            chrom_bundle = _slice_annotation_bundle(annotation_bundle, chrom)
            chromosome_results.append(
                self.compute_chromosome(
                    chrom=chrom,
                    annotation_bundle=chrom_bundle,
                    ref_panel=ref_panel,
                    ldscore_config=ldscore_config,
                    common_config=common_config,
                    regression_snps=regression_snps,
                )
            )
        result = self._aggregate_chromosome_results(chromosome_results, common_config=common_config)
        if output_spec is not None:
            summary = self.output_manager.write_outputs(result, output_spec, config_snapshot=config_snapshot)
            result = _replace_result_output_paths(result, summary.output_paths)
        return result

    def compute_chromosome(
        self,
        chrom: str,
        annotation_bundle,
        ref_panel,
        ldscore_config: LDScoreConfig,
        common_config: CommonConfig,
        regression_snps: set[str] | None = None,
    ) -> ChromLDScoreResult:
        """Compute LD scores for one chromosome."""
        backend = getattr(getattr(ref_panel, "spec", None), "backend", None)
        if backend == "parquet_r2" and ldscore_config.keep_individuals_path is not None:
            raise ValueError("keep_individuals_path/--keep is only supported for PLINK reference panels.")
        args = _namespace_from_configs(
            chrom=chrom,
            ref_panel=ref_panel,
            ldscore_config=ldscore_config,
            common_config=common_config,
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
        return self._wrap_legacy_chrom_result(legacy_result, common_config=common_config, regression_snps=regression_snps)

    def _wrap_legacy_chrom_result(
        self,
        legacy_result: _LegacyChromResult | Any,
        common_config: CommonConfig,
        regression_snps: set[str] | None = None,
    ) -> ChromLDScoreResult:
        reference_metadata = legacy_result.metadata.reset_index(drop=True).copy()
        ld_scores = pd.DataFrame(legacy_result.ld_scores, columns=list(legacy_result.ldscore_columns))
        regression_metadata = reference_metadata.copy()
        w_ld = pd.DataFrame(np.asarray(legacy_result.w_ld, dtype=np.float32), columns=["L2"])
        reference_ids = set(build_snp_id_series(reference_metadata, common_config.snp_identifier))
        retained_regression_snps = reference_ids if regression_snps is None else reference_ids.intersection(regression_snps)
        count_map = {"all_reference_snp_counts": np.asarray(legacy_result.M, dtype=np.float64)}
        if legacy_result.M_5_50 is not None:
            count_map["common_reference_snp_counts_maf_gt_0_05"] = np.asarray(legacy_result.M_5_50, dtype=np.float64)
        result = ChromLDScoreResult(
            chrom=str(legacy_result.chrom),
            reference_metadata=reference_metadata,
            ld_scores=ld_scores,
            regression_metadata=regression_metadata,
            w_ld=w_ld,
            snp_count_totals=count_map,
            baseline_columns=list(legacy_result.baseline_columns),
            query_columns=list(legacy_result.query_columns),
            reference_snps=reference_ids,
            regression_snps=retained_regression_snps,
        )
        result.validate()
        return result

    def _aggregate_chromosome_results(
        self,
        chromosome_results: Sequence[ChromLDScoreResult],
        common_config: CommonConfig,
    ) -> LDScoreResult:
        if not chromosome_results:
            raise ValueError("At least one chromosome result is required.")

        reference_metadata, ld_scores = _sort_paired_frames(
            [result.reference_metadata for result in chromosome_results],
            [result.ld_scores for result in chromosome_results],
        )
        regression_metadata, w_ld = _sort_paired_frames(
            [result.regression_metadata for result in chromosome_results],
            [result.w_ld for result in chromosome_results],
        )
        count_keys = sorted({key for result in chromosome_results for key in result.snp_count_totals})
        count_totals = {
            key: np.sum(np.vstack([result.snp_count_totals[key] for result in chromosome_results if key in result.snp_count_totals]), axis=0)
            for key in count_keys
        }
        result = LDScoreResult(
            reference_metadata=reference_metadata,
            ld_scores=ld_scores,
            regression_metadata=regression_metadata,
            w_ld=w_ld,
            snp_count_totals=count_totals,
            baseline_columns=list(chromosome_results[0].baseline_columns),
            query_columns=list(chromosome_results[0].query_columns),
            reference_snps=set().union(*(result.reference_snps for result in chromosome_results)),
            regression_snps=set().union(*(result.regression_snps for result in chromosome_results)),
            chromosome_results=list(chromosome_results),
        )
        result.validate()
        return result

    def write_outputs(
        self,
        result: LDScoreResult,
        output_spec: OutputSpec,
        config_snapshot: dict[str, Any] | None = None,
    ):
        """Write a previously computed result through the output layer."""
        return self.output_manager.write_outputs(result, output_spec, config_snapshot=config_snapshot)


def build_parser() -> argparse.ArgumentParser:
    """Build the feature parser for LD-score calculation."""
    parser = argparse.ArgumentParser(
        description="Estimate LDSC-compatible LD scores from SNP-level annotation files using PLINK or sorted parquet R2 input.",
    )
    parser.add_argument("--out", required=True, help="Output prefix.")
    parser.add_argument("--query-annot", default=None, help="Comma-separated query annotation tokens: exact paths, globs, or suite prefixes.")
    parser.add_argument("--query-annot-chr", default=None, help="Comma-separated chromosome-suite query annotation tokens; `@` is preferred.")
    parser.add_argument("--baseline-annot", default=None, help="Comma-separated baseline annotation tokens: exact paths, globs, or suite prefixes.")
    parser.add_argument("--baseline-annot-chr", default=None, help="Comma-separated chromosome-suite baseline annotation tokens; `@` is preferred.")
    parser.add_argument("--bfile", default=None, help="PLINK prefix token for the reference panel.")
    parser.add_argument("--bfile-chr", default=None, help="Chromosome-suite PLINK prefix token; `@` is preferred.")
    parser.add_argument("--r2-table", default=None, help="Comma-separated parquet R2 tokens: exact paths, globs, or suite prefixes.")
    parser.add_argument("--r2-table-chr", default=None, help="Comma-separated chromosome-suite parquet R2 tokens; `@` is preferred.")
    parser.add_argument("--snp-identifier", default="chr_pos", help="Identifier mode used to match annotations to the reference panel.")
    parser.add_argument(
        "--genome-build",
        choices=("hg19", "hg37", "GRCh37", "hg38", "GRCh38"),
        default=None,
        help="Genome build for chr_pos matching. Aliases normalize to canonical hg19/hg38.",
    )
    parser.add_argument("--r2-bias-mode", choices=("raw", "unbiased"), default=None, help="Whether parquet R2 values are raw or already unbiased.")
    parser.add_argument("--r2-sample-size", default=None, type=float, help="LD reference sample size used to correct raw parquet R2 values.")
    parser.add_argument("--regression-snps", default=None, help="Optional SNP list defining the regression SNP set for w_ld.")
    parser.add_argument("--frqfile", default=None, help="Optional frequency or metadata token for MAF and CM.")
    parser.add_argument("--frqfile-chr", default=None, help="Chromosome-suite frequency or metadata token; `@` is preferred.")
    parser.add_argument(
        "--keep",
        default=None,
        help="File with individuals to include in LD Score estimation. The file should contain one IID per row.",
    )
    parser.add_argument(
        "--print-snps",
        default=None,
        help=(
            "Only print LD scores for the SNPs listed in this one-column file. "
            "LD-score computation itself still includes SNPs not listed here."
        ),
    )
    parser.add_argument("--ld-wind-snps", default=None, type=int, help="LD window size in SNPs.")
    parser.add_argument("--ld-wind-kb", default=None, type=float, help="LD window size in kilobases.")
    parser.add_argument("--ld-wind-cm", default=None, type=float, help="LD window size in centiMorgans.")
    parser.add_argument("--maf", default=None, type=float, help="Optional MAF filter for retained SNPs when MAF is available.")
    parser.add_argument("--chunk-size", default=50, type=int, help="Chunk size for legacy PLINK block computations.")
    parser.add_argument("--per-chr-output", default=False, action="store_true", help="Emit per-chromosome outputs instead of aggregated outputs.")
    parser.add_argument("--yes-really", default=False, action="store_true", help="Allow whole-chromosome LD windows.")
    parser.add_argument("--log-level", default="INFO", choices=("DEBUG", "INFO", "WARNING", "ERROR"), help="Logging verbosity.")
    return parser


def run_ldscore_from_args(args: argparse.Namespace) -> LDScoreResult:
    """Run LD-score calculation from a parsed CLI namespace."""
    normalized_args, common_config = _normalize_run_args(args)
    kernel_ldscore.validate_args(normalized_args)
    regression_snps = _load_regression_snps(normalized_args.regression_snps, common_config)
    baseline_tokens = split_cli_path_tokens(normalized_args.baseline_annot)
    baseline_chr_tokens = split_cli_path_tokens(normalized_args.baseline_annot_chr)
    query_tokens = split_cli_path_tokens(normalized_args.query_annot)
    query_chr_tokens = split_cli_path_tokens(normalized_args.query_annot_chr)
    baseline_files = (
        []
        if not baseline_tokens
        else resolve_file_group(
            baseline_tokens,
            suffixes=ANNOTATION_SUFFIXES,
            label="baseline annotation",
            allow_chromosome_suite=True,
        )
    )
    query_files = (
        []
        if not query_tokens
        else resolve_file_group(
            query_tokens,
            suffixes=ANNOTATION_SUFFIXES,
            label="query annotation",
            allow_chromosome_suite=True,
        )
    )
    r2_tokens = split_cli_path_tokens(normalized_args.r2_table)
    r2_chr_tokens = split_cli_path_tokens(normalized_args.r2_table_chr)
    r2_files = (
        []
        if not r2_tokens
        else resolve_file_group(
            r2_tokens,
            suffixes=PARQUET_SUFFIXES,
            label="parquet R2",
            allow_chromosome_suite=True,
        )
    )
    freq_tokens = split_cli_path_tokens(normalized_args.frqfile)
    freq_chr_tokens = split_cli_path_tokens(normalized_args.frqfile_chr)
    freq_files = (
        []
        if not freq_tokens
        else resolve_file_group(
            freq_tokens,
            suffixes=FREQUENCY_SUFFIXES,
            label="frequency or metadata",
            allow_chromosome_suite=True,
        )
    )
    chromosome_results: list[ChromLDScoreResult] = []
    calculator = LDScoreCalculator()

    for chrom in _chromosome_set_from_annotation_inputs(
        baseline_files=baseline_files,
        query_files=query_files,
        baseline_chr_tokens=baseline_chr_tokens,
        query_chr_tokens=query_chr_tokens,
    ):
        chrom_baseline_files = (
            resolve_chromosome_group(
                baseline_chr_tokens,
                chrom=chrom,
                suffixes=ANNOTATION_SUFFIXES,
                label="baseline annotation",
                required=False,
            )
            if baseline_chr_tokens
            else []
        )
        chrom_baseline_files += baseline_files
        chrom_query_files = (
            resolve_chromosome_group(
                query_chr_tokens,
                chrom=chrom,
                suffixes=ANNOTATION_SUFFIXES,
                label="query annotation",
                required=False,
            )
            if query_chr_tokens
            else []
        )
        chrom_query_files += query_files
        bundle = kernel_ldscore.combine_annotation_groups(
            baseline_files=chrom_baseline_files,
            query_files=chrom_query_files,
            chrom=chrom,
            identifier_mode=normalized_args.snp_identifier,
        )
        if bundle is None:
            continue
        chrom_args = argparse.Namespace(**vars(normalized_args))
        if normalized_args.r2_table or normalized_args.r2_table_chr:
            chrom_r2_files = (
                resolve_chromosome_group(
                    r2_chr_tokens,
                    chrom=chrom,
                    suffixes=PARQUET_SUFFIXES,
                    label="parquet R2",
                    required=False,
                )
                if r2_chr_tokens
                else []
            )
            chrom_r2_files += r2_files
            chrom_freq_files = (
                resolve_chromosome_group(
                    freq_chr_tokens,
                    chrom=chrom,
                    suffixes=FREQUENCY_SUFFIXES,
                    label="frequency or metadata",
                    required=False,
                )
                if freq_chr_tokens
                else []
            )
            chrom_freq_files += freq_files
            chrom_args.r2_table = ",".join(chrom_r2_files) or None
            chrom_args.r2_table_chr = None
            chrom_args.frqfile = ",".join(chrom_freq_files) or None
            chrom_args.frqfile_chr = None
            chrom_args.bfile = None
            chrom_args.bfile_chr = None
            legacy_result = kernel_ldscore.compute_chrom_from_parquet(chrom, bundle, chrom_args, regression_snps)
        else:
            bfile_token = normalized_args.bfile if normalized_args.bfile is not None else normalized_args.bfile_chr
            chrom_args.bfile = None if bfile_token is None else resolve_plink_prefix(bfile_token, chrom=chrom)
            chrom_args.bfile_chr = None
            legacy_result = kernel_ldscore.compute_chrom_from_plink(chrom, bundle, chrom_args, regression_snps)
        chromosome_results.append(
            calculator._wrap_legacy_chrom_result(
                legacy_result,
                common_config=common_config,
                regression_snps=regression_snps,
            )
        )

    if not chromosome_results:
        raise ValueError("No chromosome results were produced.")

    result = calculator._aggregate_chromosome_results(chromosome_results, common_config=common_config)
    summary = calculator.write_outputs(result, _output_spec_from_args(normalized_args))
    return _replace_result_output_paths(result, summary.output_paths)


def run_ldscore(**kwargs) -> LDScoreResult:
    """Convenience wrapper around :func:`run_ldscore_from_args`.

    Keyword arguments are interpreted as CLI-equivalent option names without the
    leading ``--``.
    """
    parser = build_parser()
    defaults = vars(parser.parse_args(["--out", "placeholder"]))
    defaults.update(kwargs)
    args = argparse.Namespace(**defaults)
    return run_ldscore_from_args(args)


def main(argv: Sequence[str] | None = None) -> LDScoreResult:
    """Command-line entry point for the LD-score workflow."""
    parser = build_parser()
    args = parser.parse_args(argv)
    return run_ldscore_from_args(args)


def _normalize_run_args(args: argparse.Namespace) -> tuple[argparse.Namespace, CommonConfig]:
    """Normalize CLI-style args and derive the shared ``CommonConfig`` object."""
    normalized_mode = normalize_snp_identifier_mode(args.snp_identifier)
    normalized_args = argparse.Namespace(**vars(args))
    normalized_args.snp_identifier = normalized_mode
    normalized_args.out = normalize_path_token(args.out)
    normalized_args.keep = normalize_optional_path_token(getattr(args, "keep", None))
    common_config = CommonConfig(
        snp_identifier=normalized_mode,
        genome_build=getattr(args, "genome_build", None),
        log_level=getattr(args, "log_level", "INFO"),
    )
    return normalized_args, common_config


def _load_regression_snps(path: str | None, common_config: CommonConfig) -> set[str] | None:
    """Load the optional global regression SNP universe using the active identifier mode."""
    if not path:
        return None
    return read_global_snp_restriction(resolve_scalar_path(path, label="regression SNP list"), common_config.snp_identifier)


def _output_spec_from_args(args: argparse.Namespace) -> OutputSpec:
    """Translate LD-score CLI arguments into the standard output configuration."""
    out_path = Path(normalize_path_token(args.out))
    return OutputSpec(
        out_prefix=out_path.name,
        output_dir=str(out_path.parent),
        write_per_chrom=bool(getattr(args, "per_chr_output", False)),
        aggregate_across_chromosomes=not bool(getattr(args, "per_chr_output", False)),
        print_snps_path=normalize_optional_path_token(getattr(args, "print_snps", None)),
        write_summary_json=False,
        write_summary_tsv=False,
        write_run_metadata=False,
    )


def _replace_result_output_paths(result: LDScoreResult, output_paths: dict[str, str]) -> LDScoreResult:
    """Return ``result`` with updated artifact-path metadata after writing outputs."""
    return LDScoreResult(
        reference_metadata=result.reference_metadata,
        ld_scores=result.ld_scores,
        regression_metadata=result.regression_metadata,
        w_ld=result.w_ld,
        snp_count_totals=result.snp_count_totals,
        baseline_columns=result.baseline_columns,
        query_columns=result.query_columns,
        reference_snps=result.reference_snps,
        regression_snps=result.regression_snps,
        chromosome_results=result.chromosome_results,
        output_paths=dict(output_paths),
    )


def _sort_paired_frames(metadata_frames: Sequence[pd.DataFrame], value_frames: Sequence[pd.DataFrame]) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Co-sort metadata and value frames after cross-chromosome concatenation."""
    merged_frames = []
    value_columns = list(value_frames[0].columns)
    for metadata, values in zip(metadata_frames, value_frames):
        merged = pd.concat([metadata.reset_index(drop=True), values.reset_index(drop=True)], axis=1)
        merged_frames.append(merged)
    merged_frame = pd.concat(merged_frames, axis=0, ignore_index=True)
    merged_frame = kernel_ldscore.sort_frame_by_genomic_position(merged_frame)
    metadata_columns = [column for column in _LDSCORE_SUFFIX_COLUMNS if column in merged_frame.columns]
    metadata = merged_frame.loc[:, metadata_columns].copy()
    values = merged_frame.loc[:, value_columns].copy()
    return metadata.reset_index(drop=True), values.reset_index(drop=True)


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
    )


def _namespace_from_configs(chrom: str, ref_panel, ldscore_config: LDScoreConfig, common_config: CommonConfig) -> argparse.Namespace:
    """Build the legacy-kernel namespace expected by the chromosome backends."""
    spec = getattr(ref_panel, "spec", None)
    backend = getattr(spec, "backend", None)
    bfile = None
    r2_table = None
    frqfile = None
    if backend == "plink" and getattr(spec, "bfile_prefix", None) is not None:
        bfile = resolve_plink_prefix(spec.bfile_prefix, chrom=chrom)
    if backend == "parquet_r2" and getattr(spec, "r2_table_paths", ()):
        r2_table = ",".join(
            resolve_chromosome_group(
                getattr(spec, "r2_table_paths", ()),
                chrom=chrom,
                suffixes=PARQUET_SUFFIXES,
                label="parquet R2",
                required=False,
            )
        ) or None
    if getattr(spec, "maf_metadata_paths", ()):
        frqfile = ",".join(
            resolve_chromosome_group(
                getattr(spec, "maf_metadata_paths", ()),
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
        snp_identifier=common_config.snp_identifier,
        genome_build=(getattr(spec, "genome_build", None) or common_config.genome_build),
        r2_bias_mode=getattr(spec, "r2_bias_mode", None),
        r2_sample_size=getattr(spec, "sample_size", None),
        regression_snps=None,
        frqfile=frqfile,
        frqfile_chr=None,
        keep=ldscore_config.keep_individuals_path,
        ld_wind_snps=ldscore_config.ld_wind_snps,
        ld_wind_kb=ldscore_config.ld_wind_kb,
        ld_wind_cm=ldscore_config.ld_wind_cm,
        maf=ldscore_config.maf_min,
        chunk_size=ldscore_config.chunk_size,
        per_chr_output=False,
        yes_really=ldscore_config.whole_chromosome_ok,
        log_level=common_config.log_level,
    )


def _chromosome_set_from_annotation_inputs(
    *,
    baseline_files: Sequence[str],
    query_files: Sequence[str],
    baseline_chr_tokens: Sequence[str],
    query_chr_tokens: Sequence[str],
) -> list[str]:
    """Discover the chromosome set implied by workflow-resolved annotation inputs."""
    chromosomes: set[str] = set()
    for path in list(baseline_files) + list(query_files):
        df = kernel_ldscore.read_text_table(path)
        if "CHR" not in df.columns:
            raise ValueError(f"{path} is missing CHR.")
        chromosomes.update(df["CHR"].map(kernel_ldscore.normalize_chromosome).unique().tolist())

    suite_tokens = list(baseline_chr_tokens) + list(query_chr_tokens)
    if suite_tokens:
        for chrom in [str(i) for i in range(1, 23)] + ["X", "Y", "MT", "M"]:
            files = []
            if query_chr_tokens:
                files += resolve_chromosome_group(
                    query_chr_tokens,
                    chrom=chrom,
                    suffixes=ANNOTATION_SUFFIXES,
                    label="query annotation",
                    required=False,
                )
            if baseline_chr_tokens:
                files += resolve_chromosome_group(
                    baseline_chr_tokens,
                    chrom=chrom,
                    suffixes=ANNOTATION_SUFFIXES,
                    label="baseline annotation",
                    required=False,
                )
            if files:
                chromosomes.add(kernel_ldscore.normalize_chromosome(chrom))

    if not chromosomes:
        raise ValueError("No annotation chromosomes could be resolved from the supplied inputs.")
    return sorted(chromosomes, key=kernel_ldscore.chrom_sort_key)
