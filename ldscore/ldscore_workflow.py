"""Workflow-layer LD-score orchestration and result objects."""

from __future__ import annotations

import argparse
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Sequence

import numpy as np
import pandas as pd

from .config import CommonConfig, LDScoreConfig
from .identifiers import build_snp_id_series, normalize_snp_identifier_mode, read_global_snp_restriction
from .output import OutputManager, OutputSpec

try:
    from . import ldscore_new
except ImportError:  # pragma: no cover - package import fallback
    import ldscore.ldscore_new as ldscore_new


_LDSCORE_SUFFIX_COLUMNS = ("CHR", "SNP", "BP", "CM", "MAF")


@dataclass(frozen=True)
class _LegacyChromResult:
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
        if len(self.reference_metadata) != len(self.ld_scores):
            raise ValueError("reference_metadata and ld_scores must have the same number of rows.")
        if len(self.regression_metadata) != len(self.w_ld):
            raise ValueError("regression_metadata and w_ld must have the same number of rows.")

    def to_ldscore_table(self) -> pd.DataFrame:
        return pd.concat([self.reference_metadata.reset_index(drop=True), self.ld_scores.reset_index(drop=True)], axis=1)

    def to_weight_table(self) -> pd.DataFrame:
        return pd.concat([self.regression_metadata.reset_index(drop=True), self.w_ld.reset_index(drop=True)], axis=1)

    def summary(self) -> dict[str, Any]:
        return {
            "chrom": self.chrom,
            "n_reference_rows": len(self.reference_metadata),
            "n_regression_rows": len(self.regression_metadata),
            "count_keys": sorted(self.snp_count_totals.keys()),
        }


@dataclass(frozen=True)
class LDScoreResult:
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
        if len(self.reference_metadata) != len(self.ld_scores):
            raise ValueError("reference_metadata and ld_scores must have the same number of rows.")
        if len(self.regression_metadata) != len(self.w_ld):
            raise ValueError("regression_metadata and w_ld must have the same number of rows.")

    def to_ldscore_table(self) -> pd.DataFrame:
        return pd.concat([self.reference_metadata.reset_index(drop=True), self.ld_scores.reset_index(drop=True)], axis=1)

    def to_weight_table(self) -> pd.DataFrame:
        return pd.concat([self.regression_metadata.reset_index(drop=True), self.w_ld.reset_index(drop=True)], axis=1)

    def summary(self) -> dict[str, Any]:
        return {
            "n_reference_rows": len(self.reference_metadata),
            "n_regression_rows": len(self.regression_metadata),
            "chromosomes": [result.chrom for result in self.chromosome_results],
            "count_keys": sorted(self.snp_count_totals.keys()),
        }


class LDScoreCalculator:
    """Orchestrate chromosome-wise LD-score calculation around the legacy kernels."""

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
        args = _namespace_from_configs(
            chrom=chrom,
            ref_panel=ref_panel,
            ldscore_config=ldscore_config,
            common_config=common_config,
        )
        legacy_bundle = ldscore_new.AnnotationBundle(
            metadata=annotation_bundle.metadata.copy(),
            annotations=annotation_bundle.annotation_matrix(include_query=True).copy(),
            baseline_columns=list(annotation_bundle.baseline_columns),
            query_columns=list(annotation_bundle.query_columns),
        )
        if getattr(ref_panel.spec, "backend", None) == "parquet_r2":
            legacy_result = ldscore_new.compute_chrom_from_parquet(chrom, legacy_bundle, args, regression_snps)
        else:
            legacy_result = ldscore_new.compute_chrom_from_plink(chrom, legacy_bundle, args, regression_snps)
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
        return self.output_manager.write_outputs(result, output_spec, config_snapshot=config_snapshot)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Estimate LDSC-compatible LD scores from SNP-level annotation files using PLINK or sorted parquet R2 input.",
    )
    parser.add_argument("--out", required=True, help="Output prefix.")
    parser.add_argument("--query-annot", default=None, help="Comma-separated SNP-level query annotation files or prefixes.")
    parser.add_argument("--query-annot-chr", default=None, help="Comma-separated chromosome-pattern prefixes for query annotation files.")
    parser.add_argument("--baseline-annot", default=None, help="Comma-separated SNP-level baseline annotation files or prefixes.")
    parser.add_argument("--baseline-annot-chr", default=None, help="Comma-separated chromosome-pattern prefixes for baseline annotation files.")
    parser.add_argument("--bfile", default=None, help="PLINK prefix for the reference panel.")
    parser.add_argument("--bfile-chr", default=None, help="Chromosome-pattern PLINK prefix for the reference panel.")
    parser.add_argument("--r2-table", default=None, help="Comma-separated sorted parquet R2 files or prefixes.")
    parser.add_argument("--r2-table-chr", default=None, help="Comma-separated chromosome-pattern prefixes for sorted parquet R2 files.")
    parser.add_argument("--snp-identifier", default="chr_pos", help="Identifier mode used to match annotations to the reference panel.")
    parser.add_argument("--genome-build", choices=("hg19", "hg38"), default=None, help="Genome build for chr_pos matching.")
    parser.add_argument("--r2-bias-mode", choices=("raw", "unbiased"), default=None, help="Whether parquet R2 values are raw or already unbiased.")
    parser.add_argument("--r2-sample-size", default=None, type=float, help="LD reference sample size used to correct raw parquet R2 values.")
    parser.add_argument("--regression-snps", default=None, help="Optional SNP list defining the regression SNP set for w_ld.")
    parser.add_argument("--frqfile", default=None, help="Optional frequency or metadata file for MAF and CM.")
    parser.add_argument("--frqfile-chr", default=None, help="Chromosome-pattern frequency or metadata file prefix.")
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
    normalized_args, common_config = _normalize_run_args(args)
    ldscore_new.validate_args(normalized_args)
    regression_snps = _load_regression_snps(normalized_args.regression_snps, common_config)
    chromosome_results: list[ChromLDScoreResult] = []
    calculator = LDScoreCalculator()

    for chrom in ldscore_new.chromosome_set_from_annotation_inputs(normalized_args):
        baseline_files = ldscore_new.resolve_optional_chr_files(
            normalized_args.baseline_annot_chr,
            chrom,
            ("", ".annot", ".annot.gz", ".txt", ".txt.gz", ".tsv", ".tsv.gz"),
        )
        baseline_files += ldscore_new.resolve_annotation_files(normalized_args.baseline_annot)
        query_files = ldscore_new.resolve_optional_chr_files(
            normalized_args.query_annot_chr,
            chrom,
            ("", ".annot", ".annot.gz", ".txt", ".txt.gz", ".tsv", ".tsv.gz"),
        )
        query_files += ldscore_new.resolve_annotation_files(normalized_args.query_annot)
        bundle = ldscore_new.combine_annotation_groups(
            baseline_files=baseline_files,
            query_files=query_files,
            chrom=chrom,
            identifier_mode=normalized_args.snp_identifier,
        )
        if bundle is None:
            continue
        if normalized_args.r2_table or normalized_args.r2_table_chr:
            legacy_result = ldscore_new.compute_chrom_from_parquet(chrom, bundle, normalized_args, regression_snps)
        else:
            legacy_result = ldscore_new.compute_chrom_from_plink(chrom, bundle, normalized_args, regression_snps)
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
    parser = build_parser()
    defaults = vars(parser.parse_args(["--out", "placeholder"]))
    defaults.update(kwargs)
    args = argparse.Namespace(**defaults)
    return run_ldscore_from_args(args)


def main(argv: Sequence[str] | None = None) -> LDScoreResult:
    parser = build_parser()
    args = parser.parse_args(argv)
    return run_ldscore_from_args(args)


def _normalize_run_args(args: argparse.Namespace) -> tuple[argparse.Namespace, CommonConfig]:
    normalized_mode = normalize_snp_identifier_mode(args.snp_identifier)
    legacy_mode = "rsID" if normalized_mode == "rsid" else "chr_pos"
    normalized_args = argparse.Namespace(**vars(args))
    normalized_args.snp_identifier = legacy_mode
    common_config = CommonConfig(
        snp_identifier=normalized_mode,
        genome_build=getattr(args, "genome_build", None),
        log_level=getattr(args, "log_level", "INFO"),
    )
    return normalized_args, common_config


def _load_regression_snps(path: str | None, common_config: CommonConfig) -> set[str] | None:
    if not path:
        return None
    return read_global_snp_restriction(path, common_config.snp_identifier)


def _output_spec_from_args(args: argparse.Namespace) -> OutputSpec:
    out_path = Path(args.out)
    return OutputSpec(
        out_prefix=out_path.name,
        output_dir=str(out_path.parent),
        write_per_chrom=bool(getattr(args, "per_chr_output", False)),
        aggregate_across_chromosomes=not bool(getattr(args, "per_chr_output", False)),
        write_summary_json=False,
        write_summary_tsv=False,
        write_run_metadata=False,
    )


def _replace_result_output_paths(result: LDScoreResult, output_paths: dict[str, str]) -> LDScoreResult:
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
    merged_frames = []
    value_columns = list(value_frames[0].columns)
    for metadata, values in zip(metadata_frames, value_frames):
        merged = pd.concat([metadata.reset_index(drop=True), values.reset_index(drop=True)], axis=1)
        merged_frames.append(merged)
    merged_frame = pd.concat(merged_frames, axis=0, ignore_index=True)
    merged_frame = ldscore_new.sort_frame_by_genomic_position(merged_frame)
    metadata_columns = [column for column in _LDSCORE_SUFFIX_COLUMNS if column in merged_frame.columns]
    metadata = merged_frame.loc[:, metadata_columns].copy()
    values = merged_frame.loc[:, value_columns].copy()
    return metadata.reset_index(drop=True), values.reset_index(drop=True)


def _chromosomes_from_bundle(annotation_bundle) -> list[str]:
    chromosomes = getattr(annotation_bundle, "chromosomes", None)
    if chromosomes:
        return list(chromosomes)
    return sorted(annotation_bundle.metadata["CHR"].astype(str).unique().tolist())


def _slice_annotation_bundle(annotation_bundle, chrom: str):
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
    spec = getattr(ref_panel, "spec", None)
    backend = getattr(spec, "backend", None)
    return argparse.Namespace(
        out=None,
        query_annot=None,
        query_annot_chr=None,
        baseline_annot=None,
        baseline_annot_chr=None,
        bfile=getattr(spec, "bfile_prefix", None) if backend == "plink" else None,
        bfile_chr=None,
        r2_table=",".join(getattr(spec, "r2_table_paths", ())) if backend == "parquet_r2" else None,
        r2_table_chr=None,
        snp_identifier="rsID" if common_config.snp_identifier == "rsid" else "chr_pos",
        genome_build=getattr(spec, "genome_build", common_config.genome_build),
        r2_bias_mode=getattr(spec, "r2_bias_mode", None),
        r2_sample_size=getattr(spec, "sample_size", None),
        regression_snps=None,
        frqfile=",".join(getattr(spec, "maf_metadata_paths", ())) or None,
        frqfile_chr=None,
        ld_wind_snps=ldscore_config.ld_wind_snps,
        ld_wind_kb=ldscore_config.ld_wind_kb,
        ld_wind_cm=ldscore_config.ld_wind_cm,
        maf=ldscore_config.maf_min,
        chunk_size=ldscore_config.chunk_size,
        per_chr_output=False,
        yes_really=ldscore_config.whole_chromosome_ok,
        log_level=common_config.log_level,
    )
