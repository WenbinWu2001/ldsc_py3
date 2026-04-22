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
from pathlib import Path
from typing import Any, Sequence
import warnings

import numpy as np
import pandas as pd

from .column_inference import normalize_snp_identifier_mode
from .config import (
    ConfigMismatchError,
    GlobalConfig,
    LDScoreConfig,
    get_global_config,
    print_global_config_banner,
    validate_config_compatibility,
)
from .genome_build_inference import validate_auto_genome_build_mode
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
    """Normalized per-chromosome LD-score output matching the written file shape.

    Each row represents one retained regression SNP on a single chromosome. The
    table already includes the embedded ``regr_weight`` column used by the
    regression workflows; no separate weight table is required.
    """
    chrom: str
    ldscore_table: pd.DataFrame
    snp_count_totals: dict[str, np.ndarray]
    baseline_columns: list[str]
    query_columns: list[str]
    ld_reference_snps: frozenset[str]
    ld_regression_snps: frozenset[str]
    output_paths: dict[str, str] = field(default_factory=dict)
    config_snapshot: GlobalConfig | None = None

    def validate(self) -> None:
        """Check the normalized public contract for chromosome-level results."""
        required = {"CHR", "SNP", "BP", "regr_weight"}
        missing = required - set(self.ldscore_table.columns)
        if missing:
            raise ValueError(f"ldscore_table is missing required columns: {sorted(missing)}")

    def summary(self) -> dict[str, Any]:
        """Return a compact summary of chromosome-level retained rows."""
        return {
            "chrom": self.chrom,
            "n_rows": len(self.ldscore_table),
            "count_keys": sorted(self.snp_count_totals.keys()),
        }


@dataclass(frozen=True)
class LDScoreResult:
    """Aggregated cross-chromosome LD-score result in normalized row-table form.

    The public result keeps one merged ``ldscore_table`` whose rows correspond
    to the retained regression SNP set across all processed chromosomes. The
    full LD-computation reference universe is intentionally not reconstructed
    from these normalized rows and is therefore exposed as
    ``ld_reference_snps = frozenset()``.
    """
    ldscore_table: pd.DataFrame
    snp_count_totals: dict[str, np.ndarray]
    baseline_columns: list[str]
    query_columns: list[str]
    ld_reference_snps: frozenset[str]
    ld_regression_snps: frozenset[str]
    chromosome_results: list[ChromLDScoreResult]
    output_paths: dict[str, str] = field(default_factory=dict)
    config_snapshot: GlobalConfig | None = None

    def validate(self) -> None:
        """Check the normalized public contract for aggregated results."""
        required = {"CHR", "SNP", "BP", "regr_weight"}
        missing = required - set(self.ldscore_table.columns)
        if missing:
            raise ValueError(f"ldscore_table is missing required columns: {sorted(missing)}")

    def summary(self) -> dict[str, Any]:
        """Return a compact cross-chromosome summary."""
        return {
            "n_rows": len(self.ldscore_table),
            "chromosomes": [result.chrom for result in self.chromosome_results],
            "count_keys": sorted(self.snp_count_totals.keys()),
        }


class LDScoreCalculator:
    """Orchestrate chromosome-wise LD-score calculation.

    This service assembles annotation and reference-panel inputs, delegates the
    heavy computation to the internal LD-score kernel, aggregates chromosome
    outputs, and optionally hands the result to the output layer. The calculator
    treats ``ref_panel`` as the owner of the reference-panel SNP universe:
    chromosome bundles are intersected with ``ref_panel.load_metadata(chrom)``
    before the kernel runs so ``RefPanelSpec.ref_panel_snps_path`` is honored
    without leaking that setting into the calculator interface.
    """

    def __init__(self, output_manager: OutputManager | None = None) -> None:
        """Initialize the calculator with the output manager used for writes."""
        self.output_manager = output_manager or OutputManager()

    def run(
        self,
        annotation_bundle,
        ref_panel,
        ldscore_config: LDScoreConfig,
        global_config: GlobalConfig,
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
            metadata. Any ``RefPanelSpec.ref_panel_snps_path`` restriction is
            already applied when the workflow calls ``ref_panel.load_metadata``.
        ldscore_config : LDScoreConfig
            LD-window and retained-SNP settings.
        global_config : GlobalConfig
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
                    f"genome_build mismatch between RefPanelSpec ({ref_panel_build!r}) "
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
        if backend == "parquet_r2" and ldscore_config.keep_individuals_path is not None:
            raise ValueError("keep_individuals_path/--keep is only supported for PLINK reference panels.")
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
        result = ChromLDScoreResult(
            chrom=str(legacy_result.chrom),
            ldscore_table=ldscore_table,
            snp_count_totals=count_map,
            baseline_columns=list(legacy_result.baseline_columns),
            query_columns=list(legacy_result.query_columns),
            ld_reference_snps=frozenset(),
            ld_regression_snps=frozenset(build_snp_id_series(ldscore_table, global_config.snp_identifier)),
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
        ldscore_table = pd.concat(
            [result.ldscore_table for result in chromosome_results],
            axis=0,
            ignore_index=True,
        )
        ldscore_table = kernel_ldscore.sort_frame_by_genomic_position(ldscore_table)
        result = LDScoreResult(
            ldscore_table=ldscore_table,
            snp_count_totals=count_totals,
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
        output_spec: OutputSpec,
        config_snapshot: dict[str, Any] | None = None,
    ):
        """Write a previously computed result through the output layer."""
        return self.output_manager.write_outputs(result, output_spec, config_snapshot=config_snapshot)


def build_parser() -> argparse.ArgumentParser:
    """Build the feature parser for LD-score calculation."""
    parser = argparse.ArgumentParser(
        description="Estimate LDSC-compatible LD scores from SNP-level annotation files using PLINK or sorted parquet R2 input.",
        allow_abbrev=False,
    )
    parser.add_argument("--out", required=True, help="Output prefix.")
    query_group = parser.add_mutually_exclusive_group()
    query_group.add_argument(
        "--query-annot",
        default=None,
        help="Comma-separated query annotation path tokens: exact paths, globs, or explicit @ suite tokens.",
    )
    query_group.add_argument(
        "--query-annot-bed",
        default=None,
        help="Comma-separated BED file path tokens projected in memory as query annotations.",
    )
    parser.add_argument("--baseline-annot", default=None, help="Comma-separated baseline annotation path tokens: exact paths, globs, or explicit @ suite tokens.")
    parser.add_argument("--bfile", default=None, help="PLINK prefix token for the reference panel.")
    parser.add_argument("--r2-table", default=None, help="Comma-separated parquet R2 path tokens: exact paths, globs, or explicit @ suite tokens.")
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
    parser.add_argument("--frqfile", default=None, help="Optional frequency or metadata path tokens for MAF and CM.")
    parser.add_argument(
        "--keep",
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
    the normalized public ``LDScoreResult`` with a merged ``ldscore_table``.
    """
    from ._kernel.annotation import AnnotationBuilder, AnnotationSourceSpec

    normalized_args, global_config = _normalize_run_args(args)
    print_global_config_banner("run_ldscore_from_args", global_config)
    kernel_ldscore.validate_args(normalized_args)
    ldscore_config = _ldscore_config_from_args(normalized_args)
    regression_snps = _load_regression_snps(ldscore_config.regression_snps_path, global_config)
    source_spec = AnnotationSourceSpec(
        baseline_annot_paths=tuple(split_cli_path_tokens(normalized_args.baseline_annot)),
        query_annot_paths=tuple(split_cli_path_tokens(normalized_args.query_annot)),
        bed_paths=tuple(split_cli_path_tokens(getattr(normalized_args, "query_annot_bed", None))),
    )
    annotation_bundle = AnnotationBuilder(global_config).run(source_spec)
    ref_panel = _ref_panel_from_args(normalized_args, global_config)
    calculator = LDScoreCalculator()
    return calculator.run(
        annotation_bundle=annotation_bundle,
        ref_panel=ref_panel,
        ldscore_config=ldscore_config,
        global_config=global_config,
        output_spec=_output_spec_from_args(normalized_args),
        regression_snps=regression_snps,
    )


def run_ldscore(**kwargs) -> LDScoreResult:
    """Convenience wrapper around :func:`run_ldscore_from_args`.

    Keyword arguments are interpreted as CLI-equivalent option names without the
    leading ``--``. Shared runtime assumptions such as ``snp_identifier`` and
    ``genome_build`` must be supplied through ``set_global_config(...)`` first,
    while per-run controls such as ``ref_panel_snps_path`` and
    ``regression_snps_path`` remain ordinary keyword arguments here.
    """
    forbidden = sorted({"snp_identifier", "genome_build", "log_level"} & set(kwargs))
    if forbidden:
        joined = ", ".join(forbidden)
        raise ValueError(f"Python run_ldscore() no longer accepts {joined}; call set_global_config(...) first.")
    parser = build_parser()
    defaults = vars(parser.parse_args(["--out", "placeholder"]))
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
    for attr in ("query_annot_chr", "baseline_annot_chr", "bfile_chr", "r2_table_chr", "frqfile_chr", "query_annot_bed"):
        if not hasattr(normalized_args, attr):
            setattr(normalized_args, attr, None)
    for attr in ("ref_panel_snps_path", "regression_snps_path"):
        if not hasattr(normalized_args, attr):
            setattr(normalized_args, attr, None)
    normalized_args.snp_identifier = normalized_mode
    normalized_args.out = normalize_path_token(args.out)
    normalized_args.keep = normalize_optional_path_token(getattr(args, "keep", None))
    normalized_args.ref_panel_snps_path = normalize_optional_path_token(getattr(args, "ref_panel_snps_path", None))
    normalized_args.regression_snps_path = normalize_optional_path_token(getattr(args, "regression_snps_path", None))
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
    from ._kernel.ref_panel import RefPanelLoader, RefPanelSpec

    freq_tokens = split_cli_path_tokens(getattr(args, "frqfile", None))
    ref_panel_snps_path = normalize_optional_path_token(getattr(args, "ref_panel_snps_path", None))
    if getattr(args, "r2_table", None) is not None:
        r2_tokens = split_cli_path_tokens(args.r2_table)
        spec = RefPanelSpec(
            backend="parquet_r2",
            r2_table_paths=tuple(r2_tokens),
            maf_metadata_paths=tuple(freq_tokens),
            r2_bias_mode=getattr(args, "r2_bias_mode", None),
            sample_size=getattr(args, "r2_sample_size", None),
            genome_build=global_config.genome_build,
            ref_panel_snps_path=ref_panel_snps_path,
        )
    else:
        spec = RefPanelSpec(
            backend="plink",
            bfile_prefix=getattr(args, "bfile", None),
            maf_metadata_paths=tuple(freq_tokens),
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
        keep_individuals_path=getattr(args, "keep", None),
        regression_snps_path=getattr(args, "regression_snps_path", None),
        chunk_size=getattr(args, "chunk_size", 50),
        whole_chromosome_ok=getattr(args, "yes_really", False),
    )


def _output_spec_from_args(args: argparse.Namespace) -> OutputSpec:
    """Translate LD-score CLI arguments into the standard output configuration."""
    out_path = Path(normalize_path_token(args.out))
    return OutputSpec(
        out_prefix=out_path.name,
        output_dir=str(out_path.parent),
        enabled_artifacts=["ldscore", "counts"],
        write_annotation_manifest=False,
        write_summary_json=False,
        write_summary_tsv=False,
        write_run_metadata=False,
    )


def _replace_result_output_paths(result: LDScoreResult, output_paths: dict[str, str]) -> LDScoreResult:
    """Return ``result`` with updated artifact-path metadata after writing outputs."""
    return LDScoreResult(
        ldscore_table=result.ldscore_table,
        snp_count_totals=result.snp_count_totals,
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
    ``RefPanelSpec.ref_panel_snps_path``. This helper materializes the intended
    compute-time universe ``B_chrom ∩ A'_chrom`` immediately before the legacy
    kernel call.
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
        snp_identifier=global_config.snp_identifier,
        genome_build=(getattr(spec, "genome_build", None) or global_config.genome_build),
        r2_bias_mode=getattr(spec, "r2_bias_mode", None),
        r2_sample_size=getattr(spec, "sample_size", None),
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
        log_level=global_config.log_level,
    )
