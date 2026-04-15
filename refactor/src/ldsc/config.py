"""Workflow configuration dataclasses for the refactored package.

Core functionality:
    Define validated, frozen configuration objects for the public workflows.

Overview
--------
This module centralizes the user-facing configuration surface for annotation
building, reference-panel loading, LD-score calculation, summary-statistics
munging, and regression. Each dataclass is intentionally small and validation
focused so that workflow modules receive explicit, typed settings instead of a
single mutable argument namespace.

Design Notes
------------
- Config objects are frozen so the resolved runtime state stays explicit.
- Path normalization and basic value validation live here rather than being
  duplicated in the workflow modules.
- These classes describe the refactored API, not the historical script flags.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Literal


SNPIdentifierMode = Literal["rsid", "chr_pos"]
GenomeBuild = Literal["hg19", "hg38"]
LogLevel = Literal["DEBUG", "INFO", "WARNING", "ERROR"]
RefPanelBackend = Literal["auto", "plink", "parquet_r2"]
CompressionMode = Literal["auto", "gzip", "bz2", "none"]
R2BiasMode = Literal["raw", "unbiased"]


def _normalize_optional_path(path: str | Path | None) -> str | None:
    """Return ``path`` as a string or ``None``.

    This helper keeps the public config dataclasses tolerant of both ``Path``
    and ``str`` inputs while preserving ``None`` for omitted optional fields.
    """
    if path is None:
        return None
    return str(path)


def _normalize_path_tuple(values: tuple[str | Path, ...] | list[str | Path] | None) -> tuple[str, ...]:
    """Convert a sequence of path-like objects to a string tuple."""
    if not values:
        return ()
    return tuple(str(value) for value in values)


def _normalize_log_level(level: str) -> LogLevel:
    """Normalize a logging level string to the supported uppercase literal."""
    normalized = level.upper()
    allowed = {"DEBUG", "INFO", "WARNING", "ERROR"}
    if normalized not in allowed:
        raise ValueError(f"log_level must be one of {sorted(allowed)}; got {level!r}.")
    return normalized  # type: ignore[return-value]


@dataclass(frozen=True)
class CommonConfig:
    """Shared configuration used across the refactored workflows.

    Parameters
    ----------
    snp_identifier : {"rsid", "chr_pos"}, optional
        Global SNP identifier mode. Default is ``"chr_pos"``. ``"rsid"``
        expects an explicit SNP column, while ``"chr_pos"`` builds identifiers
        from chromosome and base-pair position.
    genome_build : {"hg19", "hg38"} or None, optional
        Genome-build context for ``chr_pos`` workflows. Default is ``None``.
    global_snp_restriction_path : str or None, optional
        Optional path to a SNP list or table that restricts the SNP universe used
        by annotation, reference-panel, and regression workflows. Default is
        ``None``.
    log_level : {"DEBUG", "INFO", "WARNING", "ERROR"}, optional
        Requested logging verbosity for workflow modules. Default is ``"INFO"``.
    fail_on_missing_metadata : bool, optional
        If ``True``, treat missing optional metadata as a hard error rather than
        tolerating partial metadata tables. Default is ``False``.
    """
    snp_identifier: SNPIdentifierMode = "chr_pos"
    genome_build: GenomeBuild | None = None
    global_snp_restriction_path: str | None = None
    log_level: LogLevel = "INFO"
    fail_on_missing_metadata: bool = False

    def __post_init__(self) -> None:
        if self.snp_identifier not in {"rsid", "chr_pos"}:
            raise ValueError("snp_identifier must be 'rsid' or 'chr_pos'.")
        if self.genome_build not in {None, "hg19", "hg38"}:
            raise ValueError("genome_build must be None, 'hg19', or 'hg38'.")
        object.__setattr__(self, "log_level", _normalize_log_level(self.log_level))
        object.__setattr__(
            self, "global_snp_restriction_path", _normalize_optional_path(self.global_snp_restriction_path)
        )


@dataclass(frozen=True)
class AnnotationBuildConfig:
    """Configuration for SNP-level annotation assembly and BED projection.

    Parameters
    ----------
    baseline_annotation_paths : tuple of str, optional
        Baseline annotation files to align and combine. Default is ``()``.
    query_annotation_paths : tuple of str, optional
        Query annotation files to align with the baseline. Default is ``()``.
    query_bed_paths : tuple of str, optional
        BED inputs that should be projected to SNP-level annotations. Default is
        ``()``.
    out_prefix : str or None, optional
        Output prefix used by file-writing helpers. Default is ``None``.
    batch_mode : bool, optional
        If ``True``, write one output directory per BED input in the batch
        projection workflow. Default is ``True``.
    compression : {"auto", "gzip", "bz2", "none"}, optional
        Output compression preference for generated annotation files. Default is
        ``"gzip"``.
    """
    baseline_annotation_paths: tuple[str, ...] = field(default_factory=tuple)
    query_annotation_paths: tuple[str, ...] = field(default_factory=tuple)
    query_bed_paths: tuple[str, ...] = field(default_factory=tuple)
    out_prefix: str | None = None
    batch_mode: bool = True
    compression: CompressionMode = "gzip"

    def __post_init__(self) -> None:
        object.__setattr__(self, "baseline_annotation_paths", _normalize_path_tuple(self.baseline_annotation_paths))
        object.__setattr__(self, "query_annotation_paths", _normalize_path_tuple(self.query_annotation_paths))
        object.__setattr__(self, "query_bed_paths", _normalize_path_tuple(self.query_bed_paths))
        object.__setattr__(self, "out_prefix", _normalize_optional_path(self.out_prefix))
        if self.compression not in {"auto", "gzip", "bz2", "none"}:
            raise ValueError("compression must be one of 'auto', 'gzip', 'bz2', or 'none'.")


@dataclass(frozen=True)
class RefPanelConfig:
    """Configuration for choosing and parameterizing a reference-panel backend.

    Parameters
    ----------
    backend : {"auto", "plink", "parquet_r2"}, optional
        Backend family used to supply LD reference information. Default is
        ``"auto"``.
    plink_prefix, plink_prefix_chr : str or None, optional
        PLINK ``.bed/.bim/.fam`` prefix or chromosome-pattern prefix. Default is
        ``None`` for both fields.
    parquet_r2_paths, parquet_r2_paths_chr : tuple of str, optional
        Sorted parquet R2 tables or chromosome-pattern prefixes. Default is
        ``()`` for both fields.
    frequency_paths, frequency_paths_chr : tuple of str, optional
        Sidecar frequency or metadata files aligned to the reference panel.
        Default is ``()`` for both fields.
    r2_bias_mode : {"raw", "unbiased"} or None, optional
        Whether parquet R2 values need sample-size correction. Default is
        ``None``.
    r2_sample_size : float or None, optional
        Sample size used when correcting raw parquet R2 values. Default is
        ``None``.
    """
    backend: RefPanelBackend = "auto"
    plink_prefix: str | None = None
    plink_prefix_chr: str | None = None
    parquet_r2_paths: tuple[str, ...] = field(default_factory=tuple)
    parquet_r2_paths_chr: tuple[str, ...] = field(default_factory=tuple)
    frequency_paths: tuple[str, ...] = field(default_factory=tuple)
    frequency_paths_chr: tuple[str, ...] = field(default_factory=tuple)
    r2_bias_mode: R2BiasMode | None = None
    r2_sample_size: float | None = None

    def __post_init__(self) -> None:
        if self.backend not in {"auto", "plink", "parquet_r2"}:
            raise ValueError("backend must be 'auto', 'plink', or 'parquet_r2'.")
        if self.r2_bias_mode not in {None, "raw", "unbiased"}:
            raise ValueError("r2_bias_mode must be None, 'raw', or 'unbiased'.")
        if self.r2_sample_size is not None and self.r2_sample_size <= 0:
            raise ValueError("r2_sample_size must be positive when provided.")
        object.__setattr__(self, "plink_prefix", _normalize_optional_path(self.plink_prefix))
        object.__setattr__(self, "plink_prefix_chr", _normalize_optional_path(self.plink_prefix_chr))
        object.__setattr__(self, "parquet_r2_paths", _normalize_path_tuple(self.parquet_r2_paths))
        object.__setattr__(self, "parquet_r2_paths_chr", _normalize_path_tuple(self.parquet_r2_paths_chr))
        object.__setattr__(self, "frequency_paths", _normalize_path_tuple(self.frequency_paths))
        object.__setattr__(self, "frequency_paths_chr", _normalize_path_tuple(self.frequency_paths_chr))


@dataclass(frozen=True)
class LDScoreConfig:
    """Configuration for chromosome-wise LD-score calculation.

    Exactly one LD-window field must be provided. The remaining fields control
    retained-SNP filtering and count-artifact emission.

    Parameters
    ----------
    ld_wind_snps : int or None, optional
        Window size measured in SNP count. Default is ``None``.
    ld_wind_kb : float or None, optional
        Window size measured in kilobases. Default is ``None``.
    ld_wind_cm : float or None, optional
        Window size measured in centiMorgans. Default is ``None``.
    maf_min : float or None, optional
        Minimum retained minor-allele frequency, if MAF metadata are available.
        Default is ``None``.
    chunk_size : int, optional
        Chunk size for legacy PLINK block computations. Default is ``50``.
    compute_m5_50 : bool, optional
        If ``True``, emit the common-SNP count vector used by LDSC's
        ``.M_5_50`` convention. Default is ``True``.
    whole_chromosome_ok : bool, optional
        Override the guard that rejects windows effectively spanning an entire
        chromosome. Default is ``False``.
    """
    ld_wind_snps: int | None = None
    ld_wind_kb: float | None = None
    ld_wind_cm: float | None = None
    maf_min: float | None = None
    chunk_size: int = 50
    compute_m5_50: bool = True
    whole_chromosome_ok: bool = False

    def __post_init__(self) -> None:
        windows = [self.ld_wind_snps, self.ld_wind_kb, self.ld_wind_cm]
        if sum(value is not None for value in windows) != 1:
            raise ValueError("Exactly one LD-window option must be set.")
        if self.ld_wind_snps is not None and self.ld_wind_snps <= 0:
            raise ValueError("ld_wind_snps must be positive.")
        if self.ld_wind_kb is not None and self.ld_wind_kb <= 0:
            raise ValueError("ld_wind_kb must be positive.")
        if self.ld_wind_cm is not None and self.ld_wind_cm <= 0:
            raise ValueError("ld_wind_cm must be positive.")
        if self.maf_min is not None and not 0 <= self.maf_min <= 0.5:
            raise ValueError("maf_min must lie in [0, 0.5].")
        if self.chunk_size <= 0:
            raise ValueError("chunk_size must be positive.")


@dataclass(frozen=True)
class MungeConfig:
    """Configuration for legacy-compatible summary-statistics munging.

    This dataclass preserves the established LDSC munging behavior while
    exposing the options through an explicit Python object rather than a script
    namespace.
    """
    out_prefix: str
    N: float | None = None
    N_cas: float | None = None
    N_con: float | None = None
    info_min: float = 0.9
    maf_min: float = 0.01
    n_min: float | None = None
    nstudy_min: float | None = None
    chunk_size: int = int(5e6)
    merge_alleles_path: str | None = None
    signed_sumstats_spec: str | None = None
    ignore_columns: tuple[str, ...] = field(default_factory=tuple)
    no_alleles: bool = False
    a1_inc: bool = False
    keep_maf: bool = False
    daner: bool = False
    daner_n: bool = False

    def __post_init__(self) -> None:
        if self.info_min < 0:
            raise ValueError("info_min must be non-negative.")
        if self.maf_min < 0 or self.maf_min > 0.5:
            raise ValueError("maf_min must lie in [0, 0.5].")
        if self.chunk_size <= 0:
            raise ValueError("chunk_size must be positive.")
        object.__setattr__(self, "out_prefix", str(self.out_prefix))
        object.__setattr__(self, "merge_alleles_path", _normalize_optional_path(self.merge_alleles_path))
        object.__setattr__(self, "ignore_columns", tuple(self.ignore_columns))


@dataclass(frozen=True)
class RegressionConfig:
    """Configuration for heritability and genetic-correlation regressions.

    Parameters
    ----------
    n_blocks : int, optional
        Requested number of block-jackknife partitions. Default is ``200``.
    use_m_5_50 : bool, optional
        If ``True``, prefer the common-SNP count vector when it is available.
        Default is ``True``.
    use_intercept : bool, optional
        If ``False``, constrain the intercept to the LDSC default for the model
        being fit. Default is ``True``.
    intercept_h2, intercept_gencov : float, list of float, or None, optional
        Fixed intercept values for single-trait and cross-trait models.
        Defaults are ``None``.
    two_step_cutoff : float or None, optional
        Threshold for the two-step estimator used by the regression kernel.
        Default is ``None``.
    chisq_max : float or None, optional
        Maximum allowed chi-square statistic before a row is filtered. Default is
        ``None``.
    samp_prev, pop_prev : float, list of float, or None, optional
        Liability-scale prevalence inputs. Defaults are ``None``.
    """
    n_blocks: int = 200
    use_m_5_50: bool = True
    use_intercept: bool = True
    intercept_h2: float | list[float] | None = None
    intercept_gencov: float | list[float] | None = None
    two_step_cutoff: float | None = None
    chisq_max: float | None = None
    samp_prev: float | list[float] | None = None
    pop_prev: float | list[float] | None = None

    def __post_init__(self) -> None:
        if self.n_blocks <= 1:
            raise ValueError("n_blocks must be greater than 1.")
        if self.two_step_cutoff is not None and self.two_step_cutoff <= 0:
            raise ValueError("two_step_cutoff must be positive when provided.")
        if self.chisq_max is not None and self.chisq_max <= 0:
            raise ValueError("chisq_max must be positive when provided.")
