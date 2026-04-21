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

from contextlib import contextmanager
from contextvars import ContextVar
from dataclasses import dataclass, field
from os import PathLike
from typing import Literal
import warnings

from .column_inference import normalize_genome_build
from .path_resolution import normalize_optional_path_token, normalize_path_token, normalize_path_tokens


SNPIdentifierMode = Literal["rsid", "chr_pos"]
GenomeBuild = Literal["hg19", "hg38", "auto"]
GenomeBuildInput = Literal["auto", "hg19", "hg37", "hg38", "GRCh37", "GRCh38"]
LogLevel = Literal["DEBUG", "INFO", "WARNING", "ERROR"]
RefPanelBackend = Literal["auto", "plink", "parquet_r2"]
CompressionMode = Literal["auto", "gzip", "bz2", "none"]
R2BiasMode = Literal["raw", "unbiased"]


class ConfigMismatchError(ValueError):
    """Raised when two result objects carry incompatible config snapshots."""


def _normalize_optional_path(path: str | PathLike[str] | None) -> str | None:
    """Return ``path`` as a string or ``None``.

    This helper keeps the public config dataclasses tolerant of both ``Path``
    and ``str`` inputs while preserving ``None`` for omitted optional fields.
    """
    return normalize_optional_path_token(path)


def _normalize_required_path(path: str | PathLike[str]) -> str:
    """Return a required path-like value as a plain string."""
    if path is None or not str(path).strip():
        raise ValueError("Required path value is missing.")
    return normalize_path_token(path)


def _normalize_path_tuple(values) -> tuple[str, ...]:
    """Convert one or many path-like objects to a string tuple of tokens."""
    return normalize_path_tokens(values)


def _normalize_log_level(level: str) -> LogLevel:
    """Normalize a logging level string to the supported uppercase literal."""
    normalized = level.upper()
    allowed = {"DEBUG", "INFO", "WARNING", "ERROR"}
    if normalized not in allowed:
        raise ValueError(f"log_level must be one of {sorted(allowed)}; got {level!r}.")
    return normalized  # type: ignore[return-value]


@dataclass(frozen=True)
class GlobalConfig:
    """Shared configuration used across the refactored workflows.

    Parameters
    ----------
    snp_identifier : {"rsid", "chr_pos"}, optional
        Global SNP identifier mode. Default is ``"chr_pos"``. ``"rsid"``
        expects an explicit SNP column, while ``"chr_pos"`` builds identifiers
        from chromosome and base-pair position.
    genome_build : {"auto", "hg19", "hg37", "GRCh37", "hg38", "GRCh38"} or None, optional
        Genome-build context for ``chr_pos`` workflows. Default is ``"hg38"``.
    restrict_snps_path : str or os.PathLike[str] or None, optional
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
    genome_build: GenomeBuildInput | None = "hg38"
    restrict_snps_path: str | PathLike[str] | None = None
    log_level: LogLevel = "INFO"
    fail_on_missing_metadata: bool = False

    def __post_init__(self) -> None:
        """Normalize shared path-like fields and validate common enum values."""
        if self.snp_identifier not in {"rsid", "chr_pos"}:
            raise ValueError("snp_identifier must be 'rsid' or 'chr_pos'.")
        object.__setattr__(self, "genome_build", normalize_genome_build(self.genome_build))
        object.__setattr__(self, "log_level", _normalize_log_level(self.log_level))
        object.__setattr__(
            self, "restrict_snps_path", _normalize_optional_path(self.restrict_snps_path)
        )


_GLOBAL_CONFIG: GlobalConfig = GlobalConfig()
_GLOBAL_CONFIG_BANNER_SUPPRESSION: ContextVar[int] = ContextVar(
    "ldsc_global_config_banner_suppression",
    default=0,
)


def get_global_config() -> GlobalConfig:
    """Return the process-wide default configuration used by Python workflows."""
    return _GLOBAL_CONFIG


def set_global_config(config: GlobalConfig) -> None:
    """Replace the process-wide default configuration used by Python workflows."""
    if not isinstance(config, GlobalConfig):
        raise TypeError(f"config must be a GlobalConfig instance, got {type(config).__name__}.")
    global _GLOBAL_CONFIG
    _GLOBAL_CONFIG = config


def reset_global_config() -> GlobalConfig:
    """Restore the process-wide default configuration and return it."""
    global _GLOBAL_CONFIG
    _GLOBAL_CONFIG = GlobalConfig()
    return _GLOBAL_CONFIG


def validate_config_compatibility(
    a: GlobalConfig,
    b: GlobalConfig,
    context: str = "",
) -> None:
    """Raise when ``a`` and ``b`` disagree on critical LDSC workflow settings."""
    prefix = f" when combining {context}" if context else ""
    if a.genome_build != b.genome_build:
        raise ConfigMismatchError(
            f"genome_build mismatch{prefix}: {a.genome_build!r} vs "
            f"{b.genome_build!r}. These objects were computed under different "
            "genome-build assumptions and cannot be safely merged."
        )
    if a.snp_identifier != b.snp_identifier:
        raise ConfigMismatchError(
            f"snp_identifier mismatch{prefix}: {a.snp_identifier!r} vs "
            f"{b.snp_identifier!r}. These objects were computed under different "
            "SNP-identifier modes and cannot be safely merged."
        )
    if a.restrict_snps_path != b.restrict_snps_path:
        warnings.warn(
            f"restrict_snps_path differs{prefix}: {a.restrict_snps_path!r} vs "
            f"{b.restrict_snps_path!r}. Results are combinable but SNP sets may differ.",
            UserWarning,
            stacklevel=2,
        )


@contextmanager
def suppress_global_config_banner():
    """Temporarily suppress plain-text GlobalConfig banners within this context."""
    token = _GLOBAL_CONFIG_BANNER_SUPPRESSION.set(_GLOBAL_CONFIG_BANNER_SUPPRESSION.get() + 1)
    try:
        yield
    finally:
        _GLOBAL_CONFIG_BANNER_SUPPRESSION.reset(token)


def print_global_config_banner(entrypoint: str, global_config: GlobalConfig) -> None:
    """Print the active GlobalConfig for one public workflow entrypoint."""
    if _GLOBAL_CONFIG_BANNER_SUPPRESSION.get() > 0:
        return
    print(f"{entrypoint} using {global_config!r}")


@dataclass(frozen=True)
class AnnotationBuildConfig:
    """Configuration for SNP-level annotation assembly and BED projection.

    Parameters
    ----------
    baseline_annotation_paths : str, os.PathLike[str], or sequence of those, optional
        Baseline annotation files to align and combine. Default is ``()``.
    query_annotation_paths : str, os.PathLike[str], or sequence of those, optional
        Query annotation files to align with the baseline. Default is ``()``.
    query_bed_paths : str, os.PathLike[str], or sequence of those, optional
        BED inputs that should be projected to SNP-level annotations. Default is
        ``()``.
    out_prefix : str or os.PathLike[str] or None, optional
        Output prefix used by file-writing helpers. Default is ``None``.
    batch_mode : bool, optional
        If ``True``, write one output directory per BED input in the batch
        projection workflow. Default is ``True``.
    compression : {"auto", "gzip", "bz2", "none"}, optional
        Output compression preference for generated annotation files. Default is
        ``"gzip"``.
    """
    baseline_annotation_paths: str | PathLike[str] | tuple[str | PathLike[str], ...] | list[str | PathLike[str]] = field(default_factory=tuple)
    query_annotation_paths: str | PathLike[str] | tuple[str | PathLike[str], ...] | list[str | PathLike[str]] = field(default_factory=tuple)
    query_bed_paths: str | PathLike[str] | tuple[str | PathLike[str], ...] | list[str | PathLike[str]] = field(default_factory=tuple)
    out_prefix: str | PathLike[str] | None = None
    batch_mode: bool = True
    compression: CompressionMode = "gzip"

    def __post_init__(self) -> None:
        """Normalize annotation path tokens and validate compression mode."""
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
    plink_prefix : str or os.PathLike[str] or None, optional
        PLINK ``.bed/.bim/.fam`` prefix token. Default is ``None``.
    parquet_r2_paths : str, os.PathLike[str], or sequence of those, optional
        Sorted parquet R2 path tokens. Default is ``()``.
    frequency_paths : str, os.PathLike[str], or sequence of those, optional
        Sidecar frequency or metadata path tokens aligned to the reference
        panel. Default is ``()``.
    r2_bias_mode : {"raw", "unbiased"} or None, optional
        Whether parquet R2 values need sample-size correction. Default is
        ``None``.
    r2_sample_size : float or None, optional
        Sample size used when correcting raw parquet R2 values. Default is
        ``None``.
    """
    backend: RefPanelBackend = "auto"
    plink_prefix: str | PathLike[str] | None = None
    parquet_r2_paths: str | PathLike[str] | tuple[str | PathLike[str], ...] | list[str | PathLike[str]] = field(default_factory=tuple)
    frequency_paths: str | PathLike[str] | tuple[str | PathLike[str], ...] | list[str | PathLike[str]] = field(default_factory=tuple)
    r2_bias_mode: R2BiasMode | None = None
    r2_sample_size: float | None = None

    def __post_init__(self) -> None:
        """Normalize backend path tokens and validate parquet-R2 settings."""
        if self.backend not in {"auto", "plink", "parquet_r2"}:
            raise ValueError("backend must be 'auto', 'plink', or 'parquet_r2'.")
        if self.r2_bias_mode not in {None, "raw", "unbiased"}:
            raise ValueError("r2_bias_mode must be None, 'raw', or 'unbiased'.")
        if self.r2_sample_size is not None and self.r2_sample_size <= 0:
            raise ValueError("r2_sample_size must be positive when provided.")
        object.__setattr__(self, "plink_prefix", _normalize_optional_path(self.plink_prefix))
        object.__setattr__(self, "parquet_r2_paths", _normalize_path_tuple(self.parquet_r2_paths))
        object.__setattr__(self, "frequency_paths", _normalize_path_tuple(self.frequency_paths))


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
    keep_individuals_path : str or os.PathLike[str] or None, optional
        Optional path to a one-column IID keep file applied in PLINK mode
        before SNP and MAF filtering. Default is ``None``.
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
    keep_individuals_path: str | PathLike[str] | None = None
    chunk_size: int = 50
    compute_m5_50: bool = True
    whole_chromosome_ok: bool = False

    def __post_init__(self) -> None:
        """Validate LD-window settings and normalize optional keep-file paths."""
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
        object.__setattr__(self, "keep_individuals_path", _normalize_optional_path(self.keep_individuals_path))


@dataclass(frozen=True)
class ReferencePanelBuildConfig:
    """Configuration for building standard parquet reference panels from PLINK."""

    panel_label: str
    plink_prefix: str | PathLike[str]
    source_genome_build: GenomeBuildInput
    genetic_map_hg19_path: str | PathLike[str]
    genetic_map_hg38_path: str | PathLike[str]
    output_dir: str | PathLike[str]
    liftover_chain_hg19_to_hg38_path: str | PathLike[str] | None = None
    liftover_chain_hg38_to_hg19_path: str | PathLike[str] | None = None
    ld_wind_snps: int | None = None
    ld_wind_kb: float | None = None
    ld_wind_cm: float | None = None
    maf_min: float | None = None
    restrict_snps_path: str | PathLike[str] | None = None
    keep_indivs_path: str | PathLike[str] | None = None
    chunk_size: int = 50

    def __post_init__(self) -> None:
        """Normalize build paths and validate liftover and LD-window settings."""
        if not str(self.panel_label).strip():
            raise ValueError("panel_label must be non-empty.")
        object.__setattr__(self, "panel_label", str(self.panel_label).strip())
        object.__setattr__(self, "plink_prefix", _normalize_required_path(self.plink_prefix))
        object.__setattr__(self, "source_genome_build", normalize_genome_build(self.source_genome_build))
        object.__setattr__(self, "genetic_map_hg19_path", _normalize_required_path(self.genetic_map_hg19_path))
        object.__setattr__(self, "genetic_map_hg38_path", _normalize_required_path(self.genetic_map_hg38_path))
        object.__setattr__(
            self,
            "liftover_chain_hg19_to_hg38_path",
            _normalize_optional_path(self.liftover_chain_hg19_to_hg38_path),
        )
        object.__setattr__(
            self,
            "liftover_chain_hg38_to_hg19_path",
            _normalize_optional_path(self.liftover_chain_hg38_to_hg19_path),
        )
        object.__setattr__(self, "output_dir", _normalize_required_path(self.output_dir))
        object.__setattr__(self, "restrict_snps_path", _normalize_optional_path(self.restrict_snps_path))
        object.__setattr__(self, "keep_indivs_path", _normalize_optional_path(self.keep_indivs_path))
        required_chain = (
            self.liftover_chain_hg19_to_hg38_path
            if self.source_genome_build == "hg19"
            else self.liftover_chain_hg38_to_hg19_path
        )
        if required_chain is None:
            if self.source_genome_build == "hg19":
                raise ValueError(
                    "liftover_chain_hg19_to_hg38_path is required when source_genome_build is hg19."
                )
            raise ValueError(
                "liftover_chain_hg38_to_hg19_path is required when source_genome_build is hg38."
            )
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
    out_prefix: str | PathLike[str]
    N: float | None = None
    N_cas: float | None = None
    N_con: float | None = None
    info_min: float = 0.9
    maf_min: float = 0.01
    n_min: float | None = None
    nstudy_min: float | None = None
    chunk_size: int = int(5e6)
    merge_alleles_path: str | PathLike[str] | None = None
    signed_sumstats_spec: str | None = None
    ignore_columns: tuple[str, ...] = field(default_factory=tuple)
    no_alleles: bool = False
    a1_inc: bool = False
    keep_maf: bool = False
    daner: bool = False
    daner_n: bool = False

    def __post_init__(self) -> None:
        """Validate munging thresholds and normalize optional file paths."""
        if self.info_min < 0:
            raise ValueError("info_min must be non-negative.")
        if self.maf_min < 0 or self.maf_min > 0.5:
            raise ValueError("maf_min must lie in [0, 0.5].")
        if self.chunk_size <= 0:
            raise ValueError("chunk_size must be positive.")
        object.__setattr__(self, "out_prefix", _normalize_required_path(self.out_prefix))
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
        """Validate regression hyperparameters after dataclass construction."""
        if self.n_blocks <= 1:
            raise ValueError("n_blocks must be greater than 1.")
        if self.two_step_cutoff is not None and self.two_step_cutoff <= 0:
            raise ValueError("two_step_cutoff must be positive when provided.")
        if self.chisq_max is not None and self.chisq_max <= 0:
            raise ValueError("chisq_max must be positive when provided.")
