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
import logging
from os import PathLike
from typing import Literal
import warnings

from .column_inference import normalize_genome_build
from .errors import LDSCConfigError
from .path_resolution import normalize_optional_path_token, normalize_path_token, normalize_path_tokens

LOGGER = logging.getLogger("LDSC.config")
_GENOME_BUILD_UNSET = object()


SNPIdentifierMode = Literal["rsid", "chr_pos"]
GenomeBuild = Literal["hg19", "hg38", "auto"]
GenomeBuildInput = Literal["auto", "hg19", "hg37", "hg38", "GRCh37", "GRCh38"]
LogLevel = Literal["DEBUG", "INFO", "WARNING", "ERROR"]
RefPanelBackend = Literal["auto", "plink", "parquet_r2"]
CompressionMode = Literal["auto", "gzip", "bz2", "none"]
R2BiasMode = Literal["raw", "unbiased"]


class ConfigMismatchError(LDSCConfigError, ValueError):
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


@dataclass(frozen=True, init=False)
class GlobalConfig:
    """Shared configuration used across the refactored workflows.

    Parameters
    ----------
    snp_identifier : {"rsid", "chr_pos"}, optional
        Global SNP identifier mode. Default is ``"chr_pos"``. ``"rsid"``
        expects an explicit SNP column, while ``"chr_pos"`` builds identifiers
        from chromosome and base-pair position.
    genome_build : {"auto", "hg19", "hg37", "GRCh37", "hg38", "GRCh38"} or None, optional
        Genome-build context for ``chr_pos`` workflows that require
        coordinate-build interpretation. Default is ``"auto"``. Ignored for
        ``"rsid"``.
    log_level : {"DEBUG", "INFO", "WARNING", "ERROR"}, optional
        Requested logging verbosity for workflow modules. Default is ``"INFO"``.
    fail_on_missing_metadata : bool, optional
        If ``True``, treat missing optional metadata as a hard error rather than
        tolerating partial metadata tables. Default is ``False``.

    Notes
    -----
    ``GlobalConfig`` intentionally carries only shared runtime assumptions.
    Per-run SNP-universe controls such as ``ref_panel_snps_file`` and
    ``regression_snps_file`` now live on workflow-specific configs.
    """
    snp_identifier: SNPIdentifierMode = "chr_pos"
    genome_build: GenomeBuildInput | None = "auto"
    log_level: LogLevel = "INFO"
    fail_on_missing_metadata: bool = False

    def __init__(
        self,
        snp_identifier: SNPIdentifierMode = "chr_pos",
        genome_build: GenomeBuildInput | None | object = _GENOME_BUILD_UNSET,
        log_level: LogLevel = "INFO",
        fail_on_missing_metadata: bool = False,
    ) -> None:
        """Initialize global workflow assumptions with mode-aware defaults."""
        if genome_build is _GENOME_BUILD_UNSET:
            genome_build = None if snp_identifier == "rsid" else "auto"
        object.__setattr__(self, "snp_identifier", snp_identifier)
        object.__setattr__(self, "genome_build", genome_build)
        object.__setattr__(self, "log_level", log_level)
        object.__setattr__(self, "fail_on_missing_metadata", fail_on_missing_metadata)
        self.__post_init__()

    def __post_init__(self) -> None:
        """Normalize shared path-like fields and validate common enum values."""
        if self.snp_identifier not in {"rsid", "chr_pos"}:
            raise ValueError("snp_identifier must be 'rsid' or 'chr_pos'.")
        object.__setattr__(self, "genome_build", normalize_genome_build(self.genome_build))
        object.__setattr__(self, "log_level", _normalize_log_level(self.log_level))
        if self.snp_identifier == "chr_pos" and self.genome_build is None:
            raise ValueError(
                "genome_build is required when snp_identifier='chr_pos'. "
                "Pass genome_build='auto' to infer from data, or 'hg19'/'hg38' explicitly."
            )
        if self.snp_identifier == "rsid" and self.genome_build == "auto":
            raise ValueError("genome_build='auto' is not valid for snp_identifier='rsid'.")
        if self.snp_identifier == "rsid" and self.genome_build is not None:
            warnings.warn(
                "genome_build is set but will be ignored in rsid mode.",
                UserWarning,
                stacklevel=3,
            )
            object.__setattr__(self, "genome_build", None)


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
    if a.snp_identifier != b.snp_identifier:
        raise ConfigMismatchError(
            f"snp_identifier mismatch{prefix}: {a.snp_identifier!r} vs "
            f"{b.snp_identifier!r}. These objects were computed under different "
            "SNP-identifier modes and cannot be safely merged."
        )
    if a.genome_build != b.genome_build:
        raise ConfigMismatchError(
            f"genome_build mismatch{prefix}: {a.genome_build!r} vs "
            f"{b.genome_build!r}. These objects were computed under different "
            "genome-build assumptions and cannot be safely merged."
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
    """Log the active GlobalConfig for one public workflow entrypoint."""
    if _GLOBAL_CONFIG_BANNER_SUPPRESSION.get() > 0:
        return
    LOGGER.info(f"{entrypoint} using {global_config!r}")


@dataclass(frozen=True)
class AnnotationBuildConfig:
    """Configuration for SNP-level annotation assembly and BED projection.

    Parameters
    ----------
    baseline_annot_sources : str, os.PathLike[str], or sequence of those, optional
        Baseline annotation files to align and combine. Default is ``()``.
    query_annot_sources : str, os.PathLike[str], or sequence of those, optional
        Query annotation files to align with the baseline. Default is ``()``.
    query_annot_bed_sources : str, os.PathLike[str], or sequence of those, optional
        BED inputs that should be projected to SNP-level annotations. Default is
        ``()``.
    output_dir : str or os.PathLike[str] or None, optional
        Output directory used by file-writing helpers. Default is ``None``.
    compression : {"auto", "gzip", "bz2", "none"}, optional
        Output compression preference for generated annotation files. Default is
        ``"gzip"``.
    overwrite : bool, optional
        If ``True``, replace generated annotation outputs and remove stale
        owned root-level ``query.*.annot.gz`` siblings after a successful
        projection. If ``False``, output collisions raise before writing
        starts. Default is ``False``.
    """
    baseline_annot_sources: str | PathLike[str] | tuple[str | PathLike[str], ...] | list[str | PathLike[str]] = field(default_factory=tuple)
    query_annot_sources: str | PathLike[str] | tuple[str | PathLike[str], ...] | list[str | PathLike[str]] = field(default_factory=tuple)
    query_annot_bed_sources: str | PathLike[str] | tuple[str | PathLike[str], ...] | list[str | PathLike[str]] = field(default_factory=tuple)
    output_dir: str | PathLike[str] | None = None
    compression: CompressionMode = "gzip"
    allow_missing_query: bool = True
    overwrite: bool = False

    def __post_init__(self) -> None:
        """Normalize annotation path tokens and validate compression mode."""
        object.__setattr__(self, "baseline_annot_sources", _normalize_path_tuple(self.baseline_annot_sources))
        object.__setattr__(self, "query_annot_sources", _normalize_path_tuple(self.query_annot_sources))
        object.__setattr__(self, "query_annot_bed_sources", _normalize_path_tuple(self.query_annot_bed_sources))
        object.__setattr__(self, "output_dir", _normalize_optional_path(self.output_dir))
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
    r2_dir : str or os.PathLike[str] or None, optional
        Directory containing package-built parquet R2 files named
        ``chr{chrom}_r2.parquet`` and optional metadata sidecars named
        ``chr{chrom}_meta.tsv.gz``. Default is ``None``.
    maf_min : float or None, optional
        Optional retained-reference-panel MAF threshold applied before LD-score
        computation. This is distinct from the common-SNP count threshold on
        ``LDScoreConfig``. Default is ``None``.
    keep_indivs_file : str or os.PathLike[str] or None, optional
        Optional path to a one-column IID keep file applied in PLINK mode before
        genotype-derived MAF is computed. Default is ``None``.
    r2_bias_mode : {"raw", "unbiased"} or None, optional
        Whether parquet R2 values need sample-size correction. Default is
        ``None``.
    r2_sample_size : float or None, optional
        Sample size used when correcting raw parquet R2 values. Default is
        ``None``.
    """
    backend: RefPanelBackend = "auto"
    plink_prefix: str | PathLike[str] | None = None
    r2_dir: str | PathLike[str] | None = None
    chromosomes: tuple[str, ...] | list[str] | None = None
    maf_min: float | None = None
    keep_indivs_file: str | PathLike[str] | None = None
    r2_bias_mode: R2BiasMode | None = None
    r2_sample_size: float | None = None
    sample_size: int | None = None
    ref_panel_snps_file: str | PathLike[str] | None = None

    def __post_init__(self) -> None:
        """Normalize backend path tokens and validate parquet-R2 settings."""
        if self.backend not in {"auto", "plink", "parquet_r2"}:
            raise ValueError("backend must be 'auto', 'plink', or 'parquet_r2'.")
        if self.r2_bias_mode not in {None, "raw", "unbiased"}:
            raise ValueError("r2_bias_mode must be None, 'raw', or 'unbiased'.")
        if self.r2_sample_size is not None and self.r2_sample_size <= 0:
            raise ValueError("r2_sample_size must be positive when provided.")
        if self.sample_size is not None and self.sample_size <= 0:
            raise ValueError("sample_size must be positive when provided.")
        if self.maf_min is not None and not 0 <= self.maf_min <= 0.5:
            raise ValueError("maf_min must lie in [0, 0.5].")
        if self.sample_size is None and self.r2_sample_size is not None:
            object.__setattr__(self, "sample_size", int(self.r2_sample_size))
        object.__setattr__(self, "plink_prefix", _normalize_optional_path(self.plink_prefix))
        object.__setattr__(self, "r2_dir", _normalize_optional_path(self.r2_dir))
        object.__setattr__(self, "keep_indivs_file", _normalize_optional_path(self.keep_indivs_file))
        object.__setattr__(self, "ref_panel_snps_file", _normalize_optional_path(self.ref_panel_snps_file))
        if self.chromosomes is not None:
            from .chromosome_inference import normalize_chromosome

            object.__setattr__(self, "chromosomes", tuple(normalize_chromosome(chrom) for chrom in self.chromosomes))


@dataclass(frozen=True)
class LDScoreConfig:
    """Configuration for chromosome-wise LD-score calculation.

    Exactly one LD-window field must be provided. The remaining fields control
    regression row selection and count-artifact emission.

    Parameters
    ----------
    ld_wind_snps : int or None, optional
        Window size measured in SNP count. Default is ``None``.
    ld_wind_kb : float or None, optional
        Window size measured in kilobases. Default is ``None``.
    ld_wind_cm : float or None, optional
        Window size measured in centiMorgans. Default is ``None``.
    regression_snps_file : str or os.PathLike[str] or None, optional
        Optional path to the SNP list defining the regression SNP set used for
        the persisted ``ldscore.baseline.parquet`` row set and, when query
        annotations are present, the aligned ``ldscore.query.parquet`` row set.
        Default is ``None``.
    snp_batch_size : int, optional
        Number of SNPs processed per LD-score sliding batch. Default is
        ``128``.
    common_maf_min : float, optional
        Inclusive MAF threshold used only for common-SNP count vectors
        (``MAF >= common_maf_min``). It does not change retained reference
        SNPs, LD-score rows, LD scores, or regression weights. Default is
        ``0.05``.
    whole_chromosome_ok : bool, optional
        Override the guard that rejects windows effectively spanning an entire
        chromosome. Default is ``False``.
    """
    ld_wind_snps: int | None = None
    ld_wind_kb: float | None = None
    ld_wind_cm: float | None = None
    regression_snps_file: str | PathLike[str] | None = None
    snp_batch_size: int = 128
    common_maf_min: float = 0.05
    whole_chromosome_ok: bool = False

    def __post_init__(self) -> None:
        """Validate LD-window settings and normalize optional SNP-list paths."""
        windows = [self.ld_wind_snps, self.ld_wind_kb, self.ld_wind_cm]
        if sum(value is not None for value in windows) != 1:
            raise ValueError("Exactly one LD-window option must be set.")
        if self.ld_wind_snps is not None and self.ld_wind_snps <= 0:
            raise ValueError("ld_wind_snps must be positive.")
        if self.ld_wind_kb is not None and self.ld_wind_kb <= 0:
            raise ValueError("ld_wind_kb must be positive.")
        if self.ld_wind_cm is not None and self.ld_wind_cm <= 0:
            raise ValueError("ld_wind_cm must be positive.")
        if not 0 <= self.common_maf_min <= 0.5:
            raise ValueError("common_maf_min must lie in [0, 0.5].")
        if self.snp_batch_size <= 0:
            raise ValueError("snp_batch_size must be positive.")
        object.__setattr__(self, "regression_snps_file", _normalize_optional_path(self.regression_snps_file))


@dataclass(frozen=True)
class ReferencePanelBuildConfig:
    """Configuration for building standard parquet reference panels from PLINK.

    This config owns the builder-only SNP restriction path
    ``ref_panel_snps_file``. That setting is not shared through
    ``GlobalConfig`` because it changes the retained panel rows for one build
    run rather than the global identifier or genome-build assumptions.

    Parameters
    ----------
    plink_prefix : str or os.PathLike[str]
        PLINK ``.bed/.bim/.fam`` prefix token. This may be a single prefix or an
        explicit ``@`` chromosome-suite token.
    source_genome_build : {"hg19", "hg37", "GRCh37", "hg38", "GRCh38"} or None, optional
        Genome build of the input PLINK coordinates. If ``None``, the
        build-ref-panel workflow infers the build from PLINK ``.bim`` rows
        before applying SNP restrictions.
    genetic_map_hg19_sources, genetic_map_hg38_sources : str or os.PathLike[str] or None, optional
        Genetic-map paths used to populate cM values for each emitted sidecar.
        Maps are required for every emitted build when ``ld_wind_cm`` is set.
        When maps are omitted for SNP- or kb-window builds, emitted sidecars
        store missing CM values.
    output_dir : str or os.PathLike[str]
        Output directory for the reference panel. The panel identity is
        ``Path(output_dir).name``; emitted artifacts are fixed as
        ``{build}/chr{chrom}_r2.parquet`` and
        ``{build}/chr{chrom}_meta.tsv.gz``.
    liftover_chain_hg19_to_hg38_file, liftover_chain_hg38_to_hg19_file : str or os.PathLike[str] or None, optional
        Chain files used to populate the opposite-build coordinates. If the
        chain matching the resolved ``source_genome_build`` is omitted, the
        builder emits a source-build-only panel.
    ld_wind_snps, ld_wind_kb, ld_wind_cm : int, float, float, or None, optional
        LD window specification. Exactly one must be supplied.
    maf_min : float or None, optional
        Optional retained-SNP MAF threshold. Default is ``None``.
    ref_panel_snps_file : str or os.PathLike[str] or None, optional
        Optional SNP list restricting the emitted reference-panel universe.
        Default is ``None``. The identifier mode comes from
        ``GlobalConfig.snp_identifier``. In ``chr_pos`` mode, the restriction
        file must be aligned to the resolved source reference-panel build; the
        builder never uses ``GlobalConfig.genome_build``.
    keep_indivs_file : str or os.PathLike[str] or None, optional
        Optional individual keep-file applied before R2 calculation. Default is
        ``None``.
    snp_batch_size : int, optional
        Number of SNPs loaded per pairwise-R2 computation batch. Larger values
        may improve throughput but use more memory. Default is ``128``.
    overwrite : bool, optional
        If ``True``, replace current candidate panel artifacts. This expert
        workflow does not clean stale optional target-build or ``dropped_snps``
        siblings from earlier configurations. If ``False``, output collisions
        raise before chromosome processing starts. Default is ``False``.
    """

    plink_prefix: str | PathLike[str]
    source_genome_build: GenomeBuildInput | None = None
    genetic_map_hg19_sources: str | PathLike[str] | None = None
    genetic_map_hg38_sources: str | PathLike[str] | None = None
    output_dir: str | PathLike[str] | None = None
    liftover_chain_hg19_to_hg38_file: str | PathLike[str] | None = None
    liftover_chain_hg38_to_hg19_file: str | PathLike[str] | None = None
    ld_wind_snps: int | None = None
    ld_wind_kb: float | None = None
    ld_wind_cm: float | None = None
    maf_min: float | None = None
    ref_panel_snps_file: str | PathLike[str] | None = None
    keep_indivs_file: str | PathLike[str] | None = None
    snp_batch_size: int = 128
    overwrite: bool = False
    duplicate_position_policy: str = "error"

    def __post_init__(self) -> None:
        """Normalize build paths and validate liftover and LD-window settings."""
        object.__setattr__(self, "plink_prefix", _normalize_required_path(self.plink_prefix))
        object.__setattr__(self, "source_genome_build", normalize_genome_build(self.source_genome_build))
        if self.source_genome_build == "auto":
            raise ValueError("source_genome_build must be hg19/hg38 or omitted for inference.")
        object.__setattr__(self, "genetic_map_hg19_sources", _normalize_optional_path(self.genetic_map_hg19_sources))
        object.__setattr__(self, "genetic_map_hg38_sources", _normalize_optional_path(self.genetic_map_hg38_sources))
        object.__setattr__(
            self,
            "liftover_chain_hg19_to_hg38_file",
            _normalize_optional_path(self.liftover_chain_hg19_to_hg38_file),
        )
        object.__setattr__(
            self,
            "liftover_chain_hg38_to_hg19_file",
            _normalize_optional_path(self.liftover_chain_hg38_to_hg19_file),
        )
        object.__setattr__(self, "output_dir", _normalize_required_path(self.output_dir))
        object.__setattr__(self, "ref_panel_snps_file", _normalize_optional_path(self.ref_panel_snps_file))
        object.__setattr__(self, "keep_indivs_file", _normalize_optional_path(self.keep_indivs_file))
        windows = [self.ld_wind_snps, self.ld_wind_kb, self.ld_wind_cm]
        if sum(value is not None for value in windows) != 1:
            raise ValueError("Exactly one LD-window option must be set.")
        source_map = None
        if self.source_genome_build == "hg19":
            source_map = self.genetic_map_hg19_sources
        elif self.source_genome_build == "hg38":
            source_map = self.genetic_map_hg38_sources
        if self.ld_wind_cm is not None and self.source_genome_build is not None and source_map is None:
            raise ValueError(f"{self.source_genome_build} genetic map path is required when ld_wind_cm is set.")
        if self.ld_wind_snps is not None and self.ld_wind_snps <= 0:
            raise ValueError("ld_wind_snps must be positive.")
        if self.ld_wind_kb is not None and self.ld_wind_kb <= 0:
            raise ValueError("ld_wind_kb must be positive.")
        if self.ld_wind_cm is not None and self.ld_wind_cm <= 0:
            raise ValueError("ld_wind_cm must be positive.")
        if self.maf_min is not None and not 0 <= self.maf_min <= 0.5:
            raise ValueError("maf_min must lie in [0, 0.5].")
        if self.snp_batch_size <= 0:
            raise ValueError("snp_batch_size must be positive.")
        if self.duplicate_position_policy not in {"error", "drop-all"}:
            raise ValueError(
                "duplicate_position_policy must be 'error' or 'drop-all', "
                f"got {self.duplicate_position_policy!r}."
            )


@dataclass(frozen=True)
class MungeConfig:
    """Configuration for legacy-compatible summary-statistics munging.

    This dataclass preserves the established LDSC munging behavior while
    exposing the options through an explicit Python object rather than a script
    namespace. Coordinate hints for ``CHR`` and ``POS`` may be supplied through
    ``column_hints``; when omitted, the shared column-inference registry accepts
    common aliases such as ``#CHROM``, ``CHROM``, ``CHR``, ``POS``, and ``BP``.

    Parameters
    ----------
    output_dir : str or os.PathLike[str]
        Directory that receives workflow-owned ``sumstats.parquet`` and/or
        ``sumstats.sumstats.gz``, ``sumstats.log``, and
        ``sumstats.metadata.json`` artifacts.
    raw_sumstats_file : str or os.PathLike[str] or None, optional
        Raw summary-statistics file to munge. Exact-one glob patterns are
        resolved by the workflow before entering the legacy kernel. Default is
        ``None``.
    N, N_cas, N_con : float or None, optional
        Sample-size overrides forwarded to the munging kernel. Defaults are
        ``None``.
    info_min : float, optional
        Minimum INFO score. Default is ``0.9``.
    maf_min : float, optional
        Minimum allele frequency. Default is ``0.01``.
    n_min, nstudy_min : float or None, optional
        Optional row filters for sample size and study count. Defaults are
        ``None``.
    chunk_size : int, optional
        Number of input rows processed per chunk. Default is ``1_000_000``.
    output_format : {"parquet", "tsv.gz", "both"}, optional
        Curated sumstats disk format written by the public workflow. Default is
        ``"parquet"``.
    sumstats_snps_file : str or os.PathLike[str] or None, optional
        Optional headered summary-statistics SNP keep-list path. In ``rsid``
        mode, central ``SNP`` aliases identify the keep-list column. In
        ``chr_pos`` mode, central ``CHR``/``POS`` aliases, including
        build-specific position aliases such as ``hg19_POS`` and ``hg38_POS``,
        define retained coordinates. This option restricts rows only; it does
        not allele-match, rewrite alleles, or reorder output. Default is
        ``None``.
    signed_sumstats_spec : str or None, optional
        Signed statistic specification passed through to the legacy kernel.
        Default is ``None``.
    ignore_columns : tuple of str, optional
        Source columns ignored during auto-detection. Default is ``()``.
    no_alleles, a1_inc, keep_maf, daner_old, daner_new : bool, optional
        Legacy munging switches preserved for behavior compatibility. Defaults
        are ``False``.
    overwrite : bool, optional
        If ``True``, replace current fixed sumstats outputs and remove stale
        owned ``sumstats.*`` siblings after a successful run. If ``False``, any
        owned sumstats artifact collision raises before the munging kernel
        runs. Default is ``False``.
    """
    output_dir: str | PathLike[str] | None = None
    raw_sumstats_file: str | PathLike[str] | None = None
    compression: str = "auto"
    output_format: str = "parquet"
    trait_name: str | None = None
    column_hints: dict[str, str] = field(default_factory=dict)
    N: float | None = None
    N_cas: float | None = None
    N_con: float | None = None
    info_min: float = 0.9
    maf_min: float = 0.01
    n_min: float | None = None
    nstudy_min: float | None = None
    chunk_size: int = 1_000_000
    sumstats_snps_file: str | PathLike[str] | None = None
    signed_sumstats_spec: str | None = None
    ignore_columns: tuple[str, ...] = field(default_factory=tuple)
    no_alleles: bool = False
    a1_inc: bool = False
    keep_maf: bool = False
    daner_old: bool = False
    daner_new: bool = False
    overwrite: bool = False

    def __post_init__(self) -> None:
        """Validate munging thresholds and normalize optional file paths."""
        if self.info_min < 0:
            raise ValueError("info_min must be non-negative.")
        if self.maf_min < 0 or self.maf_min > 0.5:
            raise ValueError("maf_min must lie in [0, 0.5].")
        if self.chunk_size <= 0:
            raise ValueError("chunk_size must be positive.")
        if self.output_format not in {"parquet", "tsv.gz", "both"}:
            raise ValueError("output_format must be one of 'parquet', 'tsv.gz', or 'both'.")
        object.__setattr__(self, "output_dir", _normalize_optional_path(self.output_dir))
        object.__setattr__(self, "raw_sumstats_file", _normalize_optional_path(self.raw_sumstats_file))
        object.__setattr__(self, "sumstats_snps_file", _normalize_optional_path(self.sumstats_snps_file))
        object.__setattr__(self, "ignore_columns", tuple(self.ignore_columns))
        object.__setattr__(self, "column_hints", dict(self.column_hints))


@dataclass(frozen=True)
class RegressionConfig:
    """Configuration for heritability and genetic-correlation regressions.

    Parameters
    ----------
    n_blocks : int, optional
        Requested number of block-jackknife partitions. Default is ``200``.
    use_common_counts : bool, optional
        If ``True``, prefer the manifest ``common_reference_snp_counts`` vector
        when it is available; otherwise use ``all_reference_snp_counts``.
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
    use_common_counts: bool = True
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
