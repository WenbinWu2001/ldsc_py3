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

from ._kernel.snp_identity import identity_mode_family, normalize_snp_identifier_mode
from .column_inference import normalize_genome_build
from .errors import LDSCConfigError
from .path_resolution import normalize_optional_path_token, normalize_path_token, normalize_path_tokens

LOGGER = logging.getLogger("LDSC.config")
_GENOME_BUILD_UNSET = object()


SNPIdentifierMode = Literal["rsid", "rsid_allele_aware", "chr_pos", "chr_pos_allele_aware"]
GenomeBuild = Literal["hg19", "hg38", "auto"]
GenomeBuildInput = Literal["auto", "hg19", "hg37", "hg38", "GRCh37", "GRCh38"]
LogLevel = Literal["DEBUG", "INFO", "WARNING", "ERROR"]
LOG_LEVEL_CHOICES: tuple[LogLevel, ...] = ("DEBUG", "INFO", "WARNING", "ERROR")
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
        raise LDSCConfigError(
            "Could not construct workflow config: a required path value is missing. "
            "Most likely a required input or output path argument was omitted. "
            "Provide the path before starting the workflow."
        )
    return normalize_path_token(path)


def _normalize_path_tuple(values) -> tuple[str, ...]:
    """Convert one or many path-like objects to a string tuple of tokens."""
    return normalize_path_tokens(values)


def _normalize_trait_name(value: str | None) -> str | None:
    """Return a stripped trait label, rejecting blank user-provided labels."""
    if value is None:
        return None
    normalized = str(value).strip()
    if not normalized:
        raise LDSCConfigError(
            "Could not construct MungeConfig: trait_name is blank. "
            "Most likely whitespace was passed as the trait label. "
            "Use a non-empty trait name or omit trait_name to derive it from the input."
        )
    return normalized


def _normalize_log_level(level: str) -> LogLevel:
    """Normalize a logging level string to the supported uppercase literal."""
    normalized = level.upper()
    allowed = set(LOG_LEVEL_CHOICES)
    if normalized not in allowed:
        raise LDSCConfigError(
            f"Could not construct GlobalConfig: log_level={level!r} is invalid. "
            "Most likely the log level name is misspelled. "
            f"Use one of {', '.join(LOG_LEVEL_CHOICES)}."
        )
    return normalized  # type: ignore[return-value]


def _invalid_choice_message(config_name: str, field_name: str, value, choices: str) -> str:
    """Return a self-contained message for a config enum validation failure."""
    return (
        f"Could not construct {config_name}: {field_name}={value!r} is invalid. "
        f"Most likely the option value is misspelled or from an older API. "
        f"Use one of {choices}."
    )


def _positive_number_message(config_name: str, field_name: str, value) -> str:
    """Return a self-contained message for positive numeric config fields."""
    return (
        f"Could not construct {config_name}: {field_name}={value!r} must be positive. "
        "Most likely a zero or negative value was supplied. "
        f"Set {field_name} to a value greater than 0."
    )


def _range_message(config_name: str, field_name: str, value, allowed_range: str) -> str:
    """Return a self-contained message for bounded numeric config fields."""
    return (
        f"Could not construct {config_name}: {field_name}={value!r} is outside {allowed_range}. "
        "Most likely the threshold was entered on the wrong scale. "
        f"Set {field_name} within {allowed_range}."
    )


def _mutually_exclusive_message(config_name: str, left: str, right: str) -> str:
    """Return a self-contained message for mutually exclusive config options."""
    return (
        f"Could not construct {config_name}: {left} and {right} are mutually exclusive. "
        "Most likely both a custom file and packaged shortcut were selected for the same input. "
        f"Set only one of {left} or {right}."
    )


def _one_ld_window_message(config_name: str) -> str:
    """Return a self-contained message for LD-window selector validation."""
    return (
        f"Could not construct {config_name}: exactly one LD-window option must be set. "
        "Most likely no window was supplied or multiple window units were combined. "
        "Set exactly one of ld_wind_snps, ld_wind_kb, or ld_wind_cm."
    )


def _validate_region_presets(class_name: str, names: tuple[str, ...]) -> None:
    """Raise if any preset name is outside the supported region menu."""
    from ._kernel.regions import REGION_PRESETS

    unknown = sorted(set(names) - REGION_PRESETS)
    if unknown:
        raise LDSCConfigError(
            f"Could not construct {class_name}: unknown region preset(s) {unknown}. "
            f"Most likely a name was misspelled. Valid presets are {sorted(REGION_PRESETS)}."
        )


@dataclass(frozen=True, init=False)
class GlobalConfig:
    """Shared configuration used across the refactored workflows.

    Parameters
    ----------
    snp_identifier : {"rsid", "rsid_allele_aware", "chr_pos", "chr_pos_allele_aware"}, optional
        Global SNP identifier mode. Default is ``"chr_pos_allele_aware"``.
        Mode names are exact. ``"rsid"`` and ``"chr_pos"`` are allele-blind
        base modes; allele columns are passive for identity. The allele-aware
        variants require usable ``A1/A2`` on package-written identity artifacts
        and add a normalized allele set to the merge key.
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
    snp_identifier: SNPIdentifierMode = "chr_pos_allele_aware"
    genome_build: GenomeBuildInput | None = "auto"
    log_level: LogLevel = "INFO"
    fail_on_missing_metadata: bool = False

    def __init__(
        self,
        snp_identifier: SNPIdentifierMode = "chr_pos_allele_aware",
        genome_build: GenomeBuildInput | None | object = _GENOME_BUILD_UNSET,
        log_level: LogLevel = "INFO",
        fail_on_missing_metadata: bool = False,
    ) -> None:
        """Initialize global workflow assumptions with mode-aware defaults."""
        try:
            snp_identifier = normalize_snp_identifier_mode(snp_identifier)  # type: ignore[assignment]
        except ValueError as exc:
            raise LDSCConfigError(
                f"Could not construct GlobalConfig: snp_identifier={snp_identifier!r} is invalid. "
                "Most likely the identifier mode is misspelled. "
                "Use one of rsid, rsid_allele_aware, chr_pos, or chr_pos_allele_aware."
            ) from exc
        if genome_build is _GENOME_BUILD_UNSET:
            genome_build = None if identity_mode_family(snp_identifier) == "rsid" else "auto"
        object.__setattr__(self, "snp_identifier", snp_identifier)
        object.__setattr__(self, "genome_build", genome_build)
        object.__setattr__(self, "log_level", log_level)
        object.__setattr__(self, "fail_on_missing_metadata", fail_on_missing_metadata)
        self.__post_init__()

    def __post_init__(self) -> None:
        """Normalize shared path-like fields and validate common enum values."""
        mode = normalize_snp_identifier_mode(self.snp_identifier)
        family = identity_mode_family(mode)
        object.__setattr__(self, "snp_identifier", mode)
        try:
            genome_build = normalize_genome_build(self.genome_build)
        except ValueError as exc:
            raise LDSCConfigError(
                f"Could not construct GlobalConfig: genome_build={self.genome_build!r} is invalid. "
                "Most likely the genome build name is misspelled or unsupported. "
                "Use 'auto', 'hg19', 'hg38', 'GRCh37', or 'GRCh38'."
            ) from exc
        object.__setattr__(self, "genome_build", genome_build)
        object.__setattr__(self, "log_level", _normalize_log_level(self.log_level))
        if family == "chr_pos" and self.genome_build is None:
            raise LDSCConfigError(
                "Could not construct GlobalConfig: genome_build is required for chr_pos-family "
                f"snp_identifier={mode!r}. Most likely genome_build=None was passed explicitly. "
                "Pass genome_build='auto' to infer from data, or 'hg19'/'hg38' explicitly."
            )
        if family == "rsid" and self.genome_build == "auto":
            raise LDSCConfigError(
                f"Could not construct GlobalConfig: genome_build='auto' is not valid for rsid-family "
                f"snp_identifier={mode!r}. Most likely an auto-build coordinate setting was reused with rsID mode. "
                "Omit genome_build or pass genome_build=None for rsID-family workflows."
            )
        if family == "rsid" and self.genome_build is not None:
            warnings.warn(
                "genome_build is set but will be ignored in rsid-family mode.",
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
        raise TypeError(
            f"set_global_config expected a GlobalConfig instance, got {type(config).__name__}. "
            "Most likely an unrelated config object was passed. Construct GlobalConfig(...) first."
        )
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
            f"Cannot combine LDSC results{prefix}: snp_identifier mismatch "
            f"{a.snp_identifier!r} vs {b.snp_identifier!r}. Most likely the inputs were produced "
            "with different --snp-identifier modes. Re-run the upstream workflows with the same "
            "snp_identifier before combining them."
        )
    if a.genome_build != b.genome_build:
        raise ConfigMismatchError(
            f"Cannot combine LDSC results{prefix}: genome_build mismatch "
            f"{a.genome_build!r} vs {b.genome_build!r}. Most likely the inputs were produced "
            "under different coordinate-build assumptions. Re-run or reload the inputs with the "
            "same genome_build before combining them."
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
            raise LDSCConfigError(
                _invalid_choice_message("AnnotationBuildConfig", "compression", self.compression, "'auto', 'gzip', 'bz2', or 'none'")
            )


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
    use_hm3_ref_panel_snps : bool, optional
        If ``True``, restrict the runtime reference-panel universe to the
        packaged curated HM3 SNP map. Mutually exclusive with
        ``ref_panel_snps_file``. Default is ``False``.
    exclude_regions : tuple of str, optional
        Named region presets to exclude (e.g. ``"mhc"``, ``"centromeres"``).
        Requires ``exclude_regions_build``. Default is ``()``.
    exclude_regions_bed : tuple of str, optional
        Paths to user-supplied BED files whose intervals are excluded.
        Build-agnostic; may be combined with ``exclude_regions``. Default is ``()``.
    exclude_regions_build : {"hg19", "hg38"} or None, optional
        Genome build used to resolve preset intervals. Required when
        ``exclude_regions`` is non-empty. Default is ``None``.
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
    use_hm3_ref_panel_snps: bool = False
    exclude_regions: tuple[str, ...] = ()
    exclude_regions_bed: tuple[str, ...] = ()
    exclude_regions_build: Literal["hg19", "hg38"] | None = None

    def __post_init__(self) -> None:
        """Normalize backend path tokens and validate parquet-R2 settings."""
        if self.backend not in {"auto", "plink", "parquet_r2"}:
            raise LDSCConfigError(_invalid_choice_message("RefPanelConfig", "backend", self.backend, "'auto', 'plink', or 'parquet_r2'"))
        if self.r2_bias_mode not in {None, "raw", "unbiased"}:
            raise LDSCConfigError(_invalid_choice_message("RefPanelConfig", "r2_bias_mode", self.r2_bias_mode, "None, 'raw', or 'unbiased'"))
        if self.r2_sample_size is not None and self.r2_sample_size <= 0:
            raise LDSCConfigError(_positive_number_message("RefPanelConfig", "r2_sample_size", self.r2_sample_size))
        if self.sample_size is not None and self.sample_size <= 0:
            raise LDSCConfigError(_positive_number_message("RefPanelConfig", "sample_size", self.sample_size))
        if self.maf_min is not None and not 0 <= self.maf_min <= 0.5:
            raise LDSCConfigError(_range_message("RefPanelConfig", "maf_min", self.maf_min, "[0, 0.5]"))
        if self.sample_size is None and self.r2_sample_size is not None:
            object.__setattr__(self, "sample_size", int(self.r2_sample_size))
        object.__setattr__(self, "plink_prefix", _normalize_optional_path(self.plink_prefix))
        object.__setattr__(self, "r2_dir", _normalize_optional_path(self.r2_dir))
        object.__setattr__(self, "keep_indivs_file", _normalize_optional_path(self.keep_indivs_file))
        object.__setattr__(self, "ref_panel_snps_file", _normalize_optional_path(self.ref_panel_snps_file))
        if self.ref_panel_snps_file is not None and self.use_hm3_ref_panel_snps:
            raise LDSCConfigError(_mutually_exclusive_message("RefPanelConfig", "ref_panel_snps_file", "use_hm3_ref_panel_snps"))
        if self.chromosomes is not None:
            from .chromosome_inference import normalize_chromosome

            object.__setattr__(self, "chromosomes", tuple(normalize_chromosome(chrom) for chrom in self.chromosomes))
        object.__setattr__(self, "exclude_regions", tuple(self.exclude_regions))
        object.__setattr__(
            self,
            "exclude_regions_bed",
            tuple(_normalize_optional_path(token) for token in self.exclude_regions_bed if token),
        )
        _validate_region_presets("RefPanelConfig", self.exclude_regions)
        if self.exclude_regions_build not in {None, "hg19", "hg38"}:
            raise LDSCConfigError(
                _invalid_choice_message(
                    "RefPanelConfig", "exclude_regions_build", self.exclude_regions_build, "None, 'hg19', or 'hg38'"
                )
            )
        if self.exclude_regions and self.exclude_regions_build is None:
            raise LDSCConfigError(
                "RefPanelConfig received exclude_regions presets without exclude_regions_build. "
                "Most likely --exclude-regions was passed without --exclude-regions-build. "
                "Region presets are build-specific; pass --exclude-regions-build hg19 or hg38."
            )


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
    use_hm3_regression_snps : bool, optional
        If ``True``, use the packaged curated HM3 SNP map as the persisted
        regression SNP set. Mutually exclusive with ``regression_snps_file``.
        Default is ``False``.
    snp_batch_size : int, optional
        Number of SNPs processed per LD-score sliding batch. Default is
        ``128``.
    common_maf_min : float, optional
        Inclusive MAF threshold used only for common-SNP count vectors
        (``MAF >= common_maf_min``). It does not change retained reference
        SNPs, LD-score rows, LD scores, or persisted regression-universe LD
        scores. Default is ``0.05``.
    whole_chromosome_ok : bool, optional
        Override the guard that rejects windows effectively spanning an entire
        chromosome. Default is ``False``.
    threads : int, optional
        Number of worker processes for cross-chromosome parallelism on one
        machine, using the scikit-learn/joblib ``n_jobs`` convention: ``1``
        (default) runs sequentially in-process; ``-1`` uses all available cores;
        ``-2`` uses all but one; any other negative ``-k`` uses
        ``n_cpus + 1 - k``. Core counts respect CPU affinity (SLURM/cgroup
        allocations), not the raw machine size. The effective count is capped at
        the chromosome count, and ``0`` is rejected. Output is identical
        regardless of this value.
    """
    ld_wind_snps: int | None = None
    ld_wind_kb: float | None = None
    ld_wind_cm: float | None = None
    regression_snps_file: str | PathLike[str] | None = None
    use_hm3_regression_snps: bool = False
    snp_batch_size: int = 128
    common_maf_min: float = 0.05
    whole_chromosome_ok: bool = False
    threads: int = 1

    def __post_init__(self) -> None:
        """Validate LD-window settings and normalize optional SNP-list paths."""
        windows = [self.ld_wind_snps, self.ld_wind_kb, self.ld_wind_cm]
        if sum(value is not None for value in windows) != 1:
            raise LDSCConfigError(_one_ld_window_message("LDScoreConfig"))
        if self.ld_wind_snps is not None and self.ld_wind_snps <= 0:
            raise LDSCConfigError(_positive_number_message("LDScoreConfig", "ld_wind_snps", self.ld_wind_snps))
        if self.ld_wind_kb is not None and self.ld_wind_kb <= 0:
            raise LDSCConfigError(_positive_number_message("LDScoreConfig", "ld_wind_kb", self.ld_wind_kb))
        if self.ld_wind_cm is not None and self.ld_wind_cm <= 0:
            raise LDSCConfigError(_positive_number_message("LDScoreConfig", "ld_wind_cm", self.ld_wind_cm))
        if not 0 <= self.common_maf_min <= 0.5:
            raise LDSCConfigError(_range_message("LDScoreConfig", "common_maf_min", self.common_maf_min, "[0, 0.5]"))
        if self.snp_batch_size <= 0:
            raise LDSCConfigError(_positive_number_message("LDScoreConfig", "snp_batch_size", self.snp_batch_size))
        if self.threads == 0:
            raise LDSCConfigError(
                "LDScoreConfig received threads=0, which is ambiguous. Most likely the "
                "thread count was left unset to a sentinel. Use 1 for sequential, a "
                "positive count for that many workers, -1 for all cores, or -2 for all "
                "but one."
            )
        object.__setattr__(self, "regression_snps_file", _normalize_optional_path(self.regression_snps_file))
        if self.regression_snps_file is not None and self.use_hm3_regression_snps:
            raise LDSCConfigError(_mutually_exclusive_message("LDScoreConfig", "regression_snps_file", "use_hm3_regression_snps"))


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
    source_genome_build : {"auto", "hg19", "hg37", "GRCh37", "hg38", "GRCh38"}, optional
        Genome build of the input PLINK coordinates. If ``"auto"``, the
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
        builder emits a source-build-only panel. Matching chain files are valid
        only when the active ``GlobalConfig.snp_identifier`` is in the
        ``chr_pos`` family.
    ld_wind_snps, ld_wind_kb, ld_wind_cm : int, float, float, or None, optional
        LD window specification. Exactly one must be supplied.
    maf_min : float or None, optional
        Optional retained-SNP MAF threshold. Default is ``None``.
    ref_panel_snps_file : str or os.PathLike[str] or None, optional
        Optional SNP list restricting the emitted reference-panel universe.
        Default is ``None``. The identifier mode comes from
        ``GlobalConfig.snp_identifier``. Restriction files may omit alleles and
        then match by base key. In ``chr_pos``-family modes, the restriction
        file must be aligned to the resolved source reference-panel build; the
        builder never uses ``GlobalConfig.genome_build``.
    use_hm3_snps : bool, optional
        If ``True``, restrict the emitted reference-panel universe to the
        packaged curated HM3 SNP map. Mutually exclusive with
        ``ref_panel_snps_file``. Default is ``False``.
    use_hm3_quick_liftover : bool, optional
        If ``True``, emit the opposite-build reference-panel artifacts for the
        HM3-restricted coordinate universe using the packaged curated HM3 map.
        Requires ``use_hm3_snps``, is valid only in ``chr_pos``-family modes,
        and is mutually exclusive with chain-file liftover. Default is
        ``False``.
    keep_indivs_file : str or os.PathLike[str] or None, optional
        Optional individual keep-file applied before R2 calculation. Default is
        ``None``.
    snp_batch_size : int, optional
        Number of SNPs loaded per pairwise-R2 computation batch. Larger values
        may improve throughput but use more memory. Default is ``128``.
    min_r2 : float, optional
        Opt-in unbiased-R2 floor for emitted pairs. ``0.0`` (the default, and
        any non-positive value) writes every retained pair, preserving exact
        pairwise completeness. A positive value drops pairs whose unbiased R2
        is below the threshold, reducing output size and memory at the cost of
        completeness (the read path treats absent pairs as R2=0). The applied
        value is recorded as ``ldsc:min_r2`` parquet metadata. Default is
        ``0.0``.
    overwrite : bool, optional
        If ``True``, replace current panel artifacts and remove stale
        workflow-owned parquet, metadata, dropped-SNP, or log siblings after a
        successful run. If ``False``, output collisions raise before
        chromosome processing starts. Default is ``False``.
    exclude_regions : tuple of str, optional
        Named region presets to exclude. Build is resolved from
        ``source_genome_build``. Default is ``()``.
    exclude_regions_bed : tuple of str, optional
        Paths to user-supplied BED files whose intervals are excluded.
        Default is ``()``.
    """

    plink_prefix: str | PathLike[str]
    source_genome_build: GenomeBuildInput = "auto"
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
    use_hm3_snps: bool = False
    use_hm3_quick_liftover: bool = False
    keep_indivs_file: str | PathLike[str] | None = None
    snp_batch_size: int = 128
    min_r2: float = 0.0
    overwrite: bool = False
    exclude_regions: tuple[str, ...] = ()
    exclude_regions_bed: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        """Normalize build paths and validate liftover and LD-window settings."""
        object.__setattr__(self, "plink_prefix", _normalize_required_path(self.plink_prefix))
        try:
            source_genome_build = normalize_genome_build(self.source_genome_build)
        except ValueError as exc:
            raise LDSCConfigError(
                f"Could not construct ReferencePanelBuildConfig: source_genome_build={self.source_genome_build!r} is invalid. "
                "Most likely the genome build name is misspelled or unsupported. "
                "Use 'auto', 'hg19', 'hg38', 'GRCh37', or 'GRCh38'."
            ) from exc
        object.__setattr__(self, "source_genome_build", source_genome_build)
        if self.source_genome_build is None:
            raise LDSCConfigError(
                "Could not construct ReferencePanelBuildConfig: source_genome_build=None is invalid. "
                "Most likely the build was omitted for a coordinate reference-panel build. "
                "Use source_genome_build='auto', 'hg19', or 'hg38'."
            )
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
        if self.ref_panel_snps_file is not None and self.use_hm3_snps:
            raise LDSCConfigError(_mutually_exclusive_message("ReferencePanelBuildConfig", "ref_panel_snps_file", "use_hm3_snps"))
        if self.use_hm3_quick_liftover and not self.use_hm3_snps:
            raise LDSCConfigError(
                "Could not construct ReferencePanelBuildConfig: use_hm3_quick_liftover=True requires use_hm3_snps=True. "
                "Most likely quick liftover was enabled without selecting the packaged HM3 SNP set. "
                "Enable use_hm3_snps or use a chain-file liftover option instead."
            )
        if self.use_hm3_quick_liftover and (
            self.liftover_chain_hg19_to_hg38_file is not None
            or self.liftover_chain_hg38_to_hg19_file is not None
        ):
            raise LDSCConfigError(
                "Could not construct ReferencePanelBuildConfig: chain-file liftover and use_hm3_quick_liftover are mutually exclusive. "
                "Most likely two liftover methods were selected for the same reference-panel build. "
                "Choose either the HM3 quick liftover shortcut or chain-file liftover, not both."
            )
        windows = [self.ld_wind_snps, self.ld_wind_kb, self.ld_wind_cm]
        if sum(value is not None for value in windows) != 1:
            raise LDSCConfigError(_one_ld_window_message("ReferencePanelBuildConfig"))
        source_map = None
        if self.source_genome_build == "hg19":
            source_map = self.genetic_map_hg19_sources
        elif self.source_genome_build == "hg38":
            source_map = self.genetic_map_hg38_sources
        if self.ld_wind_cm is not None and self.source_genome_build is not None and source_map is None:
            raise LDSCConfigError(
                f"Could not construct ReferencePanelBuildConfig: ld_wind_cm requires a {self.source_genome_build} genetic map path. "
                "Most likely cM-window LD scores were requested without the matching genetic map. "
                f"Set genetic_map_{self.source_genome_build}_sources or use ld_wind_snps/ld_wind_kb."
            )
        if self.ld_wind_snps is not None and self.ld_wind_snps <= 0:
            raise LDSCConfigError(_positive_number_message("ReferencePanelBuildConfig", "ld_wind_snps", self.ld_wind_snps))
        if self.ld_wind_kb is not None and self.ld_wind_kb <= 0:
            raise LDSCConfigError(_positive_number_message("ReferencePanelBuildConfig", "ld_wind_kb", self.ld_wind_kb))
        if self.ld_wind_cm is not None and self.ld_wind_cm <= 0:
            raise LDSCConfigError(_positive_number_message("ReferencePanelBuildConfig", "ld_wind_cm", self.ld_wind_cm))
        if self.maf_min is not None and not 0 <= self.maf_min <= 0.5:
            raise LDSCConfigError(_range_message("ReferencePanelBuildConfig", "maf_min", self.maf_min, "[0, 0.5]"))
        if self.snp_batch_size <= 0:
            raise LDSCConfigError(_positive_number_message("ReferencePanelBuildConfig", "snp_batch_size", self.snp_batch_size))
        if self.min_r2 > 1:
            raise LDSCConfigError(
                f"Could not construct ReferencePanelBuildConfig: min_r2={self.min_r2!r} must not exceed 1. "
                "Most likely an R2 percentage or invalid threshold was supplied. "
                "Use a decimal threshold in [0, 1]."
            )
        object.__setattr__(self, "exclude_regions", tuple(self.exclude_regions))
        object.__setattr__(
            self,
            "exclude_regions_bed",
            tuple(_normalize_optional_path(token) for token in self.exclude_regions_bed if token),
        )
        _validate_region_presets("ReferencePanelBuildConfig", self.exclude_regions)


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
        ``sumstats.sumstats.gz``, root ``metadata.json``, and
        ``diagnostics/sumstats.log`` artifacts.
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
        ``chr_pos``-family modes, central ``CHR``/``POS`` aliases, including
        build-specific position aliases such as ``hg19_POS`` and ``hg38_POS``,
        define retained coordinates. Restriction files may omit alleles even in
        allele-aware modes; allele-free restrictions match by base key before
        later artifact cleanup. This option restricts rows only; it does not
        rewrite alleles or reorder output. Default is ``None``.
    use_hm3_snps : bool, optional
        If ``True``, restrict summary-statistics rows to the packaged curated
        HM3 SNP map. Mutually exclusive with ``sumstats_snps_file``. Default is
        ``False``.
    source_genome_build : {"auto", "hg19", "hg37", "GRCh37", "hg38", "GRCh38"}, optional
        Genome build of raw ``CHR``/``POS`` coordinates. ``"auto"`` asks the
        munger to infer the source build from the raw file. Default is
        ``"auto"``.
    output_genome_build : {"hg19", "hg37", "GRCh37", "hg38", "GRCh38"} or None, optional
        Desired output coordinate build for ``chr_pos``-family munging. When
        it differs from the resolved source build, exactly one liftover method
        is required. The workflow requires this field in coordinate-family
        modes and rejects it in rsID-family modes. Default is ``None``.
    liftover_chain_file : str or os.PathLike[str] or None, optional
        Chain file used to convert ``CHR``/``POS`` from the resolved source
        build to ``output_genome_build``. Mutually exclusive with
        ``use_hm3_quick_liftover``. Sumstats liftover is valid only in
        ``chr_pos``-family modes; it updates coordinates and does not rewrite
        ``SNP``. Default is ``None``.
    use_hm3_quick_liftover : bool, optional
        If ``True``, use the packaged curated dual-build HM3 map for a
        coordinate-only quick liftover after HM3 SNP restriction. Requires
        ``use_hm3_snps`` and is mutually exclusive with ``liftover_chain_file``.
        If the resolved source build already equals ``output_genome_build``,
        the workflow warns and ignores this flag.
        Default is ``False``.
    signed_sumstats_spec : str or None, optional
        Signed statistic specification passed through to the legacy kernel.
        Default is ``None``.
    ignore_columns : tuple of str, optional
        Source columns ignored during auto-detection. Default is ``()``.
    sumstats_format : {"auto", "plain", "daner-old", "daner-new"}, optional
        Raw summary-statistics format profile. ``"auto"`` detects common
        formats from headers and leading metadata. Default is ``"auto"``.
    info_list_columns : tuple of str, optional
        INFO-like columns with comma-separated per-study values. Values are
        summarized before INFO filtering by taking the mean of numeric,
        non-missing tokens, for example ``IMPINFO=0.852,0.113,NA``. Mixed
        nonnumeric tokens are rejected with a repair suggestion. Default is
        ``()``.
    a1_inc, keep_maf, daner_old, daner_new : bool, optional
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
    use_hm3_snps: bool = False
    source_genome_build: GenomeBuildInput = "auto"
    output_genome_build: GenomeBuildInput | None = None
    liftover_chain_file: str | PathLike[str] | None = None
    use_hm3_quick_liftover: bool = False
    signed_sumstats_spec: str | None = None
    ignore_columns: tuple[str, ...] = field(default_factory=tuple)
    info_list_columns: tuple[str, ...] = field(default_factory=tuple)
    sumstats_format: str = "auto"
    a1_inc: bool = False
    keep_maf: bool = False
    daner_old: bool = False
    daner_new: bool = False
    overwrite: bool = False

    def __post_init__(self) -> None:
        """Validate munging thresholds and normalize optional file paths."""
        if self.info_min < 0:
            raise LDSCConfigError(
                f"Could not construct MungeConfig: info_min={self.info_min!r} must be non-negative. "
                "Most likely a negative INFO threshold was supplied. Set info_min to 0 or higher."
            )
        if self.maf_min < 0 or self.maf_min > 0.5:
            raise LDSCConfigError(_range_message("MungeConfig", "maf_min", self.maf_min, "[0, 0.5]"))
        if self.chunk_size <= 0:
            raise LDSCConfigError(_positive_number_message("MungeConfig", "chunk_size", self.chunk_size))
        if self.output_format not in {"parquet", "tsv.gz", "both"}:
            raise LDSCConfigError(_invalid_choice_message("MungeConfig", "output_format", self.output_format, "'parquet', 'tsv.gz', or 'both'"))
        if self.sumstats_format not in {"auto", "plain", "daner-old", "daner-new"}:
            raise LDSCConfigError(_invalid_choice_message("MungeConfig", "sumstats_format", self.sumstats_format, "'auto', 'plain', 'daner-old', or 'daner-new'"))
        try:
            source_genome_build = normalize_genome_build(self.source_genome_build)
        except ValueError as exc:
            raise LDSCConfigError(
                f"Could not construct MungeConfig: source_genome_build={self.source_genome_build!r} is invalid. "
                "Most likely the genome build name is misspelled or unsupported. "
                "Use 'auto', 'hg19', 'hg38', 'GRCh37', or 'GRCh38'."
            ) from exc
        if source_genome_build is None:
            raise LDSCConfigError(
                "Could not construct MungeConfig: source_genome_build=None is invalid. "
                "Most likely the raw coordinate build was omitted. "
                "Use source_genome_build='auto', 'hg19', or 'hg38'."
            )
        try:
            output_genome_build = normalize_genome_build(self.output_genome_build)
        except ValueError as exc:
            raise LDSCConfigError(
                f"Could not construct MungeConfig: output_genome_build={self.output_genome_build!r} is invalid. "
                "Most likely the output genome build name is misspelled or unsupported. "
                "Use 'hg19', 'hg38', 'GRCh37', or 'GRCh38', or omit it for rsID workflows."
            ) from exc
        if output_genome_build == "auto":
            raise LDSCConfigError(
                "Could not construct MungeConfig: output_genome_build='auto' is invalid. "
                "Most likely an inference setting was used for an output coordinate build. "
                "Choose a concrete output_genome_build of 'hg19' or 'hg38'."
            )
        object.__setattr__(self, "output_dir", _normalize_optional_path(self.output_dir))
        object.__setattr__(self, "raw_sumstats_file", _normalize_optional_path(self.raw_sumstats_file))
        object.__setattr__(self, "sumstats_snps_file", _normalize_optional_path(self.sumstats_snps_file))
        object.__setattr__(self, "source_genome_build", source_genome_build)
        object.__setattr__(self, "output_genome_build", output_genome_build)
        object.__setattr__(self, "liftover_chain_file", _normalize_optional_path(self.liftover_chain_file))
        object.__setattr__(self, "trait_name", _normalize_trait_name(self.trait_name))
        object.__setattr__(self, "ignore_columns", tuple(self.ignore_columns))
        object.__setattr__(self, "info_list_columns", tuple(self.info_list_columns))
        object.__setattr__(self, "column_hints", dict(self.column_hints))
        if self.sumstats_snps_file is not None and self.use_hm3_snps:
            raise LDSCConfigError(_mutually_exclusive_message("MungeConfig", "sumstats_snps_file", "use_hm3_snps"))
        if self.use_hm3_quick_liftover and not self.use_hm3_snps:
            raise LDSCConfigError(
                "Could not construct MungeConfig: use_hm3_quick_liftover=True requires use_hm3_snps=True. "
                "Most likely quick liftover was enabled without selecting the packaged HM3 SNP set. "
                "Enable use_hm3_snps or use a chain-file liftover option instead."
            )
        if self.liftover_chain_file is not None and self.use_hm3_quick_liftover:
            raise LDSCConfigError(
                "Could not construct MungeConfig: liftover_chain_file and use_hm3_quick_liftover are mutually exclusive. "
                "Most likely two liftover methods were selected for the same munge-sumstats run. "
                "Choose either a chain file or HM3 quick liftover, not both."
            )


@dataclass(frozen=True)
class RegressionConfig:
    """Configuration for heritability and genetic-correlation regressions.

    Parameters
    ----------
    n_blocks : int, optional
        Requested number of block-jackknife partitions. Default is ``200``.
    use_common_counts : bool, optional
        If ``True``, prefer the LD-score metadata ``common_reference_snp_counts`` vector
        when it is available; otherwise use ``all_reference_snp_counts``.
        Default is ``True``.
    use_intercept : bool, optional
        If ``False``, constrain the intercept to the LDSC default for the model
        being fit. Default is ``True``.
    intercept_h2, intercept_gencov : float or None, optional
        Fixed intercept values for single-trait and cross-trait models. In
        multi-trait rg runs, scalar values are broadcast to every pair.
        Defaults are ``None``.
    two_step_cutoff : float or None, optional
        Threshold for the two-step estimator used by the regression kernel.
        Default is ``None``.
    chisq_max : float or None, optional
        Maximum allowed chi-square statistic before a row is filtered. Default is
        ``None``.
    samp_prev, pop_prev : float, list of float, or None, optional
        Liability-scale prevalence inputs. Defaults are ``None``.
    allow_identity_downgrade : bool, optional
        If ``True``, same-family allele-aware/base regression inputs may run
        under the base identity mode. Cross-family mixes remain rejected.
        Original modes and dropped duplicate counts are logged by regression
        workflows, not persisted in detail metadata. Default is ``False``.
    """
    n_blocks: int = 200
    use_common_counts: bool = True
    use_intercept: bool = True
    intercept_h2: float | None = None
    intercept_gencov: float | None = None
    two_step_cutoff: float | None = None
    chisq_max: float | None = None
    samp_prev: float | list[float] | None = None
    pop_prev: float | list[float] | None = None
    allow_identity_downgrade: bool = False

    def __post_init__(self) -> None:
        """Validate regression hyperparameters after dataclass construction."""
        if self.n_blocks <= 1:
            raise LDSCConfigError(
                f"Could not construct RegressionConfig: n_blocks={self.n_blocks!r} must be greater than 1. "
                "Most likely too few jackknife blocks were requested. Set n_blocks to at least 2."
            )
        if self.two_step_cutoff is not None and self.two_step_cutoff <= 0:
            raise LDSCConfigError(_positive_number_message("RegressionConfig", "two_step_cutoff", self.two_step_cutoff))
        if self.chisq_max is not None and self.chisq_max <= 0:
            raise LDSCConfigError(_positive_number_message("RegressionConfig", "chisq_max", self.chisq_max))
