"""Shared workflow configuration dataclasses for the refactor scaffold."""

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
    if path is None:
        return None
    return str(path)


def _normalize_path_tuple(values: tuple[str | Path, ...] | list[str | Path] | None) -> tuple[str, ...]:
    if not values:
        return ()
    return tuple(str(value) for value in values)


def _normalize_log_level(level: str) -> LogLevel:
    normalized = level.upper()
    allowed = {"DEBUG", "INFO", "WARNING", "ERROR"}
    if normalized not in allowed:
        raise ValueError(f"log_level must be one of {sorted(allowed)}; got {level!r}.")
    return normalized  # type: ignore[return-value]


@dataclass(frozen=True)
class CommonConfig:
    snp_identifier: SNPIdentifierMode = "rsid"
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
