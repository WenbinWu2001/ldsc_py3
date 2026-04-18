"""Reference-panel specifications and workflow-to-kernel backend adapters.

The public reference-panel specs accept normalized token strings, including
glob patterns and chromosome-suite prefixes, but this module resolves those
tokens to concrete primitive strings before the execution kernel sees them.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from os import PathLike
from pathlib import Path
from typing import Any

import pandas as pd

from ..config import CommonConfig, RefPanelConfig, _normalize_optional_path, _normalize_path_tuple
from ..path_resolution import (
    FREQUENCY_SUFFIXES,
    PARQUET_SUFFIXES,
    resolve_chromosome_group,
    resolve_file_group,
    resolve_plink_prefix,
    resolve_plink_prefix_group,
)
from . import formats as legacy_parse
from . import ldscore as kernel_ldscore
from .identifiers import (
    build_snp_id_series,
    infer_chr_pos_columns,
    infer_snp_column,
    normalize_chromosome,
    normalize_snp_identifier_mode,
    read_global_snp_restriction,
    validate_unique_snp_ids,
)


@dataclass(frozen=True)
class RefPanelSpec:
    """Immutable description of a reference-panel backend and its sources.

    The public spec accepts path-like inputs and normalizes them to primitive
    strings. Resolution is deferred until workflow execution so callers may use
    exact paths, standard Python globs, explicit ``@`` chromosome-suite
    placeholders, or legacy bare prefixes in suite-capable fields.
    """
    backend: str
    bfile_prefix: str | PathLike[str] | None = None
    r2_table_paths: str | PathLike[str] | tuple[str | PathLike[str], ...] | list[str | PathLike[str]] = field(default_factory=tuple)
    sample_size: int | None = None
    maf_metadata_paths: str | PathLike[str] | tuple[str | PathLike[str], ...] | list[str | PathLike[str]] = field(default_factory=tuple)
    chromosomes: tuple[str, ...] | None = None
    genome_build: str | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "backend", self.backend)
        object.__setattr__(self, "bfile_prefix", _normalize_optional_path(self.bfile_prefix))
        object.__setattr__(self, "r2_table_paths", _normalize_path_tuple(self.r2_table_paths))
        object.__setattr__(self, "maf_metadata_paths", _normalize_path_tuple(self.maf_metadata_paths))
        if self.chromosomes is not None:
            object.__setattr__(self, "chromosomes", tuple(normalize_chromosome(chrom) for chrom in self.chromosomes))


class RefPanel(ABC):
    """Abstract chromosome-scoped reference-panel interface."""
    def __init__(self, common_config: CommonConfig, spec: RefPanelSpec) -> None:
        self.common_config = common_config
        self.spec = spec
        self._metadata_cache: dict[str, pd.DataFrame] = {}

    @abstractmethod
    def available_chromosomes(self) -> list[str]:
        """Return chromosomes available from the backend."""
        raise NotImplementedError

    @abstractmethod
    def load_metadata(self, chrom: str) -> pd.DataFrame:
        """Load SNP metadata for ``chrom``."""
        raise NotImplementedError

    @abstractmethod
    def build_reader(self, chrom: str, **kwargs: Any) -> Any:
        """Build a backend-specific reader for ``chrom``."""
        raise NotImplementedError

    def filter_to_snps(self, chrom: str, snps: set[str] | list[str]) -> pd.DataFrame:
        """Subset chromosome metadata to the requested SNP identifiers."""
        metadata = self.load_metadata(chrom)
        keep = build_snp_id_series(metadata, self.common_config.snp_identifier).isin(set(snps))
        return metadata.loc[keep].reset_index(drop=True)

    def summary(self) -> dict[str, Any]:
        """Summarize backend type and source paths."""
        return {
            "backend": self.spec.backend,
            "chromosomes": self.available_chromosomes(),
            "source": {
                "bfile_prefix": self.spec.bfile_prefix,
                "r2_table_paths": list(self.spec.r2_table_paths),
                "maf_metadata_paths": list(self.spec.maf_metadata_paths),
            },
        }

    def _apply_global_restriction(self, metadata: pd.DataFrame) -> pd.DataFrame:
        restrict_path = self.common_config.global_snp_restriction_path
        if restrict_path is None or len(metadata) == 0:
            return metadata
        restrict_ids = read_global_snp_restriction(restrict_path, self.common_config.snp_identifier)
        keep = build_snp_id_series(metadata, self.common_config.snp_identifier).isin(restrict_ids)
        return metadata.loc[keep].reset_index(drop=True)

    def _validate_metadata(self, metadata: pd.DataFrame, chrom: str) -> pd.DataFrame:
        metadata = metadata.reset_index(drop=True)
        validate_unique_snp_ids(metadata, self.common_config.snp_identifier, context=f"{type(self).__name__}[{chrom}]")
        return metadata


class PlinkRefPanel(RefPanel):
    """PLINK-backed reference-panel adapter."""
    def available_chromosomes(self) -> list[str]:
        df = self._read_bim_table(chrom=None)
        chromosomes = sorted(df["CHR"].astype(str).map(normalize_chromosome).unique().tolist(), key=_chrom_sort_key)
        return chromosomes

    def load_metadata(self, chrom: str) -> pd.DataFrame:
        chrom = normalize_chromosome(chrom)
        if chrom in self._metadata_cache:
            return self._metadata_cache[chrom].copy()

        df = self._read_bim_table(chrom=chrom)
        df["CHR"] = df["CHR"].map(normalize_chromosome)
        df["SNP"] = df["SNP"].astype(str)
        df["CM"] = pd.to_numeric(df["CM"], errors="coerce")
        df["POS"] = pd.to_numeric(df["POS"], errors="raise").astype(int)
        metadata = df.loc[df["CHR"] == chrom, ["CHR", "SNP", "CM", "POS"]].reset_index(drop=True)
        if len(metadata) == 0:
            raise ValueError(f"No PLINK metadata rows found for chromosome {chrom}.")
        metadata = self._apply_global_restriction(metadata)
        metadata = self._validate_metadata(metadata, chrom)
        self._metadata_cache[chrom] = metadata.copy()
        return metadata

    def build_reader(self, chrom: str, keep_snps: set[str] | list[str] | None = None, keep_indivs: list[int] | None = None, maf_min: float | None = None):
        prefix = None if self.spec.bfile_prefix is None else resolve_plink_prefix(self.spec.bfile_prefix, chrom=chrom)
        if prefix is None:
            raise ValueError("PlinkRefPanel requires bfile_prefix.")
        bim = legacy_parse.PlinkBIMFile(prefix + ".bim")
        fam = legacy_parse.PlinkFAMFile(prefix + ".fam")

        keep_snp_indices = None
        if keep_snps is not None:
            metadata = self.load_metadata(chrom)
            kept = self.filter_to_snps(chrom, set(keep_snps))
            key_to_index = dict(zip(build_snp_id_series(metadata, self.common_config.snp_identifier), metadata.index))
            keep_snp_indices = [int(key_to_index[key]) for key in build_snp_id_series(kept, self.common_config.snp_identifier)]

        return kernel_ldscore.PlinkBEDFile(
            prefix + ".bed",
            len(fam.IDList),
            bim,
            keep_snps=keep_snp_indices,
            keep_indivs=keep_indivs,
            mafMin=maf_min,
        )

    def _read_bim_table(self, chrom: str | None) -> pd.DataFrame:
        prefixes = (
            []
            if self.spec.bfile_prefix is None
            else resolve_plink_prefix_group(
                (self.spec.bfile_prefix,),
                chrom=chrom,
                allow_chromosome_suite=(chrom is None),
            )
        )
        if not prefixes:
            raise ValueError("PlinkRefPanel requires bfile_prefix.")
        frames = [
            pd.read_csv(
                prefix + ".bim",
                sep=r"\s+",
                header=None,
                names=["CHR", "SNP", "CM", "POS", "A1", "A2"],
                usecols=[0, 1, 2, 3],
            )
            for prefix in prefixes
        ]
        return pd.concat(frames, axis=0, ignore_index=True)


class ParquetR2RefPanel(RefPanel):
    """Sorted-parquet-R2-backed reference-panel adapter."""
    def available_chromosomes(self) -> list[str]:
        if self.spec.chromosomes is not None:
            return sorted(set(self.spec.chromosomes), key=_chrom_sort_key)

        chromosomes: set[str] = set()
        for path in resolve_file_group(
            self.spec.maf_metadata_paths,
            suffixes=FREQUENCY_SUFFIXES,
            label="reference metadata",
            allow_chromosome_suite=True,
        ):
            df = _read_metadata_table(path, chrom=None, common_config=self.common_config)
            chromosomes.update(df["CHR"].map(normalize_chromosome).unique().tolist())
        if chromosomes:
            return sorted(chromosomes, key=_chrom_sort_key)
        raise ImportError("Chromosome discovery for parquet R2 requires explicit chromosomes or sidecar metadata in this environment.")

    def load_metadata(self, chrom: str) -> pd.DataFrame:
        chrom = normalize_chromosome(chrom)
        if chrom in self._metadata_cache:
            return self._metadata_cache[chrom].copy()
        if not self.spec.maf_metadata_paths:
            raise ImportError("ParquetR2RefPanel.load_metadata requires metadata sidecar files.")

        resolved_paths = resolve_chromosome_group(
            self.spec.maf_metadata_paths,
            chrom=chrom,
            suffixes=FREQUENCY_SUFFIXES,
            label="reference metadata",
            required=False,
        )
        frames = [_read_metadata_table(path, chrom=chrom, common_config=self.common_config) for path in resolved_paths]
        frames = [frame for frame in frames if len(frame) > 0]
        if not frames:
            raise ValueError(f"No parquet metadata rows found for chromosome {chrom}.")
        metadata = pd.concat(frames, axis=0, ignore_index=True)
        metadata = self._apply_global_restriction(metadata)
        metadata = self._validate_metadata(metadata, chrom)
        self._metadata_cache[chrom] = metadata.copy()
        return metadata

    def build_reader(
        self,
        chrom: str,
        metadata: pd.DataFrame | None = None,
        r2_bias_mode: str | None = None,
        r2_sample_size: float | None = None,
    ):
        metadata = metadata if metadata is not None else self.load_metadata(chrom)
        return kernel_ldscore.SortedR2BlockReader(
            paths=resolve_chromosome_group(
                self.spec.r2_table_paths,
                chrom=chrom,
                suffixes=PARQUET_SUFFIXES,
                label="parquet R2",
                required=False,
            ),
            chrom=chrom,
            metadata=metadata,
            identifier_mode=self.common_config.snp_identifier,
            r2_bias_mode="unbiased" if r2_bias_mode is None else r2_bias_mode,
            r2_sample_size=r2_sample_size if r2_sample_size is not None else self.spec.sample_size,
            genome_build=self.spec.genome_build or self.common_config.genome_build,
        )


class RefPanelLoader:
    """Instantiate the backend requested by a :class:`RefPanelSpec`."""
    def __init__(self, common_config: CommonConfig, ref_panel_config: RefPanelConfig | None = None) -> None:
        self.common_config = common_config
        self.ref_panel_config = ref_panel_config or RefPanelConfig()

    def load(self, ref_panel_spec: RefPanelSpec) -> RefPanel:
        """Return a concrete reference-panel adapter for ``ref_panel_spec``."""
        backend = ref_panel_spec.backend
        if backend == "plink":
            return PlinkRefPanel(self.common_config, ref_panel_spec)
        if backend == "parquet_r2":
            return ParquetR2RefPanel(self.common_config, ref_panel_spec)
        raise ValueError(f"Unsupported reference-panel backend: {backend}")


def _chrom_sort_key(chrom: str) -> tuple[int, str]:
    """Return a stable chromosome sort key matching the package-wide ordering."""
    chrom = normalize_chromosome(chrom)
    if chrom.isdigit():
        return (int(chrom), chrom)
    special = {"X": 23, "Y": 24, "MT": 25, "M": 25}
    return (special.get(chrom, 99), chrom)


def _read_metadata_table(path: str | Path, chrom: str | None, common_config: CommonConfig) -> pd.DataFrame:
    """Read one reference-panel metadata table into the normalized column set."""
    df = pd.read_csv(path, sep=r"\s+", compression="gzip" if str(path).endswith(".gz") else None)
    snp_identifier = normalize_snp_identifier_mode(common_config.snp_identifier)

    chr_col = None
    bp_col = None
    snp_col = None
    try:
        chr_col, pos_col = infer_chr_pos_columns(df.columns)
    except ValueError:
        pass
    try:
        snp_col = infer_snp_column(df.columns)
    except ValueError:
        pass

    if snp_identifier == "rsid" and snp_col is None:
        raise ValueError(f"{path} must contain a SNP column in rsid mode.")
    if snp_identifier == "chr_pos" and (chr_col is None or pos_col is None):
        raise ValueError(f"{path} must contain CHR and POS columns in chr_pos mode.")

    out = pd.DataFrame(index=df.index)
    if chr_col is not None:
        out["CHR"] = df[chr_col].map(normalize_chromosome)
    if pos_col is not None:
        out["POS"] = pd.to_numeric(df[pos_col], errors="raise").astype(int)
    if snp_col is not None:
        out["SNP"] = df[snp_col].astype(str)

    cm_col = _find_optional_column(df.columns, ("CM", "CMBP", "CENTIMORGAN"))
    maf_col = _find_optional_column(df.columns, ("MAF", "FRQ", "FREQ", "FREQUENCY"))
    if cm_col is not None:
        out["CM"] = pd.to_numeric(df[cm_col], errors="coerce")
    else:
        out["CM"] = pd.NA
    if maf_col is not None:
        maf = pd.to_numeric(df[maf_col], errors="coerce").astype(float)
        out["MAF"] = pd.Series(maf).map(lambda value: value if pd.isna(value) else min(value, 1.0 - value))
    elif common_config.fail_on_missing_metadata:
        raise ValueError(f"{path} is missing MAF metadata.")

    if chrom is not None and "CHR" in out.columns:
        out = out.loc[out["CHR"] == normalize_chromosome(chrom)].reset_index(drop=True)
    return out


def _find_optional_column(columns, aliases: tuple[str, ...]) -> str | None:
    """Return the first optional column whose uppercased name matches ``aliases``."""
    normalized_aliases = {alias.upper() for alias in aliases}
    for column in columns:
        if column.upper() in normalized_aliases:
            return column
    return None
