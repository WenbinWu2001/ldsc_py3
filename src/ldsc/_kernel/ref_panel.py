"""Reference-panel configuration and workflow-to-kernel backend adapters.

This module resolves public reference-panel configuration into concrete
primitive strings before the execution kernel sees them.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
import logging
from pathlib import Path
import re
from typing import Any

import pandas as pd

from ..chromosome_inference import chrom_sort_key, normalize_chromosome
from ..column_inference import (
    REFERENCE_METADATA_SPEC_MAP,
    resolve_optional_column,
    resolve_required_column,
)
from ..config import GlobalConfig, RefPanelConfig
from ..genome_build_inference import validate_auto_genome_build_mode
from ..path_resolution import resolve_plink_prefix, resolve_plink_prefix_group
from . import formats as legacy_parse
from . import ldscore as kernel_ldscore
from .identifiers import (
    build_snp_id_series,
    normalize_snp_identifier_mode,
    read_global_snp_restriction,
    validate_unique_snp_ids,
)

LOGGER = logging.getLogger("LDSC.ref_panel")
_REF_PANEL_R2_RE = re.compile(r"^chr(?P<chrom>.+)_r2\.parquet$", flags=re.IGNORECASE)


class RefPanel(ABC):
    """Abstract chromosome-scoped reference-panel interface."""
    def __init__(self, global_config: GlobalConfig, spec: RefPanelConfig) -> None:
        """Store shared config and initialize the per-chromosome metadata cache."""
        validate_auto_genome_build_mode(global_config.snp_identifier, global_config.genome_build)
        self.global_config = global_config
        self.spec = spec
        self._metadata_cache: dict[str, pd.DataFrame] = {}

    @abstractmethod
    def available_chromosomes(self) -> list[str]:
        """Return chromosomes available from the backend."""
        raise NotImplementedError

    @abstractmethod
    def load_metadata(self, chrom: str) -> pd.DataFrame:
        """Load prepared SNP metadata for ``chrom`` after panel-layer filters."""
        raise NotImplementedError

    @abstractmethod
    def build_reader(self, chrom: str, **kwargs: Any) -> Any:
        """Build a backend-specific reader for ``chrom``."""
        raise NotImplementedError

    def filter_to_snps(self, chrom: str, snps: set[str] | list[str]) -> pd.DataFrame:
        """Subset chromosome metadata to the requested SNP identifiers."""
        metadata = self.load_metadata(chrom)
        keep = build_snp_id_series(metadata, self.global_config.snp_identifier).isin(set(snps))
        return metadata.loc[keep].reset_index(drop=True)

    def summary(self) -> dict[str, Any]:
        """Summarize backend type and source paths."""
        return {
            "backend": self.spec.backend,
            "chromosomes": self.available_chromosomes(),
            "source": {
                "plink_prefix": self.spec.plink_prefix,
                "r2_dir": self.spec.r2_dir,
            },
        }

    def _apply_maf_filter(self, metadata: pd.DataFrame, chrom: str) -> pd.DataFrame:
        """Apply the retained-panel MAF filter from ``RefPanelConfig`` when set."""
        maf_min = self.spec.maf_min
        if maf_min is None or len(metadata) == 0:
            return metadata.reset_index(drop=True)
        if "MAF" not in metadata.columns or metadata["MAF"].isna().all():
            LOGGER.warning(f"Cannot apply --maf-min on chromosome {chrom} because MAF metadata is unavailable.")
            return metadata.reset_index(drop=True)
        keep = metadata["MAF"] > maf_min
        removed = int((~keep).sum())
        if removed:
            LOGGER.info(f"Removed {removed} reference-panel SNPs with MAF <= {maf_min} on chromosome {chrom}.")
        return metadata.loc[keep].reset_index(drop=True)

    def _apply_snp_restriction(self, metadata: pd.DataFrame) -> pd.DataFrame:
        """Filter metadata rows to ``RefPanelConfig.ref_panel_snps_file`` when set."""
        restrict_path = self.spec.ref_panel_snps_file
        if restrict_path is None or len(metadata) == 0:
            return metadata
        restrict_ids = read_global_snp_restriction(
            restrict_path,
            self.global_config.snp_identifier,
            genome_build=self.global_config.genome_build,
            logger=LOGGER,
        )
        keep = build_snp_id_series(metadata, self.global_config.snp_identifier).isin(restrict_ids)
        return metadata.loc[keep].reset_index(drop=True)

    def _validate_metadata(self, metadata: pd.DataFrame, chrom: str) -> pd.DataFrame:
        """Reset row order and validate identifier uniqueness for one chromosome."""
        metadata = metadata.reset_index(drop=True)
        validate_unique_snp_ids(metadata, self.global_config.snp_identifier, context=f"{type(self).__name__}[{chrom}]")
        return metadata


class PlinkRefPanel(RefPanel):
    """PLINK-backed reference-panel adapter."""
    def available_chromosomes(self) -> list[str]:
        """List normalized chromosomes present in the resolved PLINK inputs."""
        df = self._read_bim_table(chrom=None)
        chromosomes = sorted(df["CHR"].astype(str).map(normalize_chromosome).unique().tolist(), key=_chrom_sort_key)
        return chromosomes

    def load_metadata(self, chrom: str) -> pd.DataFrame:
        """Load and cache normalized PLINK metadata for one chromosome."""
        chrom = normalize_chromosome(chrom)
        if chrom in self._metadata_cache:
            return self._metadata_cache[chrom].copy()

        df = self._read_bim_table(chrom=chrom).reset_index(drop=False).rename(columns={"index": "_raw_index"})
        df["CHR"] = df["CHR"].map(normalize_chromosome)
        df["SNP"] = df["SNP"].astype(str)
        df["CM"] = pd.to_numeric(df["CM"], errors="coerce")
        df["POS"] = pd.to_numeric(df["POS"], errors="raise").astype(int)
        metadata = df.loc[df["CHR"] == chrom, ["_raw_index", "CHR", "SNP", "CM", "POS"]].reset_index(drop=True)
        if len(metadata) == 0:
            raise ValueError(f"No PLINK metadata rows found for chromosome {chrom}.")
        metadata = self._apply_snp_restriction(metadata)
        validate_unique_snp_ids(metadata, self.global_config.snp_identifier, context=f"{type(self).__name__}[{chrom}]")
        metadata = self._load_genotype_metadata(chrom, metadata)
        metadata = self._validate_metadata(metadata, chrom)
        self._metadata_cache[chrom] = metadata.copy()
        return metadata

    def build_reader(self, chrom: str, keep_snps: set[str] | list[str] | None = None, keep_indivs: list[int] | None = None, maf_min: float | None = None):
        """Build a PLINK BED reader with optional SNP, sample, and MAF filters."""
        prefix = None if self.spec.plink_prefix is None else resolve_plink_prefix(self.spec.plink_prefix, chrom=chrom)
        if prefix is None:
            raise ValueError("PlinkRefPanel requires plink_prefix.")
        bim = legacy_parse.PlinkBIMFile(prefix + ".bim")
        fam = legacy_parse.PlinkFAMFile(prefix + ".fam")

        keep_snp_indices = None
        if keep_snps is not None:
            raw_metadata = self._read_bim_table(chrom=chrom)
            raw_metadata["CHR"] = raw_metadata["CHR"].map(normalize_chromosome)
            raw_metadata["SNP"] = raw_metadata["SNP"].astype(str)
            raw_metadata["POS"] = pd.to_numeric(raw_metadata["POS"], errors="raise").astype(int)
            raw_metadata = raw_metadata.loc[raw_metadata["CHR"] == normalize_chromosome(chrom)].reset_index(drop=True)
            kept = self.filter_to_snps(chrom, set(keep_snps))
            key_to_index = dict(zip(build_snp_id_series(raw_metadata, self.global_config.snp_identifier), raw_metadata.index))
            keep_snp_indices = [int(key_to_index[key]) for key in build_snp_id_series(kept, self.global_config.snp_identifier)]

        return kernel_ldscore.PlinkBEDFile(
            prefix + ".bed",
            len(fam.IDList),
            bim,
            keep_snps=keep_snp_indices,
            keep_indivs=keep_indivs,
            mafMin=maf_min,
        )

    def _load_genotype_metadata(self, chrom: str, metadata: pd.DataFrame) -> pd.DataFrame:
        """Return PLINK metadata with genotype-derived MAF when the reader is available."""
        prefix = None if self.spec.plink_prefix is None else resolve_plink_prefix(self.spec.plink_prefix, chrom=chrom)
        if prefix is None:
            raise ValueError("PlinkRefPanel requires plink_prefix.")
        bim = legacy_parse.PlinkBIMFile(prefix + ".bim")
        fam = legacy_parse.PlinkFAMFile(prefix + ".fam")
        keep_indivs = self._resolve_keep_individuals(fam)
        keep_snps = [int(index) for index in metadata["_raw_index"].tolist()]
        try:
            bed = kernel_ldscore.PlinkBEDFile(
                prefix + ".bed",
                len(fam.IDList),
                bim,
                keep_snps=keep_snps,
                keep_indivs=keep_indivs,
                mafMin=self.spec.maf_min,
            )
        except Exception as exc:
            if self.spec.maf_min is not None or self.spec.keep_indivs_file is not None:
                raise
            LOGGER.debug(f"Falling back to BIM-only PLINK metadata for chromosome {chrom}: {exc}", exc_info=True)
            return metadata.drop(columns="_raw_index", errors="ignore").reset_index(drop=True)

        out = pd.DataFrame(bed.df, columns=bed.colnames)
        if "BP" in out.columns:
            out = out.rename(columns={"BP": "POS"})
        out["CHR"] = out["CHR"].map(lambda value: normalize_chromosome(value, context=prefix + ".bim"))
        out["SNP"] = out["SNP"].astype(str)
        out["POS"] = pd.to_numeric(out["POS"], errors="raise").astype(int)
        out["CM"] = pd.to_numeric(out["CM"], errors="coerce")
        out["MAF"] = pd.to_numeric(out["MAF"], errors="coerce")
        return out.loc[out["CHR"] == normalize_chromosome(chrom)].reset_index(drop=True)

    def _resolve_keep_individuals(self, fam) -> list[int] | None:
        """Return FAM row indices retained by ``RefPanelConfig.keep_indivs_file``."""
        if self.spec.keep_indivs_file is None:
            return None
        keep_indivs = fam.loj(legacy_parse.FilterFile(self.spec.keep_indivs_file).IDList)
        if len(keep_indivs) == 0:
            raise ValueError("No individuals retained for analysis")
        return keep_indivs.tolist()

    def _read_bim_table(self, chrom: str | None) -> pd.DataFrame:
        """Read one or many `.bim` tables into a normalized DataFrame."""
        prefixes = (
            []
            if self.spec.plink_prefix is None
            else resolve_plink_prefix_group(
                (self.spec.plink_prefix,),
                chrom=chrom,
                allow_chromosome_suite=(chrom is None),
            )
        )
        if not prefixes:
            raise ValueError("PlinkRefPanel requires plink_prefix.")
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
    """
    Canonical parquet-R2 reference-panel adapter.

    The metadata sidecar is the authoritative per-SNP universe for this backend:
    ``load_metadata()`` reads A, applies ``ref_panel_snps_file`` to form A', and
    returns that restricted table to the LD-score workflow. When no sidecar is
    present, ``load_metadata()`` scans R2 parquet endpoints to synthesize
    minimal SNP metadata.
    """
    def available_chromosomes(self) -> list[str]:
        """List chromosomes from explicit config, panel directory, sidecars, or R2 files."""
        if self.spec.chromosomes is not None:
            return sorted(set(self.spec.chromosomes), key=_chrom_sort_key)

        if self.spec.r2_dir is not None:
            paths = _r2_dir_r2_paths(
                self.spec.r2_dir,
                genome_build=self.global_config.genome_build,
                chrom=None,
            )
            chromosomes = {
                chrom
                for path in paths
                if (chrom := _chromosome_from_ref_panel_r2_path(Path(path))) is not None
            }
            if chromosomes:
                return sorted(chromosomes, key=_chrom_sort_key)
            raise FileNotFoundError(f"No chr*_r2.parquet files found in R2 directory '{self.spec.r2_dir}'.")

        raise ImportError("Chromosome discovery for parquet R2 requires explicit chromosomes or r2_dir.")

    def load_metadata(self, chrom: str) -> pd.DataFrame:
        """Load, restrict, validate, and cache metadata for one chromosome."""
        chrom = normalize_chromosome(chrom)
        if chrom in self._metadata_cache:
            return self._metadata_cache[chrom].copy()
        if self.spec.keep_indivs_file is not None:
            raise ValueError("--keep-indivs-file is only supported in PLINK mode.")

        r2_paths = self.resolve_r2_paths(chrom)
        metadata_paths = self.resolve_metadata_paths(chrom)
        if metadata_paths:
            frames = [_read_metadata_table(path, chrom=chrom, global_config=self.global_config) for path in metadata_paths]
            frames = [frame for frame in frames if len(frame) > 0]
            if not frames:
                raise ValueError(f"No parquet metadata rows found for chromosome {chrom}.")
            metadata = pd.concat(frames, axis=0, ignore_index=True)
        else:
            if self.global_config.fail_on_missing_metadata:
                raise ValueError(f"Parquet reference panel metadata sidecar is missing for chromosome {chrom}.")
            if not r2_paths:
                raise ImportError("ParquetR2RefPanel.load_metadata requires metadata sidecar files or parquet R2 files.")
            LOGGER.warning(
                f"No metadata sidecar found for chromosome {chrom}; synthesizing minimal metadata from parquet R2 endpoints. "
                "SNPs without emitted R2 pairs, CM, and MAF are unavailable in this fallback mode."
            )
            metadata = _synthesize_metadata_from_r2_paths(
                r2_paths,
                chrom=chrom,
                global_config=self.global_config,
            )
            if len(metadata) == 0:
                raise ValueError(f"No parquet R2 endpoint metadata rows found for chromosome {chrom}.")
        metadata = self._apply_snp_restriction(metadata)
        metadata = self._apply_maf_filter(metadata, chrom)
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
        """Build a row-group-pruning block reader over one chromosome's R2 parquet."""
        metadata = metadata if metadata is not None else self.load_metadata(chrom)
        return kernel_ldscore.SortedR2BlockReader(
            paths=self.resolve_r2_paths(chrom),
            chrom=chrom,
            metadata=metadata,
            identifier_mode=self.global_config.snp_identifier,
            r2_bias_mode=self.spec.r2_bias_mode if r2_bias_mode is None else r2_bias_mode,
            r2_sample_size=r2_sample_size if r2_sample_size is not None else self.spec.sample_size,
            genome_build=self.global_config.genome_build,
        )

    def resolve_r2_paths(self, chrom: str, *, required: bool = True) -> list[str]:
        """Resolve parquet R2 paths for one chromosome from ``r2_dir``."""
        chrom = normalize_chromosome(chrom)
        if self.spec.r2_dir is not None:
            return _r2_dir_r2_paths(
                self.spec.r2_dir,
                genome_build=self.global_config.genome_build,
                chrom=chrom,
            )
        raise FileNotFoundError("Parquet R2 reference-panel loading requires r2_dir.")

    def resolve_metadata_paths(self, chrom: str) -> list[str]:
        """Resolve optional metadata sidecars for one chromosome."""
        chrom = normalize_chromosome(chrom)
        if self.spec.r2_dir is not None:
            return _r2_dir_metadata_paths(
                self.spec.r2_dir,
                genome_build=self.global_config.genome_build,
                chrom=chrom,
            )
        return []


def _chromosome_from_ref_panel_r2_path(path: Path) -> str | None:
    """Extract the chromosome label from ``chr{chrom}_r2.parquet``."""
    match = _REF_PANEL_R2_RE.match(path.name)
    if match is None:
        return None
    return normalize_chromosome(match.group("chrom"), context=str(path))


def _resolve_r2_build_dir(r2_dir: str | Path, genome_build: str | None) -> Path:
    """Resolve a build-specific R2 directory from user input."""
    root = Path(r2_dir)
    if not root.exists():
        raise FileNotFoundError(f"Reference panel directory does not exist: {root}")
    if not root.is_dir():
        raise NotADirectoryError(f"Reference panel path is not a directory: {root}")
    if list(root.glob("chr*_r2.parquet")):
        return root

    requested_build = genome_build if genome_build in {"hg19", "hg38"} else None
    build_dirs = {build: root / build for build in ("hg19", "hg38") if (root / build).is_dir()}
    if requested_build is not None and requested_build in build_dirs:
        return build_dirs[requested_build]
    if len(build_dirs) > 1:
        builds = ", ".join(sorted(build_dirs))
        raise ValueError(
            f"Reference panel directory '{root}' is ambiguous because it contains multiple genome-build directories ({builds}). "
            "Pass the build-specific directory or set genome_build to hg19 or hg38."
        )
    if len(build_dirs) == 1:
        return next(iter(build_dirs.values()))
    if requested_build is not None:
        return root
    return root


def _r2_dir_r2_paths(r2_dir: str | Path, *, genome_build: str | None, chrom: str | None) -> list[str]:
    """Resolve required ``chr{chrom}_r2.parquet`` files inside a panel directory."""
    build_dir = _resolve_r2_build_dir(r2_dir, genome_build)
    if chrom is not None:
        path = build_dir / f"chr{normalize_chromosome(chrom)}_r2.parquet"
        if not path.exists():
            raise FileNotFoundError(f"Required parquet R2 file is missing: {path}")
        return [str(path)]
    paths = sorted(
        build_dir.glob("chr*_r2.parquet"),
        key=lambda path: _chrom_sort_key(_chromosome_from_ref_panel_r2_path(path) or path.name),
    )
    return [str(path) for path in paths]


def _r2_dir_metadata_paths(r2_dir: str | Path, *, genome_build: str | None, chrom: str) -> list[str]:
    """Resolve the optional ``chr{chrom}_meta.tsv.gz`` sidecar from a panel directory."""
    build_dir = _resolve_r2_build_dir(r2_dir, genome_build)
    path = build_dir / f"chr{normalize_chromosome(chrom)}_meta.tsv.gz"
    return [str(path)] if path.exists() else []


def _synthesize_metadata_from_r2_paths(
    paths: list[str],
    *,
    chrom: str | None,
    global_config: GlobalConfig,
) -> pd.DataFrame:
    """Build minimal SNP metadata by scanning parquet R2 endpoint columns."""
    try:
        import pyarrow.parquet as pq
    except ImportError as exc:
        raise ImportError("pyarrow is required to synthesize metadata from parquet R2 files.") from exc

    frames = []
    for path in paths:
        pf = pq.ParquetFile(path)
        schema_names = list(pf.schema_arrow.names)
        layout = kernel_ldscore._parquet_schema_layout(schema_names)
        if layout == "canonical":
            mapping = kernel_ldscore._resolve_canonical_parquet_columns(
                schema_names,
                context=f"canonical parquet R2 schema in {path}",
            )
            table = pf.read(
                columns=[
                    mapping["CHR"],
                    mapping["POS_1"],
                    mapping["POS_2"],
                    mapping["SNP_1"],
                    mapping["SNP_2"],
                ]
            ).to_pandas()
            frames.append(
                _endpoint_frame(
                    table,
                    chr_col=mapping["CHR"],
                    left_pos_col=mapping["POS_1"],
                    right_pos_col=mapping["POS_2"],
                    left_snp_col=mapping["SNP_1"],
                    right_snp_col=mapping["SNP_2"],
                )
            )
        elif layout == "raw":
            genome_build = kernel_ldscore._require_runtime_genome_build(global_config.genome_build)
            mapping = kernel_ldscore._resolve_r2_source_columns(schema_names, context=f"raw parquet R2 schema in {path}")
            left_pos_col, right_pos_col = kernel_ldscore.get_r2_build_columns(genome_build, schema_names)
            table = pf.read(
                columns=[
                    mapping["chr"],
                    left_pos_col,
                    right_pos_col,
                    mapping["rsID_1"],
                    mapping["rsID_2"],
                ]
            ).to_pandas()
            frames.append(
                _endpoint_frame(
                    table,
                    chr_col=mapping["chr"],
                    left_pos_col=left_pos_col,
                    right_pos_col=right_pos_col,
                    left_snp_col=mapping["rsID_1"],
                    right_snp_col=mapping["rsID_2"],
                )
            )
        else:
            raise ValueError(
                "Parquet R2 metadata fallback requires canonical logical fields "
                "or the raw legacy pairwise columns."
            )
    return _finalize_endpoint_metadata(frames, chrom=chrom)


def _endpoint_frame(
    table: pd.DataFrame,
    *,
    chr_col: str,
    left_pos_col: str,
    right_pos_col: str,
    left_snp_col: str,
    right_snp_col: str,
) -> pd.DataFrame:
    """Return one SNP row per left/right endpoint in an R2 table."""
    left = pd.DataFrame(
        {
            "CHR": table[chr_col],
            "POS": table[left_pos_col],
            "SNP": table[left_snp_col],
        }
    )
    right = pd.DataFrame(
        {
            "CHR": table[chr_col],
            "POS": table[right_pos_col],
            "SNP": table[right_snp_col],
        }
    )
    return pd.concat([left, right], axis=0, ignore_index=True)


def _finalize_endpoint_metadata(frames: list[pd.DataFrame], *, chrom: str | None) -> pd.DataFrame:
    """Normalize, deduplicate, and sort synthesized endpoint metadata."""
    if not frames:
        return pd.DataFrame(columns=["CHR", "POS", "SNP", "CM"])
    metadata = pd.concat(frames, axis=0, ignore_index=True)
    if len(metadata) == 0:
        metadata["CM"] = pd.Series(pd.array([], dtype="Float64"))
        return metadata.loc[:, ["CHR", "POS", "SNP", "CM"]]
    metadata["CHR"] = metadata["CHR"].map(lambda value: normalize_chromosome(value, context="parquet R2 endpoint metadata"))
    metadata["POS"] = pd.to_numeric(metadata["POS"], errors="raise").astype(int)
    metadata["SNP"] = metadata["SNP"].astype(str)
    if chrom is not None:
        metadata = metadata.loc[metadata["CHR"] == normalize_chromosome(chrom)].reset_index(drop=True)
    metadata = metadata.drop_duplicates(subset=["CHR", "POS", "SNP"], keep="first").copy()
    metadata["_chrom_sort"] = metadata["CHR"].map(_chrom_sort_key)
    metadata = metadata.sort_values(["_chrom_sort", "POS", "SNP"], kind="mergesort").drop(columns="_chrom_sort")
    metadata["CM"] = pd.Series(pd.array([pd.NA] * len(metadata), dtype="Float64"), index=metadata.index)
    return metadata.loc[:, ["CHR", "POS", "SNP", "CM"]].reset_index(drop=True)


class RefPanelLoader:
    """Instantiate the backend requested by a :class:`RefPanelConfig`."""
    def __init__(self, global_config: GlobalConfig, ref_panel_config: RefPanelConfig | None = None) -> None:
        """Store shared configuration for future backend instantiation."""
        self.global_config = global_config
        self.ref_panel_config = ref_panel_config or RefPanelConfig()

    def load(self, ref_panel_spec: RefPanelConfig) -> RefPanel:
        """Return a concrete reference-panel adapter for ``ref_panel_spec``."""
        backend = ref_panel_spec.backend
        if backend == "plink":
            return PlinkRefPanel(self.global_config, ref_panel_spec)
        if backend == "parquet_r2":
            return ParquetR2RefPanel(self.global_config, ref_panel_spec)
        raise ValueError(f"Unsupported reference-panel backend: {backend}")


def _chrom_sort_key(chrom: str) -> tuple[int, str]:
    """Return a stable chromosome sort key matching the package-wide ordering."""
    return chrom_sort_key(chrom)


def _read_metadata_table(path: str | Path, chrom: str | None, global_config: GlobalConfig) -> pd.DataFrame:
    """Read one reference-panel metadata table into the normalized column set."""
    df = pd.read_csv(path, sep=r"\s+", compression="gzip" if str(path).endswith(".gz") else None)
    snp_identifier = normalize_snp_identifier_mode(global_config.snp_identifier)
    assert global_config.genome_build in {"hg19", "hg38", None}, (
        f"genome_build reached kernel as {global_config.genome_build!r}; "
        "should have been resolved at workflow entry."
    )
    context = str(path)

    chr_col = None
    pos_col = None
    snp_col = None
    try:
        chr_col = resolve_required_column(df.columns, REFERENCE_METADATA_SPEC_MAP["CHR"], context=context)
        pos_col = resolve_required_column(df.columns, REFERENCE_METADATA_SPEC_MAP["POS"], context=context)
    except ValueError:
        # Metadata may be rsid-only; validate the required identifier columns
        # after probing both coordinate and SNP schemas.
        pass
    try:
        snp_col = resolve_required_column(df.columns, REFERENCE_METADATA_SPEC_MAP["SNP"], context=context)
    except ValueError:
        # Metadata may be chr_pos-only; mode-specific validation below decides
        # whether the missing SNP column is an error.
        pass

    if snp_identifier == "rsid" and snp_col is None:
        raise ValueError(f"{path} must contain a SNP column in rsid mode.")
    if snp_identifier == "chr_pos" and (chr_col is None or pos_col is None):
        raise ValueError(f"{path} must contain CHR and POS columns in chr_pos mode.")

    out = pd.DataFrame(index=df.index)
    if chr_col is not None:
        out["CHR"] = df[chr_col].map(lambda value: normalize_chromosome(value, context=context))
    if pos_col is not None:
        out["POS"] = pd.to_numeric(df[pos_col], errors="raise").astype(int)
    if snp_col is not None:
        out["SNP"] = df[snp_col].astype(str)

    cm_col = resolve_optional_column(df.columns, REFERENCE_METADATA_SPEC_MAP["CM"], context=context)
    maf_col = resolve_optional_column(df.columns, REFERENCE_METADATA_SPEC_MAP["MAF"], context=context)
    if cm_col is not None:
        out["CM"] = pd.to_numeric(df[cm_col], errors="coerce")
    else:
        out["CM"] = pd.NA
    if maf_col is not None:
        maf = pd.to_numeric(df[maf_col], errors="coerce").astype(float)
        out["MAF"] = pd.Series(maf).map(lambda value: value if pd.isna(value) else min(value, 1.0 - value))
    elif global_config.fail_on_missing_metadata:
        raise ValueError(f"{path} is missing MAF metadata.")

    if chrom is not None and "CHR" in out.columns:
        out = out.loc[out["CHR"] == normalize_chromosome(chrom, context=context)].reset_index(drop=True)
    return out
