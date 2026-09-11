"""Reference-panel configuration and workflow-to-kernel backend adapters.

This module prepares chromosome-scoped reference data for LD-score projection.
The PLINK and parquet adapters own identity and SNP restrictions, sample/MAF
selection, annotation alignment, reader policy, LD windows, and file lifetime.
Numerical kernels consume the resulting aligned state without resolving paths
or reopening a panel. Metadata inspection returns caller-owned chromosome tables without caching them.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass
from types import SimpleNamespace
import gzip
import logging
from pathlib import Path
import re
from typing import Any

import numpy as np
import pandas as pd

from .._coordinates import CHR_POS_KEY_COLUMN, build_chr_pos_key_frame
from ..chromosome_inference import chrom_sort_key, normalize_chromosome
from ..column_inference import (
    A1_COLUMN_SPEC,
    A2_COLUMN_SPEC,
    REFERENCE_METADATA_SPEC_MAP,
    normalize_genome_build,
    resolve_optional_column,
    resolve_required_column,
)
from ..config import GlobalConfig, LDScoreConfig, RefPanelConfig, validate_config_compatibility
from ..errors import EmptyReferenceSNPs, LDSCConfigError, LDSCDependencyError, LDSCInputError, LDSCUsageError
from ..genome_build_inference import resolve_genome_build, validate_auto_genome_build_mode
from ..path_resolution import resolve_plink_prefix, resolve_plink_prefix_group, split_cli_path_tokens
from . import formats as legacy_parse
from . import ldscore as kernel_ldscore
from .ldscore_projection import MappedAnnotations
from .identifiers import (
    build_snp_id_series,
    normalize_snp_identifier_mode,
    read_snp_restriction_keys,
    validate_unique_snp_ids,
)
from .snp_identity import (
    clean_identity_artifact_table,
    empty_identity_drop_frame,
    identity_base_mode,
    identity_mode_family,
    is_allele_aware_mode,
    restriction_membership_mask,
    validate_identity_artifact_metadata,
)

LOGGER = logging.getLogger("LDSC.ref_panel")
_REF_PANEL_R2_RE = re.compile(r"^chr(?P<chrom>.+)_r2\.parquet$", flags=re.IGNORECASE)
_REF_PANEL_ARTIFACT_DOC = "docs/troubleshooting.md#build-r2-panel-reference-panel-artifact-is-incompatible"
_REF_PANEL_EMPTY_DOC = "docs/troubleshooting.md#build-r2-panel-no-reference-panel-artifacts-were-produced"


def _snp_id_series_for_matching(metadata: pd.DataFrame, snp_identifier: str, *, context: str) -> pd.Series:
    """Build SNP IDs for a match boundary, dropping missing CHR/POS in base coordinate mode."""
    mode = normalize_snp_identifier_mode(snp_identifier)
    if identity_mode_family(mode) == "rsid" or is_allele_aware_mode(mode):
        return build_snp_id_series(metadata, mode)
    keyed, _report = build_chr_pos_key_frame(
        metadata,
        context=context,
        drop_missing=True,
        logger=LOGGER,
    )
    keys = pd.Series(pd.NA, index=metadata.index, dtype="object")
    keys.loc[keyed.index] = keyed[CHR_POS_KEY_COLUMN].astype(str)
    return keys


@dataclass(frozen=True)
class _R2SchemaMeta:
    """Parsed LDSC R2 schema metadata from one parquet file."""

    n_samples: int | None
    r2_bias: str | None


def _read_r2_schema_meta(path: str) -> _R2SchemaMeta:
    """
    Read LDSC R2 sample-size and bias metadata from an Arrow parquet schema.

    Missing keys are returned as ``None`` so legacy parquet files keep their
    historical behavior. A file with ``ldsc:n_samples`` but no ``ldsc:r2_bias``
    is treated as malformed raw-R2 metadata and resolved to ``r2_bias="raw"``
    with a warning.
    """
    try:
        import pyarrow.parquet as pq
    except ImportError:
        return _R2SchemaMeta(n_samples=None, r2_bias=None)

    raw_meta = pq.read_schema(path).metadata or {}
    n_raw = raw_meta.get(b"ldsc:n_samples")
    bias_raw = raw_meta.get(b"ldsc:r2_bias")
    n_samples = int(n_raw.decode("utf-8")) if n_raw is not None else None
    r2_bias = bias_raw.decode("utf-8") if bias_raw is not None else None

    if n_samples is not None and r2_bias is None:
        LOGGER.warning(f"'{path}' has ldsc:n_samples but no ldsc:r2_bias; treating as 'raw'.")
        r2_bias = "raw"

    return _R2SchemaMeta(n_samples=n_samples, r2_bias=r2_bias)


def _read_identity_schema_meta(path: str, *, expected_artifact_type: str) -> dict[str, object]:
    """Read and validate LDSC identity metadata from package-written parquet schema metadata."""
    try:
        import pyarrow.parquet as pq
    except ImportError as exc:
        raise LDSCDependencyError(
            "Reading LDSC R2 identity schema metadata requires pyarrow. Most likely "
            "parquet reference-panel input was requested in an environment without pyarrow. "
            "Install pyarrow or activate the LDSC environment that includes it."
        ) from exc

    raw = pq.read_schema(path).metadata or {}
    required_keys = {
        b"ldsc:artifact_type",
        b"ldsc:snp_identifier",
        b"ldsc:genome_build",
    }
    if not required_keys.issubset(raw):
        raise LDSCInputError(
            "Could not read LDSC reference-panel R2 artifact metadata: required "
            "identity/provenance keys are missing from the parquet schema. Most likely "
            "the R2 file was written by an older LDSC version or by another tool. "
            "Regenerate the reference panel with the current `ldsc build-r2-panel`. "
            f"Other causes & fixes: {_REF_PANEL_ARTIFACT_DOC}"
        )
    metadata = {
        "artifact_type": raw[b"ldsc:artifact_type"].decode("utf-8"),
        "snp_identifier": raw[b"ldsc:snp_identifier"].decode("utf-8"),
        "genome_build": raw[b"ldsc:genome_build"].decode("utf-8"),
    }
    validate_identity_artifact_metadata(metadata, expected_artifact_type=expected_artifact_type)
    return metadata


def _r2_path_has_ldsc_package_schema(path: str) -> bool:
    """Return whether parquet schema metadata marks a package-written canonical R2 artifact."""
    try:
        import pyarrow.parquet as pq
    except ImportError:
        return False
    raw = pq.read_schema(path).metadata or {}
    package_keys = {
        b"ldsc:sorted_by_build",
        b"ldsc:row_group_size",
        b"ldsc:artifact_type",
    }
    return any(key in raw for key in package_keys)


def _validate_r2_identity_matches_config(metadata: dict[str, object], global_config: GlobalConfig) -> None:
    """Reject package R2 artifacts built under different identity assumptions."""
    artifact_config = GlobalConfig(
        snp_identifier=str(metadata["snp_identifier"]),
        genome_build=metadata.get("genome_build"),
    )
    validate_config_compatibility(
        global_config,
        artifact_config,
        context="runtime config and parquet R2 artifact",
    )


def _resolve_r2_bias_from_meta(
    r2_bias_mode: str | None,
    r2_sample_size: float | None,
    meta: _R2SchemaMeta,
) -> tuple[str, float | None]:
    """
    Resolve effective R2 bias mode and correction sample size.

    User-supplied values take precedence over parquet metadata. When no user
    mode is supplied, stored ``ldsc:r2_bias`` selects ``"unbiased"`` or
    ``"raw"``; raw mode auto-fills the sample size from ``ldsc:n_samples`` when
    available. Legacy files with neither key default to unbiased, matching the
    pre-metadata behavior.
    """
    stored_bias = meta.r2_bias
    stored_n = float(meta.n_samples) if meta.n_samples is not None else None

    if stored_bias is None and stored_n is not None:
        LOGGER.warning("R2 schema metadata has ldsc:n_samples but no ldsc:r2_bias; treating as 'raw'.")
        stored_bias = "raw"

    effective_bias = r2_bias_mode if r2_bias_mode is not None else (stored_bias or "unbiased")

    if effective_bias == "unbiased":
        if r2_sample_size is not None:
            LOGGER.warning(
                "r2_sample_size is ignored because R2 values are already unbiased "
                "(ldsc:r2_bias=unbiased in parquet schema metadata)."
            )
        return "unbiased", None

    resolved_n = r2_sample_size if r2_sample_size is not None else stored_n
    return "raw", resolved_n


class RefPanel(ABC):
    """Abstract chromosome-scoped reference-panel interface."""
    def __init__(self, global_config: GlobalConfig, spec: RefPanelConfig) -> None:
        """Store shared configuration without chromosome payload caches."""
        validate_auto_genome_build_mode(global_config.snp_identifier, global_config.genome_build)
        self.global_config = global_config
        self.spec = spec

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

    @abstractmethod
    def prepare_chromosome(
        self, chrom: str, annotations: kernel_ldscore.AnnotationBundle,
        config: LDScoreConfig, *, genetic_map: pd.DataFrame | None = None,
    ) -> kernel_ldscore.PreparedChromosome:
        """Prepare aligned reference data and an owned reader for LD-score computation.

        Parameters
        ----------
        chrom : str
            Chromosome label to resolve in this panel.
        annotations : ldsc._kernel.ldscore.AnnotationBundle
            SNP metadata in genomic order and numeric annotation columns in
            baseline-then-query order. Input tables are not modified.
        config : LDScoreConfig
            Window geometry and whole-chromosome-window policy.
        genetic_map : pandas.DataFrame or None, optional
            Already resolved PLINK genetic map for cM windows. If absent, the
            PLINK adapter resolves configured map sources when needed.

        Returns
        -------
        ldsc._kernel.ldscore.PreparedChromosome
            Reference metadata, a selected-read annotation mapping, window bounds, and
            an owned reader over the same retained SNP rows. Panel identity,
            SNP/sample/MAF restrictions, authoritative reference CM/MAF, and
            reader bias settings are resolved together. Use the result as a
            context manager to close the reader on success or failure.

        Notes
        -----
        Preparation reads inputs but writes no artifacts. Metadata inspection
        through ``load_metadata`` reads independently and retains no cache.
        """
        raise NotImplementedError

    def filter_to_snps(self, chrom: str, snps: set[str] | list[str]) -> pd.DataFrame:
        """Subset chromosome metadata to the requested SNP identifiers."""
        metadata = self.load_metadata(chrom)
        keys = _snp_id_series_for_matching(
            metadata,
            self.global_config.snp_identifier,
            context=f"reference-panel SNP filtering for chromosome {chrom}",
        )
        keep = keys.isin(set(snps))
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
        keep = metadata["MAF"] >= maf_min
        removed = int((~keep).sum())
        if removed:
            LOGGER.info(f"Removed {removed} reference-panel SNPs with MAF < {maf_min} on chromosome {chrom}.")
        return metadata.loc[keep].reset_index(drop=True)

    def _apply_snp_restriction(self, metadata: pd.DataFrame) -> pd.DataFrame:
        """Apply the explicit custom LD-reference SNP restriction when configured."""
        return self._apply_keep_restriction(metadata)

    def _apply_keep_restriction(self, metadata: pd.DataFrame) -> pd.DataFrame:
        """Filter metadata rows to ``RefPanelConfig.ref_panel_snps_file`` when set."""
        restrict_path = self.spec.ref_panel_snps_file
        if restrict_path is None or len(metadata) == 0:
            return metadata
        restriction = read_snp_restriction_keys(
            restrict_path,
            self.global_config.snp_identifier,
            genome_build=self.global_config.genome_build,
            logger=LOGGER,
        )
        keep = restriction_membership_mask(
            metadata,
            restriction,
            self.global_config.snp_identifier,
            context=f"reference-panel restriction matching for {restrict_path}",
        )
        return metadata.loc[keep].reset_index(drop=True)

    def _validate_metadata(self, metadata: pd.DataFrame, chrom: str) -> pd.DataFrame:
        """Reset row order and validate identifier uniqueness for one chromosome."""
        metadata = metadata.reset_index(drop=True)
        validate_unique_snp_ids(metadata, self.global_config.snp_identifier, context=f"{type(self).__name__}[{chrom}]")
        return metadata


class PlinkRefPanel(RefPanel):
    """PLINK-backed reference-panel adapter.

    ``chromosome_prefixes`` optionally binds normalized chromosome labels to
    complete PLINK prefixes validated by the calling workflow. This bypasses
    filename-based routing for content-authoritative query input globs.
    """

    def __init__(self, global_config: GlobalConfig, spec: RefPanelConfig, *,
                 chromosome_prefixes: dict[str, str] | None = None) -> None:
        """Store panel configuration and optional validated chromosome routing."""
        super().__init__(global_config, spec)
        self._chromosome_prefixes = None if chromosome_prefixes is None else dict(chromosome_prefixes)

    def available_chromosomes(self) -> list[str]:
        """List normalized chromosomes present in the resolved PLINK inputs."""
        chromosomes = set()
        for prefix in self._bim_prefixes(None):
            with pd.read_csv(prefix + ".bim", sep=r"\s+", header=None, usecols=[0], chunksize=65536) as reader:
                for chunk in reader:
                    chromosomes.update(chunk.iloc[:,0].map(normalize_chromosome))
        return sorted(chromosomes, key=_chrom_sort_key)

    def _load_source(self, chrom: str):
        """Resolve physical BED rows and apply panel identity/SNP restrictions once."""
        prefix = (self._chromosome_prefixes[normalize_chromosome(chrom)] if self._chromosome_prefixes is not None
                  else resolve_plink_prefix(self.spec.plink_prefix, chrom=chrom))
        bim = legacy_parse.PlinkBIMFile(prefix + ".bim")
        fam = legacy_parse.PlinkFAMFile(prefix + ".fam")
        metadata = bim.df.rename(columns={"BP": "POS"}).copy()
        metadata["_raw_index"] = np.arange(len(metadata))
        metadata["CHR"] = metadata["CHR"].map(normalize_chromosome)
        metadata["SNP"] = metadata["SNP"].astype(str)
        metadata["POS"] = pd.to_numeric(metadata["POS"], errors="raise").astype(np.int64)
        metadata["CM"] = pd.to_numeric(metadata["CM"], errors="coerce")
        metadata[["A1", "A2"]] = metadata[["A1", "A2"]].astype(str)
        metadata = metadata.loc[metadata["CHR"] == normalize_chromosome(chrom)].reset_index(drop=True)
        if metadata.empty:
            raise LDSCInputError(f"Reference-panel loading found no PLINK metadata rows for chromosome {chrom}.")
        drops = empty_identity_drop_frame()
        if self.global_config.snp_identifier in {"rsid", "chr_pos"}:
            cleanup = clean_identity_artifact_table(
                metadata, self.global_config.snp_identifier,
                context=f"PLINK reference-panel metadata chromosome {chrom}",
                stage="plink_reference_identity_cleanup", logger=LOGGER,
            )
            metadata, drops = cleanup.cleaned, cleanup.dropped
            if metadata.empty:
                raise LDSCInputError(
                    f"Reference-panel loading retained no PLINK rows on chromosome {chrom} after duplicate SNP identity cleanup."
                )
        metadata = self._apply_snp_restriction(metadata)
        metadata = self._validate_metadata(metadata, chrom)
        return prefix, bim, fam, metadata, drops

    def _load_genotypes(self, chrom, source, annotations=None, *, keep_indivs=None, maf_min=None):
        """Open one BED and align its retained rows to annotations without another read."""
        prefix, bim, fam, reference, drops = source
        if annotations is None:
            annotations = kernel_ldscore.AnnotationBundle(
                reference.drop(columns="_raw_index"), pd.DataFrame(index=reference.index), [], [],
            )
        aligned = _align_annotations(annotations, reference, chrom, self.global_config, "PLINK")
        mode = self.global_config.snp_identifier
        raw_keys = kernel_ldscore.identifier_keys(reference, mode)
        keys = kernel_ldscore.identifier_keys(aligned.metadata, mode)
        physical_rows = pd.Series(reference["_raw_index"].to_numpy(), index=raw_keys).loc[keys].tolist()
        if keep_indivs is None:
            keep_indivs = kernel_ldscore.resolve_keep_individuals(self.spec.keep_indivs_file, fam)
        geno = kernel_ldscore.PlinkBEDFile(
            prefix + ".bed", len(fam.IDList), bim, keep_snps=physical_rows,
            keep_indivs=keep_indivs, mafMin=self.spec.maf_min if maf_min is None else maf_min,
        )
        try:
            metadata = pd.DataFrame(geno.df, columns=geno.colnames).rename(columns={"BP": "POS"})
            metadata["CHR"] = metadata["CHR"].map(normalize_chromosome)
            metadata["SNP"] = metadata["SNP"].astype(str)
            metadata["POS"] = pd.to_numeric(metadata["POS"], errors="raise").astype(np.int64)
            metadata["CM"] = pd.to_numeric(metadata["CM"], errors="coerce")
            metadata["MAF"] = pd.to_numeric(metadata["MAF"], errors="coerce")
            metadata = metadata.merge(reference[["CHR", "SNP", "POS", "A1", "A2"]],
                                      on=["CHR", "SNP", "POS"], how="left", sort=False)
            retained_keys = kernel_ldscore.identifier_keys(metadata, mode)
            indices = pd.Series(np.arange(len(keys)), index=keys).loc[retained_keys].to_numpy(dtype=np.int64)
            values = MappedAnnotations(aligned.annotations, indices, aligned.annotations.columns)
            return geno, metadata, values, drops, len(aligned.metadata)
        except BaseException:
            geno.close()
            raise

    def load_metadata(self, chrom: str) -> pd.DataFrame:
        """Load one chromosome's metadata; the caller owns its lifetime."""
        chrom = normalize_chromosome(chrom)
        source = self._load_source(chrom)
        try:
            geno, metadata, _, _, _ = self._load_genotypes(chrom, source)
        except Exception:
            if self.spec.maf_min is not None or self.spec.keep_indivs_file is not None:
                raise
            # Preserve metadata-only BIM inspection when BED data is unavailable.
            metadata = source[3].drop(columns="_raw_index").reset_index(drop=True)
        else:
            geno.close()
        return metadata

    def build_reader(self, chrom: str, keep_snps=None, keep_indivs=None, maf_min=None):
        """Build an owned BED reader through the same preparation used for LD scores.

        Configured panel restrictions apply before optional ``keep_snps`` keys,
        sample-index overrides, and MAF overrides. The caller must close the
        returned reader.
        """
        source = self._load_source(chrom)
        if keep_snps is not None:
            prefix, bim, fam, metadata, drops = source
            keys = _snp_id_series_for_matching(metadata, self.global_config.snp_identifier,
                                              context=f"PLINK chromosome {chrom}")
            source = prefix, bim, fam, metadata.loc[keys.isin(keep_snps)].reset_index(drop=True), drops
        return self._load_genotypes(chrom, source, keep_indivs=keep_indivs, maf_min=maf_min)[0]

    def prepare_chromosome(self, chrom, annotations, config, *, genetic_map=None):
        """Prepare one genotype-filtered chromosome with aligned annotations and LD windows."""
        source = self._load_source(chrom)
        geno, metadata, values, drops, input_rows = self._load_genotypes(chrom, source, annotations)
        try:
            if config.ld_wind_cm is not None:
                if genetic_map is None:
                    genetic_map = _resolve_genetic_map(self, chrom, metadata)
                if genetic_map is not None:
                    from .ref_panel_builder import interpolate_genetic_map_cm
                    metadata["CM"] = interpolate_genetic_map_cm(
                        normalize_chromosome(chrom), metadata["POS"].to_numpy(dtype=np.int64), genetic_map,
                    )
                else:
                    kernel_ldscore.assert_cm_usable(metadata["CM"], chrom)
            block_left = _prepare_window(metadata, chrom, config, self.global_config.snp_identifier)
            return kernel_ldscore.PreparedChromosome(
                backend="plink", reader=geno, metadata=metadata,
                annotations=values, block_left=block_left,
                baseline_columns=list(annotations.baseline_columns), query_columns=list(annotations.query_columns),
                reference_rows_before_genotype_qc=input_rows,
                genotype_qc_removed=geno.genotype_qc_removed, maf_removed=geno.maf_removed,
                selected_individual_count=int(geno.n),
                cm_source="explicit_genetic_map" if genetic_map is not None else "bim_cm", identity_drops=drops,
            )
        except BaseException:
            geno.close()
            raise

    def _bim_prefixes(self, chrom: str | None) -> list[str]:
        """Resolve source prefixes without retaining BIM tables."""
        if self._chromosome_prefixes is not None:
            prefixes = ([self._chromosome_prefixes[normalize_chromosome(chrom)]] if chrom is not None
                        else list(dict.fromkeys(self._chromosome_prefixes.values())))
        else:
            prefixes = ([] if self.spec.plink_prefix is None else resolve_plink_prefix_group(
                (self.spec.plink_prefix,), chrom=chrom, allow_chromosome_suite=(chrom is None),
            ))
        if not prefixes:
            raise LDSCUsageError(
                "PLINK reference-panel loading requires a PLINK prefix. Most likely "
                "`RefPanelConfig(backend='plink')` was used without `plink_prefix` or "
                "the prefix token resolved to no files. Pass a valid PLINK prefix."
            )
        return prefixes


class ParquetR2RefPanel(RefPanel):
    """
    Canonical index-format parquet R2 reference-panel adapter.

    The metadata sidecar (``chrN_meta.tsv.gz``) is mandatory and authoritative:
    ``load_metadata()`` reads it, applies ``ref_panel_snps_file`` to form A',
    and returns the restricted table to the LD-score workflow. A missing sidecar
    is a hard error. The paired R2 parquet (``chrN_r2.parquet``) stores only
    sidecar-row indices and is meaningless without the exact matching sidecar.
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
            raise LDSCInputError(
                f"Reference-panel loading found no `chr*_r2.parquet` files in R2 "
                f"directory '{self.spec.r2_dir}'. Most likely the directory is not an "
                "`ldsc build-r2-panel` output directory or the wrong genome-build child "
                "directory was selected. Pass the directory containing canonical R2 parquet files. "
                f"Other causes & fixes: {_REF_PANEL_ARTIFACT_DOC}"
            )

        raise LDSCUsageError(
            "Parquet R2 chromosome discovery requires explicit chromosomes or `r2_dir`. "
            "Most likely `RefPanelConfig(backend='parquet_r2')` was used without an R2 "
            "directory. Pass `r2_dir` or an explicit chromosome list."
        )

    def load_metadata(self, chrom: str) -> pd.DataFrame:
        """Load, restrict, and validate metadata owned by this chromosome consumer."""
        chrom = normalize_chromosome(chrom)
        if self.spec.keep_indivs_file is not None:
            raise LDSCUsageError(
                "Reference-panel loading cannot apply `--keep-indivs-file` in parquet R2 mode. "
                "Most likely individual filtering was combined with precomputed parquet R2 "
                "input. Remove the keep file or use PLINK mode."
            )

        r2_paths = self.resolve_r2_paths(chrom)
        metadata_paths = self.resolve_metadata_paths(chrom)
        metadata_is_external = False
        if metadata_paths:
            require_identity_metadata = any(_r2_path_has_ldsc_package_schema(path) for path in r2_paths)
            metadata_is_external = not all(_read_metadata_sidecar_identity(path) is not None for path in metadata_paths)
            frames = [
                _read_metadata_table(
                    path,
                    chrom=chrom,
                    global_config=self.global_config,
                    require_identity_metadata=require_identity_metadata,
                )
                for path in metadata_paths
            ]
            frames = [frame for frame in frames if len(frame) > 0]
            if not frames:
                raise LDSCInputError(
                    f"Reference-panel loading found no metadata rows for chromosome {chrom}. "
                    "Most likely the metadata sidecar exists but contains no rows for this "
                    "chromosome after filtering. Check the sidecar and selected chromosome."
                )
            metadata = pd.concat(frames, axis=0, ignore_index=True)
        else:
            raise LDSCInputError(
                f"Reference-panel metadata sidecar is missing for chromosome {chrom}. "
                "Index-format R2 parquets require their `chrN_meta.tsv.gz` sidecar because "
                "the parquet stores only sidecar-row indices. Most likely the sidecar was "
                "not copied with the R2 parquet. Restore the matching sidecar or regenerate "
                "the panel. "
                f"Other causes & fixes: {_REF_PANEL_ARTIFACT_DOC}"
            )
        metadata = self._apply_snp_restriction(metadata)
        metadata = self._apply_maf_filter(metadata, chrom)
        if metadata_is_external and not metadata.empty:
            cleanup = clean_identity_artifact_table(
                metadata,
                self.global_config.snp_identifier,
                context=f"external parquet reference-panel metadata chromosome {chrom}",
                stage="parquet_metadata_identity_cleanup",
                logger=LOGGER,
            )
            metadata = cleanup.cleaned
            if len(metadata) == 0:
                raise LDSCInputError(
                    f"Reference-panel loading retained no parquet metadata rows on chromosome "
                    f"{chrom} after SNP identity cleanup. Most likely metadata rows have "
                    "missing or duplicate SNP identities under the active identifier mode. "
                    "Regenerate the reference panel or use an identifier mode compatible "
                    f"with the sidecar. Other causes & fixes: {_REF_PANEL_EMPTY_DOC}"
                )
        metadata = self._validate_metadata(metadata, chrom)
        return metadata

    def build_reader(
        self,
        chrom: str,
        metadata: pd.DataFrame | None = None,
        r2_bias_mode: str | None = None,
        r2_sample_size: float | None = None,
    ):
        """
        Build an owned streaming reader with R2 bias settings resolved.

        Runtime overrides are merged with ``RefPanelConfig`` first, then the
        first chromosome parquet's LDSC schema metadata is used to auto-fill
        missing R2 bias mode and sample size before constructing
        ``SortedR2BlockReader``. ``prepare_chromosome`` uses this same factory;
        direct callers must close the returned reader.
        """
        metadata = metadata if metadata is not None else self.load_metadata(chrom)
        paths = self.resolve_r2_paths(chrom)
        # Bias mode comes from parquet metadata (resolved below); sample size
        # falls back to any explicit RefPanelConfig.sample_size.
        effective_bias = r2_bias_mode
        effective_n = r2_sample_size if r2_sample_size is not None else self.spec.sample_size

        if paths:
            if _r2_path_has_ldsc_package_schema(paths[0]):
                identity_metadata = _read_identity_schema_meta(paths[0], expected_artifact_type="ref_panel_r2")
                _validate_r2_identity_matches_config(identity_metadata, self.global_config)
            stored = _read_r2_schema_meta(paths[0])
            effective_bias, effective_n = _resolve_r2_bias_from_meta(
                effective_bias,
                effective_n,
                stored,
            )

        return kernel_ldscore.SortedR2BlockReader(
            paths=paths,
            chrom=chrom,
            metadata=metadata,
            identifier_mode=self.global_config.snp_identifier,
            r2_bias_mode=effective_bias,
            r2_sample_size=effective_n,
            genome_build=self.global_config.genome_build,
        )

    def prepare_chromosome(self, chrom, annotations, config, *, genetic_map=None):
        """Prepare a sidecar-aligned streaming reader with one resolved R2 policy."""
        reference = self.load_metadata(chrom)
        aligned = _align_annotations(annotations, reference, chrom, self.global_config, "parquet")
        keys = kernel_ldscore.identifier_keys(aligned.metadata, self.global_config.snp_identifier)
        reference_keys = kernel_ldscore.identifier_keys(reference, self.global_config.snp_identifier)
        lookup = reference.set_index(reference_keys)
        metadata = aligned.metadata.copy()
        for column in ("CM", "MAF"):
            if column in reference:
                metadata[column] = lookup.loc[keys, column].to_numpy()
        if config.ld_wind_cm is not None:
            _resolve_genetic_map(self, chrom, reference)
        block_left = _prepare_window(metadata, chrom, config, self.global_config.snp_identifier)
        kernel_ldscore.validate_ldscore_window_within_r2_panel_window(
            config, parquet_paths=self.resolve_r2_paths(chrom), chrom=chrom,
        )
        values = aligned.annotations
        reader = self.build_reader(chrom, metadata=metadata)
        return kernel_ldscore.PreparedChromosome(
            backend="parquet_r2", reader=reader, metadata=metadata,
            annotations=values, block_left=block_left,
            baseline_columns=list(annotations.baseline_columns), query_columns=list(annotations.query_columns),
            cm_source="parquet_sidecar",
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
        raise LDSCUsageError(
            "Parquet R2 reference-panel loading requires `r2_dir`. Most likely "
            "`RefPanelConfig(backend='parquet_r2')` was used without an R2 directory. "
            "Pass the directory containing `chrN_r2.parquet` files."
        )

    def resolve_metadata_paths(self, chrom: str) -> list[str]:
        """Resolve the required metadata sidecar for one chromosome."""
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
        raise LDSCInputError(
            f"Reference-panel R2 directory does not exist: '{root}'. Most likely the "
            "path is misspelled or relative to a different working directory. Pass the "
            "existing `ldsc build-r2-panel` output directory."
        )
    if not root.is_dir():
        raise LDSCInputError(
            f"Reference-panel R2 path is not a directory: '{root}'. Most likely a file "
            "path was passed where an R2 directory was expected. Pass the directory "
            "containing `chrN_r2.parquet` files."
        )
    if list(root.glob("chr*_r2.parquet")):
        return root

    requested_build = genome_build if genome_build in {"hg19", "hg38"} else None
    build_dirs = {build: root / build for build in ("hg19", "hg38") if (root / build).is_dir()}
    if requested_build is not None and requested_build in build_dirs:
        return build_dirs[requested_build]
    if len(build_dirs) > 1:
        builds = ", ".join(sorted(build_dirs))
        raise LDSCInputError(
            f"Reference-panel directory '{root}' is ambiguous because it contains multiple "
            f"genome-build directories ({builds}). Most likely a parent output directory "
            "was passed while genome_build is unset or auto. Pass the build-specific "
            "directory or set genome_build to hg19 or hg38."
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
            raise LDSCInputError(
                f"Required parquet R2 file is missing: '{path}'. Most likely the "
                "reference panel was not built for this chromosome/build, or the wrong "
                "R2 directory was selected. Regenerate the missing chromosome or pass "
                "the correct `r2_dir`. "
                f"Other causes & fixes: {_REF_PANEL_ARTIFACT_DOC}"
            )
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


class RefPanelLoader:
    """Instantiate the backend requested by a :class:`RefPanelConfig`."""
    def __init__(self, global_config: GlobalConfig, ref_panel_config: RefPanelConfig | None = None) -> None:
        """Store shared configuration for future backend instantiation."""
        self.global_config = global_config
        self.ref_panel_config = ref_panel_config or RefPanelConfig()

    def load(self, ref_panel_spec: RefPanelConfig) -> RefPanel:
        """Return a concrete reference-panel adapter for ``ref_panel_spec``."""
        global_config = self._global_config_for_spec(ref_panel_spec)
        backend = ref_panel_spec.backend
        if backend == "plink":
            return PlinkRefPanel(global_config, ref_panel_spec)
        if backend == "parquet_r2":
            return ParquetR2RefPanel(global_config, ref_panel_spec)
        raise LDSCConfigError(
            f"Unsupported reference-panel backend {backend!r}. Most likely a Python "
            "caller supplied an invalid RefPanelConfig backend. Use `plink` or `parquet_r2`."
        )

    def _global_config_for_spec(self, ref_panel_spec: RefPanelConfig) -> GlobalConfig:
        """Resolve direct Python auto-build restrictions before backend loading."""
        mode = normalize_snp_identifier_mode(self.global_config.snp_identifier)
        if identity_mode_family(mode) != "chr_pos" or self.global_config.genome_build != "auto":
            return self.global_config
        if not ref_panel_spec.ref_panel_snps_file:
            return self.global_config
        resolved_build = _infer_restricted_ref_panel_build(ref_panel_spec)
        return GlobalConfig(
            snp_identifier=mode,
            genome_build=resolved_build,
            log_level=self.global_config.log_level,
        )


def _infer_restricted_ref_panel_build(ref_panel_spec: RefPanelConfig) -> str:
    """Infer a concrete build for direct ``chr_pos`` ref-panel restrictions."""
    if ref_panel_spec.ref_panel_snps_file:
        header_build = _infer_restriction_header_build(Path(ref_panel_spec.ref_panel_snps_file))
        if header_build is not None:
            return header_build
    try:
        if ref_panel_spec.backend == "parquet_r2":
            inferred = _infer_parquet_r2_build(ref_panel_spec)
        elif ref_panel_spec.backend == "plink":
            inferred = _infer_plink_ref_panel_build(ref_panel_spec)
        else:
            inferred = None
    except Exception as exc:
        raise LDSCInputError(
            "Reference-panel loading could not infer genome_build='auto' before applying "
            "chr_pos SNP restrictions. Most likely the restriction file, PLINK BIM, or "
            "parquet R2 metadata did not contain enough build evidence. Pass "
            "GlobalConfig(genome_build='hg19') or GlobalConfig(genome_build='hg38') explicitly."
        ) from exc
    if inferred not in {"hg19", "hg38"}:
        raise LDSCInputError(
            "Reference-panel loading could not infer genome_build='auto' before applying "
            "chr_pos SNP restrictions. Most likely the restriction file and reference panel "
            "do not expose concrete build metadata. Pass GlobalConfig(genome_build='hg19') "
            "or GlobalConfig(genome_build='hg38') explicitly."
        )
    return inferred


def _infer_restriction_header_build(path: Path) -> str | None:
    """Infer build from an unambiguous build-specific restriction POS header."""
    header, _rows, _delimiter = _parse_restriction_rows(path)
    normalized = {str(column).upper().replace("-", "_") for column in header}
    has_hg19 = any(token in normalized for token in {"HG19_POS", "HG19_BP", "HG37_POS", "GRCH37_POS"})
    has_hg38 = any(token in normalized for token in {"HG38_POS", "HG38_BP", "GRCH38_POS"})
    if has_hg19 and not has_hg38:
        return "hg19"
    if has_hg38 and not has_hg19:
        return "hg38"
    return None


def _parse_restriction_rows(path: Path):
    """Reuse the identifier parser for header inspection."""
    from .identifiers import _parse_restriction_rows as parse_rows

    return parse_rows(path)


def _infer_parquet_r2_build(ref_panel_spec: RefPanelConfig) -> str | None:
    """Infer reference-panel build from parquet R2 schema metadata."""
    if ref_panel_spec.r2_dir is None:
        return None
    root = Path(ref_panel_spec.r2_dir)
    paths: list[Path] = []
    if root.is_dir():
        paths.extend(sorted(root.glob("chr*_r2.parquet")))
        for build in ("hg19", "hg38"):
            child = root / build
            if child.is_dir():
                paths.extend(sorted(child.glob("chr*_r2.parquet")))
    builds = {_read_r2_sorted_by_build(path) for path in paths}
    builds.discard(None)
    if len(builds) == 1:
        return next(iter(builds))
    return None


def _read_r2_sorted_by_build(path: Path) -> str | None:
    """Read ``ldsc:sorted_by_build`` from a parquet R2 schema."""
    try:
        import pyarrow.parquet as pq
    except ImportError:
        return None
    raw_meta = pq.read_schema(str(path)).metadata or {}
    raw_build = raw_meta.get(b"ldsc:sorted_by_build")
    if raw_build is None:
        return None
    build = raw_build.decode("utf-8")
    return build if build in {"hg19", "hg38"} else None


def _infer_plink_ref_panel_build(ref_panel_spec: RefPanelConfig) -> str | None:
    """Infer reference-panel build from PLINK BIM coordinates."""
    if ref_panel_spec.plink_prefix is None:
        return None
    prefixes = resolve_plink_prefix_group((ref_panel_spec.plink_prefix,), allow_chromosome_suite=True)
    frames: list[pd.DataFrame] = []
    for prefix in prefixes:
        frame = pd.read_csv(
            prefix + ".bim",
            sep=r"\s+",
            header=None,
            usecols=[0, 3],
            names=["CHR", "POS"],
        )
        frames.append(frame)
    if not frames:
        return None
    sample = pd.concat(frames, ignore_index=True)
    return resolve_genome_build("auto", "chr_pos", sample, context="reference-panel SNP restriction", logger=LOGGER)


def _chrom_sort_key(chrom: str) -> tuple[int, str]:
    """Return a stable chromosome sort key matching the package-wide ordering."""
    return chrom_sort_key(chrom)


def _read_metadata_sidecar_identity(path: str | Path) -> dict[str, object] | None:
    """Read leading ``# ldsc:*`` identity metadata from a runtime sidecar."""
    opener = gzip.open if str(path).endswith(".gz") else open
    raw: dict[str, str] = {}
    with opener(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith("#"):
                break
            marker = "# ldsc:"
            if not line.startswith(marker):
                continue
            key, sep, value = line[len(marker) :].strip().partition("=")
            if sep != "":
                raw[key] = value
    if not raw:
        return None
    required = {"artifact_type", "snp_identifier", "genome_build"}
    if not required.issubset(raw):
        raise LDSCInputError(
            f"Reference-panel metadata sidecar '{path}' is missing LDSC identity "
            "metadata. Most likely it was written by an older LDSC version or edited. "
            "Regenerate the reference panel with the current `ldsc build-r2-panel`. "
            f"Other causes & fixes: {_REF_PANEL_ARTIFACT_DOC}"
        )
    return {
        "artifact_type": raw["artifact_type"],
        "snp_identifier": raw["snp_identifier"],
        "genome_build": raw["genome_build"],
    }


def _validate_metadata_identity_matches_config(metadata: dict[str, object], global_config: GlobalConfig) -> None:
    """Reject package metadata sidecars built under different identity assumptions."""
    artifact_config = GlobalConfig(
        snp_identifier=str(metadata["snp_identifier"]),
        genome_build=metadata.get("genome_build"),
    )
    validate_config_compatibility(
        artifact_config,
        global_config,
        context="runtime config and reference-panel metadata artifact",
    )


def _read_metadata_table(
    path: str | Path,
    chrom: str | None,
    global_config: GlobalConfig,
    *,
    require_identity_metadata: bool = False,
) -> pd.DataFrame:
    """Read one reference-panel metadata table into the normalized column set."""
    snp_identifier = normalize_snp_identifier_mode(global_config.snp_identifier)
    assert global_config.genome_build in {"hg19", "hg38", None}, (
        f"genome_build reached kernel as {global_config.genome_build!r}; "
        "should have been resolved at workflow entry."
    )

    identity_metadata = _read_metadata_sidecar_identity(path)
    if require_identity_metadata and identity_metadata is None:
        raise LDSCInputError(
            f"Reference-panel metadata sidecar '{path}' has no LDSC identity metadata, "
            "but its paired R2 parquet requires package provenance. Most likely the "
            "sidecar was written by an older LDSC version or copied from an external "
            "source. Regenerate the reference panel with the current `ldsc build-r2-panel`. "
            f"Other causes & fixes: {_REF_PANEL_ARTIFACT_DOC}"
        )
    if identity_metadata is not None:
        validate_identity_artifact_metadata(identity_metadata, expected_artifact_type="ref_panel_metadata")
        _validate_metadata_identity_matches_config(identity_metadata, global_config)

    df = pd.read_csv(
        path,
        sep=r"\s+",
        compression="gzip" if str(path).endswith(".gz") else None,
        comment="#",
    )
    context = str(path)

    chr_col = None
    pos_col = None
    snp_col = None
    try:
        chr_col = resolve_required_column(df.columns, REFERENCE_METADATA_SPEC_MAP["CHR"], context=context)
        pos_col = resolve_required_column(df.columns, REFERENCE_METADATA_SPEC_MAP["POS"], context=context)
    except (ValueError, LDSCInputError):
        # Metadata may be rsid-only; validate the required identifier columns
        # after probing both coordinate and SNP schemas.
        pass
    try:
        snp_col = resolve_required_column(df.columns, REFERENCE_METADATA_SPEC_MAP["SNP"], context=context)
    except (ValueError, LDSCInputError):
        # Metadata may be coordinate-only; mode-specific validation below
        # decides whether the missing SNP column is an error.
        pass

    family = identity_mode_family(snp_identifier)
    if family == "rsid" and snp_col is None:
        raise LDSCInputError(
            f"Reference-panel metadata sidecar '{path}' is incompatible with rsID-family "
            "mode: no SNP column was found. Most likely the sidecar uses CHR/POS-only "
            "identifiers or an unrecognized SNP header. Regenerate the panel or run with "
            "a chr_pos-family SNP identifier."
        )
    if family == "chr_pos" and (chr_col is None or pos_col is None):
        raise LDSCInputError(
            f"Reference-panel metadata sidecar '{path}' is incompatible with chr_pos-family "
            "mode: CHR and POS columns are required. Most likely the sidecar uses rsID-only "
            "identifiers or unrecognized CHR/POS headers. Regenerate the panel or run with "
            "an rsID-family SNP identifier."
        )

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
    if maf_col is None:
        raise LDSCInputError(
            f"Reference-panel metadata sidecar '{path}' is missing MAF metadata. MAF is "
            "required for LD-score calculation. Most likely the sidecar was produced "
            "without allele-frequency values. Regenerate the panel with MAF metadata."
        )
    maf = pd.to_numeric(df[maf_col], errors="coerce").astype(float)
    out["MAF"] = pd.Series(maf).map(lambda value: value if pd.isna(value) else min(value, 1.0 - value))

    a1_col = resolve_optional_column(df.columns, A1_COLUMN_SPEC, context=context)
    a2_col = resolve_optional_column(df.columns, A2_COLUMN_SPEC, context=context)
    if (a1_col is None) != (a2_col is None):
        raise LDSCInputError(
            f"Reference-panel metadata sidecar '{path}' has only one allele column. "
            "Most likely A1 or A2 is missing from the sidecar. Provide both allele "
            "columns, or regenerate the panel in a base SNP identifier mode."
        )
    if is_allele_aware_mode(snp_identifier) and (a1_col is None or a2_col is None):
        raise LDSCInputError(
            f"Reference-panel metadata sidecar '{path}' is malformed for allele-aware "
            "SNP identity: A1/A2 columns are required. Most likely the panel was built "
            "without allele metadata. Regenerate with A1/A2 metadata or use a base SNP "
            "identifier mode."
        )
    if a1_col is not None and a2_col is not None:
        out["A1"] = df[a1_col].astype(str)
        out["A2"] = df[a2_col].astype(str)

    if chrom is not None and "CHR" in out.columns:
        out = out.loc[out["CHR"] == normalize_chromosome(chrom, context=context)].reset_index(drop=True)
    return out


def _prepare_window(metadata, chrom, config, snp_identifier):
    """Validate retained reference metadata and derive one LD-window traversal."""
    kernel_ldscore.validate_retained_identifier_uniqueness(metadata, snp_identifier, chrom)
    kernel_ldscore.require_reference_maf(metadata, chrom)
    kernel_ldscore.validate_window_positions_sorted(metadata, chrom)
    coords, distance = kernel_ldscore.build_window_coordinates(metadata, config)
    block_left = kernel_ldscore.get_block_lefts(coords, distance)
    kernel_ldscore.check_whole_chromosome_window(
        block_left, SimpleNamespace(yes_really=config.whole_chromosome_ok), chrom,
    )
    return block_left


def _resolve_genetic_map(panel, chrom, metadata):
    """Resolve the configured PLINK map against the prepared reference coordinates."""
    spec = panel.spec
    if not (spec.genetic_map_hg19_sources or spec.genetic_map_hg38_sources):
        return None
    if spec.backend != "plink":
        LOGGER.warning("Ignoring --genetic-map-*-sources for the parquet R2 reference panel: CM is taken from the panel metadata sidecar (authoritative). Genetic-map flags apply only to PLINK panels.")
        return None
    build = normalize_genome_build(panel.global_config.genome_build)
    if build == "auto":
        build = resolve_genome_build(build, panel.global_config.snp_identifier, metadata[["CHR", "POS"]],
                                     context=f"PLINK panel chromosome {chrom} for genetic-map build selection")
    if build not in ("hg19", "hg38"):
        raise LDSCInputError("ldscore could not determine the genome build needed to select a genetic map for the PLINK panel. Pass --genome-build hg19 or hg38.")
    sources = spec.genetic_map_hg38_sources if build == "hg38" else spec.genetic_map_hg19_sources
    if sources is None:
        raise LDSCInputError(f"ldscore needs a `--genetic-map-{build}-sources` file to derive CM for the PLINK panel (resolved build {build}), but it was not supplied.")
    from .ref_panel_builder import load_genetic_map_group
    return load_genetic_map_group(split_cli_path_tokens(sources))


def _align_annotations(bundle, reference, chrom, global_config, backend):
    """Align annotation rows once to retained reference identities and recover alleles."""
    if is_allele_aware_mode(global_config.snp_identifier):
        reference = clean_identity_artifact_table(
            reference, global_config.snp_identifier,
            context=f"reference panel chromosome {chrom}",
            stage="annotation_reference_identity_cleanup", logger=LOGGER,
        ).cleaned
    mode = _annotation_reference_match_mode(bundle.metadata, global_config.snp_identifier)
    reference_keys = build_snp_id_series(reference, mode)
    annotation_keys = build_snp_id_series(bundle.metadata, mode)
    keep = annotation_keys.isin(set(reference_keys))
    if not keep.any():
        raise EmptyReferenceSNPs(
            f"ldscore retained no annotation SNPs on chromosome {chrom} after intersecting with the {backend} reference panel. "
            "Most likely the annotation SNP identifiers, genome build, or allele-aware identifier mode do not match the reference panel. "
            "Use annotation and reference-panel artifacts built with the same SNP identifier mode and genome build. "
            f"Other causes & fixes: {kernel_ldscore._LDSCORE_INTERSECTION_DOC}"
        )
    metadata = _annotation_metadata_with_reference_alleles(
        annotation_metadata=bundle.metadata.loc[keep].reset_index(drop=True),
        annotation_keys=annotation_keys.loc[keep].reset_index(drop=True),
        reference_metadata=reference, reference_keys=reference_keys,
        snp_identifier=global_config.snp_identifier,
    )
    return kernel_ldscore.AnnotationBundle(metadata,
        MappedAnnotations(bundle.annotations, np.flatnonzero(keep.to_numpy()),
                          tuple(bundle.baseline_columns + bundle.query_columns)),
        list(bundle.baseline_columns), list(bundle.query_columns))


def _annotation_reference_match_mode(
    annotation_metadata: pd.DataFrame,
    snp_identifier: str,
) -> str:
    """Return the identity mode usable by both annotation and reference metadata."""
    mode = normalize_snp_identifier_mode(snp_identifier)
    if is_allele_aware_mode(mode) and not {"A1", "A2"}.issubset(annotation_metadata.columns):
        return identity_base_mode(mode)
    return mode


def _annotation_metadata_with_reference_alleles(
    *,
    annotation_metadata: pd.DataFrame,
    annotation_keys: pd.Series,
    reference_metadata: pd.DataFrame,
    reference_keys: pd.Series,
    snp_identifier: str,
) -> pd.DataFrame:
    """Fill missing annotation alleles from unambiguous retained reference metadata."""
    mode = normalize_snp_identifier_mode(snp_identifier)
    if not is_allele_aware_mode(mode) or {"A1", "A2"}.issubset(annotation_metadata.columns):
        return annotation_metadata
    if not {"A1", "A2"}.issubset(reference_metadata.columns):
        return annotation_metadata

    reference_lookup = reference_metadata.loc[:, ["A1", "A2"]].copy()
    reference_lookup["_match_key"] = reference_keys.to_numpy()
    reference_lookup = reference_lookup.loc[reference_lookup["_match_key"].notna()].copy()
    duplicate_mask = reference_lookup["_match_key"].duplicated(keep=False)
    if bool(duplicate_mask.any()):
        duplicate_key = str(reference_lookup.loc[duplicate_mask, "_match_key"].iloc[0])
        raise LDSCInputError(
            "ldscore could not infer missing annotation alleles from reference metadata. "
            f"Base identity {duplicate_key!r} is not unique under SNP identifier mode "
            f"'{mode}'. Most likely the reference panel contains duplicate base SNP "
            "identities after dropping alleles. Provide allele-aware annotation columns "
            "A1/A2, or regenerate the reference panel with unique retained SNP identities."
        )
    reference_lookup = reference_lookup.set_index("_match_key")
    missing = ~annotation_keys.isin(reference_lookup.index)
    if bool(missing.any()):
        missing_key = str(annotation_keys.loc[missing].iloc[0])
        raise LDSCInputError(
            "ldscore could not infer missing annotation alleles from reference metadata. "
            f"Retained annotation SNP {missing_key!r} has no matching A1/A2 values in "
            "the reference metadata. Most likely the annotation and reference panel were "
            "matched with base chr_pos/rsID identities but the allele metadata is incomplete. "
            "Provide A1/A2 in the annotation input or regenerate the reference panel with allele columns."
        )

    enriched = annotation_metadata.copy()
    enriched["A1"] = annotation_keys.map(reference_lookup["A1"]).astype(str)
    enriched["A2"] = annotation_keys.map(reference_lookup["A2"]).astype(str)
    return enriched
