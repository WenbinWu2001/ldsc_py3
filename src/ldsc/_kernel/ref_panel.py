"""Reference-panel configuration and workflow-to-kernel backend adapters.

This module prepares chromosome-scoped reference data for LD-score projection.
The PLINK and parquet adapters own identity and SNP restrictions, sample/MAF
selection, annotation alignment, reader policy, LD windows, and file lifetime.
Numerical kernels consume the resulting aligned state without resolving paths
or reopening a panel. Metadata inspection returns caller-owned chromosome tables without caching them.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import replace
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
from ..errors import EmptyReferenceSNPs, LDSCConfigError, LDSCInputError, LDSCUsageError
from ..genome_build_inference import resolve_genome_build, validate_auto_genome_build_mode
from ..path_resolution import inspect_plink_inputs, split_cli_path_tokens
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
_REF_PANEL_ARTIFACT_DOC = "docs/user-guide/troubleshooting.md#build-r2-panel-reference-panel-artifact-is-incompatible"


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


def _artifact_global_config(metadata, global_config, source):
    """Validate declared coordinates before constructing an artifact snapshot."""
    build = metadata.get("genome_build")
    if identity_mode_family(global_config.snp_identifier) == "chr_pos":
        from .._genome_build_policy import GenomeBuildCheck
        check = GenomeBuildCheck(global_config.genome_build, context="reference-panel validation", logger=LOGGER)
        if build is not None:
            check.declared("reference artifact", source, build)
        build = check.resolve()
    return GlobalConfig(snp_identifier=str(metadata["snp_identifier"]), genome_build=build)


def _validate_r2_identity_matches_config(metadata: dict[str, object], global_config: GlobalConfig, *, source="R2 artifact") -> None:
    """Reject package R2 artifacts built under different identity assumptions."""
    artifact_config = _artifact_global_config(metadata, global_config, source)
    validate_config_compatibility(
        global_config,
        artifact_config,
        context="runtime config and parquet R2 artifact",
    )


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
            SNP metadata in genomic order and a selected-read annotation source
            in baseline-then-query order. The source may access packed binary
            and dense continuous stores. Input tables are not modified.
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
            the validated R2 parquet contract are resolved together. Use the result as a
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
    repeated discovery for workflow-validated inputs. Direct Python use resolves
    and caches the same mapping lazily from BIM contents.
    """

    def __init__(self, global_config: GlobalConfig, spec: RefPanelConfig, *,
                 chromosome_prefixes: dict[str, str] | None = None) -> None:
        """Store panel configuration and optional validated chromosome routing."""
        super().__init__(global_config, spec)
        self._chromosome_prefixes = None if chromosome_prefixes is None else dict(chromosome_prefixes)

    def available_chromosomes(self) -> list[str]:
        """List normalized chromosomes present in the resolved PLINK inputs."""
        return sorted(self._resolved_prefixes(), key=chrom_sort_key)

    def _resolved_prefixes(self):
        """Cache workflow resolution once for direct Python panel use."""
        if self._chromosome_prefixes is None:
            plink = inspect_plink_inputs(self.spec.plink_prefix)
            plink.require_valid()
            self._chromosome_prefixes = plink.chromosome_prefixes
        return self._chromosome_prefixes

    def _load_source(self, chrom: str):
        """Resolve physical BED rows and apply panel identity/SNP restrictions once."""
        prefix = self._bim_prefixes(chrom)[0]
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
        prefixes = self._resolved_prefixes()
        if chrom is None:
            return list(dict.fromkeys(prefixes.values()))
        chrom = normalize_chromosome(chrom)
        if chrom not in prefixes:
            raise LDSCInputError(f"No validated PLINK trio contains chromosome {chrom}.")
        return [prefixes[chrom]]


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
            return sorted(set(self.spec.chromosomes), key=chrom_sort_key)

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
                return sorted(chromosomes, key=chrom_sort_key)
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

        # A sidecar without its R2 parquet is an incomplete panel.
        self.resolve_r2_paths(chrom)
        metadata_paths = self.resolve_metadata_paths(chrom)
        if metadata_paths:
            frames = [
                _read_metadata_table(path, chrom=chrom, global_config=self.global_config)
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
        metadata = self._validate_metadata(metadata, chrom)
        return metadata

    def build_reader(self, chrom: str, metadata: pd.DataFrame | None = None):
        """
        Build an owned streaming reader over one chromosome's validated R2 parquet.

        The parquet must satisfy the ``build-r2-panel`` metadata contract and
        record the active SNP identifier and genome build; its values are read
        as already unbiased. ``prepare_chromosome`` uses this same factory;
        direct callers must close the returned reader.
        """
        metadata = metadata if metadata is not None else self.load_metadata(chrom)
        paths = self.resolve_r2_paths(chrom)
        if paths:
            contract = kernel_ldscore.read_r2_parquet_contract(paths[0])
            _validate_r2_identity_matches_config(
                {"snp_identifier": contract.snp_identifier, "genome_build": contract.genome_build},
                self.global_config, source=paths[0],
            )
        return kernel_ldscore.SortedR2BlockReader(
            paths=paths,
            chrom=chrom,
            metadata=metadata,
            identifier_mode=self.global_config.snp_identifier,
            genome_build=self.global_config.genome_build,
        )

    def prepare_chromosome(self, chrom, annotations, config, *, genetic_map=None):
        """Prepare a sidecar-aligned streaming reader over a contract-validated R2 parquet."""
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
        parquet_paths = self.resolve_r2_paths(chrom)
        kernel_ldscore.validate_ldscore_window_within_r2_panel_window(config, parquet_paths=parquet_paths, chrom=chrom)
        kernel_ldscore.reject_truncated_r2_panel(parquet_paths, chrom)
        mapped_annotations = aligned.annotations
        reader = self.build_reader(chrom, metadata=metadata)
        if np.any(np.diff(reader.panel_rows) < 0):
            # Same-position SNPs listed in another order than the sidecar: compute in
            # panel-row order, which --ld-wind-snps needs. Annotation rows move with
            # their SNPs, and outputs are re-sorted by CHR, POS, and SNP afterwards.
            order = np.argsort(reader.panel_rows, kind="stable")
            reader.close()
            metadata = metadata.iloc[order].reset_index(drop=True)
            mapped_annotations = replace(mapped_annotations, row_indices=mapped_annotations.row_indices[order])
            reader = self.build_reader(chrom, metadata=metadata)
        try:
            # --ld-wind-snps counts panel rows, so a window no wider than the
            # panel's stays inside its stored pairs when this run keeps fewer SNPs.
            block_left = _prepare_window(metadata, chrom, config, self.global_config.snp_identifier,
                                         snp_positions=reader.panel_rows)
        except BaseException:
            reader.close()
            raise
        return kernel_ldscore.PreparedChromosome(
            backend="parquet_r2", reader=reader, metadata=metadata,
            annotations=mapped_annotations, block_left=block_left,
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
        key=lambda path: chrom_sort_key(_chromosome_from_ref_panel_r2_path(path) or path.name),
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
        """Validate direct Python build declarations before backend loading."""
        mode = normalize_snp_identifier_mode(self.global_config.snp_identifier)
        if identity_mode_family(mode) != "chr_pos":
            return self.global_config
        resolved_build = _infer_restricted_ref_panel_build(ref_panel_spec, self.global_config.genome_build)
        return GlobalConfig(
            snp_identifier=mode,
            genome_build=resolved_build,
            log_level=self.global_config.log_level,
        )


def _infer_restricted_ref_panel_build(ref_panel_spec: RefPanelConfig, requested="auto") -> str:
    """Resolve all panel shards and coordinate restrictions before filtering."""
    from .._genome_build_policy import GenomeBuildCheck

    check = GenomeBuildCheck(requested, context='reference-panel loading', logger=LOGGER)
    if ref_panel_spec.backend == 'plink':
        check.plink('--plink-prefix', ref_panel_spec.plink_prefix)
    elif ref_panel_spec.r2_dir is not None:
        check.r2('--r2-dir', ref_panel_spec.r2_dir)
    check.paths('--ref-panel-snps-file', ref_panel_spec.ref_panel_snps_file, restriction_mode='chr_pos')
    return check.resolve()


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


def _validate_metadata_identity_matches_config(metadata: dict[str, object], global_config: GlobalConfig, *, source="reference metadata") -> None:
    """Reject package metadata sidecars built under different identity assumptions."""
    artifact_config = _artifact_global_config(metadata, global_config, source)
    validate_config_compatibility(
        artifact_config,
        global_config,
        context="runtime config and reference-panel metadata artifact",
    )


def _read_metadata_table(
    path: str | Path,
    chrom: str | None,
    global_config: GlobalConfig,
) -> pd.DataFrame:
    """Read one package-written reference-panel sidecar into the normalized column set."""
    snp_identifier = normalize_snp_identifier_mode(global_config.snp_identifier)
    assert global_config.genome_build in {"hg19", "hg38", None}, (
        f"genome_build reached kernel as {global_config.genome_build!r}; "
        "should have been resolved at workflow entry."
    )

    identity_metadata = _read_metadata_sidecar_identity(path)
    if identity_metadata is None:
        raise LDSCInputError(
            f"Reference-panel metadata sidecar '{path}' has no LDSC identity metadata. "
            "Most likely it was written by an older LDSC version or copied from an external "
            "source. Regenerate the reference panel with the current `ldsc build-r2-panel`. "
            f"Other causes & fixes: {_REF_PANEL_ARTIFACT_DOC}"
        )
    validate_identity_artifact_metadata(identity_metadata, expected_artifact_type="ref_panel_metadata")
    _validate_metadata_identity_matches_config(identity_metadata, global_config, source=path)

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
    # NaN fails both comparisons, so this also catches missing and non-numeric values.
    invalid_maf = ~((maf >= 0.0) & (maf <= 0.5))
    if invalid_maf.any():
        labels = out["SNP"] if "SNP" in out else out["CHR"].astype(str) + ":" + out["POS"].astype(str)
        examples = ", ".join(labels[invalid_maf].head(5))
        raise LDSCInputError(
            f"Reference-panel metadata sidecar '{path}' has missing, non-numeric, or out-of-range MAF for "
            f"{int(invalid_maf.sum())} of {len(maf)} SNP rows (first: {examples}). "
            "`ldsc build-r2-panel` records MAF as the minor-allele (A1) frequency, in [0, 0.5], "
            "for every row, so the sidecar was most likely edited or produced by another tool. "
            "Regenerate the reference panel with the current `ldsc build-r2-panel`. "
            f"Other causes & fixes: {_REF_PANEL_ARTIFACT_DOC}"
        )
    out["MAF"] = maf

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


def _prepare_window(metadata, chrom, config, snp_identifier, *, snp_positions=None):
    """Validate retained reference metadata and derive one LD-window traversal."""
    kernel_ldscore.validate_retained_identifier_uniqueness(metadata, snp_identifier, chrom)
    kernel_ldscore.validate_window_positions_sorted(metadata, chrom)
    coords, distance = kernel_ldscore.build_window_coordinates(metadata, config, snp_positions=snp_positions)
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
    from .._genome_build_policy import GenomeBuildCheck
    check = GenomeBuildCheck(build, context="PLINK genetic map", logger=LOGGER)
    check.paths(f"--genetic-map-{build}-sources", sources)
    check.resolve()
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
