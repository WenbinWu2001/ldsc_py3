"""Reference-panel configuration and workflow-to-kernel backend adapters.

This module prepares chromosome-scoped reference data for LD-score projection.
The PLINK and R2 matrix adapters own identity and SNP restrictions, sample/MAF
selection, annotation alignment, reader policy, LD windows, and file lifetime.
Numerical kernels consume the resulting aligned state without resolving paths
or reopening a panel. Metadata inspection returns caller-owned chromosome tables without caching them.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import replace
from types import SimpleNamespace
import logging
from typing import Any

import numpy as np
import pandas as pd

from .._core.coordinates import CHR_POS_KEY_COLUMN, build_chr_pos_key_frame
from .._core.chromosomes import chrom_sort_key, normalize_chromosome
from .._core.columns import normalize_genome_build
from ..config import GlobalConfig, LDScoreConfig, RefPanelConfig, validate_config_compatibility
from ..errors import EmptyReferenceSNPs, LDSCConfigError, LDSCInputError, LDSCInternalError, LDSCUsageError
from .._genome.build_inference import resolve_genome_build, validate_auto_genome_build_mode
from .._runtime.paths import inspect_plink_inputs, split_cli_path_tokens
from .._core import plink_files as legacy_parse
from .._core.identifiers import identifier_keys, validate_retained_identifier_uniqueness
from .funnel import ANNOTATION_SNP_GRID, REFERENCE_PANEL, ReferenceSNPFunnel
from .plink import PlinkBEDFile
from .prepared import ChromosomeAnnotations, PreparedChromosome
from .r2_format import (
    R2_MATRIX_TROUBLESHOOTING,
    R2MatrixChromosome,
    r2_matrix_chromosomes,
    resolve_r2_build_dir,
)
from .windows import (
    assert_cm_usable,
    check_whole_chromosome_window,
    get_block_lefts,
    validate_window_positions_sorted,
    window_coordinates,
)
from .prepared import MappedAnnotations
from .._core.identifiers import (
    build_snp_id_series,
    normalize_snp_identifier_mode,
    read_snp_restriction_keys,
    validate_unique_snp_ids,
)
from .._core.snp_identity import (
    clean_identity_artifact_table,
    empty_identity_drop_frame,
    identity_base_mode,
    identity_mode_family,
    is_allele_aware_mode,
    restriction_membership_mask,
)

LOGGER = logging.getLogger("LDSC.ref_panel")


def _snp_id_series_for_matching(metadata: pd.DataFrame, snp_identifier: str, *, context: str) -> pd.Series:
    """Build SNP identities for a match boundary, dropping missing CHR/POS in `chr_pos` mode."""
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
        from .._genome.build_policy import GenomeBuildCheck
        check = GenomeBuildCheck(global_config.genome_build, context="reference-panel validation", logger=LOGGER)
        if build is not None:
            check.declared("reference artifact", source, build)
        build = check.resolve()
    mode = str(metadata["snp_identifier"])
    # Artifacts always record a build; rsID-family identity ignores it, so drop it
    # here instead of triggering GlobalConfig's warning about a user-supplied build.
    if identity_mode_family(mode) == "rsid":
        build = None
    return GlobalConfig(snp_identifier=mode, genome_build=build)


def _validate_r2_identity_matches_config(metadata: dict[str, object], global_config: GlobalConfig, *, source="R2 matrix") -> None:
    """Reject R2 matrices built under different identity assumptions."""
    artifact_config = _artifact_global_config(metadata, global_config, source)
    validate_config_compatibility(
        global_config,
        artifact_config,
        context="runtime config and R2 matrix",
    )


class RefPanel(ABC):
    """Abstract chromosome-scoped reference-panel interface."""
    _backend: str

    def __init__(self, global_config: GlobalConfig, spec: RefPanelConfig) -> None:
        """Store shared configuration without chromosome payload caches."""
        validate_auto_genome_build_mode(global_config.snp_identifier, global_config.genome_build)
        if spec.backend != self._backend:
            raise LDSCConfigError(
                f"{type(self).__name__} cannot use RefPanelConfig(backend={spec.backend!r}). "
                "Most likely the configuration was written for the other reference-panel form. "
                f"Pass RefPanelConfig(backend={self._backend!r}, ...), or let RefPanelLoader.load choose "
                "the adapter from the configuration."
            )
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
        self, chrom: str, annotations: ChromosomeAnnotations,
        config: LDScoreConfig, *, genetic_map: pd.DataFrame | None = None,
    ) -> PreparedChromosome:
        """Prepare aligned reference data and an owned reader for LD-score computation.

        Parameters
        ----------
        chrom : str
            Chromosome label to resolve in this panel.
        annotations : ldsc._refpanel.prepared.ChromosomeAnnotations
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
        ldsc._refpanel.prepared.PreparedChromosome
            Reference metadata, a selected-read annotation mapping, window bounds, and
            an owned reader over the same retained SNP rows. Panel identity,
            SNP/sample/MAF restrictions, authoritative reference CM/MAF, and
            the validated R2 Parquet file contract are resolved together. Use the result as a
            context manager to close the reader on success or failure.
            ``reference_snp_funnel`` counts what each step removed.

        Raises
        ------
        ldsc.errors.EmptyReferenceSNPs
            No reference SNP is left on the chromosome. The error's
            ``reference_snp_funnel`` shows which step removed the last SNPs.

        Notes
        -----
        The annotation rows are the annotation SNP grid: panel SNPs that pass
        the panel's filters but have no annotation row are excluded, and the
        funnel records how many. Preparation reads inputs but writes no
        artifacts. Metadata inspection through ``load_metadata`` reads
        independently and retains no cache.
        """
        raise NotImplementedError

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
    _backend = "plink"

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
        """Resolve physical BED rows and apply panel identity/SNP restrictions once.

        Returns the prefix, BIM, FAM, restricted metadata, identity drops and
        the chromosome's reference-SNP funnel so far.
        """
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
        funnel = ReferenceSNPFunnel(normalize_chromosome(chrom), "plink")
        funnel.record(REFERENCE_PANEL, "panel_rows", len(metadata))
        drops = empty_identity_drop_frame()
        if self.global_config.snp_identifier in {"rsid", "chr_pos"}:
            cleanup = clean_identity_artifact_table(
                metadata, self.global_config.snp_identifier,
                context=f"PLINK reference-panel metadata chromosome {chrom}",
                stage="plink_reference_identity_cleanup", logger=LOGGER,
            )
            funnel.record(REFERENCE_PANEL, "identity_cleanup", len(metadata) - len(cleanup.cleaned))
            metadata, drops = cleanup.cleaned, cleanup.dropped
            if metadata.empty:
                raise LDSCInputError(
                    f"Reference-panel loading retained no PLINK rows on chromosome {chrom} after duplicate SNP identity cleanup."
                )
        if self.spec.ref_panel_snps_file is not None:
            before = len(metadata)
            metadata = self._apply_snp_restriction(metadata)
            funnel.record(REFERENCE_PANEL, "ref_panel_snps_file", before - len(metadata))
        metadata = self._validate_metadata(metadata, chrom)
        return prefix, bim, fam, metadata, drops, funnel

    def _load_genotypes(self, chrom, source, annotations=None, *, keep_indivs=None, maf_min=None):
        """Open one BED and align its retained rows to annotations without another read."""
        prefix, bim, fam, reference, drops, funnel = source
        if annotations is None:
            annotations = ChromosomeAnnotations(
                reference.drop(columns="_raw_index"), pd.DataFrame(index=reference.index), [], [],
            )
        aligned = _align_annotations(annotations, reference, chrom, self.global_config, funnel)
        mode = self.global_config.snp_identifier
        raw_keys = identifier_keys(reference, mode)
        keys = identifier_keys(aligned.metadata, mode)
        physical_rows = pd.Series(reference["_raw_index"].to_numpy(), index=raw_keys).loc[keys].tolist()
        if keep_indivs is None:
            keep_indivs = resolve_keep_individuals(self.spec.keep_indivs_file, fam)
        maf_min = self.spec.maf_min if maf_min is None else maf_min
        try:
            geno = PlinkBEDFile(
                prefix + ".bed", len(fam.IDList), bim, keep_snps=physical_rows,
                keep_indivs=keep_indivs, mafMin=maf_min,
            )
        except EmptyReferenceSNPs as exc:
            _record_genotype_qc(funnel, exc.genotype_qc_removed, exc.maf_removed, maf_min)
            exc.reference_snp_funnel = funnel.frame()
            raise
        _record_genotype_qc(funnel, geno.genotype_qc_removed, geno.maf_removed, maf_min)
        try:
            metadata = pd.DataFrame(geno.df, columns=geno.colnames).rename(columns={"BP": "POS"})
            metadata["CHR"] = metadata["CHR"].map(normalize_chromosome)
            metadata["SNP"] = metadata["SNP"].astype(str)
            metadata["POS"] = pd.to_numeric(metadata["POS"], errors="raise").astype(np.int64)
            metadata["CM"] = pd.to_numeric(metadata["CM"], errors="coerce")
            metadata["MAF"] = pd.to_numeric(metadata["MAF"], errors="coerce")
            metadata = metadata.merge(reference[["CHR", "SNP", "POS", "A1", "A2"]],
                                      on=["CHR", "SNP", "POS"], how="left", sort=False)
            retained_keys = identifier_keys(metadata, mode)
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
            prefix, bim, fam, metadata, drops, funnel = source
            keys = _snp_id_series_for_matching(metadata, self.global_config.snp_identifier,
                                              context=f"PLINK chromosome {chrom}")
            source = prefix, bim, fam, metadata.loc[keys.isin(keep_snps)].reset_index(drop=True), drops, funnel
        return self._load_genotypes(chrom, source, keep_indivs=keep_indivs, maf_min=maf_min)[0]

    def prepare_chromosome(self, chrom, annotations, config, *, genetic_map=None):
        """Prepare one genotype-filtered chromosome with aligned annotations and LD windows."""
        source = self._load_source(chrom)
        geno, metadata, values, drops, input_rows = self._load_genotypes(chrom, source, annotations)
        funnel = source[-1].frame()
        try:
            if config.ld_wind_cm is not None:
                if genetic_map is None:
                    genetic_map = _resolve_genetic_map(self, chrom, metadata)
                if genetic_map is not None:
                    from .build_kernel import interpolate_genetic_map_cm
                    metadata["CM"] = interpolate_genetic_map_cm(
                        normalize_chromosome(chrom), metadata["POS"].to_numpy(dtype=np.int64), genetic_map,
                    )
                else:
                    assert_cm_usable(metadata["CM"], chrom)
            block_left = _prepare_window(metadata, chrom, config, self.global_config.snp_identifier)
            return PreparedChromosome(
                backend="plink", reader=geno, metadata=metadata,
                annotations=values, block_left=block_left,
                baseline_columns=list(annotations.baseline_columns), query_columns=list(annotations.query_columns),
                reference_rows_before_genotype_qc=input_rows,
                genotype_qc_removed=geno.genotype_qc_removed, maf_removed=geno.maf_removed,
                selected_individual_count=int(geno.n),
                cm_source="explicit_genetic_map" if genetic_map is not None else "bim_cm", identity_drops=drops,
                reference_snp_funnel=funnel,
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


class R2MatrixRefPanel(RefPanel):
    """
    Reference-panel adapter for R2 matrices.

    Each chromosome is read through one :class:`R2MatrixChromosome`: the
    sidecar (``chrN_meta.tsv.gz``) is mandatory and authoritative for SNP
    identity, CM and MAF, and the paired R2 Parquet file (``chrN_r2.parquet``)
    stores only sidecar-row indices, meaningless without the exact matching
    sidecar. ``load_metadata()`` applies ``ref_panel_snps_file`` and
    ``maf_min`` to the sidecar and returns the restricted table.
    """
    _backend = "r2_matrix"

    def available_chromosomes(self) -> list[str]:
        """List chromosomes from explicit config or from the files in the R2 matrix directory."""
        if self.spec.chromosomes is not None:
            return sorted(set(self.spec.chromosomes), key=chrom_sort_key)

        if self.spec.r2_matrix_dir is not None:
            chromosomes = r2_matrix_chromosomes(resolve_r2_build_dir(self.spec.r2_matrix_dir,
                                                                     self.global_config.genome_build))
            if chromosomes:
                return chromosomes
            raise LDSCInputError(
                f"Reference-panel loading found no `chr*_r2.parquet` files in R2 matrix "
                f"directory '{self.spec.r2_matrix_dir}'. Most likely the directory is not an "
                "`ldsc build-r2-matrix` output directory or the wrong genome-build child "
                "directory was selected. Pass the directory containing the R2 Parquet files. "
                f"Other causes & fixes: {R2_MATRIX_TROUBLESHOOTING['ldscore']}"
            )

        raise LDSCUsageError(
            "R2 matrix chromosome discovery requires explicit chromosomes or `r2_matrix_dir`. "
            "Most likely `RefPanelConfig(backend='r2_matrix')` was used without an R2 matrix "
            "directory. Pass `r2_matrix_dir` or an explicit chromosome list."
        )

    def open_chromosome(self, chrom: str) -> R2MatrixChromosome:
        """Open one chromosome of the configured R2 matrix for ``ldscore``.

        The R2 Parquet footer must record this panel's SNP identifier mode and
        genome build. The returned handle parses its sidecar once, on first use.
        """
        assert self.global_config.genome_build in {"hg19", "hg38", None}, (
            f"genome_build reached kernel as {self.global_config.genome_build!r}; "
            "should have been resolved at workflow entry."
        )
        if self.spec.r2_matrix_dir is None:
            raise LDSCUsageError(
                "R2 matrix loading requires `r2_matrix_dir`. Most likely "
                "`RefPanelConfig(backend='r2_matrix')` was used without an R2 matrix directory. "
                "Pass the directory containing `chrN_r2.parquet` files."
            )
        build_dir = resolve_r2_build_dir(self.spec.r2_matrix_dir, self.global_config.genome_build)
        chromosome = R2MatrixChromosome.open(build_dir, chrom, command="ldscore")
        _validate_r2_identity_matches_config(
            {"snp_identifier": chromosome.contract.snp_identifier, "genome_build": chromosome.contract.genome_build},
            self.global_config, source=str(chromosome.parquet_path),
        )
        return chromosome

    def load_metadata(self, chrom: str) -> pd.DataFrame:
        """Load, restrict, and validate metadata owned by this chromosome consumer."""
        return self._load_reference(self.open_chromosome(chrom))[0]

    def _load_reference(self, chromosome: R2MatrixChromosome) -> tuple[pd.DataFrame, ReferenceSNPFunnel]:
        """Return the restricted sidecar and the chromosome's reference-SNP funnel so far."""
        if self.spec.keep_indivs_file is not None:
            raise LDSCUsageError(
                "Reference-panel loading cannot apply `--keep-indivs-file` to an R2 matrix. "
                "Most likely individual filtering was combined with a precomputed "
                "R2 matrix. Remove the keep file or use a PLINK reference panel."
            )
        metadata = chromosome.sidecar()
        if metadata.empty:
            raise LDSCInputError(
                f"Reference-panel loading found no metadata rows for chromosome {chromosome.chrom} in "
                f"'{chromosome.sidecar_path}'. Most likely the R2 matrix was built from a panel with no SNPs on "
                "this chromosome. Leave the chromosome out of the analysis or rebuild the R2 matrix."
            )
        funnel = ReferenceSNPFunnel(chromosome.chrom, "r2_matrix")
        funnel.record(REFERENCE_PANEL, "panel_rows", len(metadata))
        if self.spec.ref_panel_snps_file is not None:
            before = len(metadata)
            metadata = self._apply_snp_restriction(metadata)
            funnel.record(REFERENCE_PANEL, "ref_panel_snps_file", before - len(metadata))
        if self.spec.maf_min is not None:
            before = len(metadata)
            metadata = self._apply_maf_filter(metadata, chromosome.chrom)
            funnel.record(REFERENCE_PANEL, "maf_floor", before - len(metadata))
        metadata = self._validate_metadata(metadata, chromosome.chrom)
        return metadata, funnel

    def build_reader(self, chrom: str, metadata: pd.DataFrame | None = None):
        """
        Build an owned streaming reader over one chromosome's validated R2 Parquet file.

        The pairs among the SNPs in ``metadata`` (default: the restricted
        sidecar from :meth:`load_metadata`) are streamed; values are read as
        already unbiased. Direct callers must close the returned reader.
        """
        chromosome = self.open_chromosome(chrom)
        if metadata is None:
            metadata = self._load_reference(chromosome)[0]
        return chromosome.pair_reader(metadata, self.global_config.snp_identifier, self.global_config.genome_build)

    def prepare_chromosome(self, chrom, annotations, config, *, genetic_map=None):
        """Prepare a sidecar-aligned streaming reader over a contract-validated R2 Parquet file."""
        chromosome = self.open_chromosome(chrom)
        chromosome.check_ld_window(config)
        chromosome.reject_truncation()
        reference, funnel = self._load_reference(chromosome)
        aligned = _align_annotations(annotations, reference, chrom, self.global_config, funnel)
        keys = identifier_keys(aligned.metadata, self.global_config.snp_identifier)
        reference_keys = identifier_keys(reference, self.global_config.snp_identifier)
        lookup = reference.set_index(reference_keys)
        metadata = aligned.metadata.copy()
        for column in ("CM", "MAF"):
            metadata[column] = lookup.loc[keys, column].to_numpy()
        if config.ld_wind_cm is not None:
            _resolve_genetic_map(self, chrom, reference)
        mapped_annotations = aligned.annotations
        mode, build = self.global_config.snp_identifier, self.global_config.genome_build
        reader = chromosome.pair_reader(metadata, mode, build)
        if np.any(np.diff(reader.panel_rows) < 0):
            # Same-position SNPs listed in another order than the sidecar: compute in
            # R2 matrix row order, which --ld-wind-snps needs. Annotation rows move with
            # their SNPs, and outputs are re-sorted by CHR, POS, and SNP afterwards.
            order = np.argsort(reader.panel_rows, kind="stable")
            reader.close()
            metadata = metadata.iloc[order].reset_index(drop=True)
            mapped_annotations = replace(mapped_annotations, row_indices=mapped_annotations.row_indices[order])
            reader = chromosome.pair_reader(metadata, mode, build)
        try:
            # --ld-wind-snps counts R2 matrix rows, so a window no wider than the
            # R2 matrix's stays inside its stored pairs when this run keeps fewer SNPs.
            block_left = _prepare_window(metadata, chrom, config, self.global_config.snp_identifier,
                                         snp_positions=reader.panel_rows)
        except BaseException:
            reader.close()
            raise
        return PreparedChromosome(
            backend="r2_matrix", reader=reader, metadata=metadata,
            annotations=mapped_annotations, block_left=block_left,
            baseline_columns=list(annotations.baseline_columns), query_columns=list(annotations.query_columns),
            cm_source="parquet_sidecar", reference_snp_funnel=funnel.frame(),
        )


class RefPanelLoader:
    """Instantiate the backend requested by a :class:`RefPanelConfig`."""
    def __init__(self, global_config: GlobalConfig) -> None:
        """Store the shared configuration used for every loaded panel."""
        self.global_config = global_config

    def load(self, ref_panel_spec: RefPanelConfig) -> RefPanel:
        """Return the adapter that ``ref_panel_spec.backend`` names."""
        global_config = self._global_config_for_spec(ref_panel_spec)
        adapter = PlinkRefPanel if ref_panel_spec.backend == "plink" else R2MatrixRefPanel
        return adapter(global_config, ref_panel_spec)

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
    from .._genome.build_policy import GenomeBuildCheck

    check = GenomeBuildCheck(requested, context='reference-panel loading', logger=LOGGER)
    if ref_panel_spec.backend == 'plink':
        check.plink('--plink-prefix', ref_panel_spec.plink_prefix)
    elif ref_panel_spec.r2_matrix_dir is not None:
        check.r2('--r2-matrix-dir', ref_panel_spec.r2_matrix_dir)
    check.paths('--ref-panel-snps-file', ref_panel_spec.ref_panel_snps_file, restriction_mode='chr_pos')
    return check.resolve()


def _prepare_window(metadata, chrom, config, snp_identifier, *, snp_positions=None):
    """Validate retained reference metadata and derive one LD-window traversal."""
    validate_retained_identifier_uniqueness(metadata, snp_identifier, chrom)
    validate_window_positions_sorted(metadata, chrom)
    if config.ld_wind_cm is not None and metadata["CM"].isna().any():
        raise LDSCInputError(
            "ldscore cannot use `--ld-wind-cm` because retained SNP metadata contains "
            "missing CM values. Most likely the annotation or frequency sidecar lacks "
            "genetic-map positions for at least one retained SNP. Provide complete CM "
            "metadata or use `--ld-wind-kb` / `--ld-wind-snps`."
        )
    coords, distance = window_coordinates(
        metadata["POS"], metadata["CM"] if config.ld_wind_cm is not None else None,
        ld_wind_snps=config.ld_wind_snps, ld_wind_kb=config.ld_wind_kb, ld_wind_cm=config.ld_wind_cm,
        snp_positions=snp_positions,
    )
    block_left = get_block_lefts(coords, distance)
    check_whole_chromosome_window(
        block_left, SimpleNamespace(yes_really=config.whole_chromosome_ok), chrom,
    )
    return block_left


def _resolve_genetic_map(panel, chrom, metadata):
    """Resolve the configured PLINK map against the prepared reference coordinates."""
    spec = panel.spec
    if not (spec.genetic_map_hg19_sources or spec.genetic_map_hg38_sources):
        return None
    if spec.backend != "plink":
        LOGGER.warning("Ignoring --genetic-map-*-sources for the R2 matrix: CM is taken from the panel metadata sidecar (authoritative). Genetic-map flags apply only to PLINK panels.")
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
    from .build_kernel import load_genetic_map_group
    from .._genome.build_policy import GenomeBuildCheck
    check = GenomeBuildCheck(build, context="PLINK genetic map", logger=LOGGER)
    check.paths(f"--genetic-map-{build}-sources", sources)
    check.resolve()
    return load_genetic_map_group(split_cli_path_tokens(sources))


def _record_genotype_qc(funnel, genotype_qc_removed, maf_removed, maf_min):
    """Record PLINK genotype QC on both funnels; the MAF floor ran only when one was set."""
    for name in (REFERENCE_PANEL, ANNOTATION_SNP_GRID):
        funnel.record(name, "genotype_qc", genotype_qc_removed)
        if maf_min is not None:
            funnel.record(name, "maf_floor", maf_removed)


def _align_annotations(bundle, reference, chrom, global_config, funnel):
    """Align annotation rows once to retained reference identities and recover alleles.

    The annotation rows are the annotation SNP grid. Panel SNPs without an
    annotation row and annotation rows without a panel SNP are both dropped and
    counted in ``funnel`` at the ``annotation_snp_grid`` step.
    """
    if is_allele_aware_mode(global_config.snp_identifier):
        cleaned = clean_identity_artifact_table(
            reference, global_config.snp_identifier,
            context=f"reference panel chromosome {chrom}",
            stage="annotation_reference_identity_cleanup", logger=LOGGER,
        ).cleaned
        funnel.record(REFERENCE_PANEL, "allele_aware_identity_cleanup", len(reference) - len(cleaned))
        reference = cleaned
    mode = _annotation_reference_match_mode(bundle.metadata, global_config.snp_identifier)
    reference_keys = build_snp_id_series(reference, mode)
    annotation_keys = build_snp_id_series(bundle.metadata, mode)
    keep = annotation_keys.isin(set(reference_keys))
    on_grid = int(reference_keys.isin(set(annotation_keys)).sum())
    funnel.record(REFERENCE_PANEL, ANNOTATION_SNP_GRID, len(reference) - on_grid)
    funnel.record(ANNOTATION_SNP_GRID, "grid_rows", len(bundle.metadata))
    funnel.record(ANNOTATION_SNP_GRID, ANNOTATION_SNP_GRID, int((~keep).sum()))
    if not keep.any():
        error = EmptyReferenceSNPs(
            f"ldscore retained no reference SNPs on chromosome {chrom}: none of its {len(bundle.metadata)} "
            f"annotation rows matches the {len(reference)} {funnel.backend} reference-panel SNPs left "
            "after the panel's filters. Most likely the annotation SNP identifiers, genome build, or "
            "allele-aware identifier mode do not match the reference panel. Use annotation and reference-panel "
            "artifacts built with the same SNP identifier mode and genome build. "
            f"Other causes & fixes: {_LDSCORE_INTERSECTION_DOC}"
        )
        error.reference_snp_funnel = funnel.frame()
        raise error
    metadata = _annotation_metadata_with_reference_alleles(
        annotation_metadata=bundle.metadata.loc[keep].reset_index(drop=True),
        annotation_keys=annotation_keys.loc[keep].reset_index(drop=True),
        reference_metadata=reference, reference_keys=reference_keys,
        snp_identifier=global_config.snp_identifier,
    )
    if on_grid != len(metadata):
        raise LDSCInternalError(
            f"Reference preparation matched {on_grid} panel SNPs but {len(metadata)} annotation rows on "
            f"chromosome {chrom}; SNP identities must be unique on both sides of the match. "
            "Report the command, package version, and traceback."
        )
    return ChromosomeAnnotations(metadata,
        MappedAnnotations(bundle.annotations, np.flatnonzero(keep.to_numpy()),
                          tuple(bundle.baseline_columns + bundle.query_columns)),
        list(bundle.baseline_columns), list(bundle.query_columns))


def _annotation_reference_match_mode(
    annotation_metadata: pd.DataFrame,
    snp_identifier: str,
) -> str:
    """Return the SNP identifier mode usable by both annotation and reference metadata."""
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


_LDSCORE_INTERSECTION_DOC = (
    "docs/user-guide/troubleshooting.md#ldscore-no-reference-snps-remain-on-the-annotation-snp-grid"
)


def resolve_keep_individuals(keep_path: str | None, fam) -> list[int] | None:
    """Return FAM row indices retained by ``keep_path`` using IID matching."""
    if not keep_path:
        return None
    keep_indivs = fam.loj(legacy_parse.FilterFile(keep_path).IDList)
    if len(keep_indivs) == 0:
        raise LDSCInputError(
            f"ldscore retained no PLINK individuals after applying keep file '{keep_path}'. "
            "Most likely the keep-list IDs do not match the `.fam` file. Check IID/FID "
            "columns in the keep file and pass a keep list from the same PLINK sample set."
        )
    return keep_indivs.tolist()
