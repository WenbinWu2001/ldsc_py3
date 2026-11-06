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
import logging

import numpy as np
import pandas as pd

from .._core.chromosomes import chrom_sort_key, normalize_chromosome
from .._core.columns import normalize_genome_build
from ..config import GlobalConfig, LDScoreConfig, RefPanelConfig, validate_config_compatibility
from ..errors import EmptyReferenceSNPs, LDSCConfigError, LDSCInputError, LDSCInternalError
from .._genome.build_inference import resolve_genome_build, validate_auto_genome_build_mode
from .._runtime.paths import inspect_plink_inputs, split_cli_path_tokens
from .._runtime.preflight import InputGate
from .._core.identifiers import identifier_keys, validate_retained_identifier_uniqueness
from .funnel import ANNOTATION_SNP_GRID, REFERENCE_PANEL, ReferenceSNPFunnel
from .plink import PlinkBEDFile
from .plink_source import (
    PlinkChromosomeSource,
    build_plink_metadata_frame,
    build_restriction_mask,
    read_keep_individuals,
    select_plink_chromosome,
)
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
    RestrictionIdentityKeys,
    identity_base_mode,
    identity_mode_family,
    is_allele_aware_mode,
)
from .._genome.build_policy import GenomeBuildCheck
from .genetic_map import interpolate_genetic_map_cm, load_genetic_map_group

LOGGER = logging.getLogger("LDSC.ref_panel")


def _artifact_global_config(metadata, global_config, source):
    """Validate declared coordinates before constructing an artifact snapshot."""
    build = metadata.get("genome_build")
    if identity_mode_family(global_config.snp_identifier) == "chr_pos":
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
        self._restriction_keys: RestrictionIdentityKeys | None = None

    def _resolve_run_inputs(self) -> None:
        """Read the run-level selection inputs now, once for this panel.

        Every later chromosome reuses them, and a copy of the panel sent to a
        worker process carries them, so workers do not read them again. The
        base panel reads ``ref_panel_snps_file``.
        """
        self._restriction()

    def _restriction(self) -> RestrictionIdentityKeys | None:
        """Return the parsed ``ref_panel_snps_file`` keys, read on first use; ``None`` without one."""
        if self.spec.ref_panel_snps_file is None:
            return None
        if self._restriction_keys is None:
            self._restriction_keys = read_snp_restriction_keys(
                self.spec.ref_panel_snps_file,
                self.global_config.snp_identifier,
                genome_build=self.global_config.genome_build,
                logger=LOGGER,
            )
        return self._restriction_keys

    @abstractmethod
    def available_chromosomes(self) -> list[str]:
        """Return chromosomes available from the backend."""
        raise NotImplementedError

    @abstractmethod
    def load_metadata(self, chrom: str) -> pd.DataFrame:
        """Load prepared SNP metadata for ``chrom`` after panel-layer filters.

        Raises ``EmptyReferenceSNPs``, carrying the chromosome's
        ``reference_snp_funnel``, when those filters leave no SNP.
        """
        raise NotImplementedError

    @abstractmethod
    def check_reference(self, chromosomes, config: LDScoreConfig) -> None:
        """Check the reference files of ``chromosomes`` for ``config`` before any chromosome is computed.

        Raises ``LDSCInputError`` listing every defect found. Inputs that are
        checked only when a chromosome is prepared are not read here.
        """
        raise NotImplementedError

    @abstractmethod
    def reference_metadata_sources(self, config: LDScoreConfig) -> tuple[str, str]:
        """Return the ``(CM, MAF)`` source labels of the metadata ``prepare_chromosome`` gives under ``config``.

        PLINK panels give ``genetic_map`` or ``plink_bim`` for CM and
        ``plink_genotypes`` for MAF; R2 matrices give ``parquet_sidecar`` for both.
        """
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

    def _apply_keep_restriction(self, metadata: pd.DataFrame) -> pd.DataFrame:
        """Filter metadata rows to ``RefPanelConfig.ref_panel_snps_file`` when set.

        The restriction matches as in PLINK selection: a row whose alleles
        cannot form an allele-aware identity is kept when its base key is
        restricted, so identity cleanup, which runs afterwards, drops and
        reports it.
        """
        restriction = self._restriction()
        if restriction is None or len(metadata) == 0:
            return metadata
        keep = build_restriction_mask(metadata, restriction, self.global_config.snp_identifier)
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
        self._keep_individuals = None

    def _resolve_run_inputs(self) -> None:
        """Read ``ref_panel_snps_file`` and ``keep_indivs_file`` now, once for this panel."""
        super()._resolve_run_inputs()
        self._kept_individuals()

    def _kept_individuals(self):
        """Return the parsed ``keep_indivs_file`` IDs, read on first use; ``None`` without one."""
        if self._keep_individuals is None:
            self._keep_individuals = read_keep_individuals(self.spec.keep_indivs_file)
        return self._keep_individuals

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

    def _load_source(self, chrom: str) -> tuple[PlinkChromosomeSource, ReferenceSNPFunnel]:
        """Select one chromosome's PLINK rows and start its reference-SNP funnel.

        ``select_plink_chromosome`` applies ``ref_panel_snps_file`` first and
        runs identity cleanup on the restricted rows, as for ``build-r2-matrix``.
        When cleanup removes every row the restriction kept,
        ``EmptyReferenceSNPs`` skips the chromosome as the restriction alone
        would; when it removes every row of an unrestricted chromosome,
        ``LDSCInputError`` stops the run.
        """
        source = select_plink_chromosome(
            self._bim_prefixes(chrom)[0], chrom, self.global_config.snp_identifier,
            restriction=self._restriction(), keep_individuals=self._kept_individuals(),
            stage="plink_reference_identity_cleanup",
        )
        funnel = ReferenceSNPFunnel(normalize_chromosome(chrom), "plink")
        funnel.record(REFERENCE_PANEL, "panel_rows", source.panel_rows)
        if source.restriction_removed is not None:
            funnel.record(REFERENCE_PANEL, "ref_panel_snps_file", source.restriction_removed)
        funnel.record(REFERENCE_PANEL, "identity_cleanup", source.identity_cleanup_removed)
        removed = source.restricted_rows
        if removed and source.metadata.empty:
            # Neither outcome writes diagnostics/dropped_snps/, so the message names the drops.
            reasons = ", ".join(f"{reason}={count:,}"
                                for reason, count in source.identity_drops["reason"].value_counts().sort_index().items())
            if source.restriction_removed is not None:
                # Skipped like a chromosome the restriction alone empties.
                error = EmptyReferenceSNPs(
                    f"Reference-panel loading retained no PLINK rows on chromosome {chrom}: SNP identity cleanup "
                    f"removed all {removed:,} rows that --ref-panel-snps-file kept ({reasons})."
                )
                error.reference_snp_funnel = funnel.frame()
                raise error
            raise LDSCInputError(
                f"Reference-panel loading retained no PLINK rows on chromosome {chrom} after SNP identity cleanup "
                f"removed all {removed:,} of them ({reasons}). Under SNP identifier mode "
                f"'{self.global_config.snp_identifier}' each panel row needs a unique SNP identity, and in "
                "allele-aware modes one A/C/G/T allele pair per position. Remove the duplicate and non-SNP rows "
                "from this chromosome's panel (`plink2 --snps-only just-acgt` removes non-SNP rows), or use the "
                f"SNP identifier mode the panel was prepared for. Other causes & fixes: {_PLINK_IDENTITY_DOC}"
            )
        return replace(source, metadata=self._validate_metadata(source.metadata, chrom)), funnel

    def _load_genotypes(self, chrom, source, annotations=None):
        """Open one BED and align its retained rows to annotations without another read.

        The metadata gives each SNP's minor allele in the retained reference
        individuals as ``A1``, as ``build-r2-matrix`` does, so ``MAF`` is the
        frequency of ``A1``; at exactly 0.5 the ``.bim`` order is kept.
        """
        source, funnel = source
        reference = source.metadata
        if annotations is None:
            annotations = ChromosomeAnnotations(
                reference.drop(columns="_raw_index"), pd.DataFrame(index=reference.index), [], [],
            )
        aligned = _align_annotations(annotations, reference, chrom, self.global_config, funnel)
        mode = self.global_config.snp_identifier
        raw_keys = identifier_keys(reference, mode)
        keys = identifier_keys(aligned.metadata, mode)
        physical_rows = pd.Series(reference["_raw_index"].to_numpy(), index=raw_keys).loc[keys].tolist()
        maf_min = self.spec.maf_min
        try:
            geno = PlinkBEDFile(
                source.prefix + ".bed", len(source.fam.IDList), source.bim, keep_snps=physical_rows,
                keep_indivs=source.keep_indivs, mafMin=maf_min,
            )
        except EmptyReferenceSNPs as exc:
            _record_genotype_qc(funnel, exc.genotype_qc_removed, exc.maf_removed, maf_min)
            exc.reference_snp_funnel = funnel.frame()
            raise
        _record_genotype_qc(funnel, geno.genotype_qc_removed, geno.maf_removed, maf_min)
        try:
            metadata = build_plink_metadata_frame(bim=source.bim, kept_snps=geno.kept_snps, maf_values=geno.maf,
                                                  freq_values=geno.freq)
            retained_keys = identifier_keys(metadata, mode)
            indices = pd.Series(np.arange(len(keys)), index=keys).loc[retained_keys].to_numpy(dtype=np.int64)
            values = MappedAnnotations(aligned.annotations, indices, aligned.annotations.columns)
            return geno, metadata, values, source.identity_drops, len(aligned.metadata)
        except BaseException:
            geno.close()
            raise

    def load_metadata(self, chrom: str) -> pd.DataFrame:
        """Load one chromosome's genotype-filtered metadata; the caller owns its lifetime.

        The genotypes are read to compute ``MAF`` and orient ``A1`` to the minor
        allele, so an unreadable BED is an error. A chromosome that the SNP
        restriction, identity cleanup of the rows it kept, or genotype QC
        empties raises ``EmptyReferenceSNPs`` with its funnel.
        """
        chrom = normalize_chromosome(chrom)
        source = self._load_source(chrom)
        if source[0].metadata.empty:
            raise _panel_filter_emptied(chrom, source[1])
        geno, metadata, _, _, _ = self._load_genotypes(chrom, source)
        geno.close()
        return metadata

    def check_reference(self, chromosomes, config):
        """PLINK files are read and checked when each chromosome is prepared."""

    def reference_metadata_sources(self, config):
        """CM comes from the genetic map only for ``--ld-wind-cm``; kb and SNP windows keep the ``.bim`` CM."""
        has_map = bool(self.spec.genetic_map_hg19_sources or self.spec.genetic_map_hg38_sources)
        return ("genetic_map" if has_map and config.ld_wind_cm is not None else "plink_bim"), "plink_genotypes"

    def prepare_chromosome(self, chrom, annotations, config, *, genetic_map=None):
        """Prepare one genotype-filtered chromosome with aligned annotations and LD windows."""
        source = self._load_source(chrom)
        geno, metadata, values, drops, input_rows = self._load_genotypes(chrom, source, annotations)
        funnel = source[1].frame()
        try:
            if config.ld_wind_cm is not None:
                if genetic_map is None:
                    genetic_map = _resolve_genetic_map(self, chrom, metadata)
                if genetic_map is not None:
                    metadata["CM"] = interpolate_genetic_map_cm(
                        normalize_chromosome(chrom), metadata["POS"].to_numpy(dtype=np.int64), genetic_map,
                        command="ldscore",
                    )
                else:
                    assert_cm_usable(metadata["CM"], chrom)
            block_left = _prepare_window(metadata, chrom, config, self.global_config.snp_identifier)
            return PreparedChromosome(
                reader=geno, metadata=metadata,
                annotations=values, block_left=block_left,
                baseline_columns=list(annotations.baseline_columns), query_columns=list(annotations.query_columns),
                reference_rows_before_genotype_qc=input_rows,
                genotype_qc_removed=geno.genotype_qc_removed, maf_removed=geno.maf_removed,
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
        """List the chromosomes whose files are in the R2 matrix directory."""
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

    def open_chromosome(self, chrom: str) -> R2MatrixChromosome:
        """Open one chromosome of the configured R2 matrix for ``ldscore``.

        The R2 Parquet footer must record this panel's SNP identifier mode and
        genome build. The returned handle parses its sidecar once, on first use.
        """
        assert self.global_config.genome_build in {"hg19", "hg38", None}, (
            f"genome_build reached kernel as {self.global_config.genome_build!r}; "
            "should have been resolved at workflow entry."
        )
        build_dir = resolve_r2_build_dir(self.spec.r2_matrix_dir, self.global_config.genome_build)
        chromosome = R2MatrixChromosome.open(build_dir, chrom, command="ldscore")
        _validate_r2_identity_matches_config(
            {"snp_identifier": chromosome.contract.snp_identifier, "genome_build": chromosome.contract.genome_build},
            self.global_config, source=str(chromosome.parquet_path),
        )
        return chromosome

    def check_chromosome(self, chrom: str, *, on_pair_block=None) -> None:
        """Check one chromosome's R2 Parquet footer and its binding to the sidecar.

        Opening the chromosome checks the footer against this panel's SNP
        identifier mode and genome build; opening a pair reader over the whole
        sidecar checks that the R2 Parquet file was written for that sidecar.
        With ``on_pair_block``, every stored pair is also read, and
        ``on_pair_block()`` is called after each pair block.
        """
        chromosome = self.open_chromosome(chrom)
        reader = chromosome.pair_reader(chromosome.sidecar(), self.global_config.snp_identifier,
                                        self.global_config.genome_build)
        try:
            if on_pair_block is not None:
                for _ in reader.iter_all_pairs():
                    on_pair_block()
        finally:
            reader.close()

    def check_reference(self, chromosomes, config):
        """Check every chromosome's files, footer, LD window and sidecar binding; one error lists all issues."""
        from .preflight import inspect_r2_headers

        r2_matrix_dir = self.spec.r2_matrix_dir
        gate = InputGate("R2 matrix validation")
        gate.check("reference", r2_matrix_dir, lambda: inspect_r2_headers(
            r2_matrix_dir, self.global_config.genome_build, command="ldscore", ld_window=config.ld_window))
        for chrom in chromosomes:
            gate.check("reference", r2_matrix_dir, lambda chrom=chrom: self.check_chromosome(chrom), chrom=chrom)
        gate.finish()

    def load_metadata(self, chrom: str) -> pd.DataFrame:
        """Load, restrict, and validate metadata owned by this chromosome consumer."""
        metadata, funnel = self._load_reference(self.open_chromosome(chrom))
        if metadata.empty:
            raise _panel_filter_emptied(chrom, funnel)
        return metadata

    def _load_reference(self, chromosome: R2MatrixChromosome) -> tuple[pd.DataFrame, ReferenceSNPFunnel]:
        """Return the restricted sidecar and the chromosome's reference-SNP funnel so far."""
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
            metadata = self._apply_keep_restriction(metadata)
            funnel.record(REFERENCE_PANEL, "ref_panel_snps_file", before - len(metadata))
        if self.spec.maf_min is not None:
            before = len(metadata)
            metadata = self._apply_maf_filter(metadata, chromosome.chrom)
            funnel.record(REFERENCE_PANEL, "maf_floor", before - len(metadata))
        metadata = self._validate_metadata(metadata, chromosome.chrom)
        return metadata, funnel

    def reference_metadata_sources(self, config):
        """The sidecar is authoritative for CM and MAF under every LD window."""
        return "parquet_sidecar", "parquet_sidecar"

    def prepare_chromosome(self, chrom, annotations, config, *, genetic_map=None):
        """Prepare a sidecar-aligned streaming reader over a contract-validated R2 Parquet file."""
        chromosome = self.open_chromosome(chrom)
        chromosome.check_ld_window(config.ld_window)
        chromosome.reject_truncation()
        reference, funnel = self._load_reference(chromosome)
        aligned = _align_annotations(annotations, reference, chrom, self.global_config, funnel)
        keys = identifier_keys(aligned.metadata, self.global_config.snp_identifier)
        reference_keys = identifier_keys(reference, self.global_config.snp_identifier)
        lookup = reference.loc[:, ["CM", "MAF"]].set_index(reference_keys)
        metadata = aligned.metadata.copy()
        for column in ("CM", "MAF"):
            metadata[column] = lookup.loc[keys, column].to_numpy()
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
            reader=reader, metadata=metadata,
            annotations=mapped_annotations, block_left=block_left,
            baseline_columns=list(annotations.baseline_columns), query_columns=list(annotations.query_columns),
            cm_source="parquet_sidecar", reference_snp_funnel=funnel.frame(),
        )


class RefPanelLoader:
    """Instantiate the backend requested by a :class:`RefPanelConfig`."""
    def __init__(self, global_config: GlobalConfig) -> None:
        """Store the shared configuration used for every loaded panel."""
        self.global_config = global_config

    def load(self, ref_panel_spec: RefPanelConfig, *, chromosome_prefixes: dict[str, str] | None = None) -> RefPanel:
        """Return the adapter that ``ref_panel_spec.backend`` names.

        In a coordinate-family SNP identifier mode whose genome build is
        ``auto``, the build is resolved here from the panel and
        ``ref_panel_snps_file``. A concrete build is taken as already
        resolved, as the workflows resolve it once against all their inputs;
        each R2 matrix chromosome still compares its declared build with it
        when opened. ``chromosome_prefixes`` maps normalized chromosome labels
        to PLINK prefixes the caller has validated (``PlinkRefPanel``); an R2
        matrix ignores it.
        """
        global_config = self._global_config_for_spec(ref_panel_spec)
        if ref_panel_spec.backend == "plink":
            return PlinkRefPanel(global_config, ref_panel_spec, chromosome_prefixes=chromosome_prefixes)
        return R2MatrixRefPanel(global_config, ref_panel_spec)

    def _global_config_for_spec(self, ref_panel_spec: RefPanelConfig) -> GlobalConfig:
        """Resolve an ``auto`` build in a coordinate-family mode before backend loading."""
        mode = normalize_snp_identifier_mode(self.global_config.snp_identifier)
        if identity_mode_family(mode) != "chr_pos" or self.global_config.genome_build in {"hg19", "hg38"}:
            return self.global_config
        resolved_build = _infer_restricted_ref_panel_build(ref_panel_spec, self.global_config.genome_build)
        return GlobalConfig(
            snp_identifier=mode,
            genome_build=resolved_build,
            log_level=self.global_config.log_level,
        )


def _infer_restricted_ref_panel_build(ref_panel_spec: RefPanelConfig, requested="auto") -> str:
    """Resolve all panel shards and coordinate restrictions before filtering."""
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
        metadata["POS"], metadata["CM"] if config.ld_wind_cm is not None else None, config.ld_window,
        snp_positions=snp_positions,
    )
    block_left = get_block_lefts(coords, distance)
    check_whole_chromosome_window(block_left, chrom, allowed=config.whole_chromosome_ok)
    return block_left


def _resolve_genetic_map(panel, chrom, metadata):
    """Resolve the configured PLINK map against the prepared reference coordinates."""
    spec = panel.spec
    if not (spec.genetic_map_hg19_sources or spec.genetic_map_hg38_sources):
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
    check = GenomeBuildCheck(build, context="PLINK genetic map", logger=LOGGER)
    check.paths(f"--genetic-map-{build}-sources", sources)
    check.resolve()
    return load_genetic_map_group(split_cli_path_tokens(sources), command="ldscore")


_PANEL_FILTER_OPTIONS = {"ref_panel_snps_file": "--ref-panel-snps-file", "maf_floor": "--maf-min"}


def _panel_filter_emptied(chrom, funnel):
    """Return the error that skips a chromosome whose SNP restriction or MAF floor kept no row."""
    frame = funnel.frame()
    reference = frame.loc[frame["funnel"].eq(REFERENCE_PANEL)]
    emptied = reference.loc[(reference["removed"].fillna(0).gt(0) & reference["remaining"].eq(0)).fillna(False)]
    step = emptied.iloc[0]
    error = EmptyReferenceSNPs(
        f"Reference-panel loading retained no SNPs on chromosome {chrom}: {_PANEL_FILTER_OPTIONS[step.step]} "
        f"removed its last {int(step.removed):,} panel rows."
    )
    error.reference_snp_funnel = frame
    return error


def _record_genotype_qc(funnel, genotype_qc_removed, maf_removed, maf_min):
    """Record PLINK genotype QC on both funnels; the MAF floor ran only when one was set."""
    for name in (REFERENCE_PANEL, ANNOTATION_SNP_GRID):
        funnel.record(name, "genotype_qc", genotype_qc_removed)
        if maf_min is not None:
            funnel.record(name, "maf_floor", maf_removed)


def _align_annotations(bundle, reference, chrom, global_config, funnel):
    """Align annotation rows once to retained reference identities and take the reference alleles.

    The annotation rows are the annotation SNP grid. Panel SNPs without an
    annotation row and annotation rows without a panel SNP are both dropped and
    counted in ``funnel`` at the ``annotation_snp_grid`` step. Each kept row
    gets the ``A1``/``A2`` of its reference SNP: the R2 matrix sidecar's
    minor-allele orientation, or ``.bim`` alleles that PLINK preparation
    reorients after reading genotypes. Annotation files without alleles match
    by the base mode: identity cleanup of the reference, by ``build-r2-matrix``
    or in ``PlinkRefPanel._load_source``, already dropped multi-allelic
    positions, so base keys are unique in the reference.
    """
    mode = _annotation_reference_match_mode(bundle.metadata, global_config.snp_identifier)
    reference_keys = build_snp_id_series(reference, mode)
    annotation_keys = build_snp_id_series(bundle.metadata, mode)
    keep = annotation_keys.isin(reference_keys)
    on_grid = int(reference_keys.isin(annotation_keys).sum())
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
    rows = np.flatnonzero(keep.to_numpy())
    if on_grid != len(rows):
        raise LDSCInternalError(
            f"Reference preparation matched {on_grid} panel SNPs but {len(rows)} annotation rows on "
            f"chromosome {chrom}; SNP identities must be unique on both sides of the match. "
            "Report the command, package version, and traceback."
        )
    # Outputs carry the panel's alleles in every SNP identifier mode; annotation
    # alleles only take part in the match, which ignores their order.
    alleles = reference.loc[:, ["A1", "A2"]].set_axis(reference_keys.to_numpy()).loc[annotation_keys.take(rows)]
    metadata = bundle.metadata.take(rows).assign(A1=alleles["A1"].to_numpy(), A2=alleles["A2"].to_numpy())
    metadata.index = pd.RangeIndex(len(rows))
    return ChromosomeAnnotations(metadata,
        MappedAnnotations(bundle.annotations, rows, tuple(bundle.baseline_columns + bundle.query_columns)),
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


_LDSCORE_INTERSECTION_DOC = (
    "docs/user-guide/troubleshooting.md#ldscore-no-reference-snps-remain-on-the-annotation-snp-grid"
)
_PLINK_IDENTITY_DOC = "docs/user-guide/troubleshooting.md#ldscore-a-plink-chromosome-has-no-usable-snp-identity"
