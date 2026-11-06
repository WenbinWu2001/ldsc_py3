"""Workflow wrapper for building R2 matrices from PLINK.

Core functionality:
    Resolve public PLINK, optional genetic-map, and optional liftover inputs;
    orchestrate chromosome-wise R2 matrix builds; and emit the R2 matrix
    files consumed by LDSC3 workflows.

Overview
--------
This module is the public entry point for the `build-r2-matrix` workflow. It
keeps path resolution, logging configuration, and cross-file validation in the
workflow layer, then delegates pairwise-R2 generation and parquet serialization to
``ldsc._refpanel.build_kernel``. Before chromosome processing begins, the
workflow precomputes deterministic parquet, metadata sidecar, and dropped-SNP
audit destinations and refuses any existing workflow-owned artifact unless
``overwrite=True`` was configured. Concrete single-chromosome invocations own
only their chromosome's output artifacts so parallel array jobs can share an
output directory without deleting sibling chromosome outputs. Parsed workflow
wrappers write ``build-r2-matrix.log`` for multi-chromosome runs and
chromosome-scoped logs for concrete single-chromosome runs; direct builder calls
return data artifact paths only and do not create a log file by default. With
overwrite enabled, successful runs replace current-run artifacts and remove
stale owned parquet, metadata, audit, or log siblings from previous
incompatible configurations.
R2 matrix liftover keeps the historical source-build plus optional
opposite-build UX, but matching chain files are valid only in ``chr_pos``-family
modes and are rejected in rsID-family modes. Duplicate-position handling is
likewise coordinate-family behavior and uses drop-all for coordinate collisions.
"""

from __future__ import annotations

import argparse
import gc
import logging
import sys

from dataclasses import asdict, dataclass, field, replace as dataclass_replace
from pathlib import Path
from typing import Any, Sequence

import numpy as np
import pandas as pd

from .._runtime.serial_calculation import serial_calculation
from .._runtime.cli_help import CLIHelpFormatter, CHROMOSOME_PATH_HELP, SCALAR_PATH_HELP, add_output_dir_argument, add_log_level_argument
from .._core.snp_identity import (
    DEFAULT_SNP_IDENTIFIER,
    SNP_IDENTIFIER_MODES,
    RestrictionIdentityKeys,
    coerce_identity_drop_frame,
    concat_identity_drop_frames,
    emit_dropped_snps_notice,
    empty_identity_drop_frame,
    identity_mode_family,
    sidecar_identity_sha256,
    normalize_snp_identifier_mode,
)
from ..config import (
    DEFAULT_SNP_BATCH_SIZE, GlobalConfig, R2MatrixBuildConfig, get_global_config, print_global_config_banner,
)
from .._runtime.progress import report_phase, PhaseProgress
from .._runtime.paths import ensure_output_directory, remove_output_artifacts, resolve_file_group, resolve_scalar_path
from .._runtime.workflow_logging import (
    log_inputs,
    log_outputs,
    overwrite_failure_marker,
    workflow_logging,
)
from ..errors import LDSCConfigError, LDSCInputError, LDSCInternalError, LDSCUsageError
from .._core import identifiers as kernel_identifiers
from .._core.chromosomes import chrom_sort_key
from .plink import PlinkBEDFile
from .plink_source import (
    KeepIndividuals,
    build_plink_metadata_frame,
    inspect_plink_inputs,
    orientation_flip_sign,
    read_keep_individuals,
    select_plink_chromosome,
)
from . import build_kernel as kernel_builder
from .genetic_map import GENETIC_MAP_HEADER_CHECKS, interpolate_genetic_map_cm, load_genetic_map_group
from .outputs import BUILD_R2_MATRIX_LOG, R2MatrixBuildDirectoryWriter
from .r2_format import (
    build_runtime_metadata_table,
    r2_matrix_file_paths,
    round_to_sidecar_precision,
    write_r2_parquet,
    write_runtime_metadata_sidecar,
)
from .windows import get_block_lefts, window_coordinates
from .._genome.liftover import (
    LiftOverMappingResult,
    LiftOverTranslator,
    duplicate_coordinate_drop_result,
    liftover_drop_report,
    log_liftover_drop_report,
)
from .._runtime.input_issues import input_issue
from .._runtime.preflight import inspect_declared_inputs, InputGate
from .._genome.build_policy import GenomeBuildCheck


LOGGER = logging.getLogger("LDSC.r2_matrix_builder")
_BUILD_REF_PANEL_EMPTY_DOC = "docs/user-guide/troubleshooting.md#build-r2-matrix-no-reference-panel-artifacts-were-produced"
_BUILD_REF_PANEL_LIFTOVER_DOC = (
    "docs/user-guide/troubleshooting.md#build-r2-matrix-liftover-or-genetic-map-configuration-is-incomplete"
)


@dataclass(frozen=True)
class R2MatrixBuildResult:
    """Summary of files produced by one R2 matrix build run.

    Parameters
    ----------
    matrix_name : str
        Name inferred from ``Path(output_dir).name``.
    chromosomes : list of str
        Chromosomes for which artifacts were written.
    output_paths : dict of str to list of str, optional
        Paths grouped by artifact kind. Build-specific keys such as
        ``r2_hg19``, ``r2_hg38``, ``meta_hg19``, and ``meta_hg38`` are present
        only for emitted builds. Workflow log paths are not included.
    config_snapshot : dict, optional
        Build configuration and resolved input provenance.
    """

    matrix_name: str
    chromosomes: list[str]
    output_paths: dict[str, list[str]] = field(default_factory=dict)
    config_snapshot: dict[str, Any] = field(default_factory=dict)


@dataclass
class _BuildState:
    """Shared resources reused across chromosomes during one build run."""
    genetic_map_hg19: Any | None
    genetic_map_hg38: Any | None
    liftover_chain_paths: dict[tuple[str, str], str | None] = field(default_factory=dict)
    restriction_keys: RestrictionIdentityKeys | None = None
    keep_individuals: KeepIndividuals | None = None
    translator_cache: dict[tuple[str, str], LiftOverTranslator] = field(default_factory=dict)


_AUTOSOMES = frozenset(map(str, range(1, 23)))
_EMPTIED_AUTOSOME_REPAIR = (
    "Check that the SNP restriction file, --snp-identifier mode, --source-genome-build, and liftover chain "
    "match this PLINK input (identity, liftover, and duplicate-position drops are listed in "
    "diagnostics/dropped_snps/); an R2 matrix missing an autosome changes downstream LD scores"
)


@dataclass
class _ChromosomeSelection:
    """One chromosome's SNPs after the BIM-level selection stages.

    ``keep_indivs`` holds the retained ``.fam`` rows (``None`` keeps every
    individual). ``empty_stage`` names the stage that removed every SNP, or
    is ``None`` when SNPs remain for genotype reading.
    """
    bim: Any
    fam: Any
    keep_snps: Any
    keep_indivs: list[int] | None
    sidecar_path: Path
    n_bim_snps: int
    hg19_lookup: dict | None = None
    hg38_lookup: dict | None = None
    empty_stage: str | None = None


def _emitted_genome_builds(config: R2MatrixBuildConfig) -> list[str]:
    """Return source and optional target builds emitted by one build config."""
    source_build = config.source_genome_build
    if source_build not in {"hg19", "hg38"}:
        raise LDSCInternalError(
            "build-r2-matrix output planning failed in _emitted_genome_builds(): "
            f"source_genome_build={source_build!r} is not concrete. Most likely source "
            "build inference was skipped before output path generation. Re-run with DEBUG "
            "logging and report the traceback."
        )
    target_build = "hg38" if source_build == "hg19" else "hg19"
    matching_chain = (
        config.liftover_chain_hg19_to_hg38_file
        if source_build == "hg19"
        else config.liftover_chain_hg38_to_hg19_file
    )
    return [source_build] if matching_chain is None else [source_build, target_build]


class R2MatrixBuilder:
    """Build R2 matrices from PLINK input."""

    def __init__(self, global_config: GlobalConfig | None = None) -> None:
        """Initialize the builder with shared identifier and build defaults."""
        self.global_config = global_config or get_global_config()

    def run(self, config: R2MatrixBuildConfig) -> R2MatrixBuildResult:
        """Build fixed-name R2 matrix files for every chromosome.

        Chromosomes run serially in one worker with 1 BLAS thread.
        The caller's numerical-library settings remain unchanged.

        Parameters
        ----------
        config : R2MatrixBuildConfig
            Public build configuration containing one PLINK prefix token,
            optional genetic maps, optional liftover chain paths, output
            directory, and exactly one LD window. If
            ``config.ref_panel_snps_file`` is set, this builder interprets the restriction using
            ``self.global_config.snp_identifier`` against the resolved source
            PLINK build. ``GlobalConfig.genome_build`` is not consulted by this
            workflow. Existing workflow-owned parquet, metadata, dropped-SNP,
            or log artifacts are refused before chromosome processing unless
            ``config.overwrite`` is true. Successful overwrites remove stale
            owned siblings not produced by the current run. Concrete
            single-chromosome invocations own only that chromosome's package;
            full ``@`` chromosome-suite invocations own the all-chromosome
            output package. Failed input checks write
            ``diagnostics/input_issues.tsv`` or
            ``diagnostics/plink_input_issues.tsv``; only a run that owns the
            whole panel removes them, because parallel single-chromosome runs
            share them. Direct calls through this method do not create a
            workflow log; parsed wrappers add a build-r2-matrix workflow log
            through the shared logging context.

        Returns
        -------
        R2MatrixBuildResult
            Summary with ``matrix_name`` inferred from ``config.output_dir`` and
            build-specific artifact paths directly under ``{build}``. The
            result contains data artifacts only, not the workflow log path.

        Raises
        ------
        ValueError
            If no chromosome artifacts are produced or if multiple PLINK inputs
            resolve to the same chromosome.
        """
        output_dir, chrom_sources, resolved_prefixes = self._inspect_inputs(config)
        return self._run(config, output_dir, chrom_sources, resolved_prefixes, workflow_log_path=None)

    def _inspect_inputs(
        self,
        config: R2MatrixBuildConfig,
    ) -> tuple[Path, list[tuple[str, str]], list[str]]:
        """Check every declared input in the calling process and assign PLINK chromosomes from BIM contents.

        Returns the output directory, the ``(prefix, chromosome)`` sources, and
        the resolved PLINK prefixes. Input issues are written under the output
        directory's ``diagnostics/`` before the error is raised. Passing checks
        keep existing issue tables, which parallel chromosome runs share; only
        a whole-panel run removes them, in ``_run``.
        """
        print_global_config_banner(type(self).__name__, self.global_config)
        output_dir = ensure_output_directory(config.output_dir, label="output directory")
        declared_inputs = inspect_declared_inputs(plink=config.plink_prefix, header_checks=GENETIC_MAP_HEADER_CHECKS,
            files=[(name, getattr(config, name, None)) for name in (
                'genetic_map_hg19_sources', 'genetic_map_hg38_sources')],
            scalar_files=[(name, getattr(config, name, None)) for name in ('keep_indivs_file', 'ref_panel_snps_file')],
            issues_path=output_dir/'diagnostics/input_issues.tsv',
            plink_issues_path=output_dir/'diagnostics/plink_input_issues.tsv', remove_on_pass=False)
        plink = inspect_plink_inputs(config.plink_prefix, discovered_members=declared_inputs.plink_members)
        issues_path = output_dir / "diagnostics" / "plink_input_issues.tsv"
        try:
            plink.require_valid()
        except LDSCInputError:
            issues_path.parent.mkdir(parents=True, exist_ok=True)
            pd.DataFrame(plink.issues).to_csv(issues_path, sep="\t", index=False)
            LOGGER.error("PLINK input validation failed; complete repair audit: %s", issues_path)
            raise
        return output_dir, [(prefix, chrom) for chrom, prefix in plink.chromosome_prefixes.items()], plink.prefixes

    @serial_calculation
    def _run(
        self,
        config: R2MatrixBuildConfig,
        output_dir: Path,
        chrom_sources: list[tuple[str, str]],
        resolved_prefixes: list[str],
        *,
        workflow_log_path: Path | None,
    ) -> R2MatrixBuildResult:
        """Build from inspected inputs in one worker; shared by the direct API and the parsed workflow."""
        writer = R2MatrixBuildDirectoryWriter()
        chromosomes = [chrom for _, chrom in chrom_sources]
        scope = writer.chromosome_scope(config.plink_prefix, chromosomes)
        log_paths = [] if workflow_log_path is None else [writer.log_path(workflow_log_path, scope)]
        workflow_log_path = log_paths[0] if log_paths else None
        if log_paths:
            writer.artifact_family(output_dir, scope=scope).preflight(
                overwrite=config.overwrite, additional_paths=log_paths)
        try:
            config = self._resolve_source_genome_build(config, resolved_prefixes)
        except (LDSCInputError, LDSCConfigError):
            # Output builds are still unknown; the log alone has been preflighted.
            with workflow_logging("build-r2-matrix", workflow_log_path, log_level=self.global_config.log_level):
                raise
        snp_identifier_mode = normalize_snp_identifier_mode(self.global_config.snp_identifier)

        family = writer.artifact_family(
            output_dir, scope=scope, builds=_emitted_genome_builds(config), chromosomes=chromosomes)
        metadata_path = family.paths["metadata"]
        stale_paths = family.preflight(overwrite=config.overwrite, additional_paths=log_paths)
        (output_dir / "diagnostics").mkdir(parents=True, exist_ok=True)
        if scope is None:
            # A chromosome run leaves the shared tables to its siblings; a whole-panel run clears them.
            for name in ("input_issues.tsv", "plink_input_issues.tsv"):
                (output_dir / "diagnostics" / name).unlink(missing_ok=True)

        with workflow_logging("build-r2-matrix", workflow_log_path, log_level=self.global_config.log_level):
            build_state = self._prepare_build_state(config)
            log_inputs(
                plink_prefix=config.plink_prefix,
                output_dir=str(output_dir),
                source_genome_build=config.source_genome_build,
            )
            if identity_mode_family(snp_identifier_mode) == "rsid" and len(_emitted_genome_builds(config)) == 1:
                LOGGER.info(
                    "Coordinate duplicate filtering applies only for chr_pos-family snp_identifier modes; "
                    "keeping duplicate CHR/POS rows in rsid source-only R2 matrix builds."
                )
            chrom_sources = self._check_chromosome_selections(
                chrom_sources, config, build_state, output_dir / "diagnostics" / "input_issues.tsv"
            )
            chrom_records: list[tuple[str, dict[str, str]]] = []
            for prefix, chrom in chrom_sources:
                LOGGER.info('Computing and publishing reference chromosome %s (%d/%d).', chrom, len(chrom_records)+1, len(chrom_sources))
                chrom_paths = self._build_chromosome(prefix, chrom, config, build_state)
                if chrom_paths is None:
                    continue
                chrom_records.append((chrom, chrom_paths))

            if not chrom_records:
                raise LDSCInputError(
                    "build-r2-matrix produced no chromosome artifacts from the supplied PLINK inputs. "
                    "Most likely every chromosome was removed by the SNP restriction, identity cleanup, "
                    "liftover, duplicate-coordinate filtering, MAF filtering, or missing chromosome data. "
                    "Inspect the error/log counts and any written dropped-SNP reports to identify the failing stage. "
                    "Correct mismatched inputs or unintended settings; changing an intended restriction or filter changes the reference panel. "
                    f"Other causes & fixes: {_BUILD_REF_PANEL_EMPTY_DOC}"
                )

            chrom_records.sort(key=lambda item: chrom_sort_key(item[0]))
            output_keys = sorted({key for _, paths in chrom_records for key in paths})
            output_paths = {key: [paths[key] for _, paths in chrom_records] for key in output_keys}
            output_paths["metadata"] = [str(metadata_path)]
            writer.write_metadata(metadata_path, output_dir, output_paths, {
                "snp_identifier": self.global_config.snp_identifier,
                "genome_build": self.global_config.genome_build,
                "source_genome_build": config.source_genome_build,
                "emitted_genome_builds": _emitted_genome_builds(config),
                "chromosomes": [str(chrom) for chrom, _ in chrom_records],
                "plink_prefix": config.plink_prefix,
                "resolved_plink_prefixes": [str(prefix) for prefix in resolved_prefixes],
                "ref_panel_snps_file": None if config.ref_panel_snps_file is None else str(config.ref_panel_snps_file),
            })
            result = R2MatrixBuildResult(
                matrix_name=Path(config.output_dir).name,
                chromosomes=[chrom for chrom, _ in chrom_records],
                output_paths=output_paths,
                config_snapshot={
                    "build_config": asdict(config),
                    "global_config": asdict(self.global_config),
                    "resolved_plink_prefixes": list(resolved_prefixes),
                },
            )
            log_outputs(**{key: ", ".join(paths) for key, paths in output_paths.items()})
            remove_output_artifacts(stale_paths)
            return result

    def _resolve_source_genome_build(
        self,
        config: R2MatrixBuildConfig,
        resolved_prefixes: Sequence[str],
    ) -> R2MatrixBuildConfig:
        """Return a config whose PLINK source genome build is concrete."""
        check = GenomeBuildCheck(config.source_genome_build, context="build-r2-matrix", option="--source-genome-build", logger=LOGGER)
        for prefix in resolved_prefixes:
            check.table("--plink-prefix", str(prefix) + ".bim", bim=True)
        check.paths("--ref-panel-snps-file", config.ref_panel_snps_file,
                    restriction_mode=self.global_config.snp_identifier)
        source_build = check.resolve()
        return dataclass_replace(config, source_genome_build=source_build)

    @report_phase(LOGGER, 'validation', 'maps, restrictions, and liftover configuration')
    def _prepare_build_state(self, config: R2MatrixBuildConfig) -> _BuildState:
        """Resolve optional maps, chain paths, and source-build restrictions.

        A chain matching the resolved source build enables opposite-build
        emission in ``chr_pos``-family modes. The same matching chain is
        rejected in ``rsID``-family modes because row identity is the SNP label
        and cross-build coordinate emission would be ambiguous.
        """
        hg19_files = None
        hg38_files = None
        if config.genetic_map_hg19_sources is not None:
            hg19_files = resolve_file_group(
                (config.genetic_map_hg19_sources,),
                label="hg19 genetic map",
                allow_chromosome_suite=True,
            )
        if config.genetic_map_hg38_sources is not None:
            hg38_files = resolve_file_group(
                (config.genetic_map_hg38_sources,),
                label="hg38 genetic map",
                allow_chromosome_suite=True,
            )
        for build, paths in (("hg19", hg19_files), ("hg38", hg38_files)):
            if paths:
                check = GenomeBuildCheck(build, context=f"build-r2-matrix {build} genetic map", logger=LOGGER)
                for path in paths:
                    check.table(f"--genetic-map-{build}-sources", path)
                check.resolve()
        source_build = config.source_genome_build
        if source_build not in {"hg19", "hg38"}:
            raise LDSCInternalError(
                "build-r2-matrix shared-state preparation failed: source_genome_build is "
                f"{source_build!r}, not a concrete build. Most likely source build inference "
                "was skipped before state preparation. Re-run with DEBUG logging and report the traceback."
            )
        target_build = "hg38" if source_build == "hg19" else "hg19"
        matching_chain = (
            config.liftover_chain_hg19_to_hg38_file
            if source_build == "hg19"
            else config.liftover_chain_hg38_to_hg19_file
        )
        nonmatching_chain = (
            config.liftover_chain_hg38_to_hg19_file
            if source_build == "hg19"
            else config.liftover_chain_hg19_to_hg38_file
        )
        if matching_chain is not None and identity_mode_family(self.global_config.snp_identifier) == "rsid":
            raise LDSCUsageError(
                "build-r2-matrix cannot use chain liftover in rsID-family SNP identifier modes. "
                "Most likely a matching liftover chain was supplied with `--snp-identifier rsid` "
                "or `rsid_allele_aware`. Use a chr_pos-family SNP identifier mode, or omit "
                "the matching liftover chain."
            )
        restriction_keys = None
        restriction_source = config.ref_panel_snps_file
        if restriction_source:
            restriction_path = resolve_scalar_path(
                restriction_source,
                label="reference-panel SNP restriction",
            )
            restriction_mode = normalize_snp_identifier_mode(self.global_config.snp_identifier)
            if identity_mode_family(restriction_mode) == "chr_pos":
                LOGGER.info(
                    f"Interpreting reference-panel SNP restriction '{restriction_path}' "
                    f"against source_genome_build='{source_build}' PLINK coordinates. "
                    "GlobalConfig.genome_build is ignored by build-r2-matrix."
                )
            restriction_keys = _read_ref_panel_snp_restriction(
                restriction_path,
                restriction_mode,
                source_genome_build=source_build,
            )
            message = (
                f"Loaded {len(restriction_keys.keys)} SNP restriction identifiers "
                f"in {restriction_mode} mode from '{restriction_path}'."
            )
            LOGGER.info(message)
        if matching_chain is not None:
            if config.ld_wind_cm is not None and (
                (target_build == "hg38" and config.genetic_map_hg38_sources is None)
                or (target_build == "hg19" and config.genetic_map_hg19_sources is None)
            ):
                raise LDSCUsageError(
                    f"build-r2-matrix cannot emit {target_build} with `--ld-wind-cm` because "
                    f"the {target_build} genetic map was not supplied. Most likely liftover "
                    "enabled an opposite-build output without its matching genetic map. "
                    f"Pass `--genetic-map-{target_build}-sources <path>` or use an SNP/kb window. "
                    f"Other causes & fixes: {_BUILD_REF_PANEL_LIFTOVER_DOC}"
                )
            if (target_build == "hg38" and config.genetic_map_hg38_sources is None) or (
                target_build == "hg19" and config.genetic_map_hg19_sources is None
            ):
                LOGGER.warning(
                    f"No {target_build} genetic map was provided; "
                    f"{target_build} metadata CM values will be written as NA. "
                    "These missing genetic distances cannot support later cM-window LD-score calculation. "
                    f"Supply --genetic-map-{target_build}-sources when building an R2 matrix intended for cM windows."
                )
        elif matching_chain is None:
            if nonmatching_chain is not None:
                message = (
                    f"No usable liftover chain was provided for {source_build} -> {target_build}; "
                    "ignoring the opposite-direction chain and running source-build-only."
                )
                LOGGER.info(message)
            else:
                LOGGER.info(
                    f"No liftover chain was provided for {source_build} -> {target_build}; "
                    "running source-build-only."
                )
        if config.ld_wind_cm is not None and (
            (source_build == "hg38" and config.genetic_map_hg38_sources is None)
            or (source_build == "hg19" and config.genetic_map_hg19_sources is None)
        ):
            raise LDSCUsageError(
                f"build-r2-matrix cannot use `--ld-wind-cm` because the {source_build} "
                "genetic map was not supplied. Most likely the source-build map option "
                "was omitted. Pass the matching `--genetic-map-*-sources` file or use "
                "`--ld-wind-kb` / `--ld-wind-snps`. "
                f"Other causes & fixes: {_BUILD_REF_PANEL_LIFTOVER_DOC}"
            )
        if (source_build == "hg38" and config.genetic_map_hg38_sources is None) or (
            source_build == "hg19" and config.genetic_map_hg19_sources is None
        ):
            LOGGER.warning(
                f"No {source_build} genetic map was provided; "
                f"{source_build} metadata CM values will be written as NA. "
                "These missing genetic distances cannot support later cM-window LD-score calculation. "
                f"Supply --genetic-map-{source_build}-sources when building an R2 matrix intended for cM windows."
            )
        return _BuildState(
            genetic_map_hg19=None if hg19_files is None else load_genetic_map_group(hg19_files, command="build-r2-matrix"),
            genetic_map_hg38=None if hg38_files is None else load_genetic_map_group(hg38_files, command="build-r2-matrix"),
            liftover_chain_paths={
                ("hg19", "hg38"): config.liftover_chain_hg19_to_hg38_file,
                ("hg38", "hg19"): config.liftover_chain_hg38_to_hg19_file,
            },
            restriction_keys=restriction_keys,
            keep_individuals=None if config.keep_indivs_file is None else read_keep_individuals(
                resolve_scalar_path(config.keep_indivs_file, label="individual keep file")),
        )


    def _check_chromosome_selections(self, chrom_sources, config, build_state, issues_path):
        """Run SNP selection on every chromosome before any genotype reading.

        Autosomes left without SNPs are collected into one error, because a
        R2 matrix missing an autosome silently changes downstream LD scores. Other
        chromosomes left without SNPs are skipped with a warning. Returns the
        ``(prefix, chrom)`` pairs that keep SNPs.
        """
        issues, buildable = [], []
        with PhaseProgress(LOGGER, 'validation', 'reference-panel SNP selection', len(chrom_sources)) as progress:
            for prefix, chrom in chrom_sources:
                selection = self._select_chromosome_snps(prefix, chrom, config, build_state)
                progress.advance(object=f'chromosome {chrom}')
                if selection.empty_stage is None:
                    buildable.append((prefix, chrom))
                elif chrom in _AUTOSOMES:
                    details = f"all {selection.n_bim_snps} SNPs were removed at {selection.empty_stage}"
                    issues.append(input_issue('reference panel PLINK input', prefix, chrom, 'no_snps_remain', details,
                                              repair=_EMPTIED_AUTOSOME_REPAIR))
                else:
                    LOGGER.warning(
                        f"Skipping chromosome {chrom}: all {selection.n_bim_snps} SNPs were removed at "
                        f"{selection.empty_stage}; the R2 matrix will not include this chromosome."
                    )
        InputGate('Reference-panel SNP selection', issues_path, issues,
                  troubleshooting=_BUILD_REF_PANEL_EMPTY_DOC, remove_on_pass=False).finish()
        return buildable

    def _select_chromosome_snps(
        self,
        prefix: str,
        chrom: str,
        config: R2MatrixBuildConfig,
        build_state: _BuildState,
    ) -> _ChromosomeSelection:
        """Apply the BIM-level SNP selection stages to one chromosome.

        ``select_plink_chromosome`` applies the SNP restriction and identity
        cleanup, as for ``ldscore --plink-prefix``; liftover plus
        duplicate-position filtering follow on the rows it returns. No stage
        reads genotypes. The dropped-SNP table is always written. Reason
        rows cover source identity cleanup and liftover-stage drops;
        restriction-file drops are outside that audit vocabulary.
        """
        sidecar_path = R2MatrixBuildDirectoryWriter.dropped_snps_path(config.output_dir, chrom)
        source = select_plink_chromosome(
            prefix, chrom, self.global_config.snp_identifier,
            restriction=build_state.restriction_keys, keep_individuals=build_state.keep_individuals,
            stage="ref_panel_source_identity_cleanup",
        )
        selection = _ChromosomeSelection(source.bim, source.fam, np.empty(0, dtype=int), source.keep_indivs,
                                         sidecar_path, source.panel_rows)
        if source.restricted_rows == 0:
            R2MatrixBuildDirectoryWriter.write_dropped_snps(empty_identity_drop_frame(), sidecar_path, chrom)
            return dataclass_replace(selection, empty_stage="SNP restriction")

        dropped_frames: list[pd.DataFrame] = []
        if not source.identity_drops.empty:
            dropped_frames.append(source.identity_drops)
            dropped_examples = (
                source.identity_drops.loc[:, ["SNP", "reason", "base_key", "identity_key", "allele_set"]]
                .head(10)
                .to_dict(orient="records")
            )
            LOGGER.debug(f"R2 matrix identity cleanup drops on chromosome {chrom}: {dropped_examples}")
        # The selected rows by physical BIM row, which liftover and duplicate filtering index.
        rows = source.metadata.set_index("_raw_index").loc[:, ["SNP", "POS"]]
        keep_snps = rows.index.to_numpy(dtype=int)
        if len(keep_snps) == 0:
            R2MatrixBuildDirectoryWriter.write_dropped_snps(concat_identity_drop_frames(dropped_frames), sidecar_path, chrom)
            return dataclass_replace(selection, empty_stage="identity cleanup")

        keep_snps, hg19_lookup, hg38_lookup, liftover_drop_frame = self._resolve_mappable_snp_positions(
            build_state=build_state,
            chrom=chrom,
            source_build=config.source_genome_build,
            rows=rows,
            keep_snps=keep_snps,
            sidecar_path=sidecar_path,
        )
        if not liftover_drop_frame.empty:
            dropped_frames.append(liftover_drop_frame)
        if identity_mode_family(self.global_config.snp_identifier) == "chr_pos":
            before_liftover_collision_count = len(keep_snps)
            keep_snps, liftover_collision_df = _resolve_unique_snp_set(
                chrom=chrom,
                rows=rows,
                keep_snps=keep_snps,
                hg19_lookup=hg19_lookup,
                hg38_lookup=hg38_lookup,
                sidecar_path=sidecar_path,
            )
            if not liftover_collision_df.empty:
                dropped_frames.append(liftover_collision_df)
            synchronized_drop_count = before_liftover_collision_count - len(keep_snps)
            if synchronized_drop_count:
                LOGGER.info(
                    "R2 matrix liftover synchronized cross-build drops: dropped %d SNPs from all emitted builds "
                    "because they collided in at least one emitted build.",
                    synchronized_drop_count,
                )
        R2MatrixBuildDirectoryWriter.write_dropped_snps(concat_identity_drop_frames(dropped_frames), sidecar_path, chrom)
        return dataclass_replace(
            selection, keep_snps=keep_snps, hg19_lookup=hg19_lookup, hg38_lookup=hg38_lookup,
            empty_stage=None if len(keep_snps) else "liftover and duplicate-position filtering",
        )

    @report_phase(LOGGER, 'computation/publication', 'reference chromosome')
    def _build_chromosome(
        self,
        prefix: str,
        chrom: str,
        config: R2MatrixBuildConfig,
        build_state: _BuildState,
    ) -> dict[str, str] | None:
        """Build all parquet and metadata artifacts for one chromosome.

        SNP selection is recomputed here rather than kept from the pre-pass,
        so only one chromosome's BIM-level state is in memory at a time.
        Genotype filtering that removes every SNP stops the run for an
        autosome and skips any other chromosome with a warning.
        """
        LOGGER.info(f"Building R2 matrix artifacts for chromosome {chrom} from '{prefix}'.")
        selection = self._select_chromosome_snps(prefix, chrom, config, build_state)
        bim, fam, keep_snps = selection.bim, selection.fam, selection.keep_snps
        hg19_lookup, hg38_lookup = selection.hg19_lookup, selection.hg38_lookup
        sidecar_path, keep_indivs = selection.sidecar_path, selection.keep_indivs

        output_paths: dict[str, str] = {}
        retained_count: int | None = None
        for build in _emitted_genome_builds(config):
            build_keep_snps = _sort_retained_snps_by_build_position(
                keep_snps,
                genome_build=build,
                hg19_lookup=hg19_lookup,
                hg38_lookup=hg38_lookup,
            )
            # Unrestricted builds keep ~every SNP, so selective read cannot shrink
            # the genotype payload; stream it instead. Restricted builds stay in
            # RAM, where the kept set is already small.
            geno = PlinkBEDFile(
                prefix + ".bed",
                len(fam.IDList),
                bim,
                keep_snps=list(build_keep_snps),
                keep_indivs=keep_indivs,
                mafMin=config.maf_min,
                streaming=build_state.restriction_keys is None,
            )
            metadata = build_plink_metadata_frame(
                bim=bim,
                kept_snps=geno.kept_snps,
                maf_values=geno.maf,
                freq_values=geno.freq,
            )
            if len(metadata) == 0:
                maf_rule = f"MAF below --maf-min {config.maf_min}, " if config.maf_min is not None else ""
                if chrom in _AUTOSOMES:
                    raise LDSCInputError(
                        f"build-r2-matrix stopped at chromosome {chrom}: PLINK genotype filtering removed all "
                        f"{len(build_keep_snps)} SNPs that passed SNP selection ({maf_rule}monomorphic, or "
                        "all genotypes missing among the kept individuals). An R2 matrix without this autosome "
                        "would change downstream LD scores. Most likely --maf-min is too high or "
                        "--keep-indivs-file selects individuals without variation on this chromosome. "
                        "Correct an unintended threshold or individual selection; any revised filter "
                        f"changes the R2 matrix. Other causes & fixes: {_BUILD_REF_PANEL_EMPTY_DOC}"
                    )
                LOGGER.warning(
                    f"Skipping chromosome {chrom}: PLINK genotype filtering removed all "
                    f"{len(build_keep_snps)} SNPs that passed SNP selection."
                )
                return None
            if set(metadata["CHR"]) != {chrom}:
                raise LDSCInternalError(
                    f"build-r2-matrix PLINK filtering for chromosome {chrom} retained rows "
                    "from multiple chromosomes. Most likely the retained SNP index mapping "
                    "desynchronized from the BIM table. Re-run with DEBUG logging and report the traceback."
                )

            # _sort_retained_snps_by_build_position has already rejected a missing lookup.
            build_lookup = hg19_lookup if build == "hg19" else hg38_lookup
            build_positions = np.asarray([build_lookup[int(idx)] for idx in geno.kept_snps], dtype=int)
            _validate_emitted_build_chr_pos_uniqueness(
                metadata=metadata,
                positions=build_positions,
                genome_build=build,
                chrom=chrom,
                snp_identifier=self.global_config.snp_identifier,
            )

            cm_values = (
                interpolate_genetic_map_cm(chrom, build_positions, build_state.genetic_map_hg19, command="build-r2-matrix")
                if build == "hg19" and build_state.genetic_map_hg19 is not None
                else interpolate_genetic_map_cm(chrom, build_positions, build_state.genetic_map_hg38, command="build-r2-matrix")
                if build == "hg38" and build_state.genetic_map_hg38 is not None
                else None
            )
            if cm_values is not None:
                # Window on the CM the sidecar will store, which ldscore reads back.
                cm_values = round_to_sidecar_precision(cm_values)
            block_left = _r2_matrix_block_left(chrom, build, build_positions, cm_values, config)
            r2_path, meta_path = r2_matrix_file_paths(Path(config.output_dir) / build, chrom)

            # Build the sidecar table first so the parquet can record the binding
            # hash (CHR:POS:A1:A2) of the exact sidecar shipped alongside it. The
            # pair indices emitted below reference this table's row order.
            runtime_metadata = build_runtime_metadata_table(
                metadata=metadata,
                positions=build_positions,
                cm_values=cm_values,
            )
            identity_hash = sidecar_identity_sha256(runtime_metadata)

            # Canonical A1=minor orientation: negate standardized columns where the
            # .bim A2 allele is minor (freq < 0.5) so SIGN_R is computed in the same
            # orientation as the swapped meta sidecar. nextSNPs stays orientation-free.
            flip_sign = orientation_flip_sign(geno.freq)
            oriented_snp_getter = kernel_builder.make_oriented_snp_getter(
                lambda b: geno.nextSNPs(b, dtype=np.float32),
                flip_sign,
            )
            write_r2_parquet(
                pair_chunks=kernel_builder.yield_pairwise_r2_rows(
                    block_left=block_left,
                    snp_batch_size=config.snp_batch_size,
                    standardized_snp_getter=oriented_snp_getter,
                    m=geno.m,
                    n=geno.n,
                    r2_min=config.r2_min,
                ),
                path=r2_path,
                genome_build=build,
                snp_identifier=self.global_config.snp_identifier,
                n_samples=geno.n,
                r2_min=config.r2_min,
                n_snps=len(runtime_metadata),
                sidecar_identity_sha256=identity_hash,
                ld_window=config.ld_window,
            )
            write_runtime_metadata_sidecar(
                runtime_metadata,
                meta_path,
                genome_build=build,
                snp_identifier=self.global_config.snp_identifier,
            )
            output_paths[f"r2_{build}"] = str(r2_path)
            output_paths[f"meta_{build}"] = str(meta_path)
            retained_count = len(metadata)
            # Release the genotype bitarray before loading the next build so the
            # two builds' decoded genotypes never coexist at peak memory.
            del geno
            gc.collect()
        output_paths["dropped_snps"] = str(sidecar_path)
        LOGGER.info(f"Finished chromosome {chrom} with {retained_count} retained SNPs.")
        return output_paths

    def _resolve_mappable_snp_positions(
        self,
        *,
        build_state: _BuildState,
        chrom: str,
        source_build: str,
        rows: pd.DataFrame,
        keep_snps,
        sidecar_path: Path | None = None,
    ) -> tuple[Any, dict[int, int], dict[int, int], pd.DataFrame]:
        """Resolve retained SNP positions and liftover-stage dropped-SNP rows.

        ``rows`` holds ``SNP`` and source ``POS`` by physical BIM row. Returns
        retained PLINK row indices, source/target build position lookups, and
        a unified nullable sidecar frame for unmapped or cross-chromosome
        chain liftover drops. Source-only builds return an empty sidecar frame
        with the same schema.
        """
        keep_snps = list(keep_snps)
        candidate_positions = rows.loc[keep_snps, "POS"].to_numpy(dtype=int)
        if source_build == "hg19":
            if build_state.liftover_chain_paths.get(("hg19", "hg38")) is None:
                retained_snps = np.asarray(keep_snps, dtype=int)
                hg19_lookup = {
                    int(idx): int(pos)
                    for idx, pos in zip(retained_snps, candidate_positions)
                }
                return retained_snps, hg19_lookup, {}, empty_identity_drop_frame()
            mapping = self._map_positions(build_state, chrom, candidate_positions, "hg19", "hg38")
            retained_snps = np.asarray(keep_snps, dtype=int)[mapping.keep_mask]
            hg19_lookup = {
                int(idx): int(pos)
                for idx, pos in zip(retained_snps, candidate_positions[mapping.keep_mask])
            }
            hg38_lookup = {
                int(idx): int(pos)
                for idx, pos in zip(retained_snps, mapping.translated_positions)
            }
        else:
            if build_state.liftover_chain_paths.get(("hg38", "hg19")) is None:
                retained_snps = np.asarray(keep_snps, dtype=int)
                hg38_lookup = {
                    int(idx): int(pos)
                    for idx, pos in zip(retained_snps, candidate_positions)
                }
                return retained_snps, {}, hg38_lookup, empty_identity_drop_frame()
            mapping = self._map_positions(build_state, chrom, candidate_positions, "hg38", "hg19")
            retained_snps = np.asarray(keep_snps, dtype=int)[mapping.keep_mask]
            hg38_lookup = {
                int(idx): int(pos)
                for idx, pos in zip(retained_snps, candidate_positions[mapping.keep_mask])
            }
            hg19_lookup = {
                int(idx): int(pos)
                for idx, pos in zip(retained_snps, mapping.translated_positions)
            }

        dropped = len(keep_snps) - len(retained_snps)
        liftover_drop_frame = empty_identity_drop_frame()
        if dropped:
            drop_frame = pd.DataFrame(
                {
                    "CHR": [chrom] * len(keep_snps),
                    "SNP": rows.loc[keep_snps, "SNP"].astype(str).tolist(),
                    "source_POS": candidate_positions,
                },
                index=keep_snps,
            )
            unmapped_mask, cross_chrom_mask = mapping.unmapped_mask, mapping.cross_chrom_mask
            liftover_drop_frame = concat_identity_drop_frames(
                [
                    coerce_identity_drop_frame(
                        pd.DataFrame(
                            {
                                "CHR": [chrom] * int(unmapped_mask.sum()),
                                "SNP": drop_frame.loc[unmapped_mask, "SNP"].astype(str).tolist(),
                                "source_pos": drop_frame.loc[unmapped_mask, "source_POS"].tolist(),
                                "target_pos": [pd.NA] * int(unmapped_mask.sum()),
                                "reason": ["unmapped_liftover"] * int(unmapped_mask.sum()),
                            }
                        ),
                        default_stage="ref_panel_liftover",
                    ),
                    coerce_identity_drop_frame(
                        pd.DataFrame(
                            {
                                "CHR": [chrom] * int(cross_chrom_mask.sum()),
                                "SNP": drop_frame.loc[cross_chrom_mask, "SNP"].astype(str).tolist(),
                                "source_pos": drop_frame.loc[cross_chrom_mask, "source_POS"].tolist(),
                                "target_pos": [pd.NA] * int(cross_chrom_mask.sum()),
                                "reason": ["cross_chromosome_liftover"] * int(cross_chrom_mask.sum()),
                            }
                        ),
                        default_stage="ref_panel_liftover",
                    ),
                ]
            )
            log_liftover_drop_report(
                LOGGER,
                liftover_drop_report(
                    drop_frame,
                    unmapped_mask,
                    reason="unmapped_liftover",
                    source_pos_col="source_POS",
                ),
                workflow_label="R2 matrix liftover",
                sidecar_path=sidecar_path,
            )
            log_liftover_drop_report(
                LOGGER,
                liftover_drop_report(
                    drop_frame,
                    cross_chrom_mask,
                    reason="cross_chromosome_liftover",
                    source_pos_col="source_POS",
                ),
                workflow_label="R2 matrix liftover",
                sidecar_path=sidecar_path,
            )
            message = (
                f"Dropping {dropped} SNPs on chromosome {chrom} after liftover filtering "
                f"({mapping.unmapped_count} unmapped, {mapping.cross_chrom_count} cross-chromosome)."
            )
            LOGGER.info(message)
        return retained_snps, hg19_lookup, hg38_lookup, liftover_drop_frame

    def _map_positions(
        self,
        build_state: _BuildState,
        chrom: str,
        positions,
        source_build: str,
        target_build: str,
    ) -> LiftOverMappingResult:
        """Translate positions between builds, caching translators by build pair."""
        key = (source_build, target_build)
        translator = build_state.translator_cache.get(key)
        if translator is None:
            translator = LiftOverTranslator(
                source_build=source_build,
                target_build=target_build,
                chain_path=build_state.liftover_chain_paths.get(key),
                chain_flag_hint=f"--liftover-chain-{source_build}-to-{target_build}-file",
                workflow_label="R2 matrix liftover",
            )
            build_state.translator_cache[key] = translator
        return translator.map_positions(chrom, positions)


def _read_ref_panel_snp_restriction(
    path: str | Path,
    restriction_mode: str,
    *,
    source_genome_build: str,
) -> RestrictionIdentityKeys:
    """Load restriction keys using the resolved, declaration-consistent source build.

    Auto validation has already examined this input alongside every BIM shard.
    rsID filtering ignores coordinates; coordinate filtering prefers a matching
    build-specific position column and otherwise uses generic POS.
    """
    return kernel_identifiers.read_snp_restriction_keys(
        path, restriction_mode, genome_build=source_genome_build, logger=LOGGER,
    )


def _position_lookup_for_build(
    genome_build: str,
    *,
    hg19_lookup: dict[int, int],
    hg38_lookup: dict[int, int],
) -> dict[int, int]:
    """Return the retained-position lookup for one genome build."""
    return hg19_lookup if genome_build == "hg19" else hg38_lookup


def _r2_matrix_block_left(
    chrom: str,
    build: str,
    positions: np.ndarray,
    cm_values: np.ndarray | None,
    config: R2MatrixBuildConfig,
) -> np.ndarray:
    """Return one emitted build's LD-window block-left array."""
    if config.ld_wind_cm is not None:
        if cm_values is None:
            raise LDSCUsageError(
                f"build-r2-matrix cannot compute cM LD windows for chromosome {chrom} "
                f"in {build}: no genetic map values are available. Most likely the "
                f"{build} genetic map is missing this chromosome or was not supplied. "
                "Pass the matching genetic map or use an SNP/kb LD window. "
                f"Other causes & fixes: {_BUILD_REF_PANEL_LIFTOVER_DOC}"
            )
        if np.isnan(np.asarray(cm_values, dtype=float)).any():
            raise LDSCInputError(
                "build-r2-matrix cannot use cM LD windows because interpolated CM values "
                "are missing for retained SNPs. Most likely the genetic map does not cover "
                "all retained positions. Provide a complete genetic map or use an SNP/kb "
                f"LD window. Other causes & fixes: {_BUILD_REF_PANEL_LIFTOVER_DOC}"
            )
    coords, max_dist = window_coordinates(positions, cm_values, config.ld_window)
    return get_block_lefts(coords, max_dist)


def _metadata_with_build_positions(metadata: pd.DataFrame, positions: np.ndarray) -> pd.DataFrame:
    """Return metadata whose ``POS`` column is aligned to one emitted build."""
    build_metadata = metadata.copy()
    build_metadata["POS"] = np.asarray(positions, dtype=np.int64)
    return build_metadata


def _validate_emitted_build_chr_pos_uniqueness(
    *,
    metadata: pd.DataFrame,
    positions: np.ndarray,
    genome_build: str,
    chrom: str,
    snp_identifier: str,
) -> None:
    """Reject emitted-build coordinate collisions when building for chr_pos matching."""
    if identity_mode_family(snp_identifier) != "chr_pos":
        return
    build_metadata = _metadata_with_build_positions(metadata, positions)
    kernel_identifiers.validate_unique_snp_ids(
        build_metadata,
        "chr_pos",
        context=f"R2 matrix build output {genome_build}[{chrom}]",
    )


def _sort_retained_snps_by_build_position(
    keep_snps,
    *,
    genome_build: str,
    hg19_lookup: dict[int, int],
    hg38_lookup: dict[int, int],
) -> np.ndarray:
    """Return retained PLINK indices in one build's genomic position order."""
    # Monotone position order here makes index i -> position monotone, which lets
    # yield_pairwise_r2_rows flush cross-batch pairs with non-decreasing IDX_1.
    keep_snps = np.asarray(keep_snps, dtype=int)
    lookup = _position_lookup_for_build(
        genome_build,
        hg19_lookup=hg19_lookup,
        hg38_lookup=hg38_lookup,
    )
    if not lookup:
        raise LDSCInternalError(
            f"build-r2-matrix retained SNP sorting failed for {genome_build}: no position "
            "lookup is available. Most likely opposite-build emission reached sorting "
            "without a liftover lookup. Re-run with DEBUG logging and report the traceback."
        )
    return np.asarray(
        sorted(keep_snps.tolist(), key=lambda idx: (int(lookup[int(idx)]), int(idx))),
        dtype=int,
    )


def _resolve_unique_snp_set(
    chrom: str,
    rows: pd.DataFrame,
    keep_snps: np.ndarray,
    hg19_lookup: dict[int, int],
    hg38_lookup: dict[int, int],
    sidecar_path: Path | None = None,
) -> tuple[np.ndarray, pd.DataFrame]:
    """Drop SNPs whose emitted positions collide in a chr_pos R2 matrix build.

    ``rows`` holds ``SNP`` and source ``POS`` by physical BIM row. Source
    positions are already unique: identity cleanup drops every
    same-position group in chr_pos-family modes (``duplicate_identity`` or
    ``multi_allelic_base_key``). Target collisions are checked after the
    liftover position lookups are available, resolved with drop-all in every
    emitted build, and returned in the shared dropped-SNP sidecar schema.
    """
    drop_set: set[int] = set()
    prov_rows: list[dict[str, object]] = []

    retained = [int(idx) for idx in keep_snps.tolist()]
    for lookup in (hg19_lookup, hg38_lookup):
        if not lookup or len(retained) < 2:
            continue
        target_pos = [lookup[idx] for idx in retained]
        target_frame = pd.DataFrame(
            {
                "CHR": [chrom] * len(retained),
                "SNP": rows.loc[retained, "SNP"].astype(str).tolist(),
                "source_pos": [int(rows.loc[idx, "POS"]) for idx in retained],
                "target_pos": target_pos,
            },
            index=retained,
        )
        target_dup_result = duplicate_coordinate_drop_result(
            target_frame,
            pos_col="target_pos",
            source_pos_col="source_pos",
            target_pos_col="target_pos",
            reason="target_collision",
        )
        target_dup_mask = ~target_dup_result.keep_mask
        if not target_dup_mask.any():
            continue
        log_liftover_drop_report(
            LOGGER,
            target_dup_result.report,
            workflow_label="R2 matrix liftover",
            sidecar_path=sidecar_path,
        )
        for i, idx in enumerate(retained):
            if target_dup_mask[i] and idx not in drop_set:
                drop_set.add(int(idx))
                prov_rows.append(
                    {
                        "CHR": chrom,
                        "SNP": rows.loc[idx, "SNP"],
                        "source_pos": int(rows.loc[idx, "POS"]),
                        "target_pos": int(target_pos[i]),
                        "reason": "target_collision",
                    }
                )

    cleaned = np.asarray([idx for idx in keep_snps if idx not in drop_set], dtype=keep_snps.dtype)
    columns = ["CHR", "SNP", "source_pos", "target_pos", "reason"]
    if prov_rows:
        dropped_df = coerce_identity_drop_frame(pd.DataFrame(prov_rows, columns=columns), default_stage="ref_panel_liftover")
        dropped_df["target_pos"] = dropped_df["target_pos"].astype("Int64")
    else:
        dropped_df = empty_identity_drop_frame()
    return cleaned, dropped_df


def build_parser() -> argparse.ArgumentParser:
    """Build the feature parser for R2 matrix generation."""
    parser = argparse.ArgumentParser(allow_abbrev=False)
    parser.prog = 'ldsc build-r2-matrix'
    parser.formatter_class = CLIHelpFormatter
    parser.description = 'Build an R2 matrix of pairwise R2 values and SNP metadata from PLINK genotypes.'
    identity = parser.add_argument_group('Step 1: SNP identity and source genome build')
    panel = parser.add_argument_group('Step 2: PLINK reference panel')
    window = parser.add_argument_group('Step 3: LD window (choose one)')
    maps = parser.add_argument_group('Step 4: Genetic maps and liftover (optional)')
    output = parser.add_argument_group('Step 5: Output')
    storage = parser.add_argument_group('Optional: R2 storage threshold')
    computing = parser.add_argument_group('Computing resources')
    logs = parser.add_argument_group('Logging')

    identity.add_argument(
        '--snp-identifier', default=DEFAULT_SNP_IDENTIFIER, choices=SNP_IDENTIFIER_MODES,
        help=(
            'Match SNPs by rsID (rsid) or chromosome and position (chr_pos); allele-aware variants also use '
            'A1/A2. Default: chr_pos_allele_aware; usable alleles are required. Controls duplicate removal, '
            'SNP-list matching, and the written R2 matrix identity.'
        ),
    )
    identity.add_argument(
        '--source-genome-build', default='auto', choices=('auto', 'hg19', 'hg37', 'GRCh37', 'hg38', 'GRCh38'),
        help=(
            'Genome build of input PLINK coordinates. Default: auto, check every BIM shard and coordinate '
            'keep list before filtering. Concrete builds must agree with input declarations; cM maps follow the resolved build.'
        ),
    )

    panel.add_argument(
        '--plink-prefix', required=True, metavar='PREFIX',
        help=(
            'Required PLINK .bed/.bim/.fam filename prefix; chromosome suite stems are also accepted. '
            + CHROMOSOME_PATH_HELP
        ),
    )
    panel.add_argument(
        '--keep-indivs-file', default=None, metavar='FILE',
        help=(
            'File listing one individual ID (IID) per row to retain from PLINK genotypes; an IID listed more '
            'than once counts once. If omitted, retain all individuals. '
            + SCALAR_PATH_HELP
        ),
    )
    panel.add_argument(
        '--ref-panel-snps-file', default=None, metavar='FILE',
        help=(
            'Headered SNP identity table restricting reference SNPs. Coordinate identities must use the PLINK '
            'source build; duplicate keys collapse and extra columns are ignored. If omitted, use all otherwise '
            'eligible SNPs. '
            + SCALAR_PATH_HELP
        ),
    )
    panel.add_argument(
        '--maf-min', default=None, type=float, metavar='VALUE',
        help=(
            'Keep reference SNPs with MAF at or above this value, in [0, 0.5]. If omitted, apply no '
            'additional MAF threshold.'
        ),
    )

    window_choice = window.add_mutually_exclusive_group(required=True)
    window_choice.add_argument(
        '--ld-wind-cm', default=None, type=float, metavar='CM',
        help=(
            'Maximum pair distance on either side of each SNP, in centiMorgans. Requires '
            '--genetic-map-hg19-sources or --genetic-map-hg38-sources for each emitted build. No default; '
            'specify exactly one of --ld-wind-snps, --ld-wind-kb, or --ld-wind-cm.'
        ),
    )
    window_choice.add_argument(
        '--ld-wind-kb', default=None, type=float, metavar='KB',
        help=(
            'Maximum pair distance on either side of each SNP, in kilobases. No default; specify exactly '
            'one of --ld-wind-snps, --ld-wind-kb, or --ld-wind-cm.'
        ),
    )
    window_choice.add_argument(
        '--ld-wind-snps', default=None, type=int, metavar='N',
        help=(
            'Maximum pair distance on either side of each SNP, counted in retained reference-panel SNP '
            'positions. No default; specify exactly one of --ld-wind-snps, --ld-wind-kb, or --ld-wind-cm.'
        ),
    )

    maps.add_argument(
        '--genetic-map-hg19-sources', default=None, metavar='SOURCES',
        help=(
            "hg19 genetic-map files for deriving SNP centiMorgan coordinates. Accepts exact paths, '*' and '@' "
            'patterns; see --plink-prefix for pattern rules. Required with --ld-wind-cm when hg19 output is '
            'written; otherwise optional with no map supplied by default.'
        ),
    )
    maps.add_argument(
        '--genetic-map-hg38-sources', default=None, metavar='SOURCES',
        help=(
            "hg38 genetic-map files for deriving SNP centiMorgan coordinates. Accepts exact paths, '*' and '@' "
            'patterns; see --plink-prefix for pattern rules. Required with --ld-wind-cm when hg38 output is '
            'written; otherwise optional with no map supplied by default.'
        ),
    )
    maps.add_argument(
        '--liftover-chain-hg19-to-hg38-file', default=None, metavar='FILE',
        help=(
            'Chain file enabling additional hg38 output from hg19 coordinates. Requires coordinate-based '
            '--snp-identifier; invalid for rsID-family modes. If omitted, do not convert hg19 input to hg38. Supply an '
            'exact file path; patterns are not expanded.'
        ),
    )
    maps.add_argument(
        '--liftover-chain-hg38-to-hg19-file', default=None, metavar='FILE',
        help=(
            'Chain file enabling additional hg19 output from hg38 coordinates. Requires coordinate-based '
            '--snp-identifier; invalid for rsID-family modes. If omitted, do not convert hg38 input to hg19. Supply an '
            'exact file path; patterns are not expanded.'
        ),
    )

    add_output_dir_argument(
        output,
        help=(
            'Required destination for R2 Parquet files, chromosome SNP metadata, and diagnostics.'
        ),
    )
    output.add_argument(
        '--overwrite', action='store_true', default=False,
        help=(
            "Replace this command's existing R2 matrix files and remove obsolete outputs after success. A "
            'concrete chromosome prefix replaces only that chromosome; a suite replaces the full R2 matrix. '
            'Default: off; stop if output files exist.'
        ),
    )

    storage.add_argument(
        '--r2-min', dest='r2_min', default=0.0, type=float, metavar='VALUE',
        help=(
            'Minimum unbiased R2 to store. Default: 0; zero or negative values keep every pair within the '
            'window. A positive value discards smaller R2 values; ldscore then rejects the R2 matrix, and '
            'query-r2-matrix reports discarded pairs as NaN with status=absent.'
        ),
    )

    computing.add_argument(
        '--snp-batch-size', dest='snp_batch_size', default=DEFAULT_SNP_BATCH_SIZE, type=int, metavar='N',
        help=(
            f'Number of reference SNPs processed per genotype batch. Default: {DEFAULT_SNP_BATCH_SIZE}; larger '
            'batches may improve speed but use more memory. Stored pairs do not depend on it.'
        ),
    )

    add_log_level_argument(logs)
    return parser


def config_from_args(args: argparse.Namespace) -> R2MatrixBuildConfig:
    """Return the build configuration of parsed ``build-r2-matrix`` arguments.

    The PLINK source build remains local to this workflow and is inferred from
    ``.bim`` coordinates during the build when set to ``"auto"``.
    """
    return R2MatrixBuildConfig(
        plink_prefix=args.plink_prefix,
        source_genome_build=args.source_genome_build,
        genetic_map_hg19_sources=args.genetic_map_hg19_sources,
        genetic_map_hg38_sources=args.genetic_map_hg38_sources,
        liftover_chain_hg19_to_hg38_file=args.liftover_chain_hg19_to_hg38_file,
        liftover_chain_hg38_to_hg19_file=args.liftover_chain_hg38_to_hg19_file,
        output_dir=args.output_dir,
        overwrite=args.overwrite,
        ld_wind_snps=args.ld_wind_snps,
        ld_wind_kb=args.ld_wind_kb,
        ld_wind_cm=args.ld_wind_cm,
        maf_min=args.maf_min,
        ref_panel_snps_file=args.ref_panel_snps_file,
        keep_indivs_file=args.keep_indivs_file,
        snp_batch_size=args.snp_batch_size,
        r2_min=args.r2_min,
    )


def _global_config_from_args(args: argparse.Namespace) -> GlobalConfig:
    """Return the run's ``GlobalConfig`` from the parsed SNP identifier mode and log level.

    The SNP identifier mode interprets ``--ref-panel-snps-file`` and is the
    identity the R2 matrix records. ``build-r2-matrix`` reads no genome build
    from ``GlobalConfig``; the PLINK build is ``--source-genome-build``, and
    the recorded ``genome_build`` is ``"auto"`` in coordinate-family modes and
    ``None`` in rsID-family modes. On the Python route the namespace carries
    the registered values (see ``run_build_r2_matrix``).
    """
    snp_identifier = normalize_snp_identifier_mode(args.snp_identifier)
    return GlobalConfig(
        snp_identifier=snp_identifier,
        genome_build="auto" if identity_mode_family(snp_identifier) == "chr_pos" else None,
        log_level=args.log_level,
    )


def run_build_r2_matrix_from_args(args: argparse.Namespace) -> R2MatrixBuildResult:
    """Run R2 matrix building from parsed CLI arguments.

    The parsed workflow preflights all deterministic R2 matrix artifacts plus the
    selected build-r2-matrix log before chromosome processing. The returned
    ``R2MatrixBuildResult`` contains R2 matrix artifact paths only.

    A failed authorized overwrite leaves ``RUN_FAILED.chr<chrom>.txt`` when the
    PLINK prefix has no ``@`` and its BIM files hold one chromosome, and
    ``RUN_FAILED.txt`` otherwise. Errors raised before the chromosomes are
    known, in the options or the input checks, leave no marker: no scientific
    output has been replaced, and a marker name chosen then could not be the
    one a successful retry removes.
    """
    return _build(config_from_args(args), _global_config_from_args(args), command="run_build_r2_matrix_from_args(...)")


def _build(
    build_config: R2MatrixBuildConfig, global_config: GlobalConfig, *, command: str,
) -> R2MatrixBuildResult:
    """Run one build from its configs; ``command`` names the entry point in a failure marker."""
    builder = R2MatrixBuilder(global_config=global_config)
    output_dir, chrom_sources, resolved_prefixes = builder._inspect_inputs(build_config)
    with overwrite_failure_marker(
        output_dir,
        overwrite=build_config.overwrite,
        command=f"Python API: {command}",
        marker_name=R2MatrixBuildDirectoryWriter.failure_marker(R2MatrixBuildDirectoryWriter.chromosome_scope(
            build_config.plink_prefix, [chrom for _, chrom in chrom_sources])),
    ):
        return builder._run(
            build_config, output_dir, chrom_sources, resolved_prefixes,
            workflow_log_path=output_dir / BUILD_R2_MATRIX_LOG,
        )


def run_build_r2_matrix(**kwargs: Any) -> R2MatrixBuildResult:
    """Run the R2 matrix builder with CLI-equivalent keyword arguments.

    The wrapper accepts the same modern I/O and build names as
    ``ldsc build-r2-matrix``. The SNP identifier mode and log level come from
    the registered ``GlobalConfig``; call ``set_global_config(...)`` before
    invoking this wrapper. ``GlobalConfig.genome_build`` is ignored here; pass
    ``source_genome_build`` for the PLINK input build. Both ``plink_prefix``
    and ``output_dir`` are required because this public workflow materializes
    an R2 matrix artifact family.

    Raises
    ------
    LDSCUsageError
        Before anything is written, when a keyword is not a ``build-r2-matrix``
        option name, a required one is missing, a value is outside the
        command's choices, a flag is not a bool, or ``snp_identifier``,
        ``genome_build`` or ``log_level`` is passed.
    """
    from .._runtime.workflow_kwargs import GLOBAL_CONFIG_RUNTIME_OPTIONS, namespace_from_kwargs

    rejected = {
        **GLOBAL_CONFIG_RUNTIME_OPTIONS,
        "genome_build": "pass source_genome_build for the PLINK input build; GlobalConfig.genome_build is not used",
    }
    args = namespace_from_kwargs(
        "run_build_r2_matrix", build_parser(), kwargs, rejected=rejected, registered=get_global_config(),
    )
    return _build(config_from_args(args), _global_config_from_args(args), command="run_build_r2_matrix(...)")


def main(argv: Sequence[str] | None = None) -> R2MatrixBuildResult:
    """Command-line entry point for the R2 matrix builder."""
    from .._runtime.workflow_kwargs import parse_command_line

    args = parse_command_line(build_parser(), list(sys.argv[1:] if argv is None else argv))
    result = run_build_r2_matrix_from_args(args)
    emit_dropped_snps_notice(result.output_paths.get("dropped_snps", []), command="build-r2-matrix")
    return result
