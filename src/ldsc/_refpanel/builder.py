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
import os
import re

from dataclasses import asdict, dataclass, field, replace as dataclass_replace
from pathlib import Path
from typing import Any, Sequence

import numpy as np
import pandas as pd

from .._runtime.artifacts import atomic_write_json
from .._runtime.serial_calculation import serial_calculation
from .._runtime.cli_help import CLIHelpFormatter, CHROMOSOME_PATH_HELP, SCALAR_PATH_HELP
from .._runtime.workflow_logging import FAILURE_MARKER, LOG_LEVEL_HELP
from .._core.columns import (
    normalize_genome_build,
)
from .._core.snp_identity import (
    IDENTITY_DROP_COLUMNS,
    RestrictionIdentityKeys,
    clean_identity_artifact_table,
    coerce_identity_drop_frame,
    empty_identity_drop_frame,
    identity_mode_family,
    sidecar_identity_sha256,
    normalize_snp_identifier_mode,
)
from ..config import (
    DEFAULT_SNP_BATCH_SIZE, GlobalConfig, R2MatrixBuildConfig, get_global_config, print_global_config_banner,
)
from .._core.coordinates import normalize_chr_pos_frame
from .._runtime.progress import report_phase
from .._runtime.paths import ensure_output_directory, preflight_output_artifact_family, remove_output_artifacts, resolve_file_group, inspect_plink_inputs, resolve_scalar_path
from .._runtime.workflow_logging import (
    log_inputs,
    log_outputs,
    materializing_overwrite_guard,
    workflow_logging,
)
from ..errors import LDSCConfigError, LDSCInputError, LDSCInternalError, LDSCUsageError
from .._core import plink_files as legacy_parse
from .._core import identifiers as kernel_identifiers
from .._core.chromosomes import chrom_sort_key
from .plink import PlinkBEDFile
from . import build_kernel as kernel_builder
from .r2_format import (
    R2_MATRIX_FILE_GLOBS,
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
    mapping_reason_masks,
)


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
    restriction_mode: str | None = None
    restriction_values: set[str] | None = None
    restriction_keys: RestrictionIdentityKeys | None = None
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

    ``empty_stage`` names the stage that removed every SNP, or is ``None``
    when SNPs remain for genotype reading.
    """
    bim: Any
    fam: Any
    chrom_df: pd.DataFrame
    keep_snps: Any
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


def _resolve_build_ref_panel_log_path(
    workflow_log_path: Path | None,
    config: R2MatrixBuildConfig,
    chrom_sources: Sequence[tuple[str, str]],
) -> Path | None:
    """Return a collision-resistant workflow log path for chromosome-scoped runs."""
    if workflow_log_path is None:
        return None
    if len(chrom_sources) != 1:
        return workflow_log_path
    if "@" in str(config.plink_prefix):
        return workflow_log_path
    chrom = chrom_sources[0][1]
    return workflow_log_path.with_name(f"{workflow_log_path.stem}.chr{chrom}{workflow_log_path.suffix}")


def _is_chromosome_scoped_invocation(
    config: R2MatrixBuildConfig,
    chrom_sources: Sequence[tuple[str, str]],
) -> bool:
    """Return whether this run is scoped to one concrete chromosome input."""
    return len(chrom_sources) == 1 and "@" not in str(config.plink_prefix)


def _resolve_build_metadata_path(
    output_dir: Path,
    config: R2MatrixBuildConfig,
    chrom_sources: Sequence[tuple[str, str]],
) -> Path:
    """Return a metadata path that does not collide across parallel chromosome runs."""
    if _is_chromosome_scoped_invocation(config, chrom_sources):
        chrom = chrom_sources[0][1]
        return output_dir / "diagnostics" / f"metadata.chr{chrom}.json"
    return output_dir / "diagnostics" / "metadata.json"


def _expected_ref_panel_output_paths(
    config: R2MatrixBuildConfig,
    chromosomes: Sequence[str],
    *,
    metadata_path: Path,
) -> list[Path]:
    """Return deterministic R2 matrix artifact paths that may be written.

    The conservative preflight list includes build-specific R2 Parquet files and
    metadata sidecars for every emitted chromosome/build pair.
    """
    out_root = Path(config.output_dir)
    paths: list[Path] = [metadata_path]
    for chrom in chromosomes:
        for build in _emitted_genome_builds(config):
            paths.extend(r2_matrix_file_paths(out_root / build, chrom))
        # The per-chromosome dropped-SNP sidecar is class-1: always written
        # for every chromosome this invocation processes, even when header-only.
        paths.append(out_root / "diagnostics" / "dropped_snps" / f"chr{chrom}_dropped.tsv.gz")
    return paths


def _empty_unified_drop_frame() -> pd.DataFrame:
    """Return the canonical empty dropped-SNP sidecar frame."""
    return _coerce_unified_drop_frame(empty_identity_drop_frame())


def _coerce_unified_drop_frame(frame: pd.DataFrame | None) -> pd.DataFrame:
    """Return dropped-SNP rows with the unified nullable sidecar schema."""
    if frame is None:
        return _empty_unified_drop_frame()
    output = frame.copy()
    if "stage" not in output.columns and not output.empty:
        output["stage"] = "ref_panel_liftover"
    output = coerce_identity_drop_frame(output)
    for column in IDENTITY_DROP_COLUMNS:
        if column not in output.columns:
            output[column] = pd.NA
    output = output.loc[:, IDENTITY_DROP_COLUMNS]
    for column in ("CHR", "SNP", "reason", "base_key", "identity_key", "allele_set", "stage"):
        output[column] = output[column].astype("string")
    output["source_pos"] = pd.to_numeric(output["source_pos"], errors="coerce").astype("Int64")
    output["target_pos"] = pd.to_numeric(output["target_pos"], errors="coerce").astype("Int64")
    return output.reset_index(drop=True)


def _concat_drop_frames(frames: Sequence[pd.DataFrame]) -> pd.DataFrame:
    """Concatenate sidecar rows while preserving the unified nullable schema."""
    non_empty = [frame for frame in frames if frame is not None and not frame.empty]
    if not non_empty:
        return _empty_unified_drop_frame()
    return _coerce_unified_drop_frame(pd.concat(non_empty, ignore_index=True))


def _ref_panel_output_family(
    output_dir: Path,
    produced_paths: Sequence[Path] = (),
    *,
    chromosomes: Sequence[str] | None = None,
) -> list[Path]:
    """Return the build-r2-matrix artifact package owned by this invocation.

    When ``chromosomes`` is ``None``, the package spans the full output
    directory contract and is appropriate for a full ``@`` chromosome-suite
    build. When a concrete chromosome list is supplied, only those chromosomes'
    R2 Parquet, metadata, dropped-SNP, and chromosome-scoped diagnostic paths are owned
    so concurrent per-chromosome runs do not treat sibling chromosome outputs as
    stale.
    """
    scoped_chromosomes = None if chromosomes is None else [str(chrom) for chrom in chromosomes]
    diagnostics_dir = output_dir / "diagnostics"
    if scoped_chromosomes is None:
        paths: list[Path] = [diagnostics_dir / "metadata.json"]
        paths.extend(sorted(diagnostics_dir.glob("metadata.chr*.json")))
    else:
        paths = [diagnostics_dir / f"metadata.chr{chrom}.json" for chrom in scoped_chromosomes]
    for build_subdir_name in ("hg19", "hg38"):
        build_subdir = output_dir / build_subdir_name
        if not build_subdir.is_dir():
            continue
        if scoped_chromosomes is None:
            for pattern in R2_MATRIX_FILE_GLOBS:
                paths.extend(sorted(build_subdir.glob(pattern)))
        else:
            for chrom in scoped_chromosomes:
                paths.extend(r2_matrix_file_paths(build_subdir, chrom))
    dropped_dir = output_dir / "diagnostics" / "dropped_snps"
    if dropped_dir.is_dir():
        if scoped_chromosomes is None:
            paths.extend(sorted(dropped_dir.glob("chr*_dropped.tsv.gz")))
        else:
            paths.extend(dropped_dir / f"chr{chrom}_dropped.tsv.gz" for chrom in scoped_chromosomes)
    if scoped_chromosomes is None:
        paths.extend(sorted(diagnostics_dir.glob("build-r2-matrix*.log")))
    else:
        paths.extend(diagnostics_dir / f"build-r2-matrix.chr{chrom}.log" for chrom in scoped_chromosomes)
    paths.extend(Path(path) for path in produced_paths)
    return paths


@report_phase(LOGGER, 'publication', 'R2 matrix metadata')
def _write_build_metadata(
    path: Path,
    *,
    config: R2MatrixBuildConfig,
    global_config: GlobalConfig,
    chromosomes: Sequence[str],
    output_paths: dict[str, list[str]],
    resolved_plink_prefixes: Sequence[str],
) -> None:
    """Write diagnostic metadata for an R2 matrix output directory."""
    output_dir = Path(config.output_dir)
    files = {
        key: [_relative_to_output_dir(item, output_dir) for item in values]
        for key, values in output_paths.items()
        if key != "metadata"
    }
    payload = {
        "artifact_type": "r2_matrix_build",
        "files": files,
        "snp_identifier": global_config.snp_identifier,
        "genome_build": global_config.genome_build,
        "source_genome_build": config.source_genome_build,
        "emitted_genome_builds": _emitted_genome_builds(config),
        "chromosomes": [str(chrom) for chrom in chromosomes],
        "plink_prefix": config.plink_prefix,
        "resolved_plink_prefixes": [str(prefix) for prefix in resolved_plink_prefixes],
        "ref_panel_snps_file": None if config.ref_panel_snps_file is None else str(config.ref_panel_snps_file),
    }
    atomic_write_json(payload, path)


def _relative_to_output_dir(path: str, output_dir: Path) -> str:
    """Return output-relative path strings when possible."""
    candidate = Path(path)
    try:
        return str(candidate.relative_to(output_dir))
    except ValueError:
        return str(candidate)


class R2MatrixBuilder:
    """Build R2 matrices from PLINK input."""

    def __init__(self, global_config: GlobalConfig | None = None) -> None:
        """Initialize the builder with shared identifier and build defaults."""
        self.global_config = global_config or get_global_config()
        self._workflow_log_path: Path | None = None

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
            output package. Direct calls through this method do not create a
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
        return self._run(config, workflow_log_path=self._workflow_log_path)

    @serial_calculation
    def _run(
        self,
        config: R2MatrixBuildConfig,
        *,
        workflow_log_path: Path | None,
    ) -> R2MatrixBuildResult:
        """Shared implementation for direct API and parsed workflow entry points."""
        print_global_config_banner(type(self).__name__, self.global_config)
        output_dir = ensure_output_directory(config.output_dir, label="output directory")
        from .._runtime.preflight import inspect_declared_inputs
        from .preflight import GENETIC_MAP_HEADER_CHECKS
        declared_inputs = inspect_declared_inputs(plink=config.plink_prefix, header_checks=GENETIC_MAP_HEADER_CHECKS,
            files=[(name, getattr(config, name, None)) for name in (
                'genetic_map_hg19_sources', 'genetic_map_hg38_sources')],
            scalar_files=[(name, getattr(config, name, None)) for name in ('keep_indivs_file', 'ref_panel_snps_file')],
            issues_path=output_dir/'diagnostics/input_issues.tsv',
            plink_issues_path=output_dir/'diagnostics/plink_input_issues.tsv')
        plink = inspect_plink_inputs(config.plink_prefix, discovered_members=declared_inputs.files['plink_members'])
        issues_path = output_dir / "diagnostics" / "plink_input_issues.tsv"
        try:
            plink.require_valid()
        except LDSCInputError:
            issues_path.parent.mkdir(parents=True, exist_ok=True)
            pd.DataFrame(plink.issues).to_csv(issues_path, sep="\t", index=False)
            LOGGER.error("PLINK input validation failed; complete repair audit: %s", issues_path)
            raise
        resolved_prefixes = plink.prefixes
        chrom_sources = [(prefix, chrom) for chrom, prefix in plink.chromosome_prefixes.items()]
        workflow_log_path = _resolve_build_ref_panel_log_path(workflow_log_path, config, chrom_sources)
        if workflow_log_path is not None:
            scope = [chrom for _, chrom in chrom_sources] if _is_chromosome_scoped_invocation(config, chrom_sources) else None
            preflight_output_artifact_family(
                [workflow_log_path], _ref_panel_output_family(output_dir, [workflow_log_path], chromosomes=scope),
                overwrite=config.overwrite, label="R2 matrix output artifact",
            )
        try:
            config = self._resolve_source_genome_build(config, resolved_prefixes)
        except (LDSCInputError, LDSCConfigError):
            # Output builds are still unknown; the log alone has been preflighted.
            with workflow_logging("build-r2-matrix", workflow_log_path, log_level=self.global_config.log_level):
                raise
        snp_identifier_mode = normalize_snp_identifier_mode(self.global_config.snp_identifier)

        chromosomes = [chrom for _, chrom in chrom_sources]
        metadata_path = _resolve_build_metadata_path(output_dir, config, chrom_sources)
        produced_paths = _expected_ref_panel_output_paths(config, chromosomes, metadata_path=metadata_path)
        if workflow_log_path is not None:
            produced_paths.append(workflow_log_path)
        owned_chromosomes = chromosomes if _is_chromosome_scoped_invocation(config, chrom_sources) else None
        stale_paths = preflight_output_artifact_family(
            produced_paths,
            _ref_panel_output_family(output_dir, produced_paths, chromosomes=owned_chromosomes),
            overwrite=config.overwrite,
            label="R2 matrix output artifact",
        )
        (output_dir / "diagnostics").mkdir(parents=True, exist_ok=True)
        issues_path.unlink(missing_ok=True)

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
            _write_build_metadata(
                metadata_path,
                config=config,
                global_config=self.global_config,
                chromosomes=[chrom for chrom, _ in chrom_records],
                output_paths=output_paths,
                resolved_plink_prefixes=resolved_prefixes,
            )
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
        if config.source_genome_build not in {"auto", "hg19", "hg38"}:
            raise LDSCConfigError(
                f"build-r2-matrix received unsupported source_genome_build={config.source_genome_build!r}. "
                "Most likely a Python caller bypassed CLI choices with an invalid build token. "
                "Use `source_genome_build='auto'`, `'hg19'`, or `'hg38'`."
            )
        from .._genome.build_policy import GenomeBuildCheck

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
        from .._genome.build_policy import GenomeBuildCheck
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
        restriction_mode = None
        restriction_values = None
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
            restriction_values = restriction_keys.keys
            message = (
                f"Loaded {len(restriction_values)} SNP restriction identifiers "
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
            genetic_map_hg19=None if hg19_files is None else kernel_builder.load_genetic_map_group(hg19_files),
            genetic_map_hg38=None if hg38_files is None else kernel_builder.load_genetic_map_group(hg38_files),
            liftover_chain_paths={
                ("hg19", "hg38"): config.liftover_chain_hg19_to_hg38_file,
                ("hg38", "hg19"): config.liftover_chain_hg38_to_hg19_file,
            },
            restriction_mode=restriction_mode,
            restriction_values=restriction_values,
            restriction_keys=restriction_keys,
        )


    def _check_chromosome_selections(self, chrom_sources, config, build_state, issues_path):
        """Run SNP selection on every chromosome before any genotype reading.

        Autosomes left without SNPs are collected into one error, because a
        R2 matrix missing an autosome silently changes downstream LD scores. Other
        chromosomes left without SNPs are skipped with a warning. Returns the
        ``(prefix, chrom)`` pairs that keep SNPs.
        """
        from .._runtime.preflight import InputGate, input_issue
        from .._runtime.progress import PhaseProgress

        issues, buildable = [], []
        with PhaseProgress(LOGGER, 'validation', 'reference-panel SNP selection', len(chrom_sources)) as progress:
            for prefix, chrom in chrom_sources:
                selection = self._select_chromosome_snps(prefix, chrom, config, build_state)
                progress.advance(object=f'chromosome {chrom}')
                if selection.empty_stage is None:
                    buildable.append((prefix, chrom))
                elif chrom in _AUTOSOMES:
                    details = f"all {selection.n_bim_snps} SNPs were removed at {selection.empty_stage}"
                    issues.append(dict(input_issue('reference panel PLINK input', prefix, chrom, 'no_snps_remain', details),
                                       repair=_EMPTIED_AUTOSOME_REPAIR))
                else:
                    LOGGER.warning(
                        f"Skipping chromosome {chrom}: all {selection.n_bim_snps} SNPs were removed at "
                        f"{selection.empty_stage}; the R2 matrix will not include this chromosome."
                    )
        InputGate('Reference-panel SNP selection', issues_path, issues,
                  troubleshooting=_BUILD_REF_PANEL_EMPTY_DOC).finish()
        return buildable

    def _select_chromosome_snps(
        self,
        prefix: str,
        chrom: str,
        config: R2MatrixBuildConfig,
        build_state: _BuildState,
    ) -> _ChromosomeSelection:
        """Apply the BIM-level SNP selection stages to one chromosome.

        Stages are SNP restriction, identity cleanup, and liftover plus
        duplicate-position filtering; none reads genotypes. The dropped-SNP
        sidecar is always written. Reason rows cover source identity cleanup
        and liftover-stage drops; restriction-file drops are outside that
        audit vocabulary.
        """
        sidecar_path = Path(config.output_dir) / "diagnostics" / "dropped_snps" / f"chr{chrom}_dropped.tsv.gz"
        try:
            bim = legacy_parse.PlinkBIMFile(prefix + ".bim")
        except ValueError as exc:
            if "Usecols do not match columns" in str(exc):
                raise LDSCInputError(
                    f"build-r2-matrix could not read PLINK BIM file '{prefix}.bim': "
                    "A1 and A2 allele columns are required. Most likely the BIM file "
                    "has fewer than six standard PLINK columns. Provide a complete "
                    ".bim file with CHR/SNP/CM/BP/A1/A2 columns."
                ) from exc
            raise
        fam = legacy_parse.PlinkFAMFile(prefix + ".fam")
        panel_df, _report = normalize_chr_pos_frame(
            bim.df,
            context=f"{prefix}.bim",
            chr_col="CHR",
            pos_col="BP",
            coordinate_policy="drop",
            logger=LOGGER,
        )
        chrom_df = panel_df.loc[panel_df["CHR"] == chrom].copy()
        if len(chrom_df) == 0:
            raise LDSCInputError(
                f"build-r2-matrix found no SNPs for chromosome {chrom} in PLINK prefix "
                f"'{prefix}'. Most likely the `.bim` file does not contain that chromosome "
                "after chromosome normalization. Pass the correct chromosome-specific PLINK "
                "prefix or rebuild the PLINK panel with matching chromosome labels."
            )

        metadata_columns = ["CHR", "SNP", "BP", "A1", "A2"]
        chrom_metadata = chrom_df.loc[:, metadata_columns].copy().rename(columns={"BP": "POS"})
        chrom_metadata["POS"] = chrom_metadata["POS"].astype(int)
        chrom_metadata["_plink_row_index"] = chrom_metadata.index.astype(int)
        keep_snps = chrom_df.index.to_numpy(dtype=int)
        if (
            build_state.restriction_values is not None
            and build_state.restriction_mode is not None
        ):
            active_restriction = build_state.restriction_keys or build_state.restriction_values
            keep_mask = kernel_builder.build_restriction_mask(
                chrom_metadata.reset_index(drop=True),
                active_restriction,
                build_state.restriction_mode,
            )
            # Use _plink_row_index, not chrom_metadata.index: region exclusion may have
            # compacted the frame, so the positional index no longer maps to PLINK rows;
            # the column preserves the original indices keep_mask selects against.
            keep_snps = chrom_metadata["_plink_row_index"].to_numpy(dtype=int)[keep_mask]
            if len(keep_snps) == 0:
                _write_dropped_sidecar(_empty_unified_drop_frame(), sidecar_path, chrom)
                return _ChromosomeSelection(bim, fam, chrom_df, keep_snps, sidecar_path, len(chrom_df),
                                            empty_stage="SNP restriction")

        dropped_frames: list[pd.DataFrame] = []
        identity_input = chrom_metadata.loc[keep_snps].copy()
        identity_result = clean_identity_artifact_table(
            identity_input,
            self.global_config.snp_identifier,
            context=f"R2 matrix chromosome {chrom}",
            stage="ref_panel_source_identity_cleanup",
            logger=LOGGER,
        )
        if not identity_result.dropped.empty:
            dropped_frames.append(identity_result.dropped)
            dropped_examples = (
                identity_result.dropped.loc[:, ["SNP", "reason", "base_key", "identity_key", "allele_set"]]
                .head(10)
                .to_dict(orient="records")
            )
            message = (
                f"R2 matrix identity cleanup drops on chromosome {chrom}: {dropped_examples}"
            )
            LOGGER.debug(message)
        keep_snps = identity_result.cleaned["_plink_row_index"].to_numpy(dtype=int)
        if len(keep_snps) == 0:
            dropped_df = _concat_drop_frames(dropped_frames)
            _write_dropped_sidecar(dropped_df, sidecar_path, chrom)
            return _ChromosomeSelection(bim, fam, chrom_df, keep_snps, sidecar_path, len(chrom_df),
                                        empty_stage="identity cleanup")

        keep_snps, hg19_lookup, hg38_lookup, liftover_drop_frame = self._resolve_mappable_snp_positions(
            build_state=build_state,
            chrom=chrom,
            source_build=config.source_genome_build,
            chrom_df=chrom_df,
            keep_snps=keep_snps,
            sidecar_path=sidecar_path,
        )
        if not liftover_drop_frame.empty:
            dropped_frames.append(liftover_drop_frame)
        if identity_mode_family(self.global_config.snp_identifier) == "chr_pos":
            before_liftover_collision_count = len(keep_snps)
            keep_snps, liftover_collision_df = _resolve_unique_snp_set(
                chrom=chrom,
                chrom_df=chrom_df,
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
        dropped_df = _concat_drop_frames(dropped_frames)
        _write_dropped_sidecar(dropped_df, sidecar_path, chrom)
        return _ChromosomeSelection(
            bim, fam, chrom_df, keep_snps, sidecar_path, len(chrom_df), hg19_lookup, hg38_lookup,
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
        sidecar_path = selection.sidecar_path

        keep_indivs = None
        if config.keep_indivs_file is not None:
            keep_path = resolve_scalar_path(config.keep_indivs_file, label="individual keep file")
            keep_indivs = fam.loj(legacy_parse.FilterFile(keep_path).IDList)

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
                keep_indivs=None if keep_indivs is None else list(keep_indivs),
                mafMin=config.maf_min,
                streaming=build_state.restriction_values is None,
            )
            metadata = kernel_builder.build_plink_metadata_frame(
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

            geno_kept_snps = list(geno.kept_snps)
            hg19_positions = self._positions_from_lookup(geno_kept_snps, hg19_lookup) if hg19_lookup else None
            hg38_positions = self._positions_from_lookup(geno_kept_snps, hg38_lookup) if hg38_lookup else None
            build_positions = hg19_positions if build == "hg19" else hg38_positions
            if build_positions is None:
                raise LDSCInputError(
                    f"build-r2-matrix cannot emit {build} artifacts for chromosome {chrom}: "
                    "positions for that build are unavailable. Most likely the matching "
                    "liftover chain was omitted or dropped all retained SNPs. Provide the "
                    "matching liftover chain or disable opposite-build emission. "
                    f"Other causes & fixes: {_BUILD_REF_PANEL_LIFTOVER_DOC}"
                )
            _validate_emitted_build_chr_pos_uniqueness(
                metadata=metadata,
                positions=build_positions,
                genome_build=build,
                chrom=chrom,
                snp_identifier=self.global_config.snp_identifier,
            )

            cm_values = (
                kernel_builder.interpolate_genetic_map_cm(chrom, build_positions, build_state.genetic_map_hg19)
                if build == "hg19" and build_state.genetic_map_hg19 is not None
                else kernel_builder.interpolate_genetic_map_cm(chrom, build_positions, build_state.genetic_map_hg38)
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
            flip_sign = kernel_builder.orientation_flip_sign(geno.freq)
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
                ld_window_mode=(
                    "snps"
                    if config.ld_wind_snps is not None
                    else "kb"
                    if config.ld_wind_kb is not None
                    else "cm"
                ),
                ld_window_value=(
                    config.ld_wind_snps
                    if config.ld_wind_snps is not None
                    else config.ld_wind_kb
                    if config.ld_wind_kb is not None
                    else config.ld_wind_cm
                ),
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
        chrom_df,
        keep_snps,
        sidecar_path: Path | None = None,
    ) -> tuple[Any, dict[int, int], dict[int, int], pd.DataFrame]:
        """Resolve retained SNP positions and liftover-stage dropped-SNP rows.

        Returns retained PLINK row indices, source/target build position
        lookups, and a unified nullable sidecar frame for unmapped or
        cross-chromosome chain liftover drops. Source-only builds return an
        empty sidecar frame with the same schema.
        """
        keep_snps = list(keep_snps)
        candidate_positions = chrom_df.loc[keep_snps, "BP"].to_numpy(dtype=int)
        if source_build == "hg19":
            if build_state.liftover_chain_paths.get(("hg19", "hg38")) is None:
                retained_snps = np.asarray(keep_snps, dtype=int)
                hg19_lookup = {
                    int(idx): int(pos)
                    for idx, pos in zip(retained_snps, candidate_positions)
                }
                return retained_snps, hg19_lookup, {}, _empty_unified_drop_frame()
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
                return retained_snps, {}, hg38_lookup, _empty_unified_drop_frame()
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
        liftover_drop_frame = _empty_unified_drop_frame()
        if dropped:
            drop_frame = pd.DataFrame(
                {
                    "CHR": [chrom] * len(keep_snps),
                    "SNP": chrom_df.loc[keep_snps, "SNP"].astype(str).tolist(),
                    "source_POS": candidate_positions,
                },
                index=keep_snps,
            )
            unmapped_mask, cross_chrom_mask = mapping_reason_masks(mapping)
            liftover_drop_frame = _concat_drop_frames(
                [
                    _coerce_unified_drop_frame(
                        pd.DataFrame(
                            {
                                "CHR": [chrom] * int(unmapped_mask.sum()),
                                "SNP": drop_frame.loc[unmapped_mask, "SNP"].astype(str).tolist(),
                                "source_pos": drop_frame.loc[unmapped_mask, "source_POS"].tolist(),
                                "target_pos": [pd.NA] * int(unmapped_mask.sum()),
                                "reason": ["unmapped_liftover"] * int(unmapped_mask.sum()),
                            }
                        )
                    ),
                    _coerce_unified_drop_frame(
                        pd.DataFrame(
                            {
                                "CHR": [chrom] * int(cross_chrom_mask.sum()),
                                "SNP": drop_frame.loc[cross_chrom_mask, "SNP"].astype(str).tolist(),
                                "source_pos": drop_frame.loc[cross_chrom_mask, "source_POS"].tolist(),
                                "target_pos": [pd.NA] * int(cross_chrom_mask.sum()),
                                "reason": ["cross_chromosome_liftover"] * int(cross_chrom_mask.sum()),
                            }
                        )
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
        if source_build == target_build:
            array = np.asarray(positions, dtype=int)
            return LiftOverMappingResult(
                translated_positions=array.copy(),
                keep_mask=np.ones(len(array), dtype=bool),
                unmapped_count=0,
                cross_chrom_count=0,
            )
        key = (source_build, target_build)
        translator = build_state.translator_cache.get(key)
        if translator is None:
            translator = LiftOverTranslator(
                source_build=source_build,
                target_build=target_build,
                chain_path=build_state.liftover_chain_paths.get(key),
                chain_flag_hint=_chain_flag_hint(source_build, target_build),
                workflow_label="R2 matrix liftover",
            )
            build_state.translator_cache[key] = translator
        return translator.map_positions(chrom, positions)

    def _positions_from_lookup(self, kept_snps: Sequence[int], lookup: dict[int, int]):
        """Materialize retained positions in the same order as ``kept_snps``."""
        return np.asarray([lookup[int(idx)] for idx in kept_snps], dtype=int)


def _chain_flag_hint(source_build: str, target_build: str) -> str:
    """Return the CLI chain flag that corresponds to an R2 matrix build pair."""
    return f"--liftover-chain-{source_build}-to-{target_build}-file"


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
    coords, max_dist = window_coordinates(
        positions, cm_values,
        ld_wind_snps=config.ld_wind_snps, ld_wind_kb=config.ld_wind_kb, ld_wind_cm=config.ld_wind_cm,
    )
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
    chrom_df: pd.DataFrame,
    keep_snps: np.ndarray,
    hg19_lookup: dict[int, int],
    hg38_lookup: dict[int, int],
    sidecar_path: Path | None = None,
) -> tuple[np.ndarray, pd.DataFrame]:
    """Detect and resolve duplicate coordinate groups for chr_pos R2 matrix builds.

    Source duplicates are detectable before liftover. Target collisions are
    checked after position lookups are available. Duplicate groups are always
    resolved with drop-all, and returned rows use the shared dropped-SNP
    sidecar schema.
    """
    keep_list = keep_snps.tolist()
    drop_set: set[int] = set()
    prov_rows: list[dict[str, object]] = []

    source_pos = chrom_df.loc[keep_list, "BP"].to_numpy(dtype=int)
    source_frame = pd.DataFrame(
        {
            "CHR": [chrom] * len(keep_list),
            "SNP": chrom_df.loc[keep_list, "SNP"].astype(str).tolist(),
            "POS": source_pos,
        },
        index=keep_list,
    )
    source_dup_result = duplicate_coordinate_drop_result(
        source_frame,
        source_pos_col="POS",
        reason="source_duplicate",
    )
    source_dup_mask = ~source_dup_result.keep_mask
    if source_dup_mask.any():
        log_liftover_drop_report(
            LOGGER,
            source_dup_result.report,
            workflow_label="R2 matrix liftover",
            sidecar_path=sidecar_path,
        )
        for i, idx in enumerate(keep_list):
            if source_dup_mask[i]:
                drop_set.add(int(idx))
                prov_rows.append(
                    {
                        "CHR": chrom,
                        "SNP": chrom_df.loc[idx, "SNP"],
                        "source_pos": int(source_pos[i]),
                        "target_pos": pd.NA,
                        "reason": "source_duplicate",
                    }
                )

    source_survivors = [int(idx) for idx in keep_list if idx not in drop_set]
    for lookup in (hg19_lookup, hg38_lookup):
        if not lookup or len(source_survivors) < 2:
            continue
        target_pos = [lookup[idx] for idx in source_survivors]
        target_frame = pd.DataFrame(
            {
                "CHR": [chrom] * len(source_survivors),
                "SNP": chrom_df.loc[source_survivors, "SNP"].astype(str).tolist(),
                "source_pos": [int(chrom_df.loc[idx, "BP"]) for idx in source_survivors],
                "target_pos": target_pos,
            },
            index=source_survivors,
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
        for i, idx in enumerate(source_survivors):
            if target_dup_mask[i] and idx not in drop_set:
                drop_set.add(int(idx))
                prov_rows.append(
                    {
                        "CHR": chrom,
                        "SNP": chrom_df.loc[idx, "SNP"],
                        "source_pos": int(chrom_df.loc[idx, "BP"]),
                        "target_pos": int(target_pos[i]),
                        "reason": "target_collision",
                    }
                )

    cleaned = np.asarray([idx for idx in keep_snps if idx not in drop_set], dtype=keep_snps.dtype)
    columns = ["CHR", "SNP", "source_pos", "target_pos", "reason"]
    if prov_rows:
        dropped_df = _coerce_unified_drop_frame(pd.DataFrame(prov_rows, columns=columns))
        dropped_df["target_pos"] = dropped_df["target_pos"].astype("Int64")
    else:
        dropped_df = _empty_unified_drop_frame()
    return cleaned, dropped_df


def _write_dropped_sidecar(dropped_df: pd.DataFrame, path: Path, chrom: str) -> None:
    """Write the always-owned dropped-SNP audit sidecar for one chromosome."""
    path.parent.mkdir(parents=True, exist_ok=True)
    dropped_df = _coerce_unified_drop_frame(dropped_df)
    dropped_df.to_csv(path, sep="\t", index=False, compression="gzip", na_rep="")
    n_dropped = len(dropped_df)
    if n_dropped == 0:
        LOGGER.info(f"No SNPs dropped on chromosome {chrom}; audit sidecar at '{path}'.")
        return
    counts = dropped_df["reason"].value_counts(sort=False)
    count_text = ", ".join(f"{reason}={int(count)}" for reason, count in counts.items())
    LOGGER.info(
        f"R2 matrix dropped SNPs: {n_dropped} SNPs on chromosome {chrom} "
        f"({count_text}); audit sidecar at '{path}'."
    )


def build_parser() -> argparse.ArgumentParser:
    """Build the feature parser for R2 matrix generation."""
    parser = argparse.ArgumentParser(
        description="Build an R2 matrix from a PLINK reference input.",
        allow_abbrev=False,
    )
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
        '--snp-identifier', default='chr_pos_allele_aware', choices=('rsid', 'rsid_allele_aware', 'chr_pos', 'chr_pos_allele_aware'),
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
            'File listing one individual ID (IID) per row to retain from PLINK genotypes. If omitted, retain all '
            'individuals. '
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

    window.add_argument(
        '--ld-wind-cm', default=None, type=float, metavar='CM',
        help=(
            'Maximum pair distance on either side of each SNP, in centiMorgans. Requires '
            '--genetic-map-hg19-sources or --genetic-map-hg38-sources for each emitted build. No default; '
            'specify exactly one of --ld-wind-snps, --ld-wind-kb, or --ld-wind-cm.'
        ),
    )
    window.add_argument(
        '--ld-wind-kb', default=None, type=float, metavar='KB',
        help=(
            'Maximum pair distance on either side of each SNP, in kilobases. No default; specify exactly '
            'one of --ld-wind-snps, --ld-wind-kb, or --ld-wind-cm.'
        ),
    )
    window.add_argument(
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

    output.add_argument(
        '--output-dir', required=True, metavar='DIR',
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

    logs.add_argument(
        '--log-level', default='INFO', choices=('DEBUG', 'INFO', 'WARNING', 'ERROR'),
        help=LOG_LEVEL_HELP,
    )
    return parser


def config_from_args(args: argparse.Namespace) -> tuple[R2MatrixBuildConfig, GlobalConfig]:
    """Normalize CLI args into public config objects.

    ``--ref-panel-snps-file`` uses the parsed ``--snp-identifier`` value,
    falling back to the registered :class:`GlobalConfig` identifier mode only
    for programmatic namespaces that omit the attribute.
    The PLINK source build remains local to this workflow and is inferred from
    ``.bim`` coordinates during :meth:`R2MatrixBuilder.run` when set to
    ``"auto"``.
    ``GlobalConfig.genome_build`` is ignored by ``build-r2-matrix``.
    """
    registered_config = get_global_config()
    snp_identifier = normalize_snp_identifier_mode(args.snp_identifier or registered_config.snp_identifier)
    source_genome_build = normalize_genome_build(args.source_genome_build)
    if source_genome_build not in {"auto", "hg19", "hg38"}:
        raise LDSCConfigError(
            f"build-r2-matrix received invalid `--source-genome-build={args.source_genome_build}`. "
            "Most likely the source build was misspelled or is unsupported. Pass "
            "`--source-genome-build auto`, `hg19`, or `hg38`."
        )

    build_config = R2MatrixBuildConfig(
        plink_prefix=args.plink_prefix,
        source_genome_build=source_genome_build,
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
        r2_min=getattr(args, "r2_min", 0.0),
    )
    global_config = GlobalConfig(
        snp_identifier=snp_identifier,
        genome_build="auto" if identity_mode_family(snp_identifier) == "chr_pos" else None,
        log_level=args.log_level,
    )
    return build_config, global_config


def _reference_panel_failure_marker_scope(
    args: argparse.Namespace,
) -> tuple[str | os.PathLike[str], bool, str] | None:
    """Return the full-suite or concrete-chromosome overwrite marker scope."""
    output_dir = getattr(args, "output_dir", None)
    if output_dir is None:
        return None
    return output_dir, bool(getattr(args, "overwrite", False)), failure_marker_name(getattr(args, "plink_prefix", None))


def failure_marker_name(plink_prefix: str | None) -> str:
    """Return the failure marker of a build-r2-matrix run on ``plink_prefix``.

    A prefix whose file name names exactly one chromosome (``panel.chr22``)
    builds that chromosome alone, so its failed overwrite leaves
    ``RUN_FAILED.chr<chrom>.txt``; any other prefix, including a ``@`` suite,
    leaves ``RUN_FAILED.txt``. Inputs are not resolved or read.
    """
    prefix = str(plink_prefix or "")
    if "@" not in prefix:
        matches = re.findall(
            r"(?:^|[._-])(?:chr)?(1[0-9]|2[0-2]|[1-9]|X|Y)(?=$|[._-])",
            Path(prefix).name,
            flags=re.IGNORECASE,
        )
        chromosomes = {match.upper() for match in matches}
        if len(chromosomes) == 1:
            return f"RUN_FAILED.chr{next(iter(chromosomes))}.txt"
    return FAILURE_MARKER


@materializing_overwrite_guard(
    _reference_panel_failure_marker_scope,
    command="run_build_r2_matrix_from_args(...)",
)
def run_build_r2_matrix_from_args(args: argparse.Namespace) -> R2MatrixBuildResult:
    """Run R2 matrix building from parsed CLI arguments.

    The parsed workflow preflights all deterministic R2 matrix artifacts plus the
    selected build-r2-matrix log before chromosome processing. The returned
    ``R2MatrixBuildResult`` contains R2 matrix artifact paths only.
    """

    build_config, global_config = config_from_args(args)
    builder = R2MatrixBuilder(global_config=global_config)
    builder._workflow_log_path = Path(build_config.output_dir) / "diagnostics" / "build-r2-matrix.log"
    return builder.run(build_config)


def run_build_r2_matrix(**kwargs: Any) -> R2MatrixBuildResult:
    """Run the R2 matrix builder with CLI-equivalent keyword arguments.

    The wrapper accepts the same modern I/O and build names as
    ``ldsc build-r2-matrix``. The shared identifier mode comes from the
    registered ``GlobalConfig``; call ``set_global_config(...)`` before
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
        "run_build_r2_matrix", build_parser(), kwargs, rejected=rejected, log_level=get_global_config().log_level,
    )
    return run_build_r2_matrix_from_args(args)


def main(argv: Sequence[str] | None = None) -> R2MatrixBuildResult:
    """Command-line entry point for the R2 matrix builder."""

    parser = build_parser()
    args = parser.parse_args(argv)
    return run_build_r2_matrix_from_args(args)
