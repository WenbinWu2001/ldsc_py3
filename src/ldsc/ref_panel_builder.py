"""Workflow wrapper for building standard parquet reference panels from PLINK.

Core functionality:
    Resolve public PLINK, optional genetic-map, and optional liftover inputs;
    orchestrate chromosome-wise reference-panel builds; and emit the standard
    parquet artifacts consumed by the refactored LDSC workflows.

Overview
--------
This module is the public entry point for the `build-ref-panel` workflow. It
keeps path resolution, logging configuration, and cross-file validation in the
workflow layer, then delegates pairwise-R2 generation and parquet serialization to
``ldsc._kernel.ref_panel_builder``. Before chromosome processing begins, the
workflow precomputes deterministic parquet and metadata sidecar destinations
and refuses existing files unless ``overwrite=True`` was configured. Parsed
workflow wrappers write ``build-ref-panel.log`` for multi-chromosome runs and
chromosome-scoped logs for concrete single-chromosome runs; direct builder
calls return data artifact paths only and do not create a log file by default.

Unlike the user-facing result-directory writers, ``build-ref-panel`` keeps an
expert-oriented overwrite contract: ``overwrite=True`` permits replacing
candidate artifacts for the current chromosome/build set but does not remove
stale optional target-build or ``dropped_snps`` siblings from earlier
configurations. Use a fresh output directory when changing emitted builds,
liftover/coordinate configuration, or chromosome scope.
Reference-panel liftover keeps the historical source-build plus optional
opposite-build UX, but matching chain files are valid only when the active
identifier mode is ``chr_pos``. Duplicate-position handling is likewise
``chr_pos``-only and uses drop-all for coordinate collisions.
"""

from __future__ import annotations

import argparse
import logging
from dataclasses import asdict, dataclass, field, replace as dataclass_replace
from pathlib import Path
from typing import Any, Sequence

import numpy as np
import pandas as pd

from .column_inference import (
    RESTRICTION_CHRPOS_SPEC_MAP,
    RESTRICTION_HG19_POS_SPEC,
    RESTRICTION_HG38_POS_SPEC,
    normalize_column_token,
    normalize_snp_identifier_mode,
    resolve_required_column,
)
from ._kernel.snp_identity import identity_mode_family
from .config import GlobalConfig, ReferencePanelBuildConfig, get_global_config, print_global_config_banner
from ._coordinates import normalize_chr_pos_frame
from .genome_build_inference import resolve_genome_build
from .hm3 import packaged_hm3_curated_map_path
from .path_resolution import (
    ensure_output_directory,
    ensure_output_paths_available,
    resolve_file_group,
    resolve_plink_prefix_group,
    resolve_scalar_path,
)
from ._logging import configure_package_logging, log_inputs, log_outputs, workflow_logging
from ._kernel import formats as legacy_parse
from ._kernel import identifiers as kernel_identifiers
from ._kernel import ldscore as kernel_ldscore
from ._kernel import ref_panel_builder as kernel_builder
from ._kernel.liftover import (
    Hm3DualBuildLifter,
    LIFTOVER_DROP_COLUMNS,
    duplicate_coordinate_drop_result,
    liftover_drop_report,
    log_liftover_drop_report,
    mapping_reason_masks,
)


LOGGER = logging.getLogger("LDSC.ref_panel_builder")
_GENETIC_MAP_SUFFIXES = ("", ".txt", ".txt.gz", ".tsv", ".tsv.gz", ".csv", ".csv.gz")
_TABLE_SUFFIXES = ("", ".txt", ".txt.gz", ".tsv", ".tsv.gz", ".csv", ".csv.gz")


@dataclass(frozen=True)
class ReferencePanelBuildResult:
    """Summary of files produced by one reference-panel build run.

    Parameters
    ----------
    panel_name : str
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

    panel_name: str
    chromosomes: list[str]
    output_paths: dict[str, list[str]] = field(default_factory=dict)
    config_snapshot: dict[str, Any] = field(default_factory=dict)


@dataclass
class _BuildState:
    """Shared resources reused across chromosomes during one build run."""
    genetic_map_hg19: Any | None
    genetic_map_hg38: Any | None
    liftover_chain_paths: dict[tuple[str, str], str | None] = field(default_factory=dict)
    use_hm3_quick_liftover: bool = False
    hm3_map_file: str | None = None
    restriction_mode: str | None = None
    restriction_values: set[str] | None = None
    translator_cache: dict[tuple[str, str], kernel_builder.LiftOverTranslator] = field(default_factory=dict)


def _emitted_genome_builds(config: ReferencePanelBuildConfig) -> list[str]:
    """Return source and optional target builds emitted by one build config."""
    source_build = config.source_genome_build
    if source_build not in {"hg19", "hg38"}:
        raise ValueError("source_genome_build must be resolved before output paths are generated.")
    target_build = "hg38" if source_build == "hg19" else "hg19"
    if config.use_hm3_quick_liftover:
        return [source_build, target_build]
    matching_chain = (
        config.liftover_chain_hg19_to_hg38_file
        if source_build == "hg19"
        else config.liftover_chain_hg38_to_hg19_file
    )
    return [source_build] if matching_chain is None else [source_build, target_build]


def _resolve_build_ref_panel_log_path(
    workflow_log_path: Path | None,
    config: ReferencePanelBuildConfig,
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


def _expected_ref_panel_output_paths(config: ReferencePanelBuildConfig, chromosomes: Sequence[str]) -> list[Path]:
    """Return deterministic reference-panel artifact paths that may be written.

    The conservative preflight list includes build-specific R2 parquet and
    metadata sidecars for every emitted chromosome/build pair.
    """
    out_root = Path(config.output_dir)
    paths: list[Path] = []
    for chrom in chromosomes:
        for build in _emitted_genome_builds(config):
            paths.extend(
                [
                    out_root / build / f"chr{chrom}_r2.parquet",
                    out_root / build / f"chr{chrom}_meta.tsv.gz",
                ]
            )
        # The per-chromosome dropped-SNP sidecar is class-1: always written
        # for every chromosome this invocation processes, even when header-only.
        paths.append(out_root / "dropped_snps" / f"chr{chrom}_dropped.tsv.gz")
    return paths


def _empty_unified_drop_frame() -> pd.DataFrame:
    """Return the canonical empty dropped-SNP sidecar frame."""
    return pd.DataFrame(
        {
            "CHR": pd.Series(dtype="string"),
            "SNP": pd.Series(dtype="string"),
            "source_pos": pd.Series(dtype="Int64"),
            "target_pos": pd.Series(dtype="Int64"),
            "reason": pd.Series(dtype="string"),
        },
        columns=LIFTOVER_DROP_COLUMNS,
    )


def _coerce_unified_drop_frame(frame: pd.DataFrame | None) -> pd.DataFrame:
    """Return dropped-SNP rows with the unified nullable sidecar schema."""
    if frame is None:
        return _empty_unified_drop_frame()
    output = frame.copy()
    for column in LIFTOVER_DROP_COLUMNS:
        if column not in output.columns:
            output[column] = pd.NA
    output = output.loc[:, LIFTOVER_DROP_COLUMNS]
    output["CHR"] = output["CHR"].astype("string")
    output["SNP"] = output["SNP"].astype("string")
    output["source_pos"] = pd.to_numeric(output["source_pos"], errors="coerce").astype("Int64")
    output["target_pos"] = pd.to_numeric(output["target_pos"], errors="coerce").astype("Int64")
    output["reason"] = output["reason"].astype("string")
    return output.reset_index(drop=True)


def _concat_drop_frames(frames: Sequence[pd.DataFrame]) -> pd.DataFrame:
    """Concatenate sidecar rows while preserving the unified nullable schema."""
    non_empty = [frame for frame in frames if frame is not None and not frame.empty]
    if not non_empty:
        return _empty_unified_drop_frame()
    return _coerce_unified_drop_frame(pd.concat(non_empty, ignore_index=True))


def _warn_stale_class2_artifacts(output_dir: Path, expected_paths: Sequence[Path]) -> None:
    """Warn when prior ref-panel parquet/meta artifacts are outside this run.

    The helper only inspects workflow-owned class-2 territory under ``hg19/``
    and ``hg38/``. It never deletes files; it emits one warning so users can
    detect accidental mixed-output directories before chromosome processing.
    """
    expected = {str(Path(path)) for path in expected_paths}
    stale: list[Path] = []
    for build_subdir_name in ("hg19", "hg38"):
        build_subdir = output_dir / build_subdir_name
        if not build_subdir.is_dir():
            continue
        for pattern in ("chr*_r2.parquet", "chr*_meta.tsv.gz"):
            for path in build_subdir.glob(pattern):
                if str(path) not in expected:
                    stale.append(path)
    if not stale:
        return
    stale_sorted = sorted(stale, key=lambda path: str(path))
    shown = stale_sorted[:10]
    path_lines = "\n".join(f"  {path}" for path in shown)
    if len(stale_sorted) > 10:
        path_lines = f"{path_lines}\n  ... and {len(stale_sorted) - 10} more"
    LOGGER.warning(
        "Found %d pre-existing reference-panel artifacts in '%s' that this run "
        "will NOT touch. If those came from a previous run with different "
        "configuration (different chromosome scope, source build, liftover "
        "settings, MAF filter, etc.), downstream tools that load this directory "
        "as a complete reference panel will silently see a mix of two different "
        "builds. If this is an intentional chromosome-by-chromosome batch with "
        "consistent configuration, no action is needed.\n\n"
        "To remove the silent-mix risk, either:\n"
        "  (a) start fresh with a new --output-dir, or\n"
        "  (b) delete the listed files before re-running.\n\n"
        "Pre-existing artifacts not produced by this run:\n%s",
        len(stale_sorted),
        output_dir,
        path_lines,
    )


class ReferencePanelBuilder:
    """Build standard parquet reference-panel artifacts from PLINK input."""

    def __init__(self, global_config: GlobalConfig | None = None) -> None:
        """Initialize the builder with shared identifier and build defaults."""
        self.global_config = global_config or get_global_config()
        self._workflow_log_path: Path | None = None

    def run(self, config: ReferencePanelBuildConfig) -> ReferencePanelBuildResult:
        """Build fixed-name parquet reference artifacts for every chromosome.

        Parameters
        ----------
        config : ReferencePanelBuildConfig
            Public build configuration containing one PLINK prefix token,
            optional genetic maps, optional liftover chain paths, output
            directory, and exactly one LD window. If
            ``config.ref_panel_snps_file`` or ``config.use_hm3_snps`` is set,
            this builder interprets the restriction using
            ``self.global_config.snp_identifier`` against the resolved source
            PLINK build. ``GlobalConfig.genome_build`` is not consulted by this
            workflow. Existing deterministic output
            paths are refused before chromosome processing unless
            ``config.overwrite`` is true. Direct calls through this method do
            not create a workflow log; parsed wrappers add a build-ref-panel
            workflow log through the shared logging context.

        Returns
        -------
        ReferencePanelBuildResult
            Summary with ``panel_name`` inferred from ``config.output_dir`` and
            build-specific artifact paths directly under ``{build}``. The
            result contains data artifacts only, not the workflow log path.

        Raises
        ------
        ValueError
            If no chromosome artifacts are produced or if multiple PLINK inputs
            resolve to the same chromosome.
        """
        return self._run(config, workflow_log_path=self._workflow_log_path)

    def _run(
        self,
        config: ReferencePanelBuildConfig,
        *,
        workflow_log_path: Path | None,
    ) -> ReferencePanelBuildResult:
        """Shared implementation for direct API and parsed workflow entry points."""
        print_global_config_banner(type(self).__name__, self.global_config)
        self._configure_logging()
        output_dir = ensure_output_directory(config.output_dir, label="output directory")
        resolved_prefixes = resolve_plink_prefix_group((config.plink_prefix,), allow_chromosome_suite=True)
        config = self._resolve_source_genome_build(config, resolved_prefixes)
        build_state = self._prepare_build_state(config)
        snp_identifier_mode = normalize_snp_identifier_mode(self.global_config.snp_identifier)

        chrom_sources: list[tuple[str, str]] = []
        seen_chromosomes: set[str] = set()
        for prefix in resolved_prefixes:
            for chrom in self._discover_prefix_chromosomes(prefix):
                if chrom in seen_chromosomes:
                    raise ValueError(f"Chromosome {chrom} is present in multiple PLINK inputs. Emit one source per chromosome.")
                seen_chromosomes.add(chrom)
                chrom_sources.append((prefix, chrom))

        workflow_log_path = _resolve_build_ref_panel_log_path(workflow_log_path, config, chrom_sources)
        preflight_paths = _expected_ref_panel_output_paths(config, [chrom for _, chrom in chrom_sources])
        if workflow_log_path is not None:
            preflight_paths.append(workflow_log_path)
        ensure_output_paths_available(
            preflight_paths,
            overwrite=config.overwrite,
            label="reference-panel output artifact",
        )
        # Class-1 sidecars are always written by this run. Class-2 r2/meta
        # siblings from other chromosome/build scopes are not removed here;
        # warn early so users can choose fresh output dirs or manual cleanup.
        _warn_stale_class2_artifacts(Path(config.output_dir), preflight_paths)

        with workflow_logging("build-ref-panel", workflow_log_path, log_level=self.global_config.log_level):
            log_inputs(
                plink_prefix=config.plink_prefix,
                output_dir=str(output_dir),
                source_genome_build=config.source_genome_build,
            )
            if identity_mode_family(snp_identifier_mode) == "rsid" and len(_emitted_genome_builds(config)) == 1:
                LOGGER.info(
                    "Coordinate duplicate filtering applies only for chr_pos-family snp_identifier modes; "
                    "keeping duplicate CHR/POS rows in rsid source-only reference-panel builds."
                )
            chrom_records: list[tuple[str, dict[str, str]]] = []
            for prefix, chrom in chrom_sources:
                chrom_paths = self._build_chromosome(prefix, chrom, config, build_state)
                if chrom_paths is None:
                    continue
                chrom_records.append((chrom, chrom_paths))

            if not chrom_records:
                raise ValueError("No chromosome artifacts were produced from the supplied PLINK inputs.")

            chrom_records.sort(key=lambda item: kernel_ldscore.chrom_sort_key(item[0]))
            output_keys = sorted({key for _, paths in chrom_records for key in paths})
            output_paths = {key: [paths[key] for _, paths in chrom_records] for key in output_keys}
            result = ReferencePanelBuildResult(
                panel_name=Path(config.output_dir).name,
                chromosomes=[chrom for chrom, _ in chrom_records],
                output_paths=output_paths,
                config_snapshot={
                    "build_config": asdict(config),
                    "global_config": asdict(self.global_config),
                    "resolved_plink_prefixes": list(resolved_prefixes),
                },
            )
            log_outputs(**{key: ", ".join(paths) for key, paths in output_paths.items()})
            return result

    def _configure_logging(self) -> None:
        """Configure package logging for the current build run."""
        configure_package_logging(self.global_config.log_level)

    def _resolve_source_genome_build(
        self,
        config: ReferencePanelBuildConfig,
        resolved_prefixes: Sequence[str],
    ) -> ReferencePanelBuildConfig:
        """Return a config whose PLINK source genome build is concrete."""
        if config.source_genome_build in {"hg19", "hg38"}:
            return config
        if config.source_genome_build is not None:
            raise ValueError(f"Unsupported source_genome_build: {config.source_genome_build!r}.")
        sample_frame = _plink_bim_chr_pos_frame(resolved_prefixes)
        try:
            source_build = resolve_genome_build(
                "auto",
                "chr_pos",
                sample_frame,
                context="build-ref-panel PLINK .bim",
                logger=LOGGER,
            )
        except ValueError as exc:
            details = str(exc).replace(
                "Pass --genome-build hg19 or --genome-build hg38",
                "Pass --source-genome-build hg19 or --source-genome-build hg38",
            )
            raise ValueError(
                "Unable to infer source_genome_build from PLINK .bim coordinates. "
                f"{details}"
            ) from exc
        if source_build not in {"hg19", "hg38"}:
            raise ValueError("Unable to resolve source_genome_build from PLINK .bim coordinates.")
        LOGGER.info(f"Using source_genome_build='{source_build}' for PLINK reference-panel coordinates.")
        return dataclass_replace(config, source_genome_build=source_build)

    def _prepare_build_state(self, config: ReferencePanelBuildConfig) -> _BuildState:
        """Resolve optional maps, chain paths, and source-build restrictions.

        A chain matching the resolved source build enables opposite-build
        emission in ``chr_pos`` mode. The same matching chain is rejected in
        ``rsid`` mode because row identity is the SNP label and cross-build
        coordinate emission would be ambiguous.
        """
        hg19_files = None
        hg38_files = None
        if config.genetic_map_hg19_sources is not None:
            hg19_files = resolve_file_group(
                (config.genetic_map_hg19_sources,),
                suffixes=_GENETIC_MAP_SUFFIXES,
                label="hg19 genetic map",
                allow_chromosome_suite=True,
            )
        if config.genetic_map_hg38_sources is not None:
            hg38_files = resolve_file_group(
                (config.genetic_map_hg38_sources,),
                suffixes=_GENETIC_MAP_SUFFIXES,
                label="hg38 genetic map",
                allow_chromosome_suite=True,
            )
        source_build = config.source_genome_build
        if source_build not in {"hg19", "hg38"}:
            raise ValueError("source_genome_build must be resolved before loading build-ref-panel shared state.")
        restriction_mode = None
        restriction_values = None
        restriction_source = packaged_hm3_curated_map_path() if config.use_hm3_snps else config.ref_panel_snps_file
        if restriction_source:
            restriction_path = resolve_scalar_path(
                restriction_source,
                suffixes=_TABLE_SUFFIXES,
                label="packaged HM3 SNP map" if config.use_hm3_snps else "reference-panel SNP restriction",
            )
            restriction_mode = normalize_snp_identifier_mode(self.global_config.snp_identifier)
            if identity_mode_family(restriction_mode) == "chr_pos":
                LOGGER.info(
                    f"Interpreting reference-panel SNP restriction '{restriction_path}' "
                    f"against source_genome_build='{source_build}' PLINK coordinates. "
                    "GlobalConfig.genome_build is ignored by build-ref-panel."
                )
            restriction_values = _read_ref_panel_snp_restriction(
                restriction_path,
                restriction_mode,
                source_genome_build=source_build,
            )
            message = (
                f"Loaded {len(restriction_values)} SNP restriction identifiers "
                f"in {restriction_mode} mode from '{restriction_path}'"
                f"{' via use_hm3_snps' if config.use_hm3_snps else ''}."
            )
            LOGGER.info(message)
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
        if config.use_hm3_quick_liftover and identity_mode_family(self.global_config.snp_identifier) == "rsid":
            raise ValueError("Reference-panel HM3 quick liftover is only valid for chr_pos-family snp_identifier modes.")
        if matching_chain is not None and identity_mode_family(self.global_config.snp_identifier) == "rsid":
            raise ValueError(
                "Reference-panel chain liftover is only valid for chr_pos-family snp_identifier modes. "
                "In rsid mode, omit the matching liftover chain and build source-genome coordinates only."
            )
        if config.use_hm3_quick_liftover:
            LOGGER.info(
                f"Using packaged curated HM3 map for reference-panel quick liftover "
                f"{source_build} -> {target_build}: {packaged_hm3_curated_map_path()}."
            )
        if config.use_hm3_quick_liftover or matching_chain is not None:
            if config.ld_wind_cm is not None and (
                (target_build == "hg38" and config.genetic_map_hg38_sources is None)
                or (target_build == "hg19" and config.genetic_map_hg19_sources is None)
            ):
                raise ValueError(f"{target_build} genetic map path is required when ld_wind_cm is set and liftover emits {target_build}.")
            if (target_build == "hg38" and config.genetic_map_hg38_sources is None) or (
                target_build == "hg19" and config.genetic_map_hg19_sources is None
            ):
                LOGGER.warning(
                    f"No {target_build} genetic map was provided; "
                    f"{target_build} metadata CM values will be written as NA."
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
            raise ValueError(f"{source_build} genetic map path is required when ld_wind_cm is set.")
        if (source_build == "hg38" and config.genetic_map_hg38_sources is None) or (
            source_build == "hg19" and config.genetic_map_hg19_sources is None
        ):
            LOGGER.warning(
                f"No {source_build} genetic map was provided; "
                f"{source_build} metadata CM values will be written as NA."
            )
        return _BuildState(
            genetic_map_hg19=None if hg19_files is None else kernel_builder.load_genetic_map_group(hg19_files),
            genetic_map_hg38=None if hg38_files is None else kernel_builder.load_genetic_map_group(hg38_files),
            liftover_chain_paths={
                ("hg19", "hg38"): config.liftover_chain_hg19_to_hg38_file,
                ("hg38", "hg19"): config.liftover_chain_hg38_to_hg19_file,
            },
            use_hm3_quick_liftover=config.use_hm3_quick_liftover,
            hm3_map_file=packaged_hm3_curated_map_path() if config.use_hm3_quick_liftover else None,
            restriction_mode=restriction_mode,
            restriction_values=restriction_values,
        )

    def _discover_prefix_chromosomes(self, prefix: str) -> list[str]:
        """List the normalized chromosomes present in one PLINK prefix."""
        bim = legacy_parse.PlinkBIMFile(prefix + ".bim")
        normalized, _report = normalize_chr_pos_frame(
            bim.df,
            context=f"{prefix}.bim",
            chr_col="CHR",
            pos_col="BP",
            coordinate_policy="drop",
            logger=LOGGER,
        )
        chromosomes = sorted(
            set(normalized["CHR"]),
            key=kernel_ldscore.chrom_sort_key,
        )
        return chromosomes

    def _build_chromosome(
        self,
        prefix: str,
        chrom: str,
        config: ReferencePanelBuildConfig,
        build_state: _BuildState,
    ) -> dict[str, str] | None:
        """Build all parquet and metadata artifacts for one chromosome.

        The per-chromosome dropped-SNP sidecar is always written for any
        chromosome this method is invoked on. Reason rows cover liftover-stage
        drops only; restriction-file and PLINK filtering drops are outside that
        audit vocabulary.
        """
        LOGGER.info(f"Building reference-panel artifacts for chromosome {chrom} from '{prefix}'.")
        sidecar_path = Path(config.output_dir) / "dropped_snps" / f"chr{chrom}_dropped.tsv.gz"
        bim = legacy_parse.PlinkBIMFile(prefix + ".bim")
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
            raise ValueError(f"No SNPs found for chromosome {chrom} in {prefix}.")

        chrom_metadata = chrom_df.loc[:, ["CHR", "SNP", "BP"]].copy().rename(columns={"BP": "POS"})
        chrom_metadata["POS"] = chrom_metadata["POS"].astype(int)
        keep_snps = chrom_df.index.to_numpy(dtype=int)
        if (
            build_state.restriction_values is not None
            and build_state.restriction_mode is not None
        ):
            keep_mask = kernel_builder.build_restriction_mask(
                chrom_metadata.reset_index(drop=True),
                build_state.restriction_values,
                build_state.restriction_mode,
            )
            keep_snps = chrom_df.index[keep_mask].to_numpy(dtype=int)
            if len(keep_snps) == 0:
                _write_dropped_sidecar(_empty_unified_drop_frame(), sidecar_path, chrom)
                LOGGER.info(f"Skipping chromosome {chrom} because no SNPs remain after restriction.")
                return None

        dropped_frames: list[pd.DataFrame] = []
        duplicate_policy_applies = identity_mode_family(self.global_config.snp_identifier) == "chr_pos"
        if duplicate_policy_applies:
            keep_snps, source_duplicate_df = _resolve_unique_snp_set(
                chrom=chrom,
                chrom_df=chrom_df,
                keep_snps=keep_snps,
                hg19_lookup={},
                hg38_lookup={},
                sidecar_path=sidecar_path,
            )
            if not source_duplicate_df.empty:
                dropped_frames.append(source_duplicate_df)
            if len(keep_snps) == 0:
                dropped_df = _concat_drop_frames(dropped_frames)
                _write_dropped_sidecar(dropped_df, sidecar_path, chrom)
                LOGGER.info(f"Skipping chromosome {chrom}: no SNPs remain after duplicate-position filtering.")
                return None

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
        if duplicate_policy_applies:
            keep_snps, target_duplicate_df = _resolve_unique_snp_set(
                chrom=chrom,
                chrom_df=chrom_df,
                keep_snps=keep_snps,
                hg19_lookup=hg19_lookup,
                hg38_lookup=hg38_lookup,
                sidecar_path=sidecar_path,
            )
            if not target_duplicate_df.empty:
                dropped_frames.append(target_duplicate_df)
            dropped_df = _concat_drop_frames(dropped_frames)
        else:
            dropped_df = _concat_drop_frames(dropped_frames)
        _write_dropped_sidecar(dropped_df, sidecar_path, chrom)
        if len(keep_snps) == 0:
            LOGGER.info(f"Skipping chromosome {chrom}: no SNPs remain after duplicate-position filtering.")
            return None

        keep_indivs = None
        if config.keep_indivs_file is not None:
            keep_path = resolve_scalar_path(config.keep_indivs_file, suffixes=_TABLE_SUFFIXES, label="individual keep file")
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
            geno = kernel_ldscore.PlinkBEDFile(
                prefix + ".bed",
                len(fam.IDList),
                bim,
                keep_snps=list(build_keep_snps),
                keep_indivs=None if keep_indivs is None else list(keep_indivs),
                mafMin=config.maf_min,
            )
            metadata = kernel_builder.build_plink_metadata_frame(
                bim=bim,
                kept_snps=geno.kept_snps,
                maf_values=geno.maf,
            )
            if len(metadata) == 0:
                LOGGER.info(f"Skipping chromosome {chrom} because no SNPs remain after PLINK filtering.")
                return None
            if set(metadata["CHR"]) != {chrom}:
                raise ValueError(f"PLINK filtering for chromosome {chrom} retained rows from multiple chromosomes.")

            geno_kept_snps = list(geno.kept_snps)
            hg19_positions = self._positions_from_lookup(geno_kept_snps, hg19_lookup) if hg19_lookup else None
            hg38_positions = self._positions_from_lookup(geno_kept_snps, hg38_lookup) if hg38_lookup else None
            build_positions = hg19_positions if build == "hg19" else hg38_positions
            if build_positions is None:
                raise ValueError(f"{build} positions are unavailable for chromosome {chrom}; provide a matching liftover chain.")
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
            if cm_values is None and config.ld_wind_cm is not None:
                raise ValueError(f"{build} genetic map is required for chromosome {chrom}.")
            coords, max_dist = kernel_builder.build_window_coordinates(
                metadata=_metadata_with_build_positions(metadata, build_positions),
                cm_values=cm_values,
                ld_wind_snps=config.ld_wind_snps,
                ld_wind_kb=config.ld_wind_kb,
                ld_wind_cm=config.ld_wind_cm,
            )
            block_left = kernel_builder.compute_block_left(coords, max_dist)
            reference_snp_table = kernel_builder.build_reference_snp_table(
                metadata=metadata,
                hg19_positions=hg19_positions,
                hg38_positions=hg38_positions,
            )
            r2_path = Path(config.output_dir) / build / f"chr{chrom}_r2.parquet"
            meta_path = Path(config.output_dir) / build / f"chr{chrom}_meta.tsv.gz"
            kernel_builder.write_r2_parquet(
                pair_rows=kernel_builder.yield_pairwise_r2_rows(
                    block_left=block_left,
                    snp_batch_size=config.snp_batch_size,
                    standardized_snp_getter=geno.nextSNPs,
                    m=geno.m,
                    n=geno.n,
                ),
                reference_snp_table=reference_snp_table,
                path=r2_path,
                genome_build=build,
                snp_identifier=self.global_config.snp_identifier,
                n_samples=geno.n,
            )
            runtime_metadata = kernel_builder.build_runtime_metadata_table(
                metadata=metadata,
                positions=build_positions,
                cm_values=cm_values,
            )
            kernel_builder.write_runtime_metadata_sidecar(runtime_metadata, meta_path)
            output_paths[f"r2_{build}"] = str(r2_path)
            output_paths[f"meta_{build}"] = str(meta_path)
            retained_count = len(metadata)
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
        if build_state.use_hm3_quick_liftover:
            target_build = "hg38" if source_build == "hg19" else "hg19"
            if build_state.hm3_map_file is None:
                raise ValueError("HM3 quick liftover requires a packaged HM3 map path.")
            query = pd.DataFrame(
                {
                    "CHR": [chrom] * len(keep_snps),
                    "POS": candidate_positions,
                    "SNP": chrom_df.loc[keep_snps, "SNP"].astype(str).tolist(),
                },
                index=keep_snps,
            )
            lifted, unmapped_indices = Hm3DualBuildLifter(source_build, target_build, build_state.hm3_map_file).lift(query)
            unmapped_set = {int(idx) for idx in unmapped_indices}
            retained_snps = np.asarray([int(idx) for idx in keep_snps if int(idx) not in unmapped_set], dtype=int)
            retained_set = set(retained_snps.tolist())
            source_lookup = {
                int(idx): int(pos)
                for idx, pos in zip(keep_snps, candidate_positions)
                if int(idx) in retained_set
            }
            target_lookup = {int(idx): int(lifted.loc[idx, "POS"]) for idx in retained_snps}
            if source_build == "hg19":
                hg19_lookup, hg38_lookup = source_lookup, target_lookup
            else:
                hg38_lookup, hg19_lookup = source_lookup, target_lookup
            liftover_drop_frame = _empty_unified_drop_frame()
            if len(unmapped_indices):
                unmapped_ordered = [int(idx) for idx in keep_snps if int(idx) in unmapped_set]
                liftover_drop_frame = _coerce_unified_drop_frame(
                    pd.DataFrame(
                        {
                            "CHR": [chrom] * len(unmapped_ordered),
                            "SNP": chrom_df.loc[unmapped_ordered, "SNP"].astype(str).tolist(),
                            "source_pos": chrom_df.loc[unmapped_ordered, "BP"].astype(int).tolist(),
                            "target_pos": [pd.NA] * len(unmapped_ordered),
                            "reason": ["unmapped_liftover"] * len(unmapped_ordered),
                        }
                    )
                )
                mask = np.asarray([idx in unmapped_set for idx in keep_snps], dtype=bool)
                log_liftover_drop_report(
                    LOGGER,
                    liftover_drop_report(
                        query,
                        mask,
                        reason="unmapped_liftover",
                        source_pos_col="POS",
                    ),
                    workflow_label="Reference-panel liftover",
                    sidecar_path=sidecar_path,
                )
                LOGGER.info(
                    f"Dropping {len(unmapped_set)} SNPs on chromosome {chrom} after HM3 quick liftover filtering."
                )
            return retained_snps, hg19_lookup, hg38_lookup, liftover_drop_frame
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
                workflow_label="Reference-panel liftover",
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
                workflow_label="Reference-panel liftover",
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
    ) -> kernel_builder.LiftOverMappingResult:
        """Translate positions between builds, caching translators by build pair."""
        if source_build == target_build:
            array = np.asarray(positions, dtype=int)
            return kernel_builder.LiftOverMappingResult(
                translated_positions=array.copy(),
                keep_mask=np.ones(len(array), dtype=bool),
                unmapped_count=0,
                cross_chrom_count=0,
            )
        key = (source_build, target_build)
        translator = build_state.translator_cache.get(key)
        if translator is None:
            translator = kernel_builder.LiftOverTranslator(
                source_build=source_build,
                target_build=target_build,
                chain_path=build_state.liftover_chain_paths.get(key),
                chain_flag_hint=_chain_flag_hint(source_build, target_build),
                workflow_label="reference-panel liftover",
            )
            build_state.translator_cache[key] = translator
        return translator.map_positions(chrom, positions)

    def _translate_positions(
        self,
        build_state: _BuildState,
        chrom: str,
        positions,
        source_build: str,
        target_build: str,
    ):
        """Translate positions between builds without returning drop statistics."""
        if source_build == target_build:
            return positions.copy()
        key = (source_build, target_build)
        translator = build_state.translator_cache.get(key)
        if translator is None:
            translator = kernel_builder.LiftOverTranslator(
                source_build=source_build,
                target_build=target_build,
                chain_path=build_state.liftover_chain_paths.get(key),
                chain_flag_hint=_chain_flag_hint(source_build, target_build),
                workflow_label="reference-panel liftover",
            )
            build_state.translator_cache[key] = translator
        return translator.translate_positions(chrom, positions)

    def _positions_from_lookup(self, kept_snps: Sequence[int], lookup: dict[int, int]):
        """Materialize retained positions in the same order as ``kept_snps``."""
        return np.asarray([lookup[int(idx)] for idx in kept_snps], dtype=int)


def _chain_flag_hint(source_build: str, target_build: str) -> str:
    """Return the CLI chain flag that corresponds to a reference-panel build pair."""
    return f"--liftover-chain-{source_build}-to-{target_build}-file"


def _plink_bim_chr_pos_frame(resolved_prefixes: Sequence[str]) -> pd.DataFrame:
    """Read CHR/BP coordinates from resolved PLINK ``.bim`` files."""
    frames: list[pd.DataFrame] = []
    for prefix in resolved_prefixes:
        frame = pd.read_csv(
            prefix + ".bim",
            sep=r"\s+",
            header=None,
            usecols=[0, 3],
            names=["CHR", "POS"],
        )
        if not frame.empty:
            normalized, _report = normalize_chr_pos_frame(
                frame,
                context=f"{prefix}.bim",
                coordinate_policy="drop",
                logger=LOGGER,
                min_position=0,
            )
            frames.append(normalized)
    if not frames:
        return pd.DataFrame(columns=["CHR", "POS"])
    return pd.concat(frames, ignore_index=True)


def _read_ref_panel_snp_restriction(
    path: str | Path,
    restriction_mode: str,
    *,
    source_genome_build: str,
) -> set[str]:
    """Read a builder SNP restriction file in the source PLINK build."""
    if identity_mode_family(restriction_mode) == "rsid":
        return kernel_identifiers.read_global_snp_restriction(path, restriction_mode)
    if source_genome_build not in {"hg19", "hg38"}:
        raise ValueError("source_genome_build must be resolved before CHR/POS SNP restriction loading.")
    return _read_source_build_chr_pos_restriction(Path(path), source_genome_build)


def _source_position_spec(source_genome_build: str):
    """Return the build-specific position column spec for the source build."""
    return RESTRICTION_HG19_POS_SPEC if source_genome_build == "hg19" else RESTRICTION_HG38_POS_SPEC


def _resolve_build_specific_position_column(
    header: Sequence[str],
    *,
    source_genome_build: str,
    path: Path,
) -> str | None:
    """Resolve only explicit build-specific position aliases, not generic POS."""
    spec = _source_position_spec(source_genome_build)
    alias_tokens = {normalize_column_token(alias) for alias in spec.aliases}
    matches = [
        column
        for column in header
        if (normalized := normalize_column_token(column)) in alias_tokens
        or any(normalized.endswith(alias) for alias in alias_tokens)
    ]
    if not matches:
        return None
    if len(matches) > 1:
        raise ValueError(
            f"Ambiguous {source_genome_build} position columns in reference-panel "
            f"SNP restriction '{path}': {', '.join(map(str, matches))}"
        )
    return matches[0]


def _restriction_frame_from_columns(
    *,
    path: Path,
    header: Sequence[str],
    rows: Sequence[Sequence[str]],
    chr_col: str,
    pos_col: str,
) -> pd.DataFrame:
    """Materialize CHR/POS restriction rows from selected source columns."""
    chr_idx = list(header).index(chr_col)
    pos_idx = list(header).index(pos_col)
    values: list[tuple[object, object]] = []
    for row in rows:
        if not row:
            continue
        if len(row) <= max(chr_idx, pos_idx):
            raise ValueError(f"Restriction row in {path} is missing the CHR or POS column.")
        values.append((row[chr_idx], row[pos_idx]))
    return pd.DataFrame(values, columns=["CHR", "POS"])


def _infer_generic_restriction_build(frame: pd.DataFrame, path: Path) -> str:
    """Infer the genome build represented by a generic restriction POS column."""
    inference_frame, _report = normalize_chr_pos_frame(
        frame.loc[:, ["CHR", "POS"]].copy(),
        context=f"reference-panel SNP restriction '{path}'",
        coordinate_policy="drop",
        logger=None,
    )
    try:
        inferred_build = resolve_genome_build(
            "auto",
            "chr_pos",
            inference_frame,
            context=f"reference-panel SNP restriction '{path}'",
            logger=LOGGER,
        )
    except ValueError as exc:
        raise ValueError(
            "Unable to infer the genome build of the generic POS column in "
            f"reference-panel SNP restriction '{path}'. Provide a source-build-specific "
            "position column such as hg19_POS or hg38_POS, aligned to the PLINK source build."
        ) from exc
    if inferred_build not in {"hg19", "hg38"}:
        raise ValueError(f"Unable to infer the genome build of reference-panel SNP restriction '{path}'.")
    return inferred_build


def _read_source_build_chr_pos_restriction(path: Path, source_genome_build: str) -> set[str]:
    """Read a CHR/POS restriction file that must align to the source build."""
    header, rows, _delimiter = kernel_identifiers._parse_restriction_rows(path)
    if not header:
        return set()
    chr_col = resolve_required_column(
        header,
        RESTRICTION_CHRPOS_SPEC_MAP["CHR"],
        context=str(path),
    )
    source_pos_col = _resolve_build_specific_position_column(
        header,
        source_genome_build=source_genome_build,
        path=path,
    )
    if source_pos_col is not None:
        LOGGER.info(
            f"Using source-build-specific restriction column '{source_pos_col}' "
            f"for source genome build '{source_genome_build}'."
        )
        frame = _restriction_frame_from_columns(
            path=path,
            header=header,
            rows=rows,
            chr_col=chr_col,
            pos_col=source_pos_col,
        )
        return kernel_identifiers._finalize_chr_pos_restriction_frame(frame, path=path, logger=LOGGER)

    generic_pos_col = resolve_required_column(
        header,
        RESTRICTION_CHRPOS_SPEC_MAP["POS"],
        context=str(path),
    )
    frame = _restriction_frame_from_columns(
        path=path,
        header=header,
        rows=rows,
        chr_col=chr_col,
        pos_col=generic_pos_col,
    )
    restriction_build = _infer_generic_restriction_build(frame, path)
    if restriction_build != source_genome_build:
        raise ValueError(
            "Reference-panel SNP restriction build mismatch: generic POS column "
            f"in '{path}' was inferred as {restriction_build}, but the source "
            f"reference panel build is {source_genome_build}. Provide SNP restrictions "
            "aligned to the PLINK source build, or add a source-build-specific "
            f"{source_genome_build}_POS column."
        )
    LOGGER.info(
        f"Generic POS column in reference-panel SNP restriction '{path}' "
        f"was inferred as source_genome_build='{source_genome_build}'."
    )
    return kernel_identifiers._finalize_chr_pos_restriction_frame(frame, path=path, logger=LOGGER)


def _position_lookup_for_build(
    genome_build: str,
    *,
    hg19_lookup: dict[int, int],
    hg38_lookup: dict[int, int],
) -> dict[int, int]:
    """Return the retained-position lookup for one genome build."""
    return hg19_lookup if genome_build == "hg19" else hg38_lookup


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
        context=f"Reference-panel build output {genome_build}[{chrom}]",
    )


def _sort_retained_snps_by_build_position(
    keep_snps,
    *,
    genome_build: str,
    hg19_lookup: dict[int, int],
    hg38_lookup: dict[int, int],
) -> np.ndarray:
    """Return retained PLINK indices in one build's genomic position order."""
    # Monotone position order here makes index i → position monotone, which lets
    # yield_pairwise_r2_rows flush cross-batch pairs with non-decreasing POS_1.
    keep_snps = np.asarray(keep_snps, dtype=int)
    lookup = _position_lookup_for_build(
        genome_build,
        hg19_lookup=hg19_lookup,
        hg38_lookup=hg38_lookup,
    )
    if not lookup:
        raise ValueError(f"{genome_build} positions are unavailable for retained SNP sorting.")
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
    """Detect and resolve duplicate coordinate groups for chr_pos panel builds.

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
            workflow_label="Reference-panel liftover",
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
            workflow_label="Reference-panel liftover",
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
        dropped_df = pd.DataFrame(prov_rows, columns=columns)
        dropped_df["target_pos"] = dropped_df["target_pos"].astype("Int64")
    else:
        dropped_df = pd.DataFrame(columns=columns)
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
        f"Reference-panel liftover-stage drops: {n_dropped} SNPs on chromosome {chrom} "
        f"({count_text}); audit sidecar at '{path}'."
    )


def build_parser() -> argparse.ArgumentParser:
    """Build the feature parser for reference-panel generation."""

    parser = argparse.ArgumentParser(
        description="Build standard parquet reference-panel artifacts from a PLINK reference input.",
        allow_abbrev=False,
    )
    parser.add_argument("--plink-prefix", required=True, help="PLINK prefix token for the reference panel.")
    parser.add_argument(
        "--source-genome-build",
        default=None,
        choices=("hg19", "hg37", "GRCh37", "hg38", "GRCh38"),
        help="Genome build of the PLINK coordinates. If omitted, infer it from .bim CHR/BP rows before SNP restriction.",
    )
    parser.add_argument(
        "--genetic-map-hg19-sources",
        default=None,
        help="Optional genetic map aligned to hg19; required for cM windows when hg19 output is emitted.",
    )
    parser.add_argument(
        "--genetic-map-hg38-sources",
        default=None,
        help="Optional genetic map aligned to hg38; required for cM windows when hg38 output is emitted.",
    )
    parser.add_argument(
        "--liftover-chain-hg19-to-hg38-file",
        default=None,
        help=(
            "Optional chain file for hg19->hg38 liftover. Enables hg38 outputs "
            "for hg19 source builds in chr_pos mode; invalid in rsid mode."
        ),
    )
    parser.add_argument(
        "--liftover-chain-hg38-to-hg19-file",
        default=None,
        help=(
            "Optional chain file for hg38->hg19 liftover. Enables hg19 outputs "
            "for hg38 source builds in chr_pos mode; invalid in rsid mode."
        ),
    )
    parser.add_argument("--output-dir", required=True, help="Output root directory for emitted parquet artifacts.")
    parser.add_argument("--overwrite", action="store_true", default=False, help="Replace existing candidate panel output files.")
    parser.add_argument("--ld-wind-snps", default=None, type=int, help="LD window size in SNPs.")
    parser.add_argument("--ld-wind-kb", default=None, type=float, help="LD window size in kilobases.")
    parser.add_argument("--ld-wind-cm", default=None, type=float, help="LD window size in centiMorgans.")
    parser.add_argument("--maf-min", default=None, type=float, help="Optional MAF filter for retained SNPs.")
    parser.add_argument(
        "--ref-panel-snps-file",
        default=None,
        help=(
            "Optional SNP restriction file defining the retained reference-panel universe. "
            "In chr_pos mode, coordinates must match the PLINK source genome build."
        ),
    )
    parser.add_argument(
        "--use-hm3-snps",
        action="store_true",
        default=False,
        help="Restrict the emitted reference-panel universe to the packaged curated HM3 SNP map.",
    )
    parser.add_argument(
        "--use-hm3-quick-liftover",
        action="store_true",
        default=False,
        help="Use the packaged curated HM3 map for HM3-only coordinate liftover; requires --use-hm3-snps.",
    )
    parser.add_argument(
        "--snp-identifier",
        default="chr_pos_allele_aware",
        choices=("rsid", "rsid_allele_aware", "chr_pos", "chr_pos_allele_aware"),
        help="SNP identifier mode for --ref-panel-snps-file.",
    )
    parser.add_argument("--keep-indivs-file", default=None, help="Optional individual-keep file.")
    parser.add_argument(
        "--snp-batch-size",
        dest="snp_batch_size",
        default=128,
        type=int,
        help="Number of SNPs loaded per pairwise-R2 computation batch. Larger values may improve throughput but use more memory.",
    )
    parser.add_argument(
        "--chunk-size",
        dest="snp_batch_size",
        type=int,
        help=argparse.SUPPRESS,
    )
    parser.add_argument("--log-level", default="INFO", choices=("DEBUG", "INFO", "WARNING", "ERROR"))
    return parser


def config_from_args(args: argparse.Namespace) -> tuple[ReferencePanelBuildConfig, GlobalConfig]:
    """Normalize CLI args into public config objects.

    ``--ref-panel-snps-file`` uses the parsed ``--snp-identifier`` value,
    falling back to the registered :class:`GlobalConfig` identifier mode only
    for programmatic namespaces that omit the attribute.
    The PLINK source build remains local to this workflow and is inferred from
    ``.bim`` coordinates during :meth:`ReferencePanelBuilder.run` when omitted.
    ``GlobalConfig.genome_build`` is ignored by ``build-ref-panel``.
    """
    registered_config = get_global_config()
    snp_identifier = normalize_snp_identifier_mode(args.snp_identifier or registered_config.snp_identifier)
    source_genome_build = None
    if args.source_genome_build is not None:
        source_genome_build = resolve_genome_build(
            args.source_genome_build,
            "chr_pos",
            None,
            context="build-ref-panel source",
            logger=LOGGER,
        )

    build_config = ReferencePanelBuildConfig(
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
        use_hm3_snps=getattr(args, "use_hm3_snps", False),
        use_hm3_quick_liftover=getattr(args, "use_hm3_quick_liftover", False),
        keep_indivs_file=args.keep_indivs_file,
        snp_batch_size=args.snp_batch_size,
    )
    global_config = GlobalConfig(
        snp_identifier=snp_identifier,
        genome_build="auto" if identity_mode_family(snp_identifier) == "chr_pos" else None,
        log_level=args.log_level,
    )
    return build_config, global_config


def run_build_ref_panel_from_args(args: argparse.Namespace) -> ReferencePanelBuildResult:
    """Run reference-panel building from parsed CLI arguments.

    The parsed workflow preflights all deterministic panel artifacts plus the
    selected build-ref-panel log before chromosome processing. The returned
    ``ReferencePanelBuildResult`` contains panel artifact paths only.
    """

    build_config, global_config = config_from_args(args)
    builder = ReferencePanelBuilder(global_config=global_config)
    builder._workflow_log_path = Path(build_config.output_dir) / "build-ref-panel.log"
    return builder.run(build_config)


def run_build_ref_panel(**kwargs: Any) -> ReferencePanelBuildResult:
    """Run the reference-panel builder with CLI-equivalent keyword arguments.

    The wrapper accepts the same modern I/O and build names as
    ``ldsc build-ref-panel``. The shared identifier mode comes from the
    registered ``GlobalConfig``; call ``set_global_config(...)`` before
    invoking this wrapper. ``GlobalConfig.genome_build`` is ignored here.
    """
    forbidden = sorted({"genome_build", "log_level", "snp_identifier"} & set(kwargs))
    if forbidden:
        joined = ", ".join(forbidden)
        raise ValueError(
            f"Python run_build_ref_panel() no longer accepts {joined}; call set_global_config(...) first."
        )
    removed = sorted(
        {
            "bfile",
            "out",
            "panel_label",
            "plink_path",
            "genetic_map_hg19_path",
            "genetic_map_hg38_path",
            "liftover_chain_hg19_to_hg38_path",
            "liftover_chain_hg38_to_hg19_path",
            "ref_panel_snps_path",
            "keep_indivs_path",
            "genetic_map_hg19",
            "genetic_map_hg38",
            "liftover_chain_hg19_to_hg38",
            "liftover_chain_hg38_to_hg19",
            "keep_indivs",
            "maf",
        }
        & set(kwargs)
    )
    if removed:
        joined = ", ".join(removed)
        raise ValueError(f"Python run_build_ref_panel() no longer accepts removed IO argument(s): {joined}.")

    parser = build_parser()
    defaults = vars(
        parser.parse_args(
            [
                "--plink-prefix",
                "plink/panel.@",
                "--genetic-map-hg19-sources",
                "hg19.map",
                "--genetic-map-hg38-sources",
                "hg38.map",
                "--output-dir",
                "out",
                "--ld-wind-kb",
                "1",
            ]
        )
    )
    global_config = get_global_config()
    defaults["log_level"] = global_config.log_level
    if "chunk_size" in kwargs:
        if "snp_batch_size" in kwargs:
            raise ValueError("Pass only one of chunk_size or snp_batch_size.")
        kwargs["snp_batch_size"] = kwargs.pop("chunk_size")
    defaults.update(kwargs)
    args = argparse.Namespace(**defaults)
    return run_build_ref_panel_from_args(args)


def main(argv: Sequence[str] | None = None) -> ReferencePanelBuildResult:
    """Command-line entry point for the reference-panel builder."""

    parser = build_parser()
    args = parser.parse_args(argv)
    return run_build_ref_panel_from_args(args)
