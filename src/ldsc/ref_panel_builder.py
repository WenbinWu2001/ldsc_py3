"""Workflow wrapper for building standard parquet reference panels from PLINK.

Core functionality:
    Resolve public PLINK, optional genetic-map, and optional liftover inputs;
    orchestrate chromosome-wise reference-panel builds; and emit the standard
    parquet artifacts consumed by the refactored LDSC workflows.

Overview
--------
This module is the public entry point for the `build-ref-panel` workflow. It
keeps path resolution, logging, and cross-file validation in the workflow
layer, then delegates pairwise-LD generation and parquet serialization to
``ldsc._kernel.ref_panel_builder``.
"""

from __future__ import annotations

import argparse
import logging
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Sequence

import numpy as np

from .column_inference import normalize_snp_identifier_mode
from .config import GlobalConfig, ReferencePanelBuildConfig, get_global_config, print_global_config_banner
from .path_resolution import ensure_output_directory, resolve_file_group, resolve_plink_prefix_group, resolve_scalar_path
from ._kernel import formats as legacy_parse
from ._kernel import ldscore as kernel_ldscore
from ._kernel import ref_panel_builder as kernel_builder


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
        Paths grouped by artifact kind. ``ann`` and ``ld`` are always present;
        metadata keys such as ``meta_hg19`` and ``meta_hg38`` are present only
        for emitted builds.
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
    restriction_mode: str | None = None
    restriction_values: set[str] | None = None
    translator_cache: dict[tuple[str, str], kernel_builder.LiftOverTranslator] = field(default_factory=dict)


class ReferencePanelBuilder:
    """Build standard parquet reference-panel artifacts from PLINK input."""

    def __init__(self, global_config: GlobalConfig | None = None) -> None:
        """Initialize the builder with shared identifier and build defaults."""
        self.global_config = global_config or get_global_config()

    def run(self, config: ReferencePanelBuildConfig) -> ReferencePanelBuildResult:
        """Build fixed-name parquet reference artifacts for every chromosome.

        Parameters
        ----------
        config : ReferencePanelBuildConfig
            Public build configuration containing one PLINK prefix token,
            optional genetic maps, optional liftover chain paths, output
            directory, and exactly one LD window. If
            ``config.ref_panel_snps_path`` is set, this builder interprets that
            file using ``self.global_config.snp_identifier`` and
            ``self.global_config.genome_build``.

        Returns
        -------
        ReferencePanelBuildResult
            Summary with ``panel_name`` inferred from ``config.output_dir`` and
            artifact paths under ``parquet/ann``, ``parquet/ld``, and
            ``parquet/meta``.

        Raises
        ------
        ValueError
            If no chromosome artifacts are produced or if multiple PLINK inputs
            resolve to the same chromosome.
        """
        print_global_config_banner(type(self).__name__, self.global_config)
        self._configure_logging()
        ensure_output_directory(config.output_dir, label="output directory")
        resolved_prefixes = resolve_plink_prefix_group((config.plink_path,), allow_chromosome_suite=True)
        build_state = self._prepare_build_state(config)

        chrom_records: list[tuple[str, dict[str, str]]] = []
        seen_chromosomes: set[str] = set()
        for prefix in resolved_prefixes:
            for chrom in self._discover_prefix_chromosomes(prefix):
                if chrom in seen_chromosomes:
                    raise ValueError(f"Chromosome {chrom} is present in multiple PLINK inputs. Emit one source per chromosome.")
                seen_chromosomes.add(chrom)
                chrom_paths = self._build_chromosome(prefix, chrom, config, build_state)
                if chrom_paths is None:
                    continue
                chrom_records.append((chrom, chrom_paths))

        if not chrom_records:
            raise ValueError("No chromosome artifacts were produced from the supplied PLINK inputs.")

        chrom_records.sort(key=lambda item: kernel_ldscore.chrom_sort_key(item[0]))
        output_keys = sorted({key for _, paths in chrom_records for key in paths})
        output_paths = {key: [paths[key] for _, paths in chrom_records] for key in output_keys}
        return ReferencePanelBuildResult(
            panel_name=Path(config.output_dir).name,
            chromosomes=[chrom for chrom, _ in chrom_records],
            output_paths=output_paths,
            config_snapshot={
                "build_config": asdict(config),
                "global_config": asdict(self.global_config),
                "resolved_plink_prefixes": list(resolved_prefixes),
            },
        )

    def _configure_logging(self) -> None:
        """Configure package logging for the current build run."""
        logging.basicConfig(
            level=getattr(logging, self.global_config.log_level),
            format="%(levelname)s: %(message)s",
        )

    def _prepare_build_state(self, config: ReferencePanelBuildConfig) -> _BuildState:
        """Resolve optional shared maps, chains, and a reference-universe filter."""
        hg19_files = None
        hg38_files = None
        if config.genetic_map_hg19_path is not None:
            hg19_files = resolve_file_group(
                (config.genetic_map_hg19_path,),
                suffixes=_GENETIC_MAP_SUFFIXES,
                label="hg19 genetic map",
                allow_chromosome_suite=True,
            )
        if config.genetic_map_hg38_path is not None:
            hg38_files = resolve_file_group(
                (config.genetic_map_hg38_path,),
                suffixes=_GENETIC_MAP_SUFFIXES,
                label="hg38 genetic map",
                allow_chromosome_suite=True,
            )
        restriction_mode = None
        restriction_values = None
        if config.ref_panel_snps_path:
            restriction_path = resolve_scalar_path(
                config.ref_panel_snps_path,
                suffixes=_TABLE_SUFFIXES,
                label="reference-panel SNP restriction",
            )
            restriction_mode = normalize_snp_identifier_mode(self.global_config.snp_identifier)
            restriction_values = kernel_ldscore.read_global_snp_restriction(
                restriction_path,
                restriction_mode,
                genome_build=self.global_config.genome_build,
            )
            LOGGER.info(
                "Loaded %d SNP restriction identifiers in %s mode from %s.",
                len(restriction_values),
                restriction_mode,
                restriction_path,
            )
        source_build = config.source_genome_build
        target_build = "hg38" if source_build == "hg19" else "hg19"
        matching_chain = (
            config.liftover_chain_hg19_to_hg38_path
            if source_build == "hg19"
            else config.liftover_chain_hg38_to_hg19_path
        )
        nonmatching_chain = (
            config.liftover_chain_hg38_to_hg19_path
            if source_build == "hg19"
            else config.liftover_chain_hg19_to_hg38_path
        )
        if matching_chain is None:
            if nonmatching_chain is not None:
                LOGGER.warning(
                    "No usable liftover chain was provided for %s -> %s; ignoring the opposite-direction chain and running source-build-only.",
                    source_build,
                    target_build,
                )
            else:
                LOGGER.warning(
                    "No liftover chain was provided for %s -> %s; running source-build-only.",
                    source_build,
                    target_build,
                )
        elif (target_build == "hg38" and config.genetic_map_hg38_path is None) or (
            target_build == "hg19" and config.genetic_map_hg19_path is None
        ):
            LOGGER.warning(
                "No %s genetic map was provided; %s metadata CM values will be written as NA.",
                target_build,
                target_build,
            )
        if (source_build == "hg38" and config.genetic_map_hg38_path is None) or (
            source_build == "hg19" and config.genetic_map_hg19_path is None
        ):
            LOGGER.warning(
                "No %s genetic map was provided; %s metadata CM values will be written as NA.",
                source_build,
                source_build,
            )
        return _BuildState(
            genetic_map_hg19=None if hg19_files is None else kernel_builder.load_genetic_map_group(hg19_files),
            genetic_map_hg38=None if hg38_files is None else kernel_builder.load_genetic_map_group(hg38_files),
            liftover_chain_paths={
                ("hg19", "hg38"): config.liftover_chain_hg19_to_hg38_path,
                ("hg38", "hg19"): config.liftover_chain_hg38_to_hg19_path,
            },
            restriction_mode=restriction_mode,
            restriction_values=restriction_values,
        )

    def _discover_prefix_chromosomes(self, prefix: str) -> list[str]:
        """List the normalized chromosomes present in one PLINK prefix."""
        bim = legacy_parse.PlinkBIMFile(prefix + ".bim")
        chromosomes = sorted(
            {kernel_builder._normalize_map_chromosome(chrom) for chrom in bim.df["CHR"]},
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
        """Build all parquet and metadata artifacts for one chromosome."""
        LOGGER.info("Building reference-panel artifacts for chromosome %s from %s.", chrom, prefix)
        bim = legacy_parse.PlinkBIMFile(prefix + ".bim")
        fam = legacy_parse.PlinkFAMFile(prefix + ".fam")
        panel_df = bim.df.copy()
        panel_df["CHR"] = panel_df["CHR"].map(kernel_builder._normalize_map_chromosome)
        chrom_df = panel_df.loc[panel_df["CHR"] == chrom].copy()
        if len(chrom_df) == 0:
            raise ValueError(f"No SNPs found for chromosome {chrom} in {prefix}.")

        chrom_metadata = chrom_df.loc[:, ["CHR", "SNP", "BP"]].copy().rename(columns={"BP": "POS"})
        chrom_metadata["POS"] = chrom_metadata["POS"].astype(int)
        keep_snps = chrom_df.index.to_numpy(dtype=int)
        if build_state.restriction_values is not None and build_state.restriction_mode is not None:
            keep_mask = kernel_builder.build_restriction_mask(
                chrom_metadata.reset_index(drop=True),
                build_state.restriction_values,
                build_state.restriction_mode,
            )
            keep_snps = chrom_df.index[keep_mask].to_numpy(dtype=int)
            if len(keep_snps) == 0:
                LOGGER.info("Skipping chromosome %s because no SNPs remain after restriction.", chrom)
                return None

        keep_snps, hg19_lookup, hg38_lookup = self._resolve_mappable_snp_positions(
            build_state=build_state,
            chrom=chrom,
            source_build=config.source_genome_build,
            chrom_df=chrom_df,
            keep_snps=keep_snps,
        )
        if len(keep_snps) == 0:
            LOGGER.info("Skipping chromosome %s because no SNPs remain after liftover filtering.", chrom)
            return None
        keep_snps = _sort_retained_snps_by_source_position(
            keep_snps,
            source_build=config.source_genome_build,
            hg19_lookup=hg19_lookup,
            hg38_lookup=hg38_lookup,
        )

        keep_indivs = None
        if config.keep_indivs_path is not None:
            keep_path = resolve_scalar_path(config.keep_indivs_path, suffixes=_TABLE_SUFFIXES, label="individual keep file")
            keep_indivs = fam.loj(legacy_parse.FilterFile(keep_path).IDList)

        geno = kernel_ldscore.PlinkBEDFile(
            prefix + ".bed",
            len(fam.IDList),
            bim,
            keep_snps=list(keep_snps),
            keep_indivs=None if keep_indivs is None else list(keep_indivs),
            mafMin=config.maf_min,
        )
        metadata = kernel_builder.build_plink_metadata_frame(
            bim=bim,
            kept_snps=geno.kept_snps,
            maf_values=geno.maf,
        )
        if len(metadata) == 0:
            LOGGER.info("Skipping chromosome %s because no SNPs remain after PLINK filtering.", chrom)
            return None
        if set(metadata["CHR"]) != {chrom}:
            raise ValueError(f"PLINK filtering for chromosome {chrom} retained rows from multiple chromosomes.")
        geno_kept_snps = list(geno.kept_snps)
        hg19_positions = self._positions_from_lookup(geno_kept_snps, hg19_lookup) if hg19_lookup else None
        hg38_positions = self._positions_from_lookup(geno_kept_snps, hg38_lookup) if hg38_lookup else None

        cm_hg19 = (
            kernel_builder.interpolate_genetic_map_cm(chrom, hg19_positions, build_state.genetic_map_hg19)
            if hg19_positions is not None and build_state.genetic_map_hg19 is not None
            else None
        )
        cm_hg38 = (
            kernel_builder.interpolate_genetic_map_cm(chrom, hg38_positions, build_state.genetic_map_hg38)
            if hg38_positions is not None and build_state.genetic_map_hg38 is not None
            else None
        )
        source_cm = cm_hg19 if config.source_genome_build == "hg19" else cm_hg38
        if source_cm is None and config.ld_wind_cm is not None:
            raise ValueError(f"{config.source_genome_build} genetic map is required for chromosome {chrom}.")
        coords, max_dist = kernel_builder.build_window_coordinates(
            metadata=metadata,
            cm_values=source_cm,
            ld_wind_snps=config.ld_wind_snps,
            ld_wind_kb=config.ld_wind_kb,
            ld_wind_cm=config.ld_wind_cm,
        )
        block_left = kernel_builder.compute_block_left(coords, max_dist)

        annotation_table = kernel_builder.build_standard_annotation_table(
            metadata=metadata,
            hg19_positions=hg19_positions,
            hg38_positions=hg38_positions,
        )
        meta_hg19 = (
            kernel_builder.build_runtime_metadata_table(
                metadata=metadata,
                positions=hg19_positions,
                cm_values=cm_hg19,
            )
            if hg19_positions is not None
            else None
        )
        meta_hg38 = (
            kernel_builder.build_runtime_metadata_table(
                metadata=metadata,
                positions=hg38_positions,
                cm_values=cm_hg38,
            )
            if hg38_positions is not None
            else None
        )

        out_root = Path(config.output_dir) / "parquet"
        ann_path = out_root / "ann" / f"chr{chrom}_ann.parquet"
        ld_path = out_root / "ld" / f"chr{chrom}_LD.parquet"
        meta_hg19_path = out_root / "meta" / f"chr{chrom}_meta_hg19.tsv.gz"
        meta_hg38_path = out_root / "meta" / f"chr{chrom}_meta_hg38.tsv.gz"

        kernel_builder.write_dataframe_to_parquet(annotation_table, ann_path)
        kernel_builder.write_ld_parquet(
            pair_rows=kernel_builder.yield_pairwise_r2_rows(
                block_left=block_left,
                chunk_size=config.chunk_size,
                standardized_snp_getter=geno.nextSNPs,
                m=geno.m,
                n=geno.n,
            ),
            annotation_table=annotation_table,
            path=ld_path,
            genome_build=config.source_genome_build,
        )
        output_paths = {
            "ann": str(ann_path),
            "ld": str(ld_path),
        }
        if meta_hg19 is not None:
            kernel_builder.write_runtime_metadata_sidecar(meta_hg19, meta_hg19_path)
            output_paths["meta_hg19"] = str(meta_hg19_path)
        if meta_hg38 is not None:
            kernel_builder.write_runtime_metadata_sidecar(meta_hg38, meta_hg38_path)
            output_paths["meta_hg38"] = str(meta_hg38_path)
        LOGGER.info("Finished chromosome %s with %d retained SNPs.", chrom, len(metadata))
        return output_paths

    def _resolve_mappable_snp_positions(
        self,
        *,
        build_state: _BuildState,
        chrom: str,
        source_build: str,
        chrom_df,
        keep_snps,
    ) -> tuple[Any, dict[int, int], dict[int, int]]:
        """Resolve retained SNP positions, applying liftover only when configured."""
        keep_snps = list(keep_snps)
        candidate_positions = chrom_df.loc[keep_snps, "BP"].to_numpy(dtype=int)
        if source_build == "hg19":
            if build_state.liftover_chain_paths.get(("hg19", "hg38")) is None:
                retained_snps = np.asarray(keep_snps, dtype=int)
                hg19_lookup = {
                    int(idx): int(pos)
                    for idx, pos in zip(retained_snps, candidate_positions)
                }
                return retained_snps, hg19_lookup, {}
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
                return retained_snps, {}, hg38_lookup
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
        if dropped:
            LOGGER.warning(
                "Dropping %d SNPs on chromosome %s after liftover filtering (%d unmapped, %d cross-chromosome).",
                dropped,
                chrom,
                mapping.unmapped_count,
                mapping.cross_chrom_count,
            )
        return retained_snps, hg19_lookup, hg38_lookup

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
            )
            build_state.translator_cache[key] = translator
        return translator.translate_positions(chrom, positions)

    def _positions_from_lookup(self, kept_snps: Sequence[int], lookup: dict[int, int]):
        """Materialize retained positions in the same order as ``kept_snps``."""
        return np.asarray([lookup[int(idx)] for idx in kept_snps], dtype=int)


def _sort_retained_snps_by_source_position(
    keep_snps,
    *,
    source_build: str,
    hg19_lookup: dict[int, int],
    hg38_lookup: dict[int, int],
) -> np.ndarray:
    """Return retained PLINK indices in source-build genomic position order."""
    # Monotone position order here makes index i → position monotone, which lets
    # yield_pairwise_r2_rows flush cross-chunk pairs with non-decreasing POS_1.
    keep_snps = np.asarray(keep_snps, dtype=int)
    lookup = hg19_lookup if source_build == "hg19" else hg38_lookup
    return np.asarray(
        sorted(keep_snps.tolist(), key=lambda idx: (int(lookup[int(idx)]), int(idx))),
        dtype=int,
    )


def build_parser() -> argparse.ArgumentParser:
    """Build the feature parser for reference-panel generation."""

    parser = argparse.ArgumentParser(
        description="Build standard parquet reference-panel artifacts from a PLINK reference input.",
        allow_abbrev=False,
    )
    parser.add_argument("--plink-path", required=True, help="PLINK prefix token for the reference panel.")
    parser.add_argument(
        "--source-genome-build",
        required=True,
        choices=("hg19", "hg37", "GRCh37", "hg38", "GRCh38"),
        help="Genome build of the PLINK coordinates.",
    )
    parser.add_argument(
        "--genetic-map-hg19-path",
        default=None,
        help="Optional genetic map aligned to hg19; required only for hg19 source-build cM windows.",
    )
    parser.add_argument(
        "--genetic-map-hg38-path",
        default=None,
        help="Optional genetic map aligned to hg38; required only for hg38 source-build cM windows.",
    )
    parser.add_argument(
        "--liftover-chain-hg19-to-hg38-path",
        default=None,
        help="Optional chain file for hg19->hg38 liftover. Enables hg38 outputs for hg19 source builds.",
    )
    parser.add_argument(
        "--liftover-chain-hg38-to-hg19-path",
        default=None,
        help="Optional chain file for hg38->hg19 liftover. Enables hg19 outputs for hg38 source builds.",
    )
    parser.add_argument("--output-dir", required=True, help="Output root directory for emitted parquet artifacts.")
    parser.add_argument("--ld-wind-snps", default=None, type=int, help="LD window size in SNPs.")
    parser.add_argument("--ld-wind-kb", default=None, type=float, help="LD window size in kilobases.")
    parser.add_argument("--ld-wind-cm", default=None, type=float, help="LD window size in centiMorgans.")
    parser.add_argument("--maf", default=None, type=float, help="Optional MAF filter for retained SNPs.")
    parser.add_argument("--ref-panel-snps-path", default=None, help="Optional SNP restriction file defining the retained reference-panel universe.")
    parser.add_argument("--snp-identifier", default=None, choices=("rsid", "chr_pos"), help="SNP identifier mode for --ref-panel-snps-path.")
    parser.add_argument("--keep-indivs-path", default=None, help="Optional individual-keep file.")
    parser.add_argument("--chunk-size", default=50, type=int, help="Chunk size for block processing.")
    parser.add_argument("--log-level", default="INFO", choices=("DEBUG", "INFO", "WARNING", "ERROR"))
    return parser


def config_from_args(args: argparse.Namespace) -> tuple[ReferencePanelBuildConfig, GlobalConfig]:
    """Normalize CLI args into public config objects.

    ``--snp-identifier`` is required only when ``--ref-panel-snps-path`` is
    supplied. The derived ``GlobalConfig`` then carries the restriction-file
    identifier mode and uses ``source_genome_build`` as its genome build.
    """
    if args.ref_panel_snps_path is not None and args.snp_identifier is None:
        raise ValueError("--snp-identifier is required when --ref-panel-snps-path is set.")

    build_config = ReferencePanelBuildConfig(
        plink_path=args.plink_path,
        source_genome_build=args.source_genome_build,
        genetic_map_hg19_path=args.genetic_map_hg19_path,
        genetic_map_hg38_path=args.genetic_map_hg38_path,
        liftover_chain_hg19_to_hg38_path=args.liftover_chain_hg19_to_hg38_path,
        liftover_chain_hg38_to_hg19_path=args.liftover_chain_hg38_to_hg19_path,
        output_dir=args.output_dir,
        ld_wind_snps=args.ld_wind_snps,
        ld_wind_kb=args.ld_wind_kb,
        ld_wind_cm=args.ld_wind_cm,
        maf_min=args.maf,
        ref_panel_snps_path=args.ref_panel_snps_path,
        keep_indivs_path=args.keep_indivs_path,
        chunk_size=args.chunk_size,
    )
    global_config = GlobalConfig(
        snp_identifier=args.snp_identifier or get_global_config().snp_identifier,
        genome_build=build_config.source_genome_build,
        log_level=args.log_level,
    )
    return build_config, global_config


def run_build_ref_panel_from_args(args: argparse.Namespace) -> ReferencePanelBuildResult:
    """Run reference-panel building from parsed CLI arguments."""

    build_config, global_config = config_from_args(args)
    return ReferencePanelBuilder(global_config=global_config).run(build_config)


def run_build_ref_panel(**kwargs: Any) -> ReferencePanelBuildResult:
    """Run the reference-panel builder with CLI-equivalent keyword arguments.

    The wrapper accepts the same modern names as ``ldsc build-ref-panel``. When
    ``ref_panel_snps_path`` is supplied, callers must also pass
    ``snp_identifier`` so the SNP restriction file is interpreted explicitly.
    Without a restriction file, no identifier keyword is required.
    """
    forbidden = sorted({"log_level"} & set(kwargs))
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
            "plink_prefix",
            "genetic_map_hg19",
            "genetic_map_hg38",
            "liftover_chain_hg19_to_hg38",
            "liftover_chain_hg38_to_hg19",
            "keep_indivs",
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
                "--plink-path",
                "plink/panel.@",
                "--source-genome-build",
                "hg19",
                "--genetic-map-hg19-path",
                "hg19.map",
                "--genetic-map-hg38-path",
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
    defaults.update(kwargs)
    args = argparse.Namespace(**defaults)
    return run_build_ref_panel_from_args(args)


def main(argv: Sequence[str] | None = None) -> ReferencePanelBuildResult:
    """Command-line entry point for the reference-panel builder."""

    parser = build_parser()
    args = parser.parse_args(argv)
    return run_build_ref_panel_from_args(args)
