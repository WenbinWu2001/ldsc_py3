"""Public workflow for building standard parquet reference panels from PLINK."""

from __future__ import annotations

import argparse
import logging
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Sequence

import numpy as np

from .config import CommonConfig, ReferencePanelBuildConfig
from .path_resolution import resolve_file_group, resolve_plink_prefix_group, resolve_scalar_path
from ._kernel import formats as legacy_parse
from ._kernel import ldscore as kernel_ldscore
from ._kernel import ref_panel_builder as kernel_builder


LOGGER = logging.getLogger("LDSC.ref_panel_builder")
_GENETIC_MAP_SUFFIXES = ("", ".txt", ".txt.gz", ".tsv", ".tsv.gz", ".csv", ".csv.gz")
_TABLE_SUFFIXES = ("", ".txt", ".txt.gz", ".tsv", ".tsv.gz", ".csv", ".csv.gz")


@dataclass(frozen=True)
class ReferencePanelBuildResult:
    """Summary of files produced by one reference-panel build run."""

    panel_label: str
    chromosomes: list[str]
    output_paths: dict[str, list[str]] = field(default_factory=dict)
    config_snapshot: dict[str, Any] = field(default_factory=dict)


@dataclass
class _BuildState:
    genetic_map_hg19: Any
    genetic_map_hg38: Any
    restriction_mode: str | None = None
    restriction_values: set[str] | None = None
    translator_cache: dict[tuple[str, str], kernel_builder.LiftOverTranslator] = field(default_factory=dict)


class ReferencePanelBuilder:
    """Build standard parquet reference-panel artifacts from PLINK input."""

    def __init__(self, common_config: CommonConfig | None = None) -> None:
        self.common_config = common_config or CommonConfig()

    def run(self, config: ReferencePanelBuildConfig) -> ReferencePanelBuildResult:
        self._configure_logging()
        resolved_prefixes = resolve_plink_prefix_group((config.plink_prefix,), allow_chromosome_suite=True)
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
            panel_label=config.panel_label,
            chromosomes=[chrom for chrom, _ in chrom_records],
            output_paths=output_paths,
            config_snapshot={
                "build_config": asdict(config),
                "common_config": asdict(self.common_config),
                "resolved_plink_prefixes": list(resolved_prefixes),
            },
        )

    def _configure_logging(self) -> None:
        logging.basicConfig(
            level=getattr(logging, self.common_config.log_level),
            format="%(levelname)s: %(message)s",
        )

    def _prepare_build_state(self, config: ReferencePanelBuildConfig) -> _BuildState:
        hg19_files = resolve_file_group(
            (config.genetic_map_hg19_path,),
            suffixes=_GENETIC_MAP_SUFFIXES,
            label="hg19 genetic map",
            allow_chromosome_suite=True,
        )
        hg38_files = resolve_file_group(
            (config.genetic_map_hg38_path,),
            suffixes=_GENETIC_MAP_SUFFIXES,
            label="hg38 genetic map",
            allow_chromosome_suite=True,
        )
        restriction_mode = None
        restriction_values = None
        if self.common_config.global_snp_restriction_path:
            restriction_path = resolve_scalar_path(
                self.common_config.global_snp_restriction_path,
                suffixes=_TABLE_SUFFIXES,
                label="global SNP restriction",
            )
            restriction_mode = kernel_builder.detect_restriction_identifier_mode(restriction_path)
            restriction_values = kernel_ldscore.read_global_snp_restriction(restriction_path, restriction_mode)
            LOGGER.info(
                "Loaded %d SNP restriction identifiers in %s mode from %s.",
                len(restriction_values),
                restriction_mode,
                restriction_path,
            )
        return _BuildState(
            genetic_map_hg19=kernel_builder.load_genetic_map_group(hg19_files),
            genetic_map_hg38=kernel_builder.load_genetic_map_group(hg38_files),
            restriction_mode=restriction_mode,
            restriction_values=restriction_values,
        )

    def _discover_prefix_chromosomes(self, prefix: str) -> list[str]:
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
        hg19_positions = self._positions_from_lookup(geno_kept_snps, hg19_lookup)
        hg38_positions = self._positions_from_lookup(geno_kept_snps, hg38_lookup)

        cm_hg19 = kernel_builder.interpolate_genetic_map_cm(chrom, hg19_positions, build_state.genetic_map_hg19)
        cm_hg38 = kernel_builder.interpolate_genetic_map_cm(chrom, hg38_positions, build_state.genetic_map_hg38)
        source_cm = cm_hg19 if config.source_genome_build == "hg19" else cm_hg38
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
        meta_hg19 = kernel_builder.build_runtime_metadata_table(
            metadata=metadata,
            positions=hg19_positions,
            cm_values=cm_hg19,
        )
        meta_hg38 = kernel_builder.build_runtime_metadata_table(
            metadata=metadata,
            positions=hg38_positions,
            cm_values=cm_hg38,
        )

        out_root = Path(config.output_dir) / "parquet"
        ann_path = out_root / "ann" / f"{config.panel_label}_chr{chrom}_ann.parquet"
        ld_path = out_root / "ld" / f"{config.panel_label}_chr{chrom}_LD.parquet"
        meta_hg19_path = out_root / "meta" / f"{config.panel_label}_chr{chrom}_meta_hg19.tsv.gz"
        meta_hg38_path = out_root / "meta" / f"{config.panel_label}_chr{chrom}_meta_hg38.tsv.gz"

        kernel_builder.write_parquet_table(annotation_table, ann_path)
        kernel_builder.write_standard_ld_parquet(
            pair_rows=kernel_builder.yield_pairwise_r2_rows(
                block_left=block_left,
                chunk_size=config.chunk_size,
                standardized_snp_getter=geno.nextSNPs,
                m=geno.m,
                n=geno.n,
            ),
            annotation_table=annotation_table,
            path=ld_path,
        )
        kernel_builder.write_runtime_metadata_sidecar(meta_hg19, meta_hg19_path)
        kernel_builder.write_runtime_metadata_sidecar(meta_hg38, meta_hg38_path)
        LOGGER.info("Finished chromosome %s with %d retained SNPs.", chrom, len(metadata))
        return {
            "ann": str(ann_path),
            "ld": str(ld_path),
            "meta_hg19": str(meta_hg19_path),
            "meta_hg38": str(meta_hg38_path),
        }

    def _resolve_mappable_snp_positions(
        self,
        *,
        build_state: _BuildState,
        chrom: str,
        source_build: str,
        chrom_df,
        keep_snps,
    ) -> tuple[Any, dict[int, int], dict[int, int]]:
        keep_snps = list(keep_snps)
        candidate_positions = chrom_df.loc[keep_snps, "BP"].to_numpy(dtype=int)
        if source_build == "hg19":
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
            translator = kernel_builder.LiftOverTranslator(source_build=source_build, target_build=target_build)
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
        if source_build == target_build:
            return positions.copy()
        key = (source_build, target_build)
        translator = build_state.translator_cache.get(key)
        if translator is None:
            translator = kernel_builder.LiftOverTranslator(source_build=source_build, target_build=target_build)
            build_state.translator_cache[key] = translator
        return translator.translate_positions(chrom, positions)

    def _positions_from_lookup(self, kept_snps: Sequence[int], lookup: dict[int, int]):
        return np.asarray([lookup[int(idx)] for idx in kept_snps], dtype=int)


def build_parser() -> argparse.ArgumentParser:
    """Build the feature parser for reference-panel generation."""

    parser = argparse.ArgumentParser(
        description="Build standard parquet reference-panel artifacts from a PLINK reference input.",
    )
    parser.add_argument("--bfile", default=None, help="PLINK prefix token for the reference panel.")
    parser.add_argument("--bfile-chr", default=None, help="Chromosome-suite PLINK prefix token; `@` is preferred.")
    parser.add_argument("--panel-label", required=True, help="Label prefix used in emitted standard filenames.")
    parser.add_argument(
        "--source-genome-build",
        required=True,
        choices=("hg19", "hg37", "GRCh37", "hg38", "GRCh38"),
        help="Genome build of the PLINK coordinates.",
    )
    parser.add_argument("--genetic-map-hg19", required=True, help="Genetic map file or token aligned to hg19 coordinates.")
    parser.add_argument("--genetic-map-hg38", required=True, help="Genetic map file or token aligned to hg38 coordinates.")
    parser.add_argument("--out", required=True, help="Output root directory for emitted parquet artifacts.")
    parser.add_argument("--ld-wind-snps", default=None, type=int, help="LD window size in SNPs.")
    parser.add_argument("--ld-wind-kb", default=None, type=float, help="LD window size in kilobases.")
    parser.add_argument("--ld-wind-cm", default=None, type=float, help="LD window size in centiMorgans.")
    parser.add_argument("--maf", default=None, type=float, help="Optional MAF filter for retained SNPs.")
    parser.add_argument("--restrict-snps-path", default=None, help="Optional global SNP restriction file.")
    parser.add_argument("--keep-indivs", default=None, help="Optional individual-keep file.")
    parser.add_argument("--chunk-size", default=50, type=int, help="Chunk size for block processing.")
    parser.add_argument("--log-level", default="INFO", choices=("DEBUG", "INFO", "WARNING", "ERROR"))
    return parser


def config_from_args(args: argparse.Namespace) -> tuple[ReferencePanelBuildConfig, CommonConfig]:
    """Normalize CLI args into public config objects."""

    plink_prefix = args.bfile if args.bfile is not None else args.bfile_chr
    if plink_prefix is None:
        raise ValueError("One of --bfile or --bfile-chr is required.")
    build_config = ReferencePanelBuildConfig(
        panel_label=args.panel_label,
        plink_prefix=plink_prefix,
        source_genome_build=args.source_genome_build,
        genetic_map_hg19_path=args.genetic_map_hg19,
        genetic_map_hg38_path=args.genetic_map_hg38,
        output_dir=args.out,
        ld_wind_snps=args.ld_wind_snps,
        ld_wind_kb=args.ld_wind_kb,
        ld_wind_cm=args.ld_wind_cm,
        maf_min=args.maf,
        restrict_snps_path=args.restrict_snps_path,
        keep_indivs_path=args.keep_indivs,
        chunk_size=args.chunk_size,
    )
    common_config = CommonConfig(
        genome_build=build_config.source_genome_build,
        global_snp_restriction_path=build_config.restrict_snps_path,
        log_level=args.log_level,
    )
    return build_config, common_config


def run_build_ref_panel_from_args(args: argparse.Namespace) -> ReferencePanelBuildResult:
    """Run reference-panel building from parsed CLI arguments."""

    build_config, common_config = config_from_args(args)
    return ReferencePanelBuilder(common_config=common_config).run(build_config)


def run_build_ref_panel(**kwargs: Any) -> ReferencePanelBuildResult:
    """Convenience wrapper around :func:`run_build_ref_panel_from_args`."""

    parser = build_parser()
    defaults = vars(
        parser.parse_args(
            [
                "--panel-label",
                "placeholder",
                "--source-genome-build",
                "hg19",
                "--genetic-map-hg19",
                "hg19.map",
                "--genetic-map-hg38",
                "hg38.map",
                "--out",
                "out",
                "--ld-wind-kb",
                "1",
            ]
        )
    )
    defaults.update(kwargs)
    args = argparse.Namespace(**defaults)
    return run_build_ref_panel_from_args(args)


def main(argv: Sequence[str] | None = None) -> ReferencePanelBuildResult:
    """Command-line entry point for the reference-panel builder."""

    parser = build_parser()
    args = parser.parse_args(argv)
    return run_build_ref_panel_from_args(args)
