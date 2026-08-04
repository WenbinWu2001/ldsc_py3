"""Build, validate, and query exact disjoint-atom gene LD-score indexes.

The offline workflow computes a PLINK-backed LD-score operator for one fixed
baseline suite, reference panel, and gene-projection profile. The online
workflow resolves gene lists against the embedded catalog, assembles their
Boolean union from the stored operator, and writes an ordinary self-contained
LD-score directory. It never discovers an index or silently falls back to live
LD calculation.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
import hashlib
import json
import logging
import os
from pathlib import Path
import socket
import shutil
import sys
import tempfile
import time
from typing import Sequence

import numpy as np
import pandas as pd
from scipy import sparse

from ._kernel import ldscore as kernel_ldscore
from ._kernel.gene_ldscore_index import (
    AtomStatistics,
    ChromosomeAtomModel,
    build_disjoint_atoms,
    compute_atom_statistics,
    iter_snp_atom_blocks,
    map_snps_to_atoms,
    assemble_atom_selector,
    assemble_indexed_ld_scores,
    validate_ldscore_operator,
)
from ._kernel.regions import RegionIntervals
from ._kernel.identifiers import read_snp_restriction_keys
from ._kernel.snp_identity import RestrictionIdentityKeys
from ._kernel import regions as kernel_regions
from ._logging import log_inputs, log_outputs, set_workflow_log_path, workflow_logging
from .annotation_builder import AnnotationBuilder
from .chromosome_inference import normalize_chromosome
from .config import AnnotationBuildConfig, GeneLDScoreIndexBuildConfig, GlobalConfig
from .errors import LDSCInputError, LDSCInternalError
from .gene_list_resolver import GeneCatalog, GeneListResolution, resolve_all_protein_coding, resolve_gene_list
from .path_resolution import resolve_file_group
from .query_annotations import QueryAnnotationStatus
from .hm3 import packaged_hm3_curated_map_path
from .path_resolution import split_cli_path_tokens


LOGGER = logging.getLogger("LDSC.gene_ldscore_index")


_STRICT_IDENTITY_COLUMNS = ["CHR", "POS", "SNP"]


def _semantic_sha256(payload: dict) -> str:
    """Hash one scientific identity with canonical JSON normalization."""
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def calculate_suite_id(suite_identity: dict) -> str:
    """Return the immutable common-suite semantic identity."""
    return _semantic_sha256(suite_identity)


def calculate_profile_id(suite_id: str, profile_identity: dict) -> str:
    """Return the profile identity bound to one exact common suite."""
    return _semantic_sha256({"suite_id": suite_id, **profile_identity})


def build_parser() -> argparse.ArgumentParser:
    """Build the closed-v1 offline gene LD-score index parser."""
    parser = argparse.ArgumentParser(
        prog="ldsc build-gene-ldscore-index",
        description="Build an exact PLINK-backed disjoint-atom gene LD-score index.",
        allow_abbrev=False,
    )
    parser.add_argument("--baseline-annot-sources", required=True)
    parser.add_argument("--plink-prefix", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--genome-build", choices=("hg19",), default="hg19")
    parser.add_argument("--snp-identifier", choices=("rsid",), default="rsid")
    parser.add_argument("--padding-bp", type=int, default=100000)
    parser.add_argument("--gene-exclude-regions", choices=("none", "mhc"), default="mhc")
    parser.add_argument("--ld-wind-cm", type=float, default=1.0)
    parser.add_argument("--maf-min", type=float, default=None)
    parser.add_argument("--common-maf-min", type=float, default=0.05)
    parser.add_argument("--keep-indivs-file", default=None)
    parser.add_argument("--genetic-map-hg19-sources", default=None)
    parser.add_argument("--genetic-map-hg38-sources", default=None)
    parser.add_argument("--chromosomes", default="1-22")
    parser.add_argument("--snp-batch-size", type=int, default=128)
    parser.add_argument("--atom-batch-size", type=int, default=64)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--overwrite", action="store_true", default=False)
    parser.add_argument("--log-level", choices=("DEBUG", "INFO", "WARNING", "ERROR"), default="INFO")
    return parser


def main(argv: Sequence[str] | None = None):
    """Run the offline gene LD-score index workflow."""
    args = build_parser().parse_args(argv)
    return run_build_gene_ldscore_index_from_args(args)


def build_gene_ldscore_index(
    config: GeneLDScoreIndexBuildConfig,
    *,
    genetic_map_hg19_sources: str | None = None,
    genetic_map_hg38_sources: str | None = None,
    overwrite: bool = False,
    log_level: str = "INFO",
) -> Path:
    """Build and transactionally publish one exact gene LD-score profile.

    Parameters
    ----------
    config : GeneLDScoreIndexBuildConfig
        Validated baseline, PLINK, chromosome, projection, window, filtering,
        and batching settings. V1 supports hg19, rsID identity, and a 1 cM
        PLINK-backed workflow.
    genetic_map_hg19_sources : str, optional
        Comma-separated hg19 genetic-map paths. When omitted, informative cM
        coordinates from each PLINK BIM are used.
    genetic_map_hg38_sources : str, optional
        Reserved matching-build argument. V1 rejects hg38 maps because its
        projection build is hg19.
    overwrite : bool, optional
        Replace only the targeted profile when it already exists. Matching
        common data and sibling profiles are preserved. Default is ``False``.
    log_level : {"DEBUG", "INFO", "WARNING", "ERROR"}, optional
        Workflow logging threshold. Default is ``"INFO"``.

    Returns
    -------
    pathlib.Path
        Published ``profiles/<profile-name>`` directory.

    Raises
    ------
    TypeError
        If ``config`` is not a ``GeneLDScoreIndexBuildConfig``.
    LDSCInputError
        If an input is missing or incompatible, baseline/BIM identities do not
        match before genotype QC, or the staged artifact fails validation.
    FileExistsError
        If the targeted profile exists and ``overwrite`` is false.

    Notes
    -----
    Construction uses float64 adjusted-r-squared accumulation, retains
    negative values, and batches internal atom columns. ``threads`` controls
    chromosome workers and can multiply chromosome-local peak memory.
    """
    if not isinstance(config, GeneLDScoreIndexBuildConfig):
        raise TypeError("build_gene_ldscore_index requires GeneLDScoreIndexBuildConfig.")
    return run_build_gene_ldscore_index_from_args(
        argparse.Namespace(
            baseline_annot_sources=",".join(config.baseline_annot_sources),
            plink_prefix=config.plink_prefix,
            output_dir=config.output_dir,
            genome_build=config.genome_build,
            snp_identifier=config.snp_identifier,
            padding_bp=config.padding_bp,
            gene_exclude_regions=config.gene_exclude_regions,
            ld_wind_cm=config.ld_wind_cm,
            maf_min=config.maf_min,
            common_maf_min=config.common_maf_min,
            keep_indivs_file=config.keep_indivs_file,
            genetic_map_hg19_sources=genetic_map_hg19_sources,
            genetic_map_hg38_sources=genetic_map_hg38_sources,
            chromosomes=",".join(config.chromosomes),
            snp_batch_size=config.snp_batch_size,
            atom_batch_size=config.atom_batch_size,
            threads=config.threads,
            overwrite=overwrite,
            log_level=log_level,
        )
    )


def _parse_chromosomes(value: str | Sequence[str]) -> tuple[str, ...]:
    """Parse a canonical autosome subset such as ``1-3,22``."""
    tokens = [str(item).strip() for item in (value if not isinstance(value, str) else value.split(","))]
    selected: set[int] = set()
    for token in tokens:
        if not token:
            continue
        if "-" in token:
            parts = token.split("-", 1)
            try:
                start, end = (int(part) for part in parts)
            except ValueError as exc:
                raise LDSCInputError(f"Invalid chromosome range {token!r}.") from exc
            selected.update(range(start, end + 1))
        else:
            try:
                selected.add(int(token))
            except ValueError as exc:
                raise LDSCInputError(f"Invalid chromosome token {token!r}.") from exc
    if not selected or min(selected) < 1 or max(selected) > 22:
        raise LDSCInputError("Gene LD-score index chromosomes must be autosomes 1 through 22.")
    return tuple(str(chrom) for chrom in sorted(selected))


def _build_embedded_gene_catalog(
    catalog: GeneCatalog,
    *,
    genome_build: str,
    gene_exclude_regions: str,
) -> pd.DataFrame:
    """Build the profile's self-contained ordered catalog and inclusion policy."""
    eligible = set(
        resolve_all_protein_coding(
            catalog,
            genome_build=genome_build,
            gene_exclude_regions=gene_exclude_regions,
        ).catalog_indices
    )
    rows: list[dict] = []
    next_row: dict[str, int] = {}
    for gene_index, row in catalog.frame.iterrows():
        chrom_value = row[f"{genome_build}_chr"]
        has_coordinates = not pd.isna(chrom_value)
        chrom = normalize_chromosome(chrom_value) if has_coordinates else None
        chromosome_gene_row = -1
        if chrom is not None:
            chromosome_gene_row = next_row.get(chrom, 0)
            next_row[chrom] = chromosome_gene_row + 1
        included = int(gene_index) in eligible
        rows.append(
            {
                "gene_index": int(gene_index),
                "canonical_ensembl_id": str(row["ensgid"]),
                "gene_name": str(row["gene_name"]),
                "CHR": chrom,
                "start0": None if not has_coordinates else int(row[f"{genome_build}_start0"]),
                "end": None if not has_coordinates else int(row[f"{genome_build}_end"]),
                "included": bool(included),
                "exclusion_reason": "" if included else ("excluded_gene_region" if has_coordinates else "build_missing"),
                "chromosome_gene_row": chromosome_gene_row,
            }
        )
    return pd.DataFrame(rows)


def run_build_gene_ldscore_index_from_args(args: argparse.Namespace) -> Path:
    """Build, validate, and publish one exact v1 index profile."""
    started = time.perf_counter()
    chromosomes = _parse_chromosomes(args.chromosomes)
    baseline_sources = tuple(split_cli_path_tokens(args.baseline_annot_sources))
    config = GeneLDScoreIndexBuildConfig(
        baseline_annot_sources=baseline_sources,
        plink_prefix=args.plink_prefix,
        output_dir=args.output_dir,
        chromosomes=chromosomes,
        genome_build=args.genome_build,
        snp_identifier=args.snp_identifier,
        padding_bp=args.padding_bp,
        gene_exclude_regions=args.gene_exclude_regions,
        ld_wind_cm=args.ld_wind_cm,
        maf_min=args.maf_min,
        common_maf_min=args.common_maf_min,
        keep_indivs_file=args.keep_indivs_file,
        snp_batch_size=args.snp_batch_size,
        atom_batch_size=args.atom_batch_size,
        threads=args.threads,
    )

    profile_name = f"padding-{config.padding_bp}bp-{config.gene_exclude_regions}"
    final_log_path = (
        Path(config.output_dir)
        / "profiles"
        / profile_name
        / "diagnostics"
        / "build-gene-ldscore-index.log"
    )
    Path(config.output_dir).parent.mkdir(parents=True, exist_ok=True)
    log_stage_dir = Path(
        tempfile.mkdtemp(prefix=f".{Path(config.output_dir).name}.gene-index-log-", dir=Path(config.output_dir).parent)
    )
    staged_log_path = log_stage_dir / final_log_path.name
    try:
        with workflow_logging("build-gene-ldscore-index", staged_log_path, log_level=args.log_level):
            log_inputs(
                output_dir=str(config.output_dir),
                genome_build=config.genome_build,
                snp_identifier=config.snp_identifier,
                chromosomes=", ".join(chromosomes),
                baseline_annot_sources=", ".join(config.baseline_annot_sources),
                plink_prefix=config.plink_prefix,
                ld_wind_cm=config.ld_wind_cm,
                padding_bp=config.padding_bp,
                gene_exclude_regions=config.gene_exclude_regions,
                maf_min=config.maf_min,
                common_maf_min=config.common_maf_min,
                keep_indivs_file=config.keep_indivs_file or "all individuals",
                genetic_map="explicit hg19 map" if getattr(args, "genetic_map_hg19_sources", None) else "BIM cM fallback",
                regression_policy="bundled HM3 minus MHC and centromeres",
                snp_batch_size=config.snp_batch_size,
                atom_batch_size=config.atom_batch_size,
                threads=config.threads,
                effective_log_level=args.log_level,
                process_id=os.getpid(),
                hostname=socket.gethostname(),
            )
            return _run_gene_ldscore_index_build(
                args,
                config,
                chromosomes,
                started=started,
                profile_name=profile_name,
            )
    finally:
        _materialize_gene_index_log(staged_log_path, final_log_path, log_stage_dir)


def _run_gene_ldscore_index_build(
    args: argparse.Namespace,
    config: GeneLDScoreIndexBuildConfig,
    chromosomes: tuple[str, ...],
    *,
    started: float,
    profile_name: str,
) -> Path:
    """Run the scientific builder inside its shared workflow log context."""
    catalog = GeneCatalog.load()
    embedded_catalog = _build_embedded_gene_catalog(
        catalog,
        genome_build=config.genome_build,
        gene_exclude_regions=config.gene_exclude_regions,
    )
    global_config = GlobalConfig(snp_identifier="rsid")
    annotation_spec = AnnotationBuildConfig(baseline_annot_sources=config.baseline_annot_sources)
    regression_path = Path(packaged_hm3_curated_map_path())
    regression_keys = read_snp_restriction_keys(regression_path, "rsid", genome_build=None)
    regression_regions = kernel_regions.load_preset_intervals(("mhc", "centromeres"), "hg19")
    genetic_map = _load_builder_genetic_map(args)

    def build_one(chrom: str) -> tuple[str, IndexChromosomeData, dict]:
        chrom_started = time.perf_counter()
        LOGGER.info(f"Starting chromosome {chrom}.")
        phase = "baseline annotation resolution"
        try:
            public_bundle = AnnotationBuilder(global_config).run(annotation_spec, chrom=chrom)
            chrom_mask = public_bundle.metadata["CHR"].astype(str) == chrom
            metadata = public_bundle.metadata.loc[chrom_mask].reset_index(drop=True)
            baseline = public_bundle.baseline_annotations.loc[chrom_mask].reset_index(drop=True)
            if metadata.empty:
                raise LDSCInputError(f"No baseline annotation rows were found for chromosome {chrom}.")
            kernel_bundle = kernel_ldscore.AnnotationBundle(
                metadata=metadata,
                annotations=baseline,
                baseline_columns=list(public_bundle.baseline_columns),
                query_columns=[],
            )
            kernel_args = argparse.Namespace(
                bfile=config.plink_prefix,
                keep=config.keep_indivs_file,
                maf_min=config.maf_min,
                maf=None,
                snp_identifier="rsid",
                ld_wind_snps=None,
                ld_wind_kb=None,
                ld_wind_cm=config.ld_wind_cm,
                yes_really=False,
                snp_batch_size=config.snp_batch_size,
                common_maf_min=config.common_maf_min,
                genetic_map=genetic_map,
            )
            prefix = kernel_ldscore.resolve_bfile_prefix(kernel_args, chrom=chrom)
            if prefix is None:
                raise LDSCInputError(f"Could not resolve PLINK prefix for chromosome {chrom}.")
            bim_rows = _read_bim_identity(prefix)
            phase = "strict baseline/BIM validation"
            validate_strict_baseline_plink_identity(metadata, bim_rows, chrom=chrom)
            chrom_catalog = embedded_catalog.loc[embedded_catalog["CHR"].astype(str) == chrom].sort_values(
                "chromosome_gene_row", kind="mergesort"
            )
            intervals = chrom_catalog[["start0", "end"]].to_numpy(dtype=np.int64)
            included = chrom_catalog["included"].to_numpy(dtype=bool)
            phase = "genotype QC and atomic LD-score construction"
            record = build_plink_index_chromosome(
                chrom,
                kernel_bundle,
                kernel_args,
                regression_keys=regression_keys,
                regression_regions=regression_regions,
                gene_intervals=intervals,
                included=included,
                padding_bp=config.padding_bp,
                atom_batch_size=config.atom_batch_size,
            )
        except Exception as exc:
            LOGGER.error(
                f"Chromosome {chrom} failed during {phase}: {type(exc).__name__}: {exc}"
            )
            raise
        protein_coding_genes = int(included.sum())
        genes_with_atoms = int(np.count_nonzero(np.diff(record.atom_model.gene_to_atom.indptr)))
        rows_before_genotype_qc = int(record.reference_rows_before_genotype_qc or record.total_reference_snps_all)
        genotype_qc_removed = int(record.genotype_qc_removed)
        maf_removed = int(record.maf_removed)
        evidence = {
            "pre_qc_rows": int(len(bim_rows)),
            "reference_rows_before_genotype_qc": rows_before_genotype_qc,
            "annotation_intersection_removed": int(len(bim_rows) - rows_before_genotype_qc),
            "retained_reference_rows": int(record.total_reference_snps_all),
            "retained_common_reference_rows": int(record.total_reference_snps_common),
            "genotype_qc_removed": genotype_qc_removed,
            "maf_removed": maf_removed,
            "genotype_qc_or_maf_removed": genotype_qc_removed + maf_removed,
            "regression_rows": int(len(record.baseline_rows)),
            "baseline_bim_identity": "exact",
            "protein_coding_genes": protein_coding_genes,
            "genes_with_padded_atoms": genes_with_atoms,
            "atom_count": int(record.atom_model.n_atoms),
            "operator_nnz": int(record.operator.nnz),
            "nnz_Y": int(record.operator.nnz),
            "snp_batch_size": int(config.snp_batch_size),
            "atom_batch_size": int(config.atom_batch_size),
            "maf_filter_policy": "disabled" if config.maf_min is None else f">={config.maf_min}",
            "cm_source": record.cm_source,
            "elapsed_seconds": time.perf_counter() - chrom_started,
        }
        LOGGER.info(
            f"Finished chromosome {chrom}: protein-coding genes={protein_coding_genes}, "
            f"retained-reference={evidence['retained_reference_rows']}, "
            f"regression-rows={evidence['regression_rows']}, atoms={evidence['atom_count']}, "
            f"nnz(Y)={evidence['nnz_Y']}, elapsed={evidence['elapsed_seconds']:.3f}s."
        )
        return chrom, record, evidence

    if config.threads == 1 or len(chromosomes) == 1:
        built = [build_one(chrom) for chrom in chromosomes]
    else:
        max_workers = os.cpu_count() if config.threads == -1 else config.threads
        if max_workers is None or max_workers < 1:
            max_workers = max(1, (os.cpu_count() or 1) + 1 + config.threads)
        with ThreadPoolExecutor(max_workers=min(max_workers, len(chromosomes))) as pool:
            built = list(pool.map(build_one, chromosomes))
    chromosome_data = {chrom: record for chrom, record, _evidence in built}
    evidence_by_chrom = {chrom: evidence for chrom, _record, evidence in built}
    suite_identity = _builder_suite_identity(config, args, chromosome_data, regression_path)
    profile_identity = {
        "catalog": {
            "resource": catalog.resource,
            "release": catalog.release,
            "content_sha256": catalog.content_sha256,
        },
        "projection_build": config.genome_build,
        "padding_bp": config.padding_bp,
        "gene_exclude_regions": config.gene_exclude_regions,
    }
    suite_path = Path(config.output_dir)
    target_profile_path = suite_path / "profiles" / profile_name
    suite_existed = suite_path.exists()
    common_existed = (suite_path / "common").exists()
    target_profile_existed = target_profile_path.exists()
    sibling_profiles_before = _profile_names(suite_path, exclude=profile_name)
    suite_id = calculate_suite_id(suite_identity)
    profile_id = calculate_profile_id(suite_id, profile_identity)
    try:
        profile_path = publish_gene_ldscore_index(
            config.output_dir,
            profile_name=profile_name,
            suite_identity=suite_identity,
            profile_identity=profile_identity,
            gene_catalog=embedded_catalog,
            chromosomes=chromosome_data,
            overwrite=bool(args.overwrite),
        )
    except Exception as exc:
        LOGGER.error(
            f"Gene LD-score index build failed during staged validation/publication: "
            f"{type(exc).__name__}: {exc}"
        )
        raise
    diagnostics = profile_path / "diagnostics"
    diagnostics.mkdir(exist_ok=True)
    elapsed = time.perf_counter() - started
    import resource

    peak_rss = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    if sys.platform != "darwin":
        peak_rss *= 1024
    payload_bytes = _scientific_payload_bytes(suite_path)
    for chrom, values in evidence_by_chrom.items():
        values["common_payload_bytes"] = _directory_payload_bytes(suite_path / "common" / f"chr{chrom}")
        values["profile_payload_bytes"] = _directory_payload_bytes(profile_path / f"chr{chrom}")
    selected_individuals = suite_identity.get("selected_individuals", {})
    plink_sources = suite_identity.get("plink_sources", [])
    map_identity = suite_identity.get("genetic_map", "bim_cm")
    if isinstance(map_identity, dict):
        map_source = "explicit hg19 map"
        map_file_count = len(map_identity.get("sources", []))
    else:
        map_source = "BIM cM fallback"
        map_file_count = 0
    publication = {
        "suite_existed": suite_existed,
        "common": "reused" if common_existed else "created",
        "target_profile_existed": target_profile_existed,
        "targeted_overwrite": bool(args.overwrite),
        "staging_reload_validation": "passed",
        "atomic_publication": "completed",
        "sibling_profiles_preserved": sibling_profiles_before,
        "profile_path": str(profile_path),
    }
    diagnostic_payload = {
        "suite_id": suite_id,
        "profile_id": profile_id,
        "elapsed_seconds": elapsed,
        "peak_rss_bytes": peak_rss,
        "payload_bytes": payload_bytes,
        "snp_batch_size": config.snp_batch_size,
        "atom_batch_size": config.atom_batch_size,
        "threads": config.threads,
        "effective_log_level": args.log_level,
        "configuration": {
            "genome_build": config.genome_build,
            "snp_identifier": config.snp_identifier,
            "chromosomes": list(chromosomes),
            "baseline_annot_sources": list(config.baseline_annot_sources),
            "plink_prefix": config.plink_prefix,
            "ld_window": {"unit": "cm", "value": config.ld_wind_cm},
            "padding_bp": config.padding_bp,
            "gene_exclude_regions": config.gene_exclude_regions,
            "maf_min": config.maf_min,
            "common_maf_min": config.common_maf_min,
            "keep_individuals": "all" if config.keep_indivs_file is None else "keep_file",
            "genetic_map": map_source,
            "regression_policy": "bundled HM3 minus MHC and centromeres",
        },
        "input_resolution": {
            "baseline_annotation_files": len(suite_identity.get("baseline_sources", [])),
            "plink_bed_files": sum(source.get("kind") == "bed" for source in plink_sources),
            "plink_bim_files": sum(source.get("kind") == "bim" for source in plink_sources),
            "plink_fam_files": sum(source.get("kind") == "fam" for source in plink_sources),
            "chromosome_shards": len(chromosomes),
            "genetic_map_files": map_file_count,
            "selected_individual_count": selected_individuals.get("selected_count"),
            "keep_individual_source": selected_individuals.get("source", "all"),
            "strict_baseline_bim_identity": {
                chrom: values["baseline_bim_identity"] for chrom, values in evidence_by_chrom.items()
            },
            "map_validation": "passed and informative",
        },
        "publication": publication,
        "chromosomes": evidence_by_chrom,
        "totals": {
            "chromosome_count": len(chromosomes),
            "protein_coding_genes": sum(values["protein_coding_genes"] for values in evidence_by_chrom.values()),
            "atom_count": sum(values["atom_count"] for values in evidence_by_chrom.values()),
            "operator_nnz": sum(values["operator_nnz"] for values in evidence_by_chrom.values()),
            "regression_rows": sum(values["regression_rows"] for values in evidence_by_chrom.values()),
            "retained_reference_rows": sum(values["retained_reference_rows"] for values in evidence_by_chrom.values()),
        },
    }
    _write_json(
        diagnostics / "build-gene-ldscore-index.json",
        diagnostic_payload,
    )
    _log_gene_index_summary(diagnostic_payload)
    log_outputs(
        suite_id=suite_id,
        profile_id=profile_id,
        profile=str(profile_path),
        summary_json=str(diagnostics / "build-gene-ldscore-index.json"),
        publication="staged, reloaded, validated, and atomically replaced",
    )
    return profile_path


def _materialize_gene_index_log(staged_log_path: Path, final_log_path: Path, stage_dir: Path) -> None:
    """Copy the live workflow log into the published profile or failure location."""
    try:
        if staged_log_path.exists():
            final_log_path.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(staged_log_path, final_log_path)
            set_workflow_log_path(final_log_path)
    finally:
        if stage_dir.exists():
            shutil.rmtree(stage_dir, ignore_errors=True)


def _profile_names(suite_path: Path, *, exclude: str | None = None) -> list[str]:
    """Return existing sibling profile names without reading profile payloads."""
    profiles = suite_path / "profiles"
    if not profiles.is_dir():
        return []
    return sorted(
        path.name for path in profiles.iterdir() if path.is_dir() and path.name != exclude
    )


def _directory_payload_bytes(path: Path) -> int:
    """Return the byte size of one scientific artifact component."""
    if not path.exists():
        return 0
    return sum(file.stat().st_size for file in path.rglob("*") if file.is_file())


def _scientific_payload_bytes(suite_path: Path) -> int:
    """Return suite bytes excluding human-readable and JSON diagnostics."""
    total = 0
    for file in suite_path.rglob("*"):
        if not file.is_file() or "diagnostics" in file.parts:
            continue
        total += file.stat().st_size
    return total


def _is_abandoned_gene_index_profile(path: Path) -> bool:
    """Recognize a diagnostics-only profile left by a failed build."""
    if not path.is_dir() or (path / "metadata.json").exists():
        return False
    files = [file for file in path.rglob("*") if file.is_file()]
    return bool(files) and all(
        file.relative_to(path).parts == ("diagnostics", "build-gene-ldscore-index.log")
        for file in files
    )


def _is_abandoned_gene_index_suite(path: Path) -> bool:
    """Recognize a diagnostics-only suite left by a failed first build."""
    if not path.is_dir() or (path / "metadata.json").exists():
        return False
    files = [file for file in path.rglob("*") if file.is_file()]
    if not files:
        return False
    for file in files:
        relative = file.relative_to(path)
        if (
            len(relative.parts) != 4
            or relative.parts[0] != "profiles"
            or relative.parts[2] != "diagnostics"
            or relative.parts[3] != "build-gene-ldscore-index.log"
        ):
            return False
    return True


def _log_gene_index_summary(payload: dict) -> None:
    """Render the shared diagnostic payload as a concise operational narrative."""
    config = payload["configuration"]
    inputs = payload["input_resolution"]
    publication = payload["publication"]
    totals = payload["totals"]
    LOGGER.info("Gene LD-score index build configuration resolved.")
    LOGGER.info(
        f"Build settings: build={config['genome_build']}, identity={config['snp_identifier']}, "
        f"chromosomes={','.join(config['chromosomes'])}, ld-window={config['ld_window']['value']} "
        f"{config['ld_window']['unit']}, padding={config['padding_bp']} bp, "
        f"gene-exclusion={config['gene_exclude_regions']}, maf-min={config['maf_min']}, "
        f"common-maf-min={config['common_maf_min']}, keep-individuals={config['keep_individuals']}, "
        f"map={config['genetic_map']}."
    )
    LOGGER.info(
        f"Resolved inputs: baseline-files={inputs['baseline_annotation_files']}, "
        f"PLINK bed/bim/fam={inputs['plink_bed_files']}/{inputs['plink_bim_files']}/"
        f"{inputs['plink_fam_files']}, map-files={inputs['genetic_map_files']}, "
        f"chromosome-shards={inputs['chromosome_shards']}, "
        f"selected-individuals={inputs['selected_individual_count']} "
        f"({inputs['keep_individual_source']})."
    )
    LOGGER.info(
        "SNP universe: broad retained PLINK SNPs are LD-score contributors and count/overlap members; "
        "filtered HM3 SNPs are persisted regression rows; w_ld uses the filtered regression set as rows and contributors."
    )
    LOGGER.info("Regression-row policy: bundled HM3 minus MHC and centromeres.")
    for chrom in payload["chromosomes"]:
        values = payload["chromosomes"][chrom]
        LOGGER.info(
            f"Chromosome {chrom} metrics: pre-QC={values['pre_qc_rows']}, "
            f"strict-baseline-BIM={values['baseline_bim_identity']}, "
            f"annotation-intersection-removed={values['annotation_intersection_removed']}, "
            f"genotype-QC-removed={values['genotype_qc_removed']}, "
            f"MAF-removed={values['maf_removed']} (policy={values['maf_filter_policy']}), "
            f"retained-reference={values['retained_reference_rows']}, "
            f"common-reference={values['retained_common_reference_rows']}, "
            f"regression-rows={values['regression_rows']}, "
            f"protein-coding genes={values['protein_coding_genes']}, "
            f"genes-with-padded-atoms={values['genes_with_padded_atoms']}, "
            f"atoms={values['atom_count']}, nnz(Y)={values['nnz_Y']}, "
            f"operator_nnz={values['operator_nnz']}, "
            f"common-bytes={values['common_payload_bytes']}, "
            f"profile-bytes={values['profile_payload_bytes']}, "
            f"elapsed={values['elapsed_seconds']:.3f}s."
        )
    LOGGER.info(
        f"Publication: suite_id={payload['suite_id']}, profile_id={payload['profile_id']}, "
        f"common={publication['common']}, target-overwrite={publication['targeted_overwrite']}, "
        f"sibling-profiles-preserved={len(publication['sibling_profiles_preserved'])}, "
        f"staged-reload={publication['staging_reload_validation']}, "
        f"atomic-publication={publication['atomic_publication']}."
    )
    LOGGER.info(
        f"Gene LD-score index build completed and validated: chromosomes={totals['chromosome_count']}, "
        f"protein-coding genes={totals['protein_coding_genes']}, atoms={totals['atom_count']}, "
        f"operator_nnz={totals['operator_nnz']}, payload_bytes={payload['payload_bytes']}, "
        f"peak_rss_bytes={payload['peak_rss_bytes']}, elapsed_seconds={payload['elapsed_seconds']:.6f}, "
        f"profile={publication['profile_path']}."
    )


def _read_bim_identity(prefix: str) -> pd.DataFrame:
    frame = pd.read_csv(
        prefix + ".bim",
        sep=r"\s+",
        header=None,
        names=["CHR", "SNP", "CM", "POS", "A1", "A2"],
    )
    return frame[["CHR", "POS", "SNP"]]


def _file_sha256(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _load_builder_genetic_map(args):
    hg19 = getattr(args, "genetic_map_hg19_sources", None)
    hg38 = getattr(args, "genetic_map_hg38_sources", None)
    if hg38:
        raise LDSCInputError("The hg19 v1 index builder cannot use --genetic-map-hg38-sources.")
    if not hg19:
        return None
    from ._kernel.ref_panel_builder import load_genetic_map_group

    return load_genetic_map_group(split_cli_path_tokens(hg19))


def _genetic_map_identity(args) -> str | dict:
    """Return the content-bound scientific identity of the effective cM map."""
    sources = split_cli_path_tokens(getattr(args, "genetic_map_hg19_sources", None))
    if not sources:
        return "bim_cm"
    return {
        "kind": "explicit_hg19",
        "sources": [
            {"name": Path(path).name, "sha256": _file_sha256(path)}
            for path in sources
        ],
    }


def _builder_suite_identity(config, args, chromosomes, regression_path: Path) -> dict:
    from .path_resolution import resolve_file_group

    baseline_files = resolve_file_group(
        config.baseline_annot_sources,
        label="baseline annotation",
        allow_chromosome_suite=True,
    )
    plink_files = []
    removal_counts = {}
    first_prefix = None
    for chrom, record in chromosomes.items():
        kernel_args = argparse.Namespace(bfile=config.plink_prefix)
        prefix = kernel_ldscore.resolve_bfile_prefix(kernel_args, chrom=chrom)
        if first_prefix is None:
            first_prefix = prefix
        for suffix in (".bed", ".bim", ".fam"):
            path = str(prefix) + suffix
            plink_files.append({"chromosome": chrom, "kind": suffix[1:], "sha256": _file_sha256(path)})
        pre_qc = len(_read_bim_identity(str(prefix)))
        removal_counts[chrom] = {
            "pre_qc": pre_qc,
            "retained": int(record.total_reference_snps_all),
            "genotype_or_maf_removed": pre_qc - int(record.total_reference_snps_all),
            "regression": int(len(record.baseline_rows)),
        }
    assert first_prefix is not None
    fam = kernel_ldscore.legacy_parse.PlinkFAMFile(str(first_prefix) + ".fam")
    keep_indices = kernel_ldscore.resolve_keep_individuals(config.keep_indivs_file, fam)
    selected_indices = list(range(len(fam.IDList))) if keep_indices is None else list(keep_indices)
    selected_ids = _selected_individual_ids(fam, selected_indices)
    selected_digest = hashlib.sha256(
        json.dumps(selected_ids, separators=(",", ":"), ensure_ascii=True).encode("utf-8")
    ).hexdigest()
    return {
        "baseline_sources": [{"name": Path(path).name, "sha256": _file_sha256(path)} for path in baseline_files],
        "plink_sources": plink_files,
        "chromosomes": list(config.chromosomes),
        "genome_build": config.genome_build,
        "snp_identifier": config.snp_identifier,
        "plink_suite": "1000G_EUR_Phase3",
        "selected_individuals": {
            "source": "all" if config.keep_indivs_file is None else "keep_file",
            "selected_content_sha256": selected_digest,
            "selected_count": len(selected_ids),
            "keep_file_sha256": None if config.keep_indivs_file is None else _file_sha256(config.keep_indivs_file),
        },
        "maf_min": config.maf_min,
        "common_maf_min": config.common_maf_min,
        "ld_window": {"unit": "cm", "value": config.ld_wind_cm},
        "genetic_map": _genetic_map_identity(args),
        "regression_snps": {"kind": "bundled_hm3", "sha256": _file_sha256(regression_path)},
        "snp_exclude_regions": "mhc-and-centromeres",
        "removal_counts": removal_counts,
    }


def _selected_individual_ids(fam, selected_indices: Sequence[int]) -> list[str]:
    """Return selected PLINK individual IDs in their effective genotype order."""
    identifiers = fam.IDList
    if not isinstance(identifiers, pd.DataFrame) or identifiers.shape[1] != 1:
        raise LDSCInternalError("PLINK FAM individual identity table must contain exactly one IID column.")
    return identifiers.iloc[list(selected_indices), 0].astype(str).tolist()


@dataclass(frozen=True)
class IndexChromosomeData:
    """Validated in-memory common/profile payload for one chromosome."""

    baseline_rows: pd.DataFrame
    baseline_count_all: np.ndarray
    baseline_count_common: np.ndarray
    baseline_overlap_all: np.ndarray
    baseline_overlap_common: np.ndarray
    total_reference_snps_all: int
    total_reference_snps_common: int
    atom_model: ChromosomeAtomModel
    operator: sparse.csr_matrix
    atom_statistics: AtomStatistics
    reference_metadata: pd.DataFrame
    reference_rows_before_genotype_qc: int = 0
    genotype_qc_removed: int = 0
    maf_removed: int = 0
    selected_individual_count: int = 0
    cm_source: str = "bim_cm"


@dataclass(frozen=True)
class LoadedGeneLDScoreIndex:
    """A validated profile and its immutable common LD-score layer.

    Attributes
    ----------
    suite_id, profile_id : str
        Content-bound scientific identities for the common suite and selected
        gene profile.
    chromosomes : tuple of str
        Ordered immutable chromosome coverage.
    suite_identity, profile_identity : dict
        Scientific settings used to calculate the semantic identities.
    gene_catalog : pandas.DataFrame
        Embedded catalog with unpadded 0-based half-open gene coordinates,
        inclusion flags, and chromosome-local row ordering.
    profile_chromosomes : dict of str to IndexChromosomeData
        Validated chromosome-local baseline rows, sparse operators, atom
        geometry, and sufficient statistics.
    """

    suite_id: str
    profile_id: str
    chromosomes: tuple[str, ...]
    suite_identity: dict
    profile_identity: dict
    gene_catalog: pd.DataFrame
    profile_chromosomes: dict[str, IndexChromosomeData]


def publish_gene_ldscore_index(
    suite_dir: str | Path,
    *,
    profile_name: str,
    suite_identity: dict,
    profile_identity: dict,
    gene_catalog: pd.DataFrame,
    chromosomes: dict[str, IndexChromosomeData],
    overwrite: bool,
) -> Path:
    """Stage, validate, and atomically publish one profile without losing siblings."""
    destination = Path(suite_dir)
    destination.parent.mkdir(parents=True, exist_ok=True)
    suite_id = calculate_suite_id(suite_identity)
    profile_id = calculate_profile_id(suite_id, profile_identity)
    stage_parent = Path(tempfile.mkdtemp(prefix=f".{destination.name}.stage-", dir=destination.parent))
    staged_suite = stage_parent / destination.name
    try:
        if destination.exists() and not _is_abandoned_gene_index_suite(destination):
            root = _read_json(destination / "metadata.json", "suite root")
            _validate_metadata_identity(root, suite_id=suite_id)
            shutil.copytree(destination, staged_suite)
        else:
            if destination.exists():
                shutil.rmtree(destination)
            staged_suite.mkdir()
            _write_json(
                staged_suite / "metadata.json",
                {
                    "artifact_type": "gene_ldscore_index",
                    "suite_id": suite_id,
                    "suite_identity": suite_identity,
                    "chromosomes": list(chromosomes),
                },
            )
            _write_common_layer(staged_suite, suite_id, suite_identity, chromosomes)
        profile_dir = staged_suite / "profiles" / profile_name
        if profile_dir.exists() and not overwrite:
            if _is_abandoned_gene_index_profile(profile_dir):
                shutil.rmtree(profile_dir)
            else:
                raise FileExistsError(
                    f"Gene LD-score index profile already exists: {destination / 'profiles' / profile_name}. "
                    "Pass --overwrite to replace only this profile."
                )
        if profile_dir.exists():
            shutil.rmtree(profile_dir)
        _write_profile_layer(
            profile_dir,
            suite_id=suite_id,
            profile_id=profile_id,
            profile_identity=profile_identity,
            gene_catalog=gene_catalog,
            chromosomes=chromosomes,
        )
        load_gene_ldscore_index(profile_dir)
        backup = stage_parent / f"{destination.name}.backup"
        if destination.exists():
            os.replace(destination, backup)
        try:
            os.replace(staged_suite, destination)
        except Exception:
            if backup.exists():
                os.replace(backup, destination)
            raise
        if backup.exists():
            shutil.rmtree(backup)
        return destination / "profiles" / profile_name
    finally:
        if stage_parent.exists():
            shutil.rmtree(stage_parent)


def load_gene_ldscore_index(profile_dir: str | Path) -> LoadedGeneLDScoreIndex:
    """Load and fully validate one explicit gene LD-score index profile.

    Parameters
    ----------
    profile_dir : path-like
        Exact ``<suite>/profiles/<profile>`` directory. A suite root is not
        accepted and no profile discovery is performed.

    Returns
    -------
    LoadedGeneLDScoreIndex
        Common and profile components after validating semantic identities,
        chromosome coverage, catalog and row ordering, shapes, dtypes, finite
        values, and CSR invariants.

    Raises
    ------
    LDSCInputError
        If any required component is missing, corrupt, structurally invalid,
        or semantically incompatible.
    """
    profile_path = Path(profile_dir)
    if profile_path.parent.name != "profiles":
        raise LDSCInputError(
            f"Gene LD-score index path must name one profile directory beneath profiles/: {profile_path}"
        )
    suite_path = profile_path.parent.parent
    root = _read_json(suite_path / "metadata.json", "suite root")
    if root.get("artifact_type") != "gene_ldscore_index":
        raise LDSCInputError("Gene LD-score index root has the wrong artifact_type.")
    suite_identity = root.get("suite_identity")
    suite_id = root.get("suite_id")
    if not isinstance(suite_identity, dict) or suite_id != calculate_suite_id(suite_identity):
        raise LDSCInputError("Gene LD-score index suite_id does not match its scientific identity.")
    common_meta = _read_json(suite_path / "common" / "metadata.json", "common layer")
    _validate_metadata_identity(common_meta, suite_id=suite_id)
    root_chromosomes = tuple(str(chrom) for chrom in root.get("chromosomes", ()))
    if (
        common_meta.get("suite_identity") != suite_identity
        or tuple(str(chrom) for chrom in common_meta.get("chromosomes", ())) != root_chromosomes
    ):
        raise LDSCInputError("Gene LD-score index common chromosome coverage or identity is invalid.")
    profile_meta = _read_json(profile_path / "metadata.json", "profile")
    profile_identity = profile_meta.get("profile_identity")
    profile_id = profile_meta.get("profile_id")
    if (
        profile_meta.get("artifact_type") != "gene_ldscore_index"
        or profile_meta.get("suite_id") != suite_id
        or not isinstance(profile_identity, dict)
        or profile_id != calculate_profile_id(suite_id, profile_identity)
    ):
        raise LDSCInputError("Gene LD-score index profile semantic identity is invalid.")
    chromosomes = tuple(str(chrom) for chrom in profile_meta.get("chromosomes", ()))
    if chromosomes != root_chromosomes:
        raise LDSCInputError("Gene LD-score index profile chromosome coverage does not match the suite.")
    catalog_path = profile_path / "gene_catalog.parquet"
    if not catalog_path.exists():
        raise LDSCInputError("Gene LD-score index profile is missing gene_catalog.parquet.")
    catalog = pd.read_parquet(catalog_path)
    required_catalog = [
        "gene_index", "canonical_ensembl_id", "gene_name", "CHR", "start0", "end",
        "included", "exclusion_reason", "chromosome_gene_row",
    ]
    if list(catalog.columns) != required_catalog:
        raise LDSCInputError("Gene LD-score index catalog schema is invalid.")
    expected_gene_indices = np.arange(len(catalog), dtype=np.int64)
    if not np.array_equal(catalog["gene_index"].to_numpy(), expected_gene_indices):
        raise LDSCInputError("Gene LD-score index catalog gene_index ordering is invalid.")
    if catalog["canonical_ensembl_id"].astype(str).duplicated().any():
        raise LDSCInputError("Gene LD-score index catalog contains duplicate canonical gene identifiers.")
    if catalog["included"].dtype != np.bool_:
        raise LDSCInputError("Gene LD-score index catalog included flags must use Boolean dtype.")
    loaded: dict[str, IndexChromosomeData] = {}
    for chrom in chromosomes:
        chromosome_catalog = catalog.loc[catalog["CHR"].astype(str) == str(chrom)]
        expected_chromosome_rows = np.arange(len(chromosome_catalog), dtype=np.int64)
        if not np.array_equal(
            chromosome_catalog["chromosome_gene_row"].to_numpy(), expected_chromosome_rows
        ):
            raise LDSCInputError(
                f"Gene LD-score index catalog chromosome_gene_row ordering is invalid for chromosome {chrom}."
            )
        loaded[chrom] = _load_index_chromosome(
            suite_path,
            profile_path,
            chrom,
            suite_id=suite_id,
            profile_id=profile_id,
            expected_gene_rows=len(chromosome_catalog),
        )
    return LoadedGeneLDScoreIndex(
        suite_id=suite_id,
        profile_id=profile_id,
        chromosomes=chromosomes,
        suite_identity=suite_identity,
        profile_identity=profile_identity,
        gene_catalog=catalog,
        profile_chromosomes=loaded,
    )


def run_indexed_ldscore(
    profile_dir: str | Path,
    *,
    query_gene_list_sources: Sequence[str | Path],
    control_gene_list_source: str | Path = "all-protein-coding",
    output_dir: str | Path,
    overwrite: bool = False,
):
    """Assemble exact gene-list LD scores from an explicit index profile.

    Parameters
    ----------
    profile_dir : path-like
        Exact validated profile directory beneath ``profiles/``.
    query_gene_list_sources : sequence of path-like
        One-column gene-list files. Exact Ensembl IDs and case-sensitive gene
        names are resolved against the embedded catalog in source order.
    control_gene_list_source : path-like, {"all-protein-coding", "none"}, optional
        Fixed control source. The default adds the union of all eligible
        protein-coding genes as ``gene_control`` to the baseline block.
    output_dir : path-like
        Destination for the canonical self-contained LD-score directory.
    overwrite : bool, optional
        Replace workflow-owned output artifacts. Default is ``False``.

    Returns
    -------
    LDScoreResult
        Assembled baseline/control and focal-query tables, counts, overlaps,
        statuses, provenance, and canonical output paths.

    Raises
    ------
    LDSCInputError
        If the profile is invalid, a control is unusable, or every requested
        focal query is skipped. The all-skipped case writes diagnostics but no
        scientific LD-score tables.

    Notes
    -----
    Overlapping, nested, duplicated, and alias-selected genes are combined by
    Boolean union. All requested focal columns are assembled together in
    float64 and are narrowed only by the canonical Parquet writer.
    """
    from .config import GlobalConfig
    from .ldscore_calculator import LDScoreCalculator, LDScoreResult
    from .outputs import LDScoreDirectoryWriter, LDScoreOutputConfig
    from .overlap_matrix import LDScoreOverlap

    index = load_gene_ldscore_index(profile_dir)
    profile_catalog = dict(index.profile_identity.get("catalog", {}))
    projection_build = str(index.profile_identity.get("projection_build"))
    gene_policy = str(index.profile_identity.get("gene_exclude_regions", "none"))
    catalog = GeneCatalog.from_index_frame(
        index.gene_catalog,
        genome_build=projection_build,
        release=str(profile_catalog.get("release", "GENCODE v49")),
        content_sha256=str(profile_catalog.get("content_sha256", "")),
    )
    paths = resolve_file_group(
        query_gene_list_sources,
        label="gene-list file",
        allow_chromosome_suite=False,
    )
    resolutions = tuple(
        resolve_gene_list(
            path,
            catalog,
            genome_build=projection_build,
            source_ordinal=ordinal,
            gene_exclude_regions=gene_policy,
        )
        for ordinal, path in enumerate(paths, start=1)
    )
    statuses = tuple(
        QueryAnnotationStatus(
            resolution.query,
            resolution.source,
            "gene_list",
            resolution.status,
            resolution.reason,
            details=resolution.details,
        )
        for resolution in resolutions
    )
    usable = [resolution for resolution in resolutions if resolution.status in {"ok", "warning"}]
    if not usable:
        _write_indexed_all_skipped_diagnostics(
            statuses,
            resolutions,
            output_dir=output_dir,
            overwrite=overwrite,
        )
    control_source = str(control_gene_list_source)
    control_resolution: GeneListResolution | None
    if control_source == "none":
        control_resolution = None
    elif control_source == "all-protein-coding":
        control_resolution = resolve_all_protein_coding(
            catalog,
            genome_build=projection_build,
            gene_exclude_regions=gene_policy,
        )
    else:
        control_resolution = resolve_gene_list(
            control_source,
            catalog,
            genome_build=projection_build,
            source_ordinal=0,
            gene_exclude_regions=gene_policy,
        )
        if control_resolution.status not in {"ok", "warning"}:
            raise LDSCInputError(
                f"The requested indexed control gene list is unusable (reason={control_resolution.reason})."
            )
    if any(resolution.query == "gene_control" for resolution in usable):
        raise LDSCInputError("Focal query name 'gene_control' collides with the reserved fixed control column.")

    supplied_baseline = _baseline_columns_from_rows(next(iter(index.profile_chromosomes.values())).baseline_rows)
    if "gene_control" in supplied_baseline:
        raise LDSCInputError("Index supplied baseline columns contain reserved name 'gene_control'.")
    baseline_columns = [*supplied_baseline, *(("gene_control",) if control_resolution is not None else ())]
    query_columns = [resolution.query for resolution in usable]
    baseline_tables: list[pd.DataFrame] = []
    query_tables: list[pd.DataFrame] = []
    baseline_count_all = np.zeros(len(baseline_columns), dtype=np.float64)
    baseline_count_common = np.zeros(len(baseline_columns), dtype=np.float64)
    query_count_all = np.zeros(len(query_columns), dtype=np.float64)
    query_count_common = np.zeros(len(query_columns), dtype=np.float64)
    block_all = np.zeros((len(baseline_columns), len(baseline_columns) + len(query_columns)), dtype=np.float64)
    block_common = np.zeros_like(block_all)
    query_diag_all = np.zeros(len(query_columns), dtype=np.float64)
    query_diag_common = np.zeros(len(query_columns), dtype=np.float64)
    total_all = total_common = 0

    for chrom in index.chromosomes:
        record = index.profile_chromosomes[chrom]
        selectors = [_selector_for_resolution(index.gene_catalog, record.atom_model, chrom, resolution) for resolution in usable]
        control_selector = (
            None
            if control_resolution is None
            else _selector_for_resolution(index.gene_catalog, record.atom_model, chrom, control_resolution)
        )
        baseline_table = record.baseline_rows.copy()
        if control_selector is not None:
            baseline_table["gene_control"] = assemble_indexed_ld_scores(record.operator, control_selector)
        identity_columns = [column for column in ("CHR", "SNP", "POS", "A1", "A2") if column in baseline_table.columns]
        query_table = baseline_table.loc[:, identity_columns].copy()
        if selectors:
            selector_matrix = np.column_stack(selectors)
            scores = assemble_indexed_ld_scores(record.operator, selector_matrix)
            for idx, column in enumerate(query_columns):
                query_table[column] = scores[:, idx]
        baseline_tables.append(baseline_table)
        query_tables.append(query_table)

        b = len(supplied_baseline)
        baseline_count_all[:b] += record.baseline_count_all
        baseline_count_common[:b] += record.baseline_count_common
        block_all[:b, :b] += record.baseline_overlap_all
        block_common[:b, :b] += record.baseline_overlap_common
        atom_all = record.atom_statistics.atom_count_all.astype(np.float64)
        atom_common = record.atom_statistics.atom_count_common.astype(np.float64)
        baseline_atom_all = record.atom_statistics.baseline_atom_overlap_all
        baseline_atom_common = record.atom_statistics.baseline_atom_overlap_common
        if control_selector is not None:
            control_pos = b
            baseline_count_all[control_pos] += atom_all @ control_selector
            baseline_count_common[control_pos] += atom_common @ control_selector
            base_control_all = baseline_atom_all @ control_selector.astype(np.float64)
            base_control_common = baseline_atom_common @ control_selector.astype(np.float64)
            block_all[:b, control_pos] += base_control_all
            block_all[control_pos, :b] += base_control_all
            block_common[:b, control_pos] += base_control_common
            block_common[control_pos, :b] += base_control_common
            block_all[control_pos, control_pos] += atom_all @ control_selector
            block_common[control_pos, control_pos] += atom_common @ control_selector
        for query_pos, selector in enumerate(selectors):
            column_pos = len(baseline_columns) + query_pos
            query_count_all[query_pos] += atom_all @ selector
            query_count_common[query_pos] += atom_common @ selector
            query_diag_all[query_pos] += atom_all @ selector
            query_diag_common[query_pos] += atom_common @ selector
            block_all[:b, column_pos] += baseline_atom_all @ selector.astype(np.float64)
            block_common[:b, column_pos] += baseline_atom_common @ selector.astype(np.float64)
            if control_selector is not None:
                intersection = selector & control_selector
                block_all[b, column_pos] += atom_all @ intersection
                block_common[b, column_pos] += atom_common @ intersection
        total_all += record.total_reference_snps_all
        total_common += record.total_reference_snps_common

    baseline_table = pd.concat(baseline_tables, ignore_index=True)
    query_table = pd.concat(query_tables, ignore_index=True)
    count_records = [
        {
            "group": "baseline",
            "column": column,
            "all_reference_snp_count": float(baseline_count_all[idx]),
            "common_reference_snp_count": float(baseline_count_common[idx]),
        }
        for idx, column in enumerate(baseline_columns)
    ] + [
        {
            "group": "query",
            "column": column,
            "all_reference_snp_count": float(query_count_all[idx]),
            "common_reference_snp_count": float(query_count_common[idx]),
        }
        for idx, column in enumerate(query_columns)
    ]
    all_columns = [*baseline_columns, *query_columns]
    overlap = LDScoreOverlap(
        baseline_block_all=pd.DataFrame(block_all, index=baseline_columns, columns=all_columns),
        baseline_block_common=pd.DataFrame(block_common, index=baseline_columns, columns=all_columns),
        query_diagonal_all=pd.Series(query_diag_all, index=query_columns),
        query_diagonal_common=pd.Series(query_diag_common, index=query_columns),
        total_all_reference_snps=float(total_all),
        total_common_reference_snps=float(total_common),
    )
    result = LDScoreResult(
        baseline_table=baseline_table,
        query_table=query_table,
        count_records=count_records,
        baseline_columns=baseline_columns,
        query_columns=query_columns,
        ld_reference_snps=frozenset(),
        ld_regression_snps=frozenset(baseline_table["SNP"].astype(str)),
        chromosome_results=[],
        count_config={"common_reference_snp_maf_min": 0.05, "common_reference_snp_maf_operator": ">="},
        config_snapshot=GlobalConfig(snp_identifier="rsid"),
        overlap=overlap,
        query_statuses=statuses,
        gene_list_resolutions=resolutions,
        control_gene_list_resolution=control_resolution,
        gene_catalog_provenance={
            "resource": "gene_catalog.parquet",
            "release": catalog.release,
            "genome_build": projection_build,
            "content_sha256": catalog.content_sha256,
        },
        index_provenance={"suite_id": index.suite_id, "profile_id": index.profile_id},
    )
    result = LDScoreCalculator()._finalize_query_statuses(result, statuses)
    if not result.query_columns:
        LDScoreDirectoryWriter().write_query_diagnostics(
            result,
            LDScoreOutputConfig(output_dir=output_dir, overwrite=overwrite),
        )
        from .ldscore_calculator import _all_query_annotations_skipped_message

        raise LDSCInputError(_all_query_annotations_skipped_message(result.query_statuses))
    result.validate()
    output_paths = LDScoreDirectoryWriter().write(
        result,
        LDScoreOutputConfig(output_dir=output_dir, overwrite=overwrite),
    )
    from dataclasses import replace

    return replace(result, output_paths=output_paths)


def _write_indexed_all_skipped_diagnostics(
    statuses: tuple[QueryAnnotationStatus, ...],
    resolutions: tuple[GeneListResolution, ...],
    *,
    output_dir: str | Path,
    overwrite: bool,
) -> None:
    """Write the ordinary diagnostic-only family and raise the consolidated error."""
    from types import SimpleNamespace
    from .ldscore_calculator import _all_query_annotations_skipped_message
    from .outputs import LDScoreDirectoryWriter, LDScoreOutputConfig

    diagnostic_result = SimpleNamespace(
        query_statuses=statuses,
        gene_list_resolutions=resolutions,
        control_gene_list_resolution=None,
    )
    LDScoreDirectoryWriter().write_query_diagnostics(
        diagnostic_result,
        LDScoreOutputConfig(output_dir=output_dir, overwrite=overwrite),
    )
    raise LDSCInputError(_all_query_annotations_skipped_message(statuses))


def _baseline_columns_from_rows(frame: pd.DataFrame) -> list[str]:
    """Return ordered supplied baseline score columns from one common shard."""
    reserved = {"CHR", "SNP", "POS", "A1", "A2", "regression_ld_scores"}
    return [column for column in frame.columns if column not in reserved]


def _selector_for_resolution(
    embedded_catalog: pd.DataFrame,
    atom_model: ChromosomeAtomModel,
    chrom: str,
    resolution: GeneListResolution,
) -> np.ndarray:
    """Translate global embedded-catalog indices into one chromosome atom selector."""
    selected = set(int(index) for index in resolution.catalog_indices)
    rows = embedded_catalog.loc[
        embedded_catalog.index.isin(selected)
        & (embedded_catalog["CHR"].astype(str) == str(chrom))
        & embedded_catalog["included"].astype(bool),
        "chromosome_gene_row",
    ].astype(int)
    return assemble_atom_selector(atom_model.gene_to_atom, rows.to_numpy(dtype=np.int64))


def build_plink_index_chromosome(
    chrom: str,
    baseline_bundle: kernel_ldscore.AnnotationBundle,
    args,
    *,
    regression_keys: set[str] | RestrictionIdentityKeys,
    regression_regions: RegionIntervals | None,
    gene_intervals: np.ndarray,
    included: np.ndarray,
    padding_bp: int,
    atom_batch_size: int,
) -> IndexChromosomeData:
    """Construct exact common and atom payloads from one prepared PLINK chromosome."""
    prepared = kernel_ldscore.prepare_plink_chromosome(chrom, baseline_bundle, args)
    metadata = prepared.metadata
    baseline = np.asarray(prepared.annotation_matrix, dtype=np.float64)
    geno = prepared.geno

    geno._currentSNP = 0
    baseline_scores = np.asarray(
        geno.ldScoreVarBlocks(prepared.block_left, args.snp_batch_size, annot=baseline),
        dtype=np.float64,
    )
    persisted = kernel_ldscore.regression_mask_from_keys(
        metadata,
        regression_keys,
        args.snp_identifier,
        region_intervals=regression_regions,
    ).astype(bool)
    geno._currentSNP = 0
    regression_scores = np.asarray(
        geno.ldScoreVarBlocks(
            prepared.block_left,
            args.snp_batch_size,
            annot=persisted.reshape(-1, 1),
        ),
        dtype=np.float64,
    ).reshape(-1)

    baseline_frame = pd.DataFrame(baseline, columns=baseline_bundle.baseline_columns)
    baseline_count_all, baseline_count_common = kernel_ldscore.compute_counts(
        metadata,
        baseline_frame,
        common_maf_min=getattr(args, "common_maf_min", 0.05),
    )
    overlap = kernel_ldscore.compute_overlap(
        metadata,
        baseline_frame,
        n_baseline=len(baseline_bundle.baseline_columns),
        common_maf_min=getattr(args, "common_maf_min", 0.05),
    )
    common_mask = metadata["MAF"].to_numpy(dtype=np.float64) >= float(
        getattr(args, "common_maf_min", 0.05)
    )

    atom_model = build_disjoint_atoms(
        chrom,
        gene_intervals,
        included=included,
        padding_bp=padding_bp,
    )
    snp_atoms = map_snps_to_atoms(metadata["POS"].to_numpy(dtype=np.int64) - 1, atom_model)
    atom_statistics = compute_atom_statistics(snp_atoms, atom_model.n_atoms, baseline, common_mask)
    operator_blocks: list[sparse.csr_matrix] = []
    for _start, _end, atom_block in iter_snp_atom_blocks(
        snp_atoms,
        atom_model.n_atoms,
        atom_batch_size,
    ):
        geno._currentSNP = 0
        block_scores = np.asarray(
            geno.ldScoreVarBlocks(
                prepared.block_left,
                args.snp_batch_size,
                annot=atom_block,
            ),
            dtype=np.float64,
        )
        operator_blocks.append(sparse.csr_matrix(block_scores[persisted], dtype=np.float64))
    operator = (
        sparse.hstack(operator_blocks, format="csr", dtype=np.float64)
        if operator_blocks
        else sparse.csr_matrix((int(persisted.sum()), 0), dtype=np.float64)
    )
    operator.sort_indices()
    validate_ldscore_operator(operator, n_rows=int(persisted.sum()), n_atoms=atom_model.n_atoms)

    identity_columns = ["CHR", "SNP", "POS", *[c for c in ("A1", "A2") if c in metadata.columns]]
    baseline_rows = metadata.loc[persisted, identity_columns].reset_index(drop=True)
    baseline_rows["regression_ld_scores"] = regression_scores[persisted]
    for column_index, column in enumerate(baseline_bundle.baseline_columns):
        baseline_rows[column] = baseline_scores[persisted, column_index]
    return IndexChromosomeData(
        baseline_rows=baseline_rows,
        baseline_count_all=np.asarray(baseline_count_all, dtype=np.float64),
        baseline_count_common=np.asarray(baseline_count_common, dtype=np.float64),
        baseline_overlap_all=np.asarray(overlap.baseline_block_all, dtype=np.float64),
        baseline_overlap_common=np.asarray(overlap.baseline_block_common, dtype=np.float64),
        total_reference_snps_all=len(metadata),
        total_reference_snps_common=int(common_mask.sum()),
        atom_model=atom_model,
        operator=operator,
        atom_statistics=atom_statistics,
        reference_metadata=metadata.copy(),
        reference_rows_before_genotype_qc=prepared.reference_rows_before_genotype_qc,
        genotype_qc_removed=prepared.genotype_qc_removed,
        maf_removed=prepared.maf_removed,
        selected_individual_count=prepared.selected_individual_count,
        cm_source=prepared.cm_source,
    )


def validate_strict_baseline_plink_identity(
    baseline_rows: pd.DataFrame,
    plink_rows: pd.DataFrame,
    *,
    chrom: str,
) -> pd.DataFrame:
    """Require exact pre-genotype-QC ``CHR/POS/SNP`` identity equality.

    Input order may differ. The returned baseline frame is canonically sorted
    by position and SNP while preserving all baseline columns.
    """
    baseline = _normalize_identity_rows(baseline_rows, label="baseline", chrom=chrom)
    plink = _normalize_identity_rows(plink_rows, label="PLINK BIM", chrom=chrom)
    for label, frame in (("baseline", baseline), ("PLINK BIM", plink)):
        duplicate = frame.duplicated(_STRICT_IDENTITY_COLUMNS, keep=False)
        if duplicate.any():
            raise LDSCInputError(
                f"Gene LD-score index build found duplicate {label} CHR/POS/SNP identities "
                f"on chromosome {chrom}. Baseline identities must equal PLINK identities before genotype QC."
            )
    baseline_sorted = baseline.sort_values(["POS", "SNP"], kind="mergesort").reset_index(drop=True)
    plink_sorted = plink.sort_values(["POS", "SNP"], kind="mergesort").reset_index(drop=True)
    baseline_position_keys = set(zip(baseline_sorted["CHR"], baseline_sorted["POS"]))
    plink_position_keys = set(zip(plink_sorted["CHR"], plink_sorted["POS"]))
    if baseline_position_keys == plink_position_keys and not baseline_sorted[_STRICT_IDENTITY_COLUMNS].equals(
        plink_sorted[_STRICT_IDENTITY_COLUMNS]
    ):
        raise LDSCInputError(
            f"Gene LD-score index build found conflicting SNP identities at matching CHR/POS rows "
            f"on chromosome {chrom}. Baseline identities must equal PLINK identities before genotype QC."
        )
    baseline_keys = set(map(tuple, baseline_sorted[_STRICT_IDENTITY_COLUMNS].to_numpy()))
    plink_keys = set(map(tuple, plink_sorted[_STRICT_IDENTITY_COLUMNS].to_numpy()))
    if baseline_keys != plink_keys:
        raise LDSCInputError(
            f"Gene LD-score index build found missing or extra baseline/PLINK identities on chromosome {chrom}: "
            f"baseline_only={len(baseline_keys - plink_keys)}, plink_only={len(plink_keys - baseline_keys)}. "
            "Baseline identities must equal PLINK identities before genotype QC."
        )
    return baseline_sorted


def _normalize_identity_rows(frame: pd.DataFrame, *, label: str, chrom: str) -> pd.DataFrame:
    """Normalize one strict-identity frame without dropping caller columns."""
    missing = set(_STRICT_IDENTITY_COLUMNS) - set(frame.columns)
    if missing:
        raise LDSCInputError(
            f"Gene LD-score index build cannot compare {label} identities because columns {sorted(missing)} are missing."
        )
    normalized = frame.copy()
    normalized["CHR"] = normalized["CHR"].map(normalize_chromosome).astype(str)
    normalized["POS"] = pd.to_numeric(normalized["POS"], errors="raise").astype("int64")
    normalized["SNP"] = normalized["SNP"].astype(str)
    expected_chrom = normalize_chromosome(chrom)
    normalized = normalized.loc[normalized["CHR"] == expected_chrom].reset_index(drop=True)
    return normalized


def _write_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")


def _read_json(path: Path, label: str) -> dict:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise LDSCInputError(f"Gene LD-score index {label} metadata is missing or unreadable: {path}") from exc
    if not isinstance(payload, dict):
        raise LDSCInputError(f"Gene LD-score index {label} metadata must be a JSON object.")
    return payload


def _validate_metadata_identity(metadata: dict, *, suite_id: str, profile_id: str | None = None) -> None:
    if metadata.get("artifact_type") != "gene_ldscore_index" or metadata.get("suite_id") != suite_id:
        raise LDSCInputError("Gene LD-score index component suite identity is invalid.")
    if profile_id is not None and metadata.get("profile_id") != profile_id:
        raise LDSCInputError("Gene LD-score index component profile identity is invalid.")


def _cast_index_baseline_rows(frame: pd.DataFrame) -> pd.DataFrame:
    """Match the canonical LD-score writer's persisted floating dtypes."""
    float64_columns = [column for column in frame.columns if frame[column].dtype == np.float64]
    return frame.astype({column: np.float32 for column in float64_columns}) if float64_columns else frame


def _write_common_layer(
    suite_path: Path,
    suite_id: str,
    suite_identity: dict,
    chromosomes: dict[str, IndexChromosomeData],
) -> None:
    common_path = suite_path / "common"
    _write_json(
        common_path / "metadata.json",
        {
            "artifact_type": "gene_ldscore_index",
            "suite_id": suite_id,
            "suite_identity": suite_identity,
            "chromosomes": list(chromosomes),
        },
    )
    for chrom, record in chromosomes.items():
        chrom_path = common_path / f"chr{chrom}"
        chrom_path.mkdir(parents=True)
        baseline_columns = [
            column
            for column in record.baseline_rows.columns
            if column not in {"CHR", "SNP", "POS", "A1", "A2", "regression_ld_scores"}
        ]
        _write_json(
            chrom_path / "metadata.json",
            {
                "artifact_type": "gene_ldscore_index",
                "suite_id": suite_id,
                "chromosome": str(chrom),
                "n_rows": len(record.baseline_rows),
                "n_baseline": len(baseline_columns),
                "baseline_columns": baseline_columns,
            },
        )
        _cast_index_baseline_rows(record.baseline_rows).to_parquet(
            chrom_path / "baseline_rows.parquet", index=False
        )
        np.savez(
            chrom_path / "baseline_statistics.npz",
            baseline_count_all=np.asarray(record.baseline_count_all, dtype=np.float64),
            baseline_count_common=np.asarray(record.baseline_count_common, dtype=np.float64),
            baseline_overlap_all=np.asarray(record.baseline_overlap_all, dtype=np.float64),
            baseline_overlap_common=np.asarray(record.baseline_overlap_common, dtype=np.float64),
            total_reference_snps_all=np.asarray(record.total_reference_snps_all, dtype=np.int64),
            total_reference_snps_common=np.asarray(record.total_reference_snps_common, dtype=np.int64),
        )


def _write_profile_layer(
    profile_path: Path,
    *,
    suite_id: str,
    profile_id: str,
    profile_identity: dict,
    gene_catalog: pd.DataFrame,
    chromosomes: dict[str, IndexChromosomeData],
) -> None:
    profile_path.mkdir(parents=True)
    _write_json(
        profile_path / "metadata.json",
        {
            "artifact_type": "gene_ldscore_index",
            "suite_id": suite_id,
            "profile_id": profile_id,
            "profile_identity": profile_identity,
            "chromosomes": list(chromosomes),
        },
    )
    gene_catalog.to_parquet(profile_path / "gene_catalog.parquet", index=False)
    for chrom, record in chromosomes.items():
        chrom_path = profile_path / f"chr{chrom}"
        chrom_path.mkdir()
        atoms = pd.DataFrame(
            {
                "atom_id": np.arange(record.atom_model.n_atoms, dtype=np.int64),
                "CHR": str(chrom),
                "start0": record.atom_model.starts,
                "end": record.atom_model.ends,
            }
        )
        _write_json(
            chrom_path / "metadata.json",
            {
                "artifact_type": "gene_ldscore_index",
                "suite_id": suite_id,
                "profile_id": profile_id,
                "chromosome": str(chrom),
                "n_genes": record.atom_model.gene_to_atom.shape[0],
                "n_atoms": record.atom_model.n_atoms,
                "n_rows": record.operator.shape[0],
                "gene_to_atom_format": "csr_bool_int32",
                "ldscore_operator_format": "csr_float64_int32",
            },
        )
        atoms.to_parquet(chrom_path / "atoms.parquet", index=False)
        gene_to_atom = record.atom_model.gene_to_atom.astype(bool).tocsr()
        gene_to_atom.sort_indices()
        operator = record.operator.astype(np.float64).tocsr()
        operator.sort_indices()
        sparse.save_npz(chrom_path / "gene_to_atom.npz", gene_to_atom, compressed=True)
        sparse.save_npz(chrom_path / "ldscore_operator.npz", operator, compressed=True)
        np.savez(
            chrom_path / "atom_statistics.npz",
            atom_count_all=np.asarray(record.atom_statistics.atom_count_all, dtype=np.int64),
            atom_count_common=np.asarray(record.atom_statistics.atom_count_common, dtype=np.int64),
            baseline_atom_overlap_all=np.asarray(record.atom_statistics.baseline_atom_overlap_all, dtype=np.float64),
            baseline_atom_overlap_common=np.asarray(record.atom_statistics.baseline_atom_overlap_common, dtype=np.float64),
        )


def _load_npz_members(path: Path, required: set[str]) -> dict[str, np.ndarray]:
    try:
        with np.load(path, allow_pickle=False) as archive:
            if set(archive.files) != required:
                raise LDSCInputError(
                    f"Gene LD-score index payload {path.name} has required members {sorted(required)}, "
                    f"found {sorted(archive.files)}."
                )
            return {name: archive[name] for name in archive.files}
    except (OSError, ValueError) as exc:
        raise LDSCInputError(f"Gene LD-score index payload is unreadable: {path}") from exc


def _load_index_chromosome(
    suite_path: Path,
    profile_path: Path,
    chrom: str,
    *,
    suite_id: str,
    profile_id: str,
    expected_gene_rows: int,
) -> IndexChromosomeData:
    common_path = suite_path / "common" / f"chr{chrom}"
    component_path = profile_path / f"chr{chrom}"
    common_meta = _read_json(common_path / "metadata.json", f"common chromosome {chrom}")
    profile_meta = _read_json(component_path / "metadata.json", f"profile chromosome {chrom}")
    _validate_metadata_identity(common_meta, suite_id=suite_id)
    _validate_metadata_identity(profile_meta, suite_id=suite_id, profile_id=profile_id)
    if common_meta.get("chromosome") != str(chrom) or profile_meta.get("chromosome") != str(chrom):
        raise LDSCInputError("Gene LD-score index component chromosome identity is invalid.")
    baseline_rows = pd.read_parquet(common_path / "baseline_rows.parquet")
    if len(baseline_rows) != int(common_meta.get("n_rows", -1)):
        raise LDSCInputError("Gene LD-score index baseline row count disagrees with component metadata.")
    baseline_columns = list(common_meta.get("baseline_columns", ()))
    allele_columns = [column for column in ("A1", "A2") if column in baseline_rows.columns]
    if allele_columns not in ([], ["A1", "A2"]):
        raise LDSCInputError("Gene LD-score index baseline rows have incomplete allele identity columns.")
    required_baseline = [
        "CHR", "SNP", "POS", *allele_columns, "regression_ld_scores", *baseline_columns
    ]
    if list(baseline_rows.columns) != required_baseline:
        raise LDSCInputError("Gene LD-score index baseline rows have an invalid schema.")
    if any(baseline_rows[column].dtype != np.float32 for column in ["regression_ld_scores", *baseline_columns]):
        raise LDSCInputError("Gene LD-score index baseline LD-score columns must use float32 dtype.")
    if baseline_rows["CHR"].astype(str).ne(str(chrom)).any():
        raise LDSCInputError("Gene LD-score index baseline rows contain the wrong chromosome.")
    if baseline_rows["POS"].dtype != np.int64:
        raise LDSCInputError("Gene LD-score index baseline POS must use int64 dtype.")
    if baseline_rows.duplicated(["CHR", "POS", "SNP"]).any():
        raise LDSCInputError("Gene LD-score index baseline rows contain duplicate identities.")
    canonical_rows = baseline_rows.sort_values(["POS", "SNP"], kind="mergesort").reset_index(drop=True)
    if not baseline_rows.reset_index(drop=True).equals(canonical_rows):
        raise LDSCInputError("Gene LD-score index baseline rows are not in canonical genomic order.")
    stats = _load_npz_members(
        common_path / "baseline_statistics.npz",
        {
            "baseline_count_all", "baseline_count_common", "baseline_overlap_all",
            "baseline_overlap_common", "total_reference_snps_all", "total_reference_snps_common",
        },
    )
    b = int(common_meta.get("n_baseline", -1))
    for name in ("baseline_count_all", "baseline_count_common"):
        if stats[name].shape != (b,) or stats[name].dtype != np.float64 or not np.isfinite(stats[name]).all():
            raise LDSCInputError(f"Gene LD-score index {name} has invalid shape, dtype, or values.")
    for name in ("baseline_overlap_all", "baseline_overlap_common"):
        if stats[name].shape != (b, b) or stats[name].dtype != np.float64 or not np.isfinite(stats[name]).all():
            raise LDSCInputError(f"Gene LD-score index {name} has invalid shape, dtype, or values.")
    atoms = pd.read_parquet(component_path / "atoms.parquet")
    if list(atoms.columns) != ["atom_id", "CHR", "start0", "end"]:
        raise LDSCInputError("Gene LD-score index atoms.parquet schema is invalid.")
    if any(atoms[column].dtype != np.int64 for column in ("atom_id", "start0", "end")):
        raise LDSCInputError("Gene LD-score index atom identifiers and coordinates must use int64 dtype.")
    if atoms["CHR"].astype(str).ne(str(chrom)).any():
        raise LDSCInputError("Gene LD-score index atoms contain the wrong chromosome.")
    if not np.array_equal(atoms["atom_id"].to_numpy(), np.arange(len(atoms), dtype=np.int64)):
        raise LDSCInputError("Gene LD-score index atom_id ordering is invalid.")
    try:
        gene_to_atom = sparse.load_npz(component_path / "gene_to_atom.npz").tocsr()
        operator = sparse.load_npz(component_path / "ldscore_operator.npz").tocsr()
    except (OSError, ValueError) as exc:
        raise LDSCInputError("Gene LD-score index sparse payload is unreadable.") from exc
    model = ChromosomeAtomModel(
        chromosome=str(chrom),
        starts=atoms["start0"].to_numpy(dtype=np.int64),
        ends=atoms["end"].to_numpy(dtype=np.int64),
        gene_to_atom=gene_to_atom,
    )
    n_atoms = int(profile_meta.get("n_atoms", -1))
    n_rows = int(profile_meta.get("n_rows", -1))
    if (
        model.n_atoms != n_atoms
        or model.gene_to_atom.shape[0] != int(profile_meta.get("n_genes", -1))
        or model.gene_to_atom.shape[0] != expected_gene_rows
    ):
        raise LDSCInputError("Gene LD-score index atom dimensions disagree with metadata.")
    if n_rows != len(baseline_rows):
        raise LDSCInputError("Gene LD-score index operator row count disagrees with baseline rows.")
    try:
        model.validate()
        validate_ldscore_operator(operator, n_rows=n_rows, n_atoms=n_atoms)
    except LDSCInternalError as exc:
        raise LDSCInputError(f"Gene LD-score index sparse structure is invalid: {exc}") from exc
    atom_stats = _load_npz_members(
        component_path / "atom_statistics.npz",
        {"atom_count_all", "atom_count_common", "baseline_atom_overlap_all", "baseline_atom_overlap_common"},
    )
    for name in ("atom_count_all", "atom_count_common"):
        if atom_stats[name].shape != (n_atoms,) or atom_stats[name].dtype != np.int64:
            raise LDSCInputError(f"Gene LD-score index {name} has invalid shape or dtype.")
    for name in ("baseline_atom_overlap_all", "baseline_atom_overlap_common"):
        if atom_stats[name].shape != (b, n_atoms) or atom_stats[name].dtype != np.float64 or not np.isfinite(atom_stats[name]).all():
            raise LDSCInputError(f"Gene LD-score index {name} has invalid shape, dtype, or values.")
    for name in ("total_reference_snps_all", "total_reference_snps_common"):
        if stats[name].shape != () or stats[name].dtype != np.int64:
            raise LDSCInputError(f"Gene LD-score index {name} must be a scalar int64 value.")
    return IndexChromosomeData(
        baseline_rows=baseline_rows,
        baseline_count_all=stats["baseline_count_all"],
        baseline_count_common=stats["baseline_count_common"],
        baseline_overlap_all=stats["baseline_overlap_all"],
        baseline_overlap_common=stats["baseline_overlap_common"],
        total_reference_snps_all=int(stats["total_reference_snps_all"]),
        total_reference_snps_common=int(stats["total_reference_snps_common"]),
        atom_model=model,
        operator=operator,
        atom_statistics=AtomStatistics(
            atom_stats["atom_count_all"],
            atom_stats["atom_count_common"],
            atom_stats["baseline_atom_overlap_all"],
            atom_stats["baseline_atom_overlap_common"],
        ),
        reference_metadata=pd.DataFrame(),
    )
