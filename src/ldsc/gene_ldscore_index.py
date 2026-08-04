"""Build, validate, and query exact disjoint-atom gene LD-score indexes.

The offline workflow computes a PLINK-backed LD-score operator for one fixed
baseline, reference panel, regression-row policy, and gene projection. The online
workflow resolves gene lists against the embedded catalog, assembles their
Boolean union from the stored operator, and writes an ordinary self-contained
LD-score directory. It never discovers an index or silently falls back to live
LD calculation.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
from contextlib import contextmanager
from dataclasses import dataclass
from datetime import datetime, timezone
import fcntl
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
from ._logging import log_inputs, log_outputs, workflow_logging
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


_IDENTITY_COLUMNS = ["CHR", "POS", "SNP"]


def _semantic_sha256(payload: dict) -> str:
    """Hash one scientific identity with canonical JSON normalization."""
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def calculate_index_id(index_identity: dict) -> str:
    """Return the immutable semantic identity of one complete index."""
    return _semantic_sha256(index_identity)


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
    parser.add_argument(
        "--regression-snps-file",
        default=None,
        help=(
            "Optional identity-only SNP list defining the persisted regression/output rows. "
            "Duplicate restriction keys collapse and non-identity columns are ignored."
        ),
    )
    parser.add_argument(
        "--exclude-regions",
        choices=kernel_regions.EXCLUDE_REGIONS_CHOICES,
        default="mhc-and-centromeres",
        help=(
            "Curated regions subtracted after selecting bundled HapMap3 or "
            "--regression-snps-file SNPs."
        ),
    )
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
    """Build and transactionally publish one complete gene LD-score index.

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
        Rebuild and replace a complete existing valid index. Gene indexes do
        not support incremental updates or component reuse. Default is
        ``False``.
    log_level : {"DEBUG", "INFO", "WARNING", "ERROR"}, optional
        Workflow logging threshold. The live log is written to the sibling
        ``<output_dir>.build/build-gene-ldscore-index.log`` path so it is never
        part of the replaceable index transaction. Default is ``"INFO"``.

    Returns
    -------
    pathlib.Path
        Published index directory.

    Raises
    ------
    TypeError
        If ``config`` is not a ``GeneLDScoreIndexBuildConfig``.
    LDSCInputError
        If an input is missing or incompatible, the baseline/PLINK identifier
        intersection is empty or ambiguous, or staged validation fails.
    FileExistsError
        If a valid index exists and ``overwrite`` is false, or a nonempty
        unrecognized directory occupies the destination.

    Notes
    -----
    Construction uses float64 adjusted-r-squared accumulation, retains
    negative values, and batches internal atom columns. ``threads`` controls
    chromosome workers and can multiply chromosome-local peak memory. A missing
    or empty destination is not mutated before commit. After destination reload
    validation succeeds, failure to remove builder-owned staging/backup data is
    warned with its retained path and does not change publication success.
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
            regression_snps_file=config.regression_snps_file,
            exclude_regions=config.exclude_regions,
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
    """Build the index's self-contained ordered catalog and inclusion policy."""
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
    """Build, validate, and publish one complete exact v1 index."""
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
        regression_snps_file=getattr(args, "regression_snps_file", None),
        exclude_regions=getattr(args, "exclude_regions", "mhc-and-centromeres"),
        snp_batch_size=args.snp_batch_size,
        atom_batch_size=args.atom_batch_size,
        threads=args.threads,
    )

    index_path = Path(config.output_dir)
    with _gene_index_build_lock(index_path):
        live_log_path = _prepare_gene_index_log(index_path)
        with workflow_logging("build-gene-ldscore-index", live_log_path, log_level=args.log_level):
            _recover_gene_index_publication(index_path)
            _preflight_gene_index_output(index_path, overwrite=bool(args.overwrite))
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
                regression_snps_file=config.regression_snps_file or "bundled HapMap3",
                exclude_regions=config.exclude_regions,
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
            )


def _run_gene_ldscore_index_build(
    args: argparse.Namespace,
    config: GeneLDScoreIndexBuildConfig,
    chromosomes: tuple[str, ...],
    *,
    started: float,
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
    regression_path = Path(config.regression_snps_file or packaged_hm3_curated_map_path())
    regression_keys = read_snp_restriction_keys(regression_path, "rsid", genome_build=None)
    regression_presets = kernel_regions.exclude_regions_choice_to_presets(config.exclude_regions)
    regression_regions = (
        None
        if not regression_presets
        else kernel_regions.load_preset_intervals(regression_presets, config.genome_build)
    )
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
            phase = "baseline/PLINK identifier intersection"
            intersection = intersect_baseline_plink_by_identifier(
                metadata,
                baseline,
                bim_rows,
                snp_identifier=config.snp_identifier,
                chrom=chrom,
            )
            kernel_bundle = kernel_ldscore.AnnotationBundle(
                metadata=intersection.metadata,
                annotations=intersection.annotations,
                baseline_columns=list(public_bundle.baseline_columns),
                query_columns=[],
            )
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
            **intersection.diagnostics,
            "baseline_content_sha256": _canonical_frame_sha256(
                pd.concat([metadata.reset_index(drop=True), baseline.reset_index(drop=True)], axis=1),
                sort_by=("CHR", "POS", "SNP"),
            ),
            "reference_rows_before_genotype_qc": rows_before_genotype_qc,
            "annotation_intersection_removed": int(len(bim_rows) - rows_before_genotype_qc),
            "retained_reference_rows": int(record.total_reference_snps_all),
            "retained_common_reference_rows": int(record.total_reference_snps_common),
            "genotype_qc_removed": genotype_qc_removed,
            "maf_removed": maf_removed,
            "genotype_qc_or_maf_removed": genotype_qc_removed + maf_removed,
            "regression_rows": int(len(record.baseline_rows)),
            "baseline_plink_identity": "inner_join_by_rsid",
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
        if not len(record.baseline_rows):
            LOGGER.warning(
                f"Chromosome {chrom} has zero persisted regression SNP rows after restriction and region exclusion."
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
    total_regression_rows = sum(len(record.baseline_rows) for record in chromosome_data.values())
    if total_regression_rows == 0:
        raise LDSCInputError(
            "Gene LD-score index has zero regression SNP rows across all selected chromosomes "
            "after applying the SNP restriction and --exclude-regions policy."
        )
    index_identity = _builder_index_identity(
        config,
        args,
        chromosome_data,
        regression_keys,
        genetic_map,
        catalog,
        evidence_by_chrom,
    )
    index_path = Path(config.output_dir)
    index_id = calculate_index_id(index_identity)
    elapsed = time.perf_counter() - started
    import resource

    peak_rss = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    if sys.platform != "darwin":
        peak_rss *= 1024
    selected_individuals = index_identity.get("selected_individuals", {})
    plink_sources = index_identity.get("plink_sources", [])
    map_identity = index_identity.get("genetic_map", "bim_cm")
    map_source = "explicit hg19 map" if isinstance(map_identity, dict) else "BIM cM fallback"
    diagnostic_payload = {
        "index_id": index_id,
        "elapsed_seconds": elapsed,
        "peak_rss_bytes": peak_rss,
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
            "regression_snps": "custom" if config.regression_snps_file else "bundled_hapmap3",
            "exclude_regions": config.exclude_regions,
        },
        "input_resolution": {
            "plink_bed_files": sum(source.get("kind") == "bed" for source in plink_sources),
            "plink_bim_files": sum(source.get("kind") == "bim" for source in plink_sources),
            "chromosome_shards": len(chromosomes),
            "selected_individual_count": selected_individuals.get("selected_count"),
            "keep_individual_source": "all" if config.keep_indivs_file is None else "keep_file",
            "baseline_plink_intersection": {
                chrom: {
                    key: values[key]
                    for key in (
                        "baseline_rows", "plink_rows", "matched_rows", "baseline_only_rows",
                        "plink_only_rows", "coordinate_discordant_rows",
                    )
                }
                for chrom, values in evidence_by_chrom.items()
            },
        },
        "publication": {
            "overwrite": bool(args.overwrite),
            "staging_reload_validation": "passed",
            "atomic_publication": "completed",
            "index_path": str(index_path),
        },
        "chromosomes": evidence_by_chrom,
        "totals": {
            "chromosome_count": len(chromosomes),
            "protein_coding_genes": sum(values["protein_coding_genes"] for values in evidence_by_chrom.values()),
            "atom_count": sum(values["atom_count"] for values in evidence_by_chrom.values()),
            "operator_nnz": sum(values["operator_nnz"] for values in evidence_by_chrom.values()),
            "regression_rows": total_regression_rows,
            "retained_reference_rows": sum(values["retained_reference_rows"] for values in evidence_by_chrom.values()),
        },
    }
    try:
        published_path = publish_gene_ldscore_index(
            config.output_dir,
            index_identity=index_identity,
            gene_catalog=embedded_catalog,
            chromosomes=chromosome_data,
            overwrite=bool(args.overwrite),
            diagnostic_payload=diagnostic_payload,
        )
    except Exception as exc:
        LOGGER.error(
            f"Gene LD-score index build failed during staged validation/publication: "
            f"{type(exc).__name__}: {exc}"
        )
        raise
    _log_gene_index_summary(diagnostic_payload)
    log_outputs(
        index_id=index_id,
        index=str(published_path),
        summary_json=str(published_path / "diagnostics" / "build-gene-ldscore-index.json"),
        publication="staged, reloaded, validated, and atomically replaced",
    )
    return published_path


@contextmanager
def _gene_index_build_lock(index_path: Path):
    """Hold a nonblocking process lock for one absolute index destination."""
    destination = index_path.expanduser().resolve()
    destination.parent.mkdir(parents=True, exist_ok=True)
    lock_path = destination.parent / f".{destination.name}.build-gene-ldscore-index.lock"
    handle = lock_path.open("a+", encoding="utf-8")
    try:
        try:
            fcntl.flock(handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            handle.seek(0)
            owner = handle.read().strip() or "unknown owner"
            raise LDSCInputError(
                f"Another build-gene-ldscore-index process is active for {destination} ({owner}). "
                f"Monitor {_gene_index_build_state_dir(destination) / 'build-gene-ldscore-index.log'}."
            ) from exc
        handle.seek(0)
        handle.truncate()
        handle.write(f"pid={os.getpid()} host={socket.gethostname()} target={destination}\n")
        handle.flush()
        yield
    finally:
        try:
            fcntl.flock(handle.fileno(), fcntl.LOCK_UN)
        finally:
            handle.close()


def _preflight_gene_index_output(index_path: Path, *, overwrite: bool) -> None:
    """Validate the publication destination without creating or mutating it."""
    if index_path.exists() and not index_path.is_dir():
        raise FileExistsError(f"Gene LD-score index output is not a directory: {index_path}")
    if not index_path.exists():
        return
    if not any(index_path.iterdir()) or _is_diagnostics_only_gene_index(index_path):
        return
    try:
        load_gene_ldscore_index(index_path)
    except Exception as exc:
        raise FileExistsError(
            f"Gene LD-score index output directory is nonempty but invalid: {index_path}. "
            "Choose an empty directory; --overwrite does not replace unrecognized contents."
        ) from exc
    if not overwrite:
        raise FileExistsError(
            f"Gene LD-score index already exists: {index_path}. Pass --overwrite to rebuild it completely."
        )


def _recover_gene_index_publication(index_path: Path) -> None:
    """Recover or clean recognized interrupted publication transactions."""
    destination = index_path.expanduser().resolve()
    candidates: list[Path] = []
    for path in destination.parent.glob(f".{destination.name}.stage-*"):
        marker_path = path / ".gene-index-publication.json"
        try:
            marker = _read_json(marker_path, "publication marker")
        except LDSCInputError:
            continue
        if (
            marker.get("artifact_type") == "gene_ldscore_index_publication"
            and marker.get("target") == str(destination)
        ):
            candidates.append(path)
    if not candidates:
        return

    target_valid = False
    if destination.exists():
        try:
            load_gene_ldscore_index(destination)
        except Exception:
            target_valid = False
        else:
            target_valid = True
    if target_valid:
        try:
            (destination / ".gene-index-publication.json").unlink(missing_ok=True)
        except Exception as exc:
            LOGGER.warning(
                "The published gene LD-score index is valid, but its completed-publication marker "
                f"could not be removed: {destination / '.gene-index-publication.json'} "
                f"({type(exc).__name__}: {exc}). The marker is safe to remove after this build."
            )
        for candidate in candidates:
            _remove_gene_index_transaction(candidate, published=True)
        return

    valid_backups: list[Path] = []
    for candidate in candidates:
        backup = candidate / f"{destination.name}.backup"
        if not backup.is_dir():
            continue
        marker_path = backup / ".gene-index-publication.json"
        try:
            marker = _read_json(marker_path, "backup publication marker")
            load_gene_ldscore_index(backup)
        except Exception:
            continue
        if marker.get("target") == str(destination):
            valid_backups.append(backup)
    if len(valid_backups) != 1:
        listed = ", ".join(str(path) for path in valid_backups) or "none"
        raise LDSCInputError(
            "Gene LD-score index publication recovery is ambiguous: the target is missing or invalid "
            f"and valid matching backups are {listed}. Preserve these paths and resolve them manually."
        )
    if destination.exists():
        target_marker = destination / ".gene-index-publication.json"
        try:
            marker = _read_json(target_marker, "target publication marker")
        except LDSCInputError as exc:
            raise LDSCInputError(
                f"Cannot replace invalid unrecognized gene-index target during recovery: {destination}."
            ) from exc
        if marker.get("target") != str(destination):
            raise LDSCInputError(
                f"Cannot replace invalid unrecognized gene-index target during recovery: {destination}."
            )
        shutil.rmtree(destination)
    backup = valid_backups[0]
    os.replace(backup, destination)
    try:
        (destination / ".gene-index-publication.json").unlink(missing_ok=True)
    except Exception as exc:
        LOGGER.warning(
            "The recovered gene LD-score index is valid, but its publication marker could not be "
            f"removed: {destination / '.gene-index-publication.json'} "
            f"({type(exc).__name__}: {exc})."
        )
    for candidate in candidates:
        if candidate.exists():
            _remove_gene_index_transaction(candidate, published=True)


def _gene_index_build_state_dir(index_path: Path) -> Path:
    """Return the stable sibling directory for mutable build diagnostics."""
    destination = index_path.expanduser()
    return destination.with_name(f"{destination.name}.build")


def _prepare_gene_index_log(index_path: Path) -> Path:
    """Archive the prior sidecar log and return the stable live log path."""
    build_state = _gene_index_build_state_dir(index_path)
    build_state.mkdir(parents=True, exist_ok=True)
    log_path = build_state / "build-gene-ldscore-index.log"
    if log_path.exists():
        timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S.%fZ")
        history = build_state / "history"
        history.mkdir(exist_ok=True)
        os.replace(log_path, history / f"build-gene-ldscore-index.{timestamp}.log")
    return log_path


def _remove_gene_index_transaction(path: Path, *, published: bool) -> bool:
    """Best-effort removal for builder-owned staging and backup data.

    Cleanup after a reload-validated destination is garbage collection and
    cannot change the successful publication outcome. Before commit, cleanup
    errors likewise must not replace the primary build or rollback exception.
    """
    if not path.exists():
        return True
    try:
        shutil.rmtree(path)
    except Exception as exc:
        if published:
            LOGGER.warning(
                "The gene LD-score index was published and reload-validated, but transaction "
                f"cleanup failed ({type(exc).__name__}: {exc}). The published index is valid. "
                f"Retained builder-owned cleanup path: {path}. Remove it after confirming no "
                "build is active for this destination."
            )
        else:
            LOGGER.warning(
                "Gene LD-score index transaction cleanup also failed after an earlier build or "
                f"publication error ({type(exc).__name__}: {exc}). Retained builder-owned path: "
                f"{path}. The earlier exception remains the build failure."
            )
        return False
    return True


def _directory_payload_bytes(path: Path) -> int:
    """Return the byte size of one scientific artifact component."""
    if not path.exists():
        return 0
    return sum(file.stat().st_size for file in path.rglob("*") if file.is_file())


def _scientific_payload_bytes(index_path: Path) -> int:
    """Return index bytes excluding human-readable and JSON diagnostics."""
    total = 0
    for file in index_path.rglob("*"):
        if not file.is_file() or "diagnostics" in file.parts:
            continue
        total += file.stat().st_size
    return total


def _is_diagnostics_only_gene_index(path: Path) -> bool:
    """Recognize a legacy failed-build directory containing only owned diagnostics."""
    if not path.is_dir() or (path / "metadata.json").exists():
        return False
    files = [file for file in path.rglob("*") if file.is_file()]
    if not files:
        return False
    return all(
        file.relative_to(path).parts[0] == "diagnostics"
        and file.name.startswith("build-gene-ldscore-index")
        for file in files
    )


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
        f"Resolved inputs: PLINK bed/bim={inputs['plink_bed_files']}/{inputs['plink_bim_files']}, "
        f"chromosome-shards={inputs['chromosome_shards']}, "
        f"selected-individuals={inputs['selected_individual_count']} "
        f"({inputs['keep_individual_source']})."
    )
    LOGGER.info(
        "SNP universe: the baseline/PLINK identifier intersection defines LD-score contributors and "
        "count/overlap members; the configured regression restriction defines persisted rows."
    )
    LOGGER.info(
        f"Regression-row policy: {config['regression_snps']} with exclude-regions={config['exclude_regions']}."
    )
    for chrom in payload["chromosomes"]:
        values = payload["chromosomes"][chrom]
        LOGGER.info(
            f"Chromosome {chrom} metrics: pre-QC={values['pre_qc_rows']}, "
            f"baseline/PLINK-match={values['matched_rows']}, "
            f"baseline-only={values['baseline_only_rows']}, PLINK-only={values['plink_only_rows']}, "
            f"coordinate-discordant={values['coordinate_discordant_rows']}, "
            f"genotype-QC-removed={values['genotype_qc_removed']}, "
            f"MAF-removed={values['maf_removed']} (policy={values['maf_filter_policy']}), "
            f"retained-reference={values['retained_reference_rows']}, "
            f"common-reference={values['retained_common_reference_rows']}, "
            f"regression-rows={values['regression_rows']}, "
            f"protein-coding genes={values['protein_coding_genes']}, "
            f"genes-with-padded-atoms={values['genes_with_padded_atoms']}, "
            f"atoms={values['atom_count']}, nnz(Y)={values['nnz_Y']}, "
            f"operator_nnz={values['operator_nnz']}, "
            f"payload-bytes={values['payload_bytes']}, "
            f"elapsed={values['elapsed_seconds']:.3f}s."
        )
    LOGGER.info(
        f"Publication: index_id={payload['index_id']}, overwrite={publication['overwrite']}, "
        f"staged-reload={publication['staging_reload_validation']}, "
        f"atomic-publication={publication['atomic_publication']}."
    )
    LOGGER.info(
        f"Gene LD-score index build completed and validated: chromosomes={totals['chromosome_count']}, "
        f"protein-coding genes={totals['protein_coding_genes']}, atoms={totals['atom_count']}, "
        f"operator_nnz={totals['operator_nnz']}, payload_bytes={payload['payload_bytes']}, "
        f"peak_rss_bytes={payload['peak_rss_bytes']}, elapsed_seconds={payload['elapsed_seconds']:.6f}, "
        f"index={publication['index_path']}."
    )


def _read_bim_identity(prefix: str) -> pd.DataFrame:
    frame = pd.read_csv(
        prefix + ".bim",
        sep=r"\s+",
        header=None,
        names=["CHR", "SNP", "CM", "POS", "A1", "A2"],
    )
    return frame[["CHR", "SNP", "CM", "POS", "A1", "A2"]]


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


def _canonical_frame_sha256(
    frame: pd.DataFrame,
    *,
    sort_by: Sequence[str] = (),
) -> str:
    """Hash a normalized table without path, compression, or text-layout effects."""
    normalized = frame.copy()
    if sort_by:
        normalized = normalized.sort_values(list(sort_by), kind="mergesort").reset_index(drop=True)
    normalized.columns = [str(column) for column in normalized.columns]
    for column in normalized.columns:
        if pd.api.types.is_float_dtype(normalized[column]):
            normalized[column] = normalized[column].map(
                lambda value: "<NA>" if pd.isna(value) else format(float(value), ".17g")
            )
        elif pd.api.types.is_bool_dtype(normalized[column]):
            normalized[column] = normalized[column].astype(bool).map({True: "true", False: "false"})
        else:
            normalized[column] = normalized[column].map(
                lambda value: "<NA>" if pd.isna(value) else str(value)
            )
    payload = normalized.to_csv(index=False, sep="\t", lineterminator="\n")
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def _genetic_map_identity(args, genetic_map: pd.DataFrame | None = None) -> str | dict:
    """Return the content-bound scientific identity of the effective cM map."""
    sources = split_cli_path_tokens(getattr(args, "genetic_map_hg19_sources", None))
    if not sources:
        return "bim_cm"
    if genetic_map is None:
        genetic_map = _load_builder_genetic_map(args)
    return {
        "kind": "explicit_hg19",
        "content_sha256": _canonical_frame_sha256(genetic_map, sort_by=("CHR", "POS")),
    }


def _builder_index_identity(
    config,
    args,
    chromosomes,
    regression_keys: set[str] | RestrictionIdentityKeys,
    genetic_map: pd.DataFrame | None,
    catalog: GeneCatalog,
    evidence_by_chrom: dict[str, dict],
) -> dict:
    """Build the path-insensitive canonical scientific identity."""
    plink_files = []
    first_prefix = None
    for chrom in chromosomes:
        kernel_args = argparse.Namespace(bfile=config.plink_prefix)
        prefix = kernel_ldscore.resolve_bfile_prefix(kernel_args, chrom=chrom)
        if first_prefix is None:
            first_prefix = prefix
        plink_files.extend(
            [
                {"chromosome": chrom, "kind": "bed", "content_sha256": _file_sha256(str(prefix) + ".bed")},
                {
                    "chromosome": chrom,
                    "kind": "bim",
                    "content_sha256": _canonical_frame_sha256(_read_bim_identity(str(prefix))),
                },
            ]
        )
    assert first_prefix is not None
    fam = kernel_ldscore.legacy_parse.PlinkFAMFile(str(first_prefix) + ".fam")
    keep_indices = kernel_ldscore.resolve_keep_individuals(config.keep_indivs_file, fam)
    selected_indices = list(range(len(fam.IDList))) if keep_indices is None else list(keep_indices)
    selected_ids = _selected_individual_ids(fam, selected_indices)
    selected_digest = hashlib.sha256(
        json.dumps(selected_ids, separators=(",", ":"), ensure_ascii=True).encode("utf-8")
    ).hexdigest()
    effective_keys = regression_keys.keys if isinstance(regression_keys, RestrictionIdentityKeys) else regression_keys
    regression_digest = hashlib.sha256(
        json.dumps(sorted(map(str, effective_keys)), separators=(",", ":")).encode("utf-8")
    ).hexdigest()
    return {
        "baseline_by_chromosome": {
            chrom: evidence_by_chrom[chrom]["baseline_content_sha256"]
            for chrom in chromosomes
        },
        "plink_sources": plink_files,
        "chromosomes": list(config.chromosomes),
        "genome_build": config.genome_build,
        "snp_identifier": config.snp_identifier,
        "selected_individuals": {
            "selected_content_sha256": selected_digest,
            "selected_count": len(selected_ids),
        },
        "maf_min": config.maf_min,
        "common_maf_min": config.common_maf_min,
        "ld_window": {"unit": "cm", "value": config.ld_wind_cm},
        "genetic_map": _genetic_map_identity(args, genetic_map),
        "regression_snps": {
            "kind": "custom" if config.regression_snps_file else "bundled_hapmap3",
            "canonical_keys_sha256": regression_digest,
        },
        "exclude_regions": config.exclude_regions,
        "catalog": {
            "content_sha256": _canonical_frame_sha256(catalog.frame, sort_by=("ensgid",)),
        },
        "projection_build": config.genome_build,
        "padding_bp": config.padding_bp,
        "gene_exclude_regions": config.gene_exclude_regions,
    }


def _selected_individual_ids(fam, selected_indices: Sequence[int]) -> list[str]:
    """Return selected PLINK individual IDs in their effective genotype order."""
    identifiers = fam.IDList
    if not isinstance(identifiers, pd.DataFrame) or identifiers.shape[1] != 1:
        raise LDSCInternalError("PLINK FAM individual identity table must contain exactly one IID column.")
    return identifiers.iloc[list(selected_indices), 0].astype(str).tolist()


@dataclass(frozen=True)
class IndexChromosomeData:
    """Validated in-memory scientific payload for one index chromosome."""

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
    """One validated, self-contained gene LD-score index.

    Attributes
    ----------
    index_id : str
        Content-bound scientific identity for the complete index.
    chromosomes : tuple of str
        Ordered immutable chromosome coverage.
    index_identity : dict
        Scientific settings and canonical input-content identities used to
        calculate ``index_id``.
    gene_catalog : pandas.DataFrame
        Embedded catalog with unpadded 0-based half-open gene coordinates,
        inclusion flags, and chromosome-local row ordering.
    index_chromosomes : dict of str to IndexChromosomeData
        Validated chromosome-local baseline rows, sparse operators, atom
        geometry, and sufficient statistics.
    """

    index_id: str
    chromosomes: tuple[str, ...]
    index_identity: dict
    gene_catalog: pd.DataFrame
    index_chromosomes: dict[str, IndexChromosomeData]


def publish_gene_ldscore_index(
    index_dir: str | Path,
    *,
    index_identity: dict,
    gene_catalog: pd.DataFrame,
    chromosomes: dict[str, IndexChromosomeData],
    overwrite: bool,
    diagnostic_payload: dict | None = None,
) -> Path:
    """Stage, reload, and atomically publish one complete index directory.

    Existing valid indexes require ``overwrite=True`` and are kept loadable
    until their complete replacement has passed staged reload validation.
    Nonempty directories that are neither valid indexes nor recognized
    legacy diagnostics-only failed builds are rejected even with overwrite.
    Mutable live logs are stored in a sibling build-state directory and never
    enter the replaceable transaction. Once the destination has been replaced
    and reload-validated, failure to remove the builder-owned transaction tree
    is reported as a warning and does not change publication success.
    """
    destination = Path(index_dir)
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists():
        if not destination.is_dir():
            raise FileExistsError(f"Gene LD-score index output is not a directory: {destination}")
        if not any(destination.iterdir()):
            pass
        elif _is_diagnostics_only_gene_index(destination):
            pass
        else:
            try:
                load_gene_ldscore_index(destination)
            except Exception as exc:
                raise FileExistsError(
                    f"Gene LD-score index output directory is nonempty but invalid: {destination}. "
                    "Choose an empty directory; --overwrite does not replace unrecognized contents."
                ) from exc
            if not overwrite:
                raise FileExistsError(
                    f"Gene LD-score index already exists: {destination}. Pass --overwrite to rebuild it completely."
                )
    index_id = calculate_index_id(index_identity)
    stage_parent = Path(tempfile.mkdtemp(prefix=f".{destination.name}.stage-", dir=destination.parent))
    staged_index = stage_parent / destination.name
    publication_marker = {
        "artifact_type": "gene_ldscore_index_publication",
        "target": str(destination.expanduser().resolve()),
    }
    _write_json(stage_parent / ".gene-index-publication.json", publication_marker)
    committed = False
    try:
        _write_index_artifact(
            staged_index,
            index_id=index_id,
            index_identity=index_identity,
            gene_catalog=gene_catalog,
            chromosomes=chromosomes,
        )
        if diagnostic_payload is not None:
            diagnostic_payload["payload_bytes"] = _scientific_payload_bytes(staged_index)
            chromosome_diagnostics = diagnostic_payload.get("chromosomes", {})
            for chrom in chromosomes:
                values = chromosome_diagnostics.get(chrom)
                if isinstance(values, dict):
                    values["payload_bytes"] = _directory_payload_bytes(
                        staged_index / "chromosomes" / f"chr{chrom}"
                    )
            _write_json(
                staged_index / "diagnostics" / "build-gene-ldscore-index.json",
                diagnostic_payload,
            )
        _write_json(staged_index / ".gene-index-publication.json", publication_marker)
        load_gene_ldscore_index(staged_index)
        backup = stage_parent / f"{destination.name}.backup"
        if destination.exists():
            _write_json(destination / ".gene-index-publication.json", publication_marker)
            try:
                os.replace(destination, backup)
            except Exception:
                (destination / ".gene-index-publication.json").unlink(missing_ok=True)
                raise
        try:
            os.replace(staged_index, destination)
            load_gene_ldscore_index(destination)
        except Exception:
            if destination.exists():
                shutil.rmtree(destination)
            if backup.exists():
                os.replace(backup, destination)
                (destination / ".gene-index-publication.json").unlink(missing_ok=True)
            raise
        committed = True
        try:
            (destination / ".gene-index-publication.json").unlink(missing_ok=True)
        except Exception as exc:
            LOGGER.warning(
                "The gene LD-score index was published and reload-validated, but its publication "
                f"marker could not be removed ({type(exc).__name__}: {exc}): "
                f"{destination / '.gene-index-publication.json'}. The published index is valid."
            )
    except Exception:
        _remove_gene_index_transaction(stage_parent, published=False)
        raise
    _remove_gene_index_transaction(stage_parent, published=committed)
    return destination


def load_gene_ldscore_index(index_dir: str | Path) -> LoadedGeneLDScoreIndex:
    """Load and fully validate one explicit gene LD-score index directory.

    Parameters
    ----------
    index_dir : path-like
        Directory containing one complete index. No profile discovery or
        compatibility migration is performed.

    Returns
    -------
    LoadedGeneLDScoreIndex
        Index components after validating semantic identity,
        chromosome coverage, catalog and row ordering, shapes, dtypes, finite
        values, and CSR invariants.

    Raises
    ------
    LDSCInputError
        If any required component is missing, corrupt, structurally invalid,
        or semantically incompatible.
    """
    index_path = Path(index_dir)
    root = _read_json(index_path / "metadata.json", "root")
    if root.get("artifact_type") != "gene_ldscore_index":
        raise LDSCInputError("Gene LD-score index root has the wrong artifact_type.")
    index_identity = root.get("index_identity")
    index_id = root.get("index_id")
    if not isinstance(index_identity, dict) or index_id != calculate_index_id(index_identity):
        raise LDSCInputError("Gene LD-score index semantic identity is invalid.")
    chromosomes = tuple(str(chrom) for chrom in root.get("chromosomes", ()))
    expected_coverage = tuple(str(chrom) for chrom in index_identity.get("chromosomes", ()))
    if not chromosomes or chromosomes != expected_coverage:
        raise LDSCInputError("Gene LD-score index chromosome coverage is invalid.")
    catalog_path = index_path / "gene_catalog.parquet"
    if not catalog_path.exists():
        raise LDSCInputError("Gene LD-score index is missing gene_catalog.parquet.")
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
            index_path,
            chrom,
            index_id=index_id,
            expected_gene_rows=len(chromosome_catalog),
        )
    return LoadedGeneLDScoreIndex(
        index_id=index_id,
        chromosomes=chromosomes,
        index_identity=index_identity,
        gene_catalog=catalog,
        index_chromosomes=loaded,
    )


def run_indexed_ldscore(
    index_dir: str | Path,
    *,
    query_gene_list_sources: Sequence[str | Path],
    control_gene_list_source: str | Path = "all-protein-coding",
    output_dir: str | Path,
    overwrite: bool = False,
):
    """Assemble exact gene-list LD scores from one explicit index directory.

    Parameters
    ----------
    index_dir : path-like
        Exact validated single-index directory.
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
        If the index is invalid, a control is unusable, or every requested
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

    index = load_gene_ldscore_index(index_dir)
    catalog_identity = dict(index.index_identity.get("catalog", {}))
    projection_build = str(index.index_identity.get("projection_build"))
    gene_policy = str(index.index_identity.get("gene_exclude_regions", "none"))
    catalog = GeneCatalog.from_index_frame(
        index.gene_catalog,
        genome_build=projection_build,
        release=str(catalog_identity.get("release", "GENCODE v49")),
        content_sha256=str(catalog_identity.get("content_sha256", "")),
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

    supplied_baseline = _baseline_columns_from_rows(next(iter(index.index_chromosomes.values())).baseline_rows)
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
        record = index.index_chromosomes[chrom]
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
        index_provenance={"index_id": index.index_id},
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


@dataclass(frozen=True)
class BaselinePlinkIntersection:
    """Baseline annotations aligned to the PLINK-authoritative SNP order."""

    metadata: pd.DataFrame
    annotations: pd.DataFrame
    diagnostics: dict[str, int]


def intersect_baseline_plink_by_identifier(
    baseline_rows: pd.DataFrame,
    baseline_annotations: pd.DataFrame,
    plink_rows: pd.DataFrame,
    *,
    snp_identifier: str,
    chrom: str,
) -> BaselinePlinkIntersection:
    """Inner-join baseline annotations to PLINK by the configured SNP identity.

    PLINK supplies chromosome, position, alleles, cM coordinates, genotype
    order, and therefore the returned metadata order. Baseline-only and
    PLINK-only identities are dropped. Duplicate effective identifiers and an
    empty intersection are rejected because either condition makes scientific
    matching ambiguous or unusable.
    """
    if snp_identifier != "rsid":
        raise LDSCInputError(
            "The v1 gene LD-score index builder supports only rsid identity."
        )
    if len(baseline_rows) != len(baseline_annotations):
        raise LDSCInputError(
            "Baseline metadata and annotation row counts differ before PLINK intersection."
        )
    baseline = _normalize_identity_rows(baseline_rows, label="baseline", chrom=chrom)
    plink = _normalize_identity_rows(plink_rows, label="PLINK BIM", chrom=chrom)
    for label, frame in (("baseline", baseline), ("PLINK BIM", plink)):
        duplicate = frame["SNP"].duplicated(keep=False)
        if duplicate.any():
            raise LDSCInputError(
                f"Gene LD-score index build found duplicate {label} effective rsid identifiers "
                f"on chromosome {chrom}; identifier-key matching would be ambiguous."
            )
    baseline_keys = set(baseline["SNP"])
    plink_keys = set(plink["SNP"])
    matched = plink["SNP"].isin(baseline_keys)
    if not matched.any():
        raise LDSCInputError(
            f"Gene LD-score index baseline/PLINK rsid intersection is empty on chromosome {chrom}."
        )
    baseline_positions = baseline.set_index("SNP")["POS"]
    matched_plink = plink.loc[matched].reset_index(drop=True)
    coordinate_discordant = int(
        np.count_nonzero(
            matched_plink["POS"].to_numpy(dtype=np.int64)
            != matched_plink["SNP"].map(baseline_positions).to_numpy(dtype=np.int64)
        )
    )
    if coordinate_discordant:
        LOGGER.warning(
            f"Baseline/PLINK rsid intersection on chromosome {chrom} has "
            f"{coordinate_discordant} coordinate disagreement row(s); PLINK coordinates are authoritative."
        )
    baseline_index = pd.Series(baseline.index, index=baseline["SNP"])
    annotation_rows = baseline_index.loc[matched_plink["SNP"]].to_numpy(dtype=np.int64)
    aligned_annotations = baseline_annotations.iloc[annotation_rows].reset_index(drop=True)
    diagnostics = {
        "baseline_rows": int(len(baseline)),
        "plink_rows": int(len(plink)),
        "matched_rows": int(len(matched_plink)),
        "baseline_only_rows": int(len(baseline_keys - plink_keys)),
        "plink_only_rows": int(len(plink_keys - baseline_keys)),
        "coordinate_discordant_rows": coordinate_discordant,
    }
    return BaselinePlinkIntersection(matched_plink, aligned_annotations, diagnostics)


def _normalize_identity_rows(frame: pd.DataFrame, *, label: str, chrom: str) -> pd.DataFrame:
    """Normalize one strict-identity frame without dropping caller columns."""
    missing = set(_IDENTITY_COLUMNS) - set(frame.columns)
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


def _validate_metadata_identity(metadata: dict, *, index_id: str) -> None:
    if metadata.get("artifact_type") != "gene_ldscore_index" or metadata.get("index_id") != index_id:
        raise LDSCInputError("Gene LD-score index component identity is invalid.")


def _cast_index_baseline_rows(frame: pd.DataFrame) -> pd.DataFrame:
    """Match the canonical LD-score writer's persisted floating dtypes."""
    float64_columns = [column for column in frame.columns if frame[column].dtype == np.float64]
    return frame.astype({column: np.float32 for column in float64_columns}) if float64_columns else frame


def _write_index_artifact(
    index_path: Path,
    *,
    index_id: str,
    index_identity: dict,
    gene_catalog: pd.DataFrame,
    chromosomes: dict[str, IndexChromosomeData],
) -> None:
    """Write one complete staged index in the canonical single-directory layout."""
    index_path.mkdir(parents=True)
    _write_json(
        index_path / "metadata.json",
        {
            "artifact_type": "gene_ldscore_index",
            "index_id": index_id,
            "index_identity": index_identity,
            "chromosomes": list(chromosomes),
        },
    )
    gene_catalog.to_parquet(index_path / "gene_catalog.parquet", index=False)
    for chrom, record in chromosomes.items():
        chrom_path = index_path / "chromosomes" / f"chr{chrom}"
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
                "index_id": index_id,
                "chromosome": str(chrom),
                "n_rows": len(record.baseline_rows),
                "n_baseline": len(baseline_columns),
                "baseline_columns": baseline_columns,
                "n_genes": record.atom_model.gene_to_atom.shape[0],
                "n_atoms": record.atom_model.n_atoms,
                "gene_to_atom_format": "csr_bool_int32",
                "ldscore_operator_format": "csr_float64_int32",
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
        atoms = pd.DataFrame(
            {
                "atom_id": np.arange(record.atom_model.n_atoms, dtype=np.int64),
                "CHR": str(chrom),
                "start0": record.atom_model.starts,
                "end": record.atom_model.ends,
            }
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
    index_path: Path,
    chrom: str,
    *,
    index_id: str,
    expected_gene_rows: int,
) -> IndexChromosomeData:
    component_path = index_path / "chromosomes" / f"chr{chrom}"
    component_meta = _read_json(component_path / "metadata.json", f"chromosome {chrom}")
    _validate_metadata_identity(component_meta, index_id=index_id)
    if component_meta.get("chromosome") != str(chrom):
        raise LDSCInputError("Gene LD-score index component chromosome identity is invalid.")
    baseline_rows = pd.read_parquet(component_path / "baseline_rows.parquet")
    if len(baseline_rows) != int(component_meta.get("n_rows", -1)):
        raise LDSCInputError("Gene LD-score index baseline row count disagrees with component metadata.")
    baseline_columns = list(component_meta.get("baseline_columns", ()))
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
        component_path / "baseline_statistics.npz",
        {
            "baseline_count_all", "baseline_count_common", "baseline_overlap_all",
            "baseline_overlap_common", "total_reference_snps_all", "total_reference_snps_common",
        },
    )
    b = int(component_meta.get("n_baseline", -1))
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
    n_atoms = int(component_meta.get("n_atoms", -1))
    n_rows = int(component_meta.get("n_rows", -1))
    if (
        model.n_atoms != n_atoms
        or model.gene_to_atom.shape[0] != int(component_meta.get("n_genes", -1))
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
