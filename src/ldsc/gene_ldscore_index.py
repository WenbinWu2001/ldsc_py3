"""Build, validate, and query exact disjoint-atom gene LD-score indexes.

The offline workflow computes a PLINK-backed LD-score operator for one fixed
baseline, reference panel, regression-row policy, and gene projection. The online
workflow resolves gene lists against the embedded catalog, assembles their
Boolean union from the stored operator, and writes an ordinary self-contained
LD-score directory. It never discovers an index or silently falls back to live
LD calculation. During an offline build, completed chromosome payloads become
durable inside one private run transaction and are released from memory; only a
complete reload-validated transaction is published as an index.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
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
from types import SimpleNamespace
import warnings
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
from ._kernel.snp_identity import (
    RestrictionIdentityKeys,
    clean_identity_artifact_table,
    coerce_identity_drop_frame,
    effective_merge_key_series,
    empty_identity_drop_frame,
    normalize_snp_identifier_mode,
)
from ._kernel import regions as kernel_regions
from ._logging import (
    log_inputs,
    log_outputs,
    materializing_overwrite_guard,
    set_workflow_log_path,
    workflow_logging,
)
from .annotation_builder import AnnotationBuilder
from .chromosome_inference import normalize_chromosome, normalize_chromosome_series
from ._coordinates import positive_int_position_series
from .config import AnnotationBuildConfig, GeneLDScoreIndexBuildConfig, GlobalConfig
from .errors import LDSCInputError, LDSCInternalError
from .gene_list_resolver import (
    GeneCatalog,
    GeneCatalogValidationError,
    GeneSourceSelection,
    resolve_gene_lists,
    select_index_eligible_gene_indices,
)
from .query_annotations import QueryAnnotationStatus
from .hm3 import packaged_hm3_curated_map_path
from .path_resolution import split_cli_path_tokens


LOGGER = logging.getLogger("LDSC.gene_ldscore_index")


_IDENTITY_COLUMNS = ["CHR", "POS", "SNP"]


class _GeneIndexArgumentParser(argparse.ArgumentParser):
    """Argument parser with the builder's deliberate identity omissions."""

    def error(self, message: str) -> None:
        missing: list[str] = []
        if "required" in message and "--snp-identifier" in message:
            missing.append("--snp-identifier is required; choose rsid or chr_pos.")
        if "required" in message and "--genome-build" in message:
            missing.append("--genome-build is required; choose hg19.")
        super().error(" ".join(missing) if missing else message)


def _semantic_sha256(payload: dict) -> str:
    """Hash one scientific identity with canonical JSON normalization."""
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def calculate_index_id(index_identity: dict) -> str:
    """Return the immutable semantic identity of one complete index."""
    return _semantic_sha256(index_identity)


def build_parser() -> argparse.ArgumentParser:
    """Build the closed-v1 offline gene LD-score index parser."""
    parser = _GeneIndexArgumentParser(
        prog="ldsc build-gene-ldscore-index",
        description="Build an exact PLINK-backed disjoint-atom gene LD-score index.",
        allow_abbrev=False,
    )
    parser.add_argument("--baseline-annot-sources", required=True)
    parser.add_argument("--plink-prefix", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument(
        "--gene-coordinate-file",
        required=True,
        help="Canonical one-build gene-coordinate TSV/TSV.GZ validated before atom construction.",
    )
    parser.add_argument(
        "--genome-build",
        choices=("hg19",),
        required=True,
        help="Required explicit build assertion for all coordinate-bearing inputs; no inference or liftover.",
    )
    parser.add_argument(
        "--snp-identifier",
        choices=("rsid", "chr_pos"),
        required=True,
        help="Required effective identity mode; no default, auto mode, or column-based inference.",
    )
    parser.add_argument("--padding-bp", type=int, default=0)
    parser.add_argument("--gene-exclude-regions", choices=("none", "mhc"), default="mhc")
    parser.add_argument("--ld-wind-cm", type=float, default=1.0)
    parser.add_argument("--maf-min", type=float, default=None)
    parser.add_argument("--common-maf-min", type=float, default=0.05)
    parser.add_argument("--keep-indivs-file", default=None)
    parser.add_argument(
        "--regr-snps-file",
        default=None,
        help=(
            "Optional identity-only SNP list defining the persisted regression/output rows. "
            "Duplicate restriction keys collapse and non-identity columns are ignored."
        ),
    )
    parser.add_argument(
        "--regr-snps-exclude-regions",
        choices=kernel_regions.REGR_SNPS_EXCLUDE_REGIONS_CHOICES,
        default="mhc-and-centromeres",
        help=(
            "Curated regions subtracted after selecting bundled HapMap3 or "
            "--regr-snps-file SNPs."
        ),
    )
    parser.add_argument(
        "--exclude-regions",
        dest="regr_snps_exclude_regions",
        choices=kernel_regions.REGR_SNPS_EXCLUDE_REGIONS_CHOICES,
        default=argparse.SUPPRESS,
        help=argparse.SUPPRESS,
    )
    parser.add_argument("--genetic-map-hg19-sources", default=None)
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
        and batching settings. The builder requires explicit hg19 and either
        base rsID or base CHR/POS identity in a 1 cM PLINK-backed workflow.
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
        Workflow logging threshold. The open log uses hidden sibling
        ``.<output-name>.build-state/build-gene-ldscore-index.log`` so it is
        never part of the replaceable index transaction. After success, the
        closed log moves into ``<output_dir>/diagnostics``. Default is
        ``"INFO"``.

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
        intersection is empty after duplicate-group cleanup, an hg38 genetic
        map is supplied, or staged validation fails.
    FileExistsError
        If a valid index exists and ``overwrite`` is false, or a nonempty
        unrecognized directory occupies the destination.

    Notes
    -----
    Construction uses float64 adjusted-r-squared accumulation, retains
    negative values, and batches internal atom columns. Mutable baseline and
    PLINK duplicate effective-key groups are dropped in full with warnings and
    diagnostic rows; repeated restriction keys collapse because restrictions
    are sets. PLINK supplies the published variant labels, coordinates, and
    alleles. Each chromosome worker
    atomically persists its completed shard in a private run stage before
    releasing the in-memory record; ``threads`` therefore controls concurrent
    chromosome-local peak memory. A ``Finished chromosome N`` log record means
    that internal shard is durable, not that a partial public index exists. The
    destination remains missing, empty, or loadable at its prior version until
    the complete staged index passes reload validation and is committed. Stages
    are never resumed or reused. After destination reload validation succeeds,
    failure to remove builder-owned staging/backup data is warned with its
    retained path and does not change publication success.
    """
    if not isinstance(config, GeneLDScoreIndexBuildConfig):
        raise TypeError("build_gene_ldscore_index requires GeneLDScoreIndexBuildConfig.")
    return run_build_gene_ldscore_index_from_args(
        argparse.Namespace(
            baseline_annot_sources=",".join(config.baseline_annot_sources),
            plink_prefix=config.plink_prefix,
            output_dir=config.output_dir,
            gene_coordinate_file=config.gene_coordinate_file,
            genome_build=config.genome_build,
            snp_identifier=config.snp_identifier,
            padding_bp=config.padding_bp,
            gene_exclude_regions=config.gene_exclude_regions,
            ld_wind_cm=config.ld_wind_cm,
            maf_min=config.maf_min,
            common_maf_min=config.common_maf_min,
            keep_indivs_file=config.keep_indivs_file,
            regr_snps_file=config.regr_snps_file,
            regr_snps_exclude_regions=config.regr_snps_exclude_regions,
            genetic_map_hg19_sources=genetic_map_hg19_sources,
            genetic_map_hg38_sources=genetic_map_hg38_sources,
            _test_chromosomes=tuple(str(value) for value in range(1, 23)),
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
        select_index_eligible_gene_indices(
            catalog,
            gene_exclude_regions=gene_exclude_regions,
        )
    )
    frame = catalog.frame.copy().reset_index(drop=True)
    embedded = frame.loc[:, ["gene_id", "gene_name", "chrom", "start", "end", "genome_build", "catalog_line"]].copy()
    embedded.insert(0, "gene_index", np.arange(len(embedded), dtype=np.int64))
    embedded["source"] = catalog.source
    embedded["included"] = embedded["gene_index"].isin(eligible)
    embedded["exclusion_reason"] = np.where(embedded["included"], "", "excluded_gene_region")
    embedded["chromosome_gene_row"] = embedded.groupby("chrom", sort=False).cumcount().astype(np.int64)
    return embedded[
        [
            "gene_index",
            "gene_id",
            "gene_name",
            "chrom",
            "start",
            "end",
            "genome_build",
            "catalog_line",
            "source",
            "included",
            "exclusion_reason",
            "chromosome_gene_row",
        ]
    ]


@materializing_overwrite_guard(
    lambda args: (
        (getattr(args, "output_dir"), getattr(args, "overwrite", False), "RUN_FAILED.txt")
        if getattr(args, "output_dir", None)
        else None
    ),
    command="run_build_gene_ldscore_index_from_args(...)",
)
def run_build_gene_ldscore_index_from_args(args: argparse.Namespace) -> Path:
    """Build one exact v1 index through locked preflight and staged publication."""
    started = time.perf_counter()
    chromosomes = tuple(getattr(args, "_test_chromosomes", tuple(str(value) for value in range(1, 23))))
    allow_partial_for_tests = chromosomes != tuple(str(value) for value in range(1, 23))
    baseline_sources = tuple(split_cli_path_tokens(args.baseline_annot_sources))
    config = GeneLDScoreIndexBuildConfig(
        baseline_annot_sources=baseline_sources,
        plink_prefix=args.plink_prefix,
        output_dir=args.output_dir,
        gene_coordinate_file=args.gene_coordinate_file,
        genome_build=args.genome_build,
        snp_identifier=args.snp_identifier,
        padding_bp=args.padding_bp,
        gene_exclude_regions=args.gene_exclude_regions,
        ld_wind_cm=args.ld_wind_cm,
        maf_min=args.maf_min,
        common_maf_min=args.common_maf_min,
        keep_indivs_file=args.keep_indivs_file,
        regr_snps_file=getattr(args, "regr_snps_file", None),
        regr_snps_exclude_regions=getattr(args, "regr_snps_exclude_regions", "mhc-and-centromeres"),
        snp_batch_size=args.snp_batch_size,
        atom_batch_size=args.atom_batch_size,
        threads=args.threads,
    )

    index_path = Path(config.output_dir)
    with _gene_index_build_lock(index_path):
        live_log_path = _prepare_gene_index_log(index_path)
        with workflow_logging("build-gene-ldscore-index", live_log_path, log_level=args.log_level):
            _recover_gene_index_publication(
                index_path, _allow_partial_for_tests=allow_partial_for_tests
            )
            _preflight_gene_index_output(
                index_path,
                overwrite=bool(args.overwrite),
                _allow_partial_for_tests=allow_partial_for_tests,
            )
            log_inputs(
                output_dir=str(config.output_dir),
                genome_build=config.genome_build,
                snp_identifier=config.snp_identifier,
                chromosomes=", ".join(chromosomes),
                baseline_annot_sources=", ".join(config.baseline_annot_sources),
                gene_coordinate_file=Path(config.gene_coordinate_file).name,
                plink_prefix=config.plink_prefix,
                ld_wind_cm=config.ld_wind_cm,
                padding_bp=config.padding_bp,
                gene_exclude_regions=config.gene_exclude_regions,
                maf_min=config.maf_min,
                common_maf_min=config.common_maf_min,
                keep_indivs_file=config.keep_indivs_file or "all individuals",
                genetic_map="explicit hg19 map" if getattr(args, "genetic_map_hg19_sources", None) else "BIM cM fallback",
                regr_snps_file=config.regr_snps_file or "bundled HapMap3",
                regr_snps_exclude_regions=config.regr_snps_exclude_regions,
                snp_batch_size=config.snp_batch_size,
                atom_batch_size=config.atom_batch_size,
                threads=config.threads,
                effective_log_level=args.log_level,
                process_id=os.getpid(),
                hostname=socket.gethostname(),
            )
            try:
                catalog = GeneCatalog.load(config.gene_coordinate_file, require_canonical=True)
            except GeneCatalogValidationError as exc:
                issues_path = _write_gene_catalog_issues(index_path, exc.issues)
                for row in exc.issues.itertuples(index=False):
                    LOGGER.warning(
                        "Gene-coordinate catalog defect: line=%s field=%s reason=%s observed=%r",
                        row.catalog_line,
                        row.field,
                        row.reason,
                        row.observed_value,
                    )
                LOGGER.error(
                    "Gene-coordinate catalog validation failed with %s defect(s); full repair audit: %s",
                    len(exc.issues),
                    issues_path,
                )
                raise
            if catalog.genome_build != config.genome_build:
                raise LDSCInputError(
                    "Gene-coordinate catalog build does not match the index build: "
                    f"catalog={catalog.genome_build}, index={config.genome_build}. "
                    "Use a canonical catalog generated for the index build; no implicit liftover is performed."
                )
            stage_parent = _create_gene_index_transaction(index_path)
            try:
                published_path = _run_gene_ldscore_index_build(
                    args,
                    config,
                    chromosomes,
                    catalog=catalog,
                    stage_parent=stage_parent,
                    started=started,
                )
            except Exception:
                _remove_gene_index_transaction(stage_parent, published=False)
                raise
        _finalize_gene_index_log(live_log_path, published_path)
        return published_path


def _run_gene_ldscore_index_build(
    args: argparse.Namespace,
    config: GeneLDScoreIndexBuildConfig,
    chromosomes: tuple[str, ...],
    *,
    catalog: GeneCatalog,
    stage_parent: Path,
    started: float,
) -> Path:
    """Build chromosome shards into one transaction and publish it atomically."""
    staged_index = stage_parent / Path(config.output_dir).name
    embedded_catalog = _build_embedded_gene_catalog(
        catalog,
        genome_build=config.genome_build,
        gene_exclude_regions=config.gene_exclude_regions,
    )
    # Package-wide rsID identity deliberately carries no coordinate-build
    # semantics. The index contract independently records explicit hg19 for
    # provenance and gene/region projection in both builder modes.
    global_config = (
        GlobalConfig(snp_identifier="rsid")
        if config.snp_identifier == "rsid"
        else GlobalConfig(snp_identifier="chr_pos", genome_build=config.genome_build)
    )
    annotation_spec = AnnotationBuildConfig(baseline_annot_sources=config.baseline_annot_sources)
    regression_path = Path(config.regr_snps_file or packaged_hm3_curated_map_path())
    regression_keys = read_snp_restriction_keys(
        regression_path,
        config.snp_identifier,
        genome_build=config.genome_build,
    )
    regression_presets = kernel_regions.regr_snps_exclude_regions_choice_to_presets(config.regr_snps_exclude_regions)
    regression_regions = (
        None
        if not regression_presets
        else kernel_regions.load_preset_intervals(regression_presets, config.genome_build)
    )
    genetic_map = _load_builder_genetic_map(args)

    def build_one(chrom: str) -> StagedIndexChromosome:
        chrom_started = time.perf_counter()
        LOGGER.info(f"Starting chromosome {chrom}.")
        phase = "baseline annotation resolution"
        try:
            annotation_builder = AnnotationBuilder(global_config)
            public_bundle = annotation_builder.run(annotation_spec, chrom=chrom)
            baseline_builder_drops = coerce_identity_drop_frame(
                getattr(annotation_builder, "_identity_drop_frame", None)
            )
            if not baseline_builder_drops.empty:
                baseline_builder_drops = baseline_builder_drops.assign(
                    stage="gene_index_baseline_identity_cleanup"
                )
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
                snp_identifier=config.snp_identifier,
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
            if not baseline_builder_drops.empty:
                intersection.diagnostics["baseline_duplicate_rows_dropped"] += int(
                    (baseline_builder_drops["reason"] == "duplicate_identity").sum()
                )
            all_identity_drops = coerce_identity_drop_frame(
                pd.concat(
                    [baseline_builder_drops, intersection.dropped_rows],
                    ignore_index=True,
                )
            )
            kernel_bundle = kernel_ldscore.AnnotationBundle(
                metadata=intersection.metadata,
                annotations=intersection.annotations,
                baseline_columns=list(public_bundle.baseline_columns),
                query_columns=[],
            )
            chrom_catalog = embedded_catalog.loc[embedded_catalog["chrom"].astype(str) == chrom].sort_values(
                "chromosome_gene_row", kind="mergesort"
            )
            intervals = np.column_stack(
                [
                    chrom_catalog["start"].to_numpy(dtype=np.int64) - 1,
                    chrom_catalog["end"].to_numpy(dtype=np.int64),
                ]
            )
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
        catalog_genes = int(included.sum())
        genes_with_atoms = int(np.count_nonzero(np.diff(record.atom_model.gene_to_atom.indptr)))
        rows_before_genotype_qc = int(record.reference_rows_before_genotype_qc or record.total_reference_snps_all)
        genotype_qc_removed = int(record.genotype_qc_removed)
        maf_removed = int(record.maf_removed)
        evidence = {
            "pre_qc_rows": int(len(bim_rows)),
            **intersection.diagnostics,
            "baseline_content_sha256": _retained_baseline_content_sha256(
                intersection.metadata,
                intersection.annotations,
                config.snp_identifier,
            ),
            "reference_rows_before_genotype_qc": rows_before_genotype_qc,
            "annotation_intersection_removed": int(len(bim_rows) - rows_before_genotype_qc),
            "retained_reference_rows": int(record.total_reference_snps_all),
            "retained_common_reference_rows": int(record.total_reference_snps_common),
            "genotype_qc_removed": genotype_qc_removed,
            "maf_removed": maf_removed,
            "genotype_qc_or_maf_removed": genotype_qc_removed + maf_removed,
            "regression_rows": int(len(record.baseline_rows)),
            "baseline_plink_identity": f"inner_join_by_{config.snp_identifier}",
            "catalog_genes": catalog_genes,
            "genes_with_padded_atoms": genes_with_atoms,
            "atom_count": int(record.atom_model.n_atoms),
            "operator_nnz": int(record.operator.nnz),
            "nnz_Y": int(record.operator.nnz),
            "snp_batch_size": int(config.snp_batch_size),
            "atom_batch_size": int(config.atom_batch_size),
            "maf_filter_policy": "disabled" if config.maf_min is None else f">={config.maf_min}",
            "cm_source": record.cm_source,
        }
        try:
            component_metadata = _index_chromosome_metadata(
                chrom,
                record,
                snp_identifier=config.snp_identifier,
                genome_build=config.genome_build,
            )
            _stage_index_chromosome(staged_index, chrom, record)
            dropped_path = (
                staged_index / "diagnostics" / "dropped_snps" / f"chr{chrom}_dropped.tsv.gz"
            )
            dropped_path.parent.mkdir(parents=True, exist_ok=True)
            all_identity_drops.to_csv(
                dropped_path, sep="\t", index=False, compression="gzip", na_rep=""
            )
        except Exception as exc:
            LOGGER.error(
                f"Chromosome {chrom} failed during durable shard staging: "
                f"{type(exc).__name__}: {exc}"
            )
            raise
        evidence["elapsed_seconds"] = time.perf_counter() - chrom_started
        LOGGER.info(
            f"Finished chromosome {chrom}: catalog genes={catalog_genes}, "
            f"retained-reference={evidence['retained_reference_rows']}, "
            f"regression-rows={evidence['regression_rows']}, atoms={evidence['atom_count']}, "
            f"nnz(Y)={evidence['nnz_Y']}, elapsed={evidence['elapsed_seconds']:.3f}s."
        )
        if not len(record.baseline_rows):
            LOGGER.warning(
                f"Chromosome {chrom} has zero persisted regression SNP rows after restriction and region exclusion."
            )
        return StagedIndexChromosome(
            chromosome=chrom,
            evidence=evidence,
            component_metadata=component_metadata,
        )

    if config.threads == 1 or len(chromosomes) == 1:
        completed = {chrom: build_one(chrom) for chrom in chromosomes}
    else:
        max_workers = os.cpu_count() if config.threads == -1 else config.threads
        if max_workers is None or max_workers < 1:
            max_workers = max(1, (os.cpu_count() or 1) + 1 + config.threads)
        with ThreadPoolExecutor(max_workers=min(max_workers, len(chromosomes))) as pool:
            futures = {pool.submit(build_one, chrom): chrom for chrom in chromosomes}
            completed = {}
            try:
                for future in as_completed(futures):
                    result = future.result()
                    completed[result.chromosome] = result
            except Exception:
                for future in futures:
                    future.cancel()
                raise
    staged_chromosomes = {chrom: completed[chrom] for chrom in chromosomes}
    evidence_by_chrom = {
        chrom: staged_chromosomes[chrom].evidence for chrom in chromosomes
    }
    total_regression_rows = sum(
        values["regression_rows"] for values in evidence_by_chrom.values()
    )
    if total_regression_rows == 0:
        raise LDSCInputError(
            "Gene LD-score index has zero regression SNP rows across all selected chromosomes "
            "after applying the SNP restriction and --regr-snps-exclude-regions policy."
        )
    index_identity = _builder_index_identity(
        config,
        args,
        chromosomes,
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
            "regression_snps": "custom" if config.regr_snps_file else "bundled_hapmap3",
            "regr_snps_exclude_regions": config.regr_snps_exclude_regions,
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
            "catalog_genes": sum(values["catalog_genes"] for values in evidence_by_chrom.values()),
            "atom_count": sum(values["atom_count"] for values in evidence_by_chrom.values()),
            "operator_nnz": sum(values["operator_nnz"] for values in evidence_by_chrom.values()),
            "regression_rows": total_regression_rows,
            "retained_reference_rows": sum(values["retained_reference_rows"] for values in evidence_by_chrom.values()),
        },
    }
    try:
        _finalize_staged_index_artifact(
            staged_index,
            index_id=index_id,
            index_identity=index_identity,
            gene_catalog=embedded_catalog,
            chromosomes=staged_chromosomes,
            diagnostic_payload=diagnostic_payload,
        )
        published_path = _commit_staged_gene_ldscore_index(
            stage_parent,
            Path(config.output_dir),
            _allow_partial_for_tests=chromosomes
            != tuple(str(value) for value in range(1, 23)),
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
    """Hold nonblocking legacy and current locks for one absolute destination."""
    destination = index_path.expanduser().resolve()
    destination.parent.mkdir(parents=True, exist_ok=True)
    build_state = _gene_index_build_state_dir(destination)
    build_state.mkdir(parents=True, exist_ok=True)
    legacy_lock_path = destination.parent / f".{destination.name}.build-gene-ldscore-index.lock"
    lock_paths = (legacy_lock_path, build_state / "build-gene-ldscore-index.lock")
    handles = []
    try:
        for lock_path in lock_paths:
            handle = lock_path.open("a+", encoding="utf-8")
            try:
                fcntl.flock(handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
            except BlockingIOError as exc:
                handle.seek(0)
                owner = handle.read().strip() or "unknown owner"
                handle.close()
                raise LDSCInputError(
                    f"Another build-gene-ldscore-index process is active for {destination} ({owner}). "
                    f"Monitor {_gene_index_build_state_dir(destination) / 'build-gene-ldscore-index.log'}."
                ) from exc
            handle.seek(0)
            handle.truncate()
            handle.write(f"pid={os.getpid()} host={socket.gethostname()} target={destination}\n")
            handle.flush()
            handles.append(handle)
        yield
    finally:
        if handles:
            try:
                legacy_lock_path.unlink(missing_ok=True)
            except OSError as exc:
                warnings.warn(
                    f"The obsolete gene-index lock file could not be removed: {legacy_lock_path} "
                    f"({type(exc).__name__}: {exc}). It is safe to remove when no build is active.",
                    UserWarning,
                    stacklevel=2,
                )
        for handle in reversed(handles):
            try:
                fcntl.flock(handle.fileno(), fcntl.LOCK_UN)
            finally:
                handle.close()


def _preflight_gene_index_output(
    index_path: Path,
    *,
    overwrite: bool,
    _allow_partial_for_tests: bool = False,
) -> None:
    """Validate the publication destination without creating or mutating it."""
    if index_path.exists() and not index_path.is_dir():
        raise FileExistsError(f"Gene LD-score index output is not a directory: {index_path}")
    if not index_path.exists():
        return
    if not any(index_path.iterdir()) or _is_diagnostics_only_gene_index(index_path):
        return
    try:
        _load_gene_ldscore_index(
            index_path, _allow_partial_for_tests=_allow_partial_for_tests
        )
    except Exception as exc:
        raise FileExistsError(
            f"Gene LD-score index output directory is nonempty but invalid: {index_path}. "
            "Choose an empty directory; --overwrite does not replace unrecognized contents."
        ) from exc
    if not overwrite:
        raise FileExistsError(
            f"Gene LD-score index already exists: {index_path}. Pass --overwrite to rebuild it completely."
        )


def _recover_gene_index_publication(
    index_path: Path,
    *,
    _allow_partial_for_tests: bool = False,
) -> None:
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
            _load_gene_ldscore_index(
                destination, _allow_partial_for_tests=_allow_partial_for_tests
            )
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

    backup_paths: list[Path] = []
    valid_backups: list[Path] = []
    for candidate in candidates:
        backup = candidate / f"{destination.name}.backup"
        if not backup.is_dir():
            continue
        backup_paths.append(backup)
        marker_path = backup / ".gene-index-publication.json"
        try:
            marker = _read_json(marker_path, "backup publication marker")
            _load_gene_ldscore_index(
                backup, _allow_partial_for_tests=_allow_partial_for_tests
            )
        except Exception:
            continue
        if marker.get("target") == str(destination):
            valid_backups.append(backup)
    if not valid_backups and not backup_paths:
        for candidate in candidates:
            _remove_gene_index_transaction(candidate, published=False)
        return
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
    return destination.with_name(f".{destination.name}.build-state")


def _legacy_gene_index_build_state_dir(index_path: Path) -> Path:
    """Return the visible sidecar path used by the preceding log contract."""
    destination = index_path.expanduser()
    return destination.with_name(f"{destination.name}.build")


def _archive_gene_index_log(log_path: Path, history: Path, *, label: str = "") -> Path:
    """Move one closed operational log to a collision-resistant history path."""
    timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S.%fZ")
    suffix = f".{label}" if label else ""
    archived = history / f"build-gene-ldscore-index.{timestamp}{suffix}.log"
    os.replace(log_path, archived)
    return archived


def _migrate_legacy_gene_index_build_state(index_path: Path, build_state: Path) -> None:
    """Move recognized logs from the former visible ``<index>.build`` sidecar."""
    legacy = _legacy_gene_index_build_state_dir(index_path)
    if not legacy.is_dir():
        return
    history = build_state / "history"
    history.mkdir(parents=True, exist_ok=True)
    legacy_log = legacy / "build-gene-ldscore-index.log"
    if legacy_log.is_file():
        _archive_gene_index_log(legacy_log, history, label="legacy")
    legacy_history = legacy / "history"
    if legacy_history.is_dir():
        for source in sorted(legacy_history.iterdir()):
            if not source.is_file():
                continue
            destination = history / source.name
            if destination.exists():
                timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S.%fZ")
                destination = history / f"{source.stem}.migrated.{timestamp}{source.suffix}"
            os.replace(source, destination)
        try:
            legacy_history.rmdir()
        except OSError:
            pass
    try:
        legacy.rmdir()
    except OSError:
        warnings.warn(
            f"The former gene-index build-state directory contains unrecognized files and was "
            f"not removed: {legacy}. Recognized logs were migrated to {build_state}.",
            UserWarning,
            stacklevel=2,
        )


def _prepare_gene_index_log(index_path: Path) -> Path:
    """Archive prior attempts and return the hidden stable live-log path."""
    build_state = _gene_index_build_state_dir(index_path)
    build_state.mkdir(parents=True, exist_ok=True)
    _migrate_legacy_gene_index_build_state(index_path, build_state)
    log_path = build_state / "build-gene-ldscore-index.log"
    if log_path.exists():
        history = build_state / "history"
        history.mkdir(exist_ok=True)
        _archive_gene_index_log(log_path, history)
    issues_path = build_state / "gene_coordinate_catalog_issues.tsv.gz"
    if issues_path.exists():
        history = build_state / "history"
        history.mkdir(exist_ok=True)
        timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S.%fZ")
        os.replace(
            issues_path,
            history / f"gene_coordinate_catalog_issues.{timestamp}.tsv.gz",
        )
    return log_path


def _write_gene_catalog_issues(index_path: Path, issues: pd.DataFrame) -> Path:
    """Write the complete current-attempt catalog repair audit in build state."""
    path = _gene_index_build_state_dir(index_path) / "gene_coordinate_catalog_issues.tsv.gz"
    path.parent.mkdir(parents=True, exist_ok=True)
    issues.to_csv(path, sep="\t", index=False, na_rep="", compression="gzip")
    return path


def _finalize_gene_index_log(live_log_path: Path, published_path: Path) -> Path:
    """Move one closed successful log into the published index diagnostics.

    Log placement occurs after the workflow context writes its ``Finished``
    footer and closes the file handler. Failure is warning-only because the
    scientific index has already been published and reload-validated; the live
    log remains at its hidden path for inspection and later manual movement.
    """
    final_log_path = published_path / "diagnostics" / live_log_path.name
    try:
        final_log_path.parent.mkdir(parents=True, exist_ok=True)
        os.replace(live_log_path, final_log_path)
    except Exception as exc:
        message = (
            "The gene LD-score index was published and reload-validated, but the closed workflow "
            f"log could not be moved into {final_log_path} ({type(exc).__name__}: {exc}). The "
            f"published index remains valid; the complete log is retained at {live_log_path}."
        )
        try:
            with live_log_path.open("a", encoding="utf-8") as handle:
                handle.write(f"\n[WARNING] {message}\n")
        except OSError:
            pass
        warnings.warn(message, UserWarning, stacklevel=2)
        return live_log_path
    set_workflow_log_path(final_log_path)
    return final_log_path


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
        f"Regression-row policy: {config['regression_snps']} with "
        f"regr-snps-exclude-regions={config['regr_snps_exclude_regions']}."
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
            f"catalog genes={values['catalog_genes']}, "
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
        f"catalog genes={totals['catalog_genes']}, atoms={totals['atom_count']}, "
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


def _retained_baseline_content_sha256(
    authoritative_rows: pd.DataFrame,
    annotations: pd.DataFrame,
    snp_identifier: str,
) -> str:
    """Hash retained annotation values under only the active effective identity."""
    identity_columns = ["SNP"] if snp_identifier == "rsid" else ["CHR", "POS"]
    content = pd.concat(
        [
            authoritative_rows.loc[:, identity_columns].reset_index(drop=True),
            annotations.reset_index(drop=True),
        ],
        axis=1,
    )
    return _canonical_frame_sha256(content, sort_by=tuple(identity_columns))


def _catalog_identity_records(frame: pd.DataFrame) -> list[dict[str, object]]:
    """Return JSON-stable normalized catalog rows for the index identity."""
    columns = [
        "gene_id",
        "gene_name",
        "chrom",
        "start",
        "end",
        "genome_build",
        "catalog_line",
    ]
    normalized = frame.loc[:, columns].copy()
    normalized["gene_name"] = normalized["gene_name"].fillna("").astype(str)
    for column in ("gene_id", "chrom", "genome_build"):
        normalized[column] = normalized[column].astype(str)
    for column in ("start", "end", "catalog_line"):
        normalized[column] = pd.to_numeric(normalized[column], errors="raise").astype(int)
    return normalized.sort_values(
        ["gene_id", "catalog_line"], kind="stable"
    ).to_dict(orient="records")


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
        "chromosomes": list(chromosomes),
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
            "kind": "custom" if config.regr_snps_file else "bundled_hapmap3",
            "canonical_keys_sha256": regression_digest,
        },
        "regr_snps_exclude_regions": config.regr_snps_exclude_regions,
        "catalog": _catalog_identity_records(catalog.frame),
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
class StagedIndexChromosome:
    """Compact worker result returned after one chromosome shard is durable.

    Attributes
    ----------
    chromosome : str
        Canonical chromosome label for the staged shard.
    evidence : dict
        Scalar identity inputs and build diagnostics used by the coordinator.
        This mapping contains no chromosome payload tables or sparse matrices.
    component_metadata : dict
        Dimensions, ordered baseline columns, and sparse-format declarations
        that receive the complete ``index_id`` during coordinator finalization.
    """

    chromosome: str
    evidence: dict
    component_metadata: dict


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
    snp_identifier: str
    genome_build: str
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
    _allow_partial_for_tests: bool = False,
) -> Path:
    """Stage, reload, and atomically publish one complete in-memory index.

    Parameters
    ----------
    index_dir : path-like
        Public destination for one complete immutable index.
    index_identity : dict
        Canonical scientific identity payload used to calculate ``index_id``.
    gene_catalog : pandas.DataFrame
        Complete embedded catalog in canonical ``gene_index`` order.
    chromosomes : dict of str to IndexChromosomeData
        Ordered chromosome records to persist in the index.
    overwrite : bool
        Replace an existing valid complete index when true. Invalid nonempty
        destinations are never replaced.
    diagnostic_payload : dict, optional
        Successful-build diagnostic record. When supplied, payload-byte fields
        are added in place before the record is written.

    Returns
    -------
    pathlib.Path
        Published index directory after staged and destination reload
        validation.

    Raises
    ------
    FileExistsError
        If the destination is an existing valid index without ``overwrite``,
        or is a nonempty unrecognized directory.
    LDSCInputError
        If the complete staged or published index fails reload validation.

    Notes
    -----
    Existing valid indexes require ``overwrite=True`` and are kept loadable
    until their complete replacement has passed staged reload validation.
    Nonempty directories that are neither valid indexes nor recognized
    legacy diagnostics-only failed builds are rejected even with overwrite.
    Mutable live logs are stored in hidden sibling build state and never enter
    the replaceable transaction. The workflow wrapper moves the closed
    successful log into published diagnostics. Once the destination has been
    replaced and reload-validated, failure to remove the builder-owned
    transaction tree is reported as a warning and does not change publication
    success. This convenience seam writes supplied chromosome records into a
    new transaction; the CLI builder instead persists worker shards directly
    into its run transaction and moves that already-written tree at commit.
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
                _load_gene_ldscore_index(
                    destination, _allow_partial_for_tests=_allow_partial_for_tests
                )
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
    stage_parent = _create_gene_index_transaction(destination)
    staged_index = stage_parent / destination.name
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
        return _commit_staged_gene_ldscore_index(
            stage_parent,
            destination,
            _allow_partial_for_tests=_allow_partial_for_tests,
        )
    except Exception:
        _remove_gene_index_transaction(stage_parent, published=False)
        raise


def _create_gene_index_transaction(destination: Path) -> Path:
    """Create one marked run-specific sibling stage without touching the target."""
    destination = destination.expanduser()
    destination.parent.mkdir(parents=True, exist_ok=True)
    stage_parent = Path(
        tempfile.mkdtemp(prefix=f".{destination.name}.stage-", dir=destination.parent)
    )
    marker = {
        "artifact_type": "gene_ldscore_index_publication",
        "target": str(destination.resolve()),
    }
    _write_json(stage_parent / ".gene-index-publication.json", marker)
    (stage_parent / destination.name / "chromosomes").mkdir(parents=True)
    return stage_parent


def _commit_staged_gene_ldscore_index(
    stage_parent: Path,
    destination: Path,
    *,
    _allow_partial_for_tests: bool = False,
) -> Path:
    """Validate and install an already-written index without recopying payloads.

    The prior valid destination is held as a transaction-local backup until the
    replacement reloads successfully. Cleanup after that commit point is
    warning-only; pre-commit failure restores the backup and remains fatal.
    """
    staged_index = stage_parent / destination.name
    marker = _read_json(stage_parent / ".gene-index-publication.json", "publication marker")
    committed = False
    try:
        _write_json(staged_index / ".gene-index-publication.json", marker)
        _load_gene_ldscore_index(
            staged_index, _allow_partial_for_tests=_allow_partial_for_tests
        )
        backup = stage_parent / f"{destination.name}.backup"
        if destination.exists():
            _write_json(destination / ".gene-index-publication.json", marker)
            try:
                os.replace(destination, backup)
            except Exception:
                (destination / ".gene-index-publication.json").unlink(missing_ok=True)
                raise
        try:
            os.replace(staged_index, destination)
            _load_gene_ldscore_index(
                destination, _allow_partial_for_tests=_allow_partial_for_tests
            )
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
        raise
    _remove_gene_index_transaction(stage_parent, published=committed)
    return destination


def load_gene_ldscore_index(
    index_dir: str | Path,
) -> LoadedGeneLDScoreIndex:
    """Load one complete public autosomes 1--22 gene LD-score index."""
    return _load_gene_ldscore_index(index_dir, _allow_partial_for_tests=False)


def _load_gene_ldscore_index(
    index_dir: str | Path,
    *,
    _allow_partial_for_tests: bool,
) -> LoadedGeneLDScoreIndex:
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
    try:
        snp_identifier = normalize_snp_identifier_mode(str(index_identity.get("snp_identifier")))
    except Exception as exc:
        raise LDSCInputError("Gene LD-score index SNP identity metadata is invalid.") from exc
    if snp_identifier not in {"rsid", "chr_pos"}:
        raise LDSCInputError("Gene LD-score index supports only rsid or chr_pos identity metadata.")
    genome_build = index_identity.get("genome_build")
    if genome_build != "hg19":
        raise LDSCInputError("Gene LD-score index genome-build metadata must be explicit hg19.")
    if root.get("snp_identifier") != snp_identifier or root.get("genome_build") != genome_build:
        raise LDSCInputError("Gene LD-score index root identity metadata disagrees with index_identity.")
    chromosomes = tuple(str(chrom) for chrom in root.get("chromosomes", ()))
    expected_coverage = tuple(str(chrom) for chrom in index_identity.get("chromosomes", ()))
    if not chromosomes or chromosomes != expected_coverage:
        raise LDSCInputError("Gene LD-score index chromosome coverage is invalid.")
    full_coverage = tuple(str(value) for value in range(1, 23))
    if not _allow_partial_for_tests and chromosomes != full_coverage:
        raise LDSCInputError(
            "Gene LD-score index must cover autosomes 1 through 22; partial production indexes are unsupported."
        )
    catalog_path = index_path / "gene_catalog.parquet"
    if not catalog_path.exists():
        raise LDSCInputError("Gene LD-score index is missing gene_catalog.parquet.")
    catalog = pd.read_parquet(catalog_path)
    required_catalog = [
        "gene_index", "gene_id", "gene_name", "chrom", "start", "end", "genome_build",
        "catalog_line", "source", "included", "exclusion_reason", "chromosome_gene_row",
    ]
    if list(catalog.columns) != required_catalog:
        raise LDSCInputError("Gene LD-score index catalog schema is invalid.")
    expected_gene_indices = np.arange(len(catalog), dtype=np.int64)
    if not np.array_equal(catalog["gene_index"].to_numpy(), expected_gene_indices):
        raise LDSCInputError("Gene LD-score index catalog gene_index ordering is invalid.")
    catalog_authority = GeneCatalog.from_embedded_frame(catalog)
    if index_identity.get("catalog") != _catalog_identity_records(catalog):
        raise LDSCInputError(
            "Gene LD-score index embedded catalog disagrees with its semantic identity."
        )
    if catalog["included"].dtype != np.bool_:
        raise LDSCInputError("Gene LD-score index catalog included flags must use Boolean dtype.")
    source_values = catalog["source"].fillna("").astype(str)
    if (
        source_values.eq("").any()
        or source_values.nunique() != 1
        or source_values.str.contains(r"[/\\]", regex=True).any()
    ):
        raise LDSCInputError(
            "Gene LD-score index catalog source must be one nonempty source basename."
        )
    gene_policy = index_identity.get("gene_exclude_regions")
    if gene_policy not in {"none", "mhc"}:
        raise LDSCInputError("Gene LD-score index gene-exclusion policy is invalid.")
    eligible = select_index_eligible_gene_indices(
        catalog_authority,
        gene_exclude_regions=str(gene_policy),
    )
    expected_included = catalog.index.isin(eligible)
    expected_reasons = np.where(expected_included, "", "excluded_gene_region")
    if not np.array_equal(catalog["included"].to_numpy(dtype=bool), expected_included):
        raise LDSCInputError(
            "Gene LD-score index catalog inclusion policy disagrees with its coordinates and metadata."
        )
    if not np.array_equal(catalog["exclusion_reason"].fillna("").astype(str), expected_reasons):
        raise LDSCInputError(
            "Gene LD-score index catalog exclusion reasons disagree with its inclusion policy."
        )
    loaded: dict[str, IndexChromosomeData] = {}
    for chrom in chromosomes:
        chromosome_catalog = catalog.loc[catalog["chrom"].astype(str) == str(chrom)]
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
            snp_identifier=snp_identifier,
            genome_build=genome_build,
            expected_gene_rows=len(chromosome_catalog),
        )
    return LoadedGeneLDScoreIndex(
        index_id=index_id,
        chromosomes=chromosomes,
        snp_identifier=snp_identifier,
        genome_build=genome_build,
        index_identity=index_identity,
        gene_catalog=catalog,
        index_chromosomes=loaded,
    )


def run_indexed_ldscore(
    index_dir: str | Path,
    *,
    query_gene_list_sources: Sequence[str | Path],
    control_gene_list_file: str | Path | None = None,
    gene_list_resolution_policy: str = "strict",
    output_dir: str | Path,
    overwrite: bool = False,
    _allow_partial_for_tests: bool = False,
):
    """Assemble exact gene-list LD scores from one explicit index directory.

    Parameters
    ----------
    index_dir : path-like
        Exact validated single-index directory.
    query_gene_list_sources : sequence of path-like
        One-column gene-list files. Exact authoritative gene IDs and
        case-sensitive gene names are resolved against the embedded catalog in
        source order.
    control_gene_list_file : path-like or None, optional
        Optional fixed-control gene-list file. When omitted, no
        ``gene_control`` annotation is added.
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

    index = _load_gene_ldscore_index(
        index_dir, _allow_partial_for_tests=_allow_partial_for_tests
    )
    projection_build = str(index.index_identity.get("projection_build"))
    gene_policy = str(index.index_identity.get("gene_exclude_regions", "none"))
    catalog = GeneCatalog.from_embedded_frame(index.gene_catalog)
    if catalog.genome_build != projection_build:
        raise LDSCInputError("Embedded gene catalog build disagrees with the index projection build.")
    batch = resolve_gene_lists(
        query_gene_list_sources,
        catalog,
        control_path=control_gene_list_file,
        resolution_policy=gene_list_resolution_policy,
        gene_exclude_regions=gene_policy,
        index_chromosome_coverage=index.chromosomes,
    )
    from .ldscore_calculator import (
        _gene_list_gate_a_message,
        _log_gene_list_rejections,
        _log_gene_list_snp_support,
    )

    _log_gene_list_rejections(batch)
    if batch.has_fatal_gate_a_issues:
        LDScoreDirectoryWriter().write_gene_list_preflight(
            batch,
            LDScoreOutputConfig(output_dir=output_dir, overwrite=overwrite),
        )
        raise LDSCInputError(_gene_list_gate_a_message(batch))
    support = _indexed_gene_support(index)
    batch = batch.with_snp_support(support)
    _log_gene_list_snp_support(batch)
    statuses = _indexed_query_statuses(batch, support)
    usable = [
        selection
        for selection, status in zip(
            (item for item in batch.selections if item.input_role == "focal"),
            statuses,
            strict=True,
        )
        if status.status in {"ok", "warning"}
    ]
    control_resolution = next((item for item in batch.selections if item.input_role == "control"), None)
    if control_resolution is not None and not support.reindex(control_resolution.catalog_indices).fillna(0).gt(0).any():
        diagnostic_result = SimpleNamespace(query_statuses=statuses, gene_list_batch=batch)
        LDScoreDirectoryWriter().write_query_diagnostics(
            diagnostic_result,
            LDScoreOutputConfig(output_dir=output_dir, overwrite=overwrite),
        )
        raise LDSCInputError("The requested indexed control gene list has zero retained-SNP support.")
    if not usable:
        _write_indexed_all_skipped_diagnostics(
            statuses,
            batch,
            output_dir=output_dir,
            overwrite=overwrite,
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
        ld_regression_snps=frozenset(
            effective_merge_key_series(
                baseline_table,
                index.snp_identifier,
                context="indexed gene LD-score regression rows",
            ).astype(str)
        ),
        chromosome_results=[],
        count_config={"common_reference_snp_maf_min": 0.05, "common_reference_snp_maf_operator": ">="},
        config_snapshot=(
            GlobalConfig(snp_identifier="rsid")
            if index.snp_identifier == "rsid"
            else GlobalConfig(snp_identifier="chr_pos", genome_build=index.genome_build)
        ),
        overlap=overlap,
        query_statuses=statuses,
        gene_list_batch=batch,
        index_provenance={
            "index_id": index.index_id,
            "index_snp_identifier": index.snp_identifier,
            "index_genome_build": index.genome_build,
        },
    )
    result = LDScoreCalculator()._finalize_query_statuses(result, statuses)
    from .ldscore_calculator import _log_query_annotation_statuses

    _log_query_annotation_statuses(result.query_statuses)
    if control_resolution is not None and pd.to_numeric(
        result.baseline_table["gene_control"], errors="coerce"
    ).nunique(dropna=False) <= 1:
        LDScoreDirectoryWriter().write_query_diagnostics(
            result,
            LDScoreOutputConfig(output_dir=output_dir, overwrite=overwrite),
        )
        raise LDSCInputError("The requested indexed control gene list produced zero-variance LD scores.")
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


def _indexed_gene_support(index: LoadedGeneLDScoreIndex) -> pd.Series:
    """Return retained reference-SNP support for every embedded catalog row."""
    support = pd.Series(0, index=index.gene_catalog["gene_index"].astype(int), dtype="Int64")
    for chrom in index.chromosomes:
        record = index.index_chromosomes[chrom]
        catalog_rows = index.gene_catalog.loc[
            index.gene_catalog["chrom"].astype(str).eq(str(chrom))
        ].sort_values("chromosome_gene_row", kind="stable")
        counts = record.atom_model.gene_to_atom.astype(np.int64) @ record.atom_statistics.atom_count_all.astype(np.int64)
        support.loc[catalog_rows["gene_index"].astype(int).to_numpy()] = np.asarray(counts).reshape(-1)
    return support


def _indexed_query_statuses(batch, support: pd.Series) -> tuple[QueryAnnotationStatus, ...]:
    """Derive Gate A/B focal statuses from one shared indexed batch."""
    statuses: list[QueryAnnotationStatus] = []
    for selection in (item for item in batch.selections if item.input_role == "focal"):
        summary = batch.summary.loc[
            (batch.summary["input_role"] == "focal")
            & (batch.summary["source_ordinal"] == selection.source_ordinal)
        ].iloc[0]
        nonblank = int(summary["nonblank_input_rows"])
        rejected = int(summary["rejected_rows"])
        support_counts = support.reindex(selection.catalog_indices).fillna(0).astype(int)
        supported = int(support_counts.gt(0).sum())
        if not selection.canonical_gene_ids:
            status, reason = "skipped", "empty_gene_list" if nonblank == 0 else "zero_resolved_genes"
        elif supported == 0:
            status, reason = "skipped", "zero_annotation_snps"
        elif rejected:
            status, reason = "warning", "partial_gene_resolution"
        elif supported < len(selection.catalog_indices):
            status, reason = "warning", "partial_snp_support"
        else:
            status, reason = "ok", ""
        details_parts: list[str] = []
        if rejected:
            details_parts.append(
                f"{rejected} submitted row(s) were rejected during catalog resolution."
            )
        if selection.canonical_gene_ids and supported < len(selection.catalog_indices):
            details_parts.append(
                f"{len(selection.catalog_indices) - supported} resolved gene(s) have zero "
                "retained reference-SNP support."
            )
        if status != "ok":
            details_parts.append(
                "See diagnostics/gene_list_resolution_summary.tsv and "
                "diagnostics/gene_list_audit.tsv.gz."
            )
        statuses.append(
            QueryAnnotationStatus(
                selection.query,
                selection.source,
                "gene_list",
                status,
                reason,
                n_annotation_snps=(0.0 if reason == "zero_annotation_snps" else None),
                details=" ".join(details_parts) or None,
            )
        )
    return tuple(statuses)


def _write_indexed_all_skipped_diagnostics(
    statuses: tuple[QueryAnnotationStatus, ...],
    batch,
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
        gene_list_batch=batch,
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
    resolution: GeneSourceSelection,
) -> np.ndarray:
    """Translate global embedded-catalog indices into one chromosome atom selector."""
    selected = np.asarray(resolution.catalog_indices, dtype=np.int64)
    rows = embedded_catalog.loc[
        embedded_catalog["gene_index"].isin(selected)
        & (embedded_catalog["chrom"].astype(str) == str(chrom))
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

    persisted = kernel_ldscore.regression_mask_from_keys(
        metadata,
        regression_keys,
        args.snp_identifier,
        region_intervals=regression_regions,
    ).astype(bool)
    n_baseline_columns = baseline.shape[1]
    combined_annotation = np.column_stack([baseline, persisted])
    geno._currentSNP = 0
    combined_scores = np.asarray(
        geno.ldScoreVarBlocks(
            prepared.block_left,
            args.snp_batch_size,
            annot=combined_annotation,
        ),
        dtype=np.float64,
    )
    baseline_scores = combined_scores[:, :n_baseline_columns]
    regression_scores = combined_scores[:, n_baseline_columns]

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
    """Baseline annotations aligned to canonical PLINK-authored rows."""

    metadata: pd.DataFrame
    annotations: pd.DataFrame
    diagnostics: dict[str, int]
    dropped_rows: pd.DataFrame


def intersect_baseline_plink_by_identifier(
    baseline_rows: pd.DataFrame,
    baseline_annotations: pd.DataFrame,
    plink_rows: pd.DataFrame,
    *,
    snp_identifier: str,
    chrom: str,
) -> BaselinePlinkIntersection:
    """Drop ambiguous source groups and inner-join by one base identity mode."""
    if snp_identifier not in {"rsid", "chr_pos"}:
        raise LDSCInputError("Gene LD-score index matching supports only rsid or chr_pos identity.")
    if len(baseline_rows) != len(baseline_annotations):
        raise LDSCInputError(
            "Baseline metadata and annotation row counts differ before PLINK intersection."
        )
    baseline = _normalize_identity_rows(
        baseline_rows.assign(_source_row=np.arange(len(baseline_rows), dtype=np.int64)),
        label="baseline",
        chrom=chrom,
        snp_identifier=snp_identifier,
    )
    plink = _normalize_identity_rows(
        plink_rows,
        label="PLINK BIM",
        chrom=chrom,
        snp_identifier=snp_identifier,
    )
    baseline_cleanup = clean_identity_artifact_table(
        baseline,
        snp_identifier,
        context=f"gene-index baseline chromosome {chrom}",
        stage="gene_index_baseline_identity_cleanup",
        logger=LOGGER,
    )
    plink_cleanup = clean_identity_artifact_table(
        plink,
        snp_identifier,
        context=f"gene-index PLINK BIM chromosome {chrom}",
        stage="gene_index_plink_identity_cleanup",
        logger=LOGGER,
    )
    baseline = baseline_cleanup.cleaned
    plink = plink_cleanup.cleaned
    if baseline.empty or plink.empty:
        raise LDSCInputError(
            f"Gene LD-score index baseline/PLINK {snp_identifier} intersection is empty on chromosome {chrom} "
            "after duplicate identity cleanup."
        )
    baseline["_identity_key"] = effective_merge_key_series(
        baseline, snp_identifier, context=f"gene-index baseline chromosome {chrom}"
    )
    plink["_identity_key"] = effective_merge_key_series(
        plink, snp_identifier, context=f"gene-index PLINK BIM chromosome {chrom}"
    )
    baseline_keys = set(baseline["_identity_key"])
    plink_keys = set(plink["_identity_key"])
    matched = plink["_identity_key"].isin(baseline_keys)
    if not matched.any():
        raise LDSCInputError(
            f"Gene LD-score index baseline/PLINK {snp_identifier} intersection is empty on chromosome {chrom}."
        )
    matched_plink = (
        plink.loc[matched]
        .sort_values(["CHR", "POS", "SNP"], kind="mergesort")
        .reset_index(drop=True)
    )
    baseline_lookup = baseline.set_index("_identity_key", verify_integrity=True)
    matched_baseline = baseline_lookup.loc[matched_plink["_identity_key"]].reset_index()
    coordinate_discordant = 0
    if snp_identifier == "rsid":
        coordinate_discordant = int(
            np.count_nonzero(
                matched_plink["POS"].to_numpy(dtype=np.int64)
                != matched_baseline["POS"].to_numpy(dtype=np.int64)
            )
        )
    if coordinate_discordant:
        LOGGER.warning(
            f"Baseline/PLINK rsid intersection on chromosome {chrom} has "
            f"{coordinate_discordant} coordinate disagreement row(s); PLINK coordinates are authoritative."
        )
    snp_label_discordant = 0
    if snp_identifier == "chr_pos" and "SNP" in matched_baseline.columns:
        snp_label_discordant = int(
            np.count_nonzero(
                matched_plink["SNP"].astype(str).to_numpy()
                != matched_baseline["SNP"].astype(str).to_numpy()
            )
        )
        if snp_label_discordant:
            LOGGER.info(
                f"Baseline/PLINK chr_pos intersection on chromosome {chrom} has "
                f"{snp_label_discordant} SNP label disagreement row(s); "
                "PLINK SNP labels are authoritative and will be published."
            )
    annotation_rows = matched_baseline["_source_row"].to_numpy(dtype=np.int64)
    aligned_annotations = baseline_annotations.iloc[annotation_rows].reset_index(drop=True)
    diagnostics = {
        "baseline_rows": int(len(baseline_rows)),
        "plink_rows": int(len(plink_rows)),
        "baseline_duplicate_rows_dropped": int(len(baseline_cleanup.dropped)),
        "plink_bim_duplicate_rows_dropped": int(len(plink_cleanup.dropped)),
        "matched_rows": int(len(matched_plink)),
        "baseline_only_rows": int(len(baseline_keys - plink_keys)),
        "plink_only_rows": int(len(plink_keys - baseline_keys)),
        "coordinate_discordant_rows": coordinate_discordant,
        "snp_label_discordant_rows": snp_label_discordant,
    }
    drop_frames = [
        frame for frame in (baseline_cleanup.dropped, plink_cleanup.dropped) if not frame.empty
    ]
    dropped = pd.concat(drop_frames, ignore_index=True) if drop_frames else empty_identity_drop_frame()
    matched_plink = matched_plink.drop(columns=["_identity_key"])
    return BaselinePlinkIntersection(matched_plink, aligned_annotations, diagnostics, dropped)


def _normalize_identity_rows(
    frame: pd.DataFrame,
    *,
    label: str,
    chrom: str,
    snp_identifier: str,
) -> pd.DataFrame:
    """Normalize one strict-identity frame without dropping caller columns."""
    required = {"CHR", "POS"} | ({"SNP"} if snp_identifier == "rsid" or label == "PLINK BIM" else set())
    missing = required - set(frame.columns)
    if missing:
        raise LDSCInputError(
            f"Gene LD-score index build cannot compare {label} identities because columns {sorted(missing)} are missing."
        )
    normalized = frame.copy()
    normalized["CHR"] = normalize_chromosome_series(normalized["CHR"], context=label).astype(str)
    normalized["POS"] = positive_int_position_series(normalized["POS"], context=label).astype("int64")
    if "SNP" in normalized.columns:
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


def _index_row_digests(frame: pd.DataFrame, snp_identifier: str) -> tuple[str, str]:
    """Hash ordered effective identities and ordered published row metadata."""
    required = ["CHR", "SNP", "POS", "A1", "A2"]
    missing = [column for column in required if column not in frame.columns]
    if missing:
        raise LDSCInputError(
            "New gene LD-score indexes require PLINK-authored A1/A2 and canonical row "
            f"metadata; missing columns: {missing}."
        )
    effective_columns = ["SNP"] if snp_identifier == "rsid" else ["CHR", "POS"]
    return (
        _canonical_frame_sha256(frame.loc[:, effective_columns]),
        _canonical_frame_sha256(frame.loc[:, required]),
    )


def _index_chromosome_metadata(
    chrom: str,
    record: IndexChromosomeData,
    *,
    snp_identifier: str,
    genome_build: str,
) -> dict:
    """Return compact component metadata without the not-yet-known index ID."""
    baseline_columns = [
        column
        for column in record.baseline_rows.columns
        if column not in {"CHR", "SNP", "POS", "A1", "A2", "regression_ld_scores"}
    ]
    effective_digest, published_digest = _index_row_digests(
        _cast_index_baseline_rows(record.baseline_rows), snp_identifier
    )
    return {
        "artifact_type": "gene_ldscore_index",
        "chromosome": str(chrom),
        "snp_identifier": snp_identifier,
        "genome_build": genome_build,
        "effective_identity_sha256": effective_digest,
        "published_row_metadata_sha256": published_digest,
        "n_rows": len(record.baseline_rows),
        "n_baseline": len(baseline_columns),
        "baseline_columns": baseline_columns,
        "n_genes": record.atom_model.gene_to_atom.shape[0],
        "n_atoms": record.atom_model.n_atoms,
        "gene_to_atom_format": "csr_bool_int32",
        "ldscore_operator_format": "csr_float64_int32",
    }


def _write_index_chromosome_payload(
    chrom_path: Path,
    *,
    chrom: str,
    record: IndexChromosomeData,
) -> None:
    """Write one chromosome's large scientific payloads without shared metadata."""
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


def _stage_index_chromosome(
    staged_index: Path,
    chrom: str,
    record: IndexChromosomeData,
) -> Path:
    """Atomically install one payload-only shard in the private run stage.

    All payload files are first closed in a chromosome-private temporary
    directory. Renaming that directory to ``chrN`` is the durability boundary
    observed by the worker's subsequent ``Finished chromosome N`` log record.
    Shared and component identity metadata are intentionally deferred to the
    coordinator.
    """
    chromosomes_path = staged_index / "chromosomes"
    chromosomes_path.mkdir(parents=True, exist_ok=True)
    temporary = Path(
        tempfile.mkdtemp(prefix=f".chr{chrom}.tmp-", dir=chromosomes_path)
    )
    destination = chromosomes_path / f"chr{chrom}"
    try:
        _write_index_chromosome_payload(temporary, chrom=chrom, record=record)
        os.replace(temporary, destination)
    except Exception:
        _remove_gene_index_transaction(temporary, published=False)
        raise
    return destination


def _finalize_staged_index_artifact(
    index_path: Path,
    *,
    index_id: str,
    index_identity: dict,
    gene_catalog: pd.DataFrame,
    chromosomes: dict[str, StagedIndexChromosome],
    diagnostic_payload: dict | None = None,
) -> None:
    """Finalize coordinator-owned metadata for a complete staged index.

    ``chromosomes`` is already in canonical requested order. This function
    writes root identity, component identity, the gene catalog, and successful
    diagnostics only after every payload-only shard is durable.
    """
    _write_json(
        index_path / "metadata.json",
        {
            "artifact_type": "gene_ldscore_index",
            "index_id": index_id,
            "snp_identifier": index_identity["snp_identifier"],
            "genome_build": index_identity["genome_build"],
            "index_identity": index_identity,
            "chromosomes": list(chromosomes),
        },
    )
    gene_catalog.to_parquet(index_path / "gene_catalog.parquet", index=False)
    for chrom, result in chromosomes.items():
        _write_json(
            index_path / "chromosomes" / f"chr{chrom}" / "metadata.json",
            {**result.component_metadata, "index_id": index_id},
        )
    if diagnostic_payload is not None:
        diagnostic_payload["payload_bytes"] = _scientific_payload_bytes(index_path)
        chromosome_diagnostics = diagnostic_payload.get("chromosomes", {})
        for chrom in chromosomes:
            values = chromosome_diagnostics.get(chrom)
            if isinstance(values, dict):
                values["payload_bytes"] = _directory_payload_bytes(
                    index_path / "chromosomes" / f"chr{chrom}"
                )
        _write_json(
            index_path / "diagnostics" / "build-gene-ldscore-index.json",
            diagnostic_payload,
        )


def _write_index_artifact(
    index_path: Path,
    *,
    index_id: str,
    index_identity: dict,
    gene_catalog: pd.DataFrame,
    chromosomes: dict[str, IndexChromosomeData],
) -> None:
    """Write one complete staged index in the canonical single-directory layout."""
    index_path.mkdir(parents=True, exist_ok=True)
    _write_json(
        index_path / "metadata.json",
        {
            "artifact_type": "gene_ldscore_index",
            "index_id": index_id,
            "snp_identifier": index_identity["snp_identifier"],
            "genome_build": index_identity["genome_build"],
            "index_identity": index_identity,
            "chromosomes": list(chromosomes),
        },
    )
    gene_catalog.to_parquet(index_path / "gene_catalog.parquet", index=False)
    for chrom, record in chromosomes.items():
        chrom_path = index_path / "chromosomes" / f"chr{chrom}"
        chrom_path.mkdir(parents=True)
        _write_json(
            chrom_path / "metadata.json",
            {
                **_index_chromosome_metadata(
                    chrom,
                    record,
                    snp_identifier=str(index_identity["snp_identifier"]),
                    genome_build=str(index_identity["genome_build"]),
                ),
                "index_id": index_id,
            },
        )
        _write_index_chromosome_payload(chrom_path, chrom=chrom, record=record)


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
    snp_identifier: str,
    genome_build: str,
    expected_gene_rows: int,
) -> IndexChromosomeData:
    component_path = index_path / "chromosomes" / f"chr{chrom}"
    component_meta = _read_json(component_path / "metadata.json", f"chromosome {chrom}")
    _validate_metadata_identity(component_meta, index_id=index_id)
    if (
        component_meta.get("snp_identifier") != snp_identifier
        or component_meta.get("genome_build") != genome_build
    ):
        raise LDSCInputError("Gene LD-score index component identity metadata disagrees with the root.")
    if component_meta.get("chromosome") != str(chrom):
        raise LDSCInputError("Gene LD-score index component chromosome identity is invalid.")
    baseline_rows = pd.read_parquet(component_path / "baseline_rows.parquet")
    if len(baseline_rows) != int(component_meta.get("n_rows", -1)):
        raise LDSCInputError("Gene LD-score index baseline row count disagrees with component metadata.")
    baseline_columns = list(component_meta.get("baseline_columns", ()))
    allele_columns = [column for column in ("A1", "A2") if column in baseline_rows.columns]
    if allele_columns != ["A1", "A2"]:
        raise LDSCInputError("Gene LD-score index baseline rows require both PLINK allele columns A1 and A2.")
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
    effective_columns = ["SNP"] if snp_identifier == "rsid" else ["CHR", "POS"]
    if baseline_rows.duplicated(effective_columns, keep=False).any():
        raise LDSCInputError("Gene LD-score index baseline rows contain duplicate effective identities.")
    canonical_rows = baseline_rows.sort_values(["POS", "SNP"], kind="mergesort").reset_index(drop=True)
    if not baseline_rows.reset_index(drop=True).equals(canonical_rows):
        raise LDSCInputError("Gene LD-score index baseline rows are not in canonical genomic order.")
    effective_digest, published_digest = _index_row_digests(baseline_rows, snp_identifier)
    if component_meta.get("effective_identity_sha256") != effective_digest:
        raise LDSCInputError("Gene LD-score index effective identity digest is invalid.")
    if component_meta.get("published_row_metadata_sha256") != published_digest:
        raise LDSCInputError("Gene LD-score index published row metadata digest is invalid.")
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
