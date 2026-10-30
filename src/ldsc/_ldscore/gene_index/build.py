"""Build, validate, and query exact disjoint-atom gene LD-score indexes.

The offline workflow computes a PLINK-backed LD-score operator for one fixed
baseline, reference panel, regression-SNP policy, and gene projection. The online
workflow resolves gene lists against the embedded catalog, assembles their
Boolean union from the stored operator, and writes an ordinary self-contained
LD-score directory. It never discovers an index or silently falls back to live
LD calculation. During an offline build, completed chromosome payloads become
durable inside one private run transaction and are released from memory; only a
complete reload-validated transaction is published as an index.
"""

from __future__ import annotations

import argparse
from contextlib import closing, contextmanager
from dataclasses import dataclass, replace
from datetime import datetime, timezone
import fcntl
import gzip
import hashlib
import json
import logging
import os
from pathlib import Path
import socket
import sqlite3
import shutil
import sys
import tempfile
import time
import warnings
from typing import Sequence

import numpy as np
import pandas as pd
from scipy import sparse

from ..._runtime.worker_runtime import initialize_worker, native_thread_state
from ..._runtime.cli_help import CLIHelpFormatter, CHROMOSOME_PATH_HELP
from ..._runtime.workflow_logging import args_output_scope, FAILURE_MARKER, LOG_LEVEL_HELP
from ..._runtime.parallelism import (
    WORKERS_HELP, chromosome_pool, parse_workers, resolve_worker_count, run_chromosome_tasks,
)
from ..._runtime.progress import report_phase
from .. import compute as kernel_ldscore
from ..diagnostics import emit_reference_snp_funnel_notice, log_reference_snp_funnel_warning
from ..._annotation.overlap import compute_overlap
from ..._refpanel.adapters import resolve_keep_individuals
from ..._refpanel.funnel import concat_funnels, funnel_from_records, funnel_records, with_grid_rows
from ..._refpanel.prepared import ChromosomeAnnotations
from .atoms import (
    build_disjoint_atoms,
    compute_atom_statistics,
    iter_snp_atom_blocks,
    map_snps_to_atoms,
    validate_ldscore_operator,
)
from ..._core.chromosomes import normalize_chromosome
from ..._core.regions import RegionIntervals
from ..._core.identifiers import read_snp_restriction_keys
from ..._core.plink_files import PlinkFAMFile
from ..._core.snp_identity import (
    RestrictionIdentityKeys,
    allele_set_series,
    base_key_series,
    clean_identity_artifact_table,
    effective_merge_key_series,
    empty_identity_drop_frame,
    normalize_snp_identifier_mode,
    identity_mode_family,
    is_allele_aware_mode,
)
from ..._core import regions as kernel_regions
from ..._runtime.workflow_logging import (
    log_inputs,
    log_outputs,
    materializing_overwrite_guard,
    set_workflow_log_path,
    workflow_logging,
)
from ..._annotation.storage import AnnotationWorkspace, FrameSpool
from ..._annotation.sources import prepare_annotation_sources
from ..._core.chromosomes import normalize_chromosome, normalize_chromosome_series
from ..._core.coordinates import positive_int_position_series
from ..._core.columns import normalize_genome_build
from ...config import DEFAULT_SNP_BATCH_SIZE, GeneLDScoreIndexBuildConfig, GlobalConfig, get_global_config
from ...errors import LDSCInputError, LDSCInternalError, LDSCUsageError
from ..._annotation.gene_lists import (
    GeneCatalog,
    GeneCatalogValidationError,
    select_index_eligible_gene_indices,
)
from ..._genome.hm3 import packaged_hm3_curated_map_path
from ..._runtime.paths import inspect_plink_inputs, split_cli_path_tokens
from .artifact import (
    IndexChromosomeData,
    StagedIndexChromosome,
    _canonical_frame_sha256,
    _catalog_identity_records,
    _finalize_staged_index_artifact,
    _index_chromosome_metadata,
    _read_json,
    _remove_gene_index_transaction,
    _stage_index_chromosome,
    _write_json,
    calculate_index_id,
    index_global_config,
)
from . import artifact as index_artifact
from .scoring import load_gene_ldscore_index


LOGGER = logging.getLogger("LDSC.gene_ldscore_index")


class _GeneIndexArgumentParser(argparse.ArgumentParser):
    """Argument parser with the builder's deliberate identity omissions."""

    def error(self, message: str) -> None:
        missing: list[str] = []
        if "required" in message and "--snp-identifier" in message:
            missing.append("--snp-identifier is required; choose rsid, chr_pos, rsid_allele_aware, or chr_pos_allele_aware.")
        if "required" in message and "--genome-build" in message:
            missing.append("--genome-build is required; choose hg19 or hg38.")
        super().error(" ".join(missing) if missing else message)


def build_parser() -> argparse.ArgumentParser:
    """Build the offline gene LD-score index parser."""
    parser = _GeneIndexArgumentParser(
        prog="ldsc build-gene-ldscore-index",
        description="Build an exact PLINK-backed disjoint-atom gene LD-score index.",
        allow_abbrev=False,
    )
    parser.prog = 'ldsc build-gene-ldscore-index'
    parser.formatter_class = CLIHelpFormatter
    parser.description = 'Build a reusable gene LD-score index from hg19 or hg38 PLINK genotypes and baseline annotations.'
    parser.usage = (
        '%(prog)s --snp-identifier MODE --genome-build BUILD\n'
        '    --baseline-annot-sources SOURCES --gene-coordinate-file FILE [--padding-bp BP]\n'
        '    --plink-prefix PREFIX (--ld-wind-cm CM | --ld-wind-kb KB | --ld-wind-snps N)\n'
        '    --output-dir DIR [options]'
    )
    identity = parser.add_argument_group('Step 1: SNP identity and genome build')
    baseline = parser.add_argument_group('Step 2: Baseline annotations')
    genes = parser.add_argument_group('Step 3: Gene coordinate catalog')
    panel = parser.add_argument_group('Step 4: PLINK reference panel')
    windows = parser.add_argument_group('Step 5: LD window (choose one)')
    output = parser.add_argument_group('Step 6: Output')
    snps = parser.add_argument_group('Optional: regression SNPs and reference SNPs', description='By default, the regression SNPs are the HM3 map SNPs outside MHC and centromeres.')
    computing = parser.add_argument_group('Computing resources')
    logs = parser.add_argument_group('Logging')

    identity.add_argument(
        '--snp-identifier', choices=('rsid', 'chr_pos', 'rsid_allele_aware', 'chr_pos_allele_aware'), required=True,
        help=(
            'Required SNP matching rule: rsid uses SNP identifiers; chr_pos uses chromosome and position. '
            'The allele-aware variants also compare normalized alleles, following ldscore rules. No default.'
        ),
    )
    identity.add_argument(
        '--genome-build', type=normalize_genome_build, choices=('hg19', 'hg38'), required=True,
        help=(
            'Required build assertion for coordinate-bearing inputs: hg19 or hg38. Input build declarations '
            'must agree after alias normalization. No inference or liftover.'
        ),
    )

    baseline.add_argument(
        '--baseline-annot-sources', required=True, metavar='SOURCES',
        help=(
            'Required baseline annotation tables covering chromosomes 1-22. Accepts exact paths; this command '
            'requires full autosomal coverage even when using a glob. '
            + CHROMOSOME_PATH_HELP
        ),
    )

    genes.add_argument(
        '--gene-coordinate-file', required=True, metavar='FILE',
        help=(
            'Required headered TSV/TSV.GZ gene coordinate catalog in the declared build, with one-based inclusive intervals. '
            'Defines the full gene set stored in the index. Supply an exact file path; patterns are not expanded.'
        ),
    )
    genes.add_argument(
        '--padding-bp', type=int, default=0, metavar='BP',
        help=(
            'Extend each gene boundary by this many base pairs at both ends before finding covered SNPs. '
            'Default: 0, use gene bodies. Must be nonnegative and is fixed in the index.'
        ),
    )
    genes.add_argument(
        '--gene-exclude-regions', choices=('none', 'mhc'), default='mhc',
        help=(
            'Exclude MHC-overlapping genes before padding, or none to exclude no genes. Default: mhc. The '
            'choice is fixed in the index.'
        ),
    )

    panel.add_argument(
        '--plink-prefix', required=True, metavar='PREFIX',
        help=(
            'Required PLINK .bed/.bim/.fam filename prefix or chromosome suite covering 1-22. Accepts suite '
            "stems, '*' and '@' patterns; see --baseline-annot-sources for pattern rules."
        ),
    )
    panel.add_argument(
        '--keep-indivs-file', default=None, metavar='FILE',
        help=(
            'File listing one individual ID (IID) per row to retain from PLINK genotypes. If omitted, retain all '
            'individuals. Supply an exact file path; patterns are not expanded.'
        ),
    )
    panel.add_argument(
        '--genetic-map-hg19-sources', default=None, metavar='SOURCES',
        help=(
            'hg19 genetic-map files for deriving SNP centiMorgan coordinates. Accepts comma-separated exact file '
            "paths; '*' and '@' are not expanded. If omitted, use informative PLINK CM values; one usable "
            'coordinate source is required.'
        ),
    )
    panel.add_argument('--genetic-map-hg38-sources', default=None, metavar='SOURCES',
                      help='hg38 genetic-map files for --ld-wind-cm with --genome-build hg38. '
                           'Same path rules as --genetic-map-hg19-sources.')

    window = windows.add_mutually_exclusive_group(required=True)
    window.add_argument(
        '--ld-wind-cm', type=float, default=None, metavar='CM',
        help=(
            'Maximum LD-window distance on either side of each SNP, in centiMorgans. No default; requires '
            'informative PLINK CM values or a matching-build genetic map. Specify exactly one LD window.'
        ),
    )
    window.add_argument('--ld-wind-kb', type=float, default=None, metavar='KB',
                        help='Maximum LD-window distance in kilobases. No default; specify exactly one LD window.')
    window.add_argument('--ld-wind-snps', type=int, default=None, metavar='N',
                        help='Maximum LD-window distance in retained SNPs after identity and genotype filtering. '
                             'No default; specify exactly one LD window. Whole-chromosome windows are rejected.')

    output.add_argument(
        '--output-dir', required=True, metavar='DIR',
        help=(
            'Required destination for the complete reusable gene LD-score index and diagnostics.'
        ),
    )
    output.add_argument(
        '--overwrite', action='store_true', default=False,
        help=(
            'Replace an existing complete index after the new build succeeds. Default: off; stop if the '
            'destination already exists.'
        ),
    )

    snps.add_argument(
        '--regr-snps-file', default=None, metavar='FILE',
        help=(
            'Headered SNP identity table selecting the regression SNPs: stored output rows and regression weights. If '
            'omitted, use bundled HapMap3 SNPs. --regr-snps-exclude-regions then applies. Supply an exact file '
            'path; patterns are not expanded.'
        ),
    )
    snps.add_argument(
        '--regr-snps-exclude-regions', choices=kernel_regions.REGR_SNPS_EXCLUDE_REGIONS_CHOICES, default='mhc-and-centromeres',
        help=(
            'Exclude MHC, centromeres, both, or neither from --regr-snps-file or HapMap3 output SNPs. '
            'Default: mhc-and-centromeres. Reference SNPs and annotation counts are unchanged.'
        ),
    )
    snps.add_argument(
        '--maf-min', type=float, default=None, metavar='VALUE',
        help=(
            'Keep reference SNPs with MAF at or above this value, in [0, 0.5]. If omitted, apply no '
            'additional MAF threshold.'
        ),
    )
    snps.add_argument(
        '--common-maf-min', type=float, default=0.05, metavar='VALUE',
        help=(
            'Minimum MAF for common-SNP annotation counts, including the boundary. Default: 0.05; does not '
            'filter reference SNPs.'
        ),
    )

    computing.add_argument(
        '--workers', type=parse_workers, default=1, metavar='N',
        help='Maximum simultaneous chromosome builds. ' + WORKERS_HELP,
    )
    computing.add_argument(
        '--snp-batch-size', type=int, default=DEFAULT_SNP_BATCH_SIZE, metavar='N',
        help=(
            f'Number of reference SNPs processed per PLINK genotype batch. Default: {DEFAULT_SNP_BATCH_SIZE}; '
            'larger batches use more memory. LD scores do not depend on it.'
        ),
    )
    computing.add_argument(
        '--atom-batch-size', type=int, default=64, metavar='N',
        help=(
            'Number of disjoint gene-overlap regions processed together during index construction. Default: '
            '64; a smaller positive value reduces memory use.'
        ),
    )

    logs.add_argument(
        '--log-level', choices=('DEBUG', 'INFO', 'WARNING', 'ERROR'), default='INFO',
        help=LOG_LEVEL_HELP,
    )
    return parser


def main(argv: Sequence[str] | None = None):
    """Run the offline gene LD-score index workflow."""
    args = build_parser().parse_args(argv)
    published = run_build_gene_ldscore_index_from_args(args)
    summary = Path(published) / "diagnostics" / "build-gene-ldscore-index.json"
    emit_reference_snp_funnel_notice(_evidence_funnel(json.loads(summary.read_text(encoding="utf-8"))["chromosomes"]),
                                     summary, command="build-gene-ldscore-index")
    return published


def _evidence_funnel(evidence_by_chrom: dict) -> pd.DataFrame | None:
    """Combine the reference-SNP funnels recorded in per-chromosome build evidence."""
    return concat_funnels(funnel_from_records(values["reference_snp_funnel"]) for values in evidence_by_chrom.values())


def build_gene_ldscore_index(
    config: GeneLDScoreIndexBuildConfig,
    *,
    genetic_map_hg19_sources: str | None = None,
    genetic_map_hg38_sources: str | None = None,
    overwrite: bool = False,
    log_level: str | None = None,
) -> Path:
    """Build and transactionally publish one complete gene LD-score index.

    Parameters
    ----------
    config : GeneLDScoreIndexBuildConfig
        Validated baseline, PLINK, chromosome, projection, window, filtering,
        and batching settings. Requires an explicit hg19/hg38 build, identity
        mode, and exactly one cM, kb, or SNP-count LD window.
    genetic_map_hg19_sources : str, optional
        Comma-separated hg19 genetic-map paths for an hg19 cM window.
        When no matching map is supplied, cM windows require informative
        coordinates from each PLINK BIM. Maps are invalid for kb/SNP windows.
    genetic_map_hg38_sources : str, optional
        Comma-separated hg38 genetic-map paths for an hg38 cM window.
        Maps must match the explicit build; kb/SNP windows need no map.
    overwrite : bool, optional
        Rebuild and replace a complete existing valid index. Gene LD-score indexes do
        not support incremental updates or component reuse. Default is
        ``False``.
    log_level : {"DEBUG", "INFO", "WARNING", "ERROR"}, optional
        Workflow logging threshold. The open log uses hidden sibling
        ``.<output-name>.build-state/build-gene-ldscore-index.log`` so it is
        never part of the replaceable index transaction. After success, the
        closed log moves into ``<output_dir>/diagnostics``. Default is
        ``None``: use the registered ``GlobalConfig.log_level`` (``"INFO"``
        unless changed).

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
        intersection is empty after identity cleanup, build or SNP declarations
        conflict, a whole-chromosome window is requested, or staged validation fails.
    FileExistsError
        If a valid index exists and ``overwrite`` is false, or a nonempty
        unrecognized directory occupies the destination.

    Notes
    -----
    All autosomes 1 through 22 are required. Build declarations in applicable
    inputs must agree with the explicit build after alias normalization;
    construction performs neither build inference nor liftover.
    Allele-aware modes use shared allele normalization and identity cleanup.
    Annotations may omit both allele columns and match by base identity;
    supplying only one allele column is invalid. Shared-key coordinate or
    usable-allele contradictions stop the build with a repair table.

    Construction uses float64 adjusted-r-squared accumulation, retains
    negative values, and batches internal atom columns. Mutable baseline and
    PLINK duplicate effective-key groups are dropped in full with warnings and
    diagnostic rows; repeated restriction keys collapse because restrictions
    are sets. PLINK supplies the published variant labels, coordinates, and
    alleles. Each chromosome worker
    atomically persists its completed shard in a private run stage before
    releasing the in-memory record; ``workers`` therefore controls concurrent
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
            ld_wind_kb=config.ld_wind_kb,
            ld_wind_snps=config.ld_wind_snps,
            maf_min=config.maf_min,
            common_maf_min=config.common_maf_min,
            keep_indivs_file=config.keep_indivs_file,
            regr_snps_file=config.regr_snps_file,
            regr_snps_exclude_regions=config.regr_snps_exclude_regions,
            genetic_map_hg19_sources=genetic_map_hg19_sources,
            genetic_map_hg38_sources=genetic_map_hg38_sources,
            snp_batch_size=config.snp_batch_size,
            atom_batch_size=config.atom_batch_size,
            workers=config.workers,
            overwrite=overwrite,
            log_level=get_global_config().log_level if log_level is None else log_level,
        )
    )


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
    args_output_scope,
    command="run_build_gene_ldscore_index_from_args(...)",
)
def run_build_gene_ldscore_index_from_args(args: argparse.Namespace) -> Path:
    """Build one gene LD-score index through locked preflight and staged publication."""
    started = time.perf_counter()
    chromosomes = index_artifact.INDEX_CHROMOSOMES
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
        ld_wind_kb=getattr(args, "ld_wind_kb", None),
        ld_wind_snps=getattr(args, "ld_wind_snps", None),
        maf_min=args.maf_min,
        common_maf_min=args.common_maf_min,
        keep_indivs_file=args.keep_indivs_file,
        regr_snps_file=getattr(args, "regr_snps_file", None),
        regr_snps_exclude_regions=getattr(args, "regr_snps_exclude_regions", "mhc-and-centromeres"),
        snp_batch_size=args.snp_batch_size,
        atom_batch_size=args.atom_batch_size,
        workers=args.workers,
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
                gene_coordinate_file=Path(config.gene_coordinate_file).name,
                plink_prefix=config.plink_prefix,
                ld_wind_cm=config.ld_wind_cm,
                ld_wind_kb=config.ld_wind_kb,
                ld_wind_snps=config.ld_wind_snps,
                padding_bp=config.padding_bp,
                gene_exclude_regions=config.gene_exclude_regions,
                maf_min=config.maf_min,
                common_maf_min=config.common_maf_min,
                keep_indivs_file=config.keep_indivs_file or "all individuals",
                genetic_map=("not used" if config.ld_wind_cm is None else
                             f"explicit {config.genome_build} map" if getattr(args, f"genetic_map_{config.genome_build}_sources", None)
                             else "BIM cM fallback"),
                regr_snps_file=config.regr_snps_file or "bundled HapMap3",
                regr_snps_exclude_regions=config.regr_snps_exclude_regions,
                snp_batch_size=config.snp_batch_size,
                atom_batch_size=config.atom_batch_size,
                workers=config.workers,
                effective_log_level=args.log_level,
                process_id=os.getpid(),
                hostname=socket.gethostname(),
            )
            from ..._annotation.preflight import inspect_annotation_inputs
            from ..._refpanel.preflight import GENETIC_MAP_HEADER_CHECKS
            declared_inputs = inspect_annotation_inputs(
                header_checks=GENETIC_MAP_HEADER_CHECKS,
                baseline=config.baseline_annot_sources, plink=config.plink_prefix, chromosomes=chromosomes,
                require_complete_suite=True, mode=config.snp_identifier,
                files=[(name, getattr(args, name, None)) for name in (
                    'genetic_map_hg19_sources', 'genetic_map_hg38_sources')],
                scalar_files=[(name, getattr(args, name, None)) for name in (
                    'gene_coordinate_file', 'keep_indivs_file', 'regr_snps_file')],
                issues_path=_gene_index_build_state_dir(index_path)/'input_issues.tsv',
                plink_issues_path=_gene_index_build_state_dir(index_path)/'plink_input_issues.tsv')
            try:
                catalog = GeneCatalog.load(config.gene_coordinate_file, require_canonical=True)
            except GeneCatalogValidationError as exc:
                issues_path = _write_gene_catalog_issues(index_path, exc.issues)
                LOGGER.error(
                    "Gene coordinate catalog validation failed with %s defect(s); full repair audit: %s",
                    len(exc.issues),
                    issues_path,
                )
                raise
            from ..._genome.build_policy import GenomeBuildCheck
            check = GenomeBuildCheck(config.genome_build, context="gene LD-score index construction", logger=LOGGER)
            check.declared("gene catalog", config.gene_coordinate_file, catalog.genome_build)
            for path in declared_inputs.baseline:
                check.table("baseline annotations", path)
            check.paths("regression SNPs", config.regr_snps_file, restriction_mode=config.snp_identifier)
            check.resolve()
            stage_parent = _create_gene_index_transaction(index_path)
            try:
                published_path = _run_gene_ldscore_index_build(
                    args,
                    config,
                    chromosomes,
                    catalog=catalog,
                    stage_parent=stage_parent,
                    started=started,
                    declared_inputs=declared_inputs,
                )
            except Exception:
                _remove_gene_index_transaction(stage_parent, published=False)
                raise
        _finalize_gene_index_log(live_log_path, published_path)
        return published_path


@report_phase(LOGGER, 'computation', 'gene-index chromosome')
def _build_one_index_chromosome(chrom, shard, baseline_columns, baseline_drop, chrom_catalog,
                                config, chromosome_prefixes, genetic_map, regression_keys,
                                regression_regions, staged_index):
    """Read one staged chromosome and return only its persisted component metadata."""
    chrom_started = time.perf_counter()
    LOGGER.info(f"Starting chromosome {chrom}.")
    phase = "baseline annotation resolution"
    try:
        metadata = shard.metadata()
        baseline = pd.DataFrame(shard.read(columns=baseline_columns), columns=baseline_columns)
        if metadata.empty:
            raise LDSCInputError(f"No baseline annotation rows were found for chromosome {chrom}.")
        kernel_args = argparse.Namespace(
            bfile=chromosome_prefixes[chrom],
            keep=config.keep_indivs_file,
            maf_min=config.maf_min,
            maf=None,
            snp_identifier=config.snp_identifier,
            genome_build=config.genome_build,
            ld_wind_snps=config.ld_wind_snps,
            ld_wind_kb=config.ld_wind_kb,
            ld_wind_cm=config.ld_wind_cm,
            yes_really=False,
            snp_batch_size=config.snp_batch_size,
            common_maf_min=config.common_maf_min,
            genetic_map=genetic_map,
        )
        prefix = chromosome_prefixes[chrom]
        bim_rows = _read_bim_identity(prefix)
        phase = "baseline/PLINK identifier intersection"
        intersection = intersect_baseline_plink_by_identifier(
            metadata,
            baseline,
            bim_rows,
            snp_identifier=config.snp_identifier,
            chrom=chrom,
        )
        for frame in baseline_drop.frames():
            intersection.diagnostics["baseline_duplicate_rows_dropped"] += int(frame.reason.eq("duplicate_identity").sum())
            intersection.diagnostics["baseline_identity_rows_dropped"] += len(frame)
        kernel_bundle = ChromosomeAnnotations(
            metadata=intersection.metadata,
            annotations=intersection.annotations,
            baseline_columns=list(baseline_columns),
            query_columns=[],
        )
        # Only the intersection is needed from here on.
        grid_rows, metadata, baseline = len(metadata), None, None
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
            chromosome_prefixes=chromosome_prefixes,
            regression_keys=regression_keys,
            regression_regions=regression_regions,
            gene_intervals=intervals,
            included=included,
            padding_bp=config.padding_bp,
            atom_batch_size=config.atom_batch_size,
        )
    except Exception as exc:
        if isinstance(exc, LDSCUsageError) and "whole-chromosome" in str(exc):
            raise LDSCInputError(
                f"Gene-index construction rejects the whole-chromosome LD window on chromosome {chrom}. "
                "Specify a smaller --ld-wind-cm, --ld-wind-kb, or --ld-wind-snps value; "
                "the builder has no whole-chromosome override.") from exc
        exc.add_note(f"Chromosome {chrom} failed during {phase}.")
        LOGGER.error(
            f"Chromosome {chrom} failed during {phase}: {type(exc).__name__}: {exc}"
        )
        raise
    catalog_genes = int(included.sum())
    genes_with_atoms = int(np.count_nonzero(np.diff(record.atom_model.gene_to_atom.indptr)))
    rows_before_genotype_qc = int(record.reference_rows_before_genotype_qc or record.total_reference_snps_all)
    genotype_qc_removed = int(record.genotype_qc_removed)
    maf_removed = int(record.maf_removed)
    # Reference preparation saw only the baseline/PLINK intersection; the grid is every baseline row.
    record = replace(record, reference_snp_funnel=with_grid_rows(record.reference_snp_funnel, grid_rows))
    evidence = {
        **intersection.diagnostics,
        "baseline_content_sha256": _retained_baseline_content_sha256(
            intersection.metadata,
            intersection.annotations,
            config.snp_identifier,
        ),
        "reference_snp_funnel": funnel_records(record.reference_snp_funnel),
        "reference_rows_before_genotype_qc": rows_before_genotype_qc,
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
        with gzip.open(dropped_path, "wt") as stream:
            header = True
            for frame in baseline_drop.frames():
                frame.to_csv(stream, sep="\t", index=False, header=header, na_rep="")
                header = False
            intersection.dropped_rows.to_csv(stream, sep="\t", index=False, header=header, na_rep="")
    except Exception as exc:
        exc.add_note(f"Chromosome {chrom} failed during durable shard staging.")
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
        f"nonzero SNP-atom LD scores={evidence['nnz_Y']}, elapsed={evidence['elapsed_seconds']:.3f}s."
    )
    if not len(record.baseline_rows):
        LOGGER.warning(
            f"Chromosome {chrom} has zero persisted regression SNP rows after restriction and region exclusion. "
            "Index construction continues, but this chromosome supplies no regression SNP rows for indexed scoring. "
            "Check the requested regression-SNP restriction and region exclusions; correct unintended settings "
            "before using the index."
        )
    evidence["worker_runtime"] = native_thread_state()
    return StagedIndexChromosome(
        chromosome=chrom,
        evidence=evidence,
        component_metadata=component_metadata,
    )


@report_phase(LOGGER, 'workflow', 'gene-index build')
def _run_gene_ldscore_index_build(
    args: argparse.Namespace,
    config: GeneLDScoreIndexBuildConfig,
    chromosomes: tuple[str, ...],
    *,
    catalog: GeneCatalog,
    stage_parent: Path,
    started: float,
    declared_inputs,
) -> Path:
    """Build chromosome shards into one transaction and publish it atomically."""
    staged_index = stage_parent / Path(config.output_dir).name
    embedded_catalog = _build_embedded_gene_catalog(
        catalog,
        genome_build=config.genome_build,
        gene_exclude_regions=config.gene_exclude_regions,
    )
    # The index independently records the explicit projection build.
    global_config = index_global_config(config.snp_identifier, config.genome_build)
    regression_path = Path(config.regr_snps_file or packaged_hm3_curated_map_path())
    from ..._runtime.preflight import InputGate
    from ..._runtime.progress import PhaseProgress
    state_dir = _gene_index_build_state_dir(Path(config.output_dir))
    gate = InputGate('Reference/control content validation', state_dir/'input_issues.tsv')
    with PhaseProgress(LOGGER, 'validation', 'PLINK, regression restriction, and genetic map content'):
        regression_keys = gate.check('regression SNPs', regression_path,
            lambda: read_snp_restriction_keys(regression_path, config.snp_identifier, genome_build=config.genome_build))
        genetic_map = gate.check('genetic map', getattr(args, f'genetic_map_{config.genome_build}_sources', None),
                                 lambda: _load_builder_genetic_map(args))
        plink = inspect_plink_inputs(config.plink_prefix, chromosomes=chromosomes, require_complete_suite=True,
                                     discovered_members=declared_inputs.files['plink_members'])
        try:
            plink.require_valid(required_chromosomes=chromosomes)
        except LDSCInputError:
            pd.DataFrame(plink.issues).to_csv(state_dir/'plink_input_issues.tsv', sep='\t', index=False)
            gate.issues.extend(plink.issues)
        gate.finish()
    regression_presets = kernel_regions.regr_snps_exclude_regions_choice_to_presets(config.regr_snps_exclude_regions)
    regression_regions = (
        None
        if not regression_presets
        else kernel_regions.load_preset_intervals(regression_presets, config.genome_build)
    )
    chromosome_prefixes = plink.chromosome_prefixes

    with AnnotationWorkspace(stage_parent) as workspace:
        paths, declared, issues = declared_inputs.baseline, declared_inputs.declarations, []
        sources = prepare_annotation_sources(
            workspace, paths, [], mode=global_config.snp_identifier,
            declared_chromosomes=declared, input_issues=issues, header_widths=declared_inputs.widths,
        )
        missing = set(chromosomes).difference(sources.shards)
        if missing:
            raise LDSCInputError(f"No baseline annotation rows were found for chromosomes {sorted(missing)}.")
        baseline_drops = {chrom: FrameSpool(workspace.path / 'baseline-drops' / chrom) for chrom in chromosomes}
        for frame in sources.drops.frames():
            for chrom, records in frame.groupby(frame.CHR.astype(str), sort=False):
                if chrom in baseline_drops:
                    baseline_drops[chrom].append(records.assign(stage="gene_index_baseline_identity_cleanup"))

        _validate_baseline_plink_inputs(sources, chromosome_prefixes, config, chromosomes, state_dir)

        def task_arguments(chrom):
            chrom_catalog = embedded_catalog.loc[embedded_catalog["chrom"].astype(str) == chrom].sort_values(
                "chromosome_gene_row", kind="mergesort")
            # A task interpolates only its own chromosome's map rows.
            chrom_map = (None if genetic_map is None
                         else genetic_map.loc[genetic_map["CHR"] == normalize_chromosome(chrom)])
            return (chrom, sources.shards[chrom], sources.baseline_columns, baseline_drops[chrom],
                    chrom_catalog, config, chromosome_prefixes, chrom_map, regression_keys,
                    regression_regions, staged_index)


        max_workers = resolve_worker_count(config.workers, len(chromosomes))
        with chromosome_pool(max_workers, operation="build-gene-ldscore-index chromosome construction",
                             initializer=initialize_worker, initargs=(args.log_level,)) as pool:
            completed = {}
            def submit(chrom):
                LOGGER.info("Starting chromosome %s gene-index construction.", chrom)
                return pool.submit(_build_one_index_chromosome, *task_arguments(chrom))

            def collect(chrom, result):
                completed[result.chromosome] = result
                evidence = result.evidence
                LOGGER.info(
                    "Finished chromosome %s: catalog genes=%s, retained-reference=%s, "
                    "regression-rows=%s, atoms=%s, nonzero SNP-atom LD scores=%s, elapsed=%.3fs.",
                    result.chromosome, evidence["catalog_genes"], evidence["retained_reference_rows"],
                    evidence["regression_rows"], evidence["atom_count"], evidence["nnz_Y"],
                    evidence["elapsed_seconds"],
                )
                runtime = result.evidence["worker_runtime"]
                LOGGER.log(logging.INFO if runtime["verified"] else logging.WARNING,
                           "Finished chromosome %s; native thread policy %s: %s",
                           result.chromosome, "verified" if runtime["verified"] else "unverified or unmet (continuing)",
                           runtime)

            def on_failure(error, failed_chrom, outstanding):
                LOGGER.exception("Gene-index construction failed; affected in-flight chromosomes: %s.",
                                 list(outstanding))

            run_chromosome_tasks(chromosomes, max_workers, submit, collect, on_failure)

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
        chromosome_prefixes,
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
    map_source = (f"explicit {config.genome_build} map" if isinstance(map_identity, dict)
                  else "not used" if map_identity == "not_used" else "BIM cM fallback")
    diagnostic_payload = {
        "index_id": index_id,
        "elapsed_seconds": elapsed,
        "peak_rss_bytes": peak_rss,
        "snp_batch_size": config.snp_batch_size,
        "atom_batch_size": config.atom_batch_size,
        "workers": config.workers,
        "effective_log_level": args.log_level,
        "configuration": {
            "genome_build": config.genome_build,
            "snp_identifier": config.snp_identifier,
            "chromosomes": list(chromosomes),
            "baseline_annot_sources": list(config.baseline_annot_sources),
            "plink_prefix": config.plink_prefix,
            "ld_window": _index_ld_window(config),
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
        published_path = _commit_staged_gene_ldscore_index(stage_parent, Path(config.output_dir))
    except Exception as exc:
        LOGGER.error(
            f"Gene LD-score index build failed during staged validation/publication: "
            f"{type(exc).__name__}: {exc}"
        )
        raise
    _log_gene_index_summary(diagnostic_payload)
    log_reference_snp_funnel_warning(
        _evidence_funnel(evidence_by_chrom), published_path / "diagnostics" / "build-gene-ldscore-index.json",
        command="build-gene-ldscore-index")
    log_outputs(
        index_id=index_id,
        index=str(published_path),
        summary_json=str(published_path / "diagnostics" / "build-gene-ldscore-index.json"),
        publication="staged, reloaded, validated, and atomically replaced",
    )
    return published_path


@contextmanager
def _gene_index_build_lock(index_path: Path):
    """Hold a nonblocking build lock for one absolute destination."""
    destination = index_path.expanduser().resolve()
    destination.parent.mkdir(parents=True, exist_ok=True)
    build_state = _gene_index_build_state_dir(destination)
    build_state.mkdir(parents=True, exist_ok=True)
    handle = (build_state / "build-gene-ldscore-index.lock").open("a+", encoding="utf-8")
    try:
        try:
            fcntl.flock(handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            handle.seek(0)
            owner = handle.read().strip() or "unknown owner"
            raise LDSCInputError(
                f"Another build-gene-ldscore-index process is active for {destination} ({owner}). "
                f"Monitor {build_state / 'build-gene-ldscore-index.log'}."
            ) from exc
        handle.seek(0)
        handle.truncate()
        handle.write(f"pid={os.getpid()} host={socket.gethostname()} target={destination}\n")
        handle.flush()
        try:
            yield
        finally:
            fcntl.flock(handle.fileno(), fcntl.LOCK_UN)
    finally:
        handle.close()


def _preflight_gene_index_output(
    index_path: Path,
    *,
    overwrite: bool,
) -> None:
    """Validate the publication destination without creating or mutating it."""
    if index_path.exists() and not index_path.is_dir():
        raise FileExistsError(f"Gene LD-score index output is not a directory: {index_path}")
    if not index_path.exists():
        return
    if not any(index_path.iterdir()) or _is_failure_marker_only_gene_index(index_path):
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


def _recover_gene_index_publication(
    index_path: Path,
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
            load_gene_ldscore_index(backup)
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
    backup = valid_backups[0]
    if destination.exists():
        if _is_failure_marker_only_gene_index(destination):
            failure_marker = destination / FAILURE_MARKER
            if failure_marker.is_file():
                # Recovery is not completion of the current build; retain its marker.
                os.replace(failure_marker, backup / failure_marker.name)
        else:
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


def _archive_gene_index_log(log_path: Path, history: Path) -> Path:
    """Move one closed operational log to a collision-resistant history path."""
    timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S.%fZ")
    archived = history / f"build-gene-ldscore-index.{timestamp}.log"
    os.replace(log_path, archived)
    return archived


def _prepare_gene_index_log(index_path: Path) -> Path:
    """Archive prior attempts and return the hidden stable live-log path."""
    build_state = _gene_index_build_state_dir(index_path)
    build_state.mkdir(parents=True, exist_ok=True)
    log_path = build_state / "build-gene-ldscore-index.log"
    if log_path.exists():
        history = build_state / "history"
        history.mkdir(exist_ok=True)
        _archive_gene_index_log(log_path, history)
    plink_issues = build_state / "plink_input_issues.tsv"
    if plink_issues.exists():
        history = build_state / "history"
        history.mkdir(exist_ok=True)
        timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S.%fZ")
        os.replace(plink_issues, history / f"plink_input_issues.{timestamp}.tsv")
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


def _is_failure_marker_only_gene_index(path: Path) -> bool:
    """Recognize a failed build's destination: only the regular root ``RUN_FAILED.txt`` file.

    Builds keep their live logs in hidden sibling build state, so a failed
    build leaves at most the marker in the destination.
    """
    if not path.is_dir():
        return False
    marker = path / FAILURE_MARKER
    return list(path.iterdir()) == [marker] and marker.is_file() and not marker.is_symlink()


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
        "SNP sets: the baseline/PLINK identifier intersection defines the reference SNPs, which enter "
        "LD scores, counts and overlap; the configured regression restriction defines the regression SNPs."
    )
    LOGGER.info(
        f"Regression-row policy: {config['regression_snps']} with "
        f"regr-snps-exclude-regions={config['regr_snps_exclude_regions']}."
    )
    for chrom in payload["chromosomes"]:
        values = payload["chromosomes"][chrom]
        funnel = funnel_from_records(values["reference_snp_funnel"])
        reference = funnel.loc[funnel["funnel"].eq("reference_panel")].set_index("step")
        LOGGER.info(
            f"Chromosome {chrom} metrics: panel-rows={reference.loc['panel_rows', 'remaining']}, "
            f"off-grid-panel-SNPs={reference.loc['annotation_snp_grid', 'removed']}, "
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
            f"atoms={values['atom_count']}, nonzero SNP-atom LD scores={values['nnz_Y']}, "
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
        f"nonzero SNP-atom LD scores={totals['operator_nnz']}, payload_bytes={payload['payload_bytes']}, "
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
    build = normalize_genome_build(args.genome_build)
    other = "hg38" if build == "hg19" else "hg19"
    sources = getattr(args, f"genetic_map_{build}_sources", None)
    if getattr(args, f"genetic_map_{other}_sources", None):
        raise LDSCInputError(
            f"Gene-index genetic map declares {other}, but the provided genome build is {build}. "
            f"Supply --genetic-map-{build}-sources with verified {build} coordinates. "
            "No coordinates or build labels were changed."
        )
    if not sources:
        return None
    if args.ld_wind_cm is None:
        raise LDSCInputError("Genetic-map sources require --ld-wind-cm. Omit the map for kb or SNP-count windows.")
    from ..._genome.build_policy import GenomeBuildCheck
    from ..._refpanel.build_kernel import load_genetic_map_group
    check = GenomeBuildCheck(build, context="gene-index genetic map", logger=LOGGER)
    check.paths("genetic map", sources)
    check.resolve()
    return load_genetic_map_group(split_cli_path_tokens(sources))


def _retained_baseline_content_sha256(
    authoritative_rows: pd.DataFrame,
    annotations: pd.DataFrame,
    snp_identifier: str,
) -> str:
    """Hash retained annotation values under only the active effective identity."""
    identity_columns = ["SNP"] if identity_mode_family(snp_identifier) == "rsid" else ["CHR", "POS"]
    if is_allele_aware_mode(snp_identifier):
        authoritative_rows = authoritative_rows.assign(
            allele_set=allele_set_series(authoritative_rows, context="retained baseline identity")[0])
        identity_columns.append("allele_set")
    content = pd.concat(
        [
            authoritative_rows.loc[:, identity_columns].reset_index(drop=True),
            annotations.reset_index(drop=True),
        ],
        axis=1,
    )
    return _canonical_frame_sha256(content, sort_by=tuple(identity_columns))


def _index_ld_window(config) -> dict:
    """Describe the one validated LD window in canonical index metadata."""
    return next({"unit": unit, "value": getattr(config, f"ld_wind_{unit}")}
                for unit in ("cm", "kb", "snps") if getattr(config, f"ld_wind_{unit}") is not None)


def _genetic_map_identity(args, genetic_map: pd.DataFrame | None = None) -> str | dict:
    """Return the content-bound scientific identity of the effective cM map."""
    if args.ld_wind_cm is None:
        return "not_used"
    build = normalize_genome_build(args.genome_build)
    sources = split_cli_path_tokens(getattr(args, f"genetic_map_{build}_sources", None))
    if not sources:
        return "bim_cm"
    if genetic_map is None:
        genetic_map = _load_builder_genetic_map(args)
    return {
        "kind": f"explicit_{build}",
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
    chromosome_prefixes: dict[str, str],
) -> dict:
    """Build the path-insensitive canonical scientific identity."""
    plink_files = []
    first_prefix = None
    for chrom in chromosomes:
        prefix = chromosome_prefixes[chrom]
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
    fam = PlinkFAMFile(str(first_prefix) + ".fam")
    keep_indices = resolve_keep_individuals(config.keep_indivs_file, fam)
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
        "ld_window": _index_ld_window(config),
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
        load_gene_ldscore_index(staged_index)
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
        raise
    _remove_gene_index_transaction(stage_parent, published=committed)
    return destination


def build_plink_index_chromosome(
    chrom: str,
    baseline_bundle: ChromosomeAnnotations,
    args,
    *,
    regression_keys: set[str] | RestrictionIdentityKeys,
    regression_regions: RegionIntervals | None,
    gene_intervals: np.ndarray,
    included: np.ndarray,
    padding_bp: int,
    atom_batch_size: int,
    chromosome_prefixes: dict[str, str] | None = None,
) -> IndexChromosomeData:
    """Construct exact common and atom payloads from one PLINK chromosome.

    ``chromosome_prefixes`` is the optional workflow-validated mapping from
    normalized BIM chromosome labels to complete prefixes. Supplying it avoids
    resolving ``args.bfile`` again; direct callers may omit it to use the same
    content validation lazily through ``PlinkRefPanel``. Genotype alignment,
    windowing, regression SNPs, and atom projection are unchanged.

    Baseline and regression-weight scores share one traversal. Each atom batch
    repeats genotype traversal, bounding the active dense atom projection.
    Sparse float64 operator blocks remain live until horizontal assembly, so
    ``atom_batch_size`` does not bound the final operator or assembly copies.
    The chromosome baseline matrix is also explicitly materialized in float64.
    """
    from ...config import LDScoreConfig, RefPanelConfig
    from ..._refpanel.adapters import PlinkRefPanel

    global_options = {"snp_identifier": args.snp_identifier}
    if hasattr(args, "genome_build"):
        global_options["genome_build"] = args.genome_build
    panel = PlinkRefPanel(
        GlobalConfig(**global_options),
        RefPanelConfig(backend="plink", plink_prefix=args.bfile,
                       keep_indivs_file=getattr(args, "keep", None), maf_min=getattr(args, "maf_min", None)),
        chromosome_prefixes=chromosome_prefixes,
    )
    config = LDScoreConfig(
        ld_wind_snps=args.ld_wind_snps, ld_wind_kb=args.ld_wind_kb, ld_wind_cm=args.ld_wind_cm,
        snp_batch_size=args.snp_batch_size, common_maf_min=getattr(args, "common_maf_min", 0.05),
        whole_chromosome_ok=args.yes_really,
    )
    with panel.prepare_chromosome(chrom, baseline_bundle, config,
                                  genetic_map=getattr(args, "genetic_map", None)) as prepared:
        metadata = prepared.metadata
        baseline = np.asarray(prepared.annotations.read(), dtype=np.float64)
        geno = prepared.reader

        persisted = kernel_ldscore.regression_mask_from_keys(
            metadata,
            regression_keys,
            args.snp_identifier,
            region_intervals=regression_regions,
        ).astype(bool)
        n_baseline_columns = baseline.shape[1]
        output_rows = np.flatnonzero(persisted)
        geno._currentSNP = 0
        combined_scores = np.asarray(
            kernel_ldscore.ld_scores_from_genotypes(
                geno,
                prepared.block_left,
                args.snp_batch_size,
                annot=baseline,
                output_rows=output_rows,
                n_baseline=n_baseline_columns,
                weight_mask=persisted,
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
        overlap = compute_overlap(
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
                kernel_ldscore.ld_scores_from_genotypes(
                    geno,
                    prepared.block_left,
                    args.snp_batch_size,
                    annot=atom_block,
                    output_rows=output_rows,
                ),
                dtype=np.float64,
            )
            operator_blocks.append(sparse.csr_matrix(block_scores, dtype=np.float64))
        operator = (
            sparse.hstack(operator_blocks, format="csr", dtype=np.float64)
            if operator_blocks
            else sparse.csr_matrix((int(persisted.sum()), 0), dtype=np.float64)
        )
        operator_blocks = block_scores = None  # the assembled operator owns these values
        operator.sort_indices()
        validate_ldscore_operator(operator, n_rows=int(persisted.sum()), n_atoms=atom_model.n_atoms)

        identity_columns = ["CHR", "SNP", "POS", *[c for c in ("A1", "A2") if c in metadata.columns]]
        baseline_rows = metadata.loc[persisted, identity_columns].reset_index(drop=True)
        baseline_rows["regression_ld_scores"] = regression_scores
        for column_index, column in enumerate(baseline_bundle.baseline_columns):
            baseline_rows[column] = baseline_scores[:, column_index]
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
            reference_metadata=metadata,
            reference_rows_before_genotype_qc=prepared.reference_rows_before_genotype_qc,
            genotype_qc_removed=prepared.genotype_qc_removed,
            maf_removed=prepared.maf_removed,
            selected_individual_count=prepared.selected_individual_count,
            cm_source=prepared.cm_source,
            reference_snp_funnel=prepared.reference_snp_funnel,
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
    """Validate and align baseline annotations to consistent PLINK identities.

    Shared cleanup removes ambiguous identity groups in the selected base or
    allele-aware mode. If baseline metadata omits both allele columns,
    matching uses base keys after allele cleanup against the reference panel; one allele column
    is invalid. Shared rsIDs with different coordinates and shared base keys
    with incompatible usable normalized alleles raise ``LDSCInputError``
    carrying an ``input_issues`` repair table. Ordinary unmatched rows are
    omitted by the inner intersection.

    The returned metadata uses PLINK coordinates, labels, and alleles, in
    genomic order. Annotation values follow that same row order; diagnostics
    and dropped rows describe this intersection's identity cleanup. Neither
    input frame is modified.
    """
    snp_identifier = normalize_snp_identifier_mode(snp_identifier)
    family = identity_mode_family(snp_identifier)
    has_alleles = {"A1", "A2"}.intersection(baseline_rows.columns)
    if len(has_alleles) == 1:
        raise LDSCInputError("Baseline annotations must supply both A1 and A2, or omit both allele columns.")
    baseline_mode = snp_identifier if has_alleles else family
    if len(baseline_rows) != len(baseline_annotations):
        raise LDSCInputError(
            f"Baseline SNP metadata has {len(baseline_rows)} rows but its annotations have "
            f"{len(baseline_annotations)} rows before matching the PLINK panel. "
            "Supply baseline values and SNP metadata with one corresponding row per SNP in the same order."
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
        baseline_mode,
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
            "after duplicate identity cleanup. Check whether the baseline annotations and PLINK panel "
            "use the intended SNP identifiers and matching coordinates, and check both inputs for "
            "duplicate identities. Restore compatible inputs rather than dropping ambiguous SNPs to bypass this check."
        )
    issues = _baseline_plink_contradictions(baseline, plink, mode=snp_identifier, chrom=chrom)
    if not issues.empty:
        error = LDSCInputError(
            f"Baseline annotations and PLINK BIM contain {len(issues)} contradictory SNP records on chromosome {chrom}. "
            "Inspect the repair table; supply matching coordinates and alleles, or correct erroneous input records "
            "only after verification. No SNP labels, alleles, or coordinates were changed.")
        error.input_issues = issues
        raise error
    baseline["_identity_key"] = effective_merge_key_series(
        baseline, baseline_mode, context=f"gene-index baseline chromosome {chrom}"
    )
    plink["_identity_key"] = effective_merge_key_series(
        plink, baseline_mode, context=f"gene-index PLINK BIM chromosome {chrom}"
    )
    baseline_keys = set(baseline["_identity_key"])
    plink_keys = set(plink["_identity_key"])
    matched = plink["_identity_key"].isin(baseline_keys)
    if not matched.any():
        raise LDSCInputError(
            f"Gene LD-score index baseline/PLINK {snp_identifier} intersection is empty on chromosome {chrom}. "
            "Check that baseline annotations and the PLINK panel describe the same SNP population, "
            "use compatible identifiers, and have matching coordinates; supply compatible inputs."
        )
    matched_plink = (
        plink.loc[matched]
        .sort_values(["CHR", "POS", "SNP"], kind="mergesort")
        .reset_index(drop=True)
    )
    baseline_lookup = baseline.set_index("_identity_key", verify_integrity=True)
    matched_baseline = baseline_lookup.loc[matched_plink["_identity_key"]].reset_index()
    coordinate_discordant = 0
    snp_label_discordant = 0
    if family == "chr_pos" and "SNP" in matched_baseline.columns:
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
    aligned_annotations = baseline_annotations.take(annotation_rows)
    aligned_annotations.index = pd.RangeIndex(len(aligned_annotations))
    diagnostics = {
        "baseline_rows": int(len(baseline_rows)),
        "plink_rows": int(len(plink_rows)),
        "baseline_duplicate_rows_dropped": int(baseline_cleanup.dropped.reason.eq("duplicate_identity").sum()),
        "plink_bim_duplicate_rows_dropped": int(plink_cleanup.dropped.reason.eq("duplicate_identity").sum()),
        "baseline_identity_rows_dropped": int(len(baseline_cleanup.dropped)),
        "plink_bim_identity_rows_dropped": int(len(plink_cleanup.dropped)),
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


def _validate_baseline_plink_inputs(sources, prefixes, config, chromosomes, state_dir):
    """Collect metadata conflicts across chromosomes before starting LD workers."""
    issues = []
    causes = []
    for chrom in chromosomes:
        metadata = sources.shards[chrom].metadata()
        try:
            intersect_baseline_plink_by_identifier(
                metadata, pd.DataFrame(index=range(len(metadata))), _read_bim_identity(prefixes[chrom]),
                snp_identifier=config.snp_identifier, chrom=chrom)
        except LDSCInputError as exc:
            causes.append(exc)
            records = getattr(exc, "input_issues", pd.DataFrame([{"chrom": chrom, "reason": "invalid_intersection",
                                                                  "details": str(exc)}]))
            issues.append(records.assign(baseline_sources=",".join(config.baseline_annot_sources),
                                         plink_source=str(prefixes[chrom]) + ".bim"))
    if identity_mode_family(config.snp_identifier) == "rsid":
        cross_chromosome = _cross_chromosome_rsid_issues(sources, prefixes, config, chromosomes, state_dir)
        if not cross_chromosome.empty:
            issues.append(cross_chromosome.assign(baseline_sources=",".join(config.baseline_annot_sources),
                                                  reason="cross_chromosome_coordinate_mismatch"))
            causes.append(LDSCInputError("Shared baseline/PLINK rsIDs declare different chromosomes."))
    path = state_dir / "baseline_plink_issues.tsv"
    if issues:
        frame = pd.concat(issues, ignore_index=True)
        frame.to_csv(path, sep="\t", index=False)
        error = LDSCInputError(f"Baseline/PLINK consistency validation failed. Repair the contradictory or invalid "
                               f"SNP records listed in {path}, using matching coordinates and alleles. "
                               + " ".join(str(cause) for cause in causes[:5]) + " " +
                               "Correct metadata only after verification; no input records were changed. "
                               "See docs/user-guide/troubleshooting.md#build-gene-ldscore-index-input-conflicts.")
        error.input_issues = frame
        LOGGER.error("%s", error)
        raise error from ExceptionGroup("Baseline/PLINK input conflicts", causes)
    path.unlink(missing_ok=True)


def _cross_chromosome_rsid_issues(sources, prefixes, config, chromosomes, state_dir):
    """Join chromosome declarations on disk without retaining genome-wide metadata."""
    with tempfile.TemporaryDirectory(prefix="identity-", dir=state_dir) as temporary:
        with closing(sqlite3.connect(str(Path(temporary) / "coordinates.sqlite"))) as database:
            for table in ("baseline", "plink"):
                database.execute(f"CREATE TABLE {table} (SNP TEXT, CHR TEXT, POS INTEGER)")
            for chrom in chromosomes:
                baseline = sources.shards[chrom].metadata()
                plink = clean_identity_artifact_table(
                    _read_bim_identity(prefixes[chrom]), config.snp_identifier,
                    context=f"PLINK chromosome {chrom}", stage="gene_index_build_validation").cleaned
                for table, frame in (("baseline", baseline), ("plink", plink)):
                    frame.loc[:, ["SNP", "CHR", "POS"]].astype({"CHR": str}).to_sql(
                        table, database, if_exists="append", index=False, chunksize=1000)
            database.execute("CREATE INDEX baseline_snp ON baseline (SNP)")
            database.execute("CREATE INDEX plink_snp ON plink (SNP)")
            issues = pd.read_sql_query(
                "SELECT b.SNP AS base_key, b.CHR AS CHR_baseline, b.POS AS POS_baseline, "
                "p.CHR AS CHR_plink, p.POS AS POS_plink FROM baseline b JOIN plink p ON b.SNP=p.SNP "
                "WHERE b.CHR != p.CHR", database)
    if not issues.empty:
        issues["plink_source"] = issues.CHR_plink.map(lambda chrom: str(prefixes[chrom]) + ".bim")
    return issues


def _baseline_plink_contradictions(baseline, plink, *, mode, chrom):
    """Report conflicting coordinates or usable alleles for shared base keys."""
    left = baseline.copy()
    right = plink.copy()
    for frame in (left, right):
        frame["base_key"] = base_key_series(frame, mode, context="baseline/PLINK consistency")
        frame["allele_set"] = (allele_set_series(frame, context="baseline/PLINK consistency")[0]
                               if {"A1", "A2"}.issubset(frame.columns) else pd.NA)
    pairs = left.merge(right, on="base_key", suffixes=("_baseline", "_plink"), validate="one_to_one")
    masks = {"allele_mismatch": (pairs.allele_set_baseline.notna() & pairs.allele_set_plink.notna()
                                  & pairs.allele_set_baseline.fillna("").ne(pairs.allele_set_plink.fillna("")))}
    if identity_mode_family(mode) == "rsid":
        masks["coordinate_mismatch"] = pairs.CHR_baseline.ne(pairs.CHR_plink) | pairs.POS_baseline.ne(pairs.POS_plink)
    issues = [pairs.loc[mask].assign(reason=reason, chrom=chrom) for reason, mask in masks.items() if mask.any()]
    return pd.concat(issues, ignore_index=True) if issues else pd.DataFrame()


def _normalize_identity_rows(
    frame: pd.DataFrame,
    *,
    label: str,
    chrom: str,
    snp_identifier: str,
) -> pd.DataFrame:
    """Normalize one strict-identity frame without dropping caller columns."""
    required = {"CHR", "POS"} | ({"SNP"} if identity_mode_family(snp_identifier) == "rsid" or label == "PLINK BIM" else set())
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


