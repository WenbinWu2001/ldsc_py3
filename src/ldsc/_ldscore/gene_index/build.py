"""Build exact disjoint-atom gene LD-score indexes.

The ``build-gene-ldscore-index`` workflow computes a PLINK-backed LD-score
operator for one fixed baseline, reference panel, regression-SNP policy, and
gene projection. Completed chromosome payloads become durable in the staged
index under the hidden ``.build/`` directory and are released from memory; only
a complete reload-validated stage is published as an index. Writing and loading
the index format live in ``ldsc._ldscore.gene_index.artifact``, scoring gene
lists with an index in ``ldsc._ldscore.gene_index.scoring``.
"""

from __future__ import annotations

import argparse
from contextlib import ExitStack, closing, contextmanager
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
from ..._runtime.cli_help import CLIHelpFormatter, CHROMOSOME_PATH_HELP, add_output_dir_argument, add_log_level_argument
from ..._runtime.workflow_logging import args_output_scope, FAILURE_MARKER
from ..._runtime.parallelism import (
    WORKERS_HELP, chromosome_pool, parse_workers, resolve_worker_count, run_chromosome_tasks,
)
from ..._runtime.progress import report_phase, PhaseProgress
from .. import compute as kernel_ldscore
from ..diagnostics import chromosome_log, emit_reference_snp_funnel_notice, log_reference_snp_funnel_warning
from ..outputs import CHROMOSOME_LOG_DIR, snp_universe_policy_record
from ..source import REGRESSION_LD_SCORE_COLUMN
from ..._annotation.semantics import combine_annotation_types
from ..._refpanel.adapters import PlinkRefPanel, resolved_ref_panel
from ..._refpanel.funnel import ANNOTATION_SNP_GRID, REFERENCE_PANEL, concat_funnels, funnel_from_records, funnel_records
from ..._refpanel.plink_source import (
    inspect_plink_inputs, keep_individual_indices, read_keep_individuals, select_plink_chromosome,
)
from ..._refpanel.prepared import ChromosomeAnnotations
from ..._refpanel.windows import WholeChromosomeWindowError
from .atoms import (
    build_disjoint_atoms,
    compute_atom_statistics,
    iter_snp_atom_blocks,
    map_snps_to_atoms,
    validate_ldscore_operator,
)
from ..._core.regions import RegionIntervals
from ..._core.identifiers import read_snp_restriction_keys
from ..._core.plink_files import PlinkFAMFile
from ..._core.snp_identity import (
    RestrictionIdentityKeys,
    allele_set_series,
    base_key_series,
    emit_dropped_snps_notice,
    normalize_snp_identifier_mode,
    identity_mode_family,
    is_allele_aware_mode,
)
from ..._core import regions as kernel_regions
from ..._runtime.workflow_logging import (
    log_inputs,
    log_outputs,
    materializing_overwrite_guard,
    overwrite_failure_marker,
    set_workflow_log_path,
    workflow_logging,
)
from ..._runtime.scratch import FrameSpool, ScratchWorkspace
from ..._annotation.sources import prepare_annotation_sources
from ..._core.chromosomes import chrom_sort_key, normalize_chromosome, normalize_chromosome_series
from ..._core.coordinates import positive_int_position_series
from ..._core.columns import normalize_genome_build
from ...config import (
    DEFAULT_SNP_BATCH_SIZE, GeneLDScoreIndexBuildConfig, LDScoreConfig, RefPanelConfig, get_global_config,
)
from ...errors import LDSCInputError, LDSCInternalError
from ..._annotation.gene_lists import (
    GeneCatalog,
    GeneCatalogValidationError,
    select_index_eligible_gene_indices,
)
from ..._genome.hm3 import packaged_hm3_curated_map_path
from ..._runtime.paths import split_cli_path_tokens
from .artifact import (
    IndexChromosomeData,
    StagedIndexChromosome,
    _canonical_frame_sha256,
    _catalog_identity_records,
    _finalize_staged_index_artifact,
    _index_chromosome_metadata,
    _remove_gene_index_transaction,
    _stage_index_chromosome,
    calculate_index_id,
    index_global_config,
    load_gene_ldscore_index,
)
from . import artifact as index_artifact
from ..._annotation.preflight import inspect_annotation_inputs
from ..._genome.build_policy import GenomeBuildCheck
from ..._runtime.preflight import InputGate
from ..._refpanel.genetic_map import GENETIC_MAP_HEADER_CHECKS, GENETIC_MAP_TROUBLESHOOTING, load_genetic_map_group


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
    parser = _GeneIndexArgumentParser(allow_abbrev=False)
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
            'File listing one individual ID (IID) per row to retain from PLINK genotypes; an IID listed more '
            'than once counts once. If omitted, retain all individuals. Supply an exact file path; patterns are '
            'not expanded.'
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

    add_output_dir_argument(
        output,
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

    add_log_level_argument(logs)
    return parser


def main(argv: Sequence[str] | None = None):
    """Run the offline gene LD-score index workflow."""
    args = build_parser().parse_args(argv)
    published = run_build_gene_ldscore_index_from_args(args)
    summary = Path(published) / "diagnostics" / "build-gene-ldscore-index.json"
    evidence_by_chrom = json.loads(summary.read_text(encoding="utf-8"))["chromosomes"]
    emit_reference_snp_funnel_notice(_evidence_funnel(evidence_by_chrom), summary, command="build-gene-ldscore-index")
    emit_dropped_snps_notice(sorted((Path(published) / "diagnostics" / "dropped_snps").glob("*.tsv.gz")),
                             command="build-gene-ldscore-index")
    _emit_zero_regression_snp_notice(evidence_by_chrom)
    return published


def _emit_zero_regression_snp_notice(evidence_by_chrom: dict) -> None:
    """Print one stderr WARNING naming the chromosomes that persisted no regression SNPs.

    The build log already holds one WARNING per such chromosome; the notice
    exists because batch-scheduler users rarely read logs after a
    successful run.
    """
    empty = sorted((chrom for chrom, values in evidence_by_chrom.items() if not values["regression_rows"]), key=int)
    if not empty:
        return
    subject = f"chromosome {empty[0]} has" if len(empty) == 1 else f"chromosomes {', '.join(empty)} have"
    print(
        f"WARNING: build-gene-ldscore-index published an index in which {subject} zero regression SNPs after the "
        "regression-SNP restriction and region exclusions, so indexed scoring has no regression SNPs there. "
        "Check those settings and correct unintended ones before using the index.",
        file=sys.stderr,
    )


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
        Workflow logging threshold. The open log is
        ``<output_dir>/.build/build-gene-ldscore-index.log``, in the hidden
        build directory beside the index components; each chromosome worker
        writes ``<output_dir>/.build/chromosome-logs/chr<N>.log``. After
        success, the closed log moves into ``<output_dir>/diagnostics``, the
        chromosome logs into ``<output_dir>/diagnostics/chromosome-logs``, and
        ``.build/`` is removed. Default is ``None``: use the registered
        ``GlobalConfig.log_level`` (``"INFO"`` unless changed).

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
    that internal shard is durable, not that a partial public index exists. A
    missing destination is created when the build starts and holds only
    ``.build/`` until a build succeeds. A failure before publication begins
    keeps an earlier index complete and loadable; a failure during publication
    leaves no loadable index, because publication removes the old root
    ``metadata.json`` first and writes the new one last. Stages are never
    resumed or reused. After destination reload validation
    succeeds, failure to move the closed logs or remove ``.build/`` is warned
    with its retained path and does not change publication success.
    """
    if not isinstance(config, GeneLDScoreIndexBuildConfig):
        raise TypeError("build_gene_ldscore_index requires GeneLDScoreIndexBuildConfig.")
    with overwrite_failure_marker(
        config.output_dir, overwrite=overwrite, command="Python API: build_gene_ldscore_index(...)",
    ):
        return _build(
            config,
            genetic_map_sources={"hg19": genetic_map_hg19_sources, "hg38": genetic_map_hg38_sources},
            overwrite=overwrite,
            log_level=get_global_config().log_level if log_level is None else log_level,
        )


def _build_embedded_gene_catalog(
    catalog: GeneCatalog,
    *,
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
    """Build one gene LD-score index from parsed command-line arguments.

    The arguments become a ``GeneLDScoreIndexBuildConfig`` before anything is
    written. Every file the build touches is inside ``args.output_dir``: the
    lock, live log, per-chromosome logs, issue tables and staged index sit in
    its hidden ``.build/`` directory, which a successful build removes after
    moving the logs into ``diagnostics/`` and ``diagnostics/chromosome-logs/``.
    """
    config = GeneLDScoreIndexBuildConfig(
        baseline_annot_sources=args.baseline_annot_sources,
        plink_prefix=args.plink_prefix,
        output_dir=args.output_dir,
        gene_coordinate_file=args.gene_coordinate_file,
        genome_build=args.genome_build,
        snp_identifier=args.snp_identifier,
        padding_bp=args.padding_bp,
        gene_exclude_regions=args.gene_exclude_regions,
        ld_wind_cm=args.ld_wind_cm,
        ld_wind_kb=args.ld_wind_kb,
        ld_wind_snps=args.ld_wind_snps,
        maf_min=args.maf_min,
        common_maf_min=args.common_maf_min,
        keep_indivs_file=args.keep_indivs_file,
        regr_snps_file=args.regr_snps_file,
        regr_snps_exclude_regions=args.regr_snps_exclude_regions,
        snp_batch_size=args.snp_batch_size,
        atom_batch_size=args.atom_batch_size,
        workers=args.workers,
    )
    return _build(
        config,
        genetic_map_sources={"hg19": args.genetic_map_hg19_sources, "hg38": args.genetic_map_hg38_sources},
        overwrite=args.overwrite,
        log_level=args.log_level,
    )


def _build(
    config: GeneLDScoreIndexBuildConfig,
    *,
    genetic_map_sources: dict[str, str | None],
    overwrite: bool,
    log_level: str,
) -> Path:
    """Build one gene LD-score index through output preflight, a locked build and staged publication.

    ``genetic_map_sources`` maps ``"hg19"`` and ``"hg38"`` to the comma-separated
    genetic-map paths supplied for that build, or ``None``. The caller opens
    the failure-marker scope. The output preflight is read-only and runs
    before the lock creates the build directory, so a refused destination is
    left as it was; only a build directory that an earlier attempt left is
    cleaned for this one. The preflight, lock, build directory and
    publication all use the one destination path set here; the configuration
    has already expanded ``~`` and environment variables in it.
    """
    started = time.perf_counter()
    chromosomes = index_artifact.INDEX_CHROMOSOMES
    index_path = Path(config.output_dir)
    accepted_entries = _preflight_gene_index_output(index_path, overwrite=bool(overwrite))
    with _gene_index_build_lock(index_path):
        live_log_path = _prepare_gene_index_log(index_path, earlier_attempt=_BUILD_DIR in accepted_entries)
        with workflow_logging("build-gene-ldscore-index", live_log_path, log_level=log_level):
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
                             f"explicit {config.genome_build} map" if genetic_map_sources[config.genome_build]
                             else "BIM cM fallback"),
                regr_snps_file=config.regr_snps_file or "bundled HapMap3",
                regr_snps_exclude_regions=config.regr_snps_exclude_regions,
                snp_batch_size=config.snp_batch_size,
                atom_batch_size=config.atom_batch_size,
                workers=config.workers,
                effective_log_level=log_level,
                process_id=os.getpid(),
                hostname=socket.gethostname(),
            )
            declared_inputs = inspect_annotation_inputs(
                header_checks=GENETIC_MAP_HEADER_CHECKS,
                baseline=config.baseline_annot_sources, plink=config.plink_prefix, chromosomes=chromosomes,
                require_complete_suite=True, mode=config.snp_identifier,
                files=[(f'genetic_map_{build}_sources', genetic_map_sources[build]) for build in ('hg19', 'hg38')],
                scalar_files=[(name, getattr(config, name)) for name in (
                    'gene_coordinate_file', 'keep_indivs_file', 'regr_snps_file')],
                issues_path=_gene_index_build_dir(index_path)/'input_issues.tsv',
                plink_issues_path=_gene_index_build_dir(index_path)/'plink_input_issues.tsv')
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
            check = GenomeBuildCheck(config.genome_build, context="gene LD-score index construction", logger=LOGGER)
            check.declared("gene catalog", config.gene_coordinate_file, catalog.genome_build)
            for path in declared_inputs.baseline:
                check.table("baseline annotations", path)
            check.paths("regression SNPs", config.regr_snps_file, restriction_mode=config.snp_identifier)
            check.resolve()
            staged_index = _create_gene_index_stage(index_path)
            try:
                published_path = _run_gene_ldscore_index_build(
                    config,
                    chromosomes,
                    index_path=index_path,
                    accepted_entries=accepted_entries,
                    genetic_map_sources=genetic_map_sources,
                    overwrite=overwrite,
                    log_level=log_level,
                    catalog=catalog,
                    staged_index=staged_index,
                    started=started,
                    declared_inputs=declared_inputs,
                )
            except Exception:
                _remove_gene_index_transaction(staged_index)
                raise
        _finalize_gene_index_log(live_log_path, published_path)
        return published_path


def _index_chromosome_task(log_path, chrom, *task):
    """Construct one chromosome in a worker, its records going to its chromosome log.

    Spawn workers have no handler of their own; ``chromosome_log`` attaches the
    shared per-chromosome file handler that ``ldscore`` workers use, so worker
    records, phase timings and tracebacks reach ``log_path``.
    """
    with chromosome_log(log_path, chrom, "gene-index construction"):
        return _build_one_index_chromosome(chrom, *task)


@report_phase(LOGGER, 'computation', 'gene-index chromosome')
def _build_one_index_chromosome(chrom, shard, baseline_columns, baseline_drop, chrom_catalog,
                                config, ref_panel, genetic_map, regression_keys,
                                regression_regions, snp_label_discordant_rows, staged_index):
    """Read one staged chromosome and return only its persisted component metadata.

    The baseline shard goes to ``ref_panel.prepare_chromosome`` as the
    annotation SNP grid, as in direct ``ldscore``: the panel selects and
    cleans the PLINK rows and matches them to the grid. ``snp_label_discordant_rows``
    comes from the coordinator's baseline/PLINK consistency gate.
    """
    chrom_started = time.perf_counter()
    phase = "baseline annotation resolution"
    try:
        metadata = shard.metadata()
        baseline = pd.DataFrame(shard.read(columns=baseline_columns), columns=baseline_columns)
        if metadata.empty:
            raise LDSCInputError(f"No baseline annotation rows were found for chromosome {chrom}.")
        metadata, baseline = _in_index_row_order(metadata, baseline)
        # whole_chromosome_ok stays False: the index has no whole-chromosome override.
        ldscore_config = LDScoreConfig(
            ld_wind_snps=config.ld_wind_snps, ld_wind_kb=config.ld_wind_kb, ld_wind_cm=config.ld_wind_cm,
            snp_batch_size=config.snp_batch_size, common_maf_min=config.common_maf_min,
        )
        kernel_bundle = ChromosomeAnnotations(
            metadata=metadata,
            annotations=baseline,
            baseline_columns=list(baseline_columns),
            query_columns=[],
        )
        grid_rows = len(metadata)
        intervals = np.column_stack(
            [
                chrom_catalog["start"].to_numpy(dtype=np.int64) - 1,
                chrom_catalog["end"].to_numpy(dtype=np.int64),
            ]
        )
        included = chrom_catalog["included"].to_numpy(dtype=bool)
        phase = "genotype QC and atomic LD-score construction"
        record, build_evidence, plink_drops = build_plink_index_chromosome(
            chrom,
            kernel_bundle,
            ref_panel,
            ldscore_config,
            genetic_map=genetic_map,
            regression_keys=regression_keys,
            regression_regions=regression_regions,
            gene_intervals=intervals,
            included=included,
            padding_bp=config.padding_bp,
            atom_batch_size=config.atom_batch_size,
        )
    except Exception as exc:
        if isinstance(exc, WholeChromosomeWindowError):
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
    plink_drops = plink_drops.assign(stage="gene_index_plink_identity_cleanup")
    evidence = {
        **_baseline_plink_evidence(
            record.reference_snp_funnel, grid_rows=grid_rows, matched_rows=build_evidence["reference_rows_before_genotype_qc"],
            baseline_drops=list(baseline_drop.frames()), plink_drops=plink_drops,
            snp_label_discordant_rows=snp_label_discordant_rows,
        ),
        **build_evidence,
        "reference_snp_funnel": funnel_records(record.reference_snp_funnel),
        "retained_reference_rows": int(record.total_reference_snps_all),
        "retained_common_reference_rows": int(record.total_reference_snps_common),
        "genotype_qc_or_maf_removed": build_evidence["genotype_qc_removed"] + build_evidence["maf_removed"],
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
            plink_drops.to_csv(stream, sep="\t", index=False, header=header, na_rep="")
    except Exception as exc:
        exc.add_note(f"Chromosome {chrom} failed during durable shard staging.")
        LOGGER.error(
            f"Chromosome {chrom} failed during durable shard staging: "
            f"{type(exc).__name__}: {exc}"
        )
        raise
    evidence["elapsed_seconds"] = time.perf_counter() - chrom_started
    evidence["worker_runtime"] = native_thread_state()
    return StagedIndexChromosome(
        chromosome=chrom,
        evidence=evidence,
        component_metadata=component_metadata,
    )


@report_phase(LOGGER, 'workflow', 'gene-index build')
def _run_gene_ldscore_index_build(
    config: GeneLDScoreIndexBuildConfig,
    chromosomes: tuple[str, ...],
    *,
    index_path: Path,
    accepted_entries: frozenset[str],
    genetic_map_sources: dict[str, str | None],
    overwrite: bool,
    log_level: str,
    catalog: GeneCatalog,
    staged_index: Path,
    started: float,
    declared_inputs,
) -> Path:
    """Build chromosome shards into the staged index, then publish it with root metadata last.

    ``accepted_entries`` are the destination's top-level entries that the
    output preflight accepted; publication refuses a destination that has
    gained any other entry since.
    """
    embedded_catalog = _build_embedded_gene_catalog(
        catalog,
        gene_exclude_regions=config.gene_exclude_regions,
    )
    # The index independently records the explicit projection build.
    global_config = index_global_config(config.snp_identifier, config.genome_build)
    regression_path = Path(config.regr_snps_file or packaged_hm3_curated_map_path())
    build_dir = _gene_index_build_dir(index_path)
    gate = InputGate('Reference/control content validation', build_dir/'input_issues.tsv')
    with PhaseProgress(LOGGER, 'validation', 'PLINK, regression restriction, and genetic map content'):
        regression_keys = gate.check('regression SNPs', regression_path,
            lambda: read_snp_restriction_keys(regression_path, config.snp_identifier, genome_build=config.genome_build))
        genetic_map = gate.check('genetic map', genetic_map_sources[config.genome_build],
                                 lambda: _require_map_chromosomes(
                                     _load_builder_genetic_map(config, genetic_map_sources), chromosomes))
        plink = inspect_plink_inputs(config.plink_prefix, chromosomes=chromosomes, require_complete_suite=True,
                                     discovered_members=declared_inputs.inputs.plink_members)
        try:
            plink.require_valid(required_chromosomes=chromosomes)
        except LDSCInputError:
            pd.DataFrame(plink.issues).to_csv(build_dir/'plink_input_issues.tsv', sep='\t', index=False)
            gate.issues.extend(plink.issues)
        gate.finish()
    regression_presets = kernel_regions.regr_snps_exclude_regions_choice_to_presets(config.regr_snps_exclude_regions)
    regression_regions = (
        None
        if not regression_presets
        else kernel_regions.load_preset_intervals(regression_presets, config.genome_build)
    )
    chromosome_prefixes = plink.chromosome_prefixes

    with ScratchWorkspace(build_dir) as workspace:
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

        snp_label_discordant = _validate_baseline_plink_inputs(
            sources, chromosome_prefixes, config, chromosomes, build_dir)
        # One panel for the run: its keep list is read here, and every task's copy carries it.
        ref_panel = resolved_ref_panel(
            global_config,
            RefPanelConfig(backend="plink", plink_prefix=config.plink_prefix,
                           keep_indivs_file=config.keep_indivs_file, maf_min=config.maf_min),
            chromosome_prefixes=chromosome_prefixes,
        )
        ref_panel._resolve_run_inputs()

        def task_arguments(chrom):
            chrom_catalog = embedded_catalog.loc[embedded_catalog["chrom"].astype(str) == chrom].sort_values(
                "chromosome_gene_row", kind="mergesort")
            # A task interpolates only its own chromosome's map rows.
            chrom_map = (None if genetic_map is None
                         else genetic_map.loc[genetic_map["CHR"] == normalize_chromosome(chrom)])
            return (chrom, sources.shards[chrom], sources.baseline_columns, baseline_drops[chrom],
                    chrom_catalog, config, ref_panel, chrom_map, regression_keys,
                    regression_regions, snp_label_discordant[chrom], staged_index)


        chromosome_logs = _chromosome_log_dir(build_dir)
        chromosome_logs.mkdir()
        LOGGER.info("Chromosome logs: %s; a successful build moves them to %s.",
                    chromosome_logs, index_path / CHROMOSOME_LOG_DIR)
        max_workers = resolve_worker_count(config.workers, len(chromosomes))
        with chromosome_pool(max_workers, operation="build-gene-ldscore-index chromosome construction",
                             initializer=initialize_worker, initargs=(log_level,)) as pool:
            completed = {}
            def submit(chrom):
                LOGGER.info("Starting chromosome %s gene-index construction.", chrom)
                return pool.submit(_index_chromosome_task, chromosome_logs / f"chr{chrom}.log",
                                   *task_arguments(chrom))

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
                if not evidence["regression_rows"]:
                    LOGGER.warning(
                        f"Chromosome {result.chromosome} has zero persisted regression SNP rows after restriction "
                        "and region exclusion. Index construction continues, but this chromosome supplies no "
                        "regression SNP rows for indexed scoring. Check the requested regression-SNP restriction "
                        "and region exclusions; correct unintended settings before using the index."
                    )
                runtime = evidence["worker_runtime"]
                LOGGER.log(logging.INFO if runtime["verified"] else logging.WARNING,
                           "Chromosome %s native thread policy %s: %s",
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
        genetic_map_sources,
        chromosomes,
        regression_keys,
        genetic_map,
        catalog,
        evidence_by_chrom,
        chromosome_prefixes,
    )
    index_id = calculate_index_id(index_identity)
    # Indexed ldscore output records these exactly as a direct run on the same inputs would.
    annotation_types = combine_annotation_types(
        sources.baseline_columns, [values["annotation_types"] for values in evidence_by_chrom.values()])
    snp_universe_policy = snp_universe_policy_record(
        ref_panel_snps_file=None,
        regr_snps_file=config.regr_snps_file,
        regression_regions=regression_regions,
        retained_snp_count=sum(values["retained_reference_rows"] for values in evidence_by_chrom.values()),
        selected_snp_count_before_region_exclusion=sum(
            values["regression_selected_rows"] for values in evidence_by_chrom.values()),
        region_excluded_snp_count=sum(
            values["regression_selected_rows"] - values["regression_rows"] for values in evidence_by_chrom.values()),
        written_snp_count=total_regression_rows,
    )
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
        "effective_log_level": log_level,
        "configuration": {
            "genome_build": config.genome_build,
            "snp_identifier": config.snp_identifier,
            "chromosomes": list(chromosomes),
            "baseline_annot_sources": list(config.baseline_annot_sources),
            "plink_prefix": config.plink_prefix,
            "ld_window": config.ld_window.record(),
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
            "overwrite": bool(overwrite),
            "staging_reload_validation": "passed",
            "published_reload_validation": "passed",
            "index_path": str(index_path),
        },
        # Root metadata.json holds the combined annotation_types; repeating each chromosome's here is redundant.
        "chromosomes": {chrom: {key: value for key, value in values.items() if key != "annotation_types"}
                        for chrom, values in evidence_by_chrom.items()},
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
            annotation_types=annotation_types,
            snp_universe_policy=snp_universe_policy,
            diagnostic_payload=diagnostic_payload,
        )
        published_path = _publish_staged_gene_ldscore_index(
            staged_index, index_path, accepted_entries=accepted_entries)
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
        publication="staged, reload-validated, published with root metadata last, and reload-validated",
    )
    return published_path


@contextmanager
def _gene_index_build_lock(index_path: Path):
    """Hold a nonblocking build lock for one destination, inside its build directory.

    A successful build removes its build directory while still holding the
    lock, so a lock taken on a file that is no longer at the lock path is
    retried on a fresh lock file.
    """
    destination = index_path.absolute()
    lock_path = _gene_index_build_dir(index_path) / _LOCK_FILE
    for _attempt in range(3):
        lock_path.parent.mkdir(parents=True, exist_ok=True)
        handle = lock_path.open("a+", encoding="utf-8")
        try:
            fcntl.flock(handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            handle.seek(0)
            owner = handle.read().strip() or "unknown owner"
            handle.close()
            raise LDSCInputError(
                f"Another build-gene-ldscore-index process is active for {destination} ({owner}). "
                f"Monitor {lock_path.parent / _LOG_FILE}."
            ) from exc
        if lock_path.exists() and lock_path.stat().st_ino == os.fstat(handle.fileno()).st_ino:
            break
        handle.close()
    else:
        raise LDSCInputError(
            f"Could not hold the build lock for {destination}: another build kept replacing "
            f"{lock_path}. Retry once no build-gene-ldscore-index process targets this directory."
        )
    try:
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


# Root entries a published index owns; publication replaces them, root metadata last.
_INDEX_COMPONENTS = ("metadata.json", "gene_catalog.parquet", "chromosomes", "diagnostics")


def _preflight_gene_index_output(
    index_path: Path,
    *,
    overwrite: bool,
) -> frozenset[str]:
    """Validate the publication destination without changing anything in it.

    This runs before the build lock creates the build directory, so a build
    directory it sees was there before this run. It records an earlier
    attempt only if it holds the lock or log file that every attempt writes
    there first; an empty one, which a successful build can leave on a shared
    filesystem, is harmless; any other is refused. The build directory and a
    failure marker are the builder's own, and a directory holding nothing
    else is reused. Only a successful build removes the build directory, so
    index components beside one that records an earlier attempt, without a
    loadable index, are the remains of an interrupted build or publication
    and need ``overwrite``. Any other content is refused, including component
    names without such a build directory.

    Returns the names of the top-level entries it accepted, which publication
    checks again; ``.build`` is among them only when it is empty or records an
    earlier attempt.
    """
    if index_path.exists() and not index_path.is_dir():
        raise FileExistsError(f"Gene LD-score index output is not a directory: {index_path}")
    if not index_path.exists():
        return frozenset()
    _refuse_unusable_build_dir(index_path)
    entries = frozenset(entry.name for entry in index_path.iterdir())
    unrecognized = FileExistsError(
        f"Gene LD-score index output directory is nonempty but invalid: {index_path}. "
        "Choose an empty directory; --overwrite does not replace unrecognized contents."
    )
    build_dir = index_path / _BUILD_DIR
    earlier_attempt = _BUILD_DIR in entries and _holds_build_evidence(build_dir)
    if _BUILD_DIR in entries and not earlier_attempt and any(build_dir.iterdir()):
        raise unrecognized
    if _is_build_leftover_only(index_path):
        return entries
    try:
        load_gene_ldscore_index(index_path)
    except Exception as exc:
        if earlier_attempt and all(map(_is_owned_index_entry, index_path.iterdir())):
            if overwrite:
                return entries
            raise FileExistsError(
                f"Gene LD-score index output directory holds an incomplete gene LD-score index, "
                f"most likely from an interrupted publication: {index_path}. Pass --overwrite to rebuild it."
            ) from exc
        raise unrecognized from exc
    if not overwrite:
        raise FileExistsError(
            f"Gene LD-score index already exists: {index_path}. Pass --overwrite to rebuild it completely."
        )
    return entries


def _holds_build_evidence(build_dir: Path) -> bool:
    """Return whether a build directory holds the lock or log file that every build attempt writes there first."""
    return any((build_dir / name).is_file() for name in (_LOCK_FILE, _LOG_FILE))


def _refuse_unusable_build_dir(index_path: Path) -> None:
    """Refuse a build directory that is a symlink or not a directory, or that holds a symlink.

    The lock, log preparation and input checks create, truncate and remove
    entries of the build directory; following a link there would act on
    files outside the index directory.
    """
    build_dir = index_path / _BUILD_DIR
    if build_dir.is_symlink():
        problems = [f"{build_dir} is a symlink"]
    elif not build_dir.exists():
        return
    elif not build_dir.is_dir():
        problems = [f"{build_dir} is not a directory"]
    else:
        problems = [f"{entry} is a symlink" for entry in sorted(build_dir.iterdir()) if entry.is_symlink()]
    if problems:
        raise FileExistsError(
            f"Gene LD-score index build directory cannot be used: {'; '.join(problems)}. "
            "The builder keeps its lock, logs and staged index in the hidden .build/ directory of the "
            "output directory and never follows a link out of it. Remove or rename what is listed, "
            "or choose another output directory. "
            "See docs/user-guide/troubleshooting.md#build-gene-ldscore-index-the-build-directory-cannot-be-used."
        )


# Hidden directory inside the index for the lock, live log, issue tables and staged index.
_BUILD_DIR = ".build"
_LOCK_FILE = "build-gene-ldscore-index.lock"
_LOG_FILE = "build-gene-ldscore-index.log"


def _gene_index_build_dir(index_path: Path) -> Path:
    """Return the hidden build directory inside the index directory."""
    return index_path / _BUILD_DIR


def _chromosome_log_dir(build_dir: Path) -> Path:
    """Return the live chromosome-log directory; success moves it to ``diagnostics/chromosome-logs``."""
    return build_dir / Path(CHROMOSOME_LOG_DIR).name


def _archive_gene_index_log(log_path: Path, history: Path) -> Path:
    """Move one closed operational log to a collision-resistant history path."""
    timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S.%fZ")
    archived = history / f"build-gene-ldscore-index.{timestamp}.log"
    os.replace(log_path, archived)
    return archived


def _prepare_gene_index_log(index_path: Path, *, earlier_attempt: bool) -> Path:
    """Archive a failed attempt's logs and issue tables, discard its scratch, and return the live-log path.

    Only a failed or interrupted attempt leaves a build directory; its staged
    index and workspaces are never reused. ``earlier_attempt`` is the
    preflight's finding that a build directory, empty or recording an earlier
    attempt, was there before this run; without it nothing in the build
    directory is archived or removed.
    """
    build_dir = _gene_index_build_dir(index_path)
    build_dir.mkdir(parents=True, exist_ok=True)
    log_path = build_dir / _LOG_FILE
    if not earlier_attempt:
        return log_path
    chromosome_logs = _chromosome_log_dir(build_dir)
    if chromosome_logs.is_dir() and not chromosome_logs.is_symlink():
        history = build_dir / "history"
        history.mkdir(exist_ok=True)
        timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S.%fZ")
        os.replace(chromosome_logs, history / f"{chromosome_logs.name}.{timestamp}")
    for entry in build_dir.iterdir():
        if entry.is_dir() and not entry.is_symlink() and entry.name != "history":
            shutil.rmtree(entry)
    if log_path.exists():
        history = build_dir / "history"
        history.mkdir(exist_ok=True)
        _archive_gene_index_log(log_path, history)
    plink_issues = build_dir / "plink_input_issues.tsv"
    if plink_issues.exists():
        history = build_dir / "history"
        history.mkdir(exist_ok=True)
        timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S.%fZ")
        os.replace(plink_issues, history / f"plink_input_issues.{timestamp}.tsv")
    issues_path = build_dir / "gene_coordinate_catalog_issues.tsv.gz"
    if issues_path.exists():
        history = build_dir / "history"
        history.mkdir(exist_ok=True)
        timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S.%fZ")
        os.replace(
            issues_path,
            history / f"gene_coordinate_catalog_issues.{timestamp}.tsv.gz",
        )
    return log_path


def _write_gene_catalog_issues(index_path: Path, issues: pd.DataFrame) -> Path:
    """Write the complete current-attempt catalog repair audit in the build directory."""
    path = _gene_index_build_dir(index_path) / "gene_coordinate_catalog_issues.tsv.gz"
    path.parent.mkdir(parents=True, exist_ok=True)
    issues.to_csv(path, sep="\t", index=False, na_rep="", compression="gzip")
    return path


def _finalize_gene_index_log(live_log_path: Path, published_path: Path) -> Path:
    """Move the closed successful logs into the published index diagnostics and remove the build directory.

    Log placement occurs after the workflow context writes its ``Finished``
    footer and closes the file handler; the chromosome logs follow the main
    log into ``diagnostics/chromosome-logs``. Failure is warning-only because the
    scientific index has already been published and reload-validated; the live
    logs and their build directory then remain for inspection and later manual movement.
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
    chromosome_logs = _chromosome_log_dir(live_log_path.parent)
    try:
        if chromosome_logs.is_dir():
            os.replace(chromosome_logs, published_path / CHROMOSOME_LOG_DIR)
    except OSError as exc:
        warnings.warn(
            "The gene LD-score index was published and reload-validated, but its chromosome logs could not "
            f"be moved into {published_path / CHROMOSOME_LOG_DIR} ({type(exc).__name__}: {exc}). The published "
            f"index remains valid; the chromosome logs are retained at {chromosome_logs}.",
            UserWarning,
            stacklevel=2,
        )
        return final_log_path
    try:
        shutil.rmtree(live_log_path.parent)
    except OSError as exc:
        warnings.warn(
            "The gene LD-score index was published and reload-validated, but its build directory "
            f"could not be removed ({type(exc).__name__}: {exc}): {live_log_path.parent}. The published "
            "index remains valid; the next build of this index removes the directory.",
            UserWarning,
            stacklevel=2,
        )
    return final_log_path


def _is_owned_index_entry(entry: Path) -> bool:
    """Return whether one top-level entry is the builder's: build directory, marker file, or index component.

    The preflight has already refused a nonempty build directory without the builder's lock or log file.
    """
    if entry.is_symlink():
        return False
    if entry.name == _BUILD_DIR:
        return entry.is_dir()
    if entry.name == FAILURE_MARKER:
        return entry.is_file()
    return entry.name in _INDEX_COMPONENTS


def _is_build_leftover_only(path: Path) -> bool:
    """Recognize a failed build's destination: only the build directory and a regular ``RUN_FAILED.txt``.

    The preflight has already refused a nonempty build directory without the builder's lock or log file.
    """
    if not path.is_dir():
        return False
    for entry in path.iterdir():
        if entry.is_symlink():
            return False
        if entry.name == _BUILD_DIR and entry.is_dir():
            continue
        if entry.name == FAILURE_MARKER and entry.is_file():
            continue
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
        f"published-reload={publication['published_reload_validation']}."
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


def _load_builder_genetic_map(config: GeneLDScoreIndexBuildConfig, genetic_map_sources: dict[str, str | None]):
    """Load the genetic map supplied for the configured build, or return ``None`` when none is supplied."""
    build = config.genome_build
    other = "hg38" if build == "hg19" else "hg19"
    sources = genetic_map_sources[build]
    if genetic_map_sources[other]:
        raise LDSCInputError(
            f"Gene-index genetic map declares {other}, but the provided genome build is {build}. "
            f"Supply --genetic-map-{build}-sources with verified {build} coordinates. "
            "No coordinates or build labels were changed."
        )
    if not sources:
        return None
    if config.ld_wind_cm is None:
        raise LDSCInputError("Genetic-map sources require --ld-wind-cm. Omit the map for kb or SNP-count windows.")
    check = GenomeBuildCheck(build, context="gene-index genetic map", logger=LOGGER)
    check.paths("genetic map", sources)
    check.resolve()
    return load_genetic_map_group(split_cli_path_tokens(sources), command="build-gene-ldscore-index")


def _require_map_chromosomes(genetic_map: pd.DataFrame | None, chromosomes) -> pd.DataFrame | None:
    """Reject a loaded genetic map that has no rows for an index chromosome, before any task starts."""
    if genetic_map is None:
        return None
    missing = sorted(set(chromosomes).difference(genetic_map["CHR"]), key=chrom_sort_key)
    if missing:
        raise LDSCInputError(
            f"build-gene-ldscore-index could not use genetic map: it has no rows for chromosome"
            f"{'s' if len(missing) > 1 else ''} {', '.join(missing)} of the index, so cM positions there "
            "cannot be interpolated. Most likely the map is for the other genome build or covers only some "
            "chromosomes. Pass a map for --genome-build that covers every index chromosome. "
            f"Other causes & fixes: {GENETIC_MAP_TROUBLESHOOTING['build-gene-ldscore-index']}"
        )
    return genetic_map


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


def _genetic_map_identity(
    config: GeneLDScoreIndexBuildConfig,
    genetic_map_sources: dict[str, str | None],
    genetic_map: pd.DataFrame | None = None,
) -> str | dict:
    """Return the content-bound scientific identity of the effective cM map."""
    if config.ld_wind_cm is None:
        return "not_used"
    build = config.genome_build
    sources = split_cli_path_tokens(genetic_map_sources[build])
    if not sources:
        return "bim_cm"
    if genetic_map is None:
        genetic_map = _load_builder_genetic_map(config, genetic_map_sources)
    return {
        "kind": f"explicit_{build}",
        "content_sha256": _canonical_frame_sha256(genetic_map, sort_by=("CHR", "POS")),
    }


def _builder_index_identity(
    config,
    genetic_map_sources,
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
    keep_indices = keep_individual_indices(read_keep_individuals(config.keep_indivs_file), fam)
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
        "ld_window": config.ld_window.record(),
        "genetic_map": _genetic_map_identity(config, genetic_map_sources, genetic_map),
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


def _create_gene_index_stage(destination: Path) -> Path:
    """Create the empty staged index inside the destination's build directory."""
    staged_index = _gene_index_build_dir(destination) / "stage"
    if staged_index.exists():
        shutil.rmtree(staged_index)
    (staged_index / "chromosomes").mkdir(parents=True)
    return staged_index


def _publish_staged_gene_ldscore_index(
    staged_index: Path, destination: Path, *, accepted_entries: frozenset[str],
) -> Path:
    """Move a reload-validated staged index into its destination, root metadata last.

    The output preflight ran before the build lock, so before changing
    anything the destination is checked again: a top-level entry other than
    the build directory that ``accepted_entries`` does not name, such as an
    index another builder published meanwhile or a file a user added, stops
    publication with the destination unchanged. A regular failure marker is
    the builder's own: a concurrent ``--overwrite`` attempt refused by the
    lock writes one, and this build's successful marker scope removes it.

    The root ``metadata.json`` records ``index_id`` and the loader rejects a
    directory without it. Removing the old root metadata first and writing
    the new one last means an interrupted publication leaves an index that
    no reader accepts, never a mixture of two builds. Each staged entry
    replaces its counterpart, with ``diagnostics/`` merged entry by entry;
    payloads are moved, not copied, and other files in the destination are
    kept. The old build's log and ``diagnostics/chromosome-logs/`` are
    removed here; the caller moves the new build's logs in afterwards.
    """
    load_gene_ldscore_index(staged_index)
    gained = sorted(
        entry.name for entry in destination.iterdir()
        if entry.name not in accepted_entries | {_BUILD_DIR}
        and not (entry.name == FAILURE_MARKER and entry.is_file() and not entry.is_symlink())
    )
    if gained:
        raise FileExistsError(
            f"Gene LD-score index output directory changed during the build: {destination} gained "
            f"{', '.join(gained)}, which the output check before the build did not accept. Nothing was "
            "published and the directory's other contents were left as they were. Move those entries "
            "elsewhere or choose another output directory, then rerun."
        )
    (destination / "metadata.json").unlink(missing_ok=True)
    (destination / "diagnostics" / "build-gene-ldscore-index.log").unlink(missing_ok=True)
    _remove_published_entry(destination / CHROMOSOME_LOG_DIR)
    staged = [entry for entry in staged_index.iterdir() if entry.name not in {"metadata.json", "diagnostics"}]
    if (staged_index / "diagnostics").is_dir():
        staged.extend((staged_index / "diagnostics").iterdir())
    for source in staged:
        target = destination / source.relative_to(staged_index)
        _remove_published_entry(target)
        target.parent.mkdir(parents=True, exist_ok=True)
        os.replace(source, target)
    os.replace(staged_index / "metadata.json", destination / "metadata.json")
    load_gene_ldscore_index(destination)
    return destination


def _remove_published_entry(target: Path) -> None:
    """Remove an earlier build's file or directory; a symlink is unlinked, never followed."""
    if target.is_dir() and not target.is_symlink():
        shutil.rmtree(target)
    else:
        target.unlink(missing_ok=True)


def build_plink_index_chromosome(
    chrom: str,
    baseline_bundle: ChromosomeAnnotations,
    ref_panel: PlinkRefPanel,
    ldscore_config: LDScoreConfig,
    *,
    regression_keys: set[str] | RestrictionIdentityKeys,
    regression_regions: RegionIntervals | None,
    gene_intervals: np.ndarray,
    included: np.ndarray,
    padding_bp: int,
    atom_batch_size: int,
    genetic_map: pd.DataFrame | None = None,
) -> tuple[IndexChromosomeData, dict, pd.DataFrame]:
    """Construct exact common and atom payloads from one PLINK chromosome.

    Returns the chromosome payload the index stores, the build evidence the
    coordinator records but the payload does not store, and the PLINK rows
    that SNP identity cleanup dropped. The evidence holds reference rows
    before genotype QC, genotype-QC and MAF removals, regression SNPs before
    the region subtraction (``regression_selected_rows``), each baseline
    column's ``binary`` or ``quantitative`` classification
    (``annotation_types``), the CM source, and ``baseline_content_sha256``,
    the ``index_id`` fingerprint of the baseline values on the retained
    reference SNPs (after genotype QC).

    ``ref_panel`` prepares the chromosome with ``ldscore_config``, which gives
    the LD window, the whole-chromosome policy, the genotype batch size and
    the common-SNP MAF threshold. ``genetic_map`` is a preloaded map for cM
    windows; without one the panel resolves its own. Genotype alignment,
    windowing, regression SNPs, and atom projection are those of direct
    PLINK LD-score preparation.

    The baseline pass is direct ``ldscore``'s ``compute_chromosome``: baseline
    and regression-weight scores share one traversal, and baseline counts,
    overlaps and ``annotation_types`` classifications come from its float64
    ``annotation_statistics`` call. Each atom batch then repeats genotype
    traversal, bounding the active dense atom projection. Sparse float64
    operator blocks remain live until horizontal assembly, so
    ``atom_batch_size`` does not bound the final operator or assembly copies.
    The chromosome baseline matrix is also explicitly materialized in float64
    for the atom statistics.
    """
    snp_batch_size = ldscore_config.snp_batch_size
    common_maf_min = ldscore_config.common_maf_min
    with ref_panel.prepare_chromosome(chrom, baseline_bundle, ldscore_config, genetic_map=genetic_map) as prepared:
        computed = kernel_ldscore.compute_chromosome(
            chrom, prepared, snp_identifier=ref_panel.global_config.snp_identifier,
            snp_batch_size=snp_batch_size, common_maf_min=common_maf_min,
            regression_keys=regression_keys, regression_regions=regression_regions,
        )
        metadata = prepared.metadata
        # A SNP without MAF would silently leave the common-SNP counts the index stores.
        if computed.M_5_50 is None or metadata["MAF"].isna().any():
            raise LDSCInternalError(
                f"Gene LD-score index construction received reference metadata without MAF on chromosome "
                f"{chrom}, so the common-SNP counts the index stores cannot be computed. PLINK panel "
                "preparation always supplies MAF for every SNP; re-run with DEBUG logging and report the traceback."
            )
        common_mask = metadata["MAF"].to_numpy(dtype=np.float64) >= float(common_maf_min)
        output_rows = computed.regression_rows
        baseline = np.asarray(prepared.annotations.read(), dtype=np.float64)
        geno = prepared.reader

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
            geno.rewind()
            block_scores = np.asarray(
                kernel_ldscore.ld_scores_from_genotypes(
                    geno,
                    prepared.block_left,
                    snp_batch_size,
                    annot=atom_block,
                    output_rows=output_rows,
                ),
                dtype=np.float64,
            )
            operator_blocks.append(sparse.csr_matrix(block_scores, dtype=np.float64))
        operator = (
            sparse.hstack(operator_blocks, format="csr", dtype=np.float64)
            if operator_blocks
            else sparse.csr_matrix((len(output_rows), 0), dtype=np.float64)
        )
        operator_blocks = block_scores = None  # the assembled operator owns these values
        operator.sort_indices()
        validate_ldscore_operator(operator, n_rows=len(output_rows), n_atoms=atom_model.n_atoms)

        identity_columns = ["CHR", "SNP", "POS", *[c for c in ("A1", "A2") if c in metadata.columns]]
        baseline_rows = computed.metadata.loc[:, identity_columns].copy()
        baseline_rows[REGRESSION_LD_SCORE_COLUMN] = computed.w_ld[:, 0]
        for column_index, column in enumerate(baseline_bundle.baseline_columns):
            baseline_rows[column] = computed.ld_scores[:, column_index]
        payload = IndexChromosomeData(
            baseline_rows=baseline_rows,
            baseline_count_all=np.asarray(computed.M, dtype=np.float64),
            baseline_count_common=np.asarray(computed.M_5_50, dtype=np.float64),
            baseline_overlap_all=np.asarray(computed.overlap.baseline_block_all, dtype=np.float64),
            baseline_overlap_common=np.asarray(computed.overlap.baseline_block_common, dtype=np.float64),
            total_reference_snps_all=len(metadata),
            total_reference_snps_common=int(common_mask.sum()),
            atom_model=atom_model,
            operator=operator,
            atom_statistics=atom_statistics,
            reference_snp_funnel=prepared.reference_snp_funnel,
        )
        build_evidence = {
            "reference_rows_before_genotype_qc": int(prepared.reference_rows_before_genotype_qc or len(metadata)),
            "genotype_qc_removed": int(prepared.genotype_qc_removed),
            "maf_removed": int(prepared.maf_removed),
            "regression_selected_rows": int(computed.regression_selected_snp_count),
            "annotation_types": dict(computed.annotation_types),
            "cm_source": prepared.cm_source,
            "baseline_content_sha256": _retained_baseline_content_sha256(
                metadata, pd.DataFrame(prepared.annotations.read(), columns=baseline_bundle.baseline_columns),
                ref_panel.global_config.snp_identifier,
            ),
        }
        return payload, build_evidence, prepared.identity_drops


def _in_index_row_order(metadata: pd.DataFrame, annotations: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return baseline rows and their annotation values in ``(POS, SNP)`` order.

    Index rows are stored in that order, which the loader checks, so
    same-position rows (possible only in rsID modes) reach PLINK preparation,
    and SNP-count windows, in that order. Matched rows carry the PLINK
    ``SNP`` label in rsID modes and unique positions in ``chr_pos`` modes.
    """
    keys = ["POS", "SNP"] if "SNP" in metadata.columns else ["POS"]
    order = (metadata.loc[:, keys].astype({key: str for key in keys if key == "SNP"})
             .sort_values(keys, kind="mergesort").index.to_numpy())
    if np.array_equal(order, np.arange(len(order))):
        return metadata, annotations
    return metadata.take(order).reset_index(drop=True), annotations.take(order).reset_index(drop=True)


def _baseline_plink_evidence(funnel, *, grid_rows, matched_rows, baseline_drops, plink_drops,
                             snp_label_discordant_rows) -> dict[str, int]:
    """Return one chromosome's baseline/PLINK match counts for the build diagnostics.

    The counts come from the reference-SNP funnel of PLINK preparation, the
    baseline rows annotation preparation dropped, and the PLINK rows identity
    cleanup dropped. Matching identities are unique on both sides, so the
    rows each side has without a match are the rows the ``annotation_snp_grid``
    step removed from it. A coordinate contradiction stops the build first,
    so ``coordinate_discordant_rows`` is 0.
    """
    def step(funnel_name, step_name):
        rows = funnel.loc[funnel["funnel"].eq(funnel_name) & funnel["step"].eq(step_name)]
        return int(rows["remaining" if step_name == "panel_rows" else "removed"].iloc[0])

    baseline_duplicates = sum(int(frame.reason.eq("duplicate_identity").sum()) for frame in baseline_drops)
    return {
        "baseline_rows": int(grid_rows),
        "plink_rows": step(REFERENCE_PANEL, "panel_rows"),
        "baseline_duplicate_rows_dropped": baseline_duplicates,
        "plink_bim_duplicate_rows_dropped": int(plink_drops.reason.eq("duplicate_identity").sum()),
        "baseline_identity_rows_dropped": sum(len(frame) for frame in baseline_drops),
        "plink_bim_identity_rows_dropped": int(len(plink_drops)),
        "matched_rows": int(matched_rows),
        "baseline_only_rows": step(ANNOTATION_SNP_GRID, ANNOTATION_SNP_GRID),
        "plink_only_rows": step(REFERENCE_PANEL, ANNOTATION_SNP_GRID),
        "coordinate_discordant_rows": 0,
        "snp_label_discordant_rows": int(snp_label_discordant_rows),
    }


def _selected_plink_rows(prefix: str, chrom: str, snp_identifier: str) -> pd.DataFrame:
    """Return one chromosome's PLINK rows after identity cleanup, as PLINK panel preparation selects them."""
    return select_plink_chromosome(prefix, chrom, snp_identifier, stage="gene_index_build_validation").metadata


def _check_baseline_plink_chromosome(baseline_rows: pd.DataFrame, plink_rows: pd.DataFrame, *,
                                     snp_identifier: str, chrom: str) -> int:
    """Check one chromosome's baseline rows against its cleaned PLINK rows before workers start.

    ``baseline_rows`` are the staged baseline shard's rows, which annotation
    preparation already cleaned; ``plink_rows`` are the PLINK rows that panel
    preparation keeps after identity cleanup. Shared rsIDs with different
    coordinates and shared base keys with incompatible usable alleles raise
    ``LDSCInputError`` carrying an ``input_issues`` repair table, and so does
    an empty match. Returns the number of matched positions whose baseline
    ``SNP`` label differs from the PLINK label in ``chr_pos``-family modes;
    the PLINK label is the one published.
    """
    snp_identifier = normalize_snp_identifier_mode(snp_identifier)
    has_alleles = {"A1", "A2"}.intersection(baseline_rows.columns)
    if len(has_alleles) == 1:
        raise LDSCInputError("Baseline annotations must supply both A1 and A2, or omit both allele columns.")
    baseline = _normalize_identity_rows(baseline_rows, label="baseline", chrom=chrom, snp_identifier=snp_identifier)
    if baseline.empty or plink_rows.empty:
        raise LDSCInputError(
            f"Gene LD-score index baseline/PLINK {snp_identifier} intersection is empty on chromosome {chrom} "
            "after duplicate identity cleanup. Check whether the baseline annotations and PLINK panel "
            "use the intended SNP identifiers and matching coordinates, and check both inputs for "
            "duplicate identities. Restore compatible inputs rather than dropping ambiguous SNPs to bypass this check."
        )
    pairs = _baseline_plink_pairs(baseline, plink_rows, mode=snp_identifier)
    issues = _baseline_plink_contradictions(pairs, mode=snp_identifier, chrom=chrom)
    if not issues.empty:
        error = LDSCInputError(
            f"Baseline annotations and PLINK BIM contain {len(issues)} contradictory SNP records on chromosome {chrom}. "
            "Inspect the repair table; supply matching coordinates and alleles, or correct erroneous input records "
            "only after verification. No SNP labels, alleles, or coordinates were changed.")
        error.input_issues = issues
        raise error
    # Without a contradiction, every shared base key is a match under the run's matching mode.
    if pairs.empty:
        raise LDSCInputError(
            f"Gene LD-score index baseline/PLINK {snp_identifier} intersection is empty on chromosome {chrom}. "
            "Check that baseline annotations and the PLINK panel describe the same SNP population, "
            "use compatible identifiers, and have matching coordinates; supply compatible inputs."
        )
    discordant = 0
    if identity_mode_family(snp_identifier) == "chr_pos" and "SNP" in baseline.columns:
        discordant = int(pairs.SNP_plink.astype(str).ne(pairs.SNP_baseline.astype(str)).sum())
        if discordant:
            LOGGER.info(
                f"Baseline/PLINK chr_pos intersection on chromosome {chrom} has "
                f"{discordant} SNP label disagreement row(s); "
                "PLINK SNP labels are authoritative and will be published."
            )
    return discordant


def _validate_baseline_plink_inputs(sources, prefixes, config, chromosomes, build_dir) -> dict[str, int]:
    """Collect metadata conflicts across chromosomes before starting LD workers.

    Each chromosome's PLINK rows are selected, and their identity cleanup
    logged, once; in rsID-family modes the same rows also feed the
    cross-chromosome rsID join. Returns each chromosome's count of matched
    positions whose baseline and PLINK ``SNP`` labels differ
    (``chr_pos``-family modes; 0 otherwise).
    """
    issues = []
    causes = []
    label_discordant = {}
    rsid_family = identity_mode_family(config.snp_identifier) == "rsid"
    with ExitStack() as stack:
        coordinates = stack.enter_context(_rsid_coordinate_database(build_dir)) if rsid_family else None
        for chrom in chromosomes:
            metadata = sources.shards[chrom].metadata()
            try:
                plink_rows = _selected_plink_rows(prefixes[chrom], chrom, config.snp_identifier)
                if coordinates is not None:
                    _add_rsid_coordinates(coordinates, metadata, plink_rows)
                label_discordant[chrom] = _check_baseline_plink_chromosome(
                    metadata, plink_rows, snp_identifier=config.snp_identifier, chrom=chrom)
            except LDSCInputError as exc:
                causes.append(exc)
                records = getattr(exc, "input_issues", pd.DataFrame([{"chrom": chrom, "reason": "invalid_intersection",
                                                                      "details": str(exc)}]))
                issues.append(records.assign(baseline_sources=",".join(config.baseline_annot_sources),
                                             plink_source=str(prefixes[chrom]) + ".bim"))
        cross_chromosome = None if coordinates is None else _cross_chromosome_rsid_issues(coordinates, prefixes)
    if cross_chromosome is not None and not cross_chromosome.empty:
        issues.append(cross_chromosome.assign(baseline_sources=",".join(config.baseline_annot_sources),
                                              reason="cross_chromosome_coordinate_mismatch"))
        causes.append(LDSCInputError("Shared baseline/PLINK rsIDs declare different chromosomes."))
    path = build_dir / "baseline_plink_issues.tsv"
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
    return label_discordant


@contextmanager
def _rsid_coordinate_database(build_dir):
    """Hold the on-disk database that joins chromosome declarations without retaining genome-wide metadata.

    The database is in the build directory. Each table's ``SNP`` index exists
    before its rows are inserted, so SQLite maintains it as rows arrive and
    the join probes it directly; creating it afterwards would sort the
    genome-wide table in a temporary file in the system temporary directory.
    """
    with tempfile.TemporaryDirectory(prefix="identity-", dir=build_dir) as temporary:
        with closing(sqlite3.connect(str(Path(temporary) / "coordinates.sqlite"))) as database:
            for table in ("baseline", "plink"):
                database.execute(f"CREATE TABLE {table} (SNP TEXT, CHR TEXT, POS INTEGER)")
                database.execute(f"CREATE INDEX {table}_snp ON {table} (SNP)")
            yield database


def _add_rsid_coordinates(database, baseline, plink):
    """Append one chromosome's baseline and cleaned PLINK coordinates to the rsID join."""
    for table, frame in (("baseline", baseline), ("plink", plink)):
        frame.loc[:, ["SNP", "CHR", "POS"]].astype({"CHR": str}).to_sql(
            table, database, if_exists="append", index=False, chunksize=1000)


def _cross_chromosome_rsid_issues(database, prefixes):
    """Return shared baseline/PLINK rsIDs whose declared chromosomes differ."""
    issues = pd.read_sql_query(
        "SELECT b.SNP AS base_key, b.CHR AS CHR_baseline, b.POS AS POS_baseline, "
        "p.CHR AS CHR_plink, p.POS AS POS_plink FROM baseline b JOIN plink p ON b.SNP=p.SNP "
        "WHERE b.CHR != p.CHR", database)
    if not issues.empty:
        issues["plink_source"] = issues.CHR_plink.map(lambda chrom: str(prefixes[chrom]) + ".bim")
    return issues


def _baseline_plink_pairs(baseline, plink, *, mode):
    """Pair baseline and PLINK rows that share a base key; both sides have unique base keys after cleanup."""
    left = baseline.copy()
    right = plink.copy()
    for frame in (left, right):
        frame["base_key"] = base_key_series(frame, mode, context="baseline/PLINK consistency")
        frame["allele_set"] = (allele_set_series(frame, context="baseline/PLINK consistency")[0]
                               if {"A1", "A2"}.issubset(frame.columns) else pd.NA)
    return left.merge(right, on="base_key", suffixes=("_baseline", "_plink"), validate="one_to_one")


def _baseline_plink_contradictions(pairs, *, mode, chrom):
    """Report conflicting coordinates or usable alleles for shared base keys."""
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
    required = {"CHR", "POS"} | ({"SNP"} if identity_mode_family(snp_identifier) == "rsid" else set())
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


