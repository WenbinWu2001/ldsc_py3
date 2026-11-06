"""The ``ldscore`` command and its keyword-argument Python workflow ``run_ldscore``.

This module is the public boundary for LD-score path handling. Callers may
provide exact paths, glob patterns, or explicit chromosome-suite tokens using
``@``. ``_configs_from_args`` checks the input-mode rules and builds every
configuration of a run before the output preflight; with ``--overwrite`` a
failed check still leaves ``RUN_FAILED.txt``. A direct run then
preflights its outputs, runs the input gates, resolves the genome build once,
prepares the annotations, measures gene-list support on the prepared reference
SNPs, and hands ``LDScoreCalculator`` (``ldsc._ldscore.calculator``) a
``_RunContext`` that says what it already settled. An indexed run assembles
its query LD scores from a gene LD-score index through ``run_indexed_ldscore``.

For ordinary unpartitioned LD-score runs, callers may omit both baseline and
query annotations. In that case the workflow constructs an all-ones baseline
annotation named exactly ``base`` over the retained reference-panel metadata and
continues through the same calculator and output writer used by partitioned
runs. Query annotations remain partitioned-heritability inputs: they are accepted only
when explicit baseline annotations are supplied.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass, replace as dataclass_replace
import logging
import sys
from typing import Sequence

import numpy as np
import pandas as pd

from .calculator import LDScoreCalculator, _RunContext
from .diagnostics import (
    emit_reference_snp_funnel_notice,
    fail_gene_list_resolution,
    fail_unusable_queries,
    log_query_outcomes,
)
from .direct_batches import _chromosome_annotations, _grid_emptied_message, _slice_annotation_bundle
from .._runtime.cli_help import CLIHelpFormatter, CHROMOSOME_PATH_HELP, SCALAR_PATH_HELP, add_output_dir_argument, add_log_level_argument
from .._runtime.workflow_logging import args_output_scope
from .._runtime.parallelism import (
    WORKERS_HELP, parse_query_batch_size, parse_workers, validate_query_batch_size, validate_workers,
)
from .._runtime.progress import report_phase, PhaseProgress
from .._runtime.scratch import ScratchWorkspace
from .source import LDScoreSource
from .._core.snp_identity import DEFAULT_SNP_IDENTIFIER, SNP_IDENTIFIER_MODES, identity_mode_family
from .._core.columns import normalize_genome_build
from .._core.snp_identity import emit_dropped_snps_notice, normalize_snp_identifier_mode
from ..config import (
    DEFAULT_SNP_BATCH_SIZE,
    GlobalConfig,
    AnnotationBuildConfig,
    LDScoreConfig,
    RefPanelConfig,
    get_global_config,
    print_global_config_banner,
)
from .._genome.hm3 import packaged_hm3_curated_map_path
from .outputs import (
    INPUT_ISSUES_FILE,
    LDSCORE_LOG,
    LDScoreDirectoryWriter,
    LDScoreOutputConfig,
    QueryDiagnostics,
)
from .._runtime.paths import (
    ensure_output_directory,
    normalize_optional_path_token,
    normalize_path_token,
    output_directory_path,
    resolve_scalar_path,
    split_cli_path_tokens,
)
from .._runtime.workflow_logging import (
    log_inputs, log_outputs, materializing_overwrite_guard, overwrite_failure_marker, workflow_logging,
)
from .._refpanel.adapters import resolved_ref_panel
from .._refpanel.funnel import annotation_snp_grid_emptied, read_reference_snp_funnel
from .._core import regions as kernel_regions
from .._core.regions import REGR_SNPS_EXCLUDE_REGIONS_CHOICES, regr_snps_exclude_regions_choice_to_presets
from .._core.identifiers import read_snp_restriction_keys
from .._core.snp_identity import RestrictionIdentityKeys
from ..errors import (
    EmptyReferenceSNPs,
    LDSCConfigError,
    LDSCInputError,
    LDSCUsageError,
    LDSCUserError,
)
from .._annotation.query_outcomes import (
    gene_query_statuses, gene_viability_errors, all_query_annotations_skipped_message,
    bed_query_notice, gene_gate_b_notice, log_notices, print_notices, resolved_only_notice,
)
from .._genome.build_policy import GenomeBuildCheck
from .._annotation.gene_lists import GeneCatalog, GeneCatalogValidationError
from .._annotation.interval_projection import ChromosomeProjector
from .._annotation.gene_query_storage import resolve_gene_lists_staged, persistent_gene_diagnostics
from .._annotation.preflight import (
    annotation_columns,
    require_usable_bed_query_names,
    inspect_annotation_inputs,
)
from .._runtime.input_issues import input_issue
from .._runtime.preflight import InputGate
from .gene_index.scoring import run_indexed_ldscore


LOGGER = logging.getLogger("LDSC.ldscore_calculator")


_QUERY_REQUIRES_BASELINE_MESSAGE = (
    "ldscore cannot run query annotations without baseline annotations. Most likely "
    "a prebuilt, BED, or gene-list query source was supplied without `--baseline-annot-sources`. "
    "Pass matching baseline annotations, or create an explicit all-ones `base` "
    "baseline annotation over the same SNP rows as the query annotations before estimating "
    "partitioned heritability."
)


def build_parser() -> argparse.ArgumentParser:
    """Build the feature parser for LD-score calculation."""
    parser = argparse.ArgumentParser(allow_abbrev=False)
    parser.prog = 'ldsc ldscore'
    parser.formatter_class = CLIHelpFormatter
    parser.description = (
        'Compute LD scores from a reference panel by following Steps 1-6, or assemble gene-list LD scores from a '
        "prebuilt gene LD-score index (see 'Alternative: gene LD-score index')."
    )
    parser.usage = (
        '%(prog)s [--snp-identifier MODE] --genome-build BUILD\n'
        '    [--baseline-annot-sources SOURCES\n'
        '      [--query-annot-bed-sources SOURCES | --query-annot-gene-list-sources SOURCES]]\n'
        '    (--r2-matrix-dir DIR | --plink-prefix PREFIX)\n'
        '    (--ld-wind-cm CM | --ld-wind-kb KB | --ld-wind-snps N)\n'
        '    --output-dir DIR [options]\n'
        '   or: %(prog)s --gene-ldscore-index-dir DIR --query-annot-gene-list-sources SOURCES\n'
        '    --output-dir DIR [options]'
    )
    identity =parser.add_argument_group('Step 1: SNP identity and genome build')
    baseline = parser.add_argument_group('Step 2: Baseline annotations', description='Required when Step 3 adds query annotations. Omit Steps 2 and 3 to compute unpartitioned LD scores when you only need SNP heritability or genetic correlation.')
    queries = parser.add_argument_group('Step 3: Query annotations (optional; choose at most one source type)')
    panel = parser.add_argument_group('Step 4: Reference panel (choose one: R2 matrix or PLINK)')
    plink = parser.add_argument_group('Step 4 continued: PLINK reference panel only (optional)')
    window = parser.add_argument_group('Step 5: LD window (choose one)', description='Specify exactly one window; it must fit within the window of any R2 matrix.')
    output = parser.add_argument_group('Step 6: Output')
    snps = parser.add_argument_group('Optional: regression SNPs and reference SNPs', description='By default, the regression SNPs are the HM3 map SNPs outside MHC and centromeres.')
    indexed = parser.add_argument_group('Alternative: gene LD-score index', description='Use instead of Steps 1, 2, 4, and 5 and the optional SNP settings, which the index stores. Supply --query-annot-gene-list-sources (Step 3) and --output-dir (Step 6); omit --padding-bp, --gene-coordinate-file, and --gene-exclude-regions.')
    computing = parser.add_argument_group('Computing resources')
    logs = parser.add_argument_group('Logging')

    identity.add_argument(
        '--snp-identifier', default=DEFAULT_SNP_IDENTIFIER, choices=SNP_IDENTIFIER_MODES,
        help=(
            'Match SNPs by rsID (rsid) or chromosome and position (chr_pos); allele-aware variants also use '
            'A1/A2. Default: chr_pos_allele_aware; usable alleles are required. With '
            '--gene-ldscore-index-dir, inherit the stored identity and omit this flag.'
        ),
    )
    identity.add_argument(
        '--genome-build', choices=('auto', 'hg19', 'hg37', 'GRCh37', 'hg38', 'GRCh38'), default=None,
        help=(
            'Genome build for coordinate matching, gene intervals and the --regr-snps-exclude-regions presets. '
            'Required for coordinate-family modes; direct gene-list runs default to auto. Auto checks all '
            'coordinate inputs and shards; concrete builds must agree with input declarations. Other rsid runs '
            'need hg19 or hg38 unless --regr-snps-exclude-regions is none. With --gene-ldscore-index-dir, '
            'inherit its build and omit this flag.'
        ),
    )

    baseline.add_argument(
        '--baseline-annot-sources', default=None, metavar='SOURCES',
        help=(
            "Baseline annotation tables: comma-separated paths, '*' or '@' patterns. Required with direct query "
            'inputs. If omitted without queries, use a single all-SNP baseline. Cannot be combined with '
            '--gene-ldscore-index-dir. '
            + CHROMOSOME_PATH_HELP
        ),
    )

    query_group = queries.add_mutually_exclusive_group()
    query_group.add_argument(
        '--query-annot-bed-sources', default=None, metavar='SOURCES',
        help=(
            "BED files whose intervals define query SNPs; accepts comma-separated paths or '*' patterns. Requires "
            '--baseline-annot-sources; cannot be combined with --query-annot-sources, '
            "--query-annot-gene-list-sources, or --gene-ldscore-index-dir. If omitted, add no BED queries. '@' is "
            'not expanded; see --baseline-annot-sources for wildcard and quoting rules.'
        ),
    )
    query_group.add_argument(
        '--query-annot-gene-list-sources', default=None, metavar='SOURCES',
        help=(
            'One-column gene-list files, separated by commas. Direct mode requires --baseline-annot-sources, '
            '--gene-coordinate-file, and --padding-bp; indexed mode requires --gene-ldscore-index-dir. Cannot be '
            'combined with --query-annot-sources or --query-annot-bed-sources. If omitted, add no gene-list '
            "queries. Accepts '*' patterns (see --baseline-annot-sources); '@' is not expanded."
        ),
    )
    query_group.add_argument(
        '--query-annot-sources', default=None, metavar='SOURCES',
        help=(
            "Query annotation tables: comma-separated paths, '*' or '@' patterns (see --baseline-annot-sources for "
            'pattern rules). Requires --baseline-annot-sources; cannot be combined with --query-annot-bed-sources, '
            '--query-annot-gene-list-sources, or --gene-ldscore-index-dir. If omitted, add no table-based '
            'queries.'
        ),
    )
    queries.add_argument(
        '--padding-bp', type=int, default=None, metavar='BP',
        help=(
            'Extend BED intervals or gene boundaries by this many base pairs at both ends, then find covered '
            'SNPs. Omission means 0 for BED input; direct gene-list input requires an explicit value, including '
            '0. Cannot be used with prebuilt annotations, no query, or --gene-ldscore-index-dir. Set to 0 if your '
            'BED intervals are already padded to avoid double padding.'
        ),
    )
    queries.add_argument(
        '--gene-coordinate-file', default=None, metavar='FILE',
        help=(
            'Headered TSV/TSV.GZ of gene coordinates, using one-based inclusive intervals. Required with direct '
            '--query-annot-gene-list-sources; indexed runs use their stored catalog and reject this flag. If '
            'omitted outside gene-list mode, no catalog is used. Supply an exact file path; patterns are not '
            'expanded.'
        ),
    )
    queries.add_argument(
        '--control-gene-list-file', default=None, metavar='FILE',
        help=(
            'One-column gene list added as a fixed control alongside each gene-list query. Requires '
            '--query-annot-gene-list-sources. If omitted, add no gene control. Supply an exact file path; '
            'patterns are not expanded.'
        ),
    )
    queries.add_argument(
        '--gene-exclude-regions', choices=('none', 'mhc'), default='none',
        help=(
            'Exclude genes overlapping MHC before extending their intervals; none excludes no genes. '
            'Default: none for direct gene-list input. With --gene-ldscore-index-dir, use the stored policy '
            'and omit this flag.'
        ),
    )
    queries.add_argument(
        '--allow-unresolved-genes', action='store_const', const='resolved-only',
        dest='gene_list_resolution_policy', default='strict',
        help=(
            'Omit rejected gene identifiers and continue with the usable subset, recording omissions '
            'in diagnostics. Default: off (strict resolution; rejected identifiers stop the run). '
            'Requires --query-annot-gene-list-sources. Other validation errors remain fatal.'
        ),
    )

    panel.add_argument(
        '--r2-matrix-dir', dest='r2_matrix_dir', default=None, metavar='DIR',
        help=(
            'R2 matrix directory written by build-r2-matrix, containing chrN_r2.parquet and required chrN_meta.tsv.gz '
            'files. The LD window must use the unit the R2 matrix was built with and be no wider. Direct mode '
            'requires exactly one of --r2-matrix-dir or --plink-prefix. Cannot be combined with '
            '--gene-ldscore-index-dir; no default.'
        ),
    )
    panel.add_argument(
        '--plink-prefix', default=None, metavar='PREFIX',
        help=(
            'PLINK .bed/.bim/.fam filename prefix or chromosome suite stem. Direct mode requires exactly one of '
            '--plink-prefix or --r2-matrix-dir. Cannot be combined with --gene-ldscore-index-dir; no default. '
            "Direct query runs require every declared chromosome member. Accepts '*' and '@' patterns; see "
            '--baseline-annot-sources for pattern rules.'
        ),
    )

    plink.add_argument(
        '--keep-indivs-file', default=None, metavar='FILE',
        help=(
            'File listing one individual ID (IID) per row to retain from PLINK genotypes; an IID listed more '
            'than once counts once. Requires '
            '--plink-prefix; cannot be used with --r2-matrix-dir or --gene-ldscore-index-dir. If omitted, retain all '
            'individuals. Supply an exact file path; patterns are not expanded.'
        ),
    )
    plink.add_argument(
        '--genetic-map-hg19-sources', default=None, metavar='SOURCES',
        help=(
            'hg19 genetic-map files used to derive centiMorgan coordinates for PLINK SNPs. Accepts '
            "comma-separated exact file paths; '*' and '@' are not expanded. If omitted, use informative PLINK CM "
            'values; cM windows require one of these sources. Direct PLINK mode only.'
        ),
    )
    plink.add_argument(
        '--genetic-map-hg38-sources', default=None, metavar='SOURCES',
        help=(
            'hg38 genetic-map files used to derive centiMorgan coordinates for PLINK SNPs. Accepts '
            "comma-separated exact file paths; '*' and '@' are not expanded. If omitted, use informative PLINK CM "
            'values; cM windows require one of these sources. Direct PLINK mode only.'
        ),
    )
    plink.add_argument(
        '--export-ref-metadata', default=False, action='store_true',
        help=(
            'Write reference-panel SNP metadata, with the minor allele as A1, to ref_metadata/chrN_meta.tsv.gz '
            'under --output-dir. Applies to --plink-prefix only; rejected with --r2-matrix-dir, whose '
            'chrN_meta.tsv.gz sidecars hold the same metadata, and with --gene-ldscore-index-dir. Default: off; '
            'these additional files are not written.'
        ),
    )

    # Not required: indexed mode takes its window from the index.
    window_choice = window.add_mutually_exclusive_group()
    window_choice.add_argument(
        '--ld-wind-cm', default=None, type=float, metavar='CM',
        help=(
            'Maximum LD-window distance on either side of each SNP, in centiMorgans; requires usable '
            'genetic-map coordinates. No default; specify exactly one of --ld-wind-snps, --ld-wind-kb, or '
            '--ld-wind-cm. Direct mode only.'
        ),
    )
    window_choice.add_argument(
        '--ld-wind-kb', default=None, type=float, metavar='KB',
        help=(
            'Maximum LD-window distance on either side of each SNP, in kilobases. No default; specify '
            'exactly one of --ld-wind-snps, --ld-wind-kb, or --ld-wind-cm. Direct mode only.'
        ),
    )
    window_choice.add_argument(
        '--ld-wind-snps', default=None, type=int, metavar='N',
        help=(
            'Maximum LD-window distance on either side of each SNP, counted in retained SNPs with '
            '--plink-prefix and in R2 matrix rows with --r2-matrix-dir. No default; specify exactly one of '
            '--ld-wind-snps, --ld-wind-kb, or --ld-wind-cm. Direct mode only.'
        ),
    )
    window.add_argument(
        '--yes-really', default=False, action='store_true',
        help=(
            'Allow an LD window covering an entire chromosome, bypassing the usual guard. Default: off; '
            'reject whole-chromosome windows. Direct mode only.'
        ),
    )

    add_output_dir_argument(
        output,
        help=(
            'Required destination for LD-score tables, SNP counts, and diagnostics.'
        ),
    )
    output.add_argument(
        '--overwrite', action='store_true', default=False,
        help=(
            "Replace this command's existing output files and remove obsolete outputs from an earlier run. "
            'Default: off; stop if output files already exist.'
        ),
    )

    snps.add_argument(
        '--regr-snps-file', default=None, metavar='FILE',
        help=(
            'Headered SNP identity table selecting the regression SNPs: stored output rows and regression weights. If '
            'omitted, use bundled HapMap3 SNPs. --regr-snps-exclude-regions then applies. Direct mode only; '
            'rejected with --gene-ldscore-index-dir. '
            + SCALAR_PATH_HELP
        ),
    )
    snps.add_argument(
        '--regr-snps-exclude-regions', choices=REGR_SNPS_EXCLUDE_REGIONS_CHOICES, default='mhc-and-centromeres',
        help=(
            'Exclude MHC, centromeres, both, or neither from regression SNPs after --regr-snps-file '
            'or HapMap3 selection. Default: mhc-and-centromeres. Region coordinates follow the run\'s '
            'genome build (see --genome-build). Reference SNPs and annotation counts are unchanged. '
            'Direct mode only; indexed runs inherit this setting.'
        ),
    )
    snps.add_argument(
        '--ref-panel-snps-file', default=None, metavar='FILE',
        help=(
            'Headered SNP identity table restricting the reference SNPs contributing to LD scores. If omitted, '
            'use all otherwise eligible reference SNPs. Extra columns are ignored and duplicate keys collapse. '
            'Direct mode only; rejected with --gene-ldscore-index-dir. Supply an exact file path; patterns are '
            'not expanded.'
        ),
    )
    snps.add_argument(
        '--maf-min', default=None, type=float, metavar='VALUE',
        help=(
            'Keep reference SNPs with MAF at or above this value, in [0, 0.5]. If omitted, apply no '
            'additional MAF threshold. Direct mode only; indexed runs inherit their reference selection.'
        ),
    )
    snps.add_argument(
        '--common-maf-min', default=0.05, type=float, metavar='VALUE',
        help=(
            'Minimum MAF for common-SNP annotation counts, including the boundary. Default: 0.05; does not '
            'filter reference SNPs. Direct mode only; indexed runs inherit their count threshold.'
        ),
    )

    indexed.add_argument(
        '--gene-ldscore-index-dir', default=None, metavar='DIR',
        help=(
            'Use a complete prebuilt gene LD-score index instead of a live reference panel. Requires '
            '--query-annot-gene-list-sources; cannot be combined with direct-input or direct-computation '
            'flags. If omitted, compute LD scores directly.'
        ),
    )

    computing.add_argument(
        '--workers', default=1, type=parse_workers, metavar='N',
        help='Maximum simultaneous chromosomes in direct or indexed mode. ' + WORKERS_HELP,
    )
    computing.add_argument(
        '--query-batch-size', default=100, type=parse_query_batch_size, metavar='N',
        help=(
            'Maximum number of query annotations processed at once. Default: 100; use a smaller positive '
            'value to reduce memory. Each completed batch has its own query file; direct mode repeats reference traversal.'
        ),
    )
    computing.add_argument(
        '--snp-batch-size', default=DEFAULT_SNP_BATCH_SIZE, type=int, metavar='N',
        help=(
            f'Number of reference SNPs processed per PLINK genotype batch. Default: {DEFAULT_SNP_BATCH_SIZE}; '
            'larger batches use more memory. LD scores do not depend on it. Ignored with --r2-matrix-dir; '
            'rejected with --gene-ldscore-index-dir.'
        ),
    )

    add_log_level_argument(logs)
    return parser


@dataclass(frozen=True)
class _DirectConfigs:
    """The configuration of one direct ``ldscore`` run, built before the output preflight.

    ``global_config`` holds the requested genome build, which the run resolves
    against its inputs (``_resolve_build``): ``"auto"`` or a concrete build in
    coordinate-family modes. ``genome_build`` is ``--genome-build`` as given,
    normalized, or ``None``. rsID-family identity records no build, but gene-list
    projection and region-exclusion presets still use one.
    """

    global_config: GlobalConfig
    genome_build: str | None
    annotation_config: AnnotationBuildConfig
    ref_panel_config: RefPanelConfig
    ldscore_config: LDScoreConfig
    output_config: LDScoreOutputConfig
    regr_snps_exclude_regions: str

    @property
    def has_queries(self) -> bool:
        """Return whether the run adds prebuilt, BED, or gene-list query annotations."""
        config = self.annotation_config
        return bool(config.query_annot_sources or config.query_annot_bed_sources
                    or config.query_annot_gene_list_sources)

    @property
    def has_gene_lists(self) -> bool:
        """Return whether the run projects gene-list query annotations."""
        return bool(self.annotation_config.query_annot_gene_list_sources)


@dataclass(frozen=True)
class _IndexedConfigs:
    """One indexed ``ldscore`` run, checked before the output preflight.

    The gene LD-score index supplies the SNP identity, genome build, reference
    and regression SNPs, and LD window; the run adds its gene lists and
    execution settings.
    """

    index_dir: str
    query_gene_list_sources: tuple[str, ...]
    control_gene_list_file: str | None
    gene_list_resolution_policy: str
    query_batch_size: int
    workers: int
    output_dir: str
    overwrite: bool
    log_level: str


@materializing_overwrite_guard(
    args_output_scope,
    command="run_ldscore_from_args(...)",
)
def run_ldscore_from_args(args: argparse.Namespace) -> "LDScoreSource":
    """Run LD-score calculation from a parsed CLI namespace.

    The workflow resolves unified path tokens for baseline annotations, optional
    prebuilt, BED, or gene-list queries, and the reference panel. If no baseline
    or query annotations are supplied, it synthesizes an all-ones ``base``
    annotation over the retained reference-panel metadata. Before calculation
    it preflights the complete baseline/query/overlap artifact family, root
    metadata, applicable diagnostics, and ``diagnostics/ldscore.log`` under
    ``output_dir``. Every run writes ``diagnostics/reference_snp_funnel.tsv``;
    reference-panel SNPs without an annotation row are excluded and reported
    as one log WARNING, which a successful command-line run also prints to
    stderr. Each execution batch is privately written and released;
    successful publication records ordered files and chromosome row groups in
    ``metadata.json.query_batches``. A single batch uses
    ``ldscore.query.parquet``; multiple batches use numbered query files.
    Successful authorized overwrites remove stale owned batch files.

    Direct batches repeat reference traversal. Indexed runs reuse one
    chromosome operator per worker through its query batches. Both modes cap
    chromosome workers at the chromosome count. Scientific restrictions and
    count/overlap definitions are independent of these runtime settings.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed LD-score command arguments, including an output directory and
        the record of explicitly supplied options that ``parse_command_line``
        and ``namespace_from_kwargs`` keep (``_explicit_cli_options``).

    Returns
    -------
    LDScoreSource
        Validated saved-artifact handle with shared baseline values and
        explicit ``read_queries(names)`` access. No query LD-score tables
        remain resident in the returned object. ``output_paths`` contains
        scientific artifacts only.
    """
    configs = _configs_from_args(args)
    if isinstance(configs, _IndexedConfigs):
        return _run_explicit_indexed_ldscore(configs)
    return _run_direct_ldscore(configs)


# Range errors that a configuration words by field name, worded with the option a command-line user typed.
_OPTION_RANGE_ERRORS = {
    "maf_min": (
        "ldscore received invalid `--maf-min={value}`. Most likely the minor-allele frequency threshold was "
        "entered outside the valid [0, 0.5] range. Pass a value between 0 and 0.5, or omit the option."
    ),
    "common_maf_min": (
        "ldscore received invalid `--common-maf-min={value}`. Most likely the common-SNP MAF threshold was "
        "entered outside the valid [0, 0.5] range. Pass a value between 0 and 0.5."
    ),
    "snp_batch_size": (
        "ldscore received invalid `--snp-batch-size={value}`. Most likely the parquet query batch size was set "
        "to zero or a negative value. Pass a positive integer batch size."
    ),
}


def _config_from_options(config_type, **fields):
    """Build ``config_type(**fields)``; a range error on a field in ``_OPTION_RANGE_ERRORS`` names its option."""
    try:
        return config_type(**fields)
    except LDSCConfigError as exc:
        for name, message in _OPTION_RANGE_ERRORS.items():
            if name in fields and str(exc).startswith(f"Could not construct {config_type.__name__}: {name}="):
                raise LDSCConfigError(message.format(value=fields[name])) from None
        raise


# PLINK-only ldscore options, with why an R2 matrix run cannot use them.
_PLINK_ONLY_OPTIONS = {
    "keep_indivs_file": "individuals can be filtered only before PLINK genotype LD calculation",
    "genetic_map_hg19_sources": "the R2 matrix sidecars already hold the CM of every SNP",
    "genetic_map_hg38_sources": "the R2 matrix sidecars already hold the CM of every SNP",
}


def _reject_plink_only_options(args: argparse.Namespace) -> None:
    """Name every PLINK-only option supplied with ``--r2-matrix-dir`` in one command-worded error.

    ``RefPanelConfig`` rejects the same fields for Python callers; this words
    the rule with the options a command-line user typed.
    """
    supplied = [name for name in _PLINK_ONLY_OPTIONS if getattr(args, name)]
    if not supplied:
        return
    options = ", ".join(f"`--{name.replace('_', '-')}`" for name in supplied)
    reasons = "; ".join(dict.fromkeys(_PLINK_ONLY_OPTIONS[name] for name in supplied))
    raise LDSCUsageError(
        f"ldscore cannot use {options} with `--r2-matrix-dir`: {reasons}. "
        f"Remove {'these options' if len(supplied) > 1 else 'this option'}, "
        "or rerun in PLINK mode with `--plink-prefix`."
    )


def _configs_from_args(args: argparse.Namespace) -> "_DirectConfigs | _IndexedConfigs":
    """Check the input-mode rules and build every configuration of one ``ldscore`` run.

    Runs before the output preflight, for the command line and for
    ``run_ldscore()`` alike, but inside the failure-marker scope: with
    ``--overwrite`` (``overwrite=True``) an error raised here still leaves
    ``RUN_FAILED.txt`` in the output directory. An option counts as supplied
    when the caller named it (``_explicit_cli_options``), whatever its value.
    """
    explicit_options = args._explicit_cli_options
    _validate_padding_usage(args)
    _validate_gene_list_mode_args(args, explicit_options)
    if args.gene_ldscore_index_dir is not None:
        return _indexed_configs(args, explicit_options)
    gene_lists = tuple(split_cli_path_tokens(args.query_annot_gene_list_sources))
    if not _has_cli_tokens(args.baseline_annot_sources) and (
        gene_lists or _has_cli_tokens(args.query_annot_sources) or _has_cli_tokens(args.query_annot_bed_sources)
    ):
        raise LDSCUsageError(_QUERY_REQUIRES_BASELINE_MESSAGE)
    r2_matrix_dir = normalize_optional_path_token(args.r2_matrix_dir)
    if bool(r2_matrix_dir) == bool(args.plink_prefix):
        raise LDSCUsageError(
            "ldscore could not choose a reference-panel backend. Most likely both R2 matrix "
            "input and PLINK input were supplied, or neither was supplied. Pass exactly "
            "one of `--r2-matrix-dir` for R2 matrix mode or `--plink-prefix` for PLINK mode."
        )
    snp_identifier = normalize_snp_identifier_mode(args.snp_identifier)
    genome_build = normalize_genome_build(args.genome_build)
    if r2_matrix_dir:
        _reject_plink_only_options(args)
        if identity_mode_family(snp_identifier) == "chr_pos" and genome_build is None:
            raise LDSCUsageError(
                "ldscore cannot run R2 matrix mode with a coordinate-family SNP identifier mode "
                "without a genome build. Most likely `--snp-identifier chr_pos` or "
                "`--snp-identifier chr_pos_allele_aware` was supplied without `--genome-build`. Pass "
                "`--genome-build hg19`, `--genome-build hg38`, or `--genome-build auto`."
            )
    ref_panel_config = _config_from_options(
        RefPanelConfig,
        backend="r2_matrix" if r2_matrix_dir else "plink",
        plink_prefix=None if r2_matrix_dir else args.plink_prefix,
        r2_matrix_dir=r2_matrix_dir or None,
        ref_panel_snps_file=args.ref_panel_snps_file,
        maf_min=args.maf_min,
        keep_indivs_file=args.keep_indivs_file,
        genetic_map_hg19_sources=args.genetic_map_hg19_sources,
        genetic_map_hg38_sources=args.genetic_map_hg38_sources,
    )
    ldscore_config = _config_from_options(
        LDScoreConfig,
        ld_wind_snps=args.ld_wind_snps,
        ld_wind_kb=args.ld_wind_kb,
        ld_wind_cm=args.ld_wind_cm,
        regr_snps_file=args.regr_snps_file,
        snp_batch_size=args.snp_batch_size,
        query_batch_size=args.query_batch_size,
        common_maf_min=args.common_maf_min,
        whole_chromosome_ok=args.yes_really,
        export_ref_metadata=args.export_ref_metadata,
        workers=args.workers,
    )
    global_config = _requested_global_config(
        snp_identifier, genome_build, has_gene_lists=bool(gene_lists), log_level=args.log_level)
    annotation_config = AnnotationBuildConfig(
        baseline_annot_sources=tuple(split_cli_path_tokens(args.baseline_annot_sources)),
        query_annot_sources=tuple(split_cli_path_tokens(args.query_annot_sources)),
        query_annot_bed_sources=tuple(split_cli_path_tokens(args.query_annot_bed_sources)),
        query_annot_gene_list_sources=gene_lists,
        gene_coordinate_file=args.gene_coordinate_file,
        control_gene_list_file=args.control_gene_list_file,
        gene_list_resolution_policy=args.gene_list_resolution_policy,
        gene_exclude_regions=args.gene_exclude_regions,
        padding_bp=0 if args.padding_bp is None else args.padding_bp,
    )
    return _DirectConfigs(
        global_config=global_config,
        genome_build=genome_build,
        annotation_config=annotation_config,
        ref_panel_config=ref_panel_config,
        ldscore_config=ldscore_config,
        output_config=LDScoreOutputConfig(output_dir=normalize_path_token(args.output_dir), overwrite=args.overwrite),
        regr_snps_exclude_regions=args.regr_snps_exclude_regions,
    )


def _requested_global_config(snp_identifier: str, genome_build: str | None, *, has_gene_lists: bool,
                             log_level: str) -> GlobalConfig:
    """Return the run's ``GlobalConfig`` with the requested genome build, before it is resolved.

    Coordinate-family modes need a build; a gene-list run without one infers
    it (``auto``). rsID-family identity records no build: there
    ``--genome-build`` selects only coordinates, the gene projection build
    and the region-exclusion presets. ``GlobalConfig`` rejects ``auto`` in an
    rsID-family run without gene lists, which has nothing to infer it from.
    """
    if identity_mode_family(snp_identifier) == "chr_pos":
        if genome_build is None and not has_gene_lists:
            raise LDSCUsageError(
                "ldscore cannot resolve a coordinate-family SNP identifier mode without a genome build. "
                "Most likely `--snp-identifier` was set to a coordinate-family mode but `--genome-build` "
                "was omitted. Pass `--genome-build auto`, `--genome-build hg19`, or "
                "`--genome-build hg38`."
            )
        return GlobalConfig(snp_identifier=snp_identifier, genome_build=genome_build or "auto", log_level=log_level)
    requested = genome_build if genome_build == "auto" and not has_gene_lists else None
    return GlobalConfig(snp_identifier=snp_identifier, genome_build=requested, log_level=log_level)


def _run_direct_ldscore(configs: _DirectConfigs) -> "LDScoreSource":
    """Run one direct ``ldscore`` request: output preflight, the input gates, then the calculation."""
    from .direct_annotation import prepare_direct_annotations, prepare_synthetic_base

    annotation_config, panel_config = configs.annotation_config, configs.ref_panel_config
    ldscore_config, output_config = configs.ldscore_config, configs.output_config
    output_dir = ensure_output_directory(output_config.output_dir, label="LD-score output directory")
    log_path = output_dir / LDSCORE_LOG
    LDScoreDirectoryWriter.artifact_family(output_dir, include_logs=True).preflight(
        overwrite=output_config.overwrite, additional_paths=[log_path])
    log_path.parent.mkdir(parents=True, exist_ok=True)
    issues_path = output_dir / INPUT_ISSUES_FILE
    with ScratchWorkspace(output_dir) as workspace:
        with workflow_logging("ldscore", log_path, log_level=configs.global_config.log_level):
            from .._refpanel.preflight import inspect_r2_headers
            from .preflight import inspect_reference_inputs
            r2_matrix_dir = panel_config.r2_matrix_dir
            declared_inputs = inspect_annotation_inputs(
                baseline=annotation_config.baseline_annot_sources, query=annotation_config.query_annot_sources,
                files=[("query_annot_bed_sources", annotation_config.query_annot_bed_sources),
                       ("query_annot_gene_list_sources", annotation_config.query_annot_gene_list_sources)],
                scalar_files=[("control_gene_list_file", annotation_config.control_gene_list_file),
                              ("gene_coordinate_file", annotation_config.gene_coordinate_file),
                              ("regr_snps_file", ldscore_config.regr_snps_file),
                              ("ref_panel_snps_file", panel_config.ref_panel_snps_file),
                              ("keep_indivs_file", panel_config.keep_indivs_file)],
                plink=panel_config.plink_prefix, require_complete_suite=configs.has_queries,
                mode=configs.global_config.snp_identifier, issues_path=issues_path,
                checks=[("reference", r2_matrix_dir, lambda: inspect_r2_headers(
                    r2_matrix_dir, configs.genome_build, command="ldscore", ld_window=ldscore_config.ld_window))]
                       if r2_matrix_dir else ())
            require_usable_bed_query_names(declared_inputs.inputs.files["query_annot_bed_sources"],
                                           annotation_columns(declared_inputs), issues_path)
            batch = catalog = None
            if configs.has_gene_lists:
                try:
                    catalog = GeneCatalog.load(annotation_config.gene_coordinate_file)
                except GeneCatalogValidationError as exc:
                    # An input defect found before any calculation: earlier owned outputs stay in place.
                    LDScoreDirectoryWriter.write_input_preflight_diagnostics(
                        QueryDiagnostics({}, catalog_issues=exc.issues), output_dir)
                    raise
                batch = resolve_gene_lists_staged(list(annotation_config.query_annot_gene_list_sources), catalog,
                    workspace, existing_columns=annotation_columns(declared_inputs),
                    control_path=annotation_config.control_gene_list_file,
                    resolution_policy=annotation_config.gene_list_resolution_policy,
                    gene_exclude_regions=annotation_config.gene_exclude_regions)
                if batch.has_fatal_gate_a_issues:
                    fail_gene_list_resolution(batch, output_config)
            global_config, gene_catalog_build = _resolve_build(configs, gene_catalog=catalog)
            print_global_config_banner("run_ldscore_from_args", global_config)
            with PhaseProgress(LOGGER, "validation", "reference content and compatibility"):
                reference_inputs = inspect_reference_inputs(
                    panel_config, global_config, require_autosomes=configs.has_queries,
                    plink_members=declared_inputs.inputs.plink_members)
                gate = InputGate("Reference/control content validation", issues_path,
                                 reference_inputs.issues.to_dict("records"))
                regression_path = _regr_snps_file_from_config(ldscore_config)
                regression_snps = gate.check("regression SNPs", regression_path,
                    lambda: _load_regression_snps(regression_path, global_config,
                        label="packaged HM3 regression SNP map" if ldscore_config.regr_snps_file is None
                        else "regression SNP list"))
                regression_regions = gate.check("regression regions", configs.regr_snps_exclude_regions,
                    lambda: _regression_region_intervals(configs.regr_snps_exclude_regions, _regression_region_build(
                        global_config, gene_catalog_build, configs.genome_build)))
                gate.finish()
            chromosome_prefixes = reference_inputs.scope.get("reference_prefixes_by_chrom")
            scope, skipped_funnels, deferred_queries = {}, [], None
            if annotation_config.baseline_annot_sources:
                annotation_bundle, deferred_queries, scope = prepare_direct_annotations(
                    annotation_config, global_config, workspace, output_config, batch=batch,
                    gene_catalog_build=gene_catalog_build, declared_inputs=declared_inputs,
                    reference_inputs=reference_inputs)
            # The adapter owns the A -> A' restriction; global_config carries the resolved build.
            ref_panel = resolved_ref_panel(global_config, panel_config, chromosome_prefixes=chromosome_prefixes)
            if not annotation_config.baseline_annot_sources:
                annotation_bundle, skipped_funnels = prepare_synthetic_base(ref_panel, global_config, workspace)
            if scope.get("chromosomes"):
                LOGGER.info("Chromosomes resolved and entering the analysis: %s.", ", ".join(scope['chromosomes']))
            log_inputs(output_dir=str(output_dir), reference_mode=ref_panel.spec.backend,
                reference_panel_source=ref_panel.spec.plink_prefix or ref_panel.spec.r2_matrix_dir,
                snp_identifier=global_config.snp_identifier, genome_build=global_config.genome_build,
                ref_panel_snps_file=panel_config.ref_panel_snps_file or "none",
                regr_snps_file=ldscore_config.regr_snps_file or "packaged_hm3_default",
                regression_region_presets=", ".join(regression_regions.source_labels) if regression_regions else "none")
            if gene_catalog_build is not None:
                LOGGER.info("Gene-list catalog projection build: %s (selected from --genome-build and baseline/reference-panel evidence).", gene_catalog_build)
            if annotation_bundle.gene_list_batch is not None:
                annotation_bundle, gate_b_error, gate_b_issues = _apply_direct_gene_gate_b(
                    annotation_bundle, ref_panel, ldscore_config)
                if gate_b_error:
                    fail_unusable_queries(_bundle_diagnostics(annotation_bundle, gate_b_issues), gate_b_error, output_config)
            statuses = annotation_bundle.query_statuses
            if statuses and not annotation_bundle.query_columns:
                fail_unusable_queries(_bundle_diagnostics(annotation_bundle),
                                      all_query_annotations_skipped_message(statuses), output_config)
            context = _RunContext(global_config, ref_panel, main_log_open=True, build_resolved=True,
                                  reference_checked=True)
            result = LDScoreCalculator()._run(annotation_bundle, ldscore_config, output_config=output_config,
                regression_snps=regression_snps, regression_regions=regression_regions, context=context,
                skipped_funnels=skipped_funnels, deferred_queries=deferred_queries)
            log_query_outcomes(result.gene_list_batch, result.query_statuses)
            log_outputs(**result.output_paths)
            log_notices(_query_notices(result))
            return dataclass_replace(result, gene_list_batch=persistent_gene_diagnostics(result.gene_list_batch,
                                                                                         result.output_paths))


def _resolve_build(configs: _DirectConfigs, *, gene_catalog=None) -> tuple[GlobalConfig, str | None]:
    """Resolve the run's genome build once against every coordinate input.

    Returns the run's ``GlobalConfig`` and the gene-list projection build
    (``None`` without gene lists). A coordinate-family run carries the
    resolved build. rsID-family identity records none: a gene-list run
    resolves its projection build (``auto`` when ``--genome-build`` is
    omitted), and another run checks a concrete ``--genome-build`` against
    the input declarations.
    """
    global_config = configs.global_config
    catalog_build = getattr(gene_catalog, "genome_build", None)
    if identity_mode_family(global_config.snp_identifier) == "chr_pos":
        build = _resolve_ldscore_chr_pos_genome_build(configs, global_config.genome_build,
                                                      gene_catalog_build=catalog_build)
        return dataclass_replace(global_config, genome_build=build), build if configs.has_gene_lists else None
    if configs.has_gene_lists:
        return global_config, _resolve_ldscore_chr_pos_genome_build(
            configs, configs.genome_build or "auto", gene_catalog_build=catalog_build)
    if configs.genome_build in {"hg19", "hg38"}:
        _resolve_ldscore_chr_pos_genome_build(configs, configs.genome_build)
    return global_config, None


# Options of a direct run; an indexed run takes their values from its gene LD-score index.
_INDEX_OWNED_OPTIONS = frozenset({
    "--baseline-annot-sources",
    "--query-annot-sources",
    "--query-annot-bed-sources",
    "--padding-bp",
    "--gene-coordinate-file",
    "--gene-exclude-regions",
    "--plink-prefix",
    "--r2-matrix-dir",
    "--snp-identifier",
    "--genome-build",
    "--ref-panel-snps-file",
    "--regr-snps-file",
    "--regr-snps-exclude-regions",
    "--keep-indivs-file",
    "--ld-wind-snps",
    "--ld-wind-kb",
    "--ld-wind-cm",
    "--maf-min",
    "--common-maf-min",
    "--genetic-map-hg19-sources",
    "--genetic-map-hg38-sources",
    "--export-ref-metadata",
    "--snp-batch-size",
    "--yes-really",
})


def _indexed_configs(args: argparse.Namespace, explicit_options) -> _IndexedConfigs:
    """Check the closed indexed gene-list mode: no direct-run option, gene lists, valid execution settings."""
    supplied = sorted(explicit_options & _INDEX_OWNED_OPTIONS)
    if supplied:
        raise LDSCInputError(
            "ldscore indexed mode inherits immutable scientific inputs, SNP identity, and genome build "
            "from the index and cannot accept live overrides: "
            + ", ".join(supplied)
            + ". Remove these options, or remove --gene-ldscore-index-dir to run direct mode."
        )
    gene_lists = split_cli_path_tokens(args.query_annot_gene_list_sources)
    if not gene_lists:
        raise LDSCInputError(
            "ldscore indexed mode requires --query-annot-gene-list-sources."
        )
    return _IndexedConfigs(
        index_dir=args.gene_ldscore_index_dir,
        query_gene_list_sources=tuple(gene_lists),
        control_gene_list_file=args.control_gene_list_file,
        gene_list_resolution_policy=args.gene_list_resolution_policy,
        query_batch_size=validate_query_batch_size(args.query_batch_size),
        workers=validate_workers(args.workers, error_type=LDSCConfigError),
        output_dir=args.output_dir,
        overwrite=args.overwrite,
        log_level=args.log_level,
    )


def _run_explicit_indexed_ldscore(configs: _IndexedConfigs) -> "LDScoreSource":
    """Run one indexed ``ldscore`` request under its own workflow log."""
    output_dir = output_directory_path(configs.output_dir)
    log_path = output_dir / LDSCORE_LOG
    LDScoreDirectoryWriter.artifact_family(output_dir, include_logs=True).preflight(
        overwrite=configs.overwrite, additional_paths=[log_path],
    )
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with workflow_logging("ldscore", log_path, log_level=configs.log_level):
        log_inputs(
            output_dir=str(output_dir),
            reference_mode="gene_ldscore_index",
            gene_ldscore_index=str(configs.index_dir),
            query_gene_list_count=len(configs.query_gene_list_sources),
            control_gene_list_file=configs.control_gene_list_file,
        )
        result = run_indexed_ldscore(
            configs.index_dir,
            query_gene_list_sources=configs.query_gene_list_sources,
            control_gene_list_file=configs.control_gene_list_file,
            gene_list_resolution_policy=configs.gene_list_resolution_policy,
            query_batch_size=configs.query_batch_size,
            workers=configs.workers,
            output_dir=str(output_dir),
            overwrite=configs.overwrite,
            main_log_open=True,
        )
        log_notices(_query_notices(result))
        output_paths = getattr(result, "output_paths", None)
        if output_paths:
            log_outputs(**output_paths)
        return result


def _validate_padding_usage(args: argparse.Namespace) -> None:
    """Reject explicit padding outside live BED or gene-list projection.

    The parser's ``None`` sentinel distinguishes omission from an explicit
    zero; the direct configuration turns an allowed omission into 0.
    """
    has_gene_lists = _has_cli_tokens(args.query_annot_gene_list_sources)
    uses_index = args.gene_ldscore_index_dir is not None
    if args.padding_bp is None:
        if has_gene_lists and not uses_index:
            raise LDSCUsageError(
                "Live gene-list mode requires an explicit --padding-bp value; it must be chosen deliberately. "
                "Use --padding-bp 0 for gene bodies or a positive value for padded intervals."
            )
        return
    if (has_gene_lists or _has_cli_tokens(args.query_annot_bed_sources)) and not uses_index:
        return
    raise LDSCUsageError(
        "ldscore accepts `--padding-bp` only with live `--query-annot-bed-sources` "
        "or `--query-annot-gene-list-sources`. It cannot be used with prebuilt "
        "annotation queries, without a query, or in indexed mode. Remove "
        "`--padding-bp` from the command (or remove the Python `padding_bp` keyword)."
    )


# Options of gene-list mode; --gene-list-resolution-policy is the keyword route's --allow-unresolved-genes.
_GENE_LIST_OPTIONS = frozenset({
    "--gene-coordinate-file",
    "--gene-exclude-regions",
    "--control-gene-list-file",
    "--allow-unresolved-genes",
    "--gene-list-resolution-policy",
})


def _validate_gene_list_mode_args(args: argparse.Namespace, explicit_options) -> None:
    """Reject gene-only options outside their direct or indexed input mode."""
    if not _has_cli_tokens(args.query_annot_gene_list_sources):
        if explicit_options & _GENE_LIST_OPTIONS:
            raise LDSCUsageError(
                "Gene coordinate, exclusion, control-list, and resolution policy options are valid only with "
                "gene-list query annotations. Remove the gene-specific option or use gene-list mode."
            )
        return
    if args.gene_ldscore_index_dir is None and args.gene_coordinate_file is None:
        raise LDSCUsageError(
            "Live gene-list mode requires --gene-coordinate-file. Supply a one-build, one-based "
            "gene-coordinate TSV/TSV.GZ; there is no packaged catalog fallback."
        )


def _apply_direct_gene_gate_b(annotation_bundle, ref_panel, ldscore_config):
    """Measure interval support on the prepared reference SNPs.

    Each gene interval counts the retained reference SNPs it covers with
    ``ChromosomeProjector.counts``, the padding rule of annotation projection.
    The chromosomes are prepared one at a time in the calling process, before
    any LD score is computed.

    Return the updated bundle, the Gate B error message or ``None``, and the
    chromosomes whose reference preparation failed as an input-issue table or
    ``None``.
    """
    batch = annotation_bundle.gene_list_batch
    padding_bp = int(annotation_bundle.source_summary.get("padding_bp", 0))
    intervals_by_index = {}
    for selection in batch.selections:
        intervals_by_index.update(zip(selection.catalog_indices,selection.intervals))
    interval_frame = pd.DataFrame({'catalog_index':list(intervals_by_index), 'interval':list(intervals_by_index.values())})
    interval_coordinates = pd.DataFrame(
        interval_frame.pop("interval").tolist(),
        columns=["chrom", "start0", "end"],
        index=interval_frame.index,
    )
    interval_frame = pd.concat([interval_frame, interval_coordinates], axis=1)
    interval_frame["catalog_index"] = interval_frame["catalog_index"].astype(int)
    support = pd.Series(pd.NA, index=interval_frame["catalog_index"], dtype="Int64")
    errors = []
    input_issues = []
    grid_emptied = []
    for chrom in annotation_bundle.chromosomes:
        genes = interval_frame[interval_frame["chrom"].astype(str).eq(str(chrom))]
        chrom_bundle = _slice_annotation_bundle(annotation_bundle, chrom)
        chrom_annotations = _chromosome_annotations(chrom_bundle,chrom)
        starts, ends = genes["start0"].to_numpy(dtype=np.int64), genes["end"].to_numpy(dtype=np.int64)
        try:
            with ref_panel.prepare_chromosome(chrom, chrom_annotations, ldscore_config) as prepared:
                counts = ChromosomeProjector(prepared.metadata).counts(starts, ends, padding_bp=padding_bp)
        except EmptyReferenceSNPs as exc:
            funnel = getattr(exc, "reference_snp_funnel", None)
            if funnel is not None and annotation_snp_grid_emptied(funnel):
                grid_emptied.append(funnel)
            counts = np.zeros(len(genes), dtype=np.int64)
        except (OSError, EOFError, ValueError, LDSCUserError) as exc:
            errors.append(f"Chromosome {chrom} reference preparation failed: {exc}")
            input_issues.append(input_issue("reference", f"chromosome {chrom}", str(chrom), "invalid_required_input", exc,
                                            repair="Repair reference inputs and runtime filters before retrying."))
            continue
        support.loc[genes["catalog_index"].astype(int).to_numpy()] = counts
    if grid_emptied:
        errors.append(_grid_emptied_message(grid_emptied))
    updated_batch = batch.with_snp_support(support)
    statuses = gene_query_statuses(updated_batch)
    retained_queries = [status.query for status in statuses if status.status in {"ok", "warning"}]
    errors.extend(gene_viability_errors(updated_batch, statuses))
    changes = dict(query_columns=retained_queries, query_statuses=statuses, gene_list_batch=updated_batch)
    return (dataclass_replace(annotation_bundle, **changes), "; ".join(errors) if errors else None,
            pd.DataFrame(input_issues) if input_issues else None)


def _has_cli_tokens(value: str | Sequence[str] | None) -> bool:
    """Return whether a CLI path field contains at least one non-empty token."""
    return bool(split_cli_path_tokens(value))


def run_ldscore(**kwargs) -> "LDScoreSource":
    """Run LD-score calculation from Python using public CLI-style names.

    Keyword arguments are interpreted as CLI-equivalent option names without
    leading ``--``; for example ``baseline_annot_sources``, ``query_annot_sources``,
    ``query_annot_bed_sources``, ``query_annot_gene_list_sources``, ``plink_prefix``, ``r2_matrix_dir``,
    ``keep_indivs_file``, ``snp_batch_size``, ``common_maf_min``, ``workers``, and
    ``output_dir``. Shared runtime assumptions such as ``snp_identifier`` and
    ``genome_build`` must be supplied through ``set_global_config(...)`` first,
    while per-run controls such as ``ref_panel_snps_file`` and
    ``regr_snps_file`` remain ordinary keyword arguments here. ``output_dir``
    is required because this public workflow always writes canonical artifacts.
    ``workers`` defaults to one worker; each worker uses 1 BLAS thread
    in serial and parallel execution. Python scripts must launch workflow
    calls under ``if __name__ == "__main__":``.

    When ``baseline_annot_sources`` and query inputs are omitted, the workflow
    builds a synthetic all-ones baseline column named ``base`` from retained
    reference-panel metadata. ``query_annot_sources`` and
    ``query_annot_bed_sources`` and ``query_annot_gene_list_sources`` require
    explicit baseline annotations because query columns are interpreted
    relative to the baseline annotation rows. Gene-list inputs are projected
    from the required one-based gene coordinate catalog for the resolved build.
    ``padding_bp`` defaults to zero for live BED input but must be chosen
    explicitly for live gene lists; passing it with a prebuilt annotation
    query, no query, or a gene LD-score index raises ``LDSCUsageError`` even when its
    value is zero.

    With ``gene_ldscore_index_dir``, the gene LD-score index supplies the SNP
    identifier mode and genome build, as on the command line; the registered
    ``GlobalConfig.snp_identifier`` and ``genome_build`` are ignored, and the
    returned ``config_snapshot`` records the index values.

    Gene-list resolution defaults to ``gene_list_resolution_policy="strict"``.
    For deliberate subset analysis, pass
    ``gene_list_resolution_policy="resolved-only"`` to omit allowlisted rejected
    rows with diagnostics. Structural input, chromosome-coverage, and requested
    control validation still apply. This Python keyword retains its policy
    strings; the CLI uses the value-free ``--allow-unresolved-genes`` flag to
    select ``"resolved-only"`` and rejects ``--gene-list-resolution-policy``.

    The run's configuration is checked before the output preflight; an error
    found there writes no output, but with ``overwrite=True`` the failed call
    still leaves ``RUN_FAILED.txt`` in ``output_dir``. Keyword errors are
    found earlier and leave no marker.

    Returns
    -------
    LDScoreSource
        Saved artifacts with shared baseline, counts, overlap, diagnostics, and
        canonical paths. Use ``read_queries(names)`` to load selected queries
        across batch files. No query LD tables remain resident after return.
        For a small calculation without writes, use prepared inputs with
        ``LDScoreCalculator.run(..., output_config=None)`` instead.

    Raises
    ------
    LDSCUsageError
        Before anything is written, when a keyword is not an ``ldscore``
        option name, ``output_dir`` is missing, a value is outside the
        command's choices, a flag is not a bool, or ``snp_identifier``,
        ``genome_build`` or ``log_level`` is passed. Before the output
        preflight, when the keywords combine input modes the command rejects;
        a passed keyword counts as supplied even at its default value, so
        ``gene_coordinate_file``, ``gene_exclude_regions``,
        ``control_gene_list_file`` or ``gene_list_resolution_policy``
        without gene lists is rejected.
    LDSCConfigError
        Before the output preflight, when a value is outside its valid range,
        except ``query_batch_size`` in indexed mode.
    ValueError
        Before the output preflight, in indexed mode, when ``query_batch_size``
        is not a positive integer.
    LDSCInputError
        Before the output preflight, in indexed mode, when a direct-run option
        such as ``maf_min`` is passed (the index supplies it) or
        ``query_annot_gene_list_sources`` is missing. After the output
        preflight, when an input is missing, unreadable or inconsistent with
        the others.
    FileExistsError
        Before anything is written, when owned outputs exist and ``overwrite``
        is not ``True``.
    """
    from .._runtime.workflow_kwargs import GLOBAL_CONFIG_RUNTIME_OPTIONS, namespace_from_kwargs

    # Indexed mode ignores the registered identity: its live-override check reads explicit options only.
    args = namespace_from_kwargs(
        "run_ldscore", build_parser(), kwargs, rejected=GLOBAL_CONFIG_RUNTIME_OPTIONS, registered=get_global_config(),
    )
    # The outer marker scope makes a failed overwrite name this function, not run_ldscore_from_args.
    with overwrite_failure_marker(args.output_dir, overwrite=args.overwrite, command="Python API: run_ldscore(...)"):
        return run_ldscore_from_args(args)


def main(argv: Sequence[str] | None = None) -> LDScoreSource:
    """Command-line entry point for the LD-score workflow."""
    from .._runtime.workflow_kwargs import parse_command_line

    args = parse_command_line(build_parser(), list(sys.argv[1:] if argv is None else argv))
    result = run_ldscore_from_args(args)
    print_notices(_query_notices(result))
    _emit_reference_snp_notices(result.output_paths)
    return result


def _resolve_ldscore_chr_pos_genome_build(
    configs: _DirectConfigs,
    genome_build: str,
    *,
    gene_catalog_build: str | None = None,
) -> str:
    """Collect every coordinate input, then resolve or report all build evidence.

    A concrete build checks declarations without fresh inference. Auto trusts declarations, otherwise
    samples each file/chromosome; generic coordinate restrictions must resolve
    independently and BED inherits the result.
    """
    annotations, panel = configs.annotation_config, configs.ref_panel_config
    check = GenomeBuildCheck(genome_build, context="ldscore", logger=LOGGER)
    if gene_catalog_build is not None:
        check.declared("--gene-coordinate-file", annotations.gene_coordinate_file or "gene coordinate catalog",
                       gene_catalog_build)
    check.paths("--baseline-annot-sources", annotations.baseline_annot_sources)
    check.paths("--query-annot-sources", annotations.query_annot_sources)
    check.plink("--plink-prefix", panel.plink_prefix)
    if panel.r2_matrix_dir is not None:
        check.r2("--r2-matrix-dir", panel.r2_matrix_dir)
    mode = configs.global_config.snp_identifier
    check.paths("--ref-panel-snps-file", panel.ref_panel_snps_file, restriction_mode=mode)
    check.paths("--regr-snps-file", configs.ldscore_config.regr_snps_file, restriction_mode=mode)
    check.unchecked("--query-annot-bed-sources", annotations.query_annot_bed_sources,
                    "ordinary BED; inherits workflow build")
    resolved = check.resolve()
    if panel.plink_prefix:
        map_sources = getattr(panel, f"genetic_map_{resolved}_sources")
        if map_sources:
            maps = GenomeBuildCheck(resolved, context="ldscore genetic map", logger=LOGGER)
            maps.paths(f"--genetic-map-{resolved}-sources", map_sources)
            maps.resolve()
    return resolved


@report_phase(LOGGER, 'validation', 'regression SNP restriction')
def _load_regression_snps(path: str, global_config: GlobalConfig, *, label: str) -> RestrictionIdentityKeys:
    """Load the regression SNP list, or the packaged HM3 map, using the active identifier mode."""
    return read_snp_restriction_keys(
        resolve_scalar_path(path, label=label),
        global_config.snp_identifier,
        genome_build=global_config.genome_build,
    )


def _regr_snps_file_from_config(config: LDScoreConfig) -> str:
    """Return the explicit regression list or the canonical bundled HM3 default."""
    return config.regr_snps_file or packaged_hm3_curated_map_path()


def _regression_region_build(global_config: GlobalConfig, gene_catalog_build: str | None,
                             genome_build: str | None) -> str | None:
    """Return the genome build that selects region-exclusion preset BEDs, or ``None`` without one.

    Coordinate-family runs use their resolved build. rsID-family identity
    carries no build, so gene-list runs use their projection build and other
    runs a concrete ``genome_build``, the ``--genome-build`` given.
    """
    for build in (global_config.genome_build, gene_catalog_build, genome_build):
        if build in {"hg19", "hg38"}:
            return build
    return None


def _regression_region_intervals(
    choice: str, build: str | None
) -> kernel_regions.RegionIntervals | None:
    """Load the named intervals of ``choice`` in ``build``, applied only after selecting regression SNPs."""
    presets = regr_snps_exclude_regions_choice_to_presets(choice)
    if not presets:
        return None
    if build is None:
        raise LDSCUsageError(
            "ldscore cannot select named regression exclusion regions without a concrete genome build. "
            "Pass `--genome-build hg19` or `--genome-build hg38`, or use `--regr-snps-exclude-regions none`."
        )
    return kernel_regions.load_preset_intervals(presets, build)


def _query_notices(result) -> list[str | None]:
    """Return a run's gene-list and BED query notices in console order, ``None`` where one does not apply."""
    return [gene_gate_b_notice(result.query_statuses, result.gene_list_batch, command="ldscore"),
            resolved_only_notice(result.gene_list_batch),
            bed_query_notice(result.query_statuses, command="ldscore")]


def _emit_reference_snp_notices(output_paths: dict) -> None:
    """Print the funnel and dropped-SNP WARNINGs to stderr after a successful command-line run."""
    path = output_paths.get("reference_snp_funnel")
    if path is not None:
        emit_reference_snp_funnel_notice(read_reference_snp_funnel(path), path)
    emit_dropped_snps_notice([path for key, path in output_paths.items() if key.startswith("dropped_snps")],
                             command="ldscore")


def _bundle_diagnostics(bundle, input_issues=None) -> QueryDiagnostics:
    """Return the query diagnostics an annotation bundle holds, plus any Gate B input issues, for the writer."""
    return QueryDiagnostics(dict(bundle.source_summary.get("chromosome_scope", {})), tuple(bundle.query_statuses),
                            bundle.gene_list_batch, input_issues)
