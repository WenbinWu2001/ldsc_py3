"""Resolve annotation declarations, query names and the projection build.

``annotate``, ``ldscore`` and ``build_annotation_shards`` (``AnnotationBuilder``
and ``quantile-h2``) share these checks; each workflow keeps its own option
checks and product steps.
"""

import logging
from collections import Counter
from dataclasses import replace
from pathlib import Path

import pandas as pd

from .parsing import normalize_annotation_chunk
from .interval_queries import bed_query_name
from .query_outcomes import gene_list_gate_a_message
from .._core.chromosomes import AUTOSOMES
from .._core.columns import normalize_genome_build
from .._core.snp_identity import identity_mode_family
from .._genome.build_policy import GenomeBuildCheck
from .._runtime.preflight import InputGate, input_issue, inspect_declared_inputs
from .._runtime.paths import (
    normalize_path_token, normalize_path_tokens, resolve_file_group, split_cli_path_tokens,
)
from ..errors import LDSCInputError, LDSCUsageError, LDSCUserError

LOGGER = logging.getLogger('LDSC.annotation')


def resolve_annotation_inputs(tokens, *, role='baseline'):
    """Return concrete paths, declared chromosomes, and accumulated path issues."""
    paths, declared, issues = [], {}, []
    for raw in split_cli_path_tokens(tokens):
        token = normalize_path_token(raw)
        members = [(token.replace('@', chrom), chrom) for chrom in AUTOSOMES] if '@' in token else [(token, '')]
        for path, chrom in members:
            try:
                matches = resolve_file_group([path], label=f'{role} annotation')
                for match in matches:
                    if match not in paths:
                        paths.append(match)
                    if chrom:
                        declared[match] = chrom
            except (OSError, LDSCUserError) as exc:
                issues.append(input_issue(role, path, chrom, 'missing_required_input', exc))
    return paths, declared, issues


def inspect_annotation_inputs(*, baseline=(), query=(), mode='rsid', strict_annotations=True, **kwargs):
    """Run the declared-input gate with baseline and query annotation checks first.

    Annotation ``@`` completeness remains a workflow-specific policy. Reading
    one annotation row verifies headers and records width without a full scan.
    The other keyword arguments go to ``inspect_declared_inputs``, whose result
    gains the annotation ``declarations``, ``widths`` and ``columns``.
    """
    def check(gate, progress, resolved):
        declarations, widths, columns = {}, {}, {}
        for role, tokens in (('baseline', baseline), ('query', query)):
            if strict_annotations:
                paths, declared, issues = resolve_annotation_inputs(tokens, role=role)
                declarations.update(declared)
                gate.issues.extend(issues)
            else:
                paths = []
                for token in split_cli_path_tokens(normalize_path_tokens(tokens)):
                    matches = gate.check(role, token, lambda token=token: resolve_file_group([token], allow_chromosome_suite=True))
                    paths.extend(matches or [])
            resolved[role] = paths
            for path in paths:
                def header():
                    frame = pd.read_csv(path, sep=r'\s+', nrows=1)
                    if frame.empty:
                        raise LDSCInputError(f"Annotation file '{path}' has no SNP rows. Supply a nonempty annotation SNP grid.")
                    metadata, values = normalize_annotation_chunk(frame, path, mode, log_ignored_metadata=False)
                    member = declarations.get(path)
                    if member and set(metadata.CHR.astype(str)) != {member}:
                        raise LDSCInputError(
                            f"Annotation file '{path}' selected by @ for chromosome {member} starts with "
                            f"chromosome {metadata.CHR.iloc[0]}. "
                            "Restore the matching chromosome file or correct the path template; "
                            "renaming a file does not change its SNP coordinates."
                        )
                    columns[str(path)] = tuple(values.columns)
                    return len(frame.columns)
                width = gate.check(role, path, header, chrom=declarations.get(path, ''))
                if width is not None:
                    widths[str(path)] = width
                progress.advance(object=str(path))
            declared_headers = [(path, columns[path]) for path in paths if path in declarations and path in columns]
            if declared_headers:
                expected = declared_headers[0][1]
                for path, names in declared_headers[1:]:
                    if names != expected:
                        gate.issues.append(input_issue(role, path, declarations[path], 'invalid_required_input',
                            'Annotation columns differ across declared chromosomes. Use the same annotation header for every shard.'))
        return dict(declarations=declarations, widths=widths, columns=columns)

    return inspect_declared_inputs(declaration_check=check, **kwargs)


# SNP metadata columns and the gene-list control column; a query with one of
# these names would duplicate a column of the annotation files.
RESERVED_QUERY_NAMES = frozenset({'gene_control', 'CHR', 'POS', 'BP', 'SNP', 'CM', 'A1', 'A2', 'MAF'})


def annotation_columns(declared):
    """Return the annotation columns of the declared baseline and query files."""
    return {column for path in (*declared.baseline, *declared.query) for column in declared.columns.get(str(path), ())}


def require_usable_bed_query_names(paths, existing_columns, issues_path=None):
    """Stop before any BED file is parsed when a query name cannot become a column.

    A BED query is named by its file name without the last suffix. Each file
    whose name is reserved (``RESERVED_QUERY_NAMES``), is already an annotation
    column in ``existing_columns``, or is shared with another BED file of the
    run gets one ``annotation_name_collision`` row, and all rows are reported
    in one ``LDSCInputError`` (written to ``issues_path`` when given).
    """
    names = {Path(path): bed_query_name(path) for path in paths}
    counts = Counter(names.values())
    issues = []
    for path, name in names.items():
        if name in RESERVED_QUERY_NAMES:
            problem = 'is a reserved SNP metadata or gene-list control column name'
        elif name in existing_columns:
            problem = 'is already an annotation column of the supplied annotation files'
        elif counts[name] > 1:
            problem = 'is also the query name of another BED file of this run'
        else:
            continue
        issues.append(dict(input_issue('query BED', path, '', 'annotation_name_collision',
                                       f'its query name {name!r} {problem}'),
                           repair='Rename the BED file; its name without the last suffix becomes the query annotation name'))
    InputGate('BED query name validation', issues_path, issues,
              troubleshooting='docs/user-guide/troubleshooting.md#common-a-bed-query-name-is-not-usable').finish()


def validate_gene_query_names(batch, existing_columns):
    """Record gene-list query names that would collide with output columns (Gate A).

    A query named like one of ``existing_columns`` (the annotation columns of
    the supplied files), or a focal query named like a reserved metadata column
    or ``gene_control``, gets the ``annotation_name_collision`` reason, so the
    run stops at Gate A together with any rejected identifiers instead of after
    LD-score calculation.
    """
    summary = batch.summary.copy()
    invalid = summary['query'].isin(existing_columns) | (
        summary.input_role.eq('focal') & summary['query'].isin(RESERVED_QUERY_NAMES))
    if invalid.any():
        summary.loc[invalid, 'source_status'] = 'error'
        summary.loc[invalid, 'source_reasons'] = summary.loc[invalid, 'source_reasons'].fillna('').map(
            lambda reason: ';'.join(filter(None, (reason, 'annotation_name_collision'))))
        batch = replace(batch, summary=summary, has_fatal_gate_a_issues=True)
    return batch


def gene_list_gate_a_error(batch):
    """Return the Gate A error for a gene-list batch with fatal issues.

    The error keeps ``batch`` as ``gene_list_batch`` so a caller that writes
    diagnostics after the failure can write the gene-list audit and summary.
    """
    error = LDSCInputError(gene_list_gate_a_message(batch) +
                           ' Other causes & fixes: docs/user-guide/troubleshooting.md#annotate-gene-list-preflight')
    error.gene_list_batch = batch
    return error


def resolve_projection_build(requested, mode, catalog, *, context, baseline_option, baseline_tokens=(),
                             bim_tokens=(), query_tokens=(), bed_tokens=()):
    """Reconcile grid and catalog declarations into the projection genome build.

    An rsID-family run without a gene catalog or a requested build needs no
    projection build and returns ``None``; with a catalog an omitted build
    means ``auto``. Otherwise every labelled grid and catalog input is checked
    against the requested build, and unlabelled BED files inherit the result.
    ``context`` and ``baseline_option`` name the calling workflow and its
    baseline option in build-consistency reports.
    """
    requested = normalize_genome_build(requested)
    if catalog is None and identity_mode_family(mode) == 'rsid' and requested is None:
        return None
    if catalog is not None and requested is None:
        requested = 'auto'
    if requested is None:
        raise LDSCUsageError('Coordinate-family SNP identifier modes require a genome build: supply auto, hg19 or hg38 '
                             '(--genome-build, or genome_build in Python).')
    check = GenomeBuildCheck(requested, context=context, logger=LOGGER)
    check.paths(baseline_option, baseline_tokens)
    check.paths('--query-annot-sources', query_tokens)
    if bim_tokens:
        check.paths('--ref-snp-bim-sources', bim_tokens, bim=True)
    if catalog is not None:
        check.declared('--gene-coordinate-file', getattr(catalog, 'source', 'gene coordinate catalog'), catalog.genome_build)
    check.unchecked('--query-annot-bed-sources', bed_tokens, 'ordinary BED; inherits workflow build')
    return check.resolve()
