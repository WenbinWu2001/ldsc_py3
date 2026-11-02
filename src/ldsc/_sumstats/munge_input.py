# (C) 2014-2019 Brendan Bulik-Sullivan and Hilary Finucane
# Header parsing and aliases derive from LDSC munge_sumstats.py; GPLv3 (LICENSE).
# Python 3 adaptation by Anthony Abrantes, 2024. Modified 2026-09-10:
# extracted input resolution from the munging kernel.
# See NOTICE for attribution and dated modifications.

"""Resolve raw summary-statistics inputs for the munging workflow.

Header aliases, DANER profiles, sample-size selection and bounded source-build
inference belong here. :func:`resolve_raw_schema` reads a raw header once for
both the munging run and ``--infer-only``; the numerical kernel receives only
resolved settings.
"""

import bz2
from dataclasses import dataclass, replace
import gzip
import logging
import warnings

import pandas as pd

from .._core.columns import (
    RAW_SUMSTATS_REQUIRED_OR_OPTIONAL_SPECS,
    RAW_SUMSTATS_SIGNED_STAT_SPECS,
    build_cleaned_alias_lookup,
    clean_header,
)
from ..config import GlobalConfig, MungeConfig
from ..errors import LDSCInputError, LDSCUsageError, LDSCUserError
from .._genome.liftover import SumstatsLiftoverRequest
from .._core.snp_identity import identity_base_mode, identity_mode_family, is_allele_aware_mode
from .munge_kernel import MungeQC, ResolvedMungeInput, mean_info_list_value, prepare_sumstats_restriction

LOGGER = logging.getLogger('LDSC.sumstats_munger')
_BUILD_INFERENCE_CHUNKSIZE = 5_000


COLUMN_HINT_TARGETS = {
    'snp': 'SNP', 'chr': 'CHR', 'pos': 'POS',
    'N': 'N', 'N_col': 'N', 'N_cas': 'N_CAS', 'N_cas_col': 'N_CAS',
    'N_con': 'N_CON', 'N_con_col': 'N_CON', 'a1': 'A1', 'a2': 'A2',
    'p': 'P', 'frq': 'FRQ', 'info': 'INFO', 'nstudy': 'NSTUDY',
}


null_values = {
    'LOG_ODDS': 0,
    'BETA': 0,
    'OR': 1,
    'Z': 0
}

default_cnames = build_cleaned_alias_lookup(
    RAW_SUMSTATS_REQUIRED_OR_OPTIONAL_SPECS + RAW_SUMSTATS_SIGNED_STAT_SPECS
)

describe_cname = {
    'SNP': 'Variant ID (e.g., rs number)',
    'CHR': 'Chromosome',
    'POS': 'Base-pair position',
    'P': 'p-Value',
    'A1': 'A1, the effect allele: the allele the signed statistic refers to.',
    'A2': 'A2, the other allele.',
    'N': 'Sample size',
    'N_CAS': 'Number of cases',
    'N_CON': 'Number of controls',
    'Z': 'Z-score (0 --> no effect; above 0 --> A1 is trait/risk increasing)',
    'OR': 'Odds ratio (1 --> no effect; above 1 --> A1 is risk increasing)',
    'BETA': '[linear/logistic] regression coefficient (0 --> no effect; above 0 --> A1 is trait/risk increasing)',
    'LOG_ODDS': 'Log odds ratio (0 --> no effect; above 0 --> A1 is risk increasing)',
    'INFO': 'INFO score (imputation quality; higher --> better imputation)',
    'FRQ': 'Allele frequency',
    'SIGNED_SUMSTAT': 'Directional summary statistic as specified by --signed-sumstats.',
    'NSTUDY': 'Number of studies in which the SNP was genotyped.'
}

_LIKELY_SIGNED_SUMSTAT_NAMES = {'EFFECT_SIZE', 'EFFECTSIZE', 'LOGOR', 'LOG_OR', 'BETA_HAT'}

# Allele column names whose effect orientation depends on the upstream tool: genome
# REF/ALT, and ALLELE0/1/2, whose effect allele differs between GWAS programs. They
# are never mapped to A1/A2; a file whose only allele columns are these stops the run.
ORIENTATION_AMBIGUOUS_ALLELE_COLUMNS = (
    'REF', 'ALT', 'REFERENCE_ALLELE', 'ALLELE0', 'ALLELE1', 'ALLELE_1', 'ALLELE2', 'ALLELE_2',
)
# Allele column names that say A1 increases the trait, checked against the signed statistic.
_INCREASING_ALLELE_COLUMNS = {'A1': 'INC_ALLELE', 'A2': 'DEC_ALLELE'}
_ALLELE_TROUBLESHOOTING = "docs/user-guide/troubleshooting.md#munge-sumstats-a1-cannot-be-identified"


@dataclass(frozen=True)
class RawSchemaIssue:
    """One header-level problem that stops a munging run.

    ``message`` is the self-contained error the run raises, of class ``error``.
    ``label`` names the problem under ``Missing fields`` in the ``--infer-only``
    report, and ``note`` is the report's explanation, or ``None`` when the
    report explains the problem in its own words.
    """

    label: str
    message: str
    note: str | None
    error: type[LDSCUserError] = LDSCInputError


@dataclass(frozen=True)
class RawSchema:
    """How one raw summary-statistics file's header is read for munging.

    Built by :func:`resolve_raw_schema` for both the munging run and
    ``--infer-only``.

    Attributes
    ----------
    source_path : str
        Resolved raw summary-statistics path.
    headers : tuple of str
        Raw header names in file order.
    detected_format : str
        ``plain``, ``daner-old`` or ``daner-new``.
    inferred_hints : dict of str to str
        Column hints that the format implies; explicit hints override them.
    config : MungeConfig
        The munging config with the detected format, the merged INFO-list
        columns and any ``signed_sumstats`` hint applied.
    column_map : dict of str to str
        Raw column to canonical field; the chosen signed statistic maps to
        ``SIGNED_SUMSTAT``.
    field_descriptions : dict of str to str
        Raw column to a description of the field it supplies, for the log.
    signed_sumstat_name : str
        The raw signed-statistic column, or ``SIGNED_SUMSTATS`` when
        ``--signed-sumstats`` or ``--a1-inc`` decides it.
    signed_sumstat_null : float or None
        Null value of the signed statistic.
    qc : MungeQC
        Sample-size constants, old DANER case and control counts, and the
        effective INFO and MAF thresholds.
    notes : tuple of str
        Observations for ``--infer-only``, such as comma-separated INFO lists.
    issues : tuple of RawSchemaIssue
        Every header-level problem found, in the order the checks run.
    """

    source_path: str
    headers: tuple[str, ...]
    detected_format: str
    inferred_hints: dict[str, str]
    config: MungeConfig
    column_map: dict[str, str]
    field_descriptions: dict[str, str]
    signed_sumstat_name: str
    signed_sumstat_null: float | None
    qc: MungeQC
    notes: tuple[str, ...]
    issues: tuple[RawSchemaIssue, ...]

    def raise_issues(self) -> None:
        """Raise one error for all issues; a single issue keeps its own message and class."""
        if not self.issues:
            return
        if len(self.issues) == 1:
            issue = self.issues[0]
            raise issue.error(issue.message)
        listed = "\n".join(f"  {index}. {issue.message}" for index, issue in enumerate(self.issues, start=1))
        raise LDSCInputError(
            f"munge-sumstats found {len(self.issues)} problems with the columns of '{self.source_path}'; "
            f"fix all of them before rerunning:\n{listed}"
        )


def _decode_header_line(line):
    """Decode one raw header/comment line from plain or compressed input."""
    if isinstance(line, bytes):
        return line.decode('utf-8')
    return line


def count_leading_sumstats_comment_lines(fh):
    """Count leading raw summary-statistics metadata lines that begin with ``##``."""
    openfunc, _compression = get_compression(fh)
    count = 0
    with openfunc(fh) as handle:
        for raw_line in handle:
            line = _decode_header_line(raw_line)
            if not line.startswith('##'):
                break
            count += 1
    return count


def read_header(fh):
    '''Read the first non-metadata line of a file and return its column names.'''
    skiprows = count_leading_sumstats_comment_lines(fh)
    (openfunc, _compression) = get_compression(fh)
    with openfunc(fh) as handle:
        for _idx in range(skiprows):
            handle.readline()
        line = _decode_header_line(handle.readline())
        return [x.rstrip('\n') for x in line.split()]


def _read_first_data_row(path, headers):
    """Return the first data row keyed by raw header names."""
    openfunc, _compression = get_compression(path)
    skiprows = count_leading_sumstats_comment_lines(path)
    with openfunc(path) as handle:
        for _idx in range(skiprows + 1):
            handle.readline()
        for line in handle:
            line = _decode_header_line(line)
            if line.strip():
                values = line.split()
                return {column: values[idx] for idx, column in enumerate(headers) if idx < len(values)}
    return {}


def get_cname_map(flag, default, ignore):
    """Map cleaned headers, preferring ignores, explicit hints, then aliases."""
    clean_ignore = [clean_header(x) for x in ignore]
    cname_map = {x: flag[x] for x in flag if x not in clean_ignore}
    used_targets = set(cname_map.values())
    cname_map.update(
        {x: default[x] for x in default if x not in clean_ignore + list(flag.keys()) and default[x] not in used_targets})
    return cname_map


def get_compression(fh):
    '''
    Read filename suffixes and figure out whether it is gzipped,bzip2'ed or not compressed
    '''
    if fh.endswith('gz'):
        compression = 'gzip'
        openfunc = gzip.open
    elif fh.endswith('bz2'):
        compression = 'bz2'
        openfunc = bz2.BZ2File
    else:
        openfunc = open
        compression = None

    return openfunc, compression


def _detect_format(clean_headers, requested_format):
    """Return the requested raw format, or detect it from cleaned headers when it is ``auto``."""
    if requested_format != 'auto':
        return requested_format
    if any(column.startswith('FRQ_A_') for column in clean_headers) and any(
            column.startswith('FRQ_U_') for column in clean_headers):
        return 'daner-old'
    if {'NCA', 'NCO'}.issubset(clean_headers) or {'NCAS', 'NCON'}.issubset(clean_headers):
        return 'daner-new'
    return 'plain'


def _format_hints(headers, detected_format):
    """Return the column hints a raw format implies; allele orientation is never inferred."""
    hints = {}
    if detected_format == 'daner-new':
        frq_u = next((column for column in headers if clean_header(column).startswith('FRQ_U_')), None)
        if frq_u is not None:
            hints['frq'] = frq_u
    return hints


def _first_row_info_lists(source_path, headers):
    """Return INFO-alias columns whose first value is a comma-separated list, with a note for each."""
    sample = _read_first_data_row(source_path, headers)
    columns = [column for column in headers
               if default_cnames.get(clean_header(column)) == 'INFO' and ',' in sample.get(column, '')]
    return columns, [f"{column} appears to contain comma-separated per-study INFO values." for column in columns]


def likely_signed_sumstat_column(file_cnames):
    """Return the first column whose name suggests an unrecognized signed effect, or ``None``."""
    return next((column for column in file_cnames if clean_header(column) in _LIKELY_SIGNED_SUMSTAT_NAMES), None)


def _column_hints(columns: dict[str, str], config: MungeConfig):
    '''
    Parse flags that specify how to interpret nonstandard column names.

    Returns the flags, a dict that maps (cleaned) arguments to internal column
    names, the null value of ``--signed-sumstats``, and a usage message when
    that option cannot be parsed.
    '''
    flag_cnames = {clean_header(source): target for target, source in columns.items()}
    flag_cnames.update({clean_header(source): 'INFO' for source in config.info_list_columns})

    null_value, problem = None, None
    if config.signed_sumstats_spec:
        try:
            cname, null_value = config.signed_sumstats_spec.split(',')
            null_value = float(null_value)
            flag_cnames[clean_header(cname)] = 'SIGNED_SUMSTAT'
        except ValueError:
            null_value = None
            problem = (
                f"munge-sumstats could not parse --signed-sumstats={config.signed_sumstats_spec!r}. "
                "Most likely the value is missing '<column>,<null_value>' or the null value is not numeric. "
                "Use a value such as --signed-sumstats BETA,0 or --signed-sumstats OR,1."
            )

    return flag_cnames, null_value, problem


def _explicit_sample_size_problem(columns):
    """Return why the explicit direct-N and case/control column hints conflict, or ``None``."""
    has_direct_n = columns.get("N") is not None
    has_n_cas = columns.get("N_CAS") is not None
    has_n_con = columns.get("N_CON") is not None
    if has_n_cas != has_n_con:
        missing = '--N-con-col' if has_n_cas else '--N-cas-col'
        return (
            f"munge-sumstats --N-cas-col and --N-con-col must be provided together; {missing} is missing. "
            "Pass both case/control column names, or drop the incomplete case/control hint and use --N-col."
        )
    if has_direct_n and has_n_cas:
        return (
            "munge-sumstats --N-col cannot be combined with --N-cas-col and --N-con-col. "
            "Choose one sample-size strategy: a direct per-variant N column, or a paired case/control column strategy."
        )
    return None


def _warn_sample_size_column_suppression(message):
    """Emit one visible and logged warning for explicit sample-size precedence."""
    warnings.warn(message, UserWarning, stacklevel=3)
    LOGGER.warning(message)


def _resolve_sample_size_column_strategy(columns, cname_translation, source_path):
    """Apply explicit sample-size precedence; return why automatic inference is ambiguous, or ``None``."""
    if columns.get("N") is not None:
        suppressed = [
            source for source, target in cname_translation.items()
            if target in {'N_CAS', 'N_CON'}
        ]
        for source in suppressed:
            del cname_translation[source]
        if suppressed:
            _warn_sample_size_column_suppression(
                f"--N-col {columns.get('N')} selected the direct-N sample-size strategy; ignored automatically inferred "
                f"N_CAS/N_CON column(s): {', '.join(suppressed)}."
            )
        return None

    if columns.get("N_CAS") is not None:
        suppressed = [
            source for source, target in cname_translation.items()
            if target == 'N'
        ]
        for source in suppressed:
            del cname_translation[source]
        if suppressed:
            _warn_sample_size_column_suppression(
                f"--N-cas-col {columns.get('N_CAS')} and --N-con-col {columns.get('N_CON')} selected the case/control "
                f"sample-size strategy; ignored automatically inferred direct-N column(s): {', '.join(suppressed)}."
            )
        return None

    sources_by_target = {
        target: [source for source, mapped_target in cname_translation.items() if mapped_target == target]
        for target in ('N', 'N_CAS', 'N_CON')
    }
    if all(sources_by_target[target] for target in ('N', 'N_CAS', 'N_CON')):
        direct = sources_by_target['N'][0]
        n_cas = sources_by_target['N_CAS'][0]
        n_con = sources_by_target['N_CON'][0]
        return (
            f"munge-sumstats found multiple sample-size strategies through automatic inference in '{source_path}': "
            f"direct N column '{direct}' and case/control columns '{n_cas}'/'{n_con}'. "
            "LDSC3 will not choose one silently because the values can differ. Choose one explicitly with "
            f"--N-col {direct} or --N-cas-col {n_cas} --N-con-col {n_con}."
        )
    return None


# Legacy LDSC QC defaults; each applies only when its column is selected.
DEFAULT_INFO_MIN = 0.9
DEFAULT_MAF_MIN = 0.01
_QC_THRESHOLDS = (
    ("INFO", "info_min", "--info-min", DEFAULT_INFO_MIN, "an INFO column", "INFO",
     "imputation-quality", "--info COLUMN (or --info-list COLUMN for comma-separated per-study values)"),
    ("FRQ", "maf_min", "--maf-min", DEFAULT_MAF_MIN, "an allele-frequency column", "MAF",
     "allele-frequency", "--frq COLUMN"),
)


def resolve_qc_thresholds(config, mapped_fields):
    """Resolve effective INFO/FRQ thresholds against the selected columns.

    A ``None`` threshold applies the legacy LDSC default only when its column is
    mapped; an explicit threshold requires the column. Returns the effective
    thresholds keyed by canonical field (``None`` for a filter that will not
    run) and ``(flag, message)`` conflicts for explicit thresholds whose column
    is absent.
    """
    thresholds, conflicts = {}, []
    for field, attribute, flag, default, needs, label, kind, hint in _QC_THRESHOLDS:
        requested = getattr(config, attribute)
        if field in mapped_fields:
            thresholds[field] = default if requested is None else requested
            continue
        thresholds[field] = None
        if requested is not None:
            conflicts.append((flag, (
                f"{flag}={requested} requires {needs}, but none was selected. Most likely the "
                f"{kind} column has a name LDSC3 does not recognize or was passed to --ignore. "
                f"Select it with {hint}, or omit {flag} to run without {label} filtering."
            )))
    return thresholds, conflicts


def _old_daner_sample_sizes(file_cnames, source_path):
    """Return the case count, control count and control-frequency column of an old DANER header.

    The counts are the numeric suffixes of the first ``FRQ_A_<cases>`` and
    ``FRQ_U_<controls>`` columns, matched after header cleaning as in format
    detection, so lowercase headers are read too. The last value is a message
    naming every missing column or non-count suffix, or ``None``.
    """
    columns, counts, problems = {}, {}, []
    for prefix, role in (('FRQ_A_', 'cases'), ('FRQ_U_', 'controls')):
        column = next((name for name in file_cnames if clean_header(name).startswith(prefix)), None)
        if column is None:
            problems.append(f"the header has no {prefix}<{role}> column")
            continue
        suffix = column[len(prefix):]
        try:
            count = float(suffix)
        except ValueError:
            count = None
        if count is None or not 0 < count < float('inf'):
            problems.append(f"{suffix!r} in {column} is not a positive count of {role}")
        columns[prefix], counts[prefix] = column, count
    if problems:
        return None, None, None, (
            f"munge-sumstats reads '{source_path}' as old DANER format, which takes the case and control counts "
            f"from columns named FRQ_A_<cases> and FRQ_U_<controls>, but {'; '.join(problems)}. "
            "Most likely the file is not in old DANER format or these columns were renamed. Restore the "
            "count-named frequency columns, or pass another --input-format with --N-cas/--N-con or --N-col."
        )
    return counts['FRQ_A_'], counts['FRQ_U_'], columns['FRQ_U_'], None


def _suggest_signed_sumstat_fix(file_cnames):
    column = likely_signed_sumstat_column(file_cnames)
    if column is not None:
        return f" Try --signed-sumstats {column},0 if that column is the signed effect relative to A1."
    return ""


def _suggest_n_fix(file_cnames):
    if any(clean_header(column) == 'NEFF' for column in file_cnames):
        return " NEFF is not treated as N automatically; pass --N-col NEFF only if that is appropriate for this analysis."
    return ""


def resolve_raw_schema(
    source_path: str,
    raw_config: MungeConfig,
    config: MungeConfig,
    global_config: GlobalConfig,
) -> RawSchema:
    """Read one raw summary-statistics header the way the munging run will.

    The munging run and ``--infer-only`` both call this, so the report and the
    run cannot disagree about the header. It detects the raw format, applies the
    column hints the format implies and comma-separated INFO lists found in the
    first data row (explicit hints take priority), maps raw columns to canonical
    fields, chooses the signed statistic, reads old DANER case and control
    counts and resolves the INFO and MAF thresholds. Only the header and the
    first data row are read; nothing is written.

    Parameters
    ----------
    source_path : str
        Resolved raw summary-statistics path.
    raw_config : MungeConfig
        Supplies the explicit ``column_hints``.
    config : MungeConfig
        Munging options: input format, ignored and INFO-list columns,
        ``signed_sumstats_spec``, ``a1_inc``, sample-size constants and QC
        thresholds.
    global_config : GlobalConfig
        Supplies the SNP identifier mode; allele-aware modes require A1/A2.

    Returns
    -------
    RawSchema
        The resolved reading. Problems that would stop the run are collected in
        ``issues`` instead of being raised, so one report or one error lists
        them all; :meth:`RawSchema.raise_issues` raises them.
    """
    headers = tuple(read_header(source_path))
    clean_set = {clean_header(column) for column in headers}
    detected_format = _detect_format(clean_set, config.sumstats_format)
    inferred_hints = _format_hints(headers, detected_format)
    info_lists, notes = _first_row_info_lists(source_path, headers)
    hints = {**inferred_hints, **raw_config.column_hints}
    config = replace(
        config,
        sumstats_format=detected_format,
        signed_sumstats_spec=hints.get('signed_sumstats', config.signed_sumstats_spec),
        info_list_columns=tuple(dict.fromkeys([*config.info_list_columns, *info_lists]))
        or tuple(column for column in hints.get('info_list', '').split(',') if column),
    )
    columns = {COLUMN_HINT_TARGETS[key]: source for key, source in hints.items()
               if key in COLUMN_HINT_TARGETS and source is not None}
    qc = MungeQC(N=config.N, N_cas=config.N_cas, N_con=config.N_con,
                 N_min=config.N_min, nstudy_min=config.nstudy_min)
    issues = []

    def problem(label, message, note=None, error=LDSCInputError):
        issues.append(RawSchemaIssue(label, message, note, error))

    explicit_n_problem = _explicit_sample_size_problem(columns)
    if explicit_n_problem is not None:
        problem('sample-size options', explicit_n_problem, explicit_n_problem, LDSCUsageError)
    flag_cnames, signed_sumstat_null, signed_spec_problem = _column_hints(columns, config)
    if signed_spec_problem is not None:
        problem('--signed-sumstats', signed_spec_problem, signed_spec_problem, LDSCUsageError)
    ignore_cnames = [clean_header(x) for x in config.ignore_columns]

    # remove LOG_ODDS, BETA, Z, OR from the default list
    if config.signed_sumstats_spec is not None or config.a1_inc:
        aliases = {x: default_cnames[x] for x in default_cnames if default_cnames[x] not in null_values}
    else:
        aliases = default_cnames
    cname_map = get_cname_map(flag_cnames, aliases, ignore_cnames)
    daner_problem = None
    if detected_format == 'daner-old':
        N_cas, N_con, frq_u, daner_problem = _old_daner_sample_sizes(headers, source_path)
        if daner_problem is not None:
            problem('DANER case/control counts', daner_problem, daner_problem)
        else:
            qc = replace(qc, N_cas=N_cas, N_con=N_con)
            cname_map[clean_header(frq_u)] = 'FRQ'

    cname_translation = {x: cname_map[clean_header(x)] for x in headers if
                         clean_header(x) in cname_map}  # note keys not cleaned
    # Sample size is judged only once the explicit sample-size options are valid.
    if explicit_n_problem is None:
        ambiguous_n = _resolve_sample_size_column_strategy(columns, cname_translation, source_path)
        if ambiguous_n is not None:
            problem('conflict: sample size', ambiguous_n, ambiguous_n)
    field_descriptions = {x: describe_cname[cname_translation[x]] for x in cname_translation}
    mapped = set(cname_translation.values())
    if explicit_n_problem is None and daner_problem is None and not (
            qc.N or (qc.N_cas and qc.N_con) or 'N' in mapped or {'N_CAS', 'N_CON'}.issubset(mapped)):
        problem('N', (
            f"munge-sumstats could not determine sample size (N) for '{source_path}'. "
            "Most likely the input has no recognized N, N_CAS/N_CON, or DANER sample-size columns. "
            "Provide --N, provide both --N-cas and --N-con, or include an inferable N column."
            f"{_suggest_n_fix(headers)}"
        ), "NEFF is not treated as N automatically; pass --N-col NEFF only if appropriate." if 'NEFF' in clean_set else None)

    sign_cname = 'SIGNED_SUMSTATS'
    signed_problem = signed_spec_problem is not None
    if config.signed_sumstats_spec is None and not config.a1_inc:
        sign_cnames = [
            x for x in cname_translation if cname_translation[x] in null_values]
        if len(sign_cnames) > 1:
            message = (
                f"munge-sumstats found multiple signed statistic columns in '{source_path}': {sign_cnames}. "
                "Most likely more than one effect column is present. Pass --signed-sumstats <column>,<null_value> "
                "or use --ignore for the extra signed-statistic columns."
            )
            problem('conflict: signed statistic', message, message)
            signed_problem = True
        elif not sign_cnames:
            available = ', '.join(headers)
            accepted = ', '.join(sorted(null_values))
            problem('signed statistic', (
                f"munge-sumstats could not find a signed summary statistic column in '{source_path}'. "
                f"Available columns: {available}. Most likely the effect column has an unrecognized name. "
                f"Expected one of: {accepted}, or pass --signed-sumstats <column>,<null_value>."
                f"{_suggest_signed_sumstat_fix(headers)}"
            ))
            signed_problem = True
        else:
            sign_cname = sign_cnames[0]
            signed_sumstat_null = null_values[cname_translation[sign_cname]]
            cname_translation[sign_cname] = 'SIGNED_SUMSTAT'

    # check that we have all the columns we need
    req_cols = ['SNP', 'P'] if config.a1_inc else ['SNP', 'P', 'SIGNED_SUMSTAT']
    for c in req_cols:
        if c in cname_translation.values() or (c == 'SIGNED_SUMSTAT' and signed_problem):
            continue
        available = ', '.join(headers)
        message = (
            f"munge-sumstats could not map the required column '{c}' from the header of '{source_path}'. "
            f"Available columns: {available}. Most likely the file uses an unrecognized name for it. "
            "Pass the matching column flag or rename the column. Other causes & fixes: "
            "docs/user-guide/troubleshooting.md#munge-sumstats-could-not-map-a-required-column"
        )
        problem('signed statistic' if c == 'SIGNED_SUMSTAT' else c, message, message)

    # check aren't any duplicated column names in mapping
    for field in cname_translation:
        numk = headers.count(field)
        if numk > 1:
            message = (
                f"munge-sumstats found {numk} columns named '{field}' in '{source_path}'. "
                "Most likely the raw header contains duplicate labels. Rename or remove the duplicate column before munging."
            )
            problem(f'duplicate column: {field}', message, message)

    # check multiple different column names don't map to same data field
    targets = list(cname_translation.values())
    for head in dict.fromkeys(targets):
        numc = targets.count(head)
        if numc > 1:
            message = (
                f"munge-sumstats mapped {numc} different input columns to canonical field '{head}' in '{source_path}'. "
                "Most likely an explicit column hint conflicts with an auto-detected alias. "
                "Use --ignore for the extra column or remove the conflicting hint."
            )
            problem(f'conflict: {head}', message, message)

    thresholds, conflicts = resolve_qc_thresholds(config, set(cname_translation.values()))
    for flag, conflict in conflicts:
        problem(f"column required by {flag}",
                f"munge-sumstats cannot apply the requested QC thresholds to '{source_path}'. {conflict}", conflict)
    qc = replace(qc, info_min=thresholds['INFO'], maf_min=thresholds['FRQ'])

    if ('N' in cname_translation.values() or all(x in cname_translation.values() for x in ['N_CAS', 'N_CON']))\
            and 'NSTUDY' in cname_translation.values():
        nstudy = [
            x for x in cname_translation if cname_translation[x] == 'NSTUDY']
        for x in nstudy:
            del cname_translation[x]
    _allele_issues(problem, headers, cname_translation, ignore_cnames, source_path, global_config.snp_identifier)

    return RawSchema(
        source_path=source_path,
        headers=headers,
        detected_format=detected_format,
        inferred_hints=inferred_hints,
        config=config,
        column_map=cname_translation,
        field_descriptions=field_descriptions,
        signed_sumstat_name=sign_cname,
        signed_sumstat_null=signed_sumstat_null,
        qc=qc,
        notes=tuple(notes),
        issues=tuple(issues),
    )


def _allele_issues(problem, headers, cname_translation, ignore_cnames, source_path, mode):
    """Record allele problems: an unpaired allele column, unknown orientation, or no alleles in an allele-aware mode.

    Without alleles, a base mode munges without them unless an orientation-ambiguous
    column shows that the file has alleles whose orientation is unknown.
    """
    alleles = {target: raw for raw, target in cname_translation.items() if target in ('A1', 'A2')}
    ambiguous = [column for column in headers
                 if clean_header(column) in ORIENTATION_AMBIGUOUS_ALLELE_COLUMNS
                 and clean_header(column) not in ignore_cnames and column not in cname_translation]
    base = not is_allele_aware_mode(mode)
    effect_allele = "A1 (the effect allele, the allele the signed statistic refers to)"
    if len(alleles) == 1:
        ((target, lone),) = alleles.items()
        other = 'A2' if target == 'A1' else 'A1'
        message = (
            f"munge-sumstats found only one allele column in '{source_path}': {lone} is read as {target}, but no "
            f"column is read as {other}. {effect_allele} and A2 (the other allele) must come as a pair."
            + (f" Columns that may hold the other allele: {', '.join(ambiguous)}." if ambiguous else "")
            + f" Name the {other} column with --{other.lower()} <column>"
            + (f", or pass --ignore {lone} to munge without alleles, which rg cannot use." if base else ".")
            + f" Other causes & fixes: {_ALLELE_TROUBLESHOOTING}"
        )
        problem('A1/A2 pair', message, message)
    elif not alleles and ambiguous:
        message = (
            f"munge-sumstats cannot tell which allele column of '{source_path}' is {effect_allele}: "
            f"{', '.join(ambiguous)} do not say which allele the signed statistic refers to, and LDSC3 does not guess. "
            "Check in the source of these summary statistics which allele the signed statistic refers to, then pass "
            "--a1 <effect-allele column> --a2 <other column>"
            + (f"; or pass --ignore {','.join(ambiguous)} to munge without alleles, which rg cannot use." if base else ".")
            + f" Other causes & fixes: {_ALLELE_TROUBLESHOOTING}"
        )
        problem('A1/A2 orientation', message, message)
    elif not alleles and not base:
        problem('A1/A2', (
            f"This run is using snp_identifier={mode!r}, which requires A1/A2 allele columns. "
            f"munge-sumstats could not map usable allele columns from '{source_path}'. "
            "Most likely the file has no allele columns or names them in a way LDSC3 does not recognize. "
            f"Pass --a1 <effect-allele column> --a2 <other column>, or rerun with --snp-identifier {identity_base_mode(mode)}."
        ))


def _increasing_allele_columns(column_map):
    """Return ``"<column> as A1"`` or ``"<column> as A2"`` for each INC_ALLELE column read as A1 and DEC_ALLELE column read as A2."""
    return [f"{raw} as {target}" for raw, target in column_map.items()
            if _INCREASING_ALLELE_COLUMNS.get(target) == clean_header(raw)]


def _check_increasing_allele_signs(schema: RawSchema) -> None:
    """Stop when the signed statistic contradicts an INC_ALLELE or DEC_ALLELE column read as A1 or A2.

    Those names say A1 increases the trait, so a signed statistic below its null
    value (a negative BETA, Z or LOG_ODDS, or an OR below 1) means the allele
    columns are not oriented as named. Rows with a missing or non-numeric
    statistic are not counted. Reads the signed-statistic column in chunks.
    """
    named = _increasing_allele_columns(schema.column_map)
    signed = next((raw for raw, target in schema.column_map.items() if target == 'SIGNED_SUMSTAT'), None)
    if not named or signed is None:
        return
    snp = next(raw for raw, target in schema.column_map.items() if target == 'SNP')
    null = schema.signed_sumstat_null
    contradicting, examples = 0, []
    _, compression = get_compression(schema.source_path)
    with pd.read_csv(
        schema.source_path, sep=r"\s+", compression=compression, na_values=[".", "NA"],
        skiprows=count_leading_sumstats_comment_lines(schema.source_path),
        chunksize=schema.config.chunk_size, dtype=str, usecols=[snp, signed],
    ) as reader:
        for chunk in reader:
            values = pd.to_numeric(chunk[signed], errors='coerce')
            below = values < null
            contradicting += int(below.sum())
            examples.extend(f"{name} ({signed} = {value})" for name, value in
                            zip(chunk.loc[below, snp], chunk.loc[below, signed]))
            del examples[5:]
    if contradicting:
        raise LDSCInputError(
            f"munge-sumstats reads {' and '.join(named)} because these names say A1 (the effect allele, the allele "
            f"the signed statistic refers to) increases the trait, but {contradicting} rows of '{schema.source_path}' "
            f"have {signed} below its null value {null:g}, which says A1 decreases it; for example "
            f"{', '.join(examples)}. Most likely these columns do not hold the increasing and decreasing alleles. "
            "Check in the source of these summary statistics which allele the signed statistic refers to, then pass "
            f"--a1 <effect-allele column> --a2 <other column>. Other causes & fixes: {_ALLELE_TROUBLESHOOTING}"
        )


def prepare_munge_input(
    schema: RawSchema,
    global_config: GlobalConfig,
    liftover_request: SumstatsLiftoverRequest,
    restriction_path: str | None,
) -> ResolvedMungeInput:
    """Turn a resolved raw schema into a kernel request: source coordinates and keep-list.

    Raises every header issue of ``schema`` first. Only build inference
    performs a streaming coordinate prepass. No output directory or artifact is
    touched here; the returned request can be executed without a CLI namespace
    or any further format/build decisions.
    """
    schema.raise_issues()
    _check_increasing_allele_signs(schema)
    if schema.detected_format == 'daner-old':
        LOGGER.info(
            f"Inferred that N_cas = {schema.qc.N_cas}, N_con = {schema.qc.N_con} from the FRQ_[A/U] columns."
        )
    LOGGER.info('Interpreting column names as follows:')
    LOGGER.info('\n'.join([x + ':\t' + schema.field_descriptions[x]
                       for x in schema.field_descriptions]) + '\n')

    source_path, config, mode = schema.source_path, schema.config, global_config.snp_identifier
    _, compression = get_compression(source_path)
    metadata_skiprows = count_leading_sumstats_comment_lines(source_path)
    coordinate_metadata = resolve_source_coordinates(source_path, schema.column_map, global_config, compression, metadata_skiprows, restriction_path=restriction_path)
    restriction = prepare_sumstats_restriction(restriction_path, mode, coordinate_metadata['genome_build'])
    info_lists = {clean_header(column) for column in config.info_list_columns}
    return ResolvedMungeInput(
        source_path=source_path,
        column_map=schema.column_map,
        compression=compression,
        metadata_skiprows=metadata_skiprows,
        chunk_size=config.chunk_size,
        info_list_columns=tuple(raw for raw in schema.column_map if clean_header(raw) in info_lists),
        signed_sumstat_null=schema.signed_sumstat_null,
        signed_sumstat_name=schema.signed_sumstat_name,
        a1_inc=config.a1_inc,
        # A statistic oriented to the increasing allele is never negative, so its median is far from null by
        # design; _check_increasing_allele_signs checks every row instead.
        check_signed_median=not config.a1_inc and not _increasing_allele_columns(schema.column_map),
        snp_identifier=mode,
        genome_build=coordinate_metadata['genome_build'],
        coordinate_basis=coordinate_metadata['coordinate_basis'],
        coordinate_metadata=coordinate_metadata,
        restriction=restriction,
        liftover_request=liftover_request,
        qc=schema.qc,
    )


def all_missing_columns(schema: RawSchema) -> list[str]:
    """Return the mapped raw columns that hold no nonmissing value, in header order.

    Reads the mapped columns in chunks of ``schema.config.chunk_size`` and stops
    once every column has a value. A comma-separated INFO list counts as missing
    when none of its study values is numeric.
    """
    empty = set(schema.column_map)
    if not empty:
        return []
    info_lists = {clean_header(column) for column in schema.config.info_list_columns}
    _, compression = get_compression(schema.source_path)
    with pd.read_csv(
        schema.source_path, sep=r"\s+", compression=compression, na_values=[".", "NA"],
        skiprows=count_leading_sumstats_comment_lines(schema.source_path),
        chunksize=schema.config.chunk_size, dtype=str, usecols=list(schema.column_map),
    ) as reader:
        for chunk in reader:
            for column in tuple(empty):
                values = chunk[column]
                if clean_header(column) in info_lists:
                    values = values.map(lambda value: mean_info_list_value(value, column))
                if values.notna().any():
                    empty.remove(column)
            if not empty:
                break
    return [column for column in schema.headers if column in empty]


def resolve_source_coordinates(source_path, cname_translation, global_config, compression, metadata_skiprows, *, restriction_path=None):
    """Resolve source build and basis jointly with coordinate restriction evidence.

    Concrete source builds must agree with declarations and imply one-based coordinates.
    Auto trusts a selected build-specific position column, otherwise samples
    each chromosome before QC. It checks the keep list in the actual filtering
    mode and reports all unresolved/conflicting inputs together.
    """
    from .._genome.build_policy import GenomeBuildCheck

    mode, build = global_config.snp_identifier, global_config.genome_build
    metadata = {'snp_identifier': mode, 'genome_build': build,
                'genome_build_inferred': False,
                'coordinate_basis': '1-based' if identity_mode_family(mode) == 'chr_pos' else None}
    if identity_mode_family(mode) != 'chr_pos':
        return metadata
    check = GenomeBuildCheck(build, context='munge-sumstats', option='--source-genome-build', logger=LOGGER)
    position_columns = [raw for raw, target in cname_translation.items() if target == 'POS']
    chromosome_columns = [raw for raw, target in cname_translation.items() if target == 'CHR']
    if position_columns and chromosome_columns:
        check.table('--raw-sumstats-file', source_path,
                    coordinate_columns=(chromosome_columns[0], position_columns[0]))
    else:
        check.frames('--raw-sumstats-file', source_path, _coordinate_frames(
            source_path, cname_translation, compression=compression, metadata_skiprows=metadata_skiprows))
    check.paths('--sumstats-snps-file', restriction_path, restriction_mode=mode)
    source_inferences = [value for source, value in check.inferences if source == str(source_path)]
    bases = {value.coordinate_basis for value in source_inferences}
    if len(bases) > 1:
        check.unresolved('--raw-sumstats-file', source_path, 'inconsistent coordinate basis across chromosomes; use one-based coordinates')
    build = check.resolve()
    if not source_inferences:
        metadata['genome_build'] = build
        return metadata
    inference = source_inferences[0]
    counts = {key: sum(value.match_counts[key] for value in source_inferences) for key in inference.match_counts}
    inspected = sum(value.inspected_snp_count for value in source_inferences)
    metadata.update(genome_build=build, genome_build_inferred=True, coordinate_basis=inference.coordinate_basis,
        build_inference={'inspected_snp_count': inspected, 'match_counts': counts,
            'match_fractions': {key: count / inspected for key, count in counts.items()},
            'summary_message': '; '.join(value.summary_message for value in source_inferences)})
    if inference.coordinate_basis == '0-based':
        LOGGER.warning('Converting inferred zero-based source positions to one-based coordinates.')
    return metadata


def _coordinate_frames(source_path, cname_translation, *, compression, metadata_skiprows):
    """Stream mapped source coordinates without applying QC or restrictions."""
    raw_chr = [raw for raw, target in cname_translation.items() if target == 'CHR']
    raw_pos = [raw for raw, target in cname_translation.items() if target == 'POS']
    if not raw_chr or not raw_pos:
        raise LDSCInputError(
            f"munge-sumstats could not infer genome_build='auto' for '{source_path}' "
            "because no CHR/POS column pair was mapped. Most likely the coordinate columns use unrecognized names. "
            "Map the existing coordinates with --chr/--pos or supply an input with both columns. "
            "Declaring --source-genome-build does not supply missing chromosome/position values."
        )
    raw_columns = [raw_chr[0], raw_pos[0]]
    reader = pd.read_csv(
        source_path,
        sep=r'\s+',
        header=0,
        compression=compression,
        usecols=raw_columns,
        na_values=['.', 'NA'],
        skiprows=metadata_skiprows,
        chunksize=_BUILD_INFERENCE_CHUNKSIZE,
    )
    with reader:
        frames = (chunk.rename(columns={raw_chr[0]: 'CHR', raw_pos[0]: 'POS'}) for chunk in reader)
        yield from frames
