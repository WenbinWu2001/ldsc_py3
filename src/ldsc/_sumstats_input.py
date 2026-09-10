"""Resolve raw summary-statistics inputs for the munging workflow.

Header aliases, DANER profiles, sample-size selection and bounded source-build
inference belong here. The numerical kernel receives only resolved settings.
"""

import bz2
from dataclasses import replace
import gzip
import logging
import warnings

import pandas as pd

from .column_inference import (
    RAW_SUMSTATS_REQUIRED_OR_OPTIONAL_SPECS,
    RAW_SUMSTATS_SIGNED_STAT_SPECS,
    build_cleaned_alias_lookup,
    clean_header,
)
from .config import GlobalConfig, MungeConfig
from .errors import LDSCInputError, LDSCUsageError
from .genome_build_inference import collect_chr_pos_build_evidence_frame, resolve_chr_pos_table
from ._kernel.liftover import SumstatsLiftoverRequest
from ._kernel.snp_identity import identity_base_mode, identity_mode_family, is_allele_aware_mode
from ._kernel.sumstats_munger import MungeQC, ResolvedMungeInput, prepare_sumstats_restriction

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
    'A1': 'Allele 1; the allele that the signed statistic is relative to, usually the effect/increasing allele.',
    'A2': 'Allele 2; the counterpart allele to A1.',
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


def _decode_header_line(line):
    """Decode one raw header/comment line from plain or compressed input."""
    if isinstance(line, bytes):
        return line.decode('utf-8')
    return line


def count_leading_sumstats_comment_lines(fh):
    """Count leading raw sumstats metadata lines that begin with ``##``."""
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


def _column_hints(columns: dict[str, str], config: MungeConfig):
    '''
    Parse flags that specify how to interpret nonstandard column names.

    flag_cnames is a dict that maps (cleaned) arguments to internal column names
    '''
    flag_cnames = {clean_header(source): target for target, source in columns.items()}
    flag_cnames.update({clean_header(source): 'INFO' for source in config.info_list_columns})

    null_value = None
    if config.signed_sumstats_spec:
        try:
            cname, null_value = config.signed_sumstats_spec.split(',')
            null_value = float(null_value)
            flag_cnames[clean_header(cname)] = 'SIGNED_SUMSTAT'
        except ValueError as exc:
            raise LDSCUsageError(
                f"munge-sumstats could not parse --signed-sumstats={config.signed_sumstats_spec!r}. "
                "Most likely the value is missing '<column>,<null_value>' or the null value is not numeric. "
                "Use a value such as --signed-sumstats BETA,0 or --signed-sumstats OR,1."
            ) from exc

    return [flag_cnames, null_value]


def _validate_explicit_sample_size_column_strategy(columns):
    """Validate mutually exclusive direct-N and case/control column hints."""
    has_direct_n = columns.get("N") is not None
    has_n_cas = columns.get("N_CAS") is not None
    has_n_con = columns.get("N_CON") is not None
    if has_n_cas != has_n_con:
        missing = '--N-con-col' if has_n_cas else '--N-cas-col'
        raise LDSCUsageError(
            f"munge-sumstats --N-cas-col and --N-con-col must be provided together; {missing} is missing. "
            "Pass both case/control column names, or drop the incomplete case/control hint and use --N-col."
        )
    if has_direct_n and has_n_cas:
        raise LDSCUsageError(
            "munge-sumstats --N-col cannot be combined with --N-cas-col and --N-con-col. "
            "Choose one sample-size strategy: a direct per-variant N column, or a paired case/control column strategy."
        )


def _warn_sample_size_column_suppression(message):
    """Emit one visible and logged warning for explicit sample-size precedence."""
    warnings.warn(message, UserWarning, stacklevel=3)
    LOGGER.warning(message)


def _resolve_sample_size_column_strategy(columns, cname_translation, source_path):
    """Apply explicit sample-size precedence and reject ambiguous inference."""
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
        return

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
        return

    sources_by_target = {
        target: [source for source, mapped_target in cname_translation.items() if mapped_target == target]
        for target in ('N', 'N_CAS', 'N_CON')
    }
    if all(sources_by_target[target] for target in ('N', 'N_CAS', 'N_CON')):
        direct = sources_by_target['N'][0]
        n_cas = sources_by_target['N_CAS'][0]
        n_con = sources_by_target['N_CON'][0]
        raise LDSCInputError(
            f"munge-sumstats found multiple sample-size strategies through automatic inference in '{source_path}': "
            f"direct N column '{direct}' and case/control columns '{n_cas}'/'{n_con}'. "
            "LDSC3 will not choose one silently because the values can differ. Choose one explicitly with "
            f"--N-col {direct} or --N-cas-col {n_cas} --N-con-col {n_con}."
        )


def _suggest_allele_fix(file_cnames):
    clean_to_original = {clean_header(column): column for column in file_cnames}
    if 'REF' in clean_to_original and 'ALT' in clean_to_original:
        return f" Try --a1 {clean_to_original['REF']} --a2 {clean_to_original['ALT']} if the signed statistic is relative to REF."
    return ""


def _suggest_signed_sumstat_fix(file_cnames):
    likely = {'EFFECT_SIZE', 'EFFECTSIZE', 'LOGOR', 'LOG_OR', 'BETA_HAT'}
    for column in file_cnames:
        if clean_header(column) in likely:
            return f" Try --signed-sumstats {column},0 if that column is the signed effect relative to A1."
    return ""


def _suggest_n_fix(file_cnames):
    if any(clean_header(column) == 'NEFF' for column in file_cnames):
        return " NEFF is not treated as N automatically; pass --N-col NEFF only if that is appropriate for this analysis."
    return ""


def prepare_munge_input(
    source_path: str,
    raw_config: MungeConfig,
    config: MungeConfig,
    global_config: GlobalConfig,
    liftover_request: SumstatsLiftoverRequest,
    restriction_path: str | None,
) -> ResolvedMungeInput:
    """Resolve raw schema, sample-size choices, source coordinates and keep-list before QC.

    Only build inference performs a streaming coordinate prepass. No output
    directory or artifact is touched here; the returned request can be executed
    without a CLI namespace or any further format/build decisions.
    """
    columns = {COLUMN_HINT_TARGETS[key]: source for key, source in raw_config.column_hints.items()
               if key in COLUMN_HINT_TARGETS and source is not None}
    config = replace(
        config,
        signed_sumstats_spec=raw_config.column_hints.get('signed_sumstats', config.signed_sumstats_spec),
        info_list_columns=config.info_list_columns or tuple(column for column in raw_config.column_hints.get('info_list', '').split(',') if column),
    )
    qc = MungeQC(N=config.N, N_cas=config.N_cas, N_con=config.N_con,
                 info_min=config.info_min, maf_min=config.maf_min,
                 n_min=config.n_min, nstudy_min=config.nstudy_min, keep_maf=config.keep_maf)
    _validate_explicit_sample_size_column_strategy(columns)
    file_cnames = read_header(source_path)  # note keys not cleaned
    flag_cnames, signed_sumstat_null = _column_hints(columns, config)
    if config.ignore_columns:
        ignore_cnames = [clean_header(x) for x in config.ignore_columns]
    else:
        ignore_cnames = []

    # remove LOG_ODDS, BETA, Z, OR from the default list
    if config.signed_sumstats_spec is not None or config.a1_inc:
        mod_default_cnames = {x: default_cnames[
            x] for x in default_cnames if default_cnames[x] not in null_values}
    else:
        mod_default_cnames = default_cnames

    cname_map = get_cname_map(
        flag_cnames, mod_default_cnames, ignore_cnames)
    if config.sumstats_format == 'daner-old':
        frq_u = list(filter(lambda x: x.startswith('FRQ_U_'), file_cnames))[0]
        frq_a = list(filter(lambda x: x.startswith('FRQ_A_'), file_cnames))[0]
        N_cas = float(frq_a[6:])
        N_con = float(frq_u[6:])
        LOGGER.info(
            f"Inferred that N_cas = {N_cas}, N_con = {N_con} from the FRQ_[A/U] columns."
        )
        qc = replace(qc, N_cas=N_cas, N_con=N_con)
        cname_map[frq_u] = 'FRQ'

    if config.sumstats_format == 'daner-new':
        frq_u = list(filter(lambda x: x.startswith('FRQ_U_'), file_cnames))[0]
        cname_map[frq_u] = 'FRQ'
        try:
            dan_cas = clean_header(file_cnames[file_cnames.index('Nca')])
        except ValueError as exc:
            raise LDSCInputError(
                f"munge-sumstats could not find the Nca column required by --daner-new in '{source_path}'. "
                "Most likely the file is not in new-DANER format or uses a different case-count header. "
                "Drop --daner-new or pass the correct case-count column with --N-cas-col."
            ) from exc
        try:
            dan_con = clean_header(file_cnames[file_cnames.index('Nco')])
        except ValueError as exc:
            raise LDSCInputError(
                f"munge-sumstats could not find the Nco column required by --daner-new in '{source_path}'. "
                "Most likely the file is not in new-DANER format or uses a different control-count header. "
                "Drop --daner-new or pass the correct control-count column with --N-con-col."
            ) from exc
        cname_map[dan_cas] = 'N_CAS'
        cname_map[dan_con] = 'N_CON'

    cname_translation = {x: cname_map[clean_header(x)] for x in file_cnames if
                         clean_header(x) in cname_map}  # note keys not cleaned
    _resolve_sample_size_column_strategy(columns, cname_translation, source_path)
    cname_description = {
        x: describe_cname[cname_translation[x]] for x in cname_translation}
    if config.signed_sumstats_spec is None and not config.a1_inc:
        sign_cnames = [
            x for x in cname_translation if cname_translation[x] in null_values]
        if len(sign_cnames) > 1:
            raise LDSCInputError(
                f"munge-sumstats found multiple signed statistic columns in '{source_path}': {sign_cnames}. "
                "Most likely more than one effect column is present. Pass --signed-sumstats <column>,<null_value> "
                "or use --ignore for the extra signed-statistic columns."
            )
        if len(sign_cnames) == 0:
            available = ', '.join(file_cnames)
            accepted = ', '.join(sorted(null_values))
            raise LDSCInputError(
                f"munge-sumstats could not find a signed summary statistic column in '{source_path}'. "
                f"Available columns: {available}. Most likely the effect column has an unrecognized name. "
                f"Expected one of: {accepted}, or pass --signed-sumstats <column>,<null_value>."
                f"{_suggest_signed_sumstat_fix(file_cnames)}"
            )
        sign_cname = sign_cnames[0]
        signed_sumstat_null = null_values[cname_translation[sign_cname]]
        cname_translation[sign_cname] = 'SIGNED_SUMSTAT'
    else:
        sign_cname = 'SIGNED_SUMSTATS'

    # check that we have all the columns we need
    if not config.a1_inc:
        req_cols = ['SNP', 'P', 'SIGNED_SUMSTAT']
    else:
        req_cols = ['SNP', 'P']

    for c in req_cols:
        if c not in cname_translation.values():
            available = ', '.join(file_cnames)
            raise LDSCInputError(
                f"munge-sumstats could not map the required column '{c}' from the header of '{source_path}'. "
                f"Available columns: {available}. Most likely the file uses an unrecognized name for it. "
                "Pass the matching column flag or rename the column. Other causes & fixes: "
                "docs/troubleshooting.md#munge-sumstats-could-not-map-a-required-column"
            )

    # check aren't any duplicated column names in mapping
    for field in cname_translation:
        numk = file_cnames.count(field)
        if numk > 1:
            raise LDSCInputError(
                f"munge-sumstats found {numk} columns named '{field}' in '{source_path}'. "
                "Most likely the raw header contains duplicate labels. Rename or remove the duplicate column before munging."
            )

    # check multiple different column names don't map to same data field
    for head in cname_translation.values():
        numc = list(cname_translation.values()).count(head)
        if numc > 1:
            raise LDSCInputError(
                f"munge-sumstats mapped {numc} different input columns to canonical field '{head}' in '{source_path}'. "
                "Most likely an explicit column hint conflicts with an auto-detected alias. "
                "Use --ignore for the extra column or remove the conflicting hint."
            )

    if (not qc.N) and (not (qc.N_cas and qc.N_con)) and ('N' not in cname_translation.values()) and\
            (any(x not in cname_translation.values() for x in ['N_CAS', 'N_CON'])):
        raise LDSCInputError(
            f"munge-sumstats could not determine sample size (N) for '{source_path}'. "
            "Most likely the input has no recognized N, N_CAS/N_CON, or DANER sample-size columns. "
            "Provide --N, provide both --N-cas and --N-con, or include an inferable N column."
            f"{_suggest_n_fix(file_cnames)}"
        )
    if ('N' in cname_translation.values() or all(x in cname_translation.values() for x in ['N_CAS', 'N_CON']))\
            and 'NSTUDY' in cname_translation.values():
        nstudy = [
            x for x in cname_translation if cname_translation[x] == 'NSTUDY']
        for x in nstudy:
            del cname_translation[x]
    mode = global_config.snp_identifier
    requires_alleles = is_allele_aware_mode(mode)
    if requires_alleles and not all(x in cname_translation.values() for x in ['A1', 'A2']):
        raise LDSCInputError(
            f"This run is using snp_identifier={mode!r}, which requires A1/A2 allele columns. "
            f"munge-sumstats could not map usable allele columns from '{source_path}'. "
            "Most likely the file uses REF/ALT or other unrecognized allele headers. "
            f"Pass --a1/--a2 column hints, or rerun with --snp-identifier {identity_base_mode(mode)}."
            f"{_suggest_allele_fix(file_cnames)}"
        )

    LOGGER.info('Interpreting column names as follows:')
    LOGGER.info('\n'.join([x + ':\t' + cname_description[x]
                       for x in cname_description]) + '\n')

    (openfunc, compression) = get_compression(source_path)
    metadata_skiprows = count_leading_sumstats_comment_lines(source_path)
    coordinate_metadata = _resolve_source_coordinates(source_path, cname_translation, global_config, compression, metadata_skiprows)
    restriction = prepare_sumstats_restriction(restriction_path, mode, coordinate_metadata['genome_build'])
    return ResolvedMungeInput(
        source_path=source_path,
        column_map=cname_translation,
        compression=compression,
        metadata_skiprows=metadata_skiprows,
        chunk_size=config.chunk_size,
        info_list_columns=tuple(raw for raw in cname_translation if clean_header(raw) in {clean_header(column) for column in config.info_list_columns}),
        signed_sumstat_null=signed_sumstat_null,
        signed_sumstat_name=sign_cname,
        a1_inc=config.a1_inc,
        snp_identifier=mode,
        genome_build=coordinate_metadata['genome_build'],
        coordinate_basis=coordinate_metadata['coordinate_basis'],
        coordinate_metadata=coordinate_metadata,
        restriction=restriction,
        liftover_request=liftover_request,
        qc=qc,
    )


def _resolve_source_coordinates(source_path, cname_translation, global_config, compression, metadata_skiprows):
    """Resolve build and basis from raw coordinates before any chunk QC."""
    mode = global_config.snp_identifier
    genome_build = global_config.genome_build
    metadata = {'snp_identifier': mode, 'genome_build': genome_build,
                'genome_build_inferred': False,
                'coordinate_basis': '1-based' if identity_mode_family(mode) == 'chr_pos' else None}
    if identity_mode_family(mode) != 'chr_pos' or genome_build != 'auto':
        return metadata
    coordinate_frame = read_coordinate_evidence(
        source_path,
        cname_translation,
        compression=compression,
        metadata_skiprows=metadata_skiprows,
    )
    try:
        _normalized, inference = resolve_chr_pos_table(
            coordinate_frame,
            context=source_path,
            logger=None,
        )
    except ValueError as exc:
        raise LDSCInputError(
            f"munge-sumstats could not infer genome_build from raw CHR/POS coordinates in "
            f"'{source_path}'. Most likely too few coordinates match a known build. "
            "Pass --source-genome-build hg19 or --source-genome-build hg38 explicitly."
        ) from exc

    metadata = {
        'snp_identifier': mode,
        'genome_build': inference.genome_build,
        'genome_build_inferred': True,
        'coordinate_basis': inference.coordinate_basis,
        'build_inference': {
            'inspected_snp_count': int(inference.inspected_snp_count),
            'match_counts': dict(inference.match_counts),
            'match_fractions': dict(inference.match_fractions),
            'summary_message': inference.summary_message,
        },
    }
    if inference.coordinate_basis == '0-based':
        LOGGER.warning(inference.summary_message)
    else:
        LOGGER.info(inference.summary_message)
    LOGGER.info(f"Resolved genome_build='{inference.genome_build}' before chunk parsing for chr_pos summary statistics.")
    return metadata


def read_coordinate_evidence(source_path, cname_translation, *, compression, metadata_skiprows):
    """Read only raw CHR/POS evidence needed for early genome-build inference."""
    raw_chr = [raw for raw, target in cname_translation.items() if target == 'CHR']
    raw_pos = [raw for raw, target in cname_translation.items() if target == 'POS']
    if not raw_chr or not raw_pos:
        raise LDSCInputError(
            f"munge-sumstats could not infer genome_build='auto' for '{source_path}' "
            "because no CHR/POS column pair was mapped. Most likely the coordinate columns use unrecognized names. "
            "Pass --chr/--pos column hints, or pass --source-genome-build hg19 or hg38 explicitly."
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
        return collect_chr_pos_build_evidence_frame(frames, context=source_path)
