#!/usr/bin/env python
"""Legacy-compatible summary-statistics munging kernel.

Core functionality:
    Parse heterogeneous GWAS summary-statistics files, infer canonical LDSC
    columns, apply quality-control filters, and emit LDSC-ready outputs.

Overview
--------
This module contains the low-level munging implementation reused by the public
``ldsc.sumstats_munger`` workflow wrapper. It stays close to the historical
LDSC behavior so filtering semantics and output formats remain stable while the
rest of the refactored package gains a cleaner public interface. Public CLI
orchestration, output preflight, metadata sidecars, and log-file ownership live
in ``ldsc.sumstats_munger``; this module emits ordinary package logger records
while retaining the legacy-compatible `.sumstats.gz` writer.

The physical raw-input reader accepts plain, gzip-compressed, or bzip2-compressed
whitespace-delimited text. DANER inputs are distinguished by schema flags rather
than file suffix: ``--daner-old`` reads case/control counts from
``FRQ_A_<Ncas>`` and ``FRQ_U_<Ncon>`` headers, while ``--daner-new`` reads
per-SNP case/control counts from exact ``Nca`` and ``Nco`` columns. Optional
HM3 or ``--sumstats-snps-file`` restrictions are prepared once before chunk
parsing and applied inside ``parse_dat`` after canonical columns, coordinate
normalization, and chunk-local QC are complete, so rows outside the keep-list
are not retained for the full-table concatenation step.
"""
import pandas as pd
import numpy as np
#import os
import gzip
import bz2
import argparse
import logging
from dataclasses import dataclass
from scipy.stats import chi2
from .._coordinates import (
    coordinate_missing_mask,
    normalize_chr_pos_frame,
    positive_int_position_series,
)
from ..column_inference import (
    RAW_SUMSTATS_REQUIRED_OR_OPTIONAL_SPECS,
    RAW_SUMSTATS_SIGNED_STAT_SPECS,
    build_cleaned_alias_lookup,
    clean_header,
    normalize_genome_build,
    normalize_snp_identifier_mode,
)
from ..genome_build_inference import resolve_chr_pos_table
from . import regression as sumstats
from .identifiers import build_packed_chr_pos_series, read_global_chr_pos_restriction_key_set, read_global_snp_restriction
from .liftover import SumstatsLiftoverRequest, apply_sumstats_liftover
np.seterr(invalid='ignore')

LOGGER = logging.getLogger("LDSC.sumstats_munger.kernel")


@dataclass
class _SumstatsRestriction:
    """Prepared keep-list and whole-run counters for chunk-stage filtering."""

    path: str
    mode: str
    genome_build: str | None
    identifiers: set
    n_rows_before_filter: int = 0
    n_usable_row_identifiers: int = 0
    n_rows_kept: int = 0

    @property
    def n_rows_removed(self):
        return self.n_rows_before_filter - self.n_rows_kept

try:
    x = pd.DataFrame({'A': [1, 2, 3]})
    x.sort_values(by='A')
except AttributeError:
    raise ImportError('LDSC requires pandas version >= 0.17.0')

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

numeric_cols = ['P', 'N', 'N_CAS', 'N_CON', 'POS', 'Z', 'OR', 'BETA', 'LOG_ODDS', 'INFO', 'FRQ', 'SIGNED_SUMSTAT', 'NSTUDY']


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
    '''
    Figure out which column names to use.

    Priority is
    (1) ignore everything in ignore
    (2) use everything in flags that is not in ignore
    (3) use everything in default that is not in ignore or in flags

    The keys of flag are cleaned. The entries of ignore are not cleaned. The keys of defualt
    are cleaned. But all equality is modulo clean_header().
    
    added list wrapper to flag.keys() for py3
    '''
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


def filter_pvals(P, args):
    '''Remove out-of-bounds P-values'''
    ii = (P > 0) & (P <= 1)
    bad_p = (~ii).sum()
    if bad_p > 0:
        msg = 'WARNING: {N} SNPs had P outside of (0,1]. The P column may be mislabeled.'
        LOGGER.warning(msg.format(N=bad_p))

    return ii


def filter_info(info, args):
    '''Remove INFO < args.info_min (default 0.9) and complain about out-of-bounds INFO.'''
    if type(info) is pd.Series:  # one INFO column
        jj = ((info > 2.0) | (info < 0)) & info.notnull()
        ii = info >= args.info_min
    elif type(info) is pd.DataFrame:  # several INFO columns
        jj = (((info > 2.0) & info.notnull()).any(axis=1) | (
            (info < 0) & info.notnull()).any(axis=1))
        ii = (info.sum(axis=1) >= args.info_min * (len(info.columns)))
    else:
        raise ValueError('Expected pd.DataFrame or pd.Series.')

    bad_info = jj.sum()
    if bad_info > 0:
        msg = 'WARNING: {N} SNPs had INFO outside of [0,1.5]. The INFO column may be mislabeled.'
        LOGGER.warning(msg.format(N=bad_info))

    return ii


def filter_frq(frq, args):
    '''
    Filter on MAF. Remove MAF < args.maf_min and out-of-bounds MAF.
    '''
    jj = (frq < 0) | (frq > 1)
    bad_frq = jj.sum()
    if bad_frq > 0:
        msg = 'WARNING: {N} SNPs had FRQ outside of [0,1]. The FRQ column may be mislabeled.'
        LOGGER.warning(msg.format(N=bad_frq))

    frq = np.minimum(frq, 1 - frq)
    ii = frq > args.maf_min
    return ii & ~jj


def filter_alleles(a):
    '''Remove alleles that do not describe strand-unambiguous SNPs'''
    return a.isin(sumstats.VALID_SNPS)


def filter_signed_sumstats(x, null_value):
    '''Remove signed statistics that equal the declared null value.'''
    ii = x != null_value
    removed = (~ii).sum()
    if removed > 0:
        LOGGER.info('Removed {N} SNPs with null signed summary statistics.'.format(N=removed))
    return ii


def _looks_missing_info_token(value):
    return str(value).strip().upper() in {'', '.', 'NA', 'NAN'}


def _mean_info_list_value(value, column):
    """Return the mean of a comma-separated INFO list, ignoring missing tokens."""
    if pd.isna(value):
        return np.nan
    tokens = [token.strip() for token in str(value).split(',')]
    numeric = []
    for token in tokens:
        if _looks_missing_info_token(token):
            continue
        try:
            numeric.append(float(token))
        except ValueError as exc:
            raise ValueError(
                f"Column {column} contains comma-separated INFO values but token {token!r} is not numeric. "
                f"Use --ignore {column} to skip it, or pass --info-list {column} only for numeric/NA per-study INFO lists."
            ) from exc
    if not numeric:
        return np.nan
    return float(np.mean(numeric))


def _coerce_info_list_columns(dat, convert_colname, args):
    if not getattr(args, 'info_list', None):
        return dat
    info_list = {clean_header(column) for column in args.info_list.split(',') if column.strip()}
    for raw_col in list(dat.columns):
        if clean_header(raw_col) in info_list and convert_colname.get(raw_col) == 'INFO':
            dat[raw_col] = dat[raw_col].map(lambda value, column=raw_col: _mean_info_list_value(value, column))
    return dat


def parse_dat(dat_gen, convert_colname, args, restriction=None):
    '''Parse and filter a sumstats file chunk-wise'''
    tot_snps = 0
    dat_list = []
    msg = 'Reading sumstats from {F} into memory {N} SNPs at a time.'
    LOGGER.info(msg.format(F=args.sumstats, N=int(args.chunksize)))
    drops = {'NA': 0, 'P': 0, 'INFO': 0,
             'FRQ': 0, 'A': 0, 'SNP': 0}
    mode = normalize_snp_identifier_mode(getattr(args, 'snp_identifier', 'chr_pos'))
    if mode == 'chr_pos':
        args._coordinates_finalized_chunkwise = True
        args._coordinate_drop_counts = _empty_coordinate_drop_counts()
    for block_num, dat in enumerate(dat_gen):
        tot_snps += len(dat)
        old = len(dat)
        required_raw_cols = [
            raw_col for raw_col in dat.columns
            if convert_colname[raw_col] not in {'INFO', 'CHR', 'POS'}
        ]
        dat = dat.dropna(axis=0, how="any", subset=required_raw_cols).reset_index(drop=True)
        drops['NA'] += old - len(dat)
        dat = _coerce_info_list_columns(dat, convert_colname, args)
        dat.columns = map(lambda x: convert_colname[x], dat.columns)

        wrong_types = [
            c for c in dat.columns
            if c in numeric_cols and c != 'POS' and not np.issubdtype(dat[c].dtype, np.number)
        ]
        if len(wrong_types) > 0:
            raise ValueError('Columns {} are expected to be numeric'.format(wrong_types))

        if mode == 'chr_pos':
            dat = _normalize_chr_pos_chunk(dat, args)
            if len(dat) == 0:
                continue

        ii = np.array([True for i in range(len(dat))])
        if 'INFO' in dat.columns:
            old = ii.sum()
            ii &= filter_info(dat['INFO'], args)
            new = ii.sum()
            drops['INFO'] += old - new
            old = new

        if 'FRQ' in dat.columns:
            old = ii.sum()
            ii &= filter_frq(dat['FRQ'], args)
            new = ii.sum()
            drops['FRQ'] += old - new
            old = new

        old = ii.sum()
        if args.keep_maf:
            dat.drop(
                [x for x in ['INFO'] if x in dat.columns], inplace=True, axis=1)
        else:
            dat.drop(
                [x for x in ['INFO', 'FRQ'] if x in dat.columns], inplace=True, axis=1)
        ii &= filter_pvals(dat.P, args)
        new = ii.sum()
        drops['P'] += old - new
        old = new
        if not args.no_alleles:
            dat.A1 = dat.A1.str.upper()
            dat.A2 = dat.A2.str.upper()
            ii &= filter_alleles(dat.A1 + dat.A2)
            new = ii.sum()
            drops['A'] += old - new
            old = new

        if ii.sum() == 0:
            continue

        retained = dat[ii].reset_index(drop=True)
        retained = _filter_sumstats_chunk_to_restriction(retained, restriction, args)
        if len(retained) == 0:
            continue

        dat_list.append(retained.reset_index(drop=True))

    dat = pd.concat(dat_list, axis=0).reset_index(drop=True) if dat_list else pd.DataFrame()
    msg = 'Read {N} SNPs from --sumstats file.\n'.format(N=tot_snps)
    msg += 'Removed {N} SNPs with missing values.\n'.format(N=drops['NA'])
    msg += 'Removed {N} SNPs with INFO <= {I}.\n'.format(
        N=drops['INFO'], I=args.info_min)
    msg += 'Removed {N} SNPs with MAF <= {M}.\n'.format(
        N=drops['FRQ'], M=args.maf_min)
    msg += 'Removed {N} SNPs with out-of-bounds p-values.\n'.format(
        N=drops['P'])
    msg += 'Removed {N} variants that were not SNPs or were strand-ambiguous.\n'.format(
        N=drops['A'])
    msg += '{N} SNPs remain.'.format(N=len(dat))
    LOGGER.info(msg)
    _log_coordinate_drop_summary(args)
    if restriction is not None:
        _log_sumstats_restriction_summary(restriction)
        if restriction.n_rows_before_filter > 0 and restriction.n_rows_kept == 0:
            _raise_sumstats_restriction_empty(restriction)
    return dat


def _prepare_sumstats_restriction(args):
    """
    Load the active summary-statistics keep-list before raw chunk parsing.

    The returned object stores either rsID strings or packed uint64 CHR/POS
    keys, depending on ``args.snp_identifier``. Empty files with no usable
    identifiers fail here, before the raw sumstats iterator is consumed.
    """
    snps_path = getattr(args, 'sumstats_snps', None)
    if not snps_path:
        return None

    mode = normalize_snp_identifier_mode(getattr(args, 'snp_identifier', 'chr_pos'))
    coordinate_metadata = getattr(args, '_coordinate_metadata', {})
    genome_build = coordinate_metadata.get('genome_build', getattr(args, 'genome_build', None))
    if mode == 'chr_pos':
        identifiers = read_global_chr_pos_restriction_key_set(
            snps_path,
            genome_build=genome_build,
            logger=LOGGER,
        )
    else:
        identifiers = read_global_snp_restriction(
            snps_path,
            mode,
            genome_build=None,
            logger=LOGGER,
        )

    restriction = _SumstatsRestriction(
        path=str(snps_path),
        mode=mode,
        genome_build=genome_build if mode == 'chr_pos' else None,
        identifiers=identifiers,
    )
    if len(restriction.identifiers) == 0:
        _raise_sumstats_restriction_empty(restriction, input_rows='not inspected', usable_rows=0)
    return restriction


def _filter_sumstats_chunk_to_restriction(dat, restriction, args):
    """
    Apply a prepared keep-list to one canonical, already-QCed chunk.

    ``rsid`` mode matches the canonical ``SNP`` column. ``chr_pos`` mode packs
    normalized ``CHR`` and ``POS`` values into uint64 keys for membership tests,
    avoiding per-row ``"CHR:POS"`` strings in the munger hot path.
    """
    if restriction is None or len(dat) == 0:
        return dat

    old = len(dat)
    if restriction.mode == 'rsid':
        keys = dat['SNP'].astype(str)
        usable_mask = keys.notna()
    else:
        keys = build_packed_chr_pos_series(
            dat['CHR'].reset_index(drop=True),
            dat['POS'].reset_index(drop=True),
            context=f"--sumstats-snps-file filtering for {getattr(args, 'sumstats', 'sumstats')}",
        )
        keys.index = dat.index
        usable_mask = pd.Series(True, index=dat.index)

    keep = keys.isin(restriction.identifiers) & usable_mask
    kept = int(keep.sum())
    restriction.n_rows_before_filter += old
    restriction.n_usable_row_identifiers += int(usable_mask.sum())
    restriction.n_rows_kept += kept
    return dat.loc[keep].reset_index(drop=True)


def _log_sumstats_restriction_summary(restriction):
    build_label = restriction.genome_build if restriction.mode == 'chr_pos' else 'not used'
    LOGGER.info(
        f"Applying --sumstats-snps-file keep-list from {restriction.path} "
        f"using snp_identifier={restriction.mode}, genome_build={build_label}; "
        f"read {len(restriction.identifiers)} keep-list identifiers and found "
        f"{restriction.n_usable_row_identifiers}/{restriction.n_rows_before_filter} usable row identifiers."
    )
    LOGGER.info(
        f"Removed {restriction.n_rows_removed} SNPs not in --sumstats-snps-file "
        f"({restriction.n_rows_kept} SNPs remain; source={restriction.path})."
    )


def _raise_sumstats_restriction_empty(restriction, *, input_rows=None, usable_rows=None):
    build_label = restriction.genome_build if restriction.mode == 'chr_pos' else 'not used'
    input_rows = restriction.n_rows_before_filter if input_rows is None else input_rows
    usable_rows = restriction.n_usable_row_identifiers if usable_rows is None else usable_rows
    raise ValueError(
        "After applying --sumstats-snps-file, no SNPs remain. "
        f"Keep-list file: {restriction.path}. "
        f"snp_identifier={restriction.mode}; genome_build={build_label}; "
        f"input rows before filtering={input_rows}; usable row identifiers={usable_rows}; "
        f"keep-list identifiers={len(restriction.identifiers)}. "
        "Check that the keep-list uses the same identifier mode and genome build as the munged sumstats."
    )


def _empty_coordinate_drop_counts():
    return {
        'n_input': 0,
        'n_retained': 0,
        'n_dropped': 0,
        'n_missing_chr': 0,
        'n_missing_pos': 0,
        'n_invalid_chr': 0,
        'n_invalid_pos': 0,
        'examples': [],
    }


def _normalize_chr_pos_chunk(dat, args):
    if 'CHR' not in dat.columns:
        dat['CHR'] = pd.NA
    if 'POS' not in dat.columns:
        dat['POS'] = pd.NA

    coordinate_basis = getattr(args, '_coordinate_basis', None)
    min_position = 0 if coordinate_basis == '0-based' else 1
    normalized, report = normalize_chr_pos_frame(
        dat,
        context=getattr(args, 'sumstats', 'sumstats'),
        coordinate_policy='drop',
        logger=None,
        example_columns=('SNP', 'CHR', 'POS'),
        min_position=min_position,
    )
    _accumulate_coordinate_drop_counts(args, report)
    if coordinate_basis == '0-based' and len(normalized) > 0:
        normalized['POS'] = normalized['POS'].astype('int64') + 1
    return normalized.reset_index(drop=True)


def _accumulate_coordinate_drop_counts(args, report):
    counts = getattr(args, '_coordinate_drop_counts', None)
    if counts is None:
        counts = _empty_coordinate_drop_counts()
        args._coordinate_drop_counts = counts
    counts['n_input'] += int(report.n_input)
    counts['n_retained'] += int(report.n_retained)
    counts['n_dropped'] += int(report.n_dropped)
    counts['n_missing_chr'] += int(report.n_missing_chr)
    counts['n_missing_pos'] += int(report.n_missing_pos)
    counts['n_invalid_chr'] += int(report.n_invalid_chr)
    counts['n_invalid_pos'] += int(report.n_invalid_pos)
    if report.examples and len(counts['examples']) < 5:
        remaining = 5 - len(counts['examples'])
        counts['examples'].extend(report.examples[:remaining])


def _log_coordinate_drop_summary(args):
    counts = getattr(args, '_coordinate_drop_counts', None)
    if not counts or counts.get('n_dropped', 0) == 0:
        return
    LOGGER.warning(
        f"Dropped {counts['n_dropped']} SNPs with invalid or missing CHR/POS in "
        f"{getattr(args, 'sumstats', 'sumstats')}; {counts['n_retained']} rows remain "
        f"(missing CHR={counts['n_missing_chr']}, missing POS={counts['n_missing_pos']}, "
        f"invalid CHR={counts['n_invalid_chr']}, invalid POS={counts['n_invalid_pos']})."
    )
    if counts['examples']:
        LOGGER.warning(
            f"Example rows dropped for invalid or missing CHR/POS in "
            f"{getattr(args, 'sumstats', 'sumstats')}: {counts['examples']}"
        )


def process_n(dat, args):
    '''Determine sample size from --N* flags or N* columns. Filter out low N SNPs.s'''
    if all(i in dat.columns for i in ['N_CAS', 'N_CON']):
        N = dat.N_CAS + dat.N_CON
        P = dat.N_CAS / N
        dat['N'] = N * P / P[N == N.max()].mean()
        dat.drop(['N_CAS', 'N_CON'], inplace=True, axis=1)
        # NB no filtering on N done here -- that is done in the next code block

    if 'N' in dat.columns:
        n_min = args.n_min if args.n_min else dat.N.quantile(0.9) / 1.5
        old = len(dat)
        dat = dat[dat.N >= n_min].reset_index(drop=True)
        new = len(dat)
        LOGGER.info('Removed {M} SNPs with N < {MIN} ({N} SNPs remain).'.format(
            M=old - new, N=new, MIN=n_min))

    elif 'NSTUDY' in dat.columns and 'N' not in dat.columns:
        nstudy_min = args.nstudy_min if args.nstudy_min else dat.NSTUDY.max()
        old = len(dat)
        dat = dat[dat.NSTUDY >= nstudy_min].drop(
            ['NSTUDY'], axis=1).reset_index(drop=True)
        new = len(dat)
        LOGGER.info('Removed {M} SNPs with NSTUDY < {MIN} ({N} SNPs remain).'.format(
            M=old - new, N=new, MIN=nstudy_min))

    if 'N' not in dat.columns:
        if args.N:
            dat['N'] = args.N
            LOGGER.info('Using N = {N}'.format(N=args.N))
        elif args.N_cas and args.N_con:
            dat['N'] = args.N_cas + args.N_con
            if not args.daner_old:
                msg = 'Using N_cas = {N1}; N_con = {N2}'
                LOGGER.info(msg.format(N1=args.N_cas, N2=args.N_con))
        else:
            raise ValueError('Cannot determine N. This message indicates a bug.\n'
                             'N should have been checked earlier in the program.')

    return dat


def p_to_z(P, N):
    '''Convert P-value and N to standardized beta.'''
    return np.sqrt(chi2.isf(P, 1))


def check_median(x, expected_median, tolerance, name):
    '''Check that median(x) is within tolerance of expected_median.'''
    m = np.median(x)
    if np.abs(m - expected_median) > tolerance:
        msg = 'WARNING: median value of {F} is {V} (should be close to {M}). This column may be mislabeled.'
        raise ValueError(msg.format(F=name, M=expected_median, V=round(m, 2)))
    else:
        msg = 'Median value of {F} was {C}, which seems sensible.'.format(
            C=m, F=name)

    return msg


def parse_flag_cnames(args):
    '''
    Parse flags that specify how to interpret nonstandard column names.

    flag_cnames is a dict that maps (cleaned) arguments to internal column names
    '''
    cname_options = [
        [args.nstudy, 'NSTUDY', '--nstudy'],
        [args.snp, 'SNP', '--snp'],
        [args.chr, 'CHR', '--chr'],
        [args.pos, 'POS', '--pos'],
        [args.N_col, 'N', '--N'],
        [args.N_cas_col, 'N_CAS', '--N-cas-col'],
        [args.N_con_col, 'N_CON', '--N-con-col'],
        [args.a1, 'A1', '--a1'],
        [args.a2, 'A2', '--a2'],
        [args.p, 'P', '--P'],
        [args.frq, 'FRQ', '--nstudy'],
        [args.info, 'INFO', '--info']
    ]
    flag_cnames = {clean_header(x[0]): x[1]
                   for x in cname_options if x[0] is not None}
    if args.info_list:
        try:
            flag_cnames.update(
                {clean_header(x): 'INFO' for x in args.info_list.split(',')})
        except ValueError as exc:
            raise ValueError(
                f"Invalid --info-list value {args.info_list!r}. "
                "Expected a comma-separated list of INFO column names, for example 'INFO,INFO_SCORE'."
            ) from exc

    null_value = None
    if args.signed_sumstats:
        try:
            cname, null_value = args.signed_sumstats.split(',')
            null_value = float(null_value)
            flag_cnames[clean_header(cname)] = 'SIGNED_SUMSTAT'
        except ValueError as exc:
            raise ValueError(
                f"Invalid --signed-sumstats value {args.signed_sumstats!r}. "
                "Expected '<column>,<null_value>', for example 'BETA,0' or 'OR,1'."
            ) from exc

    return [flag_cnames, null_value]


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


def _resolve_auto_genome_build_before_chunks(args, cname_translation, compression, metadata_skiprows):
    """
    Resolve ``genome_build=auto`` for chr_pos sumstats before chunk parsing.

    The prepass reads only the raw coordinate columns selected by column
    inference, then delegates the build and coordinate-basis decision to the
    shared genome-build inference module used by ref-panel-builder.
    """
    mode = normalize_snp_identifier_mode(getattr(args, 'snp_identifier', 'chr_pos'))
    genome_build = normalize_genome_build(getattr(args, 'genome_build', 'hg38'))
    if mode != 'chr_pos' or genome_build != 'auto':
        if mode == 'chr_pos':
            args._coordinate_basis = '1-based'
        return

    coordinate_frame = _read_raw_sumstats_coordinate_frame(
        args,
        cname_translation,
        compression=compression,
        metadata_skiprows=metadata_skiprows,
    )
    try:
        _normalized, inference = resolve_chr_pos_table(
            coordinate_frame,
            context=getattr(args, 'sumstats', 'sumstats'),
            logger=None,
        )
    except ValueError as exc:
        raise ValueError(
            "Unable to infer genome_build from raw sumstats CHR/POS coordinates. "
            f"{exc} Pass --genome-build hg19 or --genome-build hg38 explicitly."
        ) from exc

    args.genome_build = inference.genome_build
    args._coordinate_basis = inference.coordinate_basis
    args._coordinate_metadata = {
        'format': 'ldsc.sumstats.v1',
        'snp_identifier': 'chr_pos',
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
    LOGGER.info(
        f"Resolved genome_build='{inference.genome_build}' before chunk parsing for chr_pos summary statistics."
    )


def _read_raw_sumstats_coordinate_frame(args, cname_translation, *, compression, metadata_skiprows):
    """Read only raw CHR/POS evidence needed for early genome-build inference."""
    raw_chr = [raw for raw, target in cname_translation.items() if target == 'CHR']
    raw_pos = [raw for raw, target in cname_translation.items() if target == 'POS']
    if not raw_chr or not raw_pos:
        raise ValueError(
            "genome_build='auto' requires CHR and POS columns in --sumstats input. "
            "Pass --chr/--pos column hints, or pass --genome-build hg19 or --genome-build hg38 explicitly."
        )
    raw_columns = [raw_chr[0], raw_pos[0]]
    frame = pd.read_csv(
        args.sumstats,
        sep=r'\s+',
        header=0,
        compression=compression,
        usecols=raw_columns,
        na_values=['.', 'NA'],
        skiprows=metadata_skiprows,
    )
    return frame.rename(columns={raw_chr[0]: 'CHR', raw_pos[0]: 'POS'})


def _finalize_coordinate_columns(dat, args):
    """Ensure canonical CHR/POS columns exist and normalize them when requested."""
    if 'CHR' not in dat.columns:
        dat['CHR'] = pd.NA
    if 'POS' not in dat.columns:
        dat['POS'] = pd.NA

    mode = normalize_snp_identifier_mode(getattr(args, 'snp_identifier', 'chr_pos'))
    genome_build = normalize_genome_build(getattr(args, 'genome_build', 'hg38'))
    existing_metadata = dict(getattr(args, '_coordinate_metadata', {}) or {})
    pre_inferred_basis = existing_metadata.get('coordinate_basis')
    source_columns = getattr(args, '_coordinate_source_columns', {})
    metadata = {
        'format': 'ldsc.sumstats.v1',
        'snp_identifier': mode,
        'genome_build': genome_build,
        'genome_build_inferred': bool(existing_metadata.get('genome_build_inferred', False)),
        'coordinate_basis': pre_inferred_basis if pre_inferred_basis in {'0-based', '1-based'} else None,
        'coordinate_columns': {
            'CHR': source_columns.get('CHR'),
            'POS': source_columns.get('POS'),
        },
        'n_rows': int(len(dat)),
        'n_retained_after_chr_pos_policy': int(len(dat)),
        'n_missing_chr_pos': 0,
        'n_invalid_chr_pos': 0,
        'n_dropped_invalid_chr_pos': 0,
    }
    if 'build_inference' in existing_metadata:
        metadata['build_inference'] = existing_metadata['build_inference']

    if mode == 'chr_pos':
        if getattr(args, '_coordinates_finalized_chunkwise', False):
            counts = getattr(args, '_coordinate_drop_counts', _empty_coordinate_drop_counts())
            metadata.update(
                {
                    'n_rows': int(counts['n_input']),
                    'n_retained_after_chr_pos_policy': int(counts['n_retained']),
                    'n_missing_chr_pos': int(counts['n_missing_chr'] + counts['n_missing_pos']),
                    'n_invalid_chr_pos': int(counts['n_invalid_chr'] + counts['n_invalid_pos']),
                    'n_dropped_invalid_chr_pos': int(counts['n_dropped']),
                    'coordinate_drop_report': {
                        'n_input': int(counts['n_input']),
                        'n_retained': int(counts['n_retained']),
                        'n_dropped': int(counts['n_dropped']),
                        'n_missing_chr': int(counts['n_missing_chr']),
                        'n_missing_pos': int(counts['n_missing_pos']),
                        'n_invalid_chr': int(counts['n_invalid_chr']),
                        'n_invalid_pos': int(counts['n_invalid_pos']),
                        'examples': list(counts['examples']),
                    },
                }
            )
            if genome_build == 'auto' and len(dat) > 0:
                normalized, inference = resolve_chr_pos_table(
                    dat,
                    context=getattr(args, 'sumstats', 'sumstats'),
                    logger=LOGGER,
                )
                dat = normalized
                metadata.update(
                    {
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
                )
            elif genome_build == 'auto':
                metadata['genome_build'] = None
                LOGGER.warning('WARNING: Cannot infer genome build because no rows have valid CHR/POS coordinates.')
            elif len(dat) > 0 and genome_build is not None and metadata['coordinate_basis'] is None:
                metadata['coordinate_basis'] = '1-based'
            args._coordinate_metadata = metadata
            return dat
        min_position = 0 if genome_build == 'auto' or pre_inferred_basis == '0-based' else 1
        dat, coordinate_report = normalize_chr_pos_frame(
            dat,
            context=getattr(args, 'sumstats', 'sumstats'),
            coordinate_policy='drop',
            logger=LOGGER,
            example_columns=('SNP', 'CHR', 'POS'),
            min_position=min_position,
        )
        metadata.update(
            {
                'n_retained_after_chr_pos_policy': int(len(dat)),
                'n_missing_chr_pos': int(
                    coordinate_report.n_missing_chr + coordinate_report.n_missing_pos
                ),
                'n_invalid_chr_pos': int(
                    coordinate_report.n_invalid_chr + coordinate_report.n_invalid_pos
                ),
                'n_dropped_invalid_chr_pos': int(coordinate_report.n_dropped),
                'coordinate_drop_report': {
                    'n_input': int(coordinate_report.n_input),
                    'n_retained': int(coordinate_report.n_retained),
                    'n_dropped': int(coordinate_report.n_dropped),
                    'n_missing_chr': int(coordinate_report.n_missing_chr),
                    'n_missing_pos': int(coordinate_report.n_missing_pos),
                    'n_invalid_chr': int(coordinate_report.n_invalid_chr),
                    'n_invalid_pos': int(coordinate_report.n_invalid_pos),
                    'examples': coordinate_report.examples,
                },
            }
        )
        if (
            pre_inferred_basis == '0-based'
            and len(dat) > 0
            and _should_shift_pre_inferred_zero_based(dat, existing_metadata)
        ):
            dat['POS'] = dat['POS'].astype('int64') + 1
        if genome_build == 'auto' and len(dat) > 0:
            normalized, inference = resolve_chr_pos_table(
                dat,
                context=getattr(args, 'sumstats', 'sumstats'),
                logger=LOGGER,
            )
            dat = normalized
            metadata.update(
                {
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
            )
        elif genome_build == 'auto':
            metadata['genome_build'] = None
            LOGGER.warning('WARNING: Cannot infer genome build because no rows have valid CHR/POS coordinates.')
        elif len(dat) > 0 and genome_build is not None and metadata['coordinate_basis'] is None:
            metadata['coordinate_basis'] = '1-based'
    else:
        complete = ~(coordinate_missing_mask(dat['CHR']) | coordinate_missing_mask(dat['POS']))
        if complete.any():
            positive_int_position_series(
                dat.loc[complete, 'POS'],
                context=getattr(args, 'sumstats', 'sumstats'),
                label='POS',
            )

    args._coordinate_metadata = metadata
    return dat


def _should_shift_pre_inferred_zero_based(dat, existing_metadata):
    """Return whether parsed rows still appear to be in raw 0-based coordinates."""
    raw_examples = existing_metadata.get('_raw_pos_examples')
    if not raw_examples:
        return True
    n = min(len(raw_examples), len(dat))
    if n == 0:
        return False
    parsed = pd.to_numeric(dat['POS'].head(n), errors='coerce').tolist()
    raw = pd.to_numeric(pd.Series(raw_examples[:n]), errors='coerce').tolist()
    return parsed == raw

parser = argparse.ArgumentParser()
parser.add_argument('--sumstats', default=None, type=str,
                    help="Input filename.")
parser.add_argument('--N', default=None, type=float,
                    help="Sample size If this option is not set, will try to infer the sample "
                    "size from the input file. If the input file contains a sample size "
                    "column, and this flag is set, the argument to this flag has priority.")
parser.add_argument('--N-cas', default=None, type=float,
                    help="Number of cases. If this option is not set, will try to infer the number "
                    "of cases from the input file. If the input file contains a number of cases "
                    "column, and this flag is set, the argument to this flag has priority.")
parser.add_argument('--N-con', default=None, type=float,
                    help="Number of controls. If this option is not set, will try to infer the number "
                    "of controls from the input file. If the input file contains a number of controls "
                    "column, and this flag is set, the argument to this flag has priority.")
parser.add_argument('--out', default=None, type=str,
                    help="Output filename prefix.")
parser.add_argument('--info-min', default=0.9, type=float,
                    help="Minimum INFO score.")
parser.add_argument('--maf-min', default=0.01, type=float,
                    help="Minimum MAF.")
parser.add_argument('--daner-old', default=False, action='store_true',
                    help="Parse old DANER format with case/control sample sizes encoded in FRQ_A_<Ncas> "
                    "and FRQ_U_<Ncon> column names.")
parser.add_argument('--daner-new', default=False, action='store_true',
                    help="Parse new DANER format with per-SNP case/control sample sizes in exact "
                    "Nca and Nco columns.")
parser.add_argument('--no-alleles', default=False, action="store_true",
                    help="Don't require alleles. Useful if only unsigned summary statistics are available "
                    "and the goal is h2 / partitioned h2 estimation rather than rg estimation.")
parser.add_argument('--n-min', default=None, type=float,
                    help='Minimum N (sample size). Default is (90th percentile N) / 2.')
parser.add_argument('--chunksize', default=1_000_000, type=int,
                    help='Chunksize.')

# optional args to specify column names
parser.add_argument('--snp', default=None, type=str,
                    help='Name of SNP column (if not a name that ldsc understands). NB: case insensitive.')
parser.add_argument('--chr', default=None, type=str,
                    help='Name of chromosome column (if not a name that ldsc understands). NB: case insensitive.')
parser.add_argument('--pos', default=None, type=str,
                    help='Name of base-pair position column (if not a name that ldsc understands). NB: case insensitive.')
parser.add_argument('--N-col', default=None, type=str,
                    help='Name of N column (if not a name that ldsc understands). NB: case insensitive.')
parser.add_argument('--N-cas-col', default=None, type=str,
                    help='Name of N column (if not a name that ldsc understands). NB: case insensitive.')
parser.add_argument('--N-con-col', default=None, type=str,
                    help='Name of N column (if not a name that ldsc understands). NB: case insensitive.')
parser.add_argument('--a1', default=None, type=str,
                    help='Name of A1 column: the allele that the signed statistic is relative to. NB: case insensitive.')
parser.add_argument('--a2', default=None, type=str,
                    help='Name of A2 column: the counterpart allele to A1. NB: case insensitive.')
parser.add_argument('--p', default=None, type=str,
                    help='Name of p-value column (if not a name that ldsc understands). NB: case insensitive.')
parser.add_argument('--frq', default=None, type=str,
                    help='Name of FRQ or MAF column (if not a name that ldsc understands). NB: case insensitive.')
parser.add_argument('--signed-sumstats', default=None, type=str,
                    help='Name of signed sumstat column, comma null value (e.g., Z,0 or OR,1), oriented relative to A1. NB: case insensitive.')
parser.add_argument('--info', default=None, type=str,
                    help='Name of INFO column (if not a name that ldsc understands). NB: case insensitive.')
parser.add_argument('--info-list', default=None, type=str,
                    help='Comma-separated list of numeric/NA per-study INFO columns. Filters on the mean, e.g. IMPINFO=0.852,0.113,NA. NB: case insensitive.')
parser.add_argument('--nstudy', default=None, type=str,
                    help='Name of NSTUDY column (if not a name that ldsc understands). NB: case insensitive.')
parser.add_argument('--nstudy-min', default=None, type=float,
                    help='Minimum # of studies. Default is to remove everything below the max, unless there is an N column,'
                    ' in which case do nothing.')
parser.add_argument('--ignore', default=None, type=str,
                    help='Comma-separated list of column names to ignore.')
parser.add_argument('--a1-inc', default=False, action='store_true',
                    help='A1 is the increasing allele.')
parser.add_argument('--keep-maf', default=False, action='store_true',
                    help='Keep the MAF column (if one exists).')
parser.add_argument('--snp-identifier', default='chr_pos', choices=('rsid', 'chr_pos'),
                    help="SNP identifier mode recorded in munged metadata.")
parser.add_argument('--genome-build', default='hg38', choices=('auto', 'hg19', 'hg37', 'GRCh37', 'hg38', 'GRCh38'),
                    help="Genome build for CHR/POS coordinates, or 'auto' to infer when possible.")


# set p = False for testing in order to prevent printing
def munge_sumstats(args, p=True):
    """Run the historical LDSC munging pipeline.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed legacy-style munging arguments.
    p : bool, optional
        Historical flag controlling whether the processed table is returned for
        programmatic reuse. Default is ``True``.

    Returns
    -------
    pandas.DataFrame
        Munged summary-statistics table when ``p`` is true.
    """
    if args.out is None:
        raise ValueError('The --out flag is required.')
    args._coordinate_metadata = {
        'format': 'ldsc.sumstats.v1',
        'snp_identifier': normalize_snp_identifier_mode(getattr(args, 'snp_identifier', 'chr_pos')),
        'genome_build': normalize_genome_build(getattr(args, 'genome_build', 'hg38')),
        'genome_build_inferred': False,
    }

    if args.sumstats is None:
        raise ValueError('The --sumstats flag is required.')
    if args.daner_old and args.daner_new:
        raise ValueError('--daner-old and --daner-new are not compatible. Use --daner-old for sample '
        'size from FRQ_A/FRQ_U headers, use --daner-new for values from Nca/Nco columns')

    file_cnames = read_header(args.sumstats)  # note keys not cleaned
    flag_cnames, signed_sumstat_null = parse_flag_cnames(args)
    if args.ignore:
        ignore_cnames = [clean_header(x) for x in args.ignore.split(',')]
    else:
        ignore_cnames = []

    # remove LOG_ODDS, BETA, Z, OR from the default list
    if args.signed_sumstats is not None or args.a1_inc:
        mod_default_cnames = {x: default_cnames[
            x] for x in default_cnames if default_cnames[x] not in null_values}
    else:
        mod_default_cnames = default_cnames

    cname_map = get_cname_map(
        flag_cnames, mod_default_cnames, ignore_cnames)
    if args.daner_old:
        frq_u = list(filter(lambda x: x.startswith('FRQ_U_'), file_cnames))[0]
        frq_a = list(filter(lambda x: x.startswith('FRQ_A_'), file_cnames))[0]
        N_cas = float(frq_a[6:])
        N_con = float(frq_u[6:])
        LOGGER.info(
            'Inferred that N_cas = {N1}, N_con = {N2} from the FRQ_[A/U] columns.'.format(N1=N_cas, N2=N_con))
        args.N_cas = N_cas
        args.N_con = N_con
        # drop any N, N_cas, N_con or FRQ columns
        for c in ['N', 'N_CAS', 'N_CON', 'FRQ']:
            for d in [x for x in cname_map if cname_map[x] == 'c']:
                del cname_map[d]
        cname_map[frq_u] = 'FRQ'

    if args.daner_new:
        frq_u = list(filter(lambda x: x.startswith('FRQ_U_'), file_cnames))[0]
        cname_map[frq_u] = 'FRQ'
        try:
            dan_cas = clean_header(file_cnames[file_cnames.index('Nca')])
        except ValueError:
            raise ValueError('Could not find Nca column expected for daner-new format')
        try:
            dan_con = clean_header(file_cnames[file_cnames.index('Nco')])
        except ValueError:
            raise ValueError('Could not find Nco column expected for daner-new format')
        cname_map[dan_cas] = 'N_CAS'
        cname_map[dan_con] = 'N_CON'

    cname_translation = {x: cname_map[clean_header(x)] for x in file_cnames if
                         clean_header(x) in cname_map}  # note keys not cleaned
    args._coordinate_source_columns = {
        target: source for source, target in cname_translation.items()
        if target in {'CHR', 'POS'}
    }
    cname_description = {
        x: describe_cname[cname_translation[x]] for x in cname_translation}
    if args.signed_sumstats is None and not args.a1_inc:
        sign_cnames = [
            x for x in cname_translation if cname_translation[x] in null_values]
        if len(sign_cnames) > 1:
            raise ValueError(
                'Too many signed sumstat columns. Specify which to ignore with the --ignore flag.')
        if len(sign_cnames) == 0:
            available = ', '.join(file_cnames)
            accepted = ', '.join(sorted(null_values))
            raise ValueError(
                f"Could not find a signed summary statistic column in --sumstats input. "
                f"Available columns: {available}. Expected one of: {accepted}, "
                "or pass --signed-sumstats <column>,<null_value>."
                f"{_suggest_signed_sumstat_fix(file_cnames)}"
            )
        sign_cname = sign_cnames[0]
        signed_sumstat_null = null_values[cname_translation[sign_cname]]
        cname_translation[sign_cname] = 'SIGNED_SUMSTAT'    
    else:
        sign_cname = 'SIGNED_SUMSTATS'

    # check that we have all the columns we need
    if not args.a1_inc:
        req_cols = ['SNP', 'P', 'SIGNED_SUMSTAT']
    else:
        req_cols = ['SNP', 'P']

    for c in req_cols:
        if c not in cname_translation.values():
            available = ', '.join(file_cnames)
            raise ValueError(
                f"Could not find {c} column in --sumstats input. "
                f"Available columns: {available}. Use the matching column flag or --ignore to resolve ambiguous headers."
            )

    # check aren't any duplicated column names in mapping
    for field in cname_translation:
        numk = file_cnames.count(field)
        if numk > 1:
            raise ValueError('Found {num} columns named {C}'.format(C=field,num=str(numk)))

    # check multiple different column names don't map to same data field
    for head in cname_translation.values():
        numc = list(cname_translation.values()).count(head)
        if numc > 1:
            raise ValueError('Found {num} different {C} columns'.format(C=head,num=str(numc)))

    if (not args.N) and (not (args.N_cas and args.N_con)) and ('N' not in cname_translation.values()) and\
            (any(x not in cname_translation.values() for x in ['N_CAS', 'N_CON'])):
        raise ValueError(
            'Could not determine N. Provide --N, provide both --N-cas and --N-con, '
            'or include an inferable N column in the input summary statistics.'
            f"{_suggest_n_fix(file_cnames)}"
        )
    if ('N' in cname_translation.values() or all(x in cname_translation.values() for x in ['N_CAS', 'N_CON']))\
            and 'NSTUDY' in cname_translation.values():
        nstudy = [
            x for x in cname_translation if cname_translation[x] == 'NSTUDY']
        for x in nstudy:
            del cname_translation[x]
    if not args.no_alleles and not all(x in cname_translation.values() for x in ['A1', 'A2']):
        raise ValueError(
            'Could not find A1/A2 columns. Provide allele columns with --a1 and --a2, '
            'or pass --no-alleles for h2/partitioned-h2 workflows that do not require allele matching.'
            f"{_suggest_allele_fix(file_cnames)}"
        )

    LOGGER.info('Interpreting column names as follows:')
    LOGGER.info('\n'.join([x + ':\t' + cname_description[x]
                       for x in cname_description]) + '\n')

    (openfunc, compression) = get_compression(args.sumstats)
    metadata_skiprows = count_leading_sumstats_comment_lines(args.sumstats)
    _resolve_auto_genome_build_before_chunks(
        args,
        cname_translation,
        compression,
        metadata_skiprows,
    )
    restriction = _prepare_sumstats_restriction(args)
    args._prepared_sumstats_restriction = restriction
    if normalize_snp_identifier_mode(getattr(args, 'snp_identifier', 'chr_pos')) == 'chr_pos':
        args._coordinates_finalized_chunkwise = True

    # figure out which columns are going to involve sign information, so we can ensure
    # they're read as floats
    signed_sumstat_cols = [k for k,v in cname_translation.items() if v=='SIGNED_SUMSTAT']
    dat_gen = pd.read_csv(args.sumstats, sep=r'\s+', header=0,
            compression=compression, usecols=cname_translation.keys(),
            na_values=['.', 'NA'], iterator=True, chunksize=args.chunksize,
            skiprows=metadata_skiprows,
            dtype={c:np.float64 for c in signed_sumstat_cols})

    dat = parse_dat(dat_gen, cname_translation, args, restriction=restriction)
    if len(dat) == 0:
        raise ValueError('After applying filters, no SNPs remain.')

    if normalize_snp_identifier_mode(getattr(args, 'snp_identifier', 'chr_pos')) == 'rsid':
        old = len(dat)
        dat = dat.drop_duplicates(subset='SNP').reset_index(drop=True)
        new = len(dat)
        LOGGER.info('Removed {M} SNPs with duplicated rs numbers ({N} SNPs remain).'.format(
            M=old - new, N=new))
    else:
        dat = dat.reset_index(drop=True)
        LOGGER.info('Retained duplicated SNP labels because snp_identifier=chr_pos uses CHR/POS as row identity.')
    median_msg = None
    if not args.a1_inc:
        median_msg = check_median(dat.SIGNED_SUMSTAT, signed_sumstat_null, 0.1, sign_cname)
    dat = _finalize_coordinate_columns(dat, args)
    # filtering on N cannot be done chunkwise
    dat = process_n(dat, args)
    dat.P = p_to_z(dat.P, dat.N)
    dat.rename(columns={'P': 'Z'}, inplace=True)
    if not args.a1_inc:
        LOGGER.info(median_msg)
        dat.Z *= (-1) ** (dat.SIGNED_SUMSTAT < signed_sumstat_null)
        dat.drop('SIGNED_SUMSTAT', inplace=True, axis=1)
    dat = _apply_liftover_if_requested(dat, args)

    out_fname = args.out + '.sumstats'
    print_colnames = [
        c for c in ['SNP', 'CHR', 'POS', 'A1', 'A2', 'Z', 'N'] if c in dat.columns]
    if args.keep_maf and 'FRQ' in dat.columns:
        print_colnames.append('FRQ')
    if p:
        msg = 'Writing summary statistics for {M} SNPs ({N} with nonmissing beta) to {F}.'
        LOGGER.info(
            msg.format(M=len(dat), F=out_fname + '.gz', N=dat.N.notnull().sum()))
        dat.to_csv(out_fname + '.gz', sep="\t", index=False,
                   columns=print_colnames, float_format='%.3f', compression = 'gzip')
    else:
        LOGGER.info(
            'Prepared summary statistics for {M} SNPs ({N} with nonmissing beta).'.format(
                M=len(dat),
                N=dat.N.notnull().sum(),
            )
        )

    LOGGER.info('\nMetadata:')
    CHISQ = (dat.Z ** 2)
    mean_chisq = CHISQ.mean()
    LOGGER.info('Mean chi^2 = ' + str(round(mean_chisq, 3)))
    if mean_chisq < 1.02:
        LOGGER.warning("WARNING: mean chi^2 may be too small.")

    LOGGER.info('Lambda GC = ' + str(round(CHISQ.median() / 0.4549, 3)))
    LOGGER.info('Max chi^2 = ' + str(round(CHISQ.max(), 3)))
    LOGGER.info('{N} Genome-wide significant SNPs (some may have been removed by filtering).'.format(N=(CHISQ
                                                                                                    > 29).sum()))
    return dat


def _apply_liftover_if_requested(dat, args):
    """Apply optional summary-statistics liftover after source-build filters."""
    coordinate_metadata = dict(getattr(args, '_coordinate_metadata', {}))
    request = getattr(args, '_liftover_request', None) or SumstatsLiftoverRequest()
    dat, liftover_report, liftover_drop_frame = apply_sumstats_liftover(
        dat,
        request,
        source_build=coordinate_metadata.get('genome_build', getattr(args, 'genome_build', None)),
        snp_identifier=coordinate_metadata.get('snp_identifier', getattr(args, 'snp_identifier', 'chr_pos')),
        logger=LOGGER,
    )
    coordinate_metadata['liftover'] = liftover_report
    coordinate_metadata['liftover_drop_frame'] = liftover_drop_frame
    if liftover_report.get('applied'):
        coordinate_metadata['genome_build'] = liftover_report['target_build']
    args._coordinate_metadata = coordinate_metadata
    return dat
