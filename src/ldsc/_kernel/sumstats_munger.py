# (C) 2014-2019 Brendan Bulik-Sullivan and Hilary Finucane
# Derived from upstream LDSC munge_sumstats.py; GPLv3 (LICENSE).
# Performance enhancement and python3 adaption by Anthony Abrantes 2024.
# Modified 2026-05 for workflow/kernel separation and 2026-09 for resolved
# inputs and QC results. Removed banner attribution restored 2026-09-10.
# See NOTICE for attribution and dated modifications.

"""Summary-statistics QC over resolved raw inputs.

``ResolvedMungeInput`` supplies the raw reader schema, concrete source build,
coordinate basis, numeric thresholds and prepared keep-list. Header inference,
DANER choices and build inference live in ``ldsc._sumstats_input``; CLI parsing
and artifact writing live in ``ldsc.sumstats_munger``.

``munge_sumstats`` streams plain/gzip/bzip2 input, normalizes coordinates and
applies chunk QC and SNP restriction before concatenation. Whole-table sample
size filtering and sign conversion precede optional liftover and global
identity cleanup. ``MungeResult`` returns the curated table, exclusive drop
counts, coordinate provenance and row-level liftover/identity audit records.
No input configuration is mutated and no artifacts are written here.
"""
import logging
from dataclasses import dataclass, replace

import numpy as np
import pandas as pd
from scipy.stats import chi2
from .._coordinates import (
    coordinate_missing_mask,
    normalize_chr_pos_frame,
    positive_int_position_series,
)
from ..errors import LDSCInputError, LDSCInternalError
from . import regression as sumstats
from .identifiers import (
    build_packed_chr_pos_series,
    read_global_chr_pos_restriction_key_set,
    read_global_snp_restriction,
    read_snp_restriction_keys,
    restriction_file_has_allele_columns,
)
from .liftover import SumstatsLiftoverRequest, apply_sumstats_liftover
from .snp_identity import (
    RestrictionIdentityKeys,
    allele_set_series,
    base_key_series,
    clean_identity_artifact_table,
    identity_mode_family,
    is_allele_aware_mode,
    restriction_membership_mask,
)

LOGGER = logging.getLogger("LDSC.sumstats_munger.kernel")


@dataclass
class _SumstatsRestriction:
    """Prepared keep-list and whole-run counters for chunk-stage filtering."""

    path: str
    mode: str
    genome_build: str | None
    identifiers: set
    identity_keys: RestrictionIdentityKeys | None = None
    n_rows_before_filter: int = 0
    n_usable_row_identifiers: int = 0
    n_rows_kept: int = 0

    @property
    def n_rows_removed(self):
        return self.n_rows_before_filter - self.n_rows_kept


@dataclass(frozen=True)
class MungeQC:
    """Resolved numeric thresholds and constant sample sizes."""

    info_min: float = 0.9
    maf_min: float = 0.01
    n_min: float | None = None
    nstudy_min: float | None = None
    N: float | None = None
    N_cas: float | None = None
    N_con: float | None = None


@dataclass(frozen=True)
class ResolvedMungeInput:
    """Raw reader and scientific settings resolved before chunk filtering.

    Coordinate-family inputs carry a concrete source build and raw coordinate
    basis. There are no CLI flags, output paths or hidden result attributes.
    The kernel does not mutate this request, including its prepared keep-list.
    """

    source_path: str
    column_map: dict[str, str]
    compression: str | None
    metadata_skiprows: int
    chunk_size: int
    info_list_columns: tuple[str, ...]
    signed_sumstat_null: float | None
    signed_sumstat_name: str
    a1_inc: bool
    snp_identifier: str
    genome_build: str | None
    coordinate_basis: str | None
    coordinate_metadata: dict
    restriction: _SumstatsRestriction | None
    liftover_request: SumstatsLiftoverRequest
    qc: MungeQC


@dataclass(frozen=True)
class MungeResult:
    """Curated rows and explicit accounting/provenance from one kernel run.

    Drop counts are exclusive, ordered stage counts; their sum equals parsed
    input rows minus curated output rows. Detailed coordinate reasons may
    overlap and are recorded separately in coordinate_metadata. Drop frames
    describe liftover and global identity cleanup, not every QC exclusion.
    """

    data: pd.DataFrame
    n_input_rows: int
    drop_counts: dict[str, int]
    used_n_rule: str
    coordinate_metadata: dict
    liftover_drop_frame: pd.DataFrame
    identity_drop_frame: pd.DataFrame


@dataclass
class _ParsedChunks:
    data: pd.DataFrame
    n_input_rows: int
    drop_counts: dict[str, int]
    coordinate_counts: dict


numeric_cols = ['P', 'N', 'N_CAS', 'N_CON', 'POS', 'Z', 'OR', 'BETA', 'LOG_ODDS', 'INFO', 'FRQ', 'SIGNED_SUMSTAT', 'NSTUDY']


def filter_pvals(P):
    '''Remove out-of-bounds P-values'''
    ii = (P > 0) & (P <= 1)
    bad_p = (~ii).sum()
    if bad_p > 0:
        LOGGER.warning(f"WARNING: {bad_p} SNPs had P outside of (0,1]. The P column may be mislabeled.")

    return ii


def filter_info(info, qc: MungeQC):
    '''Remove INFO < qc.info_min (default 0.9) and complain about out-of-bounds INFO.'''
    if type(info) is pd.Series:  # one INFO column
        jj = ((info > 2.0) | (info < 0)) & info.notnull()
        ii = info >= qc.info_min
    elif type(info) is pd.DataFrame:  # several INFO columns
        jj = (((info > 2.0) & info.notnull()).any(axis=1) | (
            (info < 0) & info.notnull()).any(axis=1))
        ii = (info.sum(axis=1) >= qc.info_min * (len(info.columns)))
    else:
        raise LDSCInternalError(
            "munge-sumstats received an invalid INFO object inside _kernel.sumstats_munger.filter_info(); "
            "expected a pandas Series or DataFrame. Most likely an internal parser contract changed. "
            "Re-run with --debug and report the traceback."
        )

    bad_info = jj.sum()
    if bad_info > 0:
        LOGGER.warning(
            f"WARNING: {bad_info} SNPs had INFO outside of [0,1.5]. The INFO column may be mislabeled."
        )

    return ii


def filter_frq(frq, qc: MungeQC):
    '''
    Filter on MAF. Remove MAF < qc.maf_min and out-of-bounds MAF.
    '''
    jj = (frq < 0) | (frq > 1)
    bad_frq = jj.sum()
    if bad_frq > 0:
        LOGGER.warning(f"WARNING: {bad_frq} SNPs had FRQ outside of [0,1]. The FRQ column may be mislabeled.")

    # Carve-out: sumstats A1 is the EFFECT allele, not the minor allele. The fold
    # below is a local mask for the --maf-min threshold only; never reorient
    # A1/A2 or overwrite FRQ here (FRQ stays freq(A1), may exceed 0.5).
    frq = np.minimum(frq, 1 - frq)
    ii = frq >= qc.maf_min
    return ii & ~jj


def filter_alleles(a):
    '''Remove alleles that do not describe strand-unambiguous SNPs'''
    return a.isin(sumstats.VALID_SNPS)


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
            raise LDSCInputError(
                f"munge-sumstats could not parse INFO-list column '{column}' in _kernel.sumstats_munger._mean_info_list_value(): "
                f"token {token!r} is not numeric or NA. Most likely --info-list was pointed at a non-INFO column. "
                f"Use --ignore {column} to skip it, or pass --info-list {column} only for numeric/NA per-study INFO lists."
            ) from exc
    if not numeric:
        return np.nan
    return float(np.mean(numeric))


def _coerce_info_list_columns(dat, columns):
    for raw_col in columns:
        dat[raw_col] = dat[raw_col].map(lambda value, column=raw_col: _mean_info_list_value(value, column))
    return dat


def parse_dat(dat_gen, request: ResolvedMungeInput) -> _ParsedChunks:
    """Parse chunks once, returning their retained rows and exclusive counts."""
    convert_colname = request.column_map
    restriction = replace(request.restriction) if request.restriction is not None else None
    coordinate_counts = _empty_coordinate_drop_counts()
    qc = request.qc
    tot_snps = 0
    dat_list = []
    LOGGER.info(f"Reading sumstats from {request.source_path} into memory {request.chunk_size} SNPs at a time.")
    drops = dict.fromkeys(('NA', 'coordinates', 'INFO', 'FRQ', 'P', 'sumstats_snps', 'N', 'NSTUDY', 'liftover', 'identity'), 0)
    mode = request.snp_identifier
    for dat in dat_gen:
        tot_snps += len(dat)
        old = len(dat)
        required_raw_cols = [
            raw_col for raw_col in dat.columns
            if convert_colname[raw_col] not in {'INFO', 'CHR', 'POS', 'A1', 'A2'}
        ]
        dat = dat.dropna(axis=0, how="any", subset=required_raw_cols).reset_index(drop=True)
        drops['NA'] += old - len(dat)
        dat = _coerce_info_list_columns(dat, request.info_list_columns)
        dat.columns = map(lambda x: convert_colname[x], dat.columns)

        wrong_types = [
            c for c in dat.columns
            if c in numeric_cols and c != 'POS' and not np.issubdtype(dat[c].dtype, np.number)
        ]
        if len(wrong_types) > 0:
            raise LDSCInputError(
                f"munge-sumstats could not parse numeric column(s) {wrong_types} from '{request.source_path}'. "
                "Most likely one of these columns contains text values or was mapped to the wrong header. "
                "Check the raw header and pass the correct column hints or clean the non-numeric values."
            )

        if identity_mode_family(mode) == 'chr_pos':
            dat = _normalize_chr_pos_chunk(dat, request, coordinate_counts)
            if len(dat) == 0:
                continue

        ii = np.array([True for i in range(len(dat))])
        if 'INFO' in dat.columns:
            old = ii.sum()
            ii &= filter_info(dat['INFO'], qc)
            new = ii.sum()
            drops['INFO'] += old - new
            old = new

        if 'FRQ' in dat.columns:
            old = ii.sum()
            ii &= filter_frq(dat['FRQ'], qc)
            new = ii.sum()
            drops['FRQ'] += old - new
            old = new

        old = ii.sum()
        dat.drop(columns=['INFO'], errors='ignore', inplace=True)
        ii &= filter_pvals(dat.P)
        new = ii.sum()
        drops['P'] += old - new
        old = new
        if is_allele_aware_mode(mode):
            dat.A1 = dat.A1.astype("string").str.upper()
            dat.A2 = dat.A2.astype("string").str.upper()

        if ii.sum() == 0:
            continue

        retained = dat[ii].reset_index(drop=True)
        retained = _filter_sumstats_chunk_to_restriction(retained, restriction, request.source_path)
        if len(retained) == 0:
            continue

        dat_list.append(retained.reset_index(drop=True))

    dat = pd.concat(dat_list, axis=0).reset_index(drop=True) if dat_list else pd.DataFrame()
    msg = (
        f"Read {tot_snps} SNPs from --sumstats file.\n"
        f"Removed {drops['NA']} SNPs with missing values.\n"
        f"Removed {drops['INFO']} SNPs with INFO < {qc.info_min}.\n"
        f"Removed {drops['FRQ']} SNPs with MAF < {qc.maf_min}.\n"
        f"Removed {drops['P']} SNPs with out-of-bounds p-values.\n"
        f"{len(dat)} SNPs remain."
    )
    LOGGER.info(msg)
    _log_coordinate_drop_summary(coordinate_counts, request.source_path)
    if restriction is not None:
        _log_sumstats_restriction_summary(restriction)
        if len(dat) == 0 and restriction.n_rows_before_filter > 0 and restriction.n_rows_kept == 0:
            _raise_sumstats_restriction_empty(restriction)
    drops['coordinates'] = coordinate_counts['n_dropped']
    drops['sumstats_snps'] = restriction.n_rows_removed if restriction is not None else 0
    return _ParsedChunks(dat, tot_snps, {key: int(value) for key, value in drops.items()}, coordinate_counts)


def prepare_sumstats_restriction(snps_path: str | None, mode: str, genome_build: str | None):
    """
    Load the active summary-statistics keep-list before raw chunk parsing.

    The returned object stores either rsID strings or packed uint64 CHR/POS
    keys, depending on the resolved identity mode. Empty files with no usable
    identifiers fail here, before the raw sumstats iterator is consumed.
    """
    if not snps_path:
        return None

    if identity_mode_family(mode) == 'chr_pos' and (
        mode == 'chr_pos' or not restriction_file_has_allele_columns(snps_path, context=str(snps_path))
    ):
        identifiers = read_global_chr_pos_restriction_key_set(
            snps_path,
            genome_build=genome_build,
            logger=LOGGER,
        )
        identity_keys = None
    elif mode in {'rsid_allele_aware', 'chr_pos_allele_aware'}:
        identity_keys = read_snp_restriction_keys(
            snps_path,
            mode,
            genome_build=genome_build if identity_mode_family(mode) == 'chr_pos' else None,
            logger=LOGGER,
        )
        identifiers = identity_keys.keys
    else:
        identifiers = read_global_snp_restriction(
            snps_path,
            mode,
            genome_build=None,
            logger=LOGGER,
        )
        identity_keys = None

    restriction = _SumstatsRestriction(
        path=str(snps_path),
        mode=mode,
        genome_build=genome_build if identity_mode_family(mode) == 'chr_pos' else None,
        identifiers=identifiers,
        identity_keys=identity_keys,
    )
    if len(restriction.identifiers) == 0:
        if restriction.identity_keys is not None:
            _raise_sumstats_restriction_empty(
                restriction,
                input_rows=restriction.identity_keys.n_input_rows,
                usable_rows=restriction.identity_keys.n_retained_keys,
            )
        _raise_sumstats_restriction_empty(restriction, input_rows='not inspected', usable_rows=0)
    return restriction


def _filter_sumstats_chunk_to_restriction(dat, restriction, source_path):
    """
    Apply a prepared keep-list to one canonical, already-QCed chunk.

    ``rsid`` mode matches the canonical ``SNP`` column. ``chr_pos`` mode packs
    normalized ``CHR`` and ``POS`` values into uint64 keys for membership tests,
    avoiding per-row ``"CHR:POS"`` strings in the munger hot path.
    """
    if restriction is None or len(dat) == 0:
        return dat

    if restriction.identity_keys is not None:
        return _filter_allele_aware_chunk_to_identity_restriction(dat, restriction, source_path)
    old = len(dat)
    if identity_mode_family(restriction.mode) == 'rsid':
        keys = dat['SNP'].astype(str)
        usable_mask = keys.notna()
        keep = keys.isin(restriction.identifiers) & usable_mask
    else:
        keys = build_packed_chr_pos_series(
            dat['CHR'].reset_index(drop=True),
            dat['POS'].reset_index(drop=True),
            context=f"--sumstats-snps-file filtering for {source_path}",
        )
        keys.index = dat.index
        usable_mask = pd.Series(True, index=dat.index)
        keep = keys.isin(restriction.identifiers) & usable_mask
    kept = int(keep.sum())
    restriction.n_rows_before_filter += old
    restriction.n_usable_row_identifiers += int(usable_mask.sum())
    restriction.n_rows_kept += kept
    return dat.loc[keep].reset_index(drop=True)


def _filter_allele_aware_chunk_to_identity_restriction(dat, restriction, source_path):
    """
    Apply an allele-aware restriction without raising on raw rows with bad alleles.

    Invalid raw allele rows must reach final global identity cleanup so the
    dropped-SNP sidecar reports identity-specific reasons.  For membership
    matching, use only the chunk rows that can produce effective identity keys.
    """
    original_row_column = '_ldsc_original_chunk_row'
    work = dat.copy()
    work[original_row_column] = range(len(work))
    _, allele_reasons = allele_set_series(
        work,
        context=f"--sumstats-snps-file filtering for {source_path}",
    )
    matchable = work.loc[allele_reasons.isna()].copy()
    if len(matchable) > 0:
        keep_matchable = restriction_membership_mask(
            matchable,
            restriction.identity_keys,
            restriction.mode,
            context=f"--sumstats-snps-file filtering for {source_path}",
        )
        kept_matchable = matchable.loc[keep_matchable].copy()
        matchable_original_rows = set(matchable[original_row_column].tolist())
    else:
        kept_matchable = matchable
        matchable_original_rows = set()

    unmatchable_mask = ~work[original_row_column].isin(matchable_original_rows)
    base_keys = base_key_series(
        work,
        restriction.mode,
        context=f"--sumstats-snps-file filtering for {source_path}",
    )
    unmatchable_base_universe = _restriction_base_key_universe(restriction.identity_keys)
    unmatchable = work.loc[unmatchable_mask & base_keys.isin(unmatchable_base_universe)].copy()
    retained = pd.concat([kept_matchable, unmatchable], axis=0, ignore_index=True)
    retained = retained.drop(columns=[original_row_column], errors='ignore')

    restriction.n_rows_before_filter += len(dat)
    restriction.n_usable_row_identifiers += len(matchable)
    # Bad allele rows retained by base key are counted as drops only when
    # global identity cleanup actually removes them.
    restriction.n_rows_kept += len(retained)
    return retained.reset_index(drop=True)


def _restriction_base_key_universe(identity_keys: RestrictionIdentityKeys) -> set[str]:
    """Return allele-blind keep-list keys from prepared restriction identity keys."""
    if identity_keys.match_kind == "base":
        return {str(key) for key in identity_keys.keys}
    return {str(key).rsplit(":", 2)[0] for key in identity_keys.keys}


def _log_sumstats_restriction_summary(restriction):
    build_label = restriction.genome_build if identity_mode_family(restriction.mode) == 'chr_pos' else 'not used'
    LOGGER.info(
        f"Applying SNP keep-list restriction from {restriction.path} "
        f"using snp_identifier={restriction.mode}, genome_build={build_label}; "
        f"read {len(restriction.identifiers)} keep-list identifiers and found "
        f"{restriction.n_usable_row_identifiers}/{restriction.n_rows_before_filter} usable row identifiers."
    )
    LOGGER.info(
        f"Removed {restriction.n_rows_removed} SNPs not in the keep-list "
        f"({restriction.n_rows_kept} SNPs remain; source={restriction.path})."
    )


def _raise_sumstats_restriction_empty(restriction, *, input_rows=None, usable_rows=None):
    build_label = restriction.genome_build if identity_mode_family(restriction.mode) == 'chr_pos' else 'not used'
    input_rows = restriction.n_rows_before_filter if input_rows is None else input_rows
    usable_rows = restriction.n_usable_row_identifiers if usable_rows is None else usable_rows
    drop_details = ""
    if restriction.identity_keys is not None and not restriction.identity_keys.dropped.empty:
        counts = restriction.identity_keys.dropped["reason"].value_counts(sort=False).to_dict()
        drop_details = f" Restriction rows dropped before matching: {counts}."
    raise LDSCInputError(
        "munge-sumstats no SNPs remain after SNP keep-list restriction. "
        "Most likely the input and keep-list have no matching identifiers or use different genome builds. "
        "Check the list and source build. To process all SNPs subject to QC, use --no-snp-restriction. "
        f"Keep-list file: {restriction.path}. snp_identifier={restriction.mode}; genome_build={build_label}; "
        f"input rows before filtering={input_rows}; usable row identifiers={usable_rows}; "
        f"keep-list identifiers={len(restriction.identifiers)}.{drop_details} "
        "Other causes & fixes: docs/troubleshooting.md#munge-sumstats-no-snps-remain-after-filtering"
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


def _normalize_chr_pos_chunk(dat, request, counts):
    if 'CHR' not in dat.columns:
        dat['CHR'] = pd.NA
    if 'POS' not in dat.columns:
        dat['POS'] = pd.NA

    coordinate_basis = request.coordinate_basis
    min_position = 0 if coordinate_basis == '0-based' else 1
    normalized, report = normalize_chr_pos_frame(
        dat,
        context=request.source_path,
        coordinate_policy='drop',
        logger=None,
        example_columns=('SNP', 'CHR', 'POS'),
        min_position=min_position,
    )
    _accumulate_coordinate_drop_counts(counts, report)
    if coordinate_basis == '0-based' and len(normalized) > 0:
        normalized['POS'] = normalized['POS'].astype('int64') + 1
    return normalized.reset_index(drop=True)


def _accumulate_coordinate_drop_counts(counts, report):
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


def _log_coordinate_drop_summary(counts, source_path):
    if not counts or counts.get('n_dropped', 0) == 0:
        return
    LOGGER.warning(
        f"Dropped {counts['n_dropped']} SNPs with invalid or missing CHR/POS in "
        f"{source_path}; {counts['n_retained']} rows remain "
        f"(missing CHR={counts['n_missing_chr']}, missing POS={counts['n_missing_pos']}, "
        f"invalid CHR={counts['n_invalid_chr']}, invalid POS={counts['n_invalid_pos']})."
    )
    if counts['examples']:
        LOGGER.warning(
            f"Example rows dropped for invalid or missing CHR/POS in "
            f"{source_path}: {counts['examples']}"
        )


def process_n(dat, qc: MungeQC):
    """Apply the resolved N strategy and whole-table N/NSTUDY filters."""
    if all(i in dat.columns for i in ['N_CAS', 'N_CON']):
        N = dat.N_CAS + dat.N_CON
        P = dat.N_CAS / N
        dat['N'] = N * P / P[N == N.max()].mean()
        dat.drop(['N_CAS', 'N_CON'], inplace=True, axis=1)
        # NB no filtering on N done here -- that is done in the next code block

    if 'N' in dat.columns:
        n_min = qc.n_min if qc.n_min else dat.N.quantile(0.9) / 1.5
        old = len(dat)
        dat = dat[dat.N >= n_min].reset_index(drop=True)
        new = len(dat)
        LOGGER.info(f"Removed {old - new} SNPs with N < {n_min} ({new} SNPs remain).")

    elif 'NSTUDY' in dat.columns and 'N' not in dat.columns:
        nstudy_min = qc.nstudy_min if qc.nstudy_min else dat.NSTUDY.max()
        old = len(dat)
        dat = dat[dat.NSTUDY >= nstudy_min].drop(
            ['NSTUDY'], axis=1).reset_index(drop=True)
        new = len(dat)
        LOGGER.info(f"Removed {old - new} SNPs with NSTUDY < {nstudy_min} ({new} SNPs remain).")

    if 'N' not in dat.columns:
        if qc.N:
            dat['N'] = qc.N
            LOGGER.info(f"Using N = {qc.N}")
        elif qc.N_cas and qc.N_con:
            dat['N'] = qc.N_cas + qc.N_con
            LOGGER.info(f"Using N_cas = {qc.N_cas}; N_con = {qc.N_con}")
        else:
            raise LDSCInternalError(
                "munge-sumstats could not derive a sample size (N) and reached a state that should be unreachable. "
                "Most likely an earlier input validation step failed to stop a missing-N file. "
                "Re-run with --debug and report the traceback."
            )

    return dat


def p_to_z(P, N):
    '''Convert P-value and N to standardized beta.'''
    return np.sqrt(chi2.isf(P, 1))


def check_median(x, expected_median, tolerance, name):
    '''Check that median(x) is within tolerance of expected_median.'''
    m = np.median(x)
    if np.abs(m - expected_median) > tolerance:
        raise LDSCInputError(
            f"munge-sumstats rejected signed statistic column '{name}' because its median is {round(m, 2)}, "
            f"but expected a value close to {expected_median}. Most likely the signed statistic column or "
            "its null value is mislabeled. Pass the correct --signed-sumstats <column>,<null_value> option."
        )
    else:
        msg = f"Median value of {name} was {m}, which seems sensible."

    return msg


def _coordinate_metadata(dat, request, counts):
    """Record the chunk coordinate policy and validate optional rsID coordinates."""
    for column in ('CHR', 'POS'):
        if column not in dat.columns:
            dat[column] = pd.NA
    metadata = {
        **request.coordinate_metadata,
        'coordinate_columns': {target: next((raw for raw, canonical in request.column_map.items() if canonical == target), None) for target in ('CHR', 'POS')},
        'n_rows': int(len(dat)), 'n_retained_after_chr_pos_policy': int(len(dat)),
        'n_missing_chr_pos': 0, 'n_invalid_chr_pos': 0, 'n_dropped_invalid_chr_pos': 0,
    }
    if identity_mode_family(request.snp_identifier) == 'chr_pos':
        metadata.update(
            n_rows=counts['n_input'], n_retained_after_chr_pos_policy=counts['n_retained'],
            n_missing_chr_pos=counts['n_missing_chr'] + counts['n_missing_pos'],
            n_invalid_chr_pos=counts['n_invalid_chr'] + counts['n_invalid_pos'],
            n_dropped_invalid_chr_pos=counts['n_dropped'], coordinate_drop_report=counts,
        )
    else:
        complete = ~(coordinate_missing_mask(dat['CHR']) | coordinate_missing_mask(dat['POS']))
        if complete.any():
            positive_int_position_series(dat.loc[complete, 'POS'], context=request.source_path, label='POS')
    return metadata


def munge_sumstats(request: ResolvedMungeInput) -> MungeResult:
    """Execute resolved chunk QC, global N filtering, sign conversion and liftover.

    The input is read once here. Source-build inference and keep-list loading
    have already completed; all data and accounting are returned explicitly.
    """
    cname_translation = request.column_map
    qc = request.qc
    signed_sumstat_null = request.signed_sumstat_null
    sign_cname = request.signed_sumstat_name
    # figure out which columns are going to involve sign information, so we can ensure
    # they're read as floats
    signed_sumstat_cols = [k for k,v in cname_translation.items() if v=='SIGNED_SUMSTAT']
    dat_gen = pd.read_csv(request.source_path, sep=r'\s+', header=0,
            compression=request.compression, usecols=cname_translation.keys(),
            na_values=['.', 'NA'], iterator=True, chunksize=request.chunk_size,
            skiprows=request.metadata_skiprows,
            dtype={c:np.float64 for c in signed_sumstat_cols})

    with dat_gen:
        parsed = parse_dat(dat_gen, request)
    dat = parsed.data
    if len(dat) == 0:
        raise LDSCInputError(
            "munge-sumstats removed every SNP during quality filtering. Most likely --info-min/--maf-min "
            "are too strict for this file. Relax the thresholds and check the dropped-SNP sidecar. "
            "Other causes & fixes: docs/troubleshooting.md#munge-sumstats-no-snps-remain-after-filtering"
        )

    dat = dat.reset_index(drop=True)
    median_msg = None
    if not request.a1_inc:
        median_msg = check_median(dat.SIGNED_SUMSTAT, signed_sumstat_null, 0.1, sign_cname)
    coordinate_metadata = _coordinate_metadata(dat, request, parsed.coordinate_counts)
    # filtering on N cannot be done chunkwise
    n_before = len(dat)
    has_n_column = 'N' in dat or {'N_CAS', 'N_CON'}.issubset(dat.columns)
    used_n_rule = 'input_columns' if has_n_column else ('fixed_N' if qc.N else 'fixed_case_control_N')
    dat = process_n(dat, qc)
    parsed.drop_counts['N' if has_n_column else 'NSTUDY'] = n_before - len(dat)
    dat.P = p_to_z(dat.P, dat.N)
    dat.rename(columns={'P': 'Z'}, inplace=True)
    if not request.a1_inc:
        LOGGER.info(median_msg)
        dat.Z *= (-1) ** (dat.SIGNED_SUMSTAT < signed_sumstat_null)
        dat.drop('SIGNED_SUMSTAT', inplace=True, axis=1)
    n_before = len(dat)
    dat, liftover_report, liftover_drop_frame = apply_sumstats_liftover(
        dat, request.liftover_request, source_build=request.genome_build,
        snp_identifier=request.snp_identifier, logger=LOGGER,
    )
    parsed.drop_counts['liftover'] = n_before - len(dat)
    cleanup = clean_identity_artifact_table(
        dat, request.snp_identifier, context="munged sumstats",
        stage="post_liftover_identity_cleanup", logger=LOGGER,
    )
    parsed.drop_counts['identity'] = len(dat) - len(cleanup.cleaned)
    dat = cleanup.cleaned
    coordinate_metadata['liftover'] = liftover_report
    if liftover_report.get('applied'):
        coordinate_metadata['genome_build'] = liftover_report['target_build']

    LOGGER.info(
        f"Prepared summary statistics for {len(dat)} SNPs ({dat.N.notnull().sum()} with nonmissing beta)."
    )

    LOGGER.info('\nMetadata:')
    CHISQ = (dat.Z ** 2)
    mean_chisq = CHISQ.mean()
    LOGGER.info(f"Mean chi^2 = {round(mean_chisq, 3)}")
    if mean_chisq < 1.02:
        LOGGER.warning("WARNING: mean chi^2 may be too small.")

    LOGGER.info(f"Lambda GC = {round(CHISQ.median() / 0.4549, 3)}")
    LOGGER.info(f"Max chi^2 = {round(CHISQ.max(), 3)}")
    LOGGER.info(f"{(CHISQ > 29).sum()} Genome-wide significant SNPs (some may have been removed by filtering).")
    return MungeResult(
        data=dat,
        n_input_rows=parsed.n_input_rows,
        drop_counts=parsed.drop_counts,
        used_n_rule=used_n_rule,
        coordinate_metadata=coordinate_metadata,
        liftover_drop_frame=liftover_drop_frame,
        identity_drop_frame=cleanup.dropped,
    )
