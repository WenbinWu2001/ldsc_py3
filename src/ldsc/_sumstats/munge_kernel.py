# (C) 2014-2019 Brendan Bulik-Sullivan and Hilary Finucane
# Derived from legacy LDSC munge_sumstats.py; GPLv3 (LICENSE).
# Performance enhancement and python3 adaption by Anthony Abrantes 2024.
# Modified 2026-05 for workflow/kernel separation and 2026-09 for resolved
# inputs and QC results. Removed banner attribution restored 2026-09-10.
# See NOTICE for attribution and dated modifications.

"""Summary-statistics QC over resolved raw inputs.

``ResolvedMungeInput`` supplies the raw reader schema, concrete source build,
coordinate basis, numeric thresholds and prepared keep-list. Header inference,
DANER choices and build inference live in ``ldsc._sumstats.munge_input``; CLI parsing
lives in ``ldsc._sumstats.munger`` and artifact writing in ``ldsc._sumstats.outputs``.

``munge_sumstats`` streams plain/gzip/bzip2 input, normalizes coordinates and
applies chunk QC and SNP restriction before concatenation. Whole-table sample
size filtering and sign conversion precede optional liftover and global
identity cleanup. ``MungeResult`` returns the munged table, exclusive drop
counts, coordinate provenance and row-level liftover/identity audit records.
No input configuration is mutated and no artifacts are written here.
"""
import logging
from collections import Counter
from contextlib import contextmanager
from dataclasses import dataclass, field, replace
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import chi2
from .._core.coordinates import (
    complete_coordinate_mask,
    coordinate_missing_mask,
    invalid_chromosome_mask,
    normalize_chr_pos_frame,
    positive_int_position_series,
)
from ..errors import LDSCInputError, LDSCInternalError
from .._runtime.progress import report_phase, advance
from .._core.identifiers import (
    build_packed_chr_pos_series,
    read_global_chr_pos_restriction_key_set,
    read_global_snp_restriction,
    read_snp_restriction_keys,
    restriction_file_has_allele_columns,
)
from .liftover import SumstatsLiftoverRequest, apply_sumstats_liftover
from .._core.snp_identity import (
    RestrictionIdentityKeys,
    _combine_base_and_allele_keys,
    allele_set_series,
    base_key_series,
    canonical_allele_series,
    clean_identity_artifact_table,
    collapse_restriction_identity_keys,
    identity_mode_family,
    is_allele_aware_mode,
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
    _key_indexes: dict = field(default_factory=dict, repr=False, compare=False)

    @property
    def n_rows_removed(self):
        return self.n_rows_before_filter - self.n_rows_kept

    def key_index(self, kind):
        """Return the allele-aware keep-list keys as an Index, built once per run.

        ``"identity"`` gives the matching keys and ``"base"`` their base keys
        without alleles; membership is ``get_indexer(...) >= 0``.
        """
        if kind not in self._key_indexes:
            keys = (self.identity_keys.keys if kind == "identity"
                    else _restriction_base_key_universe(self.identity_keys))
            self._key_indexes[kind] = pd.Index(list(keys), dtype=object)
        return self._key_indexes[kind]


@dataclass(frozen=True)
class MungeQC:
    """Resolved numeric thresholds and constant sample sizes.

    ``info_min``/``maf_min`` are ``None`` when the corresponding filter does
    not run because its column was not selected.
    """

    info_min: float | None = None
    maf_min: float | None = None
    N_min: float | None = None
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
    ``hm3_map`` is the dual-build HM3 map loaded once for HM3 quick liftover
    and the packaged keep-list; it is ``None`` for every other run.
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
    # Whether the signed statistic's median must lie near its null value.
    check_signed_median: bool
    snp_identifier: str
    genome_build: str | None
    coordinate_basis: str | None
    coordinate_metadata: dict
    restriction: _SumstatsRestriction | None
    liftover_request: SumstatsLiftoverRequest
    qc: MungeQC
    hm3_map: pd.DataFrame | None = field(default=None, repr=False, compare=False)


@dataclass(frozen=True)
class MungeResult:
    """Munged rows and explicit accounting/provenance from one kernel run.

    Drop counts are exclusive, ordered stage counts; their sum equals parsed
    input rows minus munged output rows. Detailed coordinate reasons may
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


def read_raw_sumstats(source_path, columns, *, compression, skiprows, chunk_size, dtype=None):
    """Open selected raw summary-statistics columns in chunks with the munging read rules.

    Fields are whitespace-delimited, ``.`` and ``NA`` read as missing, and
    the ``skiprows`` leading ``##`` metadata lines are skipped before the
    header. The munging pass, the source-build check, the
    ``INC_ALLELE``/``DEC_ALLELE`` sign check and the all-missing column scan
    all read the raw file through this reader, so they see the same values.
    Returns a pandas chunk reader to use as a context manager.
    """
    return pd.read_csv(source_path, sep=r'\s+', header=0, compression=compression, usecols=list(columns),
                       na_values=['.', 'NA'], chunksize=chunk_size, skiprows=skiprows, dtype=dtype)


@contextmanager
def _invalid_value_summary():
    """Keep warning totals local to this parse, including a failed partial scan."""
    counts = Counter()
    try:
        yield counts
    finally:
        if counts:
            LOGGER.warning("Invalid-value totals observed: %s. Each count is SNP rows with an invalid value in that column. "
                           "Counts can overlap; these are not exclusive QC drop counts.",
                           ", ".join(f"{column}={count}" for column, count in sorted(counts.items())))


def _record_invalid_values(column, count, accepted, counts):
    """Warn promptly on the first affected chunk and accumulate subsequent rows."""
    if not count:
        return
    if counts is None:
        LOGGER.warning(
            "%d SNP rows had %s outside of %s. The column may be mislabeled. "
            "Check the column mapping and original source values; column-specific QC still applies.",
            count, column, accepted,
        )
        return
    if not counts[column]:
        LOGGER.warning(
            "Invalid %s values outside of %s detected; the column may be mislabeled. "
            "Check the column mapping and original source values. Existing QC rules still apply; "
            "totals follow after parsing.",
            column, accepted,
        )
    counts[column] += int(count)


def filter_pvals(P, *, warning_counts=None):
    '''Remove out-of-bounds P-values'''
    ii = (P > 0) & (P <= 1)
    bad_p = (~ii).sum()
    _record_invalid_values("P", bad_p, "(0,1]", warning_counts)

    return ii


def filter_info(info, qc: MungeQC, *, warning_counts=None):
    '''Remove INFO < qc.info_min and complain about out-of-bounds INFO.'''
    if type(info) is pd.Series:  # one INFO column
        jj = ((info > 2.0) | (info < 0)) & info.notnull()
        ii = info >= qc.info_min
    elif type(info) is pd.DataFrame:  # several INFO columns
        jj = (((info > 2.0) & info.notnull()).any(axis=1) | (
            (info < 0) & info.notnull()).any(axis=1))
        ii = (info.sum(axis=1) >= qc.info_min * (len(info.columns)))
    else:
        raise LDSCInternalError(
            "munge-sumstats stopped unexpectedly while checking imputation-quality (INFO) values: "
            "the parser supplied an unsupported data object instead of a column or table. "
            "Report the command, LDSC3 version, and saved traceback; --log-level DEBUG adds troubleshooting detail."
        )

    bad_info = jj.sum()
    _record_invalid_values("INFO", bad_info, "[0,2]", warning_counts)

    return ii


def filter_frq(frq, qc: MungeQC, *, warning_counts=None):
    '''
    Filter on MAF. Remove MAF < qc.maf_min and out-of-bounds MAF.
    '''
    jj = (frq < 0) | (frq > 1)
    bad_frq = jj.sum()
    _record_invalid_values("FRQ", bad_frq, "[0,1]", warning_counts)

    # Carve-out: sumstats A1 is the EFFECT allele, not the minor allele. The fold
    # below is a local mask for the --maf-min threshold only; never reorient
    # A1/A2 or overwrite FRQ here (FRQ stays freq(A1), may exceed 0.5).
    frq = np.minimum(frq, 1 - frq)
    ii = frq >= qc.maf_min
    return ii & ~jj


def _looks_missing_info_token(value):
    return str(value).strip().upper() in {'', '.', 'NA', 'NAN'}


def mean_info_list_value(value, column):
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
                f"munge-sumstats could not parse INFO-list column '{column}': "
                f"token {token!r} is not numeric or a supported missing value. "
                f"Correct the value or map --info-list to the intended comma-separated per-study INFO column. "
                f"Use --ignore {column} only if omitting this quality field is intended for the analysis."
            ) from exc
    if not numeric:
        return np.nan
    return float(np.mean(numeric))


def _coerce_info_list_columns(dat, columns):
    for raw_col in columns:
        dat[raw_col] = dat[raw_col].map(lambda value, column=raw_col: mean_info_list_value(value, column))
    return dat


def _canonical_alleles(alleles, mode):
    """Return one chunk's ``A1`` or ``A2`` column through the canonical allele cleaner, in every SNP identifier mode.

    Legacy ``munge_sumstats.py`` uppercases both allele columns whenever they
    exist, so the letters a munged file lists, and so the allele comparison in
    ``rg``, do not depend on the mode. ``canonical_allele_series`` reads every
    parsed value as text and keeps missing alleles missing; its stripping and
    empty-to-missing rules change nothing here, because whitespace-delimited
    fields hold no blanks. Allele-aware modes store its pandas string dtype;
    base modes, which do not validate alleles, store Python objects and leave
    a column without letters (numeric allele codes, or only missing values)
    as parsed. A missing base-mode allele stays the float ``NaN`` the reader
    gave it, not ``pd.NA``: the ``rg`` allele comparison of an in-memory table
    needs comparisons that return ``False``.
    """
    if is_allele_aware_mode(mode):
        return canonical_allele_series(alleles)
    if pd.api.types.is_numeric_dtype(alleles):
        return alleles
    cleaned = canonical_allele_series(alleles)
    return cleaned.astype(object).where(cleaned.notna(), np.nan)


@report_phase(LOGGER, "preparation", "summary-statistic rows")
def parse_dat(dat_gen, request: ResolvedMungeInput) -> _ParsedChunks:
    """Parse chunks once, returning their retained rows and exclusive counts."""
    convert_colname = request.column_map
    restriction = replace(request.restriction) if request.restriction is not None else None
    coordinate_counts = _empty_coordinate_drop_counts()
    qc = request.qc
    tot_snps = 0
    dat_list = []
    LOGGER.info(f"Reading summary statistics from {request.source_path} into memory {request.chunk_size} SNPs at a time.")
    drops = dict.fromkeys(('NA', 'coordinates', 'INFO', 'FRQ', 'P', 'sumstats_snps', 'N', 'NSTUDY', 'liftover', 'identity'), 0)
    mode = request.snp_identifier
    with _invalid_value_summary() as warning_counts:
        for dat in dat_gen:
            tot_snps += len(dat)
            advance(len(dat))
            old = len(dat)
            required_raw_cols = [
                raw_col for raw_col in dat.columns
                if convert_colname[raw_col] not in {'INFO', 'CHR', 'POS', 'A1', 'A2'}
            ]
            complete = dat[required_raw_cols].notna().all(axis=1).to_numpy()
            if not complete.all():
                dat = dat.take(np.flatnonzero(complete))
            dat.index = pd.RangeIndex(len(dat))
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
            else:
                dat = _drop_unsupported_chromosome_rows(dat, request, coordinate_counts)
            if len(dat) == 0:
                continue

            ii = np.ones(len(dat), dtype=bool)
            if qc.info_min is not None:
                old = ii.sum()
                ii &= filter_info(dat['INFO'], qc, warning_counts=warning_counts)
                new = ii.sum()
                drops['INFO'] += old - new
                old = new

            if qc.maf_min is not None:
                old = ii.sum()
                ii &= filter_frq(dat['FRQ'], qc, warning_counts=warning_counts)
                new = ii.sum()
                drops['FRQ'] += old - new
                old = new

            old = ii.sum()
            if 'INFO' in dat.columns:
                del dat['INFO']
            ii &= filter_pvals(dat.P, warning_counts=warning_counts)
            new = ii.sum()
            drops['P'] += old - new
            old = new
            for allele in ('A1', 'A2'):
                if allele in dat.columns:
                    dat[allele] = _canonical_alleles(dat[allele], mode)

            if ii.sum() == 0:
                continue

            retained = dat if ii.all() else dat.take(np.flatnonzero(ii))
            retained.index = pd.RangeIndex(len(retained))
            retained = _filter_sumstats_chunk_to_restriction(retained, restriction, request.source_path)
            if len(retained) == 0:
                continue

            dat_list.append(retained)

    dat = pd.concat(dat_list, axis=0, ignore_index=True) if dat_list else pd.DataFrame()
    msg = (
        f"Read {tot_snps} input SNP rows from --raw-sumstats-file.\n"
        f"Removed {drops['NA']} SNPs with missing values.\n"
        + (f"Removed {drops['INFO']} SNPs with INFO < {qc.info_min}.\n" if qc.info_min is not None
           else "INFO filter not applied: no INFO column was selected.\n")
        + (f"Removed {drops['FRQ']} SNPs with MAF < {qc.maf_min}.\n" if qc.maf_min is not None
           else "MAF filter not applied: no allele-frequency column was selected.\n")
        + f"Removed {drops['P']} SNPs with out-of-bounds p-values.\n"
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


def prepare_sumstats_restriction(snps_path: str | None, mode: str, genome_build: str | None, *, hm3_map=None):
    """
    Load the active summary-statistics keep-list before raw chunk parsing.

    The returned object stores either rsID strings or packed uint64 CHR/POS
    keys, depending on the resolved SNP identifier mode. Empty files with no usable
    identifiers fail here, before the raw summary-statistics iterator is consumed.
    ``hm3_map`` is the already loaded dual-build HM3 map that ``snps_path``
    names, passed for HM3 quick liftover so the file is parsed once; its keys
    equal those the file readers give for a valid HM3 map.
    """
    if not snps_path:
        return None

    if hm3_map is not None:
        identifiers, identity_keys = _hm3_restriction_keys(hm3_map, snps_path, mode, genome_build)
    elif identity_mode_family(mode) == 'chr_pos' and (
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


def _hm3_restriction_keys(hm3_map, snps_path, mode, genome_build):
    """Return chr_pos-family keep-list keys from a loaded HM3 map, as the file readers build them.

    Packed ``CHR``/``POS`` keys in ``chr_pos`` mode or for a map without
    alleles; collapsed allele-aware identity keys otherwise.
    """
    positions = hm3_map[f"{genome_build}_POS"]
    if mode == 'chr_pos' or not {'A1', 'A2'}.issubset(hm3_map.columns):
        keys = build_packed_chr_pos_series(hm3_map['CHR'], positions, context=f"restriction file '{snps_path}'")
        return {int(value) for value in keys.to_numpy(dtype=np.uint64)}, None
    identity_keys = collapse_restriction_identity_keys(
        pd.DataFrame({'CHR': hm3_map['CHR'], 'POS': positions, 'A1': hm3_map['A1'], 'A2': hm3_map['A2']}),
        mode,
        context=f"SNP restriction file {Path(snps_path)}",
        has_allele_columns=True,
        logger=LOGGER,
    )
    return identity_keys.keys, identity_keys


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
    retained = dat.take(np.flatnonzero(keep.to_numpy(dtype=bool)))
    retained.index = pd.RangeIndex(len(retained))
    return retained


def _filter_allele_aware_chunk_to_identity_restriction(dat, restriction, source_path):
    """
    Apply an allele-aware restriction without raising on raw rows with bad alleles.

    Invalid raw allele rows must reach final global identity cleanup so the
    dropped-SNP sidecar reports identity-specific reasons.  For membership
    matching, use only the chunk rows that can produce effective identity keys.
    """
    context = f"--sumstats-snps-file filtering for {source_path}"
    # Allele sets and base keys are formed once per chunk; both are per row,
    # so restricting them to the matchable rows equals keying those rows alone.
    allele_set, allele_reasons = allele_set_series(dat, context=context)
    matchable = allele_reasons.isna().to_numpy()
    base_keys = base_key_series(dat, restriction.mode, context=context)
    keys = (base_keys if restriction.identity_keys.match_kind == "base"
            else _combine_base_and_allele_keys(base_keys, allele_set))
    kept_matchable = np.flatnonzero(matchable & (restriction.key_index("identity").get_indexer(keys) >= 0))
    unmatchable = np.flatnonzero(~matchable & (restriction.key_index("base").get_indexer(base_keys) >= 0))
    # Matched rows come first, then rows kept by base key, each in input order.
    retained = dat.take(np.concatenate([kept_matchable, unmatchable]))
    retained.index = pd.RangeIndex(len(retained))

    restriction.n_rows_before_filter += len(dat)
    restriction.n_usable_row_identifiers += int(matchable.sum())
    # Bad allele rows retained by base key are counted as drops only when
    # global identity cleanup actually removes them.
    restriction.n_rows_kept += len(retained)
    return retained


def _restriction_base_key_universe(identity_keys: RestrictionIdentityKeys) -> set[str]:
    """Return base keep-list keys, without alleles, from prepared restriction identity keys."""
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
        "Check that the list and input use compatible identifiers and the same coordinate build when applicable. "
        "Use --no-snp-restriction only if analyzing the broader QC-passing SNP population is intended. "
        f"Keep-list file: {restriction.path}. snp_identifier={restriction.mode}; genome_build={build_label}; "
        f"input rows before filtering={input_rows}; usable row identifiers={usable_rows}; "
        f"keep-list identifiers={len(restriction.identifiers)}.{drop_details} "
        "Other causes & fixes: docs/user-guide/troubleshooting.md#munge-sumstats-no-snps-remain-after-filtering"
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
        min_position=min_position,
    )
    _accumulate_coordinate_drop_counts(counts, report)
    if coordinate_basis == '0-based' and len(normalized) > 0:
        normalized['POS'] = normalized['POS'].astype('int64') + 1
    normalized.index = pd.RangeIndex(len(normalized))  # a new frame, never the caller's
    return normalized


def _drop_unsupported_chromosome_rows(dat, request, counts):
    """Drop rsID-family rows whose complete coordinates carry an unsupported CHR label, as coordinate modes do.

    rsID-family identity does not use coordinates, so rows without a
    complete ``CHR``/``POS`` pair stay. A complete pair whose label is not one
    of 1-22, X, Y, M or MT, such as the contig ``chrUn_gl000220``, cannot be
    written as a coordinate, so the row is dropped and counted as an invalid
    ``CHR`` among the coordinate drops, exactly as coordinate-family modes
    drop it, before any output exists.
    """
    counts['n_input'] += len(dat)
    invalid = pd.Series(False, index=dat.index)
    if {'CHR', 'POS'}.issubset(dat.columns):
        complete = complete_coordinate_mask(dat)
        invalid = invalid_chromosome_mask(dat['CHR'], missing=~complete, context=request.source_path)
    n_invalid = int(invalid.sum())
    if n_invalid:
        if len(counts['examples']) < 5:
            columns = [column for column in ('SNP', 'CHR', 'POS') if column in dat.columns]
            counts['examples'].extend(dat.loc[invalid, columns].head(5 - len(counts['examples'])).to_dict(orient='records'))
        dat = dat.take(np.flatnonzero(~invalid.to_numpy()))
        dat.index = pd.RangeIndex(len(dat))
    counts['n_invalid_chr'] += n_invalid
    counts['n_dropped'] += n_invalid
    counts['n_retained'] += len(dat)
    return dat


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
        LOGGER.debug(
            f"Example rows dropped for invalid or missing CHR/POS in "
            f"{source_path}: {counts['examples']}"
        )


def process_n(dat, qc: MungeQC):
    """Apply the resolved N strategy and whole-table N/NSTUDY filters."""
    if all(i in dat.columns for i in ['N_CAS', 'N_CON']):
        N = dat.N_CAS + dat.N_CON
        P = dat.N_CAS / N
        dat['N'] = N * P / P[N == N.max()].mean()
        del dat['N_CAS'], dat['N_CON']
        # NB no filtering on N done here -- that is done in the next code block

    if 'N' in dat.columns:
        N_min = qc.N_min if qc.N_min else dat.N.quantile(0.9) / 1.5
        old = len(dat)
        keep = (dat.N >= N_min).to_numpy()
        if not keep.all():
            dat = dat.take(np.flatnonzero(keep))
            dat.index = pd.RangeIndex(len(dat))
        new = len(dat)
        LOGGER.info(f"Removed {old - new} SNPs with N < {N_min} ({new} SNPs remain).")

    elif 'NSTUDY' in dat.columns and 'N' not in dat.columns:
        nstudy_min = qc.nstudy_min if qc.nstudy_min else dat.NSTUDY.max()
        old = len(dat)
        keep = (dat.NSTUDY >= nstudy_min).to_numpy()
        if not keep.all():
            dat = dat.take(np.flatnonzero(keep))
            dat.index = pd.RangeIndex(len(dat))
        del dat['NSTUDY']
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
                "munge-sumstats could not derive a sample size (N) after input validation; this is an unexpected internal state. "
                "Report the command, LDSC3 version, and saved traceback. "
                "--log-level DEBUG adds troubleshooting detail."
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
    """Record the chunk coordinate policy and validate optional rsID coordinates.

    Both mode families report the coordinate stage of :func:`parse_dat` from
    ``counts``; rsID-family modes drop only complete pairs with an unsupported
    ``CHR`` label there, so their missing-coordinate count stays 0.
    """
    for column in ('CHR', 'POS'):
        if column not in dat.columns:
            dat[column] = pd.NA
    metadata = {
        **request.coordinate_metadata,
        'coordinate_columns': {target: next((raw for raw, canonical in request.column_map.items() if canonical == target), None) for target in ('CHR', 'POS')},
        'n_rows': counts['n_input'], 'n_retained_after_chr_pos_policy': counts['n_retained'],
        'n_missing_chr_pos': counts['n_missing_chr'] + counts['n_missing_pos'],
        'n_invalid_chr_pos': counts['n_invalid_chr'] + counts['n_invalid_pos'],
        'n_dropped_invalid_chr_pos': counts['n_dropped'], 'coordinate_drop_report': counts,
    }
    if identity_mode_family(request.snp_identifier) != 'chr_pos':
        complete = ~(coordinate_missing_mask(dat['CHR']) | coordinate_missing_mask(dat['POS']))
        if complete.any():
            positive_int_position_series(dat.loc[complete, 'POS'], context=request.source_path, label='POS')
    return metadata


def _require_retained_rows(dat, parsed, stage, advice):
    """Stop the run when ``stage`` left no row, before anything is written: a 0-row artifact would change the analysis."""
    if len(dat):
        return
    removed = ", ".join(f"{reason}={int(count)}" for reason, count in parsed.drop_counts.items() if count) or "none recorded"
    raise LDSCInputError(
        f"munge-sumstats retained 0 of {parsed.n_input_rows} input rows after {stage}. "
        f"Recorded removal counts: {removed}. {advice} "
        "Other causes & fixes: docs/user-guide/troubleshooting.md#munge-sumstats-no-snps-remain-after-filtering"
    )


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
    dat_gen = read_raw_sumstats(request.source_path, cname_translation, compression=request.compression,
                                skiprows=request.metadata_skiprows, chunk_size=request.chunk_size,
                                dtype={c: np.float64 for c in signed_sumstat_cols})

    with dat_gen:
        parsed = parse_dat(dat_gen, request)
    dat = parsed.data
    _require_retained_rows(dat, parsed, "QC and SNP restriction", (
        "Check the mapped input columns, --info-min/--maf-min, the reported QC counts, and the selected SNP restriction. "
        "Correct unintended settings; changing an intended threshold or restriction changes the analysis. "
        "Earlier QC and restriction removals are not individual rows in the dropped-SNP sidecar."
    ))

    # Parsed chunks are concatenated with a fresh RangeIndex; the parsed
    # result is not used again, so the steps below may change it in place.
    median_msg = None
    if request.check_signed_median:
        median_msg = check_median(dat.SIGNED_SUMSTAT, signed_sumstat_null, 0.1, sign_cname)
    coordinate_metadata = _coordinate_metadata(dat, request, parsed.coordinate_counts)
    # filtering on N cannot be done chunkwise
    n_before = len(dat)
    has_n_column = 'N' in dat or {'N_CAS', 'N_CON'}.issubset(dat.columns)
    used_n_rule = 'input_columns' if has_n_column else ('fixed_N' if qc.N else 'fixed_case_control_N')
    dat = process_n(dat, qc)
    parsed.drop_counts['N' if has_n_column else 'NSTUDY'] = n_before - len(dat)
    _require_retained_rows(dat, parsed, "sample-size filtering", (
        "Check --N-min or --nstudy-min against the sample sizes in the input, and the mapped sample-size columns; "
        "changing an intended threshold changes the analysis."
    ))
    dat.P = p_to_z(dat.P, dat.N)
    dat.rename(columns={'P': 'Z'}, inplace=True)
    if not request.a1_inc:
        LOGGER.info(median_msg)
        dat.Z *= (-1) ** (dat.SIGNED_SUMSTAT < signed_sumstat_null)
        del dat['SIGNED_SUMSTAT']
    n_before = len(dat)
    dat, liftover_report, liftover_drop_frame = apply_sumstats_liftover(
        dat, request.liftover_request, source_build=request.genome_build, hm3_map=request.hm3_map, logger=LOGGER,
    )
    parsed.drop_counts['liftover'] = n_before - len(dat)
    cleanup = clean_identity_artifact_table(
        dat, request.snp_identifier, context="munged sumstats",
        stage="post_liftover_identity_cleanup", logger=LOGGER,
    )
    parsed.drop_counts['identity'] = len(dat) - len(cleanup.cleaned)
    dat = cleanup.cleaned
    reasons = ", ".join(f"{reason}={int(count)}" for reason, count in
                        cleanup.dropped["reason"].value_counts(sort=False).items())
    _require_retained_rows(dat, parsed, "SNP identity cleanup", (
        f"No remaining row has a usable SNP identity under snp_identifier='{request.snp_identifier}' "
        f"(identity drop reasons: {reasons or 'none recorded'}). Most likely the allele columns hold codes other than "
        "A/C/G/T or no values (allele-aware identifier modes need both letters), or every SNP is listed more than once. "
        "Check the --a1/--a2 columns, or rerun with a base --snp-identifier if allele-aware identity is not needed."
    ))
    coordinate_metadata['liftover'] = liftover_report
    if liftover_report.get('applied'):
        coordinate_metadata['genome_build'] = liftover_report['target_build']

    LOGGER.info(
        f"Prepared summary statistics for {len(dat)} SNPs ({dat.N.notnull().sum()} with nonmissing sample size N)."
    )

    LOGGER.info('\nMetadata:')
    CHISQ = (dat.Z ** 2)
    mean_chisq = CHISQ.mean()
    LOGGER.info(f"Mean chi^2 = {round(mean_chisq, 3)}")
    if mean_chisq < 1.02:
        LOGGER.warning(
            "Mean chi^2 is %.3f, below 1.02. Munging continues; review the input p-values, "
            "sample-size mapping, and GWAS QC before interpreting regression results. "
            "This value alone does not establish an input error.", mean_chisq
        )

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
