"""Query unbiased R2 for SNP pairs from R2 matrices.

Public surface for reading an R2 matrix written by ``ldsc build-r2-matrix`` by SNP pair:
the :class:`R2Matrix` handle, the one-shot :func:`query_r2_matrix` wrapper, and the pure
:func:`unbiased_r2_to_pearson_r` converter. See
``docs/reference/query-r2-matrix.md`` for the current input and output contract.
Stored ``SIGN_R`` values use R2 matrix alleles; returned ``sign_r`` values use the
query alleles after orientation harmonization.
"""
from __future__ import annotations

import argparse
import logging
import sys
from dataclasses import dataclass
from typing import Sequence

import numpy as np
import pandas as pd

from .._runtime.cli_help import CLIHelpFormatter
from .._runtime.workflow_logging import args_output_scope, LOG_LEVEL_HELP
from .._core.identifiers import build_snp_id_series
from .r2_format import R2MatrixChromosome, lookup_pairs_in_parquet, r2_matrix_chromosomes, resolve_r2_build_dir
from .._runtime.workflow_logging import log_inputs, log_outputs, materializing_overwrite_guard, workflow_logging
from .outputs import QueryR2MatrixDirectoryWriter, QueryR2MatrixOutputConfig
from .._runtime.paths import ensure_output_directory
from .._core.snp_identity import (
    COMPLEMENT,
    allele_set_series,
    base_key_series,
    identity_base_mode,
    identity_mode_family,
    is_allele_aware_mode,
    normalize_snp_identifier_mode,
)
from .._core.chromosomes import normalize_chromosome
from .._core.columns import (
    A1_COLUMN_SPEC,
    A2_COLUMN_SPEC,
    CHR_COLUMN_SPEC,
    ColumnSpec,
    POS_COLUMN_SPEC,
    SNP_COLUMN_SPEC,
    normalize_column_token,
    resolve_optional_column,
)
from ..errors import LDSCInputError, LDSCUsageError

LOGGER = logging.getLogger("LDSC.r2_matrix_query")


def unbiased_r2_to_pearson_r(r2_adj, n, sign=None):
    """Convert unbiased R2 to a signed Pearson correlation r.

    Inverts the unbiased correction ``r2_adj = r2_raw - (1 - r2_raw) / (n - 2)``
    to recover the biased squared correlation, then takes the square root and
    applies ``sign``. Vectorized: ``r2_adj`` and ``sign`` may be scalars or
    array-likes, and ``n`` may be a scalar or array-like reference sample size.

    Parameters
    ----------
    r2_adj : float or array-like
        Unbiased squared correlation as stored in the R2 matrix.
    n : int or array-like
        Reference sample size used to estimate R2 (``ldsc:n_samples``).
    sign : int, array-like, or None, optional
        Sign of the Pearson correlation (``+1``/``-1``). ``None`` returns the
        non-negative magnitude ``|r|``.

    Returns
    -------
    float or numpy.ndarray
        Signed Pearson r (or ``|r|`` when ``sign`` is ``None``).
    """
    r2_adj = np.asarray(r2_adj, dtype=np.float64)
    n = np.asarray(n, dtype=np.float64)
    r2_raw = (r2_adj * (n - 2.0) + 1.0) / (n - 1.0)
    r2_raw = np.clip(r2_raw, 0.0, 1.0)
    r = np.sqrt(r2_raw)
    if sign is not None:
        r = r * np.asarray(sign, dtype=np.float64)
    return r if r.ndim else r.item()


def query_r2_matrix(
    pairs: pd.DataFrame,
    *,
    r2_matrix_dir=None,
    snp_identifier: str | None = None,
    genome_build: str | None = None,
) -> pd.DataFrame:
    """Open an R2 matrix, query ``pairs`` once, and return the annotated DataFrame.

    Parameters
    ----------
    pairs : pandas.DataFrame
        SNP pairs with endpoint columns suffixed ``_1`` and ``_2``; see
        :meth:`R2Matrix.query_pairs` for identity and allele requirements.
    r2_matrix_dir : str or os.PathLike
        Required build-r2-matrix output directory or its genome-build child.
    snp_identifier : str or None, optional
        The R2 matrix's recorded SNP identifier mode or its base form.
        None uses the recorded mode.
    genome_build : {"hg19", "hg38"} or None, optional
        Select the coordinate build when the R2 matrix directory contains both.
        A concrete value must agree with normalized input declarations.
        Otherwise coordinate-panel shard
        declarations must agree. Pair endpoints use the selected R2 matrix build;
        their coordinate build is not independently inferred.

    Returns
    -------
    pandas.DataFrame
        Copy of the input with ``r2``, ``sign_r``, ``r``, and ``status``.
        Signed results use query allele orientation and are missing in base
        modes. See :meth:`R2Matrix.query_pairs` for missing-pair semantics.
    """
    matrix = R2Matrix.open(
        r2_matrix_dir,
        snp_identifier=snp_identifier,
        genome_build=genome_build,
    )
    return matrix.query_pairs(pairs)


@dataclass
class _ChromState:
    """Cached per-chromosome state for one open R2 matrix chromosome."""

    key_index: pd.Index          # identity key -> row position (unique)
    a1: np.ndarray               # R2 matrix A1 (object), aligned to sidecar rows
    a2: np.ndarray               # R2 matrix A2 (object), aligned to sidecar rows
    parquet_file: object         # pyarrow.parquet.ParquetFile
    n_snps: int
    r2_scale: float


class R2Matrix:
    """Handle over an R2 matrix written by ``build-r2-matrix``, for pair queries.

    Open once with :meth:`open`, then call :meth:`query_pairs` many times.
    :meth:`open` validates every chromosome's parquet footer up front. Each
    chromosome's sidecar (validated and bound to its parquet), identity-key map,
    and parquet handle are loaded lazily on first touch, then cached.
    """

    def __init__(self, chromosomes: dict[str, R2MatrixChromosome], snp_identifier: str) -> None:
        """Store the opened ``{chrom: R2MatrixChromosome}`` map and the active mode."""
        self._chromosomes = chromosomes
        self._mode = snp_identifier
        self._state: dict[str, _ChromState] = {}

    @classmethod
    def open(
        cls,
        r2_matrix_dir=None,
        *,
        snp_identifier: str | None = None,
        genome_build: str | None = None,
    ) -> "R2Matrix":
        """Open an R2 matrix directory and lazily cache chromosomes on first use.

        ``r2_matrix_dir`` is required and may name the output root or one build
        child. ``genome_build`` selects hg19/hg38 when needed; an omitted
        ``snp_identifier`` uses the mode recorded in R2 matrix metadata. Each
        chromosome requires a Parquet satisfying the ``build-r2-matrix``
        metadata contract and its matching SNP sidecar; every chromosome's
        files and Parquet metadata are validated here, and all defects are
        reported in one error. Each sidecar is validated against its parquet's
        identity hash and row count when its chromosome is first queried. An
        explicit ``snp_identifier`` must be the R2 matrix's recorded mode or
        that mode's base form. Returns a reusable :class:`R2Matrix` handle.
        """
        if r2_matrix_dir is None:
            raise LDSCUsageError(
                "Opening an R2 matrix requires an R2 matrix directory. "
                "Set --r2-matrix-dir (r2_matrix_dir in Python) to a build-r2-matrix output directory "
                "containing R2 Parquet files."
            )

        build_dir = resolve_r2_build_dir(r2_matrix_dir, genome_build)
        members = r2_matrix_chromosomes(build_dir)
        if not members:
            raise LDSCInputError(
                f"Could not open R2 matrix '{r2_matrix_dir}': no chr*_r2.parquet files were found. "
                "The directory may not contain build-r2-matrix outputs, or the selected genome-build "
                "subdirectory may be wrong. Set --r2-matrix-dir to the matching R2 matrix output directory "
                "and check --genome-build."
            )

        chromosomes = {chromosome.chrom: chromosome for chromosome in
                       R2MatrixChromosome.open_all(build_dir, members, command="query-r2-matrix")}
        contracts = [chromosome.contract for chromosome in chromosomes.values()]
        panel_modes = sorted({contract.snp_identifier for contract in contracts})
        sample_sizes = sorted({contract.n_samples for contract in contracts})
        mixed = []
        if len(panel_modes) > 1:
            mixed.append(f"different SNP identifier modes ({', '.join(panel_modes)})")
        if len(sample_sizes) > 1:
            mixed.append(f"different sample sizes (ldsc:n_samples {', '.join(map(str, sample_sizes))})")
        if mixed:
            raise LDSCInputError(
                f"Could not open R2 matrix '{r2_matrix_dir}': its chromosome files record "
                f"{' and '.join(mixed)}. Most likely files from different build-r2-matrix runs were "
                "combined. Use one complete R2 matrix directory or regenerate the R2 matrix."
            )
        panel_mode = panel_modes[0]
        mode = panel_mode if snp_identifier is None else normalize_snp_identifier_mode(snp_identifier)
        if mode not in {panel_mode, identity_base_mode(panel_mode)}:
            allowed = sorted({panel_mode, identity_base_mode(panel_mode)})
            raise LDSCInputError(
                f"query-r2-matrix cannot use --snp-identifier {mode} with R2 matrix '{r2_matrix_dir}', which was "
                f"built in {panel_mode} mode. An R2 matrix answers only its own mode or that mode's base "
                f"form ({' or '.join(allowed)}). Use one of those modes, or an R2 matrix built in {mode} mode."
            )
        if mode.startswith('chr_pos'):
            from .._genome.build_policy import GenomeBuildCheck
            check = GenomeBuildCheck(genome_build, context='query-r2-matrix', logger=LOGGER)
            for chromosome in chromosomes.values():
                check.table('--r2-matrix-dir', chromosome.parquet_path)
            check.resolve()
        return cls(chromosomes, mode)

    @property
    def chromosomes(self) -> list[str]:
        """Chromosomes available in this R2 matrix, in genomic order."""
        return list(self._chromosomes)

    @property
    def snp_identifier(self) -> str:
        """Active SNP identifier mode used to resolve query SNPs."""
        return self._mode

    @property
    def n_samples(self) -> int:
        """Reference sample size, which :meth:`open` requires to agree across chromosome files."""
        return next(iter(self._chromosomes.values())).contract.n_samples

    def _chrom_state(self, chrom: str) -> _ChromState:
        """Load and cache per-chromosome state on first touch."""
        chrom = normalize_chromosome(chrom)
        if chrom in self._state:
            return self._state[chrom]
        if chrom not in self._chromosomes:
            raise LDSCInputError(
                f"The R2 matrix has no chromosome {chrom}; available chromosomes: {self.chromosomes}. "
                "The query requests a chromosome absent from this R2 matrix. "
                "Check the selected R2 matrix directory and genome build, or use an R2 matrix covering "
                "the intended query chromosome."
            )
        chromosome = self._chromosomes[chrom]
        import pyarrow.parquet as pq

        sidecar = chromosome.sidecar()
        pf = pq.ParquetFile(chromosome.parquet_path)
        key = build_snp_id_series(sidecar, self._mode)
        key_index = pd.Index(key)
        if key_index.has_duplicates:
            raise LDSCInputError(
                f"R2Matrix cannot key chromosome {chrom} in mode '{self._mode}': the sidecar "
                "has duplicate SNP identity keys, so a query SNP could resolve ambiguously. "
                "Most likely the mode is too coarse for the R2 matrix (e.g. duplicate positions "
                "in chr_pos). Use an allele-aware mode."
            )
        state = _ChromState(
            key_index=key_index,
            a1=sidecar["A1"].to_numpy(dtype=object),
            a2=sidecar["A2"].to_numpy(dtype=object),
            parquet_file=pf,
            n_snps=chromosome.contract.n_snps,
            r2_scale=chromosome.contract.r2_scale,
        )
        self._state[chrom] = state
        return state

    def _resolve_endpoint(self, endpoint: pd.DataFrame):
        """Resolve one endpoint table to (chrom, row, panel_a1, panel_a2) arrays.

        ``endpoint`` carries canonical columns ``CHR/POS/SNP/A1/A2`` (subset by
        mode). Returns parallel arrays the length of ``endpoint``: chromosome
        label (object, ``None`` if unresolved), row index (int64, ``-1`` if
        unresolved), and R2 matrix A1/A2 (object, ``None`` if unresolved).
        """
        n = len(endpoint)
        chrom_out = np.full(n, None, dtype=object)
        idx_out = np.full(n, -1, dtype=np.int64)
        pa1_out = np.full(n, None, dtype=object)
        pa2_out = np.full(n, None, dtype=object)
        family = identity_mode_family(self._mode)
        keys = _query_keys(endpoint, self._mode)

        if family == "chr_pos" or "CHR" in endpoint.columns:
            chroms = endpoint["CHR"].map(lambda v: normalize_chromosome(v)).to_numpy(dtype=object)
            for chrom in {c for c in chroms if c in self._chromosomes}:
                rows = np.where(chroms == chrom)[0]
                state = self._chrom_state(chrom)
                positions = state.key_index.get_indexer(keys[rows])
                found = positions >= 0
                sel = rows[found]
                idx_out[sel] = positions[found]
                chrom_out[sel] = chrom
                pa1_out[sel] = state.a1[positions[found]]
                pa2_out[sel] = state.a2[positions[found]]
        else:
            # rsid family without CHR column: search every chromosome.
            for chrom in self._chromosomes:
                state = self._chrom_state(chrom)
                positions = state.key_index.get_indexer(keys)
                found = (positions >= 0) & (idx_out < 0)
                idx_out[found] = positions[found]
                chrom_out[found] = chrom
                pa1_out[found] = state.a1[positions[found]]
                pa2_out[found] = state.a2[positions[found]]
        return chrom_out, idx_out, pa1_out, pa2_out

    def _resolve_with_query_alleles(self, endpoint: pd.DataFrame):
        """Resolve an endpoint and also return its query A1/A2 (or None columns)."""
        chrom, idx, pa1, pa2 = self._resolve_endpoint(endpoint)
        qa1 = endpoint["A1"].to_numpy(dtype=object) if "A1" in endpoint.columns else np.full(len(endpoint), None, dtype=object)
        qa2 = endpoint["A2"].to_numpy(dtype=object) if "A2" in endpoint.columns else np.full(len(endpoint), None, dtype=object)
        return chrom, idx, qa1, qa2, pa1, pa2

    def query_pairs(self, pairs: pd.DataFrame) -> pd.DataFrame:
        """Return unbiased R2 and query-oriented signed r for each SNP pair.

        Parameters
        ----------
        pairs : pandas.DataFrame
            Endpoint columns are suffixed ``_1``/``_2``: ``CHR``/``POS`` in
            coordinate-family modes or ``SNP`` in rsID-family modes, with ``A1``/``A2`` in
            allele-aware modes. Coordinates use the selected R2 matrix build.
            Registered aliases are accepted; input rows are not reordered.

        Returns
        -------
        pandas.DataFrame
            Copy of the input with ``r2`` (float32 unbiased R2), ``sign_r``
            (nullable Int8, +1/-1), ``r`` (float64 Pearson correlation), and
            ``status`` (string). Existing result-named columns are replaced.
            ``sign_r`` and ``r`` use query A1 dosage orientation, allowing
            allele swaps and strand complements. Both are missing in base
            modes.
            Unresolved, cross-chromosome, and unstored pairs have NaN R2 and
            status ``not_in_matrix``, ``cross_chromosome``, or ``absent``.
            Resolved diagonals have R2=1; numeric results have empty status.

        Raises
        ------
        LDSCInputError
            Endpoint columns or the R2 matrix schema, identity, or binding are
            incompatible with the selected mode.
        """
        ep1 = _endpoint_frame(pairs, "1", self._mode)
        ep2 = _endpoint_frame(pairs, "2", self._mode)
        c1, i1, qa1_1, _qa2_1, pa1_1, pa2_1 = self._resolve_with_query_alleles(ep1)
        c2, i2, qa1_2, _qa2_2, pa1_2, pa2_2 = self._resolve_with_query_alleles(ep2)

        n = len(pairs)
        r2 = np.full(n, np.nan, dtype=np.float32)
        sign = np.zeros(n, dtype=np.int8)
        status = np.array([""] * n, dtype=object)

        resolved = (i1 >= 0) & (i2 >= 0)
        status[~resolved] = "not_in_matrix"
        cross = resolved & (c1 != c2)
        status[cross] = "cross_chromosome"
        same = resolved & (c1 == c2)
        diagonal = same & (i1 == i2)
        r2[diagonal] = 1.0
        sign[diagonal] = 1  # R2 matrix orientation; harmonized below
        lookup = same & (i1 != i2)

        for chrom in {c for c in c1[lookup]}:
            rows = np.where(lookup & (c1 == chrom))[0]
            state = self._chrom_state(chrom)
            lo = np.minimum(i1[rows], i2[rows])
            hi = np.maximum(i1[rows], i2[rows])
            chrom_r2, chrom_sign = lookup_pairs_in_parquet(
                state.parquet_file, lo, hi, n_snps=state.n_snps, r2_scale=state.r2_scale,
            )
            r2[rows] = chrom_r2
            sign[rows] = chrom_sign

        found_numeric = ~np.isnan(r2)
        absent = lookup & ~found_numeric
        status[absent] = "absent"

        if is_allele_aware_mode(self._mode):
            mult1 = _orientation_multiplier(qa1_1, pa1_1, pa2_1)
            mult2 = _orientation_multiplier(qa1_2, pa1_2, pa2_2)
            harmonized = sign.astype(np.int16) * mult1 * mult2
            sign = np.where(found_numeric, harmonized, 0).astype(np.int8)
        else:
            sign = np.zeros(n, dtype=np.int8)  # base mode: always NA

        out = pairs.copy()
        out["r2"] = r2
        sign_col = pd.array(sign, dtype="Int8")
        sign_col[sign == 0] = pd.NA
        out["sign_r"] = sign_col
        out["r"] = self._signed_r(r2, sign_col)
        out["status"] = pd.array(status, dtype="string")
        return out

    def _signed_r(self, r2: np.ndarray, sign_col) -> np.ndarray:
        """Recover Pearson r using the harmonized query-orientation sign_r.

        Positive r denotes positive correlation between query A1 dosages;
        query A1 need not be the R2 matrix's minor allele. Base modes carry no
        sign_r, so ``r`` is all-NaN. The required ``ldsc:n_samples`` recovers
        raw R2 from the stored unbiased value.
        """
        n_samples = self.n_samples
        sign_float = sign_col.to_numpy(dtype="float64", na_value=np.nan)
        if not is_allele_aware_mode(self._mode):
            LOGGER.warning(
                "query-r2-matrix returns R2 but all signed r values are unavailable (NaN) because the selected "
                "SNP identifier mode omits allele orientation. To obtain signed r, select an allele-aware "
                "--snp-identifier matching the R2 matrix and provide A1/A2 for both query endpoints."
            )
        return unbiased_r2_to_pearson_r(r2.astype(np.float64), n_samples, sign=sign_float)


def _query_keys(endpoint: pd.DataFrame, mode: str) -> np.ndarray:
    """Build identity keys for query endpoints, tolerating unresolvable rows.

    Unlike :func:`build_snp_id_series`, invalid/ambiguous allele-aware rows do not
    raise: they yield ``None`` so the endpoint resolves to ``not_in_matrix`` (spec
    Section 7) rather than aborting the whole query.
    """
    mode = normalize_snp_identifier_mode(mode)
    base = base_key_series(endpoint, mode, context="query pairs")
    if not is_allele_aware_mode(mode):
        return base.to_numpy(dtype=object)
    allele_set, _reasons = allele_set_series(endpoint, context="query pairs")
    keys = pd.Series(None, index=endpoint.index, dtype=object)
    valid = base.notna() & allele_set.notna()
    keys.loc[valid] = base.loc[valid].astype(str) + ":" + allele_set.loc[valid].astype(str)
    return keys.to_numpy(dtype=object)


_ENDPOINT_COLUMN_SPECS = {
    "CHR": CHR_COLUMN_SPEC,
    "POS": POS_COLUMN_SPEC,
    "SNP": SNP_COLUMN_SPEC,
    "A1": A1_COLUMN_SPEC,
    "A2": A2_COLUMN_SPEC,
}
_ENDPOINT_COLUMNS = tuple(_ENDPOINT_COLUMN_SPECS)


def _endpoint_column_spec(base: ColumnSpec, suffix: str) -> ColumnSpec:
    """Return a ``query-r2-matrix`` endpoint spec for aliases suffixed by endpoint."""
    aliases: list[str] = []
    for alias in base.aliases:
        aliases.extend((f"{alias}_{suffix}", f"{alias}{suffix}"))
    return ColumnSpec(
        base.canonical,
        tuple(aliases),
        f"endpoint-{suffix} {base.label}",
        allow_suffix_match=False,
    )


def _endpoint_frame(pairs: pd.DataFrame, suffix: str, mode: str) -> pd.DataFrame:
    """Extract one endpoint's canonical columns from suffixed pair columns.

    Resolves ``<NAME>_<suffix>`` or ``<alias>_<suffix>`` input columns and returns
    canonical ``CHR/POS/SNP/A1/A2`` fields for matching. In base modes
    (``rsid``/``chr_pos``) allele columns are dropped before use, so a base-mode
    query behaves identically with or without alleles supplied.
    """
    allele_aware = is_allele_aware_mode(mode)
    suffix_token = normalize_column_token(suffix)
    suffixed_columns = [
        column for column in pairs.columns if normalize_column_token(column).endswith(suffix_token)
    ]
    out = {}
    for canonical in _ENDPOINT_COLUMNS:
        if canonical in ("A1", "A2") and not allele_aware:
            continue
        spec = _endpoint_column_spec(_ENDPOINT_COLUMN_SPECS[canonical], suffix)
        col = resolve_optional_column(suffixed_columns, spec, context=f"query-r2-matrix endpoint {suffix}")
        if col is not None:
            out[canonical] = pairs[col].to_numpy()
    if not out:
        raise LDSCInputError(
            f"query-r2-matrix found no endpoint-{suffix} columns. Most likely the pair table is "
            f"missing `CHR_{suffix}/POS_{suffix}` or `SNP_{suffix}`. Provide the columns "
            "required by the R2 matrix's SNP identifier mode."
        )
    return pd.DataFrame(out)


def _orientation_multiplier(query_a1: np.ndarray, panel_a1: np.ndarray, panel_a2: np.ndarray) -> np.ndarray:
    """Return per-endpoint sign multiplier (+1 aligned, -1 swapped, 0 unknown).

    Aligned when the query A1 equals R2 matrix A1 (or its strand complement); swapped
    when it equals R2 matrix A2 (or its complement). Strand-ambiguous R2 matrix SNPs are
    excluded from allele-aware matching, so a matched endpoint always classifies
    cleanly; ``0`` only occurs for unresolved endpoints (R2 matrix allele ``None``).
    """
    n = query_a1.shape[0]
    out = np.zeros(n, dtype=np.int16)
    for k in range(n):
        q, p1, p2 = query_a1[k], panel_a1[k], panel_a2[k]
        if q is None or p1 is None or p2 is None or (isinstance(q, float) and np.isnan(q)):
            continue
        q = str(q).upper()
        p1u, p2u = str(p1).upper(), str(p2).upper()
        if q == p1u or q == p1u.translate(COMPLEMENT):
            out[k] = 1
        elif q == p2u or q == p2u.translate(COMPLEMENT):
            out[k] = -1
    return out


def build_parser() -> argparse.ArgumentParser:
    """Build the ``ldsc query-r2-matrix`` argument parser."""
    parser = argparse.ArgumentParser(
        prog="ldsc query-r2-matrix",
        description="Query unbiased R2 (and signed r) for SNP pairs from an R2 matrix.",
        allow_abbrev=False,
    )
    parser.prog = 'ldsc query-r2-matrix'
    parser.formatter_class = CLIHelpFormatter
    matrix = parser.add_argument_group('Step 1: R2 matrix')
    pairs = parser.add_argument_group('Step 2: SNP pairs')
    output = parser.add_argument_group('Step 3: Output')
    logs = parser.add_argument_group('Logging')

    matrix.add_argument(
        '--r2-matrix-dir', required=True, metavar='DIR',
        help=(
            'Required R2 matrix directory produced by build-r2-matrix, or its root containing '
            'genome-build subdirectories. Use --genome-build to select a subdirectory when needed.'
        ),
    )
    matrix.add_argument(
        '--genome-build', choices=['hg19', 'hg38'], default=None,
        help=(
            'Select the hg19 or hg38 subdirectory under --r2-matrix-dir. If omitted, use the supplied '
            'build-specific directory or its sole available build; ambiguous build roots require a '
            'selection.'
        ),
    )

    pairs.add_argument(
        '--pairs-file', required=True, metavar='FILE',
        help=(
            'Required TSV/CSV of SNP pairs; use - to read standard input. Endpoint columns use _1/_2 '
            'suffixes and must match --snp-identifier, including allele columns for allele-aware matching.'
        ),
    )
    pairs.add_argument(
        '--snp-identifier', default=None,
        help=(
            'SNP matching rule for pair endpoints: the mode recorded by the R2 matrix (the default) or the '
            'base form of that mode (rsid or chr_pos), which ignores alleles. Pair columns must support '
            'the selected mode.'
        ),
    )

    output.add_argument(
        '--output-dir', required=True, metavar='DIR',
        help=(
            'Required destination for query_r2_matrix.tsv (r2, sign_r, r, status plus input columns) '
            'and diagnostics. sign_r and r use query allele orientation.'
        ),
    )
    output.add_argument(
        '--overwrite', action='store_true', default=False,
        help=(
            "Replace this command's existing output files and remove obsolete outputs from an earlier run. "
            'Default: off; stop if output files already exist.'
        ),
    )

    logs.add_argument(
        '--log-level', default='INFO', choices=('DEBUG', 'INFO', 'WARNING', 'ERROR'),
        help=LOG_LEVEL_HELP,
    )
    return parser


@materializing_overwrite_guard(
    args_output_scope,
    command="run_query_r2_matrix_from_args(...)",
)
def run_query_r2_matrix_from_args(args: argparse.Namespace) -> pd.DataFrame:
    """Run ``query-r2-matrix`` from parsed CLI arguments and emit the result table.

    The command writes the canonical result directory: ``query_r2_matrix.tsv`` plus
    ``diagnostics/metadata.json`` and ``diagnostics/query-r2-matrix.log``, matching the
    other workflow commands.
    """
    return _write_query_r2_matrix_directory(args)


def _write_query_r2_matrix_directory(args: argparse.Namespace, pairs: pd.DataFrame | None = None) -> pd.DataFrame:
    """Query against the R2 matrix and write the canonical result directory."""
    output_dir = ensure_output_directory(args.output_dir, label="output directory")
    diagnostics_dir = output_dir / "diagnostics"
    log_path = diagnostics_dir / "query-r2-matrix.log"
    QueryR2MatrixDirectoryWriter.artifact_family(output_dir).preflight(
        overwrite=args.overwrite, additional_paths=[log_path],
    )
    diagnostics_dir.mkdir(parents=True, exist_ok=True)
    with workflow_logging("query-r2-matrix", log_path, log_level=args.log_level):
        from .._runtime.preflight import inspect_declared_inputs
        from .preflight import inspect_r2_headers
        from .._runtime.progress import PhaseProgress
        inspect_declared_inputs(files=[] if args.pairs_file == '-' else [('pairs', args.pairs_file)],
            checks=[('reference', args.r2_matrix_dir, lambda: inspect_r2_headers(
                args.r2_matrix_dir, args.genome_build, command='query-r2-matrix'))],
            issues_path=diagnostics_dir/'input_issues.tsv')
        if pairs is None:
            with PhaseProgress(LOGGER, 'validation', 'pair input loading'):
                pairs_handle = sys.stdin if args.pairs_file == '-' else args.pairs_file
                sep = ',' if str(args.pairs_file).endswith('.csv') else '\t'
                pairs = pd.read_csv(pairs_handle, sep=sep)
        log_inputs(r2_matrix_dir=args.r2_matrix_dir, pairs=args.pairs_file, output_dir=str(output_dir))
        LOGGER.info(f"Querying R2 for pairs from '{args.pairs_file}' against R2 matrix '{args.r2_matrix_dir}'.")
        matrix = R2Matrix.open(args.r2_matrix_dir, snp_identifier=args.snp_identifier, genome_build=args.genome_build)
        result = matrix.query_pairs(pairs)
        written = QueryR2MatrixDirectoryWriter().write(
            result,
            QueryR2MatrixOutputConfig(output_dir=output_dir, overwrite=args.overwrite),
            metadata=_query_r2_matrix_metadata(args, matrix, result),
        )
        LOGGER.info("R2 query complete: pairs=%d; status counts=%s. See query_r2_matrix.tsv and diagnostics/metadata.json.",
                    len(result), result["status"].value_counts(dropna=False).to_dict())
        log_outputs(**written)
    return result


def _query_r2_matrix_metadata(args: argparse.Namespace, matrix: "R2Matrix", result: pd.DataFrame) -> dict[str, object]:
    """Build the diagnostic metadata payload for one query-r2-matrix run."""
    counts = result["status"].value_counts(dropna=False)
    status_counts = {(str(label) or "resolved"): int(count) for label, count in counts.items()}
    return {
        "r2_matrix_dir": str(args.r2_matrix_dir),
        "snp_identifier": matrix.snp_identifier,
        "genome_build": args.genome_build,
        "n_samples": matrix.n_samples,
        "n_pairs": int(len(result)),
        "status_counts": status_counts,
    }


def main(argv: Sequence[str] | None = None) -> pd.DataFrame:
    """Command-line entry point for the R2 pair-query workflow."""
    args = build_parser().parse_args(argv)
    return run_query_r2_matrix_from_args(args)
