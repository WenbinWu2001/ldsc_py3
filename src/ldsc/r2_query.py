"""Query adjusted R² for SNP pairs from index-format reference panels.

Public surface for reading a ``ldsc build-ref-panel`` parquet panel by SNP pair:
the :class:`R2Panel` handle, the one-shot :func:`query_r2` wrapper, and the pure
:func:`unbiased_r2_to_pearson_r` converter. See
``docs/current/ref-panel-r2-query.md`` and the design spec
``docs/superpowers/specs/2026-06-06-ref-panel-r2-query-design.md``.
"""
from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from ._kernel.identifiers import build_snp_id_series
from ._kernel.ldscore import _load_full_panel_sidecar, _validate_index_binding
from ._kernel.r2_query import lookup_pairs_in_parquet
from ._kernel.ref_panel import (
    _r2_dir_metadata_paths,
    _r2_dir_r2_paths,
    _read_r2_schema_meta,
)
from ._kernel.snp_identity import (
    COMPLEMENT,
    allele_set_series,
    base_key_series,
    identity_mode_family,
    is_allele_aware_mode,
    normalize_snp_identifier_mode,
)
from .chromosome_inference import normalize_chromosome
from .errors import LDSCInputError, LDSCUsageError

LOGGER = logging.getLogger("LDSC.r2_query")


def unbiased_r2_to_pearson_r(r2_adj, n, sign=None):
    """Convert adjusted (unbiased) R² to a signed Pearson correlation r.

    Inverts the unbiased correction ``r2_adj = r2_raw - (1 - r2_raw) / (n - 2)``
    to recover the biased squared correlation, then takes the square root and
    applies ``sign``. Vectorized: ``r2_adj`` and ``sign`` may be scalars or
    array-likes, and ``n`` may be a scalar or array-like reference sample size.

    Parameters
    ----------
    r2_adj : float or array-like
        Adjusted (unbiased) squared correlation as stored in the panel.
    n : int or array-like
        Reference sample size used to estimate R² (``ldsc:n_samples``).
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


@dataclass
class _ChromState:
    """Cached per-chromosome state for one open panel chromosome."""

    sidecar: pd.DataFrame
    key_index: pd.Index          # identity key -> row position (unique)
    a1: np.ndarray               # panel A1 (object), aligned to sidecar rows
    a2: np.ndarray               # panel A2 (object), aligned to sidecar rows
    parquet_file: object         # pyarrow.parquet.ParquetFile
    n_snps: int
    r2_scale: float | None


class R2Panel:
    """Handle over a build-ref-panel index-format R² panel for pair queries.

    Open once with :meth:`open`, then call :meth:`query_pairs` many times. Each
    chromosome is loaded lazily on first touch (sidecar, schema metadata,
    sidecar-binding check, identity-key map, parquet handle), then cached.
    """

    def __init__(self, paths: dict[str, tuple[str, str]], snp_identifier: str) -> None:
        """Store the resolved ``{chrom: (meta, parquet)}`` map and active mode."""
        self._paths = paths
        self._mode = snp_identifier
        self._state: dict[str, _ChromState] = {}
        self._n_samples: int | None = None
        self._n_samples_loaded = False
        self._validate_binding = True

    @classmethod
    def open(
        cls,
        panel_dir=None,
        *,
        meta_path=None,
        parquet_path=None,
        snp_identifier: str | None = None,
        genome_build: str | None = None,
        validate_binding: bool = True,
    ) -> "R2Panel":
        """Open a panel from a directory or an explicit meta/parquet pair."""
        import pyarrow.parquet as pq

        dir_given = panel_dir is not None
        explicit_given = meta_path is not None or parquet_path is not None
        if dir_given == explicit_given:
            raise LDSCUsageError(
                "R2Panel.open requires exactly one input mode: either `panel_dir`, or "
                "both `meta_path` and `parquet_path`. Most likely neither or both were "
                "passed. Pass a build-ref-panel directory, or a single chromosome's "
                "meta+parquet pair."
            )

        paths: dict[str, tuple[str, str]] = {}
        if dir_given:
            r2_paths = _r2_dir_r2_paths(panel_dir, genome_build=genome_build, chrom=None)
            for r2 in r2_paths:
                chrom = _chrom_from_r2_path(r2)
                meta = _r2_dir_metadata_paths(panel_dir, genome_build=genome_build, chrom=chrom)
                if not meta:
                    raise LDSCInputError(
                        f"R2Panel could not find the metadata sidecar for chromosome {chrom} "
                        f"next to '{r2}'. Most likely the `chr{chrom}_meta.tsv.gz` sidecar was "
                        "not copied with the parquet. Restore the sidecar or regenerate the panel."
                    )
                paths[chrom] = (meta[0], r2)
        else:
            if meta_path is None or parquet_path is None:
                raise LDSCUsageError(
                    "R2Panel.open with explicit paths requires both `meta_path` and "
                    "`parquet_path`. Most likely only one was given. Pass both."
                )
            chrom = _chrom_from_r2_path(str(parquet_path))
            paths[chrom] = (str(meta_path), str(parquet_path))

        if not paths:
            raise LDSCInputError(
                "R2Panel.open found no `chr*_r2.parquet` files. Most likely `panel_dir` is "
                "not a build-ref-panel output directory or the wrong genome-build sub-dir "
                "was selected. Pass the directory containing canonical R2 parquet files."
            )

        if snp_identifier is None:
            first_r2 = next(iter(paths.values()))[1]
            schema_meta = pq.ParquetFile(first_r2).schema_arrow.metadata or {}
            raw = schema_meta.get(b"ldsc:snp_identifier")
            if raw is None:
                raise LDSCInputError(
                    f"R2Panel could not read `ldsc:snp_identifier` from '{first_r2}'. Most "
                    "likely the panel predates identity provenance. Pass `snp_identifier=` "
                    "explicitly (rsid, rsid_allele_aware, chr_pos, or chr_pos_allele_aware)."
                )
            snp_identifier = raw.decode("utf-8")

        panel = cls(paths, normalize_snp_identifier_mode(snp_identifier))
        panel._validate_binding = validate_binding
        return panel

    @property
    def chromosomes(self) -> list[str]:
        """Chromosomes available in this panel (in resolved order)."""
        return list(self._paths.keys())

    @property
    def snp_identifier(self) -> str:
        """Active SNP identifier mode used to resolve query SNPs."""
        return self._mode

    @property
    def n_samples(self) -> int | None:
        """Reference sample size from parquet metadata, or ``None`` if absent."""
        if not self._n_samples_loaded:
            first_r2 = next(iter(self._paths.values()))[1]
            self._n_samples = _read_r2_schema_meta(first_r2).n_samples
            self._n_samples_loaded = True
        return self._n_samples

    def _chrom_state(self, chrom: str) -> _ChromState:
        """Load and cache per-chromosome state on first touch."""
        chrom = normalize_chromosome(chrom)
        if chrom in self._state:
            return self._state[chrom]
        if chrom not in self._paths:
            raise LDSCInputError(
                f"R2Panel has no chromosome {chrom}. Available: {self.chromosomes}. Most "
                "likely the query references a chromosome not built into this panel."
            )
        meta_path, r2_path = self._paths[chrom]
        import pyarrow as pa
        import pyarrow.parquet as pq

        sidecar = _load_full_panel_sidecar(r2_path)
        pf = pq.ParquetFile(r2_path)
        schema_meta = pf.schema_arrow.metadata or {}
        n_snps = int(schema_meta[b"ldsc:n_snps"].decode("utf-8"))
        if self._validate_binding:
            identity_hash = schema_meta[b"ldsc:sidecar_identity_sha256"].decode("utf-8")
            _validate_index_binding(
                sidecar, n_snps=n_snps, identity_hash=identity_hash, context=f"R2Panel[{chrom}] {r2_path}"
            )
        key = build_snp_id_series(sidecar, self._mode)
        key_index = pd.Index(key)
        if key_index.has_duplicates:
            raise LDSCInputError(
                f"R2Panel cannot key chromosome {chrom} in mode '{self._mode}': the sidecar "
                "has duplicate SNP identity keys, so a query SNP could resolve ambiguously. "
                "Most likely the mode is too coarse for the panel (e.g. duplicate positions "
                "in chr_pos). Use an allele-aware mode."
            )
        r2_scale = None
        if pa.types.is_integer(pf.schema_arrow.field("R2").type):
            scale_raw = schema_meta.get(b"ldsc:r2_scale")
            r2_scale = float(scale_raw.decode("utf-8")) if scale_raw is not None else 32767.0
        state = _ChromState(
            sidecar=sidecar,
            key_index=key_index,
            a1=sidecar["A1"].to_numpy(dtype=object),
            a2=sidecar["A2"].to_numpy(dtype=object),
            parquet_file=pf,
            n_snps=n_snps,
            r2_scale=r2_scale,
        )
        self._state[chrom] = state
        return state

    def _resolve_endpoint(self, endpoint: pd.DataFrame):
        """Resolve one endpoint table to (chrom, row, panel_a1, panel_a2) arrays.

        ``endpoint`` carries canonical columns ``CHR/POS/SNP/A1/A2`` (subset by
        mode). Returns parallel arrays the length of ``endpoint``: chromosome
        label (object, ``None`` if unresolved), row index (int64, ``-1`` if
        unresolved), and panel A1/A2 (object, ``None`` if unresolved).
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
            for chrom in {c for c in chroms if c in self._paths}:
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
            for chrom in self._paths:
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

    def query_pairs(self, pairs: pd.DataFrame, *, with_r: bool = False,
                    strategy: str = "auto", strategy_threshold: int = 50_000) -> pd.DataFrame:
        """Return adjusted R², sign, and status for each SNP pair in ``pairs``.

        ``pairs`` has per-endpoint columns suffixed ``_1``/``_2``
        (``CHR/POS/A1/A2`` and/or ``SNP``). See the design spec for the full
        contract. ``with_r=True`` appends a signed Pearson ``r`` column.
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
        status[~resolved] = "not_in_panel"
        cross = resolved & (c1 != c2)
        status[cross] = "cross_chromosome"
        same = resolved & (c1 == c2)
        diagonal = same & (i1 == i2)
        r2[diagonal] = 1.0
        sign[diagonal] = 1  # panel orientation; harmonized below
        lookup = same & (i1 != i2)

        for chrom in {c for c in c1[lookup]}:
            rows = np.where(lookup & (c1 == chrom))[0]
            state = self._chrom_state(chrom)
            lo = np.minimum(i1[rows], i2[rows])
            hi = np.maximum(i1[rows], i2[rows])
            chrom_r2, chrom_sign = lookup_pairs_in_parquet(
                state.parquet_file, lo, hi, n_snps=state.n_snps, r2_scale=state.r2_scale,
                strategy=strategy, strategy_threshold=strategy_threshold,
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
        out["sign"] = sign_col
        out["status"] = pd.array(status, dtype="string")
        if with_r:
            out["r"] = self._signed_r(r2, sign_col)
        return out

    def _signed_r(self, r2: np.ndarray, sign_col) -> np.ndarray:
        """Compute the signed r column from adjusted r2 and the sign column."""
        n_samples = self.n_samples
        if n_samples is None:
            raise LDSCInputError(
                "R2Panel cannot compute signed r without `ldsc:n_samples` in the panel "
                "metadata. Most likely a legacy panel was used. Omit `with_r`."
            )
        sign_float = sign_col.to_numpy(dtype="float64", na_value=np.nan)
        if not is_allele_aware_mode(self._mode):
            LOGGER.warning(
                "query-r2 `with_r` in base mode produces an all-NaN r column because base "
                "modes carry no sign. Use an allele-aware panel/mode for signed r."
            )
        return unbiased_r2_to_pearson_r(r2.astype(np.float64), n_samples, sign=sign_float)


def _chrom_from_r2_path(r2_path: str) -> str:
    """Extract the chromosome label from a ``chr{N}_r2.parquet`` path."""
    name = Path(r2_path).name
    if not name.startswith("chr") or not name.endswith("_r2.parquet"):
        raise LDSCInputError(
            f"R2Panel expected a canonical `chr{{N}}_r2.parquet` filename but got '{name}'. "
            "Most likely a non-panel parquet was passed. Pass a build-ref-panel R2 file."
        )
    return normalize_chromosome(name[len("chr") : -len("_r2.parquet")])


def _query_keys(endpoint: pd.DataFrame, mode: str) -> np.ndarray:
    """Build identity keys for query endpoints, tolerating unresolvable rows.

    Unlike :func:`build_snp_id_series`, invalid/ambiguous allele-aware rows do not
    raise: they yield ``None`` so the endpoint resolves to ``not_in_panel`` (spec
    §7) rather than aborting the whole query.
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


_ENDPOINT_COLUMNS = ("CHR", "POS", "SNP", "A1", "A2")


def _endpoint_frame(pairs: pd.DataFrame, suffix: str, mode: str) -> pd.DataFrame:
    """Extract one endpoint's canonical columns from suffixed pair columns.

    Reads ``<NAME>_<suffix>`` columns and renames to canonical ``CHR/POS/SNP/A1/A2``.
    In base modes (``rsid``/``chr_pos``) allele columns are dropped before use, so
    a base-mode query behaves identically with or without alleles supplied.
    """
    allele_aware = is_allele_aware_mode(mode)
    out = {}
    for canonical in _ENDPOINT_COLUMNS:
        if canonical in ("A1", "A2") and not allele_aware:
            continue
        col = f"{canonical}_{suffix}"
        if col in pairs.columns:
            out[canonical] = pairs[col].to_numpy()
    if not out:
        raise LDSCInputError(
            f"query-r2 found no endpoint-{suffix} columns. Most likely the pair table is "
            f"missing `CHR_{suffix}/POS_{suffix}` or `SNP_{suffix}`. Provide the columns "
            "required by the panel's SNP identifier mode."
        )
    return pd.DataFrame(out)


def _orientation_multiplier(query_a1: np.ndarray, panel_a1: np.ndarray, panel_a2: np.ndarray) -> np.ndarray:
    """Return per-endpoint sign multiplier (+1 aligned, -1 swapped, 0 unknown).

    Aligned when the query A1 equals panel A1 (or its strand complement); swapped
    when it equals panel A2 (or its complement). Strand-ambiguous panel SNPs are
    excluded from allele-aware matching, so a matched endpoint always classifies
    cleanly; ``0`` only occurs for unresolved endpoints (panel allele ``None``).
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
