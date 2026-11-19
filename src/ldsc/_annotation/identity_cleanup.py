"""Exact global annotation identity cleanup: a hash prefilter with bounded SQLite bookkeeping.

A row can belong to a duplicate or multi-allelic group only when its base key
occurs on another row too. ``IdentityIndex`` therefore keeps one 64-bit hash
of each valid base key in memory, eight bytes per row, and finds the hashes
that repeat. Only rows with a repeated hash enter ``ExactKeyCounts``, a SQLite
table of exact per-key row counts and allele conflicts, so a hash collision
costs one exact comparison and never changes a result. Input in which no hash
repeats creates no database.

The table never stores SNP matrices. Its lookups and upserts use an 8 MiB page
cache. IN lookups build only a bounded key set, never an unbounded sort or
temporary table. SQLite's journal stays beside the owned database; temporary
storage is memory-only.
"""

import sqlite3

import numpy as np
import pandas as pd

from .._runtime.scratch import FrameSpool
from .._core.snp_identity import (
    IDENTITY_DROP_COLUMNS,
    IDENTITY_DROP_REASONS,
    identity_drop_rows,
    allele_set_series,
    base_key_series,
    coerce_identity_drop_frame,
    empty_identity_drop_frame,
    is_allele_aware_mode,
)


_LOOKUP_KEYS = 512
_TRANSACTION_ROWS = 262144


class IdentityDropSpool(FrameSpool):
    """Replay global cleanup records in policy-reason order, then source order."""

    def __init__(self, path):
        super().__init__(path, columns=tuple(IDENTITY_DROP_COLUMNS))
        self._parts = {reason: FrameSpool(path / reason) for reason in IDENTITY_DROP_REASONS}

    def append(self, frame):
        for reason, records in frame.groupby("reason", sort=False):
            self._parts[reason].append(records)
        self.n_rows += len(frame)

    def frames(self):
        for spool in self._parts.values():
            yield from spool.frames()


def _key_hashes(keys):
    """Return one 64-bit hash per base key; equal keys always get equal hashes."""
    return pd.util.hash_array(keys.to_numpy(dtype=object), categorize=False)


class ExactKeyCounts:
    """Exact row count and allele conflict of every base key added, in one SQLite table."""

    def __init__(self, path):
        self._pending_rows = 0
        self.connection = sqlite3.connect(path)
        self.connection.execute("PRAGMA cache_size=-8192")
        self.connection.execute("PRAGMA temp_store=MEMORY")
        self.connection.execute("PRAGMA journal_mode=DELETE")
        self.connection.execute(
            "CREATE TABLE identities (base TEXT PRIMARY KEY, allele TEXT, n INTEGER, multi INTEGER) WITHOUT ROWID"
        )

    def add(self, base, allele):
        """Count one row per base key; ``allele`` is its allele set, or ``''`` when unknown."""
        start = 0
        while start < len(base):
            stop = min(len(base), start + _TRANSACTION_ROWS - self._pending_rows)
            self.connection.executemany(
                "INSERT INTO identities VALUES (?, ?, 1, 0) "
                "ON CONFLICT(base) DO UPDATE SET n=n+1, "
                "multi=multi OR (allele!='' AND excluded.allele!='' AND allele!=excluded.allele), "
                "allele=CASE WHEN allele='' THEN excluded.allele ELSE allele END",
                zip(base.iloc[start:stop], allele.iloc[start:stop]),
            )
            self._pending_rows += stop - start
            start = stop
            if self._pending_rows == _TRANSACTION_ROWS:
                self._finish_transaction()

    def _finish_transaction(self):
        if self._pending_rows:
            self.connection.commit()
            self._pending_rows = 0

    def counts(self, keys):
        """Return ``{key: (rows, allele_conflict)}`` for the added ``keys``, after committing pending adds."""
        self._finish_transaction()
        batch_size = min(_LOOKUP_KEYS, self.connection.getlimit(sqlite3.SQLITE_LIMIT_VARIABLE_NUMBER))
        counts = {}
        for start in range(0, len(keys), batch_size):
            batch = tuple(str(key) for key in keys[start:start + batch_size])
            placeholders = ",".join("?" for _ in batch)
            counts.update((key, (n, multi)) for key, n, multi in self.connection.execute(
                f"SELECT base, n, multi FROM identities WHERE base IN ({placeholders})", batch,
            ))
        return counts

    def close(self):
        self.connection.close()


class IdentityIndex:
    """Index logical annotation rows, then select globally unique rows.

    Aligned column sources describe one logical row: pass its aligned metadata
    once, without loading annotation values. The passes run in this order:
    ``add`` every chunk of every chromosome; call ``has_repeated_keys``; when
    it returns true, replay every chunk through ``add_repeated``; then
    ``select`` each chunk. The database at ``path`` is created only when a
    base-key hash repeats.

    Invalid allele rows do not participate in duplicate or multi-allelic
    detection, matching the in-memory scientific identity policy. A logical
    shard omitting both allele columns uses base identities. Its unknown allele
    set can establish a duplicate, but cannot establish an allele conflict.
    """

    def __init__(self, path, mode):
        self.mode = mode
        self._path = path
        self._hashes = []
        self._repeated = None
        self._exact = None

    def _keys(self, frame):
        base = base_key_series(frame, self.mode, context="annotation")
        if is_allele_aware_mode(self.mode) and {"A1", "A2"}.issubset(frame.columns):
            allele, reasons = allele_set_series(frame, context="annotation")
        else:
            allele = pd.Series("", index=frame.index)
            reasons = pd.Series(None, index=frame.index, dtype=object)
        return base, allele, reasons

    def add(self, frame):
        """Record the base-key hashes of one bounded, already aligned row chunk."""
        base, _, reasons = self._keys(frame)
        self._hashes.append(_key_hashes(base.loc[base.notna() & reasons.isna()]))

    def has_repeated_keys(self):
        """End the first pass and return whether any base-key hash occurs on more than one row."""
        hashes = np.concatenate(self._hashes) if self._hashes else np.empty(0, dtype=np.uint64)
        self._hashes = None
        hashes.sort()
        self._repeated = np.unique(hashes[1:][hashes[1:] == hashes[:-1]])
        return bool(self._repeated.size)

    def _repeated_rows(self, base, valid):
        """Mark the valid rows whose base-key hash repeats somewhere in the complete input."""
        valid = valid.to_numpy(dtype=bool)
        repeated = np.zeros(len(valid), dtype=bool)
        if self._repeated.size and valid.any():
            repeated[valid] = np.isin(_key_hashes(base.loc[valid]), self._repeated)
        return repeated

    def add_repeated(self, frame):
        """Count exactly the rows of one chunk whose base-key hash repeats."""
        base, allele, reasons = self._keys(frame)
        rows = self._repeated_rows(base, base.notna() & reasons.isna())
        if rows.any():
            if self._exact is None:
                self._exact = ExactKeyCounts(self._path)
            self._exact.add(base.loc[rows].astype(str), allele.loc[rows].astype(str))

    def select(self, frame):
        """Return a keep mask and complete drop records for this bounded chunk."""
        base, allele, reasons = self._keys(frame)
        reasons = reasons.copy()
        valid = reasons.isna() & base.notna()
        multi = pd.Series(False, index=frame.index)
        repeated = self._repeated_rows(base, valid)
        if repeated.any():
            counts = self._exact.counts(pd.unique(base.loc[repeated]))
            multi = base.map({key: bool(value[1]) for key, value in counts.items()}).eq(True) & valid
            duplicate = base.map({key: value[0] for key, value in counts.items()}).gt(1) & valid & ~multi
            reasons.loc[multi] = "multi_allelic_base_key"
            reasons.loc[duplicate] = "duplicate_identity"
        if reasons.isna().all():
            return np.ones(len(frame), dtype=bool), empty_identity_drop_frame()
        work = frame.copy()
        work["_ldsc_base_key"] = base
        aware = is_allele_aware_mode(self.mode) and {"A1", "A2"}.issubset(frame.columns)
        work["_ldsc_allele_set"] = allele if aware else pd.NA
        work["_ldsc_identity_key"] = pd.NA
        keys = base.astype("string") + ":" + allele.astype("string") if aware else base
        work.loc[valid & ~multi, "_ldsc_identity_key"] = keys.loc[valid & ~multi]
        records = [
            identity_drop_rows(work.loc[reasons == reason], reason=reason, stage="annotation_identity_cleanup")
            for reason in IDENTITY_DROP_REASONS if bool((reasons == reason).any())
        ]
        return reasons.isna().to_numpy(dtype=bool), coerce_identity_drop_frame(pd.concat(records, ignore_index=True))

    def close(self):
        if self._exact is not None:
            self._exact.close()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()
