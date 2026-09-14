"""Exact global annotation identity cleanup with bounded SQLite bookkeeping.

The primary-key table stores counts and allele conflicts, never SNP matrices.
Index lookups and upserts use an 8 MiB page cache. IN lookups build only a
bounded key set, never an unbounded sort or temporary table. SQLite's journal
stays beside the owned database; temporary storage is memory-only.
"""

import sqlite3

import pandas as pd

from ._annotation_storage import FrameSpool
from ._kernel.snp_identity import (
    IDENTITY_DROP_REASONS,
    _identity_drop_rows,
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
        super().__init__(path)
        self._parts = {reason: FrameSpool(path / reason) for reason in IDENTITY_DROP_REASONS}

    def append(self, frame):
        for reason, records in frame.groupby("reason", sort=False):
            self._parts[reason].append(records)
        self.n_rows += len(frame)

    def frames(self):
        for spool in self._parts.values():
            yield from spool.frames()


class DiskIdentityIndex:
    """Index logical annotation rows once, then select globally unique rows.

    Aligned column sources describe one logical row and must be combined before
    ``add``. Invalid allele rows do not participate in duplicate or multi-allelic
    detection, matching the in-memory scientific identity policy. A logical
    shard omitting both allele columns uses base identities. Its unknown allele
    set can establish a duplicate, but cannot establish an allele conflict.
    """

    def __init__(self, path, mode):
        self.mode = mode
        self._pending_rows = 0
        self.connection = sqlite3.connect(path)
        self.connection.execute("PRAGMA cache_size=-8192")
        self.connection.execute("PRAGMA temp_store=MEMORY")
        self.connection.execute("PRAGMA journal_mode=DELETE")
        self.connection.execute(
            "CREATE TABLE identities (base TEXT PRIMARY KEY, allele TEXT, n INTEGER, multi INTEGER) WITHOUT ROWID"
        )

    def _keys(self, frame):
        base = base_key_series(frame, self.mode, context="annotation")
        if is_allele_aware_mode(self.mode) and {"A1", "A2"}.issubset(frame.columns):
            allele, reasons = allele_set_series(frame, context="annotation")
        else:
            allele = pd.Series("", index=frame.index)
            reasons = pd.Series(None, index=frame.index, dtype=object)
        return base, allele, reasons

    def add(self, frame):
        """Update global counts from one bounded, already aligned row chunk."""
        base, allele, reasons = self._keys(frame)
        valid = base.notna() & reasons.isna()
        base, allele = base.loc[valid].astype(str), allele.loc[valid].astype(str)
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

    def select(self, frame):
        """Return a keep mask and complete drop records for this bounded chunk."""
        self._finish_transaction()
        base, allele, reasons = self._keys(frame)
        reasons = reasons.copy()
        work = frame.copy()
        work["_ldsc_base_key"] = base
        aware = is_allele_aware_mode(self.mode) and {"A1", "A2"}.issubset(frame.columns)
        work["_ldsc_allele_set"] = allele if aware else pd.NA
        work["_ldsc_identity_key"] = pd.NA
        valid = reasons.isna() & base.notna()
        keys = pd.unique(base.loc[valid])
        batch_size = min(_LOOKUP_KEYS, self.connection.getlimit(sqlite3.SQLITE_LIMIT_VARIABLE_NUMBER))
        counts = {}
        for start in range(0, len(keys), batch_size):
            batch = tuple(str(key) for key in keys[start:start + batch_size])
            placeholders = ",".join("?" for _ in batch)
            counts.update((key, (n, multi)) for key, n, multi in self.connection.execute(
                f"SELECT base, n, multi FROM identities WHERE base IN ({placeholders})", batch,
            ))
        multi = base.map({key: bool(value[1]) for key, value in counts.items()}).eq(True) & valid
        duplicate = base.map({key: value[0] for key, value in counts.items()}).gt(1) & valid & ~multi
        reasons.loc[multi] = "multi_allelic_base_key"
        reasons.loc[duplicate] = "duplicate_identity"
        keys = base.astype("string") + ":" + allele.astype("string") if aware else base
        work.loc[valid & ~multi, "_ldsc_identity_key"] = keys.loc[valid & ~multi]
        records = [
            _identity_drop_rows(work.loc[reasons == reason], reason=reason, stage="annotation_identity_cleanup")
            for reason in IDENTITY_DROP_REASONS if bool((reasons == reason).any())
        ]
        dropped = coerce_identity_drop_frame(pd.concat(records, ignore_index=True)) if records else empty_identity_drop_frame()
        return reasons.isna().to_numpy(dtype=bool), dropped

    def close(self):
        self.connection.close()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()
