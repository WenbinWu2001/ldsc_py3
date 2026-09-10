"""Exact global annotation identity cleanup with bounded SQLite bookkeeping.

The primary-key table stores counts and allele conflicts, never SNP matrices.
There are no unbounded SQL sorts or temporary tables: index lookups and upserts
use an 8 MiB page cache. SQLite's journal stays beside the owned database;
temporary storage is memory-only and no statement builds a temporary B-tree.
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
    detection, matching the in-memory scientific identity policy.
    """

    def __init__(self, path, mode):
        self.mode = mode
        self.connection = sqlite3.connect(path)
        self.connection.execute("PRAGMA cache_size=-8192")
        self.connection.execute("PRAGMA temp_store=MEMORY")
        self.connection.execute("PRAGMA journal_mode=DELETE")
        self.connection.execute(
            "CREATE TABLE identities (base TEXT PRIMARY KEY, allele TEXT, n INTEGER, multi INTEGER) WITHOUT ROWID"
        )

    def _keys(self, frame):
        base = base_key_series(frame, self.mode, context="annotation")
        if is_allele_aware_mode(self.mode):
            allele, reasons = allele_set_series(frame, context="annotation")
        else:
            allele = pd.Series("", index=frame.index)
            reasons = pd.Series(None, index=frame.index, dtype=object)
        return base, allele, reasons

    def add(self, frame):
        """Update global counts from one bounded, already aligned row chunk."""
        base, allele, reasons = self._keys(frame)
        valid = base.notna() & reasons.isna()
        self.connection.executemany(
            "INSERT INTO identities VALUES (?, ?, 1, 0) "
            "ON CONFLICT(base) DO UPDATE SET n=n+1, multi=multi OR allele!=excluded.allele",
            zip(base.loc[valid].astype(str), allele.loc[valid].astype(str)),
        )
        self.connection.commit()

    def select(self, frame):
        """Return a keep mask and complete drop records for this bounded chunk."""
        base, allele, reasons = self._keys(frame)
        reasons = reasons.copy()
        work = frame.copy()
        work["_ldsc_base_key"] = base
        work["_ldsc_allele_set"] = allele if is_allele_aware_mode(self.mode) else pd.NA
        work["_ldsc_identity_key"] = pd.NA
        valid = reasons.isna() & base.notna()
        counts = {
            str(key): self.connection.execute("SELECT n, multi FROM identities WHERE base=?", (str(key),)).fetchone()
            for key in pd.unique(base.loc[valid])
        }
        multi = base.map({key: bool(value[1]) for key, value in counts.items()}).eq(True) & valid
        duplicate = base.map({key: value[0] for key, value in counts.items()}).gt(1) & valid & ~multi
        reasons.loc[multi] = "multi_allelic_base_key"
        reasons.loc[duplicate] = "duplicate_identity"
        keys = base.astype("string") + ":" + allele.astype("string") if is_allele_aware_mode(self.mode) else base
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
