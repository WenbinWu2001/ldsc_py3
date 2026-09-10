"""Bounded source parsing and exact disk-backed quantile identity alignment.

SQLite stores narrow identity/row mappings, never fitted annotation matrices.
Its indices are maintained incrementally; queries traverse existing indices
without temporary sorts. All database files belong to the analysis workspace,
and the page cache is limited to 8 MiB. Diagnostics stream to persistent gzip.
"""

import gzip
import sqlite3

import numpy as np
import pandas as pd

from .column_inference import (
    A1_COLUMN_SPEC, A2_COLUMN_SPEC, CHR_COLUMN_SPEC, POS_COLUMN_SPEC,
    SNP_COLUMN_SPEC, resolve_optional_column, resolve_required_column,
)
from .errors import LDSCInputError
from .outputs import SNP_ALIGNMENT_ISSUE_COLUMNS


ROW_CHUNK_SIZE = 65536


def target_chunks(paths, target_annotation, missing_token):
    """Yield selected target metadata, float64 values, and raw-token exclusions."""
    for path in paths:
        try:
            header = pd.read_csv(path, sep=r"\s+", nrows=0).columns
            rename = {
                resolve_required_column(header, CHR_COLUMN_SPEC, context=path): "CHR",
                resolve_required_column(header, POS_COLUMN_SPEC, context=path): "POS",
                resolve_required_column(header, SNP_COLUMN_SPEC, context=path): "SNP",
            }
            if target_annotation not in header:
                raise LDSCInputError(
                    f"Target annotation '{target_annotation}' is absent from '{path}'. Choose an existing column "
                    "with --target-annotation or supply the correct --target-annot-sources."
                )
            a1 = resolve_optional_column(header, A1_COLUMN_SPEC, context=path)
            a2 = resolve_optional_column(header, A2_COLUMN_SPEC, context=path)
            if a1 is not None and a2 is not None:
                rename.update({a1: "A1", a2: "A2"})
            with pd.read_csv(path, sep=r"\s+", dtype=str, keep_default_na=False,
                             usecols=list(dict.fromkeys([*rename, target_annotation])),
                             chunksize=ROW_CHUNK_SIZE) as reader:
                for frame in reader:
                    raw = frame[target_annotation].astype("string")
                    excluded = pd.Series(False, index=frame.index)
                    if missing_token is not None:
                        try:
                            token = float(missing_token)
                        except ValueError:
                            token = np.nan
                        excluded = (pd.to_numeric(raw, errors="coerce").eq(token) if np.isfinite(token)
                                    else raw.str.strip().eq(str(missing_token).strip())).fillna(False)
                    numeric = pd.to_numeric(raw.mask(excluded), errors="coerce").to_numpy(dtype=np.float64, na_value=np.nan)
                    invalid = ~excluded.to_numpy(dtype=bool) & ~np.isfinite(numeric)
                    if invalid.any():
                        examples = sorted(set(raw.loc[invalid].astype(str).head(5)))
                        raise LDSCInputError(
                            f"Target annotation '{target_annotation}' contains nonnumeric, NaN, or infinite value(s): {examples}. "
                            "If one token denotes missingness, pass it with --target-missing-value; otherwise clean the source."
                        )
                    metadata = frame.loc[:, list(rename)].rename(columns=rename).reset_index(drop=True)
                    metadata["CHR"] = metadata.CHR.astype(str).str.replace(r"^chr", "", regex=True)
                    metadata["POS"] = pd.to_numeric(metadata.POS, errors="raise").astype(np.int64)
                    metadata["target_value"] = numeric
                    metadata["target_excluded"] = excluded.to_numpy(dtype=bool)
                    yield metadata
        except (pd.errors.ParserError, UnicodeDecodeError, ValueError) as exc:
            raise LDSCInputError(f"quantile-h2 could not parse target annotation '{path}': {exc}") from exc


def reference_chunks(paths):
    """Yield narrow reference metadata while preserving existing MAF semantics."""
    for path in paths:
        header = pd.read_csv(path, sep=r"\s+", nrows=0).columns
        aliases = {old: new for old, new in {"BP": "POS", "FRQ": "MAF"}.items()
                   if old in header and new not in header}
        normalized = [aliases.get(name, name) for name in header]
        missing = [name for name in ("CHR", "POS", "SNP", "MAF") if name not in normalized]
        if missing:
            raise LDSCInputError(f"Reference metadata '{path}' is missing required column(s) {missing}. Expected CHR, POS (or BP), SNP, and MAF.")
        keep = [name for name in header if aliases.get(name, name) in {"CHR", "POS", "SNP", "MAF", "A1", "A2"}]
        with pd.read_csv(path, sep=r"\s+", usecols=keep, chunksize=ROW_CHUNK_SIZE) as reader:
            for frame in reader:
                frame = frame.rename(columns=aliases)
                frame["CHR"] = frame.CHR.astype(str).str.replace(r"^chr", "", regex=True)
                frame["POS"] = pd.to_numeric(frame.POS, errors="raise").astype(np.int64)
                frame["MAF"] = pd.to_numeric(frame.MAF, errors="coerce")
                yield frame


class AlignmentDiagnostics:
    """Append complete stage-ordered issue records, retaining only counts."""

    def __init__(self, path):
        self.path = path
        self.counts = {}
        self.stream = gzip.open(path, "wt")
        pd.DataFrame(columns=SNP_ALIGNMENT_ISSUE_COLUMNS).to_csv(self.stream, sep="\t", index=False)

    def append(self, frame, *, role, source, issue, action, details):
        if frame.empty:
            return
        rows = frame.reindex(columns=SNP_ALIGNMENT_ISSUE_COLUMNS).copy()
        rows["source_role"], rows["source"] = role, source
        rows["issue"], rows["action"], rows["details"] = issue, action, details
        rows.to_csv(self.stream, sep="\t", index=False, header=False, na_rep="")
        self.counts[issue] = self.counts.get(issue, 0) + len(rows)

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.stream.close()


class QuantileIdentityStore:
    """Analysis-owned exact identity joins with bounded cursor reads.

    Reference and target duplicates remain fatal across chromosomes. Fitted
    row locations point into the annotation bundle; numeric values stay there.
    """

    def __init__(self, path):
        self.db = sqlite3.connect(path)
        self.db.execute("PRAGMA cache_size=-8192")
        self.db.execute("PRAGMA temp_store=MEMORY")
        self.db.execute("PRAGMA journal_mode=DELETE")
        common = "ordinal INTEGER PRIMARY KEY, effective_snp_id TEXT, CHR TEXT, POS INTEGER, SNP TEXT, A1 TEXT, A2 TEXT"
        for table, extra in (("reference", ", base_key TEXT, MAF REAL"),
                             ("fitted", ", shard TEXT, annotation_row INTEGER"),
                             ("target", ", target_value REAL, target_excluded INTEGER")):
            self.db.execute(f"CREATE TABLE {table} ({common}{extra})")
            self.db.execute(f"CREATE INDEX {table}_identity ON {table}(effective_snp_id)")
        self.db.execute("CREATE INDEX reference_base ON reference(base_key)")
        self.db.execute("CREATE INDEX fitted_location ON fitted(shard, annotation_row)")
        self.next_rows = dict.fromkeys(("reference", "fitted", "target"), 0)

    def append(self, table, frame):
        """Insert one bounded metadata chunk in original source order."""
        frame = frame.copy()
        start = self.next_rows[table]
        frame["ordinal"] = np.arange(start, start + len(frame))
        self.next_rows[table] += len(frame)
        # Pandas' default sqlite writer builds Python rows only for this chunk.
        frame.to_sql(table, self.db, if_exists="append", index=False)

    def frames(self, sql, params=()):
        """Read a bounded cursor without collecting the complete result."""
        yield from pd.read_sql_query(sql, self.db, params=params, chunksize=ROW_CHUNK_SIZE)

    def scalar(self, sql, params=()):
        return self.db.execute(sql, params).fetchone()[0]

    def issues(self, diagnostics, table, condition, *, role, source, issue, action, details, params=()):
        """Stream matching rows in source order and return their count."""
        total = 0
        for frame in self.frames(f"SELECT r.* FROM {table} r WHERE {condition} ORDER BY r.ordinal", params):
            diagnostics.append(frame, role=role, source=source, issue=issue, action=action, details=details)
            total += len(frame)
        return total

    def duplicates(self, diagnostics, table, source, *, base=False):
        column = "base_key" if base else "effective_snp_id"
        return self.issues(diagnostics, table,
                           f"EXISTS (SELECT 1 FROM {table} s WHERE s.{column}=r.{column} AND s.ordinal!=r.ordinal)",
                           role="reference_metadata" if table == "reference" else "target_annotation", source=source,
                           issue="ambiguous_allele_inference" if base else "duplicate_identity", action="fatal",
                           details="reference base identity is not unique" if base else "effective SNP identity occurs more than once")

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.db.close()
