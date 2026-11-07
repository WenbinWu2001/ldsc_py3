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

from .._core.columns import (
    A1_COLUMN_SPEC, A2_COLUMN_SPEC, ANNOTATION_METADATA_SPEC_MAP, CHR_COLUMN_SPEC, MAF_COLUMN_SPEC,
    POS_COLUMN_SPEC, SNP_COLUMN_SPEC, lookalike_hint, resolve_optional_column, resolve_required_column,
)
from .._core.identifiers import open_identity_table
from ..errors import LDSCInputError
from .outputs import SNP_ALIGNMENT_ISSUE_COLUMNS


ROW_CHUNK_SIZE = 65536


def parse_missing_token(token):
    """Return ``--target-missing-value`` in the form :func:`target_chunks` compares target values with.

    A token that parses as a finite number becomes that float, which excludes
    numerically equal values (``-999`` excludes ``-999.0``); any other token,
    such as ``NaN``, becomes its stripped text, which excludes values whose
    stripped text equals it. ``None`` excludes no value.
    """
    if token is None:
        return None
    try:
        number = float(token)
    except ValueError:
        number = np.nan
    return number if np.isfinite(number) else str(token).strip()


def target_chunks(paths, target_annotation, missing):
    """Yield selected target metadata, float64 values, and exclusions.

    ``missing`` is the value :func:`parse_missing_token` returned for the
    ``--target-missing-value`` token. Header names resolve through the
    annotation-file aliases, so ``SNP``, ``A1`` and ``A2`` match only their
    name or an alias, as in every other annotation file.
    """
    specs = ANNOTATION_METADATA_SPEC_MAP
    for path in paths:
        try:
            header = pd.read_csv(path, sep=r"\s+", nrows=0).columns
            rename = {
                resolve_required_column(header, specs["CHR"], context=path): "CHR",
                resolve_required_column(header, specs["POS"], context=path): "POS",
                resolve_required_column(header, specs["SNP"], context=path): "SNP",
            }
            if target_annotation not in header:
                raise LDSCInputError(
                    f"Target annotation '{target_annotation}' is absent from '{path}'. Choose an existing column "
                    "with --target-annot-column or supply the correct --target-annot-sources."
                )
            a1 = resolve_optional_column(header, specs["A1"], context=path)
            a2 = resolve_optional_column(header, specs["A2"], context=path)
            if a1 is not None and a2 is not None:
                rename.update({a1: "A1", a2: "A2"})
            with pd.read_csv(path, sep=r"\s+", dtype=str, keep_default_na=False,
                             usecols=list(dict.fromkeys([*rename, target_annotation])),
                             chunksize=ROW_CHUNK_SIZE) as reader:
                for frame in reader:
                    raw = frame[target_annotation].astype("string")
                    excluded = pd.Series(False, index=frame.index)
                    if missing is not None:
                        excluded = (pd.to_numeric(raw, errors="coerce").eq(missing) if isinstance(missing, float)
                                    else raw.str.strip().eq(missing)).fillna(False)
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
            raise LDSCInputError(
                f"quantile-h2 could not parse target annotation '{path}': {exc}. "
                "Check its text encoding, table columns, CHR labels, and numeric POS coordinates; "
                "supply the intended table with --target-annot-sources. "
                "See docs/user-guide/troubleshooting.md#quantile-h2-target-values-are-invalid-or-quantiles-are-empty."
            ) from exc


_REFERENCE_TROUBLESHOOTING = "docs/user-guide/troubleshooting.md#quantile-h2-reference-metadata-cannot-be-read"


def _reference_columns(header, path):
    """Map reference-metadata header names to CHR, POS, SNP, MAF and, when both exist, A1/A2.

    Each field is taken only by its canonical name or a listed alias, so the
    frequency comes from the column named ``MAF``; a population column such as
    ``EUR_MAF`` is not read as ``MAF``, and the error names it as a column to
    rename. ``FRQ`` is never read as ``MAF``: it usually holds the unfolded A1
    frequency, which would select the wrong common SNPs.
    """
    specs = {spec.canonical: spec for spec in (CHR_COLUMN_SPEC, POS_COLUMN_SPEC, SNP_COLUMN_SPEC, MAF_COLUMN_SPEC)}
    found = {name: resolve_optional_column(header, spec, context=path) for name, spec in specs.items()}
    missing = [name for name, column in found.items() if column is None]
    if missing:
        frq = (" Its FRQ column is not read as MAF, because FRQ usually holds the unfolded A1 frequency, "
               "which would select the wrong common SNPs."
               if "MAF" in missing and any(name.upper() == "FRQ" for name in header) else "")
        raise LDSCInputError(
            f"quantile-h2 reference metadata '{path}' is missing required column(s) {', '.join(missing)}.{frq} "
            "Reference metadata needs CHR, POS (or BP), SNP and MAF, the folded minor-allele frequency. "
            "Set --ref-metadata-sources to the R2 matrix sidecars (chr<N>_meta.tsv.gz) or the files that "
            "ldscore --export-ref-metadata wrote for the panel of the fit; both carry these columns. "
            f"Do not fill missing frequencies with invented values. See {_REFERENCE_TROUBLESHOOTING}."
            + "".join(lookalike_hint(header, specs[name]) for name in missing)
        )
    a1 = resolve_optional_column(header, A1_COLUMN_SPEC, context=path)
    a2 = resolve_optional_column(header, A2_COLUMN_SPEC, context=path)
    if a1 is not None and a2 is not None:
        found.update(A1=a1, A2=a2)
    return {column: name for name, column in found.items()}


def reference_chunks(paths):
    """Yield narrow reference metadata with canonical CHR, POS, SNP, MAF and optional A1/A2 columns.

    Each path is a UTF-8 text table, plain or gzip-compressed, with a
    whitespace-delimited header. Leading ``# ldsc:`` identity lines, which
    R2 matrix sidecars carry, are skipped. Header names resolve through the
    shared column aliases; A1/A2 are kept only when both are present.
    Nonnumeric MAF values become NaN, for the caller's MAF check.
    """
    for path in paths:
        try:
            with open_identity_table(path) as (_, header, stream):
                rename = _reference_columns(header, path)
                with pd.read_csv(stream, sep=r"\s+", header=None, names=header, usecols=list(rename),
                                 chunksize=ROW_CHUNK_SIZE) as reader:
                    for frame in reader:
                        frame = frame.rename(columns=rename)
                        frame["CHR"] = frame.CHR.astype(str).str.replace(r"^chr", "", regex=True)
                        frame["POS"] = pd.to_numeric(frame.POS, errors="raise").astype(np.int64)
                        frame["MAF"] = pd.to_numeric(frame.MAF, errors="coerce")
                        yield frame
        except (pd.errors.ParserError, UnicodeDecodeError, ValueError) as exc:
            raise LDSCInputError(
                f"quantile-h2 could not parse reference metadata '{path}': {exc}. "
                "Check that it is UTF-8 text, plain or gzip-compressed, with a whitespace-delimited header "
                "and numeric POS coordinates; set --ref-metadata-sources to the R2 matrix sidecars "
                "(chr<N>_meta.tsv.gz) or the files that ldscore --export-ref-metadata wrote for the panel "
                f"of the fit. See {_REFERENCE_TROUBLESHOOTING}."
            ) from exc


class AlignmentDiagnostics:
    """Append complete stage-ordered issue records to one gzip table."""

    def __init__(self, path):
        self.path = path
        self.stream = gzip.open(path, "wt")
        pd.DataFrame(columns=SNP_ALIGNMENT_ISSUE_COLUMNS).to_csv(self.stream, sep="\t", index=False)

    def append(self, frame, *, role, source, issue, action, details):
        if frame.empty:
            return
        rows = frame.reindex(columns=SNP_ALIGNMENT_ISSUE_COLUMNS).copy()
        rows["source_role"], rows["source"] = role, source
        rows["issue"], rows["action"], rows["details"] = issue, action, details
        rows.to_csv(self.stream, sep="\t", index=False, header=False, na_rep="")

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
        frame = frame.copy(deep=False)  # the ordinal column is added to this frame only
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
