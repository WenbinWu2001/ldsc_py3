"""Selective access to aggregate LD-score directories.

Directory opening, manifest checks, and shared table/identity validation
live here, independently of regression orchestration.
Shared baseline values and SNP metadata stay resident for regression alignment.
Query reads return detached batches; the source never caches query values.
The original directory must remain available and unchanged while it is used.
"""

from dataclasses import dataclass, field
import json
import logging
import math
from pathlib import Path
import re
from typing import Any

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

from .._annotation.semantics import require_unique_annotation_names
from ..config import ConfigMismatchError, GlobalConfig
from .._core.snp_identity import DEFAULT_SNP_IDENTIFIER, normalize_snp_identifier_mode
from ..errors import LDSCInputError, LDSCInternalError
from .._core.row_alignment import assert_same_snp_rows, first_identity_mismatch
from .._core.snp_identity import (
    identity_key_series,
    is_allele_aware_mode,
    validate_identity_artifact_metadata,
)
from .._runtime.artifacts import declared_member
from .._runtime.paths import normalize_path_token

_QUERY_BATCH_FILE_RE = re.compile(r"ldscore\.query\.batch\d{5}\.parquet")

# Schema facts of an LD-score directory that the direct and indexed producers,
# the legacy converter and the readers share.
REGRESSION_LD_SCORE_COLUMN = "regression_ld_scores"
ALL_COUNT_FIELD = "all_reference_snp_count"
COMMON_COUNT_FIELD = "common_reference_snp_count"
# ``count_config`` fields: the MAF threshold and comparison that define common-SNP counts.
COMMON_MAF_MIN_FIELD = "common_reference_snp_maf_min"
COMMON_MAF_OPERATOR_FIELD = "common_reference_snp_maf_operator"
INDEX_PROVENANCE_FIELDS = ("index_id", "index_snp_identifier", "index_genome_build")


def count_record(group: str, column: str, all_count, common_count=None) -> dict[str, Any]:
    """Return one entry of the ``counts`` metadata of an LD-score directory.

    ``group`` is ``"baseline"`` or ``"query"``. An unavailable all-SNP count
    (``all_count=None``, legacy conversion without ``.l2.M``) is recorded as
    ``None``; an unavailable common-SNP count (``common_count=None``) is left
    out. Regression treats both as unavailable and never substitutes another
    count kind.
    """
    record = {"group": group, "column": column, ALL_COUNT_FIELD: None if all_count is None else float(all_count)}
    if common_count is not None:
        record[COMMON_COUNT_FIELD] = float(common_count)
    return record


def query_batch_file_name(ordinal: int, count: int) -> str:
    """Return the file name of query batch ``ordinal`` (1-based) in an LD-score directory of ``count`` batches."""
    return "ldscore.query.parquet" if count == 1 else f"ldscore.query.batch{ordinal:05d}.parquet"


def is_numbered_query_batch_file(name: str) -> bool:
    """Return whether ``name`` is a query-batch file of a multi-batch LD-score directory."""
    return _QUERY_BATCH_FILE_RE.fullmatch(name) is not None


LOGGER = logging.getLogger("LDSC.ldscore_source")
_LDSCORE_SCHEMA_DOC = "docs/user-guide/troubleshooting.md#common-ldsc3-artifact-schema-or-provenance-is-incompatible"
# Repair advice closing the reader's errors about a damaged or mixed LD-score directory.
_RESTORE_LDSCORE_DIR = (
    "Restore the complete LD-score directory from the same run, or regenerate it with `ldsc ldscore`. "
    f"Other causes & fixes: {_LDSCORE_SCHEMA_DOC}"
)
_RESTORE_QUERY_BATCHES = (
    " Restore the complete LD-score directory with metadata and query files from the same run, or "
    "regenerate it with ldsc ldscore. Do not mix query batches from different runs."
)


def load_query_manifest(root, metadata, baseline_columns, baseline_schema):
    """Validate the current complete manifest without retaining query values."""
    entries = metadata.get("query_batches")
    if not isinstance(entries, list):
        raise LDSCInputError("LD-score metadata is missing the current query_batches manifest. Regenerate the LD-score directory with the current package.")
    columns, paths = [], set()
    batches = []
    for ordinal, entry in enumerate(entries, 1):
        if not isinstance(entry, dict) or not isinstance(entry.get("query_columns"), list) or not entry["query_columns"] or not isinstance(entry.get("file"), str):
            raise LDSCInputError("LD-score query_batches contains an invalid file or query-column declaration. Regenerate the complete directory.")
        expected_name = query_batch_file_name(ordinal, len(entries))
        if entry["file"] != expected_name or not isinstance(entry.get("row_groups"), list):
            raise LDSCInputError("LD-score query_batches must declare canonical ordered filenames and chromosome row_groups. Regenerate the complete directory.")
        path = root / entry["file"]
        if path in paths:
            raise LDSCInputError(
                f"LD-score query_batches repeats file {entry['file']!r}." + _RESTORE_QUERY_BATCHES
            )
        paths.add(path)
        schema = _read_ldscore_parquet(path, pq.read_schema).names
        _validate_row_groups(
            path, entry["row_groups"],
            f"LD-score query batch {entry['file']} has inconsistent row-group metadata." + _RESTORE_QUERY_BATCHES,
        )
        validate_ldscore_schemas(baseline_columns, entry["query_columns"], baseline_schema, schema)
        columns.extend(entry["query_columns"])
        batches.append({**entry, "path": path, "schema": schema})
    if columns != metadata.get("query_columns", []):
        raise LDSCInputError("LD-score query_batches membership/order disagrees with query_columns. Regenerate the complete directory.")
    require_unique_annotation_names(baseline_columns, columns)
    validate_ldscore_schemas(baseline_columns, columns, baseline_schema,
                             ["CHR", "SNP", "POS", *columns] if batches else None)
    return tuple(batches)


def _validate_row_groups(path, row_groups, message: str) -> None:
    """Raise ``LDSCInputError(message)`` unless ``row_groups`` describes ``path``'s row groups, in order.

    Each record needs its position, its row offset in the file, its row count
    and a chromosome label; one record per chromosome row group.
    """
    parquet_metadata = _read_ldscore_parquet(path, pq.read_metadata)
    if not isinstance(row_groups, list) or len(row_groups) != parquet_metadata.num_row_groups:
        raise LDSCInputError(message)
    offset = 0
    for group_index, group in enumerate(row_groups):
        n_rows = parquet_metadata.row_group(group_index).num_rows
        if (not isinstance(group, dict) or group.get("row_group_index") != group_index
                or group.get("row_offset") != offset or group.get("n_rows") != n_rows
                or not isinstance(group.get("chrom"), str)):
            raise LDSCInputError(message)
        offset += n_rows


def validate_ldscore_schemas(baseline_columns, query_columns, baseline_schema, query_schema):
    """Check every declared column without loading its values."""
    require_unique_annotation_names(baseline_columns, query_columns)
    required = {"CHR", "SNP", "POS", REGRESSION_LD_SCORE_COLUMN, *baseline_columns}
    missing = required.difference(baseline_schema)
    if missing:
        raise LDSCInternalError(f"LD-score baseline_table is missing required columns {sorted(missing)}. Regenerate the complete LD-score directory.")
    if bool(query_columns) != (query_schema is not None):
        raise LDSCInternalError("LD-score query columns and query_table presence disagree. Regenerate the complete LD-score directory.")
    if query_schema is not None:
        missing = {"CHR", "SNP", "POS", *query_columns}.difference(query_schema)
        if missing:
            raise LDSCInternalError(f"LD-score query_table is missing required columns {sorted(missing)}. Regenerate the complete LD-score directory.")


def validate_ldscore_baseline_table(
    table: pd.DataFrame, baseline_columns: list[str], config_snapshot: GlobalConfig | None,
) -> str:
    """Validate the shared output baseline and return its SNP identifier mode.

    Query values may still be on disk. Validate this complete baseline before
    passing its matching row slices to ``validate_ldscore_query_table``.
    """
    if config_snapshot is None:
        raise LDSCInternalError(
            "LD-score output writer could not validate tables because LDScoreResult.config_snapshot is missing. "
            "Most likely an upstream LD-score workflow returned an older or incomplete result object. "
            "Re-run with DEBUG logging and report the traceback."
        )
    snp_identifier = config_snapshot.snp_identifier
    _validate_ldscore_output_columns(
        table, [REGRESSION_LD_SCORE_COLUMN, *baseline_columns],
        table_name="baseline_table", snp_identifier=snp_identifier,
    )
    return snp_identifier


def validate_ldscore_query_table(
    baseline: pd.DataFrame, query: pd.DataFrame | None, query_columns: list[str], *,
    snp_identifier: str,
) -> None:
    """Validate query columns and identity against already-validated baseline rows.

    Accept either a complete materialized query table or one chromosome
    fragment paired with its baseline slice. An absent table requires no
    query columns; a present table must declare at least one query column.
    """
    if bool(query_columns) != (query is not None):
        raise LDSCInternalError(
            "LD-score query columns and query_table presence disagree. Most likely an upstream LD-score "
            "step returned incomplete query annotation metadata. Re-run with DEBUG logging and report the traceback."
        )
    if query is None:
        return
    _validate_ldscore_output_columns(
        query, query_columns, table_name="query_table", snp_identifier=snp_identifier,
    )
    assert_same_snp_rows(
        baseline, query, context="query rows must match baseline rows on CHR/SNP/POS",
        snp_identifier=snp_identifier,
    )


def _validate_ldscore_output_columns(table, columns, *, table_name, snp_identifier):
    """Check output columns and allele identity without comparing paired rows."""
    identity_columns = ["A1", "A2"] if is_allele_aware_mode(snp_identifier) else []
    missing = [name for name in ["CHR", "POS", "SNP", *identity_columns, *columns] if name not in table.columns]
    if missing:
        if any(name in missing for name in ("A1", "A2")):
            raise LDSCInternalError(
                f"LD-score output writer could not validate {table_name} for allele-aware output because "
                f"required columns are missing: {missing}. Most likely an upstream LD-score step dropped A1/A2 "
                "from an allele-aware result table. Regenerate the LD-score directory with the current LDSC3 "
                "package; if this was produced by this run, re-run with DEBUG logging and report the traceback."
            )
        raise LDSCInternalError(
            f"LD-score output writer could not validate {table_name} because required columns are missing: "
            f"{missing}. Most likely an upstream LD-score step returned a result table with the wrong schema. "
            "Re-run with DEBUG logging and report the traceback."
        )
    validate_ldscore_allele_columns(table, table_name=table_name, snp_identifier=snp_identifier)


def validate_ldscore_allele_columns(table: pd.DataFrame, *, table_name: str, snp_identifier: str) -> None:
    """Reject allele-aware LD-score tables without usable allele identity."""
    if not is_allele_aware_mode(snp_identifier):
        return
    missing = [column for column in ("A1", "A2") if column not in table.columns]
    if missing:
        raise LDSCInputError(
            f"LD-score table {table_name} is malformed for snp_identifier='{snp_identifier}' because "
            f"allele columns are missing: {missing}. Most likely this artifact was produced by an older package "
            "version or copied without allele-aware columns. Regenerate the LD-score directory with the current "
            "LDSC3 package."
        )
    try:
        identity_key_series(table, snp_identifier, context=table_name)
    except (LDSCInputError, ValueError) as exc:
        raise LDSCInputError(
            f"LD-score table {table_name} is malformed for snp_identifier='{snp_identifier}' because "
            f"its A1/A2 columns are unusable: {exc}. Most likely the artifact contains malformed allele values "
            "or duplicate allele-aware SNP identities. Regenerate the LD-score directory with the current LDSC3 package."
        ) from exc


@dataclass(frozen=True)
class LDScoreSource:
    """Saved LD scores with resident baseline values and selective query reads.

    Obtain a source from a writing LD-score workflow or
    ``load_ldscore_from_dir``. The source retains shared baseline values, SNP
    metadata, annotation counts and overlaps, and artifact paths. It retains
    no query LD-score values or per-chromosome result tables.

    Attributes
    ----------
    baseline_table : pandas.DataFrame
        Genome-wide regression SNP rows, baseline LD-score columns, and the
        separate ``regression_ld_scores`` weight column.
    query_columns, baseline_columns : list of str
        Ordered annotation names declared by the saved artifact.
    query_batches : tuple of dict
        Validated file, column, and chromosome row-group declarations, with
        resolved paths used by ``read_queries``.
    count_records : list of dict
        All-reference and common-reference annotation counts; these count all
        retained reference SNPs, not only the regression SNPs.
    overlap : LDScoreOverlap or None
        Shared annotation cross-products and SNP totals, when applicable.
    output_paths : dict of str to str
        Paths to the published scientific artifacts.

    Notes
    -----
    The original directory must remain available and unchanged during use.
    Complete artifact schemas and allele metadata are checked when a source
    is opened, including columns that no model will select. Alignment with
    baseline rows and selected-query numerical checks remain fit-time checks.
    Explicit query reads are uncached and have no width limit; the caller is
    responsible for their memory use.
    """

    baseline_table: pd.DataFrame
    query_metadata: pd.DataFrame | None
    query_batches: tuple
    baseline_columns: list[str]
    query_columns: list[str]
    count_records: list[dict]
    config_snapshot: GlobalConfig
    count_config: dict
    overlap: object | None
    output_paths: dict[str, str]
    annotation_types: dict[str, str]
    query_statuses: tuple = ()
    gene_list_batch: object | None = None
    chromosome_scope: dict = field(default_factory=dict)
    snp_universe_policy: dict | None = None
    identity_drops_by_chrom: dict = field(default_factory=dict)
    index_provenance: dict | None = None
    legacy_ldscore_conversion: dict | None = None
    snp_identifier: str = DEFAULT_SNP_IDENTIFIER

    def read_queries(self, columns):
        """Load named query columns across saved files without caching values.

        Parameters
        ----------
        columns : sequence of str
            Query annotation names to load, in the desired output order.
            The selection may span several saved execution batches and has
            no automatic width or memory limit.

        Returns
        -------
        pandas.DataFrame or None
            Stored ``CHR``, ``SNP``, ``POS``, and available ``A1``/``A2``
            columns followed by the selected LD-score columns in requested
            order. SNP rows retain their stored order. An empty selection
            returns SNP metadata only, or ``None`` for a baseline-only source.

        Raises
        ------
        LDSCInputError
            A requested name is not a declared query annotation, or selected
            query files disagree on their effective SNP rows.

        Notes
        -----
        The returned frame owns the loaded query values. This method creates
        no files and retains no values in the source. It checks row alignment
        across selected files; alignment to baseline rows is checked by the
        consuming regression workflow.
        """
        columns = list(columns)
        missing = set(columns).difference(self.query_columns)
        if missing:
            raise LDSCInputError(
                f"Unknown query LD-score columns: {sorted(missing)}."
                " Select query names recorded in metadata.json query_columns. If query annotations were not "
                "generated, rerun ldsc ldscore with the intended baseline and query inputs."
            )
        if not self.query_batches:
            if columns:
                raise LDSCInputError(
                    "LD-score directory has no query table."
                    " Select query names recorded in metadata.json query_columns. If query annotations were "
                    "not generated, rerun ldsc ldscore with the intended baseline and query inputs."
                )
            return None
        metadata = None
        values = {}
        frames = 0
        for batch in self.query_batches:
            selected = [name for name in dict.fromkeys(columns) if name in batch["query_columns"]]
            if not selected:
                continue
            identity = [name for name in ("CHR", "SNP", "POS", "A1", "A2") if name in batch["schema"]]
            frame = pd.read_parquet(batch["path"], columns=[*identity, *selected])
            if metadata is None:
                metadata = frame.loc[:, identity]
            else:
                assert_same_snp_rows(metadata, frame.loc[:, identity],
                                     context="selected query batch files must have identical SNP rows",
                                     snp_identifier=self.snp_identifier)
            values.update((name, frame[name]) for name in selected)
            frames += 1
        if metadata is None:
            return self.query_metadata.copy()
        if frames == 1 and len(set(columns)) == len(columns):
            # One file read in the requested order is already the selection.
            return frame
        return pd.concat([metadata, *(values[name] for name in columns)], axis=1)


@dataclass(frozen=True)
class LDScoreManifest:
    """The checked ``metadata.json`` manifest of one LD-score directory, without its tables.

    Obtain it from ``read_ldscore_manifest``. Every check that needs no
    LD-score values has passed: the identity metadata, the declared files
    and their containment, the baseline and query row groups and schemas,
    the count configuration, the overlap totals of a declared overlap, the
    agreement of the count and overlap configurations, and the complete
    index provenance of an indexed run.

    Attributes
    ----------
    root : pathlib.Path
        The LD-score directory, as given.
    metadata : dict
        The parsed ``metadata.json``.
    config_snapshot : GlobalConfig
        The saved identity configuration.
    snp_identifier : str
        The saved SNP identifier mode, validated.
    genome_build : str or None
        The saved genome build, as written.
    files : dict of str to pathlib.Path
        Every declared file, under ``root``.
    baseline_columns, query_columns : list of str
        Ordered annotation names.
    query_batches : tuple of dict
        Validated query-batch declarations with resolved paths and schemas.
    count_records : list of dict
        All-reference and common-reference annotation counts.
    count_config : dict
        The common-SNP MAF threshold and comparison of the counts.
    overlap_config : dict or None
        The overlap totals and common-SNP rule, when the directory has an overlap.
    annotation_types : dict of str to str
        Per-annotation semantic classification.
    """

    root: Path
    metadata: dict[str, Any]
    config_snapshot: GlobalConfig
    snp_identifier: str
    genome_build: str | None
    files: dict[str, Path]
    baseline_columns: list[str]
    query_columns: list[str]
    query_batches: tuple
    count_records: list[dict]
    count_config: dict
    overlap_config: dict | None
    annotation_types: dict[str, str]

    def load_overlap(self):
        """Read the overlap matrix as an ``LDScoreOverlap``, or return ``None`` when none is declared.

        Raises ``LDSCInputError`` when the overlap file is unreadable or lacks
        a column of the overlap table, or ``overlap_config`` lacks usable
        reference-SNP totals.
        """
        path = self.files.get("overlap")
        if path is None:
            return None
        config = self.overlap_config
        if not _has_overlap_totals(config):
            raise _missing_overlap_totals_error(self.root)
        from .overlap_table import OVERLAP_TABLE_COLUMNS, overlap_from_long_frame
        frame = _read_ldscore_parquet(path, pd.read_parquet)
        missing = [column for column in OVERLAP_TABLE_COLUMNS if column not in frame.columns]
        if missing:
            raise LDSCInputError(
                f"Could not read LD-score file '{path}': the overlap table lacks the column(s) "
                f"{', '.join(f'`{column}`' for column in missing)}. Most likely the file was edited by hand or "
                "written by another tool, or metadata.json and the data files come from different runs. "
                + _RESTORE_LDSCORE_DIR
            )
        return overlap_from_long_frame(
            frame,
            baseline_columns=self.baseline_columns,
            query_columns=self.query_columns,
            total_all_reference_snps=config["total_all_reference_snps"],
            total_common_reference_snps=config.get("total_common_reference_snps"),
        )


def read_ldscore_manifest(ldscore_dir) -> LDScoreManifest:
    """Read and check the manifest of an LD-score directory without reading LD-score values.

    Parameters
    ----------
    ldscore_dir : str or os.PathLike[str]
        Directory written by ``ldsc ldscore``, an indexed run or
        ``convert-legacy-ldscores``.

    Returns
    -------
    LDScoreManifest
        Metadata, identity, declared files, counts, the overlap configuration
        and a lazy overlap read. Only ``metadata.json`` and Parquet footers
        are read.

    Raises
    ------
    LDSCInputError
        The directory or ``metadata.json`` is missing or unreadable, a
        declared file is missing, absolute or outside the directory, a
        Parquet footer is unreadable, the recorded row groups disagree with a
        file, or the identity, query-batch, count, overlap or index
        provenance metadata is invalid.
    LDSCInternalError
        A declared baseline/query schema lacks required columns.
    """
    root = Path(normalize_path_token(ldscore_dir))
    if not root.is_dir():
        raise LDSCInputError(
            f"Could not load LD-score directory '{root}': path is not an existing directory. "
            "Most likely `--ldscore-dir` points at the wrong location. Pass the directory written by `ldsc ldscore`."
        )
    metadata_path = root / "metadata.json"
    if not metadata_path.exists():
        raise LDSCInputError(
            f"Could not load LD-score directory '{root}': missing `metadata.json`. "
            "Most likely this is not an LD-score directory written by `ldsc ldscore` or the diagnostics file "
            "was copied without the data files. Pass the complete LD-score output directory."
        )
    try:
        metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        metadata = exc
    if not isinstance(metadata, dict) or not isinstance(metadata.get("files", {}), dict):
        detail = f": {metadata}" if isinstance(metadata, Exception) else ""
        raise LDSCInputError(
            f"Could not load LD-score directory '{root}': `metadata.json` is not a JSON object with a `files` "
            f"mapping{detail}. Most likely it was truncated, edited by hand or copied from another artifact. "
            + _RESTORE_LDSCORE_DIR
        )
    config_snapshot = _global_config_from_metadata(metadata)
    files = metadata.get("files", {})
    baseline_rel = files.get("baseline")
    if not baseline_rel:
        raise LDSCInputError(
            f"Could not load LD-score directory '{root}': metadata is missing `files.baseline`. "
            "Most likely the LD-score directory was written by an older LDSC3 version or edited by hand. "
            f"Regenerate it with the current `ldsc ldscore`. Other causes & fixes: {_LDSCORE_SCHEMA_DOC}"
        )
    # Only the data files this reader opens must exist; the others are checked for containment.
    declared = {key: _declared_file(root, key, value, required=key in {"baseline", "overlap"})
                for key, value in files.items()}
    _validate_row_groups(
        declared["baseline"], metadata.get("baseline_row_groups"),
        f"LD-score baseline file {baseline_rel} has inconsistent row-group metadata: metadata.json does not "
        "describe its row groups. Most likely an earlier LDSC3 version wrote the directory, splitting a "
        "chromosome of more than 1,048,576 rows into several row groups, or metadata.json and the data files "
        "come from different runs. Regenerate the directory with the current ldsc ldscore. "
        f"Other causes & fixes: {_LDSCORE_SCHEMA_DOC}",
    )
    baseline_columns = [str(column) for column in metadata.get("baseline_columns", [])]
    query_columns = [str(column) for column in metadata.get("query_columns", [])]
    query_batches = load_query_manifest(root, metadata, baseline_columns,
                                        _read_ldscore_parquet(declared["baseline"], pq.read_schema).names)
    count_config = metadata.get("count_config")
    threshold = count_config.get(COMMON_MAF_MIN_FIELD) if isinstance(count_config, dict) else None
    if not _is_number(threshold) or not isinstance(count_config.get(COMMON_MAF_OPERATOR_FIELD), str):
        raise LDSCInputError(
            f"Could not load LD-score directory '{root}': metadata is missing a usable `count_config`, the "
            f"common-SNP MAF threshold (`{COMMON_MAF_MIN_FIELD}`) and comparison "
            f"(`{COMMON_MAF_OPERATOR_FIELD}`) of its counts. Most likely the directory was written by an "
            "older LDSC3 version or edited by hand. Regenerate it with the current `ldsc ldscore`. "
            f"Other causes & fixes: {_LDSCORE_SCHEMA_DOC}"
        )
    _validate_count_overlap_config(metadata, root, has_overlap="overlap" in declared)
    missing_provenance = [key for key in INDEX_PROVENANCE_FIELDS if metadata.get(key) is None]
    if metadata.get("index_id") is not None and missing_provenance:
        raise LDSCInputError(
            f"Could not load LD-score directory '{root}': metadata records `index_id` but no value for "
            f"{', '.join(f'`{key}`' for key in missing_provenance)}, the rest of its gene LD-score index "
            "provenance. Most likely metadata.json was edited by hand or copied from another run. "
            + _RESTORE_LDSCORE_DIR
        )
    return LDScoreManifest(
        root=root,
        metadata=metadata,
        config_snapshot=config_snapshot,
        snp_identifier=config_snapshot.snp_identifier,
        genome_build=metadata.get("genome_build"),
        files=declared,
        baseline_columns=baseline_columns,
        query_columns=query_columns,
        query_batches=query_batches,
        count_records=[dict(record) for record in metadata.get("counts", [])],
        count_config=dict(count_config),
        overlap_config=metadata.get("overlap_config") or None,
        annotation_types={str(key): str(value) for key, value in (metadata.get("annotation_types") or {}).items()},
    )


def _declared_file(root: Path, key: str, value, *, required: bool) -> Path:
    """Return the declared file ``root / value``, raising ``LDSCInputError`` when it cannot be one of ``root``'s files.

    The name must pass ``_runtime.artifacts.declared_member``: containment is
    checked on the declared name, not on its resolved target, so a declared
    file or folder may be a symbolic link, for example to a large table on
    shared storage. ``required`` files must exist as regular files.
    """
    path = declared_member(root, value, external_links=True)
    if path is None:
        problem = ("points outside the directory" if isinstance(value, str) and value.strip()
                   else "is not a file name")
        raise LDSCInputError(
            f"Could not load LD-score directory '{root}': `files.{key}`={value!r} {problem}: declared files are "
            "named relative to the directory, without `..`. Most likely metadata.json was edited by hand or "
            "copied from another run. " + _RESTORE_LDSCORE_DIR
        )
    if (required or path.exists()) and not path.is_file():
        raise LDSCInputError(
            f"Could not load LD-score directory '{root}': `files.{key}`={value!r} declares '{path}', which is "
            "missing or is not a regular file. Most likely the file was deleted or not copied with metadata.json. "
            + _RESTORE_LDSCORE_DIR
        )
    return path


def _is_number(value) -> bool:
    """Return whether a decoded JSON value is a finite number (an ``int`` or ``float``, not a ``bool``).

    ``json.loads`` decodes the literals ``NaN`` and ``Infinity`` to floats,
    which are not usable numbers here.
    """
    return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value)


def _has_overlap_totals(overlap_config) -> bool:
    """Return whether ``overlap_config`` holds usable overlap reference-SNP totals.

    ``total_all_reference_snps`` must be a positive number and
    ``total_common_reference_snps`` null or a number from 0 to
    ``total_all_reference_snps``, since the common reference SNPs are a
    subset of all reference SNPs.
    """
    if not isinstance(overlap_config, dict):
        return False
    total_all = overlap_config.get("total_all_reference_snps")
    total_common = overlap_config.get("total_common_reference_snps")
    return (_is_number(total_all) and total_all > 0
            and (total_common is None or (_is_number(total_common) and 0 <= total_common <= total_all)))


def _missing_overlap_totals_error(root) -> LDSCInputError:
    """Return the error for a declared overlap file whose ``overlap_config`` lacks usable totals."""
    return LDSCInputError(
        f"Could not load LD-score directory '{root}': `files.overlap` is declared, but metadata has no usable "
        "`overlap_config` with the overlap's reference-SNP totals: a positive number `total_all_reference_snps` "
        "and a `total_common_reference_snps` that is null or a number from 0 to `total_all_reference_snps`. "
        "Most likely metadata.json was edited by hand or copied from another run. " + _RESTORE_LDSCORE_DIR
    )


def _read_ldscore_parquet(path, read):
    """Return ``read(path)`` for a Parquet footer or table reader, raising ``LDSCInputError`` when the file is unreadable.

    The error names the file. A missing file's ``FileNotFoundError`` is not
    chained: an input gate reads that cause as a failed path resolution and
    would replace this message with its generic wording.
    """
    try:
        return read(path)
    except (OSError, pa.ArrowException) as exc:
        raise LDSCInputError(
            f"Could not read LD-score file '{path}': {exc}. Most likely it is missing, truncated or not a "
            "Parquet file, or metadata.json and the data files come from different runs. " + _RESTORE_LDSCORE_DIR
        ) from (None if isinstance(exc, FileNotFoundError) else exc)


def load_ldscore_from_dir(
    ldscore_dir: str,
    snp_identifier: str | None = None,
) -> LDScoreSource:
    """Open an LD-score directory with query values read on demand.

    Parameters
    ----------
    ldscore_dir : str
        Directory containing ``metadata.json``, ``ldscore.baseline.parquet``,
        and the query files declared by the current ``query_batches`` manifest.
        Older directories without this manifest must be regenerated.
    snp_identifier : {"rsid", "rsid_allele_aware", "chr_pos", "chr_pos_allele_aware"} or None, optional
        Identifier mode used to check the saved allele identity. When
        omitted, the metadata value is used. An explicit value must agree
        with the saved identity configuration.

    Returns
    -------
    LDScoreSource
        Shared baseline values, SNP metadata, counts, and overlap statistics.
        Call ``read_queries(columns)`` for explicit query selections, including
        selections spanning saved files. Reads preserve requested column order,
        are uncached, and have no width limit. The caller owns their memory
        cost. The original directory must remain available and unchanged.

    Raises
    ------
    LDSCInputError
        The directory or root metadata is missing or unreadable, a declared
        baseline, query or overlap file is missing, unreadable or outside the
        directory, the recorded baseline or query row groups disagree with
        the Parquet files, or the manifest, identity, count, or overlap
        contract is invalid.
    LDSCInternalError
        A declared baseline/query schema lacks required columns.
    ConfigMismatchError
        An explicit SNP identifier mode disagrees with the saved configuration.

    Notes
    -----
    ``read_ldscore_manifest`` checks the manifest first. Opening checks all
    declared schemas and allele metadata without loading query LD-score
    values. Selected-query numerical checks and alignment to baseline rows
    are performed by the consuming regression workflow.
    """
    manifest = read_ldscore_manifest(ldscore_dir)
    root, metadata = manifest.root, manifest.metadata
    LOGGER.info(f"Loading LD-score directory from '{root}'.")
    baseline_table = pd.read_parquet(manifest.files["baseline"])
    overlap = manifest.load_overlap()
    metadata_identifier = manifest.snp_identifier
    if snp_identifier is not None:
        requested_identifier = normalize_snp_identifier_mode(snp_identifier)
        if requested_identifier != metadata_identifier:
            raise ConfigMismatchError(
                "snp_identifier mismatch when loading LD-score directory: "
                f"override {requested_identifier!r} vs metadata {metadata_identifier!r}. "
                "Use an LD-score directory written for the requested SNP identifier mode."
            )
        effective_identifier = requested_identifier
    else:
        effective_identifier = metadata_identifier
    validate_ldscore_allele_columns(
        baseline_table,
        table_name="baseline_table",
        snp_identifier=effective_identifier,
    )
    query_table = None
    identity_columns = ("CHR", "SNP", "POS", "A1", "A2")
    baseline_identity = [name for name in identity_columns if name in baseline_table.columns]
    for entry in manifest.query_batches:
        identities = pd.read_parquet(entry["path"], columns=[
            name for name in identity_columns if name in entry["schema"]])
        # Rows identical to the validated baseline rows, with the same identity
        # columns, need no second allele check.
        if (list(identities.columns) != baseline_identity
                or first_identity_mismatch(baseline_table, identities, baseline_identity) is not None):
            validate_ldscore_allele_columns(
                identities,
                table_name="query_table",
                snp_identifier=effective_identifier,
            )
        if query_table is None:
            query_table = identities
        del identities
    result = LDScoreSource(
        baseline_table=baseline_table,
        query_metadata=query_table,
        query_batches=manifest.query_batches,
        count_records=manifest.count_records,
        baseline_columns=manifest.baseline_columns,
        query_columns=manifest.query_columns,
        output_paths={
            "metadata": str(root / "metadata.json"),
            "baseline": str(manifest.files["baseline"]),
            **{key: str(path) for key, path in manifest.files.items()},
        },
        count_config=manifest.count_config,
        config_snapshot=manifest.config_snapshot,
        overlap=overlap,
        annotation_types=manifest.annotation_types,
        snp_identifier=effective_identifier,
        chromosome_scope=metadata.get("chromosome_scope") or {},
        snp_universe_policy=metadata.get("snp_universe_policy"),
        index_provenance=({key: metadata[key] for key in INDEX_PROVENANCE_FIELDS}
                          if metadata.get("index_id") is not None else None),
        legacy_ldscore_conversion=metadata.get("legacy_ldscore_conversion"),
    )
    query_rows = 0 if query_table is None else len(query_table)
    LOGGER.info(
        f"Loaded LD-score directory '{root}' with {len(baseline_table)} baseline rows, "
        f"{query_rows} query rows, and config provenance present."
    )
    return result


def _validate_count_overlap_config(metadata: dict[str, Any], root: Path, *, has_overlap: bool) -> None:
    """Require count and overlap metadata to describe the same common reference SNPs.

    ``overlap_config`` is a JSON object or null. A declared overlap file
    (``has_overlap``) stores cross-products only, so its ``overlap_config``
    must also hold the reference-SNP totals.
    """
    overlap_config = metadata.get("overlap_config")
    if has_overlap and not _has_overlap_totals(overlap_config):
        raise _missing_overlap_totals_error(root)
    if overlap_config is not None and not isinstance(overlap_config, dict):
        raise LDSCInputError(
            f"Could not load LD-score directory '{root}': metadata `overlap_config` is {overlap_config!r}, not a "
            "JSON object or null. Most likely metadata.json was edited by hand or copied from another run. "
            + _RESTORE_LDSCORE_DIR
        )
    if not overlap_config:
        return
    count_config = metadata.get("count_config") or {}
    count_threshold = count_config.get(COMMON_MAF_MIN_FIELD)
    overlap_threshold = overlap_config.get("common_maf_min")
    count_operator = count_config.get(COMMON_MAF_OPERATOR_FIELD)
    overlap_operator = overlap_config.get("common_maf_operator")
    thresholds_agree = False
    try:
        thresholds_agree = float(count_threshold) == float(overlap_threshold)
    except (TypeError, ValueError):
        thresholds_agree = False
    if not thresholds_agree or count_operator != overlap_operator:
        raise LDSCInputError(
            f"Cannot load LD-score directory '{root}': count_config and overlap_config disagree "
            f"about common-SNP semantics (count threshold/operator {count_threshold!r}/{count_operator!r}, "
            f"overlap {overlap_threshold!r}/{overlap_operator!r}). Reconvert or regenerate the directory so "
            "counts and overlap use the same threshold and operator."
        )


def _global_config_from_metadata(metadata: dict[str, Any]) -> GlobalConfig:
    """Recreate the GlobalConfig snapshot recorded in LD-score metadata."""
    if "artifact_type" not in metadata:
        raise LDSCInputError(
            "Could not read LD-score directory provenance: metadata lacks `artifact_type`. "
            "Most likely the LD-score directory was written by an older LDSC3 version. "
            f"Regenerate it with the current `ldsc ldscore`. Other causes & fixes: {_LDSCORE_SCHEMA_DOC}"
        )
    mode = validate_identity_artifact_metadata(metadata, expected_artifact_type="ldscore")
    return GlobalConfig(
        snp_identifier=mode,
        genome_build=metadata.get("genome_build"),
        log_level="INFO",
    )
