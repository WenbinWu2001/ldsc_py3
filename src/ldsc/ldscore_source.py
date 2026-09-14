"""Selective access to canonical aggregate LD-score directories.

Shared baseline values and SNP metadata stay resident for regression alignment.
Query reads return detached batches; the source never caches query values.
The original directory must remain available and unchanged while it is used.
"""

from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd
import pyarrow.parquet as pq

from .annotation_semantics import require_unique_annotation_names
from .config import GlobalConfig
from .errors import LDSCInputError, LDSCInternalError
from ._row_alignment import assert_same_snp_rows


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
        expected_name = "ldscore.query.parquet" if len(entries) == 1 else f"ldscore.query.batch{ordinal:05d}.parquet"
        if entry["file"] != expected_name or not isinstance(entry.get("row_groups"), list):
            raise LDSCInputError("LD-score query_batches must declare canonical ordered filenames and chromosome row_groups. Regenerate the complete directory.")
        path = root / entry["file"]
        if path in paths:
            raise LDSCInputError(f"LD-score query_batches repeats file {entry['file']!r}.")
        paths.add(path)
        schema = pq.read_schema(path).names
        parquet_metadata = pq.read_metadata(path)
        offset = 0
        if len(entry["row_groups"]) != parquet_metadata.num_row_groups:
            raise LDSCInputError(f"LD-score query batch {entry['file']} has inconsistent row-group metadata.")
        for group_index, group in enumerate(entry["row_groups"]):
            n_rows = parquet_metadata.row_group(group_index).num_rows
            if (not isinstance(group, dict) or group.get("row_group_index") != group_index
                    or group.get("row_offset") != offset or group.get("n_rows") != n_rows
                    or not isinstance(group.get("chrom"), str)):
                raise LDSCInputError(f"LD-score query batch {entry['file']} has inconsistent row-group metadata.")
            offset += n_rows
        validate_ldscore_schemas(baseline_columns, entry["query_columns"], baseline_schema, schema)
        columns.extend(entry["query_columns"])
        batches.append({**entry, "path": path, "schema": schema})
    if columns != metadata.get("query_columns", []):
        raise LDSCInputError("LD-score query_batches membership/order disagrees with query_columns. Regenerate the complete directory.")
    require_unique_annotation_names(baseline_columns, columns)
    validate_ldscore_schemas(baseline_columns, columns, baseline_schema,
                             ["CHR", "SNP", "POS", *columns] if batches else None)
    return tuple(batches)


def validate_ldscore_schemas(baseline_columns, query_columns, baseline_schema, query_schema):
    """Check every declared column without loading its values."""
    require_unique_annotation_names(baseline_columns, query_columns)
    required = {"CHR", "SNP", "POS", "regression_ld_scores", *baseline_columns}
    missing = required.difference(baseline_schema)
    if missing:
        raise LDSCInternalError(f"LD-score baseline_table is missing required columns {sorted(missing)}. Regenerate the complete LD-score directory.")
    if bool(query_columns) != (query_schema is not None):
        raise LDSCInternalError("LD-score query columns and query_table presence disagree. Regenerate the complete LD-score directory.")
    if query_schema is not None:
        missing = {"CHR", "SNP", "POS", *query_columns}.difference(query_schema)
        if missing:
            raise LDSCInternalError(f"LD-score query_table is missing required columns {sorted(missing)}. Regenerate the complete LD-score directory.")


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
        All-reference and common-reference annotation counts; these use the
        retained reference universe, not only the regression output rows.
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
    config_snapshot: GlobalConfig | None
    count_config: dict
    overlap: object | None
    output_paths: dict[str, str]
    annotation_types: dict[str, str]
    ld_regression_snps: frozenset[str]
    chromosome_results: tuple = ()
    ld_reference_snps: frozenset[str] = frozenset()
    query_statuses: tuple = ()
    gene_list_batch: object | None = None
    chromosome_scope: dict = field(default_factory=dict)
    snp_universe_policy: dict | None = None
    identity_drops_by_chrom: dict = field(default_factory=dict)
    index_provenance: dict | None = None
    legacy_ldsc2_import: dict | None = None
    snp_identifier: str = "chr_pos_allele_aware"

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
            raise LDSCInputError(f"Unknown query LD-score columns: {sorted(missing)}.")
        if not self.query_batches:
            if columns:
                raise LDSCInputError("LD-score directory has no query table.")
            return None
        metadata = None
        values = []
        for batch in self.query_batches:
            selected = [name for name in dict.fromkeys(columns) if name in batch["query_columns"]]
            if not selected:
                continue
            identity = [name for name in ("CHR", "SNP", "POS", "A1", "A2") if name in batch["schema"]]
            frame = pd.read_parquet(batch["path"], columns=[*identity, *selected])
            current = frame.loc[:, identity]
            if metadata is None:
                metadata = current
            else:
                assert_same_snp_rows(metadata, current, context="selected query batch files must have identical SNP rows",
                                     snp_identifier=self.snp_identifier)
            values.append(frame.loc[:, selected])
        if metadata is None:
            return self.query_metadata.copy()
        return pd.concat([metadata, pd.concat(values, axis=1).loc[:, columns]], axis=1)
