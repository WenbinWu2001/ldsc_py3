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
from pathlib import Path
from typing import Any

import pandas as pd
import pyarrow.parquet as pq

from ._annotation.semantics import require_unique_annotation_names
from .config import ConfigMismatchError, GlobalConfig
from ._core.columns import normalize_snp_identifier_mode
from .errors import LDSCInputError, LDSCInternalError
from ._core.row_alignment import assert_same_snp_rows
from ._core.identifiers import build_snp_id_series
from ._core.snp_identity import (
    effective_merge_key_series,
    is_allele_aware_mode,
    validate_identity_artifact_metadata,
)
from ._runtime.paths import normalize_path_token


LOGGER = logging.getLogger("LDSC.ldscore_source")
_LDSCORE_SCHEMA_DOC = "docs/user-guide/troubleshooting.md#common-ldsc3-artifact-schema-or-provenance-is-incompatible"


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
            raise LDSCInputError(
                f"LD-score query_batches repeats file {entry['file']!r}."
                " Restore the complete LD-score directory with metadata and query files from the same run, or "
                "regenerate it with ldsc ldscore. Do not mix query batches from different runs."
            )
        paths.add(path)
        schema = pq.read_schema(path).names
        parquet_metadata = pq.read_metadata(path)
        offset = 0
        if len(entry["row_groups"]) != parquet_metadata.num_row_groups:
            raise LDSCInputError(
                f"LD-score query batch {entry['file']} has inconsistent row-group metadata."
                " Restore the complete LD-score directory with metadata and query files from the same run, or "
                "regenerate it with ldsc ldscore. Do not mix query batches from different runs."
            )
        for group_index, group in enumerate(entry["row_groups"]):
            n_rows = parquet_metadata.row_group(group_index).num_rows
            if (not isinstance(group, dict) or group.get("row_group_index") != group_index
                    or group.get("row_offset") != offset or group.get("n_rows") != n_rows
                    or not isinstance(group.get("chrom"), str)):
                raise LDSCInputError(
                    f"LD-score query batch {entry['file']} has inconsistent row-group metadata."
                    " Restore the complete LD-score directory with metadata and query files from the same run, "
                    "or regenerate it with ldsc ldscore. Do not mix query batches from different runs."
                )
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
        table, ["regression_ld_scores", *baseline_columns],
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
        effective_merge_key_series(table, snp_identifier, context=table_name)
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
        Identifier mode used to reconstruct the public regression SNP set from
        the baseline table. When omitted, the metadata value is used. An
        explicit value must agree with the saved identity configuration.

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
        The directory or root metadata is missing, or the manifest, identity,
        count, or overlap contract is invalid.
    LDSCInternalError
        A declared baseline/query schema lacks required columns.
    ConfigMismatchError
        An explicit SNP identifier mode disagrees with the saved configuration.

    Notes
    -----
    Opening checks all declared schemas and allele metadata without loading
    query LD-score values. Selected-query numerical checks and alignment to
    baseline rows are performed by the consuming regression workflow.
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
            "Most likely this is not a canonical `ldsc ldscore` output directory or the diagnostics file "
            "was copied without the data files. Pass the complete LD-score output directory."
        )
    LOGGER.info(f"Loading LD-score directory from '{root}'.")
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    config_snapshot = _global_config_from_metadata(metadata)
    files = metadata.get("files", {})
    baseline_rel = files.get("baseline")
    if not baseline_rel:
        raise LDSCInputError(
            f"Could not load LD-score directory '{root}': metadata is missing `files.baseline`. "
            "Most likely the LD-score directory was written by an older LDSC3 version or edited by hand. "
            f"Regenerate it with the current `ldsc ldscore`. Other causes & fixes: {_LDSCORE_SCHEMA_DOC}"
        )
    baseline_table = pd.read_parquet(root / baseline_rel)
    baseline_columns = [str(column) for column in metadata.get("baseline_columns", [])]
    query_columns = [str(column) for column in metadata.get("query_columns", [])]
    count_records = [dict(record) for record in metadata.get("counts", [])]
    query_batches = load_query_manifest(root, metadata, baseline_columns, baseline_table.columns)
    _validate_count_overlap_config(metadata, root)
    overlap = None
    overlap_rel = files.get("overlap")
    if overlap_rel:
        from .overlap_matrix import overlap_from_long_frame
        overlap_config = metadata.get("overlap_config") or {}
        overlap = overlap_from_long_frame(
            pd.read_parquet(root / overlap_rel),
            baseline_columns=baseline_columns,
            query_columns=query_columns,
            total_all_reference_snps=overlap_config.get("total_all_reference_snps"),
            total_common_reference_snps=overlap_config.get("total_common_reference_snps"),
        )
    metadata_identifier = config_snapshot.snp_identifier if config_snapshot is not None else metadata.get("snp_identifier")
    if snp_identifier is not None:
        requested_identifier = normalize_snp_identifier_mode(snp_identifier)
        if metadata_identifier is not None and requested_identifier != metadata_identifier:
            raise ConfigMismatchError(
                "snp_identifier mismatch when loading LD-score directory: "
                f"override {requested_identifier!r} vs metadata {metadata_identifier!r}. "
                "Use an LD-score directory written for the requested SNP identifier mode."
            )
        effective_identifier = requested_identifier
    else:
        effective_identifier = metadata_identifier or "chr_pos_allele_aware"
    validate_ldscore_allele_columns(
        baseline_table,
        table_name="baseline_table",
        snp_identifier=effective_identifier,
    )
    query_table = None
    for entry in query_batches:
        identities = pd.read_parquet(entry["path"], columns=[
            name for name in ("CHR", "SNP", "POS", "A1", "A2") if name in entry["schema"]])
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
        query_batches=query_batches,
        count_records=count_records,
        baseline_columns=baseline_columns,
        query_columns=query_columns,
        ld_reference_snps=frozenset(),
        ld_regression_snps=frozenset(build_snp_id_series(baseline_table, effective_identifier)),
        chromosome_results=[],
        output_paths={
            "metadata": str(metadata_path),
            "baseline": str(root / baseline_rel),
            **{key: str(root / value) for key, value in files.items()},
        },
        count_config=dict(metadata.get("count_config", {})),
        config_snapshot=config_snapshot,
        overlap=overlap,
        annotation_types={str(key): str(value) for key, value in (metadata.get("annotation_types") or {}).items()},
        snp_identifier=effective_identifier,
        chromosome_scope=metadata.get("chromosome_scope") or {},
        snp_universe_policy=metadata.get("snp_universe_policy"),
        index_provenance=({key: metadata[key] for key in ("index_id", "index_snp_identifier", "index_genome_build")}
                          if metadata.get("index_id") is not None else None),
        legacy_ldsc2_import=metadata.get("legacy_ldsc2_import"),
    )
    query_rows = 0 if query_table is None else len(query_table)
    LOGGER.info(
        f"Loaded LD-score directory '{root}' with {len(baseline_table)} baseline rows, "
        f"{query_rows} query rows, and config provenance {'present' if config_snapshot is not None else 'unknown'}."
    )
    return result


def _validate_count_overlap_config(metadata: dict[str, Any], root: Path) -> None:
    """Require count and overlap metadata to describe the same common reference SNPs."""
    overlap_config = metadata.get("overlap_config")
    if not overlap_config:
        return
    count_config = metadata.get("count_config") or {}
    count_threshold = count_config.get("common_reference_snp_maf_min")
    overlap_threshold = overlap_config.get("common_maf_min")
    count_operator = count_config.get("common_reference_snp_maf_operator")
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


def _global_config_from_metadata(metadata: dict[str, Any]) -> GlobalConfig | None:
    """Recreate a GlobalConfig snapshot when LD-score metadata contains one."""
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
