"""Output writers for the refactored LDSC workflows.

The public LD-score workflow writes one canonical result directory containing
``metadata.json``, ``ldscore.baseline.parquet``, and optional
``ldscore.query.parquet``. Run identity comes from the chosen directory name;
output filenames inside that directory are fixed. The parquet payloads are
written with one row group per chromosome and matching metadata so
downstream readers can load a single chromosome without scanning the whole
table. The writer creates missing directories, reuses existing directories, and
refuses existing canonical family files unless the caller explicitly sets
``overwrite=True``; successful overwrites remove stale owned siblings that the
current result did not produce.

Partitioned-h2 regression summaries use the same directory-oriented output
policy. ``PartitionedH2DirectoryWriter`` always writes the compact aggregate
``partitioned_h2.tsv`` and can optionally stage a per-query
``diagnostics/query_annotations`` tree with ``manifest.tsv``, one-row
``partitioned_h2.tsv`` summaries, full ``partitioned_h2_full.tsv`` category
tables, and ``metadata.json`` files before moving it into place.
"""

from __future__ import annotations

import json
import os
import re
import shutil
import tempfile
import unicodedata
from dataclasses import asdict, dataclass, is_dataclass
from os import PathLike
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

from ._row_alignment import assert_same_snp_rows
from ._kernel.snp_identity import effective_merge_key_series, identity_artifact_metadata, is_allele_aware_mode
from .config import _normalize_required_path
from .path_resolution import (
    ensure_output_directory,
    preflight_output_artifact_family,
    remove_output_artifacts,
)


def _cast_parquet_floats(df: pd.DataFrame) -> pd.DataFrame:
    """Return a parquet-ready frame with float64 columns narrowed to float32."""
    float64_cols = [c for c in df.columns if df[c].dtype == np.float64]
    return df.astype({c: np.float32 for c in float64_cols}) if float64_cols else df


def _write_chromosome_aligned_parquet(
    df: pd.DataFrame, path: Path, compression: str | None
) -> list[dict]:
    """Write one parquet row group per chromosome and return row-group metadata.

    The input frame must contain a ``CHR`` column and is assumed to already be
    sorted by genomic position. ``groupby(sort=False)`` preserves that order,
    giving row groups whose offsets match the public row order in the written
    table. Float columns are narrowed through ``_cast_parquet_floats`` before
    schema construction.
    """
    df = _cast_parquet_floats(df)
    schema = pa.Schema.from_pandas(df, preserve_index=False)
    row_group_meta: list[dict] = []
    offset = 0
    with pq.ParquetWriter(path, schema, compression=compression) as writer:
        for chrom, chrom_df in df.groupby("CHR", sort=False):
            writer.write_table(pa.Table.from_pandas(chrom_df, preserve_index=False))
            row_group_meta.append(
                {
                    "chrom": str(chrom),
                    "row_group_index": len(row_group_meta),
                    "row_offset": offset,
                    "n_rows": len(chrom_df),
                }
            )
            offset += len(chrom_df)
    return row_group_meta


REGRESSION_LD_SCORE_COLUMN = "regression_ld_scores"
DEFAULT_COUNT_CONFIG = {
    "common_reference_snp_maf_min": 0.05,
    "common_reference_snp_maf_operator": ">=",
}
PARTITIONED_H2_AGGREGATE_COLUMNS = [
    "Category",
    "Prop._SNPs",
    "Prop._h2",
    "Enrichment",
    "Enrichment_p",
    "Coefficient",
    "Coefficient_p",
]
PARTITIONED_H2_FULL_COLUMNS = [
    "Category",
    "Prop._SNPs",
    "Category_h2",
    "Category_h2_std_error",
    "Prop._h2",
    "Prop._h2_std_error",
    "Enrichment",
    "Enrichment_std_error",
    "Enrichment_p",
    "Coefficient",
    "Coefficient_std_error",
    "Coefficient_p",
]
RG_CONCISE_COLUMNS = [
    "trait_1",
    "trait_2",
    "n_snps_used",
    "rg",
    "rg_se",
    "p",
    "p_fdr_bh",
    "note",
]
RG_FULL_COLUMNS = [
    "trait_1",
    "trait_2",
    "n_snps_used",
    "rg",
    "rg_se",
    "z",
    "p",
    "p_fdr_bh",
    "p_bonferroni",
    "h2_1",
    "h2_1_se",
    "h2_2",
    "h2_2_se",
    "gencov",
    "gencov_se",
    "intercept_h2_1",
    "intercept_h2_1_se",
    "intercept_h2_2",
    "intercept_h2_2_se",
    "intercept_gencov",
    "intercept_gencov_se",
    "ratio_1",
    "ratio_1_se",
    "ratio_2",
    "ratio_2_se",
    "lambda_gc_1",
    "lambda_gc_2",
    "mean_chisq_1",
    "mean_chisq_2",
    "pair_kind",
    "status",
    "error",
]


@dataclass(frozen=True)
class H2OutputConfig:
    """Directory-oriented output config for unpartitioned h2 summaries.

    Parameters
    ----------
    output_dir : str or os.PathLike[str]
        Directory that receives ``h2.tsv`` and diagnostic metadata.
    overwrite : bool, optional
        If ``True``, replace existing fixed h2 outputs. Default is ``False``.
    """

    output_dir: str | PathLike[str]
    overwrite: bool = False

    def __post_init__(self) -> None:
        """Normalize the output directory."""
        object.__setattr__(self, "output_dir", _normalize_required_path(self.output_dir))


class H2DirectoryWriter:
    """Write the unpartitioned h2 summary and diagnostic metadata sidecar."""

    def write(
        self,
        summary: pd.DataFrame,
        output_config: H2OutputConfig,
        *,
        metadata: dict[str, object],
    ) -> dict[str, str]:
        """Write ``h2.tsv`` and ``diagnostics/metadata.json`` to a result directory.

        Existing fixed h2 artifacts are checked before any output file is
        written. Replacement requires ``output_config.overwrite=True``.
        """
        output_dir = ensure_output_directory(output_config.output_dir, label="output directory")
        summary_path = output_dir / "h2.tsv"
        diagnostics_dir = output_dir / "diagnostics"
        metadata_path = diagnostics_dir / "metadata.json"
        legacy_metadata_path = output_dir / "metadata.json"
        preflight_output_artifact_family(
            [summary_path, metadata_path],
            [summary_path, metadata_path, legacy_metadata_path],
            overwrite=output_config.overwrite,
            label="h2 output artifact",
        )
        diagnostics_dir.mkdir(parents=True, exist_ok=True)
        _atomic_write_dataframe(summary, summary_path)
        _atomic_write_json(
            _result_metadata(metadata, artifact_type="h2_result", files={"summary": "h2.tsv"}),
            metadata_path,
        )
        return {
            "summary": str(summary_path),
            "metadata": str(metadata_path),
        }


@dataclass(frozen=True)
class LDScoreOutputConfig:
    """Directory-oriented output config for canonical LD-score results.

    Parameters
    ----------
    output_dir : str or os.PathLike[str]
        Directory that receives ``metadata.json``,
        ``ldscore.baseline.parquet``, and optional ``ldscore.query.parquet``.
    overwrite : bool, optional
        If ``True``, replace existing canonical LD-score files in
        ``output_dir`` and remove stale owned siblings after a successful
        write. If ``False``, any existing canonical family file raises
        ``FileExistsError`` before any parquet or metadata file is written.
        Default is ``False``.
    parquet_compression : {"snappy", "gzip", "brotli", "zstd", "none", None}, optional
        Compression codec used for parquet tables. Default is ``"snappy"``.
    """

    output_dir: str | PathLike[str]
    overwrite: bool = False
    parquet_compression: str = "snappy"

    def __post_init__(self) -> None:
        """Normalize the output directory and validate parquet compression."""
        object.__setattr__(self, "output_dir", _normalize_required_path(self.output_dir))
        allowed = {"snappy", "gzip", "brotli", "zstd", "none", None}
        if self.parquet_compression not in allowed:
            raise ValueError("parquet_compression must be one of 'snappy', 'gzip', 'brotli', 'zstd', 'none', or None.")


class LDScoreDirectoryWriter:
    """Write canonical LD-score result directories.

    The writer owns the fixed files ``metadata.json``,
    ``ldscore.baseline.parquet``, and optional ``ldscore.query.parquet``.
    Parquet files are flat files for backward compatibility, but their internal
    row groups are chromosome-aligned and described in root metadata.
    """

    def write(self, result: Any, output_config: LDScoreOutputConfig) -> dict[str, str]:
        """Write ``metadata.json`` and canonical LD-score parquet files.

        Existing canonical family files are checked before any output file is
        written. Replacement requires ``output_config.overwrite=True``;
        unrelated files in the directory are ignored. The returned paths map
        includes ``"baseline"``, ``"metadata"``, and ``"query"`` when query
        annotations were supplied.
        """
        output_dir = ensure_output_directory(output_config.output_dir, label="LD-score output directory")
        baseline_table = getattr(result, "baseline_table", None)
        query_table = getattr(result, "query_table", None)
        if baseline_table is None:
            raise ValueError("LDScoreResult is missing baseline_table.")
        self._validate_tables(result)

        paths = {
            "metadata": output_dir / "metadata.json",
            "baseline": output_dir / "ldscore.baseline.parquet",
        }
        if query_table is not None:
            paths["query"] = output_dir / "ldscore.query.parquet"
        stale_paths = preflight_output_artifact_family(
            paths.values(),
            _ldscore_output_family(output_dir),
            overwrite=output_config.overwrite,
            label="LD-score output artifact",
        )

        compression = None if output_config.parquet_compression in {None, "none"} else output_config.parquet_compression
        baseline_rg = _write_chromosome_aligned_parquet(baseline_table, paths["baseline"], compression)
        query_rg = None
        if query_table is not None:
            query_rg = _write_chromosome_aligned_parquet(query_table, paths["query"], compression)
        metadata = self.build_metadata(
            result,
            files={name: path.name for name, path in paths.items() if name != "metadata"},
            baseline_rg=baseline_rg,
            query_rg=query_rg,
        )
        paths["metadata"].write_text(json.dumps(_to_serializable(metadata), indent=2, sort_keys=True), encoding="utf-8")
        remove_output_artifacts(stale_paths)
        return {name: str(path) for name, path in paths.items()}

    def build_metadata(
        self,
        result: Any,
        files: dict[str, str],
        baseline_rg: list[dict] | None = None,
        query_rg: list[dict] | None = None,
    ) -> dict[str, Any]:
        """Build the JSON metadata payload for one LD-score result.

        The metadata always includes ``count_config`` so downstream regression
        code can report the common-SNP count threshold even when MAF metadata is
        unavailable and per-column common counts are omitted. Row-group metadata
        records the chromosome, row-group index, row offset, and row count for
        each chromosome-aligned parquet row group.
        """
        baseline_table = getattr(result, "baseline_table")
        query_table = getattr(result, "query_table", None)
        config_snapshot = getattr(result, "config_snapshot", None)
        if config_snapshot is None:
            raise ValueError("LDScoreResult.config_snapshot is required to write current LD-score artifacts.")
        identity_metadata = identity_artifact_metadata(
            artifact_type="ldscore",
            snp_identifier=config_snapshot.snp_identifier,
            genome_build=getattr(config_snapshot, "genome_build", None),
        )
        chromosomes = baseline_table["CHR"].astype(str).drop_duplicates().tolist()
        return {
            **identity_metadata,
            "files": dict(files),
            "chromosomes": chromosomes,
            "baseline_columns": list(getattr(result, "baseline_columns", [])),
            "query_columns": list(getattr(result, "query_columns", [])),
            "counts": list(getattr(result, "count_records", [])),
            "count_config": dict(getattr(result, "count_config", None) or DEFAULT_COUNT_CONFIG),
            "n_baseline_rows": int(len(baseline_table)),
            "n_query_rows": 0 if query_table is None else int(len(query_table)),
            "row_group_layout": "one_per_chromosome",
            "baseline_row_groups": baseline_rg or [],
            "query_row_groups": query_rg,
        }

    def _validate_tables(self, result: Any) -> None:
        """Validate baseline/query table shape before any files are written."""
        baseline_table = getattr(result, "baseline_table")
        query_table = getattr(result, "query_table", None)
        baseline_columns = list(getattr(result, "baseline_columns", []))
        query_columns = list(getattr(result, "query_columns", []))
        config_snapshot = getattr(result, "config_snapshot", None)
        if config_snapshot is None:
            raise ValueError("LDScoreResult.config_snapshot is required to write current LD-score artifacts.")
        snp_identifier = config_snapshot.snp_identifier
        identity_columns = ["A1", "A2"] if is_allele_aware_mode(snp_identifier) else []
        required_baseline = ["CHR", "POS", "SNP", *identity_columns, REGRESSION_LD_SCORE_COLUMN, *baseline_columns]
        missing = [column for column in required_baseline if column not in baseline_table.columns]
        if missing:
            if any(column in missing for column in ("A1", "A2")):
                raise ValueError(
                    f"baseline_table is missing required columns: {missing}; allele-aware LD-score artifacts "
                    "require A1/A2 columns. Regenerate it with the current LDSC package."
                )
            raise ValueError(f"baseline_table is missing required columns: {missing}")
        _validate_ldscore_allele_columns(baseline_table, table_name="baseline_table", snp_identifier=snp_identifier)
        if query_columns and query_table is None:
            raise ValueError("query_table is required when query_columns are present.")
        if not query_columns and query_table is not None:
            raise ValueError("query_table was provided but query_columns is empty.")
        if query_table is None:
            return
        required_query = ["CHR", "POS", "SNP", *identity_columns, *query_columns]
        missing = [column for column in required_query if column not in query_table.columns]
        if missing:
            if any(column in missing for column in ("A1", "A2")):
                raise ValueError(
                    f"query_table is missing required columns: {missing}; allele-aware LD-score artifacts "
                    "require A1/A2 columns. Regenerate it with the current LDSC package."
                )
            raise ValueError(f"query_table is missing required columns: {missing}")
        _validate_ldscore_allele_columns(query_table, table_name="query_table", snp_identifier=snp_identifier)
        assert_same_snp_rows(
            baseline_table,
            query_table,
            context="query rows must match baseline rows on CHR/SNP/POS",
            snp_identifier=snp_identifier,
        )


def _validate_ldscore_allele_columns(table: pd.DataFrame, *, table_name: str, snp_identifier: str) -> None:
    """Reject allele-aware LD-score artifacts without usable allele identity."""
    if not is_allele_aware_mode(snp_identifier):
        return
    missing = [column for column in ("A1", "A2") if column not in table.columns]
    if missing:
        raise ValueError(
            f"{table_name} is malformed for snp_identifier='{snp_identifier}': allele-aware LD-score artifacts "
            "require A1/A2 columns. Regenerate it with the current LDSC package."
        )
    try:
        effective_merge_key_series(table, snp_identifier, context=table_name)
    except ValueError as exc:
        raise ValueError(
            f"{table_name} has unusable A1/A2 columns for snp_identifier='{snp_identifier}': {exc}"
        ) from exc


@dataclass(frozen=True)
class PartitionedH2OutputConfig:
    """Directory-oriented output config for partitioned-h2 regression summaries.

    The config keeps the partitioned-h2 output contract explicit and separate
    from LD-score parquet output. The default mode writes only the stable
    aggregate summary. ``write_per_query_results=True`` enables the richer
    per-query tree without changing the aggregate file name.

    Parameters
    ----------
    output_dir : str or os.PathLike[str]
        Directory that receives ``partitioned_h2.tsv`` and, when requested,
        the optional ``diagnostics/query_annotations`` per-query result tree.
    overwrite : bool, optional
        If ``True``, replace existing fixed partitioned-h2 outputs and remove
        stale owned siblings after a successful write. Default is ``False``.
    write_per_query_results : bool, optional
        If ``True``, also write one subdirectory per query annotation under
        ``diagnostics/query_annotations``. Default is ``False``.
    """

    output_dir: str | PathLike[str]
    overwrite: bool = False
    write_per_query_results: bool = False

    def __post_init__(self) -> None:
        """Normalize the output directory."""
        object.__setattr__(self, "output_dir", _normalize_required_path(self.output_dir))


class PartitionedH2DirectoryWriter:
    """Write aggregate and optional per-query partitioned-h2 outputs.

    The writer owns the fixed regression summary layout:
    ``partitioned_h2.tsv`` at the output root uses the compact public schema,
    while optional diagnostic query folders contain the same one-row summary
    plus a full baseline-plus-query ``partitioned_h2_full.tsv`` table. The
    per-query tree
    is written to a temporary sibling directory before it is moved into the
    final location, so ordinary validation and I/O failures do not expose a
    partially populated final tree.
    """

    def write(
        self,
        summary: pd.DataFrame,
        output_config: PartitionedH2OutputConfig,
        per_query_category_tables: dict[str, pd.DataFrame] | None = None,
        metadata: dict[str, object] | None = None,
        per_query_metadata: dict[str, dict[str, object]] | None = None,
    ) -> dict[str, str]:
        """Write partitioned-h2 summary artifacts.

        Parameters
        ----------
        summary : pandas.DataFrame
            Aggregate partitioned-h2 table with the compact public columns in
            ``PARTITIONED_H2_AGGREGATE_COLUMNS``.
        output_config : PartitionedH2OutputConfig
            Output directory, overwrite policy, and per-query mode.
        per_query_category_tables : dict of str to pandas.DataFrame, optional
            Optional full baseline-plus-query category tables keyed by
            original query annotation name. Tables must follow
            ``PARTITIONED_H2_FULL_COLUMNS``; missing keys write empty
            ``partitioned_h2_full.tsv`` files.
        metadata : dict, optional
            Run-level metadata copied into every per-query ``metadata.json``.
        per_query_metadata : dict of str to dict, optional
            Query-specific metadata copied into the matching
            ``metadata.json``.

        Returns
        -------
        dict of str to str
            Written path map. Always includes ``"summary"`` and additionally
            includes ``"per_query_root"`` and ``"per_query_manifest"`` when
            per-query output is enabled.

        Raises
        ------
        ValueError
            If ``summary`` lacks the required ``Category`` column.
        FileExistsError
            If a final output path already exists and overwrite is disabled.
        """
        self._validate_summary(summary)
        output_dir = ensure_output_directory(output_config.output_dir, label="output directory")
        summary_path = output_dir / "partitioned_h2.tsv"
        diagnostics_dir = output_dir / "diagnostics"
        metadata_path = diagnostics_dir / "metadata.json"
        query_root = diagnostics_dir / "query_annotations"
        legacy_metadata_path = output_dir / "metadata.json"
        legacy_query_root = output_dir / "query_annotations"
        produced_paths = [summary_path, metadata_path]
        if output_config.write_per_query_results:
            produced_paths.append(query_root)
        stale_paths = preflight_output_artifact_family(
            produced_paths,
            [summary_path, metadata_path, query_root, legacy_metadata_path, legacy_query_root],
            overwrite=output_config.overwrite,
            label="partitioned-h2 output artifact",
        )

        root_files = {"summary": "partitioned_h2.tsv"}
        if output_config.write_per_query_results:
            root_files["query_annotations"] = "diagnostics/query_annotations"
        paths = {"summary": str(summary_path), "metadata": str(metadata_path)}
        diagnostics_dir.mkdir(parents=True, exist_ok=True)
        if not output_config.write_per_query_results:
            _atomic_write_dataframe(
                _select_columns(summary, PARTITIONED_H2_AGGREGATE_COLUMNS, label="partitioned-h2 summary"),
                summary_path,
            )
            _atomic_write_json(
                _result_metadata(metadata, artifact_type="partitioned_h2_result", files=root_files),
                metadata_path,
            )
            remove_output_artifacts(stale_paths)
            return paths

        query_records = self._query_records(summary)
        staging_dir = Path(tempfile.mkdtemp(prefix=".query_annotations.tmp.", dir=str(diagnostics_dir)))
        backup_dir: Path | None = None
        try:
            manifest_rows = self._write_staged_query_tree(
                staging_dir,
                query_records,
                summary,
                per_query_category_tables or {},
                metadata or {},
                per_query_metadata or {},
            )
            _atomic_write_dataframe(pd.DataFrame(manifest_rows), staging_dir / "manifest.tsv")
            _atomic_write_dataframe(
                _select_columns(summary, PARTITIONED_H2_AGGREGATE_COLUMNS, label="partitioned-h2 summary"),
                summary_path,
            )
            _atomic_write_json(
                _result_metadata(metadata, artifact_type="partitioned_h2_result", files=root_files),
                metadata_path,
            )
            if output_config.overwrite and query_root.exists():
                backup_dir = Path(tempfile.mkdtemp(prefix=".query_annotations.backup.", dir=str(diagnostics_dir)))
                backup_dir.rmdir()
                os.replace(query_root, backup_dir)
            os.replace(staging_dir, query_root)
        except Exception:
            if staging_dir.exists():
                shutil.rmtree(staging_dir, ignore_errors=True)
            if backup_dir is not None and backup_dir.exists() and not query_root.exists():
                os.replace(backup_dir, query_root)
            raise
        else:
            if backup_dir is not None and backup_dir.exists():
                shutil.rmtree(backup_dir, ignore_errors=True)

        remove_output_artifacts(stale_paths)
        paths["per_query_root"] = str(query_root)
        paths["per_query_manifest"] = str(query_root / "manifest.tsv")
        return paths

    def _validate_summary(self, summary: pd.DataFrame) -> None:
        """Validate the aggregate summary before writing any final files."""
        if "Category" not in summary.columns:
            raise ValueError("partitioned-h2 summary is missing required column: Category")

    def _query_records(self, summary: pd.DataFrame) -> list[dict[str, object]]:
        """Return ordered query records with deterministic safe folder names."""
        records: list[dict[str, object]] = []
        width = max(4, len(str(len(summary))))
        for idx, query_name in enumerate(summary["Category"].astype(str).tolist(), start=1):
            slug = _slugify_query_name(query_name)
            folder = f"{idx:0{width}d}_{slug}"
            records.append(
                {
                    "ordinal": idx,
                    "query_annotation": query_name,
                    "slug": slug,
                    "folder": folder,
                }
            )
        return records

    def _write_staged_query_tree(
        self,
        staging_dir: Path,
        query_records: list[dict[str, object]],
        summary: pd.DataFrame,
        per_query_category_tables: dict[str, pd.DataFrame],
        metadata: dict[str, object],
        per_query_metadata: dict[str, dict[str, object]],
    ) -> list[dict[str, object]]:
        """Populate the staged per-query result tree and return manifest rows."""
        manifest_rows: list[dict[str, object]] = []
        for record in query_records:
            query_name = str(record["query_annotation"])
            folder = str(record["folder"])
            query_dir = staging_dir / folder
            query_dir.mkdir(parents=True, exist_ok=False)
            query_summary = summary.loc[summary["Category"].astype(str) == query_name].reset_index(drop=True)
            summary_rel = f"diagnostics/query_annotations/{folder}/partitioned_h2.tsv"
            full_rel = f"diagnostics/query_annotations/{folder}/partitioned_h2_full.tsv"
            metadata_rel = f"diagnostics/query_annotations/{folder}/metadata.json"
            _atomic_write_dataframe(
                _select_columns(query_summary, PARTITIONED_H2_AGGREGATE_COLUMNS, label="query summary"),
                query_dir / "partitioned_h2.tsv",
            )
            category_table = per_query_category_tables.get(query_name)
            if category_table is None:
                category_table = pd.DataFrame()
            _atomic_write_dataframe(
                _select_columns(category_table, PARTITIONED_H2_FULL_COLUMNS, label="full partitioned-h2 summary"),
                query_dir / "partitioned_h2_full.tsv",
            )
            payload = {
                **metadata,
                **per_query_metadata.get(query_name, {}),
                "schema_version": 1,
                "artifact_type": "partitioned_h2_query_result",
                "files": {"summary": summary_rel, "full": full_rel},
                "ordinal": record["ordinal"],
                "query_annotation": query_name,
                "slug": record["slug"],
                "folder": folder,
            }
            _atomic_write_json(payload, query_dir / "metadata.json")
            manifest_rows.append(
                {
                    "ordinal": record["ordinal"],
                    "query_annotation": query_name,
                    "slug": record["slug"],
                    "folder": folder,
                    "summary_path": summary_rel,
                    "partitioned_h2_full_path": full_rel,
                    "metadata_path": metadata_rel,
                }
            )
        return manifest_rows


@dataclass(frozen=True)
class RgOutputConfig:
    """Directory-oriented output config for genetic-correlation summaries.

    Parameters
    ----------
    output_dir : str or os.PathLike[str]
        Directory that receives ``rg.tsv``, ``rg_full.tsv``, and
        ``h2_per_trait.tsv``.
    overwrite : bool, optional
        If ``True``, replace existing fixed rg outputs and remove stale owned
        siblings after a successful write. Default is ``False``.
    write_per_pair_detail : bool, optional
        If ``True``, also write one subdirectory per tested pair under
        ``diagnostics/pairs/``. Default is ``False``.
    """

    output_dir: str | PathLike[str]
    overwrite: bool = False
    write_per_pair_detail: bool = False

    def __post_init__(self) -> None:
        """Normalize the output directory."""
        object.__setattr__(self, "output_dir", _normalize_required_path(self.output_dir))


class RgDirectoryWriter:
    """Write concise, full, per-trait, and optional per-pair rg outputs.

    The workflow-level logger owns ``diagnostics/rg.log``. This writer owns the public data
    artifacts plus diagnostic metadata in the rg output family:
    ``rg.tsv``, ``rg_full.tsv``, ``h2_per_trait.tsv``, and optional
    ``diagnostics/pairs/``.
    The optional pair tree is staged before replacement so a failed write does
    not expose a partially populated final tree.
    """

    def write(self, result: Any, output_config: RgOutputConfig) -> dict[str, str]:
        """Write an ``RgResultFamily``-like object to a result directory.

        Parameters
        ----------
        result : object
            Object exposing ``rg``, ``rg_full``, ``h2_per_trait``, and
            ``per_pair_metadata`` attributes. The first two tables must contain
            the canonical rg concise and full schemas.
        output_config : RgOutputConfig
            Output directory, overwrite policy, and optional per-pair detail
            setting.

        Returns
        -------
        dict of str to str
            Written path map. Always includes ``rg``, ``rg_full``, and
            ``h2_per_trait``; also includes ``pairs_root`` and
            ``pairs_manifest`` when per-pair detail is enabled.

        Raises
        ------
        ValueError
            If required result columns are missing or the pair metadata length
            does not match ``rg_full`` when detail output is requested.
        FileExistsError
            If an owned rg artifact already exists and overwrite is disabled.
        """
        rg = _select_columns(getattr(result, "rg"), RG_CONCISE_COLUMNS, label="rg summary")
        rg_full = _select_columns(getattr(result, "rg_full"), RG_FULL_COLUMNS, label="full rg summary")
        h2_per_trait = getattr(result, "h2_per_trait")
        per_pair_metadata = list(getattr(result, "per_pair_metadata", []))
        if output_config.write_per_pair_detail and len(per_pair_metadata) != len(rg_full):
            raise ValueError("per_pair_metadata must contain one record per rg_full row when pair detail is enabled.")

        output_dir = ensure_output_directory(output_config.output_dir, label="output directory")
        diagnostics_dir = output_dir / "diagnostics"
        metadata_path = diagnostics_dir / "metadata.json"
        rg_path = output_dir / "rg.tsv"
        full_path = output_dir / "rg_full.tsv"
        h2_path = output_dir / "h2_per_trait.tsv"
        pairs_root = diagnostics_dir / "pairs"
        legacy_metadata_path = output_dir / "metadata.json"
        legacy_pairs_root = output_dir / "pairs"
        produced_paths = [metadata_path, rg_path, full_path, h2_path]
        if output_config.write_per_pair_detail:
            produced_paths.append(pairs_root)
        stale_paths = preflight_output_artifact_family(
            produced_paths,
            [metadata_path, rg_path, full_path, h2_path, pairs_root, legacy_metadata_path, legacy_pairs_root],
            overwrite=output_config.overwrite,
            label="rg output artifact",
        )

        paths = {
            "metadata": str(metadata_path),
            "rg": str(rg_path),
            "rg_full": str(full_path),
            "h2_per_trait": str(h2_path),
        }
        diagnostics_dir.mkdir(parents=True, exist_ok=True)
        if not output_config.write_per_pair_detail:
            _atomic_write_dataframe(rg, rg_path, na_rep="NaN")
            _atomic_write_dataframe(rg_full, full_path, na_rep="NaN")
            _atomic_write_dataframe(h2_per_trait, h2_path, na_rep="NaN")
            _atomic_write_json(_rg_root_metadata(result, include_pairs=False), metadata_path)
            remove_output_artifacts(stale_paths)
            return paths

        pair_records = self._pair_records(rg_full, per_pair_metadata)
        staging_dir = Path(tempfile.mkdtemp(prefix=".pairs.tmp.", dir=str(diagnostics_dir)))
        backup_dir: Path | None = None
        try:
            manifest_rows = self._write_staged_pair_tree(staging_dir, pair_records, rg_full)
            _atomic_write_dataframe(pd.DataFrame(manifest_rows), staging_dir / "manifest.tsv", na_rep="NaN")
            _atomic_write_dataframe(rg, rg_path, na_rep="NaN")
            _atomic_write_dataframe(rg_full, full_path, na_rep="NaN")
            _atomic_write_dataframe(h2_per_trait, h2_path, na_rep="NaN")
            _atomic_write_json(_rg_root_metadata(result, include_pairs=True), metadata_path)
            if output_config.overwrite and pairs_root.exists():
                backup_dir = Path(tempfile.mkdtemp(prefix=".pairs.backup.", dir=str(diagnostics_dir)))
                backup_dir.rmdir()
                os.replace(pairs_root, backup_dir)
            os.replace(staging_dir, pairs_root)
        except Exception:
            if staging_dir.exists():
                shutil.rmtree(staging_dir, ignore_errors=True)
            if backup_dir is not None and backup_dir.exists() and not pairs_root.exists():
                os.replace(backup_dir, pairs_root)
            raise
        else:
            if backup_dir is not None and backup_dir.exists():
                shutil.rmtree(backup_dir, ignore_errors=True)

        remove_output_artifacts(stale_paths)
        paths["pairs_root"] = str(pairs_root)
        paths["pairs_manifest"] = str(pairs_root / "manifest.tsv")
        return paths

    def _pair_records(self, rg_full: pd.DataFrame, per_pair_metadata: list[dict[str, object]]) -> list[dict[str, object]]:
        """Return ordered pair records with deterministic safe folder names."""
        width = max(4, len(str(len(rg_full))))
        records: list[dict[str, object]] = []
        for idx, (_, row) in enumerate(rg_full.iterrows(), start=1):
            metadata = dict(per_pair_metadata[idx - 1])
            trait_1 = str(row["trait_1"])
            trait_2 = str(row["trait_2"])
            slug = f"{_slugify_rg_trait_name(trait_1)}_vs_{_slugify_rg_trait_name(trait_2)}"
            folder = f"{idx:0{width}d}_{slug}"
            records.append(
                {
                    "ordinal": idx,
                    "trait_1": trait_1,
                    "trait_2": trait_2,
                    "slug": slug,
                    "folder": folder,
                    "metadata": metadata,
                }
            )
        return records

    def _write_staged_pair_tree(
        self,
        staging_dir: Path,
        pair_records: list[dict[str, object]],
        rg_full: pd.DataFrame,
    ) -> list[dict[str, object]]:
        """Populate the staged pair result tree and return manifest rows."""
        manifest_rows: list[dict[str, object]] = []
        for record in pair_records:
            folder = str(record["folder"])
            pair_dir = staging_dir / folder
            pair_dir.mkdir(parents=True, exist_ok=False)
            row = rg_full.iloc[[int(record["ordinal"]) - 1]].reset_index(drop=True)
            detail_rel = f"diagnostics/pairs/{folder}/rg_full.tsv"
            metadata_rel = f"diagnostics/pairs/{folder}/metadata.json"
            _atomic_write_dataframe(row, pair_dir / "rg_full.tsv", na_rep="NaN")
            payload = {
                **dict(record["metadata"]),
                "schema_version": 1,
                "artifact_type": "rg_pair_result",
                "files": {"rg_full": detail_rel},
                "ordinal": record["ordinal"],
                "trait_1": record["trait_1"],
                "trait_2": record["trait_2"],
                "slug": record["slug"],
                "folder": folder,
            }
            _atomic_write_json(payload, pair_dir / "metadata.json")
            manifest_rows.append(
                {
                    "ordinal": record["ordinal"],
                    "trait_1": record["trait_1"],
                    "trait_2": record["trait_2"],
                    "slug": record["slug"],
                    "folder": folder,
                    "rg_full_path": detail_rel,
                    "metadata_path": metadata_rel,
                }
            )
        return manifest_rows


def _ldscore_output_family(output_dir: Path) -> list[Path]:
    """Return all fixed data artifacts owned by one LD-score result directory."""
    return [
        output_dir / "metadata.json",
        output_dir / "ldscore.baseline.parquet",
        output_dir / "ldscore.query.parquet",
    ]


def _result_metadata(
    metadata: dict[str, object] | None,
    *,
    artifact_type: str,
    files: dict[str, str],
) -> dict[str, object]:
    """Return root result metadata with the canonical discriminator fields."""
    payload = dict(metadata or {})
    payload["schema_version"] = 1
    payload["artifact_type"] = artifact_type
    payload["files"] = dict(files)
    payload.pop("format", None)
    return payload


def _rg_root_metadata(result: Any, *, include_pairs: bool) -> dict[str, object]:
    """Build diagnostic metadata for an rg result directory from result tables."""
    rg_full = getattr(result, "rg_full", pd.DataFrame())
    h2_per_trait = getattr(result, "h2_per_trait", pd.DataFrame())
    files = {
        "rg": "rg.tsv",
        "rg_full": "rg_full.tsv",
        "h2_per_trait": "h2_per_trait.tsv",
    }
    if include_pairs:
        files["pairs"] = "diagnostics/pairs"
    trait_names: list[str] = []
    if "trait_name" in h2_per_trait.columns:
        trait_names = [str(value) for value in h2_per_trait["trait_name"].dropna().tolist()]
    pair_kind = None
    if "pair_kind" in rg_full.columns:
        values = [str(value) for value in rg_full["pair_kind"].dropna().unique().tolist()]
        pair_kind = values[0] if len(values) == 1 else ("mixed" if values else None)
    payload = _result_metadata(
        {
            "trait_names": trait_names,
            "pair_kind": pair_kind,
            "n_pairs": int(len(rg_full)),
            "n_traits": int(len(trait_names)),
        },
        artifact_type="rg_result",
        files=files,
    )
    ldscore_dir = getattr(result, "ldscore_dir", None)
    if ldscore_dir is not None:
        payload["ldscore_dir"] = ldscore_dir
    sumstats_files = getattr(result, "sumstats_files", None)
    if sumstats_files is not None:
        payload["sumstats_files"] = list(sumstats_files)
    return payload


def _slugify_query_name(value: str) -> str:
    """Return a filesystem-safe slug for a query annotation name."""
    normalized = unicodedata.normalize("NFKD", str(value)).encode("ascii", "ignore").decode("ascii")
    slug = re.sub(r"[^a-z0-9._-]+", "_", normalized.lower()).strip("._-")
    return slug or "annotation"


def _slugify_rg_trait_name(value: str) -> str:
    """Return a filesystem-safe slug for an rg trait name."""
    normalized = unicodedata.normalize("NFKD", str(value)).encode("ascii", "ignore").decode("ascii")
    slug = re.sub(r"[^a-z0-9._-]+", "_", normalized.lower()).strip("._-")
    return slug or "trait"


def _select_columns(df: pd.DataFrame, columns: list[str], *, label: str) -> pd.DataFrame:
    """Return ``df`` with required public output columns in canonical order."""
    if df.empty and len(df.columns) == 0:
        return df
    missing = [column for column in columns if column not in df.columns]
    if missing:
        raise ValueError(f"{label} is missing required columns: {missing}")
    return df.loc[:, columns]


def _atomic_write_dataframe(df: pd.DataFrame, path: Path, *, na_rep: str | None = None) -> None:
    """Write a dataframe through a temporary sibling file, then replace."""
    fd, tmp_name = tempfile.mkstemp(prefix=f".{path.name}.", suffix=".tmp", dir=str(path.parent))
    os.close(fd)
    tmp_path = Path(tmp_name)
    try:
        kwargs = {"sep": "\t", "index": False}
        if na_rep is not None:
            kwargs["na_rep"] = na_rep
        df.to_csv(tmp_path, **kwargs)
        os.replace(tmp_path, path)
    except Exception:
        tmp_path.unlink(missing_ok=True)
        raise


def _atomic_write_json(payload: dict[str, object], path: Path) -> None:
    """Write JSON through a temporary sibling file, then replace."""
    fd, tmp_name = tempfile.mkstemp(prefix=f".{path.name}.", suffix=".tmp", dir=str(path.parent))
    os.close(fd)
    tmp_path = Path(tmp_name)
    try:
        tmp_path.write_text(json.dumps(_to_serializable(payload), indent=2, sort_keys=True), encoding="utf-8")
        os.replace(tmp_path, path)
    except Exception:
        tmp_path.unlink(missing_ok=True)
        raise


def _to_serializable(value: Any) -> Any:
    """Recursively convert workflow objects into JSON-serializable structures."""
    if is_dataclass(value):
        return _to_serializable(asdict(value))
    if isinstance(value, pd.DataFrame):
        return value.to_dict(orient="records")
    if isinstance(value, pd.Series):
        return value.to_list()
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, dict):
        return {str(key): _to_serializable(val) for key, val in value.items()}
    if isinstance(value, (list, tuple)):
        return [_to_serializable(item) for item in value]
    return value
