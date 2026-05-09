"""Output writers for the refactored LDSC workflows.

The public LD-score workflow writes one canonical result directory containing
``manifest.json``, ``ldscore.baseline.parquet``, and optional
``ldscore.query.parquet``. Run identity comes from the chosen directory name;
output filenames inside that directory are fixed. The parquet payloads are
written with one row group per chromosome and matching manifest metadata so
downstream readers can load a single chromosome without scanning the whole
table. The writer creates missing directories, reuses existing directories, and
refuses existing canonical family files unless the caller explicitly sets
``overwrite=True``; successful overwrites remove stale owned siblings that the
current result did not produce.

Partitioned-h2 regression summaries use the same directory-oriented output
policy. ``PartitionedH2DirectoryWriter`` always writes the compact aggregate
``partitioned_h2.tsv`` and can optionally stage a per-query
``query_annotations`` tree with ``manifest.tsv``, one-row
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
    """Write one parquet row group per chromosome and return manifest metadata.

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


LDSCORE_RESULT_FORMAT = "ldsc.ldscore_result.v1"
REGRESSION_LD_SCORE_COLUMN = "regression_ld_scores"
DEFAULT_COUNT_CONFIG = {
    "common_reference_snp_maf_min": 0.05,
    "common_reference_snp_maf_operator": ">=",
}
PARTITIONED_H2_RESULT_FORMAT = "ldsc.partitioned_h2_result.v1"
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


@dataclass(frozen=True)
class LDScoreOutputConfig:
    """Directory-oriented output config for canonical LD-score results.

    Parameters
    ----------
    output_dir : str or os.PathLike[str]
        Directory that receives ``manifest.json``,
        ``ldscore.baseline.parquet``, and optional ``ldscore.query.parquet``.
    overwrite : bool, optional
        If ``True``, replace existing canonical LD-score files in
        ``output_dir`` and remove stale owned siblings after a successful
        write. If ``False``, any existing canonical family file raises
        ``FileExistsError`` before any parquet or manifest file is written.
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

    The writer owns the fixed files ``manifest.json``,
    ``ldscore.baseline.parquet``, and optional ``ldscore.query.parquet``.
    Parquet files are flat files for backward compatibility, but their internal
    row groups are chromosome-aligned and described in the manifest.
    """

    def write(self, result: Any, output_config: LDScoreOutputConfig) -> dict[str, str]:
        """Write ``manifest.json`` and canonical LD-score parquet files.

        Existing canonical family files are checked before any output file is
        written. Replacement requires ``output_config.overwrite=True``;
        unrelated files in the directory are ignored. The returned paths map
        includes ``"baseline"``, ``"manifest"``, and ``"query"`` when query
        annotations were supplied.
        """
        output_dir = ensure_output_directory(output_config.output_dir, label="LD-score output directory")
        baseline_table = getattr(result, "baseline_table", None)
        query_table = getattr(result, "query_table", None)
        if baseline_table is None:
            raise ValueError("LDScoreResult is missing baseline_table.")
        self._validate_tables(result)

        paths = {
            "manifest": output_dir / "manifest.json",
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
        manifest = self.build_manifest(
            result,
            files={name: path.name for name, path in paths.items() if name != "manifest"},
            baseline_rg=baseline_rg,
            query_rg=query_rg,
        )
        paths["manifest"].write_text(json.dumps(_to_serializable(manifest), indent=2, sort_keys=True), encoding="utf-8")
        remove_output_artifacts(stale_paths)
        return {name: str(path) for name, path in paths.items()}

    def build_manifest(
        self,
        result: Any,
        files: dict[str, str],
        baseline_rg: list[dict] | None = None,
        query_rg: list[dict] | None = None,
    ) -> dict[str, Any]:
        """Build the JSON manifest payload for one LD-score result.

        The manifest always includes ``count_config`` so downstream regression
        code can report the common-SNP count threshold even when MAF metadata is
        unavailable and per-column common counts are omitted. Row-group metadata
        records the chromosome, row-group index, row offset, and row count for
        each chromosome-aligned parquet row group.
        """
        baseline_table = getattr(result, "baseline_table")
        query_table = getattr(result, "query_table", None)
        config_snapshot = getattr(result, "config_snapshot", None)
        chromosomes = baseline_table["CHR"].astype(str).drop_duplicates().tolist()
        return {
            "format": LDSCORE_RESULT_FORMAT,
            "files": dict(files),
            "snp_identifier": getattr(config_snapshot, "snp_identifier", None),
            "genome_build": getattr(config_snapshot, "genome_build", None),
            "chromosomes": chromosomes,
            "baseline_columns": list(getattr(result, "baseline_columns", [])),
            "query_columns": list(getattr(result, "query_columns", [])),
            "counts": list(getattr(result, "count_records", [])),
            "count_config": dict(getattr(result, "count_config", None) or DEFAULT_COUNT_CONFIG),
            "config_snapshot": config_snapshot,
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
        required_baseline = ["CHR", "POS", "SNP", REGRESSION_LD_SCORE_COLUMN, *baseline_columns]
        missing = [column for column in required_baseline if column not in baseline_table.columns]
        if missing:
            raise ValueError(f"baseline_table is missing required columns: {missing}")
        if query_columns and query_table is None:
            raise ValueError("query_table is required when query_columns are present.")
        if not query_columns and query_table is not None:
            raise ValueError("query_table was provided but query_columns is empty.")
        if query_table is None:
            return
        required_query = ["CHR", "POS", "SNP", *query_columns]
        missing = [column for column in required_query if column not in query_table.columns]
        if missing:
            raise ValueError(f"query_table is missing required columns: {missing}")
        assert_same_snp_rows(
            baseline_table,
            query_table,
            context="query rows must match baseline rows on CHR/SNP/POS",
        )


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
        the optional ``query_annotations`` per-query result tree.
    overwrite : bool, optional
        If ``True``, replace existing fixed partitioned-h2 outputs and remove
        stale owned siblings after a successful write. Default is ``False``.
    write_per_query_results : bool, optional
        If ``True``, also write one subdirectory per query annotation under
        ``query_annotations``. Default is ``False``.
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
    while optional query folders contain the same one-row summary plus a full
    baseline-plus-query ``partitioned_h2_full.tsv`` table. The per-query tree
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
        query_root = output_dir / "query_annotations"
        produced_paths = [summary_path]
        if output_config.write_per_query_results:
            produced_paths.append(query_root)
        stale_paths = preflight_output_artifact_family(
            produced_paths,
            [summary_path, query_root],
            overwrite=output_config.overwrite,
            label="partitioned-h2 output artifact",
        )

        paths = {"summary": str(summary_path)}
        if not output_config.write_per_query_results:
            _atomic_write_dataframe(
                _select_columns(summary, PARTITIONED_H2_AGGREGATE_COLUMNS, label="partitioned-h2 summary"),
                summary_path,
            )
            remove_output_artifacts(stale_paths)
            return paths

        query_records = self._query_records(summary)
        staging_dir = Path(
            tempfile.mkdtemp(prefix=".query_annotations.tmp.", dir=str(output_dir))
        )
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
            if output_config.overwrite and query_root.exists():
                backup_dir = Path(tempfile.mkdtemp(prefix=".query_annotations.backup.", dir=str(output_dir)))
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
            summary_rel = f"query_annotations/{folder}/partitioned_h2.tsv"
            full_rel = f"query_annotations/{folder}/partitioned_h2_full.tsv"
            metadata_rel = f"query_annotations/{folder}/metadata.json"
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
                "format": PARTITIONED_H2_RESULT_FORMAT,
                **metadata,
                **per_query_metadata.get(query_name, {}),
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


def _ldscore_output_family(output_dir: Path) -> list[Path]:
    """Return all fixed data artifacts owned by one LD-score result directory."""
    return [
        output_dir / "manifest.json",
        output_dir / "ldscore.baseline.parquet",
        output_dir / "ldscore.query.parquet",
    ]


def _slugify_query_name(value: str) -> str:
    """Return a filesystem-safe slug for a query annotation name."""
    normalized = unicodedata.normalize("NFKD", str(value)).encode("ascii", "ignore").decode("ascii")
    slug = re.sub(r"[^a-z0-9._-]+", "_", normalized.lower()).strip("._-")
    return slug or "annotation"


def _select_columns(df: pd.DataFrame, columns: list[str], *, label: str) -> pd.DataFrame:
    """Return ``df`` with required public output columns in canonical order."""
    if df.empty and len(df.columns) == 0:
        return df
    missing = [column for column in columns if column not in df.columns]
    if missing:
        raise ValueError(f"{label} is missing required columns: {missing}")
    return df.loc[:, columns]


def _atomic_write_dataframe(df: pd.DataFrame, path: Path) -> None:
    """Write a dataframe through a temporary sibling file, then replace."""
    fd, tmp_name = tempfile.mkstemp(prefix=f".{path.name}.", suffix=".tmp", dir=str(path.parent))
    os.close(fd)
    tmp_path = Path(tmp_name)
    try:
        df.to_csv(tmp_path, sep="\t", index=False)
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
