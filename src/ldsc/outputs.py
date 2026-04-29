"""Output writers for the refactored LDSC workflows.

The public LD-score workflow writes one canonical result directory containing
``manifest.json``, ``baseline.parquet``, and optional ``query.parquet``. Run
identity comes from the chosen directory name; output filenames inside that
directory are fixed. The parquet payloads are written with one row group per
chromosome and matching manifest metadata so downstream readers can load a
single chromosome without scanning the whole table. The writer creates missing
directories, reuses existing directories, and refuses existing canonical files
unless the caller explicitly sets ``overwrite=True``.
"""

from __future__ import annotations

import json
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
from .path_resolution import ensure_output_directory, ensure_output_paths_available


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
DEFAULT_COUNT_CONFIG = {
    "common_reference_snp_maf_min": 0.05,
    "common_reference_snp_maf_operator": ">=",
}


@dataclass(frozen=True)
class LDScoreOutputConfig:
    """Directory-oriented output config for canonical LD-score results.

    Parameters
    ----------
    output_dir : str or os.PathLike[str]
        Directory that receives ``manifest.json``, ``baseline.parquet``, and
        optional ``query.parquet``.
    overwrite : bool, optional
        If ``True``, replace existing canonical LD-score files in
        ``output_dir``. If ``False``, existing canonical files raise
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

    The writer owns the fixed files ``manifest.json``, ``baseline.parquet``,
    and optional ``query.parquet``. Parquet files are flat files for backward
    compatibility, but their internal row groups are chromosome-aligned and
    described in the manifest.
    """

    def write(self, result: Any, output_config: LDScoreOutputConfig) -> dict[str, str]:
        """Write ``manifest.json``, ``baseline.parquet``, and optional ``query.parquet``.

        Existing canonical files are checked before any output file is written.
        Replacement requires ``output_config.overwrite=True``; unrelated files
        in the directory are ignored. The returned paths map includes
        ``"baseline"``, ``"manifest"``, and ``"query"`` when query annotations
        were supplied.
        """
        output_dir = ensure_output_directory(output_config.output_dir, label="LD-score output directory")
        baseline_table = getattr(result, "baseline_table", None)
        query_table = getattr(result, "query_table", None)
        if baseline_table is None:
            raise ValueError("LDScoreResult is missing baseline_table.")
        self._validate_tables(result)

        paths = {
            "manifest": output_dir / "manifest.json",
            "baseline": output_dir / "baseline.parquet",
        }
        if query_table is not None:
            paths["query"] = output_dir / "query.parquet"
        ensure_output_paths_available(
            paths.values(),
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
        required_baseline = ["CHR", "POS", "SNP", "regr_weight", *baseline_columns]
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
