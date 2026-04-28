"""Output writers for the refactored LDSC workflows.

The public LD-score workflow writes one canonical result directory containing
``manifest.json``, ``baseline.parquet``, and optional ``query.parquet``. Run
identity comes from the chosen directory name; output filenames inside that
directory are fixed.
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass, is_dataclass
from os import PathLike
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from .config import _normalize_required_path
from .path_resolution import ensure_output_directory


LDSCORE_RESULT_FORMAT = "ldsc.ldscore_result.v1"


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
        ``output_dir``. Default is ``False``.
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
    """Write the canonical LD-score result directory."""

    def write(self, result: Any, output_config: LDScoreOutputConfig) -> dict[str, str]:
        """Write ``manifest.json``, ``baseline.parquet``, and optional ``query.parquet``."""
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
        existing = [path for path in paths.values() if path.exists()]
        if existing and not output_config.overwrite:
            raise FileExistsError(_format_existing_artifact_message(existing))

        compression = None if output_config.parquet_compression in {None, "none"} else output_config.parquet_compression
        baseline_table.to_parquet(paths["baseline"], index=False, compression=compression)
        if query_table is not None:
            query_table.to_parquet(paths["query"], index=False, compression=compression)
        manifest = self.build_manifest(result, files={name: path.name for name, path in paths.items() if name != "manifest"})
        paths["manifest"].write_text(json.dumps(_to_serializable(manifest), indent=2, sort_keys=True), encoding="utf-8")
        return {name: str(path) for name, path in paths.items()}

    def build_manifest(self, result: Any, files: dict[str, str]) -> dict[str, Any]:
        """Build the JSON manifest payload for one LD-score result."""
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
            "config_snapshot": config_snapshot,
            "n_baseline_rows": int(len(baseline_table)),
            "n_query_rows": 0 if query_table is None else int(len(query_table)),
        }

    def _validate_tables(self, result: Any) -> None:
        """Validate baseline/query table shape before any files are written."""
        baseline_table = getattr(result, "baseline_table")
        query_table = getattr(result, "query_table", None)
        baseline_columns = list(getattr(result, "baseline_columns", []))
        query_columns = list(getattr(result, "query_columns", []))
        required_baseline = ["CHR", "SNP", "BP", "regr_weight", *baseline_columns]
        missing = [column for column in required_baseline if column not in baseline_table.columns]
        if missing:
            raise ValueError(f"baseline_table is missing required columns: {missing}")
        if query_columns and query_table is None:
            raise ValueError("query_table is required when query_columns are present.")
        if not query_columns and query_table is not None:
            raise ValueError("query_table was provided but query_columns is empty.")
        if query_table is None:
            return
        required_query = ["CHR", "SNP", "BP", *query_columns]
        missing = [column for column in required_query if column not in query_table.columns]
        if missing:
            raise ValueError(f"query_table is missing required columns: {missing}")
        baseline_keys = baseline_table.loc[:, ["CHR", "SNP", "BP"]].reset_index(drop=True)
        query_keys = query_table.loc[:, ["CHR", "SNP", "BP"]].reset_index(drop=True)
        if not baseline_keys.equals(query_keys):
            raise ValueError("query rows must match baseline rows on CHR/SNP/BP.")


def _format_existing_artifact_message(paths: list[Path]) -> str:
    """Build the shared overwrite-guard message for one or more destinations."""
    preview = ", ".join(str(path) for path in paths[:3])
    if len(paths) > 3:
        preview = f"{preview}, ... ({len(paths) - 3} more)"
    label = "artifact" if len(paths) == 1 else "artifacts"
    return f"Refusing to overwrite existing {label}: {preview}. Pass overwrite=True to replace them."


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
