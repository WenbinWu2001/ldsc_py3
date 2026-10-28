"""Shared declarations and atomic writes for result directories.

Each workflow writer declares the files it owns as an ``ArtifactFamily``,
which output preflight, final writes, metadata file entries and stale-file
cleanup all share. The atomic writers publish tables and JSON through a
temporary sibling file, so a reader never sees a partly written file.
"""

from __future__ import annotations

import json
import os
import tempfile
from dataclasses import asdict, dataclass, is_dataclass
from os import PathLike
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd

from .paths import preflight_output_artifact_family


@dataclass(frozen=True)
class ArtifactFamily:
    """Resolved writer declaration shared by collision checks and publication.

    Parameters
    ----------
    output_dir : pathlib.Path
        Result root used to make metadata file entries relative.
    paths : dict of str to pathlib.Path
        Outputs selected for this write, keyed by their artifact names.
    owned_paths : tuple of pathlib.Path
        All owned paths, including conditional siblings and derived roots
        eligible for cleanup after a successful replacement.
    label : str
        Workflow-specific description used in collision errors.

    Notes
    -----
    Construction only describes paths. Callers retain their existing write
    order and remove returned stale paths only after successful publication.
    """

    output_dir: Path
    paths: dict[str, Path]
    owned_paths: tuple[Path, ...]
    label: str

    def preflight(
        self, *, overwrite: bool, additional_paths: Iterable[str | PathLike[str]] = (),
    ) -> list[Path]:
        """Check ownership, including workflow-only paths such as the run log.

        ``overwrite=False`` raises ``FileExistsError`` for any existing owned
        path. With overwrite enabled, return existing owned paths that this
        declaration does not produce. ``additional_paths`` are protected
        workflow outputs; their writer controls their cleanup separately.

        Early callers discard the return value because production is not yet
        known. Writers call again with the final declaration and use that
        stale list after writing. Additional paths never enter result metadata.
        """
        additional_paths = tuple(additional_paths)
        return preflight_output_artifact_family(
            [*self.paths.values(), *additional_paths],
            [*self.owned_paths, *additional_paths],
            overwrite=overwrite, label=self.label,
        )

    def metadata_files(self, *, exclude: Iterable[str] = ()) -> dict[str, str]:
        """Return selected paths relative to the result root, excluding metadata itself."""
        omitted = {"metadata", *exclude}
        return {
            name: str(path.relative_to(self.output_dir))
            for name, path in self.paths.items() if name not in omitted
        }


def _declare_artifacts(output_dir, entries, *, label, extra_owned=()) -> ArtifactFamily:
    """Resolve each named (relative path, produced) entry exactly once."""
    output_dir = Path(output_dir)
    owned = {name: output_dir / relative for name, (relative, _) in entries.items()}
    return ArtifactFamily(
        output_dir,
        {name: owned[name] for name, (_, produced) in entries.items() if produced},
        (*owned.values(), *extra_owned), label,
    )


def _result_metadata(
    metadata: dict[str, object] | None,
    *,
    artifact_type: str,
    files: dict[str, str],
) -> dict[str, object]:
    """Return root result metadata with the canonical discriminator fields."""
    payload = dict(metadata or {})
    payload["artifact_type"] = artifact_type
    payload["files"] = dict(files)
    return payload


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
    """Publish serializable JSON through a temporary sibling, cleaning it on failure."""
    fd, tmp_name = tempfile.mkstemp(prefix=f".{path.name}.", suffix=".tmp", dir=str(path.parent))
    os.close(fd)
    tmp_path = Path(tmp_name)
    try:
        tmp_path.write_text(json.dumps(_to_serializable(payload), indent=2, sort_keys=True) + "\n", encoding="utf-8")
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
