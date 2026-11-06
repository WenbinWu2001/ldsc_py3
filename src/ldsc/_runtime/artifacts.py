"""Shared declarations and atomic writes for result directories.

Each workflow writer declares the files it owns as an ``ArtifactFamily``,
which output preflight, final writes, metadata file entries and stale-file
cleanup all share; ``declared_member`` and ``declared_result_file`` read those
``files`` entries back with one containment rule. ``atomic_replace`` publishes a file through a temporary
sibling, so a reader never sees a partly written file; the table and JSON
writers here, the failure marker and the legacy summary-statistics audit use
it.
"""

from __future__ import annotations

import json
import os
import tempfile
from dataclasses import asdict, dataclass, is_dataclass
from os import PathLike
from pathlib import Path
from typing import Any, Callable, Iterable

import numpy as np
import pandas as pd

from ..errors import LDSCInputError
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


def declared_member(root: Path, value, *, external_links: bool) -> Path | None:
    """Return the file a metadata ``files`` entry names inside ``root``, or ``None`` when it names none.

    The one reader of ``files`` entries: an entry names a file inside
    ``root`` when it is a nonblank relative path without ``..`` parts. With
    ``external_links`` the file, or a folder on its path, may be a symbolic
    link to anywhere, as LD-score directories allow for large tables on
    shared storage. Without it the resolved file must also lie inside the
    resolved ``root``, as derived-result readers and the genome-build check
    require. Existence is not checked; callers word their own errors.
    """
    token = Path(value) if isinstance(value, str) and value.strip() else None
    if token is None or token.is_absolute() or ".." in token.parts:
        return None
    path = Path(root) / token
    if not external_links and not path.resolve().is_relative_to(Path(root).resolve()):
        return None
    return path


def declared_result_file(
    result_dir: Path, metadata: dict, key: str, *, context: str, allow_missing: bool = False
) -> Path:
    """Return the file ``metadata["files"][key]`` names inside a result directory, by :func:`declared_member`.

    ``allow_missing`` permits an absent target, for example an optional ``rg``
    heritability table. It does not relax the manifest checks: a missing,
    blank, absolute or ``..`` entry, a symbolic link that leaves the
    directory, and an existing non-file raise ``LDSCInputError``.
    """
    files = metadata.get("files")
    value = files.get(key) if isinstance(files, dict) else None
    if not isinstance(value, str) or not value.strip():
        raise LDSCInputError(
            f"{context} cannot locate the saved {key} artifact: metadata has no usable files.{key} entry. "
            f"Restore the complete result directory '{result_dir}' with metadata and data from the same run, "
            "or regenerate the upstream result."
        )
    path = declared_member(result_dir, value, external_links=False)
    if Path(value).is_absolute():
        raise LDSCInputError(
            f"{context} files.{key} must be relative to the result directory; got '{value}'. "
            "Restore the original metadata and its matching data files from the same run."
        )
    if path is None:
        raise LDSCInputError(
            f"{context} files.{key} escapes the result directory '{result_dir}': '{value}'. "
            "Restore a complete result with its referenced files inside that directory; "
            "do not combine metadata and data from different runs."
        )
    if allow_missing and not path.exists():
        return path
    if not path.is_file():
        raise LDSCInputError(
            f"{context} declares files.{key}='{value}', but '{path}' is missing or is not a regular file. "
            "Restore this file and its matching metadata from the same run, or regenerate the upstream result."
        )
    return path


def declare_artifacts(output_dir, entries, *, label, extra_owned=()) -> ArtifactFamily:
    """Resolve each named (relative path, produced) entry exactly once."""
    output_dir = Path(output_dir)
    owned = {name: output_dir / relative for name, (relative, _) in entries.items()}
    return ArtifactFamily(
        output_dir,
        {name: owned[name] for name, (_, produced) in entries.items() if produced},
        (*owned.values(), *extra_owned), label,
    )


def result_metadata(
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


def atomic_replace(path: Path, write: Callable[[Path], None]) -> None:
    """Publish ``path`` atomically: ``write`` fills a temporary sibling, which then replaces ``path``.

    The temporary file ``.<name>.<random>.tmp`` lies in ``path``'s directory,
    so the replacement is one rename on the same file system. If ``write`` or
    the rename fails, or the run is interrupted, the temporary file is removed
    and ``path`` keeps its previous content.
    """
    fd, tmp_name = tempfile.mkstemp(prefix=f".{path.name}.", suffix=".tmp", dir=str(path.parent))
    os.close(fd)
    tmp_path = Path(tmp_name)
    try:
        write(tmp_path)
        os.replace(tmp_path, path)
    except BaseException:
        tmp_path.unlink(missing_ok=True)
        raise


def atomic_write_dataframe(df: pd.DataFrame, path: Path, *, na_rep: str | None = None) -> None:
    """Write a dataframe as a tab-separated table through :func:`atomic_replace`."""
    kwargs = {"sep": "\t", "index": False}
    if na_rep is not None:
        kwargs["na_rep"] = na_rep
    atomic_replace(path, lambda tmp_path: df.to_csv(tmp_path, **kwargs))


def atomic_write_json(payload: dict[str, object], path: Path) -> None:
    """Publish serializable JSON through :func:`atomic_replace`."""
    text = json.dumps(to_serializable(payload), indent=2, sort_keys=True) + "\n"
    atomic_replace(path, lambda tmp_path: tmp_path.write_text(text, encoding="utf-8"))


def to_serializable(value: Any) -> Any:
    """Recursively convert workflow objects into JSON-serializable structures."""
    if is_dataclass(value):
        return to_serializable(asdict(value))
    if isinstance(value, pd.DataFrame):
        return value.to_dict(orient="records")
    if isinstance(value, pd.Series):
        return value.to_list()
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, dict):
        return {str(key): to_serializable(val) for key, val in value.items()}
    if isinstance(value, (list, tuple)):
        return [to_serializable(item) for item in value]
    return value
