"""Small file operations shared by derived result workflows.

Workflow-specific metadata validation and scientific dispatch stay with callers.
"""

import json
import os
from pathlib import Path
import tempfile

from .errors import LDSCInputError


def declared_result_file(
    result_dir: Path, metadata: dict, key: str, *, context: str, allow_missing: bool = False
) -> Path:
    """Resolve a declared input inside the canonical result root.

    ``allow_missing`` permits an absent target, for example an optional rg
    heritability table. It does not relax metadata or containment checks:
    absolute paths, escaping symlinks, and existing non-files are rejected.
    """
    files = metadata.get("files")
    if not isinstance(files, dict) or not isinstance(files.get(key), str) or not files[key].strip():
        raise LDSCInputError(f"{context} must declare files.{key}.")
    token = Path(files[key])
    if token.is_absolute():
        raise LDSCInputError(f"{context} files.{key} must be relative to the result directory.")
    path = (result_dir / token).resolve()
    try:
        path.relative_to(result_dir)
    except ValueError as exc:
        raise LDSCInputError(f"{context} files.{key} escapes the result directory.") from exc
    if allow_missing and not path.exists():
        return path
    if not path.is_file():
        raise LDSCInputError(f"{context} declares files.{key}='{files[key]}', but it is missing.")
    return path


def atomic_write_json(payload: dict[str, object], path: Path) -> None:
    """Publish JSON through a temporary sibling, cleaning it on failure."""
    file_descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=str(path.parent)
    )
    os.close(file_descriptor)
    temporary_path = Path(temporary_name)
    try:
        temporary_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        os.replace(temporary_path, path)
    except Exception:
        temporary_path.unlink(missing_ok=True)
        raise
