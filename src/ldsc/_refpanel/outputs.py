"""Writer for the ``ldsc query-r2`` result directory.

``QueryR2DirectoryWriter`` declares the files a pair query owns and writes
the queried pairs with their metadata.
"""

from __future__ import annotations

import json
import os
from dataclasses import dataclass
from os import PathLike

import pandas as pd

from ..config import normalize_required_path
from .._runtime.artifacts import (
    ArtifactFamily,
    atomic_write_dataframe,
    atomic_write_json,
    declare_artifacts,
    result_metadata,
)
from .._runtime.paths import ensure_output_directory


@dataclass(frozen=True)
class QueryR2OutputConfig:
    """Directory-oriented output config for query-r2 pair results.

    Parameters
    ----------
    output_dir : str or os.PathLike[str]
        Directory that receives ``query_r2.tsv`` and ``diagnostics/metadata.json``.
    overwrite : bool, optional
        If ``True``, replace existing query-r2 outputs. Default is ``False``.
    """

    output_dir: str | PathLike[str]
    overwrite: bool = False

    def __post_init__(self) -> None:
        """Normalize the output directory."""
        object.__setattr__(self, "output_dir", normalize_required_path(self.output_dir))


class QueryR2DirectoryWriter:
    """Write query-r2 results and diagnostic provenance.

    The workflow supplies input columns plus ``r2``, nullable ``sign_r``,
    ``r``, and ``status``. The TSV preserves those columns and writes missing
    values as ``NaN``; it does not add compatibility aliases.
    """

    @staticmethod
    def artifact_family(output_dir) -> ArtifactFamily:
        """Declare the pair-query result table and provenance."""
        return declare_artifacts(output_dir, {
            "result": ("query_r2.tsv", True),
            "metadata": ("diagnostics/metadata.json", True),
        }, label="query-r2 output artifact")

    def write(
        self,
        result: pd.DataFrame,
        output_config: QueryR2OutputConfig,
        *,
        metadata: dict[str, object],
    ) -> dict[str, str]:
        """Write ``query_r2.tsv`` and ``diagnostics/metadata.json`` to a result directory.

        Existing fixed query-r2 artifacts are checked before any output file is
        written. Replacement requires ``output_config.overwrite=True``.
        """
        output_dir = ensure_output_directory(output_config.output_dir, label="output directory")
        family = self.artifact_family(output_dir)
        result_path = family.paths["result"]
        metadata_path = family.paths["metadata"]
        diagnostics_dir = metadata_path.parent
        family.preflight(overwrite=output_config.overwrite)
        diagnostics_dir.mkdir(parents=True, exist_ok=True)
        atomic_write_dataframe(result, result_path, na_rep="NaN")
        atomic_write_json(
            result_metadata(metadata, artifact_type="query_r2_result", files=family.metadata_files()),
            metadata_path,
        )
        return {name: str(path) for name, path in family.paths.items()}
