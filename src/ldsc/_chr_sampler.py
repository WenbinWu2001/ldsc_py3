"""Sample CHR/POS frames from chromosome-suite input tokens."""

from __future__ import annotations

import os
from collections.abc import Sequence
from pathlib import Path

import pandas as pd

from .column_inference import infer_chr_pos_columns
from .path_resolution import substitute_chromosome


def sample_frame_from_chr_pattern(
    tokens: Sequence[str],
    *,
    chromosomes: Sequence[str] | None = None,
    nrows: int = 5000,
    context: str,
) -> tuple[pd.DataFrame, str]:
    """Read a small CHR/POS sample from the first resolvable ``@`` token."""
    candidates = _candidate_chromosomes(chromosomes)
    saw_suite_token = False
    attempted: list[str] = []
    for token in tokens:
        if "@" not in token:
            continue
        saw_suite_token = True
        for chrom in candidates:
            path = substitute_chromosome(token, chrom)
            attempted.append(path)
            if not os.path.exists(path):
                continue
            frame = _read_head(path, nrows=nrows)
            chr_col, pos_col = infer_chr_pos_columns(frame.columns, context=path)
            return frame.loc[:, [chr_col, pos_col]].rename(columns={chr_col: "CHR", pos_col: "POS"}), path

    if not saw_suite_token:
        raise ValueError(
            f"Cannot infer --genome-build for {context}: at least one input path token must contain '@' "
            "so a chromosome-specific sample can be read."
        )
    attempted_text = ", ".join(attempted[:5])
    if len(attempted) > 5:
        attempted_text += ", ..."
    raise FileNotFoundError(
        f"Cannot infer --genome-build for {context}: no chromosome sample file exists for the provided "
        f"'@' token(s). Tried: {attempted_text}"
    )


def _candidate_chromosomes(chromosomes: Sequence[str] | None) -> list[str]:
    values = ["1"]
    if chromosomes is not None:
        values.extend(str(chrom) for chrom in chromosomes)
    seen: set[str] = set()
    candidates: list[str] = []
    for value in values:
        if value not in seen:
            seen.add(value)
            candidates.append(value)
    return candidates


def _read_head(path: str, *, nrows: int) -> pd.DataFrame:
    if Path(path).suffix == ".parquet":
        import pyarrow.parquet as pq

        return pq.read_table(path).to_pandas().head(nrows)
    return pd.read_csv(path, sep="\t", compression="infer", nrows=nrows)
