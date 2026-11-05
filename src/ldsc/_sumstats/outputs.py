"""Writer for the ``munge-sumstats`` output directory.

``curated_sumstats_frame`` puts munged rows into the one form that the
returned ``SumstatsTable`` and both data files share: fixed columns,
canonical coordinates and genomic row order.
``SumstatsDirectoryWriter`` declares the files a munging run owns and writes
that frame as self-describing Parquet (one row group per chromosome,
identity footer from ``ldsc._sumstats.table``) and/or the legacy-compatible
``.sumstats.gz`` text, plus the dropped-SNP table.
Every file is written into a staging directory inside the output directory
and moved into place once all of them are complete, so a reader never sees a
partly written file. The column contract of both formats is in
``docs/reference/formats/column-schema.md``.
"""

from __future__ import annotations

from dataclasses import dataclass
import os
from pathlib import Path
import re
from tempfile import TemporaryDirectory
from typing import Any

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

from .._core.chromosomes import chrom_sort_key, normalize_chromosome, normalize_chromosome_series
from .._core.coordinates import coordinate_missing_mask, positive_int_position_series
from .._runtime.artifacts import ArtifactFamily, declare_artifacts
from .._runtime.paths import ensure_output_directory, remove_output_artifacts
from ..errors import LDSCInputError, LDSCUsageError


# File names inside a munging output directory that the workflow builds paths from.
SUMSTATS_LOG = "diagnostics/sumstats.log"
DROPPED_SNPS_FILE = "diagnostics/dropped_snps/dropped.tsv.gz"
INPUT_ISSUES_FILE = "diagnostics/input_issues.tsv"
SUMSTATS_PARQUET_COMPRESSION = "snappy"
SUMSTATS_OUTPUT_FORMATS = {"parquet", "tsv.gz", "both"}
_DATA_FILES = ("tsv.gz", "parquet")
_SUMSTATS_COLUMNS = ("SNP", "CHR", "POS", "A1", "A2", "Z", "N", "FRQ")


@dataclass(frozen=True)
class WrittenSumstats:
    """Files one :meth:`SumstatsDirectoryWriter.write` call published.

    Attributes
    ----------
    paths : dict of str to str
        Written path of each file, keyed ``parquet``, ``tsv.gz`` and, for a
        munging run, ``dropped_snps``.
    primary : str
        The Parquet file when one was written, otherwise the ``.sumstats.gz``
        file.
    parquet_row_groups : list of dict
        One entry per Parquet row group (``chrom``, ``row_group_index``,
        ``row_offset``, ``n_rows``); empty without Parquet output.
    """

    paths: dict[str, str]
    primary: str
    parquet_row_groups: list[dict[str, Any]]


class SumstatsDirectoryWriter:
    """Write munged summary statistics and the dropped-SNP table."""

    @staticmethod
    def artifact_family(output_dir, *, trait_name: str | None = None, output_format: str = "parquet",
                        run: bool = False) -> ArtifactFamily:
        """Declare the munged data files of one filename stem and, for a run, its diagnostics.

        Both data formats are always owned, so a successful overwrite removes
        a stale sibling format, and only ``output_format`` is produced. The
        trait label names the files: ``<stem>.parquet`` and
        ``<stem>.sumstats.gz``, or ``sumstats.parquet`` and ``sumstats.gz``
        without a label. ``run=True`` adds what a munging run writes besides
        the data: the workflow log and the dropped-SNP table, which every
        run writes (header only when empty), and the input-issue table, owned
        but not produced because a passing input check removes it.
        ``SumstatsMunger.write_output`` uses ``run=False`` and leaves those
        diagnostics alone.
        """
        output_format = normalize_output_format(output_format)
        stem = "sumstats" if trait_name is None else re.sub(r"[^\w.-]+", "_", trait_name).strip("._-") or "trait"
        entries = {
            "parquet": (f"{stem}.parquet", output_format in {"parquet", "both"}),
            "tsv.gz": ("sumstats.gz" if trait_name is None else f"{stem}.sumstats.gz",
                       output_format in {"tsv.gz", "both"}),
        }
        if run:
            entries.update({
                "log": (SUMSTATS_LOG, True),
                "dropped_snps": (DROPPED_SNPS_FILE, True),
                "input_issues": (INPUT_ISSUES_FILE, False),
            })
        return declare_artifacts(output_dir, entries, label="munged output artifact")

    def write(
        self,
        data: pd.DataFrame,
        output_dir,
        *,
        trait_name: str | None,
        output_format: str,
        footer_metadata: dict[bytes, bytes],
        dropped_snps: pd.DataFrame | None = None,
        overwrite: bool = False,
        preflighted: bool = False,
    ) -> WrittenSumstats:
        """Write the selected munged formats and, for a run, the dropped-SNP table.

        Parameters
        ----------
        data : pandas.DataFrame
            Munged rows in the form :func:`curated_sumstats_frame` returns;
            both formats write them in this order.
        output_dir : str or os.PathLike[str]
            Output directory; created when missing.
        trait_name : str or None
            Trait label that names the files (see :meth:`artifact_family`).
        output_format : {"parquet", "tsv.gz", "both"}
            Formats to write.
        footer_metadata : dict of bytes to bytes
            ``ldsc:*`` footer keys merged into the Parquet schema metadata.
        dropped_snps : pandas.DataFrame or None, optional
            The run's dropped-SNP table rows, written to
            ``diagnostics/dropped_snps/dropped.tsv.gz``. ``SumstatsMunger.run``
            supplies them, and the declared family then also owns the log and
            the input-issue table. ``None``, the default, writes data files
            only.
        overwrite : bool, optional
            Replace existing owned files and remove stale owned siblings after
            writing. Default is ``False``.
        preflighted : bool, optional
            The workflow already checked ownership before writing its log, so
            this call only reconciles stale siblings. Default is ``False``.

        Returns
        -------
        WrittenSumstats
            Written paths, the primary data file and the Parquet row groups.

        Raises
        ------
        FileExistsError
            If an owned file exists and neither ``overwrite`` nor
            ``preflighted`` is set.
        """
        output_format = normalize_output_format(output_format)
        root = ensure_output_directory(output_dir, label="output directory")
        family = self.artifact_family(root, trait_name=trait_name, output_format=output_format,
                                      run=dropped_snps is not None)
        stale = family.preflight(overwrite=overwrite or preflighted)
        selected = [name for name in _DATA_FILES if name in family.paths]
        if dropped_snps is not None:
            selected.append("dropped_snps")
        row_groups: list[dict[str, Any]] = []
        with TemporaryDirectory(prefix=".ldsc-sumstats-", dir=root) as temporary:
            staged = {name: Path(temporary) / family.paths[name].name for name in selected}
            if "tsv.gz" in staged:
                _write_sumstats_tsv_gz(data, staged["tsv.gz"])
            if "parquet" in staged:
                row_groups = write_sumstats_parquet(data, staged["parquet"], footer_metadata)
            if "dropped_snps" in staged:
                dropped_snps.to_csv(staged["dropped_snps"], sep="\t", index=False, compression="gzip", na_rep="")
            for name, path in staged.items():
                family.paths[name].parent.mkdir(parents=True, exist_ok=True)
                os.replace(path, family.paths[name])
        remove_output_artifacts(stale)
        primary = family.paths["parquet" if "parquet" in family.paths else "tsv.gz"]
        return WrittenSumstats({name: str(family.paths[name]) for name in selected}, str(primary), row_groups)


def normalize_output_format(output_format: str) -> str:
    """Return a validated munged summary-statistics output-format token."""
    if output_format not in SUMSTATS_OUTPUT_FORMATS:
        raise LDSCUsageError(
            f"munge-sumstats does not support output_format={output_format!r}. "
            "Most likely the output format was misspelled. Use one of 'parquet', 'tsv.gz', or 'both'."
        )
    return output_format


def write_sumstats_parquet(
    frame: pd.DataFrame, path: str | os.PathLike[str], footer_metadata: dict[bytes, bytes] | None = None
) -> list[dict[str, Any]]:
    """Write snappy-compressed Parquet to ``path`` and return row-group metadata.

    ``frame`` comes from :func:`curated_sumstats_frame`, so its rows are
    already in genomic order; this writer only splits them. The Parquet
    payload keeps the munger's numeric precision. One row group is emitted
    per chromosome among complete-coordinate rows. If rows without complete
    coordinates exist, they are emitted as the last row group with ``chrom``
    recorded as ``None``. When ``footer_metadata`` is supplied, its ``ldsc:*``
    identity keys are merged into the Parquet schema footer so the file is
    self-describing downstream.
    """
    schema = pa.Schema.from_pandas(frame, preserve_index=False)
    if footer_metadata:
        schema = schema.with_metadata({**(schema.metadata or {}), **footer_metadata})
    row_groups: list[dict[str, Any]] = []
    offset = 0
    with pq.ParquetWriter(str(path), schema, compression=SUMSTATS_PARQUET_COMPRESSION) as writer:
        if frame.empty:
            writer.write_table(pa.Table.from_pandas(frame, schema=schema, preserve_index=False))
            return row_groups
        complete = ~(coordinate_missing_mask(frame["CHR"]) | coordinate_missing_mask(frame["POS"]))
        for chrom, chrom_df in frame.loc[complete].groupby("CHR", sort=False):
            writer.write_table(pa.Table.from_pandas(chrom_df, schema=schema, preserve_index=False))
            row_groups.append(
                {
                    "chrom": str(chrom),
                    "row_group_index": len(row_groups),
                    "row_offset": offset,
                    "n_rows": len(chrom_df),
                }
            )
            offset += len(chrom_df)
        missing_df = frame.loc[~complete]
        if not missing_df.empty:
            writer.write_table(pa.Table.from_pandas(missing_df, schema=schema, preserve_index=False))
            row_groups.append(
                {
                    "chrom": None,
                    "row_group_index": len(row_groups),
                    "row_offset": offset,
                    "n_rows": len(missing_df),
                }
            )
    return row_groups


def curated_sumstats_frame(data: pd.DataFrame) -> pd.DataFrame:
    """Return munged rows in the one form that the returned table and both data files share.

    Keeps the columns ``SNP, CHR, POS, A1, A2, Z, N, FRQ`` that are present,
    in that order, and adds missing ``CHR``/``POS`` columns as missing
    values. On rows with a complete ``CHR``/``POS`` pair, ``CHR`` takes its
    canonical label (no ``chr`` prefix) and ``POS`` must be a positive
    integer. A ``CHR`` value on a row without a complete pair becomes a string
    too: its canonical label where it has one, otherwise its text, so the
    column holds strings whenever any ``CHR`` value is present. A column
    without any ``CHR`` value keeps its missing values (Parquet type
    ``double`` or ``null``). Complete rows are ordered by
    chromosome sort key, ``POS`` and input order; rows without a complete pair
    follow in input order, so the Parquet writer can place them in one final
    missing-coordinate row group.
    Numeric precision is kept. ``SumstatsMunger.run`` and
    ``SumstatsMunger.write_output`` build this frame once, before any file is
    written.

    Raises
    ------
    LDSCInputError
        If a complete-coordinate row has a non-numeric or non-positive
        ``POS`` or an unsupported ``CHR`` label. Munging drops or rejects such
        rows before this step, so only a table built by hand reaches it.
    """
    columns = [col for col in _SUMSTATS_COLUMNS if col in data.columns or col in ("CHR", "POS")]
    frame = data.loc[:, [col for col in columns if col in data.columns]]
    for col in ("CHR", "POS"):
        if col not in frame.columns:
            frame.insert(columns.index(col), col, pd.NA)
    chr_missing = coordinate_missing_mask(frame["CHR"])
    pos_missing = coordinate_missing_mask(frame["POS"])
    pos_numeric = pd.to_numeric(frame["POS"], errors="coerce")
    invalid_pos = (~pos_missing) & pos_numeric.isna()
    if invalid_pos.any():
        bad_value = frame.loc[invalid_pos, "POS"].iloc[0]
        raise LDSCInputError(
            f"munge-sumstats could not write munged summary statistics because POS contains a non-numeric value "
            f"{bad_value!r}. Most likely a coordinate column was mapped incorrectly. Pass the correct --pos column or "
            "fix the POS values."
        )

    complete = ~(chr_missing | pos_missing)
    if (~chr_missing).any() and frame["CHR"].dtype.kind in "biuf":
        # Labels are strings: convert a numeric column to object first, as
        # assignment would, instead of relying on implicit upcasting.
        frame["CHR"] = frame["CHR"].astype(object)
    if complete.any():
        complete_pos = positive_int_position_series(
            frame.loc[complete, "POS"],
            context="munged summary statistics",
            label="POS",
        )
        frame.loc[complete, "CHR"] = normalize_chromosome_series(
            frame.loc[complete, "CHR"],
            context="munged summary statistics",
        ).astype(object)
        frame.loc[complete, "POS"] = complete_pos.astype("int64")
    partial = pos_missing & ~chr_missing
    if partial.any():
        values = frame.loc[partial, "CHR"]
        labels = {value: _chromosome_label(value) for value in pd.unique(values)}
        frame.loc[partial, "CHR"] = values.map(labels)

    # Sort a frame of the keys alone, then take the rows once.
    keys = pd.DataFrame({"missing": ~complete, "pos": pos_numeric.where(complete, pd.NA)})
    keys["rank"] = 10_000
    if complete.any():
        unique_chroms = pd.unique(frame.loc[complete, "CHR"])
        rank_map = {chrom: chrom_sort_key(chrom)[1] for chrom in unique_chroms}
        keys.loc[complete, "rank"] = frame.loc[complete, "CHR"].map(rank_map).astype("int64")
    keys["order"] = range(len(frame))
    keys.index = pd.RangeIndex(len(keys))
    order = keys.sort_values(by=["missing", "rank", "pos", "order"], kind="mergesort").index.to_numpy()
    frame = frame.take(order)
    frame.index = pd.RangeIndex(len(frame))
    return frame


def _chromosome_label(value: object) -> str:
    """Return the canonical label of a ``CHR`` value without a position, or its text when it has none."""
    if isinstance(value, float) and value.is_integer():
        value = int(value)  # a numeric column holding missing values was parsed as float
    try:
        return normalize_chromosome(value, context="munged summary statistics")
    except LDSCInputError:
        return str(value).strip()


def _write_sumstats_tsv_gz(frame: pd.DataFrame, path: Path) -> None:
    """Write a curated frame as gzip TSV with explicit missing fields and unrounded frequency.

    ``FRQ`` is written from its text form, so the float format does not round
    it; the caller's frame, which the returned table shares, is not changed.
    """
    if "FRQ" in frame:
        frame = frame.assign(FRQ=frame["FRQ"].astype("string"))
    frame.to_csv(path, sep="\t", index=False, float_format="%.3f", compression="gzip", na_rep="NA")
