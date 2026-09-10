"""Bounded discovery, positional alignment, and preparation of annotation files.

Each compressed source is scanned once into append spools. Aligned source rows
are joined before chromosome partitioning and global identity cleanup. Prepared
shards contain seekable float32 values and separate Parquet row metadata.
"""

from dataclasses import dataclass
from itertools import zip_longest
from pathlib import Path
import shutil

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

from ._annotation_identity import DiskIdentityIndex, IdentityDropSpool
from ._annotation_parsing import normalize_annotation_chunk
from ._annotation_storage import AnnotationShard, ColumnStore, FrameSpool
from ._kernel.annotation import _chrom_sort_key, _annotation_parse_error_message
from ._kernel.snp_identity import identity_base_mode, is_allele_aware_mode
from ._row_alignment import assert_same_snp_rows
from .annotation_semantics import require_unique_annotation_names
from .chromosome_inference import normalize_chromosome
from .errors import LDSCInputError


@dataclass
class PreparedAnnotationSources:
    """Complete source dataset described by paths and compact shared metadata."""

    shards: dict[str, AnnotationShard]
    baseline_columns: tuple[str, ...]
    query_columns: tuple[str, ...]
    drops: FrameSpool


@dataclass
class _Source:
    path: Path
    spool: FrameSpool
    columns: tuple[str, ...]
    metadata_columns: tuple[str, ...]
    chromosomes: set[str]


def _scan(path, directory, mode, chunk_rows):
    spool = FrameSpool(directory)
    chromosomes = set()
    columns = metadata_columns = ()
    try:
        with pd.read_csv(path, sep=r"\s+", compression="infer", chunksize=chunk_rows) as reader:
            for chunk in reader:
                metadata, values = normalize_annotation_chunk(chunk, path, mode)
                if metadata.empty:
                    continue
                columns, metadata_columns = tuple(values.columns), tuple(metadata.columns)
                chromosomes.update(metadata["CHR"])
                spool.append(pd.concat([metadata, values], axis=1))
    except (pd.errors.ParserError, UnicodeDecodeError, ValueError) as exc:
        raise LDSCInputError(_annotation_parse_error_message(path, details=str(exc))) from exc
    if not chromosomes:
        raise LDSCInputError(f"Annotation file '{path}' has no SNP rows. Supply a nonempty annotation grid.")
    return _Source(Path(path), spool, columns, metadata_columns, chromosomes)


def _layout(sources, label):
    single = [len(source.chromosomes) == 1 for source in sources]
    if any(single) and not all(single):
        raise LDSCInputError(f"Cannot mix whole-genome and chromosome-sharded {label} annotation files. Use one layout.")
    if not all(single):
        return None
    result = {}
    for source in sources:
        chrom = next(iter(source.chromosomes))
        if chrom in result:
            raise LDSCInputError(f"Multiple {label} annotation files describe chromosome {chrom}; supply one shard per chromosome.")
        result[chrom] = source
    return result


def _aligned_rows(sources, mode):
    for chunks in zip_longest(*(source.spool.frames() for source in sources)):
        if any(chunk is None for chunk in chunks):
            raise LDSCInputError("Annotation SNP rows do not match across files: row count mismatch. Regenerate aligned sources.")
        reference = chunks[0].loc[:, sources[0].metadata_columns].reset_index(drop=True)
        values = []
        for source, chunk in zip(sources, chunks):
            current = chunk.loc[:, source.metadata_columns].reset_index(drop=True)
            alignment_mode = mode
            if not all({"A1", "A2"}.issubset(frame.columns) for frame in (reference, current)):
                alignment_mode = identity_base_mode(mode)
            if source is not sources[0]:
                assert_same_snp_rows(reference, current, context=f"Annotation SNP rows do not match across files: {source.path}", snp_identifier=alignment_mode)
            if {"A1", "A2"}.issubset(current.columns) and "A1" not in reference:
                reference[["A1", "A2"]] = current[["A1", "A2"]]
            values.append(chunk.loc[:, source.columns].reset_index(drop=True))
        yield reference, pd.concat(values, axis=1)


def _finish_shard(spool, root, columns, metadata_columns):
    store = ColumnStore.create(root / "values.npy", spool.n_rows, columns)
    metadata_path = root / "metadata.parquet"
    offset = 0
    writer = None
    try:
        for frame in spool.frames():
            metadata = pa.Table.from_pandas(frame.loc[:, metadata_columns], preserve_index=False)
            if writer is None:
                writer = pq.ParquetWriter(metadata_path, metadata.schema)
            writer.write_table(metadata)
            store.write(offset, frame.loc[:, columns].to_numpy(dtype=np.float32))
            offset += len(frame)
    finally:
        if writer is not None:
            writer.close()
    shutil.rmtree(spool.path)
    return AnnotationShard(metadata_path, (store,), spool.n_rows)


def prepare_annotation_sources(workspace, baseline_files, query_files, *, mode, chunk_rows=None, chrom=None):
    """Scan, validate, and stage aligned whole-genome or chromosome-sharded inputs.

    ``workspace`` owns every scratch file. Automatic chunk sizes target 16 MiB
    of parsed numeric cells and at most 65,536 rows across all aligned sources.
    Global duplicate detection precedes requested chromosome selection.
    """
    workspace.require_open()
    paths = [*baseline_files, *query_files]
    if not baseline_files:
        raise LDSCInputError("Annotation preparation requires baseline annotation sources.")
    if chunk_rows is None:
        width = sum(len(pd.read_csv(path, sep=r"\s+", nrows=0).columns) for path in paths)
        chunk_rows = max(1, min(65536, (16 * 1024 * 1024) // (max(1, width) * 8)))
    sources = [_scan(path, workspace.path / f"source-{i}", mode, chunk_rows) for i, path in enumerate(paths)]
    baselines, queries = sources[:len(baseline_files)], sources[len(baseline_files):]
    baseline_layout = _layout(baselines, "baseline")
    query_layout = _layout(queries, "query") if queries else {}
    if baseline_layout is not None:
        if queries and (query_layout is None or set(query_layout) != set(baseline_layout)):
            raise LDSCInputError("Query annotation shards do not match baseline chromosome shards exactly. Supply matching shards.")
        groups = [[baseline_layout[c]] + ([query_layout[c]] if queries else []) for c in sorted(baseline_layout, key=_chrom_sort_key)]
        baseline_columns = baselines[0].columns
        query_columns = queries[0].columns if queries else ()
        if any(s.columns != baseline_columns for s in baselines) or any(s.columns != query_columns for s in queries):
            raise LDSCInputError("Annotation columns differ across chromosomes. Make every chromosome shard use the same annotation header.")
    else:
        if queries and query_layout is not None:
            raise LDSCInputError("Whole-genome baseline annotations require aligned whole-genome query inputs.")
        groups = [sources]
        baseline_columns = tuple(c for s in baselines for c in s.columns)
        query_columns = tuple(c for s in queries for c in s.columns)
    require_unique_annotation_names(baseline_columns, query_columns)
    columns = baseline_columns + query_columns
    drops = IdentityDropSpool(workspace.path / "identity-drops")
    spools, metadata_columns = {}, {}
    effective_mode = mode
    if is_allele_aware_mode(mode) and not any("A1" in s.metadata_columns for s in sources):
        effective_mode = identity_base_mode(mode)
    with DiskIdentityIndex(workspace.path / "identity.sqlite", effective_mode) as identity:
        for group in groups:
            for metadata, values in _aligned_rows(group, mode):
                identity.add(metadata)
        for group in groups:
            for metadata, values in _aligned_rows(group, mode):
                keep, dropped = identity.select(metadata)
                drops.append(dropped)
                if chrom is not None:
                    keep &= metadata["CHR"].to_numpy() == normalize_chromosome(chrom)
                retained = pd.concat([metadata.loc[keep], values.loc[keep]], axis=1)
                for label, rows in retained.groupby("CHR", sort=False):
                    root = workspace.path / f"chrom-{label}"
                    spool = spools.setdefault(label, FrameSpool(root / "spool"))
                    spool.append(rows.reset_index(drop=True))
                    metadata_columns[label] = tuple(metadata.columns)
    for source in sources:
        shutil.rmtree(source.spool.path)
    shards = {label: _finish_shard(spools[label], workspace.path / f"chrom-{label}", columns, metadata_columns[label]) for label in sorted(spools, key=_chrom_sort_key)}
    (workspace.path / "identity.sqlite").unlink()
    if not shards:
        error = LDSCInputError("annotate retained no annotation rows after SNP identity cleanup. Supply valid, unique SNP identities. Other causes & fixes: docs/troubleshooting.md#annotate-no-annotation-snp-rows-remain")
        error.annotation_drops = drops
        raise error
    return PreparedAnnotationSources(shards, baseline_columns, query_columns, drops)
