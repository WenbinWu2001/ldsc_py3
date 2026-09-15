"""Bounded discovery, positional alignment, and preparation of annotation files.

The first bounded pass classifies complete columns and stages only metadata.
After identity selection, a second input pass writes retained values directly
to packed-binary and dense-continuous NumPy stores beside Parquet metadata.
"""

from dataclasses import dataclass
from contextlib import ExitStack
from itertools import zip_longest
from pathlib import Path
from time import perf_counter
import logging
import shutil

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

from ._progress import report_phase, advance

from ._annotation_identity import DiskIdentityIndex, IdentityDropSpool
from ._annotation_parsing import normalize_annotation_chunk
from ._annotation_preflight import AUTOSOMES, INPUT_ISSUE_COLUMNS, input_issue
from ._annotation_storage import AnnotationShard, ColumnStore, FrameSpool
from ._kernel.annotation import _chrom_sort_key, _annotation_parse_error_message, _validate_annotation_values
from ._kernel.snp_identity import identity_base_mode, is_allele_aware_mode
from ._row_alignment import assert_same_snp_rows
from .annotation_semantics import require_unique_annotation_names
from .chromosome_inference import normalize_chromosome
from .errors import LDSCInputError, LDSCUserError

LOGGER = logging.getLogger("LDSC.annotation")


@dataclass
class PreparedAnnotationSources:
    """Complete source dataset described by paths and compact shared metadata."""

    shards: dict[str, AnnotationShard]
    baseline_columns: tuple[str, ...]
    query_columns: tuple[str, ...]
    drops: FrameSpool
    scope_chromosomes: tuple[str, ...] = ()
    input_chromosomes: tuple = ()


@dataclass
class _Source:
    path: Path
    spool: FrameSpool
    columns: tuple[str, ...]
    metadata_columns: tuple[str, ...]
    chromosomes: set[str]
    width: int
    binary_columns: tuple[str, ...]


def _chunk_size(width):
    return max(1, min(65536, (16 * 1024 * 1024) // (max(1, width) * 8)))


@report_phase(LOGGER, 'validation/staging', 'annotation classification and metadata')
def _scan(path, directory, mode, chunk_rows):
    advance(0, object=str(path), force=True)
    directory.mkdir(parents=True, exist_ok=True)
    spool = FrameSpool(directory / "metadata")
    chromosomes = set()
    quantitative = set()
    columns = metadata_columns = ()
    try:
        with pd.read_csv(path, sep=r"\s+", compression="infer", chunksize=chunk_rows) as reader:
            for chunk_index, chunk in enumerate(reader):
                advance(len(chunk), object=str(path))
                metadata, values = normalize_annotation_chunk(chunk, path, mode, log_ignored_metadata=chunk_index == 0,
                                                              value_dtype=None)
                if metadata.empty:
                    continue
                columns, metadata_columns = tuple(values.columns), tuple(metadata.columns)
                chromosomes.update(metadata["CHR"])
                spool.append(metadata)
                quantitative.update(values.columns[~values.isin((0, 1)).all(axis=0)])
                width = len(chunk.columns)
    except (pd.errors.ParserError, UnicodeDecodeError, ValueError) as exc:
        raise LDSCInputError(_annotation_parse_error_message(path, details=str(exc))) from exc
    if not chromosomes:
        raise LDSCInputError(f"Annotation file '{path}' has no SNP rows. Supply a nonempty annotation grid.")
    return _Source(Path(path), spool, columns, metadata_columns, chromosomes, width,
                   tuple(name for name in columns if name not in quantitative))


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


def _metadata_chunks(source, chunk_rows):
    """Reblock metadata after discovering which sources describe aligned rows."""
    parts, size = [], 0
    for frame in source.spool.frames():
        start = 0
        while start < len(frame):
            stop = min(len(frame), start + chunk_rows - size)
            parts.append(frame.iloc[start:stop])
            size += stop - start
            start = stop
            if size == chunk_rows:
                yield pd.concat(parts, ignore_index=True)
                parts, size = [], 0
    if parts:
        yield pd.concat(parts, ignore_index=True)


def _aligned_metadata(sources, mode, chunk_rows):
    for chunks in zip_longest(*(_metadata_chunks(source, chunk_rows) for source in sources)):
        if any(chunk is None for chunk in chunks):
            raise LDSCInputError("Annotation SNP rows do not match across files: row count mismatch. Regenerate aligned sources.")
        reference = chunks[0]
        for source, current in zip(sources[1:], chunks[1:]):
            alignment_mode = mode
            if not all({"A1", "A2"}.issubset(frame.columns) for frame in (reference, current)):
                alignment_mode = identity_base_mode(mode)
            assert_same_snp_rows(reference, current, context=f"Annotation SNP rows do not match across files: {source.path}", snp_identifier=alignment_mode)
            if {"A1", "A2"}.issubset(current.columns) and "A1" not in reference:
                reference[["A1", "A2"]] = current[["A1", "A2"]]
        yield reference


@report_phase(LOGGER, 'staging', 'retained annotation metadata')
def _select_rows(groups, group_rows, workspace, identity, drops, mode, chrom):
    """Write final metadata and bounded row locators, counting retained rows."""
    selections, counts, writers = [], {}, {}
    selected_chrom = normalize_chromosome(chrom) if chrom is not None else None
    with ExitStack() as stack:
        for ordinal, (group, chunk_rows) in enumerate(zip(groups, group_rows)):
            selection = FrameSpool(workspace.path / f"selection-{ordinal}")
            selections.append(selection)
            offset = 0
            for metadata in _aligned_metadata(group, mode, chunk_rows):
                advance(len(metadata), object=f'group {ordinal + 1}/{len(groups)}')
                keep, dropped = identity.select(metadata)
                drops.append(dropped)
                if selected_chrom is not None:
                    keep &= metadata["CHR"].to_numpy() == selected_chrom
                retained = metadata.loc[keep]
                selection.append(pd.DataFrame({"row": np.flatnonzero(keep) + offset,
                                               "CHR": retained["CHR"].to_numpy()}))
                offset += len(metadata)
                for label, rows in retained.groupby("CHR", sort=False):
                    table = pa.Table.from_pandas(rows, preserve_index=False)
                    if label not in writers:
                        root = workspace.path / f"chrom-{label}"
                        root.mkdir()
                        writers[label] = stack.enter_context(pq.ParquetWriter(root / "metadata.parquet", table.schema))
                    writers[label].write_table(table)
                    counts[label] = counts.get(label, 0) + len(rows)
    return selections, counts


@report_phase(LOGGER, 'staging', 'retained annotation values')
def _write_values(groups, selections, group_rows, counts, workspace, columns):
    """Reread input tiles directly into final stores using retained row locators."""
    quantitative = set()
    for group in groups:
        for source in group:
            quantitative.update(set(source.columns).difference(source.binary_columns))
    layouts = [(np.dtype(bool), "binary.npy", tuple(name for name in columns if name not in quantitative)),
               (np.dtype(np.float32), "continuous.npy", tuple(name for name in columns if name in quantitative))]
    stores = {chrom: {dtype: ColumnStore.create(workspace.path / f"chrom-{chrom}" / filename, count, names, dtype=dtype)
                      for dtype, filename, names in layouts if names} for chrom, count in counts.items()}
    for group, selection, chunk_rows in zip(groups, selections, group_rows):
        for source in group:
            offsets = dict.fromkeys(counts, 0)
            locators = iter(selection.frames())
            retained = next(locators, None)
            positions = {name: index for index, name in enumerate(source.columns)}
            selected_columns = [(dtype, [name for name in names if name in positions])
                                for dtype, _, names in layouts]
            selected_columns = [(dtype, names, [positions[name] for name in names])
                                for dtype, names in selected_columns if names]
            first = 0
            # Use the identity-selection tile width so each locator frame belongs
            # to exactly one input tile, even when earlier tiles retained no rows.
            with pd.read_csv(source.path, sep=r"\s+", compression="infer", usecols=list(source.columns),
                             chunksize=chunk_rows) as reader:
                for chunk in reader:
                    stop = first + len(chunk)
                    advance(len(chunk), object=str(source.path))
                    if retained is not None and retained["row"].iloc[0] < stop:
                        if retained["row"].iloc[0] < first or retained["row"].iloc[-1] >= stop:
                            raise LDSCInputError(f"Annotation rows changed between input passes: {source.path}. Keep source files unchanged during preparation.")
                        values = _validate_annotation_values(chunk, source.columns, path=source.path).to_numpy()
                        for chrom, rows in retained.groupby("CHR", sort=False):
                            indices = rows["row"].to_numpy() - first
                            for dtype, names, positions in selected_columns:
                                stores[chrom][dtype].write(offsets[chrom], values[np.ix_(indices, positions)], columns=names)
                            offsets[chrom] += len(rows)
                        retained = next(locators, None)
                    first = stop
            if retained is not None or first != source.spool.n_rows:
                raise LDSCInputError(f"Annotation row count changed between input passes: {source.path}. Keep source files unchanged during preparation.")
            shutil.rmtree(source.spool.path.parent)
        if selection.n_rows:
            shutil.rmtree(selection.path)
    return {label: AnnotationShard(workspace.path / f"chrom-{label}" / "metadata.parquet", tuple(stores[label].values()), counts[label])
            for label in sorted(stores, key=_chrom_sort_key)}


@report_phase(LOGGER, 'preparation', 'annotation alignment and identity validation')
def prepare_annotation_sources(workspace, baseline_files, query_files, *, mode, chunk_rows=None, chrom=None, declared_chromosomes=None, input_issues=(), autosomes_only=False, header_widths=None):
    """Scan, validate, and stage aligned whole-genome or chromosome-sharded inputs.

    ``workspace`` owns every scratch file. Automatic chunk sizes target 16 MiB
    of parsed cells and at most 65,536 rows per independently scanned file.
    Metadata is reblocked using the combined width of each discovered aligned
    group; explicit ``chunk_rows`` controls both input passes and alignment. Identity
    passes inspect metadata, and global duplicates precede requested chromosome
    selection. A second input pass writes final values using retained locators;
    no numeric staging copy is created. Binary eligibility spans every input
    chunk and chromosome and is determined before float32 narrowing.
    INFO records mark input reading, identity/shard preparation, storage column
    counts and value-file bytes, and successful completion.
    """
    workspace.require_open()
    paths = [*baseline_files, *query_files]
    if not baseline_files and not input_issues:
        raise LDSCInputError("Annotation preparation requires baseline annotation sources.")
    if header_widths is None:
        from ._input_preflight import inspect_declared_inputs
        declared_inputs = inspect_declared_inputs(baseline=baseline_files, query=query_files,
            mode=mode, initial_issues=input_issues)
        header_widths = declared_inputs.widths
    started = perf_counter()
    LOGGER.info("Reading annotation inputs: baseline files=%d, query files=%d.", len(baseline_files), len(query_files))
    issues = list(input_issues)
    declared_chromosomes = declared_chromosomes or {}
    sources = []
    for i, path in enumerate(paths):
        declared = declared_chromosomes.get(str(path), '')
        try:
            scan_rows = chunk_rows
            if scan_rows is None and header_widths is not None:
                scan_rows = _chunk_size(header_widths[str(path)])
            if scan_rows is None:
                width = 1
                try:
                    width = len(pd.read_csv(path, sep=r"\s+", nrows=0).columns)
                except (OSError, EOFError, ValueError, UnicodeError):
                    pass  # The complete scan records the existing parse diagnostic.
                scan_rows = _chunk_size(width)
            source = _scan(path, workspace.path / f"source-{i}", mode, scan_rows)
            if autosomes_only and not source.chromosomes.issubset(AUTOSOMES):
                raise LDSCInputError('Annotation contents must identify an autosomal chromosome set.')
            if declared and source.chromosomes != {declared}:
                raise LDSCInputError(f'@ member for chromosome {declared} contains chromosomes {sorted(source.chromosomes)}.')
            sources.append(source)
        except (OSError, EOFError, ValueError, UnicodeError, LDSCUserError) as exc:
            issues.append(input_issue('baseline' if i < len(baseline_files) else 'query', path, declared, 'invalid_required_input', exc))
    if issues:
        error = LDSCInputError('Annotation input preflight failed: ' + '; '.join(f"{row['source']}: {row['details']}" for row in issues[:10]) + '. Complete diagnostics: diagnostics/input_issues.tsv. Other causes & fixes: docs/troubleshooting.md#annotate-input-preflight')
        error.input_issues = pd.DataFrame(issues, columns=INPUT_ISSUE_COLUMNS)
        raise error
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
    group_rows = [chunk_rows if chunk_rows is not None else _chunk_size(sum(source.width for source in group)) for group in groups]
    effective_mode = mode
    if is_allele_aware_mode(mode) and not any("A1" in s.metadata_columns for s in sources):
        effective_mode = identity_base_mode(mode)
    LOGGER.info("Checking SNP identities and preparing chromosome annotations.")
    with DiskIdentityIndex(workspace.path / "identity.sqlite", effective_mode) as identity:
        for group, rows in zip(groups, group_rows):
            for metadata in _aligned_metadata(group, mode, rows):
                advance(len(metadata), object=','.join(str(source.path) for source in group))
                identity.add(metadata)
        selections, counts = _select_rows(groups, group_rows, workspace, identity, drops, mode, chrom)
    (workspace.path / "identity.sqlite").unlink()
    shards = _write_values(groups, selections, group_rows, counts, workspace, columns)
    if not shards:
        error = LDSCInputError("annotate retained no annotation rows after SNP identity cleanup. Supply valid, unique SNP identities. Other causes & fixes: docs/troubleshooting.md#annotate-no-annotation-snp-rows-remain")
        error.annotation_drops = drops
        raise error
    scope = tuple(sorted(set().union(*(s.chromosomes for s in baselines)), key=_chrom_sort_key))
    prepared = PreparedAnnotationSources(shards, baseline_columns, query_columns, drops, scope,
        tuple((str(s.path), "baseline" if s in baselines else "query", tuple(s.chromosomes)) for s in sources))
    binary_count = sum(len(store.columns) for store in next(iter(shards.values())).stores if store.bitorder)
    LOGGER.info("Annotation storage: binary=%d packed-bit columns, continuous=%d float32 columns; value files=%s bytes.",
                binary_count, len(columns) - binary_count,
                f"{sum(store.path.stat().st_size for shard in shards.values() for store in shard.stores):,}")
    LOGGER.info("Annotation preparation complete: chromosomes=%d, retained SNPs=%s, elapsed=%.2fs.",
                len(shards), f"{sum(counts.values()):,}", perf_counter() - started)
    return prepared
