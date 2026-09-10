"""Prepare interval queries once and construct chromosome query artifacts."""

from dataclasses import dataclass, replace
import gzip
from pathlib import Path

import numpy as np
import pandas as pd

from ._annotation_projection import ChromosomeProjector
from ._annotation_storage import ColumnStore, FrameSpool
from ._kernel.regions import iter_bed_rows
from .chromosome_inference import normalize_chromosome
from .annotation_semantics import require_unique_annotation_names
from .query_annotations import QueryAnnotationStatus, gene_query_statuses
from .errors import LDSCInputError


@dataclass
class BedQuerySource:
    """One BED source represented by chromosome interval spool descriptors."""

    query: str
    source: str
    shards: dict[str, FrameSpool]
    failure_reason: str = ''
    details: str | None = None

    @property
    def n_intervals(self):
        return sum(spool.n_rows for spool in self.shards.values())


def prepare_bed_queries(paths, workspace, *, chunk_rows=4096):
    """Parse each BED stream once, preserving header and line validation."""
    sources = []
    for ordinal, path in enumerate(map(Path, paths)):
        shards, chunk = {}, []
        def append(rows):
            frame = pd.DataFrame(rows, columns=['chrom', 'start', 'end'])
            for chrom, group in frame.groupby('chrom', sort=False):
                spool = shards.setdefault(chrom, FrameSpool(workspace.path / f'bed-{ordinal}' / chrom))
                spool.append(group.reset_index(drop=True))
        opener = gzip.open if path.suffix.lower() == '.gz' else open
        failure, details = '', None
        try:
            with opener(path, 'rt') as stream:
                for row in iter_bed_rows(stream, str(path)):
                    chunk.append((normalize_chromosome(row.chrom), row.start, row.end))
                    if len(chunk) == chunk_rows:
                        append(chunk)
                        chunk = []
                if chunk:
                    append(chunk)
        except LDSCInputError as exc:
            failure, details = 'malformed_input', str(exc).replace(str(path), path.name)
        except (OSError, UnicodeError, EOFError) as exc:
            failure, details = 'unreadable_source', str(exc).replace(str(path), path.name)
        if failure:
            shards = {}
        sources.append(BedQuerySource(path.stem, path.name, shards, failure, details))
    return sources


def build_query_shards(bundle, *, bed_sources=(), gene_batch=None, padding_bp=0, support_kind='annotation'):
    """Project one query at a time, retaining only per-query and per-gene counts.

    Zero-valued chromosome shards are retained for globally supported queries.
    Gene controls are appended to baseline columns and never enter focal batches.
    This primitive does not apply reference-panel or LD-score variance checks.
    """
    if gene_batch is None:
        statuses = tuple(QueryAnnotationStatus(source.query, source.source, 'bed', 'ok' if source.n_intervals else 'skipped', source.failure_reason or ('' if source.n_intervals else 'empty_input'), details=source.details or (None if source.n_intervals else 'BED source contains no intervals')) for source in bed_sources)
        names = [item.query for item in statuses if item.status == 'ok']
        control = False
    else:
        statuses = gene_query_statuses(gene_batch)
        names = [item.query for item in statuses if item.status in {'ok', 'warning'}]
        control = any(item['input_role'] == 'control' for item in gene_batch.declarations)
    baseline_names = [*bundle.baseline_columns, *(['gene_control'] if control else [])]
    require_unique_annotation_names(baseline_names, names)
    totals = dict.fromkeys(names, 0)
    support = None
    if gene_batch is not None:
        selected = np.zeros(len(gene_batch.catalog.frame), dtype=bool)
        for selection in gene_batch.selections:
            selected[list(selection.catalog_indices)] = True
        support = pd.Series(pd.NA, index=gene_batch.catalog.frame.index, dtype='Int64')
        support.loc[selected] = 0
    for chrom in bundle.chromosomes:
        shard = bundle.shard(chrom)
        metadata = shard.metadata()
        projector = ChromosomeProjector(metadata)
        store = ColumnStore.create(bundle.workspace.path / f'query-{chrom}.npy', shard.n_rows, [*(['gene_control'] if control else []), *names])
        if gene_batch is not None:
            for selection in gene_batch.selections:
                if selection.query not in store.columns:
                    continue
                intervals = ((start, end) for c, start, end in selection.intervals if c == chrom)
                values = projector.project(intervals, padding_bp=padding_bp)
                store.write(0, values[:, None], columns=[selection.query])
                if selection.input_role == 'focal':
                    totals[selection.query] += int(values.sum())
            genes = gene_batch.catalog.frame.loc[selected & gene_batch.catalog.frame.chrom.eq(chrom).to_numpy()]
            support.loc[genes.index] = projector.counts(genes.start0, genes.end, padding_bp=padding_bp)
        else:
            for source in bed_sources:
                if not source.n_intervals:
                    continue
                values = np.zeros(shard.n_rows, dtype=bool)
                if chrom in source.shards:
                    for intervals in source.shards[chrom].frames():
                        values |= projector.project(zip(intervals.start, intervals.end), padding_bp=padding_bp)
                store.write(0, values[:, None], columns=[source.query])
                totals[source.query] += int(values.sum())
        bundle.shards[chrom] = replace(shard, stores=(*shard.stores, store))
        del metadata, projector
        values = genes = intervals = selection = None
    if gene_batch is not None:
        gene_batch = gene_batch.with_snp_support(support, support_kind=support_kind)
        statuses = gene_query_statuses(gene_batch)
        bundle.gene_list_batch = gene_batch
    finalized = []
    for status in statuses:
        if status.status in {'ok', 'warning'}:
            count = totals[status.query]
            # BED retains its existing standalone semantics, including nonempty
            # inputs that happen not to overlap this baseline grid.
            if gene_batch is not None and count == 0:
                status = status.updated(status='skipped', reason='zero_annotation_snps')
            status = status.updated(n_annotation_snps=float(count))
        elif status.reason == 'zero_annotation_snps':
            status = status.updated(n_annotation_snps=0.0)
        finalized.append(status)
    bundle.baseline_columns = baseline_names
    bundle.query_statuses = tuple(finalized)
    bundle.query_columns = [item.query for item in finalized if item.status in {'ok', 'warning'}]
    bundle.validate()
    return bundle
