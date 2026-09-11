"""Staged gene-source resolution, selections, and complete diagnostic replay.

One catalog lookup is shared. At most one source's selected catalog indices and
one audit chunk are loaded at a time; pathway selections are never collected
into an all-pathway interval matrix.
"""

from dataclasses import dataclass, replace
import gzip
from pathlib import Path
import shutil
import zlib

import numpy as np
import pandas as pd

from ._annotation_storage import FrameSpool
from .gene_list_resolver import (
    AUDIT_COLUMNS, SUMMARY_COLUMNS, RESOLVED_ONLY_REASONS, GeneSourceSelection,
    _expand_focal_gene_list_sources, _resolve_rows, _summarize_sources,
    gene_list_query_name,
)
from .errors import LDSCInputError
from .path_resolution import normalize_path_token


@dataclass(frozen=True)
class GeneListDiagnostics:
    """Completed resolution summaries and a persistent, streamed row audit."""

    summary: pd.DataFrame
    resolution_policy: str
    audit_path: Path

    def audit_frames(self):
        """Read bounded audit chunks from the canonical output artifact."""
        yield from pd.read_csv(self.audit_path, sep='\t', chunksize=4096)

    def write_audit(self, path):
        """Copy an existing complete audit when writing another output directory."""
        if Path(path).resolve() != self.audit_path.resolve():
            shutil.copyfile(self.audit_path, path)


def persistent_gene_diagnostics(batch, output_paths):
    """Detach a completed result from private resolution storage after writing."""
    if batch is None:
        return None
    return GeneListDiagnostics(batch.summary, batch.resolution_policy, Path(output_paths['gene_list_audit']))


@dataclass(frozen=True)
class StagedGeneListBatch:
    """Compact batch result owning paths through the caller's annotation workspace."""

    catalog: object
    path: Path
    declarations: tuple[dict, ...]
    audit_spools: tuple[FrameSpool, ...]
    summary: pd.DataFrame
    resolution_policy: str
    has_fatal_gate_a_issues: bool
    coverage_scope: frozenset[str] | None = None
    support: pd.Series | None = None
    support_kind: str = "reference"

    def selection(self, input_role, source_ordinal):
        declaration = next(item for item in self.declarations if item['input_role'] == input_role and item['source_ordinal'] == source_ordinal)
        indices = np.load(self.path / f"{input_role}-{source_ordinal}.npy")
        frame = self.catalog.frame
        positions = frame.index.get_indexer(indices)
        return GeneSourceSelection(
            input_role, declaration['argument'], declaration['query'], declaration['source'], source_ordinal,
            tuple(map(str, frame.gene_id.to_numpy()[positions])), tuple(indices.tolist()),
            tuple(zip(map(str, frame.chrom.to_numpy()[positions]), map(int, frame.start0.to_numpy()[positions]),
                      map(int, frame.end.to_numpy()[positions]))),
        )

    @property
    def selections(self):
        """Iterate one source selection at a time; no shared selection cache."""
        return (self.selection(item['input_role'], item['source_ordinal']) for item in self.declarations)

    def audit_frames(self):
        """Replay every source row with applicable coverage/support annotations."""
        for declaration, spool in zip(self.declarations, self.audit_spools):
            selection = self.selection(declaration['input_role'], declaration['source_ordinal'])
            selected = set(selection.canonical_gene_ids)
            support = None if self.support is None else dict(zip(selection.canonical_gene_ids, self.support.reindex(selection.catalog_indices)))
            for original in spool.frames():
                frame = original.copy()
                eligible = frame.canonical_gene_id.isin(selected)
                if self.coverage_scope is not None:
                    frame['coverage_status'] = frame['coverage_status'].astype(object)
                    frame.loc[eligible, 'coverage_status'] = frame.loc[eligible, 'chrom'].map(lambda c: 'covered' if str(c) in self.coverage_scope else 'uncovered')
                if support is not None:
                    column = f'{self.support_kind}_snp_count'
                    frame[column] = frame.canonical_gene_id.map(support).astype('Int64')
                    unsupported = frame.disposition.eq('retained') & frame[column].eq(0).fillna(False)
                    frame.loc[unsupported, 'disposition'] = 'unsupported'
                    frame.loc[unsupported, 'reason'] = f'zero_{self.support_kind}_snp_support'
                    universe = 'annotation-grid' if self.support_kind == 'annotation' else 'retained reference-panel'
                    frame.loc[unsupported, 'details'] = f'No {universe} SNP overlaps this gene interval.'
                yield frame

    def with_snp_support(self, counts, *, support_kind='reference'):
        """Attach measured per-gene counts without retaining an enriched audit."""
        if support_kind not in {'reference', 'annotation'}:
            raise ValueError('support_kind must be reference or annotation')
        support = pd.Series(counts, dtype='Int64')
        summary = self.summary.copy()
        for selection in self.selections:
            mask = summary.input_role.eq(selection.input_role) & summary.source_ordinal.eq(selection.source_ordinal)
            values = support.reindex(selection.catalog_indices)
            summary.loc[mask, 'zero_support_genes'] = pd.NA if values.isna().any() else int(values.eq(0).sum())
            summary.loc[mask, 'genes_with_snp_support'] = pd.NA if values.isna().any() else int(values.gt(0).sum())
        return replace(self, summary=summary, support=support, support_kind=support_kind)

    def with_coverage(self, chromosomes):
        """Check selected, post-exclusion genes against validated input scope."""
        scope = frozenset(map(str, chromosomes))
        summary, errors = self.summary.copy(), []
        for column in ('coverage_status', 'missing_chromosomes', 'uncovered_gene_ids'):
            summary[column] = summary[column].astype(object)
        for selection in self.selections:
            mask = summary.input_role.eq(selection.input_role) & summary.source_ordinal.eq(selection.source_ordinal)
            missing = [(gene, interval[0]) for gene, interval in zip(selection.canonical_gene_ids, selection.intervals) if interval[0] not in scope]
            selected, uncovered = len(selection.canonical_gene_ids), len(missing)
            status = 'empty' if not selected else 'full' if not missing else 'none' if uncovered == selected else 'partial'
            missing_chroms = ','.join(sorted({c for _, c in missing}, key=int))
            summary.loc[mask, ['coverage_status', 'selected_genes', 'covered_genes', 'uncovered_genes', 'missing_chromosomes', 'uncovered_gene_ids']] = [status, selected, selected-uncovered, uncovered, missing_chroms, ','.join(g for g, _ in missing)]
            if missing:
                errors.append(f"{selection.input_role} {selection.query!r}: incomplete chromosome coverage ({selected-uncovered}/{selected} selected genes covered); missing chromosomes {missing_chroms}; affected genes {', '.join(g for g, _ in missing)}. Supply matching inputs covering these genes or explicitly revise the submitted list.")
        return replace(self, summary=summary, coverage_scope=scope), errors

    def write_audit(self, path):
        """Write the row-complete audit in deterministic source/line order."""
        with gzip.open(path, 'wt') as stream:
            header = True
            for frame in self.audit_frames():
                frame.to_csv(stream, sep='\t', index=False, header=header, na_rep='')
                header = False
            if header:
                columns = [*AUDIT_COLUMNS, *(['annotation_snp_count'] if self.support_kind == 'annotation' else [])]
                pd.DataFrame(columns=columns).to_csv(stream, sep='\t', index=False)


def _source_chunks(declaration, chunk_rows):
    path = Path(declaration['source_path'])
    opener = gzip.open if path.name.lower().endswith('.gz') else open
    with opener(path, 'rt', encoding='utf-8') as stream:
        rows = []
        for line, text in enumerate(stream, 1):
            value = text.strip()
            if not value:
                continue
            rows.append({**{k: declaration[k] for k in ('argument', 'input_role', 'query', 'source', 'source_ordinal')},
                         'line': line, 'input_gene': value, '_malformed': '\t' in text})
            if len(rows) == chunk_rows:
                yield pd.DataFrame(rows)
                rows = []
        if rows:
            yield pd.DataFrame(rows)


def resolve_gene_lists_staged(focal_paths, catalog, workspace, *, control_path=None, resolution_policy='strict', gene_exclude_regions='none', chunk_rows=4096):
    """Resolve all safely readable sources with bounded rows and staged audits.

    Deduplication retains each canonical gene's first source line, including
    duplicates crossing chunks. A late decoding/decompression failure marks
    the entire source unreadable and discards its partial resolution records.
    """
    if resolution_policy not in {'strict', 'resolved-only'} or gene_exclude_regions not in {'none', 'mhc'}:
        raise LDSCInputError('Use strict/resolved-only resolution and none/mhc gene exclusion.')
    workspace.require_open()
    path = workspace.path / 'gene-lists'
    path.mkdir()
    declarations = [{'argument': '--query-annot-gene-list-sources', 'input_role': 'focal', 'query': gene_list_query_name(source),
                     'source': Path(source).name, 'source_path': str(source), 'source_ordinal': i}
                    for i, source in enumerate(_expand_focal_gene_list_sources(focal_paths), 1)]
    if control_path is not None:
        source = normalize_path_token(control_path)
        declarations.append({'argument': '--control-gene-list-file', 'input_role': 'control', 'query': 'gene_control',
                             'source': Path(source).name, 'source_path': source, 'source_ordinal': 0})
    name_counts = pd.Series([d['query'] for d in declarations if d['input_role'] == 'focal']).value_counts()
    catalog_indices = catalog.frame.reset_index().drop_duplicates('gene_id').set_index('gene_id')['index']
    seed = pd.DataFrame([{**declaration, 'source_status': 'ok', 'source_reasons': ''} for declaration in declarations],
                        columns=SUMMARY_COLUMNS[:7])
    empty_summaries = _summarize_sources(pd.DataFrame(columns=AUDIT_COLUMNS), seed, resolution_policy, support_evaluated=False)
    summaries, spools, fatal = [], [], False
    for source_index, declaration in enumerate(declarations):
        label = f"{declaration['input_role']}-{declaration['source_ordinal']}"
        spool = FrameSpool(path / label)
        seen, selected, rejected_reasons = {}, [], set()
        counts = empty_summaries.iloc[source_index].copy()
        source_errors = []
        if declaration['input_role'] == 'focal' and name_counts[declaration['query']] > 1:
            source_errors.append('duplicate_query_name')
        try:
            for rows in _source_chunks(declaration, chunk_rows):
                # Resolve each distinct submitted token once per chunk. This
                # prevents repeated ambiguous identifiers from multiplying the
                # catalog candidate table before their audit rows are replayed.
                unique = rows.drop_duplicates(['input_gene', '_malformed']).reset_index(drop=True)
                resolved = _resolve_rows(unique, catalog, gene_exclude_regions)
                resolved['_malformed'] = unique['_malformed'].to_numpy()
                audit = rows.merge(resolved.loc[:, ['input_gene', '_malformed', *AUDIT_COLUMNS[7:]]],
                                   on=['input_gene', '_malformed'], how='left', sort=False).reindex(columns=AUDIT_COLUMNS)
                for index in audit.index[audit.canonical_gene_id.notna() & ~audit.disposition.eq('rejected')]:
                    gene = str(audit.at[index, 'canonical_gene_id'])
                    if gene in seen:
                        audit.loc[index, ['disposition', 'reason', 'details']] = ['duplicate', 'duplicate_canonical_gene', f'first occurrence at line {seen[gene]}']
                    else:
                        seen[gene] = int(audit.at[index, 'line'])
                        if audit.at[index, 'disposition'] == 'retained':
                            selected.append(int(catalog_indices[gene]))
                rejected_reasons.update(audit.loc[audit.disposition.eq('rejected'), 'reason'].astype(str))
                # Each chunk belongs to one source: direct totals avoid a
                # groupby/merge and schema reconstruction for every pathway.
                counts['nonblank_input_rows'] += len(audit)
                counts['uniquely_resolved_rows'] += int((audit.canonical_gene_id.notna() & ~audit.disposition.eq('rejected')).sum())
                for column, disposition in (('rejected_rows', 'rejected'), ('duplicate_rows', 'duplicate'), ('excluded_genes', 'excluded')):
                    counts[column] += int(audit.disposition.eq(disposition).sum())
                spool.append(audit)
        except UnicodeDecodeError:
            source_errors.append('invalid_utf8')
        except (gzip.BadGzipFile, EOFError, zlib.error):
            source_errors.append('invalid_gzip')
        except OSError:
            source_errors.append('unreadable_gene_list')
        unreadable = any(reason != 'duplicate_query_name' for reason in source_errors)
        if unreadable:
            if spool.path.exists():
                shutil.rmtree(spool.path)
            spool, selected = FrameSpool(path / label), []
        if unreadable:
            counts = empty_summaries.iloc[source_index].copy()
        counts['unique_resolved_genes'] = len(seen) if not unreadable else pd.NA
        counts['source_reasons'] = ';'.join(source_errors)
        counts['source_status'] = 'error' if source_errors else 'ok'
        if source_errors:
            for column in ('nonblank_input_rows', 'uniquely_resolved_rows', 'rejected_rows', 'unique_resolved_genes', 'duplicate_rows', 'excluded_genes'):
                counts[column] = pd.NA
            counts['resolution_fraction'] = pd.NA
        elif counts['nonblank_input_rows']:
            counts['resolution_fraction'] = counts['uniquely_resolved_rows'] / counts['nonblank_input_rows']
        fatal |= bool(source_errors) or (bool(rejected_reasons) and (resolution_policy == 'strict' or not rejected_reasons.issubset(RESOLVED_ONLY_REASONS)))
        np.save(path / f'{label}.npy', np.asarray(selected, dtype=np.int64))
        summaries.append(counts)
        spools.append(spool)
    summary = pd.DataFrame(summaries).reindex(columns=SUMMARY_COLUMNS).reset_index(drop=True)
    for column in ('nonblank_input_rows','uniquely_resolved_rows','rejected_rows','unique_resolved_genes','duplicate_rows','excluded_genes'):
        summary[column] = pd.to_numeric(summary[column],errors='coerce').astype('Int64')
    summary['resolution_fraction'] = summary['resolution_fraction'].astype(object)
    return StagedGeneListBatch(catalog, path, tuple(declarations), tuple(spools), summary, resolution_policy, fatal)
