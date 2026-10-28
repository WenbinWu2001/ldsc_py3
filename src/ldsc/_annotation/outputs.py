"""Writer for the ``ldsc annotate`` result directory.

``AnnotationDirectoryWriter`` declares the files an annotate run owns and
exports aligned baseline and query annotation shards in bounded tiles, with
the query diagnostics under ``diagnostics/``.
"""

from __future__ import annotations

import gzip
import json
from pathlib import Path

import numpy as np
import pandas as pd

from .._runtime.artifacts import ArtifactFamily, _atomic_write_json, _declare_artifacts
from .._runtime.paths import remove_output_artifacts


class AnnotationDirectoryWriter:
    """Export aligned baseline/query shards in bounded tiles and query diagnostics."""

    @staticmethod
    def artifact_family(output_dir, *, chromosomes=(), gene_lists=False, diagnostics=(), baseline=False) -> ArtifactFamily:
        output_dir = Path(output_dir)
        entries = {
            'metadata': ('diagnostics/metadata.json', True),
            'log': ('diagnostics/annotate.log', True),
            'dropped_snps': ('diagnostics/dropped_snps/dropped.tsv.gz', True),
            'query_annotation_status': ('diagnostics/query_annotation_status.tsv', True),
            'gene_list_audit': ('diagnostics/gene_list_audit.tsv.gz', gene_lists),
            'gene_list_resolution_summary': ('diagnostics/gene_list_resolution_summary.tsv', gene_lists),
            'gene_catalog_issues': ('diagnostics/gene_catalog_issues.tsv', 'gene_catalog_issues' in diagnostics),
            'input_issues': ('diagnostics/input_issues.tsv', 'input_issues' in diagnostics),
            'chromosome_scope': ('diagnostics/chromosome_scope.json', 'chromosome_scope' in diagnostics),
        }
        entries.update({f'query_{chrom}': (f'query.{chrom}.annot.gz', True) for chrom in chromosomes})
        if baseline:
            entries.update({f'baseline_{chrom}': (f'baseline.{chrom}.annot.gz', True) for chrom in chromosomes})
        return _declare_artifacts(output_dir, entries, label='annotation output artifact',
                                  extra_owned=(*output_dir.glob('query.*.annot.gz'), *output_dir.glob('baseline.*.annot.gz')))

    def write_diagnostics(self, output_dir, *, batch=None, statuses=None, drops=None, catalog_issues=None, input_issues=None, scope=None):
        """Persist every available gate result without inventing later-stage data."""
        from .._core.snp_identity import IDENTITY_DROP_COLUMNS

        root = Path(output_dir) / 'diagnostics'
        root.mkdir(parents=True, exist_ok=True)
        paths = {}
        if batch is not None:
            paths['gene_list_audit'] = root / 'gene_list_audit.tsv.gz'
            paths['gene_list_resolution_summary'] = root / 'gene_list_resolution_summary.tsv'
            batch.write_audit(paths['gene_list_audit'])
            batch.summary.to_csv(paths['gene_list_resolution_summary'], sep='\t', index=False, na_rep='')
        if statuses is not None:
            paths['query_annotation_status'] = root / 'query_annotation_status.tsv'
            columns = ['query', 'source', 'input_type', 'status', 'reason', 'n_annotation_snps', 'details']
            pd.DataFrame([item.as_dict() for item in statuses], columns=columns).to_csv(paths['query_annotation_status'], sep='\t', index=False, na_rep='')
        if drops is not None:
            paths['dropped_snps'] = root / 'dropped_snps' / 'dropped.tsv.gz'
            drops.write_tsv(paths['dropped_snps'], columns=IDENTITY_DROP_COLUMNS)
        for key, frame in [('gene_catalog_issues', catalog_issues), ('input_issues', input_issues)]:
            if frame is not None:
                paths[key] = root / f'{key}.tsv'
                frame.to_csv(paths[key], sep='\t', index=False, na_rep='')
        if scope is not None:
            paths['chromosome_scope'] = root / 'chromosome_scope.json'
            _atomic_write_json(scope, paths['chromosome_scope'])
        return paths

    def write(self, bundle, output_dir, *, overwrite=False, provenance=None, scope=None, catalog_issues=None, preflighted=False):
        """Write supplied baselines as float32 and generated queries as 0/1 text.

        Both families share cleaned metadata, rows, and chromosome order. Numeric
        tile width is bounded by the family column count; source files are never
        modified. Query-only bundles emit no baseline family.
        """
        output_dir = Path(output_dir)
        diagnostics = ['chromosome_scope'] if scope is not None else []
        if catalog_issues is not None:
            diagnostics.append('gene_catalog_issues')
        family = self.artifact_family(output_dir, chromosomes=bundle.chromosomes, gene_lists=bundle.gene_list_batch is not None,
                                      diagnostics=diagnostics, baseline=bool(bundle.baseline_columns))
        # The workflow already preflighted before writing its log/diagnostics.
        # Reconcile the final produced family for successful stale cleanup.
        stale = family.preflight(overwrite=overwrite or preflighted)
        paths = self.write_diagnostics(output_dir, batch=bundle.gene_list_batch, statuses=bundle.query_statuses,
                                       drops=bundle.identity_drops, scope=scope, catalog_issues=catalog_issues)
        groups = [('query', bundle.query_columns)]
        if bundle.baseline_columns:
            groups.insert(0, ('baseline', bundle.baseline_columns))
        written = {f'{role}_annotations': [] for role, _ in groups}
        for chrom in bundle.chromosomes:
            metadata = bundle.metadata_for_chromosome(chrom).rename(columns={'POS': 'BP'})
            preferred = [name for name in ('CHR', 'BP', 'SNP', 'CM') if name in metadata]
            metadata = metadata.loc[:, [*preferred, *[name for name in metadata if name not in preferred]]]
            for role, columns in groups:
                path = output_dir / f'{role}.{chrom}.annot.gz'
                row_batch = max(1, min(65536, 16 * 1024 * 1024 // (8 * max(1, len(columns) + len(metadata.columns)))))
                with gzip.open(path, 'wt') as stream:
                    for start in range(0, len(metadata), row_batch):
                        stop = min(start + row_batch, len(metadata))
                        values = bundle.read(chrom, rows=slice(start, stop), columns=columns)
                        annotations = pd.DataFrame(values.astype(np.uint8) if role == 'query' else values,
                                                   columns=columns)
                        frame = pd.concat([metadata.iloc[start:stop].reset_index(drop=True), annotations], axis=1)
                        frame.to_csv(stream, sep='\t', index=False, header=start == 0, na_rep='NA')
                written[f'{role}_annotations'].append(path)
            del metadata
        config = bundle.config_snapshot
        payload = {**(provenance or {}), 'artifact_type': 'annotation_projection',
                   'snp_identifier': config.snp_identifier, 'genome_build': config.genome_build,
                   'chromosomes': bundle.chromosomes, 'n_snps': bundle.n_rows,
                   'baseline_columns': bundle.baseline_columns, 'query_columns': bundle.query_columns,
                   'files': {**{key: [path.name for path in files] for key, files in written.items()},
                             **{key: str(path.relative_to(output_dir)) for key, path in paths.items()}}}
        _atomic_write_json(payload, family.paths['metadata'])
        remove_output_artifacts(stale)
        bundle.output_paths = {**{key: str(path) for key, path in paths.items()}, 'metadata': str(family.paths['metadata']),
                               **{key: [str(path) for path in files] for key, files in written.items()}}
        return bundle.output_paths
