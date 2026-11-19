"""Writer for the ``ldsc annotate`` result directory.

``AnnotationDirectoryWriter`` declares the files an annotate run owns and
exports aligned baseline and query annotation shards in bounded tiles, with
the query diagnostics under ``diagnostics/``.
"""

from __future__ import annotations

import gzip
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.csv as pacsv

from .._runtime.artifacts import ArtifactFamily, atomic_write_json, declare_artifacts
from .._runtime.paths import remove_output_artifacts
from ..errors import LDSCInputError

# Level 6 compresses annotation text several times faster than level 9, the gzip.open default, for slightly larger files.
_GZIP_LEVEL = 6
_TEXT_OPTIONS = pacsv.WriteOptions(delimiter='\t', quoting_style='none', quoting_header='none')


def _metadata_table(metadata):
    """Return SNP metadata as an Arrow table in which a missing value is the text ``NA``."""
    return pa.table({name: pa.array(column.astype('string').fillna('NA') if column.isna().any() else column)
                     for name, column in metadata.items()})


def _write_annotation_file(path, bundle, chrom, metadata, columns, dtype):
    """Write one chromosome's ``columns`` as a tab-delimited ``.annot.gz`` file in bounded row tiles.

    ``metadata`` is the chromosome's ``_metadata_table`` and ``dtype`` the type the
    values are written in. Text is unquoted, so a double quote, tab or line break
    in a SNP identifier, allele or annotation name stops the write.
    """
    row_batch = max(1, min(65536, 16 * 1024 * 1024 // (8 * max(1, len(columns) + metadata.num_columns))))
    names = [*metadata.column_names, *columns]
    try:
        with gzip.open(path, 'wb', compresslevel=_GZIP_LEVEL) as stream:
            writer = None
            for start in range(0, metadata.num_rows, row_batch):
                stop = min(start + row_batch, metadata.num_rows)
                # Column-major values give Arrow each column without a copy.
                values = np.asfortranarray(bundle.read(chrom, rows=slice(start, stop), columns=columns), dtype=dtype)
                tile = pa.Table.from_arrays([*metadata.slice(start, stop - start).columns,
                                             *(pa.array(values[:, index]) for index in range(len(columns)))], names=names)
                if writer is None:
                    writer = pacsv.CSVWriter(stream, tile.schema, write_options=_TEXT_OPTIONS)
                writer.write_table(tile)
            if writer is not None:
                writer.close()
    except pa.ArrowInvalid as exc:
        raise LDSCInputError(
            f"Cannot write annotation file '{path}': a SNP identifier, allele or annotation name contains a "
            f"double quote, tab or line break, which an annotation file cannot hold ({exc}). "
            "Rename the affected value in the input files and run again."
        ) from exc


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
        return declare_artifacts(output_dir, entries, label='annotation output artifact',
                                  extra_owned=(*output_dir.glob('query.*.annot.gz'), *output_dir.glob('baseline.*.annot.gz')))

    def write_diagnostics(self, output_dir, *, batch=None, statuses=None, drops=None, catalog_issues=None, input_issues=None, scope=None):
        """Persist every available gate result without inventing later-stage data."""
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
            drops.write_to(paths['dropped_snps'])
        for key, frame in [('gene_catalog_issues', catalog_issues), ('input_issues', input_issues)]:
            if frame is not None:
                paths[key] = root / f'{key}.tsv'
                frame.to_csv(paths[key], sep='\t', index=False, na_rep='')
        if scope is not None:
            paths['chromosome_scope'] = root / 'chromosome_scope.json'
            atomic_write_json(scope, paths['chromosome_scope'])
        return paths

    def write(self, bundle, output_dir, *, overwrite=False, provenance=None, scope=None, catalog_issues=None, preflighted=False):
        """Write supplied baseline values and generated 0/1 queries as annotation text.

        Baseline values are written as the shortest text that reads back as
        the same float32 value, so an integer value has no decimal point.
        Both families share cleaned metadata, rows, and chromosome order. Numeric
        tile width is bounded by the family column count; source files are never
        modified. Query-only bundles emit no baseline family. Return the written
        artifact paths, which the caller records on its bundle.
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
            metadata = _metadata_table(metadata.loc[:, [*preferred, *[name for name in metadata if name not in preferred]]])
            for role, columns in groups:
                path = output_dir / f'{role}.{chrom}.annot.gz'
                _write_annotation_file(path, bundle, chrom, metadata, columns, np.uint8 if role == 'query' else np.float32)
                written[f'{role}_annotations'].append(path)
            del metadata
        config = bundle.config_snapshot
        payload = {**(provenance or {}), 'artifact_type': 'annotation_projection',
                   'snp_identifier': config.snp_identifier, 'genome_build': config.genome_build,
                   'chromosomes': bundle.chromosomes, 'n_snps': bundle.n_rows,
                   'baseline_columns': bundle.baseline_columns, 'query_columns': bundle.query_columns,
                   'files': {**{key: [path.name for path in files] for key, files in written.items()},
                             **{key: str(path.relative_to(output_dir)) for key, path in paths.items()}}}
        atomic_write_json(payload, family.paths['metadata'])
        remove_output_artifacts(stale)
        return {**{key: str(path) for key, path in paths.items()}, 'metadata': str(family.paths['metadata']),
                **{key: [str(path) for path in files] for key, files in written.items()}}
