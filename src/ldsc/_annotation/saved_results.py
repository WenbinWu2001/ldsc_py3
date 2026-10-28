"""Persistent annotation results with deferred, output-owned selected reads."""

from dataclasses import dataclass, replace
from pathlib import Path

from .sources import prepare_annotation_sources
from .storage import AnnotationWorkspace, TsvDiagnostics
from .gene_query_storage import persistent_gene_diagnostics


@dataclass(frozen=True)
class SavedAnnotationShard:
    """Persistent query file and known dimensions, without loaded values."""

    path: Path
    n_rows: int
    columns: tuple[str, ...]


@dataclass(frozen=True)
class AnnotationOutputSources:
    """Normalize the exported annotation families once on first use.

    Saved baseline/query files share their cleaned rows. Query-only output
    contains real generated columns, with no synthetic baseline. The consumer
    owns bounded private stores until closure; original inputs are not needed.
    """

    baselines: tuple[Path, ...]
    queries: tuple[Path, ...]

    def prepare(self, workspace, mode):
        prepared = prepare_annotation_sources(
            workspace, self.baselines or self.queries, self.queries if self.baselines else [], mode=mode,
        )
        return prepared.shards


def persistent_annotation_bundle(bundle, output_dir):
    """Detach a completed standalone result, then release construction scratch."""
    paths = {chrom: Path(path).resolve() for chrom, path in zip(bundle.chromosomes, bundle.output_paths['query_annotations'])}
    columns = tuple(bundle.baseline_columns + bundle.query_columns)
    shards = {chrom: SavedAnnotationShard(paths[chrom], shard.n_rows, columns) for chrom, shard in bundle.shards.items()}
    loader = AnnotationOutputSources(tuple(Path(path).resolve() for path in bundle.output_paths.get('baseline_annotations', ())),
                                    tuple(paths.values()))
    result = replace(bundle, shards=shards, workspace=AnnotationWorkspace(Path(output_dir).resolve(), defer=True),
                     gene_list_batch=persistent_gene_diagnostics(bundle.gene_list_batch, bundle.output_paths),
                     identity_drops=TsvDiagnostics(Path(bundle.output_paths['dropped_snps']).resolve()),
                     _source_loader=loader)
    bundle.close()
    return result
