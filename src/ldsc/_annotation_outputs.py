"""Persistent annotation results with deferred, output-owned selected reads."""

from dataclasses import dataclass, replace
from pathlib import Path

from ._annotation_sources import prepare_annotation_sources
from ._annotation_storage import AnnotationShard, AnnotationWorkspace, TsvDiagnostics
from ._gene_query_storage import persistent_gene_diagnostics
from ._kernel.snp_identity import identity_base_mode
from ._row_alignment import assert_same_snp_rows
from .errors import LDSCInputError


@dataclass(frozen=True)
class SavedAnnotationShard:
    """Persistent query file and known dimensions, without loaded values."""

    path: Path
    n_rows: int
    columns: tuple[str, ...]


@dataclass(frozen=True)
class AnnotationOutputSources:
    """Normalize original baselines and saved query outputs once on first use.

    Keeping the original baseline dependency avoids copying baseline matrices
    into the public query artifact. A subsequent consumer owns any prepared
    private shards until it closes its bundle; no gene/BED reprojection occurs.
    """

    baselines: tuple[Path, ...]
    queries: tuple[Path, ...]

    def prepare(self, workspace, mode):
        baseline = prepare_annotation_sources(AnnotationWorkspace(workspace.path), self.baselines, [], mode=mode)
        query = prepare_annotation_sources(AnnotationWorkspace(workspace.path), self.queries, [], mode=mode)
        if set(baseline.shards) != set(query.shards):
            raise LDSCInputError("Saved annotation outputs no longer match their original baseline chromosome scope. Keep the original inputs unchanged or regenerate annotations.")
        shards = {}
        for chrom, source in query.shards.items():
            base = baseline.shards[chrom]
            left, right = base.metadata(), source.metadata()
            alignment_mode = mode if all({'A1', 'A2'}.issubset(frame.columns) for frame in (left, right)) else identity_base_mode(mode)
            assert_same_snp_rows(left, right, context='Saved annotation outputs and original baseline', snp_identifier=alignment_mode)
            shards[chrom] = AnnotationShard(source.metadata_path, base.stores + source.stores, source.n_rows)
            del left, right
        return shards


def persistent_annotation_bundle(bundle, baseline_paths, output_dir):
    """Detach a completed standalone result, then release construction scratch."""
    paths = {chrom: Path(path).resolve() for chrom, path in zip(bundle.chromosomes, bundle.output_paths['query_annotations'])}
    columns = tuple(bundle.baseline_columns + bundle.query_columns)
    shards = {chrom: SavedAnnotationShard(paths[chrom], shard.n_rows, columns) for chrom, shard in bundle.shards.items()}
    loader = AnnotationOutputSources(tuple(Path(path).resolve() for path in baseline_paths), tuple(paths.values()))
    result = replace(bundle, shards=shards, workspace=AnnotationWorkspace(Path(output_dir).resolve(), defer=True),
                     gene_list_batch=persistent_gene_diagnostics(bundle.gene_list_batch, bundle.output_paths),
                     identity_drops=TsvDiagnostics(Path(bundle.output_paths['dropped_snps']).resolve()),
                     _source_loader=loader)
    bundle.close()
    return result
