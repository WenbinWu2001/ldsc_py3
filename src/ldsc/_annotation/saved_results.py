"""The record of a completed ``run_annotate`` and its explicit preparation for reads."""

from dataclasses import dataclass, field
from pathlib import Path

from .bundle import AnnotationBundle, _ShardDimensions
from .gene_query_storage import GeneListDiagnostics, persistent_gene_diagnostics
from .sources import prepare_annotation_sources
from .._runtime.scratch import ScratchWorkspace, TsvDiagnostics
from ..config import GlobalConfig


@dataclass(frozen=True)
class SavedAnnotationShard:
    """Persistent query file and known dimensions, without loaded values.

    ``columns`` are the query annotation columns of the file at ``path``,
    which holds the SNP metadata and no baseline column.
    """

    path: Path
    n_rows: int
    columns: tuple[str, ...]


@dataclass(frozen=True)
class AnnotationOutputSources:
    """The exported annotation families, read again on explicit preparation.

    Saved baseline/query files share their cleaned rows. Query-only output
    contains real generated columns, with no synthetic baseline.
    """

    baselines: tuple[Path, ...]
    queries: tuple[Path, ...]

    def prepare(self, workspace, mode):
        prepared = prepare_annotation_sources(
            workspace, self.baselines or self.queries, self.queries if self.baselines else [], mode=mode,
        )
        return prepared.shards


@dataclass(frozen=True)
class AnnotationResult(_ShardDimensions):
    """Files a completed ``run_annotate`` wrote; values are read after explicit preparation.

    The result owns no files and needs no closing. It records the exported
    annotation files with their query statuses, gene-list diagnostics and
    dropped SNPs; the original inputs are not needed afterwards. To read
    annotation values, call ``prepared_in(directory)``.

    Attributes
    ----------
    shards : dict of str to SavedAnnotationShard
        Exported query file, row count and query columns of each chromosome.
    baseline_columns, query_columns : list of str
        Ordered annotation names; a BIM grid has no baseline columns.
    query_statuses : tuple of QueryAnnotationStatus
        Outcome of every requested query, including skipped ones.
    gene_list_batch : GeneListDiagnostics or None
        Resolution summary and streamed audit of a gene-list run.
    identity_drops : TsvDiagnostics
        SNP rows removed by identity cleanup.
    output_paths : dict
        Every published artifact path.
    source_summary : dict
        Run provenance, as written to ``diagnostics/metadata.json``.
    config_snapshot : GlobalConfig
        SNP identifier mode and resolved projection build of the run.
    """

    shards: dict[str, SavedAnnotationShard]
    baseline_columns: list[str]
    query_columns: list[str]
    query_statuses: tuple
    gene_list_batch: GeneListDiagnostics | None
    identity_drops: TsvDiagnostics
    output_paths: dict
    source_summary: dict
    config_snapshot: GlobalConfig
    exported_sources: AnnotationOutputSources = field(repr=False)

    def prepared_in(self, directory):
        """Prepare the exported files for reading in a private folder of ``directory``.

        Parameters
        ----------
        directory : path-like
            Directory, created when missing, that receives one hidden
            ``.ldsc-annotation-*`` folder of packed-binary and dense-float32
            stores. Any directory is accepted, including the annotate output
            directory; a missing or blank path is rejected.

        Returns
        -------
        AnnotationBundle
            A readable bundle that owns that folder: use it as a context
            manager or call ``close()`` to remove the folder. Each call returns
            an independent bundle; the exported files are only read.
        """
        workspace = ScratchWorkspace(directory)
        try:
            shards = self.exported_sources.prepare(workspace, self.config_snapshot.snp_identifier)
        except BaseException:
            workspace.close()
            raise
        return AnnotationBundle(shards, list(self.baseline_columns), list(self.query_columns), workspace,
                                source_summary=dict(self.source_summary), config_snapshot=self.config_snapshot,
                                query_statuses=self.query_statuses, gene_list_batch=self.gene_list_batch,
                                identity_drops=self.identity_drops, output_paths=dict(self.output_paths))


def annotation_result(bundle):
    """Record a completed standalone run, then release its construction scratch."""
    paths = {chrom: Path(path).resolve() for chrom, path in zip(bundle.chromosomes, bundle.output_paths['query_annotations'])}
    columns = tuple(bundle.query_columns)
    result = AnnotationResult(
        shards={chrom: SavedAnnotationShard(paths[chrom], shard.n_rows, columns) for chrom, shard in bundle.shards.items()},
        baseline_columns=list(bundle.baseline_columns), query_columns=list(bundle.query_columns),
        query_statuses=tuple(bundle.query_statuses),
        gene_list_batch=persistent_gene_diagnostics(bundle.gene_list_batch, bundle.output_paths),
        identity_drops=TsvDiagnostics(Path(bundle.output_paths['dropped_snps']).resolve()),
        output_paths=dict(bundle.output_paths), source_summary=dict(bundle.source_summary),
        config_snapshot=bundle.config_snapshot,
        exported_sources=AnnotationOutputSources(
            tuple(Path(path).resolve() for path in bundle.output_paths.get('baseline_annotations', ())),
            tuple(paths.values())),
    )
    bundle.close()
    return result
