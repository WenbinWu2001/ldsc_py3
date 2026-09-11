"""Resource-owning annotation dataset with explicit chromosome/column access."""

from dataclasses import dataclass, field
from ._annotation_storage import AnnotationShard, AnnotationWorkspace, FrameSpool
from .annotation_semantics import require_unique_annotation_names
from .chromosome_inference import normalize_chromosome
from .config import GlobalConfig
from .errors import LDSCInternalError


@dataclass
class AnnotationBundle:
    """Complete annotation dataset represented by on-demand chromosome shards.

    ``shards`` contains descriptors, never DataFrames or annotation arrays.
    Reads return detached chromosome/column/row selections owned by the caller.
    Use the bundle as a context manager or call ``close()`` to release private
    staging. A calculator borrowing this bundle must not close it. Closing
    invalidates subsequent reads but preserves canonical outputs and original
    sources. Freed arrays become reusable memory; process RSS may not fall.

    Attributes
    ----------
    baseline_columns, query_columns : list of str
        Ordered model and focal annotation names shared by all chromosomes.
    workspace : AnnotationWorkspace
        Owner of private, output-contained staging files.
    output_paths : dict
        Canonical persistent artifact paths, populated by writing workflows.
    """

    shards: dict[str, AnnotationShard]
    baseline_columns: list[str]
    query_columns: list[str]
    workspace: AnnotationWorkspace
    source_summary: dict = field(default_factory=dict)
    config_snapshot: GlobalConfig | None = None
    query_statuses: tuple = ()
    gene_list_batch: object | None = None
    input_issues: object | None = None
    identity_drops: FrameSpool | None = None
    chromosome_identity_drops: dict | None = None
    output_paths: dict = field(default_factory=dict)

    @property
    def chromosomes(self):
        return list(self.shards)

    @property
    def n_rows(self):
        return sum(shard.n_rows for shard in self.shards.values())

    def shard(self, chrom):
        """Return one descriptor while verifying its owner's lifetime."""
        self.workspace.require_open()
        return self.shards[normalize_chromosome(chrom)]

    def metadata_for_chromosome(self, chrom):
        """Load row metadata for exactly one chromosome; do not cache it."""
        return self.shard(chrom).metadata()

    def read(self, chrom, *, rows=None, columns=None):
        """Read detached values for selected SNP rows and annotation columns.

        Ordinary annotations are float32. Explicitly selecting a query batch
        and LD block bounds the returned matrix; omitting selectors requests
        all active columns/rows of this chromosome only.
        """
        columns = self.baseline_columns + self.query_columns if columns is None else list(columns)
        allowed = set(self.baseline_columns + self.query_columns)
        for name in columns:
            if name not in allowed:
                raise KeyError(name)
        return self.shard(chrom).read(rows=rows, columns=columns)

    def validate(self):
        """Validate descriptor/column contracts without loading annotation data."""
        self.workspace.require_open()
        require_unique_annotation_names(self.baseline_columns, self.query_columns)
        columns = set(self.baseline_columns + self.query_columns)
        for chrom, shard in self.shards.items():
            if not columns.issubset(shard.columns):
                raise LDSCInternalError(f"Annotation shard {chrom} is missing declared annotation columns.")
        if self.query_statuses:
            usable = [item.query for item in self.query_statuses if item.status in {"ok", "warning"}]
            if usable != self.query_columns:
                raise LDSCInternalError("Usable annotation query statuses do not match declared query columns.")

    def summary(self):
        """Return compact dimensions and provenance without reading any shard."""
        return {"n_rows": self.n_rows, "baseline_columns": list(self.baseline_columns),
                "query_columns": list(self.query_columns), "chromosomes": self.chromosomes,
                "source_summary": dict(self.source_summary)}

    def close(self):
        """Remove owned private staging; persistent outputs and sources survive."""
        self.workspace.close()

    def __enter__(self):
        self.workspace.require_open()
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()
