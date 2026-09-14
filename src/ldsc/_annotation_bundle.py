"""Resource-owning annotation dataset with explicit chromosome/column access."""

from dataclasses import dataclass, field
from ._annotation_storage import AnnotationShard, AnnotationWorkspace, FrameSpool, TsvDiagnostics
from .annotation_semantics import require_unique_annotation_names
from .chromosome_inference import normalize_chromosome
from .config import GlobalConfig
from .errors import LDSCInternalError


@dataclass
class AnnotationBundle:
    """Complete annotation dataset represented by on-demand chromosome shards.

    File-backed ``shards`` contain descriptors. ``from_frames`` instead owns
    arrays and diagnostics in memory for small prepared-input calculations.
    Reads return detached chromosome/column/row selections owned by the caller.
    Use the bundle as a context manager or call ``close()`` to release private
    staging. A calculator borrowing this bundle must not close it. Closing
    invalidates subsequent reads but preserves canonical outputs and original
    sources. Freed arrays become reusable memory; process RSS may not fall.

    Attributes
    ----------
    baseline_columns, query_columns : list of str
        Ordered model and focal annotation names shared by all chromosomes.
    workspace : AnnotationWorkspace or None
        Owner of private, output-contained staging files; ``None`` for the
        file-free ``from_frames`` route.
    output_paths : dict
        Canonical persistent artifact paths, populated by writing workflows.
    """

    shards: dict[str, AnnotationShard]
    baseline_columns: list[str]
    query_columns: list[str]
    workspace: AnnotationWorkspace | None
    source_summary: dict = field(default_factory=dict)
    config_snapshot: GlobalConfig | None = None
    query_statuses: tuple = ()
    gene_list_batch: object | None = None
    input_issues: object | None = None
    identity_drops: FrameSpool | TsvDiagnostics | None = None
    chromosome_identity_drops: dict | None = None
    output_paths: dict = field(default_factory=dict)
    _source_loader: object | None = field(default=None, repr=False)
    diagnostics_in_memory: bool = False
    query_preparation: object | None = field(default=None, repr=False)
    _closed: bool = field(default=False, repr=False)

    @classmethod
    def from_frames(cls, metadata, baseline_annotations, query_annotations=None, *, config_snapshot=None):
        """Prepare small, aligned DataFrames entirely in memory.

        Parameters
        ----------
        metadata : pandas.DataFrame
            SNP rows with ``CHR`` and 1-based integer ``POS`` in the configured
            build; ``SNP`` is required for rsID identity modes. Chromosome
            labels are normalized. Supply both ``A1`` and ``A2``, or neither;
            allele-free inputs use base identity matching. Optional ``CM`` and
            ``MAF`` do not override reference-panel values during LD scoring.
            Rows align positionally with the annotation frames, ignoring
            DataFrame index labels. No liftover or build inference is performed.
        baseline_annotations : pandas.DataFrame
            Nonempty set of numeric annotation columns, with one row per
            metadata row. Names must be unique across baseline and queries and
            must not overlap metadata names. Values are copied to float32 and
            must remain finite after conversion, as for file preparation.
        query_annotations : pandas.DataFrame or None, optional
            Optional query columns subject to the same row and value rules.
            Defaults to ``None`` for a baseline-only bundle.
        config_snapshot : GlobalConfig or None, optional
            SNP identity and genome-build settings; defaults to global config.

        Returns
        -------
        AnnotationBundle
            Owned in-memory chromosome values and identity-drop diagnostics.
            All members of duplicate effective-SNP-key groups are dropped
            globally before rows are ordered by position within chromosomes.
            Inputs are not mutated. This constructor creates no files; callers
            accept its resident RAM and can close it with a context manager.

        Raises
        ------
        LDSCInputError
            Required identity columns, paired alleles, annotation names, row
            counts, or finite-value requirements are violated.
        ValueError
            Coordinates or annotation values cannot be converted to their
            numeric types.
        """
        from ._annotation_memory import bundle_from_frames
        return bundle_from_frames(cls, metadata, baseline_annotations, query_annotations, config_snapshot)

    def _require_open(self):
        if self._closed:
            raise ValueError("Annotation bundle is closed.")
        if self.workspace is not None:
            self.workspace.require_open()

    @property
    def chromosomes(self):
        return list(self.shards)

    @property
    def n_rows(self):
        return sum(shard.n_rows for shard in self.shards.values())

    def shard(self, chrom):
        """Return one descriptor while verifying its owner's lifetime."""
        self._require_open()
        self._prepare_sources()
        return self.shards[normalize_chromosome(chrom)]

    def _prepare_sources(self):
        if self._source_loader is not None:
            try:
                self.shards = self._source_loader.prepare(self.workspace, self.config_snapshot.snp_identifier)
            except BaseException:
                self.workspace.close()
                raise
            self._source_loader = None

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
        """Prepare deferred sources if needed, then validate descriptor contracts."""
        self._require_open()
        self._prepare_sources()
        require_unique_annotation_names(self.baseline_columns, self.query_columns)
        columns = set(self.baseline_columns + ([] if self.query_preparation is not None else self.query_columns))
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
        if self.workspace is not None:
            self.workspace.close()
        self._closed = True

    def __enter__(self):
        self._require_open()
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()
