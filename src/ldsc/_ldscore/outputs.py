"""Writer for the LD-score directory.

``LDScoreDirectoryWriter`` declares the files an LD-score run owns and
publishes ``ldscore.baseline.parquet``, the query batch files, the optional
``ldscore.overlap.parquet`` and the root ``metadata.json`` manifest, with
one Parquet row group per chromosome. Publication consumes one query batch
at a time, stages every file inside the output directory, and moves the
staged files into place after all batches succeed, ``metadata.json`` last.
Only the shared baseline, compact counts, overlap blocks and diagnostic
descriptors survive between batches.
"""

from __future__ import annotations

import json
import os
from dataclasses import dataclass, field, replace
from os import PathLike
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import TYPE_CHECKING, Any, Iterable

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

from .._annotation.query_outcomes import gene_viability_errors
from .._annotation.semantics import require_unique_annotation_names
from .._annotation.storage import TsvDiagnostics
from .._core.snp_identity import identity_artifact_metadata
from ..config import normalize_required_path
from .source import (is_numbered_query_batch_file, load_ldscore_from_dir, query_batch_file_name,
                     validate_ldscore_baseline_table, validate_ldscore_query_table)
from ..errors import LDSCConfigError, LDSCInputError, LDSCInternalError
from .._runtime.artifacts import ArtifactFamily, atomic_write_json, declare_artifacts, to_serializable
from .._runtime.paths import ensure_output_directory, remove_output_artifacts

if TYPE_CHECKING:
    from ..config import GlobalConfig
    from .._annotation.gene_query_storage import StagedGeneListBatch
    from .._annotation.query_outcomes import QueryAnnotationStatus
    from .fragments import StagedLDScoreBatch
    from .calculator import LDScoreResult
    from .overlap_table import LDScoreOverlap
    from .source import LDScoreSource


# File names inside an LD-score directory that other modules build paths from.
BASELINE_FILE = "ldscore.baseline.parquet"
LDSCORE_LOG = "diagnostics/ldscore.log"
CHROMOSOME_LOG_DIR = "diagnostics/chromosome-logs"
REF_METADATA_GLOB = "ref_metadata/chr*_meta.tsv.gz"


def dropped_snps_file(chrom: str) -> str:
    """Return a chromosome's dropped-SNP diagnostic, relative to an LD-score directory."""
    return f"diagnostics/dropped_snps/chr{chrom}_dropped.tsv.gz"


def ref_metadata_file(chrom: str) -> str:
    """Return a chromosome's exported reference-metadata file, relative to an LD-score directory."""
    return f"ref_metadata/chr{chrom}_meta.tsv.gz"


def _cast_parquet_floats(df: pd.DataFrame) -> pd.DataFrame:
    """Return a parquet-ready frame with float64 columns narrowed to float32."""
    float64_cols = [c for c in df.columns if df[c].dtype == np.float64]
    return df.astype({c: np.float32 for c in float64_cols}, copy=False) if float64_cols else df


def _write_chromosome_aligned_parquet(
    df: pd.DataFrame, path: Path, compression: str | None
) -> list[dict]:
    """Write one parquet row group per chromosome and return row-group metadata.

    The input frame must contain a ``CHR`` column and is assumed to already be
    sorted by genomic position. ``groupby(sort=False)`` preserves that order,
    giving row groups whose offsets match the public row order in the written
    table. Float columns are narrowed through ``_cast_parquet_floats`` before
    schema construction.
    """
    df = _cast_parquet_floats(df)
    schema = pa.Schema.from_pandas(df, preserve_index=False)
    row_group_meta: list[dict] = []
    offset = 0
    with pq.ParquetWriter(path, schema, compression=compression) as writer:
        for chrom, chrom_df in df.groupby("CHR", sort=False):
            writer.write_table(pa.Table.from_pandas(chrom_df, preserve_index=False))
            row_group_meta.append(
                {
                    "chrom": str(chrom),
                    "row_group_index": len(row_group_meta),
                    "row_offset": offset,
                    "n_rows": len(chrom_df),
                }
            )
            offset += len(chrom_df)
    return row_group_meta


REGRESSION_LD_SCORE_COLUMN = "regression_ld_scores"


QUERY_STATUS_COLUMNS = ["query", "source", "input_type", "status", "reason", "n_annotation_snps", "details"]


def common_count_config(common_maf_min: float) -> dict[str, Any]:
    """Return count metadata for common-SNP counts, ``MAF >= common_maf_min``."""
    return {
        "common_reference_snp_maf_min": float(common_maf_min),
        "common_reference_snp_maf_operator": ">=",
    }


DEFAULT_COUNT_CONFIG = common_count_config(0.05)


@dataclass(frozen=True)
class QueryDiagnostics:
    """Diagnostics of an LD-score run that stopped before publication.

    The only input of the writer's diagnostics methods, which write each
    record that is present: ``chromosome_scope.json``, ``input_issues.tsv``,
    the query statuses, and the gene-list audit and resolution summary.
    """

    chromosome_scope: dict
    query_statuses: tuple[QueryAnnotationStatus, ...] = ()
    gene_list_batch: StagedGeneListBatch | None = None
    input_issues: pd.DataFrame | None = None


@dataclass(slots=True)
class _LDScorePublication:
    """Retain only shared artifacts and metadata while consuming query batches.

    Query values, the LD-score producer's SNP sets, and chromosome results
    remain owned by the active batch. Diagnostic paths here refer to publication staging.
    """

    baseline_table: pd.DataFrame
    baseline_columns: list[str]
    count_records: list[dict]
    count_config: dict
    config_snapshot: GlobalConfig
    gene_list_batch: StagedGeneListBatch | None
    chromosome_scope: dict
    index_provenance: dict | None
    annotation_types: dict[str, str]
    snp_universe_policy: dict | None
    legacy_ldscore_conversion: dict | None
    overlap: LDScoreOverlap | None
    ld_window: dict | None = None
    ld_window_snp_counts: pd.DataFrame | None = None
    reference_snp_funnel: pd.DataFrame | None = None
    query_columns: list[str] = field(default_factory=list)
    query_statuses: tuple[QueryAnnotationStatus, ...] = ()
    query_batches: list[dict] = field(default_factory=list)
    identity_drops_by_chrom: dict[str, TsvDiagnostics] = field(default_factory=dict)
    ref_metadata_chroms: tuple[str, ...] = ()


@dataclass(frozen=True)
class LDScoreOutputConfig:
    """Directory-oriented output config for LD-score directories.

    Parameters
    ----------
    output_dir : str or os.PathLike[str]
        Directory that receives ``metadata.json``,
        ``ldscore.baseline.parquet``, and optional ``ldscore.query.parquet``.
    overwrite : bool, optional
        If ``True``, replace existing LD-score directory files in
        ``output_dir`` and remove stale owned siblings after a successful
        write. If ``False``, any existing LD-score directory file raises
        ``FileExistsError`` before any parquet or metadata file is written.
        Default is ``False``.
    parquet_compression : {"snappy", "gzip", "brotli", "zstd", "none", None}, optional
        Compression codec used for parquet tables. Default is ``"snappy"``.
    """

    output_dir: str | PathLike[str]
    overwrite: bool = False
    parquet_compression: str = "snappy"

    def __post_init__(self) -> None:
        """Normalize the output directory and validate parquet compression."""
        object.__setattr__(self, "output_dir", normalize_required_path(self.output_dir))
        allowed = {"snappy", "gzip", "brotli", "zstd", "none", None}
        if self.parquet_compression not in allowed:
            raise LDSCConfigError(
                "LD-score output configuration could not use parquet_compression="
                f"{self.parquet_compression!r}. Most likely the output config was created with an unsupported "
                "compression codec. Choose one of 'snappy', 'gzip', 'brotli', 'zstd', 'none', or None."
            )


class LDScoreDirectoryWriter:
    """Write LD-score directories.

    The writer owns ``metadata.json``, one shared baseline, optional single or
    numbered query batch files, and conditional ``ldscore.overlap.parquet``.
    Annotation counts and the ordered ``query_batches`` manifest are embedded
    in root metadata. Each Parquet file has chromosome-aligned row groups.
    ``write_batches`` releases completed query values, publishes after complete
    success, and returns an ``LDScoreSource``.
    """

    # Diagnostics owned by the writer but not listed in metadata ``files``.
    METADATA_EXCLUDED_FILES = ("query_status", "gene_list_audit", "gene_list_resolution_summary",
                               "gene_catalog_issues", "chromosome_scope", "input_issues", "reference_snp_funnel")

    @staticmethod
    def artifact_family(output_dir, result: _LDScorePublication | QueryDiagnostics | None = None, *,
                        include_logs=False) -> ArtifactFamily:
        """Describe all owned LD-score directory files and the outputs selected by a result.

        With no result, describe ownership for early collision checks without
        predicting outputs. A ``QueryDiagnostics`` selects only the diagnostics
        it holds; an ``_LDScorePublication`` selects the complete directory.
        Existing chromosome drop reports remain owned even when the current
        result has no drops for that chromosome.
        """
        output_dir = Path(output_dir)
        publication = result if isinstance(result, _LDScorePublication) else None
        diagnostics = result if isinstance(result, QueryDiagnostics) else None
        scientific = publication is not None
        batch = None if result is None else result.gene_list_batch
        query_batches = publication.query_batches if scientific else ()
        entries = {
            "metadata": ("metadata.json", scientific),
            "baseline": (BASELINE_FILE, scientific),
            "query": ("ldscore.query.parquet", scientific and len(query_batches) == 1),
            "overlap": ("ldscore.overlap.parquet", scientific and publication.overlap is not None),
            "query_status": ("diagnostics/query_annotation_status.tsv",
                             result is not None and bool(result.query_statuses)),
            "chromosome_scope": ("diagnostics/chromosome_scope.json",
                                 result is not None and bool(result.chromosome_scope)),
            "input_issues": ("diagnostics/input_issues.tsv",
                             diagnostics is not None and diagnostics.input_issues is not None),
            "gene_list_audit": ("diagnostics/gene_list_audit.tsv.gz", batch is not None),
            "gene_list_resolution_summary": ("diagnostics/gene_list_resolution_summary.tsv", batch is not None),
            "gene_catalog_issues": ("diagnostics/gene_catalog_issues.tsv",
                                    batch is not None and batch.catalog is not None and not batch.catalog.issues.empty),
            "ld_window_snp_counts": ("diagnostics/ld_window_snp_counts.tsv",
                                     scientific and publication.ld_window_snp_counts is not None),
            "reference_snp_funnel": ("diagnostics/reference_snp_funnel.tsv",
                                     scientific and publication.reference_snp_funnel is not None),
        }
        if scientific:
            if len(query_batches) > 1:
                for ordinal, entry in enumerate(query_batches, 1):
                    entries[f"query_batch{ordinal:05d}"] = (entry["file"], True)
            for chrom in publication.identity_drops_by_chrom:
                entries[f"dropped_snps_chr{chrom}"] = (dropped_snps_file(chrom), True)
            for chrom in publication.ref_metadata_chroms:
                entries[f"ref_metadata_chr{chrom}"] = (ref_metadata_file(chrom), True)
        return declare_artifacts(
            output_dir, entries, label="LD-score output artifact",
            extra_owned=[*([output_dir / LDSCORE_LOG,
                            *sorted((output_dir / CHROMOSOME_LOG_DIR).glob("chr*.log"))]
                           if include_logs else []),
                         *sorted(output_dir.glob(REF_METADATA_GLOB)),
                         *sorted(output_dir.glob(dropped_snps_file("*"))),
                         *sorted(path for path in output_dir.glob("ldscore.query.batch*.parquet")
                                 if is_numbered_query_batch_file(path.name))],
        )

    def write_batches(
        self, results: Iterable[LDScoreResult | StagedLDScoreBatch],
        output_config: LDScoreOutputConfig, *, query_status_order: Iterable[str] = (),
        ref_metadata_dir: str | PathLike[str] | None = None,
    ) -> LDScoreSource:
        """Publish sequential LD-score batches and return a saved source.

        Parameters
        ----------
        results : iterable of LDScoreResult or StagedLDScoreBatch
            Nonempty sequence of genome-wide results with identical baseline
            values and SNP rows, disjoint query columns, and corresponding
            counts and overlap statistics. ``LDScoreResult`` holds materialized
            query values; private ``StagedLDScoreBatch`` descriptors hold
            ordered chromosome fragments. The producer owns each yielded
            result and its files until consumption finishes. Use a generator
            that releases the consumed batch when resumed, and shuts down any
            workers before removing their borrowed files when closed.
        output_config : LDScoreOutputConfig
            Output directory, Parquet compression, and overwrite policy.
        query_status_order : sequence of str, optional
            Original query-source order for diagnostic statuses. Empty by
            default, preserving the order produced by the batches.
        ref_metadata_dir : str or os.PathLike[str], optional
            Directory holding exported ``ref_metadata/chr*_meta.tsv.gz``
            reference-metadata files, published with the other files and
            listed in ``metadata.json`` and ``output_paths``. Default is
            ``None``: no reference metadata is exported.

        Returns
        -------
        LDScoreSource
            Saved-artifact handle with shared baseline values, compact
            statistics and diagnostics, and explicit ``read_queries`` access.
            Completed query tables are not retained.

        Raises
        ------
        FileExistsError
            An owned artifact exists and overwrite is disabled.
        LDSCInputError
            No results were yielded, a batch's query rows differ from the
            first batch's baseline rows, annotation names repeat, or the
            gene-list control or every query is unusable.

        Notes
        -----
        Each query batch is privately written under the output directory and
        released before requesting the next result. Staged batches are read
        one chromosome at a time. Retained queries arrive finalized at the
        LD-score producer's precision, direct post-kernel float32 or indexed
        float64; the writer narrows retained values to public float32 without
        changing eligibility.
        Publication retains only the shared baseline, counts, overlap, and
        provenance/diagnostics, without copying every field of the LD-score
        producer's results.
        One query batch uses ``ldscore.query.parquet``; multiple batches use
        numbered files. After
        all batches succeed, the writer publishes shared artifacts and the
        ordered ``query_batches`` manifest, with root metadata last. Handled
        failure closes the producer before removing owned staging. Staging
        adds no rollback for a later filesystem publication failure, and
        completed batches are not resumable checkpoints.
        """
        from .fragments import StagedLDScoreBatch, write_query_fragments

        root = ensure_output_directory(output_config.output_dir, label="LD-score output directory")
        self.artifact_family(root).preflight(overwrite=output_config.overwrite)
        compression = None if output_config.parquet_compression in {None, "none"} else output_config.parquet_compression
        shared = None
        query_batches, overlaps, statuses = [], [], []
        status_positions = {name: index for index, name in enumerate(query_status_order)}
        def ordered_statuses():
            return tuple(sorted(statuses, key=lambda status: status_positions.get(status.query, len(status_positions))))
        with TemporaryDirectory(prefix=".ldsc-batches-", dir=root) as temporary:
            stage = Path(temporary)
            try:
                for result in results:
                    streamed = isinstance(result, StagedLDScoreBatch)
                    if shared is None:
                        # The first batch's baseline is the run's baseline; later
                        # batches add queries whose rows must match it exactly.
                        snp_identifier = validate_ldscore_baseline_table(
                            result.baseline_table, result.baseline_columns, result.config_snapshot,
                        )
                        if not streamed:
                            validate_ldscore_query_table(
                                result.baseline_table, result.query_table, result.query_columns,
                                snp_identifier=snp_identifier,
                            )
                        shared = _LDScorePublication(
                            baseline_table=result.baseline_table,
                            baseline_columns=list(result.baseline_columns),
                            count_records=[record for record in result.count_records if record["group"] == "baseline"],
                            count_config=dict(result.count_config),
                            config_snapshot=result.config_snapshot,
                            gene_list_batch=result.gene_list_batch,
                            chromosome_scope=result.chromosome_scope,
                            index_provenance=result.index_provenance,
                            annotation_types={name: result.annotation_types.get(name, "unknown") for name in result.baseline_columns},
                            snp_universe_policy=result.snp_universe_policy,
                            legacy_ldscore_conversion=None if streamed else result.legacy_ldscore_conversion,
                            overlap=None if result.overlap is None else result.overlap.select_queries([]),
                            ld_window=result.ld_window,
                            ld_window_snp_counts=result.ld_window_snp_counts,
                            reference_snp_funnel=result.reference_snp_funnel,
                        )
                        baseline_rg = _write_chromosome_aligned_parquet(shared.baseline_table, stage / BASELINE_FILE, compression)
                        for chrom, artifact in result.identity_drops_by_chrom.items():
                            path = stage / dropped_snps_file(chrom)
                            artifact.write_to(path)
                            shared.identity_drops_by_chrom[chrom] = TsvDiagnostics(path)
                    elif not streamed:
                        validate_ldscore_query_table(
                            shared.baseline_table, result.query_table, result.query_columns,
                            snp_identifier=snp_identifier,
                        )
                    require_unique_annotation_names(shared.baseline_columns, [*shared.query_columns, *result.query_columns])
                    shared.query_columns.extend(result.query_columns)
                    shared.count_records.extend(record for record in result.count_records if record["group"] == "query")
                    shared.annotation_types.update(result.annotation_types)
                    statuses.extend(result.query_statuses)
                    if result.query_columns:
                        filename = f"query-{len(query_batches) + 1}.parquet"
                        groups = (write_query_fragments(result, shared.baseline_table, stage / filename, compression)
                                  if streamed else
                                  _write_chromosome_aligned_parquet(result.query_table, stage / filename, compression))
                        query_batches.append({"file": filename, "query_columns": list(result.query_columns), "row_groups": groups})
                    if result.overlap is not None:
                        overlaps.append(result.overlap)
                    del result
                if shared is None:
                    raise LDSCInputError("LD-score batch writing received no calculation results.")
                shared.query_statuses = ordered_statuses()
                errors = gene_viability_errors(shared.gene_list_batch, shared.query_statuses, shared)
                if errors:
                    raise LDSCInputError("; ".join(errors))
                if overlaps:
                    first = overlaps[0]
                    def blocks(name):
                        block = getattr(first, name)
                        if block is None:
                            return None
                        return pd.concat([block.loc[:, shared.baseline_columns], *[
                            getattr(overlap, name).loc[:, overlap.query_diagonal_all.index] for overlap in overlaps]], axis=1)
                    def diagonal(name):
                        if getattr(first, name) is None:
                            return None
                        return pd.concat([getattr(overlap, name) for overlap in overlaps])
                    shared.overlap = replace(first, baseline_block_all=blocks("baseline_block_all"),
                                             baseline_block_common=blocks("baseline_block_common"),
                                             query_diagonal_all=diagonal("query_diagonal_all"),
                                             query_diagonal_common=diagonal("query_diagonal_common"))
                for ordinal, entry in enumerate(query_batches, 1):
                    name = query_batch_file_name(ordinal, len(query_batches))
                    (stage / entry["file"]).rename(stage / name)
                    entry["file"] = name
                shared.query_batches = query_batches
                if ref_metadata_dir is not None:
                    chroms = []
                    for path in sorted(Path(ref_metadata_dir).glob(REF_METADATA_GLOB)):
                        chrom = path.name.removeprefix("chr").removesuffix("_meta.tsv.gz")
                        (stage / ref_metadata_file(chrom)).parent.mkdir(parents=True, exist_ok=True)
                        os.replace(path, stage / ref_metadata_file(chrom))
                        chroms.append(chrom)
                    shared.ref_metadata_chroms = tuple(chroms)
                staged_family = self.artifact_family(stage, shared)
                paths = staged_family.paths
                if shared.overlap is not None:
                    from .overlap_table import overlap_to_long_frame
                    overlap_to_long_frame(shared.overlap).to_parquet(paths["overlap"], index=False)
                self._write_query_diagnostic_files(shared, paths)
                self._write_ld_window_snp_counts(shared, paths)
                self._write_reference_snp_funnel(shared, paths)
                metadata = self._build_metadata(shared,
                    files=staged_family.metadata_files(exclude=self.METADATA_EXCLUDED_FILES),
                    baseline_rg=baseline_rg,
                    query_rg=query_batches[0]["row_groups"] if len(query_batches) == 1 else None)
                paths["metadata"].write_text(json.dumps(to_serializable(metadata), indent=2, sort_keys=True), encoding="utf-8")
                family = self.artifact_family(root, shared)
                stale = family.preflight(overwrite=output_config.overwrite)
                for key, destination in family.paths.items():
                    if key == "metadata":
                        continue
                    destination.parent.mkdir(parents=True, exist_ok=True)
                    os.replace(paths[key], destination)
                os.replace(paths["metadata"], family.paths["metadata"])
                remove_output_artifacts(stale)
            except BaseException:
                if shared is not None:
                    shared.query_statuses = ordered_statuses()
                    if shared.query_statuses or shared.gene_list_batch is not None or shared.chromosome_scope:
                        self.write_query_diagnostics(QueryDiagnostics(
                            shared.chromosome_scope, shared.query_statuses, shared.gene_list_batch), output_config)
                raise
            finally:
                close = getattr(results, "close", None)
                if close is not None:
                    close()
        from .._annotation.gene_query_storage import persistent_gene_diagnostics
        shared.baseline_table = None  # the reload below reads the published baseline
        source = load_ldscore_from_dir(str(root))
        paths = {key: str(value) for key, value in family.paths.items()}
        return replace(source, output_paths=paths, query_statuses=shared.query_statuses,
                       gene_list_batch=persistent_gene_diagnostics(shared.gene_list_batch, paths),
                       chromosome_scope=shared.chromosome_scope,
                       snp_universe_policy=shared.snp_universe_policy,
                       identity_drops_by_chrom={chrom: TsvDiagnostics(root / dropped_snps_file(chrom))
                                                for chrom in shared.identity_drops_by_chrom})

    def write_gene_list_preflight(
        self,
        batch: StagedGeneListBatch,
        output_config: LDScoreOutputConfig,
    ) -> dict[str, str]:
        """Write Gate A audit/summary artifacts without scientific outputs."""
        output_dir = ensure_output_directory(output_config.output_dir, label="LD-score output directory")
        family = self.artifact_family(output_dir, QueryDiagnostics({}, gene_list_batch=batch))
        paths = family.paths
        stale_paths = family.preflight(overwrite=output_config.overwrite)
        self._write_gene_list_diagnostic_files(batch, paths)
        remove_output_artifacts(stale_paths)
        return {name: str(path) for name, path in paths.items()}

    def write_query_diagnostics(
        self,
        diagnostics: QueryDiagnostics,
        output_config: LDScoreOutputConfig,
    ) -> dict[str, str]:
        """Write available validation diagnostics without scientific artifacts."""
        output_dir = ensure_output_directory(output_config.output_dir, label="LD-score output directory")
        family = self.artifact_family(output_dir, diagnostics)
        paths = family.paths
        if not any(name in paths for name in ("query_status", "gene_list_audit", "chromosome_scope", "input_issues")):
            raise LDSCInternalError(
                "LD-score diagnostics-only writing requires query, gene, scope, or input-issue records. "
                "Re-run with DEBUG logging and report the traceback."
            )
        stale_paths = family.preflight(overwrite=output_config.overwrite)
        self._write_query_diagnostic_files(diagnostics, paths)
        remove_output_artifacts(stale_paths)
        return {name: str(path) for name, path in paths.items()}

    @classmethod
    def write_input_preflight_diagnostics(cls, diagnostics: QueryDiagnostics, output_dir) -> dict[str, str]:
        """Write the diagnostics of a failed input preflight, keeping earlier outputs.

        The input preflight runs before any calculation, so an input mistake
        such as a mistyped path leaves an earlier run's owned outputs in place;
        unlike ``write_query_diagnostics``, this neither checks for nor removes
        them. The caller has already run the output-family preflight.
        """
        paths = cls.artifact_family(output_dir, diagnostics).paths
        cls._write_query_diagnostic_files(diagnostics, paths)
        return {name: str(path) for name, path in paths.items()}

    @staticmethod
    def _write_ld_window_snp_counts(publication: _LDScorePublication, paths: dict[str, Path]) -> None:
        """Write per-chromosome R2 matrix/kept SNP counts behind an R2 matrix ``--ld-wind-snps`` window."""
        if "ld_window_snp_counts" in paths:
            paths["ld_window_snp_counts"].parent.mkdir(parents=True, exist_ok=True)
            publication.ld_window_snp_counts.to_csv(paths["ld_window_snp_counts"], sep="\t", index=False)

    @staticmethod
    def _write_reference_snp_funnel(publication: _LDScorePublication, paths: dict[str, Path]) -> None:
        """Write every chromosome's reference-SNP funnel, including skipped chromosomes."""
        if "reference_snp_funnel" in paths:
            from .._refpanel.funnel import write_reference_snp_funnel

            write_reference_snp_funnel(publication.reference_snp_funnel, paths["reference_snp_funnel"])

    @staticmethod
    def _write_query_diagnostic_files(result: _LDScorePublication | QueryDiagnostics,
                                      paths: dict[str, Path]) -> None:
        """Serialize the query diagnostics that ``artifact_family`` selected for ``result``."""
        if "chromosome_scope" in paths:
            paths["chromosome_scope"].parent.mkdir(parents=True, exist_ok=True)
            atomic_write_json(result.chromosome_scope, paths["chromosome_scope"])
        if "input_issues" in paths:
            paths["input_issues"].parent.mkdir(parents=True, exist_ok=True)
            result.input_issues.to_csv(paths["input_issues"], sep="\t", index=False, na_rep="")
        if "query_status" in paths:
            paths["query_status"].parent.mkdir(parents=True, exist_ok=True)
            pd.DataFrame([record.as_dict() for record in result.query_statuses], columns=QUERY_STATUS_COLUMNS).to_csv(
                paths["query_status"], sep="\t", index=False, na_rep=""
            )
        if result.gene_list_batch is not None:
            LDScoreDirectoryWriter._write_gene_list_diagnostic_files(result.gene_list_batch, paths)

    @staticmethod
    def _write_gene_list_diagnostic_files(batch: StagedGeneListBatch, paths: dict[str, Path]) -> None:
        """Persist gene resolution and any latent catalog defects before logging."""
        paths["gene_list_audit"].parent.mkdir(parents=True, exist_ok=True)
        batch.write_audit(paths["gene_list_audit"])
        batch.summary.to_csv(
            paths["gene_list_resolution_summary"],
            sep="\t",
            index=False,
            na_rep="",
        )
        if "gene_catalog_issues" in paths:
            from .._annotation.query_outcomes import log_gene_catalog_issues

            batch.catalog.issues.to_csv(paths["gene_catalog_issues"], sep="\t", index=False, na_rep="")
            log_gene_catalog_issues(batch.catalog.issues)

    def _build_metadata(
        self,
        publication: _LDScorePublication,
        files: dict[str, str],
        baseline_rg: list[dict] | None = None,
        query_rg: list[dict] | None = None,
    ) -> dict[str, Any]:
        """Build the JSON metadata payload for one LD-score directory.

        The metadata always includes ``count_config`` so downstream regression
        code can report the common-SNP count threshold even when MAF metadata is
        unavailable and per-column common-SNP counts are omitted. Row-group metadata
        records the chromosome, row-group index, row offset, and row count for
        each chromosome-aligned parquet row group.
        """
        baseline_table = publication.baseline_table
        config_snapshot = publication.config_snapshot
        if config_snapshot is None:
            raise LDSCInternalError(
                "LD-score output writer could not build metadata because LDScoreResult.config_snapshot is missing. "
                "Most likely an upstream LD-score workflow returned an older or incomplete result object. "
                "Re-run with DEBUG logging and report the traceback."
            )
        identity_metadata = identity_artifact_metadata(
            artifact_type="ldscore",
            snp_identifier=config_snapshot.snp_identifier,
            genome_build=config_snapshot.genome_build,
        )
        chromosomes = baseline_table["CHR"].astype(str).drop_duplicates().tolist()
        overlap = publication.overlap
        overlap_config = None
        if overlap is not None:
            count_config = publication.count_config
            overlap_config = {
                "total_all_reference_snps": float(overlap.total_all_reference_snps),
                "total_common_reference_snps": (
                    None if overlap.total_common_reference_snps is None
                    else float(overlap.total_common_reference_snps)
                ),
                "common_maf_min": float(count_config.get("common_reference_snp_maf_min", 0.05)),
                "common_maf_operator": str(count_config.get("common_reference_snp_maf_operator", ">=")),
                "stored_block": (
                    "baseline_rows_plus_query_diagonal" if publication.query_columns else "baseline_by_baseline"
                ),
            }
        payload = {
            **identity_metadata,
            "files": dict(files),
            "chromosomes": chromosomes,
            "chromosome_scope": publication.chromosome_scope,
            "baseline_columns": list(publication.baseline_columns),
            "query_columns": list(publication.query_columns),
            "annotation_types": dict(publication.annotation_types),
            "counts": list(publication.count_records),
            "count_config": dict(publication.count_config or DEFAULT_COUNT_CONFIG),
            "overlap_config": overlap_config,
            "n_baseline_rows": int(len(baseline_table)),
            "n_query_rows": int(len(baseline_table)) if publication.query_columns else 0,
            "row_group_layout": "one_per_chromosome",
            "baseline_row_groups": baseline_rg or [],
            "query_row_groups": query_rg,
            "query_batches": list(publication.query_batches),
        }
        gene_list_batch = publication.gene_list_batch
        if publication.query_statuses:
            payload["query_diagnostics"] = {"status": "diagnostics/query_annotation_status.tsv"}
        if publication.chromosome_scope:
            payload.setdefault("query_diagnostics", {})["chromosome_scope"] = "diagnostics/chromosome_scope.json"
        if gene_list_batch is not None:
            payload["gene_list_resolution_policy"] = gene_list_batch.resolution_policy
            catalog = gene_list_batch.catalog
            if catalog is not None and not catalog.issues.empty:
                payload.setdefault("query_diagnostics", {})["gene_catalog_issues"] = "diagnostics/gene_catalog_issues.tsv"
            payload["gene_list_resolution_counts"] = {
                "nonblank_input_rows": int(gene_list_batch.summary["nonblank_input_rows"].fillna(0).sum()),
                "rejected_rows": int(gene_list_batch.summary["rejected_rows"].fillna(0).sum()),
                "unique_resolved_genes": int(gene_list_batch.summary["unique_resolved_genes"].fillna(0).sum()),
            }
            payload.setdefault("query_diagnostics", {}).update(
                {
                    "gene_list_audit": "diagnostics/gene_list_audit.tsv.gz",
                    "gene_list_resolution_summary": "diagnostics/gene_list_resolution_summary.tsv",
                }
            )
        if publication.index_provenance is not None:
            payload.update(dict(publication.index_provenance))
        if publication.snp_universe_policy is not None:
            payload["snp_universe_policy"] = dict(publication.snp_universe_policy)
        if publication.legacy_ldscore_conversion is not None:
            payload["legacy_ldscore_conversion"] = dict(publication.legacy_ldscore_conversion)
        if publication.ld_window is not None:
            payload["ld_window"] = dict(publication.ld_window)
        return payload
