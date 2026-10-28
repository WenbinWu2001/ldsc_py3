"""Writer for the LD-score directory.

``LDScoreDirectoryWriter`` declares the files an LD-score run owns and
writes ``ldscore.baseline.parquet``, the query batch files, the optional
``ldscore.overlap.parquet`` and the root ``metadata.json`` manifest, with
one Parquet row group per chromosome.
"""

from __future__ import annotations

import gzip
import json
import os
from dataclasses import dataclass
from os import PathLike
from pathlib import Path
from typing import TYPE_CHECKING, Any, Iterable

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

from .._core.snp_identity import identity_artifact_metadata
from ..config import normalize_required_path
from .source import validate_ldscore_baseline_table, validate_ldscore_query_table
from ..errors import LDSCConfigError, LDSCInputError, LDSCInternalError
from .._runtime.artifacts import ArtifactFamily, atomic_write_json, declare_artifacts, to_serializable
from .._runtime.paths import ensure_output_directory, remove_output_artifacts

if TYPE_CHECKING:
    from .fragments import StagedLDScoreBatch
    from .calculator import LDScoreResult
    from .source import LDScoreSource



def _cast_parquet_floats(df: pd.DataFrame) -> pd.DataFrame:
    """Return a parquet-ready frame with float64 columns narrowed to float32."""
    float64_cols = [c for c in df.columns if df[c].dtype == np.float64]
    return df.astype({c: np.float32 for c in float64_cols}) if float64_cols else df


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


DEFAULT_COUNT_CONFIG = {
    "common_reference_snp_maf_min": 0.05,
    "common_reference_snp_maf_operator": ">=",
}


@dataclass(frozen=True)
class LDScoreOutputConfig:
    """Directory-oriented output config for canonical LD-score results.

    Parameters
    ----------
    output_dir : str or os.PathLike[str]
        Directory that receives ``metadata.json``,
        ``ldscore.baseline.parquet``, and optional ``ldscore.query.parquet``.
    overwrite : bool, optional
        If ``True``, replace existing canonical LD-score files in
        ``output_dir`` and remove stale owned siblings after a successful
        write. If ``False``, any existing canonical family file raises
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
    success, and returns an ``LDScoreSource``. ``write`` serializes an already
    materialized single result and returns artifact paths.
    """

    # Diagnostics owned by the writer but not listed in metadata ``files``.
    METADATA_EXCLUDED_FILES = ("query_status", "gene_list_audit", "gene_list_resolution_summary",
                               "gene_catalog_issues", "chromosome_scope", "input_issues", "reference_snp_funnel")

    @staticmethod
    def artifact_family(output_dir, result=None, *, diagnostics_only=False, gene_list_batch=None,
                        include_logs=False) -> ArtifactFamily:
        """Describe all owned LD-score directory files and the outputs selected by a result.

        With no result, describe ownership for early collision checks without
        predicting scientific outputs. Diagnostic-only writes use the same
        declaration. Existing chromosome drop reports remain owned even when
        the current result has no drops for that chromosome.
        """
        output_dir = Path(output_dir)
        scientific = result is not None and not diagnostics_only
        query_batches = getattr(result, "query_batches", ())
        batch = gene_list_batch if gene_list_batch is not None else getattr(result, "gene_list_batch", None)
        entries = {
            "metadata": ("metadata.json", scientific),
            "baseline": ("ldscore.baseline.parquet", scientific),
            "query": ("ldscore.query.parquet", scientific and (getattr(result, "query_table", None) is not None or len(query_batches) == 1)),
            "overlap": ("ldscore.overlap.parquet", scientific and getattr(result, "overlap", None) is not None),
            "query_status": ("diagnostics/query_annotation_status.tsv", bool(getattr(result, "query_statuses", ()))),
            "chromosome_scope": ("diagnostics/chromosome_scope.json", bool(getattr(result, "chromosome_scope", None) or getattr(result, "source_summary", {}).get("chromosome_scope"))),
            "input_issues": ("diagnostics/input_issues.tsv", getattr(result, "input_issues", None) is not None),
            "gene_list_audit": ("diagnostics/gene_list_audit.tsv.gz", batch is not None),
            "gene_list_resolution_summary": ("diagnostics/gene_list_resolution_summary.tsv", batch is not None),
            "gene_catalog_issues": ("diagnostics/gene_catalog_issues.tsv", getattr(batch, "catalog", None) is not None and not batch.catalog.issues.empty),
            "ld_window_snp_counts": ("diagnostics/ld_window_snp_counts.tsv",
                                     scientific and getattr(result, "ld_window_snp_counts", None) is not None),
            "reference_snp_funnel": ("diagnostics/reference_snp_funnel.tsv",
                                     scientific and getattr(result, "reference_snp_funnel", None) is not None),
        }
        if scientific:
            if len(query_batches) > 1:
                for ordinal, entry in enumerate(query_batches, 1):
                    entries[f"query_batch{ordinal:05d}"] = (entry["file"], True)
            for chrom in (getattr(result, "identity_drops_by_chrom", {}) or {}):
                entries[f"dropped_snps_chr{chrom}"] = (f"diagnostics/dropped_snps/chr{chrom}_dropped.tsv.gz", True)
        return declare_artifacts(
            output_dir, entries, label="LD-score output artifact",
            extra_owned=[*([output_dir / "diagnostics/ldscore.log",
                            *sorted((output_dir / "diagnostics/chromosome-logs").glob("chr*.log"))]
                           if include_logs else []),
                         *sorted((output_dir / "ref_metadata").glob("chr*_meta.tsv.gz")),
                         *sorted((output_dir / "diagnostics/dropped_snps").glob("chr*_dropped.tsv.gz")),
                         *sorted(path for path in output_dir.glob("ldscore.query.batch*.parquet")
                                 if path.name.removeprefix("ldscore.query.batch").removesuffix(".parquet").isdigit())],
        )

    def write_batches(
        self, results: Iterable[LDScoreResult | StagedLDScoreBatch],
        output_config: LDScoreOutputConfig, *, query_status_order: Iterable[str] = (),
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
            No results were yielded, baseline values or SNP rows disagree
            across batches, or annotation names are invalid or repeated.

        Notes
        -----
        Each query batch is privately written under the output directory and
        released before requesting the next result. Staged batches are read
        one chromosome at a time. Their producers finalize retained queries
        at direct post-kernel float32 or indexed float64 precision; the writer
        narrows retained values to public float32 without changing eligibility.
        Publication retains only the shared baseline, counts, overlap, and
        provenance/diagnostics, without copying all producer result fields.
        One query batch uses ``ldscore.query.parquet``; multiple batches use
        numbered files. After
        all batches succeed, the writer publishes shared artifacts and the
        ordered ``query_batches`` manifest, with root metadata last. Handled
        failure closes the producer before removing owned staging. Staging
        adds no rollback for a later filesystem publication failure, and
        completed batches are not resumable checkpoints.
        """
        from .batch_publication import write_ldscore_batches

        return write_ldscore_batches(self, results, output_config, query_status_order)

    def write_gene_list_preflight(
        self,
        batch: Any,
        output_config: LDScoreOutputConfig,
    ) -> dict[str, str]:
        """Write Gate A audit/summary artifacts without scientific outputs."""
        output_dir = ensure_output_directory(output_config.output_dir, label="LD-score output directory")
        family = self.artifact_family(output_dir, gene_list_batch=batch)
        paths = family.paths
        stale_paths = family.preflight(overwrite=output_config.overwrite)
        self._write_gene_list_diagnostic_files(batch, paths)
        remove_output_artifacts(stale_paths)
        return {name: str(path) for name, path in paths.items()}

    def write_query_diagnostics(
        self,
        result: Any,
        output_config: LDScoreOutputConfig,
    ) -> dict[str, str]:
        """Write available validation diagnostics without scientific artifacts."""
        output_dir = ensure_output_directory(output_config.output_dir, label="LD-score output directory")
        family = self.artifact_family(output_dir, result, diagnostics_only=True)
        paths = family.paths
        if not any(name in paths for name in ("query_status", "gene_list_audit", "chromosome_scope", "input_issues")):
            raise LDSCInternalError(
                "LD-score diagnostics-only writing requires query, gene, scope, or input-issue records. "
                "Re-run with DEBUG logging and report the traceback."
            )
        stale_paths = family.preflight(overwrite=output_config.overwrite)
        self._write_query_diagnostic_files(result, paths)
        remove_output_artifacts(stale_paths)
        return {name: str(path) for name, path in paths.items()}

    def write(self, result: Any, output_config: LDScoreOutputConfig) -> dict[str, str]:
        """Write ``metadata.json`` and canonical LD-score parquet files.

        Existing canonical family files are checked before any output file is
        written. Replacement requires ``output_config.overwrite=True``;
        unrelated files in the directory are ignored. The returned paths map
        includes ``"baseline"``, ``"metadata"``, ``"query"`` when query
        annotations were supplied, and ``"overlap"`` for partitioned runs.
        """
        output_dir = ensure_output_directory(output_config.output_dir, label="LD-score output directory")
        baseline_table = getattr(result, "baseline_table", None)
        query_table = getattr(result, "query_table", None)
        if baseline_table is None:
            raise LDSCInternalError(
                "LD-score output writer could not write the result because LDScoreDirectoryWriter.write() "
                "received an LDScoreResult without baseline_table. Most likely an upstream LD-score workflow "
                "returned an incomplete result object. Re-run with DEBUG logging and report the traceback."
            )
        snp_identifier = validate_ldscore_baseline_table(
            baseline_table, result.baseline_columns, result.config_snapshot,
        )
        validate_ldscore_query_table(
            baseline_table, query_table, result.query_columns, snp_identifier=snp_identifier,
        )

        overlap = getattr(result, "overlap", None)
        family = self.artifact_family(output_dir, result)
        paths = family.paths
        identity_drops_by_chrom = dict(getattr(result, "identity_drops_by_chrom", {}) or {})
        stale_paths = family.preflight(overwrite=output_config.overwrite)

        compression = None if output_config.parquet_compression in {None, "none"} else output_config.parquet_compression
        baseline_rg = _write_chromosome_aligned_parquet(baseline_table, paths["baseline"], compression)
        query_rg = None
        if query_table is not None:
            query_rg = _write_chromosome_aligned_parquet(query_table, paths["query"], compression)
        if overlap is not None:
            from .overlap_table import overlap_to_long_frame
            overlap_to_long_frame(overlap).to_parquet(paths["overlap"], index=False)
        for chrom, artifact in identity_drops_by_chrom.items():
            artifact.write_to(paths[f"dropped_snps_chr{chrom}"])
        self._write_query_diagnostic_files(result, paths)
        self._write_ld_window_snp_counts(result, paths)
        self._write_reference_snp_funnel(result, paths)
        metadata = self.build_metadata(
            result,
            files=family.metadata_files(exclude=self.METADATA_EXCLUDED_FILES),
            baseline_rg=baseline_rg,
            query_rg=query_rg,
        )
        paths["metadata"].write_text(json.dumps(to_serializable(metadata), indent=2, sort_keys=True), encoding="utf-8")
        remove_output_artifacts(stale_paths)
        return {name: str(path) for name, path in paths.items()}

    @staticmethod
    def _write_ld_window_snp_counts(result: Any, paths: dict[str, Path]) -> None:
        """Write per-chromosome panel/kept SNP counts behind a parquet ``--ld-wind-snps`` window."""
        if "ld_window_snp_counts" in paths:
            paths["ld_window_snp_counts"].parent.mkdir(parents=True, exist_ok=True)
            result.ld_window_snp_counts.to_csv(paths["ld_window_snp_counts"], sep="\t", index=False)

    @staticmethod
    def _write_reference_snp_funnel(result: Any, paths: dict[str, Path]) -> None:
        """Write every chromosome's reference-SNP funnel, including skipped chromosomes."""
        if "reference_snp_funnel" in paths:
            from .._refpanel.funnel import write_reference_snp_funnel

            write_reference_snp_funnel(result.reference_snp_funnel, paths["reference_snp_funnel"])

    @staticmethod
    def _write_query_diagnostic_files(result: Any, paths: dict[str, Path]) -> None:
        """Serialize fixed query diagnostics after family preflight."""
        query_statuses = tuple(getattr(result, "query_statuses", ()))
        gene_list_batch = getattr(result, "gene_list_batch", None)
        if "chromosome_scope" in paths:
            paths["chromosome_scope"].parent.mkdir(parents=True, exist_ok=True)
            atomic_write_json(getattr(result, "chromosome_scope", None) or result.source_summary["chromosome_scope"], paths["chromosome_scope"])
        if "input_issues" in paths:
            paths["input_issues"].parent.mkdir(parents=True, exist_ok=True)
            result.input_issues.to_csv(paths["input_issues"], sep="\t", index=False, na_rep="")
        if "query_status" in paths:
            paths["query_status"].parent.mkdir(parents=True, exist_ok=True)
            pd.DataFrame([record.as_dict() for record in query_statuses], columns=QUERY_STATUS_COLUMNS).to_csv(
                paths["query_status"], sep="\t", index=False, na_rep=""
            )
        if gene_list_batch is not None:
            LDScoreDirectoryWriter._write_gene_list_diagnostic_files(gene_list_batch, paths)

    @staticmethod
    def _write_gene_list_diagnostic_files(batch: Any, paths: dict[str, Path]) -> None:
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

    def build_metadata(
        self,
        result: Any,
        files: dict[str, str],
        baseline_rg: list[dict] | None = None,
        query_rg: list[dict] | None = None,
    ) -> dict[str, Any]:
        """Build the JSON metadata payload for one LD-score result.

        The metadata always includes ``count_config`` so downstream regression
        code can report the common-SNP count threshold even when MAF metadata is
        unavailable and per-column common-SNP counts are omitted. Row-group metadata
        records the chromosome, row-group index, row offset, and row count for
        each chromosome-aligned parquet row group.
        """
        baseline_table = getattr(result, "baseline_table")
        query_table = getattr(result, "query_table", None)
        config_snapshot = getattr(result, "config_snapshot", None)
        if config_snapshot is None:
            raise LDSCInternalError(
                "LD-score output writer could not build metadata because LDScoreResult.config_snapshot is missing. "
                "Most likely an upstream LD-score workflow returned an older or incomplete result object. "
                "Re-run with DEBUG logging and report the traceback."
            )
        identity_metadata = identity_artifact_metadata(
            artifact_type="ldscore",
            snp_identifier=config_snapshot.snp_identifier,
            genome_build=getattr(config_snapshot, "genome_build", None),
        )
        chromosomes = baseline_table["CHR"].astype(str).drop_duplicates().tolist()
        overlap = getattr(result, "overlap", None)
        overlap_config = None
        if overlap is not None:
            count_config = dict(getattr(result, "count_config", None) or {})
            overlap_config = {
                "total_all_reference_snps": float(overlap.total_all_reference_snps),
                "total_common_reference_snps": (
                    None if overlap.total_common_reference_snps is None
                    else float(overlap.total_common_reference_snps)
                ),
                "common_maf_min": float(count_config.get("common_reference_snp_maf_min", 0.05)),
                "common_maf_operator": str(count_config.get("common_reference_snp_maf_operator", ">=")),
                "stored_block": (
                    "baseline_rows_plus_query_diagonal"
                    if list(getattr(result, "query_columns", []))
                    else "baseline_by_baseline"
                ),
            }
        payload = {
            **identity_metadata,
            "files": dict(files),
            "chromosomes": chromosomes,
            "chromosome_scope": getattr(result, "chromosome_scope", None),
            "baseline_columns": list(getattr(result, "baseline_columns", [])),
            "query_columns": list(getattr(result, "query_columns", [])),
            "annotation_types": (
                dict(getattr(result, "annotation_types", {}) or {})
                or {
                    str(name): "unknown"
                    for name in [
                        *list(getattr(result, "baseline_columns", [])),
                        *list(getattr(result, "query_columns", [])),
                    ]
                }
            ),
            "counts": list(getattr(result, "count_records", [])),
            "count_config": dict(getattr(result, "count_config", None) or DEFAULT_COUNT_CONFIG),
            "overlap_config": overlap_config,
            "n_baseline_rows": int(len(baseline_table)),
            "n_query_rows": int(len(baseline_table)) if getattr(result, "query_columns", []) else 0,
            "row_group_layout": "one_per_chromosome",
            "baseline_row_groups": baseline_rg or [],
            "query_row_groups": query_rg,
            "query_batches": getattr(result, "query_batches", None) or ([{
                "file": files["query"], "query_columns": list(result.query_columns),
                "row_groups": query_rg or [],
            }] if query_table is not None else []),
        }
        query_statuses = tuple(getattr(result, "query_statuses", ()))
        gene_list_batch = getattr(result, "gene_list_batch", None)
        if query_statuses:
            payload["query_diagnostics"] = {"status": "diagnostics/query_annotation_status.tsv"}
        if getattr(result, "chromosome_scope", None):
            payload.setdefault("query_diagnostics", {})["chromosome_scope"] = "diagnostics/chromosome_scope.json"
        input_issues = getattr(result, "input_issues", None)
        if input_issues is not None and not input_issues.empty:
            payload.setdefault("query_diagnostics", {})["input_issues"] = "diagnostics/input_issues.tsv"
        if gene_list_batch is not None:
            payload["gene_list_resolution_policy"] = gene_list_batch.resolution_policy
            catalog = getattr(gene_list_batch, "catalog", None)
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
        index_provenance = getattr(result, "index_provenance", None)
        if index_provenance is not None:
            payload.update(dict(index_provenance))
        snp_universe_policy = getattr(result, "snp_universe_policy", None)
        if snp_universe_policy is not None:
            payload["snp_universe_policy"] = dict(snp_universe_policy)
        legacy_import = getattr(result, "legacy_ldsc2_import", None)
        if legacy_import is not None:
            payload["legacy_ldsc2_import"] = dict(legacy_import)
        return payload
