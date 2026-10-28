"""Publish genome-wide query batches without retaining completed score tables.

Private Parquet files own completed values. Only the shared baseline, compact
counts, overlap blocks, and diagnostic descriptors survive between batches.
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import TYPE_CHECKING
import json
import os

import pandas as pd

from .._annotation.storage import TsvDiagnostics
from .._annotation.semantics import require_unique_annotation_names
from ..errors import LDSCInputError
from .source import validate_ldscore_baseline_table, validate_ldscore_query_table
from .._runtime.artifacts import to_serializable
from .outputs import _write_chromosome_aligned_parquet
from .._runtime.paths import ensure_output_directory, remove_output_artifacts
from .._annotation.query_outcomes import gene_viability_errors

if TYPE_CHECKING:
    from ..config import GlobalConfig
    from .._annotation.gene_query_storage import StagedGeneListBatch
    from .overlap_table import LDScoreOverlap
    from .._annotation.query_outcomes import QueryAnnotationStatus


@dataclass(slots=True)
class _LDScorePublication:
    """Retain only shared artifacts and metadata while consuming query batches.

    Query values, producer SNP sets, and chromosome results remain owned by
    the active batch. Diagnostic paths here refer to publication staging.
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
    legacy_ldsc2_import: dict | None
    overlap: LDScoreOverlap | None
    ld_window_snp_counts: pd.DataFrame | None = None
    reference_snp_funnel: pd.DataFrame | None = None
    query_columns: list[str] = field(default_factory=list)
    query_statuses: tuple[QueryAnnotationStatus, ...] = ()
    query_batches: list[dict] = field(default_factory=list)
    identity_drops_by_chrom: dict[str, TsvDiagnostics] = field(default_factory=dict)


def write_ldscore_batches(writer, results, output_config, query_status_order=()):
    """Consume one result at a time, then publish the complete artifact family."""
    from .fragments import StagedLDScoreBatch, write_query_fragments

    root = ensure_output_directory(output_config.output_dir, label="LD-score output directory")
    writer.artifact_family(root).preflight(overwrite=output_config.overwrite)
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
                snp_identifier = validate_ldscore_baseline_table(
                    result.baseline_table, result.baseline_columns, result.config_snapshot,
                )
                if not streamed:
                    validate_ldscore_query_table(
                        result.baseline_table, result.query_table, result.query_columns,
                        snp_identifier=snp_identifier,
                    )
                if shared is None:
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
                        legacy_ldsc2_import=None if streamed else result.legacy_ldsc2_import,
                        overlap=None if result.overlap is None else result.overlap.select_queries([]),
                        ld_window_snp_counts=getattr(result, "ld_window_snp_counts", None),
                        reference_snp_funnel=getattr(result, "reference_snp_funnel", None),
                    )
                    baseline_rg = _write_chromosome_aligned_parquet(shared.baseline_table, stage / "ldscore.baseline.parquet", compression)
                    for chrom, artifact in result.identity_drops_by_chrom.items():
                        path = stage / "diagnostics" / "dropped_snps" / f"chr{chrom}_dropped.tsv.gz"
                        artifact.write_to(path)
                        shared.identity_drops_by_chrom[chrom] = TsvDiagnostics(path)
                elif not shared.baseline_table.equals(result.baseline_table):
                    raise LDSCInputError("LD-score execution batches disagree on shared baseline values or SNP rows. Use the same annotation sources and reference SNPs for the complete run.")
                require_unique_annotation_names(shared.baseline_columns, [*shared.query_columns, *result.query_columns])
                shared.query_columns.extend(result.query_columns)
                shared.count_records.extend(record for record in result.count_records if record["group"] == "query")
                shared.annotation_types.update(result.annotation_types)
                statuses.extend(result.query_statuses)
                if result.query_columns:
                    filename = f"query-{len(query_batches) + 1}.parquet"
                    groups = (write_query_fragments(result, stage / filename, compression) if streamed else
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
                name = "ldscore.query.parquet" if len(query_batches) == 1 else f"ldscore.query.batch{ordinal:05d}.parquet"
                (stage / entry["file"]).rename(stage / name)
                entry["file"] = name
            shared.query_batches = query_batches
            staged_family = writer.artifact_family(stage, shared)
            paths = staged_family.paths
            if shared.overlap is not None:
                from .overlap_table import overlap_to_long_frame
                overlap_to_long_frame(shared.overlap).to_parquet(paths["overlap"], index=False)
            writer._write_query_diagnostic_files(shared, paths)
            writer._write_ld_window_snp_counts(shared, paths)
            writer._write_reference_snp_funnel(shared, paths)
            metadata = writer.build_metadata(shared,
                files=staged_family.metadata_files(exclude=writer.METADATA_EXCLUDED_FILES),
                baseline_rg=baseline_rg,
                query_rg=query_batches[0]["row_groups"] if len(query_batches) == 1 else None)
            paths["metadata"].write_text(json.dumps(to_serializable(metadata), indent=2, sort_keys=True), encoding="utf-8")
            family = writer.artifact_family(root, shared)
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
                    writer.write_query_diagnostics(shared, output_config)
            raise
        finally:
            close = getattr(results, "close", None)
            if close is not None:
                close()
    from .source import load_ldscore_from_dir
    from .._annotation.gene_query_storage import persistent_gene_diagnostics
    source = load_ldscore_from_dir(str(root))
    paths = {key: str(value) for key, value in family.paths.items()}
    return replace(source, output_paths=paths, query_statuses=shared.query_statuses,
                   gene_list_batch=persistent_gene_diagnostics(shared.gene_list_batch, paths),
                   chromosome_scope=shared.chromosome_scope,
                   snp_universe_policy=shared.snp_universe_policy,
                   identity_drops_by_chrom={chrom: TsvDiagnostics(root / "diagnostics" / "dropped_snps" / f"chr{chrom}_dropped.tsv.gz")
                                            for chrom in shared.identity_drops_by_chrom})
