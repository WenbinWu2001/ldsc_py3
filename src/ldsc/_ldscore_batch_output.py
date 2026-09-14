"""Publish genome-wide query batches without retaining completed score tables.

Private Parquet files own completed values. Only the shared baseline, compact
counts, overlap blocks, and diagnostic descriptors survive between batches.
"""

from dataclasses import replace
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import json
import os

import pandas as pd

from ._annotation_storage import TsvDiagnostics
from .annotation_semantics import require_unique_annotation_names
from .errors import LDSCInputError
from .outputs import _to_serializable, _write_chromosome_aligned_parquet
from .path_resolution import ensure_output_directory, remove_output_artifacts
from .query_annotations import gene_viability_errors


def write_ldscore_batches(writer, results, output_config, query_status_order=()):
    """Consume one result at a time, then publish the complete artifact family."""
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
                writer._validate_tables(result)
                if shared is None:
                    shared = SimpleNamespace(**vars(result))
                    shared.query_table = None
                    shared.chromosome_results = []
                    shared.query_columns = []
                    shared.count_records = [record for record in result.count_records if record["group"] == "baseline"]
                    shared.annotation_types = {name: result.annotation_types.get(name, "unknown") for name in result.baseline_columns}
                    shared.overlap = None if result.overlap is None else result.overlap.select_queries([])
                    baseline_rg = _write_chromosome_aligned_parquet(shared.baseline_table, stage / "ldscore.baseline.parquet", compression)
                    shared.identity_drops_by_chrom = {}
                    for chrom, artifact in result.identity_drops_by_chrom.items():
                        path = stage / "diagnostics" / "dropped_snps" / f"chr{chrom}_dropped.tsv.gz"
                        artifact.write_to(path)
                        shared.identity_drops_by_chrom[chrom] = TsvDiagnostics(path)
                elif not shared.baseline_table.equals(result.baseline_table):
                    raise LDSCInputError("LD-score execution batches disagree on shared baseline values or SNP rows. Use one fixed annotation/reference universe for the complete run.")
                require_unique_annotation_names(shared.baseline_columns, [*shared.query_columns, *result.query_columns])
                shared.query_columns.extend(result.query_columns)
                shared.count_records.extend(record for record in result.count_records if record["group"] == "query")
                shared.annotation_types.update(result.annotation_types)
                statuses.extend(result.query_statuses)
                if result.query_table is not None:
                    filename = f"query-{len(query_batches) + 1}.parquet"
                    groups = _write_chromosome_aligned_parquet(result.query_table, stage / filename, compression)
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
                from .overlap_matrix import overlap_to_long_frame
                overlap_to_long_frame(shared.overlap).to_parquet(paths["overlap"], index=False)
            writer._write_query_diagnostic_files(shared, paths)
            metadata = writer.build_metadata(shared,
                files=staged_family.metadata_files(exclude=("query_status", "gene_list_audit", "gene_list_resolution_summary", "chromosome_scope", "input_issues")),
                baseline_rg=baseline_rg,
                query_rg=query_batches[0]["row_groups"] if len(query_batches) == 1 else None)
            paths["metadata"].write_text(json.dumps(_to_serializable(metadata), indent=2, sort_keys=True), encoding="utf-8")
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
    from .regression_runner import load_ldscore_from_dir
    from ._gene_query_storage import persistent_gene_diagnostics
    source = load_ldscore_from_dir(str(root))
    paths = {key: str(value) for key, value in family.paths.items()}
    return replace(source, output_paths=paths, query_statuses=shared.query_statuses,
                   gene_list_batch=persistent_gene_diagnostics(shared.gene_list_batch, paths),
                   chromosome_scope=shared.chromosome_scope,
                   snp_universe_policy=getattr(shared, "snp_universe_policy", None),
                   identity_drops_by_chrom={chrom: TsvDiagnostics(root / "diagnostics" / "dropped_snps" / f"chr{chrom}_dropped.tsv.gz")
                                            for chrom in shared.identity_drops_by_chrom})
