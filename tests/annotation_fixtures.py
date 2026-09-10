"""Small, explicitly materialized scientific fixtures around the shard API.

These helpers own test-only temporary output parents and release them when a
fixture is discarded. Production APIs require callers to own their output
directories and never provide these whole-fixture inspection conveniences.
"""

from pathlib import Path
import tempfile
import weakref

import numpy as np
import pandas as pd

from ldsc._annotation_bundle import AnnotationBundle
from ldsc._annotation_storage import AnnotationShard, AnnotationWorkspace, ColumnStore
from ldsc._annotation_parsing import normalize_annotation_chunk
from ldsc._kernel.snp_identity import effective_merge_key_series, identity_base_mode, is_allele_aware_mode


def _own(bundle, temporary):
    workspace = bundle.workspace
    def cleanup():
        workspace.close()
        temporary.cleanup()
    weakref.finalize(bundle, cleanup)
    return bundle


def build_annotation_fixture(builder, source_spec=None):
    """Run the real source preparer in an independently owned test directory."""
    temporary = tempfile.TemporaryDirectory(prefix="ldsc-annotation-test-")
    try:
        return _own(builder.run(source_spec, output_dir=temporary.name), temporary)
    except BaseException:
        temporary.cleanup()
        raise


def make_annotation_bundle(*, metadata, baseline_annotations, query_annotations,
                           baseline_columns, query_columns, chromosomes, source_summary,
                           config_snapshot=None, query_statuses=(), gene_list_batch=None):
    """Write small independently specified arrays directly into test shards."""
    temporary = tempfile.TemporaryDirectory(prefix="ldsc-annotation-test-")
    workspace = AnnotationWorkspace(temporary.name)
    shards = {}
    for chrom in chromosomes:
        rows = np.flatnonzero(metadata.CHR.astype(str).eq(str(chrom)))
        path = workspace.path / str(chrom)
        path.mkdir()
        metadata.iloc[rows].reset_index(drop=True).to_parquet(path / "metadata.parquet", index=False)
        stores = []
        for role, frame in (("baseline", baseline_annotations), ("query", query_annotations)):
            store = ColumnStore.create(path / f"{role}.npy", len(rows), frame.columns)
            store.write(0, frame.iloc[rows].to_numpy(dtype=np.float32))
            stores.append(store)
        shards[str(chrom)] = AnnotationShard(path / "metadata.parquet", tuple(stores), len(rows))
    return _own(AnnotationBundle(shards, list(baseline_columns), list(query_columns), workspace,
                                source_summary=source_summary, config_snapshot=config_snapshot,
                                query_statuses=query_statuses, gene_list_batch=gene_list_batch), temporary)


def fixture_metadata(bundle):
    """Explicitly inspect all rows of a small test fixture, never cache them."""
    if not isinstance(bundle, AnnotationBundle):
        return bundle.metadata
    return pd.concat([bundle.metadata_for_chromosome(chrom) for chrom in bundle.chromosomes], ignore_index=True)


def fixture_values(bundle, role="all", *, include_query=True):
    """Materialize only test-oracle values in declared chromosome/column order."""
    columns = bundle.baseline_columns if role == "baseline" or not include_query else bundle.query_columns if role == "query" else bundle.baseline_columns + bundle.query_columns
    return pd.concat([pd.DataFrame(bundle.read(chrom, columns=columns), columns=columns) for chrom in bundle.chromosomes], ignore_index=True)


def fixture_ids(bundle, mode):
    metadata = fixture_metadata(bundle)
    if is_allele_aware_mode(mode) and not {"A1", "A2"}.issubset(metadata):
        mode = identity_base_mode(mode)
    return set(effective_merge_key_series(metadata, mode).astype(str))


def read_annotation_fixture(path, mode):
    """Exercise normalization on one deliberately small parser fixture."""
    return normalize_annotation_chunk(pd.read_csv(path, sep=r"\s+"), path, mode)


def replace_fixture_queries(bundle, *, query_annotations, query_columns=None):
    """Construct a new fixture with independently supplied query values."""
    return make_annotation_bundle(metadata=fixture_metadata(bundle), baseline_annotations=fixture_values(bundle, "baseline"),
                                  query_annotations=query_annotations, baseline_columns=bundle.baseline_columns,
                                  query_columns=list(query_annotations.columns) if query_columns is None else query_columns,
                                  chromosomes=bundle.chromosomes, source_summary=bundle.source_summary,
                                  config_snapshot=bundle.config_snapshot)


def resolve_gene_fixture(focal_paths, catalog, *, output_dir, **kwargs):
    """Stage gene fixtures under the test framework's owned output parent."""
    from ldsc._gene_query_storage import resolve_gene_lists_staged
    return resolve_gene_lists_staged(focal_paths, catalog, AnnotationWorkspace(output_dir), **kwargs)


def fixture_gene_audit(batch):
    """Explicitly inspect the complete small audit used by an independent test."""
    from ldsc.gene_list_resolver import AUDIT_COLUMNS
    return pd.concat(batch.audit_frames(), ignore_index=True) if any(spool.n_rows for spool in batch.audit_spools) else pd.DataFrame(columns=AUDIT_COLUMNS)
