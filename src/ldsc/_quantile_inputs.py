"""Reconstruct quantile universes and accumulate fitted-annotation statistics.

Only the eligible float64 target vector is genome-wide. Exact boundaries are
selected once, that vector is released, and all annotation statistics use
bounded reads from chromosome artifacts. Chan's merged centered moments avoid
the cancellation of E[x²] - E[x]² for high-offset annotations.
"""

from dataclasses import dataclass
import json

import numpy as np
import pandas as pd

from ._annotation_loading import build_annotation_shards
from ._kernel.snp_identity import effective_merge_key_series, identity_base_mode, is_allele_aware_mode
from ._quantile_storage import QuantileIdentityStore, reference_chunks, target_chunks
from .config import AnnotationBuildConfig, GlobalConfig
from .errors import LDSCInputError
from .path_resolution import ANNOTATION_SUFFIXES, resolve_file_group, split_cli_path_tokens


@dataclass(frozen=True)
class QuantileInputStatistics:
    """Compact validated sufficient statistics and original-source provenance."""

    annotation_sums: np.ndarray
    full_sums: np.ndarray
    annotation_sd: np.ndarray
    lower: np.ndarray
    upper: np.ndarray
    counts: np.ndarray
    n_common: int
    n_excluded: int
    ldscore_metadata: dict
    reference_paths: list
    target_paths: list
    common_maf_min: float


def _model_annotations(args, model, metadata, workspace):
    spec = AnnotationBuildConfig(
        **{name: tuple(split_cli_path_tokens(getattr(args, name, None))) for name in (
            "baseline_annot_sources", "query_annot_sources", "query_annot_bed_sources", "query_annot_gene_list_sources")},
        **{name: getattr(args, name, default) for name, default in (
            ("gene_coordinate_file", None), ("control_gene_list_file", None),
            ("gene_list_resolution_policy", "strict"), ("gene_exclude_regions", "none"), ("padding_bp", 0))},
    )
    mode, build = str(metadata.get("snp_identifier")), metadata.get("genome_build")
    bundle = build_annotation_shards(spec, GlobalConfig(snp_identifier=mode, genome_build=build), workspace,
                                     projection_genome_build=build)
    missing = [name for name in model.annotation_names if name not in bundle.baseline_columns + bundle.query_columns]
    if missing:
        raise LDSCInputError(f"quantile-h2 resupplied annotation sources are missing fitted annotation(s): {missing}. Resupply every original baseline/query source for this fitted model.")
    return bundle


def _infer_keys(store, frame, mode, role, diagnostics, reference_source):
    if not is_allele_aware_mode(mode) or {"A1", "A2"}.issubset(frame.columns):
        return effective_merge_key_series(frame, mode, context=role).astype(str)
    base_mode = identity_base_mode(mode)
    if not getattr(store, "inference_validated", False):
        if store.duplicates(diagnostics, "reference", reference_source, base=True):
            raise LDSCInputError(f"quantile-h2 cannot infer omitted alleles for {role}: reference metadata has non-unique {base_mode} identities. Supply A1/A2 in the annotation source.")
        store.inference_validated = True
    bases = effective_merge_key_series(frame, base_mode, context=role).astype(str)
    lookup = {}
    for base in pd.unique(bases):
        found = store.db.execute("SELECT effective_snp_id FROM reference WHERE base_key=?", (base,)).fetchone()
        lookup[base] = found[0] if found else f"__unresolved_reference_metadata__:{role}:{base}"
    return bases.map(lookup)


def _alignment(args, bundle, metadata, store, diagnostics, reference_paths, target_paths):
    mode = str(metadata.get("snp_identifier"))
    reference_source, target_source = ",".join(reference_paths), ",".join(target_paths)
    invalid_maf = 0
    for frame in reference_chunks(reference_paths):
        bad = ~np.isfinite(frame.MAF)
        invalid_maf += int(bad.sum())
        diagnostics.append(frame.loc[bad], role="reference_metadata", source=reference_source,
                           issue="missing_maf", action="fatal", details="MAF must be finite for every intersected reference SNP")
        store.append("reference", frame)
    if invalid_maf:
        raise LDSCInputError("quantile-h2 reference metadata contains missing or nonfinite MAF values.")
    for frame in store.frames("SELECT * FROM reference ORDER BY ordinal"):
        keys = effective_merge_key_series(frame, mode, context="reference metadata").astype(str)
        bases = effective_merge_key_series(frame, identity_base_mode(mode), context="reference metadata base identity").astype(str)
        store.db.executemany("UPDATE reference SET effective_snp_id=?, base_key=? WHERE ordinal=?",
                             zip(keys, bases, map(int, frame.ordinal)))
    store.db.commit()
    if store.duplicates(diagnostics, "reference", reference_source):
        raise LDSCInputError("quantile-h2 reference metadata has duplicate effective SNP identities.")
    for chrom in bundle.chromosomes:
        frame = bundle.metadata_for_chromosome(chrom)
        frame["effective_snp_id"] = _infer_keys(store, frame, mode, "fitted annotations", diagnostics, reference_source)
        frame["shard"], frame["annotation_row"] = chrom, np.arange(len(frame))
        store.append("fitted", frame.drop(columns=[name for name in frame if name not in {
            "CHR", "POS", "SNP", "A1", "A2", "effective_snp_id", "shard", "annotation_row"}]))
        del frame
    store.issues(diagnostics, "fitted", "NOT EXISTS (SELECT 1 FROM reference s WHERE s.effective_snp_id=r.effective_snp_id)",
                 role="fitted_annotations", source="resupplied annotation sources", issue="missing_reference_metadata",
                 action="excluded", details="annotation SNP is absent from reference metadata")
    store.issues(diagnostics, "reference", "NOT EXISTS (SELECT 1 FROM fitted s WHERE s.effective_snp_id=r.effective_snp_id)",
                 role="reference_metadata", source=reference_source, issue="missing_baseline_annotation",
                 action="excluded", details="reference SNP is outside the fitted annotation grid")
    threshold = float((metadata.get("count_config") or {}).get("common_reference_snp_maf_min", .05))
    operator = str((metadata.get("count_config") or {}).get("common_reference_snp_maf_operator", ">="))
    if operator != ">=":
        raise LDSCInputError(f"quantile-h2 does not support inherited common-MAF operator {operator!r}.")
    intersection = "EXISTS (SELECT 1 FROM fitted f WHERE f.effective_snp_id=r.effective_snp_id)"
    n_all = store.scalar(f"SELECT count(*) FROM reference r WHERE {intersection}")
    _validate_universe(metadata, "all", n_all)
    for frame in target_chunks(target_paths, args.target_annotation, getattr(args, "target_missing_value", None)):
        frame["effective_snp_id"] = _infer_keys(store, frame, mode, "target annotation", diagnostics, reference_source)
        store.append("target", frame)
    if store.duplicates(diagnostics, "target", target_source):
        raise LDSCInputError("quantile-h2 target annotation has duplicate effective SNP identities.")
    store.issues(diagnostics, "target", "NOT EXISTS (SELECT 1 FROM reference s WHERE s.effective_snp_id=r.effective_snp_id)",
                 role="target_annotation", source=target_source, issue="missing_reference_metadata",
                 action="excluded", details="target SNP is absent from reference metadata")
    missing = store.issues(diagnostics, "reference",
                          f"MAF>=? AND {intersection} AND NOT EXISTS (SELECT 1 FROM target t WHERE t.effective_snp_id=r.effective_snp_id)",
                          params=(threshold,), role="target_annotation", source=target_source, issue="missing_target_annotation",
                          action="fatal", details="common reference-SNP universe member is absent from target source")
    if missing:
        raise LDSCInputError("quantile-h2 target annotation does not cover the common reference-SNP universe.")
    n_common = store.scalar(f"SELECT count(*) FROM reference r WHERE MAF>=? AND {intersection}", (threshold,))
    _validate_universe(metadata, "common", n_common)
    return threshold, n_common


def _validate_universe(metadata, universe, actual):
    expected = (metadata.get("overlap_config") or {}).get(f"total_{universe}_reference_snps")
    if expected is not None and int(round(float(expected))) != actual:
        label = "all-reference SNP" if universe == "all" else "common reference-SNP"
        raise LDSCInputError(f"quantile-h2 reconstructed a different {label} universe size than the fitted LD-score artifact: {actual} vs {expected}. Check --ref-metadata-sources and annotation sources.")


def _validate_aggregates(model, metadata, sums, products):
    records = {str(row.get("column")): row for row in metadata.get("counts", [])}
    for name, actual in zip(model.annotation_names, sums):
        expected = records.get(name, {}).get("common_reference_snp_count")
        if expected is not None and not np.isclose(actual, float(expected), rtol=1e-6, atol=1e-8):
            raise LDSCInputError(f"quantile-h2 fitted annotation '{name}' has common reference-SNP universe sum {actual}, expected {expected}. Resupply the exact annotation sources used for LD-score construction.")
    overlap_path = (metadata.get("files") or {}).get("overlap")
    if overlap_path:
        from .overlap_matrix import assemble_model_overlap, overlap_from_long_frame

        config = metadata.get("overlap_config") or {}
        overlap = overlap_from_long_frame(pd.read_parquet(model.ldscore_dir / str(overlap_path)),
            baseline_columns=metadata.get("baseline_columns", []), query_columns=metadata.get("query_columns", []),
            total_all_reference_snps=config.get("total_all_reference_snps"), total_common_reference_snps=config.get("total_common_reference_snps"))
        expected = assemble_model_overlap(overlap, model.annotation_names, use_common=True)
        if not np.allclose(products, expected, rtol=1e-6, atol=1e-8):
            raise LDSCInputError("quantile-h2 reconstructed fitted-annotation cross-products that disagree with the linked LD-score overlap artifact. Resupply the exact original annotations and reference metadata.")


def prepare_quantile_statistics(args, model, workspace, diagnostics):
    """Validate sources and reduce chromosome artifacts using exact global bounds."""
    metadata = json.loads((model.ldscore_dir / "metadata.json").read_text(encoding="utf-8"))
    bundle = _model_annotations(args, model, metadata, workspace)
    references = resolve_file_group(split_cli_path_tokens(args.ref_metadata_sources),
        suffixes=("", ".tsv", ".tsv.gz", ".txt", ".txt.gz"), label="reference metadata", allow_chromosome_suite=True)
    targets = resolve_file_group(split_cli_path_tokens(args.target_annot_sources), suffixes=ANNOTATION_SUFFIXES,
        label="target annotation", allow_chromosome_suite=True)
    with QuantileIdentityStore(workspace.path / "quantile-identities.sqlite") as store:
        threshold, n_common = _alignment(args, bundle, metadata, store, diagnostics, references, targets)
        joins = "FROM fitted f JOIN reference r ON r.effective_snp_id=f.effective_snp_id JOIN target t ON t.effective_snp_id=r.effective_snp_id WHERE r.MAF>=?"
        n_eligible = store.scalar(f"SELECT count(*) {joins} AND NOT t.target_excluded", (threshold,))
        target_vector = np.empty(n_eligible, dtype=np.float64)
        offset = 0
        p = len(model.annotation_names)
        full_sums, mean, m2 = (np.zeros(p) for _ in range(3))
        products, seen = np.zeros((p, p)), 0
        same_name_mismatch = False
        target_index = model.annotation_names.index(args.target_annotation) if args.target_annotation in model.annotation_names else None
        for values, aligned in _common_values(bundle, model.annotation_names, store, joins, threshold):
            n = len(values)
            local_mean = values.mean(axis=0)
            centered = values - local_mean
            delta = local_mean - mean
            m2 += np.einsum("ij,ij->j", centered, centered) + delta * delta * (seen * n / (seen + n))
            mean += delta * (n / (seen + n))
            seen += n
            full_sums += values.sum(axis=0, dtype=np.float64)
            products += values.T @ values
            eligible = ~aligned.target_excluded.to_numpy(dtype=bool)
            target_values = aligned.target_value.to_numpy(dtype=np.float64)[eligible]
            target_vector[offset:offset + len(target_values)] = target_values
            offset += len(target_values)
            if target_index is not None:
                same_name_mismatch |= not np.array_equal(values[eligible, target_index].astype(np.float32), target_values.astype(np.float32))
            del values, centered, aligned, eligible, target_values
        _validate_aggregates(model, metadata, full_sums, products)
        if same_name_mismatch:
            raise LDSCInputError(f"Target annotation '{args.target_annotation}' matches a fitted annotation name but its values differ. Supply the same values or rename the external target annotation.")
        q = args.num_quantiles
        if int(q) != q or q < 2:
            raise LDSCInputError("quantile-h2 requires --num-quantiles to be an integer of at least 2.")
        if not n_eligible:
            raise LDSCInputError("quantile-h2 target values must be nonempty, numeric, and finite after exclusions.")
        target_vector.sort(kind="stable")
        indices = np.floor(np.arange(q + 1) * (n_eligible - 1) / q + .5).astype(int)
        bounds = target_vector[indices].copy()
        del target_vector
        quantile_sums = np.zeros((p, q))
        counts = np.zeros(q, dtype=np.int64)
        for values, aligned in _common_values(bundle, model.annotation_names, store, joins, threshold):
            eligible = ~aligned.target_excluded.to_numpy(dtype=bool)
            target_values = aligned.target_value.to_numpy(dtype=np.float64)[eligible]
            labels = np.searchsorted(bounds[1:], target_values, side="left")
            counts += np.bincount(labels, minlength=q)
            # Add indexed rows directly; no dense SNP-by-quantile indicator.
            for column in range(p):
                np.add.at(quantile_sums[column], labels, values[eligible, column])
            del values, aligned, eligible, target_values, labels
        if (counts == 0).any():
            empty = (np.flatnonzero(counts == 0) + 1).tolist()
            raise LDSCInputError(f"quantile-h2 produced an empty quantile: {empty}. Target-value ties do not provide enough distinct realized intervals; reduce --num-quantiles or use a less discrete target annotation.")
    return QuantileInputStatistics(quantile_sums, full_sums, np.sqrt(m2 / n_common), bounds[:-1], bounds[1:], counts,
                                  n_common, n_common - n_eligible, metadata, references, targets, threshold)


def _common_values(bundle, names, store, joins, threshold):
    """Yield one bounded aligned tile, releasing it before the next read."""
    for chrom in bundle.chromosomes:
        for aligned in store.frames(f"SELECT f.annotation_row,t.target_value,t.target_excluded {joins} AND f.shard=? ORDER BY f.annotation_row", (threshold, chrom)):
            if aligned.empty:
                continue
            values = bundle.read(chrom, rows=aligned.annotation_row.to_numpy(dtype=np.int64), columns=names).astype(np.float64)
            yield values, aligned
            del values, aligned
