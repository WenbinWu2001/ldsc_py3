"""Reconstruct the quantile SNP sets and accumulate fitted-annotation statistics.

Only the eligible float64 target vector is genome-wide. Exact boundaries are
selected once, that vector is released, and all annotation statistics use
bounded reads from chromosome artifacts. Chan's merged centered moments avoid
the cancellation of E[x^2] - E[x]^2 for high-offset annotations.
"""

from dataclasses import dataclass, field
import json
import logging

import numpy as np
import pandas as pd

from .._annotation.loading import build_annotation_shards
from .._core.snp_identity import effective_merge_key_series, identity_base_mode, is_allele_aware_mode
from .quantile_storage import QuantileIdentityStore, parse_missing_token, reference_chunks, target_chunks
from ..config import GlobalConfig
from ..errors import LDSCInputError
from .._runtime.paths import resolve_file_group, split_cli_path_tokens
from .overlap_model import assemble_model_overlap
from .._ldscore.overlap_table import overlap_from_long_frame
from .._genome.build_policy import GenomeBuildCheck

LOGGER = logging.getLogger("LDSC.quantile_h2")


@dataclass(frozen=True)
class QuantileTarget:
    """The target annotation that defines the quantiles, checked when it is built, before anything is written.

    ``sources`` holds the ``--target-annot-sources`` path tokens, ``column``
    the ``--target-annot-column``, ``num_quantiles`` the ``--num-quantiles``,
    and ``missing_value`` the ``--target-missing-value`` token as given, which
    the result metadata records; ``missing`` is that token parsed once by
    :func:`parse_missing_token` for the target reads.

    Raises
    ------
    LDSCInputError
        If ``num_quantiles`` is not an integer of at least 2.
    """

    sources: object
    column: str
    num_quantiles: int
    missing_value: object = None
    missing: object = field(init=False)

    def __post_init__(self):
        q = self.num_quantiles
        if int(q) != q or q < 2:
            raise LDSCInputError("quantile-h2 requires --num-quantiles to be an integer of at least 2.")
        object.__setattr__(self, "missing", parse_missing_token(self.missing_value))


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


def _model_annotations(annotation_config, model, metadata, workspace, issues_path):
    mode, build = str(metadata.get("snp_identifier")), metadata.get("genome_build")
    bundle = build_annotation_shards(annotation_config, GlobalConfig(snp_identifier=mode, genome_build=build), workspace,
                                     projection_genome_build=build, issues_path=issues_path)
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


def _alignment(target, bundle, metadata, store, diagnostics, reference_paths, target_paths):
    mode = str(metadata.get("snp_identifier"))
    reference_source, target_source = ",".join(reference_paths), ",".join(target_paths)
    invalid_maf = 0
    for path in reference_paths:
        for frame in reference_chunks([path]):
            if is_allele_aware_mode(mode) and not {"A1", "A2"}.issubset(frame.columns):
                # The reference metadata defines the SNP universe, so it cannot fall back to base identity.
                raise LDSCInputError(
                    f"quantile-h2 reference metadata '{path}' does not have both A1 and A2, but the fit is "
                    f"allele-aware (snp_identifier {mode}), so each reference SNP is identified by its alleles. "
                    "Set --ref-metadata-sources to the R2 matrix sidecars (chr<N>_meta.tsv.gz) or the files "
                    "that ldscore --export-ref-metadata wrote for the panel of the fit; both carry A1 and A2."
                )
            bad = ~np.isfinite(frame.MAF)
            invalid_maf += int(bad.sum())
            diagnostics.append(frame.loc[bad], role="reference_metadata", source=reference_source, issue="missing_maf",
                               action="fatal", details="MAF must be finite for every intersected reference SNP")
            store.append("reference", frame)
    if invalid_maf:
        raise LDSCInputError(
            f"quantile-h2 reference metadata contains missing or nonfinite MAF values in {invalid_maf} SNP rows. "
            "Restore the original SNP frequencies from the panel used for the fit and set "
            "--ref-metadata-sources to that matching metadata; do not replace unknown frequencies with zero. "
            f"See {diagnostics.path} for the affected SNP rows."
        )
    for frame in store.frames("SELECT * FROM reference ORDER BY ordinal"):
        keys = effective_merge_key_series(frame, mode, context="reference metadata").astype(str)
        bases = effective_merge_key_series(frame, identity_base_mode(mode), context="reference metadata base identity").astype(str)
        store.db.executemany("UPDATE reference SET effective_snp_id=?, base_key=? WHERE ordinal=?",
                             zip(keys, bases, map(int, frame.ordinal)))
    store.db.commit()
    if store.duplicates(diagnostics, "reference", reference_source):
        raise LDSCInputError(
            "quantile-h2 reference metadata has duplicate effective SNP identities."
            " The selected SNP identifier mode cannot distinguish these rows. Set --ref-metadata-sources "
            "to the unique SNP metadata from the original fitted panel; do not remove ambiguous rows "
            "just to bypass this check. "
            f"See {diagnostics.path} for the duplicate identities."
        )
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
                 action="excluded", details="reference SNP is outside the fitted annotation SNP grid")
    threshold = float((metadata.get("count_config") or {}).get("common_reference_snp_maf_min", .05))
    operator = str((metadata.get("count_config") or {}).get("common_reference_snp_maf_operator", ">="))
    if operator != ">=":
        raise LDSCInputError(
            f"quantile-h2 does not support inherited common-MAF operator {operator!r}."
            " This workflow requires the native inclusive MAF threshold (>=). Use a fitted result "
            "with that count policy; changing the SNP population requires regenerating LD scores "
            "and refitting partitioned-h2, not editing the saved operator. See "
            "docs/user-guide/troubleshooting.md#quantile-h2-fitted-model-or-resupplied-sources-do-not-verify."
        )
    intersection = "EXISTS (SELECT 1 FROM fitted f WHERE f.effective_snp_id=r.effective_snp_id)"
    n_all = store.scalar(f"SELECT count(*) FROM reference r WHERE {intersection}")
    _validate_universe(metadata, "all", n_all)
    for path in target_paths:
        for chunk, frame in enumerate(target_chunks([path], target.column, target.missing)):
            if chunk == 0 and is_allele_aware_mode(mode) and not {"A1", "A2"}.issubset(frame.columns):
                LOGGER.info(f"quantile-h2 target annotation '{path}' has no A1/A2 columns, so its SNPs are matched "
                            f"by {identity_base_mode(mode)}, the base mode of {mode}, without checking alleles.")
            frame["effective_snp_id"] = _infer_keys(store, frame, mode, "target annotation", diagnostics, reference_source)
            store.append("target", frame)
    if store.duplicates(diagnostics, "target", target_source):
        raise LDSCInputError(
            "quantile-h2 target annotation has duplicate effective SNP identities."
            " The selected SNP identifier mode cannot distinguish these rows. Supply the intended target "
            "source with one value for each fitted SNP identity, with the original panel's alleles "
            "where required. "
            f"See {diagnostics.path} for the duplicate identities."
        )
    store.issues(diagnostics, "target", "NOT EXISTS (SELECT 1 FROM reference s WHERE s.effective_snp_id=r.effective_snp_id)",
                 role="target_annotation", source=target_source, issue="missing_reference_metadata",
                 action="excluded", details="target SNP is absent from reference metadata")
    missing = store.issues(diagnostics, "reference",
                          f"MAF>=? AND {intersection} AND NOT EXISTS (SELECT 1 FROM target t WHERE t.effective_snp_id=r.effective_snp_id)",
                          params=(threshold,), role="target_annotation", source=target_source, issue="missing_target_annotation",
                          action="fatal", details="common reference SNP is absent from target source")
    if missing:
        raise LDSCInputError(
            "quantile-h2 target annotation does not cover the common reference SNPs; "
            f"missing SNP rows={missing}. "
            "Supply target values for the full common-SNP population used for this fit and check SNP "
            "identifiers and alleles; do not fill absent SNP values with zero. "
            f"See {diagnostics.path} for the missing SNP identities. "
            "Other causes & fixes: "
            "docs/user-guide/troubleshooting.md#quantile-h2-fitted-model-or-resupplied-sources-do-not-verify."
        )
    n_common = store.scalar(f"SELECT count(*) FROM reference r WHERE MAF>=? AND {intersection}", (threshold,))
    _validate_universe(metadata, "common", n_common)
    return threshold, n_common


def _validate_universe(metadata, universe, actual):
    expected = (metadata.get("overlap_config") or {}).get(f"total_{universe}_reference_snps")
    if expected is not None and int(round(float(expected))) != actual:
        label = "all-reference SNP" if universe == "all" else "common reference-SNP"
        raise LDSCInputError(f"quantile-h2 reconstructed a different {label} count than the fitted LD-score directory: {actual} vs {expected}. Check --ref-metadata-sources and annotation sources.")


def _validate_aggregates(model, metadata, sums, products):
    records = {str(row.get("column")): row for row in metadata.get("counts", [])}
    for name, actual in zip(model.annotation_names, sums):
        expected = records.get(name, {}).get("common_reference_snp_count")
        if expected is not None and not np.isclose(actual, float(expected), rtol=1e-6, atol=1e-8):
            raise LDSCInputError(f"quantile-h2 fitted annotation '{name}' has common reference-SNP sum {actual}, expected {expected}. Resupply the exact annotation sources used for LD-score construction.")
    overlap_path = (metadata.get("files") or {}).get("overlap")
    if overlap_path:
        config = metadata.get("overlap_config") or {}
        overlap = overlap_from_long_frame(pd.read_parquet(model.ldscore_dir / str(overlap_path)),
            baseline_columns=metadata.get("baseline_columns", []), query_columns=metadata.get("query_columns", []),
            total_all_reference_snps=config.get("total_all_reference_snps"), total_common_reference_snps=config.get("total_common_reference_snps"))
        expected = assemble_model_overlap(overlap, model.annotation_names, use_common=True)
        if not np.allclose(products, expected, rtol=1e-6, atol=1e-8):
            raise LDSCInputError("quantile-h2 reconstructed fitted-annotation cross-products that disagree with the linked LD-score overlap artifact. Resupply the exact original annotations and reference metadata.")


def prepare_quantile_statistics(annotation_config, target, ref_metadata_sources, model, workspace, diagnostics,
                                input_issues_path):
    """Validate sources and reduce chromosome artifacts using exact global bounds.

    ``annotation_config`` holds the resupplied annotation sources of the fit,
    ``target`` the checked :class:`QuantileTarget`, and
    ``ref_metadata_sources`` the reference SNP metadata path tokens. Unusable
    BED query names are written to ``input_issues_path`` before the run stops.
    """
    metadata = json.loads((model.ldscore_dir / "metadata.json").read_text(encoding="utf-8"))
    if str(metadata.get('snp_identifier', '')).startswith('chr_pos'):
        check = GenomeBuildCheck('auto', context='quantile-h2', option=None)
        check.declared('saved LD scores', model.ldscore_dir, metadata.get('genome_build'))
        for flag, tokens in (('--baseline-annot-sources', annotation_config.baseline_annot_sources),
                             ('--query-annot-sources', annotation_config.query_annot_sources),
                             ('--target-annot-sources', target.sources),
                             ('--ref-metadata-sources', ref_metadata_sources),
                             ('--gene-coordinate-file', annotation_config.gene_coordinate_file)):
            check.paths(flag, tokens)
        check.unchecked('--query-annot-bed-sources', annotation_config.query_annot_bed_sources,
                        'ordinary BED; inherits saved build')
        check.resolve()
    bundle = _model_annotations(annotation_config, model, metadata, workspace, input_issues_path)
    references = resolve_file_group(split_cli_path_tokens(ref_metadata_sources),
        label="reference metadata", allow_chromosome_suite=True)
    targets = resolve_file_group(split_cli_path_tokens(target.sources),
        label="target annotation", allow_chromosome_suite=True)
    with QuantileIdentityStore(workspace.path / "quantile-identities.sqlite") as store:
        threshold, n_common = _alignment(target, bundle, metadata, store, diagnostics, references, targets)
        joins = "FROM fitted f JOIN reference r ON r.effective_snp_id=f.effective_snp_id JOIN target t ON t.effective_snp_id=r.effective_snp_id WHERE r.MAF>=?"
        n_eligible = store.scalar(f"SELECT count(*) {joins} AND NOT t.target_excluded", (threshold,))
        target_vector = np.empty(n_eligible, dtype=np.float64)
        offset = 0
        p = len(model.annotation_names)
        full_sums, mean, m2 = (np.zeros(p) for _ in range(3))
        products, seen = np.zeros((p, p)), 0
        same_name_mismatch = False
        target_index = model.annotation_names.index(target.column) if target.column in model.annotation_names else None
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
            raise LDSCInputError(f"Target annotation '{target.column}' matches a fitted annotation name but its values differ. Supply the same values or rename the external target annotation.")
        q = target.num_quantiles
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
            eligible_values = values[eligible]
            for column in range(p):
                np.add.at(quantile_sums[column], labels, eligible_values[:, column])
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
