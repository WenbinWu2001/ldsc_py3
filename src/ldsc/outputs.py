"""Output writers for the refactored LDSC workflows.

The public LD-score workflow writes one canonical result directory containing
``metadata.json``, ``ldscore.baseline.parquet``, an optional
``ldscore.query.parquet``, and an optional ``ldscore.overlap.parquet`` (the
annotation overlap matrix required by partitioned-h2's overlap-aware category
tables and read by the shared h2 collinearity guard when there are two or more
LD-score columns; ``rg`` does not use it). The overlap sidecar is
written only when the run has two or more annotation columns; an unpartitioned
single-annotation run (e.g. the synthetic ``base``) omits it because the matrix
would collapse to a SNP count already in ``metadata.json``. Query runs also
write ``diagnostics/query_annotation_status.tsv`` and gene-list runs write a
row-complete audit plus per-source resolution summary. Run identity comes from the chosen
directory name;
output filenames inside that directory are fixed. The parquet payloads are
written with one row group per chromosome and matching metadata so
downstream readers can load a single chromosome without scanning the whole
table. The writer creates missing directories, reuses existing directories, and
refuses existing canonical family files unless the caller explicitly sets
``overwrite=True``; successful overwrites remove stale owned siblings that the
current result did not produce. Root ``metadata.json`` is part of the
downstream contract for LD-score directories. Regression
metadata emitted by this module is diagnostic provenance and is written below
``diagnostics/`` without legacy top-level ``format`` discriminators.

Each directory writer's ``artifact_family`` declares its owned paths and
conditional outputs. Workflows consult that declaration for early collision
checks, adding their own log/audit paths. Writers use the final result to
select produced paths, metadata file entries, and post-write stale cleanup;
workflows do not retain a predicted stale list across computation.

Partitioned-h2 regression summaries use the same directory-oriented output
policy. ``PartitionedH2DirectoryWriter`` writes the aggregate ``partitioned_h2.tsv``
(one ``PARTITIONED_H2_COLUMNS`` schema for both regimes) and can optionally stage a
per-query ``diagnostics/query_annotations`` tree with ``manifest.tsv``, one-row
``partitioned_h2.tsv`` summaries, full ``partitioned_h2_full.tsv`` model
tables, and ``metadata.json`` files before moving it into place. ``RgDirectoryWriter``
uses the same whole-tree staging and replacement policy for optional
``diagnostics/pairs`` details.

Unpartitioned h2 results also carry
``diagnostics/ld_score_regression_bins.tsv``. This table summarizes the exact
SNP population and final fitted model used by h2 and is the authoritative input
for the post-processing binned LD Score regression diagnostic plot. The writer
layer persists the table but never imports a plotting library.
"""

from __future__ import annotations

import gzip
import json
import os
import re
import shutil
import tempfile
import unicodedata
from dataclasses import asdict, dataclass, is_dataclass
from os import PathLike
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

from ._row_alignment import assert_same_snp_rows
from ._kernel.snp_identity import effective_merge_key_series, identity_artifact_metadata, is_allele_aware_mode
from .config import _normalize_required_path
from .errors import LDSCConfigError, LDSCInputError, LDSCInternalError
from .path_resolution import (
    ensure_output_directory,
    preflight_output_artifact_family,
    remove_output_artifacts,
)


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
H2_REGRESSION_BIN_COLUMNS = [
    "bin",
    "n_snps",
    "ld_score_min",
    "ld_score_max",
    "mean_ld_score",
    "mean_chi_square",
    "sd_chi_square",
    "mean_sample_size",
    "mean_fitted_chi_square",
    "mean_regression_weight",
]
QUERY_STATUS_COLUMNS = ["query", "source", "input_type", "status", "reason", "n_annotation_snps", "details"]
DEFAULT_COUNT_CONFIG = {
    "common_reference_snp_maf_min": 0.05,
    "common_reference_snp_maf_operator": ">=",
}
# Single partitioned-h2 schema for both regimes. The functional regime fills it
# with one row per baseline category; the cell-type regime fills it with one row
# per query annotation. Differs only in rows and default sort, never in columns.
PARTITIONED_H2_COLUMNS = [
    "category",
    "prop_snps",
    "prop_h2",
    "prop_h2_se",
    "enrichment",
    "enrichment_se",
    "enrichment_p",
    "coefficient",
    "coefficient_se",
    "coefficient_z",
    "coefficient_p",
    "overlap_annot",
    "total_h2_obs",
    "total_h2_obs_se",
    "total_h2_liab",
    "total_h2_liab_se",
    "category_h2_obs",
    "category_h2_obs_se",
    "category_h2_liab",
    "category_h2_liab_se",
    "samp_prev",
    "pop_prev",
]
QUANTILE_H2_COLUMNS = [
    "quantile", "target_value_lower", "target_value_upper", "n_snps", "prop_snps",
    "h2_obs", "h2_obs_se", "h2_liab", "h2_liab_se", "prop_h2", "prop_h2_se",
    "enrichment", "enrichment_se", "enrichment_p",
]
STANDARDIZED_COEFFICIENT_COLUMNS = [
    "annotation", "annotation_type", "annotation_sd", "tau", "tau_se", "tau_z", "tau_p",
    "tau_star", "tau_star_se", "tau_star_z", "tau_star_p",
]
SNP_ALIGNMENT_ISSUE_COLUMNS = [
    "source_role", "source", "CHR", "POS", "SNP", "A1", "A2", "effective_snp_id",
    "issue", "action", "details",
]
RG_CONCISE_COLUMNS = [
    "trait_1",
    "trait_2",
    "n_snps_used",
    "rg",
    "rg_se",
    "p",
    "note",
]
RG_FULL_COLUMNS = [
    "trait_1",
    "trait_2",
    "n_snps_used",
    "rg",
    "rg_se",
    "z",
    "p",
    "h2_1_obs",
    "h2_1_obs_se",
    "h2_1_liab",
    "h2_1_liab_se",
    "h2_2_obs",
    "h2_2_obs_se",
    "h2_2_liab",
    "h2_2_liab_se",
    "gencov_obs",
    "gencov_obs_se",
    "gencov_liab",
    "gencov_liab_se",
    "samp_prev_1",
    "pop_prev_1",
    "samp_prev_2",
    "pop_prev_2",
    "intercept_h2_1",
    "intercept_h2_1_se",
    "intercept_h2_2",
    "intercept_h2_2_se",
    "intercept_gencov",
    "intercept_gencov_se",
    "ratio_1",
    "ratio_1_se",
    "ratio_2",
    "ratio_2_se",
    "lambda_gc_1",
    "lambda_gc_2",
    "mean_chisq_1",
    "mean_chisq_2",
    "pair_kind",
    "status",
    "error",
]


@dataclass(frozen=True)
class H2OutputConfig:
    """Directory-oriented output config for unpartitioned h2 summaries.

    Parameters
    ----------
    output_dir : str or os.PathLike[str]
        Directory that receives ``h2.tsv``, the exact LD Score regression-bin
        diagnostic, and metadata.
    overwrite : bool, optional
        If ``True``, replace existing fixed h2 outputs. Default is ``False``.
    """

    output_dir: str | PathLike[str]
    overwrite: bool = False

    def __post_init__(self) -> None:
        """Normalize the output directory."""
        object.__setattr__(self, "output_dir", _normalize_required_path(self.output_dir))


@dataclass(frozen=True)
class ArtifactFamily:
    """Resolved writer declaration shared by collision checks and publication.

    Parameters
    ----------
    output_dir : pathlib.Path
        Result root used to make metadata file entries relative.
    paths : dict of str to pathlib.Path
        Outputs selected for this write, keyed by their artifact names.
    owned_paths : tuple of pathlib.Path
        All owned paths, including conditional siblings and derived roots
        eligible for cleanup after a successful replacement.
    label : str
        Workflow-specific description used in collision errors.

    Notes
    -----
    Construction only describes paths. Callers retain their existing write
    order and remove returned stale paths only after successful publication.
    """

    output_dir: Path
    paths: dict[str, Path]
    owned_paths: tuple[Path, ...]
    label: str

    def preflight(
        self, *, overwrite: bool, additional_paths: Iterable[str | PathLike[str]] = (),
    ) -> list[Path]:
        """Check ownership, including workflow-only paths such as the run log.

        ``overwrite=False`` raises ``FileExistsError`` for any existing owned
        path. With overwrite enabled, return existing owned paths that this
        declaration does not produce. ``additional_paths`` are protected
        workflow outputs; their writer controls their cleanup separately.

        Early callers discard the return value because production is not yet
        known. Writers call again with the final declaration and use that
        stale list after writing. Additional paths never enter result metadata.
        """
        additional_paths = tuple(additional_paths)
        return preflight_output_artifact_family(
            [*self.paths.values(), *additional_paths],
            [*self.owned_paths, *additional_paths],
            overwrite=overwrite, label=self.label,
        )

    def metadata_files(self, *, exclude: Iterable[str] = ()) -> dict[str, str]:
        """Return selected paths relative to the result root, excluding metadata itself."""
        omitted = {"metadata", *exclude}
        return {
            name: str(path.relative_to(self.output_dir))
            for name, path in self.paths.items() if name not in omitted
        }


def _declare_artifacts(output_dir, entries, *, label, extra_owned=()) -> ArtifactFamily:
    """Resolve each named (relative path, produced) entry exactly once."""
    output_dir = Path(output_dir)
    owned = {name: output_dir / relative for name, (relative, _) in entries.items()}
    return ArtifactFamily(
        output_dir,
        {name: owned[name] for name, (_, produced) in entries.items() if produced},
        (*owned.values(), *extra_owned), label,
    )


class AnnotationDirectoryWriter:
    """Write persistent annotation shards and stream complete query diagnostics."""

    @staticmethod
    def artifact_family(output_dir, *, chromosomes=(), gene_lists=False, diagnostics=()) -> ArtifactFamily:
        output_dir = Path(output_dir)
        entries = {
            'metadata': ('diagnostics/metadata.json', True),
            'log': ('diagnostics/annotate.log', True),
            'dropped_snps': ('diagnostics/dropped_snps/dropped.tsv.gz', True),
            'query_annotation_status': ('diagnostics/query_annotation_status.tsv', True),
            'gene_list_audit': ('diagnostics/gene_list_audit.tsv.gz', gene_lists),
            'gene_list_resolution_summary': ('diagnostics/gene_list_resolution_summary.tsv', gene_lists),
            'gene_catalog_issues': ('diagnostics/gene_catalog_issues.tsv', 'gene_catalog_issues' in diagnostics),
            'input_issues': ('diagnostics/input_issues.tsv', 'input_issues' in diagnostics),
            'chromosome_scope': ('diagnostics/chromosome_scope.json', 'chromosome_scope' in diagnostics),
        }
        entries.update({f'query_{chrom}': (f'query.{chrom}.annot.gz', True) for chrom in chromosomes})
        return _declare_artifacts(output_dir, entries, label='annotation output artifact', extra_owned=tuple(output_dir.glob('query.*.annot.gz')))

    def write_diagnostics(self, output_dir, *, batch=None, statuses=None, drops=None, catalog_issues=None, input_issues=None, scope=None):
        """Persist every available gate result without inventing later-stage data."""
        from ._kernel.snp_identity import IDENTITY_DROP_COLUMNS

        root = Path(output_dir) / 'diagnostics'
        root.mkdir(parents=True, exist_ok=True)
        paths = {}
        if batch is not None:
            paths['gene_list_audit'] = root / 'gene_list_audit.tsv.gz'
            paths['gene_list_resolution_summary'] = root / 'gene_list_resolution_summary.tsv'
            batch.write_audit(paths['gene_list_audit'])
            batch.summary.to_csv(paths['gene_list_resolution_summary'], sep='\t', index=False, na_rep='')
        if statuses is not None:
            paths['query_annotation_status'] = root / 'query_annotation_status.tsv'
            columns = ['query', 'source', 'input_type', 'status', 'reason', 'n_annotation_snps', 'details']
            pd.DataFrame([item.as_dict() for item in statuses], columns=columns).to_csv(paths['query_annotation_status'], sep='\t', index=False, na_rep='')
        if drops is not None:
            paths['dropped_snps'] = root / 'dropped_snps' / 'dropped.tsv.gz'
            drops.write_tsv(paths['dropped_snps'], columns=IDENTITY_DROP_COLUMNS)
        for key, frame in [('gene_catalog_issues', catalog_issues), ('input_issues', input_issues)]:
            if frame is not None:
                paths[key] = root / f'{key}.tsv'
                frame.to_csv(paths[key], sep='\t', index=False, na_rep='')
        if scope is not None:
            paths['chromosome_scope'] = root / 'chromosome_scope.json'
            _atomic_write_json(scope, paths['chromosome_scope'])
        return paths

    def write(self, bundle, output_dir, *, overwrite=False, provenance=None, scope=None, catalog_issues=None, preflighted=False):
        """Emit chromosome shards incrementally and return persistent path references."""
        output_dir = Path(output_dir)
        diagnostics = ['chromosome_scope'] if scope is not None else []
        if catalog_issues is not None:
            diagnostics.append('gene_catalog_issues')
        family = self.artifact_family(output_dir, chromosomes=bundle.chromosomes, gene_lists=bundle.gene_list_batch is not None, diagnostics=diagnostics)
        # The workflow already preflighted before writing its log/diagnostics.
        # Reconcile the final produced family for successful stale cleanup.
        stale = family.preflight(overwrite=overwrite or preflighted)
        paths = self.write_diagnostics(output_dir, batch=bundle.gene_list_batch, statuses=bundle.query_statuses,
                                       drops=bundle.identity_drops, scope=scope, catalog_issues=catalog_issues)
        written = []
        row_batch = max(1, min(65536, 16 * 1024 * 1024 // (8 * max(1, len(bundle.query_columns)))))
        for chrom in bundle.chromosomes:
            metadata = bundle.metadata_for_chromosome(chrom).rename(columns={'POS': 'BP'})
            preferred = [name for name in ('CHR', 'BP', 'SNP', 'CM') if name in metadata]
            metadata = metadata.loc[:, [*preferred, *[name for name in metadata if name not in preferred]]]
            path = output_dir / f'query.{chrom}.annot.gz'
            with gzip.open(path, 'wt') as stream:
                for start in range(0, len(metadata), row_batch):
                    stop = min(start + row_batch, len(metadata))
                    values = bundle.read(chrom, rows=slice(start, stop), columns=bundle.query_columns)
                    query = pd.DataFrame(values.astype(np.uint8), columns=bundle.query_columns)
                    frame = pd.concat([metadata.iloc[start:stop].reset_index(drop=True), query], axis=1)
                    frame.to_csv(stream, sep='\t', index=False, header=start == 0, na_rep='NA')
            written.append(path)
            del metadata
        config = bundle.config_snapshot
        payload = {**(provenance or {}), 'artifact_type': 'annotation_projection',
                   'snp_identifier': config.snp_identifier, 'genome_build': config.genome_build,
                   'chromosomes': bundle.chromosomes, 'n_snps': bundle.n_rows,
                   'baseline_columns': bundle.baseline_columns, 'query_columns': bundle.query_columns,
                   'files': {'query_annotations': [path.name for path in written],
                             **{key: str(path.relative_to(output_dir)) for key, path in paths.items()}}}
        _atomic_write_json(payload, family.paths['metadata'])
        remove_output_artifacts(stale)
        bundle.output_paths = {**{key: str(path) for key, path in paths.items()}, 'metadata': str(family.paths['metadata']),
                               'query_annotations': [str(path) for path in written]}
        return bundle.output_paths


class H2DirectoryWriter:
    """Write the unpartitioned h2 summary, exact fitted-bin diagnostic, and metadata."""

    @staticmethod
    def artifact_family(output_dir) -> ArtifactFamily:
        """Declare h2 outputs and derived roots invalidated by a successful rewrite."""
        return _declare_artifacts(output_dir, {
            "summary": ("h2.tsv", True),
            "ld_score_regression_bins": ("diagnostics/ld_score_regression_bins.tsv", True),
            "metadata": ("diagnostics/metadata.json", True),
            "plots": ("plots", False),
            "postprocessing": ("postprocessing", False),
        }, label="h2 output artifact")

    def write(
        self,
        summary: pd.DataFrame,
        output_config: H2OutputConfig,
        *,
        metadata: dict[str, object],
        diagnostic_bins: pd.DataFrame,
    ) -> dict[str, str]:
        """Write the h2 summary, regression bins, and metadata.

        Parameters
        ----------
        summary : pandas.DataFrame
            One-row h2 summary in the canonical public schema.
        output_config : H2OutputConfig
            Output directory and overwrite policy.
        metadata : dict
            Scientific and source provenance for the result.
        diagnostic_bins : pandas.DataFrame
            Exact post-filter LD Score regression-bin summary in
            :data:`H2_REGRESSION_BIN_COLUMNS` order.

        Returns
        -------
        dict of str to str
            Written summary, bin-diagnostic, and metadata paths.

        Raises
        ------
        ValueError
            If the summary or diagnostic table lacks required columns.
        FileExistsError
            If an owned h2 artifact exists and overwrite is disabled.

        Notes
        -----
        Existing fixed h2 artifacts are checked before any output file is
        written. Replacement requires ``output_config.overwrite=True``. A
        successful replacement also removes default plot and post-processing
        roots derived from the superseded h2 result.
        """
        output_dir = ensure_output_directory(output_config.output_dir, label="output directory")
        family = self.artifact_family(output_dir)
        paths = family.paths
        summary_path = paths["summary"]
        metadata_path = paths["metadata"]
        bins_path = paths["ld_score_regression_bins"]
        diagnostics_dir = metadata_path.parent
        stale = family.preflight(overwrite=output_config.overwrite)
        diagnostics_dir.mkdir(parents=True, exist_ok=True)
        _atomic_write_dataframe(summary, summary_path, na_rep="NaN")
        _atomic_write_dataframe(
            _select_columns(diagnostic_bins, H2_REGRESSION_BIN_COLUMNS, label="h2 regression-bin diagnostic"),
            bins_path,
            na_rep="NaN",
        )
        _atomic_write_json(
            _result_metadata(
                metadata,
                artifact_type="h2_result",
                files=family.metadata_files(),
            ),
            metadata_path,
        )
        remove_output_artifacts(stale)
        return {name: str(path) for name, path in family.paths.items()}


@dataclass(frozen=True)
class QueryR2OutputConfig:
    """Directory-oriented output config for query-r2 pair results.

    Parameters
    ----------
    output_dir : str or os.PathLike[str]
        Directory that receives ``query_r2.tsv`` and ``diagnostics/metadata.json``.
    overwrite : bool, optional
        If ``True``, replace existing query-r2 outputs. Default is ``False``.
    """

    output_dir: str | PathLike[str]
    overwrite: bool = False

    def __post_init__(self) -> None:
        """Normalize the output directory."""
        object.__setattr__(self, "output_dir", _normalize_required_path(self.output_dir))


class QueryR2DirectoryWriter:
    """Write query-r2 results and diagnostic provenance.

    The workflow supplies input columns plus ``r2``, nullable ``sign_r``,
    ``r``, and ``status``. The TSV preserves those columns and writes missing
    values as ``NaN``; it does not add compatibility aliases.
    """

    @staticmethod
    def artifact_family(output_dir) -> ArtifactFamily:
        """Declare the pair-query result table and provenance."""
        return _declare_artifacts(output_dir, {
            "result": ("query_r2.tsv", True),
            "metadata": ("diagnostics/metadata.json", True),
        }, label="query-r2 output artifact")

    def write(
        self,
        result: pd.DataFrame,
        output_config: QueryR2OutputConfig,
        *,
        metadata: dict[str, object],
    ) -> dict[str, str]:
        """Write ``query_r2.tsv`` and ``diagnostics/metadata.json`` to a result directory.

        Existing fixed query-r2 artifacts are checked before any output file is
        written. Replacement requires ``output_config.overwrite=True``.
        """
        output_dir = ensure_output_directory(output_config.output_dir, label="output directory")
        family = self.artifact_family(output_dir)
        result_path = family.paths["result"]
        metadata_path = family.paths["metadata"]
        diagnostics_dir = metadata_path.parent
        family.preflight(overwrite=output_config.overwrite)
        diagnostics_dir.mkdir(parents=True, exist_ok=True)
        _atomic_write_dataframe(result, result_path, na_rep="NaN")
        _atomic_write_json(
            _result_metadata(metadata, artifact_type="query_r2_result", files=family.metadata_files()),
            metadata_path,
        )
        return {name: str(path) for name, path in family.paths.items()}


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
        object.__setattr__(self, "output_dir", _normalize_required_path(self.output_dir))
        allowed = {"snappy", "gzip", "brotli", "zstd", "none", None}
        if self.parquet_compression not in allowed:
            raise LDSCConfigError(
                "LD-score output configuration could not use parquet_compression="
                f"{self.parquet_compression!r}. Most likely the output config was created with an unsupported "
                "compression codec. Choose one of 'snappy', 'gzip', 'brotli', 'zstd', 'none', or None."
            )


class LDScoreDirectoryWriter:
    """Write canonical LD-score result directories.

    The writer owns ``metadata.json``, one shared baseline, optional single or
    numbered query batch files, and conditional ``ldscore.overlap.parquet``.
    Annotation counts and the ordered ``query_batches`` manifest are embedded
    in root metadata. Each Parquet file has chromosome-aligned row groups.
    ``write_batches`` releases completed query values, publishes after complete
    success, and returns an ``LDScoreSource``. ``write`` serializes an already
    materialized single result and returns artifact paths.
    """

    @staticmethod
    def artifact_family(output_dir, result=None, *, diagnostics_only=False, gene_list_batch=None) -> ArtifactFamily:
        """Describe all owned LD-score artifacts and the outputs selected by a result.

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
        }
        if scientific:
            if len(query_batches) > 1:
                for ordinal, entry in enumerate(query_batches, 1):
                    entries[f"query_batch{ordinal:05d}"] = (entry["file"], True)
            for chrom in (getattr(result, "identity_drops_by_chrom", {}) or {}):
                entries[f"dropped_snps_chr{chrom}"] = (f"diagnostics/dropped_snps/chr{chrom}_dropped.tsv.gz", True)
        return _declare_artifacts(
            output_dir, entries, label="LD-score output artifact",
            extra_owned=[*sorted((output_dir / "diagnostics/dropped_snps").glob("chr*_dropped.tsv.gz")),
                         *sorted(path for path in output_dir.glob("ldscore.query.batch*.parquet")
                                 if path.name.removeprefix("ldscore.query.batch").removesuffix(".parquet").isdigit())],
        )

    def write_batches(self, results, output_config: LDScoreOutputConfig, *, query_status_order=()):
        """Publish sequential materialized batches and return a saved source.

        Parameters
        ----------
        results : iterable of LDScoreResult
            Nonempty sequence of genome-wide results with identical baseline
            values and SNP rows, disjoint query columns, and corresponding
            counts and overlap statistics. Use a generator that releases each
            yielded result when resumed to bound retained query values.
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
        Each query table is privately written under the output directory and
        released before requesting the next result. One query batch uses
        ``ldscore.query.parquet``; multiple batches use numbered files. After
        all batches succeed, the writer publishes shared artifacts and the
        ordered ``query_batches`` manifest, with root metadata last. Handled
        failure closes the producer and removes owned staging. Completed
        batches are not resumable checkpoints.
        """
        from ._ldscore_batch_output import write_ldscore_batches

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
        self._validate_tables(result)

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
            from .overlap_matrix import overlap_to_long_frame
            overlap_to_long_frame(overlap).to_parquet(paths["overlap"], index=False)
        for chrom, artifact in identity_drops_by_chrom.items():
            artifact.write_to(paths[f"dropped_snps_chr{chrom}"])
        self._write_query_diagnostic_files(result, paths)
        metadata = self.build_metadata(
            result,
            files=family.metadata_files(exclude=("query_status", "gene_list_audit", "gene_list_resolution_summary", "chromosome_scope", "input_issues")),
            baseline_rg=baseline_rg,
            query_rg=query_rg,
        )
        paths["metadata"].write_text(json.dumps(_to_serializable(metadata), indent=2, sort_keys=True), encoding="utf-8")
        remove_output_artifacts(stale_paths)
        return {name: str(path) for name, path in paths.items()}

    @staticmethod
    def _write_query_diagnostic_files(result: Any, paths: dict[str, Path]) -> None:
        """Serialize fixed query diagnostics after family preflight."""
        query_statuses = tuple(getattr(result, "query_statuses", ()))
        gene_list_batch = getattr(result, "gene_list_batch", None)
        if "chromosome_scope" in paths:
            paths["chromosome_scope"].parent.mkdir(parents=True, exist_ok=True)
            _atomic_write_json(getattr(result, "chromosome_scope", None) or result.source_summary["chromosome_scope"], paths["chromosome_scope"])
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
        """Serialize the approved row audit and per-source summary schemas."""
        paths["gene_list_audit"].parent.mkdir(parents=True, exist_ok=True)
        batch.write_audit(paths["gene_list_audit"])
        batch.summary.to_csv(
            paths["gene_list_resolution_summary"],
            sep="\t",
            index=False,
            na_rep="",
        )

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
        unavailable and per-column common counts are omitted. Row-group metadata
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

    def _validate_tables(self, result: Any) -> None:
        """Validate baseline/query table shape before any files are written."""
        baseline_table = getattr(result, "baseline_table")
        query_table = getattr(result, "query_table", None)
        baseline_columns = list(getattr(result, "baseline_columns", []))
        query_columns = list(getattr(result, "query_columns", []))
        config_snapshot = getattr(result, "config_snapshot", None)
        if config_snapshot is None:
            raise LDSCInternalError(
                "LD-score output writer could not validate tables because LDScoreResult.config_snapshot is missing. "
                "Most likely an upstream LD-score workflow returned an older or incomplete result object. "
                "Re-run with DEBUG logging and report the traceback."
            )
        snp_identifier = config_snapshot.snp_identifier
        identity_columns = ["A1", "A2"] if is_allele_aware_mode(snp_identifier) else []
        required_baseline = ["CHR", "POS", "SNP", *identity_columns, REGRESSION_LD_SCORE_COLUMN, *baseline_columns]
        missing = [column for column in required_baseline if column not in baseline_table.columns]
        if missing:
            if any(column in missing for column in ("A1", "A2")):
                raise LDSCInternalError(
                    "LD-score output writer could not validate baseline_table for allele-aware output because "
                    f"required columns are missing: {missing}. Most likely an upstream LD-score step dropped A1/A2 "
                    "from an allele-aware result table. Regenerate the LD-score artifact with the current LDSC "
                    "package; if this was produced by this run, re-run with DEBUG logging and report the traceback."
                )
            raise LDSCInternalError(
                "LD-score output writer could not validate baseline_table because required columns are missing: "
                f"{missing}. Most likely an upstream LD-score step returned a result table with the wrong schema. "
                "Re-run with DEBUG logging and report the traceback."
            )
        _validate_ldscore_allele_columns(baseline_table, table_name="baseline_table", snp_identifier=snp_identifier)
        if query_columns and query_table is None:
            raise LDSCInternalError(
                "LD-score output writer could not write query annotations because result.query_columns is non-empty "
                "but query_table is missing. Most likely an upstream LD-score step returned an incomplete query "
                "annotation result. Re-run with DEBUG logging and report the traceback."
            )
        if not query_columns and query_table is not None:
            raise LDSCInternalError(
                "LD-score output writer could not write the result because query_table was provided while "
                "result.query_columns is empty. Most likely an upstream LD-score step returned inconsistent "
                "query annotation metadata. Re-run with DEBUG logging and report the traceback."
            )
        if query_table is None:
            return
        required_query = ["CHR", "POS", "SNP", *identity_columns, *query_columns]
        missing = [column for column in required_query if column not in query_table.columns]
        if missing:
            if any(column in missing for column in ("A1", "A2")):
                raise LDSCInternalError(
                    "LD-score output writer could not validate query_table for allele-aware output because required "
                    f"columns are missing: {missing}. Most likely an upstream LD-score step dropped A1/A2 from an "
                    "allele-aware query table. Regenerate the LD-score artifact with the current LDSC package; if "
                    "this was produced by this run, re-run with DEBUG logging and report the traceback."
                )
            raise LDSCInternalError(
                "LD-score output writer could not validate query_table because required columns are missing: "
                f"{missing}. Most likely an upstream LD-score step returned a query table with the wrong schema. "
                "Re-run with DEBUG logging and report the traceback."
            )
        _validate_ldscore_allele_columns(query_table, table_name="query_table", snp_identifier=snp_identifier)
        assert_same_snp_rows(
            baseline_table,
            query_table,
            context="query rows must match baseline rows on CHR/SNP/POS",
            snp_identifier=snp_identifier,
        )


def _validate_ldscore_allele_columns(table: pd.DataFrame, *, table_name: str, snp_identifier: str) -> None:
    """Reject allele-aware LD-score artifacts without usable allele identity."""
    if not is_allele_aware_mode(snp_identifier):
        return
    missing = [column for column in ("A1", "A2") if column not in table.columns]
    if missing:
        raise LDSCInputError(
            f"LD-score artifact table {table_name} is malformed for snp_identifier='{snp_identifier}' because "
            f"allele columns are missing: {missing}. Most likely this artifact was produced by an older package "
            "version or copied without allele-aware columns. Regenerate the LD-score artifact with the current "
            "LDSC package."
        )
    try:
        effective_merge_key_series(table, snp_identifier, context=table_name)
    except (LDSCInputError, ValueError) as exc:
        raise LDSCInputError(
            f"LD-score artifact table {table_name} is malformed for snp_identifier='{snp_identifier}' because "
            f"its A1/A2 columns are unusable: {exc}. Most likely the artifact contains malformed allele values "
            "or duplicate allele-aware SNP identities. Regenerate the LD-score artifact with the current LDSC package."
        ) from exc


@dataclass(frozen=True)
class QuantileH2OutputConfig:
    """Output directory and replacement policy for ``quantile-h2`` results.

    Parameters
    ----------
    output_dir : str or os.PathLike[str]
        Directory receiving the two scientific TSVs and ``diagnostics/``.
    overwrite : bool, optional
        Replace existing owned quantile-h2 artifacts when ``True``. Default is
        ``False``.
    """

    output_dir: str | PathLike[str]
    overwrite: bool = False

    def __post_init__(self) -> None:
        """Normalize the required output directory."""
        object.__setattr__(self, "output_dir", _normalize_required_path(self.output_dir))


class QuantileH2DirectoryWriter:
    """Write post-fit quantile summaries, coefficients, and SNP diagnostics.

    A successful overwrite removes the default plot root derived from the
    superseded quantile result.
    """

    @staticmethod
    def artifact_family(output_dir) -> ArtifactFamily:
        """Declare quantile outputs, alignment diagnostics, and the derived plot root."""
        return _declare_artifacts(output_dir, {
            "quantile_h2": ("quantile_h2.tsv", True),
            "standardized_coefficients": ("standardized_coefficients.tsv", True),
            "metadata": ("diagnostics/metadata.json", True),
            "snp_alignment_issues": ("diagnostics/snp_alignment_issues.tsv.gz", True),
            "plots": ("plots", False),
        }, label="quantile-h2 output artifact")

    def write(
        self,
        quantile_h2: pd.DataFrame,
        standardized_coefficients: pd.DataFrame,
        alignment_issues_path: str | Path,
        output_config: QuantileH2OutputConfig,
        metadata: dict[str, object],
    ) -> dict[str, str]:
        """Write the complete successful ``quantile-h2`` artifact family.

        Parameters
        ----------
        quantile_h2 : pandas.DataFrame
            Low-to-high quantile rows in the stable public schema.
        standardized_coefficients : pandas.DataFrame
            One row per fitted annotation in coefficient order.
        alignment_issues_path : str or pathlib.Path
            Complete streamed gzip diagnostic table, including its header.
        output_config : QuantileH2OutputConfig
            Output directory and replacement policy.
        metadata : dict
            Scientific and source provenance added to the standard result
            metadata envelope.

        Returns
        -------
        dict of str to str
            Paths for both scientific tables, metadata, and SNP diagnostics.
        """
        output_dir = ensure_output_directory(output_config.output_dir, label="output directory")
        family = self.artifact_family(output_dir)
        paths = family.paths
        diagnostics_dir = paths["metadata"].parent
        stale = family.preflight(overwrite=output_config.overwrite)
        diagnostics_dir.mkdir(parents=True, exist_ok=True)
        _atomic_write_dataframe(
            _select_columns(quantile_h2, QUANTILE_H2_COLUMNS, label="quantile-h2 summary"),
            paths["quantile_h2"],
            na_rep="NaN",
        )
        _atomic_write_dataframe(
            _select_columns(
                standardized_coefficients,
                STANDARDIZED_COEFFICIENT_COLUMNS,
                label="standardized coefficients",
            ),
            paths["standardized_coefficients"],
            na_rep="NaN",
        )
        if Path(alignment_issues_path).resolve() != paths["snp_alignment_issues"].resolve():
            shutil.copyfile(alignment_issues_path, paths["snp_alignment_issues"])
        payload = _result_metadata(
            metadata,
            artifact_type="quantile_h2_result",
            files=family.metadata_files(),
        )
        _atomic_write_json(payload, paths["metadata"])
        remove_output_artifacts(stale)
        return {name: str(path) for name, path in paths.items()}


@dataclass(frozen=True)
class PartitionedH2OutputConfig:
    """Directory-oriented output config for partitioned-h2 regression summaries.

    The config keeps the partitioned-h2 output contract explicit and separate
    from LD-score parquet output. The default mode writes only the stable
    aggregate summary. ``write_per_query_results=True`` enables the richer
    per-query tree without changing the aggregate file name.

    Parameters
    ----------
    output_dir : str or os.PathLike[str]
        Directory that receives ``partitioned_h2.tsv`` and, when requested,
        the optional ``diagnostics/query_annotations`` per-query result tree.
    overwrite : bool, optional
        If ``True``, replace existing fixed partitioned-h2 outputs and remove
        stale owned siblings after a successful write. Default is ``False``.
    write_per_query_results : bool, optional
        If ``True``, also write one subdirectory per query annotation under
        ``diagnostics/query_annotations``. Default is ``False``.
    """

    output_dir: str | PathLike[str]
    overwrite: bool = False
    write_per_query_results: bool = False

    def __post_init__(self) -> None:
        """Normalize the output directory."""
        object.__setattr__(self, "output_dir", _normalize_required_path(self.output_dir))


@dataclass(frozen=True)
class PartitionedH2FitArtifacts:
    """Paths to one complete fitted model, owned by its containing directory."""

    full: Path
    delete_values: Path
    metadata: Path


def stage_partitioned_h2_fit(directory, category_table, delete_values, metadata):
    """Write one fit privately before releasing its estimator and matrices.

    The caller supplies a unique directory inside its owned output workspace.
    Returned paths remain private until the final summary determines folder
    ordering and the directory writer publishes the complete query tree.
    Concurrent query workers receive distinct ordinal directories from the
    coordinator. Only the coordinator assigns public folder names and publishes
    results, after applying the run's query-error policy.
    """
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=False)
    artifacts = PartitionedH2FitArtifacts(directory / "partitioned_h2_full.tsv",
                                         directory / "coefficient_delete_values.parquet",
                                         directory / "metadata.json")
    _atomic_write_dataframe(_select_columns(category_table, PARTITIONED_H2_COLUMNS, label="full partitioned-h2 summary"),
                            artifacts.full, na_rep="NaN")
    delete_values.to_parquet(artifacts.delete_values, index=False)
    _atomic_write_json({**metadata, "coefficient_delete_block_count": int(len(delete_values)),
                       "retained_ld_columns": [c for c in delete_values if c != "delete_block"]}, artifacts.metadata)
    return artifacts


class PartitionedH2DirectoryWriter:
    """Write aggregate and optional per-query partitioned-h2 outputs.

    The writer owns the fixed regression summary layout:
    ``partitioned_h2.tsv`` at the output root uses the compact public schema,
    while optional diagnostic query folders contain the same one-row summary
    plus a full baseline-plus-query ``partitioned_h2_full.tsv`` table. The
    per-query tree
    is written to a temporary sibling directory before it is moved into the
    final location, so ordinary validation and I/O failures do not expose a
    partially populated final tree. A successful overwrite also removes the
    default plot root derived from the superseded result.
    An optional ``diagnostics/query_status.tsv`` records all attempted queries,
    including failed queries absent from the summary and model tree. The
    regression runner owns the strict/continue decision; this writer publishes
    only the successful model artifacts supplied by that runner.
    """

    @staticmethod
    def artifact_family(output_dir, *, write_per_query_results=False, coefficient_delete_values=False,
                        query_status=False) -> ArtifactFamily:
        """Declare summaries, model artifacts, and an optional attempted-query ledger."""
        return _declare_artifacts(output_dir, {
            "summary": ("partitioned_h2.tsv", True),
            "metadata": ("diagnostics/metadata.json", True),
            "coefficient_delete_values": ("diagnostics/coefficient_delete_values.parquet", coefficient_delete_values),
            "query_annotations": ("diagnostics/query_annotations", write_per_query_results),
            "query_status": ("diagnostics/query_status.tsv", query_status),
            "plots": ("plots", False),
        }, label="partitioned-h2 output artifact")

    @staticmethod
    def write_query_status(output_dir, query_status):
        """Write attempted-query statuses after the caller's output preflight.

        This diagnostic can survive a failed run without publishing any
        scientific results. Successful runs write it with their result family.

        Parameters
        ----------
        output_dir : path-like
            Caller-preflighted result directory. This method does not request
            overwrite authorization or write scientific result files.
        query_status : pandas.DataFrame
            Attempted queries in input order, with ``query_annotation``,
            ``status``, ``stage``, ``error_type``, and ``error_message`` columns.

        Returns
        -------
        pathlib.Path
            Path to the written ``diagnostics/query_status.tsv``.
        """
        path = PartitionedH2DirectoryWriter.artifact_family(output_dir, query_status=True).paths["query_status"]
        path.parent.mkdir(parents=True, exist_ok=True)
        _atomic_write_dataframe(query_status, path, na_rep="NaN")
        return path

    def write(
        self,
        summary: pd.DataFrame,
        output_config: PartitionedH2OutputConfig,
        per_query_artifacts: dict[str, PartitionedH2FitArtifacts] | None = None,
        metadata: dict[str, object] | None = None,
        coefficient_delete_values: pd.DataFrame | None = None,
        query_status: pd.DataFrame | None = None,
    ) -> dict[str, str]:
        """Write partitioned-h2 summary artifacts.

        Parameters
        ----------
        summary : pandas.DataFrame
            Aggregate partitioned-h2 table with the compact public columns in
            ``PARTITIONED_H2_COLUMNS``.
        output_config : PartitionedH2OutputConfig
            Output directory, overwrite policy, and per-query mode.
        per_query_artifacts : dict of str to PartitionedH2FitArtifacts, optional
            Completed private or persistent fits keyed by original query name.
            Files are copied into the sorted publication tree without loading
            category tables or delete-value matrices. Workers may stage fits
            concurrently; the coordinator calls this writer after fitting and
            decides which successful artifacts may be published.
        metadata : dict, optional
            Run-level metadata copied into every per-query ``metadata.json``.
        coefficient_delete_values : pandas.DataFrame, optional
            Baseline-only fitted model's delete-one-block coefficient vectors.
            The first column is ``delete_block`` and remaining columns retain
            fitted annotation order.
        query_status : pandas.DataFrame, optional
            One diagnostic row per attempted query, including failures absent
            from the scientific summary and per-query artifact tree.

        Returns
        -------
        dict of str to str
            Written path map. Always includes ``"summary"`` and additionally
            includes ``"per_query_root"`` and ``"per_query_manifest"`` when
            per-query output is enabled, and ``"query_status"`` when the
            attempted-query table is supplied.

        Raises
        ------
        ValueError
            If ``summary`` lacks the required ``Category`` column.
        FileExistsError
            If a final output path already exists and overwrite is disabled.
        """
        self._validate_summary(summary)
        output_dir = ensure_output_directory(output_config.output_dir, label="output directory")
        family = self.artifact_family(
            output_dir, write_per_query_results=output_config.write_per_query_results,
            coefficient_delete_values=coefficient_delete_values is not None,
            query_status=query_status is not None,
        )
        summary_path = family.paths["summary"]
        metadata_path = family.paths["metadata"]
        diagnostics_dir = metadata_path.parent
        query_root = family.paths.get("query_annotations")
        coefficient_delete_path = family.paths.get("coefficient_delete_values")
        stale_paths = family.preflight(overwrite=output_config.overwrite)
        root_files = family.metadata_files()
        paths = {name: str(path) for name, path in family.paths.items() if name != "query_annotations"}
        metadata_payload = dict(metadata or {})
        if coefficient_delete_values is not None:
            metadata_payload["coefficient_delete_block_count"] = int(len(coefficient_delete_values))
            metadata_payload["retained_ld_columns"] = [
                column for column in coefficient_delete_values.columns if column != "delete_block"
            ]
        diagnostics_dir.mkdir(parents=True, exist_ok=True)
        if query_status is not None:
            self.write_query_status(output_dir, query_status)
        if not output_config.write_per_query_results:
            _atomic_write_dataframe(
                _select_columns(summary, PARTITIONED_H2_COLUMNS, label="partitioned-h2 summary"),
                summary_path,
                na_rep="NaN",
            )
            if coefficient_delete_values is not None:
                coefficient_delete_values.to_parquet(coefficient_delete_path, index=False)
            _atomic_write_json(
                _result_metadata(metadata_payload, artifact_type="partitioned_h2_result", files=root_files),
                metadata_path,
            )
            remove_output_artifacts(stale_paths)
            return paths

        query_records = self._query_records(summary)
        staging_dir = Path(tempfile.mkdtemp(prefix=".query_annotations.tmp.", dir=str(diagnostics_dir)))
        backup_dir: Path | None = None
        try:
            manifest_rows = self._write_staged_query_tree(
                staging_dir,
                query_records,
                summary,
                per_query_artifacts or {},
                metadata_payload,
            )
            _atomic_write_dataframe(pd.DataFrame(manifest_rows), staging_dir / "manifest.tsv", na_rep="NaN")
            _atomic_write_dataframe(
                _select_columns(summary, PARTITIONED_H2_COLUMNS, label="partitioned-h2 summary"),
                summary_path,
                na_rep="NaN",
            )
            _atomic_write_json(
                _result_metadata(metadata_payload, artifact_type="partitioned_h2_result", files=root_files),
                metadata_path,
            )
            if output_config.overwrite and query_root.exists():
                backup_dir = Path(tempfile.mkdtemp(prefix=".query_annotations.backup.", dir=str(diagnostics_dir)))
                backup_dir.rmdir()
                os.replace(query_root, backup_dir)
            os.replace(staging_dir, query_root)
        except Exception:
            if staging_dir.exists():
                shutil.rmtree(staging_dir, ignore_errors=True)
            if backup_dir is not None and backup_dir.exists() and not query_root.exists():
                os.replace(backup_dir, query_root)
            raise
        else:
            if backup_dir is not None and backup_dir.exists():
                shutil.rmtree(backup_dir, ignore_errors=True)

        remove_output_artifacts(stale_paths)
        paths["per_query_root"] = str(query_root)
        paths["per_query_manifest"] = str(query_root / "manifest.tsv")
        return paths

    def _validate_summary(self, summary: pd.DataFrame) -> None:
        """Validate the aggregate summary before writing any final files."""
        if "category" not in summary.columns:
            raise LDSCInternalError(
                "partitioned-h2 output writer could not validate the aggregate summary because the required "
                "'category' column is missing. Most likely the regression workflow returned a result table with "
                "the wrong schema. Re-run with DEBUG logging and report the traceback."
            )

    def _query_records(self, summary: pd.DataFrame) -> list[dict[str, object]]:
        """Return ordered query records with deterministic safe folder names."""
        records: list[dict[str, object]] = []
        width = max(4, len(str(len(summary))))
        for idx, query_name in enumerate(summary["category"].astype(str).tolist(), start=1):
            slug = _slugify_query_name(query_name)
            folder = f"{idx:0{width}d}_{slug}"
            records.append(
                {
                    "ordinal": idx,
                    "query_annotation": query_name,
                    "slug": slug,
                    "folder": folder,
                }
            )
        return records

    def _write_staged_query_tree(
        self,
        staging_dir: Path,
        query_records: list[dict[str, object]],
        summary: pd.DataFrame,
        per_query_artifacts: dict[str, PartitionedH2FitArtifacts],
        metadata: dict[str, object],
    ) -> list[dict[str, object]]:
        """Populate the staged per-query result tree and return manifest rows."""
        manifest_rows: list[dict[str, object]] = []
        for record in query_records:
            query_name = str(record["query_annotation"])
            folder = str(record["folder"])
            query_dir = staging_dir / folder
            query_dir.mkdir(parents=True, exist_ok=False)
            query_summary = summary.loc[summary["category"].astype(str) == query_name].reset_index(drop=True)
            summary_rel = f"diagnostics/query_annotations/{folder}/partitioned_h2.tsv"
            full_rel = f"diagnostics/query_annotations/{folder}/partitioned_h2_full.tsv"
            metadata_rel = f"diagnostics/query_annotations/{folder}/metadata.json"
            delete_rel = f"diagnostics/query_annotations/{folder}/coefficient_delete_values.parquet"
            _atomic_write_dataframe(
                _select_columns(query_summary, PARTITIONED_H2_COLUMNS, label="query summary"),
                query_dir / "partitioned_h2.tsv",
                na_rep="NaN",
            )
            artifacts = per_query_artifacts.get(query_name)
            if artifacts is None:
                raise LDSCInternalError(
                    f"partitioned-h2 output writer is missing coefficient delete values for query {query_name!r}. "
                    "Most likely the fitted result family was assembled incompletely."
                )
            shutil.copyfile(artifacts.delete_values, query_dir / "coefficient_delete_values.parquet")
            shutil.copyfile(artifacts.full, query_dir / "partitioned_h2_full.tsv")
            fit_metadata = json.loads(artifacts.metadata.read_text(encoding="utf-8"))
            payload = {
                **metadata,
                **fit_metadata,
                "artifact_type": "partitioned_h2_query_result",
                "files": {"summary": summary_rel, "full": full_rel, "coefficient_delete_values": delete_rel},
                "ordinal": record["ordinal"],
                "query_annotation": query_name,
                "slug": record["slug"],
                "folder": folder,
            }
            _atomic_write_json(payload, query_dir / "metadata.json")
            manifest_rows.append(
                {
                    "ordinal": record["ordinal"],
                    "query_annotation": query_name,
                    "slug": record["slug"],
                    "folder": folder,
                    "summary_path": summary_rel,
                    "partitioned_h2_full_path": full_rel,
                    "metadata_path": metadata_rel,
                    "coefficient_delete_values_path": delete_rel,
                }
            )
        return manifest_rows


@dataclass(frozen=True)
class RgOutputConfig:
    """Directory-oriented output config for genetic-correlation summaries.

    Parameters
    ----------
    output_dir : str or os.PathLike[str]
        Directory that receives ``rg.tsv``, ``rg_full.tsv``, and
        ``h2_per_trait.tsv``.
    overwrite : bool, optional
        If ``True``, replace existing fixed rg outputs and remove stale owned
        siblings after a successful write. Default is ``False``.
    write_per_pair_detail : bool, optional
        If ``True``, also write one subdirectory per tested pair under
        ``diagnostics/pairs/``. Default is ``False``.
    """

    output_dir: str | PathLike[str]
    overwrite: bool = False
    write_per_pair_detail: bool = False

    def __post_init__(self) -> None:
        """Normalize the output directory."""
        object.__setattr__(self, "output_dir", _normalize_required_path(self.output_dir))


class RgDirectoryWriter:
    """Write concise, full, per-trait, and optional per-pair rg outputs.

    The workflow-level logger owns ``diagnostics/rg.log``. This writer owns the public data
    artifacts plus diagnostic metadata in the rg output family:
    ``rg.tsv``, ``rg_full.tsv``, ``h2_per_trait.tsv``, and optional
    ``diagnostics/pairs/``.
    The optional pair tree is staged before replacement so a failed write does
    not expose a partially populated final tree. A successful overwrite also
    removes the default plot root derived from the superseded result.
    """

    @staticmethod
    def artifact_family(output_dir, *, include_pairs=False) -> ArtifactFamily:
        """Declare rg summaries, optional pair details, and the derived plot root."""
        return _declare_artifacts(output_dir, {
            "metadata": ("diagnostics/metadata.json", True),
            "rg": ("rg.tsv", True),
            "rg_full": ("rg_full.tsv", True),
            "h2_per_trait": ("h2_per_trait.tsv", True),
            "pairs": ("diagnostics/pairs", include_pairs),
            "plots": ("plots", False),
        }, label="rg output artifact")

    def write(self, result: Any, output_config: RgOutputConfig) -> dict[str, str]:
        """Write an ``RgResultFamily``-like object to a result directory.

        Parameters
        ----------
        result : object
            Object exposing ``rg``, ``rg_full``, ``h2_per_trait``, and
            ``per_pair_metadata`` attributes. The first two tables must contain
            the canonical rg concise and full schemas.
        output_config : RgOutputConfig
            Output directory, overwrite policy, and optional per-pair detail
            setting.

        Returns
        -------
        dict of str to str
            Written path map. Always includes ``rg``, ``rg_full``, and
            ``h2_per_trait``; also includes ``pairs_root`` and
            ``pairs_manifest`` when per-pair detail is enabled.

        Raises
        ------
        ValueError
            If required result columns are missing or the pair metadata length
            does not match ``rg_full`` when detail output is requested.
        FileExistsError
            If an owned rg artifact already exists and overwrite is disabled.
        """
        rg = _select_columns(getattr(result, "rg"), RG_CONCISE_COLUMNS, label="rg summary")
        rg_full = _select_columns(getattr(result, "rg_full"), RG_FULL_COLUMNS, label="full rg summary")
        h2_per_trait = getattr(result, "h2_per_trait")
        per_pair_metadata = list(getattr(result, "per_pair_metadata", []))
        if output_config.write_per_pair_detail and len(per_pair_metadata) != len(rg_full):
            raise LDSCInternalError(
                "rg output writer could not write per-pair diagnostics because per_pair_metadata does not contain "
                "one record per rg_full row. Most likely the rg workflow returned inconsistent pair-detail "
                "metadata. Re-run with DEBUG logging and report the traceback."
            )

        output_dir = ensure_output_directory(output_config.output_dir, label="output directory")
        family = self.artifact_family(output_dir, include_pairs=output_config.write_per_pair_detail)
        metadata_path = family.paths["metadata"]
        diagnostics_dir = metadata_path.parent
        rg_path = family.paths["rg"]
        full_path = family.paths["rg_full"]
        h2_path = family.paths["h2_per_trait"]
        pairs_root = family.paths.get("pairs")
        stale_paths = family.preflight(overwrite=output_config.overwrite)
        paths = {name: str(path) for name, path in family.paths.items() if name != "pairs"}
        diagnostics_dir.mkdir(parents=True, exist_ok=True)
        if not output_config.write_per_pair_detail:
            _atomic_write_dataframe(rg, rg_path, na_rep="NaN")
            _atomic_write_dataframe(rg_full, full_path, na_rep="NaN")
            _atomic_write_dataframe(h2_per_trait, h2_path, na_rep="NaN")
            _atomic_write_json(_rg_root_metadata(result, files=family.metadata_files()), metadata_path)
            remove_output_artifacts(stale_paths)
            return paths

        pair_records = self._pair_records(rg_full, per_pair_metadata)
        staging_dir = Path(tempfile.mkdtemp(prefix=".pairs.tmp.", dir=str(diagnostics_dir)))
        backup_dir: Path | None = None
        try:
            manifest_rows = self._write_staged_pair_tree(staging_dir, pair_records, rg_full)
            _atomic_write_dataframe(pd.DataFrame(manifest_rows), staging_dir / "manifest.tsv", na_rep="NaN")
            _atomic_write_dataframe(rg, rg_path, na_rep="NaN")
            _atomic_write_dataframe(rg_full, full_path, na_rep="NaN")
            _atomic_write_dataframe(h2_per_trait, h2_path, na_rep="NaN")
            _atomic_write_json(_rg_root_metadata(result, files=family.metadata_files()), metadata_path)
            if output_config.overwrite and pairs_root.exists():
                backup_dir = Path(tempfile.mkdtemp(prefix=".pairs.backup.", dir=str(diagnostics_dir)))
                backup_dir.rmdir()
                os.replace(pairs_root, backup_dir)
            os.replace(staging_dir, pairs_root)
        except Exception:
            if staging_dir.exists():
                shutil.rmtree(staging_dir, ignore_errors=True)
            if backup_dir is not None and backup_dir.exists() and not pairs_root.exists():
                os.replace(backup_dir, pairs_root)
            raise
        else:
            if backup_dir is not None and backup_dir.exists():
                shutil.rmtree(backup_dir, ignore_errors=True)

        remove_output_artifacts(stale_paths)
        paths["pairs_root"] = str(pairs_root)
        paths["pairs_manifest"] = str(pairs_root / "manifest.tsv")
        return paths

    def _pair_records(self, rg_full: pd.DataFrame, per_pair_metadata: list[dict[str, object]]) -> list[dict[str, object]]:
        """Return ordered pair records with deterministic safe folder names."""
        width = max(4, len(str(len(rg_full))))
        records: list[dict[str, object]] = []
        for idx, (_, row) in enumerate(rg_full.iterrows(), start=1):
            metadata = dict(per_pair_metadata[idx - 1])
            trait_1 = str(row["trait_1"])
            trait_2 = str(row["trait_2"])
            slug = f"{_slugify_rg_trait_name(trait_1)}_vs_{_slugify_rg_trait_name(trait_2)}"
            folder = f"{idx:0{width}d}_{slug}"
            records.append(
                {
                    "ordinal": idx,
                    "trait_1": trait_1,
                    "trait_2": trait_2,
                    "slug": slug,
                    "folder": folder,
                    "metadata": metadata,
                }
            )
        return records

    def _write_staged_pair_tree(
        self,
        staging_dir: Path,
        pair_records: list[dict[str, object]],
        rg_full: pd.DataFrame,
    ) -> list[dict[str, object]]:
        """Populate the staged pair result tree and return manifest rows."""
        manifest_rows: list[dict[str, object]] = []
        for record in pair_records:
            folder = str(record["folder"])
            pair_dir = staging_dir / folder
            pair_dir.mkdir(parents=True, exist_ok=False)
            row = rg_full.iloc[[int(record["ordinal"]) - 1]].reset_index(drop=True)
            detail_rel = f"diagnostics/pairs/{folder}/rg_full.tsv"
            metadata_rel = f"diagnostics/pairs/{folder}/metadata.json"
            _atomic_write_dataframe(row, pair_dir / "rg_full.tsv", na_rep="NaN")
            payload = {
                **dict(record["metadata"]),
                "artifact_type": "rg_pair_result",
                "files": {"rg_full": detail_rel},
                "ordinal": record["ordinal"],
                "trait_1": record["trait_1"],
                "trait_2": record["trait_2"],
                "slug": record["slug"],
                "folder": folder,
            }
            _atomic_write_json(payload, pair_dir / "metadata.json")
            manifest_rows.append(
                {
                    "ordinal": record["ordinal"],
                    "trait_1": record["trait_1"],
                    "trait_2": record["trait_2"],
                    "slug": record["slug"],
                    "folder": folder,
                    "rg_full_path": detail_rel,
                    "metadata_path": metadata_rel,
                }
            )
        return manifest_rows


def _result_metadata(
    metadata: dict[str, object] | None,
    *,
    artifact_type: str,
    files: dict[str, str],
) -> dict[str, object]:
    """Return root result metadata with the canonical discriminator fields."""
    payload = dict(metadata or {})
    payload["artifact_type"] = artifact_type
    payload["files"] = dict(files)
    return payload


def _rg_root_metadata(result: Any, *, files: dict[str, str]) -> dict[str, object]:
    """Build diagnostic metadata for an rg result directory from result tables."""
    rg_full = getattr(result, "rg_full", pd.DataFrame())
    h2_per_trait = getattr(result, "h2_per_trait", pd.DataFrame())
    trait_names: list[str] = []
    if "trait_name" in h2_per_trait.columns:
        trait_names = [str(value) for value in h2_per_trait["trait_name"].dropna().tolist()]
    pair_kind = None
    if "pair_kind" in rg_full.columns:
        values = [str(value) for value in rg_full["pair_kind"].dropna().unique().tolist()]
        pair_kind = values[0] if len(values) == 1 else ("mixed" if values else None)
    payload = _result_metadata(
        {
            "trait_names": trait_names,
            "pair_kind": pair_kind,
            "n_pairs": int(len(rg_full)),
            "n_traits": int(len(trait_names)),
        },
        artifact_type="rg_result",
        files=files,
    )
    ldscore_dir = getattr(result, "ldscore_dir", None)
    if ldscore_dir is not None:
        payload["ldscore_dir"] = ldscore_dir
    sumstats_files = getattr(result, "sumstats_files", None)
    if sumstats_files is not None:
        payload["sumstats_files"] = list(sumstats_files)
    return payload


def _slugify_query_name(value: str) -> str:
    """Return a filesystem-safe slug for a query annotation name."""
    normalized = unicodedata.normalize("NFKD", str(value)).encode("ascii", "ignore").decode("ascii")
    slug = re.sub(r"[^a-z0-9._-]+", "_", normalized.lower()).strip("._-")
    return slug or "annotation"


def _slugify_rg_trait_name(value: str) -> str:
    """Return a filesystem-safe slug for an rg trait name."""
    normalized = unicodedata.normalize("NFKD", str(value)).encode("ascii", "ignore").decode("ascii")
    slug = re.sub(r"[^a-z0-9._-]+", "_", normalized.lower()).strip("._-")
    return slug or "trait"


def _select_columns(df: pd.DataFrame, columns: list[str], *, label: str) -> pd.DataFrame:
    """Return ``df`` with required public output columns in canonical order."""
    if df.empty and len(df.columns) == 0:
        return df
    missing = [column for column in columns if column not in df.columns]
    if missing:
        raise LDSCInternalError(
            f"output writer could not select public columns for {label} because required columns are missing: "
            f"{missing}. Most likely the upstream workflow returned a result table with the wrong schema. "
            "Re-run with DEBUG logging and report the traceback."
        )
    return df.loc[:, columns]


def _atomic_write_dataframe(df: pd.DataFrame, path: Path, *, na_rep: str | None = None) -> None:
    """Write a dataframe through a temporary sibling file, then replace."""
    fd, tmp_name = tempfile.mkstemp(prefix=f".{path.name}.", suffix=".tmp", dir=str(path.parent))
    os.close(fd)
    tmp_path = Path(tmp_name)
    try:
        kwargs = {"sep": "\t", "index": False}
        if na_rep is not None:
            kwargs["na_rep"] = na_rep
        df.to_csv(tmp_path, **kwargs)
        os.replace(tmp_path, path)
    except Exception:
        tmp_path.unlink(missing_ok=True)
        raise


def _atomic_write_json(payload: dict[str, object], path: Path) -> None:
    """Write JSON through a temporary sibling file, then replace."""
    fd, tmp_name = tempfile.mkstemp(prefix=f".{path.name}.", suffix=".tmp", dir=str(path.parent))
    os.close(fd)
    tmp_path = Path(tmp_name)
    try:
        tmp_path.write_text(json.dumps(_to_serializable(payload), indent=2, sort_keys=True), encoding="utf-8")
        os.replace(tmp_path, path)
    except Exception:
        tmp_path.unlink(missing_ok=True)
        raise


def _to_serializable(value: Any) -> Any:
    """Recursively convert workflow objects into JSON-serializable structures."""
    if is_dataclass(value):
        return _to_serializable(asdict(value))
    if isinstance(value, pd.DataFrame):
        return value.to_dict(orient="records")
    if isinstance(value, pd.Series):
        return value.to_list()
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, dict):
        return {str(key): _to_serializable(val) for key, val in value.items()}
    if isinstance(value, (list, tuple)):
        return [_to_serializable(item) for item in value]
    return value
