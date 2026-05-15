"""Public workflow module for SNP-level annotation building.

Overview
--------
This module owns the supported annotation workflow: path-token resolution,
annotation bundle loading, BED-to-SNP projection, CLI argument parsing, genome
build resolution, output preflight, and fixed ``query.<chrom>.annot.gz``
writing. Parsed workflow entry points also write ``annotate.log`` beside the
query shards; direct in-memory calls stay log-file free. Use this module for
public annotation API calls; use
``ldsc._kernel.annotation`` only for low-level table and BED primitives.

Key Functions
-------------
AnnotationBuilder :
    Load aligned baseline/query annotation bundles and project BED inputs.
run_bed_to_annot :
    Project BED files with the registered ``GlobalConfig``.
run_annotate_from_args :
    Execute ``ldsc annotate`` from an already parsed namespace.
main :
    Standalone parser entry point for BED-to-annotation projection.

Design Notes
------------
- The workflow layer owns user-facing configuration and output policy.
- The kernel layer receives resolved concrete files and performs primitive
  table/BED operations only.
- Per-run logs are workflow audit artifacts and are not part of returned
  annotation bundles.
"""

from __future__ import annotations

import argparse
import gzip
import logging
import tempfile
from dataclasses import dataclass
from os import PathLike
from pathlib import Path
from typing import Sequence

import numpy as np
import pandas as pd

from ._chr_sampler import sample_frame_from_chr_pattern
from ._coordinates import CHR_POS_KEY_COLUMN, build_chr_pos_key_frame
from ._kernel import annotation as kernel_annotation
from ._kernel.identifiers import (
    build_snp_id_series,
    normalize_snp_identifier_mode,
    validate_unique_snp_ids,
)
from ._kernel.snp_identity import identity_base_mode, identity_mode_family, is_allele_aware_mode
from ._row_alignment import assert_same_snp_rows
from .chromosome_inference import normalize_chromosome
from .column_inference import (
    A1_COLUMN_SPEC,
    A2_COLUMN_SPEC,
    ANNOTATION_METADATA_SPEC_MAP,
    ColumnSpec,
    CHR_COLUMN_SPEC,
    CM_COLUMN_SPEC,
    POS_COLUMN_SPEC,
    SNP_COLUMN_SPEC,
    normalize_genome_build,
    resolve_optional_column,
    resolve_required_column,
)
from .config import (
    AnnotationBuildConfig,
    GlobalConfig,
    get_global_config,
    print_global_config_banner,
)
from .genome_build_inference import resolve_genome_build
from .path_resolution import (
    ANNOTATION_SUFFIXES,
    ensure_output_directory,
    preflight_output_artifact_family,
    remove_output_artifacts,
    resolve_file_group,
    split_cli_path_tokens,
)
from ._logging import configure_package_logging, log_inputs, log_outputs, workflow_logging


LOGGER = logging.getLogger("LDSC.annotation")
_ANNOTATION_A1_COLUMN_SPEC = ColumnSpec(
    A1_COLUMN_SPEC.canonical,
    A1_COLUMN_SPEC.aliases,
    A1_COLUMN_SPEC.label,
    allow_suffix_match=False,
)
_ANNOTATION_A2_COLUMN_SPEC = ColumnSpec(
    A2_COLUMN_SPEC.canonical,
    A2_COLUMN_SPEC.aliases,
    A2_COLUMN_SPEC.label,
    allow_suffix_match=False,
)


@dataclass(frozen=True)
class AnnotationBundle:
    """Aligned SNP metadata plus baseline and query annotation matrices.

    ``AnnotationBundle`` is the public result object shared by the annotation,
    LD-score, and partitioned-regression workflows. All three data frames are
    row-aligned: ``metadata`` carries the SNP keys, ``baseline_annotations``
    carries covariate annotations, and ``query_annotations`` carries optional
    query annotations projected from ``.annot`` or BED inputs.

    Attributes
    ----------
    metadata : pandas.DataFrame
        Row metadata with canonical ``CHR``, ``POS``, ``SNP``, and ``CM``
        columns, plus optional metadata such as ``MAF``.
    baseline_annotations : pandas.DataFrame
        Baseline annotation matrix with one row per metadata row.
    query_annotations : pandas.DataFrame
        Query annotation matrix with one row per metadata row. May have zero
        columns when no query annotations are present.
    baseline_columns, query_columns : list of str
        Ordered annotation column names for the corresponding matrices.
    chromosomes : list of str
        Canonical chromosome labels represented in the bundle.
    source_summary : dict
        Compact provenance summary of the paths or path tokens used to build
        the bundle.
    config_snapshot : GlobalConfig, optional
        Shared identifier/genome-build assumptions active when the bundle was
        built.
    """

    metadata: pd.DataFrame
    baseline_annotations: pd.DataFrame
    query_annotations: pd.DataFrame
    baseline_columns: list[str]
    query_columns: list[str]
    chromosomes: list[str]
    source_summary: dict[str, object]
    config_snapshot: GlobalConfig | None = None

    def validate(self, snp_identifier: str = "chr_pos_allele_aware") -> None:
        """Validate row alignment and uniqueness assumptions for the bundle.

        Parameters
        ----------
        snp_identifier : str, optional
            Canonical identifier mode used for uniqueness checks. Default is
            ``"chr_pos_allele_aware"``.
        """
        snp_identifier = self._annotation_identity_mode(snp_identifier)
        if len(self.metadata) != len(self.baseline_annotations):
            raise ValueError("metadata and baseline_annotations must have the same number of rows.")
        if len(self.metadata) != len(self.query_annotations):
            raise ValueError("metadata and query_annotations must have the same number of rows.")
        if list(self.baseline_annotations.columns) != list(self.baseline_columns):
            raise ValueError("baseline_annotations columns do not match baseline_columns.")
        if list(self.query_annotations.columns) != list(self.query_columns):
            raise ValueError("query_annotations columns do not match query_columns.")
        validate_unique_snp_ids(self.metadata, snp_identifier, context="AnnotationBundle.metadata")

    def reference_snps(self, snp_identifier: str = "chr_pos_allele_aware") -> set[str]:
        """Return the retained reference SNP universe for this bundle.

        Parameters
        ----------
        snp_identifier : str, optional
            Canonical identifier mode used to build the SNP universe. Default is
            ``"chr_pos_allele_aware"``.
        """
        mode = self._annotation_identity_mode(snp_identifier)
        if identity_mode_family(mode) == "rsid":
            return set(build_snp_id_series(self.metadata, mode))
        keyed, _report = build_chr_pos_key_frame(
            self.metadata,
            context="annotation reference SNP universe",
            drop_missing=True,
            logger=LOGGER,
        )
        return set(keyed[CHR_POS_KEY_COLUMN].astype(str))

    def _annotation_identity_mode(self, snp_identifier: str) -> str:
        """Return the effective annotation identity mode for this metadata."""
        mode = normalize_snp_identifier_mode(snp_identifier)
        if is_allele_aware_mode(mode) and not {"A1", "A2"}.issubset(self.metadata.columns):
            return identity_base_mode(mode)
        return mode

    def has_full_baseline_cover(self) -> bool:
        """Return ``True`` when baseline annotations cover every metadata row."""
        return bool(self.baseline_columns) and len(self.metadata) == len(self.baseline_annotations)

    def annotation_matrix(self, include_query: bool = True) -> pd.DataFrame:
        """Return the dense annotation matrix used by LD-score computation."""
        if include_query:
            return pd.concat([self.baseline_annotations, self.query_annotations], axis=1)
        return self.baseline_annotations.copy()

    def summary(self) -> dict[str, object]:
        """Summarize the bundle shape and source provenance for logging."""
        return {
            "n_rows": len(self.metadata),
            "baseline_columns": list(self.baseline_columns),
            "query_columns": list(self.query_columns),
            "chromosomes": list(self.chromosomes),
            "source_summary": dict(self.source_summary),
        }


class AnnotationBuilder:
    """Service object for loading SNP-level annotations and BED projections.

    The builder keeps annotation-side workflow rules in one place: resolving
    input paths, enforcing row alignment across baseline and query tables,
    aggregating chromosome-sharded annotation files, projecting optional BED
    inputs, and returning validated ``AnnotationBundle`` objects.

    Parameters
    ----------
    global_config : GlobalConfig
        Shared identifier, genome-build, and logging settings. ``genome_build``
        must already be normalized to ``None``, ``"auto"``, ``"hg19"``, or
        ``"hg38"``.
    build_config : AnnotationBuildConfig, optional
        Default annotation source configuration for ``run()``.
    """

    def __init__(self, global_config: GlobalConfig, build_config: AnnotationBuildConfig | None = None) -> None:
        """Store shared configuration for bundle loading and BED projection."""
        assert global_config.genome_build in {"auto", "hg19", "hg38", None}, (
            f"genome_build reached AnnotationBuilder as {global_config.genome_build!r}; "
            "should be normalized before annotation building."
        )
        self.global_config = global_config
        self.build_config = build_config or AnnotationBuildConfig()
        self._workflow_log_path: Path | None = None

    def run(self, source_spec: AnnotationBuildConfig | None = None, chrom: str | None = None) -> AnnotationBundle:
        """Build one aligned SNP-level annotation bundle.

        This method accepts whole-genome inputs where all files share one SNP
        row universe, and chromosome-sharded inputs such as
        ``baseline.1.annot.gz`` through ``baseline.22.annot.gz``. BED query
        inputs are projected onto the resolved baseline SNP grid in memory.

        Parameters
        ----------
        source_spec : AnnotationBuildConfig, optional
            Annotation source configuration. Defaults to the builder's
            ``build_config``.
        chrom : str, optional
            Restrict loading to one chromosome. Used by LD-score orchestration
            during chromosome-wise computation.

        Returns
        -------
        bundle : AnnotationBundle
            Validated bundle with aligned metadata, baseline annotations, query
            annotations, chromosome labels, source summary, and config
            snapshot.
        """
        source_spec = source_spec or self.build_config
        if not source_spec.baseline_annot_sources:
            raise ValueError("AnnotationBuildConfig must include at least one baseline annotation file.")
        baseline_files = resolve_file_group(
            source_spec.baseline_annot_sources,
            suffixes=ANNOTATION_SUFFIXES,
            label="baseline annotation",
            allow_chromosome_suite=True,
        )
        query_files = (
            []
            if not source_spec.query_annot_sources
            else resolve_file_group(
                source_spec.query_annot_sources,
                suffixes=ANNOTATION_SUFFIXES,
                label="query annotation",
                allow_chromosome_suite=True,
            )
        )

        sharded_baseline = self._detect_chromosome_shards(baseline_files, group_name="baseline")
        if sharded_baseline is not None:
            sharded_query = self._detect_chromosome_shards(query_files, group_name="query") if query_files else {}
            if query_files and sharded_query is None:
                raise ValueError(
                    "Query annotation chromosome shards must match the baseline chromosome shards exactly."
                )
            return self._run_sharded_inputs(
                source_spec,
                sharded_baseline,
                sharded_query,
                chrom=chrom,
                has_query_inputs=bool(query_files),
            )

        return self._run_single_universe(source_spec, baseline_files, query_files, chrom=chrom)

    def _run_single_universe(
        self,
        source_spec: AnnotationBuildConfig,
        baseline_files: Sequence[str],
        query_files: Sequence[str],
        chrom: str | None = None,
    ) -> AnnotationBundle:
        """Build one bundle when all inputs should share one SNP row universe."""
        metadata: pd.DataFrame | None = None
        baseline_blocks: list[pd.DataFrame] = []
        query_blocks: list[pd.DataFrame] = []
        baseline_columns: list[str] = []
        query_columns: list[str] = []
        seen_columns: set[str] = set()

        for group_name, files in (("baseline", baseline_files), ("query", query_files)):
            for path in files:
                file_metadata, annotations = self.parse_annotation_file(path, chrom=chrom)
                if len(file_metadata) == 0:
                    continue
                if metadata is None:
                    metadata = file_metadata
                else:
                    self._ensure_aligned_rows(metadata, file_metadata, path)
                    metadata = self._merge_missing_metadata(metadata, file_metadata)

                duplicate_columns = sorted(set(annotations.columns) & seen_columns)
                if duplicate_columns:
                    raise ValueError(f"Duplicate annotation column names detected: {duplicate_columns}")
                seen_columns.update(annotations.columns)

                if group_name == "baseline":
                    baseline_blocks.append(annotations)
                    baseline_columns.extend(annotations.columns.tolist())
                else:
                    query_blocks.append(annotations)
                    query_columns.extend(annotations.columns.tolist())

        if metadata is None:
            raise ValueError("No annotation rows were loaded from the supplied sources.")

        if source_spec.query_annot_bed_sources:
            resolved_bed_paths = [
                Path(path) for path in resolve_file_group(source_spec.query_annot_bed_sources, label="BED file")
            ]
            stems = [path.stem for path in resolved_bed_paths]
            duplicate_stems = sorted({stem for stem in stems if stems.count(stem) > 1})
            if duplicate_stems:
                raise ValueError(
                    "BED basenames must be unique because they become annotation column names. "
                    f"Duplicate names: {', '.join(duplicate_stems)}"
                )
            duplicate_columns = sorted(set(stems) & seen_columns)
            if duplicate_columns:
                raise ValueError(
                    "BED file stems clash with loaded annotation column names: "
                    f"{duplicate_columns}"
                )
            with tempfile.TemporaryDirectory(prefix="bed2annot_") as tmpdir:
                tempdir = Path(tmpdir)
                normalized_beds: list[Path] = []
                for bed_path in resolved_bed_paths:
                    normalized_path = tempdir / bed_path.name
                    kernel_annotation._write_normalized_bed(bed_path, normalized_path)
                    normalized_beds.append(normalized_path)
                bed_df = kernel_annotation._compute_bed_query_columns(metadata, normalized_beds, tempdir)
            query_blocks.append(bed_df)
            query_columns.extend(bed_df.columns.tolist())
            seen_columns.update(bed_df.columns)

        baseline_df = self._concat_or_empty(baseline_blocks, index=metadata.index)
        query_df = self._concat_or_empty(query_blocks, index=metadata.index)
        chromosomes = sorted(
            metadata["CHR"].astype(str).map(normalize_chromosome).unique().tolist(),
            key=kernel_annotation._chrom_sort_key,
        )
        return self._build_bundle(
            metadata=metadata,
            baseline_df=baseline_df,
            query_df=query_df,
            baseline_columns=baseline_columns,
            query_columns=query_columns,
            chromosomes=chromosomes,
            source_spec=source_spec,
        )

    def _run_sharded_inputs(
        self,
        source_spec: AnnotationBuildConfig,
        baseline_by_chrom: dict[str, str],
        query_by_chrom: dict[str, str],
        chrom: str | None,
        has_query_inputs: bool,
    ) -> AnnotationBundle:
        """Build and aggregate one bundle per chromosome shard."""
        if chrom is not None:
            chrom_key = normalize_chromosome(chrom)
            baseline_by_chrom = {chrom_key: baseline_by_chrom[chrom_key]} if chrom_key in baseline_by_chrom else {}
            query_by_chrom = {chrom_key: query_by_chrom[chrom_key]} if chrom_key in query_by_chrom else {}

        if not baseline_by_chrom:
            raise ValueError("No annotation rows were loaded from the supplied sources.")
        if has_query_inputs and set(query_by_chrom) != set(baseline_by_chrom):
            raise ValueError(
                "Query annotation chromosome shards must match the baseline chromosome shards exactly."
            )

        bundles: list[AnnotationBundle] = []
        for chrom_key in sorted(baseline_by_chrom, key=kernel_annotation._chrom_sort_key):
            chrom_source_spec = AnnotationBuildConfig(
                baseline_annot_sources=(baseline_by_chrom[chrom_key],),
                query_annot_sources=(() if not has_query_inputs else (query_by_chrom[chrom_key],)),
                query_annot_bed_sources=source_spec.query_annot_bed_sources,
                allow_missing_query=source_spec.allow_missing_query,
            )
            bundles.append(
                self._run_single_universe(
                    chrom_source_spec,
                    (baseline_by_chrom[chrom_key],),
                    (() if not has_query_inputs else (query_by_chrom[chrom_key],)),
                    chrom=chrom_key,
                )
            )

        baseline_columns = bundles[0].baseline_columns
        query_columns = bundles[0].query_columns
        for bundle in bundles[1:]:
            if bundle.baseline_columns != baseline_columns or bundle.query_columns != query_columns:
                raise ValueError(
                    "Per-chromosome annotation shards must expose the same baseline and query annotation columns."
                )

        metadata = pd.concat([bundle.metadata for bundle in bundles], axis=0, ignore_index=True)
        baseline_df = pd.concat([bundle.baseline_annotations for bundle in bundles], axis=0, ignore_index=True)
        query_df = pd.concat([bundle.query_annotations for bundle in bundles], axis=0, ignore_index=True)
        chromosomes = [chrom_key for chrom_key in sorted(baseline_by_chrom, key=kernel_annotation._chrom_sort_key)]
        return self._build_bundle(
            metadata=metadata,
            baseline_df=baseline_df,
            query_df=query_df,
            baseline_columns=baseline_columns,
            query_columns=query_columns,
            chromosomes=chromosomes,
            source_spec=source_spec,
        )

    def _build_bundle(
        self,
        metadata: pd.DataFrame,
        baseline_df: pd.DataFrame,
        query_df: pd.DataFrame,
        baseline_columns: list[str],
        query_columns: list[str],
        chromosomes: list[str],
        source_spec: AnnotationBuildConfig,
    ) -> AnnotationBundle:
        """Construct, validate, and return an ``AnnotationBundle``."""
        bundle = AnnotationBundle(
            metadata=metadata.reset_index(drop=True),
            baseline_annotations=baseline_df.reset_index(drop=True),
            query_annotations=query_df.reset_index(drop=True),
            baseline_columns=baseline_columns,
            query_columns=query_columns,
            chromosomes=chromosomes,
            source_summary={
                "baseline_annot_sources": list(source_spec.baseline_annot_sources),
                "query_annot_sources": list(source_spec.query_annot_sources),
                "query_annot_bed_sources": list(source_spec.query_annot_bed_sources),
            },
            config_snapshot=self.global_config,
        )
        bundle.validate(self.global_config.snp_identifier)
        return bundle

    def _detect_chromosome_shards(self, files: Sequence[str], group_name: str) -> dict[str, str] | None:
        """Return a chromosome-to-file mapping when inputs are cleanly sharded."""
        if not files:
            return None

        shard_map: dict[str, str] = {}
        saw_unsharded = False
        saw_sharded = False
        for path in files:
            chrom = kernel_annotation._annotation_shard_chromosome(path)
            if chrom is None:
                saw_unsharded = True
                continue
            saw_sharded = True
            if chrom in shard_map:
                raise ValueError(
                    f"Ambiguous {group_name} annotation chromosome shards: multiple files map to chromosome {chrom}."
                )
            shard_map[chrom] = path

        if saw_sharded and saw_unsharded:
            raise ValueError(
                f"Mixed whole-genome and per-chromosome {group_name} annotation inputs are not supported."
            )
        if saw_sharded:
            return shard_map
        return None

    def parse_annotation_file(self, path: str | Path, chrom: str | None = None) -> tuple[pd.DataFrame, pd.DataFrame]:
        """Read one SNP-level annotation file into normalized metadata and values."""
        df = kernel_annotation._read_text_table(str(path))
        context = str(path)
        chr_col = resolve_required_column(df.columns, CHR_COLUMN_SPEC, context=context)
        pos_col = resolve_required_column(df.columns, POS_COLUMN_SPEC, context=context)
        snp_col = resolve_required_column(df.columns, SNP_COLUMN_SPEC, context=context)
        cm_col = resolve_required_column(df.columns, CM_COLUMN_SPEC, context=context)
        a1_col = resolve_optional_column(df.columns, _ANNOTATION_A1_COLUMN_SPEC, context=context)
        a2_col = resolve_optional_column(df.columns, _ANNOTATION_A2_COLUMN_SPEC, context=context)
        if (a1_col is None) ^ (a2_col is None):
            raise ValueError("Annotation file has only one allele column; provide both A1 and A2 or neither.")
        metadata = pd.DataFrame(
            {
                "CHR": df[chr_col],
                "POS": df[pos_col],
                "SNP": df[snp_col],
                "CM": df[cm_col],
            }
        )
        metadata["CHR"] = metadata["CHR"].map(lambda value: normalize_chromosome(value, context=context))
        metadata["POS"] = pd.to_numeric(metadata["POS"], errors="raise").astype(np.int64)
        metadata["SNP"] = metadata["SNP"].astype(str)
        metadata["CM"] = pd.to_numeric(metadata["CM"], errors="coerce")
        if a1_col is not None and a2_col is not None:
            metadata["A1"] = df[a1_col]
            metadata["A2"] = df[a2_col]
        maf_col = resolve_optional_column(df.columns, ANNOTATION_METADATA_SPEC_MAP["MAF"], context=context)
        if maf_col is not None:
            metadata["MAF"] = pd.to_numeric(df[maf_col], errors="coerce")
        if chrom is not None:
            keep = metadata["CHR"] == normalize_chromosome(chrom, context=context)
            metadata = metadata.loc[keep].reset_index(drop=True)
            df = df.loc[keep].reset_index(drop=True)

        if len(metadata) == 0:
            return metadata, pd.DataFrame(index=metadata.index)

        metadata_source_columns = {chr_col, pos_col, snp_col, cm_col, maf_col, a1_col, a2_col}
        annotation_columns = [column for column in df.columns if column not in metadata_source_columns]
        if not annotation_columns:
            raise ValueError(f"{path} does not contain any annotation columns.")
        annotations = df.loc[:, annotation_columns].astype(np.float32).reset_index(drop=True)
        metadata = metadata.reset_index(drop=True)
        return metadata, annotations

    def project_bed_annotations(
        self,
        query_annot_bed_sources: str | PathLike[str] | Sequence[str | PathLike[str]],
        baseline_annot_sources: str | PathLike[str] | Sequence[str | PathLike[str]],
        output_dir: str | Path | None = None,
        log_level: str | None = None,
        overwrite: bool | None = None,
    ) -> AnnotationBundle:
        """Convert BED annotations into an in-memory ``AnnotationBundle``.

        When ``output_dir`` is provided, projected query annotations are also
        written as per-chromosome ``query.<chrom>.annot.gz`` files. Existing
        root-level ``query.*.annot.gz`` siblings are refused before writing
        unless overwrite is true; successful overwrites remove stale query
        shards outside the current chromosome set. Public workflow wrappers may
        attach a private log path so the same preflight also protects
        ``annotate.log``; ordinary direct calls do not create a log file.

        Parameters
        ----------
        query_annot_bed_sources : str, path-like, or sequence
            BED input token or tokens. Each resolved BED basename becomes one
            query annotation column.
        baseline_annot_sources : str, path-like, or sequence
            Baseline ``.annot(.gz)`` input token or tokens that define the SNP
            grid used for BED projection.
        output_dir : str or pathlib.Path, optional
            Directory for materialized ``query.<chrom>.annot.gz`` files. If
            omitted, projection stays in memory.
        log_level : str, optional
            Logging threshold for standalone projection calls. When the method
            is reached from the parsed workflow wrapper, the same threshold
            controls records written to ``annotate.log``.
        overwrite : bool, optional
            Whether fixed output files may be replaced and stale owned siblings
            removed. Defaults to the builder configuration when omitted.

        Returns
        -------
        bundle : AnnotationBundle
            Projected annotation bundle. The same object is returned whether
            or not query annotation files are written.
        """
        overwrite = self.build_config.overwrite if overwrite is None else overwrite
        _configure_logging(log_level or self.global_config.log_level)
        source_spec = AnnotationBuildConfig(
            baseline_annot_sources=baseline_annot_sources,
            query_annot_bed_sources=query_annot_bed_sources,
        )
        bundle = self.run(source_spec)
        if output_dir is not None:
            output_path = ensure_output_directory(output_dir, label="output directory")
            output_paths = _bundle_query_annot_output_paths(bundle, output_path)
            produced_paths: list[Path] = list(output_paths)
            owned_paths = sorted(output_path.glob("query.*.annot.gz"))
            owned_paths.extend(produced_paths)
            if self._workflow_log_path is not None:
                produced_paths.append(self._workflow_log_path)
                owned_paths.append(self._workflow_log_path)
            stale_paths = preflight_output_artifact_family(
                produced_paths,
                owned_paths,
                overwrite=overwrite,
                label="annotation output artifact",
            )
            with workflow_logging("annotate", self._workflow_log_path, log_level=log_level or self.global_config.log_level):
                log_inputs(
                    query_annot_bed_sources=query_annot_bed_sources,
                    baseline_annot_sources=baseline_annot_sources,
                    output_dir=str(output_path),
                )
                written = _write_bundle_query_as_annot_files(bundle, output_path)
                log_outputs(**{f"query_{chrom}": str(path) for chrom, path in zip(bundle.chromosomes, written)})
                remove_output_artifacts(stale_paths)
        return bundle

    def _ensure_aligned_rows(self, reference: pd.DataFrame, current: pd.DataFrame, path: str) -> None:
        """Raise if two annotation tables do not share identical SNP row order."""
        snp_identifier = self.global_config.snp_identifier
        if not ({"A1", "A2"}.issubset(reference.columns) and {"A1", "A2"}.issubset(current.columns)):
            snp_identifier = identity_base_mode(snp_identifier)
        assert_same_snp_rows(
            reference,
            current,
            context=f"Annotation SNP rows do not match across files: {path}",
            snp_identifier=snp_identifier,
        )

    def _merge_missing_metadata(self, reference: pd.DataFrame, current: pd.DataFrame) -> pd.DataFrame:
        """Backfill missing ``CM`` or ``MAF`` metadata from another aligned table."""
        merged = reference.copy()
        if "MAF" in current.columns:
            if "MAF" not in merged.columns:
                merged["MAF"] = np.nan
            missing_maf = merged["MAF"].isna() & current["MAF"].notna()
            if missing_maf.any():
                merged.loc[missing_maf, "MAF"] = current.loc[missing_maf, "MAF"].to_numpy()
        missing_cm = merged["CM"].isna() & current["CM"].notna()
        if missing_cm.any():
            merged.loc[missing_cm, "CM"] = current.loc[missing_cm, "CM"].to_numpy()
        return merged

    def _concat_or_empty(self, frames: list[pd.DataFrame], index: pd.Index) -> pd.DataFrame:
        """Concatenate annotation frames or return an empty frame on ``index``."""
        if frames:
            return pd.concat(frames, axis=1)
        return pd.DataFrame(index=index)


def run_bed_to_annot(
    query_annot_bed_sources: str | PathLike[str] | Sequence[str | PathLike[str]],
    baseline_annot_sources: str | PathLike[str] | Sequence[str | PathLike[str]],
    output_dir: str | Path | None = None,
    overwrite: bool = False,
) -> AnnotationBundle:
    """Project BED files to an ``AnnotationBundle`` using registered globals.

    Python workflows read shared identifier, genome-build, and logging settings
    from the registered ``GlobalConfig``. When ``output_dir`` is supplied,
    fixed ``query.<chrom>.annot.gz`` outputs are refused before writing unless
    ``overwrite=True``. With overwrite enabled, stale query shards from earlier
    chromosome sets are removed after successful writes. The wrapper also
    writes ``annotate.log`` in the same directory; the returned bundle remains
    an in-memory data object and does not expose the log path.

    Parameters
    ----------
    query_annot_bed_sources : str, path-like, or sequence
        BED input token or tokens. Exact paths, globs, and comma-separated CLI
        token fragments are normalized by the workflow.
    baseline_annot_sources : str, path-like, or sequence
        Baseline annotation token or tokens defining the SNP grid.
    output_dir : str or pathlib.Path, optional
        Destination directory for generated query annotation shards. Omit to
        keep the result in memory only.
    overwrite : bool, optional
        Replace existing ``query.<chrom>.annot.gz`` files when true. Default is
        ``False``.

    Returns
    -------
    bundle : AnnotationBundle
        Projected annotation bundle.
    """
    return _run_bed_to_annot_with_global_config(
        query_annot_bed_sources=query_annot_bed_sources,
        baseline_annot_sources=baseline_annot_sources,
        output_dir=output_dir,
        overwrite=overwrite,
        global_config=get_global_config(),
        entrypoint="run_bed_to_annot",
    )


def _run_bed_to_annot_with_global_config(
    query_annot_bed_sources: str | PathLike[str] | Sequence[str | PathLike[str]],
    baseline_annot_sources: str | PathLike[str] | Sequence[str | PathLike[str]],
    output_dir: str | Path | None,
    *,
    overwrite: bool,
    global_config: GlobalConfig,
    entrypoint: str,
) -> AnnotationBundle:
    """Run BED-to-annotation projection with an explicit resolved GlobalConfig."""
    print_global_config_banner(entrypoint, global_config)
    build_config = AnnotationBuildConfig(query_annot_bed_sources=query_annot_bed_sources)
    builder = AnnotationBuilder(global_config, build_config)
    if output_dir is not None:
        builder._workflow_log_path = Path(output_dir) / "annotate.log"
    return builder.project_bed_annotations(
        query_annot_bed_sources=query_annot_bed_sources,
        baseline_annot_sources=baseline_annot_sources,
        output_dir=output_dir,
        log_level=global_config.log_level,
        overwrite=overwrite,
    )


def add_annotate_arguments(parser: argparse.ArgumentParser) -> None:
    """Register ``ldsc annotate`` arguments on an existing parser.

    This is the shared argument definition used by the unified ``ldsc`` parser
    and by the standalone annotation parser built in this module.
    """
    parser.add_argument(
        "--query-annot-bed-sources",
        nargs="+",
        required=True,
        help="BED files, comma-separated lists, or glob patterns.",
    )
    parser.add_argument(
        "--baseline-annot-sources",
        nargs="+",
        required=True,
        help="Baseline annotation path tokens: exact paths, globs, or explicit @ suite tokens.",
    )
    parser.add_argument(
        "--output-dir",
        required=True,
        help="Destination directory for generated .annot.gz files.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        default=False,
        help="Replace annotation output artifacts and remove stale owned query shards.",
    )
    parser.add_argument(
        "--snp-identifier",
        default="chr_pos_allele_aware",
        choices=("rsid", "rsid_allele_aware", "chr_pos", "chr_pos_allele_aware"),
        help="Identifier mode used for bundle validation.",
    )
    parser.add_argument(
        "--genome-build",
        default=None,
        choices=("auto", "hg19", "hg37", "GRCh37", "hg38", "GRCh38"),
        help=(
            "Genome build for chr_pos-family inputs. Required when --snp-identifier is "
            "chr_pos or chr_pos_allele_aware. Use 'auto' to infer hg19/hg38 and "
            "0-based/1-based coordinates from data. Not used for rsid-family modes."
        ),
    )
    parser.add_argument("--log-level", default="INFO", choices=("DEBUG", "INFO", "WARNING", "ERROR"), help="Logging verbosity.")


def build_parser() -> argparse.ArgumentParser:
    """Build the standalone parser for the annotation workflow.

    Returns
    -------
    parser : argparse.ArgumentParser
        Parser for arguments accepted after ``ldsc annotate``.
    """
    parser = argparse.ArgumentParser(description="Convert BED annotations into LDSC .annot.gz files.", allow_abbrev=False)
    add_annotate_arguments(parser)
    return parser


def parse_bed_to_annot_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    """Parse CLI arguments for the BED-to-annotation workflow.

    Parameters
    ----------
    argv : sequence of str, optional
        Argument vector excluding the program name. ``None`` lets argparse read
        from ``sys.argv``.

    Returns
    -------
    args : argparse.Namespace
        Parsed annotation workflow arguments.
    """
    return build_parser().parse_args(argv)


def run_annotate_from_args(args: argparse.Namespace) -> AnnotationBundle:
    """Run the annotation workflow from a parsed CLI namespace.

    This function is the workflow-layer dispatch target used by the unified
    ``ldsc`` parser. It accepts the already parsed namespace, normalizes
    identifier/build settings, splits CLI path tokens, resolves
    ``--genome-build auto`` when requested, writes ``query.<chrom>.annot.gz``
    files, and returns the produced bundle without reparsing arguments.

    Parameters
    ----------
    args : argparse.Namespace
        Namespace produced by ``ldsc.cli.build_parser()`` or
        ``parse_bed_to_annot_args()``.

    Returns
    -------
    bundle : AnnotationBundle
        Produced annotation bundle.
    """
    if not getattr(args, "query_annot_bed_sources", None):
        raise ValueError("ldsc annotate requires --query-annot-bed-sources.")
    if not getattr(args, "baseline_annot_sources", None):
        raise ValueError("ldsc annotate requires --baseline-annot-sources.")
    if not getattr(args, "output_dir", None):
        raise ValueError("ldsc annotate requires --output-dir.")

    normalized_mode = normalize_snp_identifier_mode(args.snp_identifier)
    genome_build = _resolve_annotation_cli_genome_build(args, normalized_mode)
    cli_global_config = GlobalConfig(
        snp_identifier=normalized_mode,
        genome_build=genome_build,
        log_level=args.log_level,
    )
    return _run_bed_to_annot_with_global_config(
        query_annot_bed_sources=split_cli_path_tokens(args.query_annot_bed_sources),
        baseline_annot_sources=split_cli_path_tokens(args.baseline_annot_sources),
        output_dir=args.output_dir,
        overwrite=args.overwrite,
        global_config=cli_global_config,
        entrypoint="run_annotate_from_args",
    )


def main(argv: Sequence[str] | None = None) -> AnnotationBundle:
    """Parse annotation arguments, run projection, and return the bundle.

    ``main()`` replaces the removed ``main_bed_to_annot()`` public entry point.
    It is used by direct ``ldsc annotate`` dispatch and by script-style
    invocations that want parser behavior plus the produced result object.
    """
    return run_annotate_from_args(parse_bed_to_annot_args(argv))


def _resolve_annotation_cli_genome_build(args: argparse.Namespace, snp_identifier: str) -> str | None:
    """Resolve the effective genome build for annotation CLI inputs."""
    if identity_mode_family(snp_identifier) == "rsid":
        return normalize_genome_build(args.genome_build)
    genome_build = normalize_genome_build(args.genome_build)
    if genome_build is None:
        raise ValueError(
            "genome_build is required for chr_pos-family snp_identifier modes. "
            "Pass --genome-build auto, --genome-build hg19, or --genome-build hg38."
        )
    if genome_build == "auto":
        frame, sampled_path = sample_frame_from_chr_pattern(
            split_cli_path_tokens(args.baseline_annot_sources),
            context="annotation inputs",
        )
        genome_build = resolve_genome_build(
            "auto",
            "chr_pos",
            frame,
            context="annotation inputs",
            logger=LOGGER,
        )
        LOGGER.info(f"Resolved annotation genome build from '{sampled_path}'.")
    return genome_build


def _configure_logging(level: str) -> None:
    """Configure logging for standalone annotation entry points."""
    configure_package_logging(level)


def _write_bundle_query_as_annot_files(bundle: AnnotationBundle, output_dir: Path) -> list[Path]:
    """Write one ``query.<chrom>.annot.gz`` file per chromosome from ``bundle``."""
    output_paths: list[Path] = []
    for chrom, out_path in zip(bundle.chromosomes, _bundle_query_annot_output_paths(bundle, output_dir)):
        chrom_mask = bundle.metadata["CHR"].astype(str) == str(chrom)
        chrom_meta = bundle.metadata.loc[chrom_mask].reset_index(drop=True).copy()
        chrom_query = bundle.query_annotations.loc[chrom_mask].reset_index(drop=True)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        with gzip.open(out_path, "wt") as handle:
            pd.concat([chrom_meta, chrom_query], axis=1).to_csv(handle, sep="\t", index=False)
        output_paths.append(out_path)
    return output_paths


def _bundle_query_annot_output_paths(bundle: AnnotationBundle, output_dir: Path) -> list[Path]:
    """Return fixed query annotation output paths for one bundle."""
    return [output_dir / f"query.{chrom}.annot.gz" for chrom in bundle.chromosomes]


__all__ = [
    "AnnotationBuildConfig",
    "AnnotationBuilder",
    "AnnotationBundle",
    "add_annotate_arguments",
    "build_parser",
    "main",
    "parse_bed_to_annot_args",
    "run_annotate_from_args",
    "run_bed_to_annot",
]
