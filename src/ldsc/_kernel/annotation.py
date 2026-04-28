"""
annotation.py

Core functionality:
    Load SNP-level annotation tables into validated ``AnnotationBundle``
    objects and project region-level BED inputs into LDSC ``.annot.gz`` files.

Overview
--------
This module is the annotation-side execution kernel for the refactored codebase. It owns
two related jobs: reading already materialized SNP-level annotation tables into
an ``AnnotationBundle`` that the LD-score workflow can consume, and converting
BED intervals into chromosome-matched ``.annot.gz`` outputs using baseline
annotation templates. The loading path supports both single-table
inputs, where every file must share the same SNP row universe, and LDSC-style
chromosome-sharded inputs such as ``baseline.1.annot.gz`` through
``baseline.22.annot.gz``, which are validated per chromosome and then
aggregated into one bundle.

Key Functions
-------------
AnnotationBuilder.run :
    Build one aligned SNP-level annotation bundle from baseline, query, and
    in-memory BED-projected annotation sources, with automatic aggregation of
    chromosome-sharded files.
run_bed_to_annot :
    Project one or more BED files onto baseline annotation templates and
    return the resulting in-memory bundle, optionally writing `.annot.gz`
    files for callers that explicitly request them.

Design Notes
------------
- All annotation tables are normalized to the same required metadata columns:
  ``CHR``, ``POS``, ``SNP``, and ``CM``.
- When filenames cleanly encode one chromosome shard each, bundles are built
  independently per chromosome and concatenated in stable genomic order.
- BED projection relies on ``pybedtools`` for interval intersection but keeps
  row-order validation explicit so the generated masks stay aligned to the
  baseline SNP template.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import logging
import os
import re
import sys
import tempfile
from dataclasses import dataclass, field
from os import PathLike
from pathlib import Path
from typing import Iterator, Sequence

import numpy as np
import pandas as pd

from ..chromosome_inference import chrom_sort_key, normalize_chromosome
from ..column_inference import (
    ANNOTATION_METADATA_SPEC_MAP,
    CHR_COLUMN_SPEC,
    CM_COLUMN_SPEC,
    RESTRICTION_CHRPOS_SPEC_MAP,
    RESTRICTION_RSID_SPEC_MAP,
    POS_COLUMN_SPEC,
    SNP_COLUMN_SPEC,
    normalize_genome_build,
    resolve_optional_column,
    resolve_required_column,
)
from ..config import (
    AnnotationBuildConfig,
    GlobalConfig,
    get_global_config,
    print_global_config_banner,
    _normalize_path_tuple,
)
from .._chr_sampler import sample_frame_from_chr_pattern
from ..genome_build_inference import resolve_genome_build
from ..path_resolution import (
    ANNOTATION_SUFFIXES,
    ensure_output_directory,
    ensure_output_paths_available,
    ensure_output_parent_directory,
    normalize_path_token,
    resolve_file_group,
    resolve_scalar_path,
    split_cli_path_tokens,
)
from .identifiers import (
    build_snp_id_series,
    normalize_snp_identifier_mode,
    read_global_snp_restriction,
    validate_unique_snp_ids,
)


LOGGER = logging.getLogger("LDSC.annotation")
REQUIRED_ANNOT_COLUMNS = ("CHR", "POS", "SNP", "CM")


@dataclass(frozen=True)
class AnnotationBundle:
    """Aligned SNP metadata plus baseline and query annotation matrices."""
    metadata: pd.DataFrame
    baseline_annotations: pd.DataFrame
    query_annotations: pd.DataFrame
    baseline_columns: list[str]
    query_columns: list[str]
    chromosomes: list[str]
    source_summary: dict[str, object]
    config_snapshot: GlobalConfig | None = None

    def validate(self, snp_identifier: str = "chr_pos") -> None:
        """Validate row alignment and uniqueness assumptions for the bundle.

        Parameters
        ----------
        snp_identifier : str, optional
            Canonical identifier mode used for uniqueness checks. Default is ``"chr_pos"``.
        """
        snp_identifier = normalize_snp_identifier_mode(snp_identifier)
        if len(self.metadata) != len(self.baseline_annotations):
            raise ValueError("metadata and baseline_annotations must have the same number of rows.")
        if len(self.metadata) != len(self.query_annotations):
            raise ValueError("metadata and query_annotations must have the same number of rows.")
        if list(self.baseline_annotations.columns) != list(self.baseline_columns):
            raise ValueError("baseline_annotations columns do not match baseline_columns.")
        if list(self.query_annotations.columns) != list(self.query_columns):
            raise ValueError("query_annotations columns do not match query_columns.")
        validate_unique_snp_ids(self.metadata, snp_identifier, context="AnnotationBundle.metadata")

    def reference_snps(self, snp_identifier: str = "chr_pos") -> set[str]:
        """Return the retained reference SNP universe for this bundle.

        Parameters
        ----------
        snp_identifier : str, optional
            Canonical identifier mode used to build the SNP universe. Default is ``"chr_pos"``.
        """
        return set(build_snp_id_series(self.metadata, snp_identifier))

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


@dataclass(frozen=True)
class _BaselineRow:
    """Minimal metadata for one baseline SNP row during BED projection."""
    chrom: str
    pos: int
    snp: str
    cm: str


@dataclass(frozen=True)
class _RestrictResource:
    """Normalized representation of a reference-universe SNP filter input."""
    mode: str
    bed_path: Path | None = None
    snp_ids: frozenset[str] | None = None


class AnnotationBuilder:
    """
    Service object for loading SNP-level annotation inputs and BED projections.

    The builder keeps the annotation-side validation rules in one place:
    resolving input paths, enforcing row alignment across baseline and query
    tables, automatically aggregating chromosome-sharded annotation files when
    their filenames encode one chromosome each, and applying the global SNP
        the retained reference-panel SNP universe consistently.
    """

    def __init__(self, global_config: GlobalConfig, build_config: AnnotationBuildConfig | None = None) -> None:
        """Store shared configuration for bundle loading and BED projection."""
        assert global_config.genome_build in {"hg19", "hg38", None}, (
            f"genome_build reached kernel as {global_config.genome_build!r}; "
            "should have been resolved at workflow entry."
        )
        self.global_config = global_config
        self.build_config = build_config or AnnotationBuildConfig()

    def run(self, source_spec: AnnotationBuildConfig | None = None, chrom: str | None = None) -> AnnotationBundle:
        """
        Build one aligned SNP-level annotation bundle.

        This method accepts two input layouts. In the single-universe layout,
        every baseline and query file is interpreted as another annotation block
        over the same SNP rows, so all files must have identical SNP identifiers
        in identical order. In the chromosome-sharded layout, resolved filenames
        such as ``baseline.1.annot.gz`` and ``query.1.annot.gz`` are grouped by
        chromosome, each chromosome is validated independently, and the per-
        chromosome bundles are concatenated into one aggregate bundle.

        Parameters
        ----------
        source_spec : AnnotationBuildConfig
            Baseline and query annotation sources to load.
        chrom : str, optional
            Restrict loading to one chromosome before alignment. Default is
            ``None``, which keeps all chromosomes present in the input files.

        Returns
        -------
        AnnotationBundle
            Metadata plus dense baseline/query annotation blocks. For
            single-universe inputs, all rows come from one shared SNP universe.
            For chromosome-sharded inputs, rows are aggregated across the
            detected chromosome shards in stable chromosome order.

        Raises
        ------
        ValueError
            If no baseline annotation files are provided, if row alignment fails
            within one annotation universe, if shard detection reveals a mixed
            whole-genome/per-chromosome layout, if multiple files map
            ambiguously to the same chromosome shard, or if baseline and query
            chromosome shards do not cover the same chromosomes.
        """
        source_spec = source_spec or self.build_config
        if not source_spec.baseline_annot_paths:
            raise ValueError("AnnotationBuildConfig must include at least one baseline annotation file.")
        baseline_files = resolve_file_group(
            source_spec.baseline_annot_paths,
            suffixes=ANNOTATION_SUFFIXES,
            label="baseline annotation",
            allow_chromosome_suite=True,
        )
        query_files = (
            []
            if not source_spec.query_annot_paths
            else resolve_file_group(
                source_spec.query_annot_paths,
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

        if source_spec.query_annot_bed_paths:
            resolved_bed_paths = [Path(path) for path in resolve_file_group(source_spec.query_annot_bed_paths, label="BED file")]
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
                    _write_normalized_bed(bed_path, normalized_path)
                    normalized_beds.append(normalized_path)
                bed_df = _compute_bed_query_columns(metadata, normalized_beds, tempdir)
            query_blocks.append(bed_df)
            query_columns.extend(bed_df.columns.tolist())
            seen_columns.update(bed_df.columns)

        baseline_df = self._concat_or_empty(baseline_blocks, index=metadata.index)
        query_df = self._concat_or_empty(query_blocks, index=metadata.index)
        chromosomes = sorted(metadata["CHR"].astype(str).map(normalize_chromosome).unique().tolist(), key=_chrom_sort_key)
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
        for chrom_key in sorted(baseline_by_chrom, key=_chrom_sort_key):
            chrom_source_spec = AnnotationBuildConfig(
                baseline_annot_paths=(baseline_by_chrom[chrom_key],),
                query_annot_paths=(() if not has_query_inputs else (query_by_chrom[chrom_key],)),
                query_annot_bed_paths=source_spec.query_annot_bed_paths,
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
        chromosomes = [chrom_key for chrom_key in sorted(baseline_by_chrom, key=_chrom_sort_key)]
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
                "baseline_annot_paths": list(source_spec.baseline_annot_paths),
                "query_annot_paths": list(source_spec.query_annot_paths),
                "query_annot_bed_paths": list(source_spec.query_annot_bed_paths),
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
            chrom = _annotation_shard_chromosome(path)
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
        df = _read_text_table(str(path))
        context = str(path)
        chr_col = resolve_required_column(df.columns, CHR_COLUMN_SPEC, context=context)
        pos_col = resolve_required_column(df.columns, POS_COLUMN_SPEC, context=context)
        snp_col = resolve_required_column(df.columns, SNP_COLUMN_SPEC, context=context)
        cm_col = resolve_required_column(df.columns, CM_COLUMN_SPEC, context=context)
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
        maf_col = resolve_optional_column(df.columns, ANNOTATION_METADATA_SPEC_MAP["MAF"], context=context)
        if maf_col is not None:
            metadata["MAF"] = pd.to_numeric(df[maf_col], errors="coerce")
        if chrom is not None:
            keep = metadata["CHR"] == normalize_chromosome(chrom, context=context)
            metadata = metadata.loc[keep].reset_index(drop=True)
            df = df.loc[keep].reset_index(drop=True)

        if len(metadata) == 0:
            return metadata, pd.DataFrame(index=metadata.index)

        annotation_columns = [
            column for column in df.columns if column not in {chr_col, pos_col, snp_col, cm_col, "MAF"}
        ]
        if not annotation_columns:
            raise ValueError(f"{path} does not contain any annotation columns.")
        annotations = df.loc[:, annotation_columns].astype(np.float32).reset_index(drop=True)
        metadata = metadata.reset_index(drop=True)
        return metadata, annotations

    def project_bed_annotations(
        self,
        query_annot_bed_paths: str | PathLike[str] | Sequence[str | PathLike[str]],
        baseline_annot_paths: str | PathLike[str] | Sequence[str | PathLike[str]],
        output_dir: str | Path | None = None,
        batch: bool | None = None,
        log_level: str | None = None,
        overwrite: bool | None = None,
    ) -> AnnotationBundle:
        """
        Convert BED annotations into an in-memory ``AnnotationBundle``.

        When ``output_dir`` is provided, the projected query annotations are
        also written as per-chromosome ``query.<chrom>.annot.gz`` files. The
        ``batch`` argument is retained for compatibility but no longer changes
        behavior: all BED inputs become query columns on one bundle. This
        annotation-side projection does not apply any reference-panel SNP
        restriction; those filters are applied later when the reference panel
        is loaded for LD-score computation.

        Existing ``query.<chrom>.annot.gz`` files are checked after the bundle
        is assembled and before any shard is written. Replacement requires the
        ``overwrite`` argument, or ``AnnotationBuildConfig(overwrite=True)``
        when ``overwrite`` is left as ``None``.
        """
        _ = self.build_config.batch_mode if batch is None else batch
        overwrite = self.build_config.overwrite if overwrite is None else overwrite
        _configure_logging(log_level or self.global_config.log_level)
        source_spec = AnnotationBuildConfig(
            baseline_annot_paths=baseline_annot_paths,
            query_annot_bed_paths=query_annot_bed_paths,
        )
        bundle = self.run(source_spec)
        if output_dir is not None:
            output_path = ensure_output_directory(output_dir, label="output directory")
            ensure_output_paths_available(
                _bundle_query_annot_output_paths(bundle, output_path),
                overwrite=overwrite,
                label="annotation output artifact",
            )
            _write_bundle_query_as_annot_files(bundle, output_path)
        return bundle

    def make_single_annotation_file(
        self,
        bimfile: str | Path,
        annot_file: str | Path,
        bed_for_annot,
    ) -> Path:
        """Project one BED-like input onto one BIM file and write a legacy `.annot`."""
        pybedtools = _get_pybedtools()
        bimfile = Path(resolve_scalar_path(bimfile, label="PLINK BIM file"))
        annot_file = ensure_output_parent_directory(annot_file, label="annot_file")

        df_bim = pd.read_csv(
            bimfile,
            sep=r"\s+",
            usecols=[0, 1, 2, 3],
            names=["CHR", "SNP", "CM", "POS"],
            header=None,
        )
        iter_bim = [[_to_bed_chromosome(chrom), int(pos) - 1, int(pos)] for chrom, pos in np.array(df_bim[["CHR", "POS"]])]
        bim_bed = pybedtools.BedTool(iter_bim)
        annot_bed = bim_bed.intersect(bed_for_annot)
        pos = [feature.start + 1 for feature in annot_bed]
        df_int = pd.DataFrame({"POS": pos, "ANNOT": 1})
        df_annot = pd.merge(df_bim, df_int, how="left", on="POS")
        df_annot.fillna(0, inplace=True)
        df_annot = df_annot[["ANNOT"]].astype(int)
        annot_file.parent.mkdir(parents=True, exist_ok=True)
        if str(annot_file).endswith(".gz"):
            with gzip.open(annot_file, "wt") as handle:
                df_annot.to_csv(handle, sep="\t", index=False)
        else:
            df_annot.to_csv(annot_file, sep="\t", index=False)
        return annot_file

    def _resolve_annotation_path(self, path: str | Path) -> str:
        """Resolve one annotation token to exactly one concrete annotation file."""
        return resolve_scalar_path(path, suffixes=ANNOTATION_SUFFIXES, label="annotation")

    def _ensure_aligned_rows(self, reference: pd.DataFrame, current: pd.DataFrame, path: str) -> None:
        """Raise if two annotation tables do not share identical SNP row order."""
        ref_keys = build_snp_id_series(reference, self.global_config.snp_identifier)
        cur_keys = build_snp_id_series(current, self.global_config.snp_identifier)
        if len(reference) != len(current) or not ref_keys.equals(cur_keys):
            raise ValueError(f"Annotation SNP rows do not match across files: {path}")
        if not reference.loc[:, ["CHR", "POS", "SNP"]].equals(current.loc[:, ["CHR", "POS", "SNP"]]):
            raise ValueError(f"Annotation metadata mismatch across files: {path}")

    def _merge_missing_metadata(self, reference: pd.DataFrame, current: pd.DataFrame) -> pd.DataFrame:
        """Backfill missing `CM` or `MAF` metadata from another aligned table."""
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
    query_annot_bed_paths: str | PathLike[str] | Sequence[str | PathLike[str]],
    baseline_annot_paths: str | PathLike[str] | Sequence[str | PathLike[str]],
    output_dir: str | Path | None = None,
    batch: bool = True,
    overwrite: bool = False,
) -> AnnotationBundle:
    """
    Convenience wrapper for BED-to-annotation projection from Python code.

    Python workflows read shared identifier, genome-build, and logging
    settings from the registered ``GlobalConfig`` rather than from per-call
    keyword arguments. The returned bundle stays aligned to the baseline
    annotation rows; optional reference-panel SNP restrictions are applied only
    in the downstream LD-score workflow.

    When ``output_dir`` is supplied, fixed ``query.<chrom>.annot.gz`` outputs
    are refused before writing unless ``overwrite=True``.
    """
    return _run_bed_to_annot_with_global_config(
        query_annot_bed_paths=query_annot_bed_paths,
        baseline_annot_paths=baseline_annot_paths,
        output_dir=output_dir,
        batch=batch,
        overwrite=overwrite,
        global_config=get_global_config(),
        entrypoint="run_bed_to_annot",
    )


def _run_bed_to_annot_with_global_config(
    query_annot_bed_paths: str | PathLike[str] | Sequence[str | PathLike[str]],
    baseline_annot_paths: str | PathLike[str] | Sequence[str | PathLike[str]],
    output_dir: str | Path | None,
    *,
    batch: bool,
    overwrite: bool,
    global_config: GlobalConfig,
    entrypoint: str,
) -> AnnotationBundle:
    """Run BED-to-annotation projection with an explicit resolved GlobalConfig."""
    print_global_config_banner(entrypoint, global_config)
    build_config = AnnotationBuildConfig(query_annot_bed_paths=query_annot_bed_paths, batch_mode=batch)
    return AnnotationBuilder(global_config, build_config).project_bed_annotations(
        query_annot_bed_paths=query_annot_bed_paths,
        baseline_annot_paths=baseline_annot_paths,
        output_dir=output_dir,
        batch=batch,
        log_level=global_config.log_level,
        overwrite=overwrite,
    )


def parse_bed_to_annot_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    """Parse CLI arguments for the BED-to-annotation workflow."""
    parser = argparse.ArgumentParser(description="Convert BED annotations into LDSC .annot.gz files.", allow_abbrev=False)
    parser.add_argument("--query-annot-bed-paths", nargs="+", required=True, help="BED files, comma-separated lists, or glob patterns.")
    parser.add_argument(
        "--baseline-annot-paths",
        nargs="+",
        required=True,
        help="Baseline annotation path tokens: exact paths, globs, or explicit @ suite tokens.",
    )
    parser.add_argument("--output-dir", required=True, help="Destination directory for generated .annot.gz files.")
    parser.add_argument("--overwrite", action="store_true", default=False, help="Replace existing fixed output files.")
    parser.add_argument("--snp-identifier", default="chr_pos", help="Identifier mode used for bundle validation.")
    parser.add_argument(
        "--genome-build",
        default=None,
        choices=("auto", "hg19", "hg37", "GRCh37", "hg38", "GRCh38"),
        help=(
            "Genome build for chr_pos inputs. Required when --snp-identifier chr_pos "
            "(the default). Use 'auto' to infer hg19/hg38 and 0-based/1-based coordinates "
            "from data. Not used when --snp-identifier rsid."
        ),
    )
    parser.add_argument(
        "--no-batch",
        dest="batch",
        action="store_false",
        default=True,
        help="Compatibility flag retained for legacy scripts; current output is always combined.",
    )
    parser.add_argument("--log-level", default="INFO", choices=("DEBUG", "INFO", "WARNING", "ERROR"), help="Logging verbosity.")
    return parser.parse_args(argv)


def main_bed_to_annot(argv: Sequence[str] | None = None) -> int:
    """CLI entrypoint for BED-to-annotation projection."""
    args = parse_bed_to_annot_args(argv)
    normalized_mode = normalize_snp_identifier_mode(args.snp_identifier)
    genome_build = _resolve_annotation_cli_genome_build(args, normalized_mode)
    cli_global_config = GlobalConfig(
        snp_identifier=normalized_mode,
        genome_build=genome_build,
        log_level=args.log_level,
    )
    _run_bed_to_annot_with_global_config(
        query_annot_bed_paths=split_cli_path_tokens(args.query_annot_bed_paths),
        baseline_annot_paths=split_cli_path_tokens(args.baseline_annot_paths),
        output_dir=args.output_dir,
        batch=args.batch,
        overwrite=args.overwrite,
        global_config=cli_global_config,
        entrypoint="main_bed_to_annot",
    )
    return 0


def _resolve_annotation_cli_genome_build(args: argparse.Namespace, snp_identifier: str) -> str | None:
    if snp_identifier == "rsid":
        return normalize_genome_build(args.genome_build)
    genome_build = normalize_genome_build(args.genome_build)
    if genome_build is None:
        raise ValueError(
            "genome_build is required when snp_identifier='chr_pos'. "
            "Pass --genome-build auto, --genome-build hg19, or --genome-build hg38."
        )
    if genome_build == "auto":
        frame, sampled_path = sample_frame_from_chr_pattern(
            split_cli_path_tokens(args.baseline_annot_paths),
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
    """Configure logging for the standalone annotation compatibility CLIs."""
    logging.basicConfig(level=getattr(logging, level.upper()), format="%(levelname)s: %(message)s")


def _get_pybedtools():
    """Import ``pybedtools`` or raise a user-facing dependency error."""
    try:
        import pybedtools
    except ImportError as exc:  # pragma: no cover - dependency check
        raise SystemExit(
            "pybedtools is required for BED-based annotation projection. Install pybedtools and bedtools, then retry."
        ) from exc
    return pybedtools


def _read_text_table(path: str) -> pd.DataFrame:
    """Read a whitespace-delimited annotation table with optional gzip compression."""
    compression = "gzip" if str(path).endswith(".gz") else None
    return pd.read_csv(path, sep=r"\s+", compression=compression)


def _annotation_shard_chromosome(path: str | Path) -> str | None:
    """Extract the chromosome token from one LDSC-style shard filename."""
    match = re.match(
        r"^.+\.(?P<chrom>\d+|X|Y|M|MT)\.(?:annot|txt|tsv)(?:\.gz)?$",
        Path(path).name,
        flags=re.IGNORECASE,
    )
    if match is None:
        return None
    return normalize_chromosome(match.group("chrom"))


def _chrom_sort_key(chrom: str) -> tuple[int, str]:
    """Return the stable chromosome ordering used by annotation workflows."""
    return chrom_sort_key(chrom)


def _to_bed_chromosome(chrom: object) -> str:
    """Convert one chromosome label into UCSC-style BED notation."""
    return "chr" + normalize_chromosome(chrom)


def _split_delimited_line(line: str, delimiter: str | None) -> list[str]:
    """Split one line according to the inferred delimiter convention."""
    line = line.rstrip("\n")
    if delimiter == ",":
        return next(csv.reader([line]))
    if delimiter == "\t":
        return line.split("\t")
    return re.split(r"\s+", line.strip())


def _detect_delimiter(path: Path) -> str | None:
    """Infer the delimiter used by a BED or tabular text file."""
    if path.name.lower().endswith(".csv") or path.name.lower().endswith(".csv.gz"):
        return ","
    opener = gzip.open if path.suffix.lower() == ".gz" else open
    with opener(path, "rt") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            if "\t" in line:
                return "\t"
            if "," in line:
                return ","
            return None
    raise ValueError(f"{path} does not contain any readable data rows.")


def _iter_text_lines(path: Path) -> Iterator[str]:
    """Yield non-empty, non-comment text lines from a plain or gzipped file."""
    opener = gzip.open if path.suffix.lower() == ".gz" else open
    with opener(path, "rt") as handle:
        for line in handle:
            stripped = line.strip()
            if stripped and not stripped.startswith("#"):
                yield line.rstrip("\n")


def _infer_column_index(header: Sequence[str], spec, path: Path) -> int:
    """Infer the index of a named column family from a parsed header row."""
    column = resolve_required_column(header, spec, context=str(path))
    return list(header).index(column)


def _list_baseline_annots(baseline_dir: Path) -> list[Path]:
    """List baseline annotation template files available in one directory."""
    paths = sorted(
        path
        for path in baseline_dir.iterdir()
        if path.is_file() and (path.name.endswith(".annot") or path.name.endswith(".annot.gz"))
    )
    if not paths:
        raise ValueError(f"No .annot or .annot.gz files found in {baseline_dir}")
    return paths


def _read_baseline_annot(path: Path) -> list[_BaselineRow]:
    """Read the baseline SNP template rows from one ``.annot`` or ``.annot.gz`` file."""
    delimiter = _detect_delimiter(path)
    opener = gzip.open if path.suffix.lower() == ".gz" else open
    with opener(path, "rt") as handle:
        reader = csv.reader(handle, delimiter=delimiter) if delimiter else None
        if reader is None:
            header_line = next(handle).rstrip("\n")
            header = re.split(r"\s+", header_line.strip())
            rows_iter = handle
        else:
            header = next(reader)
            rows_iter = reader

        chr_idx = _infer_column_index(header, ANNOTATION_METADATA_SPEC_MAP["CHR"], path)
        pos_idx = _infer_column_index(header, ANNOTATION_METADATA_SPEC_MAP["POS"], path)
        snp_idx = _infer_column_index(header, ANNOTATION_METADATA_SPEC_MAP["SNP"], path)
        cm_idx = _infer_column_index(header, ANNOTATION_METADATA_SPEC_MAP["CM"], path)

        rows: list[_BaselineRow] = []
        if reader is None:
            for line in rows_iter:
                fields = re.split(r"\s+", line.strip())
                if not fields or len(fields) < 4:
                    continue
                rows.append(
                    _BaselineRow(
                        chrom=fields[chr_idx],
                        pos=int(fields[pos_idx]),
                        snp=fields[snp_idx],
                        cm=fields[cm_idx],
                    )
                )
        else:
            for fields in rows_iter:
                if not fields or len(fields) < 4:
                    continue
                rows.append(
                    _BaselineRow(
                        chrom=fields[chr_idx],
                        pos=int(fields[pos_idx]),
                        snp=fields[snp_idx],
                        cm=fields[cm_idx],
                    )
                )
    if not rows:
        raise ValueError(f"{path} does not contain any SNP rows.")
    return rows


def _write_baseline_bed(rows: Sequence[_BaselineRow], path: Path) -> Path:
    """Write baseline SNP rows as single-base BED intervals for overlap queries."""
    with path.open("w", encoding="utf-8") as handle:
        for row in rows:
            start = row.pos - 1
            if start < 0:
                raise ValueError(f"Invalid POS={row.pos} for SNP {row.snp} in baseline template.")
            handle.write(f"{_to_bed_chromosome(row.chrom)}\t{start}\t{row.pos}\t{row.snp}\n")
    return path


def _write_normalized_bed(in_path: Path, out_path: Path) -> Path:
    """Normalize one BED input into a plain tab-delimited UCSC-style BED file."""
    with (gzip.open(in_path, "rt") if in_path.suffix.lower() == ".gz" else open(in_path, "rt")) as src:
        with out_path.open("w", encoding="utf-8") as dst:
            for line in src:
                stripped = line.strip()
                if not stripped or stripped.startswith("#"):
                    continue
                fields = _split_delimited_line(line, None)
                if len(fields) < 3:
                    raise ValueError(f"BED file {in_path} has a row with fewer than three columns: {line!r}")
                chrom = _to_bed_chromosome(fields[0])
                dst.write("\t".join([chrom] + fields[1:]) + "\n")
    return out_path


def _is_bed_path(path: Path) -> bool:
    """Return ``True`` when a path looks like a BED or gzipped BED file."""
    name = path.name.lower()
    return name.endswith(".bed") or name.endswith(".bed.gz")


def _load_restrict_snp_ids(path: Path) -> set[str]:
    """Read a reference-universe SNP filter file interpreted in ``rsid`` mode."""
    snp_ids = read_global_snp_restriction(path, "rsid")
    if not snp_ids:
        raise ValueError(f"Restriction file {path} is empty.")
    return snp_ids


def _write_restrict_table_as_bed(in_path: Path, out_path: Path, genome_build: str | None = None) -> Path:
    """Convert a ``CHR``/``POS`` restriction table into BED intervals."""
    restriction = read_global_snp_restriction(
        in_path,
        "chr_pos",
        genome_build=genome_build,
        logger=LOGGER,
    )
    with out_path.open("w", encoding="utf-8") as handle:
        for value in sorted(
            restriction,
            key=lambda item: (_chrom_sort_key(item.split(":", 1)[0]), int(item.split(":", 1)[1])),
        ):
            chrom_token, pos_token = value.split(":", 1)
            chrom = _to_bed_chromosome(chrom_token)
            pos = int(pos_token)
            start = pos - 1
            if start < 0:
                raise ValueError(f"Invalid restriction POS={pos} in {in_path}")
            handle.write(f"{chrom}\t{start}\t{pos}\n")
    return out_path


def _build_restrict_resource(
    restrict_path: Path | None,
    snp_identifier: str,
    tempdir: Path,
    genome_build: str | None = None,
) -> _RestrictResource | None:
    """Normalize the optional restriction input into the resource used during projection."""
    if restrict_path is None:
        return None
    if normalize_snp_identifier_mode(snp_identifier) == "rsid":
        snp_ids = _load_restrict_snp_ids(restrict_path)
        return _RestrictResource(mode="rsid", snp_ids=frozenset(snp_ids))
    restrict_bed_path = tempdir / "restrict_snps.normalized.bed"
    if _is_bed_path(restrict_path):
        _write_normalized_bed(restrict_path, restrict_bed_path)
    else:
        _write_restrict_table_as_bed(restrict_path, restrict_bed_path, genome_build=genome_build)
    return _RestrictResource(mode="chr_pos", bed_path=restrict_bed_path)


def _validate_and_convert_intersection(results, baseline_rows: Sequence[_BaselineRow], sanity_mode: str) -> list[bool]:
    """Validate a pybedtools intersection result and convert it into a row-aligned mask."""
    mask: list[bool] = []
    for idx, feature in enumerate(results):
        if idx >= len(baseline_rows):
            raise ValueError("Intersection returned more rows than the baseline SNP template.")
        row = baseline_rows[idx]
        fields = feature.fields
        overlap_count = int(fields[-1])
        if fields[3] != row.snp:
            raise ValueError(
                f"Intersection row order mismatch at index {idx}: expected SNP {row.snp}, got {fields[3]}"
            )
        if sanity_mode == "chr_pos":
            feature_chr = _to_bed_chromosome(fields[0])
            feature_pos = int(fields[2])
            if feature_chr != _to_bed_chromosome(row.chrom) or feature_pos != row.pos:
                raise ValueError(
                    f"Intersection coordinate mismatch at index {idx}: expected ({row.chrom}, {row.pos}), got ({fields[0]}, {fields[2]})"
                )
        mask.append(overlap_count > 0)
    if len(mask) != len(baseline_rows):
        raise ValueError(f"Intersection returned {len(mask)} rows for a baseline template with {len(baseline_rows)} rows.")
    return mask


def _build_restrict_mask(baseline_rows: Sequence[_BaselineRow], baseline_bed, restrict_resource: _RestrictResource | None) -> list[bool]:
    """Build the boolean keep mask implied by the normalized restriction input."""
    if restrict_resource is None:
        return [True] * len(baseline_rows)
    if restrict_resource.mode == "rsid":
        assert restrict_resource.snp_ids is not None
        return [row.snp in restrict_resource.snp_ids for row in baseline_rows]
    assert restrict_resource.bed_path is not None
    results = baseline_bed.intersect(str(restrict_resource.bed_path), c=True, wa=True)
    return _validate_and_convert_intersection(results, baseline_rows, sanity_mode="chr_pos")


def _compute_bed_overlap_mask(baseline_rows: Sequence[_BaselineRow], baseline_bed, bed_path: Path) -> list[bool]:
    """Compute the overlap mask for one annotation BED against the baseline SNP template."""
    results = baseline_bed.intersect(str(bed_path), c=True, wa=True)
    return _validate_and_convert_intersection(results, baseline_rows, sanity_mode="rsid")


def _compute_bed_query_columns(
    metadata: pd.DataFrame,
    bed_paths: Sequence[Path],
    tempdir: Path,
) -> pd.DataFrame:
    """Project BED overlaps onto ``metadata`` without writing `.annot.gz` files."""
    pybedtools = _get_pybedtools()
    rows = [
        _BaselineRow(
            chrom=str(row.CHR),
            pos=int(row.POS),
            snp=str(row.SNP),
            cm=str(row.CM),
        )
        for row in metadata.itertuples(index=False)
    ]
    baseline_bed_path = _write_baseline_bed(rows, tempdir / "bed_query_baseline.bed")
    baseline_bed = pybedtools.BedTool(str(baseline_bed_path))
    try:
        return pd.DataFrame(
            {
                path.stem: np.asarray(_compute_bed_overlap_mask(rows, baseline_bed, path), dtype=np.float32)
                for path in bed_paths
            },
            index=metadata.index,
        )
    finally:
        pybedtools.cleanup(remove_all=True)


def _write_annot_file(out_path: Path, rows: Sequence[_BaselineRow], annotation_names: Sequence[str], masks: Sequence[Sequence[int]]) -> None:
    """Write one LDSC-compatible ``.annot.gz`` file from precomputed annotation masks."""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(out_path, "wt") as handle:
        handle.write("\t".join([*REQUIRED_ANNOT_COLUMNS, *annotation_names]) + "\n")
        for idx, row in enumerate(rows):
            annot_values = [str(int(mask[idx])) for mask in masks]
            handle.write("\t".join([row.chrom, str(row.pos), row.snp, row.cm, *annot_values]) + "\n")


def _write_bundle_query_as_annot_files(bundle: AnnotationBundle, output_dir: Path) -> list[Path]:
    """Write one `query.<chrom>.annot.gz` file per chromosome from ``bundle``."""
    output_paths: list[Path] = []
    for chrom, out_path in zip(bundle.chromosomes, _bundle_query_annot_output_paths(bundle, output_dir)):
        chrom_mask = bundle.metadata["CHR"].astype(str) == str(chrom)
        chrom_meta = bundle.metadata.loc[chrom_mask].reset_index(drop=True).copy()
        chrom_meta = chrom_meta.rename(columns={"POS": "BP"})
        chrom_query = bundle.query_annotations.loc[chrom_mask].reset_index(drop=True)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        with gzip.open(out_path, "wt") as handle:
            pd.concat([chrom_meta, chrom_query], axis=1).to_csv(handle, sep="\t", index=False)
        output_paths.append(out_path)
    return output_paths


def _bundle_query_annot_output_paths(bundle: AnnotationBundle, output_dir: Path) -> list[Path]:
    """Return the fixed query annotation output paths for one bundle."""
    return [output_dir / f"query.{chrom}.annot.gz" for chrom in bundle.chromosomes]


def _query_output_name(path: Path) -> str:
    """Derive the output filename for one projected query annotation shard."""
    chrom = _annotation_shard_chromosome(path)
    if chrom is not None:
        return f"query.{chrom}.annot.gz"

    if path.name.endswith(".annot.gz"):
        stem = path.name[: -len(".annot.gz")]
    elif path.name.endswith(".annot"):
        stem = path.name[: -len(".annot")]
    else:
        stem = path.stem
    return f"query.{stem}.annot.gz"


def _process_baseline_file(
    baseline_path: Path,
    query_annot_bed_paths: Sequence[Path],
    output_dir: Path,
    batch: bool,
    restrict_resource: _RestrictResource | None,
    tempdir: Path,
) -> None:
    """Project all requested BED annotations onto one baseline template file."""
    rows = _read_baseline_annot(baseline_path)
    baseline_bed_path = _write_baseline_bed(rows, tempdir / f"{baseline_path.name}.snps.bed")
    pybedtools = _get_pybedtools()
    baseline_bed = pybedtools.BedTool(str(baseline_bed_path))
    restrict_mask = _build_restrict_mask(rows, baseline_bed, restrict_resource)
    kept_rows = [row for row, keep in zip(rows, restrict_mask) if keep]
    if not kept_rows:
        pybedtools.cleanup(remove_all=True)
        return
    kept_bed_path = _write_baseline_bed(kept_rows, tempdir / f"{baseline_path.name}.restricted.snps.bed")
    kept_baseline_bed = pybedtools.BedTool(str(kept_bed_path))

    if batch:
        annotation_names = [path.stem for path in query_annot_bed_paths]
        masks = []
        for bed_path in query_annot_bed_paths:
            overlap_mask = _compute_bed_overlap_mask(kept_rows, kept_baseline_bed, bed_path)
            masks.append(overlap_mask)
        output_name = _query_output_name(baseline_path)
        _write_annot_file(output_dir / output_name, kept_rows, annotation_names, masks)
    else:
        for bed_path in query_annot_bed_paths:
            overlap_mask = _compute_bed_overlap_mask(kept_rows, kept_baseline_bed, bed_path)
            bed_output_dir = output_dir / bed_path.stem
            output_name = _query_output_name(baseline_path)
            _write_annot_file(bed_output_dir / output_name, kept_rows, [bed_path.stem], [overlap_mask])

    pybedtools.cleanup(remove_all=True)


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main_bed_to_annot())
