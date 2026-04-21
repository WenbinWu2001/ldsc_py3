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
 BED or gene-set intervals into chromosome-matched ``.annot.gz`` outputs using
baseline annotation templates. The loading path supports both single-table
inputs, where every file must share the same SNP row universe, and LDSC-style
chromosome-sharded inputs such as ``baseline.1.annot.gz`` through
``baseline.22.annot.gz``, which are validated per chromosome and then
aggregated into one bundle.

Key Functions
-------------
AnnotationBuilder.run :
    Build one aligned SNP-level annotation bundle from baseline and query
    annotation tables, with automatic aggregation of chromosome-sharded files.
run_bed_to_annot :
    Project one or more BED files onto baseline annotation templates resolved
    from strict path tokens.
make_annot_files :
    Preserve the legacy single-annotation ``make_annot.py`` behavior.

Design Notes
------------
- All annotation tables are normalized to the same required metadata columns:
  ``CHR``, ``POS``, ``SNP``, and ``CM``.
- When filenames cleanly encode one chromosome shard each, bundles are built
  independently per chromosome and concatenated in stable genomic order.
- The retained reference-panel SNP universe is applied after alignment so
  baseline, query, and metadata rows remain synchronized.
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
from ..genome_build_inference import (
    load_packaged_reference_table,
    resolve_chr_pos_table,
    validate_auto_genome_build_mode,
)
from ..path_resolution import (
    ANNOTATION_SUFFIXES,
    ensure_output_directory,
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
class AnnotationSourceSpec:
    """Describe the raw annotation input tokens needed to build a bundle.

    Annotation fields accept one token or a sequence of tokens. Each token may
    be an exact path, a standard Python glob pattern, or an explicit
    chromosome-suite placeholder using ``@``. Resolution happens inside the
    workflow layer; the stored dataclass values remain normalized strings until
    :meth:`AnnotationBuilder.run`.
    """
    baseline_annot_paths: str | PathLike[str] | tuple[str | PathLike[str], ...] | list[str | PathLike[str]] = field(default_factory=tuple)
    query_annot_paths: str | PathLike[str] | tuple[str | PathLike[str], ...] | list[str | PathLike[str]] = field(default_factory=tuple)
    bed_paths: str | PathLike[str] | tuple[str | PathLike[str], ...] | list[str | PathLike[str]] = field(default_factory=tuple)
    gene_set_paths: str | PathLike[str] | tuple[str | PathLike[str], ...] | list[str | PathLike[str]] = field(default_factory=tuple)
    allow_missing_query: bool = True

    def __post_init__(self) -> None:
        """Normalize all stored annotation path-token collections to tuples."""
        object.__setattr__(self, "baseline_annot_paths", _normalize_path_tuple(self.baseline_annot_paths))
        object.__setattr__(self, "query_annot_paths", _normalize_path_tuple(self.query_annot_paths))
        object.__setattr__(self, "bed_paths", _normalize_path_tuple(self.bed_paths))
        object.__setattr__(self, "gene_set_paths", _normalize_path_tuple(self.gene_set_paths))


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
        validate_auto_genome_build_mode(global_config.snp_identifier, global_config.genome_build)
        self.global_config = global_config
        self.build_config = build_config or AnnotationBuildConfig()

    def run(self, source_spec: AnnotationSourceSpec, chrom: str | None = None) -> AnnotationBundle:
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
        source_spec : AnnotationSourceSpec
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
        if not source_spec.baseline_annot_paths:
            raise ValueError("AnnotationSourceSpec must include at least one baseline annotation file.")
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
        source_spec: AnnotationSourceSpec,
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

        baseline_df = self._concat_or_empty(baseline_blocks, index=metadata.index)
        query_df = self._concat_or_empty(query_blocks, index=metadata.index)

        metadata, baseline_df, query_df = self._filter_aligned_tables_by_global_restriction(
            metadata, baseline_df, query_df
        )
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
        source_spec: AnnotationSourceSpec,
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
            chrom_source_spec = AnnotationSourceSpec(
                baseline_annot_paths=(baseline_by_chrom[chrom_key],),
                query_annot_paths=(() if not has_query_inputs else (query_by_chrom[chrom_key],)),
                bed_paths=source_spec.bed_paths,
                gene_set_paths=source_spec.gene_set_paths,
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
        source_spec: AnnotationSourceSpec,
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
                "bed_paths": list(source_spec.bed_paths),
                "gene_set_paths": list(source_spec.gene_set_paths),
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
        if self.global_config.snp_identifier == "chr_pos" and self.global_config.genome_build == "auto":
            metadata, _inference = resolve_chr_pos_table(
                metadata,
                context=context,
                reference_table=load_packaged_reference_table(),
                logger=LOGGER,
            )

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
        bed_files: str | PathLike[str] | Sequence[str | PathLike[str]],
        baseline_annot_paths: str | PathLike[str] | Sequence[str | PathLike[str]],
        output_dir: str | Path,
        batch: bool | None = None,
        log_level: str | None = None,
    ) -> Path:
        """
        Convert BED annotations into chromosome-matched ``.annot.gz`` outputs.

        This method uses each baseline annotation file as the SNP template,
        computes one binary mask per BED input, optionally filters baseline SNP
        rows to the configured reference-panel universe, and writes LDSC-compatible query shard
        files named ``query.<chrom>.annot.gz`` in either batch or
        one-directory-per-BED mode. ``bed_files`` and
        ``baseline_annot_paths`` both accept the public token language: exact
        paths, standard Python globs, and, for baseline annotations, explicit
        ``@`` chromosome-suite tokens. If a baseline token resolves to multiple
        files, all resolved baseline templates are processed.
        """
        batch = self.build_config.batch_mode if batch is None else batch
        _configure_logging(log_level or self.global_config.log_level)

        bed_paths = [Path(path) for path in resolve_file_group(bed_files, label="BED file")]
        stems = [path.stem for path in bed_paths]
        duplicates = sorted({stem for stem in stems if stems.count(stem) > 1})
        if duplicates:
            raise ValueError(
                "BED basenames must be unique because they become annotation names. "
                f"Duplicate names: {', '.join(duplicates)}"
            )
        baseline_paths = [
            Path(path)
            for path in resolve_file_group(
                baseline_annot_paths,
                label="baseline annotation",
                allow_chromosome_suite=True,
            )
        ]
        baseline_paths = sorted(
            [path for path in baseline_paths if path.name.endswith(".annot") or path.name.endswith(".annot.gz")],
            key=lambda path: str(path),
        )
        if not baseline_paths:
            raise ValueError("No baseline annotation files were resolved from the supplied tokens.")

        output_path = ensure_output_directory(output_dir, label="output directory")

        restrict_path = None
        if self.global_config.ref_panel_snps_path is not None:
            restrict_path = Path(
                resolve_scalar_path(
                    self.global_config.ref_panel_snps_path,
                    label="reference-panel SNP restriction",
                )
            )
            if not restrict_path.is_file():
                raise FileNotFoundError(f"Restriction file not found: {restrict_path}")

        with tempfile.TemporaryDirectory(prefix="bed2annot_") as tmpdir:
            tempdir = Path(tmpdir)
            normalized_beds = []
            for bed_path in bed_paths:
                normalized_path = tempdir / bed_path.name
                _write_normalized_bed(bed_path, normalized_path)
                normalized_beds.append(normalized_path)

            restrict_resource = _build_restrict_resource(
                restrict_path,
                self.global_config.snp_identifier,
                tempdir,
                genome_build=self.global_config.genome_build,
            )
            for baseline_path in baseline_paths:
                _process_baseline_file(
                    baseline_path=baseline_path,
                    bed_paths=normalized_beds,
                    output_dir=output_path,
                    batch=batch,
                    restrict_resource=restrict_resource,
                    tempdir=tempdir,
                )

        return output_path

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

    def _filter_aligned_tables_by_global_restriction(
        self,
        metadata: pd.DataFrame,
        baseline_df: pd.DataFrame,
        query_df: pd.DataFrame,
    ) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
        """Apply ``GlobalConfig.ref_panel_snps_path`` to aligned annotation tables."""
        restrict_path = self.global_config.ref_panel_snps_path
        if restrict_path is None:
            return metadata, baseline_df, query_df
        restriction = read_global_snp_restriction(
            restrict_path,
            self.global_config.snp_identifier,
            genome_build=self.global_config.genome_build,
            logger=LOGGER,
        )
        snp_ids = build_snp_id_series(metadata, self.global_config.snp_identifier)
        keep = snp_ids.isin(restriction)
        return (
            metadata.loc[keep].reset_index(drop=True),
            baseline_df.loc[keep].reset_index(drop=True),
            query_df.loc[keep].reset_index(drop=True),
        )


def run_bed_to_annot(
    bed_files: str | PathLike[str] | Sequence[str | PathLike[str]],
    baseline_annot_paths: str | PathLike[str] | Sequence[str | PathLike[str]],
    output_dir: str | Path,
    batch: bool = True,
) -> Path:
    """
    Convenience wrapper for BED-to-annotation projection from Python code.

    Python workflows read shared identifier, genome-build, reference-universe,
    regression-subset, and logging settings from the registered
    ``GlobalConfig`` rather than from per-call keyword arguments.
    """
    return _run_bed_to_annot_with_global_config(
        bed_files=bed_files,
        baseline_annot_paths=baseline_annot_paths,
        output_dir=output_dir,
        batch=batch,
        global_config=get_global_config(),
        entrypoint="run_bed_to_annot",
    )


def _run_bed_to_annot_with_global_config(
    bed_files: str | PathLike[str] | Sequence[str | PathLike[str]],
    baseline_annot_paths: str | PathLike[str] | Sequence[str | PathLike[str]],
    output_dir: str | Path,
    *,
    batch: bool,
    global_config: GlobalConfig,
    entrypoint: str,
) -> Path:
    """Run BED-to-annotation projection with an explicit resolved GlobalConfig."""
    print_global_config_banner(entrypoint, global_config)
    build_config = AnnotationBuildConfig(query_bed_paths=bed_files, batch_mode=batch)
    return AnnotationBuilder(global_config, build_config).project_bed_annotations(
        bed_files=bed_files,
        baseline_annot_paths=baseline_annot_paths,
        output_dir=output_dir,
        batch=batch,
        log_level=global_config.log_level,
    )


def gene_set_to_bed(args):
    """Expand a gene set plus window size into merged BED intervals."""
    pybedtools = _get_pybedtools()
    gene_set = pd.read_csv(args.gene_set_file, header=None, names=["GENE"])
    all_genes = pd.read_csv(args.gene_coord_file, delim_whitespace=True)
    df = pd.merge(gene_set, all_genes, on="GENE", how="inner")
    df["START"] = np.maximum(1, df["START"] - args.windowsize)
    df["END"] = df["END"] + args.windowsize
    intervals = [
        [_to_bed_chromosome(chrom), int(start) - 1, int(end)]
        for chrom, start, end in np.array(df[["CHR", "START", "END"]])
    ]
    return pybedtools.BedTool(intervals).sort().merge()


def make_annot_files(args, bed_for_annot):
    """Preserve the legacy single-output ``make_annot.py`` behavior."""
    builder = AnnotationBuilder(GlobalConfig(snp_identifier="rsid"))
    return builder.make_single_annotation_file(args.bimfile, args.annot_file, bed_for_annot)


def parse_bed_to_annot_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    """Parse CLI arguments for the BED-to-annotation workflow."""
    parser = argparse.ArgumentParser(description="Convert BED annotations into LDSC .annot.gz files.")
    parser.add_argument("--bed-files", nargs="+", required=True, help="BED files, comma-separated lists, or glob patterns.")
    parser.add_argument(
        "--baseline-annot",
        nargs="+",
        required=True,
        help="Baseline annotation path tokens: exact paths, globs, or explicit @ suite tokens.",
    )
    parser.add_argument("--output-dir", required=True, help="Destination directory for generated .annot.gz files.")
    parser.add_argument("--ref-panel-snps-path", default=None, help="Optional SNP restriction file defining the retained annotation/reference row universe.")
    parser.add_argument("--snp-identifier", default="chr_pos", help="How to interpret --ref-panel-snps-path when provided.")
    parser.add_argument(
        "--genome-build",
        default=None,
        choices=("auto", "hg19", "hg37", "GRCh37", "hg38", "GRCh38"),
        help="Genome build for chr_pos inputs. Use 'auto' to infer hg19/hg38 and 0-based/1-based coordinates.",
    )
    parser.add_argument("--no-batch", dest="batch", action="store_false", default=True, help="Write one output directory per BED file instead of combined output.")
    parser.add_argument("--log-level", default="INFO", choices=("DEBUG", "INFO", "WARNING", "ERROR"), help="Logging verbosity.")
    return parser.parse_args(argv)


def main_bed_to_annot(argv: Sequence[str] | None = None) -> int:
    """CLI entrypoint for BED-to-annotation projection."""
    args = parse_bed_to_annot_args(argv)
    default_config = GlobalConfig()
    cli_global_config = GlobalConfig(
        snp_identifier=normalize_snp_identifier_mode(args.snp_identifier),
        genome_build=default_config.genome_build if args.genome_build is None else args.genome_build,
        ref_panel_snps_path=None if args.ref_panel_snps_path is None else str(args.ref_panel_snps_path),
        log_level=args.log_level,
    )
    _run_bed_to_annot_with_global_config(
        bed_files=split_cli_path_tokens(args.bed_files),
        baseline_annot_paths=split_cli_path_tokens(args.baseline_annot),
        output_dir=args.output_dir,
        batch=args.batch,
        global_config=cli_global_config,
        entrypoint="main_bed_to_annot",
    )
    return 0


def parse_make_annot_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    """Parse CLI arguments for the legacy ``make_annot`` compatibility command."""
    parser = argparse.ArgumentParser()
    parser.add_argument("--gene-set-file", type=str, help="a file of gene names, one line per gene.")
    parser.add_argument(
        "--gene-coord-file",
        type=str,
        default="ENSG_coord.txt",
        help="a file with columns GENE, CHR, START, and END.",
    )
    parser.add_argument("--windowsize", type=int, help="how many base pairs to add around the transcribed region to make the annotation?")
    parser.add_argument("--bed-file", type=str, help="the UCSC bed file with the regions that make up your annotation")
    parser.add_argument(
        "--nomerge",
        action="store_true",
        default=False,
        help="don't merge the bed file before constructing the annotation.",
    )
    parser.add_argument("--bimfile", type=str, help="plink bim file for the dataset you will use to compute LD scores.")
    parser.add_argument("--annot-file", type=str, help="the name of the annot file to output.")
    return parser.parse_args(argv)


def main_make_annot(argv: Sequence[str] | None = None) -> int:
    """CLI entrypoint preserving the legacy ``make_annot.py`` interface."""
    args = parse_make_annot_args(argv)
    pybedtools = _get_pybedtools()
    if args.gene_set_file is not None:
        bed_for_annot = gene_set_to_bed(args)
    else:
        bed_for_annot = pybedtools.BedTool(args.bed_file).sort()
        if not args.nomerge:
            bed_for_annot = bed_for_annot.merge()
    make_annot_files(args, bed_for_annot)
    return 0


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


def _write_annot_file(out_path: Path, rows: Sequence[_BaselineRow], annotation_names: Sequence[str], masks: Sequence[Sequence[int]]) -> None:
    """Write one LDSC-compatible ``.annot.gz`` file from precomputed annotation masks."""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(out_path, "wt") as handle:
        handle.write("\t".join([*REQUIRED_ANNOT_COLUMNS, *annotation_names]) + "\n")
        for idx, row in enumerate(rows):
            annot_values = [str(int(mask[idx])) for mask in masks]
            handle.write("\t".join([row.chrom, str(row.pos), row.snp, row.cm, *annot_values]) + "\n")


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
    bed_paths: Sequence[Path],
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
        annotation_names = [path.stem for path in bed_paths]
        masks = []
        for bed_path in bed_paths:
            overlap_mask = _compute_bed_overlap_mask(kept_rows, kept_baseline_bed, bed_path)
            masks.append(overlap_mask)
        output_name = _query_output_name(baseline_path)
        _write_annot_file(output_dir / output_name, kept_rows, annotation_names, masks)
    else:
        for bed_path in bed_paths:
            overlap_mask = _compute_bed_overlap_mask(kept_rows, kept_baseline_bed, bed_path)
            bed_output_dir = output_dir / bed_path.stem
            output_name = _query_output_name(baseline_path)
            _write_annot_file(bed_output_dir / output_name, kept_rows, [bed_path.stem], [overlap_mask])

    pybedtools.cleanup(remove_all=True)


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main_bed_to_annot())
