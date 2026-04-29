"""Workflow wrapper for LDSC summary-statistics munging.

Core functionality:
    Expose a typed, package-level interface to the legacy-compatible munging
    implementation and provide a public loader for curated ``.sumstats.gz``
    artifacts.

Overview
--------
This module converts the historical munging script behavior into explicit
Python dataclasses and a small service object. The public workflow boundary
accepts path-like inputs, normalizes them once, and then passes primitive
values into the internal kernel so numerical behavior stays aligned with
established LDSC outputs.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass, field
import json
import logging
from os import PathLike
from pathlib import Path
from typing import Any
import warnings

import pandas as pd

from .column_inference import (
    INTERNAL_SUMSTATS_ARTIFACT_SPEC_MAP,
    infer_chr_pos_columns,
    normalize_genome_build,
    normalize_snp_identifier_mode,
    resolve_optional_column,
    resolve_required_column,
)
from .config import GlobalConfig, MungeConfig, get_global_config
from .genome_build_inference import GenomeBuildEvidenceAccumulator, resolve_genome_build
from .path_resolution import ensure_output_directory, ensure_output_paths_available, resolve_scalar_path
from ._kernel import sumstats_munger as kernel_munge


LOGGER = logging.getLogger("LDSC.sumstats_munger")
parser = kernel_munge.parser
null_values = kernel_munge.null_values
default_cnames = kernel_munge.default_cnames
describe_cname = kernel_munge.describe_cname
numeric_cols = kernel_munge.numeric_cols
read_header = kernel_munge.read_header
get_cname_map = kernel_munge.get_cname_map
get_compression = kernel_munge.get_compression
clean_header = kernel_munge.clean_header
filter_pvals = kernel_munge.filter_pvals
filter_info = kernel_munge.filter_info
filter_frq = kernel_munge.filter_frq
filter_alleles = kernel_munge.filter_alleles
parse_dat = kernel_munge.parse_dat
process_n = kernel_munge.process_n
p_to_z = kernel_munge.p_to_z
check_median = kernel_munge.check_median
parse_flag_cnames = kernel_munge.parse_flag_cnames
allele_merge = kernel_munge.allele_merge
munge_sumstats = kernel_munge.munge_sumstats


@dataclass(frozen=True)
class SumstatsTable:
    """In-memory LDSC-ready summary-statistics table.

    Parameters
    ----------
    data : pandas.DataFrame
        Munged table containing at least ``SNP``, ``Z``, and ``N``. Current
        package-written artifacts also include canonical ``CHR`` and ``POS``
        columns, filled with missing values when the raw input had no
        coordinates.
    has_alleles : bool
        Whether allele columns are expected to be present and validated.
    source_path : str or None
        Original data source path, if known.
    trait_name : str or None
        Optional trait label for regression summaries.
    provenance : dict, optional
        Lightweight run metadata retained for debugging and output summaries.
        Default is an empty dict.
    config_snapshot : GlobalConfig or None, optional
        Shared configuration captured when the table was produced in-process or
        recovered from a neighboring ``sumstats.metadata.json`` sidecar. Legacy
        disk artifacts without that sidecar use ``None`` because their original
        munge-time configuration is not recoverable from the table alone.
    """
    data: pd.DataFrame
    has_alleles: bool
    source_path: str | None
    trait_name: str | None
    provenance: dict[str, Any] = field(default_factory=dict)
    config_snapshot: GlobalConfig | None = None

    def validate(self) -> None:
        """Validate the minimum LDSC-ready table contract."""
        required = {"SNP", "Z", "N"}
        missing = required - set(self.data.columns)
        if missing:
            raise ValueError(f"SumstatsTable is missing required columns: {sorted(missing)}")
        if self.has_alleles and not {"A1", "A2"}.issubset(self.data.columns):
            raise ValueError("SumstatsTable.has_alleles=True requires A1 and A2 columns.")

    def snp_identifiers(self) -> pd.Series:
        """Return the SNP identifier column as a string series."""
        return self.data["SNP"].astype(str)

    def subset_to(self, snps: set[str] | list[str]) -> "SumstatsTable":
        """Return a copy restricted to ``snps`` while preserving metadata."""
        keep = self.data["SNP"].astype(str).isin(set(snps))
        return SumstatsTable(
            data=self.data.loc[keep].reset_index(drop=True),
            has_alleles=self.has_alleles,
            source_path=self.source_path,
            trait_name=self.trait_name,
            provenance=dict(self.provenance),
            config_snapshot=self.config_snapshot,
        )

    def align_to_metadata(self, metadata: pd.DataFrame) -> "SumstatsTable":
        """Inner-join the table to ``metadata`` on ``SNP`` and preserve order."""
        merged = pd.merge(metadata.loc[:, ["SNP"]], self.data, how="inner", on="SNP", sort=False)
        return SumstatsTable(
            data=merged.reset_index(drop=True),
            has_alleles=self.has_alleles,
            source_path=self.source_path,
            trait_name=self.trait_name,
            provenance=dict(self.provenance),
            config_snapshot=self.config_snapshot,
        )

    def summary(self) -> dict[str, Any]:
        """Summarize the retained rows and provenance fields."""
        return {
            "n_rows": len(self.data),
            "has_alleles": self.has_alleles,
            "trait_name": self.trait_name,
            "source_path": self.source_path,
        }


@dataclass(frozen=True)
class MungeRunSummary:
    """Compact summary of one munging run."""
    n_input_rows: int
    n_retained_rows: int
    drop_counts: dict[str, int]
    inferred_columns: dict[str, str]
    used_n_rule: str
    output_paths: dict[str, str]


def load_sumstats(path: str | PathLike[str], trait_name: str | None = None) -> SumstatsTable:
    """Load one curated LDSC ``.sumstats`` artifact into a ``SumstatsTable``.

    Parameters
    ----------
    path : str or os.PathLike[str]
        Path token for the curated summary-statistics artifact. This may be a
        literal path or an exact-one glob pattern. Resolution happens at the
        workflow layer before the artifact is parsed.
    trait_name : str or None, optional
        Optional trait label propagated into downstream regression summaries.
        When omitted, the resolved filename is used. Default is ``None``.

    Returns
    -------
    SumstatsTable
        Validated in-memory table with canonical LDSC columns such as ``SNP``,
        ``CHR``, ``POS``, ``N``, and ``Z`` when present in the artifact. When a
        ``*.metadata.json`` sidecar is present next to the artifact, the
        returned table also recovers its munge-time ``GlobalConfig`` snapshot.

    Raises
    ------
    ValueError
        If ``path`` does not resolve to exactly one file or if the curated
        artifact is missing required LDSC columns.

    Notes
    -----
    This loader emits a warning before returning an unknown-provenance table
    only when the metadata sidecar is absent. Regression compatibility
    validation is skipped for that sumstats side unless the caller supplies a
    table produced in-process by :meth:`SumstatsMunger.run` or a disk artifact
    with a valid sidecar.
    """
    resolved = resolve_scalar_path(path, label="munged sumstats")
    df = _read_curated_sumstats_artifact(resolved)
    resolved_columns = _resolve_curated_sumstats_columns(list(df.columns), context=resolved)
    df = df.loc[:, list(resolved_columns.values())].rename(
        columns={actual: canonical for canonical, actual in resolved_columns.items()}
    )
    metadata_path = _sumstats_metadata_path(resolved)
    metadata = _read_sumstats_metadata(metadata_path)
    config_snapshot = _global_config_from_sumstats_metadata(metadata) if metadata is not None else None
    if metadata is None:
        warnings.warn(
            "load_sumstats() cannot recover the GlobalConfig that was active when this "
            "file was originally munged. Treating config provenance as unknown. "
            "Validate manually if genome_build or snp_identifier matters here.",
            UserWarning,
            stacklevel=2,
        )
    table = SumstatsTable(
        data=df.reset_index(drop=True),
        has_alleles={"A1", "A2"}.issubset(df.columns),
        source_path=resolved,
        trait_name=trait_name or Path(resolved).name,
        provenance={
            **({"metadata_path": str(metadata_path), "metadata": metadata} if metadata is not None else {}),
        },
        config_snapshot=config_snapshot,
    )
    table.validate()
    return table


class SumstatsMunger:
    """Run legacy-compatible munging through a typed workflow interface."""

    def __init__(self) -> None:
        """Initialize the workflow wrapper and clear any cached run summary."""
        self._last_summary: MungeRunSummary | None = None

    def run(
        self,
        raw_source: MungeConfig,
        munge_config: MungeConfig | None = None,
        global_config: GlobalConfig | None = None,
    ) -> SumstatsTable:
        """Munge one raw summary-statistics file into LDSC-ready form.

        Parameters
        ----------
        raw_source : MungeConfig
            Munging config with raw file path and optional column hints. When
            ``munge_config`` is omitted, this object also supplies output and
            QC settings.
        munge_config : MungeConfig
            Munging thresholds and output directory. The kernel writes fixed
            files named ``sumstats.sumstats.gz`` and ``sumstats.log`` inside
            ``munge_config.output_dir``. Existing fixed files are refused
            before the legacy kernel runs unless ``munge_config.overwrite`` is
            true.
        global_config : GlobalConfig or None, optional
            Shared configuration snapshot to attach to the returned
            ``SumstatsTable``. When omitted, the current package-global
            configuration is captured. Default is ``None``.

        Returns
        -------
        SumstatsTable
            Validated, in-memory table suitable for the regression workflow.
            The table includes canonical ``CHR`` and ``POS`` columns, preserves
            the active or inferred ``GlobalConfig`` snapshot, and writes the
            same provenance into ``sumstats.metadata.json`` so downstream
            regression can detect incompatible LD-score results after reload.
            Output paths for the corresponding disk artifacts are available
            through :meth:`build_run_summary`.
        """
        if munge_config is None:
            munge_config = raw_source
        elif isinstance(munge_config, GlobalConfig) and global_config is None:
            global_config = munge_config
            munge_config = raw_source
        if raw_source.sumstats_path is None:
            raise ValueError("MungeConfig.sumstats_path is required to run summary-statistics munging.")
        if munge_config.output_dir is None:
            raise ValueError("MungeConfig.output_dir is required to run summary-statistics munging.")

        config_snapshot = global_config or get_global_config()
        source_path = resolve_scalar_path(raw_source.sumstats_path, label="raw sumstats")
        output_dir = ensure_output_directory(munge_config.output_dir, label="output directory")
        fixed_output_stem = str(output_dir / "sumstats")
        metadata_path = fixed_output_stem + ".metadata.json"
        LOGGER.info(
            f"Munging summary statistics from '{source_path}' into '{output_dir}' "
            f"with snp_identifier='{config_snapshot.snp_identifier}' and genome_build='{config_snapshot.genome_build}'."
        )
        ensure_output_paths_available(
            [fixed_output_stem + ".sumstats.gz", fixed_output_stem + ".log", metadata_path],
            overwrite=munge_config.overwrite,
            label="munged output artifact",
        )
        args = self._build_args(raw_source, munge_config, config_snapshot)
        data = kernel_munge.munge_sumstats(args, p=True)
        coordinate_metadata = dict(getattr(args, "_coordinate_metadata", {}))
        table_config_snapshot = _effective_sumstats_config(config_snapshot, coordinate_metadata)
        _write_sumstats_metadata(
            metadata_path,
            config_snapshot=table_config_snapshot,
            coordinate_metadata=coordinate_metadata,
            source_path=source_path,
            sumstats_path=fixed_output_stem + ".sumstats.gz",
        )
        table = SumstatsTable(
            data=data.reset_index(drop=True),
            has_alleles=(not munge_config.no_alleles),
            source_path=source_path,
            trait_name=raw_source.trait_name,
            provenance={
                "sumstats_path": source_path,
                "output_dir": str(output_dir),
                "column_hints": dict(raw_source.column_hints),
                "metadata_path": metadata_path,
                "metadata": coordinate_metadata,
            },
            config_snapshot=table_config_snapshot,
        )
        table.validate()
        self._last_summary = MungeRunSummary(
            n_input_rows=_count_data_rows(source_path),
            n_retained_rows=len(table.data),
            drop_counts={},
            inferred_columns=dict(raw_source.column_hints),
            used_n_rule=_infer_used_n_rule(args),
            output_paths={
                "sumstats_gz": fixed_output_stem + ".sumstats.gz",
                "log": fixed_output_stem + ".log",
                "metadata_json": metadata_path,
            },
        )
        LOGGER.info(
            f"Munged {self._last_summary.n_input_rows} input rows to {self._last_summary.n_retained_rows} retained rows; "
            f"wrote '{self._last_summary.output_paths['sumstats_gz']}'."
        )
        return table

    def write_output(
        self,
        sumstats: SumstatsTable,
        output_dir: str | PathLike[str],
        overwrite: bool = False,
    ) -> str:
        """Write a munged table to ``<output_dir>/sumstats.sumstats.gz``.

        Existing files are refused unless ``overwrite=True``. The helper does
        not remove unrelated files from ``output_dir``.
        """
        output_path = ensure_output_directory(output_dir, label="output directory") / "sumstats.sumstats.gz"
        metadata_path = output_path.with_name("sumstats.metadata.json")
        ensure_output_paths_available([output_path, metadata_path], overwrite=overwrite, label="munged output artifact")
        data = sumstats.data.copy()
        if "CHR" not in data.columns:
            data["CHR"] = pd.NA
        if "POS" not in data.columns:
            data["POS"] = pd.NA
        columns = [col for col in ("SNP", "CHR", "POS", "A1", "A2", "Z", "N", "FRQ") if col in data.columns]
        data.to_csv(output_path, sep="\t", index=False, columns=columns, float_format="%.3f", compression="gzip")
        if sumstats.config_snapshot is not None:
            _write_sumstats_metadata(
                str(metadata_path),
                config_snapshot=sumstats.config_snapshot,
                coordinate_metadata=sumstats.provenance.get("metadata", {}),
                source_path=sumstats.source_path,
                sumstats_path=str(output_path),
            )
        return str(output_path)

    def build_run_summary(self, _sumstats: SumstatsTable | None = None) -> MungeRunSummary:
        """Return the summary captured from the most recent call to :meth:`run`."""
        if self._last_summary is None:
            raise ValueError("No munging run has been executed yet.")
        return self._last_summary

    def _build_args(
        self,
        raw_source: MungeConfig,
        munge_config: MungeConfig,
        global_config: GlobalConfig,
    ) -> argparse.Namespace:
        """Translate dataclass configuration into the legacy parser namespace."""
        args = parser.parse_args("")
        args.sumstats = resolve_scalar_path(raw_source.sumstats_path, label="raw sumstats")
        output_dir = ensure_output_directory(munge_config.output_dir, label="output directory")
        args.out = str(output_dir / "sumstats")
        args.N = munge_config.N
        args.N_cas = munge_config.N_cas
        args.N_con = munge_config.N_con
        args.info_min = munge_config.info_min
        args.maf_min = munge_config.maf_min
        args.n_min = munge_config.n_min
        args.nstudy_min = munge_config.nstudy_min
        args.chunksize = munge_config.chunk_size
        args.merge_alleles = (
            None
            if munge_config.merge_alleles_path is None
            else resolve_scalar_path(munge_config.merge_alleles_path, label="merge-alleles file")
        )
        args.signed_sumstats = munge_config.signed_sumstats_spec
        args.ignore = ",".join(munge_config.ignore_columns) if munge_config.ignore_columns else None
        args.no_alleles = munge_config.no_alleles
        args.a1_inc = munge_config.a1_inc
        args.keep_maf = munge_config.keep_maf
        args.daner = munge_config.daner
        args.daner_n = munge_config.daner_n
        args.snp_identifier = global_config.snp_identifier
        args.genome_build = global_config.genome_build
        for key, value in raw_source.column_hints.items():
            attr = _COLUMN_HINT_ATTRS.get(key)
            if attr is not None:
                setattr(args, attr, value)
        return args


def main(argv: list[str] | None = None):
    """Run the historical munging parser and kernel entrypoint."""
    args = build_parser().parse_args(argv)
    args.sumstats = resolve_scalar_path(args.sumstats_path, label="raw sumstats")
    output_dir = ensure_output_directory(args.output_dir, label="output directory")
    args.out = str(output_dir / "sumstats")
    metadata_path = args.out + ".metadata.json"
    LOGGER.info(f"Starting munge-sumstats for '{args.sumstats}' into '{output_dir}'.")
    ensure_output_paths_available(
        [args.out + ".sumstats.gz", args.out + ".log", metadata_path],
        overwrite=getattr(args, "overwrite", False),
        label="munged output artifact",
    )
    if getattr(args, "merge_alleles_path", None):
        args.merge_alleles = resolve_scalar_path(args.merge_alleles_path, label="merge-alleles file")
    else:
        args.merge_alleles = None
    config_snapshot = _resolve_main_global_config(args)
    data = kernel_munge.munge_sumstats(args, p=True)
    coordinate_metadata = dict(getattr(args, "_coordinate_metadata", {}))
    config_snapshot = _effective_sumstats_config(config_snapshot, coordinate_metadata)
    _write_sumstats_metadata(
        metadata_path,
        config_snapshot=config_snapshot,
        coordinate_metadata=coordinate_metadata,
        source_path=args.sumstats,
        sumstats_path=args.out + ".sumstats.gz",
    )
    LOGGER.info(f"Munged {len(data)} retained summary-statistics rows; wrote '{args.out}.sumstats.gz'.")
    return data


def _resolve_main_global_config(args: argparse.Namespace) -> GlobalConfig:
    mode = normalize_snp_identifier_mode(getattr(args, "snp_identifier", "chr_pos"))
    if mode == "rsid":
        config = GlobalConfig(
            snp_identifier="rsid",
            genome_build=normalize_genome_build(getattr(args, "genome_build", None)),
        )
        args.genome_build = config.genome_build
        return config
    genome_build = normalize_genome_build(getattr(args, "genome_build", None))
    if genome_build is None:
        raise ValueError(
            "genome_build is required when snp_identifier='chr_pos'. "
            "Pass --genome-build auto, --genome-build hg19, or --genome-build hg38."
        )
    if genome_build == "auto":
        sample = _sample_raw_sumstats_chr_pos(args.sumstats)
        genome_build = resolve_genome_build(
            "auto",
            "chr_pos",
            sample,
            context="munge-sumstats raw input",
            logger=LOGGER,
        )
    args.genome_build = genome_build
    return GlobalConfig(snp_identifier="chr_pos", genome_build=genome_build)


_INFERENCE_CHUNK_SIZE = 25_000
_INFERENCE_MAX_ROWS = 1_000_000


def _sample_raw_sumstats_chr_pos(path: str) -> pd.DataFrame:
    """Read CHR/POS rows adaptively until both inference thresholds are met."""
    header = kernel_munge.read_header(path)
    chr_col, pos_col = infer_chr_pos_columns(header, context=path)
    skiprows = kernel_munge.count_leading_sumstats_comment_lines(path)
    scorer = GenomeBuildEvidenceAccumulator()

    rows_read = 0
    reader = pd.read_csv(
        path,
        sep=r"\s+",
        compression="infer",
        usecols=[chr_col, pos_col],
        skiprows=skiprows,
        chunksize=_INFERENCE_CHUNK_SIZE,
    )
    for chunk in reader:
        pos_numeric = pd.to_numeric(chunk[pos_col], errors="coerce")
        valid = chunk[chr_col].notna() & pos_numeric.notna() & (pos_numeric >= 0)
        scorer.update(
            (str(chrom), int(pos))
            for chrom, pos in zip(chunk.loc[valid, chr_col], pos_numeric.loc[valid].astype("int64"))
        )
        rows_read += len(chunk)
        if scorer.is_sufficient() or rows_read >= _INFERENCE_MAX_ROWS:
            break

    if rows_read >= _INFERENCE_MAX_ROWS and not scorer.is_sufficient():
        LOGGER.info(
            f"Genome-build inference sampled {rows_read:,} rows and found "
            f"{scorer.informative_count} informative HM3 matches. "
            "If inference fails, pass --genome-build hg19 or --genome-build hg38 explicitly."
        )
    return scorer.to_frame()


def kernel_parser():
    """Expose the public munging parser for CLI action cloning."""
    return build_parser()


def build_parser() -> argparse.ArgumentParser:
    """Build the public summary-statistics munging parser."""
    public = argparse.ArgumentParser(description=getattr(parser, "description", None), allow_abbrev=False)
    public.add_argument("--sumstats-path", required=True, help="Raw summary-statistics file path.")
    public.add_argument("--output-dir", required=True, help="Output directory for munged sumstats and logs.")
    public.add_argument("--overwrite", action="store_true", default=False, help="Replace existing fixed output files.")
    public.add_argument("--merge-alleles-path", default=None, help="Optional merge-alleles file path.")
    for action in parser._actions:
        if action.dest in {"help", "sumstats", "out", "merge_alleles"}:
            continue
        option_strings = list(action.option_strings)
        if not option_strings:
            continue
        kwargs = {
            "dest": action.dest,
            "default": action.default,
            "required": action.required,
            "help": action.help,
        }
        if getattr(action, "choices", None) is not None:
            kwargs["choices"] = action.choices
        if getattr(action, "type", None) is not None:
            kwargs["type"] = action.type
        if getattr(action, "nargs", None) is not None:
            kwargs["nargs"] = action.nargs
        if action.const is not None:
            kwargs["const"] = action.const
        if action.__class__.__name__ == "_StoreTrueAction":
            public.add_argument(*option_strings, action="store_true", default=action.default, help=action.help)
        elif action.__class__.__name__ == "_StoreFalseAction":
            public.add_argument(*option_strings, action="store_false", default=action.default, help=action.help)
        else:
            public.add_argument(*option_strings, **kwargs)
    return public


def _count_data_rows(path: str) -> int:
    """Count raw data rows after leading ``##`` metadata lines and the header."""
    openfunc, _compression = kernel_munge.get_compression(path)
    skiprows = kernel_munge.count_leading_sumstats_comment_lines(path)
    with openfunc(path) as handle:
        count = sum(1 for _ in handle)
    return max(count - skiprows - 1, 0)


def _infer_used_n_rule(args: argparse.Namespace) -> str:
    """Summarize which sample-size rule the current munging run will apply."""
    if getattr(args, "N", None) is not None:
        return "fixed_N"
    if getattr(args, "N_cas", None) is not None and getattr(args, "N_con", None) is not None:
        return "fixed_case_control_N"
    return "input_columns"


def _read_curated_sumstats_artifact(path: str) -> pd.DataFrame:
    """Read one curated whitespace-delimited ``.sumstats`` artifact."""
    compression = "gzip" if str(path).endswith(".gz") else "infer"
    return pd.read_csv(path, sep=r"\s+", compression=compression)


def _resolve_curated_sumstats_columns(columns: list[str], *, context: str) -> dict[str, str]:
    """Resolve canonical internal sumstats artifact columns from ``columns``."""
    resolved = {
        canonical: resolve_required_column(columns, INTERNAL_SUMSTATS_ARTIFACT_SPEC_MAP[canonical], context=context)
        for canonical in ("SNP", "N", "Z")
    }
    for canonical in ("CHR", "POS", "A1", "A2", "FRQ"):
        actual = resolve_optional_column(columns, INTERNAL_SUMSTATS_ARTIFACT_SPEC_MAP[canonical], context=context)
        if actual is not None:
            resolved[canonical] = actual
    return {canonical: resolved[canonical] for canonical in ("SNP", "CHR", "POS", "N", "Z", "A1", "A2", "FRQ") if canonical in resolved}


def _sumstats_metadata_path(path: str | PathLike[str]) -> Path:
    """Return the sidecar metadata path for a curated sumstats artifact."""
    token = str(path)
    if token.endswith(".sumstats.gz"):
        return Path(token[: -len(".sumstats.gz")] + ".metadata.json")
    if token.endswith(".sumstats"):
        return Path(token[: -len(".sumstats")] + ".metadata.json")
    return Path(token).with_suffix(Path(token).suffix + ".metadata.json")


def _read_sumstats_metadata(path: Path) -> dict[str, Any] | None:
    """Read a sumstats sidecar metadata file when it exists."""
    if not path.exists():
        return None
    return json.loads(path.read_text(encoding="utf-8"))


def _global_config_from_sumstats_metadata(metadata: dict[str, Any] | None) -> GlobalConfig | None:
    """Recreate a GlobalConfig snapshot from a sumstats sidecar."""
    if not isinstance(metadata, dict):
        return None
    snapshot = metadata.get("config_snapshot")
    source = snapshot if isinstance(snapshot, dict) else metadata
    try:
        return GlobalConfig(
            snp_identifier=normalize_snp_identifier_mode(source.get("snp_identifier", "chr_pos")),
            genome_build=normalize_genome_build(source.get("genome_build")),
            log_level=source.get("log_level", "INFO"),
            fail_on_missing_metadata=bool(source.get("fail_on_missing_metadata", False)),
        )
    except Exception as exc:
        raise ValueError("Sumstats metadata sidecar has invalid GlobalConfig provenance.") from exc


def _effective_sumstats_config(config: GlobalConfig, coordinate_metadata: dict[str, Any]) -> GlobalConfig:
    """Return the config snapshot implied by coordinate metadata."""
    genome_build = coordinate_metadata.get("genome_build", config.genome_build)
    return GlobalConfig(
        snp_identifier=normalize_snp_identifier_mode(coordinate_metadata.get("snp_identifier", config.snp_identifier)),
        genome_build=normalize_genome_build(genome_build),
        log_level=config.log_level,
        fail_on_missing_metadata=config.fail_on_missing_metadata,
    )


def _write_sumstats_metadata(
    path: str | PathLike[str],
    *,
    config_snapshot: GlobalConfig,
    coordinate_metadata: dict[str, Any],
    source_path: str | None,
    sumstats_path: str,
) -> None:
    """Write sidecar metadata for a curated sumstats artifact."""
    payload = {
        "format": "ldsc.sumstats.v1",
        "snp_identifier": config_snapshot.snp_identifier,
        "genome_build": config_snapshot.genome_build,
        "source_path": source_path,
        "sumstats_path": sumstats_path,
        "coordinate_metadata": dict(coordinate_metadata),
        "config_snapshot": {
            "snp_identifier": config_snapshot.snp_identifier,
            "genome_build": config_snapshot.genome_build,
            "log_level": config_snapshot.log_level,
            "fail_on_missing_metadata": config_snapshot.fail_on_missing_metadata,
        },
    }
    if "genome_build_inferred" in coordinate_metadata:
        payload["genome_build_inferred"] = coordinate_metadata["genome_build_inferred"]
    if "coordinate_basis" in coordinate_metadata:
        payload["coordinate_basis"] = coordinate_metadata["coordinate_basis"]
    Path(path).write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


_COLUMN_HINT_ATTRS = {
    "snp": "snp",
    "chr": "chr",
    "pos": "pos",
    "N": "N_col",
    "N_col": "N_col",
    "N_cas": "N_cas_col",
    "N_cas_col": "N_cas_col",
    "N_con": "N_con_col",
    "N_con_col": "N_con_col",
    "a1": "a1",
    "a2": "a2",
    "p": "p",
    "frq": "frq",
    "signed_sumstats": "signed_sumstats",
    "info": "info",
    "info_list": "info_list",
    "nstudy": "nstudy",
}
