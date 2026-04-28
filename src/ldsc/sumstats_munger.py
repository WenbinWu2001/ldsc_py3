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
from os import PathLike
from pathlib import Path
from typing import Any
import warnings

import pandas as pd

from .column_inference import (
    INTERNAL_SUMSTATS_ARTIFACT_SPEC_MAP,
    resolve_optional_column,
    resolve_required_column,
)
from .config import GlobalConfig, MungeConfig, get_global_config
from .path_resolution import ensure_output_directory, resolve_scalar_path
from ._kernel import sumstats_munger as kernel_munge


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


RawSumstatsSpec = MungeConfig


@dataclass(frozen=True)
class SumstatsTable:
    """In-memory LDSC-ready summary-statistics table.

    Parameters
    ----------
    data : pandas.DataFrame
        Munged table containing at least ``SNP``, ``Z``, and ``N``.
    has_alleles : bool
        Whether allele columns are expected to be present and validated.
    source_path : str or None
        Original data source path, if known.
    trait_name : str or None
        Optional trait label for regression summaries.
    provenance : dict, optional
        Lightweight run metadata retained for debugging and output summaries.
        Default is an empty dict.
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
        ``N``, and ``Z``.

    Raises
    ------
    ValueError
        If ``path`` does not resolve to exactly one file or if the curated
        artifact is missing required LDSC columns.
    """
    resolved = resolve_scalar_path(path, label="munged sumstats")
    df = _read_curated_sumstats_artifact(resolved)
    resolved_columns = _resolve_curated_sumstats_columns(list(df.columns), context=resolved)
    df = df.loc[:, list(resolved_columns.values())].rename(
        columns={actual: canonical for canonical, actual in resolved_columns.items()}
    )
    warnings.warn(
        "load_sumstats() cannot recover the GlobalConfig that was active when this "
        "file was originally munged. Using the current global config as a proxy. "
        "Validate manually if genome_build or snp_identifier matters here.",
        UserWarning,
        stacklevel=2,
    )
    table = SumstatsTable(
        data=df.reset_index(drop=True),
        has_alleles={"A1", "A2"}.issubset(df.columns),
        source_path=resolved,
        trait_name=trait_name or Path(resolved).name,
        config_snapshot=get_global_config(),
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
            ``munge_config.output_dir``.
        global_config : GlobalConfig or None, optional
            Reserved for future shared validation. It is currently accepted to keep
            the workflow interface consistent with the rest of the package. Default
            is ``None``.

        Returns
        -------
        SumstatsTable
            Validated, in-memory table suitable for the regression workflow.
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
        args = self._build_args(raw_source, munge_config)
        data = kernel_munge.munge_sumstats(args, p=True)
        table = SumstatsTable(
            data=data.reset_index(drop=True),
            has_alleles=(not munge_config.no_alleles),
            source_path=source_path,
            trait_name=raw_source.trait_name,
            provenance={
                "sumstats_path": source_path,
                "output_dir": str(output_dir),
                "column_hints": dict(raw_source.column_hints),
            },
            config_snapshot=config_snapshot,
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
            },
        )
        return table

    def write_output(self, sumstats: SumstatsTable, output_dir: str | PathLike[str]) -> str:
        """Write a munged table to ``<output_dir>/sumstats.sumstats.gz``."""
        output_path = ensure_output_directory(output_dir, label="output directory") / "sumstats.sumstats.gz"
        columns = [col for col in ("SNP", "N", "Z", "A1", "A2", "FRQ") if col in sumstats.data.columns]
        sumstats.data.to_csv(output_path, sep="\t", index=False, columns=columns, float_format="%.3f", compression="gzip")
        return str(output_path)

    def build_run_summary(self, _sumstats: SumstatsTable | None = None) -> MungeRunSummary:
        """Return the summary captured from the most recent call to :meth:`run`."""
        if self._last_summary is None:
            raise ValueError("No munging run has been executed yet.")
        return self._last_summary

    def _build_args(self, raw_source: MungeConfig, munge_config: MungeConfig) -> argparse.Namespace:
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
    if getattr(args, "merge_alleles_path", None):
        args.merge_alleles = resolve_scalar_path(args.merge_alleles_path, label="merge-alleles file")
    else:
        args.merge_alleles = None
    return kernel_munge.munge_sumstats(args, p=True)


def kernel_parser():
    """Expose the public munging parser for CLI action cloning."""
    return build_parser()


def build_parser() -> argparse.ArgumentParser:
    """Build the public summary-statistics munging parser."""
    public = argparse.ArgumentParser(description=getattr(parser, "description", None), allow_abbrev=False)
    public.add_argument("--sumstats-path", required=True, help="Raw summary-statistics file path.")
    public.add_argument("--output-dir", required=True, help="Output directory for munged sumstats and logs.")
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
    """Count non-header rows in an input summary-statistics file."""
    openfunc, _compression = kernel_munge.get_compression(path)
    with openfunc(path) as handle:
        count = sum(1 for _ in handle)
    return max(count - 1, 0)


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
    for canonical in ("A1", "A2", "FRQ"):
        actual = resolve_optional_column(columns, INTERNAL_SUMSTATS_ARTIFACT_SPEC_MAP[canonical], context=context)
        if actual is not None:
            resolved[canonical] = actual
    return resolved


_COLUMN_HINT_ATTRS = {
    "snp": "snp",
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
