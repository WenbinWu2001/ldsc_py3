"""Workflow-layer wrappers for legacy summary-statistics munging."""

from __future__ import annotations

import argparse
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import pandas as pd

from .config import CommonConfig, MungeConfig
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


@dataclass(frozen=True)
class RawSumstatsSpec:
    path: str
    compression: str = "auto"
    trait_name: str | None = None
    column_hints: dict[str, str] = field(default_factory=dict)


@dataclass(frozen=True)
class SumstatsTable:
    data: pd.DataFrame
    has_alleles: bool
    source_path: str | None
    trait_name: str | None
    provenance: dict[str, Any] = field(default_factory=dict)

    def validate(self) -> None:
        required = {"SNP", "Z", "N"}
        missing = required - set(self.data.columns)
        if missing:
            raise ValueError(f"SumstatsTable is missing required columns: {sorted(missing)}")
        if self.has_alleles and not {"A1", "A2"}.issubset(self.data.columns):
            raise ValueError("SumstatsTable.has_alleles=True requires A1 and A2 columns.")

    def snp_identifiers(self) -> pd.Series:
        return self.data["SNP"].astype(str)

    def subset_to(self, snps: set[str] | list[str]) -> "SumstatsTable":
        keep = self.data["SNP"].astype(str).isin(set(snps))
        return SumstatsTable(
            data=self.data.loc[keep].reset_index(drop=True),
            has_alleles=self.has_alleles,
            source_path=self.source_path,
            trait_name=self.trait_name,
            provenance=dict(self.provenance),
        )

    def align_to_metadata(self, metadata: pd.DataFrame) -> "SumstatsTable":
        merged = pd.merge(metadata.loc[:, ["SNP"]], self.data, how="inner", on="SNP", sort=False)
        return SumstatsTable(
            data=merged.reset_index(drop=True),
            has_alleles=self.has_alleles,
            source_path=self.source_path,
            trait_name=self.trait_name,
            provenance=dict(self.provenance),
        )

    def summary(self) -> dict[str, Any]:
        return {
            "n_rows": len(self.data),
            "has_alleles": self.has_alleles,
            "trait_name": self.trait_name,
            "source_path": self.source_path,
        }


@dataclass(frozen=True)
class MungeRunSummary:
    n_input_rows: int
    n_retained_rows: int
    drop_counts: dict[str, int]
    inferred_columns: dict[str, str]
    used_n_rule: str
    output_paths: dict[str, str]


class SumstatsMunger:
    """Service wrapper around the legacy munging workflow."""

    def __init__(self) -> None:
        self._last_summary: MungeRunSummary | None = None

    def run(
        self,
        raw_source: RawSumstatsSpec,
        munge_config: MungeConfig,
        common_config: CommonConfig | None = None,
    ) -> SumstatsTable:
        del common_config  # reserved for future shared validation
        args = self._build_args(raw_source, munge_config)
        data = kernel_munge.munge_sumstats(args, p=True)
        table = SumstatsTable(
            data=data.reset_index(drop=True),
            has_alleles=(not munge_config.no_alleles),
            source_path=raw_source.path,
            trait_name=raw_source.trait_name,
            provenance={
                "sumstats_path": raw_source.path,
                "out_prefix": munge_config.out_prefix,
                "column_hints": dict(raw_source.column_hints),
            },
        )
        table.validate()
        self._last_summary = MungeRunSummary(
            n_input_rows=_count_data_rows(raw_source.path),
            n_retained_rows=len(table.data),
            drop_counts={},
            inferred_columns=dict(raw_source.column_hints),
            used_n_rule=_infer_used_n_rule(args),
            output_paths={
                "sumstats_gz": munge_config.out_prefix + ".sumstats.gz",
                "log": munge_config.out_prefix + ".log",
            },
        )
        return table

    def write_output(self, sumstats: SumstatsTable, out_prefix: str) -> str:
        output_path = out_prefix + ".sumstats.gz"
        columns = [col for col in ("SNP", "N", "Z", "A1", "A2", "FRQ") if col in sumstats.data.columns]
        sumstats.data.to_csv(output_path, sep="\t", index=False, columns=columns, float_format="%.3f", compression="gzip")
        return output_path

    def build_run_summary(self, _sumstats: SumstatsTable | None = None) -> MungeRunSummary:
        if self._last_summary is None:
            raise ValueError("No munging run has been executed yet.")
        return self._last_summary

    def _build_args(self, raw_source: RawSumstatsSpec, munge_config: MungeConfig) -> argparse.Namespace:
        args = parser.parse_args("")
        args.sumstats = raw_source.path
        args.out = munge_config.out_prefix
        args.N = munge_config.N
        args.N_cas = munge_config.N_cas
        args.N_con = munge_config.N_con
        args.info_min = munge_config.info_min
        args.maf_min = munge_config.maf_min
        args.n_min = munge_config.n_min
        args.nstudy_min = munge_config.nstudy_min
        args.chunksize = munge_config.chunk_size
        args.merge_alleles = munge_config.merge_alleles_path
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
    args = parser.parse_args(argv)
    return kernel_munge.munge_sumstats(args, p=True)


def kernel_parser():
    return parser


def _count_data_rows(path: str) -> int:
    openfunc, _compression = kernel_munge.get_compression(path)
    with openfunc(path) as handle:
        count = sum(1 for _ in handle)
    return max(count - 1, 0)


def _infer_used_n_rule(args: argparse.Namespace) -> str:
    if getattr(args, "N", None) is not None:
        return "fixed_N"
    if getattr(args, "N_cas", None) is not None and getattr(args, "N_con", None) is not None:
        return "fixed_case_control_N"
    return "input_columns"


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
