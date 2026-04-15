"""Unified CLI for the refactored LDSC package.

Core functionality:
    Provide one command surface for annotation building, LD-score calculation,
    summary-statistics munging, and regression workflows.

Overview
--------
The refactored package intentionally exposes a single command, ``ldsc``, with
subcommands grouped by user task rather than by historical script name. This
module owns only argument parsing and dispatch. Scientific work remains in the
public workflow modules.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Sequence

from . import annotation_builder, ldscore_calculator, regression_runner, sumstats_munger


def build_parser() -> argparse.ArgumentParser:
    """Build the top-level ``ldsc`` command parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser whose subcommands mirror the supported public workflows.
    """
    parser = argparse.ArgumentParser(prog="ldsc", description="LDSC command line interface.")
    subparsers = parser.add_subparsers(dest="command", required=True)

    annotate_parser = subparsers.add_parser("annotate", help="Build SNP-level annotations.")
    _add_annotate_arguments(annotate_parser)

    ldscore_parser = subparsers.add_parser("ldscore", help="Compute LD scores.")
    _copy_actions(ldscore_parser, ldscore_calculator.build_parser())

    munge_parser = subparsers.add_parser("munge-sumstats", help="Munge GWAS summary statistics.")
    _copy_actions(munge_parser, sumstats_munger.kernel_parser())

    h2_parser = subparsers.add_parser("h2", help="Estimate heritability from munged sumstats and LD scores.")
    regression_runner.add_h2_arguments(h2_parser)

    partitioned_parser = subparsers.add_parser(
        "partitioned-h2",
        help="Estimate partitioned heritability, optionally looping over query annotations.",
    )
    regression_runner.add_partitioned_h2_arguments(partitioned_parser)

    rg_parser = subparsers.add_parser("rg", help="Estimate genetic correlation.")
    regression_runner.add_rg_arguments(rg_parser)
    return parser


def main(argv: Sequence[str] | None = None):
    """Run the unified command-line interface.

    Parameters
    ----------
    argv : sequence of str or None, optional
        Explicit argument vector. Default is ``None``, which reads arguments
        from ``sys.argv``.

    Returns
    -------
    object
        Workflow-specific result object returned by the dispatched subcommand.
    """
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.command == "annotate":
        return _run_annotate(args)
    if args.command == "ldscore":
        return ldscore_calculator.run_ldscore_from_args(args)
    if args.command == "munge-sumstats":
        return sumstats_munger.main(_namespace_to_argv(args, exclude={"command"}))
    if args.command == "h2":
        return regression_runner.run_h2_from_args(args)
    if args.command == "partitioned-h2":
        return regression_runner.run_partitioned_h2_from_args(args)
    if args.command == "rg":
        return regression_runner.run_rg_from_args(args)
    raise ValueError(f"Unsupported command: {args.command}")


def _run_annotate(args: argparse.Namespace):
    """Dispatch the ``annotate`` subcommand to the appropriate workflow."""
    if args.bed_files:
        return annotation_builder.main_bed_to_annot(_namespace_to_argv(args, exclude={"command", "gene_set_file", "gene_coord_file", "windowsize", "bed_file", "nomerge", "bimfile", "annot_file"}))
    return annotation_builder.main_make_annot(
        _namespace_to_argv(args, exclude={"command", "baseline_annot_dir", "output_dir", "restrict_snps_path", "snp_identifier", "batch", "bed_files", "log_level"})
    )


def _add_annotate_arguments(parser: argparse.ArgumentParser) -> None:
    """Register annotation-building arguments on a subparser."""
    parser.add_argument("--bed-files", nargs="+", default=None, help="BED files, comma-separated lists, or glob patterns.")
    parser.add_argument("--baseline-annot-dir", default=None, help="Directory containing chromosome-specific baseline .annot files.")
    parser.add_argument("--output-dir", default=None, help="Destination directory for generated .annot.gz files.")
    parser.add_argument("--restrict-snps-path", default=None, help="Optional global SNP restriction file.")
    parser.add_argument("--snp-identifier", default="chr_pos", help="How to interpret restriction SNP identifiers.")
    parser.add_argument("--no-batch", dest="batch", action="store_false", default=True, help="Write one output directory per BED file.")
    parser.add_argument("--log-level", default="INFO", choices=("DEBUG", "INFO", "WARNING", "ERROR"))
    parser.add_argument("--gene-set-file", default=None, help="A file of gene names, one line per gene.")
    parser.add_argument("--gene-coord-file", default="ENSG_coord.txt", help="Gene coordinate file with GENE, CHR, START, END.")
    parser.add_argument("--windowsize", type=int, default=None, help="Padding around the transcribed region.")
    parser.add_argument("--bed-file", default=None, help="A UCSC BED file with regions defining the annotation.")
    parser.add_argument("--nomerge", action="store_true", default=False, help="Do not merge BED intervals before projection.")
    parser.add_argument("--bimfile", default=None, help="PLINK BIM file used to define SNP rows.")
    parser.add_argument("--annot-file", default=None, help="Output .annot file path.")


def _copy_actions(target: argparse.ArgumentParser, source: argparse.ArgumentParser) -> None:
    """Clone option actions from ``source`` onto ``target``.

    This keeps the unified CLI aligned with the feature parsers defined in the
    workflow modules without duplicating every flag definition in two places.
    """
    for action in source._actions:
        if action.dest == "help":
            continue
        option_strings = list(action.option_strings)
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
            target.add_argument(*option_strings, action="store_true", default=action.default, help=action.help)
        elif action.__class__.__name__ == "_StoreFalseAction":
            target.add_argument(*option_strings, action="store_false", default=action.default, help=action.help)
        else:
            target.add_argument(*option_strings, **kwargs)


def _namespace_to_argv(args: argparse.Namespace, exclude: set[str] | None = None) -> list[str]:
    """Convert a parsed namespace back into a flat argument vector."""
    exclude = exclude or set()
    argv: list[str] = []
    for key, value in vars(args).items():
        if key in exclude or value is None or value is False:
            continue
        flag = "--" + key.replace("_", "-")
        if value is True:
            argv.append(flag)
            continue
        if isinstance(value, (list, tuple)):
            argv.append(flag)
            argv.extend(str(item) for item in value)
            continue
        argv.extend([flag, str(value)])
    return argv
