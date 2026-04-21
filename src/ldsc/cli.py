"""Unified CLI for the refactored LDSC package.

Core functionality:
    Provide one command surface for annotation building, parquet
    reference-panel generation, LD-score calculation, summary-statistics
    munging, and regression workflows.

Overview
--------
The refactored package intentionally exposes a single command, ``ldsc``, with
subcommands grouped by user task rather than by historical script name. This
module owns only argument parsing and dispatch. Scientific work remains in the
public workflow modules, including the standalone parquet reference-panel
builder added during the restructuring pass.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Sequence

from . import annotation_builder, ldscore_calculator, ref_panel_builder


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

    ref_panel_parser = subparsers.add_parser("build-ref-panel", help="Build standard parquet reference panels.")
    _copy_actions(ref_panel_parser, ref_panel_builder.build_parser())

    sumstats_munger = _load_sumstats_munger()
    munge_parser = subparsers.add_parser("munge-sumstats", help="Munge GWAS summary statistics.")
    _copy_actions(munge_parser, sumstats_munger.kernel_parser())

    regression_runner = _load_regression_runner()
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
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv:
        command = argv[0]
        subargv = argv[1:]
        if command == "annotate":
            annotate_parser = argparse.ArgumentParser(prog="ldsc annotate", description="Build SNP-level annotations.")
            _add_annotate_arguments(annotate_parser)
            return _run_annotate(annotate_parser.parse_args(subargv))
        if command == "ldscore":
            return ldscore_calculator.main(subargv)
        if command == "build-ref-panel":
            return ref_panel_builder.main(subargv)
        if command == "munge-sumstats":
            sumstats_munger = _load_sumstats_munger()
            return sumstats_munger.main(subargv)
        if command == "h2":
            regression_runner = _load_regression_runner()
            parser = argparse.ArgumentParser(prog="ldsc h2", description="Estimate heritability from munged sumstats and LD scores.")
            regression_runner.add_h2_arguments(parser)
            return regression_runner.run_h2_from_args(parser.parse_args(subargv))
        if command == "partitioned-h2":
            regression_runner = _load_regression_runner()
            parser = argparse.ArgumentParser(
                prog="ldsc partitioned-h2",
                description="Estimate partitioned heritability, optionally looping over query annotations.",
            )
            regression_runner.add_partitioned_h2_arguments(parser)
            return regression_runner.run_partitioned_h2_from_args(parser.parse_args(subargv))
        if command == "rg":
            regression_runner = _load_regression_runner()
            parser = argparse.ArgumentParser(prog="ldsc rg", description="Estimate genetic correlation.")
            regression_runner.add_rg_arguments(parser)
            return regression_runner.run_rg_from_args(parser.parse_args(subargv))
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.command == "annotate":
        return _run_annotate(args)
    if args.command == "ldscore":
        return ldscore_calculator.run_ldscore_from_args(args)
    if args.command == "build-ref-panel":
        return ref_panel_builder.run_build_ref_panel_from_args(args)
    if args.command == "munge-sumstats":
        return _load_sumstats_munger().main(_namespace_to_argv(args, exclude={"command"}))
    if args.command == "h2":
        return _load_regression_runner().run_h2_from_args(args)
    if args.command == "partitioned-h2":
        return _load_regression_runner().run_partitioned_h2_from_args(args)
    if args.command == "rg":
        return _load_regression_runner().run_rg_from_args(args)
    raise ValueError(f"Unsupported command: {args.command}")


def _run_annotate(args: argparse.Namespace):
    """Dispatch the ``annotate`` subcommand to the appropriate workflow."""
    if args.bed_files:
        return annotation_builder.main_bed_to_annot(_namespace_to_argv(args, exclude={"command", "gene_set_file", "gene_coord_file", "windowsize", "bed_file", "nomerge", "bimfile", "annot_file"}))
    return annotation_builder.main_make_annot(
        _namespace_to_argv(args, exclude={"command", "baseline_annot", "output_dir", "restrict_snps_path", "snp_identifier", "genome_build", "batch", "bed_files", "log_level"})
    )


def _add_annotate_arguments(parser: argparse.ArgumentParser) -> None:
    """Register annotation-building arguments on a subparser."""
    parser.add_argument("--bed-files", nargs="+", default=None, help="BED files, comma-separated lists, or glob patterns.")
    parser.add_argument(
        "--baseline-annot",
        nargs="+",
        default=None,
        help="Baseline annotation path tokens: exact paths, globs, or explicit @ suite tokens.",
    )
    parser.add_argument("--output-dir", default=None, help="Destination directory for generated .annot.gz files.")
    parser.add_argument("--restrict-snps-path", default=None, help="Optional global SNP restriction file.")
    parser.add_argument("--snp-identifier", default="chr_pos", help="How to interpret restriction SNP identifiers.")
    parser.add_argument(
        "--genome-build",
        default=None,
        choices=("auto", "hg19", "hg37", "GRCh37", "hg38", "GRCh38"),
        help="Genome build for chr_pos inputs. Use 'auto' to infer hg19/hg38 and 0-based/1-based coordinates.",
    )
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


def _load_regression_runner():
    """Import the regression workflow lazily for SciPy-light entry paths."""
    from . import regression_runner

    return regression_runner


def _load_sumstats_munger():
    """Import the munging workflow lazily for SciPy-light entry paths."""
    from . import sumstats_munger

    return sumstats_munger
