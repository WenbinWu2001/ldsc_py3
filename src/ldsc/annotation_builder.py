"""Public annotation preparation and standalone annotation command entry points.

``AnnotationBuilder`` prepares chromosome artifacts in an explicit output
workspace and returns one resource-owning dataset handle. It never returns or
implicitly caches a whole-genome annotation matrix. ``run_annotate`` adds the
standalone BED/gene validation and canonical incremental output contract.
"""

from __future__ import annotations

import argparse
from typing import Sequence

from ._annotation_bundle import AnnotationBundle
from ._annotation_loading import build_annotation_shards
from ._annotation_storage import AnnotationWorkspace
from ._kernel.identifiers import normalize_snp_identifier_mode
from ._kernel.snp_identity import identity_mode_family
from ._logging import materializing_overwrite_guard
from .annotate_workflow import run_annotate
from .config import AnnotationBuildConfig, GlobalConfig
from .errors import LDSCUsageError


class AnnotationBuilder:
    """Prepare chromosome-backed annotations under explicit resource ownership.

    Parameters
    ----------
    global_config : GlobalConfig
        SNP identity, genome build, and logging assumptions.
    annotation_config : AnnotationBuildConfig, optional
        Default baseline/query source configuration for ``run``.
    projection_genome_build : str, optional
        Concrete coordinate-projection build, separate from rsID identity.
    """

    def __init__(self, global_config: GlobalConfig, annotation_config: AnnotationBuildConfig | None = None,
                 *, projection_genome_build: str | None = None):
        self.global_config = global_config
        self.annotation_config = annotation_config or AnnotationBuildConfig()
        self.projection_genome_build = projection_genome_build

    def run(self, source_spec: AnnotationBuildConfig | None = None, *, output_dir) -> AnnotationBundle:
        """Prepare sources once and return an owned chromosome-shard handle.

        Parameters
        ----------
        source_spec : AnnotationBuildConfig, optional
            Baseline and optional prebuilt/BED/gene sources. Defaults to the
            constructor's configuration. Existing inputs must remain available
            and unchanged throughout their use.
        output_dir : path-like
            Required writable parent for unique private preparation files.
            Use ``run_annotate`` to publish canonical standalone query files.

        Returns
        -------
        AnnotationBundle
            Complete dataset metadata and chromosome descriptors. Use a context
            manager or call ``close()`` after the last borrower has finished.
            Read selected values with ``bundle.read(chrom, rows=..., columns=...)``.
        """
        workspace = AnnotationWorkspace(output_dir)
        try:
            return build_annotation_shards(source_spec or self.annotation_config, self.global_config, workspace,
                                           projection_genome_build=self.projection_genome_build)
        except BaseException as exc:
            from .outputs import AnnotationDirectoryWriter
            try:
                AnnotationDirectoryWriter().write_diagnostics(
                    output_dir, batch=getattr(exc, "gene_list_batch", None),
                    drops=getattr(exc, "annotation_drops", None), input_issues=getattr(exc, "input_issues", None))
            finally:
                workspace.close()
            raise


def add_annotate_arguments(parser: argparse.ArgumentParser) -> None:
    """Register ``ldsc annotate`` arguments on an existing parser.

    This is the shared argument definition used by the unified ``ldsc`` parser
    and by the standalone annotation parser built in this module.
    """
    queries = parser.add_mutually_exclusive_group(required=True)
    queries.add_argument(
        "--query-annot-bed-sources",
        nargs="+",
        help="BED files, comma-separated lists, or glob patterns.",
    )
    queries.add_argument("--query-annot-gene-list-sources", nargs="+", help="Focal one-column gene lists or globs; requires a coordinate catalog and explicit padding.")
    parser.add_argument("--gene-coordinate-file", help="One-build coordinate TSV/TSV.GZ for gene-list resolution.")
    parser.add_argument("--gene-list-resolution-policy", choices=("strict", "resolved-only"), default=None, help="Gene-only resolution policy. Default: strict.")
    parser.add_argument("--gene-exclude-regions", choices=("none", "mhc"), default=None, help="Gene-only exclusion before padding. Default: none.")
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
        "--padding-bp",
        type=int,
        default=None,
        help="Nonnegative interval padding. Required explicitly for genes; omitted BED padding is zero.",
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
            "0-based/1-based coordinates from data. Not used for rsid-family identity; gene projection records a separate catalog build."
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
    parser = argparse.ArgumentParser(description="Project BED or gene-list queries into LDSC .annot.gz files.", allow_abbrev=False)
    add_annotate_arguments(parser)
    return parser


def parse_annotate_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    """Parse CLI arguments for BED or gene-list annotation.

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


@materializing_overwrite_guard(
    lambda args: (
        (getattr(args, "output_dir"), getattr(args, "overwrite", False), "RUN_FAILED.txt")
        if getattr(args, "output_dir", None)
        else None
    ),
    command="run_annotate_from_args(...)",
)
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
        ``parse_annotate_args()``.

    Returns
    -------
    bundle : AnnotationBundle
        Produced annotation bundle.
    """
    mode = normalize_snp_identifier_mode(args.snp_identifier)
    genes = getattr(args, "query_annot_gene_list_sources", None)
    requested = getattr(args, "genome_build", None)
    if identity_mode_family(mode) == "chr_pos" and requested is None and not genes:
        raise LDSCUsageError("annotate requires --genome-build for chr_pos-family identifiers. Supply auto, hg19, or hg38.")
    config = GlobalConfig(snp_identifier=mode, genome_build=(requested or "auto") if identity_mode_family(mode) == "chr_pos" else None,
                          log_level=getattr(args, "log_level", "INFO"))
    return run_annotate(
        baseline_annot_sources=getattr(args, "baseline_annot_sources", None), output_dir=getattr(args, "output_dir", None),
        query_annot_bed_sources=getattr(args, "query_annot_bed_sources", None), query_annot_gene_list_sources=genes,
        gene_coordinate_file=getattr(args, "gene_coordinate_file", None), padding_bp=getattr(args, "padding_bp", None),
        gene_list_resolution_policy=getattr(args, "gene_list_resolution_policy", None), gene_exclude_regions=getattr(args, "gene_exclude_regions", None),
        genome_build=requested, global_config=config, overwrite=getattr(args, "overwrite", False),
    )


def main(argv: Sequence[str] | None = None) -> AnnotationBundle:
    """Parse annotation arguments, run projection, and return the bundle.

    ``main()`` replaces the removed ``main_bed_to_annot()`` public entry point.
    It is used by direct ``ldsc annotate`` dispatch and by script-style
    invocations that want parser behavior plus the produced result object.
    """
    return run_annotate_from_args(parse_annotate_args(argv))
