"""Validate direct LD-score input scope before projection or SNP filtering.

Path declarations select artifacts; their validated contents establish scope.
An ``@`` declaration requires all autosomes. Ordinary globs authorize their
actual matches, so a consistently missing chromosome cannot always be detected.
This workflow helper collects independent artifact failures without treating
unreadable inputs as scientific zero support.
"""

from dataclasses import dataclass
import logging

import pandas as pd

from .._core.chromosomes import AUTOSOMES, chrom_sort_key
from ..config import RefPanelConfig
from ..errors import LDSCInputError, LDSCUserError
from .._runtime.paths import inspect_plink_inputs, split_cli_path_tokens
from .._refpanel.adapters import RefPanelLoader
from .._refpanel.r2_format import r2_matrix_chromosomes, resolve_r2_build_dir


ISSUE_COLUMNS = ["input_role", "source", "chrom", "reason", "details", "repair"]
GLOB_CAVEAT = "Globs select their actual matches. A missing file may be undetectable when the remaining required artifacts consistently cover the same subset."
SCOPE_REPAIR = (
    "'@' requires all chromosomes 1-22. For an intentional subset, use quoted '*' patterns or exact paths "
    "to select matching baseline and PLINK chromosome sets, for example 'baseline.*.annot.gz' and 'panel.*'. "
    "For a full-suite run, restore the missing files instead. Every selected query/control gene must be "
    "covered: supply the missing chromosomes or explicitly revise the gene lists. Changing '@' to '*' "
    "does not filter genes. With --r2-matrix-dir, select a directory whose chromosome set matches the baseline; "
    "directory paths do not expand patterns."
)
R2_MATRIX_REPAIR = (
    "Set --r2-matrix-dir to an R2 matrix built with this run's --snp-identifier and --genome-build, "
    "or rebuild it with build-r2-matrix; the details name the check that failed."
)
NON_AUTOSOMAL_BED_REPAIR = (
    "Remove the BED lines on chromosomes outside autosomes 1-22 and rerun; LD scores for query annotations "
    "use autosomes 1-22 only, and the remaining intervals are used as given."
)


@dataclass(frozen=True)
class DirectInputPreflight:
    """Validated chromosome evidence and all safely discoverable input issues."""

    scope: dict
    issues: pd.DataFrame

    @property
    def chromosomes(self):
        """Return the validated shared baseline/reference chromosome set."""
        return self.scope["chromosomes"]


def inspect_reference_inputs(args, global_config, *, require_autosomes=True, plink_members=None):
    """Validate references before annotation staging; retain scope and PLINK mapping.

    R2 matrix integrity requires sidecar identity checks and a streamed pair scan.
    PLINK integrity reads BIM/FAM and BED header/size, without correlations.
    The returned evidence is reused by the later content-dependent scope gate.
    """
    from .._runtime.progress import PhaseProgress
    issues, declarations = [], []
    sets = {"reference": set()}
    reference_prefixes = {}

    def issue(role, source, chrom, reason, details, repair=SCOPE_REPAIR):
        issues.append(dict(input_role=role, source=str(source), chrom=chrom,
                           reason=reason, details=str(details), repair=repair))

    r2_matrix_dir = getattr(args, "r2_matrix_dir", None)
    if r2_matrix_dir:
        try:
            root = resolve_r2_build_dir(r2_matrix_dir, global_config.genome_build)
            members = r2_matrix_chromosomes(root)
            if not members:
                raise LDSCInputError(
                    f"R2 matrix directory '{root}' contains no chromosome R2 Parquet files or SNP sidecars. "
                    "Set --r2-matrix-dir to the matching build-r2-matrix output directory and check --genome-build; "
                    "restore or rebuild the complete R2 matrix if its files are missing."
                )
        except (OSError, EOFError, ValueError, LDSCUserError) as exc:
            issue("reference", r2_matrix_dir, "", "invalid_required_input", exc, R2_MATRIX_REPAIR)
            members = []
        panel = RefPanelLoader(global_config).load(RefPanelConfig(backend="r2_matrix", r2_matrix_dir=r2_matrix_dir)) if members else None
        for chrom in members:
            try:
                if require_autosomes and chrom not in AUTOSOMES:
                    raise LDSCInputError(
                        f"R2 matrix directory '{root}' holds chromosome {chrom}, but LD scores for query "
                        "annotations use autosomes 1 through 22 only. Select an R2 matrix directory that holds "
                        "only autosomes."
                    )
                chromosome = panel.open_chromosome(chrom)
                reader = chromosome.pair_reader(chromosome.sidecar(), global_config.snp_identifier,
                                                global_config.genome_build)
                try:
                    if require_autosomes:
                        with PhaseProgress(logging.getLogger('LDSC.preflight'), 'validation', f'R2 pairs chromosome {chrom}') as progress:
                            for pairs in reader.iter_all_pairs():
                                progress.advance(object=f'chromosome {chrom} pair block')
                finally:
                    reader.close()
                sets["reference"].add(chrom)
            except (OSError, EOFError, ValueError, LDSCUserError) as exc:
                issue("reference", r2_matrix_dir, chrom, "invalid_required_input", exc, R2_MATRIX_REPAIR)
    else:
        plink_tokens = split_cli_path_tokens(getattr(args, "plink_prefix", None))
        declarations.extend(plink_tokens)
        plink = inspect_plink_inputs(plink_tokens, chromosomes=AUTOSOMES, require_complete_suite=require_autosomes, discovered_members=plink_members)
        issues.extend(plink.issues)
        reference_prefixes = plink.chromosome_prefixes
        sets["reference"].update(chrom for chrom in reference_prefixes if chrom in AUTOSOMES)
        for chrom, prefix in reference_prefixes.items():
            if require_autosomes and chrom not in AUTOSOMES:
                issue("reference", prefix, chrom, "invalid_required_input", "PLINK contents must identify autosomal SNPs.")

    return DirectInputPreflight(
        dict(chromosomes=sorted(sets['reference'], key=lambda c: (len(c), c)),
             reference_prefixes_by_chrom=reference_prefixes, declarations=declarations),
        pd.DataFrame(issues, columns=ISSUE_COLUMNS))


def inspect_direct_inputs(args, global_config, *, bed_sources, reference_inputs, annotation_sources=None, annotation_issues=()):
    """Combine prepared annotation, BED and reference evidence into one scope decision.

    ``annotation_sources`` are the prepared baseline/query annotations, or
    ``annotation_issues`` the input issues of a failed preparation.
    ``bed_sources`` are the prepared BED queries and ``reference_inputs`` the
    result of ``inspect_reference_inputs``. Errors from independent files
    accumulate in one table. Input chromosomes come from validated contents,
    never from annotation filenames.
    """
    issues = list(annotation_issues)
    declarations = []

    def issue(role, source, chrom, reason, details):
        issues.append(dict(input_role=role, source=str(source), chrom=chrom,
                           reason=reason, details=str(details),
                           repair="Supply valid, matching baseline/reference artifacts for the declared scope. " + SCOPE_REPAIR))

    def bed_query(source, line_counts):
        """Add a BED file's autosomes to the query scope; report its other lines as one issue."""
        outside = {chrom: count for chrom, count in line_counts.items() if chrom not in AUTOSOMES}
        if outside:
            labels = sorted(outside, key=chrom_sort_key)
            total = sum(outside.values())
            issues.append(dict(
                input_role="query", source=str(source), chrom=",".join(labels), reason="non_autosomal_chromosomes",
                details=(f"{total} line{'s' if total != 1 else ''} on chromosomes outside autosomes 1-22: "
                         + ", ".join(f"{chrom} ({outside[chrom]})" for chrom in labels)
                         + ". Remove these lines from the BED file and rerun"),
                repair=NON_AUTOSOMAL_BED_REPAIR,
            ))
        query_sets.append((source, {chrom for chrom in line_counts if chrom in AUTOSOMES}))

    sets = {"baseline": set(), "reference": set()}
    query_sets = []
    declarations.extend(split_cli_path_tokens(getattr(args, 'baseline_annot_sources', None)))
    declarations.extend(split_cli_path_tokens(getattr(args, 'query_annot_sources', None)))
    if annotation_sources is not None:
        sets['baseline'].update(annotation_sources.scope_chromosomes)
        query_sets.extend((path, set(chroms)) for path,role,chroms in annotation_sources.input_chromosomes if role == 'query')

    for source in bed_sources:
        if source.failure_reason:
            issue('query',source.source,'','invalid_required_input',source.details)
        else:
            bed_query(source.source, {chrom: spool.n_rows for chrom, spool in source.shards.items()})

    issues.extend(reference_inputs.issues.to_dict('records'))
    declarations.extend(reference_inputs.scope['declarations'])
    sets['reference'].update(reference_inputs.chromosomes)
    reference_prefixes = reference_inputs.scope['reference_prefixes_by_chrom']

    requires_all = any("@" in token for token in declarations)
    if sets["baseline"] != sets["reference"] or (requires_all and sets["baseline"] != set(AUTOSOMES)):
        issue("alignment", "baseline/reference", "", "chromosome_set_mismatch",
              f"Validated baseline chromosomes={sorted(sets['baseline'], key=int)}; reference chromosomes={sorted(sets['reference'], key=int)}; @ requires autosomes 1-22={requires_all}.")
    scope = sets["baseline"] & sets["reference"]
    for path, chroms in query_sets:
        missing = chroms - scope
        if missing:
            issue("query", path, ",".join(sorted(missing, key=int)), "incomplete_chromosome_coverage", "Query contents include chromosomes outside the validated baseline/reference scope; no truncation is allowed.")
    return DirectInputPreflight(
        dict(chromosomes=sorted(scope, key=int), baseline_chromosomes=sorted(sets["baseline"], key=int),
             reference_chromosomes=sorted(sets["reference"], key=int), requires_all_autosomes=requires_all,
             selection="validated_input_contents", glob_selection_caveat=GLOB_CAVEAT,
             reference_prefixes_by_chrom=reference_prefixes,
             validation_status="failed" if issues else "passed"),
        pd.DataFrame(issues, columns=ISSUE_COLUMNS),
    )


def validate_direct_scope(args, global_config, batch, output_config, **prepared_inputs):
    """Write available batch diagnostics and fail once after scope inspection."""
    import logging

    from .outputs import LDScoreDirectoryWriter, QueryDiagnostics
    from .._annotation.query_outcomes import gene_query_statuses, gene_control_errors

    logger = logging.getLogger("LDSC.ldscore_calculator")
    evidence = inspect_direct_inputs(args, global_config, **prepared_inputs)
    errors = []
    if batch is not None:
        batch, coverage_errors = batch.with_coverage(evidence.chromosomes)
        errors.extend(coverage_errors)
        errors.extend(gene_control_errors(batch))
    for row in evidence.issues.itertuples(index=False):
        errors.append(f"{row.input_role} {row.source}: {row.reason}: {row.details}")
    logger.info("Resolved input chromosomes: %s; baseline=%s; reference=%s; validation=%s.",
                evidence.chromosomes, evidence.scope["baseline_chromosomes"],
                evidence.scope["reference_chromosomes"], "failed" if errors else "passed")
    logger.info(GLOB_CAVEAT)
    if errors:
        scope = {**evidence.scope, "validation_status": "failed", "analysis_chromosomes": []}
        for error in errors:
            logger.warning("Input preflight: %s", error)
        diagnostic = QueryDiagnostics(scope, query_statuses=gene_query_statuses(batch) if batch is not None else (),
                                      gene_list_batch=batch, input_issues=evidence.issues)
        LDScoreDirectoryWriter().write_query_diagnostics(diagnostic, output_config)
        raise LDSCInputError("LD-score input/coverage preflight failed: " + "; ".join(errors[:10]) +
                             ". Complete diagnostics: diagnostics/input_issues.tsv, diagnostics/chromosome_scope.json, and gene-list audit/summary when applicable. "
                             + SCOPE_REPAIR + " No pathways were truncated. "
                             "Other causes & fixes: docs/user-guide/troubleshooting.md#ldscore-chromosome-coverage-preflight")
    scope = {**evidence.scope, "analysis_chromosomes": evidence.chromosomes}
    logger.info("Chromosomes entering the analysis: %s.", ", ".join(evidence.chromosomes))
    return batch, scope
