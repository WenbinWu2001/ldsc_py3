"""Validate direct LD-score input scope before projection or SNP filtering.

Path declarations select artifacts; their validated contents establish scope.
An ``@`` declaration requires all autosomes. Ordinary globs authorize their
actual matches, so a consistently missing chromosome cannot always be detected.
This workflow helper collects independent artifact failures without treating
unreadable inputs as scientific zero support.
"""

from dataclasses import dataclass
import glob
import gzip
from pathlib import Path

import pandas as pd

from .annotation_builder import AnnotationBuilder
from .chromosome_inference import normalize_chromosome
from .config import RefPanelConfig
from .errors import LDSCInputError, LDSCUserError
from .path_resolution import normalize_path_token, split_cli_path_tokens
from ._kernel import formats as parse, regions
from ._kernel.ref_panel import RefPanelLoader, _read_metadata_table, _resolve_r2_build_dir


AUTOSOMES = list(map(str, range(1, 23)))
ISSUE_COLUMNS = ["input_role", "source", "chrom", "reason", "details", "repair"]
GLOB_CAVEAT = "Globs select their actual matches. A missing file may be undetectable when the remaining required artifacts consistently cover the same subset."


@dataclass(frozen=True)
class DirectInputPreflight:
    """Validated chromosome evidence and all safely discoverable input issues."""

    scope: dict
    issues: pd.DataFrame

    @property
    def chromosomes(self):
        """Return the validated shared baseline/reference chromosome set."""
        return self.scope["chromosomes"]


def inspect_direct_inputs(args, global_config):
    """Inspect all selected annotation and reference artifacts before filtering.

    Validate BED/BIM/FAM integrity or the complete R2/sidecar binding, without
    computing correlations. Errors from independent files accumulate in one
    table. Input chromosomes are never inferred from annotation filenames.
    """
    issues = []
    declarations = []

    def issue(role, source, chrom, reason, details):
        issues.append(dict(input_role=role, source=str(source), chrom=chrom,
                           reason=reason, details=str(details),
                           repair="Supply valid, matching baseline/reference artifacts for the declared scope; inspect glob matches and @ declarations."))

    def files(tokens, role, *, plink=False):
        selected = []
        for raw in split_cli_path_tokens(tokens):
            token = normalize_path_token(raw)
            declarations.append(token)
            expanded = [(token.replace("@", chrom), chrom) for chrom in AUTOSOMES] if "@" in token else [(token, "")]
            for pattern, chrom in expanded:
                if plink:
                    if pattern.endswith((".bed", ".bim", ".fam")):
                        pattern = pattern[:-4]
                    matches = sorted({path[:-4] for suffix in (".bed", ".bim", ".fam") for path in glob.glob(pattern + suffix)})
                    if not matches and not glob.has_magic(pattern) and "@" not in token:
                        matches = sorted({path[:-4] for suffix in (".bed", ".bim", ".fam") for path in glob.glob(pattern + "*" + suffix)})
                else:
                    matches = sorted(glob.glob(pattern)) if glob.has_magic(pattern) else ([pattern] if Path(pattern).is_file() else [])
                if not matches:
                    issue(role, pattern, chrom, "missing_required_input", "No selected input exists.")
                selected.extend((path, chrom) for path in matches)
        return list(dict.fromkeys(selected))

    builder = AnnotationBuilder(global_config)
    sets = {"baseline": set(), "reference": set()}
    reference_prefixes = {}
    query_sets = []
    for role, attribute in (("baseline", "baseline_annot_sources"), ("query", "query_annot_sources")):
        for path, declared_chrom in files(getattr(args, attribute, None), role):
            try:
                metadata, _ = builder.parse_annotation_file(path)
                chroms = set(metadata.CHR.astype(str))
                if not chroms or not chroms <= set(AUTOSOMES):
                    raise LDSCInputError("Annotation contents must identify a nonempty autosomal chromosome set.")
                if declared_chrom and chroms != {declared_chrom}:
                    raise LDSCInputError(f"@ member for chromosome {declared_chrom} contains chromosomes {sorted(chroms, key=int)}.")
                if role == "baseline":
                    sets[role].update(chroms)
                else:
                    query_sets.append((path, chroms))
            except (OSError, EOFError, ValueError, LDSCUserError) as exc:
                issue(role, path, declared_chrom, "invalid_required_input", exc)

    for path, _ in files(getattr(args, "query_annot_bed_sources", None), "query"):
        try:
            with (gzip.open(path, "rt") if path.endswith(".gz") else open(path)) as stream:
                intervals = regions.parse_bed_text(stream.read(), label=path)
            query_sets.append((path, {normalize_chromosome(row.chrom) for row in intervals}))
        except (OSError, EOFError, ValueError, LDSCUserError) as exc:
            issue("query", path, "", "invalid_required_input", exc)

    r2_dir = getattr(args, "r2_dir", None)
    if r2_dir:
        try:
            root = _resolve_r2_build_dir(r2_dir, global_config.genome_build)
            members = sorted({path.name.split("_")[0][3:] for pattern in ("chr*_r2.parquet", "chr*_meta.tsv.gz") for path in root.glob(pattern)}, key=lambda value: (len(value), value))
            if not members:
                raise LDSCInputError("No R2/sidecar chromosome artifacts found.")
        except (OSError, EOFError, ValueError, LDSCUserError) as exc:
            issue("reference", r2_dir, "", "invalid_required_input", exc)
            members = []
        panel = RefPanelLoader(global_config).load(RefPanelConfig(backend="parquet_r2", r2_dir=r2_dir))
        for chrom in members:
            try:
                path = root / f"chr{chrom}_meta.tsv.gz"
                metadata = _read_metadata_table(path, None, global_config)
                if set(metadata.CHR.astype(str)) != {chrom} or chrom not in AUTOSOMES:
                    raise LDSCInputError("R2 sidecar contents disagree with their chromosome member.")
                reader = panel.build_reader(chrom, metadata=metadata)
                try:
                    for _ in reader.iter_all_pairs():
                        pass
                finally:
                    reader.close()
                sets["reference"].add(chrom)
            except (OSError, EOFError, ValueError, LDSCUserError) as exc:
                issue("reference", r2_dir, chrom, "invalid_required_input", exc)
    else:
        for prefix, declared_chrom in files(getattr(args, "plink_prefix", None), "reference", plink=True):
            member_errors = False
            for suffix in (".bed", ".bim", ".fam"):
                if not Path(prefix + suffix).is_file():
                    issue("reference", prefix + suffix, declared_chrom, "missing_required_input", "PLINK requires the complete BED/BIM/FAM trio.")
                    member_errors = True
            if member_errors:
                continue
            try:
                bim, fam = parse.PlinkBIMFile(prefix + ".bim"), parse.PlinkFAMFile(prefix + ".fam")
                chroms = set(bim.df.CHR.map(normalize_chromosome))
                if not chroms or not chroms <= set(AUTOSOMES) or len(fam.IDList) == 0:
                    raise LDSCInputError("PLINK contents must identify autosomal SNPs and at least one sample.")
                if declared_chrom and chroms != {declared_chrom}:
                    raise LDSCInputError(f"@ PLINK member for chromosome {declared_chrom} contains chromosomes {sorted(chroms, key=int)}.")
                if (pd.to_numeric(bim.df.BP, errors="raise") <= 0).any():
                    raise LDSCInputError("PLINK positions must be positive.")
                with open(prefix + ".bed", "rb") as bed:
                    if bed.read(3) != b"\x6c\x1b\x01":
                        raise LDSCInputError("Invalid PLINK BED magic number or unsupported non-SNP-major format.")
                    bed.seek(0, 2)
                    if bed.tell() != 3 + len(bim.IDList) * ((len(fam.IDList) + 3) // 4):
                        raise LDSCInputError("PLINK BED size disagrees with BIM/FAM; truncated or mismatched trio.")
                for chrom in chroms:
                    if chrom in reference_prefixes and reference_prefixes[chrom] != prefix:
                        issue("reference", prefix, chrom, "ambiguous_chromosome_input",
                              f"Multiple PLINK trios contain chromosome {chrom}: {reference_prefixes[chrom]} and {prefix}. Select one complete trio per chromosome.")
                    else:
                        reference_prefixes[chrom] = prefix
                sets["reference"].update(chroms)
            except (OSError, EOFError, ValueError, LDSCUserError) as exc:
                issue("reference", prefix, declared_chrom, "invalid_required_input", exc)

    requires_all = any("@" in token for token in declarations)
    if sets["baseline"] != sets["reference"] or (requires_all and sets["baseline"] != set(AUTOSOMES)):
        issue("alignment", "baseline/reference", "", "chromosome_set_mismatch",
              f"Validated baseline chromosomes={sorted(sets['baseline'], key=int)}; reference chromosomes={sorted(sets['reference'], key=int)}; @ requires autosomes 1–22={requires_all}.")
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


def validate_direct_scope(args, global_config, batch, output_config):
    """Write available batch diagnostics and fail once after scope inspection."""
    import logging
    from types import SimpleNamespace

    from .outputs import LDScoreDirectoryWriter
    from .query_annotations import assess_gene_coverage, gene_query_statuses, gene_control_errors

    logger = logging.getLogger("LDSC.ldscore_calculator")
    evidence = inspect_direct_inputs(args, global_config)
    errors = []
    if batch is not None:
        batch, coverage_errors = assess_gene_coverage(batch, evidence.chromosomes)
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
        diagnostic = SimpleNamespace(gene_list_batch=batch, query_statuses=gene_query_statuses(batch) if batch is not None else (),
                                     chromosome_scope=scope, input_issues=evidence.issues)
        LDScoreDirectoryWriter().write_query_diagnostics(diagnostic, output_config)
        raise LDSCInputError("LD-score input/coverage preflight failed: " + "; ".join(errors[:10]) +
                             ". Complete diagnostics: diagnostics/input_issues.tsv, diagnostics/chromosome_scope.json, and gene-list audit/summary when applicable. "
                             "Repair all required inputs or explicitly revise the submitted pathways; no pathways were truncated. "
                             "Other causes & fixes: docs/troubleshooting.md#ldscore-chromosome-coverage-preflight")
    scope = {**evidence.scope, "analysis_chromosomes": evidence.chromosomes}
    logger.info("Chromosomes entering the analysis: %s.", ", ".join(evidence.chromosomes))
    return batch, scope
