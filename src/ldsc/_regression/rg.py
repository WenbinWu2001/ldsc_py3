"""Run genetic correlation over trait pairs.

``estimate_rg_pairs`` checks the traits, prevalences and anchor before any
fit, prepares each trait once and runs one joint genome-build check, then
fits every trait's own SNP heritability and every requested pair: all
unordered pairs in input order, or the anchor against every other trait.
A failed or unestimable trait or pair keeps its row with a status and an
error and does not stop the run. Inputs that share a trait label are
relabeled with the shortest source-path suffix that separates them.
"""

from __future__ import annotations

from dataclasses import replace
import json
import logging
from pathlib import Path
from typing import Sequence

import pandas as pd

from .._ldscore.results import LDScoreResult
from .._ldscore.source import LDScoreSource
from .._runtime.progress import PhaseProgress
from .._runtime.serial_calculation import serial_calculation
from .._sumstats.table import SumstatsTable
from ..config import GlobalConfig, RegressionConfig, normalize_prevalence_pair
from ..errors import LDSCInputError, LDSCUsageError, LDSCUserError
from .jackknife import JackknifeIdentifiabilityError
from .outputs import RG_CONCISE_COLUMNS, RG_FULL_COLUMNS
from .results import RgResultFamily
from .alignment import (
    _PreparedLDScores,
    _h2_dataset,
    _ldscore_identity_mode,
    _prepare_trait_sumstats,
    _regression_genome_build,
    _require_rg_alleles,
    _rg_dataset,
    _trait_label,
)
from .fitting import (
    _concise_rg_row,
    _dataset_metadata,
    _failed_h2_row,
    _failed_rg_full_row,
    _fit_h2_dataset,
    _fit_metadata,
    _fit_rg_dataset,
    _format_exception,
    _ldscore_dir_of,
    _rg_pair_metadata,
    _summarize_rg_pair,
    summarize_total_h2,
)


LOGGER = logging.getLogger("LDSC.regression_runner")


FAILED_RG_NOTE = "Failed; see the rg_full.tsv error column and diagnostics/rg.log."
UNESTIMABLE_RG_NOTE = (
    "Unestimable: a jackknife block deletion left the fit unidentifiable; "
    "see the rg_full.tsv error column and diagnostics/rg.log."
)
H2_OUT_OF_BOUNDS_RG_NOTE = (
    "Unestimable: a SNP heritability estimate is not positive (h2 out of bounds); "
    "see the rg_full.tsv error column and diagnostics/rg.log."
)


@serial_calculation
def estimate_rg_pairs(
    sumstats_tables: Sequence[SumstatsTable],
    ldscore_result: LDScoreSource | LDScoreResult,
    config: RegressionConfig,
    global_config: GlobalConfig,
    *,
    anchor_index: int | None = None,
    prevalences: Sequence[tuple[float | None, float | None]] | None = None,
) -> RgResultFamily:
    """Estimate genetic correlations for all requested trait pairs; see ``RegressionRunner.estimate_rg_pairs``.

    ``global_config`` is the SNP identifier fallback for an input without
    recorded provenance. The run fits in one controlled child.
    """
    if config.samp_prev is not None or config.pop_prev is not None:
        raise LDSCUsageError(
            f"rg regression does not read the prevalences of RegressionConfig (samp_prev={config.samp_prev!r}, "
            f"pop_prev={config.pop_prev!r}), which hold one trait's pair; ignoring them would report observed-scale "
            "estimates. Pass one (samp_prev, pop_prev) pair per trait, in input order, with prevalences=..., "
            "and leave samp_prev and pop_prev unset in the config."
        )
    tables = list(sumstats_tables)
    if len(tables) < 2:
        raise LDSCUserError(
            "rg regression requires at least two summary-statistics inputs. Most likely `--sumstats-sources` "
            "resolved to fewer than two files. Pass two or more munged `.sumstats` artifacts."
        )
    if anchor_index is not None and not 0 <= anchor_index < len(tables):
        raise LDSCUsageError(
            f"rg regression received anchor_index={anchor_index}, but valid indices are 0 through {len(tables) - 1}. "
            "Most likely the anchor trait was resolved against a different input list. "
            "Pass an anchor trait name or path that belongs to the current `--sumstats-sources` list."
        )
    if prevalences is not None and len(prevalences) != len(tables):
        raise LDSCUsageError(
            f"rg regression received {len(prevalences)} prevalence pairs for {len(tables)} traits. "
            "Most likely the prevalence list was built for a different input list. "
            "Pass one (samp_prev, pop_prev) pair per trait, in input order."
        )
    prev = ([normalize_prevalence_pair(*pair) for pair in prevalences] if prevalences is not None
            else [(None, None)] * len(tables))
    # Reject allele-free traits before any fit; a pair-level error would be
    # recorded as a failed pair instead of stopping the run.
    _require_rg_alleles(tables)
    tables = _disambiguate_trait_names(tables)
    # Run-level preparation: each trait is projected once, and one joint
    # genome-build check covers every pair, whose effective SNP identifier
    # mode stays in the LD-score directory's family.
    traits = [_prepare_trait_sumstats(table, ldscore_result) for table in tables]
    tables = [trait.table for trait in traits]
    _regression_genome_build(tables, ldscore_result,
                             _ldscore_identity_mode(ldscore_result, global_config), context="rg")
    ldscores = _PreparedLDScores(ldscore_result)

    h2_rows = []
    with PhaseProgress(LOGGER, "computation", "per-trait h2", len(tables)) as progress:
        for trait, (samp_prev, pop_prev) in zip(traits, prev):
            trait_name = _trait_label(trait.table)
            h2_dataset = None
            try:
                # As for pairs, an error fails only this trait's own fit; its pairs are still fitted.
                h2_dataset = _h2_dataset(trait, ldscores, config, global_config, build_checked=True)
                outcome = _fit_h2_dataset(h2_dataset, config, log_level=logging.DEBUG)
                row = summarize_total_h2(
                    outcome.estimator, outcome.dataset, trait_name=trait_name, n_snps_used=outcome.n_snps,
                    samp_prev=samp_prev, pop_prev=pop_prev,
                ).assign(status="ok", error="")
            except Exception as exc:
                status = "unestimable" if isinstance(exc, JackknifeIdentifiabilityError) else "failed"
                LOGGER.warning(f"SNP heritability of trait '{trait_name}' marked {status}: {_format_exception(exc)}",
                               exc_info=True)
                if status == "unestimable":
                    LOGGER.error("Trait %r jackknife diagnostics: %s", trait_name, json.dumps(exc.failures))
                row = _failed_h2_row(trait_name=trait_name, error=_format_exception(exc), status=status,
                                     n_snps=0 if h2_dataset is None else len(h2_dataset.merged),
                                     samp_prev=samp_prev, pop_prev=pop_prev)
            h2_rows.append(row)
            progress.advance(outcome=row["status"].iloc[0])
            h2_dataset = outcome = None  # release this trait before the next is built
    h2_per_trait = pd.concat(h2_rows, axis=0, ignore_index=True) if h2_rows else pd.DataFrame()

    rg_full_rows: list[dict[str, object]] = []
    rg_rows: list[dict[str, object]] = []
    per_pair_metadata: list[dict[str, object]] = []
    alignment_drops: list[tuple[int, int]] = []
    pair_kind = "anchor" if anchor_index is not None else "all_pairs"
    pair_count = len(tables) - 1 if anchor_index is not None else len(tables) * (len(tables) - 1) // 2
    with PhaseProgress(LOGGER, "computation", "rg pairs", pair_count) as progress:
        for i, j in _iter_rg_pairs(len(tables), anchor_index):
            trait_1 = _trait_label(tables[i])
            trait_2 = _trait_label(tables[j])
            drops = (0, 0)
            try:
                # Pair-specific keyed rows are built here, so an error fails only this pair.
                dataset = _rg_dataset(traits[i], traits[j], ldscores, config, global_config,
                                      build_checked=True)
                drops = (dataset.dropped_orientation_unknown_rows, dataset.dropped_incompatible_allele_rows)
                outcome = _fit_rg_dataset(dataset, config)
                fitted, dataset = outcome.estimator, outcome.dataset
                n_snps_used = outcome.n_snps
                full_row = _summarize_rg_pair(
                    fitted, dataset, trait_1=trait_1, trait_2=trait_2, pair_kind=pair_kind,
                    n_snps_used=n_snps_used,
                    samp_prev_1=prev[i][0], pop_prev_1=prev[i][1],
                    samp_prev_2=prev[j][0], pop_prev_2=prev[j][1],
                )
                note = ""
                if full_row["status"] == "unestimable":
                    note = H2_OUT_OF_BOUNDS_RG_NOTE
                    LOGGER.warning(f"rg for pair '{trait_1}' vs '{trait_2}' marked unestimable: {full_row['error']}")
                metadata = _rg_pair_metadata(tables[i], tables[j], full_row, config, pair_kind, (prev[i], prev[j]),
                                             {**_fit_metadata(outcome, "n_snps_used"), **_dataset_metadata(dataset)})
            except Exception as exc:
                drops = getattr(exc, "allele_alignment_drops", drops)
                error = _format_exception(exc)
                status = "unestimable" if isinstance(exc, JackknifeIdentifiabilityError) else "failed"
                note = UNESTIMABLE_RG_NOTE if status == "unestimable" else FAILED_RG_NOTE
                LOGGER.warning(f"rg for pair '{trait_1}' vs '{trait_2}' marked {status}: {error}", exc_info=True)
                if status == "unestimable":
                    LOGGER.error("Pair %r vs %r jackknife diagnostics: %s", trait_1, trait_2, json.dumps(exc.failures))
                fit_summary = getattr(exc, "fit_summary", None)
                full_row = _failed_rg_full_row(trait_1=trait_1, trait_2=trait_2, pair_kind=pair_kind, error=error,
                                               status=status, n_snps_used=fit_summary["n_snps"] if fit_summary else 0,
                                               samp_prev_1=prev[i][0], pop_prev_1=prev[i][1],
                                               samp_prev_2=prev[j][0], pop_prev_2=prev[j][1])
                fit = None if fit_summary is None else {
                    "n_snps_used": fit_summary["n_snps"], "n_blocks_used": fit_summary["n_blocks"],
                    **{key: value for key, value in fit_summary.items() if key not in {"n_snps", "n_blocks"}},
                }
                metadata = _rg_pair_metadata(tables[i], tables[j], full_row, config, pair_kind, (prev[i], prev[j]),
                                             fit)
            rg_full_rows.append(full_row)
            rg_rows.append(_concise_rg_row(full_row, note))
            per_pair_metadata.append(metadata)
            alignment_drops.append(drops)
            progress.advance(outcome=full_row["status"])
            dataset = outcome = fitted = None  # release this pair before the next is built

    rg_full = pd.DataFrame(rg_full_rows, columns=RG_FULL_COLUMNS)
    rg = pd.DataFrame(rg_rows, columns=RG_CONCISE_COLUMNS)
    legacy_drops = [trait.legacy_drops for trait in traits if trait.legacy_drops is not None]
    return RgResultFamily(
        rg=rg, rg_full=rg_full, h2_per_trait=h2_per_trait, per_pair_metadata=per_pair_metadata,
        sumstats_files=[table.source_path for table in tables],
        ldscore_dir=_ldscore_dir_of(ldscore_result),
        legacy_sumstats_drops=pd.concat(legacy_drops, ignore_index=True) if legacy_drops else None,
        dropped_nonfinite_statistic_rows=[trait.dropped_nonfinite_statistic_rows for trait in traits],
        dropped_orientation_unknown_rows=[unknown for unknown, _ in alignment_drops],
        dropped_incompatible_allele_rows=[incompatible for _, incompatible in alignment_drops],
    )


def _iter_rg_pairs(n_traits: int, anchor_index: int | None) -> list[tuple[int, int]]:
    """Return rg pair indices in the public output order."""
    if anchor_index is not None:
        return [(anchor_index, idx) for idx in range(n_traits) if idx != anchor_index]
    return [(i, j) for i in range(n_traits) for j in range(i + 1, n_traits)]


def _distinguishing_path_suffixes(paths: Sequence[Path]) -> list[str]:
    """Return the shortest trailing path parts that tell distinct paths apart.

    Parent directories are tried first, so files sharing a name show only their
    directories. Files in the same directory fall back to suffixes ending in the
    file name. One suffix length is used for the whole group, which keeps the
    labels independent of input order.
    """
    for parts in ([path.parent.parts for path in paths], [path.parts for path in paths]):
        for depth in range(1, max(map(len, parts)) + 1):
            suffixes = [Path(*item[-depth:]).as_posix() for item in parts]
            if len(set(suffixes)) == len(suffixes):
                return suffixes
    return [path.as_posix() for path in paths]


def _disambiguate_trait_names(tables: Sequence[SumstatsTable]) -> list[SumstatsTable]:
    """Return tables with deterministic unique trait names for rg outputs.

    A label shared by several inputs becomes ``label@<path suffix>``, using the
    shortest source-path suffix that separates that group.
    """
    labels = [_trait_label(table) for table in tables]
    unique = list(labels)
    for label in dict.fromkeys(labels):
        members = [index for index, item in enumerate(labels) if item == label]
        if len(members) == 1:
            continue
        if any(not tables[index].source_path for index in members):
            raise LDSCInputError(
                f"rg cannot give each input a unique trait label: '{label}' labels {len(members)} inputs, and "
                "at least one has no source path to tell them apart. Give each SumstatsTable a distinct "
                "trait_name."
            )
        suffixes = _distinguishing_path_suffixes([Path(tables[index].source_path) for index in members])
        for index, suffix in zip(members, suffixes):
            unique[index] = f"{label}@{suffix}"
    collisions = [label for label in dict.fromkeys(unique) if unique.count(label) > 1]
    if collisions:
        sources = ", ".join(f"'{table.source_path}'" for table, label in zip(tables, unique) if label == collisions[0])
        raise LDSCInputError(
            f"rg cannot give each input a unique trait label: '{collisions[0]}' would label {sources}. "
            "Most likely one trait name contains '@' and equals the label generated to separate "
            "inputs that share a trait name. Re-run `ldsc munge-sumstats` for one of these inputs "
            "with a distinct --trait-name."
        )

    renamed = sum(original != resolved for original, resolved in zip(labels, unique))
    if renamed:
        LOGGER.info("Renamed %d duplicate rg trait labels; resolved names and source paths are recorded in diagnostics/metadata.json.", renamed)
    for original, resolved in zip(labels, unique):
        if original != resolved:
            LOGGER.debug("Disambiguated duplicate rg trait name '%s' as '%s'.", original, resolved)
    return [replace(table, trait_name=label) for table, label in zip(tables, unique)]
