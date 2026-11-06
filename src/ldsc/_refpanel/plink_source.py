"""Select one chromosome of a PLINK reference panel, the same way in every PLINK route.

``select_plink_chromosome`` is the single BIM-level selection stage shared by
``PlinkRefPanel`` (``ldscore --plink-prefix``, ``--export-ref-metadata`` and
``build-gene-ldscore-index``) and ``R2MatrixBuilder`` (``build-r2-matrix``):
read the chromosome's ``.bim`` rows, apply ``--ref-panel-snps-file``, run SNP
identity cleanup on the rows the restriction keeps, and map
``--keep-indivs-file`` to ``.fam`` rows. The restriction applies first, so a
position is multi-allelic only when the restriction keeps more than one of
its allele pairs. Genotypes are not read here. ``build-r2-matrix`` adds
liftover and duplicate-position filtering on the returned rows, and the
runtime adapter records the reference-SNP funnel from the returned counts.

Callers read the restriction keys and the keep list once per run and pass
them to every chromosome. ``build_restriction_mask`` is the restriction match,
which the R2 matrix adapter also applies to sidecar rows, and
``build_plink_metadata_frame`` builds the retained SNP metadata with the minor
allele as ``A1`` once genotypes are read. This module is kept apart from the
GPL-derived genotype reader in ``plink.py``.
"""

from __future__ import annotations

from dataclasses import dataclass
import logging
from typing import Any, Sequence

import numpy as np
import pandas as pd

from .._core import plink_files
from .._core.chromosomes import normalize_chromosome
from .._core.coordinates import positive_int_position_series
from .._core.snp_identity import (
    RestrictionIdentityKeys,
    allele_set_series,
    base_key_series,
    clean_identity_artifact_table,
    restriction_membership_mask,
)
from ..errors import LDSCInputError, LDSCInternalError

LOGGER = logging.getLogger("LDSC.ref_panel")


@dataclass(frozen=True)
class KeepIndividuals:
    """Individual IDs read once from ``--keep-indivs-file``; matched to each ``.fam`` by IID."""

    path: str
    ids: pd.DataFrame


def read_keep_individuals(path: str | None) -> KeepIndividuals | None:
    """Read the individual keep list, or return ``None`` when no keep file is configured."""
    if not path:
        return None
    return KeepIndividuals(str(path), plink_files.FilterFile(str(path)).IDList)


@dataclass(frozen=True)
class PlinkChromosomeSource:
    """One PLINK chromosome after the restriction and SNP identity cleanup.

    ``metadata`` holds the selected ``.bim`` rows in ``.bim`` order with
    normalized ``CHR``, integer ``POS``, ``SNP``, ``CM``, ``A1``, ``A2`` and
    ``_raw_index``, the physical ``.bim``/``.bed`` row. ``identity_drops`` lists
    the rows cleanup removed. ``panel_rows`` counts the chromosome's ``.bim``
    rows, ``restriction_removed`` the rows the restriction removed (``None``
    without a restriction) and ``identity_cleanup_removed`` the rows cleanup
    removed. ``keep_indivs`` holds the retained ``.fam`` row indices, or
    ``None`` to keep every individual.
    """

    prefix: str
    bim: Any
    fam: Any
    metadata: pd.DataFrame
    identity_drops: pd.DataFrame
    keep_indivs: list[int] | None
    panel_rows: int
    restriction_removed: int | None
    identity_cleanup_removed: int

    @property
    def restricted_rows(self) -> int:
        """Rows that reached identity cleanup: the chromosome's rows the restriction kept."""
        return self.panel_rows - (self.restriction_removed or 0)


def select_plink_chromosome(
    prefix: str,
    chrom: str,
    snp_identifier: str,
    *,
    restriction: RestrictionIdentityKeys | None = None,
    keep_individuals: KeepIndividuals | None = None,
    stage: str,
) -> PlinkChromosomeSource:
    """Read and select one chromosome's PLINK SNP rows and kept individuals.

    Parameters
    ----------
    prefix : str
        PLINK prefix whose ``.bim`` contains ``chrom``; input inspection has
        already validated the trio and its positions.
    chrom : str
        Chromosome label; only the ``.bim`` rows on this chromosome are selected.
    snp_identifier : str
        SNP identifier mode used for restriction matching and identity cleanup.
    restriction : RestrictionIdentityKeys or None, optional
        Parsed ``--ref-panel-snps-file`` keys. A row without an allele-aware
        identity at a restricted rsID or position is kept for cleanup to drop
        and report.
    keep_individuals : KeepIndividuals or None, optional
        Parsed ``--keep-indivs-file`` IDs.
    stage : str
        ``stage`` label written with each identity-drop row.

    Returns
    -------
    PlinkChromosomeSource
        The selected rows, identity drops, step counts and kept individuals.
        Empty selections are returned, not raised; callers decide whether an
        emptied chromosome is skipped or stops the run.

    Raises
    ------
    ldsc.errors.LDSCInputError
        The ``.bim`` has no row on ``chrom``, or the keep list matches no
        ``.fam`` individual.
    """
    chrom = normalize_chromosome(chrom)
    bim = plink_files.PlinkBIMFile(prefix + ".bim")
    fam = plink_files.PlinkFAMFile(prefix + ".fam")
    metadata = bim.df.rename(columns={"BP": "POS"})  # a copy of the BIM table
    metadata["_raw_index"] = np.arange(len(metadata))
    metadata["CHR"] = metadata["CHR"].map({label: normalize_chromosome(label) for label in pd.unique(metadata["CHR"])})
    metadata = metadata.take(np.flatnonzero((metadata["CHR"] == chrom).to_numpy()))
    metadata.index = pd.RangeIndex(len(metadata))
    if metadata.empty:
        raise LDSCInputError(
            f"PLINK reference-panel selection found no SNPs for chromosome {chrom} in '{prefix}.bim'. "
            "Most likely the PLINK prefix does not hold that chromosome. Pass the PLINK prefix whose "
            ".bim lists the chromosome."
        )
    metadata["SNP"] = metadata["SNP"].astype(str)
    metadata["POS"] = positive_int_position_series(metadata["POS"], context=f"{prefix}.bim", label="BP")
    metadata["CM"] = pd.to_numeric(metadata["CM"], errors="coerce")
    panel_rows = len(metadata)
    restriction_removed = None
    if restriction is not None:
        metadata = metadata.loc[build_restriction_mask(metadata, restriction, snp_identifier)].reset_index(drop=True)
        restriction_removed = panel_rows - len(metadata)
    # Cleanup sees only the restricted rows. Allele-aware modes drop rows whose
    # alleles are missing, non-ACGT or identical and every row of a multi-allelic
    # position; every mode drops duplicate identities.
    cleanup = clean_identity_artifact_table(
        metadata, snp_identifier, context=f"PLINK reference panel chromosome {chrom}", stage=stage, logger=LOGGER,
    )
    return PlinkChromosomeSource(
        prefix=prefix, bim=bim, fam=fam, metadata=cleanup.cleaned, identity_drops=cleanup.dropped,
        keep_indivs=keep_individual_indices(keep_individuals, fam), panel_rows=panel_rows,
        restriction_removed=restriction_removed, identity_cleanup_removed=len(metadata) - len(cleanup.cleaned),
    )


def keep_individual_indices(keep_individuals: KeepIndividuals | None, fam) -> list[int] | None:
    """Return the ``.fam`` rows whose IID the keep list names, in ``.fam`` order; ``None`` keeps everyone."""
    if keep_individuals is None:
        return None
    indices = fam.loj(keep_individuals.ids)
    if len(indices) == 0:
        raise LDSCInputError(
            f"PLINK reference-panel selection retained no individuals after applying keep file "
            f"'{keep_individuals.path}'. Most likely the keep-list IDs do not match the `.fam` file. "
            "Check IID/FID columns in the keep file and pass a keep list from the same PLINK sample set."
        )
    return indices.tolist()


def orientation_flip_sign(freq_values: Sequence[float]) -> np.ndarray:
    """Per-SNP sign that re-orients genotypes so ``A1`` becomes the minor allele.

    ``freq_values`` is the frequency of the ``.bim`` second allele (the package's
    ``A2``). A SNP needs flipping when ``A2`` is the minor allele (``freq < 0.5``);
    at ``freq == 0.5`` PLINK order is kept (no flip), so the enforced invariant is
    ``freq(A1) <= 0.5``. Returns ``+1.0`` (keep) / ``-1.0`` (flip) as ``float32``
    so it broadcasts against ``float32`` standardized genotype columns without
    upcasting.
    """
    freq = np.asarray(freq_values, dtype=np.float64)
    return np.where(freq < 0.5, np.float32(-1.0), np.float32(1.0)).astype(np.float32)


def build_plink_metadata_frame(
    *,
    bim,
    kept_snps: Sequence[int],
    maf_values: Sequence[float],
    freq_values: Sequence[float],
) -> pd.DataFrame:
    """Build the retained SNP metadata table after PLINK filtering.

    ``A1``/``A2`` are emitted in canonical orientation: ``A1`` is the minor allele
    (``freq(A1) <= 0.5``). ``freq_values`` is the ``.bim`` ``A2`` allele frequency;
    rows with ``freq < 0.5`` have ``A1``/``A2`` swapped. ``maf_values`` is the
    folded ``min(freq, 1-freq)`` and, after the swap, equals ``freq(A1)``.
    ``build-r2-matrix`` sidecars and PLINK panel preparation for ``ldscore``,
    ``--export-ref-metadata`` and gene LD-score indexes share this table, so
    every PLINK-derived artifact gives the same ``A1`` for a SNP.
    """

    kept_snps = np.asarray(kept_snps, dtype=int)
    maf_values = np.asarray(maf_values, dtype=float)
    freq_values = np.asarray(freq_values, dtype=float)
    if not (len(kept_snps) == len(maf_values) == len(freq_values)):
        raise LDSCInternalError(
            "PLINK metadata assembly failed: kept_snps, maf_values, and "
            f"freq_values have different lengths ({len(kept_snps)} vs "
            f"{len(maf_values)} vs {len(freq_values)}). Most likely PLINK filtering "
            "bookkeeping desynchronized. Re-run with DEBUG logging and report the traceback."
        )
    kept = bim.df.iloc[kept_snps].reset_index(drop=True).copy()
    a1 = kept["A1"].astype(str).to_numpy()
    a2 = kept["A2"].astype(str).to_numpy()
    flip = orientation_flip_sign(freq_values) < 0
    canonical_a1 = np.where(flip, a2, a1)
    canonical_a2 = np.where(flip, a1, a2)
    out = pd.DataFrame(
        {
            "CHR": kept["CHR"].map(normalize_chromosome),
            "SNP": kept["SNP"].astype(str),
            "CM": pd.to_numeric(kept["CM"], errors="coerce"),
            "POS": pd.to_numeric(kept["BP"], errors="raise").astype(np.int64),
            "A1": canonical_a1,
            "A2": canonical_a2,
            "MAF": maf_values.astype(float),
        }
    )
    return out.reset_index(drop=True)


def build_restriction_mask(
    metadata: pd.DataFrame,
    restriction: RestrictionIdentityKeys,
    mode: str,
) -> np.ndarray:
    """Build the retained-SNP boolean mask for one SNP restriction set."""
    if restriction.match_kind == "identity":
        return _build_allele_aware_restriction_mask(
            metadata,
            restriction,
            mode,
            context="reference-panel SNP restriction matching",
        ).to_numpy(dtype=bool)
    return restriction_membership_mask(
        metadata,
        restriction,
        mode,
        context="reference-panel SNP restriction matching",
    ).to_numpy(dtype=bool)


def _build_allele_aware_restriction_mask(
    metadata: pd.DataFrame,
    restriction: RestrictionIdentityKeys,
    mode: str,
    *,
    context: str,
) -> pd.Series:
    """Match allele-aware restrictions without validating unrelated alleles."""
    base = base_key_series(metadata, mode, context=context)
    restricted_bases = _restriction_base_keys(restriction.keys)
    base_match = base.isin(restricted_bases)
    keep = pd.Series(False, index=metadata.index, dtype=bool)
    if not bool(base_match.any()):
        return keep

    allele_set, reasons = allele_set_series(metadata.loc[base_match], context=context)
    candidate_keys = pd.Series(pd.NA, index=allele_set.index, dtype=object)
    valid = base.loc[base_match].notna() & allele_set.notna()
    candidate_keys.loc[valid] = base.loc[base_match].loc[valid].astype(str) + ":" + allele_set.loc[valid].astype(str)
    keep.loc[base_match] = (candidate_keys.isin(restriction.keys) | reasons.notna()).to_numpy(dtype=bool)
    return keep


def _restriction_base_keys(keys: set[str]) -> set[str]:
    """Extract base keys from allele-aware identity keys."""
    base_keys: set[str] = set()
    for key in keys:
        parts = str(key).rsplit(":", 2)
        if len(parts) == 3:
            base_keys.add(parts[0])
    return base_keys
