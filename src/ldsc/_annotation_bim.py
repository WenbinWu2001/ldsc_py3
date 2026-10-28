"""Metadata-only BIM grids for standalone annotation projection.

No genotype companions or annotation values are read or synthesized. Bounded
metadata scans reuse the annotation identity index and retained-row writer;
generated values subsequently use the ordinary packed query stores.
"""

import shutil
from pathlib import Path

import numpy as np
import pandas as pd

from ._annotation_identity import DiskIdentityIndex, IdentityDropSpool
from ._annotation_preflight import AUTOSOMES, resolve_annotation_inputs
from ._annotation_sources import (
    PreparedAnnotationSources, _Source, _metadata_chunks, _select_rows, _snp_row_counts,
)
from ._annotation_storage import AnnotationShard, FrameSpool
from ._input_preflight import InputGate
from .chromosome_inference import chrom_sort_key, normalize_chromosome
from .errors import LDSCInputError


def resolve_bim_grid_inputs(tokens):
    """Resolve BIM paths/globs or complete @ suites without BED/FAM companions."""
    paths, declared, issues = resolve_annotation_inputs(
        tokens, role="BIM grid", suffixes=(".bim", ".bim.gz"),
    )
    InputGate("BIM grid declarations", issues=issues).finish()
    return paths, declared


def iter_bim_metadata(path, *, chunk_rows=65536):
    """Yield one-based SNP metadata from headerless six-column BIM chunks.

    CM is ignored and represented by the normal annotation NA placeholder.
    A1/A2 remain available to the configured identity policy and output writer.
    """
    with pd.read_csv(path, sep=r"\s+", header=None, dtype=str, chunksize=chunk_rows) as reader:
        for frame in reader:
            if len(frame.columns) != 6:
                raise LDSCInputError(f"BIM grid '{path}' must have six headerless columns: CHR SNP CM BP A1 A2.")
            frame.columns = ["CHR", "SNP", "CM", "POS", "A1", "A2"]
            positions = pd.to_numeric(frame["POS"], errors="coerce")
            if not (np.isfinite(positions) & positions.gt(0) & positions.eq(np.floor(positions))).all():
                raise LDSCInputError(f"BIM grid '{path}' requires positive integer, one-based BP positions.")
            metadata = frame.loc[:, ["CHR", "SNP", "POS", "A1", "A2"]].copy()
            metadata["CHR"] = metadata["CHR"].map(lambda value: normalize_chromosome(value, context=str(path)))
            metadata["POS"] = positions.astype(np.int64)
            metadata["CM"] = np.nan
            yield metadata.reset_index(drop=True)


def prepare_bim_grid(workspace, paths, declared, *, mode, autosomes_only=False, chunk_rows=65536):
    """Stage a cleaned grid with no value stores, preserving surviving row order."""
    workspace.require_open()
    sources = []
    gate = InputGate("BIM grid content")
    for ordinal, path in enumerate(paths):
        def scan():
            spool = FrameSpool(workspace.path / f"bim-{ordinal}")
            chromosomes = set()
            columns = ()
            for metadata in iter_bim_metadata(path, chunk_rows=chunk_rows):
                chromosomes.update(metadata["CHR"])
                if autosomes_only and not chromosomes.issubset(AUTOSOMES):
                    raise LDSCInputError(f"BIM grid '{path}' contains chromosomes outside autosomes 1 through 22.")
                if declared.get(path) and chromosomes != {declared[path]}:
                    raise LDSCInputError(f"BIM grid '{path}' does not match its declared @ chromosome {declared[path]}.")
                columns = tuple(metadata.columns)
                spool.append(metadata)
            if not spool.n_rows:
                raise LDSCInputError(f"BIM grid '{path}' has no SNP rows.")
            return _Source(Path(path), spool, (), columns, chromosomes, 6, ())
        source = gate.check("BIM grid", path, scan, chrom=declared.get(path, ""))
        if source is not None:
            sources.append(source)
    gate.finish()
    groups = [[source] for source in sources]
    drops = IdentityDropSpool(workspace.path / "identity-drops")
    input_counts = {}
    index_path = workspace.path / "identity.sqlite"
    with DiskIdentityIndex(index_path, mode) as identity:
        for source in sources:
            for metadata in _metadata_chunks(source, chunk_rows):
                identity.add(metadata)
        selections, counts = _select_rows(
            groups, [chunk_rows] * len(groups), workspace, identity, drops, mode, None,
            input_counts=input_counts,
        )
    index_path.unlink()
    for spool in [*(source.spool for source in sources), *selections]:
        if spool.path.exists():
            shutil.rmtree(spool.path)
    if not counts:
        error = LDSCInputError("annotate retained no BIM SNP rows after identity cleanup. Supply valid, unique SNP identities.")
        error.annotation_drops = drops
        raise error
    shards = {chrom: AnnotationShard(workspace.path / f"chrom-{chrom}" / "metadata.parquet", (), count)
              for chrom, count in counts.items()}
    return PreparedAnnotationSources(
        shards, (), (), drops, tuple(sorted(input_counts, key=chrom_sort_key)),
        snp_rows=_snp_row_counts(input_counts, counts),
    )


