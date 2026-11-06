"""Read genome-build evidence shared by the build check and build inspection.

``_genome.build_policy.GenomeBuildCheck`` resolves one workflow build, and
``_genome.build_inspection`` reports the build of a saved artifact. Both read
the same raw evidence: table chunks, leading ``# ldsc:`` comments of text
tables, bounded per-chromosome coordinate samples, and the R2 matrix files of
a directory. Each consumer keeps its own decision and error handling.
"""

import bz2
import gzip
from pathlib import Path

import pandas as pd
import pyarrow.parquet as pq

from .._core.coordinates import normalize_chr_pos_frame

SAMPLE_ROWS = 5000

# R2 matrix file names. ``_refpanel.r2_format.r2_matrix_file_paths`` owns them, but
# this lower tier cannot import it; tests/refpanel/test_r2_matrix_chromosome.py
# checks that the two spellings agree.
R2_PARQUET_GLOB = "chr*_r2.parquet"


def r2_sidecar_path(path) -> Path | None:
    """Return the sidecar paired with an R2 Parquet file, or None for any other file name."""
    path = Path(path)
    if not path.name.endswith("_r2.parquet"):
        return None
    return path.with_name(path.name.removesuffix("_r2.parquet") + "_meta.tsv.gz")


def r2_members(root) -> dict[str, list[Path]]:
    """Return a directory's R2 Parquet files and those of its build children.

    Keys are ``""`` for files directly in ``root`` and ``"hg19"``/``"hg38"``
    for the build children; each value is sorted and may be empty. A child
    name selects an available branch; it does not declare the build.
    """
    root = Path(root)
    members = {"": sorted(root.glob(R2_PARQUET_GLOB))}
    members.update((build, sorted((root / build).glob(R2_PARQUET_GLOB))) for build in ("hg19", "hg38"))
    return members


def sample_chromosomes(frames, *, context) -> dict:
    """Keep at most ``SAMPLE_ROWS`` usable coordinate rows per chromosome.

    Every frame is streamed, so every chromosome is discovered. Each frame
    needs ``CHR`` and ``POS`` columns; rows with unusable coordinates are
    dropped. Returns ``{chromosome: sample}`` in order of first appearance.
    Reader and normalization errors propagate to the caller.
    """
    samples, counts = {}, {}
    for frame in frames:
        frame, _ = normalize_chr_pos_frame(
            frame.loc[:, ["CHR", "POS"]], context=str(context),
            coordinate_policy="drop", min_position=0, logger=None,
        )
        for chrom, group in frame.groupby("CHR", sort=False):
            room = SAMPLE_ROWS - counts.get(chrom, 0)
            if room > 0:
                sample = group.head(room)
                samples.setdefault(chrom, []).append(sample)
                counts[chrom] = counts.get(chrom, 0) + len(sample)
    return {chrom: pd.concat(parts, ignore_index=True) for chrom, parts in samples.items()}


def table_frames(path, columns, parquet=False):
    """Yield ``columns`` of a Parquet or user-supplied text table in chunks."""
    if parquet:
        with pq.ParquetFile(path) as reader:
            for batch in reader.iter_batches(batch_size=50000, columns=columns):
                yield batch.to_pandas()
    else:
        with pd.read_csv(path, usecols=columns, chunksize=50000, **text_options(path)[0]) as reader:
            yield from reader


def text_options(path):
    """Recognize supported delimiters and leading artifact metadata comments.

    This reads user-supplied text tables: it sniffs the delimiter, skips any
    leading comment or blank line, and reads ``bz2`` files. Package-written
    identity-line tables, such as R2 matrix sidecars, are read with
    ``open_identity_table()``.
    """
    metadata = {}
    opener = gzip.open if str(path).endswith('.gz') else bz2.open if str(path).endswith('.bz2') else open
    with opener(path, 'rt', encoding='utf-8') as stream:
        for index, line in enumerate(stream):
            stripped = line.strip()
            if stripped.startswith('# ldsc:'):
                key, _, value = stripped[len('# ldsc:'):].partition('=')
                metadata[key] = value
            if not stripped or (stripped.startswith('#') and not stripped.upper().startswith('#CHROM')):
                continue
            separator = '\t' if '\t' in line else ',' if ',' in line else r'\s+'
            return {'sep': separator, 'skiprows': index}, metadata
    return {'sep': r'\s+'}, metadata
