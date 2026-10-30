"""Reference-panel input checks for the shared declared-input gate.

``inspect_r2_headers`` checks R2 matrix files before any pair scan, and
``GENETIC_MAP_HEADER_CHECKS`` gives the gate its header check for
genetic-map files.
"""

from pathlib import Path

import pandas as pd

from .._core.columns import CHR_COLUMN_SPEC, POS_COLUMN_SPEC, resolve_required_column
from .._runtime.preflight import InputGate, input_issue
from .build_kernel import GENETIC_MAP_CM_SPEC


def _check_genetic_map_header(path):
    """Require chromosome, position and centimorgan columns in a genetic-map header."""
    columns = pd.read_csv(path, sep=r'\s+', nrows=0).columns
    for spec in (CHR_COLUMN_SPEC, POS_COLUMN_SPEC, GENETIC_MAP_CM_SPEC):
        resolve_required_column(columns, spec, context=str(path))


GENETIC_MAP_HEADER_CHECKS = {
    'genetic_map_hg19_sources': _check_genetic_map_header,
    'genetic_map_hg38_sources': _check_genetic_map_header,
}


def inspect_r2_headers(directory, genome_build=None, *, command, ldscore_args=None):
    """Check every R2 matrix chromosome's files, Parquet footer and sidecar header before any pair scan.

    ``command`` (``ldscore`` or ``query-r2-matrix``) names the reading command
    in each issue. With ``ldscore_args`` (parsed ldscore options), each
    chromosome's recorded build window must also admit the requested LD window,
    and its R2 matrix must store every pair (``ldsc:r2_min`` not positive).
    """
    from .r2_format import R2MatrixChromosome, r2_matrix_chromosomes, r2_matrix_file_paths, resolve_r2_build_dir

    root = Path(directory)
    children = [root / build for build in ('hg19', 'hg38') if (root / build).is_dir()]
    if genome_build == 'auto' and len(children) > 1 and not r2_matrix_chromosomes(root):
        gate = InputGate('R2 matrix declarations and headers')
        for child in children:
            gate.check('reference', child, lambda child=child: inspect_r2_headers(
                child, child.name, command=command, ldscore_args=ldscore_args))
        gate.finish()
        return root
    root = resolve_r2_build_dir(directory, genome_build)
    members = r2_matrix_chromosomes(root)
    gate = InputGate('R2 matrix declarations and headers')
    if not members:
        gate.issues.append(input_issue('reference', root, '', 'missing_required_input', 'No R2 Parquet or sidecar chromosome files found.'))
    for chrom in members:
        source = r2_matrix_file_paths(root, chrom)[0]
        chromosome = gate.check('reference', source, lambda: R2MatrixChromosome.open(root, chrom, command=command),
                                chrom=chrom)
        if chromosome is None:
            continue
        gate.check('reference', chromosome.sidecar_path, chromosome.check_sidecar_header, chrom=chrom)
        if ldscore_args is not None:
            gate.check('reference', source, lambda: chromosome.check_ld_window(ldscore_args), chrom=chrom)
            gate.check('reference', source, chromosome.reject_truncation, chrom=chrom)
    gate.finish()
    return root
