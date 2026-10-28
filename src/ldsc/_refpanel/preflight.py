"""Reference-panel input checks for the shared declared-input gate.

``inspect_r2_headers`` checks R2 panel files before any pair scan, and
``GENETIC_MAP_HEADER_CHECKS`` gives the gate its header check for
genetic-map files.
"""

from pathlib import Path

import pandas as pd

from .._core.columns import CHR_COLUMN_SPEC, POS_COLUMN_SPEC, resolve_required_column
from .._runtime.preflight import InputGate, _readable, input_issue
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


def inspect_r2_headers(directory, genome_build=None, *, ldscore_args=None):
    """Check every selected R2 Parquet/sidecar path and Parquet metadata contract before pair scans.

    With ``ldscore_args`` (parsed ldscore options), each file's recorded build
    window must also admit the requested LD window, and the file must store
    every pair (``ldsc:min_r2`` not positive).
    """
    from .r2_format import (
        read_r2_parquet_contract,
        reject_truncated_r2_panel,
        validate_ldscore_window_within_r2_panel_window,
    )
    from .adapters import _resolve_r2_build_dir

    root = Path(directory)
    children = [root / build for build in ('hg19', 'hg38') if (root / build).is_dir()]
    if genome_build == 'auto' and len(children) > 1 and not list(root.glob('chr*_r2.parquet')):
        gate = InputGate('R2 panel declarations and headers')
        for child in children:
            gate.check('reference', child,
                       lambda child=child: inspect_r2_headers(child, child.name, ldscore_args=ldscore_args))
        gate.finish()
        return root
    root = _resolve_r2_build_dir(directory, genome_build)
    members = sorted({path.name.split('_')[0][3:] for pattern in ('chr*_r2.parquet', 'chr*_meta.tsv.gz') for path in root.glob(pattern)})
    gate = InputGate('R2 panel declarations and headers')
    if not members:
        gate.issues.append(input_issue('reference', root, '', 'missing_required_input', 'No R2 Parquet or sidecar chromosome files found.'))
    for chrom in members:
        pair = root / f'chr{chrom}_r2.parquet'
        sidecar = root / f'chr{chrom}_meta.tsv.gz'
        gate.check('reference', sidecar, lambda: _readable(sidecar), chrom=chrom)
        contract = gate.check('reference', pair, lambda: read_r2_parquet_contract(pair), chrom=chrom)
        if contract is not None and ldscore_args is not None:
            gate.check('reference', pair, lambda: validate_ldscore_window_within_r2_panel_window(
                ldscore_args, parquet_paths=[str(pair)], chrom=chrom), chrom=chrom)
            gate.check('reference', pair, lambda: reject_truncated_r2_panel([str(pair)], chrom), chrom=chrom)
    gate.finish()
    return root
