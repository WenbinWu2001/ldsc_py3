"""Reproduce the quantitative-count gap found in the completion audit.

Run from the repository root with ``PYTHONPATH=src python
docs/audits/annotation-memory/completion_probe.py OUTPUT_DIRECTORY``. Scratch
is created only beneath that directory and removed before returning. This is
an audit probe, not a replacement for a failing regression test and repair.
"""

import argparse
import json
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace

import numpy as np
import pandas as pd

from ldsc._annotation_storage import ColumnStore
from ldsc._kernel.ldscore import compute_counts
from ldsc._kernel.ldscore_projection import MappedAnnotations
from ldsc._kernel.overlap import annotation_statistics
from ldsc._quantile_inputs import _validate_aggregates
from ldsc.errors import LDSCInputError


def probe(output_dir):
    """Return observed reductions and the unchanged quantile count-gate result."""
    output_dir.mkdir(parents=True, exist_ok=True)
    n_rows, n_queries = 6000, 1000
    names = ('base', 'baseline2', *(f'q{i}' for i in range(n_queries)))
    values = np.empty((n_rows, len(names)), dtype=np.float32, order='F')
    values[:, 0] = 1
    values[:, 1:] = (np.arange(n_rows) % 2)[:, None]
    values[:, 2] = -np.random.default_rng(182).normal(size=n_rows).astype(np.float32)
    metadata = pd.DataFrame({'MAF': np.full(n_rows, 0.3)})
    exact_sum = float(values[:, 2].sum(dtype=np.float64))
    prior_all, prior_common = compute_counts(metadata, pd.DataFrame(values, columns=names))
    reductions = {'prior_formula': (prior_all, prior_common)}
    with TemporaryDirectory(prefix='count-probe-', dir=output_dir) as directory:
        store = ColumnStore.create(Path(directory) / 'annotations.npy', n_rows, names)
        store.write(0, values)
        roundtrip_equal = bool(np.array_equal(store.read(), values))
        annotations = MappedAnnotations(store, np.arange(n_rows), names)
        for batch_size in (32, 1000):
            all_counts, common_counts, _, _ = annotation_statistics(
                metadata, annotations, 2, query_batch_size=batch_size)
            reductions[f'batch_{batch_size}'] = (all_counts, common_counts)
    report = {'rows': n_rows, 'queries': n_queries, 'seed': 182,
              'exact_sum': exact_sum, 'staged_values_equal': roundtrip_equal,
              'reductions': {}}
    for label, (all_counts, common_counts) in reductions.items():
        counts = {'counts': [{'column': 'q0', 'common_reference_snp_count': float(common_counts[2])}]}
        try:
            _validate_aggregates(SimpleNamespace(annotation_names=['q0']), counts,
                                 [exact_sum], None)
            validation = 'passed'
        except LDSCInputError as error:
            validation = str(error)
        report['reductions'][label] = {
            'all_reference_count': float(all_counts[2]),
            'common_reference_count': float(common_counts[2]),
            'quantile_count_validation': validation,
        }
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output_dir', type=Path)
    print(json.dumps(probe(parser.parse_args().output_dir), indent=2))
