"""Compare scientific artifacts from matched annotation-memory benchmark runs."""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


def compare_table(left, right):
    def read(path):
        return pd.read_parquet(path) if path.suffix == '.parquet' else pd.read_csv(path, sep='\t')
    a, b = read(left), read(right)
    assert a.shape == b.shape and list(a) == list(b), (left, a.shape, b.shape)
    largest = 0.
    for column in a:
        if column in {'CHR', 'SNP', 'POS', 'BP', 'A1', 'A2', 'category', 'delete_block'} or not pd.api.types.is_numeric_dtype(a[column]):
            pd.testing.assert_series_equal(a[column], b[column], check_dtype=False)
        else:
            x, y = a[column].to_numpy(float), b[column].to_numpy(float)
            np.testing.assert_allclose(x, y, rtol=1e-6, atol=1e-8, equal_nan=True, err_msg=f'{left}: {column}')
            finite = np.isfinite(x) & np.isfinite(y)
            if finite.any():
                largest = max(largest, float(np.max(np.abs(x[finite]-y[finite]))))
    return largest


def compare(left, right, case):
    if case.startswith('direct') or case == 'indexed':
        paths = [Path(name) for name in ('ldscore.baseline.parquet', 'ldscore.query.parquet', 'ldscore.overlap.parquet')]
        a, b = (json.loads((root / 'metadata.json').read_text()) for root in (left, right))
        keys = ['annotation_types', 'artifact_type', 'baseline_columns', 'query_columns', 'counts', 'count_config',
                'overlap_config', 'chromosomes', 'snp_identifier', 'genome_build', 'snp_universe_policy',
                'baseline_row_groups', 'query_row_groups', 'n_baseline_rows', 'n_query_rows']
        for key in keys:
            assert a.get(key) == b.get(key), (case, key)
        if case == 'indexed':
            paths.extend(Path('diagnostics') / name for name in ('gene_list_audit.tsv.gz', 'gene_list_resolution_summary.tsv', 'query_annotation_status.tsv'))
    elif case == 'regression':
        paths = [path.relative_to(left) for path in left.rglob('*') if path.suffix in {'.tsv', '.parquet'}]
        for path in left.glob('diagnostics/query_annotations/*/metadata.json'):
            a = json.loads(path.read_text())
            b = json.loads((right / path.relative_to(left)).read_text())
            for key in ('n_snps', 'n_blocks', 'retained_ld_columns', 'dropped_zero_variance_ld_columns', 'count_kind', 'snp_identifier', 'genome_build'):
                assert a.get(key) == b.get(key), (path, key)
    elif case.startswith('annotate'):
        paths = [path.relative_to(left) for path in left.glob('query.*.annot.gz')]
    else:
        paths = [path.relative_to(left) for path in left.glob('*.tsv')]
    largest = max(compare_table(left / path, right / path) for path in paths)
    return dict(case=case, baseline=left.name, current=right.name, tables=len(paths), max_absolute_difference=largest)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--runs', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    results = []
    for path in sorted(args.runs.glob('current-*.json')):
        record = json.loads(path.read_text())
        if record['exit_code']:
            continue
        case = record['case']
        baseline = args.runs / ('baseline-annotate-bed' if case == 'annotate-gene' else f'baseline-{case}')
        results.append(compare(baseline, path.with_suffix(''), case))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(results, indent=2)+'\n')
    print(json.dumps(results, indent=2))


if __name__ == '__main__':
    main()
