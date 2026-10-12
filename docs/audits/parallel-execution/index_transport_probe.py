"""Focused local index-builder transport/RSS probe; run from the repository root.

Use PYTHONPATH=. python docs/audits/parallel-execution/index_transport_probe.py
with --mode threads|processes --output <new-directory>. This is synthetic
transport evidence, not a genome-wide performance or memory estimate.
"""
import argparse
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
import json
import multiprocessing as mp
import os
from pathlib import Path
import pickle
import resource
import subprocess
import sys
import threading
import time

import numpy as np
import pandas as pd

from ldsc import gene_ldscore_index as core
from ldsc._worker_runtime import initialize_worker, _accelerate_state, native_thread_state
from tests.test_plink_workflow_resolution import write_inputs, index_args


def measured(function, *args):
    _accelerate_state(configure=True)
    result = function(*args)
    result.evidence['probe_pid'] = os.getpid()
    result.evidence['probe_peak_rss_bytes'] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * (1 if sys.platform == 'darwin' else 1024)
    result.evidence['probe_native'] = native_thread_state()
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--mode', choices=['threads', 'processes'], required=True)
    parser.add_argument('--output', type=Path, required=True)
    options = parser.parse_args()
    options.output.mkdir(parents=True, exist_ok=False)
    write_inputs(options.output)
    # Large staged annotation values with a deliberately tiny genotype fixture.
    # This stresses descriptor transport without pretending to test LD scaling.
    rows = 25000
    frames = []
    for chrom in ('21', '22'):
        frame = pd.DataFrame({'CHR': chrom, 'SNP': [f'rs{chrom}_{i}' for i in range(rows)],
                              'POS': np.arange(1, rows + 1) * 100})
        rng = np.random.default_rng(int(chrom))
        values = pd.DataFrame(rng.normal(size=(rows, 64)).astype(np.float32),
                              columns=['base', *[f'annotation{i}' for i in range(63)]])
        frames.append(pd.concat([frame, values], axis=1))
    pd.concat(frames).to_csv(options.output / 'baseline.annot', sep='\t', index=False)
    del frames, frame, values
    initialize_worker()
    sizes = []
    class Executor:
        def __init__(self, max_workers, **kwargs):
            self.pool = (ProcessPoolExecutor(max_workers=max_workers, mp_context=mp.get_context('spawn'),
                                            initializer=initialize_worker)
                         if options.mode == 'processes' else ThreadPoolExecutor(max_workers=max_workers))
        def __enter__(self):
            return self
        def __exit__(self, *args):
            self.pool.shutdown(wait=True)
        def submit(self, function, *args):
            sizes.append(len(pickle.dumps(args, protocol=pickle.HIGHEST_PROTOCOL)))
            return self.pool.submit(measured, function, *args)
    core.ProcessPoolExecutor = Executor
    stop = threading.Event()
    samples = []
    def sample():
        while not stop.is_set():
            pids = [os.getpid(), *[child.pid for child in mp.active_children()]]
            output = subprocess.run(['ps', '-o', 'rss=', '-p', ','.join(map(str, pids))],
                                    text=True, capture_output=True, check=True).stdout
            samples.append(sum(int(line) * 1024 for line in output.splitlines()))
            stop.wait(.05)
    sampler = threading.Thread(target=sample)
    sampler.start()
    started = time.perf_counter()
    try:
        args = index_args(options.output, str(options.output / '1000G.EUR.QC.*'), options.output / 'index')
        path = core.run_build_gene_ldscore_index_from_args(args)
    finally:
        stop.set()
        sampler.join()
    index = core._load_gene_ldscore_index(path, _allow_partial_for_tests=True)
    evidence = json.loads((path / 'diagnostics/build_metadata.json').read_text()) if (path / 'diagnostics/build_metadata.json').exists() else {}
    result = {'mode': options.mode, 'wall_seconds': time.perf_counter() - started,
              'task_pickle_bytes': sizes, 'sampled_peak_sum_rss_bytes': max(samples, default=0),
              'annotation_value_bytes': 2 * rows * 64 * 4, 'index_id': index.index_id,
              'gene_support': index.gene_support.tolist(), 'diagnostics': evidence}
    (options.output / 'probe.json').write_text(json.dumps(result, indent=2))
    print(json.dumps({key: value for key, value in result.items() if key != 'diagnostics'}))


if __name__ == '__main__':
    main()
