"""Benchmark fresh partitioned-h2 CLI runs, including preparation/publication.

Run from an installed development environment with local canonical LD scores
and compatible summary statistics. Outputs are new directories, never replaced.
The POSIX ps sampler measures simultaneous process-tree RSS every 0.1 seconds;
shared pages can be counted in more than one process. No psutil dependency.
"""

import argparse
from collections import defaultdict
import json
import os
from pathlib import Path
import platform
import re
import subprocess
import sys
import time


def tree_rss(pid):
    """Return simultaneous parent-plus-descendant RSS and process count."""
    rows = [tuple(map(int, line.split())) for line in subprocess.check_output(
        ["ps", "-axo", "pid=,ppid=,rss="], text=True).splitlines()]
    members = {pid}
    while True:
        children = {child for child, parent, _ in rows if parent in members}
        expanded = members | children
        if expanded == members:
            break
        members = expanded
    return sum(rss * 1024 for child, _, rss in rows if child in members), len(members)


def phase_times(path):
    totals = defaultdict(float)
    counts = defaultdict(int)
    for phase, seconds in re.findall(r"Phase timing: (.*?) (?:completed in|failed after) ([0-9.]+)s\.", path.read_text()):
        totals[phase] += float(seconds)
        counts[phase] += 1
    return {name: {"calls": counts[name], "seconds": seconds} for name, seconds in totals.items()}


def compare_outputs(reference, actual):
    """Check all scientific tables, uncertainty arrays, statuses and ordering."""
    import pandas as pd

    def read_tsv(root, path):
        return pd.read_csv(root / path, sep="\t")

    for path in ("partitioned_h2.tsv", "diagnostics/query_status.tsv",
                 "diagnostics/query_annotations/manifest.tsv"):
        pd.testing.assert_frame_equal(read_tsv(reference, path), read_tsv(actual, path), rtol=1e-10, atol=1e-12)
    manifest = read_tsv(reference, "diagnostics/query_annotations/manifest.tsv")
    for row in manifest.itertuples():
        for path in (row.summary_path, row.partitioned_h2_full_path):
            pd.testing.assert_frame_equal(read_tsv(reference, path), read_tsv(actual, path), rtol=1e-10, atol=1e-12)
        path = row.coefficient_delete_values_path
        pd.testing.assert_frame_equal(pd.read_parquet(reference / path), pd.read_parquet(actual / path),
                                      rtol=1e-10, atol=1e-12)
        before, after = [json.loads((root / row.metadata_path).read_text()) for root in (reference, actual)]
        for field in ("n_snps", "n_blocks_used", "effective_chisq_max", "retained_ld_columns",
                      "dropped_zero_variance_ld_columns", "effective_snp_identifier", "identity_downgrade_applied"):
            assert before[field] == after[field], (row.query_annotation, field)
    return len(manifest)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ldscore-dir", required=True, type=Path)
    parser.add_argument("--sumstats-file", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--query-batch-size", type=int, default=100)
    parser.add_argument("--repeats", type=int, default=2)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=False)
    environment = os.environ.copy()
    for key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS",
                "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS", "BLIS_NUM_THREADS"):
        environment[key] = "1"
    report = {"platform": platform.platform(), "python": sys.version, "sample_interval_seconds": 0.1,
              "rss_policy": "simultaneous process-tree RSS; shared pages may be counted repeatedly",
              "ldscore_dir": str(args.ldscore_dir.resolve()), "sumstats_file": str(args.sumstats_file.resolve()),
              "query_batch_size": args.query_batch_size, "native_threads": 1, "runs": []}
    reference = None
    for repeat in range(args.repeats):
        for threads in ((1, 2, 4) if repeat % 2 == 0 else (4, 2, 1)):
            output = args.output_dir / f"repeat-{repeat + 1}-workers-{threads}"
            command = [sys.executable, "-m", "ldsc", "partitioned-h2", "--ldscore-dir", str(args.ldscore_dir),
                       "--sumstats-file", str(args.sumstats_file), "--output-dir", str(output),
                       "--threads", str(threads), "--query-batch-size", str(args.query_batch_size),
                       "--continue-on-query-error"]
            peak = processes = samples = 0
            started = time.perf_counter()
            with output.with_suffix(".stdout.log").open("w") as stdout, output.with_suffix(".stderr.log").open("w") as stderr:
                child = subprocess.Popen(command, env=environment, stdout=stdout, stderr=stderr)
                try:
                    while child.poll() is None:
                        rss, count = tree_rss(child.pid)
                        peak, processes = max(peak, rss), max(processes, count)
                        samples += 1
                        time.sleep(0.1)
                    elapsed = time.perf_counter() - started
                finally:
                    if child.poll() is None:
                        child.terminate()
                    child.wait()
            if child.returncode:
                raise RuntimeError(f"Benchmark failed: {output}; exit {child.returncode}")
            assert not list(output.glob(".ldsc-annotation-*"))
            metadata = json.loads((output / "diagnostics/metadata.json").read_text())
            record = {"repeat": repeat + 1, "threads": threads, "effective_workers": metadata["query_workers_effective"],
                      "elapsed_seconds": elapsed, "peak_process_tree_rss_bytes": peak, "max_processes": processes,
                      "samples": samples, "successful_queries": metadata["n_queries_successful"],
                      "failed_queries": metadata["n_queries_failed"],
                      "phases": phase_times(output / "diagnostics/partitioned-h2.log")}
            reference = reference or output
            record["models_compared"] = compare_outputs(reference, output)
            record["scratch_cleaned"] = True
            report["runs"].append(record)
            (args.output_dir / "measurements.json").write_text(json.dumps(report, indent=2) + "\n")
            print(json.dumps({key: value for key, value in record.items() if key != "phases"}), flush=True)


if __name__ == "__main__":
    main()
