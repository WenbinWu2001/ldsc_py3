"""Measure a fresh LD-score run using a preserved CLI command and new batch width.

The original inputs and outputs stay unchanged. RSS is sampled over the process
tree; private and persistent logical disk bytes are reported separately.
"""

import argparse
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import time

from annotation_memory import disk_bytes, tree_rss


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--command-file", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--batch", type=int, required=True)
    args = parser.parse_args()
    tree_rss(os.getpid())
    recorded = json.loads(args.command_file.read_text())["argv"]
    options = recorded[recorded.index("ldscore"):]
    options[options.index("--output-dir") + 1] = str(args.output)
    options[options.index("--query-batch-size") + 1] = str(args.batch)
    args.output.mkdir(parents=True, exist_ok=False)
    temporary = args.output / "library-tmp"
    temporary.mkdir()
    environment = os.environ.copy()
    environment.update(TMPDIR=str(temporary), OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1",
                       MKL_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1", NUMEXPR_NUM_THREADS="1")
    command = [sys.executable, "-m", "ldsc", *options]
    args.output.with_suffix(".command.json").write_text(json.dumps({"argv": command, "blas_threads": 1}, indent=2))
    log_path = args.output.with_suffix(".log")
    peak_rss = peak_private = samples = max_processes = 0
    started = time.monotonic()
    with log_path.open("w") as log:
        process = subprocess.Popen(["/usr/bin/time", "-l", *command], env=environment,
                                   stdout=log, stderr=subprocess.STDOUT)
        while process.poll() is None:
            rss, count = tree_rss(process.pid)
            private, _ = disk_bytes(args.output)
            peak_rss, peak_private = max(peak_rss, rss), max(peak_private, private)
            samples, max_processes = samples + 1, max(max_processes, count)
            time.sleep(.05)
    elapsed = time.monotonic() - started
    private, persistent = disk_bytes(args.output)
    maximum = re.search(r"(\d+)\s+maximum resident set size", log_path.read_text())
    result = dict(batch=args.batch, elapsed_seconds=elapsed, peak_process_tree_rss_bytes=peak_rss,
                  peak_single_process_rss_bytes=int(maximum.group(1)) if maximum else None,
                  peak_private_bytes=max(private, peak_private), remaining_private_bytes=private,
                  persistent_output_bytes=persistent, samples=samples, max_processes=max_processes,
                  exit_code=process.returncode, command=command)
    args.output.with_suffix(".metrics.json").write_text(json.dumps(result, indent=2))
    print(json.dumps({key: value for key, value in result.items() if key != "command"}), flush=True)
    raise SystemExit(process.returncode)


if __name__ == "__main__":
    main()
