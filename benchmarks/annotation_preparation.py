"""Reproducible local annotation-preparation timing and scratch measurements.

Generate inputs in a separate process, then run each measured case in a fresh
process with the desired checkout's src on PYTHONPATH. Example commands and
matched results are in docs/audits/annotation-memory/annotation-preparation.md.
Only the supplied output directory is monitored. No reference panel is needed.
"""

import argparse
from collections import Counter
from functools import wraps
import json
import os
from pathlib import Path
import platform
import resource
import sqlite3
import sys
import threading
import time


CASES = {
    "sharded": (22, 8000, 96, 0, True),
    "wide-sharded": (3, 6000, 8, 1000, True),
    "wide-whole": (3, 6000, 8, 1000, False),
    "small": (2, 100, 4, 2, True),
}


def generate(root):
    import numpy as np
    import pandas as pd

    root.mkdir(parents=True, exist_ok=False)
    for label, (chromosomes, rows, baseline, query, sharded) in CASES.items():
        directory = root / label
        directory.mkdir()
        paths = {"baseline": [], "query": []}
        for kind, width in (("baseline", baseline), ("query", query)):
            if not width:
                continue
            # A fixed formula permits exact value checks without retaining an
            # expected matrix or relying on the implementation under test.
            for chromosome in range(1, chromosomes + 1):
                snps = [f"rs{chromosome}_{i}" for i in range(rows)]
                snps[0] = "cross_chromosome_duplicate"
                metadata = pd.DataFrame({"CHR": chromosome, "SNP": snps,
                                         "POS": np.arange(rows) + 1})
                columns = tuple(f"{kind}{i}" for i in range(width))
                values = ((np.arange(rows)[:, None] * 7 + np.arange(width)) % 41 - 20) / 8
                frame = pd.concat([metadata, pd.DataFrame(values, columns=columns)], axis=1)
                path = directory / f"{kind}.{chromosome if sharded else 'all'}.annot.gz"
                frame.to_csv(path, sep="\t", index=False,
                             mode="w" if sharded or chromosome == 1 else "a",
                             header=sharded or chromosome == 1,
                             compression={"method": "gzip", "compresslevel": 1, "mtime": 0})
                if str(path) not in paths[kind]:
                    paths[kind].append(str(path))
        (directory / "inputs.json").write_text(json.dumps(paths, indent=2) + "\n")


def measure(args):
    import numpy as np
    import pandas as pd
    import pyarrow
    import ldsc
    from ldsc import _annotation_sources as sources
    from ldsc import _annotation_identity as identity_module
    from ldsc._annotation_identity import DiskIdentityIndex
    from ldsc._annotation_storage import AnnotationWorkspace

    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    paths = json.loads((args.inputs / args.case / "inputs.json").read_text())
    phases, calls, sizes, sql = Counter(), Counter(), {}, Counter()
    if args.sql_batch_size is not None:
        identity_module._LOOKUP_KEYS = args.sql_batch_size
    if args.transaction_rows is not None:
        identity_module._TRANSACTION_ROWS = args.transaction_rows

    def timed(owner, name):
        if not hasattr(owner, name):
            return
        original = getattr(owner, name)

        @wraps(original)
        def wrapped(*positional, **keywords):
            started = time.perf_counter()
            try:
                if name in ("add", "select"):
                    sizes.setdefault(name, []).append(len(positional[1]))
                return original(*positional, **keywords)
            finally:
                phases[name] += time.perf_counter() - started
                calls[name] += 1
        setattr(owner, name, wrapped)

    for owner, names in ((sources, ("_scan", "_finish_shard", "_select_rows", "_write_values")),
                         (DiskIdentityIndex, ("add", "select"))):
        for name in names:
            timed(owner, name)

    alignment_name = "_aligned_metadata" if hasattr(sources, "_aligned_metadata") else "_aligned_rows"
    original_alignment = getattr(sources, alignment_name)

    def aligned(*positional, **keywords):
        iterator = iter(original_alignment(*positional, **keywords))
        while True:
            started = time.perf_counter()
            try:
                chunk = next(iterator)
            except StopIteration:
                phases[alignment_name] += time.perf_counter() - started
                return
            phases[alignment_name] += time.perf_counter() - started
            calls[alignment_name] += 1
            yield chunk

    setattr(sources, alignment_name, aligned)

    original_connect = sqlite3.connect

    class MeasuredConnection(sqlite3.Connection):
        def execute(self, statement, *positional, **keywords):
            sql[statement.split()[0].upper()] += 1
            return super().execute(statement, *positional, **keywords)

        def commit(self):
            sql["commit_calls"] += 1
            return super().commit()

    sqlite3.connect = lambda *a, **kw: original_connect(*a, factory=MeasuredConnection, **kw)
    files, opens = set(), Counter()
    prefix = str(output) + os.sep

    def audit(event, values):
        if event != "open" or not isinstance(values[0], (str, bytes)):
            return
        path = os.fsdecode(values[0])
        if path.startswith(prefix):
            opens[str(values[1])] += 1
            if values[2] & os.O_CREAT:
                files.add(path)

    sys.addaudithook(audit)
    peak = {"files": 0, "bytes": 0}
    stop = threading.Event()

    def sample():
        while not stop.is_set():
            count = size = 0
            for directory, _, names in os.walk(output):
                for name in names:
                    try:
                        size += os.stat(os.path.join(directory, name)).st_size
                        count += 1
                    except FileNotFoundError:
                        pass
            peak["files"] = max(peak["files"], count)
            peak["bytes"] = max(peak["bytes"], size)
            stop.wait(.025)

    monitor = threading.Thread(target=sample, daemon=True)
    monitor.start()
    with AnnotationWorkspace(output) as workspace:
        started = time.perf_counter()
        prepared = sources.prepare_annotation_sources(
            workspace, paths["baseline"], paths["query"], mode="rsid", chunk_rows=args.chunk_rows,
            declared_chromosomes=({p: Path(p).name.split(".")[1]
                                  for group in paths.values() for p in group} if args.declared else None),
        )
        elapsed = time.perf_counter() - started
        # Capture the measured preparation peak before detached verification
        # arrays can increase this process's lifetime RSS.
        rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        final_bytes = sum(p.stat().st_size for p in workspace.path.rglob("*") if p.is_file())
        stop.set()
        monitor.join()
        open_counts, created = dict(opens), len(files)
        chromosomes, rows, baseline, query, _ = CASES[args.case]
        assert list(prepared.shards) == [str(i) for i in range(1, chromosomes + 1)]
        assert prepared.drops.n_rows == chromosomes
        for shard in prepared.shards.values():
            assert shard.n_rows == rows - 1
            metadata = shard.metadata()
            np.testing.assert_array_equal(metadata.POS, np.arange(2, rows + 1))
            for kind, width in (("baseline", baseline), ("query", query)):
                for start in range(0, width, 128):
                    indices = np.arange(start, min(width, start + 128))
                    values = shard.read(columns=[f"{kind}{i}" for i in indices])
                    expected = ((np.arange(1, rows)[:, None] * 7 + indices) % 41 - 20) / 8
                    np.testing.assert_array_equal(values, expected)
                    assert values.dtype == np.float32
        drops = pd.concat(prepared.drops.frames(), ignore_index=True)
        assert drops.CHR.tolist() == list(prepared.shards)
        assert drops.reason.eq("duplicate_identity").all()
    assert not list(output.glob(".ldsc-annotation-*"))
    result = {
        "case": args.case, "declared": args.declared, "chunk_rows": args.chunk_rows,
        "sql_batch_size": args.sql_batch_size, "transaction_rows": args.transaction_rows,
        "elapsed_seconds": elapsed, "phase_seconds": dict(phases), "phase_calls": dict(calls),
        "identity_chunk_rows": {k: {"min": min(v), "max": max(v), "count": len(v)} for k, v in sizes.items()},
        "sqlite_statements": dict(sql), "scratch_peak_sampled": peak,
        "scratch_files_created_python": created, "scratch_open_calls_python": open_counts,
        "prepared_bytes": final_bytes,
        "peak_rss_bytes": rss if sys.platform == "darwin" else rss * 1024,
        "verification": "all values, metadata order, dtype, cross-chromosome drops and owner cleanup passed",
        "import_path": ldsc.__file__, "platform": platform.platform(), "python": sys.version,
        "versions": {m.__name__: m.__version__ for m in (np, pd, pyarrow)},
        "sqlite_version": sqlite3.sqlite_version,
    }
    (output / "result.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    gen = sub.add_parser("generate")
    gen.add_argument("inputs", type=Path)
    run = sub.add_parser("measure")
    run.add_argument("inputs", type=Path)
    run.add_argument("output", type=Path)
    run.add_argument("--case", choices=CASES, required=True)
    run.add_argument("--chunk-rows", type=int)
    run.add_argument("--declared", action="store_true")
    run.add_argument("--sql-batch-size", type=int)
    run.add_argument("--transaction-rows", type=int)
    args = parser.parse_args()
    generate(args.inputs) if args.command == "generate" else measure(args)
