"""Reproducible local memory-refactor comparisons; see the accompanying audit.

Fixture generation is separate from measurement. Each command runs in a fresh
process, using the selected checkout's src directory. The monitor samples the
simultaneous RSS of that process and its descendants, not summed lifetime peaks.
Only the owned benchmark output tree is scanned for temporary disk usage.
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import time

SEED = 20260910


def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2) + "\n")


def prepare(root):
    import numpy as np
    import pandas as pd
    from ldsc._kernel.ref_panel_builder import write_runtime_metadata_sidecar, write_r2_parquet
    from ldsc._kernel.snp_identity import sidecar_identity_sha256

    root.mkdir(parents=True, exist_ok=False)
    rng = np.random.default_rng(SEED)
    n, q, samples = 6000, 1000, 128
    names = [f"pathway{i:04}" for i in range(q)]
    all_meta, all_queries, catalogs = [], [], []
    for chrom in range(1, 4):
        directory = root / f"chr{chrom}"
        directory.mkdir()
        meta = pd.DataFrame(dict(CHR=str(chrom), SNP=[f"rs{chrom}_{i}" for i in range(n)],
                                 POS=10 * np.arange(1, n + 1), A1="A", A2="C",
                                 CM=np.arange(1, n + 1) * .00001, MAF=.3))
        meta.assign(base=1, category=(np.arange(n) % 3 == 0).astype(np.int8)).to_csv(directory / "baseline.annot.gz", sep="\t", index=False)
        values = rng.integers(0, 2, (n, q), dtype=np.int8)
        query = pd.concat([meta, pd.DataFrame(values, columns=names)], axis=1)
        query.to_csv(directory / "query.annot.gz", sep="\t", index=False)
        panel = root / "r2"
        panel.mkdir(exist_ok=True)
        write_runtime_metadata_sidecar(meta, panel / f"chr{chrom}_meta.tsv.gz", genome_build="hg19", snp_identifier="rsid")
        pairs = [(i, i + offset) for i in range(n) for offset in range(1, 5) if i + offset < n]
        left, right = np.asarray(pairs, dtype=np.int32).T
        write_r2_parquet(pair_chunks=[(left, right, .08 / (right-left), np.ones(len(left), dtype=bool))],
                         path=panel / f"chr{chrom}_r2.parquet", genome_build="hg19", n_samples=samples,
                         snp_identifier="rsid", n_snps=n, sidecar_identity_sha256=sidecar_identity_sha256(meta))
        prefix = root / f"panel{chrom}"
        meta[["CHR", "SNP", "CM", "POS", "A1", "A2"]].to_csv(prefix.with_suffix(".bim"), sep="\t", header=False, index=False)
        prefix.with_suffix(".fam").write_text("".join(f"F{i} I{i} 0 0 0 -9\n" for i in range(samples)))
        genotype = rng.binomial(2, .3, (n, samples))
        codes = np.array([0, 2, 3], dtype=np.uint8)[genotype].reshape(n, -1, 4)
        packed = np.sum(codes << np.array([0, 2, 4, 6], dtype=np.uint8), axis=2).astype(np.uint8)
        prefix.with_suffix(".bed").write_bytes(b"\x6c\x1b\x01" + packed.tobytes())
        all_meta.append(meta)
        all_queries.append(query)
        catalogs.append(pd.DataFrame(dict(gene_id=[f"G{chrom}_{i}" for i in range(q)],
                                         gene_name=[f"GENE{chrom}_{i}" for i in range(q)], chrom=str(chrom),
                                         start=10 * np.arange(1, q+1), end=10 * np.arange(1, q+1), genome_build="hg19")))
    meta = pd.concat(all_meta, ignore_index=True)
    prepare_panel_views(root)
    meta.assign(base=1, category=(np.arange(len(meta)) % 3 == 0).astype(np.int8)).to_csv(root / "baseline.annot.gz", sep="\t", index=False)
    pd.concat(all_queries, ignore_index=True).to_csv(root / "query.annot.gz", sep="\t", index=False)
    meta.iloc[::5][["SNP"]].to_csv(root / "regression.tsv", sep="\t", index=False)
    pd.concat(catalogs, ignore_index=True).to_csv(root / "catalog.tsv", sep="\t", index=False)
    for kind in ("genes", "beds"):
        (root / kind).mkdir()
    for i, name in enumerate(names):
        (root / "genes" / f"{name}.txt").write_text("".join(f"G{c}_{i}\n" for c in range(1, 4)))
        (root / "beds" / f"{name}.bed").write_text("".join(f"{c}\t{10*(i+1)-1}\t{10*(i+1)}\n" for c in range(1, 4)))
    prepare_regression(root, q)
    prepare_quantile(root)
    write_json(root / "dimensions.json", dict(seed=SEED, chromosomes=3, snps_per_chromosome=n,
               query_count=q, plink_samples=samples, regression_rows=2000, quantile_rows=120000,
               quantile_annotations=60, standalone_queries=100))


def prepare_panel_views(root):
    import shutil

    for suffix in (".bim", ".bed"):
        with (root / ("panel" + suffix)).open("wb") as output:
            if suffix == ".bed":
                output.write(b"\x6c\x1b\x01")
            for chrom in range(1, 4):
                with (root / f"panel{chrom}{suffix}").open("rb") as source:
                    if suffix == ".bed":
                        source.seek(3)
                    shutil.copyfileobj(source, output)
    shutil.copyfile(root / "panel1.fam", root / "panel.fam")
    (root / "chr1/r2").mkdir(exist_ok=True)
    for path in (root / "r2").glob("chr1_*"):
        shutil.copyfile(path, root / "chr1/r2" / path.name)


def prepare_direct_subset(source, output, query_count):
    """Select a narrower direct-query input, sharing immutable reference files."""
    import pandas as pd

    output.mkdir(parents=True, exist_ok=False)
    columns = pd.read_csv(source / "query.annot.gz", sep="\t", nrows=0).columns
    identity = ["CHR", "SNP", "POS", "A1", "A2", "CM", "MAF"]
    queries = [column for column in columns if column not in identity][:query_count]
    if len(queries) != query_count or query_count < 1:
        raise ValueError("query_count must select an available, nonempty prefix.")
    for name in ("baseline.annot.gz", "regression.tsv", "r2", "panel.bed", "panel.bim", "panel.fam"):
        (output / name).symlink_to((source / name).resolve())
    frame = pd.read_csv(source / "query.annot.gz", sep="\t", usecols=[*identity, *queries])
    frame.to_csv(output / "query.annot.gz", sep="\t", index=False)
    dimensions = json.loads((source / "dimensions.json").read_text())
    write_json(output / "dimensions.json", {**dimensions, "query_count": query_count,
               "subset_scope": "direct three-chromosome LD scoring only"})


def prepare_regression(root, q):
    import numpy as np
    import pandas as pd
    from ldsc import GlobalConfig
    from ldsc.ldscore_calculator import LDScoreResult
    from ldsc.overlap_matrix import LDScoreOverlap
    from ldsc._kernel.overlap import OverlapContribution
    from ldsc.outputs import LDScoreDirectoryWriter, LDScoreOutputConfig
    from ldsc.sumstats_munger import _sumstats_footer_metadata, _write_sumstats_outputs

    n = 2000
    x = np.arange(1, n+1, dtype=float)
    names = [f"pathway{i:04}" for i in range(q)]
    baseline = pd.DataFrame(dict(CHR=np.repeat(["1", "2"], n//2), SNP=[f"s{i}" for i in range(n)],
                                 POS=np.tile(np.arange(1, n//2+1), 2), base=1 + x/100, regression_ld_scores=2.))
    queries = 2 + np.sin(x[:, None] * (.2 + np.arange(q)[None, :] / 997))
    query = pd.concat([baseline[["CHR", "SNP", "POS"]], pd.DataFrame(queries, columns=names)], axis=1)
    block = np.array([[10000.] + [2500.] * q])
    diag = np.full(q, 2500.)
    overlap = LDScoreOverlap.from_contribution(OverlapContribution(block, block.copy(), diag, diag.copy(), 10000, 10000),
                                              baseline_columns=["base"], query_columns=names)
    counts = [dict(group="baseline" if col == "base" else "query", column=col,
                   all_reference_snp_count=10000. if col == "base" else 2500.,
                   common_reference_snp_count=10000. if col == "base" else 2500.) for col in ["base", *names]]
    result = LDScoreResult(baseline, query, counts, ["base"], names, frozenset(), frozenset(baseline.SNP), [],
                          config_snapshot=GlobalConfig(snp_identifier="rsid"), overlap=overlap)
    LDScoreDirectoryWriter().write(result, LDScoreOutputConfig(output_dir=root / "regression-ld"))
    frame = pd.DataFrame(dict(SNP=baseline.SNP, N=10000., Z=np.sqrt(1.1 + .03 * baseline.base + .04 * np.cos(x))))
    _write_sumstats_outputs(frame, output_files={"parquet": str(root / "trait.parquet")}, output_format="parquet",
                           footer_metadata=_sumstats_footer_metadata(GlobalConfig(snp_identifier="rsid"), "synthetic"))


def prepare_quantile(root):
    import numpy as np
    import pandas as pd

    directory = root / "quantile"
    (directory / "ld").mkdir(parents=True)
    (directory / "fitted/diagnostics").mkdir(parents=True)
    n, p = 120000, 60
    rng = np.random.default_rng(SEED + 1)
    names = ["base", *[f"category{i}" for i in range(1, p)]]
    values = rng.integers(0, 2, (n, p)).astype(np.float32)
    values[:, 0] = 1
    meta = pd.DataFrame(dict(CHR=np.repeat([1, 2, 3], n//3), POS=np.tile(np.arange(1, n//3+1), 3),
                             SNP=[f"q{i}" for i in range(n)]))
    pd.concat([meta, pd.DataFrame(values, columns=names)], axis=1).to_csv(directory / "annotations.tsv.gz", sep="\t", index=False)
    meta.assign(target=np.arange(n, dtype=float)/7).to_csv(directory / "target.tsv.gz", sep="\t", index=False)
    meta.assign(MAF=.3).to_csv(directory / "ref.tsv.gz", sep="\t", index=False)
    gram = values.astype(float).T @ values.astype(float)
    pd.DataFrame(dict(row_annotation=np.repeat(names, p), col_annotation=np.tile(names, p),
                      overlap_all_snps=gram.ravel(), overlap_common_snps=gram.ravel())).to_parquet(directory / "ld/ldscore.overlap.parquet", index=False)
    write_json(directory / "ld/metadata.json", dict(artifact_type="ldscore", snp_identifier="rsid", genome_build=None,
               baseline_columns=names, query_columns=[], count_config=dict(common_reference_snp_maf_min=.05, common_reference_snp_maf_operator=">="),
               overlap_config=dict(total_all_reference_snps=n, total_common_reference_snps=n), files=dict(overlap="ldscore.overlap.parquet"),
               counts=[dict(column=c, common_reference_snp_count=float(v)) for c, v in zip(names, values.sum(axis=0))],
               annotation_types={c: "binary" for c in names}))
    write_json(directory / "fitted/diagnostics/metadata.json", dict(artifact_type="partitioned_h2_result", analysis_type="functional_category",
               ldscore_dir=str(directory / "ld"), retained_ld_columns=names))
    tau = np.arange(1, p+1) * 1e-9
    pd.DataFrame(dict(category=names, coefficient=tau, coefficient_se=tau/10, samp_prev=np.nan, pop_prev=np.nan)).to_csv(directory / "fitted/partitioned_h2.tsv", sep="\t", index=False)
    deletes = pd.DataFrame(tau[None, :] * np.linspace(.9, 1.1, 20)[:, None], columns=names)
    deletes.insert(0, "delete_block", np.arange(20))
    deletes.to_parquet(directory / "fitted/diagnostics/coefficient_delete_values.parquet", index=False)


def worker(args):
    import ldsc
    import platform
    import numpy, pandas, scipy, pyarrow
    from ldsc import cli

    source, output = args.inputs, args.output
    write_json(output / "worker.json", dict(import_path=ldsc.__file__, python=sys.version, platform=platform.platform(),
               versions={package.__name__: package.__version__ for package in (numpy, pandas, scipy, pyarrow)}))
    kind = args.case
    common = ["--output-dir", str(output), "--log-level", "WARNING"]
    if kind.startswith("direct"):
        one = kind.endswith("-one")
        base = source / "chr1" if one else source
        argv = ["ldscore", *common, "--baseline-annot-sources", str(base / "baseline.annot.gz"),
                "--query-annot-sources", str(base / "query.annot.gz"), "--snp-identifier", "rsid",
                "--regr-snps-file", str(source / "regression.tsv"), "--regr-snps-exclude-regions", "none",
                "--ld-wind-snps", "4", "--threads", str(args.workers)]
        if "plink" in kind:
            argv += ["--plink-prefix", str(source / ("panel1" if one else "panel"))]
        else:
            argv += ["--r2-dir", str(base / "r2")]
        if args.batch:
            argv += ["--query-batch-size", str(args.batch)]
    elif kind.startswith("annotate"):
        route = "gene-list" if kind.endswith("gene") else "bed"
        folder, suffix = ("genes", "txt") if route == "gene-list" else ("beds", "bed")
        paths = sorted((source / folder).glob(f"*.{suffix}"))[:100]
        argv = ["annotate", *common, "--baseline-annot-sources", str(source / "baseline.annot.gz"),
                f"--query-annot-{route}-sources", *map(str, paths), "--snp-identifier", "rsid", "--padding-bp", "0"]
        if route == "gene-list":
            argv += ["--gene-coordinate-file", str(source / "catalog.tsv"), "--genome-build", "hg19"]
    elif kind == "regression":
        argv = ["partitioned-h2", *common, "--ldscore-dir", str(source / "regression-ld"),
                "--sumstats-file", str(source / "trait.parquet"), "--n-blocks", "20"]
        if args.batch:
            argv += ["--query-batch-size", str(args.batch)]
    elif kind == "quantile":
        source = source / "quantile"
        argv = ["quantile-h2", *common, "--partitioned-h2-result-dir", str(source / "fitted"),
                "--baseline-annot-sources", str(source / "annotations.tsv.gz"), "--target-annot-sources", str(source / "target.tsv.gz"),
                "--target-annotation", "target", "--ref-metadata-sources", str(source / "ref.tsv.gz"), "--num-quantiles", "10"]
    elif kind == "indexed":
        from ldsc import gene_ldscore_index as index
        timings = {}
        for function, label in (("_load_gene_ldscore_index", "validation_seconds"),
                                ("assemble_indexed_ld_scores", "operator_multiplication_seconds")):
            original = getattr(index, function)
            def timed(*a, _original=original, _label=label, **kw):
                started = time.perf_counter()
                try:
                    return _original(*a, **kw)
                finally:
                    timings[_label] = timings.get(_label, 0.) + time.perf_counter() - started
            setattr(index, function, timed)
        kwargs = dict(output_dir=output, query_gene_list_sources=sorted((source / "genes").glob("*.txt")), _allow_partial_for_tests=True)
        if args.batch:
            kwargs["query_batch_size"] = args.batch
        started = time.perf_counter()
        index.run_indexed_ldscore(source / "index", **kwargs)
        timings["workflow_seconds"] = time.perf_counter() - started
        timings["outside_validation_and_multiplication_seconds"] = timings["workflow_seconds"] - sum(
            timings[key] for key in ("validation_seconds", "operator_multiplication_seconds"))
        record = json.loads((output / "worker.json").read_text())
        write_json(output / "worker.json", {**record, "timings": timings})
        return
    else:
        raise ValueError(kind)
    cli.main(argv)


def build_index(source):
    from ldsc import gene_ldscore_index as index
    args = index.build_parser().parse_args([
        "--output-dir", str(source / "index"), "--baseline-annot-sources", str(source / "baseline.annot.gz"),
        "--plink-prefix", str(source / "panel"), "--gene-coordinate-file", str(source / "catalog.tsv"),
        "--genome-build", "hg19", "--snp-identifier", "rsid", "--ld-wind-cm", ".00004",
        "--padding-bp", "0", "--gene-exclude-regions", "none", "--regr-snps-exclude-regions", "none",
        "--regr-snps-file", str(source / "regression.tsv"), "--log-level", "WARNING"])
    args._test_chromosomes = ("1", "2", "3")
    index.run_build_gene_ldscore_index_from_args(args)


def tree_rss(pid):
    rows = subprocess.check_output(["/bin/ps", "-axo", "pid=,ppid=,rss="], text=True)
    processes = {p: (parent, rss) for p, parent, rss in (map(int, row.split()) for row in rows.splitlines())}
    descendants = {pid}
    while True:
        expanded = descendants | {p for p, (parent, _) in processes.items() if parent in descendants}
        if expanded == descendants:
            break
        descendants = expanded
    return sum(processes.get(p, (0, 0))[1] for p in descendants) * 1024, len(descendants)


def disk_bytes(output):
    private = public = 0
    for directory, _, files in os.walk(output):
        path = Path(directory)
        is_private = any(part.startswith(".") or part == "library-tmp" for part in path.relative_to(output).parts)
        for name in files:
            try:
                size = (path / name).stat().st_size
            except FileNotFoundError:
                continue
            if is_private or name.startswith("."):
                private += size
            elif name not in {"worker.json"}:
                public += size
    return private, public


def measure(args):
    args.output.mkdir(parents=True, exist_ok=False)
    scratch = args.output / "library-tmp"
    scratch.mkdir()
    environment = os.environ.copy()
    environment.update(PYTHONPATH=str(args.checkout / "src"), TMPDIR=str(scratch),
                       OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1")
    command = [sys.executable, str(Path(__file__).resolve()), "worker", "--inputs", str(args.inputs),
               "--output", str(args.output), "--case", args.case, "--batch", str(args.batch), "--workers", str(args.workers)]
    peak_rss = peak_private = max_processes = samples = 0
    start = time.monotonic()
    with args.output.with_suffix(".log").open("w") as log:
        process = subprocess.Popen(command, cwd=args.checkout, env=environment, stdout=log, stderr=subprocess.STDOUT)
        while process.poll() is None:
            rss, count = tree_rss(process.pid)
            private, _ = disk_bytes(args.output)
            peak_rss, peak_private = max(peak_rss, rss), max(peak_private, private)
            max_processes, samples = max(max_processes, count), samples + 1
            time.sleep(.05)
    elapsed = time.monotonic() - start
    private, public = disk_bytes(args.output)
    record = dict(case=args.case, batch=args.batch or None, workers=args.workers, checkout=str(args.checkout), inputs=str(args.inputs),
                  revision=subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=args.checkout, text=True).strip(),
                  elapsed_seconds=elapsed, peak_process_tree_rss_bytes=peak_rss, peak_private_bytes=max(private, peak_private),
                  persistent_output_bytes=public, remaining_private_bytes=private, samples=samples,
                  max_processes=max_processes, exit_code=process.returncode, output=str(args.output))
    worker_record = args.output / "worker.json"
    if worker_record.exists():
        record["worker"] = json.loads(worker_record.read_text())
    write_json(args.output.with_suffix(".json"), record)
    print(json.dumps(record), flush=True)
    if process.returncode:
        raise SystemExit(process.returncode)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("prepare", "prepare-direct-subset", "build-index", "measure", "matrix", "worker"))
    parser.add_argument("--inputs", required=True, type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--checkout", type=Path)
    parser.add_argument("--case")
    parser.add_argument("--batch", type=int, default=0)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--query-count", type=int, default=100)
    args = parser.parse_args()
    if args.mode == "prepare":
        prepare(args.inputs)
    elif args.mode == "prepare-direct-subset":
        prepare_direct_subset(args.inputs, args.output, args.query_count)
    elif args.mode == "build-index":
        build_index(args.inputs)
    elif args.mode == "worker":
        worker(args)
    elif args.mode == "matrix":
        jobs = []
        for case in ("direct-r2", "direct-plink", "indexed", "regression", "annotate-bed", "quantile"):
            jobs.extend([(case, "baseline", 0, 1), (case, "current", 1000, 1)])
            if case in {"direct-r2", "direct-plink", "indexed", "regression"}:
                jobs.append((case, "current", 32, 1))
            if case.startswith("direct"):
                jobs.extend([(case, "current", 32, 2), (case + "-one", "baseline", 0, 1), (case + "-one", "current", 32, 1)])
        jobs.append(("annotate-gene", "current", 1000, 1))
        for case, revision, batch, workers in jobs:
            name = f"{revision}-{case}" + (f"-b{batch}-w{workers}" if batch else "")
            output = args.output / name
            if output.with_suffix(".json").exists():
                record = json.loads(output.with_suffix(".json").read_text())
                if record["exit_code"] == 0:
                    continue
                raise RuntimeError(f"A failed run needs inspection before retrying: {output}")
            checkout = args.checkout if revision == "current" else args.checkout / ".worktrees/memory-baseline"
            measure(argparse.Namespace(inputs=args.inputs, output=output, checkout=checkout, case=case, batch=batch, workers=workers))
    else:
        measure(args)


if __name__ == "__main__":
    main()
