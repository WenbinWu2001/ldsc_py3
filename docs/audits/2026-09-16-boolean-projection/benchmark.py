"""Isolated projection experiments; no production module is changed.

Run with the repository's ldsc3-dev interpreter. Results are written next to
this script. Every case/variant/replicate uses a fresh Python subprocess.
"""

import argparse
import json
import os
from pathlib import Path
import random
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parent
BUDGET = 16 * 1024**2


def imports():
    global np, sp, resource, tracemalloc
    import ctypes
    import resource
    import tracemalloc
    import numpy as np
    import scipy.sparse as sp
    from threadpoolctl import threadpool_limits

    threadpool_limits(1)
    if sys.platform == "darwin":
        lib = ctypes.CDLL("/System/Library/Frameworks/Accelerate.framework/Accelerate")
        lib.BLASSetThreading.argtypes = [ctypes.c_uint]
        lib.BLASSetThreading.restype = ctypes.c_int
        assert lib.BLASSetThreading(1) == 0


def highwater():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * (1 if sys.platform == "darwin" else 1024)


def trace_call(fn):
    tracemalloc.start()
    result = fn()
    current, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    return result, {"live_bytes": current, "peak_bytes": peak}


def operator(rows, cols, kind, density, seed=938):
    rng = np.random.default_rng(seed)
    if kind == "dense":
        # Limit generation scratch to one row. Include negative corrected R².
        out = np.empty((rows, cols), dtype=np.float64)
        for i in range(rows):
            raw = rng.random(cols) ** 5
            out[i] = raw - (1 - raw) / 487
        return out
    degree = max(1, round(cols * density))
    indices = np.empty(rows * degree, dtype=np.int32)
    data = np.empty(rows * degree, dtype=np.float64)
    for i in range(rows):
        sl = slice(i * degree, (i + 1) * degree)
        indices[sl] = np.sort(rng.choice(cols, degree, replace=False))
        raw = rng.random(degree, dtype=np.float32) ** 5
        # Match the Parquet reader's existing float32 normalization, then widen.
        data[sl] = raw - (np.float32(1) - raw) / np.float32(487)
    indptr = np.arange(rows + 1, dtype=np.int32) * degree
    return sp.csr_matrix((data, indices, indptr), shape=(rows, cols))


def op_bytes(op):
    return op.data.nbytes + op.indices.nbytes + op.indptr.nbytes if sp.issparse(op) else op.nbytes


def payload(n, q, density, fraction=0):
    """Create compact synthetic retained storage without a full dense staging tile."""
    rng = np.random.default_rng(174)
    nc = round(q * fraction)
    nb = q - nc
    packed = np.empty(((n + 7) // 8, nb), dtype=np.uint8)
    for j in range(nb):
        packed[:, j] = np.packbits(rng.random(n) < density, bitorder="little")
    continuous = np.empty((n, nc), dtype=np.float32)
    for j in range(nc):
        continuous[:, j] = rng.standard_normal(n).astype(np.float32) * (10 if j % 2 else 1)
    # Interleave types in logical order, so grouped operands must restore order.
    cidx = np.linspace(0, q - 1, nc, dtype=int) if nc else np.array([], dtype=int)
    bidx = np.setdiff1d(np.arange(q), cidx)
    return packed, continuous, bidx, cidx


def kernel(case):
    n, q, m = case["n"], case["q"], case["m"]
    blocks = case["blocks"]
    op = operator(m, n, case["operator"], case["density"])
    packed, continuous, bidx, cidx = payload(n, q, case["annotation_density"], case.get("continuous_fraction", 0))
    variant = case["variant"]
    scores = np.zeros((m * blocks, q), dtype=np.float64)
    scores.fill(0)

    def decode():
        binary = np.unpackbits(packed, axis=0, count=n, bitorder="little").view(bool)
        if not len(cidx):
            return binary
        out = np.empty((n, q), dtype=np.float32)
        out[:, bidx], out[:, cidx] = binary, continuous
        return out

    def read():
        if variant == "explicit64":
            out = np.empty((n, q), dtype=np.float64)
            out[:] = decode()
            return out
        if variant == "split":
            return (np.unpackbits(packed, axis=0, count=n, bitorder="little").view(bool), continuous.copy())
        return decode()

    def run():
        values = read()
        for i in range(blocks):
            rows = np.arange(i * m, (i + 1) * m)
            if variant == "split":
                scores[np.ix_(rows, bidx)] += op @ values[0]
                scores[np.ix_(rows, cidx)] += op @ values[1]
            else:
                scores[np.ix_(rows, np.arange(q))] += op @ values

    before = highwater()
    start = time.perf_counter(); run(); cold = time.perf_counter() - start
    first_peak = highwater()
    samples = []
    # Fixed number of same-tile reuses in each call; repeat trials for timing.
    for _ in range(5):
        scores.fill(0)
        start = time.perf_counter(); run(); samples.append(time.perf_counter() - start)
    timing_peak = highwater()
    _, whole_trace = trace_call(run)
    vals, read_trace = trace_call(read)
    if variant != "split":
        traces = []
        for _ in range(3):
            product, tr = trace_call(lambda: op @ vals)
            tr["result_bytes"] = product.nbytes
            tr["extra_peak_bytes"] = tr["peak_bytes"] - product.nbytes
            traces.append(tr)
            del product
        operand = vals.nbytes
    else:
        traces = []
        operand = sum(x.nbytes for x in vals)
    # Numerical comparisons occur AFTER all RSS and allocation measurements.
    expected = op @ decode().astype(np.float64)
    if variant == "split":
        actual = np.empty_like(expected)
        actual[:, bidx], actual[:, cidx] = op @ vals[0], op @ vals[1]
    else:
        actual = op @ vals
    close = np.isclose(actual, expected, rtol=1e-14, atol=1e-14)
    return dict(case=case, setup_peak_rss_bytes=before, first_peak_rss_bytes=first_peak,
                peak_rss_bytes=timing_peak, cold_seconds=cold, seconds=samples,
                operand_bytes=operand, operator_bytes=op_bytes(op), accumulator_bytes=scores.nbytes,
                result_block_bytes=m*q*8, retained_payload_bytes=packed.nbytes+continuous.nbytes,
                read_trace=read_trace, whole_trace=whole_trace, multiply_traces=traces,
                max_abs_difference=float(np.max(np.abs(actual-expected))), result_dtype=str(actual.dtype),
                tolerance_pass=bool(close.all()), tolerance_violations=int(np.count_nonzero(~close)))


def prepare_fixtures():
    imports()
    from ldsc._annotation_storage import ColumnStore
    folder = ROOT / "fixtures"
    folder.mkdir(exist_ok=True)
    for fraction in (0, 0.1):
        tag = "binary" if fraction == 0 else "mixed"
        packed, continuous, bidx, cidx = payload(17380, 1000, .1, fraction)
        bstore = ColumnStore.create(folder / f"{tag}-binary.npy", 17380, tuple(f"q{i}" for i in bidx), dtype=bool)
        # Write bounded row selections, using the actual packed writer.
        for first in range(0, 17380, 2048):
            last = min(first + 2048, 17380)
            rows = np.arange(first, last)
            values = ((packed[rows // 8] >> (rows % 8)[:, None]) & 1).astype(bool)
            bstore.write(first, values)
        stores = [bstore]
        if len(cidx):
            store = ColumnStore.create(folder / f"{tag}-continuous.npy", 17380, tuple(f"q{i}" for i in cidx))
            store.write(0, continuous)
            stores.append(store)
        import pickle
        with (folder / f"{tag}.pkl").open("wb") as out:
            pickle.dump(stores, out)


def projection(case):
    import pickle
    from ldsc._annotation_storage import AnnotationShard
    from ldsc._kernel.ldscore_projection import MappedAnnotations, ProjectionAccumulator

    with (ROOT / "fixtures" / f'{case["source"]}.pkl').open("rb") as inp:
        stores = pickle.load(inp)
    q = case["q"]
    shard = AnnotationShard(ROOT / "fixtures" / "unused-metadata.parquet", tuple(stores), 17380)
    source = MappedAnnotations(shard, np.arange(17380), tuple(f"q{i}" for i in range(q)))

    class NativeRead(ProjectionAccumulator):
        def _read(self, rows, columns):
            return self.annotations.read(rows=rows, columns=[self.annotations.columns[i] for i in columns])

    klass = ProjectionAccumulator if case["variant"] == "explicit64" else NativeRead
    acc = klass(source, n_baseline=0, query_batch_size=q)
    acc.values.fill(0)  # Touch all output pages, as diagonal initialization does.
    op = operator(case["m"], case["n"], case["operator"], case["density"])
    outputs, contributors = np.arange(case["m"]), np.arange(case["n"])
    def run():
        acc._accumulate(op, outputs, contributors)
    before = highwater()
    start = time.perf_counter(); run(); cold = time.perf_counter() - start
    first_peak = highwater()
    samples = []
    for _ in range(4):
        acc.values.fill(0)
        start=time.perf_counter(); run(); samples.append(time.perf_counter()-start)
    peak = highwater()
    acc.values.fill(0)
    _, tr = trace_call(run)
    # Full reference check is deliberately excluded from resource measurements.
    reference = ProjectionAccumulator(source, n_baseline=0, query_batch_size=q)
    reference._accumulate(op, outputs, contributors)
    expected = reference.values[outputs]
    close = np.isclose(acc.values[outputs], expected, rtol=1e-14, atol=1e-14)
    untiled = op @ source.read(rows=contributors, columns=source.columns).astype(np.float64)
    return dict(case=case, setup_peak_rss_bytes=before, first_peak_rss_bytes=first_peak,
                peak_rss_bytes=peak, cold_seconds=cold, seconds=samples, whole_trace=tr,
                operator_bytes=op_bytes(op), accumulator_bytes=acc.values.nbytes,
                max_abs_difference=float(np.max(np.abs(acc.values[outputs]-expected))),
                current_vs_untiled_max_abs=float(np.max(np.abs(expected-untiled))),
                tolerance_pass=bool(close.all()), tolerance_violations=int(np.count_nonzero(~close)),
                result_dtype=str(acc.values.dtype))


def numerical():
    import math
    from ldsc._kernel.ldscore_projection import ArrayAnnotations, ProjectionAccumulator
    rng=np.random.default_rng(1109)
    evidence=[]
    for q in (1,100,1000):
        n,m=257,19
        packed,cont,bidx,cidx=payload(n,q,.1, .1 if q>1 else 0)
        binary=np.unpackbits(packed,axis=0,count=n,bitorder="little").view(bool)
        values=np.empty((n,q),dtype=np.float32)
        values[:,bidx],values[:,cidx]=binary,cont
        for route in ("plink64", "parquet32"):
            raw=rng.random((m,n))**5
            if route=="parquet32": raw=raw.astype(np.float32)
            corrected=(raw-(1-raw)/487).astype(np.float64)
            corrected[rng.random((m,n))<.8]=0
            cols=np.unique(np.r_[0,q-1,bidx[:3],cidx[:3]])
            oracle=np.array([[math.fsum(float(corrected[i,k])*float(values[k,j]) for k in range(n)) for j in cols] for i in range(m)])
            for sparse in (False,True):
                op=sp.csr_matrix(corrected) if sparse else corrected
                explicit=op@values.astype(np.float64)
                native=op@values
                split=np.empty_like(explicit)
                split[:,bidx]=op@binary
                if len(cidx): split[:,cidx]=op@cont
                for label, result in (("explicit",explicit),("native",native),("split",split)):
                    np.testing.assert_allclose(result[:,cols],oracle,rtol=1e-14,atol=1e-14)
                    evidence.append(dict(q=q,route=route,operator="csr" if sparse else "dense",variant=label,
                        negative_r2=int(np.count_nonzero(corrected<0)), signed_continuous=int(np.count_nonzero(cont<0)),
                        max_abs_vs_fsum=float(np.max(np.abs(result[:,cols]-oracle))),
                        max_abs_vs_explicit=float(np.max(np.abs(result-explicit))),dtype=str(result.dtype)))
    # Adversarial float32 regression control: output widening cannot repair sums.
    n=100000
    op=sp.csr_matrix(np.full((1,n),.1,dtype=np.float64))
    a=np.ones((n,1),dtype=bool)
    truth=math.fsum(op.data)
    bad=float((op.astype(np.float32)@a)[0,0]); good=float((op@a)[0,0])
    # Independent diagonal, row order, weight and window oracle via existing seam.
    a=np.array([[1,1,0,-2],[1,0,3,4],[1,1,1,0],[1,0,0,5]],np.float32)
    weights=np.array([1,0,1,0],np.float32)
    r2=np.array([[1,.25,.5,0],[.25,1,.125,-.05],[.5,.125,1,.75],[0,-.05,.75,1]],np.float64)
    oracle=np.array([[math.fsum(float(r2[i,k])*float(x[k]) for k in range(4)) for x in [*a.T,weights]] for i in (2,0)])
    checks=0
    for width in (1,100,1000):
        for sparse in (False,True):
            acc=ProjectionAccumulator(ArrayAnnotations(a,('base','b','c','d')),n_baseline=1,output_rows=[2,0],weight_mask=weights,query_batch_size=width,read_budget_bytes=32)
            if sparse:
                acc.add_diagonal();i,j=np.triu_indices(4,1)
                for first in range(0,len(i),2):
                    sl=slice(first,first+2);acc.add_pairs(i[sl],j[sl],r2[i[sl],j[sl]],np.zeros(4,dtype=int))
            else: acc.add_dense(np.arange(4),np.arange(4),r2)
            np.testing.assert_allclose(acc.values,oracle,rtol=1e-14,atol=1e-14);checks+=1
    return dict(comparisons=evidence, production_oracle_checks=checks,
                negative_control=dict(float32_error=abs(bad-truth),float64_error=abs(good-truth)))


def cases():
    out=[]
    for q in (1,100,1000):
        for kind in ("dense","csr"):
            for ad in (.01,.1,.5,1):
                for blocks in (1,4):
                    n=BUDGET//(8*q)
                    # q=1 full-budget operand uses 8 output rows to cap operator at 128 MiB.
                    m=8 if q==1 else 128
                    out.append(dict(suite="kernel",name="budget",n=n,q=q,m=m,blocks=blocks,
                                    operator=kind,density=.01,annotation_density=ad))
            for n,m,name in ((128,4096,"plink-forward"),(4096,128,"plink-reverse")):
                if kind=="dense":
                    out.append(dict(suite="kernel",name=name,n=n,q=q,m=m,blocks=1,operator=kind,density=1,annotation_density=.1))
        if q>1:
            for kind in ("dense","csr"):
                for fraction in (.1,.5):
                    out.append(dict(suite="kernel",name="mixed",n=BUDGET//(8*q),q=q,m=128,blocks=4,operator=kind,density=.01,annotation_density=.1,continuous_fraction=fraction))
    for density in (.001,.1):
        out.append(dict(suite="kernel",name="csr-density",n=2097,q=1000,m=2097,blocks=4,operator="csr",density=density,annotation_density=.1))
    for q in (1,100,1000):
        for source in ("binary","mixed"):
            if q==1 and source=="mixed": continue
            out.append(dict(suite="projection",name="file-backed-csr",n=4096,m=17380,q=q,operator="csr",density=.01,source=source))
        for n,m,name in ((128,4096,"file-backed-forward"),(4096,128,"file-backed-reverse")):
            out.append(dict(suite="projection",name=name,n=n,m=m,q=q,operator="dense",density=1,source="binary"))
        for n,m,name in ((128,512,"filtered-forward"),(4096,16,"filtered-reverse")):
            out.append(dict(suite="projection",name=name,n=n,m=m,q=q,operator="dense",density=1,source="binary"))
    return out


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--worker");parser.add_argument("--prepare",action="store_true")
    parser.add_argument("--numerical",action="store_true")
    parser.add_argument("--replicates",type=int,default=3)
    parser.add_argument("--limit",type=int)
    parser.add_argument("--output",default="results.jsonl")
    parser.add_argument("--resume",action="store_true")
    args=parser.parse_args()
    if args.prepare: prepare_fixtures();return
    if args.worker:
        imports();case=json.loads(args.worker)
        print(json.dumps(kernel(case) if case["suite"]=="kernel" else projection(case)));return
    if args.numerical:
        imports();print(json.dumps(numerical(),indent=2));return
    schedule=[]
    for case in cases()[:args.limit]:
        variants=("explicit64","native","split") if case.get("continuous_fraction") else ("explicit64","native")
        for rep in range(args.replicates):
            for variant in variants: schedule.append(dict(case,variant=variant,replicate=rep))
    random.Random(431).shuffle(schedule)
    if args.resume:
        completed={json.dumps(json.loads(line)["case"],sort_keys=True) for line in (ROOT/args.output).read_text().splitlines()}
        schedule=[case for case in schedule if json.dumps(case,sort_keys=True) not in completed]
    child_env=os.environ.copy()
    for key in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","VECLIB_MAXIMUM_THREADS"):
        child_env[key]="1"
    with (ROOT/args.output).open("a" if args.resume else "w") as stream:
        for idx,case in enumerate(schedule):
            start=time.perf_counter()
            proc=subprocess.run([sys.executable,str(Path(__file__).resolve()),"--worker",json.dumps(case)],env=child_env,text=True,capture_output=True,check=True)
            result=json.loads(proc.stdout);result["subprocess_seconds"]=time.perf_counter()-start
            stream.write(json.dumps(result)+"\n");stream.flush()
            if idx%12==0: print(f"{idx+1}/{len(schedule)} completed",flush=True)


if __name__=="__main__": main()
