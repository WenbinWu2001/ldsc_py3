# Boolean annotation operands during LD-score projection

Last updated on: 2026-09-16

The evidence does **not** support a general Boolean-operand replacement as an 87.5% reduction in projection memory. Boolean operands preserve the current numerical results when the operator remains float64, but both installed multiplication implementations make a complete float64 operand copy on every call. Wide CSR projections consequently use more temporary memory. A narrower dense case does save memory by reducing the overlap between successive annotation tiles; that is a possible follow-up hypothesis, not an implemented optimization or a demonstrated whole-job benefit.

All work was local on `restructure`, HEAD `0984a2b05734806923565adcd54db00ab0688cd0`, with the user's existing working-tree changes. Workspace and repository `AGENTS.md`, relevant current documentation, the active annotation-memory plan, and the float32 sparse-projection lesson were read. Production code, current design documents, numerical tolerances, and inputs were not edited. No HPC access, full-size job, commit, or push occurred. Benchmark subclasses exist only in this audit directory.

## Current paths and allocations

| Stage | Verified behavior and source |
| --- | --- |
| Stored annotations | Packed Boolean columns and separate dense float32 continuous columns. `ColumnStore.read()` allocates a detached selection and decodes bounded row runs; row indices, sorting state and small decoding buffers are additional allocations. See [`_annotation_storage.py`](../../../src/ldsc/_annotation_storage.py), `ColumnStore.read`, and [current memory design, “Encoding, identities, and order”](../../current/annotation-memory-design.md#encoding-identities-and-order). |
| Internal file-backed reads | `AnnotationShard.read()` returns Boolean for binary-only selections. Multiple stores require grouping and assembly in requested order. A selection containing continuous columns is float32, even if most columns are binary. A grouped read may overlap its assembled result and a store selection. |
| Public and in-memory reads | `AnnotationBundle.read()` always converts to float32. `MemoryAnnotationShard.read()` also returns float32, including its internal binary-only selections. These are different paths from the file-backed Boolean reader. See [`_annotation_bundle.py`](../../../src/ldsc/_annotation_bundle.py), `AnnotationBundle.read`, and [`_annotation_memory.py`](../../../src/ldsc/_annotation_memory.py), `MemoryAnnotationShard.read`. |
| Reference alignment | [`ldscore_calculator.py`](../../../src/ldsc/ldscore_calculator.py), `_kernel_annotation_bundle`, borrows the shard through `MappedAnnotations`; reference adapters add retained-row mappings. Production projection does not need to route through the public float32 bundle read. |
| Explicit widening | [`ldscore_projection.py`](../../../src/ldsc/_kernel/ldscore_projection.py), `ProjectionAccumulator._read`, allocates float64 **before** reading the source, then assigns source values and the optional regression-weight column. The float64 destination and decoded source therefore overlap. |
| Dense operators | [`plink_bed.py`](../../../src/ldsc/_kernel/plink_bed.py), `nextSNPs`, `__corSumVarBlocks__`, and `__l2_unbiased__`, use float64 genotype/correlation calculations and retain negative unbiased R². `add_dense()` also makes an output-filtered correlation copy before `_accumulate()`. |
| CSR operators | `ProjectionAccumulator.add_pairs()` builds one bidirectional CSR operator with explicitly float64 data. The Parquet reader's existing normalized/corrected float32 R² values are widened here; changing annotation operands does not recover or alter that earlier normalization. |
| Multiplication and accumulation | `_accumulate()` reads once per contributor tile, then reuses that operand across output blocks. CSR slicing allocates another operator block; a dense basic slice is a view. Advanced indexing gathers a float64 accumulator tile **before** evaluating `block @ values`; the gather and product overlap. The full accumulator stays float64. |
| Tile transitions | During `values = self._read(...)`, the previous `values` remains alive until the new read returns. Old and new float64 tiles can overlap with the source read. The previous operator block also remains bound until replaced. |
| Output and other phases | `compute_from_parquet()` and the PLINK driver retain existing float32 publication. Pair buffering, genotype computation, counts, overlap, aggregation and batch retention are unchanged and outside this investigation. |

The current tile limit is exactly `16 * 1024**2` bytes, using `floor(budget / (8 * batch_columns))` for **both** contributor and output rows. This is not a process-memory limit. Current configuration defaults are `query_batch_size=100` and `snp_batch_size=128`; the shared baseline-plus-weight batch is separate from the query batch. See `ProjectionAccumulator.__init__/batches/_accumulate` and [`config.py`](../../../src/ldsc/config.py), `LDScoreConfig`.

| Batch columns | Maximum tile rows | Float64 operand bytes | Boolean operand bytes | Difference, MiB |
| ---: | ---: | ---: | ---: | ---: |
| 1 | 2,097,152 | 16,777,216 | 2,097,152 | 14.000 |
| 54: 53 baselines + weight | 38,836 | 16,777,152 | 2,097,144 | 14.000 |
| 100 | 20,971 | 16,776,800 | 2,097,100 | 14.000 |
| 1,000 | 2,097 | 16,776,000 | 2,097,000 | 13.999 |

Thus the preliminary 87.5% calculation is correct **for the retained operand alone at identical dimensions**. Increasing Boolean tile dimensions eightfold would also increase the hidden float64 conversion eightfold; these experiments keep the existing dimensions.

Documentation caveat: [the current annotation-memory design, “Direct LD-score computation”](../../current/annotation-memory-design.md#direct-ld-score-computation) reflects selected output rows and bounded reads. Portions of [the Parquet accumulation document, §5 and §7](../../current/ldscore-parquet-accumulation.md) describe the older full-matrix/two-product realization, assert no block tiling, and quote an approximate 0.45-GiB scatter footprint. Current production uses `ProjectionAccumulator.add_pairs()` and tiled `_accumulate()`. Those older allocation statements are not used as current RSS bounds. Mathematical linearity also does not imply bit-identical floating-point sums under different tilings.

## Library behavior: research and direct observation

The measured environment is Python 3.13.13, NumPy 2.4.5, SciPy 1.17.1, macOS 15.3.2 arm64, with NumPy linked to Apple Accelerate. Each worker set native thread limits to one and successfully called Accelerate's `BLASSetThreading(1)`. Exact configuration and installed SciPy Python methods are saved in [environment.json](environment.json). The sandbox did not permit a hardware-model query; no particular CPU model is assumed.

**NumPy:** installed `np.matmul.types` has homogeneous `dd->d` and `ff->f` loops, not a separate float64/Boolean mixed loop. The version-matched [generalized-ufunc implementation](https://raw.githubusercontent.com/numpy/numpy/v2.4.5/numpy/_core/src/umath/ufunc_object.c), `PyUFunc_GeneralizedFunctionInternal`, permits operand copying and explains why its iterator cannot use ordinary elementwise buffering. The observed Boolean multiplication allocates a complete double tile, consistent with this dispatch. This is not evidence that every BLAS implementation has identical additional workspace.

**SciPy:** the installed `_cs_matrix._matmul_multivector` allocates the promoted result and passes `other.ravel()` to `_sparsetools`; one-column operands use the vector path. The version-matched [C++ wrapper](https://raw.githubusercontent.com/scipy/scipy/v1.17.1/scipy/sparse/sparsetools/sparsetools.cxx), `call_thunk` and `c_array_from_object`, resolves a common numeric type and converts mismatching inputs to contiguous arrays of that type. The [CSR kernel](https://raw.githubusercontent.com/scipy/scipy/v1.17.1/scipy/sparse/sparsetools/csr.h), `csr_matvec/csr_matvecs`, then consumes operator, operand and output in that same numeric type. Boolean operands are converted to double before arithmetic; there is no Boolean-aware SpMM in this path. Noncontiguous operands can additionally require flattening copies; production selected reads and these probes are contiguous.

For every tested full-budget width, the allocation probe found approximately **16 MiB extra beyond the result** for `float64 @ bool`, versus only small object bookkeeping for `float64 @ float64`. The same extra allocation appeared in each of three successive calls, for dense and CSR. A returned operand can therefore remain Boolean at the Python level while the multiplication uses a temporary float64 copy. Source inspection and repeated allocation observations establish repeated conversion; the inference is not based merely on the output dtype.

Verified promotions are recorded in [context.json](context.json): float64 operators with Boolean **or float32** annotations produce float64 results; float32 operators with Boolean annotations produce float32. Widening the accumulator afterward cannot repair an already rounded float32 product. The present explicit float64 annotation conversion also protects callers supplying float32 operators; removing that protection would require an explicit float64-operator invariant.

## Benchmark design and limits

The complete experiment contains **81 workload configurations, 170 configuration/variant combinations, and 510 fresh-process measurements**, three process replicates each, in randomized order. Successful measured subprocesses took 362.8 seconds in total, including setup and subsequent validation. Runtime tables use the median of each process's warm-call median; kernel processes have five warm calls, projection processes four. Cold timings and every raw sample are retained. Three replicates describe variation; they are not confidence intervals.

- Binary annotation densities: 1%, 10%, 50%, and 100%. Dense operators are fully populated. CSR operators use fixed random row degrees approximating densities 0.1%, 1%, and 10%. The primary CSR density is 1%; dense rows have no sparsity reduction. Synthetic operator locality does not reproduce a genomic banded matrix or all 16 million buffered pairs.
- Widths: 1, 100, 1,000. Budget-stress cases fill the actual 16-MiB float64 tile. Width 1 uses only eight output rows to keep its dense operator at 128 MiB; its 2.1-million-contributor tile is a stress case, not a chr22 window.
- Dense PLINK-shaped cases: 128 contributors × 4,096 outputs and the reverse, plus output-filtered cases with 512 or 16 outputs. A read-only [BIM probe](context_probe.py) found 141,123 SNPs; rounded 1-cM left-window sizes were median 2,048, p90 4,992, p95 6,272, maximum 8,192. This precedes annotation, MAF and output filtering. Output-filtered dimensions approximate the historical 17,380/141,123 output fraction; they are not an instrumented current full job.
- Reuse: kernel experiments apply the same read once to one or four controlled output blocks. The actual `_accumulate()` experiments use production slicing, read order and scatter with file-backed `ColumnStore`/`AnnotationShard` fixtures. Their CSR operator has 17,380 outputs × 4,096 contributors and 712,580 stored entries. At width 1,000, two contributor tiles and nine output blocks cause **18 implicit casts**, versus two explicit conversions.
- Mixed selections: interleaved logical columns with 10% or 50% signed continuous values, normalized to float32. `native` retains Boolean only for binary selections, and float32 for mixed selections. An exploratory `split` variant separately multiplies binary and continuous groups and scatters them back into logical order. It is a benchmark option, not production integration.

Peak RSS is the OS `ru_maxrss` high-water mark through setup and the uninstrumented timed calls in each fresh child, in bytes on this macOS host. It includes imports, retained storage, operator generation, operator, touched accumulator and allocator history. It is captured **before** allocation tracing and correctness work, whose larger validation allocations are deliberately excluded. It cannot be reset to a phase-only high-water mark, and a lower observed RSS is not automatically a lower live-array requirement. RSS is not USS or scheduler job accounting.

Tracemalloc runs separately from timing, starting after persistent operator/storage/accumulator creation. It measures additional tracked live-allocation peaks, not RSS or cumulative allocation traffic. [NumPy documents its data-allocation tracking hooks](https://numpy.org/doc/2.4/reference/c-api/data_memory.html); private BLAS malloc/work buffers, fragmentation, resident shared libraries and system cache state can escape this accounting. Per-call probes subtract the returned product's bytes to isolate conversion overhead. No native malloc profiler was installed.

Benchmark child processes ran sequentially, with no parallel benchmark workers. Lightweight inspection commands ran during portions of the matrix; paired timing ran after it. Host load and thermals were not controlled. A 1.34-second focused test overlapped the start of the discarded preliminary run, not the validated dataset. File-backed trials use locally prepared synthetic fixtures and warm filesystem caches. Small source-free read/decode benchmarks are not end-to-end workflow timings. The fixtures preserve packed storage and use the production reader but omit full reference parsing, pair construction, genotype correlations, metadata preparation, aggregation and writers.

An initial untiled oracle exceeded `rtol=atol=1e-14` for two near-zero results from **current explicit64 production code**; the largest violating difference was 1.44e-14. This was a benchmark-comparison issue, not evidence against Boolean arithmetic. Production comparisons now use the same production tiling with explicit64 as the reference; independent `math.fsum` checks remain separate and unchanged. Preliminary runs are retained in `pilot.jsonl`, `results.jsonl`, and `oracle-diagnostic.txt`, and excluded from the final matrix. A later split-variant tolerance failure was recorded as a failed option and the remaining matrix resumed; no tolerance was relaxed. The first 210 final records predate explicit pass/fail fields and passed the same assertion before being saved.

## Measured memory and runtime

The [complete matrix](MATRIX.md) gives all cases and process ranges; [raw measurements](results-validated.jsonl) and [summary.json](summary.json) are machine-readable. The following selection covers actual reader/projection loops. Each cell is **explicit64 → native**. RSS and temporary allocations are different metrics and must not be subtracted from one another.

| Actual projection shape: contributors → outputs | Columns/source | Warm elapsed, ms | Peak RSS, MiB | Additional traced peak, MiB |
| --- | --- | ---: | ---: | ---: |
| CSR 4,096 → 17,380 | 1 binary | 1.76 → 1.59 | 136.73 → 133.78 | 8.520 → 8.524 |
| CSR 4,096 → 17,380 | 100 binary | 34.65 → 34.88 | 189.39 → 184.56 | 37.870 → 38.261 |
| CSR 4,096 → 17,380 | 1,000 binary | 542.53 → 704.49 | 318.91 → 344.92 | 48.551 → 50.551 |
| CSR 4,096 → 17,380 | 100 mixed | 35.34 → 36.84 | 193.27 → 199.09 | 37.870 → 39.433 |
| CSR 4,096 → 17,380 | 1,000 mixed | 558.69 → 541.45 | 364.81 → 361.70 | 48.551 → 56.551 |
| Dense 128 → 4,096 | 1,000 binary | 48.63 → 47.69 | 280.97 → 279.97 | 33.014 → 33.136 |
| Dense 4,096 → 128 | 1,000 binary | 50.24 → 55.85 | 301.59 → 288.34 | 33.356 → 19.992 |
| Dense 128 → 512, output-filtered | 1,000 binary | 17.05 → 18.84 | 263.41 → 263.97 | 8.829 → 8.951 |
| Dense 4,096 → 16, output-filtered | 1,000 binary | 58.05 → 46.22 | 283.47 → 269.11 | 33.356 → 18.283 |

The last case saves **14.36 MiB median measured RSS** and 15.07 MiB of traced temporaries. Explicit/native RSS ranges were 281.36–287.17 / 267.83–271.92 MiB; elapsed ranges were 47.73–87.74 / 42.58–46.47 ms. This is the strongest favorable local example. The reduction comes from avoiding the old-float64/new-float64/read-buffer overlap, not from removing the library's conversion.

The width-1,000 binary CSR case instead increases traced temporary peak by 2 MiB. Its 318.91/344.92-MiB RSS medians have broad, overlapping ranges: 303.39–343.59 / 282.70–354.97 MiB. Elapsed ranges likewise overlap: 497.12–1,162.43 / 534.43–960.37 ms. Treat its 30% median slowdown and 26-MiB RSS increase as noisy observations, not portable effect estimates. The stable allocation result is the extra retained Boolean tile during multiplication. Mixed CSR at the same width adds 8 MiB of traced temporary peak because the retained source operand is float32.

Memory accounting for that binary CSR example separates the persistent 132.60-MiB accumulator and 8.22-MiB original CSR operator from the measured extra workspace. At a full tile, explicit multiplication holds about 16 MiB annotation + 16 MiB accumulator gather + 16 MiB product, plus the sliced CSR and indices. Native multiplication holds 2 MiB Boolean + 16 MiB conversion + the same two 16-MiB output temporaries. This explains the approximately 48.55 versus 50.55 MiB traces. The 16-MiB budget never covered all those arrays.

To isolate repeated casts from process-to-process drift, [cast_cost.py](cast_cost.py) additionally ran 25 paired, randomized warm trials per operator/width/reuse combination with both inputs already resident. These timings exclude the explicit once-per-tile preparation cost and therefore describe multiplication cost, not the complete projection:

| Width | Operator | Four products: prepared float64, ms | Four products: Boolean, ms | Ratio |
| ---: | --- | ---: | ---: | ---: |
| 1 | Dense | 14.42 | 21.60 | 1.50 |
| 1 | CSR | 1.49 | 7.90 | 5.30 |
| 100 | Dense | 15.64 | 23.65 | 1.51 |
| 100 | CSR | 5.78 | 11.92 | 2.06 |
| 1,000 | Dense | 8.68 | 15.50 | 1.79 |
| 1,000 | CSR | 3.03 | 9.23 | 3.04 |

Both libraries repeat conversion for every output block. In the complete synthetic read-plus-four-products benchmark at width 1,000 and 10% binary density, dense elapsed was 16.58 → 18.03 ms and CSR 9.88 → 20.08 ms. One-product cases were much closer. Low operator work makes conversion comparatively costly; at 10% CSR operator density, arithmetic dominates and the separate matrix shows much smaller relative differences. Binary annotation density does not shrink the operand, conversion, or output; the library kernels still process converted zero entries. Runtime fluctuations across annotation densities do not establish a density-based optimization.

Mixed-group splitting sometimes lowered traced memory and improved synthetic CSR timing. For 1,000 columns/10% continuous/four CSR products, explicit64/native/split took 15.68/18.36/12.83 ms, with traced peaks 25.80/25.96/18.76 MiB. At 50% continuous the corresponding timings were 12.20/13.44/14.84 ms. Dense splitting was not consistently faster and failed a numerical gate below. These are insufficient grounds for a mixed-selection implementation.

## Numerical evidence

All **243 native-operand process measurements** matched their explicit64 counterparts exactly: maximum absolute difference zero, with float64 results and accumulation. This includes Boolean-only operands and mixed float32 operands with float64 operators, unchanged contributor/output tiling, negative corrected R², and logical interleaving of binary/continuous columns.

The independent [numerical checks](numerical.json) compare 36 combinations of width, dense/CSR backend, PLINK-like float64 correction or Parquet-like float32-normalized correction, and explicit/native/split operands against selected-column `math.fsum` dot products. All passed existing `rtol=1e-14, atol=1e-14`; maximum absolute difference was 1.42e-14, within the combined relative-plus-absolute threshold. Six additional production checks cover diagonals, signed annotations, reordered output rows, non-output contributors, weight-mask semantics and pair chunks. The focused unchanged suites passed **45 tests in 1.34 seconds**; see [focused-tests.txt](focused-tests.txt). No full suite was necessary because production behavior was not changed.

The large mixed dense split case, 20,971 contributors × 100 columns with 10% continuous, failed the unchanged tolerance at 29 of 12,800 output cells, consistently across three processes. Maximum absolute difference over all cells was **7.96e-12**. A [separate diagnostic](mixed-split-diagnostic.json) found native grouped products and explicitly float64 grouped products identical, while both differed from the full-width product. Thus grouping changed the floating-point reduction behavior even with double arithmetic; the failure is not a float32-accumulation regression. This option fails the requested criterion as tested.

The current tiled explicit implementation differed from an untiled reference by up to 2.13e-13 over the large projection experiments. This is recorded independently from the zero candidate-versus-current differences. The adversarial negative control used a CSR row with 100,000 copies of 0.1: float32 operator × Boolean lost 1.44336 against `math.fsum`, versus 1.88e-8 for float64. It illustrates the dtype hazard; it is not a reproduction of the historical ~3e-3 dataset-specific discrepancy.

Artifact verification compiled all five experiment scripts, checked local Markdown links and update dates, verified exactly three replicates for all 170 variant groups with no duplicate measurements, and confirmed that the 28 pre-existing working-tree status entries were preserved. `git diff --check` passed. See [validation.json](validation.json).

## Whole-job significance and options

The nominal 14 MiB × 22 = **308 MiB** is a simultaneous operand-size calculation, not measured RSS savings. Hidden conversions eliminate that simple saving during multiplication; previous/new-tile lifetimes can create a different saving in some dense calls. The 22-worker multiplier also assumes every worker reaches the relevant peak simultaneously. None of those assumptions was measured here.

At the current default width 100, the local chr22 BIM's largest left window, 8,192 contributors, fits one annotation tile; its float64 operand is 6.25 MiB, and the median window is only 1.56 MiB. The nominal 16-MiB operand is not routinely reached by this PLINK workload. Width 1,000 can span up to four contributor tiles in the same proxy, making old/new-tile overlap relevant. Forward PLINK projection reads only 128 contributors, about 0.98 MiB at width 1,000.

For context, the historical [sequential-query benchmark, “Resource measurements and limits”](../annotation-memory/sequential-query-batches.md#resource-measurements-and-limits) measured about 1.166 GiB peak on chr22, before these experiments and before all current storage changes. A hypothetical 14-MiB reduction would be about 1.2% of that peak **only if that projection phase determined the peak**. Its 32–48-GiB whole-job planning range would make 308 MiB about 0.6–0.9%; a 70-GiB run would make it about 0.43%. These are scale comparisons, not new job estimates. The reported 70-GB job's accounting scope remains unverified.

Even the dimension-only float64 accumulator term is 139.76 MiB for 17,380 output rows × (1,000 queries + 53 baselines + weight), and 8.91 GiB when summing 1,134,471 output rows at that width across chromosomes. These are array-size calculations, not simultaneous measured RSS. Genotypes, pair buffers, CSR construction, output conversion and aggregation can dominate a different phase. If the whole-job maximum occurs there, the projection operand change may save **zero** peak job RSS. Aggregation and batch retention were neither optimized nor benchmarked here.

Options for discussion:

1. **Keep the explicit float64 projection path.** The evidence supports this as the general default: predictable precision across operator dtypes, one conversion per annotation tile, and no demonstrated broad peak-memory benefit from native operands.
2. **If wide, output-filtered PLINK jobs matter enough, examine that specific phase next.** The 14.36-MiB local RSS reduction is a credible narrow signal, but a matched whole-job/phase profile is needed to establish practical value and runtime. Retaining explicit conversion while shortening tile lifetimes is another hypothesis to compare; it has not been implemented or measured here.
3. **A kernel that actually consumes Boolean values without making a float64 tile would be a different research project.** It could change the memory conclusion, but neither installed NumPy nor this SciPy path supplies it. Mixed-group splitting already failed one existing tolerance, so it is not ready as a substitute.

No production strategy has been selected. In particular, there is no evidence for promising a 308-MiB job saving, and no justification here for changing public reads, packed storage, scientific tolerances, output formats, aggregation or batch retention.

## Reproduction

Use the same editable environment and current production source. [environment.json](environment.json), [initial-working-tree.txt](initial-working-tree.txt) and [preexisting-projection.patch](preexisting-projection.patch) identify the measured state; the patch records the user's pre-existing default change from 1,000 to 100. The benchmark explicitly selects each tested width, so that default change does not confound comparisons.

From the repository root:

```bash
/Users/wenbinwu/miniforge3/envs/ldsc3-dev/bin/python docs/audits/2026-09-16-boolean-projection/benchmark.py --prepare
/Users/wenbinwu/miniforge3/envs/ldsc3-dev/bin/python docs/audits/2026-09-16-boolean-projection/benchmark.py --numerical > docs/audits/2026-09-16-boolean-projection/numerical.json
/Users/wenbinwu/miniforge3/envs/ldsc3-dev/bin/python docs/audits/2026-09-16-boolean-projection/benchmark.py --replicates 3 --output results-validated.jsonl
/Users/wenbinwu/miniforge3/envs/ldsc3-dev/bin/python docs/audits/2026-09-16-boolean-projection/cast_cost.py
/Users/wenbinwu/miniforge3/envs/ldsc3-dev/bin/python docs/audits/2026-09-16-boolean-projection/context_probe.py
/Users/wenbinwu/miniforge3/envs/ldsc3-dev/bin/python docs/audits/2026-09-16-boolean-projection/mixed_split_diagnostic.py
/Users/wenbinwu/miniforge3/envs/ldsc3-dev/bin/python docs/audits/2026-09-16-boolean-projection/summarize.py
/Users/wenbinwu/miniforge3/envs/ldsc3-dev/bin/python -m pytest -q tests/test_ldscore_projection.py tests/test_annotation_storage.py tests/test_annotation_format.py
```

`--prepare` recreates only this audit's synthetic fixtures; their pickle descriptors contain local paths and should be regenerated after moving the checkout. To preserve this run, copy the audit directory or choose a different `--output` filename before rerunning; the summary script reads `results-validated.jsonl`. `--resume` appends only missing case/variant/replicate combinations after an interrupted run. No network access or additional library installation is needed for reproduction.
