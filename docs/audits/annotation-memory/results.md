# Annotation memory: local resource and equivalence results

Last updated on: 2026-09-10

The completed refactor reduces avoidable memory for large pathway batches while keeping public HM3 LD tables aggregate. These local synthetic measurements show the trade-off: the 1,000-query direct and regression runs use less memory, while staging and smaller computations can cost runtime. This is evidence for the implemented ownership and batching design, not a production runtime guarantee. See the [developer design](../../current/annotation-memory-design.md), [plan](../../plans/2026-09-10-annotation-workflow-memory.md), [raw measurements](measurements.json), and [artifact comparisons](equivalence.json).

## Matched command measurements

The baseline is post-fingerprint-removal revision `a505c45d61d35b39a3fc086c2b623dd7d6936d96`. Direct, regression, and quantile measurements use `85a0466`; the final standalone and indexed measurements use `91ad7bf`, which fixes standalone scratch lifetime and reduces staged gene bookkeeping overhead. Later changes to those measured paths are documentation/type annotations only. Every record identifies the checkout, revision, imported package path, dimensions, and measurement result. Both versions use the same prepared inputs and command options; the new batching flag is applied only to the refactored version.

The table uses one chromosome worker and the approved default query batch size of 1000 where applicable. MiB means 1,048,576 bytes. Peak private disk and final persistent output size are separate from RSS and from each other.

| Workload | Old peak RSS MiB | New peak RSS MiB | Old elapsed s | New elapsed s | Old peak private MiB | New peak private MiB | New persistent MiB |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Direct R2, 18,000 SNPs, 1,000 queries | 876.2 | 560.1 | 4.86 | 5.65 | 0 | 140.17 | 6.20 |
| Direct PLINK, 18,000 SNPs, 1,000 queries | 934.7 | 533.4 | 4.75 | 10.94 | 0 | 126.43 | 19.73 |
| Indexed assembly, 1,000 gene lists | 389.5 | 346.9 | 7.58 | 27.86 | 0 | 1.59 | 2.49 |
| Partitioned regression, 1,000 separate models | 335.2 | 258.5 | 19.25 | 14.25 | 33.04 | 35.60 | 34.49 |
| Exact quantiles, 120,000 SNPs, 60 annotations | 730.6 | 460.9 | 3.03 | 4.89 | 0 | 58.02 | 0.02 |
| Standalone BED annotation, 100 queries | 150.5 | 150.1 | 14.04 | 1.45 | 0.91 | 7.44 | 0.16 |

Direct R2 and PLINK peak RSS fell about 36% and 43%, respectively, on the 1,000-query fixture. Regression fell about 23% at the default batch size and about 44% with batch size 32. Quantile reconstruction fell about 37%. The small standalone fixture is dominated by runtime/import overhead and shows essentially unchanged RSS. The new standalone gene route used 149.3 MiB RSS, 3.38 s, and 7.38 MiB private disk for the same 100 interval memberships. The old revision has no standalone gene command: gene-versus-BED is a membership-equivalence comparison, not a same-command speed comparison.

All final successful runs left **zero private scratch bytes**. Initial standalone measurements exposed a returned bundle that still depended on preparation scratch; `38bed04` replaced that return with persistent query descriptors plus deferred baseline access. Those initial records are marked superseded in the raw measurements. Persistent file sizes can differ slightly because provenance records the new settings and schemas; the scientific artifact comparisons below check contents directly.

## Batch size, workers, chromosomes, and query count

| Refactored workload | Queries | Chromosomes | Batch | Workers | Peak RSS MiB | Elapsed s |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| R2 | 1,000 | 3 | 1000 | 1 | 560.1 | 5.65 |
| R2 | 1,000 | 3 | 32 | 1 | 551.4 | 4.60 |
| R2 | 1,000 | 3 | 32 | 2 | 820.1 | 4.86 |
| R2 | 1,000 | 1 | 32 | 1 | 478.5 | 2.25 |
| R2 | 100 | 3 | 1000 | 1 | 335.0 | 1.78 |
| PLINK | 1,000 | 3 | 1000 | 1 | 533.4 | 10.94 |
| PLINK | 1,000 | 3 | 32 | 1 | 566.0 | 12.29 |
| PLINK | 1,000 | 3 | 32 | 2 | 907.5 | 9.93 |
| PLINK | 1,000 | 1 | 32 | 1 | 494.0 | 4.76 |
| Regression | 1,000 | 2 | 1000 | 1 | 258.5 | 14.25 |
| Regression | 1,000 | 2 | 32 | 1 | 188.7 | 14.21 |
| Indexed assembly | 1,000 | 3 | 1000 | 1 | 346.9 | 27.86 |
| Indexed assembly | 1,000 | 3 | 32 | 1 | 345.7 | 30.40 |

The one-versus-three-chromosome direct runs keep 6,000 SNPs per chromosome. Old one-chromosome peaks were 695.6 MiB for R2 and 770.5 MiB for PLINK. Final output rows grow from 1,200 to 3,600, so some growth in retained HM3 output is intentional. The 100-query R2 input selects the first 100 columns of the same 1,000-query fixture, preserving the reference, baseline, row count, and output selection. Its old peak was 326.1 MiB and elapsed time 1.82 s: staging overhead can offset the memory benefit at this smaller width.

At fixed output dimensions, a smaller batch reduces the active query projection or regression load, but total measured RSS need not decrease monotonically in short runs. The allowed output accumulator, parsing, allocator reuse, and sampling timing also contribute. Two workers retain two active chromosome payloads and additional process runtimes; these measurements show the expected increase in process-tree RSS. No old two-worker run was measured, so the worker rows compare the refactored configuration with itself, not old versus new parallel speed.

The indexed run remains slower on this small sparse-operator fixture. Profiling identified per-source resolution/summary construction and staged selection conversion as substantial costs. `91ad7bf` removes repeated summary groupby construction and full selection-frame copies without caching every source. Initial default-batch runtime was 42.39 s; two subsequent final default observations were 30.82 s and 27.86 s. Batch 32 took 30.40 s. This does not establish a statistically significant speed difference between the final batch settings; the approved default remains 1000.

The instrumented index comparison separates validation from sparse multiplication: old/new validation took 0.150/0.097 s and multiplication took 0.039/0.040 s. The remaining workflow time was 6.52/27.10 s, including resolution, staging, support checks, output assembly, and writing. Holding the operator through each chromosome's batches avoids repeated operator reads, but this fixture is too small to make multiplication the dominant cost. The profile itself is excluded from the resource table because profiling overhead changes runtime.

## Scientific and artifact equivalence

[`compare_annotation_memory.py`](../../../benchmarks/compare_annotation_memory.py), `compare()` and `compare_table()`, compares canonical tables, identities and order, annotation/count/overlap metadata, gene audit/status tables, and per-query regression metadata. Numeric comparisons use existing `rtol=1e-6`, `atol=1e-8`; identity/category/delete-block columns are exact. No reference goldens or tolerances were changed.

All measured current variants passed. Direct LD values, annotation memberships, counts, overlaps, and all 1,000 separate regression fits matched exactly in these fixtures. Each regression variant checks 3,002 tables, including category results and coefficient delete values, plus fitted-model metadata. The final standalone gene annotations match the independently specified BED intervals exactly. Quantile output differs by at most **1.5210055437364645e-13**, consistent with the changed stable accumulation order and well within the existing tolerance. The broader test suite separately covers membership boundaries, ties, omitted target values, both baseline layouts, identity cleanup, skipped siblings, gate failures, overwrite behavior, reloadability, bounded reads, and object lifetimes; see [progress](progress.md).

These checks compare the same backend before and after refactoring. The synthetic R2 fixture is not derived from the random PLINK genotypes, so comparing those two backends' benchmark values to each other would not be a valid scientific test. Independent backend and direct/indexed numerical oracles remain in the regular tests.

## Inputs, environment, and measurement method

[`annotation_memory.py`](../../../benchmarks/annotation_memory.py), `prepare()`, fixes seed **20260910**. Direct fixtures contain three chromosomes, 6,000 SNPs each, two baseline columns, and 1,000 binary query columns. The PLINK panel has 128 samples; the R2 panel uses deterministic short-range pairs. Whole-genome baseline and query gzip files are separate. Output selection retains every fifth SNP. The narrower input contains the first 100 query columns. Standalone projection uses 100 BEDs or matching lists from a 3,000-gene catalog. The index was built once from the same PLINK inputs and reused for both assembly versions; its three-chromosome fixture uses the existing private partial-index test seam. Production index commands still require autosomes 1–22. Index construction time is excluded from the assembly comparison.

The regression fixture contains 2,000 LD rows, 1,000 query columns, and 20 jackknife blocks. Every model is a separate baseline-plus-query fit over its full retained genome-wide SNP set. The quantile fixture has 120,000 rows, 60 fitted binary annotations, 10 exact global quantiles, and 20 coefficient delete vectors; fitted artifacts are deterministic synthetic valid inputs to the quantile post-processing command. Fixture creation is excluded from every measured command.

All measurements ran locally on macOS 15.3.2, arm64, using the existing `ldsc3-dev` environment: Python 3.13.13, NumPy 2.4.5, pandas 2.3.3, SciPy 1.17.1, and PyArrow 23.0.1. `measure()` starts a fresh worker with the selected checkout's `src` on `PYTHONPATH` and verifies the imported package path. BLAS/OpenMP thread limits are set to one in both versions; chromosome workers are varied independently. Commands run sequentially, include imports, preflight, computation, and writing, and use `TMPDIR` inside their own output directory for library scratch.

`tree_rss()` samples only PID, parent PID, and resident size using `ps`, then sums the worker and its descendants at that instant. It never sums separate child lifetime peaks as though simultaneous. Sampling sleeps 50 ms between observations, in addition to process/disk inspection overhead. `disk_bytes()` scans only the owned output tree and sums logical file sizes in private dot directories and `library-tmp`; it records public output bytes separately. The external captured console log and harness `worker.json` are excluded from persistent output size. A sample can miss a short peak, and scanning overhead varies with the number of output files. File size is not physical allocated disk usage.

Most cells are single observations with no cold-cache guarantee, confidence interval, or production-scale extrapolation. Runtime differences of a few tenths of a second should not be overinterpreted. Ownership tests, accessor/allocation guards, and numerical tests establish the retention contract; these measurements show practical consequences for the stated local inputs. No HPC work was performed.

## Reproduction

Run from the repository root with the development environment active. The baseline checkout is a separate local worktree; do not replace an existing unrelated directory. Generation and index construction happen once, outside timed commands. Output paths must be new because the harness refuses to overwrite an existing measurement.

```bash
git worktree add --detach .worktrees/memory-baseline a505c45
python benchmarks/annotation_memory.py prepare --inputs "$PWD/.worktrees/memory-benchmark/inputs"
python benchmarks/annotation_memory.py build-index --inputs "$PWD/.worktrees/memory-benchmark/inputs"
python benchmarks/annotation_memory.py matrix \
  --checkout "$PWD" \
  --inputs "$PWD/.worktrees/memory-benchmark/inputs" \
  --output "$PWD/.worktrees/memory-benchmark/runs"
python benchmarks/annotation_memory.py prepare-direct-subset \
  --inputs "$PWD/.worktrees/memory-benchmark/inputs" \
  --output "$PWD/.worktrees/memory-benchmark/inputs-100" --query-count 100
python benchmarks/annotation_memory.py measure \
  --checkout "$PWD/.worktrees/memory-baseline" \
  --inputs "$PWD/.worktrees/memory-benchmark/inputs-100" \
  --output "$PWD/.worktrees/memory-benchmark/runs/baseline-direct-r2-q100" \
  --case direct-r2-q100
python benchmarks/annotation_memory.py measure \
  --checkout "$PWD" \
  --inputs "$PWD/.worktrees/memory-benchmark/inputs-100" \
  --output "$PWD/.worktrees/memory-benchmark/runs/current-direct-r2-q100" \
  --case direct-r2-q100 --batch 1000
python benchmarks/compare_annotation_memory.py \
  --runs .worktrees/memory-benchmark/runs \
  --output docs/audits/annotation-memory/equivalence.json
```

The harness's `measure` mode also accepts `--batch` and `--workers` for individual configurations. On a restricted local host, process RSS inspection needs permission to run `ps`; it reads process sizes and parent IDs only. The committed JSON records retain superseded observations so the tuning and ownership correction remain distinguishable from final measurements.
