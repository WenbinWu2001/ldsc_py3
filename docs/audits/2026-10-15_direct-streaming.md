# Direct LD-score streaming integration

Last updated on: 2026-10-15

## Implemented boundary

This change integrates the approved direct file-writing design on `restructure`, starting from `848b103`.
[`stage_direct_chromosome()`](../../src/ldsc/_direct_ldscore_batches.py) stages existing post-kernel float32 query scores and returns baseline data with counts, overlaps, variation summaries, and diagnostic descriptors.
[`LDScoreCalculator`](../../src/ldsc/ldscore_calculator.py) combines chromosome statistics, applies global query filtering, and yields the existing `StagedLDScoreBatch` to the shared output writer.
The parent does not transport or assemble a genome-wide query matrix.
The materialized no-output API and the direct numerical kernels remain unchanged.

[`_query_variation.py`](../../src/ldsc/_query_variation.py) is shared with indexed scoring and preserves each producer's precision.
Direct filtering uses the existing float32 scores; indexed filtering still sees float64 scores before public float32 conversion.
The helper retains exact global variation for within-chromosome constants that differ between chromosomes, empty chromosomes, NaNs, infinities, signed zero, and differences below float32 precision.
Counts, overlaps, numerical accumulation, reference contributors, and separately selected regression/output rows keep their existing contracts.

The persistent pool spans all batches and publication.
Each active batch owns its fragment directory, and the shared writer validates chromosome/SNP order, releases each consumed frame, and deletes its fragment before reading the next.
Skipped batches also release their fragments.
Worker shutdown precedes removal of borrowed annotation/fragment workspaces on failure.
Scientific outputs remain private until every batch succeeds, with metadata last under the existing publication contract.
There is no new CLI option, output format, checkpoint, or resume behavior.

## Local validation

Validation used the existing editable `ldsc3-dev` environment: Python 3.13.13, NumPy 2.4.5, pandas 2.3.3, and PyArrow 23.0.1.
The first real-kernel transport test failed on the original implementation with `worker transported a query matrix` before implementation began.
It passed after integration.

[`tests/test_direct_batch_staging.py`](../../tests/test_direct_batch_staging.py) exercises these behaviors:

- Real two-chromosome parquet-R2 calculation with two workers and multiple batches matches materialized writing at the same batch width exactly for saved scores, baseline values, counts, overlaps, statuses, annotation types, SNP-universe policy, and ordered row groups; saved scores also match the small no-output API exactly.
- Batch widths 1, 2, and 10 cover zero support, global constants, chromosome constants that vary globally, skipped batches, and an empty first output chromosome while preserving full-reference counts.
- Parent-frame weak references and file-existence checks establish release and deletion before the next chromosome fragment is read.
- Final-batch worker, row-alignment, and writer failures preserve diagnostics, remove owned scratch, and leave no published scientific outputs or metadata.
- A controlled concurrent failure verifies that borrowed annotation and fragment directories remain present until worker shutdown.
- A real small PLINK gene-list workflow compares streamed and materialized writing, including exact tables/schemas, metadata, gene audits, identity-drop reports, chromosome scope, and exported reference metadata.

[`tests/test_indexed_batch_staging.py`](../../tests/test_indexed_batch_staging.py) retains the indexed lifecycle checks and checks shared variation in both float32 and float64 against pandas distinct-value semantics.
The unchanged no-output, baseline-only, overwrite, SNP-universe, and numerical-kernel paths are also covered by the existing package suite.

The full pytest run passed **1,991 tests and 133 subtests**, with one skip, in 412.62 seconds.
The skip is the existing missing-PyArrow dependency case, which skips because PyArrow is installed.
The run reported 209 warnings; these were not suppressed or reclassified by this change.
The separate unittest compatibility run passed **1,003 tests**, with one skip, in 74.94 seconds.

A final strict overlap assertion initially compared a single batch of three queries with batches of two and one.
The continuous query's all-reference diagonal differed by one float64 step: `5.8799999919533725` versus `5.879999991953373`.
Materialized writing at the same two-query batch width matched streaming exactly, including both overlap diagonals.
The final parity test therefore controls batch width and worker count when isolating the assembly change; the kernel's existing cross-batch numerical tolerances are unchanged.
After strengthening those assertions, the final direct/indexed streaming, batch, and lifecycle selection passed **67 tests** in 41.80 seconds.

```bash
PYTHONDONTWRITEBYTECODE=1 /Users/wenbinwu/miniforge3/envs/ldsc3-dev/bin/python3 -m pytest -q --tb=short
PYTHONDONTWRITEBYTECODE=1 /Users/wenbinwu/miniforge3/envs/ldsc3-dev/bin/python3 -m unittest discover -s tests -p 'test*.py' -v
PYTHONDONTWRITEBYTECODE=1 /Users/wenbinwu/miniforge3/envs/ldsc3-dev/bin/python3 -m pytest tests/test_direct_batch_staging.py tests/test_indexed_batch_staging.py tests/test_ldscore_query_batches.py tests/test_ldscore_streaming.py tests/test_gene_index_streaming.py -q --tb=short
```

Documentation validation checked 360 local links, 100 heading anchors, and 29 complete CLI examples without errors.
Both edited notebooks retain all 16 code cells, outputs, cell metadata, and notebook metadata unchanged; two explanatory Markdown cells changed.
Code cells and all eight changed Python files pass syntax checks.
Notebook execution was not repeated for these prose-only edits, and no documentation builder is configured.
`git diff --check` passed.

## Performance evidence and limits

The supplied frozen replay evidence remains unchanged under the workspace's [direct-streaming experiment directory](../../../diagnostics/20261015_direct_streaming_scale/report.md).
Its [lifecycle report](../../../diagnostics/20261015_direct_streaming_scale/lifecycle.md) records persistent workers and complete-success publication.
The [user-facing measurement table](../wiki/main-functionalities/ldscore.md#direct-assembly-measurements) summarizes the reported full-row medians and preserves node groups and measurement limits.
Those replays used precomputed scores and the frozen `03aa8ca` experiment baseline; they did not execute genotype/reference preparation or LD computation.
Their roughly 41–46% lower cgroup peaks and lower observed assembly wall/CPU times support this design but do not measure the integrated production pipeline.
The 257 MiB active-fragment figure is a total across 22 files, and sampled workload-file peaks exclude fixed inputs and measurement logs.

No production-scale memory, file-cache, scratch, CPU, or wall-time measurement was run for this integration.
No remote package, original results, frozen sources, or benchmark evidence was modified, and no new HPC job was submitted.
A real two-batch comparison remains subject to separate approval and a narrow exception for restricted input paths, as recorded in the supplied [end-to-end proposal](../../../diagnostics/20261015_direct_streaming_scale/end_to_end_proposal.md).
Do not reduce the production Slurm memory allocation based on these local tests or the assembly replays alone.
