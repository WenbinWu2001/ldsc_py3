# Parallel execution decisions

Last updated on: 2026-10-27

This record consolidates the user's finalized parallelization annotations and subsequent `workers` / `joblib.cpu_count()` decisions.
The implementation now follows these decisions; the maintained runtime contract is [parallel-execution.md](../reference/parallel-execution.md).
The earlier inline, inherited-BLAS, and `threads` interfaces have been replaced.
The final numerical policy is one BLAS thread per worker for every calculation, whether execution is serial or parallel.
This supersedes the earlier allowance for serial multithreading and user-selected higher BLAS limits.
Code and documentation are written in English.

## Reconciliation of the supplied decisions

| Statements to reconcile | Governing decision or interpretation |
| --- | --- |
| Earlier decisions allow multithreaded serial BLAS and cap it at the detected CPU count | Superseded by the final uniform `B=1` policy; no serial multithreading option remains. |
| Annotation 2 requires stopping and deleting temporary results; annotation 4 item 6 says to continue and retain successful work | Use the dedicated chromosome-failure policy in annotation 2 as the governing interpretation; the continuation/retention sentence is incompatible with it and should be discarded. |
| Earlier numerical limits differ between serial and parallel execution | Use one numerical policy at every worker count; retain the separate initial Arrow CPU/I/O defaults. |
| The fixed numerical policy requires verification; the follow-up answer permits continuing without it | Always request the single-threaded limit, but continue with a clear diagnostic when enforcement or verification is unavailable; this is a fixed policy with best-effort backend enforcement, not an unconditional runtime guarantee. |
| Earlier decisions retain `threads`; the latest request renames it | The target name is `workers` in Python/configuration and `--workers` in the CLI. |
| Annotation 6 retains current logging; annotation 8 specifies chromosome logs | Apply annotation 8 to direct and indexed LD-score logging; it does not require a new centralized worker-log transport. |
| Annotation 2 limits implementation to its items 3 and 4 | Interpret this as the scope of the failure-handling change: pending-task cancellation and reference-metadata staging, without a broader failure-system redesign. Other separately approved design topics remain in scope for the overall target. |
| Annotation 7 defers serial-phase and repeated-work improvements | Keep them in later benchmarking work; do not reorganize computation as part of this policy change. |

The failure-policy and implementation-scope rows state the interpretation used to reconcile conflicting text, rather than silently treating both alternatives as compatible.
The chromosome-failure policy does not remove regression's existing explicit `continue_on_query_error` behavior.

## User interface and CPU budget

Keep one routine concurrency setting, defaulting to `workers=1` / `--workers 1` on laptops and HPC alike.
Follow the repository's default against compatibility aliases: migrate the existing Python arguments, configuration fields, CLI flags, messages, tests, and examples together, without adding a second `threads` spelling.
This is a rename of the existing setting, not a new resource-control option.

Use `C = joblib.cpu_count()` with its default `only_physical_cores=False`.
Describe `C` as **the detected available CPU count**.
The function accounts for supported restrictions such as CPU affinity and Linux CPU quotas; it does not measure currently idle CPUs or read the user's requested scheduler allocation.
Restrictions reflected in affinity or quotas can affect the result; an allocation not exposed through those mechanisms cannot be inferred.
See the [joblib API documentation](https://joblib.readthedocs.io/en/stable/generated/joblib.cpu_count.html).

Use joblib's returned count directly as the worker budget, without a custom detection algorithm, Slurm-specific parsing, or custom CPU binding.
Record the detection source as `joblib.cpu_count()`; do not invent an attribution to affinity, quota, or a scheduler when the function only returns a count.
The user must choose concurrency consistent with any allocation that the execution environment does not expose.
Matching `--workers` to such an allocation bounds task concurrency; it does not impose a whole-job CPU or software-thread ceiling.
`C` determines the compute-worker budget only; it no longer selects or caps a configurable BLAS thread count, because the numerical target is always one.
Adding joblib for this function does not require replacing existing execution machinery with `joblib.Parallel`.

For a nonempty workload, let `R` be the requested setting and `K` the available task concurrency.
These are explanatory symbols, not additional user parameters.

| Request | Resolution before the task-count cap |
| --- | --- |
| Omitted or `1` | One worker |
| Positive integer `N` | `N` workers, provided `N <= C` |
| `-1` | `C` workers |
| `-2` | `max(1, C - 1)` workers |
| Negative integer `-k` | `max(1, C + 1 - k)` workers |
| `0`, a Boolean, or a non-integer Python value | Invalid |

Reject an explicit positive request above `C` before computation, then cap valid requests by available work: `W = min(resolved_request, K)`.
For example, `C=8`, two chromosomes, and `--workers 4` gives `W=2`; `C=2`, one chromosome, and `--workers 4` is an error before the task cap.
For batch regression, retain the existing query-batch-size cap on simultaneously available models.
Existing empty-workload validation and no-work behavior remain workflow concerns; do not start an idle child merely to satisfy a count.

## Tasks, processes, and scheduling

| Workflow | Task unit | Target scheduling |
| --- | --- | --- |
| Direct LD scoring | One chromosome within the active query batch | Reusable child processes; preserve the existing sequential query-batch boundaries |
| Indexed LD scoring | One chromosome, including its requested query batches | Reusable child processes; preserve existing genome-wide assembly boundaries |
| Gene-index construction | One chromosome | Managed child processes after the focused memory/transport benchmark |
| Batch regression | One complete query-plus-baseline model across all selected chromosomes | Reusable child processes; preserve read-only numerical memory maps |

An available worker receives the next task and executes one task at a time.
A chromosome or query does not permanently own a process.
Use one managed child even when `W=1`, reusing it for successive tasks rather than starting a process per task.
The target is a persistent pool for the workflow's task execution, including successive batches, without nested chromosome/model worker pools.
Keep the existing single-machine execution model and managed process lifecycle; do not build a distributed scheduler or a general executor framework.
The operating system schedules processes; the package does not reserve or pin a logical CPU to each worker.

The parent remains the coordinator for preparation, submission, collection, aggregation, and publication.
Each child has its own main execution thread and may also have native-library or runtime support threads.
Parent/support processes and threads are not part of the requested compute-worker count.
One worker therefore means neither one total operating-system thread nor one dedicated logical CPU, and worker count does not guarantee utilization or a whole-job CPU ceiling.
A physical core, a logical CPU exposed by the operating system, and a software thread are distinct concepts; do not label `C` as a physical-core count.

Keep task inputs and staged outputs separate, pass compact input descriptions where practical, and avoid retaining unnecessary duplicate scientific arrays in parent and children.
Preserve regression's suitable read-only memory maps rather than introducing a second sharing mechanism.
Process transport can still create copies; do not promise memory safety or memory-based automatic worker selection without credible task-memory estimates.
The user accepts the local pilot's approximately 0.8–1.0-second child startup cost.
Its runtime-only control added approximately 158 MB of sampled RSS and 108 MB of footprint accounting; these are distinct memory measurements, with further input duplication possible.
These are [local pilot observations](../../../pilot_runs/20261012_parallel_policy_local/REPORT.md), not measurements repeated in this design update or guarantees for other workloads.

## Numerical and Arrow thread policy

Let `B` denote the BLAS thread limit per compute worker.
Always use `B=1`, for both `W=1` and `W>1`, across direct LD scoring, indexed LD scoring, gene-index construction, and regression.
The single-threaded BLAS policy also applies to package-owned calculations outside chromosome/model task loops; it is not a parallel-pool-only setting.
It does not introduce additional public worker controls or change the separately agreed scope of process migration.
The target is fixed from worker startup through completion, including after task-count capping or when only one task remains.
Do not select a higher BLAS count from `C`, platform, workflow, task size, or remaining workers.
Higher inherited BLAS settings do not override this package policy inside executing workers, and no serial BLAS override or new user parameter is introduced.

| Component | Policy for all worker counts |
| --- | --- |
| BLAS | One thread per worker |
| OpenMP numerical computation | Keep a separate one-thread limit to prevent nested numerical parallelism; this control is distinct from `B` |
| Apple Accelerate | Request and verify supported single-threaded BLAS/LAPACK mode; automatic mode is not a selected execution policy |
| Arrow CPU pool | Initially 1 per worker |
| Arrow I/O pool | Initially 1 per worker |

Configure compatible numerical controls inside the executing child and verify all supported loaded numerical backends, including separate NumPy and SciPy BLAS libraries where present.
Use the existing threadpoolctl/backend support where applicable; environment variables alone do not establish that a loaded library obeyed a limit.
If the numerical limit cannot be enforced or verified, continue with the resolved worker count and clearly record the requested limit, any observed setting, and why the limit remains unverified or unmet.
Do not report a verified one-thread limit in that case; this is an explicitly accepted best-effort resource policy.
This exception concerns numerical-thread control, not computation failures, which retain their existing failure policy.
Do not change the notebook or Python caller's global numerical or Arrow settings.
When implementing the universal policy, account for numerical preparation and finalization outside task loops rather than assuming a worker initializer covers all calculation.
Keep Arrow CPU and I/O controls independent internally, without extra routine CLI options.
These are limits on individual library pools, not a claim that their sum equals total threads or simultaneous CPU demand.
Validate Arrow's initial defaults on representative HPC storage before treating their performance as established.

The [Longleaf pilot](../../../pilot_runs/20261012_parallel_policy_longleaf/REPORT.md) supports single-threaded regression and offers no convincing indexed-LD benefit from higher BLAS limits on its compact workload.
The direct-LD serial comparison does not establish the fastest setting, but the user has chosen the uniform one-thread policy for simplicity and resource efficiency.
Further serial BLAS tuning is not a prerequisite for this decision.
The reports' inline-serial recommendations do not supersede the accepted managed-child design, and their workload-specific four-worker finding does not change the default `W=1`.

## Chromosome failures, staging, and publication

Any failure of a required chromosome calculation fails the run.
After correcting the problem, the user reruns every chromosome requested for that run.
There is no partial-success, checkpoint, resume, successful-shard reuse, or automatic-retry workflow.
Preserve existing intentional skip/eligibility rules; these decisions concern failed required calculations, not a new definition of scientific eligibility.

On an observed task failure, stop submitting new tasks, attempt cancellation of work the executor has not started, allow the existing executor shutdown to finish, and then clean the run's temporary results.
Cancellation does not stop a task already considered running by the executor.
Direct scoring needs the explicit pending-task cancellation that indexed scoring and index construction already attempt.
Waiting before temporary cleanup already exists; add no new waiting or shutdown mechanism in this focused failure change.
An abnormal worker exit can break a shared pool; do not promise that other workers survive it or continue the run.
See Python's [executor shutdown and broken-process-pool documentation](https://docs.python.org/3.13/library/concurrent.futures.html).

Preserve the current cleanup behavior and tests, including removal of completed chromosome shards after another required chromosome fails.
Keep logs and error diagnostics outside disposable scientific-result staging.
Write optional `--export-ref-metadata` sidecars into the run's temporary directory and publish them only after successful computation of all required batches.
On computation failure, remove staged reference metadata with the other temporary results.
Preserve where this optional export is currently supported; staging does not add an export option to indexed scoring.

All required chromosome calculations must succeed before their cross-chromosome aggregation proceeds.
All required batches must succeed before final publication begins.
Keep existing final publication and overwrite behavior, including the gene index's existing transaction; do not add transactions, rollback, or handling for failures during LD-score final publication.
Keep existing interruption behavior, with no new special Ctrl+C or forced-termination cleanup.
Temporary files may remain if the parent is interrupted or killed.

## LD-score diagnostic logs

Use the following layout for writing direct and indexed LD-score workflows, in serial and parallel execution:

```text
<output-dir>/diagnostics/
  ldscore.log
  chromosome-logs/
    chr1.log
    chr2.log
    ...
```

Preserve the low-level `LDScoreCalculator.run(..., output_config=None)` contract: return the prepared single-batch result without filesystem writes or a workflow log.
The chromosome-file requirement applies to writing workflows and does not require an output directory for that low-level call.

Each chromosome owns its detailed log; successive query batches within the same run append to that chromosome's file.
The parent log records overall progress and each chromosome's diagnostic path.
For a Python exception, write the full traceback to the chromosome log, and identify the chromosome, current batch when applicable, error summary, and diagnostic path in the main log.
For abnormal worker termination, report the affected in-flight tasks and their log paths; a failed future alone does not identify which chromosome caused a shared pool to break.
No traceback is guaranteed after a forced exit, and no parent log update is guaranteed when the whole job is killed.

Logs survive computation failure and follow the existing `--overwrite` ownership and preflight rules on later runs.
The per-run append requirement does not imply appending a new run into a previous run's logs.
Use the existing artifact-family ownership mechanism for the new log paths and preserve unrelated files.
Do not collect every chromosome log into the main log, add log-transport infrastructure, or change regression's current log transport.

The startup summary reports the requested workers, detected CPU count and source, task unit/count, any existing batch concurrency cap, and resolved workers.
Report the numerical backend, verified native setting or supported mode, and Arrow CPU/I/O pool sizes from the executing workers.
When a numerical limit cannot be enforced or verified, explicitly record that status instead of presenting the target as an observed setting.
These diagnostics describe configuration, not measured sustained CPU utilization.

## Final decisions and remaining boundaries

The user's latest decisions on 2026-10-12 govern:

1. Always use the single-threaded BLAS policy, irrespective of worker count or workflow.
   The earlier serial CPU-count cap, higher user-configured BLAS values, and serial tuning alternatives are superseded.
2. Preserve the low-level no-file-output API; mandatory chromosome logs apply to writing LD-score workflows.
3. Continue when numerical limits cannot be enforced or verified, and clearly record that limitation.
   This supersedes the earlier recommendation to fail initialization solely for that reason.

There is no remaining serial-versus-parallel BLAS choice.
The only qualification to the word "always" is the retained best-effort enforcement rule: every worker requests one thread, but an unsupported or unverifiable backend must be reported accurately.
Worker count still differs from software-thread count and actual CPU use because coordination, native runtime helpers, and Arrow pools are separate.
A focused local transport probe completed the gene-index migration gate.
It used two synthetic chromosomes, each with four PLINK SNPs plus 25,000 supplied annotation rows and 64 quantitative columns, and two workers.
Threads took 2.01 s with 329.3 MiB sampled peak summed RSS; spawned processes took 2.29 s with 638.2 MiB.
Each task serialized about 5.1 KB of storage descriptors and catalog metadata rather than baseline arrays.
Both runs produced the same index identity and gene-support counts, and twenty exact comparisons of baseline tables, sparse operators, counts, overlaps, and atom statistics passed.
Each setting ran once, and summed RSS counts shared pages more than once, so these numbers show the direction of the memory cost rather than a stable overhead estimate.
Construction now uses chromosome processes; the measured additional RSS is documented, and full-genome memory remains unmeasured.
Broader chromosome scaling, representative memory use, and HPC Arrow I/O remain empirical questions rather than unsettled BLAS-policy choices.
Historical design/plan documents remain historical; current-runtime documentation now describes managed serial children and the uniform numerical target.

## Implementation boundaries and validation

The implementation uses `workers`, `joblib.cpu_count()`, reusable chromosome/model processes, and one BLAS thread per worker at every worker count.
The relevant seams are [`_parallelism.py`](../../src/ldsc/_parallelism.py), [`_worker_runtime.py`](../../src/ldsc/_worker_runtime.py), [`LDScoreCalculator.run`](../../src/ldsc/ldscore_calculator.py), [`_write_chromosomes`](../../src/ldsc/_indexed_ldscore_batches.py), [`build_gene_ldscore_index`](../../src/ldsc/gene_ldscore_index.py), and [`_QueryWorkers`](../../src/ldsc/_partitioned_h2_parallel.py).
Standalone numerical entry points use [`_serial_calculation.py`](../../src/ldsc/_serial_calculation.py), reusing the current child for nested numerical calls.
The existing small crash-handler correction remains: catch `concurrent.futures.process.BrokenProcessPool`, including initial submission failures, preserving the cause.
Direct scoring additionally cancels pending futures and stages reference metadata until computation succeeds.
Both writing LD modes retain per-chromosome logs, while the prepared no-output API remains file-free.

Tests cover CPU budgets before work caps, serial children, caller-setting isolation, inherited settings, runtime verification/continuation diagnostics, scientific equivalence, operator and batch ownership, chromosome logs, failure cleanup, and overwrite behavior.
The index transport probe described above records twenty exact scientific comparisons and the measured additional process RSS.
The earlier local and Longleaf pilots still do not establish LD chromosome scaling, genome-wide memory needs, or representative HPC I/O performance.
Serial-phase and repeated-work optimizations remain later benchmarking work, without a new performance promise.
