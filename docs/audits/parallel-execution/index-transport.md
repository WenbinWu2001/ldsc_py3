# Local gene-index process transport probe

Last updated on: 2026-10-12

This focused probe satisfies the agreed memory/transport check before moving chromosome construction from threads to processes.
It uses two synthetic chromosomes with four PLINK SNPs each, plus 25,000 supplied annotation rows and 64 quantitative columns per chromosome.
The annotation values occupy 12.8 MB before storage encoding.
This deliberately tests staged-input transport and process overhead, not representative LD computation or genome-wide performance.

| Executor, two workers | Wall time | Sampled peak summed RSS | Serialized arguments per chromosome |
| --- | ---: | ---: | ---: |
| Threads | 2.01 s | 329.3 MiB | 5,102 bytes |
| Spawned processes | 2.29 s | 638.2 MiB | 5,120 bytes |

RSS was sampled every 0.05 seconds across the probe parent and its compute children.
Summed RSS counts shared pages more than once; it is neither PSS nor a physical-footprint estimate.
Each setting ran once, so these timings do not establish a performance advantage or stable overhead estimate.
The first sandboxed trial could not sample process RSS and is excluded from the memory comparison.
The measured trials ran with permission to inspect only their own process IDs.

Both runs produced the same index identity and gene-support counts.
Twenty exact comparisons of baseline tables, sparse operators, counts, overlap matrices, and atom statistics passed.
Workers verified Accelerate single mode and Arrow CPU/I/O pool sizes of one during the probe.
Tasks passed storage descriptors and chromosome catalog metadata rather than baseline arrays, and returned persisted component metadata.
This supports the chosen process boundary while explicitly showing the additional memory cost.
The tiny restriction set and genotype fixture do not predict full-genome private memory or scaling; no automatic memory-safety claim follows.

The script is [index_transport_probe.py](index_transport_probe.py).
Run from the repository root with the configured development interpreter and fresh output destinations:

```bash
PYTHONPATH=. python docs/audits/parallel-execution/index_transport_probe.py \
  --mode threads --output ../pilot_runs/index-transport/threads
PYTHONPATH=. python docs/audits/parallel-execution/index_transport_probe.py \
  --mode processes --output ../pilot_runs/index-transport/processes
```

The measured artifacts are in [`pilot_runs/20261012_index_worker_transport`](../../../../pilot_runs/20261012_index_worker_transport).
Each contains `probe.json`, the generated inputs, and an index with build diagnostics.
The execution policy is documented in [parallel-execution.md](../../current/parallel-execution.md).
