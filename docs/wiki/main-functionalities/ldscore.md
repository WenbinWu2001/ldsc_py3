

Last updated on: 2026-09-10



### Effective CM Coordinates and Metadata Export

With `--ld-wind-cm`, `ldscore` uses CM coordinates interpolated from the provided genetic map when one is supplied; otherwise, it uses the `.bim` CM values, which must be informative. When `--export-ref-metadata` is requested, the exported sidecar records these effective CM coordinates—interpolated map values when a map is provided, or the original `.bim` values otherwise—without modifying the input `.bim` file.

## Memory for many pathways

This optimization targets workloads such as 1,000 pathways in one run, each later tested separately against shared baseline categories. Direct and indexed `ldscore --query-batch-size 1000` use a positive maximum query batch size; 1000 is the default. Smaller batches reduce query workspace. Direct `--threads 1` releases chromosome working data before advancing; larger counts use bounded worker processes. Freed memory can remain reserved by the allocator, so RSS may not drop immediately. Final HM3 baseline/query LD tables remain aggregate Parquet files and may stay materialized. Temporary annotation data stay below the output directory. See the [memory design](../../current/annotation-memory-design.md) and [batch regression guide](partitioned-h2.md#testing-enrichment-for-a-large-batch-of-pathways).
