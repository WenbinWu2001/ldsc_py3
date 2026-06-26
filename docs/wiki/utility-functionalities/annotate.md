## Functionality

curate binary query annotation files from bed files, given a set of baseline annotation files.

## Input

### baseline annotations

`1000G_EUR_Phase3_baseline/baseline.@.annot.gz`

Each row is a SNP, and a value of 1 for an annotation means this SNP is included in the annotation (e.g. pathway / differentially-experssed gene region).

The columns in this file are:

- `CHR`, `BP`, `SNP`, `CM`: metadata for SNPs.
- 53 baseline annotations, including a base annotation `base` of all ones and 52 other baseline annotations.

### query bed files

Four cell-type-specific (query) bed files: 

`Cerebellum_PC16.bed`,  `Cerebellum_PP1.bed`, `Cerebellum_PP3.bed`, `Hippocampus_PP1.bed`

**required format:**

TODO

## Output

per-chromosome binary annotations for all *query* annotations (i.e., the file does not contain baseline annotations), with rows aligned to the baseline annotations. You'd recycle the input baseline annotations in downstream regression by the flag `--baseline-annot-sources`.

Under output directory, you have `query.1.annot.gz`, `query.2.annot.gz`, ...,  `query.22.annot.gz`.

Each annotation file has the following columns:

- `CHR`, `BP`, `SNP`, `CM`: metadata for SNPs, copied from baseline annotations. `CM` is set as empty to indicate no-use in downstream.
- one binary annotation column for each query bed files: `Cerebellum_PC16`, `Cerebellum_PP1`  `Cerebellum_PP3`, `Hippocampus_PP1`



## Config

For four annotations, 

- **peak memory usage:** 6GB, 
- **running time:** <15 mins on longleaf 
- single-cpu, no parallel supported.



## TODO

1. why it needs such large memory?
2. benchmark time and memory usage with increasing number of annotation -- this gives practical suggestions on the users' choice of memory and running time config for slurm jobs
3. In the input and output sections, use header lines `head -n 2` as a table besides description of the columns.



​	
