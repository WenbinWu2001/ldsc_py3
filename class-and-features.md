# code structure planning

## main classes

Each class should come with proper parser (if needed) and input/output methods, and sanity checks. [Make suggestions on Should all methods be exposed to users? Where to store the helper functions?]

[we also need a unified format requirements for snp list inputs, e.g. require chr and bp when snp identifier is  `chr_bp` , or a column of snp ids when snp identifier is `snpid`]

### `AnnotationBundle` 

  - attibutes: 
    - baseline annotations
    - query annotations: .bed files containing gene lists, names (auto-infered from .bed), annotation (to create, row snps are defined by baselin annotations). 
    - others that may be helpful or from old codebase
  - methods: 
    - `run_bed_to_annot`
    - others that may be helpful or from old codebase

  - later this can be generalized to cell-type specific regression

### `RefPanel` (either from R2 table or PLINK files)

- attributes: 
  - the reference panel data (R2 table / PLINK) [Suggest whether I should use the path or the data objects]
  - (for R2 table) metadata / annotations for snps (such as MAF)
  - config for bookkeeping (e.g., filters used)
  - others that may be helpful or from old codebase

- methods: 
  - filtering (individuals [for PLINK], snps)
  - `convert_r2_table_to_sorted_parquet`
  - block reader for fast query (previously written as `SortedR2BlockReader`), will be  both used for ld score calculation and user queries
  - others that may be helpful

- `RefPanel` from R2 table and that from PLINK can be two separate classes. For PLINK format, you can reuse the classes from old codebase for PLINK files as well.

### `LDScore` 

Regarding reference snps: the SNPs used for LD score computation are effectively those present in both the reference panel and the annotation file, plus any optional SNP restriction list. But conceptually, the reference panel defines the candidate universe, and the annotation supplies values on that universe.

  - attributes:
    - reference snps (derived from reference panel and annotation file, also see methods). if no annotation is provided, set as all snps from reference panel.
    - regression snps (user-specified, used to compute `w_ld`, must be a subset of reference snps)
    - window size
    - others that may be helpful or from old codebase
  - methods:
    - calculate ld scores from `RefPanel` and  `AnnotationBundle` (here we need to produce a shared ref snps unverse by intersecting snps from baseline annotations and reference panel snps, and update reference snps and regression snps attributes (regression snps should always be a subset of reference snps).). outputs per-chromosome .M, .M_5_50, .ldscore, .w_ld.   [ suggest whether should I include these as outputs as attributes for later use in ldsc regressions]
    - others that may be helpful or from old codebase

### `config` (global configuration)

all classes and method should extract relevant configuration. For example, chunksize is used in ld score calculation)

- `snp_identifier` and  genome build (if use `chr_bp` as snp identifier)
- restrict snp space: optional restricted snp sets in global setting. These will restrict the regression snps,  reference panel snps and baseline annotation snps. For example, if provided hapmap3 snps, these will be the snp universe where R2 b/w pairs of snps will be summed over in ld score calculation.
- `chunksize` : used in processing sumstats and R2 table or PLINK in ld score calculation (num of rows processed at once)
-  verbose logging level
-  others that may be helpful or from old codebase

### `Sumstats` 

(skip for now)

## main features

### make annotations (requires)

### LD score estimation 

**inputs:** 

- requires `RefPanel`

- optional `AnnotationBundle` (if no annotation is provided, use all reference snps from reference panel)

- optional regression snps (if not set, default to reference snps in `LDScore`)

outputs: `LDScore` in proper format. See `LDScore` above

### regression suite

- heritability estimates (requires `Sumstats`, `LDScore`)
- partitioned heritability (requires `Sumstats`, `LDScore`)
- cross-trait genetic correlation (requires `Sumstats`, `LDScore`)

Remarks on regression suite:

- heritability estimates and partitioned heritability should share the same backbone for ld score regression. For partitioned heritability, it accepts from `AnnotationBundle` and could loop over all query annotations, each test against the baseline annotations. Cell-type specific regression is similar to partitioned heritability (don't implement yet).

- For partitioned heritability, if the annotation set (one query + baselines) does not cover the full SNP universe used for LD score regression, you need an annotation that represents “all SNPs,” typically a base/all-ones category. Throw a message for this.

### munge sumstats (sumstats)

(skip for now)
