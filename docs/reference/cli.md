# Command-line reference

Last updated on: 2026-10-28

LDSC3 provides one command, `ldsc`, with one subcommand per workflow; `python -m ldsc` is equivalent.
This page lists every option of every subcommand, taken from the parser (`ldsc.cli.build_parser()`) on 2026-10-28.
`ldsc <command> --help` shows the same text for the installed version, and a future documentation build will generate this page.
Each table lists options in help order, where the help groups them into numbered steps followed by optional settings, computing resources and logging.
Task guides are in the [user guide](../user-guide/index.md).

## Conventions

- Options must be spelled in full; abbreviated options are rejected.
- Commands that write results require `--output-dir`; `plot` and `convert-h2-scale` instead write to fixed folders inside the result directory they read (`<result-dir>/plots/` and `<h2-result-dir>/postprocessing/liability-scale/`).
- `--overwrite` replaces only the files the command owns; see [logs and outputs](logs-and-outputs.md).
- `--log-level` sets the detail of the workflow log, not of the console.
- Path options accept the path tokens described in [input paths](input-paths.md).
- `--snp-identifier` and `--genome-build` follow [SNP identity and genome build](snp-identity-and-genome-build.md).

## Subcommands

| Command | Purpose |
| --- | --- |
| [`infer-genome-build`](#ldsc-infer-genome-build) | Inspect an artifact's genome build without writing files. |
| [`annotate`](#ldsc-annotate) | Find SNPs covered by BED intervals or gene boundaries. |
| [`ldscore`](#ldsc-ldscore) | Compute LD scores. |
| [`build-r2-matrix`](#ldsc-build-r2-matrix) | Build pairwise R2 matrices from PLINK genotypes. |
| [`build-gene-ldscore-index`](#ldsc-build-gene-ldscore-index) | Build a reusable gene LD-score index. |
| [`convert-legacy-ldscores`](#ldsc-convert-legacy-ldscores) | Convert selected legacy LD-score suites to LDSC3 format. |
| [`munge-sumstats`](#ldsc-munge-sumstats) | Munge GWAS summary statistics. |
| [`h2`](#ldsc-h2) | Estimate heritability from munged summary statistics and LD scores. |
| [`partitioned-h2`](#ldsc-partitioned-h2) | Estimate partitioned heritability with independent query models. |
| [`quantile-h2`](#ldsc-quantile-h2) | Summarize fitted heritability across annotation quantiles. |
| [`rg`](#ldsc-rg) | Estimate genetic correlation. |
| [`query-r2-matrix`](#ldsc-query-r2-matrix) | Look up R2 for SNP pairs in an R2 matrix. |
| [`convert-h2-scale`](#ldsc-convert-h2-scale) | Convert a saved observed-scale SNP heritability estimate to liability scale. |
| [`plot`](#ldsc-plot) | Plot a saved LDSC3 result. |

## `ldsc infer-genome-build`

Inspect an artifact's genome build without writing files.

| Option | Description |
| --- | --- |
| `PATH` (positional) | Required exact file or recognized LDSC3 artifact directory. Read build declarations or infer hg19/hg38 from SNP coordinates. No patterns or multiple paths. Print unknown when evidence is insufficient; unknown exits successfully. Creates no files or logs. |

## `ldsc annotate`

Find SNPs covered by BED intervals or gene boundaries.

| Option | Description |
| --- | --- |
| `--snp-identifier` | Match SNPs by rsID (rsid) or chromosome and position (chr_pos); allele-aware variants also use A1/A2. Default: chr_pos_allele_aware; usable alleles are required. Coordinate-family modes require --genome-build. |
| `--genome-build` | Genome build for coordinate-based SNP matching. Required with chr_pos or chr_pos_allele_aware; auto checks every grid shard and the gene coordinate catalog using declarations or bounded inference. Concrete builds must agree with input declarations; BED inherits the resolved build. No fixed default. |
| `--ref-snp-annot-sources` | Baseline annotation files defining the SNPs to annotate; export aligned baseline and query files. Exclusive with --ref-snp-bim-sources. Accepts exact paths, with space- or comma-separated tokens. Gene-list runs require every declared chromosome member. Quote patterns: '*' matches filename text; '@' substitutes chromosome numbers 1-22 by default. Use '*' to select available files for a chromosome subset where allowed. |
| `--ref-snp-bim-sources` | Headerless six-column PLINK .bim or .bim.gz files defining SNP rows; no BED/FAM companions required. Use SNP metadata only and write query annotations without synthetic baselines. Exclusive with --ref-snp-annot-sources. Accepts space- or comma-separated paths. Quote patterns: '*' matches filename text; '@' substitutes chromosome numbers 1-22 by default. Use '*' to select available files for a chromosome subset where allowed. |
| `--query-annot-bed-sources` | BED intervals to find covered grid SNPs. Accepts exact paths or '*' patterns. Supply exactly one of --query-annot-bed-sources or --query-annot-gene-list-sources; no default. '@' is not expanded; see --ref-snp-annot-sources for wildcard and quoting rules. |
| `--query-annot-gene-list-sources` | One-column gene-list files to find SNPs within gene intervals. Requires --gene-coordinate-file and explicit --padding-bp. Supply exactly one of this flag or --query-annot-bed-sources; no default. Accepts '*' patterns (see --ref-snp-annot-sources); '@' is not expanded. |
| `--padding-bp` | Extend each BED interval or gene boundary by this many base pairs at both ends, then find covered SNPs. Omission means 0 for BED input; gene-list input requires an explicit nonnegative value, including 0 for gene bodies. Set to 0 if your BED intervals are already padded to avoid double padding. |
| `--gene-coordinate-file` | Headered TSV/TSV.GZ of one-based inclusive gene coordinates. Required with --query-annot-gene-list-sources; cannot be used with BED-only input. If omitted in BED mode, no gene coordinate catalog is used. Supply an exact file path; patterns are not expanded. |
| `--gene-exclude-regions` | Exclude MHC-overlapping genes before padding, or none to exclude no genes. Default: none. Requires --query-annot-gene-list-sources when explicitly supplied. |
| `--allow-unresolved-genes` | Omit rejected gene identifiers and continue with the usable subset, recording omissions in diagnostics. Default: off (strict resolution; rejected identifiers stop the run). Requires --query-annot-gene-list-sources. Other validation errors remain fatal. |
| `--output-dir` | Required. Required destination for chromosome .annot.gz files and diagnostics. |
| `--overwrite` | Replace this command's existing output files and remove obsolete outputs from an earlier run. Default: off; stop if output files already exist. |
| `--log-level` | Set detail in the workflow log. DEBUG includes troubleshooting details. INFO includes progress and summaries. WARNING includes warnings and errors. ERROR includes errors only. Each level includes more severe messages. Default: INFO. Run headers and completion status are always recorded. |

## `ldsc ldscore`

Compute LD scores.

| Option | Description |
| --- | --- |
| `--snp-identifier` | Match SNPs by rsID (rsid) or chromosome and position (chr_pos); allele-aware variants also use A1/A2. Default: chr_pos_allele_aware; usable alleles are required. With --gene-ldscore-index-dir, inherit the stored identity and omit this flag. |
| `--genome-build` | Genome build for coordinate matching and gene intervals. Required for coordinate-family modes; direct gene-list runs default to auto. Auto checks all coordinate inputs and shards; concrete builds must agree with input declarations. Other rsid runs need no build. Gene-list runs also use this build for --regr-snps-exclude-regions. With --gene-ldscore-index-dir, inherit its build and omit this flag. |
| `--baseline-annot-sources` | Baseline annotation tables: comma-separated paths, '*' or '@' patterns. Required with direct query inputs. If omitted without queries, use a single all-SNP baseline. Cannot be combined with --gene-ldscore-index-dir. Quote patterns: '*' matches filename text; '@' substitutes chromosome numbers 1-22 by default. Use '*' to select available files for a chromosome subset where allowed. |
| `--query-annot-bed-sources` | BED files whose intervals define query SNPs; accepts comma-separated paths or '*' patterns. Requires --baseline-annot-sources; cannot be combined with --query-annot-sources, --query-annot-gene-list-sources, or --gene-ldscore-index-dir. If omitted, add no BED queries. '@' is not expanded; see --baseline-annot-sources for wildcard and quoting rules. |
| `--query-annot-gene-list-sources` | One-column gene-list files, separated by commas. Direct mode requires --baseline-annot-sources, --gene-coordinate-file, and --padding-bp; indexed mode requires --gene-ldscore-index-dir. Cannot be combined with --query-annot-sources or --query-annot-bed-sources. If omitted, add no gene-list queries. Accepts '*' patterns (see --baseline-annot-sources); '@' is not expanded. |
| `--query-annot-sources` | Query annotation tables: comma-separated paths, '*' or '@' patterns (see --baseline-annot-sources for pattern rules). Requires --baseline-annot-sources; cannot be combined with --query-annot-bed-sources, --query-annot-gene-list-sources, or --gene-ldscore-index-dir. If omitted, add no table-based queries. |
| `--padding-bp` | Extend BED intervals or gene boundaries by this many base pairs at both ends, then find covered SNPs. Omission means 0 for BED input; direct gene-list input requires an explicit value, including 0. Cannot be used with prebuilt annotations, no query, or --gene-ldscore-index-dir. Set to 0 if your BED intervals are already padded to avoid double padding. |
| `--gene-coordinate-file` | Headered TSV/TSV.GZ of gene coordinates, using one-based inclusive intervals. Required with direct --query-annot-gene-list-sources; indexed runs use their stored catalog and reject this flag. If omitted outside gene-list mode, no catalog is used. Supply an exact file path; patterns are not expanded. |
| `--control-gene-list-file` | One-column gene list added as a fixed control alongside each gene-list query. Requires --query-annot-gene-list-sources. If omitted, add no gene control. Supply an exact file path; patterns are not expanded. |
| `--gene-exclude-regions` | Exclude genes overlapping MHC before extending their intervals; none excludes no genes. Default: none for direct gene-list input. With --gene-ldscore-index-dir, use the stored policy and omit this flag. |
| `--allow-unresolved-genes` | Omit rejected gene identifiers and continue with the usable subset, recording omissions in diagnostics. Default: off (strict resolution; rejected identifiers stop the run). Requires --query-annot-gene-list-sources. Other validation errors remain fatal. |
| `--r2-matrix-dir` | R2 matrix directory written by build-r2-matrix, containing chrN_r2.parquet and required chrN_meta.tsv.gz files. The LD window must use the unit the R2 matrix was built with and be no wider. Direct mode requires exactly one of --r2-matrix-dir or --plink-prefix. Cannot be combined with --gene-ldscore-index-dir; no default. |
| `--plink-prefix` | PLINK .bed/.bim/.fam filename prefix or chromosome suite stem. Direct mode requires exactly one of --plink-prefix or --r2-matrix-dir. Cannot be combined with --gene-ldscore-index-dir; no default. Direct query runs require every declared chromosome member. Accepts '*' and '@' patterns; see --baseline-annot-sources for pattern rules. |
| `--keep-indivs-file` | File listing one individual ID (IID) per row to retain from PLINK genotypes. Requires --plink-prefix; cannot be used with --r2-matrix-dir or --gene-ldscore-index-dir. If omitted, retain all individuals. Supply an exact file path; patterns are not expanded. |
| `--genetic-map-hg19-sources` | hg19 genetic-map files used to derive centiMorgan coordinates for PLINK SNPs. Accepts comma-separated exact file paths; '*' and '@' are not expanded. If omitted, use informative PLINK CM values; cM windows require one of these sources. Direct PLINK mode only. |
| `--genetic-map-hg38-sources` | hg38 genetic-map files used to derive centiMorgan coordinates for PLINK SNPs. Accepts comma-separated exact file paths; '*' and '@' are not expanded. If omitted, use informative PLINK CM values; cM windows require one of these sources. Direct PLINK mode only. |
| `--export-ref-metadata` | Write reference-panel SNP metadata to ref_metadata/chrN_meta.tsv.gz under --output-dir. Applies to PLINK input only; rejected with --gene-ldscore-index-dir. Default: off; these additional files are not written. |
| `--ld-wind-cm` | Maximum LD-window distance on either side of each SNP, in centiMorgans; requires usable genetic-map coordinates. No default; specify exactly one of --ld-wind-snps, --ld-wind-kb, or --ld-wind-cm. Direct mode only. |
| `--ld-wind-kb` | Maximum LD-window distance on either side of each SNP, in kilobases. No default; specify exactly one of --ld-wind-snps, --ld-wind-kb, or --ld-wind-cm. Direct mode only. |
| `--ld-wind-snps` | Maximum LD-window distance on either side of each SNP, counted in retained SNPs with --plink-prefix and in R2 matrix rows with --r2-matrix-dir. No default; specify exactly one of --ld-wind-snps, --ld-wind-kb, or --ld-wind-cm. Direct mode only. |
| `--yes-really` | Allow an LD window covering an entire chromosome, bypassing the usual guard. Default: off; reject whole-chromosome windows. Direct mode only. |
| `--output-dir` | Required. Required destination for LD-score tables, SNP counts, and diagnostics. |
| `--overwrite` | Replace this command's existing output files and remove obsolete outputs from an earlier run. Default: off; stop if output files already exist. |
| `--regr-snps-file` | Headered SNP identity table selecting the regression SNPs: stored output rows and regression weights. If omitted, use bundled HapMap3 SNPs. --regr-snps-exclude-regions then applies. Direct mode only; rejected with --gene-ldscore-index-dir. Accepts an exact path or a quoted '*' pattern matching exactly one file; '@' is not expanded. |
| `--regr-snps-exclude-regions` | Exclude MHC, centromeres, both, or neither from regression SNPs after --regr-snps-file or HapMap3 selection. Default: mhc-and-centromeres. Reference SNPs and annotation counts are unchanged. Direct mode only; indexed runs inherit this setting. |
| `--ref-panel-snps-file` | Headered SNP identity table restricting the reference SNPs contributing to LD scores. If omitted, use all otherwise eligible reference SNPs. Extra columns are ignored and duplicate keys collapse. Direct mode only; rejected with --gene-ldscore-index-dir. Supply an exact file path; patterns are not expanded. |
| `--maf-min` | Keep reference SNPs with MAF at or above this value, in [0, 0.5]. If omitted, apply no additional MAF threshold. Direct mode only; indexed runs inherit their reference selection. |
| `--common-maf-min` | Minimum MAF for common-SNP annotation counts, including the boundary. Default: 0.05; does not filter reference SNPs. Direct mode only; indexed runs inherit their count threshold. |
| `--gene-ldscore-index-dir` | Use a complete prebuilt gene LD-score index instead of a live reference panel. Requires --query-annot-gene-list-sources; cannot be combined with direct-input or direct-computation flags. If omitted, compute LD scores directly. |
| `--workers` | Maximum simultaneous chromosomes in direct or indexed mode. Default: 1, one worker. Positive N requests at most N workers; -1 uses the detected available CPU count from joblib.cpu_count(); -2 uses one fewer (minimum one). Reject zero and positive requests above the detected count, then cap by available tasks. Each worker uses 1 BLAS thread, in serial and parallel execution. Detection does not measure idle CPUs or read your scheduler allocation; choose a value consistent with that allocation. More workers can use more memory. |
| `--query-batch-size` | Maximum number of query annotations processed at once. Default: 100; use a smaller positive value to reduce memory. Each completed batch has its own query file; direct mode repeats reference traversal. |
| `--snp-batch-size` | Number of reference SNPs processed per PLINK genotype batch. Default: 512; larger batches use more memory. LD scores do not depend on it. Ignored with --r2-matrix-dir; rejected with --gene-ldscore-index-dir. |
| `--log-level` | Set detail in the workflow log. DEBUG includes troubleshooting details. INFO includes progress and summaries. WARNING includes warnings and errors. ERROR includes errors only. Each level includes more severe messages. Default: INFO. Run headers and completion status are always recorded. |

## `ldsc build-r2-matrix`

Build pairwise R2 matrices from PLINK genotypes.

| Option | Description |
| --- | --- |
| `--snp-identifier` | Match SNPs by rsID (rsid) or chromosome and position (chr_pos); allele-aware variants also use A1/A2. Default: chr_pos_allele_aware; usable alleles are required. Controls duplicate removal, SNP-list matching, and the written R2 matrix identity. |
| `--source-genome-build` | Genome build of input PLINK coordinates. Default: auto, check every BIM shard and coordinate keep list before filtering. Concrete builds must agree with input declarations; cM maps follow the resolved build. |
| `--plink-prefix` | Required. Required PLINK .bed/.bim/.fam filename prefix; chromosome suite stems are also accepted. Quote patterns: '*' matches filename text; '@' substitutes chromosome numbers 1-22 by default. Use '*' to select available files for a chromosome subset where allowed. |
| `--keep-indivs-file` | File listing one individual ID (IID) per row to retain from PLINK genotypes. If omitted, retain all individuals. Accepts an exact path or a quoted '*' pattern matching exactly one file; '@' is not expanded. |
| `--ref-panel-snps-file` | Headered SNP identity table restricting reference SNPs. Coordinate identities must use the PLINK source build; duplicate keys collapse and extra columns are ignored. If omitted, use all otherwise eligible SNPs. Accepts an exact path or a quoted '*' pattern matching exactly one file; '@' is not expanded. |
| `--maf-min` | Keep reference SNPs with MAF at or above this value, in [0, 0.5]. If omitted, apply no additional MAF threshold. |
| `--ld-wind-cm` | Maximum pair distance on either side of each SNP, in centiMorgans. Requires --genetic-map-hg19-sources or --genetic-map-hg38-sources for each emitted build. No default; specify exactly one of --ld-wind-snps, --ld-wind-kb, or --ld-wind-cm. |
| `--ld-wind-kb` | Maximum pair distance on either side of each SNP, in kilobases. No default; specify exactly one of --ld-wind-snps, --ld-wind-kb, or --ld-wind-cm. |
| `--ld-wind-snps` | Maximum pair distance on either side of each SNP, counted in retained reference-panel SNP positions. No default; specify exactly one of --ld-wind-snps, --ld-wind-kb, or --ld-wind-cm. |
| `--genetic-map-hg19-sources` | hg19 genetic-map files for deriving SNP centiMorgan coordinates. Accepts exact paths, '*' and '@' patterns; see --plink-prefix for pattern rules. Required with --ld-wind-cm when hg19 output is written; otherwise optional with no map supplied by default. |
| `--genetic-map-hg38-sources` | hg38 genetic-map files for deriving SNP centiMorgan coordinates. Accepts exact paths, '*' and '@' patterns; see --plink-prefix for pattern rules. Required with --ld-wind-cm when hg38 output is written; otherwise optional with no map supplied by default. |
| `--liftover-chain-hg19-to-hg38-file` | Chain file enabling additional hg38 output from hg19 coordinates. Requires coordinate-based --snp-identifier; invalid for rsID-family modes. If omitted, do not convert hg19 input to hg38. Supply an exact file path; patterns are not expanded. |
| `--liftover-chain-hg38-to-hg19-file` | Chain file enabling additional hg19 output from hg38 coordinates. Requires coordinate-based --snp-identifier; invalid for rsID-family modes. If omitted, do not convert hg38 input to hg19. Supply an exact file path; patterns are not expanded. |
| `--output-dir` | Required. Required destination for R2 Parquet files, chromosome SNP metadata, and diagnostics. |
| `--overwrite` | Replace this command's existing R2 matrix files and remove obsolete outputs after success. A concrete chromosome prefix replaces only that chromosome; a suite replaces the full R2 matrix. Default: off; stop if output files exist. |
| `--r2-min` | Minimum unbiased R2 to store. Default: 0; zero or negative values keep every pair within the window. A positive value discards smaller R2 values; ldscore then rejects the R2 matrix, and query-r2-matrix reports discarded pairs as NaN with status=absent. |
| `--snp-batch-size` | Number of reference SNPs processed per genotype batch. Default: 512; larger batches may improve speed but use more memory. Stored pairs do not depend on it. |
| `--log-level` | Set detail in the workflow log. DEBUG includes troubleshooting details. INFO includes progress and summaries. WARNING includes warnings and errors. ERROR includes errors only. Each level includes more severe messages. Default: INFO. Run headers and completion status are always recorded. |

## `ldsc build-gene-ldscore-index`

Build a reusable gene LD-score index.

| Option | Description |
| --- | --- |
| `--snp-identifier` | Required. Required SNP matching rule: rsid uses SNP identifiers; chr_pos uses chromosome and position. The allele-aware variants also compare normalized alleles, following ldscore rules. No default. |
| `--genome-build` | Required. Required build assertion for coordinate-bearing inputs: hg19 or hg38. Input build declarations must agree after alias normalization. No inference or liftover. |
| `--baseline-annot-sources` | Required. Required baseline annotation tables covering chromosomes 1-22. Accepts exact paths; this command requires full autosomal coverage even when using a glob. Quote patterns: '*' matches filename text; '@' substitutes chromosome numbers 1-22 by default. Use '*' to select available files for a chromosome subset where allowed. |
| `--gene-coordinate-file` | Required. Required headered TSV/TSV.GZ gene coordinate catalog in the declared build, with one-based inclusive intervals. Defines the full gene set stored in the index. Supply an exact file path; patterns are not expanded. |
| `--padding-bp` | Extend each gene boundary by this many base pairs at both ends before finding covered SNPs. Default: 0, use gene bodies. Must be nonnegative and is fixed in the index. |
| `--gene-exclude-regions` | Exclude MHC-overlapping genes before padding, or none to exclude no genes. Default: mhc. The choice is fixed in the index. |
| `--plink-prefix` | Required. Required PLINK .bed/.bim/.fam filename prefix or chromosome suite covering 1-22. Accepts suite stems, '*' and '@' patterns; see --baseline-annot-sources for pattern rules. |
| `--keep-indivs-file` | File listing one individual ID (IID) per row to retain from PLINK genotypes. If omitted, retain all individuals. Supply an exact file path; patterns are not expanded. |
| `--genetic-map-hg19-sources` | hg19 genetic-map files for deriving SNP centiMorgan coordinates. Accepts comma-separated exact file paths; '*' and '@' are not expanded. If omitted, use informative PLINK CM values; one usable coordinate source is required. |
| `--genetic-map-hg38-sources` | hg38 genetic-map files for --ld-wind-cm with --genome-build hg38. Same path rules as --genetic-map-hg19-sources. |
| `--ld-wind-cm` | Maximum LD-window distance on either side of each SNP, in centiMorgans. No default; requires informative PLINK CM values or a matching-build genetic map. Specify exactly one LD window. |
| `--ld-wind-kb` | Maximum LD-window distance in kilobases. No default; specify exactly one LD window. |
| `--ld-wind-snps` | Maximum LD-window distance in retained SNPs after identity and genotype filtering. No default; specify exactly one LD window. Whole-chromosome windows are rejected. |
| `--output-dir` | Required. Required destination for the complete reusable gene LD-score index and diagnostics. |
| `--overwrite` | Replace an existing complete index after the new build succeeds. Default: off; stop if the destination already exists. |
| `--regr-snps-file` | Headered SNP identity table selecting the regression SNPs: stored output rows and regression weights. If omitted, use bundled HapMap3 SNPs. --regr-snps-exclude-regions then applies. Supply an exact file path; patterns are not expanded. |
| `--regr-snps-exclude-regions` | Exclude MHC, centromeres, both, or neither from --regr-snps-file or HapMap3 output SNPs. Default: mhc-and-centromeres. Reference SNPs and annotation counts are unchanged. |
| `--maf-min` | Keep reference SNPs with MAF at or above this value, in [0, 0.5]. If omitted, apply no additional MAF threshold. |
| `--common-maf-min` | Minimum MAF for common-SNP annotation counts, including the boundary. Default: 0.05; does not filter reference SNPs. |
| `--workers` | Maximum simultaneous chromosome builds. Default: 1, one worker. Positive N requests at most N workers; -1 uses the detected available CPU count from joblib.cpu_count(); -2 uses one fewer (minimum one). Reject zero and positive requests above the detected count, then cap by available tasks. Each worker uses 1 BLAS thread, in serial and parallel execution. Detection does not measure idle CPUs or read your scheduler allocation; choose a value consistent with that allocation. More workers can use more memory. |
| `--snp-batch-size` | Number of reference SNPs processed per PLINK genotype batch. Default: 512; larger batches use more memory. LD scores do not depend on it. |
| `--atom-batch-size` | Number of disjoint gene-overlap regions processed together during index construction. Default: 64; a smaller positive value reduces memory use. |
| `--log-level` | Set detail in the workflow log. DEBUG includes troubleshooting details. INFO includes progress and summaries. WARNING includes warnings and errors. ERROR includes errors only. Each level includes more severe messages. Default: INFO. Run headers and completion status are always recorded. |

## `ldsc convert-legacy-ldscores`

Convert selected legacy LD-score suites to LDSC3 format.

| Option | Description |
| --- | --- |
| `--legacy-reference-dir` | Required. Required directory containing a complete chromosome 1-22 legacy reference LD-score suite and .l2.M_5_50 common-SNP count files. Baseline partitioned suites also require full annotations and --legacy-frequency-dir. |
| `--legacy-weight-dir` | Required. Required directory containing chromosome 1-22 single-column regression-weight LD scores. May equal --legacy-reference-dir for an unpartitioned suite. |
| `--legacy-frequency-dir` | Directory containing chromosome frequency files used to reconstruct common-SNP counts and overlap. Required for baseline partitioned --legacy-reference-dir input; optional for unpartitioned input, with no frequency suite read when omitted. |
| `--snp-identifier` | Output SNP matching rule: rsid uses SNP identifiers; chr_pos uses chromosome and position. Default: rsid. Both ignore alleles; chr_pos requires a resolved genome build. |
| `--genome-build` | Genome build of legacy reference coordinates; no liftover is performed. Default: auto, check every retained reference chromosome for chr_pos identity. A concrete build is required to agree with input declarations; rsid identity needs no coordinate-build inference. |
| `--output-dir` | Required. Required destination for the converted LD-score tables, counts, and diagnostics. Source files are not changed. |
| `--overwrite` | Replace this command's existing output files and remove obsolete outputs from an earlier run. Default: off; stop if output files already exist. |
| `--log-level` | Set detail in the workflow log. DEBUG includes troubleshooting details. INFO includes progress and summaries. WARNING includes warnings and errors. ERROR includes errors only. Each level includes more severe messages. Default: INFO. Run headers and completion status are always recorded. |

## `ldsc munge-sumstats`

Munge GWAS summary statistics.

| Option | Description |
| --- | --- |
| `--snp-identifier` | Match SNPs by rsID (rsid) or chromosome and position (chr_pos); allele-aware variants also use A1/A2. Default: chr_pos_allele_aware; usable alleles are required. Coordinate-family modes also require --output-genome-build. |
| `--source-genome-build` | Genome build of source coordinates. Default: auto, check raw chromosomes and coordinate keep lists; stop with a grouped report if evidence conflicts or is insufficient. Concrete builds are checked against declarations and cannot be used with rsid or rsid_allele_aware identity. |
| `--output-genome-build` | Required output build for chr_pos and chr_pos_allele_aware identity. If different from the source, packaged HM3 uses automatic quick liftover; custom or unrestricted SNPs require --liftover-chain-file, which overrides quick liftover when supplied. Matching builds need no liftover. Cannot be used with rsid or rsid_allele_aware identity; no default. |
| `--liftover-chain-file` | Chain file in the source-to-output direction for liftover to --output-genome-build. Requires chr_pos or chr_pos_allele_aware --snp-identifier and an explicit --output-genome-build. If omitted, the default packaged HM3 restriction automatically uses quick liftover with package-bundled reference HM3 metadata when builds differ; supply this flag to disable quick liftover and use the chain instead. Required when builds differ with --sumstats-snps-file or --no-snp-restriction; ignored when builds match. Accepts an exact path or a quoted '*' pattern matching exactly one file; '@' is not expanded. |
| `--raw-sumstats-file` | Required. Required nonempty raw GWAS summary-statistics file path, including with --infer-only. Columns are detected from their headers unless overridden. Accepts an exact path or a quoted '*' pattern matching exactly one file; '@' is not expanded. |
| `--trait-name` | Trait label stored in Parquet metadata and used for data filenames after replacing unsafe characters with underscores. Default: no explicit label; use sumstats.parquet or sumstats.gz. |
| `--input-format` | Input layout: plain reads ordinary tables, daner-old reads counts from frequency headers, and daner-new reads case/control count columns. Default: auto, detect the layout from headers. |
| `--N-col` | Input per-SNP sample-size column; matched case-insensitively. Overrides inferred case/control columns; cannot be combined with --N-cas-col or --N-con-col. If omitted, detect sample-size columns from headers. |
| `--N-cas-col` | Input per-SNP case-count column; matched case-insensitively. Requires --N-con-col and cannot be combined with --N-col. The pair overrides inferred N; if omitted, detect sample-size columns from headers. |
| `--N-con-col` | Input per-SNP control-count column; matched case-insensitively. Requires --N-cas-col and cannot be combined with --N-col. The pair overrides inferred N; if omitted, detect sample-size columns from headers. |
| `--N` | Constant sample-size fallback when per-SNP N or paired case/control columns are absent. Input columns take precedence; --N takes precedence over --N-cas/--N-con. If omitted, use input columns or both case/control constants. |
| `--N-cas` | Constant case-count fallback. Requires --N-con when used; the constants are added only when input sample-size columns and --N are absent. Default: no case-count constant. |
| `--N-con` | Constant control-count fallback. Requires --N-cas when used; the constants are added only when input sample-size columns and --N are absent. Default: no control-count constant. |
| `--output-dir` | Required. Required destination for munged summary statistics and diagnostics. Still required with --infer-only, which creates no files or directories. |
| `--overwrite` | Replace this command's existing output files and remove obsolete outputs from an earlier run. Default: off; stop if output files already exist. Has no effect with --infer-only. |
| `--output-format` | Write <trait>.parquet, <trait>.sumstats.gz, or both; without --trait-name, use sumstats.parquet or sumstats.gz. Default: parquet; output choice does not change SNP filtering or frequency preservation. |
| `--infer-only` | Inspect headers and check selected columns for all-missing values in bounded chunks, then print a suggested command without writing files. Requires --raw-sumstats-file and --output-dir. Default: off; perform the full conversion. |
| `--sumstats-snps-file` | Replace the default packaged HapMap3 restriction with this headered identity table, matched in the source genome build. Mutually exclusive with --no-snp-restriction. Requires --liftover-chain-file when source and output builds differ, even for a custom SNP keep-list. Same exact-one '*' pattern rules as --raw-sumstats-file; '@' is not expanded. |
| `--no-snp-restriction` | Disable the default packaged HapMap3 keep-list restriction. Ordinary QC still applies. Mutually exclusive with --sumstats-snps-file; requires --liftover-chain-file when source and output builds differ. Default: off; restrict to HapMap3. |
| `--N-min` | Minimum per-SNP sample size to retain. If omitted or zero, use the 90th percentile of N divided by 1.5. Constant sample sizes are assigned after this filter and are not filtered. |
| `--nstudy-min` | Minimum number of contributing studies to retain. If omitted or zero, keep the maximum observed count. Applies only when NSTUDY is present and per-SNP sample sizes are absent. |
| `--info-min` | Keep SNPs with INFO at or above this value. Default: 0.9 when an INFO column is selected; without one the default filter is skipped, and an explicit value is an error. Use --info or --info-list to select a nonstandard column. |
| `--maf-min` | Keep SNPs whose minor-allele frequency min(FRQ, 1-FRQ) is at or above this value. Default: 0.01 when a frequency column is selected; without one the default filter is skipped, and an explicit value is an error. Output FRQ preserves the original frequency. |
| `--snp` | Input column containing SNP identifiers; matched case-insensitively. If omitted, detect from headers. SNP identifiers are required for rsid and rsid_allele_aware identity. |
| `--chr` | Input chromosome column; matched case-insensitively. If omitted, detect from headers. Required in the input for coordinate-family modes; pair with a position column. |
| `--pos` | Input base-pair position column; matched case-insensitively. If omitted, detect from headers. Required in the input for coordinate-family modes; pair with a chromosome column. |
| `--a1` | Input effect-allele column, relative to which the signed statistic is defined; matched case-insensitively. If omitted, detect from headers. Allele-aware identity requires both A1 and A2 in the input. |
| `--a2` | Input other-allele column; matched case-insensitively. If omitted, detect from headers. Allele-aware identity requires both A1 and A2 in the input. |
| `--p` | Input p-value column used to calculate Z; matched case-insensitively. If omitted, detect from headers. Input p-values must be in (0, 1]. |
| `--signed-sumstats` | Signed-statistic column and null value, separated by a comma; for example Z,0 or OR,1. Column names are case-insensitive and effects are relative to A1. If omitted, detect from headers; --a1-inc bypasses signed-statistic use. |
| `--a1-inc` | Assert that A1 increases the trait for every SNP and derive positive Z from p-values, bypassing signed statistics. Default: off; use the detected or --signed-sumstats column to determine direction. |
| `--frq` | Input allele-frequency or MAF column; matched case-insensitively. If omitted, detect from headers; if absent, omit FRQ and its filter. Selected values are preserved as FRQ without conversion to MAF. |
| `--info` | Input scalar INFO column; matched case-insensitively. If omitted, detect from headers; if absent, skip INFO filtering. Select one INFO source, using --info-list for comma-separated values within a cell. |
| `--info-list` | Comma-separated names of input columns whose cells contain comma-separated INFO values, such as 0.852,0.113,NA. Average nonmissing values for --info-min filtering; names are case-insensitive. If omitted, detect supported INFO-list columns; competing INFO sources are rejected. |
| `--nstudy` | Input column counting studies contributing to each SNP; matched case-insensitively. If omitted, detect from headers. --nstudy-min filtering applies only when per-SNP sample sizes are absent. |
| `--ignore` | Comma-separated input column names to exclude from detection and reading; matched case-insensitively. Default: ignore no columns. Do not also select an ignored column with an explicit column flag. |
| `--chunksize` | Read this many input rows per chunk. Larger chunks use more memory. Default: 1,000,000. |
| `--log-level` | Set detail in the workflow log. DEBUG includes troubleshooting details. INFO includes progress and summaries. WARNING includes warnings and errors. ERROR includes errors only. Each level includes more severe messages. Default: INFO. Run headers and completion status are always recorded. The selected liftover method and mapping/drop summary are always recorded in diagnostics/sumstats.log and printed to stdout after a successful run. |

## `ldsc h2`

Estimate heritability from munged summary statistics and LD scores.

| Option | Description |
| --- | --- |
| `--sumstats-file` | Required. Required munged summary-statistics Parquet or legacy .sumstats(.gz) file. SNP identity must be compatible with --ldscore-dir. Accepts an exact path or a quoted '*' pattern matching exactly one file; '@' is not expanded. |
| `--ldscore-dir` | Required. Required LD-score directory from ldscore or convert-legacy-ldscores, containing predictors, regression weights, and SNP counts. |
| `--trait-name` | Trait label used in summaries. If omitted, use the stored trait name or derive one from the input filename. |
| `--output-dir` | Required. Required destination for regression estimates and diagnostics. |
| `--overwrite` | Replace this command's existing output files and remove obsolete outputs from an earlier run. Default: off; stop if output files already exist. |
| `--samp-prev` | Sample case fraction, strictly between 0 and 1, for liability-scale estimates. Requires --pop-prev; omit both, or pass nan for both, to report observed-scale estimates only. |
| `--pop-prev` | Population case fraction, strictly between 0 and 1, for liability-scale estimates. Requires --samp-prev; omit both, or pass nan for both, to report observed-scale estimates only. |
| `--count-kind` | Reference SNP counts used to scale regression estimates: common uses common-SNP counts, all uses all-SNP counts. Default: common. The requested counts must be in the LD-score directory; the other kind is never substituted. |
| `--intercept-h2` | Fix the heritability intercept to this value; use 1 for the standard fixed intercept. Cannot be combined with an explicit --two-step-cutoff. If omitted, estimate the intercept. |
| `--two-step-cutoff` | First-step chi-square cutoff, including the boundary. Requires a single-annotation model with a free intercept; cannot be combined with --intercept-h2. If omitted, use 30 for such models; otherwise disable two-step estimation. |
| `--chisq-max` | Keep SNPs with chi-square at or below this value. If omitted, use max(0.001 * maximum N, 80) for multi-annotation models; single-annotation models have no chi-square cap. |
| `--num-blocks` | Number of SNP blocks for jackknife standard errors; must be at least 2. Default: 200, capped at the number of retained SNPs. |
| `--allow-identity-downgrade` | Allow allele-aware and allele-unaware inputs of the same mode family by matching without alleles. rsID-family and coordinate-family modes still cannot be mixed. Default: off; require compatible SNP identifier modes. |
| `--log-level` | Set detail in the workflow log. DEBUG includes troubleshooting details. INFO includes progress and summaries. WARNING includes warnings and errors. ERROR includes errors only. Each level includes more severe messages. Default: INFO. Run headers and completion status are always recorded. |

## `ldsc partitioned-h2`

Estimate partitioned heritability with independent query models.

| Option | Description |
| --- | --- |
| `--sumstats-file` | Required. Required munged summary-statistics Parquet or legacy .sumstats(.gz) file. SNP identity must be compatible with --ldscore-dir. Accepts an exact path or a quoted '*' pattern matching exactly one file; '@' is not expanded. |
| `--ldscore-dir` | Required. Required LD-score directory from ldscore or convert-legacy-ldscores, containing predictors, regression weights, and SNP counts. |
| `--trait-name` | Trait label used in summaries. If omitted, use the stored trait name or derive one from the input filename. |
| `--output-dir` | Required. Required destination for regression estimates and diagnostics. |
| `--overwrite` | Replace this command's existing output files and remove obsolete outputs from an earlier run. Default: off; stop if output files already exist. |
| `--continue-on-query-error` | Skip and record queries whose model preparation, regression, or result calculation fails; publish successful fits and write diagnostics/query_status.tsv plus error tracebacks in the log. Default: off; any query error prevents publication. Shared input/output failures and runs with no successful queries still fail. |
| `--summary-sort-by` | Order partitioned_h2.tsv by the selected column. Default: auto, coefficient-p for per-query fits and category for baseline-only fits. P-values sort ascending; other numeric values sort descending. |
| `--samp-prev` | Sample case fraction, strictly between 0 and 1, for liability-scale estimates. Requires --pop-prev; omit both, or pass nan for both, to report observed-scale estimates only. |
| `--pop-prev` | Population case fraction, strictly between 0 and 1, for liability-scale estimates. Requires --samp-prev; omit both, or pass nan for both, to report observed-scale estimates only. |
| `--count-kind` | Reference SNP counts used to scale regression estimates: common uses common-SNP counts, all uses all-SNP counts. Default: common. The requested counts must be in the LD-score directory; the other kind is never substituted. |
| `--intercept-h2` | Fix the heritability intercept to this value; use 1 for the standard fixed intercept. Cannot be combined with an explicit --two-step-cutoff. If omitted, estimate the intercept. |
| `--two-step-cutoff` | First-step chi-square cutoff, including the boundary. Requires a single-annotation model with a free intercept; cannot be combined with --intercept-h2. If omitted, use 30 for such models; otherwise disable two-step estimation. |
| `--chisq-max` | Keep SNPs with chi-square at or below this value. If omitted, use max(0.001 * maximum N, 80) for multi-annotation models; single-annotation models have no chi-square cap. |
| `--num-blocks` | Number of SNP blocks for jackknife standard errors; must be at least 2. Default: 200, capped at the number of retained SNPs. |
| `--allow-identity-downgrade` | Allow allele-aware and allele-unaware inputs of the same mode family by matching without alleles. rsID-family and coordinate-family modes still cannot be mixed. Default: off; require compatible SNP identifier modes. |
| `--workers` | Maximum complete query-plus-baseline models running simultaneously; also capped by query batch size. Default: 1, one worker. Positive N requests at most N workers; -1 uses the detected available CPU count from joblib.cpu_count(); -2 uses one fewer (minimum one). Reject zero and positive requests above the detected count, then cap by available tasks. Each worker uses 1 BLAS thread, in serial and parallel execution. Detection does not measure idle CPUs or read your scheduler allocation; choose a value consistent with that allocation. More workers can use more memory. |
| `--query-batch-size` | Maximum number of query annotations loaded at once; each query is fitted separately with the shared baseline. Default: 100; a smaller positive value reduces memory and caps active query workers. |
| `--log-level` | Set detail in the workflow log. DEBUG includes troubleshooting details. INFO includes progress and summaries. WARNING includes warnings and errors. ERROR includes errors only. Each level includes more severe messages. Default: INFO. Run headers and completion status are always recorded. |

## `ldsc quantile-h2`

Summarize fitted heritability across annotation quantiles.

| Option | Description |
| --- | --- |
| `--partitioned-h2-result-dir` | Required. Required saved joint model: a baseline-only partitioned-h2 result or one diagnostics/query_annotations/<query>/ directory. An aggregate multi-query result root is not a single model. |
| `--target-annot-sources` | Required. Required annotation files containing --target-annot-column, used to divide SNPs into quantiles. Accepts exact paths. The target may be an external score or a fitted annotation. Quote patterns: '*' matches filename text; '@' substitutes chromosome numbers 1-22 by default. Use '*' to select available files for a chromosome subset where allowed. |
| `--target-annot-column` | Required. Required unique column name in --target-annot-sources whose values define the quantiles. |
| `--num-quantiles` | Number of low-to-high target quantiles; must be at least 2. Default: 5. Ties stay in the lower quantile; an empty quantile causes an error. |
| `--target-missing-value` | One target-value token to exclude before forming quantiles, for example NaN or -999. Default: exclude no value; unselected nonnumeric or nonfinite values cause an error. Zero remains a valid score unless explicitly selected. |
| `--baseline-annot-sources` | Required. Required original baseline annotation files used to fit the selected model. Accepts exact paths, '*' and '@' patterns; see --target-annot-sources for pattern rules. Resupply their original values and SNP coverage. |
| `--ref-metadata-sources` | Required. Required reference SNP metadata files with chromosome, position, SNP ID, and MAF; allele-aware matching also uses A1/A2. Accepts exact paths, '*' and '@' patterns; see --target-annot-sources for pattern rules. Use original R2 matrix sidecars or PLINK metadata exported by ldscore --export-ref-metadata. |
| `--query-annot-bed-sources` | Original BED query files for the selected per-query fit; accepts paths or '*' patterns. Cannot be combined with --query-annot-sources or --query-annot-gene-list-sources. Omit for other query types; reproduce the original --padding-bp. '@' is not expanded; see --target-annot-sources for wildcard and quoting rules. |
| `--query-annot-gene-list-sources` | Original one-column gene-list files for the selected per-query fit. Requires --gene-coordinate-file; cannot be combined with --query-annot-sources or --query-annot-bed-sources. Omit for other query types; reproduce the original gene settings. Accepts '*' patterns (see --target-annot-sources); '@' is not expanded. |
| `--query-annot-sources` | Original query annotation tables for the selected per-query fit; accepts paths, '*' and '@' patterns (see --target-annot-sources). Cannot be combined with --query-annot-bed-sources or --query-annot-gene-list-sources. Omit for baseline-only fits; per-query fits require their matching source. |
| `--padding-bp` | Extend BED intervals or gene boundaries by this many base pairs at both ends before finding covered SNPs. Default: 0. Use the original LD-score setting; not applicable to prebuilt annotation tables. Use 0 for BED files that already include the padding used in the fitted model; avoid double padding. |
| `--gene-coordinate-file` | Original one-based inclusive gene-coordinate TSV/TSV.GZ catalog. Requires --query-annot-gene-list-sources and is required by that flag; no default catalog. Supply an exact file path; patterns are not expanded. |
| `--control-gene-list-file` | Original fixed-control gene list, when the selected fit used one. Requires --query-annot-gene-list-sources. If omitted, add no gene control. Supply an exact file path; patterns are not expanded. |
| `--gene-exclude-regions` | Exclude MHC-overlapping genes before padding, or none to exclude no genes. Default: none. Selecting mhc requires --query-annot-gene-list-sources; reproduce the original LD-score setting. |
| `--allow-unresolved-genes` | Omit rejected gene identifiers and continue with the usable subset, recording omissions in diagnostics. Default: off (strict resolution; rejected identifiers stop the run). Requires --query-annot-gene-list-sources. Other validation errors remain fatal. Reproduce the original LD-score policy. |
| `--output-dir` | Required. Required destination for quantile heritability, standardized coefficients, and diagnostics. |
| `--overwrite` | Replace this command's existing output files and remove obsolete outputs from an earlier run. Default: off; stop if output files already exist. |
| `--log-level` | Set detail in the workflow log. DEBUG includes troubleshooting details. INFO includes progress and summaries. WARNING includes warnings and errors. ERROR includes errors only. Each level includes more severe messages. Default: INFO. Run headers and completion status are always recorded. |

## `ldsc rg`

Estimate genetic correlation.

| Option | Description |
| --- | --- |
| `--sumstats-sources` | Required. Required munged summary-statistics Parquet or legacy .sumstats(.gz) files; accepts exact paths or quoted '*' patterns matching filename text and requires at least two files. SNP identity must be compatible with --ldscore-dir. '@' is not expanded. |
| `--ldscore-dir` | Required. Required LD-score directory from ldscore or convert-legacy-ldscores, containing predictors, regression weights, and SNP counts. |
| `--anchor-trait` | Trait name or input file path to compare against every other trait. Must identify one of --sumstats-sources. If omitted, estimate every unique trait pair. |
| `--output-dir` | Required. Required destination for regression estimates and diagnostics. |
| `--overwrite` | Replace this command's existing output files and remove obsolete outputs from an earlier run. Default: off; stop if output files already exist. |
| `--write-per-pair-detail` | Also write detailed results for each tested pair under diagnostics/pairs/. Default: off; write aggregate pair and trait tables only. |
| `--samp-prev` | Comma-separated sample case fractions in resolved --sumstats-sources order, one per trait; each in (0, 1) or nan for a quantitative trait. Requires --pop-prev; cannot be combined with --prevalence-manifest. If neither form is supplied, report observed-scale estimates only. |
| `--pop-prev` | Comma-separated population case fractions in resolved --sumstats-sources order, one per trait; each in (0, 1) or nan. Requires --samp-prev; cannot be combined with --prevalence-manifest. If neither form is supplied, report observed-scale estimates only. |
| `--prevalence-manifest` | TSV with trait_name, samp_prev, and pop_prev columns, matched to exact input trait names. Cannot be combined with --samp-prev or --pop-prev. If omitted, use paired prevalence lists when supplied; otherwise report observed-scale estimates only. |
| `--count-kind` | Reference SNP counts used to scale regression estimates: common uses common-SNP counts, all uses all-SNP counts. Default: common. The requested counts must be in the LD-score directory; the other kind is never substituted. |
| `--intercept-h2` | Fix the heritability intercept to the same value for both traits in every pair; use 1 for the standard fixed intercept. Cannot be combined with an explicit --two-step-cutoff. If omitted, estimate the intercepts. |
| `--intercept-gencov` | Fix the genetic-covariance intercept to the same value for every trait pair; use 0 for the standard fixed intercept. Cannot be combined with --two-step-cutoff. Single-annotation fits also require --intercept-h2 because automatic two-step estimation otherwise conflicts with this fixed value. If omitted, estimate the covariance intercepts. |
| `--two-step-cutoff` | First-step chi-square cutoff, including the boundary. Requires a single-annotation model with free intercepts; cannot be combined with --intercept-h2 or --intercept-gencov. If omitted, use 30 with a free heritability intercept in single-annotation fits; otherwise disable it. |
| `--chisq-max` | Keep SNPs satisfying Z1^2 * Z2^2 <= VALUE^2 for each trait pair. If omitted, apply no chi-square product cap. This differs from the single-trait chi-square filter. |
| `--num-blocks` | Number of SNP blocks for jackknife standard errors; must be at least 2. Default: 200, capped at the number of retained SNPs. |
| `--allow-identity-downgrade` | Allow allele-aware and allele-unaware inputs of the same mode family by matching without alleles. rsID-family and coordinate-family modes still cannot be mixed. Default: off; require compatible SNP identifier modes. |
| `--log-level` | Set detail in the workflow log. DEBUG includes troubleshooting details. INFO includes progress and summaries. WARNING includes warnings and errors. ERROR includes errors only. Each level includes more severe messages. Default: INFO. Run headers and completion status are always recorded. |

## `ldsc query-r2-matrix`

Look up R2 for SNP pairs in an R2 matrix.

| Option | Description |
| --- | --- |
| `--r2-matrix-dir` | Required. Required R2 matrix directory produced by build-r2-matrix, or its root containing genome-build subdirectories. Use --genome-build to select a subdirectory when needed. |
| `--genome-build` | Select the hg19 or hg38 subdirectory under --r2-matrix-dir. If omitted, use the supplied build-specific directory or its sole available build; ambiguous build roots require a selection. |
| `--pairs-file` | Required. Required TSV/CSV of SNP pairs; use - to read standard input. Endpoint columns use _1/_2 suffixes and must match --snp-identifier, including allele columns for allele-aware matching. |
| `--snp-identifier` | SNP matching rule for pair endpoints: the mode recorded by the R2 matrix (the default) or the base form of that mode (rsid or chr_pos), which ignores alleles. Pair columns must support the selected mode. |
| `--output-dir` | Required. Required destination for query_r2_matrix.tsv (r2, sign_r, r, status plus input columns) and diagnostics. sign_r and r use query allele orientation. |
| `--overwrite` | Replace this command's existing output files and remove obsolete outputs from an earlier run. Default: off; stop if output files already exist. |
| `--log-level` | Set detail in the workflow log. DEBUG includes troubleshooting details. INFO includes progress and summaries. WARNING includes warnings and errors. ERROR includes errors only. Each level includes more severe messages. Default: INFO. Run headers and completion status are always recorded. |

## `ldsc convert-h2-scale`

Convert a saved observed-scale SNP heritability estimate to liability scale.

| Option | Description |
| --- | --- |
| `--h2-result-dir` | Required. Required saved unpartitioned h2 result directory. Conversion outputs are written below postprocessing/liability-scale/ in this directory. |
| `--samp-prev` | Required. Required sample case fraction, strictly between 0 and 1. Also supply exactly one of --pop-prev or --pop-prev-range. |
| `--pop-prev` | Population case fraction, strictly between 0 and 1, for one conversion. Requires --samp-prev; choose exactly one of --pop-prev or --pop-prev-range. No default. |
| `--pop-prev-range` | Minimum and maximum population case fractions for a sensitivity curve; both strictly between 0 and 1, with MIN < MAX. Requires --samp-prev; choose exactly one of this flag or --pop-prev. No default. |
| `--num-points` | Number of equally spaced prevalence values, including both range endpoints; must be at least 2. Used with --pop-prev-range; ignored for a single --pop-prev. Default: 201. |
| `--overwrite` | Replace existing liability-scale conversion files below the input result directory. Default: off; stop if these output files exist. |
| `--log-level` | Set detail in the workflow log. DEBUG includes troubleshooting details. INFO includes progress and summaries. WARNING includes warnings and errors. ERROR includes errors only. Each level includes more severe messages. Default: INFO. Run headers and completion status are always recorded. |

## `ldsc plot`

Plot a saved LDSC3 result.

| Option | Description |
| --- | --- |
| `--result-dir` | Required. Required saved LDSC3 result directory to plot. The result type selects the plot automatically; files are written below plots/ in this directory. |
| `--overwrite` | Replace existing plot files for the selected result type. Default: off; stop if these output files exist. |
| `--log-level` | Set detail in the workflow log. DEBUG includes troubleshooting details. INFO includes progress and summaries. WARNING includes warnings and errors. ERROR includes errors only. Each level includes more severe messages. Default: INFO. Run headers and completion status are always recorded. |
