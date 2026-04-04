'''
Update LD score calculation code `ldscore.py` to cater (1) R2 table as input (each row being `SNP1 SNP2 R2`), and (2) annotations are passed either in separate files or in one big annotation file (batch processed, `chr4	140087317	140301492	[binary indicator of inclusion in annot1] [binary indicator of inclusion in annot2]). Annotations are the output of `bed2annot.py` and are in the same format as the .annot files used in the original code, except that they can be batch processed (multiple annotations in one file) and can also be separated into multiple files (e.g. baseline annot vs query annot). 

Inputs:
- annotation files (allow batch input). [Need you to help decide the format. In original codebase, it is a list of file paths (`--ref-ld-chr` uses comma separated paths and  "baseline." using a glob pattern to specify the baseline annot per chr)]
- R2 table file (parquet) / reference panel : If provide an R2 table: a table with the following header: "chr, rsID_1, rsID_2, hg38_bp1, hg38_bp2, hg19_bp_1, hg19_bp_2, hg38_Uniq_ID_1, hg38_Uniq_ID_2, hg19_Uniq_ID_1, hg19_Uniq_ID_2, R2, Dprime, +/-corr", where `Uniq_ID`` is in "chr:pos:ref:alt" format.
- genome build: hg19 or hg38. Used to determine bp position to use in the R2 table. Only useful when identifier mode is `chr_pos`. If identifier mode is `rsID`, then the genome build is not relevant.

Outputs:
- .M, .M_5_50, ldscores, w_ld for each annotation, including baseline annotations (keep the format as the original code for baseline annot)

Identifier mode:
Similar to bed2annot, we will have two snp identifier modes: (1) `chr_pos` and (2) rsID. All the matching and intersection will be done based on the selected identfier mode. Read more in bed2annot.py on this.

Input reference panel:
Either an R2 table (parquet) or PLINK sets of genotyping data (same as before). It should only differ by how the chunk R2 matrix is constructed. The downstream steps of LD score calculation remain the same.

Need these helper functions (or modify or reuse existing ones):
- read from `--ref-ld-chr` the annot files and convert to proper proper format as a list (or use existing proper functions in the codebase). 
- column bind all annotations (each annotation is a column), including the baseline annotations (if annotations are in separate files). Make sure to save which columns are query annot and which are baseline annot.
- Make sure the snps are in the same order (sorted by CHR, BP).
- sanity check for each annotation file after they are normalized (e.g. snps are ordered by CHR, BP): (1) all annotation files must have the same snps in the same order in their identifier (rsID or chr_bp), i.e., the reference snp set, and throw an error if not (this is already done in the code). (2) 

Main steps:
- Process input reference panel (R2 table or genotype data in plink format) into proper format and establish indexing for window / chunk selection if needed.
- Sliding window mechanism in the original code, adapted to the new input format. 
    - first calculate the window boundaries based on the bp positions
    - fast query from R2 table to R2 matrix for a chunk: identify section in R2 table that corresponds to the window, then construct the R2 matrix for that window. 
    - normal sliding window iteration to calculate ld scores
- all annotations are processed and have their ld scores calculated all at once in this step.
- Build the .M and .M_5_50 output files. (TODO: change later. keep it like this for now.)
- Write out the output files: .M, .M_5_50, ldscores, w_ld (non-partitioned ld scores computed using the regression snps as reference snps)


Data structure:
- The chunk R2 matrix is symmetric, sparse, banded with values in [0,1]. each row and column is one SNP on a chromosome. Encode LD R2 matrix in `scipy.sparse.csr_matrix`. `dtype = float32.` 
- The annotation matrix (columns are all annotations) is tall, skinny, sparse. Each row is a SNP, each column is an annotation. For each annotation column, SNPs are usually binary coded to denote inclusion in that annotation (e.g., pathway). Could also contain continuous values. 0 denotes no contribution. Consider `scipy.sparse.csr_matrix`. `dtype = float32.` 

Later, this module will be used together with bed2annot.py in a larger wrapper function of end-to-end partitioned ldsc pipeline. In the pipeline, this script will be followed by a batch-alike partitioned ldsc wrapper function.  So you need to ensure the relevant ports are flexible.
'''


