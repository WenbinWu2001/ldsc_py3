'''
For pldsc, break annotation files input format into two sets: query annot and baseline annot. Each result is written separately.
- query annot allows multiple annotations, each tested separately against baseline annot
- baseline annot are read once and recycled to save comp time
Add safe handling when pldsc fails for one query annot and continues to the next query annot, instead of crashing the whole program
end to end pipeline for ldsc and pldsc:
Inputs are query and baseline annot files in column format, support batch input (multiple query annot in one file, multiple baseline annot in one file, or separate files for each annot). Also support R2 table as input reference panel. Query annot can be BED files specifying genomic regions (one BED file for one query), and will be converted to annot by bed2annot.py. 
Outputs are 
(1) one dataframe of summary of all query annotations, including the regression results (tau, se, p) and the enrichment results (enrichment, se, p) for each query annotation. This is the main output for pldsc.
(2) If user requests to save more outputs, make one directory for each query annotation using the query annotation name, and save the following files in each directory:
    - .annot file from bed2annot or provided by the user. This should have the same row snps as the baselin .annot file
    - ldsc results: .M, .M_5_50, ldscores, w_ld (though this is shared across all query annot, we can save a copy for each query annot for easier downstream use.)
    - pldsc results: partitioned pldsc results (regression results and enrichment results for that query annot and baseline annots). This is for users who want to do their own downstream analysis instead of using the summary dataframe.

By default, use batch processing throughout the pipeline. Namely, the .annot files for all query annots can be passed by a list of dataframe (column being an .annot). Avoid unnecessary I/O to save time.
SNP identifier mode is consistent throughout the pipeline, either `chr_pos` or `rsID`. 
genome_build is needed when the identifier mode is `chr_pos` to determine which bp position to use in the R2 table. If identifier mode is `rsID`, then genome_build is not relevant. If not provided (default), assume the user provide all files in the same build.



'''