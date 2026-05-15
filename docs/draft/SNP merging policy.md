# SNP merging policy

modes change to `rsid`, `chr_pos`, `rsid_alleles` and `chr_pos_alleles`.

`snp_identifier` can be either `chr_pos` (default) or `rsid`. When alleles info are used, the internal merge key is `rsID:REF:ALT` in rsid mode or `CHR:POS:REF:ALT` in `chr_pos` mode). When alleles info are not used, the internal merge key is `rsID` in rsid mode or `CHR:POS` in `chr_pos` mode.

Choose whether to incorporate alleles for SNP matching and merging by `GlobalConfig.allele_matching`. Default to True (use allele info for matching).

Drop-duplication policy always apply on the merge key. Duplications are always checked for the output artifact or after merge, not for the original input. 

- When performing liftover for `munge-sumstat`, we only apply the drop duplication policy on the target build (since the source build is not used in downstream analysis).
- When performing liftover for `ref-panel-builder`,  we apply the drop duplication policy on each build separately (namely, the SNP sets in R2 parquet tables for hg19 and hg38 can be different).
- We check for duplication for these scenarios (feel free to suggest others):
  - In LD score calculation `ldscore`: after merging reference panel and annotation, before LD score calculation.
  - In `munge-sumstats`: After we concatenate all chunks, before any sample size filter or conversion.
  - In `h2` or `partitioned-h2` module, when merging sumstats SNPs with regression SNPs from LD score results.
  - after we merge SNPs in `ref-panel-snps-file` and original ref panel snps.
  - after merging snps with any input SNP-restriction files (`--sumstats-snps-file`, `--ref-panel-snps-file`, `--regression-snps-file`).

## Case 1: no allele info

`GlobalConfig.allele_matching=False`

Allele does not participate in the merging. Only rsid or chr_pos is used.

**Conservative approach:** whichever identifier is used (`rsID` or `chr_pos`), no duplicate is allowed in the output artifact. All snps that share the same identifier are dropped.

**example:**

If both `rs123:A:C` and`rs123:A:G` occur, both are dropped (because both use the merge key  `rs123` when no allele info is used). Same for `1:100:A:C` and `1:100:A:G` (both use the merge key `1:100`).



**Remark:** 

- merge key is created internally during merging. Not visible in output artifact. 
- For this case, we do not require allele info in `--sumstats-snps-file`, `--ref-panel-snps-file`, `--regression-snps-file`.



## Case 2: allele info

`GlobalConfig.allele_matching=True`

This is mainly to prevent multi-allelic SNPs to be mistakenly matched. 

We use `rsID:REF:ALT` (`rsID` mode) and `CHR:POS:REF:ALT` (`chr_pos` mode) as the merging key and allow `REF:ALT` to be flipped.

Alleles columns in the input files are inferred by the column alias in central registry for column name inference.



**Example:**

- `rs123:A:C` should not be merged with `rs123:A:G`

- `1:100:A:C` should not be merged with `1:100:A:G`.

- `1:100:A:C` should be merged with `1:100:C:A` (in sumstats munger, this results in a sign flip for the summary statistic).



If duplication occurs (i.e., two rows have the same merge key `rsID:REF:ALT` in rsid mode  or `CHR:POS`in `chr_pos` mode)



When allele matching is enabled, 

- we require both sumstats and reference panel (PLINK or R2 parquet tables) to have alleles info. 

- we require allele info in `--sumstats-snps-file`, `--ref-panel-snps-file`, `--regression-snps-file`. Only exception is the annotation (see below). 
- we do NOT require alleles in the input annotations as they only encode whether a genomic position belongs to a gene / pathway annotation, and it has nothing to do with multi-allelic mismatch. Recall that these are used to build unpartitioned or partitioned LD scores from reference panel. Whether given or not, the resulting LD score artifact (containing the regression SNPs) always have the alleles (see below).
  - If not given, when merging the SNPs between annotation row SNPs and references panel SNPs, only the identifier will be used (i.e., merge key is `rsID` in rsid mode  or `chr_pos`in `chr_pos` mode). Because the reference panel itself includes the alleles, these will be borrowed, and the resulting LD score files will still contain the alleles.
  - If given, the alleles will participate in the merging normally (i.e., merge key is `rsID:REF:ALT` in rsid mode or `CHR:POS:REF:ALT`in `chr_pos` mode).

- The downstream regression module uses `rsID:REF:ALT` or `CHR:POS:REF:ALT` as merge key between the LD score results file and the munged sumstats to obtain the final regression SNPs.
- Duplication policy is the same. If there are duplications on the merge key (`rsID:REF:ALT` in `rsID` mode or `CHR:POS:REF:ALT` in `chr_pos` mode)

