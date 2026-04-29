`hm3_chr_pos_reference.tsv.gz` ships an 11,000-SNP reference panel in `hm3_chr_pos_reference.tsv.gz`, built as:

500 SNPs per autosome
filtered to common SNPs with MAF >= 0.2
autosomes only
non-strand-ambiguous REF/ALT
conflict_flag == False
unique in both (chr, hg19_pos) and (chr, hg38_pos)
hg19_pos != hg38_pos
evenly spaced by position after filtering