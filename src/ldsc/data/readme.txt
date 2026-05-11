`hm3_curated_map.tsv.gz` is the packaged curated HapMap3 SNP map used by
`ldsc.load_hm3_curated_map()` and HM3 convenience flags. It is a gzip-compressed
TSV with full packaged columns. Public loading normalizes `CHR`, `hg19_POS`,
`hg38_POS`, and `SNP` while preserving extra columns.

`hm3_chr_pos_reference.tsv.gz` is an older compact test/inference reference,
built as:

500 SNPs per autosome
filtered to common SNPs with MAF >= 0.2
autosomes only
non-strand-ambiguous REF/ALT
conflict_flag == False
unique in both (chr, hg19_pos) and (chr, hg38_pos)
hg19_pos != hg38_pos
evenly spaced by position after filtering
