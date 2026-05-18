`hm3_curated_map.tsv.gz` is the packaged curated HapMap3 SNP map used by
`ldsc.load_hm3_curated_map()` and HM3 convenience flags. It is a gzip-compressed
TSV with full packaged columns. Public loading normalizes `CHR`, `hg19_POS`,
`hg38_POS`, `SNP`, `A1`, and `A2` while preserving extra columns.

`hm3_chr_pos_reference.tsv.gz` is the compact test/inference reference. Rebuild
it after replacing the curated map with:

python -m ldsc.hm3_reference --curated-map src/ldsc/data/hm3_curated_map.tsv.gz --output src/ldsc/data/hm3_chr_pos_reference.tsv.gz

It is built as:

500 SNPs per autosome
filtered to common SNPs with MAF >= 0.2
autosomes only
non-strand-ambiguous A1/A2
unique in both (chr, hg19_pos) and (chr, hg38_pos)
hg19_pos != hg38_pos
evenly spaced by position after filtering
