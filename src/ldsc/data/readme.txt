Last updated on: 2026-10-27

See ATTRIBUTION.txt for independently reviewed data provenance and reuse
terms, including unresolved source records. The software GPLv3 declaration
does not establish the licensing of third-party datasets.

`hm3_curated_map.tsv.gz` is the packaged curated HapMap3 SNP map used by
`ldsc.load_hm3_curated_map()` and HM3 convenience flags. It is a gzip-compressed
TSV with six columns: `CHR`, `hg19_POS`, `hg38_POS`, `SNP`, `A1`, `A2`.
Positions are one-based native dbSNP build 157 placements in GRCh37 and
GRCh38; `SNP`, `A1`, and `A2` come from the LDSC `w_hm3.snplist`. Rows are
autosomal, non-strand-ambiguous, and unique in `SNP`, `(CHR, hg19_POS)`, and
`(CHR, hg38_POS)`. Rebuild it with the maintainer tool
`tools/hm3/build_hm3_curated_map.py`, whose docstring lists the pinned sources.

`hm3_chr_pos_reference.tsv.gz` is the compact test/inference reference. Rebuild
it after replacing the curated map with:

python tools/hm3/build_hm3_chr_pos_reference.py --curated-map src/ldsc/data/hm3_curated_map.tsv.gz --output src/ldsc/data/hm3_chr_pos_reference.tsv.gz

It is built as:

up to 1,000 SNPs per autosome
autosomes only
non-strand-ambiguous A1/A2
unique in both (chr, hg19_pos) and (chr, hg38_pos)
hg19_pos != hg38_pos
evenly spaced by hg38 position after filtering
