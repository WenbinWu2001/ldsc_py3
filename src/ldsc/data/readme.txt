`hm3_curated_map.tsv.gz` is the packaged curated HapMap3 SNP map used by
`ldsc.load_hm3_curated_map()` and HM3 convenience flags. It is a gzip-compressed
TSV with full packaged columns. Public loading normalizes `CHR`, `hg19_POS`,
`hg38_POS`, `SNP`, `A1`, and `A2` while preserving extra columns.

`hm3_chr_pos_reference.tsv.gz` is the compact test/inference reference. Rebuild
it after replacing the curated map with:

python tools/hm3/build_hm3_chr_pos_reference.py --curated-map src/ldsc/data/hm3_curated_map.tsv.gz --output src/ldsc/data/hm3_chr_pos_reference.tsv.gz

It is built as:

500 SNPs per autosome
filtered to common SNPs with MAF >= 0.2
autosomes only
non-strand-ambiguous A1/A2
unique in both (chr, hg19_pos) and (chr, hg38_pos)
hg19_pos != hg38_pos
evenly spaced by position after filtering

`protein_coding_genes.tsv.gz` is the packaged protein-coding gene catalog used
to resolve LD-score query gene lists. Its hg38 source is GENCODE Human release
49 (`gencode.v49.basic.annotation.gtf.gz`), curated by
`resources/gene_lists/protein-coding/12-pc-gene-list.R`. The catalog has one row
per canonical Ensembl gene ID and these columns:

ensgid gene_name hg38_chr hg38_start0 hg38_end hg38_strand hg19_chr hg19_start0 hg19_end hg19_strand

Coordinates are autosomal, 0-based, and half-open. `hg19_*` may be missing when
the source gene matrix has no GRCh37 coordinate; `hg38_*` is required. The
preparation record does not document the exact hg19 derivation, so the package
does not claim a liftover method for those fields. Runtime provenance records
the packaged filename, release, selected projection build, and decompressed
catalog checksum.
