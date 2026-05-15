# Legacy `merge-alleles` behavior

This note records the removed `munge-sumstats --merge-alleles` behavior as a
reference for future allele harmonization work.

## Inputs

`--merge-alleles` accepted a file with required columns:

- `SNP`
- `A1`
- `A2`

The reader ignored other columns, stripped whitespace from the required
columns, rejected missing required values, and accepted either tab-delimited or
whitespace-delimited files.

## Processing order

The legacy path used the merge file in two separate places:

1. During chunk parsing, it kept only raw sumstats rows whose `SNP` value was in
   the merge file.
2. After p-values were converted to signed `Z` scores, it merged the munged
   rows back onto the merge file by `SNP` with the merge file as the left table.

The second step made the output follow the merge-file SNP order and retained
merge-file SNPs that did not have matching munged data as rows with missing
values until the allele-match mask was applied.

## Allele check

The merge file's `A1` and `A2` were concatenated into an uppercase reference
pair named `MA`. For each row with nonmissing sumstats alleles, the legacy code
formed:

```text
sumstats_A1 + sumstats_A2 + reference_A1_reference_A2
```

That string was checked against `ldsc._kernel.regression.MATCH_ALLELES`.
Rows whose allele combination was not in the match table were set to missing
for all columns except `SNP`; if every compared row failed, the run raised an
error.

## Important distinction from `sumstats_snps_file`

`sumstats_snps_file` is only an identity keep-list. It filters rows by the
effective key for `GlobalConfig.snp_identifier`: `SNP` in `rsid`,
`SNP:<allele_set>` in `rsid_allele_aware`, `CHR:POS` in `chr_pos`, and
`CHR:POS:<allele_set>` in `chr_pos_allele_aware`. Headered keep-lists in
allele-aware modes may provide `A1/A2`; base modes ignore allele columns for
retention.

It does not reorder output to a reference allele file and does not perform
post-filter harmonization. Sumstats signed statistics are oriented by the
munger's own `A1` convention, and regression performs pairwise allele
orientation when matching current allele-aware sumstats and LD-score artifacts.
