# SNP merging policy

This draft has been reconciled into the package-wide allele-aware SNP identity
design in `docs/superpowers/specs/2026-05-15-allele-aware-snp-identity-design.md`.
The historical notes in this file are no longer a separate contract.

Current public SNP identifier modes are exactly:

- `rsid`
- `rsid_allele_aware`
- `chr_pos`
- `chr_pos_allele_aware`

The default mode is `chr_pos_allele_aware`.

## Effective Identity

Base modes are allele-blind:

- `rsid` uses only `SNP`.
- `chr_pos` uses only `CHR:POS`.

Allele-aware modes use normalized unordered allele sets after the base key:

- `rsid_allele_aware` uses `SNP:<allele_set>`.
- `chr_pos_allele_aware` uses `CHR:POS:<allele_set>`.

There is no separate boolean switch for allele matching. Alleles participate
in identity only by choosing one of the allele-aware modes.

## Drop Policy

The duplicate policy is package-wide: compute the effective merge key for the
active mode, then drop all rows in duplicate-key clusters.

Base modes must not inspect alleles for identity, duplicate filtering, or drop
reasons. Allele columns may be preserved as passive data, but retention is based
only on the base key and duplicate drops use reason `duplicate_identity`.

Allele-aware modes require usable `A1` and `A2` for sumstats, reference-panel
metadata, package-written R2 parquet endpoints, and LD-score artifacts. They
drop missing, invalid, non-SNP, identical, strand-ambiguous allele pairs,
package-wide multi-allelic base-key clusters, and duplicate effective-key
clusters.

Annotation rows are the one matching boundary where missing alleles are allowed:
if annotation metadata lacks `A1/A2`, annotation rows are aligned to the
reference panel by base key and the reference panel supplies the alleles carried
into the emitted LD-score artifact.

External raw R2 support remains base-mode only (`rsid` or `chr_pos`).
Allele-aware R2 workflows require current package-written artifacts with allele
metadata and identity provenance.
