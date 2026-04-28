# Test Fixtures

The default suite keeps fixture data small and purpose-specific.

- `annotation`: tiny annotation and frequency fixtures used by annotation tests.
- `formats`: tiny LDSC text/compressed file-format fixtures.
- `plink`: 5-individual/8-SNP PLINK fixture for low-level PLINK I/O.
- `minimal_external_resources`: 32-SNP chromosome 22 subset generated from the local 1KG 30x resources by intersecting source BIM hg38 positions with the packaged HapMap3 reference.

Do not add large external resources directly to `tests/fixtures`. For larger
parity checks, add a generator script and commit the smallest deterministic
subset that exercises the required code path.
