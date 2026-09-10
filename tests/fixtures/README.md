# Test Fixtures

Last updated on: 2026-09-10

The default suite keeps fixture data small and purpose-specific.

- `annotation`: tiny annotation and frequency fixtures used by annotation tests.
- `plink`: 5-individual/8-SNP PLINK fixture for low-level PLINK I/O.
- `minimal_external_resources`: 32-SNP chromosome 22 subset generated from the local 1KG 30x resources by intersecting source BIM hg38 positions with the packaged HapMap3 reference.

Do not add large external resources directly to `tests/fixtures`. For larger
parity checks, add a generator script and commit the smallest deterministic
subset that exercises the required code path.

Tests require the editable development installation (`python -m pip install -e ".[dev]"`). Run pytest and unittest from the repository root, sequentially. Individual test execution also requires that installation; no source-path bootstrap is provided.

The two `golden/*.npz` files are immutable expectations. Missing goldens fail; restore them from version control. `golden/gen_reader_golden.py` is a deliberate, reviewed regeneration tool for a known-good commit, never a test prerequisite. PLINK FAM/BIM format tests and converter-generated complete legacy suites retain format compatibility coverage; the eleven unused `formats/` fragments were retired.
