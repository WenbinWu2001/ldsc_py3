# SHA256 Usage Decisions

Last updated on: 2026-09-10

Status: the design, removal remedies, and reader consequences were confirmed on 2026-09-10. Both fingerprint removals are implemented and verified.

## Confirmed objective

Retain only necessary SHA256 fingerprints and avoid routine or redundant hashing. Decide which fingerprints to retain or remove and what remedies their removal requires before implementation.

## Evidence baseline

The inspected repository is the local `ldsc_py3_restructured` checkout on `restructure`, at `9d3fe47`. The repository's [change discipline](../../AGENTS.md#change-discipline) permits hashes only when strictly necessary for artifact semantic identity or integrity validation and requires their exact guarantees to be documented.

The [gene-index artifact contract](gene-ldscore-index.md#artifact-and-identity-contract) describes a content-derived `index_id`, while complete builds have no cache, resume, incremental-update, or append contract. The source owners are [gene_ldscore_index.py](../../src/ldsc/gene_ldscore_index.py), `calculate_index_id`, `_builder_index_identity`, `_load_gene_ldscore_index`, and `_published_row_metadata_sha256` (formerly `_index_row_digests`).

The reference-panel binding hashes ordered physical SNP identities in [snp_identity.py](../../src/ldsc/_kernel/snp_identity.py), `sidecar_identity_sha256`, and is checked by [ldscore.py](../../src/ldsc/_kernel/ldscore.py), `_validate_index_binding`. Gene-index row verification covers row identities and metadata rather than numerical payloads. At the evidence revision, the legacy converter wrote source-file hashes in [legacy_ldscore_converter.py](../../src/ldsc/legacy_ldscore_converter.py), `_build_provenance`, without a runtime comparison against expected source checksums.

## Confirmed fingerprint dispositions

The user confirmed the following decisions on 2026-09-10 and explicitly required the existing gene-index semantic-identity design to remain unchanged.

| Fingerprint purpose | Decision | Reason and scope |
| --- | --- | --- |
| Reference-panel `ldsc:sidecar_identity_sha256` | Retain | Preserve binding of parquet row indices to ordered physical SNP identity and allele orientation. Keep the current `CHR:POS:A1:A2` payload and current validation behavior. |
| Gene-index `published_row_metadata_sha256` | Retain | Preserve verification of ordered `CHR SNP POS A1 A2` metadata for published regression rows. |
| Gene-index `index_id` | Retain unchanged | Preserve the current content-derived semantic identity and root/component identity checks. Do not replace it with a per-build identifier. |
| Gene-index baseline content digest | Retain unchanged | Preserve the existing aligned baseline/PLINK-intersection identity and annotation-value payload. |
| Gene-index PLINK BED digest | Retain unchanged | Preserve hashing of complete selected BED source files. |
| Gene-index PLINK BIM digest | Retain unchanged | Preserve hashing of the existing ordered parsed BIM metadata. |
| Gene-index selected-individual digest | Retain unchanged | Preserve hashing of selected IIDs in effective genotype order. |
| Gene-index regression-key digest | Retain unchanged | Preserve hashing of the existing sorted canonical restriction keys. |
| Gene-index explicit genetic-map digest | Retain unchanged | Preserve hashing of the loaded normalized map and the existing `bim_cm` case. |
| Legacy-conversion `source_sha256` | Remove | This is routine source provenance with no runtime integrity comparison or artifact-identity consumer. |
| Gene-index `effective_identity_sha256` | Remove | Its ordered effective-key columns are already covered by the retained published-row metadata digest. |

The six gene-index input fingerprints retain their current payload scope and canonicalization. In particular, this work does not narrow them to post-QC rows, selected genotype cells, or only map/restriction entries that contribute to the final output. Existing hashes and semantic IDs must remain identical for identical current inputs and settings.

## Confirmed design boundary

The authorized change removes two unnecessary fingerprint purposes. It preserves current scientific computations, SNP selection and ordering, allele conventions, schemas of numerical payloads, semantic identity, component matching, and publication behavior. It does not add complete numerical-payload checksums or change the current identity/integrity design.

## Confirmed removal remedies

For legacy conversion, remove `source_sha256`, `_sha256`, and that module's otherwise-unused `hashlib` import. Preserve source directories, selected prefixes/files, ignored files, conversion profile, intersection counts, count origins, coordinate evidence, and diagnostics. These fields are independent of the retired hashes and require no replacement checksum.

For gene-index rows, replace the two-digest helper with `_published_row_metadata_sha256`, returning only the retained digest. Remove the effective-row fingerprint from the writer and its requirement from the loader. The retained digest contains all effective-key columns in the same row order. Preserve duplicate-effective-identity, schema, dtype, chromosome, canonical-order, semantic-ID, and component checks. Removing the redundant metadata field does not change `index_identity` or `index_id`.

The component reader accesses named fields and does not reject extra metadata keys. Otherwise-valid existing indexes can therefore be read without rewriting them, including indexes containing the retired field. Older readers that require `effective_identity_sha256` reject newly written indexes that omit it; use an updated reader for new indexes. Preserve existing artifacts and add no schema migration or compatibility shim.

## Verification contract

Exercise conversion and index publication/loading through the existing workflow tests. New output must omit the two retired fields while retaining explicit provenance, semantic identity, row metadata, and numerical payloads. Check fixed pre-removal semantic IDs and the published-row digest in both rsID and coordinate modes; altered SNP labels, positions, and either allele must still fail validation. Otherwise-valid indexes with an unused retired field must load without any migration.

The reference-panel binding and all six gene-index input fingerprint implementations remain unchanged. Full-payload integrity verification, new checksum algorithms, input-scope narrowing, and unrelated provenance changes are outside this design.

## Verification results

- The conversion test first failed because `source_sha256` was still emitted; the gene-index tests first failed because `effective_identity_sha256` was still emitted and required. Both now pass. See [conversion tests](../../tests/test_legacy_ldscore_converter.py), `test_unpartitioned_conversion_writes_reloadable_canonical_suite_with_nullable_all_count`, and [index tests](../../tests/test_gene_ldscore_index.py), `test_index_artifact_round_trip_uses_approved_tree_and_payloads` and `test_index_loader_rejects_component_identity_metadata_and_published_row_tampering`.
- The affected conversion, gene-index, reference-panel builder, and pair-query suites passed: 248 passed, one skipped, three warnings, and two subtests passed.
- The full suite passed with the `ldsc3-dev` interpreter: 1,502 passed, one skipped, 189 warnings, and 132 subtests passed in 95.05 seconds. The skipped test checks behavior when PyArrow is absent; PyArrow was installed.
- Fixed fingerprints captured from `9d3fe47` match in both identity modes. Existing metadata with the unused retired field remains readable, while altered SNP labels, positions, and either allele still fail. Source comparison confirms all retained gene-index hash/identity functions are unchanged.
- Changed Python files compile, local documentation links resolve, and `git diff --check` passes. No source datasets or existing scientific artifacts were modified.
