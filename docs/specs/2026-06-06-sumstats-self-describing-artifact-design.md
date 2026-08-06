# Self-Describing Single-File Sumstats Artifact

**Date:** 2026-06-06
**Scope:** `src/ldsc/sumstats_munger.py` (munge writer + `load_sumstats`),
`src/ldsc/regression_runner.py` (genome-build guard), and a package-wide removal
of `schema_version` from artifact metadata (`src/ldsc/_kernel/snp_identity.py`,
`src/ldsc/outputs.py`, `src/ldsc/ref_panel_builder.py`,
`src/ldsc/annotation_builder.py`, `src/ldsc/_kernel/ref_panel.py`).

**Relation to prior design:** Partially supersedes
[`2026-05-16-artifact-metadata-json-scheme.md`](2026-05-16-artifact-metadata-json-scheme.md).
The sumstats artifact no longer carries a root `metadata.json` sidecar; its
identity payload moves into the Parquet footer. `schema_version` is dropped from
every artifact type's metadata.

---

## Problem

`munge-sumstats` writes a **directory** artifact, not a single file:

```
<output-dir>/
  metadata.json                      # downstream identity payload
  sumstats.parquet                   # the data
  diagnostics/
    sumstats.log
    dropped_snps/dropped.tsv.gz
```

Downstream regression takes the data file path (`<output-dir>/sumstats.parquet`)
and the loader hunts beside it for `metadata.json` to recover identity semantics
(`snp_identifier`, `genome_build`, `trait_name`). The only downstream-required
information outside the Parquet is that small identity payload; everything in
`diagnostics/` is audit-only and never consumed by regression.

Three consequences motivate this work:

1. **The data file is not self-sufficient.** Move the parquet and you lose its
   identity. The contract is "data file *plus* colocated sidecar," not one file.
2. **Legacy sumstats are rejected.** A large, well-maintained repository of
   summary statistics munged by the legacy Python-2 LDSC codebase exists. Those
   files are rsID-keyed `.sumstats.gz` (`SNP A1 A2 Z N`) with no metadata. The
   current loader hard-rejects any artifact lacking `metadata.json`, forcing
   regeneration users will not perform. Metadata must be *optional* downstream.
3. **No genome-build guard.** Regression reconciles the `snp_identifier` family
   but never checks genome build. An hg19 sumstats run against an hg38 `chr_pos`
   panel silently under-merges with no actionable guidance.

Separately, `schema_version` (always `1`) adds noise to every artifact's
metadata without earning its keep; `artifact_type` already identifies the
artifact. The user wants it removed across the package.

---

## Design Decisions

| Question | Decision |
|---|---|
| Where does the sumstats identity payload live? | In the `sumstats.parquet` **footer** as discrete `ldsc:*` key/value entries, mirroring the existing reference-panel footer convention. |
| What about `metadata.json`? | No longer written. Removed from the munge output. |
| Output layout / CLI? | Unchanged. `--output-dir` and the `diagnostics/` subdirectory stay exactly as they are. Only `metadata.json` disappears. |
| `.sumstats.gz` metadata? | None. A plain TSV carries no footer and gets no sidecar. In `both` mode only the `.parquet` is self-describing. |
| Loader resolution when metadata is absent? | Footer-or-nothing. **No** `metadata.json` sidecar search. Absent footer → `config_snapshot=None`, identifier mode inferred downstream (as legacy `.sumstats.gz` already flow). No rejection. |
| How is the sumstats identifier mode chosen for a metadata-less file? | Unchanged regression behavior: it inherits the LD-score panel's mode via the existing fallback. No new flag, no change to regression identity resolution. |
| Genome build at regression? | For `chr_pos`-family runs only: build from footer metadata if present, else auto-infer from coordinates; mismatch against the panel (or other sumstats in `rg`) → hard error asking the user to liftover. rsID-family runs skip the check. |
| Inference inconclusive? | Hard error: re-munge the sumstats so it carries build metadata, or liftover to match the panel. Never proceed on an unverified build. |
| `schema_version`? | Removed from every writer and every reader/validator package-wide. `artifact_type` becomes the sole metadata guard. |

---

## New Contract

### Munge output

```
<output-dir>/
  sumstats.parquet          # self-describing: identity in the footer
  diagnostics/              # unchanged, audit-only
    sumstats.log
    dropped_snps/dropped.tsv.gz
```

The Parquet footer carries discrete keys, mirroring
`src/ldsc/_kernel/ref_panel.py`:

```
ldsc:artifact_type   = "sumstats"
ldsc:snp_identifier  = "<rsid|rsid_allele_aware|chr_pos|chr_pos_allele_aware>"
ldsc:genome_build    = "<hg19|hg38|...>"   # empty/None for rsID-family
ldsc:trait_name      = "<trait>"
```

No `schema_version` key. Pandas' own footer entries are preserved. `.sumstats.gz`
(when `--output-format both` or `tsv.gz`) is written with no embedded metadata.

### Downstream loading

`load_sumstats(path)`:

1. If the artifact is a Parquet with `ldsc:*` footer keys → recover the identity
   payload from the footer; run the existing mode-specific dedup/allele
   validation; return a `SumstatsTable` with a populated `config_snapshot`.
2. Otherwise (legacy `.sumstats.gz`, footer-less Parquet) → return a
   `SumstatsTable` with `config_snapshot=None`, log that the identifier mode
   will be inferred from the LD-score panel at regression, and skip
   mode-specific validation. **No rejection. No sidecar lookup.**

### Regression genome-build guard

In `build_dataset`, after the effective identity mode is resolved and **only for
`chr_pos`-family** runs:

1. sumstats build := `config_snapshot.genome_build` if known, else
   `infer_chr_pos_build(sumstats.data, ...)`.
2. Compare against the LD-score panel's build (and, for `rg`, every other
   sumstats' build).
3. Mismatch → hard error naming both builds and instructing the user to liftover
   the sumstats to the panel's build and re-run. No automatic liftover.
4. Inference inconclusive (too little reference overlap) → hard error: re-munge
   or liftover.
5. rsID-family runs skip this entirely — coordinate build is irrelevant to an
   rsID merge.
6. The effective genome build is logged alongside the already-logged effective
   `snp_identifier`.

Because a metadata-less sumstats already inherits the panel's identity *family*
through the existing fallback, the `chr_pos`-vs-rsID routing falls out without
new branching.

### `schema_version` removal

`schema_version` is dropped from:

- **Writers:** `_kernel/snp_identity.identity_artifact_metadata`,
  `outputs.py` (`_result_metadata` and the partitioned-h2 / rg per-unit
  payloads), `ref_panel_builder.py`, `regression_runner.py`,
  `annotation_builder.py`, and the `ldsc:schema_version` footer key in
  `_kernel/ref_panel.py`.
- **Readers/validators:** `_kernel/snp_identity.validate_identity_artifact_metadata`
  (validate `artifact_type` only), `_kernel/ref_panel.py` (drop
  `ldsc:schema_version` from the required-key set and the parsed dict), and the
  `regression_runner.py` LD-score provenance check (require `artifact_type`
  only).

Already-built artifacts that still carry the key remain readable — readers
simply ignore the extra entry. No migration is required.

---

## Out of Scope

- No change to `--output-dir` semantics, the `diagnostics/` layout, or any CLI
  flag.
- No new regression flag; regression identity resolution is otherwise unchanged.
- No automatic liftover at regression time.
- No back-compatibility path for already-written sumstats `metadata.json`
  sidecars — they are simply no longer read.

---

## Testing Strategy (TDD)

- **Writer:** footer round-trip — munge, then `pq.read_schema(...).metadata`
  carries the expected `ldsc:*` keys and no `schema_version`; assert
  `metadata.json` is **not** written and `diagnostics/` is still present.
- **Loader:** footer-only Parquet loads with a populated `config_snapshot`; a
  legacy `.sumstats.gz` loads with `config_snapshot=None`; a footer-less Parquet
  loads with `config_snapshot=None`; no sidecar is consulted.
- **Regression build guard:** matching build (metadata or inferred) runs;
  mismatching metadata build → liftover error; metadata-less `chr_pos` sumstats
  with inferred-matching build runs; inconclusive inference → hard error;
  rsID-family run with no build is unaffected.
- **`schema_version` removal:** every artifact type's metadata omits the key and
  still validates via `artifact_type`; an artifact that still carries the key
  loads without error.
- **End-to-end:** `h2` on a footer-only Parquet produces the expected summary.
