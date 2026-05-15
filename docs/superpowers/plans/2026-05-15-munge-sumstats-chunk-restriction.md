# Munge-Sumstats Chunk Restriction Implementation Plan

## Summary

Refactor `_kernel.sumstats_munger` so HM3 / `--sumstats-snps-file` filtering is
prepared once before parsing and applied inside `parse_dat()` before retained
chunks are appended. Preserve current aggregate logs and errors where possible.

## Code Changes

1. Add a small internal restriction state object in `_kernel.sumstats_munger`.
   It should store `path`, `mode`, `genome_build`, keep identifiers, and
   aggregate counters for rows inspected, usable row identifiers, rows kept, and
   rows removed.
2. Add a pre-parse preparation function:
   - return `None` when `args.sumstats_snps` is absent;
   - load `set[str]` once in `rsid` mode;
   - load packed `set[int]` once in `chr_pos` mode using the streaming
     restriction reader;
   - raise before raw parsing if the keep-set is empty.
3. Add early source-build resolution for every `snp_identifier=chr_pos` run with
   `genome_build=auto`. Read only the mapped raw `CHR` and `POS` columns, call
   shared `resolve_chr_pos_table()`, store inferred build and coordinate basis
   on `args._coordinate_metadata`, and set `args.genome_build` to the concrete
   source build before chunk parsing. Wrap inference failures with guidance to
   pass `--genome-build hg19` or `--genome-build hg38`.
4. Move coordinate normalization/drop behavior into the chunk path for
   `chr_pos`. Use the shared coordinate policy to drop missing or invalid
   coordinates with warnings. If early auto inference found 0-based coordinates,
   shift chunk POS to canonical 1-based coordinates before packed-key matching.
5. Apply restriction filtering after existing chunk QC and coordinate
   normalization, before `dat_list.append(...)`.
6. Remove the post-concat `filter_sumstats_snps()` call from `munge_sumstats()`.
   Replace or delete direct tests that assume post-concat filtering.
7. After chunk parsing, log aggregate keep-list messages equivalent to the
   current logs: identifiers read, usable row identifiers, removed rows, retained
   rows, and source path. Raise the whole-run no-SNPs-remain error only after
   all chunks have been processed.
8. Keep post-concat steps for duplicate rsID removal, `process_n()`, P-to-Z,
   signed-stat direction, liftover, and output writing.
9. Update `docs/current/inference-genome-build.md` so raw sumstats auto
   inference is documented as an early lightweight source-coordinate pass.

## Tests

Add or update tests in `tests/test_sumstats_munger.py`:

- `rsid` keep-list filtering happens inside chunks and preserves raw input order.
- `chr_pos` filtering uses packed keys and does not call `CHR:POS` string-key
  helpers.
- Non-kept rows do not reach post-concat stages such as duplicate-rsID removal.
- `chr_pos + genome_build=auto` resolves before parsing and uses inferred
  build/basis during chunk filtering.
- Empty chunks are skipped, but a whole-run zero-retained result raises the
  keep-list no-SNPs-remain error.
- A keep-list with zero usable identifiers fails before raw sumstats parsing.
- Missing/invalid CHR/POS rows in `chr_pos` mode are dropped with separate
  coordinate warnings.

Run targeted tests first, then the full sumstats munger test file.

## Constraints

Do not measure RSS in automated tests. Avoid unrelated refactors. Preserve
existing user changes in shared coordinate policy files and consume those helpers
instead of replacing them.
