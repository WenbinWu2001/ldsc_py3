# Pilot Reference — `munge-sumstats` self-contained messages + reference section

This file is the **gold-standard template** for the logging-refinement effort.
The follow-up session must:

1. Make every run-aborting message **self-contained** (what & where + most-likely
   cause + top remedy) — copy the in-code message shape in Part C exactly.
2. Add a `docs/troubleshooting.md` entry **only** for the few errors with many
   branching causes — copy the section/heading format in Part B exactly.

No error-ID codes. The reference doc is organized by command, with plain-language
headings; messages link to it by the heading's slug, and only when the error
genuinely has more causes than fit in one line.

All error sites, line numbers, and messages below were read from the live
codebase on 2026-06-04. Re-confirm line numbers before editing — they drift.

---

## Part A — `docs/troubleshooting.md` file header (write once)

```markdown
# Troubleshooting

This reference explains `ldsc` errors that can **abort a run** and have more than
one likely cause. It is organized by command. Each entry lists the likely causes
(ranked most-probable first), how to confirm each, and how to fix it.

Most errors are self-explanatory from their terminal message alone and are not
repeated here. When a message says `... see docs/troubleshooting.md#<section>`,
that slug is a heading below — jump to it.
```

---

## Part B — the `munge-sumstats` reference section (multi-cause errors only)

> Catalog section heading: `## munge-sumstats`
> Only the three errors below get an entry — they each have several distinct
> causes. Every other munge error (missing flag, option/mode conflict, missing
> `pyarrow`, internal N-inference bug) is fully handled by its self-contained
> message in Part C and gets **no** entry here.

```markdown
## munge-sumstats

### munge-sumstats: could not map a required column

**Raised by:** `_kernel/sumstats_munger` column resolution · `SumstatsTable.validate()`
· **Exception:** `LDSCInputError`
**Symptom:** `munge-sumstats could not map the required column '<field>' from the header...`

**Likely causes & how to check** (most probable first):

| # | Likely cause | How to check |
|---|--------------|--------------|
| 1 | The header uses a name the auto-mapper doesn't recognize | `zcat <file> \| head -1` (or `head -1`); compare against the recognized aliases in `column_inference.py` |
| 2 | The column exists under a synonym you must declare | Re-run with an explicit hint, e.g. `--snp-col MarkerName --a1-col Allele1` |
| 3 | Wrong delimiter, so the whole header parsed as one column | `head -1 <file> \| cat -A` — look for one field with embedded tabs/commas |
| 4 | Wrong `--format`, so expected columns differ | Confirm the format flag matches the file (e.g. drop `--daner-new` for a non-PGC file) |

**Remedies:**

1. Pass explicit column hints for the unmapped field(s): `--<field>-col <name>`.
2. Rename the columns to recognized names, or add the alias to `column_inference.py`
   if broadly useful.
3. Verify the delimiter is whitespace/tab as expected for `.sumstats`/`.txt` inputs.

### munge-sumstats: no SNPs remain after filtering

**Raised by:** `_kernel/sumstats_munger` (post-filter and keep-list paths)
· **Exception:** `LDSCInputError`
**Symptom:** `munge-sumstats removed every SNP...` / `...no SNPs remain after applying --sumstats-snps-file`

**Likely causes & how to check** (most probable first):

| # | Likely cause | How to check |
|---|--------------|--------------|
| 1 | `--sumstats-snps-file` keep-list uses a different SNP-id space (rsID vs chr:pos) | Compare the first IDs of each file; confirm both use `snp_identifier=<mode>` |
| 2 | Keep-list and sumstats are on different genome builds | Compare the `genome_build` in the error line against the keep-list build |
| 3 | INFO / MAF / N thresholds removed every row | Re-run with relaxed `--info-min` / `--maf-min`; inspect the dropped-SNP sidecar |
| 4 | Input is effectively empty after a delimiter/format mis-parse | `zcat <file> \| wc -l` |

**Remedies:**

1. Regenerate the keep-list in the same `snp_identifier` mode and build as the sumstats.
2. Relax quality filters and re-run; review the dropped-SNP audit sidecar to see
   which filter removed the rows.

### munge-sumstats: curated artifact is malformed or outdated

**Raised by:** `sumstats_munger.load_sumstats()` and its metadata helpers
· **Exception:** `LDSCInputError`
**Symptom:** `Cannot load curated sumstats at '<path>': ...` (sidecar missing / old
schema / bad provenance / missing A1-A2 / duplicate identity rows)

**Likely causes & how to check** (most probable first):

| # | Likely cause | How to check |
|---|--------------|--------------|
| 1 | Artifact predates the current schema (no/old `metadata.json` sidecar) | Check for `metadata.json` beside the artifact; inspect its `schema_version` |
| 2 | The `metadata.json` sidecar is missing or unreadable | `python -m json.tool metadata.json` |
| 3 | Identity provenance in the sidecar is invalid/corrupt | Inspect the sidecar's identity fields against the current contract |
| 4 | An allele-aware `snp_identifier` artifact lacks A1/A2 columns | `python -c "import pandas; print(pandas.read_parquet('<f>').columns)"` |
| 5 | Duplicate/invalid SNP-identity rows survived in the artifact | Re-munge from raw input; the loader reports the dropped-row reasons |

**Remedies:**

1. Re-run `ldsc munge-sumstats` from the **raw** GWAS file to regenerate the
   artifact with the current schema.
2. Do not hand-edit curated `.sumstats`/`.parquet` artifacts; treat them as outputs.
```

> Note: causes 1–5 of the last entry were previously **four+ identical**
> `REGENERATE_ARTIFACT_MESSAGE` raises. In code, give each call site a message that
> names *which* check failed (see Part C), and link all of them to this one
> heading. Retire `REGENERATE_ARTIFACT_MESSAGE` as a shared constant or make it a
> helper that interpolates the failing check.

---

## Part C — Pilot in-code rewrites (message text + exception class; behavior unchanged)

Every message becomes self-contained. Only the three multi-cause errors carry a
`docs/troubleshooting.md#...` link; the rest say everything inline. **Exception
class also changes** from bare `ValueError` to the matching `errors.py` subclass so
the CLI boundary classification is intentional, not reliant on the blanket
`ValueError` catch.

| File:line | Current | Replacement (self-contained message + class) | Links to ref? |
|-----------|---------|----------------------------------------------|---------------|
| `sumstats_munger.py:442` | `raise ValueError("MungeConfig.raw_sumstats_file is required to run summary-statistics munging.")` | `raise LDSCUserError("No input summary statistics given to munge-sumstats. Most likely --sumstats was omitted. Pass --sumstats <file> (and --out <stem>).")` | No (1 cause) |
| `sumstats_munger.py:178` | `raise ValueError(f"SumstatsTable is missing required columns: {sorted(missing)}")` | `raise LDSCInputError(f"munge-sumstats could not map required column(s) {sorted(missing)} from the input header. Most likely the file uses unrecognized column names. Pass explicit hints (e.g. --snp-col, --a1-col) or rename the columns. Other causes & fixes: docs/troubleshooting.md#munge-sumstats-could-not-map-a-required-column")` | **Yes** |
| `_kernel/sumstats_munger.py:1276` | `raise ValueError('After applying filters, no SNPs remain.')` | `raise LDSCInputError("munge-sumstats removed every SNP during quality filtering. Most likely --info-min/--maf-min are too strict for this file. Relax the thresholds and check the dropped-SNP sidecar. Other causes & fixes: docs/troubleshooting.md#munge-sumstats-no-snps-remain-after-filtering")` | **Yes** |
| `sumstats_munger.py:324` (and 1707/1709) | `raise ValueError(REGENERATE_ARTIFACT_MESSAGE)` | `raise LDSCInputError(f"Cannot load curated sumstats at '{resolved}': its metadata.json sidecar is missing, so its identity semantics are unknown. Most likely it predates the current schema. Re-run munge-sumstats from the raw GWAS file. Other causes & fixes: docs/troubleshooting.md#munge-sumstats-curated-artifact-is-malformed-or-outdated")` | **Yes** |
| `sumstats_munger.py:752` | `raise ValueError("--source-genome-build is only valid for chr_pos-family snp_identifier modes.")` | `raise LDSCUsageError("--source-genome-build applies only to chr_pos-family --snp-identifier modes, but the current mode is rsID-based. Either drop the build flag or switch to a chr_pos identifier.")` | No (2 causes inline) |
| `_kernel/sumstats_munger.py:682` | `raise ValueError('Cannot determine N. This message indicates a bug.\n' ...)` | `raise LDSCInternalError("munge-sumstats could not derive a sample size (N) and reached a state that should be unreachable. Re-run with --debug and report the traceback.")` | No (internal) |
| `_kernel/sumstats_munger.py:1532`→`sumstats_munger.py:1532` | `raise LDSCDependencyError("Writing sumstats parquet artifacts requires pyarrow.") from exc` | `raise LDSCDependencyError("munge-sumstats needs the 'pyarrow' package to write Parquet output, but it is not installed. Install it (pip install pyarrow) or choose --output-format tsv.gz.") from exc` | No (1 cause) |

Notes for the rewriter:
- Keep the *single* most-likely cause and *single* top remedy inline. For the three
  linked errors, the full ranked list lives in the Part B section.
- Decide "link or not" by cause count: 1–2 causes → say them inline, no link;
  ≥3 distinct causes → top cause inline + link.
- The link slug must match the GitHub auto-slug of the `###` heading exactly
  (lowercase, spaces → hyphens, punctuation dropped). The three slugs are:
  `munge-sumstats-could-not-map-a-required-column`,
  `munge-sumstats-no-snps-remain-after-filtering`,
  `munge-sumstats-curated-artifact-is-malformed-or-outdated`.
- Preserve every existing `from e` / `from None` (e.g. the parquet dependency and
  the chained metadata parse at line 1685) — never drop the chain.
