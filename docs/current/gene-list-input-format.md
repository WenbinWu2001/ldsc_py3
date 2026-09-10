# Gene-list query input

Last updated on: 2026-09-10

`ldsc ldscore` accepts one or more focal gene lists through
`--query-annot-gene-list-sources`. Direct mode requires an explicit coordinate
catalog; indexed mode uses only the catalog embedded in the index. There is no
installed-package catalog fallback or catalog overlay.

## The four query-input modes

The modes are mutually exclusive. `--gene-ldscore-index-dir` selects an indexed
backend but still requires focal gene lists.

| Mode | Required inputs | Mode-specific optional inputs | Forbidden combinations |
| --- | --- | --- | --- |
| Live gene list | `--query-annot-gene-list-sources`, `--gene-coordinate-file`, explicit `--padding-bp`, live baseline/reference | `--control-gene-list-file`, `--gene-exclude-regions`, `--gene-list-resolution-policy` | Index, BED query, prebuilt query |
| Indexed gene list | `--query-annot-gene-list-sources`, `--gene-ldscore-index-dir` | `--control-gene-list-file`, `--gene-list-resolution-policy` | Live catalog/baseline/reference/build/padding/exclusion overrides, BED query, prebuilt query |
| Prebuilt query annotations | `--query-annot-sources`, live baseline/reference | Ordinary LD-score controls | Gene-list, coordinate/control-gene/index/BED inputs and explicit padding |
| Live query BED | `--query-annot-bed-sources`, live baseline/reference | `--padding-bp` | Gene-list/catalog/control-gene/exclusion/index/prebuilt query inputs |

Baseline-only/no-query calculation remains valid.

Padding is deliberately mode-specific:

| Mode | Omitted | Explicit `0` | Explicit positive value |
| --- | --- | --- | --- |
| Live gene list | Error: choose deliberately | Gene bodies | Padded genes |
| Live BED | Effective `0` | Unpadded BED | Padded BED |
| Prebuilt annotation or no query | Accepted | Rejected | Rejected |
| Indexed gene list | Required omission; inherit index | Rejected | Rejected |

This hard-stop behavior is intentional for SLURM workflows: a successful job
must not depend on someone noticing a warning in a log file.

## Gene-list source format

Each source is a headerless plain-text or gzip-compressed file containing one
identifier per nonblank line. Whitespace around the token is removed. Blank
lines are ignored but physical line numbers are preserved. There is no comment
syntax: `#TP53` is an identifier. A nonblank row containing tab-separated extra
fields is malformed.

Focal source tokens may be exact paths or ordinary glob patterns. Glob matches
are expanded in lexical path order so source ordinals and audit order are
reproducible; a pattern that matches no file is a batched Gate A source error.
`--control-gene-list-file` accepts exactly one literal file and never a glob.
Both interfaces expand `~` and environment variables before resolution.

Matching is exact and case-sensitive. A token may be an authoritative
`gene_id` or an exact `gene_name`; IDs are strongly preferred because names can
be shared. LDSC does not strip Ensembl versions, infer synonyms, change case,
or choose among ambiguous matches. Repeated rows and different aliases that
resolve to the same canonical ID select the gene only once per source.

The query name is the source basename after removing optional final `.gz` and
then at most one final `.txt`, `.tsv`, or `.list`, case-insensitively. Query
names must be unique and must not collide with baseline columns or the reserved
control name `gene_control`.

## Coordinate catalog

`--gene-coordinate-file` is a headered TSV or TSV.GZ with required columns in
any order:

| Column | Meaning |
| --- | --- |
| `gene_id` | Nonempty authoritative identity. IDs must be unique. |
| `gene_name` | Optional exact alias; blank creates no alias. |
| `chrom` | Autosomal chromosome 1–22; `chr1` and `1` normalize identically. |
| `start` | One-based inclusive integer start, at least 1. |
| `end` | One-based inclusive integer end, at least `start`. |
| `genome_build` | One consistent build: hg19/GRCh37 or hg38/GRCh38. |

Blank and `#` comment lines are ignored in the catalog, while physical catalog
line numbers remain available in diagnostics. Extra columns are ignored
scientifically. Coordinates convert exactly once to the internal BED interval:

```text
[start - 1, end)
```

Only the start is decremented. Padding expands that internal interval on both
sides and clips the start at zero. MHC gene exclusion, when requested, tests
the unpadded interval before padding.

The catalog build must agree with explicit `--genome-build`. Under `auto`, its
build is combined with available baseline/reference evidence, and conflicting
evidence stops the run. No implicit liftover occurs. A plausible catalog whose
rows secretly mix builds cannot be detected reliably; generate it from one
authoritative single-build source.

## Strict and exploratory resolution

The default is:

```text
--gene-list-resolution-policy strict
```

Strict Gate A scans all focal lists and the optional control together and then
stops if any submitted row cannot supply one valid interval. It writes the
complete audit/summary before doing substantial LD-score work.

For a preliminary screen of many pathways, explicitly choose:

```text
--gene-list-resolution-policy resolved-only
```

This uses the resolved subset and records the exact effect in diagnostics,
metadata, logs, and a bounded console warning. It does not relax malformed or
unreadable inputs, duplicate query names, structural catalog/build failures,
corrupt indexes, or an unusable requested control. A focal query with no usable
genes, no annotation SNPs, or zero-variance LD scores is skipped while usable
siblings continue. A gene with zero retained-SNP support is nonfatal and
audited separately from failed resolution. A successful CLI run also prints one
bounded Gate B notice when zero-support genes or skipped/warning queries occur;
the detailed rows remain in the diagnostic files. The Python API returns these
statuses and diagnostic paths without writing directly to its caller's console.

## Direct query chromosome scope

For both direct PLINK and parquet-R² runs:

- `@` requires all autosomes 1–22. Missing, unreadable, or invalid required chromosome inputs fail the run.
- Ordinary globs select their actual matches; validated file contents establish the chromosome set. A filename such as `*.22.annot.gz` does not prove chr22-only contents. Additional chromosomes in matched files participate in scope.
- Baseline and reference inputs must have exactly the same validated chromosome set.
- Every selected focal/control gene must lie within that shared scope. A pathway does not need genes on every chromosome in scope.

Quote CLI glob patterns, for example `"annotations/*.22.annot.gz"`. Users own glob selection: a missing file may be undetectable if it disappears from the matches and all remaining required input groups consistently cover the same subset. `@` explicitly requires complete-autosomal coverage and cannot be paired with a chromosome-subset group.

Each focal pathway and required control must lie wholly within scope after unique gene resolution and explicit gene-region exclusions. They need not touch every covered chromosome. Incomplete coverage fails the whole batch under both resolution policies, without truncation or automatic coverage-based skipping. Prebuilt query SNP rows and BED regions must likewise lie within scope. Empty focal selections and valid zero-support/variance focal skips remain supported; unusable controls and all-focal-skipped batches fail.

The following examples describe the chromosome-coverage gate only; a passing example must still satisfy identifier resolution, input integrity, SNP support, and LD-score viability checks.

| Validated baseline/reference inputs | Selected focal/control genes | Coverage outcome |
| --- | --- | --- |
| Both contain only chr22; no `@` declaration | All genes on chr22 | Pass. A direct chromosome-only run is allowed without a new flag. |
| Both contain chr21–22; no `@` declaration | All genes on chr22 | Pass. Pathway chromosomes may be a subset of the shared scope. |
| Both contain chr21–22; no `@` declaration | Any selected gene on chr1 | Fail the entire batch under either resolution policy. |
| Baseline contains chr22; reference contains chr21–22 | All genes on chr22 | Fail: baseline/reference chromosome sets do not match exactly. |
| Any required suite uses `@`, but another required input group contains only chr22 | All genes on chr22 | Fail: `@` requires autosomes 1–22, regardless of pathway contents. |

The log names the chromosomes resolved and entering analysis. `diagnostics/chromosome_scope.json` records the scope and glob-selection caveat; `diagnostics/input_issues.tsv` records input failures. Gene coverage counts and affected genes extend the existing audit/summary. Unknown support is never zero. See [Gene-list diagnostics and repair](gene-list-diagnostics-and-repair.md#chromosome-scope-and-pathway-coverage).

## Direct example

```bash
ldsc ldscore \
  --output-dir results/gene-list-ldscores \
  --baseline-annot-sources "annotations/baseline.@.annot.gz" \
  --query-annot-gene-list-sources "gene-lists/*.txt.gz" \
  --gene-coordinate-file annotations/gene-coordinates.hg38.tsv.gz \
  --padding-bp 100000 \
  --control-gene-list-file gene-lists/assay-background.txt \
  --gene-exclude-regions mhc \
  --gene-list-resolution-policy strict \
  --r2-dir ref-panel/hg38 \
  --snp-identifier chr_pos_allele_aware \
  --genome-build auto \
  --ld-wind-cm 1.0
```

The coordinate catalog controls both focal and control resolution. It replaces,
rather than augments, any other gene universe.

## Indexed example

```bash
ldsc ldscore \
  --output-dir results/gene-list-ldscores \
  --gene-ldscore-index-dir indexes/1000G-EUR-baseline \
  --query-annot-gene-list-sources "gene-lists/*.txt.gz" \
  --control-gene-list-file gene-lists/assay-background.txt \
  --gene-list-resolution-policy strict
```

Do not pass a live catalog, padding, exclusion, baseline/reference, identity,
or build override. The index's embedded catalog and scientific configuration
are immutable authorities. Old indexes with the prior catalog schema must be
rebuilt.

## Diagnostics

Gene-list runs write `diagnostics/gene_list_audit.tsv.gz` and
`diagnostics/gene_list_resolution_summary.tsv`. Runs that pass Gate A also
write one focal row per source in
`diagnostics/query_annotation_status.tsv`. Strict failures publish only the
available diagnostics and log, not root scientific artifacts.

See [Gene-list diagnostics and repair](gene-list-diagnostics-and-repair.md) for
all columns, reasons, count identities, filters, and a list-first curation
workflow. The related index artifact/build contract is documented in
[Exact gene LD-score indexes](gene-ldscore-index.md).
