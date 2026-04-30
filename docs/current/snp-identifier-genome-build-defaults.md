# SNP Identifier / Genome Build Default Design

## Problem Statement

`GlobalConfig` has two fields that form an **invariant pair**:

| Field | Purpose |
|---|---|
| `snp_identifier` | Whether SNPs are identified by rsID or CHR:POS |
| `genome_build` | Which coordinate system is assumed for CHR:POS |

The pair is only coherent in certain combinations. The original defaults
(`snp_identifier="chr_pos"`, `genome_build=None`) formed an **invalid pair**
that `__post_init__` silently accepted — not by design, but because validation
did not explicitly forbid it. The process-wide registry worked around this by
registering `GlobalConfig(snp_identifier="rsid")` instead, creating a mismatch
between what the class advertised and what the package actually ran.

`None` was also overloaded: in `rsid` mode it meant "irrelevant"; in `chr_pos`
mode it effectively meant "missing required information" — two semantically
opposite states sharing the same literal.

---

## Resolution

### Default pair: `chr_pos + "auto"`

`chr_pos` is the preferred identifier mode because rsIDs have well-known
anomalies (retired IDs, build-crossing merges, platform-specific artefacts).
`chr_pos` makes the coordinate system explicit.

`"auto"` is the correct companion default because:
- It produces a valid `GlobalConfig()` with zero arguments.
- It expresses "I'll resolve the build from data when a workflow needs it,"
  which is the only honest statement when the user has not yet specified a build.
- It keeps the chr_pos intent without silently committing to `hg19` or `hg38`.
- It remains an opt-in signal: workflows resolve `auto` to a concrete build at
  their respective inference boundaries; they do not auto-infer by default.

### `None` is now unambiguous

After the change, `None` has exactly one meaning: **genome build is irrelevant**,
i.e., `rsid` mode. `chr_pos + None` is an explicit construction error.

| `snp_identifier` | `genome_build` | Meaning | Valid? |
|---|---|---|---|
| `rsid` | `None` | irrelevant | yes (only valid rsid state) |
| `rsid` | `"auto"` | nonsensical | raise |
| `rsid` | `"hg19"` / `"hg38"` | inconsistent | warn + coerce to `None` |
| `chr_pos` | `None` | missing required value | **raise** |
| `chr_pos` | `"auto"` | infer from data | yes (new default) |
| `chr_pos` | `"hg19"` / `"hg38"` | concrete | yes |

`__post_init__` enforces this table. The error message for `chr_pos + None`
includes a fix-it hint: "Pass genome_build='auto' to infer from data."

### CLI defaults remain stricter

The CLI `--genome-build` argument defaults to `None` in every workflow parser.
A command-line user in `chr_pos` mode must type `--genome-build auto` (or a
concrete build) explicitly. This is intentional: CLI users benefit from being
forced to make the choice visible in their shell history and scripts. The API
and CLI can have different ergonomic defaults without being inconsistent,
because they serve different user models.

---

## Process-Wide Registry

After the change:

```python
_GLOBAL_CONFIG: GlobalConfig = GlobalConfig()   # chr_pos + auto
```

`reset_global_config()` also returns `GlobalConfig()`. Workflows that previously
defaulted to `rsid + None` from the registry now default to `chr_pos + auto`.
The shift is intentional: it aligns the registry with the preferred identifier
convention and makes the class default honest.

---

## Regression Provenance Precedence

Regression workflows involve three independent `GlobalConfig` instances:

1. The LD-score manifest snapshot (`LDScoreResult.config_snapshot`)
2. The sumstats sidecar snapshot (`SumstatsTable.config_snapshot`)
3. The runner's live `self.global_config` (from the registry at construction time)

**Rule 1 — pair validation** (`build_dataset:132-137`): fires when **both**
on-disk snapshots are present; compares them against each other. The runner's
live config is not a participant.

**Rule 2 — effective identifier resolution** (`_effective_snp_identifier_mode`):
priority order is LD-score snapshot → sumstats snapshot → runner config → `"rsid"`
fallback. The runner's config is only the fallback for fully un-instrumented inputs.

### Failure modes by provenance state

| LD manifest | Sumstats snapshot | Outcome |
|---|---|---|
| `rsid` snapshot | `rsid` snapshot | Pair check passes; merge uses rsid. Runner's chr_pos+auto default is ignored. |
| `rsid` snapshot | `chr_pos` snapshot | `ConfigMismatchError("snp_identifier mismatch ... 'rsid' vs 'chr_pos'")`. Clear, actionable. |
| `rsid` snapshot | none | Pair check skipped. Mode resolution returns rsid from manifest. Merge proceeds in rsid mode. |
| none | none | Pair check skipped. Mode resolution falls through to runner → chr_pos. If sumstats lack CHR/POS columns: column inference raises. If empty merge: runtime error with active config in message. |

The third row (rsid manifest, no sumstats provenance) shows that on-disk
provenance always dominates: the registry's chr_pos+auto default is effectively
transparent to existing instrumented artifacts.

---

## Empty-Merge Diagnostic

When the merge yields no overlapping SNPs (scenario D above), the error message
includes the active `GlobalConfig`:

```
No overlapping chr_pos SNPs remain after merging sumstats 'my_trait.sumstats.gz'
with 45321 LD-score rows. Check that snp_identifier and genome_build match.
Active config: GlobalConfig(snp_identifier='chr_pos', genome_build='auto', ...).
```

This surfaces the runner's live config as the first diagnostic breadcrumb for a
user whose artifacts have no provenance.

---

## Affected Call Sites

All `GlobalConfig` constructor calls that omit `genome_build` while using
`snp_identifier="chr_pos"` (explicitly or by default) must be updated to pass
`genome_build="auto"` or a concrete build.

| Site | Change |
|---|---|
| `config.py` field defaults | `genome_build = "auto"` |
| `config.py __post_init__` | Add explicit raise for `chr_pos + None` |
| `config.py _GLOBAL_CONFIG` | `GlobalConfig()` (simplification) |
| `config.py reset_global_config()` | `GlobalConfig()` (simplification) |
| `ref_panel_builder.py config_from_args` | Add `genome_build="auto"` (build-ref-panel ignores it, but construction must be valid) |
| `sumstats_munger.py _effective_sumstats_config` | Guard against `coordinate_metadata["genome_build"] = None` overriding a valid config.genome_build |
| `regression_runner.py _global_config_from_manifest` | Fallback genome_build for chr_pos in recovered manifests |
| `regression_runner.py build_dataset` | Append `self.global_config!r` to empty-merge error |

Call sites in `ldscore_calculator.py` already resolve genome_build to a concrete
value or `"auto"` before constructing `GlobalConfig`, so they do not need changes.
Call sites in `sumstats_munger.py` that raise explicitly on `chr_pos + None`
remain correct and consistent.
