# SNP Identifier / Genome Build Defaults

This document records the defaults implemented in the current codebase. It is
descriptive, not a future design proposal.

## Current Defaults

`GlobalConfig` has two shared analysis fields:

| Field | Purpose | Current class default |
|---|---|---|
| `snp_identifier` | Whether SNPs are identified by rsID or normalized `CHR:POS` | `"chr_pos"` |
| `genome_build` | Coordinate system for `chr_pos` interpretation | `None` |

The process-wide Python registry is intentionally stricter for backward
compatibility with existing rsID-oriented LDSC workflows:

```python
_GLOBAL_CONFIG = GlobalConfig(snp_identifier="rsid")
reset_global_config()  # also returns GlobalConfig(snp_identifier="rsid")
```

As a result:

- `GlobalConfig()` constructs `snp_identifier="chr_pos", genome_build=None`.
- `get_global_config()` starts as `snp_identifier="rsid", genome_build=None`.
- `rsid` mode ignores concrete genome-build values after warning and coercing
  them to `None`.
- `rsid` mode rejects `genome_build="auto"` because auto-inference is meaningful
  only for coordinate identifiers.
- `chr_pos + None` is currently allowed at config construction time, but
  workflows that need external coordinate interpretation require a concrete
  build or `"auto"` before they proceed.

## Validity Table

| `snp_identifier` | `genome_build` | Current behavior |
|---|---|---|
| `rsid` | `None` | valid |
| `rsid` | `"auto"` | construction error |
| `rsid` | `"hg19"` / `"hg38"` | warning, then coerced to `None` |
| `chr_pos` | `None` | construction is allowed; coordinate-aware workflows may reject it |
| `chr_pos` | `"auto"` | valid; workflows infer `hg19`/`hg38` when needed |
| `chr_pos` | `"hg19"` / `"hg38"` | valid concrete coordinate mode |

## CLI Defaults

The command-line parsers expose workflow-specific defaults:

| Workflow | Identifier default | Genome-build default | Notes |
|---|---|---|---|
| `ldsc annotate` | `chr_pos` | `None` | `chr_pos` inputs require `--genome-build auto`, `hg19`, or `hg38`. |
| `ldsc ldscore` | `chr_pos` | `None` | `chr_pos` inputs require a resolved build; `rsid` mode ignores build metadata. |
| `ldsc build-ref-panel` | registered config or `--snp-identifier` | source build is separate | `--source-genome-build` is inferred from `.bim` when omitted; `GlobalConfig.genome_build` is ignored. |
| `ldsc munge-sumstats` | `chr_pos` | `hg38` | `--genome-build auto` can infer and normalize coordinate-basis from parsed raw rows. |
| Regression commands | registry at runner construction | registry at runner construction | On-disk provenance from LD-score manifests and sumstats sidecars takes precedence. |

## Regression Provenance Precedence

Regression workflows can see three independent `GlobalConfig` sources:

1. The LD-score manifest snapshot (`LDScoreResult.config_snapshot`)
2. The sumstats sidecar snapshot (`SumstatsTable.config_snapshot`)
3. The runner's live `self.global_config` from the process registry

Compatibility validation runs when both on-disk snapshots are present. The
effective SNP identifier mode is selected in this order:

```text
LD-score snapshot -> sumstats snapshot -> runner config -> "rsid"
```

The runner config is only a fallback for uninstrumented inputs. If both
artifacts carry provenance and disagree on `snp_identifier` or `genome_build`,
`validate_config_compatibility()` raises `ConfigMismatchError`.

## Failure Modes By Provenance State

| LD manifest | Sumstats snapshot | Outcome |
|---|---|---|
| `rsid` snapshot | `rsid` snapshot | Pair check passes; merge uses `SNP`. |
| `rsid` snapshot | `chr_pos` snapshot | `ConfigMismatchError` on `snp_identifier`. |
| `rsid` snapshot | none | Pair check skipped; mode resolution uses `rsid` from the manifest. |
| none | none | Pair check skipped; mode resolution falls through to the runner config and then `rsid` fallback. |

When a merge yields no overlapping SNPs, the current error names the effective
identifier mode and suggests checking `snp_identifier` and `genome_build`.
