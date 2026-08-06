# Legacy Munged Sumstats Compatibility

Last updated on: 2026-08-04

This document defines the public compatibility boundary for munged summary statistics written by the legacy LDSC Python 2 implementation and consumed by LDSC3 regression workflows.

## Scope

The compatibility input is a legacy whitespace-delimited `.sumstats` or `.sumstats.gz` artifact with no LDSC3 identity metadata. Footerless Parquet is not part of the LDSC2 compatibility contract.

Legacy munged sumstats are treated as rsID-based source inputs. Their `SNP` values are compatibility lookup keys, not assertions of the canonical LDSC3 regression identity. The canonical LDSC3 LD-score panel remains authoritative for its recorded SNP identifier mode and for panel-owned identity metadata.

This compatibility boundary applies when legacy munged sumstats are consumed by `h2`, `partitioned-h2`, or `rg`. It does not make legacy LD-score fragments part of the public regression input contract; regression continues to consume a canonical LDSC3 LD-score result directory.

Legacy coordinate columns, when present, do not participate in identity matching. LDSC3 does not infer, declare, compare, or liftover a genome build for a legacy sumstats input. After compatibility matching, the LD-score panel supplies the canonical `CHR` and `POS` values and the regression runs under the panel's recorded SNP identifier mode. Whether alleles are also projected from the panel depends on whether that panel uses an allele-aware identifier mode.

## Compatibility Matching

For each legacy sumstats input, LDSC3 performs this boundary conversion before normal regression dataset assembly:

1. Require the legacy artifact to contain `SNP`, `A1`, `A2`, `Z`, and `N`. Allele-less legacy artifacts are unsupported for all three regression commands, including `rg`. Require the canonical LD-score panel to contain `SNP`; require panel `A1` and `A2` only for an allele-aware panel.
2. Normalize allele letter case and apply the LDSC3 drop-all duplicate policy to the source. When an rsID occurs more than once in one legacy input, drop every row in that source duplicate cluster before panel matching.
3. For an allele-aware panel, find all panel rows with the source rsID, filter those candidates by allele compatibility, and accept the source row only when exactly one compatible panel row remains. A panel rsID need not be unique before allele filtering because panel alleles may safely disambiguate candidates.
4. For an allele-unaware panel, accept a source row only when its rsID identifies exactly one panel row. Multiple panel rows with that rsID cannot be disambiguated and are dropped.
5. With an allele-aware panel, accept direct, strand-complement, swapped, and swapped-strand-complement allele relationships. Drop strand-ambiguous A/T and C/G pairs, missing or invalid alleles, and relationships incompatible with every panel candidate. Orient the legacy association statistic to the accepted panel row's `A1`: keep `Z` for direct or strand-complement matches and negate `Z` for swapped or swapped-strand-complement matches.
6. With an allele-unaware panel, retain the source `A1`, `A2`, `Z`, and `FRQ` orientation. No panel allele exists to justify changing them. This is sufficient for `h2`, which consumes `Z` through `Z` squared. For `rg`, project each trait independently to panel SNP rows, then harmonize the two traits to one another; negate `Z` and transform a valid `FRQ` to `1 - FRQ` for the swapped trait as required.
7. In both cases, replace the matched row's `CHR`, `POS`, and `SNP` with the panel values and run regression using the panel's recorded identity mode. With an allele-aware panel, also replace `A1` and `A2` with the panel alleles.

If the legacy row carries `FRQ`, it remains a sumstats-owned value and is never imputed from panel frequency metadata. For a row whose allele orientation is swapped against an allele-aware panel, or for a trait swapped during pairwise `rg` harmonization under an allele-unaware panel, transform the legacy value to `1 - FRQ` so it remains the frequency of the oriented `A1`; otherwise retain it unchanged. Panel frequency metadata is likewise never imputed from sumstats.

`FRQ` is not consumed by regression. A missing, nonnumeric, or out-of-range legacy `FRQ` therefore does not drop an otherwise compatible SNP: the compatibility layer preserves it as missing and logs an aggregate warning. The `1 - FRQ` transformation applies only to valid numeric values in `[0, 1]`.

Source duplicate handling and panel candidate handling do not conflict. Source duplicates describe multiple association rows for one rsID within a single legacy trait and are dropped as an unverifiable source cluster. Multiple panel candidates describe possible canonical targets for one otherwise unique source row. They may be retained only when an allele-aware panel allows allele comparison to select exactly one target; an allele-unaware panel cannot resolve them.

## Drop and Failure Policy

Row-level failures are dropped and audited without aborting an otherwise usable input:

| Condition | Outcome / reason |
| --- | --- |
| duplicated source rsID | Drop the complete source cluster / `duplicate_source_rsid` |
| no panel row with the source rsID | Drop / `missing_panel_rsid` |
| no allele-compatible panel candidate | Drop / `incompatible_alleles` |
| more than one allele-compatible panel candidate | Drop / `ambiguous_panel_mapping` |
| missing, invalid, or strand-ambiguous source alleles | Drop with the corresponding allele reason |

Missing required artifact columns, missing required panel identity columns, or zero retained compatibility rows are file-level failures and abort the affected regression. The compatibility layer introduces no new retention-percentage threshold; ordinary regression SNP-count diagnostics continue to apply.

When an output directory is supplied and at least one legacy input is present, LDSC3 writes `diagnostics/dropped_snps/legacy_sumstats.tsv.gz`, including an empty header-only artifact when no rows were dropped. The audit table contains `trait_name`, `source_path`, `SNP`, `A1`, `A2`, `reason`, and `panel_candidate_count`. Workflow logs report aggregate counts by reason at `INFO` or `WARNING`; example SNP identifiers are `DEBUG`-only. Stdout-only runs report the aggregate counts without writing a sidecar.

## Workflow and API Behavior

Compatibility projection is automatic in both the public CLI regression workflows and the public Python regression workflow. `load_sumstats()` identifies `.sumstats` and `.sumstats.gz` as legacy LDSC2 inputs; users do not call a separate converter. `RegressionRunner` projects a marked legacy table onto the supplied canonical LD-score panel before ordinary dataset assembly.

An `rg` invocation may mix current self-describing Parquet sumstats and legacy text sumstats. Each legacy input is independently projected onto the panel identity. Current Parquet inputs retain their recorded provenance and are never reinterpreted by the compatibility layer. After projection, all inputs pass through the same ordinary compatibility checks against the canonical panel.

## Related Contracts

- [`config-design.md`](config-design.md) defines configuration provenance and compatibility checks.
- [`column-schema.md`](column-schema.md) defines canonical LDSC3 columns and accepted input aliases.
- [`io-argument-inventory.md`](io-argument-inventory.md) defines public command inputs and outputs.
