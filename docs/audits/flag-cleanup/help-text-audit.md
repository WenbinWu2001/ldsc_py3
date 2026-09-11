# CLI help text and ordering audit

Last updated on: 2026-09-11

This is a read-only audit of the current CLI and evidence for agreeing on help-writing guidance. No flag descriptions, parser order, defaults, validation, or scientific behavior were changed during this audit. The user has since approved the [CLI help guidelines](../../current/cli-help-guidelines.md); the findings below remain the pre-cleanup snapshot.

## Coverage

All 13 real `python -m ldsc COMMAND --help` invocations succeeded. The inventory counts visible command-option entries, excluding built-in help and suppressed options; a flag shared by several commands counts once per command. Every current help page has one undifferentiated options section.

| Command | Visible flag entries | Entries with no description |
| --- | ---: | ---: |
| `annotate` | 12 | 0 |
| `ldscore` | 33 | 0 |
| `build-ref-panel` | 18 | 1 |
| `build-gene-ldscore-index` | 20 | 14 |
| `convert-ldsc2-ldscores` | 8 | 8 |
| `munge-sumstats` | 39 | 0 |
| `h2` | 15 | 3 |
| `partitioned-h2` | 17 | 3 |
| `quantile-h2` | 18 | 18 |
| `rg` | 18 | 3 |
| `query-r2` | 7 | 0 |
| `convert-h2-scale` | 7 | 0 |
| `plot` | 3 | 0 |
| **Total** | **215** | **50** |

The command registry and construction paths are in [`cli.build_parser` and `cli.main`](../../../src/ldsc/cli.py). Temporary captured help pages and the full parsed inventory are in `/tmp/ldsc-help-guidance-audit/` for this local audit.

## Findings affecting the guidance

- Missing explanations are widespread, including every option in the converter and quantile commands. Having a description is a requirement for every visible option, not just complex options.
- Internal terminology obscures user actions: “in-memory projection,” “retained reference-panel universe A',” “active focal query columns,” and “stale owned siblings” need explanations in terms of SNP selection, inputs, outputs, and resource use. See [`ldscore_calculator.build_parser`](../../../src/ldsc/ldscore_calculator.py).
- Tautological descriptions such as munging's “Chunksize.” do not explain the unit or effect. Column overrides should explain automatic detection and what supplying an override changes. See [`sumstats_munger.build_parser`](../../../src/ldsc/sumstats_munger.py).
- A parser default of `None` is not a user-facing omission rule. `--padding-bp` is zero when omitted for BED input but must be supplied explicitly for direct gene-list input; `--n-min` uses a calculated threshold when omitted or zero. Required, inferred, calculated, disabled, inherited, and mode-dependent defaults need distinct wording.
- Some descriptions need factual corrections as well as clearer language. LD-score `--r2-dir` calls metadata files optional, although [the current reference-panel contract](../../current/path-specification.md) requires the chromosome metadata. `--export-ref-metadata` describes files beside the output, while [`_write_one_ref_metadata_sidecar`](../../../src/ldsc/ldscore_calculator.py) writes `ref_metadata/chrN_meta.tsv.gz` below the output directory.
- “LD window size in SNPs” omits the reference population and whether the value is a full width or distance on either side. [`build_window_coordinates` and `get_block_lefts`](../../../src/ldsc/_kernel/ldscore.py) use retained reference-panel SNP positions and a maximum distance from the focal SNP. Help should specify units and the population being counted.
- Logging descriptions need level meanings, the default `INFO`, and the destination being controlled. Ordinary progress and debug records go to the workflow log; the CLI console normally reports errors and selected explicit notices. Log lifecycle records are retained independently of the requested threshold. See [workflow logging, “Console vs File Routing”](../../current/workflow-logging.md#console-vs-file-routing) and [`_WorkflowLoggingContext`](../../../src/ldsc/_logging.py).
- Current order frequently reflects code assembly instead of user needs: regression appends its required summary-statistics input after logging and estimator controls; partitioned regression starts with batch size; munging places identity mode last. LD-score reference inputs appear after gene-specific controls. Required/common inputs need prominent placement, with related alternatives and dependencies kept together.
- Presentation needs checking through the real entry points. [`cli._copy_actions`](../../../src/ldsc/cli.py) currently copies individual actions without their help groups, so adding groups to a workflow parser alone would not preserve the same presentation everywhere. Several module help pages also show `__main__.py` in the usage line rather than the user-facing command. Guidance should cover readable command names, argument placeholders, and consistent entry-point presentation.

## Guidance decisions approved

1. Descriptions include purpose; input syntax or units when relevant; default value or behavior when omitted; and mode restrictions or dependencies. The user explicitly requires all flag dependencies, including mutual exclusion and flags that must be supplied together, to be stated.
2. Ordering prioritizes mode/use-case grouping, then practical importance, then no-default options before defaulted options within comparable importance in the same group. This keeps uncommon overrides near the end while retaining related options together.

The approved spelling decisions remain in the [current flag map](../../current/legacy-cli-flag-map.md), including retained `--threads` and `--chunksize`, `--input-format`, and automatic frequency preservation. This guidance discussion does not reopen those decisions.
