# CLI flag descriptions and ordering

Last updated on: 2026-09-11

These are the agreed rules for writing flag descriptions and organizing `--help` across all public `ldsc` commands. A user should be able to identify the inputs they need, understand what an option changes, and choose a valid combination without knowing the package implementation.

The [help audit](../audits/flag-cleanup/help-text-audit.md) records the starting problems. The [IO argument inventory](io-argument-inventory.md) describes the current options, and the [legacy flag map](legacy-cli-flag-map.md) records approved spelling and behavior decisions. These rules have been applied across the public command parsers; see the [implementation record](../audits/flag-cleanup/help-polish-followup.md).

## Required description content

Every visible flag must have a useful description. Include these four elements, incorporating input details and dependencies where applicable:

| Element | What the description must explain |
| --- | --- |
| Purpose | What supplying the flag changes for the user: selecting SNPs, filtering rows, choosing a model, extending intervals, or writing files. |
| Input | Whether the value is a file, directory, column name, choice, count, distance, or proportion; include relevant syntax, units, ranges, and the population being counted. |
| Default or omission | The default value and its effect, or what happens when the flag is omitted: required input, automatic detection, a calculated value, inherited settings, no restriction, or a disabled feature. |
| Dependencies and scope | Every applicable dependency on other flags, including required companions, mutual exclusion, mode restrictions, conditional requirements, and precedence. |

Use this as a writing pattern, not a mandatory sequence of labeled fields:

> What the option does, including its input or units. Essential restrictions or interactions. Default value and behavior when omitted.

Required inputs should say “Required,” including the applicable mode when the requirement is conditional. They do not need a fabricated default sentence.

## Writing rules

1. **Lead with the user-visible purpose.** Use concrete verbs such as “read,” “select,” “keep,” “exclude,” “extend,” and “write.” Replace implementation terms such as “in-memory projection,” “artifact provenance,” and “owned siblings” with the relevant action or output. Scientific terms such as SNP, LD, and MAF are appropriate when they convey the actual analysis.
2. **Specify the value precisely.** Distinguish a file from a directory or filename prefix, a column name from a constant, and a proportion from a percentage. State whether a count refers to reference-panel SNPs, input rows, individuals, or annotations. For windows, distinguish a distance on either side from a total width. For intervals and thresholds, state coordinate conventions or inclusive/exclusive boundaries when needed to choose a valid value.
3. **Explain omission in user terms.** Prefer `Default: 1000`, `Default: off; no additional files are written`, or `If omitted, detect the column from its header`. Never use `Default: None` as the whole explanation. Explain what `auto` detects and what a calculated or inherited default depends on. Describe special values such as zero or negative counts when they change behavior. State mode-dependent defaults separately.
4. **State all flag dependencies explicitly.** Follow the dependency rules below. A group heading or syntax diagram supplements the description; it does not replace a necessary dependency statement.
5. **Explain non-obvious choices.** A list such as `{strict,resolved-only}` does not explain the difference. Describe the effect of each choice whose meaning is not clear from its name. For logging, describe every level, the default, and the destination controlled.
6. **Keep the text concise and complete.** Usually use two to four short sentences. This is a target, not a limit that permits omitting dependencies or scientific consequences. Logging and options with several conditional cases may need more. Put detailed algorithms and full schemas in documentation, while keeping information needed to choose a valid value in help.
7. **Verify wording against behavior.** Check the parser, configuration, and workflow before stating defaults, file destinations, scope, or precedence. Use consistent wording for shared meanings, while preserving command-specific differences. A wording cleanup must not silently change validation, numerical behavior, accepted inputs, flag names, or defaults. Record behavior disagreements separately when they require a decision.
8. **Explain supported path patterns once per command, with precise cross-references.** At the first flag supporting chromosome suites, explain that quoted `*` patterns match filename text and `@` substitutes chromosome numbers 1-22 by default. Say “Use `*` to select available files for a chromosome subset where allowed.” Later compatible flags may say “Supports `*` and `@`; see `--first-flag` for pattern rules.” Keep exact-one match requirements, unsupported `@`, and command-specific coverage requirements in each affected description. Check the actual reader: a `sources` suffix alone does not establish pattern support. Directory paths do not acquire glob support from file flags.
9. **State cumulative interval expansion.** Padding extends the supplied BED coordinates at both ends. Tell users to set it to 0 when BED intervals already include the intended padding, to avoid double padding. When reconstructing a saved model, preserve its effective intervals and reproduce the original setting.

## Dependencies and interactions

Use exact flag names and make the relationship unambiguous:

- **Required together:** say “Requires `--other-flag`” or “Supply both …”. State the relationship on both affected options so either description is usable on its own.
- **Mutually exclusive:** say “Cannot be combined with …”. Distinguish “choose exactly one” from “choose at most one”; the former also requires a selection.
- **Conditional requirements:** name the triggering mode or value, then the required companion. Avoid calling a flag merely “optional” when it is required in that mode.
- **Mode restrictions:** say where a flag applies and whether another mode rejects it or ignores it. These are different behaviors and must not be conflated.
- **Precedence:** explain which supplied value, input column, or inferred setting wins. Do not describe a fallback as an override.
- **Shared group requirements:** place a concise explanation in the group introduction and retain the essential relationship in each affected flag description. Keep companions and alternatives adjacent.

For example, `--N-cas-col` must name its required `--N-con-col` companion and its incompatibility with `--N-col`. The control-column option must state the reciprocal requirement. These are explicit requirements for the help cleanup, not optional details to omit for brevity. See [`sumstats_munger.build_parser`](../../src/ldsc/sumstats_munger.py) and the [sample-size contract](munge-sumstats.md#sample-size-column-selection).

## Ordering rules

Apply these priorities in order:

1. **Group by task, mode, or use case.** Use visible, descriptive headings. Keep BED inputs, gene-list inputs, direct computation, and indexed computation recognizable as applicable. Shared inputs should have a clear common group.
2. **Order groups by practical necessity and importance.** Put essential inputs and output destinations first, common analysis settings next, specialized workflows afterward, and advanced overrides, performance controls, and logging near the end. Do not let the order of helper functions determine the help order.
3. **Within each group, put required inputs first, then commonly adjusted options.** Conditional requirements belong with the mode that needs them. Preserve the relationship between a primary option and its companion flags.
4. **Use default status as a tie-breaker.** Among similarly important options in the same group, put options needing an explicit value before those with usable defaults. A parser value of `None` can still represent automatic behavior; classify the effective user behavior, not the Python value.
5. **Keep alternatives and dependencies adjacent.** Do not separate paired case/control settings or the three LD-window alternatives merely to sort defaults. An uncommon override stays in its advanced group even if it has no fixed default.

There is no global “all no-default flags first” rule. Grouping and practical importance take precedence. Small commands need only the groups that help users; they do not need empty or artificial sections.

The following sequences are starting points for organizing the two largest relevant workflows, subject to keeping actual dependencies together:

| Command | Group sequence |
| --- | --- |
| `munge-sumstats` | Essential inputs/output; SNP identity and genome build; sample size; filtering; format-specific and column overrides; performance and logging. |
| `ldscore` | Essential inputs/output; direct computation; annotation sources with related BED/gene options together; indexed computation; advanced and resource controls. |

Within these sequences, required identity/build settings and mode-selecting inputs must remain easy to find. An essential selector can appear in the common input group while its less common tuning options stay with the corresponding mode.

Use readable argument placeholders such as `FILE`, `DIR`, `COLUMN`, `BP`, or `N`. Show the actual user-facing command in the usage line. Equivalent entry points must preserve the same groups and descriptions; verify the rendered help rather than assuming that copying parser actions preserves their groups. See [`cli.build_parser`, `cli.main`, and `cli._copy_actions`](../../src/ldsc/cli.py).

## Common mistakes and better descriptions

For path inputs, “accepts globs” alone is incomplete: identify `*`, quoting, match cardinality, and whether `@` is expanded. Do not claim that every suite reader enforces all 22 members; the full-coverage check is workflow-specific. Direct query `ldscore` and gene-index construction enforce it, whereas generic suite expansion can retain existing members. The [path specification](path-specification.md#pattern-support-in-command-help) records these distinctions.

The poor descriptions below are existing text or representative paraphrases. Improved excerpts illustrate a specific rule; a complete flag description must also include any additional dependencies or omission behavior that applies to that command.

| Mistake | Improper description | Better wording or correction |
| --- | --- | --- |
| Internal implementation jargon | “Base pairs to add … before in-memory projection.” | “Extend each BED interval by this many base pairs at both ends, then find SNPs within the expanded interval. Default for BED inputs: 0, using the original intervals.” Add the explicit-value requirement when the same flag supports gene-list input. |
| Repeating the flag name | `--chunksize`: “Chunksize.” | “Read this many input rows per chunk. Default: 1,000,000.” |
| Unspecified logging behavior | `--log-level`: “Logging verbosity.” | Explain all four levels, `INFO` as the default, and the workflow log as the destination; see the full example below. |
| Missing boolean default and output location | “Write a reference-metadata sidecar next to the LD-score output.” | “Write reference-panel SNP metadata to `ref_metadata/chrN_meta.tsv.gz` under the output directory. Applies to PLINK input only. Default: off; these additional files are not written.” |
| Ambiguous window unit | “LD window size in SNPs.” | “LD-window size on either side of each SNP, measured by SNP count in the retained reference panel. No default; specify exactly one of `--ld-wind-snps`, `--ld-wind-kb`, or `--ld-wind-cm`.” This example describes direct `ldscore`. |
| Implementation default instead of omission behavior | `--n-min`: “Default: None.” | “Minimum per-variant sample size. If omitted or zero, use the 90th percentile of N divided by 1.5. Constant sample sizes are not filtered.” |
| Hidden required companion | `--N-cas-col`: “Case-count column.” | “Name of the per-variant case-count column. Requires `--N-con-col`; cannot be combined with `--N-col`.” Also explain inference and selection precedence in the complete description. |
| Fallback incorrectly called an override | `--N`: “Override sample size.” | “Constant sample-size fallback when per-variant N or paired case/control columns are absent. Input columns take precedence; this value takes precedence over `--N-cas`/`--N-con`.” |
| Unqualified mode-dependent default | `--padding-bp`: “Default: 0.” | “For BED inputs, omission means 0. Direct gene-list input requires an explicit value, including 0 to use gene bodies.” Include the flag's purpose and other applicable mode restrictions. |
| Unclear choice meanings | `--gene-list-resolution-policy`: “Gene identifier policy.” | “Choose how to handle rejected gene identifiers: strict stops the run; resolved-only continues with the usable subset. Default: strict.” |
| Incorrect ordering | Logging, batch sizes, or advanced overrides before the main input file. | Start with essential inputs and output destinations, then apply the group and ordering priorities above. |

For a workflow that writes a log, the logging description can read:

> Set the detail recorded in the workflow log. DEBUG includes troubleshooting details. INFO includes progress and summaries. WARNING includes warnings and errors. ERROR includes errors only. Each level also includes more severe messages. Default: INFO. Run headers and completion status are always recorded.

Do not promise that increasing `--log-level` prints progress to the terminal. Ordinary progress/debug messages go to the workflow log; console error reporting and explicit scientific notices have their own routing. See [workflow logging, “Console vs File Routing”](workflow-logging.md#console-vs-file-routing).

The examples are grounded in [`sumstats_munger.build_parser`](../../src/ldsc/sumstats_munger.py), [`ldscore_calculator.build_parser` and `_write_one_ref_metadata_sidecar`](../../src/ldsc/ldscore_calculator.py), [`build_window_coordinates` and `get_block_lefts`](../../src/ldsc/_kernel/ldscore.py), and [`_WorkflowLoggingContext`](../../src/ldsc/_logging.py). Preserve these behavioral distinctions when adapting the examples to other commands.

## Review before accepting help changes

- Every visible option has a description with its purpose, relevant input details, effective default or omission behavior, and all dependencies.
- Conditional requirements, paired flags, mutual exclusion, and precedence agree with validation and workflow behavior.
- Default claims and output locations agree with the implementation; shared flag wording reflects any command-specific differences.
- Groups follow user tasks and importance, with required inputs prominent and uncommon controls near the end.
- Rendered `ldsc COMMAND --help` and `python -m ldsc COMMAND --help` remain readable, show the correct command, and preserve the intended grouping and ordering.
- Help-only changes preserve accepted flags, parser destinations, defaults, choices, dependencies, and scientific behavior. Relevant parser/CLI checks establish this; prose edits alone do not require numerical test suites.

The description contents and ordering precedence are approved. All public command parsers use explicit help groups and descriptions, with shared wrapping in [`CLIHelpFormatter`](../../src/ldsc/_cli_help.py) and consistent log-level wording in [`LOG_LEVEL_HELP`](../../src/ldsc/_logging.py). No guidance decisions remain open.
