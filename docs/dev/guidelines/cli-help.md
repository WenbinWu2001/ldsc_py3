# CLI flag descriptions and ordering

Last updated on: 2026-10-29

These are the agreed rules for writing flag descriptions and organizing `--help` across all public `ldsc` commands. A user should be able to identify the inputs they need, understand what an option changes, and choose a valid combination without knowing the package implementation.

The [IO argument inventory](../../reference/cli.md) describes the current options, and the [legacy flag map](../../reference/legacy-flag-map.md) records approved spelling and behavior decisions. These rules have been applied across the public command parsers.

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
5. **Explain non-obvious choices.** A list such as `{none,mhc}` does not explain which genes are excluded or when exclusion occurs. Describe the effect of each choice whose meaning is not clear from its name. For logging, describe every level, the default, and the destination controlled.
6. **Keep the text concise and complete.** Usually use two to four short sentences. This is a target, not a limit that permits omitting dependencies or scientific consequences. Logging and options with several conditional cases may need more. Put detailed algorithms and full schemas in documentation, while keeping information needed to choose a valid value in help.
7. **Verify wording against behavior.** Check the parser, configuration, and workflow before stating defaults, file destinations, scope, or precedence. Use consistent wording for shared meanings, while preserving command-specific differences. A wording cleanup must not silently change validation, numerical behavior, accepted inputs, flag names, or defaults. Record behavior disagreements separately when they require a decision.
8. **Explain supported path-token syntax once per command, with precise cross-references.**
   At the first flag in help order that supports chromosome suites, explain that quoted `*` patterns match filename text and `@` substitutes chromosome numbers 1-22 by default.
   Say “Use `*` to select available files for a chromosome subset where allowed.”
   Later compatible flags may say “Supports `*` and `@`; see `--first-flag` for pattern rules.”
   A cross-reference names only a flag shown earlier in the help, so reordering the help can move the full explanation to a different flag.
   Keep exact-one match requirements, unsupported `@`, and command-specific coverage requirements in each affected description.
   Check the actual reader: a `sources` suffix alone does not establish pattern support.
   Directory paths do not acquire glob support from file flags.
9. **State cumulative interval expansion.** Padding extends the supplied BED coordinates at both ends. Tell users to set it to 0 when BED intervals already include the intended padding, to avoid double padding. When reconstructing a saved model, preserve its effective intervals and reproduce the original setting.

## Dependencies and interactions

Use exact flag names and make the relationship unambiguous:

- **Required together:** say “Requires `--other-flag`” or “Supply both …”. State the relationship on both affected options so either description is usable on its own.
- **Mutually exclusive:** say “Cannot be combined with …”. Distinguish “choose exactly one” from “choose at most one”; the former also requires a selection.
- **Conditional requirements:** name the triggering mode or value, then the required companion. Avoid calling a flag merely “optional” when it is required in that mode.
- **Mode restrictions:** say where a flag applies and whether another mode rejects it or ignores it. These are different behaviors and must not be conflated.
- **Precedence:** explain which supplied value, input column, or inferred setting wins. Do not describe a fallback as an override.
- **Shared group requirements:** place a concise explanation in the group introduction and retain the essential relationship in each affected flag description. Keep companions and alternatives adjacent.

For example, `--N-cas-col` must name its required `--N-con-col` companion and its incompatibility with `--N-col`. The control-column option must state the reciprocal requirement. These are explicit requirements for the help cleanup, not optional details to omit for brevity. See [`_sumstats.munger.build_parser`](../../../src/ldsc/_sumstats/munger.py) and the [sample-size contract](../../reference/munge-sumstats.md#sample-size-column-selection).

## Ordering rules

Order `--help` as the steps a typical user follows to write the command, matching the [guided tutorial](../../user-guide/guided-tutorial.md).
The help then doubles as a guide to which flags to set.
Important and required flags and sections come first; optional and less essential ones come later.

Each command uses these section headings, in this order:

| Heading | Contents |
| --- | --- |
| `Step N: <decision>` | The decisions a typical run needs, numbered from 1 in tutorial order. Add `(choose one)`, `(choose at most one source type)` or `(optional)` when it applies. |
| `Step N continued: <mode> only` | Settings that apply to only one alternative of step N, placed directly after it. |
| `Optional: <topic>` | Refinements with usable defaults, most commonly adjusted first. |
| `Alternative: <mode>` | A different way to run the command; its description names the steps it replaces and the flags it still needs. |
| `Computing resources` | Workers, batch sizes, chunk sizes and other memory or speed controls. |
| `Logging` | `--log-level`, always the last section. |

Apply these rules within that structure:

1. **Start with SNP identity and genome build when the user chooses them.**
   List `--snp-identifier` before the genome-build flags, as in the tutorial's shared configuration.
   A command whose essential inputs come with an identity that has a usable default, such as `convert-legacy-ldscores`, starts with its inputs instead.
2. **Make output the last step.**
   It starts with `--output-dir` and `--overwrite`, followed by options that change what is written, such as `--output-format` or `--continue-on-query-error`.
   Commands without `--output-dir` put `--overwrite` in a section titled `Output`.
3. **Within a section, put required inputs first, then commonly adjusted options.**
   Conditional requirements belong with the mode that needs them.
   Preserve the relationship between a primary option and its companion flags.
4. **Use default status as a tie-breaker.**
   Among similarly important options in the same section, put options needing an explicit value before those with usable defaults.
   A parser value of `None` can still represent automatic behavior; classify the effective user behavior, not the Python value.
5. **Keep alternatives and dependencies adjacent.**
   Do not separate paired case/control settings or the three LD-window alternatives merely to sort defaults.
6. **Show each flag once.**
   An argparse option belongs to one section, so an alternative mode that reuses a flag from a step names that flag and its step in the section description.

`infer-genome-build` has one positional exact path and help only; describe read-only reporting, successful unknown results, and the absence of output and logging options.
`plot` and `convert-h2-scale` have one or two inputs and keep plain headings (`Input result`, `Case prevalences`, `Output`, `Logging`) without step numbers.
Small commands need only the sections that help users; they do not need empty or artificial sections.

The current sequences are:

| Command | Numbered steps | Later sections |
| --- | --- | --- |
| `munge-sumstats` | 1 SNP identity and genome build; 2 raw summary statistics; 3 sample size; 4 output | SNP selection and quality filters; column overrides; computing resources; logging |
| `annotate` | 1 SNP identity and genome build; 2 annotation SNP grid; 3 query annotations; 4 output | logging |
| `ldscore` | 1 SNP identity and genome build; 2 baseline annotations; 3 query annotations; 4 reference panel, then PLINK-only settings; 5 LD window; 6 output | regression and reference SNPs; gene LD-score index alternative; computing resources; logging |
| `build-gene-ldscore-index` | 1 SNP identity and genome build; 2 baseline annotations; 3 gene coordinate catalog; 4 PLINK reference panel; 5 LD window; 6 output | regression and reference SNPs; computing resources; logging |
| `build-r2-matrix` | 1 SNP identity and source genome build; 2 PLINK reference panel; 3 LD window; 4 genetic maps and liftover; 5 output | R2 storage threshold; computing resources; logging |
| `h2`, `partitioned-h2`, `rg` | 1 munged summary statistics and LD scores; 2 output | liability scale; model settings; advanced estimation and SNP matching; computing resources (`partitioned-h2` only); logging |
| `quantile-h2` | 1 fitted model; 2 target annotation and quantiles; 3 original fitted-model inputs, then original BED and gene-list settings; 4 output | logging |
| `query-r2-matrix` | 1 R2 matrix; 2 SNP pairs; 3 output | logging |
| `convert-legacy-ldscores` | 1 legacy LD-score files; 2 output SNP identity and genome build; 3 output | logging |

The example commands in the guided tutorial and the option tables in the [command-line reference](../../reference/cli.md) follow the same order.
Section titles, the step order, and the structural rules above are checked by `tests/package/test_cli_help.py`.

## Usage lines

Use readable argument placeholders such as `FILE`, `DIR`, `COLUMN`, `BP`, or `N`.
Show the actual user-facing command in the usage line.
`annotate`, `ldscore`, `build-gene-ldscore-index` and `munge-sumstats` set a hand-written `usage` that shows the typical command in step order.
It brackets optional flags, groups alternatives as `(A | B)`, lists another mode on an `or:` line, and ends with `[options]`.
Every flag it names must exist and every required flag must appear; the help tests check both.
Other commands use the usage line that argparse generates from the order in which flags are added, so add flags in display order.
Equivalent entry points must preserve the same usage, groups, and descriptions; verify the rendered help rather than assuming that copying parser actions preserves them.
See [`cli.main` and `cli._command_parser`](../../../src/ldsc/cli.py).
Do not add an examples section to `--help`; worked commands belong in the user guide.

## Common mistakes and better descriptions

For path inputs, “accepts globs” alone is incomplete: identify `*`, quoting, match cardinality, and whether `@` is expanded. Do not claim that every suite reader enforces all 22 members; the full-coverage check is workflow-specific. Direct query `ldscore` and gene-index construction enforce it, whereas generic suite expansion can retain existing members. The [path specification](../../reference/input-paths.md#pattern-support-in-command-help) records these distinctions.

The poor descriptions below are existing text or representative paraphrases. Improved excerpts illustrate a specific rule; a complete flag description must also include any additional dependencies or omission behavior that applies to that command.

| Mistake | Improper description | Better wording or correction |
| --- | --- | --- |
| Internal implementation jargon | “Base pairs to add … before in-memory projection.” | “Extend each BED interval by this many base pairs at both ends, then find SNPs within the expanded interval. Default for BED inputs: 0, using the original intervals.” Add the explicit-value requirement when the same flag supports gene-list input. |
| Repeating the flag name | `--chunksize`: “Chunksize.” | “Read this many input rows per chunk. Default: 1,000,000.” |
| Unspecified logging behavior | `--log-level`: “Logging verbosity.” | Explain all four levels, `INFO` as the default, and the workflow log as the destination; see the full example below. |
| Missing boolean default and output location | “Write a reference-metadata sidecar next to the LD-score output.” | “Write reference-panel SNP metadata to `ref_metadata/chrN_meta.tsv.gz` under the output directory. Applies to PLINK input only. Default: off; these additional files are not written.” |
| Ambiguous window unit | “LD window size in SNPs.” | “LD window on either side of each SNP, measured by SNP count in the retained reference panel. No default; specify exactly one of `--ld-wind-snps`, `--ld-wind-kb`, or `--ld-wind-cm`.” This example describes direct `ldscore`. |
| Implementation default instead of omission behavior | `--N-min`: “Default: None.” | “Minimum per-variant sample size. If omitted or zero, use the 90th percentile of N divided by 1.5. Constant sample sizes are not filtered.” |
| Hidden required companion | `--N-cas-col`: “Case-count column.” | “Name of the per-variant case-count column. Requires `--N-con-col`; cannot be combined with `--N-col`.” Also explain inference and selection precedence in the complete description. |
| Fallback incorrectly called an override | `--N`: “Override sample size.” | “Constant sample-size fallback when per-variant N or paired case/control columns are absent. Input columns take precedence; this value takes precedence over `--N-cas`/`--N-con`.” |
| Unqualified mode-dependent default | `--padding-bp`: “Default: 0.” | “For BED inputs, omission means 0. Direct gene-list input requires an explicit value, including 0 to use gene bodies.” Include the flag's purpose and other applicable mode restrictions. |
| Unclear opt-in behavior | `--allow-unresolved-genes`: “Gene identifier policy.” | “Omit rejected gene identifiers and continue with the usable subset, recording omissions in diagnostics. Default: off (strict resolution; rejected identifiers stop the run). Requires gene-list input.” |
| Incorrect ordering | Logging, batch sizes, or advanced overrides before the main input file; sections grouped by category instead of by the steps a user follows. | Follow the step sequence above: numbered steps in tutorial order with output last, then optional and alternative sections, computing resources, and logging. |

For a workflow that writes a log, the logging description can read:

> Set the detail recorded in the workflow log. DEBUG includes troubleshooting details. INFO includes progress and summaries. WARNING includes warnings and errors. ERROR includes errors only. Each level also includes more severe messages. Default: INFO. Run headers and completion status are always recorded.

Do not promise that increasing `--log-level` prints progress to the terminal. Ordinary progress/debug messages go to the workflow log; console error reporting and explicit scientific notices have their own routing. See [workflow logging, “Console vs File Routing”](../../reference/logs-and-outputs.md#console-vs-file-routing).

The examples are grounded in [`_sumstats.munger.build_parser`](../../../src/ldsc/_sumstats/munger.py), [`_ldscore.calculator.build_parser` and `_write_one_ref_metadata_sidecar`](../../../src/ldsc/_ldscore/calculator.py), [`build_window_coordinates` and `get_block_lefts`](../../../src/ldsc/_ldscore/compute.py), and [`_WorkflowLoggingContext`](../../../src/ldsc/_runtime/workflow_logging.py). Preserve these behavioral distinctions when adapting the examples to other commands.

## Naming rules

Public CLI options and Python configuration fields follow these suffix rules:

| Suffix | Meaning | Examples |
| --- | --- | --- |
| `*_file` | One file-like input; an exact-one glob is allowed where the resolver supports it. | `raw_sumstats_file`, `sumstats_file`, `sumstats_snps_file`, `keep_indivs_file` |
| `*_sources` | One logical input that may resolve to several files, through globs, comma lists or `@` chromosome suites where the option supports them. | `baseline_annot_sources`, `query_annot_bed_sources`, `query_annot_gene_list_sources` |
| `*_dir` | A directory, as input or as output location. | `ldscore_dir`, `result_dir`, `h2_result_dir`, `output_dir` |

Established names outside these rules: `--plink-prefix`, `--pairs-file`, `--prevalence-manifest` and `--target-annot-column`.

Naming decisions to keep:

- Keep `--plink-prefix` rather than forcing it into the `*_file` or `*_sources` rule, because a PLINK input is the stem of a `.bed`, `.bim` and `.fam` trio.
- Use `--r2-matrix-dir` / `r2_matrix_dir` for a package-built R² matrix; it mirrors `--plink-prefix` as the single reference-panel input for its backend.
- Do not add `--output-name` or `--panel-name`: output file names are fixed, and users name a run by naming its output directory.
- Keep `--output-dir` required for commands that write results; `plot` and `convert-h2-scale` keep their fixed derived destinations without a CLI override.

## Review before accepting help changes

- Every visible option has a description with its purpose, relevant input details, effective default or omission behavior, and all dependencies.
- Conditional requirements, paired flags, mutual exclusion, and precedence agree with validation and workflow behavior.
- Default claims and output locations agree with the implementation; shared flag wording reflects any command-specific differences.
- Sections use the step headings and order above: numbered steps with output last, then optional and alternative sections, computing resources, and logging.
- Pattern cross-references and hand-written usage lines name only existing flags, and cross-references point to flags shown earlier in the help.
- Rendered `ldsc COMMAND --help` and `python -m ldsc COMMAND --help` remain readable, show the correct command, and preserve the intended grouping and ordering.
- Help-only changes preserve accepted flags, parser destinations, defaults, choices, dependencies, and scientific behavior. Relevant parser/CLI checks establish this; prose edits alone do not require numerical test suites.

The description contents and the step-ordered layout are approved; the layout replaced the earlier inputs-and-outputs-first ordering on 2026-10-28. All public command parsers use explicit help groups and descriptions, with shared wrapping in [`CLIHelpFormatter`](../../../src/ldsc/_runtime/cli_help.py) and consistent log-level wording in [`LOG_LEVEL_HELP`](../../../src/ldsc/_runtime/workflow_logging.py). No guidance decisions remain open.
