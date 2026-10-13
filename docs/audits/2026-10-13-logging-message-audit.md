# Logging message audit

Last updated on: 2026-10-13

## Scope and contract

Inspected the local `restructure` checkout at `911eacc`, including its pre-existing modifications.
The [current logging contract](../current/workflow-logging.md) and [troubleshooting guide](../troubleshooting.md) govern this pass.
An AST inventory covers 81 Python source files outside packaged data, including 876 directly constructed exceptions and logging/printing call sites.
The inventory is a discovery aid, not proof that every private guard can be reached through the CLI.
The review covers all 13 command surfaces and their shared input, artifact, configuration, workflow, and error-reporting paths.
Existing scientific validation, exception classes, result schemas, log levels, output ownership, and continuation policies remain the remediation boundary.

## Findings before remediation

### User repair advice contradicts current interfaces or diagnostics

Munging kernel messages recommend the unsupported `--debug` flag and tell users who lost every row during QC to relax thresholds and inspect a sidecar that does not contain those QC removals.
The INFO-list error exposes a private function and recommends ignoring the column without explaining the change to filtering.
Reference-panel and PLINK empty-result messages also prescribe relaxing filters without first identifying the failing check.
These messages need the observed failure, current flags, and a repair that distinguishes a mistaken setting from an intentional change to the analysis.

### Input errors are not consistently actionable or bounded

`InputGate.finish()` renders every issue even after writing a complete repair table.
New annotation, chromosome, index-component, and result-metadata checks often report only an invalid field or mismatch.
The primary error should identify the affected input/check and repair while preserving complete tables, issue counts, and exception causes.
Only summaries with a written complete table should omit detailed issue rows.
The chosen cap is 10 representative input issues, consistent with other bounded error summaries in the package.

### Saved artifacts need repair guidance at the consumer boundary

LD-score query manifests, gene-index components, plotting inputs, quantile fitted-model inputs, and h2 scale conversion include terse metadata/schema/path errors.
These consumers should identify the saved artifact and failed requirement and guide users to restore coherent files from the same run or regenerate the appropriate upstream artifact.
Missing files must not be confused with unreadable or malformed files.
Repairs must not tell users to invent scientific fields or silently combine artifacts from different runs.

### Default warnings need plain-language outcomes

Shared query-status warnings primarily report status/reason codes without explaining omissions or retained partial support.
The munging mean-chi-square warning embeds its own severity prefix and omits the measured value or next action.
The records should explain the outcome and scientific consequence while retaining existing levels, counts, reason codes, and diagnostic pointers.
Intentional direct stderr notices and legacy formatted scientific summaries are distinct from messages passed to the shared severity formatter.

## Healthy patterns to preserve

File-authoritative workflow logging, complete captured tracebacks, original input causes, distinct failed-input counts, and output-family preflight already provide useful diagnostic boundaries.
The pending gene-list CLI notice improvement distinguishes retained-reference support from variation on regression SNP rows and must be preserved.
Existing batch/chromosome progress is bounded and most older user errors already provide specific constraints and remedies.
Private numerical invariant guards retain their developer detail rather than being reclassified as user mistakes.

## Behavioral boundaries

This pass does not narrow the CLI's existing `ValueError` classification or redesign bootstrap exception routing.
Compatibility fallbacks and scientific warning/error decisions require their own behavior review; clearer wording alone does not authorize changing their calculations.
No production or HPC calculation is required for a message audit.
Synthetic inputs and existing regression/fixture tests provide verification.

## Remediation and validation

Implemented the message refinements within the existing runtime boundaries.
The site counts below compare against the checkout at the start of this audit, preserving the already pending documentation, gene-list notice, and test changes.
They count directly changed exception/log message expressions; assembled summaries and diagnostic formatting also changed.
The implementation changes 26 source files, including 185 direct exception-message expressions and 19 directly supplied logging-message expressions.

| Area | Implemented behavior |
| --- | --- |
| Shared staged input validation | Show up to 10 representative issues only after writing the complete repair table; preserve every table row and original exception cause, and retain all primary issues when no table was written. |
| Annotation and PLINK inputs | Identify the concrete file or trio, requested and observed chromosomes, and the matching-file repair; explain that renaming files does not change SNP coordinates. |
| `munge-sumstats` | Correct obsolete flags and sample-size wording; include all-rows-lost counts and current QC options, distinguish earlier QC removals from the written drop sidecar, and explain the consequence of ignoring an INFO field. |
| Gene-list queries | Explain skip/warning reason codes in words, distinguish missing SNP support from constant regression LD scores, and identify omitted or retained query outputs. |
| Saved results and LD scores | Identify missing files, malformed metadata, inconsistent query manifests, and unsupported result contracts; restore coherent metadata/data from one run or regenerate the appropriate upstream result. |
| Gene-index construction and loading | Add actionable configuration and component repairs, preserve exact validation and atomic publication, and explain empty intersections or chromosomes without regression rows. |
| R2 panels and queries | Explain missing chromosomes, unreadable/missing panel files, missing genetic distances, absent signed-r values, and the consequences of existing bias/scale/MAF fallbacks. |
| `quantile-h2`, plotting, and scale conversion | Give supported options and compatible saved-fit repairs; point quantile identity failures to their written row-level diagnostics and explain when changing a count policy requires regeneration and refitting. |
| Unexpected internal failures | Keep the original exception evidence and saved traceback routing, identify the interrupted operation, and ask for the command/package version rather than guessing a scientific cause. |

Seven focused tests cover bounded versus complete primary input reports, all captured causes, missing saved-file repair without writes, all-QC-row loss and its captured traceback, and plain-language skipped-query outcomes.
Existing INFO-column and coarse-file-block tests now check meaningful user context, consequences, and repairs rather than obsolete advice or required implementation jargon.
All pre-existing pending gene-list notice/test changes were preserved.

### Measured validation

- The affected-workflow run passed 388 tests and 10 subtests.
- The later munging, reader, legacy-conversion, and new-message run passed 192 tests and 6 subtests.
- The final quantile/annotation run passed 91 tests and 4 subtests.
- Every modified source file compiles.
  Static comparison against the audit's starting state finds unchanged scientific and validation logic after excluding the reviewed message expressions and presentation helpers.
- Full `pytest -q --disable-warnings` passed: 1,949 tests and 133 subtests, with 1 skip and 209 warnings (368.39 seconds).
- Sequential `python -m unittest discover -s tests -p 'test*.py' -v` completed 1,002 tests successfully, with 1 skip (70.156 seconds).
  Both runners used `/Users/wenbinwu/miniforge3/envs/ldsc3-dev/bin/python` from this repository.
  The skip is the installed-PyArrow dependency test; the warning summary contains rsID-mode build notices and the existing pandas dtype warning.
- All 35 CLI options referenced by changed messages resolve against the public parser across its 13 commands, including both genome-specific genetic-map options.
- All 60 static source troubleshooting links resolve to the current guide's 31 referenced headings.
- `git diff --check` passed, and the final branch/HEAD remain `restructure`/`911eacc`.

### Remaining behavioral review

Message cleanup preserves existing compatibility behavior.
In particular, the R2 reader still assumes scale 32767 when integer encoding lacks a scale, resolves sample-size-only bias metadata as raw, and allows a requested MAF filter to go unapplied when MAF metadata is unavailable.
Their warnings now state the consequence and repair.
Whether these cases should instead fail under the scientific-validation rule needs a separate behavior change with producer/override cases and numerical validation.
The MAF case is a concrete policy gap: a requested filter can currently be omitted during a successful run.

The CLI also retains its existing broad treatment of `ValueError` as a user error and its bootstrap routing.
Changing exception classification, warning visibility, or failure policy was outside this message-focused implementation.
Private numerical guards retain useful developer detail.
The inventory and passing tests do not establish that every possible corrupt input or abrupt worker termination was exercised.
No production data or HPC run was used for validation.
