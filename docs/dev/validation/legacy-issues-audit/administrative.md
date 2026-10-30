# Administrative and insufficient-information records

Last updated on: 2026-10-28

[Audit overview](README.md) · [Workflow counts and membership](workflow-index.md) · [Verification](verification.md)

## Why these records remain in coverage

Empty reports, external-service questions, withdrawn questions and insufficient-information cases are retained so that exhaustive coverage is inspectable.
A duplicate still receives its own assessment and points to the shared concern when supported.
No undocumented reproduction or fix is inferred from closure, silence or a missing attachment.

## Recommendations by kind

- **Code:** no change is justified by an unidentifiable symptom alone.
- **Documentation:** request the command, package revision, relevant log, input schema and a small non-sensitive reproducer in future issue guidance.
- **Features:** no new issue-management subsystem is needed.
- **Resources:** distinguish external-service access from a defect in the local package.

## Category summary

15 issues and 1 FAQ entry.

| Severity | Source | Concern | LDSC3 disposition | Evidence level |
| --- | --- | --- | --- | --- |
| Low | [#14](#issue-14) | Coverage of historical annotation/frequency parsers | Partial coverage | Source inspection |
| Low | [#473](#issue-473) | Unrelated MetaXcan model schema | Intentionally unsupported | Source inspection |
| Low | [#485](#issue-485) | External summary-statistics converter | Intentionally unsupported | Source inspection |
| Low | [#347](#issue-347) | Missing reproducer | Indeterminate | Source inspection |
| Informational | [#16](#issue-16) | Internal request to review the FAQ | Intentionally unsupported | Source inspection |
| Informational | [#209](#issue-209) | Empty sub-annotation issue | Indeterminate | Source inspection |
| Informational | [#261](#issue-261) | Empty LDSC issue | Indeterminate | Source inspection |
| Informational | [#266](#issue-266) | Issue posted to wrong repository | Indeterminate | Source inspection |
| Informational | [#301](#issue-301) | Unspecified error | Indeterminate | Source inspection |
| Informational | [#324](#issue-324) | Withdrawn genetic-correlation content | Indeterminate | Source inspection |
| Informational | [#343](#issue-343) | Closed placeholder | Indeterminate | Source inspection |
| Informational | [#344](#issue-344) | Closed placeholder | Indeterminate | Source inspection |
| Informational | [#346](#issue-346) | Wrong repository | Indeterminate | Source inspection |
| Informational | [#365](#issue-365) | Unrelated pandas indexing | Indeterminate | Source inspection |
| Informational | [#455](#issue-455) | Empty allele issue | Indeterminate | Source inspection |
| Informational | [FAQ-01](#faq-01) | Forking and attribution | Apparently prevented | Source inspection |

## Individual assessments

<a id="issue-14"></a>

### #14: write tests for everything to do with annot/frq in parse.py

[Original issue](https://github.com/bulik/ldsc/issues/14) · Source state: **closed** · Severity: **low** · Accessible comments reviewed: **0**.

Workflow tags: Cross-cutting / not assigned to a requested group.

**Trigger and harm.** A developer requests tests and suggests shortening an unfamiliar parser.
Missing tests can conceal alignment and frequency/count regressions, but file length is not a defect.

**Legacy disposition.** Closed and linked to 288ae0a73432b533313af324b4f0d7ef1acae51d; referenced status alone does not establish complete test coverage.

**LDSC3 assessment — Partial coverage.** LDSC3 no longer has the old parse.annot function.
It has domain-specific annotation parsing, numeric validation, storage and overlap tests, plus legacy converter tests.
These provide targeted coverage, but no test inventory alone proves all historical edge cases are protected; issue-derived cases should be mapped to concrete tests or flagged as gaps.

**Recommended prevention.** Use an issue-to-test coverage column; do not refactor parsers solely because the legacy request mentioned their length.

**Evidence and verification.** source-inspected
[tests/annotation/test_annotation_format.py](<../../../../tests/annotation/test_annotation_format.py>); [tests/annotation/test_annotation_semantics.py](<../../../../tests/annotation/test_annotation_semantics.py>); [tests/ldscore/test_legacy_ldscore_converter.py](<../../../../tests/ldscore/test_legacy_ldscore_converter.py>).

<a id="issue-473"></a>

### #473: Unexpected error: Execution failed on sql 'SELECT gene, phi FROM extra': no such column: phi

[Original issue](https://github.com/bulik/ldsc/issues/473) · Source state: **closed** · Severity: **low** · Accessible comments reviewed: **0**.

Workflow tags: Cross-cutting / not assigned to a requested group.

**Trigger and harm.** SPrediXcan queries a phi column absent from a model database. A TWAS association tool fails; no LDSC invocation appears in the report.

**Legacy disposition.** Closed without an accessible explanation or LDSC code fix.

**LDSC3 assessment — Intentionally unsupported.** LDSC3 does not own MetaXcan SQLite model schemas or its calibration implementation.
Supplying a GWAS heritability number to another tool does not make that tool failure a regression-package defect.

**Recommended prevention.** Retain coverage, route diagnosis to the model/tool provider and avoid adding database adapters to LDSC3 for this unrelated issue.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[Source](<https://github.com/bulik/ldsc/issues/473>); [src/ldsc/cli.py](<../../../../src/ldsc/cli.py>).

<a id="issue-485"></a>

### #485: python sumstats.py csv not work

[Original issue](https://github.com/bulik/ldsc/issues/485) · Source state: **closed** · Severity: **low** · Accessible comments reviewed: **0**.

Workflow tags: `munge-sumstats`.

**Trigger and harm.** Runs python sumstats.py csv/qc from an unidentified tool and gets a tab-delimited output with a csv suffix.
Filename extensions are confused with actual delimiters; the responsible implementation is unspecified.

**Legacy disposition.** Closed without an accessible resolution or evidence that bulik/ldsc owns those subcommands.

**LDSC3 assessment — Intentionally unsupported.** LDSC3 has munge-sumstats, not the posted csv/qc interface.
Its file contracts specify content independently of an arbitrary output basename; changing its parser would not repair another tool.

**Recommended prevention.** Retain the record and request the actual converter repository/version; document accepted content/separator and avoid promising extension-driven conversion.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[src/ldsc/cli.py](<../../../../src/ldsc/cli.py>); [docs/reference/formats/column-schema.md](<../../../reference/formats/column-schema.md>).

<a id="issue-347"></a>

### #347: error in munge_sumstats.py, ValueError: invalid literal for float(): 05:14:57

[Original issue](https://github.com/bulik/ldsc/issues/347) · Source state: **closed** · Severity: **low** · Accessible comments reviewed: **0**.

Workflow tags: `munge-sumstats`.

**Trigger and harm.** The title reports converting a clock-like token to float; the body is empty. A malformed numeric field may prevent munging, but its source and column are unavailable.

**Legacy disposition.** Closed with no accessible explanation or verified code fix.

**LDSC3 assessment — Indeterminate.** LDSC3 coerces numeric fields and records invalid-value counts, but no evidence ties that behavior to this absent input.
The title alone cannot establish a runtime repair.

**Recommended prevention.** Request or document the offending row and mapped column when reporting numeric conversion failures; do not classify the historical case as fixed.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[src/ldsc/_sumstats/munge_kernel.py:148](<../../../../src/ldsc/_sumstats/munge_kernel.py#L148>); [src/ldsc/_sumstats/munge_kernel.py:271](<../../../../src/ldsc/_sumstats/munge_kernel.py#L271>).

<a id="issue-16"></a>

### #16: hilary should check out the FAQ

[Original issue](https://github.com/bulik/ldsc/issues/16) · Source state: **closed** · Severity: **informational** · Accessible comments reviewed: **0**.

Workflow tags: Cross-cutting / not assigned to a requested group.

**Trigger and harm.** A maintainer asks another maintainer to read docs/FAQ.
This is a coordination task, not a reproducible software defect.

**Legacy disposition.** Closed; no substantive technical claim or fix is recorded in the issue itself.

**LDSC3 assessment — Intentionally unsupported.** There is no current LDSC3 behavior to reproduce from this issue.
The externally linked FAQ is being audited independently, so it should supply actual conceptual questions and evidence instead of attributing invented requirements to this empty coordination record.
Count this issue in overall coverage but not in a workflow-frequency group.

**Recommended prevention.** Retain an administrative disposition and link the separate FAQ assessment.

**Evidence and verification.** source-inspected
[Source](<https://github.com/bulik/ldsc/issues/16>); [docs/user-guide/troubleshooting.md](<../../../user-guide/troubleshooting.md>).

<a id="issue-209"></a>

### #209: doubt with running sub-annotations

[Original issue](https://github.com/bulik/ldsc/issues/209) · Source state: **closed** · Severity: **informational** · Accessible comments reviewed: **0**.

Workflow tags: Cross-cutting / not assigned to a requested group.

**Trigger and harm.** Title refers to sub-annotations, but body and comments are empty.
No reproducible technical or scientific claim can be assessed.

**Legacy disposition.** Closed with no resolution details or linked fix.

**LDSC3 assessment — Indeterminate.** The title alone cannot identify whether the user intended baseline partitioning, a query annotation or a restricted SNP analysis.
Assigning a concrete persistent bug or use-case frequency from that title would overstate the available evidence.

**Recommended prevention.** Retain the issue in the coverage ledger as insufficient information; do not infer a defect or a resolved status.

**Evidence and verification.** source-inspected
[Source](<https://github.com/bulik/ldsc/issues/209>).

<a id="issue-261"></a>

### #261: LDSC

[Original issue](https://github.com/bulik/ldsc/issues/261) · Source state: **closed** · Severity: **informational** · Accessible comments reviewed: **0**.

Workflow tags: Cross-cutting / not assigned to a requested group.

**Trigger and harm.** An empty issue titled LDSC has no body or discussion.
No user scenario or defect is identifiable.

**Legacy disposition.** Closed without explanation or linked fix.

**LDSC3 assessment — Indeterminate.** There is no evidence to compare with a specific LDSC3 behavior or to assign one of the seven scientific workflows.
Keeping it in the inventory demonstrates exhaustive coverage without inventing a bug, method limitation or resolution.

**Recommended prevention.** Mark as no assessable technical content and exclude from workflow-specific counts unless new evidence becomes available.

**Evidence and verification.** source-inspected
[Source](<https://github.com/bulik/ldsc/issues/261>).

<a id="issue-266"></a>

### #266: .

[Original issue](https://github.com/bulik/ldsc/issues/266) · Source state: **closed** · Severity: **informational** · Accessible comments reviewed: **1**.

Workflow tags: Cross-cutting / not assigned to a requested group.

**Trigger and harm.** Author says the issue was posted to the wrong GitHub page and retitles it to a dot.
No LDSC behavior is implicated.

**Legacy disposition.** Author acknowledges misposting and closes; no code fix.

**LDSC3 assessment — Indeterminate.** There is no scientific or package-specific assertion to test against LDSC3.
Treat this as administrative material, keep one inventory row and do not infer use-case tags from discarded content or closure state.

**Recommended prevention.** Count it in total issue coverage and administrative totals only.

**Evidence and verification.** source-inspected
[Source](<https://github.com/bulik/ldsc/issues/266#issuecomment-770930424>).

<a id="issue-301"></a>

### #301: nothing

[Original issue](https://github.com/bulik/ldsc/issues/301) · Source state: **open** · Severity: **informational** · Accessible comments reviewed: **0**.

Workflow tags: Cross-cutting / not assigned to a requested group.

**Trigger and harm.** Issue is titled nothing and contains only error.
No trigger, impact or reproducible behavior can be inferred.

**Legacy disposition.** Open with only a rename event and no resolution.

**LDSC3 assessment — Indeterminate.** The record contains no information to map to a command, algorithm or reference resource.
It belongs in exhaustive issue coverage but should not inflate any of the seven workflow counts or be marked fixed because unrelated functionality exists.

**Recommended prevention.** Mark insufficient information and retain the source link without inventing a diagnosis.

**Evidence and verification.** source-inspected
[Source](<https://github.com/bulik/ldsc/issues/301>).

<a id="issue-324"></a>

### #324: Problems when perform LD score regression by --rg

[Original issue](https://github.com/bulik/ldsc/issues/324) · Source state: **open** · Severity: **informational** · Accessible comments reviewed: **0**.

Workflow tags: Cross-cutting / not assigned to a requested group.

**Trigger and harm.** Title mentions genetic correlation but the entire body is Done.
No remaining technical claim or resolution can be audited.

**Legacy disposition.** Open with no comments or linked fix.

**LDSC3 assessment — Indeterminate.** There is no evidence that a genetic-correlation bug persisted or was fixed, and the title does not reveal a concrete statistical issue.
Preserve the row in the exhaustive inventory, but avoid counting an inferred user problem from withdrawn content.

**Recommended prevention.** Mark no assessable technical content; do not assign a workflow frequency solely from the title.

**Evidence and verification.** source-inspected
[Source](<https://github.com/bulik/ldsc/issues/324>).

<a id="issue-343"></a>

### #343: close

[Original issue](https://github.com/bulik/ldsc/issues/343) · Source state: **closed** · Severity: **informational** · Accessible comments reviewed: **0**.

Workflow tags: Cross-cutting / not assigned to a requested group.

**Trigger and harm.** Thread contains only closure text or a wrong-repository notice.
No reproducible LDSC behavior or scientific question is supplied.

**Legacy disposition.** No technical resolution or code fix can be inferred from this record.

**LDSC3 assessment — Indeterminate.** There is no actionable condition to map to LDSC3.
The issue is retained in the complete inventory without inventing a defect or workflow tag.

**Recommended prevention.** No package change; preserve this explicit nontechnical disposition.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[Source](<https://github.com/bulik/ldsc/issues/343>).

<a id="issue-344"></a>

### #344: CLosed

[Original issue](https://github.com/bulik/ldsc/issues/344) · Source state: **closed** · Severity: **informational** · Accessible comments reviewed: **0**.

Workflow tags: Cross-cutting / not assigned to a requested group.

**Trigger and harm.** Thread contains only closure text or a wrong-repository notice.
No reproducible LDSC behavior or scientific question is supplied.

**Legacy disposition.** No technical resolution or code fix can be inferred from this record.

**LDSC3 assessment — Indeterminate.** There is no actionable condition to map to LDSC3.
The issue is retained in the complete inventory without inventing a defect or workflow tag.

**Recommended prevention.** No package change; preserve this explicit nontechnical disposition.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[Source](<https://github.com/bulik/ldsc/issues/344>).

<a id="issue-346"></a>

### #346: Sorry wrong github

[Original issue](https://github.com/bulik/ldsc/issues/346) · Source state: **open** · Severity: **informational** · Accessible comments reviewed: **0**.

Workflow tags: Cross-cutting / not assigned to a requested group.

**Trigger and harm.** Thread contains only closure text or a wrong-repository notice.
No reproducible LDSC behavior or scientific question is supplied.

**Legacy disposition.** No technical resolution or code fix can be inferred from this record.

**LDSC3 assessment — Indeterminate.** There is no actionable condition to map to LDSC3.
The issue is retained in the complete inventory without inventing a defect or workflow tag.

**Recommended prevention.** No package change; preserve this explicit nontechnical disposition.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[Source](<https://github.com/bulik/ldsc/issues/346>).

<a id="issue-365"></a>

### #365: Unalignable boolean Series provided as indexer (index of the boolean Series and of the indexed object do not match).

[Original issue](https://github.com/bulik/ldsc/issues/365) · Source state: **open** · Severity: **informational** · Accessible comments reviewed: **0**.

Workflow tags: Cross-cutting / not assigned to a requested group.

**Trigger and harm.** The entire reproducer is generic pandas CSV filtering, without an LDSC call.
No LDSC behavior is implicated.

**Legacy disposition.** Open; no package resolution is supplied.

**LDSC3 assessment — Indeterminate.** A Boolean Series over DataFrame columns is being used as a row selector in the posted snippet.
This independent pandas example does not establish a defect in LDSC3 or justify changing its indexing.

**Recommended prevention.** Keep the record in coverage and mark it outside the package audit; direct future reports to a minimal LDSC workflow reproducer.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[Source](<https://github.com/bulik/ldsc/issues/365>).

<a id="issue-455"></a>

### #455: snpsvalid alleles

[Original issue](https://github.com/bulik/ldsc/issues/455) · Source state: **closed** · Severity: **informational** · Accessible comments reviewed: **0**.

Workflow tags: Cross-cutting / not assigned to a requested group.

**Trigger and harm.** Closed issue has a title about valid alleles but no body or comments.
No failure mechanism or scientific question can be established.

**Legacy disposition.** Closed without an accessible explanation or code fix.

**LDSC3 assessment — Indeterminate.** The inventory preserves the record, but it cannot be mapped to an allele-validation defect merely from its title.

**Recommended prevention.** No package change; require a minimal input and actual diagnostic for a future actionable report.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[Source](<https://github.com/bulik/ldsc/issues/455>).

<a id="faq-01"></a>

### FAQ-01: Forking and attribution

[FAQ source](https://github.com/bulik/ldsc/wiki/FAQ) · Source-order identifier: **FAQ-01** · Severity: **informational**.

Workflow tags: Cross-cutting / not assigned to a requested group.

**Legacy disposition.** Historical FAQ guidance; assessed independently rather than treated as a code-fix record.

**LDSC3 assessment — Apparently prevented.** LDSC3 retains GPLv3 licensing, NOTICE and software/method citation instructions.
Forking permission does not settle the separate provenance or reuse terms of bundled scientific resources.

**Recommended prevention.** Retain notices and cite the methods actually used; resolve resource attribution separately.

**Evidence and verification.** Source-inspected; relevant executed checks are identified in verification.md.
[LICENSE](<../../../../LICENSE>); [NOTICE](<../../../../NOTICE>); [CITATION.cff](<../../../../CITATION.cff>); [README.md:44](<../../../../README.md#documentation>).

