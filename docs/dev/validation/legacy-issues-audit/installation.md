# Installation and dependencies

Last updated on: 2026-10-28

[Audit overview](README.md) · [Workflow counts and membership](workflow-index.md) · [Verification](verification.md)

## Installation choices worth documenting

The current README installs an editable checkout with development extras through conda.
That is appropriate for contributors, but routine users do not necessarily need editable source, pytest or every optional workflow dependency.
The distribution is called `ldsc`; that name alone does not identify this repository's LDSC3 code.
The currently visible [PyPI project](https://pypi.org/project/ldsc/) describes another distribution/release history, rather than this audited `restructure` revision.
Do not recommend an unqualified `pip install ldsc` as installation of this audited code.

| Choice | Appropriate use | Proposed recipe | Verification in this audit |
| --- | --- | --- | --- |
| Existing conda + editable `.[dev]` | Contributors changing source and running tests | Retain the documented development procedure. | Existing environment loaded this checkout; 499 focused tests and 36 subtests passed. Environment creation was not repeated. |
| Standard-library `venv` + non-editable local installation | Routine use from a chosen checkout | Use a supported Python, then `python -m pip install ".[plink]"` for PLINK workflows, or the base package for workflows not requiring extras. | Reviewed against package metadata; no clean-environment installation trial. |
| Repository-specific VCS installation | Routine use without retaining a working checkout | Pin this repository and an exact published commit; select only required extras. | Recipe structure and extras reviewed; remote installation and authentication not exercised. |
| Conda + non-editable local installation | Users who prefer conda's environment management | Create a supported Python environment and install the chosen checkout normally instead of editable mode. | Metadata-reviewed alternative; not separately installed. |

A minimal local-checkout recipe, proposed rather than clean-environment tested, is:

```bash
python3.13 -m venv .venv
. .venv/bin/activate
python -m pip install ".[plink]"
ldsc --help
python -m ldsc --help
```

Run it inside the selected LDSC3 checkout.
Python 3.11–3.13 is the declared support range; use an available supported interpreter rather than assume every machine has `python3.13`.
A VCS equivalent for this specific audited commit would be:

```bash
python -m pip install "ldsc[plink] @ git+https://github.com/WenbinWu2001/ldsc_py3.git@9d095e66a8fddff02831c644d286ed274876503b"
```

This requires access to that repository and commit; the audit did not test the remote recipe or establish a public release.
A moving branch such as `restructure` is convenient for development but less reproducible than an exact available commit.
No installation command above was executed to change the current environment.

The base requirements and extras are declared in [`setup.py`](../../../../setup.py).
`plink` adds `bitarray`; `liftover` adds `pyliftover`; `all` includes both workflow extras; `test` adds pytest and `dev` adds all workflow extras plus tests.
BED and gene-list interval projection uses NumPy and does not require installing pybedtools or a separate BED executable.
The current declared NumPy/pandas/SciPy ranges should not be replaced with Python 2-era pins suggested in old issues.
Using `python -m pip` ties installation to the selected interpreter and avoids the common `pip`/`pip3` mismatch.

## Recommendations by kind

- **Code:** address a packaging defect only when an actual supported-environment failure is reproduced; old ABI failures alone do not establish a current defect.
- **Documentation:** give separate routine-use and contributor procedures, state the selected branch/revision and extras, and explain how to confirm the imported module and command location.
- **Features:** no new installer framework, container platform or environment manager is needed to resolve these recurring problems.
- **Resources:** keep reference-data acquisition separate from package installation; a successful import does not supply scientific inputs.

## Category summary

60 issues and 0 FAQ entries.

| Severity | Source | Concern | LDSC3 disposition | Evidence level |
| --- | --- | --- | --- | --- |
| High | [#178](#issue-178) | Python2/3 and effect-column confusion | Partial coverage | Source inspection |
| High | [#373](#issue-373) | Distribution identity mismatch | Partial coverage | Source inspection |
| High | [#464](#issue-464) | External R wrapper and sparse overlap | Intentionally unsupported | Source inspection |
| High | [#91](#issue-91) | Version-dependent NaN/Inf covariance division | Indeterminate | Source inspection |
| High | [#173](#issue-173) | NumPy/pandas binary incompatibility | Apparently prevented | Source inspection |
| High | [#180](#issue-180) | Legacy NumPy pin | Apparently prevented | Source inspection |
| High | [#205](#issue-205) | macOS legacy NumPy ABI | Apparently prevented | Source inspection |
| High | [#277](#issue-277) | Pandas source build fails in conda | Apparently prevented | Source inspection |
| High | [#287](#issue-287) | macOS pandas0.20 installation failure | Apparently prevented | Source inspection |
| Moderate | [#26](#issue-26) | Installation failures and long mixed-topic support thread | Partial coverage | Source inspection |
| Moderate | [#36](#issue-36) | Tests fail across pandas versions and platforms | Partial coverage | Source inspection |
| Moderate | [#42](#issue-42) | SciPy LAPACK binary mismatch | Partial coverage | Source inspection |
| Moderate | [#55](#issue-55) | Shell commands entered into Python and mixed installations | Partial coverage | Source inspection |
| Moderate | [#86](#issue-86) | Release tags lag maintained master | Partial coverage | Source inspection |
| Moderate | [#112](#issue-112) | Munging appears hung during initial input reading | Partial coverage | Source inspection |
| Moderate | [#177](#issue-177) | Installable releases and version tags | Partial coverage | Source inspection |
| Moderate | [#358](#issue-358) | Legacy dependency slowdown | Partial coverage | Source inspection |
| Moderate | [#372](#issue-372) | Python lacks bzip2 support | Partial coverage | Source inspection |
| Moderate | [#445](#issue-445) | Missing pip in old environment | Partial coverage | Source inspection |
| Moderate | [#443](#issue-443) | MTAG count-width mismatch | Intentionally unsupported | Source inspection |
| Moderate | [#206](#issue-206) | Incomplete Git clone | Indeterminate | Source inspection |
| Moderate | [#423](#issue-423) | Native-library illegal instruction | Indeterminate | Source inspection |
| Moderate | [#446](#issue-446) | Fork-specific usecols error | Indeterminate | Source inspection |
| Moderate | [#492](#issue-492) | Native genotype-computation crash | Indeterminate | Source inspection |
| Moderate | [#30](#issue-30) | Old pandas rejects drop_duplicates subset keyword | Apparently prevented | Source inspection |
| Moderate | [#39](#issue-39) | Wrong Python interpreter causes syntax error | Apparently prevented | Source inspection |
| Moderate | [#64](#issue-64) | Legacy statistical tests have too few replicates | Apparently prevented | Source inspection |
| Moderate | [#90](#issue-90) | Legacy Python 2 code executed with Python 3 | Apparently prevented | Source inspection |
| Moderate | [#93](#issue-93) | Pandas boolean-mask compatibility in allele merge | Apparently prevented | Source inspection |
| Moderate | [#104](#issue-104) | Pandas upgrade breaks old masking path | Apparently prevented | Source inspection |
| Moderate | [#121](#issue-121) | Pandas mask regression with a huge keep list | Apparently prevented | Source inspection |
| Moderate | [#124](#issue-124) | Removed pandas DataFrame.sort API | Apparently prevented | Source inspection |
| Moderate | [#132](#issue-132) | Repeated pandas negative-integer mask error | Apparently prevented | Source inspection |
| Moderate | [#148](#issue-148) | Cell-type tutorial hits old pandas mask bug | Apparently prevented | Source inspection |
| Moderate | [#157](#issue-157) | Subset keep list triggers old pandas mask failure | Apparently prevented | Source inspection |
| Moderate | [#159](#issue-159) | Partial Python 3 patch exposes more legacy failures | Apparently prevented | Source inspection |
| Moderate | [#161](#issue-161) | Python 3 legacy allele_merge failure repeated | Apparently prevented | Source inspection |
| Moderate | [#202](#issue-202) | Mixed tabs and spaces | Apparently prevented | Source inspection |
| Moderate | [#230](#issue-230) | Legacy TabError before help | Apparently prevented | Source inspection |
| Moderate | [#276](#issue-276) | Script command not found | Apparently prevented | Source inspection |
| Moderate | [#296](#issue-296) | Python3 bytes/header and Boolean-mask porting | Apparently prevented | Source inspection |
| Moderate | [#300](#issue-300) | Private IRWLS import fails | Apparently prevented | Source inspection |
| Moderate | [#350](#issue-350) | Windows BED dependency | Apparently prevented | Source inspection |
| Moderate | [#351](#issue-351) | Wrong Python interpreter | Apparently prevented | Source inspection |
| Moderate | [#354](#issue-354) | Incomplete Python-3 port | Apparently prevented | Source inspection |
| Moderate | [#355](#issue-355) | Hand-edited port state | Apparently prevented | Source inspection |
| Moderate | [#361](#issue-361) | pybedtools scalar conversion | Apparently prevented | Source inspection |
| Moderate | [#369](#issue-369) | Import collision and unsafe repairs | Apparently prevented | Source inspection |
| Moderate | [#379](#issue-379) | Legacy pandas version gate | Apparently prevented | Source inspection |
| Moderate | [#380](#issue-380) | Legacy pandas version gate | Apparently prevented | Source inspection |
| Moderate | [#400](#issue-400) | Binary log handle in port | Apparently prevented | Source inspection |
| Moderate | [#417](#issue-417) | Terminal selects wrong Python | Apparently prevented | Source inspection |
| Moderate | [#432](#issue-432) | Obsolete ARM environment | Apparently prevented | Source inspection |
| Moderate | [#440](#issue-440) | Annotation dependency hang and BIM type | Apparently prevented | Source inspection |
| Moderate | [#442](#issue-442) | Python version contradicts assumption | Apparently prevented | Source inspection |
| Moderate | [#453](#issue-453) | Vendored package init filename | Apparently prevented | Source inspection |
| Moderate | [#490](#issue-490) | Modified allele-merge code | Apparently prevented | Source inspection |
| Moderate | [#504](#issue-504) | pip executable from wrong interpreter | Apparently prevented | Source inspection |
| Low | [#397](#issue-397) | MTAG vendored Python-2 copy | Intentionally unsupported | Source inspection |
| Low | [#398](#issue-398) | MTAG vendored Python-2 copy | Intentionally unsupported | Source inspection |

## Individual assessments

<a id="issue-178"></a>

### #178: Anaconda3 does not work, really...

[Original issue](https://github.com/bulik/ldsc/issues/178) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **7**.

Workflow tags: `munge-sumstats`.

**Trigger and harm.** Anaconda3 is mistaken for a Python3 environment; discussion also questions BETA,0 and duplicate frequency columns.
Environment failures and wrong effect direction can prevent use or mislead genetic-correlation analysis.

**Legacy disposition.** Comments distinguish conda distribution from interpreter, report NumPy workaround and external Python3 forks; no merge making legacy LDSC Python3-compatible is identified.

**LDSC3 assessment — Partial coverage.** LDSC3 runs on Python3.11-3.13 and documents that BETA,0 supplies direction around a zero null while P determines magnitude.
It still requires a single unambiguous frequency mapping and retains chunked input with whole-table accumulation.

**Recommended prevention.** Explain the active interpreter, sign/null semantics and explicit column selection; link large-file memory guidance.

**Evidence and verification.** source-inspected
[src/ldsc/_sumstats/munge_input.py:prepare_munge_input](<../../../../src/ldsc/_sumstats/munge_input.py>); [src/ldsc/_sumstats/munge_kernel.py:munge_sumstats](<../../../../src/ldsc/_sumstats/munge_kernel.py>).

Legacy implementation context: [`environment.yml`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/environment.yml), [`requirements.txt`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/requirements.txt).

<a id="issue-373"></a>

### #373: Presence of ldsc 2.0.1 on PyPI but not on GitHub

[Original issue](https://github.com/bulik/ldsc/issues/373) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **2**.

Workflow tags: Cross-cutting / not assigned to a requested group.

**Trigger and harm.** The PyPI package advertises 2.0.1 while bulik/ldsc GitHub advertises 1.0.1.
Users may unknowingly mix independently maintained implementations and interpret version-number changes as a single release lineage.

**Legacy disposition.** The linked Python-3 port PR #360 remains open and unmerged by API verification.
The PyPI implementation and later cbiit project are separate distributions, not a verified master release fix.

**LDSC3 assessment — Partial coverage.** LDSC3 also uses distribution name ldsc, so an unqualified pip install ldsc can select another implementation.
The current source-install README is safer but the package identity needs explicit explanation before recommending PyPI.

**Recommended prevention.** Recommend a repository/branch-pinned installation and display installed source/version provenance; do not promise that PyPI hosts this branch.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[setup.py](<../../../../setup.py>); [README.md](<../../../../README.md>); [Source](<https://github.com/bulik/ldsc/pull/360>).

Legacy implementation context: [`environment.yml`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/environment.yml), [`requirements.txt`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/requirements.txt).

<a id="issue-464"></a>

### #464: Heritability Estimation Error!!

[Original issue](https://github.com/bulik/ldsc/issues/464) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **0**.

Workflow tags: `h2`.

**Trigger and harm.** An R ldsc_h2 wrapper retains only 303 of 4685 SNPs before a matrix-operation error. A numerically repaired wrapper would still leave a very sparse, potentially selected regression dataset.

**Legacy disposition.** Open; no diagnosed source or fix in the accessible issue.

**LDSC3 assessment — Intentionally unsupported.** The traceback is from a separate R implementation, not LDSC3 Python.
LDSC3 input/fit validation can diagnose low row counts, but it cannot validate the external wrapper or infer whether this is genome-wide rather than selected GWAS data.

**Recommended prevention.** Route the implementation error to the R wrapper and independently check full-genome coverage, matching identifiers and minimum estimable block/model sizes before interpreting any estimate.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[src/ldsc/_regression/runner.py](<../../../../src/ldsc/_regression/runner.py>); [docs/user-guide/workflows/heritability.md](<../../../user-guide/workflows/heritability.md>); [Source](<https://github.com/bulik/ldsc/issues/464>).

Legacy implementation context: [`environment.yml`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/environment.yml), [`requirements.txt`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/requirements.txt).

<a id="issue-91"></a>

### #91: Partitioned heritability

[Original issue](https://github.com/bulik/ldsc/issues/91) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **5**.

Workflow tags: `baseline-partitioned-h2`.

**Trigger and harm.** Partitioned fits fail in coefficient covariance division under a newer NumPy/pandas stack. A late numerical exception is hard to distinguish from singular inputs or a dependency regression.

**Legacy disposition.** User confirms an older dependency environment works; maintainers did not establish the underlying defect or a precise source fix.

**LDSC3 assessment — Indeterminate.** LDSC3 has current dependency bounds, finite/schema gates and richer jackknife diagnostics, but the original numerical arrays are unavailable.
The matching division still exists conceptually; a working modern suite does not prove the exact historical NaN/Inf path is impossible.
Root V3 separately demonstrates a zero-total derived-ratio boundary, which is not proof of this dataset cause.

**Recommended prevention.** Validate finite inputs and intermediate variance outputs with informative model context; preserve the indeterminate historical cause instead of prescribing blanket downgrades.

**Evidence and verification.** source-inspected
[src/ldsc/_regression/estimators.py:LD_Score_Regression._coef](<../../../../src/ldsc/_regression/estimators.py>); [src/ldsc/_regression/jackknife.py](<../../../../src/ldsc/_regression/jackknife.py>); [docs/dev/validation/legacy-issues-audit/verification.md#v3-zero-total-derived-proportion](<verification.md#v3-zero-total-derived-proportion>).

Legacy implementation context: [`environment.yml`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/environment.yml), [`requirements.txt`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/requirements.txt).

<a id="issue-173"></a>

### #173: pandas error

[Original issue](https://github.com/bulik/ldsc/issues/173) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **2**.

Workflow tags: Cross-cutting / not assigned to a requested group.

**Trigger and harm.** Importing legacy LDSC fails with numpy.ufunc binary-size mismatch before help can run. A dependency ABI conflict makes every command unavailable.

**Legacy disposition.** Two users report numpy1.16.0 solves their environment; no merged fix is established in this thread.

**LDSC3 assessment — Apparently prevented.** LDSC3 declares supported modern Python and dependency bounds and installs through package metadata, removing the specific Python2/NumPy1.12 combination.
Arbitrarily mixed conda and pip binaries can still fail; this is environment compatibility rather than a numerical estimator defect.

**Recommended prevention.** Recommend one isolated supported environment and python -m pip; verify import and CLI after installation.

**Evidence and verification.** source-inspected
[setup.py:python_requires and install_requires](<../../../../setup.py>); [docs/user-guide/installation.md](<../../../user-guide/installation.md>).

Legacy implementation context: [`environment.yml`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/environment.yml), [`requirements.txt`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/requirements.txt).

<a id="issue-180"></a>

### #180: numpy==1.16.0 !!!

[Original issue](https://github.com/bulik/ldsc/issues/180) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **2**.

Workflow tags: Cross-cutting / not assigned to a requested group.

**Trigger and harm.** A NumPy1.12 environment produces pandas ABI failure; upgrading to1.16.0 works for several users.
The documented installation path fails across otherwise ordinary Linux environments.

**Legacy disposition.** PR197 proposed a NumPy update but closed unmerged; its referenced commit is not proof that master received the fix. [Source](<https://github.com/bulik/ldsc/pull/197>)

**LDSC3 assessment — Apparently prevented.** Modern LDSC3 dependency bounds exclude the reported legacy combination and support an installable Python3 package.
This prevents that exact old pin conflict, while platform-specific wheel availability and mixed binary channels still deserve installation smoke checks.

**Recommended prevention.** Avoid recommending historical pins; document supported environments and smoke-test help/import on supported Python versions.

**Evidence and verification.** source-inspected
[setup.py:python_requires and install_requires](<../../../../setup.py>); [docs/user-guide/installation.md](<../../../user-guide/installation.md>).

Legacy implementation context: [`environment.yml`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/environment.yml), [`requirements.txt`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/requirements.txt).

<a id="issue-205"></a>

### #205: numpy.ufunc has the wrong size (MacOS 10.12.6)

[Original issue](https://github.com/bulik/ldsc/issues/205) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **1**.

Workflow tags: Cross-cutting / not assigned to a requested group.

**Trigger and harm.** Legacy LDSC import on macOS10.12.6 fails with a NumPy ufunc-size mismatch.
All analysis commands are blocked before input processing.

**Legacy disposition.** Reporter and commenter confirm changing to NumPy1.16; no merged master fix is established by this issue.

**LDSC3 assessment — Apparently prevented.** LDSC3 avoids the obsolete Python2/NumPy1.12 dependency set and declares supported Python3 versions.
This does not imply current dependencies support macOS10.12 or that mixed binary installs cannot fail; platform availability must be tested rather than inferred from Python syntax compatibility.

**Recommended prevention.** State supported runtime/platform assumptions and prefer a clean environment with matching binary packages.

**Evidence and verification.** source-inspected
[setup.py:python_requires and install_requires](<../../../../setup.py>); [docs/user-guide/installation.md](<../../../user-guide/installation.md>).

Legacy implementation context: [`environment.yml`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/environment.yml), [`requirements.txt`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/requirements.txt).

<a id="issue-277"></a>

### #277: CondaEnvException: Pip failed

[Original issue](https://github.com/bulik/ldsc/issues/277) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **2**.

Workflow tags: Cross-cutting / not assigned to a requested group.

**Trigger and harm.** Conda environment creation fails while pip builds legacy pandas0.20 on macOS.
The documented dependency pin requires an unavailable or incompatible old binary/build toolchain.

**Legacy disposition.** Users report removing the exact pandas pin works; no verified merged package fix.

**LDSC3 assessment — Apparently prevented.** LDSC3 supports Python3.11-3.13 and modern bounded pandas/NumPy versions, avoiding this specific Python2 source-build combination. A compatible wheel or conda package should still be preferred over unrestricted version removal, which can introduce other incompatibilities.

**Recommended prevention.** Offer an isolated supported environment and package install with workflow extras; do not tell users to remove all version bounds.

**Evidence and verification.** source-inspected
[setup.py:python_requires and install_requires](<../../../../setup.py>); [docs/user-guide/installation.md](<../../../user-guide/installation.md>).

Legacy implementation context: [`environment.yml`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/environment.yml), [`requirements.txt`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/requirements.txt).

<a id="issue-287"></a>

### #287: ERROR: Command errored out with exit status 1: /Users/hongyan/opt/anaconda3/envs/ldsc/bin/python -u -c 'import sys, setuptools, tokenize; sys.argv[0] = '"'"'/private/var/folders/9t/965z56gn1lj_tknr10tt47w80000gn/T/pip-install-3eZFqg/pandas/setup.py'"'"'; __file__='"'"'/private/var/folders/9t/965z56gn1lj_tknr10tt47w80000gn/T/pip-install-3eZFqg/pandas/setup.py'"'"';f=getattr(tokenize, '"'"'open'"'"', open)(__file__);code=f.read().replace('"'"'\r\n'"'"', '"'"'\n'"'"');f.close();exec(compile(code, __file__, '"'"'exec'"'"'))' install --record /private/var/folders/9t/965z56gn1lj_tknr10tt47w80000gn/T/pip-record-ltQTy5/install-record.txt --single-version-externally-managed --compile --install-headers /Users/hongyan/opt/anaconda3/envs/ldsc/include/python2.7/pandas Check the logs for full command output.  failed  CondaEnvException: Pip failed

[Original issue](https://github.com/bulik/ldsc/issues/287) · Source state: **open** · Severity: **high** · Accessible comments reviewed: **0**.

Workflow tags: Cross-cutting / not assigned to a requested group.

**Trigger and harm.** Conda pip step fails building pandas for the legacy Python2 environment.
Installation ends before any LDSC command is available.

**Legacy disposition.** Open with no comment or reproducible build diagnostics beyond the final pip failure.

**LDSC3 assessment — Apparently prevented.** LDSC3 modern Python and bounded dependencies eliminate the specific old-pandas requirement, but the original compiler/platform cause cannot be recovered from a summary exit message.
Its setup metadata supports ordinary noneditable installation without legacy script execution.

**Recommended prevention.** Use a supported Python environment and available binary dependencies, and capture the first build error when installation fails.

**Evidence and verification.** source-inspected
[setup.py:python_requires and install_requires](<../../../../setup.py>); [docs/user-guide/installation.md](<../../../user-guide/installation.md>).

Legacy implementation context: [`environment.yml`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/environment.yml), [`requirements.txt`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/requirements.txt).

<a id="issue-26"></a>

### #26: munge_sumstats.py

[Original issue](https://github.com/bulik/ldsc/issues/26) · Source state: **closed** · Severity: **moderate** · Accessible comments reviewed: **45**.

Workflow tags: `munge-sumstats`, `rg`, `baseline-partitioned-h2`.

**Trigger and harm.** An incomplete checkout first causes ImportError; later comments report dependency incompatibilities, allele mismatches and methodological questions.
Users can conflate installation, numerical validity and interpretation failures in one thread.

**Legacy disposition.** Recloning fixed the initial missing module; maintainers later improved allele mismatch messaging and advised compatible dependencies.
No single commit resolves every subthread.

**LDSC3 assessment — Partial coverage.** LDSC3 packaging removes ad hoc script copies and declares dependencies; lazy plotting avoids the old pandas/matplotlib import chain.
Its allele harmonization and regime-specific p-value documentation address later questions.
Remaining lessons include dominant-APOE sensitivity, reference release/sample provenance, mixed-ancestry unsuitability, and interpreting padded annotations as distinct hypotheses.
These are substantive comments, not incidental pasted commands.

**Recommended prevention.** Use an isolated installed environment and a minimal smoke test; document locus sensitivity, two-sided depletion tests and exact resource provenance in their appropriate guides.

**Evidence and verification.** source-inspected
[setup.py](<../../../../setup.py>); [src/ldsc/_regression/runner.py:build_rg_dataset](<../../../../src/ldsc/_regression/runner.py>); [docs/reference/results/partitioned-h2.md](<../../../reference/results/partitioned-h2.md>); [Source](<https://github.com/bulik/ldsc/issues/26#issuecomment-75445614>); [Source](<https://github.com/bulik/ldsc/issues/26#issuecomment-145859412>).

Legacy implementation context: [`environment.yml`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/environment.yml), [`requirements.txt`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/requirements.txt).

<a id="issue-36"></a>

### #36: Some test failures

[Original issue](https://github.com/bulik/ldsc/issues/36) · Source state: **closed** · Severity: **moderate** · Accessible comments reviewed: **2**.

Workflow tags: Cross-cutting / not assigned to a requested group.

**Trigger and harm.** Pandas assertions begin checking Series names; later reports also show singular statistical fixtures and index-alignment errors. A passing tutorial does not certify installation or estimator tests, and a test-only mismatch is not itself a scientific failure.

**Legacy disposition.** Maintainer fixed assertion-name expectations; subsequent comment reopens related unresolved failures as #64.
No single verified commit resolves every report.

**LDSC3 assessment — Partial coverage.** LDSC3 uses pytest and unittest with purpose-specific fixtures rather than nose.
It has exact alignment and known-value tests, but no CI or cross-platform matrix is configured.
Root local checks can establish this environment only; Windows support and all supported dependency combinations remain unverified unless actually tested.

**Recommended prevention.** Separate metadata-only test failures from numerical failures and publish the tested environment; add a small CI matrix when release preparation is authorized.

**Evidence and verification.** source-inspected
[pyproject.toml](<../../../../pyproject.toml>); [tests/regression/test_kernel_regression.py](<../../../../tests/regression/test_kernel_regression.py>); [tests/package/test_package_layout.py](<../../../../tests/package/test_package_layout.py>); [AGENTS.md:Environment and verification](<../../../../AGENTS.md>).

Legacy implementation context: [`environment.yml`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/environment.yml), [`requirements.txt`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/requirements.txt).

<a id="issue-42"></a>

### #42: Run failure

[Original issue](https://github.com/bulik/ldsc/issues/42) · Source state: **closed** · Severity: **moderate** · Accessible comments reviewed: **6**.

Workflow tags: Cross-cutting / not assigned to a requested group.

**Trigger and harm.** SciPy import fails with undefined symbol zlacn2_ before help can run. A broken compiled scientific stack prevents all regression workflows.

**Legacy disposition.** The user confirmed a fresh Anaconda installation fixed the shared-library mismatch; no LDSC algorithm change was needed.

**LDSC3 assessment — Partial coverage.** LDSC3 declares SciPy/NumPy dependencies and supported Python versions, but cannot repair a system LAPACK mismatch. A coherent conda environment or compatible binary wheels avoids the historical manual build combination.
Successful package metadata resolution alone does not prove native libraries import on every platform.

**Recommended prevention.** Recommend one isolated conda route and one python -m pip virtual-environment route, each ending with import/help/pip-check smoke tests.

**Evidence and verification.** source-inspected
[setup.py](<../../../../setup.py>); [environment.yml](<../../../../environment.yml>); [docs/user-guide/installation.md](<../../../user-guide/installation.md>).

Legacy implementation context: [`environment.yml`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/environment.yml), [`requirements.txt`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/requirements.txt).

<a id="issue-55"></a>

### #55: starting off

[Original issue](https://github.com/bulik/ldsc/issues/55) · Source state: **closed** · Severity: **moderate** · Accessible comments reviewed: **9**.

Workflow tags: Cross-cutting / not assigned to a requested group.

**Trigger and harm.** A beginner runs git clone at the Python prompt, later installs dependencies into another interpreter.
Repeated installation attempts and copied scripts produce an incoherent environment.

**Legacy disposition.** Maintainer explains shell versus Python; a clean machine/environment eventually works.
No algorithmic fix is involved.

**LDSC3 assessment — Partial coverage.** LDSC3 supplies installed console/module entry points and separates Python API usage from CLI workflows.
Its current installation still assumes readers know which terminal/environment is active, and Unix line continuations do not work unchanged in every Windows shell.
This is a documentation and onboarding issue rather than missing computational functionality.

**Recommended prevention.** Label code blocks as terminal or Python, show a one-line smoke command, and avoid telling beginners to copy individual package files into their analysis directory.

**Evidence and verification.** source-inspected
[docs/user-guide/installation.md](<../../../user-guide/installation.md>); [docs/user-guide/python-usage.md](<../../../user-guide/python-usage.md>); [README.md](<../../../../README.md>).

Legacy implementation context: [`environment.yml`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/environment.yml), [`requirements.txt`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/requirements.txt).

<a id="issue-86"></a>

### #86: latest release

[Original issue](https://github.com/bulik/ldsc/issues/86) · Source state: **closed** · Severity: **moderate** · Accessible comments reviewed: **2**.

Workflow tags: Cross-cutting / not assigned to a requested group.

**Trigger and harm.** Cluster administrators want a stable release instead of an unversioned development checkout.
Unpinned installs are hard to reproduce and distribute.

**Legacy disposition.** Maintainer says master was maintained with fixes but no new release planned; commit identifiers were suggested as version provenance.

**LDSC3 assessment — Partial coverage.** LDSC3 is a distributable package with a beta version and release runbook, but package metadata alone does not establish a published wheel, stable tag or immutable installer.
An editable checkout is convenient for development; end users and cluster deployments need a pinned source revision or built artifact with tested dependencies.

**Recommended prevention.** For release preparation provide a pinned tag/wheel installation path and actual build/import evidence; avoid a complex release service or claiming an unpublished package is on an index.

**Evidence and verification.** source-inspected
[setup.py](<../../../../setup.py>); [docs/dev/release.md](<../../release.md>); [docs/user-guide/installation.md](<../../../user-guide/installation.md>).

Legacy implementation context: [`environment.yml`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/environment.yml), [`requirements.txt`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/requirements.txt).

<a id="issue-112"></a>

### #112: ./munge_sumstats.py is running for a long time

[Original issue](https://github.com/bulik/ldsc/issues/112) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **3**.

Workflow tags: `munge-sumstats`.

**Trigger and harm.** A large GWAS stalls after the reading line without useful progress.
Users repeatedly rebuild environments without knowing whether parsing, I/O, memory or computation is slow.

**Legacy disposition.** Comments suggest compatible conda environments; later report remains unresolved.
The issue has no proven universal cause.

**LDSC3 assessment — Partial coverage.** LDSC3 provides chunk size controls and phase logging, but file-authoritative progress may be invisible to users watching a terminal or batch job. A configured environment and retained row/chunk counts help distinguish a real stall from slow I/O; no throughput promise follows from passing a tiny fixture.

**Recommended prevention.** Document where to inspect progress and memory, show a small-input smoke run, and diagnose measured bottlenecks before changing dependencies or adding a progress framework.

**Evidence and verification.** source-inspected
[src/ldsc/_sumstats/munge_kernel.py:parse_dat](<../../../../src/ldsc/_sumstats/munge_kernel.py>); [src/ldsc/_runtime/progress.py](<../../../../src/ldsc/_runtime/progress.py>); [docs/reference/logs-and-outputs.md](<../../../reference/logs-and-outputs.md>).

Legacy implementation context: [`environment.yml`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/environment.yml), [`requirements.txt`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/requirements.txt).

<a id="issue-177"></a>

### #177: Tag a release

[Original issue](https://github.com/bulik/ldsc/issues/177) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **2**.

Workflow tags: Cross-cutting / not assigned to a requested group.

**Trigger and harm.** Users need a tagged release and setup.py for conda integration.
Unversioned source installs impair reproducibility and package-manager support.

**Legacy disposition.** Merged PR176 adds setup.py; merge3d0c4464777b2578bf6f13386f0c1c9ab7d55046.
The request for a newer release tag is separate.

**LDSC3 assessment — Partial coverage.** LDSC3 has setuptools package metadata, console entry points and a provisional2.0b0 version.
That solves source installation, but does not demonstrate a published wheel, current tag, PyPI release or maintained conda recipe; the release guide is procedural evidence only.

**Recommended prevention.** Provide a noneditable version-pinned source-install recipe now; release wheels/tags only through an explicitly authorized release.

**Evidence and verification.** source-inspected
[setup.py:python_requires and install_requires](<../../../../setup.py>); [docs/user-guide/installation.md](<../../../user-guide/installation.md>).

Legacy implementation context: [`environment.yml`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/environment.yml), [`requirements.txt`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/requirements.txt).

<a id="issue-358"></a>

### #358: munge_sumstats.py extremely slow runtime in the provided conda environment

[Original issue](https://github.com/bulik/ldsc/issues/358) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **7**.

Workflow tags: `munge-sumstats`.

**Trigger and harm.** Munging takes days in the old environment, while a Python-3 port completes the same data rapidly.
Avoidable runtime and confusing installed-script invocation obstruct routine preprocessing.

**Legacy disposition.** Commenter publishes belowlab/ldsc and a PyPI build; several users confirm that workaround.
This is a separate distribution, not an upstream commit fix.

**LDSC3 assessment — Partial coverage.** LDSC3 uses modern bounded dependencies and streams inputs, removing the reported old stack.
It still concatenates retained chunks, so genome-wide runtime and memory cannot be guaranteed from a small test.

**Recommended prevention.** Publish a small timed munging example and distinguish the intended Git branch from similarly named PyPI distributions; investigate large-file performance with measured inputs.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[setup.py](<../../../../setup.py>); [src/ldsc/_sumstats/munge_kernel.py:697](<../../../../src/ldsc/_sumstats/munge_kernel.py#L697>); [Source](<https://github.com/bulik/ldsc/issues/358#issuecomment-1292082956>).

Legacy implementation context: [`environment.yml`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/environment.yml), [`requirements.txt`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/requirements.txt).

<a id="issue-372"></a>

### #372: No module named bz2

[Original issue](https://github.com/bulik/ldsc/issues/372) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **0**.

Workflow tags: `munge-sumstats`.

**Trigger and harm.** The Python-2 interpreter was built without the bz2 standard-library extension.
Even help/import can fail before data processing.

**Legacy disposition.** Open with no verified fix; installing a package named bz2 would not restore an absent interpreter extension.

**LDSC3 assessment — Partial coverage.** LDSC3 declares supported Python versions and uses gzip/bz2 support in input resolution. A custom Python build without compression modules can still fail; modern version bounds alone do not prove these modules are present.

**Recommended prevention.** Use a complete conda or standard CPython installation and add a small compression-import/read smoke check to installation troubleshooting.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[src/ldsc/_sumstats/munge_input.py:15](<../../../../src/ldsc/_sumstats/munge_input.py#L15>); [src/ldsc/_sumstats/munge_input.py:118](<../../../../src/ldsc/_sumstats/munge_input.py#L118>); [docs/user-guide/installation.md](<../../../user-guide/installation.md>).

Legacy implementation context: [`environment.yml`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/environment.yml), [`requirements.txt`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/requirements.txt).

<a id="issue-445"></a>

### #445: "Pip subprocess error" during installation of LDSC

[Original issue](https://github.com/bulik/ldsc/issues/445) · Source state: **closed** · Severity: **moderate** · Accessible comments reviewed: **3**.

Workflow tags: Cross-cutting / not assigned to a requested group.

**Trigger and harm.** Conda finishes environment creation but its pip subprocess cannot import pip.
Users repeatedly recreate a broken environment without reaching package installation.

**Legacy disposition.** A commenter recommends interpreter-specific ensurepip; the author confirms it fixed the issue.
No upstream environment commit is identified.

**LDSC3 assessment — Partial coverage.** LDSC3 lists pip explicitly and supports modern Python, removing the old pinned stack.
An incomplete interpreter environment can still lack pip, so python -m pip is the correct diagnostic boundary.

**Recommended prevention.** Document python -m pip --version and an interpreter-specific ensurepip/conda repair when needed; provide a venv installation route without requiring development extras.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[environment.yml](<../../../../environment.yml>); [setup.py](<../../../../setup.py>); [docs/user-guide/installation.md](<../../../user-guide/installation.md>).

Legacy implementation context: [`environment.yml`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/environment.yml), [`requirements.txt`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/requirements.txt).

<a id="issue-443"></a>

### #443: Error  # terms in --M must match # of LD Scores in --ref-ld

[Original issue](https://github.com/bulik/ldsc/issues/443) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **0**.

Workflow tags: Cross-cutting / not assigned to a requested group.

**Trigger and harm.** MTAG reports more/fewer count terms than LD-score columns only for chromosome-split resources.
An external integration may pair incompatible counts and predictors.

**Legacy disposition.** Open with no traceback or verified fix; a concatenation workaround is reported by the author but does not establish scientific correctness.

**LDSC3 assessment — Intentionally unsupported.** LDSC3 validates count shape and ordered metadata in its own LD-score directory, but MTAG vendors a different interface and will not adopt those checks automatically.

**Recommended prevention.** Diagnose the MTAG resource/parser boundary separately and require matching per-chromosome annotation columns and counts; no arbitrary manual --M adjustment.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[src/ldsc/_regression/estimators.py:263](<../../../../src/ldsc/_regression/estimators.py#L263>); [src/ldsc/_ldscore/legacy_converter.py](<../../../../src/ldsc/_ldscore/legacy_converter.py>); [Source](<https://github.com/bulik/ldsc/issues/443>).

Legacy implementation context: [`environment.yml`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/environment.yml), [`requirements.txt`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/requirements.txt).

<a id="issue-206"></a>

### #206: did not install it and have some errors below.

[Original issue](https://github.com/bulik/ldsc/issues/206) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **0**.

Workflow tags: Cross-cutting / not assigned to a requested group.

**Trigger and harm.** Git clone fails with curl18, unexpected hangup and earlyEOF.
The repository cannot be obtained, so no package code executes.

**Legacy disposition.** Open without a demonstrated workaround or package fix.

**LDSC3 assessment — Indeterminate.** This is a source-download/network failure rather than evidence that LDSC3 contains a numerical or packaging bug.
LDSC3 still relies on obtaining source for the documented installation route; a versioned distribution artifact would reduce dependence on a full Git transfer.

**Recommended prevention.** Suggest retrying a verified source download or a pinned source archive; avoid arbitrary global Git buffer/network changes.

**Evidence and verification.** source-inspected
[setup.py:python_requires and install_requires](<../../../../setup.py>); [docs/user-guide/installation.md](<../../../user-guide/installation.md>).

Legacy implementation context: [`environment.yml`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/environment.yml), [`requirements.txt`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/requirements.txt).

<a id="issue-423"></a>

### #423: ./ldsc.py -h Illegal instruction (core dumped)

[Original issue](https://github.com/bulik/ldsc/issues/423) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **1**.

Workflow tags: Cross-cutting / not assigned to a requested group.

**Trigger and harm.** Even help exits with an illegal CPU instruction.
The process terminates before Python can report a normal package error.

**Legacy disposition.** A commenter suggests checking architecture and importing dependencies separately; no culprit or verified fix is established.

**LDSC3 assessment — Indeterminate.** LDSC3 removes pybedtools but still depends on compiled NumPy/SciPy/PyArrow libraries.
Those can be incompatible with an architecture or binary environment, so the refactor does not prove this class of crash impossible.

**Recommended prevention.** Diagnose architecture and individual dependency imports in the same interpreter; install compatible wheels or conda builds rather than changing numerical package code.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[setup.py](<../../../../setup.py>); [docs/user-guide/installation.md](<../../../user-guide/installation.md>).

Legacy implementation context: [`environment.yml`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/environment.yml), [`requirements.txt`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/requirements.txt).

<a id="issue-446"></a>

### #446: munge_sumstats issue: Usecols do not match names

[Original issue](https://github.com/bulik/ldsc/issues/446) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **1**.

Workflow tags: `munge-sumstats`.

**Trigger and harm.** CELLECT preprocessing adds a SNP column and a vendored munger fails with usecols mismatch. A lengthy shell pipeline obscures whether header editing, whitespace or fork code introduced the mismatch.

**Legacy disposition.** A commenter identifies the CELLECT fork and directs the question to its tracker; no root cause or fix is proven.

**LDSC3 assessment — Indeterminate.** LDSC3 has explicit column-hint validation and no dependence on CELLECT vendored modules.
Without the attached file/log content and fork revision, the exact failure cannot be declared fixed.

**Recommended prevention.** Show detected versus requested headers before usecols parsing and recommend a minimal direct input example; qualify external-fork support.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[src/ldsc/_sumstats/munge_input.py:287](<../../../../src/ldsc/_sumstats/munge_input.py#L287>); [tests/sumstats/test_munge_flag_contracts.py](<../../../../tests/sumstats/test_munge_flag_contracts.py>); [Source](<https://github.com/bulik/ldsc/issues/446>).

Legacy implementation context: [`environment.yml`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/environment.yml), [`requirements.txt`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/requirements.txt).

<a id="issue-492"></a>

### #492: Segmentation fault

[Original issue](https://github.com/bulik/ldsc/issues/492) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **0**.

Workflow tags: `ldscore`.

**Trigger and harm.** Legacy score generation on chr12 segfaults after genotype filtering. A process-level native crash loses the run and cannot be diagnosed from a Python traceback alone.

**Legacy disposition.** Open with no architecture, binary-version or minimized reproducer evidence.

**LDSC3 assessment — Indeterminate.** LDSC3 validates PLINK header/size and uses a different modern execution stack, but compiled libraries and malformed files can still crash. A refactor is not proof this original native failure is resolved.

**Recommended prevention.** Minimize the PLINK input and isolate reader/native-library versions and architecture; record process failure context, without attributing it to a specific library prematurely.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[src/ldsc/_runtime/paths.py:535](<../../../../src/ldsc/_runtime/paths.py#L535>); [src/ldsc/_refpanel/adapters.py](<../../../../src/ldsc/_refpanel/adapters.py>); [docs/reference/parallel-execution.md](<../../../reference/parallel-execution.md>).

Legacy implementation context: [`environment.yml`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/environment.yml), [`requirements.txt`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/requirements.txt).

<a id="issue-30"></a>

### #30: unexpected keyword argument 'subset'

[Original issue](https://github.com/bulik/ldsc/issues/30) · Source state: **closed** · Severity: **moderate** · Accessible comments reviewed: **2**.

Workflow tags: `h2`.

**Trigger and harm.** The installed pandas predates the API used by legacy LDSC.
Regression aborts before fitting even when inputs are valid.

**Legacy disposition.** Updating pandas fixed the user report; this is a dependency compatibility resolution rather than a scientific code fix.

**LDSC3 assessment — Apparently prevented.** LDSC3 declares a modern supported Python/dependency range and uses current pandas APIs.
Installing it through the declared package metadata prevents this particular ancient-pandas combination. A manually assembled environment can still violate requirements, so the install guide should consistently invoke pip through the selected interpreter.

**Recommended prevention.** Use python -m pip and pip check in a clean environment; document supported versions instead of advising an unbounded upgrade inside an existing analysis environment.

**Evidence and verification.** source-inspected
[setup.py](<../../../../setup.py>); [requirements.txt](<../../../../requirements.txt>); [environment.yml](<../../../../environment.yml>); [docs/user-guide/installation.md](<../../../user-guide/installation.md>).

Legacy implementation context: [`environment.yml`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/environment.yml), [`requirements.txt`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/requirements.txt).

<a id="issue-39"></a>

### #39: ./ldsc.py -h ERROR

[Original issue](https://github.com/bulik/ldsc/issues/39) · Source state: **closed** · Severity: **moderate** · Accessible comments reviewed: **11**.

Workflow tags: Cross-cutting / not assigned to a requested group.

**Trigger and harm.** Python 2.6 executes dictionary-comprehension syntax written for Python 2.7. A basic help command fails before a useful compatibility message appears.

**Legacy disposition.** The user confirms Python 2.7.8 resolved the failure.
Earlier advice to reset a checkout was unrelated to the final cause.

**LDSC3 assessment — Apparently prevented.** LDSC3 metadata requires Python >=3.11,<3.14 and console entry points bind to the installation environment.
Its own modern syntax cannot run on Python 2.6.
The transferable risk is invoking another interpreter after installation, not a need for compatibility with obsolete Python releases.

**Recommended prevention.** Show python -m pip and python -m ldsc from the same activated environment; avoid destructive git resets as an installation diagnostic.

**Evidence and verification.** source-inspected
[setup.py](<../../../../setup.py>); [src/ldsc/__main__.py](<../../../../src/ldsc/__main__.py>); [docs/user-guide/installation.md](<../../../user-guide/installation.md>).

Legacy implementation context: [`environment.yml`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/environment.yml), [`requirements.txt`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/requirements.txt).

<a id="issue-64"></a>

### #64: Errors in tests with nosetests

[Original issue](https://github.com/bulik/ldsc/issues/64) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **5**.

Workflow tags: Cross-cutting / not assigned to a requested group.

**Trigger and harm.** Nose tests fail across environments, including a simulated SE discrepancy and pandas index errors.
Users cannot distinguish a stochastic fixture failure from a production numerical defect.

**Legacy disposition.** Verified commit 346e6ce797c6ba1dd1beab0aa18497a4e1380778 (PR #113) increases N_REP from 200 to 500 in test/test_sumstats.py.
Maintainer confirms this resolves the stochastic slow-test failure; later user confirms the declared environment passes.

**LDSC3 assessment — Apparently prevented.** LDSC3 has deterministic golden/known-value tests, current pandas-compatible merge paths and pytest/unittest runners.
It does not ship the old reduced-replicate test.
Root local focused checks cover the current numerical/format paths, but cross-version/platform and full statistical calibration are separate validation claims.

**Recommended prevention.** Keep stochastic checks adequately replicated or make expectations analytic; publish actual test results instead of inferring complete equivalence from a passing tutorial.

**Evidence and verification.** source-inspected
[tests/regression/test_kernel_regression.py](<../../../../tests/regression/test_kernel_regression.py>); [tests/fixtures/README.md](<../../../../tests/fixtures/README.md>); [pyproject.toml](<../../../../pyproject.toml>); [Source](<https://github.com/bulik/ldsc/pull/113>); [Source](<https://github.com/bulik/ldsc/commit/346e6ce797c6ba1dd1beab0aa18497a4e1380778>).

Legacy implementation context: [`environment.yml`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/environment.yml), [`requirements.txt`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/requirements.txt).

<a id="issue-90"></a>

### #90: error launching munge_sumstats.py

[Original issue](https://github.com/bulik/ldsc/issues/90) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **5**.

Workflow tags: `munge-sumstats`.

**Trigger and harm.** Anaconda3 or another interpreter runs mixed-indentation legacy code.
Cryptic syntax/indentation errors prevent even help.

**Legacy disposition.** Several users resolve the failure with the declared legacy conda environment; an indentation patch alone was not a full Python 3 port.

**LDSC3 assessment — Apparently prevented.** LDSC3 explicitly targets Python 3.11-3.13 and installs console/module entry points.
It does not share the Python 2 indentation path.
The enduring risk is choosing an interpreter different from the installed environment, so environment activation and executable checks matter more than instructing users to patch indentation manually.

**Recommended prevention.** Give one clean environment command and both console/module smoke checks; never advise editing installed package indentation as a setup step.

**Evidence and verification.** source-inspected
[setup.py](<../../../../setup.py>); [environment.yml](<../../../../environment.yml>); [docs/user-guide/installation.md](<../../../user-guide/installation.md>); [tests/package/test_cli_help.py](<../../../../tests/package/test_cli_help.py>).

Legacy implementation context: [`environment.yml`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/environment.yml), [`requirements.txt`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/requirements.txt).

<a id="issue-93"></a>

### #93: Errors running

[Original issue](https://github.com/bulik/ldsc/issues/93) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **1**.

Workflow tags: `munge-sumstats`.

**Trigger and harm.** Pandas 0.21 interprets an old mask expression as integer labels -1/-2. A tutorial fails during incompatible-allele masking despite valid inputs.

**Legacy disposition.** A user reports downgrading to pandas 0.20.1 fixes it; no verified source fix is cited.

**LDSC3 assessment — Apparently prevented.** LDSC3 no longer uses the old allele_merge routine that padded HM3 rows.
It uses explicit identifier cleanup, restriction and downstream harmonization with modern pandas bounds.
The specific deprecated-mask path is absent; this does not certify every dtype combination because nullable masks and missing alleles still require focused tests.

**Recommended prevention.** Keep explicit Boolean masks and tests with missing/incompatible allele pairs; recommend the declared environment rather than historical pandas downgrades.

**Evidence and verification.** source-inspected
[src/ldsc/_sumstats/munge_kernel.py](<../../../../src/ldsc/_sumstats/munge_kernel.py>); [src/ldsc/_regression/runner.py:build_rg_dataset](<../../../../src/ldsc/_regression/runner.py>); [setup.py](<../../../../setup.py>).

Legacy implementation context: [`environment.yml`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/environment.yml), [`requirements.txt`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/requirements.txt).

<a id="issue-104"></a>

### #104: Error converting summary statistics

[Original issue](https://github.com/bulik/ldsc/issues/104) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **3**.

Workflow tags: `munge-sumstats`, `rg`.

**Trigger and harm.** A working legacy install fails after upgrading pandas.
Environment drift destroys a reproducible analysis despite unchanged data.

**Legacy disposition.** Maintainers point to pandas 0.21 and a pinned conda environment; user confirms reverting to 0.18.1 works.
No universal source repair is documented.

**LDSC3 assessment — Apparently prevented.** LDSC3 uses modern declared dependency ranges and no old allele_merge padding routine.
It remains vulnerable to untested dependency releases within those bounds if no CI matrix checks them.
The original linked attachment was not independently reproduced, so classification uses maintainer diagnosis and source-path comparison.

**Recommended prevention.** Recommend isolated pinned analysis environments and test supported dependency endpoints before release, rather than instructing users to downgrade arbitrary packages independently.

**Evidence and verification.** source-inspected
[setup.py](<../../../../setup.py>); [environment.yml](<../../../../environment.yml>); [src/ldsc/_sumstats/munge_kernel.py](<../../../../src/ldsc/_sumstats/munge_kernel.py>); [docs/user-guide/installation.md](<../../../user-guide/installation.md>).

Legacy implementation context: [`environment.yml`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/environment.yml), [`requirements.txt`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/requirements.txt).

<a id="issue-121"></a>

### #121: ERROR converting summary statistics

[Original issue](https://github.com/bulik/ldsc/issues/121) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **3**.

Workflow tags: `munge-sumstats`.

**Trigger and harm.** Legacy allele_merge writes through a mask converted to negative integer labels.
Preprocessing fails despite millions of passing SNPs.

**Legacy disposition.** Maintainer confirms known pandas-version issue and points to the newly released environment.yml; user found #104 workaround.

**LDSC3 assessment — Apparently prevented.** LDSC3 replaces padded allele_merge with explicit retained-identity processing and modern dependency requirements. A keep list called w_hm3 containing 29 million records is also not equivalent to the standard HM3 resource: filenames do not establish scientific role.
The current interface must preserve that distinction when users supply custom restrictions.

**Recommended prevention.** Use the package environment, show active restriction provenance/counts, and reject malformed masks rather than reproducing historical dependency pins.

**Evidence and verification.** source-inspected
[src/ldsc/_sumstats/munge_kernel.py](<../../../../src/ldsc/_sumstats/munge_kernel.py>); [environment.yml](<../../../../environment.yml>); [docs/reference/munge-sumstats.md](<../../../reference/munge-sumstats.md>).

Legacy implementation context: [`environment.yml`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/environment.yml), [`requirements.txt`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/requirements.txt).

<a id="issue-124"></a>

### #124: pandas dataframe.sort is now removed

[Original issue](https://github.com/bulik/ldsc/issues/124) · Source state: **closed** · Severity: **moderate** · Accessible comments reviewed: **1**.

Workflow tags: Cross-cutting / not assigned to a requested group.

**Trigger and harm.** A dependency update removes a deprecated method used for genomic ordering.
Regression cannot run, or a careless workaround can break ordered jackknife blocks.

**Legacy disposition.** Maintainers recommend a pinned conda environment; no production source fix is identified in the thread.

**LDSC3 assessment — Apparently prevented.** LDSC3 uses current sort_values/shared genomic ordering rather than DataFrame.sort.
Ordering remains scientifically important for block jackknife and LD windows, so replacing calls must preserve stable ties and chromosome order.
Current dependency bounds plus tests are stronger evidence than suggesting arbitrary old versions.

**Recommended prevention.** Keep stable genomic-order tests and validate supported pandas releases before publication; no compatibility shim for removed ancient APIs is needed.

**Evidence and verification.** source-inspected
[src/ldsc/_core/row_alignment.py](<../../../../src/ldsc/_core/row_alignment.py>); [tests/regression/test_regression_workflow.py:test_build_dataset_orders_merged_rows_by_ld_score_genomic_order](<../../../../tests/regression/test_regression_workflow.py>); [setup.py](<../../../../setup.py>).

Legacy implementation context: [`environment.yml`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/environment.yml), [`requirements.txt`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/requirements.txt).

<a id="issue-132"></a>

### #132: Error running munge_sumstats.py

[Original issue](https://github.com/bulik/ldsc/issues/132) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **1**.

Workflow tags: `munge-sumstats`.

**Trigger and harm.** Munging fails in allele_merge with labels [-1,-2].
Repeated reports obscure a shared dependency incompatibility.

**Legacy disposition.** User recognizes #104/#93 as the same issue and stops; no new fix is provided.

**LDSC3 assessment — Apparently prevented.** LDSC3 has no historical padded allele_merge path and uses modern package metadata.
This issue belongs in the same remedy cluster as #93/#104/#121 while remaining a separate frequency count.
Treating all repeated reports as independent scientific bugs would overstate distinct mechanisms.

**Recommended prevention.** Cross-link the common environment/mask diagnosis and retain one shared prevention action rather than duplicate implementation proposals.

**Evidence and verification.** source-inspected
[src/ldsc/_sumstats/munge_kernel.py](<../../../../src/ldsc/_sumstats/munge_kernel.py>); [setup.py](<../../../../setup.py>); [docs/user-guide/installation.md](<../../../user-guide/installation.md>).

Legacy implementation context: [`environment.yml`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/environment.yml), [`requirements.txt`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/requirements.txt).

<a id="issue-148"></a>

### #148: Tutorial Cell type specific analyses fails

[Original issue](https://github.com/bulik/ldsc/issues/148) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **1**.

Workflow tags: `munge-sumstats`.

**Trigger and harm.** A successful raw read fails in legacy allele_merge with [-2,-1] labels.
Users attribute a preprocessing dependency failure to the cell-type method.

**Legacy disposition.** Open; later follow-up asks for a solution.
Same established family as #93/#104, not a new demonstrated cell-type defect.

**LDSC3 assessment — Apparently prevented.** LDSC3 does not use that padded allele_merge routine and declares a Python 3 environment.
Default summary-statistics restriction and downstream allele harmonization have separate roles.
The issue is tagged munging only because the failure occurs before LD scores or query regression; the tutorial title alone does not make it an LDSC-SEG failure.

**Recommended prevention.** Cross-link the shared compatibility mechanism and run the native munger smoke test before debugging the cell-type model.

**Evidence and verification.** source-inspected
[src/ldsc/_sumstats/munge_kernel.py](<../../../../src/ldsc/_sumstats/munge_kernel.py>); [setup.py](<../../../../setup.py>); [docs/user-guide/installation.md](<../../../user-guide/installation.md>).

Legacy implementation context: [`environment.yml`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/environment.yml), [`requirements.txt`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/requirements.txt).

<a id="issue-157"></a>

### #157: KeyError in munge_sumstats step

[Original issue](https://github.com/bulik/ldsc/issues/157) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **2**.

Workflow tags: `munge-sumstats`.

**Trigger and harm.** A custom small merge list is followed by negative-index KeyError in allele_merge.
Users blame keep-list size instead of a known mask/dependency issue.

**Legacy disposition.** Comment suggests pandas below .21; no maintainer-confirmed new fix.

**LDSC3 assessment — Apparently prevented.** LDSC3 restriction retains passing SNPs and does not pad through the old mask assignment.
Explicit custom restrictions are allowed and recorded; small or selected lists still change the regression population and can undermine generalization even after the technical error is gone.

**Recommended prevention.** Group the technical remedy with #93/#104 and separately explain custom-list scientific scope; do not delete restriction simply to avoid an error.

**Evidence and verification.** source-inspected
[src/ldsc/_sumstats/munge_kernel.py:_filter_sumstats_chunk_to_restriction](<../../../../src/ldsc/_sumstats/munge_kernel.py>); [setup.py](<../../../../setup.py>); [docs/reference/munge-sumstats.md](<../../../reference/munge-sumstats.md>).

Legacy implementation context: [`environment.yml`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/environment.yml), [`requirements.txt`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/requirements.txt).

<a id="issue-159"></a>

### #159: munge_sumstats.py error

[Original issue](https://github.com/bulik/ldsc/issues/159) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **3**.

Workflow tags: `munge-sumstats`.

**Trigger and harm.** A user fixes list+dict_keys then hits Boolean masking and exception-formatting incompatibilities.
One-line porting advice creates an untested hybrid interpreter environment.

**Legacy disposition.** A commenter resolves by using the supported legacy Python 2 environment; no complete legacy Python 3 port is established.

**LDSC3 assessment — Apparently prevented.** LDSC3 is an actual Python 3 package with list(flag.keys()) where needed, updated logging, and rewritten munging.
It supports a declared modern version range rather than arbitrary Python 3.7.
No legacy patch should be applied to its private modules as an installation step.

**Recommended prevention.** Point migration users to the installed LDSC3 command and supported environment; do not normalize manually patched legacy scripts as a reliable deployment.

**Evidence and verification.** source-inspected
[src/ldsc/_sumstats/munge_input.py:get_cname_map](<../../../../src/ldsc/_sumstats/munge_input.py>); [src/ldsc/_runtime/workflow_logging.py](<../../../../src/ldsc/_runtime/workflow_logging.py>); [setup.py](<../../../../setup.py>).

Legacy implementation context: [`environment.yml`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/environment.yml), [`requirements.txt`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/requirements.txt).

<a id="issue-161"></a>

### #161: error in running munge_sumstats.py

[Original issue](https://github.com/bulik/ldsc/issues/161) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **0**.

Workflow tags: `munge-sumstats`.

**Trigger and harm.** An unported legacy script on Python 3.7 interprets mask values as integer indices.
The same environment failure is reported as a new munging bug.

**Legacy disposition.** Unanswered duplicate-like continuation of #159; no independent resolution.

**LDSC3 assessment — Apparently prevented.** LDSC3 no longer contains the Python 2 padded allele_merge implementation and requires supported Python 3 versions.
Validating the installed package path and environment is more useful than patching the old assignment.
Keep this issue in coverage/frequency counts but share one remedy cluster with #159 and the pandas-mask reports.

**Recommended prevention.** Provide an interpreter/package-location smoke test and a migration pointer; avoid instructions to edit installed source piecemeal.

**Evidence and verification.** source-inspected
[setup.py](<../../../../setup.py>); [src/ldsc/_sumstats/munge_kernel.py](<../../../../src/ldsc/_sumstats/munge_kernel.py>); [docs/user-guide/installation.md](<../../../user-guide/installation.md>).

Legacy implementation context: [`environment.yml`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/environment.yml), [`requirements.txt`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/requirements.txt).

<a id="issue-202"></a>

### #202: Tab&amp;space mix in munge_sumstats.py

[Original issue](https://github.com/bulik/ldsc/issues/202) · Source state: **closed** · Severity: **moderate** · Accessible comments reviewed: **1**.

Workflow tags: `munge-sumstats`.

**Trigger and harm.** A screenshot reports indentation failure in the legacy munging script.
The script cannot even parse, preventing preprocessing.

**Legacy disposition.** Closed as duplicate of issue90; no fixing commit is identified in this thread.

**LDSC3 assessment — Apparently prevented.** LDSC3 supplies importable Python3 modules and one installed CLI rather than requiring execution of the old Python2 script.
Existing package-import and help checks exercise parsing.
Screenshots were not needed to infer an additional numerical defect; the text identifies an interpreter-level indentation problem.

**Recommended prevention.** Retain import/help smoke checks and direct users to ldsc munge-sumstats instead of copied legacy scripts.

**Evidence and verification.** source-inspected
[setup.py:python_requires and install_requires](<../../../../setup.py>); [docs/user-guide/installation.md](<../../../user-guide/installation.md>).

Legacy implementation context: [`environment.yml`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/environment.yml), [`requirements.txt`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/requirements.txt).

<a id="issue-230"></a>

### #230: munge_sumstats.py indentation error

[Original issue](https://github.com/bulik/ldsc/issues/230) · Source state: **closed** · Severity: **moderate** · Accessible comments reviewed: **0**.

Workflow tags: `munge-sumstats`.

**Trigger and harm.** Executing munge_sumstats.py-h under the chosen interpreter raises inconsistent-tabs-and-spaces.
The script cannot parse and no preprocessing occurs.

**Legacy disposition.** Closed without comments or a fixing commit; similar to issues90and202.

**LDSC3 assessment — Apparently prevented.** LDSC3 uses maintained Python3 modules and installed command dispatch; ordinary package import/help tests exercise parseability.
This prevents reliance on the particular copied legacy script but does not justify describing an unspecified upstream closure as a repaired numerical bug.

**Recommended prevention.** Keep module import and CLI help checks on every supported interpreter; direct users to the package command.

**Evidence and verification.** source-inspected
[setup.py:python_requires and install_requires](<../../../../setup.py>); [docs/user-guide/installation.md](<../../../user-guide/installation.md>).

Legacy implementation context: [`environment.yml`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/environment.yml), [`requirements.txt`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/requirements.txt).

<a id="issue-276"></a>

### #276: munge_sumstats.py: command not found

[Original issue](https://github.com/bulik/ldsc/issues/276) · Source state: **closed** · Severity: **moderate** · Accessible comments reviewed: **2**.

Workflow tags: `munge-sumstats`.

**Trigger and harm.** A beginner invokes munge_sumstats.py by bare name under Windows Ubuntu.
Shell PATH/current-directory rules block onboarding before package code runs.

**Legacy disposition.** Closed without a confirmed solution; commenter suggests ./munge_sumstats.py.

**LDSC3 assessment — Apparently prevented.** LDSC3 installs the ldsc console script and supports python -m ldsc, so users do not need a root script on PATH.
Environment activation remains essential.
The help smoke checks provide an immediate check of interpreter/package discovery before an analysis is attempted.

**Recommended prevention.** Show conda activation and python -m ldsc munge-sumstats as a reliable alternative when the console command is unavailable.

**Evidence and verification.** source-inspected
[setup.py:python_requires and install_requires](<../../../../setup.py>); [docs/user-guide/installation.md](<../../../user-guide/installation.md>).

Legacy implementation context: [`environment.yml`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/environment.yml), [`requirements.txt`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/requirements.txt).

<a id="issue-296"></a>

### #296: munge_sumstats in python3

[Original issue](https://github.com/bulik/ldsc/issues/296) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **1**.

Workflow tags: `munge-sumstats`.

**Trigger and harm.** An external fork needs compressed-header decoding and Boolean-mask changes.
Partial Python3 conversion can make ordinary files fail or subtly change filtering.

**Legacy disposition.** The thread proposes local edits and later points to another fork; neither proves a legacy master fix.

**LDSC3 assessment — Apparently prevented.** LDSC3 is a native installable Python3 package, explicitly decodes compressed header bytes and uses Boolean matching masks.
Those repairs address the two concrete mechanisms without sys.path hacks.
Scientific equivalence should still be tested rather than inferred from successful execution after a port.

**Recommended prevention.** Keep compressed/plain equivalence and Boolean-mask coverage on supported Python versions; distinguish a runnable port from validated numerical behavior.

**Evidence and verification.** source-inspected
[src/ldsc/_sumstats/munge_input.py:_decode_header_line](<../../../../src/ldsc/_sumstats/munge_input.py>); [src/ldsc/_sumstats/munge_kernel.py:_filter_sumstats_chunk_to_restriction](<../../../../src/ldsc/_sumstats/munge_kernel.py>); [tests/sumstats/test_sumstats_munger.py](<../../../../tests/sumstats/test_sumstats_munger.py>).

Legacy implementation context: [`environment.yml`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/environment.yml), [`requirements.txt`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/requirements.txt).

<a id="issue-300"></a>

### #300: Problem with irwls

[Original issue](https://github.com/bulik/ldsc/issues/300) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **1**.

Workflow tags: Cross-cutting / not assigned to a requested group.

**Trigger and harm.** A converted regressions.py imports irwls as if it were an external top-level package.
Users search for a nonexistent dependency or alter sys.path with machine-specific directories.

**Legacy disposition.** Comment recommends a sys.path append workaround; no master fix verified.

**LDSC3 assessment — Apparently prevented.** LDSC3 uses package-relative imports between its _regression modules and exposes public objects from ldsc.
The installed layout removes dependence on the current directory and does not require pip-installing an unrelated irwls project.

**Recommended prevention.** Install the package and use public entry points; retain import tests outside the repository directory.

**Evidence and verification.** source-inspected
[src/ldsc/_regression/estimators.py](<../../../../src/ldsc/_regression/estimators.py>); [src/ldsc/_regression/irwls.py](<../../../../src/ldsc/_regression/irwls.py>); [tests/package/test_package_layout.py](<../../../../tests/package/test_package_layout.py>).

Legacy implementation context: [`environment.yml`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/environment.yml), [`requirements.txt`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/requirements.txt).

<a id="issue-350"></a>

### #350: Pybedtools not available for windows os

[Original issue](https://github.com/bulik/ldsc/issues/350) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **0**.

Workflow tags: `ldscore`.

**Trigger and harm.** Environment creation cannot install pybedtools on Windows.
Users cannot reach annotation generation.

**Legacy disposition.** Open; no workaround in this thread.

**LDSC3 assessment — Apparently prevented.** LDSC3 annotation projection uses NumPy and does not depend on pybedtools or a bedtools executable.
That removes the named blocker, but does not certify the complete package on Windows.

**Recommended prevention.** Document that BED projection has no external BED-tool requirement and add a Windows install/workflow smoke test before claiming Windows support.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[setup.py](<../../../../setup.py>); [src/ldsc/_annotation/interval_projection.py](<../../../../src/ldsc/_annotation/interval_projection.py>); [docs/user-guide/installation.md](<../../../user-guide/installation.md>).

Legacy implementation context: [`environment.yml`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/environment.yml), [`requirements.txt`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/requirements.txt).

<a id="issue-351"></a>

### #351: TabError: inconsistent use of tabs and spaces in indentation

[Original issue](https://github.com/bulik/ldsc/issues/351) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **6**.

Workflow tags: `munge-sumstats`, `rg`.

**Trigger and harm.** Python 3 executes Python-2 scripts and raises TabError and print SyntaxError.
Neither munging nor regression starts.

**Legacy disposition.** Users report activating the documented Python-2 environment as the workaround; no modern upstream fix is shown.

**LDSC3 assessment — Apparently prevented.** LDSC3 declares Python 3.11-3.13 and installs one console entry point.
This removes Python-2 syntax, although running an old cloned script can still select the wrong software.

**Recommended prevention.** Use ldsc and python -m ldsc in every example, and show interpreter/package-location checks for mixed installations.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[setup.py](<../../../../setup.py>); [src/ldsc/cli.py](<../../../../src/ldsc/cli.py>); [docs/user-guide/installation.md](<../../../user-guide/installation.md>).

Legacy implementation context: [`environment.yml`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/environment.yml), [`requirements.txt`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/requirements.txt).

<a id="issue-354"></a>

### #354: TypeError: can only concatenate list (not "dict_keys") to list

[Original issue](https://github.com/bulik/ldsc/issues/354) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **1**.

Workflow tags: `munge-sumstats`.

**Trigger and harm.** A manually converted script adds list and dict_keys and masks the first error with incorrect traceback.format_exc usage.
Munging fails before reading the data.

**Legacy disposition.** Comment supplies list(flag.keys()) as a local workaround and references #159; no merge evidence is supplied.

**LDSC3 assessment — Apparently prevented.** LDSC3 get_cname_map materializes key lists and its installed Python-3 entry point avoids hand-porting scripts.
This addresses the specific expression, not all user-edited copies.

**Recommended prevention.** Keep the supported installation route explicit; retain a column-hint test exercising dictionary-key combinations and preserve original exceptions.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[src/ldsc/_sumstats/munge_input.py:108](<../../../../src/ldsc/_sumstats/munge_input.py#L108>); [src/ldsc/_runtime/workflow_logging.py](<../../../../src/ldsc/_runtime/workflow_logging.py>).

Legacy implementation context: [`environment.yml`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/environment.yml), [`requirements.txt`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/requirements.txt).

<a id="issue-355"></a>

### #355: UnboundLocalError: local variable 'cname_translation' referenced before assignment

[Original issue](https://github.com/bulik/ldsc/issues/355) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **1**.

Workflow tags: `munge-sumstats`.

**Trigger and harm.** A user-modified Python-3 script references cname_translation before assignment.
The wrapper failure obscures the actual configuration and prevents munging.

**Legacy disposition.** A commenter suggests checking indentation; this is not a verified root cause or upstream fix.

**LDSC3 assessment — Apparently prevented.** LDSC3 resolves schema into a ResolvedMungeInput before running the kernel, where cname_translation comes from the request.
The shown local edit is absent, but the original script was not available to reproduce.

**Recommended prevention.** Direct users to the installed entry point rather than manual 2to3 edits; retain the earliest exception in logs.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[src/ldsc/_sumstats/munge_input.py:287](<../../../../src/ldsc/_sumstats/munge_input.py#L287>); [src/ldsc/_sumstats/munge_kernel.py:95](<../../../../src/ldsc/_sumstats/munge_kernel.py#L95>); [src/ldsc/_sumstats/munge_kernel.py:697](<../../../../src/ldsc/_sumstats/munge_kernel.py#L697>).

Legacy implementation context: [`environment.yml`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/environment.yml), [`requirements.txt`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/requirements.txt).

<a id="issue-361"></a>

### #361: make_annot.py_error "AttributeError: 'float' object has no attribute 'isdigit'"

[Original issue](https://github.com/bulik/ldsc/issues/361) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **1**.

Workflow tags: `ldscore`, `ldsc-seg`.

**Trigger and harm.** Gene-coordinate annotation creation sends numeric scalar fields to pybedtools, which calls isdigit on float.
Annotation generation fails even with otherwise usable gene coordinates.

**Legacy disposition.** Comment points to a pybedtools workaround, but explicitly cannot confirm output correctness.

**LDSC3 assessment — Apparently prevented.** LDSC3 uses its own gene-coordinate parsing and NumPy interval projection without pybedtools.
It also validates chromosome scope, so a genome-wide list with a chr22-only grid should fail actionably rather than silently drop genes.

**Recommended prevention.** Provide the current gene coordinate catalog schema and a small successful example; keep malformed numeric coordinates as explicit errors, rather than coercing blindly to strings.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[src/ldsc/_annotation/gene_lists.py](<../../../../src/ldsc/_annotation/gene_lists.py>); [src/ldsc/_annotation/interval_projection.py](<../../../../src/ldsc/_annotation/interval_projection.py>); [tests/annotation/test_annotate_preflight.py](<../../../../tests/annotation/test_annotate_preflight.py>).

Legacy implementation context: [`environment.yml`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/environment.yml), [`requirements.txt`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/requirements.txt).

<a id="issue-369"></a>

### #369: Getting module parse has no attibute error

[Original issue](https://github.com/bulik/ldsc/issues/369) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **31**.

Workflow tags: `h2`, `munge-sumstats`, `baseline-partitioned-h2`.

**Trigger and harm.** A Python-3 hand-port imports an unrelated parse package; later edits create circular imports and format errors.
Users install unrelated similarly named packages or edit imports repeatedly, while help output misleadingly appears to work.

**Legacy disposition.** The final commenter reports success after correcting example/input formats; no upstream code fix or single verified repair is identified.

**LDSC3 assessment — Apparently prevented.** LDSC3 uses package-relative private imports and an installed entry point, avoiding bare parse/regressions imports.
The thread also shows why help-only checks do not validate a real analysis.
Referenced screenshots/archive were not needed to infer the posted import failures.

**Recommended prevention.** Provide a tiny complete workflow after installation and a module-location diagnostic; do not advise pip-installing generic parse, jackknife or regressions packages.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[src/ldsc/cli.py](<../../../../src/ldsc/cli.py>); [src/ldsc/_regression/runner.py](<../../../../src/ldsc/_regression/runner.py>); [tests/package/test_package_layout.py](<../../../../tests/package/test_package_layout.py>).

Legacy implementation context: [`environment.yml`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/environment.yml), [`requirements.txt`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/requirements.txt).

<a id="issue-379"></a>

### #379: LDSC requires pandas version &gt;= 0.17.0

[Original issue](https://github.com/bulik/ldsc/issues/379) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **0**.

Workflow tags: Cross-cutting / not assigned to a requested group.

**Trigger and harm.** An import reports that pandas must be at least 0.17; no full environment or traceback is supplied.
Users cannot start the legacy program and may downgrade unrelated environments.

**Legacy disposition.** Open; no verified root cause or fix in the accessible record.

**LDSC3 assessment — Apparently prevented.** LDSC3 declares modern pandas bounds through package metadata rather than relying on the historical version gate.
That removes the old minimum-version check but does not establish the state of an existing user environment.

**Recommended prevention.** Use python -m pip in an isolated supported environment and report actual interpreter, pandas version and package source in troubleshooting.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[setup.py](<../../../../setup.py>); [requirements.txt](<../../../../requirements.txt>); [environment.yml](<../../../../environment.yml>).

Legacy implementation context: [`environment.yml`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/environment.yml), [`requirements.txt`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/requirements.txt).

<a id="issue-380"></a>

### #380: LDSC requires pandas version &gt;= 0.17.0

[Original issue](https://github.com/bulik/ldsc/issues/380) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **0**.

Workflow tags: Cross-cutting / not assigned to a requested group.

**Trigger and harm.** An import reports that pandas must be at least 0.17; no full environment or traceback is supplied.
Users cannot start the legacy program and may downgrade unrelated environments.

**Legacy disposition.** Open; no verified root cause or fix in the accessible record.

**LDSC3 assessment — Apparently prevented.** LDSC3 declares modern pandas bounds through package metadata rather than relying on the historical version gate.
That removes the old minimum-version check but does not establish the state of an existing user environment.

**Recommended prevention.** Use python -m pip in an isolated supported environment and report actual interpreter, pandas version and package source in troubleshooting.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[setup.py](<../../../../setup.py>); [requirements.txt](<../../../../requirements.txt>); [environment.yml](<../../../../environment.yml>).

Legacy implementation context: [`environment.yml`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/environment.yml), [`requirements.txt`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/requirements.txt).

<a id="issue-400"></a>

### #400: TypeError: a bytes-like object is required, not 'str'

[Original issue](https://github.com/bulik/ldsc/issues/400) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **1**.

Workflow tags: `munge-sumstats`.

**Trigger and harm.** A Python-3 conversion writes text to a legacy binary log handle.
Munging fails while writing its header, before its scientific logic runs.

**Legacy disposition.** A commenter notes Python-3 support is absent and suggests 2to3; no verified full repair is supplied.

**LDSC3 assessment — Apparently prevented.** LDSC3 uses Python logging and workflow-managed text files instead of the shown print-to-binary logger.
The obsolete Python 3.6/dependency combination is outside supported bounds.

**Recommended prevention.** Keep text encoding explicit in log writers and recommend the supported package installation rather than manual interpreter conversion.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[src/ldsc/_runtime/workflow_logging.py](<../../../../src/ldsc/_runtime/workflow_logging.py>); [setup.py](<../../../../setup.py>); [docs/reference/logs-and-outputs.md](<../../../reference/logs-and-outputs.md>).

Legacy implementation context: [`environment.yml`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/environment.yml), [`requirements.txt`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/requirements.txt).

<a id="issue-417"></a>

### #417: Can't run ldsc under screen at Linux

[Original issue](https://github.com/bulik/ldsc/issues/417) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **2**.

Workflow tags: Cross-cutting / not assigned to a requested group.

**Trigger and harm.** The same legacy script works in one shell but fails with Python-3 print syntax under screen.
Environment mismatch prevents batch execution and encourages unsafe source edits.

**Legacy disposition.** A commenter identifies PATH/interpreter mismatch; changing print syntax is proposed by another but is not a complete port.

**LDSC3 assessment — Apparently prevented.** LDSC3 requires Python 3 and offers python -m ldsc, which binds the chosen interpreter.
Different shells can still activate different installations, so interpreter provenance remains relevant.

**Recommended prevention.** Show activation inside the batch/session script plus python --version and package-source checks; do not prescribe editing print statements.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[setup.py](<../../../../setup.py>); [src/ldsc/__main__.py](<../../../../src/ldsc/__main__.py>); [docs/user-guide/installation.md](<../../../user-guide/installation.md>).

Legacy implementation context: [`environment.yml`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/environment.yml), [`requirements.txt`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/requirements.txt).

<a id="issue-432"></a>

### #432: Is python 2.7 no longer available to download using code in ldsc tutorial?

[Original issue](https://github.com/bulik/ldsc/issues/432) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **5**.

Workflow tags: Cross-cutting / not assigned to a requested group.

**Trigger and harm.** An M3 Mac cannot resolve Python2.7, pybedtools0.7 and bitarray0.8 for osx-arm64.
Installation is impossible through the old environment recipe on that platform.

**Legacy disposition.** Comments explain old Python support ended and point to conversion/another port; no bulik master environment fix is proven.

**LDSC3 assessment — Apparently prevented.** LDSC3 uses supported Python3 and modern compiled dependencies without pybedtools.
That removes the obsolete recipe, while a clean ARM install remains a separate verification task.

**Recommended prevention.** Provide a modern isolated environment and a small installed workflow check; avoid emulation or source conversion unless a demonstrated platform need remains.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[setup.py](<../../../../setup.py>); [environment.yml](<../../../../environment.yml>); [docs/user-guide/installation.md](<../../../user-guide/installation.md>).

Legacy implementation context: [`environment.yml`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/environment.yml), [`requirements.txt`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/requirements.txt).

<a id="issue-440"></a>

### #440: make_annot.py does not run

[Original issue](https://github.com/bulik/ldsc/issues/440) · Source state: **closed** · Severity: **moderate** · Accessible comments reviewed: **2**.

Workflow tags: `ldscore`.

**Trigger and harm.** Legacy make_annot hangs importing pybedtools, then another interpreter fails subtracting one from a string BP. A help-only install check misses annotation-specific dependencies and malformed coordinates.

**Legacy disposition.** Author confirms their BIM was malformed; no fix for the separate pybedtools import hang is established.

**LDSC3 assessment — Apparently prevented.** LDSC3 removes pybedtools and validates BIM positions before interval projection.
This addresses both named pathways structurally, while the original cluster environment was not reproduced.

**Recommended prevention.** Keep a complete small annotation smoke test and actionable numeric-position errors; treat the two failures as separate causes.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[src/ldsc/_annotation/bim_grid.py](<../../../../src/ldsc/_annotation/bim_grid.py>); [src/ldsc/_runtime/paths.py:529](<../../../../src/ldsc/_runtime/paths.py#L529>); [src/ldsc/_annotation/interval_projection.py](<../../../../src/ldsc/_annotation/interval_projection.py>).

Legacy implementation context: [`environment.yml`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/environment.yml), [`requirements.txt`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/requirements.txt).

<a id="issue-442"></a>

### #442: Syntax Error

[Original issue](https://github.com/bulik/ldsc/issues/442) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **4**.

Workflow tags: Cross-cutting / not assigned to a requested group.

**Trigger and harm.** User believes Python3 is active but posted python --version is 2.7 and the source has Python-3 print syntax.
Mixed source/interpreter versions prevent startup.

**Legacy disposition.** A commenter identifies the mismatch; no final successful run or upstream fix is shown.

**LDSC3 assessment — Apparently prevented.** LDSC3 requires supported Python3 and installs its entry point from that interpreter.
An old environment can still invoke the wrong executable, so the actual version matters more than the prompt name.

**Recommended prevention.** Use python -m ldsc and python -m pip from one environment, and report interpreter/source locations before advising source edits.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[setup.py](<../../../../setup.py>); [src/ldsc/__main__.py](<../../../../src/ldsc/__main__.py>); [docs/user-guide/installation.md](<../../../user-guide/installation.md>).

Legacy implementation context: [`environment.yml`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/environment.yml), [`requirements.txt`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/requirements.txt).

<a id="issue-453"></a>

### #453: Compute LD scores fails: no module named ldscore.ldscore

[Original issue](https://github.com/bulik/ldsc/issues/453) · Source state: **closed** · Severity: **moderate** · Accessible comments reviewed: **1**.

Workflow tags: `ldscore`.

**Trigger and harm.** CELLECT works on one OS but fails in a container importing ldscore. A broken vendored package layout prevents score computation.

**Legacy disposition.** Author reports the cause was a misnamed init file and closes the issue; no upstream code fix is involved.

**LDSC3 assessment — Apparently prevented.** LDSC3 installs a src-layout package with tested package contents and does not rely on copying a sibling ldscore directory.
Containers can still omit files if built from an incomplete source copy.

**Recommended prevention.** Install the built package inside the container and run a tiny score workflow; avoid advising a pip package named ldscore as a fix for a missing vendored module.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[setup.py](<../../../../setup.py>); [MANIFEST.in](<../../../../MANIFEST.in>); [tests/package/test_package_layout.py](<../../../../tests/package/test_package_layout.py>).

Legacy implementation context: [`environment.yml`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/environment.yml), [`requirements.txt`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/requirements.txt).

<a id="issue-490"></a>

### #490: Error running munge_sumstats.py

[Original issue](https://github.com/bulik/ldsc/issues/490) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **2**.

Workflow tags: `munge-sumstats`.

**Trigger and harm.** A local munger contains a malformed j-for-i assignment absent from the tracked legacy source. A bug in an unknown modified copy is mistaken for a current package failure and masks its error with traceback handling.

**Legacy disposition.** A commenter identifies the non-upstream line; the reporter plans to ask administrators for an update, without a confirmed final run.

**LDSC3 assessment — Apparently prevented.** LDSC3 uses installed package modules and identity-based filtering instead of this NaN-padding allele_merge assignment.
That does not certify unknown site-modified installations.

**Recommended prevention.** Include source revision and module location in bug reports; reproduce on an unmodified installation before proposing package code changes.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[src/ldsc/_sumstats/munge_kernel.py:458](<../../../../src/ldsc/_sumstats/munge_kernel.py#L458>); [setup.py](<../../../../setup.py>); [Source](<https://github.com/bulik/ldsc/issues/490#issuecomment-3144874418>).

Legacy implementation context: [`environment.yml`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/environment.yml), [`requirements.txt`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/requirements.txt).

<a id="issue-504"></a>

### #504: LDSC installation fails due to incorrect pip executable in Python 2.7 environment

[Original issue](https://github.com/bulik/ldsc/issues/504) · Source state: **open** · Severity: **moderate** · Accessible comments reviewed: **0**.

Workflow tags: Cross-cutting / not assigned to a requested group.

**Trigger and harm.** pip3 tries to install legacy SciPy into the wrong interpreter while the intended environment is Python2.7.
Users pursue an imp-module error instead of fixing interpreter/pip selection.

**Legacy disposition.** Author reports interpreter-specific ensurepip and correct pip selection resolved it; the open state does not negate that workaround.

**LDSC3 assessment — Apparently prevented.** LDSC3 uses modern supported dependencies and its README invokes python -m pip, binding installation to the active interpreter.
Mixed environments can still happen, but the named obsolete SciPy recipe is absent.

**Recommended prevention.** Keep python -m pip in every recipe and include a short interpreter/pip/source check; offer a standard venv path without requiring manual dependency downgrades.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[README.md](<../../../../README.md>); [setup.py](<../../../../setup.py>); [environment.yml](<../../../../environment.yml>).

Legacy implementation context: [`environment.yml`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/environment.yml), [`requirements.txt`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/requirements.txt).

<a id="issue-397"></a>

### #397: Error with 'python mtag.py -h'

[Original issue](https://github.com/bulik/ldsc/issues/397) · Source state: **closed** · Severity: **low** · Accessible comments reviewed: **0**.

Workflow tags: Cross-cutting / not assigned to a requested group.

**Trigger and harm.** MTAG imports its vendored ldsc_mod and fails on a Python-2 print statement under Python 3.
MTAG cannot start, but changing this package would not update the copy it vendors.

**Legacy disposition.** Closed without an accessible technical resolution; the two issues reproduce the same external-tool traceback.

**LDSC3 assessment — Intentionally unsupported.** LDSC3 supports Python 3, but is not a drop-in replacement for MTAG private modules and does not control its bundled code.
This is an interoperability request rather than a confirmed defect in the current installed package.

**Recommended prevention.** Route the fix to the MTAG distribution or an explicitly tested integration; do not recommend copying LDSC3 private modules into MTAG.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[setup.py](<../../../../setup.py>); [src/ldsc/__init__.py](<../../../../src/ldsc/__init__.py>); [Source](<https://github.com/bulik/ldsc/issues/397>).

Legacy implementation context: [`environment.yml`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/environment.yml), [`requirements.txt`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/requirements.txt).

<a id="issue-398"></a>

### #398: Error with 'python mtag.py -h'

[Original issue](https://github.com/bulik/ldsc/issues/398) · Source state: **closed** · Severity: **low** · Accessible comments reviewed: **0**.

Workflow tags: Cross-cutting / not assigned to a requested group.

**Trigger and harm.** MTAG imports its vendored ldsc_mod and fails on a Python-2 print statement under Python 3.
MTAG cannot start, but changing this package would not update the copy it vendors.

**Legacy disposition.** Closed without an accessible technical resolution; the two issues reproduce the same external-tool traceback.

**LDSC3 assessment — Intentionally unsupported.** LDSC3 supports Python 3, but is not a drop-in replacement for MTAG private modules and does not control its bundled code.
This is an interoperability request rather than a confirmed defect in the current installed package.

**Recommended prevention.** Route the fix to the MTAG distribution or an explicitly tested integration; do not recommend copying LDSC3 private modules into MTAG.

**Evidence and verification.** Source-inspected; original inputs and platform not reproduced.
[setup.py](<../../../../setup.py>); [src/ldsc/__init__.py](<../../../../src/ldsc/__init__.py>); [Source](<https://github.com/bulik/ldsc/issues/398>).

Legacy implementation context: [`environment.yml`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/environment.yml), [`requirements.txt`](https://github.com/bulik/ldsc/blob/2fdeeb3b44379408794154993dbd6101b8946b7e/requirements.txt).

