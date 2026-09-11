# Package attribution and metadata audit

Last updated on: 2026-09-10

Audited `ldsc_py3_restructured`, branch `restructure`, commit `69081b2e0b2235e2ed81a781ba7e32ca3d7ad61b`. This is an audit, not an authorship assignment or a licensing change. Package code and metadata were left unchanged.

**Historical baseline:** the findings and metadata below describe the initial audit. The subsequent user-authorized corrections and remaining verification gaps are recorded in [the remediation record](remediation.md). The current author list and contact are in [setup.py](../../../setup.py) and [CITATION.cff](../../../CITATION.cff); the historical email below is not the current contact.

The main problems are missing distributed license text, a removed upstream copyright notice, inherited contact information, and no central citation or contributor record. Address the first two before the next distribution. The current author's preferred public identity, contact address, citation author order, and exact license expression remain decisions to resolve.

## Verified package metadata

The values below come from [setup.py](../../../setup.py), lines 32–43, and freshly built wheel `METADATA` and source-distribution `PKG-INFO`. Both builds agree; see [build evidence](build-evidence.json).

| Field | Observed value |
| --- | --- |
| Distribution / import / command | `ldsc` / `ldsc` / `ldsc` |
| README title | `ldsc3_Jerry` |
| Version | Source `2.0b`; built artifacts `2.0b0` |
| Summary | `LD Score Regression (LDSC)` |
| Homepage | `http://github.com/abrantesas/ldsc_py3` |
| Authors | Anthony Abrantes, Brendan Bulik-Sullivan and Hilary Finucane |
| Author email | `antshaabr@gmail.com` |
| Maintainer / maintainer email | Absent |
| Declared license | Legacy free-text field `GPLv3` |
| License expression / license files | Absent |
| Project links / long description | Absent |
| README installation repository | `https://github.com/WenbinWu2001/ldsc_py3` |

## Findings

### A1 — High: neither source tree nor distributions contain license text

The tracked tree contains no `LICENSE`, `COPYING`, or equivalent license document. Both generated archives contain zero license or credit files, despite declaring `GPLv3`. `include_package_data=True` does not supply missing documents. Evidence: [setup.py](../../../setup.py), lines 39–43; [build evidence](build-evidence.json), both artifact records.

The local legacy repository contains the full GPLv3 text at `/Users/wenbinwu/Documents_local/Research/SullivanLab/LDSC/repos/ldsc_py2_Bulik_workspace/ldsc_py2_Bulik/LICENSE`. Its sections 4–5 address preserving notices, supplying the license, and identifying modifications; see also the [GPLv3 text published by SPDX](https://spdx.github.io/license-list-data/GPL-3.0.html). Restore the applicable license text and verify inclusion in both distribution formats. This is a concrete distribution gap, not a conclusion about every aspect of license compliance.

`GPLv3` also leaves the `only` versus `or-later` choice unstated. Establish the applicable upstream grant before selecting an SPDX expression; the generic “How to Apply” example inside the GPL document does not establish the project's choice. Modern metadata supports `License-Expression` and `License-File`; the former replaces the legacy `License` field. [PyPA core metadata specification](https://packaging.python.org/en/latest/specifications/core-metadata/#license-expression).

### A2 — High: restructuring removed a copyright notice from retained legacy code

[`src/ldsc/_kernel/formats.py`](../../../src/ldsc/_kernel/formats.py), lines 1–17, has a functional module description but no attribution. Its `get_compression` and `__ID_List_Factory__` retain legacy parsing code at lines 31 and 43. The upstream equivalents are in `/Users/wenbinwu/Documents_local/Research/SullivanLab/LDSC/repos/ldsc_py2_Bulik_workspace/ldsc_py2_Bulik/ldscore/parse.py`, lines 58 and 246; that file credits Brendan Bulik-Sullivan and Hilary Finucane at line 2.

History confirms the omission: commit `d7c8777b9334aecfea671ad37e75912134df9e8c` renamed `ldscore/parse.py` to `src/ldsc/_kernel/formats.py` and replaced the header, deleting the 2014 copyright line. Restore that attribution and add a concise modification notice grounded in the actual history.

Only three of 70 tracked source Python modules contain copyright notices:

| Module | Retained notice |
| --- | --- |
| [`_kernel/_irwls.py`](../../../src/ldsc/_kernel/_irwls.py), line 2 | 2015 Brendan Bulik-Sullivan and Hilary Finucane |
| [`_kernel/_jackknife.py`](../../../src/ldsc/_kernel/_jackknife.py), line 2 | 2014 Brendan Bulik-Sullivan and Hilary Finucane |
| [`_kernel/regression.py`](../../../src/ldsc/_kernel/regression.py), line 2 | 2014 Brendan Bulik-Sullivan and Hilary Finucane |

The other 67 files are not automatically defective: some are new, and some legacy modules never had individual notices. Use provenance to distinguish new and derived files; do not mechanically assign every file to one person or replace historical years with the current year. A complete derivation map remains follow-up work.

### A3 — Medium: current maintainership is missing and support metadata points to an earlier fork

The metadata homepage and email refer to the earlier Python 3 fork, while [README installation instructions](../../../README.md), lines 32 and 72, use the current repository. There is no separate current maintainer, support URL, or contributor roster. The published author list also omits the current development contributor visible in Git history. Evidence: [setup.py](../../../setup.py), lines 36–38.

`git shortlog -sn HEAD` reports 634 commits as `WenbinWu2001`, five as `Wenbin Wu (wenbinwu)`, and 12 as `abrantes`. These identities establish contribution history, not citation order, ownership, current affiliations, or permission to publish a preferred email. Preserve upstream credit and distinguish original LDSC authors, Python 3 port contributors, and current refactoring/maintenance roles. Confirm the preferred public name and contact before changing those fields; do not infer a contact address from Git history. The existing email was inspected but not contacted or verified as active.

### A4 — Medium: users have no central instructions for citing this software and its methods

The README has no citation, authors, acknowledgments, or license section. The tracked tree contains no `CITATION.cff`, bibliography file, or equivalent software citation record. Scientific references exist in individual documents, but they do not tell users how to cite this fork.

Restore a README citation section and add a software citation record with the chosen title, credited authors, repository URL, and release identifier. GitHub supports `CITATION.cff` and distinguishes a software citation from an optional preferred article citation. [GitHub citation-file documentation](https://docs.github.com/en/repositories/managing-your-repositorys-settings-and-features/customizing-your-repository/about-citation-files).

Keep method citations conditional on the analysis performed. The local legacy `README.md`, “Citation” section, lines 82–107, distinguishes basic LDSC/intercept, genetic correlation, partitioned heritability, continuous annotations, the HE relationship, and LD Hub. The current [guided tutorial](../../wiki/guided-tutorial.md), line 278, also cites Finucane et al. 2018 for cell-type analysis. These are source-supported starting points; LD Hub and other optional references should not become universal citation requirements. No software DOI, ORCID, or publication author order was established by this audit.

### A5 — Medium: bundled data have incomplete attribution and their README is not shipped

Both archives contain the two compressed HM3 tables and six BED files, but omit [`src/ldsc/data/readme.txt`](../../../src/ldsc/data/readme.txt). The packaging patterns at [setup.py](../../../setup.py), line 43, select `.tsv.gz` and `.bed` only.

The data README explains schema and how to derive the compact reference, but does not identify the curated map's original datasets, curation author, source releases, or reuse terms. Region provenance is better documented: [region presets](../../current/region-exclusion-presets.md), “Source of the coordinates,” and [the builder](../../../tools/regions/build_region_beds.py), “Sources,” identify UCSC tracks and Alkes-group genetic maps. The [HM3 validation document](../../current/bundled-hm3-map-coordinate-validation.md), “Official External Sources,” records coordinate-validation sources; this does not establish the provenance of every bundled column.

Create a small resource attribution record covering source datasets/releases, curation, relevant citations, and applicable redistribution terms, and package it. This finding identifies missing evidence; it does not assert that redistribution is prohibited or that the software license automatically governs upstream data.

### A6 — Low: software identity and release identification are weak for citations

The README title, repository name, and distribution name differ without an explanation. The built version `2.0b0` is valid normalization of `2.0b`, not a build error. However, users have no canonical software citation title, and the current CLI parser and shared logging header expose no package-version banner. See [README](../../../README.md), line 1; [setup.py](../../../setup.py), lines 33–35; [`_logging.py`](../../../src/ldsc/_logging.py), `_WorkflowLoggingContext._write_header`, line 420; and [`cli.py`](../../../src/ldsc/cli.py), `build_parser`.

Choose a display/citation title and explain its relationship to `ldsc`; changing the import or distribution name is not necessary for this repair. Define release/version updates and a way to record the software version used for an analysis. [The release runbook](../../release.md), “Procedure” and “Guards,” currently addresses branch synchronization but contains no citation/license/credit verification step.

### A7 — Low: historical contributor instructions conflict with the current authorship rule

The [2026-06-11 implementation plan](../../plans/2026-06-11-overlap-aware-partitioned-h2-plan.md), line 16, prescribes a tool-generated coauthor credit, conflicting with [AGENTS.md](../../../AGENTS.md), “Change discipline.” Git history contains 41 `Co-Authored-By` trailer lines; the count is not a verified roster of people or publication authors. Retire the stale instruction in a future documentation cleanup and avoid deriving formal authorship automatically from trailers. Rewriting historical commits is unnecessary for this audit.

## Recommended cleanup and decisions

1. Restore the applicable license document and the confirmed missing upstream notice; verify both archives include required attribution.
2. Confirm the current maintainer's public name and contact, the software citation title and author order, and the applicable SPDX license expression. Keep legal copyright, software authorship, maintainer contact, and paper authorship distinct.
3. Update metadata links and contact fields; add README license, credit, support, and citation sections plus a `CITATION.cff` and a contributor/provenance record.
4. Document and ship bundled-data attribution; review other derived files against history before adopting a consistent header policy.
5. Add release checks for metadata, citation validity, license/credit inclusion, and version consistency.

## Coverage and validation

Inventoried 363 tracked files, scanned all 70 source Python modules and relevant metadata/documentation/notebook text, inspected selected attribution history, and compared notices and citation instructions with the local legacy repository. No HPC access or inspection of the original local `main` worktree was needed. The pre-existing untracked `docs/audits/flag-cleanup/` directory was excluded and preserved.

Built fresh source and wheel distributions from a complete tracked working-tree snapshot in a scratch directory, offline, using separate `setuptools.build_meta.build_sdist` and `build_wheel` calls with setuptools 82.0.1 and wheel 0.47.0. Both builds succeeded. Examined `PKG-INFO`, `METADATA`, archive members, and bundled documentation; results are saved in [build-evidence.json](build-evidence.json). These are local builds, not an audit of any already published release. No numerical test suite was run because package behavior was not changed. Only this report and its evidence file were added.
