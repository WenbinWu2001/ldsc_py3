# Package attribution corrections

Last updated on: 2026-09-10

Scope: user-authorized changes to `ldsc_py3_restructured` on `restructure`, based on commit `69081b2`. No sibling codebase or HPC workspace was accessed during this change. Upstream public declarations were read to verify the requested license restoration. The [initial audit](README.md) and [initial build evidence](build-evidence.json) remain historical snapshots, not current author/contact metadata.

## Implemented

- [setup.py](../../../setup.py), `setup()`, declares Wenbin Wu, Anthony Abrantes, Brendan Bulik-Sullivan, and Hilary Finucane in that order. Wenbin Wu is also the maintainer. Only `drkwu999@gmail.com` is used as the author/maintainer contact. Repository and issue links now identify the current fork.
- [LICENSE](../../../LICENSE) contains the complete 35,142-byte GPLv3 document downloaded from [upstream LDSC](https://raw.githubusercontent.com/bulik/ldsc/master/LICENSE), without editorial changes. Explicit packaging declarations include `LICENSE` and [NOTICE](../../../NOTICE) in source distributions and wheels.
- [README](../../../README.md), “Authors and maintainer,” “License,” and “Citation,” separates maintainer contact, legal terms, upstream credit, and academic citation. Its existing software title `ldsc3_Jerry` is retained; package/import/CLI names and version are unchanged. [CITATION.cff](../../../CITATION.cff) records the same four software authors, only Wenbin Wu's email, and normalized package version `2.0b0`. No publication date, DOI, or ORCID was invented.
- Existing 2014/2015 notices remain in `_irwls.py`, `_jackknife.py`, and `regression.py`. The parser's 2014 notice is restored in [`formats.py`](../../../src/ldsc/_kernel/formats.py). The 2014–2019 upstream banner credit and Anthony Abrantes's 2024 Python 3 credit are restored as comments in [`sumstats_munger.py`](../../../src/ldsc/_kernel/sumstats_munger.py), without restoring obsolete runtime banner code.
- Nine derived source modules now record upstream paths and dated modifications: those five modules plus `_kernel/ldscore.py`, `_kernel/plink_bed.py`, `_sumstats_input.py`, and `column_inference.py`. Credit follows the extracted munging helpers and header-alias logic. Source notices distinguish upstream work from subsequent changes; no current-maintainer copyright replaces an upstream notice.
- [Resource attribution](../../../src/ldsc/data/ATTRIBUTION.txt) and the existing [data README](../../../src/ldsc/data/readme.txt) now ship in both archive types. The [release runbook](../../release.md), “Procedure,” includes metadata, attribution, citation, and license checks.

## History supporting restored notices

The parser notice was removed from `refactor/src/ldsc/_kernel/formats.py` by `18c2d31` on 2026-04-15; `d7c8777` on 2026-04-29 also shows the notice replacement when comparing the original parser with its final package location. Munging commit `c9d9a9c` on 2026-05-02 removed the runtime banner containing upstream copyright and Anthony Abrantes's adaptation credit. `8e513b9` on 2026-09-10 extracted raw input handling into `_sumstats_input.py`. These facts were inspected in this repository's own Git history.

## Verification still open

**Exact SPDX variant:** upstream [README](https://github.com/bulik/ldsc/blob/master/README.md), “License,” [setup.py](https://github.com/bulik/ldsc/blob/master/setup.py), `setup()`, and [ldsc.py](https://github.com/bulik/ldsc/blob/master/ldsc.py), `MASTHEAD`, declare GPLv3 without explicitly stating `only` or `or-later`. The generic application example in the license document does not establish a project-specific choice. Metadata therefore retains `license="GPLv3"`; `CITATION.cff` omits its optional SPDX license field. Further upstream grant evidence is needed before choosing either identifier. See [NOTICE](../../../NOTICE), “License verification.”

**Bundled datasets:** the [resource record](../../../src/ldsc/data/ATTRIBUTION.txt), sections 1–5, distinguishes verified sources from missing provenance. UCSC's [data policy](https://genome.ucsc.edu/license/) permits reuse from UCSC's perspective but notes possible original-producer restrictions. The original HM3 curation inputs and their column-level lineage are not established here. The Alkes-group map filenames are recorded, but their exact releases, URLs, and data-specific terms are not. Neither the software's GPL declaration nor coordinate validation resolves those gaps. The maintainer was asked for the missing source records; no dataset was relicensed or changed.

## Validation

- `python -m pytest -q tests/test_package_layout.py`: **32 passed, 20 subtests passed**, 12.48 seconds, using `ldsc3-dev`.
- All nine edited source modules compiled, and each full executable AST matched its `HEAD` version exactly. Changes are comments only; `setup.py` compiled separately.
- `CITATION.cff` passed YAML parsing and JSON Schema validation against the official CFF 1.2.0 schema. Author order, contact, and absence of other authors' email fields were checked explicitly. Validation dependencies were installed only under `/private/tmp`, outside the package environment.
- A fresh source distribution was built from a snapshot containing tracked working-tree files and the new distribution files. A wheel was then built from that extracted source distribution using setuptools 82.0.1 and wheel 0.47.0, offline. Both archives passed checks for author and maintainer metadata, the sole contact address, version, GPLv3 declaration, absence of an assumed SPDX expression, exact license/notice bytes, and included data documentation. The source distribution also contains `CITATION.cff`. See [remediation build evidence](remediation-build-evidence.json).
- All eight bundled data files matched their source bytes in both archives and matched their original Git contents. No numerical code, resource values, dependencies, or command contracts changed. No publication, commit, or release was performed.
