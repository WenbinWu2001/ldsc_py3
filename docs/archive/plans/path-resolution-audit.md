# Path Resolution Audit

Date: 2026-04-20

Purpose: capture path-resolution and implicit file-discovery behaviors that should be revisited when finalizing the package. This note focuses on user-facing runtime code, not test-only helpers.

## Summary

The explicit liftover-chain fix removed the main repo-layout-dependent bundled-resource lookup from the reference-panel builder. I did not find another active runtime module that still walks `__file__` parents or auto-discovers bundled resources from `resources/`.

There are still several other implicit path-resolution behaviors in the package. Most are legacy compatibility features rather than hidden bundled-file lookup, but they can still mask mistakes or make behavior less explicit than desired.

## Findings

### 1. Implicit default gene-coordinate file

Files:
- [src/ldsc/cli.py](/Users/wenbinwu/Documents_local/Research/SullivanLab/LDSC/repos/ldsc_py3_Jerry_workspace/ldsc_py3_restructured/src/ldsc/cli.py:149)
- [src/ldsc/_kernel/annotation.py](/Users/wenbinwu/Documents_local/Research/SullivanLab/LDSC/repos/ldsc_py3_Jerry_workspace/ldsc_py3_restructured/src/ldsc/_kernel/annotation.py:730)

Behavior:
- `--gene-coord-file` still defaults to `ENSG_coord.txt`.
- That file is not resolved from an explicit user path.
- It is looked up relative to the current working directory at runtime.

Risk:
- silent dependence on caller working directory
- easy to break in batch jobs, worktrees, and installed package usage

Recommendation:
- require `--gene-coord-file` explicitly for this workflow
- if a default is kept, make it a documented package asset with explicit resolution logic rather than a bare filename

### 2. General path token expansion in workflow path resolution

File:
- [src/ldsc/path_resolution.py](/Users/wenbinwu/Documents_local/Research/SullivanLab/LDSC/repos/ldsc_py3_Jerry_workspace/ldsc_py3_restructured/src/ldsc/path_resolution.py:74)

Behavior:
- `substitute_chromosome()` appends `@` if missing
- `resolve_scalar_path()` and `resolve_file_group()` accept:
  - exact paths
  - globs
  - suffix expansion
  - chromosome-suite token expansion when enabled

Examples:
- `token -> token + ".txt"` or similar suffixes
- `token -> token@ -> token1`, `token2`, ...

Risk:
- typo-tolerant behavior can hide user mistakes
- exact input is not always the file that gets opened

Recommendation:
- decide package policy explicitly:
  - keep current token language as a compatibility feature, or
  - add a strict mode, or
  - gradually require exact paths in new workflows

### 3. Legacy filename and compression guessing in formats kernel

File:
- [src/ldsc/_kernel/formats.py](/Users/wenbinwu/Documents_local/Research/SullivanLab/LDSC/repos/ldsc_py3_Jerry_workspace/ldsc_py3_restructured/src/ldsc/_kernel/formats.py:52)

Behavior:
- `sub_chr()` appends `@` if missing
- `get_present_chrs()` scans chromosomes by globbing `sub_chr(fh, chrom) + '.*'`
- `which_compression()` tries `.bz2`, then `.gz`, then plain

Risk:
- legacy-compatible but not explicit
- can make debugging harder when the wrong file variant is present

Recommendation:
- keep for legacy regression if needed
- consider documenting these as legacy-only semantics
- for newer workflows, prefer explicit file paths and explicit compression

### 4. Prefix and optional chromosome-file resolution in LD-score kernel

Files:
- [src/ldsc/_kernel/ldscore.py](/Users/wenbinwu/Documents_local/Research/SullivanLab/LDSC/repos/ldsc_py3_Jerry_workspace/ldsc_py3_restructured/src/ldsc/_kernel/ldscore.py:611)
- [src/ldsc/_kernel/ldscore.py](/Users/wenbinwu/Documents_local/Research/SullivanLab/LDSC/repos/ldsc_py3_Jerry_workspace/ldsc_py3_restructured/src/ldsc/_kernel/ldscore.py:636)
- [src/ldsc/_kernel/ldscore.py](/Users/wenbinwu/Documents_local/Research/SullivanLab/LDSC/repos/ldsc_py3_Jerry_workspace/ldsc_py3_restructured/src/ldsc/_kernel/ldscore.py:645)

Behavior:
- `resolve_prefixed_file()` tries:
  - exact path
  - glob
  - suffix expansion
- `resolve_chr_files()` applies implicit chromosome substitution
- `resolve_optional_chr_files()` silently ignores missing chromosome-specific files

Risk:
- optional missing-file behavior can hide typos or incomplete inputs
- suffix and prefix expansion reduce explicitness

Recommendation:
- review whether missing optional chromosome files should emit warnings
- consider reserving silent skip behavior for truly optional query inputs only

### 5. Reference-panel builder still allows non-exact token resolution for maps and PLINK prefixes

File:
- [src/ldsc/ref_panel_builder.py](/Users/wenbinwu/Documents_local/Research/SullivanLab/LDSC/repos/ldsc_py3_Jerry_workspace/ldsc_py3_restructured/src/ldsc/ref_panel_builder.py:53)

Behavior:
- PLINK prefixes are resolved through `resolve_plink_prefix_group(...)`
- genetic maps are resolved through `resolve_file_group(...)`
- both allow token patterns such as globs and chromosome-suite placeholders

Note:
- this is no longer bundled-resource guessing
- however, it still means a user-provided token may not be an exact path

Recommendation:
- decide whether this workflow should remain flexible or move to strict exact-path requirements
- if flexibility is kept, document it clearly in the tutorial and CLI help

## What Was Fixed

Files:
- [src/ldsc/config.py](/Users/wenbinwu/Documents_local/Research/SullivanLab/LDSC/repos/ldsc_py3_Jerry_workspace/ldsc_py3_restructured/src/ldsc/config.py)
- [src/ldsc/ref_panel_builder.py](/Users/wenbinwu/Documents_local/Research/SullivanLab/LDSC/repos/ldsc_py3_Jerry_workspace/ldsc_py3_restructured/src/ldsc/ref_panel_builder.py)
- [src/ldsc/_kernel/ref_panel_builder.py](/Users/wenbinwu/Documents_local/Research/SullivanLab/LDSC/repos/ldsc_py3_Jerry_workspace/ldsc_py3_restructured/src/ldsc/_kernel/ref_panel_builder.py)

Completed change:
- reference-panel building now requires explicit liftover chain paths
- there is no fallback to guessed bundled chain files

Required flags:
- `--liftover-chain-hg19-to-hg38` when `--source-genome-build hg19`
- `--liftover-chain-hg38-to-hg19` when `--source-genome-build hg38`

## Suggested Follow-up During Package Finalization

1. Decide whether new workflows should allow token expansion or require exact paths.
2. Remove the `ENSG_coord.txt` default or replace it with explicit documented asset handling.
3. Review silent optional-file skipping in LD-score query/frequency path resolution.
4. Separate strict modern workflow behavior from legacy compatibility behavior in docs and CLI help.
