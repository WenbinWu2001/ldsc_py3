# Release: squash-merge `restructure` into `main`

Last updated on: 2026-09-10

`main` keeps a squash-only history. `restructure` is the authoritative branch:
every release is a single `Squash merge restructure` commit whose tree is
identical to the `restructure` tip.

## Why not `git merge --squash`

`main` is a chain of independent squash commits and shares no real ancestry with
`restructure`. A plain `git merge --squash restructure` therefore uses a stale
3-way merge base: it raises spurious conflicts and re-adds files that
`restructure` deleted. Reconstruct the tree deterministically instead.

## Procedure

Before branch synchronization, verify the distribution metadata in `setup.py` against the README and `CITATION.cff`: author order, Wenbin Wu as maintainer, only his email as the author contact, repository URL, and release version. Check that source distributions and wheels contain the complete `LICENSE`, `NOTICE`, and packaged `data/readme.txt` and `data/ATTRIBUTION.txt`; the source distribution must also contain `CITATION.cff`. Review derived-code notices and the unresolved upstream SPDX and resource-provenance questions in `NOTICE` and `src/ldsc/data/ATTRIBUTION.txt`. Do not infer an SPDX variant or dataset license to satisfy a packaging check.

Run in the `main` worktree (`../ldsc_py3_Jerry`). Commit all work on
`restructure` first, and make sure the `main` worktree is clean.

```bash
cd /Users/wenbinwu/Documents_local/Research/SullivanLab/LDSC/repos/ldsc_py3_Jerry_workspace/ldsc_py3_Jerry
git fetch origin
git tag -f backup/main-pre-release main                  # safety net

git reset --hard origin/main                             # restore clean squash history
git diff --binary origin/main restructure | git apply --index --whitespace=nowarn

git diff restructure                                     # MUST be empty: tree == restructure
git commit -m "Squash merge restructure"
git push origin main                                     # clean fast-forward
```

## Guards

- If `git diff restructure` is non-empty, STOP -- do not commit. The tree must
  match `restructure` exactly.
- The push is a fast-forward of `origin/main`; never force-push.
- After verifying the release, drop the backup tag: `git tag -d backup/main-pre-release`.
