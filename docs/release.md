# Release: squash-merge `restructure` into `main`

`main` keeps a squash-only history. `restructure` is the authoritative branch:
every release is a single `Squash merge restructure` commit whose tree is
identical to the `restructure` tip.

## Why not `git merge --squash`

`main` is a chain of independent squash commits and shares no real ancestry with
`restructure`. A plain `git merge --squash restructure` therefore uses a stale
3-way merge base: it raises spurious conflicts and re-adds files that
`restructure` deleted. Reconstruct the tree deterministically instead.

## Procedure

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
