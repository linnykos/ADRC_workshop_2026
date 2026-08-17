# Git hooks (public/ repo)

`public/` is its own git repository, separate from the `ADRC_workshop_2026/` folder that
contains it, and it is the **only** part of the workshop tree that is version-controlled
and pushed to GitHub (`linnykos/ADRC_workshop_2026`). The parent folder is not a git repo,
so its `.githooks/` and `.gitignore` have no effect in here — git only reads those at or
below its own repo root. That is why this directory exists as a near-copy of the parent's.

`pre-commit` is **byte-identical to `../../.githooks/pre-commit`** on purpose: `diff` the two
to detect drift, and apply any fix to both.

## What's here

- **`pre-commit`** — blocks any staged file **≥ 50 MB**. Large files bloat git history
  permanently (git keeps every blob forever, even after a later delete) and GitHub warns at
  50 MB and hard-rejects at 100 MB.

## Install (run once per clone, from `public/`)

```sh
git config core.hooksPath .githooks
```

This is stored in `public/.git/config`, which is **not** version-controlled — so it does not
travel with a clone or a `git pull`. Anyone who clones this repo, and Kevin on any second
machine, must run it again. Verify with `git config --get core.hooksPath`.

If the hook doesn't fire, check it is executable: `chmod +x .githooks/pre-commit`.

## Why this repo in particular

The realistic hazards here are all much larger than 50 MB, so the hook is a backstop behind
`.gitignore`, not the first line of defense:

- **`*_cache/`** — knitr caches. Week 4's two tutorials produce **3.0 GB** and **1.6 GB**.
- **`seaad_microglia.RData`** — 525 MB, over five times GitHub's hard limit. Tutorials
  download it into `tempdir()` at runtime; it must never be committed.
- **`*_files/`** — extracted figures, small, but pointless in git since the knitted HTML is
  self-contained.

All of these are covered by `.gitignore`. The hook catches the case where something slips
past those patterns — a renamed data file, a `git add -f`, a new artifact type from a future
week — before it reaches a push that GitHub would reject anyway.

## Usage

- **Raise the limit for one commit:** `MAX_MB=100 git commit ...`
- **Keep the file out of git:** add the pattern to `.gitignore`, then
  `git restore --staged <file>`
- **Bypass once (sparingly):** `git commit --no-verify`

## Verify it works

```sh
dd if=/dev/zero of=bigfile.bin bs=1m count=60
git add bigfile.bin
git commit -m "test"      # should be BLOCKED
git restore --staged bigfile.bin && rm bigfile.bin
```
