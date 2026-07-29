# Releasing the `gamblr-collective` container

Maintainer-only process for building and publishing the Docker image that
mirrors `envs/gamblr_collective.yaml` plus all six GAMBLR packages. Not
something package users ever need to run themselves -- this is what backs
LCR-modules' Apptainer install mode, alongside the conda install mode that
`envs/gamblr_collective.yaml` already covers directly.

## Background

The image is built entirely from files already in this repo:

- `envs/gamblr_collective.yaml` provisions every CRAN/Bioconductor
  dependency as a conda binary (see that file's own header for why system
  libraries are deliberately left unpinned).
- `docker/install_gamblr.R` then installs the six GAMBLR packages from
  GitHub on top, via `remotes::install_github(..., dependencies = FALSE)`.
- `docker/Dockerfile` wires the two together.

Because the image is built *from* the conda env file rather than a
separately maintained package list, it can't drift out of sync with it --
any change to `envs/gamblr_collective.yaml` is automatically picked up by
the next image build.

Apptainer/Singularity can pull this image directly (`apptainer pull
docker://ghcr.io/morinlab/gamblr-collective:<tag>`) -- there's no separate
Apptainer-specific build step or definition file.

## Steps

1. **Make sure PRs are merged to `master`** on all six repos first. A
   normal release build installs every package from `master` by default --
   this is not the mechanism for testing open PR branches (see Testing,
   below).

2. **Bump `DESCRIPTION`'s `Version`** in this repo if it hasn't already
   been bumped as part of the release the container is meant to
   accompany. The image tag defaults to this version.

3. **Commit and push the version bump.** `release_container.sh` refuses to
   run otherwise -- same reasoning as `GAMBLR.data`'s
   `data-raw/release_mutations_db.R`: the tag is derived from the version
   on disk, so an uncommitted bump would tag an image nobody else's
   checkout can reproduce.

4. **Run the build+push via GitHub Actions** (recommended, and the
   verified path -- see below). No local Docker/disk-space/credentials
   needed, since it runs on GitHub's own runners and authenticates to GHCR
   with the built-in `GITHUB_TOKEN`. From the Actions tab, run "Release
   gamblr-collective container" (leave all ref inputs at their `master`
   default for a normal release), or via the CLI:
   ```bash
   gh workflow run release_container.yaml --repo morinlab/GAMBLR.open
   ```
   This builds `ghcr.io/morinlab/gamblr-collective:<version>`, pushes it,
   and also moves the floating `latest` tag to it (since a plain run with
   no ref overrides is, by definition, a master-based release build).

   *Verified end-to-end* against both currently-open PRs (#128, #16) --
   see Testing, below, for exactly how. A first attempt failed on an
   `install_github()` timeout (fixed in `install_gamblr.R`); the retry
   succeeded in 12m36s.

   **Fallback: running `docker/release_container.sh` locally** is the
   same script, unchanged, so it should work identically -- useful for
   quick local Dockerfile iteration without pushing a real build, or if
   Actions isn't available. Just not independently re-verified the way the
   Actions path now has been. Log in to GHCR with a GitHub PAT that has
   `write:packages` scope first:
   ```bash
   echo "$GHCR_PAT" | docker login ghcr.io -u <your-github-username> --password-stdin
   ./docker/release_container.sh
   ```

5. **Verify** by pulling the image fresh and running the same smoke test
   used to validate the conda environment:
   ```bash
   docker run --rm ghcr.io/morinlab/gamblr-collective:<version> \
     Rscript -e 'library(GAMBLR.open); meta <- get_gambl_metadata(); muts <- get_coding_ssm(these_samples_metadata = head(meta, 5)); stopifnot(nrow(muts) > 0); cat("OK\n")'
   ```
   Or the Apptainer equivalent:
   ```bash
   apptainer exec docker://ghcr.io/morinlab/gamblr-collective:<version> \
     Rscript -e '...'
   ```

## Testing pre-merge PR branches

Set the relevant ref inputs when dispatching the workflow, e.g. to test
both currently-open PRs together (`gamblr_data_ref` and `gamblr_open_ref`
set to `rmorin-dev`, tag set to something like `test-pr128-pr16`), either
from the Actions tab's "Run workflow" form, or:
```bash
gh workflow run release_container.yaml --repo morinlab/GAMBLR.open \
  -f tag=test-pr128-pr16 -f gamblr_data_ref=rmorin-dev -f gamblr_open_ref=rmorin-dev
```
(The equivalent local invocation is `GAMBLR_DATA_REF=rmorin-dev
GAMBLR_OPEN_REF=rmorin-dev ./docker/release_container.sh
test-pr128-pr16`, but this is the exact command that was actually used to
validate the current Dockerfile/install script -- prefer it unless you
have a specific reason to run locally instead.)

This pushes only the custom tag (`test-pr128-pr16` here) and never
touches `latest`, so it can't be mistaken for a real release by anyone
pulling the image without an explicit tag.

## Troubleshooting

- **"DESCRIPTION has uncommitted changes"** -- commit and push the version
  bump first (step 3).
- **`docker push` auth error** -- your GHCR login isn't set up, or the PAT
  is missing `write:packages` scope (step 4).
- **Install fails inside the image but works via `envs/gamblr_collective.yaml`
  alone** -- almost certainly one of the six packages' own `Imports:`
  version constraints isn't satisfied by what's on `master` yet; check each
  repo's `DESCRIPTION` for pins added against packages that haven't been
  bumped/merged yet.
- **`install_github()` fails with a download timeout** -- GitHub's
  on-demand tarball generation can occasionally exceed R's default 60s
  `download.file()` timeout (hit this on the very first test build).
  `install_gamblr.R` already sets `options(timeout = 600)` and retries each
  package with backoff, so a single flaky attempt shouldn't fail the
  build; a *persistent* timeout suggests a real GitHub-side issue, not
  this script.
- **The image builds and pushes fine, but nobody else can pull it** -- the
  package (not the repo) has its own separate visibility setting, private
  by default for anything pushed via `GITHUB_TOKEN`, regardless of the
  linked repo's own visibility. Check "Package settings" on the package's
  GitHub page (`github.com/morinlab/GAMBLR.open/pkgs/container/gamblr-collective`)
  for a visibility toggle. If it's grayed out, the `morinlab` org itself
  has public-package creation disabled under
  `github.com/organizations/morinlab/settings/packages` ("Package
  creation" section) -- an org owner needs to enable it there first, one
  time, before the per-package toggle becomes available. Confirm anonymous
  pullability afterward with:
  ```bash
  curl -s "https://ghcr.io/token?scope=repository:morinlab/gamblr-collective:pull"
  ```
  A real token back (not a 403) means it's genuinely public.
