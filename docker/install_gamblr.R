# install_gamblr.R
#
# Installs the six GAMBLR packages on top of the gamblr_collective conda
# environment (already provisioned by the Dockerfile via
# envs/gamblr_collective.yaml). dependencies = FALSE is required: every
# CRAN/Bioconductor dependency is already provided as a conda binary, and
# letting remotes re-resolve them would mean recompiling from source inside
# the image, defeating the point of building it this way.
#
# Order matters: GAMBLR.open (and, to a lesser extent, GAMBLR.utils/
# GAMBLR.viz/GAMBLR.predict) declare version-pinned Imports on the others,
# which R checks at install time -- so each package must already be present
# before anything that depends on it installs.
#
# Refs default to "master" (a normal release build, once PRs are merged);
# override via env vars for pre-merge testing of open PR branches.

# R's download.file() defaults to a 60s timeout, which GitHub's on-demand
# tarball generation for install_github() can occasionally exceed (seen in
# practice: a GAMBLR.data@rmorin-dev fetch failed at exactly 60.7s during
# an early test build). 600s gives ample headroom.
options(timeout = 600)

refs <- list(
  GAMBLR.data    = Sys.getenv("GAMBLR_DATA_REF", "master"),
  GAMBLR.helpers = Sys.getenv("GAMBLR_HELPERS_REF", "master"),
  GAMBLR.utils   = Sys.getenv("GAMBLR_UTILS_REF", "master"),
  GAMBLR.viz     = Sys.getenv("GAMBLR_VIZ_REF", "master"),
  GAMBLR.predict = Sys.getenv("GAMBLR_PREDICT_REF", "master"),
  GAMBLR.open    = Sys.getenv("GAMBLR_OPEN_REF", "master")
)

install_with_retry <- function(repo, max_tries = 3) {
  for (i in seq_len(max_tries)) {
    ok <- tryCatch({
      remotes::install_github(repo, dependencies = FALSE, upgrade = "never")
      TRUE
    }, error = function(e) {
      if (i < max_tries) {
        message("Install of ", repo, " failed (attempt ", i, "/", max_tries,
                "): ", conditionMessage(e), " -- retrying in ", 5 * i, "s...")
        Sys.sleep(5 * i)
        FALSE
      } else {
        stop(e)
      }
    })
    if (isTRUE(ok)) return(invisible())
  }
}

# Installed one at a time, in dependency order (the order `refs` is
# defined in), rather than as a single install_github() vector call --
# GAMBLR.open and several others declare version-pinned Imports on the
# rest, checked at install time, so a failure needs to stop the build
# immediately rather than cascading into five more doomed downstream
# installs (which is what a whole-vector call did in practice: GAMBLR.data
# failing still let remotes plow ahead into GAMBLR.helpers/utils/viz/
# predict/open, all of which then failed too, for a less useful log).
for (pkg in names(refs)) {
  repo <- paste0("morinlab/", pkg, "@", refs[[pkg]])
  message("Installing ", repo)
  install_with_retry(repo)
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package failed to install: ", pkg, call. = FALSE)
  }
}
message("All six GAMBLR packages installed successfully.")
