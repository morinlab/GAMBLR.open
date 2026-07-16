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

refs <- list(
  GAMBLR.data    = Sys.getenv("GAMBLR_DATA_REF", "master"),
  GAMBLR.helpers = Sys.getenv("GAMBLR_HELPERS_REF", "master"),
  GAMBLR.utils   = Sys.getenv("GAMBLR_UTILS_REF", "master"),
  GAMBLR.viz     = Sys.getenv("GAMBLR_VIZ_REF", "master"),
  GAMBLR.predict = Sys.getenv("GAMBLR_PREDICT_REF", "master"),
  GAMBLR.open    = Sys.getenv("GAMBLR_OPEN_REF", "master")
)

repos <- paste0("morinlab/", names(refs), "@", unlist(refs))
message("Installing: ", paste(repos, collapse = ", "))

remotes::install_github(repos, dependencies = FALSE, upgrade = "never")

# Fail the build loudly if anything didn't actually install, rather than
# only discovering it the first time someone tries to use the image.
for (pkg in names(refs)) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package failed to install: ", pkg, call. = FALSE)
  }
}
message("All six GAMBLR packages installed successfully.")
