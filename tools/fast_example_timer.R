# fast_example_timer.R
#
# Lighter-weight, per-example-timed replacement for tools/logExampleOutputs.R's
# devtools::run_examples() call. devtools::run_examples() has document = TRUE
# by default -- EVERY run calls devtools::document() (regenerating every .Rd
# file + NAMESPACE via roxygen2) and reloads the package via load_all() TWICE
# (once explicitly, once again via on.exit()), all before a single example
# actually runs. None of that is needed here: man/*.Rd is already committed
# and up to date, and this assumes the normal reinstall-then-test workflow
# (the package is already freshly installed when this runs).
#
# Usage (run from the package root, e.g. ~/GAMBLR.open):
#   Rscript tools/fast_example_timer.R [package_name] [--load-all]
#     package_name  defaults to GAMBLR.open
#     --load-all    use pkgload::load_all(".") instead of library(package_name)
#                   -- only needed to pick up uncommitted source changes
#                   without a full reinstall; slower than a plain library()
#                   call, so omit it when testing an already-installed build.
#
# Output:
#   GAMBLR_examples_output_fast.log  -- full example output (like the original)
#   GAMBLR_examples_timing.tsv       -- one row per example file, seconds + status,
#                                        slowest first

args <- commandArgs(trailingOnly = TRUE)
pkg <- if (length(args) >= 1 && !startsWith(args[1], "--")) args[1] else "GAMBLR.open"
use_load_all <- "--load-all" %in% args

log_file <- "GAMBLR_examples_output_fast.log"
timing_file <- "GAMBLR_examples_timing.tsv"

sink(log_file)
cat(paste("=== STARTED AT", Sys.time(), "===\n"))

t_load <- system.time({
  if (use_load_all) {
    suppressMessages(pkgload::load_all(".", export_all = FALSE, helpers = FALSE, quiet = TRUE))
  } else {
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }
})
cat(sprintf("[TIMING] package load (%s): %.2fs\n",
            if (use_load_all) "load_all" else "library", t_load[["elapsed"]]))

rd_files <- sort(list.files("man", pattern = "\\.Rd$", full.names = TRUE))
cat(sprintf("Found %d .Rd files\n", length(rd_files)))

timings <- data.frame(file = character(0), seconds = numeric(0), status = character(0),
                       stringsAsFactors = FALSE)

for (f in rd_files) {
  rd <- tryCatch(tools::parse_Rd(f), error = function(e) NULL)
  if (is.null(rd)) next
  tags <- vapply(rd, function(x) attr(x, "Rd_tag"), character(1))
  if (!"\\examples" %in% tags) next

  ex_file <- tempfile(fileext = ".R")
  tools::Rd2ex(rd, out = ex_file, commentDontrun = TRUE, commentDonttest = TRUE)
  if (!file.exists(ex_file) || !length(readLines(ex_file, warn = FALSE))) next

  fname <- basename(f)
  cat(sprintf("\n--- %s ---\n", fname))
  t0 <- Sys.time()
  status <- "OK"
  tryCatch(
    source(ex_file, echo = TRUE, local = new.env(parent = globalenv())),
    error = function(e) status <<- paste("ERROR:", conditionMessage(e))
  )
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  cat(sprintf("[TIMING] %-45s %8.2fs  %s\n", fname, elapsed, status))
  timings <- rbind(timings, data.frame(file = fname, seconds = elapsed, status = status,
                                        stringsAsFactors = FALSE))
}

cat(paste("=== COMPLETED AT", Sys.time(), "===\n"))
sink()

timings <- timings[order(-timings$seconds), ]
write.table(timings, timing_file, sep = "\t", quote = FALSE, row.names = FALSE)

cat(sprintf("\nTotal example time: %.1fs across %d files (package load: %.1fs)\n",
            sum(timings$seconds), nrow(timings), t_load[["elapsed"]]))
cat("Slowest 10:\n")
print(head(timings, 10), row.names = FALSE)
cat(sprintf("\nFull output: %s\nPer-file timings: %s\n", log_file, timing_file))
