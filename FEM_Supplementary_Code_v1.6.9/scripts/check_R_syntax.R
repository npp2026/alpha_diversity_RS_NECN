#!/usr/bin/env Rscript
# Parse only; do not source files or load their dependencies.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: Rscript check_R_syntax.R PACKAGE_ROOT", call. = FALSE)
root <- normalizePath(args[[1L]], mustWork = TRUE)
files <- list.files(root, pattern = "\\.[Rr]$", recursive = TRUE, full.names = TRUE)
failed <- character()
for (file in files) {
  err <- tryCatch({parse(file = file, keep.source = FALSE); NULL}, error = function(e) conditionMessage(e))
  if (!is.null(err)) {
    failed <- c(failed, file)
    cat("FAIL ", file, ": ", err, "\n", sep = "")
  }
}
cat("R parse: ", length(files) - length(failed), "/", length(files), " passed\n", sep = "")
if (length(failed)) quit(save = "no", status = 1L)
