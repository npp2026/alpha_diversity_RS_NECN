#!/usr/bin/env Rscript
# Quick checks only. Write new evidence outside this source package.
args <- commandArgs(trailingOnly = TRUE)
if (any(args %in% c("--help", "-h"))) {
  cat("Usage: Rscript scripts/run_MS_validation.R [PACKAGE_ROOT] [OUTPUT_DIR]\n",
      "OUTPUT_DIR must be empty and outside PACKAGE_ROOT.\n",
      "The package self-check is listed in scripts/quick_check_suites.csv.\n", sep = "")
  quit(save = "no", status = 0L)
}
if (length(args) > 2L) stop("Expected at most PACKAGE_ROOT and OUTPUT_DIR.")
root <- normalizePath(if (length(args)) args[[1L]] else getwd(),
                      winslash = "/", mustWork = TRUE)
if (!file.exists(file.path(root, "VERSION"))) stop("Invalid package root: ", root)
output <- if (length(args) == 2L) args[[2L]] else file.path(
  dirname(root), paste0("FEM_validation_", format(Sys.time(), "%Y%m%dT%H%M%S"), "_", Sys.getpid())
)
output <- if (dir.exists(output)) {
  normalizePath(output, winslash = "/", mustWork = TRUE)
} else {
  file.path(normalizePath(dirname(output), winslash = "/", mustWork = TRUE), basename(output))
}
if (identical(output, root) || startsWith(output, paste0(root, "/"))) {
  stop("Choose an output directory outside the source package.")
}
if (length(list.files(output, all.files = TRUE, no.. = TRUE))) {
  stop("Use a fresh output directory: ", output)
}
if (!dir.exists(output) && !dir.create(output)) stop("Cannot create output directory: ", output)
Sys.setenv(FEM_PACKAGE_ROOT = root)
inventory <- read.csv(file.path(root, "scripts", "quick_check_suites.csv"), stringsAsFactors = FALSE)
if (!identical(names(inventory), c("suite", "script")) || !nrow(inventory) ||
    anyNA(inventory) || anyDuplicated(inventory$suite) ||
    any(!grepl("^[A-Za-z0-9_]+$", inventory$suite))) stop("Invalid quick-check inventory")
suite_paths <- normalizePath(file.path(root, inventory$script), winslash = "/", mustWork = TRUE)
if (any(!startsWith(suite_paths, paste0(root, "/"))) || any(dir.exists(suite_paths))) {
  stop("Quick-check scripts must be files inside the package")
}
files <- list.files(root, pattern = "\\.[Rr]$", recursive = TRUE, full.names = TRUE)
syntax <- do.call(rbind, lapply(files, function(file) {
  error <- tryCatch({parse(file, keep.source = FALSE); ""},
                    error = function(e) conditionMessage(e))
  data.frame(file = substring(file, nchar(root) + 2L),
             status = if (nzchar(error)) "FAIL" else "PASS", detail = error)
}))
write.csv(syntax, file.path(output, "R_syntax.csv"), row.names = FALSE)
if (any(syntax$status == "FAIL")) stop("R syntax failures; see ", output)

rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
status <- data.frame(suite = character(), status = character(), exit_code = integer(), log = character())
for (i in seq_len(nrow(inventory))) {
  logfile <- file.path(output, paste0(inventory$suite[i], ".log"))
  code <- system2(rscript, c("--vanilla", shQuote(suite_paths[i])), stdout = logfile, stderr = logfile)
  status <- rbind(status, data.frame(suite = inventory$suite[i],
    status = if (code == 0L) "PASS" else "FAIL", exit_code = as.integer(code), log = logfile))
  write.csv(status, file.path(output, "quick_check_results.csv"), row.names = FALSE)
  if (code != 0L) stop("Quick-check suite failed: ", inventory$suite[i], "; see ", logfile)
}
capture.output(sessionInfo(), file = file.path(output, "R_sessionInfo.txt"))
packages <- c("ranger", "VSURF", "xgboost", "terra", "sf", "yaml",
              "digest", "jsonlite", "dplyr", "tidyr", "CAST")
versions <- vapply(packages, function(package) {
  if (requireNamespace(package, quietly = TRUE)) as.character(packageVersion(package)) else NA_character_
}, character(1))
write.csv(data.frame(package = packages, version = versions),
          file.path(output, "R_packages.csv"), row.names = FALSE)
writeLines(c(paste("Version:", readLines(file.path(root, "VERSION"))),
             "PASS: R syntax and the suites in scripts/quick_check_suites.csv.",
             "The supplied self-check is a reduced computational check, not full integration validation.",
             "Integration and original-data reproduction are separate checks."),
           file.path(output, "RUN_COMPLETE.txt"))
message("Quick checks passed. Evidence: ", output)
