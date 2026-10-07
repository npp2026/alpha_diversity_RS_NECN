#!/usr/bin/env Rscript
# Report installed dependencies; do not install or modify the environment.
args <- commandArgs(trailingOnly = TRUE)
if (any(args %in% c("--help", "-h"))) {
  cat("Usage: Rscript scripts/check_environment.R\n",
      "Checks analytical dependencies.\n",
      "Python dependencies are listed in DEPENDENCIES.md.\n", sep = "")
  quit(save = "no", status = 0L)
}
if (length(args)) stop("Unknown option; see --help.")
if (getRversion() < "4.3.0") stop("Use R >= 4.3.0 for this package.")
packages <- c("ranger", "VSURF", "xgboost", "CAST", "terra", "sf",
              "dplyr", "tidyr", "yaml", "digest", "jsonlite")
versions <- vapply(packages, function(package) {
  if (requireNamespace(package, quietly = TRUE)) as.character(packageVersion(package)) else NA_character_
}, character(1))
print(data.frame(package = packages, version = versions), row.names = FALSE)
if (anyNA(versions)) stop("Missing packages: ", paste(packages[is.na(versions)], collapse = ", "))
optional <- c("trend", "modifiedmk", "mutoss", "FNN")
optional_versions <- vapply(optional, function(package) {
  if (requireNamespace(package, quietly = TRUE)) as.character(packageVersion(package)) else NA_character_
}, character(1))
print(data.frame(optional_package = optional, version = optional_versions), row.names = FALSE)
cat("Required R dependencies are available. Optional package-backed diagnostics may use project-local fallbacks. This is not an analysis execution.\n")
