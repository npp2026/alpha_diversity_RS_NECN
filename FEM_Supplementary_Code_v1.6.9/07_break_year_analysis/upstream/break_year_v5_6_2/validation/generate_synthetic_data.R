#!/usr/bin/env Rscript
args0 <- commandArgs(FALSE); farg <- grep("^--file=", args0, value = TRUE)
PROJECT_DIR <- if (length(farg)) normalizePath(dirname(dirname(sub("^--file=", "", farg[[1L]]))), winslash = "/", mustWork = TRUE) else normalizePath(getwd(), winslash = "/", mustWork = TRUE)
source(file.path(PROJECT_DIR, "R", "load_all.R")); load_v56_modules(PROJECT_DIR)
args <- commandArgs(trailingOnly = TRUE); reuse <- "--reuse" %in% args
root <- file.path(PROJECT_DIR, "validation", "synthetic_data")
if (reuse && file.exists(file.path(root, "synthetic_manifest.csv"))) {
  cat("Reusing existing synthetic dataset in: ", root, "\n", sep = "")
} else {
  z <- generate_synthetic_dataset(PROJECT_DIR, overwrite = TRUE)
  cat("Synthetic dataset generated in: ", z$root, "\n", sep = "")
}
