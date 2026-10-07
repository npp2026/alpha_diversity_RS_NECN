#!/usr/bin/env Rscript
args0 <- commandArgs(FALSE); farg <- grep("^--file=", args0, value = TRUE)
PROJECT_DIR <- if (length(farg)) normalizePath(dirname(sub("^--file=", "", farg[[1L]])), winslash = "/", mustWork = TRUE) else normalizePath(getwd(), winslash = "/", mustWork = TRUE)
source(file.path(PROJECT_DIR, "R", "load_all.R")); load_v56_modules(PROJECT_DIR)
source(file.path(PROJECT_DIR, "validation", "synthetic_e2e.R"))
args <- commandArgs(trailingOnly = TRUE); reuse <- "--reuse-data" %in% args
z <- run_synthetic_e2e_v56(PROJECT_DIR, reuse_data = reuse, fail_on_hard = TRUE)
cat("Synthetic v5.6 E2E status: ", z$evaluation$summary$status, "\n", sep = "")
cat("Synthetic run output: ", z$run$ctx$root, "\n", sep = "")
