#!/usr/bin/env Rscript
args0 <- commandArgs(FALSE); farg <- grep("^--file=", args0, value = TRUE)
PROJECT_DIR <- if (length(farg)) normalizePath(dirname(sub("^--file=", "", farg[[1L]])), winslash = "/", mustWork = TRUE) else normalizePath(getwd(), winslash = "/", mustWork = TRUE)
source(file.path(PROJECT_DIR, "R", "load_all.R")); load_v56_modules(PROJECT_DIR)
source(file.path(PROJECT_DIR, "validation", "real_data_subset_audit.R"))
z <- run_real_data_subset_audit_v56(PROJECT_DIR)
cat("Real-data subset audit status: ", z$evaluation$status, "\n", sep = "")
cat("Audit report: ", z$evaluation$report, "\n", sep = "")
cat("Subset run output: ", z$run$ctx$root, "\n", sep = "")
