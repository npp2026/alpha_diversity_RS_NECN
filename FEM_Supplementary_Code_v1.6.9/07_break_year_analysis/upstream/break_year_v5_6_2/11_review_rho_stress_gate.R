#!/usr/bin/env Rscript
# Re-evaluate an already completed rho stress gate using the current audit
# decision policy only. This does not rerun raster/SupF/FDR analysis.
args0 <- commandArgs(FALSE); farg <- grep("^--file=", args0, value = TRUE)
PROJECT_DIR <- if (length(farg)) normalizePath(dirname(sub("^--file=", "", farg[[1L]])), winslash = "/", mustWork = TRUE) else normalizePath(getwd(), winslash = "/", mustWork = TRUE)
source(file.path(PROJECT_DIR, "R", "load_all.R")); load_v56_modules(PROJECT_DIR)
source(file.path(PROJECT_DIR, "validation", "rho_stress_gate.R"))
res <- review_rho_stress_gate_v56(PROJECT_DIR)
cat("v5.6 rho stress gate review completed for: ", res$gate_run, "\n", sep = "")
cat("Diagnostic gate: ", res$classification$overall_gate, "\n", sep = "")
cat("Sensitivity severity: ", res$classification$severity, "\n", sep = "")
cat("Report: ", res$report, "\n", sep = "")
