#!/usr/bin/env Rscript
args0 <- commandArgs(FALSE); farg <- grep("^--file=", args0, value=TRUE)
PROJECT_DIR <- if(length(farg)) normalizePath(dirname(sub("^--file=", "", farg[[1L]])), winslash="/", mustWork=TRUE) else normalizePath(getwd(), winslash="/", mustWork=TRUE)
source(file.path(PROJECT_DIR,"R","load_all.R")); load_v56_modules(PROJECT_DIR)
res <- run_pipeline_until("validation", project_dir=PROJECT_DIR)
cat("v5.6 stage 'validation' completed in: ", res$ctx$root %||% "see output root", "\n", sep="")
