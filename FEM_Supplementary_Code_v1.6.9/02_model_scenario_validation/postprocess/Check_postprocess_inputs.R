# =============================================================================
# Check_postprocess_inputs.R
# Validate paired out-of-fold inputs before computing summaries and contrasts.
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L || !nzchar(args[1])) {
  stop("Usage: Rscript Check_postprocess_inputs.R RESULTS_DIR", call. = FALSE)
}
output_dir <- args[1]

cat("[check] Output/results directory:", output_dir, "\n")
if (!dir.exists(output_dir)) {
  stop("[check] Directory does not exist: ", output_dir)
}

rds_path <- file.path(output_dir, "all_results.rds")
cat("[check] all_results.rds:", rds_path, "\n")
if (!file.exists(rds_path)) {
  stop("[check] Missing all_results.rds. Pass the directory written by model_fit/run_nested_MS.R.")
}

all_results <- readRDS(rds_path)
# Validate before any paired bootstrap or publication output.
.fem_sources <- unlist(lapply(sys.frames(),function(e)e$ofile),use.names=FALSE)
.fem_script <- if(length(.fem_sources))tail(.fem_sources,1L) else sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1])
source(file.path(dirname(normalizePath(.fem_script)),"..","..","R","oof_contracts.R"))
rm(.fem_sources,.fem_script)
all_results <- fem_use_raw_rf_oof(all_results)
cat("[check] Prediction type: RF_raw (uncalibrated outer-test OOF).\n")

cat("[check] all_results entries:", length(all_results), "\n")

scenario_keep <- c(
  "S0_Env", "S1a_Prod", "S1b_Het", "S1c_Temp",
  "S2a_Prod_Het", "S2b_Prod_Temp", "S2c_Het_Temp",
  "S3_RS_Full", "S4_Integrated"
)

summary_rows <- list()
missing_coords <- 0L
missing_pred <- 0L

for (key in names(all_results)) {
  res <- all_results[[key]]
  if (is.null(res$model_results)) next
  target <- if (!is.null(res$target_var)) res$target_var else NA_character_
  pft <- if (!is.null(res$pft)) res$pft else NA_character_
  cvm <- if (!is.null(res$cv_method)) res$cv_method else NA_character_
  dt <- if (!is.null(res$data_type)) res$data_type else NA_character_
  scns <- intersect(scenario_keep, names(res$model_results))
  ok_scns <- 0L
  for (scn in scns) {
    sr <- res$model_results[[scn]]
    if (!isTRUE(sr$success)) next
    ok_scns <- ok_scns + 1L
    if (is.null(sr$predictions) || is.null(sr$observed)) missing_pred <- missing_pred + 1L
    if (is.null(sr$coords) && is.null(sr$coordinates)) missing_coords <- missing_coords + 1L
  }
  summary_rows[[length(summary_rows)+1L]] <- data.frame(
    Key = key, Data_Type = dt, Target = target, PFT = pft, CV_Method = cvm,
    Successful_Scenarios = ok_scns, stringsAsFactors = FALSE
  )
}

summary_df <- if (length(summary_rows)) do.call(rbind, summary_rows) else data.frame()
if (nrow(summary_df) > 0) {
  cat("[check] Result strata:\n")
  print(summary_df)
  cat("[check] CV methods found:", paste(sort(unique(summary_df$CV_Method)), collapse = ", "), "\n")
  cat("[check] Targets found:", paste(sort(unique(summary_df$Target)), collapse = ", "), "\n")
} else {
  warning("[check] No model_results found in all_results.rds.")
}

if (missing_pred > 0) {
  warning("[check] Some successful scenario results lack predictions/observed: ", missing_pred,
          ". Compute_fig23_stats.R may skip or fail for those cells.")
}
if (missing_coords > 0) {
  warning("[check] Some successful scenario results lack coordinates: ", missing_coords,
          ". observation-level bootstrap can still run, but spatial-block bootstrap needs coords.")
}

cat("[check] OK. Run 02_model_scenario_validation/run_postprocess_submission.sh with this results directory.\n")
