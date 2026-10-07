load_v56_modules <- function(project_dir = getwd()) {
  files <- c("utils.R", "parallel_runtime.R", "core_break.R", "ar1_calibration.R", "fdr_discrete.R", "io.R", "rho_estimation.R",
             "provenance.R", "synthetic_data.R", "raster_pipeline.R", "sensitivity.R", "bootstrap_pixel.R", "bootstrap_trajectory.R",
             "aggregation.R", "spatial_validation.R", "tables.R", "qa_report.R", "pipeline.R")
  for (f in files) source(file.path(project_dir, "R", f), local = .GlobalEnv)
  invisible(TRUE)
}
