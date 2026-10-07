#!/usr/bin/env Rscript
# Portable orchestration; analytical calculations remain in the supplied scripts.
args <- commandArgs(trailingOnly = TRUE)
help_text <- c(
  "FEM potential-model workflow",
  "Set FEM_POTENTIAL_DATA_DIR to the directory containing input data.",
  "Usage: Rscript run_potential_workflow.R [--check] [--steps=tune,predict,oldage]",
  "--check          Parse requested scripts and check inputs/packages; do not fit models.",
  "All requested step failures return nonzero. Fresh analysis output directories are required."
)
if (any(args %in% c("--help", "-h"))) {
  cat(paste(help_text, collapse = "\n"), "\n")
  quit(save = "no", status = 0L)
}
bad <- args[!grepl("^(--check|--steps=(tune|predict|oldage)(,(tune|predict|oldage))*)$", args)]
if (length(bad)) stop("Unknown/invalid option(s): ", paste(bad, collapse = ", "), call. = FALSE)
for (prefix in "--steps=") {
  if (sum(startsWith(args, prefix)) > 1L) stop("Repeated option: ", prefix, call. = FALSE)
}
option <- function(prefix, default) {
  found <- args[startsWith(args, prefix)]
  if (length(found)) substring(found, nchar(prefix) + 1L) else default
}
steps <- strsplit(option("--steps=", "predict,oldage"), ",", fixed = TRUE)[[1L]]
canonical <- c("tune", "predict", "oldage")
if (any(!steps %in% canonical) || anyDuplicated(steps) || !identical(steps, canonical[canonical %in% steps])) {
  stop("Steps must be unique and in tune,predict,oldage order.", call. = FALSE)
}

file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (!length(file_arg)) stop("Run this entry point with Rscript.", call. = FALSE)
script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]), winslash = "/", mustWork = TRUE))
data_dir <- Sys.getenv("FEM_POTENTIAL_DATA_DIR", unset = "")
if (!nzchar(data_dir)) stop("Set FEM_POTENTIAL_DATA_DIR explicitly. See README.md.", call. = FALSE)
data_dir <- normalizePath(path.expand(data_dir), winslash = "/", mustWork = TRUE)
if (!dir.exists(data_dir)) stop("FEM_POTENTIAL_DATA_DIR is not a directory.", call. = FALSE)
Sys.setenv(FEM_POTENTIAL_DATA_DIR = data_dir)
if (getRversion() < "4.3.0") stop("R >= 4.3.0 is required by this archive.", call. = FALSE)

script_names <- c(tune = "tune_qrf_params.R", predict = "predict_age_sensitivity.R", oldage = "oldage_validation_and_aoa_coverage.R")
paths <- setNames(file.path(script_dir, unname(script_names[steps])), steps)
for (path in paths) {
  if (!file.exists(path)) stop("Missing script: ", path, call. = FALSE)
  parse(file = path, keep.source = FALSE)
}
pkgs <- c("terra", "ranger", "CAST", "sf")
missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install required R packages: ", paste(missing, collapse = ", "), call. = FALSE)
# v1.6: SI S5.4 parameters are the default; a fresh tuning result is opt-in.
qrf_parameter_mode <- Sys.getenv("FEM_QRF_PARAMETER_MODE", "manuscript")
if (!qrf_parameter_mode %in% c("manuscript", "retuned")) stop("Invalid FEM_QRF_PARAMETER_MODE")
if ("predict" %in% steps && qrf_parameter_mode == "manuscript") {
  paths<-c(paths,static_export=file.path(script_dir,"export_static_q95.R"))
  invisible(parse(paths[["static_export"]],keep.source=FALSE))
}
output_names<-c(tune="qrf_parameter_tuning_outputs",predict="quantile_sensitivity_outputs",oldage="oldage_validation_outputs",static_export="static_q95_age100")
for(nm in intersect(names(paths),names(output_names))){
 dest<-file.path(data_dir,output_names[[nm]])
 if(dir.exists(dest)&&length(list.files(dest,all.files=TRUE,no..=TRUE)))stop("Use a fresh output directory: ",dest)
}
snippet <- file.path(data_dir, "qrf_parameter_tuning_outputs", "apply_best_qrf_params_to_final_mapping.R")
if (length(steps)) {
  inputs <- file.path(data_dir, c("train4pot.csv", "ENV.tif", "templ_1km.tif"))
  missing_inputs <- inputs[!file.exists(inputs)]
  if (length(missing_inputs)) stop("Missing input(s): ", paste(missing_inputs, collapse = ", "), call. = FALSE)
  header <- names(utils::read.csv(inputs[[1L]], nrows = 0L, check.names = FALSE))
  required_columns <- c("long", "lat", "Forest_age", "Rich_tree", "Shannon_wiener", "DEMc", "bio6_wc", "bio10_wc", "bio12_wc", "bio17_wc", "bio15_wc", "bio4_wc")
  missing_columns <- setdiff(required_columns, header)
  if (length(missing_columns)) stop("CSV missing required column(s): ", paste(missing_columns, collapse = ", "), call. = FALSE)
  env <- terra::rast(inputs[[2L]])
  expected_layers <- c("DEM", "BIO6", "BIO10", "BIO12", "BIO17", "BIO15", "BIO4")
  if (anyDuplicated(names(env))) stop("ENV.tif contains duplicate layer names.", call. = FALSE)
  if (length(setdiff(expected_layers, names(env)))) stop("ENV.tif missing named predictor layer(s): ", paste(setdiff(expected_layers, names(env)), collapse = ", "), call. = FALSE)
  template <- terra::rast(inputs[[3L]])
  if (!nzchar(terra::crs(env)) || !nzchar(terra::crs(template))) stop("Both rasters must have a CRS.", call. = FALSE)
  if (terra::nlyr(template) != 1L) stop("templ_1km.tif must have one mask layer.", call. = FALSE)
  unit_name <- tolower(sf::st_crs(terra::crs(template))$units_gdal)
  if (terra::is.lonlat(template) || length(unit_name) != 1L || is.na(unit_name) || !unit_name %in% c("metre", "meter", "metres", "meters", "m") || any(abs(terra::res(template) - 1000) > 0.01)) {
    stop("The publication template must use a projected metre CRS with 1000 x 1000 m resolution.", call. = FALSE)
  }
  if (qrf_parameter_mode == "retuned" && !"tune" %in% steps && any(c("predict", "oldage") %in% steps) && !file.exists(snippet)) {
    stop("Missing tuned parameter snippet. Run --steps=tune first, or restore the verified production snippet: ", snippet, call. = FALSE)
  }
  if ("oldage" %in% steps && !"predict" %in% steps) {
    aoa_summary <- file.path(data_dir, "quantile_sensitivity_outputs", "AOA_summary_all_responses_ages.csv")
    if (!file.exists(aoa_summary)) stop("Run prediction/AOA before the old-age reporting step: ", aoa_summary, call. = FALSE)
  }
}
cat("Preflight passed for: ", paste(names(paths), collapse = ", "), "\n", sep = "")
cat("Data/output directory: ", data_dir, "\n", sep = "")
if ("--check" %in% args) quit(save = "no", status = 0L)

rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
if (!file.exists(rscript)) rscript <- Sys.which("Rscript")
if (!nzchar(rscript)) stop("Rscript executable not found.", call. = FALSE)
run_dir <- file.path(data_dir, "logs", paste0("potential_", format(Sys.time(), "%Y%m%dT%H%M%S"), "_", Sys.getpid()))
dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
writeLines(capture.output(utils::sessionInfo()), file.path(run_dir, "sessionInfo.txt"))
utils::write.csv(data.frame(package = pkgs, version = vapply(pkgs, function(p) as.character(utils::packageVersion(p)), character(1))), file.path(run_dir, "package_versions.csv"), row.names = FALSE)
writeLines(c(paste0("data_directory=", data_dir), paste0("steps=", paste(steps, collapse = ",")), "oldage_folds=10; repeats=3; training_subset=age>=100"), file.path(run_dir, "requested_workflow.txt"))
status <- data.frame(step = character(), exit_code = integer(), log = character(), stringsAsFactors = FALSE)
for (step in names(paths)) {
  if (qrf_parameter_mode == "retuned" && step %in% c("predict", "oldage") && !file.exists(snippet)) stop("Required tuned parameter snippet was not produced: ", snippet, call. = FALSE)
  logfile <- file.path(run_dir, paste0(step, ".log"))
  cat("Running ", step, "; log: ", logfile, "\n", sep = "")
  code <- system2(rscript, args = c("--vanilla", shQuote(unname(paths[[step]]))), stdout = logfile, stderr = logfile)
  status <- rbind(status, data.frame(step = step, exit_code = as.integer(code), log = logfile))
  utils::write.csv(status, file.path(run_dir, "step_status.csv"), row.names = FALSE)
  if (!identical(as.integer(code), 0L)) {
    cat("FAILED: ", step, ". Inspect ", logfile, "\n", sep = "", file = stderr())
    quit(save = "no", status = 1L)
  }
}
cat("Requested steps finished. Logs: ", run_dir, "\n", sep = "")
cat("Compare outputs and saved parameters with the manuscript before claiming reproduction.\n")

if ("predict" %in% steps && qrf_parameter_mode == "retuned") message("Retuned run completed; automatic MS static export is intentionally not produced.")
