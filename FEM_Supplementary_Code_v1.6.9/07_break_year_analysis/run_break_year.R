#!/usr/bin/env Rscript
# Portable adapter for the supplied figure-free derivative; see SOURCE_PROVENANCE.csv.
main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  if ("--help" %in% args) {
    cat(paste(
      "Usage: Rscript 07_break_year_analysis/run_break_year.R [--check] [--synthetic] [--diagnostic]",
      "       [--stage=primary|sensitivity|uncertainty|spatial|tables|validation]",
      "       [--config=/absolute/path/config.yml]",
      "Default: manuscript-method primary and sensitivity stages.",
      "Real inputs: set FEM_BREAK_DATA_DIR and FEM_BREAK_OUTPUT_DIR.",
      "Synthetic mode: uses the bundled generated fixture, not manuscript data.",
      "--diagnostic allows non-publication calibration settings and labels the run.",
      "--check verifies dependencies, annual files, geometry and region fields; no inference.",
      "Stages are cumulative; each run starts from primary inference.",
      "This package runs analysis and diagnostics; historical publication workflows are archived.",
      sep = "\n"), "\n")
    return(invisible(NULL))
  }
  valid <- args %in% c("--check", "--synthetic", "--diagnostic") | grepl("^--(stage|config)=.+$", args)
  if (any(!valid)) stop("Unknown/empty argument: ", paste(args[!valid], collapse = ", "), call. = FALSE)
  for (key in c("stage", "config")) {
    if (sum(grepl(paste0("^--", key, "="), args)) > 1L) stop("Duplicate option: ", key, call. = FALSE)
  }
  opt <- function(key, default) {
    z <- grep(paste0("^--", key, "="), args, value = TRUE)
    if (length(z)) sub(paste0("^--", key, "="), "", z[[1L]]) else default
  }
  stage <- opt("stage", "sensitivity")
  if (!stage %in% c("primary", "sensitivity", "uncertainty", "spatial", "tables", "validation")) stop("Invalid stage", call. = FALSE)
  synthetic <- "--synthetic" %in% args
  diagnostic <- synthetic || "--diagnostic" %in% args
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (!length(file_arg)) stop("Run this adapter with Rscript, not source().", call. = FALSE)
  module <- dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]), winslash = "/", mustWork = TRUE))
  engine <- normalizePath(file.path(module, "upstream", "break_year_v5_6_2"), winslash = "/", mustWork = TRUE)
  if (synthetic && any(grepl("^--config=", args))) stop("--synthetic and --config cannot be combined.", call. = FALSE)
  config <- if (synthetic) file.path(engine, "config", "v5_6_synthetic.yml") else opt("config", file.path(module, "config", "fem_break_year.yml"))
  config <- normalizePath(path.expand(config), winslash = "/", mustWork = TRUE)
  data_dir <- if (synthetic) engine else Sys.getenv("FEM_BREAK_DATA_DIR", unset = "")
  if (!nzchar(data_dir)) stop("Set FEM_BREAK_DATA_DIR to the real annual-raster input directory.", call. = FALSE)
  data_dir <- normalizePath(path.expand(data_dir), winslash = "/", mustWork = TRUE)
  output <- Sys.getenv("FEM_BREAK_OUTPUT_DIR", unset = "")
  if (!nzchar(output)) stop("Set FEM_BREAK_OUTPUT_DIR to a separate output directory.", call. = FALSE)
  output <- normalizePath(path.expand(output), winslash = "/", mustWork = FALSE)
  package <- dirname(module)
  if (identical(output, data_dir) || startsWith(output, paste0(package, "/")) || identical(output, package)) {
    stop("Choose an output directory separate from the input root and outside the source package.", call. = FALSE)
  }
  # Explicit FEM settings take precedence over unrelated ambient data/output settings.
  Sys.setenv(DATA_ROOT = data_dir, OUTPUT_ROOT = output, BREAK_CONFIG = config)
  source(file.path(engine, "R", "load_all.R"), local = .GlobalEnv)
  load_v56_modules(engine)
  assert_packages(c("terra", "yaml", "digest", "jsonlite"))
  cfg <- read_config(config, engine)
  source(file.path(module, "method_profile.R"), local = .GlobalEnv)
  validate_fem_method_profile(cfg, diagnostic)
  years <- as.integer(unlist(cfg$data$years))
  template <- NULL
  for (response in names(cfg$data$responses)) {
    files <- discover_response_files(cfg$data$responses[[response]], years)
    for (f in files$file) {
      raster <- terra::rast(f)
      if (terra::nlyr(raster) != 1L) stop("Expected one realized-diversity layer: ", f, call. = FALSE)
      if (!nzchar(terra::crs(raster))) stop("Missing CRS: ", f, call. = FALSE)
      if (terra::is.lonlat(raster) || any(abs(terra::res(raster) - 1000) > 1e-6)) {
        stop("Expected the documented projected 1-km input grid: ", f, call. = FALSE)
      }
      if (is.null(template)) template <- raster else terra::compareGeom(template, raster, stopOnError = TRUE)
    }
    cat(response, ": ", nrow(files), " aligned annual files\n", sep = "")
  }
  regions <- load_regions(cfg, template)
  if (is.null(regions)) stop("Region vector is required for this FEM integration.", call. = FALSE)
  meta <- region_metadata(regions, cfg)
  cat("Region features: ", nrow(meta), "; verify grouping against manuscript regions I-V.\n", sep = "")
  cat("Preflight passed. CRS units, OOD-mask provenance and numerical values still require review.\n")
  cat("SupF MC B: ", supf_mc_B_description(cfg), "; mode: ", if (diagnostic) "DIAGNOSTIC ONLY" else "manuscript-method replication", "\n", sep = "")
  if ("--check" %in% args) return(invisible(NULL))
  result <- run_pipeline_until(stage, config_path = config, project_dir = engine)
  capture.output(sessionInfo(), file = file.path(result$ctx$root, "logs", "FEM_sessionInfo.txt"))
  extras <- c(file.path(module, "run_break_year.R"), config)
  write_csv(file_manifest(extras, role = c("FEM_adapter", "FEM_config")), file.path(result$ctx$root, "config", "FEM_adapter_manifest.csv"))
  jsonlite::write_json(list(software_version = trimws(readLines(file.path(package, "VERSION"), n = 1L, warn = FALSE)), stage = stage,
    diagnostic_only = diagnostic, publication_freeze_certified = FALSE,
    method = "unconstrained_SSE_supF_AR1_MC_response_BH_then_direction",
    mc_B_by_response = as.list(supf_mc_B_by_response(cfg))),
    file.path(result$ctx$root, "FEM_RUN_METADATA.json"), auto_unbox = TRUE, pretty = TRUE)
  cat("Completed: ", result$ctx$root, "\n", sep = "")
}
tryCatch(main(), error = function(e) {
  message("Break-year workflow failed: ", conditionMessage(e))
  quit(save = "no", status = 1L)
})
