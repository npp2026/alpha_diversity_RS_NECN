# The adapter protects the manuscript profile; engine derivation is recorded in SOURCE_PROVENANCE.csv.
# Reduced calibration is permitted only with an explicit diagnostic label.
validate_fem_method_profile <- function(cfg, diagnostic = FALSE) {
  if (diagnostic) {
    warning("DIAGNOSTIC ONLY: settings differ from the manuscript method profile.", call. = FALSE)
    return(invisible(TRUE))
  }
  same <- function(actual, expected, field) {
    if (!isTRUE(all.equal(unlist(actual, use.names = FALSE), expected, check.attributes = FALSE))) {
      stop("Manuscript method profile mismatch: ", field,
           ". Restore publication settings or explicitly use --diagnostic for tests.", call. = FALSE)
    }
  }
  same(cfg$data$years, 2001:2020, "years")
  same(names(cfg$data$responses), c("SR", "Shannon"), "response names/order")
  same(cfg$data$require_complete_series, TRUE, "complete annual series")
  same(cfg$data$mask_to_regions, TRUE, "region mask")
  same(cfg$analysis$primary$min_segment, 5, "primary segment length")
  same(cfg$analysis$primary$rho_source, "empirical_corrected", "rho source")
  same(cfg$analysis$primary$direction, "decline_recovery", "direction filter")
  same(cfg$analysis$slope_tol[c("SR", "Shannon")], c(0, 0), "slope tolerance")
  same(cfg$analysis$tie_rule, "earliest", "SSE tie rule")
  same(cfg$multiple_testing$family, "per_response", "FDR family")
  same(cfg$analysis$primary$fdr, "BH", "primary FDR")
  same(cfg$analysis$primary$alpha, 0.05, "alpha")
  same(supf_mc_B_by_response(cfg), c(200000, 1000000), "SR/Shannon Monte Carlo B")
  same(cfg$supf$p_rule, "greater_equal", "Monte Carlo tie rule")
  r <- cfg$rho_estimation
  same(c(r$enabled, r$separate_by_response, r$bias_correction$enabled), rep(TRUE, 3), "rho estimation")
  same(r$residual_model, "hinge", "rho residual model")
  same(r$hinge_min_segment, 3, "rho hinge segment length")
  same(c(r$sample_n, r$oversample_factor, r$max_abs_rho), c(10000, 4, 0.8), "rho sample/clipping")
  b <- r$bias_correction
  same(c(b$simulation_B, b$grid_min, b$grid_max, b$grid_step), c(10000, -0.8, 0.8, 0.02), "rho bias calibration")
  same(r$diagnostic_models, c("linear", "hinge"), "rho diagnostics")
  same(c(cfg$sensitivity$boundary$enabled, cfg$sensitivity$rho$enabled,
         cfg$sensitivity$multiplicity$enabled), rep(TRUE, 3), "sensitivity stages")
  same(cfg$sensitivity$boundary$min_segments, c(5, 4, 3, 2), "boundary settings")
  same(as.character(unlist(cfg$sensitivity$rho$values)),
       c("0", "0.2", "empirical", "linear_diagnostic", "0.4", "0.6"), "rho stress settings")
  same(cfg$sensitivity$multiplicity$methods, c("BH", "BY"), "multiplicity sensitivity")
  same(cfg$rng$kind, "L'Ecuyer-CMRG", "RNG kind")
  same(cfg$rng$master_seed, 12345, "master seed")
  invisible(TRUE)
}
