regional_annual_means <- function(cfg, inventory, ctx) {
  assert_packages("terra")
  regions <- load_regions(cfg)
  if (is.null(regions)) stop("Regional aggregation requires data.regions_shapefile.", call. = FALSE)
  yrs <- as.numeric(unlist(cfg$data$years)); rows <- list()
  for (response in names(inventory$files)) {
    st <- load_response_stack(inventory$files[[response]], cfg, mask_to_regions = FALSE)
    reg <- regions; if (!terra::same.crs(reg, st)) reg <- terra::project(reg, terra::crs(st))
    meta <- region_metadata(reg, cfg, require_unique_id = TRUE)
    reg$.__zone_id__ <- seq_len(nrow(reg))
    zones <- terra::rasterize(reg, st[[1]], field = ".__zone_id__", touches = FALSE)
    area <- terra::cellSize(st[[1]], unit = "km")
    for (j in seq_along(yrs)) {
      x <- st[[j]]
      num <- terra::zonal(x * area, zones, fun = "sum", na.rm = TRUE)
      den <- terra::zonal(terra::ifel(is.na(x), NA, area), zones, fun = "sum", na.rm = TRUE)
      names(num) <- c("zone_id", "weighted_sum"); names(den) <- c("zone_id", "valid_area_km2")
      z <- merge(num, den, by = "zone_id", all = TRUE)
      z$value <- ifelse(is.finite(z$valid_area_km2) & z$valid_area_km2 > 0, z$weighted_sum / z$valid_area_km2, NA_real_)
      idx <- match(as.integer(z$zone_id), meta$feature)
      rows[[length(rows) + 1L]] <- data.frame(region = meta$region[idx], region_label = meta$region_label[idx],
        region_group = meta$region_group[idx], region_group_label = meta$region_group_label[idx],
        response = response, year = yrs[[j]], value = z$value, valid_area_km2 = z$valid_area_km2,
        stringsAsFactors = FALSE)
    }
  }
  tab <- do.call(rbind, rows); tab <- tab[order(tab$response, tab$region, tab$year), ]
  write_csv(tab, file.path(ctx$root, "uncertainty", "trajectory", "regional_annual_means.csv")); tab
}

run_regional_trajectory_analysis <- function(cfg, inventory, ctx) {
  if (!isTRUE(cfg$bootstrap$trajectory$enabled)) return(NULL)
  root <- ensure_dir(file.path(ctx$root, "uncertainty", "trajectory"))
  dat <- regional_annual_means(cfg, inventory, ctx); yrs <- as.numeric(unlist(cfg$data$years)); out <- list(); fitted_rows <- list()
  # Short-series rho correction is built once for linear-null and hinge residuals.
  rcfg <- cfg$rho_estimation; grid <- seq(as.numeric(rcfg$bias_correction$grid_min), as.numeric(rcfg$bias_correction$grid_max), by=as.numeric(rcfg$bias_correction$grid_step))
  rb <- as.integer(cfg$bootstrap$trajectory$rho_bias_sim_B %||% min(5000L, as.integer(rcfg$bias_correction$simulation_B)))
  hms <- as.integer(rcfg$hinge_min_segment %||% 3L)
  map_lin <- if (isTRUE(rcfg$bias_correction$enabled)) rho_bias_mapping_cached(cfg, yrs, grid, rb, "linear", seed_from_key(cfg$rng$master_seed,"trajectory_rho_bias","linear"), hms) else NULL
  map_hinge <- if (isTRUE(rcfg$bias_correction$enabled)) rho_bias_mapping_cached(cfg, yrs, grid, rb, "hinge", seed_from_key(cfg$rng$master_seed,"trajectory_rho_bias","hinge"), hms) else NULL
  groups <- split(dat, interaction(dat$region, dat$response, drop = TRUE))
  for (g in groups) {
    g <- g[order(g$year), ]; y <- g$value; response <- as.character(g$response[[1L]])
    if (length(y) != length(yrs) || any(!is.finite(y))) next
    slope_tol <- as.numeric(cfg$analysis$slope_tol[[response]] %||% 0); m <- as.integer(cfg$analysis$primary$min_segment)
    fit <- fit_break_core(y, yrs, m, m, slope_tol, cfg$analysis$tie_rule); if (!isTRUE(fit$testable)) next
    lin <- fit_ols_sse(cbind(1, yrs - yrs[[1L]]), y)
    rho_null_raw <- estimate_lag1_single(lin$resid, rcfg$max_abs_rho); rho_boot_raw <- estimate_lag1_single(fit$resid, rcfg$max_abs_rho)
    rho_null <- if (!is.null(map_lin)) bias_correct_rho(rho_null_raw, map_lin, rcfg$max_abs_rho) else rho_null_raw
    rho_boot <- if (!is.null(map_hinge)) bias_correct_rho(rho_boot_raw, map_hinge, rcfg$max_abs_rho) else rho_boot_raw
    key <- paste(g$region[[1L]], response, sep="|"); seed <- seed_from_key(cfg$rng$master_seed, "regional_supf", key)
    null <- simulate_supf_null(yrs, m, m, rho_null, as.integer(cfg$bootstrap$trajectory$supf_mc_B), seed,
                               min(10000L, as.integer(cfg$bootstrap$trajectory$supf_mc_B)), cfg$rng$kind)
    p <- mc_pvalue(fit$supF, null)
    boot <- bootstrap_regional_trajectory(y, yrs, fit, rho_boot, null, cfg, key, response)
    area_frac <- if (max(g$valid_area_km2, na.rm=TRUE) > 0) min(g$valid_area_km2, na.rm=TRUE) / max(g$valid_area_km2, na.rm=TRUE) else NA_real_
    row <- cbind(data.frame(region=g$region[[1L]], region_label=g$region_label[[1L]],
      region_group=g$region_group[[1L]], region_group_label=g$region_group_label[[1L]],
      response=response, break_year=fit$break_year, slope_before=fit$slope_before,
      slope_after=fit$slope_after, supF=fit$supF, p_mc=p, pointwise_alpha=as.numeric(cfg$analysis$primary$alpha), direction_ok=fit$direction_ok,
      valid_area_min_fraction=area_frac, coverage_warning=is.finite(area_frac) && area_frac < as.numeric(cfg$bootstrap$trajectory$min_valid_area_fraction %||% 0.95),
      rho_null_raw=rho_null_raw, rho_null_corrected=rho_null, rho_bootstrap_raw=rho_boot_raw, rho_bootstrap_corrected=rho_boot,
      inference_scope="secondary regional trajectory; pointwise Monte Carlo p-value, not raster FDR", stringsAsFactors=FALSE), boot)
    out[[length(out)+1L]] <- row
    fitted_rows[[length(fitted_rows)+1L]] <- data.frame(region=g$region[[1L]], region_label=g$region_label[[1L]],
      region_group=g$region_group[[1L]], region_group_label=g$region_group_label[[1L]],
      response=response, year=yrs, observed=y,
      fitted=fit$fitted, residual=fit$resid, break_year=fit$break_year,
      bootstrap_break_q025=boot$bootstrap_break_q025, bootstrap_break_median=boot$bootstrap_break_median,
      bootstrap_break_q975=boot$bootstrap_break_q975, stringsAsFactors=FALSE)
  }
  tab <- if(length(out)) do.call(rbind,out) else data.frame(); fit_tab <- if(length(fitted_rows)) do.call(rbind,fitted_rows) else data.frame()
  write_csv(tab,file.path(root,"regional_trajectory_breaks.csv")); write_csv(fit_tab,file.path(root,"regional_trajectory_fitted.csv"))
  list(summary=tab, fitted=fit_tab)
}
