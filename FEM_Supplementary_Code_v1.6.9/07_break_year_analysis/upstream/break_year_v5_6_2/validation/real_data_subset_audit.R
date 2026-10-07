# Real-data subset audit for v5.6.
#
# This module performs two distinct checks:
#   1) full-input preflight on the actual rasters/vector metadata and coverage;
#   2) a spatially aligned cropped subset run through the complete v5.6 pipeline.
#
# Subset FDR/trajectory results are AUDIT diagnostics only. They are not a
# substitute for the full-study inferential family or full-region trajectories.

real_audit_metric <- function(category, check, status = c("PASS", "WARN", "FAIL", "INFO", "NOT_TESTED"),
                              value = NA, threshold = NA, note = "") {
  status <- match.arg(status)
  data.frame(category = category, check = check, status = status,
             value = as.character(value), threshold = as.character(threshold),
             note = as.character(note), stringsAsFactors = FALSE)
}

real_audit_read_spec <- function(project_dir) {
  assert_packages("yaml")
  path <- file.path(project_dir, "validation", "real_data_subset_spec.yml")
  if (!file.exists(path)) stop("Real subset audit spec not found: ", path, call. = FALSE)
  yaml::read_yaml(path)
}

real_audit_mode <- function(spec) {
  z <- tolower(Sys.getenv("REAL_SUBSET_MODE", unset = as.character(spec$mode$default %||% "standard")))
  if (!z %in% c("quick", "standard", "full")) stop("REAL_SUBSET_MODE must be quick, standard, or full.", call. = FALSE)
  z
}

real_audit_side_cells <- function(spec) {
  x <- suppressWarnings(as.integer(Sys.getenv("REAL_SUBSET_SIDE", unset = as.character(spec$subset$side_cells %||% 256L))))
  if (!is.finite(x) || x < 32L) stop("REAL_SUBSET_SIDE must be >= 32 cells.", call. = FALSE)
  x
}

real_audit_geom_signature <- function(r) {
  e <- terra::ext(r); rr <- terra::res(r)
  data.frame(nrow = terra::nrow(r), ncol = terra::ncol(r), ncell = terra::ncell(r),
             xmin = e$xmin, xmax = e$xmax, ymin = e$ymin, ymax = e$ymax,
             res_x = rr[[1L]], res_y = rr[[2L]], crs = terra::crs(r),
             lonlat = terra::is.lonlat(r), stringsAsFactors = FALSE)
}

real_audit_compare_response_geometry <- function(files_table, response) {
  assert_packages("terra")
  ref <- terra::rast(files_table$file[[1L]])
  ok <- vapply(files_table$file, function(f) {
    z <- terra::rast(f)
    isTRUE(terra::compareGeom(ref, z, stopOnError = FALSE))
  }, logical(1L))
  data.frame(response = response, year = files_table$year, file = files_table$file,
             geometry_matches_first_year = ok, stringsAsFactors = FALSE)
}

real_audit_layer_stats <- function(files_table, response, cfg) {
  assert_packages("terra")
  st <- load_response_stack(files_table, cfg)
  g <- terra::global(st, c("min", "max", "mean"), na.rm = TRUE)
  valid <- terra::global(!is.na(st), "sum", na.rm = TRUE)
  ncell <- terra::ncell(st)
  data.frame(response = response, year = files_table$year,
             min = as.numeric(g[, "min"]), max = as.numeric(g[, "max"]), mean = as.numeric(g[, "mean"]),
             valid_cells = as.numeric(valid[[1L]]), total_raster_cells = ncell,
             valid_fraction_raster = as.numeric(valid[[1L]]) / ncell,
             stringsAsFactors = FALSE)
}

real_audit_region_diagnostics <- function(cfg, template, thresholds) {
  assert_packages("terra")
  shp <- cfg$data$regions_shapefile
  if (!nzchar(shp) || !file.exists(shp)) {
    return(list(vector = NULL, table = data.frame(), candidates = data.frame(), metrics = list(
      real_audit_metric("regions", "regions_shapefile_present", if (isTRUE(cfg$data$mask_to_regions)) "FAIL" else "WARN", 0, 1,
                        "Regional/spatial interpretation is unavailable without the configured vector."))))
  }
  v <- terra::vect(shp)
  if (!terra::same.crs(v, template)) v <- terra::project(v, terra::crs(template))
  attrs <- as.data.frame(v)
  fields <- region_field_config(cfg)
  id_present <- nzchar(fields$id) && fields$id %in% names(attrs)
  ids <- if (id_present) as.character(attrs[[fields$id]]) else character()
  ids_complete <- id_present && !any(is.na(ids) | !nzchar(trimws(ids)))
  ids_unique <- ids_complete && !anyDuplicated(ids)
  candidates <- region_unique_character_candidates(v)
  candidate_names <- candidates$field[candidates$unique_character_candidate]

  optional_field_metric <- function(check, field, role) {
    if (!nzchar(field)) return(real_audit_metric("regions", check, "INFO", "not_configured", NA, paste(role, "is optional.")))
    ok <- field %in% names(attrs)
    real_audit_metric("regions", check, if (ok) "PASS" else "WARN", as.integer(ok), 1,
                      paste0(role, ": ", field, if (!ok) "; downstream output falls back where possible." else ""))
  }

  metrics <- list(
    real_audit_metric("regions", "regions_shapefile_present", "PASS", 1, 1, shp),
    real_audit_metric("regions", "region_available_fields", "INFO", paste(names(attrs), collapse = ","), NA,
                      "Fields available in the configured region vector."),
    real_audit_metric("regions", "region_unique_character_candidates", "INFO",
                      if (length(candidate_names)) paste(candidate_names, collapse = ",") else "none", NA,
                      "Candidate machine IDs: character/factor, complete, nonblank, and unique per polygon feature."),
    real_audit_metric("regions", "region_id_field_present", if (id_present) "PASS" else "FAIL",
                      as.integer(id_present), 1, fields$id),
    real_audit_metric("regions", "region_ids_complete",
                      if (!id_present) "NOT_TESTED" else if (ids_complete) "PASS" else "FAIL",
                      if (!id_present) NA else sum(!(is.na(ids) | !nzchar(trimws(ids)))), if (!id_present) NA else length(ids),
                      if (!id_present) "Cannot test completeness until region_id_field exists." else "Region IDs must be nonmissing and nonblank."),
    real_audit_metric("regions", "region_ids_unique",
                      if (!id_present) "NOT_TESTED" else if (!ids_complete) "NOT_TESTED" else if (ids_unique) "PASS" else "FAIL",
                      if (!id_present || !ids_complete) NA else length(unique(ids)), if (!id_present || !ids_complete) NA else length(ids),
                      if (!id_present) "Cannot test uniqueness until region_id_field exists." else if (!ids_complete) "Cannot test uniqueness until region IDs are complete." else "Duplicate region IDs can break/merge regional trajectory groups."),
    optional_field_metric("region_label_field_present", fields$label, "Publication/display label field"),
    optional_field_metric("region_group_field_present", fields$group, "Broader regional grouping field"),
    optional_field_metric("region_group_label_field_present", fields$group_label, "Broader grouping display label field"),
    real_audit_metric("regions", "projected_crs_for_spatial_blocks",
                      { wkt <- terra::crs(template); metric_ok <- !terra::is.lonlat(template) && grepl('CS[Cartesian', wkt, fixed = TRUE) && grepl('LENGTHUNIT["metre",1', wkt, fixed = TRUE);
                        if (!isTRUE(cfg$spatial$block_jackknife$enabled) || metric_ok) "PASS" else "FAIL" },
                      if (terra::is.lonlat(template)) "lonlat" else "projected", "projected Cartesian metre CRS",
                      "50-km block jackknife requires a projected CRS whose Cartesian axis unit is metre.")
  )
  area_parts <- tryCatch(sum(terra::expanse(v, unit = "km"), na.rm = TRUE), error = function(e) NA_real_)
  vu <- tryCatch(terra::aggregate(v), error = function(e) NULL)
  area_union <- if (!is.null(vu)) tryCatch(sum(terra::expanse(vu, unit = "km"), na.rm = TRUE), error = function(e) NA_real_) else NA_real_
  overlap_fraction <- if (is.finite(area_parts) && area_parts > 0 && is.finite(area_union)) max(0, (area_parts - area_union) / area_parts) else NA_real_
  max_overlap <- as.numeric(thresholds$max_region_overlap_fraction %||% 1e-6)
  metrics[[length(metrics) + 1L]] <- real_audit_metric(
    "regions", "region_polygon_overlap_fraction",
    if (!is.finite(overlap_fraction)) "WARN" else if (overlap_fraction <= max_overlap) "PASS" else "FAIL",
    signif(overlap_fraction, 6), max_overlap,
    "Overlapping regions make cell-center zone assignment order-dependent."
  )

  std <- data.frame(feature = seq_len(nrow(v)), stringsAsFactors = FALSE)
  if (id_present && ids_complete) {
    meta <- region_metadata(v, cfg, require_unique_id = FALSE)
    std <- meta
  } else {
    std$region <- NA_character_; std$region_label <- NA_character_
    std$region_group <- NA_character_; std$region_group_label <- NA_character_
  }
  tab <- cbind(std, attrs)
  list(vector = v, table = tab, candidates = candidates, metrics = metrics)
}


real_audit_region_mask <- function(cfg, template, regions = NULL) {
  assert_packages("terra")
  if (is.null(regions)) {
    z <- terra::rast(template); terra::values(z) <- 1
    return(z)
  }
  u <- terra::aggregate(regions); u$.__audit_mask__ <- 1L
  terra::rasterize(u, template, field = ".__audit_mask__", background = NA, touches = FALSE)
}

real_audit_complete_mask <- function(stack, filename = "") {
  assert_packages("terra")
  nyr <- terra::nlyr(stack)
  fun <- function(v) {
    if (is.matrix(v)) return(as.integer(rowSums(is.finite(v)) == nyr))
    as.integer(sum(is.finite(v)) == nyr)
  }
  terra::app(stack, fun = fun, cores = 1L, filename = filename,
             overwrite = TRUE, wopt = list(datatype = "INT1U"))
}

real_audit_extent_from_rows_cols <- function(template, row1, row2, col1, col2) {
  e <- terra::ext(template); rr <- terra::res(template)
  left <- e$xmin + (col1 - 1L) * rr[[1L]]
  right <- e$xmin + col2 * rr[[1L]]
  top <- e$ymax - (row1 - 1L) * rr[[2L]]
  bottom <- e$ymax - row2 * rr[[2L]]
  terra::ext(left, right, bottom, top)
}

real_audit_candidate_windows <- function(template, complete_combined, zone_raster, side_cells, candidate_grid = 5L) {
  nr <- terra::nrow(template); nc <- terra::ncol(template)
  side_r <- min(as.integer(side_cells), nr); side_c <- min(as.integer(side_cells), nc)
  row_max <- max(1L, nr - side_r + 1L); col_max <- max(1L, nc - side_c + 1L)
  rows <- unique(as.integer(round(seq(1, row_max, length.out = max(1L, as.integer(candidate_grid))))))
  cols <- unique(as.integer(round(seq(1, col_max, length.out = max(1L, as.integer(candidate_grid))))))
  out <- list(); k <- 0L
  for (r0 in rows) for (c0 in cols) {
    k <- k + 1L; r1 <- min(nr, r0 + side_r - 1L); c1 <- min(nc, c0 + side_c - 1L)
    ex <- real_audit_extent_from_rows_cols(template, r0, r1, c0, c1)
    cm <- terra::crop(complete_combined, ex, snap = "near")
    vv <- terra::values(cm, mat = FALSE)
    complete_fraction <- mean(is.finite(vv) & vv == 1)
    region_count <- 0L
    if (!is.null(zone_raster)) {
      zv <- terra::values(terra::crop(zone_raster, ex, snap = "near"), mat = FALSE)
      region_count <- length(unique(zv[is.finite(zv)]))
    }
    # Coverage dominates. Region diversity is a small secondary preference.
    score <- complete_fraction + 0.03 * min(region_count, 4L)
    out[[k]] <- data.frame(candidate = k, row1 = r0, row2 = r1, col1 = c0, col2 = c1,
                           complete_fraction = complete_fraction, region_count = region_count,
                           score = score, xmin = ex$xmin, xmax = ex$xmax, ymin = ex$ymin, ymax = ex$ymax,
                           stringsAsFactors = FALSE)
  }
  d <- do.call(rbind, out)
  d[order(-d$score, -d$complete_fraction, -d$region_count, d$candidate), , drop = FALSE]
}

real_audit_prepare_subset <- function(cfg, inventory, spec, project_dir, preflight_root, reuse = TRUE) {
  assert_packages(c("terra", "yaml"))
  responses <- names(inventory$files); years <- as.integer(unlist(cfg$data$years))
  stacks <- lapply(responses, function(nm) load_response_stack(inventory$files[[nm]], cfg)); names(stacks) <- responses
  template <- stacks[[1L]][[1L]]
  for (nm in responses[-1L]) if (!terra::compareGeom(template, stacks[[nm]][[1L]], stopOnError = FALSE))
    stop("Cross-response geometry differs after configured region masking: ", responses[[1L]], " vs ", nm, call. = FALSE)

  reg_diag <- real_audit_region_diagnostics(cfg, template, spec$thresholds %||% list())
  region_mask <- real_audit_region_mask(cfg, template, reg_diag$vector)
  study_cells <- as.numeric(terra::global(!is.na(region_mask), "sum", na.rm = TRUE)[[1L, 1L]])

  mask_dir <- ensure_dir(file.path(preflight_root, "complete_masks"))
  complete <- list(); coverage_rows <- list()
  for (nm in responses) {
    f <- file.path(mask_dir, paste0(nm, "_complete.tif"))
    complete[[nm]] <- real_audit_complete_mask(stacks[[nm]], f)
    complete[[nm]] <- terra::mask(complete[[nm]], region_mask)
    n_complete <- as.numeric(terra::global(complete[[nm]] == 1, "sum", na.rm = TRUE)[[1L, 1L]])
    coverage_rows[[nm]] <- data.frame(response = nm, study_cells = study_cells, complete_cells = n_complete,
                                      complete_fraction = if (study_cells > 0) n_complete / study_cells else NA_real_)
  }
  coverage <- do.call(rbind, coverage_rows)
  combined <- complete[[1L]]
  if (length(complete) > 1L) for (i in 2:length(complete)) combined <- combined * complete[[i]]
  combined <- terra::mask(combined, region_mask)
  combined_file <- file.path(mask_dir, "combined_complete.tif")
  combined <- terra::writeRaster(combined, combined_file, overwrite = TRUE, datatype = "INT1U")
  combined_complete <- as.numeric(terra::global(combined == 1, "sum", na.rm = TRUE)[[1L, 1L]])
  combined_fraction <- if (study_cells > 0) combined_complete / study_cells else NA_real_

  zones <- NULL
  if (!is.null(reg_diag$vector) && nzchar(region_field_config(cfg)$id) && region_field_config(cfg)$id %in% names(reg_diag$vector)) {
    zz <- reg_diag$vector; zz$.__audit_zone__ <- seq_len(nrow(zz))
    zones <- terra::rasterize(zz, template, field = ".__audit_zone__", background = NA, touches = FALSE)
  }

  side <- real_audit_side_cells(spec); cg <- as.integer(spec$subset$candidate_grid %||% 5L)
  candidates <- real_audit_candidate_windows(template, combined, zones, side, cg)
  if (!nrow(candidates)) stop("No subset candidate windows could be constructed.", call. = FALSE)
  chosen <- candidates[1L, , drop = FALSE]
  ex <- terra::ext(chosen$xmin, chosen$xmax, chosen$ymin, chosen$ymax)

  subset_key <- paste0(inventory$fingerprint, "_s", side, "_r", chosen$row1, "_c", chosen$col1)
  data_root <- ensure_dir(file.path(project_dir, "validation", "real_subset_data", subset_key))
  marker <- file.path(data_root, "subset_manifest.csv")
  if (!(isTRUE(reuse) && file.exists(marker))) {
    if (dir.exists(data_root)) unlink(data_root, recursive = TRUE, force = TRUE)
    ensure_dir(data_root)
    manifest_rows <- list()
    for (nm in responses) {
      od <- ensure_dir(file.path(data_root, nm))
      tab <- inventory$files[[nm]]
      for (i in seq_len(nrow(tab))) {
        src <- terra::rast(tab$file[[i]])
        z <- terra::crop(src, ex, snap = "near")
        out <- file.path(od, basename(tab$file[[i]]))
        terra::writeRaster(z, out, overwrite = TRUE)
        manifest_rows[[length(manifest_rows) + 1L]] <- data.frame(response = nm, year = tab$year[[i]], source = tab$file[[i]], subset_file = out, stringsAsFactors = FALSE)
      }
    }
    subset_regions <- ""
    if (!is.null(reg_diag$vector)) {
      rv <- terra::crop(reg_diag$vector, ex)
      if (nrow(rv) < 1L && isTRUE(cfg$data$mask_to_regions)) stop("Chosen subset does not intersect configured regions.", call. = FALSE)
      rd <- ensure_dir(file.path(data_root, "regions")); subset_regions <- file.path(rd, "subset_regions.shp")
      if (nrow(rv)) terra::writeVector(rv, subset_regions, overwrite = TRUE, filetype = "ESRI Shapefile")
    }
    manifest <- do.call(rbind, manifest_rows)
    manifest$subset_regions <- subset_regions
    write_csv(manifest, marker)
  }

  manifest <- utils::read.csv(marker, stringsAsFactors = FALSE)
  list(data_root = data_root, manifest = manifest, subset_regions = unique(manifest$subset_regions)[1L],
       extent = ex, chosen = chosen, candidates = candidates, coverage = coverage,
       combined_fraction = combined_fraction, study_cells = study_cells, region_diagnostics = reg_diag,
       source_fingerprint = inventory$fingerprint)
}

real_audit_apply_mode_overrides <- function(cfg_clean, spec, mode) {
  if (mode == "full") return(cfg_clean)
  ov <- if (mode == "quick") spec$quick_overrides else spec$standard_overrides
  if (mode == "quick") {
    cfg_clean$rho_estimation$sample_n <- as.integer(ov$rho_sample_n)
    cfg_clean$rho_estimation$bias_correction$simulation_B <- as.integer(ov$rho_bias_sim_B)
    cfg_clean$supf$mc_B <- as.integer(ov$supf_mc_B)
    cfg_clean$supf$mc_chunk <- as.integer(ov$supf_mc_chunk)
  }
  cfg_clean$bootstrap$pixel$B <- as.integer(ov$pixel_bootstrap_B)
  cfg_clean$bootstrap$pixel$min_valid_replicates_for_ci <- as.integer(ov$pixel_bootstrap_min_ci)
  cfg_clean$bootstrap$trajectory$B <- as.integer(ov$trajectory_B)
  cfg_clean$bootstrap$trajectory$supf_mc_B <- as.integer(ov$trajectory_supf_mc_B)
  cfg_clean$bootstrap$trajectory$rho_bias_sim_B <- as.integer(ov$trajectory_rho_bias_sim_B)
  cfg_clean$bootstrap$trajectory$min_valid_replicates_for_ci <- as.integer(ov$trajectory_min_ci)
  cfg_clean
}

real_audit_write_subset_config <- function(cfg, subset, spec, mode, project_dir) {
  assert_packages("yaml")
  clean <- cfg[!grepl("^\\.", names(cfg))]
  responses <- names(cfg$data$responses)
  for (nm in responses) clean$data$responses[[nm]]$folder <- file.path(subset$data_root, nm)
  if (nzchar(subset$subset_regions %||% "") && file.exists(subset$subset_regions)) clean$data$regions_shapefile <- subset$subset_regions
  clean$data$input_version <- paste0(cfg$data$input_version %||% "", "|REAL_SUBSET_AUDIT|", basename(subset$data_root))
  clean$output$root <- file.path(project_dir, "real_subset_audit_outputs_v5_6")
  clean$supf$cache_dir <- "_cache/supf"
  clean$runtime$overwrite <- FALSE
  subset_cells <- (as.integer(subset$chosen$row2[[1L]]) - as.integer(subset$chosen$row1[[1L]]) + 1L) *
                  (as.integer(subset$chosen$col2[[1L]]) - as.integer(subset$chosen$col1[[1L]]) + 1L)
  clean$bootstrap$pixel$max_target_pixels <- max(as.numeric(clean$bootstrap$pixel$max_target_pixels %||% 0), subset_cells)
  clean$validation$run_simulation_null <- FALSE
  clean$validation$run_simulation_power <- FALSE
  clean$validation$run_bootstrap_coverage <- FALSE
  clean <- real_audit_apply_mode_overrides(clean, spec, mode)
  cfg_dir <- ensure_dir(file.path(project_dir, "validation", "real_subset_generated"))
  path <- file.path(cfg_dir, paste0(basename(subset$data_root), "_", mode, ".yml"))
  yaml::write_yaml(clean, path)
  path
}

real_audit_preflight <- function(cfg, inventory, subset, spec, preflight_root) {
  assert_packages("terra")
  metrics <- list(); add <- function(x) metrics[[length(metrics) + 1L]] <<- x
  responses <- names(inventory$files)
  geom <- do.call(rbind, lapply(responses, function(nm) real_audit_compare_response_geometry(inventory$files[[nm]], nm)))
  write_csv(geom, file.path(preflight_root, "full_input_geometry.csv"))
  for (nm in responses) add(real_audit_metric("geometry", paste0(nm, "_all_years_geometry_match"),
    if (all(geom$geometry_matches_first_year[geom$response == nm])) "PASS" else "FAIL",
    sum(geom$geometry_matches_first_year[geom$response == nm]), nrow(inventory$files[[nm]])))

  refs <- lapply(responses, function(nm) terra::rast(inventory$files[[nm]]$file[[1L]])); names(refs) <- responses
  cross_ok <- all(vapply(responses[-1L], function(nm) terra::compareGeom(refs[[1L]], refs[[nm]], stopOnError = FALSE), logical(1L)))
  add(real_audit_metric("geometry", "cross_response_geometry_match", if (cross_ok) "PASS" else "FAIL", as.integer(cross_ok), 1))
  sig <- do.call(rbind, lapply(responses, function(nm) cbind(response = nm, real_audit_geom_signature(refs[[nm]]))))
  write_csv(sig, file.path(preflight_root, "full_input_geometry_signature.csv"))
  for (nm in responses) {
    has_crs <- nzchar(terra::crs(refs[[nm]]))
    add(real_audit_metric("geometry", paste0(nm, "_crs_present"), if (has_crs) "PASS" else "FAIL", as.integer(has_crs), 1))
  }

  stats <- do.call(rbind, lapply(responses, function(nm) real_audit_layer_stats(inventory$files[[nm]], nm, cfg)))
  write_csv(stats, file.path(preflight_root, "full_input_layer_stats.csv"))
  input_rows <- inventory$manifest[inventory$manifest$response %in% responses, , drop = FALSE]
  add(real_audit_metric("performance", "full_input_raster_size_gb", "INFO", signif(sum(input_rows$size, na.rm = TRUE) / 1024^3, 6), NA))
  for (nm in responses) {
    d <- stats[stats$response == nm, , drop = FALSE]
    finite_stats <- all(is.finite(d$min) & is.finite(d$max) & is.finite(d$mean))
    add(real_audit_metric("values", paste0(nm, "_finite_annual_summary"), if (finite_stats) "PASS" else "WARN",
                          sum(is.finite(d$mean)), nrow(d), "Annual global min/max/mean after NA removal."))
    neg_years <- d$year[is.finite(d$min) & d$min < 0]
    if (tolower(nm) %in% c("sr", "shannon"))
      add(real_audit_metric("values", paste0(nm, "_negative_values"), if (!length(neg_years)) "PASS" else "WARN",
                            length(neg_years), 0, if (length(neg_years)) paste("Negative annual minima in", paste(neg_years, collapse = ",")) else ""))
    vr <- diff(range(d$valid_fraction_raster, na.rm = TRUE)); vlim <- as.numeric(spec$thresholds$annual_valid_fraction_range_warning %||% 0.05)
    add(real_audit_metric("coverage", paste0(nm, "_annual_valid_fraction_range"),
                          if (!is.finite(vr) || vr <= vlim) "PASS" else "WARN", signif(vr, 5), vlim,
                          "Large annual changes in valid-cell coverage can mimic temporal change."))
    mm <- stats::median(d$mean, na.rm = TRUE); md <- stats::mad(d$mean, center = mm, constant = 1, na.rm = TRUE)
    mz <- if (is.finite(md) && md > 0) max(abs(d$mean - mm) / md, na.rm = TRUE) else 0
    mzlim <- as.numeric(spec$thresholds$annual_mean_mad_z_warning %||% 8)
    add(real_audit_metric("values", paste0(nm, "_annual_mean_max_robust_z"),
                          if (!is.finite(mz) || mz <= mzlim) "PASS" else "WARN", signif(mz, 5), mzlim,
                          "Robust diagnostic for a single anomalous annual raster; ecological shifts still require scientific interpretation."))
    duplicate_md5 <- anyDuplicated(inventory$manifest$md5[inventory$manifest$response == nm]) > 0L
    add(real_audit_metric("inputs", paste0(nm, "_duplicate_annual_file_md5"), if (!duplicate_md5) "PASS" else "WARN",
                          as.integer(duplicate_md5), 0, "Identical annual raster bytes can be legitimate but should be inspected."))
  }

  thr_each <- as.numeric(spec$thresholds$min_complete_fraction_each_response %||% 0.5)
  for (i in seq_len(nrow(subset$coverage))) {
    d <- subset$coverage[i, ]
    add(real_audit_metric("coverage", paste0(d$response, "_complete_series_fraction_full_study"),
                          if (is.finite(d$complete_fraction) && d$complete_fraction >= thr_each) "PASS" else "WARN",
                          signif(d$complete_fraction, 5), thr_each,
                          "Fraction of study-region cells with all configured years finite."))
  }
  thr_comb <- as.numeric(spec$thresholds$min_combined_complete_fraction %||% 0.4)
  add(real_audit_metric("coverage", "combined_response_complete_fraction_full_study",
                        if (is.finite(subset$combined_fraction) && subset$combined_fraction >= thr_comb) "PASS" else "WARN",
                        signif(subset$combined_fraction, 5), thr_comb,
                        "Cells complete for every year in every configured response."))
  preferred <- as.numeric(spec$subset$min_preferred_complete_fraction %||% 0.50)
  chosen_cf <- as.numeric(subset$chosen$complete_fraction[[1L]])
  add(real_audit_metric("subset", "chosen_window_combined_complete_fraction",
                        if (is.finite(chosen_cf) && chosen_cf >= preferred) "PASS" else "WARN",
                        signif(chosen_cf, 5), preferred, "The selector maximizes complete-series coverage first and region diversity second."))

  for (m in subset$region_diagnostics$metrics) add(m)
  write_csv(subset$region_diagnostics$table, file.path(preflight_root, "region_features.csv"))
  write_csv(subset$region_diagnostics$candidates, file.path(preflight_root, "region_id_candidates.csv"))
  write_csv(subset$candidates, file.path(preflight_root, "subset_selection_candidates.csv"))
  write_csv(subset$chosen, file.path(preflight_root, "subset_selection_chosen.csv"))
  write_csv(subset$coverage, file.path(preflight_root, "full_input_complete_series_coverage.csv"))

  tab <- do.call(rbind, metrics)
  write_csv(tab, file.path(preflight_root, "preflight_metrics.csv"))
  tab
}

real_audit_capture_pipeline <- function(config_path, project_dir) {
  warnings_seen <- character()
  value <- withCallingHandlers(
    run_pipeline_until("validation", config_path = config_path, project_dir = project_dir),
    warning = function(w) { warnings_seen <<- c(warnings_seen, conditionMessage(w)) }
  )
  list(run = value, warnings = unique(warnings_seen))
}

real_audit_evaluate_run <- function(res, warnings_seen, preflight, subset, spec, mode, preflight_seconds = NA_real_, pipeline_seconds = NA_real_) {
  metrics <- list(); add <- function(x) metrics[[length(metrics) + 1L]] <<- x
  for (i in seq_len(nrow(preflight))) metrics[[length(metrics) + 1L]] <- preflight[i, , drop = FALSE]
  add(real_audit_metric("performance", "full_preflight_and_subset_preparation_seconds", "INFO", signif(preflight_seconds, 6), NA))
  add(real_audit_metric("performance", "subset_pipeline_seconds", "INFO", signif(pipeline_seconds, 6), NA))
  out_files <- list.files(res$ctx$root, recursive = TRUE, full.names = TRUE); out_files <- out_files[file.exists(out_files) & !dir.exists(out_files)]
  out_gb <- if (length(out_files)) sum(file.info(out_files)$size, na.rm = TRUE) / 1024^3 else 0
  add(real_audit_metric("performance", "subset_run_output_size_gb", "INFO", signif(out_gb, 6), NA))

  s <- res$primary$summary
  if (!is.null(s) && nrow(s)) for (i in seq_len(nrow(s))) {
    d <- s[i, ]; r <- as.character(d$response)
    add(real_audit_metric("primary", paste0(r, "_eligible_pixels"), if (d$eligible_pixels > 0) "PASS" else "FAIL",
                          d$eligible_pixels, ">0"))
    add(real_audit_metric("fdr", paste0(r, "_mc_floor_min_rank"), if (d$min_rank_at_mc_floor <= 1) "PASS" else "INFO",
                          d$min_rank_at_mc_floor, 1,
                          paste0("MC p-min=", signif(d$mc_min_p, 5), "; subset FDR family only.")))
    bb <- c(d$lower_boundary_share, d$upper_boundary_share); bb <- bb[is.finite(bb)]
    bmax <- if (length(bb)) max(bb) else NA_real_
    bw <- as.numeric(spec$thresholds$boundary_share_warning %||% 0.25)
    add(real_audit_metric("science", paste0(r, "_primary_boundary_share_max"),
                          if (!is.finite(bmax) || bmax <= bw) "PASS" else "WARN",
                          signif(bmax, 5), bw,
                          "A high share at candidate-window boundaries requires 5/4/3/2 interpretation."))
  }

  rho <- res$rho
  if (!is.null(rho) && nrow(rho)) for (i in seq_len(nrow(rho))) {
    d <- rho[i, ]; add(real_audit_metric("rho", paste0(d$response, "_corrected_rho_hits_bound"),
      if (!isTRUE(d$correction_hit_bound)) "PASS" else "WARN", as.integer(isTRUE(d$correction_hit_bound)), 0))
  }
  comp_file <- file.path(res$ctx$root, "rho", "rho_model_comparison.csv")
  if (file.exists(comp_file)) {
    comp <- utils::read.csv(comp_file, stringsAsFactors = FALSE)
    lim <- as.numeric(spec$thresholds$rho_model_difference_warning %||% 0.10)
    if ("linear_minus_hinge" %in% names(comp)) for (i in seq_len(nrow(comp))) {
      dd <- abs(comp$linear_minus_hinge[[i]])
      add(real_audit_metric("rho", paste0(comp$response[[i]], "_linear_vs_hinge_abs_difference"),
                            if (is.finite(dd) && dd <= lim) "PASS" else "WARN", signif(dd, 5), lim,
                            "Inspect rho sensitivity before publication if this is large."))
    }
  }

  # Exact response-specific linear-residual rho is a stress test only. Report
  # its impact explicitly so a large model disagreement cannot be hidden by
  # the coarser fixed-rho sensitivity grid.
  linear_stress <- res$rho_sensitivity[["linear_diagnostic"]] %||% NULL
  if (!is.null(linear_stress) && !is.null(linear_stress$summary) && nrow(linear_stress$summary)) {
    ls <- linear_stress$summary
    ps <- res$primary$summary
    for (response in intersect(as.character(ls$response), as.character(ps$response))) {
      a <- ls[ls$response == response, , drop = FALSE][1L, ]
      b <- ps[ps$response == response, , drop = FALSE][1L, ]
      retention <- if (is.finite(b$recovery_breaks) && b$recovery_breaks > 0) a$recovery_breaks / b$recovery_breaks else NA_real_
      add(real_audit_metric("rho_stress", paste0(response, "_linear_diagnostic_rho"), "INFO",
                            signif(a$rho, 7), NA,
                            "Response-specific linear-residual corrected rho; diagnostic stress test, not the primary estimator."))
      add(real_audit_metric("rho_stress", paste0(response, "_linear_diagnostic_recovery_breaks"), "INFO",
                            a$recovery_breaks, NA,
                            "Subset FDR family only; use to quantify rho-model sensitivity, not as a paper result."))
      add(real_audit_metric("rho_stress", paste0(response, "_linear_diagnostic_recovery_retention_vs_primary"), "INFO",
                            signif(retention, 5), NA,
                            "Recovery count under exact linear-diagnostic rho divided by primary hinge-empirical recovery count."))
    }
  }

  if (!is.null(res$pixel) && length(res$pixel)) for (nm in names(res$pixel)) {
    d <- res$pixel[[nm]]; if (is.null(d) || !nrow(d)) next
    unstable <- mean(d$ci_status != "ok")
    lim <- as.numeric(spec$thresholds$bootstrap_unstable_fraction_warning %||% 0.50)
    add(real_audit_metric("bootstrap", paste0(nm, "_pixel_ci_unstable_fraction"),
                          if (unstable <= lim) "PASS" else "WARN", signif(unstable, 5), lim,
                          paste0("Audit mode=", mode, "; reduced bootstrap B may inflate this fraction.")))
    add(real_audit_metric("bootstrap", paste0(nm, "_median_recovery_frequency"), "INFO",
                          signif(stats::median(d$bootstrap_recovery_frequency, na.rm = TRUE), 5), NA,
                          "Selection-conditioned stability frequency, not a posterior probability."))
  }

  if (!is.null(res$trajectory) && nrow(res$trajectory$summary)) {
    cw <- sum(res$trajectory$summary$coverage_warning, na.rm = TRUE)
    add(real_audit_metric("trajectory", "regional_coverage_warnings", if (cw == 0) "PASS" else "WARN", cw, 0,
                          "Subset regional trajectories are audit diagnostics only; clipping changes region composition."))
  }

  if (!is.null(res$spatial)) for (nm in names(res$spatial)) {
    d <- res$spatial[[nm]]$scenario_iou
    if (!is.null(d) && nrow(d) && "within1_year_agreement" %in% names(d)) {
      x <- stats::median(d$within1_year_agreement, na.rm = TRUE); lim <- as.numeric(spec$thresholds$spatial_within1_year_warning %||% 0.50)
      add(real_audit_metric("spatial", paste0(nm, "_median_boundary_scenario_within1yr_agreement"),
                            if (!is.finite(x)) "INFO" else if (x >= lim) "PASS" else "WARN", signif(x, 5), lim,
                            "Separates temporal break-year stability from mask IoU."))
    }
  }

  expected <- c(file.path(res$ctx$root, "tables", "Table1A_primary_overall_summary.csv"),
                file.path(res$ctx$root, "tables", "Table2A_boundary_sensitivity_overall.csv"),
                file.path(res$ctx$root, "rho", "rho_estimates.csv"),
                file.path(res$ctx$root, "validation", "validation_report.md"))
  add(real_audit_metric("workflow", "essential_artifacts_present", if (all(file.exists(expected))) "PASS" else "FAIL",
                        sum(file.exists(expected)), length(expected), paste(basename(expected[!file.exists(expected)]), collapse = ";")))

  if (length(warnings_seen)) {
    for (i in seq_along(warnings_seen)) {
      known_rho <- grepl("Linear- and hinge-residual corrected rho differ materially", warnings_seen[[i]], fixed = TRUE)
      add(real_audit_metric("runtime_warning", paste0("warning_", i), if (known_rho) "INFO" else "WARN", warnings_seen[[i]], NA,
        if (known_rho) "Expected scientific QA warning; quantified separately in rho metrics." else "Warnings are surfaced, not suppressed; classify before full-data analysis."))
    }
  } else add(real_audit_metric("runtime_warning", "runtime_warning_count", "PASS", 0, 0))

  do.call(rbind, metrics)
}

real_audit_write_report <- function(res, metrics, subset, config_path, mode, source_config, preflight_root) {
  root <- file.path(res$ctx$root, "validation")
  write_csv(metrics, file.path(root, "real_data_subset_audit_metrics.csv"))
  # Preserve the preflight/selection evidence with the run itself.
  pf <- list.files(preflight_root, full.names = TRUE)
  if (length(pf)) file.copy(pf[file.info(pf)$isdir %in% FALSE], root, overwrite = TRUE)
  fail <- metrics$status == "FAIL"; warn <- metrics$status == "WARN"
  status <- if (any(fail)) "FAIL" else if (any(warn)) "PASS_WITH_WARNINGS" else "PASS"
  summary <- list(status = status, failures = sum(fail), warnings = sum(warn), mode = mode,
                  source_config = source_config, subset_config = config_path,
                  subset_data_root = subset$data_root, subset_extent = as.vector(subset$extent),
                  run_root = res$ctx$root)
  write_json(summary, file.path(root, "real_data_subset_audit_summary.json"))

  failtab <- metrics[metrics$status == "FAIL", , drop = FALSE]
  warntab <- metrics[metrics$status == "WARN", , drop = FALSE]
  infotab <- metrics[metrics$status == "INFO", , drop = FALSE]
  nottested <- metrics[metrics$status == "NOT_TESTED", , drop = FALSE]
  rf <- region_field_config(res$cfg)
  lines <- c(
    "# v5.6 real-data subset audit", "",
    paste0("- Status: **", status, "**"),
    paste0("- Audit mode: `", mode, "`"),
    paste0("- Source config: `", source_config, "`"),
    paste0("- Generated subset config: `", config_path, "`"),
    paste0("- Full-input fingerprint: `", subset$source_fingerprint, "`"),
    paste0("- Subset-input fingerprint: `", res$inventory$fingerprint, "`"),
    paste0("- Subset data root: `", subset$data_root, "`"),
    paste0("- Subset run root: `", res$ctx$root, "`"),
    paste0("- Chosen subset complete fraction (all responses/all years): ", signif(subset$chosen$complete_fraction[[1L]], 5)),
    paste0("- Regions represented in chosen window: ", subset$chosen$region_count[[1L]]),
    paste0("- Region feature ID: `", rf$id, "`; label: `", rf$label, "`; group: `", rf$group, "`; group label: `", rf$group_label, "`"),
    "",
    "## Scope and interpretation", "",
    "This audit first inspects the complete real input collection, then runs a spatial crop through the complete v5.6 workflow.",
    "**Subset BH/BY results are not the full-study FDR family and must not be reported as final scientific inference.**",
    "Likewise, trajectories from clipped region fragments are implementation/data diagnostics, not final regional trajectories.",
    "The audit is intended to expose geometry, CRS, NoData, value-range, rho, boundary, bootstrap, spatial, RAM/I/O, and wiring problems before the full production run.",
    "",
    "## Failures", "",
    if (nrow(failtab)) paste(capture.output(print(failtab, row.names = FALSE)), collapse = "\n") else "None.",
    "", "## Warnings requiring review", "",
    if (nrow(warntab)) paste(capture.output(print(warntab, row.names = FALSE)), collapse = "\n") else "None.",
    "", "## Not-tested dependent checks", "",
    if (nrow(nottested)) paste(capture.output(print(nottested, row.names = FALSE)), collapse = "\n") else "None.",
    "", "## Informational diagnostics", "",
    if (nrow(infotab)) paste(capture.output(print(infotab, row.names = FALSE)), collapse = "\n") else "None.",
    "", "## Primary subset summary", "",
    paste(capture.output(print(res$primary$summary, row.names = FALSE)), collapse = "\n"),
    "", "## Response-specific rho", "",
    paste(capture.output(print(res$rho, row.names = FALSE)), collapse = "\n"),
    "", "## Exact linear-diagnostic rho stress test", "",
    if (!is.null(res$rho_sensitivity[["linear_diagnostic"]]))
      paste(capture.output(print(res$rho_sensitivity[["linear_diagnostic"]]$summary, row.names = FALSE)), collapse = "\n")
    else "Not available (linear_diagnostic was not configured or diagnostic rho values were unavailable).",
    "", "## Next gate", "",
    "If there are no FAIL items, inspect every WARN item, the generated maps, rho_model_comparison.csv, the exact linear_diagnostic rho stress scenario, boundary sensitivity tables, MC floor diagnostics, and pixel-bootstrap stability before running the full study area."
  )
  report <- file.path(root, "REAL_DATA_SUBSET_AUDIT_REPORT.md")
  writeLines(lines, report, useBytes = TRUE)
  list(status = status, report = report, summary = file.path(root, "real_data_subset_audit_summary.json"), metrics = metrics)
}

run_real_data_subset_audit_v56 <- function(project_dir = getwd(), source_config = Sys.getenv("BREAK_CONFIG", "config/v5_6.yml"),
                                           reuse_subset = NULL) {
  project_dir <- normalizePath(project_dir, winslash = "/", mustWork = TRUE)
  spec <- real_audit_read_spec(project_dir); mode <- real_audit_mode(spec)
  if (is.null(reuse_subset)) {
    env_reuse <- tolower(Sys.getenv("REAL_SUBSET_REUSE", unset = as.character(spec$subset$reuse_existing_subset %||% TRUE)))
    reuse_subset <- env_reuse %in% c("true", "1", "yes", "y")
  }
  t_pre <- proc.time()[["elapsed"]]
  cfg <- read_config(source_config, project_dir)
  inventory <- build_input_inventory(cfg)
  preflight_root <- ensure_dir(file.path(project_dir, "validation", "real_subset_preflight", inventory$fingerprint))
  log_msg("Real-data audit input fingerprint: ", inventory$fingerprint)
  subset <- real_audit_prepare_subset(cfg, inventory, spec, project_dir, preflight_root, reuse = reuse_subset)
  preflight <- real_audit_preflight(cfg, inventory, subset, spec, preflight_root)
  preflight_seconds <- proc.time()[["elapsed"]] - t_pre
  if (any(preflight$status == "FAIL")) {
    write_csv(preflight, file.path(preflight_root, "preflight_metrics.csv"))
    stop("Real-data full-input preflight has FAIL item(s). Inspect: ", file.path(preflight_root, "preflight_metrics.csv"), call. = FALSE)
  }
  subset_config <- real_audit_write_subset_config(cfg, subset, spec, mode, project_dir)
  log_msg("Chosen real-data subset: side~", real_audit_side_cells(spec), " cells; complete fraction=",
          signif(subset$chosen$complete_fraction[[1L]], 5), "; regions=", subset$chosen$region_count[[1L]])
  log_msg("Running subset E2E in audit mode '", mode, "' with config: ", subset_config)
  t_pipe <- proc.time()[["elapsed"]]
  cap <- real_audit_capture_pipeline(subset_config, project_dir)
  pipeline_seconds <- proc.time()[["elapsed"]] - t_pipe
  metrics <- real_audit_evaluate_run(cap$run, cap$warnings, preflight, subset, spec, mode, preflight_seconds, pipeline_seconds)
  out <- real_audit_write_report(cap$run, metrics, subset, subset_config, mode, cfg$.config_path, preflight_root)
  invisible(list(run = cap$run, subset = subset, preflight = preflight, warnings = cap$warnings, evaluation = out))
}
