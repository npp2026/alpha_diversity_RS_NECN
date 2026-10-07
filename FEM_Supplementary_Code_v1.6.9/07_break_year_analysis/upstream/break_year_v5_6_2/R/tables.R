regional_pixel_summary <- function(cfg, scenario) {
  assert_packages("terra")
  reg0 <- load_regions(cfg)
  if (is.null(reg0)) return(data.frame())
  rows <- list()
  for (response in names(scenario$results)) {
    pr <- scenario$results[[response]]; template <- pr$layers$p_mc
    reg <- reg0
    if (!terra::same.crs(reg, template)) reg <- terra::project(reg, terra::crs(template))
    meta <- region_metadata(reg, cfg, require_unique_id = TRUE)
    reg$.__zone_id__ <- seq_len(nrow(reg))
    zones <- terra::rasterize(reg, template, field = ".__zone_id__", touches = FALSE)
    elig <- terra::zonal(terra::ifel(!is.na(pr$layers$p_mc), 1, NA), zones, fun = "sum", na.rm = TRUE)
    sig <- terra::zonal(pr$significant, zones, fun = "sum", na.rm = TRUE)
    rec <- terra::zonal(pr$recovery, zones, fun = "sum", na.rm = TRUE)
    ar <- terra::zonal(terra::ifel(pr$recovery == 1, terra::cellSize(template, unit = "km"), NA), zones, fun = "sum", na.rm = TRUE)
    names(elig) <- c("zone_id", "eligible_pixels"); names(sig) <- c("zone_id", "significant_breaks")
    names(rec) <- c("zone_id", "recovery_breaks"); names(ar) <- c("zone_id", "recovery_area_km2")
    z <- Reduce(function(x, y) merge(x, y, by = "zone_id", all = TRUE), list(elig, sig, rec, ar))
    br <- terra::mask(pr$layers$break_year, terra::ifel(pr$recovery == 1, 1, NA))
    ct <- terra::crosstab(c(zones, br), long = TRUE); bst <- data.frame()
    if (!is.null(ct) && nrow(ct)) {
      names(ct)[1:3] <- c("zone_id", "break_year", "count")
      ct <- ct[is.finite(ct$zone_id) & is.finite(ct$break_year) & ct$count > 0, , drop = FALSE]
      if (nrow(ct)) {
        ss <- split(ct, ct$zone_id)
        bst <- do.call(rbind, lapply(ss, function(d) {
          mx <- max(d$count)
          data.frame(zone_id = d$zone_id[[1L]], modal_break_year = min(d$break_year[d$count == mx]),
                     median_break_year = weighted_discrete_quantile_from_counts(d$break_year, d$count, .5))
        }))
      }
    }
    if (nrow(bst)) z <- merge(z, bst, by = "zone_id", all = TRUE)
    else { z$modal_break_year <- NA_real_; z$median_break_year <- NA_real_ }
    idx <- match(as.integer(z$zone_id), meta$feature)
    z$region <- meta$region[idx]; z$region_label <- meta$region_label[idx]
    z$region_group <- meta$region_group[idx]; z$region_group_label <- meta$region_group_label[idx]
    z$response <- response; z$scenario_id <- scenario$id; z$min_segment <- pr$summary$min_segment[[1L]]
    z$rho <- pr$summary$rho[[1L]]
    z$recovery_pct_eligible <- ifelse(z$eligible_pixels > 0, 100 * z$recovery_breaks / z$eligible_pixels, NA_real_)
    rows[[response]] <- z[, c("scenario_id", "min_segment", "response", "region", "region_label",
      "region_group", "region_group_label", "eligible_pixels", "significant_breaks", "recovery_breaks",
      "recovery_pct_eligible", "recovery_area_km2", "modal_break_year", "median_break_year", "rho")]
  }
  if (!length(rows)) data.frame() else do.call(rbind, rows)
}


build_publication_tables <- function(cfg,primary,boundary,multiplicity,rho_tab,pixel,trajectory,spatial,ctx){
  root<-ensure_dir(file.path(ctx$root,"tables"));write_csv(primary$summary,file.path(root,"Table1A_primary_overall_summary.csv"));rp<-regional_pixel_summary(cfg,primary);if(nrow(rp))write_csv(rp,file.path(root,"Table1B_primary_regional_summary.csv"))
  if(length(boundary)){d<-do.call(rbind,lapply(boundary,function(x)x$summary));write_csv(d,file.path(root,"Table2A_boundary_sensitivity_overall.csv"));dr<-do.call(rbind,lapply(boundary,function(x)regional_pixel_summary(cfg,x)));if(nrow(dr))write_csv(dr,file.path(root,"Table2B_boundary_sensitivity_regional.csv"))}
  if(!is.null(trajectory)&&nrow(trajectory$summary))write_csv(trajectory$summary,file.path(root,"Table3_regional_trajectory_uncertainty.csv"))
  rows<-list();if(!is.null(spatial))for(response in names(spatial)){cl<-spatial[[response]]$clusters;if(!is.null(cl)&&nrow(cl$summary)){d<-cl$summary;d$response<-response;rows[[response]]<-d}}
  if(length(rows))write_csv(do.call(rbind,rows),file.path(root,"Table4_spatial_clusters.csv"))
  write_csv(rho_tab,file.path(root,"Table5_rho_calibration.csv"));if(!is.null(multiplicity)&&nrow(multiplicity))write_csv(multiplicity,file.path(root,"TableS_multiplicity_sensitivity.csv"));invisible(TRUE)
}
