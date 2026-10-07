read_synthetic_spec <- function(path, project_dir = getwd()) {
  assert_packages("yaml")
  path <- resolve_path(path, project_dir)
  if (!file.exists(path)) stop("Synthetic spec not found: ", path, call. = FALSE)
  spec <- yaml::read_yaml(path)
  spec$.path <- path
  validate_synthetic_spec(spec)
  spec
}

validate_synthetic_spec <- function(spec) {
  yrs <- as.integer(unlist(spec$years))
  if (length(yrs) < 8L || anyDuplicated(yrs) || is.unsorted(yrs)) stop("Synthetic years must be sorted and unique.", call. = FALSE)
  nr <- as.integer(spec$raster$nrow); nc <- as.integer(spec$raster$ncol)
  if (nr < 4L || nc < 4L) stop("Synthetic raster is too small.", call. = FALSE)
  cls <- vapply(spec$classes, function(z) as.integer(z$code), integer(1L))
  nms <- vapply(spec$classes, function(z) as.character(z$name), character(1L))
  if (anyDuplicated(cls) || anyDuplicated(nms) || !0L %in% cls) stop("Synthetic classes need unique code/name and code 0.", call. = FALSE)
  for (z in spec$layout) {
    if (!as.character(z$class) %in% nms) stop("Unknown class in synthetic layout: ", z$class, call. = FALSE)
    rr <- as.integer(unlist(z$rows)); cc <- as.integer(unlist(z$cols))
    if (length(rr) != 2L || length(cc) != 2L || rr[[1L]] < 1L || rr[[2L]] > nr || cc[[1L]] < 1L || cc[[2L]] > nc) stop("Invalid synthetic layout bounds.", call. = FALSE)
  }
  invisible(TRUE)
}

synthetic_class_table <- function(spec) {
  do.call(rbind, lapply(spec$classes, function(z) {
    data.frame(code = as.integer(z$code), name = as.character(z$name),
      slope_before = as.numeric(z$slope_before), slope_after = as.numeric(z$slope_after),
      break_year = if (is.null(z$break_year)) NA_integer_ else as.integer(z$break_year),
      noise = isTRUE(z$noise), missing_year = if (is.null(z$missing_year)) NA_integer_ else as.integer(z$missing_year),
      stringsAsFactors = FALSE)
  }))
}

synthetic_layout_matrix <- function(spec) {
  nr <- as.integer(spec$raster$nrow); nc <- as.integer(spec$raster$ncol)
  m <- matrix(0L, nrow = nr, ncol = nc)
  ct <- synthetic_class_table(spec); by_name <- setNames(ct$code, ct$name)
  occupied <- matrix(FALSE, nrow = nr, ncol = nc)
  for (z in spec$layout) {
    rr <- seq.int(as.integer(z$rows[[1L]]), as.integer(z$rows[[2L]]))
    cc <- seq.int(as.integer(z$cols[[1L]]), as.integer(z$cols[[2L]]))
    if (any(occupied[rr, cc])) stop("Synthetic layout rectangles overlap; make classes mutually exclusive.", call. = FALSE)
    m[rr, cc] <- as.integer(by_name[[as.character(z$class)]])
    occupied[rr, cc] <- TRUE
  }
  m
}

synthetic_region_matrix <- function(spec) {
  nr <- as.integer(spec$raster$nrow); nc <- as.integer(spec$raster$ncol)
  m <- matrix(NA_integer_, nrow = nr, ncol = nc)
  for (z in spec$regions) {
    rr <- seq.int(as.integer(z$rows[[1L]]), as.integer(z$rows[[2L]]))
    cc <- seq.int(as.integer(z$cols[[1L]]), as.integer(z$cols[[2L]]))
    if (any(is.finite(m[rr, cc]))) stop("Synthetic regions overlap.", call. = FALSE)
    m[rr, cc] <- as.integer(z$id)
  }
  if (any(!is.finite(m))) stop("Synthetic regions must cover the entire raster.", call. = FALSE)
  m
}

simulate_ar1_vector_keyed <- function(n, rho, marginal_sd, seed, rng_kind = "L'Ecuyer-CMRG") {
  set_rng(rng_kind, seed)
  as.numeric(simulate_ar1_matrix(1L, n, rho, marginal_sd))
}

synthetic_signal_for_cell <- function(years, class_row, response_cfg, row, col, master_seed, response, rng_kind = "L'Ecuyer-CMRG") {
  scale <- as.numeric(response_cfg$slope_scale %||% 1)
  sb <- as.numeric(class_row$slope_before) * scale
  sa <- as.numeric(class_row$slope_after) * scale
  tau <- as.numeric(class_row$break_year)
  center <- mean(years)
  base <- as.numeric(response_cfg$baseline)
  # Stable spatial intercept heterogeneity exercises regional aggregation without
  # changing the known temporal model class.
  int_seed <- seed_from_key(master_seed, "synthetic_intercept", response, row, col)
  set_rng(rng_kind, int_seed)
  intercept <- base + stats::rnorm(1L, sd = as.numeric(response_cfg$intercept_sd))
  y <- intercept + sb * (years - center)
  if (is.finite(tau) && abs(sa - sb) > 0) y <- y + (sa - sb) * pmax(0, years - tau)
  if (isTRUE(class_row$noise)) {
    noise_seed <- seed_from_key(master_seed, "synthetic_noise", response, row, col)
    y <- y + simulate_ar1_vector_keyed(length(years), as.numeric(response_cfg$rho), as.numeric(response_cfg$noise_sd), noise_seed, rng_kind)
  }
  miss <- as.integer(class_row$missing_year)
  if (is.finite(miss)) y[years == miss] <- NA_real_
  y
}

write_synthetic_regions <- function(template, region_matrix, spec, out_dir, overwrite = TRUE) {
  assert_packages("terra"); ensure_dir(out_dir)
  zr <- template; terra::values(zr) <- as.vector(t(region_matrix))
  names(zr) <- "region_id"
  reg <- terra::as.polygons(zr, dissolve = TRUE, values = TRUE, na.rm = TRUE)
  region_defs <- do.call(rbind, lapply(spec$regions, function(z) data.frame(region_id = as.integer(z$id), NAME = as.character(z$name), stringsAsFactors = FALSE)))
  reg$NAME <- region_defs$NAME[match(as.integer(reg$region_id), region_defs$region_id)]
  shp <- file.path(out_dir, "synthetic_regions.shp")
  terra::writeVector(reg, shp, overwrite = overwrite, filetype = "ESRI Shapefile")
  shp
}

generate_synthetic_dataset <- function(project_dir = getwd(),
                                       spec_path = "validation/synthetic_spec.yml",
                                       out_dir = "validation/synthetic_data",
                                       overwrite = TRUE) {
  assert_packages(c("terra", "yaml", "digest"))
  project_dir <- normalizePath(project_dir, winslash = "/", mustWork = TRUE)
  spec <- read_synthetic_spec(spec_path, project_dir)
  out_dir <- resolve_path(out_dir, project_dir); ensure_dir(out_dir)
  years <- as.integer(unlist(spec$years)); nr <- as.integer(spec$raster$nrow); nc <- as.integer(spec$raster$ncol)
  res <- as.numeric(spec$raster$resolution_m); xmin <- as.numeric(spec$raster$xmin); ymin <- as.numeric(spec$raster$ymin)
  template <- terra::rast(nrows = nr, ncols = nc, xmin = xmin, xmax = xmin + nc * res,
                          ymin = ymin, ymax = ymin + nr * res, crs = as.character(spec$raster$crs))
  class_m <- synthetic_layout_matrix(spec); region_m <- synthetic_region_matrix(spec); ct <- synthetic_class_table(spec)
  truth_r <- template; terra::values(truth_r) <- as.vector(t(class_m)); names(truth_r) <- "truth_class"
  truth_tif <- file.path(out_dir, "synthetic_truth_class.tif")
  terra::writeRaster(truth_r, truth_tif, overwrite = overwrite, datatype = "INT1U")
  region_r <- template; terra::values(region_r) <- as.vector(t(region_m)); names(region_r) <- "region_id"
  terra::writeRaster(region_r, file.path(out_dir, "synthetic_region_id.tif"), overwrite = overwrite, datatype = "INT1U")
  shp <- write_synthetic_regions(template, region_m, spec, file.path(out_dir, "regions"), overwrite)

  xy <- terra::xyFromCell(template, seq_len(terra::ncell(template)))
  rc <- expand.grid(col = seq_len(nc), row = seq_len(nr)); rc <- rc[order(rc$row, rc$col), ]
  truth <- data.frame(cell = seq_len(terra::ncell(template)), row = rc$row, col = rc$col,
                      x = xy[,1], y = xy[,2], class_code = as.vector(t(class_m)), region_id = as.vector(t(region_m)), stringsAsFactors = FALSE)
  truth <- merge(truth, ct, by.x = "class_code", by.y = "code", all.x = TRUE, sort = FALSE)
  region_defs <- do.call(rbind, lapply(spec$regions, function(z) data.frame(region_id = as.integer(z$id), region = as.character(z$name), stringsAsFactors = FALSE)))
  truth <- merge(truth, region_defs, by = "region_id", all.x = TRUE, sort = FALSE); truth <- truth[order(truth$cell), ]
  write_csv(truth, file.path(out_dir, "synthetic_truth.csv"))
  write_csv(ct, file.path(out_dir, "synthetic_classes.csv"))

  files <- list(); master_seed <- as.integer(spec$seed); rng_kind <- as.character(spec$rng_kind %||% "L'Ecuyer-CMRG")
  for (response in names(spec$responses)) {
    rdir <- ensure_dir(file.path(out_dir, response)); vals <- matrix(NA_real_, nrow = terra::ncell(template), ncol = length(years))
    for (cell in seq_len(terra::ncell(template))) {
      row <- truth$row[[cell]]; col <- truth$col[[cell]]; cr <- ct[ct$code == truth$class_code[[cell]], , drop = FALSE]
      vals[cell, ] <- synthetic_signal_for_cell(years, cr, spec$responses[[response]], row, col, master_seed, response, rng_kind)
    }
    rf <- character(length(years))
    for (j in seq_along(years)) {
      rr <- template; terra::values(rr) <- vals[, j]; names(rr) <- paste0(response, "_", years[[j]])
      rf[[j]] <- file.path(rdir, sprintf("synthetic_%s_%d_1km.tif", response, years[[j]]))
      terra::writeRaster(rr, rf[[j]], overwrite = overwrite, datatype = "FLT8S", NAflag = -9999)
    }
    files[[response]] <- data.frame(year = years, file = rf, stringsAsFactors = FALSE)
  }
  spec_copy <- file.path(out_dir, "synthetic_spec_resolved.yml"); file.copy(spec$.path, spec_copy, overwrite = TRUE)
  manifest_files <- c(unlist(lapply(files, function(z) z$file)), truth_tif, file.path(out_dir, "synthetic_region_id.tif"),
                      list.files(file.path(out_dir, "regions"), full.names = TRUE), file.path(out_dir, "synthetic_truth.csv"), file.path(out_dir, "synthetic_classes.csv"), spec_copy)
  manifest_files <- manifest_files[file.exists(manifest_files)]
  write_csv(file_manifest(manifest_files, role = paste0("synthetic_", basename(manifest_files))), file.path(out_dir, "synthetic_manifest.csv"))
  metadata <- list(version = spec$version, seed = master_seed, rng_kind = rng_kind, years = years, nrow = nr, ncol = nc,
                   responses = lapply(spec$responses, function(z) list(rho = z$rho, noise_sd = z$noise_sd)),
                   regions_shapefile = shp, generated_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"))
  write_json(metadata, file.path(out_dir, "synthetic_metadata.json"))
  invisible(list(root = out_dir, truth = truth, files = files, regions = shp, spec = spec))
}
