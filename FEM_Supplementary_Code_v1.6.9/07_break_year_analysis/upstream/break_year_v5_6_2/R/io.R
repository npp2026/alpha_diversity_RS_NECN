discover_response_files <- function(response_cfg, years) {
  folder <- response_cfg$folder; pattern <- response_cfg$pattern
  if (!dir.exists(folder)) stop("Response folder not found: ", folder, call. = FALSE)
  files <- list.files(folder, pattern = pattern, full.names = TRUE)
  if (!length(files)) stop("No files matched ", pattern, " in ", folder, call. = FALSE)
  yrs <- suppressWarnings(as.integer(sub(pattern, "\\1", basename(files), perl = TRUE)))
  if (any(!is.finite(yrs))) stop("Could not extract a 4-digit year using response pattern: ", pattern, call. = FALSE)
  if (anyDuplicated(yrs)) stop("Duplicate years found in ", folder, call. = FALSE)
  missing <- setdiff(years, yrs); if (length(missing)) stop("Missing required years in ", folder, ": ", paste(missing, collapse = ", "), call. = FALSE)
  data.frame(year = years, file = files[match(years, yrs)], stringsAsFactors = FALSE)
}

load_regions <- function(cfg, template = NULL) {
  assert_packages("terra"); shp <- cfg$data$regions_shapefile
  if (!nzchar(shp) || !file.exists(shp)) return(NULL)
  v <- terra::vect(shp)
  if (!is.null(template) && !terra::same.crs(v, template)) v <- terra::project(v, terra::crs(template))
  v
}

load_response_stack <- function(files_table, cfg, mask_to_regions = isTRUE(cfg$data$mask_to_regions)) {
  assert_packages("terra")
  x <- terra::rast(files_table$file); names(x) <- paste0("y", files_table$year)
  if (mask_to_regions) {
    reg <- load_regions(cfg, x)
    if (!is.null(reg)) { u <- terra::aggregate(reg); x <- terra::crop(x, u); x <- terra::mask(x, u) }
  }
  x
}

build_input_inventory <- function(cfg) {
  yrs <- as.integer(unlist(cfg$data$years))
  if (isTRUE(cfg$data$mask_to_regions) && (!nzchar(cfg$data$regions_shapefile) || !file.exists(cfg$data$regions_shapefile))) stop("mask_to_regions=true but regions_shapefile is missing: ", cfg$data$regions_shapefile, call. = FALSE)
  tabs <- lapply(names(cfg$data$responses), function(nm) discover_response_files(cfg$data$responses[[nm]], yrs)); names(tabs) <- names(cfg$data$responses)
  rows <- do.call(rbind, lapply(names(tabs), function(nm) {
    z <- file_manifest(tabs[[nm]]$file, role = paste0(nm, "_", tabs[[nm]]$year)); z$response <- nm; z$year <- tabs[[nm]]$year; z
  }))
  if (file.exists(cfg$data$regions_shapefile)) {
    stem <- tools::file_path_sans_ext(cfg$data$regions_shapefile)
    side <- paste0(stem, c(".shp", ".shx", ".dbf", ".prj", ".cpg", ".qpj")); side <- side[file.exists(side)]
    if (length(side)) rows <- rbind(rows, transform(file_manifest(side, role = paste0("regions_", basename(side))), response = "regions", year = NA_integer_))
  }
  fp <- substr(hash_text(c(paste0("input_version=", cfg$data$input_version %||% ""), apply(rows[, c("role", "path", "size", "md5")], 1L, paste, collapse = "|"))), 1L, 12L)
  list(files = tabs, manifest = rows, fingerprint = fp)
}

write_resolved_config <- function(cfg, ctx) {
  assert_packages("yaml")
  clean <- cfg[!grepl("^\\.", names(cfg))]
  yaml::write_yaml(clean, file.path(ctx$root, "config", "config_resolved.yml"))
  write_csv(cfg$.code_manifest, file.path(ctx$root, "config", "code_manifest.csv"))
  invisible(TRUE)
}
