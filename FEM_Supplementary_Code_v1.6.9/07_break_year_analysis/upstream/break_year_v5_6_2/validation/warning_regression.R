run_warning_regression_v56 <- function(project_dir = getwd()) {
  project_dir <- normalizePath(project_dir, winslash = "/", mustWork = TRUE)
  source(file.path(project_dir, "R", "load_all.R"), local = .GlobalEnv)
  load_v56_modules(project_dir)
  assert_packages("terra")

  captured <- character()
  capture_target_warnings <- function(expr) {
    withCallingHandlers(
      expr,
      warning = function(w) {
        msg <- conditionMessage(w)
        if (grepl("\\[readStart\\] source already open for reading|uneven horizontal intervals", msg)) {
          captured <<- c(captured, msg)
        }
        invokeRestart("muffleWarning")
      }
    )
  }

  td <- tempfile("v56_warning_regression_")
  dir.create(td, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(td, recursive = TRUE, force = TRUE), add = TRUE)

  r <- terra::rast(nrows = 4, ncols = 5, xmin = 0, xmax = 10, ymin = 0, ymax = 8,
                   crs = "EPSG:3857")
  terra::values(r) <- c(0, 1, NA, 2, 3,
                        1, NA, 2, 3, 4,
                        0, 1, 2, NA, 3,
                        1, 2, 3, 4, NA)
  by <- terra::rast(r); terra::values(by) <- c(2005, 2006, NA, 2008, 2009,
                                               2005, NA, 2007, 2008, 2009,
                                               2006, 2007, 2008, NA, 2010,
                                               2007, 2008, 2009, 2010, NA)
  mask <- terra::ifel(r >= 2, 1, 0)

  capture_target_warnings({
    invisible(count_mc_exceed_raster(r, B = 4L))
    invisible(map_fdr_lookup_raster(r, lookup = seq(0.1, 0.5, length.out = 5L),
                                    filename = file.path(td, "q.tif"), overwrite = TRUE))
    invisible(recovery_break_counts(by, mask))
    invisible(cells_from_mask(mask))
    invisible(binary_iou_rasters(mask, mask))
    invisible(break_year_agreement_rasters(mask, by, mask, by))
  })

  pass <- length(captured) == 0L
  if (!pass) stop("Target warning regression failed: ", paste(unique(captured), collapse = " | "), call. = FALSE)
  cat("Warning regression tests passed (terra block reads).\n")
  data.frame(check = "terra_read_connection_warning", pass = TRUE)
}
