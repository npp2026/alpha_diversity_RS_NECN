run_region_metadata_tests_v56 <- function(project_dir = getwd()) {
  project_dir <- normalizePath(project_dir, winslash = "/", mustWork = TRUE)
  source(file.path(project_dir, "R", "load_all.R"), local = .GlobalEnv)
  load_v56_modules(project_dir)
  assert_packages("terra")

  v <- terra::vect(data.frame(
    x = c(0, 2, 4, 6, 8), y = c(0, 0, 0, 0, 0),
    L2_code = c("III4", "I1", "I2", "I2", "I2"),
    L2_name = c("DXAL", "XXAL", "CBS", "CBS", "CBS"),
    Reg_CN = c("A", "B", "C", "D", "E"),
    Reg_EN = c("A_en", "B_en", "C_en", "D_en", "E_en")
  ), geom = c("x", "y"), crs = "EPSG:3857")
  cfg <- list(data = list(region_id_field = "Reg_EN", region_label_field = "Reg_CN",
                          region_group_field = "L2_code", region_group_label_field = "L2_name"))
  m <- region_metadata(v, cfg, require_unique_id = TRUE)
  stopifnot(nrow(m) == 5L, length(unique(m$region)) == 5L, length(unique(m$region_group)) == 3L,
            identical(m$region, c("A_en", "B_en", "C_en", "D_en", "E_en")))
  cand <- region_unique_character_candidates(v)
  stopifnot(isTRUE(cand$unique_character_candidate[cand$field == "Reg_EN"]),
            !isTRUE(cand$unique_character_candidate[cand$field == "L2_code"]))

  cfg_bad <- cfg; cfg_bad$data$region_id_field <- "L2_code"
  err <- try(region_metadata(v, cfg_bad, require_unique_id = TRUE), silent = TRUE)
  stopifnot(inherits(err, "try-error"))
  cat("Region metadata tests passed (feature ID vs mountain-system grouping).\n")
  data.frame(check = c("unique_feature_id", "nonunique_group_rejected_as_id", "group_metadata_preserved"), pass = TRUE)
}
