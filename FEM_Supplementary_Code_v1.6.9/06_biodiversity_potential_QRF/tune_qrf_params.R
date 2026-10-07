# ============================================================
# tune_qrf_params.R
#
# Purpose
# -------
# Tune ranger Quantile Regression Forest (QRF) parameters FIRST,
# WITHOUT producing final GeoTIFF prediction maps.
#
# Workflow
#   1) Read ENV.tif predictors: DEM + BIO6/BIO10/BIO12/BIO17/BIO15/BIO4
#   2) Read train4pot.csv WGS84 long/lat, project points to ENV CRS
#   3) Extract ENV predictors at training points
#   4) Calibrated fill of ENV-extraction NAs from the CSV's own
#      bio{x}_wc / DEMc columns (DEM + climate only; never Forest_age)
#   5) Build response-specific training data for Rich_tree / Shannon_wiener
#   6) Tune the QRF parameter grid with repeated hold-out
#      (validation splits are fixed per repeat and shared across all
#       candidates, so score differences reflect parameters, not split luck)
#   7) Run kNNDM spatial CV only for the top hold-out candidates
#   8) Save ranked parameter tables and a ready-to-source snippet
#      (apply_best_qrf_params_to_final_mapping.R) for predict_age_sensitivity.R
#
# Outputs are CSV/RDS only. No quantile surfaces or sensitivity GeoTIFFs are written.
#
# Pairs with: predict_age_sensitivity.R (final validation + mapping).
# ============================================================

suppressPackageStartupMessages({
  library(terra)
})

# ============================================================
# 0. CONFIG
# ============================================================

# Portable paths: keep input/output paths relative to an explicit data root.
# Shared input/cache/output contracts are resolved relative to this installed script.
.fem_script<-sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1])
source(file.path(dirname(normalizePath(.fem_script)),"..","R","qrf_contracts.R"))
rm(.fem_script)

base_dir <- Sys.getenv("FEM_POTENTIAL_DATA_DIR", unset = getwd())
base_dir <- normalizePath(path.expand(base_dir), winslash = "/", mustWork = TRUE)
setwd(base_dir)

env_file <- "ENV.tif"
csv_file <- "train4pot.csv"
template_file <- "templ_1km.tif"
# TRUE: templ_1km.tif 中 0 与 NA 都视为模板外；若 0 是有效像元值，请改为 FALSE。
mask_zero_as_na <- TRUE

out_dir <- "qrf_parameter_tuning_outputs"
fem_qrf_new_output(out_dir)

# ENV.tif 只读；DEM/BIO 会派生为共享的 templ_1km 对齐+mask 缓存，供调参/预测/oldage 复用。
aligned_env_file <- file.path(base_dir, "ENV_predictors_aligned_to_templ_1km_masked.tif")
reuse_aligned_env_cache <- TRUE

terra_tempdir <- file.path(base_dir, "terra_tmp")
dir.create(terra_tempdir, showWarnings = FALSE, recursive = TRUE)
terra::terraOptions(tempdir = terra_tempdir, memfrac = 0.70, progress = 0)
try(terra::terraOptions(memmax = 24), silent = TRUE)

response_vars <- c("Rich_tree", "Shannon_wiener")
climate_vars <- c("BIO6", "BIO10", "BIO12", "BIO17", "BIO15", "BIO4")
raster_predictor_vars <- c("DEM", climate_vars)
predictor_vars <- c(raster_predictor_vars, "Forest_age")

quantiles <- c(0.90, 0.95)
quantiles <- sort(unique(quantiles))
q_names <- paste0("q", round(quantiles * 100))

nonnegative_responses <- c("Rich_tree", "Shannon_wiener")
clip_negative_predictions <- TRUE

seed <- 123
set.seed(seed)

t_start <- Sys.time()

# Threads
hardware_cores <- fem_qrf_threads(8L)
reserve_cores_for_os <- 0
detected_cores <- parallel::detectCores(logical = TRUE)
if (is.na(detected_cores)) detected_cores <- hardware_cores
n_cores <- max(1, min(hardware_cores, detected_cores) - reserve_cores_for_os)
Sys.setenv(
  OMP_NUM_THREADS = as.character(n_cores),
  MKL_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  NUMEXPR_NUM_THREADS = "1"
)

# Packages
if (!requireNamespace("ranger", quietly = TRUE)) {
  stop("需要安装 ranger 包：install.packages('ranger')")
}

# -------------------------
# Tuning design
# -------------------------
# Stage 1: repeated hold-out over all candidates.
run_holdout_tuning <- TRUE
validation_fraction <- 0.20
validation_repeats <- 5
qrf_num_trees_tuning <- 500
# Final mapping tree-count candidates for optional stability checks in the production script.
# Do NOT multiply the full mtry/node/age tuning grid by these values unless you have ample time.
# More trees mainly stabilize q90/q95 surfaces; they usually do not need selection like mtry/node size.
qrf_num_trees_final_candidates <- c(800, 1200, 1500, 2000)
qrf_num_trees_validation_default <- 800
qrf_num_trees_final_default <- 1500

# Stage 2: kNNDM only for top hold-out candidates.
# This avoids running kNNDM for every candidate, which can be very slow.
run_knndm_tuning <- TRUE
top_n_for_knndm <- 5
knndm_k <- 10
knndm_samplesize <- 10000
knndm_clustering <- "kmeans"
knndm_maxp <- 0.5
knndm_dist_space <- "geographical"

# Scoring weights. Lower score is better.
# coverage error is the primary objective for quantile models.
score_weight_coverage_abs <- 1.00
score_weight_undercoverage <- 0.50
score_weight_pinball_scaled <- 0.10

# Parameter grid. Keep it small first.
tuning_grid <- expand.grid(
  qrf_mtry = c(2, 3, 4),
  qrf_min_node_size = c(3, 5, 10),
  qrf_always_split_age = c(TRUE, FALSE),
  stringsAsFactors = FALSE
)
tuning_grid$param_id <- sprintf(
  "mtry%s_node%s_age%s",
  tuning_grid$qrf_mtry,
  tuning_grid$qrf_min_node_size,
  ifelse(tuning_grid$qrf_always_split_age, "Y", "N")
)
tuning_grid <- tuning_grid[, c("param_id", "qrf_mtry", "qrf_min_node_size", "qrf_always_split_age")]

write.csv(tuning_grid, file.path(out_dir, "qrf_tuning_parameter_grid.csv"), row.names = FALSE)

# ============================================================
# 1. Utility functions
# ============================================================

safe_as_numeric <- function(x) fem_qrf_numeric(x)


safe_sample <- function(x, size, replace = FALSE, prob = NULL) {
  if (length(x) == 1L) {
    if (isTRUE(replace)) return(rep(x, size))
    if (size > 1L) {
      stop("safe_sample: cannot sample size > 1 from a length-1 vector without replacement")
    }
    return(x)
  }
  x[sample.int(length(x), size = size, replace = replace, prob = prob)]
}


stop_if_missing <- function(required, available, object_name) {
  miss <- setdiff(required, available)
  if (length(miss) > 0) {
    stop(
      object_name, " 缺少以下字段/图层：", paste(miss, collapse = ", "),
      "\n当前可用名称：", paste(available, collapse = ", ")
    )
  }
}

pinball_loss <- function(y, pred, q) {
  err <- y - pred
  mean(ifelse(err >= 0, q * err, (q - 1) * err), na.rm = TRUE)
}

enforce_quantile_order_df <- function(p, quantiles) {
  ord <- order(quantiles)
  p <- p[, ord, drop = FALSE]
  names(p) <- paste0("q", round(sort(quantiles) * 100))
  if (ncol(p) <= 1) return(p)
  for (i in 2:ncol(p)) {
    p[[i]] <- pmax(p[[i]], p[[i - 1]], na.rm = FALSE)
  }
  p
}

clip_negative_df <- function(p) {
  for (nm in names(p)) p[[nm]] <- pmax(p[[nm]], 0, na.rm = FALSE)
  p
}

postprocess_prediction_df <- function(p, response_var) {
  p <- enforce_quantile_order_df(p, quantiles)
  if (clip_negative_predictions && response_var %in% nonnegative_responses) {
    p <- clip_negative_df(p)
  }
  p
}

train_qrf <- function(train_dat,
                      qrf_mtry,
                      qrf_min_node_size,
                      qrf_always_split_age,
                      num_trees = qrf_num_trees_tuning,
                      model_seed = seed) {
  d <- train_dat[, c(predictor_vars, "response"), drop = FALSE]

  args <- list(
    formula = response ~ .,
    data = d,
    num.trees = num_trees,
    mtry = qrf_mtry,
    min.node.size = qrf_min_node_size,
    quantreg = TRUE,
    importance = "none",
    oob.error = FALSE,
    num.threads = n_cores,
    seed = model_seed
  )

  if (isTRUE(qrf_always_split_age)) {
    args$always.split.variables <- "Forest_age"
  }

  do.call(ranger::ranger, args)
}

predict_qrf <- function(model, newdata) {
  newdata <- newdata[, predictor_vars, drop = FALSE]
  p <- predict(
    model,
    data = newdata,
    type = "quantiles",
    quantiles = quantiles,
    num.threads = n_cores
  )$predictions

  if (is.vector(p)) p <- matrix(p, ncol = length(quantiles))
  p <- as.data.frame(p)
  names(p) <- q_names
  p
}

evaluate_prediction <- function(y, pred_df, response_var, param_row,
                                validation_method, repeat_id = NA_integer_,
                                fold_id = NA_integer_, n_train = NA_integer_) {
  rows <- list()
  y_sd <- stats::sd(y, na.rm = TRUE)
  if (!is.finite(y_sd) || y_sd <= 0) y_sd <- 1

  for (i in seq_along(quantiles)) {
    q <- quantiles[i]
    qn <- q_names[i]
    pred <- pred_df[[qn]]
    cov <- mean(y <= pred, na.rm = TRUE)
    pb <- pinball_loss(y, pred, q)
    rows[[i]] <- data.frame(
      response = response_var,
      validation_method = validation_method,
      param_id = param_row$param_id,
      qrf_mtry = param_row$qrf_mtry,
      qrf_min_node_size = param_row$qrf_min_node_size,
      qrf_always_split_age = param_row$qrf_always_split_age,
      repeat_id = repeat_id,
      fold_id = fold_id,
      quantile = q,
      n_train = n_train,
      n_validation = length(y),
      pinball_loss = pb,
      pinball_loss_scaled = pb / y_sd,
      empirical_coverage = cov,
      coverage_error = cov - q,
      abs_coverage_error = abs(cov - q),
      undercoverage = max(0, q - cov),
      prediction_min = min(pred, na.rm = TRUE),
      prediction_max = max(pred, na.rm = TRUE)
    )
  }
  do.call(rbind, rows)
}

summarize_and_score <- function(metrics_df, method_label) {
  if (nrow(metrics_df) == 0) return(data.frame())

  agg <- aggregate(
    cbind(pinball_loss, pinball_loss_scaled, empirical_coverage,
          abs_coverage_error, undercoverage) ~
      response + param_id + qrf_mtry + qrf_min_node_size + qrf_always_split_age + quantile,
    data = metrics_df,
    FUN = mean,
    na.rm = TRUE
  )

  score <- aggregate(
    cbind(abs_coverage_error, undercoverage, pinball_loss_scaled) ~
      response + param_id + qrf_mtry + qrf_min_node_size + qrf_always_split_age,
    data = agg,
    FUN = mean,
    na.rm = TRUE
  )

  score$score <-
    score_weight_coverage_abs * score$abs_coverage_error +
    score_weight_undercoverage * score$undercoverage +
    score_weight_pinball_scaled * score$pinball_loss_scaled

  score$validation_method <- method_label
  score <- score[order(score$response, score$score), ]

  list(by_quantile = agg, score = score)
}

# ============================================================
# templ_1km 强制对齐工具
#   - 原始 ENV.tif 只读，不覆盖、不改写。
#   - 生成派生缓存 ENV_predictors_aligned_to_templ_1km_masked.tif。
#   - 注意：调参脚本不使用 templ_1km 过滤训练点；该工具仅保留给正式预测/AOA脚本复用。
# ============================================================
read_templ_1km_mask <- function(template_file) {
  if (!file.exists(template_file)) {
    stop("缺少 templ_1km.tif：", template_file,
         "\n请把 templ_1km.tif 放在 base_dir 下，或修改 template_file。")
  }
  templ <- terra::rast(template_file)
  if (terra::nlyr(templ) > 1) templ <- templ[[1]]
  names(templ) <- "templ_1km_mask_raw"

  # 关键修正：terra::mask() 默认只按 NA 遮罩；若模板无效区用 0 表示，必须先把 0 转成 NA。
  # 如果你的 templ_1km.tif 中 0 是有效值，可在 CONFIG 中设置 mask_zero_as_na <- FALSE。
  if (!exists("mask_zero_as_na", inherits = TRUE)) mask_zero_as_na <- TRUE
  if (isTRUE(mask_zero_as_na)) {
    templ <- terra::ifel(!is.finite(templ) | templ == 0, NA, 1)
  } else {
    templ <- terra::ifel(!is.finite(templ), NA, 1)
  }
  names(templ) <- "templ_1km_mask"
  templ
}

file_fingerprint <- function(paths) fem_qrf_fingerprint(paths,mask_zero_as_na)

cache_metadata_matches <- function(meta_file,source_files) fem_qrf_cache_matches(meta_file,source_files,mask_zero_as_na)

write_cache_metadata <- function(meta_file, source_files) {
  if (is.null(source_files)) return(invisible(FALSE))
  utils::write.csv(file_fingerprint(source_files), meta_file, row.names = FALSE)
  invisible(TRUE)
}

align_env_to_templ_1km_mask <- function(env_pred, templ,
                                        aligned_file = NULL,
                                        reuse_cache = TRUE,
                                        source_files = NULL) {
  meta_file <- if (!is.null(aligned_file)) paste0(aligned_file, ".meta.csv") else NULL

  if (!is.null(aligned_file) && isTRUE(reuse_cache) && file.exists(aligned_file)) {
    metadata_ok <- if (is.null(source_files)) TRUE else cache_metadata_matches(meta_file, source_files)
    cached <- terra::rast(aligned_file)
    if (metadata_ok && all(names(env_pred) %in% names(cached)) &&
        isTRUE(terra::compareGeom(cached[[1]], templ, stopOnError = FALSE,
                                  crs = TRUE, ext = TRUE, rowcol = TRUE, res = TRUE))) {
      message("Reusing templ_1km-aligned masked ENV cache: ", aligned_file)
      cached <- cached[[names(env_pred)]]
      cached <- terra::mask(cached, templ)
      names(cached) <- names(env_pred)
      return(cached)
    }
    message("Existing aligned ENV cache is stale, missing layers, or not aligned to templ_1km; rebuilding.")
  }

  if (isTRUE(terra::compareGeom(env_pred[[1]], templ, stopOnError = FALSE,
                                crs = TRUE, ext = TRUE, rowcol = TRUE, res = TRUE))) {
    message("ENV predictors already match templ_1km geometry; applying templ_1km mask only.")
    env_aligned <- env_pred
  } else if (!isTRUE(terra::same.crs(env_pred, templ))) {
    message("ENV.tif CRS differs from templ_1km.tif; projecting ENV predictors to templ_1km grid with bilinear method.")
    env_aligned <- terra::project(env_pred, templ, method = "bilinear")
  } else {
    message("Resampling ENV predictors to templ_1km grid with bilinear method.")
    env_aligned <- terra::resample(env_pred, templ, method = "bilinear")
  }

  names(env_aligned) <- names(env_pred)
  env_aligned <- terra::mask(env_aligned, templ)
  names(env_aligned) <- names(env_pred)

  if (!isTRUE(terra::compareGeom(env_aligned[[1]], templ, stopOnError = FALSE,
                                 crs = TRUE, ext = TRUE, rowcol = TRUE, res = TRUE))) {
    stop("ENV predictors failed to align to templ_1km.tif after project/resample.")
  }

  if (!is.null(aligned_file)) {
    terra::writeRaster(
      env_aligned,
      aligned_file,
      overwrite = TRUE,
      datatype = "FLT4S",
      gdal = c("COMPRESS=LZW", "TILED=YES", "BIGTIFF=IF_SAFER", "NUM_THREADS=ALL_CPUS", "PREDICTOR=3")
    )
    write_cache_metadata(meta_file, source_files)
    message("Wrote templ_1km-aligned masked ENV predictors: ", aligned_file)
    env_aligned <- terra::rast(aligned_file)[[names(env_pred)]]
  }

  env_aligned
}

align_raster_to_templ_1km_mask <- function(r, templ, method = "near") {
  if (!isTRUE(terra::same.crs(r, templ))) {
    r2 <- terra::project(r, templ, method = method)
  } else if (!isTRUE(terra::compareGeom(r[[1]], templ, stopOnError = FALSE,
                                        crs = TRUE, ext = TRUE, rowcol = TRUE, res = TRUE))) {
    r2 <- terra::resample(r, templ, method = method)
  } else {
    r2 <- r
  }
  r2 <- terra::mask(r2, templ)
  names(r2) <- names(r)
  r2
}


# ============================================================
# 2. Read ENV and CSV; extract predictors
# ============================================================

message("==================================================")
message("Reading ENV and CSV...")

ENV <- rast(env_file)
message("Raw ENV raster is read-only and will NOT be overwritten.")
stop_if_missing(raster_predictor_vars, names(ENV), "ENV.tif")

# Training extraction retains raw-ENV support; kNNDM prediction points use
# the same aligned template domain as production/AOA.
ENV_model_pred <- ENV[[raster_predictor_vars]]
ENV_model_pred <- terra::ifel(is.finite(ENV_model_pred)&ENV_model_pred!=-9999,ENV_model_pred,NA)
template_1km <- read_templ_1km_mask(template_file)
ENV_pred <- align_env_to_templ_1km_mask(ENV_model_pred,template_1km,
  aligned_file=aligned_env_file,reuse_cache=reuse_aligned_env_cache,
  source_files=c(env_file,template_file))

na_env <- as.data.frame(global(is.na(ENV_pred), "sum", na.rm = TRUE))
na_env$layer <- rownames(na_env)
rownames(na_env) <- NULL
names(na_env)[1] <- "NA_count"
write.csv(na_env[, c("layer", "NA_count")], file.path(out_dir, "ENV_predictor_NA_summary.csv"), row.names = FALSE)

dat_raw <- read.csv(csv_file, stringsAsFactors = FALSE)
required_cols <- c("long", "lat", "Forest_age", response_vars)
stop_if_missing(required_cols, names(dat_raw), "CSV")

for (cc in required_cols) dat_raw[[cc]] <- safe_as_numeric(dat_raw[[cc]])

coord_ok <- !is.na(dat_raw$long) & !is.na(dat_raw$lat) &
  is.finite(dat_raw$long) & is.finite(dat_raw$lat)
dat0 <- dat_raw[coord_ok, , drop = FALSE]

if (nrow(dat0) == 0) stop("坐标有效样点数为 0，请检查 long/lat。")
if (any(dat0$long < -180 | dat0$long > 180 | dat0$lat < -90 | dat0$lat > 90, na.rm = TRUE)) {
  warning("发现超出常规经纬度范围的 long/lat，请确认 CSV 坐标确实为 WGS84 经纬度。")
}

pts_ll <- vect(dat0, geom = c("long", "lat"), crs = "EPSG:4326")
pts_env <- project(pts_ll, crs(ENV_pred))
# 不按 templ_1km mask 剔除训练/调参样点；模板外样点仍可用于学习环境-林龄-响应关系。
xy_env <- crds(pts_env)
dat0$x_env <- xy_env[, 1]
dat0$y_env <- xy_env[, 2]

env_extract <- extract(ENV_model_pred,project(pts_ll,crs(ENV_model_pred)),ID=FALSE)
keep_cols <- setdiff(names(dat0), raster_predictor_vars)
dat0 <- dat0[, keep_cols, drop = FALSE]
dat0 <- cbind(dat0, env_extract)

# ============================================================
# 5b. 用 CSV 中的 bio{x}_wc / DEMc 对 ENV 提取 NA 做“订正后填补”
#     只填补 raster_predictor_vars 中的 DEM + BIO 气候变量，不填补 Forest_age。
#
#     对每个变量单独建立订正模型：
#       ENV变量 ~ CSV对应变量 + long + lat
#     其中订正模型只使用 ENV 非 NA 且 CSV 对应变量非 NA 的重叠样点。
#     然后仅对 ENV 提取为 NA、但 CSV 对应变量非 NA 的样点预测填补。
# ============================================================

fill_env_na_from_csv_calibrated <- function(dat,
                                            target_vars,
                                            min_pairs = 30,
                                            min_r2_for_model = 0.70,
                                            clamp_to_observed_range = TRUE,
                                            direct_fill_if_model_poor = TRUE) {
  fallback_map <- setNames(rep(NA_character_, length(target_vars)), target_vars)

  if ("DEM" %in% target_vars) {
    fallback_map["DEM"] <- "DEMc"
  }

  bio_targets <- grep("^BIO[0-9]+$", target_vars, value = TRUE)
  fallback_map[bio_targets] <- paste0("bio", sub("^BIO", "", bio_targets), "_wc")
  fallback_map <- fallback_map[!is.na(fallback_map)]

  missing_fallback <- setdiff(unname(fallback_map), names(dat))
  if (length(missing_fallback) > 0) {
    stop(
      "CSV 缺少用于订正填补 ENV NA 的字段：",
      paste(missing_fallback, collapse = ", "),
      "\n请确认 csv_file 指向的是包含 bio{x}_wc 和 DEMc 的样表。"
    )
  }

  if (!all(c("long", "lat") %in% names(dat))) {
    stop("订正填补需要 CSV 中存在 long 和 lat 字段。")
  }

  dat$long <- safe_as_numeric(dat$long)
  dat$lat  <- safe_as_numeric(dat$lat)

  fill_summary <- data.frame()

  for (target in names(fallback_map)) {
    source <- fallback_map[[target]]

    dat[[target]] <- safe_as_numeric(dat[[target]])
    dat[[source]] <- safe_as_numeric(dat[[source]])

    target_bad_before <- is.na(dat[[target]]) | !is.finite(dat[[target]])
    target_ok_before  <- !target_bad_before
    source_ok <- !is.na(dat[[source]]) & is.finite(dat[[source]])
    xy_ok <- !is.na(dat$long) & is.finite(dat$long) &
      !is.na(dat$lat) & is.finite(dat$lat)

    # 订正模型训练数据：ENV 提取值和 CSV 原始值均非 NA 的重叠样点
    pair_idx <- target_ok_before & source_ok & xy_ok

    # 待填补数据：ENV 提取值 NA，但 CSV 原始值非 NA
    fill_idx <- target_bad_before & source_ok & xy_ok

    n_pair <- sum(pair_idx)
    n_to_fill <- sum(fill_idx)
    n_filled_actual <- 0L
    method_used <- "none_needed_or_no_csv_value"
    r2 <- NA_real_
    rmse <- NA_real_
    fit_status <- "not_run"

    if (n_to_fill > 0) {
      if (n_pair >= min_pairs) {
        cal_dat <- data.frame(
          y = dat[[target]][pair_idx],
          x = dat[[source]][pair_idx],
          long = dat$long[pair_idx],
          lat = dat$lat[pair_idx]
        )

        pred_dat <- data.frame(
          x = dat[[source]][fill_idx],
          long = dat$long[fill_idx],
          lat = dat$lat[fill_idx]
        )

        fit <- try(stats::lm(y ~ x + long + lat, data = cal_dat), silent = TRUE)

        if (!inherits(fit, "try-error")) {
          fit_status <- "ok"
          fit_pred <- try(stats::predict(fit, newdata = cal_dat), silent = TRUE)

          if (!inherits(fit_pred, "try-error")) {
            if (stats::sd(cal_dat$y, na.rm = TRUE) > 0 &&
                stats::sd(fit_pred, na.rm = TRUE) > 0) {
              r2 <- suppressWarnings(stats::cor(cal_dat$y, fit_pred, use = "complete.obs")^2)
            }
            rmse <- sqrt(mean((cal_dat$y - fit_pred)^2, na.rm = TRUE))
          }

          if (is.finite(r2) && r2 >= min_r2_for_model) {
            fill_values <- try(stats::predict(fit, newdata = pred_dat), silent = TRUE)

            if (!inherits(fill_values, "try-error")) {
              fill_values <- as.numeric(fill_values)

              if (clamp_to_observed_range) {
                obs_min <- min(cal_dat$y, na.rm = TRUE)
                obs_max <- max(cal_dat$y, na.rm = TRUE)
                fill_values <- pmax(pmin(fill_values, obs_max), obs_min)
              }

              valid_fill <- is.finite(fill_values)
              idx_all <- which(fill_idx)
              dat[[target]][idx_all[valid_fill]] <- fill_values[valid_fill]
              n_filled_actual <- sum(valid_fill)
              method_used <- "lm_calibrated_ENV_on_CSV_long_lat"
            } else if (direct_fill_if_model_poor) {
              dat[[target]][fill_idx] <- dat[[source]][fill_idx]
              n_filled_actual <- n_to_fill
              method_used <- "direct_fill_predict_failed"
              fit_status <- "predict_failed"
            } else {
              method_used <- "not_filled_predict_failed"
              fit_status <- "predict_failed"
            }
          } else if (direct_fill_if_model_poor) {
            dat[[target]][fill_idx] <- dat[[source]][fill_idx]
            n_filled_actual <- n_to_fill
            method_used <- "direct_fill_low_calibration_R2"
          } else {
            method_used <- "not_filled_low_calibration_R2"
          }
        } else if (direct_fill_if_model_poor) {
          dat[[target]][fill_idx] <- dat[[source]][fill_idx]
          n_filled_actual <- n_to_fill
          method_used <- "direct_fill_lm_failed"
          fit_status <- "lm_failed"
        } else {
          method_used <- "not_filled_lm_failed"
          fit_status <- "lm_failed"
        }
      } else if (direct_fill_if_model_poor) {
        dat[[target]][fill_idx] <- dat[[source]][fill_idx]
        n_filled_actual <- n_to_fill
        method_used <- "direct_fill_too_few_calibration_pairs"
        fit_status <- "too_few_pairs"
      } else {
        method_used <- "not_filled_too_few_calibration_pairs"
        fit_status <- "too_few_pairs"
      }
    }

    remaining_bad <- is.na(dat[[target]]) | !is.finite(dat[[target]])

    fill_summary <- rbind(
      fill_summary,
      data.frame(
        variable = target,
        fallback_variable = source,
        n_calibration_pairs = n_pair,
        extracted_NA_before_fill = sum(target_bad_before),
        candidate_NA_with_csv_value = n_to_fill,
        filled_NA = n_filled_actual,
        remaining_NA_after_fill = sum(remaining_bad),
        calibration_R2 = r2,
        calibration_RMSE = rmse,
        fit_status = fit_status,
        method = method_used,
        stringsAsFactors = FALSE
      )
    )
  }

  return(list(data = dat, summary = fill_summary))
}

fill_res <- fill_env_na_from_csv_calibrated(
  dat = dat0,
  target_vars = raster_predictor_vars,
  min_pairs = 30,
  min_r2_for_model = 0.70,
  clamp_to_observed_range = TRUE,
  direct_fill_if_model_poor = TRUE
)

dat0 <- fill_res$data
fill_summary <- fill_res$summary

write.csv(
  fill_summary,
  file.path(out_dir, "ENV_predictor_NA_filled_from_csv_calibrated_summary.csv"),
  row.names = FALSE
)

print(fill_summary)

write.csv(
  dat0,
  file.path(out_dir, "training_points_with_ENV_extracted_calibrated_filled.csv"),
  row.names = FALSE
)

extract_na_summary <- data.frame(
  variable = raster_predictor_vars,
  extracted_NA = sapply(raster_predictor_vars, function(v) sum(is.na(dat0[[v]]))),
  extracted_nonNA = sapply(raster_predictor_vars, function(v) sum(!is.na(dat0[[v]])))
)
write.csv(extract_na_summary, file.path(out_dir, "extracted_ENV_predictor_NA_summary.csv"), row.names = FALSE)

make_training_data <- function(dat, response_var) {
  use_cols <- c(predictor_vars, "x_env", "y_env", response_var)
  d <- dat[, use_cols, drop = FALSE]
  for (cc in names(d)) d[[cc]] <- safe_as_numeric(d[[cc]])
  ok <- complete.cases(d)
  for (cc in names(d)) ok <- ok & is.finite(d[[cc]])
  d <- d[ok, , drop = FALSE]
  names(d)[names(d) == response_var] <- "response"
  d
}

train_list <- list()
summary_rows <- list()
for (resp in response_vars) {
  d <- make_training_data(dat0, resp)
  message("Training rows for ", resp, ": ", nrow(d))
  if (nrow(d) < 30) stop("有效训练样本过少：", resp, "; n = ", nrow(d))
  train_list[[resp]] <- d

  row <- data.frame(
    response = resp,
    n_train = nrow(d),
    Forest_age_min = min(d$Forest_age, na.rm = TRUE),
    Forest_age_max = max(d$Forest_age, na.rm = TRUE)
  )
  for (v in raster_predictor_vars) {
    row[[paste0(v, "_min")]] <- min(d[[v]], na.rm = TRUE)
    row[[paste0(v, "_max")]] <- max(d[[v]], na.rm = TRUE)
  }
  summary_rows[[resp]] <- row
}
training_summary <- do.call(rbind, summary_rows)
write.csv(training_summary, file.path(out_dir, "model_training_summary.csv"), row.names = FALSE)

# ============================================================
# 3. Stage 1: repeated hold-out tuning over full grid
# ============================================================

holdout_metrics <- data.frame()

if (run_holdout_tuning) {
  message("==================================================")
  message("Stage 1: repeated hold-out tuning over ", nrow(tuning_grid), " parameter sets")

  holdout_rows <- list()

  for (resp in response_vars) {
    d <- train_list[[resp]]
    n <- nrow(d)
    val_n <- max(1, floor(n * validation_fraction))

    # Fix the validation splits per repeat ONCE per response and reuse them for
    # EVERY parameter candidate. This makes the comparison paired: differences in
    # the resulting score reflect the parameters, not which rows happened to land
    # in the validation set. (Previously the split depended on the grid index g,
    # so each candidate was scored on a different split, adding avoidable noise.)
    val_idx_list <- lapply(seq_len(validation_repeats), function(rep_i) {
      set.seed(seed + rep_i)
      safe_sample(seq_len(n), size = val_n)
    })

    for (g in seq_len(nrow(tuning_grid))) {
      param <- tuning_grid[g, ]
      message("Hold-out: ", resp, " | ", param$param_id, " (", g, "/", nrow(tuning_grid), ")")

      for (rep_i in seq_len(validation_repeats)) {
        val_idx <- val_idx_list[[rep_i]]
        d_train <- d[-val_idx, , drop = FALSE]
        d_val <- d[val_idx, , drop = FALSE]

        m <- train_qrf(
          d_train,
          qrf_mtry = param$qrf_mtry,
          qrf_min_node_size = param$qrf_min_node_size,
          qrf_always_split_age = param$qrf_always_split_age,
          num_trees = qrf_num_trees_tuning,
          model_seed = seed + 1000 * g + rep_i
        )
        p <- predict_qrf(m, d_val)
        p <- postprocess_prediction_df(p, resp)

        holdout_rows[[length(holdout_rows) + 1L]] <- evaluate_prediction(
          y = d_val$response,
          pred_df = p,
          response_var = resp,
          param_row = param,
          validation_method = "repeated_holdout",
          repeat_id = rep_i,
          n_train = nrow(d_train)
        )
      }
      gc()
    }
  }

  holdout_metrics <- if (length(holdout_rows)) do.call(rbind, holdout_rows) else data.frame()

  write.csv(holdout_metrics, file.path(out_dir, "qrf_tuning_holdout_metrics.csv"), row.names = FALSE)

  hsum <- summarize_and_score(holdout_metrics, "repeated_holdout")
  write.csv(hsum$by_quantile, file.path(out_dir, "qrf_tuning_holdout_summary_by_quantile.csv"), row.names = FALSE)
  write.csv(hsum$score, file.path(out_dir, "qrf_tuning_holdout_ranked_parameters.csv"), row.names = FALSE)
}

# ============================================================
# 4. Stage 2: kNNDM tuning for top hold-out candidates
# ============================================================

make_prediction_points_sf <- function(env_pred,sample_size,seed=123) {
  if(!requireNamespace("sf",quietly=TRUE))stop("Install sf")
  set.seed(seed)
  vals<-terra::values(env_pred,mat=TRUE)
  cells<-which(rowSums(!is.finite(vals)|vals==-9999)==0L)
  if(!length(cells))stop("No jointly finite prediction-domain cells")
  cells<-safe_sample(cells,min(as.integer(sample_size),length(cells)))
  xy<-terra::xyFromCell(env_pred[[1]],cells)
  sf::st_as_sf(data.frame(cell=cells,x=xy[,1],y=xy[,2]),coords=c("x","y"),crs=sf::st_crs(terra::crs(env_pred)))
}

call_knndm <- function(tpoints, predpoints, k, maxp, clustering, dist_space) {
  args_new <- list(
    tpoints = tpoints,
    predpoints = predpoints,
    dist_space = dist_space,
    k = k,
    maxp = maxp,
    clustering = clustering
  )
  out <- try(do.call(CAST::knndm, args_new), silent = TRUE)
  if (!inherits(out, "try-error")) return(out)

  args_old <- args_new
  args_old$space <- args_old$dist_space
  args_old$dist_space <- NULL
  do.call(CAST::knndm, args_old)
}

make_knndm_folds <- function(train_dat, env_pred, predpoints_sf,
                             k = 10, clustering = "kmeans", maxp = 0.5,
                             dist_space = "geographical", seed = 123) {
  if (!requireNamespace("CAST", quietly = TRUE)) stop("需要安装 CAST 包：install.packages('CAST')")
  if (!requireNamespace("sf", quietly = TRUE)) stop("需要安装 sf 包：install.packages('sf')")
  set.seed(seed)

  tpoints <- sf::st_as_sf(
    train_dat,
    coords = c("x_env", "y_env"),
    crs = sf::st_crs(terra::crs(env_pred))
  )

  folds <- tryCatch(
    call_knndm(tpoints, predpoints_sf, k, maxp, clustering, dist_space),
    error = function(e) {
      if (!identical(clustering, "hierarchical")) {
        message("CAST::knndm with clustering='", clustering, "' failed: ", conditionMessage(e))
        message("Retrying with clustering='hierarchical'...")
        call_knndm(tpoints, predpoints_sf, k, maxp, "hierarchical", dist_space)
      } else {
        stop(e)
      }
    }
  )

  if (is.null(folds$indx_train) || is.null(folds$indx_test)) {
    stop("CAST::knndm 返回对象缺少 indx_train / indx_test。")
  }
  folds
}

validate_fold_indices <- function(folds, n) {
  if (length(folds$indx_train) != length(folds$indx_test)) stop("kNNDM train/test fold 数不一致。")
  for (i in seq_along(folds$indx_test)) {
    tr <- folds$indx_train[[i]]
    te <- folds$indx_test[[i]]
    if (length(tr) == 0 || length(te) == 0) stop("kNNDM 第 ", i, " 折为空。")
    if (any(tr < 1 | tr > n) || any(te < 1 | te > n)) stop("kNNDM 第 ", i, " 折索引越界。")
    if (length(intersect(tr, te)) > 0) stop("kNNDM 第 ", i, " 折 train/test 有重叠。")
  }
  all_test <- unlist(folds$indx_test, use.names = FALSE)
  if (length(all_test) != length(unique(all_test))) stop("kNNDM 测试索引在不同 folds 间重复。")
  if (!setequal(all_test, seq_len(n))) stop("kNNDM 测试索引没有覆盖全部训练样本。")
  invisible(TRUE)
}

knndm_metrics <- data.frame()

if (run_knndm_tuning) {
  if (!requireNamespace("CAST", quietly = TRUE) || !requireNamespace("sf", quietly = TRUE)) {
    warning("CAST 或 sf 未安装，跳过 kNNDM 调参；只使用 hold-out 排名。")
  } else {
    message("==================================================")
    message("Stage 2: kNNDM tuning for top ", top_n_for_knndm, " hold-out candidates per response")

    predpoints_sf <- make_prediction_points_sf(ENV_pred, knndm_samplesize, seed = seed)
    saveRDS(predpoints_sf, file.path(out_dir, "qrf_tuning_knndm_predpoints.rds"))

    # Prefer the hold-out ranking already in memory; otherwise read it from disk.
    # If neither is available (e.g. run_holdout_tuning was FALSE and no prior run
    # produced the file), fail with a clear message instead of an obscure read error.
    if (exists("hsum") && is.list(hsum) && !is.null(hsum$score) && nrow(hsum$score) > 0) {
      holdout_rank <- hsum$score
    } else {
      holdout_rank_file <- file.path(out_dir, "qrf_tuning_holdout_ranked_parameters.csv")
      if (!file.exists(holdout_rank_file)) {
        stop("kNNDM 调参需要 hold-out 排名，但内存与磁盘上都没有。\n",
             "请先运行 Stage 1 (run_holdout_tuning <- TRUE)，或提供文件：\n  ",
             holdout_rank_file)
      }
      holdout_rank <- read.csv(holdout_rank_file, stringsAsFactors = FALSE)
    }

    knndm_rows <- list()

    for (resp in response_vars) {
      d <- train_list[[resp]]
      hr <- holdout_rank[holdout_rank$response == resp, , drop = FALSE]
      hr <- hr[order(hr$score), , drop = FALSE]
      top_ids <- head(hr$param_id, top_n_for_knndm)
      candidate_grid <- tuning_grid[tuning_grid$param_id %in% top_ids, , drop = FALSE]

      message("Making kNNDM folds for ", resp)
      folds <- make_knndm_folds(
        train_dat = d,
        env_pred = ENV_pred,
        predpoints_sf = predpoints_sf,
        k = knndm_k,
        clustering = knndm_clustering,
        maxp = knndm_maxp,
        dist_space = knndm_dist_space,
        seed = seed
      )
      validate_fold_indices(folds, nrow(d))
      saveRDS(folds, file.path(out_dir, paste0("qrf_tuning_knndm_folds_", resp, ".rds")))

      fold_info <- data.frame(
        response = resp,
        k = knndm_k,
        W = if (!is.null(folds$W)) as.numeric(folds$W) else NA_real_,
        n_predpoints = length(predpoints_sf),
        clustering = knndm_clustering,
        maxp = knndm_maxp,
        dist_space = knndm_dist_space
      )
      write.csv(fold_info, file.path(out_dir, paste0("qrf_tuning_knndm_fold_info_", resp, ".csv")), row.names = FALSE)

      for (g in seq_len(nrow(candidate_grid))) {
        param <- candidate_grid[g, ]
        message("kNNDM: ", resp, " | ", param$param_id, " (", g, "/", nrow(candidate_grid), ")")

        pred_rows <- list()
        for (fold_i in seq_along(folds$indx_test)) {
          train_idx <- folds$indx_train[[fold_i]]
          test_idx <- folds$indx_test[[fold_i]]
          d_train <- d[train_idx, , drop = FALSE]
          d_test <- d[test_idx, , drop = FALSE]

          m <- train_qrf(
            d_train,
            qrf_mtry = param$qrf_mtry,
            qrf_min_node_size = param$qrf_min_node_size,
            qrf_always_split_age = param$qrf_always_split_age,
            num_trees = qrf_num_trees_tuning,
            model_seed = seed + 5000 * g + fold_i
          )
          p <- predict_qrf(m, d_test)
          p <- postprocess_prediction_df(p, resp)

          pred_rows[[fold_i]] <- data.frame(row_id = test_idx, y = d_test$response, p)
          gc()
        }

        pred_all <- do.call(rbind, pred_rows)
        pred_all <- pred_all[order(pred_all$row_id), ]

        knndm_rows[[length(knndm_rows) + 1L]] <- evaluate_prediction(
          y = pred_all$y,
          pred_df = pred_all[, q_names, drop = FALSE],
          response_var = resp,
          param_row = param,
          validation_method = "knndm_global",
          n_train = NA_integer_
        )
      }
    }

    knndm_metrics <- if (length(knndm_rows)) do.call(rbind, knndm_rows) else data.frame()

    write.csv(knndm_metrics, file.path(out_dir, "qrf_tuning_knndm_global_metrics.csv"), row.names = FALSE)
    ksum <- summarize_and_score(knndm_metrics, "knndm_global")
    write.csv(ksum$by_quantile, file.path(out_dir, "qrf_tuning_knndm_summary_by_quantile.csv"), row.names = FALSE)
    write.csv(ksum$score, file.path(out_dir, "qrf_tuning_knndm_ranked_parameters.csv"), row.names = FALSE)
  }
}

# ============================================================
# 5. Select best parameters and write final snippets
# ============================================================

select_best <- function() {
  knndm_file <- file.path(out_dir, "qrf_tuning_knndm_ranked_parameters.csv")
  holdout_file <- file.path(out_dir, "qrf_tuning_holdout_ranked_parameters.csv")

  if (file.exists(knndm_file)) {
    ranked <- read.csv(knndm_file, stringsAsFactors = FALSE)
    method <- "knndm_global"
  } else if (file.exists(holdout_file)) {
    ranked <- read.csv(holdout_file, stringsAsFactors = FALSE)
    method <- "repeated_holdout"
  } else {
    stop("没有找到可用于选择最优参数的排名文件。")
  }

  by_resp <- do.call(rbind, lapply(split(ranked, ranked$response), function(z) {
    z[order(z$score), ][1, , drop = FALSE]
  }))
  rownames(by_resp) <- NULL
  by_resp$selection_basis <- method

  # Fair global ranking: only compare parameters that were evaluated for ALL response variables.
  # In two-stage tuning, kNNDM is run only for the top candidates per response; therefore some
  # param_id values may appear for Rich_tree but not for Shannon_wiener, or vice versa.
  # Averaging over an incomplete set would unfairly favor parameters evaluated for only one response.
  n_resp_total <- length(unique(ranked$response))

  global_all <- aggregate(
    score ~ param_id + qrf_mtry + qrf_min_node_size + qrf_always_split_age,
    data = ranked,
    FUN = mean,
    na.rm = TRUE
  )

  coverage_n <- aggregate(
    response ~ param_id + qrf_mtry + qrf_min_node_size + qrf_always_split_age,
    data = ranked,
    FUN = function(z) length(unique(z))
  )
  names(coverage_n)[names(coverage_n) == "response"] <- "n_responses_evaluated"

  global <- merge(
    global_all,
    coverage_n,
    by = c("param_id", "qrf_mtry", "qrf_min_node_size", "qrf_always_split_age"),
    all.x = TRUE
  )
  global$complete_response_coverage <- global$n_responses_evaluated == n_resp_total
  global$selection_basis <- method

  global_fair <- global[global$complete_response_coverage, , drop = FALSE]
  if (nrow(global_fair) == 0) {
    warning("没有任何参数同时覆盖全部响应变量；全局参数退回为不完整聚合结果。建议提高 top_n_for_knndm。")
    global_fair <- global
  }

  global <- global[order(!global$complete_response_coverage, global$score), ]
  global_fair <- global_fair[order(global_fair$score), ]

  list(by_response = by_resp, global = global, global_fair = global_fair)
}

best <- select_best()
write.csv(best$by_response, file.path(out_dir, "qrf_tuning_best_parameters_by_response.csv"), row.names = FALSE)
write.csv(best$global, file.path(out_dir, "qrf_tuning_best_global_parameters.csv"), row.names = FALSE)
write.csv(best$global_fair, file.path(out_dir, "qrf_tuning_best_global_parameters_fair.csv"), row.names = FALSE)

best_global <- best$global_fair[1, ]

# Write a production snippet that lets predict_age_sensitivity use response-specific parameters.
# The old scalar qrf_mtry/qrf_min_node_size values are also written as a fallback.
make_response_param_snippet <- function(best_by_response, best_global_row) {
  lines <- c(
    "# Auto-generated by tune_qrf_params.R",
    "# Source this file from predict_age_sensitivity.R.",
    "# This snippet intentionally writes only structural QRF parameters.",
    "# Tree counts are controlled in the prediction script to avoid silent overrides.",
    "qrf_params_by_response <- list("
  )

  item_lines <- character(0)
  for (i in seq_len(nrow(best_by_response))) {
    row <- best_by_response[i, ]
    resp <- row$response
    always_value <- if (isTRUE(row$qrf_always_split_age)) "\"Forest_age\"" else "NULL"
    comma <- if (i < nrow(best_by_response)) "," else ""
    item_lines <- c(
      item_lines,
      sprintf(
        "  %s = list(qrf_mtry = %s, qrf_min_node_size = %s, qrf_always_split_variables = %s)%s",
        resp, row$qrf_mtry, row$qrf_min_node_size, always_value, comma
      )
    )
  }

  c(
    lines,
    item_lines,
    ")",
    "",
    "# Scalar fallback for older scripts:",
    sprintf("qrf_mtry <- %s", best_global_row$qrf_mtry),
    sprintf("qrf_min_node_size <- %s", best_global_row$qrf_min_node_size),
    sprintf("qrf_always_split_variables <- %s", if (isTRUE(best_global_row$qrf_always_split_age)) "\"Forest_age\"" else "NULL")
  )
}

snippet <- make_response_param_snippet(best$by_response, best_global)
writeLines(snippet, file.path(out_dir, "apply_best_qrf_params_to_final_mapping.R"))

run_info <- list(
  env_file = normalizePath(env_file),
  csv_file = normalizePath(csv_file),
  out_dir = normalizePath(out_dir),
  response_vars = response_vars,
  raster_predictor_vars = raster_predictor_vars,
  predictor_vars = predictor_vars,
  quantiles = quantiles,
  tuning_grid = tuning_grid,
  qrf_num_trees_tuning = qrf_num_trees_tuning,
  qrf_num_trees_final_candidates = qrf_num_trees_final_candidates,
  qrf_num_trees_validation_default = qrf_num_trees_validation_default,
  qrf_num_trees_final_default = qrf_num_trees_final_default,
  validation_repeats = validation_repeats,
  validation_fraction = validation_fraction,
  run_knndm_tuning = run_knndm_tuning,
  top_n_for_knndm = top_n_for_knndm,
  knndm_k = knndm_k,
  knndm_samplesize = knndm_samplesize,
  seed = seed,
  n_cores = n_cores,
  best_by_response = best$by_response,
  best_global = best$global_fair[1, ],
  best_global_all_candidates = best$global,
  terra_version = as.character(utils::packageVersion("terra")),
  ranger_version = as.character(utils::packageVersion("ranger")),
  CAST_version = if (requireNamespace("CAST", quietly = TRUE)) as.character(utils::packageVersion("CAST")) else NA_character_,
  sf_version = if (requireNamespace("sf", quietly = TRUE)) as.character(utils::packageVersion("sf")) else NA_character_
)
saveRDS(run_info, file.path(out_dir, "qrf_tuning_run_info.rds"))

message("==================================================")
message("QRF parameter tuning complete. No maps were produced.")
message("Outputs written to: ", normalizePath(out_dir))
message("Best fair global parameters (complete response coverage preferred):")
print(best$global_fair[1, ])
message("Tuned-parameter snippet: ",
        normalizePath(file.path(out_dir, "apply_best_qrf_params_to_final_mapping.R")))
message("predict_age_sensitivity.R will source this snippet automatically ",
        "(auto_load_tuned_params <- TRUE).")
message("Elapsed: ", format(round(difftime(Sys.time(), t_start, units = "mins"), 2)))
