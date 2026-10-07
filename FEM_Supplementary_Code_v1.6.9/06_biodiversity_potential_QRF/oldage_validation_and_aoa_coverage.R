# ============================================================
# oldage_validation_and_aoa_coverage.R
#
# 目的（对应讨论中未做的针对性检验）
# ------------------------------------------------------------
#   Part A. 留出"最老那批"样点验证：模型能否复原老龄样点的多样性？
#           直接检验"年龄反事实（固定 Forest_age=120）"在老龄端是否可靠。
#           不是普通空间 CV，而是沿【年龄轴】的外推检验。
#   Part B. 汇总 Forest_age_max 与各年龄 AOA 覆盖率，并补算管线没有产出的
#           "年龄外推区"（环境-only AOA 内、但 age=120 完整 AOA 外 = B 型区）。
#
# 设计：独立后处理，忠实复刻 predict_age_sensitivity.R 的
#   提取 → 订正填补 → 训练矩阵 → ranger QRF → 分位数预测 → 后处理，
#   但【不 source 那个 2900 行脚本】（否则会重跑整条生产管线）。
#   直接读同样的输入（ENV.tif / train4pot.csv）与已有输出（AOA/DI、调参片段）。
#
# 运行：把本文件放在 <FEM_POTENTIAL_DATA_DIR> 下，在生产管线跑完后执行
#   Rscript oldage_validation_and_aoa_coverage.R
# 依赖：terra, ranger（Part B 读 AOA 栅格用 terra）
# ============================================================

suppressPackageStartupMessages({
  library(terra)
})

# ============================================================
# 0. CONFIG（与 predict_age_sensitivity.R 对齐）
# ============================================================
# Portable paths: keep input/output paths relative to an explicit data root.
# Shared input/cache/output contracts are resolved relative to this installed script.
.fem_script<-sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1])
source(file.path(dirname(normalizePath(.fem_script)),"..","R","qrf_contracts.R"))
rm(.fem_script)

base_dir <- Sys.getenv("FEM_POTENTIAL_DATA_DIR", unset = getwd())
base_dir <- normalizePath(path.expand(base_dir), winslash = "/", mustWork = TRUE)
setwd(base_dir)

env_file        <- "ENV.tif"
csv_file        <- "train4pot.csv"
template_file   <- "templ_1km.tif"
# TRUE: templ_1km.tif 中 0 与 NA 都视为模板外；若 0 是有效像元值，请改为 FALSE。
mask_zero_as_na <- TRUE
prod_out_dir    <- "quantile_sensitivity_outputs"          # 生产脚本输出目录（读 AOA/DI）
tuning_out_dir  <- "qrf_parameter_tuning_outputs"          # 调参输出（读最优参数片段）
out_dir         <- "oldage_validation_outputs"             # 本脚本输出目录
fem_qrf_new_output(out_dir)

# ENV.tif 只读；DEM/BIO 会派生为共享的 templ_1km 对齐+mask 缓存，供调参/预测/oldage 复用。
aligned_env_file <- file.path(base_dir, "ENV_predictors_aligned_to_templ_1km_masked.tif")
reuse_aligned_env_cache <- TRUE

response_vars <- c("Rich_tree", "Shannon_wiener")
climate_vars  <- c("BIO6", "BIO10", "BIO12", "BIO17", "BIO15", "BIO4")
raster_predictor_vars <- c("DEM", climate_vars)
predictor_vars <- c(raster_predictor_vars, "Forest_age")

quantiles <- sort(unique(c(0.90, 0.95)))
q_names   <- paste0("q", round(quantiles * 100))

nonnegative_responses <- c("Rich_tree", "Shannon_wiener")
clip_negative_predictions <- TRUE

# ---- Part A 参数 ----
# 老龄留出阈值：同时用"分位数阈值(留出最老 10%/5%)"和"绝对阈值(100/120)"
oldage_percentile_thresholds <- numeric(0)   # 留出 age >= 该分位数的点
oldage_absolute_thresholds   <- c(100, 120, 148, 152)     # 留出 age >= 100 / >= 120 的点（点数够才用）
min_holdout_n <- 20                              # 留出点数下限，不足则跳过并记为"无法验证"
min_train_n   <- 100                             # 训练点数下限
oldage_within_kfold_min_age <- 100               # Test2：老龄内部 k 折的"老龄"定义
oldage_within_kfold_k <- 10
oldage_within_kfold_repeats <- 3

# ---- Part B 参数 ----
btype_reference_age <- 120                        # B 型区用哪个年龄的完整 AOA
compute_btype_zone  <- TRUE

# QRF 参数（默认；若存在调参片段则被其覆盖）
qrf_num_trees   <- 1500L
qrf_min_node_size <- 5L
qrf_mtry        <- max(2L, floor(sqrt(length(predictor_vars))))
qrf_always_split_variables <- "Forest_age"
qrf_prediction_chunk_size <- 200000L

seed <- 123
set.seed(seed)

# 线程
hardware_cores <- fem_qrf_threads(8L)
detected_cores <- parallel::detectCores(logical = TRUE)
if (is.na(detected_cores)) detected_cores <- hardware_cores
n_cores <- max(1L, min(hardware_cores, detected_cores))
Sys.setenv(OMP_NUM_THREADS = as.character(n_cores),
           MKL_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           VECLIB_MAXIMUM_THREADS = "1", NUMEXPR_NUM_THREADS = "1")

if (!requireNamespace("ranger", quietly = TRUE)) stop("需要 ranger 包：install.packages('ranger')")

t_start <- Sys.time()

# 可选：载入调参片段（只定义 qrf_params_by_response / 标量，安全）
# v1.6: SI S5.4 parameters are the default; a fresh tuning result is opt-in.
qrf_parameter_mode <- Sys.getenv("FEM_QRF_PARAMETER_MODE", "manuscript")
if (!qrf_parameter_mode %in% c("manuscript", "retuned")) stop("Invalid FEM_QRF_PARAMETER_MODE")
qrf_params_by_response <- NULL
tuned_snippet <- file.path(tuning_out_dir, "apply_best_qrf_params_to_final_mapping.R")
if (qrf_parameter_mode == "retuned" && !file.exists(tuned_snippet)) stop("Retuned profile requires: ", tuned_snippet)
if (qrf_parameter_mode == "retuned") {
  message("Sourcing tuned QRF params: ", tuned_snippet)
  local_trees <- qrf_num_trees
  source(tuned_snippet, local = FALSE)
  qrf_num_trees <- local_trees  # 片段可能改标量树数，这里保持本脚本设定
  if (exists("qrf_mtry")) qrf_mtry <- get("qrf_mtry")
  if (exists("qrf_min_node_size")) qrf_min_node_size <- get("qrf_min_node_size")
  if (exists("qrf_always_split_variables")) qrf_always_split_variables <- get("qrf_always_split_variables")
} else {
  message("Using fixed manuscript QRF parameters; no tuning snippet is applied.")
}

if (qrf_parameter_mode == "manuscript") {
  qrf_params_by_response <- list(
    Rich_tree = list(qrf_mtry=4L, qrf_min_node_size=10L,
      qrf_always_split_variables="Forest_age", qrf_num_trees_validation=1500L, qrf_num_trees_final=1500L),
    Shannon_wiener = list(qrf_mtry=2L, qrf_min_node_size=10L,
      qrf_always_split_variables=NULL, qrf_num_trees_validation=1500L, qrf_num_trees_final=1500L))
}

get_qrf_params <- function(response_var = NULL) {
  if (!is.null(qrf_params_by_response) &&
      !is.null(response_var) && response_var %in% names(qrf_params_by_response)) {
    pp <- qrf_params_by_response[[response_var]]
    return(list(
      mtry = if (!is.null(pp$qrf_mtry)) pp$qrf_mtry else qrf_mtry,
      min_node = if (!is.null(pp$qrf_min_node_size)) pp$qrf_min_node_size else qrf_min_node_size,
      always = if (!is.null(pp$qrf_always_split_variables)) pp$qrf_always_split_variables else NULL
    ))
  }
  list(mtry = qrf_mtry, min_node = qrf_min_node_size,
       always = if (identical(qrf_always_split_variables, "Forest_age")) "Forest_age" else NULL)
}

# ============================================================
# Validate the completed prediction before combining its AOA with new validation.
prod_run<-fem_qrf_completed_run(file.path(prod_out_dir,"run_info.rds"),
                               c(env_file,csv_file,template_file),qrf_parameter_mode)
for(resp in response_vars){
 actual<-get_qrf_params(resp);saved<-prod_run$qrf_params_by_response[[resp]]
 expected<-list(mtry=saved$qrf_mtry,min_node=saved$qrf_min_node_size,always=saved$qrf_always_split_variables)
 if(!isTRUE(all.equal(actual,expected,check.attributes=FALSE)))stop("Old-age parameters differ from saved prediction: ",resp)
}

# 1. 工具函数（复刻生产脚本口径）
# ============================================================
safe_as_numeric <- function(x) fem_qrf_numeric(x)

pinball_loss <- function(y, pred, q) {
  err <- y - pred
  mean(ifelse(err >= 0, q * err, (q - 1) * err), na.rm = TRUE)
}

enforce_quantile_order <- function(p) {
  # p: matrix/data.frame n x nq，列顺序对应 sort(quantiles)
  p <- as.matrix(p)
  if (ncol(p) > 1) for (i in 2:ncol(p)) p[, i] <- pmax(p[, i], p[, i - 1])
  p
}

postprocess_pred <- function(p, response_var) {
  p <- enforce_quantile_order(p)
  if (clip_negative_predictions && response_var %in% nonnegative_responses) p[p < 0] <- 0
  colnames(p) <- q_names
  p
}

# --- 订正填补：DEM->DEMc, BIO{x}->bio{x}_wc；lm(y ~ x + long + lat) + clamp（与生产一致）---
fill_env_na_from_csv_calibrated <- function(dat, target_vars,
                                            min_pairs = 30, min_r2 = 0.70,
                                            clamp = TRUE, direct_if_poor = TRUE) {
  fmap <- setNames(rep(NA_character_, length(target_vars)), target_vars)
  if ("DEM" %in% target_vars) fmap["DEM"] <- "DEMc"
  bio_t <- grep("^BIO[0-9]+$", target_vars, value = TRUE)
  fmap[bio_t] <- paste0("bio", sub("^BIO", "", bio_t), "_wc")
  fmap <- fmap[!is.na(fmap)]
  miss <- setdiff(unname(fmap), names(dat))
  if (length(miss)) stop("CSV 缺少订正填补字段：", paste(miss, collapse = ", "))
  dat$long <- safe_as_numeric(dat$long); dat$lat <- safe_as_numeric(dat$lat)

  for (target in names(fmap)) {
    src <- fmap[[target]]
    dat[[target]] <- safe_as_numeric(dat[[target]]); dat[[src]] <- safe_as_numeric(dat[[src]])
    bad <- is.na(dat[[target]]) | !is.finite(dat[[target]]); ok <- !bad
    sok <- !is.na(dat[[src]]) & is.finite(dat[[src]])
    xyok <- is.finite(dat$long) & is.finite(dat$lat)
    pair <- ok & sok & xyok; fill <- bad & sok & xyok
    if (sum(fill) == 0) next
    if (sum(pair) >= min_pairs) {
      cal <- data.frame(y = dat[[target]][pair], x = dat[[src]][pair],
                        long = dat$long[pair], lat = dat$lat[pair])
      nd  <- data.frame(x = dat[[src]][fill], long = dat$long[fill], lat = dat$lat[fill])
      fit <- try(stats::lm(y ~ x + long + lat, data = cal), silent = TRUE)
      done <- FALSE
      if (!inherits(fit, "try-error")) {
        fp <- try(stats::predict(fit, cal), silent = TRUE)
        r2 <- if (!inherits(fp, "try-error") && stats::sd(cal$y) > 0 && stats::sd(fp) > 0)
          suppressWarnings(stats::cor(cal$y, fp)^2) else NA_real_
        if (is.finite(r2) && r2 >= min_r2) {
          fv <- try(as.numeric(stats::predict(fit, nd)), silent = TRUE)
          if (!inherits(fv, "try-error")) {
            if (clamp) { fv <- pmax(pmin(fv, max(cal$y)), min(cal$y)) }
            idx <- which(fill); v <- is.finite(fv)
            dat[[target]][idx[v]] <- fv[v]; done <- TRUE
          }
        }
      }
      if (!done && direct_if_poor) dat[[target]][fill] <- dat[[src]][fill]
    } else if (direct_if_poor) {
      dat[[target]][fill] <- dat[[src]][fill]
    }
  }
  dat
}

make_training_data <- function(dat, response_var) {
  use_cols <- c(predictor_vars, response_var)
  miss <- setdiff(use_cols, names(dat))
  if (length(miss)) stop(response_var, " 训练缺列：", paste(miss, collapse = ", "))
  d <- dat[, use_cols, drop = FALSE]
  names(d)[names(d) == response_var] <- "response"
  for (nm in names(d)) d[[nm]] <- safe_as_numeric(d[[nm]])
  d <- d[rowSums(!is.finite(as.matrix(d)))==0L, , drop = FALSE]
  d
}

train_qrf <- function(train_dat, response_var, num_trees = qrf_num_trees, model_seed = seed) {
  pp <- get_qrf_params(response_var)
  d <- train_dat[, c(predictor_vars, "response"), drop = FALSE]
  args <- list(formula = response ~ ., data = d, num.trees = num_trees,
               mtry = pp$mtry, min.node.size = pp$min_node,
               quantreg = TRUE, importance = "none", oob.error = FALSE,
               num.threads = n_cores, seed = model_seed)
  if (!is.null(pp$always)) args$always.split.variables <- pp$always
  do.call(ranger::ranger, args)
}

predict_qrf <- function(model, newdata) {
  nd <- newdata[, predictor_vars, drop = FALSE]
  p <- predict(model, data = nd, type = "quantiles", quantiles = quantiles,
               num.threads = n_cores)$predictions
  if (is.vector(p)) p <- matrix(p, ncol = length(quantiles))
  p
}

# 评价一组分位数预测（覆盖率为主）
eval_quantile_pred <- function(y, pred_mat, label = "") {
  rows <- list()
  for (j in seq_along(quantiles)) {
    q <- quantiles[j]; pj <- pred_mat[, j]
    cov <- mean(y <= pj, na.rm = TRUE)                 # 经验覆盖率，目标≈q
    rows[[j]] <- data.frame(
      label = label, quantile = q,
      n = sum(is.finite(y) & is.finite(pj)),
      empirical_coverage = cov,
      coverage_gap = cov - q,                          # <0 表示欠覆盖（分位数偏低）
      pinball = pinball_loss(y, pj, q),
      mean_obs = mean(y, na.rm = TRUE),
      mean_pred = mean(pj, na.rm = TRUE),
      bias_pred_minus_obs = mean(pj - y, na.rm = TRUE),
      pct_obs_above_pred = 100 * mean(y > pj, na.rm = TRUE),  # 观测超过该分位数预测的比例
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

# ============================================================
# templ_1km 强制对齐工具
#   - 原始 ENV.tif 只读，不覆盖、不改写。
#   - 生成派生缓存 ENV_predictors_aligned_to_templ_1km_masked.tif。
#   - Part A 训练点提取使用原始 ENV.tif；Part B AOA/B型区叠加使用对齐并 mask 的 ENV_pred。
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
# 2. 读入 ENV + 样点，提取 + 订正填补（一次，供 Part A 与年龄统计）
# ============================================================
message("== 读取 ENV 与样点，提取环境变量 ==")
ENV_raw <- terra::rast(env_file)
message("Raw ENV.tif is read-only and will NOT be overwritten.")
miss_layer <- setdiff(raster_predictor_vars, names(ENV_raw))
if (length(miss_layer)) stop("ENV.tif 缺图层：", paste(miss_layer, collapse = ", "))
ENV_model_pred <- ENV_raw[[raster_predictor_vars]]
ENV_model_pred <- terra::ifel(is.finite(ENV_model_pred)&ENV_model_pred != -9999,ENV_model_pred,NA)  # Part A 训练/留出验证使用原始 ENV.tif，不受 templ_1km mask 限制

template_1km <- read_templ_1km_mask(template_file)
ENV_pred <- align_env_to_templ_1km_mask(
  ENV_model_pred,
  template_1km,
  aligned_file = aligned_env_file,
  reuse_cache = reuse_aligned_env_cache,
  source_files = c(env_file, template_file)
)
# 数据流：Part A 老龄留出验证用 ENV_model_pred；Part B AOA/B型区叠加用 ENV_pred/template_1km。

dat_raw <- utils::read.csv(csv_file, stringsAsFactors = FALSE)
for (cc in intersect(c("long","lat","Forest_age", response_vars,
                       "DEMc", paste0("bio", sub("^BIO","",climate_vars), "_wc")), names(dat_raw)))
  dat_raw[[cc]] <- safe_as_numeric(dat_raw[[cc]])

coord_ok <- is.finite(dat_raw$long) & is.finite(dat_raw$lat)
dat0 <- dat_raw[coord_ok, , drop = FALSE]

pts_ll  <- terra::vect(dat0, geom = c("long","lat"), crs = "EPSG:4326")
# 老龄验证的训练点环境变量提取使用原始 ENV.tif，不按 templ_1km mask 剔除样点。
pts_env <- terra::project(pts_ll, terra::crs(ENV_model_pred))
env_extract <- terra::extract(ENV_model_pred, pts_env, ID = FALSE)
dat0 <- dat0[, setdiff(names(dat0), raster_predictor_vars), drop = FALSE]
dat0 <- cbind(dat0, env_extract)
dat0 <- fill_env_na_from_csv_calibrated(dat0, raster_predictor_vars)

# 每个响应的建模训练矩阵（complete.cases，与生产一致）
train_list <- lapply(response_vars, function(r) make_training_data(dat0, r))
names(train_list) <- response_vars

# ============================================================
# 3. Part A —— 老龄留出验证
# ============================================================
message("== Part A：老龄留出验证 ==")

age_stats_rows <- list()
test1_rows <- list()   # 严格年龄外推留出
test2_rows <- list()   # 老龄内部随机 k 折

for (resp in response_vars) {
  d <- train_list[[resp]]
  age <- d$Forest_age
  amin <- min(age); amax <- max(age)

  age_stats_rows[[resp]] <- data.frame(
    response = resp, n_train = nrow(d),
    Forest_age_min = amin, Forest_age_max = amax,
    Forest_age_q10 = as.numeric(quantile(age, 0.10)),
    Forest_age_q50 = as.numeric(quantile(age, 0.50)),
    Forest_age_q90 = as.numeric(quantile(age, 0.90)),
    Forest_age_q95 = as.numeric(quantile(age, 0.95)),
    n_ge_100 = sum(age >= 100), n_ge_120 = sum(age >= 120),
    max_below_120 = if (any(age < 120)) max(age[age < 120]) else NA_real_,
    stringsAsFactors = FALSE
  )

  # ---- 组装阈值集合 ----
  thr_pct <- as.numeric(quantile(age, oldage_percentile_thresholds))
  thr_abs <- oldage_absolute_thresholds
  thresholds <- sort(unique(round(c(thr_pct, thr_abs), 3)))

  for (T in thresholds) {
    hold_idx  <- which(age >= T)
    train_idx <- which(age <  T)
    n_hold <- length(hold_idx); n_tr <- length(train_idx)
    if (n_hold < min_holdout_n || n_tr < min_train_n) {
      test1_rows[[paste(resp, T)]] <- data.frame(
        response = resp, threshold = T, n_holdout = n_hold, n_train = n_tr,
        status = "SKIPPED_too_few", quantile = NA_real_, empirical_coverage = NA_real_,
        coverage_gap = NA_real_, pinball = NA_real_, mean_obs = NA_real_, mean_pred = NA_real_,
        bias_pred_minus_obs = NA_real_, pct_obs_above_pred = NA_real_,
        train_max_age = if (n_tr > 0) max(age[train_idx]) else NA_real_,
        clamp_engaged = NA, pred_mode = NA_character_, stringsAsFactors = FALSE
      )
      next
    }

    d_tr <- d[train_idx, , drop = FALSE]
    d_ho <- d[hold_idx,  , drop = FALSE]
    train_max_age <- max(d_tr$Forest_age)
    m <- train_qrf(d_tr, resp, num_trees = qrf_num_trees, model_seed = seed + round(T))

    # (i) 生产口径：输入年龄 clamp 到训练最大年龄（模拟"预测没见过的老龄"）
    d_ho_clamp <- d_ho
    d_ho_clamp$Forest_age <- pmin(d_ho_clamp$Forest_age, train_max_age)
    p_clamp <- postprocess_pred(predict_qrf(m, d_ho_clamp), resp)
    e_clamp <- eval_quantile_pred(d_ho$response, p_clamp,
                                  label = sprintf("%s_T%.1f_clamped", resp, T))
    e_clamp$response <- resp; e_clamp$threshold <- T
    e_clamp$n_holdout <- n_hold; e_clamp$n_train <- n_tr
    e_clamp$train_max_age <- train_max_age
    e_clamp$clamp_engaged <- TRUE; e_clamp$pred_mode <- "clamped_to_train_max"
    e_clamp$status <- "OK"

    # (ii) 不夹紧输入年龄（保留真实老龄）：显示 RF 内部对响应外推的天花板
    p_raw <- postprocess_pred(predict_qrf(m, d_ho), resp)
    e_raw <- eval_quantile_pred(d_ho$response, p_raw,
                                label = sprintf("%s_T%.1f_uncapped", resp, T))
    e_raw$response <- resp; e_raw$threshold <- T
    e_raw$n_holdout <- n_hold; e_raw$n_train <- n_tr
    e_raw$train_max_age <- train_max_age
    e_raw$clamp_engaged <- FALSE; e_raw$pred_mode <- "uncapped_input_age"
    e_raw$status <- "OK"

    test1_rows[[paste(resp, T)]] <- rbind(
      e_clamp[, c("response","threshold","n_holdout","n_train","train_max_age",
                  "clamp_engaged","pred_mode","status","quantile","empirical_coverage",
                  "coverage_gap","pinball","mean_obs","mean_pred","bias_pred_minus_obs","pct_obs_above_pred")],
      e_raw[,  c("response","threshold","n_holdout","n_train","train_max_age",
                 "clamp_engaged","pred_mode","status","quantile","empirical_coverage",
                 "coverage_gap","pinball","mean_obs","mean_pred","bias_pred_minus_obs","pct_obs_above_pred")]
    )
    message(sprintf("  [Test1] %s T=%.1f: holdout n=%d, train n=%d, train_max_age=%.1f",
                    resp, T, n_hold, n_tr, train_max_age))
  }

  # ---- Test2：老龄内部随机 k 折（模型仍见到部分老龄点，作为对照）----
  old_idx <- which(age >= oldage_within_kfold_min_age)
  if (length(old_idx) >= max(min_holdout_n, oldage_within_kfold_k * 5)) {
    for (rep_i in seq_len(oldage_within_kfold_repeats)) {
      set.seed(seed + 777 * rep_i)
      folds <- sample(rep(seq_len(oldage_within_kfold_k), length.out = length(old_idx)))
      preds <- matrix(NA_real_, nrow = length(old_idx), ncol = length(quantiles))
      yy <- d$response[old_idx]
      for (k in seq_len(oldage_within_kfold_k)) {
        te <- old_idx[folds == k]
        # 训练和测试都限定在 age >=100；排除当前测试折。
        tr <- setdiff(old_idx, te)  # SI: train and test within age >= 100 subset
        mk <- train_qrf(d[tr, , drop = FALSE], resp,
                        num_trees = qrf_num_trees, model_seed = seed + 100 * rep_i + k)
        pk <- postprocess_pred(predict_qrf(mk, d[te, , drop = FALSE]), resp)
        preds[folds == k, ] <- pk
      }
      e2 <- eval_quantile_pred(yy, preds, label = sprintf("%s_within_old_rep%d", resp, rep_i))
      e2$response <- resp; e2$repeat_id <- rep_i; e2$n_old <- length(old_idx)
      e2$min_age_def <- oldage_within_kfold_min_age
      test2_rows[[paste(resp, rep_i)]] <- e2
    }
    message(sprintf("  [Test2] %s: 老龄内部 %d 折 x %d 重复 (n_old=%d)",
                    resp, oldage_within_kfold_k, oldage_within_kfold_repeats, length(old_idx)))
  } else {
    message(sprintf("  [Test2] %s: 老龄点不足 (n=%d)，跳过内部 k 折。", resp, length(old_idx)))
  }
}

age_stats <- do.call(rbind, age_stats_rows)
test1 <- do.call(rbind, test1_rows); rownames(test1) <- NULL
test2 <- if (length(test2_rows)) do.call(rbind, test2_rows) else data.frame()

write.csv(age_stats, file.path(out_dir, "training_forest_age_distribution.csv"), row.names = FALSE)
write.csv(test1, file.path(out_dir, "oldage_holdout_extrapolation_metrics.csv"), row.names = FALSE)
if (nrow(test2)) write.csv(test2, file.path(out_dir, "oldage_within_random_kfold_metrics.csv"), row.names = FALSE)

# ---- Part A 解读表（把数字翻成结论）----
interp_rows <- list()
for (resp in response_vars) {
  sub <- test1[test1$response == resp & test1$status == "OK" &
               test1$pred_mode == "clamped_to_train_max" & test1$quantile == 0.95, , drop = FALSE]
  if (nrow(sub) == 0) next
  # 取最严格阈值（最高 threshold）的一行作代表
  rep_row <- sub[which.max(sub$threshold), , drop = FALSE]
  amax <- age_stats$Forest_age_max[age_stats$response == resp]
  verdict <- with(rep_row, {
    if (amax < 120) "训练最大年龄 < 120：120 面实为 clamp 到 age_max 的面，不能作 120 年结论"
    else if (empirical_coverage < 0.85 && pct_obs_above_pred > 10)
      "老龄端严重欠覆盖且大量实测超过 q95：年龄外推不可靠，潜力面是下界"
    else if (empirical_coverage < 0.90)
      "老龄端轻度欠覆盖：仅在 AOA 内可谨慎解释老龄潜力"
    else "老龄端覆盖接近目标：年龄反事实在此端相对可信"
  })
  interp_rows[[resp]] <- data.frame(
    response = resp, Forest_age_max = amax,
    strictest_threshold = rep_row$threshold,
    holdout_n = rep_row$n_holdout, train_max_age = rep_row$train_max_age,
    q95_empirical_coverage = rep_row$empirical_coverage,
    q95_pct_obs_above_pred = rep_row$pct_obs_above_pred,
    verdict = verdict, stringsAsFactors = FALSE
  )
}
if (length(interp_rows)) {
  interp <- do.call(rbind, interp_rows)
  write.csv(interp, file.path(out_dir, "oldage_validation_interpretation.csv"), row.names = FALSE)
}

# ============================================================
# 4. Part B —— AOA 覆盖率汇总 + 年龄外推区（B 型）
# ============================================================
message("== Part B：AOA 覆盖率 + 年龄外推区 ==")

# 4a. 汇总生产已算的 AOA 覆盖率
prod_aoa_summary <- file.path(prod_out_dir, "AOA_summary_all_responses_ages.csv")
if (file.exists(prod_aoa_summary)) {
  aoa_cov <- utils::read.csv(prod_aoa_summary, stringsAsFactors = FALSE)
  write.csv(aoa_cov, file.path(out_dir, "AOA_coverage_by_age_ECHO.csv"), row.names = FALSE)
  message("  已回显生产 AOA 覆盖率：", prod_aoa_summary)
} else {
  aoa_cov <- NULL
  message("  未找到 ", prod_aoa_summary, "，跳过 AOA 覆盖率回显（可改为直接读 AOA 栅格）。")
}

# 4b. 年龄外推区：env-only AOA==1 且 age=120 完整 AOA==0
btype_rows <- list()
if (compute_btype_zone) {
  env_only_file <- file.path(prod_out_dir, "AOA_environment_only_AOA.tif")
  if (!file.exists(env_only_file)) {
    message("  未找到环境-only AOA (", env_only_file, ")，跳过 B 型区计算。")
  } else {
    env_aoa <- align_raster_to_templ_1km_mask(terra::rast(env_only_file), template_1km, method = "near")
    for (resp in response_vars) {
      full_file <- file.path(prod_out_dir,
                             paste0("AOA_", resp, "_ForestAge_", btype_reference_age, "_AOA.tif"))
      if (!file.exists(full_file)) {
        message("  未找到 ", full_file, "，跳过 ", resp, " 的 B 型区。")
        next
      }
      full_aoa <- align_raster_to_templ_1km_mask(terra::rast(full_file), template_1km, method = "near")
      # env-only AOA 和 full-age AOA 均已强制对齐到 templ_1km 并使用 templ_1km mask。
      valid <- !is.na(env_aoa) & !is.na(full_aoa)
      # B 型：气候没问题(env-only 内=1) 但 该年龄外推(full=0)
      btype <- env_aoa == 1 & full_aoa == 0
      btype[!valid] <- NA
      out_tif <- file.path(out_dir, paste0("age_extrapolation_zone_", resp,
                                           "_age", btype_reference_age, ".tif"))
      terra::writeRaster(btype, out_tif, overwrite = TRUE,
                         datatype = "INT1U", gdal = c("COMPRESS=LZW"))
      nv <- terra::global(valid, "sum", na.rm = TRUE)[1, 1]
      nin_env <- terra::global(env_aoa == 1 & valid, "sum", na.rm = TRUE)[1, 1]
      nin_full <- terra::global(full_aoa == 1 & valid, "sum", na.rm = TRUE)[1, 1]
      nb <- terra::global(btype, "sum", na.rm = TRUE)[1, 1]
      btype_rows[[resp]] <- data.frame(
        response = resp, reference_age = btype_reference_age,
        valid_cells = nv,
        env_only_inside_pct = 100 * nin_env / nv,
        full_age_inside_pct  = 100 * nin_full / nv,
        age_extrapolation_zone_pct = 100 * nb / nv,   # ★ 决策核心：气候可、但该年龄外推的比例
        note = "B型: env-only AOA内 且 age完整AOA外 = 气候有支撑但该林龄未见过",
        raster = normalizePath(out_tif), stringsAsFactors = FALSE
      )
      message(sprintf("  [B型] %s @age%d：外推区占有效像元 %.1f%%",
                      resp, btype_reference_age, 100 * nb / nv))
    }
  }
}
if (length(btype_rows)) {
  btype_tab <- do.call(rbind, btype_rows)
  write.csv(btype_tab, file.path(out_dir, "age_extrapolation_zone_summary.csv"), row.names = FALSE)
}

# 4c. 生产训练年龄摘要回显（若存在）
prod_train_summary <- file.path(prod_out_dir, "model_training_summary.csv")
if (file.exists(prod_train_summary)) {
  file.copy(prod_train_summary,
            file.path(out_dir, "model_training_summary_ECHO.csv"), overwrite = TRUE)
}

# ============================================================
# 5. 控制台总判决
# ============================================================
message("\n==================== 决策摘要 ====================")
for (resp in response_vars) {
  amax <- age_stats$Forest_age_max[age_stats$response == resp]
  n120 <- age_stats$n_ge_120[age_stats$response == resp]
  msg <- sprintf("[%s] 训练 Forest_age_max=%.1f, n(age>=120)=%d", resp, amax, n120)
  if (length(btype_rows) && resp %in% names(btype_rows)) {
    msg <- paste0(msg, sprintf(", age120 完整AOA覆盖=%.1f%%, 年龄外推区=%.1f%%",
                               btype_rows[[resp]]$full_age_inside_pct,
                               btype_rows[[resp]]$age_extrapolation_zone_pct))
  }
  message(msg)
}
message("详见 ", normalizePath(out_dir))
message("关键文件：")
message("  training_forest_age_distribution.csv        —— age_max / n>=120（决定120面是否成立）")
message("  oldage_holdout_extrapolation_metrics.csv     —— 老龄端 q90/q95 覆盖率与欠覆盖")
message("  oldage_validation_interpretation.csv         —— 数字→结论")
message("  age_extrapolation_zone_summary.csv + *.tif    —— B型区（气候可但年龄外推）")
message("Elapsed: ", format(round(difftime(Sys.time(), t_start, units = "mins"), 2)))
