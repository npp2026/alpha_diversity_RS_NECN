#!/usr/bin/env Rscript
# =============================================================================
# 01_fit_final_RFM_and_QM.R   （推荐主方案：production + 保留 Shannon=0 + 完整 QM；无残差校正）
# Stage-1 训练 + 上尾压缩诊断 + QM 去衰减标定，一次运行全部产出：
#   final_models_v1/[<sub>/]<target>_Final_Model.rds          模型
#   final_models_v1/[<sub>/]<target>_Importance.csv           变量重要性
#   final_models_v1/[<sub>/]<target>_Training_Metadata.csv    训练元数据
#   train_data/[<sub>/]<target>_train.rds                     清洗后训练表(含 oob，可复用/分层)
#   qm_diag/[<sub>/]<target>_metrics.csv                      订正前/后指标
#   qm_diag/[<sub>/]<target>_oob_corrected.csv                逐样本 obs/raw/corrected
#   qm_calibration/[<sub>/]qm_<target>.rds                    QM 校正（预测脚本 load_qm_corrector 直接读）
#
# 诊断用【训练后内存里的 model$predictions(OOB) 与清洗后观测】，obs↔oob 天然行对齐。
# qm_<target>.rds is consumed by 03_predict_annual_30m_maps.R.
# =============================================================================

suppressPackageStartupMessages({
  library(ranger)
  library(dplyr)
})
options(stringsAsFactors = FALSE)

# =============================================================================
# 1. CONFIGURATION
# =============================================================================

module_dir <- dirname(normalizePath(sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1])))
source(file.path(module_dir,"..","R","qm_core.R"))
source(file.path(module_dir,"..","R","input_contracts.R"))

DATA_FILE <- Sys.getenv("DATA_FILE", unset = "landsat_v12_1D_VRT_filtered.csv")

GLOBAL <- list(
  CORES          = suppressWarnings(as.integer(Sys.getenv("RFM_CORES", unset = "12"))),
  SEED           = 42L,
  OUTPUT_DIR     = Sys.getenv("MODEL_OUTPUT_DIR", unset = "final_models_v1"),   # 与预测脚本 CONFIG$MODEL_DIR 一致
  TRAIN_DATA_DIR = Sys.getenv("TRAIN_DATA_DIR", unset = "train_data"),
  DIAG_DIR       = Sys.getenv("QM_DIAG_DIR", unset = "qm_diag"),
  CALIB_DIR      = Sys.getenv("CALIB_DIR", unset = "qm_calibration"),    # 与预测脚本 load_qm_corrector 一致
  N_KNOTS        = 1000L,
  TAIL_PROBS     = c(0.90, 0.95, 0.99),
  CLAMP_NONNEG   = TRUE,
  REMOVE_SHANNON_ZERO = FALSE,        # 推荐主方案：保留 Shannon=0 单物种样点
  ENABLE_RESIDUAL_CORRECTION = FALSE, # 推荐主方案：不启用 OOB 残差校正
  QM_MODE        = "full",            # 推荐主方案：完整经验分位数映射，不做 shrinkage
  BUILD_QM       = TRUE,                # 训练后做诊断+QM标定
  RUN_CONFIGS    = c("production")      # 默认只跑生产配置；可加 "rGRVI","noNIRStd","rLANDSAT"
)

# 每目标超参（与三套原脚本一致）+ QM 物理上界 max
BASE_PARAMS <- list(
  Rich_tree      = list(num.trees = 2000, mtry = 5, min.node.size = 5, sample.fraction = 0.65, max = 150),
  Shannon_wiener = list(num.trees = 2000, mtry = 5, min.node.size = 5, sample.fraction = 0.65, max = 4.5)
)

# Production uses the author-confirmed fixed 15 predictors below for both targets.
# Module 02 VSURF screening belongs to validation only; its output is never read
# here and never automatically replaces this production list.
# 四套特征配置（subdir=NULL => 写 flat = 生产配置，喂预测/QM）
CONFIGS <- list(
  production = list(
    subdir = NULL,
    predictors = c("NIR_Std_GS", "GRVI_Med_Aut", "NDMI_CV_900m_Sum", "NIR_TSD_Spr",
                   "EVI_mean_w7_Spr", "Prod_GPP_CV", "Prod_GPP_Mean", "Prod_LAI_Max",
                   "Pheno_Amp", "Pheno_GSL", "BIO12_baseline_mean", "BIO4_baseline_mean",
                   "ELEV", "SPEI_12_annual_baseline_mean", "STN")
  ),
  rGRVI = list(
    subdir = "rGRVI",
    predictors = c("NIR_Std_GS", "NDMI_CV_900m_Sum", "NIR_TSD_Spr",
                   "EVI_mean_w7_Spr", "Prod_GPP_CV", "Prod_GPP_Mean", "Prod_LAI_Max",
                   "Pheno_Amp", "Pheno_GSL", "BIO12_baseline_mean", "BIO4_baseline_mean",
                   "ELEV", "SPEI_12_annual_baseline_mean", "STN")
  ),
  noNIRStd = list(
    subdir = "noNIRStd",
    predictors = c("NDMI_CV_900m_Sum", "NIR_TSD_Spr",
                   "EVI_mean_w7_Spr", "Prod_GPP_CV", "Prod_GPP_Mean", "Prod_LAI_Max",
                   "Pheno_Amp", "Pheno_GSL", "BIO12_baseline_mean", "BIO4_baseline_mean",
                   "ELEV", "SPEI_12_annual_baseline_mean", "STN")
  ),
  rLANDSAT = list(
    subdir = "rLANDSAT",
    predictors = c("Prod_GPP_CV", "Prod_GPP_Mean", "Prod_LAI_Max",
                   "Pheno_Amp", "Pheno_GSL", "BIO12_baseline_mean", "BIO4_baseline_mean",
                   "ELEV", "SPEI_12_annual_baseline_mean", "STN")
  )
)

TARGETS <- c("Rich_tree", "Shannon_wiener")
if (length(CONFIGS$production$predictors) != 15L || anyDuplicated(CONFIGS$production$predictors))
  stop("Production requires exactly 15 distinct fixed predictors")

# =============================================================================
# 2. HELPERS
# =============================================================================

log_msg <- function(...) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste(...))); flush.console()
}

check_columns <- function(data, cols, config_name, target) {
  miss <- setdiff(cols, names(data))
  if (length(miss) > 0)
    stop(sprintf("[%s - %s] DATA_FILE 缺列: %s", config_name, target, paste(miss, collapse = ", ")))
}

# =============================================================================
# 3. QM 核心 + 诊断（使用 ../R/qm_core.R 共享实现）
# =============================================================================

# 全套压缩诊断指标
diag_metrics <- function(obs, pred, tail_probs) {
  ok <- is.finite(obs) & is.finite(pred); obs <- obs[ok]; pred <- pred[ok]
  res <- obs - pred; sst <- sum((obs - mean(obs))^2)
  fit <- lm(obs ~ pred)
  thr <- quantile(obs, min(tail_probs), names = FALSE); ts <- obs >= thr
  tslope <- if (sum(ts) >= 30 && sd(pred[ts]) > 0) unname(coef(lm(obs[ts] ~ pred[ts]))[2]) else NA_real_
  qo <- quantile(obs, tail_probs, names = FALSE); qp <- quantile(pred, tail_probs, names = FALSE)
  list(n = length(obs), r = cor(obs, pred), R2 = 1 - sum(res^2) / sst,
       RMSE = sqrt(mean(res^2)), MAE = mean(abs(res)), ME = mean(res),
       sd_obs = sd(obs), sd_pred = sd(pred), lambda_disp = sd(pred) / sd(obs),
       calib_slope = unname(coef(fit)[2]), calib_int = unname(coef(fit)[1]),
       max_obs = max(obs), max_pred = max(pred),
       tail_probs = tail_probs, q_obs = qo, q_pred = qp, tail_ratio = qp / qo, tail_slope = tslope)
}

# 单调经验分位映射（predicted -> corrected）
build_qm <- function(obs, oob, n_knots, ymin, ymax, response) {
  q <- fem_fit_qm(oob,obs,response,GLOBAL$SEED,n_knots)
  if (ymin != q$ymin || ymax != q$ymax) stop("QM bounds differ from the response contract")
  q$apply <- function(v) fem_apply_qm(q,v)
  q
}

# 诊断 + 构建/保存 QM（用内存中的 obs/oob；写 metrics/oob_corrected + qm_<target>.rds）
run_qm_diagnostic <- function(name, obs, oob, ymax, diag_dir, calib_dir) {
  obs <- as.numeric(obs); oob <- as.numeric(oob)
  ok <- is.finite(obs) & is.finite(oob); obs <- obs[ok]; oob <- oob[ok]
  ymin <- if (isTRUE(GLOBAL$CLAMP_NONNEG)) 0 else min(obs)
  tp <- GLOBAL$TAIL_PROBS

  d0 <- diag_metrics(obs, oob, tp)
  qm <- build_qm(obs, oob, GLOBAL$N_KNOTS, ymin, ymax, name)
  oob_c <- qm$apply(oob)
  d1 <- diag_metrics(obs, oob_c, tp)

  # ---- 控制台报告 ----
  cat(sprintf("[QM] %s  n=%d\n", name, d0$n))
  cat("                        before    after\n")
  cat(sprintf("  Pearson r            %8.3f %8.3f\n", d0$r, d1$r))
  cat(sprintf("  R^2                  %8.3f %8.3f   (QM 后可能略降，正常)\n", d0$R2, d1$R2))
  cat(sprintf("  RMSE                 %8.3f %8.3f   (QM 后可能略升:偏差-方差取舍)\n", d0$RMSE, d1$RMSE))
  cat(sprintf("  MAE                  %8.3f %8.3f\n", d0$MAE, d1$MAE))
  cat(sprintf("  lambda_disp(sdP/sdO) %8.3f %8.3f   (-> ~1: 压缩消除)\n", d0$lambda_disp, d1$lambda_disp))
  cat(sprintf("  calib slope obs~pred %8.3f %8.3f\n", d0$calib_slope, d1$calib_slope))
  cat(sprintf("  max(obs)=%.3f  max(pred raw)=%.3f   (RF 不超训练上界)\n", d0$max_obs, d0$max_pred))
  for (i in seq_along(tp))
    cat(sprintf("  上尾 q%.2f Qobs=%.3f  pred/obs(lambda_hat): %.3f -> %.3f\n",
                tp[i], d0$q_obs[i], d0$tail_ratio[i], d1$tail_ratio[i]))

  dir.create(diag_dir, recursive = TRUE, showWarnings = FALSE)

  # ---- 指标 CSV ----
  flat <- function(d, tag) data.frame(
    stage = tag, n = d$n, r = d$r, R2 = d$R2, RMSE = d$RMSE, MAE = d$MAE, ME = d$ME,
    sd_obs = d$sd_obs, sd_pred = d$sd_pred, lambda_disp = d$lambda_disp,
    calib_slope = d$calib_slope, calib_int = d$calib_int, max_obs = d$max_obs, max_pred = d$max_pred,
    q90_obs = d$q_obs[1], q90_pred = d$q_pred[1], tail_ratio_q90 = d$tail_ratio[1],
    q95_obs = d$q_obs[2], q95_pred = d$q_pred[2], tail_ratio_q95 = d$tail_ratio[2],
    q99_obs = d$q_obs[3], q99_pred = d$q_pred[3], tail_ratio_q99 = d$tail_ratio[3],
    tail_slope = d$tail_slope)
  write.csv(rbind(flat(d0, "before"), flat(d1, "after")),
            file.path(diag_dir, sprintf("%s_metrics.csv", name)), row.names = FALSE)
  write.csv(data.frame(obs = obs, oob_raw = oob, oob_corrected = oob_c),
            file.path(diag_dir, sprintf("%s_oob_corrected.csv", name)), row.names = FALSE)

  # ---- 保存与预测脚本兼容的 QM 校正 ----
  dir.create(calib_dir, recursive = TRUE, showWarnings = FALSE)
  calib <- list(target = name, x_knots = qm$x_knots, y_knots = qm$y_knots, ymin = ymin, ymax = ymax,
                diag_before = d0, diag_after = d1, n = d0$n, built_at = Sys.time(),
                method=qm$method, seed=qm$seed, continuity_correction=qm$continuity_correction)
  saveRDS(calib, file.path(calib_dir, sprintf("qm_%s.rds", name)))
  log_msg(sprintf("QM saved: %s | diag: %s/%s_{metrics.csv,oob_corrected.csv}",
                  file.path(calib_dir, sprintf("qm_%s.rds", name)), diag_dir, name))
  invisible(calib)
}

# =============================================================================
# 4. TRAIN （与原脚本一致；返回 model + 清洗后 df）
# =============================================================================

train_final_model <- function(data, target, predictors, params, config_name) {
  log_msg(sprintf("--- Training %s | target=%s ---", config_name, target))
  log_msg(sprintf("Params: trees=%d, mtry=%d, min_node=%d, sample_frac=%.2f",
                  params$num.trees, params$mtry, params$min.node.size, params$sample.fraction))
  log_msg(sprintf("Predictors (%d): %s", length(predictors), paste(predictors, collapse = ", ")))

  check_columns(data, c(predictors, target), config_name, target)

  # 推荐主方案说明：production 特征集；Shannon=0 默认保留；不启用残差校正；训练后做完整经验 QM。
  model_cols <- c(predictors, target)
  model_df <- data[, model_cols, drop = FALSE]
  complete_idx <- fem_finite_rows(model_df, model_cols)
  df <- model_df[complete_idx, , drop = FALSE]
  if (nrow(df) == 0) stop(sprintf("[%s - %s] complete.cases 后无完整样本。", config_name, target))

  if (target == "Shannon_wiener") {
    zero_idx <- is.finite(df[[target]]) & df[[target]] == 0
    n_zero <- sum(zero_idx, na.rm = TRUE)
    msg <- sprintf("Shannon=0 complete training rows before optional removal: %d", n_zero)
    if ("Rich_tree" %in% names(data)) {
      rich_vals <- data$Rich_tree[complete_idx][zero_idx]
      n_rich1 <- sum(is.finite(rich_vals) & rich_vals == 1, na.rm = TRUE)
      msg <- sprintf("%s; among them Rich_tree==1: %d/%d", msg, n_rich1, n_zero)
    }
    log_msg(msg)
    if (isTRUE(GLOBAL$REMOVE_SHANNON_ZERO) && n_zero > 0) {
      df <- df[!zero_idx, , drop = FALSE]
      log_msg(sprintf("REMOVE_SHANNON_ZERO=TRUE: removed %d Shannon=0 rows", n_zero))
    } else {
      log_msg("REMOVE_SHANNON_ZERO=FALSE: keeping Shannon=0 rows in training data")
    }
    if (nrow(df) == 0) stop(sprintf("[%s - %s] Shannon=0 处理后无完整样本。", config_name, target))
  }

  # SI does not prescribe response-tail trimming. Retain valid high diversity
  # observations; these also define the OOD training reference and QM tails.
  if (any(df[[target]] < 0)) stop("Negative diversity response")

  set.seed(GLOBAL$SEED)
  t0 <- Sys.time()
  model <- ranger::ranger(
    formula = as.formula(paste(target, "~ .")), data = df,
    num.trees = params$num.trees, mtry = params$mtry,
    min.node.size = params$min.node.size, sample.fraction = params$sample.fraction,
    importance = "permutation", write.forest = TRUE, quantreg = TRUE,
    keep.inbag = TRUE, num.threads = GLOBAL$CORES, verbose = TRUE
  )
  log_msg(sprintf("Done in %.2fs. OOB R2: %.4f | OOB RMSE: %.4f",
                  as.numeric(Sys.time() - t0, units = "secs"),
                  model$r.squared, sqrt(model$prediction.error)))
  list(model = model, data = df)
}

# =============================================================================
# 5. SAVE OUTPUTS + 内联 QM 诊断/标定
# =============================================================================

save_outputs <- function(res, target, predictors, params, cfg, config_name) {
  is_flat <- is.null(cfg$subdir)
  out_dir   <- if (is_flat) GLOBAL$OUTPUT_DIR     else file.path(GLOBAL$OUTPUT_DIR, cfg$subdir)
  train_dir <- if (is_flat) GLOBAL$TRAIN_DATA_DIR else file.path(GLOBAL$TRAIN_DATA_DIR, cfg$subdir)
  diag_dir  <- if (is_flat) GLOBAL$DIAG_DIR       else file.path(GLOBAL$DIAG_DIR, cfg$subdir)
  calib_dir <- if (is_flat) GLOBAL$CALIB_DIR      else file.path(GLOBAL$CALIB_DIR, cfg$subdir)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  # --- 模型 ---
  saveRDS(res$model, file.path(out_dir, sprintf("%s_Final_Model.rds", target)))
  # --- 变量重要性 ---
  imp <- data.frame(Variable = names(ranger::importance(res$model)),
                    Importance = as.numeric(ranger::importance(res$model)), row.names = NULL) %>%
    dplyr::arrange(dplyr::desc(Importance))
  write.csv(imp, file.path(out_dir, sprintf("%s_Importance.csv", target)), row.names = FALSE)
  # --- 元数据 ---
  meta <- data.frame(config = config_name, target = target,
                     predictor_selection = if (config_name == "production") "fixed_15" else "predefined_sensitivity",
                     vsurf_screening_used = FALSE,
                     qm_calibration_source = "training_OOB",
                     n_predictors = length(predictors), predictors = paste(predictors, collapse = ";"),
                     num.trees = params$num.trees, mtry = params$mtry,
                     min.node.size = params$min.node.size, sample.fraction = params$sample.fraction,
                     n_training_rows_after_cleaning = nrow(res$data),
                     remove_shannon_zero = GLOBAL$REMOVE_SHANNON_ZERO,
                     enable_residual_correction = GLOBAL$ENABLE_RESIDUAL_CORRECTION,
                     qm_mode = GLOBAL$QM_MODE,
                     oob_r2 = res$model$r.squared, oob_rmse = sqrt(res$model$prediction.error))
  write.csv(meta, file.path(out_dir, sprintf("%s_Training_Metadata.csv", target)), row.names = FALSE)

  # --- 训练表 + OOB（行序与 model$predictions 对齐；供复用/分层 QM）---
  oob <- res$model$predictions
  if (is.matrix(oob)) oob <- rowMeans(oob, na.rm = TRUE)
  oob <- as.numeric(oob)
  dir.create(train_dir, recursive = TRUE, showWarnings = FALSE)
  train_tab <- res$data; train_tab$oob <- oob
  saveRDS(train_tab, file.path(train_dir, sprintf("%s_train.rds", target)))
  log_msg(sprintf("Saved model + importance + metadata + train table (n=%d)", nrow(res$data)))

  # --- 内联 QM 诊断 + 标定（用内存中的 obs/oob）---
  if (isTRUE(GLOBAL$BUILD_QM)) {
    run_qm_diagnostic(name = target, obs = res$data[[target]], oob = oob,
                      ymax = params$max, diag_dir = diag_dir, calib_dir = calib_dir)
  }
}

# =============================================================================
# 6. EXECUTION
# =============================================================================

if (!file.exists(DATA_FILE)) stop("Data file not found: ", DATA_FILE)
log_msg("Reading data:", DATA_FILE)
full_data <- read.csv(DATA_FILE, check.names = FALSE, colClasses=c(plot_id="character"))
log_msg(sprintf("Loaded %d rows × %d columns", nrow(full_data), ncol(full_data)))
# Fresh output directories prevent model/train/QM objects from different runs mixing.
for (dest in c(GLOBAL$OUTPUT_DIR,GLOBAL$TRAIN_DATA_DIR,GLOBAL$DIAG_DIR,GLOBAL$CALIB_DIR)) {
  if (dir.exists(dest) && length(list.files(dest,all.files=TRUE,no..=TRUE))) stop("Use empty final-model output directories: ",dest)
  dir.create(dest,recursive=TRUE,showWarnings=FALSE)
}
matched <- fem_match_plots(full_data)
all_model_cols <- unique(c(unlist(lapply(CONFIGS[GLOBAL$RUN_CONFIGS],`[[`,"predictors")),TARGETS))
finite <- fem_finite_rows(matched$data,all_model_cols)
included <- matched$included
included[which(included)] <- finite
write.csv(data.frame(input_row=seq_len(nrow(full_data)),plot_id=full_data$plot_id,
                     included=included),file.path(GLOBAL$OUTPUT_DIR,"training_inclusion.csv"),row.names=FALSE)
full_data <- full_data[included,,drop=FALSE]
if(nrow(full_data)<10L) stop("Fewer than 10 paired finite matched training rows")
if(any(full_data$Rich_tree<0 | full_data$Shannon_wiener<0)) stop("Negative diversity response")
if(!is.finite(GLOBAL$CORES)||GLOBAL$CORES<1L) stop("RFM_CORES must be positive")
write.csv(full_data[,c("plot_id","plot_year")],file.path(GLOBAL$OUTPUT_DIR,"training_plot_ids.csv"),row.names=FALSE)


for (config_name in GLOBAL$RUN_CONFIGS) {
  cfg <- CONFIGS[[config_name]]
  if (is.null(cfg)) stop("未知配置: ", config_name)
  cat(sprintf("\n################ CONFIG: %s ################\n", config_name))
  for (target in TARGETS) {
    params <- BASE_PARAMS[[target]]
    res <- train_final_model(full_data, target, cfg$predictors, params, config_name)
    save_outputs(res, target, cfg$predictors, params, cfg, config_name)
    rm(res); gc(verbose = FALSE)
  }
}

cat("\n========================================\n")
cat("PIPELINE COMPLETE\n")
cat("Models       -> ", GLOBAL$OUTPUT_DIR, "/[<sub>/]\n", sep = "")
cat("Train tables -> ", GLOBAL$TRAIN_DATA_DIR, "/[<sub>/]<target>_train.rds (含 oob)\n", sep = "")
cat("QM diagnostics-> ", GLOBAL$DIAG_DIR, "/[<sub>/]<target>_{metrics.csv,oob_corrected.csv}\n", sep = "")
cat("QM calibration-> ", GLOBAL$CALIB_DIR, "/[<sub>/]qm_<target>.rds (used by 03_predict_annual_30m_maps.R)\n", sep = "")
cat("========================================\n")

writeLines(c("Final RF and shared OOB QM completed.",
             "Production predictors=fixed_15; VSURF validation selections do not replace this list.",
             "Annual mapping calibration=training_OOB_QM."),file.path(GLOBAL$OUTPUT_DIR,"RUN_COMPLETE.txt"))
