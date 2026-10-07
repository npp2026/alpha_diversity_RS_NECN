# ============================================================
# predict_age_sensitivity.R
#
# Forest-age sensitivity of biodiversity quantiles (q90/q95)
# ------------------------------------------------------------
# 流程:
#   1) 从 ENV.tif 读取栅格预测变量 (DEM + 气候 BIO*)
#   2) 训练点 CSV (WGS84 long/lat) 投影到 ENV CRS 并提取环境变量
#   2b) 对 ENV 提取 NA 用 CSV 自带的 bio{x}_wc / DEMc 做“订正后填补”
#       (只填 DEM + 气候，不填 Forest_age)
#   3) 用 DEM + 气候 + Forest_age 训练 Rich_tree / Shannon_wiener
#      的分位数回归森林 (ranger QRF; 可选 XGBoost reg:quantileerror)
#   4) 固定 Forest_age = 80/100/120 预测 q90/q95 曲面
#   5) 以 Forest_age = 100 为参考输出敏感性差值与百分比变化
#   6) 验证: 重复 random hold-out (快速对照) + kNNDM 空间 CV (主验证)
#      kNNDM predpoints 来自 ENV 预测域的环境分层抽样
#   7) 计算完整 AOA（DEM/BIO + Forest_age）与环境-only AOA（DEM/BIO）
#   8) Forest_age 偏依赖 (PDP) sanity check
#
# 输入 (相对 base_dir):
#   ENV.tif        多波段栅格, 须包含 raster_predictor_vars 列出的图层
#   train4pot.csv  须含 long, lat, Forest_age, response_vars, 以及
#                  bio{x}_wc / DEMc (用于订正填补 ENV 提取 NA)
# 输出:
#   out_dir/ 下的分位面、敏感性面、验证指标、AOA、PDP、诊断 CSV 与 run_info.rds
#
# 依赖: terra, ranger (+ CAST, sf 用于 kNNDM; xgboost 仅在 XGBOOST 模式;
#       FNN 可选, 用于加速最近邻诊断)
# 硬件目标: 约 32 GB RAM / 8 核 (见下方 CONFIG, 可自动探测核数并回退)
#
# 配对脚本: tune_qrf_params.R 先调参, 本脚本可自动 source 其输出的
#           apply_best_qrf_params_to_final_mapping.R (auto_load_tuned_params)。
# ============================================================

suppressPackageStartupMessages({
  library(terra)
})

# ============================================================
# 0. 参数设置
# ============================================================

# 唯一需要按机器修改的根路径; 其余路径都从它派生。
# Portable paths: keep input/output paths relative to an explicit data root.
# Shared input/cache/output contracts are resolved relative to this installed script.
.fem_script<-sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1])
source(file.path(dirname(normalizePath(.fem_script)),"..","R","qrf_contracts.R"))
.fem_software_version <- trimws(readLines(file.path(dirname(normalizePath(.fem_script)), "..", "VERSION"), n = 1L, warn = FALSE))
rm(.fem_script)

base_dir <- Sys.getenv("FEM_POTENTIAL_DATA_DIR", unset = getwd())
base_dir <- normalizePath(path.expand(base_dir), winslash = "/", mustWork = TRUE)
setwd(base_dir)

env_file <- "ENV.tif"
csv_file <- "train4pot.csv"
template_file <- "templ_1km.tif"
# TRUE: templ_1km.tif 中 0 与 NA 都视为模板外；若 0 是有效像元值，请改为 FALSE。
mask_zero_as_na <- TRUE

out_dir <- "quantile_sensitivity_outputs"
fem_qrf_new_output(out_dir)

# ENV.tif 只读；DEM/BIO 会派生为共享的 templ_1km 对齐+mask 缓存，供调参/预测/oldage 复用。
aligned_env_file <- file.path(base_dir, "ENV_predictors_aligned_to_templ_1km_masked.tif")
# 复用已对齐缓存；缓存会用 ENV.tif 与 templ_1km.tif 的文件指纹校验，源文件变化会自动重建。
reuse_aligned_env_cache <- TRUE

# terra 临时目录与内存策略
# 32G 内存机器建议给 terra 使用约 60-70% RAM，其余留给 ranger / R 对象 / 系统。
terra_tempdir <- file.path(base_dir, "terra_tmp")
dir.create(terra_tempdir, showWarnings = FALSE, recursive = TRUE)

# tempdir / memfrac / progress 各版本都支持，单独设置，确保即使下面的
# memmax 不被旧版 terra 接受，这些关键项也不会被一起丢弃。
terra::terraOptions(
  tempdir = terra_tempdir,
  memfrac = 0.70,
  progress = 0
)
# memmax (GB 上限) 仅较新版 terra 支持；老版本忽略即可。
try(terra::terraOptions(memmax = 24), silent = TRUE)

# 模型类型：
# "RANGER_QRF" = ranger quantile regression forest，推荐
# "XGBOOST"    = xgboost native quantile objective
model_type <- "RANGER_QRF"

# 响应变量
response_vars <- c("Rich_tree", "Shannon_wiener")

# ENV.tif 中必须存在的空间预测变量
climate_vars <- c("BIO6", "BIO10", "BIO12", "BIO17", "BIO15", "BIO4")
raster_predictor_vars <- c("DEM", climate_vars)

# 训练/预测模型使用的全部变量
predictor_vars <- c(raster_predictor_vars, "Forest_age")

# 预测分位数
quantiles <- c(0.90, 0.95)
quantiles <- sort(unique(quantiles))
q_names <- paste0("q", round(quantiles * 100))

# Forest age 敏感性分析固定值
age_values <- c(80, 100, 120)
ref_age <- 100

# 结果后处理
nonnegative_responses <- c("Rich_tree", "Shannon_wiener")
clip_negative_predictions <- TRUE

# 适用域/外推检查
# range mask 是一维范围筛查：可识别单变量超范围，但不能识别多变量组合外推；
# 若需要严格适用域，建议使用 CAST::aoa()。
make_range_mask <- TRUE
apply_range_mask_to_predictions <- FALSE
# 如果希望外推区直接置 NA，改成 TRUE；建议先 FALSE，输出 mask 后自行检查

# CAST Area of Applicability (AOA)：多变量适用域，比一维 range mask 更严格。
# 每个响应变量、每个固定 Forest_age 单独输出 AOA 和 DI 栅格。
make_aoa <- TRUE
apply_aoa_mask_to_predictions <- FALSE
aoa_use_variable_importance <- TRUE
aoa_use_knndm_folds <- TRUE
aoa_lpd <- FALSE
aoa_method <- "L2"
aoa_skip_existing <- FALSE # 已有 AOA/DI tif 不直接复用；确保 AOA 基于 templ_1km 对齐/mask 后的 ENV 重算
# TRUE = 每个 Forest_age 分别计算 AOA/DI，最严格但较慢；FALSE = 只算 ref_age 并复用，速度更快。
aoa_compute_per_age <- TRUE

# 额外计算“环境变量-only AOA”：只使用 DEM + BIO 气候/地形变量，
# 不包含 Forest_age。因此该 AOA 不随 80/100/120 年龄情景变化。
# 推荐保留 TRUE，用于区分“基础环境外推风险”和“林龄-环境组合外推风险”。
make_environment_only_aoa <- TRUE
# FALSE = 环境变量等权，输出一套共同的环境 AOA；TRUE = 使用参考响应变量的 permutation importance 加权。
environment_only_aoa_use_variable_importance <- FALSE
# 环境-only AOA 使用哪个响应变量的训练行和 kNNDM folds。当前两个响应变量训练行一致，默认第一个即可。
environment_only_aoa_reference_response <- response_vars[1]

# 树数稳定性检查：只在每个响应变量的“最优 mtry/node/always_split”上比较，
# 不把 800-2000 乘到完整参数网格里，避免计算量爆炸。
run_tree_count_sensitivity <- TRUE
qrf_num_trees_final_candidates <- c(800, 1200, 1500, 2000)
# TRUE: allow the tuning snippet to overwrite validation/final tree counts.
# FALSE: keep tree counts controlled by this prediction script, avoiding silent source() overrides.
allow_tuned_tree_counts <- FALSE

# 快速预测模式：一次性读取 ENV 预测矩阵，适合当前 1655 x 1077 规模
# 如果内存不足，改成 FALSE 使用 terra::predict 分块预测
use_fast_matrix_prediction <- TRUE

# 验证设置
# random hold-out 作为快速对照；kNNDM CV 作为主空间验证
run_holdout_validation <- TRUE
validation_fraction <- 0.20
validation_repeats <- 5

run_knndm_validation <- TRUE
# kNNDM 折数。全局验证 (堆叠所有折外预测) 下,总测试样本数恒为 N,
# 与 k 无关,因此覆盖率/pinball 的点估计精度不随 k 变化。提高 k 的主要好处是
# 每折训练集更大 (10 折时用 90%, 5 折时用 80%),更接近最终全量模型,
# 验证更不悲观。代价是计算量约翻倍。判优依据是 folds$W (越小越匹配),
# 见输出 validation_knndm_fold_assignment_info.csv;若 k=10 的 W 不比 k=5 差即可用。
knndm_k <- 10

# kNNDM predpoints 抽样：推荐使用预测域 ENV 的环境分层抽样
# 逻辑：ENV 有效像元 -> 标准化 -> PCA -> k-means 环境簇 -> 按簇分层抽样 -> CAST::knndm
knndm_predpoints_sampling <- "environmental_stratified"
knndm_samplesize <- 10000
knndm_env_clusters <- 50
knndm_min_per_cluster <- 50
knndm_pca_components <- 4
knndm_candidate_multiplier <- 3
knndm_kmeans_nstart <- 5
knndm_kmeans_iter_max <- 100

# CAST::knndm fold 生成参数
# clustering: "kmeans"（默认，点数大时更快）或 "hierarchical"（对重复坐标稳健）。
# 若训练点含大量重复坐标（如同一样地多年重复测量），kmeans 可能报
# "more cluster centers than distinct data points"；脚本会自动回退到 hierarchical。
# 也可直接把下行设为 "hierarchical" 跳过这次回退。
knndm_clustering <- "kmeans"
knndm_maxp <- 0.5
# kNNDM 默认在地理空间匹配最近邻距离；环境分层只用于 predpoints 的代表性抽样。
# 注意：本脚本的 predpoints 是 sf 几何点，因此这里只实现 space = "geographical"。
# 若要做 feature-space kNNDM，需要另行构建同名特征 data.frame 传给 tpoints/predpoints。
knndm_space <- "geographical"

# 输出 kNNDM predpoints 采样点与诊断
write_knndm_predpoints_gpkg <- TRUE

min_train_rows <- 30

# 随机种子与线程
seed <- 123
set.seed(seed)

t_start <- Sys.time()

# ============================================================
# 32G RAM / 8-core CPU 并行策略
# ============================================================
# 默认认为当前机器是专用计算环境，尽量用满 8 核。
# 如果你还要同时做其它工作，可把 reserve_cores_for_os 改为 1。
hardware_cores <- fem_qrf_threads(8L)
hardware_memory_gb <- 32
reserve_cores_for_os <- 0

detected_cores <- parallel::detectCores(logical = TRUE)
if (is.na(detected_cores)) detected_cores <- hardware_cores

n_model_threads <- max(1, min(hardware_cores, detected_cores) - reserve_cores_for_os)
n_cores <- n_model_threads

# 避免 ranger/xgboost 多线程时再被 BLAS/OpenMP 额外嵌套放大。
# ranger/xgboost 的线程数由 num.threads / nthread 显式控制。
Sys.setenv(
  OMP_NUM_THREADS = as.character(n_model_threads),
  MKL_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  NUMEXPR_NUM_THREADS = "1",
  GDAL_NUM_THREADS = "ALL_CPUS"
)

# ranger 参数
# 分位数森林估计 q90/q95 这类尾部分位,比估计均值需要更多树才稳定,
# v1.6 MS profile uses 1500 trees in both validation and final mapping.
# 提示: 800 -> 1500 通常只让尾部分位面更平滑,增益有限;
# 可对比 validation 指标在 500 vs 1500 下是否明显变化来判断是否值得。
qrf_num_trees_validation <- 1500L
qrf_num_trees_final <- 1500
# Chunk size for full-raster ranger quantile prediction. Lower it if RAM is tight.
qrf_prediction_chunk_size <- 200000L
qrf_min_node_size <- 5
# 默认使用 sqrt(p)，保留随机森林的变量子采样，避免退化成 bagging；
# 当前 predictor_vars = DEM + 6 个气候变量 + Forest_age，共 8 个变量，因此默认 mtry = 2。
qrf_mtry <- max(2, floor(sqrt(length(predictor_vars))))
# Forest_age 是敏感性分析的核心变量；让它在每个 split 中额外进入候选集，
# 同时保留其他变量的 mtry 随机子采样。若想对比，可设为 NULL。
qrf_always_split_variables <- "Forest_age"

# 可选：自动载入 tune_qrf_params.R 生成的最优参数片段
# (apply_best_qrf_params_to_final_mapping.R)，避免人工复制粘贴出错。
# 新版片段优先提供 qrf_params_by_response，使 Rich_tree / Shannon_wiener
# 可以分别使用各自最优的 mtry / min.node.size / always.split.variables。
# 若只找到旧版标量参数，本脚本会自动复制成每个响应变量共用的参数。
# 设为 FALSE 则始终使用本文件中的默认值。
# v1.6: SI S5.4 parameters are the default; a fresh tuning result is opt-in.
qrf_parameter_mode <- Sys.getenv("FEM_QRF_PARAMETER_MODE", "manuscript")
if (!qrf_parameter_mode %in% c("manuscript", "retuned")) stop("Invalid FEM_QRF_PARAMETER_MODE")
auto_load_tuned_params <- qrf_parameter_mode == "retuned"
tuned_params_file <- file.path(
  base_dir, "qrf_parameter_tuning_outputs",
  "apply_best_qrf_params_to_final_mapping.R"
)
# Preserve local tree-count settings before source(), because older snippets may contain
# qrf_num_trees_validation/final and qrf_num_trees_final_candidates assignments.
local_qrf_num_trees_validation <- qrf_num_trees_validation
local_qrf_num_trees_final <- qrf_num_trees_final
local_qrf_num_trees_final_candidates <- qrf_num_trees_final_candidates
if (isTRUE(auto_load_tuned_params) && file.exists(tuned_params_file)) {
  message("Loading tuned QRF parameters from: ", tuned_params_file)
  source(tuned_params_file, local = FALSE)
} else if (isTRUE(auto_load_tuned_params)) {
  stop("Retuned profile requires parameter snippet: ", tuned_params_file)
}

if (!isTRUE(allow_tuned_tree_counts)) {
  qrf_num_trees_validation <- local_qrf_num_trees_validation
  qrf_num_trees_final <- local_qrf_num_trees_final
  qrf_num_trees_final_candidates <- local_qrf_num_trees_final_candidates
}

if (qrf_parameter_mode == "manuscript") {
  qrf_params_by_response <- list(
    Rich_tree = list(qrf_mtry=4L, qrf_min_node_size=10L,
      qrf_always_split_variables="Forest_age", qrf_num_trees_validation=1500L, qrf_num_trees_final=1500L),
    Shannon_wiener = list(qrf_mtry=2L, qrf_min_node_size=10L,
      qrf_always_split_variables=NULL, qrf_num_trees_validation=1500L, qrf_num_trees_final=1500L))
}

# 如果调参片段没有提供 qrf_params_by_response，则使用当前标量参数作为所有响应变量的 fallback。
if (!exists("qrf_params_by_response")) {
  qrf_params_by_response <- setNames(
    lapply(response_vars, function(resp) {
      list(
        qrf_mtry = qrf_mtry,
        qrf_min_node_size = qrf_min_node_size,
        qrf_always_split_variables = qrf_always_split_variables,
        qrf_num_trees_validation = qrf_num_trees_validation,
        qrf_num_trees_final = qrf_num_trees_final
      )
    }),
    response_vars
  )
}

# 补齐列表中缺失的字段，避免手工改片段时遗漏。
for (resp in response_vars) {
  if (is.null(qrf_params_by_response[[resp]])) qrf_params_by_response[[resp]] <- list()
  if (is.null(qrf_params_by_response[[resp]]$qrf_mtry)) qrf_params_by_response[[resp]]$qrf_mtry <- qrf_mtry
  if (is.null(qrf_params_by_response[[resp]]$qrf_min_node_size)) qrf_params_by_response[[resp]]$qrf_min_node_size <- qrf_min_node_size
  if (!("qrf_always_split_variables" %in% names(qrf_params_by_response[[resp]]))) {
    qrf_params_by_response[[resp]]$qrf_always_split_variables <- qrf_always_split_variables
  }
  if (!isTRUE(allow_tuned_tree_counts)) {
    qrf_params_by_response[[resp]]$qrf_num_trees_validation <- qrf_num_trees_validation
    qrf_params_by_response[[resp]]$qrf_num_trees_final <- qrf_num_trees_final
  } else {
    if (is.null(qrf_params_by_response[[resp]]$qrf_num_trees_validation)) qrf_params_by_response[[resp]]$qrf_num_trees_validation <- qrf_num_trees_validation
    if (is.null(qrf_params_by_response[[resp]]$qrf_num_trees_final)) qrf_params_by_response[[resp]]$qrf_num_trees_final <- qrf_num_trees_final
  }
}

message("QRF parameters by response:")
for (resp in response_vars) {
  pp <- qrf_params_by_response[[resp]]
  message("  ", resp, ": mtry=", pp$qrf_mtry,
          ", min.node.size=", pp$qrf_min_node_size,
          ", always.split.variables=",
          if (is.null(pp$qrf_always_split_variables)) "NULL" else paste(pp$qrf_always_split_variables, collapse = ","),
          ", trees(validation)=", pp$qrf_num_trees_validation,
          ", trees(final)=", pp$qrf_num_trees_final)
}

# Forest_age 响应 sanity check：输出 PDP 曲线，帮助判断 80-120 年敏感性是否被模型学到。
run_forest_age_pdp <- TRUE
pdp_age_grid <- sort(unique(c(seq(20, 160, by = 5), age_values)))
pdp_sample_size <- 2000
pdp_background <- "training"  # 当前使用训练样点环境背景；不要用响应变量分层。

# XGBoost 参数
xgb_nrounds <- 500
xgb_eta <- 0.03
xgb_max_depth <- 4
xgb_min_child_weight <- 3
xgb_subsample <- 0.8
xgb_colsample_bytree <- 0.8

# 写出选项
write_datatype <- "FLT4S"
# GTiff 写出优化：
# - NUM_THREADS=ALL_CPUS：GDAL 写压缩 GeoTIFF 时尽量多线程；
# - PREDICTOR=3：对 FLT4S 浮点栅格的 LZW 压缩通常更有效。
write_gdal <- c("COMPRESS=LZW", "TILED=YES", "BIGTIFF=IF_SAFER", "NUM_THREADS=ALL_CPUS", "PREDICTOR=3")

# ============================================================
# 1. 基础检查与包加载
# ============================================================

if (!(ref_age %in% age_values)) {
  stop("ref_age 必须包含在 age_values 中。当前 ref_age = ", ref_age)
}

if (model_type == "RANGER_QRF") {
  if (!requireNamespace("ranger", quietly = TRUE)) {
    stop("需要安装 ranger 包：install.packages('ranger')")
  }
}

if (run_knndm_validation || make_aoa || make_environment_only_aoa) {
  if (!requireNamespace("CAST", quietly = TRUE)) {
    stop("kNNDM CV / AOA 需要 CAST 包：install.packages('CAST')；或设置 run_knndm_validation <- FALSE 且 make_aoa <- FALSE 且 make_environment_only_aoa <- FALSE")
  }
  if (run_knndm_validation && !requireNamespace("sf", quietly = TRUE)) {
    stop("kNNDM CV 需要 sf 包：install.packages('sf')；或设置 run_knndm_validation <- FALSE")
  }
}

if (model_type == "XGBOOST") {
  if (!requireNamespace("xgboost", quietly = TRUE)) {
    stop("需要安装 xgboost 包：install.packages('xgboost')")
  }
  if (utils::packageVersion("xgboost") < "2.0.0") {
    stop("XGBoost 分位数目标 reg:quantileerror 需要 xgboost >= 2.0.0。请升级 xgboost。")
  }
}

message("Model type: ", model_type)
message("Detected cores: ", detected_cores)
message("Hardware target cores: ", hardware_cores, "; model threads used: ", n_model_threads)
message("Hardware target RAM: ", hardware_memory_gb, " GB; terra memmax target: 24 GB")
message("Predictors: ", paste(predictor_vars, collapse = ", "))
message("Quantiles: ", paste(quantiles, collapse = ", "))


# ============================================================
# 2. 工具函数
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
      object_name, " 缺少以下字段/图层：",
      paste(miss, collapse = ", "),
      "\n当前可用名称：", paste(available, collapse = ", ")
    )
  }
}


# ============================================================
# QRF response-specific parameter helpers
# ============================================================

normalize_qrf_param_list <- function(pp) {
  if (is.null(pp)) pp <- list()
  if (is.null(pp$qrf_mtry)) pp$qrf_mtry <- qrf_mtry
  if (is.null(pp$qrf_min_node_size)) pp$qrf_min_node_size <- qrf_min_node_size
  if (!("qrf_always_split_variables" %in% names(pp))) {
    pp$qrf_always_split_variables <- qrf_always_split_variables
  }
  if (!isTRUE(allow_tuned_tree_counts)) {
    pp$qrf_num_trees_validation <- qrf_num_trees_validation
    pp$qrf_num_trees_final <- qrf_num_trees_final
  } else {
    if (is.null(pp$qrf_num_trees_validation)) {
      pp$qrf_num_trees_validation <- qrf_num_trees_validation
    }
    if (is.null(pp$qrf_num_trees_final)) {
      pp$qrf_num_trees_final <- qrf_num_trees_final
    }
  }
  pp$qrf_mtry <- as.integer(pp$qrf_mtry)
  pp$qrf_min_node_size <- as.integer(pp$qrf_min_node_size)
  pp$qrf_num_trees_validation <- as.integer(pp$qrf_num_trees_validation)
  pp$qrf_num_trees_final <- as.integer(pp$qrf_num_trees_final)
  if (!is.null(pp$qrf_always_split_variables) && length(pp$qrf_always_split_variables) == 0) {
    pp$qrf_always_split_variables <- NULL
  }
  pp
}

get_qrf_params <- function(response_var = NULL) {
  if (is.null(response_var) || is.na(response_var) || !(response_var %in% names(qrf_params_by_response))) {
    return(normalize_qrf_param_list(list(
      qrf_mtry = qrf_mtry,
      qrf_min_node_size = qrf_min_node_size,
      qrf_always_split_variables = qrf_always_split_variables,
      qrf_num_trees_validation = qrf_num_trees_validation,
      qrf_num_trees_final = qrf_num_trees_final
    )))
  }
  normalize_qrf_param_list(qrf_params_by_response[[response_var]])
}

get_qrf_param_value <- function(response_var, field) {
  pp <- get_qrf_params(response_var)
  if (!(field %in% names(pp))) {
    stop("QRF 参数字段不存在: ", field, "；response = ", response_var)
  }
  pp[[field]]
}

qrf_params_to_df <- function() {
  rows <- lapply(response_vars, function(resp) {
    pp <- get_qrf_params(resp)
    data.frame(
      response = resp,
      qrf_mtry = pp$qrf_mtry,
      qrf_min_node_size = pp$qrf_min_node_size,
      qrf_always_split_variables = if (is.null(pp$qrf_always_split_variables)) "NULL" else paste(pp$qrf_always_split_variables, collapse = ";"),
      qrf_num_trees_validation = pp$qrf_num_trees_validation,
      qrf_num_trees_final = pp$qrf_num_trees_final,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

# 参数一致性检查
for (resp in response_vars) {
  pp <- get_qrf_params(resp)
  if (!is.null(pp$qrf_always_split_variables) && length(pp$qrf_always_split_variables) > 0) {
    stop_if_missing(pp$qrf_always_split_variables, predictor_vars, paste0("qrf_always_split_variables for ", resp))
  }
}
write.csv(qrf_params_to_df(), file.path(out_dir, "qrf_final_parameters_by_response.csv"), row.names = FALSE)

if (run_knndm_validation && !identical(knndm_space, "geographical")) {
  stop(
    "当前脚本的 kNNDM predpoints 是 sf 几何点，仅支持 knndm_space = 'geographical'。",
    "如需 feature-space kNNDM，请改为传入训练/预测特征 data.frame，且列名必须一致。"
  )
}

if (run_forest_age_pdp && !identical(pdp_background, "training")) {
  stop(
    "当前脚本仅实现 pdp_background = 'training'。",
    "如需 prediction-domain PDP，请用 knndm_predpoints_obj$sample_values 或 fast_data$base_df 另建背景样本。"
  )
}


pinball_loss <- function(y, pred, q) {
  err <- y - pred
  mean(ifelse(err >= 0, q * err, (q - 1) * err), na.rm = TRUE)
}

enforce_quantile_order_df <- function(p, quantiles) {
  ord <- order(quantiles)
  p <- p[, ord, drop = FALSE]
  names(p) <- paste0("q", round(sort(quantiles) * 100))

  # 防止 quantiles 只有一个值时 2:ncol(p) 变成 2:1，导致访问 p[[2]] 报错
  if (ncol(p) <= 1) {
    return(p)
  }

  for (i in 2:ncol(p)) {
    p[[i]] <- pmax(p[[i]], p[[i - 1]], na.rm = FALSE)
  }

  return(p)
}

enforce_quantile_order_raster <- function(r) {
  if (nlyr(r) <= 1) return(r)
  for (i in 2:nlyr(r)) {
    r[[i]] <- max(r[[i]], r[[i - 1]], na.rm = FALSE)
  }
  return(r)
}

clip_negative_df <- function(p) {
  for (nm in names(p)) {
    p[[nm]] <- pmax(p[[nm]], 0, na.rm = FALSE)
  }
  return(p)
}

clip_negative_raster <- function(r) {
  r <- ifel(r < 0, 0, r)
  return(r)
}

postprocess_prediction_df <- function(p, response_var) {
  p <- enforce_quantile_order_df(p, quantiles)
  if (clip_negative_predictions && response_var %in% nonnegative_responses) {
    p <- clip_negative_df(p)
  }
  return(p)
}

postprocess_prediction_raster <- function(r, response_var) {
  names(r) <- q_names
  r <- enforce_quantile_order_raster(r)
  if (clip_negative_predictions && response_var %in% nonnegative_responses) {
    r <- clip_negative_raster(r)
  }
  names(r) <- q_names
  return(r)
}

write_raster_safe <- function(r, filename) {
  writeRaster(
    r,
    filename,
    overwrite = TRUE,
    datatype = write_datatype,
    gdal = write_gdal
  )
  return(rast(filename))
}

# ============================================================
# templ_1km 强制对齐工具
#   - 原始 ENV.tif 只读，不覆盖、不改写。
#   - 生成派生缓存 ENV_predictors_aligned_to_templ_1km_masked.tif。
#   - 训练点提取使用原始 ENV.tif；整图预测、kNNDM predpoints、AOA/DI 使用对齐并 mask 的 ENV_pred。
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
# 3. 读取 ENV raster，并检查图层名
# ============================================================

ENV_raw <- rast(env_file)

message("==================================================")
message("Raw ENV raster (read-only; will NOT be overwritten):")
print(ENV_raw)
message("Raw ENV names: ", paste(names(ENV_raw), collapse = ", "))

stop_if_missing(raster_predictor_vars, names(ENV_raw), "ENV.tif")

ENV_model_pred <- ENV_raw[[raster_predictor_vars]]
ENV_model_pred <- terra::ifel(is.finite(ENV_model_pred)&ENV_model_pred != -9999,ENV_model_pred,NA)  # 训练点环境提取使用原始 ENV.tif，不受 templ_1km mask 限制

template_1km <- read_templ_1km_mask(template_file)
message("templ_1km mask:")
print(template_1km)

ENV_pred <- align_env_to_templ_1km_mask(
  ENV_model_pred,
  template_1km,
  aligned_file = aligned_env_file,
  reuse_cache = reuse_aligned_env_cache,
  source_files = c(env_file, template_file)
)
# 重要数据流：
#   - 训练点提取 / 模型训练：使用 ENV_model_pred（原始 ENV.tif，保留模板外样点）。
#   - 整图预测 / kNNDM predpoints / AOA/DI：使用 ENV_pred（对齐到 templ_1km 并 masked）。

# 检查 templ_1km 对齐并 mask 后的 ENV 预测变量内部 NA 概况
na_env <- as.data.frame(global(is.na(ENV_pred), "sum", na.rm = TRUE))
na_env$layer <- rownames(na_env)
rownames(na_env) <- NULL
names(na_env)[1] <- "NA_count"
write.csv(
  na_env[, c("layer", "NA_count")],
  file.path(out_dir, "ENV_predictor_NA_summary.csv"),
  row.names = FALSE
)
write.csv(
  na_env[, c("layer", "NA_count")],
  file.path(out_dir, "ENV_predictor_aligned_to_templ_1km_masked_NA_summary.csv"),
  row.names = FALSE
)

# ============================================================
# 4. 读取 CSV，处理字段类型和坐标缺失
# ============================================================

dat_raw <- read.csv(csv_file, stringsAsFactors = FALSE)

message("==================================================")
message("CSV columns: ", paste(names(dat_raw), collapse = ", "))
message("Original CSV rows: ", nrow(dat_raw))

required_cols <- c("long", "lat", "Forest_age", response_vars)
stop_if_missing(required_cols, names(dat_raw), "CSV")

# 转 numeric，避免 CSV 读成字符
numeric_cols <- required_cols
for (cc in numeric_cols) {
  dat_raw[[cc]] <- safe_as_numeric(dat_raw[[cc]])
}

# 坐标缺失不能提取栅格，先剔除
coord_ok <- !is.na(dat_raw$long) & !is.na(dat_raw$lat) &
  is.finite(dat_raw$long) & is.finite(dat_raw$lat)

dat0 <- dat_raw[coord_ok, , drop = FALSE]
message("Rows after removing missing coordinates: ", nrow(dat0))

if (nrow(dat0) == 0) {
  stop("坐标有效样点数为 0，请检查 long/lat 字段。")
}

# 可选坐标范围提示
if (any(dat0$long < -180 | dat0$long > 180 | dat0$lat < -90 | dat0$lat > 90, na.rm = TRUE)) {
  warning("发现超出常规经纬度范围的 long/lat，请确认 CSV 坐标确实为 WGS84 经纬度。")
}

# ============================================================
# 5. CSV 经纬度点投影到 ENV CRS，并提取 DEM + climate
# ============================================================

message("==================================================")
message("Projecting CSV points and extracting ENV predictors for model training...")

pts_ll <- vect(
  dat0,
  geom = c("long", "lat"),
  crs = "EPSG:4326"
)

# 训练点环境变量提取使用原始 ENV.tif，不用 templ_1km mask 剔除样点。
pts_model_env <- project(pts_ll, crs(ENV_model_pred))

# x_env/y_env 用于空间验证/kNNDM，应与预测域 ENV_pred 的 CRS 保持一致；但不按模板 mask 过滤训练点。
pts_pred_crs <- project(pts_ll, crs(ENV_pred))
xy_env <- crds(pts_pred_crs)
dat0$x_env <- xy_env[, 1]
dat0$y_env <- xy_env[, 2]

env_extract <- extract(
  ENV_model_pred,
  pts_model_env,
  ID = FALSE
)

# 如果 CSV 原本已有同名列，先移除，避免 cbind 后重名
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

write.csv(
  extract_na_summary,
  file.path(out_dir, "extracted_ENV_predictor_NA_summary.csv"),
  row.names = FALSE
)

print(extract_na_summary)

# ============================================================
# 6. 构建每个响应变量的训练数据
# ============================================================

make_training_data <- function(dat, response_var) {
  # x_env/y_env 不进入模型，仅用于空间结构验证，如 kNNDM CV
  use_cols <- c(predictor_vars, "x_env", "y_env", response_var)
  d <- dat[, use_cols, drop = FALSE]

  # 全部转 numeric
  for (cc in names(d)) {
    d[[cc]] <- safe_as_numeric(d[[cc]])
  }

  ok <- complete.cases(d)
  for (cc in names(d)) {
    ok <- ok & is.finite(d[[cc]])
  }

  d <- d[ok, , drop = FALSE]
  names(d)[names(d) == response_var] <- "response"
  return(d)
}

train_list <- list()

for (resp in response_vars) {
  d <- make_training_data(dat0, resp)
  message("--------------------------------------------------")
  message("Response: ", resp)
  message("Training rows after removing NA: ", nrow(d))

  if (nrow(d) < min_train_rows) {
    stop("有效训练样本过少：", resp, "；当前 n = ", nrow(d))
  }
  train_list[[resp]] <- d
}

# 每个响应变量单独做 age 外推检查
for (resp in response_vars) {
  d <- train_list[[resp]]
  age_min <- min(d$Forest_age, na.rm = TRUE)
  age_max <- max(d$Forest_age, na.rm = TRUE)
  if (any(age_values < age_min | age_values > age_max)) {
    warning(
      resp, ": 部分 age_values 超出该响应变量训练样本 Forest_age 范围 [",
      age_min, ", ", age_max, "]。RF/QRF 不能可靠外推。"
    )
  }
}

# ============================================================
# 6b. 训练点空间聚集诊断
#     - 多少样点落在同一个 1km 像元 (坐标重合程度)
#     - 最近邻距离分布 (空间聚集尺度的参考)
#     坐标为 ENV CRS (Albers, 单位米)。点数很大时对最近邻用随机子样,避免 O(n^2) 内存。
# ============================================================

diagnose_point_clustering <- function(coords_xy, env_pred, out_dir,
                                      label = "training",
                                      max_full_nn = 12000,
                                      seed = 123) {
  xy <- as.matrix(coords_xy)
  xy <- xy[stats::complete.cases(xy), , drop = FALSE]
  n <- nrow(xy)
  if (n < 2) {
    warning("点数过少，跳过空间聚集诊断。")
    return(invisible(NULL))
  }

  res_m <- terra::res(env_pred)[1]

  # --- 1km 像元重合 ---
  cells <- terra::cellFromXY(env_pred[[1]], xy)
  cells <- cells[!is.na(cells)]
  cell_tab <- table(cells)
  pts_sharing <- sum(cell_tab[cell_tab > 1])

  coincidence <- data.frame(
    dataset = label,
    n_points = n,
    n_points_in_raster = length(cells),
    n_unique_cells = length(cell_tab),
    n_cells_with_multiple = sum(cell_tab > 1),
    n_points_sharing_a_cell = pts_sharing,
    pct_points_sharing = round(100 * pts_sharing / n, 2),
    max_points_per_cell = as.integer(max(cell_tab)),
    raster_res_m = res_m
  )

  # --- 最近邻距离 (米) ---
  # 优先用 FNN 做精确最近邻 (O(n log n) k-d tree)，对全部点直接计算，
  # 避免 stats::dist() 的 O(n^2) 稠密矩阵 (n=12000 时约 1.1 GB)。
  # 没装 FNN 时回退到“子样 + 稠密矩阵”，行为与旧版一致。
  if (requireNamespace("FNN", quietly = TRUE)) {
    xy_nn <- xy
    nn_note <- paste0("all ", n, " (FNN k-d tree)")
    nn_dist <- as.numeric(FNN::knn.dist(xy_nn, k = 1)[, 1])
  } else {
    if (n > max_full_nn) {
      set.seed(seed)
      xy_nn <- xy[safe_sample(seq_len(n), max_full_nn), , drop = FALSE]
      nn_note <- paste0("subsample ", max_full_nn, " / ", n, " (dense dist)")
    } else {
      xy_nn <- xy
      nn_note <- paste0("all ", n, " (dense dist)")
    }
    dmat <- as.matrix(stats::dist(xy_nn))
    diag(dmat) <- NA
    nn_dist <- apply(dmat, 1, min, na.rm = TRUE)
  }
  qs <- stats::quantile(nn_dist, c(0, .05, .25, .5, .75, .95, 1), names = FALSE)

  nn_summary <- data.frame(
    dataset = label,
    note = nn_note,
    n_used = nrow(xy_nn),
    nn_min_m = round(qs[1], 1),
    nn_q05_m = round(qs[2], 1),
    nn_q25_m = round(qs[3], 1),
    nn_median_m = round(qs[4], 1),
    nn_q75_m = round(qs[5], 1),
    nn_q95_m = round(qs[6], 1),
    nn_max_m = round(qs[7], 1),
    nn_mean_m = round(mean(nn_dist), 1),
    pct_nn_within_1cell = round(100 * mean(nn_dist <= res_m), 2)
  )

  write.csv(coincidence,
            file.path(out_dir, paste0("point_coincidence_", label, ".csv")),
            row.names = FALSE)
  write.csv(nn_summary,
            file.path(out_dir, paste0("nearest_neighbor_distance_summary_", label, ".csv")),
            row.names = FALSE)
  write.csv(data.frame(nn_distance_m = round(nn_dist, 1)),
            file.path(out_dir, paste0("nearest_neighbor_distances_", label, ".csv")),
            row.names = FALSE)

  message("空间聚集诊断 (", label, "):")
  message("  ", n, " 点落在 ", coincidence$n_unique_cells, " 个 1km 像元; ",
          coincidence$n_cells_with_multiple, " 个像元含多点; ",
          pts_sharing, " 点 (", coincidence$pct_points_sharing, "%) 与他点同格; ",
          "单格最多 ", coincidence$max_points_per_cell, " 点。")
  message("  最近邻距离(m): 中位 ", nn_summary$nn_median_m,
          ", 5% ", nn_summary$nn_q05_m,
          ", 95% ", nn_summary$nn_q95_m,
          "; ", nn_summary$pct_nn_within_1cell, "% 的点最近邻在 1 像元内 (", nn_note, ")。")

  invisible(list(coincidence = coincidence, nn_summary = nn_summary, nn_dist = nn_dist))
}

# 两个响应的训练坐标几乎一致,用第一个响应的训练点做一次诊断即可代表建模/验证点集。
diag_coords <- train_list[[response_vars[1]]][, c("x_env", "y_env")]
invisible(diagnose_point_clustering(
  coords_xy = diag_coords,
  env_pred = ENV_pred,
  out_dir = out_dir,
  label = "training",
  seed = seed
))

# ============================================================
# 7. 训练数据摘要
# ============================================================

summary_rows <- list()
for (resp in response_vars) {
  d <- train_list[[resp]]
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

summary_table <- do.call(rbind, summary_rows)
write.csv(
  summary_table,
  file.path(out_dir, "model_training_summary.csv"),
  row.names = FALSE
)
print(summary_table)

# ============================================================
# 8. 简单一维 range mask（宽松筛查，不等同于多维 AOA）
# ============================================================

make_univariate_range_mask <- function(env_pred, train_dat, filename) {
  message("Creating univariate range applicability mask: ", basename(filename))
  message("Note: this mask checks one variable at a time; it does not detect multivariate extrapolation.")

  m <- !is.na(env_pred[[1]])

  range_table <- data.frame(
    variable = raster_predictor_vars,
    train_min = NA_real_,
    train_max = NA_real_
  )

  for (i in seq_along(raster_predictor_vars)) {
    v <- raster_predictor_vars[i]
    mn <- min(train_dat[[v]], na.rm = TRUE)
    mx <- max(train_dat[[v]], na.rm = TRUE)

    range_table$train_min[i] <- mn
    range_table$train_max[i] <- mx

    m <- m & !is.na(env_pred[[v]]) & env_pred[[v]] >= mn & env_pred[[v]] <= mx
  }

  m <- ifel(m, 1, NA)
  names(m) <- "range_applicability"

  write_raster_safe(m, filename)
  return(list(mask = rast(filename), range_table = range_table))
}

range_masks <- list()

if (make_range_mask) {
  for (resp in response_vars) {
    out_mask <- file.path(out_dir, paste0("range_applicability_mask_", resp, ".tif"))
    rr <- make_univariate_range_mask(ENV_pred, train_list[[resp]], out_mask)
    range_masks[[resp]] <- rr$mask

    write.csv(
      rr$range_table,
      file.path(out_dir, paste0("range_applicability_training_ranges_", resp, ".csv")),
      row.names = FALSE
    )
  }
}

# ============================================================
# 9. 模型训练函数
# ============================================================

train_ranger_qrf_model <- function(train_dat,
                                   response_var = NULL,
                                   importance = "permutation",
                                   num_trees = NULL,
                                   model_seed = seed) {
  d <- train_dat[, c(predictor_vars, "response"), drop = FALSE]

  pp <- get_qrf_params(response_var)
  mtry_use <- pp$qrf_mtry
  if (is.null(mtry_use)) {
    mtry_use <- max(2, floor(sqrt(length(predictor_vars))))
  }
  min_node_use <- pp$qrf_min_node_size
  always_split_use <- pp$qrf_always_split_variables
  if (is.null(num_trees)) {
    num_trees <- pp$qrf_num_trees_final
  }

  ranger_args <- list(
    formula = response ~ .,
    data = d,
    num.trees = num_trees,
    mtry = mtry_use,
    min.node.size = min_node_use,
    quantreg = TRUE,
    importance = importance,
    num.threads = n_cores,
    seed = model_seed,
    # Validation uses external hold-out / kNNDM metrics; permutation importance does not
    # require prediction.error, so avoid the extra OOB error aggregation.
    oob.error = FALSE
  )

  if (!is.null(always_split_use) && length(always_split_use) > 0) {
    ranger_args$always.split.variables <- always_split_use
  }

  do.call(ranger::ranger, ranger_args)
}

predict_ranger_qrf_dataframe <- function(model, newdata, quantiles,
                                         chunk_size = qrf_prediction_chunk_size) {
  newdata <- newdata[, predictor_vars, drop = FALSE]
  n <- nrow(newdata)
  nq <- length(quantiles)

  # ranger quantile prediction can materialize an n x num.trees intermediate object.
  # Chunking keeps peak memory low while producing row-wise identical predictions.
  if (n <= chunk_size) {
    p <- predict(
      model,
      data = newdata,
      type = "quantiles",
      quantiles = quantiles,
      num.threads = n_cores
    )$predictions
    if (is.vector(p)) {
      p <- matrix(p, ncol = nq)
    }
  } else {
    p <- matrix(NA_real_, nrow = n, ncol = nq)
    for (s in seq.int(1L, n, by = chunk_size)) {
      e <- min(s + chunk_size - 1L, n)
      message("    QRF predict chunk ", s, "-", e, " / ", n)
      pp <- predict(
        model,
        data = newdata[s:e, , drop = FALSE],
        type = "quantiles",
        quantiles = quantiles,
        num.threads = n_cores
      )$predictions
      if (is.vector(pp)) {
        pp <- matrix(pp, ncol = nq)
      }
      p[s:e, ] <- pp
    }
  }

  p <- as.data.frame(p)
  names(p) <- q_names
  return(p)
}

train_xgb_quantile_models <- function(train_dat,
                                      quantiles,
                                      nrounds = xgb_nrounds,
                                      model_seed = seed) {
  x <- as.matrix(train_dat[, predictor_vars, drop = FALSE])
  y <- train_dat$response

  dtrain <- xgboost::xgb.DMatrix(
    data = x,
    label = y,
    missing = NA
  )

  models <- list()

  for (q in quantiles) {
    q_name <- paste0("q", round(q * 100))
    message("Training XGBoost quantile model: ", q_name)

    models[[q_name]] <- xgboost::xgb.train(
      data = dtrain,
      nrounds = nrounds,
      params = list(
        objective = "reg:quantileerror",
        quantile_alpha = q,
        tree_method = "hist",
        max_depth = xgb_max_depth,
        eta = xgb_eta,
        subsample = xgb_subsample,
        colsample_bytree = xgb_colsample_bytree,
        min_child_weight = xgb_min_child_weight,
        nthread = n_cores,
        seed = model_seed
      ),
      verbose = 0
    )
  }

  return(models)
}

predict_xgb_dataframe <- function(model_list, newdata) {
  newdata <- newdata[, predictor_vars, drop = FALSE]
  x <- as.matrix(newdata)

  out <- matrix(
    NA_real_,
    nrow = nrow(x),
    ncol = length(model_list)
  )

  for (i in seq_along(model_list)) {
    out[, i] <- predict(model_list[[i]], newdata = x)
  }

  out <- as.data.frame(out)
  names(out) <- names(model_list)
  out <- out[, q_names, drop = FALSE]
  return(out)
}

# ============================================================
# 10. 重复 hold-out 验证
#     默认重复 5 次随机 80/20 划分，比单次划分更稳定
# ============================================================

validation_metrics <- data.frame()
validation_metric_rows <- list()

if (run_holdout_validation) {
  message("==================================================")
  message("Running repeated hold-out validation...")

  for (resp in response_vars) {
    d <- train_list[[resp]]
    n <- nrow(d)

    if (n < 50) {
      warning(resp, ": 样本量 < 50，跳过 hold-out 验证。")
      next
    }

    val_n <- max(1, floor(n * validation_fraction))

    for (rep_i in seq_len(validation_repeats)) {
      set.seed(seed + rep_i)
      val_idx <- safe_sample(seq_len(n), size = val_n)

      d_train <- d[-val_idx, , drop = FALSE]
      d_val <- d[val_idx, , drop = FALSE]

      if (model_type == "RANGER_QRF") {
        # 验证阶段不计算 permutation importance，节省时间
        m_val <- train_ranger_qrf_model(d_train, response_var = resp, importance = "none", num_trees = get_qrf_param_value(resp, "qrf_num_trees_validation"), model_seed = seed + rep_i)
        p_val <- predict_ranger_qrf_dataframe(m_val, d_val, quantiles)
      } else if (model_type == "XGBOOST") {
        m_val <- train_xgb_quantile_models(d_train, quantiles, nrounds = max(100, floor(xgb_nrounds * 0.6)), model_seed = seed + rep_i)
        p_val <- predict_xgb_dataframe(m_val, d_val)
      } else {
        stop("model_type 只能是 'RANGER_QRF' 或 'XGBOOST'。")
      }

      p_val <- postprocess_prediction_df(p_val, resp)

      for (i in seq_along(quantiles)) {
        q <- quantiles[i]
        qn <- q_names[i]
        y <- d_val$response
        pred <- p_val[[qn]]

        validation_metric_rows[[length(validation_metric_rows) + 1L]] <- data.frame(
          response = resp,
          model_type = model_type,
          repeat_id = rep_i,
          quantile = q,
          n_train = nrow(d_train),
          n_validation = nrow(d_val),
          pinball_loss = pinball_loss(y, pred, q),
          empirical_coverage = mean(y <= pred, na.rm = TRUE),
          prediction_min = min(pred, na.rm = TRUE),
          prediction_max = max(pred, na.rm = TRUE)
        )
      }
    }
  }

  validation_metrics <- if (length(validation_metric_rows) > 0) do.call(rbind, validation_metric_rows) else data.frame()

  write.csv(
    validation_metrics,
    file.path(out_dir, "validation_metrics_repeated_holdout.csv"),
    row.names = FALSE
  )

  if (nrow(validation_metrics) > 0) {
    validation_summary <- aggregate(
      cbind(pinball_loss, empirical_coverage) ~ response + model_type + quantile,
      data = validation_metrics,
      FUN = function(z) c(mean = mean(z, na.rm = TRUE), sd = stats::sd(z, na.rm = TRUE))
    )

    # 展开 aggregate 产生的矩阵列
    validation_summary <- do.call(
      data.frame,
      validation_summary
    )

    names(validation_summary) <- gsub("\\.", "_", names(validation_summary))

    write.csv(
      validation_summary,
      file.path(out_dir, "validation_metrics_repeated_holdout_summary.csv"),
      row.names = FALSE
    )

    print(validation_summary)
  }
}


# ============================================================
# 11. kNNDM CV 空间验证
#     主验证：fold 仍分配给训练样点；但 fold 设计基于预测域 predpoints。
#     predpoints 使用 DEM + BIO 气候变量的环境分层抽样，而不是简单 regular/random。
# ============================================================

summarize_predictor_dataframe <- function(df, vars, dataset_name) {
  rows <- vector("list", length(vars))
  row_i <- 0L

  for (v in vars) {
    z <- safe_as_numeric(df[[v]])
    z <- z[is.finite(z)]

    if (length(z) == 0) {
      row <- data.frame(
        dataset = dataset_name,
        variable = v,
        n = 0,
        min = NA_real_,
        q05 = NA_real_,
        mean = NA_real_,
        median = NA_real_,
        q95 = NA_real_,
        max = NA_real_,
        sd = NA_real_
      )
    } else {
      row <- data.frame(
        dataset = dataset_name,
        variable = v,
        n = length(z),
        min = min(z, na.rm = TRUE),
        q05 = as.numeric(stats::quantile(z, 0.05, na.rm = TRUE, names = FALSE)),
        mean = mean(z, na.rm = TRUE),
        median = stats::median(z, na.rm = TRUE),
        q95 = as.numeric(stats::quantile(z, 0.95, na.rm = TRUE, names = FALSE)),
        max = max(z, na.rm = TRUE),
        sd = stats::sd(z, na.rm = TRUE)
      )
    }

    row_i <- row_i + 1L
    rows[[row_i]] <- row
  }

  return(do.call(rbind, rows[seq_len(row_i)]))
}

allocate_cluster_samples <- function(cluster_counts, target_total, min_per_cluster = 200) {
  cluster_counts <- as.integer(cluster_counts)
  target_total <- as.integer(min(target_total, sum(cluster_counts)))

  if (target_total <= 0 || sum(cluster_counts) <= 0) {
    return(rep(0L, length(cluster_counts)))
  }

  n <- rep(0L, length(cluster_counts))
  active <- which(cluster_counts > 0)

  # 如果样本数足够，先保证每个环境簇至少 1 个样本；否则优先覆盖最大簇
  if (target_total >= length(active)) {
    n[active] <- 1L
    remaining <- target_total - sum(n)
  } else {
    ord <- active[order(cluster_counts[active], decreasing = TRUE)]
    n[ord[seq_len(target_total)]] <- 1L
    remaining <- 0L
  }

  # 再尽量补到 min_per_cluster，但不超过各簇可用像元数
  if (remaining > 0 && min_per_cluster > 1) {
    ord <- active[order(cluster_counts[active], decreasing = FALSE)]
    for (idx in ord) {
      add <- min(
        remaining,
        cluster_counts[idx] - n[idx],
        min_per_cluster - n[idx]
      )
      if (add > 0) {
        n[idx] <- n[idx] + add
        remaining <- target_total - sum(n)
      }
      if (remaining <= 0) break
    }
  }

  # 剩余样本按各簇剩余容量近似比例分配
  while (remaining > 0) {
    capacity <- cluster_counts - n
    if (all(capacity <= 0)) break

    add <- floor(remaining * capacity / sum(capacity))
    add <- pmin(add, capacity)

    if (sum(add) == 0) {
      ord <- order(capacity, decreasing = TRUE)
      ord <- ord[capacity[ord] > 0]
      n_add <- min(remaining, length(ord))
      add[ord[seq_len(n_add)]] <- 1L
    }

    n <- n + add
    remaining <- target_total - sum(n)
  }

  return(as.integer(n))
}

make_environmental_stratified_predpoints <- function(env_pred,
                                                      sample_size = 100000,
                                                      n_clusters = 50,
                                                      min_per_cluster = 200,
                                                      pca_components = 4,
                                                      candidate_multiplier = 3,
                                                      kmeans_nstart = 5,
                                                      kmeans_iter_max = 100,
                                                      seed = 123,
                                                      out_dir = ".") {
  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("需要安装 sf 包：install.packages('sf')")
  }

  message("Preparing kNNDM predpoints using environmental stratified sampling...")
  message("  sample_size=", sample_size,
          ", n_clusters=", n_clusters,
          ", min_per_cluster=", min_per_cluster,
          ", pca_components=", pca_components)

  set.seed(seed)

  vals <- terra::values(env_pred, mat = TRUE)
  colnames(vals) <- names(env_pred)

  ok <- rowSums(!is.finite(vals)|vals == -9999)==0L
  ok_idx <- which(ok)

  if (length(ok_idx) == 0) {
    stop("ENV_pred 中没有完整有效像元，无法构建 kNNDM predpoints。")
  }

  # 完整预测域环境摘要，用于诊断抽样代表性
  full_df <- as.data.frame(vals[ok_idx, , drop = FALSE])
  names(full_df) <- names(env_pred)
  full_summary <- summarize_predictor_dataframe(
    full_df,
    vars = raster_predictor_vars,
    dataset_name = "full_ENV_valid_cells"
  )

  target_total <- min(as.integer(sample_size), length(ok_idx))
  candidate_size <- min(
    length(ok_idx),
    max(target_total, as.integer(target_total * candidate_multiplier), n_clusters * 10)
  )

  if (length(ok_idx) > candidate_size) {
    cand_cells <- safe_sample(ok_idx, candidate_size)
  } else {
    cand_cells <- ok_idx
  }

  cand_vals <- vals[cand_cells, , drop = FALSE]
  cand_df <- as.data.frame(cand_vals)
  names(cand_df) <- names(env_pred)

  # 去掉零方差变量，避免 PCA 出错；这些变量仍保留在输出诊断中
  sds <- apply(cand_vals, 2, stats::sd, na.rm = TRUE)
  keep_var <- is.finite(sds) & sds > 0

  if (sum(keep_var) == 0) {
    stop("kNNDM environmental stratified sampling 失败：候选预测变量全部为零方差。")
  }

  cand_scaled <- scale(cand_vals[, keep_var, drop = FALSE])
  cand_scaled[!is.finite(cand_scaled)] <- 0

  npc <- min(pca_components, ncol(cand_scaled), nrow(cand_scaled) - 1)

  if (npc >= 1 && ncol(cand_scaled) > 1) {
    pca <- stats::prcomp(cand_scaled, center = FALSE, scale. = FALSE)
    pca_x <- pca$x[, seq_len(npc), drop = FALSE]
  } else {
    pca_x <- cand_scaled[, seq_len(min(1, ncol(cand_scaled))), drop = FALSE]
  }

  # kmeans 要求聚类中心数不能超过唯一环境组合数
  distinct_env_n <- nrow(unique(as.data.frame(pca_x)))
  n_clusters_use <- min(as.integer(n_clusters), nrow(pca_x), distinct_env_n)

  cluster_counts <- NULL
  n_by_cluster <- NULL

  if (n_clusters_use < 2) {
    warning("kNNDM predpoints 候选像元过少或环境空间不足，退化为简单随机抽样。")
    sampled_pos <- safe_sample(seq_len(nrow(cand_vals)), size = target_total)
    cluster_id <- rep(1L, nrow(cand_vals))
    n_by_cluster <- target_total
  } else {
    km <- stats::kmeans(
      pca_x,
      centers = n_clusters_use,
      iter.max = kmeans_iter_max,
      nstart = kmeans_nstart
    )

    cluster_id <- km$cluster
    cluster_counts <- tabulate(cluster_id, nbins = n_clusters_use)

    n_by_cluster <- allocate_cluster_samples(
      cluster_counts = cluster_counts,
      target_total = target_total,
      min_per_cluster = min_per_cluster
    )

    sampled_pos <- integer(0)

    for (cl in seq_len(n_clusters_use)) {
      idx <- which(cluster_id == cl)
      n_i <- min(length(idx), n_by_cluster[cl])
      if (n_i > 0) {
        sampled_pos <- c(sampled_pos, safe_sample(idx, size = n_i))
      }
    }

    # 理论上不会不足；如果因为极端小簇/四舍五入不足，从候选池补齐
    if (length(sampled_pos) < target_total) {
      rest <- setdiff(seq_len(nrow(cand_vals)), sampled_pos)
      n_add <- min(length(rest), target_total - length(sampled_pos))
      if (n_add > 0) {
        sampled_pos <- c(sampled_pos, safe_sample(rest, n_add))
      }
    }
  }

  sampled_pos <- sampled_pos[seq_len(min(length(sampled_pos), target_total))]
  sampled_cells <- cand_cells[sampled_pos]
  sampled_cluster <- cluster_id[sampled_pos]

  xy <- terra::xyFromCell(env_pred[[1]], sampled_cells)
  sample_values <- as.data.frame(cand_vals[sampled_pos, , drop = FALSE])
  names(sample_values) <- names(env_pred)

  sample_df <- data.frame(
    cell = sampled_cells,
    x = xy[, 1],
    y = xy[, 2],
    env_cluster = sampled_cluster
  )
  sample_df <- cbind(sample_df, sample_values)

  sample_summary <- summarize_predictor_dataframe(
    sample_values,
    vars = raster_predictor_vars,
    dataset_name = "knndm_predpoints_environmental_stratified_sample"
  )

  cluster_summary <- aggregate(
    cell ~ env_cluster,
    data = sample_df,
    FUN = length
  )
  names(cluster_summary)[names(cluster_summary) == "cell"] <- "sampled_cells"

  if (!is.null(cluster_counts)) {
    candidate_cluster_summary <- data.frame(
      env_cluster = seq_along(cluster_counts),
      candidate_cells = as.integer(cluster_counts),
      target_sample_cells = as.integer(n_by_cluster)
    )
    cluster_summary <- merge(
      candidate_cluster_summary,
      cluster_summary,
      by = "env_cluster",
      all.x = TRUE
    )
    cluster_summary$sampled_cells[is.na(cluster_summary$sampled_cells)] <- 0
  }

  # 写出 predpoints 抽样诊断
  write.csv(
    sample_df,
    file.path(out_dir, "knndm_predpoints_sample_points.csv"),
    row.names = FALSE
  )

  write.csv(
    rbind(full_summary, sample_summary),
    file.path(out_dir, "knndm_predpoints_env_summary_full_vs_sample.csv"),
    row.names = FALSE
  )

  write.csv(
    cluster_summary,
    file.path(out_dir, "knndm_predpoints_cluster_summary.csv"),
    row.names = FALSE
  )

  predpoints_sf <- sf::st_as_sf(
    sample_df[, c("cell", "env_cluster", "x", "y")],
    coords = c("x", "y"),
    crs = sf::st_crs(terra::crs(env_pred))
  )

  if (write_knndm_predpoints_gpkg) {
    gpkg_file <- file.path(out_dir, "knndm_predpoints_sample_points.gpkg")
    if (file.exists(gpkg_file)) file.remove(gpkg_file)
    sf::st_write(predpoints_sf, gpkg_file, quiet = TRUE)
  }

  rm(vals, full_df, cand_vals, cand_df)
  gc()

  return(list(
    predpoints = predpoints_sf,
    sample_df = sample_df,
    sample_values = sample_values,
    full_summary = full_summary,
    sample_summary = sample_summary,
    cluster_summary = cluster_summary
  ))
}

write_training_vs_predpoints_summary <- function(train_dat,
                                                  predpoints_values,
                                                  response_var,
                                                  out_dir) {
  train_summary <- summarize_predictor_dataframe(
    train_dat,
    vars = raster_predictor_vars,
    dataset_name = paste0("training_points_", response_var)
  )

  predpoints_summary <- summarize_predictor_dataframe(
    predpoints_values,
    vars = raster_predictor_vars,
    dataset_name = "knndm_predpoints_sample"
  )

  out <- rbind(train_summary, predpoints_summary)

  write.csv(
    out,
    file.path(out_dir, paste0("knndm_training_vs_predpoints_env_summary_", response_var, ".csv")),
    row.names = FALSE
  )

  return(out)
}

make_knndm_folds <- function(train_dat,
                             env_pred,
                             predpoints_sample,
                             k = 5,
                             clustering = "hierarchical",
                             maxp = 0.5,
                             space = "geographical",
                             seed = 123) {
  if (!requireNamespace("CAST", quietly = TRUE)) {
    stop("需要安装 CAST 包：install.packages('CAST')")
  }
  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("需要安装 sf 包：install.packages('sf')")
  }

  set.seed(seed)

  tpoints <- sf::st_as_sf(
    train_dat,
    coords = c("x_env", "y_env"),
    crs = sf::st_crs(terra::crs(env_pred))
  )

  run_knndm_once <- function(cl) {
    base_args <- list(
      tpoints = tpoints,
      predpoints = predpoints_sample,
      k = k,
      maxp = maxp,
      clustering = cl
    )
    # CAST 不同版本中地理/特征空间的参数名不同：部分版本用 dist_space，部分版本用 space。
    # 与 tune_qrf_params 脚本保持一致：先试 dist_space，失败再回退 space。
    out <- try(do.call(CAST::knndm, c(base_args, list(dist_space = space))), silent = TRUE)
    if (!inherits(out, "try-error")) return(out)
    do.call(CAST::knndm, c(base_args, list(space = space)))
  }

  # kmeans 在训练点存在重复坐标（例如同一样地多年重复测量，long/lat 完全相同）时，
  # 会报 "more cluster centers than distinct data points"。此处自动改用 hierarchical
  # 重试；hierarchical 对重复点稳健（仅在训练点数非常大时更耗内存）。
  folds <- tryCatch(
    run_knndm_once(clustering),
    error = function(e) {
      msg <- conditionMessage(e)
      if (!identical(clustering, "hierarchical")) {
        message("CAST::knndm 用 clustering='", clustering, "' 失败（", msg,
                "）。自动改用 clustering='hierarchical' 重试 ...")
        tryCatch(
          run_knndm_once("hierarchical"),
          error = function(e2) {
            stop("CAST::knndm 在 hierarchical 下仍失败。可尝试降低 knndm_k 或增大 knndm_maxp。",
                 "原始错误：", conditionMessage(e2))
          }
        )
      } else {
        stop("CAST::knndm 生成 folds 失败。可尝试降低 knndm_k 或增大 knndm_maxp。",
             "原始错误：", msg)
      }
    }
  )

  if (is.null(folds$indx_train) || is.null(folds$indx_test)) {
    stop("CAST::knndm 返回对象中缺少 indx_train / indx_test，请检查 CAST 版本。")
  }

  return(folds)
}

validate_fold_indices <- function(folds, n) {
  if (length(folds$indx_train) != length(folds$indx_test)) {
    stop("kNNDM folds 的 indx_train 与 indx_test 长度不一致。")
  }

  for (i in seq_along(folds$indx_train)) {
    tr <- folds$indx_train[[i]]
    te <- folds$indx_test[[i]]

    if (length(tr) == 0 || length(te) == 0) {
      stop("kNNDM 第 ", i, " 折训练集或测试集为空。")
    }
    if (any(tr < 1 | tr > n) || any(te < 1 | te > n)) {
      stop("kNNDM 第 ", i, " 折索引超出训练数据行数。")
    }
    if (length(intersect(tr, te)) > 0) {
      stop("kNNDM 第 ", i, " 折训练/测试索引有重叠。")
    }
  }

  all_test <- unlist(folds$indx_test, use.names = FALSE)
  if (length(all_test) != length(unique(all_test))) {
    stop("kNNDM 测试索引在不同 folds 间重复，无法进行严格 k-fold OOS 验证。")
  }
  if (!setequal(all_test, seq_len(n))) {
    stop("kNNDM 测试索引没有覆盖全部训练样本；请检查 CAST::knndm 输出或 fold 参数。")
  }

  invisible(TRUE)
}

knndm_metrics <- data.frame()
knndm_metric_rows <- list()
knndm_info_rows <- list()
knndm_folds_by_response <- list()
knndm_predpoints_obj <- NULL

if (run_knndm_validation) {
  message("==================================================")
  message("Running kNNDM CV spatial validation...")
  message("Note: default clustering is kmeans for projected CRS efficiency; switch to hierarchical if kmeans fails.")
  message("kNNDM predpoints sampling: ", knndm_predpoints_sampling)
  message("kNNDM parameters: k=", knndm_k,
          ", samplesize=", knndm_samplesize,
          ", env_clusters=", knndm_env_clusters,
          ", min_per_cluster=", knndm_min_per_cluster,
          ", pca_components=", knndm_pca_components,
          ", clustering=", knndm_clustering,
          ", maxp=", knndm_maxp,
          ", space=", knndm_space)

  if (knndm_predpoints_sampling != "environmental_stratified") {
    stop("当前 v8 脚本只实现 knndm_predpoints_sampling = 'environmental_stratified'。")
  }

  # predpoints 由预测域 ENV 有效像元做环境分层抽样得到；所有响应变量共用同一个预测域样本。
  knndm_predpoints_obj <- make_environmental_stratified_predpoints(
    env_pred = ENV_pred,
    sample_size = knndm_samplesize,
    n_clusters = knndm_env_clusters,
    min_per_cluster = knndm_min_per_cluster,
    pca_components = knndm_pca_components,
    candidate_multiplier = knndm_candidate_multiplier,
    kmeans_nstart = knndm_kmeans_nstart,
    kmeans_iter_max = knndm_kmeans_iter_max,
    seed = seed,
    out_dir = out_dir
  )

  saveRDS(
    knndm_predpoints_obj,
    file.path(out_dir, "knndm_predpoints_environmental_stratified_object.rds")
  )

  for (resp in response_vars) {
    message("--------------------------------------------------")
    message("kNNDM CV for response: ", resp)

    d <- train_list[[resp]]
    n <- nrow(d)

    write_training_vs_predpoints_summary(
      train_dat = d,
      predpoints_values = knndm_predpoints_obj$sample_values,
      response_var = resp,
      out_dir = out_dir
    )

    if (n < max(min_train_rows, knndm_k * 2)) {
      warning(resp, ": 样本量过少，跳过 kNNDM CV。n = ", n)
      next
    }

    folds <- make_knndm_folds(
      train_dat = d,
      env_pred = ENV_pred,
      predpoints_sample = knndm_predpoints_obj$predpoints,
      k = knndm_k,
      clustering = knndm_clustering,
      maxp = knndm_maxp,
      space = knndm_space,
      seed = seed
    )

    validate_fold_indices(folds, n)
    knndm_folds_by_response[[resp]] <- folds

    # 保存 fold 索引，便于复现和诊断
    saveRDS(
      folds,
      file.path(out_dir, paste0("knndm_folds_", resp, ".rds"))
    )

    if (!is.null(folds$W)) {
      knndm_info_rows[[resp]] <- data.frame(
        response = resp,
        k = knndm_k,
        predpoints_sampling = knndm_predpoints_sampling,
        samplesize = nrow(knndm_predpoints_obj$sample_df),
        env_clusters = knndm_env_clusters,
        min_per_cluster = knndm_min_per_cluster,
        pca_components = knndm_pca_components,
        clustering = knndm_clustering,
        maxp = knndm_maxp,
        space = knndm_space,
        W = as.numeric(folds$W)
      )
    }

    fold_pred_rows <- list()

    for (fold_i in seq_along(folds$indx_test)) {
      message("  Fold ", fold_i, " / ", length(folds$indx_test))

      train_idx <- folds$indx_train[[fold_i]]
      test_idx <- folds$indx_test[[fold_i]]

      d_train <- d[train_idx, , drop = FALSE]
      d_test <- d[test_idx, , drop = FALSE]

      if (model_type == "RANGER_QRF") {
        m_cv <- train_ranger_qrf_model(d_train, response_var = resp, importance = "none", num_trees = get_qrf_param_value(resp, "qrf_num_trees_validation"), model_seed = seed + fold_i)
        p_test <- predict_ranger_qrf_dataframe(m_cv, d_test, quantiles)
      } else if (model_type == "XGBOOST") {
        m_cv <- train_xgb_quantile_models(d_train, quantiles, nrounds = max(100, floor(xgb_nrounds * 0.6)), model_seed = seed + fold_i)
        p_test <- predict_xgb_dataframe(m_cv, d_test)
      } else {
        stop("model_type 只能是 'RANGER_QRF' 或 'XGBOOST'。")
      }

      p_test <- postprocess_prediction_df(p_test, resp)

      pred_df <- data.frame(
        response = resp,
        model_type = model_type,
        validation_method = "knndm",
        fold_id = fold_i,
        row_id = test_idx,
        y = d_test$response
      )

      pred_df <- cbind(pred_df, p_test)
      fold_pred_rows[[fold_i]] <- pred_df

      # fold-level 指标：仅用于诊断，主报告建议看全局汇总
      for (qi in seq_along(quantiles)) {
        q <- quantiles[qi]
        qn <- q_names[qi]
        y <- d_test$response
        pred <- p_test[[qn]]

        knndm_metric_rows[[length(knndm_metric_rows) + 1L]] <- data.frame(
          response = resp,
          model_type = model_type,
          validation_method = "knndm_fold",
          fold_id = fold_i,
          quantile = q,
          n_train = nrow(d_train),
          n_validation = nrow(d_test),
          pinball_loss = pinball_loss(y, pred, q),
          empirical_coverage = mean(y <= pred, na.rm = TRUE),
          prediction_min = min(pred, na.rm = TRUE),
          prediction_max = max(pred, na.rm = TRUE)
        )
      }
    }

    pred_all <- do.call(rbind, fold_pred_rows)

    write.csv(
      pred_all,
      file.path(out_dir, paste0("validation_predictions_knndm_", resp, ".csv")),
      row.names = FALSE
    )

    # global validation：kNNDM folds 可能不平衡，推荐把所有 OOS 预测合并后统一计算
    for (qi in seq_along(quantiles)) {
      q <- quantiles[qi]
      qn <- q_names[qi]
      y <- pred_all$y
      pred <- pred_all[[qn]]

      knndm_metric_rows[[length(knndm_metric_rows) + 1L]] <- data.frame(
        response = resp,
        model_type = model_type,
        validation_method = "knndm_global",
        fold_id = NA_integer_,
        quantile = q,
        n_train = NA_integer_,
        n_validation = nrow(pred_all),
        pinball_loss = pinball_loss(y, pred, q),
        empirical_coverage = mean(y <= pred, na.rm = TRUE),
        prediction_min = min(pred, na.rm = TRUE),
        prediction_max = max(pred, na.rm = TRUE)
      )
    }
  }

  knndm_metrics <- if (length(knndm_metric_rows) > 0) do.call(rbind, knndm_metric_rows) else data.frame()

  if (nrow(knndm_metrics) > 0) {
    write.csv(
      knndm_metrics,
      file.path(out_dir, "validation_metrics_knndm.csv"),
      row.names = FALSE
    )

    knndm_global <- knndm_metrics[knndm_metrics$validation_method == "knndm_global", , drop = FALSE]
    write.csv(
      knndm_global,
      file.path(out_dir, "validation_metrics_knndm_global_summary.csv"),
      row.names = FALSE
    )

    print(knndm_global)
  }

  if (length(knndm_info_rows) > 0) {
    knndm_info <- do.call(rbind, knndm_info_rows)
    write.csv(
      knndm_info,
      file.path(out_dir, "validation_knndm_fold_assignment_info.csv"),
      row.names = FALSE
    )
  }
}

# ============================================================
# 11b. 输出统一验证报告：random hold-out + kNNDM
# ============================================================

summarize_validation_metrics <- function(df, method_label) {
  if (is.null(df) || nrow(df) == 0) return(data.frame())
  rows <- list()
  ii <- 0L
  for (resp in unique(df$response)) {
    for (q in sort(unique(df$quantile))) {
      z <- df[df$response == resp & df$quantile == q, , drop = FALSE]
      if (nrow(z) == 0) next
      ii <- ii + 1L
      rows[[ii]] <- data.frame(
        validation_method = method_label,
        response = resp,
        model_type = unique(z$model_type)[1],
        quantile = q,
        target_coverage = q,
        empirical_coverage_mean = mean(z$empirical_coverage, na.rm = TRUE),
        empirical_coverage_sd = if (nrow(z) > 1) stats::sd(z$empirical_coverage, na.rm = TRUE) else NA_real_,
        coverage_error_mean = mean(z$empirical_coverage - q, na.rm = TRUE),
        abs_coverage_error_mean = mean(abs(z$empirical_coverage - q), na.rm = TRUE),
        undercoverage_mean = mean(pmax(0, q - z$empirical_coverage), na.rm = TRUE),
        pinball_loss_mean = mean(z$pinball_loss, na.rm = TRUE),
        pinball_loss_sd = if (nrow(z) > 1) stats::sd(z$pinball_loss, na.rm = TRUE) else NA_real_,
        n_validation_total = sum(z$n_validation, na.rm = TRUE),
        n_rows = nrow(z),
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

validation_report_rows <- list()
if (exists("validation_metrics") && nrow(validation_metrics) > 0) {
  validation_report_rows[["repeated_holdout"]] <- summarize_validation_metrics(validation_metrics, "repeated_holdout")
}
if (exists("knndm_metrics") && nrow(knndm_metrics) > 0) {
  knndm_global_for_report <- knndm_metrics[knndm_metrics$validation_method == "knndm_global", , drop = FALSE]
  validation_report_rows[["knndm_global"]] <- summarize_validation_metrics(knndm_global_for_report, "knndm_global")
}
validation_report <- if (length(validation_report_rows) > 0) do.call(rbind, validation_report_rows) else data.frame()
if (nrow(validation_report) > 0) {
  write.csv(validation_report, file.path(out_dir, "validation_report_holdout_and_knndm.csv"), row.names = FALSE)
  print(validation_report)
}

# ============================================================
# 11c. 可选：树数稳定性检查 800/1200/1500/2000
#      只比较每个响应变量已选最优参数下的树数，不重新搜索 mtry/node。
# ============================================================

if (isTRUE(run_tree_count_sensitivity) && model_type == "RANGER_QRF") {
  message("==================================================")
  message("Running optional tree-count sensitivity check: ", paste(qrf_num_trees_final_candidates, collapse = ", "))
  tree_metrics <- data.frame()
  tree_metric_rows <- list()

  for (resp in response_vars) {
    d <- train_list[[resp]]
    n <- nrow(d)
    if (n < 50) next
    val_n <- max(1, floor(n * validation_fraction))

    for (ntree in qrf_num_trees_final_candidates) {
      for (rep_i in seq_len(validation_repeats)) {
        set.seed(seed + 7000 + rep_i)
        val_idx <- safe_sample(seq_len(n), size = val_n)
        d_train <- d[-val_idx, , drop = FALSE]
        d_val <- d[val_idx, , drop = FALSE]

        m_tree <- train_ranger_qrf_model(
          d_train,
          response_var = resp,
          importance = "none",
          num_trees = ntree,
          model_seed = seed + 9000 + ntree + rep_i
        )
        p_tree <- predict_ranger_qrf_dataframe(m_tree, d_val, quantiles)
        p_tree <- postprocess_prediction_df(p_tree, resp)

        for (qi in seq_along(quantiles)) {
          q <- quantiles[qi]
          qn <- q_names[qi]
          y <- d_val$response
          pred <- p_tree[[qn]]
          tree_metric_rows[[length(tree_metric_rows) + 1L]] <- data.frame(
            response = resp,
            model_type = model_type,
            validation_method = "tree_count_repeated_holdout",
            num_trees = ntree,
            repeat_id = rep_i,
            quantile = q,
            n_train = nrow(d_train),
            n_validation = nrow(d_val),
            pinball_loss = pinball_loss(y, pred, q),
            empirical_coverage = mean(y <= pred, na.rm = TRUE),
            stringsAsFactors = FALSE
          )
        }
      }
      gc()
    }
  }

  tree_metrics <- if (length(tree_metric_rows) > 0) do.call(rbind, tree_metric_rows) else data.frame()
  write.csv(tree_metrics, file.path(out_dir, "qrf_tree_count_sensitivity_metrics.csv"), row.names = FALSE)

  if (nrow(tree_metrics) > 0) {
    tree_summary <- aggregate(
      cbind(pinball_loss, empirical_coverage) ~ response + model_type + validation_method + num_trees + quantile,
      data = tree_metrics,
      FUN = function(z) c(mean = mean(z, na.rm = TRUE), sd = stats::sd(z, na.rm = TRUE))
    )
    tree_summary <- do.call(data.frame, tree_summary)
    names(tree_summary) <- gsub("\\.", "_", names(tree_summary))
    tree_summary$coverage_error_mean <- tree_summary$empirical_coverage_mean - tree_summary$quantile
    tree_summary$abs_coverage_error_mean <- abs(tree_summary$coverage_error_mean)
    write.csv(tree_summary, file.path(out_dir, "qrf_tree_count_sensitivity_summary.csv"), row.names = FALSE)
    print(tree_summary)
  }
}

# ============================================================
# 12. 使用全部训练数据训练最终模型
# ============================================================

models <- list()

for (resp in response_vars) {
  message("==================================================")
  message("Training final model for response: ", resp)

  if (model_type == "RANGER_QRF") {
    models[[resp]] <- train_ranger_qrf_model(train_list[[resp]], response_var = resp, importance = "permutation", num_trees = get_qrf_param_value(resp, "qrf_num_trees_final"), model_seed = seed)

    vi <- models[[resp]]$variable.importance
    vi_df <- data.frame(
      variable = names(vi),
      importance = as.numeric(vi)
    )
    vi_df <- vi_df[order(vi_df$importance, decreasing = TRUE), ]

    write.csv(
      vi_df,
      file.path(out_dir, paste0(model_type, "_", resp, "_variable_importance.csv")),
      row.names = FALSE
    )

  } else if (model_type == "XGBOOST") {
    models[[resp]] <- train_xgb_quantile_models(train_list[[resp]], quantiles, nrounds = xgb_nrounds, model_seed = seed)

    for (qn in names(models[[resp]])) {
      imp <- xgboost::xgb.importance(
        feature_names = predictor_vars,
        model = models[[resp]][[qn]]
      )
      write.csv(
        imp,
        file.path(out_dir, paste0(model_type, "_", resp, "_", qn, "_variable_importance.csv")),
        row.names = FALSE
      )
    }
  } else {
    stop("model_type 只能是 'RANGER_QRF' 或 'XGBOOST'。")
  }
}

# ============================================================
# 12b. CAST Area of Applicability (AOA) 多变量适用域
#      每个响应变量、每个固定 Forest_age 单独计算 AOA/DI。
# ============================================================

# 供 AOA、快速矩阵预测和 terra 分块预测共同使用。
make_age_raster <- function(template_raster, age_value) {
  age_r <- terra::init(template_raster, age_value)
  names(age_r) <- "Forest_age"
  return(age_r)
}

make_predictor_stack <- function(env_pred, age_value) {
  age_r <- make_age_raster(env_pred[[1]], age_value)
  pred_stack <- c(env_pred, age_r)
  names(pred_stack) <- predictor_vars
  return(pred_stack)
}

make_aoa_weight_for_variables <- function(model = NULL, variables, use_variable_importance = TRUE) {
  if (!isTRUE(use_variable_importance) || is.null(model) || is.null(model$variable.importance)) {
    w <- rep(1, length(variables))
    names(w) <- variables
  } else {
    vi <- model$variable.importance
    w <- vi[variables]
    w[is.na(w) | !is.finite(w) | w < 0] <- 0
    if (sum(w, na.rm = TRUE) <= 0) {
      w <- rep(1, length(variables))
      names(w) <- variables
    }
  }

  w <- w / mean(w, na.rm = TRUE)
  w <- as.numeric(w)
  names(w) <- variables

  # CAST::aoa/trainDI expects user weights as a one-row data.frame with
  # predictor variables as column names, not a named numeric vector.
  weight_df <- as.data.frame(as.list(w), check.names = FALSE)
  names(weight_df) <- variables
  return(weight_df)
}

make_aoa_weight <- function(model, response_var) {
  make_aoa_weight_for_variables(
    model = model,
    variables = predictor_vars,
    use_variable_importance = aoa_use_variable_importance
  )
}

default_aoa_weight <- function(variables = predictor_vars) {
  make_aoa_weight_for_variables(
    model = NULL,
    variables = variables,
    use_variable_importance = FALSE
  )
}

write_aoa_summary <- function(aoa_r, di_r, response_var, age_value, filename,
                              predictor_set = "full_predictor_space") {
  # AOA is a 0/1 raster: inside = sum of 1-valued cells; valid = non-NA cells.
  inside <- as.numeric(terra::global(aoa_r, "sum", na.rm = TRUE)[1, 1])
  valid <- as.numeric(terra::global(!is.na(aoa_r), "sum", na.rm = TRUE)[1, 1])
  outside <- valid - inside
  di_mean <- as.numeric(terra::global(di_r, "mean", na.rm = TRUE)[1, 1])
  di_max <- as.numeric(terra::global(di_r, "max", na.rm = TRUE)[1, 1])

  out <- data.frame(
    predictor_set = predictor_set,
    response = response_var,
    Forest_age = age_value,
    valid_cells = valid,
    inside_AOA_cells = inside,
    outside_AOA_cells = outside,
    inside_AOA_percent = ifelse(valid > 0, 100 * inside / valid, NA_real_),
    outside_AOA_percent = ifelse(valid > 0, 100 * outside / valid, NA_real_),
    DI_mean = di_mean,
    DI_max = di_max,
    stringsAsFactors = FALSE
  )
  write.csv(out, filename, row.names = FALSE)
  out
}

aoa_results <- list()
aoa_summary_rows <- list()

if (isTRUE(make_aoa)) {
  message("==================================================")
  message("Calculating CAST Area of Applicability (AOA)...")

  for (resp in response_vars) {
    d <- train_list[[resp]]
    train_aoa <- d[, predictor_vars, drop = FALSE]
    weight_aoa <- if (model_type == "RANGER_QRF") {
      make_aoa_weight(models[[resp]], resp)
    } else {
      default_aoa_weight()
    }
    cvtest <- NULL
    cvtrain <- NULL
    if (isTRUE(aoa_use_knndm_folds) && !is.null(knndm_folds_by_response[[resp]])) {
      # CAST::aoa() 与 CAST::CreateSpacetimeFolds/knndm 输出兼容：
      # CVtest/CVtrain 应传入“每折索引列表”，不要传入折号向量。
      cvtest <- knndm_folds_by_response[[resp]]$indx_test
      cvtrain <- knndm_folds_by_response[[resp]]$indx_train
    } else if (isTRUE(aoa_use_knndm_folds)) {
      message("  AOA for ", resp, ": no kNNDM folds available; using CAST::aoa without CVtest/CVtrain.")
    }

    aoa_results[[resp]] <- list()

    ref_aoa_obj <- NULL
    aoa_loop_values <- if (isTRUE(aoa_compute_per_age)) age_values else unique(c(ref_age, age_values))

    for (age in aoa_loop_values) {
      if (!isTRUE(aoa_compute_per_age) && !identical(as.numeric(age), as.numeric(ref_age))) {
        message("  AOA for ", resp, ", Forest_age=", age, " reused from ref_age=", ref_age)
        if (is.null(ref_aoa_obj)) {
          ref_aoa_file <- file.path(out_dir, paste0("AOA_", resp, "_ForestAge_", ref_age, "_AOA.tif"))
          ref_di_file <- file.path(out_dir, paste0("AOA_", resp, "_ForestAge_", ref_age, "_DI.tif"))
          if (file.exists(ref_aoa_file) && file.exists(ref_di_file)) {
            ref_aoa_obj <- list(AOA = rast(ref_aoa_file), DI = rast(ref_di_file))
          } else {
            stop("aoa_compute_per_age=FALSE 需要先生成 ref_age 的 AOA/DI: ", ref_age)
          }
        }
        # 写出同名 AOA/DI 文件，保持不同 Forest_age 输出文件齐全；内容复用 ref_age。
        aoa_file_reuse <- file.path(out_dir, paste0("AOA_", resp, "_ForestAge_", age, "_AOA.tif"))
        di_file_reuse <- file.path(out_dir, paste0("AOA_", resp, "_ForestAge_", age, "_DI.tif"))
        if (!file.exists(aoa_file_reuse)) write_raster_safe(ref_aoa_obj$AOA, aoa_file_reuse)
        if (!file.exists(di_file_reuse)) write_raster_safe(ref_aoa_obj$DI, di_file_reuse)
        aoa_results[[resp]][[as.character(age)]] <- list(AOA = rast(aoa_file_reuse), DI = rast(di_file_reuse))
        aoa_summary_rows[[paste(resp, age, sep = "_")]] <- write_aoa_summary(
          aoa_r = ref_aoa_obj$AOA, di_r = ref_aoa_obj$DI, response_var = resp, age_value = age,
          filename = file.path(out_dir, paste0("AOA_", resp, "_ForestAge_", age, "_summary.csv"))
        )
        next
      }
      message("  AOA for ", resp, ", Forest_age=", age)
      aoa_file <- file.path(out_dir, paste0("AOA_", resp, "_ForestAge_", age, "_AOA.tif"))
      di_file <- file.path(out_dir, paste0("AOA_", resp, "_ForestAge_", age, "_DI.tif"))
      if (isTRUE(aoa_skip_existing) && file.exists(aoa_file) && file.exists(di_file)) {
        message("    Existing AOA/DI found; loading instead of recalculating.")
        aoa_obj <- list(AOA = rast(aoa_file), DI = rast(di_file))
        if (!isTRUE(aoa_compute_per_age) && identical(as.numeric(age), as.numeric(ref_age))) {
          ref_aoa_obj <- aoa_obj
        }
        aoa_results[[resp]][[as.character(age)]] <- aoa_obj
        aoa_summary_rows[[paste(resp, age, sep = "_")]] <- write_aoa_summary(
          aoa_r = aoa_obj$AOA, di_r = aoa_obj$DI, response_var = resp, age_value = age,
          filename = file.path(out_dir, paste0("AOA_", resp, "_ForestAge_", age, "_summary.csv"))
        )
        next
      }
      pred_stack <- make_predictor_stack(ENV_pred, age)

      aoa_obj <- CAST::aoa(
        newdata = pred_stack,
        train = train_aoa,
        variables = predictor_vars,
        weight = weight_aoa,
        CVtest = cvtest,
        CVtrain = cvtrain,
        method = aoa_method,
        LPD = aoa_lpd,
        verbose = TRUE
      )

      if (!isTRUE(aoa_compute_per_age) && identical(as.numeric(age), as.numeric(ref_age))) {
        ref_aoa_obj <- aoa_obj
      }
      aoa_results[[resp]][[as.character(age)]] <- aoa_obj

      write_raster_safe(aoa_obj$AOA, aoa_file)
      write_raster_safe(aoa_obj$DI, di_file)

      aoa_summary_rows[[paste(resp, age, sep = "_")]] <- write_aoa_summary(
        aoa_r = aoa_obj$AOA,
        di_r = aoa_obj$DI,
        response_var = resp,
        age_value = age,
        filename = file.path(out_dir, paste0("AOA_", resp, "_ForestAge_", age, "_summary.csv"))
      )
    }
  }

  if (length(aoa_summary_rows) > 0) {
    aoa_summary <- do.call(rbind, aoa_summary_rows)
    write.csv(aoa_summary, file.path(out_dir, "AOA_summary_all_responses_ages.csv"), row.names = FALSE)
    print(aoa_summary)
  }

}

# ============================================================
# 12c. 环境变量-only AOA：只看 DEM + BIO，不包含 Forest_age
#      该结果不随 80/100/120 林龄变化，可用于报告基础环境适用域。
# ============================================================

environment_only_aoa_result <- NULL

if (isTRUE(make_environment_only_aoa)) {
  message("==================================================")
  message("Calculating environment-only CAST Area of Applicability (AOA)...")
  message("Environment-only AOA variables: ", paste(raster_predictor_vars, collapse = ", "))

  env_aoa_response <- environment_only_aoa_reference_response
  if (is.null(env_aoa_response) || !(env_aoa_response %in% names(train_list))) {
    env_aoa_response <- response_vars[1]
    message("  environment_only_aoa_reference_response not found; using ", env_aoa_response)
  }

  d_env_aoa <- train_list[[env_aoa_response]]
  train_aoa_env <- d_env_aoa[, raster_predictor_vars, drop = FALSE]

  if (isTRUE(environment_only_aoa_use_variable_importance) && model_type == "RANGER_QRF") {
    weight_env_aoa <- make_aoa_weight_for_variables(
      model = models[[env_aoa_response]],
      variables = raster_predictor_vars,
      use_variable_importance = TRUE
    )
  } else {
    weight_env_aoa <- default_aoa_weight(raster_predictor_vars)
  }

  cvtest_env <- NULL
  cvtrain_env <- NULL
  if (isTRUE(aoa_use_knndm_folds) && !is.null(knndm_folds_by_response[[env_aoa_response]])) {
    cvtest_env <- knndm_folds_by_response[[env_aoa_response]]$indx_test
    cvtrain_env <- knndm_folds_by_response[[env_aoa_response]]$indx_train
  } else if (isTRUE(aoa_use_knndm_folds)) {
    message("  Environment-only AOA: no kNNDM folds available; using CAST::aoa without CVtest/CVtrain.")
  }

  env_aoa_file <- file.path(out_dir, "AOA_environment_only_AOA.tif")
  env_di_file <- file.path(out_dir, "AOA_environment_only_DI.tif")

  if (isTRUE(aoa_skip_existing) && file.exists(env_aoa_file) && file.exists(env_di_file)) {
    message("  Existing environment-only AOA/DI found; loading instead of recalculating.")
    environment_only_aoa_result <- list(AOA = rast(env_aoa_file), DI = rast(env_di_file))
  } else {
    environment_only_aoa_result <- CAST::aoa(
      newdata = ENV_pred,
      train = train_aoa_env,
      variables = raster_predictor_vars,
      weight = weight_env_aoa,
      CVtest = cvtest_env,
      CVtrain = cvtrain_env,
      method = aoa_method,
      LPD = aoa_lpd,
      verbose = TRUE
    )
    write_raster_safe(environment_only_aoa_result$AOA, env_aoa_file)
    write_raster_safe(environment_only_aoa_result$DI, env_di_file)
  }

  env_aoa_summary <- write_aoa_summary(
    aoa_r = environment_only_aoa_result$AOA,
    di_r = environment_only_aoa_result$DI,
    response_var = env_aoa_response,
    age_value = NA_real_,
    filename = file.path(out_dir, "AOA_environment_only_summary.csv"),
    predictor_set = "environment_only"
  )
  env_aoa_summary$variables_used <- paste(raster_predictor_vars, collapse = ",")
  env_aoa_summary$weight_mode <- ifelse(isTRUE(environment_only_aoa_use_variable_importance),
                                        "reference_response_variable_importance", "equal_weight")
  write.csv(env_aoa_summary, file.path(out_dir, "AOA_environment_only_summary.csv"), row.names = FALSE)
  print(env_aoa_summary)
}

# ============================================================
# 13. Forest_age PDP sanity check
#     用最终模型检查 80-120 年龄敏感性是否被模型学到
# ============================================================

make_forest_age_pdp <- function(model,
                                response_var,
                                background_dat,
                                age_grid,
                                sample_size = 2000,
                                seed = 123) {
  set.seed(seed)

  bg <- background_dat[, raster_predictor_vars, drop = FALSE]
  bg <- bg[complete.cases(bg), , drop = FALSE]

  if (nrow(bg) == 0) {
    warning(response_var, ": PDP 背景样本为空，跳过。")
    return(data.frame())
  }

  if (nrow(bg) > sample_size) {
    bg <- bg[safe_sample(seq_len(nrow(bg)), sample_size), , drop = FALSE]
  }

  rows <- list()

  for (age in age_grid) {
    nd <- bg
    nd$Forest_age <- age
    nd <- nd[, predictor_vars, drop = FALSE]

    if (model_type == "RANGER_QRF") {
      p <- predict_ranger_qrf_dataframe(model, nd, quantiles)
    } else if (model_type == "XGBOOST") {
      p <- predict_xgb_dataframe(model, nd)
    } else {
      stop("model_type 只能是 'RANGER_QRF' 或 'XGBOOST'。")
    }

    p <- postprocess_prediction_df(p, response_var)

    for (qn in names(p)) {
      z <- p[[qn]]
      rows[[length(rows) + 1]] <- data.frame(
        response = response_var,
        model_type = model_type,
        background = pdp_background,
        Forest_age = age,
        quantile_layer = qn,
        n_background = nrow(bg),
        pred_mean = mean(z, na.rm = TRUE),
        pred_median = stats::median(z, na.rm = TRUE),
        pred_q05 = as.numeric(stats::quantile(z, 0.05, na.rm = TRUE, names = FALSE)),
        pred_q95 = as.numeric(stats::quantile(z, 0.95, na.rm = TRUE, names = FALSE)),
        pred_min = min(z, na.rm = TRUE),
        pred_max = max(z, na.rm = TRUE)
      )
    }
  }

  do.call(rbind, rows)
}

if (run_forest_age_pdp) {
  message("==================================================")
  message("Running Forest_age PDP sanity check...")

  pdp_all <- list()

  for (resp in response_vars) {
    bg_dat <- train_list[[resp]]

    pdp_df <- make_forest_age_pdp(
      model = models[[resp]],
      response_var = resp,
      background_dat = bg_dat,
      age_grid = pdp_age_grid,
      sample_size = pdp_sample_size,
      seed = seed
    )

    pdp_all[[resp]] <- pdp_df

    write.csv(
      pdp_df,
      file.path(out_dir, paste0("pdp_Forest_age_", resp, ".csv")),
      row.names = FALSE
    )
  }

  pdp_combined <- do.call(rbind, pdp_all)
  write.csv(
    pdp_combined,
    file.path(out_dir, "pdp_Forest_age_all_responses.csv"),
    row.names = FALSE
  )
}

# ============================================================
# 14. 快速矩阵预测工具
# ============================================================

prepare_fast_prediction_data <- function(env_pred) {
  message("==================================================")
  message("Preparing fast prediction matrix from ENV raster...")

  # NOTE: when run_knndm_validation=TRUE, make_environmental_stratified_predpoints()
  # has already read the full ENV grid once to sample predpoints. We re-read it here.
  # For the target size (~1.78M cells x 7 layers, ~100 MB) this duplicate read is a
  # few seconds and negligible next to QRF fitting/prediction over the same cells,
  # so it is intentionally left as two independent, self-contained reads rather than
  # threading a shared cache through index-sensitive sampling code.
  vals <- values(env_pred, mat = TRUE)
  colnames(vals) <- names(env_pred)

  ok <- rowSums(!is.finite(vals)|vals == -9999)==0L
  ok_idx <- which(ok)

  base_df <- as.data.frame(vals[ok_idx, , drop = FALSE])
  names(base_df) <- names(env_pred)
  base_df$Forest_age <- NA_real_
  base_df <- base_df[, predictor_vars, drop = FALSE]

  rm(vals)
  gc()

  message("Valid raster cells for prediction: ", length(ok_idx), " / ", ncell(env_pred))

  return(list(
    ok_idx = ok_idx,
    base_df = base_df,
    template = env_pred[[1]]
  ))
}

prediction_df_to_raster <- function(pred_df,
                                    ok_idx,
                                    template,
                                    response_var,
                                    filename,
                                    range_mask = NULL,
                                    aoa_mask = NULL,
                                    already_postprocessed = FALSE) {
  nq <- ncol(pred_df)
  nc <- terra::ncell(template)

  M <- matrix(NA_real_, nrow = nc, ncol = nq)
  for (j in seq_len(nq)) {
    M[ok_idx, j] <- pred_df[[j]]
  }

  r <- terra::rast(template, nlyrs = nq)
  terra::values(r) <- M
  names(r) <- names(pred_df)

  # fast 路径中 pred_df 已经做过分位数顺序修正和非负裁剪，
  # 因此这里默认不重复做栅格级 max/ifel，减少一次无效计算。
  if (!already_postprocessed) {
    r <- postprocess_prediction_raster(r, response_var)
  } else {
    names(r) <- q_names
  }

  if (!is.null(range_mask) && apply_range_mask_to_predictions) {
    r <- mask(r, range_mask)
  }
  if (!is.null(aoa_mask) && apply_aoa_mask_to_predictions) {
    r <- mask(r, aoa_mask, maskvalues = 0)
  }

  write_raster_safe(r, filename)
}

predict_one_age_fast <- function(model, response_var, age_value, fast_data, filename, range_mask = NULL, aoa_mask = NULL) {
  nd <- fast_data$base_df
  nd$Forest_age <- age_value

  if (model_type == "RANGER_QRF") {
    p <- predict_ranger_qrf_dataframe(model, nd, quantiles)
  } else if (model_type == "XGBOOST") {
    p <- predict_xgb_dataframe(model, nd)
  }

  p <- postprocess_prediction_df(p, response_var)

  prediction_df_to_raster(
    pred_df = p,
    ok_idx = fast_data$ok_idx,
    template = fast_data$template,
    response_var = response_var,
    filename = filename,
    range_mask = range_mask,
    aoa_mask = aoa_mask,
    already_postprocessed = TRUE
  )
}

# ============================================================
# 15. terra::predict 分块预测兜底工具
#     注意：先写临时文件，再后处理写最终文件，避免 Windows 同名自覆盖
# ============================================================

# make_age_raster() and make_predictor_stack() are defined once above and reused here.

predict_one_age_terra <- function(model, response_var, age_value, filename, range_mask = NULL, aoa_mask = NULL) {
  pred_stack <- make_predictor_stack(ENV_pred, age_value)
  tmp <- tempfile(pattern = "terra_predict_", fileext = ".tif")
  on.exit(unlink(tmp), add = TRUE)

  if (model_type == "RANGER_QRF") {
    pred_fun <- function(model, data) {
      data <- as.data.frame(data)
      data <- data[, predictor_vars, drop = FALSE]
      p <- predict_ranger_qrf_dataframe(model, data, quantiles)
      return(p)
    }

    r <- terra::predict(
      pred_stack,
      model,
      fun = pred_fun,
      na.rm = TRUE,
      filename = tmp,
      overwrite = TRUE,
      datatype = write_datatype,
      gdal = write_gdal
    )

  } else if (model_type == "XGBOOST") {
    pred_fun <- function(model, data) {
      data <- as.data.frame(data)
      data <- data[, predictor_vars, drop = FALSE]
      p <- predict_xgb_dataframe(model, data)
      return(p)
    }

    r <- terra::predict(
      pred_stack,
      model,
      fun = pred_fun,
      na.rm = TRUE,
      filename = tmp,
      overwrite = TRUE,
      datatype = write_datatype,
      gdal = write_gdal
    )
  }

  names(r) <- q_names
  r <- postprocess_prediction_raster(r, response_var)

  if (!is.null(range_mask) && apply_range_mask_to_predictions) {
    r <- mask(r, range_mask)
  }
  if (!is.null(aoa_mask) && apply_aoa_mask_to_predictions) {
    r <- mask(r, aoa_mask, maskvalues = 0)
  }

  write_raster_safe(r, filename)
}

# ============================================================
# 16. 不同 Forest_age 固定值下预测分位数曲面
# ============================================================

pred_rasters <- list()
fast_data <- NULL

if (use_fast_matrix_prediction) {
  fast_data <- prepare_fast_prediction_data(ENV_pred)
}

for (resp in response_vars) {
  pred_rasters[[resp]] <- list()

  for (age in age_values) {
    message("==================================================")
    message("Predicting response = ", resp, ", Forest_age = ", age)

    out_name <- file.path(
      out_dir,
      paste0(model_type, "_", resp, "_ForestAge_", age, "_quantile_surface.tif")
    )

    range_mask <- NULL
    if (make_range_mask && resp %in% names(range_masks)) {
      range_mask <- range_masks[[resp]]
    }

    aoa_mask <- NULL
    if (make_aoa && exists("aoa_results") && resp %in% names(aoa_results) &&
        as.character(age) %in% names(aoa_results[[resp]])) {
      aoa_mask <- aoa_results[[resp]][[as.character(age)]]$AOA
    }

    if (use_fast_matrix_prediction) {
      pred_r <- predict_one_age_fast(
        model = models[[resp]],
        response_var = resp,
        age_value = age,
        fast_data = fast_data,
        filename = out_name,
        range_mask = range_mask,
        aoa_mask = aoa_mask
      )
    } else {
      pred_r <- predict_one_age_terra(
        model = models[[resp]],
        response_var = resp,
        age_value = age,
        filename = out_name,
        range_mask = range_mask,
        aoa_mask = aoa_mask
      )
    }

    pred_rasters[[resp]][[as.character(age)]] <- pred_r
  }
}

# ============================================================
# 17. 以 Forest_age = 100 为参考，计算敏感性曲面
# ============================================================

for (resp in response_vars) {
  ref_r <- pred_rasters[[resp]][[as.character(ref_age)]]

  if (is.null(ref_r)) {
    stop("找不到参考年龄 ref_age = ", ref_age, " 的预测结果：", resp)
  }

  for (age in age_values) {
    if (age == ref_age) next

    message("==================================================")
    message("Calculating sensitivity: ", resp, ", age ", age, " minus age ", ref_age)

    target_r <- pred_rasters[[resp]][[as.character(age)]]

    sens_abs <- target_r - ref_r
    names(sens_abs) <- paste0(
      names(target_r),
      "_age", age,
      "_minus_age", ref_age
    )

    out_abs <- file.path(
      out_dir,
      paste0(model_type, "_", resp, "_sensitivity_age", age, "_minus_age", ref_age, ".tif")
    )

    write_raster_safe(sens_abs, out_abs)

    # 百分比变化，避免除以 0 或极小参考值
    eps <- 1e-6
    ref_safe <- ifel(abs(ref_r) > eps, ref_r, NA)
    sens_pct <- 100 * sens_abs / ref_safe

    names(sens_pct) <- paste0(
      names(target_r),
      "_pct_age", age,
      "_vs_age", ref_age
    )

    out_pct <- file.path(
      out_dir,
      paste0(model_type, "_", resp, "_sensitivity_percent_age", age, "_vs_age", ref_age, ".tif")
    )

    write_raster_safe(sens_pct, out_pct)
  }
}

# ============================================================
# 18. 保存模型对象和运行信息
# ============================================================

saveRDS(
  models,
  file.path(out_dir, paste0(model_type, "_final_models.rds"))
)

run_info <- list(
  package_profile = paste0("FEM_v", .fem_software_version),
  qrf_parameter_mode = qrf_parameter_mode,
  model_type = model_type,
  env_file = normalizePath(env_file),
  csv_file = normalizePath(csv_file),
  out_dir = normalizePath(out_dir),
  response_vars = response_vars,
  raster_predictor_vars = raster_predictor_vars,
  predictor_vars = predictor_vars,
  quantiles = quantiles,
  age_values = age_values,
  ref_age = ref_age,
  seed = seed,
  detected_cores = detected_cores,
  hardware_cores = hardware_cores,
  hardware_memory_gb = hardware_memory_gb,
  reserve_cores_for_os = reserve_cores_for_os,
  n_model_threads = n_model_threads,
  n_cores = n_cores,
  terra_tempdir = terra_tempdir,
  use_fast_matrix_prediction = use_fast_matrix_prediction,
  make_range_mask = make_range_mask,
  apply_range_mask_to_predictions = apply_range_mask_to_predictions,
  make_aoa = make_aoa,
  apply_aoa_mask_to_predictions = apply_aoa_mask_to_predictions,
  aoa_use_variable_importance = aoa_use_variable_importance,
  aoa_use_knndm_folds = aoa_use_knndm_folds,
  aoa_lpd = aoa_lpd,
  aoa_method = aoa_method,
  aoa_skip_existing = aoa_skip_existing,
  aoa_compute_per_age = aoa_compute_per_age,
  make_environment_only_aoa = make_environment_only_aoa,
  environment_only_aoa_use_variable_importance = environment_only_aoa_use_variable_importance,
  environment_only_aoa_reference_response = environment_only_aoa_reference_response,
  run_holdout_validation = run_holdout_validation,
  validation_fraction = validation_fraction,
  validation_repeats = validation_repeats,
  run_knndm_validation = run_knndm_validation,
  knndm_k = knndm_k,
  knndm_predpoints_sampling = knndm_predpoints_sampling,
  knndm_samplesize = knndm_samplesize,
  knndm_env_clusters = knndm_env_clusters,
  knndm_min_per_cluster = knndm_min_per_cluster,
  knndm_pca_components = knndm_pca_components,
  knndm_candidate_multiplier = knndm_candidate_multiplier,
  knndm_kmeans_nstart = knndm_kmeans_nstart,
  knndm_kmeans_iter_max = knndm_kmeans_iter_max,
  knndm_clustering = knndm_clustering,
  knndm_maxp = knndm_maxp,
  knndm_space = knndm_space,
  write_knndm_predpoints_gpkg = write_knndm_predpoints_gpkg,
  qrf_num_trees_validation = qrf_num_trees_validation,
  qrf_num_trees_final = qrf_num_trees_final,
  qrf_num_trees_final_candidates = qrf_num_trees_final_candidates,
  run_tree_count_sensitivity = run_tree_count_sensitivity,
  qrf_min_node_size = qrf_min_node_size,
  qrf_mtry = qrf_mtry,
  qrf_always_split_variables = qrf_always_split_variables,
  qrf_params_by_response = qrf_params_by_response,
  auto_load_tuned_params = auto_load_tuned_params,
  tuned_params_file = if (isTRUE(auto_load_tuned_params) && file.exists(tuned_params_file)) normalizePath(tuned_params_file) else NA_character_,
  run_forest_age_pdp = run_forest_age_pdp,
  pdp_age_grid = pdp_age_grid,
  pdp_sample_size = pdp_sample_size,
  pdp_background = pdp_background,
  terra_version = as.character(utils::packageVersion("terra")),
  ranger_version = if (requireNamespace("ranger", quietly = TRUE)) as.character(utils::packageVersion("ranger")) else NA_character_,
  CAST_version = if (requireNamespace("CAST", quietly = TRUE)) as.character(utils::packageVersion("CAST")) else NA_character_,
  sf_version = if (requireNamespace("sf", quietly = TRUE)) as.character(utils::packageVersion("sf")) else NA_character_,
  xgboost_version = if (requireNamespace("xgboost", quietly = TRUE)) as.character(utils::packageVersion("xgboost")) else NA_character_
)

surface_files <- file.path(out_dir,paste0("RANGER_QRF_",response_vars,"_ForestAge_100_quantile_surface.tif"))
if(any(!file.exists(surface_files)))stop("Age-100 surface missing at QRF completion")
run_info$input_md5 <- setNames(unname(tools::md5sum(c(env_file,csv_file,template_file))),c("ENV.tif","train4pot.csv","templ_1km.tif"))
run_info$age100_surface_md5 <- setNames(unname(tools::md5sum(surface_files)),response_vars)
saveRDS(run_info, file.path(out_dir, "run_info.rds"))

message("==================================================")
message("Done.")
message("Outputs written to: ", normalizePath(out_dir))
message("Main prediction variables used: ", paste(predictor_vars, collapse = ", "))
message("Elapsed: ", format(round(difftime(Sys.time(), t_start, units = "mins"), 2)))
