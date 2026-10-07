# =============================================================================
# Compute_fig23_stats.R
#
# COMPUTATION-ONLY module of the Fig. 2 / Fig. 3 pipeline.
#
# Reads the analysis results object (all_results.rds) and produces two
# statistical summary CSVs. For each
# (target, cross-validation, scenario) cell, BOTH the observation-level
# bootstrap and the 50 km spatial-block bootstrap are computed on the
# same paired observation-level out-of-fold predictions, so the two
# bootstrap methods use the same observations.
#
# Inputs
# ------
#   all_results.rds (in `output_dir`), with each entry providing:
#     $target_var, $cv_method, $pft, $data_type
#     $model_results[[scenario]]$raw_predictions (raw RF OOF, length n)
#     or $predictions with $prediction_type == "RF_raw"
#     $model_results[[scenario]]$observed        (length n)
#     $model_results[[scenario]]$coords          (or $coordinates) with lon/lat
#     $model_results[[scenario]]$success         (TRUE)
#
#   Coordinates are required for the spatial-block bootstrap. The obs-level
#   bootstrap path runs even if coords are missing.
#
# Outputs (written to `output_dir/`)
# -----------------------------------
#   fig2_panelAB_R2_summary.csv
#       Long format. One row per (target × CV × scenario × method).
#       Columns:
#         Data_Type, Target, PFT, CV_Method, Scenario, Model,
#         Method                 (obs_level / spatial_block_50km),
#         estimate, ci_lo, ci_hi, n_boot_used, n_obs, n_blocks,
#         Block_Size_km, CI_Method
#
#   fig3_panelCD_DeltaR2_summary.csv
#       Long format. One row per (target × CV × contrast × method).
#       Columns:
#         Data_Type, Target, PFT, CV_Method, Source, Removal,
#         Full, Reduced, Model, Method,
#         estimate, ci_lo, ci_hi, p, p_BH, n_boot_used, n_obs, n_blocks,
#         Block_Size_km, Test_Method
#
#   fig23_compute_log.txt
#       Free-form text summary of strata processed, n_blocks per stratum,
#       and any warnings.
#
# Both bootstrap methods use:
#   * B = N_BOOT (default 2000)
#   * Plus-one achieved-significance two-sided p-values (capped at 1)
#   * BH adjustment within (Target × PFT × CV_Method × Method × Contrast_Family)
#     where Contrast_Family is "Block_removal_8" for Panel C/D
#
# Usage (command line)
# --------------------
#   Rscript Compute_fig23_stats.R OUTPUT_DIR [N_BOOT] [BLOCK_KM]
#   Rscript Compute_fig23_stats.R /path/to/scenario_results 2000 50
#
# Usage (interactive R)
# ---------------------
#   FIG23_OUTPUT_DIR <- "/path/to/scenario_results"
#   N_BOOT_FIG23     <- 2000   # optional
#   BLOCK_KM_FIG23   <- 50     # optional
#   source("Compute_fig23_stats.R")
#
# Output filenames are retained for compatibility with existing result tables.
# =============================================================================

# -- Configuration ------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
output_dir <- if (exists("FIG23_OUTPUT_DIR", envir = .GlobalEnv)) {
  get("FIG23_OUTPUT_DIR", envir = .GlobalEnv)
} else if (length(args) >= 1) {
  args[1]
} else {
  stop("Supply OUTPUT_DIR as the first argument or set FIG23_OUTPUT_DIR in interactive R.", call. = FALSE)
}

n_boot <- if (exists("N_BOOT_FIG23", envir = .GlobalEnv)) {
  as.integer(get("N_BOOT_FIG23", envir = .GlobalEnv))
} else if (length(args) >= 2) {
  as.integer(args[2])
} else {
  2000L
}

block_size_km <- if (exists("BLOCK_KM_FIG23", envir = .GlobalEnv)) {
  as.numeric(get("BLOCK_KM_FIG23", envir = .GlobalEnv))
} else if (length(args) >= 3) {
  as.numeric(args[3])
} else {
  50
}

ci_level  <- 0.95
boot_seed <- 42L

# -- Constants (mirrored from manuscript scenario set) ------------------------

scenario_keep <- c(
  "S0_Env", "S1a_Prod", "S1b_Het", "S1c_Temp",
  "S2a_Prod_Het", "S2b_Prod_Temp", "S2c_Het_Temp",
  "S3_RS_Full", "S4_Integrated"
)

removal_def <- data.frame(
  Source  = c("RS",       "Env",
              "ProdHet",  "ProdTemp", "HetTemp",
              "Prod",     "Het",      "Temp"),
  Removal = c(
    "Remove RS from S4_Integrated",
    "Remove Env from S4_Integrated",
    "Remove Prod+Het from S3_RS_Full",
    "Remove Prod+Temp from S3_RS_Full",
    "Remove Het+Temp from S3_RS_Full",
    "Remove Prod from S3_RS_Full",
    "Remove Het from S3_RS_Full",
    "Remove Temp from S3_RS_Full"
  ),
  Full    = c("S4_Integrated", "S4_Integrated",
              "S3_RS_Full",    "S3_RS_Full",    "S3_RS_Full",
              "S3_RS_Full",    "S3_RS_Full",    "S3_RS_Full"),
  Reduced = c("S0_Env",        "S3_RS_Full",
              "S1c_Temp",      "S1b_Het",        "S1a_Prod",
              "S2c_Het_Temp",  "S2b_Prod_Temp",  "S2a_Prod_Het"),
  stringsAsFactors = FALSE
)

METHOD_OBS   <- "obs_level"
METHOD_BLOCK <- "spatial_block_50km"

# -- Tiny helpers -------------------------------------------------------------

`%||%` <- function(a, b) if (!is.null(a)) a else b
safe_numeric <- function(x) suppressWarnings(as.numeric(x))

extract_coords <- function(scn_res) {
  coords <- scn_res$coords %||% scn_res$coordinates
  if (is.null(coords)) return(NULL)
  coords <- as.data.frame(coords)
  if (ncol(coords) < 2) return(NULL)
  nms <- names(coords)
  lon_col <- if ("Lon_Export" %in% nms) "Lon_Export" else if ("lon" %in% tolower(nms)) nms[which(tolower(nms) == "lon")[1]] else nms[1]
  lat_col <- if ("Lat_Export" %in% nms) "Lat_Export" else if ("lat" %in% tolower(nms)) nms[which(tolower(nms) == "lat")[1]] else nms[2]
  data.frame(
    lon = safe_numeric(coords[[lon_col]]),
    lat = safe_numeric(coords[[lat_col]])
  )
}

# Privacy-preserving public release support.  If exact coordinates were removed
# from all_results_public_no_coords.rds, spatial-block bootstrap can still be
# reproduced using anonymized precomputed block labels.  Only group membership
# is required by the bootstrap; exact lon/lat are not required at this stage.
extract_block_id <- function(scn_res, block_size_km = 50) {
  candidates <- unique(c(
    sprintf("block_id_%dkm", as.integer(block_size_km)),
    sprintf("spatial_block_%dkm", as.integer(block_size_km)),
    "block_id_50km",
    "spatial_block_50km",
    "spatial_block_id_anonymized",
    "block_id"
  ))
  for (nm in candidates) {
    if (!is.null(scn_res[[nm]])) {
      x <- as.character(scn_res[[nm]])
      if (length(x) > 0 && any(!is.na(x))) return(x)
    }
  }
  NULL
}

make_spatial_blocks <- function(lon, lat, block_size_km = 50) {
  if (length(lon) != length(lat)) stop("lon and lat lengths do not match.")
  if (any(!is.finite(lon)) || any(!is.finite(lat))) {
    stop("Non-finite coordinates after complete-case filtering.")
  }
  lon0 <- stats::median(lon, na.rm = TRUE)
  lat0 <- stats::median(lat, na.rm = TRUE)
  x_km <- (lon - lon0) * 111.32 * cos(lat0 * pi / 180)
  y_km <- (lat - lat0) * 110.57
  paste(floor(x_km / block_size_km), floor(y_km / block_size_km), sep = "_")
}

pooled_r2 <- function(obs, pred) {
  ss_tot <- sum((obs - mean(obs))^2)
  if (!is.finite(ss_tot) || ss_tot <= 0) return(NA_real_)
  1 - sum((obs - pred)^2) / ss_tot
}

summarise_boot <- function(estimate, boot_values, ci_level = 0.95,
                           include_p = FALSE) {
  boot_values <- boot_values[is.finite(boot_values)]
  B <- length(boot_values)
  out <- list(estimate = estimate, ci_lo = NA_real_, ci_hi = NA_real_,
              n_boot_used = B)
  if (include_p) out$p <- NA_real_
  if (B < 100) return(out)
  alpha <- 1 - ci_level
  ci <- stats::quantile(boot_values, c(alpha / 2, 1 - alpha / 2), na.rm = TRUE)
  out$ci_lo <- unname(ci[1])
  out$ci_hi <- unname(ci[2])
  if (include_p) {
    p_left  <- (sum(boot_values <= 0) + 1) / (B + 1)
    p_right <- (sum(boot_values >= 0) + 1) / (B + 1)
    out$p <- min(2 * min(p_left, p_right), 1)
  }
  out
}

# -- Bootstrap primitives -----------------------------------------------------

# Single-scenario obs-level bootstrap of pooled R^2.
boot_r2_obs <- function(obs, pred, n_boot = 2000, ci_level = 0.95, seed = NULL) {
  ok <- is.finite(obs) & is.finite(pred); obs <- obs[ok]; pred <- pred[ok]
  n <- length(obs)
  if (n < 10) return(list(estimate = NA_real_, ci_lo = NA_real_,
                          ci_hi = NA_real_, n_boot_used = 0L,
                          n_obs = n, n_blocks = NA_integer_))
  estimate <- pooled_r2(obs, pred)
  if (!is.null(seed)) set.seed(seed)
  boot_idx <- matrix(sample.int(n, n * n_boot, replace = TRUE), nrow = n_boot)
  obs_b <- matrix(obs[boot_idx],  nrow = n_boot)
  prd_b <- matrix(pred[boot_idx], nrow = n_boot)
  mn <- rowMeans(obs_b)
  ss_tot <- rowSums((obs_b - mn)^2)
  r2_b <- ifelse(ss_tot > 0, 1 - rowSums((obs_b - prd_b)^2) / ss_tot, NA_real_)
  s <- summarise_boot(estimate, r2_b, ci_level)
  c(s, list(n_obs = n, n_blocks = NA_integer_))
}

# Single-scenario 50 km spatial-block bootstrap of pooled R^2.
boot_r2_block <- function(obs, pred, block_id, n_boot = 2000, ci_level = 0.95,
                          seed = NULL) {
  # Bug 1 fix: validate block_id length matches obs length BEFORE filtering.
  # Mismatch indicates that the caller passed a block_id derived from a
  # different scenario whose observation order or count does not match this
  # scenario's. Continuing would silently produce wrong CIs via R recycling.
  if (length(block_id) != length(obs)) {
    stop(sprintf("boot_r2_block: block_id length (%d) does not match obs length (%d). ",
                 length(block_id), length(obs)),
         "Caller must pass block_id matched to this scenario's observation order.")
  }
  ok <- is.finite(obs) & is.finite(pred); obs <- obs[ok]; pred <- pred[ok]
  block_id <- block_id[ok]
  n <- length(obs)
  if (n < 10) return(list(estimate = NA_real_, ci_lo = NA_real_,
                          ci_hi = NA_real_, n_boot_used = 0L,
                          n_obs = n, n_blocks = 0L))
  estimate <- pooled_r2(obs, pred)
  ub <- sort(unique(block_id))
  nb <- length(ub)
  # Bug 2 fix: warn when block count is too small for stable bootstrap.
  if (nb < 5) {
    warning(sprintf("boot_r2_block: only %d spatial blocks; bootstrap CI may be unstable.",
                    nb), call. = FALSE)
  }
  bidx <- split(seq_len(n), block_id)
  if (!is.null(seed)) set.seed(seed)
  r2_b <- numeric(n_boot)
  for (b in seq_len(n_boot)) {
    samp <- sample(ub, nb, replace = TRUE)
    idx <- unlist(bidx[samp], use.names = FALSE)
    obs_b <- obs[idx]; prd_b <- pred[idx]
    sst <- sum((obs_b - mean(obs_b))^2)
    r2_b[b] <- if (!is.finite(sst) || sst <= 0) NA_real_
               else 1 - sum((obs_b - prd_b)^2) / sst
  }
  s <- summarise_boot(estimate, r2_b, ci_level)
  c(s, list(n_obs = n, n_blocks = nb))
}

# Paired obs-level bootstrap of (R^2_full - R^2_reduced).
boot_diff_obs <- function(obs, pred_full, pred_red, n_boot = 2000,
                          ci_level = 0.95, seed = NULL) {
  ok <- is.finite(obs) & is.finite(pred_full) & is.finite(pred_red)
  obs <- obs[ok]; pred_full <- pred_full[ok]; pred_red <- pred_red[ok]
  n <- length(obs)
  if (n < 10) return(list(estimate = NA_real_, ci_lo = NA_real_,
                          ci_hi = NA_real_, p = NA_real_,
                          n_boot_used = 0L, n_obs = n, n_blocks = NA_integer_))
  sst <- sum((obs - mean(obs))^2)
  if (!is.finite(sst) || sst <= 0)
    return(list(estimate = NA_real_, ci_lo = NA_real_, ci_hi = NA_real_,
                p = NA_real_, n_boot_used = 0L, n_obs = n,
                n_blocks = NA_integer_))
  estimate <- (sum((obs - pred_red)^2) - sum((obs - pred_full)^2)) / sst
  if (!is.null(seed)) set.seed(seed)
  bi <- matrix(sample.int(n, n * n_boot, replace = TRUE), nrow = n_boot)
  ob <- matrix(obs[bi], nrow = n_boot)
  fb <- matrix(pred_full[bi], nrow = n_boot)
  rb <- matrix(pred_red[bi],  nrow = n_boot)
  mn <- rowMeans(ob)
  sst_b <- rowSums((ob - mn)^2)
  diffs <- ifelse(sst_b > 0,
                  (rowSums((ob - rb)^2) - rowSums((ob - fb)^2)) / sst_b,
                  NA_real_)
  s <- summarise_boot(estimate, diffs, ci_level, include_p = TRUE)
  c(s, list(n_obs = n, n_blocks = NA_integer_))
}

# Paired 50 km spatial-block bootstrap of (R^2_full - R^2_reduced).
boot_diff_block <- function(obs, pred_full, pred_red, block_id,
                            n_boot = 2000, ci_level = 0.95, seed = NULL) {
  # Bug 1 fix: validate block_id length matches obs length BEFORE filtering.
  if (length(block_id) != length(obs)) {
    stop(sprintf("boot_diff_block: block_id length (%d) does not match obs length (%d). ",
                 length(block_id), length(obs)),
         "Caller must pass block_id matched to the contrast's observation order.")
  }
  ok <- is.finite(obs) & is.finite(pred_full) & is.finite(pred_red)
  obs <- obs[ok]; pred_full <- pred_full[ok]; pred_red <- pred_red[ok]
  block_id <- block_id[ok]
  n <- length(obs)
  if (n < 10) return(list(estimate = NA_real_, ci_lo = NA_real_,
                          ci_hi = NA_real_, p = NA_real_,
                          n_boot_used = 0L, n_obs = n, n_blocks = 0L))
  sst <- sum((obs - mean(obs))^2)
  if (!is.finite(sst) || sst <= 0)
    return(list(estimate = NA_real_, ci_lo = NA_real_, ci_hi = NA_real_,
                p = NA_real_, n_boot_used = 0L, n_obs = n, n_blocks = 0L))
  estimate <- (sum((obs - pred_red)^2) - sum((obs - pred_full)^2)) / sst
  ub <- sort(unique(block_id))
  nb <- length(ub)
  # Bug 2 fix: warn when block count is too small for stable bootstrap.
  if (nb < 5) {
    warning(sprintf("boot_diff_block: only %d spatial blocks; bootstrap CI may be unstable.",
                    nb), call. = FALSE)
  }
  bidx <- split(seq_len(n), block_id)
  if (!is.null(seed)) set.seed(seed)
  diffs <- numeric(n_boot)
  for (b in seq_len(n_boot)) {
    samp <- sample(ub, nb, replace = TRUE)
    idx <- unlist(bidx[samp], use.names = FALSE)
    ob <- obs[idx]; fb <- pred_full[idx]; rb <- pred_red[idx]
    sst_b <- sum((ob - mean(ob))^2)
    diffs[b] <- if (!is.finite(sst_b) || sst_b <= 0) NA_real_
                else (sum((ob - rb)^2) - sum((ob - fb)^2)) / sst_b
  }
  s <- summarise_boot(estimate, diffs, ci_level, include_p = TRUE)
  c(s, list(n_obs = n, n_blocks = nb))
}

# -- Main entry point ---------------------------------------------------------

rds_path <- file.path(output_dir, "all_results.rds")
if (!file.exists(rds_path)) {
  stop("Could not find all_results.rds in: ", output_dir)
}
log_lines <- character()
add_log <- function(...) {
  msg <- paste0(...)
  cat(msg, "\n", sep = "")
  log_lines[length(log_lines) + 1L] <<- msg
}

add_log("[compute_fig23] Loading: ", rds_path)
all_results <- readRDS(rds_path)
# Validate before any paired bootstrap or publication output.
.fem_sources <- unlist(lapply(sys.frames(),function(e)e$ofile),use.names=FALSE)
.fem_script <- if(length(.fem_sources))tail(.fem_sources,1L) else sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1])
source(file.path(dirname(normalizePath(.fem_script)),"..","..","R","oof_contracts.R"))
rm(.fem_sources,.fem_script)
all_results <- fem_use_raw_rf_oof(all_results)
add_log("[compute_fig23] Prediction type: RF_raw (uncalibrated outer-test OOF); QM excluded.")


panel_ab_rows <- list()
panel_cd_rows <- list()
ab_ctr <- 0L

for (key in names(all_results)) {
  res <- all_results[[key]]
  if (is.null(res$model_results)) next
  target    <- res$target_var
  pft       <- res$pft
  cv_method <- res$cv_method %||% "kNNDM"
  data_type <- res$data_type %||% "match_only"

  add_log(sprintf("[compute_fig23] stratum: %s | %s | %s", target, pft, cv_method))

  scenarios_here <- intersect(scenario_keep, names(res$model_results))

  # Spatial-block groups once per stratum.
  # Prefer anonymized precomputed block IDs when available, so a public
  # coordinate-free all_results object can reproduce the spatial-block bootstrap.
  block_id <- NULL
  coords <- NULL
  for (scn0 in scenarios_here) {
    sr0 <- res$model_results[[scn0]]
    if (!isTRUE(sr0$success)) next
    block_id <- extract_block_id(sr0, block_size_km = block_size_km)
    if (!is.null(block_id) && length(block_id) >= 10) {
      add_log(sprintf("  Using precomputed anonymized block_id for spatial-block bootstrap."))
      break
    }
  }
  if (is.null(block_id)) {
    for (scn0 in scenarios_here) {
      sr0 <- res$model_results[[scn0]]
      if (!isTRUE(sr0$success)) next
      coords <- extract_coords(sr0)
      if (!is.null(coords) && nrow(coords) >= 10) break
    }
    block_id <- if (!is.null(coords) && nrow(coords) >= 10) {
      make_spatial_blocks(coords$lon, coords$lat, block_size_km = block_size_km)
    } else {
      add_log(sprintf("  WARNING: no coords or block_id; spatial-block bootstrap will be skipped for this stratum."))
      NULL
    }
  }

  # ---------- Panel A/B: 9 scenarios ----------
  for (scn in scenarios_here) {
    sr <- res$model_results[[scn]]
    if (!isTRUE(sr$success))     next
    if (is.null(sr$predictions)) next
    if (is.null(sr$observed))    next

    ab_ctr <- ab_ctr + 1L
    obs_seed   <- boot_seed + ab_ctr
    block_seed <- boot_seed + 100000L + ab_ctr

    # obs-level
    s_obs <- boot_r2_obs(sr$observed, sr$predictions,
                         n_boot = n_boot, ci_level = ci_level, seed = obs_seed)
    panel_ab_rows[[length(panel_ab_rows) + 1]] <- data.frame(
      Data_Type = data_type, Target = target, PFT = pft, CV_Method = cv_method,
      Scenario = scn, Model = "RF", Method = METHOD_OBS,
      estimate = s_obs$estimate, ci_lo = s_obs$ci_lo, ci_hi = s_obs$ci_hi,
      n_boot_used = s_obs$n_boot_used, n_obs = s_obs$n_obs,
      n_blocks = NA_integer_, Block_Size_km = NA_real_,
      CI_Method = sprintf("Observation-level percentile bootstrap (B=%d)", n_boot),
      stringsAsFactors = FALSE
    )

    # spatial-block (only if coords available)
    if (!is.null(block_id)) {
      s_blk <- boot_r2_block(sr$observed, sr$predictions, block_id,
                             n_boot = n_boot, ci_level = ci_level, seed = block_seed)
      panel_ab_rows[[length(panel_ab_rows) + 1]] <- data.frame(
        Data_Type = data_type, Target = target, PFT = pft, CV_Method = cv_method,
        Scenario = scn, Model = "RF", Method = METHOD_BLOCK,
        estimate = s_blk$estimate, ci_lo = s_blk$ci_lo, ci_hi = s_blk$ci_hi,
        n_boot_used = s_blk$n_boot_used, n_obs = s_blk$n_obs,
        n_blocks = s_blk$n_blocks, Block_Size_km = block_size_km,
        CI_Method = sprintf("50 km spatial-block percentile bootstrap (B=%d)", n_boot),
        stringsAsFactors = FALSE
      )
    }
  }

  # ---------- Panel C/D: 8 source-removal contrasts ----------
  for (i in seq_len(nrow(removal_def))) {
    cd <- removal_def[i, ]
    fr <- res$model_results[[cd$Full]]
    rr <- res$model_results[[cd$Reduced]]
    if (is.null(fr) || !isTRUE(fr$success)) next
    if (is.null(rr) || !isTRUE(rr$success)) next
    if (is.null(fr$predictions) || is.null(rr$predictions)) next
    if (is.null(fr$observed)    || is.null(rr$observed))    next
    if (length(fr$observed) != length(rr$observed) ||
        !isTRUE(all.equal(fr$observed, rr$observed, tolerance = 1e-8))) {
      add_log(sprintf("  Skipping contrast %s: observed values misaligned.", cd$Source))
      next
    }
    obs_seed   <- boot_seed + 1000L + i
    block_seed <- boot_seed + 200000L + i

    # obs-level
    s_obs <- boot_diff_obs(fr$observed, fr$predictions, rr$predictions,
                           n_boot = n_boot, ci_level = ci_level, seed = obs_seed)
    panel_cd_rows[[length(panel_cd_rows) + 1]] <- data.frame(
      Data_Type = data_type, Target = target, PFT = pft, CV_Method = cv_method,
      Source = cd$Source, Removal = cd$Removal, Full = cd$Full, Reduced = cd$Reduced,
      Model = "RF", Method = METHOD_OBS,
      estimate = s_obs$estimate, ci_lo = s_obs$ci_lo, ci_hi = s_obs$ci_hi,
      p = s_obs$p, n_boot_used = s_obs$n_boot_used, n_obs = s_obs$n_obs,
      n_blocks = NA_integer_, Block_Size_km = NA_real_,
      Test_Method = sprintf("Paired observation-level bootstrap (B=%d)", n_boot),
      stringsAsFactors = FALSE
    )

    # spatial-block
    if (!is.null(block_id)) {
      s_blk <- boot_diff_block(fr$observed, fr$predictions, rr$predictions,
                               block_id, n_boot = n_boot, ci_level = ci_level,
                               seed = block_seed)
      panel_cd_rows[[length(panel_cd_rows) + 1]] <- data.frame(
        Data_Type = data_type, Target = target, PFT = pft, CV_Method = cv_method,
        Source = cd$Source, Removal = cd$Removal, Full = cd$Full, Reduced = cd$Reduced,
        Model = "RF", Method = METHOD_BLOCK,
        estimate = s_blk$estimate, ci_lo = s_blk$ci_lo, ci_hi = s_blk$ci_hi,
        p = s_blk$p, n_boot_used = s_blk$n_boot_used, n_obs = s_blk$n_obs,
        n_blocks = s_blk$n_blocks, Block_Size_km = block_size_km,
        Test_Method = sprintf("Paired 50 km spatial-block bootstrap (B=%d)", n_boot),
        stringsAsFactors = FALSE
      )
    }
  }
}

if (length(panel_ab_rows) == 0) stop("No panel A/B rows produced.")
if (length(panel_cd_rows) == 0) stop("No panel C/D rows produced.")

panel_ab <- do.call(rbind, panel_ab_rows)
panel_cd <- do.call(rbind, panel_cd_rows)

# BH adjustment within (Target × PFT × CV_Method × Method) for the
# 8 contrast Block_removal family. p_BH is meaningful only for panel CD.
panel_cd$p_BH <- ave(
  panel_cd$p,
  panel_cd$Target, panel_cd$PFT, panel_cd$CV_Method, panel_cd$Method,
  FUN = function(x) stats::p.adjust(x, method = "BH")
)

# Reorder columns for readability.
panel_ab$Prediction_Type <- "RF_raw"
panel_cd$Prediction_Type <- "RF_raw"
ab_cols <- c("Data_Type", "Target", "PFT", "CV_Method", "Scenario", "Model", "Prediction_Type",
             "Method", "estimate", "ci_lo", "ci_hi",
             "n_boot_used", "n_obs", "n_blocks", "Block_Size_km", "CI_Method")
cd_cols <- c("Data_Type", "Target", "PFT", "CV_Method", "Source", "Removal",
             "Full", "Reduced", "Model", "Prediction_Type", "Method",
             "estimate", "ci_lo", "ci_hi", "p", "p_BH",
             "n_boot_used", "n_obs", "n_blocks", "Block_Size_km", "Test_Method")
panel_ab <- panel_ab[, ab_cols, drop = FALSE]
panel_cd <- panel_cd[, cd_cols, drop = FALSE]

# Outputs
out_ab  <- file.path(output_dir, "fig2_panelAB_R2_summary.csv")
out_cd  <- file.path(output_dir, "fig3_panelCD_DeltaR2_summary.csv")
out_log <- file.path(output_dir, "fig23_compute_log.txt")
write.csv(panel_ab, out_ab, row.names = FALSE, fileEncoding = "UTF-8")
write.csv(panel_cd, out_cd, row.names = FALSE, fileEncoding = "UTF-8")

add_log("")
add_log("[compute_fig23] Done.")
add_log(sprintf("  Panel A/B rows : %d (= n_strata x n_scenarios x n_methods)", nrow(panel_ab)))
add_log(sprintf("  Panel C/D rows : %d (= n_strata x 8 contrasts x n_methods)", nrow(panel_cd)))
add_log(sprintf("  fig2 CSV       : %s", out_ab))
add_log(sprintf("  fig3 CSV       : %s", out_cd))
writeLines(log_lines, con = out_log, useBytes = TRUE)
add_log(sprintf("  log            : %s", out_log))

# Quick diagnostic: show centerpiece kNNDM rows for both methods.
cat("\n=== Centerpiece check: source removal under kNNDM ===\n")
key_subset <- subset(panel_cd,
                     CV_Method == "kNNDM" & Source %in% c("RS", "Prod", "HetTemp"))
print(key_subset[, c("Target", "Source", "Method", "estimate", "ci_lo",
                     "ci_hi", "p_BH")],
      row.names = FALSE)
