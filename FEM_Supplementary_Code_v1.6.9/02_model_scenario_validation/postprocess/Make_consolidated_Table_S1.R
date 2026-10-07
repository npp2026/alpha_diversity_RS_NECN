# =============================================================================
# Make_consolidated_Table_S1.R
#
# Purpose
#   Summarize paired model contrasts for the supplementary tables.
#   Legacy output filenames are retained. Block-removal outputs
#   support current Table S3; substitution outputs support current Table S4.
#
#     Section 0: Block_removal_8
#                 3 single-block + 3 pair-block + 1 all-RS + 1 abiotic removal.
#
#     Section A: RS_marginal_9 + Env_RS_complementarity_2
#                 M(Block | Context) = R2(Context + Block) - R2(Context).
#
#     Section B: RS_substitution_3
#                 S(A, B | C) = M(B | C) - M(B | A + C)
#                             = [R2(C+B)-R2(C)] - [R2(A+B+C)-R2(A+C)].
#
#   This script reads all_results.rds directly. It does NOT retrain models.
#   It uses raw RF paired out-of-fold predictions from the strict 9-scenario factorial
#   design and computes pooled OOF R2 under a shared observation set.
#
# Main outputs
#   1. Table_S1_consolidated_block_removal.csv       (Section 0; 8 × 4 = 32 rows)
#   2. Table_S1_consolidated_marginal.csv            (Section A; 44 rows)
#   3. Table_S1_consolidated_substitution.csv        (Section B; 12 rows)
#   4. Table_S1_consolidated_all.csv                 (combined; 88 rows)
#
# Statistical inference
#   Paired spatial-block percentile bootstrap, B = 2000, block size = 50 km.
#   Blocks are resampled with replacement. All scenarios are evaluated on the
#   same resampled observations within each bootstrap replicate.
#   BH adjustment is done within
#     Data_Type x Target x PFT x CV_Method x Contrast_Family
#   with FOUR contrast families (Block_removal_8, RS_marginal_9,
#   Env_RS_complementarity_2, RS_substitution_3) adjusted independently.
#
# Usage
#   Rscript Make_consolidated_Table_S1.R /path/to/scenario_results 2000 50
#
# Optional in interactive R before source():
#   FIG4_OUTPUT_DIR <- "/path/to/scenario_results"
#   N_BOOT_TABLE_S1 <- 2000
#   BLOCK_SIZE_KM_TABLE_S1 <- 50
#   SEED_TABLE_S1 <- 42
#
# =============================================================================

options(stringsAsFactors = FALSE)

`%||%` <- function(a, b) if (!is.null(a)) a else b

# -- Configuration -------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

output_dir <- if (exists("FIG4_OUTPUT_DIR", envir = .GlobalEnv)) {
  get("FIG4_OUTPUT_DIR", envir = .GlobalEnv)
} else if (length(args) >= 1) {
  args[1]
} else {
  stop("Supply OUTPUT_DIR as the first argument or set FIG4_OUTPUT_DIR in interactive R.", call. = FALSE)
}

.get_opt <- function(name_s1, name_s2_legacy, fallback) {
  if (exists(name_s1, envir = .GlobalEnv)) return(get(name_s1, envir = .GlobalEnv))
  if (exists(name_s2_legacy, envir = .GlobalEnv)) return(get(name_s2_legacy, envir = .GlobalEnv))
  fallback
}

n_boot <- as.integer(.get_opt("N_BOOT_TABLE_S1", "N_BOOT_TABLE_S2",
                              if (length(args) >= 2) as.integer(args[2]) else 2000L))

block_size_km <- as.numeric(.get_opt("BLOCK_SIZE_KM_TABLE_S1", "BLOCK_SIZE_KM_TABLE_S2",
                                     if (length(args) >= 3) as.numeric(args[3]) else 50))

ci_level <- as.numeric(.get_opt("CI_LEVEL_TABLE_S1", "CI_LEVEL_TABLE_S2", 0.95))

seed <- as.integer(.get_opt("SEED_TABLE_S1", "SEED_TABLE_S2", 42L))

allow_observation_bootstrap_fallback <- if (exists("ALLOW_OBS_BOOT_FALLBACK", envir = .GlobalEnv)) {
  isTRUE(get("ALLOW_OBS_BOOT_FALLBACK", envir = .GlobalEnv))
} else {
  FALSE
}

output_subdir <- file.path(output_dir, "Table_S1_consolidated")
dir.create(output_subdir, showWarnings = FALSE, recursive = TRUE)

scenario_keep <- c(
  "S0_Env", "S1a_Prod", "S1b_Het", "S1c_Temp",
  "S2a_Prod_Het", "S2b_Prod_Temp", "S2c_Het_Temp",
  "S3_RS_Full", "S4_Integrated"
)

# -- Contrast definitions ------------------------------------------------------

# Section 0 (NEW in v2.6): Block-removal contrasts.
# Mathematically these are the same form as marginal contributions
# (R2_full - R2_reduced), but framed as "removing a block from a larger model".
removal_def <- data.frame(
  Table_Panel = "S1_section_0",
  Contrast_Family = "Block_removal_8",
  Contrast_ID = c(
    "R_AllRS_from_Integrated",
    "R_Env_from_Integrated",
    "R_Prod_from_RSFull",
    "R_Het_from_RSFull",
    "R_Temp_from_RSFull",
    "R_ProdHet_from_RSFull",
    "R_ProdTemp_from_RSFull",
    "R_HetTemp_from_RSFull"
  ),
  Block_Removed = c(
    "All_RS",
    "Env",
    "Prod",
    "Het",
    "Temp",
    "Prod+Het",
    "Prod+Temp",
    "Het+Temp"
  ),
  Context_Retained = c(
    "Env_only",
    "RS_Full_only",
    "Het+Temp",
    "Prod+Temp",
    "Prod+Het",
    "Temp",
    "Het",
    "Prod"
  ),
  Reduced = c(
    "S0_Env",
    "S3_RS_Full",
    "S2c_Het_Temp",
    "S2b_Prod_Temp",
    "S2a_Prod_Het",
    "S1c_Temp",
    "S1b_Het",
    "S1a_Prod"
  ),
  Full = c(
    "S4_Integrated",
    "S4_Integrated",
    "S3_RS_Full",
    "S3_RS_Full",
    "S3_RS_Full",
    "S3_RS_Full",
    "S3_RS_Full",
    "S3_RS_Full"
  ),
  Interpretation = c(
    "All remote-sensing block removed from the integrated model",
    "Abiotic block removed from the integrated model",
    "Productivity block removed from the full RS model",
    "Heterogeneity block removed from the full RS model",
    "Temporal-dynamics block removed from the full RS model",
    "Productivity + Heterogeneity pair removed from the full RS model",
    "Productivity + Temporal pair removed from the full RS model",
    "Heterogeneity + Temporal pair removed from the full RS model"
  ),
  stringsAsFactors = FALSE
)

marginal_def <- data.frame(
  Table_Panel = "S1_section_A",
  Contrast_Family = c(rep("RS_marginal_9", 9), rep("Env_RS_complementarity_2", 2)),
  Contrast_ID = c(
    "M_Prod_given_Het",
    "M_Prod_given_Temp",
    "M_Prod_given_HetTemp",
    "M_Het_given_Prod",
    "M_Het_given_Temp",
    "M_Het_given_ProdTemp",
    "M_Temp_given_Prod",
    "M_Temp_given_Het",
    "M_Temp_given_ProdHet",
    "M_Env_given_RSFull",
    "M_RSFull_given_Env"
  ),
  Block_Added = c(
    "Prod", "Prod", "Prod",
    "Het", "Het", "Het",
    "Temp", "Temp", "Temp",
    "Env", "RS_Full"
  ),
  Context = c(
    "Het", "Temp", "Het+Temp",
    "Prod", "Temp", "Prod+Temp",
    "Prod", "Het", "Prod+Het",
    "Prod+Het+Temp", "Env"
  ),
  Reduced = c(
    "S1b_Het",
    "S1c_Temp",
    "S2c_Het_Temp",
    "S1a_Prod",
    "S1c_Temp",
    "S2b_Prod_Temp",
    "S1a_Prod",
    "S1b_Het",
    "S2a_Prod_Het",
    "S3_RS_Full",
    "S0_Env"
  ),
  Full = c(
    "S2a_Prod_Het",
    "S2b_Prod_Temp",
    "S3_RS_Full",
    "S2a_Prod_Het",
    "S2c_Het_Temp",
    "S3_RS_Full",
    "S2b_Prod_Temp",
    "S2c_Het_Temp",
    "S3_RS_Full",
    "S4_Integrated",
    "S4_Integrated"
  ),
  Interpretation = c(
    "Productivity marginal contribution when heterogeneity is present",
    "Productivity marginal contribution when temporal dynamics is present",
    "Productivity marginal contribution when heterogeneity and temporal dynamics are present",
    "Heterogeneity marginal contribution when productivity is present",
    "Heterogeneity marginal contribution when temporal dynamics is present",
    "Heterogeneity marginal contribution when productivity and temporal dynamics are present",
    "Temporal-dynamics marginal contribution when productivity is present",
    "Temporal-dynamics marginal contribution when heterogeneity is present",
    "Temporal-dynamics marginal contribution when productivity and heterogeneity are present",
    "Abiotic marginal contribution when full remote sensing is present",
    "Full remote-sensing marginal contribution when abiotic variables are present"
  )
)

substitution_def <- data.frame(
  Table_Panel = "S1_section_B",
  Contrast_Family = "RS_substitution_3",
  Contrast_ID = c(
    "S_Prod_Temp_given_Het",
    "S_Prod_Het_given_Temp",
    "S_Het_Temp_given_Prod"
  ),
  Substitution_Pair = c("Prod-Temp", "Prod-Het", "Het-Temp"),
  Present_Block_A = c("Prod", "Prod", "Het"),
  Marginal_Block_B = c("Temp", "Het", "Temp"),
  Conditioning_Context_C = c("Het", "Temp", "Prod"),
  Scenario_C = c("S1b_Het", "S1c_Temp", "S1a_Prod"),
  Scenario_C_plus_B = c("S2c_Het_Temp", "S2c_Het_Temp", "S2b_Prod_Temp"),
  Scenario_A_plus_C = c("S2a_Prod_Het", "S2b_Prod_Temp", "S2a_Prod_Het"),
  Scenario_A_plus_B_plus_C = c("S3_RS_Full", "S3_RS_Full", "S3_RS_Full"),
  Interpretation = c(
    "Reduction in Temp marginal value when Prod is also present, conditioned on Het",
    "Reduction in Het marginal value when Prod is also present, conditioned on Temp",
    "Reduction in Temp marginal value when Het is also present, conditioned on Prod"
  )
)

# -- Helper functions ----------------------------------------------------------

safe_numeric <- function(x) {
  suppressWarnings(as.numeric(x))
}

pooled_r2 <- function(obs, pred) {
  ok <- is.finite(obs) & is.finite(pred)
  obs <- obs[ok]
  pred <- pred[ok]
  if (length(obs) < 10) return(NA_real_)
  ss_tot <- sum((obs - mean(obs))^2)
  if (!is.finite(ss_tot) || ss_tot <= 0) return(NA_real_)
  1 - sum((obs - pred)^2) / ss_tot
}

extract_predictions <- function(scn_res) {
  if (!is.null(scn_res$predictions)) return(safe_numeric(scn_res$predictions))
  if (!is.null(scn_res$cv_predictions)) return(safe_numeric(scn_res$cv_predictions))
  NULL
}

extract_coords <- function(scn_res) {
  coords <- scn_res$coords %||% scn_res$coordinates
  if (is.null(coords)) return(NULL)
  coords <- as.data.frame(coords)
  if (ncol(coords) < 2) return(NULL)

  nms <- names(coords)
  lon_col <- if ("Lon_Export" %in% nms) "Lon_Export" else if ("lon" %in% tolower(nms)) nms[which(tolower(nms) == "lon")[1]] else nms[1]
  lat_col <- if ("Lat_Export" %in% nms) "Lat_Export" else if ("lat" %in% tolower(nms)) nms[which(tolower(nms) == "lat")[1]] else nms[2]

  out <- data.frame(
    lon = safe_numeric(coords[[lon_col]]),
    lat = safe_numeric(coords[[lat_col]])
  )
  out
}

make_spatial_blocks <- function(lon, lat, block_size_km = 50) {
  if (length(lon) != length(lat)) stop("lon and lat lengths do not match.")
  if (any(!is.finite(lon)) || any(!is.finite(lat))) {
    stop("Non-finite coordinates found after complete-case filtering.")
  }

  lon0 <- stats::median(lon, na.rm = TRUE)
  lat0 <- stats::median(lat, na.rm = TRUE)

  # Equirectangular approximation around study-area median latitude.
  # This is adequate for 50 km blocking over Northeast China.
  x_km <- (lon - lon0) * 111.32 * cos(lat0 * pi / 180)
  y_km <- (lat - lat0) * 110.57

  bx <- floor(x_km / block_size_km)
  by <- floor(y_km / block_size_km)
  paste(bx, by, sep = "_")
}

build_prediction_table <- function(res, scenario_keep) {
  model_results <- res$model_results
  if (is.null(model_results)) stop("Result has no model_results.")

  missing <- scenario_keep[!scenario_keep %in% names(model_results)]
  if (length(missing) > 0) {
    stop(sprintf("Missing required scenarios: %s", paste(missing, collapse = ", ")))
  }

  success <- vapply(scenario_keep, function(scn) {
    isTRUE(model_results[[scn]]$success)
  }, logical(1))
  if (!all(success)) {
    stop(sprintf("Unsuccessful required scenarios: %s",
                 paste(scenario_keep[!success], collapse = ", ")))
  }

  base_scn <- scenario_keep[1]
  base_res <- model_results[[base_scn]]
  obs0 <- safe_numeric(base_res$observed)
  if (is.null(obs0)) stop("Base scenario has no observed values.")

  coords0 <- extract_coords(base_res)
  if (is.null(coords0)) {
    if (!allow_observation_bootstrap_fallback) {
      stop("No coordinates found. Spatial-block bootstrap cannot be computed.")
    }
    coords0 <- data.frame(lon = seq_along(obs0), lat = rep(0, length(obs0)))
  }

  fold0 <- base_res$fold_assignment %||% base_res$fold_indices %||% rep(NA_integer_, length(obs0))
  fold0 <- safe_numeric(fold0)

  pred_list <- list()

  for (scn in scenario_keep) {
    sr <- model_results[[scn]]
    obs <- safe_numeric(sr$observed)
    pred <- extract_predictions(sr)

    if (is.null(obs) || is.null(pred)) {
      stop(sprintf("Scenario %s lacks observed or predictions.", scn))
    }
    if (length(obs) != length(obs0) || length(pred) != length(obs0)) {
      stop(sprintf("Scenario %s length mismatch relative to %s.", scn, base_scn))
    }
    if (!isTRUE(all.equal(obs, obs0, tolerance = 1e-8))) {
      stop(sprintf("Observed values do not align between %s and %s.", base_scn, scn))
    }

    pred_list[[scn]] <- pred
  }

  pred_df <- as.data.frame(pred_list, check.names = FALSE)

  df <- cbind(
    data.frame(
      obs = obs0,
      lon = coords0$lon,
      lat = coords0$lat,
      fold = fold0
    ),
    pred_df
  )

  needed <- c("obs", "lon", "lat", scenario_keep)
  ok <- stats::complete.cases(df[, needed, drop = FALSE])
  n_removed <- sum(!ok)
  if (n_removed > 0) {
    message(sprintf("    Removed %d incomplete rows before Table S2 bootstrap.", n_removed))
  }
  df <- df[ok, , drop = FALSE]

  if (nrow(df) < 30) {
    stop(sprintf("Too few complete observations after alignment: n = %d", nrow(df)))
  }

  if (allow_observation_bootstrap_fallback && all(df$lat == 0)) {
    df$block_id <- paste0("obs_", seq_len(nrow(df)))
  } else {
    df$block_id <- make_spatial_blocks(df$lon, df$lat, block_size_km)
  }

  df
}

compute_r2_vector <- function(df, scenario_keep, idx = NULL) {
  if (is.null(idx)) {
    d <- df
  } else {
    d <- df[idx, , drop = FALSE]
  }
  obs <- d$obs
  r2 <- vapply(scenario_keep, function(scn) {
    pooled_r2(obs, d[[scn]])
  }, numeric(1))
  r2
}

marginal_from_r2 <- function(r2, def_row) {
  full <- def_row$Full
  reduced <- def_row$Reduced
  if (!all(c(full, reduced) %in% names(r2))) return(NA_real_)
  r2[[full]] - r2[[reduced]]
}

substitution_from_r2 <- function(r2, def_row) {
  need <- c(def_row$Scenario_C,
            def_row$Scenario_C_plus_B,
            def_row$Scenario_A_plus_C,
            def_row$Scenario_A_plus_B_plus_C)
  if (!all(need %in% names(r2))) return(NA_real_)

  m_b_given_c <- r2[[def_row$Scenario_C_plus_B]] - r2[[def_row$Scenario_C]]
  m_b_given_ac <- r2[[def_row$Scenario_A_plus_B_plus_C]] - r2[[def_row$Scenario_A_plus_C]]
  m_b_given_c - m_b_given_ac
}

summarise_bootstrap <- function(estimate, boot_values, ci_level = 0.95) {
  boot_values <- boot_values[is.finite(boot_values)]
  n_used <- length(boot_values)

  if (n_used < 100) {
    return(data.frame(
      estimate = estimate,
      ci_lo = NA_real_,
      ci_hi = NA_real_,
      p = NA_real_,
      n_boot_used = n_used
    ))
  }

  alpha <- 1 - ci_level
  ci <- stats::quantile(boot_values, c(alpha / 2, 1 - alpha / 2), na.rm = TRUE)

  # Plus-one achieved-significance p-value.
  p_left  <- (sum(boot_values <= 0, na.rm = TRUE) + 1) / (n_used + 1)
  p_right <- (sum(boot_values >= 0, na.rm = TRUE) + 1) / (n_used + 1)
  p <- min(1, 2 * min(p_left, p_right))

  data.frame(
    estimate = estimate,
    ci_lo = unname(ci[1]),
    ci_hi = unname(ci[2]),
    p = p,
    n_boot_used = n_used
  )
}

fmt_pct <- function(x) {
  ifelse(is.na(x), NA_character_, sprintf("%.2f", 100 * x))
}

add_bh <- function(df) {
  if (is.null(df) || nrow(df) == 0) return(df)
  df$p_BH <- NA_real_

  fam <- interaction(
    df$Data_Type, df$Target, df$PFT, df$CV_Method, df$Contrast_Family,
    drop = TRUE, lex.order = TRUE
  )

  df$p_BH <- ave(df$p, fam, FUN = function(x) {
    stats::p.adjust(x, method = "BH")
  })
  df$Sig_BH <- ifelse(is.na(df$p_BH), NA, df$p_BH < 0.05)
  df
}

run_bootstrap_for_result <- function(res, key, scenario_keep,
                                     n_boot = 2000,
                                     block_size_km = 50,
                                     ci_level = 0.95,
                                     seed = 42) {
  data_type <- res$data_type %||% "match_only"
  target <- res$target_var
  pft <- res$pft
  cv_method <- res$cv_method %||% "kNNDM"

  cat(sprintf("\n[Table S2] %s | %s | %s | %s\n", data_type, target, pft, cv_method))

  df <- build_prediction_table(res, scenario_keep)
  n_obs <- nrow(df)
  block_ids <- sort(unique(df$block_id))
  n_blocks <- length(block_ids)

  if (n_blocks < 5 && !allow_observation_bootstrap_fallback) {
    warning(sprintf("Only %d spatial blocks found for %s. CI may be unstable.",
                    n_blocks, key), call. = FALSE)
  }

  block_to_idx <- split(seq_len(nrow(df)), df$block_id)

  r2_actual <- compute_r2_vector(df, scenario_keep)
  if (any(!is.finite(r2_actual))) {
    stop(sprintf("Non-finite actual pooled R2 for key: %s", key))
  }

  # Point estimates
  marginal_est <- vapply(seq_len(nrow(marginal_def)), function(i) {
    marginal_from_r2(r2_actual, marginal_def[i, ])
  }, numeric(1))

  # Section 0 (NEW v2.6): block-removal contrasts use the same math.
  removal_est <- vapply(seq_len(nrow(removal_def)), function(i) {
    marginal_from_r2(r2_actual, removal_def[i, ])
  }, numeric(1))

  substitution_est <- vapply(seq_len(nrow(substitution_def)), function(i) {
    substitution_from_r2(r2_actual, substitution_def[i, ])
  }, numeric(1))

  # Bootstrap
  set.seed(seed)
  boot_marg <- matrix(NA_real_, nrow = n_boot, ncol = nrow(marginal_def))
  boot_rem  <- matrix(NA_real_, nrow = n_boot, ncol = nrow(removal_def))
  boot_sub  <- matrix(NA_real_, nrow = n_boot, ncol = nrow(substitution_def))

  for (b in seq_len(n_boot)) {
    sampled_blocks <- sample(block_ids, size = n_blocks, replace = TRUE)
    idx <- unlist(block_to_idx[sampled_blocks], use.names = FALSE)

    r2_b <- compute_r2_vector(df, scenario_keep, idx = idx)

    if (any(!is.finite(r2_b))) next

    boot_marg[b, ] <- vapply(seq_len(nrow(marginal_def)), function(i) {
      marginal_from_r2(r2_b, marginal_def[i, ])
    }, numeric(1))

    boot_rem[b, ] <- vapply(seq_len(nrow(removal_def)), function(i) {
      marginal_from_r2(r2_b, removal_def[i, ])
    }, numeric(1))

    boot_sub[b, ] <- vapply(seq_len(nrow(substitution_def)), function(i) {
      substitution_from_r2(r2_b, substitution_def[i, ])
    }, numeric(1))

    if (b %% max(1L, floor(n_boot / 10L)) == 0L) {
      cat(sprintf("  bootstrap %d/%d\n", b, n_boot))
    }
  }

  # Summarize marginal contrasts
  marginal_rows <- list()
  for (i in seq_len(nrow(marginal_def))) {
    stat <- summarise_bootstrap(marginal_est[i], boot_marg[, i], ci_level)
    marginal_rows[[i]] <- cbind(
      data.frame(
        Data_Type = data_type,
        Target = target,
        PFT = pft,
        CV_Method = cv_method,
        Model = "RF",
        Table_Panel = marginal_def$Table_Panel[i],
        Contrast_Family = marginal_def$Contrast_Family[i],
        Contrast_Type = "marginal_contribution",
        Contrast_ID = marginal_def$Contrast_ID[i],
        Block_Added = marginal_def$Block_Added[i],
        Context = marginal_def$Context[i],
        Reduced = marginal_def$Reduced[i],
        Full = marginal_def$Full[i],
        Substitution_Pair = NA_character_,
        Present_Block_A = NA_character_,
        Marginal_Block_B = NA_character_,
        Conditioning_Context_C = NA_character_,
        Interpretation = marginal_def$Interpretation[i],
        N_Obs = n_obs,
        N_Blocks = n_blocks,
        Block_Size_km = block_size_km,
        stringsAsFactors = FALSE
      ),
      stat
    )
  }
  marginal_out <- do.call(rbind, marginal_rows)

  # Summarize Section 0 block-removal contrasts (NEW v2.6)
  removal_rows <- list()
  for (i in seq_len(nrow(removal_def))) {
    stat <- summarise_bootstrap(removal_est[i], boot_rem[, i], ci_level)
    removal_rows[[i]] <- cbind(
      data.frame(
        Data_Type = data_type,
        Target = target,
        PFT = pft,
        CV_Method = cv_method,
        Model = "RF",
        Table_Panel = removal_def$Table_Panel[i],
        Contrast_Family = removal_def$Contrast_Family[i],
        Contrast_Type = "block_removal",
        Contrast_ID = removal_def$Contrast_ID[i],
        Block_Removed = removal_def$Block_Removed[i],
        Context_Retained = removal_def$Context_Retained[i],
        Reduced = removal_def$Reduced[i],
        Full = removal_def$Full[i],
        Interpretation = removal_def$Interpretation[i],
        N_Obs = n_obs,
        N_Blocks = n_blocks,
        Block_Size_km = block_size_km,
        stringsAsFactors = FALSE
      ),
      stat
    )
  }
  removal_out <- do.call(rbind, removal_rows)

  # Summarize substitution contrasts
  substitution_rows <- list()
  for (i in seq_len(nrow(substitution_def))) {
    stat <- summarise_bootstrap(substitution_est[i], boot_sub[, i], ci_level)
    substitution_rows[[i]] <- cbind(
      data.frame(
        Data_Type = data_type,
        Target = target,
        PFT = pft,
        CV_Method = cv_method,
        Model = "RF",
        Table_Panel = substitution_def$Table_Panel[i],
        Contrast_Family = substitution_def$Contrast_Family[i],
        Contrast_Type = "substitution_contrast",
        Contrast_ID = substitution_def$Contrast_ID[i],
        Block_Added = NA_character_,
        Context = substitution_def$Conditioning_Context_C[i],
        Reduced = NA_character_,
        Full = NA_character_,
        Substitution_Pair = substitution_def$Substitution_Pair[i],
        Present_Block_A = substitution_def$Present_Block_A[i],
        Marginal_Block_B = substitution_def$Marginal_Block_B[i],
        Conditioning_Context_C = substitution_def$Conditioning_Context_C[i],
        Interpretation = substitution_def$Interpretation[i],
        Scenario_C = substitution_def$Scenario_C[i],
        Scenario_C_plus_B = substitution_def$Scenario_C_plus_B[i],
        Scenario_A_plus_C = substitution_def$Scenario_A_plus_C[i],
        Scenario_A_plus_B_plus_C = substitution_def$Scenario_A_plus_B_plus_C[i],
        N_Obs = n_obs,
        N_Blocks = n_blocks,
        Block_Size_km = block_size_km,
        stringsAsFactors = FALSE
      ),
      stat
    )
  }
  substitution_out <- do.call(rbind, substitution_rows)

  # Scenario R2 summary for traceability
  scenario_r2 <- data.frame(
    Data_Type = data_type,
    Target = target,
    PFT = pft,
    CV_Method = cv_method,
    Model = "RF",
    Scenario = names(r2_actual),
    R2_pooled = as.numeric(r2_actual),
    N_Obs = n_obs,
    N_Blocks = n_blocks,
    stringsAsFactors = FALSE
  )

  list(
    marginal = marginal_out,
    removal = removal_out,
    substitution = substitution_out,
    scenario_r2 = scenario_r2
  )
}

# -- Main ----------------------------------------------------------------------

rds_path <- file.path(output_dir, "all_results.rds")
if (!file.exists(rds_path)) {
  stop("Could not find all_results.rds in: ", output_dir,
       "\nThis script needs the full results object with predictions, observed values, and coordinates.")
}

cat("====================================================================\n")
cat("Table S2 spatial-block marginal contribution and substitution script\n")
cat("====================================================================\n")
cat(sprintf("Input all_results: %s\n", rds_path))
cat(sprintf("Output directory: %s\n", output_subdir))
cat(sprintf("Bootstrap replicates: %d\n", n_boot))
cat(sprintf("Spatial block size: %.1f km\n", block_size_km))
cat(sprintf("CI level: %.3f\n", ci_level))
cat(sprintf("Seed: %d\n", seed))

all_results <- readRDS(rds_path)
# Validate before any paired bootstrap or publication output.
.fem_sources <- unlist(lapply(sys.frames(),function(e)e$ofile),use.names=FALSE)
.fem_script <- if(length(.fem_sources))tail(.fem_sources,1L) else sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1])
source(file.path(dirname(normalizePath(.fem_script)),"..","..","R","oof_contracts.R"))
rm(.fem_sources,.fem_script)
all_results <- fem_use_raw_rf_oof(all_results)
cat("Prediction type: RF_raw (uncalibrated outer-test OOF); QM excluded.\n")

if (is.null(all_results) || length(all_results) == 0) {
  stop("all_results.rds is empty or invalid.")
}

all_marginal <- list()
all_removal <- list()
all_substitution <- list()
all_scenario_r2 <- list()
skipped <- list()

key_counter <- 0L

for (key in names(all_results)) {
  key_counter <- key_counter + 1L
  res <- all_results[[key]]

  if (is.null(res$model_results)) {
    skipped[[length(skipped) + 1]] <- data.frame(Key = key, Reason = "No model_results")
    next
  }

  out <- tryCatch({
    run_bootstrap_for_result(
      res = res,
      key = key,
      scenario_keep = scenario_keep,
      n_boot = n_boot,
      block_size_km = block_size_km,
      ci_level = ci_level,
      seed = seed + 10000L * key_counter
    )
  }, error = function(e) {
    skipped[[length(skipped) + 1]] <<- data.frame(Key = key, Reason = e$message)
    NULL
  })

  if (is.null(out)) next
  all_marginal[[length(all_marginal) + 1]] <- out$marginal
  all_removal[[length(all_removal) + 1]] <- out$removal
  all_substitution[[length(all_substitution) + 1]] <- out$substitution
  all_scenario_r2[[length(all_scenario_r2) + 1]] <- out$scenario_r2
}

if (length(skipped)) {
  write.csv(do.call(rbind,skipped),file.path(output_subdir,"Table_S1_skipped_results.csv"),row.names=FALSE)
  stop("Incomplete result strata; refusing a partial publication table. See Table_S1_skipped_results.csv")
}

if (length(all_marginal) == 0 || length(all_substitution) == 0 || length(all_removal) == 0) {
  if (length(skipped) > 0) {
    skipped_df <- do.call(rbind, skipped)
    write.csv(skipped_df, file.path(output_subdir, "Table_S1_skipped_results.csv"),
              row.names = FALSE, fileEncoding = "UTF-8")
  }
  stop("No Table S1 rows generated. See Table_S1_skipped_results.csv if present.")
}

tab_s1_removal <- do.call(rbind, all_removal)
tab_s1_marginal <- do.call(rbind, all_marginal)
tab_s1_substitution <- do.call(rbind, all_substitution)
scenario_r2 <- do.call(rbind, all_scenario_r2)

# Apply BH adjustment within each Contrast_Family
tab_s1_removal <- add_bh(tab_s1_removal)
tab_s1_marginal <- add_bh(tab_s1_marginal)
tab_s1_substitution <- add_bh(tab_s1_substitution)

# Ensure columns missing from one panel exist before rbind for combined output.
all_cols <- unique(c(names(tab_s1_removal), names(tab_s1_marginal), names(tab_s1_substitution)))
for (cc in setdiff(all_cols, names(tab_s1_removal))) tab_s1_removal[[cc]] <- NA
for (cc in setdiff(all_cols, names(tab_s1_marginal))) tab_s1_marginal[[cc]] <- NA
for (cc in setdiff(all_cols, names(tab_s1_substitution))) tab_s1_substitution[[cc]] <- NA
combined <- rbind(
  tab_s1_removal[, all_cols, drop = FALSE],
  tab_s1_marginal[, all_cols, drop = FALSE],
  tab_s1_substitution[, all_cols, drop = FALSE]
)

# Add percent columns for easy manuscript/SI formatting.
for (df_name in c("tab_s1_removal", "tab_s1_marginal", "tab_s1_substitution", "combined")) {
  df <- get(df_name)
  df$estimate_pct <- fmt_pct(df$estimate)
  df$ci_lo_pct <- fmt_pct(df$ci_lo)
  df$ci_hi_pct <- fmt_pct(df$ci_hi)
  df$CI_pct <- ifelse(is.na(df$ci_lo_pct) | is.na(df$ci_hi_pct),
                      NA_character_,
                      paste0("[", df$ci_lo_pct, ", ", df$ci_hi_pct, "]"))
  assign(df_name, df)
}

# Write outputs (v2.6: 4 CSVs).
out_s1_0 <- file.path(output_subdir, "Table_S1_consolidated_block_removal.csv")
out_s1_a <- file.path(output_subdir, "Table_S1_consolidated_marginal.csv")
out_s1_b <- file.path(output_subdir, "Table_S1_consolidated_substitution.csv")
out_combined <- file.path(output_subdir, "Table_S1_consolidated_all.csv")
out_r2 <- file.path(output_subdir, "Table_S1_scenario_pooled_R2_trace.csv")

tab_s1_removal$Prediction_Type <- "RF_raw"
tab_s1_marginal$Prediction_Type <- "RF_raw"
tab_s1_substitution$Prediction_Type <- "RF_raw"
combined$Prediction_Type <- "RF_raw"
scenario_r2$Prediction_Type <- "RF_raw"
write.csv(tab_s1_removal, out_s1_0, row.names = FALSE, fileEncoding = "UTF-8")
write.csv(tab_s1_marginal, out_s1_a, row.names = FALSE, fileEncoding = "UTF-8")
write.csv(tab_s1_substitution, out_s1_b, row.names = FALSE, fileEncoding = "UTF-8")
write.csv(combined, out_combined, row.names = FALSE, fileEncoding = "UTF-8")
write.csv(scenario_r2, out_r2, row.names = FALSE, fileEncoding = "UTF-8")

if (length(skipped) > 0) {
  skipped_df <- do.call(rbind, skipped)
  write.csv(skipped_df, file.path(output_subdir, "Table_S1_skipped_results.csv"),
            row.names = FALSE, fileEncoding = "UTF-8")
}

# Methods / caption note.
methods_note <- c(
  "Table S1 (consolidated) methods note (v2.6)",
  "Prediction type: RF_raw. Scenario R2, source contrasts and every bootstrap replicate use uncalibrated outer-test RF OOF predictions; RF_QM is excluded.",
  "",
  paste0("This table merges three contrast families into a single SI table: (Section 0) block-removal contrasts, (Section A) context-dependent marginal contributions, and (Section B) block-level substitution contrasts. All contrasts share the same paired observation-level pooled OOF R2 framework. Marginal contributions M(Block | Context) = R2(Context + Block) - R2(Context). Substitution contrasts S(A, B | C) = M(B | C) - M(B | A + C). Block-removal contrasts are mathematically of the same form (R2_full - R2_reduced) but framed as removing a block from a larger model."),
  "",
  paste0("Uncertainty was estimated with a paired spatial-block percentile bootstrap using B = ", n_boot, " replicates and ", block_size_km, " km blocks. Within each bootstrap replicate, spatial blocks were resampled with replacement and all scenarios were evaluated on the same resampled observations. Confidence intervals are percentile intervals at the ", 100 * ci_level, "% level. Two-sided p-values were computed using an achieved-significance approach with a plus-one correction and were adjusted using the Benjamini-Hochberg procedure within Data_Type x Target x PFT x CV_Method x Contrast_Family. The four contrast families (Block_removal_8, RS_marginal_9, Env_RS_complementarity_2, RS_substitution_3) were adjusted independently."),
  "",
  "Primary interpretation target for Discussion 4.3 §2:",
  "  S_Prod_Temp_given_Het tests whether the marginal value of Temp is reduced when Prod is present under the Het-conditioned context. This directly addresses the apparent mismatch between high individual importance of temporal predictors and non-significant Temp removal in the full RS model."
)
writeLines(methods_note, con = file.path(output_subdir, "Table_S1_methods_note.txt"), useBytes = TRUE)

cat("\n====================================================================\n")
cat("Saved Table S1 (consolidated) outputs (v2.6)\n")
cat("====================================================================\n")
cat(sprintf("Section 0 block-removal:    %s\n", out_s1_0))
cat(sprintf("Section A marginal contrib: %s\n", out_s1_a))
cat(sprintf("Section B substitution:     %s\n", out_s1_b))
cat(sprintf("Combined (88 rows):         %s\n", out_combined))
cat(sprintf("Scenario pooled R2 trace:   %s\n", out_r2))
cat(sprintf("Methods note:               %s\n", file.path(output_subdir, "Table_S1_methods_note.txt")))

if (length(skipped) > 0) {
  cat(sprintf("Skipped result keys: %d. See Table_S1_skipped_results.csv\n", length(skipped)))
}

# Print key lines for Discussion 4.3.
cat("\n\n=== Section 0 block-removal contrasts (replaces previous Table S1) ===\n")
print(tab_s1_removal[, c("Target", "PFT", "CV_Method", "Contrast_ID",
                          "Block_Removed", "Context_Retained",
                          "estimate_pct", "CI_pct", "p", "p_BH", "Sig_BH",
                          "N_Obs", "N_Blocks")],
      row.names = FALSE)

cat("\n\n=== Key substitution target for Discussion 4.3 ===\n")
key_sub <- subset(tab_s1_substitution, Contrast_ID == "S_Prod_Temp_given_Het")
if (nrow(key_sub) > 0) {
  print(key_sub[, c("Target", "PFT", "CV_Method", "Contrast_ID",
                    "estimate_pct", "CI_pct", "p", "p_BH", "Sig_BH",
                    "N_Obs", "N_Blocks")],
        row.names = FALSE)
}

cat("\n\n=== Temp marginal contribution by context ===\n")
key_temp <- subset(tab_s1_marginal, Block_Added == "Temp" & Contrast_Family == "RS_marginal_9")
if (nrow(key_temp) > 0) {
  print(key_temp[, c("Target", "PFT", "CV_Method", "Block_Added", "Context",
                     "estimate_pct", "CI_pct", "p", "p_BH", "Sig_BH",
                     "N_Obs", "N_Blocks")],
        row.names = FALSE)
}

cat("\nDone.\n")
