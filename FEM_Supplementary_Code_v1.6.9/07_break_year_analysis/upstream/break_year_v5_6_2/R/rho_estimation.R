linear_residual_matrix <- function(Y, years) {
  Y <- as.matrix(Y); t <- years - mean(years); M <- residual_maker(cbind(1, t)); Y %*% M
}

hinge_residual_matrix <- function(Y, years, min_segment = 3L) {
  Y <- as.matrix(Y); n <- ncol(Y); t <- years - years[[1L]]; cand <- candidate_indices(n, min_segment, min_segment)
  best <- rep(Inf, nrow(Y)); residual <- matrix(NA_real_, nrow(Y), n)
  for (idx in cand) {
    X <- cbind(1, t, pmax(0, years - years[[idx]])); M <- residual_maker(X); R <- Y %*% M; sse <- rowSums(R^2)
    take <- sse < best
    if (any(take)) { residual[take, ] <- R[take, , drop = FALSE]; best[take] <- sse[take] }
  }
  residual
}

lag1_correlation_rows <- function(R) {
  R <- as.matrix(R); left <- R[, -ncol(R), drop = FALSE]; right <- R[, -1L, drop = FALSE]
  den <- sqrt(rowSums(left^2) * rowSums(right^2)); z <- rowSums(left * right) / den; z[is.finite(z)]
}

estimate_residual_rho_values <- function(Y, years, residual_model = c("linear", "hinge"), hinge_min_segment = 3L) {
  residual_model <- match.arg(residual_model); Y <- as.matrix(Y)
  keep <- rowSums(is.finite(Y)) == ncol(Y); Y <- Y[keep, , drop = FALSE]
  if (!nrow(Y)) return(numeric())
  R <- if (residual_model == "linear") linear_residual_matrix(Y, years) else hinge_residual_matrix(Y, years, hinge_min_segment)
  lag1_correlation_rows(R)
}

rho_bias_mapping <- function(years, grid, B, residual_model, seed, hinge_min_segment = 3L, rng_kind = "L'Ecuyer-CMRG") {
  set_rng(rng_kind, seed); out <- vector("list", length(grid))
  for (i in seq_along(grid)) {
    Y <- simulate_ar1_matrix(B, length(years), grid[[i]], 1)
    rr <- estimate_residual_rho_values(Y, years, residual_model, hinge_min_segment)
    out[[i]] <- data.frame(true_rho = grid[[i]], median_estimated_rho = stats::median(rr),
                           q10_estimated_rho = safe_quantile_discrete(rr, 0.10), q90_estimated_rho = safe_quantile_discrete(rr, 0.90),
                           n_valid = length(rr), stringsAsFactors = FALSE)
  }
  tab <- do.call(rbind, out)
  iso <- stats::isoreg(tab$true_rho, tab$median_estimated_rho)
  tab$median_isotonic <- iso$yf
  tab
}


rho_bias_mapping_cached <- function(cfg, years, grid, B, residual_model, seed, hinge_min_segment = 3L) {
  cache_dir <- ensure_dir(file.path(cfg$output$root, "_cache", "rho_bias"))
  r_version <- paste(R.version$major, R.version$minor, sep = ".")
  key <- substr(hash_text(c("v5.6-rho-bias", cfg$.code_hash %||% "no-code-hash", paste(years, collapse = ","),
                          paste(grid, collapse = ","), B, residual_model, hinge_min_segment,
                          cfg$rng$kind, seed, r_version)), 1L, 16L)
  path <- file.path(cache_dir, paste0("rho_bias_", residual_model, "_", key, ".rds"))
  if (file.exists(path)) { z <- safe_read_rds_cache(path); if (!is.null(z)) return(z) }
  z <- rho_bias_mapping(years, grid, B, residual_model, seed, hinge_min_segment, cfg$rng$kind)
  atomic_save_rds(z, path); z
}

bias_correct_rho <- function(raw, mapping, max_abs) {
  if (!is.finite(raw) || is.null(mapping) || !nrow(mapping)) return(max(-max_abs, min(max_abs, raw)))
  z <- mapping[is.finite(mapping$true_rho) & is.finite(mapping$median_isotonic), c("true_rho", "median_isotonic")]
  # Collapse isotonic plateaus before inverse interpolation.
  z <- stats::aggregate(true_rho ~ median_isotonic, data = z, FUN = median); z <- z[order(z$median_isotonic), ]
  corrected <- if (nrow(z) == 1L) z$true_rho[[1L]] else stats::approx(z$median_isotonic, z$true_rho, xout = raw, rule = 2, ties = mean)$y
  max(-max_abs, min(max_abs, corrected))
}

sample_complete_series <- function(stack, n, oversample_factor, seed, rng_kind = "L'Ecuyer-CMRG") {
  assert_packages("terra"); set_rng(rng_kind, seed)
  target <- as.integer(n); draw <- min(terra::ncell(stack), max(target, target * as.integer(oversample_factor)))
  cells <- sample.int(terra::ncell(stack), draw, replace = FALSE)
  v <- as.matrix(terra::extract(stack, cells, raw = TRUE))
  v <- v[rowSums(is.finite(v)) == ncol(v), , drop = FALSE]
  if (nrow(v) > target) v <- v[seq_len(target), , drop = FALSE]
  if (nrow(v) < min(100L, target)) stop("Too few complete sampled series for rho estimation.", call. = FALSE)
  if (nrow(v) < target) warning("rho estimation obtained ", nrow(v), " complete series versus requested ", target, "; consider increasing oversample_factor.", call. = FALSE)
  v
}

rho_response_task <- function(task, state) {
  nm <- task$response; rcfg <- state$rcfg
  st <- load_response_stack(state$inventory$files[[nm]], state$cfg)
  Y <- sample_complete_series(st, rcfg$sample_n, rcfg$oversample_factor,
                              seed_from_key(state$cfg$rng$master_seed, "rho_sample", nm), state$cfg$rng$kind)
  primary_row <- NULL; diag_rows <- list(); dist_rows <- list()
  for (model in state$diagnostic_models) {
    vals <- estimate_residual_rho_values(Y, state$years, model, state$hinge_min_segment)
    raw <- stats::median(vals)
    corrected <- if (isTRUE(rcfg$bias_correction$enabled))
      bias_correct_rho(raw, state$mappings[[model]], as.numeric(rcfg$max_abs_rho)) else raw
    row <- data.frame(response = nm, residual_model = model, n_sampled_complete = nrow(Y), n_valid_rho = length(vals),
                      residual_rho_raw_median = raw, residual_rho_raw_q25 = safe_quantile_discrete(vals, 0.25),
                      residual_rho_raw_q75 = safe_quantile_discrete(vals, 0.75), bias_corrected_rho = corrected,
                      correction_hit_bound = abs(corrected) >= as.numeric(rcfg$max_abs_rho) - 1e-12,
                      is_primary = model == state$primary_model, input_fingerprint = state$inventory$fingerprint,
                      stringsAsFactors = FALSE)
    diag_rows[[model]] <- row
    dist_rows[[model]] <- data.frame(response = nm, residual_model = model, rho = vals)
    if (model == state$primary_model) primary_row <- transform(row, rho_for_calibration = corrected,
      correction_method = if (isTRUE(rcfg$bias_correction$enabled)) "isotonic simulation-inversion of median residual lag-1 correlation" else "none")
  }
  list(primary = primary_row, diagnostic = do.call(rbind, diag_rows), distribution = do.call(rbind, dist_rows))
}

estimate_rho_all <- function(cfg, inventory, ctx) {
  yrs <- as.numeric(unlist(cfg$data$years)); rcfg <- cfg$rho_estimation
  primary_model <- tolower(rcfg$residual_model)
  if (!primary_model %in% c("linear", "hinge")) stop("rho_estimation.residual_model must be linear or hinge.", call. = FALSE)
  diagnostic_models <- unique(tolower(c(primary_model, unlist(rcfg$diagnostic_models %||% c("linear", "hinge")))))
  diagnostic_models <- diagnostic_models[diagnostic_models %in% c("linear", "hinge")]
  grid <- seq(as.numeric(rcfg$bias_correction$grid_min), as.numeric(rcfg$bias_correction$grid_max), by = as.numeric(rcfg$bias_correction$grid_step))
  hms <- as.integer(rcfg$hinge_min_segment %||% 3L)

  mappings <- list()
  for (model in diagnostic_models) {
    mp <- if (isTRUE(rcfg$bias_correction$enabled)) rho_bias_mapping_cached(cfg, yrs, grid, as.integer(rcfg$bias_correction$simulation_B), model,
           seed_from_key(cfg$rng$master_seed, "rho_bias", model), hms) else NULL
    mappings[[model]] <- mp
    if (!is.null(mp)) write_csv(mp, file.path(ctx$root, "rho", paste0("rho_bias_mapping_", model, ".csv")))
  }
  if (!is.null(mappings[[primary_model]])) write_csv(mappings[[primary_model]], file.path(ctx$root, "rho", "rho_bias_mapping.csv"))

  tasks <- lapply(names(inventory$files), function(nm) list(response = nm))
  state <- list(cfg = cfg, inventory = inventory, rcfg = rcfg, years = yrs, primary_model = primary_model,
                diagnostic_models = diagnostic_models, mappings = mappings, hinge_min_segment = hms)
  pieces <- if (parallel_stage_enabled(cfg, "response") && length(tasks) > 1L) {
    psock_task_lapply(cfg, tasks, "rho_response_task", state, stage = "response", log_file = ctx$log)
  } else lapply(tasks, rho_response_task, state = state)

  primary <- do.call(rbind, lapply(pieces, function(z) z$primary))
  diagnostic <- do.call(rbind, lapply(pieces, function(z) z$diagnostic))
  dist <- do.call(rbind, lapply(pieces, function(z) z$distribution))
  if (any(primary$correction_hit_bound, na.rm=TRUE)) warning("At least one primary corrected rho hit rho_estimation.max_abs_rho; expand the bias grid / inspect rho sensitivity before publication.", call.=FALSE)
  write_csv(primary, file.path(ctx$root, "rho", "rho_estimates.csv")); write_csv(diagnostic, file.path(ctx$root, "rho", "rho_diagnostics.csv")); write_csv(dist, file.path(ctx$root, "rho", "rho_sample_distribution.csv"))
  # Large disagreement is not automatically resolved; it is surfaced for sensitivity interpretation.
  wide <- reshape(diagnostic[, c("response", "residual_model", "bias_corrected_rho")], idvar = "response", timevar = "residual_model", direction = "wide")
  if (all(c("bias_corrected_rho.linear", "bias_corrected_rho.hinge") %in% names(wide))) {
    wide$linear_minus_hinge <- wide$bias_corrected_rho.linear - wide$bias_corrected_rho.hinge
    write_csv(wide, file.path(ctx$root, "rho", "rho_model_comparison.csv"))
    if (any(abs(wide$linear_minus_hinge) > as.numeric(rcfg$diagnostic_warning_difference %||% 0.10), na.rm = TRUE))
      warning("Linear- and hinge-residual corrected rho differ materially for at least one response; inspect rho_model_comparison.csv and rho sensitivity results.", call. = FALSE)
  }
  primary
}

rho_list_from_table <- function(tab) stats::setNames(as.list(tab$rho_for_calibration), tab$response)
