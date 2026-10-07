simulate_ar1_residual_gaussian <- function(resid, rho, n = length(resid)) {
  sig <- sqrt(sum(resid^2, na.rm = TRUE) / max(sum(is.finite(resid)) - 3L, 1L))
  simulate_ar1_with_marginal_sd(n, rho, sig)
}

simulate_ar1_residual_empirical <- function(resid, rho, n = length(resid)) {
  resid <- resid[is.finite(resid)]; m <- length(resid)
  if (m < 4L) stop("Too few residuals for empirical AR bootstrap.", call. = FALSE)
  innov <- resid[-1L] - rho * resid[-m]; innov <- innov - mean(innov)
  e <- numeric(n); e[[1L]] <- sample(resid, 1L)
  if (n > 1L) { z <- sample(innov, n - 1L, replace = TRUE); for (tt in 2:n) e[[tt]] <- rho * e[[tt - 1L]] + z[[tt - 1L]] }
  e
}

estimate_lag1_single <- function(resid, max_abs = 0.95) {
  resid <- resid[is.finite(resid)]; if (length(resid) < 3L) return(0)
  a <- resid[-length(resid)]; b <- resid[-1L]; den <- sqrt(sum(a^2) * sum(b^2)); r <- if (den > 0) sum(a*b)/den else 0
  max(-max_abs, min(max_abs, r))
}

bootstrap_regional_trajectory <- function(y, years, fit, rho, null_sorted, cfg, key, response) {
  tcfg <- cfg$bootstrap$trajectory; B <- as.integer(tcfg$B)
  seed <- seed_from_key(cfg$rng$master_seed, "trajectory_boot", key); set_rng(cfg$rng$kind, seed)
  gen <- if (tolower(tcfg$innovation) == "empirical_ar1") simulate_ar1_residual_empirical else simulate_ar1_residual_gaussian
  alpha <- as.numeric(cfg$analysis$primary$alpha); tau <- rep(NA_real_, B); rec <- logical(B); det <- logical(B)
  slope_tol <- as.numeric(cfg$analysis$slope_tol[[response]] %||% 0); m <- as.integer(cfg$analysis$primary$min_segment)
  for (b in seq_len(B)) {
    ys <- fit$fitted + gen(fit$resid, rho, length(y))
    z <- fit_break_core(ys, years, m, m, slope_tol, cfg$analysis$tie_rule)
    if (isTRUE(z$testable)) {
      pb <- mc_pvalue(z$supF, null_sorted); det[[b]] <- is.finite(pb) && pb <= alpha
      rec[[b]] <- det[[b]] && isTRUE(z$direction_ok); if (rec[[b]]) tau[[b]] <- z$break_year
    }
  }
  tt <- tau[is.finite(tau)]; minrep <- as.integer(tcfg$min_valid_replicates_for_ci)
  q <- if (length(tt) >= minrep) safe_quantile_discrete(tt, c(.025,.5,.975)) else c(NA,NA,NA)
  data.frame(bootstrap_detection_frequency = mean(det), bootstrap_recovery_frequency = mean(rec), n_recovery_boot = sum(rec),
             bootstrap_break_q025 = q[[1L]], bootstrap_break_median = q[[2L]], bootstrap_break_q975 = q[[3L]],
             bootstrap_break_sd = if(length(tt)>1) stats::sd(tt) else NA_real_,
             bootstrap_break_CI_width = if (all(is.finite(q[c(1,3)]))) q[[3L]]-q[[1L]] else NA_real_,
             ci_status = if(length(tt)>=minrep) "ok" else "unstable", stringsAsFactors=FALSE)
}
