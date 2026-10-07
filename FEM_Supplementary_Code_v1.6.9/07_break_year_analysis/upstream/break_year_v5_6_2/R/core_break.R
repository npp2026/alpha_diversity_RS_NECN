candidate_indices <- function(n, min_left, min_right) {
  lo <- as.integer(min_left); hi <- n - as.integer(min_right)
  if (lo > hi) integer() else seq.int(lo, hi)
}

fit_ols_sse <- function(X, y) {
  q <- qr(X)
  if (q$rank < ncol(X)) return(list(ok = FALSE, coef = rep(NA_real_, ncol(X)), fitted = rep(NA_real_, length(y)), resid = rep(NA_real_, length(y)), sse = NA_real_))
  b <- qr.coef(q, y); fitted <- as.numeric(X %*% b); resid <- y - fitted
  list(ok = TRUE, coef = b, fitted = fitted, resid = resid, sse = sum(resid^2))
}

fit_hinge_at_tau <- function(y, years, tau) {
  t <- as.numeric(years - years[[1L]]); z <- pmax(0, as.numeric(years - tau))
  fit <- fit_ols_sse(cbind(intercept = 1, time = t, hinge = z), y)
  if (!fit$ok) return(c(fit, list(tau = tau, slope_before = NA_real_, slope_after = NA_real_)))
  fit$tau <- tau; fit$slope_before <- unname(fit$coef[[2L]]); fit$slope_after <- unname(fit$coef[[2L]] + fit$coef[[3L]])
  fit
}

choose_tie <- function(indices, candidate_years, rule = c("earliest", "latest", "center")) {
  rule <- match.arg(rule)
  if (length(indices) == 1L) return(indices)
  if (rule == "earliest") return(indices[[which.min(candidate_years[indices])]])
  if (rule == "latest") return(indices[[which.max(candidate_years[indices])]])
  center <- mean(range(candidate_years)); z <- indices[abs(candidate_years[indices] - center) == min(abs(candidate_years[indices] - center))]
  z[[which.min(candidate_years[z])]]
}

sse_tolerance <- function(y) {
  # SSE has squared units. Use a relative tolerance on centered variation plus
  # a squared floating-point roundoff floor; do not use an O(eps) absolute SSE
  # floor, which would misclassify small-magnitude ecological indices as exact.
  scale <- sum((y - mean(y))^2)
  amp <- max(1, max(abs(y), na.rm = TRUE))
  roundoff <- length(y) * (.Machine$double.eps * amp)^2 * 1000
  max(100 * .Machine$double.eps * scale, roundoff, .Machine$double.xmin)
}

fit_break_core <- function(y, years, min_left, min_right, slope_tol = 0, tie_rule = "earliest") {
  y <- as.numeric(y); years <- as.numeric(years)
  finite <- is.finite(y) & is.finite(years); y <- y[finite]; years <- years[finite]
  ord <- order(years); y <- y[ord]; years <- years[ord]
  n <- length(y)
  out <- list(n_obs = n, testable = FALSE, fit_status = "insufficient_data", linear_sse = NA_real_,
              break_year = NA_real_, break_index = NA_integer_, hinge_sse = NA_real_, supF = NA_real_,
              slope_before = NA_real_, slope_after = NA_real_, direction_ok = FALSE, tie_count = 0L,
              fitted = rep(NA_real_, n), resid = rep(NA_real_, n))
  if (anyDuplicated(years)) { out$fit_status <- "duplicate_years"; return(out) }
  cand <- candidate_indices(n, min_left, min_right)
  if (!length(cand) || n <= 3L) return(out)

  t <- years - years[[1L]]
  lin <- fit_ols_sse(cbind(1, t), y)
  if (!lin$ok || !is.finite(lin$sse)) { out$fit_status <- "linear_singular"; return(out) }
  out$linear_sse <- lin$sse

  fits <- vector("list", length(cand)); sses <- rep(Inf, length(cand))
  for (j in seq_along(cand)) {
    fits[[j]] <- fit_hinge_at_tau(y, years, years[cand[[j]]])
    if (isTRUE(fits[[j]]$ok) && is.finite(fits[[j]]$sse)) sses[[j]] <- fits[[j]]$sse
  }
  if (!any(is.finite(sses))) { out$fit_status <- "hinge_singular"; return(out) }

  best <- min(sses, na.rm = TRUE)
  tie_tol <- max(sse_tolerance(y), 100 * .Machine$double.eps * abs(best))
  ties <- which(is.finite(sses) & abs(sses - best) <= tie_tol)
  pick <- choose_tie(ties, years[cand], tie_rule)
  hf <- fits[[pick]]; numer <- lin$sse - hf$sse; denom <- hf$sse; eps_sse <- sse_tolerance(y)
  if (numer < 0 && abs(numer) <= eps_sse) numer <- 0

  if (denom <= eps_sse) {
    if (numer <= eps_sse) { Fstat <- 0; status <- "perfect_linear" }
    else { Fstat <- Inf; status <- "exact_hinge" }
  } else {
    Fstat <- (n - 3L) * numer / denom
    if (Fstat < 0 && abs(Fstat) <= 100 * .Machine$double.eps) Fstat <- 0
    if (!is.finite(Fstat) || Fstat < 0) { out$fit_status <- "numerical_failure"; return(out) }
    status <- "normal"
  }

  out$testable <- TRUE; out$fit_status <- status; out$break_year <- hf$tau; out$break_index <- cand[[pick]]
  out$hinge_sse <- hf$sse; out$supF <- Fstat; out$slope_before <- hf$slope_before; out$slope_after <- hf$slope_after
  tol <- abs(as.numeric(slope_tol)); out$direction_ok <- is.finite(hf$slope_before) && is.finite(hf$slope_after) && hf$slope_before < -tol && hf$slope_after > tol
  out$tie_count <- length(ties); out$fitted <- hf$fitted; out$resid <- hf$resid
  out
}

fit_break_vector <- function(y, years, min_left, min_right, slope_tol, tie_rule, null_sorted = NULL) {
  z <- fit_break_core(y, years, min_left, min_right, slope_tol, tie_rule)
  if (!z$testable) return(c(n_obs = z$n_obs, supF = NA, break_year = NA, slope_before = NA, slope_after = NA,
                            direction_ok = 0, mc_exceed = NA, p_mc = NA, fit_status = 0, tie_count = z$tie_count))
  k <- if (is.null(null_sorted)) NA_integer_ else mc_exceed_ge(z$supF, null_sorted)
  p <- if (is.null(null_sorted)) NA_real_ else mc_pvalue_from_exceed(k, length(null_sorted))
  status_code <- switch(z$fit_status, normal = 1, perfect_linear = 2, exact_hinge = 3, 0)
  c(n_obs = z$n_obs, supF = z$supF, break_year = z$break_year, slope_before = z$slope_before,
    slope_after = z$slope_after, direction_ok = as.numeric(z$direction_ok), mc_exceed = k,
    p_mc = p, fit_status = status_code, tie_count = z$tie_count)
}
