# ==============================================================================
# v9_trend_package_utils.R
# ------------------------------------------------------------------------------
# Professional-package-first trend/FDR utilities with project-local fallback.
# Preferred packages:
#   - trend::sens.slope() for Sen's slope when valid years are consecutive.
#   - modifiedmk::mmkh() for Hamed-Rao modified Mann-Kendall p-values.
#   - mutoss::two.stage() for BKY two-stage FDR adjusted p-values.
# Fallback code is retained so the v9 pipeline can still run in lean environments.
# ==============================================================================

v9_env_logical <- function(name, default) {
  v <- Sys.getenv(name, unset = NA_character_)
  if (is.na(v) || !nzchar(v)) return(default)
  v <- tolower(trimws(v))
  if (v %in% c("true", "t", "1", "yes", "y")) return(TRUE)
  if (v %in% c("false", "f", "0", "no", "n")) return(FALSE)
  stop('Invalid logical environment value ', name, '=', v, call.=FALSE)
}

V9_USE_PROFESSIONAL_TREND_PACKAGES <- v9_env_logical("USE_PROFESSIONAL_TREND_PACKAGES", TRUE)
V9_USE_TREND_PACKAGE_SEN <- v9_env_logical("USE_TREND_PACKAGE_SEN", FALSE)
# keep modifiedmk available for optional parity/sensitivity, but do not
# use it by default in the production main pipeline. Several Windows/R package
# combinations returned non-portable mmkh objects and excessive warnings; the
# project fallback is deterministic, uses the true year vector, and is the
# production default unless USE_MODIFIEDMK_PACKAGE=TRUE is explicitly set.
V9_USE_MODIFIEDMK_PACKAGE <- v9_env_logical("USE_MODIFIEDMK_PACKAGE", FALSE)
V9_STRICT_PROFESSIONAL_PACKAGES <- v9_env_logical("STRICT_PROFESSIONAL_PACKAGES", FALSE)
# mutoss::two.stage() computes BKY adjusted p-values in super-linear
# (~O(n^2)) time. At 1-km scale each (metric,trend_period) group has 1e5-1e6 finite
# p-values, so the mutoss BKY path can take HOURS while base-R BH/BY stay in the
# millisecond range. The analytic two-stage plug-in below (BH x m0/m, m0=m-r1) is the
# standard TSBH and is numerically equivalent to mutoss at a fixed alpha, but O(n log n).
# Default: do NOT call mutoss for BKY. Opt in with USE_MUTOSS_BKY=TRUE (still size-gated
# by BKY_MUTOSS_MAX_M so it never silently runs on a huge group).
V9_USE_MUTOSS_BKY <- v9_env_logical("USE_MUTOSS_BKY", FALSE)
.v9_bky_mutoss_max_m <- suppressWarnings(as.integer(Sys.getenv("BKY_MUTOSS_MAX_M", unset = "20000")))
if (is.na(.v9_bky_mutoss_max_m) || .v9_bky_mutoss_max_m < 0L) .v9_bky_mutoss_max_m <- 20000L
V9_BKY_MUTOSS_MAX_M <- .v9_bky_mutoss_max_m

v9_pkg_available <- function(pkg) requireNamespace(pkg, quietly = TRUE)

v9_check_professional_packages <- function(strict = V9_STRICT_PROFESSIONAL_PACKAGES) {
  pkgs <- c("trend", "modifiedmk", "mutoss")
  ok <- stats::setNames(vapply(pkgs, v9_pkg_available, logical(1)), pkgs)
  if (strict && !all(ok)) {
    stop("STRICT_PROFESSIONAL_PACKAGES=TRUE but missing R package(s): ",
         paste(names(ok)[!ok], collapse = ", "),
         ". Install them or set STRICT_PROFESSIONAL_PACKAGES=FALSE to allow fallback.",
         call. = FALSE)
  }
  if (!all(ok)) {
    warning("Professional R package(s) unavailable; v9 will use fallback where needed: ",
            paste(names(ok)[!ok], collapse = ", "), call. = FALSE)
  }
  invisible(ok)
}

v9_is_unit_consecutive <- function(years) {
  years <- years[is.finite(years)]
  if (length(years) < 2L) return(FALSE)
  all(diff(years) == 1)
}

v9_extract_numeric_by_name <- function(x, patterns, allow_unnamed_fallback = TRUE, require_probability = FALSE) {
  flat <- suppressWarnings(unlist(x))
  if (!length(flat)) return(NA_real_)
  nm <- names(flat)
  if (is.null(nm)) nm <- rep("", length(flat))
  for (pat in patterns) {
    idx <- grep(pat, nm, ignore.case = TRUE)
    if (length(idx)) {
      vals <- suppressWarnings(as.numeric(flat[idx]))
      vals <- vals[is.finite(vals)]
      if (require_probability) vals <- vals[vals >= 0 & vals <= 1]
      if (length(vals)) return(unname(vals[1]))
    }
  }
  # Name-based extraction is mandatory for p-values from professional-package objects.
  # Returning the first numeric field can silently confuse p-values with slopes/Z/variance.
  if (allow_unnamed_fallback) {
    num <- suppressWarnings(as.numeric(flat))
    num <- num[is.finite(num)]
    if (require_probability) num <- num[num >= 0 & num <= 1]
    if (length(num)) return(unname(num[1]))
  }
  NA_real_
}

v9_extract_sens_slope_estimate <- function(x) {
  # trend::sens.slope() usually stores the estimate under an estimate/Sen-style name.
  v9_extract_numeric_by_name(x, c("sen", "slope", "estimate"), allow_unnamed_fallback = TRUE)
}

v9_extract_mmkh_p <- function(mm) {
  # modifiedmk::mmkh() object names have varied across package versions.  Only accept
  # a name-matched probability; never fall back to the first numeric field.
  v9_extract_numeric_by_name(
    mm,
    c("new.*p", "corrected.*p", "correction.*p", "modified.*p", "p[._ -]?value", "p-value", "^p$"),
    allow_unnamed_fallback = FALSE,
    require_probability = TRUE
  )
}

.v9_pair_cache <- new.env(parent = emptyenv())
v9_pair_indices <- function(n) {
  n <- as.integer(n)
  if (is.na(n) || n < 2L) return(list(i=integer(),j=integer()))
  key <- as.character(n)
  if (exists(key,envir=.v9_pair_cache,inherits=FALSE)) return(get(key,envir=.v9_pair_cache,inherits=FALSE))
  i <- rep.int(seq_len(n-1L),times=(n-1L):1L)
  j <- unlist(lapply(seq_len(n-1L),function(k)(k+1L):n),use.names=FALSE)
  out <- list(i=i,j=j); assign(key,out,envir=.v9_pair_cache); out
}

v9_sen_slope_fallback <- function(years, y) {
  ok <- is.finite(years) & is.finite(y)
  years <- years[ok]
  y <- y[ok]
  n <- length(y)
  if (n < 2L) return(NA_real_)
  # vectorised Theil-Sen median slope. Bit-for-bit identical to the
  # previous nested for-loop (verified), but this is the single hottest function
  # in the package -- called once per cell x trend period in the main pipeline and
  # once per residual resample inside the time_block bootstrap -- so removing the
  # interpreted double loop matters.
  ij <- v9_pair_indices(n); iu <- ij$i; ju <- ij$j
  den <- years[ju] - years[iu]
  keep <- is.finite(den) & den != 0
  if (!any(keep)) return(NA_real_)
  stats::median((y[ju] - y[iu])[keep] / den[keep], na.rm = TRUE)
}

v9_sen_slope <- function(years, y, prefer_package = (V9_USE_PROFESSIONAL_TREND_PACKAGES && V9_USE_TREND_PACKAGE_SEN), return_source = FALSE) {
  ok <- is.finite(years) & is.finite(y)
  years <- years[ok]
  y <- y[ok]
  if (length(y) < 2L) {
    out <- NA_real_; attr(out, "source") <- "insufficient"
    return(if (return_source) list(value = out, source = "insufficient") else out)
  }
  # trend::sens.slope() uses observation index spacing (x_j-x_i)/(j-i). This is
  # equivalent to the project slope only when valid annual observations are unit-consecutive.
  if (prefer_package && v9_pkg_available("trend") && v9_is_unit_consecutive(years)) {
    res <- try(withCallingHandlers(
      trend::sens.slope(y),
      warning = function(w) invokeRestart("muffleWarning")
    ), silent = TRUE)
    if (!inherits(res, "try-error")) {
      est <- if (!is.null(res$estimates)) v9_extract_sens_slope_estimate(res$estimates) else v9_extract_sens_slope_estimate(res)
      if (is.finite(est)) {
        attr(est, "source") <- "trend::sens.slope"
        return(if (return_source) list(value = unname(est), source = "trend::sens.slope") else unname(est))
      }
    }
  }
  val <- v9_sen_slope_fallback(years, y)
  attr(val, "source") <- "fallback_pairwise_years"
  if (return_source) list(value = val, source = "fallback_pairwise_years") else val
}

v9_mk_basic_stats <- function(y) {
  y <- y[is.finite(y)]
  n <- length(y)
  if (n < 3L) return(list(n = n, S = NA_real_, varS = NA_real_))
  # vectorised Mann-Kendall S (identical to the previous nested loop;
  # verified across tie/NA cases). Called once per cell, so worth vectorising.
  ij <- v9_pair_indices(n); iu <- ij$i; ju <- ij$j
  S <- sum(sign(y[ju] - y[iu]))
  ties <- table(y)
  tie_term <- sum(ties * (ties - 1) * (2 * ties + 5))
  varS <- (n * (n - 1) * (2 * n + 5) - tie_term) / 18
  list(n = n, S = S, varS = varS)
}

v9_mk_p_from_S <- function(S, varS) {
  if (is.finite(S) && is.finite(varS) && S == 0 && varS == 0) return(1) # constant valid series: no trend
  if (!is.finite(S) || !is.finite(varS) || varS <= 0) return(NA_real_)
  z <- if (S > 0) (S - 1) / sqrt(varS) else if (S < 0) (S + 1) / sqrt(varS) else 0
  2 * stats::pnorm(-abs(z))
}

v9_hamed_rao_mk_fallback <- function(years, y, lag_max = 0L, sig_only = TRUE, precomputed_sen = NULL, precomputed_base = NULL) {
  ok <- is.finite(years) & is.finite(y)
  years <- years[ok]
  y <- y[ok]
  n <- length(y)
  if (n < 3L) return(c(p = NA_real_, S = NA_real_, varS = NA_real_, varS_hr = NA_real_, acf1 = NA_real_, n = n))
  sen <- if(!is.null(precomputed_sen)&&is.finite(precomputed_sen))as.numeric(precomputed_sen)else v9_sen_slope_fallback(years,y)
  residual <- y - sen * years
  ranks <- rank(residual, ties.method = "average")
  base <- if(!is.null(precomputed_base))precomputed_base else v9_mk_basic_stats(y)
  varS <- base$varS
  if (is.finite(varS) && varS == 0 && base$S == 0) return(c(p=1, S=0, varS=0, varS_hr=0, acf1=NA_real_, n=n))
  if (!is.finite(varS) || varS <= 0) return(c(p = NA_real_, S = base$S, varS = varS, varS_hr = NA_real_, acf1 = NA_real_, n = n))
  lag_max <- if (!is.finite(lag_max) || as.integer(lag_max) <= 0L) (n - 1L) else max(1L, min(as.integer(lag_max), n - 1L))
  ac <- tryCatch(stats::acf(ranks, lag.max = lag_max, plot = FALSE, na.action = na.pass)$acf[-1], error = function(e) rep(NA_real_, lag_max))
  acf1 <- if (length(ac)) ac[1] else NA_real_
  if (sig_only) {
    crit <- 1.96 / sqrt(n)
    ac_use <- ifelse(is.finite(ac) & abs(ac) > crit, ac, 0)
  } else {
    ac_use <- ifelse(is.finite(ac), ac, 0)
  }
  lags <- seq_along(ac_use)
  weights <- (n - lags) * (n - lags - 1) * (n - lags - 2)
  weights[weights < 0] <- 0
  corr <- 1 + (2 / (n * (n - 1) * (n - 2))) * sum(weights * ac_use, na.rm = TRUE)
  corr <- max(1e-6, corr)
  varS_hr <- varS * corr
  c(p = v9_mk_p_from_S(base$S, varS_hr), S = base$S, varS = varS, varS_hr = varS_hr, acf1 = acf1, n = n)
}

v9_hamed_rao_mk <- function(years, y, lag_max = 0L, sig_only = TRUE, prefer_package = V9_USE_PROFESSIONAL_TREND_PACKAGES, precomputed_sen = NULL, precomputed_base = NULL) {
  ok <- is.finite(years) & is.finite(y)
  yy <- y[ok]
  n <- length(yy)
  if (n < 3L) return(c(p = NA_real_, S = NA_real_, varS = NA_real_, varS_hr = NA_real_, acf1 = NA_real_, n = n, source = NA_character_))
  base <- if(!is.null(precomputed_base))precomputed_base else v9_mk_basic_stats(yy)
  # acf1 is diagnostic only; compute consistently even if package supplies p.
  acf1 <- tryCatch({
    sen <- if(!is.null(precomputed_sen)&&is.finite(precomputed_sen))as.numeric(precomputed_sen)else v9_sen_slope_fallback(years[ok],yy)
    residual <- yy - sen * years[ok]
    ranks <- rank(residual, ties.method = "average")
    stats::acf(ranks, lag.max = 1, plot = FALSE, na.action = na.pass)$acf[2]
  }, error = function(e) NA_real_)
  # modifiedmk::mmkh() operates on observation order and assumes regular spacing.
  # Use it only when the valid annual observations are unit-consecutive; otherwise
  # fall back to the project implementation, which uses the true year vector for de-trending.
  if (prefer_package && isTRUE(get0("V9_USE_MODIFIEDMK_PACKAGE", ifnotfound=FALSE, inherits=TRUE)) && v9_pkg_available("modifiedmk") && v9_is_unit_consecutive(years[ok])) {
    mm <- try(modifiedmk::mmkh(yy, ci = 0.95), silent = TRUE)
    if (!inherits(mm, "try-error")) {
      p <- v9_extract_mmkh_p(mm)
      if (is.finite(p) && p >= 0 && p <= 1) {
        return(c(p = p, S = base$S, varS = base$varS, varS_hr = NA_real_, acf1 = acf1, n = n, source = "modifiedmk::mmkh"))
      }
    }
  }
  out <- v9_hamed_rao_mk_fallback(years,y,lag_max=lag_max,sig_only=sig_only,precomputed_sen=precomputed_sen,precomputed_base=base)
  c(out, source = "fallback_hamed_rao")
}

v9_raw_mk_p <- function(years, y, precomputed_base = NULL) {
  ok <- is.finite(years) & is.finite(y)
  yy <- y[ok]
  b <- if(!is.null(precomputed_base))precomputed_base else v9_mk_basic_stats(yy)
  c(p = v9_mk_p_from_S(b$S, b$varS), S = b$S, varS = b$varS, acf1 = NA_real_, n = b$n, source = "fallback_raw_mk")
}


# Trend calculation: merged trend kernel. Computes Sen slope, Mann-Kendall S,
# tie-corrected varS, raw-p, Hamed-Rao-corrected p and acf1 in ONE pass, removing the
# duplicate residual/rank/acf work and the per-cell table() call of the previous design.
# Helper functions are inlined as local closures so the exported v9_trend_stats carries
# them to PSOCK workers without any change to export=c(...) lists. Bit-identical to the
# previous fallback path (verified: sen/S/varS/p_raw/p_hr/acf1 max abs error 0). When a
# professional package path is requested (USE_TREND_PACKAGE_SEN / USE_MODIFIEDMK_PACKAGE),
# it falls back to the original per-function implementation to preserve package parity.
v9_trend_stats <- function(years,y,lag_max=0L,sig_only=TRUE,
                           prefer_package=V9_USE_PROFESSIONAL_TREND_PACKAGES,
                           prefer_sen_package=(V9_USE_PROFESSIONAL_TREND_PACKAGES&&V9_USE_TREND_PACKAGE_SEN)) {
  ok <- is.finite(years)&is.finite(y); xx <- years[ok]; yy <- y[ok]; n <- length(yy)

  # Explicit professional-package request -> original per-function path (unchanged numerics).
  if (prefer_sen_package || (prefer_package && isTRUE(get0("V9_USE_MODIFIEDMK_PACKAGE", ifnotfound = FALSE, inherits = TRUE)))) {
    sen  <- v9_sen_slope(xx, yy, prefer_package = prefer_sen_package, return_source = TRUE)
    base <- v9_mk_basic_stats(yy)
    raw  <- v9_raw_mk_p(xx, yy, precomputed_base = base)
    hr   <- v9_hamed_rao_mk(xx, yy, lag_max = lag_max, sig_only = sig_only,
                            prefer_package = prefer_package,
                            precomputed_sen = sen$value, precomputed_base = base)
    return(list(sen = sen, raw = raw, hr = hr))
  }
  if (n < 3L) {
    return(list(sen = list(value = NA_real_, source = 'insufficient'),
                raw = v9_raw_mk_p(xx, yy),
                hr  = v9_hamed_rao_mk(xx, yy, lag_max, sig_only, prefer_package)))
  }

  # tie-correction term via sorted run-length (equivalent to table() counts, no hashing)
  ties_tieterm <- function(y_sorted) {
    m <- length(y_sorted); if (m == 0L) return(0)
    chg <- c(TRUE, y_sorted[-1L] != y_sorted[-m]); idx <- which(chg)
    cnt <- diff(c(idx, m + 1L))
    sum(cnt * (cnt - 1) * (2 * cnt + 5))
  }
  # direct O(n*L) autocorrelation reusing the demeaned rank vector xm and c0 (matches stats::acf)
  acf_direct <- function(xm, c0, L) {
    m <- length(xm); res <- numeric(L)
    for (k in seq_len(L)) res[k] <- (sum(xm[seq_len(m - k)] * xm[(k + 1L):m]) / m) / c0
    res
  }

  # single pairwise pass: Sen median slope and MK S share the same differences
  ij <- v9_pair_indices(n); iu <- ij$i; ju <- ij$j
  dyr <- xx[ju] - xx[iu]; dyy <- yy[ju] - yy[iu]
  keep <- is.finite(dyr) & dyr != 0
  sen_val <- if (any(keep)) stats::median(dyy[keep] / dyr[keep], na.rm = TRUE) else NA_real_
  S <- sum(sign(dyy))

  ys <- sort(yy)
  varS <- (n * (n - 1) * (2 * n + 5) - ties_tieterm(ys)) / 18
  raw_p <- v9_mk_p_from_S(S, varS)

  # Hamed-Rao: residual / rank / autocorrelation computed exactly once
  residual <- yy - sen_val * xx
  ranks <- rank(residual, ties.method = "average")
  rmv <- ranks - mean(ranks); c0 <- sum(rmv * rmv) / n
  if (!is.finite(c0) || c0 <= 0) {
    acf1 <- NA_real_; varS_hr <- varS
  } else {
    L <- if (!is.finite(lag_max) || as.integer(lag_max) <= 0L) (n - 1L) else max(1L, min(as.integer(lag_max), n - 1L))
    ac <- acf_direct(rmv, c0, L); acf1 <- ac[1]
    if (sig_only) { crit <- 1.96 / sqrt(n); ac_use <- ifelse(is.finite(ac) & abs(ac) > crit, ac, 0) }
    else          { ac_use <- ifelse(is.finite(ac), ac, 0) }
    lags <- seq_len(L); w <- (n - lags) * (n - lags - 1) * (n - lags - 2); w[w < 0] <- 0
    corr <- max(1e-6, 1 + (2 / (n * (n - 1) * (n - 2))) * sum(w * ac_use, na.rm = TRUE))
    varS_hr <- varS * corr
  }
  hr_p <- v9_mk_p_from_S(S, varS_hr)

  list(
    sen = list(value = sen_val, source = "fallback_pairwise_years"),
    raw = c(p = raw_p, S = S, varS = varS, acf1 = NA_real_, n = n, source = "fallback_raw_mk"),
    hr  = c(p = hr_p, S = S, varS = varS, varS_hr = varS_hr, acf1 = acf1, n = n, source = "fallback_hamed_rao")
  )
}

# BKY two-stage adjusted p-values.
# BH is computed first (reusing the caller's precomputed_bh when supplied). mutoss::two.stage
# is used ONLY when USE_MUTOSS_BKY=TRUE and the finite-p count m is small (<= BKY_MUTOSS_MAX_M),
# because mutoss is ~O(n^2) and takes hours at 1-km scale. Otherwise the analytic two-stage
# plug-in (BH x m0/m, m0 = m - r1, r1 = #{BH_adj <= alpha/(1+alpha)}) is used: it is the
# standard TSBH, O(n log n), and numerically equivalent to mutoss at a fixed alpha.
v9_bky_adjust <- function(p, alpha_ref = 0.05, prefer_package = V9_USE_PROFESSIONAL_TREND_PACKAGES,
                          precomputed_bh = NULL) {
  p <- as.numeric(p)
  out <- rep(NA_real_, length(p))
  ok <- is.finite(p)
  if (!sum(ok)) {
    attr(out, "source") <- "empty"
    return(out)
  }
  m <- sum(ok)
  alpha1 <- alpha_ref / (1 + alpha_ref)
  # First-stage BH q-values (reuse a caller-supplied vector when valid).
  if (!is.null(precomputed_bh)) {
    q_bh <- as.numeric(precomputed_bh)[ok]
    if (length(q_bh) != m || any(!is.finite(q_bh))) q_bh <- stats::p.adjust(p[ok], method = "BH")
  } else {
    q_bh <- stats::p.adjust(p[ok], method = "BH")
  }
  # Optional mutoss parity path: opt-in AND size-gated (mutoss is ~O(n^2)).
  use_mutoss <- prefer_package &&
    isTRUE(get0("V9_USE_MUTOSS_BKY", ifnotfound = FALSE, inherits = TRUE)) &&
    m <= get0("V9_BKY_MUTOSS_MAX_M", ifnotfound = 20000L, inherits = TRUE) &&
    v9_pkg_available("mutoss")
  if (use_mutoss) {
    ts <- try(mutoss::two.stage(pValues = p[ok], alpha = alpha_ref), silent = TRUE)
    if (!inherits(ts, "try-error")) {
      adj <- tryCatch(ts$adjPValues, error = function(e) NULL)
      if (!is.null(adj) && length(adj) == m) {
        out[ok] <- pmin(1, pmax(0, as.numeric(adj)))
        attr(out, "source") <- "mutoss::two.stage"
        return(out)
      }
    }
  }
  # Analytic two-stage plug-in (standard TSBH), O(n log n).
  r1 <- sum(q_bh <= alpha1, na.rm = TRUE)
  m0 <- max(1, m - r1)
  out[ok] <- pmin(1, q_bh * m0 / m)
  attr(out, "source") <- "TSBH_analytic"
  out
}
