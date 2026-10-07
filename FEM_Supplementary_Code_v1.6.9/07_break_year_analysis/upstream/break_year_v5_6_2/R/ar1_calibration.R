simulate_ar1_matrix <- function(B, n, rho, marginal_sd = 1) {
  if (!is.finite(rho) || abs(rho) >= 1) stop("rho must satisfy abs(rho) < 1.", call. = FALSE)
  e <- matrix(0, nrow = B, ncol = n)
  e[, 1L] <- stats::rnorm(B, sd = marginal_sd)
  if (n > 1L) {
    innov_sd <- marginal_sd * sqrt(max(1 - rho^2, 0))
    for (tt in 2:n) e[, tt] <- rho * e[, tt - 1L] + stats::rnorm(B, sd = innov_sd)
  }
  e
}

residual_maker <- function(X) {
  q <- qr(X)
  if (q$rank < ncol(X)) stop("Singular design matrix in residual_maker().", call. = FALSE)
  diag(nrow(X)) - X %*% solve(crossprod(X), t(X))
}

supf_for_matrix <- function(Y, years, min_left, min_right) {
  Y <- as.matrix(Y); n <- ncol(Y)
  if (length(years) != n) stop("Y columns and years length differ.", call. = FALSE)
  if (any(!is.finite(Y))) stop("supf_for_matrix requires finite complete series.", call. = FALSE)
  t <- as.numeric(years - years[[1L]])
  M0 <- residual_maker(cbind(1, t))
  sse0 <- rowSums((Y %*% M0) * Y)
  cand <- candidate_indices(n, min_left, min_right)
  if (!length(cand)) stop("No candidate break years.", call. = FALSE)
  best <- rep(Inf, nrow(Y))
  for (idx in cand) {
    tau <- years[[idx]]
    X <- cbind(1, t, pmax(0, years - tau))
    M <- residual_maker(X)
    sse <- rowSums((Y %*% M) * Y)
    best <- pmin(best, sse)
  }
  numer <- sse0 - best
  # Row-wise SSE tolerance mirrors fit_break_core(). A fixed O(sqrt(eps))
  # threshold is not scale invariant and can misclassify small-magnitude data.
  centered <- Y - rowMeans(Y)
  scale <- rowSums(centered^2)
  amp <- pmax(1, apply(abs(Y), 1L, max))
  tol <- pmax(100 * .Machine$double.eps * scale,
              n * (.Machine$double.eps * amp)^2 * 1000,
              .Machine$double.xmin)
  numer[numer < 0 & abs(numer) <= tol] <- 0
  out <- numeric(nrow(Y))
  exact <- best <= tol
  out[exact & numer <= tol] <- 0
  out[exact & numer > tol] <- Inf
  regular <- !exact
  out[regular] <- (n - 3L) * numer[regular] / best[regular]
  tiny_neg <- out < 0 & abs(out) <= 100 * .Machine$double.eps
  out[tiny_neg] <- 0
  out[!is.finite(out) & !is.infinite(out)] <- NA_real_
  out
}

supf_mc_chunk_task <- function(task, state) {
  chunk_seed <- seed_from_key(state$base_seed, "mc_chunk", task$chunk_id)
  set_rng(state$rng_kind, chunk_seed)
  Y <- simulate_ar1_matrix(task$n, length(state$years), state$rho, marginal_sd = 1)
  supf_for_matrix(Y, state$years, state$min_left, state$min_right)
}

simulate_supf_null <- function(years, min_left, min_right, rho, B = 200000L,
                               seed = 12345L, chunk_size = 10000L, rng_kind = "L'Ecuyer-CMRG",
                               cfg = NULL, log_file = NULL) {
  years <- as.numeric(years); B <- as.integer(B); chunk_size <- as.integer(chunk_size)
  if (B < 1L || chunk_size < 1L) stop("B and chunk_size must be positive integers.", call. = FALSE)
  starts <- seq.int(1L, B, by = chunk_size)
  tasks <- lapply(seq_along(starts), function(i) list(
    chunk_id = as.integer(i),
    n = as.integer(min(chunk_size, B - starts[[i]] + 1L))
  ))
  state <- list(years = years, min_left = as.integer(min_left), min_right = as.integer(min_right),
                rho = as.numeric(rho), base_seed = as.integer(seed), rng_kind = as.character(rng_kind))
  pieces <- if (!is.null(cfg) && parallel_stage_enabled(cfg, "supf_mc")) {
    psock_task_lapply(cfg, tasks, "supf_mc_chunk_task", state, stage = "supf_mc", log_file = log_file)
  } else {
    lapply(tasks, supf_mc_chunk_task, state = state)
  }
  ans <- as.numeric(unlist(pieces, use.names = FALSE))
  if (length(ans) != B) stop("SupF MC chunk assembly returned the wrong number of draws.", call. = FALSE)
  sort(ans)
}

# Number of Monte Carlo statistics >= observed. With left.open=TRUE,
# findInterval returns the count strictly below x, so ties are included.
mc_exceed_ge <- function(observed, null_sorted) {
  if (!length(null_sorted) || is.na(observed)) return(NA_integer_)
  if (is.infinite(observed)) {
    if (observed < 0) return(as.integer(length(null_sorted)))
    # Count +Inf null values without a full scan; finite SupF values are <= .Machine$double.xmax.
    return(as.integer(length(null_sorted) - findInterval(.Machine$double.xmax, null_sorted)))
  }
  as.integer(length(null_sorted) - findInterval(observed, null_sorted, left.open = TRUE))
}

mc_pvalue_from_exceed <- function(k, B) {
  if (!is.finite(k)) return(NA_real_)
  (1 + as.numeric(k)) / (as.numeric(B) + 1)
}

mc_pvalue <- function(observed, null_sorted) mc_pvalue_from_exceed(mc_exceed_ge(observed, null_sorted), length(null_sorted))

calibration_summary <- function(null_sorted, years, min_left, min_right, rho, B, seed) {
  qs <- safe_quantile_discrete(null_sorted, c(0.90, 0.95, 0.99))
  data.frame(
    min_left = min_left, min_right = min_right, rho = rho, B = B, seed = seed,
    candidate_first = years[[min_left]], candidate_last = years[[length(years) - min_right]],
    critical_90 = qs[[1L]], critical_95 = qs[[2L]], critical_99 = qs[[3L]],
    null_mean = mean(null_sorted[is.finite(null_sorted)]), null_median = stats::median(null_sorted[is.finite(null_sorted)]),
    inf_count = sum(is.infinite(null_sorted)), stringsAsFactors = FALSE
  )
}

get_or_create_calibration <- function(cfg, years, min_segment, rho, response, ctx) {
  ensure_dir(cfg$supf$cache_dir)
  B <- supf_mc_B_for_response(cfg, response)
  chunk <- as.integer(cfg$supf$mc_chunk); rng_kind <- as.character(cfg$rng$kind); master_seed <- as.integer(cfg$rng$master_seed)
  r_version <- paste(R.version$major, R.version$minor, sep = ".")
  key <- substr(hash_text(c("v5.6-supf-chunk-key-v1", cfg$.code_hash %||% "no-code-hash", paste(years, collapse = ","), min_segment, rho, B,
                          chunk, cfg$supf$p_rule, rng_kind, master_seed, r_version)), 1L, 16L)
  cache <- file.path(cfg$supf$cache_dir, paste0("supf_", response, "_", key, ".rds"))
  caldir <- ensure_dir(file.path(ctx$root, "calibration"))
  if (isTRUE(cfg$supf$cache) && file.exists(cache)) {
    obj <- safe_read_rds_cache(cache)
    if (!is.null(obj)) {
      write_csv(obj$summary, file.path(caldir, paste0("supf_", response, "_L", min_segment, "_", key, ".csv")))
      return(obj)
    }
  }
  seed <- seed_from_key(cfg$rng$master_seed, "supf", response, min_segment, sprintf("%.8f", rho))
  null <- simulate_supf_null(years, min_segment, min_segment, rho, B, seed, chunk, rng_kind, cfg = cfg, log_file = ctx$log)
  obj <- list(null = null, summary = calibration_summary(null, years, min_segment, min_segment, rho, B, seed), key = key)
  if (isTRUE(cfg$supf$cache)) atomic_save_rds(obj, cache)
  write_csv(obj$summary, file.path(caldir, paste0("supf_", response, "_L", min_segment, "_", key, ".csv")))
  obj
}
