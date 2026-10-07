harmonic_number <- function(m) {
  m <- as.numeric(m)
  if (!is.finite(m) || m < 1 || m != floor(m)) stop("m must be a positive integer.", call. = FALSE)
  digamma(m + 1) - digamma(1)
}

discrete_fdr_lookup <- function(counts, B, method = c("BH", "BY")) {
  method <- match.arg(toupper(method), c("BH", "BY")); counts <- as.numeric(counts)
  if (length(counts) != B + 1L) stop("counts must have length B+1 for k=0..B.", call. = FALSE)
  if (any(counts < 0 | !is.finite(counts))) stop("counts must be finite nonnegative.", call. = FALSE)
  m <- sum(counts); q <- rep(NA_real_, B + 1L)
  if (m == 0) return(q)
  p <- (1 + 0:B) / (B + 1); present <- which(counts > 0); rank_max <- cumsum(counts)[present]
  factor <- if (method == "BY") harmonic_number(m) else 1
  raw <- factor * m * p[present] / rank_max
  q[present] <- pmin(1, rev(cummin(rev(raw))))
  # Fill absent levels from the next observed level only for diagnostics; no raster cell has an absent k.
  last <- NA_real_
  for (i in rev(seq_len(B + 1L))) { if (is.finite(q[[i]])) last <- q[[i]] else if (is.finite(last)) q[[i]] <- last }
  q
}

count_mc_exceed_raster <- function(r, B) {
  assert_packages("terra")
  counts <- numeric(B + 1L); bs <- terra::blocks(r)
  for (i in seq_len(bs$n)) {
    v <- terra::values(r, row = bs$row[[i]], nrows = bs$nrows[[i]], mat = FALSE)
    v <- as.integer(round(v[is.finite(v)])); v <- v[v >= 0L & v <= B]
    if (length(v)) counts <- counts + tabulate(v + 1L, nbins = B + 1L)
  }
  counts
}

map_fdr_lookup_raster <- function(mc_exceed_raster, lookup, filename, overwrite = FALSE) {
  assert_packages("terra")
  r <- mc_exceed_raster; out <- terra::rast(r); names(out) <- "q"
  wb <- terra::writeStart(out, filename = filename, overwrite = overwrite, wopt = list(datatype = "FLT8S"))
  write_open <- TRUE; on.exit(if (write_open) try(terra::writeStop(out), silent = TRUE), add = TRUE)
  for (i in seq_len(wb$n)) {
    v <- terra::values(r, row = wb$row[[i]], nrows = wb$nrows[[i]], mat = FALSE)
    z <- rep(NA_real_, length(v)); ok <- is.finite(v); kk <- as.integer(round(v[ok])); valid <- kk >= 0L & kk < length(lookup)
    zz <- rep(NA_real_, length(kk)); zz[valid] <- lookup[kk[valid] + 1L]; z[ok] <- zz
    terra::writeValues(out, z, wb$row[[i]], wb$nrows[[i]])
  }
  out <- terra::writeStop(out); write_open <- FALSE; out
}

fdr_raster_from_exceed <- function(mc_exceed_raster, B, method, filename, overwrite = FALSE) {
  counts <- count_mc_exceed_raster(mc_exceed_raster, B); lookup <- discrete_fdr_lookup(counts, B, method)
  q <- map_fdr_lookup_raster(mc_exceed_raster, lookup, filename, overwrite)
  list(raster = q, counts = counts, lookup = lookup, n_tests = sum(counts))
}

fdr_p_cutoff <- function(counts, lookup, B, alpha) {
  k <- which(counts > 0 & is.finite(lookup) & lookup <= alpha) - 1L
  if (!length(k)) return(NA_real_)
  max((1 + k) / (B + 1))
}

validate_discrete_fdr <- function(k, B, method = c("BH", "BY"), tol = 1e-12) {
  method <- match.arg(method); counts <- tabulate(k + 1L, nbins = B + 1L); lookup <- discrete_fdr_lookup(counts, B, method)
  p <- (1 + k) / (B + 1); got <- lookup[k + 1L]; ref <- stats::p.adjust(p, method = method)
  max(abs(got - ref)) <= tol
}
