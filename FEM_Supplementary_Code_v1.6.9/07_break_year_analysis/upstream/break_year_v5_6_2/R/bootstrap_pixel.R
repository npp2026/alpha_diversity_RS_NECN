simulate_ar1_with_marginal_sd <- function(n, rho, sd_marginal) {
  if (!is.finite(sd_marginal) || sd_marginal < 0) stop("sd_marginal must be finite and >= 0.", call. = FALSE)
  if (sd_marginal == 0) return(rep(0, n))
  as.numeric(simulate_ar1_matrix(1L, n, rho, sd_marginal))
}

bootstrap_one_pixel <- function(y, years, observed_fit, null_sorted, rho, fdr_p_cutoff, cfg, response, cell_id) {
  pcfg <- cfg$bootstrap$pixel; B <- as.integer(pcfg$B)
  # Marginal residual SD with hinge residual degrees of freedom. This is the
  # Gaussian AR(1) bootstrap scale, not the raw sample SD of residuals.
  sig <- sqrt(max(observed_fit$hinge_sse, 0) / max(length(y) - 3L, 1L))
  if (!is.finite(sig)) return(NULL)
  seed <- seed_from_key(cfg$rng$master_seed, "pixel_boot", response, cell_id)
  set_rng(cfg$rng$kind, seed)
  tau <- rep(NA_real_, B); rejected <- logical(B); recovery <- logical(B); testable <- logical(B)
  m <- as.integer(cfg$analysis$primary$min_segment)
  slope_tol <- as.numeric(cfg$analysis$slope_tol[[response]] %||% 0)
  candidates <- years[candidate_indices(length(years), m, m)]
  for (b in seq_len(B)) {
    ys <- observed_fit$fitted + simulate_ar1_with_marginal_sd(length(y), rho, sig)
    ft <- fit_break_core(ys, years, m, m, slope_tol, cfg$analysis$tie_rule)
    if (isTRUE(ft$testable)) {
      testable[[b]] <- TRUE
      pb <- mc_pvalue(ft$supF, null_sorted)
      # Freeze the primary whole-map BH/BY p-value cutoff. This is a stability
      # frequency under the fitted model, not a posterior probability.
      rejected[[b]] <- is.finite(pb) && is.finite(fdr_p_cutoff) && pb <= fdr_p_cutoff
      recovery[[b]] <- rejected[[b]] && isTRUE(ft$direction_ok)
      if (recovery[[b]]) tau[[b]] <- ft$break_year
    }
  }
  valid_tau <- tau[is.finite(tau)]; min_rep <- as.integer(pcfg$min_valid_replicates_for_ci)
  q <- if (length(valid_tau) >= min_rep) safe_quantile_discrete(valid_tau, c(0.025, 0.5, 0.975)) else c(NA, NA, NA)
  data.frame(
    cell = as.integer(cell_id), bootstrap_testable_frequency = mean(testable),
    bootstrap_rejection_frequency = mean(rejected), bootstrap_recovery_frequency = mean(recovery),
    n_recovery_boot = sum(recovery), break_year_q025 = q[[1L]], break_year_median = q[[2L]], break_year_q975 = q[[3L]],
    break_year_sd = if (length(valid_tau) > 1L) stats::sd(valid_tau) else NA_real_,
    break_year_IQR = if (length(valid_tau)) diff(safe_quantile_discrete(valid_tau, c(0.25, 0.75))) else NA_real_,
    break_year_CI_width = if (all(is.finite(q[c(1,3)]))) q[[3L]] - q[[1L]] else NA_real_,
    boundary_hit_frequency = if (length(valid_tau)) mean(valid_tau %in% range(candidates)) else NA_real_,
    ci_status = if (length(valid_tau) >= min_rep) "ok" else "unstable", stringsAsFactors = FALSE)
}

cells_from_mask <- function(mask) {
  assert_packages("terra")
  bs <- terra::blocks(mask); nc <- terra::ncol(mask); out <- vector("list", bs$n)
  for (i in seq_len(bs$n)) {
    v <- terra::values(mask, row = bs$row[[i]], nrows = bs$nrows[[i]], mat = FALSE)
    hit <- which(is.finite(v) & v == 1)
    first_cell <- (bs$row[[i]] - 1L) * nc + 1L
    out[[i]] <- if (length(hit)) first_cell + hit - 1L else integer()
  }
  as.integer(unlist(out, use.names = FALSE))
}

metric_raster_from_cells <- function(template, cells, values, filename, overwrite = FALSE) {
  assert_packages("terra")
  if (length(cells) != length(values)) stop("cells and values differ in length.", call. = FALSE)
  ord <- order(cells); cells <- as.integer(cells[ord]); values <- as.numeric(values[ord])
  out <- terra::rast(template); names(out) <- tools::file_path_sans_ext(basename(filename))
  wb <- terra::writeStart(out, filename = filename, overwrite = overwrite, wopt = list(datatype = "FLT8S"))
  write_open <- TRUE; on.exit(if (write_open) try(terra::writeStop(out), silent = TRUE), add = TRUE)
  nc <- terra::ncol(out); pos <- 1L
  for (i in seq_len(wb$n)) {
    ncell_block <- wb$nrows[[i]] * nc; first <- (wb$row[[i]] - 1L) * nc + 1L; last <- first + ncell_block - 1L
    z <- rep(NA_real_, ncell_block)
    while (pos <= length(cells) && cells[[pos]] < first) pos <- pos + 1L
    j <- pos
    while (j <= length(cells) && cells[[j]] <= last) { z[cells[[j]] - first + 1L] <- values[[j]]; j <- j + 1L }
    pos <- j; terra::writeValues(out, z, wb$row[[i]], wb$nrows[[i]])
  }
  out <- terra::writeStop(out); write_open <- FALSE; out
}

pixel_bootstrap_batch_task <- function(task, state) {
  response <- task$response; cc <- as.integer(task$cells)
  # A PSOCK worker may process several batches. Cache its read-only response
  # stack locally so we do not reopen 20 GeoTIFF layers for every batch.
  if (!exists(".V56_PIXEL_STACK_CACHE", envir = .GlobalEnv, inherits = FALSE))
    assign(".V56_PIXEL_STACK_CACHE", new.env(parent = emptyenv()), envir = .GlobalEnv)
  cache <- get(".V56_PIXEL_STACK_CACHE", envir = .GlobalEnv, inherits = FALSE)
  if (!exists(response, envir = cache, inherits = FALSE))
    assign(response, load_response_stack(state$inventory$files[[response]], state$cfg), envir = cache)
  st <- get(response, envir = cache, inherits = FALSE)
  vals0 <- terra::extract(st, cc, raw = TRUE)
  vals <- if (is.null(dim(vals0))) matrix(as.numeric(vals0), nrow = 1L) else as.matrix(vals0)
  if (nrow(vals) != length(cc) || ncol(vals) != length(state$years)) {
    stop("Unexpected terra::extract shape in pixel bootstrap batch.", call. = FALSE)
  }
  rows <- vector("list", length(cc))
  for (i in seq_along(cc)) {
    y <- as.numeric(vals[i, ])
    obs <- fit_break_core(y, state$years, state$cfg$analysis$primary$min_segment, state$cfg$analysis$primary$min_segment,
                          state$cfg$analysis$slope_tol[[response]] %||% 0, state$cfg$analysis$tie_rule)
    rows[[i]] <- bootstrap_one_pixel(y, state$years, obs, state$null_by_response[[response]],
                                     as.numeric(state$rho_by_response[[response]]),
                                     state$cutoff_by_response[[response]], state$cfg, response, cc[[i]])
  }
  z <- do.call(rbind, rows)
  if (!is.null(z) && nrow(z)) z$response <- response
  z
}

run_pixel_bootstrap <- function(cfg, inventory, primary, rho_by_response, ctx) {
  if (!isTRUE(cfg$bootstrap$pixel$enabled)) return(NULL)
  yrs <- as.numeric(unlist(cfg$data$years)); root <- ensure_dir(file.path(ctx$root, "uncertainty", "pixel"))
  batch_size <- as.integer(cfg$bootstrap$pixel$batch_size %||% 2000L)
  max_target <- as.numeric(cfg$bootstrap$pixel$max_target_pixels %||% Inf)

  tasks <- list(); response_cells <- list(); task_id <- 0L
  null_by_response <- list(); cutoff_by_response <- list()
  for (response in names(primary$results)) {
    res <- primary$results[[response]]; cells <- cells_from_mask(res$recovery)
    response_cells[[response]] <- cells
    if (!length(cells)) next
    if (length(cells) > max_target) stop("Pixel bootstrap target exceeds bootstrap.pixel.max_target_pixels for ", response,
                                         ": ", length(cells), " > ", max_target, call. = FALSE)
    null_by_response[[response]] <- res$calibration$null
    cutoff_by_response[[response]] <- res$frozen_fdr_p_cutoff
    for (from in seq(1L, length(cells), by = batch_size)) {
      task_id <- task_id + 1L
      idx <- from:min(length(cells), from + batch_size - 1L)
      tasks[[task_id]] <- list(task_id = task_id, response = response, cells = as.integer(cells[idx]))
    }
  }
  if (!length(tasks)) return(list())

  state <- list(cfg = cfg, inventory = inventory, years = yrs, null_by_response = null_by_response,
                cutoff_by_response = cutoff_by_response, rho_by_response = rho_by_response)
  pieces <- if (parallel_stage_enabled(cfg, "pixel_bootstrap") && length(tasks) > 1L) {
    psock_task_lapply(cfg, tasks, "pixel_bootstrap_batch_task", state, stage = "pixel_bootstrap", log_file = ctx$log)
  } else lapply(tasks, pixel_bootstrap_batch_task, state = state)

  all_rows <- list()
  for (response in names(primary$results)) {
    idx <- which(vapply(tasks, function(t) identical(t$response, response), logical(1L)))
    if (!length(idx)) next
    tab <- do.call(rbind, pieces[idx]); if (is.null(tab) || !nrow(tab)) next
    # PSOCK preserves task order; sorting by cell additionally makes output
    # invariant to future scheduler changes and batch_size choices.
    tab <- tab[order(tab$cell), , drop = FALSE]
    write_csv(tab, file.path(root, paste0(response, "_pixel_bootstrap.csv")))
    template <- primary$results[[response]]$layers$break_year
    for (nm in c("bootstrap_rejection_frequency", "bootstrap_recovery_frequency", "break_year_median", "break_year_CI_width", "boundary_hit_frequency"))
      metric_raster_from_cells(template, tab$cell, tab[[nm]], file.path(root, paste0(response, "_", nm, ".tif")), isTRUE(cfg$runtime$overwrite))
    all_rows[[response]] <- tab
  }
  if (length(all_rows)) write_csv(do.call(rbind, all_rows), file.path(root, "pixel_bootstrap_all.csv"))
  all_rows
}
