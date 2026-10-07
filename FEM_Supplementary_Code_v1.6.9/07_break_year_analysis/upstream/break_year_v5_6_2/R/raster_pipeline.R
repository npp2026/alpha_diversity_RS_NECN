write_named_layers <- function(raw, out_dir, overwrite = FALSE) {
  assert_packages("terra"); ensure_dir(out_dir)
  specs <- list(n_obs = "INT2S", supF = "FLT8S", break_year = "INT2S", slope_before = "FLT8S", slope_after = "FLT8S",
                direction_ok = "INT1U", mc_exceed = "INT4S", p_mc = "FLT8S", fit_status = "INT1U", tie_count = "INT2S")
  out <- list()
  for (nm in names(specs)) {
    f <- file.path(out_dir, paste0(nm, ".tif"))
    out[[nm]] <- terra::writeRaster(raw[[nm]], f, overwrite = overwrite, datatype = specs[[nm]])
  }
  out
}

raster_true_count <- function(r) {
  z <- terra::global(r == 1, "sum", na.rm = TRUE); as.numeric(z[[1L, 1L]])
}

recovery_break_counts <- function(break_year, recovery) {
  assert_packages("terra")
  if (!terra::compareGeom(break_year, recovery, stopOnError = FALSE)) stop("break/recovery raster geometry differs.", call. = FALSE)
  bs <- terra::blocks(break_year); counts <- new.env(parent = emptyenv())
  for (i in seq_len(bs$n)) {
    by <- terra::values(break_year, row = bs$row[[i]], nrows = bs$nrows[[i]], mat = FALSE)
    rc <- terra::values(recovery, row = bs$row[[i]], nrows = bs$nrows[[i]], mat = FALSE)
    ok <- is.finite(by) & is.finite(rc) & rc == 1
    z <- as.integer(round(by[ok]))
    if (length(z)) {
      tb <- table(z)
      for (nm in names(tb)) assign(nm, (if (exists(nm, envir = counts, inherits = FALSE)) get(nm, envir = counts) else 0) + as.numeric(tb[[nm]]), envir = counts)
    }
  }
  nms <- ls(counts, all.names = TRUE)
  if (!length(nms)) return(data.frame(break_year = integer(), count = numeric()))
  d <- data.frame(break_year = as.integer(nms), count = vapply(nms, function(nm) get(nm, envir = counts), numeric(1L)))
  d[order(d$break_year), , drop = FALSE]
}

run_response_analysis <- function(stack, response, years, min_segment, rho, cfg, out_dir, ctx) {
  assert_packages("terra"); ensure_dir(out_dir)
  cal <- get_or_create_calibration(cfg, years, min_segment, rho, response, ctx); null <- cal$null
  slope_tol <- as.numeric(cfg$analysis$slope_tol[[response]] %||% 0); tie_rule <- cfg$analysis$tie_rule %||% "earliest"
  nyr <- length(years); require_complete <- isTRUE(cfg$data$require_complete_series)
  one_cell <- function(v) {
    if (require_complete && sum(is.finite(v)) != nyr) return(c(n_obs = sum(is.finite(v)), supF = NA, break_year = NA, slope_before = NA, slope_after = NA, direction_ok = 0, mc_exceed = NA, p_mc = NA, fit_status = 0, tie_count = 0))
    fit_break_vector(v, years, min_segment, min_segment, slope_tol, tie_rule, null)
  }
  # terra::app may call a custom function with either one cell vector or a
  # matrix whose rows are cells and columns are layers. Support both contracts.
  fun <- function(v) {
    if (is.matrix(v)) return(t(vapply(seq_len(nrow(v)), function(i) one_cell(v[i, ]), numeric(10L))))
    one_cell(v)
  }
  # Deliberately keep terra::app single-worker. v5.6 parallelizes independent
  # responses/scenarios outside terra to avoid nested worker pools and helper
  # export problems in custom terra callbacks.
  raw_path <- file.path(out_dir, "core_multiband.tif")
  raw <- terra::app(stack, fun = fun, cores = 1L, filename = raw_path,
                    overwrite = isTRUE(cfg$runtime$overwrite), wopt = list(datatype = "FLT8S"))
  names(raw) <- c("n_obs", "supF", "break_year", "slope_before", "slope_after", "direction_ok", "mc_exceed", "p_mc", "fit_status", "tie_count")
  layers <- write_named_layers(raw, out_dir, isTRUE(cfg$runtime$overwrite))

  B <- length(null); alpha <- as.numeric(cfg$multiple_testing$alpha); q <- list(); fdr <- list()
  methods <- unique(toupper(c(cfg$multiple_testing$primary_method, unlist(cfg$multiple_testing$sensitivity_methods), unlist(cfg$sensitivity$multiplicity$methods))))
  for (method in methods) {
    qfile <- file.path(out_dir, paste0("q_", method, ".tif"))
    fq <- fdr_raster_from_exceed(layers$mc_exceed, B, method, qfile, isTRUE(cfg$runtime$overwrite))
    q[[method]] <- fq$raster; fdr[[method]] <- fq
    write_csv(data.frame(k = 0:B, p_mc = (1 + 0:B)/(B+1), count = fq$counts, q = fq$lookup), file.path(out_dir, paste0("fdr_lookup_", method, ".csv")))
  }
  primary_method <- toupper(cfg$multiple_testing$primary_method)
  sig <- terra::ifel(q[[primary_method]] <= alpha, 1, 0); sig <- terra::mask(sig, layers$p_mc)
  rec <- terra::ifel(sig == 1 & layers$direction_ok == 1, 1, 0); rec <- terra::mask(rec, layers$p_mc)
  sig <- terra::writeRaster(sig, file.path(out_dir, "significant_break.tif"), overwrite = isTRUE(cfg$runtime$overwrite), datatype = "INT1U")
  rec <- terra::writeRaster(rec, file.path(out_dir, "recovery_break.tif"), overwrite = isTRUE(cfg$runtime$overwrite), datatype = "INT1U")

  n_test <- as.numeric(terra::global(!is.na(layers$p_mc), "sum", na.rm = TRUE)[[1L, 1L]])
  n_sig <- raster_true_count(sig); n_rec <- raster_true_count(rec)
  dist <- recovery_break_counts(layers$break_year, rec)
  candidates <- years[candidate_indices(length(years), min_segment, min_segment)]
  if (nrow(dist) && n_rec > 0) {
    dist$proportion <- dist$count / n_rec
    dist$equal_candidate_ratio <- dist$proportion / (1 / length(candidates))
  } else { dist$proportion <- numeric(nrow(dist)); dist$equal_candidate_ratio <- numeric(nrow(dist)) }
  write_csv(dist, file.path(out_dir, "break_year_distribution.csv"))
  cutoff <- fdr_p_cutoff(fdr[[primary_method]]$counts, fdr[[primary_method]]$lookup, B, alpha)
  lower_share <- if (n_rec > 0) sum(dist$count[dist$break_year == min(candidates)]) / n_rec else NA_real_
  upper_share <- if (n_rec > 0) sum(dist$count[dist$break_year == max(candidates)]) / n_rec else NA_real_
  mc_min_p <- 1/(B+1); mt_method <- toupper(primary_method); dep_factor <- if (mt_method == "BY" && n_test > 0) harmonic_number(n_test) else 1
  min_rank_at_mc_floor <- if (n_test > 0) ceiling(mc_min_p * n_test * dep_factor / alpha) else NA_real_
  summary <- data.frame(response = response, min_segment = min_segment, rho = rho, B = B,
                        eligible_pixels = n_test, significant_breaks = n_sig, recovery_breaks = n_rec,
                        mc_min_p = mc_min_p, min_rank_at_mc_floor = min_rank_at_mc_floor,
                        recovery_pct_eligible = if (n_test > 0) 100 * n_rec / n_test else NA_real_,
                        modal_break_year = if (nrow(dist)) dist$break_year[[which.max(dist$count)]] else NA_real_,
                        lower_boundary_share = lower_share, upper_boundary_share = upper_share,
                        primary_fdr_method = primary_method, alpha = alpha, frozen_fdr_p_cutoff = cutoff,
                        stringsAsFactors = FALSE)
  write_csv(summary, file.path(out_dir, "summary.csv"))
  list(out_dir = out_dir, layers = layers, q = q, fdr = fdr, significant = sig, recovery = rec,
       summary = summary, calibration = cal, frozen_fdr_p_cutoff = cutoff, break_distribution = dist)
}

response_result_descriptor <- function(res) {
  list(
    out_dir = res$out_dir,
    q_methods = names(res$q),
    fdr = lapply(res$fdr, function(z) list(counts = z$counts, lookup = z$lookup, n_tests = z$n_tests)),
    summary = res$summary,
    calibration = res$calibration,
    frozen_fdr_p_cutoff = res$frozen_fdr_p_cutoff,
    break_distribution = res$break_distribution
  )
}

hydrate_response_result <- function(desc) {
  assert_packages("terra")
  out_dir <- desc$out_dir
  layer_names <- c("n_obs", "supF", "break_year", "slope_before", "slope_after", "direction_ok", "mc_exceed", "p_mc", "fit_status", "tie_count")
  layers <- stats::setNames(lapply(layer_names, function(nm) terra::rast(file.path(out_dir, paste0(nm, ".tif")))), layer_names)
  q <- stats::setNames(lapply(desc$q_methods, function(method) terra::rast(file.path(out_dir, paste0("q_", method, ".tif")))), desc$q_methods)
  fdr <- desc$fdr
  for (method in names(fdr)) fdr[[method]]$raster <- q[[method]]
  list(
    out_dir = out_dir, layers = layers, q = q, fdr = fdr,
    significant = terra::rast(file.path(out_dir, "significant_break.tif")),
    recovery = terra::rast(file.path(out_dir, "recovery_break.tif")),
    summary = desc$summary, calibration = desc$calibration,
    frozen_fdr_p_cutoff = desc$frozen_fdr_p_cutoff,
    break_distribution = desc$break_distribution
  )
}

run_response_scenario_task <- function(task, state) {
  response <- task$response
  st <- load_response_stack(state$inventory$files[[response]], state$cfg)
  rdir <- ensure_dir(file.path(state$sdir, response))
  res <- run_response_analysis(st, response, state$years, state$min_segment,
                               as.numeric(state$rho_by_response[[response]]),
                               state$cfg, rdir, state$ctx)
  response_result_descriptor(res)
}

scenario_result_descriptor <- function(x) {
  list(id = x$id, dir = x$dir, results = x$results, summary = x$summary)
}

hydrate_scenario_descriptor <- function(desc) {
  list(id = desc$id, dir = desc$dir,
       results = lapply(desc$results, hydrate_response_result),
       summary = desc$summary)
}

run_scenario <- function(cfg, inventory, rho_by_response, min_segment, family, root_dir, ctx,
                         fdr_label = NULL, return_descriptor = FALSE) {
  yrs <- as.numeric(unlist(cfg$data$years)); fdr_label <- fdr_label %||% cfg$multiple_testing$primary_method
  sid <- scenario_id(family, min_segment, rho_by_response, fdr_label); sdir <- ensure_dir(file.path(root_dir, sid))
  responses <- names(inventory$files)
  tasks <- lapply(responses, function(nm) list(response = nm))
  # When response-level outer parallelism is selected, pre-warm calibration
  # caches on the master. This lets deterministic MC chunks use the full MC
  # worker pool without creating nested workers inside response workers.
  if (length(tasks) > 1L && parallel_stage_enabled(cfg, "response") &&
      parallel_stage_enabled(cfg, "supf_mc") && isTRUE(cfg$supf$cache)) {
    for (nm in responses) invisible(get_or_create_calibration(cfg, yrs, min_segment, as.numeric(rho_by_response[[nm]]), nm, ctx))
  }
  state <- list(cfg = cfg, inventory = inventory, rho_by_response = rho_by_response,
                min_segment = as.integer(min_segment), years = yrs, sdir = sdir, ctx = ctx)
  descs <- if (parallel_stage_enabled(cfg, "response") && length(tasks) > 1L) {
    psock_task_lapply(cfg, tasks, "run_response_scenario_task", state, stage = "response", log_file = ctx$log)
  } else {
    lapply(tasks, run_response_scenario_task, state = state)
  }
  names(descs) <- responses
  summaries <- lapply(descs, function(z) z$summary)
  sm <- do.call(rbind, summaries); sm$scenario_id <- sid; sm$family <- family
  write_csv(sm, file.path(sdir, "scenario_summary.csv"))
  write_json(list(scenario_id = sid, family = family, min_segment = min_segment, rho_by_response = rho_by_response,
                  fdr = fdr_label, alpha = cfg$multiple_testing$alpha, supf_mc_B_by_response = as.list(supf_mc_B_by_response(cfg)), input_fingerprint = inventory$fingerprint,
                  config_hash = ctx$config_hash, code_hash = ctx$code_hash,
                  parallel = list(response_workers = resolve_worker_count(cfg, length(tasks), "response"), nested_parallelism = parallel_settings(cfg)$nested_parallelism)),
             file.path(sdir, "scenario_metadata.json"))
  desc <- list(id = sid, dir = sdir, results = descs, summary = sm)
  if (isTRUE(return_descriptor)) desc else hydrate_scenario_descriptor(desc)
}
