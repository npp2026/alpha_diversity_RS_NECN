run_sensitivity_scenario_task <- function(task, state) {
  run_scenario(state$cfg, state$inventory, task$rho_by_response, task$min_segment,
               task$family, state$root, state$ctx, return_descriptor = TRUE)
}

run_boundary_sensitivity <- function(cfg, inventory, rho_by_response, ctx, primary) {
  if (!isTRUE(cfg$sensitivity$boundary$enabled)) return(list())
  vals <- unique(as.integer(unlist(cfg$sensitivity$boundary$min_segments)))
  root <- ensure_dir(file.path(ctx$root, "sensitivity", "boundary"))
  primary_m <- as.integer(cfg$analysis$primary$min_segment)
  out <- list()
  if (primary_m %in% vals) out[[paste0("L", primary_m)]] <- primary

  todo <- vals[vals != primary_m]
  tasks <- lapply(todo, function(m) list(label = paste0("L", m), min_segment = m,
                                         rho_by_response = rho_by_response, family = "BOUND"))
  if (length(tasks)) {
    state <- list(cfg = cfg, inventory = inventory, root = root, ctx = ctx)
    descs <- if (parallel_stage_enabled(cfg, "scenario") && length(tasks) > 1L) {
      psock_task_lapply(cfg, tasks, "run_sensitivity_scenario_task", state, stage = "scenario", log_file = ctx$log)
    } else lapply(tasks, run_sensitivity_scenario_task, state = state)
    for (i in seq_along(tasks)) out[[tasks[[i]]$label]] <- hydrate_scenario_descriptor(descs[[i]])
  }
  # Preserve configured window order in tables and downstream figures.
  out <- out[paste0("L", vals)]
  tab <- do.call(rbind, lapply(out, function(x) x$summary))
  write_csv(tab, file.path(root, "boundary_sensitivity_summary.csv"))
  out
}

rho_sensitivity_sources <- function(empirical, ctx) {
  responses <- names(empirical)
  sources <- list(empirical = empirical)
  comp_file <- file.path(ctx$root, "rho", "rho_model_comparison.csv")
  if (!file.exists(comp_file)) return(sources)
  comp <- utils::read.csv(comp_file, stringsAsFactors = FALSE, check.names = FALSE)
  if (!"response" %in% names(comp)) return(sources)
  for (model in c("linear", "hinge")) {
    col <- paste0("bias_corrected_rho.", model)
    if (!col %in% names(comp)) next
    idx <- match(responses, comp$response)
    z <- suppressWarnings(as.numeric(comp[[col]][idx]))
    if (length(z) == length(responses) && all(is.finite(z)) && all(abs(z) < 1)) {
      sources[[paste0(model, "_diagnostic")]] <- stats::setNames(as.list(z), responses)
    }
  }
  sources
}

parse_rho_sensitivity <- function(values, empirical, sources = list(empirical = empirical)) {
  if (is.null(sources$empirical)) sources$empirical <- empirical
  out <- list()
  for (v in values) {
    if (is.character(v)) {
      key <- tolower(trimws(v))
      if (key %in% names(sources)) {
        out[[key]] <- sources[[key]]
        next
      }
    }
    z <- suppressWarnings(as.numeric(v))
    if (!is.finite(z) || abs(z) >= 1) {
      stop("Invalid rho sensitivity value/token: ", v,
           ". Supported data-driven tokens: ", paste(names(sources), collapse = ", "), call. = FALSE)
    }
    out[[paste0("rho", format(z, trim = TRUE))]] <- stats::setNames(as.list(rep(z, length(empirical))), names(empirical))
  }
  out[!duplicated(vapply(out, function(x) paste(round(unlist(x), 8), collapse = "|"), character(1L)))]
}

write_rho_sensitivity_values <- function(vals, root) {
  rows <- do.call(rbind, lapply(names(vals), function(label) {
    z <- unlist(vals[[label]])
    data.frame(rho_sensitivity_label = label, response = names(z), rho = as.numeric(z),
               source = if (grepl("_diagnostic$", label)) "rho_model_comparison" else if (label == "empirical") "primary_empirical" else "fixed_numeric",
               stringsAsFactors = FALSE)
  }))
  write_csv(rows, file.path(root, "rho_sensitivity_values.csv"))
  invisible(rows)
}

run_rho_sensitivity <- function(cfg, inventory, empirical_rho, ctx, primary) {
  if (!isTRUE(cfg$sensitivity$rho$enabled)) return(list())
  root <- ensure_dir(file.path(ctx$root, "sensitivity", "rho"))
  sources <- rho_sensitivity_sources(empirical_rho, ctx)
  vals <- parse_rho_sensitivity(unlist(cfg$sensitivity$rho$values), empirical_rho, sources)
  write_rho_sensitivity_values(vals, root)
  out <- list(); m <- as.integer(cfg$analysis$primary$min_segment)
  emp_token <- paste(round(unlist(empirical_rho), 8), collapse = "|")
  tasks <- list()
  for (nm in names(vals)) {
    token <- paste(round(unlist(vals[[nm]]), 8), collapse = "|")
    if (identical(token, emp_token)) out[[nm]] <- primary
    else tasks[[length(tasks) + 1L]] <- list(label = nm, min_segment = m, rho_by_response = vals[[nm]], family = "RHO")
  }
  if (length(tasks)) {
    state <- list(cfg = cfg, inventory = inventory, root = root, ctx = ctx)
    descs <- if (parallel_stage_enabled(cfg, "scenario") && length(tasks) > 1L) {
      psock_task_lapply(cfg, tasks, "run_sensitivity_scenario_task", state, stage = "scenario", log_file = ctx$log)
    } else lapply(tasks, run_sensitivity_scenario_task, state = state)
    for (i in seq_along(tasks)) out[[tasks[[i]]$label]] <- hydrate_scenario_descriptor(descs[[i]])
  }
  out <- out[names(vals)]
  tab <- do.call(rbind, lapply(names(out), function(nm) transform(out[[nm]]$summary, rho_sensitivity_label = nm)))
  write_csv(tab, file.path(root, "rho_sensitivity_summary.csv"))
  out
}

build_multiplicity_summary <- function(cfg, primary, ctx) {
  if (!isTRUE(cfg$sensitivity$multiplicity$enabled)) return(NULL)
  root <- ensure_dir(file.path(ctx$root, "sensitivity", "multiplicity"))
  alpha <- as.numeric(cfg$multiple_testing$alpha)
  rows <- list()
  for (response in names(primary$results)) {
    res <- primary$results[[response]]
    for (method in toupper(unlist(cfg$sensitivity$multiplicity$methods))) {
      q <- res$q[[method]]
      if (is.null(q)) next
      sig <- terra::ifel(q <= alpha, 1, 0)
      rec <- terra::ifel(sig == 1 & res$layers$direction_ok == 1, 1, 0)
      n <- as.numeric(terra::global(!is.na(res$layers$p_mc), "sum", na.rm = TRUE)[[1L, 1L]])
      nr <- raster_true_count(rec)
      rows[[length(rows) + 1L]] <- data.frame(response = response, method = method, eligible_pixels = n, recovery_breaks = nr,
                                               recovery_pct_eligible = if (n > 0) 100 * nr/n else NA_real_)
    }
  }
  tab <- do.call(rbind, rows); write_csv(tab, file.path(root, "multiplicity_sensitivity_summary.csv")); tab
}
