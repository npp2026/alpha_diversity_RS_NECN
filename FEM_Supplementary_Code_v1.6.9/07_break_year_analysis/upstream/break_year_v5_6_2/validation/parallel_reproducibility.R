run_parallel_reproducibility_v56 <- function(project_dir = getwd()) {
  project_dir <- normalizePath(project_dir, winslash = "/", mustWork = TRUE)
  source(file.path(project_dir, "R", "load_all.R"), local = .GlobalEnv)
  load_v56_modules(project_dir)

  years <- 2001:2020
  B <- 240L; chunk <- 40L; seed <- 76001L
  serial <- simulate_supf_null(years, 5L, 5L, 0.30, B, seed, chunk, "L'Ecuyer-CMRG")

  cfgp <- list(
    .project_dir = project_dir,
    runtime = list(cores = 2L, parallel = list(
      enabled = TRUE, backend = "psock", workers = 2L, max_workers = 2L,
      nested_parallelism = FALSE, blas_threads_per_worker = 1L,
      strategy = list(response = FALSE, scenario = FALSE, supf_mc = TRUE,
                      pixel_bootstrap = FALSE, validation_simulation = FALSE)
    ))
  )
  resolved <- resolve_worker_count(cfgp, ceiling(B/chunk), "supf_mc")
  parallel_result <- simulate_supf_null(years, 5L, 5L, 0.30, B, seed, chunk,
                                        "L'Ecuyer-CMRG", cfg = cfgp)
  stopifnot(length(serial) == length(parallel_result),
            isTRUE(all.equal(serial, parallel_result, tolerance = 0, check.attributes = FALSE)))

  # Scheduling order must not alter the sorted null because every chunk owns a
  # deterministic chunk-key seed.
  tasks <- lapply(seq_len(6L), function(i) list(chunk_id = i, n = 40L))
  state <- list(years = as.numeric(years), min_left = 5L, min_right = 5L,
                rho = 0.30, base_seed = seed, rng_kind = "L'Ecuyer-CMRG")
  forward <- sort(as.numeric(unlist(lapply(tasks, supf_mc_chunk_task, state = state), use.names = FALSE)))
  reverse_order <- sort(as.numeric(unlist(lapply(rev(tasks), supf_mc_chunk_task, state = state), use.names = FALSE)))
  stopifnot(isTRUE(all.equal(forward, reverse_order, tolerance = 0, check.attributes = FALSE)))

  # Pixel bootstrap RNG is keyed by response + cell_id, so cell scheduling and
  # batch partition must not alter any cell's bootstrap result.
  cfgb <- list(
    rng = list(kind = "L'Ecuyer-CMRG", master_seed = 99123L),
    bootstrap = list(pixel = list(B = 60L, min_valid_replicates_for_ci = 5L)),
    analysis = list(primary = list(min_segment = 5L), slope_tol = list(SR = 0), tie_rule = "earliest")
  )
  x <- years - mean(years); tx <- 2010 - mean(years)
  ys <- list(
    `101` = 5 - 0.25*x + 0.65*pmax(0, x-tx) + sin(seq_along(x))*0.015,
    `202` = 6 - 0.20*x + 0.55*pmax(0, x-tx) + cos(seq_along(x))*0.020,
    `303` = 4 - 0.18*x + 0.50*pmax(0, x-tx) + sin(seq_along(x)*0.7)*0.018
  )
  eval_cells <- function(ids) {
    z <- lapply(ids, function(id) {
      y <- ys[[as.character(id)]]
      fit <- fit_break_core(y, years, 5L, 5L, 0, "earliest")
      bootstrap_one_pixel(y, years, fit, serial, 0.30, 0.10, cfgb, "SR", id)
    })
    out <- do.call(rbind, z); out[order(out$cell), , drop = FALSE]
  }
  p1 <- eval_cells(c(101L, 202L, 303L)); p2 <- eval_cells(c(303L, 101L, 202L))
  stopifnot(identical(names(p1), names(p2)))
  for (nm in names(p1)) {
    if (is.numeric(p1[[nm]])) stopifnot(isTRUE(all.equal(p1[[nm]], p2[[nm]], tolerance = 0, check.attributes = FALSE)))
    else stopifnot(identical(as.character(p1[[nm]]), as.character(p2[[nm]])))
  }

  # Nested pools are disabled inside workers by default.
  old <- Sys.getenv("V56_PARALLEL_WORKER", unset = NA_character_)
  on.exit({ if (is.na(old)) Sys.unsetenv("V56_PARALLEL_WORKER") else Sys.setenv(V56_PARALLEL_WORKER = old) }, add = TRUE)
  Sys.setenv(V56_PARALLEL_WORKER = "1")
  stopifnot(resolve_worker_count(cfgp, 10L, "supf_mc") == 1L)

  data.frame(
    check = c("mc_serial_vs_psock", "mc_task_order_invariance", "pixel_cell_order_invariance", "nested_parallelism_guard"),
    pass = TRUE,
    detail = c(paste0("resolved PSOCK workers = ", resolved), "chunk-key RNG", "cell-key RNG", "worker resolves to one inner worker"),
    stringsAsFactors = FALSE
  )
}
