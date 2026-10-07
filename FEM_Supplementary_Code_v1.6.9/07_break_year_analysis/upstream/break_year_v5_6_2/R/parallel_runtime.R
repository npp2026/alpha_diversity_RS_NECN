# Cross-platform PSOCK parallel runtime for v5.6.
#
# Design rules:
#   * outer-parallel / inner-serial by default (no nested worker pools)
#   * worker count never changes deterministic MC chunk or pixel-bootstrap RNG
#   * workers never write the same raster; raster writes remain single-writer
#   * base R parallel only (Windows/Linux/macOS compatible)

is_parallel_worker <- function() identical(Sys.getenv("V56_PARALLEL_WORKER", unset = "0"), "1")

parallel_settings <- function(cfg) {
  p <- cfg$runtime$parallel %||% list()
  list(
    enabled = isTRUE(p$enabled %||% FALSE),
    backend = tolower(as.character(p$backend %||% "psock")),
    workers = p$workers %||% cfg$runtime$cores %||% 1L,
    max_workers = as.integer(p$max_workers %||% cfg$runtime$cores %||% 1L),
    response = isTRUE(p$strategy$response %||% FALSE),
    scenario = isTRUE(p$strategy$scenario %||% FALSE),
    supf_mc = isTRUE(p$strategy$supf_mc %||% FALSE),
    pixel_bootstrap = isTRUE(p$strategy$pixel_bootstrap %||% FALSE),
    validation_simulation = isTRUE(p$strategy$validation_simulation %||% FALSE),
    nested_parallelism = isTRUE(p$nested_parallelism %||% FALSE),
    blas_threads_per_worker = as.integer(p$blas_threads_per_worker %||% 1L)
  )
}

parallel_stage_enabled <- function(cfg, stage) {
  p <- parallel_settings(cfg)
  if (!p$enabled) return(FALSE)
  if (is_parallel_worker() && !p$nested_parallelism) return(FALSE)
  isTRUE(p[[stage]] %||% FALSE)
}

resolve_worker_count <- function(cfg, n_tasks = Inf, stage = NULL) {
  p <- parallel_settings(cfg)
  if (!p$enabled) return(1L)
  if (is_parallel_worker() && !p$nested_parallelism) return(1L)
  if (!is.null(stage) && !isTRUE(p[[stage]] %||% FALSE)) return(1L)

  detected <- suppressWarnings(parallel::detectCores(logical = TRUE))
  if (!is.finite(detected) || detected < 1L) detected <- 1L
  auto_workers <- max(1L, as.integer(detected) - 1L)

  requested <- p$workers
  if (is.character(requested) && length(requested) == 1L && tolower(requested) == "auto") requested <- auto_workers
  requested <- suppressWarnings(as.integer(requested))
  if (!is.finite(requested) || requested < 1L) requested <- 1L

  cap <- suppressWarnings(as.integer(p$max_workers))
  if (!is.finite(cap) || cap < 1L) cap <- requested
  w <- min(requested, cap, auto_workers)
  if (is.finite(n_tasks)) w <- min(w, max(1L, as.integer(n_tasks)))
  max(1L, as.integer(w))
}

parallel_worker_init <- function(project_dir, state, blas_threads = 1L) {
  Sys.setenv(
    V56_PARALLEL_WORKER = "1",
    OMP_NUM_THREADS = as.character(max(1L, as.integer(blas_threads))),
    OPENBLAS_NUM_THREADS = as.character(max(1L, as.integer(blas_threads))),
    MKL_NUM_THREADS = as.character(max(1L, as.integer(blas_threads))),
    VECLIB_MAXIMUM_THREADS = as.character(max(1L, as.integer(blas_threads))),
    NUMEXPR_NUM_THREADS = as.character(max(1L, as.integer(blas_threads)))
  )
  source(file.path(project_dir, "R", "load_all.R"), local = .GlobalEnv)
  load_v56_modules(project_dir)
  assign(".V56_PARALLEL_STATE", state, envir = .GlobalEnv)
  TRUE
}

# Run globally named task_fun(task, state) across a PSOCK cluster. The state is
# copied once per worker, not once per task. Task functions must return ordinary
# serializable R objects (never live terra external pointers).
psock_task_lapply <- function(cfg, tasks, task_fun, state = list(), stage, log_file = NULL) {
  if (!length(tasks)) return(list())
  workers <- resolve_worker_count(cfg, length(tasks), stage)
  fn <- get(task_fun, envir = .GlobalEnv, inherits = TRUE)
  if (workers <= 1L) return(lapply(tasks, fn, state = state))

  p <- parallel_settings(cfg)
  if (!identical(p$backend, "psock")) stop("v5.6 currently supports runtime.parallel.backend = psock only.", call. = FALSE)
  log_msg("Parallel stage '", stage, "': ", length(tasks), " task(s), ", workers, " PSOCK worker(s).", .file = log_file)
  cl <- parallel::makePSOCKcluster(workers, outfile = "")
  on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
  parallel::clusterCall(cl, parallel_worker_init, cfg$.project_dir, state, p$blas_threads_per_worker)
  parallel::parLapply(cl, tasks, function(task, task_fun) {
    fn <- get(task_fun, envir = .GlobalEnv, inherits = TRUE)
    fn(task, get(".V56_PARALLEL_STATE", envir = .GlobalEnv, inherits = FALSE))
  }, task_fun = task_fun)
}

parallel_plan_table <- function(cfg) {
  p <- parallel_settings(cfg)
  stages <- c("response", "scenario", "supf_mc", "pixel_bootstrap", "validation_simulation")
  data.frame(
    stage = stages,
    enabled = vapply(stages, function(s) isTRUE(p[[s]]), logical(1L)),
    resolved_worker_cap = vapply(stages, function(s) resolve_worker_count(cfg, Inf, s), integer(1L)),
    backend = p$backend,
    nested_parallelism = p$nested_parallelism,
    stringsAsFactors = FALSE
  )
}
