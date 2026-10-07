# ==============================================================================
# v9_parallel_utils.R  (optimised)
# ------------------------------------------------------------------------------
# Shared base-R parallel utilities for the independent RF_6class v9 pipeline.
# Uses PSOCK clusters so the same scripts work under Windows .bat and Linux/macOS.
# Set N_WORKERS in the shell to control parallelism. N_WORKERS=1 keeps serial mode.
#
# ==============================================================================

v9_resolve_workers <- function(n_tasks, requested = NULL, reserve_cores = 1L) {
  n_tasks <- as.integer(max(0L, n_tasks))
  if (n_tasks <= 1L) return(1L)
  if (is.null(requested)) {
    requested <- suppressWarnings(as.integer(Sys.getenv("N_WORKERS", unset = "1")))
  }
  if (length(requested) != 1L || !is.finite(requested) || is.na(requested) || requested < 1L) requested <- 1L
  cores <- suppressWarnings(parallel::detectCores(logical = TRUE))
  if (!is.finite(cores) || is.na(cores) || cores < 1L) cores <- 1L
  max_cores <- max(1L, cores - as.integer(reserve_cores))
  as.integer(max(1L, min(requested, max_cores, n_tasks)))
}

# ---- internal cluster helpers ------------------------------------------------
# One place that builds a PSOCK cluster, loads packages and exports objects.
.v9_make_cluster <- function(n_workers) {
  cl <- tryCatch(parallel::makeCluster(n_workers, type = "PSOCK", outfile = ""),
                 error = function(e) e)
  cl
}
.v9_cluster_load_packages <- function(cl, packages) {
  if (length(packages)) {
    parallel::clusterCall(cl, function(pkgs) {
      invisible(lapply(pkgs, function(p) suppressPackageStartupMessages(library(p, character.only = TRUE))))
    }, unique(packages))
  }
  invisible(NULL)
}
.v9_cluster_export <- function(cl, export, export_env) {
  if (length(export)) parallel::clusterExport(cl, unique(export), envir = export_env)
  invisible(NULL)
}
.v9_cluster_map <- function(cl, X, FUN, dots, load_balance, n_workers, chunk_factor = 4L) {
  if (isTRUE(load_balance)) {
    # parLapplyLB hands out work dynamically; a modest chunk size keeps dispatch
    # overhead low while still balancing uneven tasks across workers.
    chunk_size <- max(1L, as.integer(ceiling(length(X) / (n_workers * chunk_factor))))
    do.call(parallel::parLapplyLB, c(list(cl = cl, X = X, fun = FUN, chunk.size = chunk_size), dots))
  } else {
    do.call(parallel::parLapply, c(list(cl = cl, X = X, fun = FUN), dots))
  }
}

# ---- reusable cluster pool ---------------------------------------------------
# A single per-session cluster, keyed by worker count. Reused across calls so a
# loop of v9_cluster_lapply(..., reuse_pool = TRUE) pays makeCluster() only once.
.v9_pool <- new.env(parent = emptyenv())
.v9_pool$cl <- NULL
.v9_pool$n_workers <- NA_integer_
.v9_pool$packages <- character(0)

v9_pool_cluster <- function(n_workers, packages = character()) {
  n_workers <- as.integer(n_workers)
  cl <- .v9_pool$cl
  # (Re)create when there is no live cluster or the requested width changed.
  if (is.null(cl) || !identical(.v9_pool$n_workers, n_workers)) {
    if (!is.null(cl)) try(parallel::stopCluster(cl), silent = TRUE)
    cl <- .v9_make_cluster(n_workers)
    if (inherits(cl, "error")) { .v9_pool$cl <- NULL; return(cl) }  # caller falls back
    .v9_pool$cl <- cl
    .v9_pool$n_workers <- n_workers
    .v9_pool$packages <- character(0)
  }
  # Load only packages not already loaded on this pooled cluster (library is idempotent).
  new_pkgs <- setdiff(unique(packages), .v9_pool$packages)
  if (length(new_pkgs)) {
    .v9_cluster_load_packages(cl, new_pkgs)
    .v9_pool$packages <- union(.v9_pool$packages, new_pkgs)
  }
  cl
}

v9_close_pool <- function() {
  if (!is.null(.v9_pool$cl)) {
    try(parallel::stopCluster(.v9_pool$cl), silent = TRUE)
  }
  .v9_pool$cl <- NULL
  .v9_pool$n_workers <- NA_integer_
  .v9_pool$packages <- character(0)
  invisible(NULL)
}
# Best-effort cleanup if the script forgets to close the pool explicitly.
reg.finalizer(.v9_pool, function(e) { if (!is.null(e$cl)) try(parallel::stopCluster(e$cl), silent = TRUE) }, onexit = TRUE)

v9_cluster_lapply <- function(X, FUN, ..., workers = NULL, packages = character(),
                              export = character(), export_env = parent.frame(),
                              seed = NULL, task_label = "task",
                              reuse_pool = FALSE, load_balance = FALSE) {
  n_tasks <- length(X)
  n_workers <- v9_resolve_workers(n_tasks, workers)
  dots <- list(...)
  if (n_workers <= 1L) {
    # Serial path also honours `seed`, so serial and parallel runs agree.
    if (!is.null(seed) && length(seed) == 1L && is.finite(seed)) set.seed(as.integer(seed))
    return(lapply(X, FUN, ...))
  }
  message(sprintf("[parallel] %s: %d tasks on %d PSOCK workers%s%s", task_label, n_tasks, n_workers,
                  if (reuse_pool) " (pooled)" else "", if (load_balance) " (LB)" else ""))

  if (isTRUE(reuse_pool)) {
    cl <- v9_pool_cluster(n_workers, packages)
    if (inherits(cl, "error")) {
      warning("Could not start/reuse pooled PSOCK cluster; falling back to serial. Reason: ",
              conditionMessage(cl), call. = FALSE)
      if (!is.null(seed) && length(seed) == 1L && is.finite(seed)) set.seed(as.integer(seed))
      return(lapply(X, FUN, ...))
    }
    # NOTE: no on.exit(stopCluster) -- the pool owns the cluster's lifetime.
  } else {
    cl <- .v9_make_cluster(n_workers)
    if (inherits(cl, "error")) {
      warning("Could not start PSOCK cluster; falling back to serial. Reason: ",
              conditionMessage(cl), call. = FALSE)
      if (!is.null(seed) && length(seed) == 1L && is.finite(seed)) set.seed(as.integer(seed))
      return(lapply(X, FUN, ...))
    }
    on.exit(parallel::stopCluster(cl), add = TRUE)
    .v9_cluster_load_packages(cl, packages)
  }

  # Export every call (cheap relative to worker startup; objects like `base`/`prep`
  # change between calls, so a pooled cluster must be refreshed each time).
  .v9_cluster_export(cl, export, export_env)
  if (!is.null(seed) && length(seed) == 1L && is.finite(seed)) {
    parallel::clusterSetRNGStream(cl, iseed = as.integer(seed))
  }
  .v9_cluster_map(cl, X, FUN, dots, load_balance, n_workers)
}

# Split seq_len(n) into ~ workers*chunk_factor contiguous chunks (order-preserving).
v9_chunk_indices <- function(n, workers = NULL, chunk_factor = 4L) {
  n <- as.integer(n)
  if (n <= 0L) return(list(integer(0)))
  w <- v9_resolve_workers(n, workers)
  n_chunks <- min(n, max(1L, as.integer(w * chunk_factor)))
  unname(split(seq_len(n), cut(seq_len(n), breaks = n_chunks, labels = FALSE)))
}

v9_log_parallel_settings <- function(n_tasks = NA_integer_, requested = NULL, label = "parallel setting") {
  if (is.null(requested)) requested <- suppressWarnings(as.integer(Sys.getenv("N_WORKERS", unset = "1")))
  nw <- if (is.finite(n_tasks)) v9_resolve_workers(n_tasks, requested) else requested
  message(sprintf("[parallel] %s: requested N_WORKERS=%s; effective workers=%s", label, requested, nw))
  invisible(nw)
}
