`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L) y else x

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", args, value = TRUE)
  if (length(f)) dirname(normalizePath(sub("^--file=", "", f[[1L]]), winslash = "/", mustWork = TRUE))
  else normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

is_abs_path <- function(x) grepl("^(/|[A-Za-z]:[/\\\\]|\\\\\\\\)", x)
resolve_path <- function(x, base = getwd()) {
  if (is.null(x) || !nzchar(x)) return("")
  if (!is_abs_path(x)) x <- file.path(base, x)
  normalizePath(x, winslash = "/", mustWork = FALSE)
}
ensure_dir <- function(x) { dir.create(x, recursive = TRUE, showWarnings = FALSE); x }

region_field_config <- function(cfg) {
  id <- as.character(cfg$data$region_id_field %||% "")
  label <- as.character(cfg$data$region_label_field %||% id)
  group <- as.character(cfg$data$region_group_field %||% "")
  group_label <- as.character(cfg$data$region_group_label_field %||% group)
  list(id = id, label = label, group = group, group_label = group_label)
}

region_metadata <- function(regions, cfg, require_unique_id = TRUE) {
  attrs <- as.data.frame(regions)
  f <- region_field_config(cfg)
  if (!nzchar(f$id) || !f$id %in% names(attrs))
    stop("region_id_field not found in shapefile: ", f$id, call. = FALSE)
  ids <- as.character(attrs[[f$id]])
  bad_id <- is.na(ids) | !nzchar(trimws(ids))
  if (any(bad_id)) stop("region_id_field contains missing/blank values: ", f$id, call. = FALSE)
  if (isTRUE(require_unique_id) && anyDuplicated(ids))
    stop("region_id_field must uniquely identify each polygon feature: ", f$id, call. = FALSE)

  take_optional <- function(field, fallback = NA_character_) {
    if (!nzchar(field) || !field %in% names(attrs)) return(rep(fallback, nrow(attrs)))
    as.character(attrs[[field]])
  }
  labels <- take_optional(f$label)
  labels[is.na(labels) | !nzchar(trimws(labels))] <- ids[is.na(labels) | !nzchar(trimws(labels))]
  groups <- take_optional(f$group)
  group_labels <- take_optional(f$group_label)
  if (nzchar(f$group) && nzchar(f$group_label)) {
    missing_gl <- is.na(group_labels) | !nzchar(trimws(group_labels))
    group_labels[missing_gl] <- groups[missing_gl]
  }
  data.frame(feature = seq_len(nrow(attrs)), region = ids, region_label = labels,
             region_group = groups, region_group_label = group_labels,
             stringsAsFactors = FALSE)
}

region_unique_character_candidates <- function(regions) {
  attrs <- as.data.frame(regions)
  rows <- lapply(names(attrs), function(nm) {
    x <- attrs[[nm]]
    is_char <- is.character(x) || is.factor(x)
    xc <- if (is_char) as.character(x) else rep(NA_character_, length(x))
    blank <- if (is_char) is.na(xc) | !nzchar(trimws(xc)) else rep(TRUE, length(x))
    data.frame(field = nm, class = class(x)[1L], n_rows = length(x),
               n_missing = sum(is.na(x)),
               n_blank = if (is_char) sum(!is.na(xc) & !nzchar(trimws(xc))) else NA_integer_,
               n_unique = length(unique(x[!is.na(x)])),
               unique_character_candidate = is_char && !any(blank) && length(unique(xc)) == length(xc),
               stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

assert_packages <- function(pkgs) {
  miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1L), quietly = TRUE)]
  if (length(miss)) stop("Missing R package(s): ", paste(miss, collapse = ", "), call. = FALSE)
  invisible(TRUE)
}

log_msg <- function(..., .file = NULL) {
  z <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste0(..., collapse = ""))
  cat(z, "\n")
  if (!is.null(.file)) cat(z, "\n", file = .file, append = TRUE)
  invisible(z)
}

hash_text <- function(x, algo = "xxhash64") {
  assert_packages("digest")
  digest::digest(paste(as.character(x), collapse = "\n"), algo = algo, serialize = FALSE)
}

seed_from_key <- function(master_seed, ...) {
  assert_packages("digest")
  key <- paste(master_seed, ..., sep = "|")
  h <- digest::digest(key, algo = "xxhash32", serialize = FALSE)
  s <- suppressWarnings(strtoi(substr(h, 1L, 7L), base = 16L))
  if (!is.finite(s) || s <= 0L) s <- abs(as.integer(master_seed)) %% 100000000L + 1L
  as.integer(s)
}

set_rng <- function(kind = "L'Ecuyer-CMRG", seed = 12345L) {
  RNGkind(kind); set.seed(as.integer(seed)); invisible(seed)
}

write_json <- function(x, path) {
  assert_packages("jsonlite"); ensure_dir(dirname(path))
  jsonlite::write_json(x, path, auto_unbox = TRUE, pretty = TRUE, null = "null", digits = 16)
  invisible(path)
}
write_csv <- function(x, path) { ensure_dir(dirname(path)); utils::write.csv(x, path, row.names = FALSE, na = ""); invisible(path) }

atomic_save_rds <- function(x, path) {
  ensure_dir(dirname(path))
  tmp <- tempfile(pattern = paste0(basename(path), ".tmp_"), tmpdir = dirname(path))
  on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)
  saveRDS(x, tmp)
  if (file.exists(path)) unlink(path)
  if (!file.rename(tmp, path)) stop("Atomic RDS cache rename failed: ", path, call. = FALSE)
  invisible(path)
}

safe_read_rds_cache <- function(path) {
  if (!file.exists(path)) return(NULL)
  tryCatch(readRDS(path), error = function(e) {
    warning("Ignoring unreadable cache and recomputing: ", path, " (", conditionMessage(e), ")", call. = FALSE)
    try(unlink(path), silent = TRUE)
    NULL
  })
}

file_manifest <- function(files, role = NULL) {
  files <- normalizePath(files, winslash = "/", mustWork = TRUE); fi <- file.info(files)
  out <- data.frame(path = files, size = as.numeric(fi$size),
                    mtime_utc = format(fi$mtime, "%Y-%m-%dT%H:%M:%OS6Z", tz = "UTC"),
                    md5 = unname(tools::md5sum(files)), stringsAsFactors = FALSE)
  if (!is.null(role)) out$role <- role
  out
}

project_code_manifest <- function(project_dir) {
  files <- c(list.files(file.path(project_dir, "R"), pattern = "\\.R$", recursive = TRUE, full.names = TRUE),
             list.files(project_dir, pattern = "^[0-9]{2}.*\\.R$", full.names = TRUE),
             list.files(file.path(project_dir, "validation"), pattern = "\\.(R|py)$", full.names = TRUE))
  files <- sort(unique(files[file.exists(files)]))
  file_manifest(files, role = paste0("code_", substring(normalizePath(files, winslash="/"), nchar(normalizePath(project_dir, winslash="/")) + 2L)))
}

modal_value <- function(x) {
  x <- x[is.finite(x)]; if (!length(x)) return(NA_real_)
  tb <- table(x); as.numeric(names(tb)[which.max(tb)])
}
safe_quantile_discrete <- function(x, probs) {
  x <- x[!is.na(x)]; if (!length(x)) return(rep(NA_real_, length(probs)))
  as.numeric(stats::quantile(x, probs = probs, type = 1, names = FALSE))
}

read_config <- function(path = Sys.getenv("BREAK_CONFIG", "config/v5_6.yml"), project_dir = getwd()) {
  # Retired modes must fail explicitly, including direct upstream launchers.
  for (flag in c("PUBLICATION_FREEZE", "MC_RESOLUTION_REMEDIATION")) {
    if (tolower(trimws(Sys.getenv(flag, unset = "false"))) %in% c("true", "1", "yes", "y")) {
      stop("Unsupported historical mode: ", flag, ". This package supports analysis and diagnostics only.", call. = FALSE)
    }
  }
  assert_packages(c("yaml", "digest"))
  path <- resolve_path(path, project_dir)
  if (!file.exists(path)) stop("Config not found: ", path, call. = FALSE)
  cfg <- yaml::read_yaml(path)
  cfg$.config_path <- path
  cfg$.project_dir <- normalizePath(project_dir, winslash = "/", mustWork = TRUE)

  data_root <- Sys.getenv("DATA_ROOT", unset = "")
  if (!nzchar(data_root)) data_root <- cfg$.project_dir
  cfg$.data_root <- resolve_path(data_root, cfg$.project_dir)
  for (nm in names(cfg$data$responses)) cfg$data$responses[[nm]]$folder <- resolve_path(cfg$data$responses[[nm]]$folder, cfg$.data_root)
  cfg$data$regions_shapefile <- resolve_path(cfg$data$regions_shapefile, cfg$.data_root)

  out_override <- Sys.getenv("OUTPUT_ROOT", unset = "")
  if (nzchar(out_override)) cfg$output$root <- out_override
  cfg$output$root <- resolve_path(cfg$output$root, cfg$.project_dir)

  cfg$supf$cache_dir <- resolve_path(cfg$supf$cache_dir, cfg$output$root)

  core_override <- suppressWarnings(as.integer(Sys.getenv("MAX_CORES", unset = "")))
  if (is.finite(core_override) && core_override > 0L) {
    cfg$runtime$cores <- core_override
    if (is.null(cfg$runtime$parallel)) cfg$runtime$parallel <- list()
    cfg$runtime$parallel$max_workers <- core_override
    # MAX_CORES is an operational cap, not a request to disable auto sizing.
  }

  # analysis.primary is the single source of truth for primary multiplicity settings.
  cfg$multiple_testing$primary_method <- toupper(cfg$analysis$primary$fdr)
  cfg$multiple_testing$alpha <- as.numeric(cfg$analysis$primary$alpha)
  cfg$.code_manifest <- project_code_manifest(cfg$.project_dir)
  cfg$.code_hash <- substr(hash_text(apply(cfg$.code_manifest[, c("path", "size", "md5")], 1L, paste, collapse = "|")), 1L, 16L)
  validate_config(cfg)
  cfg
}

validate_config <- function(cfg) {
  yrs <- as.integer(unlist(cfg$data$years))
  if (length(yrs) < 8L || anyDuplicated(yrs) || is.unsorted(yrs)) stop("data.years must be sorted unique years.", call. = FALSE)
  if (!length(cfg$data$responses)) stop("At least one response is required.", call. = FALSE)
  rf <- region_field_config(cfg)
  if (isTRUE(cfg$data$mask_to_regions) && !nzchar(rf$id)) stop("data.region_id_field must be configured when mask_to_regions=true.", call. = FALSE)
  m <- as.integer(cfg$analysis$primary$min_segment)
  if (m < 2L || 2L * m >= length(yrs)) stop("Primary min_segment is incompatible with series length.", call. = FALSE)
  if (!isTRUE(cfg$data$require_complete_series)) stop("v5.6 MC inference requires complete annual series; pattern-specific missing-year calibration is not implemented.", call. = FALSE)
  if (!identical(tolower(cfg$supf$p_rule), "greater_equal")) stop("v5.6 requires supf.p_rule = greater_equal.", call. = FALSE)
  base_B <- as.integer(cfg$supf$mc_B); if (!is.finite(base_B) || base_B < 1L) stop("supf.mc_B must be a positive integer.", call. = FALSE)
  by_B <- cfg$supf$mc_B_by_response %||% NULL
  if (!is.null(by_B)) {
    if (is.null(names(by_B)) || any(!names(by_B) %in% names(cfg$data$responses))) stop("supf.mc_B_by_response names must be configured responses.", call. = FALSE)
    bz <- suppressWarnings(as.integer(unlist(by_B)))
    if (length(bz) != length(by_B) || any(!is.finite(bz) | bz < 1L)) stop("supf.mc_B_by_response values must be positive integers.", call. = FALSE)
  }
  if (!tolower(cfg$multiple_testing$family) %in% "per_response") stop("Only per_response FDR family is supported.", call. = FALSE)
  if (!toupper(cfg$multiple_testing$primary_method) %in% c("BH", "BY")) stop("Primary FDR must be BH or BY.", call. = FALSE)
  a <- as.numeric(cfg$multiple_testing$alpha); if (!is.finite(a) || a <= 0 || a >= 1) stop("Primary alpha must be in (0,1).", call. = FALSE)
  rr <- as.numeric(cfg$rho_estimation$max_abs_rho)
  gmin <- as.numeric(cfg$rho_estimation$bias_correction$grid_min); gmax <- as.numeric(cfg$rho_estimation$bias_correction$grid_max)
  if (!is.finite(rr) || rr <= 0 || rr >= 1) stop("rho_estimation.max_abs_rho must be in (0,1).", call. = FALSE)
  if (gmin >= 0 || gmax <= 0 || max(abs(c(gmin, gmax))) + 1e-12 < rr) stop("rho bias grid must span +/- max_abs_rho.", call. = FALSE)
  if (!tolower(cfg$analysis$primary$direction) %in% "decline_recovery") stop("Only decline_recovery direction is implemented.", call. = FALSE)
  gp <- cfg$validation$rho_stress_gate %||% NULL
  if (!is.null(gp)) {
    rrmin <- as.numeric(gp$robust_retention_min %||% NA_real_)
    crmax <- as.numeric(gp$collapse_retention_max %||% NA_real_)
    mcm <- as.numeric(gp$mc_floor_clear_multiplier %||% NA_real_)
    ysw <- as.numeric(gp$modal_year_shift_warn %||% NA_real_)
    if (!is.finite(rrmin) || rrmin <= 0 || rrmin > 1) stop("validation.rho_stress_gate.robust_retention_min must be in (0,1].", call. = FALSE)
    if (!is.finite(crmax) || crmax < 0 || crmax >= rrmin) stop("validation.rho_stress_gate.collapse_retention_max must be >=0 and below robust_retention_min.", call. = FALSE)
    if (!is.finite(mcm) || mcm < 1) stop("validation.rho_stress_gate.mc_floor_clear_multiplier must be >=1.", call. = FALSE)
    if (!is.finite(ysw) || ysw < 0) stop("validation.rho_stress_gate.modal_year_shift_warn must be >=0.", call. = FALSE)
  }
  sm <- as.integer(unlist(cfg$sensitivity$boundary$min_segments %||% m))
  if (any(sm < 2L | 2L * sm >= length(yrs))) stop("At least one boundary sensitivity min_segment is incompatible with series length.", call. = FALSE)
  if (!tolower(cfg$bootstrap$pixel$target %||% "primary_recovery_pixels") %in% "primary_recovery_pixels") stop("Only bootstrap.pixel.target=primary_recovery_pixels is implemented.", call. = FALSE)
  p <- cfg$runtime$parallel %||% list(enabled = FALSE)
  if (isTRUE(p$enabled)) {
    backend <- tolower(as.character(p$backend %||% "psock"))
    if (!backend %in% "psock") stop("runtime.parallel.backend must be psock.", call. = FALSE)
    wreq <- p$workers %||% cfg$runtime$cores %||% 1L
    if (is.character(wreq)) {
      if (length(wreq) != 1L || tolower(wreq) != "auto") stop("runtime.parallel.workers must be 'auto' or a positive integer.", call. = FALSE)
    } else {
      wi <- suppressWarnings(as.integer(wreq))
      if (!is.finite(wi) || wi < 1L) stop("runtime.parallel.workers must be 'auto' or a positive integer.", call. = FALSE)
    }
    mw <- suppressWarnings(as.integer(p$max_workers %||% cfg$runtime$cores %||% 1L))
    if (!is.finite(mw) || mw < 1L) stop("runtime.parallel.max_workers must be >= 1.", call. = FALSE)
    bt <- suppressWarnings(as.integer(p$blas_threads_per_worker %||% 1L))
    if (!is.finite(bt) || bt < 1L) stop("runtime.parallel.blas_threads_per_worker must be >= 1.", call. = FALSE)
    if (isTRUE(p$nested_parallelism)) warning("Nested parallelism is enabled. This is not recommended for raster/MC workloads and may oversubscribe CPU/RAM.", call. = FALSE)
  }
  invisible(TRUE)
}

supf_mc_B_for_response <- function(cfg, response = NULL) {
  base_B <- as.integer(cfg$supf$mc_B)
  if (!is.finite(base_B) || base_B < 1L) stop("supf.mc_B must be a positive integer.", call. = FALSE)
  by_response <- cfg$supf$mc_B_by_response %||% NULL
  if (!is.null(response) && !is.null(by_response) && response %in% names(by_response)) {
    z <- as.integer(by_response[[response]])
    if (!is.finite(z) || z < 1L) stop("Invalid response-specific SupF MC B for ", response, ".", call. = FALSE)
    return(z)
  }
  base_B
}

supf_mc_B_by_response <- function(cfg) {
  responses <- names(cfg$data$responses)
  stats::setNames(vapply(responses, function(response) supf_mc_B_for_response(cfg, response), integer(1L)), responses)
}

supf_mc_B_description <- function(cfg) {
  z <- supf_mc_B_by_response(cfg)
  if (length(unique(z)) == 1L) return(as.character(z[[1L]]))
  paste(paste0(names(z), "=", z), collapse = "; ")
}

config_hash <- function(cfg) {
  clean <- cfg[!grepl("^\\.", names(cfg))]
  hash_text(paste(capture.output(str(clean, give.attr = FALSE)), collapse = "\n"))
}

scenario_id <- function(family, min_segment, rho_by_response, fdr = "BH") {
  rho_token <- paste(vapply(names(rho_by_response), function(nm) paste0(nm, sprintf("%.3f", as.numeric(rho_by_response[[nm]]))), character(1L)), collapse = "_")
  gsub("[^A-Za-z0-9_.-]", "-", paste(family, paste0("L", min_segment, "R", min_segment), rho_token, toupper(fdr), sep = "_"))
}

create_run_context <- function(cfg) {
  ch <- substr(config_hash(cfg), 1L, 10L)
  run_id <- paste0("run_", format(Sys.time(), "%Y%m%d_%H%M%S"), "_", ch)
  run_root <- ensure_dir(file.path(cfg$output$root, run_id))
  dirs <- c("config", "rho", "calibration", "primary", "sensitivity/boundary", "sensitivity/rho", "sensitivity/multiplicity",
            "uncertainty/pixel", "uncertainty/trajectory", "spatial", "tables", "validation", "logs")
  invisible(vapply(dirs, function(d) ensure_dir(file.path(run_root, d)), character(1L)))
  ctx <- list(run_id = run_id, root = run_root, log = file.path(run_root, "logs", "pipeline.log"), config_hash = ch, code_hash = cfg$.code_hash)
  if (exists("parallel_plan_table", mode = "function")) write_csv(parallel_plan_table(cfg), file.path(run_root, "config", "parallel_plan.csv"))
  ctx
}
