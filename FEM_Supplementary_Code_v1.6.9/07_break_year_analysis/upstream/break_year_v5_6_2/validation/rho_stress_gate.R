# Targeted real-data rho stress gate for an already completed v5.6 run.
#
# This does NOT modify the source run and does NOT recompute primary, boundary,
# bootstrap, trajectory, spatial, tables, or figures. It reuses the source
# run's resolved input configuration and exact corrected linear-residual rho
# diagnostics, then runs one response-specific SupF/FDR scenario in a new,
# separately provenance-tracked output directory.
#
# The gate is intentionally asymmetric: the linear-residual rho remains a
# diagnostic stress model, not an alternative primary estimator. Large loss of
# recovery abundance therefore produces CONDITIONAL_GO (mandatory methods and
# interpretation caveat), not an automatic estimator switch. NO_GO is reserved
# for an uninterpretable gate (invalid/missing inputs, inconsistent responses,
# or other hard integrity failures).

rho_stress_gate_required_files <- function(source_run) {
  c(
    config = file.path(source_run, "config", "config_resolved.yml"),
    manifest = file.path(source_run, "config", "run_manifest.json"),
    rho_comparison = file.path(source_run, "rho", "rho_model_comparison.csv"),
    primary_summary = file.path(source_run, "tables", "Table1A_primary_overall_summary.csv")
  )
}

rho_stress_gate_validate_source <- function(source_run) {
  if (!nzchar(source_run)) {
    stop("Set RHO_STRESS_SOURCE_RUN to a completed v5.6 run directory.", call. = FALSE)
  }
  source_run <- normalizePath(source_run, winslash = "/", mustWork = TRUE)
  req <- rho_stress_gate_required_files(source_run)
  missing <- req[!file.exists(req)]
  if (length(missing)) stop("Rho stress source run is missing required file(s): ", paste(missing, collapse = "; "), call. = FALSE)
  source_run
}

rho_stress_gate_linear_rho <- function(source_run, responses) {
  f <- rho_stress_gate_required_files(source_run)[["rho_comparison"]]
  d <- utils::read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
  col <- "bias_corrected_rho.linear"
  if (!all(c("response", col) %in% names(d))) stop("rho_model_comparison.csv lacks response-specific linear corrected rho.", call. = FALSE)
  if (anyDuplicated(as.character(d$response))) stop("rho_model_comparison.csv has duplicate response rows.", call. = FALSE)
  idx <- match(responses, d$response)
  if (anyNA(idx)) stop("rho_model_comparison.csv is missing one or more configured responses: ", paste(responses[is.na(idx)], collapse = ", "), call. = FALSE)
  z <- suppressWarnings(as.numeric(d[[col]][idx]))
  if (any(!is.finite(z)) || any(abs(z) >= 1)) stop("Invalid/missing linear diagnostic rho for one or more responses.", call. = FALSE)
  stats::setNames(as.list(z), responses)
}

rho_stress_gate_policy <- function(project_dir) {
  # These are audit/decision thresholds, not inferential parameters. They are
  # loaded from the current project config so an old source run can be reviewed
  # without changing any of its statistical settings.
  defaults <- list(
    robust_retention_min = 0.50,
    collapse_retention_max = 0.10,
    mc_floor_clear_multiplier = 10,
    modal_year_shift_warn = 1
  )
  f <- file.path(project_dir, "config", "v5_6.yml")
  if (!file.exists(f)) return(defaults)
  x <- tryCatch(yaml::read_yaml(f), error = function(e) NULL)
  p <- if (!is.null(x)) x$validation$rho_stress_gate else NULL
  if (is.null(p)) return(defaults)
  for (nm in names(defaults)) {
    z <- suppressWarnings(as.numeric(p[[nm]]))
    if (length(z) == 1L && is.finite(z)) defaults[[nm]] <- z
  }
  if (defaults$robust_retention_min <= 0 || defaults$robust_retention_min > 1) stop("rho_stress_gate robust_retention_min must be in (0,1].", call. = FALSE)
  if (defaults$collapse_retention_max < 0 || defaults$collapse_retention_max >= defaults$robust_retention_min) stop("rho_stress_gate collapse_retention_max must be >=0 and below robust_retention_min.", call. = FALSE)
  if (defaults$mc_floor_clear_multiplier < 1) stop("rho_stress_gate mc_floor_clear_multiplier must be >=1.", call. = FALSE)
  if (defaults$modal_year_shift_warn < 0) stop("rho_stress_gate modal_year_shift_warn must be >=0.", call. = FALSE)
  defaults
}

rho_stress_gate_compare <- function(primary_summary, stress_summary) {
  needed <- c("response", "rho", "recovery_breaks", "recovery_pct_eligible", "modal_break_year")
  if (!all(needed %in% names(primary_summary))) stop("Primary summary lacks required rho stress comparison columns.", call. = FALSE)
  stress_needed <- c(needed, "significant_breaks", "mc_min_p", "min_rank_at_mc_floor", "frozen_fdr_p_cutoff")
  if (!all(stress_needed %in% names(stress_summary))) stop("Stress summary lacks required rho stress comparison columns.", call. = FALSE)
  pr <- as.character(primary_summary$response); sr <- as.character(stress_summary$response)
  if (anyDuplicated(pr) || anyDuplicated(sr)) stop("Primary/stress summaries must contain one row per response.", call. = FALSE)
  if (!setequal(pr, sr)) stop("Primary and stress summaries contain different response sets.", call. = FALSE)
  rows <- lapply(pr, function(response) {
    p <- primary_summary[primary_summary$response == response, , drop = FALSE][1L, ]
    s <- stress_summary[stress_summary$response == response, , drop = FALSE][1L, ]
    retention <- if (is.finite(p$recovery_breaks) && p$recovery_breaks > 0) s$recovery_breaks / p$recovery_breaks else NA_real_
    cutoff_ratio <- if (is.finite(s$frozen_fdr_p_cutoff) && is.finite(s$mc_min_p) && s$mc_min_p > 0) s$frozen_fdr_p_cutoff / s$mc_min_p else NA_real_
    data.frame(
      response = response,
      primary_rho = as.numeric(p$rho),
      linear_diagnostic_rho = as.numeric(s$rho),
      rho_linear_minus_primary = as.numeric(s$rho) - as.numeric(p$rho),
      rho_abs_difference = abs(as.numeric(s$rho) - as.numeric(p$rho)),
      primary_recovery_breaks = as.numeric(p$recovery_breaks),
      stress_recovery_breaks = as.numeric(s$recovery_breaks),
      recovery_retention_vs_primary = retention,
      primary_recovery_pct_eligible = as.numeric(p$recovery_pct_eligible),
      stress_recovery_pct_eligible = as.numeric(s$recovery_pct_eligible),
      stress_significant_breaks = as.numeric(s$significant_breaks),
      primary_modal_break_year = as.numeric(p$modal_break_year),
      stress_modal_break_year = as.numeric(s$modal_break_year),
      modal_break_year_shift = if (is.finite(p$modal_break_year) && is.finite(s$modal_break_year)) as.numeric(s$modal_break_year) - as.numeric(p$modal_break_year) else NA_real_,
      stress_mc_min_p = as.numeric(s$mc_min_p),
      stress_min_rank_at_mc_floor = as.numeric(s$min_rank_at_mc_floor),
      stress_fdr_p_cutoff = as.numeric(s$frozen_fdr_p_cutoff),
      stress_fdr_cutoff_to_mc_floor = cutoff_ratio,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

rho_stress_gate_add_floor_counts <- function(cmp, stress, cfg) {
  method <- toupper(as.character(cfg$multiple_testing$primary_method))
  counts <- vapply(as.character(cmp$response), function(response) {
    z <- stress$results[[response]]$fdr[[method]]$counts
    if (is.null(z) || !length(z)) return(NA_real_)
    as.numeric(z[[1L]])
  }, numeric(1L))
  cmp$stress_mc_floor_count <- counts
  cmp$stress_mc_floor_count_to_required_rank <- ifelse(
    is.finite(cmp$stress_min_rank_at_mc_floor) & cmp$stress_min_rank_at_mc_floor > 0,
    cmp$stress_mc_floor_count / cmp$stress_min_rank_at_mc_floor,
    NA_real_
  )
  cmp
}

rho_stress_gate_classify <- function(cmp, policy, context = c("subset", "full")) {
  context <- match.arg(context)
  invalid <- !is.finite(cmp$primary_rho) | !is.finite(cmp$linear_diagnostic_rho) |
    !is.finite(cmp$primary_recovery_breaks) | cmp$primary_recovery_breaks < 0 |
    !is.finite(cmp$stress_recovery_breaks) | cmp$stress_recovery_breaks < 0 |
    (cmp$primary_recovery_breaks > 0 & !is.finite(cmp$recovery_retention_vs_primary))

  retention_class <- rep("NOT_EVALUABLE", nrow(cmp))
  ok <- !invalid & is.finite(cmp$recovery_retention_vs_primary)
  retention_class[ok & cmp$recovery_retention_vs_primary >= policy$robust_retention_min] <- "ROBUST"
  retention_class[ok & cmp$recovery_retention_vs_primary < policy$robust_retention_min & cmp$recovery_retention_vs_primary >= policy$collapse_retention_max] <- "SENSITIVE"
  retention_class[ok & cmp$recovery_retention_vs_primary < policy$collapse_retention_max] <- "COLLAPSE"
  # If primary had zero recovery, abundance retention is undefined rather than a failure.
  retention_class[!invalid & cmp$primary_recovery_breaks == 0] <- "PRIMARY_ZERO_RECOVERY"

  mc_resolution_class <- rep("NO_FDR_CUTOFF", nrow(cmp))
  mc_ok <- is.finite(cmp$stress_fdr_cutoff_to_mc_floor)
  mc_resolution_class[mc_ok & cmp$stress_fdr_cutoff_to_mc_floor >= policy$mc_floor_clear_multiplier] <- "CLEAR_OF_MC_FLOOR"
  mc_resolution_class[mc_ok & cmp$stress_fdr_cutoff_to_mc_floor < policy$mc_floor_clear_multiplier] <- "NEAR_MC_FLOOR"

  year_class <- rep("NOT_EVALUABLE", nrow(cmp))
  yr_ok <- is.finite(cmp$modal_break_year_shift)
  year_class[yr_ok & abs(cmp$modal_break_year_shift) <= policy$modal_year_shift_warn] <- "MODAL_YEAR_STABLE"
  year_class[yr_ok & abs(cmp$modal_break_year_shift) > policy$modal_year_shift_warn] <- "MODAL_YEAR_SHIFT"

  cmp$retention_class <- retention_class
  cmp$mc_resolution_class <- mc_resolution_class
  cmp$modal_year_class <- year_class
  cmp$row_gate <- ifelse(invalid, "NO_GO", ifelse(retention_class %in% c("ROBUST", "PRIMARY_ZERO_RECOVERY"), "GO", "CONDITIONAL_GO"))

  overall <- if (any(invalid)) {
    "NO_GO"
  } else if (all(retention_class %in% c("ROBUST", "PRIMARY_ZERO_RECOVERY"))) {
    "GO"
  } else {
    "CONDITIONAL_GO"
  }

  severity <- if (any(retention_class == "COLLAPSE")) "SEVERE_RHO_MODEL_SENSITIVITY" else if (any(retention_class == "SENSITIVE")) "MATERIAL_RHO_MODEL_SENSITIVITY" else "RHO_STRESS_ROBUST"
  if (overall == "NO_GO") severity <- "GATE_INTEGRITY_FAILURE"

  conditions <- character()
  if (overall == "CONDITIONAL_GO") {
    conditions <- c(conditions,
      "Keep hinge-residual bias-corrected rho as the specified primary estimator; do not switch estimators post hoc.",
      "Report that detection abundance is sensitive to the temporal-dependence nuisance model, with exact linear-diagnostic stress results in validation/supplementary material.",
      if (context == "subset") "Interpret subset recovery percentages/counts only as methodological audit diagnostics, never as paper Results." else "Treat the full-data linear-diagnostic counts as sensitivity evidence; formal study inference remains the specified hinge-residual primary analysis.")
    if (any(retention_class == "COLLAPSE")) {
      conditions <- c(conditions,
        "Because at least one response nearly collapses under the linear-diagnostic rho stress, explicitly justify why linear residual autocorrelation is signal-contaminated and why the prevalidated hinge-residual estimator is primary.")
    }
  }
  if (any(mc_resolution_class == "NEAR_MC_FLOOR")) {
    conditions <- c(conditions,
      if (context == "subset") "The stress scenario has an FDR cutoff near the finite-MC floor; re-check MC resolution on the full-data primary run before reporting full-data results." else "The full-data linear-diagnostic stress scenario has an FDR cutoff near the finite-MC floor; report its finite-MC limitation separately from the primary MC-resolution audit.")
  }
  if (any(year_class == "MODAL_YEAR_SHIFT")) {
    conditions <- c(conditions,
      if (context == "subset") {
        "In the source real-data subset stress audit, the stress scenario shifts the modal break year beyond the audit tolerance; treat this as subset-only methodological sensitivity evidence and avoid model-independent exact-year claims."
      } else {
        "In the full-data exact linear-diagnostic rho stress, at least one response shifts the modal break year beyond the audit tolerance; avoid temporal-dependence-model-independent exact-year ecological claims unless supported by boundary/bootstrap uncertainty."
      })
  }
  if (overall == "GO" && !length(conditions)) {
    conditions <- "No additional rho-stress condition beyond the existing v5.6 interpretation guards."
  }
  if (overall == "NO_GO") {
    conditions <- "Do not use the diagnostic decision until the gate integrity failure is resolved and the targeted gate is rerun."
  }

  list(
    comparison = cmp,
    overall_gate = overall,
    severity = severity,
    conditions = unique(conditions),
    policy = policy
  )
}

rho_stress_gate_decision_table <- function(classification) {
  cmp <- classification$comparison
  data.frame(
    gate = "exact_linear_diagnostic_rho",
    overall_gate = classification$overall_gate,
    severity = classification$severity,
    robust_retention_min = classification$policy$robust_retention_min,
    collapse_retention_max = classification$policy$collapse_retention_max,
    mc_floor_clear_multiplier = classification$policy$mc_floor_clear_multiplier,
    modal_year_shift_warn = classification$policy$modal_year_shift_warn,
    responses = paste(cmp$response, collapse = ";"),
    min_recovery_retention = if (any(is.finite(cmp$recovery_retention_vs_primary))) min(cmp$recovery_retention_vs_primary, na.rm = TRUE) else NA_real_,
    max_rho_abs_difference = if (any(is.finite(cmp$rho_abs_difference))) max(cmp$rho_abs_difference, na.rm = TRUE) else NA_real_,
    any_near_mc_floor = any(cmp$mc_resolution_class == "NEAR_MC_FLOOR"),
    any_modal_year_shift = any(cmp$modal_year_class == "MODAL_YEAR_SHIFT"),
    conditions = paste(classification$conditions, collapse = " | "),
    stringsAsFactors = FALSE
  )
}

rho_stress_gate_write_decision <- function(classification, validation_dir) {
  cmp <- classification$comparison
  decision <- rho_stress_gate_decision_table(classification)
  write_csv(cmp, file.path(validation_dir, "rho_stress_gate_summary.csv"))
  write_csv(decision, file.path(validation_dir, "rho_stress_gate_decision.csv"))
  write_json(list(
    gate = "exact_linear_diagnostic_rho",
    overall_gate = classification$overall_gate,
    severity = classification$severity,
    policy = classification$policy,
    conditions = classification$conditions,
    by_response = lapply(seq_len(nrow(cmp)), function(i) as.list(cmp[i, , drop = FALSE]))
  ), file.path(validation_dir, "rho_stress_gate_decision.json"))
  invisible(decision)
}

rho_stress_gate_report_lines <- function(source_run, source_fingerprint, ctx_root, scenario_id, classification) {
  cmp <- classification$comparison
  c(
    "# v5.6 exact linear-diagnostic rho stress gate", "",
    paste0("- Source run: `", source_run, "`"),
    paste0("- Source input fingerprint: `", source_fingerprint, "`"),
    paste0("- Gate output: `", ctx_root, "`"),
    paste0("- Scenario: `", scenario_id, "`"),
    paste0("- Automated diagnostic gate: **", classification$overall_gate, "**"),
    paste0("- Rho sensitivity severity: **", classification$severity, "**"), "",
    "## Interpretation guard", "",
    "This is a response-specific stress test at the corrected linear-residual rho values from the completed source run.",
    "It does not replace the prevalidated hinge-residual empirical primary estimator.",
    "Because the source run is a real-data subset, all FDR recovery counts/percentages below remain methodological audit diagnostics and are not paper Results.", "",
    "## Primary versus exact linear-diagnostic rho", "",
    "```text",
    paste(capture.output(print(cmp, row.names = FALSE)), collapse = "\n"),
    "```", "",
    "## Diagnostic decision", "",
    paste0("**", classification$overall_gate, " — ", classification$severity, "**"), "",
    paste0("Policy: retention >= ", format(classification$policy$robust_retention_min, trim = TRUE),
           " is ROBUST; retention < ", format(classification$policy$collapse_retention_max, trim = TRUE),
           " is COLLAPSE; intermediate values are SENSITIVE. These are audit decision thresholds, not statistical significance thresholds."), "",
    "Conditions:",
    paste0("- ", classification$conditions), "",
    "## Gate use", "",
    "These labels summarize this diagnostic run; they do not certify a release or authorize publication.",
    "Use recovery_retention_vs_primary together with the pre-existing rho-model disagreement warning to document how strongly inferential abundance depends on the nuisance temporal-dependence model before full-data analysis.",
    "GO means the exact linear-diagnostic stress is comparatively robust. CONDITIONAL_GO means full-data production may proceed only with the listed, predeclared rho-model interpretation conditions. NO_GO is reserved for an uninterpretable/integrity-failed gate and requires repair plus rerun."
  )
}

review_rho_stress_gate_v56 <- function(project_dir = getwd(), gate_run = Sys.getenv("RHO_STRESS_GATE_RUN", unset = "")) {
  assert_packages(c("yaml", "jsonlite"))
  if (!nzchar(gate_run)) stop("Set RHO_STRESS_GATE_RUN to a completed rho stress gate run directory.", call. = FALSE)
  gate_run <- normalizePath(gate_run, winslash = "/", mustWork = TRUE)
  validation_dir <- file.path(gate_run, "validation")
  summary_file <- file.path(validation_dir, "rho_stress_gate_summary.csv")
  source_file <- file.path(gate_run, "config", "rho_stress_gate_source.json")
  if (!file.exists(summary_file) || !file.exists(source_file)) stop("Gate run lacks rho_stress_gate_summary.csv or rho_stress_gate_source.json.", call. = FALSE)
  cmp <- utils::read.csv(summary_file, stringsAsFactors = FALSE, check.names = FALSE)
  # Re-review may be applied to an older gate summary. Recompute only fields that
  # do not require raster results; if new required fields are absent, fail rather
  # than fabricate a decision.
  required <- c("response", "primary_rho", "linear_diagnostic_rho", "primary_recovery_breaks", "stress_recovery_breaks", "recovery_retention_vs_primary")
  if (!all(required %in% names(cmp))) stop("Existing gate summary is too old/incomplete for automated review; rerun 10_run_rho_stress_gate.R.", call. = FALSE)
  if (!"rho_abs_difference" %in% names(cmp)) cmp$rho_abs_difference <- abs(cmp$linear_diagnostic_rho - cmp$primary_rho)
  if (!"stress_fdr_cutoff_to_mc_floor" %in% names(cmp)) cmp$stress_fdr_cutoff_to_mc_floor <- NA_real_
  if (!"modal_break_year_shift" %in% names(cmp)) {
    if (all(c("primary_modal_break_year", "stress_modal_break_year") %in% names(cmp)))
      cmp$modal_break_year_shift <- cmp$stress_modal_break_year - cmp$primary_modal_break_year
    else cmp$modal_break_year_shift <- NA_real_
  }
  policy <- rho_stress_gate_policy(project_dir)
  classification <- rho_stress_gate_classify(cmp, policy)
  rho_stress_gate_write_decision(classification, validation_dir)
  source_info <- jsonlite::fromJSON(source_file, simplifyVector = TRUE)
  lines <- rho_stress_gate_report_lines(
    source_run = source_info$source_run %||% "unknown",
    source_fingerprint = source_info$source_input_fingerprint %||% "unknown",
    ctx_root = gate_run,
    scenario_id = source_info$scenario %||% "RHO_LINEAR_DIAGNOSTIC",
    classification = classification
  )
  report <- file.path(validation_dir, "RHO_STRESS_GATE_REPORT.md")
  writeLines(lines, report, useBytes = TRUE)
  invisible(list(gate_run = gate_run, classification = classification, report = report))
}

run_rho_stress_gate_v56 <- function(project_dir = getwd(), source_run = Sys.getenv("RHO_STRESS_SOURCE_RUN", unset = "")) {
  assert_packages(c("yaml", "terra", "digest", "jsonlite"))
  source_run <- rho_stress_gate_validate_source(source_run)
  req <- rho_stress_gate_required_files(source_run)

  cfg <- read_config(req[["config"]], project_dir)
  source_manifest <- jsonlite::fromJSON(req[["manifest"]], simplifyVector = TRUE)
  inventory <- build_input_inventory(cfg)
  if (!is.null(source_manifest$input_fingerprint) && nzchar(source_manifest$input_fingerprint) &&
      !identical(as.character(source_manifest$input_fingerprint), as.character(inventory$fingerprint))) {
    stop("Current subset input fingerprint does not match source run. Refusing a cross-input rho stress comparison.", call. = FALSE)
  }

  rho_linear <- rho_stress_gate_linear_rho(source_run, names(inventory$files))
  primary_summary <- utils::read.csv(req[["primary_summary"]], stringsAsFactors = FALSE, check.names = FALSE)

  # Keep the source statistical configuration, but isolate new outputs and cache
  # so the completed source run is never mutated by a later code version.
  cfg$output$root <- ensure_dir(file.path(project_dir, "rho_stress_gate_outputs_v5_6"))
  cfg$supf$cache_dir <- ensure_dir(file.path(cfg$output$root, "_cache", "supf"))
  cfg$runtime$overwrite <- FALSE
  ctx <- create_run_context(cfg)
  write_resolved_config(cfg, ctx)

  policy <- rho_stress_gate_policy(project_dir)
  source_meta <- list(
    gate = "exact_linear_diagnostic_rho",
    source_run = source_run,
    source_run_id = source_manifest$run_id %||% basename(source_run),
    source_config_hash = source_manifest$config_hash %||% NA_character_,
    source_code_hash = source_manifest$code_hash %||% NA_character_,
    source_input_fingerprint = source_manifest$input_fingerprint %||% inventory$fingerprint,
    gate_code_hash = ctx$code_hash,
    rho_by_response = rho_linear,
    gate_policy = policy,
    interpretation = "Scientific stress test only; hinge-residual corrected rho remains primary. Subset FDR family is audit-only."
  )
  source_meta_file <- file.path(ctx$root, "config", "rho_stress_gate_source.json")
  write_json(source_meta, source_meta_file)

  m <- as.integer(cfg$analysis$primary$min_segment)
  stress <- run_scenario(cfg, inventory, rho_linear, m, "RHO_LINEAR_DIAGNOSTIC",
                         file.path(ctx$root, "sensitivity", "rho"), ctx)
  source_meta$scenario <- stress$id
  write_json(source_meta, source_meta_file)
  cmp <- rho_stress_gate_compare(primary_summary, stress$summary)
  cmp <- rho_stress_gate_add_floor_counts(cmp, stress, cfg)
  classification <- rho_stress_gate_classify(cmp, policy)
  rho_stress_gate_write_decision(classification, file.path(ctx$root, "validation"))

  source_fingerprint <- source_manifest$input_fingerprint %||% inventory$fingerprint
  lines <- rho_stress_gate_report_lines(source_run, source_fingerprint, ctx$root, stress$id, classification)
  report <- file.path(ctx$root, "validation", "RHO_STRESS_GATE_REPORT.md")
  writeLines(lines, report, useBytes = TRUE)
  write_json(list(
    status = "COMPLETED_SCIENTIFIC_STRESS_TEST",
    diagnostic_gate = classification$overall_gate,
    severity = classification$severity,
    source_run = source_run,
    scenario = stress$id,
    summary_csv = file.path(ctx$root, "validation", "rho_stress_gate_summary.csv"),
    decision_csv = file.path(ctx$root, "validation", "rho_stress_gate_decision.csv"),
    decision_json = file.path(ctx$root, "validation", "rho_stress_gate_decision.json"),
    report = report
  ), file.path(ctx$root, "validation", "rho_stress_gate_status.json"))
  write_run_manifest(cfg, inventory, rho_tab = NULL, primary = NULL, ctx = ctx,
                     status = "completed", stage = paste0("rho_stress_gate_", tolower(classification$overall_gate)))
  log_msg("Exact linear-diagnostic rho stress gate completed: ", ctx$root,
          "; diagnostic gate = ", classification$overall_gate,
          "; severity = ", classification$severity, .file = ctx$log)
  invisible(list(cfg = cfg, ctx = ctx, source_run = source_run, rho = rho_linear,
                 stress = stress, comparison = classification$comparison,
                 classification = classification, report = report))
}
