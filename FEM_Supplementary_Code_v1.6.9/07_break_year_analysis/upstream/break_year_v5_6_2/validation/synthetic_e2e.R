synthetic_extract_cells <- function(r, cells) {
  if (!length(cells)) return(numeric())
  as.numeric(as.matrix(terra::extract(r, as.integer(cells), raw = TRUE))[, 1L])
}

synthetic_metric_row <- function(response, check, value, threshold = NA_real_, relation = "", severity = "hard", pass = TRUE, note = "") {
  data.frame(response = response, check = check, value = as.numeric(value), threshold = as.numeric(threshold),
             relation = relation, severity = severity, pass = isTRUE(pass), note = note, stringsAsFactors = FALSE)
}

evaluate_synthetic_run_v56 <- function(res, project_dir = getwd(), fail_on_hard = TRUE) {
  assert_packages(c("terra", "yaml"))
  truth_file <- file.path(project_dir, "validation", "synthetic_data", "synthetic_truth.csv")
  spec_file <- file.path(project_dir, "validation", "synthetic_spec.yml")
  if (!file.exists(truth_file) || !file.exists(spec_file)) stop("Synthetic truth/spec are missing. Run generate_synthetic_dataset() first.", call. = FALSE)
  truth <- utils::read.csv(truth_file, stringsAsFactors = FALSE); spec <- yaml::read_yaml(spec_file)
  thr <- res$cfg$validation$synthetic_e2e %||% list()
  recall_min <- as.numeric(thr$strong_recovery_recall_min %||% 0.60)
  null_max <- as.numeric(thr$null_recovery_fpr_max %||% 0.05)
  reverse_max <- as.numeric(thr$reverse_recovery_rate_max %||% 0.05)
  rho_warn <- as.numeric(thr$rho_abs_error_warning %||% 0.25)
  tau_tol <- as.numeric(thr$boundary_break_tolerance_years %||% 1)
  rows <- list(); add <- function(x) rows[[length(rows) + 1L]] <<- x
  mask_cells <- function(name) truth$cell[truth$name == name]
  strong_names <- c("recovery_2008", "recovery_2012", "boundary_recovery_2005", "boundary_recovery_2015", "exact_hinge_2010")

  for (response in names(res$primary$results)) {
    pr <- res$primary$results[[response]]
    exact_linear <- mask_cells("exact_linear"); exact_hinge <- mask_cells("exact_hinge_2010"); missing <- mask_cells("missing_2010")
    reverse <- mask_cells("reverse_2010"); null <- mask_cells("null_linear"); strong <- truth$cell[truth$name %in% strong_names]
    st_lin <- synthetic_extract_cells(pr$layers$fit_status, exact_linear); st_hinge <- synthetic_extract_cells(pr$layers$fit_status, exact_hinge)
    tau_hinge <- synthetic_extract_cells(pr$layers$break_year, exact_hinge); p_missing <- synthetic_extract_cells(pr$layers$p_mc, missing)
    rec_strong <- synthetic_extract_cells(pr$recovery, strong); rec_null <- synthetic_extract_cells(pr$recovery, null); rec_reverse <- synthetic_extract_cells(pr$recovery, reverse)
    recall <- mean(rec_strong == 1, na.rm = TRUE); null_fpr <- mean(rec_null == 1, na.rm = TRUE); reverse_rate <- mean(rec_reverse == 1, na.rm = TRUE)
    add(synthetic_metric_row(response, "exact_linear_fit_status", mean(st_lin == 2), 1, "==", "hard", all(st_lin == 2), "fit_status code 2 = perfect_linear"))
    add(synthetic_metric_row(response, "exact_hinge_fit_status", mean(st_hinge == 3), 1, "==", "hard", all(st_hinge == 3), "fit_status code 3 = exact_hinge"))
    add(synthetic_metric_row(response, "exact_hinge_break_year", max(abs(tau_hinge - 2010), na.rm = TRUE), 0, "==", "hard", all(tau_hinge == 2010), "noiseless exact hinge"))
    add(synthetic_metric_row(response, "missing_series_excluded", mean(!is.finite(p_missing)), 1, "==", "hard", all(!is.finite(p_missing)), "complete-series inference must exclude missing_2010"))
    add(synthetic_metric_row(response, "strong_recovery_recall", recall, recall_min, ">=", "hard", is.finite(recall) && recall >= recall_min))
    add(synthetic_metric_row(response, "null_recovery_fpr", null_fpr, null_max, "<=", "hard", is.finite(null_fpr) && null_fpr <= null_max))
    add(synthetic_metric_row(response, "reverse_recovery_rate", reverse_rate, reverse_max, "<=", "hard", is.finite(reverse_rate) && reverse_rate <= reverse_max,
                             "reverse breaks may be structurally significant but must fail decline-to-recovery direction"))

    for (cname in c("early_recovery_2003", "late_recovery_2017")) {
      true_tau <- if (cname == "early_recovery_2003") 2003 else 2017; cc <- mask_cells(cname)
      for (m in c(5L, 3L)) {
        sc <- res$boundary[[paste0("L", m)]]
        if (is.null(sc)) next
        z <- synthetic_extract_cells(sc$results[[response]]$layers$break_year, cc); med <- stats::median(z, na.rm = TRUE)
        add(synthetic_metric_row(response, paste0(cname, "_median_break_L", m), med, true_tau, "target", "diagnostic", TRUE,
                                 "raw break estimates; significance is intentionally not required for window-identifiability check"))
        if (m == 3L) add(synthetic_metric_row(response, paste0(cname, "_L3_absolute_error"), abs(med - true_tau), tau_tol, "<=", "hard", abs(med - true_tau) <= tau_tol))
      }
    }

    rrow <- res$rho[res$rho$response == response, , drop = FALSE]
    true_rho <- as.numeric(spec$responses[[response]]$rho)
    if (nrow(rrow)) {
      err <- abs(as.numeric(rrow$rho_for_calibration[[1L]]) - true_rho)
      add(synthetic_metric_row(response, "rho_absolute_error", err, rho_warn, "<= warning threshold", "warning", err <= rho_warn,
                               paste0("true synthetic rho = ", true_rho, "; mixed structural breaks intentionally stress linear-residual rho estimation")))
    }
  }

  # Workflow-level artifact checks catch integration failures outside the core math.
  expected_tables <- c("Table1A_primary_overall_summary.csv", "Table1B_primary_regional_summary.csv", "Table2A_boundary_sensitivity_overall.csv",
                       "Table2B_boundary_sensitivity_regional.csv", "Table3_regional_trajectory_uncertainty.csv", "Table4_spatial_clusters.csv", "Table5_rho_calibration.csv")
  have_tables <- file.exists(file.path(res$ctx$root, "tables", expected_tables))
  add(synthetic_metric_row("workflow", "publication_tables_present", mean(have_tables), 1, "==", "hard", all(have_tables), paste(expected_tables[!have_tables], collapse = ";")))
  parallel_plan <- file.path(res$ctx$root, "config", "parallel_plan.csv")
  add(synthetic_metric_row("workflow", "parallel_plan_present", as.numeric(file.exists(parallel_plan)), 1, "==", "hard", file.exists(parallel_plan),
                           "resolved response/scenario/MC/pixel parallel plan must be recorded"))
  spatial_csv <- list.files(file.path(res$ctx$root, "spatial"), pattern = "(cluster_summary|scenario_iou)\\.csv$", recursive = TRUE, full.names = TRUE)
  add(synthetic_metric_row("workflow", "spatial_outputs_present", length(spatial_csv), 2, ">=", "hard", length(spatial_csv) >= 2L))
  temporal_ok <- !is.null(res$spatial) && all(vapply(names(res$spatial), function(response) {
    d <- res$spatial[[response]]$scenario_iou
    !is.null(d) && nrow(d) > 0 && all(c("exact_year_agreement", "within1_year_agreement", "mean_abs_year_shift", "median_abs_year_shift") %in% names(d)) && any(is.finite(d$within1_year_agreement))
  }, logical(1L)))
  add(synthetic_metric_row("workflow", "spatial_temporal_agreement_present", as.numeric(temporal_ok), 1, "==", "hard", temporal_ok))

  tab <- do.call(rbind, rows)
  out_csv <- file.path(res$ctx$root, "validation", "synthetic_e2e_metrics.csv"); write_csv(tab, out_csv)
  hard_failed <- tab$severity == "hard" & !tab$pass
  warn_failed <- tab$severity == "warning" & !tab$pass
  summary <- list(total_checks = nrow(tab), hard_failures = sum(hard_failed), warnings = sum(warn_failed), metrics_file = out_csv,
                  status = if (any(hard_failed)) "FAIL" else if (any(warn_failed)) "PASS_WITH_WARNINGS" else "PASS")
  write_json(summary, file.path(res$ctx$root, "validation", "synthetic_e2e_summary.json"))
  lines <- c("# Synthetic E2E validation", "", paste0("- Status: **", summary$status, "**"),
             paste0("- Hard failures: ", summary$hard_failures), paste0("- Warnings: ", summary$warnings), "",
             "## Metrics", "", paste(capture.output(print(tab, row.names = FALSE)), collapse = "\n"), "",
             "The synthetic fixture contains null, decline-to-recovery, reverse-break, out-of-primary-window, exact-fit, missing-year, AR(1), and spatial-cluster cases.")
  writeLines(lines, file.path(res$ctx$root, "validation", "synthetic_e2e_report.md"), useBytes = TRUE)
  if (fail_on_hard && any(hard_failed)) stop("Synthetic E2E validation failed: ", paste(tab$check[hard_failed], collapse = ", "), call. = FALSE)
  invisible(list(metrics = tab, summary = summary))
}

run_synthetic_e2e_v56 <- function(project_dir = getwd(), reuse_data = FALSE, fail_on_hard = TRUE) {
  project_dir <- normalizePath(project_dir, winslash = "/", mustWork = TRUE)
  synthetic_root <- file.path(project_dir, "validation", "synthetic_data")
  if (!(isTRUE(reuse_data) && file.exists(file.path(synthetic_root, "synthetic_manifest.csv"))))
    generate_synthetic_dataset(project_dir, overwrite = TRUE)
  config_path <- file.path(project_dir, "config", "v5_6_synthetic.yml")
  res <- run_pipeline_until("validation", config_path = config_path, project_dir = project_dir)
  ev <- evaluate_synthetic_run_v56(res, project_dir, fail_on_hard = fail_on_hard)
  invisible(list(run = res, evaluation = ev))
}
