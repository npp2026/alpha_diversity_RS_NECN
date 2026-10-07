# Paired scenario comparisons require the same independent plots in the same order.
fem_validate_oof <- function(all_results) {
  scenarios <- c("S0_Env","S1a_Prod","S1b_Het","S1c_Temp","S2a_Prod_Het",
                 "S2b_Prod_Temp","S2c_Het_Temp","S3_RS_Full","S4_Integrated")
  if (!is.list(all_results) || !length(all_results) || is.null(names(all_results)) ||
      anyDuplicated(names(all_results))) stop("Invalid or empty named all_results object")
  for (key in names(all_results)) {
    m <- all_results[[key]]$model_results
    if (!all(scenarios %in% names(m))) stop(key, ": missing required scenarios")
    base <- m[[scenarios[1]]]; n <- length(base$observed)
    if (n < 2L) stop(key, ": no observations")
    ids <- base$plot_id
    if (!is.null(ids) && (length(ids)!=n || anyNA(ids) || anyDuplicated(ids) ||
                         any(!nzchar(trimws(as.character(ids)))))) stop(key, ": invalid plot IDs")
    same <- function(a,b) isTRUE(all.equal(unname(a),unname(b),check.attributes=FALSE,tolerance=0))
    fields <- c("plot_id","fold_assignment","fold_indices","coords","coordinates",
                "block_id_50km","spatial_block_50km","spatial_block_id_anonymized","block_id")
    if (is.null(ids) && is.null(base$coords) && is.null(base$coordinates))
      stop(key, ": cannot verify pairing without plot IDs or coordinates")
    for (sc in scenarios) {
      z <- m[[sc]]
      if (!isTRUE(z$success) || !is.numeric(z$observed) || !is.numeric(z$predictions) ||
          length(z$observed)!=n || length(z$predictions)!=n ||
          any(!is.finite(z$observed)) || any(!is.finite(z$predictions)))
        stop(key, "/", sc, ": unsuccessful, missing or non-finite OOF values")
      if (!same(base$observed,z$observed)) stop(key, "/", sc, ": observed order mismatch")
      for (field in fields) {
        a<-base[[field]]; b<-z[[field]]
        if (!is.null(a) || !is.null(b)) {
          if (is.null(a) || is.null(b) || !same(a,b)) stop(key, "/", sc, ": ", field, " order mismatch")
          if (if(is.null(dim(a))) length(a)!=n else nrow(a)!=n) stop(key, ": invalid ",field," length")
          if (anyNA(a)) stop(key, ": missing ",field)
        }
      }
    }
  }
  invisible(TRUE)
}

# Figure2, source contrasts and their bootstrap must use uncalibrated RF OOF.
# Original v1.6.9 stores QM in predictions and RF OOF in raw_predictions.
# Select raw explicitly on read; never guess the meaning of an unlabeled vector.
# This returns an in-memory copy and does not rewrite the input all_results.rds.
fem_use_raw_rf_oof <- function(all_results) {
  if (!is.list(all_results) || !length(all_results) || is.null(names(all_results)) ||
      anyDuplicated(names(all_results))) stop("Invalid or empty named all_results object")
  for (key in names(all_results)) {
    models <- all_results[[key]]$model_results
    if (!is.list(models) || !length(models)) stop(key, ": missing model_results")
    for (sc in names(models)) {
      z <- models[[sc]]
      raw <- z$raw_predictions
      if (is.null(raw)) {
        if (!identical(z$prediction_type, "RF_raw"))
          stop(key, "/", sc, ": raw RF OOF predictions are required. Supply raw_predictions ",
               "or predictions explicitly labeled prediction_type='RF_raw'; QM/unlabeled fallback is disabled.")
        raw <- z$predictions
      }
      if (!is.numeric(raw) || length(raw) != length(z$observed) || any(!is.finite(raw)))
        stop(key, "/", sc, ": missing, invalid or non-finite raw RF OOF predictions")
      z$predictions <- raw
      z$prediction_type <- "RF_raw"
      models[[sc]] <- z
    }
    all_results[[key]]$model_results <- models
  }
  fem_validate_oof(all_results)
  all_results
}
