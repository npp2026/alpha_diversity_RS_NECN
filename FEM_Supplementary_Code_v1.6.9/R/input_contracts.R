# Input contracts used before model fitting; no response-dependent filtering.
fem_match_plots <- function(d) {
  if (anyDuplicated(names(d))) stop("Duplicate input column names")
  req <- c("plot_id","plot_year")
  if (!all(req %in% names(d))) stop("Training CSV requires plot_id and plot_year")
  if (!is.numeric(d$plot_year)) stop("plot_year must be numeric")
  keep <- is.finite(d$plot_year) & d$plot_year %in% c(2008:2012,2014:2017)
  if ("image_year" %in% names(d)) {
    if (!is.numeric(d$image_year)) stop("image_year must be numeric")
    keep <- keep & is.finite(d$image_year) & d$image_year == d$plot_year
  }
  keep[is.na(keep)] <- FALSE
  if (anyNA(d$plot_id[keep]) || any(!nzchar(trimws(as.character(d$plot_id[keep])))) ||
      anyDuplicated(d$plot_id[keep])) stop("Need one matched row per nonempty independent plot_id")
  list(data=d[keep,,drop=FALSE],included=keep)
}
fem_finite_rows <- function(d,columns) {
  if (!all(columns %in% names(d))) stop("Missing input columns: ",paste(setdiff(columns,names(d)),collapse=", "))
  if (any(!vapply(d[,columns,drop=FALSE],is.numeric,logical(1)))) stop("Model variables must be numeric")
  x <- as.matrix(d[,columns,drop=FALSE])
  # GEE's historical numeric missing sentinel is never a valid predictor.
  rowSums(!is.finite(x) | x == -9999) == 0L
}
fem_validate_groups <- function(groups) {
  if (!setequal(names(groups),c("Env","Prod","Het","Temp")) || any(lengths(groups)==0L))
    stop("Every source group requires predictors")
  all_features <- unlist(groups,use.names=FALSE)
  forbidden <- c("plot_id","plot_year","image_year","year","aug_offset","Lon_Export","Lat_Export",
                 "Rich_tree","Shannon_wiener","Forest_age","PFT","AGB","Forest_type","random_fold","knndm_fold")
  if (any(all_features %in% forbidden)) stop("Metadata/response cannot be model predictors: ",
    paste(intersect(all_features,forbidden),collapse=", "))
  duplicated_features <- unique(all_features[duplicated(all_features)])
  if (length(duplicated_features)) stop("Predictor groups overlap: ",paste(duplicated_features,collapse=", "))
  invisible(TRUE)
}
