# ==============================================================================
# Shared post-processing utilities for RF_6class state-trajectory.
# - strict/safe environment parsing
# - canonical six-class/focal-four definitions
# - common-support paired estimands
# - support-aware, vectorised spatial block bootstrap
# - reusable block layouts for RF/ratio sensitivity scans
# - target transition/agreement summaries
# ==============================================================================
options(stringsAsFactors = FALSE)

V12_CLASS6_LEVELS <- c('H+','H0','H-','L+','L0','L-')
V12_FOCAL4_LEVELS <- c('H+','H0','L+','L0')
V12_REG_ORDER <- c('I','II','III','IV','V','outside','Overall')
V12_METRIC_LEVELS <- c('Richness','Shannon')
V12_CLASS6_CODE <- setNames(seq_along(V12_CLASS6_LEVELS), V12_CLASS6_LEVELS)

v12_env_str <- function(name, default = '') {
  v <- Sys.getenv(name, unset = NA_character_)
  if (is.na(v) || !nzchar(trimws(v))) default else trimws(v)
}

v12_env_bool <- function(name, default = FALSE) {
  raw <- v12_env_str(name, '')
  if (!nzchar(raw)) return(isTRUE(default))
  v <- tolower(raw)
  if (v %in% c('true','t','1','yes','y')) return(TRUE)
  if (v %in% c('false','f','0','no','n')) return(FALSE)
  stop('Invalid logical environment value ', name, '=', raw, call. = FALSE)
}

v12_env_num <- function(name, default, lower = -Inf, upper = Inf,
                        lower_inclusive = TRUE, upper_inclusive = TRUE) {
  raw <- v12_env_str(name, '')
  if (!nzchar(raw)) x <- as.numeric(default) else x <- suppressWarnings(as.numeric(raw))
  ok <- length(x) == 1L && is.finite(x)
  if (ok) ok <- if (lower_inclusive) x >= lower else x > lower
  if (ok) ok <- if (upper_inclusive) x <= upper else x < upper
  if (!ok) stop('Invalid numeric environment value ', name, '=', ifelse(nzchar(raw), raw, '<default>'), call. = FALSE)
  x
}

v12_env_int <- function(name, default, lower = -Inf, upper = Inf) {
  x <- v12_env_num(name, default, lower, upper)
  xi <- suppressWarnings(as.integer(round(x)))
  if (is.na(xi) || !isTRUE(all.equal(as.numeric(xi), as.numeric(x), tolerance = 1e-9))) {
    stop('Environment value ', name, ' must be a representable integer; got ', x, call. = FALSE)
  }
  xi
}

# Validate an environment expression before evaluating it.  Production settings only
# need literals, c(...), integer ranges (a:b), parentheses, and unary +/-; arbitrary
# function calls such as system(...) are rejected.
.v12_safe_expr <- function(x) {
  if (is.atomic(x) || is.null(x)) return(TRUE)
  if (is.name(x)) return(as.character(x) %in% c('TRUE','FALSE','NA','NULL','Inf','NaN','NA_integer_','NA_real_','NA_character_'))
  if (!is.call(x)) return(FALSE)
  fn <- as.character(x[[1]])
  if (!fn %in% c('c',':','(', '+','-')) return(FALSE)
  all(vapply(as.list(x)[-1], .v12_safe_expr, logical(1)))
}

v12_env_expr <- function(name, default) {
  raw <- v12_env_str(name, '')
  if (!nzchar(raw)) return(default)
  parsed <- try(parse(text = raw, keep.source = FALSE), silent = TRUE)
  if (inherits(parsed, 'try-error') || length(parsed) != 1L || !.v12_safe_expr(parsed[[1]])) {
    stop('Invalid or unsafe vector expression in ', name, ': ', raw,
         '. Allowed forms include 2001:2020, c(50,75,100), and c("H+","H0").', call. = FALSE)
  }
  safe <- new.env(parent = emptyenv())
  safe$c <- base::c; safe$`:` <- base::`:`; safe$`(` <- base::`(`
  safe$`+` <- base::`+`; safe$`-` <- base::`-`
  safe$`TRUE` <- TRUE; safe$`FALSE` <- FALSE; safe$`NA` <- NA
  safe$`NULL` <- NULL; safe$`Inf` <- Inf; safe$`NaN` <- NaN
  safe$`NA_integer_` <- NA_integer_; safe$`NA_real_` <- NA_real_; safe$`NA_character_` <- NA_character_
  out <- try(eval(parsed[[1]], envir = safe), silent = TRUE)
  if (inherits(out, 'try-error')) stop('Could not evaluate ', name, ': ', raw, call. = FALSE)
  out
}

v12_make_abs_path <- function(path, base_dir) {
  path <- gsub('\\\\', '/', path)
  if (grepl('^[A-Za-z]:/|^/', path)) path else file.path(base_dir, path)
}

v12_msg <- function(...) cat(sprintf(...), '\n')

v12_validate_year_vector <- function(x, label = 'year vector', min_n = 1L) {
  y <- suppressWarnings(as.integer(x))
  if (length(y) < min_n || any(!is.finite(y)) || any(is.na(y)) || anyDuplicated(y)) {
    stop(label, ' must contain at least ', min_n, ' unique finite integer years.', call. = FALSE)
  }
  if (!identical(as.numeric(y), as.numeric(x))) stop(label, ' must contain integer years.', call. = FALSE)
  y
}

v12_zone_value <- function(z) {
  z <- trimws(as.character(z))
  z[is.na(z) | !nzchar(z)] <- 'outside'
  bad <- !z %in% V12_REG_ORDER[V12_REG_ORDER != 'Overall']
  z[bad] <- 'outside'
  z
}

v12_validate_cell_table <- function(d, require_rf = FALSE, path = '<cell table>') {
  need <- c('cell','x','y','rich_class6','shannon_class6')
  if (require_rf) need <- c(need, 'rich_RF','rich_trend_state','shannon_RF','shannon_trend_state')
  miss <- setdiff(need, names(d))
  if (length(miss)) stop('Missing required columns in ', path, ': ', paste(miss, collapse=', '), call. = FALSE)
  if (!nrow(d)) stop('Cell table is empty: ', path, call. = FALSE)
  d$cell <- suppressWarnings(as.integer(d$cell))
  if (any(is.na(d$cell)) || any(d$cell <= 0L) || anyDuplicated(d$cell)) stop('Invalid or duplicate cell IDs in ', path, call. = FALSE)
  d$x <- suppressWarnings(as.numeric(d$x)); d$y <- suppressWarnings(as.numeric(d$y))
  if (any(!is.finite(d$x)) || any(!is.finite(d$y))) stop('Non-finite x/y coordinates in ', path, call. = FALSE)
  bad_r <- unique(as.character(d$rich_class6[!is.na(d$rich_class6) & !d$rich_class6 %in% V12_CLASS6_LEVELS]))
  bad_s <- unique(as.character(d$shannon_class6[!is.na(d$shannon_class6) & !d$shannon_class6 %in% V12_CLASS6_LEVELS]))
  if (length(bad_r) || length(bad_s)) {
    stop('Non-canonical class labels in ', path, ': ', paste(unique(c(bad_r,bad_s)), collapse=', '), call. = FALSE)
  }
  if (require_rf) {
    for (nm in c('rich_RF','shannon_RF')) {
      d[[nm]] <- suppressWarnings(as.numeric(d[[nm]]))
      bad <- !is.na(d[[nm]]) & !is.finite(d[[nm]])
      if (any(bad)) stop('Non-finite non-NA values in ', nm, ' of ', path, call. = FALSE)
    }
    for (nm in c('rich_trend_state','shannon_trend_state')) {
      bad <- !is.na(d[[nm]]) & !as.character(d[[nm]]) %in% c('T+','T0','T-')
      if (any(bad)) stop('Invalid trend-state values in ', nm, ' of ', path, call. = FALSE)
    }
  }
  if (!'roman' %in% names(d)) d$roman <- NA_character_
  d$roman <- v12_zone_value(d$roman)
  d
}

v12_ensure_area <- function(d, cell_area_default = 25) {
  if (!'area_km2' %in% names(d)) d$area_km2 <- cell_area_default
  d$area_km2 <- suppressWarnings(as.numeric(d$area_km2))
  if (any(!is.finite(d$area_km2) | d$area_km2 <= 0)) {
    stop('area_km2 must be finite and positive for every cell.', call. = FALSE)
  }
  d
}

v12_classify6 <- function(rf, trend_state, rf_threshold) {
  rf <- suppressWarnings(as.numeric(rf))
  trend_state <- as.character(trend_state)
  rf_state <- ifelse(is.finite(rf) & rf >= rf_threshold, 'H',
                     ifelse(is.finite(rf) & rf < rf_threshold, 'L', NA_character_))
  out <- rep(NA_character_, length(rf))
  out[rf_state == 'H' & trend_state == 'T+'] <- 'H+'
  out[rf_state == 'H' & trend_state == 'T0'] <- 'H0'
  out[rf_state == 'H' & trend_state == 'T-'] <- 'H-'
  out[rf_state == 'L' & trend_state == 'T+'] <- 'L+'
  out[rf_state == 'L' & trend_state == 'T0'] <- 'L0'
  out[rf_state == 'L' & trend_state == 'T-'] <- 'L-'
  out
}

v12_normalise_support_mode <- function(x) {
  x <- tolower(trimws(as.character(x)[1]))
  aliases <- c(common='common', paired='common', intersection='common',
               target_specific='target_specific', targetspecific='target_specific', separate='target_specific')
  if (!x %in% names(aliases)) stop('TARGET_CONTRAST_SUPPORT must be common or target_specific; got ', x, call. = FALSE)
  unname(aliases[x])
}

v12_normalise_block_universe <- function(x) {
  x <- tolower(trimws(as.character(x)[1]))
  aliases <- c(analysis_support='analysis_support', support='analysis_support', active='analysis_support',
               all_rows='all_rows', legacy='all_rows', global='all_rows')
  if (!x %in% names(aliases)) stop('BOOT_BLOCK_UNIVERSE must be analysis_support or all_rows; got ', x, call. = FALSE)
  unname(aliases[x])
}

v12_block_id <- function(d, block_size_km, origin_x = 0, origin_y = 0) {
  if (!is.finite(block_size_km) || block_size_km <= 0) stop('block_size_km must be positive.', call. = FALSE)
  if (!is.finite(origin_x) || !is.finite(origin_y)) stop('Block origins must be finite.', call. = FALSE)
  bs <- as.numeric(block_size_km) * 1000
  paste0(floor((d$x - origin_x) / bs), '_', floor((d$y - origin_y) / bs))
}

v12_prepare_block_layout <- function(d, block_size_km, cell_area_default = 25,
                                     origin_x = 0, origin_y = 0) {
  d <- v12_ensure_area(d, cell_area_default)
  bid <- v12_block_id(d, block_size_km, origin_x, origin_y)
  blocks <- sort(unique(bid))
  list(
    n_rows = nrow(d), block_size_km = as.numeric(block_size_km),
    origin_x = as.numeric(origin_x), origin_y = as.numeric(origin_y),
    blocks = blocks, block_idx = match(bid, blocks),
    zone_idx = match(v12_zone_value(d$roman), V12_REG_ORDER),
    area = d$area_km2
  )
}

.v12_boot_cache <- new.env(parent = emptyenv())

v12_iter_seeds <- function(base_seed, n) {
  base_seed <- suppressWarnings(as.integer(base_seed)); n <- suppressWarnings(as.integer(n))
  if (is.na(base_seed)) stop('Invalid bootstrap seed.', call. = FALSE)
  if (is.na(n) || n <= 0L) stop('N_BOOT must be positive.', call. = FALSE)
  old_seed <- if (exists('.Random.seed', envir=.GlobalEnv, inherits=FALSE)) get('.Random.seed', envir=.GlobalEnv) else NULL
  on.exit({
    if (!is.null(old_seed)) assign('.Random.seed', old_seed, envir=.GlobalEnv)
    else if (exists('.Random.seed', envir=.GlobalEnv, inherits=FALSE)) rm('.Random.seed', envir=.GlobalEnv)
  }, add=TRUE)
  set.seed(base_seed)
  sample.int(.Machine$integer.max, n, replace = (n > .Machine$integer.max))
}

v12_seed_offset <- function(base_seed, ...) {
  vals <- suppressWarnings(as.numeric(unlist(list(...))))
  vals[!is.finite(vals)] <- 0
  # Keep within R's positive integer seed range while remaining deterministic.
  z <- (as.numeric(base_seed) + sum((seq_along(vals) * 1000003) * abs(vals))) %% (.Machine$integer.max - 1)
  as.integer(max(1, floor(z)))
}

v12_bootstrap_weights <- function(n_blocks, n_boot, seed) {
  n_blocks <- as.integer(n_blocks); n_boot <- as.integer(n_boot); seed <- as.integer(seed)
  if (is.na(n_blocks) || is.na(n_boot) || n_blocks <= 0L || n_boot <= 0L) stop('Bootstrap needs positive numbers of blocks and replicates.', call. = FALSE)
  key <- paste(n_blocks, n_boot, seed, sep='|')
  if (exists(key, envir=.v12_boot_cache, inherits=FALSE)) return(get(key, envir=.v12_boot_cache, inherits=FALSE))
  seeds <- v12_iter_seeds(seed, n_boot)
  W <- matrix(0L, nrow=n_boot, ncol=n_blocks)
  for (i in seq_len(n_boot)) {
    set.seed(seeds[i])
    W[i, ] <- tabulate(sample.int(n_blocks, n_blocks, replace=TRUE), nbins=n_blocks)
  }
  assign(key, W, envir=.v12_boot_cache)
  W
}

# Build a block x zone x class area cube using a precomputed spatial layout.
# rowsum on a linear array index avoids repeated data.frame + aggregate overhead.
v12_area_cube <- function(d, class_col, block_size_km, support_mask = NULL,
                          cell_area_default = 25, layout = NULL,
                          origin_x = 0, origin_y = 0) {
  d <- v12_ensure_area(d, cell_area_default)
  if (!class_col %in% names(d)) stop('Missing class column: ', class_col, call. = FALSE)
  if (is.null(layout)) layout <- v12_prepare_block_layout(d, block_size_km, cell_area_default, origin_x, origin_y)
  if (!is.list(layout) || layout$n_rows != nrow(d) || length(layout$block_idx) != nrow(d)) stop('Invalid/reused block layout for current cell table.', call. = FALSE)
  if (abs(layout$block_size_km - block_size_km) > 1e-12) stop('Block layout size mismatch.', call. = FALSE)
  if (is.null(support_mask)) support_mask <- rep(TRUE, nrow(d))
  support_mask <- as.logical(support_mask); support_mask[is.na(support_mask)] <- FALSE
  cls <- as.character(d[[class_col]])
  ci <- match(cls, V12_CLASS6_LEVELS)
  keep <- support_mask & !is.na(ci)
  nb <- length(layout$blocks); nz <- length(V12_REG_ORDER); nc <- length(V12_CLASS6_LEVELS)
  flat <- numeric(nb * nz * nc)
  if (any(keep)) {
    b <- layout$block_idx[keep]; z <- layout$zone_idx[keep]; cidx <- ci[keep]; a <- layout$area[keep]
    overall <- match('Overall', V12_REG_ORDER)
    lin_zone <- b + (z - 1L) * nb + (cidx - 1L) * nb * nz
    lin_all <- b + (overall - 1L) * nb + (cidx - 1L) * nb * nz
    lin <- c(lin_zone, lin_all); val <- c(a, a)
    sm <- rowsum(val, group = lin, reorder = FALSE)
    flat[as.integer(rownames(sm))] <- sm[,1]
  }
  cube <- array(flat, dim=c(nb,nz,nc),
                dimnames=list(block=layout$blocks, zone=V12_REG_ORDER, class6=V12_CLASS6_LEVELS))
  list(cube=cube, blocks=layout$blocks, layout=layout)
}

v12_cube_block_area <- function(cube, zone = 'Overall') {
  zi <- match(zone, dimnames(cube)$zone)
  if (is.na(zi)) stop('Unknown cube zone: ', zone, call. = FALSE)
  mat <- matrix(cube[,zi,], nrow=dim(cube)[1], ncol=dim(cube)[3])
  rowSums(mat)
}

v12_select_cube_blocks <- function(cube, active) {
  active <- as.logical(active); active[is.na(active)] <- FALSE
  if (length(active) != dim(cube)[1]) stop('Active-block mask length mismatch.', call. = FALSE)
  cube[active,,,drop=FALSE]
}

v12_cube_zone_block_counts <- function(cube) {
  out <- setNames(integer(dim(cube)[2]), dimnames(cube)$zone)
  for (zi in seq_len(dim(cube)[2])) {
    mat <- matrix(cube[,zi,], nrow=dim(cube)[1], ncol=dim(cube)[3])
    out[zi] <- sum(rowSums(mat) > 0)
  }
  out
}

# Block bootstrap optimization: batch the per-zone GEMMs into a single large GEMM.
# The zone x class slices are flattened to one [n_blocks x (nz*nc)] matrix so the
# bootstrap propagation is one W %*% BIG instead of nz separate small products, and
# the point props use one colSums. Bit-identical to the per-zone loop (verified: max
# abs error 0, NaN pattern identical) and still passes the 1e-12 self-test tolerance.
v12_cube_point_props <- function(cube) {
  nblk <- dim(cube)[1]; nz <- dim(cube)[2]; nc <- dim(cube)[3]
  out <- matrix(NA_real_, nrow=nz, ncol=nc,
                dimnames=list(zone=dimnames(cube)$zone, class6=dimnames(cube)$class6))
  if (nblk == 0L) return(out)
  mat  <- matrix(cube, nrow=nblk, ncol=nz*nc)     # column order: zone fastest, class next
  area <- matrix(colSums(mat), nrow=nz, ncol=nc)  # [zone, class]
  totals <- rowSums(area)
  ok <- is.finite(totals) & totals > 0
  if (any(ok)) out[ok,] <- area[ok,,drop=FALSE] / totals[ok]
  out
}

v12_cube_boot_props <- function(cube, W) {
  if (ncol(W) != dim(cube)[1]) stop('Bootstrap weight/cube block mismatch.', call. = FALSE)
  nb <- nrow(W); nblk <- dim(cube)[1]; nz <- dim(cube)[2]; nc <- dim(cube)[3]
  out <- array(NA_real_, dim=c(nb,nz,nc),
               dimnames=list(iter=seq_len(nb), zone=dimnames(cube)$zone, class6=dimnames(cube)$class6))
  if (nblk == 0L) return(out)
  big   <- matrix(cube, nrow=nblk, ncol=nz*nc)    # [nblk, nz*nc]
  areas <- W %*% big                              # single large GEMM -> [nb, nz*nc]
  arr   <- array(areas, dim=c(nb,nz,nc))          # back to [nb, nz, nc]
  den <- matrix(0, nb, nz)
  for (ci in seq_len(nc)) den <- den + arr[,,ci]  # sum over class (nc=6, cheap)
  ok <- is.finite(den) & den > 0
  for (ci in seq_len(nc)) {
    slc <- arr[,,ci]
    slc[!ok] <- NA_real_
    slc[ok]  <- slc[ok] / den[ok]
    out[,,ci] <- slc
  }
  out
}

v12_safe_quantile <- function(x, p) {
  x <- as.numeric(x); x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  as.numeric(stats::quantile(x, probs=p, names=FALSE, type=7, na.rm=TRUE))
}

v12_props_ci_long <- function(point, boot, value_name='prop', ci_level=0.95,
                              zone_block_counts = NULL, n_blocks_universe = NA_integer_) {
  if (!is.finite(ci_level) || ci_level <= 0 || ci_level >= 1) stop('CI_LEVEL must be in (0,1).', call. = FALSE)
  alpha <- (1-ci_level)/2
  rows <- vector('list', length(V12_REG_ORDER) * length(V12_CLASS6_LEVELS)); k <- 0L
  for (z in V12_REG_ORDER) for (cl in V12_CLASS6_LEVELS) {
    k <- k + 1L; x <- boot[,z,cl]
    rows[[k]] <- data.frame(zone=z, class6=cl,
      point=point[z,cl], median=v12_safe_quantile(x,0.5),
      low=v12_safe_quantile(x,alpha), high=v12_safe_quantile(x,1-alpha),
      n_boot=sum(is.finite(x)), bootstrap_valid_fraction=mean(is.finite(x)),
      n_blocks_universe=as.integer(n_blocks_universe),
      n_blocks_zone=if (is.null(zone_block_counts)) NA_integer_ else as.integer(zone_block_counts[z]),
      ci_level=ci_level, stringsAsFactors=FALSE)
  }
  out <- do.call(rbind, rows)
  names(out)[names(out)=='point'] <- value_name
  names(out)[names(out)=='median'] <- paste0(value_name,'_median')
  names(out)[names(out)=='low'] <- paste0(value_name,'_low')
  names(out)[names(out)=='high'] <- paste0(value_name,'_high')
  out
}

# SI S6.3: resample blocks WITHIN each reporting unit. Overall is its own
# reporting unit. Shared support and seed preserve metric pairing in each unit.
v161_boot_props_by_zone <- function(cube,n_boot,seed,support=cube) {
  if(!identical(dim(cube),dim(support)))stop('Paired cube dimensions differ')
  ans<-array(NA_real_,dim=c(n_boot,dim(cube)[2],dim(cube)[3]),
             dimnames=list(iter=seq_len(n_boot),zone=dimnames(cube)[[2]],class6=dimnames(cube)[[3]]))
  for(g in seq_len(dim(cube)[2])) {
    support_mat<-matrix(support[,g,],nrow=dim(cube)[1],ncol=dim(cube)[3])
    active<-rowSums(support_mat)>0
    if(!any(active))next
    values<-matrix(cube[,g,],nrow=dim(cube)[1],ncol=dim(cube)[3])[active,,drop=FALSE]
    W<-v12_bootstrap_weights(nrow(values),n_boot,v12_seed_offset(seed,g))
    area<-W %*% values;den<-rowSums(area);ok<-den>0
    ans[ok,g,]<-area[ok,,drop=FALSE]/den[ok]
  }
  ans
}

v12_area_ci_metric <- function(d, class_col, metric_label, block_size_km, n_boot,
                               seed, cell_area_default=25, layout=NULL,
                               block_universe='analysis_support', ci_level=0.95,
                               origin_x=0, origin_y=0) {
  d <- v12_ensure_area(d, cell_area_default)
  block_universe <- v12_normalise_block_universe(block_universe)
  if (is.null(layout)) layout <- v12_prepare_block_layout(d, block_size_km, cell_area_default, origin_x, origin_y)
  valid <- as.character(d[[class_col]]) %in% V12_CLASS6_LEVELS
  obj <- v12_area_cube(d, class_col, block_size_km, support_mask=valid,
                       cell_area_default=cell_area_default, layout=layout)
  active <- if (block_universe=='analysis_support') v12_cube_block_area(obj$cube)>0 else rep(TRUE,length(obj$blocks))
  cube <- v12_select_cube_blocks(obj$cube, active)
  point <- v12_cube_point_props(cube)
  if (dim(cube)[1] > 0L) {
    if(block_universe=='all_rows') {
      W <- v12_bootstrap_weights(dim(cube)[1], n_boot, seed)
      boot <- v12_cube_boot_props(cube,W) # explicit historical profile only
    } else boot <- v161_boot_props_by_zone(cube,n_boot,seed)
  } else {
    boot <- array(NA_real_, dim=c(n_boot,length(V12_REG_ORDER),length(V12_CLASS6_LEVELS)),
                  dimnames=list(iter=seq_len(n_boot),zone=V12_REG_ORDER,class6=V12_CLASS6_LEVELS))
  }
  zbc <- v12_cube_zone_block_counts(cube)
  out <- v12_props_ci_long(point, boot, 'prop', ci_level, zbc, dim(cube)[1])
  dd <- d[valid,,drop=FALSE]
  rows <- list(); k <- 0L
  for (z in V12_REG_ORDER) {
    dz <- if (z=='Overall') dd else dd[v12_zone_value(dd$roman)==z,,drop=FALSE]
    for (cl in V12_CLASS6_LEVELS) {
      k <- k+1L; hit <- as.character(dz[[class_col]])==cl
      rows[[k]] <- data.frame(zone=z,class6=cl,n_cells=sum(hit,na.rm=TRUE),
                              area_km2=sum(dz$area_km2[hit],na.rm=TRUE), stringsAsFactors=FALSE)
    }
  }
  counts <- do.call(rbind, rows)
  merge(out, counts, by=c('zone','class6'), all.x=TRUE, sort=FALSE) |>
    transform(metric=metric_label, block_size_km=block_size_km,
              block_universe_mode=block_universe,
              block_origin_x=origin_x, block_origin_y=origin_y)
}

v12_metric_contrast <- function(d, block_size_km, n_boot, seed,
                                support_mode='common', cell_area_default=25,
                                layout=NULL, block_universe='analysis_support',
                                ci_level=0.95, origin_x=0, origin_y=0) {
  d <- v12_ensure_area(d, cell_area_default)
  support_mode <- v12_normalise_support_mode(support_mode)
  block_universe <- v12_normalise_block_universe(block_universe)
  if (is.null(layout)) layout <- v12_prepare_block_layout(d, block_size_km, cell_area_default, origin_x, origin_y)
  rich_valid <- as.character(d$rich_class6) %in% V12_CLASS6_LEVELS
  shan_valid <- as.character(d$shannon_class6) %in% V12_CLASS6_LEVELS
  common <- rich_valid & shan_valid
  rich_support <- if (support_mode=='common') common else rich_valid
  shan_support <- if (support_mode=='common') common else shan_valid
  r <- v12_area_cube(d, 'rich_class6', block_size_km, rich_support, cell_area_default, layout=layout)
  s <- v12_area_cube(d, 'shannon_class6', block_size_km, shan_support, cell_area_default, layout=layout)
  if (!identical(r$blocks,s$blocks)) stop('Internal block-universe mismatch.', call. = FALSE)
  rba <- v12_cube_block_area(r$cube); sba <- v12_cube_block_area(s$cube)
  active <- if (block_universe=='analysis_support') (rba>0 | sba>0) else rep(TRUE,length(r$blocks))
  rc <- v12_select_cube_blocks(r$cube,active); sc <- v12_select_cube_blocks(s$cube,active)
  pr <- v12_cube_point_props(rc); ps <- v12_cube_point_props(sc)
  if (dim(rc)[1] > 0L) {
    if(block_universe=='all_rows') {
      W <- v12_bootstrap_weights(dim(rc)[1], n_boot, seed)
      br <- v12_cube_boot_props(rc,W); bs <- v12_cube_boot_props(sc,W)
    } else {
      support<-rc+sc
      br<-v161_boot_props_by_zone(rc,n_boot,seed,support)
      bs<-v161_boot_props_by_zone(sc,n_boot,seed,support)
    }
    delta <- bs - br
  } else {
    delta <- array(NA_real_,dim=c(n_boot,length(V12_REG_ORDER),length(V12_CLASS6_LEVELS)),
                   dimnames=list(iter=seq_len(n_boot),zone=V12_REG_ORDER,class6=V12_CLASS6_LEVELS))
  }
  alpha <- (1-ci_level)/2
  rows <- vector('list', length(V12_REG_ORDER)*length(V12_CLASS6_LEVELS)); k <- 0L
  zones <- v12_zone_value(d$roman)
  rzc <- v12_cube_zone_block_counts(rc); szc <- v12_cube_zone_block_counts(sc)
  for (z in V12_REG_ORDER) for (cl in V12_CLASS6_LEVELS) {
    k <- k+1L; x <- delta[,z,cl]
    zm <- if (z=='Overall') rep(TRUE,nrow(d)) else zones==z
    q50 <- v12_safe_quantile(x,0.5); qlo <- v12_safe_quantile(x,alpha); qhi <- v12_safe_quantile(x,1-alpha)
    rows[[k]] <- data.frame(
      zone=z,class6=cl,support_mode=support_mode,
      prop_shannon=ps[z,cl],prop_richness=pr[z,cl],
      delta_metric=ps[z,cl]-pr[z,cl],
      delta_median=q50,delta_low=qlo,delta_high=qhi,
      ci_excludes_zero=is.finite(qlo) && is.finite(qhi) && (qlo>0 || qhi<0),
      n_boot=sum(is.finite(x)),bootstrap_valid_fraction=mean(is.finite(x)),ci_level=ci_level,
      n_blocks_universe=dim(rc)[1],n_blocks_rich_zone=rzc[z],n_blocks_shannon_zone=szc[z],
      block_universe_mode=block_universe,block_origin_x=origin_x,block_origin_y=origin_y,
      n_union_cells=sum(zm & (rich_valid | shan_valid)),
      n_common_cells=sum(zm & common),
      common_area_km2=sum(d$area_km2[zm & common],na.rm=TRUE),
      rich_classified_area_km2=sum(d$area_km2[zm & rich_valid],na.rm=TRUE),
      shannon_classified_area_km2=sum(d$area_km2[zm & shan_valid],na.rm=TRUE),
      stringsAsFactors=FALSE)
  }
  do.call(rbind,rows)
}

# Area-weighted (not ordinal-category-weighted) Cohen's kappa.
v12_weighted_kappa <- function(tab) {
  tab <- as.matrix(tab); total <- sum(tab)
  if (!is.finite(total) || total <= 0) return(NA_real_)
  p <- tab/total; po <- sum(diag(p)); pe <- sum(rowSums(p)*colSums(p))
  if (!is.finite(pe) || abs(1-pe)<1e-15) return(NA_real_)
  (po-pe)/(1-pe)
}

v12_target_transition <- function(d, cell_area_default=25) {
  d <- v12_ensure_area(d,cell_area_default)
  rv <- as.character(d$rich_class6) %in% V12_CLASS6_LEVELS
  sv <- as.character(d$shannon_class6) %in% V12_CLASS6_LEVELS
  common <- rv & sv; zones <- v12_zone_value(d$roman)
  trans <- list(); agree <- list(); coverage <- list(); kt <- 0L; ka <- 0L; kc <- 0L
  for (z in V12_REG_ORDER) {
    zm <- if (z=='Overall') zones %in% c('I','II','III','IV','V') else zones==z
    ok <- zm & common
    tab <- matrix(0,nrow=6,ncol=6,dimnames=list(rich_class6=V12_CLASS6_LEVELS,shannon_class6=V12_CLASS6_LEVELS))
    if (any(ok)) {
      ag <- stats::aggregate(d$area_km2[ok], by=list(rich=as.character(d$rich_class6[ok]),shannon=as.character(d$shannon_class6[ok])), FUN=sum)
      tab[cbind(match(ag$rich,V12_CLASS6_LEVELS),match(ag$shannon,V12_CLASS6_LEVELS))] <- ag$x
    }
    total <- sum(tab)
    for (i in seq_along(V12_CLASS6_LEVELS)) for (j in seq_along(V12_CLASS6_LEVELS)) {
      kt<-kt+1L; trans[[kt]]<-data.frame(zone=z,rich_class6=V12_CLASS6_LEVELS[i],shannon_class6=V12_CLASS6_LEVELS[j],
        area_km2=tab[i,j],prop_common_area=if(total>0) tab[i,j]/total else NA_real_,stringsAsFactors=FALSE)
    }
    if (any(ok)) {
      rr <- as.character(d$rich_class6[ok]); ss <- as.character(d$shannon_class6[ok]); w <- d$area_km2[ok]
      wmean <- function(x) sum(w*as.numeric(x),na.rm=TRUE)/sum(w)
      exact <- wmean(rr==ss); hl <- wmean(substr(rr,1,1)==substr(ss,1,1)); tr <- wmean(substr(rr,2,2)==substr(ss,2,2))
    } else exact<-hl<-tr<-NA_real_
    kap <- v12_weighted_kappa(tab)
    ka<-ka+1L; agree[[ka]]<-data.frame(zone=z,n_common_cells=sum(ok),common_area_km2=total,
      exact_6class_agreement=exact,hl_state_agreement=hl,trend_state_agreement=tr,
      area_weighted_cohen_kappa=kap,weighted_cohen_kappa=kap,stringsAsFactors=FALSE)
    for (metric in c('Richness','Shannon')) {
      vv <- if(metric=='Richness') rv else sv; cc <- if(metric=='Richness') as.character(d$rich_class6) else as.character(d$shannon_class6)
      den <- sum(d$area_km2[zm & vv],na.rm=TRUE); foc <- sum(d$area_km2[zm & vv & cc %in% V12_FOCAL4_LEVELS],na.rm=TRUE)
      kc<-kc+1L; coverage[[kc]]<-data.frame(zone=z,metric=metric,n_classified_cells=sum(zm&vv),classified_area_km2=den,
        focal4_area_km2=foc,focal4_share=if(den>0) foc/den else NA_real_,stringsAsFactors=FALSE)
    }
  }
  list(transition=do.call(rbind,trans),agreement=do.call(rbind,agree),coverage=do.call(rbind,coverage))
}
