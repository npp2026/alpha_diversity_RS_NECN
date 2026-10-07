#!/usr/bin/env Rscript
# ==============================================================================
# Publication-facing 1-km state-trajectory workflow (curated from compute_RF_6class_5km_trend_v9.R)
# ------------------------------------------------------------------------------
# Target-grid (default 1 km) annual time series -> trend-period-specific
# Sen/Hamed-Rao MK -> BH-FDR -> fixed RF period x trend-state six-class classification.
# v1.6 MS: RF = mean observed diversity in 2016-2020 / one static q95 surface at age 100.
#
# This script intentionally does NOT use 30 m significant-pixel fractions, SFM,
# or within-grid field significance for main classification.
# ==============================================================================
options(stringsAsFactors = FALSE)
if (getRversion() < "4.1.0") stop("State-trajectory analysis requires R >= 4.1.0 because scripts use the base pipe operator |>. Please upgrade R.", call. = FALSE)

suppressPackageStartupMessages({
  library(terra)
  library(sf)
  library(dplyr)
  library(tidyr)
})

# 1-km runs create many more cells than the original 5-km workflow.  Keep terra's
# in-memory cache conservative by default and allow users to redirect temporary
# raster blocks to a fast local disk via TERRA_TMPDIR.
terra_tmpdir <- Sys.getenv('TERRA_TMPDIR', unset = '')
if (nzchar(terra_tmpdir)) dir.create(terra_tmpdir, recursive = TRUE, showWarnings = FALSE)
terra_memfrac <- suppressWarnings(as.numeric(Sys.getenv('TERRA_MEMFRAC', unset = '0.55')))
if (!is.finite(terra_memfrac) || terra_memfrac <= 0 || terra_memfrac > 0.9) terra_memfrac <- 0.55
try({
  if (nzchar(terra_tmpdir)) terra::terraOptions(tempdir = terra_tmpdir)
  terra::terraOptions(memfrac = terra_memfrac)
}, silent = TRUE)

# ---- environment helpers -----------------------------------------------------
# Configuration parsers are bound to the shared v12 helpers below.
msg <- function(...) cat(sprintf(...), "\n")

get_script_dir <- function() {
  args <- commandArgs(FALSE)
  file_arg <- args[grep('^--file=', args)]
  if (length(file_arg)) return(dirname(normalizePath(sub('^--file=', '', file_arg[1]), winslash='/', mustWork=FALSE)))
  getwd()
}
SCRIPT_DIR <- get_script_dir()
FEM_SOFTWARE_VERSION <- trimws(readLines(file.path(SCRIPT_DIR, "..", "VERSION"), n = 1L, warn = FALSE))
UTILS_FILE <- file.path(SCRIPT_DIR, 'v9_parallel_utils.R')
if (file.exists(UTILS_FILE)) source(UTILS_FILE)
V12_UTILS_FILE <- file.path(SCRIPT_DIR, 'v12_posthoc_utils.R')
if (!file.exists(V12_UTILS_FILE)) stop('Missing shared utility file: ', V12_UTILS_FILE, call.=FALSE)
source(V12_UTILS_FILE)
# Use strict parsers for production configuration: malformed values now stop rather
# than silently reverting to defaults.
env_str <- v12_env_str
env_logical <- v12_env_bool
env_num <- function(name, default) v12_env_num(name, default)
env_int <- function(name, default) v12_env_int(name, default)
env_expr <- v12_env_expr
make_abs_path <- v12_make_abs_path
TREND_UTILS_FILE <- file.path(SCRIPT_DIR, 'v9_trend_package_utils.R')
if (file.exists(TREND_UTILS_FILE)) {
  source(TREND_UTILS_FILE)
  v9_check_professional_packages(strict = V9_STRICT_PROFESSIONAL_PACKAGES)
}


# ---- parameters --------------------------------------------------------------
DATA_DIR <- env_str("DATA_DIR", ".")
DATA_DIR_ABS <- normalizePath(DATA_DIR, winslash = "/", mustWork = TRUE)
OUT_ROOT_NAME <- env_str("OUT_ROOT_NAME", "outputs_RF_6class_fig6_minimal")
OUT_ROOT <- make_abs_path(OUT_ROOT_NAME, DATA_DIR_ABS)
CLEAN_EXISTING <- env_logical("CLEAN_EXISTING", FALSE)
if (CLEAN_EXISTING && dir.exists(OUT_ROOT)) {
  out_norm <- normalizePath(OUT_ROOT, winslash="/", mustWork=TRUE)
  data_norm <- normalizePath(DATA_DIR_ABS, winslash="/", mustWork=TRUE)
  if (identical(out_norm, data_norm) || identical(out_norm,"/") || grepl("^[A-Za-z]:/$",out_norm)) {
    stop("Refusing to clean unsafe OUT_ROOT: ", out_norm, call.=FALSE)
  }
  unlink(OUT_ROOT, recursive=TRUE, force=TRUE)
  if (dir.exists(OUT_ROOT)) stop("Failed to clean existing OUT_ROOT: ", OUT_ROOT, call.=FALSE)
}
dir.create(OUT_ROOT, recursive = TRUE, showWarnings = FALSE)

CELL_SIZE_M <- env_int("CELL_SIZE_M", 1000)
YEARS <- v12_validate_year_vector(env_expr("YEARS_R", 2001:2020), "YEARS_R", min_n=3L)
MIN_VALID_FRAC <- env_num("MIN_VALID_FRAC", 0.80)
# Trend validity is now proportional to each trend-period length so 2005-2020
# and 2001-2020 do not use asymmetric 100% vs 80% completeness rules.
TREND_MIN_VALID_YEAR_FRAC <- env_num("TREND_MIN_VALID_YEAR_FRAC", 0.80)
RF_THRESH <- env_num("RF_THRESH", 0.80)
Q_MAIN <- env_num("Q_MAIN", 0.05)
Q_WILKS_SENS <- env_num("Q_WILKS_SENS", 0.10)
MIN_DENOM <- env_num("MIN_DENOM", 1e-6)
# Annual-potential overrides are incompatible with the static age-100 contract.
POT_REF_YEARS_OVERRIDE <- env_str("POT_REF_YEARS_R", "")
if (nzchar(POT_REF_YEARS_OVERRIDE)) stop("POT_REF_YEARS_R is incompatible with the static age-100 profile.")
ALLOW_PARTIAL_OBS_YEARS <- env_logical("ALLOW_PARTIAL_OBS_YEARS", TRUE)
OBS_YEARS_MIN_FRAC <- env_num("OBS_YEARS_MIN_FRAC", 0.80)
HR_LAG_MAX <- env_int("HR_LAG_MAX", 0)  # 0 means use all available lags (n-1) in fallback, closer to modifiedmk::mmkh
HR_USE_SIGNIFICANT_ACF_ONLY <- env_logical("HR_USE_SIGNIFICANT_ACF_ONLY", TRUE)
WRITE_RASTERS <- env_logical("WRITE_RASTERS", TRUE)
WRITE_ANNUAL_TABLES <- env_logical("WRITE_ANNUAL_TABLES", FALSE)
WRITE_COMBINED_CELL_TABLE <- env_logical("WRITE_COMBINED_CELL_TABLE", FALSE)
ALLOW_AMBIGUOUS_INPUT_FILES <- env_logical("ALLOW_AMBIGUOUS_INPUT_FILES", FALSE)
ALLOW_REGION_ORDER_FALLBACK <- env_logical("ALLOW_REGION_ORDER_FALLBACK", FALSE)
INPUT_NEGATIVE_TOL <- env_num("INPUT_NEGATIVE_TOL", -1e-10)
N_WORKERS <- env_int("N_WORKERS", 1)
# Step-02 parallel mode (no effect on results; lets you profile both on your CPU/disk):
#   TRUE  (default): metric-outer parallel -- the 2 metrics run concurrently on 2 PSOCK
#                    workers, each fitting trends serially. Caps at ~2 cores but overlaps
#                    the (often I/O-bound) terra raster reads of the two metrics.
#   FALSE          : serial-per-metric -- metrics run one at a time, but each metric's
#                    trend fit parallelises across N_WORKERS cores. Reads I/O sequentially
#                    but uses all cores for the CPU-bound Sen/MK fitting (the original v9
#                    design). Faster when trend fitting dominates rather than I/O.
# Results are bit-identical in both modes (trend fits are deterministic; row order is
# preserved). Try FALSE with a high N_WORKERS if step 02 is CPU-bound on your machine.
MAIN_METRIC_PARALLEL <- env_logical("MAIN_METRIC_PARALLEL", TRUE)
USE_PROFESSIONAL_TREND_PACKAGES <- env_logical("USE_PROFESSIONAL_TREND_PACKAGES", TRUE)
# trend::sens.slope() can emit many sqrt(VS) warnings for tied/near-constant series.
# The project fallback implements the same Sen median slope and uses true year spacing.
# Default FALSE keeps main runs quiet and deterministic; set TRUE for package-parity checks.
USE_TREND_PACKAGE_SEN <- env_logical("USE_TREND_PACKAGE_SEN", FALSE)
USE_MODIFIEDMK_PACKAGE <- env_logical("USE_MODIFIEDMK_PACKAGE", FALSE)
STRICT_PROFESSIONAL_PACKAGES <- env_logical("STRICT_PROFESSIONAL_PACKAGES", FALSE)

if (!is.finite(CELL_SIZE_M) || CELL_SIZE_M < 1L) stop("CELL_SIZE_M must be a positive integer number of meters; got ", CELL_SIZE_M, call.=FALSE)
GRID_TAG <- sprintf("%dm", CELL_SIZE_M)
GRID_LABEL <- if (CELL_SIZE_M %% 1000L == 0L) sprintf("%g km", CELL_SIZE_M / 1000) else sprintf("%d m", CELL_SIZE_M)
TEMPLATE_FILE <- sprintf("template_%s.tif", GRID_TAG)
REGION_FILE <- sprintf("region_id_%s.tif", GRID_TAG)
MAP_FILE <- sprintf("RF_6class_%s.tif", GRID_TAG)
if (!is.finite(MIN_VALID_FRAC) || MIN_VALID_FRAC <= 0 || MIN_VALID_FRAC > 1) stop("MIN_VALID_FRAC must be in (0,1].", call.=FALSE)
if (!is.finite(OBS_YEARS_MIN_FRAC) || OBS_YEARS_MIN_FRAC <= 0 || OBS_YEARS_MIN_FRAC > 1) stop("OBS_YEARS_MIN_FRAC must be in (0,1].", call.=FALSE)
if (!is.finite(RF_THRESH) || RF_THRESH <= 0) stop("RF_THRESH must be positive.", call.=FALSE)
if (!is.finite(Q_MAIN) || Q_MAIN <= 0 || Q_MAIN >= 1) stop("Q_MAIN must be in (0,1).", call.=FALSE)
if (!is.finite(MIN_DENOM) || MIN_DENOM < 0) stop("MIN_DENOM must be non-negative.", call.=FALSE)
if (!is.finite(INPUT_NEGATIVE_TOL) || INPUT_NEGATIVE_TOL > 0) stop("INPUT_NEGATIVE_TOL must be finite and <= 0.", call.=FALSE)
if (N_WORKERS < 1L) stop("N_WORKERS must be >=1.", call.=FALSE)

POTENTIAL_DIR <- v12_make_abs_path(env_str("POTENTIAL_DIR", "potential_q95"), DATA_DIR_ABS)
MOUNTAIN_SHP <- file.path(DATA_DIR_ABS, "NE_Mountain_Output", "NE_Mountain_Regions_All.shp")

METRICS <- list(
  list(pref = "rich", stem = "Rich_tree", label = "Richness"),
  list(pref = "shannon", stem = "Shannon_wiener", label = "Shannon")
)

# RF period is now fixed to the late-stage observed mean by default.
# Trend periods are separate analysis axes and are NOT RF periods.
PERIODS <- list(
  `2016_2020_mean` = list(stage_label = "2016-2020 mean", obs_years = 2016:2020)
)
RUN_PERIODS_RAW <- unique(as.character(env_expr("RUN_PERIODS_R", "2016_2020_mean")))
BAD_RUN_PERIODS <- setdiff(RUN_PERIODS_RAW,names(PERIODS))
if(length(BAD_RUN_PERIODS))stop("Unknown RUN_PERIODS_R value(s): ",paste(BAD_RUN_PERIODS,collapse=","),". Allowed: ",paste(names(PERIODS),collapse=","),call.=FALSE)
RUN_PERIODS <- RUN_PERIODS_RAW
if (!length(RUN_PERIODS)) stop("RUN_PERIODS_R resolved to an empty vector.", call. = FALSE)

TREND_PERIODS <- list(
  trend_2005_2020 = list(label = "2005-2020 trend", years = 2005:2020),
  trend_2001_2020 = list(label = "2001-2020 trend", years = 2001:2020)
)
RUN_TREND_PERIODS_RAW <- unique(as.character(env_expr("RUN_TREND_PERIODS_R", "trend_2005_2020")))
BAD_TREND_PERIODS <- setdiff(RUN_TREND_PERIODS_RAW,names(TREND_PERIODS))
if(length(BAD_TREND_PERIODS))stop("Unknown RUN_TREND_PERIODS_R value(s): ",paste(BAD_TREND_PERIODS,collapse=","),". Allowed: ",paste(names(TREND_PERIODS),collapse=","),call.=FALSE)
RUN_TREND_PERIODS <- RUN_TREND_PERIODS_RAW
if (!length(RUN_TREND_PERIODS)) stop("RUN_TREND_PERIODS_R resolved to an empty vector.", call. = FALSE)
needed_trend_years <- sort(unique(unlist(lapply(TREND_PERIODS[RUN_TREND_PERIODS], `[[`, "years"))))
if (!all(needed_trend_years %in% YEARS)) {
  stop("YEARS_R must include all requested trend-period years. Missing: ", paste(setdiff(needed_trend_years, YEARS), collapse=","), call. = FALSE)
}
trend_min_valid_years <- function(n_total) {
  if (!is.finite(TREND_MIN_VALID_YEAR_FRAC) || TREND_MIN_VALID_YEAR_FRAC <= 0 || TREND_MIN_VALID_YEAR_FRAC > 1) {
    stop("TREND_MIN_VALID_YEAR_FRAC must be in (0, 1].", call. = FALSE)
  }
  # ceiling() returns a double and max(3L, <double>) promotes to double, so coerce
  # back to a true integer. This keeps the function's contract (a count of years)
  # honest and lets callers use vapply(..., integer(1)) safely.
  as.integer(max(3L, ceiling(TREND_MIN_VALID_YEAR_FRAC * n_total)))
}

CLASS6_LEVELS <- c("H+", "H0", "H-", "L+", "L0", "L-")
CLASS6_CODE <- c("H+" = 1L, "H0" = 2L, "H-" = 3L, "L+" = 4L, "L0" = 5L, "L-" = 6L)
CLASS6_LABEL <- c(
  "H+" = "High RF - significant increase 高实现-持续提升",
  "H0" = "High RF - not significant 高实现-无显著变化证据",
  "H-" = "High RF - significant decline 高实现-退化风险",
  "L+" = "Low RF - significant increase 低实现-恢复提升",
  "L0" = "Low RF - not significant 低实现-恢复不足",
  "L-" = "Low RF - significant decline 低实现-持续退化"
)
TREND_CODE <- c("T+" = 1L, "T0" = 0L, "T-" = -1L)
MOUNTAIN_LEVELS <- c("HLJDXAL", "HLJXXAL", "HLJCBS", "JLSCBS", "LNCBS")
REGION_CN <- c(HLJDXAL = "黑龙江大兴安岭", HLJXXAL = "黑龙江小兴安岭", HLJCBS  = "黑龙江长白山", JLSCBS  = "吉林省长白山", LNCBS   = "辽宁省长白山")
ROMAN_ID <- c(HLJDXAL = "I", HLJXXAL = "II", HLJCBS = "III", JLSCBS = "IV", LNCBS = "V")
REG_ORDER <- c("I", "II", "III", "IV", "V", "outside", "Overall")
zone_clean <- function(z) {
  z <- as.character(z)
  z[is.na(z) | !nzchar(z)] <- "outside"
  z
}

# ---- file discovery ----------------------------------------------------------
observed_annual_search_dirs <- function() {
  env_dirs <- env_str("OBS_ANNUAL_DIRS", "")
  if (nzchar(env_dirs)) {
    parts <- unlist(strsplit(env_dirs, ";|,", perl = TRUE))
    parts <- trimws(parts[nzchar(trimws(parts))])
  } else {
    parts <- c("mean", "1km", "obs", "obs1km", "annual", ".")
  }
  dirs <- vapply(parts, make_abs_path, character(1), base_dir = DATA_DIR_ABS)
  unique(dirs[file.exists(dirs) & dir.exists(dirs)])
}

is_safe_observed_annual_candidate <- function(path, stem, year) {
  b <- basename(path)
  if (!grepl(stem, b, fixed = TRUE)) return(FALSE)
  if (!grepl(as.character(year), b, fixed = TRUE)) return(FALSE)
  bad_pat <- "support_|quantiles|uncertainty|OOD_stack|valid_fraction|Pred_|q0\\.99|Trend|Mann|Slope|MK|2001_2005|2016_2020|2001_2020|2005_2020"
  if (grepl(bad_pat, b, ignore.case = TRUE)) return(FALSE)
  TRUE
}

find_observed_annual_file <- function(stem, year) {
  dirs <- observed_annual_search_dirs()
  exact_names <- c(sprintf("%s_%d_1km.tif",stem,year),sprintf("%s_%d.tif",stem,year),
    sprintf("%s_%d_mean_1km.tif",stem,year),sprintf("%s_%d_mean.tif",stem,year),
    sprintf("%s_mean_%d_1km.tif",stem,year),sprintf("%s_mean_%d.tif",stem,year))
  exact <- sort(unique(unlist(lapply(dirs,function(d)file.path(d,exact_names)[file.exists(file.path(d,exact_names))]),use.names=FALSE)))
  fuzzy <- sort(unique(unlist(lapply(dirs,function(d)list.files(d,pattern=paste0(".*",stem,".*",year,".*\\.tif$"),full.names=TRUE,recursive=FALSE)),use.names=FALSE)))
  fuzzy <- fuzzy[vapply(fuzzy,is_safe_observed_annual_candidate,logical(1),stem=stem,year=year)]
  # Consider every safe candidate, including non-standard names even when an
  # exact-name file exists. Otherwise a stale duplicate can be silently ignored.
  cand <- sort(unique(c(exact, fuzzy)))
  cand <- sort(unique(normalizePath(cand,winslash="/",mustWork=TRUE)))
  if(length(cand)>1L&&!ALLOW_AMBIGUOUS_INPUT_FILES)stop("Multiple annual candidates for ",stem," ",year,": ",paste(cand,collapse="; "),
    ". Remove duplicates or set ALLOW_AMBIGUOUS_INPUT_FILES=TRUE only after verifying equivalence.",call.=FALSE)
  if(length(cand)>1L)warning("Multiple annual candidates for ",stem," ",year,"; using deterministic first path: ",cand[1],call.=FALSE)
  if(length(cand))cand[1] else NA_character_
}

observed_annual_files <- function(stem, years) {
  files <- vapply(years, function(y) find_observed_annual_file(stem, y), character(1))
  names(files) <- as.character(years)
  files
}

# ---- spatial helpers ---------------------------------------------------------
make_region_table <- function(mountains_vect) {
  n <- nrow(mountains_vect); vals <- terra::values(mountains_vect); nm <- names(vals)
  pick <- function(target){i<-match(tolower(target),tolower(nm));if(is.na(i))NULL else vals[[i]]}
  code0<-pick('region_code');roman0<-pick('Roman_ID');cn0<-pick('Reg_CN')
  if(!is.null(code0)) code<-trimws(as.character(code0))
  else if(!is.null(roman0)) code<-names(ROMAN_ID)[match(trimws(as.character(roman0)),ROMAN_ID)]
  else if(!is.null(cn0)){rev_cn<-setNames(names(REGION_CN),REGION_CN);code<-unname(rev_cn[trimws(as.character(cn0))])}
  else if(n==5L&&ALLOW_REGION_ORDER_FALLBACK){code<-MOUNTAIN_LEVELS;warning('Using feature-order region fallback because ALLOW_REGION_ORDER_FALLBACK=TRUE.',call.=FALSE)}
  else stop('Mountain vector must contain region_code, Roman_ID, or Reg_CN with canonical I-V mapping. Feature-order fallback is disabled.',call.=FALSE)
  if(length(code)!=n||any(is.na(code))||any(!code%in%MOUNTAIN_LEVELS)||anyDuplicated(code))
    stop('Mountain-region attributes do not form a one-to-one canonical set: ',paste(MOUNTAIN_LEVELS,collapse=','),call.=FALSE)
  roman<-unname(ROMAN_ID[code])
  if(!setequal(roman,c('I','II','III','IV','V')))stop('Mountain-region Roman IDs must be exactly I-V.',call.=FALSE)
  data.frame(region_code=code,region_cn=unname(REGION_CN[code]),roman_id=roman,stringsAsFactors=FALSE)
}

vector_crs_is_empty <- function(v) {
  cr <- tryCatch(terra::crs(v), error = function(e) "")
  is.null(cr) || !nzchar(cr)
}

extent_overlaps <- function(a, b) {
  ea <- as.vector(terra::ext(a)); eb <- as.vector(terra::ext(b))
  isTRUE(ea[1] <= eb[2] && ea[2] >= eb[1] && ea[3] <= eb[4] && ea[4] >= eb[3])
}

project_vector_to_template <- function(v, template, label = "vector") {
  if (is.null(v) || nrow(v) == 0) return(NULL)
  if (vector_crs_is_empty(v)) stop(label, " has no CRS; cannot align safely.", call. = FALSE)
  if (!identical(terra::crs(v), terra::crs(template))) v <- terra::project(v, terra::crs(template))
  if (!extent_overlaps(v, template)) warning(label, " extent does not overlap template after CRS alignment.", call. = FALSE)
  v
}

prepare_mountains <- function(path, template) {
  if (!file.exists(path)) stop("Missing mountain-region vector: ", path, call. = FALSE)
  mtn <- project_vector_to_template(terra::vect(path), template, "Mountain-region vector")
  tab <- make_region_table(mtn)
  mtn$region_code <- tab$region_code
  mtn$Reg_CN <- tab$region_cn
  mtn$Roman_ID <- tab$roman_id
  mtn$rid <- match(tab$roman_id, c("I", "II", "III", "IV", "V"))
  mtn
}

stack_spatraster_layers <- function(layers, layer_names = NULL, context = "SpatRaster stack") {
  # Avoid do.call(c, layers): in some terra/R combinations it can dispatch to
  # base::c and return a plain list, causing terra::values(x=list) failures.
  layers <- layers[!vapply(layers, is.null, logical(1))]
  if (!length(layers)) stop("No raster layers supplied for ", context, ".", call. = FALSE)
  ok <- vapply(layers, inherits, logical(1), what = "SpatRaster")
  if (!all(ok)) {
    stop("Non-SpatRaster object supplied for ", context, ": layer index ",
         paste(which(!ok), collapse = ","), call. = FALSE)
  }
  out <- layers[[1]]
  if (length(layers) > 1L) {
    for (ii in 2:length(layers)) out <- c(out, layers[[ii]])
  }
  if (!is.null(layer_names)) {
    if (length(layer_names) != terra::nlyr(out)) {
      stop("Layer name count does not match stacked raster layer count for ", context, ".", call. = FALSE)
    }
    names(out) <- layer_names
  }
  out
}

build_target_template <- function() {
  f <- NA_character_
  for (mt in METRICS) {
    fs <- observed_annual_files(mt$stem, YEARS)
    ok <- fs[!is.na(fs) & file.exists(fs)]
    if (length(ok)) { f <- ok[1]; break }
  }
  if (is.na(f)) stop("No observed annual raster found; cannot build target-grid template.", call. = FALSE)
  r <- terra::rast(f)
  if (isTRUE(terra::is.lonlat(r))) stop("Observed annual raster appears to be lon/lat. v9 requires an equal-area/projected CRS before building a target-grid template.", call. = FALSE)
  rr <- terra::res(r)
  if (length(rr) < 2L || any(!is.finite(rr)) || rr[1] <= 0 || rr[2] <= 0) stop("Invalid observed raster resolution.", call. = FALSE)
  if (abs(rr[1] - rr[2]) > max(1e-6, 1e-6 * max(rr))) stop("Observed annual raster has non-square pixels; target-grid processing expects square cells.", call. = FALSE)
  fact_raw <- CELL_SIZE_M / rr[1]
  if (!is.finite(fact_raw) || fact_raw < 1 || abs(fact_raw - round(fact_raw)) > 1e-6) {
    stop(sprintf("CELL_SIZE_M=%s must be an integer multiple of observed raster resolution %.6f and cannot be finer than the source grid. Reproject/resample annual rasters before running.", CELL_SIZE_M, rr[1]), call. = FALSE)
  }
  fact <- max(1L, as.integer(round(fact_raw)))
  template <- if (fact > 1L) terra::aggregate(r, fact = fact, fun = mean, na.rm = TRUE) else r
  template[] <- NA
  names(template) <- paste0("template_", GRID_TAG)
  attr(template, "source_resolution") <- rr[1]
  attr(template, "source_factor") <- fact
  template
}

aggregate_mean_count_to_template <- function(file, template) {
  r <- terra::rast(file)
  if(terra::nlyr(r)!=1L)stop('Expected a single-layer annual raster: ',file,'; found ',terra::nlyr(r),' layers.',call.=FALSE)
  if(isTRUE(terra::is.lonlat(r)))stop('Input raster must use projected metric coordinates: ',file,call.=FALSE)
  same_crs<-terra::compareGeom(r,template,stopOnError=FALSE,crs=TRUE,ext=FALSE,rowcol=FALSE,res=FALSE)
  if(!same_crs)stop('Annual raster CRS differs from the analysis-grid template: ',file,'. Reproject before running; support counts must not be approximated.',call.=FALSE)
  g<-terra::global(r,c('min','max'),na.rm=TRUE);mn<-as.numeric(g[1,'min']);mx<-as.numeric(g[1,'max'])
  if(any(!is.finite(c(mn,mx))))stop('Raster contains no finite values or non-finite range: ',file,call.=FALSE)
  if(mn<INPUT_NEGATIVE_TOL)stop('Diversity/potential raster contains negative values below tolerance (min=',format(mn,digits=8),'): ',file,call.=FALSE)
  if(mn<0)r<-terra::ifel(r<0,0,r)
  fact_raw<-terra::res(template)[1]/terra::res(r)[1]
  if(!is.finite(fact_raw)||fact_raw<1||abs(fact_raw-round(fact_raw))>1e-6)stop(sprintf('Template resolution %.6f must be an integer multiple of input resolution %.6f for %s.',terra::res(template)[1],terra::res(r)[1],basename(file)),call.=FALSE)
  fact<-max(1L,as.integer(round(fact_raw)))
  mean_r<-if(fact>1L)terra::aggregate(r,fact=fact,fun=mean,na.rm=TRUE)else r
  valid<-terra::ifel(is.na(r),0,1)
  count_r<-if(fact>1L)terra::aggregate(valid,fact=fact,fun=sum,na.rm=TRUE)else valid
  # Geometric possible support is a raster, not fact^2: edge target cells can contain fewer source cells.
  possible_src<-terra::init(r,1)
  possible_r<-if(fact>1L)terra::aggregate(possible_src,fact=fact,fun=sum,na.rm=TRUE)else possible_src
  for(obj in list(mean_r,count_r,possible_r))if(!terra::compareGeom(obj,template,stopOnError=FALSE))stop('Aggregated raster does not align exactly with the analysis-grid template for ',basename(file),call.=FALSE)
  names(mean_r)<-'value';names(count_r)<-'n_valid_source';names(possible_r)<-'n_possible_source'
  list(mean=mean_r,count=count_r,possible=possible_r,
       possible_min=min(terra::values(possible_r,mat=FALSE),na.rm=TRUE),possible_max=max(terra::values(possible_r,mat=FALSE),na.rm=TRUE))
}

mk_rast_from_cells <- function(template, cells, vals, lyr_name = "layer") {
  r <- template
  terra::values(r) <- NA
  r[cells] <- vals
  names(r) <- lyr_name
  r
}

# ---- trend functions ---------------------------------------------------------
# Professional-package-first wrappers. Custom implementations remain available as
# fallbacks through v9_trend_package_utils.R. If that helper is unavailable for any
# reason, stop early rather than silently reverting to outdated local definitions.
if (!exists("v9_sen_slope") || !exists("v9_hamed_rao_mk") || !exists("v9_bky_adjust")) {
  stop("Missing v9_trend_package_utils.R; cannot initialise professional-package/fallback trend functions.", call. = FALSE)
}

bky_adjust_professional <- function(p, alpha_ref = Q_MAIN, precomputed_bh = NULL) {
  # pass a caller-computed BH q-vector so the fallback two-stage BKY
  # does not recompute p.adjust("BH"). Unchanged numerics and mutoss path.
  v9_bky_adjust(p, alpha_ref = alpha_ref, prefer_package = USE_PROFESSIONAL_TREND_PACKAGES,
                precomputed_bh = precomputed_bh)
}

# Keep FDR families consistent across BH/BY/BKY and with the pre-filtered trend table:
# only finite p-values are tested; NA rows are restored as NA afterwards. Base R
# p.adjust() otherwise counts NA values in its default n=length(p), which makes q-values
# depend on how many non-testable rows happened to remain in the table.
p_adjust_finite <- function(p, method) {
  p <- as.numeric(p)
  out <- rep(NA_real_, length(p))
  ok <- is.finite(p)
  if (any(ok)) out[ok] <- stats::p.adjust(p[ok], method = method)
  out
}

compute_metric_1km_series_and_trend <- function(metric, template, region_id = NULL, inner_workers = NULL) {
  if (is.null(inner_workers)) inner_workers <- N_WORKERS
  files <- observed_annual_files(metric$stem, YEARS)
  ok_files <- !is.na(files) & file.exists(files)
  if(any(!ok_files)&&!isTRUE(ALLOW_PARTIAL_OBS_YEARS))stop("ALLOW_PARTIAL_OBS_YEARS=FALSE but missing annual files for ",metric$stem,": ",paste(YEARS[!ok_files],collapse=","),call.=FALSE)
  if (sum(ok_files) < ceiling(length(YEARS) * OBS_YEARS_MIN_FRAC)) {
    stop("Too few annual files for ", metric$stem, ": ", sum(ok_files), "/", length(YEARS), call. = FALSE)
  }
  if (any(!ok_files)) warning(metric$stem, ": missing annual years: ", paste(YEARS[!ok_files], collapse = ","), call. = FALSE)

  value_layers <- list(); vf_layers <- list()
  count_layers <- if (WRITE_ANNUAL_TABLES) list() else NULL
  possible_n <- NA_real_
  for (yy in YEARS) {
    f <- files[as.character(yy)]
    if (!is.na(f) && file.exists(f)) {
      ag <- aggregate_mean_count_to_template(f, template)
      possible_n <- as.numeric(terra::values(ag$possible,mat=FALSE))
      val <- ag$mean
      cnt <- ag$count
      vf <- cnt / ag$possible
      val[vf < MIN_VALID_FRAC] <- NA
    } else {
      val <- template; val[] <- NA
      cnt <- template; cnt[] <- NA
      vf <- template; vf[] <- NA
    }
    names(val) <- paste0("value_", yy)
    names(cnt) <- paste0("n_valid_", yy)
    names(vf) <- paste0("valid_frac_", yy)
    value_layers[[as.character(yy)]] <- val
    if (WRITE_ANNUAL_TABLES) count_layers[[as.character(yy)]] <- cnt
    vf_layers[[as.character(yy)]] <- vf
  }
  value_stack <- stack_spatraster_layers(value_layers, context = paste0(metric$stem, " annual value layers"))
  count_stack <- if (WRITE_ANNUAL_TABLES) stack_spatraster_layers(count_layers, context = paste0(metric$stem, " annual count layers")) else NULL
  vf_stack <- stack_spatraster_layers(vf_layers, context = paste0(metric$stem, " annual valid-frac layers"))
  names(value_stack) <- paste0(metric$pref, "_", YEARS)
  names(vf_stack) <- paste0(metric$pref, "_valid_frac_", YEARS)

  vals <- terra::values(value_stack, mat = TRUE)
  vfs <- terra::values(vf_stack, mat = TRUE)
  cnts <- if (WRITE_ANNUAL_TABLES) terra::values(count_stack, mat = TRUE) else NULL
  cells <- seq_len(terra::ncell(template))
  xy <- terra::xyFromCell(template, cells)
  region_vals <- if (!is.null(region_id)) as.integer(terra::values(region_id)[cells]) else rep(NA_integer_, length(cells))
  roman_vals <- c("I", "II", "III", "IV", "V")[region_vals]
  roman_vals[is.na(roman_vals)] <- "outside"

  n_valid_years <- rowSums(is.finite(vals))

  # Optional long table of the actual analysis-grid annual series used for Sen/MK.
  # It is large and not required by the Fig.6 workflow, so it is built only
  # when WRITE_ANNUAL_TABLES=TRUE.
  annual_df <- if (WRITE_ANNUAL_TABLES) {
    dplyr::bind_rows(lapply(seq_along(YEARS), function(j) {
      data.frame(
        cell = cells,
        x = xy[, 1],
        y = xy[, 2],
        metric = metric$stem,
        metric_label = metric$label,
        year = YEARS[j],
        value_target_grid = as.numeric(vals[, j]),
        n_valid_source = as.numeric(cnts[, j]),
        possible_n_source = possible_n,
        valid_frac = as.numeric(vfs[, j]),
        region_id = region_vals,
        roman = roman_vals,
        mountain_zone = roman_vals,
        stringsAsFactors = FALSE
      )
    }))
  } else data.frame()


  trend_year_idx <- lapply(TREND_PERIODS[RUN_TREND_PERIODS], function(tp) which(YEARS %in% tp$years))

  # ---- parallel processing -----------------------------------------
  # Most 5 km template cells (ocean / outside the study mask) are all-NA and can
  # never yield a trend in any run period. Computing -- and serialising -- an NA
  # row for each of them dominates the cost of the original per-cell dispatch.
  # classify_period() left-joins trend results onto the cell table by `cell`, so a
  # cell that is simply ABSENT from trend_df receives identical NA trend columns to
  # one that produced an explicit NA row. We can therefore drop hopeless cells from
  # trend fitting with no change in results. The annual long table (annual_df above)
  # still retains every cell for reproducibility / temporal-correction.
  union_trend_idx <- sort(unique(unlist(trend_year_idx, use.names = FALSE)))
  global_min_years <- min(vapply(RUN_TREND_PERIODS,
                                 function(tp) trend_min_valid_years(length(trend_year_idx[[tp]])),
                                 integer(1)))
  cell_finite_union <- rowSums(is.finite(vals[, union_trend_idx, drop = FALSE]))
  keep_cells <- which(cell_finite_union >= global_min_years)
  if (!length(keep_cells)) keep_cells <- integer(0)
  # subset the heavy per-cell objects to candidate cells BEFORE the parallel
  # export. Only keep_cells are ever read by the worker, but the original code exported
  # the full ncell-row vals/vfs matrices (plus cells/xy) to every PSOCK worker -- often
  # ~5x larger than needed. We slice them here and reindex trend_jobs$cell_i to the
  # SUBSET (1..length(keep_cells)); cells_k preserves the true cell IDs for output.
  total_template_cells <- nrow(vals)
  vals_k  <- vals[keep_cells, , drop = FALSE]
  vfs_k   <- vfs[keep_cells, , drop = FALSE]
  cells_k <- cells[keep_cells]
  xy_k    <- xy[keep_cells, , drop = FALSE]

  # After subsetting to trend candidates, the full ncell x years matrices and
  # raster stacks are no longer used unless the optional annual table is enabled.
  # Releasing them is important for 1-km runs and for metric-parallel execution.
  if (!WRITE_ANNUAL_TABLES) {
    rm(vals, vfs, cnts, value_stack, vf_stack, count_stack, value_layers, vf_layers, count_layers)
    gc(verbose = FALSE)
  }
  trend_jobs <- expand.grid(cell_i = seq_along(keep_cells), trend_period = RUN_TREND_PERIODS, stringsAsFactors = FALSE)
  msg("  %s: fitting trends on %d/%d candidate cells x %d trend period(s) = %d tasks",
      metric$label, length(keep_cells), total_template_cells, length(RUN_TREND_PERIODS), nrow(trend_jobs))

  # Chunk worker: process a contiguous block of job indices and assemble ONE
  # data.frame from pre-allocated column vectors. This removes both the per-task
  # PSOCK dispatch overhead and the per-row data.frame + bind_rows anti-pattern
  # (tens of thousands of 1-row frames) of the original design.
  trend_chunk_worker <- function(job_idx) {
    m <- length(job_idx)
    cell_v <- integer(m); x_v <- numeric(m); y_v <- numeric(m)
    tp_v <- character(m); tpl_v <- character(m); ty_v <- character(m)
    n_years_v <- numeric(m); sen_v <- numeric(m); mkpraw_v <- numeric(m); phr_v <- numeric(m)
    S_v <- numeric(m); varS_v <- numeric(m); varShr_v <- numeric(m); acf1_v <- numeric(m)
    mksrc_v <- character(m); sensrc_v <- character(m); mvf_v <- numeric(m)
    for (k in seq_len(m)) {
      job_i <- job_idx[k]
      i <- trend_jobs$cell_i[job_i]          # row index into the keep_cells subset
      tp_name <- trend_jobs$trend_period[job_i]
      idx <- trend_year_idx[[tp_name]]
      yy <- YEARS[idx]
      y <- as.numeric(vals_k[i, idx])
      min_years <- trend_min_valid_years(length(yy))
      cell_v[k] <- cells_k[i]; x_v[k] <- xy_k[i, 1]; y_v[k] <- xy_k[i, 2]
      tp_v[k] <- tp_name; tpl_v[k] <- TREND_PERIODS[[tp_name]]$label
      ty_v[k] <- paste(yy, collapse = ",")
      mvf_v[k] <- mean(as.numeric(vfs_k[i, idx]), na.rm = TRUE)
      if (sum(is.finite(y)) < min_years) {
        n_years_v[k] <- sum(is.finite(y))
        sen_v[k] <- NA_real_; mkpraw_v[k] <- NA_real_; phr_v[k] <- NA_real_
        S_v[k] <- NA_real_; varS_v[k] <- NA_real_; varShr_v[k] <- NA_real_; acf1_v[k] <- NA_real_
        mksrc_v[k] <- NA_character_; sensrc_v[k] <- NA_character_
        next
      }
      stats_obj <- v9_trend_stats(yy,y,lag_max=HR_LAG_MAX,sig_only=HR_USE_SIGNIFICANT_ACF_ONLY,
        prefer_package=USE_PROFESSIONAL_TREND_PACKAGES,prefer_sen_package=(USE_PROFESSIONAL_TREND_PACKAGES&&USE_TREND_PACKAGE_SEN))
      sen_obj <- stats_obj$sen; mk_raw <- stats_obj$raw; mk <- stats_obj$hr
      n_years_v[k] <- as.numeric(unname(mk["n"])); sen_v[k] <- sen_obj$value
      mkpraw_v[k] <- as.numeric(unname(mk_raw["p"])); phr_v[k] <- as.numeric(unname(mk["p"]))
      S_v[k] <- as.numeric(unname(mk["S"])); varS_v[k] <- as.numeric(unname(mk["varS"]))
      varShr_v[k] <- as.numeric(unname(mk["varS_hr"])); acf1_v[k] <- as.numeric(unname(mk["acf1"]))
      mksrc_v[k] <- as.character(unname(mk["source"])); sensrc_v[k] <- sen_obj$source
    }
    data.frame(cell = cell_v, x = x_v, y = y_v,
               metric = rep(metric$stem, m), metric_label = rep(metric$label, m),
               trend_period = tp_v, trend_period_label = tpl_v, trend_years = ty_v,
               n_years = n_years_v, sen_slope = sen_v, mk_p_raw = mkpraw_v, p_hr = phr_v,
               S = S_v, varS = varS_v, varS_hr = varShr_v, acf1 = acf1_v,
               mk_source = mksrc_v, sen_source = sensrc_v, mean_valid_frac = mvf_v,
               stringsAsFactors = FALSE)
  }

  if (nrow(trend_jobs) == 0L) {
    trend_df <- trend_chunk_worker(integer(0))
  } else if (inner_workers > 1L && exists('v9_cluster_lapply')) {
    # 1-km FAST mode: use many smaller chunks and load-balanced scheduling.
    # Static parLapply with only workers*4 chunks often leaves one long "tail"
    # worker running after the other 7 workers have finished.  The two knobs below
    # can be overridden from .bat/PowerShell if needed.
    trend_chunk_factor <- env_int("TREND_CHUNK_FACTOR", 24L)
    trend_load_balance <- env_logical("TREND_LOAD_BALANCE", TRUE)
    if (!is.finite(trend_chunk_factor) || trend_chunk_factor < 4L) trend_chunk_factor <- 4L
    job_chunks <- if (exists('v9_chunk_indices')) {
      v9_chunk_indices(nrow(trend_jobs), workers = inner_workers, chunk_factor = trend_chunk_factor)
    } else list(seq_len(nrow(trend_jobs)))
    chunk_results <- v9_cluster_lapply(job_chunks, trend_chunk_worker,
      workers = inner_workers,
      export = c('vals_k','vfs_k','cells_k','xy_k','metric','YEARS','HR_LAG_MAX','HR_USE_SIGNIFICANT_ACF_ONLY',
                 'USE_PROFESSIONAL_TREND_PACKAGES','USE_TREND_PACKAGE_SEN','USE_MODIFIEDMK_PACKAGE','V9_USE_PROFESSIONAL_TREND_PACKAGES','V9_USE_TREND_PACKAGE_SEN','V9_USE_MODIFIEDMK_PACKAGE','TREND_PERIODS','RUN_TREND_PERIODS','TREND_MIN_VALID_YEAR_FRAC','trend_min_valid_years','trend_year_idx','trend_jobs',
                 'v9_pair_indices','.v9_pair_cache','v9_trend_stats','v9_sen_slope','v9_sen_slope_fallback','v9_hamed_rao_mk','v9_hamed_rao_mk_fallback','v9_raw_mk_p','v9_mk_basic_stats','v9_mk_p_from_S','v9_extract_numeric_by_name','v9_extract_sens_slope_estimate','v9_extract_mmkh_p','v9_pkg_available','v9_is_unit_consecutive'),
      export_env = environment(), task_label = paste0(GRID_LABEL, ' trend-period rows: ', metric$label),
      load_balance = trend_load_balance)
    trend_df <- dplyr::bind_rows(chunk_results)
  } else {
    # Serial trend fit. Used when inner_workers<=1, e.g. when the two METRICS are
    # themselves being run on separate workers (no nested PSOCK clusters).
    trend_df <- trend_chunk_worker(seq_len(nrow(trend_jobs)))
  }
  # NOTE: the per-year SpatRaster stacks (value/count/valid_frac) are intentionally
  # NOT returned: nothing downstream of the per-metric step uses them, and excluding
  # them keeps the result a plain data.frame list (so the metric step can run on a
  # PSOCK worker, whose results must serialise back) and avoids holding 3 extra
  # ncell x nyears rasters per metric in memory.
  list(trend = trend_df, annual = annual_df)
}

apply_fdr_by_metric <- function(trend_df) {
  trend_df %>%
    dplyr::group_by(metric, trend_period) %>%
    dplyr::group_modify(function(.x, .y) {
      # compute BH once and reuse it for the two-stage BKY fallback.
      q_bh_vec <- p_adjust_finite(.x$p_hr, method = "BH")
      q_bky_vec <- bky_adjust_professional(.x$p_hr, alpha_ref = Q_MAIN, precomputed_bh = q_bh_vec)
      .x$q_bh <- q_bh_vec
      .x$q_by <- p_adjust_finite(.x$p_hr, method = "BY")
      .x$q_bky <- q_bky_vec
      .x$bky_source <- attr(q_bky_vec, "source")
      .x
    }) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      trend_state = dplyr::case_when(
        is.finite(sen_slope) & is.finite(q_bh) & sen_slope > 0 & q_bh < Q_MAIN ~ "T+",
        is.finite(sen_slope) & is.finite(q_bh) & sen_slope < 0 & q_bh < Q_MAIN ~ "T-",
        is.finite(q_bh) ~ "T0",
        TRUE ~ NA_character_
      ),
      trend_state_q010 = dplyr::case_when(
        is.finite(sen_slope) & is.finite(q_bh) & sen_slope > 0 & q_bh < Q_WILKS_SENS ~ "T+",
        is.finite(sen_slope) & is.finite(q_bh) & sen_slope < 0 & q_bh < Q_WILKS_SENS ~ "T-",
        is.finite(q_bh) ~ "T0",
        TRUE ~ NA_character_
      )
    ) %>%
    dplyr::ungroup()
}

potential_ceiling_1km_stable <- function(metric, template) {
  # MS 2.6 / SI S5.4: exactly ONE static, age-100 q95 surface per response.
  f <- file.path(POTENTIAL_DIR, sprintf("Q95_%s_age100_1km.tif", metric$stem))
  if (!file.exists(f)) stop("Missing static age-100 q95: ", f)
  z <- terra::rast(f)
  if (terra::nlyr(z) != 1L || !terra::compareGeom(z, template, stopOnError=FALSE))
    stop("Static potential must be a single layer on the common 1-km template: ", f)
  z[z <= 0] <- NA
  n <- terra::ifel(is.finite(z), 1, 0)
  names(z) <- paste0(metric$pref, "_potential_q95_age100")
  names(n) <- paste0(metric$pref, "_potential_n_valid_surfaces")
  list(mean=z, n_valid_years=n, mean_valid_frac=n,
       diagnostics=data.frame(metric=metric$stem, reference_age=100, quantile=0.95,
         source=normalizePath(f), method_used="static_q95_age100_common_1km",
         potential_file_exists=TRUE,
         final_target_grid_finite=sum(is.finite(terra::values(z,mat=FALSE)))))
}

build_stage_mean_1km <- function(metric, years, template) {
  files <- observed_annual_files(metric$stem, years)
  ok <- !is.na(files) & file.exists(files)
  min_n <- max(1L, ceiling(length(years) * OBS_YEARS_MIN_FRAC))
  if (sum(ok) < min_n) stop("Too few observed annual rasters for stage mean ", metric$stem, ": ", sum(ok), "/", length(years), call. = FALSE)
  if (any(!ok)) warning("Stage mean ", metric$stem, " missing years: ", paste(years[!ok], collapse = ","), call. = FALSE)

  vals_list <- list(); vf_list <- list(); cnt_list <- list(); possible_n <- NA_real_
  for (yy in years) {
    f <- files[as.character(yy)]
    if (!is.na(f) && file.exists(f)) {
      ag <- aggregate_mean_count_to_template(f, template)
      possible_n <- as.numeric(terra::values(ag$possible,mat=FALSE))
      val <- ag$mean
      cnt <- ag$count
      vf <- cnt / ag$possible
      val[vf < MIN_VALID_FRAC] <- NA
    } else {
      val <- template; val[] <- NA
      cnt <- template; cnt[] <- NA
      vf <- template; vf[] <- NA
    }
    vals_list[[as.character(yy)]] <- val
    vf_list[[as.character(yy)]] <- vf
    cnt_list[[as.character(yy)]] <- cnt
  }
  value_stack <- stack_spatraster_layers(vals_list, context = paste0(metric$stem, " stage mean value layers"))
  vf_stack <- stack_spatraster_layers(vf_list, context = paste0(metric$stem, " stage mean valid-frac layers"))
  vals <- terra::values(value_stack, mat = TRUE)
  vfs <- terra::values(vf_stack, mat = TRUE)
  n_valid_years <- rowSums(is.finite(vals))
  mean_valid_frac <- rowMeans(vfs, na.rm = TRUE)
  stage_mean <- rowMeans(vals, na.rm = TRUE)
  stage_mean[!is.finite(stage_mean)] <- NA_real_
  stage_mean[n_valid_years < min_n] <- NA_real_
  mean_valid_frac[n_valid_years < min_n] <- NA_real_

  r_mean <- template; terra::values(r_mean) <- stage_mean; names(r_mean) <- paste0(metric$pref, "_stage_mean")
  r_n <- template; terra::values(r_n) <- n_valid_years; names(r_n) <- paste0(metric$pref, "_stage_n_valid_years")
  r_vf <- template; terra::values(r_vf) <- mean_valid_frac; names(r_vf) <- paste0(metric$pref, "_stage_mean_valid_frac")
  list(mean = r_mean, n_valid_years = r_n, mean_valid_frac = r_vf, min_stage_valid_years = min_n)
}

summarise_class_area <- function(df, metric_pref, period_name, trend_period_name, analysis_group) {
  class_col <- paste0(metric_pref, "_class6")
  metric_label <- if (metric_pref == "rich") "Richness" else "Shannon"
  d <- df[!is.na(df[[class_col]]), , drop = FALSE]
  if (!"area_km2" %in% names(d)) d$area_km2 <- prod(terra::res(template)) / 1e6
  summarise_zone <- function(dd) {
    if (!nrow(dd)) return(data.frame(class6 = character(), n_cells = integer(), area_km2 = numeric(), stringsAsFactors = FALSE))
    dd %>%
      dplyr::group_by(class6 = .data[[class_col]]) %>%
      dplyr::summarise(n_cells = dplyr::n(), area_km2 = sum(.data$area_km2, na.rm = TRUE), .groups = "drop") %>%
      dplyr::mutate(class6 = as.character(.data$class6), n_cells = as.integer(.data$n_cells), area_km2 = as.numeric(.data$area_km2))
  }
  overall <- summarise_zone(d) %>%
    dplyr::mutate(analysis_group = analysis_group, period = period_name, trend_period = trend_period_name, metric = metric_pref, metric_label = metric_label, zone = "Overall")
  zone <- d %>%
    dplyr::mutate(zone = zone_clean(.data$roman)) %>%
    dplyr::group_by(zone, class6 = .data[[class_col]]) %>%
    dplyr::summarise(n_cells = dplyr::n(), area_km2 = sum(.data$area_km2, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(zone = as.character(.data$zone), class6 = as.character(.data$class6), n_cells = as.integer(.data$n_cells), area_km2 = as.numeric(.data$area_km2), analysis_group = analysis_group, period = period_name, trend_period = trend_period_name, metric = metric_pref, metric_label = metric_label)
  out <- dplyr::bind_rows(overall, zone) %>%
    dplyr::mutate(zone = as.character(.data$zone), class6 = as.character(.data$class6)) %>%
    tidyr::complete(analysis_group, period, trend_period, metric, metric_label, zone = REG_ORDER, class6 = CLASS6_LEVELS, fill = list(n_cells = 0L, area_km2 = 0)) %>%
    dplyr::mutate(zone = as.character(.data$zone), class6 = as.character(.data$class6), n_cells = as.integer(.data$n_cells), area_km2 = as.numeric(.data$area_km2)) %>%
    dplyr::group_by(analysis_group, period, trend_period, metric, metric_label, zone) %>%
    dplyr::mutate(prop = if (sum(.data$area_km2, na.rm = TRUE) > 0) .data$area_km2 / sum(.data$area_km2, na.rm = TRUE) else NA_real_, pct = 100 * prop) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(class_label = CLASS6_LABEL[class6])
  out
}

write_zero_classification_diagnostics <- function(df, period_root, period_name, trend_period_name, analysis_group) {
  diag <- data.frame(
    analysis_group = analysis_group,
    period = period_name,
    trend_period = trend_period_name,
    total_cells = nrow(df),
    rich_obs_finite = sum(is.finite(df$rich_obs)),
    shannon_obs_finite = sum(is.finite(df$shannon_obs)),
    rich_potential_finite = sum(is.finite(df$rich_potential)),
    shannon_potential_finite = sum(is.finite(df$shannon_potential)),
    rich_RF_finite = sum(is.finite(df$rich_RF)),
    shannon_RF_finite = sum(is.finite(df$shannon_RF)),
    rich_trend_state_nonNA = sum(!is.na(df$rich_trend_state)),
    shannon_trend_state_nonNA = sum(!is.na(df$shannon_trend_state)),
    rich_class6_nonNA = sum(!is.na(df$rich_class6)),
    shannon_class6_nonNA = sum(!is.na(df$shannon_class6)),
    stringsAsFactors = FALSE
  )
  write.csv(diag, file.path(period_root, "RF_6class_zero_classification_diagnostics.csv"), row.names = FALSE)
  diag
}

classify_period <- function(period_name, cfg, trend_period_name, trend_cfg, analysis_group, template, region_id, trend_wide, pot_rich_obj, pot_shannon_obj) {
  msg("[%s | %s] Building RF and six classes ...", period_name, trend_period_name)
  period_root <- file.path(OUT_ROOT, analysis_group)
  dir.create(period_root, recursive = TRUE, showWarnings = FALSE)

  rich_stage <- build_stage_mean_1km(METRICS[[1]], cfg$obs_years, template)
  shan_stage <- build_stage_mean_1km(METRICS[[2]], cfg$obs_years, template)
  rich_obs <- rich_stage$mean
  shan_obs <- shan_stage$mean
  pot_rich <- pot_rich_obj$mean
  pot_shannon <- pot_shannon_obj$mean

  S <- c(rich_obs, shan_obs, pot_rich, pot_shannon, region_id,
         rich_stage$n_valid_years, rich_stage$mean_valid_frac,
         shan_stage$n_valid_years, shan_stage$mean_valid_frac,
         pot_rich_obj$n_valid_years, pot_rich_obj$mean_valid_frac,
         pot_shannon_obj$n_valid_years, pot_shannon_obj$mean_valid_frac)
  names(S) <- c("rich_obs", "shannon_obs", "rich_potential", "shannon_potential", "region_id",
                "rich_stage_n_valid_years", "rich_stage_mean_valid_frac",
                "shannon_stage_n_valid_years", "shannon_stage_mean_valid_frac",
                "rich_potential_n_valid_surfaces", "rich_potential_mean_valid_frac",
                "shannon_potential_n_valid_surfaces", "shannon_potential_mean_valid_frac")
  # Keep rows that contain at least one non-NA layer.  This preserves all potentially
  # classifiable cells, while dropping all-empty background cells from the 1-km template.
  df <- as.data.frame(S, cells = TRUE, xy = TRUE, na.rm = NA)
  df$analysis_group <- analysis_group
  df$period <- period_name
  df$rf_period <- period_name
  df$stage_label <- cfg$stage_label
  df$obs_years <- paste(cfg$obs_years, collapse = ",")
  df$potential_reference_age <- 100L
  df$potential_quantile <- 0.95
  df$rf_formula <- "mean_obs_2016_2020_divided_by_static_q95_age100"
  df$trend_period <- trend_period_name
  df$trend_period_label <- trend_cfg$label
  df$trend_years <- paste(trend_cfg$years, collapse = ",")
  df$roman <- c("I", "II", "III", "IV", "V")[df$region_id]
  df$area_km2 <- prod(terra::res(template)) / 1e6

  # RF denominator guard: the previous pmax(potential, MIN_DENOM) floor
  # turned a vanishing/zero potential into obs/1e-6, i.e. a RF of ~1e6 that ALWAYS
  # crossed RF_THRESH and spuriously classified the cell as "High RF". A potential
  # at or below MIN_DENOM carries no usable denominator signal, so RF is left
  # undefined (NA) there. Cells with a genuinely positive potential (potential >>
  # MIN_DENOM, i.e. essentially all real cells) are unchanged.
  rich_den <- df$rich_potential
  rich_den[!is.finite(rich_den) | rich_den <= MIN_DENOM] <- NA_real_
  shan_den <- df$shannon_potential
  shan_den[!is.finite(shan_den) | shan_den <= MIN_DENOM] <- NA_real_
  df$rich_RF <- df$rich_obs / rich_den
  df$shannon_RF <- df$shannon_obs / shan_den
  df$rich_gap <- pmax(0, df$rich_potential - df$rich_obs)
  df$shannon_gap <- pmax(0, df$shannon_potential - df$shannon_obs)

  trend_join <- trend_wide %>% dplyr::select(-dplyr::any_of("trend_period"))
  df <- df %>% dplyr::left_join(trend_join, by = "cell", suffix = c("", "_trend"))

  add_class <- function(d, pref) {
    rf <- d[[paste0(pref, "_RF")]]
    sen <- d[[paste0(pref, "_sen_slope")]]
    q <- d[[paste0(pref, "_q_bh")]]
    trend_state <- d[[paste0(pref, "_trend_state")]]
    rf_state <- ifelse(is.finite(rf) & rf >= RF_THRESH, "H", ifelse(is.finite(rf) & rf < RF_THRESH, "L", NA_character_))
    class6 <- dplyr::case_when(
      rf_state == "H" & trend_state == "T+" ~ "H+",
      rf_state == "H" & trend_state == "T0" ~ "H0",
      rf_state == "H" & trend_state == "T-" ~ "H-",
      rf_state == "L" & trend_state == "T+" ~ "L+",
      rf_state == "L" & trend_state == "T0" ~ "L0",
      rf_state == "L" & trend_state == "T-" ~ "L-",
      TRUE ~ NA_character_
    )
    d[[paste0(pref, "_rf_state")]] <- rf_state
    d[[paste0(pref, "_class6")]] <- class6
    d[[paste0(pref, "_class6_code")]] <- unname(CLASS6_CODE[class6])
    d[[paste0(pref, "_trend_code")]] <- unname(TREND_CODE[trend_state])
    d
  }
  df <- add_class(df, "rich")
  df <- add_class(df, "shannon")
  df$cross_trend_conflict <- as.integer((df$rich_trend_code * df$shannon_trend_code) < 0)

  # Keep only cells with enough valid information for at least one class; retain full columns.
  # Reporting-domain cells must be assigned to the five mountain regions.
  df_out <- df %>% dplyr::filter(roman %in% c("I","II","III","IV","V"),
                               !is.na(rich_class6) | !is.na(shannon_class6))
  if (!nrow(df_out)) {
    diag <- write_zero_classification_diagnostics(df, period_root, period_name, trend_period_name, analysis_group)
    stop("No classified cells were produced for ", analysis_group, ". Diagnostics written to ",
         file.path(period_root, "RF_6class_zero_classification_diagnostics.csv"),
         ". Key counts: rich_RF_finite=", diag$rich_RF_finite,
         ", shannon_RF_finite=", diag$shannon_RF_finite,
         ", rich_trend_state_nonNA=", diag$rich_trend_state_nonNA,
         ", shannon_trend_state_nonNA=", diag$shannon_trend_state_nonNA,
         call. = FALSE)
  }

  msg("[%s | %s] Writing RF class cell table ...", period_name, trend_period_name)
  write.csv(df_out, file.path(period_root, "RF_6class_cell_table.csv"), row.names = FALSE)
  msg("[%s | %s] Writing RF class area summary ...", period_name, trend_period_name)
  write.csv(dplyr::bind_rows(summarise_class_area(df_out, "rich", period_name, trend_period_name, analysis_group),
                             summarise_class_area(df_out, "shannon", period_name, trend_period_name, analysis_group)),
            file.path(period_root, "RF_6class_area_summary.csv"), row.names = FALSE)

  if (WRITE_RASTERS) {
    msg("[%s | %s] Building output raster stack ...", period_name, trend_period_name)
    out_stack <- c(
      mk_rast_from_cells(template, df_out$cell, df_out$rich_class6_code, "rich_class6_code"),
      mk_rast_from_cells(template, df_out$cell, df_out$shannon_class6_code, "shannon_class6_code"),
      mk_rast_from_cells(template, df_out$cell, df_out$rich_trend_code, "rich_trend_code"),
      mk_rast_from_cells(template, df_out$cell, df_out$shannon_trend_code, "shannon_trend_code"),
      mk_rast_from_cells(template, df_out$cell, df_out$rich_RF, "rich_RF"),
      mk_rast_from_cells(template, df_out$cell, df_out$shannon_RF, "shannon_RF"),
      region_id
    )
    names(out_stack) <- c("rich_class6_code", "shannon_class6_code", "rich_trend_code", "shannon_trend_code", "rich_RF", "shannon_RF", "region_id")
    msg("[%s | %s] Writing output raster: %s ...", period_name, trend_period_name, MAP_FILE)
    terra::writeRaster(out_stack, file.path(period_root, MAP_FILE), overwrite = TRUE, datatype = "FLT4S")
  }

  msg("[%s | %s] Writing RF class bundle RDS without compression ...", period_name, trend_period_name)
  saveRDS(list(params = list(analysis_group = analysis_group, period = period_name, trend_period = trend_period_name, RF_formula = "mean_obs_2016_2020 / static_q95_age100", potential_reference_age = 100L, potential_quantile = 0.95, RF_THRESH = RF_THRESH, Q_MAIN = Q_MAIN, MIN_VALID_FRAC = MIN_VALID_FRAC, TREND_MIN_VALID_YEAR_FRAC = TREND_MIN_VALID_YEAR_FRAC),
               cell_table = df_out), file.path(period_root, "RF_6class_bundle.rds"), compress = FALSE)
  msg("[%s | %s] Done. classified cells=%d", period_name, trend_period_name, nrow(df_out))
  invisible(df_out)
}

# ---- main --------------------------------------------------------------------
msg("RF_6class FEM %s started in: %s", FEM_SOFTWARE_VERSION, DATA_DIR_ABS)
msg("Output root: %s", OUT_ROOT)
msg("Main choices: analysis=%s; trend=Sen + Hamed-Rao MK; FDR=BH q<%.3f; RF_THRESH=%.2f", GRID_LABEL, Q_MAIN, RF_THRESH)
if (exists('v9_log_parallel_settings')) v9_log_parallel_settings(n_tasks = NA_integer_, requested = N_WORKERS, label = 'main pipeline')


template <- build_target_template()
terra::writeRaster(template, file.path(OUT_ROOT, TEMPLATE_FILE), overwrite = TRUE)
msg("Template built: ncell=%d, res=%.1f x %.1f", terra::ncell(template), terra::res(template)[1], terra::res(template)[2])

msg("Rasterising mountain regions ...")
mountains <- prepare_mountains(MOUNTAIN_SHP, template)
region_id <- terra::rasterize(mountains, template, field = "rid")
# Persist region_id so metric workers can rebuild it from disk (SpatRaster objects
# cannot be serialised across PSOCK workers).
region_id_path <- file.path(OUT_ROOT, REGION_FILE)
terra::writeRaster(region_id, region_id_path, overwrite = TRUE, datatype = "INT2S")

msg("Computing %s annual series and grid-level trends ...", GRID_LABEL)
# Metric-level parallelism: the per-metric work is dominated by SERIAL terra
# raster aggregation + value extraction + annual-table building; the inner Sen/MK fit
# is only a short parallel burst. The two metrics (Richness, Shannon) are independent,
# so we run them on separate PSOCK workers (each running its inner trend fit SERIALLY to
# avoid nested clusters). Workers rebuild template/region_id from disk and return only
# data.frames. If anything in the parallel path fails, we fall back to the proven serial
# path so the pipeline never breaks.
metric_outer_workers <- if (isTRUE(MAIN_METRIC_PARALLEL)) v9_resolve_workers(length(METRICS), N_WORKERS) else 1L
trend_objects <- NULL
if (metric_outer_workers > 1L && exists('v9_cluster_lapply')) {
  metric_worker <- function(metric) {
    template_w <- terra::rast(file.path(OUT_ROOT, TEMPLATE_FILE))
    region_w   <- terra::rast(file.path(OUT_ROOT, REGION_FILE))
    res <- compute_metric_1km_series_and_trend(metric, template_w, region_w, inner_workers = 1L)
    list(trend = res$trend, annual = res$annual)   # data.frames only -> serialise cleanly
  }
  trend_objects <- tryCatch(
    v9_cluster_lapply(METRICS, metric_worker,
      workers = metric_outer_workers,
      packages = c('terra','dplyr','tidyr'),
      export = c('OUT_ROOT','YEARS','MIN_VALID_FRAC','TREND_MIN_VALID_YEAR_FRAC',
                 'ALLOW_PARTIAL_OBS_YEARS','OBS_YEARS_MIN_FRAC','POTENTIAL_DIR','DATA_DIR_ABS',
                 'CELL_SIZE_M','GRID_TAG','GRID_LABEL','TEMPLATE_FILE','REGION_FILE','MAP_FILE','HR_LAG_MAX','HR_USE_SIGNIFICANT_ACF_ONLY','USE_PROFESSIONAL_TREND_PACKAGES',
                 'USE_TREND_PACKAGE_SEN','USE_MODIFIEDMK_PACKAGE','V9_USE_PROFESSIONAL_TREND_PACKAGES',
                 'V9_USE_TREND_PACKAGE_SEN','V9_USE_MODIFIEDMK_PACKAGE','METRICS','TREND_PERIODS','RUN_TREND_PERIODS',
                 'N_WORKERS','MIN_DENOM','WRITE_ANNUAL_TABLES','ALLOW_AMBIGUOUS_INPUT_FILES','INPUT_NEGATIVE_TOL',
                 'compute_metric_1km_series_and_trend','observed_annual_files','find_observed_annual_file',
                 'observed_annual_search_dirs','is_safe_observed_annual_candidate','aggregate_mean_count_to_template',
                 'stack_spatraster_layers','msg','make_abs_path','env_str','trend_min_valid_years',
                 'v9_chunk_indices','v9_cluster_lapply','v9_resolve_workers',
                 'v9_pair_indices','.v9_pair_cache','v9_trend_stats','v9_sen_slope','v9_sen_slope_fallback','v9_hamed_rao_mk','v9_hamed_rao_mk_fallback','v9_raw_mk_p',
                 'v9_mk_basic_stats','v9_mk_p_from_S','v9_extract_numeric_by_name','v9_extract_sens_slope_estimate',
                 'v9_extract_mmkh_p','v9_pkg_available','v9_is_unit_consecutive'),
      export_env = environment(),
      task_label = paste0('per-metric compute (', length(METRICS), ' metrics)')),
    error = function(e) { warning('Metric-level parallel path failed (', conditionMessage(e), '); falling back to serial per-metric compute.', call. = FALSE); NULL }
  )
  # Validate the parallel result; if malformed, discard and fall back.
  if (!is.null(trend_objects)) {
    ok_par <- is.list(trend_objects) && length(trend_objects) == length(METRICS) &&
      all(vapply(trend_objects, function(o) is.list(o) && all(c('trend','annual') %in% names(o)) && is.data.frame(o$trend) && is.data.frame(o$annual), logical(1)))
    if (!ok_par) { warning('Metric-level parallel result was malformed; falling back to serial per-metric compute.', call. = FALSE); trend_objects <- NULL }
  }
}
if (is.null(trend_objects)) {
  # Serial per-metric path (each metric still parallelises its own inner trend fit).
  trend_objects <- lapply(METRICS, compute_metric_1km_series_and_trend, template = template, region_id = region_id)
}
names(trend_objects) <- vapply(METRICS, function(x) x$pref, character(1))
msg("Merging annual metric series ...")
annual_1km_df <- dplyr::bind_rows(lapply(trend_objects, `[[`, "annual"))
msg("Annual metric rows merged: %s", nrow(annual_1km_df))
if (WRITE_ANNUAL_TABLES) {
  msg("Writing annual metric CSV ...")
  write.csv(annual_1km_df, file.path(OUT_ROOT, sprintf("annual_%s_metric_series_v9.csv", GRID_TAG)), row.names = FALSE)
  msg("Writing annual metric RDS without compression ...")
  saveRDS(annual_1km_df, file.path(OUT_ROOT, sprintf("annual_%s_metric_series_v9.rds", GRID_TAG)), compress = FALSE)
}
msg("Merging per-metric trend results ...")
trend_df <- dplyr::bind_rows(lapply(trend_objects, `[[`, "trend"))
msg("Trend rows merged: %s", nrow(trend_df))
msg("Applying FDR correction ...")
trend_df <- apply_fdr_by_metric(trend_df)
msg("FDR correction finished.")
msg("Writing trend CSV ...")
write.csv(trend_df, file.path(OUT_ROOT, sprintf("trend_%s_sen_HRmk_FDR_v9.csv", GRID_TAG)), row.names = FALSE)
msg("Trend CSV written.")
msg("Writing trend RDS without compression ...")
saveRDS(trend_df, file.path(OUT_ROOT, sprintf("trend_%s_sen_HRmk_FDR_v9.rds", GRID_TAG)), compress = FALSE)
msg("Trend RDS written.")

msg("Building wide trend table ...")
trend_wide <- trend_df %>%
  dplyr::mutate(pref = ifelse(metric == "Rich_tree", "rich", "shannon")) %>%
  dplyr::select(cell, trend_period, pref, n_years, sen_slope, mk_p_raw, p_hr, q_bh, q_by, q_bky, trend_state, trend_state_q010, acf1, mean_valid_frac) %>%
  tidyr::pivot_wider(names_from = pref, values_from = c(n_years, sen_slope, mk_p_raw, p_hr, q_bh, q_by, q_bky, trend_state, trend_state_q010, acf1, mean_valid_frac), names_glue = "{pref}_{.value}")

# Ensure names are pref-first.
names(trend_wide) <- gsub("^n_years_(rich|shannon)$", "\\1_n_years", names(trend_wide))
names(trend_wide) <- gsub("^sen_slope_(rich|shannon)$", "\\1_sen_slope", names(trend_wide))
names(trend_wide) <- gsub("^mk_p_raw_(rich|shannon)$", "\\1_mk_p_raw", names(trend_wide))
names(trend_wide) <- gsub("^p_hr_(rich|shannon)$", "\\1_p_hr", names(trend_wide))
names(trend_wide) <- gsub("^q_bh_(rich|shannon)$", "\\1_q_bh", names(trend_wide))
names(trend_wide) <- gsub("^q_by_(rich|shannon)$", "\\1_q_by", names(trend_wide))
names(trend_wide) <- gsub("^q_bky_(rich|shannon)$", "\\1_q_bky", names(trend_wide))
names(trend_wide) <- gsub("^trend_state_(rich|shannon)$", "\\1_trend_state", names(trend_wide))
names(trend_wide) <- gsub("^trend_state_q010_(rich|shannon)$", "\\1_trend_state_q010", names(trend_wide))
names(trend_wide) <- gsub("^acf1_(rich|shannon)$", "\\1_acf1", names(trend_wide))
names(trend_wide) <- gsub("^mean_valid_frac_(rich|shannon)$", "\\1_mean_valid_frac", names(trend_wide))
msg("Wide trend table built: %s rows x %s cols", nrow(trend_wide), ncol(trend_wide))

msg("Loading static age-100 q95 denominators on the common 1-km grid ...")
potential_by_period <- lapply(RUN_PERIODS, function(pn) {
  msg("  Static age-100 potential for state period %s", pn)
  msg("  Potential denominator period %s: Richness ...", pn)
  rich_pot <- potential_ceiling_1km_stable(METRICS[[1]], template = template)
  msg("  Potential denominator period %s: Shannon ...", pn)
  shannon_pot <- potential_ceiling_1km_stable(METRICS[[2]], template = template)
  msg("  Potential denominator period %s: finished.", pn)
  list(
    rich = rich_pot,
    shannon = shannon_pot
  )
})
names(potential_by_period) <- RUN_PERIODS

# Write RF potential denominator support diagnostics for reproducibility.
pot_diag <- dplyr::bind_rows(lapply(names(potential_by_period), function(pn) {
  pr <- potential_by_period[[pn]]
  dplyr::bind_rows(
    dplyr::mutate(pr$rich$diagnostics, period = pn, potential_reference_age = 100L),
    dplyr::mutate(pr$shannon$diagnostics, period = pn, potential_reference_age = 100L)
  )
}))
msg("Writing RF potential support diagnostics ...")
write.csv(pot_diag, file.path(OUT_ROOT, "RF_potential_support_diagnostics.csv"), row.names = FALSE)
msg("RF potential support diagnostics written.")
period_outputs <- list()
msg("Building RF six-class outputs ...")
for (tp in RUN_TREND_PERIODS) {
  tw <- trend_wide %>% dplyr::filter(.data$trend_period == tp)
  for (pn in RUN_PERIODS) {
    ag <- tp
    pot_pair <- potential_by_period[[pn]]
    period_outputs[[ag]] <- classify_period(pn, PERIODS[[pn]], tp, TREND_PERIODS[[tp]], ag, template, region_id, tw, pot_pair$rich, pot_pair$shannon)
  }
}

if (WRITE_COMBINED_CELL_TABLE) {
  msg("Combining period outputs ...")
  combined <- dplyr::bind_rows(period_outputs, .id = "analysis_group_key")
  msg("Combined period output rows: %s", nrow(combined))
  write.csv(combined, file.path(OUT_ROOT, "RF_6class_all_periods_cell_table.csv"), row.names = FALSE)
  write.csv(combined, file.path(OUT_ROOT, "RF_6class_all_trend_periods_cell_table.csv"), row.names = FALSE)
}

manifest <- data.frame(
  version = FEM_SOFTWARE_VERSION,
  data_dir = DATA_DIR_ABS,
  out_root = OUT_ROOT,
  years = paste(YEARS, collapse = ","),
  run_rf_periods = paste(RUN_PERIODS, collapse = ","),
  run_trend_periods = paste(RUN_TREND_PERIODS, collapse = ","),
  trend_periods = paste(vapply(RUN_TREND_PERIODS, function(tp) paste0(tp, ":", paste(TREND_PERIODS[[tp]]$years, collapse = "-")), character(1)), collapse = ";"),
  rf_thresh = RF_THRESH,
  q_main = Q_MAIN,
  q_wilks_sensitivity = Q_WILKS_SENS,
  min_valid_frac = MIN_VALID_FRAC,
  observed_years_min_frac = OBS_YEARS_MIN_FRAC,
  potential_reference_age = 100L,
  potential_quantile = 0.95,
  trend_min_valid_year_frac = TREND_MIN_VALID_YEAR_FRAC,
  rf_formula = "mean(Obs 2016-2020) / static q95(Forest_age=100)",
  potential_time_reference = "static (not an annual series)",
  trend_method = "Sen slope + Hamed-Rao modified MK",
  cell_size_m = CELL_SIZE_M,
  use_professional_trend_packages = USE_PROFESSIONAL_TREND_PACKAGES,
  use_trend_package_sen = USE_TREND_PACKAGE_SEN,
  use_modifiedmk_package = USE_MODIFIEDMK_PACKAGE,
  strict_professional_packages = STRICT_PROFESSIONAL_PACKAGES,
  fdr_main = "BH",
  fdr_bky_column = "q_bky (analytic two-stage TSBH by default; mutoss::two.stage only when USE_MUTOSS_BKY=TRUE and m<=BKY_MUTOSS_MAX_M)",
  write_annual_tables = WRITE_ANNUAL_TABLES,
  write_combined_cell_table = WRITE_COMBINED_CELL_TABLE,
  allow_ambiguous_input_files = ALLOW_AMBIGUOUS_INPUT_FILES,
  clean_existing = CLEAN_EXISTING,
  potential_support = "static_q95_age100_common_1km",
  stringsAsFactors = FALSE
)
write.csv(manifest, file.path(OUT_ROOT, "RF_6class_v9_run_manifest.csv"), row.names = FALSE)

msg("RF_6class FEM %s finished successfully.", FEM_SOFTWARE_VERSION)
