#!/usr/bin/env Rscript
# ==============================================================================
# Fast, support-aware spatial block-bootstrap CIs for six-class area composition.
# Reuses block layouts, samples only estimand-support blocks by default, and
# writes explicit bootstrap/block diagnostics. Set BOOT_BLOCK_UNIVERSE=all_rows
# only to reproduce the legacy global-block sampling frame.
# ==============================================================================
options(stringsAsFactors = FALSE)
if (getRversion() < '4.1.0') stop('R >= 4.1.0 is required.', call.=FALSE)
suppressPackageStartupMessages({library(dplyr); library(tidyr)})
get_script_dir <- function(){a<-commandArgs(FALSE);f<-a[grep('^--file=',a)];if(length(f))dirname(normalizePath(sub('^--file=','',f[1]),winslash='/',mustWork=FALSE))else getwd()}
SCRIPT_DIR <- get_script_dir(); source(file.path(SCRIPT_DIR,'v12_posthoc_utils.R'))

DATA_DIR <- normalizePath(v12_env_str('DATA_DIR','.'),winslash='/',mustWork=TRUE)
OUT_ROOT <- v12_make_abs_path(v12_env_str('OUT_ROOT_NAME','outputs_RF_6class_fig6_minimal'),DATA_DIR)
RUN_GROUPS <- unique(as.character(v12_env_expr('RUN_OUTPUT_GROUPS_R',v12_env_expr('RUN_TREND_PERIODS_R',c('trend_2005_2020')))))
GROUP_SINGLE <- v12_env_str('ANALYSIS_GROUP',''); if(nzchar(GROUP_SINGLE)) RUN_GROUPS <- GROUP_SINGLE
N_BOOT <- v12_env_int('N_BOOT',2000,1)
BLOCK_SIZES_KM <- sort(unique(as.numeric(v12_env_expr('BLOCK_SIZES_KM_R',c(50,75,100)))))
if(!length(BLOCK_SIZES_KM)||any(!is.finite(BLOCK_SIZES_KM)|BLOCK_SIZES_KM<=0))stop('BLOCK_SIZES_KM_R must contain positive values.',call.=FALSE)
CELL_SIZE_M <- v12_env_int('CELL_SIZE_M',1000,1); CELL_AREA_DEFAULT_KM2 <- (CELL_SIZE_M/1000)^2
CELL_AREA_KM2 <- v12_env_num('CELL_AREA_KM2',CELL_AREA_DEFAULT_KM2,0,Inf,FALSE,TRUE)
BOOT_SEED <- v12_env_int('BOOT_SEED',20260609,1)
CI_LEVEL <- v12_env_num('CI_LEVEL',0.95,0,1,FALSE,FALSE)
BLOCK_UNIVERSE <- v12_normalise_block_universe(v12_env_str('BOOT_BLOCK_UNIVERSE','analysis_support'))
BLOCK_ORIGIN_X <- v12_env_num('BLOCK_ORIGIN_X',0,-Inf,Inf,TRUE,TRUE)
BLOCK_ORIGIN_Y <- v12_env_num('BLOCK_ORIGIN_Y',0,-Inf,Inf,TRUE,TRUE)
MIN_BOOT_VALID_FRACTION <- v12_env_num('MIN_BOOT_VALID_FRACTION',0.95,0,1,TRUE,TRUE)
MIN_ZONE_BLOCKS_WARN <- v12_env_int('MIN_ZONE_BLOCKS_WARN',5,1)

process_group <- function(group, gi) {
  f <- file.path(OUT_ROOT,group,'RF_6class_cell_table.csv')
  if(!file.exists(f)){warning('Missing cell table: ',f,call.=FALSE);return(NULL)}
  d <- v12_validate_cell_table(read.csv(f,check.names=FALSE),FALSE,f)
  d <- v12_ensure_area(d,CELL_AREA_KM2)
  period <- if('period'%in%names(d))as.character(d$period[1])else NA_character_
  trend_period <- if('trend_period'%in%names(d))as.character(d$trend_period[1])else group
  layouts <- setNames(lapply(BLOCK_SIZES_KM,function(bk)v12_prepare_block_layout(d,bk,CELL_AREA_KM2,BLOCK_ORIGIN_X,BLOCK_ORIGIN_Y)),as.character(BLOCK_SIZES_KM))
  rows <- list(); k <- 0L
  for(bk in BLOCK_SIZES_KM) for(metric in c('Richness','Shannon')) {
    col <- if(metric=='Richness')'rich_class6'else'shannon_class6'
    seed <- v12_seed_offset(BOOT_SEED,gi,bk)
    r <- v12_area_ci_metric(d,col,metric,bk,N_BOOT,seed,CELL_AREA_KM2,
      layout=layouts[[as.character(bk)]],block_universe=BLOCK_UNIVERSE,ci_level=CI_LEVEL,
      origin_x=BLOCK_ORIGIN_X,origin_y=BLOCK_ORIGIN_Y)
    k <- k+1L; rows[[k]] <- r
  }
  ci <- dplyr::bind_rows(rows) |>
    dplyr::mutate(analysis_group=group,period=period,trend_period=trend_period,
      class6=factor(.data$class6,levels=V12_CLASS6_LEVELS),zone=factor(.data$zone,levels=V12_REG_ORDER),
      metric=factor(.data$metric,levels=V12_METRIC_LEVELS)) |>
    dplyr::arrange(.data$block_size_km,.data$metric,.data$zone,.data$class6) |>
    dplyr::select(.data$analysis_group,.data$period,.data$trend_period,.data$block_size_km,
      .data$metric,.data$zone,.data$class6,.data$prop_median,.data$prop_low,.data$prop_high,
      .data$prop,.data$n_cells,.data$area_km2,.data$n_boot,.data$bootstrap_valid_fraction,.data$ci_level,
      .data$n_blocks_universe,.data$n_blocks_zone,.data$block_universe_mode,.data$block_origin_x,.data$block_origin_y)
  weak <- ci |> dplyr::filter(as.character(.data$zone)!='outside',.data$n_cells>0,.data$n_blocks_zone<MIN_ZONE_BLOCKS_WARN)
  if(nrow(weak)) warning(group,': ',nrow(weak),' class/zone rows have < ',MIN_ZONE_BLOCKS_WARN,' active spatial blocks; interpret CIs cautiously.',call.=FALSE)
  bad_boot <- ci |> dplyr::filter(.data$n_cells>0,.data$bootstrap_valid_fraction<MIN_BOOT_VALID_FRACTION)
  if(nrow(bad_boot)) warning(group,': ',nrow(bad_boot),' rows have bootstrap_valid_fraction < ',MIN_BOOT_VALID_FRACTION,'.',call.=FALSE)
  point <- ci |> dplyr::filter(.data$block_size_km==BLOCK_SIZES_KM[1]) |>
    dplyr::select(.data$analysis_group,.data$period,.data$trend_period,.data$metric,.data$zone,.data$class6,.data$n_cells,.data$area_km2,.data$prop)
  write.csv(point,file.path(OUT_ROOT,group,'RF_6class_area_point_estimates_v9.csv'),row.names=FALSE)
  write.csv(ci,file.path(OUT_ROOT,group,'RF_6class_area_block_bootstrap_CI_v9.csv'),row.names=FALSE)
  v12_msg('[area_ci] %s: rows=%d, boot=%d, universe=%s, blocks=%s',group,nrow(ci),N_BOOT,BLOCK_UNIVERSE,paste(BLOCK_SIZES_KM,collapse=','))
  list(point=point,ci=ci)
}

res <- lapply(seq_along(RUN_GROUPS),function(i)process_group(RUN_GROUPS[i],i)); res <- res[!vapply(res,is.null,logical(1))]
if(length(res)) {
  all_point <- dplyr::bind_rows(lapply(res,`[[`,'point')); all_ci <- dplyr::bind_rows(lapply(res,`[[`,'ci'))
  for(nm in c('all_trend_periods_CI_v9','all_periods_CI_v9')) {
    od <- file.path(OUT_ROOT,nm); dir.create(od,recursive=TRUE,showWarnings=FALSE)
    write.csv(all_point,file.path(od,if(nm=='all_trend_periods_CI_v9')'RF_6class_area_point_estimates_all_trend_periods_v9.csv'else'RF_6class_area_point_estimates_all_periods_v9.csv'),row.names=FALSE)
    write.csv(all_ci,file.path(od,if(nm=='all_trend_periods_CI_v9')'RF_6class_area_block_bootstrap_CI_all_trend_periods_v9.csv'else'RF_6class_area_block_bootstrap_CI_all_periods_v9.csv'),row.names=FALSE)
  }
}
