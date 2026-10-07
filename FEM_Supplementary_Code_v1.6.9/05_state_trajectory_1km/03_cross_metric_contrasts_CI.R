#!/usr/bin/env Rscript
# ==============================================================================
# Paired Shannon-minus-Richness six-class contrasts with support-aware block CI.
# Primary estimand is common classified support; target-specific denominators are
# retained as sensitivity. Delta is always p_Shannon - p_Richness.
# ==============================================================================
options(stringsAsFactors=FALSE)
if(getRversion()<'4.1.0')stop('R >= 4.1.0 is required.',call.=FALSE)
suppressPackageStartupMessages({library(dplyr);library(tidyr)})
get_script_dir<-function(){a<-commandArgs(FALSE);f<-a[grep('^--file=',a)];if(length(f))dirname(normalizePath(sub('^--file=','',f[1]),winslash='/',mustWork=FALSE))else getwd()}
SCRIPT_DIR<-get_script_dir();source(file.path(SCRIPT_DIR,'v12_posthoc_utils.R'))

DATA_DIR<-normalizePath(v12_env_str('DATA_DIR','.'),winslash='/',mustWork=TRUE)
OUT_ROOT<-v12_make_abs_path(v12_env_str('OUT_ROOT_NAME','outputs_RF_6class_fig6_minimal'),DATA_DIR)
RUN_GROUPS<-unique(as.character(v12_env_expr('RUN_OUTPUT_GROUPS_R',v12_env_expr('RUN_TREND_PERIODS_R',c('trend_2005_2020')))))
GROUP_SINGLE<-v12_env_str('ANALYSIS_GROUP','');if(nzchar(GROUP_SINGLE))RUN_GROUPS<-GROUP_SINGLE
N_BOOT<-v12_env_int('N_BOOT',2000,1)
BLOCK_SIZES_KM<-sort(unique(as.numeric(v12_env_expr('BLOCK_SIZES_KM_R',c(50,75,100)))))
if(!length(BLOCK_SIZES_KM)||any(!is.finite(BLOCK_SIZES_KM)|BLOCK_SIZES_KM<=0))stop('BLOCK_SIZES_KM_R must contain positive values.',call.=FALSE)
CELL_SIZE_M<-v12_env_int('CELL_SIZE_M',1000,1);CELL_AREA_DEFAULT_KM2<-(CELL_SIZE_M/1000)^2
CELL_AREA_KM2<-v12_env_num('CELL_AREA_KM2',CELL_AREA_DEFAULT_KM2,0,Inf,FALSE,TRUE);BOOT_SEED<-v12_env_int('BOOT_SEED',20260609,1)
CI_LEVEL<-v12_env_num('CI_LEVEL',0.95,0,1,FALSE,FALSE)
BLOCK_UNIVERSE<-v12_normalise_block_universe(v12_env_str('BOOT_BLOCK_UNIVERSE','analysis_support'))
BLOCK_ORIGIN_X<-v12_env_num('BLOCK_ORIGIN_X',0,-Inf,Inf,TRUE,TRUE);BLOCK_ORIGIN_Y<-v12_env_num('BLOCK_ORIGIN_Y',0,-Inf,Inf,TRUE,TRUE)
PRIMARY_SUPPORT<-v12_normalise_support_mode(v12_env_str('TARGET_CONTRAST_SUPPORT','common'))
WRITE_BOTH_SUPPORTS<-v12_env_bool('WRITE_BOTH_SUPPORT_MODES',TRUE);CLEAN_POSTHOC_OUTPUT<-v12_env_bool('CLEAN_POSTHOC_OUTPUT',FALSE)
MIN_BOOT_VALID_FRACTION<-v12_env_num('MIN_BOOT_VALID_FRACTION',0.95,0,1,TRUE,TRUE)

read_group<-function(group){f<-file.path(OUT_ROOT,group,'RF_6class_cell_table.csv');if(!file.exists(f))return(NULL);d<-v12_validate_cell_table(read.csv(f,check.names=FALSE),FALSE,f);v12_ensure_area(d,CELL_AREA_KM2)}
rows<-list();k<-0L
for(gi in seq_along(RUN_GROUPS)){
 group<-RUN_GROUPS[gi];d<-read_group(group);if(is.null(d)){warning('Missing cell table for ',group,call.=FALSE);next}
 layouts<-setNames(lapply(BLOCK_SIZES_KM,function(bk)v12_prepare_block_layout(d,bk,CELL_AREA_KM2,BLOCK_ORIGIN_X,BLOCK_ORIGIN_Y)),as.character(BLOCK_SIZES_KM))
 supports<-unique(c(PRIMARY_SUPPORT,if(WRITE_BOTH_SUPPORTS)c('common','target_specific')else character()))
 for(support in supports)for(bk in BLOCK_SIZES_KM){
   seed<-v12_seed_offset(BOOT_SEED,gi,bk)
   r<-v12_metric_contrast(d,bk,N_BOOT,seed,support,CELL_AREA_KM2,layout=layouts[[as.character(bk)]],
     block_universe=BLOCK_UNIVERSE,ci_level=CI_LEVEL,origin_x=BLOCK_ORIGIN_X,origin_y=BLOCK_ORIGIN_Y)
   r$trend_period<-group;r$block_size_km<-bk;r$n_cells<-nrow(d);r$contrast_definition<-'Shannon_minus_Richness'
   k<-k+1L;rows[[k]]<-r
 }
}
if(!length(rows))stop('No usable cell tables for metric-difference analysis.',call.=FALSE)
all_support<-dplyr::bind_rows(rows)|>
 dplyr::mutate(class6=factor(.data$class6,levels=V12_CLASS6_LEVELS),zone=factor(.data$zone,levels=V12_REG_ORDER))|>
 dplyr::arrange(.data$trend_period,.data$support_mode,.data$block_size_km,.data$zone,.data$class6)|>
 dplyr::select(.data$trend_period,.data$zone,.data$class6,.data$block_size_km,.data$support_mode,.data$contrast_definition,
   .data$prop_shannon,.data$prop_richness,.data$delta_metric,.data$delta_median,.data$delta_low,.data$delta_high,
   .data$ci_excludes_zero,.data$n_boot,.data$bootstrap_valid_fraction,.data$ci_level,.data$n_cells,
   .data$n_blocks_universe,.data$n_blocks_rich_zone,.data$n_blocks_shannon_zone,.data$block_universe_mode,.data$block_origin_x,.data$block_origin_y,
   .data$n_union_cells,.data$n_common_cells,.data$common_area_km2,.data$rich_classified_area_km2,.data$shannon_classified_area_km2)
if(any(all_support$n_common_cells>0 & all_support$bootstrap_valid_fraction<MIN_BOOT_VALID_FRACTION,na.rm=TRUE))warning('Some metric contrasts have low bootstrap-valid fractions; inspect output diagnostics.',call.=FALSE)
out_dir<-file.path(OUT_ROOT,'metric_diff_CI_v9');if(CLEAN_POSTHOC_OUTPUT&&dir.exists(out_dir))unlink(out_dir,recursive=TRUE,force=TRUE);dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)
write.csv(all_support,file.path(out_dir,'RF_6class_metric_diff_within_period_all_support_modes_v13.csv'),row.names=FALSE)
primary<-all_support|>dplyr::filter(as.character(.data$support_mode)==PRIMARY_SUPPORT)
write.csv(primary,file.path(out_dir,'RF_6class_metric_diff_within_period_Shannon_minus_Richness_v9.csv'),row.names=FALSE)
write.csv(all_support|>dplyr::filter(as.character(.data$support_mode)=='common'),file.path(out_dir,'RF_6class_metric_diff_within_period_common_support_v13.csv'),row.names=FALSE)
write.csv(all_support|>dplyr::filter(as.character(.data$support_mode)=='target_specific'),file.path(out_dir,'RF_6class_metric_diff_within_period_target_specific_v13.csv'),row.names=FALSE)
qc<-primary|>dplyr::filter(.data$block_size_km==BLOCK_SIZES_KM[1],as.character(.data$class6)==V12_CLASS6_LEVELS[1])|>
 dplyr::transmute(trend_period=.data$trend_period,zone=.data$zone,support_mode=.data$support_mode,n_union_cells=.data$n_union_cells,n_common_cells=.data$n_common_cells,
 common_support_fraction=ifelse(.data$n_union_cells>0,.data$n_common_cells/.data$n_union_cells,NA_real_),common_area_km2=.data$common_area_km2,
 rich_classified_area_km2=.data$rich_classified_area_km2,shannon_classified_area_km2=.data$shannon_classified_area_km2,n_blocks_universe=.data$n_blocks_universe)
write.csv(qc,file.path(out_dir,'RF_6class_target_support_QC_v13.csv'),row.names=FALSE)
v12_msg('[metric_diff] support=%s; universe=%s; primary rows=%d; all rows=%d',PRIMARY_SUPPORT,BLOCK_UNIVERSE,nrow(primary),nrow(all_support))
