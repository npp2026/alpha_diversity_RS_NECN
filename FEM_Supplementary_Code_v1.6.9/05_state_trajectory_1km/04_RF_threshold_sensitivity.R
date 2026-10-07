#!/usr/bin/env Rscript
# ==============================================================================
# RF-threshold sensitivity for paired Shannon-minus-Richness class contrasts.
# Reuses spatial layouts/weights, defaults to support-aware block universes,
# and computes the baseline sign independently even when RF_SWEEP_R omits RF_THRESH.
# ==============================================================================
options(stringsAsFactors=FALSE)
if(getRversion()<'4.1.0')stop('R >= 4.1.0 is required.',call.=FALSE)
suppressPackageStartupMessages({library(dplyr);library(tidyr)})
get_script_dir<-function(){a<-commandArgs(FALSE);f<-a[grep('^--file=',a)];if(length(f))dirname(normalizePath(sub('^--file=','',f[1]),winslash='/',mustWork=FALSE))else getwd()}
SCRIPT_DIR<-get_script_dir();source(file.path(SCRIPT_DIR,'v12_posthoc_utils.R'))
DATA_DIR<-normalizePath(v12_env_str('DATA_DIR','.'),winslash='/',mustWork=TRUE)
OUT_ROOT<-v12_make_abs_path(v12_env_str('OUT_ROOT_NAME','outputs_RF_6class_fig6_minimal'),DATA_DIR)
RUN_GROUPS<-unique(as.character(v12_env_expr('RUN_OUTPUT_GROUPS_R',v12_env_expr('RUN_TREND_PERIODS_R',c('trend_2005_2020')))))
RF_SWEEP<-sort(unique(as.numeric(v12_env_expr('RF_SWEEP_R',c(0.70,0.75,0.80,0.85)))))
if(!length(RF_SWEEP)||any(!is.finite(RF_SWEEP)|RF_SWEEP<=0))stop('RF_SWEEP_R must contain positive finite thresholds.',call.=FALSE)
BASE_THRESHOLD<-v12_env_num('RF_THRESH',0.80,0,Inf,FALSE,TRUE);N_BOOT<-v12_env_int('N_BOOT',2000,1)
BLOCK_SIZES_KM<-sort(unique(as.numeric(v12_env_expr('BLOCK_SIZES_KM_R',c(50,75,100)))))
if(!length(BLOCK_SIZES_KM)||any(!is.finite(BLOCK_SIZES_KM)|BLOCK_SIZES_KM<=0))stop('BLOCK_SIZES_KM_R must contain positive values.',call.=FALSE)
CELL_SIZE_M<-v12_env_int('CELL_SIZE_M',1000,1);CELL_AREA_DEFAULT_KM2<-(CELL_SIZE_M/1000)^2
CELL_AREA_KM2<-v12_env_num('CELL_AREA_KM2',CELL_AREA_DEFAULT_KM2,0,Inf,FALSE,TRUE);BOOT_SEED<-v12_env_int('BOOT_SEED',20260609,1)
CI_LEVEL<-v12_env_num('CI_LEVEL',0.95,0,1,FALSE,FALSE);BLOCK_UNIVERSE<-v12_normalise_block_universe(v12_env_str('BOOT_BLOCK_UNIVERSE','analysis_support'))
BLOCK_ORIGIN_X<-v12_env_num('BLOCK_ORIGIN_X',0,-Inf,Inf,TRUE,TRUE);BLOCK_ORIGIN_Y<-v12_env_num('BLOCK_ORIGIN_Y',0,-Inf,Inf,TRUE,TRUE)
PRIMARY_SUPPORT<-v12_normalise_support_mode(v12_env_str('TARGET_CONTRAST_SUPPORT','common'));WRITE_BOTH_SUPPORTS<-v12_env_bool('WRITE_BOTH_SUPPORT_MODES',TRUE)

read_group<-function(group){f<-file.path(OUT_ROOT,group,'RF_6class_cell_table.csv');if(!file.exists(f))return(NULL);d<-v12_validate_cell_table(read.csv(f,check.names=FALSE),TRUE,f);v12_ensure_area(d,CELL_AREA_KM2)}
baseline_qc<-list();baseline_rows<-list();bq<-kb<-0L;rows<-list();k<-0L
for(gi in seq_along(RUN_GROUPS)){
 group<-RUN_GROUPS[gi];d0<-read_group(group);if(is.null(d0)){warning('Missing cell table for ',group,call.=FALSE);next}
 layouts<-setNames(lapply(BLOCK_SIZES_KM,function(bk)v12_prepare_block_layout(d0,bk,CELL_AREA_KM2,BLOCK_ORIGIN_X,BLOCK_ORIGIN_Y)),as.character(BLOCK_SIZES_KM))
 for(metric in c('rich','shannon')){
   pred<-v12_classify6(d0[[paste0(metric,'_RF')]],d0[[paste0(metric,'_trend_state')]],BASE_THRESHOLD);stored<-as.character(d0[[paste0(metric,'_class6')]])
   sv<-stored%in%V12_CLASS6_LEVELS;pv<-pred%in%V12_CLASS6_LEVELS;match_cls<-if(any(sv))mean(pred[sv]==stored[sv])else NA_real_;match_support<-mean(sv==pv)
   bq<-bq+1L;baseline_qc[[bq]]<-data.frame(trend_period=group,metric=metric,rf_threshold=BASE_THRESHOLD,n_rows=nrow(d0),classified_match_fraction=match_cls,support_match_fraction=match_support,stringsAsFactors=FALSE)
   if((is.finite(match_cls)&&match_cls<0.999999)||match_support<0.999999)stop('RF-threshold baseline reconstruction failed for ',group,'/',metric,'. Check RF_THRESH.',call.=FALSE)
 }
 supports<-unique(c(PRIMARY_SUPPORT,if(WRITE_BOTH_SUPPORTS)c('common','target_specific')else character()))
 # Independent baseline contrasts make custom sweeps well-defined even without RF_THRESH.
 for(support in supports)for(bk in BLOCK_SIZES_KM){
   x<-v12_metric_contrast(d0,bk,N_BOOT,v12_seed_offset(BOOT_SEED,gi,bk),support,CELL_AREA_KM2,
      layout=layouts[[as.character(bk)]],block_universe=BLOCK_UNIVERSE,ci_level=CI_LEVEL,origin_x=BLOCK_ORIGIN_X,origin_y=BLOCK_ORIGIN_Y)
   x$trend_period<-group;x$rf_threshold<-BASE_THRESHOLD;x$block_size_km<-bk
   kb<-kb+1L;baseline_rows[[kb]]<-x
 }
 for(thr in RF_SWEEP){
   d<-d0;d$rich_class6<-v12_classify6(d$rich_RF,d$rich_trend_state,thr);d$shannon_class6<-v12_classify6(d$shannon_RF,d$shannon_trend_state,thr)
   for(support in supports)for(bk in BLOCK_SIZES_KM){
     r<-v12_metric_contrast(d,bk,N_BOOT,v12_seed_offset(BOOT_SEED,gi,bk),support,CELL_AREA_KM2,
       layout=layouts[[as.character(bk)]],block_universe=BLOCK_UNIVERSE,ci_level=CI_LEVEL,origin_x=BLOCK_ORIGIN_X,origin_y=BLOCK_ORIGIN_Y)
     r$trend_period<-group;r$rf_threshold<-thr;r$block_size_km<-bk;r$n_cells<-nrow(d0);r$contrast_definition<-'Shannon_minus_Richness'
     k<-k+1L;rows[[k]]<-r
   }
 }
}
if(!length(rows))stop('No usable cell tables for RF-threshold sensitivity.',call.=FALSE)
grid_all<-dplyr::bind_rows(rows)|>dplyr::mutate(class6=factor(.data$class6,levels=V12_CLASS6_LEVELS),zone=factor(.data$zone,levels=V12_REG_ORDER))|>
 dplyr::arrange(.data$trend_period,.data$support_mode,.data$rf_threshold,.data$block_size_km,.data$zone,.data$class6)
base_all<-dplyr::bind_rows(baseline_rows)
primary<-grid_all|>dplyr::filter(as.character(.data$support_mode)==PRIMARY_SUPPORT)
out_dir<-file.path(OUT_ROOT,'metric_diff_CI_v9');dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)
write.csv(primary,file.path(out_dir,'RF_6class_metric_diff_RFthreshold_grid_v9.csv'),row.names=FALSE)
write.csv(grid_all,file.path(out_dir,'RF_6class_metric_diff_RFthreshold_grid_all_support_modes_v13.csv'),row.names=FALSE)
write.csv(dplyr::bind_rows(baseline_qc),file.path(out_dir,'RF_6class_RFthreshold_baseline_reproduction_QC_v13.csv'),row.names=FALSE)
write.csv(base_all,file.path(out_dir,'RF_6class_RFthreshold_baseline_contrasts_v13.csv'),row.names=FALSE)
safe_min<-function(x){x<-x[is.finite(x)];if(length(x))min(x)else NA_real_};safe_max<-function(x){x<-x[is.finite(x)];if(length(x))max(x)else NA_real_}
rob<-grid_all|>dplyr::group_by(.data$trend_period,.data$zone,.data$class6,.data$support_mode)|>dplyr::group_modify(function(.x,.y){
 n<-nrow(.x);sig<-.x$ci_excludes_zero%in%TRUE;pos<-sig&is.finite(.x$delta_low)&.x$delta_low>0;neg<-sig&is.finite(.x$delta_high)&.x$delta_high<0
 sign_when_sig<-if(!sum(sig))'ns'else if(any(pos)&&!any(neg))'+'else if(any(neg)&&!any(pos))'-'else'mixed'
 b<-base_all|>dplyr::filter(.data$trend_period==.y$trend_period,.data$zone==.y$zone,.data$class6==.y$class6,.data$support_mode==.y$support_mode)
 bv<-b$delta_metric[is.finite(b$delta_metric)];baseline_sign<-if(!length(bv))'0/NA'else if(mean(bv)>0)'+'else if(mean(bv)<0)'-'else'0/NA'
 finite<-is.finite(.x$delta_metric);same<-if(baseline_sign=='+').x$delta_metric>0 else if(baseline_sign=='-').x$delta_metric<0 else rep(FALSE,n)
 data.frame(n_thresholds=dplyr::n_distinct(.x$rf_threshold),n_block_sizes=dplyr::n_distinct(.x$block_size_km),n_combos=n,n_sig=sum(sig),n_sig_pos=sum(pos),n_sig_neg=sum(neg),
 sig_fraction=if(n>0)sum(sig)/n else NA_real_,baseline_sign=baseline_sign,sign_when_sig=sign_when_sig,sign_stable_fraction=if(any(finite))mean(same[finite])else NA_real_,
 sign_stable_all=if(any(finite))all(same[finite])else NA,direction_robust=n>0&&sum(sig)==n&&sign_when_sig%in%c('+','-')&&sign_when_sig==baseline_sign,
 mostly_direction_robust=n>0&&sum(sig)/n>=0.8&&sign_when_sig%in%c('+','-')&&sign_when_sig==baseline_sign,delta_min=safe_min(.x$delta_metric),delta_max=safe_max(.x$delta_metric),stringsAsFactors=FALSE)
})|>dplyr::ungroup()|>dplyr::arrange(.data$trend_period,.data$support_mode,.data$zone,.data$class6)
write.csv(rob|>dplyr::filter(as.character(.data$support_mode)==PRIMARY_SUPPORT),file.path(out_dir,'RF_6class_metric_diff_RFthreshold_direction_robustness_v9.csv'),row.names=FALSE)
write.csv(rob,file.path(out_dir,'RF_6class_metric_diff_RFthreshold_direction_robustness_all_support_modes_v13.csv'),row.names=FALSE)
v12_msg('[rf_sens] support=%s; universe=%s; rows=%d; thresholds=%s',PRIMARY_SUPPORT,BLOCK_UNIVERSE,nrow(primary),paste(RF_SWEEP,collapse=','))
