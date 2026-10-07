#!/usr/bin/env Rscript
# ==============================================================================
# Ratio-component robustness v13 (R implementation for the minimal package).
# RF* = (m_D / m_B) * RF. Expensive classification/bootstrap is performed once
# per distinct ratio, then expanded back to the requested (m_D,m_B) grid.
# Primary target contrast uses common classified support by default.
# ==============================================================================
options(stringsAsFactors=FALSE)
if(getRversion()<'4.1.0')stop('R >= 4.1.0 is required.',call.=FALSE)
suppressPackageStartupMessages({library(dplyr);library(tidyr)})
get_script_dir<-function(){a<-commandArgs(FALSE);f<-a[grep('^--file=',a)];if(length(f))dirname(normalizePath(sub('^--file=','',f[1]),winslash='/',mustWork=FALSE))else getwd()}
SCRIPT_DIR<-get_script_dir();source(file.path(SCRIPT_DIR,'v12_posthoc_utils.R'))

DATA_DIR<-normalizePath(v12_env_str('DATA_DIR','.'),winslash='/',mustWork=TRUE)
OUT_ROOT<-v12_make_abs_path(v12_env_str('OUT_ROOT_NAME','outputs_RF_6class_fig6_minimal'),DATA_DIR)
RUN_GROUPS<-as.character(v12_env_expr('RUN_OUTPUT_GROUPS_R',v12_env_expr('RUN_TREND_PERIODS_R',c('trend_2005_2020'))))
MD_GRID<-as.numeric(v12_env_expr('MD_GRID_R',c(0.8,0.9,1.0,1.1,1.2)))
MB_GRID<-as.numeric(v12_env_expr('MB_GRID_R',c(0.8,0.9,1.0,1.1,1.2)))
if(any(!is.finite(c(MD_GRID,MB_GRID))|c(MD_GRID,MB_GRID)<=0))stop('MD_GRID_R and MB_GRID_R must contain positive finite values.',call.=FALSE)
RF_THRESH<-v12_env_num('RF_THRESH',0.80,0,Inf,FALSE,TRUE)
N_BOOT<-v12_env_int('N_BOOT',2000,1)
BLOCK_SIZES_KM<-sort(unique(as.numeric(v12_env_expr('BLOCK_SIZES_KM_R',c(50,75,100)))))
if(!length(BLOCK_SIZES_KM)||any(!is.finite(BLOCK_SIZES_KM)|BLOCK_SIZES_KM<=0))stop('BLOCK_SIZES_KM_R must contain positive values.',call.=FALSE)
CELL_SIZE_M<-v12_env_int('CELL_SIZE_M',1000,1);CELL_AREA_DEFAULT_KM2<-(CELL_SIZE_M/1000)^2
CELL_AREA_KM2<-v12_env_num('CELL_AREA_KM2',CELL_AREA_DEFAULT_KM2,0,Inf,FALSE,TRUE)
CI_LEVEL<-v12_env_num('CI_LEVEL',0.95,0,1,FALSE,FALSE)
BLOCK_UNIVERSE<-v12_normalise_block_universe(v12_env_str('BOOT_BLOCK_UNIVERSE','analysis_support'))
BLOCK_ORIGIN_X<-v12_env_num('BLOCK_ORIGIN_X',0,-Inf,Inf,TRUE,TRUE)
BLOCK_ORIGIN_Y<-v12_env_num('BLOCK_ORIGIN_Y',0,-Inf,Inf,TRUE,TRUE)
BOOT_SEED<-v12_env_int('BOOT_SEED',20260609,1)
SUPPORT<-v12_normalise_support_mode(v12_env_str('TARGET_CONTRAST_SUPPORT','common'))
WRITE_BOTH_SUPPORTS<-v12_env_bool('WRITE_BOTH_SUPPORT_MODES',TRUE)
DO_SHARES<-v12_env_bool('RATIO_SHARE_CI',TRUE)
MOSTLY_MIN<-v12_env_num('RATIO_MOSTLY_ROBUST_MIN',0.80,0,1,TRUE,TRUE)
CLEAN_RATIO_OUTPUT<-v12_env_bool('CLEAN_RATIO_OUTPUT',TRUE)
FOCUS<-as.character(v12_env_expr('RATIO_FOCUS_CLASSES',V12_FOCAL4_LEVELS));FOCUS<-intersect(FOCUS,V12_CLASS6_LEVELS)
if(!length(FOCUS))stop('RATIO_FOCUS_CLASSES resolved to no canonical classes.',call.=FALSE)

out_dir<-file.path(OUT_ROOT,'ratio_robustness_v13');if(CLEAN_RATIO_OUTPUT&&dir.exists(out_dir))unlink(out_dir,recursive=TRUE,force=TRUE);dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)
mult_grid<-expand.grid(m_D=MD_GRID,m_B=MB_GRID,KEEP.OUT.ATTRS=FALSE,stringsAsFactors=FALSE)
mult_grid$ratio<-mult_grid$m_D/mult_grid$m_B
mult_grid$ratio_key<-sprintf('%.17g',mult_grid$ratio)
distinct<-mult_grid|>dplyr::distinct(.data$ratio_key,.keep_all=TRUE)|>dplyr::arrange(.data$ratio)
# Always evaluate ratio=1 internally so baseline direction is defined even for custom grids.
if(!any(abs(distinct$ratio-1)<1e-12))distinct<-dplyr::bind_rows(distinct,data.frame(m_D=1,m_B=1,ratio=1,ratio_key=sprintf('%.17g',1)))|>dplyr::arrange(.data$ratio)
v12_msg('[ratio_robust] %d (mD,mB) cells -> %d distinct ratios; T=%.3f',nrow(mult_grid),nrow(distinct),RF_THRESH)

read_group<-function(group){
 f<-file.path(OUT_ROOT,group,'RF_6class_cell_table.csv');if(!file.exists(f))return(NULL)
 d<-read.csv(f,check.names=FALSE);d<-v12_validate_cell_table(d,TRUE,f);v12_ensure_area(d,CELL_AREA_KM2)
}

agreement_rows<-function(star,base,group,ratio){
 out<-list();k<-0L;zones<-v12_zone_value(star$roman)
 for(metric in c('Richness','Shannon')){
  col<-if(metric=='Richness')'rich_class6'else'shannon_class6'
  a<-as.character(star[[col]]);b<-as.character(base[[col]]);valid<-a%in%V12_CLASS6_LEVELS & b%in%V12_CLASS6_LEVELS
  for(z in V12_REG_ORDER){
   zm<-if(z=='Overall')rep(TRUE,nrow(star))else zones==z;ok<-valid&zm;w<-star$area_km2[ok]
   wm<-function(x)if(length(w)&&sum(w)>0)sum(w*as.numeric(x))/sum(w)else NA_real_
   star_valid<-a%in%V12_CLASS6_LEVELS;base_valid<-b%in%V12_CLASS6_LEVELS
   den_base<-sum(star$area_km2[zm & base_valid],na.rm=TRUE);den_star<-sum(star$area_km2[zm & star_valid],na.rm=TRUE)
   k<-k+1L;out[[k]]<-data.frame(trend_period=group,ratio=ratio,metric=metric,zone=z,n_common_classified_cells=sum(ok),
    common_area_km2=sum(w),baseline_classified_area_km2=den_base,perturbed_classified_area_km2=den_star,
    support_retention_fraction=if(den_base>0)den_star/den_base else NA_real_,
    hl_agreement=wm(substr(a[ok],1,1)==substr(b[ok],1,1)),class_agreement=wm(a[ok]==b[ok]),stringsAsFactors=FALSE)
  }
 }
 dplyr::bind_rows(out)
}

contrast_distinct<-list();shares_distinct<-list();agree_all<-list();qc_all<-list();kc<-ks<-ka<-kq<-0L
for(gi in seq_along(RUN_GROUPS)){
 group<-RUN_GROUPS[gi];d0<-read_group(group);if(is.null(d0)){warning('Missing cell table for ',group,call.=FALSE);next}
 layouts<-setNames(lapply(BLOCK_SIZES_KM,function(bk)v12_prepare_block_layout(d0,bk,CELL_AREA_KM2,BLOCK_ORIGIN_X,BLOCK_ORIGIN_Y)),as.character(BLOCK_SIZES_KM))
 for(metric in c('rich','shannon')){
  pred<-v12_classify6(d0[[paste0(metric,'_RF')]],d0[[paste0(metric,'_trend_state')]],RF_THRESH);stored<-as.character(d0[[paste0(metric,'_class6')]])
  sv<-stored%in%V12_CLASS6_LEVELS;pv<-pred%in%V12_CLASS6_LEVELS;mf<-if(any(sv))mean(pred[sv]==stored[sv])else NA_real_;sf<-mean(sv==pv)
  kq<-kq+1L;qc_all[[kq]]<-data.frame(trend_period=group,metric=metric,rf_threshold=RF_THRESH,n_rows=nrow(d0),classified_match_fraction=mf,support_match_fraction=sf,stringsAsFactors=FALSE)
  if((is.finite(mf)&&mf<0.999999)||sf<0.999999)stop('Baseline ratio reconstruction failed for ',group,'/',metric,'.',call.=FALSE)
 }
 supports<-if(WRITE_BOTH_SUPPORTS)c('common','target_specific')else SUPPORT;supports<-unique(c(SUPPORT,supports))
 for(ri in seq_len(nrow(distinct))){
  ratio<-distinct$ratio[ri];rkey<-distinct$ratio_key[ri];d<-d0
  d$rich_class6<-v12_classify6(ratio*d$rich_RF,d$rich_trend_state,RF_THRESH)
  d$shannon_class6<-v12_classify6(ratio*d$shannon_RF,d$shannon_trend_state,RF_THRESH)
  ka<-ka+1L;agree_all[[ka]]<-agreement_rows(d,d0,group,ratio)|>dplyr::mutate(ratio_key=rkey)
  for(support in supports)for(bk in BLOCK_SIZES_KM){
   seed<-v12_seed_offset(BOOT_SEED,gi,bk)
   x<-v12_metric_contrast(d,bk,N_BOOT,seed,support,CELL_AREA_KM2,layout=layouts[[as.character(bk)]],block_universe=BLOCK_UNIVERSE,ci_level=CI_LEVEL,origin_x=BLOCK_ORIGIN_X,origin_y=BLOCK_ORIGIN_Y)
   x$trend_period<-group;x$ratio<-ratio;x$ratio_key<-rkey;x$block_size_km<-bk;x$rf_threshold<-RF_THRESH;x$contrast_definition<-'Shannon_minus_Richness'
   kc<-kc+1L;contrast_distinct[[kc]]<-x
  }
  if(DO_SHARES)for(bk in BLOCK_SIZES_KM)for(metric in c('Richness','Shannon')){
   col<-if(metric=='Richness')'rich_class6'else'shannon_class6';seed<-v12_seed_offset(BOOT_SEED,gi,bk,if(metric=='Richness')7L else 13L)
   x<-v12_area_ci_metric(d,col,metric,bk,N_BOOT,seed,CELL_AREA_KM2,layout=layouts[[as.character(bk)]],block_universe=BLOCK_UNIVERSE,ci_level=CI_LEVEL,origin_x=BLOCK_ORIGIN_X,origin_y=BLOCK_ORIGIN_Y)|>dplyr::filter(as.character(.data$class6)%in%FOCUS)
   x$trend_period<-group;x$ratio<-ratio;x$ratio_key<-rkey;x$rf_threshold<-RF_THRESH
   ks<-ks+1L;shares_distinct[[ks]]<-x
  }
 }
}
if(!length(contrast_distinct))stop('No ratio robustness results were produced.',call.=FALSE)
dist_con<-dplyr::bind_rows(contrast_distinct)
# Expand equivalent ratios back to requested multiplier grid.
grid<-dplyr::inner_join(mult_grid,dist_con,by='ratio_key',suffix=c('_requested',''))|>
 dplyr::mutate(ratio=.data$ratio_requested)|>dplyr::select(-dplyr::any_of('ratio_requested'))|>
 dplyr::arrange(.data$trend_period,.data$m_D,.data$m_B,.data$support_mode,.data$block_size_km,.data$zone,.data$class6)
write.csv(grid,file.path(out_dir,'RF_6class_ratio_perturb_grid_v13.csv'),row.names=FALSE)
write.csv(dist_con,file.path(out_dir,'RF_6class_ratio_perturb_distinct_ratios_v13.csv'),row.names=FALSE)
write.csv(dplyr::bind_rows(qc_all),file.path(out_dir,'RF_6class_ratio_baseline_reproduction_QC_v13.csv'),row.names=FALSE)
write.csv(dplyr::bind_rows(agree_all),file.path(out_dir,'RF_6class_ratio_perturb_HL_agreement_v13.csv'),row.names=FALSE)
if(length(shares_distinct)){
 sh<-dplyr::inner_join(mult_grid,dplyr::bind_rows(shares_distinct),by='ratio_key',suffix=c('_requested',''))|>dplyr::mutate(ratio=.data$ratio_requested)|>dplyr::select(-dplyr::any_of('ratio_requested'))
 write.csv(sh,file.path(out_dir,'RF_6class_ratio_perturb_share_CI_v13.csv'),row.names=FALSE)
}

summarise_robust<-function(dat,weighting){
 if(weighting=='distinct_ratios')dat<-dat|>dplyr::arrange(.data$m_D,.data$m_B)|>dplyr::distinct(.data$trend_period,.data$zone,.data$class6,.data$support_mode,.data$ratio_key,.data$block_size_km,.keep_all=TRUE)
 dat|>dplyr::group_by(.data$trend_period,.data$zone,.data$class6,.data$support_mode)|>
  dplyr::group_modify(function(.x,.y){
   n<-nrow(.x);sig<-.x$ci_excludes_zero%in%TRUE;pos<-sig & is.finite(.x$delta_low) & .x$delta_low>0;neg<-sig & is.finite(.x$delta_high) & .x$delta_high<0
   sign<-if(!sum(sig))'ns'else if(any(pos)&&!any(neg))'+'else if(any(neg)&&!any(pos))'-'else'mixed'
   b<-dist_con|>dplyr::filter(.data$trend_period==.y$trend_period,.data$zone==.y$zone,.data$class6==.y$class6,.data$support_mode==.y$support_mode,abs(.data$ratio-1)<1e-12);bv<-b$delta_metric[is.finite(b$delta_metric)];bs<-if(!length(bv))'0/NA'else if(mean(bv)>0)'+'else if(mean(bv)<0)'-'else'0/NA'
   finite_delta<-is.finite(.x$delta_metric)
   same<-if(bs=='+').x$delta_metric>0 else if(bs=='-').x$delta_metric<0 else rep(FALSE,n)
   same_valid<-same[finite_delta]
   dv<-.x$delta_metric[finite_delta]
   data.frame(robustness_weighting=weighting,n_distinct_ratios=dplyr::n_distinct(.x$ratio_key),n_grid_cells=nrow(dplyr::distinct(.x,.data$m_D,.data$m_B)),
    n_block_sizes=dplyr::n_distinct(.x$block_size_km),n_tests=n,n_sig=sum(sig),n_sig_pos=sum(pos),n_sig_neg=sum(neg),sig_fraction=if(n>0)sum(sig)/n else NA_real_,
    baseline_sign=bs,sign_when_sig=sign,sign_matches_baseline=if(bs=='0/NA')NA else sign==bs,
    sign_stable_fraction=if(length(same_valid))mean(same_valid)else NA_real_,
    sign_stable_all=if(length(same_valid))all(same_valid)else NA,
    direction_robust=n>0 && sum(sig)==n && sign%in%c('+','-') && sign==bs,
    mostly_direction_robust=n>0 && sum(sig)/n>=MOSTLY_MIN && sign%in%c('+','-') && sign==bs,
    delta_min=if(length(dv))min(dv)else NA_real_,delta_max=if(length(dv))max(dv)else NA_real_,stringsAsFactors=FALSE)
  })|>dplyr::ungroup()
}
rob_grid<-summarise_robust(grid,'grid_cells');rob_dist<-summarise_robust(grid,'distinct_ratios')
write.csv(rob_grid,file.path(out_dir,'RF_6class_ratio_perturb_contrast_robustness_v13.csv'),row.names=FALSE)
write.csv(rob_dist,file.path(out_dir,'RF_6class_ratio_perturb_contrast_robustness_distinct_ratios_v13.csv'),row.names=FALSE)
v12_msg('[ratio_robust] done: %s; universe=%s; distinct ratios=%d',out_dir,BLOCK_UNIVERSE,nrow(distinct))
