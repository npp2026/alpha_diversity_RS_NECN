#!/usr/bin/env Rscript
# ==============================================================================
# Post-run acceptance checks for RF_6class state-trajectory.
# Cross-validates cell table, class raster, area summaries, bootstrap diagnostics,
# exact six-class closure, baseline reconstruction, final-data tables and ratio run.
# ==============================================================================
options(stringsAsFactors=FALSE)
if(getRversion()<'4.1.0')stop('R >= 4.1.0 is required.',call.=FALSE)
suppressPackageStartupMessages({library(terra);library(dplyr)})
get_script_dir<-function(){a<-commandArgs(FALSE);f<-a[grep('^--file=',a)];if(length(f))dirname(normalizePath(sub('^--file=','',f[1]),winslash='/',mustWork=FALSE))else getwd()}
SCRIPT_DIR<-get_script_dir();source(file.path(SCRIPT_DIR,'v12_posthoc_utils.R'))
DATA_DIR<-normalizePath(v12_env_str('DATA_DIR','.'),winslash='/',mustWork=TRUE)
OUT_ROOT<-v12_make_abs_path(v12_env_str('OUT_ROOT_NAME','outputs_RF_6class_fig6_minimal'),DATA_DIR)
GROUP<-v12_env_str('ANALYSIS_GROUP','trend_2005_2020');RUN_RATIO<-v12_env_bool('RUN_RATIO_ROBUSTNESS',TRUE)
CELL_SIZE_M<-v12_env_int('CELL_SIZE_M',1000,1);if(!is.finite(CELL_SIZE_M)||CELL_SIZE_M<1L)stop('CELL_SIZE_M must be a positive integer number of meters.',call.=FALSE)
CELL_AREA<-v12_env_num('CELL_AREA_KM2',(CELL_SIZE_M/1000)^2,0,Inf,FALSE,TRUE);MIN_BOOT_VALID<-v12_env_num('MIN_BOOT_VALID_FRACTION',0.95,0,1,TRUE,TRUE)
MAP_FILE<-sprintf('RF_6class_%dm.tif',CELL_SIZE_M)
need<-function(p){if(!file.exists(p))stop('Acceptance check missing required file: ',p,call.=FALSE);p}
check_cols<-function(d,cols,label){m<-setdiff(cols,names(d));if(length(m))stop(label,' missing columns: ',paste(m,collapse=', '),call.=FALSE)}
check_exact_classes<-function(x,groups,label){key<-do.call(interaction,c(x[groups],list(drop=TRUE,lex.order=TRUE)));sp<-split(as.character(x$class6),key);bad<-vapply(sp,function(z)length(z)!=length(V12_CLASS6_LEVELS)||anyDuplicated(z)>0L||!setequal(z,V12_CLASS6_LEVELS),logical(1));if(any(bad))stop(label,' has incomplete or duplicated class sets in ',sum(bad),' groups.',call.=FALSE)}

cell_file<-need(file.path(OUT_ROOT,GROUP,'RF_6class_cell_table.csv'));d<-v12_validate_cell_table(read.csv(cell_file,check.names=FALSE),TRUE,cell_file);d<-v12_ensure_area(d,CELL_AREA)
map_file<-need(file.path(OUT_ROOT,GROUP,MAP_FILE));r<-terra::rast(map_file);rr<-terra::res(r)
if(length(rr)<2L||any(!is.finite(rr))||any(abs(rr-CELL_SIZE_M)>1e-6))stop('Main raster is not ',CELL_SIZE_M,' m: ',paste(rr,collapse=' x '),call.=FALSE)
if(isTRUE(terra::is.lonlat(r))||!nzchar(terra::crs(r)))stop('Main raster requires a projected CRS.',call.=FALSE)
req_bands<-c('rich_class6_code','shannon_class6_code');if(!all(req_bands%in%names(r)))stop('Main raster lacks class-code bands.',call.=FALSE)
if(any(d$cell>terra::ncell(r)))stop('Cell table contains cell IDs outside raster.',call.=FALSE)
rv<-terra::values(r[[req_bands]],mat=TRUE);dom<-sort(unique(as.numeric(rv[is.finite(rv)])))
if(any(!dom%in%1:6))stop('Raster class bands contain codes outside 1..6: ',paste(setdiff(dom,1:6),collapse=','),call.=FALSE)
for(metric in c('rich','shannon')){
 col<-paste0(metric,'_class6');band<-paste0(metric,'_class6_code');expected<-unname(V12_CLASS6_CODE[as.character(d[[col]])]);actual<-as.numeric(rv[d$cell,band])
 if(any((is.na(expected)!=is.na(actual))|(!is.na(expected)&expected!=actual)))stop('Cell-table/raster class mismatch for ',metric,': ',sum((is.na(expected)!=is.na(actual))|(!is.na(expected)&expected!=actual)),' rows.',call.=FALSE)
}

summary_file<-need(file.path(OUT_ROOT,GROUP,'RF_6class_area_summary.csv'));asum<-read.csv(summary_file,check.names=FALSE)
check_cols(asum,c('metric','zone','class6','n_cells','area_km2','prop'),'Area summary');check_exact_classes(asum,c('metric','zone'),'Area summary')
for(metric in c('Richness','Shannon'))for(zone in V12_REG_ORDER){
 dz<-if(zone=='Overall')d else d[v12_zone_value(d$roman)==zone,,drop=FALSE];col<-if(metric=='Richness')'rich_class6'else'shannon_class6';metric_key<-if(metric=='Richness')'rich'else'shannon';valid<-as.character(dz[[col]])%in%V12_CLASS6_LEVELS
 for(cl in V12_CLASS6_LEVELS){a<-asum[as.character(asum$metric)==metric_key&as.character(asum$zone)==zone&as.character(asum$class6)==cl,,drop=FALSE];if(nrow(a)!=1L)stop('Area summary row uniqueness failure.',call.=FALSE);hit<-valid&as.character(dz[[col]])==cl
   if(a$n_cells!=sum(hit)||abs(a$area_km2-sum(dz$area_km2[hit],na.rm=TRUE))>1e-7)stop('Area summary disagrees with cell table: ',metric,'/',zone,'/',cl,call.=FALSE)
 }
}

area<-read.csv(need(file.path(OUT_ROOT,GROUP,'RF_6class_area_block_bootstrap_CI_v9.csv')),check.names=FALSE)
check_cols(area,c('metric','zone','class6','block_size_km','prop','prop_low','prop_high','bootstrap_valid_fraction','n_blocks_universe'),'Area CI')
check_exact_classes(area,c('metric','zone','block_size_km'),'Area CI')
if(!50%in%area$block_size_km)stop('Area CI lacks required 50-km result.',call.=FALSE)
if(any(area$n_cells>0&area$bootstrap_valid_fraction<MIN_BOOT_VALID,na.rm=TRUE))stop('Area CI bootstrap-valid fraction below threshold.',call.=FALSE)

md<-read.csv(need(file.path(OUT_ROOT,'metric_diff_CI_v9','RF_6class_metric_diff_within_period_Shannon_minus_Richness_v9.csv')),check.names=FALSE)
check_cols(md,c('trend_period','zone','class6','block_size_km','support_mode','prop_shannon','prop_richness','delta_metric','delta_low','delta_high','bootstrap_valid_fraction'),'Metric contrast')
mdg<-md[as.character(md$trend_period)==GROUP,,drop=FALSE];if(!nrow(mdg))stop('Metric contrast has no rows for ',GROUP,call.=FALSE)
check_exact_classes(mdg,c('trend_period','zone','block_size_km','support_mode'),'Metric contrast')
if(!all(c(50,75,100)%in%unique(mdg$block_size_km)))stop('Metric contrast lacks 50/75/100-km results.',call.=FALSE)
closure<-mdg|>dplyr::group_by(.data$zone,.data$block_size_km,.data$support_mode)|>dplyr::summarise(n=dplyr::n(),n_common=max(.data$n_common_cells,na.rm=TRUE),rich_area=max(.data$rich_classified_area_km2,na.rm=TRUE),shannon_area=max(.data$shannon_classified_area_km2,na.rm=TRUE),ds=sum(.data$delta_metric),ps=sum(.data$prop_shannon),pr=sum(.data$prop_richness),.groups='drop')
app<-ifelse(as.character(closure$support_mode)=='common',is.finite(closure$n_common)&closure$n_common>0,is.finite(closure$rich_area)&closure$rich_area>0&is.finite(closure$shannon_area)&closure$shannon_area>0)
if(any(closure$n!=6L)||any(app&(abs(closure$ds)>1e-9|abs(closure$ps-1)>1e-9|abs(closure$pr-1)>1e-9)))stop('Six-class contrast completeness/closure failed.',call.=FALSE)
if(any(mdg$n_common_cells>0&mdg$bootstrap_valid_fraction<MIN_BOOT_VALID,na.rm=TRUE))stop('Metric bootstrap-valid fraction below threshold.',call.=FALSE)

rfqc<-read.csv(need(file.path(OUT_ROOT,'metric_diff_CI_v9','RF_6class_RFthreshold_baseline_reproduction_QC_v13.csv')),check.names=FALSE)
if(any(rfqc$classified_match_fraction<0.999999|rfqc$support_match_fraction<0.999999,na.rm=TRUE))stop('RF-threshold baseline reconstruction failed.',call.=FALSE)
final<-file.path(OUT_ROOT,'fig6_final_data');required_final<-c('fig_source/Fig6C_regional_6class_stackedbar_source.csv','fig_source/Fig6D_focal4_paired_contrasts_50km_source.csv','fig_source/FigS10_RFthreshold_sensitivity_focal4_50km_source.csv','tables/TableS11A_cross_classification_matrix_long.csv','tables/TableS11B_cross_metric_agreement_summary.csv','tables/TableS0a_focal4_coverage_QC.csv','tables/TableS0b_class_contrast_closure_QC.csv',file.path('maps',MAP_FILE))
for(f in required_final)need(file.path(final,f))
if(!identical(unname(tools::md5sum(map_file)),unname(tools::md5sum(file.path(final,'maps',MAP_FILE)))))stop('Final-data map is stale/different from analysis map.',call.=FALSE)
cl<-read.csv(file.path(final,'tables','TableS0b_class_contrast_closure_QC.csv'),check.names=FALSE);if(!all(cl$closure_ok%in%TRUE))stop('Final-data closure QC failed.',call.=FALSE)
ag<-read.csv(file.path(final,'tables','TableS11B_cross_metric_agreement_summary.csv'),check.names=FALSE);for(nm in c('exact_6class_agreement','hl_state_agreement','trend_state_agreement','area_weighted_cohen_kappa'))if(any(ag[[nm]]< -1|ag[[nm]]>1,na.rm=TRUE))stop('Agreement/kappa outside valid range: ',nm,call.=FALSE)

if(RUN_RATIO){
 rq<-read.csv(need(file.path(OUT_ROOT,'ratio_robustness_v13','RF_6class_ratio_baseline_reproduction_QC_v13.csv')),check.names=FALSE)
 if(any(rq$classified_match_fraction<0.999999|rq$support_match_fraction<0.999999,na.rm=TRUE))stop('Ratio baseline reconstruction failed.',call.=FALSE)
 need(file.path(OUT_ROOT,'ratio_robustness_v13','RF_6class_ratio_perturb_contrast_robustness_distinct_ratios_v13.csv'))
}
summary<-data.frame(check='RF_6class_acceptance',status='PASS',analysis_group=GROUP,n_cell_rows=nrow(d),rich_classified=sum(as.character(d$rich_class6)%in%V12_CLASS6_LEVELS),shannon_classified=sum(as.character(d$shannon_class6)%in%V12_CLASS6_LEVELS),n_metric_contrast_rows=nrow(mdg),ratio_checked=RUN_RATIO,stringsAsFactors=FALSE)
write.csv(summary,file.path(OUT_ROOT,'RF_6class_v13_acceptance_summary.csv'),row.names=FALSE)
v12_msg('[acceptance] PASS: cells=%d; contrasts=%d; ratio=%s',nrow(d),nrow(mdg),RUN_RATIO)
