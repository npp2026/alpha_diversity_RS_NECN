#!/usr/bin/env Rscript
# ==============================================================================
# Build final Fig.6 / Fig.S10 / SI data bundle.
# Complete six-class outputs are retained; manuscript-facing contrasts focus on
# H+, H0, L+, and L0. All target contrasts use Shannon minus Richness and the
# common classified support by default.
# ==============================================================================
options(stringsAsFactors=FALSE)
if(getRversion()<'4.1.0')stop('R >= 4.1.0 is required.',call.=FALSE)
suppressPackageStartupMessages({library(dplyr);library(tidyr)})
get_script_dir<-function(){a<-commandArgs(FALSE);f<-a[grep('^--file=',a)];if(length(f))dirname(normalizePath(sub('^--file=','',f[1]),winslash='/',mustWork=FALSE))else getwd()}
SCRIPT_DIR<-get_script_dir();source(file.path(SCRIPT_DIR,'v12_posthoc_utils.R'))

DATA_DIR<-normalizePath(v12_env_str('DATA_DIR','.'),winslash='/',mustWork=TRUE)
OUT_ROOT<-v12_make_abs_path(v12_env_str('OUT_ROOT_NAME','outputs_RF_6class_fig6_minimal'),DATA_DIR)
GROUP<-v12_env_str('ANALYSIS_GROUP','trend_2005_2020')
CELL_SIZE_M<-v12_env_int('CELL_SIZE_M',1000,1)
if(!is.finite(CELL_SIZE_M)||CELL_SIZE_M<1L)stop('CELL_SIZE_M must be a positive integer number of meters; got ',CELL_SIZE_M,call.=FALSE)
CELL_AREA_DEFAULT_KM2<-(CELL_SIZE_M/1000)^2
PRIMARY_SUPPORT<-v12_normalise_support_mode(v12_env_str('TARGET_CONTRAST_SUPPORT','common'))
RUN_RATIO<-v12_env_bool('RUN_RATIO_ROBUSTNESS',TRUE)
CLEAN_FINAL_DATA<-v12_env_bool('CLEAN_FINAL_DATA',TRUE)
OUT<-file.path(OUT_ROOT,'fig6_final_data')
if(CLEAN_FINAL_DATA&&dir.exists(OUT))unlink(OUT,recursive=TRUE,force=TRUE)
if(CLEAN_FINAL_DATA&&dir.exists(OUT))stop('Failed to clean derived final-data directory: ',OUT,call.=FALSE)
for(d in c(OUT,file.path(OUT,'fig_source'),file.path(OUT,'tables'),file.path(OUT,'maps')))dir.create(d,recursive=TRUE,showWarnings=FALSE)
need<-function(p){if(!file.exists(p))stop('Missing required file: ',p,call.=FALSE);p}
zone_order<-c('Overall','I','II','III','IV','V','outside')

# Fig.6C and Table S9: full six-class composition.
area_sum<-read.csv(need(file.path(OUT_ROOT,GROUP,'RF_6class_area_summary.csv')),check.names=FALSE)
if(!'metric_label'%in%names(area_sum))area_sum$metric_label<-ifelse(area_sum$metric=='rich','Richness','Shannon')
area_sum<-area_sum|>dplyr::mutate(prop_pct=100*.data$prop,class6=factor(.data$class6,levels=V12_CLASS6_LEVELS),
 zone=factor(.data$zone,levels=zone_order),metric_label=factor(.data$metric_label,levels=V12_METRIC_LEVELS))|>
 dplyr::arrange(.data$zone,.data$metric_label,.data$class6)
write.csv(area_sum,file.path(OUT,'fig_source','Fig6C_regional_6class_stackedbar_source.csv'),row.names=FALSE)
write.csv(area_sum|>dplyr::filter(as.character(.data$class6)%in%V12_FOCAL4_LEVELS),file.path(OUT,'fig_source','Fig6C_regional_focal4_stackedbar_source.csv'),row.names=FALSE)

area_ci<-read.csv(need(file.path(OUT_ROOT,GROUP,'RF_6class_area_block_bootstrap_CI_v9.csv')),check.names=FALSE)|>
 dplyr::filter(.data$block_size_km==50)|>
 dplyr::mutate(prop_pct=100*.data$prop,prop_low_pct=100*.data$prop_low,prop_high_pct=100*.data$prop_high,
 class6=factor(.data$class6,levels=V12_CLASS6_LEVELS),zone=factor(.data$zone,levels=zone_order),metric=factor(.data$metric,levels=V12_METRIC_LEVELS))|>
 dplyr::arrange(.data$zone,.data$metric,.data$class6)
write.csv(area_ci,file.path(OUT,'tables','TableS9_regional_6class_area_shares_50km_CI.csv'),row.names=FALSE)
write.csv(area_ci|>dplyr::filter(as.character(.data$class6)%in%V12_FOCAL4_LEVELS),file.path(OUT,'tables','TableS9_focal4_regional_area_shares_50km_CI.csv'),row.names=FALSE)

# Table S10 / Fig.6D: same direction for every class, common support by default.
md<-read.csv(need(file.path(OUT_ROOT,'metric_diff_CI_v9','RF_6class_metric_diff_within_period_Shannon_minus_Richness_v9.csv')),check.names=FALSE)|>
 dplyr::filter(.data$trend_period==GROUP)
if('support_mode'%in%names(md))md<-md|>dplyr::filter(as.character(.data$support_mode)==PRIMARY_SUPPORT)
md<-md|>dplyr::mutate(delta_SminusR_pp=100*.data$delta_metric,delta_low_pp=100*.data$delta_low,delta_high_pp=100*.data$delta_high,
 prop_shannon_pct=100*.data$prop_shannon,prop_richness_pct=100*.data$prop_richness,
 class6=factor(.data$class6,levels=V12_CLASS6_LEVELS),zone=factor(.data$zone,levels=zone_order))
s2<-md|>dplyr::filter(.data$block_size_km==50)|>dplyr::arrange(.data$zone,.data$class6)
write.csv(s2,file.path(OUT,'tables','TableS10_class_area_share_differences_50km.csv'),row.names=FALSE)
write.csv(s2|>dplyr::filter(as.character(.data$class6)%in%V12_FOCAL4_LEVELS),file.path(OUT,'tables','TableS10_focal4_area_share_differences_50km.csv'),row.names=FALSE)
all_support_path<-file.path(OUT_ROOT,'metric_diff_CI_v9','RF_6class_metric_diff_within_period_all_support_modes_v13.csv')
if(file.exists(all_support_path)){
  md_all<-read.csv(all_support_path,check.names=FALSE)|>dplyr::filter(.data$trend_period==GROUP,.data$block_size_km==50)
  write.csv(md_all,file.path(OUT,'tables','QC_complete_6class_metric_diff_all_support_modes_v13.csv'),row.names=FALSE)
  write.csv(md_all|>dplyr::filter(as.character(.data$class6)%in%V12_FOCAL4_LEVELS),file.path(OUT,'tables','QC_focal4_metric_diff_all_support_modes_v13.csv'),row.names=FALSE)
}
support_qc_path<-file.path(OUT_ROOT,'metric_diff_CI_v9','RF_6class_target_support_QC_v13.csv')
if(file.exists(support_qc_path))file.copy(support_qc_path,file.path(OUT,'tables',basename(support_qc_path)),overwrite=TRUE)

make_focal_contrast<-function(x){
 x|>dplyr::mutate(contrast=as.character(.data$class6),formula=paste0('p_S(',as.character(.data$class6),') - p_R(',as.character(.data$class6),')'),
 estimate=.data$delta_metric,low=.data$delta_low,high=.data$delta_high,estimate_pp=100*.data$delta_metric,low_pp=100*.data$delta_low,high_pp=100*.data$delta_high,
 contrast_definition='Shannon_minus_Richness')
}
figD<-md|>dplyr::filter(.data$block_size_km==50,as.character(.data$class6)%in%V12_FOCAL4_LEVELS,as.character(.data$zone)%in%c('Overall','I','II','III','IV','V'))|>
 make_focal_contrast()|>dplyr::arrange(factor(.data$zone,levels=c('Overall','I','II','III','IV','V')),factor(.data$class6,levels=V12_FOCAL4_LEVELS))
write.csv(figD,file.path(OUT,'fig_source','Fig6D_focal4_paired_contrasts_50km_source.csv'),row.names=FALSE)

s3<-md|>dplyr::filter(as.character(.data$class6)%in%V12_FOCAL4_LEVELS,as.character(.data$zone)%in%c('Overall','I','II','III','IV','V'))|>
 make_focal_contrast()|>dplyr::arrange(factor(.data$zone,levels=c('Overall','I','II','III','IV','V')),factor(.data$class6,levels=V12_FOCAL4_LEVELS),.data$block_size_km)
write.csv(s3,file.path(OUT,'tables','FigS9_block_size_threshold_support_data.csv'),row.names=FALSE)

rf<-read.csv(need(file.path(OUT_ROOT,'metric_diff_CI_v9','RF_6class_metric_diff_RFthreshold_grid_v9.csv')),check.names=FALSE)|>
 dplyr::filter(.data$trend_period==GROUP,as.character(.data$class6)%in%V12_FOCAL4_LEVELS,as.character(.data$zone)%in%c('Overall','I','II','III','IV','V'))
if('support_mode'%in%names(rf))rf<-rf|>dplyr::filter(as.character(.data$support_mode)==PRIMARY_SUPPORT)
rf<-make_focal_contrast(rf)|>dplyr::arrange(factor(.data$zone,levels=c('Overall','I','II','III','IV','V')),factor(.data$class6,levels=V12_FOCAL4_LEVELS),.data$rf_threshold,.data$block_size_km)
write.csv(rf,file.path(OUT,'tables','FigS9_S10_RFthreshold_sensitivity_all_blocks.csv'),row.names=FALSE)
write.csv(rf|>dplyr::filter(.data$block_size_km==50),file.path(OUT,'fig_source','FigS10_RFthreshold_sensitivity_focal4_50km_source.csv'),row.names=FALSE)
rf_all_path<-file.path(OUT_ROOT,'metric_diff_CI_v9','RF_6class_metric_diff_RFthreshold_grid_all_support_modes_v13.csv')
if(file.exists(rf_all_path)){
  rf_all<-read.csv(rf_all_path,check.names=FALSE)|>dplyr::filter(.data$trend_period==GROUP,as.character(.data$class6)%in%V12_FOCAL4_LEVELS,as.character(.data$zone)%in%c('Overall','I','II','III','IV','V'))
  write.csv(rf_all,file.path(OUT,'tables','FigS9_S10_RFthreshold_sensitivity_all_support_modes.csv'),row.names=FALSE)
}
robp<-file.path(OUT_ROOT,'metric_diff_CI_v9','RF_6class_metric_diff_RFthreshold_direction_robustness_v9.csv')
if(file.exists(robp)){
 rob<-read.csv(robp,check.names=FALSE)|>dplyr::filter(.data$trend_period==GROUP,as.character(.data$class6)%in%V12_FOCAL4_LEVELS,as.character(.data$zone)%in%c('Overall','I','II','III','IV','V'))
 write.csv(rob,file.path(OUT,'tables','RFthreshold_direction_robustness_focal4_Shannon_minus_Richness.csv'),row.names=FALSE)
}

# Same-cell target transition, agreement, focal-four coverage, and closure QC.
cell_file<-need(file.path(OUT_ROOT,GROUP,'RF_6class_cell_table.csv'))
cell<-read.csv(cell_file,check.names=FALSE);cell<-v12_validate_cell_table(cell,FALSE,cell_file);cell_area_km2<-v12_env_num('CELL_AREA_KM2',CELL_AREA_DEFAULT_KM2,0,Inf,FALSE,TRUE);cell<-v12_ensure_area(cell,cell_area_km2)
ta<-v12_target_transition(cell,cell_area_km2)
write.csv(ta$transition,file.path(OUT,'tables','TableS11A_cross_classification_matrix_long.csv'),row.names=FALSE)
write.csv(ta$agreement,file.path(OUT,'tables','TableS11B_cross_metric_agreement_summary.csv'),row.names=FALSE)
write.csv(ta$coverage,file.path(OUT,'tables','TableS0a_focal4_coverage_QC.csv'),row.names=FALSE)
closure<-md|>dplyr::group_by(.data$trend_period,.data$zone,.data$block_size_km,.data$support_mode)|>
 dplyr::summarise(n_classes=dplyr::n_distinct(as.character(.data$class6)),class_set=paste(sort(unique(as.character(.data$class6))),collapse='|'),n_common=max(.data$n_common_cells,na.rm=TRUE),rich_area=max(.data$rich_classified_area_km2,na.rm=TRUE),shannon_area=max(.data$shannon_classified_area_km2,na.rm=TRUE),sum_delta=sum(.data$delta_metric),sum_prop_shannon=sum(.data$prop_shannon),sum_prop_richness=sum(.data$prop_richness),.groups='drop')|>
 dplyr::mutate(applicable=ifelse(as.character(.data$support_mode)=='common',is.finite(.data$n_common)&.data$n_common>0,is.finite(.data$rich_area)&.data$rich_area>0&is.finite(.data$shannon_area)&.data$shannon_area>0),closure_ok=.data$n_classes==6L&.data$class_set==paste(sort(V12_CLASS6_LEVELS),collapse='|')&(!.data$applicable|(abs(.data$sum_delta)<1e-10&abs(.data$sum_prop_shannon-1)<1e-10&abs(.data$sum_prop_richness-1)<1e-10)))
write.csv(closure,file.path(OUT,'tables','TableS0b_class_contrast_closure_QC.csv'),row.names=FALSE)
if(any(!closure$closure_ok))stop('Class-composition closure QC failed; inspect TableS0b_class_contrast_closure_QC.csv.',call.=FALSE)

# Optional ratio-robustness tables are copied into the final table bundle.
ratio_dir<-file.path(OUT_ROOT,'ratio_robustness_v13')
if(RUN_RATIO&&dir.exists(ratio_dir))for(f in list.files(ratio_dir,pattern='\\.csv$',full.names=TRUE))file.copy(f,file.path(OUT,'tables',basename(f)),overwrite=TRUE)

# Map sources: hard guard against the historical 1-km/5-km mismatch.
map_name<-sprintf('RF_6class_%dm.tif',CELL_SIZE_M)
map_src<-need(file.path(OUT_ROOT,GROUP,map_name))
manifest<-data.frame(file=c(map_name,'RF_6class_cell_table.csv'),source=c(map_src,cell_file),exists=TRUE,
 purpose=c(sprintf('Fig.6A/B %dm raster stack',CELL_SIZE_M),'same-cell transition and custom map rebuilding'),stringsAsFactors=FALSE)
write.csv(manifest,file.path(OUT,'maps','Fig6AB_map_file_manifest.csv'),row.names=FALSE)
file.copy(map_src,file.path(OUT,'maps',basename(map_src)),overwrite=TRUE);file.copy(cell_file,file.path(OUT,'maps',basename(cell_file)),overwrite=TRUE)

readme<-c('RF_6class state-trajectory final data','',paste0('Analysis group: ',GROUP),paste0('Target contrast support: ',PRIMARY_SUPPORT),'Contrast direction: Shannon minus Richness for every class.','',
 'Main files:','- fig_source/Fig6C_regional_6class_stackedbar_source.csv','- fig_source/Fig6D_focal4_paired_contrasts_50km_source.csv','- fig_source/FigS10_RFthreshold_sensitivity_focal4_50km_source.csv', '- tables/TableS11A_cross_classification_matrix_long.csv','- tables/TableS11B_cross_metric_agreement_summary.csv','- tables/TableS0a_focal4_coverage_QC.csv','- tables/TableS0b_class_contrast_closure_QC.csv')
writeLines(readme,file.path(OUT,'README_final_data.txt'))
v12_msg('[Fig6 data] written: %s',OUT)
