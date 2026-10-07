#!/usr/bin/env Rscript
# ==============================================================================
# RF_6class state-trajectory preflight validation.
# Fails early on malformed configuration, ambiguous input discovery, incomplete
# trend/RF support, raster geometry mismatch, missing CRS, or invalid region data.
# ==============================================================================
options(stringsAsFactors=FALSE)
if(getRversion()<'4.1.0')stop('R >= 4.1.0 is required.',call.=FALSE)
core_pkgs<-c('terra','sf','dplyr','tidyr')
required_pkgs<-core_pkgs
missing_pkgs<-required_pkgs[!vapply(required_pkgs,requireNamespace,logical(1),quietly=TRUE)]
if(length(missing_pkgs))stop('Missing required R package(s): ',paste(missing_pkgs,collapse=', '),call.=FALSE)
suppressPackageStartupMessages({library(terra);library(sf)})
get_script_dir<-function(){a<-commandArgs(FALSE);f<-a[grep('^--file=',a)];if(length(f))dirname(normalizePath(sub('^--file=','',f[1]),winslash='/',mustWork=FALSE))else getwd()}
SCRIPT_DIR<-get_script_dir();source(file.path(SCRIPT_DIR,'v12_posthoc_utils.R'))

needed_scripts<-c('selftest.R','01_compute_state_trajectory_1km.R','02_area_shares_CI.R','03_cross_metric_contrasts_CI.R',
 '04_RF_threshold_sensitivity.R','05_relative_bias_sensitivity.R','06_cross_metric_agreement_TableS11.R',
 '07_build_publication_data.R','validate_outputs.R','v9_parallel_utils.R','v9_trend_package_utils.R','v12_posthoc_utils.R')
miss_scripts<-needed_scripts[!file.exists(file.path(SCRIPT_DIR,needed_scripts))]
if(length(miss_scripts))stop('Package is incomplete; missing script(s): ',paste(miss_scripts,collapse=', '),call.=FALSE)

DATA_DIR<-normalizePath(v12_env_str('DATA_DIR','.'),winslash='/',mustWork=TRUE)
OUT_ROOT<-v12_make_abs_path(v12_env_str('OUT_ROOT_NAME','outputs_RF_6class_fig6_minimal'),DATA_DIR)
CELL_SIZE_M<-v12_env_int('CELL_SIZE_M',1000,1)
if(!is.finite(CELL_SIZE_M)||CELL_SIZE_M<1L)stop('CELL_SIZE_M must be a positive integer number of meters.',call.=FALSE)
YEARS<-v12_validate_year_vector(v12_env_expr('YEARS_R',2001:2020),'YEARS_R',3L)
TREND_PERIODS<-list(trend_2005_2020=2005:2020,trend_2001_2020=2001:2020)
RUN_TREND<-as.character(v12_env_expr('RUN_TREND_PERIODS_R',c('trend_2005_2020')))
if(any(!RUN_TREND%in%names(TREND_PERIODS)))stop('Unknown RUN_TREND_PERIODS_R value(s): ',paste(setdiff(RUN_TREND,names(TREND_PERIODS)),collapse=','),call.=FALSE)
TREND_FRAC<-v12_env_num('TREND_MIN_VALID_YEAR_FRAC',0.80,0,1,FALSE,TRUE)
MIN_VALID_FRAC<-v12_env_num('MIN_VALID_FRAC',0.80,0,1,FALSE,TRUE)
OBS_MIN_FRAC<-v12_env_num('OBS_YEARS_MIN_FRAC',0.80,0,1,FALSE,TRUE)
RF_THRESH<-v12_env_num('RF_THRESH',0.80,0,Inf,FALSE,TRUE)
ALLOW_PARTIAL_OBS<-v12_env_bool('ALLOW_PARTIAL_OBS_YEARS',TRUE)
INPUT_NEGATIVE_TOL<-v12_env_num('INPUT_NEGATIVE_TOL',-1e-10,-Inf,0,TRUE,TRUE)
ALLOW_AMBIG<-v12_env_bool('ALLOW_AMBIGUOUS_INPUT_FILES',FALSE)
RUN_PERIODS<-as.character(v12_env_expr('RUN_PERIODS_R','2016_2020_mean'))
if(!identical(RUN_PERIODS,'2016_2020_mean'))stop('Minimal package supports only RUN_PERIODS_R=2016_2020_mean.',call.=FALSE)
if(nzchar(v12_env_str('POT_REF_YEARS_R','')))stop('v1.6 requires static age-100 q95; remove POT_REF_YEARS_R.')
PRED_DIR<-v12_make_abs_path(v12_env_str('POTENTIAL_DIR','potential_q95'),DATA_DIR);MOUNTAIN<-file.path(DATA_DIR,'NE_Mountain_Output','NE_Mountain_Regions_All.shp')
METRICS<-c('Rich_tree','Shannon_wiener')

use_prof<-v12_env_bool('USE_PROFESSIONAL_TREND_PACKAGES',TRUE);strict_prof<-v12_env_bool('STRICT_PROFESSIONAL_PACKAGES',FALSE)
prof<-c('trend','modifiedmk','mutoss');prof_ok<-vapply(prof,requireNamespace,logical(1),quietly=TRUE)
if(use_prof){message('Optional package availability: ',paste(sprintf('%s=%s',prof,prof_ok),collapse=', '));if(strict_prof&&!all(prof_ok))stop('STRICT_PROFESSIONAL_PACKAGES=TRUE but missing: ',paste(prof[!prof_ok],collapse=', '),call.=FALSE)}

obs_dirs<-function(){raw<-v12_env_str('OBS_ANNUAL_DIRS','');parts<-if(nzchar(raw))trimws(unlist(strsplit(raw,';|,')))else c('mean','1km','obs','obs1km','annual','.')
 d<-vapply(parts[nzchar(parts)],v12_make_abs_path,character(1),base_dir=DATA_DIR);unique(d[dir.exists(d)]) }
find_obs<-function(stem,year){
 dirs<-obs_dirs();nms<-c(sprintf('%s_%d_1km.tif',stem,year),sprintf('%s_%d.tif',stem,year),sprintf('%s_%d_mean_1km.tif',stem,year),sprintf('%s_%d_mean.tif',stem,year),sprintf('%s_mean_%d_1km.tif',stem,year),sprintf('%s_mean_%d.tif',stem,year))
 exact<-unique(unlist(lapply(dirs,function(d)file.path(d,nms)),use.names=FALSE));exact<-exact[file.exists(exact)]
 fuzzy<-unlist(lapply(dirs,function(d)list.files(d,pattern=paste0('.*',stem,'.*',year,'.*\\.tif$'),full.names=TRUE,recursive=FALSE)),use.names=FALSE)
 fuzzy<-fuzzy[!grepl('support_|quantiles|uncertainty|OOD_stack|valid_fraction|Pred_|q0\\.(95|99)|Trend|Mann|Slope|MK|2001_2005|2016_2020|2001_2020|2005_2020',basename(fuzzy),ignore.case=TRUE)]
 cand<-sort(unique(normalizePath(c(exact,fuzzy),winslash='/',mustWork=TRUE)))
 if(length(cand)>1L&&!ALLOW_AMBIG)stop('Multiple observed candidates for ',stem,' ',year,': ',paste(cand,collapse='; '),call.=FALSE)
 if(length(cand)>1L)warning('Multiple observed candidates for ',stem,' ',year,'; using deterministic first path: ',cand[1],call.=FALSE)
 if(length(cand))cand[1]else NA_character_
}

errs<-character()
if(!dir.exists(PRED_DIR))errs<-c(errs,paste('Missing q95 potential folder:',PRED_DIR))
if(!file.exists(MOUNTAIN))errs<-c(errs,paste('Missing mountain vector:',MOUNTAIN))
obs_by_metric<-list()
for(stem in METRICS){
 fs<-vapply(YEARS,function(y)find_obs(stem,y),character(1));names(fs)<-YEARS;obs_by_metric[[stem]]<-fs
 if(any(is.na(fs))&&!ALLOW_PARTIAL_OBS)errs<-c(errs,paste0(stem,' is missing observed years while ALLOW_PARTIAL_OBS_YEARS=FALSE: ',paste(YEARS[is.na(fs)],collapse=',')))
 for(tp in RUN_TREND){yrs<-TREND_PERIODS[[tp]];n<-sum(!is.na(fs[as.character(yrs)]));need<-max(3L,ceiling(length(yrs)*TREND_FRAC));if(n<need)errs<-c(errs,sprintf('%s has %d/%d files for %s; need at least %d.',stem,n,length(yrs),tp,need))}
 stage_n<-sum(!is.na(fs[as.character(2016:2020)]));stage_need<-max(1L,ceiling(5*OBS_MIN_FRAC));if(stage_n<stage_need)errs<-c(errs,sprintf('%s has %d/5 observed files for RF stage mean; need %d.',stem,stage_n,stage_need))
 pf<-file.path(PRED_DIR,sprintf('Q95_%s_age100_1km.tif',stem))
 if(!file.exists(pf))errs<-c(errs,paste('Missing static q95:',pf))
 message(stem,': observed=',sum(!is.na(fs)),'/',length(YEARS),'; static potential=',file.exists(pf))
}
if(length(errs))stop('Preflight failed before geometry checks:\n',paste(errs,collapse='\n'),call.=FALSE)

all_obs<-unlist(obs_by_metric,use.names=FALSE);all_obs<-all_obs[!is.na(all_obs)]
ref_file<-all_obs[1];ref<-terra::rast(ref_file);rr<-terra::res(ref)
if(isTRUE(terra::is.lonlat(ref)))errs<-c(errs,'Reference observed raster is lon/lat; projected meter CRS is required.')
if(any(!is.finite(rr))||any(rr<=0)||abs(rr[1]-rr[2])>max(1e-6,1e-6*max(rr)))errs<-c(errs,'Reference raster has invalid/non-square resolution.')
if(any(abs(rr-1000)>0.01)||CELL_SIZE_M!=1000L)stop('v1.6 requires pre-aligned 1-km observed inputs; use module 08 for 30-m aggregation.')
fact<-CELL_SIZE_M/rr[1];if(!is.finite(fact)||abs(fact-round(fact))>1e-6||fact<1)errs<-c(errs,sprintf('CELL_SIZE_M=%d must be an integer multiple of source resolution %.6f and cannot be finer than the source grid.',CELL_SIZE_M,rr[1]))
check_geom<-function(file,label){r<-tryCatch(terra::rast(file),error=function(e)e);if(inherits(r,'error'))return(paste(label,'cannot be read:',file,r$message));if(terra::nlyr(r)!=1L)return(paste(label,'must be single-layer:',file));if(isTRUE(terra::is.lonlat(r)))return(paste(label,'is lon/lat:',file));if(!terra::compareGeom(r,ref,stopOnError=FALSE,crs=TRUE,rowcol=TRUE,ext=TRUE,res=TRUE))return(paste(label,'is not aligned with reference:',file));g<-tryCatch(terra::global(r,c('min','max'),na.rm=TRUE),error=function(e)NULL);if(is.null(g)||any(!is.finite(as.numeric(g[1,]))))return(paste(label,'has no finite range:',file));if(as.numeric(g[1,'min'])<INPUT_NEGATIVE_TOL)return(paste(label,'contains negative values below tolerance:',file));NULL}
for(stem in METRICS){for(ff in obs_by_metric[[stem]][!is.na(obs_by_metric[[stem]])]){e<-check_geom(ff,paste(stem,'observed'));if(length(e))errs<-c(errs,e)};pf<-file.path(PRED_DIR,sprintf('Q95_%s_age100_1km.tif',stem));for(ff in pf[file.exists(pf)]){e<-check_geom(ff,paste(stem,'potential'));if(length(e))errs<-c(errs,e)}}

mv<-tryCatch(terra::vect(MOUNTAIN),error=function(e)e)
if(inherits(mv,'error'))errs<-c(errs,paste('Mountain vector cannot be read:',mv$message))else{
 if(nrow(mv)!=5L)errs<-c(errs,paste('Mountain vector must contain exactly five region polygons; found',nrow(mv)))
 cr<-tryCatch(terra::crs(mv),error=function(e)'');if(is.null(cr)||!nzchar(cr))errs<-c(errs,'Mountain vector has no CRS.')else{
  mv2<-tryCatch(if(!identical(terra::crs(mv),terra::crs(ref)))terra::project(mv,terra::crs(ref))else mv,error=function(e)e)
  if(inherits(mv2,'error'))errs<-c(errs,paste('Mountain projection failed:',mv2$message))else{ea<-as.vector(terra::ext(mv2));eb<-as.vector(terra::ext(ref));if(!isTRUE(ea[1]<=eb[2]&&ea[2]>=eb[1]&&ea[3]<=eb[4]&&ea[4]>=eb[3]))errs<-c(errs,'Mountain vector does not overlap raster extent.')}}
 vals<-terra::values(mv);nm<-names(vals);pick<-function(x){i<-match(tolower(x),tolower(nm));if(is.na(i))NULL else vals[[i]]};code<-pick('region_code');rom<-pick('Roman_ID');cn<-pick('Reg_CN');if(is.null(code)&&is.null(rom)&&is.null(cn))errs<-c(errs,'Mountain vector needs region_code, Roman_ID, or Reg_CN; feature-order assignment is disabled.')else{mp<-c(HLJDXAL='I',HLJXXAL='II',HLJCBS='III',JLSCBS='IV',LNCBS='V');cnmp<-c('黑龙江大兴安岭'='I','黑龙江小兴安岭'='II','黑龙江长白山'='III','吉林省长白山'='IV','辽宁省长白山'='V');z<-if(!is.null(code))unname(mp[trimws(as.character(code))])else if(!is.null(rom))trimws(as.character(rom))else unname(cnmp[trimws(as.character(cn))]);if(length(z)!=5L||any(is.na(z))||anyDuplicated(z)||!setequal(z,c('I','II','III','IV','V')))errs<-c(errs,paste('Mountain attributes do not map one-to-one to I-V:',paste(z,collapse=','))) }
}

# Verify output location is creatable/writable without deleting existing results.
dir.create(OUT_ROOT,recursive=TRUE,showWarnings=FALSE);probe<-file.path(OUT_ROOT,paste0('.write_test_',Sys.getpid()));ok<-tryCatch({writeLines('ok',probe);unlink(probe);TRUE},error=function(e)FALSE);if(!ok)errs<-c(errs,paste('Output root is not writable:',OUT_ROOT))
message('Reference observed raster: ',ref_file);message('Reference resolution: ',paste(rr,collapse=' x '),'; aggregation factor=',round(fact));message('RF threshold=',RF_THRESH,'; MIN_VALID_FRAC=',MIN_VALID_FRAC,'; potential=static age100 q95')
if(length(errs))stop('RF_6class state-trajectory setup validation failed:\n',paste(errs,collapse='\n'),call.=FALSE)
message('RF_6class state-trajectory setup validation passed.')
