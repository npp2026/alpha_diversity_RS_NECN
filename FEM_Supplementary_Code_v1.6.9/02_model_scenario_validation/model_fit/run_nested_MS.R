#!/usr/bin/env Rscript
args<-commandArgs(TRUE)
if(length(args)<3L)stop("Usage: Rscript run_nested_MS.R DATA.csv OUT_DIR SHARED_KNNDM_FOLDS.csv [FEATURE_GROUPS.csv]")
module<-dirname(normalizePath(sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1])))
source(file.path(module,"nested_core_MS.R"));source(file.path(module,"feature_groups_MS.R"))
source(file.path(module,"..","..","R","input_contracts.R"))
for(p in c("ranger","VSURF","xgboost"))if(!requireNamespace(p,quietly=TRUE))stop("Install ",p)
out<-args[2];if(dir.exists(out)&&length(list.files(out)))stop("Use an empty output directory to avoid mixing versions")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
d<-read.csv(args[1],check.names=FALSE,stringsAsFactors=FALSE,colClasses=c(plot_id="character"))
req<-c("plot_id","plot_year","Lon_Export","Lat_Export","Rich_tree","Shannon_wiener")
if(!all(req%in%names(d)))stop("Missing required columns: ",paste(setdiff(req,names(d)),collapse=","))
matched <- fem_match_plots(d)
write.csv(data.frame(input_row=seq_len(nrow(d)),plot_id=d$plot_id,included=matched$included),
          file.path(out,"survey_year_inclusion.csv"),row.names=FALSE)
d <- matched$data
if(length(args)>=4L&&nzchar(args[4])) {
  g<-read.csv(args[4],stringsAsFactors=FALSE)
  if(!all(c("feature","group")%in%names(g))||anyDuplicated(g$feature)||any(!g$group%in%c("Env","Prod","Het","Temp")))stop("Invalid group dictionary")
  if(!all(g$feature%in%names(d)))stop("Dictionary contains unavailable columns")
  groups<-split(g$feature,g$group)
}else {
  groups<-define_feature_groups(d)[c("Env","Prod","Het","Temp")]
  groups$Env<-unique(c(groups$Env,grep("^(ELEV$|BIO[0-9]+_|SPEI_)",names(d),value=TRUE)))
}
for(g in names(groups))groups[[g]]<-setdiff(groups[[g]],c(req,"Forest_age","PFT","AGB","Forest_type"))
features<-unique(unlist(groups,use.names=FALSE))
fem_validate_groups(groups)
if(any(!vapply(d[,features,drop=FALSE],is.numeric,logical(1))))stop("All predictors must be numeric")
finite<-fem_finite_rows(d,c(features,req[-1]))
write.csv(data.frame(plot_id=d$plot_id,included=finite),file.path(out,"sample_inclusion.csv"),row.names=FALSE)
d<-d[finite,,drop=FALSE]
if(nrow(d)!=3066L)warning("MS has 3066 plots; this input retains ",nrow(d),". Results are a new run.")
write.csv(do.call(rbind,lapply(names(groups),function(g)data.frame(feature=groups[[g]],group=g))),file.path(out,"feature_groups_used.csv"),row.names=FALSE)
f<-read.csv(args[3],stringsAsFactors=FALSE,colClasses=c(plot_id="character"))
if(!all(c("plot_id","knndm_fold")%in%names(f))||anyDuplicated(f$plot_id))stop("Invalid shared kNNDM assignment")
fold<-f$knndm_fold[match(d$plot_id,f$plot_id)]
if(anyNA(fold)||!setequal(unique(fold),1:10))stop("Need complete shared kNNDM folds 1:10; no random fallback")
folds<-list(random=ms_folds(d$plot_id,10L,42L),kNNDM=as.integer(fold))
write.csv(data.frame(plot_id=d$plot_id,random_fold=folds$random,knndm_fold=folds$kNNDM),file.path(out,"shared_outer_folds.csv"),row.names=FALSE)
all_results<-list();records<-list();metrics<-list();counter<-0L
for(resp in c("Rich_tree","Shannon_wiener"))for(design in names(folds)) {
  models<-list()
  for(sc in names(ms_scenarios)) {
    fs<-unique(unlist(groups[ms_scenarios[[sc]]],use.names=FALSE))
    rf<-xgb<-qm<-rep(NA_real_,nrow(d))
    for(k in 1:10) {
      tr<-which(folds[[design]]!=k);te<-which(folds[[design]]==k)
      seed<-42L+10000L*match(resp,c("Rich_tree","Shannon_wiener"))+1000L*match(design,names(folds))+k
      message(resp," ",design," ",sc," outer fold ",k)
      z<-ms_outer(d[tr,,drop=FALSE],d[te,,drop=FALSE],resp,fs,seed)
      rf[te]<-z$predictions$RF$raw;qm[te]<-z$predictions$RF$calibrated;xgb[te]<-z$predictions$XGBoost$raw
      saveRDS(z$audit,file.path(out,paste(resp,design,sc,paste0("fold",k),"audit.rds",sep="_")))
    }
    # Figure2 and all source contrasts use raw outer-test RF OOF predictions.
    # Fold-local QM remains a separately labeled calibration diagnostic.
    models[[sc]]<-list(success=TRUE,predictions=rf,raw_predictions=rf,
      prediction_type="RF_raw",calibrated_predictions=qm,calibrated_prediction_type="RF_QM",observed=d[[resp]],
      coords=d[,c("Lon_Export","Lat_Export")],fold_assignment=folds[[design]],plot_id=d$plot_id)
    for(alg in c("RF_raw","RF_QM","XGBoost_raw")) {
      pr<-switch(alg,RF_raw=rf,RF_QM=qm,XGBoost_raw=xgb);counter<-counter+1L
      records[[counter]]<-data.frame(response=resp,design=design,scenario=sc,algorithm=alg,
        plot_id=d$plot_id,fold=folds[[design]],observed=d[[resp]],prediction=pr)
      metrics[[counter]]<-data.frame(response=resp,design=design,scenario=sc,algorithm=alg,t(ms_metric(d[[resp]],pr)))
    }
  }
  all_results[[paste(resp,design,sep="_")]]<-list(target_var=resp,pft="ALL",cv_method=design,data_type="match_only",model_results=models)
  saveRDS(all_results,file.path(out,"all_results.rds"))
  write.csv(do.call(rbind,records),file.path(out,"nested_OOF_predictions.csv"),row.names=FALSE)
  write.csv(do.call(rbind,metrics),file.path(out,"nested_metrics.csv"),row.names=FALSE)
}
# Independent early -> late cohort evaluation; no test observations used for selection/tuning/QM.
tr<-which(d$plot_year%in%2008:2012);te<-which(d$plot_year%in%2014:2017)
if(length(tr)<50L||length(te)<10L)stop("Temporal validation requires both 2008-2012 and 2014-2017 cohorts")
if(length(intersect(d$plot_id[tr],d$plot_id[te])))stop("Cohort plot IDs overlap")
temporal<-list()
for(resp in c("Rich_tree","Shannon_wiener")) {
  z<-ms_outer(d[tr,,drop=FALSE],d[te,,drop=FALSE],resp,features,17042L)
  saveRDS(z,file.path(out,paste0(resp,"_temporal_cohort_audit.rds")))
  temporal[[resp]]<-data.frame(response=resp,plot_id=d$plot_id[te],year=d$plot_year[te],
    observed=d[[resp]][te],RF_raw=z$predictions$RF$raw,RF_QM=z$predictions$RF$calibrated,
    XGBoost_raw=z$predictions$XGBoost$raw)
}
write.csv(do.call(rbind,temporal),file.path(out,"temporal_cohort_predictions.csv"),row.names=FALSE)
# Table S1 compares uncalibrated algorithms on the same fold-local features.
mt<-do.call(rbind,metrics);raw<-mt[mt$algorithm%in%c("RF_raw","XGBoost_raw"),]
s1<-aggregate(cbind(R2,RMSE,MAE)~response+design+algorithm,raw,mean)
write.csv(s1,file.path(out,"Table_S1_algorithm_means.csv"),row.names=FALSE)
tm<-list()
for(resp in names(temporal))for(alg in c("RF_raw","RF_QM","XGBoost_raw")) {
 z<-temporal[[resp]];tm[[paste(resp,alg)]]<-data.frame(response=resp,algorithm=alg,n_train=length(tr),n_test=length(te),t(ms_metric(z$observed,z[[alg]])))
}
write.csv(do.call(rbind,tm),file.path(out,"temporal_cohort_metrics.csv"),row.names=FALSE)
capture.output(sessionInfo(),file=file.path(out,"sessionInfo.txt"))
writeLines(c("profile=MS_v1.6.9_rawRF_fixed15","inner_folds=5 (implementation choice)","algorithm_comparison=RF_raw versus XGBoost_raw; identical fold-local feature sets",
 "RF_scenario_outputs=RF_raw_outer_test_OOF","RF_QM_outputs=separate_calibration_diagnostics",
 "VSURF_scope=validation_only; production_uses_independent_fixed_15_predictors",
 "results_are_not_historical_manuscript_numbers"),file.path(out,"RUN_COMPLETE.txt"))
