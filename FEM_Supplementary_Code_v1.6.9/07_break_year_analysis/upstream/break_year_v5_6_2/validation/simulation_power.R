run_power_simulation_v56 <- function(project_dir=getwd(),B_cal=9999L,B_sim=2000L,rho=.3,min_segment=5L,alpha=.05,break_years=c(2005,2008,2010,2012,2015),slope_changes=c(.2,.4,.6,.8)){
  source(file.path(project_dir,"R","utils.R"),local=TRUE);source(file.path(project_dir,"R","core_break.R"),local=TRUE);source(file.path(project_dir,"R","ar1_calibration.R"),local=TRUE);years<-2001:2020;x<-years-mean(years);null<-simulate_supf_null(years,min_segment,min_segment,rho,B_cal,313,min(2000L,B_cal));rows<-list()
  for(tau in break_years)for(delta in slope_changes){tx<-tau-mean(years);det<-rec<-logical(B_sim);err<-rep(NA_real_,B_sim);set_rng("L'Ecuyer-CMRG",seed_from_key(314,tau,delta))
    for(b in seq_len(B_sim)){e<-as.numeric(simulate_ar1_matrix(1,length(years),rho,.3));y<-2-.2*x+delta*pmax(0,x-tx)+e;f<-fit_break_core(y,years,min_segment,min_segment);p<-mc_pvalue(f$supF,null);det[b]<-p<=alpha;rec[b]<-det[b]&&isTRUE(f$direction_ok);if(det[b])err[b]<-f$break_year-tau}
    rows[[length(rows)+1L]]<-data.frame(true_break=tau,slope_change=delta,detection_power=mean(det),recovery_power=mean(rec),break_rmse=if(any(is.finite(err)))sqrt(mean(err^2,na.rm=TRUE))else NA_real_)}
  do.call(rbind,rows)
}
