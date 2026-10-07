run_bootstrap_coverage_v56 <- function(project_dir=getwd(),B_cal=9999L,B_outer=200L,B_boot=400L,rho=.3,min_segment=5L,alpha=.05,true_break=2010,noise_sd=.3){
  source(file.path(project_dir,"R","utils.R"),local=TRUE);source(file.path(project_dir,"R","core_break.R"),local=TRUE);source(file.path(project_dir,"R","ar1_calibration.R"),local=TRUE)
  years<-2001:2020;x<-years-mean(years);tx<-true_break-mean(years);null<-simulate_supf_null(years,min_segment,min_segment,rho,B_cal,919,min(2000L,B_cal));set_rng("L'Ecuyer-CMRG",920)
  selected<-covered<-logical(B_outer);ci_width<-rep(NA_real_,B_outer)
  for(o in seq_len(B_outer)){
    e<-as.numeric(simulate_ar1_matrix(1,length(years),rho,noise_sd));y<-2-.2*x+.6*pmax(0,x-tx)+e;fit<-fit_break_core(y,years,min_segment,min_segment)
    selected[o]<-isTRUE(fit$testable)&&mc_pvalue(fit$supF,null)<=alpha&&isTRUE(fit$direction_ok);if(!selected[o])next
    sig<-sqrt(max(fit$hinge_sse,0)/max(length(years)-3L,1L));tau<-numeric()
    for(b in seq_len(B_boot)){yb<-fit$fitted+as.numeric(simulate_ar1_matrix(1,length(years),rho,sig));fb<-fit_break_core(yb,years,min_segment,min_segment);if(isTRUE(fb$testable)&&mc_pvalue(fb$supF,null)<=alpha&&isTRUE(fb$direction_ok))tau<-c(tau,fb$break_year)}
    if(length(tau)>=max(30L,ceiling(.1*B_boot))){q<-safe_quantile_discrete(tau,c(.025,.975));covered[o]<-q[1]<=true_break&&q[2]>=true_break;ci_width[o]<-q[2]-q[1]}
  }
  data.frame(true_break=true_break,rho=rho,noise_sd=noise_sd,B_outer=B_outer,B_boot=B_boot,outer_selection_rate=mean(selected),conditional_CI_coverage=if(any(selected&is.finite(ci_width)))mean(covered[selected&is.finite(ci_width)])else NA_real_,median_CI_width=if(any(is.finite(ci_width)))median(ci_width,na.rm=TRUE)else NA_real_,n_CI_evaluated=sum(selected&is.finite(ci_width)),stringsAsFactors=FALSE)
}
