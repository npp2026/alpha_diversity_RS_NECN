run_smoke_test_v56 <- function(project_dir=getwd(),B=1999L){
  source(file.path(project_dir,"R","utils.R"),local=TRUE);source(file.path(project_dir,"R","core_break.R"),local=TRUE);source(file.path(project_dir,"R","ar1_calibration.R"),local=TRUE)
  years<-2001:2020;null<-simulate_supf_null(years,5,5,.3,B,99,min(500L,B));set_rng("L'Ecuyer-CMRG",7);x<-years-mean(years);tau<-2010;tx<-tau-mean(years);e<-as.numeric(simulate_ar1_matrix(1,length(years),.3,.15));y<-5-.2*x+.55*pmax(0,x-tx)+e
  f<-fit_break_core(y,years,5,5);p<-mc_pvalue(f$supF,null);stopifnot(isTRUE(f$testable),is.finite(p),p>=1/(B+1),p<=1);cat("Smoke test: break=",f$break_year,", SupF=",signif(f$supF,4),", p=",signif(p,4),"\n",sep="");TRUE
}
