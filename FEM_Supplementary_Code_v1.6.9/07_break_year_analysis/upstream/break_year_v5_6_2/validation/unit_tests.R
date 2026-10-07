run_unit_tests_v56 <- function(project_dir=getwd()){
  source(file.path(project_dir,"R","utils.R"),local=TRUE);source(file.path(project_dir,"R","parallel_runtime.R"),local=TRUE);source(file.path(project_dir,"R","core_break.R"),local=TRUE);source(file.path(project_dir,"R","ar1_calibration.R"),local=TRUE);source(file.path(project_dir,"R","fdr_discrete.R"),local=TRUE);source(file.path(project_dir,"R","spatial_validation.R"),local=TRUE);source(file.path(project_dir,"R","sensitivity.R"),local=TRUE);source(file.path(project_dir,"validation","rho_stress_gate.R"),local=TRUE)
  years<-2001:2020;ms<-5L;cb<-function(m)years[candidate_indices(length(years),m,m)]
  ylin<-2+0.3*(years-mean(years));f1<-fit_break_core(ylin,years,ms,ms);stopifnot(isTRUE(f1$testable),identical(f1$fit_status,"perfect_linear"),isTRUE(all.equal(f1$supF,0)))
  tau<-2010;x<-years-mean(years);tx<-tau-mean(years);yhinge<-10-.5*x+1.2*pmax(0,x-tx);f2<-fit_break_core(yhinge,years,ms,ms);stopifnot(isTRUE(f2$testable),is.infinite(f2$supF),f2$break_year==tau,isTRUE(f2$direction_ok))
  f3<-fit_break_core(rep(1,length(years)),years,ms,ms,tie_rule="earliest");stopifnot(f3$break_year==min(cb(ms)),f3$tie_count>1L)
  null<-sort(c(1,2,2,3));stopifnot(mc_exceed_ge(2,null)==3L,abs(mc_pvalue(2,null)-0.8)<1e-15,mc_exceed_ge(Inf,c(1,2,Inf))==1L)
  # Response-specific MC B is a numerical-resolution selector only; the frozen
  # baseline remains the fallback for every response without an explicit map.
  cfgB<-list(supf=list(mc_B=200000L,mc_B_by_response=list(Shannon=1000000L)),data=list(responses=list(SR=list(),Shannon=list())))
  stopifnot(supf_mc_B_for_response(cfgB,"SR")==200000L,supf_mc_B_for_response(cfgB,"Shannon")==1000000L,
            identical(unname(supf_mc_B_by_response(cfgB)),c(200000L,1000000L)))
  k<-c(0L,0L,1L,1L,2L,5L,9L);stopifnot(validate_discrete_fdr(k,9,"BH"),validate_discrete_fdr(k,9,"BY"))
  stopifnot(seed_from_key(1,"A",22)==seed_from_key(1,"A",22),seed_from_key(1,"A",22)!=seed_from_key(1,"A",23))
  a<-c(1,1,0,0,1);b<-c(1,0,0,1,1);stopifnot(abs(cluster_iou(a,b)-0.5)<1e-12,is.finite(cluster_iou_adjusted(a,b)))
  stopifnot(identical(cb(5),2005:2015),identical(cb(4),2004:2016),identical(cb(3),2003:2017),identical(cb(2),2002:2018))
  # Rho sensitivity can include a response-specific diagnostic vector without
  # hard-coding real-subset values into the publication configuration.
  emp<-list(SR=0.07,Shannon=0.21);src<-list(empirical=emp,linear_diagnostic=list(SR=0.31,Shannon=0.405))
  rv<-parse_rho_sensitivity(c(0,"empirical","linear_diagnostic",0.4),emp,src)
  stopifnot(all(c("rho0","empirical","linear_diagnostic","rho0.4") %in% names(rv)),
            abs(rv$linear_diagnostic$SR-0.31)<1e-12,abs(rv$linear_diagnostic$Shannon-0.405)<1e-12)
  # Rho stress diagnostic policy is intentionally audit-only: robust stress
  # can GO, abundance sensitivity yields CONDITIONAL_GO, and integrity failure
  # yields NO_GO rather than silently switching the primary rho estimator.
  gate_policy<-list(robust_retention_min=.5,collapse_retention_max=.1,mc_floor_clear_multiplier=10,modal_year_shift_warn=1)
  gate_cmp<-data.frame(response=c("SR","Shannon"),primary_rho=c(.07,.21),linear_diagnostic_rho=c(.31,.405),
    rho_abs_difference=c(.24,.195),primary_recovery_breaks=c(100,100),stress_recovery_breaks=c(70,55),
    recovery_retention_vs_primary=c(.70,.55),stress_fdr_cutoff_to_mc_floor=c(100,80),modal_break_year_shift=c(0,1))
  g1<-rho_stress_gate_classify(gate_cmp,gate_policy);stopifnot(g1$overall_gate=="GO",g1$severity=="RHO_STRESS_ROBUST")
  gate_cmp$recovery_retention_vs_primary<-c(.30,.05);gate_cmp$stress_recovery_breaks<-c(30,5)
  g2<-rho_stress_gate_classify(gate_cmp,gate_policy);stopifnot(g2$overall_gate=="CONDITIONAL_GO",g2$severity=="SEVERE_RHO_MODEL_SENSITIVITY",any(grepl("nearly collapses",g2$conditions,fixed=TRUE)))
  gate_cmp$linear_diagnostic_rho[[1L]]<-NA_real_;g3<-rho_stress_gate_classify(gate_cmp,gate_policy);stopifnot(g3$overall_gate=="NO_GO",g3$severity=="GATE_INTEGRITY_FAILURE")
  gate_cmp$linear_diagnostic_rho[[1L]]<-.31;g4<-rho_stress_gate_classify(gate_cmp,gate_policy,context="full");stopifnot(any(grepl("formal study inference",g4$conditions,fixed=TRUE)),!any(grepl("subset recovery",g4$conditions,fixed=TRUE)))
  # Scale invariance: shrinking y must not convert a noisy finite fit into an exact fit.
  set.seed(4);yn<-2+0.2*(years-2001)+rnorm(length(years),sd=.03);fa<-fit_break_core(yn,years,ms,ms);fb<-fit_break_core(yn*1e-6,years,ms,ms)
  stopifnot(is.finite(fa$supF),is.finite(fb$supF),abs(fa$supF-fb$supF)<1e-6)
  Ym<-rbind(yn,yn+seq_along(yn)*1e-4);Fm1<-supf_for_matrix(Ym,years,ms,ms);Fm2<-supf_for_matrix(Ym*1e-6,years,ms,ms);stopifnot(all(is.finite(Fm1)),all(is.finite(Fm2)),max(abs(Fm1-Fm2))<1e-5)
  cat("v5.6 unit tests passed.\n");TRUE
}
