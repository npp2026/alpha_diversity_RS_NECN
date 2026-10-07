run_null_simulation_v56 <- function(project_dir=getwd(),B_cal=9999L,B_sim=5000L,rho_grid=c(0,.2,.4,.6),min_segments=c(5,4,3,2),alpha=.05){
  source(file.path(project_dir,"R","utils.R"),local=TRUE);source(file.path(project_dir,"R","core_break.R"),local=TRUE);source(file.path(project_dir,"R","ar1_calibration.R"),local=TRUE);years<-2001:2020;rows<-list()
  for(ms in min_segments)for(rho in rho_grid){null<-simulate_supf_null(years,ms,ms,rho,B_cal,seed_from_key(11,ms,rho),min(2000L,B_cal));set_rng("L'Ecuyer-CMRG",seed_from_key(12,ms,rho));p<-numeric(B_sim)
    for(b in seq_len(B_sim)){y<-as.numeric(simulate_ar1_matrix(1,length(years),rho,1));f<-fit_break_core(y,years,ms,ms);p[b]<-mc_pvalue(f$supF,null)}
    rows[[length(rows)+1L]]<-data.frame(min_segment=ms,rho=rho,alpha=alpha,empirical_type1=mean(p<=alpha),mc_se=sqrt(alpha*(1-alpha)/B_sim))}
  do.call(rbind,rows)
}
