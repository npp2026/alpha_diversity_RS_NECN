#!/usr/bin/env Rscript
# Uses full native 30-m maps for annual means and pixel-wise trends.
a<-commandArgs(TRUE)
if(length(a)<2L)stop("Usage: Rscript run_regional_30m.R OOD_OUTPUT_DIR OUT_DIR [--profile=diagnostic --sample-per-region=N --bootstrap=B]")
opts<-if(length(a)>2L)a[-c(1L,2L)]else character()
if(any(!grepl("^--(profile=(manuscript|diagnostic)|sample-per-region=[0-9]+|bootstrap=[0-9]+)$",opts)))stop("Unknown regional option")
get_opt<-function(key,default){z<-grep(paste0("^--",key,"="),opts,value=TRUE);if(length(z)>1L)stop("Repeated option: ",key);if(length(z))sub(paste0("^--",key,"="),"",z)else default}
profile<-get_opt("profile","manuscript")
if(profile=="manuscript"&&any(grepl("^--(sample-per-region|bootstrap)=",opts)))stop("Nonstandard counts require --profile=diagnostic")
sample_per_region<-as.integer(get_opt("sample-per-region","200000"));boot_B<-as.integer(get_opt("bootstrap","2000"))
if(!is.finite(sample_per_region)||sample_per_region<2L||!is.finite(boot_B)||boot_B<2L)stop("Counts must be integers >=2")
message("Regional profile=",profile,"; B=",boot_B,"; sample/region=",sample_per_region)
if(profile=="diagnostic")message("DIAGNOSTIC RUN: counts differ from manuscript; do not label as production results")
module<-dirname(normalizePath(sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1])))
source(file.path(module,"core_MS.R"))
source(file.path(module,"..","05_state_trajectory_1km","v9_trend_package_utils.R"))
if(!requireNamespace("terra",quietly=TRUE))stop("Install terra")
input<-normalizePath(a[1],mustWork=TRUE);out<-a[2]
if(dir.exists(out)&&length(list.files(out)))stop("Use an empty output directory")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
write.csv(data.frame(profile=profile,bootstrap_B=boot_B,sample_per_region=sample_per_region),file.path(out,"regional_run_config.csv"),row.names=FALSE)
m<-read.csv(file.path(input,"masked_annual_manifest.csv"),stringsAsFactors=FALSE)
zone<-terra::rast(file.path(input,"blocks50km_region_nested.tif"));region<-terra::rast(file.path(input,"regions_30m.tif"))
domain<-terra::rast(file.path(input,"forest_domain_30m.tif"))
# Establish common block rows, even if a response/year has no usable pixels.
zz<-terra::zonal(region,zone,fun="min",na.rm=TRUE);names(zz)<-c("block","region")
zz<-zz[order(zz$block),];years<-2001:2020;windows<-list(`2005_2020`=2005:2020,`2001_2020`=2001:2020)
annual_rows<-slope_rows<-contrast_rows<-change_rows<-window_rows<-list();means<-list();boots<-list()
for(resp in c("Rich_tree","Shannon_wiener")) {
 mm<-m[m$response==resp,];if(nrow(mm)!=20L||anyDuplicated(mm$year)||!setequal(mm$year,years))stop("Need exactly 20 masked annual maps per response")
 st<-terra::rast(mm$prediction[match(years,mm$year)])
 if(!terra::compareGeom(st,zone,stopOnError=FALSE)||any(abs(terra::res(st)-30)>.01))stop("Native 30-m geometry mismatch")
 colkeys<-c(paste0("y",years),paste0("sen_",names(windows)),"early","late")
 num<-den<-matrix(0,nrow(zz),length(colkeys),dimnames=list(NULL,colkeys))
 add_stat<-function(r,key) {
   z<-terra::zonal(c(r,terra::ifel(is.finite(r),1,0)),zone,fun="sum",na.rm=TRUE)
   j<-match(zz$block,z[[1]]);s<-z[j,2];n<-z[j,3];s[!is.finite(s)]<-0;n[!is.finite(n)]<-0
   num[,key]<<-s;den[,key]<<-n
 }
 for(j in seq_along(years))add_stat(st[[j]],paste0("y",years[j]))
 period_mean<-function(yrs) {
   z<-st[[match(yrs,years)]];n<-terra::app(is.finite(z),sum)
   terra::ifel(n>=ceiling(.8*length(yrs)),terra::app(z,fun="mean",na.rm=TRUE),NA)
 }
 early<-period_mean(2005:2010);late<-period_mean(2016:2020)
 # Relative changes use identical finite spatial support for the two periods.
 paired<-is.finite(early)&is.finite(late)
 add_stat(terra::ifel(paired,early,NA),"early");add_stat(terra::ifel(paired,late,NA),"late")
 means[[resp]]<-terra::writeRaster(late,file.path(out,paste0(resp,"_2016_2020_mean_30m.tif")))
 for(win in names(windows)) {
  yy<-windows[[win]]
  one<-function(v) {
   if(sum(is.finite(v))<ceiling(.8*length(yy)))return(c(sen=NA_real_,p_hr=NA_real_,n_valid=sum(is.finite(v))))
   z<-v9_trend_stats(yy,v,prefer_package=FALSE,prefer_sen_package=FALSE)
   c(sen=unname(z$sen$value),p_hr=as.numeric(z$hr["p"]),n_valid=sum(is.finite(v)))
  }
  fun<-function(v)if(is.matrix(v))t(vapply(seq_len(nrow(v)),function(i)one(v[i,]),numeric(3L)))else one(v)
  z<-terra::app(st[[match(yy,years)]],fun,cores=1L,filename=file.path(out,paste0(resp,"_",win,"_native_trend.tif")))
  names(z)<-c("sen","p_hr","n_valid")
  # Exact BH within response x window. Memory cost is proportional to valid 30-m cells.
  pv<-terra::values(z[["p_hr"]],mat=FALSE);ok<-is.finite(pv);q<-rep(NA_real_,length(pv));q[ok]<-p.adjust(pv[ok],"BH")
  qr<-z[["p_hr"]];terra::values(qr)<-q;terra::writeRaster(qr,file.path(out,paste0(resp,"_",win,"_q_BH.tif")))
  rm(pv,q,ok,qr);gc(verbose=FALSE)
  terra::writeRaster(z[["sen"]],file.path(out,paste0(resp,"_",win,"_Sen_30m.tif")))
  add_stat(z[["sen"]],paste0("sen_",win))
 }
 tab<-data.frame(zz,num,den,check.names=FALSE)
 names(tab)<-c("block","region",paste0("sum_",colkeys),paste0("n_",colkeys))
 write.csv(tab,file.path(out,paste0(resp,"_50km_block_tabulation.csv")),row.names=FALSE)
 b<-regional_block_bootstrap(num,den,zz$region,boot_B,42L);boots[[resp]]<-b
 saveRDS(b,file.path(out,paste0(resp,"_paired_block_bootstrap.rds")))
 for(g in seq_along(b$labels)) {
  label<-b$labels[g]
  for(j in seq_along(years)) {
   key<-paste0("y",years[j]);ci<-percentile_ci(b$replicates[,key,g])
   annual_rows[[paste(resp,g,j)]]<-data.frame(response=resp,region=label,year=years[j],mean=b$estimates[g,key],lower=ci[1],upper=ci[2])
  }
  for(win in names(windows)) {
   key<-paste0("sen_",win);ci<-percentile_ci(b$replicates[,key,g])
   slope_rows[[paste(resp,g,win)]]<-data.frame(response=resp,region=label,window=win,mean_sen=b$estimates[g,key],lower=ci[1],upper=ci[2])
  }
  delta<-b$replicates[,"sen_2005_2020",g]-b$replicates[,"sen_2001_2020",g]
  delta_ci<-percentile_ci(delta)
  window_rows[[paste(resp,g)]]<-data.frame(response=resp,region=label,
    difference=b$estimates[g,"sen_2005_2020"]-b$estimates[g,"sen_2001_2020"],
    lower=delta_ci[1],upper=delta_ci[2],definition="2005_2020 minus 2001_2020; shared block draws")
  v<-100*(b$replicates[,"late",g]/b$replicates[,"early",g]-1);v[b$replicates[,"early",g]<=0]<-NA
  ci<-percentile_ci(v);estimate<-if(is.finite(b$estimates[g,"early"])&&b$estimates[g,"early"]>0)100*(b$estimates[g,"late"]/b$estimates[g,"early"]-1)else NA_real_
  change_rows[[paste(resp,g)]]<-data.frame(response=resp,region=label,percent_change=estimate,lower=ci[1],upper=ci[2],definition="100*(late_mean/early_mean-1); common support")
 }
 for(win in names(windows))for(pair in combn(2:6,2,simplify=FALSE)) {
  key<-paste0("sen_",win);v<-b$replicates[,key,pair[1]]-b$replicates[,key,pair[2]];ci<-percentile_ci(v)
  # Centred two-sided bootstrap test of zero difference.
  est<-b$estimates[pair[1],key]-b$estimates[pair[2],key]
  p<-(1+sum(abs(v-est)>=abs(est),na.rm=TRUE))/(1+sum(is.finite(v)))
  contrast_rows[[paste(resp,win,pair[1],pair[2])]]<-data.frame(response=resp,window=win,region_A=b$labels[pair[1]],region_B=b$labels[pair[2]],difference=est,lower=ci[1],upper=ci[2],p_raw=p)
 }
}
annual<-do.call(rbind,annual_rows);write.csv(annual,file.path(out,"annual_30m_means_pointwise95CI.csv"),row.names=FALSE)
write.csv(do.call(rbind,slope_rows),file.path(out,"regional_mean_pixel_Sen_CI.csv"),row.names=FALSE)
write.csv(do.call(rbind,change_rows),file.path(out,"regional_relative_changes_CI.csv"),row.names=FALSE)
write.csv(do.call(rbind,window_rows),file.path(out,"Table_S6_paired_window_Sen_differences.csv"),row.names=FALSE)
ct<-do.call(rbind,contrast_rows);ct$p_BH<-ave(ct$p_raw,interaction(ct$response,ct$window),FUN=function(x)p.adjust(x,"BH"))
write.csv(ct,file.path(out,"regional_Sen_pairwise_contrasts.csv"),row.names=FALSE)
# Table S5 / Figure 4: exact 200,000 common-support pixels per region.
stack<-c(means$Rich_tree,means$Shannon_wiener,region,domain)
terra::readStart(stack);chunks<-terra::blocks(stack);samples<-vector("list",5L);set.seed(142L)
for(k in seq_len(chunks$n)) {
 v<-terra::readValues(stack,row=chunks$row[k],nrows=chunks$nrows[k],mat=TRUE)
 cell<-seq_len(nrow(v))+(chunks$row[k]-1L)*terra::ncol(stack)
 good<-is.finite(v[,1])&is.finite(v[,2])&is.finite(v[,3])&v[,4]==1
 good[is.na(good)]<-FALSE
 for(g in 1:5) {
   ix<-which(good&v[,3]==g)
   if(length(ix))samples[[g]]<-priority_sample_update(samples[[g]],data.frame(cell=cell[ix],region=g,Rich_tree=v[ix,1],Shannon_wiener=v[ix,2]),sample_per_region)
 }
}
terra::readStop(stack)
if(any(vapply(samples,function(z)if(is.null(z))0L else nrow(z),integer(1))!=sample_per_region))stop("Fewer than the requested valid paired pixels in a region; sample size is never silently reduced")
sample<-do.call(rbind,samples);sample$.priority<-NULL
write.csv(sample,file.path(out,paste0("balanced_",sample_per_region,"_per_region_sample.csv")),row.names=FALSE)
summary<-list()
for(resp in c("Rich_tree","Shannon_wiener"))for(g in 0:5) {
 v<-if(g==0)sample[[resp]]else sample[[resp]][sample$region==g]
 summary[[paste(resp,g)]]<-data.frame(response=resp,region=if(g==0)"Overall"else as.character(g),n=length(v),mean=mean(v),sd=sd(v),median=median(v),q25=unname(quantile(v,.25)),q75=unname(quantile(v,.75)))
}
write.csv(do.call(rbind,summary),file.path(out,"Table_S5_contemporary_balanced_summary.csv"),row.names=FALSE)
capture.output(sessionInfo(),file=file.path(out,"sessionInfo.txt"))
writeLines(c(paste0("Completed native 30-m regional analysis; profile=",profile),paste0("Pointwise ribbons: ",boot_B," region-nested 50-km block resamples."),paste0("Overall annual/block summaries pool block sums; Table S5 pools equal ",sample_per_region,"-pixel regional samples.")),file.path(out,"RUN_COMPLETE.txt"))
