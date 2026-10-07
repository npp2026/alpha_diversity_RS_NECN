# Pure statistical helpers; explicit denominator and resampling contracts.
fit_ood_reference <- function(training,predictors,ridge=1e-4,prob=.975) {
  x<-as.matrix(training[,predictors,drop=FALSE]);storage.mode(x)<-"double"
  x<-x[apply(x,1,function(v)all(is.finite(v) & v != -9999)),,drop=FALSE]
  if(nrow(x)<max(3L,ncol(x)+1L))stop("Insufficient complete training rows for OOD")
  center<-colMeans(x);s<-apply(x,2,sd);s[!is.finite(s)|s==0]<-1
  z<-sweep(sweep(x,2,center,"-"),2,s,"/")
  cv<-stats::cov(z)+diag(ridge,ncol(z));condition<-kappa(cv,exact=TRUE)
  if(!is.finite(condition)||condition>1e8) {
    ee<-eigen(cv,symmetric=TRUE);tol<-max(ee$values)*.Machine$double.eps*ncol(cv)
    iv<-ifelse(ee$values>tol,1/ee$values,0)
    inverse<-ee$vectors%*%(iv*t(ee$vectors));inverse_method<-"pseudoinverse"
  }else {inverse<-solve(cv);inverse_method<-"ordinary"}
  md2<-rowSums((z%*%inverse)*z)
  list(predictors=predictors,center=center,scale=s,min=apply(x,2,min),max=apply(x,2,max),
    inverse=inverse,threshold_md2=unname(quantile(md2,prob)),n_training=nrow(x),
    ridge=ridge,probability=prob,condition_number=condition,inverse_method=inverse_method)
}
predict_ood <- function(reference,x) {
  x<-as.matrix(x)
  if(!is.null(colnames(x))) {
    if(anyDuplicated(colnames(x))||!setequal(colnames(x),reference$predictors))stop("Predictor names mismatch")
    x<-x[,reference$predictors,drop=FALSE]
  }
  if(ncol(x)!=length(reference$predictors))stop("Predictor dimension mismatch")
  valid<-apply(x,1,function(v)all(is.finite(v) & v != -9999))
  range_pass<-rep(FALSE,nrow(x));md2<-rep(NA_real_,nrow(x))
  if(any(valid))range_pass[valid]<-apply(sweep(x[valid,,drop=FALSE],2,reference$min,">="),1,all)&
    apply(sweep(x[valid,,drop=FALSE],2,reference$max,"<="),1,all)
  if(any(range_pass)) {
    z<-sweep(sweep(x[range_pass,,drop=FALSE],2,reference$center,"-"),2,reference$scale,"/")
    md2[range_pass]<-rowSums((z%*%reference$inverse)*z)
  }
  data.frame(valid=valid,range_pass=range_pass,md2=md2,
    md_pass=range_pass&is.finite(md2)&md2<=reference$threshold_md2)
}
# Each row is a region-nested spatial block. Every column uses the SAME draws.
block_ratio_bootstrap <- function(numer,denom,region,B=2000L,seed=42L) {
  numer<-as.matrix(numer);denom<-as.matrix(denom)
  stopifnot(identical(dim(numer),dim(denom)),length(region)==nrow(numer),B>=2L)
  if(anyNA(region)||any(!nzchar(as.character(region))))stop("Missing region labels")
  if(length(B)!=1L||!is.finite(B)||B<2||B!=as.integer(B))stop("B must be an integer >=2")
  if(any(!is.finite(numer))||any(!is.finite(denom))||any(denom<0))stop("Invalid block sufficient statistics")
  groups<-split(seq_len(nrow(numer)),region)
  if(any(lengths(groups)<2L))stop("Need at least two spatial clusters per region for bootstrap")
  set.seed(seed);out<-matrix(NA_real_,B,ncol(numer))
  for(b in seq_len(B)) {
    ix<-unlist(lapply(groups,function(z)z[sample.int(length(z),length(z),replace=TRUE)]),use.names=FALSE)
    den<-colSums(denom[ix,,drop=FALSE]);v<-colSums(numer[ix,,drop=FALSE]);out[b,den>0]<-v[den>0]/den[den>0]
  }
  out
}
percentile_ci <- function(x) {
  x<-x[is.finite(x)];if(length(x)<2L)return(c(lower=NA_real_,upper=NA_real_))
  setNames(as.numeric(quantile(x,c(.025,.975))),c("lower","upper"))
}
regional_block_bootstrap <- function(numer,denom,region,B=2000L,seed=42L) {
  numer<-as.matrix(numer);denom<-as.matrix(denom)
  stopifnot(identical(dim(numer),dim(denom)),length(region)==nrow(numer))
  if(anyNA(region)||any(!nzchar(as.character(region))))stop("Missing region labels")
  if(length(B)!=1L||!is.finite(B)||B<2||B!=as.integer(B))stop("B must be an integer >=2")
  if(any(!is.finite(numer))||any(!is.finite(denom))||any(denom<0))stop("Invalid block sufficient statistics")
  groups<-split(seq_len(nrow(numer)),region)
  if(length(groups)!=5L||any(lengths(groups)<2L))stop("Five regions with >=2 blocks each are required")
  labels<-c("Overall",names(groups));p<-ncol(numer)
  out<-array(NA_real_,c(B,p,6L),dimnames=list(NULL,colnames(numer),labels))
  ratio<-function(ix){d<-colSums(denom[ix,,drop=FALSE]);s<-colSums(numer[ix,,drop=FALSE]);ifelse(d>0,s/d,NA_real_)}
  set.seed(seed)
  for(b in seq_len(B)) {
    sampled<-lapply(groups,function(z)z[sample.int(length(z),length(z),replace=TRUE)])
    out[b,,1]<-ratio(unlist(sampled,use.names=FALSE))
    for(g in seq_along(groups))out[b,,g+1L]<-ratio(sampled[[g]])
  }
  estimates<-rbind(Overall=ratio(seq_len(nrow(numer))),do.call(rbind,lapply(groups,ratio)))
  colnames(estimates)<-colnames(numer)
  list(estimates=estimates,replicates=out,labels=labels)
}
# Uniform sample without replacement, streamed via independent priority keys.
priority_sample_update <- function(old,chunk,n) {
  if(!nrow(chunk))return(old)
  chunk$.priority<-runif(nrow(chunk))
  z<-if(is.null(old))chunk else rbind(old,chunk)
  z[order(z$.priority)[seq_len(min(n,nrow(z)))],,drop=FALSE]
}
prepare_mean_1km <- function(masked,domain,template,min_valid=.8) {
  # Padding is required when a target-cell centre lies outside a partially
  # overlapping native extent; otherwise GDAL can drop the valid edge overlap.
  if(!terra::same.crs(masked,template)||!terra::compareGeom(masked,domain,stopOnError=FALSE))stop("Grid/CRS mismatch")
  if(length(min_valid)!=1L||!is.finite(min_valid)||min_valid<=0||min_valid>1)stop("min_valid must be in (0,1]")
  if(any(c(terra::nlyr(masked),terra::nlyr(domain),terra::nlyr(template))!=1L))stop("Single-layer rasters required")
  template<-terra::ifel(is.finite(template)&template!=0,1,NA)
  # Values outside the forest domain cannot count as valid forest support.
  domain<-terra::ifel(is.finite(domain)&domain>0,1,0)
  masked<-terra::ifel(domain==1&is.finite(masked),masked,NA)
  padded<-terra::extend(masked,template)
  dom<-terra::extend(domain,template,fill=0)
  valid<-terra::ifel(is.finite(padded),1,0)
  mean<-terra::resample(padded,template,method="average")
  fraction<-terra::resample(valid,template,method="average")/terra::resample(dom,template,method="average")
  fraction<-terra::ifel(!is.na(template),fraction,NA)
  value<-terra::ifel(is.finite(fraction)&fraction>=min_valid&!is.na(template),mean,NA)
  list(mean=value,valid_fraction=fraction)
}
