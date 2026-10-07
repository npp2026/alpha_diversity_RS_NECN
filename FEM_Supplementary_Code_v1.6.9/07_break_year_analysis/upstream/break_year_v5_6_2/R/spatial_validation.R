cluster_iou <- function(a,b){a<-as.logical(a);b<-as.logical(b);ok<-!is.na(a)&!is.na(b);u<-sum(a[ok]|b[ok]);if(!u)return(NA_real_);sum(a[ok]&b[ok])/u}
cluster_iou_expected <- function(a,b){a<-as.logical(a);b<-as.logical(b);ok<-!is.na(a)&!is.na(b);if(!any(ok))return(NA_real_);pa<-mean(a[ok]);pb<-mean(b[ok]);d<-pa+pb-pa*pb;if(!is.finite(d)||d<=0)return(NA_real_);pa*pb/d}
cluster_iou_adjusted <- function(a,b){obs<-cluster_iou(a,b);ex<-cluster_iou_expected(a,b);if(!is.finite(obs)||!is.finite(ex)||ex>=1)return(NA_real_);(obs-ex)/(1-ex)}

binary_iou_rasters <- function(a,b){
  assert_packages("terra"); if(!terra::compareGeom(a,b,stopOnError=FALSE))stop("IoU rasters have different geometry.",call.=FALSE)
  bs<-terra::blocks(a);counts<-c(n_ok=0,n_a=0,n_b=0,n_intersection=0,n_union=0)
  for(i in seq_len(bs$n)){va<-terra::values(a,row=bs$row[[i]],nrows=bs$nrows[[i]],mat=FALSE);vb<-terra::values(b,row=bs$row[[i]],nrows=bs$nrows[[i]],mat=FALSE)
    ok<-is.finite(va)&is.finite(vb);aa<-va[ok]==1;bb<-vb[ok]==1;counts["n_ok"]<-counts["n_ok"]+length(aa);counts["n_a"]<-counts["n_a"]+sum(aa);counts["n_b"]<-counts["n_b"]+sum(bb);counts["n_intersection"]<-counts["n_intersection"]+sum(aa&bb);counts["n_union"]<-counts["n_union"]+sum(aa|bb)}
  if(counts["n_ok"]<=0||counts["n_union"]<=0)return(data.frame(raw_iou=NA_real_,expected_iou=NA_real_,adjusted_iou=NA_real_,n_common=counts["n_ok"],n_primary_recovery=counts["n_a"],n_other_recovery=counts["n_b"],n_intersection=counts["n_intersection"],n_union=counts["n_union"]))
  raw<-counts["n_intersection"]/counts["n_union"];pa<-counts["n_a"]/counts["n_ok"];pb<-counts["n_b"]/counts["n_ok"];den<-pa+pb-pa*pb;ex<-if(den>0)pa*pb/den else NA_real_;adj<-if(is.finite(ex)&&ex<1)(raw-ex)/(1-ex)else NA_real_
  data.frame(raw_iou=unname(raw),expected_iou=unname(ex),adjusted_iou=unname(adj),n_common=unname(counts["n_ok"]),n_primary_recovery=unname(counts["n_a"]),n_other_recovery=unname(counts["n_b"]),n_intersection=unname(counts["n_intersection"]),n_union=unname(counts["n_union"]))
}

break_year_agreement_rasters <- function(primary_recovery, primary_break, other_recovery, other_break) {
  assert_packages("terra")
  xs <- list(primary_recovery, primary_break, other_recovery, other_break)
  for (i in 2:length(xs)) if (!terra::compareGeom(xs[[1L]], xs[[i]], stopOnError = FALSE)) stop("Temporal agreement rasters have different geometry.", call. = FALSE)
  bs <- terra::blocks(primary_recovery); hist_abs <- numeric(); n <- 0; exact <- 0; within1 <- 0; abs_sum <- 0
  for (i in seq_len(bs$n)) {
    a <- terra::values(primary_recovery,row=bs$row[[i]],nrows=bs$nrows[[i]],mat=FALSE)
    ay <- terra::values(primary_break,row=bs$row[[i]],nrows=bs$nrows[[i]],mat=FALSE)
    b <- terra::values(other_recovery,row=bs$row[[i]],nrows=bs$nrows[[i]],mat=FALSE)
    by <- terra::values(other_break,row=bs$row[[i]],nrows=bs$nrows[[i]],mat=FALSE)
    ok <- is.finite(a)&is.finite(b)&a==1&b==1&is.finite(ay)&is.finite(by)
    if (!any(ok)) next
    d <- abs(as.integer(round(by[ok])) - as.integer(round(ay[ok]))); n <- n + length(d); exact <- exact + sum(d==0L); within1 <- within1 + sum(d<=1L); abs_sum <- abs_sum + sum(d)
    tb <- table(d); need <- max(as.integer(names(tb))) + 1L; if (length(hist_abs) < need) hist_abs <- c(hist_abs, rep(0, need - length(hist_abs)))
    for (nm in names(tb)) hist_abs[as.integer(nm)+1L] <- hist_abs[as.integer(nm)+1L] + as.numeric(tb[[nm]])
  }
  if (!n) return(data.frame(n_common_recovery_year=0,exact_year_agreement=NA_real_,within1_year_agreement=NA_real_,mean_abs_year_shift=NA_real_,median_abs_year_shift=NA_real_))
  vals <- seq_along(hist_abs)-1L; keep <- hist_abs>0
  data.frame(n_common_recovery_year=n,exact_year_agreement=exact/n,within1_year_agreement=within1/n,mean_abs_year_shift=abs_sum/n,
             median_abs_year_shift=weighted_discrete_quantile_from_counts(vals[keep],hist_abs[keep],.5))
}

weighted_discrete_quantile_from_counts <- function(values, counts, prob) {
  o <- order(values); values <- values[o]; counts <- counts[o]; cs <- cumsum(counts); target <- prob * sum(counts)
  values[[which(cs >= target)[[1L]]]]
}

connected_cluster_summary <- function(recovery,break_year,cfg,out_dir,stability_raster=NULL){
  assert_packages("terra");ensure_dir(out_dir);p<-terra::patches(terra::ifel(recovery==1,1,NA),directions=as.integer(cfg$spatial$connectivity),zeroAsNA=TRUE)
  p<-terra::writeRaster(p,file.path(out_dir,"primary_clusters.tif"),overwrite=isTRUE(cfg$runtime$overwrite),datatype="INT4S")
  area<-terra::cellSize(p,unit="km");a<-terra::zonal(area,p,fun="sum",na.rm=TRUE);n<-terra::zonal(terra::ifel(!is.na(p),1,NA),p,fun="sum",na.rm=TRUE)
  names(a)<-c("cluster_id","area_km2");names(n)<-c("cluster_id","n_pixels")
  tab<-merge(a,n,by="cluster_id",all=TRUE)
  ct<-terra::crosstab(c(p,break_year),long=TRUE)
  if(!is.null(ct)&&nrow(ct)){
    names(ct)[1:3]<-c("cluster_id","break_year","count");ct<-ct[is.finite(ct$cluster_id)&is.finite(ct$break_year)&ct$count>0,,drop=FALSE]
    if(nrow(ct)){
      ss<-split(ct,ct$cluster_id);stats<-do.call(rbind,lapply(ss,function(z){mx<-max(z$count);mode<-min(z$break_year[z$count==mx]);q25<-weighted_discrete_quantile_from_counts(z$break_year,z$count,.25);q50<-weighted_discrete_quantile_from_counts(z$break_year,z$count,.50);q75<-weighted_discrete_quantile_from_counts(z$break_year,z$count,.75);data.frame(cluster_id=z$cluster_id[[1L]],modal_break_year=mode,median_break_year=q50,q25_break_year=q25,q75_break_year=q75,IQR_break_year=q75-q25)}))
      tab<-merge(tab,stats,by="cluster_id",all=TRUE)
    }
  }
  if(!is.null(stability_raster)&&terra::compareGeom(stability_raster,p,stopOnError=FALSE)){
    st<-terra::zonal(stability_raster,p,fun="mean",na.rm=TRUE);names(st)<-c("cluster_id","mean_bootstrap_recovery_frequency");tab<-merge(tab,st,by="cluster_id",all=TRUE)
  }
  write_csv(tab,file.path(out_dir,"cluster_summary.csv"));list(patches=p,summary=tab)
}

block_jackknife_summary <- function(recovery,break_year,block_size_m,out_file){
  assert_packages("terra");if(terra::is.lonlat(recovery))stop("block_jackknife requires a projected metric CRS.",call.=FALSE)
  wkt<-terra::crs(recovery);if(!grepl('CS[Cartesian',wkt,fixed=TRUE)||!grepl('LENGTHUNIT["metre",1',wkt,fixed=TRUE))stop("block_jackknife block_size_m requires a projected CRS whose Cartesian axis unit is metre.",call.=FALSE)
  cells<-cells_from_mask(recovery);if(!length(cells)){write_csv(data.frame(),out_file);return(data.frame())}
  yy<-as.numeric(terra::extract(break_year,cells,raw=TRUE));xy<-terra::xyFromCell(recovery,cells);keep<-is.finite(yy);yy<-yy[keep];xy<-xy[keep,,drop=FALSE]
  if(!length(yy)){write_csv(data.frame(),out_file);return(data.frame())};block<-paste(floor(xy[,1]/block_size_m),floor(xy[,2]/block_size_m),sep="_");blocks<-unique(block);rows<-vector("list",length(blocks))
  full_mode<-modal_value(yy);full_median<-stats::median(yy)
  for(i in seq_along(blocks)){z<-yy[block!=blocks[[i]]];md<-modal_value(z);me<-if(length(z))stats::median(z)else NA_real_;rows[[i]]<-data.frame(removed_block=blocks[[i]],n_remaining=length(z),modal_break_year=md,median_break_year=me,full_modal_break_year=full_mode,full_median_break_year=full_median,modal_changed=is.finite(md)&&md!=full_mode,median_shift=me-full_median,stringsAsFactors=FALSE)}
  tab<-do.call(rbind,rows);write_csv(tab,out_file);tab
}

freq_sizes <- function(patches, id_name, n_name){
  f<-terra::freq(patches);if(is.null(f)||!nrow(f))return(data.frame());value_col<-if("value"%in%names(f))"value" else names(f)[ncol(f)-1L];count_col<-if("count"%in%names(f))"count" else names(f)[ncol(f)]
  z<-data.frame(id=as.numeric(f[[value_col]]),n=as.numeric(f[[count_col]]));names(z)<-c(id_name,n_name);z[is.finite(z[[id_name]]),,drop=FALSE]
}

match_clusters_iou <- function(primary_patches,other_patches){
  assert_packages("terra");z<-terra::crosstab(c(primary_patches,other_patches),long=TRUE);if(is.null(z)||!nrow(z))return(data.frame())
  names(z)[1:3]<-c("primary_cluster","other_cluster","intersection");z<-z[is.finite(z$primary_cluster)&is.finite(z$other_cluster),];if(!nrow(z))return(data.frame())
  # Denominator uses FULL patch sizes, including non-overlap. The older implementation
  # summed only intersections and therefore biased IoU upward.
  pa<-freq_sizes(primary_patches,"primary_cluster","primary_n");oa<-freq_sizes(other_patches,"other_cluster","other_n")
  z<-merge(merge(z,pa,by="primary_cluster"),oa,by="other_cluster");z$iou<-z$intersection/(z$primary_n+z$other_n-z$intersection);z<-z[order(z$primary_cluster,-z$iou),]
  z[!duplicated(z$primary_cluster),c("primary_cluster","other_cluster","intersection","primary_n","other_n","iou")]
}

run_spatial_validation <- function(cfg,primary,boundary,ctx){
  if(!isTRUE(cfg$spatial$connected_clusters)&&!isTRUE(cfg$spatial$scenario_iou$enabled))return(NULL);root<-ensure_dir(file.path(ctx$root,"spatial"));out<-list()
  for(response in names(primary$results)){rdir<-ensure_dir(file.path(root,response));pr<-primary$results[[response]];sf<-file.path(ctx$root,"uncertainty","pixel",paste0(response,"_bootstrap_recovery_frequency.tif"));stab<-if(file.exists(sf))terra::rast(sf)else NULL;cl<-if(isTRUE(cfg$spatial$connected_clusters))connected_cluster_summary(pr$recovery,pr$layers$break_year,cfg,rdir,stab)else NULL
    jack<-if(isTRUE(cfg$spatial$block_jackknife$enabled))block_jackknife_summary(pr$recovery,pr$layers$break_year,as.numeric(cfg$spatial$block_jackknife$block_size_m),file.path(rdir,"block_jackknife.csv"))else NULL
    ioutab<-NULL;matchtab<-NULL
    if(isTRUE(cfg$spatial$scenario_iou$enabled)&&length(boundary)){rows<-list();matches<-list();for(nm in names(boundary)){oth<-boundary[[nm]]$results[[response]];x<-binary_iou_rasters(pr$recovery,oth$recovery);ta<-break_year_agreement_rasters(pr$recovery,pr$layers$break_year,oth$recovery,oth$layers$break_year);x<-cbind(x,ta);x$scenario<-boundary[[nm]]$id;x$min_segment<-oth$summary$min_segment[[1L]];rows[[length(rows)+1L]]<-x
        if(isTRUE(cfg$spatial$scenario_iou$cluster_matching)&&!is.null(cl)){op<-terra::patches(terra::ifel(oth$recovery==1,1,NA),directions=as.integer(cfg$spatial$connectivity),zeroAsNA=TRUE);m<-match_clusters_iou(cl$patches,op);if(nrow(m)){m$scenario<-boundary[[nm]]$id;m$min_segment<-oth$summary$min_segment[[1L]];matches[[length(matches)+1L]]<-m}}}
      ioutab<-do.call(rbind,rows);write_csv(ioutab,file.path(rdir,"scenario_iou.csv"));if(length(matches)){matchtab<-do.call(rbind,matches);write_csv(matchtab,file.path(rdir,"cluster_scenario_matching.csv"))}}
    out[[response]]<-list(clusters=cl,jackknife=jack,scenario_iou=ioutab,cluster_matching=matchtab)}
  out
}
