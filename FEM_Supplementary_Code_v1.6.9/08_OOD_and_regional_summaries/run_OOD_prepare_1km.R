#!/usr/bin/env Rscript
# Input manifest: response,year,prediction,predictors,train_rds,model_rds.
# All paths absolute. Rasters must be aligned at 30 m in the same metre CRS.
a<-commandArgs(TRUE)
if(length(a)!=5L)stop("Usage: Rscript run_OOD_prepare_1km.R MANIFEST.csv FOREST_MASK.tif REGIONS.gpkg TEMPLATE_1KM.tif OUT_DIR")
module<-dirname(normalizePath(sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1])))
source(file.path(module,"core_MS.R"))
for(p in c("terra","sf","ranger"))if(!requireNamespace(p,quietly=TRUE))stop("Install ",p)
out<-a[5];if(dir.exists(out)&&length(list.files(out)))stop("Use an empty output directory")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
m<-read.csv(a[1],stringsAsFactors=FALSE)
if(!all(c("response","year","prediction","predictors","train_rds","model_rds")%in%names(m)))stop("Manifest columns missing")
if(anyDuplicated(paste(m$response,m$year))||nrow(m)!=40L||!setequal(m$response,c("Rich_tree","Shannon_wiener")))stop("Manifest must contain 2 responses x 20 unique years")
for(resp in unique(m$response))if(!setequal(m$year[m$response==resp],2001:2020))stop("Need years 2001:2020")
for(nm in c("prediction","predictors","train_rds","model_rds"))if(any(!file.exists(m[[nm]])))stop("Missing manifest path in ",nm)
forest<-terra::rast(a[2]);template<-terra::rast(a[4])
unit<-tolower(sf::st_crs(terra::crs(forest))$units_gdal)
if(terra::nlyr(forest)!=1L||terra::is.lonlat(forest)||!unit%in%c("metre","meter","metres","meters","m")||any(abs(terra::res(forest)-30)>.01))stop("Forest mask must be a single 30-m layer in metre CRS")
if(terra::nlyr(template)!=1L||!terra::same.crs(forest,template)||any(abs(terra::res(template)-1000)>.01))stop("Template must be a single 1-km layer in the SAME projected CRS")
template<-terra::ifel(is.finite(template)&template!=0,1,NA)
regions<-terra::vect(a[3]);if(!terra::same.crs(regions,forest))regions<-terra::project(regions,terra::crs(forest))
if(!"region_code"%in%names(regions))stop("Regions need region_code: HLJDXAL,HLJXXAL,HLJCBS,JLSCBS,LNCBS")
codes<-c("HLJDXAL","HLJXXAL","HLJCBS","JLSCBS","LNCBS")
regions$region_id<-match(regions$region_code,codes)
regions$Reg_EN<-as.character(regions$region_code)
regions$Roman_ID<-c("I","II","III","IV","V")[regions$region_id]
regions$Reg_CN<-c("黑龙江大兴安岭","黑龙江小兴安岭","黑龙江长白山","吉林省长白山","辽宁省长白山")[regions$region_id]
regions$L2_code<-c("Greater_Khingan","Lesser_Khingan","Changbai","Changbai","Changbai")[regions$region_id]
regions$L2_name<-regions$L2_code
if(nrow(regions)!=5L||anyNA(regions$region_id)||anyDuplicated(regions$region_id))stop("Need exactly one polygon feature per region (multipart is allowed)")
region<-terra::rasterize(regions,forest,field="region_id",filename=file.path(out,"regions_30m.tif"))
domain<-terra::ifel(is.finite(forest)&forest>0&!is.na(region),1,0)
terra::writeRaster(domain,file.path(out,"forest_domain_30m.tif"))
# Fixed projected-coordinate origin shared by all years and responses.
bsize<-50000;nx<-ceiling((terra::xmax(forest)-terra::xmin(forest))/bsize)+1
ny<-ceiling((terra::ymax(forest)-terra::ymin(forest))/bsize)+1
bid<-floor((terra::init(forest,"x")-terra::xmin(forest))/bsize)+nx*floor((terra::init(forest,"y")-terra::ymin(forest))/bsize)+1
zone<-terra::ifel(domain==1,bid+(region-1)*nx*ny,NA)
terra::writeRaster(zone,file.path(out,"blocks50km_region_nested.tif"))
weight<-terra::cellSize(forest,unit="m");refs<-list();ood<-list();manifest<-list()
for(i in seq_len(nrow(m))) {
 row<-m[i,];resp<-row$response;yr<-row$year;message(resp," ",yr)
 if(is.null(refs[[resp]])) {
   if(length(unique(m$train_rds[m$response==resp]))!=1L||length(unique(m$model_rds[m$response==resp]))!=1L)stop("Use one final fitted model and training table per response")
   model<-readRDS(row$model_rds);predictors<-model$forest$independent.variable.names
   if(!length(predictors))stop("Cannot read ranger independent variable names")
   refs[[resp]]<-fit_ood_reference(readRDS(row$train_rds),predictors)
   saveRDS(refs[[resp]],file.path(out,paste0(resp,"_OOD_reference.rds")))
 }
 ref<-refs[[resp]];x<-terra::rast(row$predictors);pred<-terra::rast(row$prediction)
 if(anyDuplicated(names(x))||!all(ref$predictors%in%names(x)))stop("Predictor raster bands do not match model")
 x<-x[[ref$predictors]]
 if(terra::nlyr(pred)!=1L||!terra::compareGeom(x,forest,stopOnError=FALSE)||!terra::compareGeom(pred,forest,stopOnError=FALSE))stop("30-m raster geometry mismatch")
 p<-length(ref$predictors)
 fun<-function(v) {
  if(!is.matrix(v))v<-matrix(v,nrow=1)
  z<-predict_ood(ref,v[,seq_len(p),drop=FALSE])
  valid<-z$valid&is.finite(v[,p+1])&v[,p+2]==1&is.finite(v[,p+3])&v[,p+3]>0
  valid[is.na(valid)]<-FALSE;w<-ifelse(valid,v[,p+3],0)
  passed<-valid&z$md_pass
  cbind(masked=ifelse(passed,v[,p+1],NA_real_),valid_weight=w,
    range_fail_weight=ifelse(valid&!z$range_pass,w,0),
    md_eligible_weight=ifelse(valid&z$range_pass&is.finite(z$md2),w,0),
    md_fail_weight=ifelse(valid&z$range_pass&is.finite(z$md2)&!z$md_pass,w,0))
 }
 dest<-file.path(out,paste0(resp,"_",yr,"_OOD_stack.tif"))
 z<-terra::app(c(x,pred,domain,weight),fun,cores=1,filename=dest,wopt=list(datatype="FLT8S"))
 names(z)<-c("masked","valid_weight","range_fail_weight","md_eligible_weight","md_fail_weight")
 native<-file.path(out,paste0(resp,"_",yr,"_OOD_30m.tif"));terra::writeRaster(z[[1]],native)
 sums<-terra::zonal(z[[2:5]],zone,fun="sum",na.rm=TRUE)
 names(sums)[1]<-"block";sums$region<-((sums$block-1)%/%(nx*ny))+1L
 sums<-sums[is.finite(sums$valid_weight)&sums$valid_weight>0,,drop=FALSE]
 if(!setequal(unique(sums$region),1:5))stop("All five regions require valid prediction support")
 sums$total_ood_weight<-sums$range_fail_weight+sums$md_fail_weight
 boot<-block_ratio_bootstrap(sums[,c("range_fail_weight","md_fail_weight","total_ood_weight")],
   sums[,c("valid_weight","md_eligible_weight","valid_weight")],sums$region,2000L,42L)
 vals<-c(sum(sums$range_fail_weight)/sum(sums$valid_weight),sum(sums$md_fail_weight)/sum(sums$md_eligible_weight),sum(sums$total_ood_weight)/sum(sums$valid_weight))
 for(j in 1:3)ood[[paste(resp,yr,j)]]<-data.frame(response=resp,year=yr,gate=c("range","MD_conditional_on_range","total_OOD")[j],
   percentage=100*vals[j],lower=100*percentile_ci(boot[,j])[1],upper=100*percentile_ci(boot[,j])[2],
   denominator_weight=sum(sums[[c("valid_weight","md_eligible_weight","valid_weight")[j]]]),B=2000L,blocks=nrow(sums))
 write.csv(sums,file.path(out,paste0(resp,"_",yr,"_OOD_block_sums.csv")),row.names=FALSE)
 # Overlap-weighted averaging supports 30 -> 1000 m, a non-integer ratio.
 ag<-prepare_mean_1km(z[[1]],domain,template,.8)
 one<-ag$mean;vf<-ag$valid_fraction
 d1<-file.path(out,"prepared_1km",resp);dir.create(d1,recursive=TRUE,showWarnings=FALSE)
 path1<-file.path(d1,paste0(resp,"_",yr,"_1km.tif"));terra::writeRaster(one,path1)
 support_dir<-file.path(d1,"support");dir.create(support_dir,showWarnings=FALSE)
 terra::writeRaster(vf,file.path(support_dir,paste0("support_",resp,"_",yr,"_1km.tif")))
 manifest[[i]]<-data.frame(response=resp,year=yr,prediction=normalizePath(native),prepared_1km=normalizePath(path1))
}
write.csv(do.call(rbind,ood),file.path(out,"OOD_annual_percentages_CI.csv"),row.names=FALSE)
write.csv(do.call(rbind,manifest),file.path(out,"masked_annual_manifest.csv"),row.names=FALSE)
write.csv(m,file.path(out,"input_manifest_used.csv"),row.names=FALSE)
terra::writeVector(regions,file.path(out,"regions_used.gpkg"),overwrite=FALSE)
vector_dir<-file.path(out,"prepared_1km","NE_Mountain_Output");dir.create(vector_dir,recursive=TRUE,showWarnings=FALSE)
terra::writeVector(regions,file.path(vector_dir,"NE_Mountain_Regions_All.shp"),overwrite=FALSE)
capture.output(sessionInfo(),file=file.path(out,"sessionInfo.txt"))
writeLines("Completed OOD screening and common-grid preparation; native 30-m outputs feed regional analysis.",file.path(out,"RUN_COMPLETE.txt"))
