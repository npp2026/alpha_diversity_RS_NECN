#!/usr/bin/env Rscript
# No annual replication: export one q95 band at ForestAge=100 for each response.
if (!requireNamespace("terra", quietly=TRUE)) stop("Install terra")
root <- Sys.getenv("FEM_POTENTIAL_DATA_DIR", "")
if (!nzchar(root)) stop("Set FEM_POTENTIAL_DATA_DIR")
mode <- Sys.getenv("FEM_QRF_PARAMETER_MODE", "manuscript")
if (mode != "manuscript") stop("Retuned output is not the MS profile; export explicitly after scientific review.")
run_path<-file.path(root,"quantile_sensitivity_outputs","run_info.rds")
if(!file.exists(run_path))stop("Missing QRF run provenance; regenerate prediction with this package")
.fem_script<-sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1])
source(file.path(dirname(normalizePath(.fem_script)),"..","R","qrf_contracts.R"))
run<-fem_qrf_completed_run(run_path,file.path(root,c("ENV.tif","train4pot.csv","templ_1km.tif")),mode)
if(!identical(run$qrf_parameter_mode,"manuscript")||is.null(run$age100_surface_md5))
  stop("Source QRF run is unverified or retuned; environment flags cannot relabel it")
out <- file.path(root,"static_q95_age100")
if(dir.exists(out)&&length(list.files(out)))stop("Use an empty static potential output directory")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
template <- terra::rast(file.path(root,"templ_1km.tif"))
rows <- list()
for (resp in c("Rich_tree","Shannon_wiener")) {
  f <- file.path(root,"quantile_sensitivity_outputs",paste0("RANGER_QRF_",resp,"_ForestAge_100_quantile_surface.tif"))
  if (!file.exists(f)) stop("Missing prediction: ",f)
  if(!identical(unname(tools::md5sum(f)),unname(run$age100_surface_md5[resp])))stop("QRF surface differs from completed run: ",f)
  pp<-run$qrf_params_by_response[[resp]]
  if(is.null(pp)||pp$qrf_num_trees_final!=1500L||pp$qrf_min_node_size!=10L||
     pp$qrf_mtry!=if(resp=="Rich_tree")4L else 2L)stop("QRF parameters do not match manuscript profile")
  if(resp=="Rich_tree"&&!identical(pp$qrf_always_split_variables,"Forest_age"))stop("Missing forced age predictor")
  if(resp=="Shannon_wiener"&&length(pp$qrf_always_split_variables))stop("Unexpected forced Shannon predictor")
  z <- terra::rast(f)
  if (sum(names(z)=="q95")!=1L) stop("Expected exactly one q95 band: ",f)
  z <- z[["q95"]]
  if (!terra::compareGeom(z,template,stopOnError=FALSE)) stop("QRF/template geometry mismatch")
  if (terra::is.lonlat(z) || any(abs(terra::res(z)-1000)>0.01)) stop("Expected projected 1-km grid")
  dest <- file.path(out,paste0("Q95_",resp,"_age100_1km.tif"))
  z<-terra::ifel(is.finite(z)&z>0&is.finite(template)&template!=0,z,NA)
  terra::writeRaster(z,dest,overwrite=FALSE)
  rows[[resp]] <- data.frame(response=resp,reference_age=100,quantile=.95,
    parameter_profile=mode,source=normalizePath(f),source_md5=unname(tools::md5sum(f)),
    output=normalizePath(dest),output_md5=unname(tools::md5sum(dest)))
}
write.csv(do.call(rbind,rows),file.path(out,"static_potential_manifest.csv"),row.names=FALSE)
