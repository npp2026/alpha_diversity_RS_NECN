#!/usr/bin/env Rscript
# v1.6.1: fresh-output, fail-visible annual RF/QM mapping.
# Each task validates and commits mean/quantile/interval-width outputs together.
map_int <- function(name,default) {
  raw<-Sys.getenv(name,as.character(default));x<-suppressWarnings(as.numeric(raw))
  if(length(x)!=1L||!is.finite(x)||x<1||x!=as.integer(x))stop(name," must be a positive integer")
  as.integer(x)
}
map_config <- function() list(
  VRT_DIR=Sys.getenv("VRT_DIR","VRT_Combined/LND_1D"),
  MODEL_DIR=Sys.getenv("MODEL_DIR","final_models_v1"),
  TRAIN_DIR=Sys.getenv("TRAIN_DATA_DIR","train_data"),
  CALIB_DIR=Sys.getenv("CALIB_DIR","qm_calibration"),
  OUTPUT_DIR=Sys.getenv("PREDICTION_OUTPUT_DIR","predictions_final_production_fullQM"),
  TEMP_PARENT=Sys.getenv("TERRA_TMPDIR",tempdir()),
  TARGETS=c("Rich_tree","Shannon_wiener"),YEARS=2001:2020,
  QUANTILES=c(.05,.5,.95),WORKERS=map_int("N_WORKERS",1L),THREADS=map_int("N_THREADS",1L))
load_qm_corrector <- function(target,calib_dir) {
  path<-file.path(calib_dir,paste0("qm_",target,".rds"))
  if(!file.exists(path))stop("Missing QM calibration: ",path)
  q<-readRDS(path)
  if(!identical(q$target,target)||length(q$x_knots)<2L||length(q$x_knots)!=length(q$y_knots)||
     any(!is.finite(q$x_knots))||any(!is.finite(q$y_knots))||any(diff(q$x_knots)<=0)||
     any(diff(q$y_knots)<0)||!is.finite(q$ymin)||!is.finite(q$ymax)||q$ymin>=q$ymax)
    stop("Invalid QM calibration for ",target)
  fn<-stats::approxfun(q$x_knots,q$y_knots,rule=2,ties="ordered")
  function(v){p<-pmin(q$ymax,pmax(q$ymin,fn(v)));p[!is.finite(v)]<-NA_real_;p}
}
make_predict_mean <- function(num_threads,qm) {
  force(num_threads);force(qm)
  function(model,data,...)qm(predict(model,data=data,type="response",num.threads=num_threads)$predictions)
}
make_predict_quantiles <- function(quantiles,num_threads,qm) {
  force(quantiles);force(num_threads);force(qm)
  function(model,data,...) {
    p<-predict(model,data=data,type="quantiles",quantiles=quantiles,num.threads=num_threads)$predictions
    p<-matrix(qm(as.numeric(p)),nrow=nrow(data),ncol=length(quantiles))
    colnames(p)<-paste0("q",format(100*quantiles,trim=TRUE));p
  }
}
map_predictor_stack <- function(vrt,model) {
  x<-terra::rast(vrt);features<-model$forest$independent.variable.names
  if(!length(features)||anyDuplicated(names(x))||!all(features%in%names(x)))
    stop("VRT predictor bands do not match fitted model: ",vrt)
  # Missing values in UNUSED VRT bands must not remove otherwise valid pixels.
  x<-x[[features]]
  terra::ifel(is.finite(x) & x != -9999,x,NA)
}
process_task <- function(task,config,temp_root) {
  target<-task$target;year<-task$year
  worker_tmp<-file.path(temp_root,paste0(target,"_",year));dir.create(worker_tmp)
  on.exit(unlink(worker_tmp,recursive=TRUE),add=TRUE)
  terra::terraOptions(tempdir=worker_tmp,memfrac=max(.02,min(.2,.6/config$WORKERS)),progress=0)
  tryCatch({
    model_path<-file.path(config$MODEL_DIR,paste0(target,"_Final_Model.rds"))
    model<-readRDS(model_path)
    # ranger stores quantile support in random.node.values, not model$quantreg.
    if(is.null(model$random.node.values))stop("Model lacks quantile-forest support: ",target)
    qm<-load_qm_corrector(target,config$CALIB_DIR)
    vrt<-file.path(config$VRT_DIR,sprintf("Combined_%d.vrt",year))
    x<-map_predictor_stack(vrt,model)
    dest<-file.path(config$OUTPUT_DIR,target);dir.create(dest,showWarnings=FALSE)
    final<-file.path(dest,c(sprintf("%s_%d.tif",target,year),
      sprintf("%s_%d_quantiles.tif",target,year),sprintf("%s_%d_uncertainty.tif",target,year)))
    if(any(file.exists(final)))stop("Existing task output; use an empty output directory")
    stage<-file.path(dest,paste0(".partial_",Sys.getpid(),"_",basename(final)))
    on.exit(unlink(stage),add=TRUE)
    wopt<-list(datatype="FLT4S",gdal=c("COMPRESS=LZW","BIGTIFF=IF_SAFER","TILED=YES"))
    mean<-terra::predict(x,model,fun=make_predict_mean(config$THREADS,qm),na.rm=TRUE,
                       cores=1,filename=stage[1],wopt=wopt)
    quant<-terra::predict(x,model,fun=make_predict_quantiles(config$QUANTILES,config$THREADS,qm),
                        na.rm=TRUE,cores=1,filename=stage[2],wopt=wopt)
    if(terra::nlyr(mean)!=1L||terra::nlyr(quant)!=3L)stop("Unexpected output layer count")
    width<-quant[[3]]-quant[[1]];names(width)<-paste0(target,"_q95_minus_q05")
    terra::writeRaster(width,stage[3],wopt=wopt)
    for(j in 1:3)if(!terra::compareGeom(terra::rast(stage[j]),x,stopOnError=FALSE))stop("Output geometry mismatch")
    rm(mean,quant,width,x);gc(verbose=FALSE)
    committed<-file.rename(stage,final)
    if(!all(committed)){unlink(final[committed]);stop("Failed to commit complete output set")}
    list(target=target,year=year,success=TRUE,error="",prediction=normalizePath(final[1]),
         predictors=normalizePath(vrt),model_rds=normalizePath(model_path),
         train_rds=normalizePath(file.path(config$TRAIN_DIR,paste0(target,"_train.rds"))))
  },error=function(e)list(target=target,year=year,success=FALSE,error=conditionMessage(e)))
}
run_prediction <- function(config=map_config()) {
  for(p in c("terra","ranger"))if(!requireNamespace(p,quietly=TRUE))stop("Install ",p)
  # Register ranger's predict method in parent and forked processes.
  loadNamespace("ranger")
  if(dir.exists(config$OUTPUT_DIR)&&length(list.files(config$OUTPUT_DIR,all.files=TRUE,no..=TRUE)))
    stop("Use an empty prediction output directory; existing rasters are never accepted by file size")
  if(!dir.exists(config$TEMP_PARENT))dir.create(config$TEMP_PARENT,recursive=TRUE)
  tmp<-tempfile("FEM_mapping_",tmpdir=config$TEMP_PARENT);dir.create(tmp)
  on.exit(unlink(tmp,recursive=TRUE),add=TRUE)
  required<-c(file.path(config$VRT_DIR,sprintf("Combined_%d.vrt",config$YEARS)),
    file.path(config$MODEL_DIR,paste0(config$TARGETS,"_Final_Model.rds")),
    file.path(config$TRAIN_DIR,paste0(config$TARGETS,"_train.rds")),
    file.path(config$CALIB_DIR,paste0("qm_",config$TARGETS,".rds")))
  if(any(!file.exists(required)))stop("Missing mapping inputs: ",paste(required[!file.exists(required)],collapse=", "))
  dir.create(config$OUTPUT_DIR,recursive=TRUE,showWarnings=FALSE)
  tasks<-unlist(lapply(config$TARGETS,function(t)lapply(config$YEARS,function(y)list(target=t,year=y))),recursive=FALSE)
  workers<-min(config$WORKERS,length(tasks))
  if(.Platform$OS.type=="windows"&&workers>1L){message("Windows: using serial tasks; N_THREADS still controls ranger");workers<-1L}
  results<-if(workers==1L)lapply(tasks,process_task,config=config,temp_root=tmp) else
    parallel::mclapply(tasks,process_task,config=config,temp_root=tmp,mc.cores=workers,mc.preschedule=FALSE)
  # Worker crashes, NULL results and try-errors are failures, never silent success.
  good<-vapply(results,function(x)is.list(x)&&isTRUE(x$success),logical(1))
  status<-do.call(rbind,lapply(seq_along(tasks),function(i)data.frame(response=tasks[[i]]$target,
    year=tasks[[i]]$year,success=good[i],error=if(good[i])""else if(is.list(results[[i]]))results[[i]]$error else "Worker failed")))
  write.csv(status,file.path(config$OUTPUT_DIR,"mapping_task_status.csv"),row.names=FALSE)
  if(!all(good))stop(sum(!good)," mapping tasks failed; inspect mapping_task_status.csv")
  manifest<-do.call(rbind,lapply(results,function(x)data.frame(response=x$target,year=x$year,
    prediction=x$prediction,predictors=x$predictors,train_rds=x$train_rds,model_rds=x$model_rds)))
  write.csv(manifest,file.path(config$OUTPUT_DIR,"annual_manifest.csv"),row.names=FALSE)
  write.csv(data.frame(path=normalizePath(required),md5=unname(tools::md5sum(required))),
            file.path(config$OUTPUT_DIR,"mapping_input_hashes.csv"),row.names=FALSE)
  capture.output(sessionInfo(),file=file.path(config$OUTPUT_DIR,"sessionInfo.txt"))
  writeLines(c("All requested mapping tasks completed: v1.6.1",
    "Uncertainty output is QM-transformed QRF q95-q05 width; it is not a calibrated confidence interval.",
    "Input hashes cover VRT definitions, models, training tables and QM, not externally referenced VRT source pixels."),
    file.path(config$OUTPUT_DIR,"RUN_COMPLETE.txt"))
  invisible(results)
}
if(sys.nframe()==0L)run_prediction()
