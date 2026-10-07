# QRF input/cache/output contracts; no fitted scientific parameters are changed.
fem_qrf_numeric <- function(x) {
  y<-if(is.numeric(x))x else suppressWarnings(as.numeric(x))
  y[!is.finite(y)|y == -9999]<-NA_real_;y
}
fem_qrf_new_output <- function(path) {
  if(dir.exists(path)&&length(list.files(path,all.files=TRUE,no..=TRUE)))
    stop("Use a fresh QRF output directory; existing results cannot be mixed: ",path)
  dir.create(path,recursive=TRUE,showWarnings=FALSE);invisible(path)
}
fem_qrf_threads <- function(default=8L) {
  x<-suppressWarnings(as.numeric(Sys.getenv("FEM_QRF_THREADS",as.character(default))))
  if(length(x)!=1L||!is.finite(x)||x<1||x!=as.integer(x))stop("FEM_QRF_THREADS must be a positive integer")
  as.integer(x)
}
fem_qrf_fingerprint <- function(paths,mask_zero_as_na=TRUE) {
  if(!length(paths)||any(!file.exists(paths)))stop("Missing QRF cache source")
  info<-file.info(paths)
  rows<-data.frame(path=normalizePath(paths,winslash="/",mustWork=TRUE),size=as.numeric(info$size),
    mtime=as.character(info$mtime),md5=unname(tools::md5sum(paths)),stringsAsFactors=FALSE)
  rbind(rows,data.frame(path="__config_mask_zero_as_na__",size=as.numeric(isTRUE(mask_zero_as_na)),mtime="config",md5="config"))
}
fem_qrf_cache_matches <- function(meta_file,source_files,mask_zero_as_na=TRUE) {
  if(is.null(source_files)||is.null(meta_file)||!file.exists(meta_file))return(FALSE)
  old<-try(utils::read.csv(meta_file,stringsAsFactors=FALSE),silent=TRUE)
  if(inherits(old,"try-error"))return(FALSE)
  cur<-fem_qrf_fingerprint(source_files,mask_zero_as_na)
  if(!all(names(cur)%in%names(old))||nrow(old)!=nrow(cur))return(FALSE)
  old<-old[order(old$path),names(cur)];cur<-cur[order(cur$path),];rownames(old)<-rownames(cur)<-NULL
  isTRUE(all.equal(old,cur,check.attributes=FALSE))
}
fem_qrf_completed_run <- function(run_path,input_paths,mode) {
  if(!file.exists(run_path))stop("Completed QRF run provenance required: ",run_path)
  z<-readRDS(run_path)
  if(!identical(z$qrf_parameter_mode,mode)||is.null(z$input_md5))stop("QRF source profile/provenance mismatch; rerun predict_age_sensitivity.R from this package")
  current<-setNames(unname(tools::md5sum(input_paths)),basename(input_paths))
  if(anyNA(current)||!identical(current,z$input_md5[names(current)]))stop("QRF inputs changed after prediction; refusing stale AOA/validation")
  z
}
