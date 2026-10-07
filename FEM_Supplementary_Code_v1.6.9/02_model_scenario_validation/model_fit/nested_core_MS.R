# MS 2.4, SI S3.1-S3.2. No file IO or model fitting on source().
ms_scenarios <- list(S0_Env="Env", S1a_Prod="Prod", S1b_Het="Het", S1c_Temp="Temp",
  S2a_Prod_Het=c("Prod","Het"), S2b_Prod_Temp=c("Prod","Temp"),
  S2c_Het_Temp=c("Het","Temp"), S3_RS_Full=c("Prod","Het","Temp"),
  S4_Integrated=c("Env","Prod","Het","Temp"))
ms_folds <- function(id,k=5L,seed=42L) {
  u <- unique(as.character(id)); if(length(u)<k)stop("Too few independent plots for folds")
  set.seed(seed); f <- sample(rep(seq_len(k),length.out=length(u)))
  f[match(as.character(id),u)]
}
ms_lhs <- function(n,p,seed) {
  set.seed(seed)
  vapply(seq_len(p),function(j)(sample.int(n)-runif(n))/n,numeric(n))
}
ms_candidates <- function(algorithm,seed,n=30L) {
  u <- ms_lhs(n,if(algorithm=="RF")4L else 7L,seed)
  if(algorithm=="RF")return(data.frame(num.trees=round(300+1700*u[,1]),
    mtry_fraction=.25+.55*u[,2],min.node.size=round(1+14*u[,3]),sample.fraction=.5+.4*u[,4]))
  data.frame(eta=exp(log(.01)+log(12)*u[,1]),max_depth=round(3+5*u[,2]),
    min_child_weight=1+7*u[,3],subsample=.6+.35*u[,4],colsample_bytree=.5+.45*u[,5],
    lambda=exp(log(.1)+log(30)*u[,6]),alpha=u[,7])
}
ms_screen <- function(x,y,seed) {
  # Only the caller's training rows reach Spearman, VIF and VSURF.
  x <- as.data.frame(x,check.names=FALSE)
  good <- vapply(x,function(z)is.numeric(z)&&all(is.finite(z))&&stats::sd(z)>0,logical(1))
  x <- x[,good,drop=FALSE]; if(!ncol(x))stop("No variable predictors in training fold")
  rho <- vapply(x,function(z)abs(stats::cor(z,y,method="spearman")),numeric(1))
  rho[!is.finite(rho)] <- 0
  base <- gsub("_(w[0-9]+|[0-9]+m|[0-9]+km)(?=_|$)","",names(x),perl=TRUE)
  take <- vapply(split(seq_along(base),base),function(ix)ix[order(-rho[ix],names(x)[ix])[1]],integer(1))
  x <- x[,sort(take),drop=FALSE]; scale_selected <- names(x)
  repeat {
    if(ncol(x)<=1L)break
    vif <- vapply(seq_len(ncol(x)),function(j){
      z <- stats::lm.fit(cbind(1,as.matrix(x[,-j,drop=FALSE])),x[[j]])
      rss <- sum(z$residuals^2); tss <- sum((x[[j]]-mean(x[[j]]))^2)
      if(rss <= .Machine$double.eps*tss) Inf else tss/rss
    },numeric(1))
    if(max(vif)<=10)break
    x <- x[,-which.max(vif),drop=FALSE]
  }
  vif_selected <- names(x)
  set.seed(seed)
  selection_stage <- "single_available_feature"
  if(ncol(x)==1L) sel <- names(x) else {
    v <- VSURF::VSURF(x=x,y=y,parallel=FALSE,verbose=FALSE)
    if(length(v$varselect.pred)) {
      sel <- names(x)[v$varselect.pred]; selection_stage <- "VSURF_prediction"
    } else if(length(v$varselect.interp)) {
      sel <- names(x)[v$varselect.interp]; selection_stage <- "VSURF_interpretation_prediction_unavailable"
      warning("VSURF prediction step unavailable; retaining interpretation-step features. See fold audit.",call.=FALSE)
    } else stop("VSURF selected no usable predictors in this training fold")
  }
  list(features=sel,scale_selected=scale_selected,vif_selected=vif_selected,
       training_n=length(y),seed=seed,selection_stage=selection_stage)
}
ms_metric <- function(y,p) {
  stopifnot(length(y)==length(p),all(is.finite(y)),all(is.finite(p)))
  c(R2=if(sum((y-mean(y))^2)>0)1-sum((y-p)^2)/sum((y-mean(y))^2)else NA_real_,
    RMSE=sqrt(mean((y-p)^2)),MAE=mean(abs(y-p)))
}
# Source relative to this file, including when loaded by the validation runner.
.fem_source_files <- unlist(lapply(sys.frames(),function(e)e$ofile),use.names=FALSE)
source(file.path(dirname(normalizePath(tail(.fem_source_files,1L))),"..","..","R","qm_core.R"))
rm(.fem_source_files)
ms_qm <- function(oob,y,response,seed=42L) {
  q <- fem_fit_qm(oob,y,response,seed)
  list(x=q$x_knots,y=q$y_knots,upper=q$ymax)
}
ms_apply_qm <- function(q,p) fem_apply_qm(list(x_knots=q$x,y_knots=q$y,ymin=0,ymax=q$upper),p)
ms_rf <- function(x,y,h,seed) {
  ranger::ranger(x=x,y=y,num.trees=as.integer(h$num.trees),
    mtry=max(1L,min(ncol(x),round(h$mtry_fraction*ncol(x)))),
    min.node.size=as.integer(h$min.node.size),sample.fraction=h$sample.fraction,
    importance="none",num.threads=1L,seed=seed)
}
ms_xgb <- function(x,y,h,seed,nrounds=800L,valid=NULL) {
  params<-c(as.list(h),list(objective="reg:squarederror",eval_metric="rmse",nthread=1L,seed=seed))
  d<-xgboost::xgb.DMatrix(as.matrix(x),label=y)
  a<-list(params=params,data=d,nrounds=as.integer(nrounds),verbose=0)
  if(!is.null(valid)) {
    dv<-xgboost::xgb.DMatrix(as.matrix(valid$x),label=valid$y)
    # xgboost renamed watchlist to evals; select an explicitly advertised formal.
    arg<-if("evals"%in%names(formals(xgboost::xgb.train)))"evals" else "watchlist"
    a[[arg]]<-list(valid=dv);a$early_stopping_rounds<-30L
  }
  do.call(xgboost::xgb.train,a)
}
ms_best_round <- function(model) {
  history<-attr(model,"evaluation_log",exact=TRUE)
  if(is.null(history))history<-model$evaluation_log
  log<-as.data.frame(history)
  nm<-grep("valid.*rmse",names(log),value=TRUE)
  if(length(nm)!=1L||!nrow(log))stop("Cannot identify XGBoost validation RMSE history")
  which.min(log[[nm]])  # number of rounds is row index, independent of API iteration indexing
}
ms_tune <- function(x,y,algorithm,inner_screens,inner_fold,seed) {
  candidates<-ms_candidates(algorithm,seed)
  scores<-matrix(NA_real_,nrow(candidates),max(inner_fold)); rounds<-scores
  for(cand in seq_len(nrow(candidates)))for(k in sort(unique(inner_fold))) {
    tr<-which(inner_fold!=k);va<-which(inner_fold==k); fs<-inner_screens[[k]]$features
    h<-candidates[cand,,drop=FALSE]
    if(algorithm=="RF") {
      m<-ms_rf(x[tr,fs,drop=FALSE],y[tr],h,seed+1000*cand+k)
      pr<-predict(m,data=x[va,fs,drop=FALSE])$predictions
    }else {
      m<-ms_xgb(x[tr,fs,drop=FALSE],y[tr],h,seed+1000*cand+k,
                valid=list(x=x[va,fs,drop=FALSE],y=y[va]))
      best<-ms_best_round(m); rounds[cand,k]<-best
      # Refit exactly best rounds: prediction cannot accidentally use post-best trees.
      m<-ms_xgb(x[tr,fs,drop=FALSE],y[tr],h,seed+1000*cand+k,nrounds=best)
      pr<-as.numeric(predict(m,xgboost::xgb.DMatrix(as.matrix(x[va,fs,drop=FALSE]))))
    }
    scores[cand,k]<-ms_metric(y[va],pr)["RMSE"]
  }
  chosen<-which.min(rowMeans(scores)); h<-candidates[chosen,,drop=FALSE]
  list(parameters=h,selected_candidate=chosen,candidates=candidates,
       inner_RMSE=scores,rounds=rounds,
       nrounds=if(algorithm=="XGBoost")max(1L,round(median(rounds[chosen,])))else NA_integer_)
}
ms_outer <- function(train,test,response,features,seed) {
  x<-train[,features,drop=FALSE]; y<-train[[response]]
  inner_fold<-ms_folds(train$plot_id,5L,seed)
  inner_screens<-lapply(1:5,function(k)ms_screen(x[inner_fold!=k,,drop=FALSE],y[inner_fold!=k],seed+k))
  outer_screen<-ms_screen(x,y,seed+10L);fs<-outer_screen$features
  ans<-list(); audits<-list()
  for(algorithm in c("RF","XGBoost")) {
    tuning<-ms_tune(x,y,algorithm,inner_screens,inner_fold,seed+20L)
    if(algorithm=="RF") {
      m<-ms_rf(x[,fs,drop=FALSE],y,tuning$parameters,seed+30L)
      raw<-predict(m,data=test[,fs,drop=FALSE])$predictions
      qm<-ms_qm(m$predictions,y,response,seed+40L)
      calibrated<-ms_apply_qm(qm,raw)
    }else {
      m<-ms_xgb(x[,fs,drop=FALSE],y,tuning$parameters,seed+30L,tuning$nrounds)
      raw<-as.numeric(predict(m,xgboost::xgb.DMatrix(as.matrix(test[,fs,drop=FALSE]))))
      calibrated<-rep(NA_real_,length(raw)); qm<-NULL
    }
    ans[[algorithm]]<-list(raw=raw,calibrated=calibrated)
    audits[[algorithm]]<-list(tuning=tuning,QM=qm)
  }
  list(predictions=ans,audit=list(train_ids=train$plot_id,test_ids=test$plot_id,
    inner_fold=inner_fold,inner_screens=inner_screens,outer_screen=outer_screen,algorithms=audits))
}
