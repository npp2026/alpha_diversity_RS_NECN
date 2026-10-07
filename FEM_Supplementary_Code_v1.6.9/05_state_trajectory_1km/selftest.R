#!/usr/bin/env Rscript
# Fast internal self-test for the state-trajectory statistical backend. No project data needed.
options(stringsAsFactors=FALSE)
if(getRversion()<'4.1.0')stop('R >= 4.1.0 is required.',call.=FALSE)
get_script_dir<-function(){a<-commandArgs(FALSE);f<-a[grep('^--file=',a)];if(length(f))dirname(normalizePath(sub('^--file=','',f[1]),winslash='/',mustWork=FALSE))else getwd()}
SCRIPT_DIR<-get_script_dir();source(file.path(SCRIPT_DIR,'v12_posthoc_utils.R'));source(file.path(SCRIPT_DIR,'v9_trend_package_utils.R'))
set.seed(20260628);n<-2400L
d<-data.frame(cell=seq_len(n),x=runif(n,-2e5,2e5),y=runif(n,-2e5,2e5),roman=sample(c('I','II','III','IV','V','outside'),n,TRUE),area_km2=runif(n,20,25),stringsAsFactors=FALSE)
d$rich_class6<-sample(c(V12_CLASS6_LEVELS,NA_character_),n,TRUE,prob=c(.07,.22,.01,.14,.54,.01,.01));d$shannon_class6<-sample(c(V12_CLASS6_LEVELS,NA_character_),n,TRUE,prob=c(.02,.09,.01,.14,.70,.01,.03))
reference_boot<-function(d,class_col,block_size_km,support,W,blocks){bid<-v12_block_id(d,block_size_km);bidx<-match(bid,blocks);cls<-as.character(d[[class_col]]);zv<-v12_zone_value(d$roman);out<-array(NA_real_,dim=c(nrow(W),length(V12_REG_ORDER),length(V12_CLASS6_LEVELS)),dimnames=list(NULL,V12_REG_ORDER,V12_CLASS6_LEVELS));for(i in seq_len(nrow(W))){rw<-W[i,bidx];for(z in V12_REG_ORDER){zm<-if(z=='Overall')rep(TRUE,nrow(d))else zv==z;keep<-support&zm&cls%in%V12_CLASS6_LEVELS;a<-vapply(V12_CLASS6_LEVELS,function(cl)sum(d$area_km2[keep&cls==cl]*rw[keep&cls==cl]),numeric(1));if(sum(a)>0)out[i,z,]<-a/sum(a)}};out}
rv<-as.character(d$rich_class6)%in%V12_CLASS6_LEVELS;sv<-as.character(d$shannon_class6)%in%V12_CLASS6_LEVELS;layout<-v12_prepare_block_layout(d,50,25)
for(mode in c('common','target_specific'))for(universe in c('analysis_support','all_rows')){
 common<-rv&sv;rs<-if(mode=='common')common else rv;ss<-if(mode=='common')common else sv
 ro<-v12_area_cube(d,'rich_class6',50,rs,25,layout=layout);so<-v12_area_cube(d,'shannon_class6',50,ss,25,layout=layout)
 active<-if(universe=='analysis_support')(v12_cube_block_area(ro$cube)>0|v12_cube_block_area(so$cube)>0)else rep(TRUE,length(ro$blocks))
 rc<-v12_select_cube_blocks(ro$cube,active);sc<-v12_select_cube_blocks(so$cube,active);blocks<-ro$blocks[active]
 W<-v12_bootstrap_weights(length(blocks),60,20260628);nr<-v12_cube_boot_props(rc,W);ns<-v12_cube_boot_props(sc,W)
 or<-reference_boot(d,'rich_class6',50,rs,W,blocks);os<-reference_boot(d,'shannon_class6',50,ss,W,blocks)
 err<-max(abs((ns-nr)-(os-or)),na.rm=TRUE);if(!is.finite(err)||err>1e-12)stop('Bootstrap equivalence failed: ',mode,'/',universe,' err=',err,call.=FALSE)
}
obs<-v12_classify6(c(.8,.7999,.2,NA,1,0),c('T+','T0','T-','T+','T-','T0'),.8);if(!identical(obs,c('H+','L0','L-',NA,'H-','L0')))stop('Classifier threshold-boundary test failed.',call.=FALSE)
ta<-v12_target_transition(d,25);if(abs(sum(d$area_km2[rv&sv&d$roman%in%c("I","II","III","IV","V")])-sum(ta$transition$area_km2[ta$transition$zone=='Overall']))>1e-8)stop('Transition area conservation failed.',call.=FALSE)
# Restricted parser must reject arbitrary calls.
if(.v12_safe_expr(parse(text='system("echo bad")')[[1]]))stop('Safe expression parser accepted an arbitrary function call.',call.=FALSE)
if(!.v12_safe_expr(parse(text='c(50,75,100)')[[1]]))stop('Safe expression parser rejected a supported vector.',call.=FALSE)
# Cached combined trend statistics must reproduce separate fallback calculations.
years<-2005:2020;y<-cumsum(rnorm(length(years)))
s1<-v9_sen_slope(years,y,prefer_package=FALSE,return_source=TRUE);r1<-v9_raw_mk_p(years,y);h1<-v9_hamed_rao_mk(years,y,prefer_package=FALSE)
st<-v9_trend_stats(years,y,prefer_package=FALSE,prefer_sen_package=FALSE)
if(max(abs(c(s1$value-as.numeric(st$sen$value),as.numeric(r1['p'])-as.numeric(st$raw['p']),as.numeric(h1['p'])-as.numeric(st$hr['p']))),na.rm=TRUE)>1e-14)stop('Combined trend-stat cache changed estimates.',call.=FALSE)
cat('[selftest] PASS: support-aware bootstrap, legacy equivalence, safe parser, classifier, transition, and cached trend statistics.\n')
