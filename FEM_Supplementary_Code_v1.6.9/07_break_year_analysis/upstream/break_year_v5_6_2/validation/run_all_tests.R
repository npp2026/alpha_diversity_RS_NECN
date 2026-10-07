args0<-commandArgs(FALSE);farg<-grep("^--file=",args0,value=TRUE);PROJECT_DIR<-if(length(farg))normalizePath(dirname(dirname(sub("^--file=","",farg[[1L]]))),winslash="/",mustWork=TRUE)else normalizePath(getwd(),winslash="/",mustWork=TRUE)
source(file.path(PROJECT_DIR,"validation","unit_tests.R"));run_unit_tests_v56(PROJECT_DIR)
source(file.path(PROJECT_DIR,"validation","smoke_test.R"));run_smoke_test_v56(PROJECT_DIR)
source(file.path(PROJECT_DIR,"validation","parallel_reproducibility.R"));pr<-run_parallel_reproducibility_v56(PROJECT_DIR);stopifnot(all(pr$pass));cat("Parallel reproducibility tests passed.\n")
source(file.path(PROJECT_DIR,"validation","warning_regression.R"));wr<-run_warning_regression_v56(PROJECT_DIR);stopifnot(all(wr$pass))
source(file.path(PROJECT_DIR,"validation","region_metadata_tests.R"));rmtest<-run_region_metadata_tests_v56(PROJECT_DIR);stopifnot(all(rmtest$pass))
cat("All algorithm-only v5.6 tests passed.\n")
args<-commandArgs(trailingOnly=TRUE)
if("--synthetic" %in% args){
  source(file.path(PROJECT_DIR,"R","load_all.R"));load_v56_modules(PROJECT_DIR)
  source(file.path(PROJECT_DIR,"validation","synthetic_e2e.R"))
  z<-run_synthetic_e2e_v56(PROJECT_DIR,reuse_data="--reuse-data" %in% args,fail_on_hard=TRUE)
  cat("Synthetic E2E status: ",z$evaluation$summary$status,"\n",sep="")
}
