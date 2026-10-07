#!/usr/bin/env Rscript
args0<-commandArgs(FALSE);farg<-grep("^--file=",args0,value=TRUE);PROJECT_DIR<-if(length(farg))normalizePath(dirname(sub("^--file=","",farg[[1L]])),winslash="/",mustWork=TRUE)else normalizePath(getwd(),winslash="/",mustWork=TRUE)
source(file.path(PROJECT_DIR,"R","load_all.R"));load_v56_modules(PROJECT_DIR)
source(file.path(PROJECT_DIR,"validation","unit_tests.R"));run_unit_tests_v56(PROJECT_DIR)
source(file.path(PROJECT_DIR,"validation","smoke_test.R"));run_smoke_test_v56(PROJECT_DIR)
cat("Algorithm-only validation complete. Data-dependent validation is performed by 00_run_all.R when enabled.\n")
