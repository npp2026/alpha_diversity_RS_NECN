#!/usr/bin/env Rscript
# Generate current Table S11: pixel-wise richness x Shannon six-class cross-classification
# and agreement diagnostics on their common valid 1-km support.
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({library(dplyr)})

get_script_dir <- function() {
  a <- commandArgs(FALSE); f <- a[grep('^--file=', a)]
  if (length(f)) dirname(normalizePath(sub('^--file=', '', f[1]), winslash='/', mustWork=FALSE)) else getwd()
}
SCRIPT_DIR <- get_script_dir()
source(file.path(SCRIPT_DIR, 'v12_posthoc_utils.R'))

DATA_DIR <- normalizePath(v12_env_str('DATA_DIR','.'), winslash='/', mustWork=TRUE)
OUT_ROOT <- v12_make_abs_path(v12_env_str('OUT_ROOT_NAME','outputs_RF_6class_fig6_minimal'), DATA_DIR)
GROUP <- v12_env_str('ANALYSIS_GROUP','trend_2005_2020')
CELL_SIZE_M <- v12_env_int('CELL_SIZE_M',1000,1)
CELL_AREA_KM2 <- v12_env_num('CELL_AREA_KM2',(CELL_SIZE_M/1000)^2,0,Inf,FALSE,TRUE)

infile <- file.path(OUT_ROOT, GROUP, 'RF_6class_cell_table.csv')
if (!file.exists(infile)) stop('Missing cell table: ', infile, call.=FALSE)
outdir <- file.path(OUT_ROOT, GROUP, 'TableS11_cross_metric_agreement')
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

d <- read.csv(infile, check.names=FALSE)
d <- v12_validate_cell_table(d, require_rf=FALSE, path=infile)
d <- v12_ensure_area(d, CELL_AREA_KM2)
ta <- v12_target_transition(d, CELL_AREA_KM2)

# Panel A: 6 x 6 matrix, normalized within each reporting unit's common support.
write.csv(ta$transition,
          file.path(outdir,'TableS11A_cross_classification_matrix_long.csv'),
          row.names=FALSE)

# Panel B: exact six-class, H/L-state, trend-state agreement, and area-weighted Cohen kappa.
write.csv(ta$agreement,
          file.path(outdir,'TableS11B_cross_metric_agreement_summary.csv'),
          row.names=FALSE)

# Optional audit table: focal-four coverage by metric/reporting unit.
write.csv(ta$coverage,
          file.path(outdir,'TableS11_support_QC.csv'),
          row.names=FALSE)

v12_msg('[Table S11] written: %s', outdir)
