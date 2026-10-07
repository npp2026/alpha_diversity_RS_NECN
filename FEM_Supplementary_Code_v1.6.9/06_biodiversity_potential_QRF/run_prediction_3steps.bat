@echo off
setlocal EnableExtensions
REM Portable launcher; configure FEM_POTENTIAL_DATA_DIR before running.
REM Rscript run_potential_workflow.R --help describes all options.
if not defined RSCRIPT set "RSCRIPT=Rscript"
"%RSCRIPT%" "%~dp0run_potential_workflow.R" %*
exit /b %errorlevel%
