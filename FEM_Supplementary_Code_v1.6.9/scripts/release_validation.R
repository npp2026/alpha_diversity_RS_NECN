# Delegate to the actual standard-library Python integrity/syntax validator.
file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (!length(file_arg)) stop("Run this script with Rscript.", call. = FALSE)
here <- dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE))
python <- Sys.getenv("FEM_PYTHON", unset = "")
if (!nzchar(python)) {
  python_candidates <- Sys.which(c("python3", "python"))
  python_candidates <- python_candidates[nzchar(python_candidates)]
  if (!length(python_candidates)) stop("Install Python >= 3.9 or set FEM_PYTHON.", call. = FALSE)
  python <- unname(python_candidates[[1L]])
}
code <- system2(python, c(shQuote(file.path(here, "validate_release.py")), shQuote(commandArgs(trailingOnly = TRUE))))
quit(save = "no", status = as.integer(code))
