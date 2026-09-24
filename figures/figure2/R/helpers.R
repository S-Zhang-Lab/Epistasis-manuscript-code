log_session <- function(slug) {
  if (!exists("OUT_LOGS")) source(here::here("R", "paths.R"), local = FALSE)
  out <- log_path(sprintf("%s_session.txt", slug))
  con <- file(out, open = "wt")
  on.exit(close(con), add = TRUE)
  cat(sprintf("# %s -- %s\n\n", slug, format(Sys.time(), tz = "UTC", usetz = TRUE)),
      file = con)
  capture.output(sessionInfo(), file = con, append = TRUE)
  invisible(out)
}

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) {
    return(normalizePath(dirname(sub("^--file=", "", file_arg[1]))))
  }
  if (!is.null(sys.frames()) && length(sys.frames())) {
    sf <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) NULL)
    if (!is.null(sf)) return(dirname(sf))
  }
  getwd()
}
