log_session <- function(slug) {
  if (!exists("OUTPUT_LOG_DIR")) source(here::here("R", "paths.R"), local = FALSE)
  ts  <- format(Sys.time(), "%Y%m%dT%H%M%S")
  out <- file.path(OUTPUT_LOG_DIR, sprintf("%s_session_%s.txt", slug, ts))
  con <- file(out, open = "wt")
  on.exit(close(con), add = TRUE)
  cat(sprintf("# %s — sessionInfo — %s\n\n",
              slug, format(Sys.time(), tz = "UTC", usetz = TRUE)),
      file = con)
  capture.output(sessionInfo(), file = con, append = TRUE)
  invisible(out)
}

require_file <- function(path, hint = "") {
  if (!file.exists(path)) {
    stop("Required input not found:\n  ", path,
         if (nzchar(hint)) paste0("\n", hint) else "",
         call. = FALSE)
  }
  invisible(path)
}

require_columns <- function(df, cols, what = "data frame") {
  miss <- setdiff(cols, colnames(df))
  if (length(miss)) {
    stop(what, " is missing required column(s): ",
         paste(miss, collapse = ", "), call. = FALSE)
  }
  invisible(df)
}
