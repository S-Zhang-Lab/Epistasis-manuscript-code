log_session <- function(slug) {
  source(here::here("R", "paths.R"), local = FALSE)
  out <- log_path(sprintf("%s_session.txt", slug))
  con <- file(out, open = "wt")
  on.exit(close(con), add = TRUE)
  cat(sprintf("# %s -- %s\n\n", slug, format(Sys.time(), tz = "UTC", usetz = TRUE)),
      file = con)
  capture.output(sessionInfo(), file = con, append = TRUE)
  invisible(out)
}

require_cols <- function(df, cols, where = "input") {
  missing <- setdiff(cols, names(df))
  if (length(missing) > 0L) {
    stop(sprintf("%s missing required columns: %s",
                 where, paste(missing, collapse = ", ")),
         call. = FALSE)
  }
  invisible(df)
}

sym_lim <- function(vals, q = 0.97) {
  stats::quantile(abs(vals), q, na.rm = TRUE)
}

to_mouse_token <- function(g) {
  g <- trimws(g)
  if (g == "NT") return("NT")
  paste0(toupper(substr(g, 1L, 1L)),
         tolower(substr(g, 2L, nchar(g))))
}

convert_gene_pair <- function(pair) {
  parts <- strsplit(pair, "*", fixed = TRUE)[[1L]]
  paste(vapply(parts, to_mouse_token, character(1L)), collapse = "*")
}

convert_sgrna_id <- function(sgrna) {
  parts <- strsplit(sgrna, "*", fixed = TRUE)[[1L]]
  out <- vapply(parts, function(part) {
    idx <- regexpr("_", part, fixed = TRUE)
    if (idx == -1L) {
      to_mouse_token(part)
    } else {
      paste0(to_mouse_token(substr(part, 1L, idx - 1L)),
             substr(part, idx, nchar(part)))
    }
  }, character(1L))
  paste(out, collapse = "*")
}

pair_type <- function(gene_pair) {
  parts <- strsplit(gene_pair, "*", fixed = TRUE)[[1L]]
  n_nt <- sum(parts == "NT")
  if (n_nt == 2L) "NT_NT" else if (n_nt == 1L) "single_NT" else "dual"
}
