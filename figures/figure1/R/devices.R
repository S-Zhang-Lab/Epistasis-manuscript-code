PDF_ENCODING <- "WinAnsi"

pdf_device <- function(filename, width, height, ...) {
  dots <- list(...)
  keep <- names(dots) %in% names(formals(grDevices::pdf))
  do.call(grDevices::pdf,
          c(list(file = filename, width = width, height = height,
                 useDingbats = FALSE, encoding = PDF_ENCODING),
            dots[keep]))
}

check_pdf_text <- function(..., .context = "panel text") {
  x <- unlist(list(...), use.names = FALSE)
  x <- x[!is.na(x)]
  if (!length(x)) return(invisible(TRUE))
  conv <- iconv(x, from = "UTF-8", to = "CP1252")
  bad  <- x[is.na(conv)]
  if (length(bad)) {
    offending <- unique(unlist(lapply(bad, function(s) {
      ch <- strsplit(s, "")[[1]]
      ch[is.na(iconv(ch, "UTF-8", "CP1252"))]
    })))
    stop(sprintf(
      paste0("[devices] %s contains %d character(s) the PDF encoding cannot write: %s\n",
             "  The base pdf() device maps text through %s. Write Greek and math\n",
             "  relations as plotmath expressions, or spell them in ASCII.\n",
             "  Offending string(s): %s"),
      .context, length(offending),
      paste(sprintf("U+%04X (%s)", utf8ToInt(paste(offending, collapse = ""))[seq_along(offending)],
                    offending), collapse = ", "),
      PDF_ENCODING,
      paste(shQuote(utils::head(bad, 3)), collapse = ", ")),
      call. = FALSE)
  }
  invisible(TRUE)
}
