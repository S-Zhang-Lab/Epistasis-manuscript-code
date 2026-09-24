source(here::here("R", "paths.R"))

REQUIRED_RAW <- c(
  "Cx3cr1_raw_counts.csv",
  "Mrc1_raw_counts.csv"
)

OPTIONAL_RAW <- c(
  "Mrc1_Pten_x_Cdh1_Image/"
)

REQUIRED_DERIVED <- c(
  "geneLFC_Cx3cr1.csv",
  "geneLFC_Mrc1.csv"
)

REQUIRED_REFERENCE <- character()

status_table <- function(root, files) {
  paths <- file.path(root, files)
  data.frame(
    file    = files,
    path    = paths,
    present = file.exists(paths) | dir.exists(paths),
    stringsAsFactors = FALSE
  )
}

report_missing <- function(label, tbl) {
  missing <- tbl[!tbl$present, , drop = FALSE]
  if (nrow(missing) == 0L) {
    message(sprintf("[00_check_input_data] OK %s: %d file(s)", label, nrow(tbl)))
    return(invisible(character()))
  }
  message(sprintf("[00_check_input_data] MISSING %s:", label))
  for (p in missing$path) message("  - ", p)
  missing$path
}

verify_manifest <- function() {
  manifest <- raw("MANIFEST.sha256")
  if (!file.exists(manifest)) {
    message("[00_check_input_data] No data/raw/MANIFEST.sha256 found; checksum verification skipped.")
    return(invisible(TRUE))
  }
  if (!requireNamespace("digest", quietly = TRUE)) {
    stop("[00_check_input_data] Package 'digest' is required to verify MANIFEST.sha256.", call. = FALSE)
  }

  expected <- read.table(manifest, header = FALSE, col.names = c("sha", "file"),
                         stringsAsFactors = FALSE)
  bad <- character()
  missing <- character()
  for (i in seq_len(nrow(expected))) {
    f <- raw(expected$file[i])
    if (!file.exists(f)) {
      missing <- c(missing, expected$file[i])
      next
    }
    got <- digest::digest(file = f, algo = "sha256")
    if (!identical(tolower(got), tolower(expected$sha[i]))) {
      bad <- c(bad, expected$file[i])
    }
  }

  if (length(bad)) {
    stop("[00_check_input_data] MANIFEST.sha256 verification failed.\nSHA-256 mismatch:\n  ",
         paste(bad, collapse = "\n  "), call. = FALSE)
  }

  if (length(missing)) {
    message(sprintf(
      "[00_check_input_data] MANIFEST.sha256: %d verified, %d not present locally (skipped — not redistributed)",
      nrow(expected) - length(missing), length(missing)))
  } else {
    message("[00_check_input_data] OK MANIFEST.sha256")
  }
  invisible(TRUE)
}

validate_raw_counts <- function() {
  fatal <- character()
  for (fn in REQUIRED_RAW) {
    f <- raw(fn)
    if (!file.exists(f)) next
    df <- utils::read.csv(f, check.names = FALSE, stringsAsFactors = FALSE)
    if (ncol(df) < 3L) { fatal <- c(fatal, sprintf("%s: <3 columns", fn)); next }
    id <- df[[1L]]
    M  <- suppressWarnings(as.matrix(vapply(df[, -(1:2), drop = FALSE],
              function(x) as.numeric(as.character(x)), numeric(nrow(df)))))

    n_bad <- sum(is.na(M))
    if (n_bad) fatal <- c(fatal, sprintf("%s: %d non-numeric/NA count cell(s)", fn, n_bad))
    dup <- unique(id[duplicated(id)])
    if (length(dup)) fatal <- c(fatal, sprintf("%s: %d duplicate sgRNA id(s), e.g. %s",
              fn, length(dup), paste(utils::head(dup, 3L), collapse = ", ")))
    lib <- colSums(M, na.rm = TRUE)
    if (any(lib <= 0)) fatal <- c(fatal, sprintf("%s: empty library column(s): %s",
              fn, paste(names(lib)[lib <= 0], collapse = ", ")))

    neg <- which(M < 0, arr.ind = TRUE)
    if (nrow(neg)) {
      ex <- apply(utils::head(neg, 4L), 1L, function(r)
        sprintf("%s/%s=%s", id[r[[1]]], colnames(M)[r[[2]]], M[r[[1]], r[[2]]]))
      fatal <- c(fatal, sprintf("%s: %d negative count(s) -- invalid as sequencing counts (%s)",
              fn, nrow(neg), paste(ex, collapse = ", ")))
    }

    nint <- sum(M != round(M), na.rm = TRUE)
    if (nint) message(sprintf("[00_check_input_data] WARNING %s: %d non-integer count value(s).", fn, nint))

    for (s in colnames(M)) {
      x <- M[, s]; x[is.na(x) | x < 0] <- 0; tot <- sum(x)
      if (tot <= 0) next
      p <- x[x > 0] / tot; eff <- exp(-sum(p * log(p)))
      if (eff < 20)
        message(sprintf("[00_check_input_data] NOTE %s/%s: low complexity -- effGuides=%.0f, zeros=%.1f%%, top guide=%.1f%% of reads",
                fn, s, eff, mean(x == 0) * 100, max(x) / tot * 100))
    }
  }
  if (length(fatal))
    stop("[00_check_input_data] Raw-count validation failed:\n  ",
         paste(fatal, collapse = "\n  "), call. = FALSE)
  message("[00_check_input_data] OK raw-count value validation (structural checks passed).")
  invisible(TRUE)
}

raw_status       <- status_table(DATA_RAW,     REQUIRED_RAW)
derived_status   <- status_table(DATA_DERIVED, REQUIRED_DERIVED)
optional_status  <- status_table(DATA_RAW,     OPTIONAL_RAW)

missing <- c(
  report_missing("data/raw committed sgRNA count inputs", raw_status),
  report_missing("data/derived curated inputs",           derived_status)
)

for (i in seq_len(nrow(optional_status))) {
  if (optional_status$present[i]) {
    message("[00_check_input_data] OK (optional) ", optional_status$file[i])
  } else {
    message("[00_check_input_data] NOTE (optional, not redistributed) ",
            optional_status$file[i], " absent — the Panel I photographic montage ",
            "cannot be rebuilt, but every numerical result, including the Panel I ",
            "burden statistics, reproduces without it.")
  }
}

verify_manifest()

if (!length(missing)) validate_raw_counts()

if (length(missing)) {
  stop("[00_check_input_data] Required Figure 5 inputs are missing. ",
       "These are committed to the repository, so a complete clone should already ",
       "contain them — check for a partial checkout or deleted files.",
       call. = FALSE)
}

message("[00_check_input_data] OK — all required Figure 5 inputs are present.")
