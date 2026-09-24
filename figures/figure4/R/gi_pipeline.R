suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(readr)
})

TOP_DUAL  <- 5L
TOP_NT    <- 18L
LFC_COLS  <- c("WT_vs_Cell", "Depleted_vs_Cell")

load_and_normalize <- function(counts_path) {
  if (!file.exists(counts_path)) {
    stop("Raw count file not found: ", counts_path, call. = FALSE)
  }
  df_raw <- readr::read_csv(counts_path, show_col_types = FALSE,
                            progress = FALSE)
  require_cols(df_raw, c("sgRNA", "gene"), where = counts_path)

  sample_cols <- setdiff(names(df_raw), c("sgRNA", "gene"))
  if (length(sample_cols) == 0L) {
    stop("Raw count file has no sample columns beyond sgRNA/gene.",
         call. = FALSE)
  }

  by_sgrna <- df_raw %>%
    dplyr::group_by(sgRNA) %>%
    dplyr::summarise(dplyr::across(dplyr::all_of(sample_cols), sum),
                     .groups = "drop")
  gene_map <- df_raw %>%
    dplyr::group_by(sgRNA) %>%
    dplyr::summarise(gene = dplyr::first(gene), .groups = "drop") %>%
    tibble::deframe()

  counts <- as.matrix(by_sgrna[, sample_cols, drop = FALSE])
  rownames(counts) <- by_sgrna$sgRNA
  storage.mode(counts) <- "double"

  if (!isTRUE(Sys.getenv("FIG4_ALLOW_BAD_COUNTS") == "1")) {
    bad <- !is.finite(counts) | counts < 0 | counts != round(counts)
    bad[is.na(bad)] <- TRUE
    if (any(bad)) {
      rc <- which(bad, arr.ind = TRUE)
      msg <- utils::head(sprintf("    %s [%s]  = %s",
                                 rownames(counts)[rc[, 1L]],
                                 colnames(counts)[rc[, 2L]],
                                 counts[bad]), 10L)
      stop(sprintf(paste0(
        "Invalid raw counts: %d cell(s) in %s\n%s%s\n",
        "  A count must be a finite, non-negative whole number."),
        sum(bad), basename(counts_path), paste(msg, collapse = "\n"),
        if (sum(bad) > 10L) "\n    ... (truncated)" else ""),
        call. = FALSE)
    }
  }

  lib_sizes <- colSums(counts)
  if (any(lib_sizes == 0)) {
    bad <- names(lib_sizes)[lib_sizes == 0]
    stop("Zero library size in sample(s): ", paste(bad, collapse = ", "),
         call. = FALSE)
  }
  cpm <- sweep(counts, 2L, lib_sizes, "/") * 1e6
  log2cpm <- log2(cpm + 1)

  list(counts = counts, gene_map = gene_map[rownames(counts)],
       log2cpm = log2cpm, sample_cols = sample_cols)
}

classify_samples <- function(sample_cols) {
  list(
    cell     = sample_cols[startsWith(sample_cols, "Cell_")],
    wt       = sample_cols[startsWith(sample_cols, "WT_")],
    depleted = sample_cols[startsWith(sample_cols, "Depleted_")]
  )
}

calc_guide_lfc <- function(log2cpm, sample_cols, gene_map, out_path = NULL) {
  groups <- classify_samples(sample_cols)
  stopifnot(length(groups$cell) > 0L,
            length(groups$wt) > 0L,
            length(groups$depleted) > 0L)

  mean_cell <- rowMeans(log2cpm[, groups$cell,     drop = FALSE])
  mean_wt   <- rowMeans(log2cpm[, groups$wt,       drop = FALSE])
  mean_dep  <- rowMeans(log2cpm[, groups$depleted, drop = FALSE])

  lfc <- tibble::tibble(
    sgRNA_raw            = rownames(log2cpm),
    LFC_WT_vs_Cell       = round(mean_wt  - mean_cell, 4L),
    LFC_Depleted_vs_Cell = round(mean_dep - mean_cell, 4L)
  ) %>%
    dplyr::mutate(
      sgRNA = vapply(sgRNA_raw, convert_sgrna_id, character(1L),
                     USE.NAMES = FALSE),
      gene  = vapply(unname(gene_map[sgRNA_raw]), convert_gene_pair,
                     character(1L), USE.NAMES = FALSE)
    ) %>%
    dplyr::select(sgRNA, gene, LFC_WT_vs_Cell, LFC_Depleted_vs_Cell)

  if (!is.null(out_path)) {
    dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)
    readr::write_csv(lfc, out_path)
    message(sprintf("    Wrote %s (%d rows)", out_path, nrow(lfc)))
  }
  lfc
}

calc_gene_lfc <- function(guide_lfc, out_path = NULL,
                          top_dual = TOP_DUAL, top_nt = TOP_NT) {
  require_cols(guide_lfc, c("sgRNA", "gene",
                            "LFC_WT_vs_Cell", "LFC_Depleted_vs_Cell"),
               where = "guide-level LFC table")

  df <- guide_lfc %>%
    dplyr::mutate(type = vapply(gene, pair_type, character(1L),
                                USE.NAMES = FALSE))

  lfc_cols <- paste0("LFC_", LFC_COLS)

  ntnt <- df %>% dplyr::filter(type == "NT_NT")
  ntnt_mean <- vapply(lfc_cols, function(col) mean(ntnt[[col]]),
                      numeric(1L))

  grouped <- df %>%
    dplyr::filter(type != "NT_NT") %>%
    dplyr::group_split(gene)

  rows <- lapply(grouped, function(grp) {
    gene_pair <- grp$gene[[1L]]
    ptype     <- grp$type[[1L]]
    top_n     <- if (ptype == "dual") top_dual else top_nt
    n_used    <- min(top_n, nrow(grp))

    out <- tibble::tibble(
      gene       = gene_pair,
      type       = ptype,
      n_guides   = nrow(grp),
      n_top_used = n_used
    )
    for (col in lfc_cols) {
      vals     <- sort(grp[[col]], decreasing = TRUE)[seq_len(n_used)]
      raw      <- mean(vals)
      corr     <- raw - ntnt_mean[[col]]
      label    <- sub("^LFC_", "", col)
      out[[paste0("mean_top_LFC_", label)]]               <- round(raw,  4L)
      out[[paste0("mean_top_LFC_", label, "_corrected")]] <- round(corr, 4L)
    }
    out
  })

  out <- dplyr::bind_rows(rows) %>%
    dplyr::arrange(type, gene)

  if (!is.null(out_path)) {
    dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)
    readr::write_csv(out, out_path)
    n_dual <- sum(out$type == "dual")
    n_snt  <- sum(out$type == "single_NT")
    message(sprintf("    Wrote %s (%d dual + %d single+NT)",
                    out_path, n_dual, n_snt))
  }
  attr(out, "ntnt_mean") <- ntnt_mean
  out
}

calc_gi <- function(gene_lfc, out_path = NULL) {
  require_cols(gene_lfc, c("gene", "type",
                           paste0("mean_top_LFC_", LFC_COLS, "_corrected")),
               where = "gene-pair LFC table")

  dual <- gene_lfc %>% dplyr::filter(type == "dual")
  snt  <- gene_lfc %>% dplyr::filter(type == "single_NT")

  single_lfc <- list()
  for (i in seq_len(nrow(snt))) {
    parts <- strsplit(snt$gene[[i]], "*", fixed = TRUE)[[1L]]
    gene <- if (parts[[2L]] == "NT") parts[[1L]] else parts[[2L]]
    vals <- vapply(LFC_COLS,
                   function(c) snt[[paste0("mean_top_LFC_", c, "_corrected")]][[i]],
                   numeric(1L))
    names(vals) <- LFC_COLS
    single_lfc[[gene]] <- vals
  }
  snt_gene_set <- snt$gene

  records <- vector("list", nrow(dual))
  missing <- character(0L)
  k <- 0L
  for (i in seq_len(nrow(dual))) {
    parts <- strsplit(dual$gene[[i]], "*", fixed = TRUE)[[1L]]
    gA <- parts[[1L]]; gB <- parts[[2L]]
    if (is.null(single_lfc[[gA]]) || is.null(single_lfc[[gB]])) {
      missing <- c(missing, dual$gene[[i]])
      next
    }
    rec <- tibble::tibble(
      gene  = dual$gene[[i]],
      geneA = gA,
      geneB = gB,
      singleA_pair = if (paste0(gA, "*NT") %in% snt_gene_set)
                       paste0(gA, "*NT") else paste0("NT*", gA),
      singleB_pair = if (paste0("NT*", gB) %in% snt_gene_set)
                       paste0("NT*", gB) else paste0(gB, "*NT")
    )
    for (c in LFC_COLS) {
      d_c  <- dual[[paste0("mean_top_LFC_", c, "_corrected")]][[i]]
      sA_c <- single_lfc[[gA]][[c]]
      sB_c <- single_lfc[[gB]][[c]]
      rec[[paste0("LFC_dual_",    c, "_corrected")]] <- round(d_c,  4L)
      rec[[paste0("LFC_singleA_", c, "_corrected")]] <- round(sA_c, 4L)
      rec[[paste0("LFC_singleB_", c, "_corrected")]] <- round(sB_c, 4L)
      rec[[paste0("GI_", c)]]                        <- round(d_c - sA_c - sB_c, 4L)
    }
    k <- k + 1L
    records[[k]] <- rec
  }
  length(records) <- k

  out <- dplyr::bind_rows(records) %>%
    dplyr::arrange(gene) %>%
    dplyr::select(gene, geneA, geneB, singleA_pair, singleB_pair,
                  paste0("LFC_dual_WT_vs_Cell_corrected"),
                  paste0("LFC_singleA_WT_vs_Cell_corrected"),
                  paste0("LFC_singleB_WT_vs_Cell_corrected"),
                  GI_WT_vs_Cell,
                  paste0("LFC_dual_Depleted_vs_Cell_corrected"),
                  paste0("LFC_singleA_Depleted_vs_Cell_corrected"),
                  paste0("LFC_singleB_Depleted_vs_Cell_corrected"),
                  GI_Depleted_vs_Cell)

  if (length(missing) > 0L) {
    message(sprintf("    WARNING: %d dual pairs lacked a single-NT match",
                    length(missing)))
  }
  if (!is.null(out_path)) {
    dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)
    readr::write_csv(out, out_path)
    message(sprintf("    Wrote %s (%d dual pairs with GI)",
                    out_path, nrow(out)))
  }
  attr(out, "missing") <- missing
  out
}
