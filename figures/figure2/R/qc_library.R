R1_WHITELIST <- 233684L

qc_gini <- function(x) {
  x <- sort(x[x > 0])
  n <- length(x)
  if (n == 0) return(NA_real_)
  (2 * sum(seq_len(n) * x) / (n * sum(x))) - (n + 1) / n
}

qc_lorenz <- function(x) {
  x <- sort(x[x > 0])
  data.frame(
    p_constructs = c(0, seq_along(x) / length(x)),
    p_reads      = c(0, cumsum(x) / sum(x))
  )
}

qc_gene_pair <- function(construct) {
  parts <- strsplit(construct, "*", fixed = TRUE)
  vapply(parts, function(p) {
    if (length(p) < 2) return(NA_character_)
    paste(sub("_\\d+$", "", p[1]), sub("_\\d+$", "", p[2]), sep = "*")
  }, character(1))
}

qc_gene_side <- function(construct, side) {
  parts <- strsplit(construct, "*", fixed = TRUE)
  vapply(parts, function(p) {
    if (length(p) < 2) return(NA_character_)
    sub("_\\d+$", "", p[side])
  }, character(1))
}

qc_design_genes <- function(count_csv = raw("R1_construct_counts.csv")) {
  ids <- read.csv(count_csv, check.names = FALSE, stringsAsFactors = FALSE)$sgRNA
  ids <- ids[grepl("*", ids, fixed = TRUE)]
  list(pos1 = unique(qc_gene_side(ids, 1)),
       pos2 = unique(qc_gene_side(ids, 2)))
}

qc_classify <- function(construct, design) {
  a <- qc_gene_side(construct, 1)
  b <- qc_gene_side(construct, 2)
  ifelse(a %in% design$pos1 & b %in% design$pos2, "canonical",
         ifelse(a %in% design$pos2 & b %in% design$pos1, "reciprocal",
                "within_program"))
}

qc_metrics <- function(reads, pairs, whitelist = R1_WHITELIST) {
  n <- length(reads)
  q <- quantile(reads, c(0.10, 0.90), names = FALSE)
  ord <- sort(reads, decreasing = TRUE)
  top10 <- sum(ord[seq_len(max(1, floor(n * 0.10)))]) / sum(ord)
  u <- unique(pairs[!is.na(pairs)])

  is_nt  <- function(g) toupper(g) == "NT"
  sides  <- strsplit(u, "*", fixed = TRUE)
  n_nt   <- vapply(sides, function(p) sum(is_nt(p)), integer(1))
  dual   <- n_nt == 0L
  single <- n_nt == 1L
  ntnt   <- n_nt == 2L
  list(
    n_constructs = n,
    pct_detected = 100 * n / whitelist,
    n_pairs      = length(u),
    n_dual_pairs  = sum(dual),
    n_single_ctrl = sum(single),
    n_ntnt        = sum(ntnt),
    median_reads = median(reads),
    gini         = qc_gini(reads),
    skew_ratio   = if (q[1] > 0) q[2] / q[1] else NA_real_,
    top10_share  = 100 * top10,
    total_reads  = sum(reads)
  )
}
