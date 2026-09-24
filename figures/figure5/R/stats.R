suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
})

.boot_pval <- function(M) {
  apply(M, 1L, function(b) {
    b <- b[!is.na(b)]
    if (length(b) == 0L) return(NA_real_)
    n  <- length(b)
    p_above <- (0.5 + sum(b >= 0)) / (n + 1L)
    p_below <- (0.5 + sum(b <= 0)) / (n + 1L)

    min(1, 2 * min(p_above, p_below))
  })
}

permutation_gi <- function(guide_lfc, obs_gi, n_perm = 200L, seed = 123L,
                           top_dual = TOP_DUAL, top_nt = TOP_NT) {
  set.seed(seed)
  type <- vapply(guide_lfc$gene, pair_type, character(1L), USE.NAMES = FALSE)
  dual_idx <- which(type == "dual")

  obs_gi <- obs_gi %>%
    dplyr::mutate(delta_GI = GI_Depleted_vs_Cell - GI_WT_vs_Cell)

  obs_abs_WT    <- abs(obs_gi$GI_WT_vs_Cell)
  obs_abs_Dep   <- abs(obs_gi$GI_Depleted_vs_Cell)
  obs_abs_delta <- abs(obs_gi$delta_GI)

  count_WT    <- integer(nrow(obs_gi))
  count_Dep   <- integer(nrow(obs_gi))
  count_delta <- integer(nrow(obs_gi))

  for (b in seq_len(n_perm)) {
    perm <- guide_lfc
    perm$gene[dual_idx] <- sample(perm$gene[dual_idx])

    perm_gene_lfc <- calc_gene_lfc(perm,    top_dual = top_dual, top_nt = top_nt)
    perm_gi       <- calc_gi(perm_gene_lfc)

    m <- match(obs_gi$gene, perm_gi$gene)
    null_WT    <- abs(perm_gi$GI_WT_vs_Cell[m])
    null_Dep   <- abs(perm_gi$GI_Depleted_vs_Cell[m])
    null_delta <- abs(perm_gi$GI_Depleted_vs_Cell[m] - perm_gi$GI_WT_vs_Cell[m])

    null_WT[is.na(null_WT)]       <- 0
    null_Dep[is.na(null_Dep)]     <- 0
    null_delta[is.na(null_delta)] <- 0

    count_WT    <- count_WT    + (null_WT    >= obs_abs_WT)
    count_Dep   <- count_Dep   + (null_Dep   >= obs_abs_Dep)
    count_delta <- count_delta + (null_delta >= obs_abs_delta)
  }

  p_WT    <- (count_WT    + 1L) / (n_perm + 1L)
  p_Dep   <- (count_Dep   + 1L) / (n_perm + 1L)
  p_delta <- (count_delta + 1L) / (n_perm + 1L)

  tibble::tibble(
    gene                = obs_gi$gene,
    GI_WT_pvalue        = p_WT,
    GI_WT_qvalue        = stats::p.adjust(p_WT,    method = "BH"),
    GI_Depleted_pvalue  = p_Dep,
    GI_Depleted_qvalue  = stats::p.adjust(p_Dep,   method = "BH"),
    delta_GI_pvalue     = p_delta,
    delta_GI_qvalue     = stats::p.adjust(p_delta, method = "BH")
  )
}

bootstrap_gi_ci <- function(counts, gene_map, obs_gi,
                            n_boot = 200L, conf = 0.95, seed = 123L,
                            top_dual = TOP_DUAL, top_nt = TOP_NT) {
  set.seed(seed)
  sample_cols <- colnames(counts)
  groups <- classify_samples(sample_cols)
  stopifnot(length(groups$cell) > 0L,
            length(groups$wt) > 0L,
            length(groups$depleted) > 0L)

  n_pairs <- nrow(obs_gi)
  boot_WT    <- matrix(NA_real_, nrow = n_pairs, ncol = n_boot)
  boot_Dep   <- matrix(NA_real_, nrow = n_pairs, ncol = n_boot)
  boot_delta <- matrix(NA_real_, nrow = n_pairs, ncol = n_boot)

  for (b in seq_len(n_boot)) {
    resampled_cols <- c(
      sample(groups$cell,     length(groups$cell),     replace = TRUE),
      sample(groups$wt,       length(groups$wt),       replace = TRUE),
      sample(groups$depleted, length(groups$depleted), replace = TRUE)
    )
    cnt_b <- counts[, resampled_cols, drop = FALSE]
    lib <- colSums(cnt_b)
    cpm_b <- sweep(cnt_b, 2L, lib, "/") * 1e6
    log2cpm_b <- log2(cpm_b + 1)
    colnames(log2cpm_b) <- resampled_cols

    guide_lfc_b <- calc_guide_lfc(log2cpm_b, resampled_cols, gene_map)
    gene_lfc_b  <- calc_gene_lfc(guide_lfc_b, top_dual = top_dual,
                                 top_nt = top_nt)
    gi_b        <- calc_gi(gene_lfc_b)

    m <- match(obs_gi$gene, gi_b$gene)
    boot_WT[,    b] <- gi_b$GI_WT_vs_Cell[m]
    boot_Dep[,   b] <- gi_b$GI_Depleted_vs_Cell[m]
    boot_delta[, b] <- gi_b$GI_Depleted_vs_Cell[m] - gi_b$GI_WT_vs_Cell[m]
  }

  lo <- (1 - conf) / 2
  hi <- 1 - lo
  qfun <- function(M) apply(M, 1L, stats::quantile,
                            probs = c(lo, hi), na.rm = TRUE)

  q_WT    <- qfun(boot_WT)
  q_Dep   <- qfun(boot_Dep)
  q_delta <- qfun(boot_delta)

  p_WT    <- .boot_pval(boot_WT)
  p_Dep   <- .boot_pval(boot_Dep)
  p_delta <- .boot_pval(boot_delta)

  tibble::tibble(
    gene                  = obs_gi$gene,
    GI_WT_ci_lo           = q_WT[1L, ],
    GI_WT_ci_hi           = q_WT[2L, ],
    GI_Depleted_ci_lo     = q_Dep[1L, ],
    GI_Depleted_ci_hi     = q_Dep[2L, ],
    delta_GI_ci_lo        = q_delta[1L, ],
    delta_GI_ci_hi        = q_delta[2L, ],
    GI_WT_pvalue          = p_WT,
    GI_WT_qvalue          = stats::p.adjust(p_WT,    method = "BH"),
    GI_Depleted_pvalue    = p_Dep,
    GI_Depleted_qvalue    = stats::p.adjust(p_Dep,   method = "BH"),
    delta_GI_pvalue       = p_delta,
    delta_GI_qvalue       = stats::p.adjust(p_delta, method = "BH")
  )
}

.boot_one_screen_delta <- function(counts, gene_map, genes_ref,
                                    top_dual, top_nt) {
  cols   <- colnames(counts)
  groups <- classify_samples(cols)
  rs <- c(
    sample(groups$cell,     length(groups$cell),     replace = TRUE),
    sample(groups$wt,       length(groups$wt),       replace = TRUE),
    sample(groups$depleted, length(groups$depleted), replace = TRUE)
  )
  cb  <- counts[, rs, drop = FALSE]
  lib <- colSums(cb)
  l2  <- log2(sweep(cb, 2L, lib, "/") * 1e6 + 1)
  colnames(l2) <- rs
  g   <- calc_guide_lfc(l2, rs, gene_map)
  gl  <- calc_gene_lfc(g, top_dual = top_dual, top_nt = top_nt)
  gi  <- calc_gi(gl)
  d   <- gi$GI_Depleted_vs_Cell - gi$GI_WT_vs_Cell
  d[match(genes_ref, gi$gene)]
}

bootstrap_ddgi <- function(counts_cx, gene_map_cx,
                           counts_mrc, gene_map_mrc, obs,
                           n_boot = 200L, conf = 0.95, seed = 123L,
                           top_dual = TOP_DUAL, top_nt = TOP_NT) {
  set.seed(seed)
  genes_ref <- obs$gene
  np <- length(genes_ref)
  bc <- matrix(NA_real_, np, n_boot)
  bm <- matrix(NA_real_, np, n_boot)

  for (b in seq_len(n_boot)) {
    bc[, b] <- .boot_one_screen_delta(counts_cx, gene_map_cx, genes_ref,
                                      top_dual, top_nt)
    bm[, b] <- .boot_one_screen_delta(counts_mrc, gene_map_mrc, genes_ref,
                                      top_dual, top_nt)
  }
  bd <- bc - bm

  lo <- (1 - conf) / 2
  hi <- 1 - lo
  qfun <- function(M) apply(M, 1L, stats::quantile,
                            probs = c(lo, hi), na.rm = TRUE)
  q_dCx <- qfun(bc)
  q_dMrc <- qfun(bm)
  q_dd  <- qfun(bd)

  p_dd <- .boot_pval(bd)

  tibble::tibble(
    gene             = genes_ref,
    dGI_Cx_ci_lo     = q_dCx[1L, ],
    dGI_Cx_ci_hi     = q_dCx[2L, ],
    dGI_Mrc_ci_lo    = q_dMrc[1L, ],
    dGI_Mrc_ci_hi    = q_dMrc[2L, ],
    ddGI_ci_lo       = q_dd[1L, ],
    ddGI_ci_hi       = q_dd[2L, ],
    ddGI_pvalue      = p_dd,
    ddGI_qvalue      = stats::p.adjust(p_dd, method = "BH")
  )
}
