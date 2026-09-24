suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
})

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

  pval <- function(M) {
    apply(M, 1L, function(b) {
      b <- b[!is.na(b)]
      if (length(b) == 0L) return(NA_real_)
      n  <- length(b)
      p_above <- (0.5 + sum(b >= 0)) / (n + 1L)
      p_below <- (0.5 + sum(b <= 0)) / (n + 1L)
      2 * min(p_above, p_below)
    })
  }
  p_WT    <- pval(boot_WT)
  p_Dep   <- pval(boot_Dep)
  p_delta <- pval(boot_delta)

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

topN_sensitivity <- function(guide_lfc,
                             top_dual_grid = c(3L, 5L, 10L),
                             top_nt_grid   = c(10L, 18L, 25L)) {
  combos <- expand.grid(top_dual = top_dual_grid,
                        top_nt   = top_nt_grid,
                        KEEP.OUT.ATTRS = FALSE,
                        stringsAsFactors = FALSE)
  out <- vector("list", nrow(combos))
  for (i in seq_len(nrow(combos))) {
    gene_lfc_i <- calc_gene_lfc(guide_lfc,
                                top_dual = combos$top_dual[i],
                                top_nt   = combos$top_nt[i])
    gi_i <- calc_gi(gene_lfc_i)
    out[[i]] <- tibble::tibble(
      top_dual = combos$top_dual[i],
      top_nt   = combos$top_nt[i],
      gene     = gi_i$gene,
      GI_WT_vs_Cell       = gi_i$GI_WT_vs_Cell,
      GI_Depleted_vs_Cell = gi_i$GI_Depleted_vs_Cell,
      delta_GI            = gi_i$GI_Depleted_vs_Cell - gi_i$GI_WT_vs_Cell
    )
  }
  dplyr::bind_rows(out)
}

topN_rank_correlation <- function(sensitivity, ref_top_dual = 5L,
                                  ref_top_nt = 18L) {
  ref <- sensitivity %>%
    dplyr::filter(top_dual == ref_top_dual, top_nt == ref_top_nt)
  if (nrow(ref) == 0L) {
    stop("Reference (top_dual=", ref_top_dual, ", top_nt=", ref_top_nt,
         ") not present in the sensitivity table.", call. = FALSE)
  }

  sensitivity %>%
    dplyr::group_by(top_dual, top_nt) %>%
    dplyr::summarise(
      rho_GI_WT       = stats::cor(GI_WT_vs_Cell[match(ref$gene, gene)],
                                   ref$GI_WT_vs_Cell, method = "spearman"),
      rho_GI_Depleted = stats::cor(GI_Depleted_vs_Cell[match(ref$gene, gene)],
                                   ref$GI_Depleted_vs_Cell, method = "spearman"),
      rho_delta_GI    = stats::cor(delta_GI[match(ref$gene, gene)],
                                   ref$delta_GI, method = "spearman"),
      .groups = "drop"
    )
}
