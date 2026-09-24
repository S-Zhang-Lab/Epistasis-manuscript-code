suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
})

source(here::here("R", "paths.R"))
source(here::here("R", "helpers.R"))

source(here::here("R", "theme_pub.R"))

OUT      <- OUT_TABLES

FIG_DIR  <- OUT_SUPP

FDR <- 0.05
ENRICH_FDR <- 0.10
MODULE_COLORS <- c("0" = "#E41A1C", "1" = "#377EB8", "2" = "#4DAF4A")

set.seed(42)

message("Loading GI and single-KO data")
gi <- read.csv(file.path(OUT, "GI_scores_all_dualKO.csv"),
               stringsAsFactors = FALSE)
sk <- read.csv(file.path(OUT, "singleKO_reference_table.csv"),
               stringsAsFactors = FALSE)

sig <- gi |> filter(fdr_GI < FDR)
message(sprintf("Significant GI pairs: %d (synergy=%d, buffering=%d)",
                nrow(sig), sum(sig$GI_score > 0), sum(sig$GI_score < 0)))

genes_in_net <- sort(unique(c(toupper(sig$gene_A), toupper(sig$gene_B))))
n_genes <- length(genes_in_net)
gene_idx <- setNames(seq_along(genes_in_net), genes_in_net)
sk_lfc <- setNames(sk$meta_lfc, toupper(sk$gene))

W <- matrix(0, n_genes, n_genes, dimnames = list(genes_in_net, genes_in_net))
for (i in seq_len(nrow(sig))) {
  a <- toupper(sig$gene_A[i])
  b <- toupper(sig$gene_B[i])
  if (a %in% genes_in_net && b %in% genes_in_net) {
    w <- abs(sig$GI_Z[i])
    W[a, b] <- max(W[a, b], w)
    W[b, a] <- W[a, b]
  }
}
degree_w <- rowSums(W)
message(sprintf("Genes in network: %d; edges: %d; isolated: %d",
                n_genes, sum(W > 0) / 2, sum(degree_w == 0)))

spectral_embedding <- function(W, k) {
  d <- rowSums(W)
  d_inv_sqrt <- ifelse(d > 0, 1 / sqrt(d), 0)
  L_norm <- diag(nrow(W)) - diag(d_inv_sqrt) %*% W %*% diag(d_inv_sqrt)
  eig <- eigen(L_norm, symmetric = TRUE)
  ord <- order(eig$values)
  eigvals <- eig$values[ord]
  eigvecs <- eig$vectors[, ord, drop = FALSE]
  U <- eigvecs[, seq_len(k), drop = FALSE]

  for (j in seq_len(ncol(U))) {
    lead <- which.max(abs(U[, j]))
    if (U[lead, j] < 0) U[, j] <- -U[, j]
  }
  nr <- sqrt(rowSums(U^2))
  nr[nr == 0] <- 1
  U <- U / nr
  list(eigvals = eigvals, U = U)
}

eig_all <- spectral_embedding(W, min(20, n_genes))
gaps <- diff(eig_all$eigvals[seq_len(min(15, length(eig_all$eigvals)))])
k_opt <- if (length(gaps) > 1) which.max(gaps[-1]) + 1 else 3
k_opt <- max(3, min(k_opt, 8))
message("Eigengap-selected k: ", k_opt)

emb <- spectral_embedding(W, k_opt)
in_net_mask <- degree_w > 0
km <- kmeans(emb$U[in_net_mask, , drop = FALSE], centers = k_opt,
             nstart = 10, iter.max = 300)

labels_full <- rep(-1L, n_genes)
labels_full[in_net_mask] <- km$cluster - 1L
module_sizes <- sort(table(labels_full[labels_full >= 0]), decreasing = TRUE)
remap <- setNames(seq_along(module_sizes) - 1L, names(module_sizes))
labels_final <- ifelse(labels_full >= 0, as.integer(remap[as.character(labels_full)]), -1L)
n_modules <- length(module_sizes)

for (m in seq_len(n_modules) - 1L) {
  mem <- genes_in_net[labels_final == m]
  message(sprintf("  Module %d: n=%d  %s", m, length(mem),
                  paste(head(mem, 8), collapse = ", ")))
}

mod_df <- data.frame(
  gene = genes_in_net,
  module = labels_final,
  meta_lfc = unname(sk_lfc[genes_in_net]),
  stringsAsFactors = FALSE
)
mod_df$meta_lfc[is.na(mod_df$meta_lfc)] <- 0
write.csv(mod_df, file.path(OUT, "module_assignments.csv"),
          row.names = FALSE, quote = FALSE)

GENE_SETS <- list(
  "Tumor suppressors /\nDNA damage" = c(
    "PTEN", "BRCA1", "BRCA2", "RB1", "SMAD4", "TRP53", "NF2",
    "APC", "KMT2C", "CHEK2", "STK11", "TSC2"
  ),
  "ECM remodeling /\nInvasion" = c(
    "MMP3", "MMP12", "CDH1", "HAS3", "VCAN", "COL17A1",
    "ADAM12", "ADAM20", "LYVE1", "CX3CL1", "PAPPA2"
  ),
  "Immune checkpoint /\nEvasion" = c(
    "PDCD1LG2", "BTLA", "TIMD4", "CLEC4G", "FCRL5", "FCRL6",
    "TLR1", "TLR4", "TLR7", "MPEG1", "CD38"
  ),
  "Cytokines /\nImmune signaling" = c(
    "IL6", "IL6ST", "CXCL1", "CXCL2", "CXCL5", "CXCR2", "CXCR5",
    "NRG1", "TGFA", "HGF", "EREG", "ANGPTL3", "TNFSF15", "TNFSF8",
    "CCR2", "KIT"
  ),
  "Lipid metabolism /\nAdipokines" = c(
    "PLIN1", "GPD1", "THRSP", "ADH4", "PCK1", "UPP2",
    "CIDEC", "KLB", "LEP", "AADAC"
  ),
  "Epithelial identity /\nDifferentiation" = c(
    "KRT14", "KRT15", "KRT6B", "KRT23", "DSC3", "DSG3", "CDH1",
    "CLDN8", "SERPINB5", "EHF", "CLCA2", "KLK7", "CALML3"
  ),
  "Ion channels /\nMembrane excitability" = c(
    "KCNJ13", "KCNJ16", "KCNS2", "SCN11A", "CACNB4",
    "GRIK1", "GRM7", "TRPM3", "KCNRG"
  ),
  "PI3K-AKT-mTOR" = c(
    "PTEN", "PIK3C2G", "APC", "TSC2", "STK11", "BRCA1", "BRCA2"
  )
)
GENE_SETS <- lapply(GENE_SETS, toupper)

universe <- genes_in_net
M <- length(universe)
enrich_rows <- list()
for (m in seq_len(n_modules) - 1L) {
  mod_genes <- genes_in_net[labels_final == m]
  N <- length(mod_genes)
  pvals <- numeric(length(GENE_SETS))
  module_rows <- vector("list", length(GENE_SETS))
  for (i in seq_along(GENE_SETS)) {
    gs_name <- names(GENE_SETS)[i]
    gs <- intersect(GENE_SETS[[i]], universe)
    n_gs <- length(gs)
    k <- length(intersect(mod_genes, gs))
    pval <- if (k > 0) phyper(k - 1, n_gs, M - n_gs, N,
                              lower.tail = FALSE) else 1
    odds <- if (n_gs > 0 && N > 0) (k / N) / (n_gs / M) else 0
    pvals[i] <- pval
    module_rows[[i]] <- data.frame(
      module = m, gene_set = gs_name, n_module = N, n_geneset = n_gs,
      n_overlap = k, odds_ratio = odds, pval = pval,
      stringsAsFactors = FALSE
    )
  }
  module_df <- bind_rows(module_rows)
  module_df$fdr <- p.adjust(pvals, method = "BH")
  enrich_rows[[as.character(m)]] <- module_df
}
enrich_df <- bind_rows(enrich_rows)
write.csv(enrich_df, file.path(OUT, "module_enrichment.csv"),
          row.names = FALSE, quote = FALSE)

emb3 <- spectral_embedding(W, 3)
n_in_net <- sum(in_net_mask & labels_final %in% 0:2)
message(sprintf(
  "Spectral: n=%d genes, %d significant pairs at FDR<%.2f; lambda2=%.3f lambda3=%.3f lambda4=%.3f; k=%d by eigengap",
  n_in_net, sum(W > 0) / 2, FDR,
  emb3$eigvals[2], emb3$eigvals[3], emb3$eigvals[4], k_opt))

message("Computing within/cross-module GI stratification (for Fig2H forest)")
mod_lookup <- setNames(mod_df$module, toupper(mod_df$gene))
gi_mod <- gi |>
  mutate(gene_A_u = toupper(gene_A), gene_B_u = toupper(gene_B),
         mod_A = mod_lookup[gene_A_u],
         mod_B = mod_lookup[gene_B_u]) |>
  filter(!is.na(mod_A), !is.na(mod_B), mod_A >= 0, mod_B >= 0) |>
  mutate(mod_A = as.integer(mod_A), mod_B = as.integer(mod_B),
         within = mod_A == mod_B)

within_z <- gi_mod$GI_Z[gi_mod$within]
cross_z <- gi_mod$GI_Z[!gi_mod$within]
mw <- wilcox.test(within_z, cross_z, exact = FALSE)
message(sprintf("Within mean %.3f (n=%d); cross mean %.3f (n=%d); Wilcoxon p=%.2e",
                mean(within_z), length(within_z), mean(cross_z), length(cross_z),
                mw$p.value))

message("Module assignment + within/cross-module GI computation complete")

log_session("02_r1_modules")
