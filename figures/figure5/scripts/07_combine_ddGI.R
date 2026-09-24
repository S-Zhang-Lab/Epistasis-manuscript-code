suppressPackageStartupMessages({
  source(here::here("R", "paths.R"))
  source(here::here("R", "helpers.R"))
  source(here::here("R", "gi_pipeline.R"))
  source(here::here("R", "stats.R"))
  for (pkg in c("dplyr", "readr")) library(pkg, character.only = TRUE)
})

SEED   <- 123L
N_BOOT <- 200L
Q_SIG  <- 0.05
set.seed(SEED)

message("\n=== 07_combine_ddGI ===\n")

load_screen_gi <- function(screen) {
  g <- readr::read_csv(tbl(sprintf("GI_%s.csv", screen)),
                       show_col_types = FALSE, progress = FALSE)
  has_q <- "delta_GI_qvalue" %in% names(g)
  tibble::tibble(
    gene  = g$gene,
    geneA = g$geneA,
    geneB = g$geneB,
    dGI   = g$GI_Depleted_vs_Cell - g$GI_WT_vs_Cell,
    dGI_q = if (has_q) g$delta_GI_qvalue else NA_real_
  )
}

cx  <- load_screen_gi("Cx3cr1")
mrc <- load_screen_gi("Mrc1")

comb <- dplyr::inner_join(
  dplyr::rename(cx,  dGI_Cx = dGI,  q_dGI_Cx = dGI_q),
  dplyr::select(mrc, gene, dGI_Mrc = dGI, q_dGI_Mrc = dGI_q),
  by = "gene"
) %>%
  dplyr::mutate(ddGI = dGI_Cx - dGI_Mrc)

message(sprintf("[07_combine_ddGI] %d pairs present in both screens", nrow(comb)))

DD_COLS <- c("gene", "dGI_Cx_ci_lo", "dGI_Cx_ci_hi", "dGI_Mrc_ci_lo",
             "dGI_Mrc_ci_hi", "ddGI_ci_lo", "ddGI_ci_hi", "ddGI_pvalue",
             "ddGI_qvalue")
cache_path <- tbl("ddGI_combined.csv")
prov_path  <- tbl("ddGI_bootstrap_provenance.txt")

.dd_fingerprint <- function() {
  tf <- tempfile(); on.exit(unlink(tf), add = TRUE)
  saveRDS(list(
    raw_md5 = unname(tools::md5sum(c(raw("Cx3cr1_raw_counts.csv"),
                                     raw("Mrc1_raw_counts.csv")))),
    gene    = comb$gene,
    dGI_Cx  = comb$dGI_Cx,
    dGI_Mrc = comb$dGI_Mrc,
    n_boot  = N_BOOT,
    seed    = SEED
  ), tf)
  unname(tools::md5sum(tf))
}
fp <- .dd_fingerprint()

dd <- NULL
if (file.exists(cache_path) && file.exists(prov_path)) {
  prev    <- readr::read_csv(cache_path, show_col_types = FALSE, progress = FALSE)
  prev_fp <- tryCatch(trimws(readLines(prov_path, warn = FALSE))[1],
                      error = function(e) NA_character_)
  if (all(DD_COLS %in% names(prev)) && setequal(prev$gene, comb$gene) &&
      identical(prev_fp, fp)) {
    dd <- prev[, DD_COLS]
    message("[07_combine_ddGI] reusing cached ΔΔGI bootstrap (provenance fingerprint matches)")
  } else {
    message("[07_combine_ddGI] cache present but inputs/settings changed -> recomputing")
  }
}
if (is.null(dd)) {
  message(sprintf("[07_combine_ddGI] joint ΔΔGI bootstrap (%d iters, both screens) ...",
                  N_BOOT))
  nx <- load_and_normalize(raw("Cx3cr1_raw_counts.csv"))
  nm <- load_and_normalize(raw("Mrc1_raw_counts.csv"))
  t0 <- Sys.time()
  dd <- bootstrap_ddgi(nx$counts, nx$gene_map, nm$counts, nm$gene_map,
                       obs = comb, n_boot = N_BOOT, seed = SEED)
  message(sprintf("    -> %.1fs", as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  writeLines(fp, prov_path)
}

comb <- comb %>% dplyr::left_join(dd, by = "gene")

comb <- comb %>%
  dplyr::mutate(
    sig_Cx  = !is.na(q_dGI_Cx)  & q_dGI_Cx  < Q_SIG,
    sig_Mrc = !is.na(q_dGI_Mrc) & q_dGI_Mrc < Q_SIG,
    sig_dd  = !is.na(ddGI_qvalue) & ddGI_qvalue < Q_SIG,
    genuine = sig_Cx | sig_Mrc,

    category = dplyr::case_when(
      !genuine                                          ~ "Not niche-conditional",
      !sig_dd                                           ~ "Conserved across cell type",
      sig_Cx & sig_Mrc & sign(dGI_Cx) != sign(dGI_Mrc) ~ "Sign-divergent",
      ddGI < 0                                          ~ "Stronger in Mrc1",
      TRUE                                              ~ "Stronger in Cx3cr1"),
    gene_label = gsub("\\*", " x ", gene)
  )

readr::write_csv(comb, tbl("ddGI_combined.csv"))
message(sprintf("    Wrote %s (%d rows, %d cols)",
                tbl("ddGI_combined.csv"), nrow(comb), ncol(comb)))

n_total    <- nrow(comb)
n_genuine  <- sum(comb$genuine)
n_celltype <- sum(comb$genuine & comb$sig_dd)
n_dd_all   <- sum(comb$sig_dd)
n_signdiv  <- sum(comb$category == "Sign-divergent")
pct_among_genuine <- if (n_genuine > 0) 100 * n_celltype / n_genuine else NA

cat(sprintf(paste0(
  "\n------------------- ΔΔGI HEADLINE -------------------\n",
  "  pairs (both screens)                : %d\n",
  "  genuine niche-conditional (ΔGI q<%.2f in >=1 screen): %d\n",
  "  cell-type-specific (ΔΔGI q<%.2f)    : %d  (all pairs)\n",
  "  cell-type-specific AND genuine      : %d\n",
  "  => pervasiveness (cell-type-specific / genuine): %.0f%%\n",
  "  sign-divergent (opposite by cell type): %d\n",
  "-----------------------------------------------------\n"),
  n_total, Q_SIG, n_genuine, Q_SIG, n_dd_all, n_celltype, pct_among_genuine,
  n_signdiv))

cat("\n  Category breakdown:\n")
print(comb %>% dplyr::count(category, sort = TRUE))

hero <- comb %>% dplyr::filter(gene == "Pten*Cdh1")
if (nrow(hero) == 1) {
  cat(sprintf(paste0(
    "\n  HERO  Pten x Cdh1:\n",
    "    dGI_Cx  = %+.2f  [%.2f, %.2f]  q=%.3f\n",
    "    dGI_Mrc = %+.2f  [%.2f, %.2f]  q=%.3f\n",
    "    ddGI    = %+.2f  [%.2f, %.2f]  q=%.3f   -> %s\n\n"),
    hero$dGI_Cx,  hero$dGI_Cx_ci_lo,  hero$dGI_Cx_ci_hi,  hero$q_dGI_Cx,
    hero$dGI_Mrc, hero$dGI_Mrc_ci_lo, hero$dGI_Mrc_ci_hi, hero$q_dGI_Mrc,
    hero$ddGI,    hero$ddGI_ci_lo,    hero$ddGI_ci_hi,    hero$ddGI_qvalue,
    hero$category))
} else {
  cat("\n  HERO  Pten x Cdh1 not found in combined table.\n\n")
}

if (!interactive()) log_session("07_combine_ddGI")
