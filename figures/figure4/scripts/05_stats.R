suppressPackageStartupMessages({
  source(here::here("R", "paths.R"))
  source(here::here("R", "helpers.R"))
  source(here::here("R", "gi_pipeline.R"))
  source(here::here("R", "stats.R"))
})

SEED      <- 123L
N_PERM    <- 200L
N_BOOT    <- 200L
CONF      <- 0.95
TOP_DUAL_GRID <- c(3L, 5L, 10L)
TOP_NT_GRID   <- c(10L, 18L, 25L)

set.seed(SEED)

message("\n=== 05_stats ===\n")

guide_lfc <- readr::read_csv(tbl("guideLFC_Cx3cr1.csv"),
                             show_col_types = FALSE, progress = FALSE)
gi_obs    <- readr::read_csv(tbl("GI_Cx3cr1.csv"),
                             show_col_types = FALSE, progress = FALSE)
require_cols(guide_lfc, c("sgRNA", "gene",
                          "LFC_WT_vs_Cell", "LFC_Depleted_vs_Cell"),
             where = "output/tables/guideLFC_Cx3cr1.csv")
require_cols(gi_obs, c("gene", "GI_WT_vs_Cell", "GI_Depleted_vs_Cell"),
             where = "output/tables/GI_Cx3cr1.csv")

stat_cols <- c(
  "GI_WT_pvalue",      "GI_WT_qvalue",
  "GI_Depleted_pvalue","GI_Depleted_qvalue",
  "delta_GI_pvalue",   "delta_GI_qvalue",
  "GI_WT_ci_lo",       "GI_WT_ci_hi",
  "GI_Depleted_ci_lo", "GI_Depleted_ci_hi",
  "delta_GI_ci_lo",    "delta_GI_ci_hi"
)
gi_obs <- gi_obs[, setdiff(names(gi_obs), stat_cols), drop = FALSE]

message(sprintf("[05_stats] Replicate bootstrap (n_boot = %d, conf = %.2f) ...",
                N_BOOT, CONF))
t0 <- Sys.time()
normd <- load_and_normalize(raw("Cx3cr1_raw_counts.csv"))
boot <- bootstrap_gi_ci(normd$counts, normd$gene_map, gi_obs,
                        n_boot = N_BOOT, conf = CONF, seed = SEED)
message(sprintf("    -> %.1fs", as.numeric(difftime(Sys.time(), t0, units = "secs"))))

n_sig_05_WT  <- sum(boot$GI_WT_qvalue       < 0.05, na.rm = TRUE)
n_sig_10_WT  <- sum(boot$GI_WT_qvalue       < 0.10, na.rm = TRUE)
n_sig_05_Dep <- sum(boot$GI_Depleted_qvalue < 0.05, na.rm = TRUE)
n_sig_10_Dep <- sum(boot$GI_Depleted_qvalue < 0.10, na.rm = TRUE)
n_sig_05_dlt <- sum(boot$delta_GI_qvalue    < 0.05, na.rm = TRUE)
n_sig_10_dlt <- sum(boot$delta_GI_qvalue    < 0.10, na.rm = TRUE)
message(sprintf("    GI_WT       FDR<0.05: %3d   FDR<0.10: %3d", n_sig_05_WT,  n_sig_10_WT))
message(sprintf("    GI_Depleted FDR<0.05: %3d   FDR<0.10: %3d", n_sig_05_Dep, n_sig_10_Dep))
message(sprintf("    delta_GI    FDR<0.05: %3d   FDR<0.10: %3d", n_sig_05_dlt, n_sig_10_dlt))

gi_aug <- gi_obs %>%
  dplyr::left_join(boot, by = "gene") %>%
  dplyr::mutate(dplyr::across(dplyr::all_of(stat_cols),
                              ~ ifelse(is.na(.x), NA_real_, round(.x, 4L))))

readr::write_csv(gi_aug, tbl("GI_Cx3cr1.csv"))
message(sprintf("    Wrote %s (%d rows, %d cols)",
                tbl("GI_Cx3cr1.csv"), nrow(gi_aug), ncol(gi_aug)))

message(sprintf("[05_stats] (Optional) permutation null (n_perm = %d) ...",
                N_PERM))
t0 <- Sys.time()
perm <- permutation_gi(guide_lfc, gi_obs, n_perm = N_PERM, seed = SEED)
message(sprintf("    -> %.1fs", as.numeric(difftime(Sys.time(), t0, units = "secs"))))
readr::write_csv(perm, tbl("permutation_GI_Cx3cr1.csv"))
message(sprintf("    Wrote %s", tbl("permutation_GI_Cx3cr1.csv")))

message(sprintf("[05_stats] Top-N sensitivity sweep (%d x %d settings) ...",
                length(TOP_DUAL_GRID), length(TOP_NT_GRID)))
t0 <- Sys.time()
sens <- topN_sensitivity(guide_lfc,
                         top_dual_grid = TOP_DUAL_GRID,
                         top_nt_grid   = TOP_NT_GRID)
rho  <- topN_rank_correlation(sens, ref_top_dual = TOP_DUAL,
                              ref_top_nt = TOP_NT)
message(sprintf("    -> %.1fs", as.numeric(difftime(Sys.time(), t0, units = "secs"))))
print(rho, n = nrow(rho))

readr::write_csv(sens, tbl("topN_sensitivity.csv"))
readr::write_csv(rho,  tbl("topN_rank_correlation.csv"))
message(sprintf("    Wrote %s (%d rows)",
                tbl("topN_sensitivity.csv"), nrow(sens)))
message(sprintf("    Wrote %s (%d rows)",
                tbl("topN_rank_correlation.csv"), nrow(rho)))

if (!interactive()) log_session("05_stats")
invisible(list(gi = gi_aug, boot = boot, perm = perm,
               sensitivity = sens, rho = rho))
