suppressPackageStartupMessages({
  source(here::here("R", "paths.R"))
  source(here::here("R", "helpers.R"))
  source(here::here("R", "gi_pipeline.R"))
  source(here::here("R", "stats.R"))
})

SEED   <- 123L
N_BOOT <- 200L
CONF   <- 0.95
DO_PERM <- FALSE
N_PERM <- 200L
set.seed(SEED)

message("\n=== 06_stats_GI ===\n")

STAT_COLS <- c(
  "GI_WT_pvalue", "GI_WT_qvalue",
  "GI_Depleted_pvalue", "GI_Depleted_qvalue",
  "delta_GI_pvalue", "delta_GI_qvalue",
  "GI_WT_ci_lo", "GI_WT_ci_hi",
  "GI_Depleted_ci_lo", "GI_Depleted_ci_hi",
  "delta_GI_ci_lo", "delta_GI_ci_hi"
)

run_stats <- function(screen) {
  message(sprintf("\n--- %s ---", screen))
  gi_obs <- readr::read_csv(tbl(sprintf("GI_%s.csv", screen)),
                            show_col_types = FALSE, progress = FALSE)
  require_cols(gi_obs, c("gene", "GI_WT_vs_Cell", "GI_Depleted_vs_Cell"),
               where = sprintf("GI_%s.csv", screen))
  gi_obs <- gi_obs[, setdiff(names(gi_obs), STAT_COLS), drop = FALSE]

  normd <- load_and_normalize(raw(sprintf("%s_raw_counts.csv", screen)))
  t0 <- Sys.time()
  boot <- bootstrap_gi_ci(normd$counts, normd$gene_map, gi_obs,
                          n_boot = N_BOOT, conf = CONF, seed = SEED)
  message(sprintf("    bootstrap (%d iters): %.1fs", N_BOOT,
                  as.numeric(difftime(Sys.time(), t0, units = "secs"))))

  message(sprintf("    GI_WT       FDR<0.05: %3d", sum(boot$GI_WT_qvalue       < 0.05, na.rm = TRUE)))
  message(sprintf("    GI_Depleted FDR<0.05: %3d", sum(boot$GI_Depleted_qvalue < 0.05, na.rm = TRUE)))
  message(sprintf("    delta_GI    FDR<0.05: %3d", sum(boot$delta_GI_qvalue    < 0.05, na.rm = TRUE)))

  gi_aug <- gi_obs %>%
    dplyr::left_join(boot, by = "gene") %>%
    dplyr::mutate(dplyr::across(dplyr::all_of(STAT_COLS),
                                ~ ifelse(is.na(.x), NA_real_, round(.x, 4L))))
  readr::write_csv(gi_aug, tbl(sprintf("GI_%s.csv", screen)))
  message(sprintf("    Wrote %s (%d rows, %d cols)",
                  tbl(sprintf("GI_%s.csv", screen)), nrow(gi_aug), ncol(gi_aug)))

  if (DO_PERM) {
    guide_lfc <- readr::read_csv(tbl(sprintf("guideLFC_%s.csv", screen)),
                                 show_col_types = FALSE, progress = FALSE)
    perm <- permutation_gi(guide_lfc, gi_obs, n_perm = N_PERM, seed = SEED)
    readr::write_csv(perm, tbl(sprintf("permutation_GI_%s.csv", screen)))
    message(sprintf("    Wrote %s", tbl(sprintf("permutation_GI_%s.csv", screen))))
  }
  invisible(gi_aug)
}

invisible(lapply(c("Cx3cr1", "Mrc1"), run_stats))

if (!interactive()) log_session("06_stats_GI")
