suppressPackageStartupMessages({
  source(here::here("R", "paths.R"))
  source(here::here("R", "helpers.R"))
  source(here::here("R", "gi_pipeline.R"))
})

SEED <- 123L
set.seed(SEED)

message("\n=== 02_compute_GI ===\n")

SCREENS <- c("Cx3cr1", "Mrc1")

run_screen <- function(screen) {
  message(sprintf("\n--- %s ---", screen))
  normd <- load_and_normalize(raw(sprintf("%s_raw_counts.csv", screen)))

  guide_lfc <- calc_guide_lfc(normd$log2cpm, normd$sample_cols, normd$gene_map,
                              out_path = tbl(sprintf("guideLFC_%s.csv", screen)))
  gene_lfc  <- calc_gene_lfc(guide_lfc,
                             out_path = tbl(sprintf("geneLFC_%s.csv", screen)))
  gi        <- calc_gi(gene_lfc,
                       out_path = tbl(sprintf("GI_%s.csv", screen)))

  message(sprintf("[02_compute_GI] %s: %d guides -> %d gene pairs -> %d dual GI rows",
                  screen, nrow(guide_lfc), nrow(gene_lfc), nrow(gi)))

  ref_path <- derived(sprintf("geneLFC_%s.csv", screen))
  if (file.exists(ref_path)) {
    ref <- readr::read_csv(ref_path, show_col_types = FALSE, progress = FALSE)
    j <- dplyr::inner_join(
      dplyr::select(gene_lfc, gene,
                    new_WT  = mean_top_LFC_WT_vs_Cell_corrected,
                    new_Dep = mean_top_LFC_Depleted_vs_Cell_corrected),
      dplyr::select(ref, gene,
                    ref_WT  = mean_top_LFC_WT_vs_Cell_corrected,
                    ref_Dep = mean_top_LFC_Depleted_vs_Cell_corrected),
      by = "gene")
    max_d <- max(abs(c(j$new_WT - j$ref_WT, j$new_Dep - j$ref_Dep)), na.rm = TRUE)
    message(sprintf("[02_compute_GI] %s geneLFC vs derived ref: %d/%d pairs matched, max |Δ| = %.2e",
                    screen, nrow(j), nrow(gene_lfc), max_d))
    if (max_d > 1e-3) {
      message(sprintf("    NOTE: max |Δ| exceeds 1e-3 for %s -- inspect before trusting downstream stats.",
                      screen))
    }
  } else {
    message(sprintf("[02_compute_GI] %s: no derived reference to validate against.",
                    screen))
  }

  invisible(gi)
}

invisible(lapply(SCREENS, run_screen))

if (!interactive()) log_session("02_compute_GI")
