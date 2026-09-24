suppressPackageStartupMessages({
  source(here::here("R", "paths.R"))
  source(here::here("R", "helpers.R"))
  source(here::here("R", "gi_pipeline.R"))
})

SEED <- 123L
set.seed(SEED)

message("\n=== 01_compute_GI ===\n")

normd <- load_and_normalize(raw("Cx3cr1_raw_counts.csv"))

guide_lfc <- calc_guide_lfc(normd$log2cpm, normd$sample_cols, normd$gene_map,
                            out_path = tbl("guideLFC_Cx3cr1.csv"))
gene_lfc  <- calc_gene_lfc(guide_lfc, out_path = tbl("geneLFC_Cx3cr1.csv"))
gi        <- calc_gi(gene_lfc,        out_path = tbl("GI_Cx3cr1.csv"))

message(sprintf("[01_compute_GI] Summary: %d guides -> %d gene pairs -> %d dual GI rows",
                nrow(guide_lfc), nrow(gene_lfc), nrow(gi)))

if (!interactive()) log_session("01_compute_GI")

invisible(list(guide_lfc = guide_lfc, gene_lfc = gene_lfc, gi = gi))
