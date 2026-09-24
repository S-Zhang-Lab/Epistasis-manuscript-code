# Generate manuscript panels for Figure 3.
steps <- c("00_check_input_data.R",
  "00_setup.R",
  "cache/10_ucell_scores.R",
  "cache/11_aucell_hallmark.R",
  "cache/12_ucell_hallmark.R",
  "cache/13_cluster_identity.R",
  "panels/Fig3C_condition_pca.R",
  "panels/Fig3D_S3F_niche_program.R",
  "panels/Fig3E_guide_within_cluster.R",
  "panels/Fig3F_edist_dumbbell9_stats.R",
  "panels/Fig3F_edist_dumbbell9_plot.R",
  "panels/Fig3G_umap_split.R",
  "panels/Fig3H_permouse_residual.R",
  "panels/Fig3J_rank_concordance.R",
  "panels/Fig3K_seed_trajectories.R",
  "panels/FigS3C_cutting_efficiency.R",
  "panels/FigS3D_global_umap.R",
  "panels/FigS3E_cluster_dotplot.R",
  "panels/FigS3G_coefficient_flip_stats.R",
  "panels/FigS3G_coefficient_flip_plot.R",
  "panels/FigS3H_edist_bars.R",
  "panels/FigS3I_dko_density_cache.R",
  "panels/FigS3I_dko_density_plot.R",
  "panels/FigS3J_permouse_pointrange.R")
if ("--light" %in% commandArgs(trailingOnly = TRUE)) steps <- c("00_setup.R", "panels/Fig3J_rank_concordance.R", "panels/Fig3K_seed_trajectories.R", "panels/FigS3C_cutting_efficiency.R")
stopifnot(file.exists(".here"), dir.exists("scripts"))
here::i_am("scripts/run.R")
BiocParallel::register(BiocParallel::SerialParam())
env <- new.env(parent = globalenv())
for (step in steps) {
  message("Running ", step)
  source(file.path("scripts", step), local = env)
}
message("Figure 3 panel workflow completed.")
