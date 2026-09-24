# Generate manuscript panels for Figure 1.
steps <- c("00_check_input_data.R",
  "00_setup.R",
  "01_load_GSE110590.R",
  "Fig1_B_pancancer_chord.R",
  "Fig1_msk_met_panels.R",
  "Fig1_F_PCA.R",
  "Fig1_G_pathway_rewiring.R",
  "Fig1_H0_paired_DEG_modules.R",
  "Fig1_H_module_DEG_bar.R",
  "Fig1_I_niche_engagement.R",
  "Fig1_J_intrinsic_niche_coactivation.R",
  "supplementary/FigS1_cohort_and_exposure.R",
  "supplementary/FigS1_divergence.R")
stopifnot(file.exists(".here"), dir.exists("scripts"))
here::i_am("scripts/run.R")
env <- new.env(parent = globalenv())
for (step in steps) {
  message("Running ", step)
  source(file.path("scripts", step), local = env)
}
message("Figure 1 panel workflow completed.")
