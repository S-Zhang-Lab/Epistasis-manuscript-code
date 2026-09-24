# Generate manuscript panels for Figure 5.
steps <- c("00_check_input_data.R",
  "02_compute_GI.R",
  "06_stats_GI.R",
  "07_combine_ddGI.R",
  "09_filtered_objective_analysis.R",
  "Fig5_first_order_scatter.R",
  "Fig5_first_order_bubbles.R",
  "Fig5_additivity.R",
  "Fig5_ddGI_map.R",
  "Fig5_pervasiveness_waterfall.R",
  "Fig5_PtenCdh1_hero.R",
  "Fig5_Pten_partner_profile.R",
  "Fig5_invivo_burden.R",
  "Fig5_pooled_burden_specificity.R",
  "10_panel_A_pca.R",
  "Fig5S_S1_wt_filter_qc.R",
  "Fig5S_S11_scale_sensitivity.R",
  "Fig5S_S6_first_order_lollipop.R",
  "Fig5S_S2_sign_distribution.R",
  "Fig5S_S7_stability.R",
  "Fig5S_S12_PtenCdh1_detail.R")
stopifnot(file.exists(".here"), dir.exists("scripts"))
here::i_am("scripts/run.R")
env <- new.env(parent = globalenv())
for (step in steps) {
  message("Running ", step)
  status <- system2(file.path(R.home("bin"), "Rscript"), shQuote(file.path("scripts", step)))
  if (status != 0L) stop("Failed: ", step, " (exit ", status, ")")
}
message("Figure 5 panel workflow completed.")
