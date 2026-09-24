# Generate manuscript panels for Figure 4.
steps <- c("00_check_input_data.R",
  "01_compute_GI.R",
  "05_stats.R",
  "10_panel_DEF_heatmaps.R",
  "proto_panelA_hostconditional.R",
  "proto_panelB_driverpathway.R",
  "proto_CX3CL1_fitness_dependency.R",
  "30_panel_HIJ_pathway.R",
  "proto_PTEN_CHEK2_programs.R",
  "proto_A2slope_Bstructure.R",
  "proto_CHEK2_estrogen_cascade.R")
stopifnot(file.exists(".here"), dir.exists("scripts"))
here::i_am("scripts/run.R")
env <- new.env(parent = globalenv())
for (step in steps) {
  message("Running ", step)
  status <- system2(file.path(R.home("bin"), "Rscript"), shQuote(file.path("scripts", step)))
  if (status != 0L) stop("Failed: ", step, " (exit ", status, ")")
}
message("Figure 4 panel workflow completed.")
