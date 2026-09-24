# Generate manuscript panels for Figure 2.
steps <- c("00_check_input_data.R",
  "01_prepare_r1_gi_tables.R",
  "02_r1_modules.R",
  "04_r2_derive_gi.R",
  "03_main_panels.R",
  "05_fig2F_gi_z_histogram.R",
  "06_fig2I_pten_ego.R",
  "supp_01_plasmid_library_qc.R",
  "supp_02_day0_library_qc.R",
  "supp_04_enrichment_summary.R",
  "supp_08_bipartite_gi.R",
  "07_string_vs_gi.R",
  "supp_09_r1_plasmid_qc.R",
  "supp_10_r1_cell_qc.R")
stopifnot(file.exists(".here"), dir.exists("scripts"))
here::i_am("scripts/run.R")
env <- new.env(parent = globalenv())
for (step in steps) {
  message("Running ", step)
  status <- system2(file.path(R.home("bin"), "Rscript"), shQuote(file.path("scripts", step)))
  if (status != 0L) stop("Failed: ", step, " (exit ", status, ")")
}
message("Figure 2 panel workflow completed.")
