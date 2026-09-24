log_session <- function(slug) {
  source(here::here("R", "paths.R"), local = FALSE)
  out <- log_path(sprintf("%s_session.txt", slug))
  con <- file(out, open = "wt")
  on.exit(close(con), add = TRUE)
  cat(sprintf("# %s — %s\n\n", slug, format(Sys.time(), tz = "UTC", usetz = TRUE)),
      file = con)
  capture.output(sessionInfo(), file = con, append = TRUE)
  invisible(out)
}

load_merged_data <- function() {
  source(here::here("R", "paths.R"), local = FALSE)
  if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
  if (!requireNamespace("readr", quietly = TRUE)) install.packages("readr")

  cx3cr1 <- readr::read_csv(derived("geneLFC_Cx3cr1.csv"), show_col_types = FALSE)
  mrc1   <- readr::read_csv(derived("geneLFC_Mrc1.csv"),   show_col_types = FALSE)

  dplyr::inner_join(
    cx3cr1 |> dplyr::select(gene,
                            WT_Cx  = mean_top_LFC_WT_vs_Cell_corrected,
                            Dep_Cx = mean_top_LFC_Depleted_vs_Cell_corrected),
    mrc1   |> dplyr::select(gene,
                            WT_Mrc  = mean_top_LFC_WT_vs_Cell_corrected,
                            Dep_Mrc = mean_top_LFC_Depleted_vs_Cell_corrected),
    by = "gene"
  ) |>
    dplyr::mutate(
      dLFC_Cx  = Dep_Cx  - WT_Cx,
      dLFC_Mrc = Dep_Mrc - WT_Mrc,
      diff_niche_dep  = dLFC_Mrc - dLFC_Cx,
      abs_wt_diff  = abs(WT_Mrc - WT_Cx),
      sign_flip    = (WT_Cx * WT_Mrc) < 0
    ) |>
    dplyr::filter(!(sign_flip & abs_wt_diff > 2.0))
}

compute_synergy_table <- function() {
  source(here::here("R", "paths.R"), local = FALSE)
  if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
  if (!requireNamespace("tidyr", quietly = TRUE)) install.packages("tidyr")
  if (!requireNamespace("readr", quietly = TRUE)) install.packages("readr")

  cx3cr1 <- readr::read_csv(derived("geneLFC_Cx3cr1.csv"), show_col_types = FALSE)
  mrc1   <- readr::read_csv(derived("geneLFC_Mrc1.csv"),   show_col_types = FALSE)

  get_effect_df <- function(df, suffix) {
    df |>
      dplyr::mutate(Effect = mean_top_LFC_Depleted_vs_Cell_corrected -
                              mean_top_LFC_WT_vs_Cell_corrected) |>
      dplyr::select(gene, !!paste0("Effect_", suffix) := Effect)
  }

  all_eff <- dplyr::inner_join(
    get_effect_df(cx3cr1, "Cx"),
    get_effect_df(mrc1,   "Mrc"),
    by = "gene"
  )

  extract_singles <- function(eff_df, suffix) {
    col_name <- paste0("Effect_", suffix)
    list(
      A = eff_df |>
        dplyr::filter(grepl("\\*NT$", gene)) |>
        dplyr::mutate(gene_name = sub("\\*NT$", "", gene)) |>
        dplyr::select(gene_name, !!paste0("SingleA_", suffix) := !!rlang::sym(col_name)),
      B = eff_df |>
        dplyr::filter(grepl("^NT\\*", gene)) |>
        dplyr::mutate(gene_name = sub("^NT\\*", "", gene)) |>
        dplyr::select(gene_name, !!paste0("SingleB_", suffix) := !!rlang::sym(col_name))
    )
  }

  cx_s  <- extract_singles(all_eff, "Cx")
  mrc_s <- extract_singles(all_eff, "Mrc")

  all_eff |>
    dplyr::filter(!grepl("NT", gene)) |>
    tidyr::separate(gene, into = c("GeneA", "GeneB"), sep = "\\*", remove = FALSE) |>
    dplyr::left_join(cx_s$A,  by = c("GeneA" = "gene_name")) |>
    dplyr::left_join(cx_s$B,  by = c("GeneB" = "gene_name")) |>
    dplyr::left_join(mrc_s$A, by = c("GeneA" = "gene_name")) |>
    dplyr::left_join(mrc_s$B, by = c("GeneB" = "gene_name")) |>
    dplyr::mutate(
      Expected_Cx   = SingleA_Cx + SingleB_Cx,
      Synergy_Cx    = Effect_Cx  - Expected_Cx,
      Expected_Mrc  = SingleA_Mrc + SingleB_Mrc,
      Synergy_Mrc   = Effect_Mrc - Expected_Mrc,
      Delta_Synergy = Synergy_Cx - Synergy_Mrc
    ) |>
    dplyr::filter(!is.na(Delta_Synergy))
}
