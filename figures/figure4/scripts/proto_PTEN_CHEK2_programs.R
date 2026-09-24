suppressPackageStartupMessages({
  source(here::here("R", "paths.R"))
  source(here::here("R", "theme_pub.R"))
  for (pkg in c("ggplot2", "dplyr", "readr", "ggrepel")) {
    if (!requireNamespace(pkg, quietly = TRUE))
      install.packages(pkg, repos = "https://cloud.r-project.org")
    library(pkg, character.only = TRUE)
  }
})
SEED <- 123L

rd <- function(f) {
  d <- readr::read_csv(f, show_col_types = FALSE, progress = FALSE)
  names(d)[1] <- "Pathway"
  dplyr::transmute(d,
    key = sub("^HALLMARK_", "", gsub("-", "_", toupper(as.character(Pathway)))),
    GI  = as.numeric(GI))
}

pten <- rd(derived("pathwayGI_PTEN_CX3CL1.csv"))  |> dplyr::rename(PTEN = GI)
chek <- rd(derived("pathwayGI_CHEK2_CX3CL1.csv")) |> dplyr::rename(CHEK2 = GI)

m <- dplyr::inner_join(pten, chek, by = "key") |>
  dplyr::mutate(short = gsub("_", " ", key), divergence = PTEN - CHEK2)
r <- cor(m$PTEN, m$CHEK2)

t_lean <- 0.5 * stats::sd(m$divergence)
m <- m |> dplyr::mutate(leaning = dplyr::case_when(
  divergence >  t_lean ~ "Pten-specific (host-interfacing)",
  divergence < -t_lean ~ "Chek2-specific (cell-autonomous)",
  TRUE                 ~ "shared"))

cat(sprintf("\n[PTEN/CHEK2] %d pathways; cor(PTEN, CHEK2) = %.2f; t_lean = %.4f\n",
            nrow(m), r, t_lean))
cat("  Pten-specific:", sum(m$leaning == "Pten-specific (host-interfacing)"),
    "| Chek2-specific:", sum(m$leaning == "Chek2-specific (cell-autonomous)"),
    "| shared:", sum(m$leaning == "shared"), "\n")
cat("  Pten-specific programs:\n    ",
    paste(gsub("_"," ", m$key[m$leaning == "Pten-specific (host-interfacing)"]), collapse = " | "), "\n")
cat("  Chek2-specific programs:\n    ",
    paste(gsub("_"," ", m$key[m$leaning == "Chek2-specific (cell-autonomous)"]), collapse = " | "), "\n")

lim <- max(abs(c(m$PTEN, m$CHEK2))) * 1.05
lab <- m |> dplyr::slice_max(abs(divergence), n = 10)
pad <- 0.18 * lim

PTEN_FILL <- "#D6604D"; PTEN_LAB <- "#B2182B"
CHEK_FILL <- "#7570B3"; CHEK_LAB <- "#4B468F"

p <- ggplot(m, aes(PTEN, CHEK2)) +

  annotate("rect", xmin = pad, xmax = lim, ymin = -lim, ymax = -pad,
           fill = PTEN_FILL, alpha = 0.07) +

  annotate("rect", xmin = -lim, xmax = -pad, ymin = pad, ymax = lim,
           fill = CHEK_FILL, alpha = 0.07) +
  geom_abline(slope = 1, intercept = 0, color = "grey55", linewidth = 0.4) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70", linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey70", linewidth = 0.3) +
  geom_point(aes(fill = leaning), shape = 21, size = 2.1, color = "white",
             stroke = 0.25, alpha = 0.9) +
  ggrepel::geom_text_repel(data = lab, aes(label = short, color = leaning),
             size = 2.0, max.overlaps = 25, segment.size = 0.2,
             min.segment.length = 0, box.padding = 0.3, seed = SEED,
             show.legend = FALSE) +
  scale_fill_manual(values = c("Pten-specific (host-interfacing)" = PTEN_FILL,
                               "Chek2-specific (cell-autonomous)" = CHEK_FILL,
                               "shared" = "grey70"), name = NULL) +
  scale_color_manual(values = c("Pten-specific (host-interfacing)" = PTEN_LAB,
                                "Chek2-specific (cell-autonomous)" = CHEK_LAB,
                                "shared" = "grey50")) +
  coord_equal(xlim = c(-lim, lim), ylim = c(-lim, lim)) +
  labs(
    title = "The niche-enabled pair engages host-interfacing programs",
    subtitle = sprintf("within-lane raw-AUC pathway-GI (NT*NT-anchored) · Pten (UTSW37 Synergy) vs Chek2 (UTSW33) · r = %.2f", r),
    x = expression("Pten_Cx3cl1  pathway GI ("*Delta*"AUC)"),
    y = expression("Chek2_Cx3cl1  pathway GI ("*Delta*"AUC)")) +
  theme_pub(base_size = 7) +
  theme(legend.position = "right",
        plot.subtitle = element_text(size = 5.2, colour = "grey25"),
        legend.key.size = unit(7, "pt"),
        legend.text = element_text(size = 5.6))

save_pdf(p, fig("Fig4H_PTEN_vs_CHEK2_programs.pdf"), 121, 110)
message("  wrote Fig4H_PTEN_vs_CHEK2_programs.pdf\n")
if (!interactive()) tryCatch(log_session("proto_PTEN_CHEK2_programs"), error = function(e) NULL)
