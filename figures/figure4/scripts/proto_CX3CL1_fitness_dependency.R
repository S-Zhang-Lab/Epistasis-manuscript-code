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

gi <- readr::read_csv(tbl("GI_Cx3cr1.csv"), show_col_types = FALSE,
                      progress = FALSE) |>
  dplyr::mutate(dGI = GI_Depleted_vs_Cell - GI_WT_vs_Cell)

cx <- gi |>
  dplyr::filter(geneA == "Cx3cl1" | geneB == "Cx3cl1") |>
  dplyr::mutate(
    driver = ifelse(geneA == "Cx3cl1", geneB, geneA),
    host_cond = (GI_WT_qvalue < 0.05 | GI_Depleted_qvalue < 0.05) &
                delta_GI_qvalue < 0.05,
    perturbseq = driver %in% c("Pten", "Apc", "Chek2"),
    lab = host_cond | perturbseq)

cat(sprintf("\n[cx3cl1] %d driver x Cx3cl1 pairs; host-conditional (q<.05): %d\n",
            nrow(cx), sum(cx$host_cond)))
cx |> dplyr::arrange(dGI) |>
  dplyr::transmute(driver, dGI = round(dGI, 2), q = round(delta_GI_qvalue, 3),
                   host_cond) |> as.data.frame() |> print(row.names = FALSE)

lim <- max(abs(cx$dGI)) * 1.15
p <- ggplot(cx, aes(dGI, reorder(driver, dGI))) +
  geom_vline(xintercept = 0, color = "grey55", linewidth = 0.4) +
  geom_segment(aes(x = 0, xend = dGI, yend = driver,
                   color = host_cond), linewidth = 0.6) +
  geom_point(aes(color = host_cond, size = abs(dGI))) +

  geom_point(data = dplyr::filter(cx, perturbseq), aes(dGI, driver),
             shape = 5, size = 3.8, stroke = 0.7, color = "grey15") +
  scale_color_manual(values = c(`TRUE` = "#B2182B", `FALSE` = "grey65"),
                     labels = c(`TRUE` = "host-conditional",
                                `FALSE` = "not host-conditional"),
                     name = NULL) +
  scale_size_area(max_size = 3.2, guide = "none") +
  scale_x_continuous(limits = c(-lim, lim)) +
  labs(
    title = "Host-dependency of CX3CL1 interactions",
    subtitle = "8/17 host-conditional; diamond = Perturb-seq",
    x = expression("host-"*Delta*"GI  (depleted - WT)"), y = NULL) +
  theme_pub(base_size = 7) +
  theme(legend.position = "top",
        plot.subtitle = element_text(size = 5.8, colour = "grey25"),
        axis.text.y = element_text(size = 6.5),
        legend.key.size = unit(7, "pt"))

save_pdf(p, fig("Fig4F_CX3CL1_fitness_dependency.pdf"), 78, 95, snap_width = FALSE)
message("  wrote Fig4F_CX3CL1_fitness_dependency.pdf\n")

if (!interactive()) tryCatch(log_session("proto_CX3CL1_fitness_dependency"), error = function(e) NULL)
