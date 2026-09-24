suppressPackageStartupMessages({
  source(here::here("R", "paths.R"))
  source(here::here("R", "theme_pub.R"))
  source(here::here("R", "palettes.R"))
  source(here::here("R", "helpers.R"))
  for (pkg in c("ggplot2", "dplyr", "readr")) library(pkg, character.only = TRUE)
})

a <- readr::read_csv(tbl("niche_anchor_analysis.csv"), show_col_types = FALSE) |>
  dplyr::filter(layer == "first_order") |>
  dplyr::mutate(dir = dplyr::case_when(q < 0.05 & mean > 0 ~ "MRC1-biased (q<0.05)",
                                       q < 0.05 & mean < 0 ~ "Cx3cr1-biased (q<0.05)",
                                       TRUE                ~ "n.s."))

ntot <- nrow(a)
ord  <- dplyr::arrange(a, dplyr::desc(mean))
a    <- dplyr::distinct(dplyr::bind_rows(head(ord, 14), tail(ord, 2)), gene, .keep_all = TRUE)
nsig <- sum(ord$q < 0.05)

cdh1_x <- a$mean[a$gene == "Cdh1"]
COL <- c("MRC1-biased (q<0.05)" = HOST_COLORS[["Mrc1"]],
         "Cx3cr1-biased (q<0.05)" = HOST_COLORS[["Cx3cr1"]], "n.s." = "gray78")

p <- ggplot(a, aes(mean, reorder(gene, mean), color = dir)) +
  geom_vline(xintercept = 0, linewidth = 0.4, color = "gray45") +
  geom_segment(aes(x = 0, xend = mean, yend = reorder(gene, mean)), linewidth = 0.55) +
  geom_point(aes(size = -log10(q))) +
  scale_colour_manual(values = COL, name = NULL) +
  scale_size_continuous(range = c(1.2, 3.2), guide = "none") +
  labs(title = "First-order anchor ranking",
       subtitle = sprintf("%d most niche-biased of %d anchors shown\n%d reach q < 0.05; CDH1 top (q < 0.001)",
                          nrow(a), ntot, nsig),
       x = expression("diff. niche-dep. (MRC1 - CX3CR1)"), y = NULL) +
  theme_pub() +
  theme(axis.text.y = element_text(size = 6),
        legend.position = "top", legend.text = element_text(size = 5.4),
        plot.subtitle = element_text(size = 5.2, color = "gray25", lineheight = 1.15))

save_pdf(p, fig("Fig5S_S6_first_order_lollipop.pdf"), 57, 95)
message(sprintf("[Fig5S_S6] first-order anchors; CDH1 mean=%.2f", cdh1_x))
log_session("Fig5S_S6_first_order_lollipop")
