suppressPackageStartupMessages({
  source(here::here("R", "paths.R"))
  source(here::here("R", "theme_pub.R"))
  source(here::here("R", "palettes.R"))
  source(here::here("R", "helpers.R"))
  for (pkg in c("ggplot2", "dplyr", "readr")) {
    if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
})

comb <- readr::read_csv(tbl("ddGI_combined.csv"), show_col_types = FALSE, progress = FALSE)
gen  <- comb |> dplyr::filter(genuine) |> dplyr::arrange(ddGI) |>
  dplyr::mutate(rank = dplyr::row_number())

PAL <- c("Stronger in Mrc1" = "#4393C3", "Stronger in Cx3cr1" = "#C51B7D",
         "Sign-divergent" = "#E08214", "Conserved across cell type" = "gray78")
cnt  <- gen |> dplyr::count(category)
getn <- function(x) { v <- cnt$n[cnt$category == x]; if (length(v)) v else 0L }
labs_leg <- c(
  "Stronger in Mrc1"           = sprintf("Stronger in Mrc1⁺ cell depleted host (n=%d)",   getn("Stronger in Mrc1")),
  "Stronger in Cx3cr1"         = sprintf("Stronger in Cx3cr1⁺ cell depleted host (n=%d)", getn("Stronger in Cx3cr1")),
  "Sign-divergent"             = sprintf("GI reverses between niches (n=%d)", getn("Sign-divergent")),
  "Conserved across cell type" = sprintf("Not significant (FDR > 0.05, n=%d)", getn("Conserved across cell type")))
gen   <- gen |> dplyr::mutate(category = factor(category, levels = names(PAL)))
ntot  <- nrow(gen); nspec <- sum(gen$category != "Conserved across cell type"); pct <- round(100 * nspec / ntot)
hero  <- dplyr::filter(gen, gene == "Pten*Cdh1")

p <- ggplot(gen, aes(rank, ddGI, fill = category)) +
  geom_col(width = 1) +
  geom_hline(yintercept = 0, color = "gray40", linewidth = 0.3) +
  geom_point(data = hero, aes(rank, ddGI - 0.5), shape = 23, size = 2.6, fill = "#FFD400",
             color = "black", stroke = 0.5, inherit.aes = FALSE) +
  annotate("text", x = hero$rank + 4, y = hero$ddGI - 1.05, label = "Pten×Cdh1",
           size = 2.1, fontface = "bold", hjust = 0) +
  annotate("text", x = ntot * 0.74, y = -2.6,
           label = sprintf("%d of %d genuine interactions\n(%d%%) differ between hosts\nΔΔGI significant, BH-FDR q<0.05", nspec, ntot, pct),
           size = 2.35, fontface = "bold", color = "#1a3a6b", hjust = 0.5, lineheight = 1.08) +
  scale_fill_manual(values = PAL, labels = labs_leg, name = NULL,
                    breaks = c("Stronger in Mrc1", "Stronger in Cx3cr1", "Sign-divergent", "Conserved across cell type")) +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.01))) +
  coord_cartesian(ylim = c(-4.9, 4.9), clip = "off") +
  labs(title = "Niche-specific epistasis is pervasive",
       subtitle = "every genuine interaction (real GI in ≥ 1 host) ranked by ΔΔGI; colored = FDR-significant difference between hosts",
       x = "gene pairs (ranked by ΔΔGI)",
       y = expression(Delta*Delta*"GI  ("*Delta*"GI"[Cx3cr1]*" − "*Delta*"GI"[Mrc1]*")")) +
  theme_pub() +
  theme(legend.position = c(0.015, 0.99), legend.justification = c(0, 1),
        legend.text = element_text(size = 5), legend.key.size = unit(7, "pt"),
        legend.background = element_rect(fill = alpha("white", 0.85), color = "gray85", linewidth = 0.3),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        plot.subtitle = element_text(size = 5.5, color = "gray30"))

save_pdf(p, fig("Fig5_pervasiveness_waterfall.pdf"), 121, 84)
message(sprintf("[Fig5_pervasiveness_waterfall] %d/%d genuine = %d%% cell-type-specific (BH q<0.05)", nspec, ntot, pct))
log_session("Fig5_pervasiveness_waterfall")
