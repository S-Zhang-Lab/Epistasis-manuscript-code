suppressPackageStartupMessages({
  source(here::here("R", "paths.R"))
  source(here::here("R", "theme_pub.R"))
  source(here::here("R", "palettes.R"))
  source(here::here("R", "helpers.R"))
  for (pkg in c("ggplot2", "dplyr", "readr", "tidyr")) library(pkg, character.only = TRUE)
})

ANCHOR <- "Pten"

d <- readr::read_csv(tbl("ddGI_combined.csv"), show_col_types = FALSE) |>
  dplyr::filter(geneA == ANCHOR | geneB == ANCHOR) |>
  dplyr::mutate(
    partner = ifelse(geneA == ANCHOR, geneB, geneA),

    cat = dplyr::case_when(
      !genuine                        ~ "No interaction in either host",
      category == "Sign-divergent"    ~ "GI reverses between hosts",
      category == "Stronger in Mrc1"  ~ "Stronger in MRC1+ depleted host",
      category == "Stronger in Cx3cr1" ~ "Stronger in CX3CR1+ depleted host",
      TRUE                            ~ "Not significant (FDR > 0.05)"),
    sig = ddGI_qvalue < 0.05) |>
  dplyr::arrange(ddGI) |>
  dplyr::mutate(partner = factor(partner, levels = partner))

PAL <- c("Stronger in MRC1+ depleted host"   = HOST_COLORS[["Mrc1"]],
         "Stronger in CX3CR1+ depleted host" = HOST_COLORS[["Cx3cr1"]],
         "GI reverses between hosts"         = "#E08214",
         "Not significant (FDR > 0.05)"      = "gray72",
         "No interaction in either host"     = "gray88")

pts <- d |>
  dplyr::select(partner, cat, sig, dGI_Cx, dGI_Mrc) |>
  tidyr::pivot_longer(c(dGI_Cx, dGI_Mrc), names_to = "host", values_to = "dGI") |>
  dplyr::mutate(host = ifelse(host == "dGI_Cx", "Cx3cr1", "Mrc1"))

TX <- 4.55
lab <- d |>
  dplyr::mutate(stat = sprintf("%+.2f%s", ddGI, ifelse(sig, "  *", "")))

p <- ggplot() +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray45", linewidth = 0.35) +

  geom_segment(data = d,
               aes(x = dGI_Cx, xend = dGI_Mrc, y = partner, yend = partner, color = cat),
               linewidth = 0.85, lineend = "round") +

  geom_point(data = pts, aes(dGI, partner, fill = host),
             shape = 21, size = 1.05, stroke = 0.22, color = "gray20") +
  geom_text(data = lab, aes(x = TX, y = partner, label = stat,
                            fontface = ifelse(sig, "bold", "plain")),
            hjust = 1, size = 1.45, color = "gray15") +
  annotate("text", x = TX, y = 18.1, label = "Delta*Delta*\"GI\"", parse = TRUE,
           hjust = 1, size = 1.45, fontface = "bold", color = "gray30") +
  scale_colour_manual(values = PAL, name = NULL,
                      breaks = c("Stronger in MRC1+ depleted host",
                                 "GI reverses between hosts",
                                 "Not significant (FDR > 0.05)",
                                 "No interaction in either host"),

                      labels = parse(text = c(
                        '"Stronger in MRC1"^"+"*" host"',
                        '"Sign reverses"',
                        '"Not significant"',
                        '"No interaction"'))) +
  scale_fill_manual(values = c(Cx3cr1 = HOST_COLORS[["Cx3cr1"]], Mrc1 = HOST_COLORS[["Mrc1"]]),
                    name = NULL,
                    labels = parse(text = c('Delta*"GI, CX3CR1"^"+"*" dep."',
                                            'Delta*"GI, MRC1"^"+"*" dep."'))) +
  guides(color = guide_legend(order = 2, ncol = 2, byrow = TRUE,
                               override.aes = list(linewidth = 1.1)),
         fill   = guide_legend(order = 1, ncol = 2, byrow = TRUE,
                               override.aes = list(size = 1.5))) +
  coord_cartesian(xlim = c(-3.8, 3.2), ylim = c(0.5, 18.4), clip = "off") +
  labs(title = "Pten's host dependence is set by the partner",
       x = expression(Delta*"GI   (depleted - wild-type host)"),
       y = "Pten partner") +
  theme_pub() +
  theme(axis.text.y        = element_text(size = 4.0, face = "italic"),
        axis.text.x        = element_text(size = 4.0),
        axis.title         = element_text(size = 4.4),
        axis.ticks         = element_line(linewidth = 0.2),
        plot.title         = element_text(size = 4.8),
        plot.title.position = "plot",
        panel.grid.major.y = element_line(color = "gray93", linewidth = 0.18),
        legend.key.height  = unit(1.9, "mm"),
        legend.key.width   = unit(2.8, "mm"),
        legend.text        = element_text(size = 4.0),
        legend.position    = "bottom",
        legend.box         = "vertical",
        legend.box.just    = "left",
        legend.margin      = margin(0, 0, 0, 0),
        legend.box.spacing = unit(1, "mm"),
        legend.spacing.y   = unit(0.4, "mm"),
        plot.margin        = margin(1.5, 5, 1.5, 1.5, "mm"))

save_pdf(p, fig("Fig5_Pten_partner_profile.pdf"), 45.87, 52.92, snap = FALSE)

readr::write_csv(
  d[, c("partner", "dGI_Cx", "q_dGI_Cx", "dGI_Mrc", "q_dGI_Mrc",
        "ddGI", "ddGI_ci_lo", "ddGI_ci_hi", "ddGI_qvalue", "genuine", "category")],
  tbl("Pten_partner_profile.csv"))

message(sprintf("[Fig5_Pten_partner_profile] %d partners | %d genuine | %d with ddGI q<0.05: %s",
                nrow(d), sum(d$genuine), sum(d$sig),
                paste(d$partner[d$sig], collapse = ", ")))
log_session("Fig5_Pten_partner_profile")
