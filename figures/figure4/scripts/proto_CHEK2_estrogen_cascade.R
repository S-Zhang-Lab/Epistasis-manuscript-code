suppressPackageStartupMessages({
  source(here::here("R", "paths.R"))
  source(here::here("R", "theme_pub.R"))
  for (p in c("ggplot2", "dplyr", "tidyr", "readr")) {
    if (!requireNamespace(p, quietly = TRUE)) install.packages(p, repos = "https://cloud.r-project.org")
    library(p, character.only = TRUE)
  }
})

GREEN <- "#1a9850"; SYN <- "#b2182b"

d <- readr::read_csv(derived("pathwayGI_CHEK2_CX3CL1.csv"), show_col_types = FALSE)
names(d)[1] <- "Pathway"
er <- d |>
  filter(grepl("ESTROGEN", Pathway)) |>
  mutate(lab = ifelse(grepl("EARLY", Pathway),
                      "Estrogen response (early)", "Estrogen response (late)"))

CASCADE <- c("Chek2_NT", "NT_Cx3cl1", "Additive\nexpectation", "Chek2_Cx3cl1\nobserved")
o2 <- er |>
  transmute(lab,
            `Chek2_NT`               = singleA  - NTNT,
            `NT_Cx3cl1`              = singleB  - NTNT,
            `Additive\nexpectation`  = Expected - NTNT,
            `Chek2_Cx3cl1\nobserved` = dual     - NTNT) |>
  pivot_longer(-lab, names_to = "step", values_to = "dev") |>
  mutate(step = factor(step, levels = CASCADE),
         kind = ifelse(grepl("expectation", step), "expected",
                ifelse(grepl("observed", step), "observed", "single")))

p_cascade <- ggplot(o2, aes(step, dev, fill = kind)) +
  geom_hline(yintercept = 0, colour = "grey35", linewidth = 0.4) +
  geom_col(width = 0.6) +
  scale_fill_manual(values = c(single = "grey72", expected = SYN, observed = GREEN),
                    breaks = c("single", "expected", "observed"),
                    labels = c("single knockout", "additive expectation", "observed double KO"),
                    name = NULL) +
  geom_text(aes(label = sprintf("%+.4f", dev),
                vjust = ifelse(dev < 0, 1.45, -0.55)), size = 1.8, colour = "grey20") +
  facet_wrap(~lab, scales = "free_y") +
  labs(title = "Each single knockout lowers the estrogen program; losing both restores it",
       subtitle = "deviation in AUCell from the NT_NT baseline; the gap between expectation and observation is the interaction",
       x = NULL, y = expression(Delta*" AUCell vs NT_NT")) +
  theme_pub(base_size = 6.5) +
  theme(axis.text.x     = element_text(size = 5.4, lineheight = 0.85),
        plot.title      = element_text(size = 7.4, face = "bold"),
        plot.subtitle   = element_text(size = 5, colour = "grey30"),
        strip.text      = element_text(size = 6.8, face = "bold"),
        legend.position = "bottom", legend.key.size = unit(6, "pt"),
        legend.text     = element_text(size = 5.4))
save_pdf(p_cascade, fig("SuppFig4_CHEK2_estrogen_cascade.pdf"), 150, 80, snap_width = FALSE)
message("  wrote SuppFig4_CHEK2_estrogen_cascade.pdf")

if (!interactive()) tryCatch(log_session("proto_CHEK2_estrogen_cascade"), error = function(e) NULL)
