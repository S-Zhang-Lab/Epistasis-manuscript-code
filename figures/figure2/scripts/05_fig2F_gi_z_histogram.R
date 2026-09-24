pkgs <- c("ggplot2")
inst <- rownames(installed.packages())
if (length(setdiff(pkgs, inst)))
  install.packages(setdiff(pkgs, inst), repos = "https://cloud.r-project.org")
suppressPackageStartupMessages(library(ggplot2))

source(here::here("R", "paths.R"))
source(here::here("R", "helpers.R"))
source(here::here("R", "theme_pub.R"))

EDGE_CUT <- 1.5
COL_BUF  <- "#2166AC"
COL_SYN  <- "#B2182B"

gi <- read.csv(tbl("R2_GI_derived.csv"), check.names = FALSE)
z  <- gi$GI_Z_Score[gi$GI_Z_Include == 1]
stopifnot(length(z) == 252L)

n_buf <- sum(z < -EDGE_CUT); n_syn <- sum(z > EDGE_CUT)
rng   <- range(z)
message(sprintf("[Fig2F] n=%d  range %.2f..%.2f | |Z|>%.1f: %d buffering, %d synergistic",
                length(z), rng[1], rng[2], EDGE_CUT, n_buf, n_syn))

p <- ggplot(data.frame(z = z), aes(z)) +
  annotate("rect", xmin = -Inf, xmax = -EDGE_CUT, ymin = 0, ymax = Inf,
           fill = COL_BUF, alpha = 0.08) +
  annotate("rect", xmin = EDGE_CUT, xmax = Inf, ymin = 0, ymax = Inf,
           fill = COL_SYN, alpha = 0.08) +
  geom_histogram(binwidth = 0.25, boundary = 0,
                 fill = "#3F6E63", colour = "white", linewidth = 0.18) +
  geom_vline(xintercept = 0, colour = "grey35", linewidth = LW_THIN) +

  annotate("text", x = -Inf, y = Inf, vjust = 1.3, hjust = -0.08,
           size = pt(FS_LABEL), colour = COL_BUF, lineheight = 0.95,
           label = sprintf("Buffering\n%d pairs", n_buf)) +
  annotate("text", x = Inf, y = Inf, vjust = 1.3, hjust = 1.08,
           size = pt(FS_LABEL), colour = COL_SYN, lineheight = 0.95,
           label = sprintf("Synergistic\n%d pairs", n_syn)) +
  scale_x_continuous(breaks = seq(-2, 6, 2),
                     expand = expansion(mult = c(0.02, 0.02))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(x = "GI Z-score", y = "Number of gene pairs") +
  theme_panel()

save_panel(p, fig("Fig2F_GI_Z_histogram.pdf"))
log_session("05_fig2F_gi_z_histogram")
message("[Fig2F] Saved Fig2F_GI_Z_histogram.{pdf,png} for placement inside panel F")
