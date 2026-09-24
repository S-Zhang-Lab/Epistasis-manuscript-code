suppressPackageStartupMessages({
  source(here::here("R", "paths.R"))
  source(here::here("R", "theme_pub.R"))
  source(here::here("R", "palettes.R"))
  source(here::here("R", "helpers.R"))
  library(ggplot2); library(dplyr)
})

POOLDIR <- here::here("data", "Pooled_Burden_Quantification")
mrc <- read.csv(file.path(POOLDIR, "Mrc1_B6_DTR_nodule_counts.csv"))
cx  <- read.csv(file.path(POOLDIR, "Cx3cr1_B6_DTR_nodule_counts.csv"))
cl  <- read.csv(tbl("tumor_burden_foci.csv"))

p_exact <- function(a, b) {
  n1 <- length(a); n2 <- length(b); v <- c(a, b); Umax <- n1 * n2
  U  <- function(x, y) sum(outer(x, y, ">") + 0.5 * outer(x, y, "=="))
  T0 <- { o <- U(a, b); min(o, Umax - o) }
  idx <- utils::combn(n1 + n2, n1)
  mean(apply(idx, 2, function(i) { u <- U(v[i], v[-i]); min(u, Umax - u) <= T0 + 1e-9 }))
}
SEED <- 123L; N_BOOT <- 20000L
set.seed(SEED)
est <- function(dep, wt) {
  r <- replicate(N_BOOT, mean(sample(dep, length(dep), TRUE)) / mean(sample(wt, length(wt), TRUE)))
  data.frame(fold = mean(dep) / mean(wt),
             lo = unname(quantile(r, .025)), hi = unname(quantile(r, .975)),
             n_wt = length(wt), n_dep = length(dep), p = p_exact(dep, wt))
}

d <- bind_rows(
  cbind(key = "cx_pool",  est(cx$nodule_count[cx$group == "Cx3cr1_DTR"], cx$nodule_count[cx$group == "Cx3cr1_B6"])),
  cbind(key = "mrc_pool", est(mrc$nodule_count[mrc$group == "Mrc1_DTR"], mrc$nodule_count[mrc$group == "Mrc1_B6"])),
  cbind(key = "clone",    est(cl$foci[cl$host == "Mrc1-DTR"],            cl$foci[cl$host == "Mrc1-B6 (WT)"]))
) |>
  mutate(grp  = ifelse(key == "clone", "clone", "pool"),
         stat = sprintf("%.2f  (%.2f-%.2f)   P = %s", fold, lo, hi,
                        ifelse(p < 0.01, sprintf("%.4f", p), sprintf("%.2f", p))),
         key  = factor(key, levels = c("clone", "mrc_pool", "cx_pool")))

norm_pts <- function(key, dep, wt) {
  m <- mean(wt)
  rbind(data.frame(key = key, host = "wild-type", val = wt  / m),
        data.frame(key = key, host = "depleted",  val = dep / m))
}
pts <- rbind(
  norm_pts("cx_pool",  cx$nodule_count[cx$group == "Cx3cr1_DTR"], cx$nodule_count[cx$group == "Cx3cr1_B6"]),
  norm_pts("mrc_pool", mrc$nodule_count[mrc$group == "Mrc1_DTR"], mrc$nodule_count[mrc$group == "Mrc1_B6"]),
  norm_pts("clone",    cl$foci[cl$host == "Mrc1-DTR"],            cl$foci[cl$host == "Mrc1-B6 (WT)"])
)
pts$key  <- factor(pts$key, levels = c("clone", "mrc_pool", "cx_pool"))
pts$host <- factor(pts$host, levels = c("wild-type", "depleted"))
pts$y   <- as.numeric(pts$key) + ifelse(pts$host == "depleted", 0.21, -0.21)

plotlab <- c(
  cx_pool  = '"Pooled library,  Cx3cr1"^"+"*" depleted"',
  mrc_pool = '"Pooled library,  Mrc1"^"+"*" depleted"',
  clone    = '"Pten" %*% "Cdh1 clone,  Mrc1"^"+"*" depleted"')

TX <- 2^5.85
p <- ggplot(d, aes(fold, key)) +
  annotate("rect", xmin = 0.26, xmax = 40, ymin = 0.58, ymax = 1.42, fill = "#FFD400", alpha = 0.13) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray45", linewidth = 0.35) +
  geom_point(data = pts, aes(x = val, y = y, fill = host), shape = 21,
             size = 1.15, stroke = 0.28, color = "gray40", alpha = 0.9,
             inherit.aes = FALSE) +

  geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y", width = 0,
                color = "black", linewidth = 0.6) +
  geom_point(aes(shape = "est"), fill = "black", color = "black",
             size = 2.2, stroke = 0.4) +
  scale_shape_manual(values = c(est = 22), name = NULL,
                     labels = "Ratio of means, 95% CI") +
  scale_fill_manual(values = c("wild-type" = "white", "depleted" = "gray55"),
                    name = NULL, labels = c("Wild-type host", "Depleted host")) +
  guides(shape = guide_legend(order = 1),
         fill  = guide_legend(order = 2,
                              override.aes = list(size = 2.2, stroke = 0.35))) +
  geom_text(aes(x = TX, label = stat), hjust = 0, size = 2.05, color = "gray15") +
  annotate("text", x = TX, y = 3.70, label = "fold change (95% CI)", hjust = 0,
           size = 2.0, fontface = "bold", color = "gray30") +
  scale_y_discrete(labels = function(x) parse(text = plotlab[x])) +
  scale_x_continuous(trans = "log2", breaks = c(0.5, 1, 2, 4, 8, 16, 32),
                     labels = c("0.5", "1", "2", "4", "8", "16", "32")) +
  coord_cartesian(xlim = c(0.33, 32), ylim = c(0.62, 3.55), clip = "off") +
  labs(title    = expression("The Mrc1-dependent burden increase is specific to Pten" %*% "Cdh1"),
       subtitle = "Ratio of means, 95% bootstrap CI, exact permutation test. Each small dot is one\nmouse, divided by its own assay's wild-type mean.\nn = 7 vs 7, 5 vs 3, 6 vs 6, top to bottom.",
       x = "fold change in lung metastatic burden", y = NULL) +
  theme_pub() +
  theme(axis.text.y      = element_text(size = 6.0, hjust = 1),
        plot.title       = element_text(size = 7.2),
        plot.subtitle    = element_text(size = 4.8, color = "gray30", lineheight = 1.15),
        panel.grid.major.y = element_blank(),
        legend.position  = "bottom",
        legend.box       = "horizontal",
        legend.spacing.x = unit(2.5, "mm"),
        legend.justification = "left",
        legend.margin    = margin(-1, 0, 0, 0),
        legend.key.height = unit(3.2, "mm"),
        plot.margin      = margin(3, 40, 3, 3, "mm"))

save_pdf(p, fig("Fig5_pooled_burden_specificity.pdf"), 121, 62)

raw <- bind_rows(
  transform(mrc, screen = "Mrc1",   host = ifelse(group == "Mrc1_DTR",   "depleted", "wild-type")),
  transform(cx,  screen = "Cx3cr1", host = ifelse(group == "Cx3cr1_DTR", "depleted", "wild-type"))
)[, c("screen", "host", "panel", "nodule_count")]
readr::write_csv(raw, tbl("pooled_burden_nodule_counts.csv"))
readr::write_csv(d[, c("key", "fold", "lo", "hi", "n_wt", "n_dep", "p")],
                 tbl("burden_specificity_summary.csv"))

message(sprintf("[Fig5_pooled_burden_specificity] pooled Mrc1 %.2fx (P=%.2f) | pooled Cx3cr1 %.2fx (P=%.2f) | PtenxCdh1 %.2fx (P=%.4f)",
                d$fold[2], d$p[2], d$fold[1], d$p[1], d$fold[3], d$p[3]))
log_session("Fig5_pooled_burden_specificity")
