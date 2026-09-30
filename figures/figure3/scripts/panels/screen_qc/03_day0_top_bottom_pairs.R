# =====================================================================
# Day 0 (in vivo) gene-pair representation - ranked "enrichment" view
# i.e. which pairs are most / least abundant in the starting pool.
# Day0 = sum of the 3 Cell replicates; reported as reads and % of pool.
# =====================================================================

if (!exists("save_panel", mode = "function")) source(here::here("scripts", "00_setup.R"))
need <- c("readr","dplyr","ggplot2","scales","forcats","patchwork")
invisible(lapply(need, library, character.only = TRUE))

TOPN <- 10

raw <- read_csv(derived("screens","2nd_DC_counts.csv"), show_col_types = FALSE)
DAY0_COLS <- c("Cell_1","Cell_2","Cell_3")

pair <- raw |>
  transmute(gene, d0 = rowSums(across(all_of(DAY0_COLS), as.numeric), na.rm = TRUE)) |>
  group_by(gene) |> summarise(reads = sum(d0), .groups = "drop") |>
  mutate(pct = 100 * reads / sum(reads)) |>
  arrange(desc(reads)) |>
  mutate(rank = row_number())

write_csv(pair, tbl("invivo_day0_ranked_pairs.csv"))
cat(sprintf("Total pairs=%d | total Day0 reads=%s | median=%.0f\n",
            nrow(pair), comma(sum(pair$reads)), median(pair$reads)))

# ---- top & bottom for plotting ----
top  <- pair |> slice_head(n = TOPN) |> mutate(grp = "Top 10 (most abundant)")
bot  <- pair |> slice_tail(n = TOPN) |> mutate(grp = "Bottom 10 (least abundant)")
show <- bind_rows(top, bot) |>
  mutate(grp = factor(grp, levels = c("Top 10 (most abundant)",
                                      "Bottom 10 (least abundant)")),
         gene = fct_reorder(gene, reads))

p <- ggplot(show, aes(reads, gene, fill = grp)) +
  geom_col(width = 0.72) +
  geom_text(aes(label = comma(reads)), hjust = -0.1, size = 2.7) +
  facet_wrap(~ grp, scales = "free", ncol = 1) +
  scale_x_continuous(labels = comma, expand = expansion(mult = c(0, 0.18))) +
  scale_fill_manual(values = c("Top 10 (most abundant)" = "#CC3311",
                               "Bottom 10 (least abundant)" = "#3B6FB6"),
                    guide = "none") +
  labs(x = "Day 0 reads (sum of 3 Cell replicates)", y = NULL,
       title = "Day 0 gene-pair representation (in vivo starting pool)",
       subtitle = sprintf("324 pairs; median %.0f reads; most-abundant pair = %.1fx the median",
                          median(pair$reads), max(pair$reads)/median(pair$reads))) +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank())

ggsave(fig_supp("invivo_day0_top_pairs.pdf"), p, width = 8.5, height = 9)
ggsave(fig_supp("invivo_day0_top_pairs.png"), p, width = 8.5, height = 9, dpi = 300)
message("Saved invivo_day0_top_pairs.pdf / .png and invivo_day0_ranked_pairs.csv")

cat("\n=== Top 15 Day0 pairs ===\n");  print(head(pair, 15),  n = 15)
cat("\n=== Bottom 5 Day0 pairs ===\n"); print(tail(pair, 5),  n = 5)