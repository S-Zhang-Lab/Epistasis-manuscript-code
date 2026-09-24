pkgs <- c("ggplot2", "dplyr", "tibble", "igraph", "ggrepel")
inst <- rownames(installed.packages())
if (length(setdiff(pkgs, inst)))
  install.packages(setdiff(pkgs, inst), repos = "https://cloud.r-project.org")
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tibble); library(igraph)
})

source(here::here("R", "paths.R"))
source(here::here("R", "helpers.R"))
source(here::here("R", "theme_pub.R"))
source(here::here("R", "palettes.R"))

STRING_MED <- 0.40
GI_CUT     <- 1.5
FR_SEED    <- 42L

tc <- function(x) paste0(toupper(substring(x, 1, 1)), tolower(substring(x, 2)))

tumor <- c("PTEN","APC","RB1","BRCA1","BRCA2","CHEK2","TSC2","TGFBR2","SMAD4",
           "STK11","CALML3","CLCA2","SOX10","KRT15","SERPINB5","CHIT1","KIT")
niche <- c("BTLA","CXCR2","CXCR5","FCRL6","IL6","MPEG1","PDCD1LG2","TLR7",
           "TNFSF15","CX3CL1","MFAP5","MUC16","PLA2G2A","PLA2R1","HGF","CDH1",
           "MYSM1")

mod_map <- c(
  PTEN="Growth_Control", APC="Growth_Control", RB1="Growth_Control",
  SMAD4="Growth_Control", STK11="Growth_Control", TSC2="Growth_Control",
  BRCA1="DNA_Repair", BRCA2="DNA_Repair", CHEK2="DNA_Repair",
  CALML3="Barrier_ECM", SERPINB5="Barrier_ECM",
  SOX10="Cell_Identity", KRT15="Cell_Identity", CLCA2="Cell_Identity",
  CDH1="Cell_Identity", TGFBR2="Unknown", KIT="Unknown", CHIT1="Unknown",
  BTLA="Immune_Evasion", CXCR2="Immune_Evasion", CXCR5="Immune_Evasion",
  FCRL6="Immune_Evasion", IL6="Immune_Evasion", MPEG1="Immune_Evasion",
  PDCD1LG2="Immune_Evasion", TLR7="Immune_Evasion", TNFSF15="Immune_Evasion",
  CX3CL1="Barrier_ECM", MFAP5="Barrier_ECM", MUC16="Barrier_ECM",
  PLA2G2A="Barrier_ECM", PLA2R1="Barrier_ECM",
  HGF="Growth_Signaling", MYSM1="Unknown"
)

st <- read.delim(ref("R2_STRING_network.tsv"), check.names = FALSE)
gi_strong <- read.csv(tbl("Fig2I_strong_GI_edges.csv"), stringsAsFactors = FALSE)
nd <- read.csv(tbl("Fig2I_node_strong_GI_partner_counts.csv"),
               stringsAsFactors = FALSE)
gi_all <- read.csv(tbl("R2_GI_derived.csv"), stringsAsFactors = FALSE,
                   check.names = FALSE)
is_nt <- function(x) is.na(x) | toupper(trimws(x)) %in% c("NT", "NTC", "")
tp <- gi_all[!is_nt(gi_all$Gene1) & !is_nt(gi_all$Gene2), ]

st$A <- toupper(st$preferredName_A)
st$B <- toupper(st$preferredName_B)
cross <- (st$A %in% tumor & st$B %in% niche) |
         (st$A %in% niche & st$B %in% tumor)
st_edges <- st[cross & st$score >= STRING_MED, ]

gi_edges <- gi_strong

k <- function(a, b) paste(pmin(toupper(a), toupper(b)),
                          pmax(toupper(a), toupper(b)), sep = "|")
st_keys <- k(st_edges$A, st_edges$B)
gi_keys <- k(gi_edges$gene_A, gi_edges$gene_B)
shared  <- intersect(st_keys, gi_keys)
both_n  <- length(shared)

ft <- fisher.test(table(
  k(toupper(tp$Gene1), toupper(tp$Gene2)) %in% st_keys,
  k(toupper(tp$Gene1), toupper(tp$Gene2)) %in% gi_keys
))
write.csv(tibble(
  metric = c("assayed pairs", "STRING cross-program >= 0.40", "strong GI",
             "both", "shared pairs", "Fisher OR", "Fisher P"),
  value  = c(nrow(tp), nrow(st_edges), nrow(gi_edges), both_n,
             paste(gsub("\\|", "-", shared), collapse = "; "),
             round(ft$estimate, 3), round(ft$p.value, 3))
), tbl("SuppFig_STRING_vs_GI_stats.csv"), row.names = FALSE)

all_genes <- toupper(nd$gene)
stopifnot(length(all_genes) == 34L)

g_gi <- graph_from_data_frame(
  gi_edges[, c("gene_A", "gene_B")], directed = FALSE,
  vertices = data.frame(name = all_genes)
)
set.seed(FR_SEED)
lay <- layout_with_fr(g_gi, weights = gi_edges$abs_GI_Z)
rownames(lay) <- V(g_gi)$name

nodes <- tibble(
  name = all_genes, x = lay[all_genes, 1], y = lay[all_genes, 2],
  lab  = tc(all_genes),
  module = factor(mod_map[all_genes], levels = names(R2_MODULE_COLORS))
)

st_ed <- st_edges |>
  left_join(nodes |> select(name, x, y) |> rename(xa = x, ya = y),
            by = c("A" = "name")) |>
  left_join(nodes |> select(name, x, y) |> rename(xb = x, yb = y),
            by = c("B" = "name"))

pA <- ggplot() +
  geom_segment(data = st_ed, aes(xa, ya, xend = xb, yend = yb),
               colour = "grey45", linewidth = 0.4, alpha = 0.55) +
  geom_point(data = nodes, aes(x, y, fill = module),
             shape = 21, colour = "white", size = 2.0, stroke = 0.3) +
  scale_fill_manual(values = R2_MODULE_COLORS, name = "Functional\nmodule") +
  ggrepel::geom_text_repel(
    data = nodes, aes(x, y, label = lab),
    size = pt(FS_LABEL), seed = FR_SEED, max.overlaps = Inf,
    box.padding = 0.16, segment.size = 0.12, segment.alpha = 0.4,
    min.segment.length = 0.2
  ) +
  coord_equal(clip = "off") +
  labs(title = sprintf("STRING network (%d edges)", nrow(st_ed))) +
  theme_void(base_size = MIN_FONT_PT) +
  theme(plot.title    = element_text(size = FS_TITLE, hjust = 0),
        legend.title  = element_text(size = FS_LEGEND),
        legend.text   = element_text(size = FS_LEGEND),
        legend.key.size = unit(5, "pt"))

save_panel(pA, supp_fig("SuppFig2P_STRING_network.pdf"))
log_session("07_string_vs_gi")
message(sprintf(
  "[STRING vs GI] STRING %d edges | GI %d edges | shared %d (%s) | Fisher P = %.2f",
  nrow(st_ed), nrow(gi_edges), both_n,
  paste(gsub("\\|", "-", shared), collapse = ", "), ft$p.value))
