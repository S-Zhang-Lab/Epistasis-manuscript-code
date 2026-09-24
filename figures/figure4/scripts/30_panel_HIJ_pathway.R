suppressPackageStartupMessages({
  source(here::here("R", "paths.R"))
  source(here::here("R", "theme_pub.R"))
  source(here::here("R", "palettes.R"))
  source(here::here("R", "helpers.R"))
  for (pkg in c("ggplot2", "dplyr", "tidyr", "tibble", "readr", "readxl",
                "scales", "ggrepel", "igraph", "uwot")) {
    if (!requireNamespace(pkg, quietly = TRUE))
      install.packages(pkg, repos = "https://cloud.r-project.org")
    library(pkg, character.only = TRUE)
  }
})

SEED            <- 123L

LAYOUT_SEED     <- 39L
KEEP_QUANTILE   <- 0.80
LABEL_THRESHOLD <- 1.5
LABEL_NODES     <- FALSE
FR_NITER        <- 2000L
PULL            <- 0.20
PULL_ITERS      <- 2L

CLUSTER_SEP     <- suppressWarnings(as.numeric(Sys.getenv("FIG5_CLUSTER_SEP", "1.5")))
if (is.na(CLUSTER_SEP)) CLUSTER_SEP <- 1.0

COMMUNITY_LABELS <- c(
  "1" = "Lipid & oxidative metabolism",
  "2" = "Immune & inflammation",
  "3" = "Hormone response",
  "4" = "EMT & developmental",
  "5" = "Apical surface",
  "6" = "Glycolysis/mTOR",
  "7" = "Cell cycle/DNA repair",
  "8" = "KRAS signaling (down)",
  "9" = "Pancreas beta-cells")
COMMUNITY_LABELS_CANVAS <- c(
  "1" = "Metabolism",
  "2" = "Immune",
  "3" = "ER hormone",
  "4" = "EMT &\ndevelopment",
  "6" = "Glycolysis",
  "7" = "Cell cycle/DNA repair")

CLUSTER_LABEL_PUSH <- 0.55

set.seed(SEED)

message("\n=== 30_panel_HIJ_pathway ===\n")

dash_to_uscore <- function(x) gsub("-", "_", x, fixed = TRUE)
uscore_to_dash <- function(x) gsub("_", "-", x, fixed = TRUE)

build_similarity_graph <- function(S, keep_quantile = 0.8, weighted = TRUE) {
  diag(S) <- 0
  thr <- stats::quantile(S[S > 0], keep_quantile, na.rm = TRUE)
  A   <- if (weighted) S * (S >= thr) else (S >= thr) * 1
  igraph::graph_from_adjacency_matrix(A, mode = "undirected",
                                      weighted = weighted, diag = FALSE)
}

compute_layout_comm_cohesion <- function(g, pull = 0.20, pull_iters = 2L,
                                         niter_fr = 2000L, seed = 123L) {
  set.seed(seed)
  wts <- if ("weight" %in% igraph::edge_attr_names(g)) igraph::E(g)$weight else NULL
  coords <- igraph::layout_with_fr(g, weights = wts, grid = "nogrid",
                                   niter = niter_fr)
  rownames(coords) <- igraph::V(g)$name
  if (!"community" %in% igraph::vertex_attr_names(g)) return(coords)
  comm <- as.factor(igraph::V(g)$community)
  for (it in seq_len(pull_iters)) {
    for (k in levels(comm)) {
      idx <- which(comm == k)
      if (length(idx) > 1L) {
        cx <- mean(coords[idx, 1L]); cy <- mean(coords[idx, 2L])
        coords[idx, 1L] <- coords[idx, 1L] * (1 - pull) + cx * pull
        coords[idx, 2L] <- coords[idx, 2L] * (1 - pull) + cy * pull
      }
    }
  }
  coords
}

explode_communities <- function(coords, community, sep = 1.0) {
  if (!is.finite(sep) || sep == 1) return(coords)
  gc   <- colMeans(coords)
  comm <- as.factor(community)
  for (k in levels(comm)) {
    idx   <- which(comm == k)
    cc    <- colMeans(coords[idx, , drop = FALSE])
    shift <- (cc - gc) * (sep - 1)
    coords[idx, 1L] <- coords[idx, 1L] + shift[1L]
    coords[idx, 2L] <- coords[idx, 2L] + shift[2L]
  }
  coords
}

build_net_plot <- function(g, gi_attr, title_main, title_sub, coords,
                           label_threshold = LABEL_THRESHOLD,
                           ring_expand = 1.21, ring_stroke = 1.0,
                           node_stroke = 0.15, label_nodes = LABEL_NODES,
                           color_cap = NULL) {

  comm_levels  <- levels(as.factor(igraph::V(g)$community))
  comm_palette <- setNames(scales::hue_pal()(length(comm_levels)), comm_levels)
  gi_vals <- igraph::vertex_attr(g, gi_attr)

  n_lab <- if (label_nodes) 12L else 0L
  .ft <- sort(abs(gi_vals), decreasing = TRUE)
  label_threshold <- if (n_lab > 0L) .ft[min(n_lab, length(.ft))] else Inf

  nodes <- tibble::tibble(
    name      = igraph::V(g)$name,
    x         = coords[igraph::V(g)$name, 1L],
    y         = coords[igraph::V(g)$name, 2L],
    community = factor(igraph::V(g)$community, levels = comm_levels),
    GI        = gi_vals,
    abs_GI    = abs(gi_vals),
    label     = ifelse(abs(gi_vals) >= label_threshold,
                       paste0(uscore_to_dash(name), " (", signif(gi_vals, 2), ")"),
                       "")
  )
  e <- igraph::as_data_frame(g, what = "edges")
  edges <- e %>%
    dplyr::transmute(
      x1 = coords[from, 1L], y1 = coords[from, 2L],
      x2 = coords[to,   1L], y2 = coords[to,   2L],
      weight = if (!is.null(weight)) weight else 1
    )
  lim <- if (!is.null(color_cap)) color_cap else max(nodes$abs_GI, na.rm = TRUE)
  if (!is.finite(lim) || lim == 0) lim <- 1

  gx <- mean(nodes$x); gy <- mean(nodes$y)
  cl <- nodes %>%
    dplyr::group_by(community) %>%
    dplyr::summarise(cx = mean(x), cy = mean(y), .groups = "drop") %>%
    dplyr::mutate(
      label = unname(COMMUNITY_LABELS_CANVAS[as.character(community)]),
      lx = cx + (cx - gx) * CLUSTER_LABEL_PUSH,
      ly = cy + (cy - gy) * CLUSTER_LABEL_PUSH) %>%
    dplyr::filter(!is.na(label) & nzchar(label)) %>%
    dplyr::mutate(lab_colour = unname(comm_palette[as.character(community)]))

  ggplot() +
    geom_segment(data = edges,
                 aes(x = x1, y = y1, xend = x2, yend = y2, alpha = weight),
                 colour = "grey60", linewidth = 0.12, show.legend = FALSE) +
    geom_point(data = nodes,
               aes(x = x, y = y, size = abs_GI * ring_expand, colour = GI),
               shape = 21, fill = NA, stroke = ring_stroke, alpha = 0.9) +
    geom_point(data = nodes,
               aes(x = x, y = y, size = abs_GI, fill = community),
               shape = 21, colour = "grey20", stroke = node_stroke, alpha = 0.95) +
    ggrepel::geom_text_repel(
      data = dplyr::filter(nodes, abs_GI >= label_threshold),
      aes(x = x, y = y, label = label),
      size = 1.35, max.overlaps = 30, segment.size = 0.1,
      point.padding = 0.1, seed = SEED) +
    geom_text(data = cl, aes(x = lx, y = ly, label = label),
              colour = cl$lab_colour, fontface = "bold.italic",
              size = 1.65, lineheight = 0.9, show.legend = FALSE,
              inherit.aes = FALSE) +
    scale_size_continuous(name = "| GI |", range = c(0.5, 2.3),
      guide = guide_legend(keyheight = unit(2.2, "mm"), order = 2)) +
    scale_fill_manual(name = "Pathway cluster", values = comm_palette,
      labels = function(b) { lab <- COMMUNITY_LABELS[b]; ifelse(is.na(lab), b, lab) },
      guide = guide_legend(override.aes = list(shape = 21, size = 1.5, alpha = 1),
                           keyheight = unit(1.75, "mm"), order = 3)) +
    scale_colour_gradient2(name = "GI (ring)",
      low = "#2166ac", mid = "#f7f7f7", high = "#b2182b",
      midpoint = 0, limits = c(-lim, lim), oob = scales::squish,
      breaks = c(-1, 0, 1) * signif(lim * 0.8, 1),
      guide = guide_colourbar(barheight = unit(5.5, "mm"), barwidth = unit(1.4, "mm"),
                              ticks.linewidth = 0.3, order = 1)) +
    theme_void(base_size = 5) +
    theme(plot.margin      = margin(1.5, 1.5, 1.5, 1.5, "mm"),
          plot.title       = element_text(size = 5.4, hjust = 0),
          plot.subtitle    = element_text(size = 4.2, colour = "grey30", hjust = 0),
          legend.title     = element_text(size = 4.6, face = "bold"),
          legend.text      = element_text(size = 4.2),
          legend.key.size  = unit(1.9, "mm"),
          legend.spacing.y = unit(0.05, "mm"),
          legend.box.spacing = unit(0.6, "mm"),
          legend.margin    = margin(0, 0, 0, 0.5, "mm")) +
    coord_equal(clip = "off") +
    labs(title = title_main, subtitle = title_sub)
}

read_gi <- function(file, pair_name) {
  f <- derived(file)
  if (!file.exists(f)) stop("Missing pathway input: ", f, call. = FALSE)
  df <- readr::read_csv(f, show_col_types = FALSE, progress = FALSE)
  names(df)[1L] <- "Pathway"
  df %>% dplyr::transmute(Pathway = as.character(Pathway),
                          GI = as.numeric(GI), pair = pair_name)
}

gi_long_all <- read_gi("pathwayGI_PTEN_CX3CL1.csv", "PTEN_CX3CL1")

unique_pw  <- unique(gi_long_all$Pathway)
has_prefix <- grepl("^HALLMARK", unique_pw)
pathways_us <- ifelse(has_prefix, dash_to_uscore(unique_pw),
                      paste0("HALLMARK_", dash_to_uscore(unique_pw)))
rownames_map <- tibble::tibble(Pathway_dash = unique_pw,
                               Pathway_us   = pathways_us)

scaffold_nodes_f <- derived("network_scaffold_Cx3cr1_nodes.csv")
scaffold_edges_f <- derived("network_scaffold_Cx3cr1_edges.csv")

if (file.exists(scaffold_nodes_f) && file.exists(scaffold_edges_f)) {
  message("[30_panel_HIJ] FROZEN scaffold: ", basename(scaffold_nodes_f),
          " (no msigdbr call)")
  nd <- readr::read_csv(scaffold_nodes_f, show_col_types = FALSE, progress = FALSE)
  ed <- readr::read_csv(scaffold_edges_f, show_col_types = FALSE, progress = FALSE)
  g  <- igraph::graph_from_data_frame(ed[, c("from", "to", "weight")],
          directed = FALSE, vertices = nd[, "name", drop = FALSE])
  igraph::V(g)$community <- factor(nd$community[match(igraph::V(g)$name, nd$name)])
  coords <- as.matrix(nd[, c("x", "y")]); rownames(coords) <- nd$name

  .relayout_seed <- suppressWarnings(as.integer(Sys.getenv("FIG5_RELAYOUT_SEED", "")))
  if (!is.na(.relayout_seed)) {
    message("[30_panel_HIJ] RE-LAYOUT on frozen graph, FR seed ", .relayout_seed)
    coords <- compute_layout_comm_cohesion(g, pull = PULL, pull_iters = PULL_ITERS,
                                           niter_fr = FR_NITER, seed = .relayout_seed)
  }
  message(sprintf("[30_panel_HIJ] Network (frozen): %d nodes, %d edges",
                  igraph::vcount(g), igraph::ecount(g)))
} else {
  if (!requireNamespace("msigdbr", quietly = TRUE))
    install.packages("msigdbr", repos = "https://cloud.r-project.org")
  message("[30_panel_HIJ] LIVE scaffold from msigdbr ",
          as.character(utils::packageVersion("msigdbr")), " (freezing to CSV)")
  hallmark <- tryCatch(
    msigdbr::msigdbr(species = "Homo sapiens", collection = "H"),
    error = function(e) msigdbr::msigdbr(species = "Homo sapiens", category = "H"))
  msig_list <- split(hallmark$gene_symbol, hallmark$gs_name)
  msig_list <- msig_list[names(msig_list) %in% pathways_us]
  if (length(msig_list) == 0L)
    stop("No pathways matched msigdbr -- check Hallmark naming.", call. = FALSE)
  message(sprintf("[30_panel_HIJ] %d Hallmark pathways matched (of %d in inputs)",
                  length(msig_list), length(unique_pw)))

  pw <- names(msig_list); n <- length(pw)
  S  <- matrix(0, n, n, dimnames = list(pw, pw))
  for (i in seq_len(n)) {
    S[i, i] <- 1
    if (i < n) for (j in (i + 1L):n) {
      jc <- length(intersect(msig_list[[i]], msig_list[[j]])) /
            length(union    (msig_list[[i]], msig_list[[j]]))
      S[i, j] <- S[j, i] <- jc
    }
  }
  g <- build_similarity_graph(S, keep_quantile = KEEP_QUANTILE, weighted = TRUE)
  clu <- igraph::cluster_louvain(g)
  igraph::V(g)$community <- factor(igraph::membership(clu))
  message(sprintf("[30_panel_HIJ] Network: %d nodes, %d edges; Louvain %d (mod %.3f)",
                  igraph::vcount(g), igraph::ecount(g),
                  length(unique(igraph::V(g)$community)), igraph::modularity(clu)))
  coords <- compute_layout_comm_cohesion(g, pull = PULL, pull_iters = PULL_ITERS,
                                         niter_fr = FR_NITER, seed = LAYOUT_SEED)

  nd <- tibble::tibble(
    name      = igraph::V(g)$name,
    community = as.integer(as.character(igraph::V(g)$community)),
    x = coords[igraph::V(g)$name, 1L], y = coords[igraph::V(g)$name, 2L])
  ed <- igraph::as_data_frame(g, what = "edges")
  if (is.null(ed$weight)) ed$weight <- 1
  readr::write_csv(nd, scaffold_nodes_f)
  readr::write_csv(ed[, c("from", "to", "weight")], scaffold_edges_f)
  message("[30_panel_HIJ] froze scaffold -> ", basename(scaffold_nodes_f),
          " + ", basename(scaffold_edges_f))
}

attach_gi <- function(g, gi_df, attr_name) {
  df_us <- gi_df %>% dplyr::left_join(rownames_map,
                                      by = c("Pathway" = "Pathway_dash"))
  vals  <- df_us$GI[match(names(igraph::V(g)), df_us$Pathway_us)]
  vals[is.na(vals)] <- 0
  igraph::set_vertex_attr(g, attr_name, value = vals)
}
for (pair in unique(gi_long_all$pair)) {
  safe <- gsub("[*]", "_", pair)
  g <- attach_gi(g, gi_long_all %>% dplyr::filter(pair == !!pair), safe)
}

coords <- explode_communities(coords, igraph::V(g)$community, CLUSTER_SEP)
if (CLUSTER_SEP != 1)
  message(sprintf("[30_panel_HIJ] cluster separation x%.2f", CLUSTER_SEP))

comm_tbl <- tibble::tibble(
  Pathway   = igraph::V(g)$name,
  community = as.integer(as.character(igraph::V(g)$community))
)
readr::write_csv(comm_tbl, tbl("pathway_network_communities.csv"))
message(sprintf("    Wrote %s", tbl("pathway_network_communities.csv")))

PANEL_W_PT <- 190; PANEL_H_PT <- 120
panel_w <- PANEL_W_PT * 25.4 / 72
panel_h <- PANEL_H_PT * 25.4 / 72

PTEN_COLOR_CAP <- 0.005
if ("PTEN_CX3CL1" %in% igraph::vertex_attr_names(g)) {
  p_G <- build_net_plot(g, "PTEN_CX3CL1",
    title_main = "Pten_Cx3cl1 epistatic transcriptional manifold",
    title_sub  = "Hallmark pathway-GI (Perturb-seq AUCell)",
    coords     = coords, color_cap = PTEN_COLOR_CAP)
  save_pdf(p_G, fig("Fig4G_PTEN_CX3CL1_network.pdf"),
           panel_w, panel_h, snap_width = FALSE)
}

if (!interactive()) log_session("30_panel_HIJ_pathway")
invisible(g)
