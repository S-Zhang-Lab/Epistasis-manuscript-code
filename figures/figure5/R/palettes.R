HOST_COLORS <- c(
  Cx3cr1 = "#C51B7D",
  Mrc1   = "#4393C3",
  WT     = "#7F8C8D"
)

DIV_COLORS <- grDevices::colorRampPalette(
  c("#2166ac", "#4393c3", "#92c5de", "#d1e5f0",
    "#f7f7f7",
    "#fddbc7", "#f4a582", "#d6604d", "#b2182b")
)(100)

DIV_RAMP_BREAKS <- c(-2, -1, 0, 1, 2)
DIV_RAMP_COLORS <- c("#2166ac", "#92c5de", "#f7f7f7", "#f4a582", "#b2182b")

NETWORK_TIER_COLORS <- c(
  "Q4: Top 25% (Super Hub)"   = "#800026",
  "Q3: 50-75% (Hub)"          = "#E31A1C",
  "Q2: 25-50% (Peripheral)"   = "#FD8D3C",
  "Q1: 0-25% (Isolated)"      = "#FECC5C"
)

QUADRANT_COLORS <- c(
  `WT(+) Dep(+)` = "#1b7837",
  `WT(+) Dep(-)` = "#fdae61",
  `WT(-) Dep(+)` = "#2166ac",
  `WT(-) Dep(-)` = "#7b3294"
)

SWITCH_COLORS <- c(
  "Cx: Lethal | Mrc: Antagonistic" = "#E41A1C",
  "Cx: Antagonistic | Mrc: Lethal" = "#377EB8"
)

DIV_RB_COLORS <- c(
  "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "white",
  "#FDDBC7", "#F4A582", "#D6604D", "#B2182B"
)

DIV_PB_COLORS <- c(
  "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "white",
  "#FDE0EF", "#F1B6DA", "#DE77AE", "#C51B7D"
)

DIV_ANCHOR_PROPS <- c(-1, -0.6, -0.25, -0.05, 0, 0.05, 0.25, 0.6, 1)

CELLTYPE_COLORS <- c(
  "Stronger in Mrc1"              = "#4393C3",
  "Stronger in Cx3cr1"            = "#C51B7D",
  "Sign-divergent"             = "#E08214",
  "Conserved across cell type" = "#7B3294",
  "Not niche-conditional"      = "gray75"
)
