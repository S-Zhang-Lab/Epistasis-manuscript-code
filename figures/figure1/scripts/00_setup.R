if (!requireNamespace("here", quietly = TRUE)) {
  install.packages("here", repos = "https://cloud.r-project.org")
}

source(here::here("R", "paths.R"))
source(here::here("R", "palettes.R"))
source(here::here("R", "devices.R"))
source(here::here("R", "theme_pub.R"))
source(here::here("R", "helpers.R"))

suppressPackageStartupMessages({

  library(maftools)
  library(circlize)

  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(ggrepel)

  library(ComplexHeatmap)
  library(grid)

  library(limma)
  library(clusterProfiler)
  library(enrichplot)
  library(msigdbr)
  library(patchwork)

  library(readxl)
})

set.seed(42)

source(here::here("scripts", "utils_variants.R"))

message("[00_setup] OK — repo root: ", REPO_ROOT)
message("[00_setup] R: ", R.version.string)
