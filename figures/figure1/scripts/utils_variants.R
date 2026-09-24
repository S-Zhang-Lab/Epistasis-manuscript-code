.read_variants_manifest <- function(path) {
  if (!file.exists(path)) {
    stop("selected-main variants manifest not found: ", path)
  }
  lines <- readLines(path, warn = FALSE)

  lines <- sub("\\s*#.*$", "", lines)
  lines <- trimws(lines)
  lines <- lines[nzchar(lines)]

  kv <- regmatches(lines, regexec("^([A-Za-z0-9_]+)\\s*:\\s*(\\S+)\\s*$", lines))
  ok <- vapply(kv, function(m) length(m) == 3, logical(1))
  if (!all(ok)) {
    stop("Malformed line(s) in ", path, ":\n  ",
         paste(lines[!ok], collapse = "\n  "))
  }
  out <- vapply(kv, function(m) m[[3]], character(1))
  names(out) <- vapply(kv, function(m) m[[2]], character(1))
  out
}

SELECTED_VARIANTS <- .read_variants_manifest(
  file.path(BASE_DIR, "docs", "selected_main_variants.yml")
)

read_main_variant <- function(panel) {
  panel <- toupper(panel)
  if (!(panel %in% names(SELECTED_VARIANTS))) {
    stop("Panel '", panel, "' not declared in docs/selected_main_variants.yml.\n",
         "  Declared panels: ", paste(names(SELECTED_VARIANTS), collapse = ", "))
  }
  unname(SELECTED_VARIANTS[panel])
}

panel_output_path <- function(panel, variant_id, basename, placed = FALSE) {
  selected <- read_main_variant(panel)
  dir <- if (identical(variant_id, selected)) OUTPUT_FIG_DIR
         else if (isTRUE(placed))             OUTPUT_SUP_DIR
         else                                 OUTPUT_SUP_VAR_DIR
  file.path(dir, basename)
}
