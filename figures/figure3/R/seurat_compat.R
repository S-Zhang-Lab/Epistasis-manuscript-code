safe_variable_features <- function(obj, n = 2000L, assay = "SCT") {
  vf <- tryCatch(VariableFeatures(obj, assay = assay), error = function(e) character(0))
  if (length(vf) > 0) return(vf)
  mat <- safe_get_assay_data(obj, assay = assay, slot_or_layer = "data")
  vars <- Matrix::rowMeans((mat - Matrix::rowMeans(mat))^2)
  names(sort(vars, decreasing = TRUE))[seq_len(min(n, length(vars)))]
}

safe_get_assay_data <- function(obj, assay, slot_or_layer = "data") {
  if (utils::packageVersion("Seurat") >= "5.0.0") {
    GetAssayData(obj, assay = assay, layer = slot_or_layer)
  } else {
    GetAssayData(obj, assay = assay, slot = slot_or_layer)
  }
}
