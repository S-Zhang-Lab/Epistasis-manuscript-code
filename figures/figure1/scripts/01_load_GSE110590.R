message("[01_load_GSE110590] Loading expression matrix + paired subset metadata...")

exprs_file <- file.path(INPUT_DIR, "GSE110590",
                        "GSE110590_RAP_A16_log2.sne.tsv.gz")
exprs_data <- read.delim(gzfile(exprs_file), check.names = FALSE)

gene_symbols <- vapply(strsplit(as.character(exprs_data$NAME), "\\|"),
  function(x) {
    if (length(x) == 1)    return(x)
    if (x[1] == "?")       return(paste0("ENTREZ_", x[2]))
    x[1]
  }, character(1))

exprs_mat <- as.matrix(exprs_data[, -1])
rownames(exprs_mat) <- make.unique(gene_symbols)

normalize_sample_name <- function(x) sub("-RNA-RNA$", "-RNA", as.character(x))
colnames(exprs_mat) <- normalize_sample_name(colnames(exprs_mat))

custom_meta <- read.csv(
  file.path(DATA_DERIVED, "Sample_Metadata_Table.csv"),
  check.names = FALSE, stringsAsFactors = FALSE)
custom_meta$sample_id  <- trimws(custom_meta$sample_id)
custom_meta$patient_id <- trimws(custom_meta$patient_id)

exprs_sub <- exprs_mat[, custom_meta$sample_id, drop = FALSE]
exprs_sub <- apply(exprs_sub, 2, as.numeric)
rownames(exprs_sub) <- rownames(exprs_mat)
exprs_sub <- exprs_sub[complete.cases(exprs_sub), ]
exprs_sub <- exprs_sub[apply(exprs_sub, 1, var) > 1e-6, ]

group_sub   <- factor(ifelse(custom_meta$tissue_type == "primary",
                             "Primary", "LungMet"),
                      levels = c("Primary", "LungMet"))
patient_sub <- factor(custom_meta$patient_id)
design_sub  <- model.matrix(~ group_sub)
