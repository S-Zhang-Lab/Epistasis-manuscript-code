suppressPackageStartupMessages({ library(data.table) })

.args <- commandArgs(trailingOnly = FALSE)
.fa   <- sub("^--file=", "", .args[grepl("^--file=", .args)])
REPO_ROOT <- if (length(.fa) > 0) normalizePath(file.path(dirname(.fa), "..")) else normalizePath(getwd())
VERIFY <- "--verify" %in% commandArgs(trailingOnly = TRUE)

RAW <- file.path(REPO_ROOT, "data", "raw", "msk_met_2021_raw")
OUT <- if (VERIFY) file.path(tempdir(), "msk_met_verify") else
                   file.path(REPO_ROOT, "data", "raw", "MSK_MET_2021_breast")
SHIPPED <- file.path(REPO_ROOT, "data", "raw", "MSK_MET_2021_breast")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

SCHEMA <- list(
  clinical  = c("sampleId","CANCER_TYPE","CANCER_TYPE_DETAILED","SUBTYPE",
                "SUBTYPE_ABBREVIATION","SAMPLE_TYPE","METASTATIC_SITE","MET_SITE_COUNT",
                "PRIMARY_SITE","DMETS_DX_BONE","DMETS_DX_LIVER","DMETS_DX_LUNG",
                "DMETS_DX_CNS_BRAIN","DMETS_DX_DIST_LN","DMETS_DX_PLEURA","DMETS_DX_SKIN"),
  mutations = c("sampleId","hugo","entrez","mutationType","proteinChange"),
  cna       = c("sampleId","hugo","entrez","alteration"))

need <- function(f) {
  p <- file.path(RAW, f)
  if (!file.exists(p))
    stop("Missing raw input: ", p, "\n",
         "  Fetch the upstream study first, from the pinned cBioPortal DataHub\n",
         "  revision recorded in DATA.md:\n",
         "    mkdir -p data/raw/msk_met_2021_raw\n",
         "    for f in data_clinical_sample.txt data_clinical_patient.txt \\\n",
         "             data_mutations.txt data_cna.txt; do curl -fsSL -o \\\n",
         "      data/raw/msk_met_2021_raw/$f \\\n",
         "      https://media.githubusercontent.com/media/cBioPortal/datahub/",
         "db2f8008a119f6008fba5a99102d085a863995a1/public/msk_met_2021/$f; done\n",
         "  (not normally needed: the derived CSVs ship with the repository)",
         call. = FALSE)
  p
}

read_clin <- function(f) fread(need(f), sep = "\t", header = TRUE, na.strings = c("", "NA"),
                              skip = "PATIENT_ID", showProgress = FALSE)

cat("[etl] reading clinical tables ...\n")
smp <- read_clin("data_clinical_sample.txt")
pat <- read_clin("data_clinical_patient.txt")
setnames(smp, toupper(names(smp))); setnames(pat, toupper(names(pat)))

shared <- setdiff(intersect(names(smp), names(pat)), "PATIENT_ID")
if (length(shared)) pat[, (shared) := NULL]
clin <- merge(smp, pat, by = "PATIENT_ID", all.x = TRUE)

cat(sprintf("[etl] %d specimens total; filtering to breast ...\n", nrow(clin)))
clin <- clin[CANCER_TYPE == "Breast Cancer"]
if (!nrow(clin)) stop("No rows with CANCER_TYPE == 'Breast Cancer'; check the study release.", call. = FALSE)

clin[, sampleId := SAMPLE_ID]
missing_cols <- setdiff(SCHEMA$clinical, names(clin))
if (length(missing_cols))
  stop("Clinical schema mismatch. Missing column(s): ", paste(missing_cols, collapse = ", "),
       "\n  The upstream study layout may have changed; update SCHEMA in this script.", call. = FALSE)
clin_cols <- SCHEMA$clinical
clinical  <- clin[, ..clin_cols]
setcolorder(clinical, SCHEMA$clinical)
setorder(clinical, sampleId)
breast_ids <- clinical$sampleId
cat(sprintf("[etl] breast specimens: %d\n", length(breast_ids)))

cat("[etl] reading mutations ...\n")
mut <- fread(need("data_mutations.txt"), sep = "\t", header = TRUE,
             na.strings = c("", "NA"), showProgress = FALSE)
mutations <- mut[Tumor_Sample_Barcode %in% breast_ids,
                 .(sampleId = Tumor_Sample_Barcode, hugo = Hugo_Symbol,
                   entrez = Entrez_Gene_Id, mutationType = Variant_Classification,
                   proteinChange = HGVSp_Short)]
setorder(mutations, sampleId, hugo)
cat(sprintf("[etl] mutation rows (breast): %d\n", nrow(mutations)))

cat("[etl] reading and melting copy-number matrix ...\n")
cna_wide <- fread(need("data_cna.txt"), sep = "\t", header = TRUE, showProgress = FALSE)
id_cols  <- intersect(c("Hugo_Symbol","Entrez_Gene_Id"), names(cna_wide))
keep     <- intersect(breast_ids, setdiff(names(cna_wide), id_cols))
cna_long <- melt(cna_wide[, c(id_cols, keep), with = FALSE],
                 id.vars = id_cols, variable.name = "sampleId",
                 value.name = "alteration", variable.factor = FALSE)
setnames(cna_long, c("Hugo_Symbol","Entrez_Gene_Id"), c("hugo","entrez"), skip_absent = TRUE)
cna <- cna_long[!is.na(alteration), .(sampleId, hugo, entrez, alteration = as.integer(alteration))]
setorder(cna, sampleId, hugo)
cat(sprintf("[etl] copy-number rows (breast): %d\n", nrow(cna)))

for (nm in c("clinical","mutations","cna")) {
  obj <- get(nm)
  stopifnot(identical(names(obj), SCHEMA[[nm]]))
  fwrite(obj, file.path(OUT, paste0(nm, ".csv")))
}
cat(sprintf("[etl] wrote 3 CSVs to %s\n", OUT))

if (VERIFY) {
  cat("\n[etl --verify] comparing against the shipped CSVs\n")
  ok <- TRUE
  for (nm in c("clinical","mutations","cna")) {
    a <- file.path(OUT, paste0(nm, ".csv")); b <- file.path(SHIPPED, paste0(nm, ".csv"))
    if (!file.exists(b)) { cat(sprintf("  %-10s shipped copy absent, skipped\n", nm)); next }
    ha <- tools::md5sum(a); hb <- tools::md5sum(b)
    same <- identical(unname(ha), unname(hb))
    if (!same) {
      da <- fread(a); db <- fread(b)
      cat(sprintf("  %-10s DIFFERS  (derived %d rows vs shipped %d rows)\n", nm, nrow(da), nrow(db)))
      ok <- FALSE
    } else cat(sprintf("  %-10s identical\n", nm))
  }
  cat(if (ok) "[etl --verify] PASS\n" else "[etl --verify] FAIL: derived CSVs differ from the shipped copies\n")
  if (!ok) quit(status = 1)
}
