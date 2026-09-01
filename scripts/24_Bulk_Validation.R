# ============================================================
# bulk Validation
# ============================================================

library(dplyr)
library(data.table)


# ------------------------------------------------------------
# 0. EDIT THESE PATHS to match where you saved the downloaded files
# ------------------------------------------------------------

expr_path     <- "data/metabric/metabric/data_mrna_illumina_microarray.txt"
clinical_path <- "data/metabric/metabric/data_clinical_patient.txt"
sample_path   <- "data/metabric/metabric/data_clinical_sample.txt"   # optional, for subtype


# ------------------------------------------------------------
# 1. Load the gene signature from Step 8
# ------------------------------------------------------------

signature <- read.csv("results/TME/gene_signature_FINAL.csv")$gene
cat("Signature genes loaded:", length(signature), "\n")
print(signature)


# ------------------------------------------------------------
# 2. Load bulk expression matrix
#    cBioPortal files: first column = Hugo_Symbol (gene name),
#    second column often = Entrez_Gene_Id, rest = sample columns
# ------------------------------------------------------------

expr <- fread(expr_path, sep = "\t", header = TRUE, data.table = FALSE)

cat("\nExpression matrix dimensions:", dim(expr), "\n")
cat("First few column names:\n")
print(head(colnames(expr), 10))

# Identify the gene-symbol column (usually "Hugo_Symbol")
gene_col <- colnames(expr)[1]
cat("\nUsing gene ID column:", gene_col, "\n")

# Drop non-sample columns (Hugo_Symbol, Entrez_Gene_Id if present)
non_sample_cols <- intersect(colnames(expr), c("Hugo_Symbol", "Entrez_Gene_Id"))
sample_cols <- setdiff(colnames(expr), non_sample_cols)

expr_mat <- as.matrix(expr[, sample_cols])
rownames(expr_mat) <- expr[[gene_col]]

# Remove duplicate gene symbols (keep first occurrence) and NAs
expr_mat <- expr_mat[!is.na(rownames(expr_mat)) & !duplicated(rownames(expr_mat)), ]

cat("\nCleaned expression matrix:", dim(expr_mat), "\n")


# ------------------------------------------------------------
# 3. Check signature gene coverage in the bulk data
# ------------------------------------------------------------

genes_found <- intersect(signature, rownames(expr_mat))
genes_missing <- setdiff(signature, rownames(expr_mat))

cat("\nSignature genes found in bulk data:", length(genes_found), "/", length(signature), "\n")
if (length(genes_missing) > 0) {
  cat("Missing genes (not in bulk platform):\n")
  print(genes_missing)
}

stopifnot(length(genes_found) >= 5)  # need a reasonable number to proceed


# ------------------------------------------------------------
# 4. Compute per-sample signature score (mean z-score method)
# ------------------------------------------------------------

sig_expr <- expr_mat[genes_found, , drop = FALSE]
sig_expr <- apply(sig_expr, 2, as.numeric)
rownames(sig_expr) <- genes_found

# Z-score each gene across samples (row-wise)
sig_expr_z <- t(scale(t(sig_expr)))

# Signature score = mean z-score per sample, across all signature genes
signature_score <- colMeans(sig_expr_z, na.rm = TRUE)

score_df <- data.frame(
  sample_id = names(signature_score),
  signature_score = signature_score
)

cat("\nSignature score distribution:\n")
print(summary(score_df$signature_score))

hist(score_df$signature_score, breaks = 50,
     main = "Bulk Signature Score Distribution", xlab = "Mean z-score")


# ------------------------------------------------------------
# 5. Stratify into High / Low groups (median split — simple,
#    standard, and avoids p-hacking via cutpoint optimization)
# ------------------------------------------------------------

median_score <- median(score_df$signature_score, na.rm = TRUE)
score_df$group <- ifelse(score_df$signature_score > median_score, "High", "Low")

cat("\nGroup sizes:\n")
print(table(score_df$group))


# ------------------------------------------------------------
# 6. Merge with clinical data (survival info + subtype/stage)
# ------------------------------------------------------------

clinical <- fread(clinical_path, sep = "\t", header = TRUE,
                  data.table = FALSE, skip = "PATIENT_ID")

cat("\nClinical data columns:\n")
print(colnames(clinical))

# cBioPortal patient IDs vs sample IDs can differ slightly in
# formatting — inspect both before merging
cat("\nSample score IDs (first 5):\n"); print(head(score_df$sample_id, 5))
cat("Clinical PATIENT_ID (first 5):\n"); print(head(clinical$PATIENT_ID, 5))

# METABRIC sample IDs are usually same as PATIENT_ID; if your
# expression matrix uses slightly different formatting (e.g.
# dashes vs dots), fix it here before merging:
# score_df$sample_id <- gsub("\\.", "-", score_df$sample_id)

merged <- merge(score_df, clinical, by.x = "sample_id", by.y = "PATIENT_ID")

cat("\nMerged dataset:", nrow(merged), "patients matched\n")

write.csv(merged, "results/TME/bulk_signature_scores.csv", row.names = FALSE)

cat("\n============================================\n")
cat("Saved: results/TME/bulk_signature_scores.csv\n")
cat("============================================\n")