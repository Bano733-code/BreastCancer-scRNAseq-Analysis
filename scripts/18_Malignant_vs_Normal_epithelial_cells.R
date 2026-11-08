# ============================================================
# Malignant vs Normal Epithelial Cell Classification
# Manual chromosome-level CNV-score method (chromosome arms are NOT resolved)
# (lightweight alternative to inferCNV / CopyKAT — used here
#  due to local machine RAM constraints)
# Input : seurat_TME.rds  (94,195-cell annotated Seurat object)
# Output: epithelial_cnv_classified.rds, cnv_scores_per_cell.csv
#
# [REVISION] Sections 10-13 (end of this script) add the analyses requested
# during peer review: composition of the reference set, a cleaned-reference
# sensitivity analysis, a held-out null, the benchmark against the original
# atlas labels, and the patient-level CNV table (manuscript Table 2).
# They do not change the primary classification produced by sections 1-9.
# ============================================================

# ------------------------------------------------------------
# 0. Setup
# ------------------------------------------------------------

library(Seurat)
library(Matrix)
library(data.table)
library(dplyr)

set.seed(123)


# ------------------------------------------------------------
# 1. Load the full annotated Seurat object
# ------------------------------------------------------------

seurat_obj <- readRDS("results/TME/seurat_TME.rds")

cat("Full object dimensions:\n")
print(dim(seurat_obj))
print(table(seurat_obj$SingleR))


# ------------------------------------------------------------
# 2. Define target (epithelial) and reference (immune+stromal)
#    populations using SingleR labels
# ------------------------------------------------------------

epithelial_cells <- WhichCells(
  seurat_obj,
  expression = SingleR == "Epithelial_cells"
)

reference_labels <- c(
  "T_cells", "B_cell", "NK_cell", "Macrophage",
  "Fibroblasts", "Endothelial_cells",
  "Smooth_muscle_cells", "MSC"
)

# Stratified sampling — max 625 cells per reference population,
# no duplication if a population is smaller than 625
cells_per_group <- 625

reference_cells_list <- lapply(reference_labels, function(label) {
  available_cells <- WhichCells(seurat_obj, expression = SingleR == label)
  n_select <- min(cells_per_group, length(available_cells))
  sample(available_cells, size = n_select, replace = FALSE)
})
names(reference_cells_list) <- reference_labels
reference_cells <- unlist(reference_cells_list, use.names = FALSE)

# [REVISION] store the reference cells for reproducibility (no random number
# is drawn here, so the sampling above is unaffected)
dir.create("results/TME", recursive = TRUE, showWarnings = FALSE)
writeLines(reference_cells, "results/TME/reference_cells_script18.txt")

cat("Epithelial (target) cells:", length(epithelial_cells), "\n")
cat("Reference cells:", length(reference_cells), "\n")


# ------------------------------------------------------------
# 3. Build the CNV working subset (epithelial + reference only)
#    This keeps everything downstream small and RAM-safe
# ------------------------------------------------------------

cnv_cells <- unique(c(epithelial_cells, reference_cells))
cnv_obj <- subset(seurat_obj, cells = cnv_cells)

cnv_obj$InferCNV_group <- "reference"
cnv_obj$InferCNV_group[colnames(cnv_obj) %in% epithelial_cells] <- "epithelial"

cat("CNV working object:\n")
print(dim(cnv_obj))
print(table(cnv_obj$InferCNV_group))

# Free the big object — no longer needed
rm(seurat_obj)
gc()

saveRDS(cnv_obj, "results/TME/infercnv_cnv_subset.rds")


# ------------------------------------------------------------
# 4. Build / load the gene order file (chromosome positions)
#    NOTE: if gene_order.txt already exists from a previous run,
#    load it instead of re-downloading the ~90MB GTF file.
# ------------------------------------------------------------

gene_order_path <- "gene_order.txt"

if (file.exists(gene_order_path)) {
  
  cat("Loading cached gene order file...\n")
  gene_order_final <- fread(
    gene_order_path,
    header = FALSE,
    col.names = c("gene", "chromosome", "start", "end")
  )
  
} else {
  
  cat("No cached gene order found — downloading GENCODE GTF...\n")
  
  options(timeout = 600)
  gencode_url <- paste0(
    "https://ftp.ebi.ac.uk/pub/databases/gencode/",
    "Gencode_human/release_49/gencode.v49.annotation.gtf.gz"
  )
  gtf_file <- "gencode.v49.annotation.gtf.gz"
  download.file(gencode_url, destfile = gtf_file, mode = "wb", method = "libcurl")
  
  gtf <- fread(
    gtf_file, sep = "\t", header = FALSE,
    comment.char = "#", quote = "", data.table = TRUE
  )
  setnames(gtf, c("chromosome", "source", "feature", "start", "end",
                  "score", "strand", "frame", "attribute"))
  
  gtf_genes <- gtf[feature == "gene"]
  gtf_genes[, gene_name := sub('.*gene_name "([^"]+)".*', '\\1', attribute)]
  
  genes <- rownames(LayerData(cnv_obj, assay = "RNA", layer = "counts"))
  
  gene_order <- gtf_genes[
    gene_name %in% genes,
    .(gene = gene_name, chromosome = chromosome, start = start, end = end)
  ]
  gene_order <- gene_order[!duplicated(gene)]
  
  standard_chr <- c(paste0("chr", 1:22), "chrX", "chrY")
  gene_order <- gene_order[chromosome %in% standard_chr]
  
  gene_order[, chromosome := factor(chromosome, levels = standard_chr)]
  setorder(gene_order, chromosome, start)
  gene_order[, chromosome := as.character(chromosome)]
  
  gene_order_final <- gene_order[gene %in% genes]
  
  write.table(gene_order_final, gene_order_path, sep = "\t",
              quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  rm(gtf, gtf_genes)
  gc()
}

cat("Gene order entries:", nrow(gene_order_final), "\n")


# ------------------------------------------------------------
# 5. Extract normalized expression for epithelial + reference
#    (Seurat v5 uses `layer =`, NOT the old `slot =` argument)
# ------------------------------------------------------------

ref_cells <- colnames(cnv_obj)[cnv_obj$InferCNV_group == "reference"]
epi_cells <- colnames(cnv_obj)[cnv_obj$InferCNV_group == "epithelial"]

ref_expr <- GetAssayData(cnv_obj[, ref_cells], assay = "RNA", layer = "data")
epi_expr <- GetAssayData(cnv_obj[, epi_cells], assay = "RNA", layer = "data")

common_genes <- intersect(rownames(epi_expr), gene_order_final$gene)
epi_expr <- epi_expr[common_genes, ]
ref_expr <- ref_expr[common_genes, ]

ref_mean <- Matrix::rowMeans(ref_expr)

cat("Common genes used for CNV scoring:", length(common_genes), "\n")


# ------------------------------------------------------------
# 6. Chromosome-level CNV scoring (mean reference-adjusted expression per chromosome)
#    IMPORTANT: processed ONE chromosome at a time to avoid
#    ever densifying the full sparse matrix (this is what
#    caused the "cannot allocate vector of size 3.2 Gb" crash)
# ------------------------------------------------------------

go <- gene_order_final %>% filter(gene %in% common_genes)
chr_list <- unique(go$chromosome)

chr_score_mat     <- matrix(NA, nrow = ncol(epi_expr), ncol = length(chr_list),
                            dimnames = list(colnames(epi_expr), chr_list))
ref_chr_score_mat <- matrix(NA, nrow = ncol(ref_expr), ncol = length(chr_list),
                            dimnames = list(colnames(ref_expr), chr_list))

for (ch in chr_list) {
  
  genes_in_chr <- go$gene[go$chromosome == ch]
  genes_in_chr <- intersect(genes_in_chr, rownames(epi_expr))
  
  if (length(genes_in_chr) >= 5) {
    
    epi_sub  <- as.matrix(epi_expr[genes_in_chr, , drop = FALSE])
    ref_sub  <- as.matrix(ref_expr[genes_in_chr, , drop = FALSE])
    mean_sub <- ref_mean[genes_in_chr]
    
    epi_rel <- epi_sub - mean_sub
    ref_rel <- ref_sub - mean_sub
    
    chr_score_mat[, ch]     <- colMeans(epi_rel)
    ref_chr_score_mat[, ch] <- colMeans(ref_rel)
    
    rm(epi_sub, ref_sub, epi_rel, ref_rel)
  }
  
  cat("Done:", ch, "\n")
}

gc()


# ------------------------------------------------------------
# 7. Per-cell CNV score (SD across chromosomes) + cutoff
#    [REVISION] A chromosome is scored only if >= 5 of its genes are present
#    (chrY does not qualify, so 23 chromosomes contribute: chr1-22 and chrX).
#    The reference cells define both the per-gene baseline and the cutoff;
#    see sections 10-12 for the sensitivity analyses of this choice.
#    Cutoff = 95th percentile of the REFERENCE score
#    distribution (data-driven, NOT an arbitrary 50/50 split)
# ------------------------------------------------------------

cnv_score     <- apply(chr_score_mat, 1, function(x) sd(x, na.rm = TRUE))
ref_cnv_score <- apply(ref_chr_score_mat, 1, function(x) sd(x, na.rm = TRUE))

cnv_cutoff <- quantile(ref_cnv_score, 0.95, na.rm = TRUE)
cat("CNV cutoff (95th pct of reference):", cnv_cutoff, "\n")

epi_meta <- data.frame(
  barcode   = names(cnv_score),
  cnv_score = cnv_score
)
epi_meta$cnv_call <- ifelse(
  epi_meta$cnv_score > cnv_cutoff,
  "CNV-high (malignant candidate)",
  "CNV-low"
)

print(table(epi_meta$cnv_call))


# ------------------------------------------------------------
# 8. Sanity checks — ALWAYS inspect before trusting the split
# ------------------------------------------------------------

hist(ref_cnv_score, breaks = 50, main = "CNV Score Distribution (Reference)",
     xlab = "cnv_score", col = "grey80")

hist(epi_meta$cnv_score, breaks = 100, main = "CNV Score Distribution (Epithelial)",
     xlab = "cnv_score", col = "grey80")
abline(v = cnv_cutoff, col = "red", lwd = 2)

# Per-sample breakdown — malignant proportion should vary
# meaningfully across samples, not be uniform everywhere
epi_meta$sample <- sub("_[ACGT]+$", "", epi_meta$barcode)
print(table(epi_meta$sample, epi_meta$cnv_call))


# ------------------------------------------------------------
# 9. Attach labels back to the epithelial Seurat object and save
# ------------------------------------------------------------

epi <- subset(cnv_obj, cells = epi_cells)

epi$cnv_score <- epi_meta$cnv_score[match(colnames(epi), epi_meta$barcode)]
epi$cnv_call  <- epi_meta$cnv_call[match(colnames(epi), epi_meta$barcode)]

print(table(epi$cnv_call))

saveRDS(epi, "results/TME/epithelial_cnv_classified.rds")
write.csv(epi_meta, "results/TME/cnv_scores_per_cell.csv", row.names = FALSE)

cat("\n============================================\n")
cat("Saved: results/TME/epithelial_cnv_classified.rds\n")
cat("Saved: results/TME/cnv_scores_per_cell.csv\n")
cat("============================================\n")


# ============================================================
# NOTE FOR METHODS SECTION:
# This is a simplified chromosome-level expression-averaging
# proxy for CNV signal (conceptually related to inferCNV),
# used due to computational resource constraints. Cells are
# labeled "CNV-high (malignant candidate)" / "CNV-low" rather
# than definitively "malignant" / "normal", since this method
# has not been validated against a gold-standard tool
# (inferCNV / CopyKAT) on this dataset.
# ============================================================


# ============================================================
# [REVISION] 10. Reference composition according to the ORIGINAL atlas labels
#     (peer-review addition)
#     The atlas labels show that part of the reference (mainly the cells that
#     SingleR called "MSC" and some "Smooth_muscle_cells") are cancer-epithelial
#     cells, which lowers the contrast of the baseline and raises the cutoff.
# ============================================================

library(tidyr)
wu <- read.csv("data/GSE176078/GSE176078/metadata.csv",
               row.names = 1, check.names = FALSE)

sr_ref    <- cnv_obj$SingleR[ref_cells]
ref_atlas <- wu[ref_cells, "celltype_major"]

ref_tab <- data.frame(SingleR = sr_ref, atlas = ref_atlas) %>%
  group_by(SingleR) %>%
  summarise(n = n(),
            atlas_cancer = sum(atlas == "Cancer Epithelial", na.rm = TRUE),
            atlas_cancer_pct = round(100 * mean(atlas == "Cancer Epithelial", na.rm = TRUE), 1),
            .groups = "drop") %>%
  arrange(desc(atlas_cancer))
print(as.data.frame(ref_tab), row.names = FALSE)
cat("Atlas cancer-epithelial cells inside the reference:",
    sum(ref_tab$atlas_cancer), "of", length(ref_cells), "\n")   # 474 of 4743 (10.0%)
write.csv(ref_tab, "results/TME/CNV_reference_composition_atlas.csv", row.names = FALSE)


# ============================================================
# [REVISION] 11. Sensitivity analysis: cleaned reference
#     MSC- and smooth-muscle-labelled cells are removed from the reference
#     using the SingleR labels only (the atlas labels are used only afterwards
#     for the benchmark, so the benchmark is not circular).
# ============================================================

score_cells <- function(expr_mat, baseline) {
  # SD across chromosomes of the mean reference-adjusted expression per chromosome
  m <- sapply(chr_list, function(ch) {
    g <- intersect(go$gene[go$chromosome == ch], rownames(expr_mat))
    if (length(g) >= 5) {
      colMeans(as.matrix(expr_mat[g, , drop = FALSE]) - baseline[g])
    } else {
      rep(NA_real_, ncol(expr_mat))
    }
  })
  out <- apply(m, 1, sd, na.rm = TRUE)
  names(out) <- colnames(expr_mat)
  out
}

ref_clean_cells <- ref_cells[!sr_ref %in% c("MSC", "Smooth_muscle_cells")]
ref_mean_clean  <- Matrix::rowMeans(ref_expr[, ref_clean_cells])

score_epi_clean <- score_cells(epi_expr, ref_mean_clean)
score_ref_clean <- score_cells(ref_expr[, ref_clean_cells], ref_mean_clean)
cutoff_clean    <- quantile(score_ref_clean, 0.95, na.rm = TRUE)
cat("Cleaned reference:", length(ref_clean_cells), "cells | new cutoff:",
    round(cutoff_clean, 4), "| original cutoff:", round(cnv_cutoff, 4), "\n")


# ============================================================
# [REVISION] 12. Held-out null: immune/stromal cells that are NOT part of the
#     reference are scored against the cleaned baseline, so the null is not
#     derived from the cells that define the baseline.
# ============================================================

obj_all <- readRDS("results/TME/seurat_TME.rds")
set.seed(123)
heldout_labels <- c("T_cells", "B_cell", "NK_cell", "Macrophage",
                    "Fibroblasts", "Endothelial_cells")
heldout_cells <- unlist(lapply(heldout_labels, function(l) {
  cand <- setdiff(WhichCells(obj_all, expression = SingleR == l), reference_cells)
  sample(cand, min(300, length(cand)))
}))
heldout_expr <- GetAssayData(obj_all, assay = "RNA", layer = "data")[common_genes, heldout_cells]
rm(obj_all); gc()

score_heldout  <- score_cells(heldout_expr, ref_mean_clean)
cutoff_heldout <- quantile(score_heldout, 0.95, na.rm = TRUE)
cat("Held-out cells:", length(heldout_cells), "| held-out cutoff:",
    round(cutoff_heldout, 4), "\n")                       # ~0.0209 vs 0.0211


# ============================================================
# [REVISION] 13. Benchmark against the cancer/normal epithelial labels of the
#     original atlas (a reference method, not ground truth) and Table 2
# ============================================================

epi_md <- epi@meta.data
epi_md$barcode   <- rownames(epi_md)
epi_md$new_score <- score_epi_clean[epi_md$barcode]
epi_md$old_high  <- epi_md$cnv_call == "CNV-high (malignant candidate)"
epi_md$new_high  <- epi_md$new_score > cutoff_clean
epi_md$wu        <- wu[epi_md$barcode, "celltype_major"]
epi_md$subtype   <- wu[epi_md$barcode, "subtype"]
cat("CNV-high: original", sum(epi_md$old_high), "| cleaned reference", sum(epi_md$new_high), "\n")
cat("Spearman correlation of original and cleaned scores:",
    round(cor(epi_md$cnv_score, epi_md$new_score, method = "spearman"), 3), "\n")

bench_cells <- epi_md %>%
  filter(wu %in% c("Cancer Epithelial", "Normal Epithelial")) %>%
  mutate(wu_cancer = wu == "Cancer Epithelial")

if (!requireNamespace("pROC", quietly = TRUE)) install.packages("pROC")
bench <- function(call, score) {
  TP <- sum(bench_cells$wu_cancer & call);  FN <- sum(bench_cells$wu_cancer & !call)
  TN <- sum(!bench_cells$wu_cancer & !call); FP <- sum(!bench_cells$wu_cancer & call)
  N  <- nrow(bench_cells)
  po <- (TP + TN) / N
  pe <- ((TP + FN) * (TP + FP) + (TN + FP) * (TN + FN)) / N^2
  data.frame(n_high = sum(call), sensitivity = TP / (TP + FN), specificity = TN / (TN + FP),
             precision = TP / (TP + FP), accuracy = po, kappa = (po - pe) / (1 - pe),
             AUC = as.numeric(pROC::auc(pROC::roc(bench_cells$wu_cancer, score, quiet = TRUE))))
}
bench_tab <- rbind(original_reference = bench(bench_cells$old_high, bench_cells$cnv_score),
                   cleaned_reference  = bench(bench_cells$new_high, bench_cells$new_score))
print(round(bench_tab, 3))
write.csv(cbind(analysis = rownames(bench_tab), bench_tab),
          "results/TME/CNV_benchmark_summary.csv", row.names = FALSE)

# held-out threshold applied to the cleaned-reference scores
heldout_tab <- do.call(rbind, lapply(
  list(`in-sample cleaned reference` = cutoff_clean, `held-out cells` = cutoff_heldout),
  function(cut) {
    call <- bench_cells$new_score > cut
    TP <- sum(bench_cells$wu_cancer & call);  FN <- sum(bench_cells$wu_cancer & !call)
    TN <- sum(!bench_cells$wu_cancer & !call); FP <- sum(!bench_cells$wu_cancer & call)
    data.frame(threshold = round(cut, 4), CNV_high_cells = sum(epi_md$new_score > cut),
               sensitivity = round(TP / (TP + FN), 3), specificity = round(TN / (TN + FP), 3))
  }))
heldout_tab <- cbind(null = rownames(heldout_tab), heldout_tab)
print(heldout_tab)
write.csv(heldout_tab, "results/TME/CNV_heldout_null.csv", row.names = FALSE)

# per-patient agreement with the atlas labels (original call)
bench_patient <- bench_cells %>%
  group_by(Patient = orig.ident) %>%
  summarise(n = n(), atlas_cancer_pct = round(100 * mean(wu_cancer), 1),
            ours_high_pct = round(100 * mean(old_high), 1),
            agree_pct = round(100 * mean(wu_cancer == old_high), 1), .groups = "drop")
write.csv(bench_patient, "results/TME/CNV_benchmark_per_patient.csv", row.names = FALSE)

# Table 2: patient-level CNV table for all 26 patients
all_patients <- wu %>% distinct(Patient = orig.ident, Subtype = subtype)
table2 <- epi_md %>%
  group_by(Patient = orig.ident) %>%
  summarise(Epithelial = n(), CNV_high = sum(old_high), .groups = "drop") %>%
  right_join(all_patients, by = "Patient") %>%
  mutate(across(c(Epithelial, CNV_high), ~ replace_na(.x, 0)),
         CNV_low = Epithelial - CNV_high,
         Pct_CNV_high = ifelse(Epithelial > 0, round(100 * CNV_high / Epithelial, 2), NA)) %>%
  left_join(epi_md %>% group_by(Patient = orig.ident) %>%
              summarise(Pct_CNV_high_cleanRef = round(100 * mean(new_high), 1), .groups = "drop"),
            by = "Patient") %>%
  arrange(Subtype, Patient)
print(as.data.frame(table2), row.names = FALSE)
cat("Patients:", nrow(table2), "| epithelial:", sum(table2$Epithelial),
    "| CNV-high:", sum(table2$CNV_high), "\n")           # 26 | 22729 | 10548
write.csv(table2, "results/TME/Table2_per_patient_CNV.csv", row.names = FALSE)

# per-cell scores of the sensitivity analysis (supplementary table S4)
write.csv(epi_md[, c("barcode", "orig.ident", "subtype", "cnv_score", "new_score",
                     "old_high", "new_high", "wu")],
          "results/TME/CNV_sensitivity_cleanRef_cells.csv", row.names = FALSE)
