# ============================================================
# Malignant vs Normal Epithelial Cell Classification
# Manual chromosome-arm CNV-score method
# (lightweight alternative to inferCNV / CopyKAT — used here
#  due to local machine RAM constraints)
# Input : seurat_TME.rds  (94,195-cell annotated Seurat object)
# Output: epithelial_cnv_classified.rds, cnv_scores_per_cell.csv
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
# 6. Chromosome-arm CNV scoring
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
# 7. Per-cell CNV score (SD across chromosome arms) + cutoff
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
# This is a simplified chromosome-arm expression-averaging
# proxy for CNV signal (conceptually related to inferCNV),
# used due to computational resource constraints. Cells are
# labeled "CNV-high (malignant candidate)" / "CNV-low" rather
# than definitively "malignant" / "normal", since this method
# has not been validated against a gold-standard tool
# (inferCNV / CopyKAT) on this dataset.
# ============================================================