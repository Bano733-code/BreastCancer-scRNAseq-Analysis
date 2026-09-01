# ============================================================
# Malignant Cell State Characterization
# ============================================================

library(Seurat)
library(dplyr)
library(ggplot2)

set.seed(123)

# ------------------------------------------------------------
# 1. Load Step 2 output and subset malignant candidates only
# ------------------------------------------------------------

epi <- readRDS("results/TME/epithelial_cnv_classified.rds")

table(epi$cnv_call)

malignant_cells <- colnames(epi)[epi$cnv_call == "CNV-high (malignant candidate)"]
malignant_obj <- subset(epi, cells = malignant_cells)

cat("Malignant candidate cells:", ncol(malignant_obj), "\n")
print(dim(malignant_obj))


# ------------------------------------------------------------
# 2. Re-normalize and re-cluster the malignant subset alone
#    (higher resolution to capture finer substructure)
# ------------------------------------------------------------

malignant_obj <- NormalizeData(malignant_obj)
malignant_obj <- FindVariableFeatures(malignant_obj, nfeatures = 2000)
malignant_obj <- ScaleData(malignant_obj)
malignant_obj <- RunPCA(malignant_obj, npcs = 30)

# Check how many PCs to use
ElbowPlot(malignant_obj, ndims = 30)

malignant_obj <- FindNeighbors(malignant_obj, dims = 1:20)
malignant_obj <- FindClusters(malignant_obj, resolution = 1.0)
malignant_obj <- RunUMAP(malignant_obj, dims = 1:20)

DimPlot(malignant_obj, label = TRUE) +
  ggtitle("Malignant Epithelial Substructure (Sub-clusters)")

gc()


# ------------------------------------------------------------
# 3. Define breast cancer malignant-state gene signatures
#
#    NOTE: These are established literature markers for these
#    states (widely used across breast cancer scRNA-seq papers).
# ------------------------------------------------------------

basal_like_markers <- c(
  "KRT5", "KRT14", "KRT17", "TP63", "MYLK",
  "ACTA2", "OXTR", "DKK3", "KRT6B"
)

luminal_progenitor_markers <- c(
  "KIT", "ELF5", "KRT23", "ALDH1A3", "SLPI",
  "LTF", "KRT15", "PROM1"
)

mature_luminal_markers <- c(
  "ESR1", "FOXA1", "GATA3", "XBP1", "AR",
  "ANKRD30A", "AGR2", "TFF1"
)

emt_mesenchymal_markers <- c(
  "VIM", "ZEB1", "ZEB2", "SNAI2", "CDH2",
  "FN1", "TWIST1", "MMP2", "S100A4"
)

cycling_proliferative_markers <- c(
  "MKI67", "TOP2A", "CCNB1", "CCNB2", "PCNA",
  "CDK1", "BIRC5", "TYMS", "UBE2C"
)

# Keep only genes present in the dataset
gene_present <- rownames(malignant_obj)

sig_list <- list(
  Basal        = intersect(basal_like_markers, gene_present),
  LumProg      = intersect(luminal_progenitor_markers, gene_present),
  MatureLum    = intersect(mature_luminal_markers, gene_present),
  EMT          = intersect(emt_mesenchymal_markers, gene_present),
  Cycling      = intersect(cycling_proliferative_markers, gene_present)
)

# Sanity check — make sure each signature retained enough genes
print(sapply(sig_list, length))


# ------------------------------------------------------------
# 4. Module scoring (AddModuleScore)
# ------------------------------------------------------------

malignant_obj <- AddModuleScore(
  malignant_obj, features = list(sig_list$Basal),
  name = "Basal_score", ctrl = 50
)
malignant_obj <- AddModuleScore(
  malignant_obj, features = list(sig_list$LumProg),
  name = "LumProg_score", ctrl = 50
)
malignant_obj <- AddModuleScore(
  malignant_obj, features = list(sig_list$MatureLum),
  name = "MatureLum_score", ctrl = 50
)
malignant_obj <- AddModuleScore(
  malignant_obj, features = list(sig_list$EMT),
  name = "EMT_score", ctrl = 50
)
malignant_obj <- AddModuleScore(
  malignant_obj, features = list(sig_list$Cycling),
  name = "Cycling_score", ctrl = 50
)

# AddModuleScore appends "1" to the name by default — rename cleanly
meta_cols <- c("Basal_score1", "LumProg_score1", "MatureLum_score1",
               "EMT_score1", "Cycling_score1")
colnames(malignant_obj@meta.data)[
  colnames(malignant_obj@meta.data) %in% meta_cols
] <- gsub("1$", "", meta_cols)


# ------------------------------------------------------------
# 5. Assign each cell to its dominant state
#    (highest module score wins)
# ------------------------------------------------------------

score_mat <- malignant_obj@meta.data[, c("Basal_score", "LumProg_score",
                                         "MatureLum_score", "EMT_score",
                                         "Cycling_score")]

state_names <- c("Basal-like", "Luminal Progenitor-like",
                 "Mature Luminal-like", "EMT/Mesenchymal",
                 "Cycling/Proliferative")

malignant_obj$malignant_state <- state_names[apply(score_mat, 1, which.max)]

print(table(malignant_obj$malignant_state))


# ------------------------------------------------------------
# 6. Visualize
# ------------------------------------------------------------

DimPlot(malignant_obj, group.by = "malignant_state") +
  ggtitle("Malignant Cell States")

FeaturePlot(malignant_obj, features = c(
  "Basal_score", "LumProg_score", "MatureLum_score",
  "EMT_score", "Cycling_score"
), ncol = 3)

# Cross-tab: does malignant_state align with the Seurat sub-clusters?
print(table(malignant_obj$seurat_clusters, malignant_obj$malignant_state))


# ------------------------------------------------------------
# 7. Save results
# ------------------------------------------------------------

saveRDS(malignant_obj, "results/TME/malignant_states_classified.rds")

state_summary <- malignant_obj@meta.data %>%
  select(orig.ident, seurat_clusters, malignant_state,
         Basal_score, LumProg_score, MatureLum_score,
         EMT_score, Cycling_score)

write.csv(state_summary, "results/TME/malignant_state_scores_per_cell.csv",
          row.names = FALSE)

cat("\n============================================\n")
cat("Saved: results/TME/malignant_states_classified.rds\n")
cat("Saved: results/TME/malignant_state_scores_per_cell.csv\n")
cat("============================================\n")
