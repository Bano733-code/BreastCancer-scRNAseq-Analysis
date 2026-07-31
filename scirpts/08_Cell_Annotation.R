# ==========================================
# Step 8 : Cell Annotation
# ==========================================

library(Seurat)

library(SingleR)

library(celldex)

library(SingleCellExperiment)

library(dplyr)
seurat_obj <- readRDS(
  "data/processed/seurat_umap.rds"
)
sce <- as.SingleCellExperiment(seurat_obj)
ref <- HumanPrimaryCellAtlasData()
pred <- SingleR(
  
  test = sce,
  
  ref = ref,
  
  labels = ref$label.main
  
)
seurat_obj$SingleR <- pred$labels
table(seurat_obj$SingleR)
saveRDS(
  
  seurat_obj,
  
  "data/processed/seurat_annotated.rds"
  
)

cat("Cell annotation completed!\n")