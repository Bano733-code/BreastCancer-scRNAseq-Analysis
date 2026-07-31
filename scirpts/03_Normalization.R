# ==========================================
# Step 3: Normalization
# ==========================================

library(Seurat)

# Load filtered object
seurat_obj <- readRDS("data/processed/seurat_filtered.rds")

# Normalize data
seurat_obj <- NormalizeData(
  seurat_obj,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)

# Check assay
DefaultAssay(seurat_obj)

# Save normalized object
saveRDS(
  seurat_obj,
  file = "data/processed/seurat_normalized.rds"
)

cat("Normalization completed successfully!\n")
