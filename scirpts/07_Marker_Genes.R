# ==========================================
# Step 7: Marker Gene Detection
# ==========================================

library(Seurat)
library(dplyr)

# Load UMAP object
seurat_obj <- readRDS(
  "data/processed/seurat_umap.rds"
)

# Find marker genes
markers <- FindAllMarkers(
  seurat_obj,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

# Show first few rows
head(markers)

# Save markers
write.csv(
  markers,
  "results/Marker_Genes/All_Markers.csv",
  row.names = FALSE
)

# Save object
saveRDS(
  seurat_obj,
  "data/processed/seurat_markers.rds"
)

cat("Marker gene detection completed!\n")