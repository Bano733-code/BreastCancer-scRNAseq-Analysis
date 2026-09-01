# ==========================================
# Step 4: Highly Variable Gene Selection
# ==========================================

library(Seurat)
library(ggplot2)

# Load normalized object
seurat_obj <- readRDS(
  "data/processed/seurat_normalized.rds"
)


# Find highly variable genes

seurat_obj <- FindVariableFeatures(
  seurat_obj,
  selection.method = "vst",
  nfeatures = 2000
)


# Show top variable genes

top10 <- head(
  VariableFeatures(seurat_obj),
  10
)

print(top10)


# Plot variable genes

hvg_plot <- VariableFeaturePlot(seurat_obj)

ggsave(
  filename = "figures/QC/HVG_selection.png",
  plot = hvg_plot,
  width = 8,
  height = 6,
  dpi = 300,
  bg="white"
)


# Save object

saveRDS(
  seurat_obj,
  "data/processed/seurat_hvg.rds"
)


cat("HVG selection completed successfully!\n")
