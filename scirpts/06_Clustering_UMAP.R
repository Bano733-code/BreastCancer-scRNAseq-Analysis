# ==========================================
# Step 6.1: Find Neighbors
# ==========================================

library(Seurat)

# Load PCA object
seurat_obj <- readRDS(
  "data/processed/seurat_pca.rds"
)

# Build nearest-neighbor graph
seurat_obj <- FindNeighbors(
  seurat_obj,
  dims = 1:20
)

cat("Neighbor graph created successfully!\n")

# Save object
saveRDS(
  seurat_obj,
  "data/processed/seurat_neighbors.rds"
)
# ==========================================
# Step 6.2: Cell Clustering
# ==========================================

library(Seurat)

# Load neighbor object (if restarting R)
seurat_obj <- readRDS(
  "data/processed/seurat_neighbors.rds"
)

# Find clusters
seurat_obj <- FindClusters(
  seurat_obj,
  resolution = 0.5
)

# Show cluster sizes
table(Idents(seurat_obj))

# Save clustered object
saveRDS(
  seurat_obj,
  "data/processed/seurat_clustered.rds"
)

cat("Clustering completed successfully!\n")
# ==========================================
# Step 6.3: UMAP Visualization
# ==========================================

library(Seurat)
library(ggplot2)

# Load clustered object (if restarting R)
seurat_obj <- readRDS(
  "data/processed/seurat_clustered.rds"
)

# Run UMAP
seurat_obj <- RunUMAP(
  seurat_obj,
  dims = 1:20
)

# Create UMAP plot
umap_plot <- DimPlot(
  seurat_obj,
  reduction = "umap",
  label = TRUE,
  repel = TRUE
) +
  ggtitle("UMAP of Breast Cancer scRNA-seq") +
  theme_classic()

# Display
umap_plot

# Create folder
dir.create(
  "figures/UMAP",
  recursive = TRUE,
  showWarnings = FALSE
)

# Save
ggsave(
  "figures/UMAP/UMAP.png",
  plot = umap_plot,
  width = 9,
  height = 7,
  dpi = 300
)

# Save object
saveRDS(
  seurat_obj,
  "data/processed/seurat_umap.rds"
)

cat("UMAP completed successfully!\n")