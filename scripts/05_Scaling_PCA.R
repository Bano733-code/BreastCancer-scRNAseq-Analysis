# ==========================================
# Step 5.1: Scaling Data
# ==========================================

library(Seurat)

# Load HVG object
seurat_obj <- readRDS(
  "data/processed/seurat_hvg.rds"
)

# Scale only variable genes
seurat_obj <- ScaleData(
  seurat_obj,
  features = VariableFeatures(seurat_obj)
)

cat("Scaling completed successfully!\n")

# Save object
saveRDS(
  seurat_obj,
  "data/processed/seurat_scaled.rds"
)
# ==========================================
# Step 5.2: Principal Component Analysis
# ==========================================

# Run PCA using only highly variable genes
seurat_obj <- RunPCA(
  seurat_obj,
  features = VariableFeatures(seurat_obj)
)

# Print PCA summary
print(seurat_obj[["pca"]])

# Save object
saveRDS(
  seurat_obj,
  "data/processed/seurat_pca.rds"
)

cat("PCA completed successfully!\n")
# ==========================================
# Step 5.3: PCA Plot
# ==========================================

library(ggplot2)

# Create folder if it doesn't exist
dir.create("figures/PCA", recursive = TRUE, showWarnings = FALSE)

# PCA scatter plot
pca_plot <- DimPlot(
  seurat_obj,
  reduction = "pca"
) +
  ggtitle("Principal Component Analysis (PCA)") +
  theme_classic()

# Display plot
pca_plot

# Save figure
ggsave(
  filename = "figures/PCA/PCA_plot.png",
  plot = pca_plot,
  width = 8,
  height = 6,
  dpi = 300
)

cat("PCA plot saved successfully!\n")
# ==========================================
# Step 5.4: Elbow Plot
# ==========================================

library(Seurat)
library(ggplot2)

# Load PCA object (if restarting R)
seurat_obj <- readRDS(
  "data/processed/seurat_pca.rds"
)


# Generate Elbow Plot

elbow_plot <- ElbowPlot(
  seurat_obj,
  ndims = 50
) +
  ggtitle("Elbow Plot - PCA Dimensions") +
  theme_classic()


# Display plot
elbow_plot


# Create folder
dir.create(
  "figures/PCA",
  recursive = TRUE,
  showWarnings = FALSE
)


# Save figure

ggsave(
  filename = "figures/PCA/Elbow_plot.png",
  plot = elbow_plot,
  width = 8,
  height = 6,
  dpi = 300
)


cat("Elbow plot saved successfully!\n")
