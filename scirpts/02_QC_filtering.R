# ==========================================
# Step 2: Quality Control
# ==========================================

library(Seurat)
library(ggplot2)
library(patchwork)

# Load Seurat object
seurat_obj <- readRDS("data/processed/seurat_raw.rds")

# Calculate mitochondrial percentage
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(
  seurat_obj,
  pattern = "^MT-"
)

# Basic information
head(seurat_obj@meta.data)

# Summary statistics
summary(seurat_obj$nFeature_RNA)
summary(seurat_obj$nCount_RNA)
summary(seurat_obj$percent.mt)
# ==========================================
# QC Plots
# ==========================================

library(patchwork)

# Violin plots
vln_plot <- VlnPlot(
  seurat_obj,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3,
  pt.size = 0.1
)

# Save violin plot
ggsave(
  filename = "figures/QC/QC_Violin.png",
  plot = vln_plot,
  width = 10,
  height = 5,
  dpi = 300
)

# Scatter plot 1
scatter1 <- FeatureScatter(
  seurat_obj,
  feature1 = "nCount_RNA",
  feature2 = "percent.mt"
)

# Scatter plot 2
scatter2 <- FeatureScatter(
  seurat_obj,
  feature1 = "nCount_RNA",
  feature2 = "nFeature_RNA"
)

# Combine plots
qc_scatter <- scatter1 + scatter2

# Save scatter plots
ggsave(
  filename = "figures/QC/QC_Scatter.png",
  plot = qc_scatter,
  width = 10,
  height = 5,
  dpi = 300
)
# ==========================================
# Step 2.3: Filter low-quality cells
# ==========================================

seurat_filtered <- subset(
  seurat_obj,
  subset =
    nFeature_RNA > 300 &
    nFeature_RNA < 7000 &
    percent.mt < 15
)

# Compare before vs after
cat("Cells before filtering :", ncol(seurat_obj), "\n")
cat("Cells after filtering  :", ncol(seurat_filtered), "\n")

# Save filtered object
saveRDS(
  seurat_filtered,
  file = "data/processed/seurat_filtered.rds"
)