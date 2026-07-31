# ============================
# Step 1: Load Dataset
# ============================

library(Seurat)
library(Matrix)
library(dplyr)

# Dataset path
data_path <- "C:/Users/banor/Downloads/BreastCancer_scRNAseq/Single_cell_RNAseq/data/GSE176078/GSE176078/"

# Read matrix
counts <- ReadMtx(
  mtx = paste0(data_path, "matrix.mtx"),
  cells = paste0(data_path, "barcodes.tsv"),
  features = paste0(data_path, "genes.tsv"),
  feature.column = 1
)

# Check dimensions
dim(counts)

# Create Seurat object
seurat_obj <- CreateSeuratObject(
  counts = counts,
  project = "GSE176078_BreastCancer",
  min.cells = 3,
  min.features = 200
)

# Show object summary
seurat_obj