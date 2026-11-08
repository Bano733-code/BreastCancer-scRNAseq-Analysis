# ==========================================
# Step 9: Differential Expression Analysis
# ==========================================

library(Seurat)
library(dplyr)

# Load annotated object
seurat_obj <- readRDS(
  "data/processed/seurat_annotated.rds"
)

# Use SingleR labels as identities
Idents(seurat_obj) <- seurat_obj$SingleR

# Differential expression
DEGs <- FindAllMarkers(
  seurat_obj,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

# Create output folder
dir.create(
  "results/Differential_Expression",
  recursive = TRUE,
  showWarnings = FALSE
)

# Save results
write.csv(
  DEGs,
  "results/Differential_Expression/All_DEGs.csv",
  row.names = FALSE
)

cat("Differential Expression completed!\n")