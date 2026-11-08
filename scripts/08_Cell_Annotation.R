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

# ==========================================
# [REVISION] Concordance of the SingleR labels with the original atlas
# annotation (added during peer review; no label used downstream changes).
# Compares the Human Primary Cell Atlas labels (SingleR) with the
# `celltype_major` labels distributed with GSE176078 (Wu et al., 2021).
# Output: results/Annotation/SingleR_vs_atlas_crosstab.csv
# ==========================================

library(dplyr)
dir.create("results/Annotation", recursive = TRUE, showWarnings = FALSE)

atlas <- read.csv("data/GSE176078/GSE176078/metadata.csv",
                  row.names = 1, check.names = FALSE)
md <- seurat_obj@meta.data
md$atlas_major <- atlas[rownames(md), "celltype_major"]
cat("Cells:", nrow(md), "| matched to atlas labels:", sum(!is.na(md$atlas_major)), "\n")

# 1. Full cross-tabulation (SingleR label x atlas label)
tab <- table(SingleR = md$SingleR, Atlas = md$atlas_major)
write.csv(as.data.frame.matrix(tab), "results/Annotation/SingleR_vs_atlas_crosstab.csv")
pc  <- round(100 * prop.table(tab, 1), 1)
big <- rowSums(tab) >= 100
print(cbind(n = rowSums(tab)[big], pc[big, ]))

# 2. Concordance for labels that have a direct counterpart in the atlas
label_map <- c(T_cells = "T-cells", NK_cell = "T-cells", B_cell = "B-cells",
               Macrophage = "Myeloid", Monocyte = "Myeloid", DC = "Myeloid",
               Neutrophils = "Myeloid", Endothelial_cells = "Endothelial",
               Fibroblasts = "CAFs", Smooth_muscle_cells = "PVL",
               Epithelial_cells = "EPI")
md$expected  <- unname(label_map[md$SingleR])
md$atlas_grp <- ifelse(md$atlas_major %in% c("Cancer Epithelial", "Normal Epithelial"),
                       "EPI", md$atlas_major)
mapped <- md[!is.na(md$expected) & !is.na(md$atlas_grp), ]
conc <- mapped %>% group_by(SingleR) %>%
  summarise(n = n(), concordant_pct = round(100 * mean(expected == atlas_grp), 1))
print(as.data.frame(conc), row.names = FALSE)
cat("\nOverall concordance (mapped labels):",
    round(100 * mean(mapped$expected == mapped$atlas_grp), 1), "% of", nrow(mapped), "cells\n")

# 3. Labels without a direct counterpart (e.g. tissue stem cells, chondrocytes,
#    iPS cells, MSC): what are they according to the atlas?
unmapped <- md[is.na(md$expected), ]
cat("\nUnmapped SingleR labels:", nrow(unmapped), "cells (",
    round(100 * nrow(unmapped) / nrow(md), 1), "%)\n")
print(head(sort(table(unmapped$SingleR), decreasing = TRUE), 15))
cat("\nAtlas major label of these cells (%):\n")
print(round(100 * prop.table(table(unmapped$atlas_major, useNA = "ifany")), 1))
