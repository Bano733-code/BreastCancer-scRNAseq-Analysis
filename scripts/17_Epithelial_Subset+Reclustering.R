###############################################################
# EPITHELIAL CELL SUBSETTING + RE-CLUSTERING
# Dataset: GSE176078 Breast Cancer scRNA-seq
###############################################################

library(Seurat)
library(dplyr)
library(ggplot2)

###############################################################
# 1. Create output directories
###############################################################

dir.create(
  "results/Malignant",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "figures/Malignant",
  recursive = TRUE,
  showWarnings = FALSE
)

###############################################################
# 2. Check original Seurat object
###############################################################

cat("\n========================================\n")
cat("ORIGINAL SEURAT OBJECT\n")
cat("========================================\n")

cat(
  "Total cells:",
  ncol(seurat_obj),
  "\n"
)

cat(
  "Total genes:",
  nrow(seurat_obj),
  "\n"
)

###############################################################
# 3. Check SingleR annotations
###############################################################

cat("\n========================================\n")
cat("SingleR CELL-TYPE DISTRIBUTION\n")
cat("========================================\n")

print(
  sort(
    table(seurat_obj$SingleR),
    decreasing = TRUE
  )
)

###############################################################
# 4. Check epithelial label
###############################################################

cat("\nChecking epithelial annotation...\n")

if (!"Epithelial_cells" %in% unique(seurat_obj$SingleR)) {
  
  stop(
    paste0(
      "\n'Epithelial_cells' was not found in SingleR annotations.\n",
      "Please check the output of table(seurat_obj$SingleR)\n",
      "and replace 'Epithelial_cells' with the correct label."
    )
  )
  
}

cat("Epithelial label found: Epithelial_cells\n")

###############################################################
# 5. Subset epithelial cells
###############################################################

cat("\nSubsetting epithelial cells...\n")

epithelial_obj <- subset(
  seurat_obj,
  subset = SingleR == "Epithelial_cells"
)

cat(
  "Epithelial cells:",
  ncol(epithelial_obj),
  "\n"
)

cat(
  "Percentage of total cells:",
  round(
    ncol(epithelial_obj) /
      ncol(seurat_obj) * 100,
    2
  ),
  "%\n"
)

###############################################################
# 6. Save epithelial subset before re-clustering
###############################################################

saveRDS(
  epithelial_obj,
  "results/Malignant/epithelial_raw.rds"
)

###############################################################
# 7. Epithelial cells on original UMAP
###############################################################

if ("umap" %in% names(seurat_obj@reductions)) {
  
  epithelial_umap <- DimPlot(
    seurat_obj,
    cells.highlight = Cells(epithelial_obj),
    cols.highlight = "red",
    cols = "grey90",
    pt.size = 0.3
  ) +
    ggtitle(
      "Epithelial Cells Highlighted on Original UMAP"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold"
      )
    )
  
  print(epithelial_umap)
  
  ggsave(
    "figures/Malignant/Epithelial_Cells_Original_UMAP.png",
    epithelial_umap,
    width = 9,
    height = 7,
    dpi = 300
  )
  
}

###############################################################
# 8. Epithelial cells by sample
###############################################################

if ("orig.ident" %in% colnames(epithelial_obj@meta.data)) {
  
  sample_counts <- as.data.frame(
    table(
      epithelial_obj$orig.ident
    )
  )
  
  colnames(sample_counts) <- c(
    "Sample",
    "Cells"
  )
  
  sample_counts$Percentage <- round(
    sample_counts$Cells /
      sum(sample_counts$Cells) * 100,
    2
  )
  
  print(sample_counts)
  
  write.csv(
    sample_counts,
    "results/Malignant/Epithelial_Cell_Counts_By_Sample.csv",
    row.names = FALSE
  )
  
  sample_plot <- ggplot(
    sample_counts,
    aes(
      x = reorder(Sample, Cells),
      y = Cells
    )
  ) +
    geom_col() +
    coord_flip() +
    theme_bw(base_size = 13) +
    labs(
      title = "Epithelial Cells by Sample",
      x = "Sample",
      y = "Number of Epithelial Cells"
    )
  
  ggsave(
    "figures/Malignant/Epithelial_Cells_By_Sample.png",
    sample_plot,
    width = 9,
    height = 7,
    dpi = 300
  )
  
}

###############################################################
# 9. Save epithelial metadata
###############################################################

write.csv(
  epithelial_obj@meta.data,
  "results/Malignant/Epithelial_Cell_Metadata.csv",
  row.names = TRUE
)

###############################################################
# 10. Set RNA assay
###############################################################

DefaultAssay(epithelial_obj) <- "RNA"

###############################################################
# 11. Re-normalize epithelial cells
###############################################################

cat("\n========================================\n")
cat("EPITHELIAL RE-CLUSTERING\n")
cat("========================================\n")

cat("\nRunning normalization...\n")

epithelial_obj <- NormalizeData(
  epithelial_obj,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

###############################################################
# 12. Find highly variable genes
###############################################################

cat("\nFinding highly variable genes...\n")

epithelial_obj <- FindVariableFeatures(
  epithelial_obj,
  selection.method = "vst",
  nfeatures = 3000,
  verbose = FALSE
)

cat(
  "Variable genes:",
  length(VariableFeatures(epithelial_obj)),
  "\n"
)

###############################################################
# 13. Variable feature plot
###############################################################

p_variable <- VariableFeaturePlot(
  epithelial_obj
) +
  ggtitle(
    "Epithelial Cells - Highly Variable Genes"
  ) +
  theme_bw()

print(p_variable)

ggsave(
  "figures/Malignant/Epithelial_Variable_Features.png",
  p_variable,
  width = 9,
  height = 7,
  dpi = 300
)

###############################################################
# 14. Scale data
###############################################################

cat("\nScaling data...\n")

# Only regress percent.mt if the metadata column exists

if ("percent.mt" %in% colnames(epithelial_obj@meta.data)) {
  
  epithelial_obj <- ScaleData(
    epithelial_obj,
    vars.to.regress = "percent.mt",
    verbose = FALSE
  )
  
} else {
  
  epithelial_obj <- ScaleData(
    epithelial_obj,
    verbose = FALSE
  )
  
}

###############################################################
# 15. PCA
###############################################################

cat("\nRunning PCA...\n")

epithelial_obj <- RunPCA(
  epithelial_obj,
  features = VariableFeatures(epithelial_obj),
  npcs = 40,
  verbose = FALSE
)

###############################################################
# 16. PCA elbow plot
###############################################################

p_elbow <- ElbowPlot(
  epithelial_obj,
  ndims = 40
) +
  ggtitle(
    "Epithelial Cells - PCA Elbow Plot"
  ) +
  theme_bw()

print(p_elbow)

ggsave(
  "figures/Malignant/Epithelial_PCA_Elbow.png",
  p_elbow,
  width = 8,
  height = 6,
  dpi = 300
)

###############################################################
# 17. Select PCs
###############################################################

# Initial analysis uses 30 PCs.
# This should be reviewed using the elbow plot.

use_pcs <- 1:30

cat(
  "\nPCs used for clustering:",
  paste(range(use_pcs), collapse = "-"),
  "\n"
)

###############################################################
# 18. Find neighbors
###############################################################

cat("\nFinding nearest neighbors...\n")

epithelial_obj <- FindNeighbors(
  epithelial_obj,
  dims = use_pcs,
  verbose = FALSE
)

###############################################################
# 19. Higher-resolution clustering
###############################################################

cat("\nRunning epithelial clustering...\n")

epithelial_obj <- FindClusters(
  epithelial_obj,
  resolution = 1.0,
  algorithm = 1,
  verbose = FALSE
)

###############################################################
# 20. Cluster statistics
###############################################################

cluster_counts <- sort(
  table(
    Idents(epithelial_obj)
  ),
  decreasing = TRUE
)

cat("\n========================================\n")
cat("EPITHELIAL CLUSTER COUNTS\n")
cat("========================================\n")

print(cluster_counts)

write.csv(
  as.data.frame(cluster_counts),
  "results/Malignant/Epithelial_Cluster_Cell_Counts.csv",
  row.names = FALSE
)

###############################################################
# 21. Run UMAP
###############################################################

cat("\nRunning UMAP...\n")

epithelial_obj <- RunUMAP(
  epithelial_obj,
  dims = use_pcs,
  verbose = FALSE
)

###############################################################
# 22. UMAP
###############################################################

p_umap <- DimPlot(
  epithelial_obj,
  reduction = "umap",
  label = TRUE,
  repel = TRUE,
  pt.size = 0.4
) +
  ggtitle(
    "Epithelial Cell Subclusters"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    )
  )

print(p_umap)

ggsave(
  "figures/Malignant/Epithelial_Subclusters_UMAP.png",
  p_umap,
  width = 10,
  height = 8,
  dpi = 300
)

###############################################################
# 23. Save re-clustered object
###############################################################

saveRDS(
  epithelial_obj,
  "results/Malignant/epithelial_reclustered.rds"
)

###############################################################
# 24. Save metadata
###############################################################

write.csv(
  epithelial_obj@meta.data,
  "results/Malignant/Epithelial_Reclustered_Metadata.csv",
  row.names = TRUE
)

###############################################################
# 25. Final summary
###############################################################

cat("\n========================================\n")
cat("PHASE 1 - STEP 1 COMPLETED\n")
cat("========================================\n")

cat(
  "Original cells:",
  ncol(seurat_obj),
  "\n"
)

cat(
  "Epithelial cells:",
  ncol(epithelial_obj),
  "\n"
)

cat(
  "Epithelial percentage:",
  round(
    ncol(epithelial_obj) /
      ncol(seurat_obj) * 100,
    2
  ),
  "%\n"
)

cat(
  "Number of epithelial clusters:",
  length(
    unique(
      Idents(epithelial_obj)
    )
  ),
  "\n"
)

cat(
  "\nSaved object:\n",
  "results/Malignant/epithelial_reclustered.rds\n"
)
