# ============================================================
# MANUSCRIPT FIGURE 2 (Cellular landscape of the GSE176078 dataset)
#   Fig. 2a  UMAP of 94,195 cells by cluster identity (31 clusters) -> this script
#   Fig. 2b  UMAP of the four compartments (Immune, Stromal, Tumor, Other)
#            -> script 16 (figures/TME/TME_UMAP.png)
#   Fig. 2c  cell counts per annotated cell type -> script 14
#            (figures/Cell_Composition/)
#   Fig. 2d  canonical marker dot plot -> script 15 (figures/Marker_DotPlot/)
# Run from the project root (no setwd() is used).
# ============================================================

library(Seurat)
library(ggplot2)
library(dplyr)

dir.create("figures/Manuscript_Figures", recursive = TRUE, showWarnings = FALSE)

seu <- readRDS("data/processed/seurat_annotated.rds")
ncol(seu); length(unique(Idents(seu)))      # expected: 94195 cells and 31 clusters

# ------------------------------------------------------------
# Fig. 2a (final): clusters labelled on the plot, legend removed
# ------------------------------------------------------------
p2a <- DimPlot(seu, reduction = "umap", label = TRUE, repel = TRUE) +
  NoLegend() +
  ggtitle("UMAP of Breast Cancer scRNA-seq") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

ggsave("figures/Manuscript_Figures/Fig2a.png", p2a,
       width = 8, height = 7, dpi = 300, bg = "white")
ggsave("figures/UMAP/Fig2a_clusters_nolegend.png", p2a,
       width = 8, height = 7, dpi = 300, bg = "white")

# ------------------------------------------------------------
# Optional: cluster labels that include the dominant SingleR cell type
# ------------------------------------------------------------
lab <- seu@meta.data |>
  dplyr::count(cl = Idents(seu), SingleR) |>
  dplyr::group_by(cl) |>
  dplyr::slice_max(n, n = 1, with_ties = FALSE)

cluster_names <- setNames(paste0(lab$cl, ": ", lab$SingleR), lab$cl)
seu$cl_label <- cluster_names[as.character(Idents(seu))]

p2a_named <- DimPlot(seu, reduction = "umap", group.by = "cl_label", label = FALSE) +
  ggtitle("UMAP of Breast Cancer scRNA-seq") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        legend.text = element_text(size = 7))
ggsave("figures/UMAP/Fig2a_clusters_named.png", p2a_named,
       width = 11, height = 7, dpi = 300, bg = "white")
