###############################################################
## Marker Heatmap
###############################################################

library(Seurat)
library(dplyr)
library(ggplot2)

dir.create(
  "figures/Marker_Heatmap",
  recursive = TRUE,
  showWarnings = FALSE
)
DEGs <- read.csv("results/Differential_Expression/All_DEGs.csv")

top10 <- DEGs %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10)
###############################################################
# Marker Heatmap
###############################################################

heatmap_plot <- DoHeatmap(
  seurat_obj,
  features = unique(top10$gene),
  group.by = "SingleR",
  size = 3
) +
  NoLegend() +
  ggtitle("Top 10 Marker Genes of Each Cell Type")
ggsave(
  filename = "figures/Marker_Heatmap/Marker_Heatmap.png",
  plot = heatmap_plot,
  width = 14,
  height = 12,
  dpi = 300
)

cat("=====================================\n")
cat("Marker Heatmap Created Successfully!\n")
cat("=====================================\n")
