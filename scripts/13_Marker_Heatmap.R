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

# REQUIRED: load the object — this was missing, causing your error
seurat_obj <- readRDS("data/processed/seurat_annotated.rds")

DEGs <- read.csv("results/Differential_Expression/All_DEGs.csv")

top5 <- DEGs %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 5)

###############################################################
# Marker Heatmap
###############################################################
heatmap_plot <- DoHeatmap(
  seurat_obj,
  features = unique(top5$gene),
  group.by = "SingleR",
  size = 3,
  angle = 45,     # prevents label overlap
  hjust = 0
) +
  scale_fill_gradientn(colors = c("#2c115f", "black", "#f1c40f")) +  # better contrast than NoLegend()
  theme(
    axis.text.y = element_text(size = 7),
    legend.position = "right"
  ) +
  ggtitle("Top 5 Marker Genes of Each Cell Type")

ggsave(
  filename = "figures/Marker_Heatmap/Marker_Heatmap.png",
  plot = heatmap_plot,
  width = 18,
  height = 12,
  dpi = 300
)

cat("=====================================\n")
cat("Marker Heatmap Created Successfully!\n")
cat("=====================================\n")
