###############################################################
## Marker DotPlot
###############################################################

library(Seurat)
library(ggplot2)

dir.create(
  "figures/Marker_DotPlot",
  recursive = TRUE,
  showWarnings = FALSE
)

###############################################################
# Canonical markers
###############################################################

markers <- c(
  
  # T cells
  "CD3D","CD3E","IL7R",
  
  # B cells
  "MS4A1","CD79A","CD79B",
  
  # NK cells
  "NKG7","GNLY","KLRD1",
  
  # Macrophages
  "LYZ","C1QA","FCER1G",
  
  # Endothelial
  "VWF","PLVAP","EMCN",
  
  # Fibroblasts
  "COL1A1","COL3A1","DCN",
  
  # Epithelial
  "EPCAM","KRT8","KRT18"
  
)

###############################################################

dot <- DotPlot(
  
  seurat_obj,
  
  features = markers,
  
  group.by = "SingleR"
  
)+
  
  RotatedAxis()+
  
  ggtitle("Canonical Cell-Type Marker Expression")

###############################################################

ggsave(
  
  "figures/Marker_DotPlot/Marker_DotPlot.png",
  
  dot,
  
  width=14,
  
  height=7,
  
  dpi=300,
  
  bg="white"
  
)

cat("Marker DotPlot Finished!\n")