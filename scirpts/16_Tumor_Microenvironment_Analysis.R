###############################################################
## Tumor Microenvironment Analysis
###############################################################

library(Seurat)
library(ggplot2)
library(dplyr)

dir.create(
  "results/TME",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "figures/TME",
  recursive = TRUE,
  showWarnings = FALSE
)
###############################################################
## Assign Major Groups
###############################################################

immune <- c(
  "T_cells",
  "B_cell",
  "NK_cell",
  "Macrophage",
  "Monocyte",
  "DC",
  "Neutrophils"
)

stromal <- c(
  "Fibroblasts",
  "Endothelial_cells",
  "Smooth_muscle_cells",
  "MSC"
)

tumor <- c(
  "Epithelial_cells"
)

celltype <- as.character(seurat_obj$SingleR)

seurat_obj$TME <- ifelse(
  celltype %in% immune,
  "Immune",
  ifelse(
    celltype %in% stromal,
    "Stromal",
    ifelse(
      celltype %in% tumor,
      "Tumor",
      "Other"
    )
  )
)

table(seurat_obj$TME)
###############################################################
## UMAP
###############################################################

tme_umap <- DimPlot(
  seurat_obj,
  group.by = "TME",
  label = TRUE,
  repel = TRUE,
  pt.size = 0.3
) +
  ggtitle("Tumor Microenvironment")

ggsave(
  "figures/TME/TME_UMAP.png",
  tme_umap,
  width = 10,
  height = 8,
  dpi = 300
)
###############################################################
## Summary Table
###############################################################

tme_table <- as.data.frame(table(seurat_obj$TME))

colnames(tme_table) <- c("Group","Cells")

tme_table$Percentage <-
  round(
    100*tme_table$Cells/sum(tme_table$Cells),
    2
  )

print(tme_table)

write.csv(
  tme_table,
  "results/TME/TME_Composition.csv",
  row.names = FALSE
)
###############################################################
## Bar Plot
###############################################################

bar <- ggplot(
  tme_table,
  aes(
    x=reorder(Group,Cells),
    y=Cells,
    fill=Group
  )
)+
  
  geom_col(width=0.7)+
  
  coord_flip()+
  
  theme_bw(base_size=14)+
  
  labs(
    title="Tumor Microenvironment Composition",
    x="",
    y="Cells"
  )

ggsave(
  "figures/TME/TME_Barplot.png",
  bar,
  width=8,
  height=6,
  dpi=300
)
###############################################################
## Pie Chart
###############################################################

pie <- ggplot(
  tme_table,
  aes(
    x="",
    y=Percentage,
    fill=Group
  )
)+
  
  geom_col(width=1)+
  
  coord_polar("y")+
  
  theme_void()+
  
  labs(title="Tumor Microenvironment")

ggsave(
  "figures/TME/TME_PieChart.png",
  pie,
  width=7,
  height=7,
  dpi=300,
  bg="white"
)
###############################################################
## Immune Fraction
###############################################################

immune_only <- subset(
  seurat_obj,
  subset = TME == "Immune"
)

immune_fraction <-
  as.data.frame(
    table(
      immune_only$SingleR
    )
  )

colnames(
  immune_fraction
)<-
  c(
    "CellType",
    "Cells"
  )

immune_fraction<-
  immune_fraction%>%
  arrange(desc(Cells))

print(immune_fraction)

write.csv(
  immune_fraction,
  "results/TME/Immune_Cell_Fraction.csv",
  row.names=FALSE
)
###############################################################
## Immune Landscape
###############################################################

immune_plot<-
  
  ggplot(
    
    immune_fraction,
    
    aes(
      
      x=reorder(CellType,Cells),
      
      y=Cells
      
    )
    
  )+
  
  geom_col(fill="steelblue")+
  
  coord_flip()+
  
  theme_bw(base_size=13)+
  
  labs(
    
    title="Immune Cell Landscape",
    
    x="",
    
    y="Cells"
    
  )

ggsave(
  
  "figures/TME/Immune_Landscape.png",
  
  immune_plot,
  
  width=8,
  
  height=6,
  
  dpi=300
)
###############################################################
## Save Updated Object
###############################################################

saveRDS(
  seurat_obj,
  "results/TME/seurat_TME.rds"
)

cat("\n=====================================\n")
cat("Tumor Microenvironment Analysis Complete!\n")
cat("=====================================\n")