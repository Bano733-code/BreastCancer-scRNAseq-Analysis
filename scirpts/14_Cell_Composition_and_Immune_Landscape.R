###############################################################
## Cell Composition & Immune Landscape
###############################################################

library(Seurat)
library(dplyr)
library(ggplot2)

dir.create(
  "results/Immune_Landscape",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "figures/Immune_Landscape",
  recursive = TRUE,
  showWarnings = FALSE
)
###############################################################
## Major Cell Categories
###############################################################

immune_cells <- c(
  
  "T_cells",
  "B_cell",
  "NK_cell",
  "Macrophage",
  "Monocyte",
  "DC",
  "Neutrophils"
  
)

stromal_cells <- c(
  
  "Fibroblasts",
  "Endothelial_cells",
  "Smooth_muscle_cells",
  "MSC"
  
)

tumor_cells <- c(
  
  "Epithelial_cells"
  
)

###############################################################

celltype <- as.character(seurat_obj$SingleR)

MajorGroup <- ifelse(
  
  celltype %in% immune_cells,
  
  "Immune",
  
  ifelse(
    
    celltype %in% stromal_cells,
    
    "Stromal",
    
    ifelse(
      
      celltype %in% tumor_cells,
      
      "Tumor",
      
      "Other"
      
    )
    
  )
  
)

seurat_obj$MajorGroup <- MajorGroup
###############################################################
## Cell Counts
###############################################################

composition <-
  
  as.data.frame(
    
    table(seurat_obj$MajorGroup)
    
  )

colnames(composition) <-
  
  c(
    
    "Group",
    
    "Cells"
    
  )

composition$Percentage <-
  
  round(
    
    composition$Cells/
      
      sum(composition$Cells)*100,
    
    2
    
  )

print(composition)

write.csv(
  
  composition,
  
  "results/Immune_Landscape/Major_Cell_Composition.csv",
  
  row.names=FALSE
  
)
###############################################################
## Bar Plot
###############################################################

bar <- ggplot(
  
  composition,
  
  aes(
    
    x=reorder(Group,Cells),
    
    y=Cells,
    
    fill=Group
    
  )
  
)+
  
  geom_col()+
  
  coord_flip()+
  
  theme_bw(base_size=14)+
  
  labs(
    
    title="Major Cell Populations",
    
    x="",
    
    y="Number of Cells"
    
  )

ggsave(
  
  "figures/Immune_Landscape/Major_Cell_Composition.png",
  
  bar,
  
  width=8,
  
  height=6,
  
  dpi=300
  
)
###############################################################
## Pie Chart
###############################################################

pie <- ggplot(
  
  composition,
  
  aes(
    
    x="",
    
    y=Percentage,
    
    fill=Group
    
  )
  
)+
  
  geom_col(width=1)+
  
  coord_polar("y")+
  
  theme_void()+
  
  labs(
    
    title="Tumor Microenvironment Composition"
    
  )

ggsave(
  
  "figures/Immune_Landscape/TME_PieChart.png",
  
  pie,
  
  width=7,
  
  height=7,
  
  dpi=300,
  
  bg="white"
  
)
###############################################################
## Immune Cell Types
###############################################################

immune_summary <-
  
  as.data.frame(
    
    table(
      
      seurat_obj$SingleR[
        
        seurat_obj$MajorGroup=="Immune"
        
      ]
      
    )
    
  )

colnames(
  
  immune_summary
  
)<-
  
  c(
    
    "Immune_Cell",
    
    "Cells"
    
  )

immune_summary<-
  
  immune_summary %>%
  
  arrange(desc(Cells))

print(immune_summary)

write.csv(
  
  immune_summary,
  
  "results/Immune_Landscape/Immune_Cell_Subtypes.csv",
  
  row.names=FALSE
  
)
###############################################################
## Immune Landscape
###############################################################

immune_plot<-
  
  ggplot(
    
    immune_summary,
    
    aes(
      
      x=reorder(
        
        Immune_Cell,
        
        Cells
        
      ),
      
      y=Cells
      
    )
    
  )+
  
  geom_col(fill="steelblue")+
  
  coord_flip()+
  
  theme_bw(base_size=13)+
  
  labs(
    
    title="Immune Landscape",
    
    x="",
    
    y="Cells"
    
  )

ggsave(
  
  "figures/Immune_Landscape/Immune_Landscape.png",
  
  immune_plot,
  
  width=8,
  
  height=6,
  
  dpi=300
  
)
cat("Cell Composition Completed!\n")
cat("Immune Landscape Completed!\n")
