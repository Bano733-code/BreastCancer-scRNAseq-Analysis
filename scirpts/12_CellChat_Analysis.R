###############################################################
## Cell Composition Analysis
###############################################################

library(Seurat)
library(ggplot2)
library(dplyr)
library(scales)

###############################################################
# Create folders
###############################################################

dir.create("results/Cell_Composition",
           recursive = TRUE,
           showWarnings = FALSE)

dir.create("figures/Cell_Composition",
           recursive = TRUE,
           showWarnings = FALSE)

###############################################################
# Count each cell type
###############################################################

cell_counts <- as.data.frame(table(seurat_obj$SingleR))

colnames(cell_counts) <- c("CellType","Count")

###############################################################
# Calculate percentages
###############################################################

cell_counts$Percentage <-
  round(
    cell_counts$Count /
      sum(cell_counts$Count) * 100,
    2
  )

###############################################################
# Sort from largest to smallest
###############################################################

cell_counts <-
  cell_counts %>%
  arrange(desc(Count))

###############################################################
# Save table
###############################################################

write.csv(
  cell_counts,
  "results/Cell_Composition/Cell_Composition.csv",
  row.names = FALSE
)

###############################################################
# Print table
###############################################################

print(cell_counts)
###############################################################
# Publication-quality Bar Plot
###############################################################

bar_plot <-
  
  ggplot(
    cell_counts,
    aes(
      x = reorder(CellType, Count),
      y = Count
    )
  )+
  
  geom_col() +
  
  coord_flip() +
  
  theme_bw(base_size = 14)+
  
  labs(
    title="Cell Composition",
    x="Cell Type",
    y="Number of Cells"
  )

ggsave(
  "figures/Cell_Composition/Cell_Composition_Barplot.png",
  bar_plot,
  width=9,
  height=7,
  dpi=300
)
###############################################################
# Pie Chart
###############################################################

pie_plot <-
  
  ggplot(
    cell_counts,
    aes(
      x="",
      y=Percentage,
      fill=CellType
    )
  )+
  
  geom_col(width=1)+
  
  coord_polar(theta="y")+
  
  theme_void()+
  
  labs(title="Cell Type Composition")

ggsave(
  "figures/Cell_Composition/Cell_Composition_PieChart.png",
  pie_plot,
  width=8,
  height=8,
  dpi=300,
  bg="white"
)
###############################################################
# Percentage Plot
###############################################################

percent_plot <-
  
  ggplot(
    cell_counts,
    aes(
      x=reorder(CellType,Percentage),
      y=Percentage
    )
  )+
  
  geom_col()+
  
  coord_flip()+
  
  theme_bw(base_size=14)+
  
  labs(
    title="Cell Type Percentage",
    x="Cell Type",
    y="Percentage (%)"
  )

ggsave(
  "figures/Cell_Composition/Cell_Composition_Percentage.png",
  percent_plot,
  width=9,
  height=7,
  dpi=300
)

cat("\n")
cat("=====================================\n")
cat("Cell Composition Analysis Completed!\n")
cat("=====================================\n")