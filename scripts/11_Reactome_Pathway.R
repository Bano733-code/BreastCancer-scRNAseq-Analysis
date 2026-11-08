##############################
# Load Libraries
##############################

library(clusterProfiler)
library(ReactomePA)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(dplyr)

##############################
# Create Output Folder
##############################

dir.create(
  "results/Reactome",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "figures/Reactome",
  recursive = TRUE,
  showWarnings = FALSE
)

##############################
# Read Differential Expression Results
##############################

DEGs <- read.csv(
  "results/Differential_Expression/All_DEGs.csv"
)

##############################
# Filter Significant Genes
##############################

genes <- DEGs %>%
  filter(
    p_val_adj < 0.05,
    avg_log2FC > 0.25
  )

cat("Significant genes:", nrow(genes), "\n")

##############################
# Convert Gene Symbols
# to Entrez IDs
##############################

gene.df <- bitr(
  genes$gene,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

gene.df <- gene.df[!duplicated(gene.df$ENTREZID), ]

cat("Mapped genes:", nrow(gene.df), "\n")

##############################
# Reactome Enrichment
##############################

reactome <- enrichPathway(
  
  gene = gene.df$ENTREZID,
  
  organism = "human",
  
  readable = TRUE,
  
  pvalueCutoff = 0.05,
  
  qvalueCutoff = 0.20
  
)

##############################
# Check Results
##############################

cat("Number of enriched pathways:",
    nrow(as.data.frame(reactome)),
    "\n")

##############################
# Save Results
##############################

write.csv(
  
  as.data.frame(reactome),
  
  "results/Reactome/Reactome_Pathways.csv",
  
  row.names = FALSE
  
)

##############################
# Dot Plot
##############################

dot <- dotplot(
  
  reactome,
  
  showCategory = 15,
  
  font.size = 12,
  
  title = "Reactome Pathway Enrichment"
  
)

ggsave(
  
  filename = "figures/Reactome/Reactome_Dotplot.png",
  
  plot = dot,
  
  width = 10,
  
  height = 8,
  
  dpi = 300
  
)

##############################
# Bar Plot
##############################

bar <- barplot(
  
  reactome,
  
  showCategory = 15,
  
  font.size = 12,
  
  title = "Reactome Pathway Enrichment"
  
)

ggsave(
  
  filename = "figures/Reactome/Reactome_Barplot.png",
  
  plot = bar,
  
  width = 10,
  
  height = 8,
  
  dpi = 300
  
)

##############################
# Enrichment Map
##############################

emap <- pairwise_termsim(reactome)

emap_plot <- emapplot(
  
  emap,
  
  showCategory = 20
  
)

ggsave(
  
  filename = "figures/Reactome/Reactome_EnrichmentMap.png",
  
  plot = emap_plot,
  
  width = 12,
  
  height = 10,
  
  dpi = 300,
  
  bg="white"
  
)

##############################
# Display Top Pathways
##############################

print(
  
  head(
    
    as.data.frame(reactome),
    
    20
    
  )
  
)

##############################
# Completion Message
##############################


cat("Reactome Analysis Completed!\n")

cat("Results saved in:\n")

cat("results/Reactome/\n")

cat("Figures saved in:\n")

cat("figures/Reactome/\n")