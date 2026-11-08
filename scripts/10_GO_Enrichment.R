library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(enrichplot)
library(ggplot2)

# Read DEG file
DEGs <- read.csv(
  "results/Differential_Expression/All_DEGs.csv"
)

######################################################
# Choose ONE cell type
######################################################

celltype <- "T_cells"

genes <- DEGs %>%
  filter(
    cluster == celltype,
    p_val_adj < 0.05,
    avg_log2FC > 0.25
  )

######################################################
# Convert IDs
######################################################

gene.df <- bitr(
  genes$gene,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

######################################################
# GO
######################################################

GO <- enrichGO(
  gene = gene.df$ENTREZID,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  readable = TRUE
)
if (nrow(as.data.frame(GO)) == 0) {
  stop("No enriched GO terms found.")
}
######################################################
# Save
######################################################

dir.create(
  "results/GO_KEGG",
  recursive = TRUE,
  showWarnings = FALSE
)

write.csv(
  as.data.frame(GO),
  paste0("results/GO_KEGG/", celltype, "_GO.csv"),
  row.names = FALSE
)

######################################################
# Figures
######################################################
dir.create("figures/GO_KEGG", recursive = TRUE, showWarnings = FALSE)

pdf(
  paste0("figures/GO_KEGG/", celltype, "_GO_dotplot.pdf"),
  width = 8,
  height = 6
)

dotplot(GO, showCategory = 15)

dev.off()


p1 <- dotplot(GO, showCategory = 15)

ggsave(
  "figures/GO_KEGG/GO_Dotplot.png",
  p1,
  width = 10,
  height = 8,
  dpi = 300
)

p2 <- barplot(GO, showCategory = 15)

ggsave(
  "figures/GO_KEGG/GO_Barplot.png",
  p2,
  width = 10,
  height = 8,
  dpi = 300
)
cat("Finished!")