# ============================================================
# Derive a gene signature from:
#   (A) the malignant state most strongly associated with subtype
#       (from Phase 1 Step 3/4 results)
#   (B) the CellChat pathway most differentially active across
#       subtypes (from Phase 2 Step 7 results)
# ============================================================

library(Seurat)
library(CellChat)
library(dplyr)


# ------------------------------------------------------------
# 1. Identify the malignant state most associated with subtype
#    (lowest adjusted p-value from the Kruskal-Wallis results)
# ------------------------------------------------------------

kw_df <- read.csv("results/TME/malignant_state_kruskal_BHadjusted.csv")
kw_df <- kw_df %>% arrange(p_adj)

cat("Malignant states ranked by subtype association:\n")
print(kw_df)

top_state <- kw_df$state[1]
cat("\n>>> Top subtype-associated malignant state:", top_state, "<<<\n")


# ------------------------------------------------------------
# 2. Identify the pathway most differentially active across
#    subtypes (from the properly-sorted Step 7 output)
# ------------------------------------------------------------

sig_pathways <- read.csv("results/TME/cellchat_significant_pathways_sorted.csv")

cat("\nTop significant pathways by contribution:\n")
print(head(sig_pathways, 10))

top_pathway <- sig_pathways$name[1]
cat("\n>>> Top differentially-active pathway:", top_pathway, "<<<\n")


# ============================================================
# SIGNATURE A — Marker genes of the top malignant state
# (standard approach: differential expression, this state vs.
#  all other malignant cells)
# ============================================================

malignant_obj <- readRDS("results/TME/malignant_states_with_subtype.rds")

Idents(malignant_obj) <- malignant_obj$malignant_state

state_markers <- FindMarkers(
  malignant_obj,
  ident.1 = top_state,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

state_markers$gene <- rownames(state_markers)
state_markers <- state_markers %>%
  filter(p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC))

cat("\nTop 20 markers for state:", top_state, "\n")
print(head(state_markers[, c("gene", "avg_log2FC", "p_val_adj")], 20))

write.csv(state_markers, "results/TME/gene_signature_state_based.csv",
          row.names = FALSE)

# Take top 50 genes by log fold-change for the signature
signature_state_genes <- head(state_markers$gene, 50)


# ============================================================
# SIGNATURE B — Genes participating in the top pathway,
# restricted to those actually expressed in malignant cells
# ============================================================

# NOTE: mergeCellChat() does NOT store the individual objects
# inside the merged object (no @object.list slot exists on it).
# Load one of the individual subtype objects instead — all three
# used the same CellChatDB, so the gene membership of a pathway
# is identical regardless of which one we use here.
cellchat_single <- readRDS("results/TME/cellchat_TNBC.rds")

# Extract ligand-receptor gene members of the top pathway from
# the CellChatDB interaction table used in the analysis
pathway_genes <- extractEnrichedLR(
  cellchat_single,
  signaling = top_pathway,
  geneLR.return = TRUE
)$geneLR

cat("\nGenes involved in pathway", top_pathway, ":\n")
print(pathway_genes)

rm(cellchat_single); gc()

# Keep only genes actually present/expressed in the malignant object
pathway_genes_present <- intersect(pathway_genes, rownames(malignant_obj))

write.csv(data.frame(gene = pathway_genes_present),
          "results/TME/gene_signature_pathway_based.csv",
          row.names = FALSE)


final_signature <- unique(c(signature_state_genes, pathway_genes_present))

cat("\n============================================\n")
cat("FINAL GENE SIGNATURE\n")
cat("State-based genes (top", length(signature_state_genes), "):", top_state, "\n")
cat("Pathway-based genes (", length(pathway_genes_present), "):", top_pathway, "\n")
cat("Combined signature size:", length(final_signature), "genes\n")
cat("============================================\n")

write.csv(data.frame(gene = final_signature),
          "results/TME/gene_signature_FINAL.csv", row.names = FALSE)

cat("\nSaved: results/TME/gene_signature_FINAL.csv\n")
cat("This file will be used as input for (bulk validation).\n")