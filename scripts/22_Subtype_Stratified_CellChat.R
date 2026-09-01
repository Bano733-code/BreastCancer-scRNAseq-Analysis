# ============================================================
# Subtype-Stratified CellChat Comparison (TNBC / HER2+ / ER+)
# ============================================================

library(Seurat)
library(CellChat)
library(patchwork)
library(dplyr)

set.seed(123)

candidate_pathways <- c("TGFb", "SPP1", "CXCL", "CCL", "MIF")


# ------------------------------------------------------------
# 1. Reusable function: run CellChat for one subtype
#    (only called if the cached .rds doesn't already exist)
# ------------------------------------------------------------

run_cellchat_subtype <- function(seurat_obj, subtype_label,
                                 min_group_size = 15,
                                 max_cells_per_group = 300,
                                 nboot = 25) {
  
  cat("\nRunning CellChat for subtype:", subtype_label, "\n")
  
  sub_obj <- subset(seurat_obj, subset = subtype == subtype_label)
  
  celltype_counts <- table(sub_obj$final_celltype)
  keep_types <- names(celltype_counts)[celltype_counts >= min_group_size]
  sub_obj <- subset(sub_obj, subset = final_celltype %in% keep_types)
  
  set.seed(123)
  cells_keep <- unlist(lapply(unique(sub_obj$final_celltype), function(ct) {
    cells_in_type <- colnames(sub_obj)[sub_obj$final_celltype == ct]
    if (length(cells_in_type) > max_cells_per_group) {
      sample(cells_in_type, max_cells_per_group)
    } else {
      cells_in_type
    }
  }))
  sub_obj <- subset(sub_obj, cells = cells_keep)
  
  data_input <- GetAssayData(sub_obj, assay = "RNA", layer = "data")
  meta <- data.frame(labels = sub_obj$final_celltype, row.names = colnames(sub_obj))
  rm(sub_obj); gc()
  
  cc <- createCellChat(object = data_input, meta = meta, group.by = "labels")
  cc@DB <- CellChatDB.human
  
  cc <- subsetData(cc)
  cc <- identifyOverExpressedGenes(cc)
  cc <- identifyOverExpressedInteractions(cc)
  cc <- computeCommunProb(cc, raw.use = TRUE, nboot = nboot)
  cc <- filterCommunication(cc, min.cells = 10)
  cc <- computeCommunProbPathway(cc)
  cc <- aggregateNet(cc)
  cc <- netAnalysis_computeCentrality(cc, slot.name = "netP")
  
  gc()
  return(cc)
}


# ------------------------------------------------------------
# 2. Load cached results if present, otherwise compute fresh
# ------------------------------------------------------------

get_or_compute <- function(path, subtype_label, seurat_obj) {
  if (file.exists(path)) {
    cat("Loading cached:", path, "\n")
    return(readRDS(path))
  } else {
    cc <- run_cellchat_subtype(seurat_obj, subtype_label)
    saveRDS(cc, path)
    return(cc)
  }
}

paths <- c(
  TNBC = "results/TME/cellchat_TNBC.rds",
  HER2 = "results/TME/cellchat_HER2.rds",
  ER   = "results/TME/cellchat_ER.rds"
)

if (!all(file.exists(paths))) {
  seurat_obj <- readRDS("results/TME/full_annotated_for_cellchat.rds")
  subtype_meta <- read.csv("data/GSE176078/GSE176078/metadata.csv")
  sample_subtype <- unique(subtype_meta[, c("orig.ident", "subtype")])
  seurat_obj$subtype <- sample_subtype$subtype[
    match(seurat_obj$orig.ident, sample_subtype$orig.ident)
  ]
} else {
  seurat_obj <- NULL   # not needed — everything cached
}

cellchat_TNBC <- get_or_compute(paths["TNBC"], "TNBC", seurat_obj)
cellchat_HER2 <- get_or_compute(paths["HER2"], "HER2+", seurat_obj)
cellchat_ER   <- get_or_compute(paths["ER"],   "ER+",   seurat_obj)

rm(seurat_obj); gc()


# ------------------------------------------------------------
# 3. Harmonize cell-type groups across all 3 objects (THE FIX)
#    — required before merging/comparing, since the three
#    subtypes don't contain identical sets of rare cell types
# ------------------------------------------------------------

group.new <- Reduce(union, list(
  levels(cellchat_TNBC@idents),
  levels(cellchat_HER2@idents),
  levels(cellchat_ER@idents)
))

cellchat_TNBC <- liftCellChat(cellchat_TNBC, group.new)
cellchat_HER2 <- liftCellChat(cellchat_HER2, group.new)
cellchat_ER   <- liftCellChat(cellchat_ER, group.new)


# ------------------------------------------------------------
# 4. Merge
# ------------------------------------------------------------

object.list <- list(TNBC = cellchat_TNBC, HER2 = cellchat_HER2, ER = cellchat_ER)
cellchat_merged <- mergeCellChat(object.list, add.names = names(object.list))

saveRDS(cellchat_merged, "results/TME/cellchat_merged_comparison_FIXED.rds")


# ------------------------------------------------------------
# 5. Overall interaction strength (bar chart)
# ------------------------------------------------------------

pdf("results/TME/cellchat_subtype_interaction_strength.pdf", width = 8, height = 5)
print(compareInteractions(cellchat_merged, show.legend = FALSE, group = c(1, 2, 3)))
dev.off()


# ------------------------------------------------------------
# 6. rankNet — full pathway comparison
# ------------------------------------------------------------

rank_data <- rankNet(cellchat_merged, mode = "comparison", stacked = TRUE,
                     do.stat = TRUE, return.data = TRUE)

write.csv(rank_data$signaling.contribution,
          "results/TME/cellchat_subtype_pathway_rank_table.csv",
          row.names = FALSE)

pdf("results/TME/cellchat_subtype_pathway_rank_FULL.pdf", width = 9, height = 14)
print(rankNet(cellchat_merged, mode = "comparison", stacked = TRUE, do.stat = TRUE))
dev.off()


# ------------------------------------------------------------
# 7. PROPERLY interpret rank_data — sort by actual contribution,
#    not the raw unsorted table (head() alone is misleading,
#    since the raw table isn't ordered by importance)
# ------------------------------------------------------------

rc <- rank_data$signaling.contribution

cat("\n=== TOP 15 pathways by contribution (descending) ===\n")
print(rc %>% arrange(desc(contribution)) %>% head(15))

cat("\n=== Significant (p<0.05) pathways with contribution > 0, sorted ===\n")
sig_pathways <- rc %>%
  filter(pvalues < 0.05, contribution > 0) %>%
  arrange(desc(contribution))
print(sig_pathways)

write.csv(sig_pathways, "results/TME/cellchat_significant_pathways_sorted.csv",
          row.names = FALSE)

cat("\n=== Candidate pathways (TGFb/SPP1/CXCL/CCL/MIF) across subtypes ===\n")
candidate_table <- rc %>% filter(name %in% candidate_pathways) %>%
  arrange(name, group)
print(candidate_table)

write.csv(candidate_table, "results/TME/cellchat_candidate_pathways_table.csv",
          row.names = FALSE)


# ------------------------------------------------------------
# 8. Bubble plot — TME senders -> malignant states, per subtype
# ------------------------------------------------------------

common_groups <- Reduce(intersect, lapply(object.list, function(x) levels(x@idents)))
malignant_states <- intersect(
  c("Basal-like", "Luminal Progenitor-like", "Mature Luminal-like",
    "EMT/Mesenchymal", "Cycling/Proliferative"),
  common_groups
)
non_malignant <- setdiff(common_groups, c(malignant_states, "Normal_epithelial"))

pdf("results/TME/cellchat_subtype_bubble_comparison.pdf", width = 12, height = 10)
print(netVisual_bubble(
  cellchat_merged,
  sources.use = non_malignant,
  targets.use = malignant_states,
  comparison = c(1, 2, 3),
  angle.x = 45
))
dev.off()


cat("\n============================================\n")
cat("Key result files:\n")
cat("  - cellchat_significant_pathways_sorted.csv   <- properly sorted,\n")
cat("      significant pathways only. TOP of this file = your\n")
cat("      strongest candidate for the 'mechanistic axis'.\n")
cat("  - cellchat_candidate_pathways_table.csv      <- TGFb/SPP1/CXCL/\n")
cat("      CCL/MIF values across TNBC/HER2+/ER+, side by side.\n")
cat("  - cellchat_subtype_bubble_comparison.pdf\n")
cat("  - cellchat_subtype_pathway_rank_FULL.pdf\n")
cat("============================================\n")