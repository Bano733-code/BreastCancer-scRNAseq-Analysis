# ============================================================
# Build the final annotated object (all 94,195 cells) with:
#   - immune/stromal cell types (from SingleR)
#   - normal epithelial cells (CNV-low)
#   - malignant epithelial cells labeled by STATE (from Step 3/4)
#     instead of a generic "Epithelial" label
#
# This is required before CellChat, since CellChat needs one
# cell-type label per cell across the WHOLE object, not just
# the malignant subset.
# ============================================================

library(Seurat)
library(dplyr)

# ------------------------------------------------------------
# 1. Load the full 94,195-cell object
# ------------------------------------------------------------

seurat_obj <- readRDS("results/TME/seurat_TME.rds")
cat("Full object:\n"); print(dim(seurat_obj))


# ------------------------------------------------------------
# 2. Load Step 2 (CNV call) and Step 3/4 (malignant state) results
# ------------------------------------------------------------

epi          <- readRDS("results/TME/epithelial_cnv_classified.rds")
malignant_obj <- readRDS("results/TME/malignant_states_with_subtype.rds")

# barcode -> cnv_call lookup (all 22,729 epithelial cells)
cnv_lookup <- setNames(epi$cnv_call, colnames(epi))

# barcode -> malignant_state lookup (only the 10,548 CNV-high cells)
state_lookup <- setNames(malignant_obj$malignant_state, colnames(malignant_obj))


# ------------------------------------------------------------
# 3. Build the final cell-type label for EVERY cell in seurat_obj
#
#    Logic:
#      - non-epithelial cells  -> keep their SingleR label
#      - epithelial + CNV-low  -> "Normal_epithelial"
#      - epithelial + CNV-high -> its malignant_state
#                                  (Basal-like, EMT/Mesenchymal, etc.)
# ------------------------------------------------------------

all_barcodes <- colnames(seurat_obj)

final_celltype <- as.character(seurat_obj$SingleR)   # start from SingleR labels
names(final_celltype) <- all_barcodes

is_epithelial <- all_barcodes %in% names(cnv_lookup)

# Normal epithelial (CNV-low)
cnv_low_barcodes <- names(cnv_lookup)[cnv_lookup == "CNV-low"]
final_celltype[cnv_low_barcodes] <- "Normal_epithelial"

# Malignant epithelial -> specific state
malignant_barcodes <- names(state_lookup)
final_celltype[malignant_barcodes] <- state_lookup[malignant_barcodes]

seurat_obj$final_celltype <- final_celltype

cat("\nFinal cell-type composition:\n")
print(table(seurat_obj$final_celltype))


# ------------------------------------------------------------
# 4. Sanity checks
# ------------------------------------------------------------

# Every epithelial cell should now be either Normal_epithelial
# or one of the 5 malignant states — never left as raw "Epithelial_cells"
stopifnot(
  sum(seurat_obj$final_celltype == "Epithelial_cells") == 0
)

# Total counts should still add up
cat("\nTotal cells:", ncol(seurat_obj), "\n")
cat("Sum check (should match):",
    sum(table(seurat_obj$final_celltype)), "\n")


# ------------------------------------------------------------
# 5. Simplify rare/ambiguous SingleR categories (optional but
#    recommended) — CellChat works better with fewer, cleaner
#    groups. Cells with very small group sizes contribute little
#    and can slow down / destabilize the permutation testing.
# ------------------------------------------------------------

celltype_counts <- table(seurat_obj$final_celltype)
print(celltype_counts[order(celltype_counts)])

# Cells belonging to groups smaller than this will be dropped
# from the CellChat object entirely (not deleted from seurat_obj)
min_group_size <- 20
rare_groups <- names(celltype_counts)[celltype_counts < min_group_size]
cat("\nRare groups (<", min_group_size, "cells), will be excluded from CellChat:\n")
print(rare_groups)


# ------------------------------------------------------------
# 6. Save
# ------------------------------------------------------------

saveRDS(seurat_obj, "results/TME/full_annotated_for_cellchat.rds")

cat("\n============================================\n")
cat("FINAL ANNOTATED OBJECT READY FOR CELLCHAT\n")
cat("Saved: results/TME/full_annotated_for_cellchat.rds\n")
cat("============================================\n")# ============================================================
# PHASE 2 — STEP 5 & STEP 6
# Run CellChat on the full annotated object (immune + stromal +
# normal epithelial + malignant states) and identify the
# dominant signaling axis.
#
# Input : results/TME/full_annotated_for_cellchat.rds
# Output: results/TME/cellchat_full.rds + diagnostic plots
# ============================================================

library(Seurat)
library(CellChat)
library(patchwork)
library(dplyr)

setwd("C:/Users/banor/Downloads/BreastCancer_scRNAseq/Single_cell_RNAseq")
set.seed(123)


# ------------------------------------------------------------
# 1. Load the fully annotated object
# ------------------------------------------------------------

seurat_obj <- readRDS("results/TME/full_annotated_for_cellchat.rds")
cat("Full object:\n"); print(dim(seurat_obj))
print(table(seurat_obj$final_celltype))


# ------------------------------------------------------------
# 2. Drop very rare groups + downsample large groups
#
#    WHY DOWNSAMPLE: CellChat's permutation testing scales with
#    cell numbers, and your machine has already hit RAM limits
#    on smaller tasks. Capping each group at max_cells keeps the
#    run memory-safe while preserving every cell TYPE (we are
#    not removing any group, just balancing sizes within them).
# ------------------------------------------------------------

min_group_size <- 20
max_cells_per_group <- 800   # lower this further (e.g. 400) if it still crashes

celltype_counts <- table(seurat_obj$final_celltype)
keep_types <- names(celltype_counts)[celltype_counts >= min_group_size]

seurat_sub <- subset(seurat_obj, subset = final_celltype %in% keep_types)

set.seed(123)
cells_keep <- unlist(lapply(unique(seurat_sub$final_celltype), function(ct) {
  cells_in_type <- colnames(seurat_sub)[seurat_sub$final_celltype == ct]
  if (length(cells_in_type) > max_cells_per_group) {
    sample(cells_in_type, max_cells_per_group)
  } else {
    cells_in_type
  }
}))

seurat_cc <- subset(seurat_sub, cells = cells_keep)

cat("\nCellChat input object:\n")
print(dim(seurat_cc))
print(table(seurat_cc$final_celltype))

rm(seurat_obj, seurat_sub); gc()


# ------------------------------------------------------------
# 3. Create the CellChat object
# ------------------------------------------------------------

data_input <- GetAssayData(seurat_cc, assay = "RNA", layer = "data")
meta <- data.frame(
  labels = seurat_cc$final_celltype,
  row.names = colnames(seurat_cc)
)

cellchat <- createCellChat(object = data_input, meta = meta, group.by = "labels")

rm(seurat_cc, data_input); gc()


# ------------------------------------------------------------
# 4. Set the ligand-receptor database (human)
# ------------------------------------------------------------

CellChatDB <- CellChatDB.human
cellchat@DB <- CellChatDB


# ------------------------------------------------------------
# 5. Preprocessing
# ------------------------------------------------------------

cellchat <- subsetData(cellchat)   # subsets to signaling-related genes only (fast + memory-light)

# Parallel processing (optional — comment out if it causes issues)
# future::plan("multisession", workers = 4)

cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

gc()


# ------------------------------------------------------------
# 6. Compute communication probabilities
#    (this is the slow/heavy step — be patient)
# ------------------------------------------------------------

cellchat <- computeCommunProb(cellchat, raw.use = TRUE)

# Filter out cell-cell pairs with too few cells to be reliable
cellchat <- filterCommunication(cellchat, min.cells = 10)

gc()


# ------------------------------------------------------------
# 7. Pathway-level aggregation
# ------------------------------------------------------------

cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

saveRDS(cellchat, "results/TME/cellchat_full.rds")

cat("\n============================================\n")
cat("PHASE 2 STEP 5 COMPLETE — CellChat object saved\n")
cat("============================================\n")


# ------------------------------------------------------------
# 8. STEP 6 — Identify the dominant signaling pathways
# ------------------------------------------------------------

# Which pathways are active at all?
active_pathways <- cellchat@netP$pathways
cat("\nActive signaling pathways detected:\n")
print(active_pathways)

# Overall communication strength heatmap (all pathways combined)
pdf("results/TME/cellchat_overall_heatmap.pdf", width = 8, height = 10)
netVisual_heatmap(cellchat, measure = "weight", color.heatmap = "Reds")
dev.off()

# Outgoing / incoming signaling role heatmap — shows which cell
# types SEND vs RECEIVE the most signal, per pathway
pdf("results/TME/cellchat_signaling_role_heatmap.pdf", width = 10, height = 12)
netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
dev.off()

# Bubble plot: malignant states as RECEIVERS, everything else as SENDERS
# (this is the key "who signals TO the tumor" view)
malignant_states <- c("Basal-like", "Luminal Progenitor-like",
                      "Mature Luminal-like", "EMT/Mesenchymal",
                      "Cycling/Proliferative")

all_groups <- levels(cellchat@idents)
malignant_states <- intersect(malignant_states, all_groups)
non_malignant <- setdiff(all_groups, c(malignant_states, "Normal_epithelial"))

pdf("results/TME/cellchat_bubble_TME_to_malignant.pdf", width = 10, height = 14)
netVisual_bubble(
  cellchat,
  sources.use = non_malignant,
  targets.use = malignant_states,
  remove.isolate = TRUE
)
dev.off()

# Common breast-TME pathways worth checking explicitly, if present:
candidate_pathways <- c("TGFb", "SPP1", "CXCL", "CCL", "MIF")
present_candidates <- intersect(candidate_pathways, active_pathways)

cat("\nCandidate breast-TME pathways found in your data:\n")
print(present_candidates)

if (length(present_candidates) > 0) {
  pdf("results/TME/cellchat_candidate_pathways_bubble.pdf", width = 10, height = 8)
  netVisual_bubble(
    cellchat,
    sources.use = non_malignant,
    targets.use = malignant_states,
    signaling = present_candidates,
    remove.isolate = TRUE
  )
  dev.off()
}

cat("\n============================================\n")
cat("  - cellchat_overall_heatmap.pdf\n")
cat("  - cellchat_signaling_role_heatmap.pdf\n")
cat("  - cellchat_bubble_TME_to_malignant.pdf\n")
cat("  - cellchat_candidate_pathways_bubble.pdf (if any matched)\n")
cat("============================================\n")