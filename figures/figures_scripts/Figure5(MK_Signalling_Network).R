# ============================================================
# FIGURE 5D — MIDKINE (MK) SIGNALING NETWORK
# ============================================================

library(CellChat)

# ------------------------------------------------------------
# 1. Load full CellChat object
# ------------------------------------------------------------

cellchat_full <- readRDS(
  "results/TME/cellchat_full.rds"
)

cat("\nCellChat object loaded successfully.\n")

# ------------------------------------------------------------
# 2. Check available signaling pathways
# ------------------------------------------------------------

cat("\nAvailable pathways containing 'MK':\n")

mk_pathways <- cellchat_full@netP$pathways[
  grepl("^MK$|Midkine|midkine",
        cellchat_full@netP$pathways,
        ignore.case = TRUE)
]

print(mk_pathways)

# ------------------------------------------------------------
# 3. Safety check
# ------------------------------------------------------------

if (length(mk_pathways) == 0) {
  stop(
    "MK signaling pathway was not found in cellchat_full.rds"
  )
}

# CellChat normally uses "MK"
mk_pathway <- mk_pathways[1]

cat("\nUsing pathway:", mk_pathway, "\n")

# ------------------------------------------------------------
# 4. Create output directory
# ------------------------------------------------------------

dir.create(
  "figures/TME",
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# 5. PDF output
# ------------------------------------------------------------

pdf(
  "figures/TME/Figure5D_MK_Signaling_Network.pdf",
  width = 8,
  height = 8
)

netVisual_aggregate(
  cellchat_full,
  signaling = mk_pathway,
  layout = "circle"
)

dev.off()

# ------------------------------------------------------------
# 6. PNG output
# ------------------------------------------------------------

png(
  "figures/TME/Figure5D_MK_Signaling_Network.png",
  width = 2400,
  height = 2400,
  res = 300
)

netVisual_aggregate(
  cellchat_full,
  signaling = mk_pathway,
  layout = "circle"
)

dev.off()

cat(
  "\nFigure 5D generated successfully:\n",
  "figures/TME/Figure5D_MK_Signaling_Network.png\n",
  sep = ""
)