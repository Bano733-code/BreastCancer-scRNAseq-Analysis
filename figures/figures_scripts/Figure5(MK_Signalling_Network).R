# ============================================================
# MANUSCRIPT FIGURE 5 (Cell-cell communication and midkine signaling)
#   Fig. 5a  overall communication landscape (interaction strength heatmap)
#            -> script 21 (results/TME/cellchat_overall_heatmap.pdf)
#   Fig. 5b  relative information flow across subtypes (rankNet)
#            -> script 22 (results/TME/cellchat_subtype_pathway_rank_FULL.pdf)
#   Fig. 5c  MK ligand-receptor bubble plot per subtype -> section "Fig. 5c" at the END of this file
#   Fig. 5d  MK signaling networks per subtype -> this script (main part)
# Run from the project root (no setwd() is used).
# ============================================================

# ============================================================
# FIGURE 5D — SINGLE COMPOSITE PNG
# TNBC vs HER2+ vs ER+
# ============================================================

library(CellChat)
library(magick)

# ------------------------------------------------------------
# 1. Load CellChat objects
# ------------------------------------------------------------

cc_tnbc <- readRDS("results/TME/cellchat_TNBC.rds")
cc_her2 <- readRDS("results/TME/cellchat_HER2.rds")
cc_er   <- readRDS("results/TME/cellchat_ER.rds")

# ------------------------------------------------------------
# 2. Relevant breast-TME groups
# ------------------------------------------------------------

relevant_groups <- c(
  "Basal-like",
  "Luminal Progenitor-like",
  "Mature Luminal-like",
  "EMT/Mesenchymal",
  "Cycling/Proliferative",
  "Normal_epithelial",
  "T_cells",
  "B_cell",
  "NK_cell",
  "Macrophage",
  "Monocyte",
  "DC",
  "Fibroblasts",
  "Endothelial_cells",
  "Smooth_muscle_cells",
  "MSC"
)

# ------------------------------------------------------------
# 3. Groups common to all 3 objects
# ------------------------------------------------------------

common_groups <- Reduce(
  intersect,
  list(
    levels(cc_tnbc@idents),
    levels(cc_her2@idents),
    levels(cc_er@idents)
  )
)

plot_groups <- intersect(
  relevant_groups,
  common_groups
)

cat("\nGroups included:\n")
print(plot_groups)

# ------------------------------------------------------------
# 4. Output directory
# ------------------------------------------------------------

dir.create(
  "figures/TME",
  recursive = TRUE,
  showWarnings = FALSE
)

# ============================================================
# 5. Function: make individual PNG
# ============================================================

make_mk_png <- function(
    cc,
    filename,
    panel_title
) {

  png(
    filename,
    width = 1800,
    height = 1800,
    res = 300,
    bg = "white"
  )

  netVisual_aggregate(
    cc,
    signaling = "MK",
    signaling.name = panel_title,
    layout = "circle",
    sources.use = plot_groups,
    targets.use = plot_groups,
    remove.isolate = TRUE,
    vertex.label.cex = 0.55
  )

  dev.off()
}

# ============================================================
# 6. Generate temporary subtype PNGs
# ============================================================

tnbc_file <- "figures/TME/_MK_TNBC_temp.png"
her2_file <- "figures/TME/_MK_HER2_temp.png"
er_file   <- "figures/TME/_MK_ER_temp.png"

make_mk_png(
  cc_tnbc,
  tnbc_file,
  "TNBC"
)

make_mk_png(
  cc_her2,
  her2_file,
  "HER2+"
)

make_mk_png(
  cc_er,
  er_file,
  "ER+"
)

# ============================================================
# 7. Read the three images
# ============================================================

img_tnbc <- image_read(tnbc_file)
img_her2 <- image_read(her2_file)
img_er   <- image_read(er_file)

# ============================================================
# 8. Make all images same size
# ============================================================

img_tnbc <- image_resize(
  img_tnbc,
  "1800x1800!"
)

img_her2 <- image_resize(
  img_her2,
  "1800x1800!"
)

img_er <- image_resize(
  img_er,
  "1800x1800!"
)

# ============================================================
# 9. Combine side-by-side
# ============================================================

combined <- image_append(
  c(
    img_tnbc,
    img_her2,
    img_er
  ),
  stack = FALSE
)

# ============================================================
# 10. Add overall title area
# ============================================================

title_area <- image_blank(
  width = 5400,
  height = 300,
  color = "white"
)

title_area <- image_annotate(
  title_area,
  "Midkine (MK) signaling across breast cancer subtypes",
  gravity = "center",
  location = "+0+0",
  size = 70,
  weight = 700,
  color = "black"
)

# ============================================================
# 11. Combine title + networks
# ============================================================

final_image <- image_append(
  c(
    title_area,
    combined
  ),
  stack = TRUE
)

# ============================================================
# 12. Add panel labels
# ============================================================

# Add labels using image annotation
final_image <- image_annotate(
  final_image,
  "A. TNBC",
  location = "+220+330",
  size = 45,
  weight = 700,
  color = "black"
)

final_image <- image_annotate(
  final_image,
  "B. HER2+",
  location = "+2020+330",
  size = 45,
  weight = 700,
  color = "black"
)

final_image <- image_annotate(
  final_image,
  "C. ER+",
  location = "+3820+330",
  size = 45,
  weight = 700,
  color = "black"
)

# ============================================================
# 13. Save final ONE PNG
# ============================================================

final_file <- "figures/TME/Figure5D_MK_Subtype_Networks.png"

image_write(
  final_image,
  path = final_file,
  format = "png"
)

# ============================================================
# 14. Delete temporary PNGs
# ============================================================

file.remove(
  tnbc_file,
  her2_file,
  er_file
)

cat(
  "\n============================================\n",
  "Figure 5D generated successfully!\n",
  "============================================\n"
)

cat(
  "\nFinal file:\n",
  final_file,
  "\n"
)


# ============================================================
# [REVISION] Fig. 5c (final): MK bubble plot (TNBC / HER2 / ER), bold axis text
# NOTE: MSC-labelled cells are almost entirely cancer-epithelial cells in the
# original atlas (97.3%), so the "MSC" sender in TNBC is best read as malignant cells
# (probable autocrine MK signaling); the group is kept to match the manuscript figure
# and is explained in the caption.
# ============================================================

library(CellChat)
library(ggplot2)
library(dplyr)

dir.create("figures/Manuscript_Figures", recursive = TRUE, showWarnings = FALSE)
dir.create("figures/TME", recursive = TRUE, showWarnings = FALSE)

cellchat_merged <- readRDS("results/TME/cellchat_merged_comparison_FIXED.rds")

malignant_states <- c("Basal-like", "Luminal Progenitor-like", "Mature Luminal-like",
                      "EMT/Mesenchymal", "Cycling/Proliferative")

senders <- c("Fibroblasts", "Endothelial_cells", "MSC", "Macrophage",
             "NK_cell", "Smooth_muscle_cells", "Normal_epithelial")

p0 <- netVisual_bubble(cellchat_merged,
                       sources.use = senders, targets.use = malignant_states,
                       comparison = c(1, 2, 3), signaling = "MK",
                       remove.isolate = FALSE)

d <- p0$data
d$pair    <- gsub(" \\((TNBC|HER2|ER)\\)$", "", as.character(d$source.target))
d$subtype <- factor(d$dataset, levels = c("TNBC", "HER2", "ER"))

p5c <- ggplot(d, aes(x = interaction_name_2, y = pair, colour = prob, size = pval)) +
  geom_point() +
  facet_grid(. ~ subtype) +
  scale_colour_distiller(palette = "Spectral", na.value = "grey90") +
  scale_size(range = c(1.5, 4.5), guide = "none") +
  labs(x = "Ligand - receptor", y = "Sender -> target", colour = "Commun.\nProb.") +
  theme_bw(base_size = 11) +
  theme(axis.text.x  = element_text(angle = 45, hjust = 1, size = 9, face = "bold"),
        axis.text.y  = element_text(size = 8, face = "bold"),
        axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        strip.background = element_rect(fill = "grey90"),
        strip.text = element_text(face = "bold"))

ggsave("figures/TME/Figure5C_MK_bubble_swapped_bold.png", p5c,
       width = 12, height = 10, dpi = 300, bg = "white")
ggsave("figures/Manuscript_Figures/Fig5c.png", p5c,
       width = 12, height = 10, dpi = 300, bg = "white")
print(p5c)

# ------------------------------------------------------------
# Optional: bold bubble plots for the full-object CellChat (cellchat_full.rds)
# ------------------------------------------------------------

cellchat <- readRDS("results/TME/cellchat_full.rds")

bold_theme <- theme(
  axis.text.x  = element_text(face = "bold"),
  axis.text.y  = element_text(face = "bold"),
  axis.title.x = element_text(face = "bold"),
  axis.title.y = element_text(face = "bold"),
  strip.text   = element_text(face = "bold")
)

all_groups       <- levels(cellchat@idents)
malignant_states <- intersect(malignant_states, all_groups)
non_malignant    <- setdiff(all_groups, c(malignant_states, "Normal_epithelial"))

# 1) all pathways
p1 <- netVisual_bubble(cellchat, sources.use = non_malignant,
                       targets.use = malignant_states, remove.isolate = TRUE) + bold_theme
ggsave("results/TME/cellchat_bubble_TME_to_malignant_bold.png", p1,
       width = 10, height = 14, dpi = 300, bg = "white")

# 2) candidate pathways
present_candidates <- intersect(c("TGFb", "SPP1", "CXCL", "CCL", "MIF"), cellchat@netP$pathways)
if (length(present_candidates) > 0) {
  p2 <- netVisual_bubble(cellchat, sources.use = non_malignant,
                         targets.use = malignant_states,
                         signaling = present_candidates, remove.isolate = TRUE) + bold_theme
  ggsave("results/TME/cellchat_candidate_pathways_bubble_bold.png", p2,
         width = 10, height = 8, dpi = 300, bg = "white")
}
