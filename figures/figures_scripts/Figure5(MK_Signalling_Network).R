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
