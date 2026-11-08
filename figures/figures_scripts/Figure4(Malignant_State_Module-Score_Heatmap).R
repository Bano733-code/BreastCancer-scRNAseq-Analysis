# ============================================================
# MANUSCRIPT FIGURE 4 (Malignant transcriptional-state heterogeneity)
#   Fig. 4a  UMAP of the 10,548 CNV-high cells coloured by malignant state
#            -> section "Fig. 4a" at the END of this file
#   Fig. 4b  module-score heatmap -> this script (main part)
#   Fig. 4c  stacked bar of state composition by subtype -> script 20
#            (results/TME/malignant_state_by_subtype_stackedbar.png)
#   Fig. 4d  per-sample state proportions -> script 20
#            (results/TME/malignant_state_by_subtype_persample_boxplot.png)
# Run from the project root (no setwd() is used).
# ============================================================

# ============================================================
# FIGURE 4B — MALIGNANT STATE MODULE-SCORE HEATMAP
# ============================================================

library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)

# ------------------------------------------------------------
# 1. Load data
# ------------------------------------------------------------

scores <- read_csv(
  "results/TME/malignant_state_scores_per_cell.csv",
  show_col_types = FALSE
)

# ------------------------------------------------------------
# 2. Exact score columns from YOUR file
# ------------------------------------------------------------

score_cols <- c(
  "Basal_score",
  "LumProg_score",
  "MatureLum_score",
  "EMT_score",
  "Cycling_score"
)

# Check that all columns exist
missing_cols <- setdiff(score_cols, colnames(scores))

if (length(missing_cols) > 0) {
  stop(
    paste(
      "These columns are missing:",
      paste(missing_cols, collapse = ", ")
    )
  )
}

# ------------------------------------------------------------
# 3. Calculate mean module score for each malignant state
# ------------------------------------------------------------

mean_scores <- scores %>%
  summarise(
    Basal_score      = mean(Basal_score, na.rm = TRUE),
    LumProg_score    = mean(LumProg_score, na.rm = TRUE),
    MatureLum_score  = mean(MatureLum_score, na.rm = TRUE),
    EMT_score        = mean(EMT_score, na.rm = TRUE),
    Cycling_score    = mean(Cycling_score, na.rm = TRUE)
  )

print(mean_scores)

# ------------------------------------------------------------
# 4. Convert to long format
# ------------------------------------------------------------

heatmap_df <- mean_scores %>%
  pivot_longer(
    cols = everything(),
    names_to = "State",
    values_to = "MeanScore"
  )

# Human-readable labels
heatmap_df$State <- factor(
  heatmap_df$State,
  levels = score_cols,
  labels = c(
    "Basal-like",
    "Luminal progenitor-like",
    "Mature luminal-like",
    "EMT/mesenchymal",
    "Cycling/proliferative"
  )
)

# ------------------------------------------------------------
# 5. IMPORTANT:
#    We need a 5 x 5 heatmap if we want
#    malignant state x module score.
#
#    So calculate average score FOR EACH malignant state
# ------------------------------------------------------------

state_heatmap <- scores %>%
  group_by(malignant_state) %>%
  summarise(
    Basal_score     = mean(Basal_score, na.rm = TRUE),
    LumProg_score   = mean(LumProg_score, na.rm = TRUE),
    MatureLum_score = mean(MatureLum_score, na.rm = TRUE),
    EMT_score       = mean(EMT_score, na.rm = TRUE),
    Cycling_score   = mean(Cycling_score, na.rm = TRUE),
    .groups = "drop"
  )

# ------------------------------------------------------------
# 6. Convert to long format
# ------------------------------------------------------------

state_heatmap_long <- state_heatmap %>%
  pivot_longer(
    cols = all_of(score_cols),
    names_to = "Module",
    values_to = "Score"
  )

# ------------------------------------------------------------
# 7. Clean labels
# ------------------------------------------------------------

state_heatmap_long$Module <- factor(
  state_heatmap_long$Module,
  levels = score_cols,
  labels = c(
    "Basal-like",
    "Luminal progenitor-like",
    "Mature luminal-like",
    "EMT/mesenchymal",
    "Cycling/proliferative"
  )
)

# Clean malignant-state labels
state_heatmap_long$malignant_state <- factor(
  state_heatmap_long$malignant_state,
  levels = c(
    "Basal-like",
    "Luminal Progenitor-like",
    "Mature Luminal-like",
    "EMT/Mesenchymal",
    "Cycling/Proliferative"
  ),
  labels = c(
    "Basal-like",
    "Luminal progenitor-like",
    "Mature luminal-like",
    "EMT/mesenchymal",
    "Cycling/proliferative"
  )
)

# ------------------------------------------------------------
# 8. Plot heatmap
# ------------------------------------------------------------

p4B <- ggplot(
  state_heatmap_long,
  aes(
    x = Module,
    y = malignant_state,
    fill = Score
  )
) +
  
  geom_tile(
    color = "white",
    linewidth = 0.6
  ) +
  
  geom_text(
    aes(label = sprintf("%.2f", Score)),
    size = 3.5
  ) +
  
  scale_fill_gradient2(
    midpoint = 0,
    name = "Mean\nmodule score"
  ) +
  
  labs(
    title = "Malignant-state module scores",
    x = "Malignant transcriptional program",
    y = "Assigned malignant state"
  ) +
  
  theme_classic(base_size = 12) +
  
  theme(
    axis.text.x = element_text(
      angle = 35,
      hjust = 1
    ),
    axis.text.y = element_text(
      face = "bold"
    ),
    plot.title = element_text(
      face = "bold"
    )
  )

print(p4B)

# ------------------------------------------------------------
# 9. Save figure
# ------------------------------------------------------------

dir.create(
  "figures/Malignant",
  recursive = TRUE,
  showWarnings = FALSE
)

ggsave(
  "figures/Malignant/Figure4B_Malignant_State_ModuleScore_Heatmap.png",
  p4B,
  width = 8,
  height = 6,
  dpi = 600
)

cat("\nFigure 4B generated successfully.\n")


# ============================================================
# [REVISION] Fig. 4a (final): CNV-high malignant cells only (n = 10,548)
# The earlier panel clustered ALL 22,729 epithelial cells and therefore did not match
# the caption. The malignant subset was already re-clustered in script 19
# (30 PCs computed, first 20 used for the neighbour graph and UMAP, resolution 1.0),
# so the stored embedding is reused; it is recomputed with the same settings only if missing.
# ============================================================

library(Seurat)
library(ggplot2)
library(dplyr)

dir.create("figures/Manuscript_Figures", recursive = TRUE, showWarnings = FALSE)

mal <- readRDS("results/TME/malignant_states_with_subtype.rds")
ncol(mal)                                  # should be 10,548 cells

if (!"umap" %in% names(mal@reductions)) {
  DefaultAssay(mal) <- "RNA"
  mal <- NormalizeData(mal, verbose = FALSE)
  mal <- FindVariableFeatures(mal, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  mal <- ScaleData(mal, verbose = FALSE)
  mal <- RunPCA(mal, npcs = 30, verbose = FALSE)
  mal <- FindNeighbors(mal, dims = 1:20, verbose = FALSE)
  mal <- FindClusters(mal, resolution = 1.0, verbose = FALSE)
  mal <- RunUMAP(mal, dims = 1:20, verbose = FALSE)
}

p4a <- DimPlot(mal, reduction = "umap", group.by = "malignant_state") +
  ggtitle("CNV-high malignant cells (n = 10,548)") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

ggsave("figures/Manuscript_Figures/Fig4a.png", p4a,
       width = 8, height = 6, dpi = 300, bg = "white")
print(p4a)

# ------------------------------------------------------------
# Optional supporting plots (not in the final manuscript figure)
# ------------------------------------------------------------

# (1) Sub-clusters of ALL 22,729 epithelial cells, legend removed
epithelial_obj <- readRDS("results/Malignant/epithelial_reclustered.rds")
p_all_epi <- DimPlot(epithelial_obj, reduction = "umap",
                     label = TRUE, repel = TRUE, pt.size = 0.4) +
  NoLegend() +                              # the side legend is removed
  ggtitle("Epithelial Cell Subclusters") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
ggsave("figures/Malignant/Epithelial_Subclusters_UMAP_nolegend.png",
       p_all_epi, width = 8, height = 7, dpi = 300, bg = "white")

# (2) Dominant malignant state of each sub-cluster
print(table(Idents(mal))[1:3])             # check: are these the expected cluster IDs?
print(colnames(mal@meta.data))             # inspect the cluster column and malignant_state

cluster_state_tab <- mal@meta.data |>
  count(cluster = Idents(mal), malignant_state) |>
  group_by(cluster) |>
  mutate(pct = round(100 * n / sum(n), 1)) |>
  slice_max(pct, n = 1, with_ties = FALSE) |>
  arrange(as.numeric(as.character(cluster)))
print(cluster_state_tab, n = 40)
